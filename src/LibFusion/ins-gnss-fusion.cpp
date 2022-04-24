#include "rtklib.h"

#define UPD_INS_PAR  1

/* Simplfiy the filter usage */
#define KF_FILTER(kf, dz, Qdz)                                               \
    filter(kf->x, kf->P, kf->H, (dz), (Qdz), kf->nx, kf->ny, 0, 0, NULL, 0)

/* set kf->ny and set kf->H to zero */
#define KF_HINIT(kf, _ny)                                                    \
{                                                                            \
    kf->ny = _ny;                                                            \
    for(unsigned int i = 0; i < kf->ny*kf->nx; ++i)                          \
        kf->H[i] = 0.0;                                                      \
}

extern void traceins(const solins_t* sol_ins, int post, const insopt_t* opt)
{
	int week;
	double sow = time2gpst(sol_ins->time, &week);
	v3_t pos, vel, att;
	pos = sol_ins->pos; vel = sol_ins->vel;
	if (opt->mech_coord == INSMECH_ECEF) {
		//        ecef2llh(&pos,&vel, nullptr,opt->local_coord);
		quat2att(&sol_ins->quat, &att, opt->imu_coord);
	}
	fprintf(stdout, "%s(%9.4f) INS MECH(%s)\n", time_str(sol_ins->time, 3), sow, post ? "+" : "-");
	fprintf(stdout, "attitue(quat): %26.16f %26.16f %26.16f %26.16f\n", sol_ins->quat.q0, sol_ins->quat.q1, sol_ins->quat.q2, sol_ins->quat.q3);
	fprintf(stdout, "velocity(m/s): %26.16f %26.16f %26.16f\n", vel.x, vel.y, vel.z);
	fprintf(stdout, "position(llh): %26.16f %26.16f %26.16f\n", pos.x, pos.y, pos.z);
	fflush(stdout);
}

static void save2ygm(const pva_t* pva, const imu_t* imu)
{
	FILE* out_gps, * out_imu;
	int week = 0;
	double sow = 0.0, dt = 0.0, sow0 = 108143.0;

	char path_gps[1024];
	char path_imu[1024];

	strcpy(path_gps, "E:/DataProcess/GNSS_DATA/2019/364/obs_lc/100C.gps");
	strcpy(path_imu, "E:/DataProcess/GNSS_DATA/2019/364/obs_lc/100C.imu");

	out_gps = fopen(path_gps, "w");

	for (int i = 0; i < pva->n; i++) {
		sow = time2gpst(pva->time[i], nullptr);
		dt = sow - sow0;
		fprintf(out_gps, "%20.16G %20.16G %20.16G %20.16G %20.16G %20.16G %20.16G\n", pva->vel[i].x, pva->vel[i].y, pva->vel[i].z,
			pva->pos[i].x, pva->pos[i].y, pva->pos[i].z, dt);
	}

	out_imu = fopen(path_imu, "w");
	for (int i = 0; i < imu->n; i++) {
		sow = time2gpst(imu->data[i].time, nullptr);
		dt = sow - sow0;
		fprintf(out_imu, "%20.16G %26.16G %20.16G %20.16G %20.16G %20.16G %26.16G\n", imu->data[i].gyro.x, imu->data[i].gyro.y, imu->data[i].gyro.z,
			imu->data[i].accel.x, imu->data[i].accel.y, imu->data[i].accel.z, dt);
	}

	fclose(out_gps);
	fclose(out_imu);
}

extern int ins_use_rtk_init(rtk_t* rtk, prcopt_t* popt, const imup_t* imup) {
	int stat = SOLQ_FIX;
	popt->mode = PMODE_KINEMA;
	popt->modear = ARMODE_CONT;
	popt->sateph = EPHOPT_BRDC;
	popt->tropopt = TROPOPT_SAAS;
	popt->ionoopt = IONOOPT_BRDC;
	popt->tropmap = TROPMAP_NMF;
	popt->elmin = 15.0 * D2R;
	popt->nf = 2;
	popt->thresar[0] = 3.0;
	popt->thresar[1] = 0.995;
	popt->thresar[2] = 1.0;
	popt->par = 0;
	popt->navsys = SYS_GPS | SYS_CMP;

	rtkinit(rtk, popt, imup);
	return stat;
}

extern int ins_use_rtk(const prcopt_t* popt, imup_t* imup, const imu_t* imu, const int kimu_idx, int* rover_idx, int* base_idx,
	const obs_t* obs, const nav_t* nav, rtk_t* rtk) {
	obsd_t epoch_obs[MAXOBS * 2] = { {0} };
	int i, nobs = 0, n = 0;
	int match_gnss;
	match_gnss = ig_sync_obs(popt, rover_idx, obs, imu, kimu_idx, popt->insopt.back);

	if (match_gnss) {
		rtk->epoch++;
		nobs = inputgnss(epoch_obs, 0, popt, obs, rover_idx, base_idx, 0);
		if (nobs > 4) {
			for (i = 0; i < nobs; i++) {
				if ((satsys(epoch_obs[i].sat, nullptr) & popt->navsys) && popt->exsats[epoch_obs[i].sat - 1] != 1) {
					epoch_obs[n++] = epoch_obs[i];
				}
			}
			if (!rtkpos(rtk, epoch_obs, n, nav)) {
				return 0;
			}
		}
	}
	else {
		return 0;
	}
	return 1;
}

static int wrtsol(rtk_t* rtk, const solopt_t* sopt, stream_t* m) {
	unsigned char buff[1024];
	static int c = 0;
	int n = 0;

	if (rtk->sol.time.time = 0.0) return 0;

	if (m && m->port) {
		n = outsols(buff, &rtk->sol, rtk->rb, sopt, &rtk->opt, nullptr);
		strwrite(m, buff, n);
	}

	return n;
}

static int tcfilter(const prcopt_t* popt, const obs_t* obss, const nav_t* navs, int* krover_idx, int* kbase_idx, rtk_t* ins_rtk)
{
	int nobs, n = 0, stat = 0;
	obsd_t epoch_obs[MAXOBS * 2] = { 0 };
	nobs = inputgnss(epoch_obs, 0, popt, obss, krover_idx, kbase_idx, 0);
	if (nobs) {
		if (!(stat = rtkpos(ins_rtk, epoch_obs, nobs, navs))) {
			trace(2, "%s(%d): tightly coupled update failed\n", time_str(obss->data[*kbase_idx].time, 1), ins_rtk->ins_kf->couple_epoch);
		}
		if (popt->insopt.velopt) {
			doppler(ins_rtk, epoch_obs, nobs, navs, popt);
		}
		return stat;
	}
	else {
		trace(2, "%s%d: no found GNSS obs for tightly coupled\n", time_str(ins_rtk->sol.time, 1), ins_rtk->ins_kf->couple_epoch);
		return 0;
	}
}

extern int couplepos(prcopt_t* popt, const filopt_t* fopt, solopt_t* solopt, stream_t* moni)
{
	///* Global variable definition ----------------------------------------------------*///
	static nav_t navs = { 0 };   /* navigation ephemeris */
	static obs_t obss = { 0 };   /* observation data */
	static rtk_t ins_rtk = { 0 };
	static int kimu_idx = 0, krover_idx = 0, kbase_idx = 0, kpva_idx = 0;

	///* Local variable definition ----------------------------------------------------*///
	imu_t imu;    /* imu data*/
	pva_t pva;    /* pva information for lc slove*/
	od_t od;      /* Odometer data*/
	kf_t inskf;   /* filter struct */
	gtime_t ts = { 0 }, te = { 0 };

	int ins_init = 0, zero_flag = 0, ep_count = 0, vflag_count = 0, match_gnss = 0;
	int back_id[5];
	double dt;
	FILE* fp_rts = nullptr;

	/* Determine the solution mode  */
	int ppk = (popt->mode == PMODE_TC_PPK || popt->mode == PMODE_TC_DGPS) ||
		(popt->mode == PMODE_LC_DGPS || popt->mode == PMODE_LC_PPK) || popt->mode == PMODE_STC_PPK ||
		(popt->insopt.imu_align == INS_ALIGN_GNSS_PPK || popt->insopt.imu_align == INS_ALIGN_GNSS_DGPS);
	int ppp = (popt->mode == PMODE_TC_PPP || popt->insopt.imu_align == INS_ALIGN_GNSS_PPP) || (popt->mode == PMODE_LC_PPP) || popt->mode == PMODE_STC_PPP;
	int tc = (popt->mode >= PMODE_TC_SPP && popt->mode <= PMODE_TC_PPP);
	int stc = popt->mode == PMODE_STC_PPK || popt->mode == PMODE_STC_PPP;
	int lc_obs = popt->mode >= PMODE_LC_SPP && popt->mode <= PMODE_LC_PPP;
	int lc = popt->mode >= PMODE_LC_POS && popt->mode <= PMODE_LC_PPP;
	int mode = ppp ? 2 : ppk ? 1 : 0;

	///* Code start ----------------------------------------------------------------- *///
	/* prepare for debug trace*/
	traceopen("trace.txt");
	tracelevel(solopt->trace);

	/* structure initialization and read configuration */
	imu_init(&imu);
	if (!ins_readf(fopt->imupf, STRFMT_IMU_PROPERTY, &imu, NULL, NULL, &popt->insopt)) {   // IMU config
		return 0;
	}
	popt->insopt.imup = *imu.property;
	if (!ins_readf(fopt->imuf, imu.property->strfmt, &imu, NULL, NULL, &popt->insopt)) {    // IMU data
		return 0;
	}
	imud_t* imuz = (imud_t*)malloc(sizeof(imud_t) * popt->insopt.zvopt.ws), imu_data = { 0 };   // zero vel window

	pva_init(&pva, 0);
	if (popt->mode == PMODE_LC_POS) {
		if (popt->insopt.posfmt == LCPOSFMT_POS) {
			if (!ins_readf(fopt->posf, STRFMT_POS, NULL, &pva, NULL, &popt->insopt)) {
				return 0;
			}
		}
		else if (popt->insopt.posfmt == LCPOSFMT_YGM_PV) {
			if (!ins_readf(fopt->posf, STRFMT_YGM_PV, &imu, &pva, NULL, &popt->insopt)) {
				return 0;
			}
		}
	}

	od_init(&od);
	if (popt->insopt.odo) {
		if (!ins_readf(fopt->strpath[FT_IDX_OD], fopt->strfmt[FT_IDX_OD], NULL, NULL, &od, &popt->insopt)) {
			return 0;
		}
	}

	/* read gnss file */
	if (tc || lc_obs || stc) {
		if (!read_gnss_file(popt, fopt, mode, ppk ? 1 : 0, ts, te, &obss, &navs)) {
			return 0;
		}
	}

#if 0
	save2ygm(&pva, &imu);
#endif

	/* write header to output file */
	if (solopt->outhead && !outhead(fopt->solf, nullptr, 0, popt, solopt)) {
		freeobsnav(&obss, &navs);
		return 0;
	}

	/*open solution file*/
	FILE* fp_sol = fopen(fopt->solf, "a");

	/// * start slove
	while (inputimu(popt->insopt.zvopt.ws, &kimu_idx, &imu, &imu_data, imuz)) {
		if (popt->ts.time != 0.0 && timediff(imu.data[kimu_idx].time, popt->ts) < 0) {
			kimu_idx++;
			continue;
		}
		if (popt->te.time != 0.0 && timediff(imu.data[kimu_idx].time, popt->te) > 0) break;

		/* align ins and initialization */
		if (ins_init == 0) {
			ins_init = ins_align(popt, &obss, &navs, &pva, &imu, &kimu_idx, &krover_idx, &kbase_idx, &kpva_idx, &ins_rtk);
			if (ins_init) {
				popt->insopt.gnss_opt = popt;
				ins_rtk.ins_kf = &inskf;
				//rtkfree(&ins_rtk);
				rtkinit(&ins_rtk, popt, imu.property);
				inskf.last_couple_time = inskf.time = imu.data[kimu_idx].time;

				if (popt->insopt.rts && lc) {
					if (!(fp_rts = fopen(fopt->rts_ins_fw, "wb"))) return 0;
				}
				sol_t gnss_sol = { 0 };
				inskf.dthetap = imu.data[kimu_idx - 1].gyro;
				inskf.dvp = imu.data[kimu_idx - 1].accel;
				solins2sol(&popt->insopt, inskf.insstate, &gnss_sol, 0);
				//                outsol(fp_sol,&gnss_sol,nullptr,solopt,popt,inskf.insstate);

								/* open solution statistics */
				char statfile[MAXSTRPATH] = { '\0' };
				if (solopt->sstat > 0) {
					strcpy(statfile, fopt->solf);
					strcat(statfile, ".stat");
					rtkclosestat();
					rtkopenstat(statfile, solopt->sstat);
				}
			}
			if (!ins_init) {
				kimu_idx++;
				continue;
			}
		}

		/* ins mechanical */
		if (inskf.couple_epoch >= 366000)
			traceins(inskf.sol, 0, &popt->insopt);
		inskf_udstate(&inskf, &imu, &kimu_idx, imu.property, &popt->insopt);
		if (inskf.couple_epoch >= 366000)
			traceins(inskf.sol, 1, &popt->insopt);

		inskf.ins_epoch++;
		double sow = time2gpst(inskf.insstate->time, nullptr);
		ep_count++;

		inskf.imu_obs = &imu.data[kimu_idx];

		if (popt->mode > PMODE_INS_MECH) {
			/* match gnss information */
			match_gnss = ig_sync(popt, &inskf, &kpva_idx, &pva, &krover_idx, &kbase_idx, &obss, &navs, &ins_rtk, &imu, kimu_idx);
			/* ins-gnss coupled */
			if (match_gnss) {
				inskf.couple_epoch++;
				if (lc) {
					lcfilter(popt, kpva_idx, &pva, &imu, kimu_idx, &inskf);
				}
				else if (tc || stc) {
					if (!tcfilter(popt, &obss, &navs, &krover_idx, &kbase_idx, &ins_rtk)) {
						match_gnss = 0;
						/* next epoch */
						if (popt->insopt.ms >= 2 && popt->insopt.ms <= 5) {
							kimu_idx += popt->insopt.ms;
						}
						else {
							kimu_idx += 1;
						}
						//                        continue;
					}
					ins_rtk.ins_kf->last_couple_time = ins_rtk.sol.time;
				}
			}
		}

		/* Non-integrity constraints */
		if (vflag_count > popt->insopt.imup.freq_imu) {
			vflag_count = 0;
			dt = timediff(imuz[1].time, inskf.time);

			if (popt->insopt.zupt || popt->insopt.zaru) {
				zero_flag = detimustatic(&popt->insopt, imuz, inskf.insstate->pos, dt);
			}

			if (zero_flag && popt->insopt.zupt) {
				inskf_ZUPT(&inskf, &popt->insopt);
			}
			if (zero_flag && popt->insopt.zaru) {
				inskf_ZARU(&inskf, &popt->insopt);
			}

			if (popt->insopt.nhc) {
				v3_t wib_b = v3_scalar(1 / dt, imu.data[kimu_idx].gyro);
				inskf_NHC(&inskf, &wib_b, &popt->insopt);
			}
		}

		sol_t gnss_sol = ins_rtk.sol;
		solins2sol(&popt->insopt, inskf.sol, &gnss_sol, 1);

		/* write solution */
		if (solopt->outfrq == -1) {
			if (match_gnss == 1 || ins_rtk.ins_kf->gnss_loss_flag) {
				outsol(fp_sol, &gnss_sol, nullptr, solopt, popt, inskf.sol);
			}
		}
		else {
			if (solopt->outfrq == 0) {
				/*out all solution*/
				outsol(fp_sol, &gnss_sol, nullptr, solopt, popt, inskf.sol);
			}
			else {
				if ((ep_count >= solopt->outfrq) || match_gnss || popt->mode == PMODE_INS_MECH) {
					if (ep_count >= solopt->outfrq) ep_count = 0;
					outsol(fp_sol, &gnss_sol, nullptr, solopt, popt, inskf.sol);
				}
			}
		}

		/* write forward solution for rts*/
		if (lc && popt->insopt.rts) {
			rts_save_forward(fp_rts, inskf.insstate, &popt->insopt);
		}
		/* next epoch */
		if (popt->insopt.ms >= 2 && popt->insopt.ms <= 5) {
			kimu_idx += popt->insopt.ms;
		}
		else {
			kimu_idx += 1;
		}
		ins_rtk.sol.ns = ins_rtk.sol.stat = 0;
		ins_rtk.sol.dop[1] = 0; ins_rtk.sol.adop = 0;
		match_gnss = 0;
	}

	if (popt->modear >= ARMODE_OFF && (popt->mode > PMODE_INS_MECH && popt->mode != PMODE_LC_POS && popt->mode != PMODE_LC_SPP &&
		popt->mode != PMODE_TC_SPP)) {
		trace(2, "Total coupled epoch: %d, Cannot fix epoch: %d(%4.1f%s), Fix epoch: %d(%4.1f%s), Float epoch: %d(%4.2f%s), Abnormal epoch: %d\n",
			ins_rtk.ins_kf->couple_epoch,
			ins_rtk.cannotfix, (ins_rtk.cannotfix + 1E-09) / ins_rtk.ins_kf->couple_epoch * 100, "%",
			ins_rtk.fix_epoch, (ins_rtk.fix_epoch / (ins_rtk.ins_kf->couple_epoch + 1E-9)) * 100.0, "%",
			ins_rtk.ins_kf->couple_epoch - ins_rtk.fix_epoch, (ins_rtk.ins_kf->couple_epoch - ins_rtk.fix_epoch) / (ins_rtk.ins_kf->couple_epoch + 1E-09) * 100.0, "%",
			ins_rtk.abnormal_epoch);
	}

	if (lc && popt->insopt.rts) {  /// rts
		fclose(fp_rts);
		popt->insopt.back = 1;
		back_id[0] = kimu_idx; back_id[1] = kpva_idx; back_id[2] = krover_idx; back_id[3] = kbase_idx;
		rts_ins_gnss(popt, fopt, solopt, &ins_rtk, &pva, &imu, back_id, &obss, &navs);
	}

	///* Free memory ------------------------------------------------------------- */
	fclose(fp_sol);
	rtkfree(&ins_rtk);
	imu_free(&imu);
	pva_free(&pva);
	od_free(&od);
	free(imuz);
	rtkclosestat();
	traceclose();
	return 1;
}