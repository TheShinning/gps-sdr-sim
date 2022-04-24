/*-----------------------------------------------------------------------------
* ins-gnss-lc.cc : loose combination code
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/22.
*----------------------------------------------------------------------------*/
#include "rtklib.h"

#define IV  3
#define IP  0

/* Simplfiy the filter usage */
#define KF_FILTER(kf, dz, Qdz)                                               \
    filter(kf->x, kf->P, kf->H, (dz), (Qdz), kf->nx, kf->ny,0, 0, NULL,0)

/* set kf->ny and set kf->H to zero */
#define KF_HINIT(kf, _ny)                                                    \
{                                                                            \
    kf->ny = _ny;                                                            \
    for(unsigned int i = 0; i < kf->ny*kf->nx; ++i)                          \
        kf->H[i] = 0.0;                                                      \
}

/* loosely coupled update by position ---------------------------------------*/
extern int inskf_udmeasr(kf_t* inskf, const v3_t* reg_e, const m3_t* Qreg_e, const v3_t* lever_arm_gps, const insopt_t* opt) {
	KF_HINIT(inskf, 3);
	m3_t Cbe = inskf->insstate->dcm;

	//    fprintf(stdout,"%d\n",inskf->couple_epoch);
		/* reb_e measments update */
	v3_t reb_e_gnss = v3_add(*reg_e, m3_mul_v3(Cbe, *lever_arm_gps));
	v3_t dz_r = v3_del(inskf->insstate->pos, reb_e_gnss); // v

	m3_t m3H;
	jacobi_lc_meas_p2a_ECEF(&m3H, &Cbe, lever_arm_gps);
	m3H = m3_T(m3H);
	m3_paste((inskf->H + xiA(opt)), inskf->nx, &m3H);
	inskf->H[xiP(opt)] = 1.0;
	inskf->H[xiP(opt) + 1 + 1 * inskf->nx] = 1.0;
	inskf->H[xiP(opt) + 2 + 2 * inskf->nx] = 1.0;
	//    matprint(1,inskf->H,15,3,15,6);
	//    matprint(0,inskf->P,15,15,15,6);

		/* innovtaion check */
	double ndz_r[3];
	inskf_norm_innov(inskf, (const double*)&dz_r, (const double*)Qreg_e, ndz_r);
	for (unsigned int i = 0; i < inskf->ny; ++i) {
		if (fabs(ndz_r[i]) > 5.0) {
			const double* pdz_r = (const double*)&dz_r;
			trace(3, "%s too small noise(%f m) or too large innovations(%f m)\n", time_str(inskf->time, 3),
				pdz_r[i] / ndz_r[i], pdz_r[i]);
		}
	}

	//    matprint(1,inskf->x,15,1,20,15);
	//    matprint(0,inskf->F,15,15,20,15);
	//    matprint(0,inskf->P,15,15,20,15);

		/* note: Qreg_e is sysmmetry matrix */
	KF_FILTER(inskf, (const double*)&dz_r, (const double*)Qreg_e);
	//    matprint(1,inskf->x,15,1,20,15);

	return 0;
}

/* loosely coupled update by velocity ------------------------------------------ */
extern int inskf_udmeasv(kf_t* inskf, const v3_t* veg_e, const m3_t* Qveg_e,
	const v3_t* lever_arm, const v3_t* wib_b, const insopt_t* opt) {
	KF_HINIT(inskf, 3);
	m3_t Cbe = inskf->insstate->dcm;

	v3_t v1 = m3_mul_v3(Cbe, v3_cross(*wib_b, *lever_arm));
	v3_t wie_e = { 0.0, 0.0, wgs84.wie };
	v3_t v2 = v3_cross(wie_e, m3_mul_v3(Cbe, *lever_arm));
	/*todo:check*/
	v3_t veb_e_gnss = v3_del(v3_del(*veg_e, v1), v2);
	v3_t dz_v = v3_del(veb_e_gnss, inskf->insstate->vel);

	m3_t m3H;
	jacobi_lc_meas_v2a_ECEF(&m3H, &Cbe, wib_b, lever_arm);
	m3H = m3_T(m3H);
	m3_paste((inskf->H + xiA(opt)), inskf->nx, &m3H);
	if (opt->est_bg) {
		jacobi_lc_meas_v2bg_ECEF(&m3H, &Cbe, lever_arm);
		m3H = m3_T(m3H);
		m3_paste((inskf->H + xiBg(opt)), inskf->nx, &m3H);
	}
	inskf->H[xiV(opt)] = -1.0;
	inskf->H[xiV(opt) + 1 + 1 * inskf->nx] = -1.0;
	inskf->H[xiV(opt) + 2 + 2 * inskf->nx] = -1.0;

	double ndz_v[3];
	inskf_norm_innov(inskf, (const double*)&dz_v, (const double*)Qveg_e, ndz_v);
	for (unsigned int i = 0; i < inskf->ny; ++i) {
		if (fabs(ndz_v[i]) > 5.0) {
			const double* pdz_v = (const double*)&dz_v;
			trace(3, "too small noise(%f m) or too large innovations(%f m)\n", pdz_v[i] / ndz_v[i], pdz_v[i]);
		}
	}

	KF_FILTER(inskf, (const double*)&dz_v, (const double*)Qveg_e);
	//    matprint(inskf->x, 1, inskf->nx, 8, 5);
	//    fflush(stdout);

	return 0;
}

static v3_t vn2pos(solins_t* inssol, v3_t vr, double dt)
{
	v3_t dpos = { 0 };
	dpos.x = dt * vr.y / inssol->eth.RMh;
	dpos.y = dt * vr.x / inssol->eth.clRNh;
	dpos.z = dt * vr.z;

	return dpos;
}

static int lc_meas_update_pv(const insopt_t* opt, kf_t* inskf, v3_t* arm, v3_t* re, m3_t* Qre, v3_t* ve, m3_t* Qve, int v_aiding) {
	int i, j, tra = 0, ny = 6;
	double zk[6] = { 0 }, * R;
	auto* p = (prcopt_t*)opt->gnss_opt;

	if (!v_aiding) {
		ny -= 3;
	}

	v3_t gr = { 0 }, gv = { 0 };
	v3_t delay_pos = inskf->insstate->delay_pos;
	KF_HINIT(inskf, ny);
	R = zeros(ny, ny);

	ins2gnss(inskf->insstate, arm, nullptr, gr.v, gv.v, opt->mech_coord);

	if (opt->est_igdt) {
		delay_pos = vn2pos(inskf->insstate, gv, inskf->insstate->t_delay);
		inskf->insstate->delay_pos = delay_pos;
		gr = v3_add(gr, delay_pos);
	}

	asi_blk_mat(R, ny, ny, Qre->v, 3, 3, IP, IP);

	for (i = IP; i < IP + 3; i++) {
		zk[i] = gr.v[i - IP] - re->v[i - IP];
	}
	inskf->H[xiP(opt) + 0 + IP * inskf->nx] = 1.0;
	inskf->H[xiP(opt) + 1 + (IP + 1) * inskf->nx] = 1.0;
	inskf->H[xiP(opt) + 2 + (IP + 2) * inskf->nx] = 1.0;

	if (opt->est_armgps) {
		m3_t dpdl = m3_scalar(-1.0, inskf->insstate->MpvCnb);
		asi_blk_mat(inskf->H, inskf->nx, ny, dpdl.v, 3, 3, xiArmGps(opt), IP);
	}
	if (opt->est_igdt) {
		v3_t dpdt = v3_scalar(-1.0, inskf->insstate->Mpvvn);
		asi_blk_mat(inskf->H, inskf->nx, ny, dpdt.v, 1, 3, xiDt(opt), IP);
	}

	if (v_aiding) {
		gv = v3_add(gv, v3_scalar(inskf->insstate->t_delay, inskf->insstate->an));
		for (i = IV; i < IV + 3; i++) {
			zk[i] = gv.v[i - IV] - ve->v[i - IV];
		}
		asi_blk_mat(R, ny, ny, Qve->v, 3, 3, IV, IV);

		inskf->H[xiV(opt) + 0 + IV * inskf->nx] = 1.0;
		inskf->H[xiV(opt) + 1 + (IV + 1) * inskf->nx] = 1.0;
		inskf->H[xiV(opt) + 2 + (IV + 2) * inskf->nx] = 1.0;
		if (opt->est_armgps) {
			m3_t dvdl = m3_scalar(-1.0, inskf->insstate->CW);
			asi_blk_mat(inskf->H, inskf->nx, ny, dvdl.v, 3, 3, xiArmGps(opt), IV);
		}
		if (opt->est_igdt) {
			v3_t dvdt = v3_scalar(-1.0, inskf->insstate->an);
			asi_blk_mat(inskf->H, inskf->nx, ny, dvdt.v, 1, 3, xiDt(opt), IV);
		}
	}

#if 0
	matprint(1, inskf->H, inskf->nx, inskf->ny, 15, 6);
	matprint(0, R, inskf->ny, inskf->ny, 15, 6);
	matprint(1, zk, inskf->ny, 1, 15, 6);
#endif
	if (inskf->couple_epoch >= 1) {
		tra = 0;
	}
	filter(inskf->x, inskf->P, inskf->H, zk, R, inskf->nx, ny, 0, p->kfopt, nullptr, tra);

	free(R);
	return 1;
}

extern int lcfilter(const prcopt_t* popt, int kpva_idx, const pva_t* pva, const imu_t* imu, int kimu_idx, kf_t* inskf) {
	int back = popt->insopt.back ? 1 : 0;
	int exist_pos = 0, exist_vel = 0;

	v3_t gnss_pos = { 0, 0, 0 }, gnss_vel = { 0, 0, 0 };
	m3_t R_pos = { 0 }, R_vel = { 0 };

	if (pva->status[kpva_idx] != SOLQ_NONE) {
		if (popt->mode == PMODE_LC_POS || (popt->mode >= PMODE_LC_SPP && popt->mode <= PMODE_LC_PPP)) {
			gnss_pos = pva->pos[kpva_idx];
			R_pos = pva->Qpos[kpva_idx];
			if (R_pos.m11 * R_pos.m22 * R_pos.m33 == 0.0) {
				R_pos.m11 = R_pos.m22 = R_pos.m33 == 0.01;
			}
			if (popt->insopt.mech_coord == INSMECH_LLH) {
				v3_t v = { 0.1,0.1,0.1 };
				R_pos = v3_diag(v3_pow(v, 2));
				R_pos.m11 /= SQR(wgs84.R0); R_pos.m22 /= SQR(wgs84.R0);
			}
			gnss_vel = pva->vel[kpva_idx];

			if (v3_norm(gnss_vel) != 0.0) {
				R_vel = pva->Qvel[kpva_idx];
				exist_vel = 1;
				if (R_vel.m11 * R_vel.m22 * R_vel.m33 == 0.0) {
					v3_t v = { 0.1,0.1,0.1 };
					R_vel = v3_diag(v3_pow(v, 2));
				}
			}
			else exist_vel = 0;
		}

		if (popt->insopt.mech_coord == INSMECH_LLH) {
			m3_t A;
			asymmetric_mat(&inskf->insstate->web, &A);
			inskf->insstate->CW = m3_mul(m3_T(inskf->insstate->dcm), A);
			inskf->insstate->MpvCnb = m3_mul(inskf->insstate->eth.Mpv, m3_T(inskf->insstate->dcm));
			inskf->insstate->Mpvvn = m3_mul_v3(inskf->insstate->eth.Mpv, inskf->insstate->vel);
			lc_meas_update_pv(&popt->insopt, inskf, popt->insopt.est_armgps ? &inskf->insstate->arm_gps : &imu->property->lever_arm_gps,
				&gnss_pos, &R_pos, &gnss_vel, &R_vel, exist_vel);
		}
		else if (popt->insopt.mech_coord == INSMECH_ECEF) {
			inskf_udmeasr(inskf, &gnss_pos, &R_pos, &popt->insopt.imup.lever_arm_gps, &popt->insopt);
			/* loosely coupled update by velocity */
			if (exist_vel) {
				v3_t wib_b = v3_scalar(imu->property->freq_imu, imu->data[kimu_idx].gyro);
				inskf_udmeasv(inskf, &gnss_vel, &R_vel, &popt->insopt.imup.lever_arm_gps, &wib_b, &popt->insopt);
			}
			if (back) {
				for (int i = 0; i < 3; i++) inskf->x[xiV(&popt->insopt) + i] = -inskf->x[xiV(&popt->insopt) + i];
				for (int i = 0; i < 3; i++) inskf->x[xiBg(&popt->insopt) + i] = -inskf->x[xiBg(&popt->insopt) + i];
			}
		}
		inskf_feedback(pva->time[kpva_idx], SOLQ_IGLC, &popt->insopt, inskf->x, inskf->insstate);
		inskf->insstate->ns = 0;
		inskf->sol = inskf->insstate;
	}
	return 1;
}