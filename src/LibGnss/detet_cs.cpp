#include "rtklib.h"
//#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
//#define NP(opt)     ((opt)->dynamics==0?3:9)
//#define NI(opt)     ((opt)->ionoopt!=IONOOPT_UC?0:MAXSAT)
//#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
//#define NL(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
//#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
//#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
//#define NX(opt)     (NR(opt)+NB(opt))
//
///* state variable index */
//#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
//#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
//#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
//#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

#define NF(opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf))
#define NF_SYS(s,opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf_sys[s]))
#define NP(opt)      ((opt)->dynamics?9:3)
#define NC(opt)      ((opt)->sdopt?0:((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS))
#define NRDCB(opt)   (((opt)->ionoopt==IONOOPT_UC_CONS&&(opt)->nf>=2)?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
#define NIFCB(opt)   ((opt)->nf>=3?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
#define NT(opt)      ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)      ((opt)->ionoopt==IONOOPT_UC||(opt)->ionoopt==IONOOPT_UC_CONS?MAXSAT:0)
#define NR(opt)      (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt)+NT(opt)+NI(opt))
#define NB(opt)      (NF(opt)*MAXSAT)
#define NX(opt)      (NR(opt)+NB(opt))
#define IC(s,opt)    (NP(opt)+(s))
#define IRDCB(s,opt) (NP(opt)+NC(opt)+(s))
#define IIFCB(s,opt) (NP(opt)+NC(opt)+NRDCB(opt)+(s))
#define IT(opt)      (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt))
#define II(s,opt)    (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt)+NT(opt)+(s)-1)
#define IB(s,f,opt)  (NR(opt)+MAXSAT*(f)+(s)-1)

extern int myRound(const double dNum)
{
	int iNum;
	if (dNum >= 0)
		iNum = (int)(dNum + 0.5);
	else
		iNum = (int)(dNum - 0.5);

	return iNum;
}
extern bool calCsThres(rtk_t* rtk, const double sample)
{
	bool b = false;
	// jxy
	if (sample > 0.0)
	{
		if (true == rtk->opt_ex.bUsed_gfCs && fabs(rtk->opt_ex.csThresGF) < 0.01)
		{
			if (sample <= 1.0)         rtk->opt_ex.csThresGF = 0.05;
			else if (sample <= 20.0)   rtk->opt_ex.csThresGF = (0.10) / (20.0) * sample + 0.05;
			else if (sample <= 60.0)   rtk->opt_ex.csThresGF = 0.15;
			else if (sample <= 100.0)  rtk->opt_ex.csThresGF = 0.25;
			else					   rtk->opt_ex.csThresGF = 0.35;
			b = true;
		}
		if (true == rtk->opt_ex.bUsed_mwCs && fabs(rtk->opt_ex.csThresMW) < 0.01)
		{
			if (sample <= 1.0)         rtk->opt_ex.csThresMW = 3.0;
			else if (sample <= 20.0)   rtk->opt_ex.csThresMW = (2.5) / (20.0 - 0.0) * sample + 2.5;
			else if (sample <= 60.0)   rtk->opt_ex.csThresMW = 5.0;
			else                     rtk->opt_ex.csThresMW = 7.5;
			b = true;
		}
		return b;
	}
	else
	{
		printf("sample<=0.0! wait.");
		//getchar();
		rtk->opt_ex.csThresGF = 0.10;
		rtk->opt_ex.csThresMW = 4.5;
		b = false;
	}

	/*
	if ( sample > 0.0 )
	{
		if ( true ==  rtk->opt_ex.bUsed_gfCs && fabs( rtk->prcopt_ex.csThresGF) < 0.01 )
		{
			if (sample <= 1.0)         rtk->opt_ex.csThresGF = 0.05;
			else if (sample <= 20.0)   rtk->opt_ex.csThresGF = (0.10) / (20.0) * sample + 0.05;
			else if (sample <= 60.0)   rtk->opt_ex.csThresGF = 0.15;
			else if (sample <= 100.0)  rtk->opt_ex.csThresGF = 0.25;
			else                     rtk->opt_ex.csThresGF = 0.35;

			b = true;
		}

		if ( true ==  rtk->opt_ex.bUsed_mwCs && fabs( rtk->opt_ex.csThresMW) < 0.01 )
		{
			if (sample <= 1.0)         rtk->opt_ex.csThresMW = 3.0;
			else if (sample <= 20.0)   rtk->opt_ex.csThresMW = (2.5) / (20.0 - 0.0) * sample + 2.5;
			else if (sample <= 60.0)   rtk->opt_ex.csThresMW = 5.0;
			else                     rtk->opt_ex.csThresMW = 7.5;

			b = true;
		}

		return b;
	}
	else
	{
		printf("sample<=0.0! wait.");
		getchar();
		 rtk->opt_ex.csThresGF = 0.10;
		 rtk->opt_ex.csThresMW = 4.5;
		b = false;
	}
	*/
	return b;
}

// jxy 20180105
// type 1 : mw
//      2 : gf
extern double get_slip_thres(const rtk_t* rtk, const int sat, double dt, double elev, int type, double* thres_dd)
{
	//以下阈值适用于30s采样率的情况
	double thres = 0;
	const double rad_20(20 * D2R);
	const double rad_15(15 * D2R);
	const double factor(1.0 / 15.0);
	if (type == 1)
	{
		if (dt <= 1.0)        thres = 3.0;
		else if (dt <= 20.0)  thres = (2.5) / (20.0 - 0.0) * dt + 2.5;
		else if (dt <= 60.0)  thres = 5.0;
		else                  thres = 7.5;
	}
	else if (type == 2)
	{
		if (dt <= 1.0)        thres = 0.05;
		else if (dt <= 20.0)  thres = (0.10) / (20.0) * dt + 0.05;
		else if (dt <= 60.0)  thres = 0.15;
		else if (dt <= 100.0) thres = 0.25;
		else                  thres = 0.35;
	}
	else
	{
		return 0.0;
	}

	if (elev < rtk->opt.elmin)
		elev = rtk->opt.elmin;

	if (type == 1)
	{
		if (elev < rad_20)
			thres = -thres * 0.1 * elev * R2D + 3 * thres;
	}
	else if (type == 2)
	{
		if (elev < rad_15)
		{
			thres = -thres * factor * elev * R2D + 2 * thres;
			(*thres_dd) = -(*thres_dd) * elev * R2D * factor + 2 * (*thres_dd);
		}
		int prn = 0; satsys(sat, &prn);
		if (SYS_CMP == satsys(sat, NULL) && prn <= 5)
		{
			thres = 0.02;
			*thres_dd = 0.025;
		}
	}
	return thres;
}

// jxy 2017-12-30
/* detect cycle slip by widelane jump -----------------------------*/
// thres computed from time difference

// 获取最近的gf
extern void getGF_Nearby(rtk_t* rtk, const int sat, double* gf_1, double* gf_2, gtime_t* gt_1, gtime_t* gt_2, bool* b_1, bool* b_2)
{
	int i;
	double dt;
	gtime_t gt, gt0 = { 0 };

	*gf_1 = *gf_2 = 0.0;
	*gt_1 = *gt_2 = gt0;

	//vector<ion_base_t>::iterator it;
	kf_ion_t* pki = &rtk->ssat_ex[sat - 1].kf_ion;

	// jxy 20171230 kf_ion
	//int isize = (int)(pki->vec.size());
	int isize = pki->n;

	if (isize <= 0) return;

	double elminDeg = rtk->opt_ex.elMin * R2D;

	for (i = isize - 1; i >= 0; i--)
	{
		if (pki->vec[i].elev_deg < elminDeg) continue;

		dt = pki->vec[i].gf;
		gt = pki->vec[i].gt;
		bool j = OBSSTATE_CS == pki->vec[i].iCs;

		if (0.0 != dt)
		{
			if (0.0 == *gf_1)
			{
				*gf_1 = dt;
				*gt_1 = gt;
				*b_1 = j;
			}
			else if (0.0 == *gf_2)
			{
				*gf_2 = dt;
				*gt_2 = gt;
				*b_2 = j;
			}
			else break;
		}
	}
}

//Panda 20151107
//考虑最近观测值周跳，无第二近观测值的情况，只用一次差判断；其余情况可以用二次差考虑
extern int detslp_ddgf0(rtk_t* rtk, const int sat, const double g0, const double g_1, const double g_2, const gtime_t gt0,
	const gtime_t gt_1, const gtime_t gt_2, const bool bcs_1, const int nout, const double thres, const double thres_dd, double* dg01, double* dg12)
{
	bool bSlip = false;

	double dg01_tmp = g0 - g_1;
	double dg12_tmp = g_1 - g_2;
	double dd02_tmp = dg01_tmp - dg12_tmp;

	double dg01fabs = fabs(dg01_tmp);
	double dd02fabs = fabs(dd02_tmp);

	double dt01 = fabs(timediff(gt0, gt_1));
	double dt12 = fabs(timediff(gt_1, gt_2));
	int    it01 = myRound(dt01 / rtk->ppp_glo.sample);
	int    it12 = myRound(dt12 / rtk->ppp_glo.sample);

	double gflimit = thres;

	if (it01 >= 3 && dt01 >= 90.0)
	{
		gflimit = MIN(thres * it01, 0.35);
	}
#if 1	// jxy
	//if (0.0 == g_2 || bcs_1 || it12 >= nout) {	//最近观测值周跳，无第二近观测值的情况，只用一次差判断
	//	dg12_tmp = dd02_tmp = dd02fabs = 0.0;
	if (dg01fabs < gflimit) bSlip = false;
	else                    bSlip = true;
	//}
#else
	if (0.0 == g_2 || bcs_1 || it12 >= nout)  	//最近观测值周跳，无第二近观测值的情况，只用一次差判断
	{
		dg12_tmp = dd02_tmp = dd02fabs = 0.0;
		if (dg01fabs < gflimit) bSlip = false;
		else                    bSlip = true;
	}
	else
	{
		if (dd02fabs > thres_dd * MIN(it01, it12))
		{
			//if ( dg01fabs>gflimit )
			//	bSlip=true;
			if (dg01fabs > 0.03)
			{
				bSlip = true;
			}
			else
			{
				sprintf(msgbuf[rtk->is], "%s gf_dd=%.3f  _sd=%.3f\n", rtk->ppp_glo.sFlag[sat - 1].id, dd02fabs, dg01fabs);
				outDebug(rtk, true, rtk->pf.pfdebug);
			}
		}

		else if (dg01fabs > gflimit + 0.15) bSlip = true;
		else
		{
			if (dd02fabs > thres_dd) bSlip = true;
			else                   bSlip = false;
		}
	}
#endif
	if (bSlip)
		rtk->ssat_ex[sat - 1].iCycle_gf = OBSSTATE_CS;
	else if (dg01fabs > 0.05)
		rtk->ssat_ex[sat - 1].iCycle_gf = OBSSTATE_CS_PB;

	*dg01 = dg01_tmp;
	*dg12 = dg12_tmp;

	return bSlip ? 1 : 0;
}

//static bool solEWLAmb(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const prcopt_t *opt, double &ambEWL)
//{
//	double lam[6] = { 0 };
//	double L2, L5, P2, P5, P1_C1, P2_C2;
//	int i = 0, j = 1, i1, i2, i5, sys,prn;
//	sys = satsys(obs->sat, &prn);
//	int sys_idx = satsysidx(obs->sat);
//	for (i = 0; i < NFREQ + NEXOBS;i++)lam[i] = sat2freq(obs->sat, obs->code[i], nav);
//	ambEWL = 0.0;
//
//	//sys=satsys(obs->sat,NULL);
//	//sys = satsys(obs->sat, NULL);
//
//	i1 = 0;
//	i2 = 1;
//	i5 = 2;
//
//	/* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
//	if (NFREQ >= 3 && (sys & (SYS_GAL | SYS_SBS)))
//	{
//		i1 = 0;
//		i2 = 1;
//		i5 = 4;
//	}
//
//	if (NFREQ < 2 || lam[i2] == 0.0 || lam[i5] == 0.0) return 0;
//
//	L2 = obs->L[i2] * lam[i2];        /* cycle -> m */
//	L5 = obs->L[i5] * lam[i5];
//	P2 = obs->P[i2];
//	P5 = obs->P[i5];
//
//	P1_C1 = nav->cbias[obs->sat - 1][1];
//	P2_C2 = nav->cbias[obs->sat - 1][2];
//
//	if (L2 == 0.0 || L5 == 0.0 || P5 == 0.0 || P2 == 0.0) return 0;
//
//	if (sys & (SYS_CMP))
//	{
//		ambEWL = (FREQ2_CMP * L2 - FREQ3_CMP * L5) / (FREQ2_CMP - FREQ3_CMP) - (FREQ3_CMP * P5 + FREQ2_CMP * P2) / (FREQ3_CMP + FREQ2_CMP);
//		ambEWL = ambEWL / CLIGHT * (FREQ2_CMP - FREQ3_CMP);
//	}
//	else
//	{
//		ambEWL = (FREQ2 * L2 - FREQ5 * L5) / (FREQ2 - FREQ5) - (FREQ2 * P2 + FREQ5 * P5) / (FREQ2 + FREQ5);
//		ambEWL = ambEWL / CLIGHT * (FREQ2 - FREQ5);
//	}
//
//	return 1;
//}
//
//extern void detslp_ewl(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
//{
//	int i, j, sat,prn;
//	const prcopt_t *opt = &rtk->opt;
//	double ewl0, ewl1, elev, el, thres = 0.0;
//
//	for (i = 0; i < n && i < MAXOBS; i++)
//	{
//		sat = obs[i].sat;
//		elev = rtk->ssat[sat - 1].azel[1] * R2D;
//		satsys(sat, &prn);
//		ewl0 = rtk->ssat_ex[sat - 1].eWl;
//		solEWLAmb(rtk, obs + i, nav, opt, ewl1);
//		rtk->ssat_ex[sat - 1].eWl = ewl1;
//
//		if (ewl0 * ewl1 == 0.0) continue;
//
//		el = elev;
//
//		if (el < rtk->opt.elmin * R2D)
//			el = rtk->opt.elmin * R2D;
//
//
//		if (el < 15.0) thres = 0.85;
//		else		 thres = 0.75;
//
//		if (fabs(ewl1 - ewl0) > MIN(thres, 1.5))
//		{
//			for (j = 0; j < rtk->opt.nf; j++)
//				rtk->ssat[sat - 1].slip[j] |= 1;
//
//			if (elev >= rtk->opt.elmin * R2D)
//			{
//				if (1)
//					trace(3, "prn%d EWL slip.  diff=%10.3f/%-3d/%d  .  elev=%4.1f\n",  prn,
//						ewl1 - ewl0, rtk->ssat[sat - 1].lock[0], rtk->ssat[sat - 1].outc[0], elev);
//				//printf("%s %s EWL周跳.  diff=%10.3f/%-3d/%d  .  elev=%4.1f\n", rtk->ppp_glo.chTime_s, rtk->ppp_glo.sFlag[sat - 1].id,
//				//	ewl1 - ewl0, rtk->ssat[sat - 1].lock[0], rtk->ssat[sat - 1].outc[0], elev);
//				//outDebug(rtk, true, rtk->pf.pfdebug);
//			}
//		}
//	}
//}
void statRecord1(rtk_t* rtk, const obsd_t* obs, const int nobs, gtime_t currt_time)
{
	int i, sat;
	int nf = NF(&rtk->opt);
	//printf(" %02d:%02d:%02d ", (int)rtk->ppp_glo.ctNow[3], (int)rtk->ppp_glo.ctNow[4], (int)rtk->ppp_glo.ctNow[5]);
	for (i = 0; i < nobs; i++)
	{
		sat = obs[i].sat;

		//rtk->ppp_glo.obsLast[sat - 1][0] = obs[i].P[0];
		//rtk->ppp_glo.obsLast[sat - 1][1] = obs[i].P[1];
		//rtk->ppp_glo.obsLast[sat - 1][2] = obs[i].P[2];
		//rtk->ppp_glo.obsLast[sat - 1][3] = obs[i].L[0];
		//rtk->ppp_glo.obsLast[sat - 1][4] = obs[i].L[1];
		//rtk->ppp_glo.obsLast[sat - 1][5] = obs[i].L[2];

#if 0
		// for re-filter jxy 20200406
		continue;
#endif

		kf_ion_t* pki = &rtk->ssat_ex[sat - 1].kf_ion;
		/*
		int isize = (int)(pki->vec.size());
		*/
		// jxy 20171230
		int isize = pki->n;

		int k1 = 0, k2 = 1;
		//        if ( SYS_GAL == rtk->ppp_glo.sFlag[sat - 1].sys ) k2 = 2;
		//#ifdef RTPPP_BDS23_B1B3
		//		if (SYS_CMP == rtk->ppp_glo.sFlag[sat - 1].sys) k2 = 2;
		//#endif

		ion_base_t ib;
		ib.gt = obs[i].time;
		ib.iCs = rtk->ssat_ex[sat - 1].iCycle;
		ib.iCs_gf = rtk->ssat_ex[sat - 1].iCycle_gf;
		ib.elev_deg = rtk->ssat[sat - 1].azel[1] * R2D;
		ib.gf = rtk->ssat[sat - 1].gf[0];
		ib.gf_pr = obs[i].P[k1] - obs[i].P[k2];
		k1 = II(sat, &rtk->opt);
		ib.ion_val = rtk->x[k1];
		ib.ion_var = rtk->P[k1 + k1 * rtk->nx];

		if (0 == isize)
		{
			init_kf_ion(&ib, pki);
		}
		else
		{
#if 0
			// lipan
			vector<ion_base_t>::iterator it;
			for (it = pki->vec.begin(); it != pki->vec.end(); )
			{
				double dt = fabs(timediff(rtk->ppp_glo.tNow, it->gt));
				if (fabs(dt) > 600)
					it = pki->vec.erase(it);
				else
					it++;
			}
			pki->vec.push_back(ib);
#else
			// jxy 20171230 kf_ion
			// 对存储的历史电离层处理，删掉超时的
			// 如果采用反向滤波的话，pki->vec[pki->n - 1] = ib;会出问题
			int b_del = 1;
			while (b_del)
			{
				b_del = 0;
				for (int k = 0; k < pki->n; k++)
				{
					double dt = fabs(timediff(obs[i].time, pki->vec[k].gt));
					if (fabs(dt) > 600)
					{
						// 如果只有一个元素的话，就不要平移了
						if (pki->n > 1)
						{
							for (int x = k; x < pki->n - 1; x++)
								pki->vec[x] = pki->vec[x + 1];
						}
						//memmove(pki->vec + k, pki->vec + k + 1, pki->n - k - 1);
						(pki->n)--;
						b_del = 1;
						break;
					}
				}
			}
			(pki->n)++;
			pki->vec[pki->n - 1] = ib;
			//printf(" %2d-%2d",sat, pki->n);
		}
	}
	//printf("\n");
#endif
}