//
// Created by chenc on 2021/3/22.
//

#include "rtklib.h"

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
//#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
//#define NP(opt)     ((opt)->dynamics==0?3:9)
//#define NI(opt)     ((opt)->ionoopt!=IONOOPT_UC?0:MAXSAT)
//#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
//#define NL(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
//#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
//#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
//#define NX(opt)     (NR(opt)+NB(opt))

//#define NF(opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf))
//#define NF_SYS(s,opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf_sys[s]))
//#define NP(opt)      ((opt)->dynamics?9:3)
//#define NC(opt)      ((opt)->sdopt?0:((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS))
//#define NRDCB(opt)   (((opt)->ionoopt==IONOOPT_UC_CONS&&(opt)->nf>=2)?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
//#define NIFCB(opt)   ((opt)->nf>=3?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
//#define NT(opt)      ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
//#define NI(opt)      ((opt)->ionoopt==IONOOPT_UC||(opt)->ionoopt==IONOOPT_UC_CONS?MAXSAT:0)
//#define NR(opt)      (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt)+NT(opt)+NI(opt))
//#define NB(opt)      (NF(opt)*MAXSAT)
//#define NX(opt)      (NR(opt)+NB(opt))
///* state variable index */
////#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
////#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
////#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
////#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */
//
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

#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */

#define  VAR_AMB_INIT  10000.0

#define NUM_SYS NSYS+1

#define SWAP_I(x,y) do {int    _tmp=x; x=y; y=_tmp;} while (0)
#define SWAP_D(x,y) do {double _tmp=x; x=y; y=_tmp;} while (0)
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t* rtk, double xi, double* P, double var, int i)
{
	int j;

	rtk->x[i] = xi;
	for (j = 0; j < rtk->nx; j++) {
		P[i + j * rtk->nx] = P[j + i * rtk->nx] = i == j ? var : 0.0;
	}
}

static void res_class(res_t* res, int pri)
{
	int i, j = 0, k = 0, type;
	if (pri) {
		for (i = 0; i < res->nv; i++) {
			type = (res->vflag[i] >> 4) & 0xF;
			if (type == 1) {
				res->npr++;
			}
			else if (type == 0) {
				res->ncp++;
			}
		}
		res->pri_pr = mat(res->npr, 1);
		res->pri_cp = mat(res->ncp, 1);
		res->post_pr = mat(res->npr, 1);
		res->post_cp = mat(res->ncp, 1);
		res->pr_idx = imat(res->npr, 1);
		res->cp_idx = imat(res->ncp, 1);
		res->norm_pr = mat(res->npr, 1);
		res->norm_cp = mat(res->ncp, 1);
	}

	for (i = 0; i < res->nv; i++) {
		type = (res->vflag[i] >> 4) & 0xF;
		if (type == 1) {
			if (pri) {
				res->pr_idx[j] = i;
				res->pri_pr[j++] = fabs(res->pri_v[i]);
			}
			else {
				res->norm_pr[j] = res->post_v[i] / (res->sigma0 * SQRT(res->R[i + i * res->nv]));
				res->post_pr[j++] = res->post_v[i];
			}
		}
		else if (type == 0) {
			if (pri) {
				res->cp_idx[k] = i;
				res->pri_cp[k++] = fabs(res->pri_v[i]);
			}
			else {
				res->norm_cp[k] = res->post_v[i] / (res->sigma0 * SQRT(res->R[i + i * res->nv]));
				res->post_cp[k++] = res->post_v[i];
			}
		}
	}
}

extern void init_prires(const double* v, const int* vflag, int nv, res_t* res)
{
	res->pri_v = mat(nv, 1);
	res->vflag = imat(nv, 1);
	matcpy(res->pri_v, v, nv, 1);
	for (int i = 0; i < nv; i++) {
		res->vflag[i] = vflag[i];
	}
	res->nv = nv;
	res_class(res, 1);
}

extern void init_postres(rtk_t* rtk, const double* post_v, res_t* res, const double* R)
{
	res->post_v = mat(res->nv, 1);
	res->R = mat(res->nv, res->nv);
	matcpy(res->post_v, post_v, res->nv, 1);
	matcpy(res->R, R, res->nv, res->nv);

	res_class(res, 0);

	int sat, frq;
	for (int i = 0; i < res->ncp; i++) {
		sat = (res->vflag[res->cp_idx[i]] >> 8) & 0xFF;
		frq = (res->vflag[res->cp_idx[i]] & 0xF);
		rtk->ssat[sat - 1].norm_v[0][frq] = res->norm_cp[i];
		rtk->ssat[sat - 1].norm_v[1][frq] = res->norm_pr[i];
	}
}

extern void freeres(res_t* res)
{
	res->npr = res->ncp = res->nv = 0;
	if (res->vflag) {
		free(res->vflag); res->vflag = nullptr;
	}
	if (res->pri_v) {
		free(res->pri_v); res->pri_v = nullptr;
	}
	if (res->post_v) {
		free(res->post_v); res->post_v = nullptr;
	}
	if (res->pr_idx) {
		free(res->pr_idx); res->pr_idx = nullptr;
	}
	if (res->cp_idx) {
		free(res->cp_idx); res->cp_idx = nullptr;
	}
	if (res->pri_pr) {
		free(res->pri_pr); res->pri_pr = nullptr;
	}
	if (res->pri_cp) {
		free(res->pri_cp); res->pri_cp = nullptr;
	}
	if (res->post_pr) {
		free(res->post_pr); res->post_pr = nullptr;
	}
	if (res->post_cp) {
		free(res->post_cp); res->post_cp = nullptr;
	}
	if (res->norm_pr) {
		free(res->norm_pr); res->norm_pr = nullptr;
	}
	if (res->norm_cp) {
		free(res->norm_cp); res->norm_cp = nullptr;
	}
	if (res->R) {
		free(res->R); res->R = nullptr;
	}
	if (res->Qvv) {
		free(res->Qvv); res->Qvv = nullptr;
	}
}

extern int resqc_igg_pr(rtk_t* rtk, res_t* res, int* exc, int ppp) {
	/*only check norm post pseudorange residual*/
	int frq, sat, qc_flag = 0;
	double el, k0 = 2.80, k1 = 4.13, fact;

	if (rtk->opt.igg_pr_k0 != 0.0) {
		k0 = rtk->opt.igg_pr_k0;
	}
	if (rtk->opt.igg_pr_k1 != 0.0) {
		k1 = rtk->opt.igg_pr_k1;
	}

	double max_n_pr;
	int max_n_pr_idx = findmax(res->norm_pr, res->npr, &max_n_pr);
	sat = (res->vflag[res->pr_idx[max_n_pr_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->pr_idx[max_n_pr_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	if (max_n_pr > k1) {
		rtk->ssat[sat - 1].var_fact[1][frq] = 100000.0;
		//exc[sat-1]=1;
		qc_flag = 1;
		trace(2, "%s(%d): %s P%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
			time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
			res->post_v[res->pr_idx[max_n_pr_idx]], max_n_pr, SQRT(res->R[res->pr_idx[max_n_pr_idx] + res->pr_idx[max_n_pr_idx] * res->nv]));
	}
	else if (max_n_pr >= k0 && max_n_pr <= k1) {
		fact = (max_n_pr / k0) * SQR((k1 - k0) / (k1 - max_n_pr));
		rtk->ssat[sat - 1].var_fact[1][frq] *= fact;
		qc_flag = 1;
		trace(3, "%s(%d): %s P%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
			time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
			res->post_v[res->pr_idx[max_n_pr_idx]], max_n_pr, SQRT(res->R[res->pr_idx[max_n_pr_idx] + res->pr_idx[max_n_pr_idx] * res->nv]), fact);
	}
	else {
		rtk->ssat[sat - 1].var_fact[1][frq] = 1.0;
	}

	return qc_flag;
}

extern int resqc_igg_cp(rtk_t* rtk, res_t* res, int* exc, int ppp) {
	/*only check norm post phase residual*/
	int frq, sat, qc_flag = 0;
	double el, k0 = 2.80, k1 = 4.13, fact;

	if (rtk->opt.igg_cp_k0 != 0.0) {
		k0 = rtk->opt.igg_cp_k0;
	}
	if (rtk->opt.igg_cp_k1 != 0.0) {
		k1 = rtk->opt.igg_cp_k1;
	}

#if 1
	double max_n_cp;
	int max_n_cp_idx = findmax(res->norm_cp, res->ncp, &max_n_cp);
	sat = (res->vflag[res->cp_idx[max_n_cp_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->cp_idx[max_n_cp_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	if (fabs(max_n_cp) > k1) {
		//rtk->ssat[sat - 1].var_fact[0][frq] = 100000.0;
		//exc[sat - 1] = 1;
		int iamb = (ppp ? iamb_ppp(&rtk->opt, sat, frq) : iamb_ppk(&rtk->opt, sat, frq));
		double amb = rtk->ssat[sat - 1].amb_from_pr[frq];
		initx(rtk, amb, rtk->P, SQR(rtk->opt.std[0]), iamb);

		qc_flag = 1;
		trace(2, "%s(%d): %s L%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
			time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
			res->post_v[res->cp_idx[max_n_cp_idx]], max_n_cp, SQRT(res->R[res->cp_idx[max_n_cp_idx] + res->cp_idx[max_n_cp_idx] * res->nv]));
	}
	else if (fabs(max_n_cp) >= k0 && fabs(max_n_cp) <= k1) {
		fact = (max_n_cp / k0) * SQR((k1 - k0) / (k1 - max_n_cp));
		rtk->ssat[sat - 1].var_fact[0][frq] *= fact;
		qc_flag = 1;
		trace(3, "%s(%d): %s L%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
			time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
			res->post_v[res->cp_idx[max_n_cp_idx]], max_n_cp, SQRT(res->R[res->cp_idx[max_n_cp_idx] + res->cp_idx[max_n_cp_idx] * res->nv]), fact);
	}
	else {
		rtk->ssat[sat - 1].var_fact[0][frq] = 1.0;
	}
#endif
	return qc_flag;
}

extern int resqc_igg(rtk_t* rtk, res_t* res, int* exc, int ppp) {
	int frq, sat, qc_flag = 0;
	double el, k0 = 2.80, k1 = 4.13, fact;

	if (rtk->opt.igg_k0 != 0.0) {
		k0 = rtk->opt.igg_k0;
	}
	if (rtk->opt.igg_k1 != 0.0) {
		k1 = rtk->opt.igg_k1;
	}

	double max_n_pr;
	int max_n_pr_idx = findmax(res->norm_pr, res->npr, &max_n_pr);
	sat = (res->vflag[res->pr_idx[max_n_pr_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->pr_idx[max_n_pr_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	if (max_n_pr > k1) {
		//rtk->ssat[sat - 1].var_fact[1][frq] = 100000.0;
		//rtk->ssat[sat - 1].init_amb[frq] = 1;
		exc[sat - 1] = 1;
		qc_flag = 1;
		trace(2, "%s(%d): %s P%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
			time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
			res->post_v[res->pr_idx[max_n_pr_idx]], max_n_pr, SQRT(res->R[res->pr_idx[max_n_pr_idx] + res->pr_idx[max_n_pr_idx] * res->nv]));
	}
	else {
		double max_n_cp;
		int max_n_cp_idx = findmax(res->norm_cp, res->ncp, &max_n_cp);
		sat = (res->vflag[res->cp_idx[max_n_cp_idx]] >> 8) & 0xFF;
		frq = (res->vflag[res->cp_idx[max_n_cp_idx]] & 0xF);
		el = rtk->ssat[sat - 1].azel[1] * R2D;
		if (fabs(max_n_cp) > k1) {
			rtk->ssat[sat - 1].var_fact[0][frq] = 100000.0;

			//int iamb = (ppp ? iamb_ppp(&rtk->opt, sat, frq) : iamb_ppk(&rtk->opt, sat, frq));
			//double amb = rtk->ssat[sat - 1].amb_from_pr[frq];
			//initx(rtk, amb, rtk->P, SQR(rtk->opt.std[0]), iamb);
			//rtk->ssat[sat - 1].init_amb[frq] = 1;
			//exc[sat - 1] = 1;
			qc_flag = 1;
			trace(2, "%s(%d): %s L%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
				time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
				res->post_v[res->cp_idx[max_n_cp_idx]], max_n_cp, SQRT(res->R[res->cp_idx[max_n_cp_idx] + res->cp_idx[max_n_cp_idx] * res->nv]));
		}
		else if (fabs(max_n_cp) >= k0 && fabs(max_n_cp) <= k1) {
			fact = (max_n_cp / k0) * SQR((k1 - k0) / (k1 - max_n_cp));
			rtk->ssat[sat - 1].var_fact[0][frq] *= fact;
			//exc[sat - 1] = 1;
			qc_flag = 1;
			trace(3, "%s(%d): %s L%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
				time_str(rtk->sol.time, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), frq + 1, el,
				res->post_v[res->cp_idx[max_n_cp_idx]], max_n_cp, SQRT(res->R[res->cp_idx[max_n_cp_idx] + res->cp_idx[max_n_cp_idx] * res->nv]), fact);
		}
		else {
			rtk->ssat[sat - 1].var_fact[0][frq] = 1.0;
		}
	}
	return qc_flag;
}

extern int resqc_shi(rtk_t* rtk, res_t* res, int* exc, int ppp) {
	int frq, sat, qc_flag = 0, iamb;
	double el;

	/*STEP 1: check post pseudorange residual*/
	double max_pr;
	int max_pr_idx = findmax(res->post_pr, res->npr, &max_pr);
	sat = (res->vflag[res->pr_idx[max_pr_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->pr_idx[max_pr_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	if (fabs(max_pr) > 3.0 / sin(el * D2R)) {
		exc[sat - 1] = 1;
		qc_flag = 1;
		return qc_flag;
	}

	/*STEP 2: check norm pseudorange residual*/
	double max_n_pr;
	int max_n_pr_idx = findmax(res->norm_pr, res->npr, &max_n_pr);
	sat = (res->vflag[res->pr_idx[max_n_pr_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->pr_idx[max_n_pr_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	if (fabs(max_n_pr) > 2.0) {
		exc[sat - 1] = 1;
		qc_flag = 1;
		return qc_flag;
	}

	/*STEP 3: check post phase residual*/
	double max_cp;
	int max_cp_idx = findmax(res->post_cp, res->ncp, &max_cp);
	sat = (res->vflag[res->cp_idx[max_cp_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->cp_idx[max_cp_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	iamb = (ppp ? iamb_ppp(&rtk->opt, sat, frq) : iamb_ppk(&rtk->opt, sat, frq));
	if (fabs(max_cp) > 0.03 / sin(el * D2R)) {
		initx(rtk, rtk->x[iamb], rtk->P, SQR(rtk->opt.std[0]), iamb);
		rtk->ssat[sat - 1].init_amb[frq] = 1;
		qc_flag = 1;
		return qc_flag;
	}

	/*STEP 4: check norm post phase residual*/
	double max_n_cp;
	int max_n_cp_idx = findmax(res->norm_cp, res->ncp, &max_n_cp);
	sat = (res->vflag[res->cp_idx[max_n_cp_idx]] >> 8) & 0xFF;
	frq = (res->vflag[res->cp_idx[max_n_cp_idx]] & 0xF);
	el = rtk->ssat[sat - 1].azel[1] * R2D;
	iamb = (ppp ? iamb_ppp(&rtk->opt, sat, frq) : iamb_ppk(&rtk->opt, sat, frq));
	if (fabs(max_n_cp) > 2.0) {
		initx(rtk, rtk->x[iamb], rtk->P, SQR(rtk->opt.std[0]), iamb);
		rtk->ssat[sat - 1].init_amb[frq] = 1;
		qc_flag = 1;
		return qc_flag;
	}
	return 0;
}

extern int resqc_zhao(gtime_t t, rtk_t* rtk, const double* v, const double* norm_v, const int* v_flag, int nv, int* exc) {
	int i, j, k, m = 0, n = 0, type, frq, nf = NF(&rtk->opt), satj;
	double vv_p[MAXOBS * 2] = { 0 }, vv_c[MAXOBS * 2] = { 0 }, vv_np[MAXOBS * 2] = { 0 }, vv_nc[MAXOBS * 2] = { 0 };
	int vv_sat_p[MAXOBS] = { 0 }, vv_sat_c[MAXOBS] = { 0 }, vv_frq_p[MAXOBS] = { 0 }, vv_frq_c[MAXOBS] = { 0 };
	double vv_el_p[MAXOBS] = { 0 }, vv_el_c[MAXOBS] = { 0 };
	double max_vp, max_vc, max_np, max_nc, max;
	int max_idx_vp, max_idx_vc, max_idx_np, max_idx_nc;
	double K0 = 1.5, K1 = 2.5;
	int qc_flag = 0;

	for (i = 0, j = 0, k = 0; i < nv; i++) {
		satj = (v_flag[i] >> 8) & 0xFF;
		type = (v_flag[i] >> 4) & 0xF;
		frq = (v_flag[i] & 0xF);

		if (type == 1) {
			vv_p[j] = fabs(v[i]);
			vv_np[j] = fabs(norm_v[i]);
			vv_frq_p[j] = frq;
			vv_el_p[j] = rtk->ssat[satj - 1].azel[1];
			vv_sat_p[j++] = satj;
		}
		else if (type == 0) {
			vv_c[k] = fabs(v[i]);
			vv_nc[k] = fabs(norm_v[i]);
			vv_frq_c[k] = frq;
			vv_el_c[k] = rtk->ssat[satj - 1].azel[1];
			vv_sat_c[k++] = satj;
		}
	}
	int sat = 0;
	double el = 0.0;
	double fact = 1.0;
	/* pseudorange t-test*/
	double avp_vnp = 0.0, delta_avp_vnp = 0.0;
	for (i = 0; i < j; i++) {
		avp_vnp += vv_np[i];
	}
	avp_vnp /= j;
	for (i = 0; i < j; i++) {
		delta_avp_vnp += SQR(vv_np[i] - avp_vnp);
	}

	double thres1 = tdistb_0050[j], thres2 = tdistb_0010[j];  ///0050 and 0010 for cpt
	double T_p[MAXOBS * 2] = { 0 };
	for (i = 0; i < j; i++) {
		sat = vv_sat_p[i];
		el = vv_el_p[i];
		T_p[i] = fabs(vv_np[i] - avp_vnp) / SQRT(delta_avp_vnp / j);
		if (T_p[i] > thres1 && T_p[i] < thres2) {
			fact = thres1 / T_p[i] * SQR((thres2 - T_p[i]) / (thres2 - thres1));
			rtk->ssat[sat - 1].var_fact[0][frq] = 1.0 / fact;
			qc_flag = 1;
			trace(2, "%s(%d): %s down weight by larger pseudorange residual el=%4.1f Tp=%5.2f thres1=%5.2f thres2=%5.2f\n",
				time_str(t, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), el * R2D, T_p[i], thres1, thres2);
		}
		else if (T_p[i] > thres2) {
			//            rtk->ssat[sat-1].var_fact=10000.0;
			exc[sat - 1] = 1;
			qc_flag = 1;
			trace(2, "%s(%d): %s exclude by larger pseudorange residual el=%4.1f Tp=%5.2f thres=%5.2f\n",
				time_str(t, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), el * R2D, T_p[i], thres2);
		}
		else {
			rtk->ssat[sat - 1].var_fact[0][frq] = 1.0;
			qc_flag = 0;
		}
	}

	/* phase t-test */
	double avp_vnc = 0.0, delta_avp_vnc = 0.0;
	for (i = 0; i < k; i++) {
		avp_vnc += vv_nc[i];
	}
	avp_vnc /= k;
	for (i = 0; i < k; i++) {
		delta_avp_vnc += SQR(vv_nc[i] - avp_vnc);
	}
	thres1 = tdistb_0050[k], thres2 = tdistb_0010[k];
	double T_c[MAXOBS * 2] = { 0 };
	for (i = 0; i < k; i++) {
		sat = vv_sat_c[i];
		el = vv_el_c[i];
		T_c[i] = fabs(vv_nc[i] - avp_vnc) / SQRT(delta_avp_vnc / k);
		if (T_c[i] > thres1 && T_c[i] < thres2) {
			fact = thres1 / T_c[i] * SQR((thres2 - T_c[i]) / (thres2 - thres1));
			rtk->ssat[sat - 1].var_fact[0][frq] = 1.0 / fact;
			qc_flag = 1;
			trace(2, "%s(%d): %s down weight by larger phase residual el=%4.1f Tc=%5.2f thres1=%5.2f thres2=%5.2f\n",
				time_str(t, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), el * R2D, T_c[i], thres1, thres2);
		}
		else if (T_c[i] > thres2) {
			exc[sat - 1] = 1;
			qc_flag = 1;
			trace(2, "%s(%d): %s exclude by larger phase residual el=%4.1f Tp=%5.2f thres=%5.2f\n",
				time_str(t, 1), rtk->tc ? rtk->ins_kf->couple_epoch : rtk->epoch, sat_id(sat), el * R2D, T_c[i], thres2);
		}
	}
	return qc_flag;
}
extern int resqc(gtime_t t, rtk_t* rtk, res_t* res, int* exc, int ppp) {
	if (rtk->opt.robust == ROBUST_QC_OFF) {
		return 0;
	}
	else if (rtk->opt.robust == ROBUST_QC_IGG_PR) {
		return resqc_igg_pr(rtk, res, exc, ppp);
	}
	else if (rtk->opt.robust == ROBUST_QC_IGG_CP) {
		return resqc_igg_cp(rtk, res, exc, ppp);
	}
	else if (rtk->opt.robust == ROBUST_QC_IGG) {
		return resqc_igg(rtk, res, exc, ppp);
	}
	else if (rtk->opt.robust == ROBUST_QC_SHI) {
		return resqc_shi(rtk, res, exc, ppp);
	}
	else return 0;
}

extern int pri_res_check(gtime_t t, rtk_t* rtk, const double* pri_v, const int* vflag, int nv, int* exc)
{
	double v_SYS[NUM_SYS][MAXOBS] = { 0 }, v_copy[MAXOBS] = { 0 }, mean, thres = 20.0;
	int nv_SYS[NUM_SYS] = { 0 }, sat_SYS[NUM_SYS][MAXOBS] = { 0 };
	int i, j, sat, sys_idx, qc_flag = 0;
	int ppk = 0, ppp = 0;
	prcopt_t* popt = &rtk->opt;

	ppk = ((popt->mode >= PMODE_DGPS && popt->mode <= PMODE_STATIC_START) ||
		(popt->mode == PMODE_TC_DGPS || popt->mode == PMODE_TC_PPK) ||
		(popt->mode == PMODE_LC_DGPS || popt->mode == PMODE_LC_PPK) ||
		(popt->insopt.imu_align == INS_ALIGN_GNSS_PPK || popt->insopt.imu_align == INS_ALIGN_GNSS_DGPS));
	ppp = ((popt->mode >= PMODE_PPP_KINEMA && popt->mode <= PMODE_PPP_FIXED) || (popt->mode == PMODE_TC_PPP || popt->mode == PMODE_LC_PPP));

	if (ppk) {
		thres = 3.0;
	}
	if (ppp) {
		thres = 100.0;
	}

	for (j = 0; j < nv; j++) {
		sat = (vflag[j] >> 8) & 0xFF;
		if (j == 64)
		{
			j = j;
		}
		sys_idx = satsysidx(sat);
		if (sys_idx == -1)continue;
		sat_SYS[sys_idx][nv_SYS[sys_idx]] = sat;
		v_SYS[sys_idx][nv_SYS[sys_idx]++] = pri_v[j];
	}

	for (j = 0; j < NUM_SYS; j++) {
		if (nv_SYS[j] <= 0) continue;
		matcpy(v_copy, v_SYS[j], nv_SYS[j], 1);
		mean = median(v_copy, nv_SYS[j]);
		for (i = 0; i < nv_SYS[j]; i++) {
			if (fabs(v_SYS[j][i] - mean) > thres) {
				sat = sat_SYS[j][i];
				exc[sat - 1] = 1;
				qc_flag = 1;
			}
		}
	}

	return qc_flag;
}
//JXY qc

static bool sysErrPrc_v1(rtk_t* rtk, const nav_t* nav, const int sat, const int fqu, const statIndicator_t* sI, int* exc)
{
	bool bPrc = false;
	int  nmsg = 0, j = 0,sys_idx = satsysidx(sat),fre;
	char chmsg[20][MAXCHARS] = { '\0' };
	double delev, dave, dstd, dsr, dre;

	int index; //, ib = rtk->para.idAmb_sat2ix[sat - 1][fq] IB_PPP(sat,fq,&rtk->opt);
	bool bInitRe = false, bVarInf = false;
	char msg[MAXCHARS];
	msg[0] = '\0';
	int prn = 0; satsys(sat, &prn);
	//////////////////////////////////////////////////////////////////////////
	for (int fq = 0; fq < NFREQ; fq++)
	{
		index = sI[sat - 1].iStResi_PR[fq];
		if (index>0&&index<4)
		{
			if (rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_cur[fq] <= 4.0)
			{
				rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq] *= 4.0;
			}
			else if (rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_cur[fq] >= 36.0)
			{
				rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq] *= 36.0;
			}
			else
			{
				rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq] *= rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_cur[fq];
			}
			if (rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq] > 5000.0)
			{
				rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq] = 5000.0;
			}

			dave = rtk->ppp_glo.stResiAve[OTYPE_PR_LP];
			dstd = rtk->ppp_glo.stResiStd[OTYPE_PR_LP];
			delev = rtk->ssat[sat - 1].azel[1] * R2D;
			dsr = rtk->ppp_glo.equS[sat - 1].eS_pr.stResi[fq];
			dre = rtk->ppp_glo.equS[sat - 1].eS_pr.resi_1[fq];

			if (satsys(sat, &prn) == SYS_GLO)exc[sat - 1] = 1;
			trace(3, "prn%d SR_P%d %6.2f %6.2f %6.2f %5.2f el=%4.1f %4.1f \0", prn, fq + 1,
				dsr, dre, dave, dstd, delev, rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_fin[fq]);
			bPrc = true;
		}
		if (index >= 4) {
			exc[sat - 1] = 1;
			return true;
		}

		//////////////////////////////////////////////////////////////////////////

		bInitRe = false;
		index = sI[sat - 1].iResi[fq];

		if (1 < index)
		{
			if (3 == index)
			{
				int nf = NF(&rtk->opt);
				for (int f = 0; f < nf; f++)
				{

					fre = rtk->opt.gnss_frq_idx[sys_idx][f] - 1;
					rtk->ssat[sat - 1].slip[fre] = true;

					/* reinitialize phase-bias if detecting cycle slip */
					if (rtk->opt.ionoopt == IONOOPT_UC)
					{
						j = rtk->tc ? xiIon(&rtk->opt.insopt, sat) : II(sat, &rtk->opt);
						initx(rtk, rtk->ssat[sat - 1].ion, rtk->P, VAR_IONO, j);
					}
					int iamb = rtk->tc ? xiAmb(&rtk->opt.insopt, sat, f) : IB(sat, f, &rtk->opt);
					initx(rtk, rtk->ssat[sat - 1].amb_from_pr[f], rtk->P, VAR_BIAS, iamb);

					rtk->ssat[sat - 1].lock[fre] = -rtk->opt.minlock;
					if (satsys(sat, &prn) == SYS_GLO)exc[sat - 1] = 1;
					//exc[sat - 1] = 1;
					bInitRe = true;
					trace(3, "Init sat %c%d amb=%2f(L%d)", satno2sys(sat), prn, f, rtk->ssat[sat - 1].amb_from_pr[f]);
				}

			}
			if (index >= 4) {
				exc[sat - 1] = 1;
				return true;
			}


			if (bVarInf || bInitRe)
			{
				dave = rtk->ppp_glo.resiAve_C;
				dstd = rtk->ppp_glo.resiStd_C;
				delev = rtk->ssat[sat - 1].azel[1] * R2D;
				dsr = rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[fq]; // 相位标准化残差
				dre = rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[fq]; // 验后相位残差

				trace(3, "prn%d Re_L%d %6.2f %6.2f %6.2f %5.2f el=%4.1f\0", prn, fq + 1,
					dsr, dre, dave, dstd, delev);
				bPrc = true;
			}
		}

		///*
		if (bInitRe)
		{
			int nf = NF(&rtk->opt);
			for (int f = 0; f < nf; f++)
			{
				//rtk->ssat[sat-1].slip[f]=true;
				rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_cur[f] = 1.0;
				rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[f] = 1.0;
			}
		}//*/
		else
		{
			index = sI[sat - 1].iStResi[fq];
			if (index>0&&index<4)
			{
				if (rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_cur[fq] <= 5.0)
				{
					rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq] *= 5.0;
				}
				else if (rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_cur[fq] >= 25.0)
				{
					rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq] *= 25.0;
				}
				else
				{
					rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq] *= rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_cur[fq];
				}
				if (rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq] > 225.0)
				{
					rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq] = 225.0;
				}

				dave = rtk->ppp_glo.stResiAve[OTYPE_CP_LP];
				dstd = rtk->ppp_glo.stResiStd[OTYPE_CP_LP];
				delev = rtk->ssat[sat - 1].azel[1] * R2D;
				dsr = rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[fq];
				dre = rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[fq];

				trace(3, "prn%d SR_L%d %6.2f %6.2f %6.2f %5.2f el=%4.1f %4.1f \0", prn, fq + 1,
					dsr, dre, dave, dstd, delev, rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_fin[fq]);
				//exc[sat - 1] = 1;
				bPrc = true;
			}
			if(index>=4)
			{
				exc[sat - 1] = 1;
				return true;
			}
		}
	}
	return bPrc;
}
static int cmpres(const void* p1, const void* p2)
{
	double* q1 = (double*)p1, * q2 = (double*)p2;
	double delta = *q1 - *q2;
	return delta < -0.0 ? -1 : (delta > 0.0 ? 1 : 0);
}

static int analyseModelIndex_v1(rtk_t* rtk, const nav_t* navs, const statIndicator_t* sI, int* exc)
{
	bool bPrc = false;
	int i, j = 0, k = 0, sat, worstSat, worstSat2 = -1, nBad = 0, satBad[MAXSAT] = { 0 }, fqBad[MAXSAT] = { 0 };
	double dBadMax = -1.0, dBadAmounts[MAXSAT] = { 0 }, dsat_bad_counts[MAXSAT] = { 0 }, res_phase[MAXSAT] = { 0 };
	char msg[MAXCHARS] = { '\0' };
	//////////////////////////////////////////////////////////////////////////

	int prn = 0,nsat=0;

	for (i = 0; i < MAXSAT; i++)
	{
		dBadAmounts[i] = 0.0;
		fqBad[i] = -1;

		if (!rtk->ppp_glo.bUsed_sat[i]) continue;
		satsys(i, &prn);
		double sum1 = 0.0,res_ph=0.0;
		double sum2 = 0.0;
		double dmax = -999.9;

		for (int f = 0; f < NFREQ; f++)
		{
			double dtmp = 0.0;
			dtmp += sI[i].iResi[f];
			dtmp += sI[i].iStResi[f];
			dtmp += sI[i].iStResi_PR[f] * 0.5;

			if (dtmp > dmax)
			{
				dmax = dtmp;
				fqBad[i] = f;
			}
			if (res_ph < rtk->ppp_glo.equS[i].eS_cp.resi_1[f])
			{
				res_ph = rtk->ppp_glo.equS[i].eS_cp.resi_1[f];
			}
			sum1 += sI[i].iResi[f];
			sum1 += sI[i].iStResi[f];
			sum1 += sI[i].iStResi_PR[f] * 0.5;

			sum2 += sI[i].iResi[f];
			sum2 += sI[i].iStResi[f];
		}

		dBadAmounts[i] += sum1;
		nsat++;
		// jxy 20200406
		//if (dBadAmounts[i] > 0.0)
		//	dBadAmounts[i]  += rtk->ssat[i].outc[0] - 1;

		int outc = rtk->ssat[i].outc[0] - 1;
		if (dBadAmounts[i] > 0.0)dBadAmounts[i] += outc;

		if (dBadAmounts[i] > dBadMax)
		{
			dBadMax = dBadAmounts[i];
			worstSat = i + 1;
		}

		if (0.0 < dBadAmounts[i])
		{
			satBad[nBad] = i + 1;
			dsat_bad_counts[nBad] = sum1;
			res_phase[nBad] = res_ph;
			nBad += 1;
		}
	}
	if (0 == nBad|| dBadMax==1)
		return 0;

	if (1 == nBad)
	{
		bPrc = sysErrPrc_v1(rtk, navs, worstSat, fqBad[worstSat - 1], sI, exc);
		satsys(worstSat, &prn);
		trace(3, "\n>>>>>>>>>>>>>>>>>>>>>> Sat %s%d observation data quality is poor\n", satno2sys(worstSat), prn);

		if (!bPrc)
			return 0;
	}
	else
	{
		//Panda 20141018
		//(周跳修复+钟差建模)处理grcb0020.12o 00:31:00时刻，G25伪距存在粗差，导致12 15 22号卫星相位残差异常，模型分析到这里时
		//并没有真正找到问题所在
		//
		for (j = 0; j < nBad - 1; j++) for (k = j + 1; k < nBad; k++) {
			if (dsat_bad_counts[j] >= dsat_bad_counts[k]) continue;
			SWAP_I(satBad[j], satBad[k]);
			SWAP_I(dsat_bad_counts[j], dsat_bad_counts[k]);
			SWAP_D(res_phase[j], res_phase[k]);
		}

		int  nbadsat_same = 0, bad_delsatno[MAXSAT] = { 0 };
		for (i = 0; i < nBad; i++)
		{
			if (dsat_bad_counts[i] == dsat_bad_counts[0]|| dsat_bad_counts[i]>=6)
			{
				nbadsat_same++;
				bad_delsatno[i] = satBad[i];
			}
		}
		int prn = 0;
		if (nbadsat_same >= 1)
		{		
			if (nbadsat_same >= nsat/3)
			{
				trace(3, ">>>>>>>>>>>>>>>>>>>>>> Urgent!!! Urgent!!! %d of the %d satellites were of poor quality! \n", nbadsat_same,nsat);
				for (j = 0; j < nbadsat_same; j++)
				{
					if (res_phase[j] < 0.03)continue;
					bPrc = sysErrPrc_v1(rtk, navs, bad_delsatno[j], fqBad[bad_delsatno[j] - 1], sI, exc);
					//bPrc = sysErrPrc_v1(rtk, navs, bad_delsatno[j], 1, sI, exc);
					satsys(bad_delsatno[j], &prn);
					trace(3, "\n>>>>>>>>>>>>>>>>>>>>>> Sat %s%d observation data quality is poor\n", satno2sys(bad_delsatno[j]), prn);
				}
			}
			else {
				trace(3, "\n>>>>>>>>>>>>>>>>>>>>>> Poor observation data more than 1,count is %d\n", nbadsat_same);
				for (j = 0; j < nbadsat_same; j++)
				{
					bPrc = sysErrPrc_v1(rtk, navs, bad_delsatno[j], fqBad[bad_delsatno[j] - 1], sI, exc);
					//bPrc = sysErrPrc_v1(rtk, navs, bad_delsatno[j], 1, sI, exc);
					satsys(bad_delsatno[j], &prn);
					trace(3, "\n>>>>>>>>>>>>>>>>>>>>>> Sat %s%d observation data quality is poor\n", satno2sys(bad_delsatno[j]), prn);
				}
			}
		}
		else {
			worstSat = bad_delsatno[0];
			//bPrc = sysErrPrc_v1(rtk, navs, worstSat, fqBad[worstSat - 1], sI, exc);
			bPrc = sysErrPrc_v1(rtk, navs, worstSat, fqBad[worstSat - 1], sI, exc);
			//bPrc = sysErrPrc_v1(rtk, navs, worstSat, 1, sI, exc);
			satsys(worstSat, &prn);
			trace(3, "\n>>>>>>>>>>>>>>>>>>>>>> Sat %s%d observation data quality is poor\n", satno2sys(worstSat), prn);

		}
	}

	if (!bPrc)
		return 0;

	return nBad;
}

static bool qcPost_standResi_1sat(const double dfabsmax, const double dfac, const double std, int* ibad)
{
	bool bBad = false;
	double faclimited = 1.0;

	*ibad = 0;

	if (dfabsmax <= 1.25)
		return bBad;
	else if (dfabsmax <= 1.5)
	{
		if (std > 0.35)		faclimited = 4.5;
		else if (std > 0.22)  faclimited = 6.0;
		else                faclimited = 8.0;

		if (dfac > faclimited)
		{
			bBad = true;
			*ibad = 4;
		}
	}
	else if (dfabsmax <= 2.5)
	{
		if (std > 0.5)        faclimited = 2.5;
		else if (std > 0.35)  faclimited = 3.0;
		else if (std > 0.22)  faclimited = 4.0;
		else                faclimited = 5.0;

		if (dfac > faclimited)
		{
			bBad = true;
			*ibad = 3;
		}
	}
	else if (dfabsmax <= 5.0)
	{
		bBad = true;
		*ibad = 2;
	}
	else
	{
		bBad = true;
		*ibad = 1;
	}

	return bBad;
}

static bool qcPost_standResi_PR_1sat(const double dfabsmax, const double dfac, const double std, int* ibad)
{
	bool bBad = false;
	double facLimited = 1.0;

	*ibad = 0;

	if (dfabsmax < 1.25)
		return bBad;
	else if (dfabsmax <= 1.5)
	{
		if (std > 0.35)		facLimited = 3.5;
		else if (std > 0.22)  facLimited = 5.0;
		else                facLimited = 7.0;

		if (dfac > facLimited)
		{
			*ibad = 4;
			bBad = true;
		}
	}
	else if (dfabsmax <= 2.0)
	{
		if (std > 0.35)		facLimited = 2.0;
		else if (std > 0.22)  facLimited = 3.0;
		else                facLimited = 4.0;

		if (dfac > facLimited)
		{
			*ibad = 3;
			bBad = true;
		}
	}
	else if (dfabsmax <= 3.0)
	{
		facLimited = 1.75;

		if (dfac > facLimited)
		{
			*ibad = 2;
			bBad = true;
		}
	}
	else
	{
		*ibad = 1;
		bBad = true;
	}

	return bBad;
}

static bool qcPost_ModelError_1sat(const double dfabsmax, const double dfac, int* ibad)
{
	*ibad = 0;
	bool bBad = false;
	double factor_limit = 1.0;

	if (dfabsmax <= 1.0)
	{
		return bBad;
	}
	else if (dfabsmax <= 1.5)
	{
		factor_limit = 7.5;
		if (dfac > factor_limit)
		{
			*ibad = 4;
			bBad = true;
		}
	}
	else if (dfabsmax <= 2.0)
	{
		factor_limit = 4.5;
		if (dfac > factor_limit)
		{
			*ibad = 3;
			bBad = true;
		}
	}
	else if (dfabsmax <= 3.0)
	{
		factor_limit = 3.0;

		if (dfac > factor_limit)
		{
			*ibad = 2;
			bBad = true;
		}
	}
	else
	{
		*ibad = 1;
		bBad = true;
	}

	return bBad;
}

static bool calEQW(const bool& bPR, const double& v, const double& sig, double* standResi, double* eqw)
{//                     伪距识别符     伪距验后残差    伪距验后残差方差    标准化伪距残差  伪距残差降权因子
	bool bBiased = false;
	double dd = 0.0, k0 = 1.25;;
	const double k1 = 5.0;
	const double maxfactor = 64.0;

	if (fabs(sig) < 1.0e-10)
		return bBiased;

	*standResi = dd = v / sig;//标准化残差
	dd = fabs(dd);//取绝对值

	if (bPR)k0 = 1.0;//伪距残差阈值
	//标准化残差判断，经验值
	if (dd <= k0)
	{
		(*eqw) = 1.0;//残差无异常
	}
	else
	{
		//残差异常
		bBiased = true;

		if (dd <= k1)
		{
			double dtmp = 0.0;
			dtmp = dd * (k1 - k0) * (k1 - k0) / (k0 * (k1 - dd) * (k1 - dd));//放大因子
			if (dtmp > maxfactor)(*eqw) = maxfactor;
			else(*eqw) = dtmp;
		}
		else
		{
			(*eqw) = maxfactor;
		}
	}

	return bBiased;
}
extern double calStd_ex(const double* v, const int& n, bool* bused, double* ave, double* rms)
{
	int i = 0, m = 0;
	double dave = 0.0, std = 0.0;

	for (i = m = 0; i < n; i++)
	{
		if (bused) if (false == bused[i])
			continue;

		dave += v[i];
		m++;
	}

	if (m <= 0)
		return -1.0;

	dave /= m;

	for (i = 0; i < n; i++)
	{
		if (bused) if (false == bused[i])
			continue;

		std += (v[i] - dave) * (v[i] - dave);
	}

	std = sqrt(std / m);

	if (ave)
		*ave = dave;

	if (rms)
	{
		double drms = 0.0;
		for (i = 0; i < n; i++)
		{
			if (bused) if (false == bused[i])
				continue;

			drms += v[i] * v[i];
		}
		drms = SQRT(drms / m);
		*rms = drms;
	}

	return std;
}
//计算标准化残差                +=标准化  -残差    +相位-AVE      +相位-STD        +伪距-AVE      +伪距-STD
static void calStResi_AveStd_v1(rtk_t* rtk, double* aveExC, double* stdExC, double* aveExP, double* stdExP)
{
	int i, f, nC, nP, iWorstP, iWorstC;
	double dvC[MAXSAT] = { 0 }, dvP[MAXSAT] = { 0 };
	bool bUsedC[MAXSAT] = { 0 }, bUsedP[MAXSAT] = { 0 };
	double aveTmpC, stdTmpC, aveExTmpC, stdExTmpC, dMaxC;
	double aveTmpP, stdTmpP, aveExTmpP, stdExTmpP, dMaxP;

	nC = nP = 0;
	aveTmpC = stdTmpC = aveExTmpC = stdExTmpC = dMaxC = 0.0;
	aveTmpP = stdTmpP = aveExTmpP = stdExTmpP = dMaxP = 0.0;
	dMaxP = dMaxC = -1.0;

	for (i = 0; i < MAXSAT; i++)
	{
		for (f = 0; f < NFREQ; f++)
		{
			rtk->ppp_glo.equS[i].eS_pr.stResi[f] = 0.0;
			rtk->ppp_glo.equS[i].eS_cp.stResi[f] = 0.0;
		}
	}

	for (i = 0; i < rtk->ppp_glo.nOEI; i++)
	{
		int sat = rtk->ppp_glo.oEI[i].sat;//卫星索引
		int ot = rtk->ppp_glo.oEI[i].ty; //残差类型 0 相位； 1 伪距
		int f = rtk->ppp_glo.oEI[i].fq;  //频率索引
		if (f == 2)continue;//待添加
		//for ( f=0;f<NFREQ;f++ ) {
		if (OTYPE_PR_LP == ot)//伪距
		{
			if (rtk->ppp_glo.equS[sat - 1].eS_pr.bUsed[f])//sat的f频率
			{
				// Panda 20140930
				// 计算伪距标准化残差不检查vsat。在L超限的情况下，vsat为false，但P可能临近超限
											   //sat的f频率验后伪距残差                   sat的f频率验后伪距残差方差
				calEQW(true, rtk->ppp_glo.equS[sat - 1].eS_pr.resi_1[f], rtk->ppp_glo.equS[sat - 1].eS_pr.sig_vp[f],
					//            sat的f频率标准化验后伪距残差               sat的f频率的降权因子
					&rtk->ppp_glo.equS[sat - 1].eS_pr.stResi[f], &rtk->ppp_glo.equS[sat - 1].eS_pr.eWF_cur[f]);

				dvP[nP] = rtk->ppp_glo.equS[sat - 1].eS_pr.stResi[f];
				aveTmpP += dvP[nP];       //标准化验后伪距残差之和

				if (fabs(dvP[nP]) > dMaxP)
				{
					dMaxP = fabs(dvP[nP]);//记录最大标准化验后伪距残差最大值
					iWorstP = nP;         //记录最大标准化验后伪距残差最大值索引
				}

				nP++;//标准化验后伪距残差数量
			}
		}
		else if (OTYPE_CP_LP == ot)//相位
		{
			if (rtk->ppp_glo.equS[sat - 1].eS_cp.bUsed[f])
			{
				// 计算标准化残差，倒数第二个参数
												//sat的f频率验后相位残差                   sat的f频率验后相位残差方差
				calEQW(false, rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[f], rtk->ppp_glo.equS[sat - 1].eS_cp.sig_vp[f],
					//            sat的f频率标准化验后相位残差               sat的f频率的降权因子
					&rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[f], &rtk->ppp_glo.equS[sat - 1].eS_cp.eWF_cur[f]);
				// 计算标准化残差
				dvC[nC] = rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[f];
				aveTmpC += dvC[nC];//标准化验后相位残差之和

				if (fabs(dvC[nC]) > dMaxC)
				{
					dMaxC = fabs(dvC[nC]);//记录最大标准化验后相位残差最大值
					iWorstC = nC;		  //记录最大标准化验后相位残差最大值索引
				}

				nC++;                     //标准化验后相位残差数量
			}
		}
	}

	//伪距
	if (1 < nP)
	{
		aveTmpP /= nP;                 //标准化验后伪距残差AVE
		for (i = 0; i < nP; i++)
		{
			bUsedP[i] = true;
			stdTmpP += (dvP[i] - aveTmpP) * (dvP[i] - aveTmpP);
		}
		stdTmpP = sqrt(stdTmpP / nP); //标准化验后伪距残差STD

		if (aveExP || stdExP)
		{
			if (nP >= 4)
			{
				bUsedP[iWorstP] = false;
				stdExTmpP = calStd_ex(dvP, nP, bUsedP, &aveExTmpP, 0);//stdExTmpP-排除最大残差后的标准化验后伪距残差STD；aveExTmpP-.....均值
			}
			else
			{
				aveExTmpP = aveTmpP;
				stdExTmpP = stdTmpP;
			}
		}
	}
	//相位  算法同伪距
	if (1 < nC)
	{
		aveTmpC /= nC;
		for (i = 0; i < nC; i++)
		{
			bUsedC[i] = true;
			stdTmpC += (dvC[i] - aveTmpC) * (dvC[i] - aveTmpC);
		}

		stdTmpC = sqrt(stdTmpC / nC);

		if (aveExC || stdExC)
		{
			if (nC >= 4)
			{
				bUsedC[iWorstC] = false;
				stdExTmpC = calStd_ex(dvC, nC, bUsedC, &aveExTmpC, 0);
			}
			else
			{
				aveExTmpC = aveTmpC;
				stdExTmpC = stdTmpC;
			}
		}
	}

	if (aveExC) *aveExC = aveExTmpC;
	if (stdExC) *stdExC = stdExTmpC;
	if (aveExP) *aveExP = aveExTmpP;
	if (stdExP) *stdExP = stdExTmpP;
}

// 计算验后标准化残差                                 验后相位残差均值 STD          验后伪距残差AVE   STD
static void calPostResi_AveStd_v1(rtk_t* rtk, double* aveExC, double* stdExC, double* aveExP, double* stdExP/*, statIndicator_t *sI*/)
{          //伪距数量 载波  最大残差索引
	int i, nC, nP, iCWorst = 0, iPWorst = 0; // jxy 20180310
	double dvC[MAXSAT] = { 0 }, dvP[MAXSAT] = { 0 };//伪距/相位残差和
	bool bUsedC[MAXSAT] = { 0 }, bUsedP[MAXSAT] = { 0 };//C/P可用索引？
	double aveTmpC = 0, stdTmpC = 0, aveExTmpC = 0, stdExTmpC = 0, dMaxC = 0;//C 平均值 STD Ex(exclude?) 最大值
	double aveTmpP = 0, stdTmpP = 0, aveExTmpP = 0, stdExTmpP = 0, dMaxP = 0;//P

	//初始化
	nC = nP = 0;
	aveTmpC = stdTmpC = aveExTmpC = stdExTmpC = dMaxC = 0.0;
	aveTmpP = stdTmpP = aveExTmpP = stdExTmpP = dMaxP = 0.0;
	//观测值数量
	for (i = 0; i < rtk->ppp_glo.nOEI; i++)
	{
		int sat = rtk->ppp_glo.oEI[i].sat;//索引
		int ot = rtk->ppp_glo.oEI[i].ty;//观测值类型 0-载波 1-伪距
		int fq = rtk->ppp_glo.oEI[i].fq;//频率索引
		if (fq == 2)continue;//暂不处理第三频率（待改进）
		//sI[sat-1].iResi=0;
		if (OTYPE_PR_LP == ot)// 伪距
		{
			dvP[nP] = rtk->ppp_glo.equS[sat - 1].eS_pr.resi_1[fq];//验后伪距残差
			aveTmpP += dvP[nP];//验后伪距残差之和  P1+P2

			if (fabs(dvP[nP]) > dMaxP)//记录最大验后伪距残差
			{
				dMaxP = fabs(dvP[nP]);  //最大验后伪距残差
				iPWorst = nP;           //最大验后伪距残差索引
			}

			nP++;                       //伪距观测值数量
		}
		else if (OTYPE_CP_LP == ot)	// 载波
		{
			dvC[nC] = rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[fq];//验后相位残差
			aveTmpC += dvC[nC];//验后相位残差之和 L1+L2

			if (fabs(dvC[nC]) > dMaxC)//记录最大相位残差
			{
				dMaxC = fabs(dvC[nC]);//最大相位残差
				iCWorst = nC;        //最大相位残差索引
			}
			nC++;                    //相位观测值数量
		}
	}
	//伪距
	if (1 < nP)
	{
		aveTmpP /= nP;//验后伪距残差均值
		for (i = 0; i < nP; i++)
		{
			bUsedP[i] = true;
			stdTmpP += (dvP[i] - aveTmpP) * (dvP[i] - aveTmpP);
		}
		stdTmpP = sqrt(stdTmpP / nP);//验后伪距残差STD

		if (stdExP || aveExP)
		{
			if (nP >= 4)
			{
				bUsedP[iPWorst] = false;
				stdExTmpP = calStd_ex(dvP, nP, bUsedP, &aveExTmpP, 0);//stdExTmpP-排除最大残差后的验后伪距残差STD，aveExTmpP排除最大残差后的验后伪距残差AVE
			}
			else
			{
				aveExTmpP = aveTmpP;
				stdExTmpP = stdTmpP;
			}
		}
	}
	//相位
	if (1 < nC)
	{
		aveTmpC /= nC;//验后相位残差均值
		for (i = 0; i < nC; i++)
		{
			bUsedC[i] = true;
			stdTmpC += (dvC[i] - aveTmpC) * (dvC[i] - aveTmpC);
		}
		stdTmpC = sqrt(stdTmpC / nC);//验后相位残差STD

		if (stdExC || aveExC)
		{
			if (nC >= 4)
			{
				bUsedC[iCWorst] = false;
				stdExTmpC = calStd_ex(dvC, nC, bUsedC, &aveExTmpC, 0);//stdExTmpC-排除最大残差后的验后相位残差STD，aveExTmpC排除最大残差后的验后相位残差AVE
			}
			else
			{
				aveExTmpC = aveTmpC;
				stdExTmpC = stdTmpC;
			}
		}
	}

	if (aveExC) *aveExC = aveExTmpC;
	if (stdExC) *stdExC = stdExTmpC;
	if (aveExP) *aveExP = aveExTmpP;
	if (stdExP) *stdExP = stdExTmpP;
}

static bool qcPost_standResi_PR_1sat_v1(const double dfabsmax, const double dfac, const double std, int* ibad)
{
	bool bBad = false;
	double facLimited = 1.0;

	//return bBad;

	*ibad = 0;

	if (dfabsmax < 1.0)
		return bBad;
	else if (dfabsmax <= 1.33)
	{
		if (std > 0.35)		facLimited = 5.0;
		else if (std > 0.22)  facLimited = 6.5;
		else                facLimited = 8.0;

		if (dfac > facLimited)
		{
			*ibad = 1;
			bBad = true;
		}
	}
	else if (dfabsmax <= 2.0)
	{
		if (std > 0.35)		facLimited = 4.0;
		else if (std > 0.22)  facLimited = 5.0;
		else                facLimited = 6.0;

		if (dfac > facLimited)
		{
			*ibad = 2;
			bBad = true;
		}
	}
	else if (dfabsmax <= 3.0)
	{
		facLimited = 2.5;

		if (dfac > facLimited)
		{
			*ibad = 3;
			bBad = true;
		}
	}
	else
	{
		*ibad = 4;
		bBad = true;
	}

	return bBad;
}

//Panda 20151027
//验后伪距标准化残差质量控制
static bool qcPost_standResi_PR_v1(rtk_t* rtk, const double* ave, const double* std, int* iWorstSat, statIndicator_t* sI)
{                                     //     标准化验后残差ave     验后标准化残差std
	int i;
	bool bBad = false, bRes;
	double dtmp, dfabs, dfac, dRecipStd, dmax = 0.0;

	*iWorstSat = -1;//初始化
	if (*std < 1.0e-20) return false;
	if (fabs(*ave) < 1.0e-20) return false;

	dmax = 0.0;
	dRecipStd = 1.0 / (*std);//临界值

	for (i = 0; i < rtk->ppp_glo.nOEI; i++)
	{
		int sat = rtk->ppp_glo.oEI[i].sat;//卫星索引
		int ot = rtk->ppp_glo.oEI[i].ty;//观测值类型
		int fq = rtk->ppp_glo.oEI[i].fq;//观测值频率
		if (fq == 2)continue;

		if (!rtk->ppp_glo.equS[sat - 1].eS_pr.bUsed[fq]) continue;

		dtmp = rtk->ppp_glo.equS[sat - 1].eS_pr.stResi[fq];//标准伪距化验后残差
		dfabs = fabs(dtmp);                                //取绝对值
		dfac = fabs(dtmp - (*ave)) * dRecipStd;              //高斯缩放法

		//if ( fabs(rtk->ppp_glo.equS[sat-1].eS_pr.resi_1[fq])>10.0 ) {
		//	getchar();
		//}

		bRes = qcPost_standResi_PR_1sat_v1(dfabs, dfac, *std, sI[sat - 1].iStResi_PR + fq);//sat的fq频率的观测值质量量化

		if (!bRes)continue;

		bBad = true;

		if (dfabs > dmax)
		{
			dmax = dfabs;
			*iWorstSat = sat;
		}
	}

	return bBad;
}

static bool qcPost_standResi_1sat_v1(const double dfabsmax, const double dfac, const double std, int* ibad)
{
	bool bBad = false;
	double faclimited = 1.0;

	*ibad = 0;

	if (dfabsmax <= 1.0)
		return bBad;
	else if (dfabsmax <= 1.5)
	{
		if (std > 0.35)		faclimited = 5.0;
		else if (std > 0.22)  faclimited = 7.0;
		else                faclimited = 9.0;

		if (dfac > faclimited)
		{
			bBad = true;
			*ibad = 1;
		}
	}
	else if (dfabsmax <= 2.5)
	{
		if (std > 0.5)        faclimited = 2.5;
		else if (std > 0.35)  faclimited = 3.0;
		else if (std > 0.22)  faclimited = 4.0;
		else                faclimited = 5.0;

		if (dfac > faclimited)
		{
			bBad = true;
			*ibad = 2;
		}
	}
	else if (dfabsmax <= 5.0)
	{
		bBad = true;
		*ibad = 3;
	}
	else
	{
		bBad = true;
		*ibad = 4;
	}

	return bBad;
}
//算法同伪距
static bool qcPost_standResi_v1(rtk_t* rtk, const double* ave, const double* std, int* iWorstSat, statIndicator_t* sI)
{
	bool badv = false;
	int i;
	double dtmp, dfabs, dfac, dRecipStd, dmax = 0.0; //,dMax,dMaxEx;

	*iWorstSat = -1;
	if (*std < 1.0e-20)       return false;
	if (fabs(*ave) < 1.0e-20) return false;

	dmax = 0.0;
	dRecipStd = 1.0 / (*std);

	for (i = 0; i < rtk->ppp_glo.nOEI; i++)
	{
		int sat = rtk->ppp_glo.oEI[i].sat;
		int ot = rtk->ppp_glo.oEI[i].ty;
		int fq = rtk->ppp_glo.oEI[i].fq;
		if (fq == 2)continue;

		if (OTYPE_CP_LP != ot) continue;

		sI[sat - 1].iStResi[fq] = 0;

		dtmp = fabs(rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[fq]);
		if (dtmp < 0.0001) continue;
		else if (dtmp < 0.003)
		{
			if (rtk->ppp_glo.equS[sat - 1].eS_cp.sig_vp[fq] < 0.001)
				continue;
		}

		dtmp = rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[fq];
		dfabs = fabs(dtmp);
		dfac = fabs(dtmp - (*ave)) * dRecipStd;//标准化de标准化  （x-ave）/方差

		if (qcPost_standResi_1sat_v1(dfabs, dfac, (*std), sI[sat - 1].iStResi + fq))
		{
			badv = true;

			if (dfabs > dmax)
			{
				dmax = dfabs;
				*iWorstSat = sat;
			}
		}
	}

	return badv;
}

static bool qcPost_Resi_1sat_v1(const double lam, const double _limitedMax_, const double _limitB_,
	const double _limitC_, const double dFac, const double dFabs, const double dFabsSResi, int* ibad)
{
	bool badv = false;
	double factor = 1.0;

	*ibad = 0;

	if (dFabs <= 0.005)
		return badv;
	else if (dFabs <= 0.03)
	{
		if (dFabsSResi > 1.5)		factor = 4.0;
		else if (dFabsSResi > 1.15)	factor = 5.0;
		else if (dFabsSResi > 0.8)  factor = 6.0;
		else                        factor = 9.0;

		if (dFac > factor)
		{
			*ibad = 1;
			badv = true;
		}
	}
	else if (dFabs <= 0.06)  		//如果最大的卫星残差大于阈值的一半，且大于2.5倍的ave或rms0，将其方差斜方差增大9倍
	{
		if (dFabsSResi > 1.5)		factor = 3.0;
		else if (dFabsSResi > 1.15)	factor = 4.0;
		else if (dFabsSResi > 0.8)	factor = 5.0;
		else						factor = 8.0;

		if (dFac > factor)
		{
			*ibad = 2;
			badv = true;
		}
	}
	else  			//如果最大的残差大于0.25m，直接初始化
	{
		bool bSlip = false;

		/*if (  rtk->opt_ex.bElevLmtEx ) {
		if ( rtk->ssat[sat-1].azel[1]<opt->elmin ) {
		bool bVarBad=rtk->ppp_glo.sigma0>chisqr_0p75[rtk->ppp_glo.nst-1];
		bool bSResiBad=dFabsSResi>0.6;

		if ( false==bVarBad || false==bSResiBad )
		bSlip=true;
		}
		}*/

		if (false == bSlip)
		{
			*ibad = 3;
			badv = true;
		}
	}

	return badv;
}

///*
static bool qcPost_Resi_v1(rtk_t* rtk, const nav_t* nav, const double* ave, const double* std, int* iWorstSat, statIndicator_t* sI)
{
	int i;
	bool bBad;
	double dTmp, dFac, dFabs, dFabsSResi, dRecipStd, dMax = 0.0;

	//////////////////////////////////////////////////////////////////////////
	const double limitedMax = 0.10;
	const double limitB = 0.04;
	const double limitC = 0.01;
	//////////////////////////////////////////////////////////////////////////

	bBad = false;
	*iWorstSat = -1;

	if (*std < 1.0e-20) return false;
	if (fabs(*ave) < 1.0e-20) return false;

	dMax = 0.0;
	dRecipStd = 1.0 / (*std);

	for (i = 0; i < rtk->ppp_glo.nOEI; i++)
	{
		int sat = rtk->ppp_glo.oEI[i].sat;
		int ot = rtk->ppp_glo.oEI[i].ty;
		int fq = rtk->ppp_glo.oEI[i].fq;
		if (fq == 2)continue;

		if (OTYPE_CP_LP != ot) continue;

		dTmp = rtk->ppp_glo.equS[sat - 1].eS_cp.resi_1[fq];
		dFabs = fabs(dTmp);
		dFac = fabs(dTmp - (ave[0])) * dRecipStd;
		dFabsSResi = fabs(rtk->ppp_glo.equS[sat - 1].eS_cp.stResi[fq]);

		//if ( qcPost_Resi_1sat_v1(rtk->ppp_glo.lam[i][fq], limitedMax, limitB, limitC, dFac, dFabs, dFabsSResi, sI[sat - 1].iResi + fq ) ) // jxy
		if (qcPost_Resi_1sat_v1(0/*nav->lam[i][fq]*/, limitedMax, limitB, limitC, dFac, dFabs, dFabsSResi, sI[sat - 1].iResi + fq))
		{
			if (dFabs > dMax)
			{
				dMax = dFabs;
				*iWorstSat = sat;
			}
			bBad = true;
		}
	}

	return bBad;
}

extern bool qualityControl_v1(rtk_t* rtk, const nav_t* navs, const int iter, const  int* v_flag, const int nv, int* exc)
{
	int i = 0, sat1 = 0;
	bool bbadqc = false;
	statIndicator_t sI[MAXSAT] = { 0 };
	char msg[MAXCHARS];
	//for (i = 0; i < MAXSAT; i++)rtk->ppp_glo.bUsed_sat[i]=false;//排除未使用卫星
	//for (i = 0; i < nv; i++)
	//{
	//	sat1 = (v_flag[i] >> 8) & 0xFF;
	//	rtk->ppp_glo.bUsed_sat[sat1] = true;
	//}
	rtk->ppp_glo.bQcNewCsFound = false;
	rtk->ppp_glo.iWorstSat_Resi_C = -1;
	rtk->ppp_glo.iWorstSat_stResi[OTYPE_CP_LP] = -1;
	rtk->ppp_glo.iWorstSat_stResi[OTYPE_PR_LP] = -1;

	for (i = 0; i < MAXSAT; i++)
	{
		if (!rtk->ppp_glo.bUsed_sat[i]) //排除未使用卫星
		{
			continue;
		}
		for (int f = 0; f < NFREQ; f++)
		{
			sI[i].iMErr[f] = 0;     //模型化误差指标	    0：正常
			sI[i].iResi[f] = 0;     //相位残差指标			0：正常		越大表示越严重
			sI[i].iStResi[f] = 0;   //相位标准化残差指标	0：正常		越大表示越严重
			sI[i].iStResi_PR[f] = 0;//伪距标准化残差指标	0：正常		越大表示越严重
		}
		sI[i].iter = iter;
	}
	//计算验后残差均值和STD               验后相位残差均值        验后相位残差std
	calPostResi_AveStd_v1(rtk, &rtk->ppp_glo.resiAve_C, &rtk->ppp_glo.resiStd_C, NULL, NULL);
	// 计算标准化验后残差均值和STD          标准化验后相位残差均值               标准化验后相位残差std                 标准化验后伪距残差ave                 标准化验后伪距残差std
	calStResi_AveStd_v1(rtk, &rtk->ppp_glo.stResiAve[OTYPE_CP_LP], &rtk->ppp_glo.stResiStd[OTYPE_CP_LP], &rtk->ppp_glo.stResiAve[OTYPE_PR_LP], &rtk->ppp_glo.stResiStd[OTYPE_PR_LP]);
	if (rtk->epoch == 1239)
	{
		i = i;
	}
	if (1/*rtk->opt.robust!=0*/)
	{
		bool bBad = false;

		if (2 >= iter)
		{//                                     标准化验后残差ave           伪距          标准化验后残差std     伪距              伪距质量检核失败卫星索引             伪距质量指标
			bBad = qcPost_standResi_PR_v1(rtk, &rtk->ppp_glo.stResiAve[OTYPE_PR_LP], &rtk->ppp_glo.stResiStd[OTYPE_PR_LP], &rtk->ppp_glo.iWorstSat_stResi[OTYPE_PR_LP], sI);
			bbadqc = bbadqc || bBad;

			if (!bBad)
			{//                                      标准化验后残差ave      相位          标准化验后残差std     相位              相位质量检核失败卫星索引             伪距质量指标
				bBad = qcPost_standResi_v1(rtk, &rtk->ppp_glo.stResiAve[OTYPE_CP_LP], &rtk->ppp_glo.stResiStd[OTYPE_CP_LP], &rtk->ppp_glo.iWorstSat_stResi[OTYPE_CP_LP], sI);
				bbadqc = bbadqc || bBad;
				//                                      相位验后残差ave           相位验后残差std        相位质量检核失败卫星索引
				bBad = qcPost_Resi_v1(rtk, navs, &rtk->ppp_glo.resiAve_C, &rtk->ppp_glo.resiStd_C, &rtk->ppp_glo.iWorstSat_Resi_C, sI);
				bbadqc = bbadqc || bBad;
			}
		}
		else if (iter < 10)
		{
			bBad = qcPost_standResi_PR_v1(rtk, &rtk->ppp_glo.stResiAve[OTYPE_PR_LP], &rtk->ppp_glo.stResiStd[OTYPE_PR_LP], &rtk->ppp_glo.iWorstSat_stResi[OTYPE_PR_LP], sI);
			bbadqc = bbadqc || bBad;
			bBad = qcPost_standResi_v1(rtk, &rtk->ppp_glo.stResiAve[OTYPE_CP_LP], &rtk->ppp_glo.stResiStd[OTYPE_CP_LP], &rtk->ppp_glo.iWorstSat_stResi[OTYPE_CP_LP], sI);
			bbadqc = bbadqc || bBad;
			bBad = qcPost_Resi_v1(rtk, navs, &rtk->ppp_glo.resiAve_C, &rtk->ppp_glo.resiStd_C, &rtk->ppp_glo.iWorstSat_Resi_C, sI);
			bbadqc = bbadqc || bBad;
		}

		if (bbadqc)
		{
			if (analyseModelIndex_v1(rtk, navs, sI, exc))
			{
				trace(3, "\nobservation data quality check  iter=%d\n", iter + 1);
			}
			else {
				bbadqc = false;
			}
		}
	}

#ifdef KFPPP_BETA
	if (!bbadqc)
	{
		if (0 != rtk->ppp_glo.iModelStat)
		{
			sprintf(msgbuf[rtk->is], "%11s iModelStat=%d. xxxxxxxxx\n", "", rtk->ppp_glo.iModelStat);
			outDebug(rtk, true, rtk->pf.pfdebug);
		}
	}
#endif

	return bbadqc;
}
extern double NormDistribution(const double u) {
	if (u < -5.0) return 0.0;
	if (u > 5.0) return 1.0;

	double y = fabs(u) / sqrt(2.0);

	double p = 1.0 + y * (0.0705230784 + y * (0.0422820123 + y * (0.0092705272 +
		y * (0.0001520143 + y * (0.0002765672 + y * 0.0000430638)))));
	double er = 1.0 - pow(p, -16.0);
	p = (u < 0.0) ? 0.5 - 0.5 * er : 0.5 + 0.5 * er;
	return p;
}

extern double ReNorm(double p) {
	if (p == 0.5) return 0.0;
	if (p > 0.9999997) return 5.0;
	if (p < 0.0000003) return -5.0;
	if (p < 0.5) return -ReNorm(1.0 - p);

	double y = -log(4.0 * p * (1.0 - p));
	y = y * (1.570796288 + y * (0.3706987906e-1
		+ y * (-0.8364353589e-3 + y * (-0.2250947176e-3
			+ y * (0.6841218299e-5 + y * (0.5824238515e-5
				+ y * (-0.1045274970e-5 + y * (0.8360937017e-7
					+ y * (-0.3231081277e-8 + y * (0.3657763036e-10
						+ y * 0.6936233982e-12))))))))));
	return sqrt(y);
}

extern void MatMul(const char* tr, int n, int k, int m, double alpha,
	const double* A, const double* B, double beta, double* C)
{
	double d;
	int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m]; break;
		case 2: for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k]; break;
		case 3: for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m]; break;
		case 4: for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k]; break;
		}
		if (beta == 0.0) C[i + j * n] = alpha * d; else C[i + j * n] = alpha * d + beta * C[i + j * n];
	}
}
extern int max_element(const int n, const double* v, double* max_v) {
	int i = 0, inx = 0;
	for (i = 0; i < n; i++)if (*max_v < fabs(v[i]))*max_v = fabs(v[i]), inx = i;
	return inx;
}
extern int valins(const prcopt_t* opt, const double* x)
{
	const insopt_t* insopt = &opt->insopt;

	int na, ia, nba, iba, nbg, ibg;

	na = xnA(insopt); ia = xiA(insopt);
	nba = xnBa(insopt); iba = xiBa(insopt);
	nbg = xnBg(insopt); ibg = xiBg(insopt);

	if (norm(x + ia, na) > 5.0 * D2R || (nba ? norm(x + iba, 3) > 1E4 * MG2MPS2:0) ||
		(nbg ? norm(x + ibg, 3) > 5.0 * D2R:0)) {
		trace(2, "ins solutions validate failed\n");
		return 0;
	}
	return 1;
}