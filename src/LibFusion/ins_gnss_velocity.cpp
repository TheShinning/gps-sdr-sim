/*-----------------------------------------------------------------------------
* ins-gnss-velocity.cc :
*
* version :
* history : Created by lizhen on 2021/4/17.
*----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <algorithm>
using namespace std;

/* constants -----------------------------------------------------------------*/
#define MAXIONV     1E3        /* max innovations for doppler measurement */
#define STDOPP      SQR(0.20)        /* doppler measurement standard variance */
#define UNC_CLKR    (60.0)     /* default initial receiver clock drift uncertainty (s/s) */
#define Mg2M        9.80665E-6          /* micro-g to meters per second squared */

/* functions for initial error covariance matrix ----------------------------*/
extern void initP(int is, int ni, int nx, double unc, double unc0, double* P0) {
    int i, j;
    for (i = is; i < is + ni; i++)
        for (j = 0; j < nx; j++) {
            if (j == i) P0[j + i * nx] = SQR(unc == 0.0 ? unc0 : unc);
            else P0[j + i * nx] = P0[i + j * nx] = 0.0;
        }
}

/* jacobian of doppler measurement by ins velocity---------------------------*/
static void jaco_dop_dv(const double* e, const double* rs, double* dopdv) {
    dopdv[0] = e[0];//-OMGE/CLIGHT*rs[1];
    dopdv[1] = e[1];//+OMGE/CLIGHT*rs[0];
    dopdv[2] = e[2];
}

/* sensitive matrix for doppler measurement----------------------------------*/
static int doppHVR(const prcopt_t* opt, const rtk_t* rtk, solins_t* sol, const double* rs,
    const double* dts, const double* var, const int* svh, int n, const obsd_t* obs,
    const nav_t* nav, double* v, double* H, double* R, double* fact) {
    double rr[3], vr[3], vs[3], e[3], azel[2], pos[3], rate, lam, dtr, * r, freq, cosel;
    double dopdv[3], a[3], E[9];
    int i, j, nv = 0, nx = rtk->nx;
    int irr, IV, NV;
    const insopt_t* opt_ = &opt->insopt;
    v3_t wbib;
    trace(3, "doppHVR:\n");

    r = mat(n, 1);
    wbib = v3_scalar(1 / rtk->ins_kf->idt, rtk->ins_kf->omgb);
    ins2gnss(sol, &rtk->ins_kf->insstate->arm_gps, &wbib, rr, vr, opt->insopt.mech_coord);
    ecef2pos(rr, pos); xyz2enu(pos, E);

    dtr = rtk->dtrr;  /* receiver clock drift */

    irr = xiRr(opt_); /* estimated states index */
    IV = xiV(opt_); NV = xnV(opt_);

    /* sensitive matrix for doppler measurement */
    for (i = 0; i < n && i < MAXOBS; i++) {
        freq = sat2freq(obs[i].sat, obs[i].code[0], nav);

        if (obs[i].D[0] == 0.0 || freq == 0.0 || norm(rs + 3 + i * 6, 3) <= 0.0) {
            continue;
        }
        if (!satsys(obs[i].sat, NULL)) continue;

        lam = CLIGHT / freq;
        /* geometric distance/azimuth/elevation angle */
        if (geodist(rs + i * 6, rr, e) <= 0.0 || satazel(pos, e, azel) < opt->elmin) continue;

        /* excluded satellite */
        if (satexclude(obs[i].sat, var[i], svh[i], opt)) continue;

        /* satellite velocity relative to receiver in ecef */
        for (j = 0; j < 3; j++) vs[j] = rs[j + 3 + i * 6] - vr[j];

        /* range rate with earth rotation correction */
        rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * vr[0] -
            rs[3 + i * 6] * rr[1] - rs[i * 6] * vr[1]);

        /* doppler residual */
        v[nv] = -lam * obs[i].D[0] - (rate + dtr - CLIGHT * dts[1 + i * 2]);

        if (fabs(v[nv]) < MAXIONV) {

            /* jacobian matrix */
            if (H) {
                jaco_dop_dv(e, rs + 6 * i, dopdv);
                for (j = IV; j < IV + NV; j++) H[j + nv * nx] = dopdv[j - IV];
                H[irr + nv * nx] = 1.0; /* receiver clock drift */
            }
            /* measurement variance */
            r[nv++] = STDOPP / sin(azel[1]);
        }
    }

    if (R) {
        for (i = 0; i < nv; i++) R[i + i * nv] = r[i] * fact[i];
        trace(3, "R=\n");
        tracemat(5, R, nv, nv, 15, 6);
    }
    if (nv && v) {
        trace(3, "v=\n");
        tracemat(5, v, nv, 1, 15, 6);
    }
    if (nv && H) {
        trace(3, "H=\n");
        tracemat(5, H, nx, nv, 15, 6);
    }
    free(r);
    return nv;
}

/* doppler updates solutions valid-------------------------------------------*/
static int valdopp(const insopt_t* opt, const double* x, const double* R,
    const double* v, int nv, double thres) {
    double fact = SQR(thres);
    int iA, nba, nbg, iba, ibg, i;

    trace(3, "valdoop:\n");

    iA = xiA(opt);
    nba = xnBa(opt); iba = xiBa(opt);
    nbg = xnBg(opt); ibg = xiBg(opt);

    /* check estimated states */
    if (norm(x + iA, 3) > 15.0 * D2R || (nba ? norm(x + iba, 3) > 1E5 * Mg2M : false)
        || (nbg ? norm(x + ibg, 3) > 30.0 * D2R : false)) {
        trace(2, "too large estimated state error\n");
        return 0;
    }
    /* post-fit residual test */
    for (i = 0; i < nv; i++) {
        if (v[i] * v[i] < fact * R[i + i * nv]) continue;
        trace(2, "large residual (v=%6.3f sig=%.3f)\n", v[i], SQRT(R[i + i * nv]));
        return 0;
    }
    return 1;
}

static int qc_vel(double* v, int nv, double* R, double* fact, int n) {
    int qc_flag = 0, i;
#if 0  // IGG3
    double max_num, * norm_v;
    double k0 = 2.0, k1 = 3.5;
    norm_v = mat(nv, 1);

    for (i = 0; i < nv; i++) norm_v[i] = v[i] / sqrt(R[i + i * n]);
    int max_id = findmax(norm_v, nv, &max_num);

    if (max_num > k1) {
        fact[max_id] = 1e9;
        qc_flag = 1;
        trace(3, "Doppler velocity: large norm residual in rejected segment %d v=%7.3f norm_v = %7.3f var=%7.3f\n",
            max_id, v[max_id], max_num, R[max_id + max_id * n]);
        return qc_flag;
    }
    else if (max_num >= k0 && max_num <= k1) {
        fact[max_id] = (max_num / k0) * SQR((k1 - k0) / (k1 - max_num));
        qc_flag = 1;
        trace(3, "Doppler velocity: large norm residual in reduced segment %d v=%7.3f norm_v = %7.3f var=%7.3f\n",
            max_id, v[max_id], max_num, R[max_id + max_id * nv]);
        return qc_flag;
    }
#endif
#if 1  // median absolute deviation
    double med, mad;
    double* vs;
    int k0 = 5, k1 = 10;

    vs = mat(nv, 1); matcpy(vs, v, nv, 1);
    sort(vs, vs + nv);
    if (nv % 2 == 1) {   // median
        med = vs[nv / 2];
    }
    else {
        med = (vs[nv / 2] + vs[nv / 2 - 1]) / 2;
    }
    matcpy(vs, v, nv, 1);
    for (i = 0; i < nv; i++) vs[i] = fabs(vs[i] - med) / 0.6745;
    sort(vs, vs + nv);
    if (nv % 2 == 1) {   // median absolute deviation
        mad = vs[nv / 2];
    }
    else {
        mad = (vs[nv / 2] + vs[nv / 2 - 1]) / 2;
    }

    if (mad < 0.02) mad = 0.02;

    for (i = 0; i < nv; i++) {
        if (fabs(v[i]) > k1 * mad && fact[i] < 2) {
            fact[i] = 1e3;
            qc_flag = 1;
            trace(3, "Doppler velocity: residual larger than median absolute deviation v=%7.3f mad=%7.3f\n", v[i], mad);
            return qc_flag;
        }
        else if (fabs(v[i]) >= k0 * mad && fact[i] && fabs(v[i]) <= k1 * mad && fact[i] && fact[i] < 2) {
            fact[i] = 5;
            qc_flag = 1;
            trace(3, "Doppler velocity: residual larger than median absolute deviation v=%7.3f mad=%7.3f\n", v[i], mad);
            return qc_flag;
        }
    }
#endif

    return qc_flag;
}

/* doppler filter for ins states updates-------------------------------------*/
static int doppfilt(rtk_t* rtk, const prcopt_t* opt, const obsd_t* obs, int n,
    const nav_t* nav, const double* rs, const double* dts,
    const double* var, const int* svh) {
    int i, j, qc_flag = 0;
    int nx = rtk->nx, na = 15, nv, info = 0, irr = 0, nrr = 0;
    double* v, * H, * R, * x, * P, * fact;
    static insstate_t inss;
    solins_t sol_copy = { 0 };

    trace(3, "doppfilt:\n");

    irr = xiRr(&opt->insopt);
    nrr = xnRr(&opt->insopt);

    x = zeros(nx, 1);  P = mat(nx, nx);
    v = mat(n, 1); R = zeros(n, n); H = zeros(n, nx);
    fact = mat(MAXOBS, 1);
    for (i = 0; i < MAXOBS; i++) fact[i] = 1.0;

    for (int iter = 0; iter < 5; iter++) {

        matcpy(x, rtk->x, rtk->na, 1);
        matcpy(P, rtk->P, nx, nx);
        /* initialize every epoch for clock drift (white noise) */
        initP(irr, nrr, nx, 0, UNC_CLKR, P);

        /* sensitive matrix for doppler measurement */
        nv = doppHVR(opt, rtk, rtk->ins_kf->insstate, rs, dts, var, svh, n, obs, nav, v, H, R, fact);

        if (nv <= 3 || norm(v, nv) / nv > 0.2) break;
        if (var[0] < 0.001) {
            int qq = 1;
        }

        /* receiver clock drift and non-zero due to close-loop */
        x[irr] = rtk->dtrr == 0 ? 1e-32 : rtk->dtrr;

        /* ekf filter */
        if (filter(x, P, H, v, R, nx, nv, 0, 0, nullptr, 0)) {
            trace(2, "filter error\n");
        }
        else {
            sol_copy = *rtk->ins_kf->insstate;
            /* postfit doppler residuals */
            if (doppHVR(opt, rtk, &sol_copy, rs, dts, var, svh, n, obs, nav, v, nullptr, nullptr, fact)) {
                qc_flag = qc_vel(v, nv, R, fact, n);
                if (qc_flag) continue;
                /* valid solutions */
                info = valdopp(&opt->insopt, x, R, v, nv, 4.0);

                /* close loop */
                if (info) {
                    inskf_feedback(rtk->ins_kf->time, SOLQ_DOP_AID, &opt->insopt, x, &sol_copy);
                    *rtk->ins_kf->insstate = sol_copy;
                    //    matcpy(rtk->P, P, nx, nx);
                    rtk->dtrr = x[irr];
                }
                if (!qc_flag) break;
            }
        }

    }

    free(x); free(P);
    free(v); free(H); free(R); free(fact);
    return info;
}


extern int doppler(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, const prcopt_t* opt) {
    double* rs, * dts, * var;
    int i, no, stat, svh[MAXOBS];
    obsd_t obsd[MAXOBS];

    trace(3, "doppler:\n");
    trace(4, "obs=\n");
    traceobs(4, obs, n);
    trace(5, "nav=\n");
    tracenav(5, nav);

    /* collect rover station observation */
    for (i = 0, no = 0; i < n; i++) {
        if (obs[i].rcv == 1) obsd[no++] = obs[i];
    }
    if (no <= 3) {
        trace(2, "no observation data\n");
        return 0;
    }

    rs = mat(6, n);
    dts = mat(2, n);
    var = mat(1, n);

    /* satellite positons,velocities and clocks */
    satposs(&rtk->opt, obs[0].time, obsd, no, nav, opt->sateph, rs, dts,NULL, var, svh);

    /* doppler measurement aid ins states updates */
    stat = doppfilt(rtk, opt, obsd, no, nav, rs, dts, var, svh);

    if (stat) {
        trace(3, "doppler aid ins ok\n");
    }
    else {
        trace(2, "doppler aid ins fail\n");
        stat = 0;
    }
    free(rs);
    free(dts);
    free(var);
    return stat;
}

