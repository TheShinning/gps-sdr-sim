//
// Created by chenc on 2021/7/27.
//

#include "rtklib.h"

static int ig_stc_pr(const rtk_t* rtk, const double* v, const int nv, const double* H, const int* vflag, const double* R)
{
    int i, j, ip, nx, npr = 0, idx_pr[MAXOBS * 2] = { 0 }, vflag_[MAXOBS * 2] = { 0 }, type, info;
    double* H_, * v_, * R_, * xp, * Pp;
    solins_t insstate_copy = { 0 };

    /*extract pseudorange residuals*/
    for (i = 0; i < nv; i++) {
        type = (vflag[i] >> 4) & 0xF;
        if (type == 0) {
            idx_pr[npr] = i;
            vflag_[npr] = vflag[i];
            npr++;
        }
    }

    nx = rtk->ins_kf->nx;
    xp = mat(nx, 1); Pp = mat(nx, nx);
    insstate_copy = *rtk->ins_kf->insstate;

    /*construct matrix*/
    matcpy(xp, rtk->ins_kf->x, nx, 1);
    matcpy(Pp, rtk->ins_kf->P, nx, nx);
    v_ = mat(npr, 1); H_ = mat(nx, npr); R_ = mat(npr, npr);

    for (i = 0; i < npr; i++) {
        v_[i] = v[idx_pr[i]];
        for (j = 0; j < nx; j++) H_[j + i * nx] = 0.0;
        for (j = 0; j < npr; j++) R_[j + i * npr] = R[j + npr + (i + npr) * nv];

        ip = xiP(&rtk->opt.insopt);
        for (j = ip; j < ip + 3; j++) H_[j + i * nx] = -H[j - ip + rtk->nx * idx_pr[i]];
    }

    int tra = 0;
    if (rtk->ins_kf->couple_epoch >= 367000) {
        tra = 1;
    }

    /*kalman filter*/
    if ((info = filter(xp, Pp, H_, v_, R_, nx, npr, rtk->opt.pos_res_qc, rtk->opt.kfopt, nullptr, tra))) {
        trace(2, "%s(%d): filter error (info=%d)\n", time_str(rtk->ins_kf->sol->time, 1), rtk->ins_kf->couple_epoch, info);
        return 0;
    }


    /*update ins state by pr-only, float solution*/
    inskf_feedback(rtk->sol.time, SOLQ_IGTC, &rtk->opt.insopt, xp, &insstate_copy);
    insstate_copy.g_status = SOLQ_FLOAT;

    *rtk->ins_kf->insstate = insstate_copy;
    *rtk->ins_kf->sol = insstate_copy;
    matcpy(rtk->ins_kf->x, xp, nx, 1);
    matcpy(rtk->ins_kf->P, Pp, nx, nx);

    free(xp); free(Pp); free(H_); free(R_); free(v_);
    return 1;
}

static int ig_stc_cp(const rtk_t* rtk, const double* v, const int nv, const double* H, const int* vflag, const double* R)
{
    int i, j, ip, nx, ncp = 0, idx_cp[MAXOBS * 2] = { 0 }, vflag_[MAXOBS * 2] = { 0 }, type, info;
    double* H_, * v_, * R_, * xp, * Pp;
    solins_t insstate_copy = { 0 };

    /*extract pseudorange residuals*/
    for (i = 0; i < nv; i++) {
        type = (vflag[i] >> 4) & 0xF;
        if (type == 0) {
            idx_cp[ncp] = i;
            vflag_[ncp] = vflag[i];
            ncp++;
        }
    }

    nx = rtk->ins_kf->nx;
    xp = mat(nx, 1); Pp = mat(nx, nx);
    insstate_copy = *rtk->ins_kf->insstate;

    /*construct matrix*/
    matcpy(xp, rtk->ins_kf->x, nx, 1);
    matcpy(Pp, rtk->ins_kf->P, nx, nx);
    v_ = mat(ncp, 1); H_ = mat(nx, ncp); R_ = mat(ncp, ncp);

    for (i = 0; i < ncp; i++) {
        v_[i] = v[idx_cp[i]];
        for (j = 0; j < nx; j++) H_[j + i * nx] = 0.0;
        for (j = 0; j < ncp; j++) R_[j + i * ncp] = R[j + i * (nv)];

        ip = xiP(&rtk->opt.insopt);
        for (j = ip; j < ip + 3; j++) H_[j + i * nx] = -H[j - ip + rtk->nx * idx_cp[i]];
    }

    int tra = 0;
    if (rtk->ins_kf->couple_epoch >= 17000) {
        tra = 0;
    }

    /*kalman filter*/
    if ((info = filter(xp, Pp, H_, v_, R_, nx, ncp, rtk->opt.pos_res_qc, rtk->opt.kfopt, nullptr, tra))) {
        trace(2, "%s(%d): filter error (info=%d)\n", time_str(rtk->ins_kf->sol->time, 1), rtk->ins_kf->couple_epoch, info);
        return 0;
    }

    /*update ins state by pr-only, float solution*/
    inskf_feedback(rtk->sol.time, SOLQ_IGTC, &rtk->opt.insopt, xp, &insstate_copy);
    insstate_copy.g_status = SOLQ_FIX;

    *rtk->ins_kf->insstate = insstate_copy;
    *rtk->ins_kf->sol = insstate_copy;
    matcpy(rtk->ins_kf->x, xp, nx, 1);

    /*dont update variance matrix by fix solution*/
#if 1
    matcpy(rtk->ins_kf->P, Pp, nx, nx);
#endif

    free(xp); free(Pp); free(H_); free(R_); free(v_);
    return 1;
}

extern int ig_stc(const rtk_t* rtk, const double* v, const int nv, const double* H, const int* vflag, const double* R, const int upd_type)
{
    int stat = 0;
    if (upd_type == 0) {
        stat = ig_stc_pr(rtk, v, nv, H, vflag, R);
    }
    else {
        stat = ig_stc_cp(rtk, v, nv, H, vflag, R);
    }

    return stat;
}