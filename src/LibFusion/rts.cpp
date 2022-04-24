/*------------------------------------------------------------------------------
* rts.cc : forward and backward smoothing function
* reference : [1] 严恭敏 捷联惯导算法与组合导航原理讲义, 2016.9
* version :
* history : Created by lizhen on 2021/5/20.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <vector>

using namespace std;

#define MAXSTATE 26


/* expand P mat for initialization */
extern void expand_pmat(double* P, int m, int n, double scale) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            P[i + j * m] *= scale;
        }
    }
}

/* sx=(fp+bp)-1*(fp*fx+bp*bx), sp=(fp^-1+bp^-1)^-1) */
extern void smooth(double* fx, double* fp, double* bx, double* bp, double* sx, double* sp, int mode) {
    int i;
    if (mode == 1) for (i = 0; i < 3; i++) bx[i] = -bx[i];
    else if (mode == 2) {   // slove euler angle
        if (fx[2] * bx[2] < 0.0) { // sign opposite
            if (bx[2] > 0.95 * PI) { bx[2] = bx[2] - 2 * PI; }
            else if (bx[2] < -0.95 * PI) { bx[2] = 2 * PI + bx[2]; }
        }
    }
    for (i = 0; i < 3; i++) {
        sp[i] = 1.0 / fp[i] + 1.0 / bp[i];
        sp[i] = 1.0 / sp[i];
        sx[i] = fp[i] * fx[i] + bp[i] * bx[i];
        sx[i] = sx[i] / (fp[i] + bp[i]);
    }

}

/* save forward result to temporary bin files */
extern void rts_save_forward(FILE* fp, const solins_t* sol, const insopt_t* opt) {
    double time;
    v3_t Enb, Qpos, Qvel, Qatt;
    vector<v3_t> x, P;

    time = (time2gpst(sol->time, NULL));

    m3_diag_v3(sol->Qpos, &Qpos);
    m3_diag_v3(sol->Qvel, &Qvel);

    solins_t solins_copy = *sol;
    //    attsync(&solins_copy,opt);   //todo cheak
    Enb.x = solins_copy.rpy.x;
    Enb.y = solins_copy.rpy.y;
    Enb.z = solins_copy.rpy.z;
    Qatt.x = (float)solins_copy.Qatt.m11;
    Qatt.y = (float)solins_copy.Qatt.m22;
    Qatt.z = (float)solins_copy.Qatt.m33;

    x.push_back(sol->pos); P.push_back(Qpos);
    x.push_back(sol->vel); P.push_back(Qvel);
    x.push_back(Enb);      P.push_back(Qatt);
    if (opt->est_ba) {
        x.push_back(sol->ba);
        P.push_back(sol->ba_std);
    }
    if (opt->est_bg) {
        x.push_back(sol->bg);
        P.push_back(sol->bg_std);
    }

    if (opt->est_sa) {
        x.push_back(sol->sa);
        P.push_back(sol->sa_std);
    }
    if (opt->est_sg) {
        x.push_back(sol->sg);
        P.push_back(sol->sg_std);
    }

    // format [x1*3, P1*3, x2*3, P2*3...]
    fwrite(&time, sizeof(double), 1, fp);
    for (size_t i = 0; i < x.size(); i++) {
        fwrite(&x[i], sizeof(v3_t), 1, fp);
        fwrite(&P[i], sizeof(v3_t), 1, fp);
    }
}

/* match forward result and variance matrix in bin file*/
extern int rts_match_forward(FILE* fp, const insopt_t* opt, double time, sol_t* sol) {
    int offset = 3 * 2; // pva num
    double tt;
    v3_t buff[MAXSTATE];

    if (opt->est_ba) offset += 1 * 2;
    if (opt->est_bg) offset += 1 * 2;
    if (opt->est_sg) offset += 1 * 2;
    if (opt->est_sa) offset += 1 * 2;

    if (time < 0.0) return 0;
    while (fread(&tt, sizeof(double), 1, fp)) {
        fread(buff, sizeof(v3_t), offset, fp);
        if (tt - time > 0.0025) {  // earlier than first forward result epoch
            fseek(fp, -sizeof(v3_t) * offset - sizeof(double), SEEK_CUR);
            break;
        }
        if (abs(time - tt) < 1e-4) {   // found match
            sol->time.time = time;
            for (int i = 0; i < 3; i++) {
                sol->rr[i] = buff[xiP(opt) / 3 * 2].v[i];
                sol->qr[i] = buff[xiP(opt) / 3 * 2 + 1].v[i];
                sol->rr[i + 3] = buff[xiV(opt) / 3 * 2].v[i];
                sol->qv[i] = buff[xiV(opt) / 3 * 2 + 1].v[i];
                sol->rpy[i] = buff[xiA(opt) / 3 * 2].v[i];
                sol->qa[i] = buff[xiA(opt) / 3 * 2 + 1].v[i];
                if (opt->est_ba) {
                    sol->ba[i] = buff[xiBa(opt) / 3 * 2].v[i];
                    sol->qba[i] = buff[xiBa(opt) / 3 * 2 + 1].v[i];
                }
                if (opt->est_bg) {
                    sol->bg[i] = buff[xiBg(opt) / 3 * 2].v[i];
                    sol->qbg[i] = buff[xiBg(opt) / 3 * 2 + 1].v[i];
                }
                if (opt->est_sa) {
                    sol->sa[i] = buff[xiSa(opt) / 3 * 2].v[i];
                    sol->qsa[i] = buff[xiSa(opt) / 3 * 2 + 1].v[i];
                }
                if (opt->est_sg) {
                    sol->sg[i] = buff[xiSg(opt) / 3 * 2].v[i];
                    sol->qsg[i] = buff[xiSg(opt) / 3 * 2 + 1].v[i];
                }
            }
            return 1;
        }
    }
    return 0;
}

/* smooth each variable */
extern void rts_smooth(const insopt_t* opt, sol_t* fwd, sol_t* bwd) {
    /// PVA
    smooth(fwd->rr, fwd->qr, bwd->rr, bwd->qr, bwd->rr, bwd->qr, 0);
    smooth(&fwd->rr[3], fwd->qv, &bwd->rr[3], bwd->qv, &bwd->rr[3], bwd->qv, 0);
    smooth(fwd->rpy, fwd->qa, bwd->rpy, bwd->qa, bwd->rpy, bwd->qa, 2);

    if (opt->est_ba) smooth(fwd->ba, fwd->qba, bwd->ba, bwd->qba, bwd->ba, bwd->qba, 1);
    if (opt->est_bg) smooth(fwd->bg, fwd->qbg, bwd->bg, bwd->qbg, bwd->bg, bwd->qbg, 1);

    //todo check
    if (opt->est_sa) smooth(fwd->sa, fwd->qsa, bwd->sa, bwd->qsa, bwd->sa, bwd->qsa, 1);
    if (opt->est_sg) smooth(fwd->sg, fwd->qsg, bwd->sg, bwd->qsg, bwd->sg, bwd->qsg, 1);

}


extern int rts_ins_gnss(const prcopt_t* popt, const filopt_t* fopt, const solopt_t* solopt, rtk_t* ins_rtk, pva_t* pva,
    const imu_t* imu, int* id, const obs_t* obss, const nav_t* navs) {
    /// back_id[0] = kimu_idx, kpva_idx, krover_idx, kbase_idx;

    /* Local variable definition ----------------------------------------------------- */
    int match_gnss, match_forward, ep_count = 0;
    double time;
    vector<sol_t> result;
    kf_t* inskf = ins_rtk->ins_kf;  /// for lc
    sol_t fwd_sol = { 0 };
    insopt_t opt = popt->insopt;
    imud_t* imuz = (imud_t*)malloc(sizeof(imud_t) * opt.zvopt.ws), imu_data = { 0 };
    int tc = popt->mode >= PMODE_TC_SPP && popt->mode <= PMODE_TC_PPP;
    int lc_obs = popt->mode >= PMODE_LC_SPP && popt->mode <= PMODE_LC_PPP;
    int lc = popt->mode >= PMODE_LC_POS && popt->mode <= PMODE_LC_PPP;

    int imu_idx = id[0] - 1, pva_idx = id[1] - 1;
    int rover_idx = id[2], base_idx = id[3];   /// todo check for tc

    /* Code start --------------------------------------------------------------------- */
    /* write header to output file */
    if (!outhead(fopt->rtsfile, nullptr, 0, popt, solopt)) return 0;

    FILE* fp_sol = fopen(fopt->rtsfile, "a");
    FILE* fp_fwd = fopen(fopt->rts_ins_fw, "rb");

    expand_pmat(inskf->P, inskf->nx, inskf->nx, 100.0);
    inskf->insstate->ba = v3_scalar(-1.0, inskf->insstate->ba);
    inskf->insstate->bg = v3_scalar(-1.0, inskf->insstate->bg);

    while (inputimu(opt.zvopt.ws, &imu_idx, imu, &imu_data, imuz)) {
        if (popt->te.time != 0.0 && timediff(popt->te, imu->data[imu_idx].time) < 0) {
            imu_idx--;
            continue;
        }
        if (popt->te.time != 0.0 && timediff(imu->data[imu_idx].time, popt->ts) < 0) break;

        /* ins mechanical */
        ep_count++;
        inskf_udstate(inskf, imu, &imu_idx, imu->property, &popt->insopt);
        inskf->imu_obs = &imu->data[imu_idx];

        /* match gnss information */
        match_gnss = matchgnssinfo(popt, inskf, &pva_idx, pva, &rover_idx, &base_idx, obss, navs, ins_rtk, imu);

        /* ins-gnss coupled */
        if (match_gnss) {
            if (lc) lcfilter(popt, pva_idx, pva, imu, imu_idx, inskf);
            else if (1) {
                /*                if(!tcfilter(popt,&obss,&navs,&krover_idx,&kbase_idx,&ins_rtk)){
                                    match_gnss=0;
                                    continue;*/
            }
        }

        //   solins2sol(&popt->insopt,inskf->sol,&ins_rtk->sol);  //todo check
           /* write solution */
        if (match_gnss) {
            result.push_back(ins_rtk->sol);
            //   outsol(fp_sol,&ins_rtk->sol,nullptr,solopt,popt);
        }
        else {
            /*            if(ep_count>=solopt->outfrq){
                            ep_count=0;
                            outsol(fp_sol,&ins_rtk->sol,nullptr,solopt,popt);*/

        }

        /* next epoch */
        imu_idx -= 1;
        match_gnss = 0;
    }


    for (size_t i = result.size(); i > 0; i--) {
        time = time2gpst(result[i].time, NULL);
        match_forward = rts_match_forward(fp_fwd, &opt, time, &fwd_sol);
        if (match_forward) {
            rts_smooth(&opt, &fwd_sol, &result[i]);
            // outsol(fp_sol, &result[i], nullptr, solopt, popt);   // todo check
        }
    }
    fclose(fp_fwd); remove(fopt->rts_ins_fw);
    fclose(fp_sol);
    return 0;
}

