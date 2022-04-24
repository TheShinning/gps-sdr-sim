/*------------------------------------------------------------------------------
* inputobs.cc :
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/4/19.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"


extern int inputimu(int ws, int* imu_idx, const imu_t* imus, imud_t* imud, imud_t* imuz) {
    int i, j, k;
    if (*imu_idx < 0 || *imu_idx + 1 > imus->n) return 0;
    for (k = 0, i = *imu_idx; k < ws && i >= 0 && i < imus->n && imuz; k++, i++) {
        imuz[k] = imus->data[i];
    }
    *imud = imus->data[*imu_idx];
    return 1;
}


extern int inputgnss(obsd_t* obs, int solq, const prcopt_t* popt, const obs_t* obss, int* rover_idx, int* base_idx, int revs) {
    gtime_t time = { 0 };
    int i, nu, nr, n = 0;
    int base = (popt->mode == PMODE_DGPS || popt->mode == PMODE_KINEMA || popt->mode == PMODE_LC_DGPS ||
        popt->mode == PMODE_LC_PPK ||
        popt->mode == PMODE_TC_DGPS || popt->mode == PMODE_TC_PPK || popt->mode == PMODE_STC_PPK);

    if (*rover_idx < 0 || *rover_idx > obss->n) {
        return 0;
    }
    if (!revs) { /* input forward data */
        if ((nu = nextobsf(obss, rover_idx, 1)) <= 0) return -1;
        if (popt->intpref) {
            for (; (nr = nextobsf(obss, base_idx, 2)) > 0; base_idx += nr)
                if (timediff(obss->data[*base_idx].time, obss->data[*rover_idx].time) > -DTTOL) break;
        }
        else {
            for (i = *base_idx; (nr = nextobsf(obss, &i, 2)) > 0; *base_idx = i, i += nr)
                if (timediff(obss->data[i].time, obss->data[*rover_idx].time) > DTTOL) break;
        }
        nr = nextobsf(obss, base_idx, 2);
        if (nr <= 0) {
            nr = nextobsf(obss, base_idx, 2);
        }
        for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss->data[*rover_idx + i];
        if (base) {
            for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss->data[*base_idx + i];
        }
        *rover_idx += nu;
    }
    else { /* input backward data */
        if ((nu = nextobsb(obss, rover_idx, 1)) <= 0) return -1;
        if (popt->intpref) {
            for (; (nr = nextobsb(obss, base_idx, 2)) > 0; base_idx -= nr)
                if (timediff(obss->data[*base_idx].time, obss->data[*rover_idx].time) < DTTOL) break;
        }
        else {
            for (i = *base_idx; (nr = nextobsb(obss, &i, 2)) > 0; *base_idx = i, i -= nr)
                if (timediff(obss->data[i].time, obss->data[*rover_idx].time) < -DTTOL) break;
        }
        nr = nextobsb(obss, base_idx, 2);
        for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss->data[*rover_idx - nu + 1 + i];
        if (base) for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss->data[*base_idx - nr + 1 + i];
        *rover_idx -= nu;
    }
    return n;
}

extern int read_gnss_file(prcopt_t* popt, const filopt_t* fopt, int mode, int need_base, gtime_t ts, gtime_t te, obs_t* obss, nav_t* navs) {
    /// mode 0:other 1:ppk 2:ppp
    int n = 0, ti = 0;
    int index[MAXFILE] = { 0 };
    char* infile[MAXFILE] = { nullptr };
    static pcvs_t pcvss = { 0 }; /* antenna correction data */
    static sta_t stas[2];    /* station parameter */
    static sbs_t sbss = { 0 };   /* sbss messages */
    static lex_t lexs = { 0 };   /* QZSS LEX messages */

    for (int i = 0; i < MAXFILE; i++) {
        if (!(infile[i] = (char*)malloc(1024))) {
            for (; i >= 0; i--) free(infile[i]);
            return -1;
        }
        infile[i][0] = '\0';
    }

    /// rover
    strcpy(infile[n], fopt->robsf);
    index[n++] = n;
    /// base
    if (need_base) {
        strcpy(infile[n], fopt->bobsf);
        index[n++] = n;
    }


    /// brdc
    for (int i = 0; i < 3; i++) {
        if (strcmp(fopt->navf[i], "")) {
            strcpy(infile[n], fopt->navf[i]);
            index[n++] = n;
        }
    }

    if (!readobsnav_a(fopt,ts, te, ti, infile, index, n, popt, obss, navs, stas)) {
        /* free obs and nav data */
        freeobsnav(obss, navs);
        return 0;
    }

    /// sp3 and clk
    if (mode == 2) {  // PPP
        for (int i = 0; i < 3; i++) {
            if (strcmp(fopt->sp3f[i], "")) {
                strcpy(infile[n++], fopt->sp3f[i]);
            }
        }
        for (int i = 0; i < 3; i++) {
            if (strcmp(fopt->clkf[i], "")) {
                strcpy(infile[n++], fopt->clkf[i]);
            }
        }

        readpreceph(infile, n, popt, fopt, navs, &sbss, &lexs);

        /* read satellite antenna parameters */
        if (*fopt->atx && !(readpcv(fopt->atx, &pcvss))) {
            trace(1, "sat antenna pcv read error: %s\n", fopt->atx);
            return 0;
        }

        char path[1024];
        if (*fopt->eop) {
            free(navs->erp.data);
            navs->erp.data = nullptr;
            navs->erp.n = navs->erp.nmax = 0;
            reppath(fopt->eop, path, ts, "", "");
            if (!readerp(path, &navs->erp)) {
                trace(2, "no erp data %s\n", path);
            }
        }

        /* read dcb parameters */
        if (*fopt->dcb) {
            reppath(fopt->dcb, path, ts, "", "");
            readdcb(path, navs, stas);
        }
        /* set antenna parameters */
        if (popt->mode != PMODE_SINGLE && *fopt->atx) {
            setpcv(obss->n > 0 ? obss->data[0].time : timeget(), popt, navs, &pcvss, NULL, stas);
        }
        /* read ocean tide loading parameters */
        if (popt->mode > PMODE_SINGLE && *fopt->blq) {
            readotl(popt, fopt->blq, stas);
        }
    }

    for (int i = 0; i < MAXFILE; i++) free(infile[i]);

    return 1;
}

extern int matchgnsssol(gtime_t imu_time, int* sol_idx, const pva_t* pva, double dt, int mode) {
    /// mode = 0, forward; 1, backward
    int find_sol = 0, i;

    if (*sol_idx < 0 || *sol_idx > pva->n) return 0;

    double sow_ins, sow_sol;

    sow_ins = time2gpst(imu_time, NULL);

    if (mode == 0) {  // forward
        for (i = *sol_idx > 5 ? *sol_idx - 5 : 0; i < pva->n; i++) {
            sow_sol = time2gpst(pva->time[i], NULL);
            if (fabs(sow_ins - sow_sol) < dt) {
                find_sol = 1;
                *sol_idx = i;
                break;
            }
            else if (sow_sol - sow_ins > 2.0 * dt) {
                find_sol = 0;
                break;
            }
        }
    }
    else { // backward
        for (i = *sol_idx; i > 0; i--) {
            sow_sol = time2gpst(pva->time[i], NULL);
            if (fabs(sow_ins - sow_sol) < dt) {
                find_sol = 1;
                *sol_idx = i;
                break;
            }
            else if (sow_ins - sow_sol > 2.0 * dt) {
                find_sol = 0;
                break;
            }
        }
    }
    return find_sol;
}

extern int matchgnsssol1(gtime_t imu_time1, gtime_t imu_time2, int* sol_idx, const pva_t* pva, double dt, int mode) {
    /// mode = 0, forward; 1, backward
    int find_sol = 0, i;

    if (*sol_idx < 0 || *sol_idx > pva->n) return 0;

    double sow_ins1, sow_ins2, sow_sol;

    sow_ins1 = time2gpst(imu_time1, NULL);
    sow_ins2 = time2gpst(imu_time2, NULL);

    if (mode == 0) {  // forward
        for (i = *sol_idx > 5 ? *sol_idx - 5 : 0; i < pva->n; i++) {
            sow_sol = time2gpst(pva->time[i], NULL);
            if ((sow_sol - sow_ins1) >= 0.0 && (sow_sol - sow_ins2) <= 0.0) {
                find_sol = 1;
                *sol_idx = i;
                break;
            }
            else if (sow_sol - sow_ins1 > 10.0 * dt) {
                find_sol = 0;
                break;
            }
        }
    }
    else { // backward
        for (i = *sol_idx; i > 0; i--) {
            sow_sol = time2gpst(pva->time[i], NULL);
            if (fabs(sow_ins1 - sow_sol) < dt) {
                find_sol = 1;
                *sol_idx = i;
                break;
            }
            else if (sow_ins1 - sow_sol > 2.0 * dt) {
                find_sol = 0;
                break;
            }
        }
    }
    return find_sol;
}

extern int matchgnssobs(gtime_t imu_time, int* obs_idx, const obs_t* obs, double dt) {
    double sow_ins, sow_obs;
    int find_obs = 0, i;

    if (*obs_idx < 0 || *obs_idx > obs->n) return 0;
    for (i = *obs_idx - 100 < 0 ? 0 : *obs_idx - 10; i < obs->n; i++) {
        if (obs->data[i].rcv != 1) continue;
        sow_obs = time2gpst(obs->data[i].time, NULL);
        sow_ins = time2gpst(imu_time, NULL);
        if (fabs(sow_obs - sow_ins) < dt) {
            *obs_idx = i;
            find_obs = 1;
            break;
        }
        else if (sow_obs - sow_ins > 2.0 * dt) {
            find_obs = 0;
            break;
        }
    }
    return find_obs;
}

extern int matchgnssinfo(const prcopt_t* popt, kf_t* inskf, int* kpva_idx, pva_t* pva, int* krover_idx, int* kbase_idx,
    const obs_t* obss, const nav_t* navs, rtk_t* ins_rtk, const imu_t* imu)
{
    int match_gnss = 0, match_pva = 0, lc, tc;
    static gtime_t gps_loss_start = { 0 }, gps_loss_end = { 0 };
    tc = popt->mode >= PMODE_TC_SPP && popt->mode <= PMODE_TC_PPP;
    lc = popt->mode >= PMODE_LC_POS && popt->mode <= PMODE_LC_PPP;

    if (lc && inskf->ins_epoch != 1) {
        if (popt->mode == PMODE_LC_POS) {
            //            match_pva=matchgnsssol(inskf->imu_obs[inskf->nsample-1].time, kpva_idx, pva, 0.5/imu->property->freq_imu, 0);
            match_pva = matchgnsssol1(timeadd(inskf->imu_obs->time, -inskf->idt), inskf->imu_obs[inskf->nsample - 1].time, kpva_idx, pva, 0.5 / imu->property->freq_imu, 0);
            match_gnss = match_pva;
            if (match_pva) {
                inskf->last_couple_time = pva->time[*kpva_idx];
                inskf->insstate->g_status = pva->status[*kpva_idx];
            }
        }
        else if (popt->mode >= PMODE_LC_SPP && popt->mode <= PMODE_LC_PPP) {
            //            if (ins_use_rtk(popt, imu->property, inskf->imu_obs, krover_idx, kbase_idx, obss, navs, ins_rtk)) {
            //                match_gnss=1;
            //                *kpva_idx=0;
            //                sol2pva(&popt->insopt,*kpva_idx,&ins_rtk->sol,pva);
            //                inskf->last_couple_time=ins_rtk->sol.time;
            //                inskf->sol->g_status=ins_rtk->sol.stat;
            //            } else match_gnss=0;
        }
    }
    else if (tc) {
        match_gnss = matchgnssobs(inskf->imu_obs->time, krover_idx, obss, 0.5 / imu->property->freq_imu);
    }

    if (popt->insopt.gps_loss_s > 0.0 && popt->insopt.gps_loss_last > 0.0 && match_gnss == 1 && inskf->couple_epoch == 0) {
        gps_loss_start = timeadd(inskf->imu_obs->time, popt->insopt.gps_loss_s);
        gps_loss_end = timeadd(gps_loss_start, popt->insopt.gps_loss_last);
    }

    if (match_gnss && inskf->couple_epoch > 0 && popt->insopt.gps_loss_s > 0.0 && popt->insopt.gps_loss_last > 0.0) {
        if (timediff(inskf->imu_obs->time, gps_loss_start) > 0.0 && timediff(inskf->imu_obs->time, gps_loss_end) < 0.0) {
            match_gnss = 0;
        }
    }
    return match_gnss;
}

extern int ig_sync_obs(const prcopt_t* popt, int* krover_idx, const obs_t* obs, const imu_t* imu, const int kimu, int dir)
{
    static int first = 1;
    int k0 = kimu, k1 = 0, i, sync = 0;
    double imut0 = 0.0, imut1 = 0.0, gpst = 0.0;
    if (popt->insopt.ms >= 2 && popt->insopt.ms <= 5) {
        k1 = k0 + popt->insopt.ms - 1;
    }
    else k1 = k0;

    k0 -= 1;
    if (first) {
        k0 = k0 + 1;
        first = 0;
    }

    if (*krover_idx<0 || *krover_idx>obs->n) return 0;

    imut0 = time2gpst(imu->data[k0].time, nullptr);
    imut1 = time2gpst(imu->data[k1].time, nullptr);
    if (dir == 0) { /*Forward*/
        for (i = *krover_idx > 100 ? *krover_idx - 100 : 0; i < obs->n; i++) {
            if (obs->data[i].rcv != 1) continue;
            gpst = time2gpst(obs->data[i].time, nullptr);
            if (gpst > imut0 && gpst <= imut1) {
                sync = 1;
#if 0
                if (k0 >= 578982) {
                    fprintf(stdout, "imu index: %6d\n", kimu + 1 - 69999);
                    fprintf(stdout, "imu tag  : %8.3f\n", imut1);
                    fprintf(stdout, "imu epoch: %6d \n", kimu);

                    fprintf(stdout, "gps index: %6d\n", i + 1);
                    fprintf(stdout, "gps tag  : %8.3f %20s\n\n", gpst, time_str(obs->data[i].time, 2));
                    fflush(stdout);

                }

#endif
                * krover_idx = i;
                break;
            }
            else if (gpst - imut1 > 1.0) {
                sync = 0;
                break;
            }
        }
    }
    else { /*Backward*/
        for (i = *krover_idx > (*krover_idx - 100) ? *krover_idx : *krover_idx - 100; i > 0; i--) {
            gpst = time2gpst(obs->data[i].time, nullptr);
            if (gpst >= imut0 && gpst < imut1) {
                sync = 1;
                *krover_idx = i;
                break;
            }
            else if (gpst - imut0 < 1.0) {
                sync = 0;
                break;
            }
        }
    }

    return sync;
}

static int ig_sync_sol(const prcopt_t* popt, int* kpva_idx, const pva_t* pva, const imu_t* imu, const int kimu, int dir)
{
    int k0 = kimu, k1 = 0, i, sync = 0;
    double imut0 = 0.0, imut1 = 0.0, gpst = 0.0;
    static int first = 1;

    if (popt->insopt.ms >= 2 && popt->insopt.ms <= 5) {
        k1 = k0 + popt->insopt.ms - 1;
    }

    k0 -= 1;
    if (first) {
        k0 = k0 + 1;
        first = 0;
    }

    if (*kpva_idx<0 || *kpva_idx>pva->n) return 0;

    imut0 = time2gpst(imu->data[k0].time, nullptr);
    imut1 = time2gpst(imu->data[k1].time, nullptr);
    if (dir == 0) { /*Forward*/
        for (i = *kpva_idx > 5 ? *kpva_idx - 5 : 0; i < pva->n; i++) {
            gpst = time2gpst(pva->time[i], nullptr);
            if (gpst > imut0 && gpst <= imut1) {
                sync = 1;
#if 0
                fprintf(stdout, "imu index: %6d\n", kimu + 1 - 69999);
                fprintf(stdout, "imu tag  : %8.3f\n", imut1);
                fprintf(stdout, "imu epoch: %6d \n", kimu);

                fprintf(stdout, "gps index: %6d\n", i + 1);
                fprintf(stdout, "gps tag  : %8.3f\n\n", gpst);
                fflush(stdout);
#endif
                * kpva_idx = i;
                break;
            }
            else if (gpst - imut1 > 1.0) {
                sync = 0;
                break;
            }
        }
    }
    else { /*Backward*/
        for (i = *kpva_idx > (*kpva_idx - 5) ? *kpva_idx : *kpva_idx - 5; i > 0; i--) {
            gpst = time2gpst(pva->time[i], nullptr);
            if (gpst >= imut0 && gpst < imut1) {
                sync = 1;
                *kpva_idx = i;
                break;
            }
            else if (gpst - imut0 < 1.0) {
                sync = 0;
                break;
            }
        }
    }

    return sync;
}

/*-------------------------------------
 *                    k0          k1
 * imu_t: ------------|------*----|----
 *                           <-dt->   Forward
 *                    <--dt-->        Backward
 * gps_t: -------------------|---------
 *                           kg
---------------------------------------*/
extern int ig_sync(const prcopt_t* popt, kf_t* inskf, int* kpva_idx, pva_t* pva, int* krover_idx, int* kbase_idx,
    const obs_t* obss, const nav_t* navs, rtk_t* ins_rtk, const imu_t* imu, const int kimu)
{
    int match_gnss = 0;
    int stc = popt->mode == PMODE_STC_PPK || popt->mode == PMODE_STC_PPP;
    static gtime_t gps_loss_start = { 0 }, gps_loss_end = { 0 };

    if (popt->mode == PMODE_LC_POS) {
        match_gnss = ig_sync_sol(popt, kpva_idx, pva, imu, kimu, popt->insopt.back);
    }
    else if ((popt->mode >= PMODE_LC_SPP && popt->mode <= PMODE_LC_PPP)) {
        if (ins_use_rtk(popt, imu->property, imu, kimu, krover_idx, kbase_idx, obss, navs, ins_rtk)) {
            match_gnss = 1;
            *kpva_idx = 0;

            sol2pva(&popt->insopt, *kpva_idx, &ins_rtk->sol, pva);
            inskf->last_couple_time = ins_rtk->sol.time;
            if (ins_rtk->sol.stat == SOLQ_FIX) {
                inskf->insstate->g_status = GNSS_STATUS_PPK_AR;
            }
            else if (ins_rtk->sol.stat == SOLQ_FLOAT) {
                inskf->insstate->g_status = GNSS_STATUS_PPK;
            }
            else if (ins_rtk->sol.stat == SOLQ_PPP) {
                inskf->insstate->g_status = GNSS_STATUS_PPP;
            }
            else if (ins_rtk->sol.stat == SOLQ_PPP_AR) {
                inskf->insstate->g_status = GNSS_STATUS_PPP_AR;
            }
        }
        else match_gnss = 0;
    }
    else if ((popt->mode >= PMODE_TC_SPP && popt->mode <= PMODE_TC_PPP) || stc) {
        match_gnss = ig_sync_obs(popt, krover_idx, obss, imu, kimu, popt->insopt.back);
    }
    else return 0;

    if (popt->insopt.gps_loss_s > 0.0 && popt->insopt.gps_loss_last > 0.0 && match_gnss == 1 && inskf->couple_epoch == 0) {
        gps_loss_start = timeadd(imu->property->tstart, popt->insopt.gps_loss_s);
        gps_loss_end = timeadd(gps_loss_start, popt->insopt.gps_loss_last);
    }

    if (match_gnss && inskf->couple_epoch > 0 && popt->insopt.gps_loss_s > 0.0 && popt->insopt.gps_loss_last > 0.0) {
        gtime_t gt = { 0 };
        if (popt->mode >= PMODE_LC_POS && popt->mode <= PMODE_LC_PPP) {
            gt = pva->time[*kpva_idx];
        }
        else if (popt->mode >= PMODE_TC_SPP) {
            gt = obss->data[*krover_idx].time;
        }
        if (timediff(gt, gps_loss_start) > 0.0 && timediff(gt, gps_loss_end) < 0.0) {
            ins_rtk->nfix = 0;
            match_gnss = 0;
        }
    }

    /*ÅÐ¶ÏÊÇ·ñGNSSÈ±Ê§*/
    double dt;
    dt = timediff(inskf->sol->time, inskf->last_couple_time);
    if (dt > popt->sample * (inskf->gnss_loss_flag + 1) + 0.01) {
        inskf->gnss_loss_flag++;
        ins_rtk->nfix = 0;
        trace(2, "%s(%d): GNSS signal loss (%d)\n", time_str(inskf->sol->time, 1), inskf->couple_epoch, inskf->gnss_loss_flag);
        if (ins_rtk->tc) {
            for (int i = 0; i < MAXSAT; i++) {
                for (int f = 0; f < NFREQ; f++) {
                    ins_rtk->ssat[i].lock[f] = 0;
                    ins_rtk->ssat[i].outc[f]++;
                }
            }
        }
    }
    if (dt < popt->sample) {
        inskf->gnss_loss_flag = 0;
    }
    return match_gnss;
}