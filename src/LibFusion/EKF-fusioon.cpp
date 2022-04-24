/*-----------------------------------------------------------------------------
* EKF-fusion.cc : Extended kalman filter of multi sensor fusion
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/22.
*----------------------------------------------------------------------------*/
#include "rtklib.h"

#define PMIN_RP        0.2*D2R
#define PMIN_YAW       1.0*D2R
#define PMIN_VEL       0.01
#define PMIN_POS_H       0.01
#define PMIN_POS_V       0.05
#define PMIN_GYR_BIA   0.1*DPH2RPS
#define PMIN_ACC_BIA   500*UG2MPS2
#define PMIN_ARMGPS    0.01
#define PMIN_IGDT      0.001

// form F mat in ECEF
extern void etm_ECEF(kf_t* inskf, const imup_t* imup, const insopt_t* opt, double dt, const v3_t dv, const v3_t dtheta)
{
    m3_t m3F;
    m3_t Cbe = inskf->insstate->dcm;
    jacobi_trans_v2a_ECEF(&m3F, &Cbe, &dv);
    m3_paste(inskf->F + xiV(opt) + xiA(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_v2p_ECEF(&m3F, dt, &inskf->insstate->pos, opt->local_coord);
    m3_paste(inskf->F + xiV(opt) + xiP(opt) * inskf->nx, inskf->nx, &m3F);
    if (opt->est_ba) {
        jacobi_trans_v2ba_ECEF(&m3F, dt, &Cbe);
        m3_paste(inskf->F + xiV(opt) + xiBa(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_bg) {
        jacobi_trans_a2bg_ECEF(&m3F, dt, &Cbe);
        m3_paste(inskf->F + xiA(opt) + xiBg(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_sa) {  // F26
        jacobi_trans_v2sa_ECEF(&m3F, &Cbe, &dv);
        m3_paste(inskf->F + xiV(opt) + xiSa(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_sg) {  // F37
        jacobi_trans_a2sg_ECEF(&m3F, &Cbe, &dtheta);
        m3_paste(inskf->F + xiA(opt) + xiSg(opt) * inskf->nx, inskf->nx, &m3F);
    }
}

// form F mat in BLH
extern void etm_LLH_NED(kf_t* inskf, const imup_t* imup, const insopt_t* opt, double dt, const v3_t dv, const v3_t dtheta)
{
    /* Local variable definition ------------------------------------------- */
    int upd_nx = xnX(opt);
    m3_t m3F;
    m3_t Cbn = inskf->insstate->dcm;

    /* Code start ---------------------------------------------------------- */
    jacobi_trans_p2p_NED(&m3F, inskf->insstate, dt);
    m3_paste(inskf->F + xiP(opt) + xiP(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_p2v_NED(&m3F, dt);
    m3_paste(inskf->F + xiP(opt) + xiV(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_v2p_NED(&m3F, inskf->insstate, dt);
    m3_paste(inskf->F + xiV(opt) + xiP(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_v2v_NED(&m3F, inskf->insstate, dt);
    m3_paste(inskf->F + xiV(opt) + xiV(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_v2a_NED(&m3F, &Cbn, &dv);
    m3_paste(inskf->F + xiV(opt) + xiA(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_a2p_NED(&m3F, inskf->insstate, dt);
    m3_paste(inskf->F + xiA(opt) + xiP(opt) * inskf->nx, inskf->nx, &m3F);
    jacobi_trans_a2v_NED(&m3F, inskf->insstate, dt);
    m3_paste(inskf->F + xiA(opt) + xiV(opt) * inskf->nx, inskf->nx, &m3F);

    if (opt->est_ba) {
        jacobi_trans_v2ba_NED(&m3F, &Cbn, dt);
        m3_paste(inskf->F + xiV(opt) + xiBa(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_bg) {
        jacobi_trans_a2bg_NED(&m3F, &Cbn, dt);
        m3_paste(inskf->F + xiA(opt) + xiBg(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_sa) {
        jacobi_trans_v2sa_NED(&m3F, &Cbn, &dv);
        m3_paste(inskf->F + xiV(opt) + xiSa(opt) * inskf->nx, inskf->nx, &m3F);
    }
    if (opt->est_sg) {
        jacobi_trans_a2sg_NED(&m3F, &Cbn, &dtheta);
        m3_paste(inskf->F + xiA(opt) + xiSg(opt) * inskf->nx, inskf->nx, &m3F);
    }
}

EXPORT void etm_LLH_ENU(kf_t* inskf, const imup_t* imup, const insopt_t* opt, double dt, const v3_t dv, const v3_t dtheta)
{
    int upd_nx = inskf->nx;
    double tl = inskf->insstate->eth.tl;
    double secl = 1.0 / inskf->insstate->eth.cl;
    double f_RMh = 1.0 / inskf->insstate->eth.RMh;
    double f_RNh = 1.0 / inskf->insstate->eth.RNh;
    double f_clRNh = 1.0 / inskf->insstate->eth.clRNh;
    double f_RMh2 = f_RMh * f_RMh;
    double f_RNh2 = f_RNh * f_RNh;
    v3_t vn = inskf->insstate->vel;
    double vE_clRNh = vn.x * f_clRNh;
    double vE_RNh2 = vn.x * f_RNh2;
    double vN_RMh2 = vn.y * f_RMh2;
    m3_t Mp1 = O33;
    Mp1.m21 = -inskf->insstate->eth.wnie.z;
    Mp1.m31 = inskf->insstate->eth.wnie.y;

    m3_t Mp2 = O33;
    Mp2.m31 = vE_clRNh * secl;
    Mp2.m13 = vN_RMh2;
    Mp2.m23 = -vE_RNh2;
    Mp2.m33 = -vE_RNh2 * tl;

    m3_t Avn = O33;
    Avn = v3_askew(vn);

    m3_t Awn = O33;
    Awn = v3_askew(inskf->insstate->eth.wnien);

    m3_t Maa = O33;
    Maa = v3_askew(inskf->insstate->eth.wnin);

    m3_t Mav = O33;
    Mav.m21 = f_RNh;
    Mav.m31 = f_RNh * tl;
    Mav.m12 = -f_RMh;

    m3_t Map = m3_add(Mp1, Mp2);

    m3_t Mva = v3_askew(inskf->insstate->fn);
    m3_t Mvv = m3_del(m3_mul(Avn, Mav), Awn);
    m3_t Mvp = m3_mul(Avn, m3_add(Mp1, Map));
    double scl = inskf->insstate->eth.sl * inskf->insstate->eth.cl;
    Mvp.m31 -= G2MPS2 * (5.27094e-3 * 2 + 2.32718e-5 * 4 * inskf->insstate->eth.sl2) * scl;
    Mvp.m33 += 3.086e-6;

    m3_t Mpv = inskf->insstate->eth.Mpv;
    m3_t Mpp = O33;
    Mpp.m21 = vE_clRNh * tl;
    Mpp.m13 = -vN_RMh2;
    Mpp.m23 = -vE_RNh2 * secl;

    Maa = m3_scalar(-dt, Maa);
    Maa = m3_diag_add(Maa, 1.0);
    Mav = m3_scalar(dt, Mav);
    Map = m3_scalar(dt, Map);
    Mva = m3_scalar(dt, Mva);
    Mvv = m3_scalar(dt, Mvv);
    Mvv = m3_diag_add(Mvv, 1.0);
    Mvp = m3_scalar(dt, Mvp);
    Mpv = m3_scalar(dt, Mpv);
    Mpp = m3_scalar(dt, Mpp);
    Mpp = m3_diag_add(Mpp, 1.0);

    m3_paste(inskf->F + xiA(opt) + xiA(opt) * upd_nx, upd_nx, &Maa);
    m3_paste(inskf->F + xiA(opt) + xiV(opt) * upd_nx, upd_nx, &Mav);
    m3_paste(inskf->F + xiA(opt) + xiP(opt) * upd_nx, upd_nx, &Map);

    m3_paste(inskf->F + xiV(opt) + xiA(opt) * upd_nx, upd_nx, &Mva);
    m3_paste(inskf->F + xiV(opt) + xiV(opt) * upd_nx, upd_nx, &Mvv);
    m3_paste(inskf->F + xiV(opt) + xiP(opt) * upd_nx, upd_nx, &Mvp);

    m3_paste(inskf->F + xiP(opt) + xiV(opt) * upd_nx, upd_nx, &Mpv);
    m3_paste(inskf->F + xiP(opt) + xiP(opt) * upd_nx, upd_nx, &Mpp);
    if (opt->est_bg) {
        m3_t Cbn = m3_T(m3_scalar(-dt, inskf->insstate->dcm));
        m3_paste(inskf->F + xiA(opt) + xiBg(opt) * upd_nx, upd_nx, &Cbn);
        m3_t m3F;
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Tg);
        m3_paste(inskf->F + xiBg(opt) + xiBg(opt) * upd_nx, upd_nx, &m3F);
    }
    if (opt->est_ba) {
        m3_t Cbn = m3_T(m3_scalar(dt, inskf->insstate->dcm));
        m3_paste(inskf->F + xiV(opt) + xiBa(opt) * upd_nx, upd_nx, &Cbn);
        m3_t m3F;
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Ta);
        m3_paste(inskf->F + xiBa(opt) + xiBa(opt) * upd_nx, upd_nx, &m3F);
    }
}

/**
 * @brief inskf_norm_innov
 * @param[in]   inskf ins kalman filter struct
 * @param[in]   dz    innovation(measment - model, inskf->ny x 1)
 * @param[in]   R     measment noise matrix(inskf->ny x inskf->ny)
 * @param[out]  ndz   ouput normialized innovations(inskf->ny x 1)
 * @return  0: OK
 */
extern int inskf_norm_innov(kf_t* inskf, const double* dz, const double* R, double* ndz) {
    /* H'*PH */
    double* tmpm = zeros(inskf->ny, inskf->nx);
    double* Rm = zeros(inskf->ny, inskf->ny);
    matmul("TN", inskf->ny, inskf->nx, inskf->nx, 1.0, inskf->H, inskf->P, 0.0, tmpm);
    matmul("NN", inskf->ny, inskf->ny, inskf->nx, 1.0, tmpm, inskf->H, 0.0, Rm);

    for (int i = 0; i < inskf->ny; ++i)
        ndz[i] = dz[i] / sqrt(Rm[i + i * inskf->ny] + R[i + i * inskf->ny]);

    free(tmpm);
    free(Rm);
    return 0;
}

/**
 * @brief form transfer noise matrix of kalman matrix
 *          P = F * P * F' + Q
 * @param[out]  Q           transer noise matrix
 * @param[in]   Cbe         Attitude transformation from b to e
 * @param[in]   dt_zoom     time interval[s] x zoom
 * @param[in]   imup        imu property struct
 * @return status(0: OK)
 */
extern int kf_formQ(double* Q, const m3_t* Cbe, double dt_zoom, const imup_t* imup, const insopt_t* opt) {
    int upd_nx = 0;
    upd_nx = xnX(opt);

    for (int i = 0; i < upd_nx * upd_nx; i++)
        Q[i] = 0.0;

    /* Angular white noise */
    m3_t Qphi = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->gyrnd, 2.0)));
    if (!opt->is_imu_samenoise && Cbe != NULL)
        Qphi = m3_mul(m3_mul(*Cbe, Qphi), m3_T(*Cbe));
    m3_paste((Q + xiA(opt) + xiA(opt) * upd_nx), upd_nx, &Qphi);

    /* Velocity white noise */
    m3_t Qv = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->accnd, 2.0)));
    if (!opt->is_imu_samenoise && Cbe != NULL)
        Qv = m3_mul(m3_mul(*Cbe, Qv), m3_T(*Cbe));
    m3_paste((Q + xiV(opt) + xiV(opt) * upd_nx), upd_nx, &Qv);

    if (opt->est_ba) {
        /* Accel bias random walk */
        m3_t Qba = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->abrw, 2.0)));
        m3_paste((Q + xiBa(opt) + xiBa(opt) * upd_nx), upd_nx, &Qba);
    }
    if (opt->est_bg) {
        /* Gyro bias random walk */
        m3_t Qbg = v3_diag(v3_scalar(dt_zoom, v3_pow(imup->gbrw, 2.0)));
        m3_paste((Q + xiBg(opt) + xiBg(opt) * upd_nx), upd_nx, &Qbg);
    }

    return 0;
}

extern void inskf_proP(double* x, double* P, kf_t* inskf, int all_nx, int upd_nx, int tra) {
    double* P_, * Pp_, * F_, * Q_;
    int i, j, k, * ix;

    ix = imat(upd_nx, 1);
    for (i = k = 0; i < upd_nx; i++) if (x[i] != 0.0 && P[i + i * all_nx] > 0.0) ix[k++] = i;
    P_ = mat(k, k);
    Pp_ = mat(k, k);
    F_ = mat(k, k);
    Q_ = mat(k, k);

    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) P_[i + j * k] = P[ix[i] + ix[j] * all_nx];
        for (j = 0; j < k; j++) F_[i + j * k] = inskf->F[ix[i] + ix[j] * all_nx];
        for (j = 0; j < k; j++) Q_[i + j * k] = inskf->Q[ix[i] + ix[j] * all_nx];
    }

    if (tra) {
        matprint(0, F_, k, k, 20, 10);
        matprint(1, P_, k, k, 20, 10);
        matprint(1, Q_, k, k, 20, 10);
    }

    double* FP = mat(k, k);
    matmul("NN", k, k, k, 1.0, F_, P_, 0.0, FP);
    matmul("NT", k, k, k, 1.0, FP, F_, 0.0, Pp_);

    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            Pp_[i + j * k] += Q_[i + j * k];
        }
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) P[ix[i] + ix[j] * all_nx] = Pp_[i + j * k];
    }

    free(P_);
    free(Pp_);
    free(F_);
    free(FP);
    free(Q_);
    free(ix);
}

static void init_inssol(solins_t* solins)
{
    gtime_t t0 = { 0 };
    solins->time = t0;
    solins->status = SOLQ_INS_ALIGN;
    solins->an = V0;
    solins->ba = V0;
    solins->bg = V0;
    solins->sa = V0;
    solins->sg = V0;
    solins->ba_std = V0;
    solins->bg_std = V0;
    solins->sa_std = V0;
    solins->sg_std = V0;
    solins->arm_gps = V0;
    solins->arm_gps_std = V0;
    solins->delay_pos = V0;
    solins->t_delay = 0.0;
}

extern int inskf_init(kf_t* inskf, const imup_t* imup, insopt_t* opt) {
    int upd_nx = 0;
    auto* p = (prcopt_t*)opt->gnss_opt;

    inskf->couple_epoch = 0;
    inskf->ins_epoch = 0;
    inskf->imup = const_cast<imup_t*>(imup);
    inskf->time = imup->tstart;
    inskf->gnss_loss_flag = 0;

    inskf->nix = xnIns(opt);
    if (opt->gnss_opt && ((prcopt_t*)opt->gnss_opt)->mode >= PMODE_TC_SPP &&
        ((prcopt_t*)opt->gnss_opt)->mode <= PMODE_TC_PPP)
        inskf->ngx = xnGnss(opt);
    inskf->nx = xnX(opt);
    inskf->na = xnR(opt);

    upd_nx = inskf->nx;
    if (imup->freq_imu > 0) {
        int N = (opt->ms == -2 || opt->ms == 0 || opt->ms > 5) ? 1 : opt->ms;
        inskf->idt = (1.0 / imup->freq_imu) * N;
    }
    if (imup->freq_od > 0) {
        inskf->odt = 1.0 / (double)imup->freq_od;
    }
    else inskf->odt = 0.0;

    if (!(inskf->x = (double*)malloc(inskf->nx * sizeof(double)))) {
        trace(1, "inskf->x memory allocation error\n");
    }
    if (!(inskf->insstate = (solins_t*)malloc(sizeof(solins_t)))) {
        trace(1, "inskf->sol memory allocation error\n");
    }
    if (!(inskf->sol = (solins_t*)malloc(sizeof(solins_t)))) {
        trace(1, "inskf->sol memory allocation error\n");
    }
    if (!(inskf->P = (double*)calloc(inskf->nx * inskf->nx, sizeof(double)))) {
        trace(1, "inskf->P memory allocation error\n");
    }
    if (!(inskf->F = (double*)calloc(upd_nx * upd_nx, sizeof(double)))) {
        trace(1, "inskf->F memory allocation error\n");
    }
    if (!(inskf->Q = (double*)calloc(upd_nx * upd_nx, sizeof(double)))) {
        trace(1, "inskf->Q memory allocation error\n");
    }

    if (p && p->mode >= PMODE_TC_SPP && p->mode <= PMODE_TC_PPP) {
        if (!(inskf->xa = (double*)calloc(inskf->na, sizeof(double)))) {
            trace(1, "inskf->Xa memory allocation error\n");
        }
        if (!(inskf->Pa = (double*)calloc(inskf->na * inskf->na, sizeof(double)))) {
            trace(1, "inskf->Pa memory allocation error\n");
        }
    }

    int i, j;
    /* init x */
    for (i = 0; i < inskf->nx; i++) inskf->x[i] = i < inskf->nix ? 1E-32 : 0.0;
    for (i = 0; i < upd_nx; i++) {
        for (j = 0; j < upd_nx; j++) {
            inskf->F[i + j * upd_nx] = i == j ? 1.0 : 0.0;
        }
    }

    /* init sol*/
    init_inssol(inskf->insstate); init_inssol(inskf->sol);
    inskf->sol->g_status = inskf->insstate->g_status = GNSS_STATUS_NONE;
    inskf->sol->time = inskf->insstate->time = imup->init_tag;
    v3_t pos = imup->initr;
    if (opt->mech_coord == INSMECH_ECEF) {
        ecef2llh(&pos, nullptr, nullptr, opt->local_coord);
    }

    /*init ins*/
    inskf->insstate->rpy = imup->inita;
    inskf->insstate->ns = 0;
    if (opt->mech_coord == INSMECH_ECEF) {
        inskf->insstate->dcm = att2Cbe(&pos, &imup->inita, opt->imu_coord, opt->local_coord);
        dcm2quat(&inskf->insstate->dcm, &inskf->insstate->quat);
    }
    else if (opt->mech_coord == INSMECH_LLH) {
        v3_t rpy = imup->inita;
        attsync(&inskf->insstate->pos, &rpy, &inskf->insstate->dcm, &inskf->insstate->quat, opt->imu_coord, opt->local_coord, opt->mech_coord, 0);
    }
    inskf->insstate->vel = imup->initv;
    inskf->insstate->pos = imup->initr;
    if (opt->mech_coord == INSMECH_LLH) {
        inskf->insstate->eth = wgs84;
        updearth(&inskf->insstate->eth, inskf->insstate->pos, inskf->insstate->vel, opt->local_coord, 1);
    }

    /* init P */
    m3_t Qr, Qv, Qa;
    Qr = v3_diag(v3_pow(imup->initr_err, 2));
    Qv = v3_diag(v3_pow(imup->initv_err, 2));
    /*TODO: check*/
    Qa = v3_diag(v3_pow(imup->inita_err, 2));
    if (opt->mech_coord == INSMECH_ECEF) {   // note here
        llh2ecefQ(&pos, nullptr, nullptr, &Qa, opt->local_coord);
    }
    m3_paste(inskf->P + xiA(opt) + xiA(opt) * inskf->nx, inskf->nx, &Qa);
    m3_paste(inskf->P + xiV(opt) + xiV(opt) * inskf->nx, inskf->nx, &Qv);
    m3_paste(inskf->P + xiP(opt) + xiP(opt) * inskf->nx, inskf->nx, &Qr);

    /* init F */
    m3_t m3F;
    if (opt->mech_coord == INSMECH_ECEF) {
        jacobi_trans_a2a_ECEF(&m3F, inskf->idt);
        m3_paste(inskf->F + xiA(opt) + xiA(opt) * upd_nx, upd_nx, &m3F);
        jacobi_trans_v2v_ECEF(&m3F, inskf->idt);
        m3_paste(inskf->F + xiV(opt) + xiV(opt) * upd_nx, upd_nx, &m3F);
        jacobi_trans_p2v_ECEF(&m3F, inskf->idt);
        m3_paste(inskf->F + xiP(opt) + xiV(opt) * upd_nx, upd_nx, &m3F);
        for (i = 0; i < xiP(opt) + xnP(opt); i++) {
            inskf->F[i + i * upd_nx] = 1.0;
        }
    }

    if (opt->est_ba) {
        inskf->insstate->ba = V0;
        m3_t initQba = v3_diag(v3_pow(imup->ba_err, 2.0));
        m3_paste(inskf->P + xiBa(opt) + xiBa(opt) * inskf->nx, inskf->nx, &initQba);
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Ta); /*I+F*dt*/
        m3_paste(inskf->F + xiBa(opt) + xiBa(opt) * upd_nx, upd_nx, &m3F);
    }
    if (opt->est_bg) {
        inskf->insstate->bg = V0;
        m3_t initQbg = v3_diag(v3_pow(imup->bg_err, 2.0));
        m3_paste(inskf->P + xiBg(opt) + xiBg(opt) * inskf->nx, inskf->nx, &initQbg);
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Tg);
        m3_paste(inskf->F + xiBg(opt) + xiBg(opt) * upd_nx, upd_nx, &m3F);
    }
    if (opt->est_sa) {
        inskf->insstate->sa = V1;
        m3_t initQba = v3_diag(v3_pow(imup->sa_err, 2.0));
        m3_paste(inskf->P + xiSa(opt) + xiSa(opt) * inskf->nx, inskf->nx, &initQba);
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Ta);
        m3_paste(inskf->F + xiSa(opt) + xiSa(opt) * upd_nx, upd_nx, &m3F);
    }
    if (opt->est_sg) {
        inskf->insstate->sg = V1;
        m3_t initQbg = v3_diag(v3_pow(imup->sg_err, 2.0));
        m3_paste(inskf->P + xiSg(opt) + xiSg(opt) * inskf->nx, inskf->nx, &initQbg);
        jacobi_trans_markov(&m3F, inskf->idt, &imup->Tg);
        m3_paste(inskf->F + xiSg(opt) + xiSg(opt) * upd_nx, upd_nx, &m3F);
    }

    if (opt->est_armgps) {
        inskf->insstate->arm_gps = V0;
    }
    else {
        inskf->insstate->arm_gps = imup->lever_arm_gps;
    }
    if (opt->est_armgps) {
        v3_t arm_std = imup->lever_arm_gps_std;
        if (norm(imup->lever_arm_gps_std.v, 3) == 0.0) {
            arm_std = V1;
        }
        m3_t initQarm = v3_diag(v3_pow(arm_std, 2.0));
        m3_paste(inskf->P + xiArmGps(opt) + xiArmGps(opt) * inskf->nx, inskf->nx, &initQarm);
    }
    if (opt->est_igdt) {
        inskf->P[xiDt(opt) + xiDt(opt) * inskf->nx] = SQR(0.1);
    }

    if (imup->kod != 0.0) {
        inskf->insstate->kod = imup->kod;
    }
    else {
        inskf->insstate->kod = 1.0;
    }

    if (opt->est_kod) {
        inskf->F[xiKod(opt) + xiKod(opt) * upd_nx] = 1.0;
        inskf->P[xiKod(opt) + xiKod(opt) * inskf->nx] = SQR(imup->kod_err);
    }

    inskf->ny = 6;
    if (!(inskf->H = (double*)calloc(inskf->ny * inskf->nx, sizeof(double)))) {
        trace(1, "inskf->H memory allocation error\n");
    }

    if (!(inskf->R = (double*)calloc(inskf->ny, sizeof(double)))) {
        trace(1, "inskf->R memory allocation error\n");
    }

    if (opt->zvopt.ws > 0) {
        if (!(inskf->imud = (imud_t*)malloc(opt->zvopt.ws * sizeof(imud_t)))) {
            trace(1, "inskf->imud memory allocation error\n");
        }
        inskf->nimud = 0;
        inskf->imudend = opt->zvopt.ws - 1;
    }

    if (!opt->is_imu_samenoise) {
        if (fabs(imup->gyrnd.x - imup->gyrnd.y) < 1e-32 &&
            fabs(imup->gyrnd.x - imup->gyrnd.z) < 1e-32 &&
            fabs(imup->accnd.x - imup->accnd.y) < 1e-32 &&
            fabs(imup->accnd.x - imup->accnd.z) < 1e-32)
            opt->is_imu_samenoise = 1;
    }
    if (opt->is_imu_samenoise) {
        kf_formQ(inskf->Q, nullptr, inskf->idt, imup, opt);
    }

    kf_copyQ2sol(inskf, opt);

#if 0
    matprint(0, inskf->F, inskf->nx, inskf->nx, 20, 10);
    matprint(0, inskf->P, inskf->nx, inskf->nx, 20, 10);
    matprint(0, inskf->Q, inskf->nx, inskf->nx, 20, 10);
#endif

    return 1;
}

static void inskf_P_constrain(const insopt_t* opt, kf_t* inskf)
{
    int i, nx = inskf->nx;
    /*Position constrain*/
    if (opt->mech_coord == INSMECH_ECEF) {
        //        for(i=xiP(opt);i<xiP(opt)+xnP(opt);i++){
        //            if(inskf->P[i+i*nx]<=SQR(PMIN_POS)){
        //                inskf->P[i+i*nx]=SQR(PMIN_POS);
        //            }
        //        }
        if (inskf->P[xiP(opt) + xiP(opt) * nx] <= SQR(PMIN_POS_H)) {
            inskf->P[xiP(opt) + xiP(opt) * nx] = SQR(PMIN_POS_H);
        }
        if (inskf->P[xiP(opt) + 1 + (xiP(opt) + 1) * nx] <= SQR(PMIN_POS_H)) {
            inskf->P[xiP(opt) + 1 + (xiP(opt) + 1) * nx] = SQR(PMIN_POS_H);
        }
        if (inskf->P[xiP(opt) + 2 + (xiP(opt) + 2) * nx] <= SQR(PMIN_POS_V)) {
            inskf->P[xiP(opt) + 2 + (xiP(opt) + 2) * nx] = SQR(PMIN_POS_V);
        }
    }
    else if (opt->mech_coord == INSMECH_LLH) {
        /*lat*/
        inskf->P[xiP(opt) + xiP(opt) * nx] = SQR(PMIN_POS_H / wgs84.R0);
        /*lon*/
        inskf->P[(xiP(opt) + 1) + (xiP(opt) + 1) * nx] = SQR(PMIN_POS_H / wgs84.R0);
        /*hgt*/
        inskf->P[(xiP(opt) + 2) + (xiP(opt) + 2) * nx] = SQR(PMIN_POS_V);
    }
    /*Velocity constrain*/
//    for(i=xiV(opt);i<xiV(opt)+xnV(opt);i++){
//        if(inskf->P[i+i*nx]<=SQR(PMIN_VEL)){
//            inskf->P[i+i*nx]=SQR(PMIN_VEL);
//        }
//    }
//    /*Attitude roll and pitch constrain*/
//    if(inskf->P[xiA(opt)+xiA(opt)*nx]<=SQR(PMIN_RP)){
//        inskf->P[xiA(opt)+xiA(opt)*nx]=PMIN_RP;
//    }
//    if(inskf->P[xiA(opt)+1+(xiA(opt)+1)*nx]<=SQR(PMIN_RP)){
//        inskf->P[xiA(opt)+1+(xiA(opt)+1)*nx]=PMIN_RP;
//    }
//    /*Attitude yaw constrain*/
//    if(inskf->P[xiA(opt)+2+(xiA(opt)+2)*nx]<=SQR(PMIN_YAW)){
//        inskf->P[xiA(opt)+2+(xiA(opt)+2)*nx]=PMIN_YAW;
//    }
    if (opt->est_bg) {
        /*gyro bias constrain*/
        for (i = xiBg(opt); i < xiBg(opt) + xnBg(opt); i++) {
            if (inskf->P[i + i * nx] <= SQR(PMIN_GYR_BIA)) {
                inskf->P[i + i * nx] = SQR(PMIN_GYR_BIA);
            }
        }
    }
    if (opt->est_ba) {
        /*acc bias constrain*/
        for (i = xiBa(opt); i < xiBa(opt) + xnBa(opt); i++) {
            if (inskf->P[i + i * nx] <= SQR(PMIN_ACC_BIA)) {
                inskf->P[i + i * nx] = SQR(PMIN_ACC_BIA);
            }
        }
    }
    if (opt->est_armgps) {
        /*arm gps*/
        for (i = xiArmGps(opt); i < xiArmGps(opt) + xnArmGps(opt); i++) {
            if (inskf->P[i + i * nx] <= SQR(PMIN_ARMGPS)) {
                inskf->P[i + i * nx] = SQR(PMIN_ARMGPS);
            }
        }
    }
    if (opt->est_igdt) {
        /*igdt*/
        if (inskf->P[xiDt(opt) + xiDt(opt) * nx] <= SQR(PMIN_IGDT)) {
            inskf->P[xiDt(opt) + xiDt(opt) * nx] = SQR(PMIN_IGDT);
        }
    }
}

/* time update include ins mechanization */
extern int inskf_udstate(kf_t* inskf, const imu_t* imus, int* imu_idx, const imup_t* imup, const insopt_t* opt) {
    int upd_nx = inskf->nx, nsample = 1, tra = 0;
    gtime_t cur_time;
    v3_t dv, dtheta;

    /* conning&sculling error compensation */
    int N = 2;
    N = (opt->ms == -2 || opt->ms == 0 || opt->ms > 5) ? 2 : opt->ms;
    if (opt->ms >= 2 && opt->ms <= 5) inskf->nsample = opt->ms;
    else inskf->nsample = 1;

    v3_t* dtheta_list;
    v3_t* dv_list;
    dtheta_list = (v3_t*)malloc(N * sizeof(v3_t));
    dv_list = (v3_t*)malloc(N * sizeof(v3_t));
    if (opt->ms == -2) {/*one-plus-previous*/
        dtheta_list[0] = inskf->dthetap;
        dtheta_list[1] = imus->data[*imu_idx].gyro;
        dv_list[0] = inskf->dvp;
        dv_list[1] = imus->data[*imu_idx].accel;
        cur_time = imus->data[*imu_idx].time;
        if (multisample(dtheta_list, dv_list, opt->ms, &dtheta, &dv)) {
            trace(3, "failed to conning&sculling error compensation\n");
        }
        if (opt->ms == -2) {
            inskf->dthetap = dtheta;
            inskf->dvp = dv;
        }
    }
    else if (opt->ms >= 2 && opt->ms <= 5) {/*1-5 sample*/
        for (int i = 0; i < N; i++) {
            dtheta_list[i] = imus->data[*imu_idx + i].gyro;
            dv_list[i] = imus->data[*imu_idx + i].accel;
        }
        cur_time = imus->data[*imu_idx + N - 1].time;
        if (multisample(dtheta_list, dv_list, opt->ms, &dtheta, &dv)) {
            trace(3, "failed to conning&sculling error compensation\n");
        }
    }
    else { /*no*/
        cur_time = imus->data[*imu_idx].time;
        dv = imus->data[*imu_idx].accel;
        dtheta = imus->data[*imu_idx].gyro;
    }

    inskf->omgb = dtheta;
    inskf->fb = dv;
    double dt = inskf->idt;

    /* corrent imu output, true = (output - bias)*factor */
    if (opt->est_sa) { dv = v3_mul_cxc(inskf->insstate->sa, dv); }
    if (opt->est_sg) { dtheta = v3_mul_cxc(inskf->insstate->sg, dtheta); }
    dv = v3_del(dv, v3_scalar(dt, inskf->insstate->ba));
    dtheta = v3_del(dtheta, v3_scalar(dt, inskf->insstate->bg));
    inskf->time = cur_time;

    int mech_stat = 0;
    if (opt->mech_coord == INSMECH_ECEF) {
        if (opt->back) {
            mech_stat = ins_nav_ecef_back(dt, &dtheta, &dv, &inskf->insstate->pos, &inskf->insstate->vel, &inskf->insstate->quat, opt->local_coord);
        }
        else {
            mech_stat = ins_nav_ecef(dt, &dtheta, &dv, &inskf->insstate->pos, &inskf->insstate->vel, &inskf->insstate->quat, opt->local_coord);
        }
    }
    else if (opt->mech_coord == INSMECH_LLH) {
        if (opt->back) {
            mech_stat = ins_nav_llh_back(dt, &dtheta, &dv, inskf->insstate, opt->local_coord);
        }
        else {
            mech_stat = ins_nav_llh(dt, &dtheta, &dv, inskf->insstate, opt->local_coord);
        }
    }

    if (mech_stat == 1) {
        inskf->insstate->time = inskf->time;
        inskf->insstate->status = SOLQ_INS;
        quat2dcm(&inskf->insstate->quat, &inskf->insstate->dcm, opt->imu_coord);
        //        quat2att(&inskf->sol->quat,&inskf->sol->rpy,opt->imu_coord);
        kf_copyQ2sol(inskf, opt);
    }

    /* Update F */
    if (opt->mech_coord == INSMECH_ECEF && ((prcopt_t*)opt->gnss_opt)->mode > PMODE_INS_MECH) {
        etm_ECEF(inskf, inskf->imup, opt, dt, dv, dtheta);
    }
    else if (opt->mech_coord == INSMECH_LLH && ((prcopt_t*)opt->gnss_opt)->mode > PMODE_INS_MECH) {
        if (opt->local_coord == INSLOCAL_NED) {
            etm_LLH_NED(inskf, imup, opt, dt, dv, dtheta);
        }
        else if (opt->local_coord == INSLOCAL_ENU) {
            etm_LLH_ENU(inskf, imup, opt, dt, dv, dtheta);
        }
    }

    if (inskf->couple_epoch == 0) {
        tra = 0;
    }
    if (((prcopt_t*)opt->gnss_opt)->mode > PMODE_INS_MECH) {
        inskf_proP(inskf->x, inskf->P, inskf, inskf->nx, upd_nx, tra);
#if 0
        inskf_P_constrain(opt, inskf);
#endif
    }

    inskf->insstate->zero_flag = 0;
    inskf->insstate->g_status = GNSS_STATUS_NONE;
    memcpy(inskf->sol, inskf->insstate, sizeof(solins_t));
    free(dtheta_list);
    free(dv_list);
    return 0;
}

extern int inskf_feedback(gtime_t t, unsigned int SOL_TYPE, const insopt_t* opt, double* x, solins_t* solins) {
    if (opt->feedratio < 1)
        trace(2, "do not support opt->feedratio < 1\n");

    /* attitude feedback */
    solins->time = t;
    solins->status = SOL_TYPE;

    if (norm(x + xiA(opt), 3) > 10.0 * D2R) {
        trace(2, "%s: too large residuals for attitude update\n", time_str(t, 1));
    }

    if (opt->fb_att) {
        if (opt->mech_coord == INSMECH_LLH) {
            quatdelphi(&solins->quat, (v3_t*)(x + xiA(opt)));
            quat2dcm(&solins->quat, &solins->dcm, opt->imu_coord);
            quat2att(&solins->quat, &solins->rpy, opt->imu_coord);
        }
        else if (opt->mech_coord == INSMECH_ECEF) {
            m3_t phix = v3_askew(*(v3_t*)(x + xiA(opt)));
            solins->dcm = m3_mul(m3_del(I33, phix), solins->dcm);
            dcm2quat(&solins->dcm, &solins->quat);
        }

        for (int i = xiA(opt); i < xnA(opt) + xiA(opt); ++i)
            x[i] = 1e-32;
    }

    if (opt->fb_vel) {
        solins->vel = v3_del(solins->vel, *(v3_t*)(x + xiV(opt)));
        for (int i = xiV(opt); i < xnV(opt) + xiV(opt); ++i)
            x[i] = 1e-32;
    }
    if (opt->fb_pos) {
        solins->pos = v3_del(solins->pos, *(v3_t*)(x + xiP(opt)));
        for (int i = xiP(opt); i < 3 + xiP(opt); ++i) x[i] = 1e-32;
    }
    double feedratio = opt->feedratio;
    if (feedratio > 1.0 || feedratio < 0.0) {
        feedratio = 1.0;
    }

    v3_t x3;
    if (opt->est_ba && opt->fb_ba) {
        x3 = v3_scalar(feedratio, *(v3_t*)(x + xiBa(opt)));
        if (opt->mech_coord == INSMECH_LLH) {
            solins->ba = v3_add(solins->ba, x3);
        }
        else if (opt->mech_coord == INSMECH_ECEF) {
            solins->ba = v3_del(solins->ba, x3);
        }
        for (int i = xiBa(opt); i < 3 + xiBa(opt); ++i)
            x[i] = feedratio == 1.0 ? 1E-32 : (1.0 - feedratio) * x[i];
    }
    if (opt->est_bg && opt->fb_bg) {
        x3 = v3_scalar(feedratio, *(v3_t*)(x + xiBg(opt)));
        if (opt->mech_coord == INSMECH_LLH) {
            solins->bg = v3_add(solins->bg, x3);
        }
        else if (opt->mech_coord == INSMECH_ECEF) {
            solins->bg = v3_del(solins->bg, x3);
        }
        for (int i = xiBg(opt); i < 3 + xiBg(opt); ++i) {
            x[i] = feedratio == 1.0 ? 1E-32 : (1.0 - feedratio) * x[i];
        }
    }
    if (opt->est_sa && opt->fb_sa) {
        x3 = v3_scalar(feedratio, *(v3_t*)(x + xiSa(opt)));
        if (opt->mech_coord == INSMECH_LLH) {
            solins->sa = v3_add(solins->sa, x3);
        }
        else if (opt->mech_coord == INSMECH_ECEF) {
            solins->sa = v3_del(solins->sa, x3);
        }
        for (int i = xiSa(opt); i < 3 + xiSa(opt); ++i)
            x[i] = feedratio == 1.0 ? 1E-32 : (1.0 - feedratio) * x[i];
    }
    if (opt->est_sg && opt->fb_sg) {
        x3 = v3_scalar(feedratio, *(v3_t*)(x + xiSg(opt)));
        if (opt->mech_coord == INSMECH_LLH) {
            solins->sg = v3_add(solins->sg, x3);
        }
        else if (opt->mech_coord == INSMECH_ECEF) {
            solins->sg = v3_del(solins->sg, x3);
        }
        for (int i = xiSg(opt); i < 3 + xiSg(opt); ++i)
            x[i] = feedratio == 1.0 ? 1E-32 : (1.0 - feedratio) * x[i];
    }
    if (opt->est_armgps && opt->fb_armgps) {
        v3_t arm = *(v3_t*)(x + xiArmGps(opt));
        solins->arm_gps = v3_add(solins->arm_gps, arm);
        for (int i = xiArmGps(opt); i < xiArmGps(opt) + xnArmGps(opt); ++i) {
            x[i] = 1E-32;
        }
    }
    if (opt->est_igdt && opt->fb_igdt) {
        double dt = x[xiDt(opt)];
        solins->t_delay += dt;
        x[xiDt(opt)] = 1E-32;
    }
    if (opt->est_kod) {
        solins->kod -= x[xiKod(opt)];
        x[xiKod(opt)] = 1e-32;
    }

    return 0;
}
