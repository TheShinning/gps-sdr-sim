/*-----------------------------------------------------------------------------
* ins-tc-jacobi.cc : form jacobi matrix during combine solution
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/29.
*----------------------------------------------------------------------------*/
#include "rtklib.h"

/**
 * @brief Form state transformation jacobi matrix from Ebe to Ebe
 * @param[out]  F   jacobi matrix
 * @param[in]   dt  time interval[s]
 * @return O: OK
 * @note This function relate to the earth rotation rate, use wgs84 parameter.
 *
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_a2a_ECEF(m3_t* F, double dt) {
    /* I3 - OMGie_e * dt */
    v3_t omg_ie = { 0.0, 0.0, wgs84.wie };
    omg_ie = v3_scalar(-dt, omg_ie);
    asymmetric_mat(&omg_ie, F);
    F->m11 += 1;
    F->m22 += 1;
    F->m33 += 1;
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from Ebe to gryo bias
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval
 * @param[in]   Cbe Current average attitude
 * @return 0: OK
 * @note Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_a2bg_ECEF(m3_t* F, double dt, const m3_t* Cbe) {
    /* - Cbe * dt */
    *F = m3_scalar(-dt, *Cbe);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to Ebe
 * @param[out]  F   Jacobi matrix
 * @param[in]   Cbe Current average attitude.
 * @param[in]   dv  Velocity increment[m/s]
 * @return 0: OK
 * @note Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_v2a_ECEF(m3_t* F, const m3_t* Cbe, const v3_t* dv) {
    /* [- Cbe * fib_b ]x * dt */
    v3_t Tib_e = v3_scalar(-1.0, m3_mul_v3(*Cbe, *dv));
    asymmetric_mat(&Tib_e, F);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to veb_e
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval[s]
 * @return 0: OK
 * @note Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_v2v_ECEF(m3_t* F, double dt) {
    /* I3 - 2 * OMGie_e * dt */
    v3_t omg_ie = { 0.0, 0.0, wgs84.wie };
    omg_ie = v3_scalar(-2 * dt, omg_ie);
    asymmetric_mat(&omg_ie, F);
    F->m11 += 1;
    F->m22 += 1;
    F->m33 += 1;
    return 0;
}

/**
 * @brief Form state transforamtion jacobi matrix from veb_e to reb_e
 * @param[out]  F       Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   reb_e   ecef position[m]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_v2p_ECEF(m3_t* F, double dt, const v3_t* reb_e, int opt) {
    /* F23 * dt */
    /* Ref Paul P71:2.137 */
    v3_t pos = *reb_e;
    ecef2llh(&pos, NULL, NULL, opt);
    double RE = earth_RE(&wgs84, pos.x);
    double grtmp = SQR(cos(pos.x)) + SQR(1 - SQR(wgs84.e)) * SQR(sin(pos.x));
    double geocentric_radius = RE * sqrt(grtmp);
    v3_t ge;
    gravity_ecef(reb_e, &ge);
    double fac = -dt * 2 / geocentric_radius / v3_norm(*reb_e);
    ge = v3_scalar(fac, ge);
    *F = v3_mul_cxr(ge, *reb_e);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from veb_e to ba
 * @param[out]  F       Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   Cbe     Current average attitude over time[m]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_v2ba_ECEF(m3_t* F, double dt, const m3_t* Cbe) {
    *F = m3_scalar(-dt, *Cbe);
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix from reb_e to veb_e
 * @param[out]  F   Jacobi matrix
 * @param[in]   dt  Time interval[s]
 * @return 0: OK
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P583, 14.48
 */
extern int jacobi_trans_p2v_ECEF(m3_t* F, double dt) {
    F->m11 = dt;
    F->m12 = 0.0;
    F->m13 = 0.0;
    F->m21 = 0.0;
    F->m22 = dt;
    F->m23 = 0.0;
    F->m31 = 0.0;
    F->m32 = 0.0;
    F->m33 = dt;
    return 0.0;
}

/**
 * @brief Form state transformation jacobi matrix of
 *      1st order markov/random walk process
 * @param[out]   F      Jacobi matrix
 * @param[in]   dt      Time interval[s]
 * @param[in]   T       correlation time[s]
 * @warning dt should far less than T(dt << T)
 * @note if T == V0, this function will return random walk process factor(I3)
 *  markov process should not used if you can not comprehend what really it is.
 */
extern int jacobi_trans_markov(m3_t* F, double dt, const v3_t* T) {
    /*TODO: check*/
    F->m11 = T->x == 0.0 ? 1.0 : 1.0 - dt / T->x;
    F->m22 = T->y == 0.0 ? 1.0 : 1.0 - dt / T->y;
    F->m33 = T->z == 0.0 ? 1.0 : 1.0 - dt / T->z;
    F->m12 = 0.0;
    F->m13 = 0.0;
    F->m21 = 0.0;
    F->m23 = 0.0;
    F->m31 = 0.0;
    F->m32 = 0.0;
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from reb_e to Ebe
 * @param[out]  H           Jacobi matrix
 * @param[in]   Cbe         Current average attitude
 * @param[in]   lever_arm   gnss to imu lever arm[m]
 * @return status(0: OK)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P600, 14.111
 */
extern int jacobi_lc_meas_p2a_ECEF(m3_t* H, const m3_t* Cbe, const v3_t* lever_arm) {
    v3_t v = m3_mul_v3(*Cbe, *lever_arm);
    asymmetric_mat(&v, H);
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from veb_e to Ebe
 * @param[out]  H			Jacobi matrix
 * @param[in]   Cbe			Current average attitude
 * @param[in]   wib_b		imu output angluar rate[rad/s]
 * @param[in]   lever_arm	gnss to imu lever arm[m]
 * @return status(0: OK)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P543, 14.111
 */
extern int jacobi_lc_meas_v2a_ECEF(m3_t* H, const m3_t* Cbe, const v3_t* wib_b, const v3_t* lever_arm) {
    v3_t v1 = m3_mul_v3(*Cbe, v3_cross(*wib_b, *lever_arm));
    v3_t wie_e = { 0.0, 0.0, wgs84.wie };
    v3_t v2 = v3_cross(wie_e, m3_mul_v3(*Cbe, *lever_arm));
    v3_t v = v3_del(v1, v2);
    asymmetric_mat(&v, H);
    return 0;
}

/**
 * @brief Form measurement jacobi matrix from veb_e to bg
 * @param[out] 		H			Jacobi matrix
 * @param[in] 		Cbe			Current average attitude
 * @param[in] 		lever_arm	gnss lever arm[m]
 * @return status(0: 0K)
 * @note
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P600, 14.111
 */
extern int jacobi_lc_meas_v2bg_ECEF(m3_t* H, const m3_t* Cbe, const v3_t* lever_arm) {
    m3_t lx;
    asymmetric_mat(lever_arm, &lx);
    *H = m3_mul(*Cbe, lx);
    return 0;
}

extern int jacobi_tc_meas2a_ECEF(v3_t* dpda, const m3_t* Cbe, const v3_t* lever_arm, double* ei, double* ej) {
    m3_t T;
    double sl[3] = { 0 };

    *dpda = m3_mul_v3(*Cbe, *lever_arm);
    T = v3_askew(*dpda);

    if (ei) {
        for (int i = 0; i < 3; i++) sl[i] = -ei[i] + ej[i];
    }
    else {
        for (int i = 0; i < 3; i++) sl[i] = ej[i];
    }

    matmul("NN", 1, 3, 3, 1.0, sl, T.v, 0.0, dpda->v);
    return 1;
}

EXPORT int jacobi_tc_meas2p_LLH(v3_t* dpda, const v3_t* pos, double* ei, double* ej)
{
    double sl[3] = { 0 };
    m3_t C;
    dpos2decef(pos->v, C.v);

    if (ei) {
        for (int i = 0; i < 3; i++) sl[i] = -ei[i] + ej[i];
    }
    else {
        for (int i = 0; i < 3; i++) sl[i] = ej[i];
    }

    matmul("NN", 1, 3, 3, 1.0, sl, C.v, 0.0, dpda->v);
    return 1;
}

EXPORT int jacobi_tc_meas2a_LLH(v3_t* dpda, const m3_t* Cbn, const v3_t* pos, const v3_t* lever_arm, double* ei, double* ej)
{
    m3_t C = { 0 };
    dpos2decef(pos->v, C.v);

    m3_t T;
    double sl[3] = { 0 };

    *dpda = m3_mul_v3(*Cbn, *lever_arm);
    T = v3_askew(*dpda);
    T = m3_mul(C, T);

    if (ei) {
        for (int i = 0; i < 3; i++) sl[i] = -ei[i] + ej[i];
    }
    else {
        for (int i = 0; i < 3; i++) sl[i] = ej[i];
    }

    matmul("NN", 1, 3, 3, 1.0, sl, T.v, 0.0, dpda->v);

    return 1;
}


/**
 * @brief Form state transformation jacobi matrix from
 * @param[out]  F37   Jacobi matrix
 * @param[in]   Cbe Current average attitude.
 * @param[in]   dtheta  Angle increment[rad]
 * @return 0: OK
 * @note Ref:
 */
extern int jacobi_trans_a2sg_ECEF(m3_t* F, const m3_t* Cbe, const v3_t* dtheta) {
    m3_t T = m3_scalar(-1, m3_mul(*Cbe, v3_diag(*dtheta)));
    F = &T;
    return 0;
}

/**
 * @brief Form state transformation jacobi matrix
 * @param[out]  F26   Jacobi matrix
 * @param[in]   Cbe Current average attitude.
 * @param[in]   dv  Velocity increment[m/s]
 * @return 0: OK
 * @note Ref:
 */
extern int jacobi_trans_v2sa_ECEF(m3_t* F, const m3_t* Cbe, const v3_t* dv) {
    m3_t T = m3_scalar(-1, m3_mul(*Cbe, v3_diag(*dv)));
    F = &T;
    return 0;
}


/** LLH Jacobi matrix form-------------------------------------------- */
/** @note Ref: 陈起金. GNSS/INS松组合算法设计 */

extern void jacobi_trans_p2p_NED(m3_t* F, solins_t* sol, double dt) {
    m3_t T = O33;
    v3_t vel = sol->vel;
    earth_t eth = sol->eth;

    T.m11 = -vel.z / eth.RMh;
    T.m13 = vel.x / eth.RMh;
    T.m21 = vel.y * eth.tl / eth.RNh;
    T.m22 = -(vel.z + vel.x * eth.tl) / eth.RNh;
    T.m23 = vel.y / eth.RNh;

    T = m3_scalar(dt, T);
    T = m3_del(I33, T);
    F = &T;
}

extern void jacobi_trans_p2v_NED(m3_t* F, double dt) {
    m3_t T = m3_scalar(dt, I33);
    F = &T;
}

extern void jacobi_trans_v2p_NED(m3_t* F, solins_t* sol, double dt) {
    m3_t T = O33;
    v3_t vel = sol->vel;
    earth_t eth = sol->eth;
    double ve2, vn2, secB2, RMh2, RNh2;

    vn2 = vel.x * vel.x;
    ve2 = vel.y * vel.y;
    secB2 = (eth.cl / eth.sl) * (eth.cl / eth.sl);
    RMh2 = eth.RMh * eth.RMh;
    RNh2 = eth.RNh * eth.RNh;

    T.m11 = -2 * vel.y * wgs84.e * eth.cl / eth.RMh - ve2 * secB2 / eth.RMh / eth.RNh;
    T.m13 = vel.x * vel.z / RMh2 - ve2 * eth.tl / RNh2;
    T.m21 = 2 * wgs84.e * (vel.x * eth.cl - vel.z * eth.sl) / eth.RMh +
        vel.x * vel.y * secB2 / eth.RMh / eth.RNh;
    T.m23 = vel.y * vel.z + vel.x * vel.y * eth.tl / RNh2;
    T.m31 = 2 * wgs84.e * vel.y * eth.sl / eth.RMh;
    T.m33 = -ve2 / RNh2 - vn2 / RMh2 + 2 * eth.g / (sqrt(eth.RM * eth.RN) + sol->pos.z);

    T = m3_scalar(dt, T);
    F = &T;
}

extern void jacobi_trans_v2v_NED(m3_t* F, solins_t* sol, double dt) {
    m3_t T = O33;
    v3_t vel = sol->vel;
    earth_t eth = sol->eth;

    T.m11 = vel.z / eth.RMh;
    T.m12 = -2.0 * (wgs84.e * eth.sl + vel.y * eth.tl / eth.RNh);
    T.m13 = vel.x / eth.RMh;
    T.m21 = 2.0 * wgs84.e * eth.sl + vel.y * eth.tl / eth.RNh;
    T.m22 = (vel.z + vel.x * eth.tl) / eth.RNh;
    T.m23 = 2.0 * wgs84.e * eth.cl + vel.y / eth.RNh;
    T.m31 = -2.0 * vel.x / eth.RMh;
    T.m32 = -2.0 * (wgs84.e * eth.cl + vel.y / eth.RNh);

    T = m3_scalar(dt, T);
    T = m3_del(I33, T);
    F = &T;
}

extern void jacobi_trans_v2a_NED(m3_t* F, const m3_t* Cbn, const v3_t* dv) {
    /* [- Cbn * fib_b ]x * dt */
    v3_t Tib_e = v3_scalar(-1.0, m3_mul_v3(*Cbn, *dv));
    asymmetric_mat(&Tib_e, F);
}

extern void jacobi_trans_a2p_NED(m3_t* F, solins_t* sol, double dt) {
    m3_t T = O33;
    v3_t vel = sol->vel;
    earth_t eth = sol->eth;
    double RMh2, RNh2, secB2;

    RMh2 = eth.RMh * eth.RMh;
    RNh2 = eth.RNh * eth.RNh;
    secB2 = (eth.cl / eth.sl) * (eth.cl / eth.sl);

    T.m11 = -wgs84.e * eth.sl / eth.RMh;
    T.m13 = vel.y / RNh2;
    T.m23 = -vel.x / RMh2;
    T.m31 = -wgs84.e * eth.cl / eth.RMh - vel.y * secB2 / eth.RMh / eth.RNh;
    T.m33 = -vel.y * eth.tl / RNh2;

    T = m3_scalar(dt, T);
    F = &T;
}

extern void jacobi_trans_a2v_NED(m3_t* F, solins_t* sol, double dt) {
    m3_t T = O33;
    earth_t eth = sol->eth;

    T.m12 = 1.0 / eth.RNh;
    T.m21 = -1.0 / eth.RMh;
    T.m32 = -eth.tl / eth.RNh;

    T = m3_scalar(dt, T);
    F = &T;
}

extern void jacobi_trans_v2ba_NED(m3_t* F, const m3_t* Cbn, double dt) {
    *F = m3_scalar(dt, *Cbn);
}

extern void jacobi_trans_a2bg_NED(m3_t* F, const m3_t* Cbn, double dt) {
    *F = m3_scalar(-dt, *Cbn);
}

extern void jacobi_trans_v2sa_NED(m3_t* F, const m3_t* Cbn, const v3_t* dv) {
    *F = m3_mul(*Cbn, v3_diag(*dv));
}

extern void jacobi_trans_a2sg_NED(m3_t* F, const m3_t* Cbn, const v3_t* dtheta) {
    *F = m3_scalar(-1.0, m3_mul(*Cbn, v3_diag(*dtheta)));
}
