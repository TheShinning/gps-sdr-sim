/*------------------------------------------------------------------------------
* PVA-transform.cc :
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/4/19.
*-----------------------------------------------------------------------------*/

#include "rtklib.h"


/**
 * @brief Convert attitude angle(Enb) to DCM transform from b-frame to
 *      e-frame(Cbe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return DCM transform from b-frame to e-frame, Cbe
 */
m3_t att2Cbe(const v3_t* pos, const v3_t* Enb, int imu_opt, int local_opt) {
    m3_t Cbn;
    att2dcm(Enb, &Cbn, imu_opt);
    m3_t Cen = formCen_llh(pos->x, pos->y, local_opt);
    return m3_mul(m3_T(Cen), Cbn);
}

/**
 * @brief Convert attitude angle(Enb) to quaternion transform from b-frame to
 *      e-frame(Qbe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return quaternion transform from b-frame to e-frame, Qbe
 */
quat_t att2Qbe(const v3_t* pos, const v3_t* Enb, int imu_opt, int local_opt) {
    m3_t Cbe = att2Cbe(pos, Enb, imu_opt, local_opt);
    quat_t Qbe;
    dcm2quat(&Cbe, &Qbe);
    return Qbe;
}

/**
 * @brief Convert attitude angle(Enb) to Euler transform from b-frame to
 *      e-frame(Ebe)
 * @param[in] pos   Current geodetic position, BLH[rad,m]
 * @param[in] Enb   Current attitue angle(roll, pitch, yaw)[rad]
 * @return Ebe, Euler transform from b-frame to e-frame
 */
v3_t att2Ebe(const v3_t* pos, const v3_t* Enb, int imu_coord_opt, int local_opt) {
    m3_t Cbe = att2Cbe(pos, Enb, imu_coord_opt, local_opt);
    v3_t Ebe;
    dcm2euler(&Cbe, &Ebe, imu_coord_opt);
    return Ebe;
}

/**
 * @brief Convert quaternion transform from b-frame to e-frame(Qbe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Cbe   Current DCM transform from b-frame to e-frame, Cbe
 * @return attitude angle(roll, pitch, yaw)[Enb]
 */
v3_t Cbe2att(const v3_t* xyz, const m3_t* Cbe, int imu_opt, int local_opt) {
    v3_t pos = *xyz;
    m3_t Cbn = *Cbe;
    ecef2llh(&pos, nullptr, &Cbn, local_opt);
    v3_t att;
    dcm2att(&Cbn, &att, imu_opt);
    return att;
}

/**
 * @brief Convert quaternion transform from b-frame to e-frame(Qbe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Qbe   Current quaternion transform from b-frame to e-frame, Qbe
 * @return Enb, attitude angle(roll, pitch, yaw) [rad]
 */
v3_t Qbe2att(const v3_t* xyz, const quat_t* Qbe, int imu_opt, int local_opt) {
    m3_t Cbe;
    quat2dcm(Qbe, &Cbe, imu_opt);
    return Cbe2att(xyz, &Cbe, imu_opt, local_opt);
}

/**
 * @brief Convert Euler transform from b-frame to e-frame(Ebe) to attitude
 *      angle(Enb)
 * @param[in] xyz   Current geodetic position, BLH[rad,m]
 * @param[in] Ebe   Current Euler transform from b-frame to e-frame, Ebe[rad]
 * @return Enb, attitude angle(roll, pitch, yaw)[rad]
 */
v3_t Ebe2att(const v3_t* xyz, const v3_t* Ebe, int imu_opt, int local_opt) {
    m3_t Cbe;
    euler2dcm(Ebe, &Cbe, imu_opt);
    return Cbe2att(xyz, &Cbe, imu_opt, local_opt);
}

extern void pos2vel(const v3_t* pos1, const v3_t* pos2, v3_t* v, double dt) {
    *v = v3_scalar(1.0 / dt, v3_del(*pos1, *pos2));
}

static void vel2att(const v3_t vn, imup_t* imup, int local_opt) {
    imup->inita.x = imup->inita.y = 0.0;
    if (local_opt == INSLOCAL_NED) {
        imup->inita.z = atan2(vn.y, fabs(vn.x) < 1E-4 ? 1E-4 : vn.x);
    }
    else if (local_opt == INSLOCAL_ENU) {
        imup->inita.z = atan2(vn.x, fabs(vn.y) < 1E-4 ? 1E-4 : vn.y);
    }
}

extern void ins2gnss(const solins_t* ins_sol, const v3_t* arm, const v3_t* wib_b, double* gr, double* gv, int mech_coord) {
    m3_t dcm = ins_sol->dcm;

    if (mech_coord == INSMECH_ECEF) {
        if (gr != nullptr) {
            v3_t re = v3_del(ins_sol->pos, m3_mul_v3(dcm, *arm));
            gr[0] = re.x;
            gr[1] = re.y;
            gr[2] = re.z;
        }

        if (gv != nullptr && wib_b != nullptr) {
            v3_t ve_t = m3_mul_v3(dcm, v3_cross(*wib_b, *arm));
            v3_t ve = v3_add(ins_sol->vel, ve_t);
            gv[0] = ve.x;
            gv[1] = ve.y;
            gv[2] = ve.z;
        }
    }
    else if (mech_coord == INSMECH_LLH) {
        if (gr != nullptr) {
            v3_t re = v3_add(ins_sol->pos, m3_mul_v3(ins_sol->MpvCnb, *arm));
            gr[0] = re.x;
            gr[1] = re.y;
            gr[2] = re.z;
        }
        if (gv != nullptr) {
            v3_t ve = v3_add(ins_sol->vel, m3_mul_v3(ins_sol->CW, *arm));
            gv[0] = ve.x;
            gv[1] = ve.y;
            gv[2] = ve.z;
        }
    }
}

extern void gnss2ins(const v3_t* gr, const v3_t* gv, const v3_t* arm, const v3_t* wib_b, solins_t* ins_sol) {
    m3_t Cbe = ins_sol->dcm;

    if (gr != nullptr) {
        ins_sol->pos = v3_add(*gr, m3_mul_v3(Cbe, *arm));
    }

    if (gv != nullptr && wib_b != nullptr) {
        v3_t w = { 0 };
        for (int i = 0; i < 3; i++) w.v[i] = wib_b->v[i];
        ins_sol->vel = m3_mul_v3(Cbe, v3_cross(w, *arm));
    }
}

extern int gnsspva2ins(const v3_t re, const v3_t ve, const imud_t* imu_data, imup_t* imup, int imu_coord_opt, int local_opt) {
    v3_t wiee = { 0, 0, wgs84.wie }, rn = re, vn = ve;
    m3_t Cbe = O33;

    ecef2llh(&rn, &vn, nullptr, local_opt);
    vel2att(vn, imup, local_opt);

    /*TODO: check*/
#if 1
    Cbe = att2Cbe(&rn, &imup->inita, imu_coord_opt, local_opt); /*�˴�����̬��Ϊnϵ*/
    imup->initr = v3_add(re, m3_mul_v3(Cbe, imup->lever_arm_gps));
    //    for(int i=0;i<3;i++) imup->initr.v[i]+=3.0;
    m3_t T = v3_askew(v3_scalar(1.0 / imup->freq_imu, imu_data->gyro));
    m3_t omge = v3_askew(wiee);
    imup->initv = v3_del(v3_del(ve, m3_mul_v3(m3_mul(omge, T), imup->lever_arm_gps)),
        m3_mul_v3(m3_mul(omge, Cbe), imup->lever_arm_gps));
#else
    imup->initr = re;
    imup->initv = ve;
#endif
    return 1;
}

extern int vel2yaw(const v3_t* veb_n, double* yaw,
    const m3_t* Qveb_n, double* Qyaw) {
    double SQR_vel = SQR(veb_n->x) + SQR(veb_n->y);
    if (sqrt(SQR_vel) < 1.0) {
        return 1;
    }
    *yaw = atan2(veb_n->y, veb_n->x);
    if (Qveb_n != NULL && Qyaw != NULL) {
        v3_t A = { -veb_n->x / SQR_vel, veb_n->y / SQR_vel, 0.0 };
        *Qyaw = v3_mul_rxc(A, m3_mul_v3(*Qveb_n, A));
    }
    return 0;
}

//EXPORT void attsync(solins_t *sol_ins,const insopt_t *opt)
//{
//    if(opt->mech_coord==INSMECH_ECEF){  // convert ecef to llh
//        ecef2llhQ(&sol_ins->pos, &sol_ins->Qpos, &sol_ins->Qvel, &sol_ins->Qatt,opt->local_coord);
//        ecef2llh(&sol_ins->pos, &sol_ins->vel, &sol_ins->dcm,opt->local_coord);
//    }
//    dcm2att(&sol_ins->dcm,&sol_ins->rpy);
//}


/**
 * @brief Convert n-frame(NED) position/velocity/attitude to e-frame(ECEF)
 * @param[in,out] pos   Input (lat,lon,hgt)[rad,m] / Output ECEF position xyz[m]
 * @param[in,out] vel   Input velocity(vN, vE, vD)[m/s] / Output velocity ECEF[m/s]
 * @param[in,out] dcm   Input Cbn attitude / Ouput Cbe attitude
 * @return 0: OK
 * @see ecef2ned()
 * @note  vel and dcm could set to NULL pointer when do NOT interested.
 */
extern int llh2ecef(v3_t* pos, v3_t* vel, m3_t* dcm, int opt) {
    double lat = pos->x, lon = pos->y, hgt = pos->z;
    double coslat = cos(lat), sinlat = sin(lat);
    double coslon = cos(lon), sinlon = sin(lon);

    double tmp = wgs84.e * sinlat;
    double Re = wgs84.R0 / sqrt(1 - SQR(tmp));

    pos->x = (Re + hgt) * coslat * coslon;
    pos->y = (Re + hgt) * coslat * sinlon;
    pos->z = ((1 - wgs84.e * wgs84.e) * Re + hgt) * sinlat;

    if (vel != NULL || dcm != NULL) {
        m3_t Cne;
        if (opt == INSLOCAL_ENU) { /*ENU*/
            Cne.m11 = -sinlon;
            Cne.m12 = coslon;
            Cne.m13 = 0.0;
            Cne.m21 = -sinlat * coslon;
            Cne.m22 = -sinlat * sinlon;
            Cne.m23 = coslat;
            Cne.m31 = coslat * coslon;
            Cne.m32 = coslat * sinlon;
            Cne.m33 = sinlat;
            Cne = m3_T(Cne);
        }
        else if (opt == INSLOCAL_NED) {    /*NED*/
            Cne.m11 = -sinlat * coslon;
            Cne.m21 = -sinlat * sinlon;
            Cne.m31 = coslat;
            Cne.m12 = -sinlon;
            Cne.m22 = coslon;
            Cne.m32 = 0;
            Cne.m13 = -coslat * coslon;
            Cne.m23 = -coslat * sinlon;
            Cne.m33 = -sinlat;
        }
        if (vel != NULL)
            *vel = m3_mul_v3(Cne, *vel);    /* Veb_n => Veb_e */
        if (dcm != NULL)
            *dcm = m3_mul(Cne, *dcm);       /* Cb_n => Cb_e */
    }
    return 0;
}

/**
 * @brief Convert e-frame(ECEF) position/velocity/attitude to n-frame(NED)
 * @param[in,out] pos   Input ECEF position xyz[m] / Output (lat,lon,hgt)[rad,m]
 * @param[in,out] vel   Input velocity ECEF[m/s] / Output velocity(vN, vE, vD)[m/s]
 * @param[in,out] dcm   Input Cbe attitude / Output Cbn attitude
 * @see ned2ecef()
 * @return 0: OK
 * @note vel and dcm could set to NULL pointer when do NOT interested.
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, C.29-C.38
 */
extern int ecef2llh(v3_t* pos, v3_t* vel, m3_t* dcm, int opt) {
    double lon = atan2(pos->y, pos->x);

    double e2 = wgs84.e * wgs84.e;
    double k1 = sqrt(1.0 - e2) * fabs(pos->z);
    double k2 = e2 * wgs84.R0;
    double beta = sqrt(pos->x * pos->x + pos->y * pos->y);
    double E = (k1 - k2) / beta, F = (k1 + k2) / beta;
    double P = 4.0 / 3.0 * (E * F + 1.0);
    double Q = 2.0 * (E * E - F * F);
    double D = P * P * P + Q * Q;
    double V = pow(sqrt(D) - Q, 1.0 / 3.0) - pow(sqrt(D) + Q, 1.0 / 3.0);
    double G = 0.5 * (sqrt(E * E + V) + E);
    double T = sqrt(G * G + (F - V * G) / (2.0 * G - E)) - G;
    double signz = pos->z > 0 ? 1.0 : -1.0;
    double lat = signz * atan((1 - T * T) / (2 * T * sqrt(1 - e2)));

    double coslat = cos(lat), sinlat = sin(lat);
    double hgt = (beta - wgs84.R0 * T) * coslat
        + (pos->z - signz * wgs84.R0 * sqrt(1.0 - e2)) * sinlat;

    pos->x = lat;
    pos->y = lon;
    pos->z = hgt;

    if (vel != NULL || dcm != NULL) {
        double coslon = cos(lon), sinlon = sin(lon);
        m3_t C;
        if (opt == INSLOCAL_ENU) {
            C.m11 = -sinlon;
            C.m12 = -sinlat * coslon;
            C.m13 = coslat * coslon;
            C.m21 = coslon;
            C.m22 = -sinlat * sinlon;
            C.m23 = coslat * sinlon;
            C.m31 = 0.0;
            C.m32 = coslat;
            C.m33 = sinlat;
            C = m3_T(C);
        }
        else if (opt == INSLOCAL_NED) {
            C.m11 = -sinlat * coslon;
            C.m12 = -sinlat * sinlon;
            C.m13 = coslat;
            C.m21 = -sinlon;
            C.m22 = coslon;
            C.m23 = 0;
            C.m31 = -coslat * coslon;
            C.m32 = -coslat * sinlon;
            C.m33 = -sinlat;
        }
        if (vel != nullptr)
            *vel = m3_mul_v3(C, *vel); /* Veb_e => Veb_n */
        if (dcm != nullptr)
            *dcm = m3_mul(C, *dcm); /* Cbe => Cbn */
    }
    return 0;
}

/**
 * @brief form the DCM transform from e-frame to n-frame(Cbe)
 * @param[in] lat   Current postion's latitude[rad]
 * @param[in] lon   Current postion's longitude[rad]
 * @return Cen, DCM attitude transform from e-frame to n-frame
 */
extern m3_t formCen_llh(double lat, double lon, int opt) {
    double coslat = cos(lat), sinlat = sin(lat);
    double coslon = cos(lon), sinlon = sin(lon);
    m3_t Cen;
    if (opt == INSLOCAL_ENU) {
        Cen.m11 = -sinlon;
        Cen.m12 = -sinlat * coslon;
        Cen.m13 = coslat * coslon;
        Cen.m21 = coslon;
        Cen.m22 = -sinlat * sinlon;
        Cen.m23 = coslat * sinlon;
        Cen.m31 = 0.0;
        Cen.m32 = coslat;
        Cen.m33 = sinlat;
    }
    else if (opt == INSLOCAL_NED) {
        Cen.m11 = -sinlat * coslon;
        Cen.m12 = -sinlat * sinlon;
        Cen.m13 = coslat;
        Cen.m21 = -sinlon;
        Cen.m22 = coslon;
        Cen.m23 = 0;
        Cen.m31 = -coslat * coslon;
        Cen.m32 = -coslat * sinlon;
        Cen.m33 = -sinlat;
    }

    return Cen;
}


///* format trasform *///

extern void sol2pva(const insopt_t* opt, int pva_idx, sol_t* sol, pva_t* pva) {
    pva->time[pva_idx] = sol->time;
    pva->pos[pva_idx] = { sol->rr[0], sol->rr[1], sol->rr[2] };
    pva->vel[pva_idx] = { sol->rr[3], sol->rr[4], sol->rr[5] };
    pva->att[pva_idx] ={ 0, 0, 0 };

    pva->Qpos[pva_idx].m11 = ((double)sol->qr[0]);
    pva->Qpos[pva_idx].m22 = ((double)sol->qr[1]);
    pva->Qpos[pva_idx].m33 = ((double)sol->qr[2]);
    pva->Qpos[pva_idx].m21 = pva->Qpos[pva_idx].m12 = ((double)sol->qr[3]);
    pva->Qpos[pva_idx].m32 = pva->Qpos[pva_idx].m23 = ((double)sol->qr[4]);
    pva->Qpos[pva_idx].m31 = pva->Qpos[pva_idx].m13 = ((double)sol->qr[5]);

    pva->Qvel[pva_idx].m11 = ((double)sol->qv[0]);
    pva->Qvel[pva_idx].m22 = ((double)sol->qv[1]);
    pva->Qvel[pva_idx].m33 = ((double)sol->qv[2]);
    pva->Qvel[pva_idx].m21 = pva->Qvel[pva_idx].m12 = (double)sol->qv[3];
    pva->Qvel[pva_idx].m32 = pva->Qvel[pva_idx].m23 = (double)sol->qv[4];
    pva->Qvel[pva_idx].m31 = pva->Qvel[pva_idx].m13 = (double)sol->qv[5];

    if (opt->mech_coord == INSMECH_LLH && fabs(pva->pos[pva_idx].x) > 100 && fabs(pva->pos[pva_idx].y > 200)) {
        v3_t re = pva->pos[pva_idx];
        ecef2llh(&pva->pos[pva_idx], &pva->vel[pva_idx], nullptr, opt->local_coord);
        ecef2llhQ(&re, &pva->Qpos[pva_idx], &pva->Qvel[pva_idx], nullptr, opt->local_coord);
    }

    pva->status[pva_idx] = sol->stat;
}

extern void pva2sol(const pva_t* pva, sol_t* sol, int pva_idx) {
    sol->time = pva->time[pva_idx];
    sol->stat = pva->status[pva_idx];
    sol->rr[0] = pva->pos[pva_idx].x;
    sol->rr[1] = pva->pos[pva_idx].y;
    sol->rr[2] = pva->pos[pva_idx].z;
    sol->rr[3] = pva->vel[pva_idx].x;
    sol->rr[4] = pva->vel[pva_idx].y;
    sol->rr[5] = pva->vel[pva_idx].z;
    sol->qr[0] = (float)pva->Qpos[pva_idx].m11;
    sol->qr[1] = (float)pva->Qpos[pva_idx].m22;
    sol->qr[2] = (float)pva->Qpos[pva_idx].m33;
    sol->qr[3] = (float)pva->Qpos[pva_idx].m21;
    sol->qr[4] = (float)pva->Qpos[pva_idx].m32;
    sol->qr[5] = (float)pva->Qpos[pva_idx].m31;
    sol->qv[0] = (float)pva->Qvel[pva_idx].m11;
    sol->qv[1] = (float)pva->Qvel[pva_idx].m22;
    sol->qv[2] = (float)pva->Qvel[pva_idx].m33;
    sol->qv[3] = (float)pva->Qvel[pva_idx].m21;
    sol->qv[4] = (float)pva->Qvel[pva_idx].m32;
    sol->qv[5] = (float)pva->Qvel[pva_idx].m31;
}

extern void sol2solins(const sol_t* sol, solins_t* solins) {
    solins->time = sol->time;
    solins->status = sol->stat;
    solins->pos = { sol->rr[0], sol->rr[1], sol->rr[2] };
    solins->vel = { sol->rr[3], sol->rr[4], sol->rr[5] };
    solins->Qpos.m11 = (double)(sol->qr[0]);
    solins->Qpos.m22 = (double)(sol->qr[1]);
    solins->Qpos.m33 = (double)(sol->qr[2]);
    solins->Qpos.m21 = solins->Qpos.m12 = (double)(sol->qr[3]);
    solins->Qpos.m32 = solins->Qpos.m32 = (double)(sol->qr[4]);
    solins->Qpos.m31 = solins->Qpos.m13 = (double)(sol->qr[5]);

    solins->Qvel.m11 = (double)(sol->qv[0]);
    solins->Qvel.m22 = (double)(sol->qv[1]);
    solins->Qvel.m33 = (double)(sol->qv[2]);
    solins->Qvel.m21 = solins->Qvel.m12 = (double)(sol->qv[3]);
    solins->Qvel.m32 = solins->Qvel.m23 = (double)(sol->qv[4]);
    solins->Qvel.m31 = solins->Qvel.m13 = (double)(sol->qv[5]);
}

extern void solins2sol(const insopt_t* opt, const solins_t* solins, sol_t* sol, int att_fmt) {
    sol->time = solins->time;
    sol->stat = solins->status;
    solins_t solins_copy = *solins;

    if (opt->mech_coord == INSMECH_LLH) {
        llh2ecef(&solins_copy.pos, &solins_copy.vel, nullptr, opt->local_coord);
        llh2ecefQ(&solins->pos, &solins_copy.Qpos, &solins_copy.Qvel, &solins_copy.Qatt, opt->local_coord);
    }
    sol->rr[0] = solins_copy.pos.x;
    sol->rr[1] = solins_copy.pos.y;
    sol->rr[2] = solins_copy.pos.z;
    sol->rr[3] = solins_copy.vel.x;
    sol->rr[4] = solins_copy.vel.y;
    sol->rr[5] = solins_copy.vel.z;

    sol->qr[0] = (float)solins_copy.Qpos.m11;
    sol->qr[1] = (float)solins_copy.Qpos.m22;
    sol->qr[2] = (float)solins_copy.Qpos.m33;
    sol->qr[3] = (float)solins_copy.Qpos.m21;
    sol->qr[4] = (float)solins_copy.Qpos.m32;
    sol->qr[5] = (float)solins_copy.Qpos.m31;
    sol->qv[0] = (float)solins_copy.Qvel.m11;
    sol->qv[1] = (float)solins_copy.Qvel.m22;
    sol->qv[2] = (float)solins_copy.Qvel.m33;
    sol->qv[3] = (float)solins_copy.Qvel.m21;
    sol->qv[4] = (float)solins_copy.Qvel.m32;
    sol->qv[5] = (float)solins_copy.Qvel.m31;

    attsync(&solins_copy.pos, &solins_copy.rpy, &solins_copy.dcm, &solins_copy.quat, opt->imu_coord, opt->local_coord, opt->mech_coord, att_fmt);
    sol->rpy[0] = solins_copy.rpy.x;
    sol->rpy[1] = solins_copy.rpy.y;
    sol->rpy[2] = solins_copy.rpy.z;
    sol->qa[0] = (float)solins_copy.Qatt.m11;
    sol->qa[1] = (float)solins_copy.Qatt.m22;
    sol->qa[2] = (float)solins_copy.Qatt.m33;

    sol->ba[0] = solins_copy.ba.x;
    sol->ba[1] = solins_copy.ba.y;
    sol->ba[2] = solins_copy.ba.z;
    sol->bg[0] = solins_copy.bg.x;
    sol->bg[1] = solins_copy.bg.y;
    sol->bg[2] = solins_copy.bg.z;
    sol->qba[0] = (float)solins_copy.ba_std.x;
    sol->qba[1] = (float)solins_copy.ba_std.y;
    sol->qba[2] = (float)solins_copy.ba_std.z;
    sol->qbg[0] = (float)solins_copy.bg_std.x;
    sol->qbg[1] = (float)solins_copy.bg_std.y;
    sol->qbg[2] = (float)solins_copy.bg_std.z;

    sol->sa[0] = solins->sa.x;
    sol->sa[1] = solins->sa.y;
    sol->sa[2] = solins->sa.z;
    sol->sg[0] = solins->sg.x;
    sol->sg[1] = solins->sg.y;
    sol->sg[2] = solins->sg.z;
    sol->qsa[0] = (float)solins->sa_std.x;
    sol->qsa[1] = (float)solins->sa_std.y;
    sol->qsa[2] = (float)solins->sa_std.z;
    sol->qsg[0] = (float)solins->sg_std.x;
    sol->qsg[1] = (float)solins->sg_std.y;
    sol->qsg[2] = (float)solins->sg_std.z;

    sol->arm_gps[0] = solins->arm_gps.x;
    sol->arm_gps[1] = solins->arm_gps.y;
    sol->arm_gps[2] = solins->arm_gps.z;
    sol->qarm_gps[0] = solins->arm_gps_std.x;
    sol->qarm_gps[1] = solins->arm_gps_std.y;
    sol->qarm_gps[2] = solins->arm_gps_std.z;
}

/**
 * @brief copy varainace(stanadard error) from inskf->P to inskf->sol.
 * @param   inskf   kalman filter struct
 */
extern void kf_copyQ2sol(kf_t* inskf, const insopt_t* opt) {
#define GETSTD(x)   sqrt(inskf->P[(x)+(x)*inskf->nx])
    m3_copy(&inskf->insstate->Qpos, &inskf->P[xiP(opt) + xiP(opt) * inskf->nx], inskf->nx);
    m3_copy(&inskf->insstate->Qvel, &inskf->P[xiV(opt) + xiV(opt) * inskf->nx], inskf->nx);
    m3_copy(&inskf->insstate->Qatt, &inskf->P[xiA(opt) + xiA(opt) * inskf->nx], inskf->nx);
    if (opt->est_ba) {
        inskf->insstate->ba_std.x = GETSTD(xiBa(opt));
        inskf->insstate->ba_std.y = GETSTD(xiBa(opt) + 1);
        inskf->insstate->ba_std.z = GETSTD(xiBa(opt) + 2);
    }
    if (opt->est_bg) {
        inskf->insstate->bg_std.x = GETSTD(xiBg(opt));
        inskf->insstate->bg_std.y = GETSTD(xiBg(opt) + 1);
        inskf->insstate->bg_std.z = GETSTD(xiBg(opt) + 2);
    }
    if (opt->est_sa) {
        inskf->insstate->sa_std.x = GETSTD(xiSa(opt));
        inskf->insstate->sa_std.y = GETSTD(xiSa(opt) + 1);
        inskf->insstate->sa_std.z = GETSTD(xiSa(opt) + 2);
    }
    if (opt->est_sg) {
        inskf->insstate->sg_std.x = GETSTD(xiSg(opt));
        inskf->insstate->sg_std.y = GETSTD(xiSg(opt) + 1);
        inskf->insstate->sg_std.z = GETSTD(xiSg(opt) + 2);
    }

    if (opt->est_kod) inskf->insstate->std_kod = GETSTD(xiKod(opt));
#undef GETSTD
}
