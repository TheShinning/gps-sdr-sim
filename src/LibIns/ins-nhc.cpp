/*------------------------------------------------------------------------------
* ins-nhc.cc : Non Holonomic Constraint update for ins navigation
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/28.
*-----------------------------------------------------------------------------*/

#include "rtklib.h"

/* Simplfiy the filter usage */
#define KF_FILTER(kf, dz, Qdz)                                               \
    filter(kf->x, kf->P, kf->H, (dz), (Qdz), kf->nx, kf->ny, 0, 0, NULL, 0)

/* set kf->ny and set kf->H to zero */
#define KF_HINIT(kf, _ny)                                                    \
{                                                                            \
    kf->ny = _ny;                                                            \
    for(unsigned int i = 0; i < kf->ny*kf->nx; ++i)                          \
        kf->H[i] = 0.0;                                                      \
}

extern int inskf_NHC(kf_t* inskf, const v3_t* wib_b, const insopt_t* opt) {
	KF_HINIT(inskf, 2);
	v3_t veb_b = m3_mul_v3(m3_T(inskf->insstate->dcm), inskf->insstate->vel);
	v3_t vec_b = v3_del(veb_b, v3_cross(*wib_b, opt->imup.lever_arm_car));
	v3_t vec_c = m3_mul_v3(inskf->insstate->Cbc, vec_b);

	double dz[2];
	dz[0] = vec_c.y;
	dz[1] = vec_c.z;

	m3_t Cec = m3_mul(inskf->insstate->Cbc, m3_T(inskf->insstate->dcm));
	m3_t Mphi = m3_mul(Cec, v3_askew(inskf->insstate->vel));

	inskf->H[xiA(opt)] = Mphi.m21;
	inskf->H[xiA(opt) + 1] = Mphi.m22;
	inskf->H[xiA(opt) + 2] = Mphi.m23;
	inskf->H[xiA(opt) + inskf->nx] = Mphi.m31;
	inskf->H[xiA(opt) + 1 + inskf->nx] = Mphi.m32;
	inskf->H[xiA(opt) + 2 + inskf->nx] = Mphi.m33;

	inskf->H[xiV(opt)] = Cec.m21;
	inskf->H[xiV(opt) + 1] = Cec.m22;
	inskf->H[xiV(opt) + 2] = Cec.m23;
	inskf->H[xiV(opt) + inskf->nx] = Cec.m31;
	inskf->H[xiV(opt) + 1 + inskf->nx] = Cec.m32;
	inskf->H[xiV(opt) + 2 + inskf->nx] = Cec.m33;

	//    matprint(inskf->H, inskf->nx, ny, 12, 6);
	//    fflush(stdout);

	//    std_vy = std_vy * (fabs(wib_b->z) > 0.01 ? fabs(wib_b->z)/0.01 : 1.0 );
	//    LOG_INFO("std_vy: %f, ratio %f", std_vy, fabs(wib_b->z)/0.01);

	double Qdz[4] = { SQR(opt->nhc_var), 0.0, 0.0, SQR(opt->nhc_var) };
	double ndz[2];
	inskf_norm_innov(inskf, dz, Qdz, ndz);  // check norm
	if (fabs(ndz[0]) > 10.0 || fabs(ndz[1]) > 10.0) {
		return 1;
	}

	KF_FILTER(inskf, dz, Qdz);
	inskf_feedback(inskf->time, SOLQ_NHC_AID, opt, inskf->x, inskf->insstate);
	return 0;
}