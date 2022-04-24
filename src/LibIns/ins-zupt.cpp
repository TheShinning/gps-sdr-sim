/*------------------------------------------------------------------------------
* ins-zupt.cc : zero velocity and zero angular rate update for ins navigation
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

/* ins detect static imu measurement functions--------------------------------*/
/* static detector for per accl axis,this function check accl variance--------*/
static int detstaticaxis(const imud_t* imu, int n, int axis, double thres) {
	int i;
	double mean, var, * a = mat(n, 1);

	for (mean = 0.0, i = 0; i < n; i++) mean += imu[i].accel.v[axis];
	mean /= n;
	for (i = 0; i < n; i++) a[i] = imu[i].accel.v[axis] - mean;

	matmul("NT", 1, 1, n, 1.0 / (n - 1), a, a, 0.0, &var);

	/* detect static for imu measurement data */
	for (i = 0; i < n; i++) {
		if (fabs(a[i] / SQRT(var)) < thres) continue;
		free(a);
		return 0;
	}
	free(a);
	return 1;
}

/* static detector by gyro measurement --------------------------------------*/
static int detstaticgyro(const imud_t* imu, int n, const double* thres) {
	int i;
	for (i = 0; i < n; i++) {
		if (fabs(imu[i].gyro.v[0]) < thres[0] &&
			fabs(imu[i].gyro.v[1]) < thres[1] &&
			fabs(imu[i].gyro.v[2]) < thres[2])
			continue;
		return 0;
	}
	return 1;
}

/* static detector for gravity in ned-frame----------------------------------*/
static int detstaticg(const double* pos, const imud_t* imu, int n, double thres) {
	int i;
	v3_t gn;

	gravity_ned(pos[0], pos[2], &gn);
	for (i = 0; i < n; i++) {
		if (fabs(norm(imu[i].accel.v, 3) - norm(gn.v, 3)) < thres) continue;
		return 0;
	}
	return 1;
}

extern int detstatic_MEAS(const insopt_t* opt, const imud_t* imus, int ws, const double* pos) {
	return detstaticaxis(imus, ws, 0, opt->zvopt.athres[0]) &&
		detstaticaxis(imus, ws, 1, opt->zvopt.athres[1]) &&
		detstaticaxis(imus, ws, 2, opt->zvopt.athres[2]) &&
		detstaticgyro(imus, ws, opt->zvopt.gyrothres);
}

/* runs the generalized likelihood test for detect static imu measurement----
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 *           double *pos     I  ins position (lat,lon,h)
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_GLRT(const imud_t* imu, int n, const insopt_t* opt,
	const double* pos, double dt, double* z) {
	int i, j;
	double ym[3] = { 0 }, tmp[3];
	double T = 0, sg = opt->zvopt.sig_g, sa = opt->zvopt.sig_a;
	v3_t gn;

	trace(3, "detstatic_GLRT:\n");

	gravity_ned(pos[0], pos[2], &gn);

	for (i = 0; i < 3; i++) {
		for (j = 0; j < n; j++) ym[i] += imu[j].accel.v[i] / dt;
		ym[i] /= n;
	}
	double gyro[3];
	/* detect zero velocity */
	for (i = 0; i < n; i++) {
		v3_2_array(v3_scalar(1 / dt, imu[i].gyro), gyro);
		for (j = 0; j < 3; j++) tmp[j] = imu[i].accel.v[j] / dt - norm(gn.v, 3) / norm(ym, 3) * ym[j];
		T += SQR(norm(gyro, 3)) / SQR(sg) + SQR(norm(tmp, 3)) / SQR(sa);
	}
	T /= n;

	trace(3, "T=%6.4lf,gamma=%6.4lf\n", T, opt->zvopt.gamma[0]);
	//printf("T=%6.4lf,gamma=%6.4lf\n", T, opt->zvopt.gamma[0]);
	if (z) {
		*z = T;
	}

	return T < opt->zvopt.gamma[0];
	// *z = T;
}

/* runs the acceleration moving variance detector ---------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_MV(const imud_t* imu, int n, const insopt_t* opt, double dt, double* z) {
	int i, j;
	double ym[3], tmp[3], T = 0.0;

	trace(3, "detstatic_MV:\n");

	for (i = 0; i < 3; i++) {
		for (j = 0; j < n; j++) ym[i] += imu[j].accel.v[i] / dt;
		ym[i] /= n;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < 3; j++) tmp[j] = imu[i].accel.v[j] / dt - ym[j];
		T += SQR(norm(tmp, 3));
	}
	T /= (SQR(opt->zvopt.sig_a) * n);

	trace(3, "T=%6.4lf,gamma=%6.4lf\n", T, opt->zvopt.gamma[1]);
	if (z) {
		*z = T;
	}
	return T < opt->zvopt.gamma[1];
}

/* runs the acceleration magnitude detector ---------------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 *           double *pos     I  ins position (lat,lon,h)
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_MAG(const imud_t* imu, int n, const insopt_t* opt,
	const double* pos, double dt, double* z) {
	int i, j;
	double sa2 = SQR(opt->zvopt.sig_a), T;
	v3_t gn;

	trace(3, "detstatic_MAG:\n");

	gravity_ned(pos[0], pos[2], &gn);
	double acc[3];
	for (T = 0.0, i = 0; i < n; i++) {
		v3_2_array(v3_scalar(1 / dt, imu[i].accel), acc);
		T += SQR(norm(gn.v, 3) - norm(acc, 3));
	}
	T /= (sa2 * n);

	trace(3, "T=%6.4lf,gamma=%6.4lf\n", T, opt->zvopt.gamma[2]);
	if (z) {
		*z = T;
	}
	return T < opt->zvopt.gamma[2];
}

/* runs the angular rate energy detector -------------------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * --------------------------------------------------------------------------*/
extern int detstatic_ARE(const imud_t* imu, int n, const insopt_t* opt, double dt, double* z) {
	int i;
	double T = 0.0, sg2 = SQR(opt->zvopt.sig_g);

	trace(3, "detstatic_ARE:\n");
	double gyro[3];
	for (i = 0; i < n; i++) {
		v3_2_array(v3_scalar(1 / dt, imu[i].gyro), gyro);
		T += SQR(norm(gyro, 3));
	}
	T /= (sg2 * n);

	trace(3, "T=%6.4lf,gamma=%6.4lf\n", T, opt->zvopt.gamma[3]);
	if (z) {
		*z = T;
	}
	return T < opt->zvopt.gamma[3];
}

extern int detimustatic(const insopt_t* opt, const imud_t* imu, const v3_t solpos, double dt) {
	int det = opt->zvopt.det, ws = opt->zvopt.ws;
	double xyz[3], pos[3];

	v3_2_array(solpos, xyz);
	ecef2pos(xyz, pos);

	switch (det) {
	case IMUDETST_GLRT:
		return detstatic_GLRT(imu, ws, opt, pos, dt, NULL);
	case IMUDETST_MV:
		return detstatic_MV(imu, ws, opt, dt, NULL);
	case IMUDETST_MAG:
		return detstatic_MAG(imu, ws, opt, pos, dt, NULL);
	case IMUDETST_ARE:
		return detstatic_ARE(imu, ws, opt, dt, NULL);
	case IMUDETST_ALL:
		return detstatic_GLRT(imu, ws, opt, pos, dt, NULL) && detstatic_MV(imu, ws, opt, dt, NULL) &&
			detstatic_MAG(imu, ws, opt, pos, dt, NULL) && detstatic_ARE(imu, ws, opt, dt, NULL);
	default:
		return detstatic_GLRT(imu, ws, opt, pos, dt, NULL);
	}
}

/* zero velocity constraint, using in static status */
extern int inskf_ZUPT(kf_t* inskf, const insopt_t* opt) {
	double* Qvel;
	double ndz[3];
	KF_HINIT(inskf, 3);
	Qvel = zeros(3, 3);

	for (int i = 0; i < 3; i++) Qvel[i + i * 3] = opt->zupt_var;
	v3_t dz_v = inskf->insstate->vel;

	inskf->H[xiV(opt)] = 1.0;
	inskf->H[xiV(opt) + 1 + 1 * inskf->nx] = 1.0;
	inskf->H[xiV(opt) + 2 + 2 * inskf->nx] = 1.0;

	inskf_norm_innov(inskf, (const double*)&dz_v, Qvel, ndz);  // check norm
	if (fabs(ndz[0]) > 10.0 || fabs(ndz[1]) > 10.0 || fabs(ndz[2]) > 10.0) {
		return 1;
	}

	KF_FILTER(inskf, (const double*)&dz_v, Qvel);

	inskf_feedback(inskf->time, SOLQ_ZV_AID, opt, inskf->x, inskf->insstate);
	return 0;
}

/* zero angle constraint, using in static status or low cost mems no change in driving direction */
extern int inskf_ZARU(kf_t* inskf, const insopt_t* opt) {
	double* Qvel;
	double ndz[3];
	KF_HINIT(inskf, 3);
	Qvel = zeros(3, 3);

	for (int i = 0; i < 3; i++) Qvel[i + i * 3] = opt->zaru_var;
	v3_t dz_v = inskf->dthetap;

	inskf->H[xiBg(opt)] = 1.0;
	inskf->H[xiBg(opt) + 1 + 1 * inskf->nx] = 1.0;
	inskf->H[xiBg(opt) + 2 + 2 * inskf->nx] = 1.0;

	inskf_norm_innov(inskf, (const double*)&dz_v, Qvel, ndz);  // check norm
	if (fabs(ndz[0]) > 10.0 || fabs(ndz[1]) > 10.0 || fabs(ndz[2]) > 10.0) {
		return 1;
	}

	KF_FILTER(inskf, (const double*)&dz_v, Qvel);
	inskf_feedback(inskf->time, SOLQ_ZV_AID, opt, inskf->x, inskf->insstate);
	return 0;
}

/* zero yaw angle constraint, using in low cost mems no change in driving direction */
// todo check
extern int inskf_ZAU(kf_t* inskf, double yaw, double Qyaw, const insopt_t* opt) {
	KF_HINIT(inskf, 1);
	v3_t att = Cbe2att(&inskf->insstate->pos, &inskf->insstate->dcm, opt->imu_coord, opt->local_coord);
	double dz_yaw = yaw_del(yaw, att.z);

	v3_t pos = inskf->insstate->pos;
	ecef2llh(&pos, NULL, NULL, opt->local_coord);
	m3_t Cen = formCen_llh(pos.x, pos.y, opt->local_coord);

	inskf->H[xiA(opt)] = -Cen.m31;
	inskf->H[xiA(opt) + 1] = -Cen.m32;
	inskf->H[xiA(opt) + 2] = -Cen.m33;

	double ndz_yaw;
	inskf_norm_innov(inskf, &dz_yaw, &Qyaw, &ndz_yaw);

	if (fabs(ndz_yaw) > 5.0) {
		return 1;
	}

	KF_FILTER(inskf, &dz_yaw, &Qyaw);
	return 0;
}