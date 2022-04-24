/*-----------------------------------------------------------------------------
* ins-tmechanization.cc : nertial navigation system mechanization function
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/29.
*----------------------------------------------------------------------------*/
#include "rtklib.h"

/**
 * @brief Use multi-subsample to compensate the conning&scull error
 * @param[in]   dtheta_list     Angular increment list,order by time,length:abs(N)
 * @param[in]   dv_list         Velocity increment list,order by time,length:abs(N)
 * @param[in]   N               Subsample number(1=<N<=5, -2: one-plus-previous)
 * @param[out]  dtheta          Sum of angular increment with conning error compensation
 * @param[out]  dv              Sum of velocity incrment with scull error compensation
 * @return 0: OK    1: Error
 * @note Ref: \n
 *  1. Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition) \n
 *  2. Yan Gongming, 捷联惯导算法与组合导航讲义, 2016
 */
extern int multisample(const v3_t* dtheta_list, const v3_t* dv_list, int N, v3_t* dtheta, v3_t* dv) {
	if (abs(N) == 1) {
		*dtheta = dtheta_list[0];
		*dv = dv_list[0];
		return 0;
	}
	/* Multi-subsample coning compensation */
	if (N >= 2 && N <= 5) {
		// coning error confficients, ref Yan2016(P38:table 2.5-2)
		static double conefactors[5][4] = {
				{2. / 3},                                           // 2
				{9. / 20,    27. / 20},                             // 3
				{54. / 105,  92. / 105,  214. / 105},               // 4
				{250. / 504, 525. / 504, 650. / 504, 1375. / 504}   // 5
		};
		double* pcf = conefactors[N - 2];
		v3_t sum_c = { 0.0, 0.0, 0.0 }, sum_s = { 0.0, 0.0, 0.0 };
		v3_t sum_w = { 0.0, 0.0, 0.0 }, sum_v = { 0.0, 0.0, 0.0 };
		int i = 0;
		for (; i < N - 1; i++) {
			/* sum of cross_product factor */
			sum_c = v3_add(sum_c, v3_scalar(pcf[i], dtheta_list[i]));
			sum_s = v3_add(sum_s, v3_scalar(pcf[i], dv_list[i]));
			/* sum of angular increment and velocity increment */
			sum_w = v3_add(sum_w, dtheta_list[i]);
			sum_v = v3_add(sum_v, dv_list[i]);
		}
		sum_w = v3_add(sum_w, dtheta_list[i]);
		sum_v = v3_add(sum_v, dv_list[i]);
		/* coning error compensation for angular increment, ref Paul2013(5.97) */
		v3_t dphim = v3_cross(sum_c, dtheta_list[i]);
		*dtheta = v3_add(sum_w, v3_cross(sum_c, dtheta_list[i]));
		/* sculling error compensation for velocity increment, ref Paul2013(5.98)*/
		v3_t scul = v3_add(
			v3_cross(sum_c, dv_list[i]), v3_cross(sum_s, dtheta_list[i]));
		v3_t rotm = v3_scalar(0.5, v3_cross(sum_w, sum_v));
		*dv = v3_add(v3_add(sum_v, scul), rotm);
		return 0;
	}
	if (N == -2) {
		v3_t sum_c = v3_scalar(1.0 / 12, dtheta_list[0]);
		v3_t sum_s = v3_scalar(1.0 / 12, dv_list[0]);
		/* ref Yan2016(P31:2.5-37) */
		*dtheta = v3_add(dtheta_list[1], v3_cross(sum_c, dtheta_list[1]));
		/* ref Yan2016(P73:4.1-36,P76:4.1-55) */
		v3_t scul = v3_add(v3_cross(sum_c, dv_list[1]),
			v3_cross(sum_s, dtheta_list[1]));
		*dv = v3_add(dv_list[1], scul);
		return 0;
	}
	return 1;
}

extern int ins_nav_llh(double dt, const v3_t* dtheta, const v3_t* dv, solins_t* sol, int opt) {
	if (opt == INSLOCAL_NED) {
		/* Local variable definition ------------------------------------------- */
		m3_t theta;
		v3_t zeta, vng_cor, vbf, vnf;
		v3_t old_v = sol->vel;
		quat_t q_nk_nk1, q_bk1_bk;

		/* Code start ---------------------------------------------------------- */
		updearth(&sol->eth, sol->pos, sol->vel, opt, 0);

		/* Attitude update */
		zeta = v3_scalar(-dt, sol->eth.wnin);
		rv2quat(&zeta, &q_nk_nk1);    //update n
		rv2quat(dtheta, &q_bk1_bk);   //update b

		sol->quat = quat_mul(quat_mul(q_nk_nk1, sol->quat), q_bk1_bk);
		quat_normalize(&sol->quat);

		/* Velocity update */
		vng_cor = v3_scalar(dt, v3_del(sol->eth.gn, v3_cross(sol->eth.wnien, sol->vel)));  // gravity and Coriolis
		vbf = v3_add(*dv, v3_scalar(0.5, v3_cross(*dtheta, *dv)));
		theta = v3_askew(zeta);
		vnf = m3_mul_v3(m3_mul(m3_del(I33, m3_scalar(0.5, theta)), sol->dcm), vbf); // specific force
		sol->vel = v3_add(v3_add(sol->vel, vng_cor), vnf);

		quat2dcm(&sol->quat, &sol->dcm, IMUCOOR_FRD);

		/* Position update */
		sol->pos.x = sol->pos.x + 0.5 * (sol->vel.x + old_v.x) / (sol->eth.RMh) * dt;
		sol->pos.y = sol->pos.y + 0.5 * (sol->vel.y + old_v.y) / (sol->eth.RNh * sol->eth.cl) * dt;
		sol->pos.z = sol->pos.z - 0.5 * (sol->vel.z + old_v.z) * dt;

		return 1;
	}
	else if (opt == INSLOCAL_ENU) {
		v3_t vn01, pos01;
		vn01 = v3_add(sol->vel, v3_scalar(dt / 2.0, sol->an));
		pos01 = v3_add(sol->pos, m3_mul_v3(sol->eth.Mpv, v3_scalar(dt / 2.0, vn01)));
		updearth(&sol->eth, pos01, vn01, opt, 0);

		sol->wib = v3_scalar(1 / dt, *dtheta);
		sol->fb = v3_scalar(1 / dt, *dv);
		v3_t wie = m3_mul_v3(sol->dcm, sol->eth.wnie);
		sol->web = v3_del(sol->wib, wie);
		v3_t dtheta2 = v3_scalar(0.5, *dtheta);
		m3_t m1, m2, m3;
		rv2dcm(&dtheta2, &m1);
		m2 = m3_mul(sol->dcm, m1);
		m3 = m3_T(m2);
		sol->wnb = v3_del(sol->wib, m3_mul_v3(m3, sol->eth.wnin));

		sol->fn = quat_mul_v3(sol->quat, sol->fb);
		v3_t wnin2 = v3_scalar(-dt / 2, sol->eth.wnin);
		v3_t a = sol->fn;
		rotvec(&wnin2, &a);
		sol->an = v3_add(a, sol->eth.gcc);
		v3_t vn1 = v3_add(sol->vel, v3_scalar(dt, sol->an));

		sol->eth.Mpv.m11 = 0.0;
		sol->eth.Mpv.m12 = 1.0 / sol->eth.RMh;
		sol->eth.Mpv.m13 = 0.0;
		sol->eth.Mpv.m21 = 1.0 / sol->eth.clRNh;
		sol->eth.Mpv.m22 = 0.0;
		sol->eth.Mpv.m23 = 0.0;
		sol->eth.Mpv.m31 = 0.0;
		sol->eth.Mpv.m32 = 0.0;
		sol->eth.Mpv.m33 = 1.0;
		sol->Mpvvn = m3_mul_v3(sol->eth.Mpv, v3_scalar(0.5, v3_add(sol->vel, vn1)));
		sol->pos = v3_add(sol->pos, v3_scalar(dt, sol->Mpvvn));
		sol->vel = vn1;

		rvupdquat(*dtheta, v3_scalar(dt, sol->eth.wnin), &sol->quat);
		return 1;
	}
	return 0;
}

// FRD
extern int ins_nav_llh_back(double dt, const v3_t* dtheta, const v3_t* dv, solins_t* sol, int opt) {
	/* Local variable definition ------------------------------------------- */
	m3_t theta;
	v3_t zeta, vng_cor, vbf, vnf, ndtheta = *dtheta;
	v3_t old_v = sol->vel;
	quat_t q_nk_nk1, q_bk1_bk;

	/* Code start ---------------------------------------------------------- */
	updearth(&sol->eth, sol->pos, sol->vel, opt, 0);

	/* Velocity update */
	vng_cor = v3_scalar(dt, v3_del(sol->eth.gn, v3_cross(sol->eth.wnien, sol->vel)));  // gravity and Coriolis
	vbf = v3_add(*dv, v3_scalar(0.5, v3_cross(*dtheta, *dv)));
	theta = v3_askew(zeta);
	vnf = m3_mul_v3(m3_mul(m3_del(I33, m3_scalar(0.5, theta)), sol->dcm), vbf); // specific force
	sol->vel = v3_del(v3_del(sol->vel, vng_cor), vnf);

	/* Position update */
	sol->pos.x = sol->pos.x - 0.5 * (sol->vel.x + old_v.x) / (sol->eth.RMh) * dt;
	sol->pos.y = sol->pos.y - 0.5 * (sol->vel.y + old_v.y) / (sol->eth.RNh * sol->eth.cl) * dt;
	sol->pos.z = sol->pos.z + 0.5 * (sol->vel.z + old_v.z) * dt;

	/* Attitude update */
	zeta = v3_scalar(dt, sol->eth.wnin);
	rv2quat(&zeta, &q_nk_nk1);      //update n
	ndtheta = v3_scalar(-1.0, ndtheta);
	rv2quat(&ndtheta, &q_bk1_bk);   //update b
	sol->quat = quat_mul(quat_mul(q_nk_nk1, sol->quat), q_bk1_bk);
	quat_normalize(&sol->quat);
	quat2dcm(&sol->quat, &sol->dcm, IMUCOOR_FRD);

	return 1;
}

/**
 * @brief Strapdown-INS equations under ECEF frame
 * @param[in]   dt      Time interval [s]
 * @param[in]   dtheta  Angular increment [rad]
 * @param[in]   dv      Velocity increment [m/s]
 * @param[in,out]   r   Start/End postion in ECEF(reb_e) [m]
 * @param[in,out]   v   Start/End velocity in ECEF(veb_e) [m]
 * @param[in,out]   q   Start/End attitude trans express by quaternion(qbe)
 * @return 0: OK
 * @see multisample()
 * @note This function do not contain conning&sculling error compensation,
 *      so it should work with multisample() function.
 *
 *      Ref: Paul. D. Groves. Principles of GNSS, Inertial, and Multisensor
 *      Integrated Navigation Systems(2nd Edition), P163
 */
extern int ins_nav_ecef(double dt, const v3_t* dtheta, const v3_t* dv, v3_t* r, v3_t* v, quat_t* q, int opt) {
	quat_t q_earth, q_body, old_q = *q;

	/* Attitude update */
	v3_t dtheta_ie = { 0.0, 0.0, -wgs84.wie * dt };
	rv2quat(&dtheta_ie, &q_earth);   // update e
	rv2quat(dtheta, &q_body);        // update b
	*q = quat_mul(quat_mul(q_earth, old_q), q_body);
	quat_normalize(q);

	/* Specific force transform(velocity form) */
	v3_t dtheta_ie_half = { 0.0, 0.0, -wgs84.wie * dt / 2.0 };
	rv2quat(&dtheta_ie_half, &q_earth);
	v3_t dv_rot = v3_scalar(0.5, v3_cross(*dtheta, *dv));  // rotation compensation
	v3_t dv_e = quat_mul_v3(quat_mul(q_earth, old_q), v3_add(*dv, dv_rot));  // vef

	/* Velocity update */
	v3_t ge, gn, gee;
	gravity_ecef(r, &ge);
	double pos[3];
	ecef2pos(r->v, pos);
	gravity_ned(pos[0], pos[1], &gn);
	double Cne[9];
	ned2xyz(pos, Cne);
	matmul3v("N", Cne, gn.v, gee.v);

	v3_t old_v = *v;
	v->x = old_v.x + dv_e.x + dt * (gee.x + 2.0 * wgs84.wie * old_v.y);  //gravity and Coriolis
	v->y = old_v.y + dv_e.y + dt * (gee.y - 2.0 * wgs84.wie * old_v.x);
	v->z = old_v.z + dv_e.z + dt * (gee.z);  /// Grove, p154

	/* Position update */
	v3_t old_r = *r;
	*r = v3_add(old_r, v3_scalar(0.5 * dt, v3_add(old_v, *v)));
	return 1;
}

extern int ins_nav_ecef_back(double dt, const v3_t* dtheta, const v3_t* dv, v3_t* r, v3_t* v, quat_t* q, int opt) {
	dt = fabs(dt);
	v3_t ge, ndtheta = *dtheta, old_r = *r, old_v = *v;
	quat_t q_earth, q_body, old_q = *q;

	/* Specific force transform(velocity form) */
	v3_t dtheta_ie_half = { 0.0, 0.0, -wgs84.wie * dt / 2.0 };
	rv2quat(&dtheta_ie_half, &q_earth);
	v3_t dv_rot = v3_scalar(0.5, v3_cross(*dtheta, *dv));  // rotation compensation
	v3_t dv_e = quat_mul_v3(quat_mul(q_earth, old_q), v3_add(*dv, dv_rot));  // vef

	/* Velocity update */
	gravity_ecef(r, &ge);
	v->x = old_v.x - dv_e.x - dt * (ge.x + 2.0 * wgs84.wie * old_v.y);
	v->y = old_v.y - dv_e.y - dt * (ge.y - 2.0 * wgs84.wie * old_v.x);
	v->z = old_v.z - dv_e.z - dt * ge.z;  // gravity and Coriolis

	/* Position update */
	*r = v3_del(old_r, v3_scalar(0.5 * dt, v3_add(old_v, *v)));

	/* Attitude update */
	ndtheta = v3_scalar(-1.0, ndtheta);
	v3_t dtheta_ie = { 0.0, 0.0, wgs84.wie * dt };
	rv2quat(&dtheta_ie, &q_earth);     // update e
	rv2quat(&ndtheta, &q_body);        // update b
	*q = quat_mul(quat_mul(q_earth, old_q), q_body);
	quat_normalize(q);

	return 1;
}