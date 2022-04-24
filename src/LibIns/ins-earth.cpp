/*------------------------------------------------------------------------------
 * ins-com.cc : ins common functions
 *
 * version : reorganized form chen chao's code
 * history : Created by lizhen on 2021/5/28.
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

 /**
  * @brief Caculate meridian radius of curvature(from the north to south)
  * @param[in] eth   earth parameters struct
  * @param[in] lat   latitude [rad]
  * @return meridian radius of curvature [m]
  * @see earth_RE()
  * @note
  *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
  *          Integrated Navigation Systems(2nd Edition), 2013, P59
  */
extern double earth_RN(const earth_t* eth, double lat) {
	double e2 = SQR(eth->e);
	return eth->R0 * (1 - e2) * pow(1 - e2 * SQR(sin(lat)), -1.5);
}

/**
 * @brief Caculate transverse radius of curvature(from the east to west)
 * @param[in] eth   earth parameters struct
 * @param[in] lat   latitude [rad]
 * @return transverse radius of curvature [m]
 * @see earth_RN()
 * @note Also called nomal radius of curvature, or prime vertical radius of
 *  curvature.
 *
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, P59
 */
extern double earth_RE(const earth_t* eth, double lat) {
	return eth->R0 / sqrt(1 - SQR(eth->e) * SQR(sin(lat)));
}

extern double updearth(earth_t* eth, const v3_t pos, const v3_t vel, int opt, int init)
{   // in ned
	eth->sl = sin(pos.x); eth->cl = cos(pos.x); eth->tl = eth->sl / eth->cl;
	eth->sl2 = eth->sl * eth->sl;
	double sl4 = eth->sl2 * eth->sl2;
	double sq = 1.0 - wgs84.e * wgs84.e * eth->sl2;
	double RN = wgs84.R0 / sqrt(sq);
	eth->RN = RN;
	eth->RNh = RN + pos.z;
	eth->clRNh = eth->cl * eth->RNh;
	eth->RM = RN * (1.0 - SQR(wgs84.e)) / sq;
	eth->RMh = eth->RM + pos.z;

	if (opt == INSLOCAL_NED) {
		eth->wnie.x = wgs84.wie * eth->cl;  ///£¡note now ned
		eth->wnie.y = 0.0;
		eth->wnie.z = -wgs84.wie * eth->sl;

		eth->wnen.x = vel.y / eth->RNh;   ///£¡note now ned
		eth->wnen.y = -vel.x / eth->RMh;
		eth->wnen.z = -vel.y * eth->tl / eth->RNh;

		if (init) {
			eth->Mpv.m11 = 0.0; eth->Mpv.m12 = 1.0 / eth->RMh; eth->Mpv.m13 = 0.0;
			eth->Mpv.m21 = 1.0 / eth->clRNh; eth->Mpv.m22 = 0.0; eth->Mpv.m23 = 0.0;
			eth->Mpv.m31 = 0.0; eth->Mpv.m32 = 0.0; eth->Mpv.m33 = 1.0;
		}
	}
	else if (opt == INSLOCAL_ENU) {
		eth->wnie.x = 0.0;
		eth->wnie.y = wgs84.wie * eth->cl;
		eth->wnie.z = wgs84.wie * eth->sl;

		eth->wnen.x = -vel.y / eth->RMh;
		eth->wnen.y = vel.x / eth->RNh;
		eth->wnen.z = (vel.x / eth->RNh) * eth->tl;

		if (init) {
			eth->Mpv.m11 = 0.0; eth->Mpv.m12 = 1.0 / eth->RMh; eth->Mpv.m13 = 0.0;
			eth->Mpv.m21 = 1.0 / eth->clRNh; eth->Mpv.m22 = 0.0; eth->Mpv.m23 = 0.0;
			eth->Mpv.m31 = 0.0; eth->Mpv.m32 = 0.0; eth->Mpv.m33 = 1.0;
		}
	}

	eth->wnin = v3_add(eth->wnie, eth->wnen);
	eth->wnien = v3_add(eth->wnie, eth->wnin);
	eth->g = G2MPS2 * (1.0 + 5.27094e-3 * eth->sl2 + 2.32718e-5 * sl4) - 3.086e-6 * pos.z;
	//    eth->g = G2MPS2 * (1.0 + 5.2790414e-3 * eth->sl2 + 2.32718e-5 * sl4) +
	//              (-3.087691089e-6 + 4.397731e-9 * eth->sl2) * pos.z + 7.21e-13 * pos.z * pos.z;

	eth->gn.x = eth->gn.y = 0.0; eth->gn.z = opt == INSLOCAL_NED ? eth->g : -eth->g;

	eth->gcc = v3_del(eth->gn, v3_cross(eth->wnien, vel));  // gravity and Coriolis

	return 0;
}

/**
 * @brief Gravitational acceleration of Earth project to e-frame
 * @param[in]   r   position under e-axis
 * @param[out]  ge  Gravitational acceleration under ECEF(m s^-2)
 * @return 0: OK
 * @see gravity_ned()
 * @note ge do NOT contain centrifugal force
 */
extern int gravity_ecef(const v3_t* r, v3_t* ge) {
	double mag_r = sqrt(r->x * r->x + r->y * r->y + r->z * r->z);
	if (fabs(mag_r) < 1e-32) {
		ge->x = 0.0;
		ge->y = 0.0;
		ge->z = 0.0;
	}
	else {
		/* Calculate gravitational acceleration using (2.142) */
		double f1 = -wgs84.mu / pow(mag_r, 3);
		double f2 = 1.5 * wgs84.J2 * pow(wgs84.R0 / mag_r, 2.0);
		double z_scale = 5.0 * pow((r->z / mag_r), 2.0);
		double g1 = f1 * (r->x + f2 * (1.0 - z_scale) * r->x);
		double g2 = f1 * (r->y + f2 * (1.0 - z_scale) * r->y);
		double g3 = f1 * (r->z + f2 * (3.0 - z_scale) * r->z);
		/* Add centripetal acceleration using (2.133) */
		ge->x = g1 + wgs84.wie * wgs84.wie * r->x;
		ge->y = g2 + wgs84.wie * wgs84.wie * r->y;
		ge->z = g3;
	}
	return 0;
}

/**
 * @brief Acceleration of gravity under n-axis(NED)
 * @param[in]   lat     latitude [rad]
 * @param[in]   hgt     ellipsoidal height [m]
 * @param[out]  gn      acceleration of gravity [m s^-2]
 * @return 0: OK
 * @see gravity_ecef()
 * @note gravity contains two part: gravitational and centrifugal accelration
 */
extern int gravity_ned(double lat, double hgt, v3_t* gn) {
	double sinlat2 = sin(lat) * sin(lat);
	double e2 = SQR(wgs84.e);

	/* Calculate surface gravity using Somigliana model */
	double g0 = 9.7803253359 * (1.0 + 0.001931853 * sinlat2)
		/ sqrt(1.0 - e2 * sinlat2);

	gn->x = -8.08E-9 * hgt * sin(2.0 * lat); /* North */
	gn->y = 0.0;                             /* East */
	/* Down */
	double tmp = 1.0 + wgs84.f * (1.0 - 2.0 * sinlat2)
		+ (SQR(wgs84.wie) * SQR(wgs84.R0) * wgs84.RP / wgs84.mu);
	gn->z = g0 *
		(1.0 - (2.0 / wgs84.R0) * tmp * hgt + (3.0 * SQR(hgt) / SQR(wgs84.R0)));

	return 0;
}