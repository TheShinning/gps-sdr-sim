/*------------------------------------------------------------------------------
 * ins-com.cc : ins common functions
 *
 * version : reorganized form chen chao's code
 * history : Created by lizhen on 2021/3/21.
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

extern const m3_t I33 = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }; /**< unit 3D matrix */
extern const m3_t O33 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; /**< zero 3D matrix */
extern const v3_t V0 = { 0.0, 0.0, 0.0 };      /**< zero 3D vector */
extern const v3_t V1 = { 1.0, 1.0, 1.0 };      /**< unit 3D vector */
extern const double Crf[9] = { 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0 };

earth_t wgs84 = {
	/*.wie =*/ 7.2921151467E-5, /**< WGS84 earth parameters */
	/*.R0 =*/ 6378137,
	/*.RP =*/ 6356752.31425,
	/*.mu =*/ 3.986004418E14,
	/*.J2 =*/ 1.082627E-3,
	/*.e =*/ 0.08181919084255230,
	/*.f = */1.0 / 298.257223563
};

/**
 * @brief Convert Euler attitude to Quaternion attitude(Eab => Qab)
 * @param[in] euler Input Euler attitude
 * @param[out] quat Ouput Quaternion attitude
 * @return 0: OK
 */
extern int euler2quat(const v3_t* euler, quat_t* quat, int opt) {
	double si = sin(euler->x / 2), ci = cos(euler->x / 2);
	double sj = sin(euler->y / 2), cj = cos(euler->y / 2);
	double sk = sin(euler->z / 2), ck = cos(euler->z / 2);

	if (opt == IMUCOOR_FRD) {
		quat->q0 = ci * cj * ck + si * sj * sk;
		quat->q1 = si * cj * ck - ci * sj * sk;
		quat->q2 = ci * sj * ck + si * cj * sk;
		quat->q3 = ci * cj * sk - si * sj * ck;
	}
	else if (opt == IMUCOOR_RFU) {
		quat->q0 = ci * cj * ck - si * sj * sk;
		quat->q1 = si * cj * ck - ci * sj * sk;
		quat->q2 = ci * sj * ck + si * cj * sk;
		quat->q3 = ci * cj * sk + si * sj * ck;
	}

	quat_conj(quat);
	return 0;
}

/**
 * @brief Convert Quaternion attitude to Euler attitude(Qab => Eab)
 * @param[in]  quat     Input Quaternion attitude
 * @param[out] euler    Output Euler attitue
 * @return  0: OK
 */
extern int quat2euler(const quat_t* quat, v3_t* euler, int imu_opt) {
	m3_t dcm;
	quat2dcm(quat, &dcm, imu_opt);
	dcm2att(&dcm, euler, imu_opt);
	return 0;
}

/**
 * @brief Convert DCM attitude to Euler attitude
 * @param[in]   dcm     Input DCM attitude
 * @param[out]  euler   Ouput Euler attitude
 * @return  0: OK
 */
extern int dcm2euler(const m3_t* dcm, v3_t* euler, int imu_opt) {
	if (imu_opt == IMUCOOR_FRD) {
		euler->x = atan2(dcm->m23, dcm->m33);
		euler->y = -asin(dcm->m13);
		euler->z = atan2(dcm->m12, dcm->m11);
	}
	else if (imu_opt == IMUCOOR_RFU) {
		euler->x = asin(dcm->m32);
		euler->y = atan2(-dcm->m31, dcm->m33);
		euler->z = atan2(-dcm->m12, dcm->m22);
	}

	return 0;
}

/**
 * @brief Convert Euler to DCM attitude
 * @param[in] euler     Input Euler attitude
 * @param[out] dcm      Ouput DCM attitude
 * @return  0: OK
 */
extern int euler2dcm(const v3_t* euler, m3_t* dcm, int opt) {
	double sin_phi = sin(euler->x);
	double cos_phi = cos(euler->x);
	double sin_theta = sin(euler->y);
	double cos_theta = cos(euler->y);
	double sin_psi = sin(euler->z);
	double cos_psi = cos(euler->z);

	if (opt == IMUCOOR_FRD) {
		/* Calculate coordinate transformation matrix using (2.22) */
		dcm->m11 = cos_theta * cos_psi;
		dcm->m12 = cos_theta * sin_psi;
		dcm->m13 = -sin_theta;
		dcm->m21 = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
		dcm->m22 = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
		dcm->m23 = sin_phi * cos_theta;
		dcm->m31 = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
		dcm->m32 = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
		dcm->m33 = cos_phi * cos_theta;
	}
	else if (opt == IMUCOOR_RFU) {
		dcm->m11 = cos_theta * cos_psi - sin_phi * sin_theta * sin_psi;
		dcm->m12 = -cos_phi * sin_psi;
		dcm->m13 = sin_theta * cos_psi + sin_phi * cos_theta * sin_psi;
		dcm->m21 = cos_theta * sin_psi + sin_phi * sin_theta * cos_psi;
		dcm->m22 = cos_phi * cos_psi;
		dcm->m23 = sin_theta * sin_psi - sin_phi * cos_theta * cos_psi;
		dcm->m31 = -cos_phi * sin_theta;
		dcm->m32 = sin_phi;
		dcm->m33 = cos_phi * cos_theta;
	}
	return 0;
}

/**
 * @brief Convert DCM attitude to Quaternion attitude
 * @param[in]   dcm     DCM attitude
 * @param[out]  quat    Quaternion attitude
 * @return  O: OK
 */
extern int dcm2quat(const m3_t* dcm, quat_t* quat) {
	double qq4;
	if (dcm->m11 >= dcm->m22 + dcm->m33) {
		quat->q1 = 0.5 * sqrt(1 + dcm->m11 - dcm->m22 - dcm->m33);
		qq4 = 4 * quat->q1;
		quat->q0 = (dcm->m32 - dcm->m23) / qq4;
		quat->q2 = (dcm->m12 + dcm->m21) / qq4;
		quat->q3 = (dcm->m13 + dcm->m31) / qq4;
	}
	else if (dcm->m22 >= dcm->m11 + dcm->m33) {
		quat->q2 = 0.5 * sqrt(1 - dcm->m11 + dcm->m22 - dcm->m33);
		qq4 = 4 * quat->q2;
		quat->q0 = (dcm->m13 - dcm->m31) / qq4;
		quat->q1 = (dcm->m12 + dcm->m21) / qq4;
		quat->q3 = (dcm->m23 + dcm->m32) / qq4;
	}
	else if (dcm->m33 >= dcm->m11 + dcm->m22) {
		quat->q3 = 0.5 * sqrt(1 - dcm->m11 - dcm->m22 + dcm->m33);
		qq4 = 4 * quat->q3;
		quat->q0 = (dcm->m21 - dcm->m12) / qq4;
		quat->q1 = (dcm->m13 + dcm->m31) / qq4;
		quat->q2 = (dcm->m23 + dcm->m32) / qq4;
	}
	else {
		quat->q0 = 0.5 * sqrt(1 + dcm->m11 + dcm->m22 + dcm->m33);
		qq4 = 4 * quat->q0;
		quat->q1 = (dcm->m32 - dcm->m23) / qq4;
		quat->q2 = (dcm->m13 - dcm->m31) / qq4;
		quat->q3 = (dcm->m21 - dcm->m12) / qq4;
	}
	quat_normalize(quat);
	return 0;
}

/**
 * @brief Convert Quaternion to DCM attitude
 * @param[in]   quat    Input Quaternion attitude
 * @param[out]  dcm     Oput DCM attitude
 * @return 0: OK
 */
extern int quat2dcm(const quat_t* quat, m3_t* dcm, int opt) {
	double q11 = quat->q0 * quat->q0, q12 = quat->q0 * quat->q1,
		q13 = quat->q0 * quat->q2, q14 = quat->q0 * quat->q3,
		q22 = quat->q1 * quat->q1, q23 = quat->q1 * quat->q2,
		q24 = quat->q1 * quat->q3, q33 = quat->q2 * quat->q2,
		q34 = quat->q2 * quat->q3, q44 = quat->q3 * quat->q3;

	if (opt == IMUCOOR_FRD) {
		dcm->m11 = q11 + q22 - q33 - q44;
		dcm->m12 = 2 * (q23 - q14);
		dcm->m13 = 2 * (q24 + q13);
		dcm->m21 = 2 * (q23 + q14);
		dcm->m22 = q11 - q22 + q33 - q44;
		dcm->m23 = 2 * (q34 - q12);
		dcm->m31 = 2 * (q24 - q13);
		dcm->m32 = 2 * (q34 + q12);
		dcm->m33 = q11 - q22 - q33 + q44;
	}
	else if (opt == IMUCOOR_RFU) {
		dcm->m11 = q11 + q22 - q33 - q44;
		dcm->m12 = 2 * (q23 + q14);
		dcm->m13 = 2 * (q24 - q13);
		dcm->m21 = 2 * (q23 - q14);
		dcm->m22 = q11 - q22 + q33 - q44;
		dcm->m23 = 2 * (q34 + q12);
		dcm->m31 = 2 * (q24 + q13);
		dcm->m32 = 2 * (q34 - q12);
		dcm->m33 = q11 - q22 - q33 + q44;
	}

	return 0;
}

/**
 * @brief Convert Euler attitude(roll, pitch, yaw) to DCM
 * @param[in]   att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @param[out]  dcm     DCM, Cbn
 * @see euler2dcm() dcm2att()
 * @return 0: OK
 */
extern int att2dcm(const v3_t* att, m3_t* dcm, int opt) {
	euler2dcm(att, dcm, opt);
	*dcm = m3_T(*dcm);
	return 0;
}

/**
 * @brief Convert DCM to Euler attitude(roll, pitch, yaw)
 * @param[in]   dcm     DCM, Cbn
 * @param[out]  att     Euler attitude[roll, pitch, yaw], Enb[rad]
 * @see dcm2euler() att2dcm()
 * @return 0: OK
 */
extern int dcm2att(const m3_t* dcm, v3_t* att, int imu_opt) {
	m3_t M = m3_T(*dcm);
	dcm2euler(&M, att, imu_opt);
	return 0;
}

/**
 * @brief Convert Euler attitude to Quaternion
 * @param[in]   att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @param[out]  quat    Quaternion, Qbn
 * @see euler2quat() quat2att()
 * @return O: OK
 */
extern int att2quat(const v3_t* att, quat_t* quat, int opt) {
	euler2quat(att, quat, opt);
	quat_conj(quat);
	return 0;
}

/**
 * @brief Convert Quaternion to Euler attitude
 * @param[in]   quat    Quaternion, Qbn
 * @param[out]  att     Euler attitude(roll, pitch, yaw), Enb[rad]
 * @see quat2euler() att2quat()
 * @return 0: OK
 */
extern int quat2att(const quat_t* quat, v3_t* att, int opt) {
	quat_t Q = *quat;
	quat_conj(&Q);
	quat2euler(&Q, att, opt);
	return 0;
}

extern int rotvec(const v3_t* rv, v3_t* vi) {
	double n, n2 = SQR(v3_norm(*rv)), n_2, q1, s;
	if (n2 < 1.0e-8) {
		q1 = 1.0 - n2 * (1.0 / 8 - n2 / 384.0);
		s = 1.0 / 2.0 - n2 * (1.0 / 48.0 - n2 / 3840.0);
	}
	else {
		n = sqrt(n2);
		n_2 = n / 2;
		q1 = cos(n_2);
		s = sin(n_2) / n;
	}
	double q2, q3, q4;
	q2 = s * rv->x;
	q3 = s * rv->y;
	q4 = s * rv->z;
	double qo1, qo2, qo3, qo4;
	qo1 = -q2 * vi->x - q3 * vi->y - q4 * vi->z;
	qo2 = q1 * vi->x + q3 * vi->z - q4 * vi->y;
	qo3 = q1 * vi->y + q4 * vi->x - q2 * vi->z;
	qo4 = q1 * vi->z + q2 * vi->y - q3 * vi->x;
	vi->x = -qo1 * q2 + qo2 * q1 - qo3 * q4 + qo4 * q3;
	vi->y = -qo1 * q3 + qo3 * q1 - qo4 * q2 + qo2 * q4;
	vi->z = -qo1 * q4 + qo4 * q1 - qo2 * q3 + qo3 * q2;

	return 0;
}

extern int rvupdquat(const v3_t rv_ib, const v3_t rv_in, quat_t* q) {
	v3_t rv1, rv2;
	/*rv2q(rv_ib)*/
	double n, n2 = SQR(v3_norm(rv_ib)), n_2, rv_ib0, rv_in0, s;
	if (n2 < 1.0e-8) {
		rv_ib0 = 1.0 - n2 * (1.0 / 8.0 - n2 / 384.0);
		s = 1.0 / 2.0 - n2 * (1.0 / 48.0 - n2 / 3840.0);
	}
	else {
		n = sqrt(n2);
		n_2 = n / 2.0;
		rv_ib0 = cos(n_2);
		s = sin(n_2) / n;
	}
	rv1 = v3_scalar(s, rv_ib);
	double qb1, qb2, qb3, qb4;
	qb1 = q->q0 * rv_ib0 - q->q1 * rv1.x - q->q2 * rv1.y - q->q3 * rv1.z;
	qb2 = q->q0 * rv1.x + q->q1 * rv_ib0 + q->q2 * rv1.z - q->q3 * rv1.y;
	qb3 = q->q0 * rv1.y + q->q2 * rv_ib0 + q->q3 * rv1.x - q->q1 * rv1.z;
	qb4 = q->q0 * rv1.z + q->q3 * rv_ib0 + q->q1 * rv1.y - q->q2 * rv1.x;

	/*rv2q(-rv_in)*/
	n2 = SQR(v3_norm(rv_in));
	if (n2 < 1.0e-8) {
		rv_in0 = 1.0 - n2 * (1.0 / 8.0 - n2 / 384.0);
		s = -1.0 / 2.0 + n2 * (1.0 / 48.0 - n2 / 3840.0);
	}
	else {
		n = sqrt(n2);
		n_2 = n / 2;
		rv_in0 = cos(n_2);
		s = -sin(n_2) / n;
	}
	rv2 = v3_scalar(s, rv_in);
	quat_t q1 = *q;
	q1.q0 = rv_in0 * qb1 - rv2.x * qb2 - rv2.y * qb3 - rv2.z * qb4;
	q1.q1 = rv_in0 * qb2 + rv2.x * qb1 + rv2.y * qb4 - rv2.z * qb3;
	q1.q2 = rv_in0 * qb3 + rv2.y * qb1 + rv2.z * qb2 - rv2.x * qb4;
	q1.q3 = rv_in0 * qb4 + rv2.z * qb1 + rv2.x * qb3 - rv2.y * qb2;

	n2 = SQR(quat_norm(q1));
	if (n2 > 1.000001 || n2 < 0.999999) {
		double nq = 1 / sqrt(n2);
		q->q0 = q1.q0 * nq;
		q->q1 = q1.q1 * nq;
		q->q2 = q1.q2 * nq;
		q->q3 = q1.q3 * nq;
	}
	else {
		*q = q1;
	}
	return 0;
}

extern void quatdelphi(quat_t* q, const v3_t* phi) {
	quat_t q0;
	rv2quat(phi, &q0);
	*q = quat_mul(q0, *q);
}

extern void attsync(const v3_t* rr, v3_t* rpy, m3_t* dcm, quat_t* quat, int imu_opt, int local_opt, int mech_opt, int att_fmt)
{
	if (rpy && att_fmt == 0) {
		att2dcm(rpy, dcm, imu_opt);
		att2quat(rpy, quat, imu_opt);
	}
	if (dcm && att_fmt == 1) {
		if (mech_opt == INSMECH_ECEF) {
			v3_t re = *rr;
			ecef2llh(&re, nullptr, dcm, local_opt);
		}
		dcm2att(dcm, rpy, imu_opt);
		dcm2quat(dcm, quat);
	}
	if (quat && att_fmt == 2) {
		quat2att(quat, rpy, imu_opt);
		quat2dcm(quat, dcm, imu_opt);
	}
}

/**
 * @brief Quaternion normalization and make first element large than zero
 * @param[in,out] quat  Quaternion/unit Quaternion
 * @return  O: OK
 */
extern int quat_normalize(quat_t* quat) {
	double nq = 1.0 / quat_norm(*quat);
	quat->q0 *= nq;
	quat->q1 *= nq;
	quat->q2 *= nq;
	quat->q3 *= nq;
	if (quat->q0 < 0) {
		quat->q0 = -quat->q0;
		quat->q1 = -quat->q1;
		quat->q2 = -quat->q2;
		quat->q3 = -quat->q3;
	}
	return 0;
}

/**
 * @brief Quaternion conjugation
 * @param[in,out] quat  Quaternion/Quaternion conjugation
 * @return 0: OK
 * @see quat_inv()
 * @note Unit Quaternion Conjugation is the same as inversion, but conjugation
 *  is much faster
 */
extern int quat_conj(quat_t* quat) {
	quat->q1 = -quat->q1;
	quat->q2 = -quat->q2;
	quat->q3 = -quat->q3;
	return 0;
}

/**
 * @brief Quaternion inversion
 * @param[in,out] quat  Input Quaternion/Ouput Quaternion invertion
 * @see quat_dot()
 * @see quat_conj()
 * @return 0: OK, 1: input Quaternion norm too small, less than sqrt(EPS)
 * @note Unit Quaternion Conjugation is the same as inversion, but conjugation
 *  is much faster
 */
extern int quat_inv(quat_t* quat) {
	double inv = 1.0 / quat_dot(*quat, *quat);
	if (fabs(inv) < EPS) return 1;

	quat_conj(quat);
	quat->q0 *= inv;
	quat->q1 *= inv;
	quat->q2 *= inv;
	quat->q3 *= inv;
	return 0;
}

/**
 * @brief Quaterion dot production, sum of corresponding elements product
 * @param[in] P     Fisrt Quaternion
 * @param[in] Q     Second Quaternion
 * @return dot product result number
 * @note dot production also called inner porduction
 */
extern double quat_dot(quat_t P, quat_t Q) {
	return P.q0 * Q.q0 + P.q1 * Q.q1 + P.q2 * Q.q2 + P.q3 * Q.q3;
}

/**
 * @brief Quaternion cross product, result of the vector cross procdut of two
 *      quaternions' imaginary part(q1,q2,q3)
 * @param[in] P     Fisrt quaternion
 * @param[in] Q     Second quatenion
 * @return Quaternion cross product result vector
 */
extern v3_t quat_cross(quat_t P, quat_t Q) {
	v3_t P_t, Q_t;
	P_t.x = P.q1; P_t.y = P.q2; P_t.z = P.q3;
	Q_t.x = Q.q1; Q_t.y = Q.q2; Q_t.z = Q.q3;

	//return v3_cross((v3_t) { P.q1, P.q2, P.q3 }, (v3_t) { Q.q1, Q.q2, Q.q3 });
	return v3_cross(P_t, Q_t);
}

/**
 * @brief L2-norm of quaternion(also call modulus of quaternion)
 * @param[in] P Input quaternion
 * @return L2-norm of quaternion
 */
extern double quat_norm(quat_t P) {
	return sqrt(P.q0 * P.q0 + P.q1 * P.q1 + P.q2 * P.q2 + P.q3 * P.q3);
}

/**
 * @brief Check two quaternions are equal or not
 * @param[in] P     Fisrt quaternion
 * @param[in] Q     Second quaternion
 * @param[in] eps   Zero threshold, e.g. 1E-10
 * @return true: P eqaul to Q, false: P not equal to Q
 * @warning minus quaternion do NOT equal to origin quaternion
 * @note    two quaternions' corresponding elements difference should be less
 *  than zero threshold, then function return true, however, minus quaterions
 *  are the same as origin if quaternion represent attitude, but this function
 *  can not identify.
 */
extern int quat_equal(const quat_t* P, const quat_t* Q, double eps) {
	double diff;
	const double* pP = (const double*)P;
	const double* pQ = (const double*)Q;
	for (int i = 0; i < 4; ++i) {
		if ((diff = pP[i] - pQ[i]) > eps || diff < -eps)
			return 0;
	}
	return 1;
}

extern quat_t qaddphi(const quat_t* quat, v3_t phi) {
	quat_t q, qpb;
	phi = v3_scalar(-1.0, phi);
	rv2quat(&phi, &q);
	qpb = quat_mul(q, *quat);
	return qpb;
}

/**
 * @brief Rotation vector(Angular increment) to quaternion attitude
 *      transformation(New to Old)
 * @param[in]   dtheta  Rotation vector(Angular increment) [rad]
 * @param[out]  quat    Quaternion attitude transformation
 * @return O: OK
 * @see rv2dcm()
 * @note Quaternion attitude transforman can be expressed as Qbi+_bi-, And
 *      fellow the compluting: (Qbi+  = Qbi-  * Qbi+_bi- )
 *      Ref: Yan Gongming, 捷联惯导算法与组合导航原理讲义, 2016, P29
 */
extern int rv2quat(const v3_t* dtheta, quat_t* quat) {
	const double F1 = 2 * 1; // define Fk = 2^k * k!
	const double F2 = F1 * 2 * 2;
	const double F3 = F2 * 2 * 3;
	const double F4 = F3 * 2 * 4;
	const double F5 = F4 * 2 * 5;
	double n2 = dtheta->x * dtheta->x + dtheta->y * dtheta->y + dtheta->z * dtheta->z;
	double f;
	if (n2 < 1.0e-08) {
		double n4 = n2 * n2;
		quat->q0 = 1.0 - n2 * (1.0 / F2) + n4 * (1.0 / F4);
		f = 0.5 - n2 * (1.0 / F3) + n4 * (1.0 / F5);
	}
	else {
		double n_2 = sqrt(n2) / 2.0;
		quat->q0 = cos(n_2);
		f = sin(n_2) / n_2 * 0.5;
	}
	quat->q1 = f * dtheta->x;
	quat->q2 = f * dtheta->y;
	quat->q3 = f * dtheta->z;
	return 0;
}

/**
 * @brief Rotation vector(Angular increment) to DCM attitude transformation
 *      (New to Old)
 * @param[in]   dtheta  Rotation vector(Angular increment) [rad]
 * @param[out]  dcm     DCM attitude transformation
 * @return  0: OK
 * @see rv2quat()
 * @note DCM attitude transforman can be expressed as Cbi+_bi-, And  fellow the
 *      compluting: (Cbi+  = Cbi-  * Cbi+_bi- )
 *      The DCM can also be expressed in "exp(skew-symmetric(dtheta))".
 *
 *      Ref: Paul D. Groves, Principles of GNSS, Inertial, and Multisensor
 *          Integrated Navigation Systems(2nd Edition), 2013, P184
 */
extern int rv2dcm(const v3_t* dtheta, m3_t* dcm) {
	double norm_dtheta = v3_norm(*dtheta);
	m3_t Ax;
	asymmetric_mat(dtheta, &Ax);
	if (norm_dtheta > 1e-8) {
		m3_t m1 = m3_scalar(sin(norm_dtheta) / norm_dtheta, Ax);
		m3_t m2 = m3_scalar((1 - cos(norm_dtheta)) / SQR(norm_dtheta), m3_mul(Ax, Ax));
		*dcm = m3_add(m3_add(I33, m1), m2);
	}
	else {
		*dcm = m3_add(I33, Ax);
	}
	return 0;
}

/**
 * @brief Convert Quaternion attitude to Rotation vector(Angular increment)
 * @param[in]  quat     Quaternion attitude trasnformation
 * @param[out] dtheta   Rotation vector(Angular increment) [rad]
 * @return 0: OK
 * @see dcm2rv()
 * @note  q = [ cos(|rv|/2); sin(|rv|/2)/|rv|*rv ]
 */
extern int quat2rv(const quat_t* quat, v3_t* dtheta) {
	double n2 = acos(fabs(quat->q0));
	double k;
	if (n2 > 1e-40) k = 2 * n2 / sin(n2); else k = 2.0;
	if (quat->q0 < 0.0) k = -k;
	dtheta->x = k * quat->q1;
	dtheta->y = k * quat->q2;
	dtheta->z = k * quat->q3;
	return 0;
}

/**
 * @brief Convert DCM attitude to Rotation vector(Angular increment)
 * @param[in]   dcm     DCM attitude transformation
 * @param[out]  dtheta  Rotation vector(Angular increment) [rad]
 * @see quat2rv()
 * @return 0: OK
 * @note
 *      DCM = I + sin(|rv|)/|rv|*(rvx) + [1-cos(|rv|)]/|rv|^2*(rvx)^2
 *      rvx: askew matrix of rv
 */
extern int dcm2rv(const m3_t* dcm, v3_t* dtheta) {
	quat_t quat;
	dcm2quat(dcm, &quat);
	quat2rv(&quat, dtheta);
	/* TODO: maybe need to iterate to get more accurate solution */
	return 0;
}

/**
 * @brief Hamiton product, is determined by the products of the basis elements
 *  and the ditributive law. Hamition product could represent the ratation
 *  quaternions
 * @param[in] P Fisrt quaternion
 * @param[in] Q Second quaternion
 * @return Hamiton product result quaternion (P * Q)
 */
extern quat_t quat_mul(quat_t P, quat_t Q) {
	quat_t qtmp;
	qtmp.q0 = P.q0 * Q.q0 - P.q1 * Q.q1 - P.q2 * Q.q2 - P.q3 * Q.q3;
	qtmp.q1 = P.q0 * Q.q1 + P.q1 * Q.q0 + P.q2 * Q.q3 - P.q3 * Q.q2;
	qtmp.q2 = P.q0 * Q.q2 + P.q2 * Q.q0 + P.q3 * Q.q1 - P.q1 * Q.q3;
	qtmp.q3 = P.q0 * Q.q3 + P.q3 * Q.q0 + P.q1 * Q.q2 - P.q2 * Q.q1;
	return qtmp;
}

/**
 * @brief Rotate a 3D vector, to make (vn = Qbn * vb) work
 * @param[in] quat  Input Quaternion
 * @param[in] vec   Input 3D vector
 * @return result of rotation 3D vector ( quat * vec )
 */
extern v3_t quat_mul_v3(quat_t quat, v3_t vec) {
	quat_t qtmp;
	v3_t vtmp;
	qtmp.q0 = -quat.q1 * vec.x - quat.q2 * vec.y - quat.q3 * vec.z;
	qtmp.q1 = quat.q0 * vec.x + quat.q2 * vec.z - quat.q3 * vec.y;
	qtmp.q2 = quat.q0 * vec.y + quat.q3 * vec.x - quat.q1 * vec.z;
	qtmp.q3 = quat.q0 * vec.z + quat.q1 * vec.y - quat.q2 * vec.x;
	vtmp.x = -qtmp.q0 * quat.q1 + qtmp.q1 * quat.q0 - qtmp.q2 * quat.q3
		+ qtmp.q3 * quat.q2;
	vtmp.y = -qtmp.q0 * quat.q2 + qtmp.q2 * quat.q0 - qtmp.q3 * quat.q1
		+ qtmp.q1 * quat.q3;
	vtmp.z = -qtmp.q0 * quat.q3 + qtmp.q3 * quat.q0 - qtmp.q1 * quat.q2
		+ qtmp.q2 * quat.q1;
	return vtmp;
}

/**
 * @brief vector outer product/cross product
 * @param[in] v1    First vector
 * @param[in] v2    Second vector
 * @return cross product result of v1 and v2
 */
extern v3_t v3_cross(v3_t v1, v3_t v2) {
	v3_t v;
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}

/**
 * @brief Sum of the correspoding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return sum of v1 and v2 (v1 + v2)
 */
extern v3_t v3_add(v3_t v1, v3_t v2) {
	v3_t temp_t;
	temp_t.x = v1.x + v2.x;
	temp_t.y = v1.y + v2.y;
	temp_t.z = v1.z + v2.z;
	//return (v3_t) { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	return temp_t;
}

/**
 * @brief Substraction of corresponding elements of two vector
 * @param[in] v1    Fisrt vector
 * @param[in] v2    Second vector
 * @return Substraction of v1 and v2 ( v1 - v2 )
 */
extern v3_t v3_del(v3_t v1, v3_t v2) {
	//return (v3_t) { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	v3_t temp_t;
	temp_t.x = v1.x - v2.x;
	temp_t.y = v1.y - v2.y;
	temp_t.z = v1.z - v2.z;
	return temp_t;
}

/**
 * @brief Scalar multiplication between number and vector
 * @param[in] s     Input number
 * @param[in] v     Input vector
 * @return Scalar muliplication result of s and v ( s x v )
 */
extern v3_t v3_scalar(double s, v3_t v) {
	//return (v3_t) { s* v.x, s* v.y, s* v.z };
	v3_t temp_t;
	temp_t.x = s * v.x;
	temp_t.y = s * v.y;
	temp_t.z = s * v.z;
	return temp_t;
}

/**
 * @brief 3D vector normalize(to an unit vector)
 * @param[in,out] v input & output vector
 * @return  0: OK, 1: Input vector are zero vector
 * @see v3_norm()
 * @warning  if input vector length less than EPS(Zero determination threshold),
 *      v3_normalize() will do nothing.
 */
extern int v3_normalize(v3_t* v) {
	double norm = v3_norm(*v);
	if (fabs(norm) < EPS) return 1;
	double fac = 1.0 / norm;
	v->x *= fac;
	v->y *= fac;
	v->z *= fac;
	return 0;
}

/**
 * @brief L2-norm/Euclidean Metric/Euclidean Distance of vector
 * @param[in] v Input vector
 * @return L2-norm of v
 * @note Ref: http://mathworld.wolfram.com/VectorNorm.html
 */
extern double v3_norm(v3_t v) {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief row vector multiply by column vector, also called vector dot product
 * @param[in] v1    row vector
 * @param[in] v2    column vector
 * @see v3_mul_cxr() v3_mul_cxc() v3_cross()
 * @return product result number of v1 and v2
 */
extern double v3_mul_rxc(v3_t v1, v3_t v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/**
 * @brief column vector multiply by row vector
 * @param[in] v1    column vector
 * @param[in] v2    row vector
 * @return product result matrix of v1 and v2
 * @see v3_mul_rxc() v3_mul_cxc() v3_cross()
 */
extern m3_t v3_mul_cxr(v3_t v1, v3_t v2) {
	//return (m3_t) {
	//    v1.x* v2.x,
	//    v1.x* v2.y,
	//    v1.x* v2.z,
	//    v1.y* v2.x,
	//    v1.y* v2.y,
	//    v1.y* v2.z,
	//    v1.z* v2.x,
	//    v1.z* v2.y,
	//    v1.z* v2.z
	//};

	m3_t temp_t;
	temp_t.m11 = v1.x * v2.x;
	temp_t.m12 = v1.x * v2.y;
	temp_t.m13 = v1.x * v2.z;
	temp_t.m21 = v1.y * v2.x;
	temp_t.m22 = v1.y * v2.y;
	temp_t.m23 = v1.y * v2.z;
	temp_t.m31 = v1.z * v2.x;
	temp_t.m32 = v1.z * v2.y;
	temp_t.m33 = v1.z * v2.z;
	return temp_t;
}

/**
 * @brief multiply coresponding elements of two vector
 * @param[in] v1    vector 1
 * @param[in] v2    vector 2
 * @return result vector
 * @see v3_mul_cxr() v3_mul_rxc() v3_cross()
 */
extern v3_t v3_mul_cxc(v3_t v1, v3_t v2) {
	v3_t temp_t;
	temp_t.x = v1.x * v2.x;
	temp_t.y = v1.y * v2.y;
	temp_t.z = v1.z * v2.z;
	return temp_t;

	//return (v3_t) { v1.x* v2.x, v1.y* v2.y, v1.z* v2.z };
}

/**
 * @brief Form 3D diagonal matrix by using 3D vector
 * @param[in] v Input vector
 * @return  3D diagonal matrix
 */
extern m3_t v3_diag(v3_t v) {
	//return (m3_t) { v.x, 0.0, 0.0, 0.0, v.y, 0.0, 0.0, 0.0, v.z };
	m3_t temp_t;
	temp_t.m11 = v.x;
	temp_t.m12 = 0.0;
	temp_t.m13 = 0.0;
	temp_t.m21 = 0.0;
	temp_t.m22 = v.y;
	temp_t.m23 = 0.0;
	temp_t.m31 = 0.0;
	temp_t.m32 = 0.0;
	temp_t.m33 = v.z;
	return temp_t;
}

/**
 * @brief Power exponent of every elements in vector
 * @param[in] v     Input vector
 * @param[in] order Power order, e.g. 2.0 means square, 0.5 means root square
 * @return Power exponent result vector
 */
extern v3_t v3_pow(v3_t v, double order) {
	//return (v3_t) { pow(v.x, order), pow(v.y, order), pow(v.z, order) };
	v3_t temp_t;
	temp_t.x = pow(v.x, order);
	temp_t.y = pow(v.y, order);
	temp_t.z = pow(v.z, order);
	return temp_t;
}

/**
 * @brief Check two 3D vector are equal or not
 * @param[in] v1    Fisrt 3D vector
 * @param[in] v2    Second 3D vector
 * @param[in] eps   Zero threshold, e.g. 1E-10
 * @return true: v1 eqaul to v2, false: v1 not equal to v2
 * @note    two vector's corresponding elements difference should be less than
 *  zero threshold, then function return true
 */
extern int v3_equal(const v3_t* v1, const v3_t* v2, double eps) {
	double diff;
	const double* pv1 = (const double*)v1;
	const double* pv2 = (const double*)v2;
	for (int i = 0; i < 3; ++i) {
		if ((diff = pv1[i] - pv2[i]) > eps || diff < -eps)
			return 0;
	}
	return 1;
}

/**
 * @brief Square root of every elements in vector
 * @param[in] v     Input vector
 * @return Square root result vector
 */
extern v3_t v3_sqrt(v3_t v) {
	//return (v3_t) { sqrt(v.x), sqrt(v.y), sqrt(v.z) };

	v3_t temp_t;
	temp_t.x = sqrt(v.x);
	temp_t.y = sqrt(v.y);
	temp_t.z = sqrt(v.z);
	return temp_t;
}

/**
 * @brief sum of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  sum vector
 */
extern v3_t v3_sum(const v3_t* v3_list, int n) {
	v3_t vsum = V0;
	for (int i = 0; i < n; ++i) {
		vsum.x += v3_list[i].x;
		vsum.y += v3_list[i].y;
		vsum.z += v3_list[i].z;
	}
	return vsum;
}

/**
 * @brief mean value of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return   mean value vector
 */
extern v3_t v3_mean(const v3_t* v3_list, int n) {
	return v3_scalar(1.0 / n, v3_sum(v3_list, n));
}

/**
 * @brief root-mean-square value of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  root-mean-square vector
 */
extern v3_t v3_rms(const v3_t* v3_list, int n) {
	v3_t vsum = V0;
	for (int i = 0; i < n; ++i) {
		vsum.x += SQR(v3_list[i].x);
		vsum.y += SQR(v3_list[i].y);
		vsum.z += SQR(v3_list[i].z);
	}
	return v3_sqrt(v3_scalar(1.0 / n, vsum));
}

/**
 * @brief standard deviation of 3D vector list's corresponding elements
 * @param[in] v3_list   Input vector list
 * @param[in] n         number of vector
 * @return  standard deviation vector
 */
extern v3_t v3_std(const v3_t* v3_list, int n) {
	v3_t vmean = v3_mean(v3_list, n);
	v3_t vsum = V0;
	for (int i = 0; i < n; ++i) {
		vsum.x += SQR(v3_list[i].x - vmean.x);
		vsum.y += SQR(v3_list[i].y - vmean.y);
		vsum.z += SQR(v3_list[i].z - vmean.z);
	}
	return v3_sqrt(v3_scalar(1.0 / (n - 1), vsum));
}

extern void v3_2_array(const v3_t v, double* a) {
	*(a + 0) = v.x;
	*(a + 1) = v.y;
	*(a + 2) = v.z;
}
/**
 * @brief  get anti-symmetric matrix of a 3D vector
 * @param[in] v3    Input 3D vector
 * @param[in] mat   result anti-symmetric matrix
 * @return 0: OK
 */
extern int asymmetric_mat(const v3_t* v3, m3_t* mat) {
	mat->m11 = 0.0;
	mat->m12 = -v3->z;
	mat->m13 = v3->y;
	mat->m21 = v3->z;
	mat->m22 = 0.0;
	mat->m23 = -v3->x;
	mat->m31 = -v3->y;
	mat->m32 = v3->x;
	mat->m33 = 0.0;
	return 0;
}

/**
 * @brief Convert 3D vector to 3D skew symmetric matrix.
 * @param[in] v     3D vector
 * @return skew symmetric matrix
 * @see m3_iaskew()
 */
extern m3_t v3_askew(v3_t v) {
	m3_t Vx = O33;
	Vx.m12 = -v.z;
	Vx.m13 = v.y;
	Vx.m21 = v.z;
	Vx.m23 = -v.x;
	Vx.m31 = -v.y;
	Vx.m32 = v.x;
	return Vx;
}

/**
 * @brief Convert 3D skew symmetric matrix to 3D vector, the inverse precedure
	of function v3_askew().
 * @param[in] m     3D skew symmetric matrix
 * @return correspoding 3D vector
 * @see v3_askew()
 */
extern v3_t m3_iaskew(m3_t m) {
	v3_t v;
	v.x = (m.m32 - m.m23) / 2.0;
	v.y = (m.m13 - m.m31) / 2.0;
	v.z = (m.m21 - m.m12) / 2.0;
	return v;
}

/**
 * @brief 3D matrix transposition
 * @param[in] A Input 3D matrix
 * @return transposition matrix of A
 */
extern m3_t m3_T(m3_t A) {
	m3_t dcm;
	dcm.m11 = A.m11;
	dcm.m12 = A.m21;
	dcm.m13 = A.m31;
	dcm.m21 = A.m12;
	dcm.m22 = A.m22;
	dcm.m23 = A.m32;
	dcm.m31 = A.m13;
	dcm.m32 = A.m23;
	dcm.m33 = A.m33;
	return dcm;
}

/**
 * @brief 3D matrix inverse
 * @param[in,out] A Input 3D matrix, output inverse matrix of A
 * @return 0: OK
 */
extern int m3_inv(m3_t* A) {
	/* PLU decomposition */
	m3_t P, L, U;
	if (!m3_LU(A, &L, &U, &P)) {
		/* inverse of L */
		L.m21 = L.m21 / (L.m22 * L.m11);
		L.m31 = L.m31 / (L.m33 * L.m11) - L.m32 * L.m21 / L.m33;
		L.m22 = -1.0 / L.m22;
		L.m32 = -L.m32 * L.m22 / L.m33;
		L.m11 = -1.0 / L.m11;
		L.m33 = -1.0 / L.m33;

		/* inverse of U */
		U.m12 = U.m12 / (U.m22 * U.m11);
		U.m13 = U.m13 / (U.m33 * U.m11) - U.m23 * U.m12 / U.m33;
		U.m22 = -1.0 / U.m22;
		U.m23 = -U.m23 * U.m22 / U.m33;
		U.m11 = -1.0 / U.m11;
		U.m33 = -1.0 / U.m33;
	}
	*A = m3_mul(U, m3_mul(L, P));
	return 0;
}

/**
 * @brief Sum of corresponding elements of two 3D matrix
 * @param[in] A Fisrt matrix
 * @param[in] B Second matrix
 * @return Sum of A and B (A + B)
 */
extern m3_t m3_add(m3_t A, m3_t B) {
	//return (m3_t) {
	//    A.m11 + B.m11, A.m12 + B.m12, A.m13 + B.m13,
	//        A.m21 + B.m21, A.m22 + B.m22, A.m23 + B.m23,
	//        A.m31 + B.m31, A.m32 + B.m32, A.m33 + B.m33
	//};

	m3_t temp_t;
	temp_t.m11 = A.m11 + B.m11;
	temp_t.m12 = A.m12 + B.m12;
	temp_t.m13 = A.m13 + B.m13;
	temp_t.m21 = A.m21 + B.m21;
	temp_t.m22 = A.m22 + B.m22;
	temp_t.m23 = A.m23 + B.m23;
	temp_t.m31 = A.m31 + B.m31;
	temp_t.m32 = A.m32 + B.m32;
	temp_t.m33 = A.m33 + B.m33;
	return temp_t;
}

/**
 * @brief Substraction of corresponding elements of two matrix
 * @param[in] A Frist 3D matrix
 * @param[in] B Second 3D matrix
 * @return A minus B result (A - B)
 */
extern m3_t m3_del(m3_t A, m3_t B) {
	//return (m3_t) {
	//    A.m11 - B.m11, A.m12 - B.m12, A.m13 - B.m13,
	//        A.m21 - B.m21, A.m22 - B.m22, A.m23 - B.m23,
	//        A.m31 - B.m31, A.m32 - B.m32, A.m33 - B.m33
	//};

	m3_t temp_t;
	temp_t.m11 = A.m11 - B.m11;
	temp_t.m12 = A.m12 - B.m12;
	temp_t.m13 = A.m13 - B.m13;
	temp_t.m21 = A.m21 - B.m21;
	temp_t.m22 = A.m22 - B.m22;
	temp_t.m23 = A.m23 - B.m23;
	temp_t.m31 = A.m31 - B.m31;
	temp_t.m32 = A.m32 - B.m32;
	temp_t.m33 = A.m33 - B.m33;
	return temp_t;
}

/**
 * @brief Scalar multiplication between number and 3D matrix
 * @param[in] alpha Input number
 * @param[in] A     Input 3D matrix
 * @return Scalar muliplication result of s and A ( s x v )
 */
extern m3_t m3_scalar(double alpha, m3_t A) {
	m3_t mat;
	mat.m11 = alpha * A.m11;
	mat.m12 = alpha * A.m12;
	mat.m13 = alpha * A.m13;
	mat.m21 = alpha * A.m21;
	mat.m22 = alpha * A.m22;
	mat.m23 = alpha * A.m23;
	mat.m31 = alpha * A.m31;
	mat.m32 = alpha * A.m32;
	mat.m33 = alpha * A.m33;
	return mat;
}

/**
 * @brief matrix multiplication
 * @param[in] A Fisrt 3D matrix
 * @param[in] B Second 3D matrix
 * @return matrix multiplication result of A and B (A * B)
 */
extern m3_t m3_mul(m3_t A, m3_t B) {
	m3_t C;
	C.m11 = A.m11 * B.m11 + A.m12 * B.m21 + A.m13 * B.m31;
	C.m12 = A.m11 * B.m12 + A.m12 * B.m22 + A.m13 * B.m32;
	C.m13 = A.m11 * B.m13 + A.m12 * B.m23 + A.m13 * B.m33;
	C.m21 = A.m21 * B.m11 + A.m22 * B.m21 + A.m23 * B.m31;
	C.m22 = A.m21 * B.m12 + A.m22 * B.m22 + A.m23 * B.m32;
	C.m23 = A.m21 * B.m13 + A.m22 * B.m23 + A.m23 * B.m33;
	C.m31 = A.m31 * B.m11 + A.m32 * B.m21 + A.m33 * B.m31;
	C.m32 = A.m31 * B.m12 + A.m32 * B.m22 + A.m33 * B.m32;
	C.m33 = A.m31 * B.m13 + A.m32 * B.m23 + A.m33 * B.m33;
	return C;
}

/**
 * @brief 3D matrix multiply 3D column vector
 * @param[in] A     Input 3D matrix
 * @param[in] B     Input 3D column vector
 * @return result 3D vector (A * B)
 */
extern v3_t m3_mul_v3(m3_t A, v3_t B) {
	v3_t C;
	C.x = A.m11 * B.x + A.m12 * B.y + A.m13 * B.z;
	C.y = A.m21 * B.x + A.m22 * B.y + A.m23 * B.z;
	C.z = A.m31 * B.x + A.m32 * B.y + A.m33 * B.z;
	return C;
}

/**
 * @biref main diagonal elements of a 3D matrix
 * @param[in] A     Input 3D matrix
 * @return 3D vector of main diagonal elements
 */
extern v3_t m3_diag(m3_t A) {
	//return (v3_t) { A.m11, A.m22, A.m33 };
	v3_t temp_t;
	temp_t.x = A.m11;
	temp_t.y = A.m22;
	temp_t.z = A.m33;
	return temp_t;
}

extern m3_t m3_diag_add(m3_t A, double a)
{
	m3_t B = A;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) B.v[i + j * 3] += a;
		}
	}
	return B;
}

/**
 * @brief Power exponent of every elements in vector
 * @param[in] A     Input vector
 * @param[in] order Power order, e.g. 2.0 means square, 0.5 means root square
 * @return Power exponent result vector
 */
extern m3_t m3_pow(m3_t A, double order) {
	//return (m3_t) {
	//    pow(A.m11, order), pow(A.m12, order), pow(A.m13, order),
	//        pow(A.m21, order), pow(A.m22, order), pow(A.m23, order),
	//        pow(A.m31, order), pow(A.m32, order), pow(A.m33, order)
	//};
	m3_t temp_t;
	temp_t.m11 = pow(A.m11, order);
	temp_t.m12 = pow(A.m12, order);
	temp_t.m13 = pow(A.m13, order);
	temp_t.m21 = pow(A.m21, order);
	temp_t.m22 = pow(A.m22, order);
	temp_t.m23 = pow(A.m23, order);
	temp_t.m31 = pow(A.m31, order);
	temp_t.m32 = pow(A.m32, order);
	temp_t.m33 = pow(A.m33, order);
	return temp_t;
}

/**
 * @brief judge the two matrix if equal or not
 * @param[in] A     First matrix
 * @param[in] B     Second matrix
 * @param[in] eps   the zero threshold e.g. 1e-20
 * @return true: A == B, false: A != B
 */
extern int m3_equal(const m3_t* A, const m3_t* B, double eps) {
	double diff;
	const double* pA = (const double*)A;
	const double* pB = (const double*)B;
	for (int i = 0; i < 9; ++i) {
		if ((diff = pA[i] - pB[i]) > eps || diff < -eps)
			return 0;
	}
	return 1;
}

/**
 * @brief Singular Value Decomposition(SVD) of 3D Matrix
 *          A = U * diag(D) * V'
 * @param[in] A   3D matrix
 * @param[out] U  3D unit othogonal matrix
 * @param[out] D  3D vector
 * @param[out] V  3D unit othogonal matrix
 * @return 0: OK
 */
extern int m3_SVD(const m3_t* A, m3_t* U, v3_t* D, m3_t* V) {
	m3_t B = m3_mul(m3_T(*A), *A);
	/* Jacobi iteration method to solve eigenvalue and eigenvector of B */
	const double thres = 1E-15;
	m3_t eVal = B, eVec = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
	double afa, sa, ca;
	m3_t Upq;
	for (int i = 0; i < 100; i++) {
		double d12 = fabs(eVal.m12), d13 = fabs(eVal.m13), d23 = fabs(eVal.m23);
		if (d12 > thres && d12 >= d13 && d12 >= d23) {
			afa = atan2(2 * eVal.m12, eVal.m11 - eVal.m22) / 2.0;
			sa = sin(afa);
			ca = cos(afa);
			//Upq = (m3_t){ ca, -sa, 0.0, sa, ca, 0.0, 0.0, 0.0, 1.0 };
			Upq.m11 = ca; Upq.m12 = -sa; Upq.m13 = 0.0;
			Upq.m21 = sa; Upq.m22 = ca; Upq.m23 = 0.0;
			Upq.m31 = 0.0; Upq.m32 = 0.0; Upq.m33 = 1.0;
		}
		else if (d13 > thres && d13 >= d23 && d13 >= d12) {
			afa = atan2(2 * eVal.m13, eVal.m33 - eVal.m11) / 2.0;
			sa = sin(afa);
			ca = cos(afa);
			//Upq = (m3_t){ ca, 0.0, sa, 0.0, 1.0, 0.0, -sa, 0.0, ca };

			Upq.m11 = ca;  Upq.m12 = 0.0; Upq.m13 = sa;
			Upq.m21 = 0.0; Upq.m22 = 1.0; Upq.m23 = 0.0;
			Upq.m31 = -sa; Upq.m32 = 0.0; Upq.m33 = ca;
		}
		else if (d23 > thres && d23 >= d12 && d23 >= d13) {
			afa = atan2(2 * eVal.m23, eVal.m22 - eVal.m33) / 2.0;
			sa = sin(afa);
			ca = cos(afa);
			//Upq = (m3_t){ 1.0, 0.0, 0.0, 0.0, ca, -sa, 0.0, sa, ca };
			Upq.m11 = 1.0; Upq.m12 = 0.0; Upq.m13 = 0.0;
			Upq.m21 = 0.0; Upq.m22 = ca;  Upq.m23 = -sa;
			Upq.m31 = 0.0; Upq.m32 = sa;  Upq.m33 = ca;
		}
		else {
			break;
		}
		eVec = m3_mul(eVec, Upq);
		eVal = m3_mul(m3_T(eVec), m3_mul(B, eVec));
	}
	/* Using eigenvalue&eigenvecotor to calculate U D V */
	*D = v3_pow(m3_diag(eVal), 0.5);
	*V = eVec;
	*U = m3_mul(*A, m3_mul(eVec, v3_diag(v3_pow(*D, -1.0))));
	return 0;
}

/**
 * @brief swap 3D matrix row
 * @param[in,out]   A   Input and swap row matrix(swap row1 and row2)
 * @param[in]       r1  row1, range: 1~3
 * @param[in]       r2  row2, range: 1~3
 */
extern void m3_swap_row(m3_t* A, int r1, int r2) {
	if (r1 == r2) return;
	double* pA = (double*)A;
	double tmp;
	int rr1 = r1 - 1, rr2 = r2 - 1;
	for (int i = 0; i < 3; ++i) {
		tmp = pA[rr1 * 3 + i];
		pA[rr1 * 3 + i] = pA[rr2 * 3 + i];
		pA[rr2 * 3 + i] = tmp;
	}
}

/**
 * @brief swap 3D matrix column
 * @param[in,out]   A   Input and swap column matrix(swap column1 and column2)
 * @param[in]       c1  column1, range: 1~3
 * @param[in]       c2 column2, range: 1~3
 */
extern void m3_swap_clm(m3_t* A, int c1, int c2) {
	if (c1 == c2) return;
	double* pA = (double*)A;
	double tmp;
	int cc1 = c1 - 1, cc2 = c2 - 1;
	for (int i = 0; i < 3; ++i) {
		tmp = pA[i * 3 + cc1];
		pA[i * 3 + cc1] = pA[i * 3 + cc2];
		pA[i * 3 + cc2] = tmp;
	}
}

/* trasform m3 diag to v3 format */
extern void m3_diag_v3(m3_t m3, v3_t* v3) {
	v3->x = m3.m11;
	v3->y = m3.m22;
	v3->z = m3.m33;
}

/**
 * @brief determinant value of a 3D matrix
 * @param[in] A Input matrix
 * @return determinant value
 */
extern double m3_det(const m3_t* A) {
	return A->m11 * A->m22 * A->m33 + A->m12 * A->m23 * A->m31
		+ A->m13 * A->m21 * A->m32 - A->m31 * A->m22 * A->m13
		- A->m21 * A->m12 * A->m33 - A->m11 * A->m32 * A->m23;
}

/**
 * @brief 3D matrix LU decomposition, express matrix A as the product of two
 *      essentially triangular matrices, and satisfying the equation:
 *      PA = LU, where L is lower triangular matrix, and U is Upper trianglur
 *      matrix, P is the row permutation matrix.
 * @param[in]   A   Input matrix
 * @param[out]  L   Lower triangular matrix
 * @param[out]  U   Upper triangular matrix
 * @param[out]  P   Row permutation matrix
 * @return 0: OK
 * @note  LU decompostion by Elimination with Maximal Column Pivoting
 */
extern int m3_LU(const m3_t* A, m3_t* L, m3_t* U, m3_t* P) {
	*L = I33;
	*U = *A;
	*P = I33;
	double* pU = (double*)U, * pL = (double*)L;
	for (int i = 0; i < 2; ++i) {
		/* choose the abs max column element as pivot */
		if (i == 0) {
			if (fabs(U->m11) >= fabs(U->m21)) {
				if (fabs(U->m11) < fabs(U->m31)) {
					m3_swap_row(U, 1, 3);
					m3_swap_row(P, 1, 3);
				}
			}
			else {
				if (fabs(U->m21) >= fabs(U->m31)) {
					m3_swap_row(U, 1, 2);
					m3_swap_row(P, 1, 2);
				}
				else {
					m3_swap_row(U, 1, 3);
					m3_swap_row(P, 1, 3);
				}
			}
		}
		else if (i == 1) {
			if (fabs(U->m22) < fabs(U->m32)) {
				m3_swap_row(U, 2, 3);
				m3_swap_row(P, 2, 3);
				/* Should exchange the m31 and m21 of L */
				double tmp = L->m31;
				L->m31 = L->m21;
				L->m21 = tmp;
			}
		}
		/* solve pivot element equal to zero */
		if (fabs(pU[i * 3 + i]) < 1E-64)
			continue;
		/* Gaussian elimination */
		for (int j = i + 1; j < 3; ++j) {
			pL[j * 3 + i] = pU[j * 3 + i] / pU[i * 3 + i];
			for (int k = 0; k < 3; ++k) {
				pU[j * 3 + k] -= pL[j * 3 + i] * pU[i * 3 + k];
			}
		}
	}
	return 0;
}

/**
 * @brief Euler attitude add operation, the same as left-multiply DCM, and
 *      similar to vector add operation by superscript&subscript
 * @param[in] Eab   Euler attitude transform from a-frame to b-frame[rad]
 * @param[in] Ebc   Euler attitude transform from b-frame to c-frame[rad]
 * @return Eac, Euler attitude transofrom from a-frame to c-frame[rad]
 */
extern v3_t euler_add(v3_t Eab, v3_t Ebc, int opt) {
	m3_t Cab;
	euler2dcm(&Eab, &Cab, opt);
	m3_t Cbc;
	euler2dcm(&Ebc, &Cbc, opt);
	m3_t Cac = m3_mul(Cbc, Cab);
	v3_t Eac;
	dcm2euler(&Cac, &Eac, opt);
	return Eac;
}

/**
 * @brief Eulear attitue delete operation, get two Euler angel difference, and
 *      similar to vector delete operation by superscript&subscript
 * @param[in] Eac   Euler attitude transform from a-frame to c-frame
 * @param[in] Eab   Euler attitude transform from a-frame to b-frame
 * @return Eac, Euler attitude transform from a-frame to c-frame
 */
extern v3_t euler_del(v3_t Eac, v3_t Eab, int opt) {
	m3_t Cab;
	euler2dcm(&Eab, &Cab, opt);
	m3_t Cac;
	euler2dcm(&Eac, &Cac, opt);
	m3_t Cbc = m3_mul(Cac, m3_T(Cab));
	v3_t Ebc;
	dcm2euler(&Cbc, &Ebc, opt);
	return Ebc;
}

/**
 * @brief generate stanard normal distribution random number
 * @return radnom number
 */
extern double ins_randn(void) {
	double N = 25.0, sum_x = 0.0;
	for (int i = 0; i < N; ++i) {
		sum_x += ((double)(rand())) / RAND_MAX;
	}
	/* standard normal distrubution */
	return (sum_x - 0.5 * N) / sqrt(N / 12.0);
}

/**
 * @brief generate an 3D vector filled with normal distribution random
 *      number
 * @param[in] mean      mean of normal distribution
 * @param[in] sigma     standard error of normal distribution
 * @return 3D random vector number
 * @see ins_randn()
 */
extern v3_t v3_randn(double mean, double sigma) {
	//return (v3_t) {
	//    mean + sigma * ins_randn(), mean + sigma * ins_randn(),
	//        mean + sigma * ins_randn()
	//};
	v3_t temp_t;
	temp_t.x = mean + sigma * ins_randn();
	temp_t.y = mean + sigma * ins_randn();
	temp_t.z = mean + sigma * ins_randn();
	return temp_t;
}

/**
 * @brief copy double array elements to 3D vector
 * @param[out] v3_dest  3D vector
 * @param[in]  src      first pointer of source double array
 */
extern void v3_copy(v3_t* v3_dest, const double* src) {
	v3_dest->x = src[0];
	v3_dest->y = src[1];
	v3_dest->z = src[2];
}

/**
 * @brief copy double matrix(column-major) elements to 3D matrix
 * @param[out]  m3_dest 3D matrix
 * @param[in]   src     first pointer of source double matrix
 * @param[in]   n       row number of source double matrix
 */
extern void m3_copy(m3_t* m3_dest, const double* src, int n) {
	double* imat = (double*)m3_dest;
	for (int j = 0; j < 3; ++j)
		for (int i = 0; i < 3; ++i)
			imat[i * 3 + j] = src[j * n + i];
}

/**
 * @brief paste matrix to double array
 * @param[out]  dest    first pointer of double array
 * @param[in]   v3      3D vector
 * @see m3_copy() v3_copy() m3_paste()
 */
extern void v3_paste(double* dest, const v3_t* v3) {
	dest[0] = v3->x;
	dest[1] = v3->y;
	dest[2] = v3->z;
}

/**
 * @brief paste matrix to double array matrix
 * @param[out]  dest    double array first address
 * @param[in]   n       matrix row number(n >= 3)
 * @param[in]   m3      3D matrix
 * @see m3_copy() v3_copy() v3_paste()
 */
extern void m3_paste(double* dest, int n, const m3_t* m3) {
	const double* imat = (const double*)m3;
	for (int j = 0; j < 3; ++j)
		for (int i = 0; i < 3; ++i)
			dest[j * n + i] = imat[i * 3 + j];
}

/**
 * @brief limit angle to (-pi, pi]
 * @param[in] 	angle	input angle[rad]
 * @return limited angle[rad]
 */
extern double angle_to180(double angle) {
	if (angle > PI) {
		angle -= 2.0 * PI;
		if (angle > PI) return angle_to180(angle);
	}
	else if (angle <= -PI) {
		angle += 2.0 * PI;
		if (angle < -PI) return angle_to180(angle);
	}
	return angle;
}

/**
 * @brief limit angle to [0, 2pi)
 * @param[in] 	angle	input angle[rad]
 * @return limit angle [rad]
 */
extern double angle_to360(double angle) {
	if (angle < 0.0) {
		angle += 2.0 * PI;
		if (angle < 0.0) return angle_to360(angle);
	}
	else if (angle >= 2.0 * PI) {
		angle -= 2.0 * PI;
		if (angle > 2.0 * PI) return angle_to360(angle);
	}
	return angle;
}