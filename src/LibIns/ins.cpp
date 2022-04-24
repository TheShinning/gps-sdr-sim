//
// Created by cc on 11/25/20.
//

#include "rtklib.h"

/**
 * @brief Euler angle attitude add a samll error angle
 * @param[in,out]   E       Euler angle attitude[rad]
 * @param[in]       phi     Small error angle at three axis[rad]
 * @return (I + phi x)E
 * @note phi must keep small, otherwise this function would lead to large
 *  calcuate error
 */
extern int euler_addphi(v3_t* E, const v3_t* phi, int opt) {
	if (v3_norm(*phi) > 1.0 * D2R) {
	}
	m3_t dcm, phix;
	euler2dcm(E, &dcm, opt);
	asymmetric_mat(phi, &phix);
	phix.m11 += 1.0;
	phix.m22 += 1.0;
	phix.m33 += 1.0;
	dcm = m3_mul(phix, dcm);
	dcm2euler(&dcm, E, opt);
	return 0;
}

/**
 * @brief Euler angle attitude add a samll error angle(minus)
 * @param[in,out]   E       Euler angle attitude[rad]
 * @param[in]       phi     Small error angle at three axis[rad]
 * @return (I - phi x)E
 * @note phi must keep small, otherwise this function would lead to large
 *  calcuate error
 */
extern int euler_delphi(v3_t* E, const v3_t* phi, int opt) {
	if (v3_norm(*phi) > 1.0 * D2R) {
	}
	m3_t dcm, phix;
	euler2dcm(E, &dcm, opt);
	asymmetric_mat(phi, &phix);
	phix = m3_scalar(-1.0, phix);
	phix.m11 += 1.0;
	phix.m22 += 1.0;
	phix.m33 += 1.0;
	dcm = m3_mul(phix, dcm);
	dcm2euler(&dcm, E, opt);
	return 0;
}

/**
 * @brief yaw angel delete, get two yaw difference
 * @param[in] yaw1  Fisrt yaw angel[rad]
 * @param[in] yaw2  Second yaw angel[rad]
 * @return yaw1 - yaw2, angle differnece, range from -PI to PI[rad]
 */
extern double yaw_del(double yaw1, double yaw2) {
	return angle_to180(yaw1 - yaw2);
}

/**
 * @brief Convert Var-covariance matrix from ned frame to ecef frame
 * @param[in]      pos  Current geodetic position, BLH[rad,m]
 * @param[in,out]  Qpos pos covariance matrix, Q_ned to Q_xyz[m^2]
 * @param[in,out]  Qvel vel covariance matrix, Q_ned to Q_xyz[m^2/s^2]
 * @param[in,out]  Qatt att covariance matrix, Q_ned to Q_xyz[rad^2]
 * @return 0: OK
 * @note Qpos, Qvel and Qatt could be NULL pointer when NOT interested
 */
extern int llh2ecefQ(const v3_t* pos, m3_t* Qpos, m3_t* Qvel, m3_t* Qatt, int opt) {
	m3_t Cen = formCen_llh(pos->x, pos->y, opt);
	m3_t Cne = m3_T(Cen);
	if (Qpos != NULL) *Qpos = m3_mul(Cne, m3_mul(*Qpos, Cen));
	if (Qvel != NULL) *Qvel = m3_mul(Cne, m3_mul(*Qvel, Cen));
	if (Qatt != NULL) *Qatt = m3_mul(Cne, m3_mul(*Qatt, Cen));
	return 0;
}

/**
 * @brief Convert var-covariance matrix from ecef frame to ned frame
 * @param[in]       xyz  Current ecef position, xyz[m]
 * @param[in,out]   Qxyz pos covaraince matrix, Q_xyz to Q_ned[m^2]
 * @param[in,out]   Qvel vel covariance matrix, Q_xyz to Q_ned[m^2/s^2]
 * @param[in,out]   Qatt att covariance matrix, Q_xyz to Q_ned[rad^2]
 * @return 0: OK
 * @note Qpos, Qvel and Qatt could be NULL pointer when NOT interested
 */
extern int ecef2llhQ(const v3_t* xyz, m3_t* Qxyz, m3_t* Qvel, m3_t* Qatt, int opt) {
	v3_t pos = *xyz;
	ecef2llh(&pos, NULL, NULL, opt);
	m3_t Cen = formCen_llh(pos.x, pos.y, opt);
	m3_t Cne = m3_T(Cen);
	if (Qxyz != NULL) *Qxyz = m3_mul(Cen, m3_mul(*Qxyz, Cne));
	if (Qvel != NULL) *Qvel = m3_mul(Cen, m3_mul(*Qvel, Cne));
	if (Qatt != NULL) *Qatt = m3_mul(Cen, m3_mul(*Qatt, Cne));
	return 0;
}

extern int checkpva(const pva_t* pva, int idx) {
	if (pva->status[idx] == SOLQ_NONE) return 0;
	if (pva->time[idx].time == 0.0) return 0;
	//    if(norm(&pva->pos[idx],3)==0.0)  return 0;
	return 1;
}

/** FB UD LR */
void imu_orientation_adjust(imud_t* imud, const char* or1, const char* or2) {
	imud_t old_imud = *imud;
	const char OR[3][3] = { "FB", "UD", "LR" };
	double* a = (double*)&imud->accel;
	double* old_a = (double*)&old_imud.accel;
	double* g = (double*)&imud->gyro;
	double* old_g = (double*)&old_imud.gyro;

	for (int i = 0; i < 3; i++) {
		if (or1[i] == OR[i][0] || or1[i] == OR[i][1]) {
			for (int j = 0; j < 3; j++) {
				if (or2[j] == OR[i][0] || or2[j] == OR[i][1]) {
					if (or1[i] == or2[j] && i != j) {
						a[i] = old_a[j];
						g[i] = old_g[j];
					}
					if (or1[i] != or2[j]) {
						a[i] = -old_a[j];
						g[i] = -old_g[j];
					}
					break;
				}
			}
		}
	}
}

/**
 * @brief Add solution type to status
 * @param[in] status        old status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  added SOL_TYPE new status
 */
extern unsigned int soltype_add(unsigned int status,
	unsigned int SOL_TYPE) {
	return status | SOL_TYPE;
}

/**
 * @brief remove solution type to status
 * @param[in] status        old status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  removed SOL_TYPE new status
 */
extern unsigned int soltype_remove(unsigned int status,
	unsigned int SOL_TYPE) {
	return status & (~(SOL_TYPE));
}

/**
 * @brief judge if status contains SOL_TYPE or not
 * @param[in] status        current solution status
 * @param[in] SOL_TYPE      solution type, see enum SOL
 * @return  removed SOL_TYPE new status
 */
extern int is_soltype(unsigned int status, unsigned int SOL_TYPE) {
	return (SOL_TYPE == (status & SOL_TYPE));
}

/**
 * @brief judge current position type either geodetic Coordinates or ECEF
 * Cartesian Coordinates by position data features
 * @param[in] 	 pos  	input pos, [m or rad,m]
 * @retval true  geodetic coordinates(Lat, lon, hgt)
 * @retval false ECEF Cartesian coordiante, xyz
 * @warning this function only work with earth surface position
 */
extern int is_blh(const v3_t* pos) {
	if (fabs(pos->x) < PI + 1e-6 && fabs(pos->y) < 2 * PI + 1e-6) {
		if (fabs(pos->z) < 1e6) return 1;
		else if (fabs(v3_norm(*pos) - wgs84.R0) < 1e6) return 0;
		else {
			return 0;
		}
	}
	return 0;
}

/**
 * @brief Convert Markov random process standard error(often named unstability)
 *      to it's driven white noise(random walk coefficient)
 * @param[in] std   standard error      [SI]
 * @param[in] T     correlative time    [SI]
 * @return random walk coefficient [SI]
 * @see markov_rw2std()
 * @note  ref Yan, IMU test & data analysis, 2012, P152
 */
extern double markov_std2rw(double std, double T) { return sqrt(2.0 * SQR(std)) / T; }

/**
 * @brief Convert Markov random process driven white noise(random walk
 *  coefficient to standard error
 * @param[in] rw    random warlk coeffcient [SI]
 * @param[in] T     correlative time   [SI]
 * @return standard error
 * @see markov_std2rw()
 * @note  ref Yan, IMU test & data analysis, 2012, P152
 */
extern double markov_rw2std(double rw, double T) { return sqrt(SQR(rw * T) / 2.0); }

extern int imu_add(imu_t* imu, const imud_t* data) {
	if (imu->n == 0 && imu->nmax == 0) {
		imu->nmax = MAXIMUOBS;
		MALLOC(imu->data, imud_t, imu->nmax);
	}
	else if (imu->nmax <= imu->n) {
		imud_t* imu_data;
		imu->nmax *= 2;
		if (!(imu_data
			= (imud_t*)realloc(imu->data, sizeof(imud_t) * imu->nmax))) {
			FREE(imu->data);
			imu->n = imu->nmax = 0;
			return -1;
		}
		imu->data = imu_data;
	}
	imu->data[imu->n++] = *data;
	return 0;
}

extern int imu_init(imu_t* imu) {
	MALLOC(imu->property, imup_t, 1);

	imu->n = 0;
	imu->nmax = MAXIMUOBS;
	MALLOC(imu->data, imud_t, imu->nmax);
	return 0;
}

extern void imu_free(imu_t* imu) {
	if (imu->data != NULL) {
		FREE(imu->data);
	}
	if (imu->property != NULL) {
		FREE(imu->property);
	}
}

extern void od_init(od_t* od) {
	od->n = 0;
	od->nmax = MAXODOBS;

	MALLOC(od->time, gtime_t, od->nmax);
	MALLOC(od->dS, double, od->nmax);
}

extern void od_init1(od_t* od) {
	od->n = 1;
	od->nmax = 1;
	MALLOC(od->time, gtime_t, 1);
	MALLOC(od->dS, double, 1);
}

extern void od_add(od_t* od, const gtime_t* time, const double* dS) {
	if (od->n == 0 && od->nmax == 0)
		od_init(od);
	else if (od->n + 1 > od->nmax) {
		od->nmax *= 2;
		REALLOC(od->time, gtime_t, od->nmax);
		REALLOC(od->dS, double, od->nmax);
	}
	od->time[od->n] = *time;
	od->dS[od->n] = *dS;
	od->n++;
}

extern void od_free(od_t* od) {
	FREE(od->time);
	FREE(od->dS);
	od->n = 0;
	od->nmax = 0;
}

extern void pva_init(pva_t* pva, int n) {
	pva->is_cov = 1;
	pva->is_ext = 1;
	pva->n = 0;
	if (n != 0) {
		pva->nmax = n;
	}
	else {
		pva->nmax = MAXPVAOBS;
	}
	MALLOC(pva->time, gtime_t, pva->nmax);
	MALLOC(pva->pos, v3_t, pva->nmax);
	MALLOC(pva->vel, v3_t, pva->nmax);
	MALLOC(pva->att, v3_t, pva->nmax);
	MALLOC(pva->status, unsigned int, pva->nmax);

	if (pva->is_cov) {
		MALLOC(pva->Qpos, m3_t, pva->nmax);
		MALLOC(pva->Qvel, m3_t, pva->nmax);
		MALLOC(pva->Qatt, m3_t, pva->nmax);
	}
	else {
	}

	if (pva->is_ext) {
		MALLOC(pva->yaw2, double, pva->nmax);
		MALLOC(pva->std_yaw2, double, pva->nmax);
		MALLOC(pva->pitch2, double, pva->nmax);
		MALLOC(pva->std_pitch2, double, pva->nmax);
		MALLOC(pva->ext_status, unsigned int, pva->nmax);
	}
	else {
	}
}

extern void pva_resize(pva_t* pva, unsigned int nmax) {
	pva->nmax = nmax;
	REALLOC(pva->time, gtime_t, nmax);
	REALLOC(pva->pos, v3_t, nmax);
	REALLOC(pva->vel, v3_t, nmax);
	REALLOC(pva->att, v3_t, nmax);
	REALLOC(pva->status, unsigned int, nmax);
	if (pva->is_cov) {
		REALLOC(pva->Qpos, m3_t, nmax);
		REALLOC(pva->Qvel, m3_t, nmax);
		REALLOC(pva->Qatt, m3_t, nmax);
	}
	if (pva->is_ext) {
		REALLOC(pva->yaw2, double, nmax);
		REALLOC(pva->std_yaw2, double, nmax);
		REALLOC(pva->pitch2, double, nmax);
		REALLOC(pva->std_pitch2, double, nmax);
		REALLOC(pva->ext_status, unsigned int, nmax);
	}
}

extern void pva_free(pva_t* pva) {
	pva->n = 0;
	pva->nmax = 0;
	FREE(pva->time);
	FREE(pva->att);
	FREE(pva->att);
	FREE(pva->vel);
	FREE(pva->pos);
	FREE(pva->status);
	if (pva->is_cov) {
		FREE(pva->Qatt);
		FREE(pva->Qvel);
		FREE(pva->Qpos);
	}
	if (pva->is_ext) {
		FREE(pva->yaw2);
		FREE(pva->std_yaw2);
		FREE(pva->pitch2);
		FREE(pva->std_pitch2);
		FREE(pva->ext_status);
	}
}