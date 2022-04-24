//
// Created by lizhen on 2021/5/31.
//
#include "rtklib.h"

/* todo check */
extern int inskf_udmeas_yaw(kf_t* inskf, double yaw, double Qyaw, const insopt_t* opt) {
	//    KF_HINIT(inskf, 1);
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

	//    KF_FILTER(inskf, &dz_yaw, &Qyaw);
	return 0;
}

static double dyaw_swerve(double veb_b_x, double web_n_z, double dL) {
	if (fabs(veb_b_x) < 0.05) {
		return 0.0;
	}
	double sinphi = dL * web_n_z / veb_b_x;
	if (fabs(sinphi) > 1.0) {
		return 0.0;
	}
	return asin(sinphi);
}

static double dyaw_swerve_Q(double dyaw, double veb_b_x, double web_n_z,
	double dL, double Qv, double Qw, double QL) {
	if (fabs(veb_b_x) < 0.05) {
		return 999.9;
	}
	double tmp = 1.0 / (veb_b_x * cos(dyaw));
	v3_t H = { -web_n_z * tmp, -dL * tmp, dL * web_n_z * tmp / veb_b_x };
	m3_t Q = I33;
	Q.m11 = QL;
	Q.m22 = Qw;
	Q.m33 = Qv;
	return v3_mul_rxc(H, m3_mul_v3(Q, H));
}

extern int inskf_constraint_yaw(kf_t* inskf, const v3_t* wib_b, const imup_t* imup, const insopt_t* opt) {
	v3_t we = { 0, 0, wgs84.wie };
	v3_t pos = inskf->insstate->pos, vel = inskf->insstate->vel;
	m3_t dcm = inskf->insstate->dcm;

	m3_t Qveb_n;
	m3_copy(&Qveb_n, &inskf->P[xiV(opt) + xiV(opt) * inskf->nx], inskf->nx);
	ecef2llhQ(&pos, NULL, &Qveb_n, NULL, opt->local_coord);
	ecef2llh(&pos, &vel, &dcm, opt->local_coord); /* Cbe => Cbn */

	double yaw = 0.0, Qyaw = 0.0;
	if (!vel2yaw(&vel, &yaw, &Qveb_n, &Qyaw)) {
		m3_t Cen = formCen_llh(pos.x, pos.y, opt->local_coord);
		v3_t web_n = v3_del(m3_mul_v3(dcm, *wib_b), m3_mul_v3(Cen, we));
		v3_t veb_b = m3_mul_v3(m3_T(dcm), vel);

		if (veb_b.x < -0.1) {
			return 1;
		}

		/* adust vecocity */
		m3_t Qveb_b = m3_mul(m3_T(dcm), m3_mul(Qveb_n, dcm));
		double Qv = Qveb_b.m11 < SQR(0.01) ? SQR(0.01) : Qveb_b.m11 * 2.0;

		double dL = imup->lever_arm_car.x;
		double dyaw_sw = dyaw_swerve(veb_b.x, web_n.z, dL);
		double Qyaw_sw = dyaw_swerve_Q(dyaw_sw, veb_b.x, web_n.z * 10.0, dL,
			Qv, SQR(imup->gyro_noise.z), SQR(imup->lever_arm_car_std.x));

		yaw = yaw - imup->err_angle_imu.z - dyaw_sw;
		Qyaw = Qyaw + SQR(imup->err_angle_imu_std.z) + Qyaw_sw;

		inskf_udmeas_yaw(inskf, yaw, Qyaw, opt);
	}
	return 0;
}

/* todo check */
static int imu_mean_var_dv(const kf_t* inskf, double* mean_dv, double* var_dv) {
	/* calculate varianace of accelration */
	double* acc = static_cast<double*>(malloc(inskf->nimud * sizeof(double)));
	double sum_acc = 0.0;
	for (unsigned short i = 0; i < inskf->nimud; ++i) {
		acc[i] = v3_norm(inskf->imud[i].accel);
		sum_acc += acc[i];
	}
	*mean_dv = sum_acc / inskf->nimud;

	sum_acc = 0.0;
	for (unsigned short i = 0; i < inskf->nimud; ++i) {
		sum_acc += SQR(acc[i] - *mean_dv);
	}
	*var_dv = sum_acc / (inskf->nimud - 1);
	free(acc);
	return 0;
}

extern int inskf_ZST_imu(kf_t* inskf, double mean_dv, double std_dv, const insopt_t* opt) {
	if (opt->zvopt.ws == 0) {
		trace(3, "opt.zvu_size = 0, Zero Speed Test disabled\n");
		inskf->ZST_count = 0;
		return false;
	}
	if (inskf->nimud < opt->zvopt.ws) {
		trace(3, "Do not have enough imu data\n");
		inskf->ZST_count = 0;
		return false;
	}
	if (fabs(mean_dv) < 1e-32 || std_dv < 1e-32) {
		inskf->ZST_count = 0;
		return false;
	}

	double mean_dv_cur, var_dv_cur;
	imu_mean_var_dv(inskf, &mean_dv_cur, &var_dv_cur);

	/* F-test */
	double Ftest = var_dv_cur / SQR(std_dv);
	/* t-test */
	if (Ftest < 20.0) {
		double t = (mean_dv_cur - mean_dv) / sqrt(var_dv_cur / (opt->zvopt.ws - 1));
		if (fabs(t) < 5.0) {
			inskf->ZST_count++;
			if (inskf->ZST_count > 5) return true; else return false;
		}
	}

	inskf->ZST_count = 0;
	return false;
}

/* todo check */
extern int inskf_ZST_pos(kf_t* inskf, const solins_t* last_sol) {
	double dt = timediff(inskf->insstate->time, last_sol->time);
	double meanv = v3_norm(v3_del(inskf->insstate->pos, last_sol->pos)) / dt;

	v3_t fac1 = v3_scalar(1.0 / v3_norm(inskf->insstate->pos), inskf->insstate->pos);
	v3_t fac2 = v3_scalar(1.0 / v3_norm(last_sol->pos), last_sol->pos);
	double Qpos1 = v3_mul_rxc(fac1, m3_mul_v3(inskf->insstate->Qpos, fac1));
	double Qpos2 = v3_mul_rxc(fac2, m3_mul_v3(last_sol->Qpos, fac2));
	double meanv_std = sqrt(Qpos1 + Qpos2) / dt;

	trace(4, "meanv %f, meanv_std %f\n", meanv, meanv_std);
	if (meanv < 0.04 && meanv_std < 0.5)
		return 1;
	else if (meanv > 0.04 && meanv > 3 * meanv_std)
		return 0;
	else
		return -1;
}