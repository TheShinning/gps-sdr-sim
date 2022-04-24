/*-----------------------------------------------------------------------------
* ins-align.cc : ins init align
*
* history : Created by chenchao on 2021/3/29.
*----------------------------------------------------------------------------*/

#include "rtklib.h"

extern int dblvec2att(const v3_t* vn1, const v3_t* vn2, const v3_t* vb1,
	const v3_t* vb2, m3_t* Cnb) {
	/**
	 * @brief Determinate attitude by using double vector
	 * @param[in] vn1   First vector under n-frame
	 * @param[in] vn2   Second vector under n-frame
	 * @param[in] vb1   First vector under b-frame
	 * @param[in] vb2   Second vector under b-frame
	 * @param[out] Cnb  Attitude transform from n-frame to b-frame(or Attitude
	 *      b-frame to n-frame)
	 * @return 0: OK
	 * @note Ref: Yan Gongming, 捷联惯导算法与组合导航讲义, 2016, P143:6.1-7
	 */
	v3_t vx = v3_scalar(1.0 / v3_norm(*vb1), *vb1);
	v3_t vtmp = v3_cross(*vb1, *vb2);
	v3_t vy = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
	vtmp = v3_cross(vtmp, *vb1);
	v3_t vz = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
	m3_t m3_b = { vx.x, vy.x, vz.x, vx.y, vy.y, vz.y, vx.z, vy.z, vz.z };

	vx = v3_scalar(1.0 / v3_norm(*vn1), *vn1);
	vtmp = v3_cross(*vn1, *vn2);
	vy = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
	vtmp = v3_cross(vtmp, *vn1);
	vz = v3_scalar(1.0 / v3_norm(vtmp), vtmp);
	m3_t m3_n = { vx.x, vx.y, vx.z, vy.x, vy.y, vy.z, vz.x, vz.y, vz.z };

	*Cnb = m3_mul(m3_b, m3_n);
	return 0;
}

static int ins_align_manual(gtime_t imu_tag, imup_t* imup, v3_t* r, v3_t* v, v3_t* rpy, gtime_t t, const insopt_t* opt, kf_t* ins_kf) {
	if (fabs(timediff(imu_tag, t)) < 0.5 / imup->freq_imu) {
		imup->initr = *r;
		imup->initv = *v;
		imup->inita = *rpy;

		char time_str[64];
		time2str(t, time_str, 3);

		if (opt->mech_coord == INSMECH_LLH) {
			v3_t re = imup->initr, ve = imup->initv;
			llh2ecef(&re, &ve, nullptr, opt->local_coord);
			fprintf(stdout, "%s: ins align ok\n", time_str);
			fprintf(stdout, "ins init position(e): %15.4f %15.4f %15.4f\n", re.x, re.y, re.z);
			fprintf(stdout, "ins init velocity(e): %15.4f %15.4f %15.4f\n", ve.x, ve.y, ve.z);
			fprintf(stdout, "ins init attitude(e): %15.4f %15.4f %15.4f\n", imup->inita.x * R2D, imup->inita.y * R2D, imup->inita.z * R2D);
			fprintf(stdout, "ins init position(n): %15.4f %15.4f %15.4f\n", imup->initr.x, imup->initr.y, imup->initr.z);
			fprintf(stdout, "ins init velocity(n): %15.4f %15.4f %15.4f\n", imup->initv.x, imup->initv.y, imup->initv.z);
			fflush(stdout);
		}
		else if (opt->mech_coord == INSMECH_ECEF) {
			v3_t rn = imup->initr, vn = imup->initv;
			//            ecef2llh(&rn, &vn, nullptr,opt->local_coord);
			fprintf(stdout, "%s: ins align ok\n", time_str);
			fprintf(stdout, "ins init position(e): %15.4f %15.4f %15.4f\n", imup->initr.x, imup->initr.y, imup->initr.z);
			fprintf(stdout, "ins init velocity(e): %15.4f %15.4f %15.4f\n", imup->initv.x, imup->initv.y, imup->initv.z);
			fprintf(stdout, "ins init attitude(e): %15.4f %15.4f %15.4f\n", imup->inita.x * R2D, imup->inita.y * R2D, imup->inita.z * R2D);
			fprintf(stdout, "ins init position(n): %15.4f %15.4f %15.4f\n", rn.x, rn.y, rn.z);
			fprintf(stdout, "ins init velocity(n): %15.4f %15.4f %15.4f\n", vn.x, vn.y, vn.z);
			fflush(stdout);
		}
		return 1;
	}
	else {
		return 0;
	}
}

extern int align_coarse_static_base(const imu_t* imu, double lat, m3_t* Cnb) {
	/**
 * @brief Coarse align on the static base
 * @param[in]   imu     static imu data
 * @param[in]   lat     imu latitude [rad]
 * @param[out]  Cnb     Output DCM attitude(a-axis refer to n-axis) on average
 * @return 0: Ok
 * @see align_coarse_inertial()
 * @see dblvec2att()
 * @note ref: Yan Gongming, 捷联惯导算法与组合导航讲义, 2016, P144
 */
 /* Fetch mean mean angular rate and specific force */
	v3_t mean_wib_b = { 0 }, mean_fib_b = { 0 };
	/* Start from the second(first epoch should not add in)*/
	for (unsigned int i = 1; i < imu->n; ++i) {
		mean_wib_b = v3_add(mean_wib_b, imu->data[i].gyro);
		mean_fib_b = v3_add(mean_fib_b, imu->data[i].accel);
	}
	double T = timediff(imu->data[imu->n - 1].time, imu->data[0].time);
	mean_wib_b = v3_scalar(1.0 / T, mean_wib_b);
	mean_fib_b = v3_scalar(1.0 / T, mean_fib_b);

	/* Double Vector to Attitude */
	v3_t gn;
	gravity_ned(lat, 0.0, &gn);
	v3_t wie_n = { wgs84.wie * cos(lat), 0.0, -wgs84.wie * sin(lat) };
	v3_t gib_b = v3_scalar(-1.0, mean_fib_b);

	dblvec2att(&gn, &wie_n, &gib_b, &mean_wib_b, Cnb);
	return 0;
}

extern int align_coarse_inertial(const imu_t* imu, double lat, m3_t* Cnb) {
	/**
	 * @brief Coarse align under inertial fame(anti-vibration method)
	 * @param[in]   imu     (quasi-)static imu data(recommend imu->n/8 = 0)
	 * @param[in]   lat     imu latitude [rad]
	 * @param[out]  Cnb     Output DCM attitude(b-axis refer to n-axis) at last moment
	 * @return 0: OK
	 * @see align_coarse_static_base()
	 * @see dblvec2att()
	 * @note ref: Qin Yunyuan, 惯性导航(2nd Edition), P361
	 */
	 /* I-frame: inertial frame of n-frame at start moment */
	 /* B-frame: inertial frame of b-frame at start moment */

	unsigned int sample_N = 4;
	double ts = timediff(imu->data[1].time, imu->data[0].time);
	/* double nts = sample_N * ts; */

	double sin_lat = sin(lat), cos_lat = cos(lat);
	v3_t gn;
	gravity_ned(lat, 0.0, &gn);

	/* Calculate vib_B1,vib_B2*/
	v3_t dtheta[4], dv[4], sum_dtheta, sum_dv;
	v3_t vib_B = { 0.0, 0.0, 0.0 }, vib_B1 = { 0.0, 0.0, 0.0 };
	quat_t qb_B = { 1.0, 0.0, 0.0, 0.0 };   /* initial attitde*/
	quat_t qk_k1; /* trans from k to k+1 */
	unsigned int ind_mid = (imu->n / sample_N) / 2 * sample_N;
	for (unsigned int i = 0; i <= imu->n - sample_N; i += sample_N) {
		for (unsigned int j = 0; j < sample_N; ++j) {
			dtheta[j] = imu->data[i + j].gyro;
			dv[j] = imu->data[i + j].accel;
		}
		/* Calculate current fib_B */
		multisample(dtheta, dv, (int)sample_N, &sum_dtheta, &sum_dv);
		vib_B = v3_add(vib_B, quat_mul_v3(qb_B, sum_dv));
		/* qb_B attitude update uner inertial frame */
		rv2quat(&sum_dtheta, &qk_k1);
		qb_B = quat_mul(qb_B, qk_k1);

		/* record middle vib_B */
		if (i == ind_mid - sample_N) vib_B1 = vib_B;
	}
	/* Calculate vib_I1, vib_I2 */
	double total_t = ts * ind_mid * 2;
	double wie_dtheta = wgs84.wie * total_t;
	double gcl_wie = gn.z * cos_lat / wgs84.wie;
	v3_t vib_I1 = { gcl_wie * sin(wie_dtheta / 2.0),
				   gcl_wie * (1 - cos(wie_dtheta / 2.0)), total_t / 2.0 * gn.z * sin_lat };
	v3_t vib_I2 = { gcl_wie * sin(wie_dtheta),
				   gcl_wie * (1 - cos(wie_dtheta)), total_t * gn.z * sin_lat };

	/* double vector to attitude */
	m3_t CB_I;
	dblvec2att(&vib_B1, &vib_B, &vib_I1, &vib_I2, &CB_I);

	/* Calculate Cnb */
	double cos_wie = cos(wie_dtheta), sin_wie = sin(wie_dtheta);
	m3_t CI_n = { -sin_lat * cos_wie, -sin_lat * sin_wie, cos_lat,
				 -sin_wie, cos_wie, 0.0,
				 -cos_lat * cos_wie, -cos_lat * sin_wie, -sin_lat };
	m3_t Cb_B;
	quat2dcm(&qb_B, &Cb_B, IMUCOOR_FRD);
	*Cnb = m3_T(m3_mul(m3_mul(CI_n, CB_I), Cb_B));

	return 0;
}

extern int align_coarse_wuhba(const imu_t* imu, double lat, const v3_t* veb_n,
	unsigned int Nveb_n, m3_t* Cnb) {
	/**
	 * @brief Coarse Alignment by solving Wuhba problem under inertial frame
	 * @param[in]   imu     IMU data
	 * @param[in]   lat     Imu latitude(Average Latitude) [rad]
	 * @param[in]   veb_n   Imu velocity uner n-frame, length: Nveb_n [m/s]
	 * @param[in]   Nveb_n  Imu velocity number
	 * @param[out]  Cnb     Output DCM attitude at last moment
	 * @return 0: Ok
	 * @warning     1. (imu->n - 1)/(Nveb_n - 1) should be an interger
	 *              2. Nveb_n >= 4(shouldn't be too small)
	 *              3. imu data should be uniform sampling
	 * @see align_coarse_static_base()
	 * @see align_coarse_inertial()
	 * @note Ref: \n
	 *  1. Peter M.G. Silson, Coarse Align of Ship's Strapdown Inertial Attitude
	 *  Reference System Using Velocity Loci, 2011 \n
	 *  2. F. Landis Markley, Attitude Determination using Vector Observations and the
	 *  Singular Value Decompostion, 1988
	 */
	 /* N-frame: inertial frame of n-frame at start moment */
	 /* B-frame: inertial frame of b-frame at start moment */
	double ts = timediff(imu->data[1].time, imu->data[0].time);
	unsigned int len_dv = (imu->n - 1) / (unsigned int)(Nveb_n - 1);
	double nts = len_dv * ts;

	double sin_lat = sin(lat), cos_lat = cos(lat);
	v3_t gn;
	gravity_ned(lat, 0.0, &gn);
	v3_t wie_n = { wgs84.wie * cos_lat, 0.0, -wgs84.wie * sin_lat };
	v3_t dtheta_ie_n = v3_scalar(nts, wie_n);

	/* Calculate vib_B1,vib_B2 */
	v3_t vib_B = { 0.0, 0.0, 0.0 };
	m3_t CbB = I33, CnN = I33;   /* initial attitde at start moment*/
	m3_t Ck_k1; /* trans from k to k+1 */
	/* int ind_mid = (imu->n) / 2; */

	v3_t* dv_N = (v3_t*)malloc(sizeof(v3_t) * (Nveb_n - 1));
	v3_t* dv_B = (v3_t*)malloc(sizeof(v3_t) * (Nveb_n - 1));

	/* TN, TN_last: wie_N x veb_N - gN */
	v3_t veb_N, veb_N_last = V0, TN, TN_last = V0;
	v3_t mean_v, wen_n, dtheta_Nn_n;

	for (unsigned int i = 1, n = 0; i < imu->n; ++i) {
		/* i start with 1: becase imu->data[0](first obs) should NOT be counted */
		/* fib_B integration, ref: Peter, 2011, eq.7 */
		vib_B = v3_add(vib_B, m3_mul_v3(CbB, imu->data[i].accel));
		/* CbB attitude update uner inertial frame */
		rv2dcm(&imu->data[i].gyro, &Ck_k1);
		CbB = m3_mul(CbB, Ck_k1);

		if (i % len_dv == 0) {  /* Save vib_B and vib_N, then reset */
			dv_B[n] = vib_B;
			vib_B = /*(v3_t)*/{ 0.0, 0.0, 0.0 };

			if (n == 0) {
				veb_N_last = veb_n[0];
				TN_last = v3_del(v3_cross(wie_n, veb_n[0]), gn);
			}
			mean_v = v3_scalar(0.5, v3_add(veb_n[n], veb_n[n + 1]));
			wen_n = /*(v3_t)*/{ mean_v.y / wgs84.R0, -mean_v.x / wgs84.R0,
							-mean_v.y * tan(lat) / wgs84.R0 };
			/* update CnN */
			dtheta_Nn_n = v3_add(dtheta_ie_n, v3_scalar(nts, wen_n));
			rv2dcm(&dtheta_Nn_n, &Ck_k1);
			CnN = m3_mul(CnN, Ck_k1);
			/* Calculate  dv_N  and save */
			veb_N = m3_mul_v3(CnN, veb_n[n + 1]);
			TN = v3_del(v3_cross(m3_mul_v3(CnN, wie_n), m3_mul_v3(CnN, veb_n[n + 1])),
				m3_mul_v3(CnN, gn));
			dv_N[n] = v3_add(v3_del(veb_N, veb_N_last),
				v3_scalar(0.5 * nts, v3_add(TN, TN_last)));

			TN_last = TN;
			veb_N_last = veb_N;
			n++;
		}
	}

	/* interleave accumulate */
	len_dv = Nveb_n / 2;
	unsigned int N_dvsum = Nveb_n - len_dv;  /* Number of vectors to determiate attitude */
	v3_t* dv_N_sum = (v3_t*)malloc(sizeof(v3_t) * N_dvsum);
	v3_t* dv_B_sum = (v3_t*)malloc(sizeof(v3_t) * N_dvsum);
	for (unsigned int i = 0; i < N_dvsum; ++i) {
		//dv_N_sum[i] = /*(v3_t)*/{ 0.0, 0.0, 0.0 };
		//dv_B_sum[i] = /*(v3_t)*/{ 0.0, 0.0, 0.0 };
		dv_N_sum[i].x = 0.0; dv_N_sum[i].y = 0.0; dv_N_sum[i].z = 0.0;
		dv_B_sum[i].x = 0.0; dv_B_sum[i].y = 0.0; dv_B_sum[i].z = 0.0;

		for (unsigned int j = 0; j < len_dv; ++j) {
			dv_N_sum[i] = v3_add(dv_N_sum[i], dv_N[i + j]);
			dv_B_sum[i] = v3_add(dv_B_sum[i], dv_B[i + j]);
		}
		/* Normalize to unit vector */
		v3_normalize(&dv_N_sum[i]);
		v3_normalize(&dv_B_sum[i]);
	}

	/* Solve wuhba problem by using SVD solution */
	m3_t B = O33;
	for (unsigned int i = 0; i < N_dvsum; ++i) {
		B = m3_add(B, v3_mul_cxr(dv_N_sum[i], dv_B_sum[i]));
	}
	m3_t U, V;
	v3_t D, temp_t = { 0 };
	m3_SVD(&B, &U, &D, &V);
	double d = m3_det(&U) * m3_det(&V); /* d = +-1 */
	temp_t.x = 1.0; temp_t.y = 1.0; temp_t.z = d;
	m3_t CBN_opt = m3_mul(U, m3_mul(v3_diag(temp_t), m3_T(V)));

	/* Calculate variance-covariance */
	/* Beasue of sensors' bias, difference between CBN_opt and CBN_true is a
	 * systemic bias, could not evaluating by RMSE */
	 /* if(Q_Enb != NULL){
		 v3_t diff_dv = {0.0}, SS_diff_dv = {0.0};
		 m3_t AN = {0.0};
		 for(int i = 0; i < N_dvsum; ++i){
			 diff_dv = v3_del(dv_N_sum[i],m3_mul_v3(CBN_opt,dv_B_sum[i]));
			 SS_diff_dv = v3_add(SS_diff_dv, v3_pow(diff_dv,2.0));
			 AN = m3_add(AN,v3_mul_cxr(dv_N_sum[i],dv_N_sum[i]));
		 }
		 v3_t var_dv_v3 = v3_scalar(1.0/(Nveb_n - len_dv),SS_diff_dv);
		 double var_dv = SQR(sqrt(var_dv_v3.i) + sqrt(var_dv_v3.j) + sqrt(var_dv_v3.k)) / 9.0;
		 AN = m3_scalar(1.0/N_dvsum, AN);
		 AN.m11 = 1 - AN.m11; AN.m22 = 1 - AN.m22; AN.m33 = 1 - AN.m33;
		 m3_inv(&AN);
		 *Q_Enb = m3_scalar(var_dv/N_dvsum, AN);
		 *Q_Enb = m3_transpose(m3_mul(m3_mul(m3_transpose(CnN),*Q_Enb),CnN));
	 }*/

	 /* Caludate Cnb at last moment: Cbn = CNn * CBN * CbB */
	* Cnb = m3_T(m3_mul(m3_mul(m3_T(CnN), CBN_opt), CbB));

	free(dv_N);
	free(dv_B);
	free(dv_N_sum);
	free(dv_B_sum);

	return 0;
}

static int ins_align_gnss_pv(imup_t* imup, const imud_t* imu_data, int* pva_idx, const pva_t* pva, const insopt_t iopt, kf_t* ins_kf) {
	int align_ok = 0, i, nsol = 5;
	static int idxs[5] = { 0 };

	if (!matchgnsssol(imu_data->time, pva_idx, pva, 0.5 / imup->freq_imu, 0)) {
		return align_ok;
	}

	if (!checkpva(pva, *pva_idx)) {
		return align_ok;
	}

	for (i = 0; i < nsol - 1; i++) idxs[i] = idxs[i + 1];
	idxs[i] = *pva_idx;
	for (i = 0; i < nsol; i++) {
		if (pva->status[idxs[i]] == SOLQ_NONE) {
			return 0;
		}
	}
	if (idxs[0] == 0) return 0;

	v3_t re = pva->pos[idxs[nsol - 1]];
	v3_t ve;
#if 0
	if (!v3_norm(pva->vel[idxs[nsol - 1]])) {
		pos2vel(&pva->pos[idxs[nsol - 1]], &pva->pos[idxs[nsol - 2]], &ve, timediff(pva->time[idxs[nsol - 1]], pva->time[idxs[nsol - 2]]));
	}
	else {
		ve = pva->vel[idxs[nsol - 1]];
	}
#endif
	pos2vel(&pva->pos[idxs[nsol - 1]], &pva->pos[idxs[nsol - 2]], &ve,
		timediff(pva->time[idxs[nsol - 1]], pva->time[idxs[nsol - 2]]));

	if (v3_norm(ve) < 5.0 || v3_norm(imu_data->gyro) > 30.0 * D2R) {
		return 0;
	}

	if (gnsspva2ins(re, ve, imu_data, imup, iopt.imu_coord, iopt.local_coord)) {
		char time_str[64];
		time2str(imu_data->time, time_str, 3);
		imup->init_tag = imu_data->time;
		v3_t rn = imup->initr, vn = imup->initv;
		ecef2llh(&rn, &vn, nullptr, iopt.local_coord);
		fprintf(stdout, "%s: ins align ok\n", time_str);
		fprintf(stdout, "ins inint position(e): %10.3f %10.3f %10.3f\n", imup->initr.x, imup->initr.y, imup->initr.z);
		fprintf(stdout, "ins inint velocity(e): %10.3f %10.3f %10.3f\n", imup->initv.x, imup->initv.y, imup->initv.z);
		fprintf(stdout, "ins inint attitude(e): %10.3f %10.3f %10.3f\n", imup->inita.x * R2D, imup->inita.y * R2D,
			imup->inita.z * R2D);
		fprintf(stdout, "ins inint position(n): %10.3f %10.3f %10.3f\n", rn.x, rn.y, rn.z);
		fprintf(stdout, "ins inint velocity(n): %10.3f %10.3f %10.3f\n", vn.x, vn.y, vn.z);
		fflush(stdout);
		if (iopt.mech_coord == INSMECH_LLH) {
			imup->initr = rn;
			imup->initv = vn;
		}
		return 1;
	}
	else {
		return 0;
	}
}

static int ins_align_gnss_obs(const prcopt_t* popt, imup_t* imup, const imu_t* imu, const int kimu_idx, int* rover_idx, int* base_idx,
	const obs_t* obs, const nav_t* nav, rtk_t* rtk) {
	int min_sol = 5, i;
	static sol_t sols[5] = { 0 };
	static int first = 1;
	unsigned char stat = SOLQ_FIX;
	static prcopt_t opt = *popt;

	if (first) {
		stat = ins_use_rtk_init(rtk, &opt, nullptr);
		first = 0;
	}

	if (!ins_use_rtk(&opt, imup, imu, kimu_idx, rover_idx, base_idx, obs, nav, rtk)) {
		return 0;
	}
	else {
		for (i = 0; i < min_sol - 1; i++) sols[i] = sols[i + 1];
		sols[i] = rtk->sol;
		for (i = 0; i < min_sol; i++) {
			//if (/*sols[i].stat > stat || */sols[i].stat == SOLQ_NONE) return 0;
			if (sols[i].stat != SOLQ_FIX) return 0;
		}
		for (i = 0; i < min_sol - 1; i++) {
			if (timediff(sols[i + 1].time, sols[i].time) > 3.0) return 0;
		}

		//v3_t re = (v3_t){ sols[min_sol - 1].rr[0], sols[min_sol - 1].rr[1], sols[min_sol - 1].rr[2] };
		//v3_t re2 = (v3_t){ sols[min_sol - 2].rr[0], sols[min_sol - 2].rr[1], sols[min_sol - 2].rr[2] };
		v3_t ve;

		v3_t re, re2;
		re.x = sols[min_sol - 1].rr[0]; re2.x = sols[min_sol - 2].rr[0];
		re.y = sols[min_sol - 1].rr[1]; re2.y = sols[min_sol - 2].rr[1];
		re.z = sols[min_sol - 1].rr[2]; re2.z = sols[min_sol - 2].rr[2];

		double dt = timediff(sols[min_sol - 1].time, sols[min_sol - 2].time);
		if (fabs(dt) < 1E-03) return 0;
		pos2vel(&re, &re2, &ve, dt);

		if (v3_norm(ve) < 5.0 || v3_norm(imu->data[kimu_idx].gyro) > 30.0 * D2R) {
			return 0;
		}

		if (gnsspva2ins(re, ve, &imu->data[kimu_idx], imup, opt.insopt.imu_coord, opt.insopt.local_coord)) {
			char time_str[64];
			time2str(imu->data[kimu_idx].time, time_str, 3);
			imup->init_tag = imu->data[kimu_idx].time;
			double sow = time2gpst(imu->data[kimu_idx].time, nullptr);
			v3_t rn = imup->initr, vn = imup->initv;
			ecef2llh(&rn, &vn, nullptr, opt.insopt.local_coord);
			fprintf(stdout, "%s(%9.2f): ins align ok\n", time_str, sow);
			fprintf(stdout, "ins inint position(e): %10.3f %10.3f %10.3f\n", imup->initr.x, imup->initr.y,
				imup->initr.z);
			fprintf(stdout, "ins inint velocity(e): %10.3f %10.3f %10.3f\n", imup->initv.x, imup->initv.y,
				imup->initv.z);
			fprintf(stdout, "ins inint attitude(e): %10.3f %10.3f %10.3f\n", imup->inita.x * R2D, imup->inita.y * R2D,
				imup->inita.z * R2D);
			fprintf(stdout, "ins inint position(n): %10.3f %10.3f %10.3f\n", rn.x * R2D, rn.y * R2D, rn.z);
			fprintf(stdout, "ins inint velocity(n): %10.3f %10.3f %10.3f\n", vn.x, vn.y, vn.z);
			fflush(stdout);

			if (opt.insopt.mech_coord == INSMECH_LLH) {
				imup->initr = rn;
				imup->initv = vn;
			}
			rtkfree(rtk);
			return 1;
		}
		else {
			return 0;
		}
	}
}

extern int ins_align(const prcopt_t* popt, const obs_t* obs, const nav_t* nav, const pva_t* pva, imu_t* imus, int* imu_idx,
	int* rover_idx, int* base_idx, int* pva_idx, rtk_t* rtk) {
	int init_ok = 0;
	switch (popt->insopt.imu_align) {
	case INS_ALIGN_MANUAL:
		init_ok = ins_align_manual(imus->data[*imu_idx].time, imus->property, &imus->property->initr, &imus->property->initv, &imus->property->inita, imus->property->init_tag, &popt->insopt, rtk->ins_kf);
		break;
	case INS_ALIGN_CORSE:
		init_ok = align_coarse_static_base(imus, 0.0, NULL);
		break;
	case INS_ALIGN_GNSS_PV:
		init_ok = ins_align_gnss_pv(imus->property, &imus->data[*imu_idx], pva_idx, pva, popt->insopt, rtk->ins_kf);
		break;
	case INS_ALIGN_GNSS_SPP:
	case INS_ALIGN_GNSS_DGPS:
	case INS_ALIGN_GNSS_PPK:
	case INS_ALIGN_GNSS_PPP:
		init_ok = ins_align_gnss_obs(popt, imus->property, imus, *imu_idx, rover_idx, base_idx, obs, nav, rtk);
		break;
	}
	if (init_ok) imus->property->tstart = imus->data[*imu_idx].time;

	return init_ok;
}