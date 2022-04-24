//
// Created by chenc on 2021/3/22.
//

#include "rtklib.h"
//#include"../Eigen/Dense"
#define NF_SYS(s,opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf_sys[s]))
static void cal_Qvv(const double *R, const double *H, int n, int m, double *Qvv)
{
	double *T1 = mat(n, m), *T2 = mat(n, n), *T3 = mat(m, n), *R_ = mat(m, m);
	matcpy(R_, R, m, m);
	matcpy(Qvv, R, m, m);

	matprint(1, H, n, m, 15, 6);
	matprint(1, R, m, m, 15, 6);


	matinv(R_, m);

	matprint(0, R_, m, m, 15, 6);
	matmul("NN", n, m, m, 1.0, H, R_, 0.0, T1);
	matmul("NT", n, n, m, 1.0, T1, H, 0.0, T2);

	matprint(0, T2, n, n, 15, 6);

	if (!matinv(T2, n)) {
		matmul("TN", m, n, n, 1.0, H, T2, 0.0, T3);

		matmul("NN", m, m, n, -1.0, T3, H, 0.0, Qvv);
		matprint(0, Qvv, m, m, 15, 9);
	}

	free(T1); free(T2); free(T3); free(R_);
}

/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
static int filter_(const double *x, const double *P, const double *H,
	const double *v, const double *R, int n, int m,
	double *xp, double *Pp, int qc, res_t *res)
{
	double *F, *Q, *K, *I, *dx;
	double *P1, *P2, *R1;
	F = mat(n, m), Q = mat(m, m), K = mat(n, m), I = eye(n), dx = mat(n, 1);
	P1 = mat(n, n), P2 = mat(n, n), R1 = mat(n, m);

	int info;

	matcpy(Q, R, m, m);
	matcpy(xp, x, n, 1);
	matmul("NN", n, m, n, 1.0, P, H, 0.0, F);       /* Q=H'*P*H+R */
	matmul("TN", m, m, n, 1.0, H, F, 1.0, Q);
	//Pp=(I-K*H')*P*(I-K*H')'+K*R*K' 公式比 Pp=(I-K*H')*P 更稳健
	if (!(info = matinv(Q, m))) {
		matmul("NN", n, m, m, 1.0, F, Q, 0.0, K);   /* K=P*H*Q^-1 */
		matmul("NN", n, 1, m, 1.0, K, v, 1.0, xp);  /* xp=x+K*v */
		matmul("NT", n, n, m, -1.0, K, H, 1.0, I);  /* I=(I-K*H') */
		matmul("NN", n, n, n, 1.0, I, P, 0.0, P1);  /* P1=(I-K*H')*P */
		matmul("NT", n, n, n, 1.0, P1, I, 0.0, P2); /* P2=(I-K*H')*P*(I-K*H')' */
		matcpy(Pp, P2, n, n);
		matmul("NN", n, m, m, 1.0, K, R, 0.0, R1);  /* R1=K*R */
		matmul("NT", n, n, m, 1.0, R1, K, 1.0, Pp); /* Pp=P2+K*R*K'=(I-K*H')*P*(I-K*H')'+K*R*K' */

		if (res) {
			double *RQ = mat(m, m);
			res->post_v = mat(m, 1); res->Qvv = mat(m, m);

			matmul("TN", m, m, m, 1.0, R, Q, 0.0, RQ);
			matmul("NN", m, m, m, 1.0, RQ, R, 0.0, res->Qvv);
			matmul("NN", m, 1, m, -1.0, RQ, v, 0.0, res->post_v);

			/*观测值的单位权中误差*/
			double *Qv = mat(m, 1), sigma0 = 0.0;
			matmul("NT", m, 1, m, 1.0, Q, v, 0.0, Qv);
			matmul("NN", 1, 1, m, 1.0, v, Qv, 0.0, &sigma0);
			res->sigma0 = SQRT(sigma0 / m);

			free(Qv); free(RQ);
		}
	}
	free(F); free(Q); free(K); free(I); free(dx);
	return info;
}
static int filter_uc(const double* x, const double* P, const double* H,
	const double* v, const double* R, int n, int m,
	double* xp, double* Pp)
{
	double* F, * Q, * K, * I, * dx;
	double* P1, * P2, * R1;
	F = mat(n, m), Q = mat(m, m), K = mat(n, m), I = eye(n), dx = mat(n, 1);
	P1 = mat(n, n), P2 = mat(n, n), R1 = mat(n, m);

	int info;

	matcpy(Q, R, m, m);
	matcpy(xp, x, n, 1);
	matmul("NN", n, m, n, 1.0, P, H, 0.0, F);       /* Q=H'*P*H+R */
	matmul("TN", m, m, n, 1.0, H, F, 1.0, Q);
	//Pp=(I-K*H')*P*(I-K*H')'+K*R*K' 公式比 Pp=(I-K*H')*P 更稳健
	if (!(info = matinv(Q, m))) {
		matmul("NN", n, m, m, 1.0, F, Q, 0.0, K);   /* K=P*H*Q^-1 */
		matmul("NN", n, 1, m, 1.0, K, v, 1.0, xp);  /* xp=x+K*v */
		matmul("NT", n, n, m, -1.0, K, H, 1.0, I);  /* I=(I-K*H') */
		matmul("NN", n, n, n, 1.0, I, P, 0.0, P1);  /* P1=(I-K*H')*P */
		matmul("NT", n, n, n, 1.0, P1, I, 0.0, P2); /* P2=(I-K*H')*P*(I-K*H')' */
		matcpy(Pp, P2, n, n);
		matmul("NN", n, m, m, 1.0, K, R, 0.0, R1);  /* R1=K*R */
		matmul("NT", n, n, m, 1.0, R1, K, 1.0, Pp); /* Pp=P2+K*R*K'=(I-K*H')*P*(I-K*H')'+K*R*K' */

	}
	free(F); free(Q); free(K); free(I); free(dx);
	return info;
}

static int filter0_(rtk_t *rtk, const double *x, const double *P, const double *H,
	const double *v, const double *R, int n, int m,
	double *xp, double *Pp, int qc, res_t *res)
{
	double *F, *Q, *K, *I, *dx,*x1,*Pp1;
	double *P1, *P2, *R1;
	F = mat(n, m), Q = mat(m, m), K = mat(n, m), I = eye(n), dx = mat(n, 1);
	P1 = mat(n, n), P2 = mat(n, n), R1 = mat(n, m);
	int info;

	matcpy(Q, R, m, m);
	matcpy(xp, x, n, 1);                                
	matmul("NN", n, m, n, 1.0, P, H, 0.0, F);       /* Q=H'*P*H+R */
	matmul("TN", m, m, n, 1.0, H, F, 1.0, Q);
	//Pp=(I-K*H')*P*(I-K*H')'+K*R*K' 公式比 Pp=(I-K*H')*P 更稳健
	if (!(info = matinv(Q, m))) {
		matmul("NN", n, m, m, 1.0, F, Q, 0.0, K);   /* K=P*H*Q^-1 */
		matmul("NN", n, 1, m, 1.0, K, v, 1.0, xp);  /* xp=x+K*v */
		matmul("NT", n, n, m, -1.0, K, H, 1.0, I);  /* I=(I-K*H') */
		matmul("NN", n, n, n, 1.0, I, P, 0.0, P1);  /* P1=(I-K*H')*P */
		matmul("NT", n, n, n, 1.0, P1, I, 0.0, P2); /* P2=(I-K*H')*P*(I-K*H')' */
		matcpy(Pp, P2, n, n);
		matmul("NN", n, m, m, 1.0, K, R, 0.0, R1);  /* R1=K*R */
		matmul("NT", n, n, m, 1.0, R1, K, 1.0, Pp); /* Pp=P2+K*R*K'=(I-K*H')*P*(I-K*H')'+K*R*K' */
#if 1
		if (1)
		{
			double *Dx = mat(n, 1);
			double *Vp = mat(m, 1);
			double *Qv = mat(m, m);
			double *tmp = mat(m, n);

			// arma::mat matDx = matxp - matx;// Dx = xp - x
			for (size_t i = 0; i < n; i++)	Dx[i] = xp[i] - x[i];
			//arma::mat matVp = matv - matH * matDx;// Vp = v - H' * Dx
			matcpy(Vp, v, m, 1);
			matmul("TN", m, 1, n, -1.0, H, Dx, 1.0, Vp);
			//arma::mat matQv = matR - matH * matPp * matH0;// Qv = R - H' * Pp * H
			matcpy(Qv, R, m, m);
			matmul("TN", m, n, n, 1.0, H, Pp, 0.0, tmp);
			matmul("NN", m, m, n, -1.0, tmp, H, 1.0, Qv);

			int i, j, s, ty, fq;

			// 主要是输出残差，但是浮点解和固定解不一致的情况下，这里不应该输出
			// 暂时还未想好怎么处理
			if (1)
			{
				for (i = 0; i < rtk->ppp_glo.nOEI; i++)
				{
					s = rtk->ppp_glo.oEI[i].sat;
					ty = rtk->ppp_glo.oEI[i].ty;
					fq = rtk->ppp_glo.oEI[i].fq;

					if (OTYPE_CP_LP == ty)//相位
					{
						//rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(matQv(i, i));
						//rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = matVp[i];
						rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(Qv[i + m * i]);//验后方差
						rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = Vp[i];//验后残差
					}
					else if (OTYPE_PR_LP == ty)//伪距
					{
						//rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(matQv(i, i));
						//rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = matVp[i];
						rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(Qv[i + m * i]);
						rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = Vp[i];
					}
				}
			}
			free(Dx);
			free(Vp);
			free(Qv);
			free(tmp);
		}
#endif

		if (res) {
			double *RQ = mat(m, m);
			res->post_v = mat(m, 1); res->Qvv = mat(m, m);

			matmul("TN", m, m, m, 1.0, R, Q, 0.0, RQ);
			matmul("NN", m, m, m, 1.0, RQ, R, 0.0, res->Qvv);
			matmul("NN", m, 1, m, -1.0, RQ, v, 0.0, res->post_v);
			/*观测值的单位权中误差*/
			double *Qv = mat(m, 1), sigma0 = 0.0;
			matmul("NT", m, 1, m, 1.0, Q, v, 0.0, Qv);
			matmul("NN", 1, 1, m, 1.0, v, Qv, 0.0, &sigma0);
			res->sigma0 = SQRT(sigma0 / m);
			if (0)
			{
				int i, j, s, ty, fq;
				for (i = 0; i < rtk->ppp_glo.nOEI; i++)
				{
					s = rtk->ppp_glo.oEI[i].sat;
					ty = rtk->ppp_glo.oEI[i].ty;
					fq = rtk->ppp_glo.oEI[i].fq;

					if (OTYPE_CP_LP == ty)//相位
					{
						//rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(matQv(i, i));
						//rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = matVp[i];
						rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(res->Qvv[i + m * i])/** res->sigma0*/;//验后方差
						rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = res->post_v[i];//验后残差
					}
					else if (OTYPE_PR_LP == ty)//伪距
					{
						//rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(matQv(i, i));
						//rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = matVp[i];
						rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(res->Qvv[i + m * i]);
						rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = res->post_v[i];
					}
					i = i;
				}
			}

			free(Qv); free(RQ);
		}
	}
	free(F); free(Q); free(K); free(I); free(dx);
	free(P1); free(P2); free(R1); 
	return info;
}

/*ref to "A Variational Bayesian-Based Robust Adaptive Filtering for Precise Point Positioning Using Undifferenced and Uncombined Observations"*/
static int vbakf_(const double *x, const double *P, const double *H, const double *v,
	const double *R, int n, int m, double *xp, double *Pp)
{
	double *x_1 = mat(n, 1), *xk1k = mat(n, 1), *Pk1k1 = mat(n, n), *Pk1k = mat(n, n);
	double *Tk1k = mat(n, n), *Ak = mat(n, n), *Tkk = mat(n, n), *E_i_Pk1k = mat(n, n);
	double *P_0 = mat(m, n), *P_1 = mat(m, m), *zk1k = mat(m, 1), *E_i_Pk1k_0 = mat(m, n);
	double *Pzzk1k = mat(m, m), *Pxzk1k = mat(n, m), *Kk = mat(n, m), *Kk_ = mat(n, n);
	double *P1 = mat(n, n), *P2 = mat(n, n), *R1 = mat(n, m), *I = eye(n), *v_post_ = mat(m, m), *v_post = mat(m, 1);
	double tao_P = 3.0, v_all1 = 0.0, v_all2 = 0.0, V_all1 = 0.0, V_all2 = 0.0;
	double *F = mat(n, m), *Q = mat(m, m), *v_N = mat(m, 1), *T = mat(m, 1), *RI = mat(m, m);
	int N = 10, tk1k, tkk, info, i, j, jj;
	float fabs_v;

	matcpy(Q, R, m, m);
	matmul("NN", n, m, n, 1.0, P, H, 0.0, F);                //F=P*H       F(n,m),P(n,n),H(n,m)
	matmul("TN", m, m, n, 1.0, H, F, 1.0, Q);                //Q=H'*F+R    Q(m,m),H(n,m),F(n,m)
	if (!(info = matinv(Q, m))) {
		matmul("NN", m, m, m, 1.0, R, Q, 0.0, v_post_);
	}
	matmul("NN", m, 1, m, -1.0, v_post_, v, 0.0, v_post);    //v_postres=-R*inv(Q)*v_prires;

	/*robust,基于后验残差*/
	matcpy(RI, R, m, m);
	for (j = 0; j < m; j++) {
		fabs_v = fabs(v_post[j]);
		v_N[j] = fabs_v / sqrt(v_post_[j*m + j] * RI[j*m + j]);
		if (j % 2 == 0) v_all1 = v_all1 + v_N[j];
		if (j % 2 == 1) v_all2 = v_all2 + v_N[j];
	}
	for (jj = 0; jj < (m / 2); jj++) {/*phase*/
		V_all1 = V_all1 + SQR(v_N[2 * jj] - v_all1 / (m / 2));
	}
	for (jj = 0; jj < (m / 2); jj++) { /*pseudorange*/
		V_all2 = V_all2 + SQR(v_N[2 * jj + 1] - v_all2 / (m / 2));
	}
	for (j = 0; j < m; j++) {
		if (j % 2 == 0) { /*phase*/
			T[j] = fabs(v_N[j] - v_all1 / (m / 2)) / SQRT(V_all1 / (m / 2));
			if (T[j] > tdistb_0250[m / 2] && T[j] < tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * T[j] / tdistb_0250[m / 2] * SQR((tdistb_0005[m / 2] - tdistb_0250[m / 2]) / (tdistb_0005[m / 2] - T[j])); //down weight
			}
			if (T[j] > tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * 10000000.0; //rejected
			}
		}
		if (j % 2 == 1) { /*pseudorange*/
			T[j] = fabs(v_N[j] - v_all2 / (m / 2)) / SQRT(V_all2 / (m / 2));
			if (T[j] > tdistb_0250[m / 2] && T[j] < tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * T[j] / tdistb_0250[m / 2] * SQR((tdistb_0005[m / 2] - tdistb_0250[m / 2]) / (tdistb_0005[m / 2] - T[j]));
			}
			if (T[j] > tdistb_0005[m / 2 - 1]) {
				RI[j*m + j] = RI[j*m + j] * 100000000.0;
			}
		}
	}

	/*adaptive*/
	matcpy(Pk1k, P, n, n);
	tk1k = n + 1 + tao_P;                                                                     /*tk1k = (nx + 1 + tao_P)*/
	matmul("NN", n, n, n, 1.0, mat_scale(n, tao_P), Pk1k, 0.0, Tk1k);    /*Tk1k=tao_P*Pk1k*/
	matcpy(xp, x, n, 1);                                                             /*xkk=xk1k*/
	matcpy(xk1k, x, n, 1);
	matcpy(Pp, Pk1k, n, n);                                                              /*Pkk=Pk1k*/
	for (i = 0; i < N; i++) {
		matcpy(Ak, Pp, n, n);
		matcpy(x_1, xp, n, 1);
		matmul("NN", n, 1, n, -1.0, eye(n), xk1k, 1.0, x_1);
		matmul("NT", n, n, 1, 1.0, x_1, x_1, 1.0, Ak);                      /*Ak=(xkk-xk1k)*(xkk-xk1k)'+Pp*/
		tkk = tk1k + 1;                                                                        /*tkk=tk1k+1*/
		matcpy(Tkk, Ak, n, n);                                                                 /*Tkk=Tk1k+Ak*/
		matmul("NN", n, n, n, 1.0, eye(n), Tk1k, 1.0, Tkk);
		if (!(info = matinv(Tkk, n))) {                                                            /*E_i_Pk1k=(tkk-nx-1)*inv(Tkk)*/
			matmul("NN", n, n, n, 1.0, mat_scale(n, (tkk - n - 1)*1.0), Tkk, 0.0, E_i_Pk1k);
		}
		matcpy(Pzzk1k, RI, m, m);
		if (!(info = matinv(E_i_Pk1k, n))) {                                                    /*D_Pk1k = inv(E_i_Pk1k)*/
			matmul("TN", m, n, n, 1.0, H, E_i_Pk1k, 0.0, E_i_Pk1k_0);
			matmul("NN", m, m, n, 1.0, E_i_Pk1k_0, H, 1.0, Pzzk1k);             /*Pzzk1k = H*D_Pk1k*H'+D_R   Pzzk1k=H*Pk1k*H'+R*/
			matmul("NN", n, m, n, 1.0, E_i_Pk1k, H, 0.0, Pxzk1k);               /*Pxzk1k = D_Pk1k*H'         Pxzk1k=Pk1k*H'*/
			if (!(info = matinv(Pzzk1k, m))) {
				matmul("NN", n, m, m, 1.0, Pxzk1k, Pzzk1k, 0.0, Kk);            /*Kk=Pxzk1k*inv(Pzzk1k)      Kk=Pxzk1k*inv(Pzzk1k)*/
				matcpy(xp, xk1k, n, 1);
				matmul("NN", n, 1, m, 1.0, Kk, v, 1.0, xp);//v                /*xp=xk1k+Kk*(z-H*xk1k)      xkk=xk1k+Kk*(z-H*xk1k)*/
				matmul("NT", n, n, m, 1.0, Kk, H, 0.0, Kk_);
				matcpy(Pp, E_i_Pk1k, n, n);
				matmul("NN", n, n, n, -1.0, Kk_, E_i_Pk1k, 1.0, Pp);            /*Pp=D_Pk1k-Kk*H*D_Pk1k      Pkk=Pk1k-Kk*H*Pk1k*/
			}
		}
	}
	free(x_1);
	free(xk1k); free(Pk1k1); free(Pk1k);
	free(Tk1k); free(Ak); free(Tkk); free(E_i_Pk1k);
	free(P_0); free(P_1); free(zk1k); free(E_i_Pk1k_0);
	free(Pzzk1k); free(Pxzk1k); free(Kk); free(Kk_);
	free(P1); free(P2); free(R1); free(I);
	free(F); free(RI); free(v_N); free(T); free(Q);
	free(v_post); free(v_post_);
	return info;
}
static int vbakf0_(rtk_t *rtk,const double *x, const double *P, const double *H, const double *v,
	const double *R, int n, int m, double *xp, double *Pp, res_t *res)
{
	double *x_1 = mat(n, 1), *xk1k = mat(n, 1), *Pk1k1 = mat(n, n), *Pk1k = mat(n, n);
	double *Tk1k = mat(n, n), *Ak = mat(n, n), *Tkk = mat(n, n), *E_i_Pk1k = mat(n, n);
	double *P_0 = mat(m, n), *P_1 = mat(m, m), *zk1k = mat(m, 1), *E_i_Pk1k_0 = mat(m, n);
	double *Pzzk1k = mat(m, m), *Pxzk1k = mat(n, m), *Kk = mat(n, m), *Kk_ = mat(n, n);
	double *P1 = mat(n, n), *P2 = mat(n, n), *R1 = mat(n, m), *I = eye(n), *v_post_ = mat(m, m), *v_post = mat(m, 1);
	double tao_P = 3.0, v_all1 = 0.0, v_all2 = 0.0, V_all1 = 0.0, V_all2 = 0.0;
	double *F = mat(n, m), *Q = mat(m, m), *v_N = mat(m, 1), *T = mat(m, 1), *RI = mat(m, m);
	int N = 10, tk1k, tkk, info, i, j, jj;
	float fabs_v;

	matcpy(Q, R, m, m);
	matmul("NN", n, m, n, 1.0, P, H, 0.0, F);                //F=P*H       F(n,m),P(n,n),H(n,m)
	matmul("TN", m, m, n, 1.0, H, F, 1.0, Q);                //Q=H'*F+R    Q(m,m),H(n,m),F(n,m)
	if (!(info = matinv(Q, m))) {
		matmul("NN", m, m, m, 1.0, R, Q, 0.0, v_post_);
	}
	matmul("NN", m, 1, m, -1.0, v_post_, v, 0.0, v_post);    //v_postres=-R*inv(Q)*v_prires;

	/*robust,基于后验残差*/
	matcpy(RI, R, m, m);
	for (j = 0; j < m; j++) {
		fabs_v = fabs(v_post[j]);
		v_N[j] = fabs_v / sqrt(v_post_[j*m + j] * RI[j*m + j]);
		if (j % 2 == 0) v_all1 = v_all1 + v_N[j];
		if (j % 2 == 1) v_all2 = v_all2 + v_N[j];
	}
	for (jj = 0; jj < (m / 2); jj++) {/*phase*/
		V_all1 = V_all1 + SQR(v_N[2 * jj] - v_all1 / (m / 2));
	}
	for (jj = 0; jj < (m / 2); jj++) { /*pseudorange*/
		V_all2 = V_all2 + SQR(v_N[2 * jj + 1] - v_all2 / (m / 2));
	}
	for (j = 0; j < m; j++) {
		if (j % 2 == 0) { /*phase*/
			T[j] = fabs(v_N[j] - v_all1 / (m / 2)) / SQRT(V_all1 / (m / 2));
			if (T[j] > tdistb_0250[m / 2] && T[j] < tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * T[j] / tdistb_0250[m / 2] * SQR((tdistb_0005[m / 2] - tdistb_0250[m / 2]) / (tdistb_0005[m / 2] - T[j])); //down weight
			}
			if (T[j] > tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * 10000000.0; //rejected
			}
		}
		if (j % 2 == 1) { /*pseudorange*/
			T[j] = fabs(v_N[j] - v_all2 / (m / 2)) / SQRT(V_all2 / (m / 2));
			if (T[j] > tdistb_0250[m / 2] && T[j] < tdistb_0005[m / 2]) {
				RI[j*m + j] = RI[j*m + j] * T[j] / tdistb_0250[m / 2] * SQR((tdistb_0005[m / 2] - tdistb_0250[m / 2]) / (tdistb_0005[m / 2] - T[j]));
			}
			if (T[j] > tdistb_0005[m / 2 - 1]) {
				RI[j*m + j] = RI[j*m + j] * 100000000.0;
			}
		}
	}

	/*adaptive*/
	matcpy(Pk1k, P, n, n);
	tk1k = n + 1 + tao_P;                                                                     /*tk1k = (nx + 1 + tao_P)*/
	matmul("NN", n, n, n, 1.0, mat_scale(n, tao_P), Pk1k, 0.0, Tk1k);    /*Tk1k=tao_P*Pk1k*/
	matcpy(xp, x, n, 1);                                                             /*xkk=xk1k*/
	matcpy(xk1k, x, n, 1);
	matcpy(Pp, Pk1k, n, n);                                                              /*Pkk=Pk1k*/
	for (i = 0; i < N; i++) {
		matcpy(Ak, Pp, n, n);
		matcpy(x_1, xp, n, 1);
		matmul("NN", n, 1, n, -1.0, eye(n), xk1k, 1.0, x_1);
		matmul("NT", n, n, 1, 1.0, x_1, x_1, 1.0, Ak);                      /*Ak=(xkk-xk1k)*(xkk-xk1k)'+Pp*/
		tkk = tk1k + 1;                                                                        /*tkk=tk1k+1*/
		matcpy(Tkk, Ak, n, n);                                                                 /*Tkk=Tk1k+Ak*/
		matmul("NN", n, n, n, 1.0, eye(n), Tk1k, 1.0, Tkk);
		if (!(info = matinv(Tkk, n))) {                                                            /*E_i_Pk1k=(tkk-nx-1)*inv(Tkk)*/
			matmul("NN", n, n, n, 1.0, mat_scale(n, (tkk - n - 1)*1.0), Tkk, 0.0, E_i_Pk1k);
		}
		matcpy(Pzzk1k, RI, m, m);
		if (!(info = matinv(E_i_Pk1k, n))) {                                                    /*D_Pk1k = inv(E_i_Pk1k)*/
			matmul("TN", m, n, n, 1.0, H, E_i_Pk1k, 0.0, E_i_Pk1k_0);
			matmul("NN", m, m, n, 1.0, E_i_Pk1k_0, H, 1.0, Pzzk1k);             /*Pzzk1k = H*D_Pk1k*H'+D_R   Pzzk1k=H*Pk1k*H'+R*/
			matmul("NN", n, m, n, 1.0, E_i_Pk1k, H, 0.0, Pxzk1k);               /*Pxzk1k = D_Pk1k*H'         Pxzk1k=Pk1k*H'*/
			if (!(info = matinv(Pzzk1k, m))) {
				matmul("NN", n, m, m, 1.0, Pxzk1k, Pzzk1k, 0.0, Kk);            /*Kk=Pxzk1k*inv(Pzzk1k)      Kk=Pxzk1k*inv(Pzzk1k)*/
				matcpy(xp, xk1k, n, 1);
				matmul("NN", n, 1, m, 1.0, Kk, v, 1.0, xp);//v                /*xp=xk1k+Kk*(z-H*xk1k)      xkk=xk1k+Kk*(z-H*xk1k)*/
				matmul("NT", n, n, m, 1.0, Kk, H, 0.0, Kk_);
				matcpy(Pp, E_i_Pk1k, n, n);
				matmul("NN", n, n, n, -1.0, Kk_, E_i_Pk1k, 1.0, Pp);            /*Pp=D_Pk1k-Kk*H*D_Pk1k      Pkk=Pk1k-Kk*H*Pk1k*/
			}
		}
	}
#if 1
	if (1)
	{
		double *Dx = mat(n, 1);
		double *Vp = mat(m, 1);
		double *Qv = mat(m, m);
		double *tmp = mat(m, n);

		// arma::mat matDx = matxp - matx;// Dx = xp - x
		for (size_t i = 0; i < n; i++)	Dx[i] = xp[i] - x[i];
		//arma::mat matVp = matv - matH * matDx;// Vp = v - H' * Dx
		matcpy(Vp, v, m, 1);
		matmul("TN", m, 1, n, -1.0, H, Dx, 1.0, Vp);
		//arma::mat matQv = matR - matH * matPp * matH0;// Qv = R - H' * Pp * H
		matcpy(Qv, R, m, m);
		matmul("TN", m, n, n, 1.0, H, Pp, 0.0, tmp);
		matmul("NN", m, m, n, -1.0, tmp, H, 1.0, Qv);

		int i, j, s, ty, fq;

		// 主要是输出残差，但是浮点解和固定解不一致的情况下，这里不应该输出
		// 暂时还未想好怎么处理
		if (1)
		{
			for (i = 0; i < rtk->ppp_glo.nOEI; i++)
			{
				s = rtk->ppp_glo.oEI[i].sat;
				ty = rtk->ppp_glo.oEI[i].ty;
				fq = rtk->ppp_glo.oEI[i].fq;

				if (OTYPE_CP_LP == ty)
				{
					//rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(matQv(i, i));
					//rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = matVp[i];
					rtk->ppp_glo.equS[s - 1].eS_cp.sig_vp[fq] = SQRT(Qv[i + m * i]);
					rtk->ppp_glo.equS[s - 1].eS_cp.resi_1[fq] = Vp[i];
				}
				else if (OTYPE_PR_LP == ty)
				{
					//rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(matQv(i, i));
					//rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = matVp[i];
					rtk->ppp_glo.equS[s - 1].eS_pr.sig_vp[fq] = SQRT(Qv[i + m * i]);
					rtk->ppp_glo.equS[s - 1].eS_pr.resi_1[fq] = Vp[i];
				}
			}
		}
		free(Dx);
		free(Vp);
		free(Qv);
		free(tmp);
	}
#endif
	if (res) {
		double *RQ = mat(m, m);
		res->post_v = mat(m, 1); res->Qvv = mat(m, m);

		matmul("TN", m, m, m, 1.0, R, Q, 0.0, RQ);
		matmul("NN", m, m, m, 1.0, RQ, R, 0.0, res->Qvv);
		matmul("NN", m, 1, m, -1.0, RQ, v, 0.0, res->post_v);

		/*观测值的单位权中误差*/
		double *Qv = mat(m, 1), sigma0 = 0.0;
		matmul("NT", m, 1, m, 1.0, Q, v, 0.0, Qv);
		matmul("NN", 1, 1, m, 1.0, v, Qv, 0.0, &sigma0);
		res->sigma0 = SQRT(sigma0 / m);

		free(Qv); free(RQ);
	}

	free(x_1);
	free(xk1k); free(Pk1k1); free(Pk1k);
	free(Tk1k); free(Ak); free(Tkk); free(E_i_Pk1k);
	free(P_0); free(P_1); free(zk1k); free(E_i_Pk1k_0);
	free(Pzzk1k); free(Pxzk1k); free(Kk); free(Kk_);
	free(P1); free(P2); free(R1); free(I);
	free(F); free(RI); free(v_N); free(T); free(Q);
	free(v_post); free(v_post_);
	return info;
}

static int sage_husa_(const double *x, const double *P, const double *H, const double *v,
	const double *R, int n, int m, double *xp, double *Pp)
{
	int i, j, info;
	double *Pxykk_1, *Py0, *ykk_1, *rk, ry = 1.0, *R_;
	static double beta = 1.0;

	Pxykk_1 = mat(n, m); Py0 = mat(m, m); ykk_1 = mat(m, 1); rk = mat(m, 1); R_ = mat(m, m);
	matcpy(R_, R, m, m);
	matcpy(Pp, P, n, n);

	matprint(0, P, n, n, 15, 6);
	matprint(0, H, n, m, 15, 6);
	matprint(0, R, m, m, 15, 6);


	matmul("NN", n, m, n, 1.0, P, H, 0.0, Pxykk_1);
	matmul("TN", m, m, n, 1.0, H, Pxykk_1, 0.0, Py0);
	matmul("TN", m, 1, n, 1.0, H, x, 0.0, ykk_1);
#if 1
	matprint(0, Pxykk_1, n, m, 15, 6);
	matprint(0, Py0, m, m, 15, 6);
	matprint(0, ykk_1, m, 1, 15, 6);

#endif
	for (i = 0; i < m; i++) {
		rk[i] = v[i] - ykk_1[i];
	}

	for (i = 0; i < m; i++) {
		if (v[i] > 1E10) continue;
		ry = SQR(rk[i]) - Py0[i + i * m];
		double bbb = R_[i + i * m];
		double aaa = (1.0 - beta)*R_[i + i * m] + beta * ry;
		R_[i + i * m] = (1.0 - beta)*R_[i + i * m] + beta * ry;
	}
	beta = beta / (beta + 0.5);
	matprint(0, R_, m, m, 15, 6);

	double *Pykk_1, *Pykk_1_, *Kk;
	Pykk_1 = mat(m, m), Pykk_1_ = mat(m, m);
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			Pykk_1[i + j * m] = Py0[i + j * m] + R_[i + j * m];
		}
	}
	matcpy(Pykk_1_, Pykk_1, m, m);
	matprint(0, Pykk_1, m, m, 15, 6);

	Kk = mat(n, m);
	if (!(info = matinv(Pykk_1, m))) {
		matmul("NN", n, m, m, 1.0, Pxykk_1, Pykk_1, 0.0, Kk);
		matmul("NN", n, 1, m, 1.0, Kk, rk, 1.0, xp);
	}
	matprint(0, Kk, n, m, 15, 6);
	matprint(1, xp, n, 1, 15, 6);

	double *Pk = mat(n, m);
	matmul("NN", n, m, m, 1.0, Kk, Pykk_1_, 0.0, Pk);
	matprint(0, Pp, n, n, 15, 6);
	matprint(0, Pk, n, m, 15, 6);

	matmul("NT", n, n, m, -1.0, Pk, Kk, 1.0, Pp);
	matprint(0, Pp, n, n, 15, 6);

	//    for(i=0;i<m;i++){
	//        for(j=0;j<m;j++){
	//            Pp[i+j*m]=(Pp[i+j*m]+Pp[j+i*m])/2.0;
	//        }
	//    }

	matprint(0, Pp, n, n, 15, 6);


	free(Pxykk_1); free(Py0); free(ykk_1); free(rk); free(R_);
	free(Pykk_1); free(Pykk_1_); free(Kk); free(Pk);
	return info;
}

//extern int filter(double *x, double *P, const double *H, const double *v,
//                  const double *R, int n, int m,int qc, const int kf_type,res_t *res,int tc)
//{
//    double *x_,*xp_,*P_,*Pp_,*H_;
//    int i,j,k,info,*ix;
//
//    /* create list of non-zero states */
//    ix=imat(n,1); for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
//    x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
//    /* compress array by removing zero elements to save computation time */
//    for (i=0;i<k;i++) {
//        x_[i]=x[ix[i]];
//        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
//        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
//    }
//    if(tc){
//#if 1
//        fprintf(stdout,"coupled prior:\n");
//        matprint(0,v,1,m,20,10);
//        matprint(1,H_,k,m,20,10);
//        matprint(0,x_,1,k,20,10);
//        matprint(0,P_,k,k,20,10);
//#endif
//    }
//
//    /* do kalman filter state update on compressed arrays */
//	switch (kf_type) {
//	case KFOPT_VBKF:
//		info = vbakf_(x_, P_, H_, v, R, k, m, xp_, Pp_);
//		break;
//	case KFOPT_SAGE_HUSA:
//		info = sage_husa_(x_, P_, H_, v, R, k, m, xp_, Pp_);
//		break;
//	default:
//		info = filter_(x_, P_, H_, v, R, k, m, xp_, Pp_, qc, res);
//		break;
//	}
//
//    if(tc){
//#if 1
//        fprintf(stdout,"coupled post:\n");
//        matprint(0,R,m,m,15,6);
//        matprint(0,xp_,1,k,20,10);
//        matprint(0,Pp_,k,k,20,10);
//#endif
//    }
//
//    /* copy values from compressed arrays back to full arrays */
//    for (i=0;i<k;i++) {
//        x[ix[i]]=xp_[i];
//        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
//    }
//    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
//    return info;
//}
extern int filter(double *x, double *P, const double *H, const double *v,
	const double *R, int n, int m, int qc, const int kf_type, res_t *res, int tc)
{
	double *x_, *xp_, *P_, *Pp_, *H_;
	int i, j, k, info, *ix;

	/* create list of non-zero states */
	ix = imat(n, 1); for (i = k = 0; i < n; i++) if (x[i] != 0.0&&P[i + i * n] > 0.0) ix[k++] = i;
	x_ = mat(k, 1); xp_ = mat(k, 1); P_ = mat(k, k); Pp_ = mat(k, k); H_ = mat(k, m);
	/* compress array by removing zero elements to save computation time */
	for (i = 0; i < k; i++) {
		x_[i] = x[ix[i]];
		for (j = 0; j < k; j++) P_[i + j * k] = P[ix[i] + ix[j] * n];
		for (j = 0; j < m; j++) H_[i + j * k] = H[ix[i] + j * n];
	}
	trace(3, "befor filter\n P(%d)=\n", i + 1); tracemat(3, P_, k, k, 13, 4);
	trace(3, "H(%d)=\n", i + 1); tracemat(3, H_, k, m, 13, 4);
	trace(3, "X(%d)=\n", i + 1); tracemat(3, x_, 1, k, 13, 4);

	if (tc) {
#if 1
		fprintf(stdout, "coupled prior:\n");
		matprint(0, v, 1, m, 20, 10);
		matprint(1, H_, k, m, 20, 10);
		matprint(0, x_, 1, k, 20, 10);
		matprint(0, P_, k, k, 20, 10);
#endif
	}

	/* do kalman filter state update on compressed arrays */
	switch (kf_type) {
	case KFOPT_VBKF:
		info = vbakf_(x_, P_, H_, v, R, k, m, xp_, Pp_);
		break;
	case KFOPT_SAGE_HUSA:
		info = sage_husa_(x_, P_, H_, v, R, k, m, xp_, Pp_);
		break;
	default:
		info = filter_(x_, P_, H_, v, R, k, m, xp_, Pp_, qc, res);
		break;
	}

	if (tc) {
#if 1
		fprintf(stdout, "coupled post:\n");
		matprint(0, R, m, m, 15, 6);
		matprint(0, xp_, 1, k, 20, 10);
		matprint(0, Pp_, k, k, 20, 10);
#endif
	}
	trace(3, "Xp(%d)=\n", i + 1); tracemat(3, xp_, 1, k, 13, 4);
	/* copy values from compressed arrays back to full arrays */
	for (i = 0; i < k; i++) {
		x[ix[i]] = xp_[i];
		for (j = 0; j < k; j++) P[ix[i] + ix[j] * n] = Pp_[i + j * k];
	}
	free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
	return info;
}
//rtk->ppp_glo.nSat;

extern int filter0(rtk_t* rtk, double* x, double* P, const double* H, const double* v,
	const double* R, int n, int m, int qc, const int kf_type, res_t* res, int tc)
{
	double* x_, * xp_, * P_, * Pp_, * H_;
	int i, j, k, info, * ix;
	/* create list of non-zero states */
	ix = imat(n, 1);
	for (i = k = 0; i < n; i++)
	{
		if (x[i] != 0.0 && P[i + i * n] > 0.0)
		{
			ix[k++] = i;
		}
	}
	x_ = mat(k, 1); xp_ = mat(k, 1); P_ = mat(k, k); Pp_ = mat(k, k); H_ = mat(k, m);
	/* compress array by removing zero elements to save computation time */
	for (i = 0; i < k; i++) {
		x_[i] = x[ix[i]];
		for (j = 0; j < k; j++) P_[i + j * k] = P[ix[i] + ix[j] * n];
		for (j = 0; j < m; j++) H_[i + j * k] = H[ix[i] + j * n];
	}
	trace(3, "befor filter\n P(%d)=\n", i + 1); tracemat(3, P_, k, k, 13, 4);
	trace(3, "H(%d)=\n", i + 1); tracemat(3, H_, k, m, 13, 4);
	trace(3, "X(%d)=\n", i + 1); tracemat(3, x_, 1, k, 13, 4);
	trace(3, "V(%d)=\n", i + 1); tracemat(3, v, 1, m, 13, 4);
	trace(3, "R(%d)=\n", i + 1); tracemat(3, R, m, m, 13, 6);
#if 1
	if (tc) {
		fprintf(stdout, "coupled prior:\n");
		matprint(0, v, 1, m, 20, 10);
		matprint(1, H_, k, m, 20, 10);
		matprint(0, x_, 1, k, 20, 10);
		matprint(0, P_, k, k, 20, 10);
	}
#endif
	/* do kalman filter state update on compressed arrays */
	switch (kf_type) {
	case KFOPT_VBKF:
		//info = vbakf_(x_, P_, H_, v, R, k, m, xp_, Pp_);
		info = vbakf0_(rtk,x_, P_, H_, v, R, k, m, xp_, Pp_,res);
		break;
	case KFOPT_SAGE_HUSA:
		info = sage_husa_(x_, P_, H_, v, R, k, m, xp_, Pp_);
		break;
	default:
		info = filter0_(rtk,x_, P_, H_, v, R, k, m, xp_, Pp_, qc, res);
		break;
	}

	trace(3, "after filter\n", i + 1);
	trace(3, "X(%d)=\n", i + 1); tracemat(3, xp_, 1, k, 13, 4);
	trace(3, "Pp(%d)=\n", i + 1); tracemat(0, Pp_, k, k, 20, 10);

	if (tc) {
#if 0
		fprintf(stdout, "coupled post:\n");
		matprint(0, R, m, m, 15, 6);
		matprint(0, xp_, 1, k, 20, 10);
		matprint(0, Pp_, k, k, 20, 10);
#endif
	}

	/* copy values from compressed arrays back to full arrays */
	for (i = 0; i < k; i++) {
		x[ix[i]] = xp_[i];
		for (j = 0; j < k; j++) P[ix[i] + ix[j] * n] = Pp_[i + j * k];
	}
	free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
	return info;
}

extern void MatrixMult(const double *a, const double *b, double *c, const int m, const int n, const int k)
{
	int i, j, l;
	for (i = 0; i < m; i++) {
		for (j = 0; j < k; j++) {
			c[i * k + j] = 0.0;
			for (l = 0; l < n; l++) {
				c[i * k + j] = c[i * k + j] + a[i * n + l] * b[l * k + j];
			}
		}
	}
}
/* purpose : matrix translate.														*/
/* 矩阵转置a[m][n]*/
extern void MatrixTranslate(const double *a, double *b, const int m, const int n)
{
	int i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			b[j * m + i] = a[i * n + j];
	}
}


