/*------------------------------------------------------------------------------
* zupt-threshold.cc : tool for test threshold used in zupt update
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/28.
*-----------------------------------------------------------------------------*/

#include "rtklib.h"

extern int detimustatic_test(const insopt_t* opt, const imud_t* imu, const double* pos, double dt, int det, int ws) {
	int n = opt->zvopt.ws; //opt->zvopt.det
	double* z1, * z2, * z3, * z4;
	z1 = new double;
	z2 = new double;
	z3 = new double;
	z4 = new double;
	detstatic_GLRT(imu, n, opt, pos, dt, z1);
	detstatic_MV(imu, n, opt, dt, z2);
	detstatic_MAG(imu, n, opt, pos, dt, z3);
	detstatic_ARE(imu, n, opt, dt, z4);
	fprintf(stdout, "%12.5f  ,%12.5f  ,%12.5f  ,%12.5f", *z1, *z2, *z3, *z4);
	fprintf(stdout, "\n");
	fflush(stdout);
	delete z1;
	delete z2;
	delete z3;
	delete z4;
	return 1;
}

int test_zero_threshold(kf_t* inskf, FILE* fp_sol, const insopt_t* opt, imud_t* imuz) {
	double ep[6];
	time2epoch(inskf->time, ep);
	fprintf(fp_sol, "%5.0f,%2.0f,%2.0f,%2.0f,%2.0f,%7.4f,", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
	double ep1[6] = { 2020, 12, 10, 3, 31, 0.0 };
	double ep2[6] = { 2020, 12, 10, 3, 34, 55.0 };
	gtime_t start = epoch2time(ep1);
	gtime_t end = epoch2time(ep2);
	int vflag = 0;
	double xyz[3], pos[3];
	v3_2_array(inskf->insstate->pos, xyz);
	ecef2pos(xyz, pos);
	gtime_t cur_time = imuz[1].time;
	double dt = timediff(cur_time, inskf->time);
	//    vflag = detimustatic_test(opt, imuz, pos, dt, fp_sol);
		//fprintf(fp, "%12.5f  ,%12.5f  ,%12.5f  ,%12.5f,%12.5f,%12.5f,", outsol.pos.x, outsol.pos.y, outsol.pos.z, outsol.vel.x, outsol.vel.y,outsol.vel.z);
	fprintf(fp_sol, "%12.5f,%12.5f,%12.5f,", inskf->insstate->vel.x, inskf->insstate->vel.y, inskf->insstate->vel.z);
	fprintf(fp_sol, "\n");
	fflush(fp_sol);
	return 0;
}