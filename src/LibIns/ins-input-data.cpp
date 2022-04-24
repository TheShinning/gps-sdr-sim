/*-----------------------------------------------------------------------------
* ins-input-data.cc : read imu data from file and reform data to structure
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/5/23.
*----------------------------------------------------------------------------*/

#include "rtklib.h"
#include <iostream>
#include <fstream>

/* constants/macros ----------------------------------------------------------*/
#define NUMBYTES_GI310  43                      /* numbers of bytes of gi310 imu raw data */
#define MAXDIFFTIME     10.0                    /* max time difference to reset  */
#define FREQOCXO    1E8                 /* crystal frequency (100Mhz) */

using namespace std;

const unsigned char gi310_head[2] = { 0x55,0xAA };  /* imu message header */

/* get fields (little-endian) ------------------------------------------------*/
#define U1(p) (*((unsigned char *)(p)))
#define I1(p) (*((char *)(p)))

/* get fields (little-endian) ------------------------------------------------*/
static unsigned short U2(unsigned char* p) { unsigned short u; memcpy(&u, p, 2); return u; }
static unsigned int   U4(unsigned char* p) { unsigned int   u; memcpy(&u, p, 4); return u; }
static short          I2(unsigned char* p) { short          i; memcpy(&i, p, 2); return i; }
static int            I4(unsigned char* p) { int            i; memcpy(&i, p, 4); return i; }
static float          R4(unsigned char* p) { float          r; memcpy(&r, p, 4); return r; }
static double         R8(unsigned char* p) { double         r; memcpy(&r, p, 8); return r; }

static int cmpimu(const void* p1, const void* p2)
{
	imud_t* q1 = (imud_t*)p1, * q2 = (imud_t*)p2;
	double tt = timediff(q1->time, q2->time);
	if (fabs(tt) > INS_DTTOL) return tt < 0 ? -1 : 1;
	return 0;
}

extern int sortimudata(imu_t* imu)
{
	int i, j, n;

	trace(3, "sortobs: nobs=%d\n", imu->n);

	if (imu->n <= 0) return 0;

	qsort(imu->data, imu->n, sizeof(imud_t), cmpimu);

	/* delete duplicated data */
	for (i = j = 0; i < imu->n; i++) {
		if (timediff(imu->data[i].time, imu->data[j].time) != 0.0) {
			imu->data[++j] = imu->data[i];
		}
	}
	imu->n = j + 1;

	for (i = n = 0; i < imu->n; i = j, n++) {
		for (j = i + 1; j < imu->n; j++) {
			if (timediff(imu->data[j].time, imu->data[i].time) > INS_DTTOL) break;
		}
	}
	return n;
}

static void adjustimu(imud_t* imu, const imup_t* imup)
{
	double gyro[3], accl[3], dt = 1.0 / imup->freq_imu;
	int i;

	trace(5, "adjustimu:\n");

	if (imup->imucoors == IMUCOOR_RFU) { /* convert to frd-ned-frame */
		matcpy(gyro, imu->gyro.v, 1, 3);
		matcpy(accl, imu->accel.v, 1, 3);
		matmul("NN", 3, 1, 3, 1.0, Crf, gyro, 0.0, imu->gyro.v);
		matmul("NN", 3, 1, 3, 1.0, Crf, accl, 0.0, imu->accel.v);
	}
	if (imup->imudecfmt == IMUDECFMT_INCR) {
		for (i = 0; i < 3; i++) {
			imu->gyro.v[i] /= dt;
			imu->accel.v[i] /= dt; /* convert to rate/acceleration */
		}
	}
	if (imup->imuvalfmt == IMUVALFMT_DEG) {
		for (i = 0; i < 3; i++) imu->gyro.v[i] *= D2R; /* convert to rad */
	}
}

static void adj_imudata(imu_t* imu)
{
	int i;
	for (i = 0; i < imu->n; i++) adjustimu(&imu->data[i], imu->property);
}

/* check header --------------------------------------------------------------*/
static int chkhead(const unsigned char* buff, const unsigned char* head)
{
	return (buff[0] == head[0]) && (buff[1] == head[1]);
}
/* checksum ------------------------------------------------------------------*/
static unsigned char chksum(const unsigned char* buff, int len)
{
	int i;
	unsigned char sum = 0;
	for (i = 0; i < len; i++) sum += buff[i]; return sum;
}
/* decode imu time------------------------------------------------------------*/
static void decode_sow_time(raw_t* raw, double* sow, int* start)
{
	static unsigned int pps = 0;

	if (pps == 0) {
		pps = U4(raw->buff + 6); *sow = U4(raw->buff + 2);
		return;
	}
	if (pps != U4(raw->buff + 6)) { (*sow)++; *start = 1; }
	pps = U4(raw->buff + 6);
}
/* adjust gps seconds of week and imu time------------------------------------*/
static void adjtime(raw_t* raw, const double sowi, double* sowo, double* timu,
	int* week, unsigned int* dcc)
{
	static unsigned int imuc = 0, dc = 0;
	int d = 0;

	if (imuc == 0) {
		imuc = U4(raw->buff + 10); *sowo = sowi;
		return;
	}
	*dcc = dc = U4(raw->buff + 10) - imuc < 0 ? UINT_MAX + U4(raw->buff + 10) - imuc : U4(raw->buff + 10) - imuc;
	imuc = U4(raw->buff + 10);

	/* increase week */
	if ((*sowo = (int)(1.0 / FREQOCXO * dc + sowi)) >= 604800.0) {
		*sowo -= 604800.0; (*week)++;
	}
	d = U4(raw->buff + 10) - U4(raw->buff + 06);
	*timu = *sowo + 1.0 / FREQOCXO * d;
}

static int addimudata(imu_t* imu, const imud_t* data)
{
	imud_t* obs_data;

	if (imu->nmax <= imu->n) {
		if (imu->nmax <= 0) imu->nmax = 64; else imu->nmax *= 2;
		if (!(obs_data = (imud_t*)realloc(imu->data, sizeof(imud_t) * imu->nmax))) {
			trace(1, "addimudata: memalloc error n=%dx%d\n", sizeof(imud_t), imu->nmax);
			free(imu->data); imu->data = NULL; imu->n = imu->nmax = 0;
			return -1;
		}
		imu->data = obs_data;
	}
	imu->data[imu->n++] = *data;
	return 1;
}

static int decode_odo_data(raw_t* raw, double dt)
{
	static int dc;

	if (raw->imu.time.time == 0 || dt == 0.0) {
		raw->imu.odoc = I2(raw->buff + 38); return 0;
	}
	dc = I2(raw->buff + 38) - raw->imu.odoc <= SHRT_MIN ?
		I2(raw->buff + 38) - raw->imu.odoc + USHRT_MAX :
		I2(raw->buff + 38) - raw->imu.odoc >= SHRT_MAX ?
		I2(raw->buff + 38) - raw->imu.odoc - USHRT_MAX :
		I2(raw->buff + 38) - raw->imu.odoc;

	raw->imu.odoc = I2(raw->buff + 38);
	//    raw->imu.odo.time=raw->imu.time;
	//    raw->imu.odo.dt=dt;
	//    raw->imu.odo.dr=dc/res*PI*d;
	return 1;
}

static int decode_imu_data(raw_t* raw)
{
	int i, week = 0;
	unsigned int dc = 0;
	static int start = 0;
	static double sow = 0.0, timu = 0.0;
	gtime_t t0;

	raw->imut.n = 0;

	/* decode GPS sow (s) */
	decode_sow_time(raw, &sow, &start);

	/* start decode imu time */
	if (start) {
		adjtime(raw, sow, &sow, &timu, &week, &dc);
	}
	else return 0;

	/* current and precious time difference is too large */
	if (dc * 1.0 / FREQOCXO > MAXDIFFTIME) return 0;

	t0 = gpst2time(week, timu);
	raw->imu.pps = U4(raw->buff + 06);
	raw->imu.imuc = U4(raw->buff + 10);

	for (i = 0; i < 3; i++) {
		raw->imu.gyro.v[i] = R4(raw->buff + 14 + i * 4);
		raw->imu.accel.v[i] = R4(raw->buff + 26 + i * 4);
	}
	//    decode_odo_data(raw,dc*1.0/FREQOCXO);
	raw->imu.time = t0;

	/* add imu measurement data */
	addimudata(&raw->imut, &raw->imu);
	return timu > 0.0 ? 4 : 0;
}

static int nextimub(const imu_t* imu, imu_t* data, int* index)
{
	int i; data->n = 0;
	for (i = 0; i < 100; i++) addimudata(data, &imu->data[(*index)--]);
	if (data->n) return 4; else return 0;
}

static int decode_imu_m39b(raw_t* raw, unsigned char data)
{
	if (raw->imub.n == 0) {
		prcopt_t* opt = (prcopt_t*)raw->optp;
		stream_t* str = (stream_t*)raw->strp;
		readimu_m39(str->path, &raw->imub, NULL);

		raw->curb = raw->imub.n - 1;
	}
	raw->imut.n = 0;
	return nextimub(&raw->imub, &raw->imut, &raw->curb);
}

static int decode_imu_m39(raw_t* raw, unsigned char data)
{
	raw->buff[raw->nbyte++] = data;

	if (raw->nbyte < 2) return 0; /* synchronize frame */

	if (!chkhead(raw->buff, gi310_head)) { raw->nbyte = 0; return 0; }
	if (raw->nbyte < NUMBYTES_GI310) return 0;

	if (chksum(raw->buff + 2, 40) == data) { /* checksum */
		raw->nbyte = 0;
		return decode_imu_data(raw);
	}
	else {
		raw->nbyte = 0; return -1; /* checksum fail */
	}
}

static int input_m39(raw_t* raw, unsigned char data)
{
	trace(5, "input_m39: data=%02X\n", data);
	raw->len = NUMBYTES_GI310;
	if (raw->dire) return decode_imu_m39b(raw, data);
	else           return decode_imu_m39(raw, data);
}

extern int readimu_m39(const char* file, imu_t* imu, FILE* fp)
{
	raw_t raw = { 0 };
	int data, siz, nimu;

	trace(3, "readimub:\n");

	imu->n = imu->nmax = 0; imu->data = NULL;

	/* read imu message from file */
	while (1) {
		if ((data = fgetc(fp)) == EOF) break;
		if ((input_m39(&raw, (unsigned char)data))) {
			if (imu->n >= imu->nmax) {
				trace(5, "readimub: imu->n=%d nmax=%d\n", imu->n, imu->nmax);

				siz = sizeof(imud_t) * (imu->nmax += 4096);
				if (!(imu->data = (imud_t*)realloc(imu->data, siz))) {
					fprintf(stderr, "readimub :memory allocation error\n");
					free(imu->data); imu->n = imu->nmax = 0;
					break;
				}
			}
			imu->data[imu->n++] = raw.imu;
		}
	}
	nimu = sortimudata(imu);
	adj_imudata(imu);
	return imu->n;
}

static int readf_imu_nvt(FILE* fp, imu_t* imu, const insopt_t* opt)
{
	if (!(imu->data = (imud_t*)malloc(sizeof(imud_t) * MAXIMUOBS))) {
		trace(1, "imu->data memory alloction error\n");
		return 0;
	}
	imu->n = 0;
	imu->nmax = MAXIMUOBS;

	char buff[MAXLINELEN];
	imud_t data;
	char seps[] = ",;*";
	char* token;
	int gpsw;

	double a_scale = 1.52587890625E-06;     /* 0.05*2^-15 */
	double g_scale = 1.085069444444444E-07; /* 0.1/(3600*256) */
	while (fgets(buff, MAXLINELEN, fp)) {
		token = strtok(buff, seps);
		if (!strcmp(token, "%RAWIMUSA")) {
			for (int i = 0; i < 3; ++i)
				token = strtok(NULL, seps); /* skip header */
			gpsw = atoi(token);
			token = strtok(NULL, seps);
			data.time = gpst2time(gpsw, atof(token));
			for (int i = 0; i < 2; ++i)
				token = strtok(NULL, seps); /* skip flag */

			if (opt->local_coord == INSLOCAL_NED) {
				/* velocity increment, m/s */
				data.accel.z = -atof(token) * a_scale;
				token = strtok(NULL, seps);
				data.accel.x = -atof(token) * a_scale;
				token = strtok(NULL, seps);
				data.accel.y = atof(token) * a_scale;
				token = strtok(NULL, seps);
				/* angle increment, rad */
				data.gyro.z = -atof(token) * g_scale;
				token = strtok(NULL, seps);
				data.gyro.x = -atof(token) * g_scale;
				token = strtok(NULL, seps);
				data.gyro.y = atof(token) * g_scale;
			}
			else if (opt->local_coord == INSLOCAL_ENU) {
				/* velocity increment, m/s */
				data.accel.z = atof(token) * a_scale;
				token = strtok(NULL, seps);
				data.accel.y = -atof(token) * a_scale;
				token = strtok(NULL, seps);
				data.accel.x = atof(token) * a_scale;
				token = strtok(NULL, seps);
				/* angle increment, rad */
				data.gyro.z = atof(token) * g_scale;
				token = strtok(NULL, seps);
				data.gyro.y = -atof(token) * g_scale;
				token = strtok(NULL, seps);
				data.gyro.x = atof(token) * g_scale;
			}

			imu_add(imu, &data);
		}
	}
	return 1;
}

static int readf_nvt_kvh(FILE* fp, imu_t* imu, const insopt_t* opt)
{  // 速率
	char buff[MAXLINELEN];
	imud_t imud;
	int week;
	double sec, gx, gy, gz, ax, ay, az;
	double dt = 1.0 / opt->imup.freq_imu;

	fgets(buff, MAXLINELEN, fp);
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &week, &sec, &gx, &gy, &gz, &ax, &ay, &az);
		if (opt->mech_coord == INSLOCAL_NED) {
			//imud.gyro = (v3_t){ gy * dt,gx * dt,-gz * dt };
			//imud.accel = (v3_t){ ay * dt,ax * dt,-az * dt };

			imud.gyro.x = gy * dt; imud.gyro.y = gx * dt; imud.gyro.z = -gz * dt;
			imud.accel.x = ay * dt; imud.accel.y = ax * dt; imud.accel.z = -az * dt;
		}
		else if (opt->mech_coord == INSLOCAL_ENU) {
			//imud.gyro = (v3_t){ gx * dt,gy * dt,gz * dt };
			//imud.accel = (v3_t){ ax * dt,ay * dt,az * dt };
			imud.gyro.x = gx * dt;   imud.accel.x = ax * dt;
			imud.gyro.y = gy * dt;   imud.accel.y = ay * dt;
			imud.gyro.z = gz * dt;   imud.accel.z = az * dt;
		}
		imud.time = gpst2time(week, sec);
		imu_add(imu, &imud);
	}
	return 0;
}

static int readf_imu_increment(const char* file, imu_t* imu, const insopt_t* opt)
{  // 增量
	unsigned char buff[MAXLINELEN];
	imud_t imud;
	int week;
	double sec, gx, gy, gz, ax, ay, az;
	double data[8];
	ifstream imuFile(file, ios::in | ios::binary);

	while (1) {
		imuFile.read((char*)&data, sizeof(double) * 8);
		week = int(data[0]); sec = data[1];
		gx = data[2]; gy = data[3]; gz = data[4];
		ax = data[5]; ay = data[6]; az = data[7];
		if (imuFile.eof()) return 0;
		if (opt->local_coord == INSLOCAL_NED) {
			//imud.gyro = (v3_t){ gx,gy,gz };
			//imud.accel = (v3_t){ ax,ay,az };
			imud.gyro.x = gx;   imud.accel.x = ax;
			imud.gyro.y = gy;   imud.accel.y = ay;
			imud.gyro.z = gz;   imud.accel.z = az;
		}
		else if (opt->local_coord == INSLOCAL_ENU) {
			//imud.gyro = (v3_t){ gy,gx,-gz };
			//imud.accel = (v3_t){ ay,ax,-az };
			imud.gyro.x = gy;    imud.accel.x = ay;
			imud.gyro.y = gx;    imud.accel.y = ax;
			imud.gyro.z = -gz;   imud.accel.z = -az;
		}
		imud.time = gpst2time(week, sec);
		imu_add(imu, &imud);
	}
	imuFile.close();
	return 0;
}

static int readf_nvt(FILE* fp, imu_t* imu, const insopt_t* opt)
{
	char buff[MAXLINELEN];
	imud_t imud;
	int week = imu->property->week;
	double sec, gx, gy, gz, ax, ay, az;
	double dt = 1.0 / opt->imup.freq_imu;

	fgets(buff, MAXLINELEN, fp);
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf", &sec, &gx, &gy, &gz, &ax, &ay, &az);
		if (opt->imu_coord == IMUCOOR_FRD) {
			//imud.gyro = (v3_t){ gy * dt * D2R,gx * dt * D2R,-gz * dt * D2R };
			//imud.accel = (v3_t){ ay * dt,ax * dt,-az * dt };
			imud.gyro.x = gy * dt * D2R;   imud.accel.x = ay * dt;
			imud.gyro.y = gx * dt * D2R;   imud.accel.y = ax * dt;
			imud.gyro.z = -gz * dt * D2R;   imud.accel.z = -az * dt;
		}
		else if (opt->imu_coord == IMUCOOR_RFU) {
			//imud.gyro = (v3_t){ gx * dt * D2R,gy * dt * D2R,gz * dt * D2R };
			//imud.accel = (v3_t){ ax * dt,ay * dt,az * dt };
			imud.gyro.x = gx * dt * D2R;   imud.accel.x = ax * dt;
			imud.gyro.y = gy * dt * D2R;   imud.accel.y = ay * dt;
			imud.gyro.z = gz * dt * D2R;   imud.accel.z = az * dt;
		}
		imud.time = gpst2time(week, sec);
		imu_add(imu, &imud);
	}
	return 0;
}

static int readf_by(FILE* fp, imu_t* imu, const insopt_t* opt)
{
	char buff[MAXLINELEN];
	imud_t imud;
	int week = imu->property->week;
	double sec, gx, gy, gz, ax, ay, az;
	double dt = 1.0 / opt->imup.freq_imu;

	fgets(buff, MAXLINELEN, fp);
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf", &sec, &gx, &gy, &gz, &ax, &ay, &az);
		if (opt->imu_coord == IMUCOOR_FRD) {
			//imud.gyro = (v3_t){ gz * dt * D2R,gy * dt * D2R,-gx * dt * D2R };
			//imud.accel = (v3_t){ az * dt,ay * dt,-ax * dt };
			imud.gyro.x = gz * dt * D2R;   imud.accel.x = az * dt;
			imud.gyro.y = gy * dt * D2R;   imud.accel.y = ay * dt;
			imud.gyro.z = -gx * dt * D2R;   imud.accel.z = -ax * dt;
		}
		else if (opt->imu_coord == IMUCOOR_RFU) {
			//imud.gyro = (v3_t){ gy * dt * D2R,gz * dt * D2R,gx * dt * D2R };
			//imud.accel = (v3_t){ ay * dt,az * dt,ax * dt };
			imud.gyro.x = gy * dt * D2R;   imud.accel.x = ay * dt;
			imud.gyro.y = gz * dt * D2R;   imud.accel.y = az * dt;
			imud.gyro.z = gx * dt * D2R;   imud.accel.z = ax * dt;
		}
		imud.time = gpst2time(week, sec);
		imu_add(imu, &imud);
	}
	return 0;
}

static int str2time_atlan(const char* s, int i, int n, gtime_t* t)
{
	double ep[6];
	char str[256], * p = str;

	if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < i) return -1;
	for (s += i; *s && --n >= 0;) *p++ = *s++;
	*p = '\0';
	if (sscanf(str, "%lf/%lf/%lf %lf:%lf:%lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5) < 6)
		return -1;
	if (ep[0] < 100.0) ep[0] += ep[0] < 80.0 ? 2000.0 : 1900.0;
	*t = epoch2time(ep);
	return 0;
}

static int readf_atlan(FILE* fp, imu_t* imu, int opt)
{
	char buff[MAXLINELEN];
	imud_t imud;
	double gx, gy, gz, ax, ay, az;
	char t1[40] = { '\0' }, t2[40] = { '\0' };
	int frame;
	double scale_gyro = 10E-3, scale_acc = 10E-7;
	while (fgets(buff, MAXLINELEN, fp)) {
		if (strstr(buff, "#")) {
			continue;
		}
		sscanf(buff, "%s %s %d %lf %lf %lf %lf %lf %lf",
			t1, t2, &frame, &gx, &gy, &gz, &ax, &ay, &az); /*RFU*/
		strcat(t1, " ");
		strcat(t1, t2);
		str2time_atlan(t1, 0, 28, &imud.time);
		if (opt == IMUCOOR_FRD) {
			//imud.gyro = (v3_t){ gy * scale_gyro / 3600 * D2R, gx * scale_gyro / 3600 * D2R, -gz * scale_gyro / 3600 * D2R };
			//imud.accel = (v3_t){ ay * scale_acc, ax * scale_acc, -az * scale_acc };

			imud.gyro.x = gy * scale_gyro / 3600 * D2R; imud.gyro.y = gx * scale_gyro / 3600 * D2R; imud.gyro.z = -gz * scale_gyro / 3600 * D2R;
			imud.accel.x = ay * scale_acc;              imud.accel.y = ax * scale_acc;                 imud.accel.z = -az * scale_acc;
		}
		else if (opt == IMUCOOR_RFU) {
			//imud.gyro = (v3_t){ gx * scale_gyro / 3600 * D2R, gy * scale_gyro / 3600 * D2R, gz * scale_gyro / 3600 * D2R };
			//imud.accel = (v3_t){ ax * scale_acc, ay * scale_acc,az * scale_acc };
			imud.gyro.x = gx * scale_gyro / 3600 * D2R; imud.gyro.y = gy * scale_gyro / 3600 * D2R;    imud.gyro.z = gz * scale_gyro / 3600 * D2R;
			imud.accel.x = ax * scale_acc;              imud.accel.y = ay * scale_acc;                 imud.accel.z = az * scale_acc;
		}
		imu_add(imu, &imud);
	}
	return 0;
}

static int readf_ygm_imu(FILE* fp, imu_t* imu, int opt)
{
	char buff[MAXLINELEN];
	imud_t imud;
	double gx, gy, gz, ax, ay, az, t;
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf",
			&gx, &gy, &gz, &ax, &ay, &az, &t);

		imud.time = gpst2time(imu->property->week, t + imu->property->sow);
		if (opt == IMUCOOR_FRD) {
			//imud.gyro = (v3_t){ gy, gx, -gz };
			//imud.accel = (v3_t){ ay, ax, -az };
			imud.gyro.x = gy; imud.gyro.y = gx; imud.gyro.z = -gz;
			imud.accel.x = ay; imud.accel.y = ax; imud.accel.z = -az;
		}
		else if (opt == IMUCOOR_RFU) {
			//imud.gyro = (v3_t){ gx, gy, gz };
			//imud.accel = (v3_t){ ax, ay,az };
			imud.gyro.x = gx;  imud.gyro.y = gy; imud.gyro.z = gz;
			imud.accel.x = ax; imud.accel.y = ay; imud.accel.z = az;
		}
		imu_add(imu, &imud);
	}
	return 0;
}

static int readf_ygm_avp(FILE* fp, pva_t* pva, const insopt_t* opt, int week)
{
	char buff[MAXLINELEN];
	double lat = 0.0, lon = 0.0, hgt = 0.0;
	double vE = 0.0, vN = 0.0, vU = 0.0;
	double pitch = 0.0, roll = 0.0, yaw = 0.0;
	double t;
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf",
			&vE, &vN, &vU, &lat, &lon, &hgt, &t);

		if (pva->n > pva->nmax) {
			pva->nmax *= 2;
			pva_resize(pva, pva->nmax);
		}

		pva->time[pva->n] = gpst2time(week, t);

		//pva->pos[pva->n] = (v3_t){ lat,lon,hgt };
		pva->pos[pva->n].x = lat;
		pva->pos[pva->n].y = lon;
		pva->pos[pva->n].z = hgt;
		if (opt->local_coord == INSLOCAL_ENU) {
			//pva->vel[pva->n] = (v3_t){ vE,vN,vU };
			//pva->att[pva->n] = (v3_t){ roll,pitch,yaw };

			pva->vel[pva->n].x = vE; pva->vel[pva->n].y = vN; pva->vel[pva->n].z = vU;
			pva->att[pva->n].x = roll; pva->att[pva->n].y = pitch; pva->att[pva->n].z = yaw;
		}
		else if (opt->local_coord == INSLOCAL_NED) {
			//pva->vel[pva->n] = (v3_t){ vN,vE,-vU };
			//pva->att[pva->n] = (v3_t){ pitch,roll,yaw };

			pva->vel[pva->n].x = vN; pva->vel[pva->n].y = vE; pva->vel[pva->n].z = -vU;
			pva->att[pva->n].x = pitch; pva->att[pva->n].y = roll; pva->att[pva->n].z = yaw;
		}

		if (opt->mech_coord == INSMECH_ECEF) {
			llh2ecef(&pva->pos[pva->n], &pva->vel[pva->n], nullptr, opt->local_coord);
		}

		pva->status[pva->n] = GNSS_STATUS_POS;
		pva->n++;
	}
	return 0;
}

static int readf_ygm_od(FILE* fp, od_t* od)
{
	char buff[MAXLINELEN];
	double dS, t;
	gtime_t time; time.time = 0.0;
	while (fgets(buff, MAXLINELEN, fp)) {
		sscanf(buff, "%lf %lf", &dS, &t);
		time.sec = t;
		od_add(od, &time, &dS);
	}
	return 0;
}

static int readf_pos(const insopt_t* iopt, FILE* fp, pva_t* pva)
{
	solopt_t opt = solopt_default;
	solbuf_t sols = { 0 };
	gtime_t ts = { 0 }, te = { 0 };

	readsolopt(fp, &opt);
	rewind(fp);

	if (!readsoldata(fp, ts, te, 0.0, 0, &opt, &sols, iopt->mech_coord)) {
		trace(2, "readsolt: no solution\n");
		return 0;
	}

	for (int i = 0; i < sols.n; i++) {
		if (pva->n > pva->nmax) {
			pva->nmax *= 2;
			pva_resize(pva, pva->nmax);
		}

		sol2pva(iopt, pva->n, &sols.data[i], pva);
		pva->n++;
	}

	return 1;
}

inline static void tokenconvert(char* token, v3_t* imup)
{
	if ((token = strtok(NULL, ","))) imup->x = atof(token);
	if ((token = strtok(NULL, ","))) imup->y = atof(token);
	if ((token = strtok(NULL, ","))) imup->z = atof(token);
}

static int readf_imup(FILE* fp, imup_t* imup, const insopt_t* opt)
{
	char buff[MAXLINELEN];
	v3_t temp;
	while (fgets(buff, MAXLINELEN, fp)) {
		if (buff[0] == '>') break;
		if (buff[0] == '#') continue;

		/* remove comment part from buff(after #) */
		for (unsigned int i = 0; i < strlen(buff); ++i) {
			if (buff[i] == '#') buff[i] = '\0';
		}
		char* token = strtok(buff, ":");
		if (imup != nullptr) {
			if (!strcmp(token, "imup.freq_imu")) {
				if ((token = strtok(nullptr, ",")))
					imup->freq_imu = (unsigned int)atoi(token);
			}
			else if (!strcmp(token, "imup.strfmt")) {
				if ((token = strtok(nullptr, ",")))
					imup->strfmt = (int)atoi(token) + STRFMT_IMU_YGM_SIM;
			}
			else if (!strcmp(token, "imup.decfmt")) {
				if ((token = strtok(nullptr, ",")))
					imup->imudecfmt = (unsigned int)atoi(token);
			}
			else if (!strcmp(token, "imup.coords")) {
				if ((token = strtok(nullptr, ",")))
					imup->imucoors = (unsigned int)atoi(token);
			}
			else if (!strcmp(token, "imup.valfmt")) {
				if ((token = strtok(nullptr, ",")))
					imup->imuvalfmt = (unsigned int)atoi(token);
			}
			else if (!strcmp(token, "imup.week")) {
				if ((token = strtok(nullptr, ",")))
					imup->week = atoi(token);
			}
			else if (!strcmp(token, "imup.sow")) {
				if ((token = strtok(nullptr, ",")))
					imup->sow = atof(token);
			}
			else if (!strcmp(token, "imup.freq_od")) {
				if ((token = strtok(nullptr, ",")))
					imup->freq_od = (unsigned int)atoi(token);
			}
			else if (!strcmp(token, "imup.accel_noise")) {
				tokenconvert(token, &imup->accel_noise);
			}
			else if (!strcmp(token, "imup.gyro_noise")) {    /* deg/s */
				tokenconvert(token, &imup->gyro_noise);
				imup->gyro_noise = v3_scalar(D2R, imup->gyro_noise);
			}
			else if (!strcmp(token, "imup.init_tag")) {
				int week;
				double sow;
				if ((token = strtok(nullptr, ","))) {
					week = atoi(token);
					if ((token = strtok(nullptr, ","))) {
						sow = atof(token);
					}
				}
				imup->init_tag = gpst2time(week, sow);
			}
			else if (!strcmp(token, "imup.init_pos")) {
				int flag = 0;
				if ((token = strtok(nullptr, ","))) {
					flag = atoi(token);
				}
				tokenconvert(token, &temp);
				if (flag == 0) {
					temp.x *= D2R;
					temp.y *= D2R;
				}
				if (opt->mech_coord == INSMECH_ECEF && flag != 2) {
					pos2ecef(temp.v, imup->initr.v);
				}
				else {
					imup->initr = temp;
				}
			}
			else if (!strcmp(token, "imup.init_vel")) {
				tokenconvert(token, &imup->initv);
				if (opt->mech_coord == INSMECH_ECEF) {
					enu2ecef(imup->initr.v, imup->initv.v, imup->initv.v);
				}
			}
			else if (!strcmp(token, "imup.init_att")) {
				int flag = 0;
				if ((token = strtok(nullptr, ","))) {
					flag = atoi(token);
				}
				tokenconvert(token, &imup->inita);
				if (flag == 0) {
					imup->inita = v3_scalar(D2R, imup->inita);
				}
			}
			else if (!strcmp(token, "imup.pos_err")) {
				tokenconvert(token, &temp);
				if (opt->mech_coord == INSMECH_ECEF) {
					imup->initr_err = temp;
				}
				else {
					temp.x /= RE_WGS84;
					temp.y /= RE_WGS84;
					imup->initr_err = temp;
				}
			}
			else if (!strcmp(token, "imup.vel_err")) {
				tokenconvert(token, &imup->initv_err);
			}
			else if (!strcmp(token, "imup.att_err")) {
				tokenconvert(token, &temp);
				temp.x *= D2R;
				temp.y *= D2R;
				temp.z *= D2R;
				imup->inita_err = temp;
			}
			else if (!strcmp(token, "imup.ba")) {             /* mg */
				tokenconvert(token, &imup->ba);
				imup->ba = v3_scalar(MG2MPS2, imup->ba);
			}
			else if (!strcmp(token, "imup.ba_err")) {        /* mg */
				tokenconvert(token, &imup->ba_err);
				imup->ba_err = v3_scalar(MG2MPS2, imup->ba_err);
			}
			else if (!strcmp(token, "imup.bg")) {            /* deg/h */
				tokenconvert(token, &imup->bg);
				imup->bg = v3_scalar(DPH2RPS, imup->bg);
			}
			else if (!strcmp(token, "imup.bg_err")) {        /* deg/h */
				tokenconvert(token, &imup->bg_err);
				imup->bg_err = v3_scalar(DPH2RPS, imup->bg_err);
			}
			else if (!strcmp(token, "imup.sa")) {
				tokenconvert(token, &imup->sa);
			}
			else if (!strcmp(token, "imup.sa_err")) {
				tokenconvert(token, &imup->sa_err);
			}
			else if (!strcmp(token, "imup.sg")) {
				tokenconvert(token, &imup->sg);
			}
			else if (!strcmp(token, "imup.sg_err")) {
				tokenconvert(token, &imup->sg_err);
			}
			else if (!strcmp(token, "imup.gyrnd")) {        /* deg/sqrt(h) */
				tokenconvert(token, &imup->gyrnd);
				imup->gyrnd = v3_scalar(DPSH2RPSS, imup->gyrnd);
			}
			else if (!strcmp(token, "imup.gbrw")) {
				tokenconvert(token, &imup->gbrw);
			}
			else if (!strcmp(token, "imup.accnd")) {       /* mg/sqrt(Hz) */
				tokenconvert(token, &imup->accnd);
				imup->accnd = v3_scalar(MG2MPS2, imup->accnd);
			}
			else if (!strcmp(token, "imup.abrw")) {
				tokenconvert(token, &imup->abrw);
			}
			else if (!strcmp(token, "imup.Ta")) {
				tokenconvert(token, &imup->Ta);
			}
			else if (!strcmp(token, "imup.Tg")) {
				tokenconvert(token, &imup->Tg);
			}
			else if (!strcmp(token, "imup.kod")) {
				if ((token = strtok(nullptr, ","))) imup->kod = atof(token);
			}
			else if (!strcmp(token, "imup.kod_err")) {
				if ((token = strtok(nullptr, ","))) imup->kod_err = atof(token);
			}
			else if (!strcmp(token, "imup.lever_arm_gps")) {
				tokenconvert(token, &temp);
				if (opt->imu_coord == IMUCOOR_RFU) {
					imup->lever_arm_gps.x = -temp.y;   /*imu to gnss*/
					imup->lever_arm_gps.y = -temp.x;
					imup->lever_arm_gps.z = temp.z;
				}
				else if (opt->imu_coord == IMUCOOR_FRD) {
					imup->lever_arm_gps = temp;
				}
			}
			else if (!strcmp(token, "imup.lever_arm_gps_std")) {
				tokenconvert(token, &imup->lever_arm_gps_std);
			}
			else if (!strcmp(token, "imup.lever_arm_od")) {
				tokenconvert(token, &imup->lever_arm_od);
			}
			else if (!strcmp(token, "imup.lever_arm_od_std")) {
				tokenconvert(token, &imup->lever_arm_od_std);
			}
			else if (!strcmp(token, "imup.lever_arm_car")) {
				tokenconvert(token, &imup->lever_arm_car);
			}
			else if (!strcmp(token, "imup.lever_arm_car_std")) {
				tokenconvert(token, &imup->lever_arm_car_std);
			}
			else if (!strcmp(token, "imup.err_angle_imu")) {     /* deg */
				tokenconvert(token, &imup->err_angle_imu);
				imup->err_angle_imu = v3_scalar(D2R, imup->err_angle_imu);
			}
			else if (!strcmp(token, "imup.err_angle_imu_std")) {   /* deg */
				tokenconvert(token, &imup->err_angle_imu_std);
				imup->err_angle_imu_std =
					v3_scalar(D2R, imup->err_angle_imu_std);
			}
			else if (!strcmp(token, "imup.err_angle_imu_rw")) {
				tokenconvert(token, &imup->err_angle_imu_rw);
			}
			else if (!strcmp(token, "imup.Terr_angle_imu")) {
				tokenconvert(token, &imup->Terr_angle_imu);
			}
			else if (!strcmp(token, "imup.err_angle_gps")) {
				tokenconvert(token, &imup->err_angle_gps);
				imup->err_angle_gps = v3_scalar(D2R, imup->err_angle_gps);
			}
			else if (!strcmp(token, "imup.err_angle_gps_std")) {
				tokenconvert(token, &imup->err_angle_gps_std);
				imup->err_angle_gps_std =
					v3_scalar(D2R, imup->err_angle_gps_std);
			}
			else if (!strcmp(token, "imup.ref_point")) {
				tokenconvert(token, &imup->ref_point);
			}
		}
	}
	return 1;
}

extern int ins_readf(const char* fname, int ft, imu_t* imu, pva_t* pva, od_t* od, const insopt_t* opt)
{
	FILE* fp;
	if (!(fp = fopen(fname, "r"))) {
		trace(1, "%s open error\n", fname);
		return 0;
	}

	switch (ft) {
	case STRFMT_IMU_PROPERTY:
		readf_imup(fp, imu->property, opt); break;
	case STRFMT_IMU_YGM_SIM:
		readf_ygm_imu(fp, imu, opt->imu_coord); break;
	case STRFMT_IMU_M39:
		readimu_m39(fname, imu, fp); break;
	case STRFMT_IMU_CPT_TOKYO:
		readf_nvt_kvh(fp, imu, opt); break;           // vel
	case STRFMT_IMU_A15:
		readf_imu_increment(fname, imu, opt); break;  // increment
	case STRFMT_IMU_NVT_IGM:
	case STRFMT_IMU_NVT_KVH:
	case STRFMT_IMU_NVT_CPT:
		readf_nvt(fp, imu, opt); break;
	case STRFMT_IMU_BY:
		readf_by(fp, imu, opt); break;
	case STRFMT_IMU_ATLAN:
		readf_atlan(fp, imu, opt->imu_coord); break;
	case STRFMT_YGM_PV:
		readf_ygm_avp(fp, pva, opt, imu->property->week); break;
	case STRFMT_POS:
		readf_pos(opt, fp, pva); break;
	case STRFMT_OD_YGM_SIM:
		readf_ygm_od(fp, od); break;

	default:
		trace(1, "file type error %s\n", fname);
		fclose(fp);
		return 0;
	}

	fclose(fp);
	return 1;
}