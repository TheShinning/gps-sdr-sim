#include <stdlib.h>
#include <signal.h>
//#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>//#include <sys/time.h>

#include  <winsock2.h>//#include <sys/socket.h>

//#include <netinet/in.h>

//#include <netinet/tcp.h>

#include<windows.h>// <arpa/inet.h>

//#include <netdb.h>
//#include"json.h"

#include <iostream>
#include <fstream>

//#include "vt.h"
#include <errno.h>
#include "rtklib.h"

//#include "pthread.h"
#include "WS2tcpip.h"

#define PRGNAME     "rtkrcv"            /* program name */
#define CMDPROMPT   "rtkrcv> "          /* command prompt */
#define MAXCON      32                  /* max number of consoles */
#define MAXARG      10                  /* max number of args in a command */
#define MAXCMD      256                 /* max length of a command */
#define MAXSTR      1024                /* max length of a stream */
#define OPTSDIR     "."                 /* default config directory */
#define OPTSFILE    "rtkrcv.conf"       /* default config file */
#define NAVIFILE    "rtkrcv.nav"        /* navigation save file */
#define STATFILE    "rtkrcv_%Y%m%d%h%M.stat"  /* solution status file */
#define TRACEFILE   "rtkrcv_%Y%m%d%h%M.trace" /* debug trace file */
#define INTKEEPALIVE 1000               /* keep alive interval (ms) */

#define ESC_CLEAR   "\033[H\033[2J"     /* ansi/vt100 escape: erase screen */
#define ESC_RESET   "\033[0m"           /* ansi/vt100: reset attribute */
#define ESC_BOLD    "\033[1m"           /* ansi/vt100: bold */

//#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static rtksvr_t svr;                    /* rtk server struct */
static stream_t moni;                   /* monitor stream */

static int intflg = 0;             /* interrupt flag (2:shtdown) */

static char passwd[MAXSTR] = "admin";     /* login password */
static int timetype = 0;             /* time format (0:gpst,1:utc,2:jst,3:tow) */
static int soltype = 0;             /* sol format (0:dms,1:xyz,2:enu,3:pyl,7:deg,) *///0:llh, 1 : xyz, 2 : enu, 3 : nmea, 4 : stat, 5 : gsif, 6 : ins
static int solflag = 2;             /* sol flag (1:std+2:age/ratio/ns) */
static int strtype[] = {                  /* stream types */
	STR_SERIAL,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,
	STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE
};
static char strpath[18][MAXSTR] = { "","","","","","","","","","","","","","","","","","" }; /* stream paths */
static int strfmt[] = {                   /* stream formats */
	//STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,SOLF_LLH,SOLF_NMEA,STRFMT_NONE,STRFMT_NONE,STRFMT_NONE,STRFMT_NONE,STRFMT_NONE
	STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,SOLF_LLH,SOLF_NMEA
	,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3
};
static int stream_number = 4;       /* the total number of stream. Note:1:base 2:eph or other type 3:output stream  4 to end: rover stream */
static int svrcycle = 10;            /* server cycle (ms) */
static int timeout = 10000;         /* timeout time (ms) */
static int reconnect = 10000;         /* reconnect interval (ms) */
static int nmeacycle = 5000;          /* nmea request cycle (ms) */
static int buffsize = 32768;         /* input buffer size (bytes) */
static int navmsgsel = 0;             /* navigation mesaage select */
static char proxyaddr[256] = "";          /* http/ntrip proxy */
static int nmeareq = 0;             /* nmea request type (0:off,1:lat/lon,2:single) */
static double nmeapos[] = { 0,0,0 };       /* nmea position (lat/lon/height) (deg,m) */
static char rcvcmds[10][MAXSTR] = { "" };    /* receiver commands files */
static char startcmd[MAXSTR] = "";        /* start command */
static char stopcmd[MAXSTR] = "";        /* stop command */
static int modflgr[256] = { 0 };           /* modified flags of receiver options */
static int modflgs[256] = { 0 };           /* modified flags of system options */
static int moniport = 0;             /* monitor port */
static int keepalive = 0;             /* keep alive flag */
static int start = 0;             /* auto start */
static int fswapmargin = 30;            /* file swap margin (s) */
static char sta_name[256] = "";           /* station name */

static prcopt_t prcopt;                 /* processing options */
static solopt_t solopt[2] = { {0} };        /* solution options */
static filopt_t filopt = { "" };          /* file options */

static void prtime(gtime_t time)
{
	double tow;
	int week;
	char tstr[64] = "";

	if (timetype == 1) {
		time2str(gpst2utc(time), tstr, 2);
	}
	else if (timetype == 2) {
		time2str(timeadd(gpst2utc(time), 9 * 3600.0), tstr, 2);
	}
	else if (timetype == 3) {
		tow = time2gpst(time, &week); sprintf(tstr, "  %04d %9.2f", week, tow);
	}
	else time2str(time, tstr, 1);
	printf("%s\n", tstr);
}

static void prstream()
{
	//const char *ch[] = {
	//	"input rover","input base","input corr","output sol1","output sol2",
	//	"log rover","log base","log corr","monitor"
	//};
	const char* ch[] = {
	"Rover","Ssr","Eph","output1","log rover","monitor"
	};
	const char* type[] = {
		"-","serial","file","tcpsvr","tcpcli","udp","ntrips","ntripc","ftp",
		"http","ntripcas"
	};
	const char* fmt[] = { "rtcm2","rtcm3","oem4","","ubx","swift","hemis","skytreq",
					   "javad","nvs","binex","rt17","sbf","","","sp3","" };
	const char* sol[] = { "llh","xyz","enu","nmea","stat","-" };
	stream_t stream[10];
	int i, format[10] = { 0 };

	trace(4, "prstream:\n");

	rtksvrlock(&svr);
	for (i = 0; i < 3; i++) stream[i] = svr.stream[i];
	for (i = 0; i < 3; i++) format[i] = svr.format[i];
	stream[3] = svr.stream[6];
	format[3] = svr.format[6];
	stream[4] = svr.stream[7];
	format[4] = svr.format[7];

	///*for (i = 6; i < 9; i++)*/ format[2] = svr.solopt[0].posf;
	//stream[8] = moni;
	//format[8] = SOLF_LLH;
	rtksvrunlock(&svr);

	printf("\n%s%-12s %-8s %-5s %s %10s %7s %10s %7s %-24s %s%s\n", ESC_BOLD,
		"Stream", "Type", "Fmt", "S", "In-byte", "In-bps", "Out-byte", "Out-bps",
		"Path", "Message", ESC_RESET);
	for (i = 0; i < 5; i++) {
		//vt_printf(vt, "%-12s %-8s %-5s %s %10d %7d %10d %7d %-24.24s %s\n",
		//	ch[i], type[stream[i].type], i < 3 ? fmt[format[i]] : (i < 5 || i == 8 ? sol[format[i]] : "-"),
		//	stream[i].state < 0 ? "E" : (stream[i].state ? "C" : "-"),
		//	stream[i].inb, stream[i].inr, stream[i].outb, stream[i].outr,
		//	stream[i].path, stream[i].msg);
		printf("%-12s %-8s %-5s %s %10d %7d %10d %7d %-24.24s %s\n",
			ch[i], type[stream[i].type], fmt[format[i]],
			stream[i].state < 0 ? "E" : (stream[i].state ? "C" : "-"),
			stream[i].inb, stream[i].inr, stream[i].outb, stream[i].outr,
			stream[i].path, stream[i].msg);
	}
}

static void prsolution(const sol_t* sol, const double* rb)
{
	const char* solstr[] = { "------","FIX","FLOAT","SBAS","DGPS","SINGLE","PPP","PPP-AR" };
	double pos[3] = { 0 }, Qr[9], Qe[9] = { 0 }, dms1[3] = { 0 }, dms2[3] = { 0 }, bl[3] = { 0 };
	double enu[3] = { 0 }, pitch = 0.0, yaw = 0.0, len;
	int i;

	trace(4, "prsolution:\n");
	prtime(sol->time);

	if (sol->time.time == 0 || !sol->stat) return;
	printf("\n(%-6s)", solstr[sol->stat]);

	if (norm(sol->rr, 3) > 0.0 && norm(rb, 3) >= 0.0) {
		for (i = 0; i < 3; i++) bl[i] = sol->rr[i] - rb[i];
	}
	len = norm(bl, 3);
	Qr[0] = sol->qr[0];
	Qr[4] = sol->qr[1];
	Qr[8] = sol->qr[2];
	Qr[1] = Qr[3] = sol->qr[3];
	Qr[5] = Qr[7] = sol->qr[4];
	Qr[2] = Qr[6] = sol->qr[5];

	if (soltype == 0) {
		if (norm(sol->rr, 3) > 0.0) {
			ecef2pos(sol->rr, pos);
			covenu(pos, Qr, Qe);
			deg2dms(pos[0] * R2D, dms1, 4);
			deg2dms(pos[1] * R2D, dms2, 4);
			if (solopt[0].height == 1) pos[2] -= geoidh(pos); /* geodetic */
		}
		printf(" %s:%2.0f %02.0f %07.4f", pos[0] < 0 ? "S" : "N", fabs(dms1[0]), dms1[1], dms1[2]);

		printf(" %s:%3.0f %02.0f %07.4f", pos[1] < 0 ? "W" : "E", fabs(dms2[0]), dms2[1], dms2[2]);

		printf(" H:%8.3f", pos[2]);
		if (solflag & 1) {
			printf(" (N:%6.3f E:%6.3f U:%6.3f)", SQRT(Qe[4]), SQRT(Qe[0]), SQRT(Qe[8]));
		}
	}
	else if (soltype == 4) {
		if (norm(sol->rr, 3) > 0.0) {
			ecef2pos(sol->rr, pos);
			covenu(pos, Qr, Qe);
			if (solopt[0].height == 1) pos[2] -= geoidh(pos); /* geodetic */
		}

		printf(" %s:%11.8f", pos[0] < 0.0 ? "S" : "N", fabs(pos[0]) * R2D);
		printf(" %s:%12.8f", pos[1] < 0.0 ? "W" : "E", fabs(pos[1]) * R2D);
		printf(" H:%8.3f", pos[2]);

		if (solflag & 1) {
			printf(" (E:%6.3f N:%6.3f U:%6.3fm)", SQRT(Qe[0]), SQRT(Qe[4]), SQRT(Qe[8]));
		}
	}
	else if (soltype == 1) {
		printf(" X:%12.3f", sol->rr[0]);
		printf(" Y:%12.3f", sol->rr[1]);
		printf(" Z:%12.3f", sol->rr[2]);

		if (solflag & 1) {
			printf(" (X:%6.3f Y:%6.3f Z:%6.3f)", SQRT(Qr[0]), SQRT(Qr[4]), SQRT(Qr[8]));
		}
	}
	else if (soltype == 2) {
		if (len > 0.0) {
			ecef2pos(rb, pos);
			ecef2enu(pos, bl, enu);
			covenu(pos, Qr, Qe);
		}

		printf(" E:%12.3f", enu[0]);
		printf(" N:%12.3f", enu[1]);
		printf(" U:%12.3f", enu[2]);

		if (solflag & 1) {
			printf(" (E:%6.3f N:%6.3f U:%6.3f)", SQRT(Qe[0]), SQRT(Qe[4]), SQRT(Qe[8]));
		}
	}
	else if (soltype == 3) {
		if (len > 0.0) {
			ecef2pos(rb, pos);
			ecef2enu(pos, bl, enu);
			covenu(pos, Qr, Qe);
			pitch = asin(enu[2] / len);
			yaw = atan2(enu[0], enu[1]); if (yaw < 0.0) yaw += 2.0 * PI;
		}
		printf(" P:%12.3f", pitch * R2D);
		printf(" Y:%12.3f", yaw * R2D);
		printf(" L:%12.3f", len);
		if (solflag & 1) {
			printf(" (E:%6.3f N:%6.3f U:%6.3f)", SQRT(Qe[0]), SQRT(Qe[4]), SQRT(Qe[8]));
		}
	}
	if (solflag & 2) {
		printf(" A:%4.1f R:%5.1f N:%2d", sol->age, sol->ratio, sol->ns);
	}
	printf("\n");
}

static void extra_info(const char* path, char* addr, char* port, char* user,
	char* passwd, char* mntpnt, char* str)
{
	char buff[MAXSTRPATH], * p, * q;

	tracet(4, "decodetcpepath: path=%s\n", path);

	if (port) *port = '\0';
	if (user) *user = '\0';
	if (passwd) *passwd = '\0';
	if (mntpnt) *mntpnt = '\0';
	if (str) *str = '\0';

	strcpy(buff, path);

	if (!(p = strrchr(buff, '@'))) p = buff;

	if ((p = strchr(p, '/'))) {
		if ((q = strchr(p + 1, ':'))) {
			*q = '\0'; if (str) sprintf(str, "%.*s", 256 - 1, q + 1);
		}
		*p = '\0'; if (mntpnt) sprintf(mntpnt, "%.255s", p + 1);
		mntpnt[4] = '\0';
	}
	if ((p = strrchr(buff, '@'))) {
		*p++ = '\0';
		if ((q = strchr(buff, ':'))) {
			*q = '\0'; if (passwd) sprintf(passwd, "%.255s", q + 1);
		}
		if (user) sprintf(user, "%.255s", buff);
	}
	else p = buff;

	if ((q = strchr(p, ':'))) {
		*q = '\0'; if (port) sprintf(port, "%.255s", q + 1);
	}
	if (addr) sprintf(addr, "%.255s", p);
}
static int get_site_info(const char* path, prcopt_t* popt, sta_t* sta)
{
	FILE* fp = NULL;
	int i = 0;
	if (!(fp = fopen(path, "r"))) {
		trace(2, "IGSNetwork file open error: %s\n", path);
		return 0;
	}
	char buff[1024], /*temp0[1024],*/* p, * q;
	char name[5] = { '\0' };
	while (fgets(buff, sizeof(buff), fp)) {
		strncpy(name, buff + 15, 3);

		if ((strstr(name, ": {"))) {
			strncpy(name, buff + 5, 4);
			if (!(strstr(name, sta[0].name))) continue;
			i = i;

			////sprintf(p, "%.255s", p + 8);

			while (fgets(buff, sizeof(buff), fp)) {
				char temp0[1024] = { '\0' };
				if (buff[0] == '}')
				{
					fclose(fp);
					return 0;
				}
				if (strstr(buff + 9, "X\":"))
				{
					p = strchr(buff, ':');
					strncpy(temp0, buff + 12, 20);
					temp0[strlen(temp0) - 1] = '\0';
					temp0[strlen(temp0) - 1] = '\0';
					sta[0].pos[0] = atof(temp0);
				}
				if (strstr(buff + 9, "Y\":"))
				{
					p = strchr(buff, ':');
					strncpy(temp0, buff + 12, 20);
					temp0[strlen(temp0) - 1] = '\0';
					temp0[strlen(temp0) - 1] = '\0';
					sta[0].pos[1] = atof(temp0);
				}
				if (strstr(buff + 9, "Z\":"))
				{
					p = strchr(buff, ':');
					strncpy(temp0, buff + 12, 20);
					temp0[strlen(temp0) - 1] = '\0';
					temp0[strlen(temp0) - 1] = '\0';
					sta[0].pos[2] = atof(temp0);
				}

				if (strstr(buff + 9, "Receiver\":"))
				{
					fgets(buff, sizeof(buff), fp);
					strncpy(temp0, buff + 21, 20);
					int len_str = strlen(temp0);
					temp0[len_str - 1] = temp0[len_str - 2] = temp0[len_str - 3] = '\0';
					strncpy(sta[0].rectype, temp0, strlen(temp0));
				}
				if (strstr(buff + 9, "Antenna\":"))
				{
					fgets(buff, sizeof(buff), fp);
					strncpy(temp0, buff + 21, 20);
					int len_str = strlen(temp0);
					temp0[len_str - 1] = temp0[len_str - 2] = temp0[len_str - 3] = '\0';
					strncpy(sta[0].antdes, temp0, strlen(temp0));

					while (fgets(buff, sizeof(buff), fp)) {
						if (strstr(buff + 13, "MarkerUp\":"))
						{
							strncpy(temp0, buff + 25, 20);
							int len_str = strlen(temp0);
							temp0[len_str - 1] = temp0[len_str - 2] = temp0[len_str - 3] = '\0';
							sta[0].del[2] = atof(temp0);
						}
						if (strstr(buff + 13, "MarkerNorth\":"))
						{
							strncpy(temp0, buff + 28, 20);
							int len_str = strlen(temp0);
							temp0[len_str - 1] = temp0[len_str - 2] = temp0[len_str - 3] = '\0';
							sta[0].del[1] = atof(temp0);
						}
						if (strstr(buff + 13, "MarkerEast\":"))
						{
							strncpy(temp0, buff + 27, 20);
							int len_str = strlen(temp0);
							temp0[len_str - 1] = temp0[len_str - 2] = temp0[len_str - 3] = '\0';
							sta[0].del[0] = atof(temp0);
							if (norm(sta[0].del, 3) > 0)sta[0].deltype = 1;
							fclose(fp);
							return 1;
						}
					}
				}

				//if ((strstr(buff, "    },")))break;
				//return 1;
			}
		}
	}

	sta[0].pos[0];

	fclose(fp);
	return 0;
}
static void matchoutfile(prcopt_t* popt, solopt_t* sopt, filopt_t* fopt,char* outfile,char* tracefile,char* statfile) {
	popt->ts = timeget();
	matchout(popt, popt->prcdir, fopt, sopt);
	/* open debug trace */
	if (sopt->trace <= 5 && sopt->trace >= 1) {
		if (*outfile) {
			strcpy(tracefile, outfile);
			strcat(tracefile, ".trace");
		}
		else {
			strcpy(tracefile, fopt->trace);
		}
		traceclose();
		traceopen(tracefile);
		tracelevel(sopt->trace);
	}
	/* open solution statistics */
	if (sopt->sstat > 0) {
		strcpy(statfile, fopt->solf);
		strcat(statfile, ".stat");
		rtkclosestat();
		rtkopenstat(statfile, sopt->sstat);
	}
	if (sopt->ambres && (popt->modear == ARMODE_PPPAR_ILS || popt->modear == ARMODE_CONT)) {
		rtkopenfcbstat(fopt->wl_amb, fopt->nl_amb, fopt->lc_amb);
	}

}
int main(int argc, char** argv)
{
	///*int argc = 3;*/ char** argv = NULL;

	int i, port = 0, outstat = 1, trace = 0, sock = 0;
	char* infile[16] = { 0 }, * outfile = { (char *)"" }, * p;
	int n = 2;
	char tracefile[1024]={0}, statfile[1024] = { 0 }, path[1024] = { 0 }, * ext, outfiletm[1024] = {0};
	gtime_t time = timeget();
	pcvs_t pcvr = { 0 }, pcvs = { 0 };
	pcv_t* pcv = { 0 };
	sta_t stas[MAXRCV] = { {0} };      /* station information */

//	int argc = 8;
//	char* argv[] = { "-C","H:/winmove/project/pppar/pppar_main/conf/PPP/ppp_mgex_2021_rt.conf",
//		"-S", "GREC", "-M", "PPP-KINE", "-A", "0", "L", "0" };
	long t1 = clock();
	prcopt_t popt = prcopt_default/*, popt_*/;
	prcopt_t* popt_;
	solopt_t sopt = solopt_default, sopt_;
	filopt_t fopt = { 0 }, fopt_;
	int nsta;
	port = 0;

	/* initialize rtk server and monitor port */
	rtksvrinit(&svr);
	strinit(&moni);

	/* load options file */
	if (!parsecmd(argc, argv, &popt, &sopt, &fopt, &port)) return 0;

	/*static file read*/
	/*set name*/
	extra_info(fopt.strpath[0], NULL, NULL, NULL, NULL, stas[0].name, NULL);
	//sprintf(stas[0].name, "BRUX");
	soltype = fopt.strfmt[7];
	get_site_info(fopt.rcvantp, &popt, stas);
	/*read atx file*/
	if (popt.mode != PMODE_SINGLE && *fopt.atx) {
		if (*fopt.atx && !(readpcv(fopt.atx, &pcvs))) {
			showmsg("error : no sat ant pcv in %s", fopt.atx);
			//return 0;
		}
		setpcv(timeget(), &popt, &svr.nav, &pcvs, &pcvr, stas);
	}
	if (norm(popt.rb, 3) == 0) {
		popt.rb[0] = stas[0].pos[0];
		popt.rb[1] = stas[0].pos[1];
		popt.rb[2] = stas[0].pos[2];
	}

	//read blq
	if (popt.mode > PMODE_SINGLE && *fopt.blq) {
		//Ã·«∞≤‚’æ√˚
		//fopt.strpath[0];
		//*stas[0].name = *popt.site_list;
		readotl(&popt, fopt.blq, stas);
	}
	/* read erp data */
	if (*fopt.eop) {
		free(svr.nav.erp.data); svr.nav.erp.data = NULL; svr.nav.erp.n = svr.nav.erp.nmax = 0;
		reppath(fopt.eop, fopt.eop, timeget(), "", "");
		if (!readerp(fopt.eop, &svr.nav.erp)) {
			showmsg("error : no erp data %s", fopt.eop);
		}
	}

	matchoutfile(&popt, &sopt, &fopt,outfile,tracefile,statfile);
	static sta_t sta[MAXRCV] = { {""} };
	double pos[3], npos[3] = { 0 };
	char s1[18][MAXRCVCMD] = { "","","","","","","","","","","","","","","","","" }, * cmds[] = { NULL,NULL,NULL,NULL,NULL,NULL ,NULL,NULL,NULL,NULL,NULL ,NULL,NULL,NULL ,NULL,NULL,NULL };
	char s2[18][MAXRCVCMD] = { "","","" ,"","","","","","","","","","","" }, * cmds_periodic[] = { NULL,NULL,NULL ,NULL,NULL,NULL ,NULL,NULL,NULL ,NULL,NULL ,NULL,NULL,NULL ,NULL,NULL,NULL };
	char* ropts[] = { "","","" ,"","","" ,"","","","","","","","","","","" };
	char* paths[] = { "","","" ,"","","" ,"","","","","","","","","","","" };
	char errmsg[2048] = "";
	int  ret, stropt[8] = { 0 };

	/* set stream options */
	stropt[0] = timeout;
	stropt[1] = reconnect;
	stropt[2] = 1000;
	stropt[3] = buffsize;
	stropt[4] = fswapmargin;
	strsetopt(stropt);

	/* set ftp/http directory and proxy */
	//strsetdir(filopt.tempdir);
	strsetproxy(proxyaddr);

	soltype = sopt.posf;

	//paths[4] = paths[3];
	for (i = 0; i < 11; i++)
	{
		paths[i] = fopt.strpath[i];
		strtype[i] = fopt.strtype[i];
		strfmt[i] = fopt.strfmt[i];
	}
	if (fopt.strtype[6] == 2)paths[6] = fopt.solf;
	if (fopt.strtype[8] == 2)
	{		
		strcpy(tracefile, fopt.solf);
		strcat(tracefile, ".ssr");
		paths[8] = tracefile;
	}
	solopt[0] = sopt;
	solopt[1] = sopt;
	solopt[0].posf = fopt.strfmt[6];
	solopt[1].posf = fopt.strfmt[7];

	/* start rtk server */
	if (!rtksvrstart(&svr, svrcycle, buffsize, strtype, paths, strfmt, navmsgsel,
		cmds, cmds_periodic, ropts, nmeacycle, nmeareq, npos, &popt,
		solopt, &moni, errmsg)) {
		traceclose();
		return 0;
	}
	return 0;
}

void TimerTimer(char** in_str) {
	// TODO: Add your implementation code here.
	static int n = 0, inactive = 0;
	sol_t* sol;
	int i, update = 0;

	trace(4, "TimerTimer\n");

	rtksvrlock(&svr);
	//in_str += svr.solbuf->time.sec;
	//in_str += svr.rtk->rb[0].c_str();

	svr.nsol = 0;
	int SolCurrentStat = svr.state ? svr.rtk.sol.stat : 0;

	rtksvrunlock(&svr);

	gtime_t time = timeget();

	//in_str += xxx.c_str() + ' ';
	//vt_t* vt = NULL;

	//prtime(vt, time);
	prstream();
	prsolution(&svr.rtk.sol, svr.rtk.rb);
	//prstream(vt);
	//printf("\nRover 1£∫\n");
	//prsolution(NULL, &svr.nrtk[svr.nrtk[1].thread_num[0]].sol, &svr.nrtk[svr.nrtk[1].thread_num[0]].opt.ru);
	//printf("\nRover 2£∫\n");
	//prsolution(NULL, &svr.nrtk[svr.nrtk[2].thread_num[0]].sol, &svr.nrtk[svr.nrtk[2].thread_num[0]].opt.ru);

	//prsatellite(NULL,2);
	//probserv(vt, 2);
	//prnavidata(NULL);

	// keep alive for monitor port
	//if (!(++n % (1000 / 100)) && 52001) {
	//	unsigned char buf[1];
	//	buf[0] = '\r';
	//	strwrite(&moni, buf, 1);
	//}
};