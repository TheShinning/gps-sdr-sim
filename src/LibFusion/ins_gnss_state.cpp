/*-----------------------------------------------------------------------------
* ins-gnss-state.cc : ins-gnss coupled estimated states interface
* get number and index for each parameter
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/21.
*----------------------------------------------------------------------------*/
#include "rtklib.h"

/* get number of position states---------------------------------------------*/
extern int xnA(const insopt_t* opt) { return 3; }

/* get number of velocity states---------------------------------------------*/
extern int xnV(const insopt_t* opt) { return 3; }

/* get number of attitude states---------------------------------------------*/
extern int xnP(const insopt_t* opt) { return 3; }

/* get number of accl/gyro bias states---------------------------------------*/
extern int xnBg(const insopt_t* opt) { return (opt->est_bg ? 3 : 0); }
extern int xnBa(const insopt_t* opt) { return (opt->est_ba ? 3 : 0); }

/* get number of residual scale factors of accl/gyro------------------------*/
extern int xnSg(const insopt_t* opt) { return (opt->est_sg ? 3 : 0); }
extern int xnSa(const insopt_t* opt) { return (opt->est_sa ? 3 : 0); }

/* get number of lever arm for ins body to ant.------------------------------*/
extern int xnArmGps(const insopt_t* opt) { return (opt->est_armgps ? 3 : 0); }

extern int xnDt(const insopt_t* opt) { return (opt->est_igdt ? 1 : 0); }

/* get number of odometry scale factor---------------------------------------*/
extern int xnKod(const insopt_t* opt) { return (opt->est_kod ? 1 : 0); }

/* get number of GNSS frequency ---------------------------------------------*/
extern int xnFrq(const insopt_t* opt) {
    auto* p = (prcopt_t*)opt->gnss_opt;
    return opt->gnss_opt ? (p->ionoopt == IONOOPT_IFLC ? 1 : p->nf) : 0;
}

/* get number of GNSS clock parameter ----------------------------------------*/
extern int xnClk(const insopt_t* opt) {
    auto* p = (prcopt_t*)opt->gnss_opt;
    int no_clk =
        (p->mode >= PMODE_LC_POS && p->mode <= PMODE_LC_PPP) || p->mode == PMODE_TC_PPK ||
        p->mode == PMODE_TC_DGPS || p->mode == PMODE_STC_PPP || p->mode == PMODE_STC_PPK || p->sdopt == 1;
    return opt->gnss_opt ? (no_clk ? 0 : (p->bd3opt >= BD3OPT_BD2_3 ? NSYS + 1 : NSYS)) : 0;
}
/* get number of receiver clock drift (non-close-loop correction states) ----*/
extern int xnRr(const insopt_t* opt)
{
    return (opt->velopt ? 1 : 0);
}

/* get number of troposphere parameter----------------------------------------*/
extern int xnTrp(const insopt_t* opt) {
    auto* p = (prcopt_t*)opt->gnss_opt;
    int no_trp = (p->mode >= PMODE_INS_MECH && p->mode <= PMODE_LC_PPP) || (p->mode == PMODE_STC_PPP || p->mode == PMODE_STC_PPK);
    if (no_trp) return 0;
    return (opt->gnss_opt ? (p->tropopt < TROPOPT_EST ? 0 : (p->tropopt < TROPOPT_ESTG ? 1 : 2)) : 0) *
        (p->mode == PMODE_TC_PPK ? 2 : 1);
}

/* get number of ionosphere parameter----------------------------------------*/
extern int xnIon(const insopt_t* opt) {
    auto* p = (prcopt_t*)opt->gnss_opt;
    int no_ion = (p->mode >= PMODE_INS_MECH && p->mode <= PMODE_LC_PPP) || (p->mode == PMODE_STC_PPK || p->mode == PMODE_STC_PPP);
    if (no_ion) return 0;
    return opt->gnss_opt ? (p->ionoopt != IONOOPT_UC ? 0 : MAXSAT) : 0;
}

/* get number of ambiguity parameter-----------------------------------------*/
extern int xnAmb(const insopt_t* opt) {
    auto* p = (prcopt_t*)opt->gnss_opt;
    int no_amb = (p->mode >= PMODE_INS_MECH && p->mode <= PMODE_LC_PPP) || (p->mode == PMODE_STC_PPK || p->mode == PMODE_STC_PPP);
    if (no_amb) return 0;
    return opt->gnss_opt ? ((p->mode == PMODE_TC_PPK || p->mode == PMODE_TC_PPP) ? xnFrq(opt) * MAXSAT : 0) : 0;
}

/* get number of all ins parameter ------------------------------------------*/
extern int xnIns(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt)
        + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt);
}

/* get number of close-loop correction states--------------------------------*/
extern int xnClp(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt)
        + xnSa(opt) + xnSg(opt);
}

/* get number of all ins parameter ------------------------------------------*/
extern int xnGnss(const insopt_t* opt) {
    double a = xnClk(opt);
    double b = xnTrp(opt);
    double c = xnRr(opt);
    return xnClk(opt) + xnRr(opt) + xnTrp(opt) + xnIon(opt) + xnAmb(opt);
}

extern int xnR(const insopt_t* opt) {
    if (xnAmb(opt) == 0) return 0;
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt)
        + xnSa(opt) + xnSg(opt)
        + xnKod(opt)
        + xnArmGps(opt)
        + xnClk(opt) + xnTrp(opt) + xnIon(opt) + xnRr(opt);
}

/* get number of all states--------------------------------------------------*/
extern int xnX(const insopt_t* opt) { return xnIns(opt) + xnGnss(opt); }

/* position states index ----------------------------------------------------*/
extern int xiP(const insopt_t* opt) { return 0; }

/* velocity states index ----------------------------------------------------*/
extern int xiV(const insopt_t* opt) { return xnP(opt); }

/* attitude  states index ---------------------------------------------------*/
extern int xiA(const insopt_t* opt) { return xnP(opt) + xnV(opt); }

/* accl bias state index ----------------------------------------------------*/
extern int xiBg(const insopt_t* opt) { return xnP(opt) + xnV(opt) + xnA(opt); }

/* gyro bias state index ----------------------------------------------------*/
extern int xiBa(const insopt_t* opt) { return xnP(opt) + xnV(opt) + xnA(opt) + xnBg(opt); }

/* get index of residual scale factors of accl-------------------------------*/
extern int xiSg(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt) + xnBa(opt) + xnBg(opt);
}

/* get index of residual scale factors of gyroscopes-------------------------*/
extern int xiSa(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt) + xnBa(opt) + xnBg(opt) + xnSg(opt);
}

/* get index of lever arm for ins body to ant.-------------------------------*/
extern int xiArmGps(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt);
}

EXPORT int xiDt(const insopt_t* opt)
{
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt);
}

/* get index of odometry scale factor states---------------------------------*/
extern int xiKod(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt) + xnA(opt) + xnDt(opt);
}

/* get index of GNSS clock parameter ---------------------------------------*/
extern int xiClk(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt);
}

/* get index of troposphere parameter --------------------------------------*/
extern int xiTrp(const insopt_t* opt) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt)
        + xnClk(opt);
}

/* get index of ionosphere parameter ---------------------------------------*/
extern int xiIon(const insopt_t* opt, int sat) {
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt)
        + xnClk(opt) + xnTrp(opt) + sat - 1;
}

/* get index of receiver clock drift-----------------------------------------*/
extern int xiRr(const insopt_t* opt)
{
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt)
        + xnClk(opt) + xnTrp(opt) + xnIon(opt);
}

/* get index of ambiguity parameter ---------------------------------------*/
extern int xiAmb(const insopt_t* opt, int sat, int f) {
    int a = xnP(opt);
    int b = xnV(opt);
    int c = xnA(opt);
    int d = xnBa(opt);
    int e = xnBg(opt);
    int i = xnSa(opt);
    int j = xnSg(opt);
    int k = xnKod(opt);
    int l = xnArmGps(opt);
    int m = xnClk(opt);
    int n = xnTrp(opt);
    int o = xnIon(opt);
    int p = xnRr(opt);
    int q = MAXSAT;
    int cc = a + b + c + d + e + i + j + k + l + m + n + o + p + q * f + sat - 1;
    return xnP(opt) + xnV(opt) + xnA(opt)
        + xnBa(opt) + xnBg(opt) + xnSa(opt) + xnSg(opt)
        + xnArmGps(opt) + xnDt(opt)
        + xnKod(opt)
        + xnClk(opt) + xnTrp(opt) + xnIon(opt) + xnRr(opt)
        + MAXSAT * f + sat - 1;
}


