
static char *identifier = "FeWadsleyite_slb_em.emml:9f0377218f0e21090036fa7c395cbfa56d93d020:Thu Feb 10 16:52:21 2022";



#include <math.h>
#include "coder_error.h"

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static double chebvalat(double x) {
    double c[17] = {
        2.707737068327440945 / 2.0, 0.340068135211091751, -0.12945150184440869e-01, 0.7963755380173816e-03,
        -0.546360009590824e-04, 0.39243019598805e-05, -0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08,
        0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, -0.607797e-13, 0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17
    };
    double x2 = 2 * x;
    double c0 = c[17-2];
    double c1 = c[17-1];
    for (int i=3; i<18; i++) {
        double tmp = c0;
        c0 = c[17-i] - c1;
        c1 = tmp + c1 * x2;
    }
    return c0 + c1 * x;
}

static double Debye(double x) {
    //
    // returns D_3(x) = 3/x^3\int_0^x t^3/(e^t - 1) dt
    //
    
    double val_infinity = 19.4818182068004875;
    double sqrt_eps = sqrt(DBL_EPSILON);
    double log_eps = log(DBL_EPSILON);
    double xcut = -log_eps;

    //Check for negative x (was returning zero)
    assert(x >= 0.);

    if (x < (2.0*sqrt(2.0)*sqrt_eps)) return 1.0 - 3.0*x/8.0 + x*x/20.0;
    else if (x <= 4.0) {
        double t = x*x/8.0 - 1.0;
        double c = chebvalat(t);
        return c - 0.375*x;
    } else if (x < -(log(2.0)+log_eps)) {
        int nexp = (int)(floor(xcut / x));
        double ex = exp(-x);
        double xk = nexp * x;
        double rk = nexp;
        double sum = 0.0;
        for (int i=nexp; i>0; i--) {
            double xk_inv = 1.0/xk;
            sum *= ex;
            sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0)/rk;
            rk -= 1.0;
            xk -= x;
        }
        return val_infinity / (x * x * x) - 3.0 * sum * ex;
    } else if (x < xcut) {
        double x3 = x*x*x;
        double sum = 6.0 + 6.0*x + 3.0*x*x + x3;
        return (val_infinity - 3.0*sum*exp(-x))/x3;
    } else return ((val_infinity/x)/x)/x;
}

double born_B(double t, double p);
double born_Q(double t, double p);
double born_N(double t, double p);
double born_U(double t, double p);
double born_Y(double t, double p);
double born_X(double t, double p);
double born_dUdT(double t, double p);
double born_dUdP(double t, double p);
double born_dNdT(double t, double p);
double born_dNdP(double t, double p);
double born_dXdT(double t, double p);
double gSolvent(double t, double p);
double DgSolventDt(double t, double p);
double DgSolventDp(double t, double p);
double D2gSolventDt2(double t, double p);
double D2gSolventDtDp(double t, double p);
double D2gSolventDp2(double t, double p);
double D3gSolventDt3(double t, double p);
double D3gSolventDt2Dp(double t, double p);
double D3gSolventDtDp2(double t, double p);
double D3gSolventDp3(double t, double p);
double D4gSolventDt4(double t, double p);

static double coder_a(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = sqrt(13.447928624308624*x1 - 5.1532061946641798*x2 - 3.3598705850000004);
    double x4 = 2.2181640000000002*x3;
    double x5 = 665.44920000000002*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 26.76322271834383*T - 32433372.848412335*x1 + 29087789.962084591*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 5442392.9838575013 + 24008015.650664754/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-5.1532061946641798*pow(x0, 4.0/3.0) + 13.447928624308624*pow(x0, 2.0/3.0) - 3.3598705850000004);
    double x2 = x1/T;
    double x3 = 665.44920000000002*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -116189.90245127954*x2*x5/x6 + 38729.96748375985*x2*(-0.0045082329349858709*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(13.447928624308624*x1 - 5.1532061946641798*x3 - 3.3598705850000004);
    double x6 = 2.2181640000000002*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-4.4826428747695406*x2 + 3.4354707964427864*x4);
    double x10 = 665.44920000000002*x5/T;
    double x11 = exp(-x10);

    result += 116189.90245127954*x11*x9/(-x11 + 1) + 21622248.565608222*x2 - 38783719.949446119*x4 - 116189.90245127956*x7*x9/(-x7 + 1) + 38729.967483759858*x9*(-1.3524698804957613*x8*Debye(x6) + 3/(exp(x6) - 1)) - 38729.96748375985*x9*(-0.0045082329349858709*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 48016031.301329508/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(13.447928624308624*x2 - 5.1532061946641798*x3 - 3.3598705850000004);
    double x5 = x0*x4;
    double x6 = 665.44920000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1039773368.5660272*x2 + 398438057.90698588*x3 + 259780078.68040454);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1330.8984*x5)/((x8)*(x8)) - 38729.96748375985*x4*(x0*(0.013524698804957613*T*x11 - 9.0/x13) + 0.0045082329349858709*x11 - 1996.3476000000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 13.447928624308624*x2;
    double x4 = 5.1532061946641798*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 3.3598705850000004;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 665.44920000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 77318477.634282008*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(3.4354707964427864*x2 - 4.4826428747695406)*(-38729.96748375985*x6*(1996.3476000000001*x0*x13*x14/((x15)*(x15)) - 0.0045082329349858709*x12/pow(x5, 3.0/2.0) + (0.013524698804957613*x12*x13 - 9.0/x15)/(-x3 + x4 + 3.3598705850000004)) + x11*x9/x10 + x11*exp(-1330.8984*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 13.447928624308624*x2;
    double x5 = 5.1532061946641798*x3;
    double x6 = x4 - x5 - 3.3598705850000004;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(8.0160985250331684*x2 - 7.4710714579492343);
    double x10 = 2.2181640000000002*x7;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = x11/x12;
    double x14 = 116189.90245127956*x13;
    double x15 = 1.0/(-x4 + x5 + 3.3598705850000004);
    double x16 = x3*((3.4354707964427864*x2 - 4.4826428747695406)*(3.4354707964427864*x2 - 4.4826428747695406));
    double x17 = x15*x16;
    double x18 = 257728.25878094009*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 1.0/T;
    double x22 = x21*x7;
    double x23 = 665.44920000000002*x22;
    double x24 = exp(-x23);
    double x25 = -x24 + 1;
    double x26 = x24/x25;
    double x27 = 116189.90245127954*x26;
    double x28 = 77318477.634282008*x17*x21;
    double x29 = exp(x10);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x10);
    double x33 = x32*x8;
    double x34 = -116189.90245127957*x31 + 52381.114494365414*x33;
    double x35 = exp(x23);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x23);
    double x39 = x38*x8;
    double x40 = -116189.90245127954*x37 + 174.60371498121802*x39;
    double x41 = x16*x8;

    result += x0*(144048093.90398854*x0 - x13*x18 + x14*x20 + x14*x9 - 36037080.942680374*x2 - x20*x27 + x20*x34 - x20*x40 + x26*x28 - x27*x9 + 90495346.548707604*x3 + x34*x9 - x40*x9 - 38729.967483759858*x41*(x15*(-9.0000000000000018*x31 + 4.0574096414872844*x33) - 1.3524698804957613*x19*x32 + 6.6544920000000012*x29*x8/((x30)*(x30))) + 38729.96748375985*x41*(x15*(-9.0*x37 + 0.013524698804957613*x39) - 0.0045082329349858709*x19*x38 + 1996.3476000000001*x21*x35*x8/((x36)*(x36))) + x28*exp(-1330.8984*x22)/((x25)*(x25)) - x18*exp(-4.4363280000000005*x7)/((x12)*(x12)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-3119320105.6980815*x2 + 1195314173.7209578*x3 + 779340236.04121363);
    double x5 = 1.0/T;
    double x6 = 13.447928624308624*x2;
    double x7 = 5.1532061946641798*x3;
    double x8 = x6 - x7 - 3.3598705850000004;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 665.44920000000002*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1330.8984*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -1996.3476000000001*x0*x22*x9;
    double x24 = 0.013524698804957613*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 38729.96748375985*x9;
    double x27 = x17*(-x6 + x7 + 3.3598705850000004);

    result += x0*(-154354557260.85257*x15*x18 - x15*x4 - 51451519086.950859*x16*x18 - x16*x4 + x26*(0.0045082329349858709*x19 + x23 + x25) - x26*(-1328467.9133419201*x22*x27 + x23 + x24 + 3.0*x25 + 2656935.8266838402*x27*exp(x14)/((x21)*(x21)*(x21))) - 102903038173.90172*x18*exp(-1996.3476000000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(531250743.87598115*x2 - 693182245.71068466);
    double x5 = 13.447928624308624*x2;
    double x6 = 5.1532061946641798*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 3.3598705850000004;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 665.44920000000002*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1330.8984*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 3.4354707964427864*x2 - 4.4826428747695406;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0045082329349858709*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 1996.3476000000001*x3;
    double x26 = 0.013524698804957613*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 38729.96748375985*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 3.3598705850000004);

    result += x0*x1*x2*(154354557260.85257*x14*x18 - x14*x4 + 51451519086.950859*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(5989.0428000000002*x0*x31 - x26*x30 + 3.0*x27*x32) + 1328467.9133419201*x16*x24 - 2656935.8266838402*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 102903038173.90172*x18*exp(-1996.3476000000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(619792534.52197802*x2 - 577651871.42557049);
    double x4 = 13.447928624308624*x2;
    double x5 = 5.1532061946641798*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 3.3598705850000004;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 665.44920000000002*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1330.8984*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 3.4354707964427864*x2 - 4.4826428747695406;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0045082329349858709*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 1996.3476000000001*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 3.3598705850000004;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.013524698804957613*x33));
    double x35 = 8.0160985250331684*x2 - 7.4710714579492343;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0045082329349858709*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(154354557260.85257*x13*x20 + x13*x3 + 51451519086.950859*x14*x20 + x14*x3 + 38729.96748375985*x15*x34 + 38729.96748375985*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(10.30641238932836*x2 - 13.447928624308622)/pow(x6, 5.0/2.0) + x24*x35 - 1328467.9133419201*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(6.8709415928855728*x2 - 8.9652857495390812) + 2656935.8266838402*x37*exp(x12)/((x26)*(x26)*(x26))) + 102903038173.90172*x20*exp(-1996.3476000000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 13.447928624308624*x2;
    double x5 = 5.1532061946641798*x3;
    double x6 = x4 - x5 - 3.3598705850000004;
    double x7 = sqrt(x6);
    double x8 = 2.2181640000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 3.4354707964427864*x2 - 4.4826428747695406;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.4363280000000005*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 3.3598705850000004;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 257728.25878094009*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = x2*(26.720328416777228*x2 - 19.922857221197958);
    double x24 = x22*x23;
    double x25 = 116189.90245127956*x21;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 665.44920000000002*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 116189.90245127954*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1330.8984*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 77318477.634282008*x26;
    double x38 = x19*x37;
    double x39 = x18*(6.8709415928855728*x2 - 8.9652857495390812);
    double x40 = ((x12)*(x12));
    double x41 = x0*x40;
    double x42 = x39*x41;
    double x43 = 257728.25878094009*x42;
    double x44 = (10.30641238932836*x2 - 13.447928624308622)/pow(x6, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 1.0/x17;
    double x47 = 8.0160985250331684*x2 - 7.4710714579492343;
    double x48 = x46*x47;
    double x49 = x12*x3;
    double x50 = 257728.25878094009*x49;
    double x51 = x16*x50;
    double x52 = x46*(16.032197050066337*x2 - 14.942142915898469);
    double x53 = x47*x49;
    double x54 = x11*x53;
    double x55 = x21*x50;
    double x56 = x37*x42;
    double x57 = x37*x49;
    double x58 = x36*x57;
    double x59 = x31*x57;
    double x60 = exp(x8);
    double x61 = x60 - 1;
    double x62 = 1.0/x61;
    double x63 = Debye(x8);
    double x64 = x22*x63;
    double x65 = -3*x62 + 1.3524698804957613*x64;
    double x66 = 38729.967483759858*x22;
    double x67 = exp(x28);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x28);
    double x71 = T*x22*x70;
    double x72 = -3*x69 + 0.0045082329349858709*x71;
    double x73 = 38729.96748375985*x22;
    double x74 = 1.3524698804957613*x63;
    double x75 = x11*x74;
    double x76 = x60/((x61)*(x61));
    double x77 = 6.6544920000000012*x76;
    double x78 = x22*x77;
    double x79 = x46*(-9.0000000000000018*x62 + 4.0574096414872844*x64) - x75 + x78;
    double x80 = 77459.934967519715*x79;
    double x81 = x22*x53;
    double x82 = 0.0045082329349858709*T*x70;
    double x83 = x11*x82;
    double x84 = x67/((x68)*(x68));
    double x85 = 1996.3476000000001*x26*x84;
    double x86 = x22*x85;
    double x87 = x46*(-9.0*x69 + 0.013524698804957613*x71) - x83 + x86;
    double x88 = 77459.934967519701*x87;
    double x89 = x12*x2;
    double x90 = x44*x89;
    double x91 = x2*x40;
    double x92 = x46*x91;
    double x93 = x11*x91;
    double x94 = 3.0000000000000004*x65;
    double x95 = x18*x91;
    double x96 = x39*x89;
    double x97 = x33*x92;
    double x98 = 3.0*x72;

    result += (-576192375.61595416*x0 - 1715050.6362316958*x14*x16 - 571683.54541056522*x14*x21 + x14*x80 - x14*x88 - x16*x20 - x16*x43 + 96098882.51381433*x2 - x20*x21 - x21*x43 - 348569.70735383866*x21*x54 - x23*x65*x66 + x23*x72*x73 - x24*x25 + x24*x32 - x25*x45 - 301651155.16235864*x3 + 51451519086.950859*x31*x34 + x31*x38 + 348569.7073538386*x31*x54 + x31*x56 + x32*x45 + 154354557260.85257*x34*x36 + x36*x38 + x36*x56 - 38729.967483759858*x45*x65 + 38729.96748375985*x45*x72 + x48*x51 + x48*x55 - x48*x58 - x48*x59 - x49*x66*(x47*x75 - x47*x78 - x48*x94 + x74*x90 - 14.760754592688004*x76*x92 - x77*x93 + 3.0000000000000004*x79*x92 + x94*x95 + x94*x96 + 29.521509185376008*x92*exp(x15)/((x61)*(x61)*(x61))) + x49*x73*(x47*x83 - x47*x86 - x48*x98 + x82*x90 - 1328467.9133419201*x84*x97 - x85*x93 + 3.0*x87*x92 + x95*x98 + x96*x98 + 2656935.8266838402*x97*exp(x35)/((x68)*(x68)*(x68))) + x51*x52 + x52*x55 - x52*x58 - x52*x59 - 116189.90245127957*x54*x65 + 116189.90245127954*x54*x72 + x80*x81 - x81*x88 + 102903038173.90172*x34*exp(-1996.3476000000001*x27)/((x30)*(x30)*(x30)) - 1143367.0908211304*x14*exp(-6.6544920000000012*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*4.28;
static const double Vmax = 1.15*4.28;
static double V = 0.9*4.28;

static void coder_solve_V(double T, double P) {
    // Newtsafe routine, newton + bisection with bracket check
    // Update if *either* T or P changes (was only if both changed)
    if ((T != Told) || (P != Pold)) {
        // check bracket
        double fa = -coder_dadv(T, Vmin) - P;
        double fb = -coder_dadv(T, Vmax) - P;
        if ( isnan(fa) ) {
            fprintf(stderr, "Error: lower bracket isnan for Vmin=%g\n",Vmin);
            coder_error_flag = CODER_ERR_NAN;
            return;
        }
        if ( isnan(fb) ) {
            fprintf(stderr, "Error: upper bracket isnan for Vmax=%g\n",Vmax);
            coder_error_flag = CODER_ERR_NAN;
            return;
        }
        if ( fa*fb > 0.) 
        {
            fprintf(stderr, "Error: improper  initial bracket in solve_V\n");
            coder_error_flag = CODER_ERR_BOUNDS;
            return;
        }
        //set up bisection parameters
        double a = Vmin;
        double b = Vmax;
        double c = 0.;
        double fc = 0.;
        Told = T;
        Pold = P;
        //start Newton step
        double f = 0.0;
        int iter = 0;
        do {
            //newton step
            f = -coder_dadv(T, V) - P;
            double df = -coder_d2adv2(T, V);
            V -= f/df;
            if (V < a || V > b || df == 0. || isnan(V) ) {
                //take a bisection step
                c = 0.5*(a + b);
                fc = -coder_dadv(T, c) - P;
                if ( fa*fc < 0 ) { // left bracket
                    b = c;
                    fb = fc;
                }
                else if ( fb*fc < 0 ) { //right bracket
                    a = c;
                    fa = fc;
                }
                //reset V to middle of new bracket and clear f
                V = 0.5*( a + b);
                f = 1.;
            }
            iter++;
        } while ((fabs(f) > tol) && (iter < MAX_ITS));
        if ( iter == MAX_ITS) {
            fprintf(stderr, "Error: max iterations exceeded in solve_V\n");
            coder_error_flag = CODER_ERR_MAXITS;
            return;
        }
    }
}

static double coder_g(double T, double P) {
    coder_solve_V(T, P);
    double A = coder_a(T, V);
    return A + P*V;
}

static double coder_dgdt(double T, double P) {
    coder_solve_V(T, P);
    double dAdT = coder_dadt(T, V);
    return dAdT;
}

static double coder_dgdp(double T, double P) {
    coder_solve_V(T, P);
    return V;
}

static double coder_d2gdt2(double T, double P) {
    coder_solve_V(T, P);
    double d2AdT2 = coder_d2adt2(T, V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    return d2AdT2 - d2AdTdV*d2AdTdV/d2AdV2;
}

static double coder_d2gdtdp(double T, double P) {
    coder_solve_V(T, P);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    return - d2AdTdV/d2AdV2;
}

static double coder_d2gdp2(double T, double P) {
    coder_solve_V(T, P);
    double d2AdV2 = coder_d2adv2(T, V);
    return - 1.0/d2AdV2;
}

static double coder_d3gdt3(double T, double P) {
    coder_solve_V(T, P);
    double d3AdT3 = coder_d3adt3(T, V);
    double d3AdT2dV = coder_d3adt2dv(T, V);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    double d2VdT2 = (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)/d2AdV2;
    return d3AdT3 + 2.0*d3AdT2dV*dVdT + d3AdTdV2*dVdT*dVdT + d2AdTdV*d2VdT2;
}

static double coder_d3gdt2dp(double T, double P) {
    coder_solve_V(T, P);
    double d3AdT2dV = coder_d3adt2dv(T, V);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdP = -1.0/d2AdV2;
    double d2VdTdP = (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)/d2AdV2;
    return d3AdT2dV*dVdP + d3AdTdV2*dVdT*dVdP + d2AdTdV*d2VdTdP;
}

static double coder_d3gdtdp2(double T, double P) {
    coder_solve_V(T, P);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    return (d3AdTdV2 + d3AdV3*dVdT)/d2AdV2/d2AdV2;
}

static double coder_d3gdp3(double T, double P) {
    coder_solve_V(T, P);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdV2 = coder_d2adv2(T, V);
    double dVdP = -1.0/d2AdV2;
    return d3AdV3*dVdP/d2AdV2/d2AdV2;
}

static double coder_s(double T, double P) {
    double result = -coder_dgdt(T, P);
    return result;
}

static double coder_v(double T, double P) {
    double result = coder_dgdp(T, P);
    return result;
}

static double coder_cv(double T, double P) {
    double result = -T*coder_d2gdt2(T, P);
    double dvdt = coder_d2gdtdp(T, P);
    double dvdp = coder_d2gdp2(T, P);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P) {
    double result = -T*coder_d2gdt2(T, P);
    return result;
}

static double coder_dcpdt(double T, double P) {
    double result = -T*coder_d3gdt3(T, P) - coder_d2gdt2(T, P);
    return result;
}

static double coder_alpha(double T, double P) {
    double result = coder_d2gdtdp(T, P)/coder_dgdp(T, P);
    return result;
}

static double coder_beta(double T, double P) {
    double result = -coder_d2gdp2(T, P)/coder_dgdp(T, P);
    return result;
}

static double coder_K(double T, double P) {
    double result = -coder_dgdp(T, P)/coder_d2gdp2(T, P);
    return result;
}

static double coder_Kp(double T, double P) {
    double result = coder_dgdp(T, P);
    result *= coder_d3gdp3(T, P);
    result /= pow(coder_d2gdp2(T, P), 2.0);
    return result - 1.0;
}


#include <math.h>

static double coder_dparam_a(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_dadt(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_dadv(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2adt2(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2adtdv(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2adv2(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3adt3(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3adt2dv(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3adtdv2(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3adv3(double T, double V, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}


static double coder_dparam_g(double T, double P, int index) {
    coder_solve_V(T, P);
    double dAdz = coder_dparam_a(T, V, index);
    return dAdz + P*V;
}

static double coder_dparam_dgdt(double T, double P, int index) {
    coder_solve_V(T, P);
    double dAdTdz = coder_dparam_dadt(T, V, index);
    return dAdTdz;
}

static double coder_dparam_dgdp(double T, double P, int index) {
    coder_solve_V(T, P);
    return 0.0; /* V; */
}

static double coder_dparam_d2gdt2(double T, double P, int index) {
    coder_solve_V(T, P);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d2AdT2dz = coder_dparam_d2adt2(T, V, index);
    double d2AdTdVdz = coder_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    /* return d2AdT2 - d2AdTdV*d2AdTdV/d2AdV2; */
    return d2AdT2dz - 2.0*d2AdTdV*d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
}

static double coder_dparam_d2gdtdp(double T, double P, int index) {
    coder_solve_V(T, P);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d2AdTdVdz = coder_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    /* return - d2AdTdV/d2AdV2; */
    return - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
}

static double coder_dparam_d2gdp2(double T, double P, int index) {
    coder_solve_V(T, P);
    double d2AdV2 = coder_d2adv2(T, V);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    /* return - 1.0/d2AdV2; */
    return d2AdV2dz/d2AdV2/d2AdV2;
}

static double coder_dparam_d3gdt3(double T, double P, int index) {
    coder_solve_V(T, P);
    double d3AdT2dV = coder_d3adt2dv(T, V);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d3AdT3dz = coder_dparam_d3adt3(T, V, index);
    double d3AdT2dVdz = coder_dparam_d3adt2dv(T, V, index);
    double d3AdTdV2dz = coder_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = coder_dparam_d3adv3(T, V, index);
    double d2AdTdVdz = coder_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    double dVdT = - d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    double d2VdT2 = (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)/d2AdV2;
    double d2VdT2dz = (-d3AdT2dVdz - 2.0*d3AdTdV2dz*dVdT - 2.0*d3AdTdV2*dVdTdz
                        - d3AdV3dz*dVdT*dVdT - 2.0*d3AdV3*dVdT*dVdTdz)/d2AdV2
                    - (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)*d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdT3 + 2.0*d3AdT2dV*dVdT + d3AdTdV2*dVdT*dVdT + d2AdTdV*d2VdT2; */
    return d3AdT3dz + 2.0*d3AdT2dVdz*dVdT + 2.0*d3AdT2dV*dVdTdz + d3AdTdV2dz*dVdT*dVdT
            + 2.0*d3AdTdV2*dVdT*dVdTdz + d2AdTdVdz*d2VdT2 + d2AdTdV*d2VdT2dz;
}

static double coder_dparam_d3gdt2dp(double T, double P, int index) {
    coder_solve_V(T, P);
    double d3AdT2dV = coder_d3adt2dv(T, V);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d3AdT2dVdz = coder_dparam_d3adt2dv(T, V, index);
    double d3AdTdV2dz = coder_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = coder_dparam_d3adv3(T,V, index);
    double d2AdTdVdz = coder_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    double dVdP = -1.0/d2AdV2;
    double dVdPdz = d2AdV2dz/d2AdV2/d2AdV2;
    double d2VdTdP = (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)/d2AdV2;
    double d2VdTdPdz = (-d3AdTdV2dz*dVdP -d3AdTdV2*dVdPdz
            - d3AdV3dz*dVdT*dVdP - d3AdV3*dVdTdz*dVdP - d3AdV3*dVdT*dVdPdz)/d2AdV2
            - (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)*d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdT2dV*dVdP + d3AdTdV2*dVdT*dVdP + d2AdTdV*d2VdTdP; */
    return d3AdT2dVdz*dVdP + d3AdT2dV*dVdPdz + d3AdTdV2dz*dVdT*dVdP + d3AdTdV2*dVdTdz*dVdP
        + d3AdTdV2*dVdT*dVdPdz + d2AdTdVdz*d2VdTdP + d2AdTdV*d2VdTdPdz;
}

static double coder_dparam_d3gdtdp2(double T, double P, int index) {
    coder_solve_V(T, P);
    double d3AdTdV2 = coder_d3adtdv2(T, V);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdTdV = coder_d2adtdv(T, V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d3AdTdV2dz = coder_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = coder_dparam_d3adv3(T, V, index);
    double d2AdTdVdz = coder_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    /* return (d3AdTdV2 + d3AdV3*dVdT)/d2AdV2/d2AdV2; */
    return (d3AdTdV2dz + d3AdV3dz*dVdT + d3AdV3*dVdTdz)/d2AdV2/d2AdV2
        - 2.0*(d3AdTdV2 + d3AdV3*dVdT)*d2AdV2dz/d2AdV2/d2AdV2/d2AdV2;
}

static double coder_dparam_d3gdp3(double T, double P, int index) {
    coder_solve_V(T, P);
    double d3AdV3 = coder_d3adv3(T,V);
    double d2AdV2 = coder_d2adv2(T, V);
    double d3AdV3dz = coder_dparam_d3adv3(T, V, index);
    double d2AdV2dz = coder_dparam_d2adv2(T, V, index);
    double dVdP = -1.0/d2AdV2;
    double dVdPdz = d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdV3*dVdP/d2AdV2/d2AdV2; */
    return d3AdV3dz*dVdP/d2AdV2/d2AdV2 + d3AdV3*dVdPdz/d2AdV2/d2AdV2
        - 2.0*d3AdV3*dVdP*d2AdV2dz/d2AdV2/d2AdV2/d2AdV2;
}

static int coder_get_param_number(void) {
    return 0;
}

static const char *paramNames[0] = {  };

static const char *paramUnits[0] = {  };

static const char **coder_get_param_names(void) {
    return paramNames;
}

static const char **coder_get_param_units(void) {
    return paramUnits;
}

static void coder_get_param_values(double **values) {
}

static int coder_set_param_values(double *values) {
    return 1;
}

static double coder_get_param_value(int index) {
    double result = 0.0;
    switch (index) {
     default:
         break;
    }
    return result;
}

static int coder_set_param_value(int index, double value) {
    int result = 1;
    switch (index) {
     default:
         break;
    }
    return result;
}



const char *FeWadsleyite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *FeWadsleyite_slb_em_coder_calib_name(void) {
    return "FeWadsleyite_slb_em";
}

const char *FeWadsleyite_slb_em_coder_calib_formula(void) {
    return "Fe2SiO4";
}

const double FeWadsleyite_slb_em_coder_calib_mw(void) {
    return 203.77710000000002;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,4.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,2.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0
    };

const double *FeWadsleyite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double FeWadsleyite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double FeWadsleyite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double FeWadsleyite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double FeWadsleyite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double FeWadsleyite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double FeWadsleyite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double FeWadsleyite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double FeWadsleyite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double FeWadsleyite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int FeWadsleyite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **FeWadsleyite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **FeWadsleyite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void FeWadsleyite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int FeWadsleyite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double FeWadsleyite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int FeWadsleyite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double FeWadsleyite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

