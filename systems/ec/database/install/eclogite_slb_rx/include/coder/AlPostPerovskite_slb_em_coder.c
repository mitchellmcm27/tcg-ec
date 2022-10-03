
static char *identifier = "AlPostPerovskite_slb_em.emml:933c0ee6afb091d644ed9c4e7bbdeca91a80183b:Thu Feb 10 16:51:27 2022";



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
    double x3 = sqrt(-11.464926975874729*x1 + 18.097006708817457*x2 + 1.7429853801249982);
    double x4 = 2.5406503333333337*x3;
    double x5 = 762.19510000000002*x3/T;

    result += 124.7169392722986*T*log(1 - exp(-x5)) - 41.572313090766201*T*Debye(x5) - 23847239.755976487*x1 + 21282891.132854495*x2 - 37415.08178168958*log(1 - exp(-x4)) + 12471.69392722986*Debye(x4) + 5302558.875;
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(18.097006708817457*pow(x0, 4.0/3.0) - 11.464926975874729*pow(x0, 2.0/3.0) + 1.7429853801249982);
    double x2 = x1/T;
    double x3 = 762.19510000000002*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -95058.64000034357*x2*x5/x6 + 31686.213333447853*x2*(-0.0039360001133568034*T*x4/x1 + 3/(exp(x3) - 1)) - 41.572313090766201*x4 + 124.7169392722986*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-11.464926975874729*x1 + 18.097006708817457*x3 + 1.7429853801249982);
    double x6 = 2.5406503333333337*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(3.8216423252915761*x2 - 12.064671139211637*x4);
    double x10 = 95058.64000034357*x9;
    double x11 = 762.19510000000002*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 15898159.837317657*x2 - 28377188.177139327*x4 + 31686.213333447857*x9*(-1.1808000340070408*x8*Debye(x6) + 3/(exp(x6) - 1)) - 31686.213333447853*x9*(-0.0039360001133568034*T*x8*Debye(x11) + 3/(exp(x11) - 1));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-11.464926975874729*x2 + 18.097006708817457*x3 + 1.7429853801249982);
    double x5 = x0*x4;
    double x6 = 762.19510000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-830670986.77019894*x2 + 1311186582.5253873*x3 + 126284919.97211327);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1524.3902*x5)/((x8)*(x8)) + 31686.213333447853*x4*(x0*(0.01180800034007041*T*x11 - 9.0/x13) + 0.0039360001133568034*x11 - 2286.5853000000002*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 18.097006708817457*pow(x1, 4.0/3.0) - 11.464926975874729*x2 + 1.7429853801249982;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 762.19510000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 72453229.620925874*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(12.064671139211637*x2 - 3.8216423252915761)*(31686.213333447853*x4*(-2286.5853000000002*x0*x11*x12/((x13)*(x13)) + 0.0039360001133568034*x10/pow(x3, 3.0/2.0) + (0.01180800034007041*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1524.3902*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = 18.097006708817457*pow(x0, 4.0/3.0) - 11.464926975874729*x1 + 1.7429853801249982;
    double x3 = sqrt(x2);
    double x4 = 1.0/x3;
    double x5 = x4*(2675986.2046046592*x1 - 605466.87001661293);
    double x6 = 2.5406503333333337*x3;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = x7/x8;
    double x10 = 1.0/T;
    double x11 = x10*x3;
    double x12 = 762.19510000000002*x11;
    double x13 = exp(-x12);
    double x14 = -x13 + 1;
    double x15 = x13/x14;
    double x16 = 1.0/x2;
    double x17 = x1*((12.064671139211637*x1 - 3.8216423252915761)*(12.064671139211637*x1 - 3.8216423252915761));
    double x18 = x16*x17;
    double x19 = 241510.76540308626*x18;
    double x20 = pow(x2, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 95058.64000034357*x21;
    double x23 = 72453229.620925874*x10*x18;
    double x24 = x4*(28.150899324827151*x1 - 6.3694038754859594);
    double x25 = exp(x6);
    double x26 = x25 - 1;
    double x27 = 1.0/x26;
    double x28 = Debye(x6);
    double x29 = x28*x4;
    double x30 = -95058.64000034357*x27 + 37415.08178168958*x29;
    double x31 = exp(x12);
    double x32 = x31 - 1;
    double x33 = 1.0/x32;
    double x34 = Debye(x12);
    double x35 = T*x34*x4;
    double x36 = -95058.640000343556*x33 + 124.7169392722986*x35;
    double x37 = x17*x4;

    result += x1*(66213439.079991758*x1 - x15*x22 - x15*x23 + x15*x5 + x19*x9 + x19*exp(-5.0813006666666674*x3)/((x8)*(x8)) + x21*x30 - x21*x36 + x22*x9 - x24*x30 + x24*x36 + 31686.213333447857*x37*(x16*(-9.0*x27 + 3.5424001020211224*x29) + 1.1808000340070408*x20*x28 - 7.621951000000001*x25*x4/((x26)*(x26))) - 31686.213333447853*x37*(0.0039360001133568034*T*x20*x34 - 2286.5853000000002*x10*x31*x4/((x32)*(x32)) + x16*(-9.0*x33 + 0.01180800034007041*x35)) - x5*x9 - 26496933.062196095 - x23*exp(-1524.3902*x11)/((x14)*(x14)))/((V)*(V));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-2492012960.3105969*x2 + 3933559747.5761614*x3 + 378854759.91633976);
    double x5 = 1.0/T;
    double x6 = -11.464926975874729*x2 + 18.097006708817457*x3 + 1.7429853801249982;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 762.19510000000002*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1524.3902*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.01180800034007041*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 31686.213333447853*x7;
    double x25 = x15*x6;

    result += x0*(-165670489788.73367*x13*x16 + x13*x4 - 55223496596.24456*x14*x16 + x14*x4 + x24*(0.0039360001133568034*x17 - 2286.5853000000002*x21 + x23) - x24*(1742824.1113920303*x20*x25 - 2286.5852999999997*x21 + x22 + 3.0*x23 - 3485648.2227840605*x25*exp(x12)/((x19)*(x19)*(x19))) - 110446993192.48912*x16*exp(-2286.5853000000002*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(1748248776.700516*x2 - 553780657.84679925);
    double x5 = 18.097006708817457*pow(x1, 4.0/3.0) - 11.464926975874729*x2 + 1.7429853801249982;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 762.19510000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1524.3902*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 12.064671139211637*x2 - 3.8216423252915761;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0039360001133568034*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2286.5853000000002*x22*x3;
    double x24 = 0.01180800034007041*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 31686.213333447853*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-165670489788.73367*x12*x16 + x12*x4 - 55223496596.24456*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-6859.7559000000001*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 1742824.1113920303*x15*x22 - 3485648.2227840605*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 110446993192.48912*x16*exp(-2286.5853000000002*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(2039623572.8172686*x2 - 461483881.53899938);
    double x4 = 18.097006708817457*pow(x1, 4.0/3.0) - 11.464926975874729*x2 + 1.7429853801249982;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 762.19510000000002*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1524.3902*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 12.064671139211637*x2 - 3.8216423252915761;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0039360001133568034*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2286.5853000000002*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.01180800034007041*x30));
    double x32 = 28.150899324827151*x2 - 6.3694038754859594;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0039360001133568034*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-165670489788.73367*x11*x18 + x11*x3 - 55223496596.24456*x12*x18 + x12*x3 + 31686.213333447853*x13*x31 - 110446993192.48912*x18*exp(-2286.5853000000002*x6)/((x9)*(x9)*(x9)) - 31686.213333447853*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(36.194013417634906*x2 - 11.464926975874729)/pow(x4, 5.0/2.0) - x22*x32 + 1742824.1113920303*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(24.129342278423273*x2 - 7.6432846505831522) - 3485648.2227840605*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = -11.464926975874729*x1 + 18.097006708817457*x2 + 1.7429853801249982;
    double x4 = sqrt(x3);
    double x5 = 2.5406503333333337*x4;
    double x6 = exp(-x5);
    double x7 = -x6 + 1;
    double x8 = pow(x3, -3.0/2.0);
    double x9 = pow(V, -2);
    double x10 = 12.064671139211637*x1 - 3.8216423252915761;
    double x11 = ((x10)*(x10)*(x10))*x9;
    double x12 = x11*x8;
    double x13 = 5.0813006666666674*x4;
    double x14 = exp(-x13)/((x7)*(x7));
    double x15 = pow(x3, -2);
    double x16 = 241510.76540308626*x15;
    double x17 = x11*x16;
    double x18 = x6/x7;
    double x19 = 95058.64000034357*x18;
    double x20 = 1.0/x4;
    double x21 = 93.836331082757169*x1 - 16.985077001295892;
    double x22 = x1*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x4;
    double x25 = 762.19510000000002*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 95058.64000034357*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1524.3902*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = x11*x15;
    double x35 = 72453229.620925874*x23;
    double x36 = x33*x35;
    double x37 = x28*x35;
    double x38 = 24.129342278423273*x1 - 7.6432846505831522;
    double x39 = ((x10)*(x10));
    double x40 = x39*x9;
    double x41 = x38*x40;
    double x42 = x16*x41;
    double x43 = (36.194013417634906*x1 - 11.464926975874729)/pow(x3, 5.0/2.0);
    double x44 = x40*x43;
    double x45 = 28.150899324827151*x1 - 6.3694038754859594;
    double x46 = x10*x2;
    double x47 = x45*x46;
    double x48 = 1.0/x3;
    double x49 = 241510.76540308626*x48;
    double x50 = x14*x49;
    double x51 = x46*(56.301798649654302*x1 - 12.738807750971919);
    double x52 = x46*x8;
    double x53 = 285175.92000103072*x45*x52;
    double x54 = x18*x49;
    double x55 = x15*x41;
    double x56 = x36*x48;
    double x57 = x37*x48;
    double x58 = exp(x5);
    double x59 = x58 - 1;
    double x60 = 1.0/x59;
    double x61 = Debye(x5);
    double x62 = x20*x61;
    double x63 = -3*x60 + 1.1808000340070408*x62;
    double x64 = x1*x63;
    double x65 = 31686.213333447857*x20;
    double x66 = 31686.213333447853*x20;
    double x67 = exp(x25);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x25);
    double x71 = T*x20*x70;
    double x72 = -3*x69 + 0.0039360001133568034*x71;
    double x73 = x1*x72;
    double x74 = x45*x63;
    double x75 = x45*x72;
    double x76 = 1.1808000340070408*x61;
    double x77 = x76*x8;
    double x78 = x58/((x59)*(x59));
    double x79 = 7.621951000000001*x78;
    double x80 = x20*x79;
    double x81 = x48*(-9.0*x60 + 3.5424001020211224*x62) + x77 - x80;
    double x82 = 63372.426666895713*x81;
    double x83 = x20*x47;
    double x84 = 0.0039360001133568034*T*x70;
    double x85 = x8*x84;
    double x86 = x67/((x68)*(x68));
    double x87 = 2286.5853000000002*x23*x86;
    double x88 = x20*x87;
    double x89 = x48*(-9.0*x69 + 0.01180800034007041*x71) + x85 - x88;
    double x90 = 63372.426666895706*x89;
    double x91 = x1*x10*x43;
    double x92 = x1*x39;
    double x93 = x48*x92;
    double x94 = x8*x92;
    double x95 = 3.0*x48;
    double x96 = 3.0*x15;
    double x97 = x64*x96;
    double x98 = x10*x38;
    double x99 = x92*x95;
    double x100 = x30*x93;

    result += (70658488.165856242*x1 + 1840783.2198748188*x12*x14 + 613594.40662493964*x12*x18 + x12*x82 - x12*x90 + 1227188.8132498793*x12*exp(-7.621951000000001*x4)/((x7)*(x7)*(x7)) + x14*x17 + x14*x42 + x17*x18 + x18*x42 - x18*x53 + x19*x22 + x19*x44 - 220711463.59997252*x2 + x21*x64*x65 - x21*x66*x73 - x22*x29 - 55223496596.24456*x28*x31 + x28*x53 - x29*x44 - 165670489788.73367*x31*x33 - x34*x36 - x34*x37 - x36*x55 - x37*x55 + 31686.213333447857*x44*x63 - 31686.213333447853*x44*x72 + x46*x65*(x39*x97 - x45*x77 + x45*x80 - x74*x95 + x76*x91 + 19.364712348800339*x78*x93 - x79*x94 + x81*x99 + x97*x98 - 38.729424697600678*x93*exp(x13)/((x59)*(x59)*(x59))) - x46*x66*(1742824.1113920303*x100*x86 - 3485648.2227840605*x100*exp(x32)/((x68)*(x68)*(x68)) - x45*x85 + x45*x88 + x72*x92*x96 + x73*x96*x98 - x75*x95 + x84*x91 - x87*x94 + x89*x99) - x47*x50 - x47*x54 + x47*x56 + x47*x57 - x50*x51 - x51*x54 + x51*x56 + x51*x57 - 95058.64000034357*x52*x74 + 95058.640000343556*x52*x75 - x82*x83 + x83*x90 - 110446993192.48912*x31*exp(-2286.5853000000002*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.3847;
static const double Vmax = 1.15*2.3847;
static double V = 0.9*2.3847;

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



const char *AlPostPerovskite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *AlPostPerovskite_slb_em_coder_calib_name(void) {
    return "AlPostPerovskite_slb_em";
}

const char *AlPostPerovskite_slb_em_coder_calib_formula(void) {
    return "Al2O3";
}

const double AlPostPerovskite_slb_em_coder_calib_mw(void) {
    return 101.96127999999999;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,3.0,0.0,0.0,0.0,
        0.0,2.0,0.0,0.0,0.0,0.0,
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
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0
    };

const double *AlPostPerovskite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double AlPostPerovskite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int AlPostPerovskite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **AlPostPerovskite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **AlPostPerovskite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void AlPostPerovskite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int AlPostPerovskite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double AlPostPerovskite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int AlPostPerovskite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double AlPostPerovskite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

