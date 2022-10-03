
static char *identifier = "Jadeite_slb_em.emml:6cc8135f0217fbf233338a740711ccd94531b672:Thu Feb 10 16:52:48 2022";



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
    double x3 = sqrt(-1.083848401829469*x1 + 16.734165453935571*x2 - 0.19129629500000012);
    double x4 = 2.7358743333333333*x3;
    double x5 = 820.76229999999998*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) - 4558127.371876426*x1 - 91658972.795326322*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) + 830947.06115761958 + 219657372.85005447/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(16.734165453935571*pow(x0, 4.0/3.0) - 1.083848401829469*pow(x0, 2.0/3.0) - 0.19129629500000012);
    double x2 = x1/T;
    double x3 = 820.76229999999998*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -204725.92385218426*x2*x5/x6 + 68241.97461739475*x2*(-0.0036551386436730828*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-1.083848401829469*x1 + 16.734165453935571*x3 - 0.19129629500000012);
    double x6 = 2.7358743333333333*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(0.36128280060982298*x2 - 11.156110302623713*x4);
    double x10 = 820.76229999999998*x5/T;
    double x11 = exp(-x10);
    double x12 = 68241.97461739475*x9;

    result += 204725.92385218426*x11*x9/(-x11 + 1) + x12*(-1.0965415931019249*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.0036551386436730828*T*x8*Debye(x10) + 3/(exp(x10) - 1)) + 3038751.5812509507*x2 + 122211963.72710176*x4 - 204725.92385218423*x7*x9/(-x7 + 1) - 439314745.70010895/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-1.083848401829469*x2 + 16.734165453935571*x3 - 0.19129629500000012);
    double x5 = x0*x4;
    double x6 = 820.76229999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(182120477.78078559*x2 - 2811863912.5077319*x3 + 32143768.984931931);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1641.5246*x5)/((x8)*(x8)) - 68241.97461739475*x4*(x0*(0.010965415931019249*T*x11 - 9.0/x13) + 0.0036551386436730828*x11 - 2462.2869000000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 1.083848401829469*x2;
    double x4 = 16.734165453935571*pow(x1, 4.0/3.0);
    double x5 = -x3 + x4 - 0.19129629500000012;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 820.76229999999998*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 168031320.13054362*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(11.156110302623713*x2 - 0.36128280060982298)*(68241.97461739475*x6*(2462.2869000000001*x0*x13*x14/((x15)*(x15)) - 0.0036551386436730828*x12/pow(x5, 3.0/2.0) + (0.010965415931019249*x12*x13 - 9.0/x15)/(x3 - x4 + 0.19129629500000012)) - x11*x9/x10 - x11*exp(-1641.5246*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 1.083848401829469*x2;
    double x5 = 16.734165453935571*x3;
    double x6 = -x4 + x5 - 0.19129629500000012;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(26.030924039455329*x2 - 0.6021380010163716);
    double x10 = x8*x9;
    double x11 = 2.7358743333333333*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 204725.92385218423*x14;
    double x16 = 1.0/(x4 - x5 + 0.19129629500000012);
    double x17 = x3*((11.156110302623713*x2 - 0.36128280060982298)*(11.156110302623713*x2 - 0.36128280060982298));
    double x18 = x16*x17;
    double x19 = 560104.40043514525*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 820.76229999999998*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 204725.92385218426*x27;
    double x29 = 168031320.13054362*x18*x22;
    double x30 = exp(x11);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x11);
    double x34 = x33*x8;
    double x35 = -3*x32 + 1.0965415931019249*x34;
    double x36 = 68241.97461739475*x8;
    double x37 = x36*x9;
    double x38 = 68241.97461739475*x21;
    double x39 = exp(x24);
    double x40 = x39 - 1;
    double x41 = 1.0/x40;
    double x42 = Debye(x24);
    double x43 = T*x42*x8;
    double x44 = -3*x41 + 0.0036551386436730828*x43;
    double x45 = x17*x36;

    result += x0*(1317944237.1003268*x0 - x10*x15 + x10*x28 - x14*x19 + x15*x21 - 5064585.9687515842*x2 - x21*x28 + x27*x29 - 285161248.69657075*x3 - x35*x37 + x35*x38 + x37*x44 - x38*x44 - x45*(x16*(-9.0*x32 + 3.2896247793057745*x34) - 1.0965415931019249*x20*x33 + 8.2076229999999999*x30*x8/((x31)*(x31))) + x45*(-0.0036551386436730828*T*x20*x42 + x16*(-9.0*x41 + 0.010965415931019249*x43) + 2462.2869000000001*x22*x39*x8/((x40)*(x40))) + x29*exp(-1641.5246*x23)/((x26)*(x26)) - x19*exp(-5.4717486666666666*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(546361433.3423568*x2 - 8435591737.5231953*x3 + 96431306.954795793);
    double x5 = 1.0/T;
    double x6 = 1.083848401829469*x2;
    double x7 = 16.734165453935571*x3;
    double x8 = -x6 + x7 - 0.19129629500000012;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 820.76229999999998*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1641.5246*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2462.2869000000001*x0*x22*x9;
    double x24 = 0.010965415931019249*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 68241.97461739475*x9;
    double x27 = x17*(x6 - x7 + 0.19129629500000012);

    result += x0*(-413741318347.14386*x15*x18 - x15*x4 - 137913772782.38129*x16*x18 - x16*x4 + x26*(0.0036551386436730828*x19 + x23 + x25) - x26*(-2020952.2593038699*x22*x27 + x23 + x24 + 3.0*x25 + 4041904.5186077398*x27*exp(x14)/((x21)*(x21)*(x21))) - 275827545564.76257*x18*exp(-2462.2869000000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(3749151883.3436418*x2 - 121413651.85385706);
    double x5 = 1.083848401829469*x2;
    double x6 = 16.734165453935571*pow(x1, 4.0/3.0);
    double x7 = -x5 + x6 - 0.19129629500000012;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 820.76229999999998*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1641.5246*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = 11.156110302623713*x2 - 0.36128280060982298;
    double x17 = pow(T, -3);
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0036551386436730828*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2462.2869000000001*x24*x3;
    double x26 = 0.010965415931019249*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 68241.97461739475*x16;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = 1.0/(x5 - x6 + 0.19129629500000012);

    result += x0*x1*x2*(-413741318347.14386*x14*x18 + x14*x4 - 137913772782.38129*x15*x18 + x15*x4 + x19*x29*(x19*x21 - x25*x8 + x28) - x29*x8*(x0*(-7386.8607000000002*x0*x19*x24 + x26*x30 - 3.0*x27*x31) + 2020952.2593038699*x17*x24 - 4041904.5186077398*x17*exp(x13)/((x23)*(x23)*(x23)) + x19*x25 + x21*x30 - x28*x31) - 275827545564.76257*x18*exp(-2462.2869000000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(4374010530.5675821*x2 - 101178043.21154754);
    double x4 = 1.083848401829469*x2;
    double x5 = 16.734165453935571*pow(x1, 4.0/3.0);
    double x6 = -x4 + x5 - 0.19129629500000012;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 820.76229999999998*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1641.5246*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -2);
    double x16 = 1.0/x7;
    double x17 = 11.156110302623713*x2 - 0.36128280060982298;
    double x18 = ((x17)*(x17));
    double x19 = x18*x2;
    double x20 = x16*x19;
    double x21 = x15*x20;
    double x22 = pow(x6, -3.0/2.0);
    double x23 = Debye(x9);
    double x24 = 0.0036551386436730828*T*x23;
    double x25 = x22*x24;
    double x26 = exp(x9);
    double x27 = x26 - 1;
    double x28 = x26/((x27)*(x27));
    double x29 = 2462.2869000000001*x28;
    double x30 = x0*x16*x29;
    double x31 = x4 - x5 + 0.19129629500000012;
    double x32 = 1.0/x31;
    double x33 = 1.0/x27;
    double x34 = T*x16*x23;
    double x35 = x32*(-9.0*x33 + 0.010965415931019249*x34);
    double x36 = 26.030924039455329*x2 - 0.6021380010163716;
    double x37 = x17*x2;
    double x38 = x19*x32;
    double x39 = x15*x38;
    double x40 = x0*x2;
    double x41 = -9.0*x33 + 0.010965415931019249*x34;
    double x42 = x41/((x31)*(x31));

    result += x40*(-413741318347.14386*x13*x21 + x13*x3 - 137913772782.38129*x14*x21 + x14*x3 - 68241.97461739475*x20*(-x25 + x30 + x35) - 68241.97461739475*x7*(-x18*x22*x29*x40 + x19*x42 + x24*x37*(33.468330907871135*x2 - 1.083848401829469)/pow(x6, 5.0/2.0) - x25*x36 - 2020952.2593038699*x28*x39 + x30*x36 + x32*x36*x41 + x37*x42*(22.312220605247425*x2 - 0.72256560121964597) - 3.0*x38*(x25 - x30 - x35) + 4041904.5186077398*x39*exp(x12)/((x27)*(x27)*(x27))) - 275827545564.76257*x21*exp(-2462.2869000000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 1.083848401829469*x2;
    double x5 = 16.734165453935571*x3;
    double x6 = -x4 + x5 - 0.19129629500000012;
    double x7 = sqrt(x6);
    double x8 = 2.7358743333333333*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 11.156110302623713*x2 - 0.36128280060982298;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = x4 - x5 + 0.19129629500000012;
    double x16 = pow(x15, -2);
    double x17 = 560104.40043514525*x16;
    double x18 = x13*x17;
    double x19 = 5.4717486666666666*x7;
    double x20 = exp(-x19)/((x10)*(x10));
    double x21 = x9/x10;
    double x22 = 204725.92385218423*x21;
    double x23 = 1.0/x7;
    double x24 = 86.769746798184428*x2 - 1.6057013360436576;
    double x25 = x2*x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 820.76229999999998*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 204725.92385218426*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1641.5246*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 168031320.13054362*x26;
    double x38 = x16*x37;
    double x39 = x13*x38;
    double x40 = 22.312220605247425*x2 - 0.72256560121964597;
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x40*x42;
    double x44 = x17*x43;
    double x45 = (33.468330907871135*x2 - 1.083848401829469)/pow(x6, 5.0/2.0);
    double x46 = x42*x45;
    double x47 = 26.030924039455329*x2 - 0.6021380010163716;
    double x48 = x12*x3;
    double x49 = x47*x48;
    double x50 = 1.0/x15;
    double x51 = 560104.40043514525*x50;
    double x52 = x20*x51;
    double x53 = x48*(52.061848078910657*x2 - 1.2042760020327432);
    double x54 = x21*x51;
    double x55 = x11*x48;
    double x56 = 614177.77155655273*x47*x55;
    double x57 = x38*x43;
    double x58 = x37*x50;
    double x59 = x36*x58;
    double x60 = x31*x58;
    double x61 = exp(x8);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x8);
    double x65 = x23*x64;
    double x66 = -3*x63 + 1.0965415931019249*x65;
    double x67 = x2*x66;
    double x68 = 68241.97461739475*x23;
    double x69 = x24*x68;
    double x70 = exp(x28);
    double x71 = x70 - 1;
    double x72 = 1.0/x71;
    double x73 = T*Debye(x28);
    double x74 = x23*x73;
    double x75 = -3*x72 + 0.0036551386436730828*x74;
    double x76 = x2*x75;
    double x77 = 68241.97461739475*x46;
    double x78 = x47*x66;
    double x79 = 204725.92385218426*x55;
    double x80 = x47*x75;
    double x81 = 1.0965415931019249*x64;
    double x82 = x11*x81;
    double x83 = x61/((x62)*(x62));
    double x84 = 8.2076229999999999*x83;
    double x85 = x23*x84;
    double x86 = x50*(-9.0*x63 + 3.2896247793057745*x65) - x82 + x85;
    double x87 = 136483.9492347895*x86;
    double x88 = x23*x49;
    double x89 = 0.0036551386436730828*x73;
    double x90 = x11*x89;
    double x91 = x70/((x71)*(x71));
    double x92 = 2462.2869000000001*x26*x91;
    double x93 = x23*x92;
    double x94 = x50*(-9.0*x72 + 0.010965415931019249*x74) - x90 + x93;
    double x95 = 136483.9492347895*x94;
    double x96 = x12*x2*x45;
    double x97 = x2*x41;
    double x98 = x50*x97;
    double x99 = x11*x97;
    double x100 = 3.0*x50;
    double x101 = 3.0*x16;
    double x102 = x101*x67;
    double x103 = x12*x40;
    double x104 = x100*x97;
    double x105 = x48*x68;
    double x106 = x33*x98;

    result += (-5271776948.4013071*x0 + x105*(x100*x78 + x102*x103 + x102*x41 + x104*x86 - x47*x82 + x47*x85 + x81*x96 - 22.455025103376332*x83*x98 - x84*x99 + 44.910050206752665*x98*exp(x19)/((x62)*(x62)*(x62))) - x105*(x100*x80 + x101*x103*x76 + x101*x75*x97 + x104*x94 - 2020952.2593038699*x106*x91 + 4041904.5186077398*x106*exp(x35)/((x71)*(x71)*(x71)) - x47*x90 + x47*x93 + x89*x96 - x92*x99) + 4597125.7594127078*x14*x20 + 1532375.2531375694*x14*x21 - x14*x87 + x14*x95 + x18*x20 + x18*x21 + 13505562.583337557*x2 + x20*x44 + x21*x44 - x21*x56 + x22*x25 + x22*x46 - x25*x32 + 950537495.65523577*x3 - 137913772782.38129*x31*x34 - x31*x39 + x31*x56 - x31*x57 - x32*x46 - 413741318347.14386*x34*x36 - x36*x39 - x36*x57 + x49*x52 + x49*x54 - x49*x59 - x49*x60 + x52*x53 + x53*x54 - x53*x59 - x53*x60 + x66*x77 + x67*x69 - x69*x76 - x75*x77 - x78*x79 + x79*x80 + x87*x88 - x88*x95 - 275827545564.76257*x34*exp(-2462.2869000000001*x27)/((x30)*(x30)*(x30)) + 3064750.5062751388*x14*exp(-8.2076229999999999*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.0508;
static const double Vmax = 1.15*6.0508;
static double V = 0.9*6.0508;

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



const char *Jadeite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Jadeite_slb_em_coder_calib_name(void) {
    return "Jadeite_slb_em";
}

const char *Jadeite_slb_em_coder_calib_formula(void) {
    return "NaAlSi2O6";
}

const double Jadeite_slb_em_coder_calib_mw(void) {
    return 202.13871;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,6.0,0.0,0.0,1.0,
        0.0,1.0,2.0,0.0,0.0,0.0,
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

const double *Jadeite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Jadeite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Jadeite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Jadeite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Jadeite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Jadeite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Jadeite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Jadeite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Jadeite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Jadeite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Jadeite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Jadeite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Jadeite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Jadeite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Jadeite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Jadeite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Jadeite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Jadeite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Jadeite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Jadeite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Jadeite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Jadeite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Jadeite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Jadeite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Jadeite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Jadeite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Jadeite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Jadeite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Jadeite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Jadeite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Jadeite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

