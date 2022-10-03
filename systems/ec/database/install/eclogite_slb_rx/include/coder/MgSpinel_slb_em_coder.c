
static char *identifier = "MgSpinel_slb_em.emml:b50c5cd1431d18e1e149688ee9ca72fa6fd19f0d:Thu Feb 10 16:53:14 2022";



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
    double x3 = sqrt(58.211824175067086*x1 - 122.71057094199776*x2 - 5.1364123293499997);
    double x4 = 2.8093680000000001*x3;
    double x5 = 842.81039999999996*x3/T;

    result += 698.41485992487208*T*log(1 - exp(-x5)) - 232.80495330829069*T*Debye(x5) + 116831922.09757045*x1 - 2148369432.3589268*x2 - 209524.45797746163*log(1 - exp(-x4)) + 69841.485992487214*Debye(x4) - 3079050.1555562988 + 7500408834.8521709/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-122.71057094199776*pow(x0, 4.0/3.0) + 58.211824175067086*pow(x0, 2.0/3.0) - 5.1364123293499997);
    double x2 = x1/T;
    double x3 = 842.81039999999996*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -588631.30745922541*x2*x5/x6 + 196210.4358197418*x2*(-0.0035595194363999307*T*x4/x1 + 3/(exp(x3) - 1)) - 232.80495330829069*x4 + 698.41485992487208*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(58.211824175067086*x1 - 122.71057094199776*x3 - 5.1364123293499997);
    double x6 = 2.8093680000000001*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-19.403941391689028*x2 + 81.807047294665168*x4);
    double x10 = 588631.30745922541*x9;
    double x11 = 842.81039999999996*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) - 77887948.065046966*x2 + 2864492576.478569*x4 + 196210.43581974183*x9*(-1.0678558309199793*x8*Debye(x6) + 3/(exp(x6) - 1)) - 196210.4358197418*x9*(-0.0035595194363999307*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 15000817669.704342/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(58.211824175067086*x2 - 122.71057094199776*x3 - 5.1364123293499997);
    double x5 = x0*x4;
    double x6 = 842.81039999999996*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-28879153031.184402*x2 + 60877277202.658272*x3 + 2548197720.8694825);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1685.6207999999999*x5)/((x8)*(x8)) - 196210.4358197418*x4*(x0*(0.010678558309199792*T*x11 - 9.0/x13) + 0.0035595194363999307*x11 - 2528.4312*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 58.211824175067086*x2;
    double x4 = 122.71057094199776*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 5.1364123293499997;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 842.81039999999996*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 496104587.69223273*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(81.807047294665168*x2 - 19.403941391689028)*(-196210.4358197418*x6*(2528.4312*x0*x13*x14/((x15)*(x15)) - 0.0035595194363999307*x12/pow(x5, 3.0/2.0) + (0.010678558309199792*x12*x13 - 9.0/x15)/(-x3 + x4 + 5.1364123293499997)) + x11*x9/x10 + x11*exp(-1685.6207999999999*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 58.211824175067086*x2;
    double x5 = 122.71057094199776*x3;
    double x6 = x4 - x5 - 5.1364123293499997;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(190.88311035421873*x2 - 32.339902319481709);
    double x10 = x8*x9;
    double x11 = 588631.30745922541*x10;
    double x12 = 2.8093680000000001*x7;
    double x13 = exp(-x12);
    double x14 = -x13 + 1;
    double x15 = x13/x14;
    double x16 = 1.0/(-x4 + x5 + 5.1364123293499997);
    double x17 = x3*((81.807047294665168*x2 - 19.403941391689028)*(81.807047294665168*x2 - 19.403941391689028));
    double x18 = x16*x17;
    double x19 = 1653681.9589741093*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 588631.30745922541*x21;
    double x23 = 1.0/T;
    double x24 = x23*x7;
    double x25 = 842.81039999999996*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 496104587.69223273*x18*x23;
    double x30 = exp(x12);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x12);
    double x34 = x33*x8;
    double x35 = -588631.30745922553*x32 + 209524.45797746169*x34;
    double x36 = exp(x25);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = T*Debye(x25);
    double x40 = x39*x8;
    double x41 = -3*x38 + 0.0035595194363999307*x40;
    double x42 = 196210.4358197418*x8;

    result += x0*(45002453009.113022*x0 + x10*x35 + x11*x15 - x11*x28 - x15*x19 + x15*x22 + x17*x42*(x16*(-9.0*x38 + 0.010678558309199792*x40) - 0.0035595194363999307*x20*x39 + 2528.4312*x23*x36*x8/((x37)*(x37))) - 196210.43581974183*x17*x8*(x16*(-9.0000000000000018*x32 + 3.2035674927599382*x34) - 1.0678558309199793*x20*x33 + 8.4281040000000012*x30*x8/((x31)*(x31))) + 129813246.77507827*x2 + x21*x35 - 196210.4358197418*x21*x41 - x22*x28 + x28*x29 - 6683816011.7833271*x3 - x41*x42*x9 + x29*exp(-1685.6207999999999*x24)/((x27)*(x27)) - x19*exp(-5.6187360000000002*x7)/((x14)*(x14)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-86637459093.553207*x2 + 182631831607.97482*x3 + 7644593162.6084471);
    double x5 = 1.0/T;
    double x6 = 58.211824175067086*x2;
    double x7 = 122.71057094199776*x3;
    double x8 = x6 - x7 - 5.1364123293499997;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 842.81039999999996*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1685.6207999999999*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2528.4312*x0*x22*x9;
    double x24 = 0.010678558309199792*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 196210.4358197418*x9;
    double x27 = x17*(-x6 + x7 + 5.1364123293499997);

    result += x0*(-1254366317984.1772*x15*x18 - x15*x4 - 418122105994.72571*x16*x18 - x16*x4 + x26*(0.0035595194363999307*x19 + x23 + x25) - x26*(-2130988.11104448*x22*x27 + x23 + x24 + 3.0*x25 + 4261976.22208896*x27*exp(x14)/((x21)*(x21)*(x21))) - 836244211989.45142*x18*exp(-2528.4312*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(81169702936.877686*x2 - 19252768687.456268);
    double x5 = 58.211824175067086*x2;
    double x6 = 122.71057094199776*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 5.1364123293499997;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 842.81039999999996*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1685.6207999999999*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 81.807047294665168*x2 - 19.403941391689028;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0035595194363999307*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2528.4312*x3;
    double x26 = 0.010678558309199792*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 196210.4358197418*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 5.1364123293499997);

    result += x0*x1*x2*(1254366317984.1772*x14*x18 - x14*x4 + 418122105994.72571*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7585.2936*x0*x31 - x26*x30 + 3.0*x27*x32) + 2130988.11104448*x16*x24 - 4261976.22208896*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 836244211989.45142*x18*exp(-2528.4312*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(94697986759.690643*x2 - 16043973906.213554);
    double x4 = 58.211824175067086*x2;
    double x5 = 122.71057094199776*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 5.1364123293499997;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 842.81039999999996*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1685.6207999999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 81.807047294665168*x2 - 19.403941391689028;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0035595194363999307*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2528.4312*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 5.1364123293499997;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.010678558309199792*x33));
    double x35 = 190.88311035421873*x2 - 32.339902319481709;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0035595194363999307*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(1254366317984.1772*x13*x20 + x13*x3 + 418122105994.72571*x14*x20 + x14*x3 + 196210.4358197418*x15*x34 + 196210.4358197418*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(245.42114188399552*x2 - 58.211824175067079)/pow(x6, 5.0/2.0) + x24*x35 - 2130988.11104448*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(163.61409458933034*x2 - 38.807882783378055) + 4261976.22208896*x37*exp(x12)/((x26)*(x26)*(x26))) + 836244211989.45142*x20*exp(-2528.4312*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 58.211824175067086*x2;
    double x5 = 122.71057094199776*x3;
    double x6 = x4 - x5 - 5.1364123293499997;
    double x7 = sqrt(x6);
    double x8 = 2.8093680000000001*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 81.807047294665168*x2 - 19.403941391689028;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.6187360000000002*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 5.1364123293499997;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1653681.9589741093*x16;
    double x21 = x9/x10;
    double x22 = 1653681.9589741093*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(636.27703451406239*x2 - 86.239739518617881);
    double x25 = 588631.30745922541*x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 842.81039999999996*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = pow(T, -2);
    double x33 = x14*x32;
    double x34 = 1685.6207999999999*x27;
    double x35 = exp(-x34)/((x30)*(x30));
    double x36 = 496104587.69223273*x26;
    double x37 = x19*x36;
    double x38 = x18*(163.61409458933034*x2 - 38.807882783378055);
    double x39 = ((x12)*(x12));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = (245.42114188399552*x2 - 58.211824175067079)/pow(x6, 5.0/2.0);
    double x43 = x40*x42;
    double x44 = 588631.30745922541*x43;
    double x45 = 1.0/x17;
    double x46 = 190.88311035421873*x2 - 32.339902319481709;
    double x47 = x45*x46;
    double x48 = x12*x3;
    double x49 = x20*x48;
    double x50 = x45*(381.76622070843746*x2 - 64.679804638963418);
    double x51 = x46*x48;
    double x52 = x11*x51;
    double x53 = 1765893.9223776762*x52;
    double x54 = x22*x48;
    double x55 = x36*x41;
    double x56 = x36*x48;
    double x57 = x35*x56;
    double x58 = x31*x56;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.0678558309199793*x63;
    double x65 = 196210.43581974183*x23;
    double x66 = exp(x28);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x28);
    double x70 = T*x23*x69;
    double x71 = -3*x68 + 0.0035595194363999307*x70;
    double x72 = 196210.4358197418*x23;
    double x73 = 1.0678558309199793*x62;
    double x74 = x11*x73;
    double x75 = x59/((x60)*(x60));
    double x76 = 8.4281040000000012*x75;
    double x77 = x23*x76;
    double x78 = x45*(-9.0000000000000018*x61 + 3.2035674927599382*x63) - x74 + x77;
    double x79 = 392420.87163948367*x78;
    double x80 = x23*x51;
    double x81 = 0.0035595194363999307*T*x69;
    double x82 = x11*x81;
    double x83 = x66/((x67)*(x67));
    double x84 = 2528.4312*x26*x83;
    double x85 = x23*x84;
    double x86 = x45*(-9.0*x68 + 0.010678558309199792*x70) - x82 + x85;
    double x87 = 392420.87163948361*x86;
    double x88 = x12*x2;
    double x89 = x42*x88;
    double x90 = x2*x39;
    double x91 = x45*x90;
    double x92 = x11*x90;
    double x93 = 3.0000000000000004*x64;
    double x94 = x18*x90;
    double x95 = x38*x88;
    double x96 = x32*x91;
    double x97 = 3.0*x71;

    result += (-180009812036.45209*x0 - 13937403.533157527*x14*x16 - 4645801.1777191758*x14*x21 + x14*x79 - x14*x87 - x19*x20 - x19*x22 - 346168658.06687534*x2 - x20*x41 - x21*x25 - x21*x44 - x21*x53 - x22*x41 - x24*x64*x65 + x24*x71*x72 + x25*x31 + 22279386705.944424*x3 + 418122105994.72571*x31*x33 + x31*x37 + x31*x44 + x31*x53 + x31*x55 + 1254366317984.1772*x33*x35 + x35*x37 + x35*x55 - 196210.43581974183*x43*x64 + 196210.4358197418*x43*x71 + x47*x49 + x47*x54 - x47*x57 - x47*x58 - x48*x65*(x46*x74 - x46*x77 - x47*x93 + x73*x89 - 23.677645678272004*x75*x91 - x76*x92 + 3.0000000000000004*x78*x91 + x93*x94 + x93*x95 + 47.355291356544008*x91*exp(x15)/((x60)*(x60)*(x60))) + x48*x72*(x46*x82 - x46*x85 - x47*x97 + x81*x89 - 2130988.11104448*x83*x96 - x84*x92 + 3.0*x86*x91 + x94*x97 + x95*x97 + 4261976.22208896*x96*exp(x34)/((x67)*(x67)*(x67))) + x49*x50 + x50*x54 - x50*x57 - x50*x58 - 588631.30745922553*x52*x64 + 588631.30745922541*x52*x71 + x79*x80 - x80*x87 + 836244211989.45142*x33*exp(-2528.4312*x27)/((x30)*(x30)*(x30)) - 9291602.3554383516*x14*exp(-8.4281040000000012*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*15.9048;
static const double Vmax = 1.15*15.9048;
static double V = 0.9*15.9048;

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



const char *MgSpinel_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *MgSpinel_slb_em_coder_calib_name(void) {
    return "MgSpinel_slb_em";
}

const char *MgSpinel_slb_em_coder_calib_formula(void) {
    return "Mg4Al8O16";
}

const double MgSpinel_slb_em_coder_calib_mw(void) {
    return 569.06272;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,16.0,0.0,0.0,0.0,
        4.0,8.0,0.0,0.0,0.0,0.0,
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

const double *MgSpinel_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double MgSpinel_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double MgSpinel_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double MgSpinel_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double MgSpinel_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double MgSpinel_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double MgSpinel_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double MgSpinel_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double MgSpinel_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double MgSpinel_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double MgSpinel_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double MgSpinel_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double MgSpinel_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int MgSpinel_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **MgSpinel_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **MgSpinel_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void MgSpinel_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int MgSpinel_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double MgSpinel_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int MgSpinel_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double MgSpinel_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

