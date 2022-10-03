
static char *identifier = "MgCaFerrite_slb_em.emml:94ce0b160a4722f550f7a6936009855c506b48e4:Thu Feb 10 16:52:58 2022";



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
    double x3 = sqrt(-4.0305790459603035*x1 + 15.675131736893716*x2 - 0.11217664879999978);
    double x4 = 2.7954303333333335*x3;
    double x5 = 838.62909999999999*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 38810626.521710575*x1 + 43845212.956688762*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 6225413.2173550297 + 2962437.2956025009/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(15.675131736893716*pow(x0, 4.0/3.0) - 4.0305790459603035*pow(x0, 2.0/3.0) - 0.11217664879999978);
    double x2 = x1/T;
    double x3 = 838.62909999999999*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -146427.75635135538*x2*x5/x6 + 48809.252117118464*x2*(-0.0035772667559472952*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-4.0305790459603035*x1 + 15.675131736893716*x3 - 0.11217664879999978);
    double x6 = 2.7954303333333335*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(1.3435263486534343*x2 - 10.45008782459581*x4);
    double x10 = 838.62909999999999*x5/T;
    double x11 = exp(-x10);

    result += 146427.75635135538*x11*x9/(-x11 + 1) + 25873751.014473714*x2 - 58460283.942251682*x4 - 146427.75635135541*x7*x9/(-x7 + 1) + 48809.252117118471*x9*(-1.0731800267841887*x8*Debye(x6) + 3/(exp(x6) - 1)) - 48809.252117118464*x9*(-0.0035772667559472952*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 5924874.5912050018/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-4.0305790459603035*x2 + 15.675131736893716*x3 - 0.11217664879999978);
    double x5 = x0*x4;
    double x6 = 838.62909999999999*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(494949373.44179076*x2 - 1924883879.791173*x3 + 13775132.904044408);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1677.2582*x5)/((x8)*(x8)) - 48809.252117118464*x4*(x0*(0.010731800267841884*T*x11 - 8.9999999999999982/x13) + 0.0035772667559472952*x11 - 2515.8872999999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 4.0305790459603035*x2;
    double x4 = 15.675131736893716*pow(x1, 4.0/3.0);
    double x5 = -x3 + x4 - 0.11217664879999978;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 838.62909999999999*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 122798577.52395645*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(10.45008782459581*x2 - 1.3435263486534343)*(48809.252117118464*x6*(2515.8872999999999*x0*x13*x14/((x15)*(x15)) - 0.0035772667559472952*x12/pow(x5, 3.0/2.0) + (0.010731800267841884*x12*x13 - 8.9999999999999982/x15)/(x3 - x4 + 0.11217664879999978)) - x11*x9/x10 - x11*exp(-1677.2582*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 4.0305790459603035*x2;
    double x5 = 15.675131736893716*x3;
    double x6 = -x4 + x5 - 0.11217664879999978;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(24.383538257390221*x2 - 2.2392105810890572);
    double x10 = x8*x9;
    double x11 = 2.7954303333333335*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 146427.75635135541*x14;
    double x16 = 1.0/(x4 - x5 + 0.11217664879999978);
    double x17 = x3*((10.45008782459581*x2 - 1.3435263486534343)*(10.45008782459581*x2 - 1.3435263486534343));
    double x18 = x16*x17;
    double x19 = 409328.59174652159*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 838.62909999999999*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 146427.75635135538*x27;
    double x29 = 122798577.52395645*x18*x22;
    double x30 = exp(x11);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x11);
    double x34 = x33*x8;
    double x35 = -146427.75635135541*x32 + 52381.114494365414*x34;
    double x36 = exp(x24);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = Debye(x24);
    double x40 = T*x39*x8;
    double x41 = -3*x38 + 0.0035772667559472952*x40;
    double x42 = 48809.252117118464*x8;

    result += x0*(17774623.773615006*x0 - x10*x15 + x10*x28 - x10*x35 - x14*x19 + x15*x21 + x17*x42*(-0.0035772667559472952*T*x20*x39 + x16*(-8.9999999999999982*x38 + 0.010731800267841884*x40) + 2515.8872999999999*x22*x36*x8/((x37)*(x37))) - 48809.252117118471*x17*x8*(x16*(-9.0000000000000018*x32 + 3.2195400803525667*x34) - 1.0731800267841887*x20*x33 + 8.3862909999999999*x30*x8/((x31)*(x31))) - 43122918.357456192*x2 - x21*x28 + x21*x35 - 48809.252117118464*x21*x41 + x27*x29 + 136407329.19858724*x3 + x41*x42*x9 + x29*exp(-1677.2582*x23)/((x26)*(x26)) - x19*exp(-5.5908606666666669*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(1484848120.3253725*x2 - 5774651639.3735199*x3 + 41325398.712133229);
    double x5 = 1.0/T;
    double x6 = 4.0305790459603035*x2;
    double x7 = 15.675131736893716*x3;
    double x8 = -x6 + x7 - 0.11217664879999978;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 838.62909999999999*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1677.2582*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.010731800267841884*x19;
    double x25 = x5*(T*x24 - 8.9999999999999982/x21);
    double x26 = 48809.252117118464*x9;
    double x27 = x17*(x6 - x7 + 0.11217664879999978);

    result += x0*(-308947381650.58746*x15*x18 - x15*x4 - 102982460550.19582*x16*x18 - x16*x4 + x26*(0.0035772667559472952*x19 - 2515.8872999999999*x23 + x25) - x26*(-2109896.3021004298*x22*x27 - 2515.8872999999985*x23 + x24 + 2.9999999999999996*x25 + 4219792.6042008596*x27*exp(x14)/((x21)*(x21)*(x21))) - 205964921100.39163*x18*exp(-2515.8872999999999*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(2566511839.7215638*x2 - 329966248.9611938);
    double x5 = 4.0305790459603035*x2;
    double x6 = 15.675131736893716*pow(x1, 4.0/3.0);
    double x7 = -x5 + x6 - 0.11217664879999978;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 838.62909999999999*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1677.2582*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = 10.45008782459581*x2 - 1.3435263486534343;
    double x17 = pow(T, -3);
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0035772667559472952*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2515.8872999999999*x24*x3;
    double x26 = 0.010731800267841884*T*x20;
    double x27 = x19*x26 - 8.9999999999999982/x23;
    double x28 = x0*x27;
    double x29 = 48809.252117118464*x16;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = 1.0/(x5 - x6 + 0.11217664879999978);

    result += x0*x1*x2*(-308947381650.58746*x14*x18 + x14*x4 - 102982460550.19582*x15*x18 + x15*x4 + x19*x29*(x19*x21 - x25*x8 + x28) - x29*x8*(x0*(-7547.6618999999982*x0*x19*x24 + x26*x30 - 2.9999999999999996*x27*x31) + 2109896.3021004298*x17*x24 - 4219792.6042008596*x17*exp(x13)/((x23)*(x23)*(x23)) + x19*x25 + x21*x30 - x28*x31) - 205964921100.39163*x18*exp(-2515.8872999999999*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(2994263813.008491*x2 - 274971874.13432819);
    double x4 = 4.0305790459603035*x2;
    double x5 = 15.675131736893716*pow(x1, 4.0/3.0);
    double x6 = -x4 + x5 - 0.11217664879999978;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 838.62909999999999*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1677.2582*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -2);
    double x16 = 1.0/x7;
    double x17 = 10.45008782459581*x2 - 1.3435263486534343;
    double x18 = ((x17)*(x17));
    double x19 = x18*x2;
    double x20 = x16*x19;
    double x21 = x15*x20;
    double x22 = pow(x6, -3.0/2.0);
    double x23 = Debye(x9);
    double x24 = 0.0035772667559472952*T*x23;
    double x25 = x22*x24;
    double x26 = exp(x9);
    double x27 = x26 - 1;
    double x28 = x26/((x27)*(x27));
    double x29 = 2515.8872999999999*x28;
    double x30 = x0*x16*x29;
    double x31 = x4 - x5 + 0.11217664879999978;
    double x32 = 1.0/x31;
    double x33 = 1.0/x27;
    double x34 = T*x16*x23;
    double x35 = x32*(-8.9999999999999982*x33 + 0.010731800267841884*x34);
    double x36 = 24.383538257390221*x2 - 2.2392105810890572;
    double x37 = x17*x2;
    double x38 = x19*x32;
    double x39 = x15*x38;
    double x40 = x0*x2;
    double x41 = -8.9999999999999982*x33 + 0.010731800267841884*x34;
    double x42 = x41/((x31)*(x31));

    result += x40*(-308947381650.58746*x13*x21 + x13*x3 - 102982460550.19582*x14*x21 + x14*x3 - 48809.252117118464*x20*(-x25 + x30 + x35) - 48809.252117118464*x7*(-x18*x22*x29*x40 + x19*x42 + x24*x37*(31.350263473787429*x2 - 4.0305790459603035)/pow(x6, 5.0/2.0) - x25*x36 - 2109896.3021004298*x28*x39 + x30*x36 + x32*x36*x41 + x37*x42*(20.90017564919162*x2 - 2.6870526973068687) - 2.9999999999999996*x38*(x25 - x30 - x35) + 4219792.6042008596*x39*exp(x12)/((x27)*(x27)*(x27))) - 205964921100.39163*x21*exp(-2515.8872999999999*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 4.0305790459603035*x2;
    double x5 = 15.675131736893716*x3;
    double x6 = -x4 + x5 - 0.11217664879999978;
    double x7 = sqrt(x6);
    double x8 = 2.7954303333333335*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 10.45008782459581*x2 - 1.3435263486534343;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = x4 - x5 + 0.11217664879999978;
    double x16 = pow(x15, -2);
    double x17 = x13*x16;
    double x18 = 5.5908606666666669*x7;
    double x19 = exp(-x18)/((x10)*(x10));
    double x20 = 409328.59174652159*x19;
    double x21 = x9/x10;
    double x22 = 409328.59174652159*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(81.278460857967403*x2 - 5.9712282162374857);
    double x25 = x23*x24;
    double x26 = 146427.75635135541*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 838.62909999999999*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 146427.75635135538*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1677.2582*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 122798577.52395645*x27;
    double x39 = x17*x38;
    double x40 = x16*(20.90017564919162*x2 - 2.6870526973068687);
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x40*x42;
    double x44 = (31.350263473787429*x2 - 4.0305790459603035)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 1.0/x15;
    double x47 = 24.383538257390221*x2 - 2.2392105810890572;
    double x48 = x46*x47;
    double x49 = x12*x3;
    double x50 = x20*x49;
    double x51 = x46*(48.767076514780442*x2 - 4.4784211621781145);
    double x52 = x22*x49;
    double x53 = x47*x49;
    double x54 = x11*x53;
    double x55 = x38*x43;
    double x56 = x38*x49;
    double x57 = x37*x56;
    double x58 = x32*x56;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.0731800267841887*x63;
    double x65 = 48809.252117118471*x23;
    double x66 = exp(x29);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = T*Debye(x29);
    double x70 = x23*x69;
    double x71 = -3*x68 + 0.0035772667559472952*x70;
    double x72 = 48809.252117118464*x23;
    double x73 = 146427.75635135541*x54;
    double x74 = 1.0731800267841887*x62;
    double x75 = x11*x74;
    double x76 = x59/((x60)*(x60));
    double x77 = 8.3862909999999999*x76;
    double x78 = x23*x77;
    double x79 = x46*(-9.0000000000000018*x61 + 3.2195400803525667*x63);
    double x80 = -97618.504234236942*x75 + 97618.504234236942*x78 + 97618.504234236942*x79;
    double x81 = x23*x53;
    double x82 = 0.0035772667559472952*x69;
    double x83 = x11*x82;
    double x84 = x66/((x67)*(x67));
    double x85 = 2515.8872999999999*x27*x84;
    double x86 = x23*x85;
    double x87 = x46*(-8.9999999999999982*x68 + 0.010731800267841884*x70);
    double x88 = -97618.504234236927*x83 + 97618.504234236927*x86 + 97618.504234236927*x87;
    double x89 = x12*x2;
    double x90 = x44*x89;
    double x91 = x2*x41;
    double x92 = x46*x91;
    double x93 = x11*x91;
    double x94 = 3.0000000000000004*x64;
    double x95 = x16*x91;
    double x96 = x40*x89;
    double x97 = x34*x92;
    double x98 = 2.9999999999999996*x71;

    result += (-71098495.094460025*x0 + 3432748.6850065282*x14*x19 + 1144249.5616688428*x14*x21 - x14*x80 + x14*x88 + x17*x20 + x17*x22 + 114994448.95321651*x2 + x20*x43 - 439283.26905406622*x21*x54 + x22*x43 + x24*x64*x65 - x24*x71*x72 + x25*x26 - x25*x33 + x26*x45 - 454691097.32862413*x3 - 102982460550.19582*x32*x35 - x32*x39 + 439283.26905406616*x32*x54 - x32*x55 - x33*x45 - 308947381650.58746*x35*x37 - x37*x39 - x37*x55 + 48809.252117118471*x45*x64 - 48809.252117118464*x45*x71 + x48*x50 + x48*x52 - x48*x57 - x48*x58 + x49*x65*(-x47*x75 + x47*x78 + x48*x94 + x74*x90 - 23.443292245560333*x76*x92 - x77*x93 - 3.0000000000000004*x92*(x75 - x78 - x79) + x94*x95 + x94*x96 + 46.886584491120665*x92*exp(x18)/((x60)*(x60)*(x60))) - x49*x72*(-x47*x83 + x47*x86 + x48*x98 + x82*x90 - 2109896.3021004298*x84*x97 - x85*x93 - 2.9999999999999996*x92*(x83 - x86 - x87) + x95*x98 + x96*x98 + 4219792.6042008596*x97*exp(x36)/((x67)*(x67)*(x67))) + x50*x51 + x51*x52 - x51*x57 - x51*x58 - x64*x73 + x71*x73 + x80*x81 - x81*x88 - 205964921100.39163*x35*exp(-2515.8872999999999*x28)/((x31)*(x31)*(x31)) + 2288499.1233376856*x14*exp(-8.3862909999999999*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*3.6177;
static const double Vmax = 1.15*3.6177;
static double V = 0.9*3.6177;

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



const char *MgCaFerrite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *MgCaFerrite_slb_em_coder_calib_name(void) {
    return "MgCaFerrite_slb_em";
}

const char *MgCaFerrite_slb_em_coder_calib_formula(void) {
    return "MgAl2O4";
}

const double MgCaFerrite_slb_em_coder_calib_mw(void) {
    return 142.26568;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,4.0,0.0,0.0,0.0,
        1.0,2.0,0.0,0.0,0.0,0.0,
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

const double *MgCaFerrite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double MgCaFerrite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double MgCaFerrite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double MgCaFerrite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double MgCaFerrite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double MgCaFerrite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double MgCaFerrite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double MgCaFerrite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double MgCaFerrite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double MgCaFerrite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double MgCaFerrite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double MgCaFerrite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double MgCaFerrite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double MgCaFerrite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int MgCaFerrite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **MgCaFerrite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **MgCaFerrite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void MgCaFerrite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int MgCaFerrite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double MgCaFerrite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int MgCaFerrite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double MgCaFerrite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double MgCaFerrite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

