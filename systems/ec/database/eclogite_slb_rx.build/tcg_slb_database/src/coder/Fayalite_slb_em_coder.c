
static char *identifier = "Fayalite_slb_em.emml:346cf78f2814cadf31d6c1c950b3ee277761b771:Thu Feb 10 16:52:01 2022";



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
    double x3 = sqrt(37.892806648613799*x1 - 40.355170896226504*x2 - 7.4116686774499998);
    double x4 = 2.0623356666666668*x3;
    double x5 = 618.70069999999998*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 26.76322271834383*T - 32660399.211665526*x1 + 36494034.179121077*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 4891859.3800706053 + 16412434.93426137/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-40.355170896226504*pow(x0, 4.0/3.0) + 37.892806648613799*pow(x0, 2.0/3.0) - 7.4116686774499998);
    double x2 = x1/T;
    double x3 = 618.70069999999998*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -108027.44068148007*x2*x5/x6 + 36009.146893826692*x2*(-0.0048488711908682182*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(37.892806648613799*x1 - 40.355170896226504*x3 - 7.4116686774499998);
    double x6 = 2.0623356666666668*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-12.630935549537933*x2 + 26.903447264151001*x4);
    double x10 = 618.70069999999998*x5/T;
    double x11 = exp(-x10);

    result += 108027.44068148007*x11*x9/(-x11 + 1) + 21773599.474443682*x2 - 48658712.2388281*x4 - 108027.44068148009*x7*x9/(-x7 + 1) + 36009.146893826699*x9*(-1.4546613572604652*x8*Debye(x6) + 3/(exp(x6) - 1)) - 36009.146893826692*x9*(-0.0048488711908682182*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 32824869.868522741/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(37.892806648613799*x2 - 40.355170896226504*x3 - 7.4116686774499998);
    double x5 = x0*x4;
    double x6 = 618.70069999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2532628375.5673223*x2 + 2697204560.760365*x3 + 495371128.79708219);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1237.4014*x5)/((x8)*(x8)) - 36009.146893826692*x4*(x0*(0.014546613572604655*T*x11 - 9.0/x13) + 0.0048488711908682182*x11 - 1856.1021000000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 37.892806648613799*x2;
    double x4 = 40.355170896226504*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 7.4116686774499998;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 618.70069999999998*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 66836653.1688402*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(26.903447264151001*x2 - 12.630935549537933)*(-36009.146893826692*x6*(1856.1021000000001*x0*x13*x14/((x15)*(x15)) - 0.0048488711908682182*x12/pow(x5, 3.0/2.0) + (0.014546613572604655*x12*x13 - 9.0/x15)/(-x3 + x4 + 7.4116686774499998)) + x11*x9/x10 + x11*exp(-1237.4014*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 37.892806648613799*x2;
    double x5 = 40.355170896226504*x3;
    double x6 = x4 - x5 - 7.4116686774499998;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(62.774710283019004*x2 - 21.051559249229889);
    double x10 = x8*x9;
    double x11 = 2.0623356666666668*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 108027.44068148009*x14;
    double x16 = 1.0/(-x4 + x5 + 7.4116686774499998);
    double x17 = x3*((26.903447264151001*x2 - 12.630935549537933)*(26.903447264151001*x2 - 12.630935549537933));
    double x18 = x16*x17;
    double x19 = 222788.84389613403*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 618.70069999999998*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 108027.44068148007*x27;
    double x29 = 66836653.1688402*x18*x22;
    double x30 = exp(x11);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x11);
    double x34 = x33*x8;
    double x35 = -108027.4406814801*x32 + 52381.114494365407*x34;
    double x36 = exp(x24);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = T*Debye(x24);
    double x40 = x39*x8;
    double x41 = -3*x38 + 0.0048488711908682182*x40;
    double x42 = 36009.146893826692*x8;

    result += x0*(98474609.60556823*x0 + x10*x15 - x10*x28 + x10*x35 - x14*x19 + x15*x21 + x17*x42*(x16*(-9.0*x38 + 0.014546613572604655*x40) - 0.0048488711908682182*x20*x39 + 1856.1021000000001*x22*x36*x8/((x37)*(x37))) - 36009.146893826699*x17*x8*(x16*(-9.0*x32 + 4.3639840717813954*x34) - 1.4546613572604652*x20*x33 + 6.1870070000000004*x30*x8/((x31)*(x31))) - 36289332.457406133*x2 - x21*x28 + x21*x35 - 36009.146893826692*x21*x41 + x27*x29 + 113536995.22393224*x3 - x41*x42*x9 + x29*exp(-1237.4014*x23)/((x26)*(x26)) - x19*exp(-4.1246713333333336*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-7597885126.7019672*x2 + 8091613682.2810946*x3 + 1486113386.3912466);
    double x5 = 1.0/T;
    double x6 = 37.892806648613799*x2;
    double x7 = 40.355170896226504*x3;
    double x8 = x6 - x7 - 7.4116686774499998;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 618.70069999999998*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1237.4014*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -1856.1021000000001*x0*x22*x9;
    double x24 = 0.014546613572604655*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 36009.146893826692*x9;
    double x27 = x17*(-x6 + x7 + 7.4116686774499998);

    result += x0*(-124055652303.65594*x15*x18 - x15*x4 - 41351884101.218651*x16*x18 - x16*x4 + x26*(0.0048488711908682182*x19 + x23 + x25) - x26*(-1148371.6685414701*x22*x27 + x23 + x24 + 3.0*x25 + 2296743.3370829402*x27*exp(x14)/((x21)*(x21)*(x21))) - 82703768202.437302*x18*exp(-1856.1021000000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(3596272747.6804867*x2 - 1688418917.0448816);
    double x5 = 37.892806648613799*x2;
    double x6 = 40.355170896226504*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 7.4116686774499998;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 618.70069999999998*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1237.4014*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 26.903447264151001*x2 - 12.630935549537933;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0048488711908682182*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 1856.1021000000001*x3;
    double x26 = 0.014546613572604655*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 36009.146893826692*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 7.4116686774499998);

    result += x0*x1*x2*(124055652303.65594*x14*x18 - x14*x4 + 41351884101.218651*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(5568.3063000000002*x0*x31 - x26*x30 + 3.0*x27*x32) + 1148371.6685414701*x16*x24 - 2296743.3370829402*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 82703768202.437302*x18*exp(-1856.1021000000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(4195651538.9605675*x2 - 1407015764.2040682);
    double x4 = 37.892806648613799*x2;
    double x5 = 40.355170896226504*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 7.4116686774499998;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 618.70069999999998*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1237.4014*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 26.903447264151001*x2 - 12.630935549537933;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0048488711908682182*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 1856.1021000000001*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 7.4116686774499998;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.014546613572604655*x33));
    double x35 = 62.774710283019004*x2 - 21.051559249229889;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0048488711908682182*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(124055652303.65594*x13*x20 + x13*x3 + 41351884101.218651*x14*x20 + x14*x3 + 36009.146893826692*x15*x34 + 36009.146893826692*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(80.710341792453008*x2 - 37.892806648613799)/pow(x6, 5.0/2.0) + x24*x35 - 1148371.6685414701*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(53.806894528302003*x2 - 25.261871099075865) + 2296743.3370829402*x37*exp(x12)/((x26)*(x26)*(x26))) + 82703768202.437302*x20*exp(-1856.1021000000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 37.892806648613799*x2;
    double x5 = 40.355170896226504*x3;
    double x6 = x4 - x5 - 7.4116686774499998;
    double x7 = sqrt(x6);
    double x8 = 2.0623356666666668*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 26.903447264151001*x2 - 12.630935549537933;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.1246713333333336*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 7.4116686774499998;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 222788.84389613403*x16;
    double x21 = x9/x10;
    double x22 = 222788.84389613403*x21;
    double x23 = 1.0/x7;
    double x24 = 209.24903427673001*x2 - 56.137491331279705;
    double x25 = x2*x23*x24;
    double x26 = 108027.44068148009*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 618.70069999999998*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 108027.44068148007*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1237.4014*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 66836653.1688402*x27;
    double x39 = x19*x38;
    double x40 = 53.806894528302003*x2 - 25.261871099075865;
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x18*x40*x42;
    double x44 = (80.710341792453008*x2 - 37.892806648613799)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 62.774710283019004*x2 - 21.051559249229889;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = x20*x49;
    double x51 = x47*(125.54942056603801*x2 - 42.103118498459779);
    double x52 = x11*x47;
    double x53 = x46*x52;
    double x54 = x22*x49;
    double x55 = x38*x43;
    double x56 = x38*x49;
    double x57 = x37*x56;
    double x58 = x32*x56;
    double x59 = 36009.146893826699*x23;
    double x60 = exp(x8);
    double x61 = x60 - 1;
    double x62 = 1.0/x61;
    double x63 = Debye(x8);
    double x64 = x23*x63;
    double x65 = -3*x62 + 1.4546613572604652*x64;
    double x66 = x2*x65;
    double x67 = exp(x29);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x29);
    double x71 = T*x23*x70;
    double x72 = -3*x69 + 0.0048488711908682182*x71;
    double x73 = x2*x72;
    double x74 = 36009.146893826692*x23;
    double x75 = x46*x65;
    double x76 = x46*x72;
    double x77 = 1.4546613572604652*x63;
    double x78 = x11*x77;
    double x79 = x60/((x61)*(x61));
    double x80 = 6.1870070000000004*x79;
    double x81 = x23*x80;
    double x82 = x49*(-9.0*x62 + 4.3639840717813954*x64) - x78 + x81;
    double x83 = 72018.293787653398*x82;
    double x84 = x23*x48;
    double x85 = 0.0048488711908682182*T*x70;
    double x86 = x11*x85;
    double x87 = x67/((x68)*(x68));
    double x88 = 1856.1021000000001*x27*x87;
    double x89 = x23*x88;
    double x90 = x49*(-9.0*x69 + 0.014546613572604655*x71) - x86 + x89;
    double x91 = 72018.293787653383*x90;
    double x92 = x12*x2*x44;
    double x93 = x2*x41;
    double x94 = x49*x93;
    double x95 = x11*x93;
    double x96 = 3.0*x49;
    double x97 = 3.0*x18;
    double x98 = x12*x40;
    double x99 = x93*x96;
    double x100 = x34*x94;
    double x101 = x73*x97;

    result += (-393898438.42227292*x0 - 1378396.1367072887*x14*x16 - 459465.37890242954*x14*x21 + x14*x83 - x14*x91 - x19*x20 - x19*x22 + 96771553.219749689*x2 - x20*x43 - 324082.32204444025*x21*x53 - x22*x43 - x24*x59*x66 + x24*x73*x74 - x25*x26 + x25*x33 - x26*x45 - 378456650.74644077*x3 + 41351884101.218651*x32*x35 + x32*x39 + 324082.3220444402*x32*x53 + x32*x55 + x33*x45 + 124055652303.65594*x35*x37 + x37*x39 + x37*x55 - 36009.146893826699*x45*x65 + 36009.146893826692*x45*x72 - x47*x59*(x46*x78 - x46*x81 + x65*x93*x97 + x66*x97*x98 - x75*x96 + x77*x92 - 12.759685206016336*x79*x94 - x80*x95 + x82*x99 + 25.519370412032671*x94*exp(x15)/((x61)*(x61)*(x61))) + x47*x74*(-1148371.6685414701*x100*x87 + 2296743.3370829402*x100*exp(x36)/((x68)*(x68)*(x68)) + x101*x41 + x101*x98 + x46*x86 - x46*x89 - x76*x96 + x85*x92 - x88*x95 + x90*x99) + x48*x50 + x48*x54 - x48*x57 - x48*x58 + x50*x51 + x51*x54 - x51*x57 - x51*x58 - 108027.4406814801*x52*x75 + 108027.44068148007*x52*x76 + x83*x84 - x84*x91 + 82703768202.437302*x35*exp(-1856.1021000000001*x28)/((x31)*(x31)*(x31)) - 918930.75780485908*x14*exp(-6.1870070000000004*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*4.629;
static const double Vmax = 1.15*4.629;
static double V = 0.9*4.629;

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



const char *Fayalite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Fayalite_slb_em_coder_calib_name(void) {
    return "Fayalite_slb_em";
}

const char *Fayalite_slb_em_coder_calib_formula(void) {
    return "Fe2SiO4";
}

const double Fayalite_slb_em_coder_calib_mw(void) {
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

const double *Fayalite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Fayalite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Fayalite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Fayalite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Fayalite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Fayalite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Fayalite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Fayalite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Fayalite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Fayalite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Fayalite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Fayalite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Fayalite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Fayalite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Fayalite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Fayalite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Fayalite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Fayalite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Fayalite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Fayalite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Fayalite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Fayalite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Fayalite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Fayalite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Fayalite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Fayalite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Fayalite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Fayalite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Fayalite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Fayalite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Fayalite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

