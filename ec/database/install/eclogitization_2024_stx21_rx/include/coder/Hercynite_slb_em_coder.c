
static char *identifier = "Hercynite_slb_em.emml:3c9e7bbaa007a1612579b21f2f3502bbc8703c8d:Tue Aug 22 19:33:41 2023";



#include <math.h>

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
    double x3 = sqrt(103.63546402744704*x1 - 259.42629510728*x2 - 8.8380791260999985);
    double x4 = 2.6473728000000003*x3;
    double x5 = 794.21184000000005*x3/T;

    result += 698.41485992487208*T*log(1 - exp(-x5)) - 232.80495330829069*T*Debye(x5) - 53.526445436687659*T - 262166826.43358806*x1 + 95717796.575061202*x2 - 209524.45797746163*log(1 - exp(-x4)) + 69841.485992487214*Debye(x4) + 19019872.603473764 + 3211799223.8258071/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-259.42629510728*pow(x0, 4.0/3.0) + 103.63546402744704*pow(x0, 2.0/3.0) - 8.8380791260999985);
    double x2 = x1/T;
    double x3 = 794.21184000000005*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -554689.35098427499*x2*x5/x6 + 184896.45032809165*x2*(-0.0037773297361066785*T*x4/x1 + 3/(exp(x3) - 1)) - 232.80495330829069*x4 + 698.41485992487208*log(x6) - 53.526445436687659;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(103.63546402744704*x1 - 259.42629510728*x3 - 8.8380791260999985);
    double x6 = 2.6473728000000003*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-34.545154675815681*x2 + 172.95086340485332*x4);
    double x10 = 554689.35098427499*x9;
    double x11 = 794.21184000000005*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 174777884.28905869*x2 - 127623728.76674826*x4 + 184896.45032809168*x9*(-1.1331989208320037*x8*Debye(x6) + 3/(exp(x6) - 1)) - 184896.45032809165*x9*(-0.0037773297361066785*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 6423598447.6516142/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(103.63546402744704*x2 - 259.42629510728*x3 - 8.8380791260999985);
    double x5 = x0*x4;
    double x6 = 794.21184000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-45655655420.4263*x2 + 114287880578.01273*x3 + 3893534891.2300706);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1588.4236800000001*x5)/((x8)*(x8)) - 184896.45032809165*x4*(x0*(0.011331989208320034*T*x11 - 8.9999999999999982/x13) + 0.0037773297361066785*x11 - 2382.6355200000003*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 103.63546402744704*x2;
    double x4 = 259.42629510728*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 8.8380791260999985;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 794.21184000000005*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 440540850.07362688*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(172.95086340485332*x2 - 34.545154675815681)*(-184896.45032809165*x6*(2382.6355200000003*x0*x13*x14/((x15)*(x15)) - 0.0037773297361066785*x12/pow(x5, 3.0/2.0) + (0.011331989208320034*x12*x13 - 8.9999999999999982/x15)/(-x3 + x4 + 8.8380791260999985)) + x11*x9/x10 + x11*exp(-1588.4236800000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 103.63546402744704*x2;
    double x5 = 259.42629510728*x3;
    double x6 = x4 - x5 - 8.8380791260999985;
    double x7 = sqrt(x6);
    double x8 = 2.6473728000000003*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x12*x2*(403.55201461132441*x2 - 57.575257793026132);
    double x14 = 554689.35098427499*x13;
    double x15 = 1.0/(-x4 + x5 + 8.8380791260999985);
    double x16 = x3*((172.95086340485332*x2 - 34.545154675815681)*(172.95086340485332*x2 - 34.545154675815681));
    double x17 = x15*x16;
    double x18 = 1468469.5002454231*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 554689.35098427499*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 794.21184000000005*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 440540850.07362688*x17*x22;
    double x29 = exp(x8);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x8);
    double x33 = x12*x32;
    double x34 = -554689.35098427511*x31 + 209524.45797746166*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x12*x38;
    double x40 = -554689.35098427499*x37 + 698.41485992487208*x39;
    double x41 = x12*x16;

    result += x0*(19270795342.954842*x0 + x11*x14 - x11*x18 + x11*x21 + x13*x34 - x13*x40 - x14*x27 - 291296473.81509781*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 + 297788700.45574594*x3 - 184896.45032809168*x41*(7.9421184000000009*x12*x29/((x30)*(x30)) + x15*(-9.0000000000000018*x31 + 3.3995967624960115*x33) - 1.1331989208320037*x19*x32) + 184896.45032809165*x41*(2382.6355200000003*x12*x22*x35/((x36)*(x36)) + x15*(-8.9999999999999982*x37 + 0.011331989208320034*x39) - 0.0037773297361066785*x19*x38) + x28*exp(-1588.4236800000001*x23)/((x26)*(x26)) - x18*exp(-5.2947456000000006*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-136966966261.27888*x2 + 342863641734.03815*x3 + 11680604673.69021);
    double x5 = 1.0/T;
    double x6 = 103.63546402744704*x2;
    double x7 = 259.42629510728*x3;
    double x8 = x6 - x7 - 8.8380791260999985;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 794.21184000000005*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1588.4236800000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.011331989208320034*x19;
    double x25 = x5*(T*x24 - 8.9999999999999982/x21);
    double x26 = 184896.45032809165*x9;
    double x27 = x17*(-x6 + x7 + 8.8380791260999985);

    result += x0*(-1049648277396.418*x15*x18 - x15*x4 - 349882759132.13934*x16*x18 - x16*x4 + x26*(0.0037773297361066785*x19 - 2382.6355200000003*x23 + x25) - x26*(-1892317.3403885572*x22*x27 - 2382.6355199999989*x23 + x24 + 2.9999999999999996*x25 + 3784634.6807771144*x27*exp(x14)/((x21)*(x21)*(x21))) - 699765518264.27869*x18*exp(-2382.6355200000003*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(152383840770.68362*x2 - 30437103613.617531);
    double x5 = 103.63546402744704*x2;
    double x6 = 259.42629510728*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 8.8380791260999985;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 794.21184000000005*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1588.4236800000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 172.95086340485332*x2 - 34.545154675815681;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0037773297361066785*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2382.6355200000003*x3;
    double x26 = 0.011331989208320034*T*x20;
    double x27 = x19*x26 - 8.9999999999999982/x23;
    double x28 = x0*x27;
    double x29 = 184896.45032809165*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 8.8380791260999985);

    result += x0*x1*x2*(1049648277396.418*x14*x18 - x14*x4 + 349882759132.13934*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7147.9065599999994*x0*x31 - x26*x30 + 2.9999999999999996*x27*x32) + 1892317.3403885572*x16*x24 - 3784634.6807771144*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 699765518264.27869*x18*exp(-2382.6355200000003*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(177781147565.79755*x2 - 25364253011.347942);
    double x4 = 103.63546402744704*x2;
    double x5 = 259.42629510728*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 8.8380791260999985;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 794.21184000000005*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1588.4236800000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 172.95086340485332*x2 - 34.545154675815681;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0037773297361066785*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2382.6355200000003*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 8.8380791260999985;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-8.9999999999999982*x32 + 0.011331989208320034*x33));
    double x35 = 403.55201461132441*x2 - 57.575257793026132;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0037773297361066785*x33;
    double x40 = 2.9999999999999996*x31;
    double x41 = 2.9999999999999996*x39/((x30)*(x30));

    result += -x38*(1049648277396.418*x13*x20 + x13*x3 + 349882759132.13934*x14*x20 + x14*x3 + 184896.45032809165*x15*x34 + 184896.45032809165*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(518.85259021456*x2 - 103.63546402744704)/pow(x6, 5.0/2.0) + x24*x35 - 1892317.3403885572*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(345.90172680970664*x2 - 69.090309351631362) + 3784634.6807771144*x37*exp(x12)/((x26)*(x26)*(x26))) + 699765518264.27869*x20*exp(-2382.6355200000003*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 103.63546402744704*x2;
    double x5 = 259.42629510728*x3;
    double x6 = x4 - x5 - 8.8380791260999985;
    double x7 = sqrt(x6);
    double x8 = 2.6473728000000003*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 172.95086340485332*x2 - 34.545154675815681;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.2947456000000006*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 8.8380791260999985;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1468469.5002454231*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = x2*(1345.1733820377481*x2 - 153.534020781403);
    double x24 = 554689.35098427499*x22*x23;
    double x25 = 1.0/T;
    double x26 = x25*x7;
    double x27 = 794.21184000000005*x26;
    double x28 = exp(-x27);
    double x29 = -x28 + 1;
    double x30 = x28/x29;
    double x31 = pow(T, -2);
    double x32 = x14*x31;
    double x33 = 1588.4236800000001*x26;
    double x34 = exp(-x33)/((x29)*(x29));
    double x35 = 440540850.07362688*x25;
    double x36 = x19*x35;
    double x37 = x18*(345.90172680970664*x2 - 69.090309351631362);
    double x38 = ((x12)*(x12));
    double x39 = x0*x38;
    double x40 = x37*x39;
    double x41 = 1468469.5002454231*x40;
    double x42 = (518.85259021456*x2 - 103.63546402744704)/pow(x6, 5.0/2.0);
    double x43 = x39*x42;
    double x44 = 554689.35098427499*x43;
    double x45 = 1.0/x17;
    double x46 = 403.55201461132441*x2 - 57.575257793026132;
    double x47 = x45*x46;
    double x48 = x12*x3;
    double x49 = 1468469.5002454231*x48;
    double x50 = x16*x49;
    double x51 = x45*(807.10402922264882*x2 - 115.15051558605226);
    double x52 = x46*x48;
    double x53 = x11*x52;
    double x54 = 1664068.0529528251*x53;
    double x55 = x21*x49;
    double x56 = x35*x40;
    double x57 = x35*x48;
    double x58 = x34*x57;
    double x59 = x30*x57;
    double x60 = exp(x8);
    double x61 = x60 - 1;
    double x62 = 1.0/x61;
    double x63 = Debye(x8);
    double x64 = x22*x63;
    double x65 = -3*x62 + 1.1331989208320037*x64;
    double x66 = 184896.45032809168*x22;
    double x67 = exp(x27);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x27);
    double x71 = T*x22*x70;
    double x72 = -3*x69 + 0.0037773297361066785*x71;
    double x73 = 184896.45032809165*x22;
    double x74 = 1.1331989208320037*x63;
    double x75 = x11*x74;
    double x76 = x60/((x61)*(x61));
    double x77 = 7.9421184000000009*x76;
    double x78 = x22*x77;
    double x79 = x45*(-9.0000000000000018*x62 + 3.3995967624960115*x64) - x75 + x78;
    double x80 = 369792.90065618337*x79;
    double x81 = x22*x52;
    double x82 = 0.0037773297361066785*T*x70;
    double x83 = x11*x82;
    double x84 = x67/((x68)*(x68));
    double x85 = 2382.6355200000003*x25*x84;
    double x86 = x22*x85;
    double x87 = x45*(-8.9999999999999982*x69 + 0.011331989208320034*x71) - x83 + x86;
    double x88 = 369792.90065618331*x87;
    double x89 = x12*x2;
    double x90 = x42*x89;
    double x91 = x2*x38;
    double x92 = x45*x91;
    double x93 = x11*x91;
    double x94 = 3.0000000000000004*x65;
    double x95 = x18*x91;
    double x96 = x37*x89;
    double x97 = x31*x92;
    double x98 = 2.9999999999999996*x72;

    result += (-77083181371.819366*x0 - 11662758.63773798*x14*x16 - 3887586.2125793267*x14*x21 + x14*x80 - x14*x88 - x16*x20 - x16*x41 + 776790596.84026074*x2 - x20*x21 - x21*x24 - x21*x41 - x21*x44 - x21*x54 - x23*x65*x66 + x23*x72*x73 + x24*x30 - 992629001.51915312*x3 + 349882759132.13934*x30*x32 + x30*x36 + x30*x44 + x30*x54 + x30*x56 + 1049648277396.418*x32*x34 + x34*x36 + x34*x56 - 184896.45032809168*x43*x65 + 184896.45032809165*x43*x72 + x47*x50 + x47*x55 - x47*x58 - x47*x59 - x48*x66*(x46*x75 - x46*x78 - x47*x94 + x74*x90 - 21.025748226539523*x76*x92 - x77*x93 + 3.0000000000000004*x79*x92 + x94*x95 + x94*x96 + 42.051496453079046*x92*exp(x15)/((x61)*(x61)*(x61))) + x48*x73*(x46*x83 - x46*x86 - x47*x98 + x82*x90 - 1892317.3403885572*x84*x97 - x85*x93 + 2.9999999999999996*x87*x92 + x95*x98 + x96*x98 + 3784634.6807771144*x97*exp(x33)/((x68)*(x68)*(x68))) + x50*x51 + x51*x55 - x51*x58 - x51*x59 - 554689.35098427511*x53*x65 + 554689.35098427499*x53*x72 + x80*x81 - x81*x88 + 699765518264.27869*x32*exp(-2382.6355200000003*x26)/((x29)*(x29)*(x29)) - 7775172.4251586534*x14*exp(-7.9421184000000009*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*16.3372;
static const double Vmax = 1.15*16.3372;
static double V = 0.9*16.3372;

static void coder_solve_V(double T, double P) {
    // Newtsafe routine, newton + bisection with bracket check
    // Update if *either* T or P changes (was only if both changed)
    if ((T != Told) || (P != Pold)) {
        // check bracket
        double fa = -coder_dadv(T, Vmin) - P;
        double fb = -coder_dadv(T, Vmax) - P;
        if ( isnan(fa) ) {
            printf("Error: lower bracket isnan for Vmin=%g\n",Vmin);
        }
        if ( isnan(fb) ) {
            printf("Error: upper bracket isnan for Vmax=%g\n",Vmax);
        }
        if ( isnan(fa) || isnan(fb) || fa*fb > 0.) 
        {
            printf("Error: improper  initial bracket in solve_V\n");
            exit(EXIT_FAILURE); 
        }
        // assert(fa*fb < 0.);
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
            printf("Error: max iterations exceeded in solve_V\n");
            exit(EXIT_FAILURE);
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



const char *Hercynite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Hercynite_slb_em_coder_calib_name(void) {
    return "Hercynite_slb_em";
}

const char *Hercynite_slb_em_coder_calib_formula(void) {
    return "Fe4Al8O16";
}

const double Hercynite_slb_em_coder_calib_mw(void) {
    return 695.23072;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,16.0,0.0,0.0,0.0,
        0.0,8.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,4.0,0.0,0.0,0.0,
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

const double *Hercynite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Hercynite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Hercynite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Hercynite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Hercynite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Hercynite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Hercynite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Hercynite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Hercynite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Hercynite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Hercynite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Hercynite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Hercynite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Hercynite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Hercynite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Hercynite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Hercynite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Hercynite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Hercynite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Hercynite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Hercynite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Hercynite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Hercynite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Hercynite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Hercynite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Hercynite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Hercynite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Hercynite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Hercynite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Hercynite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Hercynite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

