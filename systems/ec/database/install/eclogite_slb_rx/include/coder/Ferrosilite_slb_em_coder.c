
static char *identifier = "Ferrosilite_slb_em.emml:224b1e32b1db209aba4d4377cb41276ee9650c1e:Thu Feb 10 16:52:24 2022";



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
    double x3 = sqrt(37.926780180236925*x1 - 53.349677172598199*x2 - 5.4711099017000011);
    double x4 = 2.2482563333333334*x3;
    double x5 = 674.4769*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) - 26.76322271834383*T + 66649136.311499238*x1 - 326592858.44354093*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) - 6057466.061557347 + 490916391.62581241/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-53.349677172598199*pow(x0, 4.0/3.0) + 37.926780180236925*pow(x0, 2.0/3.0) - 5.4711099017000011);
    double x2 = x1/T;
    double x3 = 674.4769*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -168237.38915573643*x2*x5/x6 + 56079.129718578813*x2*(-0.0044478913955392687*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(37.926780180236925*x1 - 53.349677172598199*x3 - 5.4711099017000011);
    double x6 = 2.2482563333333334*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-12.642260060078975*x2 + 35.566451448398794*x4);
    double x10 = 168237.38915573643*x9;
    double x11 = 674.4769*x5/T;
    double x12 = exp(-x11);
    double x13 = 56079.129718578813*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.3343674186617807*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0044478913955392687*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 44432757.540999487*x2 + 435457144.59138787*x4 - 981832783.25162482/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(37.926780180236925*x2 - 53.349677172598199*x3 - 5.4711099017000011);
    double x5 = x0*x4;
    double x6 = 674.4769*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-4303636426.2439356*x2 + 6053706982.6978893*x3 + 620819055.90312397);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1348.9538*x5)/((x8)*(x8)) - 56079.129718578813*x4*(x0*(0.013343674186617806*T*x11 - 9.0/x13) + 0.0044478913955392687*x11 - 2023.4306999999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 37.926780180236925*x2;
    double x4 = 53.349677172598199*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 5.4711099017000011;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 674.4769*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 113472232.70185472*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(35.566451448398794*x2 - 12.642260060078975)*(-56079.129718578813*x6*(2023.4306999999999*x0*x13*x14/((x15)*(x15)) - 0.0044478913955392687*x12/pow(x5, 3.0/2.0) + (0.013343674186617806*x12*x13 - 9.0/x15)/(-x3 + x4 + 5.4711099017000011)) + x11*x9/x10 + x11*exp(-1348.9538*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 37.926780180236925*x2;
    double x5 = 53.349677172598199*x3;
    double x6 = x4 - x5 - 5.4711099017000011;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(82.98838671293052*x2 - 21.070433433464956);
    double x10 = 168237.38915573643*x9;
    double x11 = 2.2482563333333334*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 1.0/(-x4 + x5 + 5.4711099017000011);
    double x16 = x3*((35.566451448398794*x2 - 12.642260060078975)*(35.566451448398794*x2 - 12.642260060078975));
    double x17 = x15*x16;
    double x18 = 378240.77567284909*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 168237.38915573643*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 674.4769*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 113472232.70185472*x17*x22;
    double x29 = exp(x11);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x11);
    double x33 = x32*x8;
    double x34 = -168237.38915573643*x31 + 74830.163563379159*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x38*x8;
    double x40 = -168237.38915573643*x37 + 249.43387854459721*x39;
    double x41 = 56079.129718578813*x16*x8;

    result += x0*(2945498349.7548742*x0 + x10*x14 - x10*x27 - x14*x18 + x14*x21 + 74054595.901665807*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 - 1016066670.7132384*x3 + x34*x9 - x40*x9 - x41*(x15*(-9.0*x31 + 4.0031022559853415*x33) - 1.3343674186617807*x19*x32 + 6.7447689999999998*x29*x8/((x30)*(x30))) + x41*(x15*(-9.0*x37 + 0.013343674186617806*x39) - 0.0044478913955392687*x19*x38 + 2023.4306999999999*x22*x35*x8/((x36)*(x36))) + x28*exp(-1348.9538*x23)/((x26)*(x26)) - x18*exp(-4.4965126666666668*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-12910909278.731808*x2 + 18161120948.09367*x3 + 1862457167.7093723);
    double x5 = 1.0/T;
    double x6 = 37.926780180236925*x2;
    double x7 = 53.349677172598199*x3;
    double x8 = x6 - x7 - 5.4711099017000011;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 674.4769*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1348.9538*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2023.4306999999999*x0*x22*x9;
    double x24 = 0.013343674186617806*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 56079.129718578813*x9;
    double x27 = x17*(-x6 + x7 + 5.4711099017000011);

    result += x0*(-229603199246.47678*x15*x18 - x15*x4 - 76534399748.825592*x16*x18 - x16*x4 + x26*(0.0044478913955392687*x19 + x23 + x25) - x26*(-1364757.26590083*x22*x27 + x23 + x24 + 3.0*x25 + 2729514.5318016601*x27*exp(x14)/((x21)*(x21)*(x21))) - 153068799497.65118*x18*exp(-2023.4306999999999*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(8071609310.2638521*x2 - 2869090950.8292904);
    double x5 = 37.926780180236925*x2;
    double x6 = 53.349677172598199*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 5.4711099017000011;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 674.4769*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1348.9538*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 35.566451448398794*x2 - 12.642260060078975;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0044478913955392687*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2023.4306999999999*x3;
    double x26 = 0.013343674186617806*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 56079.129718578813*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 5.4711099017000011);

    result += x0*x1*x2*(229603199246.47678*x14*x18 - x14*x4 + 76534399748.825592*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(6070.2920999999997*x0*x31 - x26*x30 + 3.0*x27*x32) + 1364757.26590083*x16*x24 - 2729514.5318016601*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 153068799497.65118*x18*exp(-2023.4306999999999*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(9416877528.641161*x2 - 2390909125.6910753);
    double x4 = 37.926780180236925*x2;
    double x5 = 53.349677172598199*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 5.4711099017000011;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 674.4769*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1348.9538*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 35.566451448398794*x2 - 12.642260060078975;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0044478913955392687*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2023.4306999999999*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 5.4711099017000011;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.013343674186617806*x33));
    double x35 = 82.98838671293052*x2 - 21.070433433464956;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0044478913955392687*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(229603199246.47678*x13*x20 + x13*x3 + 76534399748.825592*x14*x20 + x14*x3 + 56079.129718578813*x15*x34 + 56079.129718578813*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(106.69935434519638*x2 - 37.926780180236925)/pow(x6, 5.0/2.0) + x24*x35 - 1364757.26590083*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(71.132902896797589*x2 - 25.28452012015795) + 2729514.5318016601*x37*exp(x12)/((x26)*(x26)*(x26))) + 153068799497.65118*x20*exp(-2023.4306999999999*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 37.926780180236925*x2;
    double x5 = 53.349677172598199*x3;
    double x6 = x4 - x5 - 5.4711099017000011;
    double x7 = sqrt(x6);
    double x8 = 2.2482563333333334*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 35.566451448398794*x2 - 12.642260060078975;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.4965126666666668*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 5.4711099017000011;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 378240.77567284909*x16;
    double x21 = x9/x10;
    double x22 = 378240.77567284909*x21;
    double x23 = 1.0/x7;
    double x24 = 276.6279557097684*x2 - 56.187822489239878;
    double x25 = x2*x23*x24;
    double x26 = 168237.38915573643*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 674.4769*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 168237.38915573643*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1348.9538*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 113472232.70185472*x27;
    double x39 = x19*x38;
    double x40 = 71.132902896797589*x2 - 25.28452012015795;
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x18*x40*x42;
    double x44 = (106.69935434519638*x2 - 37.926780180236925)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 82.98838671293052*x2 - 21.070433433464956;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = x20*x49;
    double x51 = x47*(165.97677342586104*x2 - 42.140866866929912);
    double x52 = x11*x47;
    double x53 = 504712.16746720928*x46*x52;
    double x54 = x22*x49;
    double x55 = x38*x43;
    double x56 = x38*x49;
    double x57 = x37*x56;
    double x58 = x32*x56;
    double x59 = 56079.129718578813*x23;
    double x60 = x24*x59;
    double x61 = exp(x8);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x8);
    double x65 = x23*x64;
    double x66 = -3*x63 + 1.3343674186617807*x65;
    double x67 = x2*x66;
    double x68 = exp(x29);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x29);
    double x72 = T*x23*x71;
    double x73 = -3*x70 + 0.0044478913955392687*x72;
    double x74 = x2*x73;
    double x75 = 56079.129718578813*x45;
    double x76 = x46*x66;
    double x77 = 168237.38915573643*x52;
    double x78 = x46*x73;
    double x79 = 1.3343674186617807*x64;
    double x80 = x11*x79;
    double x81 = x61/((x62)*(x62));
    double x82 = 6.7447689999999998*x81;
    double x83 = x23*x82;
    double x84 = x49*(-9.0*x63 + 4.0031022559853415*x65) - x80 + x83;
    double x85 = 112158.25943715763*x84;
    double x86 = x23*x48;
    double x87 = 0.0044478913955392687*T*x71;
    double x88 = x11*x87;
    double x89 = x68/((x69)*(x69));
    double x90 = 2023.4306999999999*x27*x89;
    double x91 = x23*x90;
    double x92 = x49*(-9.0*x70 + 0.013343674186617806*x72) - x88 + x91;
    double x93 = 112158.25943715763*x92;
    double x94 = x12*x2*x44;
    double x95 = x2*x41;
    double x96 = x49*x95;
    double x97 = x11*x95;
    double x98 = 3.0*x49;
    double x99 = 3.0*x18;
    double x100 = x12*x40;
    double x101 = x95*x98;
    double x102 = x47*x59;
    double x103 = x34*x96;
    double x104 = x74*x99;

    result += (-11781993399.019497*x0 + x102*(x100*x104 + x101*x92 - 1364757.26590083*x103*x89 + 2729514.5318016601*x103*exp(x36)/((x69)*(x69)*(x69)) + x104*x41 + x46*x88 - x46*x91 - x78*x98 + x87*x94 - x90*x97) - x102*(x100*x67*x99 + x101*x84 + x46*x80 - x46*x83 + x66*x95*x99 - x76*x98 + x79*x94 - 15.163969621120334*x81*x96 - x82*x97 + 30.327939242240667*x96*exp(x15)/((x62)*(x62)*(x62))) - 2551146.6582941869*x14*x16 - 850382.21943139564*x14*x21 + x14*x85 - x14*x93 - x19*x20 - x19*x22 - 197478922.40444213*x2 - x20*x43 - x21*x53 - x22*x43 - x25*x26 + x25*x33 - x26*x45 + 3386888902.3774614*x3 + 76534399748.825592*x32*x35 + x32*x39 + x32*x53 + x32*x55 + x33*x45 + 229603199246.47678*x35*x37 + x37*x39 + x37*x55 + x48*x50 + x48*x54 - x48*x57 - x48*x58 + x50*x51 + x51*x54 - x51*x57 - x51*x58 - x60*x67 + x60*x74 - x66*x75 + x73*x75 - x76*x77 + x77*x78 + x85*x86 - x86*x93 + 153068799497.65118*x35*exp(-2023.4306999999999*x28)/((x31)*(x31)*(x31)) - 1700764.4388627913*x14*exp(-6.7447689999999998*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.5941;
static const double Vmax = 1.15*6.5941;
static double V = 0.9*6.5941;

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



const char *Ferrosilite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Ferrosilite_slb_em_coder_calib_name(void) {
    return "Ferrosilite_slb_em";
}

const char *Ferrosilite_slb_em_coder_calib_formula(void) {
    return "Fe2Si2O6";
}

const double Ferrosilite_slb_em_coder_calib_mw(void) {
    return 263.8614;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,6.0,0.0,0.0,0.0,
        0.0,0.0,2.0,0.0,0.0,0.0,
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

const double *Ferrosilite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Ferrosilite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Ferrosilite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Ferrosilite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Ferrosilite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Ferrosilite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Ferrosilite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Ferrosilite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Ferrosilite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Ferrosilite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Ferrosilite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Ferrosilite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Ferrosilite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Ferrosilite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Ferrosilite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Ferrosilite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Ferrosilite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Ferrosilite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Ferrosilite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Ferrosilite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Ferrosilite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Ferrosilite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Ferrosilite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Ferrosilite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Ferrosilite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Ferrosilite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Ferrosilite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Ferrosilite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Ferrosilite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

