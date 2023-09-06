
static char *identifier = "MgMajorite_slb_em.emml:b7bc07aee8a257a14dc3e591fa640f25ac93af95:Fri May 26 02:49:51 2023";



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
    double x3 = sqrt(20.422879221122948*x1 - 14.086132733762462*x2 - 2.4773727812750002);
    double x4 = 2.7416187333333335*x3;
    double x5 = 822.48562000000004*x3/T;

    result += 498.86775708919441*T*log(1 - exp(-x5)) - 166.2892523630648*T*Debye(x5) - 181301791.31923041*x1 + 373148702.01547527*x2 - 149660.32712675832*log(1 - exp(-x4)) + 49886.775708919442*Debye(x4) + 13294048.561349168 + 293979887.61160558/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-14.086132733762462*pow(x0, 4.0/3.0) + 20.422879221122948*pow(x0, 2.0/3.0) - 2.4773727812750002);
    double x2 = x1/T;
    double x3 = 822.48562000000004*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -410311.55648751551*x2*x5/x6 + 136770.51882917184*x2*(-0.003647480183300955*T*x4/x1 + 3/(exp(x3) - 1)) - 166.2892523630648*x4 + 498.86775708919441*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(20.422879221122948*x1 - 14.086132733762462*x3 - 2.4773727812750002);
    double x6 = 2.7416187333333335*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-6.8076264070409822*x2 + 9.3907551558416404*x4);
    double x10 = 822.48562000000004*x5/T;
    double x11 = exp(-x10);
    double x12 = 136770.51882917184*x9;

    result += 410311.55648751551*x11*x9/(-x11 + 1) + x12*(-1.0942440549902865*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.003647480183300955*T*x8*Debye(x10) + 3/(exp(x10) - 1)) + 120867860.87948693*x2 - 497531602.68730032*x4 - 410311.55648751545*x7*x9/(-x7 + 1) - 587959775.22321117/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(20.422879221122948*x2 - 14.086132733762462*x3 - 2.4773727812750002);
    double x5 = x0*x4;
    double x6 = 822.48562000000004*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-6892218413.8573112*x2 + 4753722643.9288359*x3 + 836052258.65668201);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1644.9712400000001*x5)/((x8)*(x8)) - 136770.51882917184*x4*(x0*(0.010942440549902864*T*x11 - 8.9999999999999982/x13) + 0.003647480183300955*x11 - 2467.4568600000002*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 20.422879221122948*x2;
    double x4 = 14.086132733762462*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 2.4773727812750002;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 822.48562000000004*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 337475354.93079925*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(9.3907551558416404*x2 - 6.8076264070409822)*(-136770.51882917184*x6*(2467.4568600000002*x0*x13*x14/((x15)*(x15)) - 0.003647480183300955*x12/pow(x5, 3.0/2.0) + (0.010942440549902864*x12*x13 - 8.9999999999999982/x15)/(-x3 + x4 + 2.4773727812750002)) + x11*x9/x10 + x11*exp(-1644.9712400000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 20.422879221122948*x2;
    double x5 = 14.086132733762462*x3;
    double x6 = x4 - x5 - 2.4773727812750002;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(21.911762030297162*x2 - 11.346044011734971);
    double x10 = 2.7416187333333335*x7;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = x11/x12;
    double x14 = 410311.55648751545*x13;
    double x15 = 1.0/(-x4 + x5 + 2.4773727812750002);
    double x16 = x3*((9.3907551558416404*x2 - 6.8076264070409822)*(9.3907551558416404*x2 - 6.8076264070409822));
    double x17 = x15*x16;
    double x18 = 1124917.8497693306*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 1.0/T;
    double x22 = x21*x7;
    double x23 = 822.48562000000004*x22;
    double x24 = exp(-x23);
    double x25 = -x24 + 1;
    double x26 = x24/x25;
    double x27 = 410311.55648751551*x26;
    double x28 = 337475354.93079925*x17*x21;
    double x29 = exp(x10);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x10);
    double x33 = x32*x8;
    double x34 = -410311.55648751551*x31 + 149660.32712675832*x33;
    double x35 = exp(x23);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x23);
    double x39 = x38*x8;
    double x40 = -410311.55648751551*x37 + 498.86775708919441*x39;
    double x41 = 136770.51882917184*x16*x8;

    result += x0*(1763879325.6696334*x0 - x13*x18 + x14*x20 + x14*x9 - 201446434.79914486*x2 - x20*x27 + x20*x34 - x20*x40 + x26*x28 - x27*x9 + 1160907072.9370341*x3 + x34*x9 - x40*x9 - x41*(x15*(-9.0*x31 + 3.2827321649708594*x33) - 1.0942440549902865*x19*x32 + 8.2248562000000014*x29*x8/((x30)*(x30))) + x41*(x15*(-8.9999999999999982*x37 + 0.010942440549902864*x39) - 0.003647480183300955*x19*x38 + 2467.4568600000002*x21*x35*x8/((x36)*(x36))) + x28*exp(-1644.9712400000001*x22)/((x25)*(x25)) - x18*exp(-5.483237466666667*x7)/((x12)*(x12)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-20676655241.571934*x2 + 14261167931.786509*x3 + 2508156775.970046);
    double x5 = 1.0/T;
    double x6 = 20.422879221122948*x2;
    double x7 = 14.086132733762462*x3;
    double x8 = x6 - x7 - 2.4773727812750002;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 822.48562000000004*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1644.9712400000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.010942440549902864*x19;
    double x25 = x5*(T*x24 - 8.9999999999999982/x21);
    double x26 = 136770.51882917184*x9;
    double x27 = x17*(-x6 + x7 + 2.4773727812750002);

    result += x0*(-832705879604.93555*x15*x18 - x15*x4 - 277568626534.97852*x16*x18 - x16*x4 + x26*(0.003647480183300955*x19 - 2467.4568600000002*x23 + x25) - x26*(-2029447.7853203535*x22*x27 - 2467.4568599999993*x23 + x24 + 2.9999999999999996*x25 + 4058895.5706407069*x27*exp(x14)/((x21)*(x21)*(x21))) - 555137253069.95703*x18*exp(-2467.4568600000002*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(6338296858.5717812*x2 - 4594812275.9048738);
    double x5 = 20.422879221122948*x2;
    double x6 = 14.086132733762462*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 2.4773727812750002;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 822.48562000000004*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1644.9712400000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 9.3907551558416404*x2 - 6.8076264070409822;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.003647480183300955*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2467.4568600000002*x3;
    double x26 = 0.010942440549902864*T*x20;
    double x27 = x19*x26 - 8.9999999999999982/x23;
    double x28 = x0*x27;
    double x29 = 136770.51882917184*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 2.4773727812750002);

    result += x0*x1*x2*(832705879604.93555*x14*x18 - x14*x4 + 277568626534.97852*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7402.3705799999998*x0*x31 - x26*x30 + 2.9999999999999996*x27*x32) + 2029447.7853203535*x16*x24 - 4058895.5706407069*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 555137253069.95703*x18*exp(-2467.4568600000002*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(7394679668.333745*x2 - 3829010229.9207287);
    double x4 = 20.422879221122948*x2;
    double x5 = 14.086132733762462*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 2.4773727812750002;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 822.48562000000004*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1644.9712400000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 9.3907551558416404*x2 - 6.8076264070409822;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.003647480183300955*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2467.4568600000002*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 2.4773727812750002;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-8.9999999999999982*x32 + 0.010942440549902864*x33));
    double x35 = 21.911762030297162*x2 - 11.346044011734971;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.003647480183300955*x33;
    double x40 = 2.9999999999999996*x31;
    double x41 = 2.9999999999999996*x39/((x30)*(x30));

    result += -x38*(832705879604.93555*x13*x20 + x13*x3 + 277568626534.97852*x14*x20 + x14*x3 + 136770.51882917184*x15*x34 + 136770.51882917184*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(28.172265467524923*x2 - 20.422879221122948)/pow(x6, 5.0/2.0) + x24*x35 - 2029447.7853203535*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(18.781510311683281*x2 - 13.615252814081964) + 4058895.5706407069*x37*exp(x12)/((x26)*(x26)*(x26))) + 555137253069.95703*x20*exp(-2467.4568600000002*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 20.422879221122948*x2;
    double x5 = 14.086132733762462*x3;
    double x6 = x4 - x5 - 2.4773727812750002;
    double x7 = sqrt(x6);
    double x8 = 2.7416187333333335*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 9.3907551558416404*x2 - 6.8076264070409822;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.483237466666667*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 2.4773727812750002;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1124917.8497693306*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = x2*(73.039206767657205*x2 - 30.25611736462659);
    double x24 = x22*x23;
    double x25 = 410311.55648751545*x21;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 822.48562000000004*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 410311.55648751551*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1644.9712400000001*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 337475354.93079925*x26;
    double x38 = x19*x37;
    double x39 = x18*(18.781510311683281*x2 - 13.615252814081964);
    double x40 = ((x12)*(x12));
    double x41 = x0*x40;
    double x42 = x39*x41;
    double x43 = 1124917.8497693306*x42;
    double x44 = (28.172265467524923*x2 - 20.422879221122948)/pow(x6, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 1.0/x17;
    double x47 = 21.911762030297162*x2 - 11.346044011734971;
    double x48 = x46*x47;
    double x49 = x12*x3;
    double x50 = 1124917.8497693306*x49;
    double x51 = x16*x50;
    double x52 = x46*(43.823524060594323*x2 - 22.692088023469942);
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
    double x65 = -3*x62 + 1.0942440549902865*x64;
    double x66 = 136770.51882917184*x22;
    double x67 = x23*x66;
    double x68 = exp(x28);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x28);
    double x72 = T*x22*x71;
    double x73 = -3*x70 + 0.003647480183300955*x72;
    double x74 = 136770.51882917184*x45;
    double x75 = 410311.55648751551*x54;
    double x76 = 1.0942440549902865*x63;
    double x77 = x11*x76;
    double x78 = x60/((x61)*(x61));
    double x79 = 8.2248562000000014*x78;
    double x80 = x22*x79;
    double x81 = x46*(-9.0*x62 + 3.2827321649708594*x64) - x77 + x80;
    double x82 = 273541.03765834367*x81;
    double x83 = x22*x53;
    double x84 = 0.003647480183300955*T*x71;
    double x85 = x11*x84;
    double x86 = x68/((x69)*(x69));
    double x87 = 2467.4568600000002*x26*x86;
    double x88 = x22*x87;
    double x89 = x46*(-8.9999999999999982*x70 + 0.010942440549902864*x72) - x85 + x88;
    double x90 = 273541.03765834367*x89;
    double x91 = x12*x2;
    double x92 = x44*x91;
    double x93 = x2*x40;
    double x94 = x46*x93;
    double x95 = x11*x93;
    double x96 = 3.0*x65;
    double x97 = x18*x93;
    double x98 = x39*x91;
    double x99 = x49*x66;
    double x100 = x33*x94;
    double x101 = 2.9999999999999996*x73;

    result += (-7055517302.6785336*x0 - 9252287.5511659477*x14*x16 - 3084095.8503886494*x14*x21 + x14*x82 - x14*x90 - x16*x20 - x16*x43 + 537190492.7977196*x2 - x20*x21 - x21*x43 - 1230934.6694625462*x21*x54 - x24*x25 + x24*x32 - x25*x45 - 3869690243.1234469*x3 + 277568626534.97852*x31*x34 + x31*x38 + 1230934.6694625465*x31*x54 + x31*x56 + x32*x45 + 832705879604.93555*x34*x36 + x36*x38 + x36*x56 + x48*x51 + x48*x55 - x48*x58 - x48*x59 + x51*x52 + x52*x55 - x52*x58 - x52*x59 - x65*x67 - x65*x74 - x65*x75 + x67*x73 + x73*x74 + x73*x75 + x82*x83 - x83*x90 + x99*(-2029447.7853203535*x100*x86 + 4058895.5706407069*x100*exp(x35)/((x69)*(x69)*(x69)) - x101*x48 + x101*x97 + x101*x98 + x47*x85 - x47*x88 + x84*x92 - x87*x95 + 2.9999999999999996*x89*x94) - x99*(x47*x77 - x47*x80 - x48*x96 + x76*x92 - 22.54941983689282*x78*x94 - x79*x95 + 3.0*x81*x94 + x96*x97 + x96*x98 + 45.098839673785641*x94*exp(x15)/((x61)*(x61)*(x61))) + 555137253069.95703*x34*exp(-2467.4568600000002*x27)/((x30)*(x30)*(x30)) - 6168191.7007772988*x14*exp(-8.2248562000000014*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*11.4324;
static const double Vmax = 1.15*11.4324;
static double V = 0.9*11.4324;

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



const char *MgMajorite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *MgMajorite_slb_em_coder_calib_name(void) {
    return "MgMajorite_slb_em";
}

const char *MgMajorite_slb_em_coder_calib_formula(void) {
    return "Mg4Si4O12";
}

const double MgMajorite_slb_em_coder_calib_mw(void) {
    return 401.5548;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,12.0,0.0,0.0,0.0,
        4.0,0.0,4.0,0.0,0.0,0.0,
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

const double *MgMajorite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double MgMajorite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double MgMajorite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double MgMajorite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double MgMajorite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double MgMajorite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double MgMajorite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double MgMajorite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double MgMajorite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double MgMajorite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double MgMajorite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double MgMajorite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double MgMajorite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double MgMajorite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double MgMajorite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double MgMajorite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double MgMajorite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double MgMajorite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double MgMajorite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double MgMajorite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int MgMajorite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **MgMajorite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **MgMajorite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void MgMajorite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int MgMajorite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double MgMajorite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int MgMajorite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double MgMajorite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double MgMajorite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

