
static char *identifier = "Grossular_slb_em.emml:f919d813f8ca5062a1c5ad7dccde44610c5d8fa4:Thu Feb 10 16:52:31 2022";



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
    double x3 = sqrt(28.480921258115991*x1 - 30.824277637309436*x2 - 3.2233037561);
    double x4 = 2.7424766666666671*x3;
    double x5 = 822.74300000000005*x3/T;

    result += 498.86775708919441*T*log(1 - exp(-x5)) - 166.2892523630648*T*Debye(x5) - 269552541.60092473*x1 + 769700951.84614527*x2 - 149660.32712675832*log(1 - exp(-x4)) + 49886.775708919442*Debye(x4) + 18231983.022500161 - 155648845.76175645/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-30.824277637309436*pow(x0, 4.0/3.0) + 28.480921258115991*pow(x0, 2.0/3.0) - 3.2233037561);
    double x2 = x1/T;
    double x3 = 822.74300000000005*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -410439.95507083513*x2*x5/x6 + 136813.31835694503*x2*(-0.0036463391362794944*T*x4/x1 + 3/(exp(x3) - 1)) - 166.2892523630648*x4 + 498.86775708919441*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(28.480921258115991*x1 - 30.824277637309436*x3 - 3.2233037561);
    double x6 = 2.7424766666666671*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-9.4936404193719959*x2 + 20.549518424872957*x4);
    double x10 = 410439.95507083513*x9;
    double x11 = 822.74300000000005*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 179701694.40061647*x2 - 1026267935.7948604*x4 + 136813.31835694506*x9*(-1.0939017408838483*x8*Debye(x6) + 3/(exp(x6) - 1)) - 136813.31835694503*x9*(-0.0036463391362794944*T*x8*Debye(x11) + 3/(exp(x11) - 1)) + 311297691.5235129/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(28.480921258115991*x2 - 30.824277637309436*x3 - 3.2233037561);
    double x5 = x0*x4;
    double x6 = 822.74300000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-9617625463.2348309*x2 + 10408945511.40716*x3 + 1088466486.0190871);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1645.4860000000001*x5)/((x8)*(x8)) - 136813.31835694503*x4*(x0*(0.010939017408838484*T*x11 - 9.0000000000000018/x13) + 0.0036463391362794944*x11 - 2468.2290000000003*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 28.480921258115991*x2;
    double x4 = 30.824277637309436*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 3.2233037561;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 822.74300000000005*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 337686599.95484412*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(20.549518424872957*x2 - 9.4936404193719959)*(-136813.31835694503*x6*(2468.2290000000003*x0*x13*x14/((x15)*(x15)) - 0.0036463391362794944*x12/pow(x5, 3.0/2.0) + (0.010939017408838484*x12*x13 - 9.0000000000000018/x15)/(-x3 + x4 + 3.2233037561)) + x11*x9/x10 + x11*exp(-1645.4860000000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 28.480921258115991*x2;
    double x5 = 30.824277637309436*x3;
    double x6 = x4 - x5 - 3.2233037561;
    double x7 = sqrt(x6);
    double x8 = 2.7424766666666671*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x12*x2*(47.948876324703562*x2 - 15.822734032286659);
    double x14 = 410439.95507083513*x13;
    double x15 = 1.0/(-x4 + x5 + 3.2233037561);
    double x16 = x3*((20.549518424872957*x2 - 9.4936404193719959)*(20.549518424872957*x2 - 9.4936404193719959));
    double x17 = x15*x16;
    double x18 = 1125621.9998494806*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 410439.95507083513*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 822.74300000000005*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 337686599.95484412*x17*x22;
    double x29 = exp(x8);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x8);
    double x33 = x12*x32;
    double x34 = -410439.95507083519*x31 + 149660.32712675835*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x12*x38;
    double x40 = -410439.95507083507*x37 + 498.86775708919447*x39;
    double x41 = x12*x16;

    result += x0*(-933893074.57053876*x0 + x11*x14 - x11*x18 + x11*x21 + x13*x34 - x13*x40 - x14*x27 - 299502824.00102746*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 + 2394625183.5213408*x3 - 136813.31835694506*x41*(8.2274300000000018*x12*x29/((x30)*(x30)) + x15*(-9.0000000000000018*x31 + 3.2817052226515453*x33) - 1.0939017408838483*x19*x32) + 136813.31835694503*x41*(2468.2290000000003*x12*x22*x35/((x36)*(x36)) + x15*(-9.0000000000000018*x37 + 0.010939017408838484*x39) - 0.0036463391362794944*x19*x38) + x28*exp(-1645.4860000000001*x23)/((x26)*(x26)) - x18*exp(-5.4849533333333342*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-28852876389.704491*x2 + 31226836534.221478*x3 + 3265399458.0572615);
    double x5 = 1.0/T;
    double x6 = 28.480921258115991*x2;
    double x7 = 30.824277637309436*x3;
    double x8 = x6 - x7 - 3.2233037561;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 822.74300000000005*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1645.4860000000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.010939017408838484*x19;
    double x25 = x5*(T*x24 - 9.0000000000000018/x21);
    double x26 = 136813.31835694503*x9;
    double x27 = x17*(-x6 + x7 + 3.2233037561);

    result += x0*(-833487858919.94495*x15*x18 - x15*x4 - 277829286306.64832*x16*x18 - x16*x4 + x26*(0.0036463391362794944*x19 - 2468.2290000000003*x23 + x25) - x26*(-2030718.1321470004*x22*x27 - 2468.2290000000012*x23 + x24 + 3.0000000000000004*x25 + 4061436.2642940008*x27*exp(x14)/((x21)*(x21)*(x21))) - 555658572613.29663*x18*exp(-2468.2290000000003*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(13878594015.209545*x2 - 6411750308.8232193);
    double x5 = 28.480921258115991*x2;
    double x6 = 30.824277637309436*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 3.2233037561;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 822.74300000000005*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1645.4860000000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 20.549518424872957*x2 - 9.4936404193719959;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0036463391362794944*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2468.2290000000003*x3;
    double x26 = 0.010939017408838484*T*x20;
    double x27 = x19*x26 - 9.0000000000000018/x23;
    double x28 = x0*x27;
    double x29 = 136813.31835694503*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 3.2233037561);

    result += x0*x1*x2*(833487858919.94495*x14*x18 - x14*x4 + 277829286306.64832*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7404.6870000000017*x0*x31 - x26*x30 + 3.0000000000000004*x27*x32) + 2030718.1321470004*x16*x24 - 4061436.2642940008*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 555658572613.29663*x18*exp(-2468.2290000000003*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(16191693017.744469*x2 - 5343125257.3526831);
    double x4 = 28.480921258115991*x2;
    double x5 = 30.824277637309436*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 3.2233037561;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 822.74300000000005*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1645.4860000000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 20.549518424872957*x2 - 9.4936404193719959;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0036463391362794944*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2468.2290000000003*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 3.2233037561;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0000000000000018*x32 + 0.010939017408838484*x33));
    double x35 = 47.948876324703562*x2 - 15.822734032286659;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0036463391362794944*x33;
    double x40 = 3.0000000000000004*x31;
    double x41 = 3.0000000000000004*x39/((x30)*(x30));

    result += -x38*(833487858919.94495*x13*x20 + x13*x3 + 277829286306.64832*x14*x20 + x14*x3 + 136813.31835694503*x15*x34 + 136813.31835694503*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(61.648555274618872*x2 - 28.480921258115988)/pow(x6, 5.0/2.0) + x24*x35 - 2030718.1321470004*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(41.099036849745914*x2 - 18.987280838743992) + 4061436.2642940008*x37*exp(x12)/((x26)*(x26)*(x26))) + 555658572613.29663*x20*exp(-2468.2290000000003*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 28.480921258115991*x2;
    double x5 = 30.824277637309436*x3;
    double x6 = x4 - x5 - 3.2233037561;
    double x7 = sqrt(x6);
    double x8 = 2.7424766666666671*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 20.549518424872957*x2 - 9.4936404193719959;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.4849533333333342*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 3.2233037561;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1125621.9998494806*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = 159.82958774901186*x2 - 42.193957419431086;
    double x24 = x2*x22*x23;
    double x25 = 410439.95507083513*x21;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 822.74300000000005*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 410439.95507083513*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1645.4860000000001*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 337686599.95484412*x26;
    double x38 = x19*x37;
    double x39 = 41.099036849745914*x2 - 18.987280838743992;
    double x40 = ((x12)*(x12));
    double x41 = x0*x40;
    double x42 = x18*x39*x41;
    double x43 = 1125621.9998494806*x42;
    double x44 = (61.648555274618872*x2 - 28.480921258115988)/pow(x6, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 47.948876324703562*x2 - 15.822734032286659;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = 1125621.9998494806*x49;
    double x51 = x16*x50;
    double x52 = x47*(95.897752649407124*x2 - 31.645468064573318);
    double x53 = x11*x47;
    double x54 = 1231319.8652125055*x46*x53;
    double x55 = x21*x50;
    double x56 = x37*x42;
    double x57 = x37*x49;
    double x58 = x36*x57;
    double x59 = x31*x57;
    double x60 = 136813.31835694506*x22;
    double x61 = exp(x8);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x8);
    double x65 = x22*x64;
    double x66 = -3*x63 + 1.0939017408838483*x65;
    double x67 = x2*x66;
    double x68 = exp(x28);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x28);
    double x72 = T*x22*x71;
    double x73 = -3*x70 + 0.0036463391362794944*x72;
    double x74 = x2*x73;
    double x75 = 136813.31835694503*x22;
    double x76 = x46*x66;
    double x77 = x46*x73;
    double x78 = 1.0939017408838483*x64;
    double x79 = x11*x78;
    double x80 = x61/((x62)*(x62));
    double x81 = 8.2274300000000018*x80;
    double x82 = x22*x81;
    double x83 = x49*(-9.0000000000000018*x63 + 3.2817052226515453*x65) - x79 + x82;
    double x84 = 273626.63671389013*x83;
    double x85 = x22*x48;
    double x86 = 0.0036463391362794944*T*x71;
    double x87 = x11*x86;
    double x88 = x68/((x69)*(x69));
    double x89 = 2468.2290000000003*x26*x88;
    double x90 = x22*x89;
    double x91 = x49*(-9.0000000000000018*x70 + 0.010939017408838484*x72) - x87 + x90;
    double x92 = 273626.63671389007*x91;
    double x93 = x12*x2*x44;
    double x94 = x2*x40;
    double x95 = x49*x94;
    double x96 = x11*x94;
    double x97 = 3.0000000000000004*x49;
    double x98 = 3.0000000000000004*x18;
    double x99 = x12*x39;
    double x100 = x94*x97;
    double x101 = x33*x95;
    double x102 = x74*x98;

    result += (3735572298.282155*x0 - 9260976.2102216128*x14*x16 - 3086992.0700738709*x14*x21 + x14*x84 - x14*x92 - x16*x20 - x16*x43 + 798674197.33607316*x2 - x20*x21 - x21*x43 - x21*x54 - x23*x60*x67 + x23*x74*x75 - x24*x25 + x24*x32 - x25*x45 - 7982083945.0711365*x3 + 277829286306.64832*x31*x34 + x31*x38 + x31*x54 + x31*x56 + x32*x45 + 833487858919.94495*x34*x36 + x36*x38 + x36*x56 - 136813.31835694506*x45*x66 + 136813.31835694503*x45*x73 - x47*x60*(x100*x83 + x46*x79 - x46*x82 + x66*x94*x98 + x67*x98*x99 - x76*x97 + x78*x93 - 22.563534801633342*x80*x95 - x81*x96 + 45.127069603266683*x95*exp(x15)/((x62)*(x62)*(x62))) + x47*x75*(x100*x91 - 2030718.1321470004*x101*x88 + 4061436.2642940008*x101*exp(x35)/((x69)*(x69)*(x69)) + x102*x40 + x102*x99 + x46*x87 - x46*x90 - x77*x97 + x86*x93 - x89*x96) + x48*x51 + x48*x55 - x48*x58 - x48*x59 + x51*x52 + x52*x55 - x52*x58 - x52*x59 - 410439.95507083519*x53*x76 + 410439.95507083507*x53*x77 + x84*x85 - x85*x92 + 555658572613.29663*x34*exp(-2468.2290000000003*x27)/((x30)*(x30)*(x30)) - 6173984.1401477419*x14*exp(-8.2274300000000018*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*12.512;
static const double Vmax = 1.15*12.512;
static double V = 0.9*12.512;

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



const char *Grossular_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Grossular_slb_em_coder_calib_name(void) {
    return "Grossular_slb_em";
}

const char *Grossular_slb_em_coder_calib_formula(void) {
    return "Ca3Al2Si3O12";
}

const double Grossular_slb_em_coder_calib_mw(void) {
    return 450.45238;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,12.0,0.0,0.0,0.0,
        0.0,2.0,3.0,0.0,0.0,0.0,
        0.0,0.0,3.0,0.0,0.0,0.0,
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

const double *Grossular_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Grossular_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Grossular_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Grossular_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Grossular_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Grossular_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Grossular_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Grossular_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Grossular_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Grossular_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Grossular_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Grossular_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Grossular_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Grossular_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Grossular_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Grossular_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Grossular_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Grossular_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Grossular_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Grossular_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Grossular_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Grossular_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Grossular_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Grossular_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Grossular_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Grossular_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Grossular_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Grossular_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Grossular_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Grossular_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Grossular_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Grossular_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Grossular_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Grossular_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Grossular_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Grossular_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Grossular_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

