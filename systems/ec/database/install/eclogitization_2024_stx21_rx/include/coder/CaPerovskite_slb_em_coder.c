
static char *identifier = "CaPerovskite_slb_em.emml:2beec7d215a9c1f08440bacb3919138b05bbe470:Tue Aug 22 19:32:48 2023";



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
    double x3 = sqrt(-25.818448177904244*x1 + 36.201168247903183*x2 + 4.7505941272000003);
    double x4 = 2.6626193666666667*x3;
    double x5 = 798.78580999999997*x3/T;

    result += 124.7169392722986*T*log(1 - exp(-x5)) - 41.572313090766201*T*Debye(x5) - 30718987.695086863*x1 + 32212753.277169228*x2 - 37415.08178168958*log(1 - exp(-x4)) + 12471.69392722986*Debye(x4) + 6192463.5700000003 - 2745753.6912187519/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(36.201168247903183*pow(x0, 4.0/3.0) - 25.818448177904244*pow(x0, 2.0/3.0) + 4.7505941272000003);
    double x2 = x1/T;
    double x3 = 798.78580999999997*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -99622.121357343844*x2*x5/x6 + 33207.373785781281*x2*(-0.0037557001669821853*T*x4/x1 + 3/(exp(x3) - 1)) - 41.572313090766201*x4 + 124.7169392722986*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-25.818448177904244*x1 + 36.201168247903183*x3 + 4.7505941272000003);
    double x6 = 2.6626193666666667*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(8.6061493926347481*x2 - 24.134112165268789*x4);
    double x10 = 99622.121357343844*x9;
    double x11 = 798.78580999999997*x5/T;
    double x12 = exp(-x11);
    double x13 = 33207.373785781281*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.1267100500946554*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0037557001669821853*T*x8*Debye(x11) + 3/(exp(x11) - 1)) + 20479325.130057909*x2 - 42950337.702892303*x4 + 5491507.3824375039/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-25.818448177904244*x2 + 36.201168247903183*x3 + 4.7505941272000003);
    double x5 = x0*x4;
    double x6 = 798.78580999999997*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2054547857.879894*x2 + 2880770841.2208881*x3 + 378036778.99001586);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1597.5716199999999*x5)/((x8)*(x8)) + 33207.373785781281*x4*(x0*(0.011267100500946555*T*x11 - 9.0/x13) + 0.0037557001669821853*x11 - 2396.35743*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 36.201168247903183*pow(x1, 4.0/3.0) - 25.818448177904244*x2 + 4.7505941272000003;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 798.78580999999997*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 79576736.902344197*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(24.134112165268789*x2 - 8.6061493926347481)*(33207.373785781281*x4*(-2396.35743*x0*x11*x12/((x13)*(x13)) + 0.0037557001669821853*x10/pow(x3, 3.0/2.0) + (0.011267100500946555*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1597.5716199999999*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -25.818448177904244*x2 + 36.201168247903183*x3 + 4.7505941272000003;
    double x5 = sqrt(x4);
    double x6 = 2.6626193666666667*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = x7/x8;
    double x10 = 1.0/x5;
    double x11 = x2*(56.312928385627174*x2 - 14.343582321057912);
    double x12 = 99622.121357343844*x10*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((24.134112165268789*x2 - 8.6061493926347481)*(24.134112165268789*x2 - 8.6061493926347481));
    double x15 = x13*x14;
    double x16 = 265255.78967448068*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 99622.121357343844*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 798.78580999999997*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 79576736.902344197*x15*x20;
    double x27 = exp(x6);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x6);
    double x31 = x10*x30;
    double x32 = -3*x29 + 1.1267100500946554*x31;
    double x33 = 33207.373785781281*x10;
    double x34 = x11*x33;
    double x35 = 33207.373785781281*x18;
    double x36 = exp(x22);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = Debye(x22);
    double x40 = T*x10*x39;
    double x41 = -3*x38 + 0.0037557001669821853*x40;
    double x42 = x14*x33;

    result += x0*(-16474522.147312511*x0 + x12*x25 - x12*x9 + x16*x9 + x16*exp(-5.3252387333333333*x5)/((x8)*(x8)) - x19*x25 + x19*x9 - 34132208.550096512*x2 - x25*x26 + 100217454.64008203*x3 - x32*x34 + x32*x35 + x34*x41 - x35*x41 - x42*(0.0037557001669821853*T*x17*x39 - 2396.35743*x10*x20*x36/((x37)*(x37)) + x13*(-9.0*x38 + 0.011267100500946555*x40)) + x42*(-7.9878581000000004*x10*x27/((x28)*(x28)) + x13*(-8.9999999999999982*x29 + 3.3801301502839656*x31) + 1.1267100500946554*x17*x30) - x26*exp(-1597.5716199999999*x21)/((x24)*(x24)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-6163643573.6396818*x2 + 8642312523.6626644*x3 + 1134110336.9700477);
    double x5 = 1.0/T;
    double x6 = -25.818448177904244*x2 + 36.201168247903183*x3 + 4.7505941272000003;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 798.78580999999997*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1597.5716199999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = -2396.35743*x0*x20*x7;
    double x22 = 0.011267100500946555*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 33207.373785781281*x7;
    double x25 = x15*x6;

    result += x0*(-190694304731.08771*x13*x16 + x13*x4 - 63564768243.6959*x14*x16 + x14*x4 + x24*(0.0037557001669821853*x17 + x21 + x23) - x24*(1914176.3107720683*x20*x25 + x21 + x22 + 3.0*x23 - 3828352.6215441367*x25*exp(x12)/((x19)*(x19)*(x19))) - 127129536487.3918*x16*exp(-2396.35743*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(3841027788.2945175*x2 - 1369698571.9199293);
    double x5 = 36.201168247903183*pow(x1, 4.0/3.0) - 25.818448177904244*x2 + 4.7505941272000003;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 798.78580999999997*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1597.5716199999999*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 24.134112165268789*x2 - 8.6061493926347481;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0037557001669821853*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2396.35743*x22*x3;
    double x24 = 0.011267100500946555*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 33207.373785781281*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-190694304731.08771*x12*x16 + x12*x4 - 63564768243.6959*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7189.0722900000001*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 1914176.3107720683*x15*x22 - 3828352.6215441367*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 127129536487.3918*x16*exp(-2396.35743*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(4481199086.3436041*x2 - 1141415476.599941);
    double x4 = 36.201168247903183*pow(x1, 4.0/3.0) - 25.818448177904244*x2 + 4.7505941272000003;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 798.78580999999997*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1597.5716199999999*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 24.134112165268789*x2 - 8.6061493926347481;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0037557001669821853*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2396.35743*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.011267100500946555*x30));
    double x32 = 56.312928385627174*x2 - 14.343582321057912;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0037557001669821853*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-190694304731.08771*x11*x18 + x11*x3 - 63564768243.6959*x12*x18 + x12*x3 + 33207.373785781281*x13*x31 - 127129536487.3918*x18*exp(-2396.35743*x6)/((x9)*(x9)*(x9)) - 33207.373785781281*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(72.402336495806367*x2 - 25.818448177904244)/pow(x4, 5.0/2.0) - x22*x32 + 1914176.3107720683*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(48.268224330537578*x2 - 17.212298785269496) - 3828352.6215441367*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -25.818448177904244*x2 + 36.201168247903183*x3 + 4.7505941272000003;
    double x5 = sqrt(x4);
    double x6 = 2.6626193666666667*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 24.134112165268789*x2 - 8.6061493926347481;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 265255.78967448068*x14;
    double x16 = 5.3252387333333333*x5;
    double x17 = exp(-x16)/((x8)*(x8));
    double x18 = x7/x8;
    double x19 = 99622.121357343844*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(187.70976128542389*x2 - 38.24955285615443);
    double x22 = x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 798.78580999999997*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 99622.121357343844*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1597.5716199999999*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = 79576736.902344197*x23;
    double x35 = x14*x34;
    double x36 = x13*(48.268224330537578*x2 - 17.212298785269496);
    double x37 = ((x10)*(x10));
    double x38 = x0*x37;
    double x39 = x36*x38;
    double x40 = 265255.78967448068*x17;
    double x41 = 265255.78967448068*x18;
    double x42 = (72.402336495806367*x2 - 25.818448177904244)/pow(x4, 5.0/2.0);
    double x43 = x38*x42;
    double x44 = 1.0/x4;
    double x45 = 56.312928385627174*x2 - 14.343582321057912;
    double x46 = x44*x45;
    double x47 = x10*x3;
    double x48 = x40*x47;
    double x49 = x44*(112.62585677125435*x2 - 28.687164642115825);
    double x50 = x45*x47;
    double x51 = x50*x9;
    double x52 = 298866.36407203152*x51;
    double x53 = x41*x47;
    double x54 = x34*x39;
    double x55 = x34*x47;
    double x56 = x33*x55;
    double x57 = x28*x55;
    double x58 = exp(x6);
    double x59 = x58 - 1;
    double x60 = 1.0/x59;
    double x61 = Debye(x6);
    double x62 = x20*x61;
    double x63 = -3*x60 + 1.1267100500946554*x62;
    double x64 = 33207.373785781281*x20;
    double x65 = x21*x64;
    double x66 = exp(x25);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x25);
    double x70 = T*x20*x69;
    double x71 = -3*x68 + 0.0037557001669821853*x70;
    double x72 = 33207.373785781281*x43;
    double x73 = 99622.121357343844*x51;
    double x74 = 1.1267100500946554*x61;
    double x75 = x74*x9;
    double x76 = x58/((x59)*(x59));
    double x77 = 7.9878581000000004*x76;
    double x78 = x20*x77;
    double x79 = x44*(-8.9999999999999982*x60 + 3.3801301502839656*x62) + x75 - x78;
    double x80 = 66414.747571562562*x79;
    double x81 = x20*x50;
    double x82 = 0.0037557001669821853*T*x69;
    double x83 = x82*x9;
    double x84 = x66/((x67)*(x67));
    double x85 = 2396.35743*x23*x84;
    double x86 = x20*x85;
    double x87 = x44*(-9.0*x68 + 0.011267100500946555*x70) + x83 - x86;
    double x88 = 66414.747571562562*x87;
    double x89 = x10*x2;
    double x90 = x42*x89;
    double x91 = x2*x37;
    double x92 = x44*x91;
    double x93 = x9*x91;
    double x94 = 2.9999999999999996*x63;
    double x95 = x13*x91;
    double x96 = x36*x89;
    double x97 = x47*x64;
    double x98 = x30*x92;
    double x99 = 3.0*x71;

    result += (65898088.589250043*x0 + 2118825.6081231972*x12*x17 + 706275.20270773233*x12*x18 + x12*x80 - x12*x88 + 1412550.4054154647*x12*exp(-7.9878581000000004*x5)/((x8)*(x8)*(x8)) + x15*x17 + x15*x18 - x18*x52 + x19*x22 + x19*x43 + 91019222.800257355*x2 - x22*x29 - 63564768243.6959*x28*x31 - x28*x35 + x28*x52 - x28*x54 - x29*x43 - 334058182.13360679*x3 - 190694304731.08771*x31*x33 - x33*x35 - x33*x54 + x39*x40 + x39*x41 - x46*x48 - x46*x53 + x46*x56 + x46*x57 - x48*x49 - x49*x53 + x49*x56 + x49*x57 + x63*x65 + x63*x72 - x63*x73 - x65*x71 - x71*x72 + x71*x73 - x80*x81 + x81*x88 + x97*(-x45*x75 + x45*x78 - x46*x94 + x74*x90 + 21.268625675245204*x76*x92 - x77*x93 + 2.9999999999999996*x79*x92 + x94*x95 + x94*x96 - 42.537251350490408*x92*exp(x16)/((x59)*(x59)*(x59))) - x97*(-x45*x83 + x45*x86 - x46*x99 + x82*x90 + 1914176.3107720683*x84*x98 - x85*x93 + 3.0*x87*x92 + x95*x99 + x96*x99 - 3828352.6215441367*x98*exp(x32)/((x67)*(x67)*(x67))) - 127129536487.3918*x31*exp(-2396.35743*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.745;
static const double Vmax = 1.15*2.745;
static double V = 0.9*2.745;

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



const char *CaPerovskite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *CaPerovskite_slb_em_coder_calib_name(void) {
    return "CaPerovskite_slb_em";
}

const char *CaPerovskite_slb_em_coder_calib_formula(void) {
    return "CaSiO3";
}

const double CaPerovskite_slb_em_coder_calib_mw(void) {
    return 116.16369999999999;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,3.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
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

const double *CaPerovskite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double CaPerovskite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double CaPerovskite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double CaPerovskite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double CaPerovskite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double CaPerovskite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double CaPerovskite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double CaPerovskite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double CaPerovskite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double CaPerovskite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double CaPerovskite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double CaPerovskite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double CaPerovskite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double CaPerovskite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double CaPerovskite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double CaPerovskite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double CaPerovskite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double CaPerovskite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double CaPerovskite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double CaPerovskite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int CaPerovskite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **CaPerovskite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **CaPerovskite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void CaPerovskite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int CaPerovskite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double CaPerovskite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int CaPerovskite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double CaPerovskite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double CaPerovskite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

