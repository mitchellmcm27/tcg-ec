
static char *identifier = "MgCaFerrite_slb_em.emml:b1761122b6211beee5251ee7be74df5fd6c2fa19:Tue Aug 22 19:33:51 2023";



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
    double x3 = sqrt(-13.283234193829802*x1 + 28.670877818546344*x2 + 1.4703920127999988);
    double x4 = 2.7704740333333335*x3;
    double x5 = 831.14220999999998*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 37721065.590121299*x1 + 40811559.401290312*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 6103033.2162500024 + 5653095.3067274103/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(28.670877818546344*pow(x0, 4.0/3.0) - 13.283234193829802*pow(x0, 2.0/3.0) + 1.4703920127999988);
    double x2 = x1/T;
    double x3 = 831.14220999999998*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -145120.51754369965*x2*x5/x6 + 48373.505847899884*x2*(-0.0036094906069082931*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-13.283234193829802*x1 + 28.670877818546344*x3 + 1.4703920127999988);
    double x6 = 2.7704740333333335*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(4.4277447312766007*x2 - 19.113918545697562*x4);
    double x10 = 831.14220999999998*x5/T;
    double x11 = exp(-x10);

    result += 145120.51754369965*x11*x9/(-x11 + 1) + 25147377.060080864*x2 - 54415412.535053745*x4 - 145120.51754369968*x7*x9/(-x7 + 1) + 48373.505847899891*x9*(-1.0828471820724879*x8*Debye(x6) + 3/(exp(x6) - 1)) - 48373.505847899884*x9*(-0.0036094906069082931*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 11306190.613454821/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-13.283234193829802*x2 + 28.670877818546344*x3 + 1.4703920127999988);
    double x5 = x0*x4;
    double x6 = 831.14220999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1602167755.0621691*x2 + 3458160511.2058983*x3 + 177352490.80404064);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1662.28442*x5)/((x8)*(x8)) + 48373.505847899884*x4*(x0*(0.01082847182072488*T*x11 - 9.0/x13) + 0.0036094906069082931*x11 - 2493.4266299999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 28.670877818546344*pow(x1, 4.0/3.0) - 13.283234193829802*x2 + 1.4703920127999988;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 831.14220999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 120615787.6676143*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(19.113918545697562*x2 - 4.4277447312766007)*(48373.505847899884*x4*(-2493.4266299999999*x0*x11*x12/((x13)*(x13)) + 0.0036094906069082931*x10/pow(x3, 3.0/2.0) + (0.01082847182072488*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1662.28442*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -13.283234193829802*x2 + 28.670877818546344*x3 + 1.4703920127999988;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x2*x6*(44.59914327329431*x2 - 7.3795745521276679);
    double x8 = 2.7704740333333335*x5;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 145120.51754369968*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((19.113918545697562*x2 - 4.4277447312766007)*(19.113918545697562*x2 - 4.4277447312766007));
    double x15 = x13*x14;
    double x16 = 402052.62555871444*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 1.0/T;
    double x20 = x19*x5;
    double x21 = 831.14220999999998*x20;
    double x22 = exp(-x21);
    double x23 = -x22 + 1;
    double x24 = x22/x23;
    double x25 = 145120.51754369965*x24;
    double x26 = 120615787.6676143*x15*x19;
    double x27 = exp(x8);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x8);
    double x31 = x30*x6;
    double x32 = -145120.51754369968*x29 + 52381.114494365407*x31;
    double x33 = exp(x21);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x21);
    double x37 = T*x36*x6;
    double x38 = -145120.51754369965*x35 + 174.60371498121802*x37;
    double x39 = x14*x6;

    result += x0*(33918571.840364464*x0 + x11*x16 + x12*x18 - x12*x7 - x18*x25 + x18*x32 - x18*x38 - 41912295.100134775*x2 - x24*x26 + x25*x7 + 126969295.9151254*x3 - x32*x7 + x38*x7 + 48373.505847899891*x39*(x13*(-9.0*x29 + 3.2485415462174636*x31) + 1.0828471820724879*x17*x30 - 8.3114221000000015*x27*x6/((x28)*(x28))) - 48373.505847899884*x39*(0.0036094906069082931*T*x17*x36 + x13*(-9.0*x35 + 0.01082847182072488*x37) - 2493.4266299999999*x19*x33*x6/((x34)*(x34))) - x26*exp(-1662.28442*x20)/((x23)*(x23)) + x16*exp(-5.5409480666666671*x5)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-4806503265.1865082*x2 + 10374481533.617695*x3 + 532057472.41212195);
    double x5 = 1.0/T;
    double x6 = -13.283234193829802*x2 + 28.670877818546344*x3 + 1.4703920127999988;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 831.14220999999998*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1662.28442*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = -2493.4266299999999*x0*x20*x7;
    double x22 = 0.01082847182072488*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 48373.505847899884*x7;
    double x25 = x15*x6;

    result += x0*(-300746616968.8551*x13*x16 + x13*x4 - 100248872322.95169*x14*x16 + x14*x4 + x24*(0.0036094906069082931*x17 + x21 + x23) - x24*(2072392.1197310521*x20*x25 + x21 + x22 + 3.0*x23 - 4144784.2394621042*x25*exp(x12)/((x19)*(x19)*(x19))) - 200497744645.90338*x16*exp(-2493.4266299999999*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(4610880681.6078644*x2 - 1068111836.7081128);
    double x5 = 28.670877818546344*pow(x1, 4.0/3.0) - 13.283234193829802*x2 + 1.4703920127999988;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 831.14220999999998*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1662.28442*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 19.113918545697562*x2 - 4.4277447312766007;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0036094906069082931*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2493.4266299999999*x22*x3;
    double x24 = 0.01082847182072488*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 48373.505847899884*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-300746616968.8551*x12*x16 + x12*x4 - 100248872322.95169*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7480.2798899999998*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 2072392.1197310521*x15*x22 - 4144784.2394621042*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 200497744645.90338*x16*exp(-2493.4266299999999*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(5379360795.2091751*x2 - 890093197.2567606);
    double x4 = 28.670877818546344*pow(x1, 4.0/3.0) - 13.283234193829802*x2 + 1.4703920127999988;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 831.14220999999998*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1662.28442*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 19.113918545697562*x2 - 4.4277447312766007;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0036094906069082931*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2493.4266299999999*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.01082847182072488*x30));
    double x32 = 44.59914327329431*x2 - 7.3795745521276679;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0036094906069082931*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-300746616968.8551*x11*x18 + x11*x3 - 100248872322.95169*x12*x18 + x12*x3 + 48373.505847899884*x13*x31 - 200497744645.90338*x18*exp(-2493.4266299999999*x6)/((x9)*(x9)*(x9)) - 48373.505847899884*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(57.341755637092689*x2 - 13.283234193829802)/pow(x4, 5.0/2.0) - x22*x32 + 2072392.1197310521*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(38.227837091395124*x2 - 8.8554894625532015) - 4144784.2394621042*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -13.283234193829802*x2 + 28.670877818546344*x3 + 1.4703920127999988;
    double x5 = sqrt(x4);
    double x6 = 2.7704740333333335*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 19.113918545697562*x2 - 4.4277447312766007;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = 5.5409480666666671*x5;
    double x14 = exp(-x13)/((x8)*(x8));
    double x15 = pow(x4, -2);
    double x16 = 402052.62555871444*x15;
    double x17 = x11*x16;
    double x18 = x7/x8;
    double x19 = 145120.51754369968*x18;
    double x20 = 1.0/x5;
    double x21 = 148.66381091098103*x2 - 19.678865472340448;
    double x22 = x2*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 831.14220999999998*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 145120.51754369965*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1662.28442*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = x11*x15;
    double x35 = 120615787.6676143*x23;
    double x36 = x33*x35;
    double x37 = x28*x35;
    double x38 = 38.227837091395124*x2 - 8.8554894625532015;
    double x39 = ((x10)*(x10));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = x16*x41;
    double x43 = (57.341755637092689*x2 - 13.283234193829802)/pow(x4, 5.0/2.0);
    double x44 = x40*x43;
    double x45 = 44.59914327329431*x2 - 7.3795745521276679;
    double x46 = x10*x3;
    double x47 = x45*x46;
    double x48 = 1.0/x4;
    double x49 = 402052.62555871444*x48;
    double x50 = x14*x49;
    double x51 = x46*(89.198286546588619*x2 - 14.759149104255336);
    double x52 = x47*x9;
    double x53 = x18*x49;
    double x54 = x15*x41;
    double x55 = x36*x48;
    double x56 = x37*x48;
    double x57 = exp(x6);
    double x58 = x57 - 1;
    double x59 = 1.0/x58;
    double x60 = Debye(x6);
    double x61 = x20*x60;
    double x62 = -3*x59 + 1.0828471820724879*x61;
    double x63 = x2*x62;
    double x64 = 48373.505847899891*x20;
    double x65 = 48373.505847899884*x20;
    double x66 = exp(x25);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x25);
    double x70 = T*x20*x69;
    double x71 = -3*x68 + 0.0036094906069082931*x70;
    double x72 = x2*x71;
    double x73 = x45*x62;
    double x74 = x46*x9;
    double x75 = x45*x71;
    double x76 = 1.0828471820724879*x60;
    double x77 = x76*x9;
    double x78 = x57/((x58)*(x58));
    double x79 = 8.3114221000000015*x78;
    double x80 = x20*x79;
    double x81 = x48*(-9.0*x59 + 3.2485415462174636*x61) + x77 - x80;
    double x82 = 96747.011695799782*x81;
    double x83 = x20*x47;
    double x84 = 0.0036094906069082931*T*x69;
    double x85 = x84*x9;
    double x86 = x66/((x67)*(x67));
    double x87 = 2493.4266299999999*x23*x86;
    double x88 = x20*x87;
    double x89 = x48*(-9.0*x68 + 0.01082847182072488*x70) + x85 - x88;
    double x90 = 96747.011695799767*x89;
    double x91 = x10*x2*x43;
    double x92 = x2*x39;
    double x93 = x48*x92;
    double x94 = x9*x92;
    double x95 = 3.0*x48;
    double x96 = 3.0*x15;
    double x97 = x63*x96;
    double x98 = x10*x38;
    double x99 = x92*x95;
    double x100 = x30*x93;

    result += (-135674287.36145785*x0 + 3341629.0774317244*x12*x14 + 1113876.3591439081*x12*x18 + x12*x82 - x12*x90 + 2227752.7182878163*x12*exp(-8.3114221000000015*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x42 + x17*x18 + x18*x42 - 435361.55263109901*x18*x52 + x19*x22 + x19*x44 + 111766120.26702607*x2 + x21*x63*x64 - x21*x65*x72 - x22*x29 - 100248872322.95169*x28*x31 + 435361.55263109895*x28*x52 - x29*x44 - 423230986.38375133*x3 - 300746616968.8551*x31*x33 - x34*x36 - x34*x37 - x36*x54 - x37*x54 + 48373.505847899891*x44*x62 - 48373.505847899884*x44*x71 + x46*x64*(x39*x97 - x45*x77 + x45*x80 - x73*x95 + x76*x91 + 23.02657910812281*x78*x93 - x79*x94 + x81*x99 + x97*x98 - 46.05315821624562*x93*exp(x13)/((x58)*(x58)*(x58))) - x46*x65*(2072392.1197310521*x100*x86 - 4144784.2394621042*x100*exp(x32)/((x67)*(x67)*(x67)) - x45*x85 + x45*x88 + x71*x92*x96 + x72*x96*x98 - x75*x95 + x84*x91 - x87*x94 + x89*x99) - x47*x50 - x47*x53 + x47*x55 + x47*x56 - x50*x51 - x51*x53 + x51*x55 + x51*x56 - 145120.51754369968*x73*x74 + 145120.51754369965*x74*x75 - x82*x83 + x83*x90 - 200497744645.90338*x31*exp(-2493.4266299999999*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*3.6135;
static const double Vmax = 1.15*3.6135;
static double V = 0.9*3.6135;

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

