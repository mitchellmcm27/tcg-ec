
static char *identifier = "HPClinoferrosilite_slb_em.emml:8367a7349510c97d4ff66df55011d7ab4c8777be:Thu Feb 10 16:52:38 2022";



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
    double x3 = sqrt(-12.405223829538892*x1 + 41.334448283499526*x2 + 1.1150177016250002);
    double x4 = 2.3052133333333336*x3;
    double x5 = 691.56399999999996*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) - 26.76322271834383*T + 38878963.430427521*x1 - 232552779.15056384*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) - 3209230.3395718634 + 380083869.91161942/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(41.334448283499526*pow(x0, 4.0/3.0) - 12.405223829538892*pow(x0, 2.0/3.0) + 1.1150177016250002);
    double x2 = x1/T;
    double x3 = 691.56399999999996*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -172499.49078181581*x2*x5/x6 + 57499.830260605275*x2*(-0.0043379933021383413*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-12.405223829538892*x1 + 41.334448283499526*x3 + 1.1150177016250002);
    double x6 = 2.3052133333333336*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(4.1350746098462974*x2 - 27.556298855666348*x4);
    double x10 = 691.56399999999996*x5/T;
    double x11 = exp(-x10);

    result += 172499.49078181581*x11*x9/(-x11 + 1) - 25919308.953618348*x2 + 310070372.20075178*x4 - 172499.49078181584*x7*x9/(-x7 + 1) + 57499.830260605282*x9*(-1.3013979906415023*x8*Debye(x6) + 3/(exp(x6) - 1)) - 57499.830260605275*x9*(-0.0043379933021383413*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 760167739.82323885/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-12.405223829538892*x2 + 41.334448283499526*x3 + 1.1150177016250002);
    double x5 = x0*x4;
    double x6 = 691.56399999999996*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1479874203.0618722*x2 + 4930969771.5321064*x3 + 133015409.90038808);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1383.1279999999999*x5)/((x8)*(x8)) + 57499.830260605275*x4*(x0*(0.013013979906415021*T*x11 - 8.9999999999999982/x13) + 0.0043379933021383413*x11 - 2074.692*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 41.334448283499526*pow(x1, 4.0/3.0) - 12.405223829538892*x2 + 1.1150177016250002;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 691.56399999999996*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 119294437.84303567*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(27.556298855666348*x2 - 4.1350746098462974)*(57499.830260605275*x4*(-2074.692*x0*x11*x12/((x13)*(x13)) + 0.0043379933021383413*x10/pow(x3, 3.0/2.0) + (0.013013979906415021*x10*x11 - 8.9999999999999982/x13)/x3) + x7*x9/x8 + x9*exp(-1383.1279999999999*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -12.405223829538892*x2 + 41.334448283499526*x3 + 1.1150177016250002;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x2*x6*(64.298030663221482*x2 - 6.8917910164104956);
    double x8 = 2.3052133333333336*x5;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 172499.49078181584*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((27.556298855666348*x2 - 4.1350746098462974)*(27.556298855666348*x2 - 4.1350746098462974));
    double x15 = x13*x14;
    double x16 = 397648.12614345236*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 1.0/T;
    double x20 = x19*x5;
    double x21 = 691.56399999999996*x20;
    double x22 = exp(-x21);
    double x23 = -x22 + 1;
    double x24 = x22/x23;
    double x25 = 172499.49078181581*x24;
    double x26 = 119294437.84303567*x15*x19;
    double x27 = exp(x8);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x8);
    double x31 = x30*x6;
    double x32 = -172499.49078181584*x29 + 74830.163563379159*x31;
    double x33 = exp(x21);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x21);
    double x37 = T*x36*x6;
    double x38 = -172499.49078181584*x35 + 249.43387854459721*x37;
    double x39 = x14*x6;

    result += x0*(2280503219.4697165*x0 + x11*x16 + x12*x18 - x12*x7 - x18*x25 + x18*x32 - x18*x38 + 43198848.256030574*x2 - x24*x26 + x25*x7 - 723497535.13508749*x3 - x32*x7 + x38*x7 + 57499.830260605282*x39*(x13*(-9.0*x29 + 3.9041939719245069*x31) + 1.3013979906415023*x17*x30 - 6.9156400000000007*x27*x6/((x28)*(x28))) - 57499.830260605275*x39*(0.0043379933021383413*T*x17*x36 + x13*(-8.9999999999999982*x35 + 0.013013979906415021*x37) - 2074.692*x19*x33*x6/((x34)*(x34))) - x26*exp(-1383.1279999999999*x20)/((x23)*(x23)) + x16*exp(-4.6104266666666671*x5)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-4439622609.1856165*x2 + 14792909314.596319*x3 + 399046229.70116419);
    double x5 = 1.0/T;
    double x6 = -12.405223829538892*x2 + 41.334448283499526*x3 + 1.1150177016250002;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 691.56399999999996*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1383.1279999999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.013013979906415021*x17;
    double x23 = x5*(T*x22 - 8.9999999999999982/x19);
    double x24 = 57499.830260605275*x7;
    double x25 = x15*x6;

    result += x0*(-247499215837.44333*x13*x16 + x13*x4 - 82499738612.48111*x14*x16 + x14*x4 + x24*(0.0043379933021383413*x17 - 2074.692*x21 + x23) - x24*(1434782.2982879998*x20*x25 - 2074.6919999999991*x21 + x22 + 2.9999999999999996*x23 - 2869564.5965759996*x25*exp(x12)/((x19)*(x19)*(x19))) - 164999477224.96222*x16*exp(-2074.692*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(6574626362.0428085*x2 - 986582802.0412482);
    double x5 = 41.334448283499526*pow(x1, 4.0/3.0) - 12.405223829538892*x2 + 1.1150177016250002;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 691.56399999999996*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1383.1279999999999*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 27.556298855666348*x2 - 4.1350746098462974;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0043379933021383413*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2074.692*x22*x3;
    double x24 = 0.013013979906415021*T*x18;
    double x25 = x17*x24 - 8.9999999999999982/x21;
    double x26 = x0*x25;
    double x27 = 57499.830260605275*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-247499215837.44333*x12*x16 + x12*x4 - 82499738612.48111*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-6224.0759999999991*x0*x17*x22 + x24*x28 + 2.9999999999999996*x25*x29) + 1434782.2982879998*x15*x22 - 2869564.5965759996*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 164999477224.96222*x16*exp(-2074.692*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(7670397422.3832769*x2 - 822152335.03437352);
    double x4 = 41.334448283499526*pow(x1, 4.0/3.0) - 12.405223829538892*x2 + 1.1150177016250002;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 691.56399999999996*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1383.1279999999999*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 27.556298855666348*x2 - 4.1350746098462974;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0043379933021383413*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2074.692*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-8.9999999999999982*x29 + 0.013013979906415021*x30));
    double x32 = 64.298030663221482*x2 - 6.8917910164104956;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0043379933021383413*x30;
    double x37 = 2.9999999999999996*x28;
    double x38 = 2.9999999999999996*x36/((x4)*(x4));

    result += x35*(-247499215837.44333*x11*x18 + x11*x3 - 82499738612.48111*x12*x18 + x12*x3 + 57499.830260605275*x13*x31 - 164999477224.96222*x18*exp(-2074.692*x6)/((x9)*(x9)*(x9)) - 57499.830260605275*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(82.668896566999052*x2 - 12.405223829538892)/pow(x4, 5.0/2.0) - x22*x32 + 1434782.2982879998*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(55.112597711332697*x2 - 8.2701492196925948) - 2869564.5965759996*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -12.405223829538892*x2 + 41.334448283499526*x3 + 1.1150177016250002;
    double x5 = sqrt(x4);
    double x6 = 2.3052133333333336*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 27.556298855666348*x2 - 4.1350746098462974;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 4.6104266666666671*x5;
    double x16 = exp(-x15)/((x8)*(x8));
    double x17 = 397648.12614345236*x16;
    double x18 = x7/x8;
    double x19 = 397648.12614345236*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(214.32676887740493*x2 - 18.378109377094653);
    double x22 = x20*x21;
    double x23 = 172499.49078181584*x18;
    double x24 = 1.0/T;
    double x25 = x24*x5;
    double x26 = 691.56399999999996*x25;
    double x27 = exp(-x26);
    double x28 = -x27 + 1;
    double x29 = x27/x28;
    double x30 = 172499.49078181581*x29;
    double x31 = pow(T, -2);
    double x32 = x12*x31;
    double x33 = 1383.1279999999999*x25;
    double x34 = exp(-x33)/((x28)*(x28));
    double x35 = 119294437.84303567*x24;
    double x36 = x14*x35;
    double x37 = x13*(55.112597711332697*x2 - 8.2701492196925948);
    double x38 = ((x10)*(x10));
    double x39 = x0*x38;
    double x40 = x37*x39;
    double x41 = (82.668896566999052*x2 - 12.405223829538892)/pow(x4, 5.0/2.0);
    double x42 = x39*x41;
    double x43 = 1.0/x4;
    double x44 = 64.298030663221482*x2 - 6.8917910164104956;
    double x45 = x43*x44;
    double x46 = x10*x3;
    double x47 = x17*x46;
    double x48 = x43*(128.59606132644296*x2 - 13.783582032820991);
    double x49 = x44*x46;
    double x50 = x49*x9;
    double x51 = x19*x46;
    double x52 = x35*x40;
    double x53 = x35*x46;
    double x54 = x34*x53;
    double x55 = x29*x53;
    double x56 = exp(x6);
    double x57 = x56 - 1;
    double x58 = 1.0/x57;
    double x59 = Debye(x6);
    double x60 = x20*x59;
    double x61 = -3*x58 + 1.3013979906415023*x60;
    double x62 = 57499.830260605282*x20;
    double x63 = exp(x26);
    double x64 = x63 - 1;
    double x65 = 1.0/x64;
    double x66 = Debye(x26);
    double x67 = T*x20*x66;
    double x68 = -3*x65 + 0.0043379933021383413*x67;
    double x69 = 57499.830260605275*x20;
    double x70 = 172499.49078181584*x50;
    double x71 = 1.3013979906415023*x59;
    double x72 = x71*x9;
    double x73 = x56/((x57)*(x57));
    double x74 = 6.9156400000000007*x73;
    double x75 = x20*x74;
    double x76 = x43*(-9.0*x58 + 3.9041939719245069*x60) + x72 - x75;
    double x77 = 114999.66052121056*x76;
    double x78 = x20*x49;
    double x79 = 0.0043379933021383413*T*x66;
    double x80 = x79*x9;
    double x81 = x63/((x64)*(x64));
    double x82 = 2074.692*x24*x81;
    double x83 = x20*x82;
    double x84 = x43*(-8.9999999999999982*x65 + 0.013013979906415021*x67) + x80 - x83;
    double x85 = 114999.66052121055*x84;
    double x86 = x10*x2;
    double x87 = x41*x86;
    double x88 = x2*x38;
    double x89 = x43*x88;
    double x90 = x88*x9;
    double x91 = 3.0*x61;
    double x92 = x13*x88;
    double x93 = x37*x86;
    double x94 = x31*x89;
    double x95 = 2.9999999999999996*x68;

    result += (-9122012877.8788662*x0 + 2749991.2870827052*x12*x16 + 916663.76236090169*x12*x18 + x12*x77 - x12*x85 + 1833327.5247218034*x12*exp(-6.9156400000000007*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x19 + x17*x40 - 517498.47234544752*x18*x50 + x19*x40 - 115196928.6827482*x2 + x21*x61*x62 - x21*x68*x69 + x22*x23 - x22*x30 + x23*x42 - 82499738612.48111*x29*x32 - x29*x36 + 517498.4723454474*x29*x50 - x29*x52 + 2411658450.4502916*x3 - x30*x42 - 247499215837.44333*x32*x34 - x34*x36 - x34*x52 + 57499.830260605282*x42*x61 - 57499.830260605275*x42*x68 - x45*x47 - x45*x51 + x45*x54 + x45*x55 + x46*x62*(-x44*x72 + x44*x75 - x45*x91 + x71*x87 + 15.942025536533336*x73*x89 - x74*x90 + 3.0*x76*x89 + x91*x92 + x91*x93 - 31.884051073066672*x89*exp(x15)/((x57)*(x57)*(x57))) - x46*x69*(-x44*x80 + x44*x83 - x45*x95 + x79*x87 + 1434782.2982879998*x81*x94 - x82*x90 + 2.9999999999999996*x84*x89 + x92*x95 + x93*x95 - 2869564.5965759996*x94*exp(x33)/((x64)*(x64)*(x64))) - x47*x48 - x48*x51 + x48*x54 + x48*x55 - x61*x70 + x68*x70 - x77*x78 + x78*x85 - 164999477224.96222*x32*exp(-2074.692*x25)/((x28)*(x28)*(x28)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.385413;
static const double Vmax = 1.15*6.385413;
static double V = 0.9*6.385413;

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



const char *HPClinoferrosilite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *HPClinoferrosilite_slb_em_coder_calib_name(void) {
    return "HPClinoferrosilite_slb_em";
}

const char *HPClinoferrosilite_slb_em_coder_calib_formula(void) {
    return "Fe2Si2O6";
}

const double HPClinoferrosilite_slb_em_coder_calib_mw(void) {
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

const double *HPClinoferrosilite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double HPClinoferrosilite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double HPClinoferrosilite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int HPClinoferrosilite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **HPClinoferrosilite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **HPClinoferrosilite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void HPClinoferrosilite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int HPClinoferrosilite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double HPClinoferrosilite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int HPClinoferrosilite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double HPClinoferrosilite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double HPClinoferrosilite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

