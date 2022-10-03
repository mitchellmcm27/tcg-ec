
static char *identifier = "Anorthite_slb_em.emml:58f035c76f748b0eb85e1a6930205619df53aeb6:Thu Feb 10 16:51:37 2022";



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
    double x3 = sqrt(12.743683058578839*x1 - 16.911048197110393*x2 - 0.95583226354999995);
    double x4 = 2.5079703333333336*x3;
    double x5 = 752.39110000000005*x3/T;

    result += 324.26404210797637*T*log(1 - exp(-x5)) - 108.08801403599212*T*Debye(x5) - 88713775.470475569*x1 + 206722857.38080168*x2 - 97279.212632392911*log(1 - exp(-x4)) + 32426.404210797635*Debye(x4) + 5503116.5541875008;
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-16.911048197110393*pow(x0, 4.0/3.0) + 12.743683058578839*pow(x0, 2.0/3.0) - 0.95583226354999995);
    double x2 = x1/T;
    double x3 = 752.39110000000005*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -243973.37933206666*x2*x5/x6 + 81324.459777355558*x2*(-0.0039872879942359765*T*x4/x1 + 3/(exp(x3) - 1)) - 108.08801403599212*x4 + 324.26404210797637*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(12.743683058578839*x1 - 16.911048197110393*x3 - 0.95583226354999995);
    double x6 = 2.5079703333333336*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-4.2478943528596123*x2 + 11.274032131406928*x4);
    double x10 = 752.39110000000005*x5/T;
    double x11 = exp(-x10);
    double x12 = 81324.459777355558*x9;

    result += 243973.37933206666*x11*x9/(-x11 + 1) + x12*(-1.1961863982707928*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.0039872879942359765*T*x8*Debye(x10) + 3/(exp(x10) - 1)) + 59142516.980317041*x2 - 275630476.50773555*x4 - 243973.37933206669*x7*x9/(-x7 + 1);
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(12.743683058578839*x2 - 16.911048197110393*x3 - 0.95583226354999995);
    double x5 = x0*x4;
    double x6 = 752.39110000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2339273781.1511207*x2 + 3104249491.880796*x3 + 175455819.40659106);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1504.7822000000001*x5)/((x8)*(x8)) - 81324.459777355558*x4*(x0*(0.011961863982707929*T*x11 - 9.0/x13) + 0.0039872879942359765*x11 - 2257.1733000000004*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 12.743683058578839*x2;
    double x4 = 16.911048197110393*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 0.95583226354999995;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 752.39110000000005*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 183563399.24637091*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(11.274032131406928*x2 - 4.2478943528596123)*(-81324.459777355558*x6*(2257.1733000000004*x0*x13*x14/((x15)*(x15)) - 0.0039872879942359765*x12/pow(x5, 3.0/2.0) + (0.011961863982707929*x12*x13 - 9.0/x15)/(-x3 + x4 + 0.95583226354999995)) + x11*x9/x10 + x11*exp(-1504.7822000000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = 12.743683058578839*x1;
    double x3 = 16.911048197110393*pow(x0, 4.0/3.0);
    double x4 = x2 - x3 - 0.95583226354999995;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = 6417982.0081945183*x1;
    double x8 = 2.5079703333333336*x5;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/T;
    double x13 = x12*x5;
    double x14 = 752.39110000000005*x13;
    double x15 = exp(-x14);
    double x16 = -x15 + 1;
    double x17 = x15/x16;
    double x18 = 1.0/(-x2 + x3 + 0.95583226354999995);
    double x19 = x1*((11.274032131406928*x1 - 4.2478943528596123)*(11.274032131406928*x1 - 4.2478943528596123));
    double x20 = x18*x19;
    double x21 = 611877.99748790311*x20;
    double x22 = pow(x4, -3.0/2.0);
    double x23 = x19*x22;
    double x24 = 183563399.24637091*x12*x20;
    double x25 = exp(x8);
    double x26 = x25 - 1;
    double x27 = 1.0/x26;
    double x28 = Debye(x8);
    double x29 = x28*x6;
    double x30 = -243973.37933206669*x27 + 97279.212632392897*x29;
    double x31 = x6*(26.306074973282833*x1 - 7.0798239214326868);
    double x32 = exp(x14);
    double x33 = x32 - 1;
    double x34 = 1.0/x33;
    double x35 = Debye(x14);
    double x36 = T*x35*x6;
    double x37 = -243973.37933206669*x34 + 324.26404210797637*x36;
    double x38 = 81324.459777355558*x19*x6;

    result += x1*(643137778.5180496*x1 - x11*x21 + 243973.37933206669*x11*x23 + x11*x6*(x7 - 1727288.5671879367) - 243973.37933206666*x17*x23 + x17*x24 - x17*x6*(x7 - 1727288.5671879365) + x23*x30 - x23*x37 + x30*x31 - x31*x37 - x38*(x18*(-9.0*x27 + 3.5885591948123783*x29) - 1.1961863982707928*x22*x28 + 7.5239110000000009*x25*x6/((x26)*(x26))) + x38*(-0.0039872879942359765*T*x22*x35 + 2257.1733000000004*x12*x32*x6/((x33)*(x33)) + x18*(-9.0*x34 + 0.011961863982707929*x36)) - 98570861.633861735 + x24*exp(-1504.7822000000001*x13)/((x16)*(x16)) - x21*exp(-5.0159406666666673*x5)/((x10)*(x10)))/((V)*(V));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-7017821343.4533615*x2 + 9312748475.6423874*x3 + 526367458.21977317);
    double x5 = 1.0/T;
    double x6 = 12.743683058578839*x2;
    double x7 = 16.911048197110393*x3;
    double x8 = x6 - x7 - 0.95583226354999995;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 752.39110000000005*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1504.7822000000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2257.1733000000004*x0*x22*x9;
    double x24 = 0.011961863982707929*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 81324.459777355558*x9;
    double x27 = x17*(-x6 + x7 + 0.95583226354999995);

    result += x0*(-414334403636.14856*x15*x18 - x15*x4 - 138111467878.71619*x16*x18 - x16*x4 + x26*(0.0039872879942359765*x19 + x23 + x25) - x26*(-1698277.1020776303*x22*x27 + x23 + x24 + 3.0*x25 + 3396554.2041552607*x27*exp(x14)/((x21)*(x21)*(x21))) - 276222935757.43237*x18*exp(-2257.1733000000004*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(4138999322.5077281*x2 - 1559515854.1007469);
    double x5 = 12.743683058578839*x2;
    double x6 = 16.911048197110393*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 0.95583226354999995;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 752.39110000000005*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1504.7822000000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 11.274032131406928*x2 - 4.2478943528596123;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0039872879942359765*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2257.1733000000004*x3;
    double x26 = 0.011961863982707929*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 81324.459777355558*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 0.95583226354999995);

    result += x0*x1*x2*(414334403636.14856*x14*x18 - x14*x4 + 138111467878.71619*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(6771.5199000000011*x0*x31 - x26*x30 + 3.0*x27*x32) + 1698277.1020776303*x16*x24 - 3396554.2041552607*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 276222935757.43237*x18*exp(-2257.1733000000004*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(4828832542.925683*x2 - 1299596545.0839555);
    double x4 = 12.743683058578839*x2;
    double x5 = 16.911048197110393*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 0.95583226354999995;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 752.39110000000005*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1504.7822000000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 11.274032131406928*x2 - 4.2478943528596123;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0039872879942359765*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2257.1733000000004*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 0.95583226354999995;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.011961863982707929*x33));
    double x35 = 26.306074973282833*x2 - 7.0798239214326868;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0039872879942359765*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(414334403636.14856*x13*x20 + x13*x3 + 138111467878.71619*x14*x20 + x14*x3 + 81324.459777355558*x15*x34 + 81324.459777355558*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(33.822096394220786*x2 - 12.743683058578837)/pow(x6, 5.0/2.0) + x24*x35 - 1698277.1020776303*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(22.548064262813856*x2 - 8.4957887057192245) + 3396554.2041552607*x37*exp(x12)/((x26)*(x26)*(x26))) + 276222935757.43237*x20*exp(-2257.1733000000004*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = 12.743683058578839*x1;
    double x4 = 16.911048197110393*x2;
    double x5 = x3 - x4 - 0.95583226354999995;
    double x6 = sqrt(x5);
    double x7 = 2.5079703333333336*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = pow(x5, -3.0/2.0);
    double x11 = pow(V, -2);
    double x12 = 11.274032131406928*x1 - 4.2478943528596123;
    double x13 = x11*((x12)*(x12)*(x12));
    double x14 = x10*x13;
    double x15 = 5.0159406666666673*x6;
    double x16 = exp(-x15)/((x9)*(x9));
    double x17 = -x3 + x4 + 0.95583226354999995;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 611877.99748790311*x19;
    double x21 = x8/x9;
    double x22 = 1.0/x6;
    double x23 = 87.686916577609438*x1 - 18.87953045715383;
    double x24 = x1*x22*x23;
    double x25 = 243973.37933206669*x21;
    double x26 = 1.0/T;
    double x27 = x26*x6;
    double x28 = 752.39110000000005*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 243973.37933206666*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1504.7822000000001*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 183563399.24637091*x26;
    double x38 = x19*x37;
    double x39 = 22.548064262813856*x1 - 8.4957887057192245;
    double x40 = ((x12)*(x12));
    double x41 = x11*x40;
    double x42 = x18*x39*x41;
    double x43 = 611877.99748790311*x42;
    double x44 = (33.822096394220786*x1 - 12.743683058578837)/pow(x5, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 26.306074973282833*x1 - 7.0798239214326868;
    double x47 = x12*x2;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = 611877.99748790311*x49;
    double x51 = x16*x50;
    double x52 = x47*(52.612149946565665*x1 - 14.159647842865374);
    double x53 = x10*x47;
    double x54 = 731920.13799620001*x46*x53;
    double x55 = x21*x50;
    double x56 = x37*x42;
    double x57 = x37*x49;
    double x58 = x36*x57;
    double x59 = x31*x57;
    double x60 = 81324.459777355558*x22;
    double x61 = x23*x60;
    double x62 = exp(x7);
    double x63 = x62 - 1;
    double x64 = 1.0/x63;
    double x65 = Debye(x7);
    double x66 = x22*x65;
    double x67 = -3*x64 + 1.1961863982707928*x66;
    double x68 = x1*x67;
    double x69 = exp(x28);
    double x70 = x69 - 1;
    double x71 = 1.0/x70;
    double x72 = Debye(x28);
    double x73 = T*x22*x72;
    double x74 = -3*x71 + 0.0039872879942359765*x73;
    double x75 = x1*x74;
    double x76 = 81324.459777355558*x45;
    double x77 = x46*x67;
    double x78 = 243973.37933206669*x53;
    double x79 = x46*x74;
    double x80 = 1.1961863982707928*x65;
    double x81 = x10*x80;
    double x82 = x62/((x63)*(x63));
    double x83 = 7.5239110000000009*x82;
    double x84 = x22*x83;
    double x85 = x49*(-9.0*x64 + 3.5885591948123783*x66) - x81 + x84;
    double x86 = 162648.91955471112*x85;
    double x87 = x22*x48;
    double x88 = 0.0039872879942359765*T*x72;
    double x89 = x10*x88;
    double x90 = x69/((x70)*(x70));
    double x91 = 2257.1733000000004*x26*x90;
    double x92 = x22*x91;
    double x93 = x49*(-9.0*x71 + 0.011961863982707929*x73) - x89 + x92;
    double x94 = 162648.91955471112*x93;
    double x95 = x1*x12*x44;
    double x96 = x1*x40;
    double x97 = x49*x96;
    double x98 = x10*x96;
    double x99 = 3.0*x49;
    double x100 = 3.0*x18;
    double x101 = x12*x39;
    double x102 = x96*x99;
    double x103 = x47*x60;
    double x104 = x33*x97;
    double x105 = x100*x75;

    result += (262855631.02363127*x1 + x103*(x101*x105 + x102*x93 - 1698277.1020776303*x104*x90 + 3396554.2041552607*x104*exp(x35)/((x70)*(x70)*(x70)) + x105*x40 + x46*x89 - x46*x92 - x79*x99 + x88*x95 - x91*x98) - x103*(x100*x101*x68 + x100*x67*x96 + x102*x85 + x46*x81 - x46*x84 - x77*x99 + x80*x95 - 18.869745578640337*x82*x97 - x83*x98 + 37.739491157280675*x97*exp(x15)/((x63)*(x63)*(x63))) - 4603715.5959572066*x14*x16 - 1534571.865319069*x14*x21 + x14*x86 - x14*x94 - 3069143.730638138*x14*exp(-7.5239110000000009*x6)/((x9)*(x9)*(x9)) - x16*x20 - x16*x43 - 2143792595.0601654*x2 - x20*x21 - x21*x43 - x21*x54 - x24*x25 + x24*x32 - x25*x45 + 138111467878.71619*x31*x34 + x31*x38 + x31*x54 + x31*x56 + x32*x45 + 414334403636.14856*x34*x36 + x36*x38 + x36*x56 + x48*x51 + x48*x55 - x48*x58 - x48*x59 + x51*x52 + x52*x55 - x52*x58 - x52*x59 - x61*x68 + x61*x75 - x67*x76 + x74*x76 - x77*x78 + x78*x79 + x86*x87 - x87*x94 + 276222935757.43237*x34*exp(-2257.1733000000004*x27)/((x30)*(x30)*(x30)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*10.061;
static const double Vmax = 1.15*10.061;
static double V = 0.9*10.061;

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



const char *Anorthite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Anorthite_slb_em_coder_calib_name(void) {
    return "Anorthite_slb_em";
}

const char *Anorthite_slb_em_coder_calib_formula(void) {
    return "CaAl2Si2O8";
}

const double Anorthite_slb_em_coder_calib_mw(void) {
    return 278.20928;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,8.0,0.0,0.0,0.0,
        0.0,2.0,2.0,0.0,0.0,0.0,
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

const double *Anorthite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Anorthite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Anorthite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Anorthite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Anorthite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Anorthite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Anorthite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Anorthite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Anorthite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Anorthite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Anorthite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Anorthite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Anorthite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Anorthite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Anorthite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Anorthite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Anorthite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Anorthite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Anorthite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Anorthite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Anorthite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Anorthite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Anorthite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Anorthite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Anorthite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Anorthite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Anorthite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Anorthite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Anorthite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Anorthite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Anorthite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

