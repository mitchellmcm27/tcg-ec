
static char *identifier = "Seifertite_slb_em.emml:6aa8147b013f3f58dd237f7b16c71faa88fd92f1:Fri May 26 02:51:14 2023";



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
    double x3 = sqrt(3.70558619230939*x1 + 1.2605567001450875*x2 - 2.8393411597999991);
    double x4 = 3.763153*x3;
    double x5 = 1128.9458999999999*x3/T;

    result += 74.830163563379159*T*log(1 - exp(-x5)) - 24.943387854459719*T*Debye(x5) - 12239700.060537552*x1 + 7442801.2293260405*x2 - 22449.049069013745*log(1 - exp(-x4)) + 7483.0163563379156*Debye(x4) + 4196073.9572256776 + 78129.310714971871/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(1.2605567001450875*pow(x0, 4.0/3.0) + 3.70558619230939*pow(x0, 2.0/3.0) - 2.8393411597999991);
    double x2 = x1/T;
    double x3 = 1128.9458999999999*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -84479.206351206289*x2*x5/x6 + 28159.735450402095*x2*(-0.0026573461137508895*T*x4/x1 + 3/(exp(x3) - 1)) - 24.943387854459719*x4 + 74.830163563379159*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(3.70558619230939*x1 + 1.2605567001450875*x3 - 2.8393411597999991);
    double x6 = 3.763153*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-1.2351953974364633*x2 - 0.84037113343005831*x4);
    double x10 = 1128.9458999999999*x5/T;
    double x11 = exp(-x10);
    double x12 = 28159.735450402095*x9;

    result += 84479.206351206289*x11*x9/(-x11 + 1) + x12*(-0.79720383412526674*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.0026573461137508895*T*x8*Debye(x10) + 3/(exp(x10) - 1)) + 8159800.0403583683*x2 - 9923734.97243472*x4 - 84479.206351206274*x7*x9/(-x7 + 1) - 156258.62142994374/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(3.70558619230939*x2 + 1.2605567001450875*x3 - 2.8393411597999991);
    double x5 = x0*x4;
    double x6 = 1128.9458999999999*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(353410847.35524058*x2 + 120222385.45204662*x3 - 270794933.14663881);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-2257.8917999999999*x5)/((x8)*(x8)) + 28159.735450402095*x4*(x0*(0.0079720383412526692*T*x11 - 9.0/x13) + 0.0026573461137508895*x11 - 3386.8377*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 1.2605567001450875*pow(x1, 4.0/3.0) + 3.70558619230939*x2 - 2.8393411597999991;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 1128.9458999999999*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 95372453.645448297*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(0.84037113343005831*x2 + 1.2351953974364633)*(28159.735450402095*x4*(-3386.8377*x0*x11*x12/((x13)*(x13)) + 0.0026573461137508895*x10/pow(x3, 3.0/2.0) + (0.0079720383412526692*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-2257.8917999999999*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 3.70558619230939*x2 + 1.2605567001450875*x3 - 2.8393411597999991;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x2*(1.9608659780034694*x2 + 2.0586589957274386);
    double x8 = x6*x7;
    double x9 = 3.763153*x5;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = x10/x11;
    double x13 = 84479.206351206274*x12;
    double x14 = 1.0/x4;
    double x15 = x3*((0.84037113343005831*x2 + 1.2351953974364633)*(0.84037113343005831*x2 + 1.2351953974364633));
    double x16 = x14*x15;
    double x17 = 317908.17881816096*x16;
    double x18 = pow(x4, -3.0/2.0);
    double x19 = x15*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 1128.9458999999999*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 84479.206351206289*x25;
    double x27 = 95372453.645448297*x16*x20;
    double x28 = exp(x9);
    double x29 = x28 - 1;
    double x30 = 1.0/x29;
    double x31 = Debye(x9);
    double x32 = x31*x6;
    double x33 = -3*x30 + 0.79720383412526674*x32;
    double x34 = 28159.735450402095*x6;
    double x35 = x34*x7;
    double x36 = 28159.735450402095*x19;
    double x37 = exp(x22);
    double x38 = x37 - 1;
    double x39 = 1.0/x38;
    double x40 = Debye(x22);
    double x41 = T*x40*x6;
    double x42 = -3*x39 + 0.0026573461137508895*x41;
    double x43 = x15*x34;

    result += x0*(468775.86428983125*x0 + x12*x17 + x13*x19 - x13*x8 - x19*x26 - 13599666.733930614*x2 - x25*x27 + x26*x8 + 23155381.602347679*x3 - x33*x35 + x33*x36 + x35*x42 - x36*x42 + x43*(x14*(-9.0*x30 + 2.3916115023758002*x32) + 0.79720383412526674*x18*x31 - 11.289459000000001*x28*x6/((x29)*(x29))) - x43*(0.0026573461137508895*T*x18*x40 + x14*(-9.0*x39 + 0.0079720383412526692*x41) - 3386.8377*x20*x37*x6/((x38)*(x38))) - x27*exp(-2257.8917999999999*x21)/((x24)*(x24)) + x17*exp(-7.5263059999999999*x5)/((x11)*(x11)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(1060232542.0657215*x2 + 360667156.35613984*x3 - 812384799.43991637);
    double x5 = 1.0/T;
    double x6 = 3.70558619230939*x2 + 1.2605567001450875*x3 - 2.8393411597999991;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 1128.9458999999999*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 2257.8917999999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = -3386.8377*x0*x20*x7;
    double x22 = 0.0079720383412526692*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 28159.735450402095*x7;
    double x25 = x15*x6;

    result += x0*(-323011021547.90674*x13*x16 + x13*x4 - 107670340515.9689*x14*x16 + x14*x4 + x24*(0.0026573461137508895*x17 + x21 + x23) - x24*(3823556.5353804301*x20*x25 + x21 + x22 + 3.0*x23 - 7647113.0707608601*x25*exp(x12)/((x19)*(x19)*(x19))) - 215340681031.93781*x16*exp(-3386.8377*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(160296513.93606216*x2 + 235607231.57016036);
    double x5 = 1.2605567001450875*pow(x1, 4.0/3.0) + 3.70558619230939*x2 - 2.8393411597999991;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 1128.9458999999999*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 2257.8917999999999*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 0.84037113343005831*x2 + 1.2351953974364633;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0026573461137508895*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 3386.8377*x22*x3;
    double x24 = 0.0079720383412526692*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 28159.735450402095*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-323011021547.90674*x12*x16 + x12*x4 - 107670340515.9689*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-10160.5131*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 3823556.5353804301*x15*x22 - 7647113.0707608601*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 215340681031.93781*x16*exp(-3386.8377*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(187012599.59207252*x2 + 196339359.64180028);
    double x4 = 1.2605567001450875*pow(x1, 4.0/3.0) + 3.70558619230939*x2 - 2.8393411597999991;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 1128.9458999999999*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 2257.8917999999999*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 0.84037113343005831*x2 + 1.2351953974364633;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0026573461137508895*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 3386.8377*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.0079720383412526692*x30));
    double x32 = 1.9608659780034694*x2 + 2.0586589957274386;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0026573461137508895*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-323011021547.90674*x11*x18 + x11*x3 - 107670340515.9689*x12*x18 + x12*x3 + 28159.735450402095*x13*x31 - 215340681031.93781*x18*exp(-3386.8377*x6)/((x9)*(x9)*(x9)) - 28159.735450402095*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(2.5211134002901749*x2 + 3.70558619230939)/pow(x4, 5.0/2.0) - x22*x32 + 3823556.5353804301*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(1.6807422668601166*x2 + 2.4703907948729267) - 7647113.0707608601*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 3.70558619230939*x2 + 1.2605567001450875*x3 - 2.8393411597999991;
    double x5 = sqrt(x4);
    double x6 = 3.763153*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 0.84037113343005831*x2 + 1.2351953974364633;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = 7.5263059999999999*x5;
    double x14 = exp(-x13)/((x8)*(x8));
    double x15 = pow(x4, -2);
    double x16 = 317908.17881816096*x15;
    double x17 = x11*x16;
    double x18 = x7/x8;
    double x19 = 84479.206351206274*x18;
    double x20 = 1.0/x5;
    double x21 = 6.536219926678231*x2 + 5.4897573219398357;
    double x22 = x2*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 1128.9458999999999*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 84479.206351206289*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 2257.8917999999999*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = x11*x15;
    double x35 = 95372453.645448297*x23;
    double x36 = x33*x35;
    double x37 = x28*x35;
    double x38 = 1.6807422668601166*x2 + 2.4703907948729267;
    double x39 = ((x10)*(x10));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = x16*x41;
    double x43 = (2.5211134002901749*x2 + 3.70558619230939)/pow(x4, 5.0/2.0);
    double x44 = x40*x43;
    double x45 = 1.9608659780034694*x2 + 2.0586589957274386;
    double x46 = x10*x3;
    double x47 = x45*x46;
    double x48 = 1.0/x4;
    double x49 = 317908.17881816096*x48;
    double x50 = x14*x49;
    double x51 = x46*(3.9217319560069388*x2 + 4.1173179914548772);
    double x52 = x46*x9;
    double x53 = x45*x52;
    double x54 = x18*x49;
    double x55 = x15*x41;
    double x56 = x36*x48;
    double x57 = x37*x48;
    double x58 = exp(x6);
    double x59 = x58 - 1;
    double x60 = 1.0/x59;
    double x61 = Debye(x6);
    double x62 = x20*x61;
    double x63 = -3*x60 + 0.79720383412526674*x62;
    double x64 = x2*x63;
    double x65 = 28159.735450402095*x20;
    double x66 = x21*x65;
    double x67 = exp(x25);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x25);
    double x71 = T*x20*x70;
    double x72 = -3*x69 + 0.0026573461137508895*x71;
    double x73 = x2*x72;
    double x74 = 28159.735450402095*x44;
    double x75 = x45*x63;
    double x76 = 84479.206351206289*x52;
    double x77 = x45*x72;
    double x78 = 0.79720383412526674*x61;
    double x79 = x78*x9;
    double x80 = x58/((x59)*(x59));
    double x81 = 11.289459000000001*x80;
    double x82 = x20*x81;
    double x83 = x48*(-9.0*x60 + 2.3916115023758002*x62) + x79 - x82;
    double x84 = 56319.47090080419*x83;
    double x85 = x20*x47;
    double x86 = 0.0026573461137508895*T*x70;
    double x87 = x86*x9;
    double x88 = x67/((x68)*(x68));
    double x89 = 3386.8377*x23*x88;
    double x90 = x20*x89;
    double x91 = x48*(-9.0*x69 + 0.0079720383412526692*x71) + x87 - x90;
    double x92 = 56319.47090080419*x91;
    double x93 = x10*x2*x43;
    double x94 = x2*x39;
    double x95 = x48*x94;
    double x96 = x9*x94;
    double x97 = 3.0*x48;
    double x98 = 3.0*x15;
    double x99 = x64*x98;
    double x100 = x10*x38;
    double x101 = x94*x97;
    double x102 = x46*x65;
    double x103 = x30*x95;

    result += (-1875103.457159325*x0 + x102*(x100*x99 + x101*x83 + x39*x99 - x45*x79 + x45*x82 - x75*x97 + x78*x93 + 42.483961504227004*x80*x95 - x81*x96 - 84.967923008454008*x95*exp(x13)/((x59)*(x59)*(x59))) - x102*(x100*x73*x98 + x101*x91 + 3823556.5353804301*x103*x88 - 7647113.0707608601*x103*exp(x32)/((x68)*(x68)*(x68)) - x45*x87 + x45*x90 + x72*x94*x98 - x77*x97 + x86*x93 - x89*x96) + 3589011.3505322961*x12*x14 + 1196337.1168440988*x12*x18 + x12*x84 - x12*x92 + 2392674.2336881976*x12*exp(-11.289459000000001*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x42 + x17*x18 + x18*x42 - 253437.61905361881*x18*x53 + x19*x22 + x19*x44 + 36265777.957148299*x2 - x22*x29 - 107670340515.9689*x28*x31 + 253437.61905361887*x28*x53 - x29*x44 - 77184605.341158926*x3 - 323011021547.90674*x31*x33 - x34*x36 - x34*x37 - x36*x55 - x37*x55 - x47*x50 - x47*x54 + x47*x56 + x47*x57 - x50*x51 - x51*x54 + x51*x56 + x51*x57 + x63*x74 + x64*x66 - x66*x73 - x72*x74 - x75*x76 + x76*x77 - x84*x85 + x85*x92 - 215340681031.93781*x31*exp(-3386.8377*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*1.367;
static const double Vmax = 1.15*1.367;
static double V = 0.9*1.367;

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



const char *Seifertite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Seifertite_slb_em_coder_calib_name(void) {
    return "Seifertite_slb_em";
}

const char *Seifertite_slb_em_coder_calib_formula(void) {
    return "SiO2";
}

const double Seifertite_slb_em_coder_calib_mw(void) {
    return 60.0843;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,2.0,0.0,0.0,0.0,
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
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0
    };

const double *Seifertite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Seifertite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Seifertite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Seifertite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Seifertite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Seifertite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Seifertite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Seifertite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Seifertite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Seifertite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Seifertite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Seifertite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Seifertite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Seifertite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Seifertite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Seifertite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Seifertite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Seifertite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Seifertite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Seifertite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Seifertite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Seifertite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Seifertite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

