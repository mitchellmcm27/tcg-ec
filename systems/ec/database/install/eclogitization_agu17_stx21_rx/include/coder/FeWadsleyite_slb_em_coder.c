
static char *identifier = "FeWadsleyite_slb_em.emml:dbc67b38acba30aeb19279e965d0adf8de37feb6:Tue Aug 22 19:33:24 2023";



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
    double x3 = sqrt(16.218107763452309*x1 - 8.8275864332619456*x2 - 3.8819711662999996);
    double x4 = 2.1578186333333336*x3;
    double x5 = 647.34559000000002*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 26.76322271834383*T - 38846178.987580091*x1 + 45996557.700691968*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 6259402.8940333277 + 9146793.4508254919/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-8.8275864332619456*pow(x0, 4.0/3.0) + 16.218107763452309*pow(x0, 2.0/3.0) - 3.8819711662999996);
    double x2 = x1/T;
    double x3 = 647.34559000000002*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -113028.94489070842*x2*x5/x6 + 37676.314963569472*x2*(-0.0046343097818894539*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(16.218107763452309*x1 - 8.8275864332619456*x3 - 3.8819711662999996);
    double x6 = 2.1578186333333336*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-5.4060359211507691*x2 + 5.8850576221746298*x4);
    double x10 = 647.34559000000002*x5/T;
    double x11 = exp(-x10);

    result += 113028.94489070842*x11*x9/(-x11 + 1) + 25897452.658386726*x2 - 61328743.600922622*x4 - 113028.94489070843*x7*x9/(-x7 + 1) + 37676.314963569479*x9*(-1.390292934566836*x8*Debye(x6) + 3/(exp(x6) - 1)) - 37676.314963569472*x9*(-0.0046343097818894539*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 18293586.901650984/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(16.218107763452309*x2 - 8.8275864332619456*x3 - 3.8819711662999996);
    double x5 = x0*x4;
    double x6 = 647.34559000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1186659305.2047389*x2 + 645903809.26779211*x3 + 284039129.23845297);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1294.69118*x5)/((x8)*(x8)) - 37676.314963569472*x4*(x0*(0.013902929345668362*T*x11 - 9.0/x13) + 0.0046343097818894539*x11 - 1942.0367700000002*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 16.218107763452309*x2;
    double x4 = 8.8275864332619456*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 3.8819711662999996;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 647.34559000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 73168789.017353132*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(5.8850576221746298*x2 - 5.4060359211507691)*(-37676.314963569472*x6*(1942.0367700000002*x0*x13*x14/((x15)*(x15)) - 0.0046343097818894539*x12/pow(x5, 3.0/2.0) + (0.013902929345668362*x12*x13 - 9.0/x15)/(-x3 + x4 + 3.8819711662999996)) + x11*x9/x10 + x11*exp(-1294.69118*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 16.218107763452309*x2;
    double x5 = 8.8275864332619456*x3;
    double x6 = x4 - x5 - 3.8819711662999996;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(13.73180111840747*x2 - 9.0100598685846158);
    double x10 = x8*x9;
    double x11 = 2.1578186333333336*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 113028.94489070843*x14;
    double x16 = 1.0/(-x4 + x5 + 3.8819711662999996);
    double x17 = x3*((5.8850576221746298*x2 - 5.4060359211507691)*(5.8850576221746298*x2 - 5.4060359211507691));
    double x18 = x16*x17;
    double x19 = 243895.96339117715*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 647.34559000000002*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 113028.94489070842*x27;
    double x29 = 73168789.017353132*x18*x22;
    double x30 = exp(x11);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x11);
    double x34 = x33*x8;
    double x35 = -113028.94489070843*x32 + 52381.114494365407*x34;
    double x36 = exp(x24);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = T*Debye(x24);
    double x40 = x39*x8;
    double x41 = -3*x38 + 0.0046343097818894539*x40;
    double x42 = 37676.314963569472*x8;

    result += x0*(54880760.704952955*x0 + x10*x15 - x10*x28 + x10*x35 - x14*x19 + x15*x21 + x17*x42*(x16*(-9.0*x38 + 0.013902929345668362*x40) - 0.0046343097818894539*x20*x39 + 1942.0367700000002*x22*x36*x8/((x37)*(x37))) - 37676.314963569479*x17*x8*(x16*(-9.0*x32 + 4.170878803700508*x34) - 1.390292934566836*x20*x33 + 6.4734559000000012*x30*x8/((x31)*(x31))) - 43162421.097311214*x2 - x21*x28 + x21*x35 - 37676.314963569472*x21*x41 + x27*x29 + 143100401.73548612*x3 - x41*x42*x9 + x29*exp(-1294.69118*x23)/((x26)*(x26)) - x19*exp(-4.3156372666666671*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-3559977915.6142168*x2 + 1937711427.8033767*x3 + 852117387.71535885);
    double x5 = 1.0/T;
    double x6 = 16.218107763452309*x2;
    double x7 = 8.8275864332619456*x3;
    double x8 = x6 - x7 - 3.8819711662999996;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 647.34559000000002*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1294.69118*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.013902929345668362*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 37676.314963569472*x9;
    double x27 = x17*(-x6 + x7 + 3.8819711662999996);

    result += x0*(-142096478688.07196*x15*x18 - x15*x4 - 47365492896.023987*x16*x18 - x16*x4 + x26*(0.0046343097818894539*x19 - 1942.0367700000002*x23 + x25) - x26*(-1257168.9386773445*x22*x27 - 1942.0367699999997*x23 + x24 + 3.0*x25 + 2514337.8773546889*x27*exp(x14)/((x21)*(x21)*(x21))) - 94730985792.047974*x18*exp(-1942.0367700000002*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(861205079.02372277*x2 - 791106203.46982586);
    double x5 = 16.218107763452309*x2;
    double x6 = 8.8275864332619456*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 3.8819711662999996;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 647.34559000000002*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1294.69118*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 5.8850576221746298*x2 - 5.4060359211507691;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0046343097818894539*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 1942.0367700000002*x3;
    double x26 = 0.013902929345668362*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 37676.314963569472*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 3.8819711662999996);

    result += x0*x1*x2*(142096478688.07196*x14*x18 - x14*x4 + 47365492896.023987*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(5826.11031*x0*x31 - x26*x30 + 3.0*x27*x32) + 1257168.9386773445*x16*x24 - 2514337.8773546889*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 94730985792.047974*x18*exp(-1942.0367700000002*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(1004739258.86101*x2 - 659255169.5581882);
    double x4 = 16.218107763452309*x2;
    double x5 = 8.8275864332619456*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 3.8819711662999996;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 647.34559000000002*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1294.69118*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 5.8850576221746298*x2 - 5.4060359211507691;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0046343097818894539*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 1942.0367700000002*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 3.8819711662999996;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.013902929345668362*x33));
    double x35 = 13.73180111840747*x2 - 9.0100598685846158;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0046343097818894539*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(142096478688.07196*x13*x20 + x13*x3 + 47365492896.023987*x14*x20 + x14*x3 + 37676.314963569472*x15*x34 + 37676.314963569472*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(17.655172866523891*x2 - 16.218107763452309)/pow(x6, 5.0/2.0) + x24*x35 - 1257168.9386773445*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(11.77011524434926*x2 - 10.812071842301538) + 2514337.8773546889*x37*exp(x12)/((x26)*(x26)*(x26))) + 94730985792.047974*x20*exp(-1942.0367700000002*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 16.218107763452309*x2;
    double x5 = 8.8275864332619456*x3;
    double x6 = x4 - x5 - 3.8819711662999996;
    double x7 = sqrt(x6);
    double x8 = 2.1578186333333336*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 5.8850576221746298*x2 - 5.4060359211507691;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.3156372666666671*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 3.8819711662999996;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 243895.96339117715*x16;
    double x21 = x9/x10;
    double x22 = 243895.96339117715*x21;
    double x23 = 1.0/x7;
    double x24 = 45.772670394691559*x2 - 24.026826316225641;
    double x25 = x2*x23*x24;
    double x26 = 113028.94489070843*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 647.34559000000002*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 113028.94489070842*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1294.69118*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 73168789.017353132*x27;
    double x39 = x19*x38;
    double x40 = 11.77011524434926*x2 - 10.812071842301538;
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x18*x40*x42;
    double x44 = (17.655172866523891*x2 - 16.218107763452309)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 13.73180111840747*x2 - 9.0100598685846158;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = x20*x49;
    double x51 = x47*(27.463602236814939*x2 - 18.020119737169232);
    double x52 = x11*x48;
    double x53 = x22*x49;
    double x54 = x38*x43;
    double x55 = x38*x49;
    double x56 = x37*x55;
    double x57 = x32*x55;
    double x58 = 37676.314963569479*x23;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.390292934566836*x63;
    double x65 = x2*x64;
    double x66 = exp(x29);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x29);
    double x70 = T*x23*x69;
    double x71 = -3*x68 + 0.0046343097818894539*x70;
    double x72 = x2*x71;
    double x73 = 37676.314963569472*x23;
    double x74 = x46*x64;
    double x75 = x11*x47;
    double x76 = x46*x71;
    double x77 = 1.390292934566836*x62;
    double x78 = x11*x77;
    double x79 = x59/((x60)*(x60));
    double x80 = 6.4734559000000012*x79;
    double x81 = x23*x80;
    double x82 = x49*(-9.0*x61 + 4.170878803700508*x63) - x78 + x81;
    double x83 = 75352.629927138958*x82;
    double x84 = x23*x48;
    double x85 = 0.0046343097818894539*T*x69;
    double x86 = x11*x85;
    double x87 = x66/((x67)*(x67));
    double x88 = 1942.0367700000002*x27*x87;
    double x89 = x23*x88;
    double x90 = x49*(-9.0*x68 + 0.013902929345668362*x70) - x86 + x89;
    double x91 = 75352.629927138943*x90;
    double x92 = x12*x2*x44;
    double x93 = x2*x41;
    double x94 = x49*x93;
    double x95 = x11*x93;
    double x96 = 3.0*x49;
    double x97 = 3.0*x18;
    double x98 = x12*x40;
    double x99 = x93*x96;
    double x100 = x34*x94;
    double x101 = x72*x97;

    result += (-219523042.81981182*x0 - 1578849.7632007999*x14*x16 - 526283.25440026668*x14*x21 + x14*x83 - x14*x91 - x19*x20 - x19*x22 + 115099789.5928299*x2 - x20*x43 - 339086.83467212529*x21*x52 - x22*x43 - x24*x58*x65 + x24*x72*x73 - x25*x26 + x25*x33 - x26*x45 - 477001339.11828709*x3 + 47365492896.023987*x32*x35 + x32*x39 + 339086.83467212523*x32*x52 + x32*x54 + x33*x45 + 142096478688.07196*x35*x37 + x37*x39 + x37*x54 - 37676.314963569479*x45*x64 + 37676.314963569472*x45*x71 - x47*x58*(x46*x78 - x46*x81 + x64*x93*x97 + x65*x97*x98 - x74*x96 + x77*x92 - 13.968543763081607*x79*x94 - x80*x95 + x82*x99 + 27.937087526163214*x94*exp(x15)/((x60)*(x60)*(x60))) + x47*x73*(-1257168.9386773445*x100*x87 + 2514337.8773546889*x100*exp(x36)/((x67)*(x67)*(x67)) + x101*x41 + x101*x98 + x46*x86 - x46*x89 - x76*x96 + x85*x92 - x88*x95 + x90*x99) + x48*x50 + x48*x53 - x48*x56 - x48*x57 + x50*x51 + x51*x53 - x51*x56 - x51*x57 - 113028.94489070843*x74*x75 + 113028.94489070842*x75*x76 + x83*x84 - x84*x91 + 94730985792.047974*x35*exp(-1942.0367700000002*x28)/((x31)*(x31)*(x31)) - 1052566.5088005334*x14*exp(-6.4734559000000012*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*4.28;
static const double Vmax = 1.15*4.28;
static double V = 0.9*4.28;

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



const char *FeWadsleyite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *FeWadsleyite_slb_em_coder_calib_name(void) {
    return "FeWadsleyite_slb_em";
}

const char *FeWadsleyite_slb_em_coder_calib_formula(void) {
    return "Fe2SiO4";
}

const double FeWadsleyite_slb_em_coder_calib_mw(void) {
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

const double *FeWadsleyite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double FeWadsleyite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double FeWadsleyite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double FeWadsleyite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double FeWadsleyite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double FeWadsleyite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double FeWadsleyite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double FeWadsleyite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double FeWadsleyite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double FeWadsleyite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double FeWadsleyite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double FeWadsleyite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int FeWadsleyite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **FeWadsleyite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **FeWadsleyite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void FeWadsleyite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int FeWadsleyite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double FeWadsleyite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int FeWadsleyite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double FeWadsleyite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double FeWadsleyite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

