
static char *identifier = "NaNAL_stx21_em.emml:0f6ad14ea9e0af49482b51e3e20cdf33ea13764b:Sat Mar  9 00:08:44 2024";



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
    double x3 = sqrt(-16.680099032862614*x1 + 93.184318685221172*x2 + 0.54777557845000002);
    double x4 = 2.8344277666666668*x3;
    double x5 = 850.32833000000005*x3/T;

    result += 523.81114494365409*T*log(1 - exp(-x5)) - 174.60371498121802*T*Debye(x5) - 188284546.21692681*x1 + 318133361.59967774*x2 - 157143.34348309625*log(1 - exp(-x4)) + 52381.114494365414*Debye(x4) + 15537243.005017532 + 479047685.95908189/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(93.184318685221172*pow(x0, 4.0/3.0) - 16.680099032862614*pow(x0, 2.0/3.0) + 0.54777557845000002);
    double x2 = x1/T;
    double x3 = 850.32833000000005*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -445411.45611532533*x2*x5/x6 + 148470.48537177511*x2*(-0.0035280489831498382*T*x4/x1 + 3/(exp(x3) - 1)) - 174.60371498121802*x4 + 523.81114494365409*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-16.680099032862614*x1 + 93.184318685221172*x3 + 0.54777557845000002);
    double x6 = 2.8344277666666668*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(5.5600330109542044*x2 - 62.122879123480779*x4);
    double x10 = 850.32833000000005*x5/T;
    double x11 = exp(-x10);

    result += 445411.45611532533*x11*x9/(-x11 + 1) + 125523030.81128454*x2 - 424177815.46623695*x4 - 445411.45611532545*x7*x9/(-x7 + 1) + 148470.48537177514*x9*(-1.0584146949449513*x8*Debye(x6) + 3/(exp(x6) - 1)) - 148470.48537177511*x9*(-0.0035280489831498382*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 958095371.91816378/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-16.680099032862614*x2 + 93.184318685221172*x3 + 0.54777557845000002);
    double x5 = x0*x4;
    double x6 = 850.32833000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-6317520448.7173347*x2 + 35293186067.651711*x3 + 207467798.08368689);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1700.6566600000001*x5)/((x8)*(x8)) + 148470.48537177511*x4*(x0*(0.010584146949449516*T*x11 - 9.0000000000000018/x13) + 0.0035280489831498382*x11 - 2550.9849899999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 93.184318685221172*pow(x1, 4.0/3.0) - 16.680099032862614*x2 + 0.54777557845000002;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 850.32833000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 378745979.64141291*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(62.122879123480779*x2 - 5.5600330109542044)*(148470.48537177511*x4*(-2550.9849899999999*x0*x11*x12/((x13)*(x13)) + 0.0035280489831498382*x10/pow(x3, 3.0/2.0) + (0.010584146949449516*x10*x11 - 9.0000000000000018/x13)/x3) + x7*x9/x8 + x9*exp(-1700.6566600000001*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -16.680099032862614*x2 + 93.184318685221172*x3 + 0.54777557845000002;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x2*x6*(144.95338462145514*x2 - 9.2667216849236738);
    double x8 = 2.8344277666666668*x5;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 445411.45611532545*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((62.122879123480779*x2 - 5.5600330109542044)*(62.122879123480779*x2 - 5.5600330109542044));
    double x15 = x13*x14;
    double x16 = 1262486.59880471*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 1.0/T;
    double x20 = x19*x5;
    double x21 = 850.32833000000005*x20;
    double x22 = exp(-x21);
    double x23 = -x22 + 1;
    double x24 = x22/x23;
    double x25 = 445411.45611532533*x24;
    double x26 = 378745979.64141291*x15*x19;
    double x27 = exp(x8);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x8);
    double x31 = x30*x6;
    double x32 = -445411.45611532545*x29 + 157143.34348309625*x31;
    double x33 = exp(x21);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x21);
    double x37 = T*x36*x6;
    double x38 = -445411.45611532533*x35 + 523.81114494365409*x37;
    double x39 = x14*x6;

    result += x0*(2874286115.7544913*x0 + x11*x16 + x12*x18 - x12*x7 - x18*x25 + x18*x32 - x18*x38 - 209205051.3521409*x2 - x24*x26 + x25*x7 + 989748236.08788621*x3 - x32*x7 + x38*x7 + 148470.48537177514*x39*(x13*(-9.0*x29 + 3.175244084834854*x31) + 1.0584146949449513*x17*x30 - 8.5032832999999997*x27*x6/((x28)*(x28))) - 148470.48537177511*x39*(0.0035280489831498382*T*x17*x36 + x13*(-9.0000000000000018*x35 + 0.010584146949449516*x37) - 2550.9849899999999*x19*x33*x6/((x34)*(x34))) - x26*exp(-1700.6566600000001*x20)/((x23)*(x23)) + x16*exp(-5.6688555333333337*x5)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-18952561346.152004*x2 + 105879558202.95512*x3 + 622403394.25106061);
    double x5 = 1.0/T;
    double x6 = -16.680099032862614*x2 + 93.184318685221172*x3 + 0.54777557845000002;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 850.32833000000005*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1700.6566600000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.010584146949449516*x17;
    double x23 = x5*(T*x22 - 9.0000000000000018/x19);
    double x24 = 148470.48537177511*x7;
    double x25 = x15*x6;

    result += x0*(-966175309088.08997*x13*x16 + x13*x4 - 322058436362.69666*x14*x16 + x14*x4 + x24*(0.0035280489831498382*x17 - 2550.9849899999999*x21 + x23) - x24*(2169174.8064017668*x20*x25 - 2550.9849900000008*x21 + x22 + 3.0000000000000004*x23 - 4338349.6128035337*x25*exp(x12)/((x19)*(x19)*(x19))) - 644116872725.39331*x16*exp(-2550.9849899999999*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(47057581423.535614*x2 - 4211680299.1448898);
    double x5 = 93.184318685221172*pow(x1, 4.0/3.0) - 16.680099032862614*x2 + 0.54777557845000002;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 850.32833000000005*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1700.6566600000001*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 62.122879123480779*x2 - 5.5600330109542044;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0035280489831498382*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2550.9849899999999*x22*x3;
    double x24 = 0.010584146949449516*T*x18;
    double x25 = x17*x24 - 9.0000000000000018/x21;
    double x26 = x0*x25;
    double x27 = 148470.48537177511*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-966175309088.08997*x12*x16 + x12*x4 - 322058436362.69666*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7652.9549700000007*x0*x17*x22 + x24*x28 + 3.0000000000000004*x25*x29) + 2169174.8064017668*x15*x22 - 4338349.6128035337*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 644116872725.39331*x16*exp(-2550.9849899999999*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(54900511660.791542*x2 - 3509733582.6207414);
    double x4 = 93.184318685221172*pow(x1, 4.0/3.0) - 16.680099032862614*x2 + 0.54777557845000002;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 850.32833000000005*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1700.6566600000001*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 62.122879123480779*x2 - 5.5600330109542044;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0035280489831498382*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2550.9849899999999*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0000000000000018*x29 + 0.010584146949449516*x30));
    double x32 = 144.95338462145514*x2 - 9.2667216849236738;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0035280489831498382*x30;
    double x37 = 3.0000000000000004*x28;
    double x38 = 3.0000000000000004*x36/((x4)*(x4));

    result += x35*(-966175309088.08997*x11*x18 + x11*x3 - 322058436362.69666*x12*x18 + x12*x3 + 148470.48537177511*x13*x31 - 644116872725.39331*x18*exp(-2550.9849899999999*x6)/((x9)*(x9)*(x9)) - 148470.48537177511*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(186.36863737044234*x2 - 16.680099032862614)/pow(x4, 5.0/2.0) - x22*x32 + 2169174.8064017668*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(124.24575824696156*x2 - 11.120066021908409) - 4338349.6128035337*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -16.680099032862614*x2 + 93.184318685221172*x3 + 0.54777557845000002;
    double x5 = sqrt(x4);
    double x6 = 2.8344277666666668*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 62.122879123480779*x2 - 5.5600330109542044;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 5.6688555333333337*x5;
    double x16 = exp(-x15)/((x8)*(x8));
    double x17 = 1262486.59880471*x16;
    double x18 = x7/x8;
    double x19 = 1262486.59880471*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(483.17794873818377*x2 - 24.711257826463129);
    double x22 = x20*x21;
    double x23 = 445411.45611532545*x18;
    double x24 = 1.0/T;
    double x25 = x24*x5;
    double x26 = 850.32833000000005*x25;
    double x27 = exp(-x26);
    double x28 = -x27 + 1;
    double x29 = x27/x28;
    double x30 = 445411.45611532533*x29;
    double x31 = pow(T, -2);
    double x32 = x12*x31;
    double x33 = 1700.6566600000001*x25;
    double x34 = exp(-x33)/((x28)*(x28));
    double x35 = 378745979.64141291*x24;
    double x36 = x14*x35;
    double x37 = x13*(124.24575824696156*x2 - 11.120066021908409);
    double x38 = ((x10)*(x10));
    double x39 = x0*x38;
    double x40 = x37*x39;
    double x41 = (186.36863737044234*x2 - 16.680099032862614)/pow(x4, 5.0/2.0);
    double x42 = x39*x41;
    double x43 = 1.0/x4;
    double x44 = 144.95338462145514*x2 - 9.2667216849236738;
    double x45 = x43*x44;
    double x46 = x10*x3;
    double x47 = x17*x46;
    double x48 = x43*(289.90676924291029*x2 - 18.533443369847348);
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
    double x61 = -3*x58 + 1.0584146949449513*x60;
    double x62 = 148470.48537177514*x20;
    double x63 = exp(x26);
    double x64 = x63 - 1;
    double x65 = 1.0/x64;
    double x66 = Debye(x26);
    double x67 = T*x20*x66;
    double x68 = -3*x65 + 0.0035280489831498382*x67;
    double x69 = 148470.48537177511*x20;
    double x70 = 1.0584146949449513*x59;
    double x71 = x70*x9;
    double x72 = x56/((x57)*(x57));
    double x73 = 8.5032832999999997*x72;
    double x74 = x20*x73;
    double x75 = x43*(-9.0*x58 + 3.175244084834854*x60) + x71 - x74;
    double x76 = 296940.97074355028*x75;
    double x77 = x20*x49;
    double x78 = 0.0035280489831498382*T*x66;
    double x79 = x78*x9;
    double x80 = x63/((x64)*(x64));
    double x81 = 2550.9849899999999*x24*x80;
    double x82 = x20*x81;
    double x83 = x43*(-9.0000000000000018*x65 + 0.010584146949449516*x67) + x79 - x82;
    double x84 = 296940.97074355022*x83;
    double x85 = x10*x2;
    double x86 = x41*x85;
    double x87 = x2*x38;
    double x88 = x43*x87;
    double x89 = x87*x9;
    double x90 = 3.0*x61;
    double x91 = x13*x87;
    double x92 = x37*x85;
    double x93 = x31*x88;
    double x94 = 3.0000000000000004*x68;

    result += (-11497144463.017965*x0 + 10735281.212089892*x12*x16 + 3578427.0706966305*x12*x18 + x12*x76 - x12*x84 + 7156854.141393261*x12*exp(-8.5032832999999997*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x19 + x17*x40 - 1336234.3683459763*x18*x50 + x19*x40 + 557880136.93904233*x2 + x21*x61*x62 - x21*x68*x69 + x22*x23 - x22*x30 + x23*x42 - 322058436362.69666*x29*x32 - x29*x36 + 1336234.3683459759*x29*x50 - x29*x52 - 3299160786.9596205*x3 - x30*x42 - 966175309088.08997*x32*x34 - x34*x36 - x34*x52 + 148470.48537177514*x42*x61 - 148470.48537177511*x42*x68 - x45*x47 - x45*x51 + x45*x54 + x45*x55 + x46*x62*(-x44*x71 + x44*x74 - x45*x90 + x70*x86 + 24.101942293352963*x72*x88 - x73*x89 + 3.0*x75*x88 + x90*x91 + x90*x92 - 48.203884586705925*x88*exp(x15)/((x57)*(x57)*(x57))) - x46*x69*(-x44*x79 + x44*x82 - x45*x94 + x78*x86 + 2169174.8064017668*x80*x93 - x81*x89 + 3.0000000000000004*x83*x88 + x91*x94 + x92*x94 - 4338349.6128035337*x93*exp(x33)/((x64)*(x64)*(x64))) - x47*x48 - x48*x51 + x48*x54 + x48*x55 - 445411.45611532545*x50*x61 + 445411.45611532533*x50*x68 - x76*x77 + x77*x84 - 644116872725.39331*x32*exp(-2550.9849899999999*x25)/((x28)*(x28)*(x28)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*10.9401;
static const double Vmax = 1.15*10.9401;
static double V = 0.9*10.9401;

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



const char *NaNAL_stx21_em_coder_calib_identifier(void) {
    return identifier;
}

const char *NaNAL_stx21_em_coder_calib_name(void) {
    return "NaNAL_stx21_em";
}

const char *NaNAL_stx21_em_coder_calib_formula(void) {
    return "Na3Al3Si3O12";
}

const double NaNAL_stx21_em_coder_calib_mw(void) {
    return 426.16323;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,12.0,0.0,0.0,3.0,
        0.0,3.0,3.0,0.0,0.0,0.0,
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

const double *NaNAL_stx21_em_coder_calib_elements(void) {
    return elmformula;
}

double NaNAL_stx21_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double NaNAL_stx21_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double NaNAL_stx21_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double NaNAL_stx21_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double NaNAL_stx21_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double NaNAL_stx21_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double NaNAL_stx21_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double NaNAL_stx21_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double NaNAL_stx21_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double NaNAL_stx21_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double NaNAL_stx21_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double NaNAL_stx21_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double NaNAL_stx21_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double NaNAL_stx21_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double NaNAL_stx21_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double NaNAL_stx21_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double NaNAL_stx21_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double NaNAL_stx21_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double NaNAL_stx21_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int NaNAL_stx21_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **NaNAL_stx21_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **NaNAL_stx21_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void NaNAL_stx21_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int NaNAL_stx21_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double NaNAL_stx21_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int NaNAL_stx21_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double NaNAL_stx21_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double NaNAL_stx21_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

