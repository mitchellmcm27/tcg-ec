
static char *identifier = "Quartz_stx21_em.emml:58925e2a960ec0569cd08d5fb7f5d6fbdee0987c:Sat Mar  9 00:08:57 2024";



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
    double x3 = sqrt(-0.73607897226924934*x1 + 0.45624634927040242*x2 + 1.2742145938);
    double x4 = 2.9473493666666668*x3;
    double x5 = 884.20480999999995*x3/T;

    result += 74.830163563379159*T*log(1 - exp(-x5)) - 24.943387854459719*T*Debye(x5) + 57515842.554625735*x1 - 103074191.02015355*x2 - 22449.049069013745*log(1 - exp(-x4)) + 7483.0163563379156*Debye(x4) - 11535365.979120873 + 61453308.643032715/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(0.45624634927040242*pow(x0, 4.0/3.0) - 0.73607897226924934*pow(x0, 2.0/3.0) + 1.2742145938);
    double x2 = x1/T;
    double x3 = 884.20480999999995*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -66165.190555826586*x2*x5/x6 + 22055.063518608862*x2*(-0.0033928790774164644*T*x4/x1 + 3/(exp(x3) - 1)) - 24.943387854459719*x4 + 74.830163563379159*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-0.73607897226924934*x1 + 0.45624634927040242*x3 + 1.2742145938);
    double x6 = 2.9473493666666668*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(0.24535965742308311*x2 - 0.30416423284693495*x4);
    double x10 = 66165.190555826586*x9;
    double x11 = 884.20480999999995*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) - 38343895.036417156*x2 + 137432254.69353807*x4 + 22055.063518608866*x9*(-1.0178637232249392*x8*Debye(x6) + 3/(exp(x6) - 1)) - 22055.063518608862*x9*(-0.0033928790774164644*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 122906617.28606543/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-0.73607897226924934*x2 + 0.45624634927040242*x3 + 1.2742145938);
    double x5 = x0*x4;
    double x6 = 884.20480999999995*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-43063254.852056526*x2 + 26692044.677462842*x3 + 74546115.099383116);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1768.4096199999999*x5)/((x8)*(x8)) + 22055.063518608862*x4*(x0*(0.010178637232249394*T*x11 - 9.0/x13) + 0.0033928790774164644*x11 - 2652.6144299999996*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 0.45624634927040242*pow(x1, 4.0/3.0) - 0.73607897226924934*x2 + 1.2742145938;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 884.20480999999995*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 58503579.744028442*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(0.30416423284693495*x2 - 0.24535965742308311)*(22055.063518608862*x4*(-2652.6144299999996*x0*x11*x12/((x13)*(x13)) + 0.0033928790774164644*x10/pow(x3, 3.0/2.0) + (0.010178637232249394*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1768.4096199999999*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -0.73607897226924934*x2 + 0.45624634927040242*x3 + 1.2742145938;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x2*x6*(0.70971654330951484*x2 - 0.4089327623718052);
    double x8 = 66165.190555826586*x7;
    double x9 = 2.9473493666666668*x5;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = x10/x11;
    double x13 = 1.0/x4;
    double x14 = x3*((0.30416423284693495*x2 - 0.24535965742308311)*(0.30416423284693495*x2 - 0.24535965742308311));
    double x15 = x13*x14;
    double x16 = 195011.93248009481*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 66165.190555826586*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 884.20480999999995*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 58503579.744028442*x15*x20;
    double x27 = exp(x9);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x9);
    double x31 = x30*x6;
    double x32 = -66165.190555826601*x29 + 22449.049069013749*x31;
    double x33 = exp(x22);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x22);
    double x37 = T*x36*x6;
    double x38 = -66165.190555826586*x35 + 74.830163563379159*x37;
    double x39 = x14*x6;

    result += x0*(368719851.85819626*x0 + x12*x16 + x12*x19 - x12*x8 + x18*x32 - x18*x38 - x19*x25 + 63906491.727361925*x2 - x25*x26 + x25*x8 - 320675260.95158881*x3 - x32*x7 + x38*x7 + 22055.063518608866*x39*(x13*(-9.0*x29 + 3.0535911696748177*x31) + 1.0178637232249392*x17*x30 - 8.8420480999999995*x27*x6/((x28)*(x28))) - 22055.063518608862*x39*(0.0033928790774164644*T*x17*x36 + x13*(-9.0*x35 + 0.010178637232249394*x37) - 2652.6144299999996*x20*x33*x6/((x34)*(x34))) - x26*exp(-1768.4096199999999*x21)/((x24)*(x24)) + x16*exp(-5.8946987333333336*x5)/((x11)*(x11)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-129189764.55616958*x2 + 80076134.032388523*x3 + 223638345.29814932);
    double x5 = 1.0/T;
    double x6 = -0.73607897226924934*x2 + 0.45624634927040242*x3 + 1.2742145938;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 884.20480999999995*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1768.4096199999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = -2652.6144299999996*x0*x20*x7;
    double x22 = 0.010178637232249394*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 22055.063518608862*x7;
    double x25 = x15*x6;

    result += x0*(-155187439835.66553*x13*x16 + x13*x4 - 51729146611.888512*x14*x16 + x14*x4 + x24*(0.0033928790774164644*x17 + x21 + x23) - x24*(2345454.4380814079*x20*x25 + x21 + x22 + 3.0*x23 - 4690908.8761628158*x25*exp(x12)/((x19)*(x19)*(x19))) - 103458293223.77702*x16*exp(-2652.6144299999996*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(35589392.90328379*x2 - 28708836.568037685);
    double x5 = 0.45624634927040242*pow(x1, 4.0/3.0) - 0.73607897226924934*x2 + 1.2742145938;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 884.20480999999995*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1768.4096199999999*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 0.30416423284693495*x2 - 0.24535965742308311;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0033928790774164644*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2652.6144299999996*x22*x3;
    double x24 = 0.010178637232249394*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 22055.063518608862*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-155187439835.66553*x12*x16 + x12*x4 - 51729146611.888512*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7957.8432899999989*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 2345454.4380814079*x15*x22 - 4690908.8761628158*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 103458293223.77702*x16*exp(-2652.6144299999996*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(41520958.387164414*x2 - 23924030.473364741);
    double x4 = 0.45624634927040242*pow(x1, 4.0/3.0) - 0.73607897226924934*x2 + 1.2742145938;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 884.20480999999995*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1768.4096199999999*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 0.30416423284693495*x2 - 0.24535965742308311;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0033928790774164644*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2652.6144299999996*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.010178637232249394*x30));
    double x32 = 0.70971654330951484*x2 - 0.4089327623718052;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0033928790774164644*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-155187439835.66553*x11*x18 + x11*x3 - 51729146611.888512*x12*x18 + x12*x3 + 22055.063518608862*x13*x31 - 103458293223.77702*x18*exp(-2652.6144299999996*x6)/((x9)*(x9)*(x9)) - 22055.063518608862*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(0.91249269854080484*x2 - 0.73607897226924934)/pow(x4, 5.0/2.0) - x22*x32 + 2345454.4380814079*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(0.60832846569386989*x2 - 0.49071931484616621) - 4690908.8761628158*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -0.73607897226924934*x2 + 0.45624634927040242*x3 + 1.2742145938;
    double x5 = sqrt(x4);
    double x6 = 2.9473493666666668*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 0.30416423284693495*x2 - 0.24535965742308311;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = 5.8946987333333336*x5;
    double x14 = exp(-x13)/((x8)*(x8));
    double x15 = pow(x4, -2);
    double x16 = 195011.93248009481*x15;
    double x17 = x11*x16;
    double x18 = x7/x8;
    double x19 = 66165.190555826586*x18;
    double x20 = 1.0/x5;
    double x21 = 2.3657218110317162*x2 - 1.0904873663248138;
    double x22 = x2*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 884.20480999999995*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 66165.190555826586*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1768.4096199999999*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = x11*x15;
    double x35 = 58503579.744028442*x23;
    double x36 = x33*x35;
    double x37 = x28*x35;
    double x38 = 0.60832846569386989*x2 - 0.49071931484616621;
    double x39 = ((x10)*(x10));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = x16*x41;
    double x43 = (0.91249269854080484*x2 - 0.73607897226924934)/pow(x4, 5.0/2.0);
    double x44 = x40*x43;
    double x45 = 0.70971654330951484*x2 - 0.4089327623718052;
    double x46 = x10*x3;
    double x47 = x45*x46;
    double x48 = 1.0/x4;
    double x49 = 195011.93248009481*x48;
    double x50 = x46*(1.4194330866190297*x2 - 0.81786552474361041);
    double x51 = x49*x50;
    double x52 = 198495.57166747976*x9;
    double x53 = x18*x47;
    double x54 = x28*x47;
    double x55 = x15*x41;
    double x56 = x48*x50;
    double x57 = exp(x6);
    double x58 = x57 - 1;
    double x59 = 1.0/x58;
    double x60 = Debye(x6);
    double x61 = x20*x60;
    double x62 = -3*x59 + 1.0178637232249392*x61;
    double x63 = x2*x62;
    double x64 = 22055.063518608866*x20;
    double x65 = 22055.063518608862*x20;
    double x66 = exp(x25);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x25);
    double x70 = T*x20*x69;
    double x71 = -3*x68 + 0.0033928790774164644*x70;
    double x72 = x2*x71;
    double x73 = x45*x62;
    double x74 = x46*x9;
    double x75 = x45*x71;
    double x76 = 1.0178637232249392*x60;
    double x77 = x76*x9;
    double x78 = x57/((x58)*(x58));
    double x79 = 8.8420480999999995*x78;
    double x80 = x20*x79;
    double x81 = x48*(-9.0*x59 + 3.0535911696748177*x61) + x77 - x80;
    double x82 = 44110.127037217731*x81;
    double x83 = x20*x47;
    double x84 = 0.0033928790774164644*T*x69;
    double x85 = x84*x9;
    double x86 = x66/((x67)*(x67));
    double x87 = 2652.6144299999996*x23*x86;
    double x88 = x20*x87;
    double x89 = x48*(-9.0*x68 + 0.010178637232249394*x70) + x85 - x88;
    double x90 = 44110.127037217724*x89;
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

    result += (-1474879407.432785*x0 + 1724304.8870629505*x12*x14 + 574768.29568765021*x12*x18 + x12*x82 - x12*x90 + 1149536.5913753004*x12*exp(-8.8420480999999995*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x42 - x14*x47*x49 - x14*x51 + x17*x18 + x18*x42 - x18*x51 + x19*x22 + x19*x44 - 170417311.27296513*x2 + x21*x63*x64 - x21*x65*x72 - x22*x29 - 51729146611.888512*x28*x31 - x29*x44 + 1068917536.505296*x3 - 155187439835.66553*x31*x33 - x34*x36 - x34*x37 + x35*x48*x54 + x36*x47*x48 - x36*x55 + x36*x56 - x37*x55 + x37*x56 + 22055.063518608866*x44*x62 - 22055.063518608862*x44*x71 + x46*x64*(x39*x97 - x45*x77 + x45*x80 - x73*x95 + x76*x91 + 26.060604867571204*x78*x93 - x79*x94 + x81*x99 + x97*x98 - 52.121209735142408*x93*exp(x13)/((x58)*(x58)*(x58))) - x46*x65*(2345454.4380814079*x100*x86 - 4690908.8761628158*x100*exp(x32)/((x67)*(x67)*(x67)) - x45*x85 + x45*x88 + x71*x92*x96 + x72*x96*x98 - x75*x95 + x84*x91 - x87*x94 + x89*x99) - x49*x53 - x52*x53 + x52*x54 - 66165.190555826601*x73*x74 + 66165.190555826586*x74*x75 - x82*x83 + x83*x90 - 103458293223.77702*x31*exp(-2652.6144299999996*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.2421;
static const double Vmax = 1.15*2.2421;
static double V = 0.9*2.2421;

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



const char *Quartz_stx21_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Quartz_stx21_em_coder_calib_name(void) {
    return "Quartz_stx21_em";
}

const char *Quartz_stx21_em_coder_calib_formula(void) {
    return "SiO2";
}

const double Quartz_stx21_em_coder_calib_mw(void) {
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

const double *Quartz_stx21_em_coder_calib_elements(void) {
    return elmformula;
}

double Quartz_stx21_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Quartz_stx21_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Quartz_stx21_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Quartz_stx21_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Quartz_stx21_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Quartz_stx21_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Quartz_stx21_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Quartz_stx21_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Quartz_stx21_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Quartz_stx21_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Quartz_stx21_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Quartz_stx21_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Quartz_stx21_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Quartz_stx21_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Quartz_stx21_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Quartz_stx21_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Quartz_stx21_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Quartz_stx21_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Quartz_stx21_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Quartz_stx21_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Quartz_stx21_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Quartz_stx21_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Quartz_stx21_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Quartz_stx21_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Quartz_stx21_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Quartz_stx21_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Quartz_stx21_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Quartz_stx21_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Quartz_stx21_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Quartz_stx21_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

