
static char *identifier = "Diopside_slb_em.emml:e57d33c1e586a6de659f2f3e9d23bf5a05e72051:Fri May 26 02:48:56 2023";



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
    double x3 = sqrt(-1.3374580978259498*x1 + 21.109873037677488*x2 - 0.32383121749999999);
    double x4 = 2.6085768666666671*x3;
    double x5 = 782.57306000000005*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) - 23527186.23224432*x1 - 21901174.174610846*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) + 2015612.3198061325 + 148560076.5378395/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(21.109873037677488*pow(x0, 4.0/3.0) - 1.3374580978259498*pow(x0, 2.0/3.0) - 0.32383121749999999);
    double x2 = x1/T;
    double x3 = 782.57306000000005*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -195200.2336003138*x2*x5/x6 + 65066.74453343793*x2*(-0.0038335078899853769*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-1.3374580978259498*x1 + 21.109873037677488*x3 - 0.32383121749999999);
    double x6 = 2.6085768666666671*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(0.44581936594198324*x2 - 14.073248691784991*x4);
    double x10 = 195200.2336003138*x9;
    double x11 = 782.57306000000005*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 15684790.821496213*x2 + 29201565.566147793*x4 + 65066.744533437937*x9*(-1.1500523669956129*x8*Debye(x6) + 3/(exp(x6) - 1)) - 65066.74453343793*x9*(-0.0038335078899853769*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 297120153.075679/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-1.3374580978259498*x2 + 21.109873037677488*x3 - 0.32383121749999999);
    double x5 = x0*x4;
    double x6 = 782.57306000000005*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(204308018.10134214*x2 - 3224711360.8340559*x3 + 49467952.943210311);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1565.1461200000001*x5)/((x8)*(x8)) - 65066.74453343793*x4*(x0*(0.011500523669956131*T*x11 - 9.0/x13) + 0.0038335078899853769*x11 - 2347.7191800000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 1.3374580978259498*x2;
    double x4 = 21.109873037677488*pow(x1, 4.0/3.0);
    double x5 = -x3 + x4 - 0.32383121749999999;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 782.57306000000005*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 152758444.12131241*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(14.073248691784991*x2 - 0.44581936594198324)*(65066.74453343793*x6*(2347.7191800000001*x0*x13*x14/((x15)*(x15)) - 0.0038335078899853769*x12/pow(x5, 3.0/2.0) + (0.011500523669956131*x12*x13 - 9.0/x15)/(x3 - x4 + 0.32383121749999999)) - x11*x9/x10 - x11*exp(-1565.1461200000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 1.3374580978259498*x2;
    double x5 = 21.109873037677488*x3;
    double x6 = -x4 + x5 - 0.32383121749999999;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*(32.837580280831645*x2 - 0.74303227656997206);
    double x10 = x8*x9;
    double x11 = 195200.2336003138*x10;
    double x12 = 2.6085768666666671*x7;
    double x13 = exp(-x12);
    double x14 = -x13 + 1;
    double x15 = x13/x14;
    double x16 = 1.0/(x4 - x5 + 0.32383121749999999);
    double x17 = x3*((14.073248691784991*x2 - 0.44581936594198324)*(14.073248691784991*x2 - 0.44581936594198324));
    double x18 = x16*x17;
    double x19 = 509194.81373770803*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 195200.2336003138*x21;
    double x23 = 1.0/T;
    double x24 = x23*x7;
    double x25 = 782.57306000000005*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 152758444.12131241*x18*x23;
    double x30 = exp(x12);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x12);
    double x34 = x33*x8;
    double x35 = -195200.2336003138*x32 + 74830.163563379159*x34;
    double x36 = exp(x25);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = Debye(x25);
    double x40 = T*x39*x8;
    double x41 = -3*x38 + 0.0038335078899853769*x40;
    double x42 = 65066.74453343793*x8;

    result += x0*(891360459.22703695*x0 - x10*x35 - x11*x15 + x11*x28 - x15*x19 + x15*x22 + x17*x42*(-0.0038335078899853769*T*x20*x39 + x16*(-9.0*x38 + 0.011500523669956131*x40) + 2347.7191800000001*x23*x36*x8/((x37)*(x37))) - 65066.744533437937*x17*x8*(x16*(-8.9999999999999982*x32 + 3.4501571009868379*x34) - 1.1500523669956129*x20*x33 + 7.8257306000000018*x30*x8/((x31)*(x31))) - 26141318.035827018*x2 + x21*x35 - 65066.74453343793*x21*x41 - x22*x28 + x28*x29 - 68136986.321011513*x3 + x41*x42*x9 + x29*exp(-1565.1461200000001*x24)/((x27)*(x27)) - x19*exp(-5.2171537333333342*x7)/((x14)*(x14)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(612924054.30402648*x2 - 9674134082.5021687*x3 + 148403858.82963094);
    double x5 = 1.0/T;
    double x6 = 1.3374580978259498*x2;
    double x7 = 21.109873037677488*x3;
    double x8 = -x6 + x7 - 0.32383121749999999;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 782.57306000000005*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1565.1461200000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2347.7191800000001*x0*x22*x9;
    double x24 = 0.011500523669956131*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 65066.74453343793*x9;
    double x27 = x17*(x6 - x7 + 0.32383121749999999);

    result += x0*(-358633929170.56342*x15*x18 - x15*x4 - 119544643056.85448*x16*x18 - x16*x4 + x26*(0.0038335078899853769*x19 + x23 + x25) - x26*(-1837261.782713291*x22*x27 + x23 + x24 + 3.0*x25 + 3674523.565426582*x27*exp(x14)/((x21)*(x21)*(x21))) - 239089286113.70895*x18*exp(-2347.7191800000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(4299615147.7787409*x2 - 136205345.40089476);
    double x5 = 1.3374580978259498*x2;
    double x6 = 21.109873037677488*pow(x1, 4.0/3.0);
    double x7 = -x5 + x6 - 0.32383121749999999;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 782.57306000000005*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1565.1461200000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = 14.073248691784991*x2 - 0.44581936594198324;
    double x17 = pow(T, -3);
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0038335078899853769*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2347.7191800000001*x24*x3;
    double x26 = 0.011500523669956131*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 65066.74453343793*x16;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = 1.0/(x5 - x6 + 0.32383121749999999);

    result += x0*x1*x2*(-358633929170.56342*x14*x18 + x14*x4 - 119544643056.85448*x15*x18 + x15*x4 + x19*x29*(x19*x21 - x25*x8 + x28) - x29*x8*(x0*(-7043.1575400000002*x0*x19*x24 + x26*x30 - 3.0*x27*x31) + 1837261.782713291*x17*x24 - 3674523.565426582*x17*exp(x13)/((x23)*(x23)*(x23)) + x19*x25 + x21*x30 - x28*x31) - 239089286113.70895*x18*exp(-2347.7191800000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(5016217672.4085312*x2 - 113504454.50074562);
    double x4 = 1.3374580978259498*x2;
    double x5 = 21.109873037677488*pow(x1, 4.0/3.0);
    double x6 = -x4 + x5 - 0.32383121749999999;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 782.57306000000005*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1565.1461200000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -2);
    double x16 = 1.0/x7;
    double x17 = 14.073248691784991*x2 - 0.44581936594198324;
    double x18 = ((x17)*(x17));
    double x19 = x18*x2;
    double x20 = x16*x19;
    double x21 = x15*x20;
    double x22 = pow(x6, -3.0/2.0);
    double x23 = Debye(x9);
    double x24 = 0.0038335078899853769*T*x23;
    double x25 = x22*x24;
    double x26 = exp(x9);
    double x27 = x26 - 1;
    double x28 = x26/((x27)*(x27));
    double x29 = 2347.7191800000001*x28;
    double x30 = x0*x16*x29;
    double x31 = x4 - x5 + 0.32383121749999999;
    double x32 = 1.0/x31;
    double x33 = 1.0/x27;
    double x34 = T*x16*x23;
    double x35 = x32*(-9.0*x33 + 0.011500523669956131*x34);
    double x36 = 32.837580280831645*x2 - 0.74303227656997206;
    double x37 = x17*x2;
    double x38 = x19*x32;
    double x39 = x15*x38;
    double x40 = x0*x2;
    double x41 = -9.0*x33 + 0.011500523669956131*x34;
    double x42 = x41/((x31)*(x31));

    result += x40*(-358633929170.56342*x13*x21 + x13*x3 - 119544643056.85448*x14*x21 + x14*x3 - 65066.74453343793*x20*(-x25 + x30 + x35) - 65066.74453343793*x7*(-x18*x22*x29*x40 + x19*x42 + x24*x37*(42.219746075354976*x2 - 1.3374580978259498)/pow(x6, 5.0/2.0) - x25*x36 - 1837261.782713291*x28*x39 + x30*x36 + x32*x36*x41 + x37*x42*(28.146497383569983*x2 - 0.89163873188396647) - 3.0*x38*(x25 - x30 - x35) + 3674523.565426582*x39*exp(x12)/((x27)*(x27)*(x27))) - 239089286113.70895*x21*exp(-2347.7191800000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 1.3374580978259498*x2;
    double x5 = 21.109873037677488*x3;
    double x6 = -x4 + x5 - 0.32383121749999999;
    double x7 = sqrt(x6);
    double x8 = 2.6085768666666671*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 14.073248691784991*x2 - 0.44581936594198324;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.2171537333333342*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = x4 - x5 + 0.32383121749999999;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 509194.81373770803*x19;
    double x21 = x9/x10;
    double x22 = 195200.2336003138*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(109.45860093610548*x2 - 1.981419404186592);
    double x25 = x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 782.57306000000005*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 195200.2336003138*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1565.1461200000001*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 152758444.12131241*x26;
    double x38 = x19*x37;
    double x39 = x18*(28.146497383569983*x2 - 0.89163873188396647);
    double x40 = ((x12)*(x12));
    double x41 = x0*x40;
    double x42 = x39*x41;
    double x43 = 509194.81373770803*x16;
    double x44 = 509194.81373770803*x21;
    double x45 = (42.219746075354976*x2 - 1.3374580978259498)/pow(x6, 5.0/2.0);
    double x46 = x41*x45;
    double x47 = 1.0/x17;
    double x48 = 32.837580280831645*x2 - 0.74303227656997206;
    double x49 = x47*x48;
    double x50 = x12*x3;
    double x51 = x43*x50;
    double x52 = x47*(65.67516056166329*x2 - 1.4860645531399441);
    double x53 = x44*x50;
    double x54 = x48*x50;
    double x55 = x11*x54;
    double x56 = 585600.70080094144*x55;
    double x57 = x37*x42;
    double x58 = x37*x50;
    double x59 = x36*x58;
    double x60 = x31*x58;
    double x61 = exp(x8);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x8);
    double x65 = x23*x64;
    double x66 = -3*x63 + 1.1500523669956129*x65;
    double x67 = 65066.744533437937*x23;
    double x68 = exp(x28);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = T*Debye(x28);
    double x72 = x23*x71;
    double x73 = -3*x70 + 0.0038335078899853769*x72;
    double x74 = 65066.74453343793*x23;
    double x75 = 1.1500523669956129*x64;
    double x76 = x11*x75;
    double x77 = x61/((x62)*(x62));
    double x78 = 7.8257306000000018*x77;
    double x79 = x23*x78;
    double x80 = x47*(-8.9999999999999982*x63 + 3.4501571009868379*x65) - x76 + x79;
    double x81 = 130133.48906687587*x80;
    double x82 = x23*x54;
    double x83 = 0.0038335078899853769*x71;
    double x84 = x11*x83;
    double x85 = x68/((x69)*(x69));
    double x86 = 2347.7191800000001*x26*x85;
    double x87 = x23*x86;
    double x88 = x47*(-9.0*x70 + 0.011500523669956131*x72) - x84 + x87;
    double x89 = 130133.48906687586*x88;
    double x90 = x12*x2;
    double x91 = x45*x90;
    double x92 = x2*x40;
    double x93 = x47*x92;
    double x94 = x11*x92;
    double x95 = 2.9999999999999996*x66;
    double x96 = x18*x92;
    double x97 = x39*x90;
    double x98 = x33*x93;
    double x99 = 3.0*x73;

    result += (-3565441836.9081478*x0 + 3984821.4352284828*x14*x16 + 1328273.8117428275*x14*x21 - x14*x81 + x14*x89 + x16*x20 + 69710181.428872049*x2 + x20*x21 - x21*x56 + x22*x25 + x22*x46 + x24*x66*x67 - x24*x73*x74 - x25*x32 + 227123287.73670503*x3 - 119544643056.85448*x31*x34 - x31*x38 + x31*x56 - x31*x57 - x32*x46 - 358633929170.56342*x34*x36 - x36*x38 - x36*x57 + x42*x43 + x42*x44 + 65066.744533437937*x46*x66 - 65066.74453343793*x46*x73 + x49*x51 + x49*x53 - x49*x59 - x49*x60 + x50*x67*(-x48*x76 + x48*x79 + x49*x95 + x75*x91 - 20.414019807925463*x77*x93 - x78*x94 + 2.9999999999999996*x80*x93 + x95*x96 + x95*x97 + 40.828039615850926*x93*exp(x15)/((x62)*(x62)*(x62))) - x50*x74*(-x48*x84 + x48*x87 + x49*x99 + x83*x91 - 1837261.782713291*x85*x98 - x86*x94 + 3.0*x88*x93 + x96*x99 + x97*x99 + 3674523.565426582*x98*exp(x35)/((x69)*(x69)*(x69))) + x51*x52 + x52*x53 - x52*x59 - x52*x60 - 195200.2336003138*x55*x66 + 195200.23360031378*x55*x73 + x81*x82 - x82*x89 - 239089286113.70895*x34*exp(-2347.7191800000001*x27)/((x30)*(x30)*(x30)) + 2656547.6234856551*x14*exp(-7.8257306000000018*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.6039;
static const double Vmax = 1.15*6.6039;
static double V = 0.9*6.6039;

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



const char *Diopside_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Diopside_slb_em_coder_calib_name(void) {
    return "Diopside_slb_em";
}

const char *Diopside_slb_em_coder_calib_formula(void) {
    return "CaMgSi2O6";
}

const double Diopside_slb_em_coder_calib_mw(void) {
    return 216.55239999999998;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,6.0,0.0,0.0,0.0,
        1.0,0.0,2.0,0.0,0.0,0.0,
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

const double *Diopside_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Diopside_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Diopside_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Diopside_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Diopside_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Diopside_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Diopside_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Diopside_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Diopside_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Diopside_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Diopside_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Diopside_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Diopside_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Diopside_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Diopside_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Diopside_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Diopside_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Diopside_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Diopside_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Diopside_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Diopside_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Diopside_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Diopside_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Diopside_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Diopside_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Diopside_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Diopside_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Diopside_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Diopside_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Diopside_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Diopside_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Diopside_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Diopside_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Diopside_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Diopside_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Diopside_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Diopside_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

