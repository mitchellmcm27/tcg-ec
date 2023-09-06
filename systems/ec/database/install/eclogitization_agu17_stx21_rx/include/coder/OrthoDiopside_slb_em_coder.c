
static char *identifier = "OrthoDiopside_slb_em.emml:cdbecfc30ddad5508c6656fb89e3f2f4e4cc4f44:Fri May 26 02:50:43 2023";



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
    double x3 = sqrt(40.611377456322202*x1 - 57.739979479979205*x2 - 5.831495441225);
    double x4 = 2.4818297666666669*x3;
    double x5 = 744.54893000000004*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) + 74814767.359745055*x1 - 374396321.41549832*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) - 7227436.4221247602 + 574730847.58657849/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-57.739979479979205*pow(x0, 4.0/3.0) + 40.611377456322202*pow(x0, 2.0/3.0) - 5.831495441225);
    double x2 = x1/T;
    double x3 = 744.54893000000004*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -185715.72737612983*x2*x5/x6 + 61905.242458709938*x2*(-0.0040292852210532352*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(40.611377456322202*x1 - 57.739979479979205*x3 - 5.831495441225);
    double x6 = 2.4818297666666669*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-13.537125818774067*x2 + 38.493319653319467*x4);
    double x10 = 185715.72737612983*x9;
    double x11 = 744.54893000000004*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) - 49876511.573163368*x2 + 499195095.22066438*x4 + 61905.242458709945*x9*(-1.2087855663159706*x8*Debye(x6) + 3/(exp(x6) - 1)) - 61905.242458709938*x9*(-0.0040292852210532352*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 1149461695.173157/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(40.611377456322202*x2 - 57.739979479979205*x3 - 5.831495441225);
    double x5 = x0*x4;
    double x6 = 744.54893000000004*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-5615515723.2150116*x2 + 7983963680.5389643*x3 + 806346802.08212829);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1489.0978600000001*x5)/((x8)*(x8)) - 61905.242458709938*x4*(x0*(0.012087855663159706*T*x11 - 9.0/x13) + 0.0040292852210532352*x11 - 2233.6467900000002*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 40.611377456322202*x2;
    double x4 = 57.739979479979205*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 5.831495441225;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 744.54893000000004*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 138274446.10206917*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(38.493319653319467*x2 - 13.537125818774067)*(-61905.242458709938*x6*(2233.6467900000002*x0*x13*x14/((x15)*(x15)) - 0.0040292852210532352*x12/pow(x5, 3.0/2.0) + (0.012087855663159706*x12*x13 - 9.0/x15)/(-x3 + x4 + 5.831495441225)) + x11*x9/x10 + x11*exp(-1489.0978600000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 40.611377456322202*x2;
    double x5 = 57.739979479979205*x3;
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = 2.4818297666666669*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x2*(89.817745857745422*x2 - 22.561876364623444);
    double x14 = x12*x13;
    double x15 = 185715.72737612983*x14;
    double x16 = 1.0/(-x4 + x5 + 5.831495441225);
    double x17 = x3*((38.493319653319467*x2 - 13.537125818774067)*(38.493319653319467*x2 - 13.537125818774067));
    double x18 = x16*x17;
    double x19 = 460914.8203402306*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 185715.72737612983*x21;
    double x23 = 1.0/T;
    double x24 = x23*x7;
    double x25 = 744.54893000000004*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 138274446.10206917*x18*x23;
    double x30 = exp(x8);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x8);
    double x34 = x12*x33;
    double x35 = -185715.72737612983*x32 + 74830.163563379159*x34;
    double x36 = exp(x25);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = T*Debye(x25);
    double x40 = x12*x39;
    double x41 = -3*x38 + 0.0040292852210532352*x40;
    double x42 = 61905.242458709938*x12;

    result += x0*(3448385085.5194712*x0 + x11*x15 - x11*x19 + x11*x22 - 61905.242458709945*x12*x17*(7.4454893000000002*x12*x30/((x31)*(x31)) + x16*(-9.0*x32 + 3.6263566989479115*x34) - 1.2087855663159706*x20*x33) - x13*x41*x42 + x14*x35 - x15*x28 + x17*x42*(2233.6467900000002*x12*x23*x36/((x37)*(x37)) + x16*(-9.0*x38 + 0.012087855663159706*x40) - 0.0040292852210532352*x20*x39) + 83127519.288605615*x2 + x21*x35 - 61905.242458709938*x21*x41 - x22*x28 + x28*x29 - 1164788555.5148835*x3 + x29*exp(-1489.0978600000001*x24)/((x27)*(x27)) - x19*exp(-4.9636595333333338*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-16846547169.645035*x2 + 23951891041.616894*x3 + 2419040406.2463851);
    double x5 = 1.0/T;
    double x6 = 40.611377456322202*x2;
    double x7 = 57.739979479979205*x3;
    double x8 = x6 - x7 - 5.831495441225;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 744.54893000000004*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1489.0978600000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.012087855663159706*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 61905.242458709938*x9;
    double x27 = x17*(-x6 + x7 + 5.831495441225);

    result += x0*(-308856272674.91479*x15*x18 - x15*x4 - 102952090891.63828*x16*x18 - x16*x4 + x26*(0.0040292852210532352*x19 - 2233.6467900000002*x23 + x25) - x26*(-1663059.327492435*x22*x27 - 2233.6467899999998*x23 + x24 + 3.0*x25 + 3326118.65498487*x27*exp(x14)/((x21)*(x21)*(x21))) - 205904181783.27655*x18*exp(-2233.6467900000002*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(10645284907.385284*x2 - 3743677148.8100076);
    double x5 = 40.611377456322202*x2;
    double x6 = 57.739979479979205*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 5.831495441225;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 744.54893000000004*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1489.0978600000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 38.493319653319467*x2 - 13.537125818774067;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0040292852210532352*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2233.6467900000002*x3;
    double x26 = 0.012087855663159706*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 61905.242458709938*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 5.831495441225);

    result += x0*x1*x2*(308856272674.91479*x14*x18 - x14*x4 + 102952090891.63828*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(6700.9403700000003*x0*x31 - x26*x30 + 3.0*x27*x32) + 1663059.327492435*x16*x24 - 3326118.65498487*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 205904181783.27655*x18*exp(-2233.6467900000002*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(12419499058.616165*x2 - 3119730957.3416729);
    double x4 = 40.611377456322202*x2;
    double x5 = 57.739979479979205*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 744.54893000000004*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1489.0978600000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 38.493319653319467*x2 - 13.537125818774067;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0040292852210532352*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2233.6467900000002*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 5.831495441225;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.012087855663159706*x33));
    double x35 = 89.817745857745422*x2 - 22.561876364623444;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0040292852210532352*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(308856272674.91479*x13*x20 + x13*x3 + 102952090891.63828*x14*x20 + x14*x3 + 61905.242458709938*x15*x34 + 61905.242458709938*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(115.4799589599584*x2 - 40.611377456322202)/pow(x6, 5.0/2.0) + x24*x35 - 1663059.327492435*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(76.986639306638935*x2 - 27.074251637548134) + 3326118.65498487*x37*exp(x12)/((x26)*(x26)*(x26))) + 205904181783.27655*x20*exp(-2233.6467900000002*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 40.611377456322202*x2;
    double x5 = 57.739979479979205*x3;
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = 2.4818297666666669*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 38.493319653319467*x2 - 13.537125818774067;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.9636595333333338*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 5.831495441225;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 460914.8203402306*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = 299.39248619248474*x2 - 60.165003638995849;
    double x24 = x2*x22*x23;
    double x25 = 185715.72737612983*x21;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 744.54893000000004*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 185715.72737612983*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1489.0978600000001*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 138274446.10206917*x26;
    double x38 = x19*x37;
    double x39 = 76.986639306638935*x2 - 27.074251637548134;
    double x40 = ((x12)*(x12));
    double x41 = x0*x40;
    double x42 = x18*x39*x41;
    double x43 = 460914.8203402306*x42;
    double x44 = (115.4799589599584*x2 - 40.611377456322202)/pow(x6, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 89.817745857745422*x2 - 22.561876364623444;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = 460914.8203402306*x49;
    double x51 = x16*x50;
    double x52 = x47*(179.63549171549084*x2 - 45.123752729246888);
    double x53 = x11*x47;
    double x54 = 557147.18212838948*x46*x53;
    double x55 = x21*x50;
    double x56 = x37*x42;
    double x57 = x37*x49;
    double x58 = x36*x57;
    double x59 = x31*x57;
    double x60 = 61905.242458709945*x22;
    double x61 = exp(x8);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x8);
    double x65 = x22*x64;
    double x66 = -3*x63 + 1.2087855663159706*x65;
    double x67 = x2*x66;
    double x68 = exp(x28);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x28);
    double x72 = T*x22*x71;
    double x73 = -3*x70 + 0.0040292852210532352*x72;
    double x74 = x2*x73;
    double x75 = 61905.242458709938*x22;
    double x76 = x46*x66;
    double x77 = 185715.72737612983*x53;
    double x78 = x46*x73;
    double x79 = 1.2087855663159706*x64;
    double x80 = x11*x79;
    double x81 = x61/((x62)*(x62));
    double x82 = 7.4454893000000002*x81;
    double x83 = x22*x82;
    double x84 = x49*(-9.0*x63 + 3.6263566989479115*x65) - x80 + x83;
    double x85 = 123810.48491741989*x84;
    double x86 = x22*x48;
    double x87 = 0.0040292852210532352*T*x71;
    double x88 = x11*x87;
    double x89 = x68/((x69)*(x69));
    double x90 = 2233.6467900000002*x26*x89;
    double x91 = x22*x90;
    double x92 = x49*(-9.0*x70 + 0.012087855663159706*x72) - x88 + x91;
    double x93 = 123810.48491741988*x92;
    double x94 = x12*x2*x44;
    double x95 = x2*x40;
    double x96 = x49*x95;
    double x97 = x11*x95;
    double x98 = 3.0*x49;
    double x99 = 3.0*x18;
    double x100 = x12*x39;
    double x101 = x95*x98;
    double x102 = x33*x96;
    double x103 = x74*x99;

    result += (-13793540342.077885*x0 - 3431736.3630546099*x14*x16 - 1143912.1210182032*x14*x21 + x14*x85 - x14*x93 - x16*x20 - x16*x43 - 221673384.76961496*x2 - x20*x21 - x21*x43 - x21*x54 - x23*x60*x67 + x23*x74*x75 - x24*x25 + x24*x32 - x25*x45 + 3882628518.3829451*x3 + 102952090891.63828*x31*x34 + x31*x38 + x31*x54 + x31*x56 + x32*x45 + 308856272674.91479*x34*x36 + x36*x38 + x36*x56 - 61905.242458709945*x45*x66 + 61905.242458709938*x45*x73 - x47*x60*(x100*x67*x99 + x101*x84 + x46*x80 - x46*x83 + x66*x95*x99 - x76*x98 + x79*x94 - 18.478436972138166*x81*x96 - x82*x97 + 36.956873944276332*x96*exp(x15)/((x62)*(x62)*(x62))) + x47*x75*(x100*x103 + x101*x92 - 1663059.327492435*x102*x89 + 3326118.65498487*x102*exp(x35)/((x69)*(x69)*(x69)) + x103*x40 + x46*x88 - x46*x91 - x78*x98 + x87*x94 - x90*x97) + x48*x51 + x48*x55 - x48*x58 - x48*x59 + x51*x52 + x52*x55 - x52*x58 - x52*x59 - x76*x77 + x77*x78 + x85*x86 - x86*x93 + 205904181783.27655*x34*exp(-2233.6467900000002*x27)/((x30)*(x30)*(x30)) - 2287824.2420364064*x14*exp(-7.4454893000000002*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.8054;
static const double Vmax = 1.15*6.8054;
static double V = 0.9*6.8054;

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



const char *OrthoDiopside_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *OrthoDiopside_slb_em_coder_calib_name(void) {
    return "OrthoDiopside_slb_em";
}

const char *OrthoDiopside_slb_em_coder_calib_formula(void) {
    return "CaMgSi2O6";
}

const double OrthoDiopside_slb_em_coder_calib_mw(void) {
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

const double *OrthoDiopside_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double OrthoDiopside_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double OrthoDiopside_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double OrthoDiopside_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double OrthoDiopside_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double OrthoDiopside_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double OrthoDiopside_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double OrthoDiopside_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double OrthoDiopside_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double OrthoDiopside_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double OrthoDiopside_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double OrthoDiopside_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double OrthoDiopside_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double OrthoDiopside_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int OrthoDiopside_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **OrthoDiopside_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **OrthoDiopside_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void OrthoDiopside_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int OrthoDiopside_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double OrthoDiopside_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int OrthoDiopside_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double OrthoDiopside_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double OrthoDiopside_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

