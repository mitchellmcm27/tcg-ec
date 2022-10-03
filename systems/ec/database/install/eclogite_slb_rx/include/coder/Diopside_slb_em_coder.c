
static char *identifier = "Diopside_slb_em.emml:e14730139f27616c46c70bed48ef5633c5bc8c43:Thu Feb 10 16:51:54 2022";



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
    double x3 = sqrt(14.341527253702218*x1 - 7.4227312191405694*x2 - 2.4752859960500002);
    double x4 = 2.6053820000000001*x3;
    double x5 = 781.6146*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) - 4159928.0909287706*x1 - 88674649.21226716*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) + 144023.45621367078 + 225265874.2487379/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-7.4227312191405694*pow(x0, 4.0/3.0) + 14.341527253702218*pow(x0, 2.0/3.0) - 2.4752859960500002);
    double x2 = x1/T;
    double x3 = 781.6146*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -194961.16120508392*x2*x5/x6 + 64987.053735027977*x2*(-0.0038382087540329976*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(14.341527253702218*x1 - 7.4227312191405694*x3 - 2.4752859960500002);
    double x6 = 2.6053820000000001*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-4.7805090845674059*x2 + 4.9484874794270457*x4);
    double x10 = 194961.16120508392*x9;
    double x11 = 781.6146*x5/T;
    double x12 = exp(-x11);
    double x13 = 64987.053735027977*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.1514626262098995*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0038382087540329976*T*x8*Debye(x11) + 3/(exp(x11) - 1)) + 2773285.3939525136*x2 + 118232865.61635621*x4 - 450531748.4974758/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(14.341527253702218*x2 - 7.4227312191405694*x3 - 2.4752859960500002);
    double x5 = x0*x4;
    double x6 = 781.6146*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2185426316.8189087*x2 + 1131109111.4647844*x3 + 377195194.18857694);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1563.2292*x5)/((x8)*(x8)) - 64987.053735027977*x4*(x0*(0.011514626262098992*T*x11 - 8.9999999999999982/x13) + 0.0038382087540329976*x11 - 2344.8438000000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 14.341527253702218*x2;
    double x4 = 7.4227312191405694*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 2.4752859960500002;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 781.6146*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 152384490.03084719*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(4.9484874794270457*x2 - 4.7805090845674059)*(-64987.053735027977*x6*(2344.8438000000001*x0*x13*x14/((x15)*(x15)) - 0.0038382087540329976*x12/pow(x5, 3.0/2.0) + (0.011514626262098992*x12*x13 - 8.9999999999999982/x15)/(-x3 + x4 + 2.4752859960500002)) + x11*x9/x10 + x11*exp(-1563.2292*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 14.341527253702218*x2;
    double x5 = 7.4227312191405694*x3;
    double x6 = x4 - x5 - 2.4752859960500002;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(11.546470785329774*x2 - 7.9675151409456761);
    double x10 = 194961.16120508392*x9;
    double x11 = 2.6053820000000001*x7;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 1.0/(-x4 + x5 + 2.4752859960500002);
    double x16 = x3*((4.9484874794270457*x2 - 4.7805090845674059)*(4.9484874794270457*x2 - 4.7805090845674059));
    double x17 = x15*x16;
    double x18 = 507948.30010282394*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 194961.16120508392*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 781.6146*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 152384490.03084719*x17*x22;
    double x29 = exp(x11);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x11);
    double x33 = x32*x8;
    double x34 = -194961.16120508395*x31 + 74830.163563379174*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x38*x8;
    double x40 = -194961.16120508395*x37 + 249.43387854459721*x39;
    double x41 = 64987.053735027977*x16*x8;

    result += x0*(1351595245.4924273*x0 + x10*x14 - x10*x27 - x14*x18 + x14*x21 - 4622142.3232541895*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 - 275876686.43816447*x3 + x34*x9 - x40*x9 - x41*(x15*(-9.0000000000000018*x31 + 3.4543878786296989*x33) - 1.1514626262098995*x19*x32 + 7.8161459999999998*x29*x8/((x30)*(x30))) + x41*(x15*(-8.9999999999999982*x37 + 0.011514626262098992*x39) - 0.0038382087540329976*x19*x38 + 2344.8438000000001*x22*x35*x8/((x36)*(x36))) + x28*exp(-1563.2292*x23)/((x26)*(x26)) - x18*exp(-5.2107640000000002*x7)/((x13)*(x13)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-6556278950.456727*x2 + 3393327334.3943529*x3 + 1131585582.5657308);
    double x5 = 1.0/T;
    double x6 = 14.341527253702218*x2;
    double x7 = 7.4227312191405694*x3;
    double x8 = x6 - x7 - 2.4752859960500002;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 781.6146*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1563.2292*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.011514626262098992*x19;
    double x25 = x5*(T*x24 - 8.9999999999999982/x21);
    double x26 = 64987.053735027977*x9;
    double x27 = x17*(-x6 + x7 + 2.4752859960500002);

    result += x0*(-357317826664.99384*x15*x18 - x15*x4 - 119105942221.66461*x16*x18 - x16*x4 + x26*(0.0038382087540329976*x19 - 2344.8438000000001*x23 + x25) - x26*(-1832764.1487994802*x22*x27 - 2344.8437999999987*x23 + x24 + 2.9999999999999996*x25 + 3665528.2975989603*x27*exp(x14)/((x21)*(x21)*(x21))) - 238211884443.32922*x18*exp(-2344.8438000000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(1508145481.9530456*x2 - 1456950877.8792727);
    double x5 = 14.341527253702218*x2;
    double x6 = 7.4227312191405694*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 2.4752859960500002;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 781.6146*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1563.2292*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 4.9484874794270457*x2 - 4.7805090845674059;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0038382087540329976*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2344.8438000000001*x3;
    double x26 = 0.011514626262098992*T*x20;
    double x27 = x19*x26 - 8.9999999999999982/x23;
    double x28 = x0*x27;
    double x29 = 64987.053735027977*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 2.4752859960500002);

    result += x0*x1*x2*(357317826664.99384*x14*x18 - x14*x4 + 119105942221.66461*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7034.5313999999989*x0*x31 - x26*x30 + 2.9999999999999996*x27*x32) + 1832764.1487994802*x16*x24 - 3665528.2975989603*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 238211884443.32922*x18*exp(-2344.8438000000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(1759503062.2785532*x2 - 1214125731.5660605);
    double x4 = 14.341527253702218*x2;
    double x5 = 7.4227312191405694*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 2.4752859960500002;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 781.6146*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1563.2292*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 4.9484874794270457*x2 - 4.7805090845674059;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0038382087540329976*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2344.8438000000001*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 2.4752859960500002;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-8.9999999999999982*x32 + 0.011514626262098992*x33));
    double x35 = 11.546470785329774*x2 - 7.9675151409456761;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0038382087540329976*x33;
    double x40 = 2.9999999999999996*x31;
    double x41 = 2.9999999999999996*x39/((x30)*(x30));

    result += -x38*(357317826664.99384*x13*x20 + x13*x3 + 119105942221.66461*x14*x20 + x14*x3 + 64987.053735027977*x15*x34 + 64987.053735027977*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(14.845462438281137*x2 - 14.341527253702218)/pow(x6, 5.0/2.0) + x24*x35 - 1832764.1487994802*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(9.8969749588540914*x2 - 9.5610181691348117) + 3665528.2975989603*x37*exp(x12)/((x26)*(x26)*(x26))) + 238211884443.32922*x20*exp(-2344.8438000000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 14.341527253702218*x2;
    double x5 = 7.4227312191405694*x3;
    double x6 = x4 - x5 - 2.4752859960500002;
    double x7 = sqrt(x6);
    double x8 = 2.6053820000000001*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 4.9484874794270457*x2 - 4.7805090845674059;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.2107640000000002*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 2.4752859960500002;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 507948.30010282394*x16;
    double x21 = x9/x10;
    double x22 = 507948.30010282394*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(38.488235951099242*x2 - 21.246707042521802);
    double x25 = 194961.16120508392*x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 781.6146*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = pow(T, -2);
    double x33 = x14*x32;
    double x34 = 1563.2292*x27;
    double x35 = exp(-x34)/((x30)*(x30));
    double x36 = 152384490.03084719*x26;
    double x37 = x19*x36;
    double x38 = x18*(9.8969749588540914*x2 - 9.5610181691348117);
    double x39 = ((x12)*(x12));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = (14.845462438281137*x2 - 14.341527253702218)/pow(x6, 5.0/2.0);
    double x43 = x40*x42;
    double x44 = 194961.16120508392*x43;
    double x45 = 1.0/x17;
    double x46 = 11.546470785329774*x2 - 7.9675151409456761;
    double x47 = x45*x46;
    double x48 = x12*x3;
    double x49 = x20*x48;
    double x50 = x45*(23.092941570659548*x2 - 15.935030281891352);
    double x51 = x46*x48;
    double x52 = x11*x51;
    double x53 = 584883.48361525172*x52;
    double x54 = x22*x48;
    double x55 = x36*x41;
    double x56 = x36*x48;
    double x57 = x35*x56;
    double x58 = x31*x56;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.1514626262098995*x63;
    double x65 = 64987.053735027977*x23;
    double x66 = x24*x65;
    double x67 = exp(x28);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x28);
    double x71 = T*x23*x70;
    double x72 = -3*x69 + 0.0038382087540329976*x71;
    double x73 = 64987.053735027977*x43;
    double x74 = 194961.16120508395*x52;
    double x75 = 1.1514626262098995*x62;
    double x76 = x11*x75;
    double x77 = x59/((x60)*(x60));
    double x78 = 7.8161459999999998*x77;
    double x79 = x23*x78;
    double x80 = x45*(-9.0000000000000018*x61 + 3.4543878786296989*x63) - x76 + x79;
    double x81 = 129974.10747005595*x80;
    double x82 = x23*x51;
    double x83 = 0.0038382087540329976*T*x70;
    double x84 = x11*x83;
    double x85 = x67/((x68)*(x68));
    double x86 = 2344.8438000000001*x26*x85;
    double x87 = x23*x86;
    double x88 = x45*(-8.9999999999999982*x69 + 0.011514626262098992*x71) - x84 + x87;
    double x89 = 129974.10747005595*x88;
    double x90 = x12*x2;
    double x91 = x42*x90;
    double x92 = x2*x39;
    double x93 = x45*x92;
    double x94 = x11*x92;
    double x95 = 3.0000000000000004*x64;
    double x96 = x18*x92;
    double x97 = x38*x90;
    double x98 = x48*x65;
    double x99 = x32*x93;
    double x100 = 2.9999999999999996*x72;

    result += (-5406380981.9697094*x0 - 3970198.0740554868*x14*x16 - 1323399.3580184956*x14*x21 + x14*x81 - x14*x89 - x19*x20 - x19*x22 + 12325712.862011172*x2 - x20*x41 - x21*x25 - x21*x44 - x21*x53 - x22*x41 + x25*x31 + 919588954.79388154*x3 + 119105942221.66461*x31*x33 + x31*x37 + x31*x44 + x31*x53 + x31*x55 + 357317826664.99384*x33*x35 + x35*x37 + x35*x55 + x47*x49 + x47*x54 - x47*x57 - x47*x58 + x49*x50 + x50*x54 - x50*x57 - x50*x58 - x64*x66 - x64*x73 - x64*x74 + x66*x72 + x72*x73 + x72*x74 + x81*x82 - x82*x89 + x98*(-x100*x47 + x100*x96 + x100*x97 + x46*x84 - x46*x87 + x83*x91 - 1832764.1487994802*x85*x99 - x86*x94 + 2.9999999999999996*x88*x93 + 3665528.2975989603*x99*exp(x34)/((x68)*(x68)*(x68))) - x98*(x46*x76 - x46*x79 - x47*x95 + x75*x91 - 20.364046097772*x77*x93 - x78*x94 + 3.0000000000000004*x80*x93 + x95*x96 + x95*x97 + 40.728092195544001*x93*exp(x15)/((x60)*(x60)*(x60))) + 238211884443.32922*x33*exp(-2344.8438000000001*x27)/((x30)*(x30)*(x30)) - 2646798.7160369912*x14*exp(-7.8161459999999998*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
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

