
static char *identifier = "Nepheline_stx21_em.emml:f9daf0c488ceeebc1844cb15f6a6c757419660fb:Sat Mar  9 00:08:46 2024";



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
    double x3 = sqrt(9.0541789153976104*x1 - 4.040278199555674*x2 - 1.5185617549999999);
    double x4 = 2.4785995000000001*x3;
    double x5 = 743.57984999999996*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 19760488.757031649*x1 + 30361385.600990448*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 1220773.3132499999;
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-4.040278199555674*pow(x0, 4.0/3.0) + 9.0541789153976104*pow(x0, 2.0/3.0) - 1.5185617549999999);
    double x2 = x1/T;
    double x3 = 743.57984999999996*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -129831.80419517684*x2*x5/x6 + 43277.268065058946*x2*(-0.0040345364388236181*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(9.0541789153976104*x1 - 4.040278199555674*x3 - 1.5185617549999999);
    double x6 = 2.4785995000000001*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-3.01805963846587*x2 + 2.6935187997037824*x4);
    double x10 = 743.57984999999996*x5/T;
    double x11 = exp(-x10);

    result += 129831.80419517684*x11*x9/(-x11 + 1) + 13173659.171354432*x2 - 40481847.467987262*x4 - 129831.80419517685*x7*x9/(-x7 + 1) + 43277.268065058954*x9*(-1.2103609316470854*x8*Debye(x6) + 3/(exp(x6) - 1)) - 43277.268065058946*x9*(-0.0040345364388236181*T*x8*Debye(x10) + 3/(exp(x10) - 1));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(9.0541789153976104*x2 - 4.040278199555674*x3 - 1.5185617549999999);
    double x5 = x0*x4;
    double x6 = 743.57984999999996*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-874093270.8750726*x2 + 390049723.96658021*x3 + 146602427.8796185);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1487.1596999999999*x5)/((x8)*(x8)) - 43277.268065058946*x4*(x0*(0.012103609316470854*T*x11 - 9.0/x13) + 0.0040345364388236181*x11 - 2230.7395499999998*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 9.0541789153976104*x2;
    double x4 = 4.040278199555674*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 1.5185617549999999;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 743.57984999999996*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 96540313.488678962*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(2.6935187997037824*x2 - 3.01805963846587)*(-43277.268065058946*x6*(2230.7395499999998*x0*x13*x14/((x15)*(x15)) - 0.0040345364388236181*x12/pow(x5, 3.0/2.0) + (0.012103609316470854*x12*x13 - 9.0/x15)/(-x3 + x4 + 1.5185617549999999)) + x11*x9/x10 + x11*exp(-1487.1596999999999*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = 9.0541789153976104*x1;
    double x3 = 4.040278199555674*pow(x0, 4.0/3.0);
    double x4 = x2 - x3 - 1.5185617549999999;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = 2.4785995000000001*x5;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = x8/x9;
    double x11 = 1.0/T;
    double x12 = x11*x5;
    double x13 = 743.57984999999996*x12;
    double x14 = exp(-x13);
    double x15 = -x14 + 1;
    double x16 = x14/x15;
    double x17 = 1.0/(-x2 + x3 + 1.5185617549999999);
    double x18 = x1*((2.6935187997037824*x1 - 3.01805963846587)*(2.6935187997037824*x1 - 3.01805963846587));
    double x19 = x17*x18;
    double x20 = 321801.04496226326*x19;
    double x21 = pow(x4, -3.0/2.0);
    double x22 = x18*x21;
    double x23 = 96540313.488678962*x11*x19;
    double x24 = x6*(6.2848771993088253*x1 - 5.0300993974431165);
    double x25 = exp(x7);
    double x26 = x25 - 1;
    double x27 = 1.0/x26;
    double x28 = Debye(x7);
    double x29 = x28*x6;
    double x30 = -129831.80419517687*x27 + 52381.114494365407*x29;
    double x31 = exp(x13);
    double x32 = x31 - 1;
    double x33 = 1.0/x32;
    double x34 = Debye(x13);
    double x35 = T*x34*x6;
    double x36 = -129831.80419517684*x33 + 174.60371498121802*x35;
    double x37 = x18*x6;

    result += x1*(94457644.091970265*x1 - x10*x20 + 129831.80419517685*x10*x22 + x10*x6*(815976.94593139493*x1 - 653066.8800511118) - 129831.80419517684*x16*x22 + x16*x23 - x16*x6*(815976.94593139482*x1 - 653066.88005111169) - x20*exp(-4.9571990000000001*x5)/((x9)*(x9)) + x22*x30 - x22*x36 + x24*x30 - x24*x36 - 43277.268065058954*x37*(x17*(-9.0*x27 + 3.6310827949412561*x29) - 1.2103609316470854*x21*x28 + 7.4357985000000006*x25*x6/((x26)*(x26))) + 43277.268065058946*x37*(-0.0040345364388236181*T*x21*x34 + 2230.7395499999998*x11*x31*x6/((x32)*(x32)) + x17*(-9.0*x33 + 0.012103609316470854*x35)) - 21956098.618924052 + x23*exp(-1487.1596999999999*x12)/((x15)*(x15)))/((V)*(V));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-2622279812.6252179*x2 + 1170149171.8997407*x3 + 439807283.63885552);
    double x5 = 1.0/T;
    double x6 = 9.0541789153976104*x2;
    double x7 = 4.040278199555674*x3;
    double x8 = x6 - x7 - 1.5185617549999999;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 743.57984999999996*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1487.1596999999999*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.012103609316470854*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 43277.268065058946*x9;
    double x27 = x17*(-x6 + x7 + 1.5185617549999999);

    result += x0*(-215356295468.59467*x15*x18 - x15*x4 - 71785431822.864883*x16*x18 - x16*x4 + x26*(0.0040345364388236181*x19 - 2230.7395499999998*x23 + x25) - x26*(-1658732.9799780673*x22*x27 - 2230.7395499999993*x23 + x24 + 3.0*x25 + 3317465.9599561347*x27*exp(x14)/((x21)*(x21)*(x21))) - 143570863645.72977*x18*exp(-2230.7395499999998*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(520066298.62210685*x2 - 582728847.2500484);
    double x5 = 9.0541789153976104*x2;
    double x6 = 4.040278199555674*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 1.5185617549999999;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 743.57984999999996*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1487.1596999999999*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 2.6935187997037824*x2 - 3.01805963846587;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0040345364388236181*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2230.7395499999998*x3;
    double x26 = 0.012103609316470854*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 43277.268065058946*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 1.5185617549999999);

    result += x0*x1*x2*(215356295468.59467*x14*x18 - x14*x4 + 71785431822.864883*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(6692.2186499999989*x0*x31 - x26*x30 + 3.0*x27*x32) + 1658732.9799780673*x16*x24 - 3317465.9599561347*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 143570863645.72977*x18*exp(-2230.7395499999998*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(606744015.05912459*x2 - 485607372.70837361);
    double x4 = 9.0541789153976104*x2;
    double x5 = 4.040278199555674*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 1.5185617549999999;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 743.57984999999996*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1487.1596999999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 2.6935187997037824*x2 - 3.01805963846587;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0040345364388236181*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2230.7395499999998*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 1.5185617549999999;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.012103609316470854*x33));
    double x35 = 6.2848771993088253*x2 - 5.0300993974431165;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0040345364388236181*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(215356295468.59467*x13*x20 + x13*x3 + 71785431822.864883*x14*x20 + x14*x3 + 43277.268065058946*x15*x34 + 43277.268065058946*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(8.080556399111348*x2 - 9.0541789153976104)/pow(x6, 5.0/2.0) + x24*x35 - 1658732.9799780673*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(5.3870375994075648*x2 - 6.03611927693174) + 3317465.9599561347*x37*exp(x12)/((x26)*(x26)*(x26))) + 143570863645.72977*x20*exp(-2230.7395499999998*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = 9.0541789153976104*x1;
    double x4 = 4.040278199555674*x2;
    double x5 = x3 - x4 - 1.5185617549999999;
    double x6 = sqrt(x5);
    double x7 = 2.4785995000000001*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = pow(x5, -3.0/2.0);
    double x11 = pow(V, -2);
    double x12 = 2.6935187997037824*x1 - 3.01805963846587;
    double x13 = x11*((x12)*(x12)*(x12));
    double x14 = x10*x13;
    double x15 = 4.9571990000000001*x6;
    double x16 = exp(-x15)/((x9)*(x9));
    double x17 = -x3 + x4 + 1.5185617549999999;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 321801.04496226326*x19;
    double x21 = x8/x9;
    double x22 = 1.0/x6;
    double x23 = 20.949590664362752*x1 - 13.413598393181644;
    double x24 = x1*x22*x23;
    double x25 = 129831.80419517685*x21;
    double x26 = 1.0/T;
    double x27 = x26*x6;
    double x28 = 743.57984999999996*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = 129831.80419517684*x31;
    double x33 = pow(T, -2);
    double x34 = x14*x33;
    double x35 = 1487.1596999999999*x27;
    double x36 = exp(-x35)/((x30)*(x30));
    double x37 = 96540313.488678962*x26;
    double x38 = x19*x37;
    double x39 = 5.3870375994075648*x1 - 6.03611927693174;
    double x40 = ((x12)*(x12));
    double x41 = x11*x40;
    double x42 = x18*x39*x41;
    double x43 = 321801.04496226326*x42;
    double x44 = (8.080556399111348*x1 - 9.0541789153976104)/pow(x5, 5.0/2.0);
    double x45 = x41*x44;
    double x46 = 6.2848771993088253*x1 - 5.0300993974431165;
    double x47 = x12*x2;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = 321801.04496226326*x49;
    double x51 = x16*x50;
    double x52 = x47*(12.569754398617651*x1 - 10.060198794886233);
    double x53 = x10*x47;
    double x54 = 389495.41258553055*x46*x53;
    double x55 = x21*x50;
    double x56 = x37*x42;
    double x57 = x37*x49;
    double x58 = x36*x57;
    double x59 = x31*x57;
    double x60 = 43277.268065058954*x22;
    double x61 = exp(x7);
    double x62 = x61 - 1;
    double x63 = 1.0/x62;
    double x64 = Debye(x7);
    double x65 = x22*x64;
    double x66 = -3*x63 + 1.2103609316470854*x65;
    double x67 = x1*x66;
    double x68 = exp(x28);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x28);
    double x72 = T*x22*x71;
    double x73 = -3*x70 + 0.0040345364388236181*x72;
    double x74 = x1*x73;
    double x75 = 43277.268065058946*x22;
    double x76 = x46*x66;
    double x77 = x46*x73;
    double x78 = 1.2103609316470854*x64;
    double x79 = x10*x78;
    double x80 = x61/((x62)*(x62));
    double x81 = 7.4357985000000006*x80;
    double x82 = x22*x81;
    double x83 = x49*(-9.0*x63 + 3.6310827949412561*x65) - x79 + x82;
    double x84 = 86554.536130117907*x83;
    double x85 = x22*x48;
    double x86 = 0.0040345364388236181*T*x71;
    double x87 = x10*x86;
    double x88 = x68/((x69)*(x69));
    double x89 = 2230.7395499999998*x26*x88;
    double x90 = x22*x89;
    double x91 = x49*(-9.0*x70 + 0.012103609316470854*x72) - x87 + x90;
    double x92 = 86554.536130117893*x91;
    double x93 = x1*x12*x44;
    double x94 = x1*x40;
    double x95 = x49*x94;
    double x96 = x10*x94;
    double x97 = 3.0*x49;
    double x98 = 3.0*x18;
    double x99 = x12*x39;
    double x100 = x94*x97;
    double x101 = x33*x95;
    double x102 = x74*x98;

    result += (58549596.317130804*x1 - 2392847.7274288298*x14*x16 - 797615.90914294322*x14*x21 + x14*x84 - x14*x92 - 1595231.8182858864*x14*exp(-7.4357985000000006*x6)/((x9)*(x9)*(x9)) - x16*x20 - x16*x43 - 314858813.63990086*x2 - x20*x21 - x21*x43 - x21*x54 - x23*x60*x67 + x23*x74*x75 - x24*x25 + x24*x32 - x25*x45 + 71785431822.864883*x31*x34 + x31*x38 + x31*x54 + x31*x56 + x32*x45 + 215356295468.59467*x34*x36 + x36*x38 + x36*x56 - 43277.268065058954*x45*x66 + 43277.268065058946*x45*x73 - x47*x60*(x100*x83 + x46*x79 - x46*x82 + x66*x94*x98 + x67*x98*x99 - x76*x97 + x78*x93 - 18.430366444200754*x80*x95 - x81*x96 + 36.860732888401508*x95*exp(x15)/((x62)*(x62)*(x62))) + x47*x75*(x100*x91 - 1658732.9799780673*x101*x88 + 3317465.9599561347*x101*exp(x35)/((x69)*(x69)*(x69)) + x102*x40 + x102*x99 + x46*x87 - x46*x90 - x77*x97 + x86*x93 - x89*x96) + x48*x51 + x48*x55 - x48*x58 - x48*x59 + x51*x52 + x52*x55 - x52*x58 - x52*x59 - 129831.80419517687*x53*x76 + 129831.80419517684*x53*x77 + x84*x85 - x85*x92 + 143570863645.72977*x34*exp(-2230.7395499999998*x27)/((x30)*(x30)*(x30)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*5.3868;
static const double Vmax = 1.15*5.3868;
static double V = 0.9*5.3868;

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



const char *Nepheline_stx21_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Nepheline_stx21_em_coder_calib_name(void) {
    return "Nepheline_stx21_em";
}

const char *Nepheline_stx21_em_coder_calib_formula(void) {
    return "NaAlSiO4";
}

const double Nepheline_stx21_em_coder_calib_mw(void) {
    return 142.05441;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,4.0,0.0,0.0,1.0,
        0.0,1.0,1.0,0.0,0.0,0.0,
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

const double *Nepheline_stx21_em_coder_calib_elements(void) {
    return elmformula;
}

double Nepheline_stx21_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Nepheline_stx21_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Nepheline_stx21_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Nepheline_stx21_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Nepheline_stx21_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Nepheline_stx21_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Nepheline_stx21_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Nepheline_stx21_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Nepheline_stx21_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Nepheline_stx21_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Nepheline_stx21_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Nepheline_stx21_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Nepheline_stx21_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Nepheline_stx21_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Nepheline_stx21_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Nepheline_stx21_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Nepheline_stx21_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Nepheline_stx21_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Nepheline_stx21_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Nepheline_stx21_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Nepheline_stx21_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Nepheline_stx21_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Nepheline_stx21_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Nepheline_stx21_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Nepheline_stx21_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Nepheline_stx21_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Nepheline_stx21_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Nepheline_stx21_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

