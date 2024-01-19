
static char *identifier = "AlPostPerovskite_slb_em.emml:90340eb0ef433418519592fd89f54a804c208fa0:Tue Aug 22 19:32:36 2023";



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
    double x3 = sqrt(-5.9966012499579247*x1 + 14.372501747429684*x2 - 0.15158896324999915);
    double x4 = 2.4097945000000003*x3;
    double x5 = 722.93835000000001*x3/T;

    result += 124.7169392722986*T*log(1 - exp(-x5)) - 41.572313090766201*T*Debye(x5) - 23847239.755976487*x1 + 21282891.132854495*x2 - 37415.08178168958*log(1 - exp(-x4)) + 12471.69392722986*Debye(x4) + 5343676.1449999996;
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(14.372501747429684*pow(x0, 4.0/3.0) - 5.9966012499579247*pow(x0, 2.0/3.0) - 0.15158896324999915);
    double x2 = x1/T;
    double x3 = 722.93835000000001*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -90162.658294565757*x2*x5/x6 + 30054.219431521917*x2*(-0.0041497314397555482*T*x4/x1 + 3/(exp(x3) - 1)) - 41.572313090766201*x4 + 124.7169392722986*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-5.9966012499579247*x1 + 14.372501747429684*x3 - 0.15158896324999915);
    double x6 = 2.4097945000000003*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(1.9988670833193081*x2 - 9.5816678316197894*x4);
    double x10 = 90162.658294565757*x9;
    double x11 = 722.93835000000001*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 15898159.837317657*x2 - 28377188.177139327*x4 + 30054.21943152192*x9*(-1.2449194319266641*x8*Debye(x6) + 3/(exp(x6) - 1)) - 30054.219431521917*x9*(-0.0041497314397555482*T*x8*Debye(x11) + 3/(exp(x11) - 1));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-5.9966012499579247*x2 + 14.372501747429684*x3 - 0.15158896324999915);
    double x5 = x0*x4;
    double x6 = 722.93835000000001*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(390870723.04170996*x2 - 936829032.94186819*x3 + 9880878.3844158575);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1445.8767*x5)/((x8)*(x8)) - 30054.219431521917*x4*(x0*(0.012449194319266646*T*x11 - 9.0000000000000018/x13) + 0.0041497314397555482*x11 - 2168.8150500000002*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 5.9966012499579247*x2;
    double x4 = 14.372501747429684*pow(x1, 4.0/3.0);
    double x5 = -x3 + x4 - 0.15158896324999915;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 722.93835000000001*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 65182043.419087186*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(9.5816678316197894*x2 - 1.9988670833193081)*(30054.219431521917*x6*(2168.8150500000002*x0*x13*x14/((x15)*(x15)) - 0.0041497314397555482*x12/pow(x5, 3.0/2.0) + (0.012449194319266646*x12*x13 - 9.0000000000000018/x15)/(x3 - x4 + 0.15158896324999915)) - x11*x9/x10 - x11*exp(-1445.8767*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = 5.9966012499579247*x1;
    double x3 = 14.372501747429684*pow(x0, 4.0/3.0);
    double x4 = -x2 + x3 - 0.15158896324999915;
    double x5 = sqrt(x4);
    double x6 = 1.0/x5;
    double x7 = x6*(2015786.8327201915*x1 - 300371.94968262344);
    double x8 = 2.4097945000000003*x5;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/T;
    double x13 = x12*x5;
    double x14 = 722.93835000000001*x13;
    double x15 = exp(-x14);
    double x16 = -x15 + 1;
    double x17 = x15/x16;
    double x18 = 1.0/(x2 - x3 + 0.15158896324999915);
    double x19 = x1*((9.5816678316197894*x1 - 1.9988670833193081)*(9.5816678316197894*x1 - 1.9988670833193081));
    double x20 = x18*x19;
    double x21 = 217273.47806362397*x20;
    double x22 = pow(x4, -3.0/2.0);
    double x23 = x19*x22;
    double x24 = 90162.658294565757*x23;
    double x25 = 65182043.419087186*x12*x20;
    double x26 = x6*(22.357224940446173*x1 - 3.3314451388655133);
    double x27 = exp(x8);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x8);
    double x31 = x30*x6;
    double x32 = -90162.658294565757*x29 + 37415.08178168958*x31;
    double x33 = exp(x14);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x14);
    double x37 = T*x36*x6;
    double x38 = -90162.658294565743*x35 + 124.71693927229862*x37;
    double x39 = x19*x6;

    result += x1*(66213439.079991758*x1 - x11*x21 + x11*x24 - x11*x7 - x17*x24 + x17*x25 + x17*x7 + x23*x32 - x23*x38 - x26*x32 + x26*x38 - 30054.21943152192*x39*(x18*(-9.0*x29 + 3.7347582957799923*x31) - 1.2449194319266641*x22*x30 + 7.2293835000000009*x27*x6/((x28)*(x28))) + 30054.219431521917*x39*(-0.0041497314397555482*T*x22*x36 + 2168.8150500000002*x12*x33*x6/((x34)*(x34)) + x18*(-9.0000000000000018*x35 + 0.012449194319266646*x37)) - 26496933.062196095 + x25*exp(-1445.8767*x13)/((x16)*(x16)) - x21*exp(-4.8195890000000006*x5)/((x10)*(x10)))/((V)*(V));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(1172612169.1251299*x2 - 2810487098.8256044*x3 + 29642635.153247572);
    double x5 = 1.0/T;
    double x6 = 5.9966012499579247*x2;
    double x7 = 14.372501747429684*x3;
    double x8 = -x6 + x7 - 0.15158896324999915;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 722.93835000000001*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1445.8767*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.012449194319266646*x19;
    double x25 = x5*(T*x24 - 9.0000000000000018/x21);
    double x26 = 30054.219431521917*x9;
    double x27 = x17*(x6 - x7 + 0.15158896324999915);

    result += x0*(-141367796757.06973*x15*x18 - x15*x4 - 47122598919.023247*x16*x18 - x16*x4 + x26*(0.0041497314397555482*x19 - 2168.8150500000002*x23 + x25) - x26*(-1567919.5737021677*x22*x27 - 2168.8150500000011*x23 + x24 + 3.0000000000000004*x25 + 3135839.1474043354*x27*exp(x14)/((x21)*(x21)*(x21))) - 94245197838.046494*x18*exp(-2168.8150500000002*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(1249105377.2558241*x2 - 260580482.02780661);
    double x5 = 5.9966012499579247*x2;
    double x6 = 14.372501747429684*pow(x1, 4.0/3.0);
    double x7 = -x5 + x6 - 0.15158896324999915;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 722.93835000000001*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1445.8767*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = 9.5816678316197894*x2 - 1.9988670833193081;
    double x17 = pow(T, -3);
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0041497314397555482*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2168.8150500000002*x24*x3;
    double x26 = 0.012449194319266646*T*x20;
    double x27 = x19*x26 - 9.0000000000000018/x23;
    double x28 = x0*x27;
    double x29 = 30054.219431521917*x16;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = 1.0/(x5 - x6 + 0.15158896324999915);

    result += x0*x1*x2*(-141367796757.06973*x14*x18 + x14*x4 - 47122598919.023247*x15*x18 + x15*x4 + x19*x29*(x19*x21 - x25*x8 + x28) - x29*x8*(x0*(-6506.4451500000014*x0*x19*x24 + x26*x30 - 3.0000000000000004*x27*x31) + 1567919.5737021677*x17*x24 - 3135839.1474043354*x17*exp(x13)/((x23)*(x23)*(x23)) + x19*x25 + x21*x30 - x28*x31) - 94245197838.046494*x18*exp(-2168.8150500000002*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(1457289606.7984614*x2 - 217150401.68983883);
    double x4 = 5.9966012499579247*x2;
    double x5 = 14.372501747429684*pow(x1, 4.0/3.0);
    double x6 = -x4 + x5 - 0.15158896324999915;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 722.93835000000001*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1445.8767*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -2);
    double x16 = 1.0/x7;
    double x17 = 9.5816678316197894*x2 - 1.9988670833193081;
    double x18 = ((x17)*(x17));
    double x19 = x18*x2;
    double x20 = x16*x19;
    double x21 = x15*x20;
    double x22 = pow(x6, -3.0/2.0);
    double x23 = Debye(x9);
    double x24 = 0.0041497314397555482*T*x23;
    double x25 = x22*x24;
    double x26 = exp(x9);
    double x27 = x26 - 1;
    double x28 = x26/((x27)*(x27));
    double x29 = 2168.8150500000002*x28;
    double x30 = x0*x16*x29;
    double x31 = x4 - x5 + 0.15158896324999915;
    double x32 = 1.0/x31;
    double x33 = 1.0/x27;
    double x34 = T*x16*x23;
    double x35 = x32*(-9.0000000000000018*x33 + 0.012449194319266646*x34);
    double x36 = 22.357224940446173*x2 - 3.3314451388655133;
    double x37 = x17*x2;
    double x38 = x19*x32;
    double x39 = x15*x38;
    double x40 = x0*x2;
    double x41 = -9.0000000000000018*x33 + 0.012449194319266646*x34;
    double x42 = x41/((x31)*(x31));

    result += x40*(-141367796757.06973*x13*x21 + x13*x3 - 47122598919.023247*x14*x21 + x14*x3 - 30054.219431521917*x20*(-x25 + x30 + x35) - 30054.219431521917*x7*(-x18*x22*x29*x40 + x19*x42 + x24*x37*(28.745003494859368*x2 - 5.9966012499579247)/pow(x6, 5.0/2.0) - x25*x36 - 1567919.5737021677*x28*x39 + x30*x36 + x32*x36*x41 + x37*x42*(19.163335663239579*x2 - 3.9977341666386161) - 3.0000000000000004*x38*(x25 - x30 - x35) + 3135839.1474043354*x39*exp(x12)/((x27)*(x27)*(x27))) - 94245197838.046494*x21*exp(-2168.8150500000002*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = 5.9966012499579247*x1;
    double x4 = 14.372501747429684*x2;
    double x5 = -x3 + x4 - 0.15158896324999915;
    double x6 = sqrt(x5);
    double x7 = 2.4097945000000003*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = pow(x5, -3.0/2.0);
    double x11 = pow(V, -2);
    double x12 = 9.5816678316197894*x1 - 1.9988670833193081;
    double x13 = x11*((x12)*(x12)*(x12));
    double x14 = x10*x13;
    double x15 = x3 - x4 + 0.15158896324999915;
    double x16 = pow(x15, -2);
    double x17 = x13*x16;
    double x18 = 4.8195890000000006*x6;
    double x19 = exp(-x18)/((x9)*(x9));
    double x20 = 217273.47806362397*x19;
    double x21 = x8/x9;
    double x22 = 217273.47806362397*x21;
    double x23 = 1.0/x6;
    double x24 = x1*(74.524083134820572*x1 - 8.8838537036413676);
    double x25 = 90162.658294565757*x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x6;
    double x28 = 722.93835000000001*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = pow(T, -2);
    double x33 = x14*x32;
    double x34 = 1445.8767*x27;
    double x35 = exp(-x34)/((x30)*(x30));
    double x36 = 65182043.419087186*x26;
    double x37 = x17*x36;
    double x38 = x16*(19.163335663239579*x1 - 3.9977341666386161);
    double x39 = ((x12)*(x12));
    double x40 = x11*x39;
    double x41 = x38*x40;
    double x42 = (28.745003494859368*x1 - 5.9966012499579247)/pow(x5, 5.0/2.0);
    double x43 = x40*x42;
    double x44 = 90162.658294565757*x43;
    double x45 = 1.0/x15;
    double x46 = 22.357224940446173*x1 - 3.3314451388655133;
    double x47 = x45*x46;
    double x48 = x12*x2;
    double x49 = x20*x48;
    double x50 = x45*(44.714449880892346*x1 - 6.6628902777310266);
    double x51 = x22*x48;
    double x52 = x46*x48;
    double x53 = x10*x52;
    double x54 = 270487.97488369729*x53;
    double x55 = x36*x41;
    double x56 = x36*x48;
    double x57 = x35*x56;
    double x58 = x31*x56;
    double x59 = exp(x7);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x7);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.2449194319266641*x63;
    double x65 = 30054.21943152192*x23;
    double x66 = exp(x28);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = T*Debye(x28);
    double x70 = x23*x69;
    double x71 = -3*x68 + 0.0041497314397555482*x70;
    double x72 = 30054.219431521917*x23;
    double x73 = 1.2449194319266641*x62;
    double x74 = x10*x73;
    double x75 = x59/((x60)*(x60));
    double x76 = 7.2293835000000009*x75;
    double x77 = x23*x76;
    double x78 = x45*(-9.0*x61 + 3.7347582957799923*x63) - x74 + x77;
    double x79 = 60108.438863043841*x78;
    double x80 = x23*x52;
    double x81 = 0.0041497314397555482*x69;
    double x82 = x10*x81;
    double x83 = x66/((x67)*(x67));
    double x84 = 2168.8150500000002*x26*x83;
    double x85 = x23*x84;
    double x86 = x45*(-9.0000000000000018*x68 + 0.012449194319266646*x70) - x82 + x85;
    double x87 = 60108.438863043833*x86;
    double x88 = x1*x12;
    double x89 = x42*x88;
    double x90 = x1*x39;
    double x91 = x45*x90;
    double x92 = x10*x90;
    double x93 = 3.0*x64;
    double x94 = x16*x90;
    double x95 = x38*x88;
    double x96 = x32*x91;
    double x97 = 3.0000000000000004*x71;

    result += (70658488.165856242*x1 + 1570753.2973007753*x14*x19 + 523584.43243359175*x14*x21 - x14*x79 + x14*x87 + 1047168.8648671835*x14*exp(-7.2293835000000009*x6)/((x9)*(x9)*(x9)) + x17*x20 + x17*x22 - 220711463.59997252*x2 + x20*x41 + x21*x25 + x21*x44 - x21*x54 + x22*x41 + x24*x64*x65 - x24*x71*x72 - x25*x31 - 47122598919.023247*x31*x33 - x31*x37 - x31*x44 + x31*x54 - x31*x55 - 141367796757.06973*x33*x35 - x35*x37 - x35*x55 + 30054.21943152192*x43*x64 - 30054.219431521917*x43*x71 + x47*x49 + x47*x51 - x47*x57 - x47*x58 + x48*x65*(-x46*x74 + x46*x77 + x47*x93 + x73*x89 - 17.421328596690753*x75*x91 - x76*x92 + 3.0*x78*x91 + x93*x94 + x93*x95 + 34.842657193381505*x91*exp(x18)/((x60)*(x60)*(x60))) - x48*x72*(-x46*x82 + x46*x85 + x47*x97 + x81*x89 - 1567919.5737021677*x83*x96 - x84*x92 + 3.0000000000000004*x86*x91 + x94*x97 + x95*x97 + 3135839.1474043354*x96*exp(x34)/((x67)*(x67)*(x67))) + x49*x50 + x50*x51 - x50*x57 - x50*x58 - 90162.658294565757*x53*x64 + 90162.658294565743*x53*x71 + x79*x80 - x80*x87 - 94245197838.046494*x33*exp(-2168.8150500000002*x27)/((x30)*(x30)*(x30)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.3847;
static const double Vmax = 1.15*2.3847;
static double V = 0.9*2.3847;

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



const char *AlPostPerovskite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *AlPostPerovskite_slb_em_coder_calib_name(void) {
    return "AlPostPerovskite_slb_em";
}

const char *AlPostPerovskite_slb_em_coder_calib_formula(void) {
    return "Al2O3";
}

const double AlPostPerovskite_slb_em_coder_calib_mw(void) {
    return 101.96127999999999;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,3.0,0.0,0.0,0.0,
        0.0,2.0,0.0,0.0,0.0,0.0,
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

const double *AlPostPerovskite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double AlPostPerovskite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double AlPostPerovskite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int AlPostPerovskite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **AlPostPerovskite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **AlPostPerovskite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void AlPostPerovskite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int AlPostPerovskite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double AlPostPerovskite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int AlPostPerovskite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double AlPostPerovskite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double AlPostPerovskite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

