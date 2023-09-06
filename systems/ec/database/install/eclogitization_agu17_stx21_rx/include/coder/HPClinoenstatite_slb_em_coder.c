
static char *identifier = "HPClinoenstatite_slb_em.emml:680cd5c84af4a1467f49e3d389afaf540081f630:Fri May 26 02:49:30 2023";



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
    double x3 = sqrt(-11.993992636550246*x1 + 38.668590279711047*x2 + 1.1144041129000004);
    double x4 = 2.7429867333333338*x3;
    double x5 = 822.89602000000002*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) + 35788734.403792508*x1 - 207099238.45293558*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) - 3844571.492297194 + 327459852.87336564/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(38.668590279711047*pow(x0, 4.0/3.0) - 11.993992636550246*pow(x0, 2.0/3.0) + 1.1144041129000004);
    double x2 = x1/T;
    double x3 = 822.89602000000002*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -205258.14590751243*x2*x5/x6 + 68419.38196917082*x2*(-0.003645661088505447*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-11.993992636550246*x1 + 38.668590279711047*x3 + 1.1144041129000004);
    double x6 = 2.7429867333333338*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(3.9979975455167485*x2 - 25.77906018647403*x4);
    double x10 = 822.89602000000002*x5/T;
    double x11 = exp(-x10);
    double x12 = 68419.38196917082*x9;

    result += 205258.14590751243*x11*x9/(-x11 + 1) + x12*(-1.093698326551634*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.003645661088505447*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 23859156.269195005*x2 + 276132317.9372474*x4 - 205258.14590751246*x7*x9/(-x7 + 1) - 654919705.74673128/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-11.993992636550246*x2 + 38.668590279711047*x3 + 1.1144041129000004);
    double x5 = x0*x4;
    double x6 = 822.89602000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2025858655.6787519*x2 + 6531361215.1407375*x3 + 188229665.17109793);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1645.79204*x5)/((x8)*(x8)) + 68419.38196917082*x4*(x0*(0.010936983265516341*T*x11 - 9.0/x13) + 0.003645661088505447*x11 - 2468.68806*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 38.668590279711047*pow(x1, 4.0/3.0) - 11.993992636550246*x2 + 1.1144041129000004;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 822.89602000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 168906111.33987126*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(25.77906018647403*x2 - 3.9979975455167485)*(68419.38196917082*x4*(-2468.68806*x0*x11*x12/((x13)*(x13)) + 0.003645661088505447*x10/pow(x3, 3.0/2.0) + (0.010936983265516341*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1645.79204*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -11.993992636550246*x2 + 38.668590279711047*x3 + 1.1144041129000004;
    double x5 = sqrt(x4);
    double x6 = 2.7429867333333338*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = x7/x8;
    double x10 = 205258.14590751246*x9;
    double x11 = 1.0/x5;
    double x12 = x2*(60.151140435106072*x2 - 6.6633292425279134);
    double x13 = x11*x12;
    double x14 = 1.0/x4;
    double x15 = x3*((25.77906018647403*x2 - 3.9979975455167485)*(25.77906018647403*x2 - 3.9979975455167485));
    double x16 = x14*x15;
    double x17 = 563020.37113290443*x16;
    double x18 = pow(x4, -3.0/2.0);
    double x19 = x15*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 822.89602000000002*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 205258.14590751243*x25;
    double x27 = 168906111.33987126*x16*x20;
    double x28 = exp(x6);
    double x29 = x28 - 1;
    double x30 = 1.0/x29;
    double x31 = Debye(x6);
    double x32 = x11*x31;
    double x33 = -3*x30 + 1.093698326551634*x32;
    double x34 = 68419.38196917082*x11;
    double x35 = x12*x34;
    double x36 = 68419.38196917082*x19;
    double x37 = exp(x22);
    double x38 = x37 - 1;
    double x39 = 1.0/x38;
    double x40 = Debye(x22);
    double x41 = T*x11*x40;
    double x42 = -3*x39 + 0.003645661088505447*x41;
    double x43 = x15*x34;

    result += x0*(1964759117.2401938*x0 - x10*x13 + x10*x19 + x13*x26 + x17*x9 + x17*exp(-5.4859734666666675*x5)/((x8)*(x8)) - x19*x26 + 39765260.44865834*x2 - x25*x27 - 644308741.85357726*x3 - x33*x35 + x33*x36 + x35*x42 - x36*x42 - x43*(0.003645661088505447*T*x18*x40 - 2468.68806*x11*x20*x37/((x38)*(x38)) + x14*(-9.0*x39 + 0.010936983265516341*x41)) + x43*(-8.2289602000000013*x11*x28/((x29)*(x29)) + x14*(-9.0000000000000018*x30 + 3.2810949796549025*x32) + 1.093698326551634*x18*x31) - x27*exp(-1645.79204*x21)/((x24)*(x24)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-6077575967.0362549*x2 + 19594083645.422211*x3 + 564688995.51329374);
    double x5 = 1.0/T;
    double x6 = -11.993992636550246*x2 + 38.668590279711047*x3 + 1.1144041129000004;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 822.89602000000002*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1645.79204*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.010936983265516341*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 68419.38196917082*x7;
    double x25 = x15*x6;

    result += x0*(-416976500325.77075*x13*x16 + x13*x4 - 138992166775.25693*x14*x16 + x14*x4 + x24*(0.003645661088505447*x17 - 2468.68806*x21 + x23) - x24*(2031473.5791955213*x20*x25 - 2468.6880599999995*x21 + x22 + 3.0*x23 - 4062947.1583910426*x25*exp(x12)/((x19)*(x19)*(x19))) - 277984333550.51385*x16*exp(-2468.68806*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(8708481620.1876507*x2 - 1350572437.1191678);
    double x5 = 38.668590279711047*pow(x1, 4.0/3.0) - 11.993992636550246*x2 + 1.1144041129000004;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 822.89602000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1645.79204*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 25.77906018647403*x2 - 3.9979975455167485;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.003645661088505447*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2468.68806*x22*x3;
    double x24 = 0.010936983265516341*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 68419.38196917082*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-416976500325.77075*x12*x16 + x12*x4 - 138992166775.25693*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7406.0641799999994*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 2031473.5791955213*x15*x22 - 4062947.1583910426*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 277984333550.51385*x16*exp(-2468.68806*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(10159895223.552258*x2 - 1125477030.9326398);
    double x4 = 38.668590279711047*pow(x1, 4.0/3.0) - 11.993992636550246*x2 + 1.1144041129000004;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 822.89602000000002*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1645.79204*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 25.77906018647403*x2 - 3.9979975455167485;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.003645661088505447*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2468.68806*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.010936983265516341*x30));
    double x32 = 60.151140435106072*x2 - 6.6633292425279134;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.003645661088505447*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-416976500325.77075*x11*x18 + x11*x3 - 138992166775.25693*x12*x18 + x12*x3 + 68419.38196917082*x13*x31 - 277984333550.51385*x18*exp(-2468.68806*x6)/((x9)*(x9)*(x9)) - 68419.38196917082*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(77.337180559422094*x2 - 11.993992636550246)/pow(x4, 5.0/2.0) - x22*x32 + 2031473.5791955213*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(51.558120372948061*x2 - 7.9959950910334969) - 4062947.1583910426*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -11.993992636550246*x2 + 38.668590279711047*x3 + 1.1144041129000004;
    double x5 = sqrt(x4);
    double x6 = 2.7429867333333338*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 25.77906018647403*x2 - 3.9979975455167485;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 5.4859734666666675*x5;
    double x16 = exp(-x15)/((x8)*(x8));
    double x17 = 563020.37113290443*x16;
    double x18 = x7/x8;
    double x19 = 563020.37113290443*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(200.50380145035356*x2 - 17.768877980074436);
    double x22 = x20*x21;
    double x23 = 205258.14590751246*x18;
    double x24 = 1.0/T;
    double x25 = x24*x5;
    double x26 = 822.89602000000002*x25;
    double x27 = exp(-x26);
    double x28 = -x27 + 1;
    double x29 = x27/x28;
    double x30 = 205258.14590751243*x29;
    double x31 = pow(T, -2);
    double x32 = x12*x31;
    double x33 = 1645.79204*x25;
    double x34 = exp(-x33)/((x28)*(x28));
    double x35 = 168906111.33987126*x24;
    double x36 = x14*x35;
    double x37 = x13*(51.558120372948061*x2 - 7.9959950910334969);
    double x38 = ((x10)*(x10));
    double x39 = x0*x38;
    double x40 = x37*x39;
    double x41 = (77.337180559422094*x2 - 11.993992636550246)/pow(x4, 5.0/2.0);
    double x42 = x39*x41;
    double x43 = 1.0/x4;
    double x44 = 60.151140435106072*x2 - 6.6633292425279134;
    double x45 = x43*x44;
    double x46 = x10*x3;
    double x47 = x17*x46;
    double x48 = x43*(120.30228087021214*x2 - 13.326658485055827);
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
    double x61 = -3*x58 + 1.093698326551634*x60;
    double x62 = 68419.38196917082*x20;
    double x63 = x21*x62;
    double x64 = exp(x26);
    double x65 = x64 - 1;
    double x66 = 1.0/x65;
    double x67 = Debye(x26);
    double x68 = T*x20*x67;
    double x69 = -3*x66 + 0.003645661088505447*x68;
    double x70 = 68419.38196917082*x42;
    double x71 = 205258.14590751246*x50;
    double x72 = 1.093698326551634*x59;
    double x73 = x72*x9;
    double x74 = x56/((x57)*(x57));
    double x75 = 8.2289602000000013*x74;
    double x76 = x20*x75;
    double x77 = x43*(-9.0000000000000018*x58 + 3.2810949796549025*x60) + x73 - x76;
    double x78 = 136838.76393834164*x77;
    double x79 = x20*x49;
    double x80 = 0.003645661088505447*T*x67;
    double x81 = x80*x9;
    double x82 = x64/((x65)*(x65));
    double x83 = 2468.68806*x24*x82;
    double x84 = x20*x83;
    double x85 = x43*(-9.0*x66 + 0.010936983265516341*x68) + x81 - x84;
    double x86 = 136838.76393834164*x85;
    double x87 = x10*x2;
    double x88 = x41*x87;
    double x89 = x2*x38;
    double x90 = x43*x89;
    double x91 = x89*x9;
    double x92 = 3.0000000000000004*x61;
    double x93 = x13*x89;
    double x94 = x37*x87;
    double x95 = x46*x62;
    double x96 = x31*x90;
    double x97 = 3.0*x69;

    result += (-7859036468.9607754*x0 + 4633072.2258419003*x12*x16 + 1544357.4086139668*x12*x18 + x12*x78 - x12*x86 + 3088714.8172279336*x12*exp(-8.2289602000000013*x5)/((x8)*(x8)*(x8)) + x14*x17 + x14*x19 + x17*x40 - 615774.43772253743*x18*x50 + x19*x40 - 106040694.52975556*x2 + x22*x23 - x22*x30 + x23*x42 - 138992166775.25693*x29*x32 - x29*x36 + 615774.43772253732*x29*x50 - x29*x52 + 2147695806.1785908*x3 - x30*x42 - 416976500325.77075*x32*x34 - x34*x36 - x34*x52 - x45*x47 - x45*x51 + x45*x54 + x45*x55 - x47*x48 - x48*x51 + x48*x54 + x48*x55 + x61*x63 + x61*x70 - x61*x71 - x63*x69 - x69*x70 + x69*x71 - x78*x79 + x79*x86 + x95*(-x44*x73 + x44*x76 - x45*x92 + x72*x88 + 22.571928657728019*x74*x90 - x75*x91 + 3.0000000000000004*x77*x90 + x92*x93 + x92*x94 - 45.143857315456039*x90*exp(x15)/((x57)*(x57)*(x57))) - x95*(-x44*x81 + x44*x84 - x45*x97 + x80*x88 + 2031473.5791955213*x82*x96 - x83*x91 + 3.0*x85*x90 + x93*x97 + x94*x97 - 4062947.1583910426*x96*exp(x33)/((x65)*(x65)*(x65))) - 277984333550.51385*x32*exp(-2468.68806*x25)/((x28)*(x28)*(x28)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.076;
static const double Vmax = 1.15*6.076;
static double V = 0.9*6.076;

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



const char *HPClinoenstatite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *HPClinoenstatite_slb_em_coder_calib_name(void) {
    return "HPClinoenstatite_slb_em";
}

const char *HPClinoenstatite_slb_em_coder_calib_formula(void) {
    return "Mg2Si2O6";
}

const double HPClinoenstatite_slb_em_coder_calib_mw(void) {
    return 200.7774;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,6.0,0.0,0.0,0.0,
        2.0,0.0,2.0,0.0,0.0,0.0,
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

const double *HPClinoenstatite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double HPClinoenstatite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double HPClinoenstatite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int HPClinoenstatite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **HPClinoenstatite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **HPClinoenstatite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void HPClinoenstatite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int HPClinoenstatite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double HPClinoenstatite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int HPClinoenstatite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double HPClinoenstatite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double HPClinoenstatite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

