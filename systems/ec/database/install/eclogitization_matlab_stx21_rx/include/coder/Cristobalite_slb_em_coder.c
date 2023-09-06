
static char *identifier = "Cristobalite_slb_em.emml:257d40464dd2113401b38015c97b50484b48960f:Fri May 26 02:48:53 2023";



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
    double x3 = sqrt(0.30624107406733952*x1 - 0.21220319221710396*x2 + 0.89877035125000004);
    double x4 = 2.6491484000000001*x3;
    double x5 = 794.74451999999997*x3/T;

    result += 74.830163563379159*T*log(1 - exp(-x5)) - 24.943387854459719*T*Debye(x5) + 3780578.18812097*x1 - 8714590.0180236064*x2 - 22449.049069013745*log(1 - exp(-x4)) + 7483.0163563379156*Debye(x4) - 1381165.6800000002 + 6537499.6886999998/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-0.21220319221710396*pow(x0, 4.0/3.0) + 0.30624107406733952*pow(x0, 2.0/3.0) + 0.89877035125000004);
    double x2 = x1/T;
    double x3 = 794.74451999999997*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -59470.862422699254*x2*x5/x6 + 19823.620807566418*x2*(-0.0037747979690378993*T*x4/x1 + 3/(exp(x3) - 1)) - 24.943387854459719*x4 + 74.830163563379159*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(0.30624107406733952*x1 - 0.21220319221710396*x3 + 0.89877035125000004);
    double x6 = 2.6491484000000001*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-0.1020803580224465*x2 + 0.14146879481140262*x4);
    double x10 = 59470.862422699254*x9;
    double x11 = 794.74451999999997*x5/T;
    double x12 = exp(-x11);
    double x13 = 19823.620807566418*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.1324393907113697*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0037747979690378993*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 2520385.4587473134*x2 + 11619453.357364807*x4 - 13074999.3774/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(0.30624107406733952*x2 - 0.21220319221710396*x3 + 0.89877035125000004);
    double x5 = x0*x4;
    double x6 = 794.74451999999997*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(14474221.614048623*x2 - 10029601.811948752*x3 + 42479609.515960179);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1589.4890399999999*x5)/((x8)*(x8)) + 19823.620807566418*x4*(x0*(0.011324393907113699*T*x11 - 9.0/x13) + 0.0037747979690378993*x11 - 2384.2335599999997*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = -0.21220319221710396*pow(x1, 4.0/3.0) + 0.30624107406733952*x2 + 0.89877035125000004;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 794.74451999999997*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 47264142.010114156*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*x1*x2*(0.14146879481140262*x2 - 0.1020803580224465)*(19823.620807566418*x4*(-2384.2335599999997*x0*x11*x12/((x13)*(x13)) + 0.0037747979690378993*x10/pow(x3, 3.0/2.0) + (0.011324393907113699*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1589.4890399999999*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 0.30624107406733952*x2 - 0.21220319221710396*x3 + 0.89877035125000004;
    double x5 = sqrt(x4);
    double x6 = 2.6491484000000001*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = x7/x8;
    double x10 = 1.0/x5;
    double x11 = x10*x2*(0.33009385455993945*x2 - 0.17013393003741084);
    double x12 = 59470.862422699254*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((0.14146879481140262*x2 - 0.1020803580224465)*(0.14146879481140262*x2 - 0.1020803580224465));
    double x15 = x13*x14;
    double x16 = 157547.14003371386*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 59470.862422699254*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 794.74451999999997*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 47264142.010114156*x15*x20;
    double x27 = exp(x6);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x6);
    double x31 = x10*x30;
    double x32 = -59470.862422699254*x29 + 22449.049069013745*x31;
    double x33 = exp(x22);
    double x34 = x33 - 1;
    double x35 = 1.0/x34;
    double x36 = Debye(x22);
    double x37 = T*x10*x36;
    double x38 = -59470.862422699254*x35 + 74.830163563379159*x37;
    double x39 = 19823.620807566418*x10*x14;

    result += x0*(39224998.132200003*x0 + x11*x32 - x11*x38 - x12*x25 + x12*x9 + x16*x9 + x16*exp(-5.2982968000000001*x5)/((x8)*(x8)) + x18*x32 - x18*x38 - x19*x25 + x19*x9 + 4200642.4312455226*x2 - x25*x26 - 27112057.833851218*x3 - x39*(0.0037747979690378993*T*x17*x36 - 2384.2335599999997*x10*x20*x33/((x34)*(x34)) + x13*(-9.0*x35 + 0.011324393907113699*x37)) + x39*(-7.9474452000000007*x10*x27/((x28)*(x28)) + x13*(-9.0*x29 + 3.397318172134109*x31) + 1.1324393907113697*x17*x30) - x26*exp(-1589.4890399999999*x21)/((x24)*(x24)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(43422664.842145868*x2 - 30088805.435846254*x3 + 127438828.54788055);
    double x5 = 1.0/T;
    double x6 = 0.30624107406733952*x2 - 0.21220319221710396*x3 + 0.89877035125000004;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 794.74451999999997*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1589.4890399999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = -2384.2335599999997*x0*x20*x7;
    double x22 = 0.011324393907113699*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 19823.620807566418*x7;
    double x25 = x15*x6;

    result += x0*(-112688753565.12003*x13*x16 + x13*x4 - 37562917855.040009*x14*x16 + x14*x4 + x24*(0.0037747979690378993*x17 + x21 + x23) - x24*(1894856.5562100909*x20*x25 + x21 + x22 + 3.0*x23 - 3789713.1124201817*x25*exp(x12)/((x19)*(x19)*(x19))) - 75125835710.080017*x16*exp(-2384.2335599999997*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(13372802.415931668*x2 - 9649481.076032415);
    double x5 = -0.21220319221710396*pow(x1, 4.0/3.0) + 0.30624107406733952*x2 + 0.89877035125000004;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 794.74451999999997*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1589.4890399999999*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = pow(T, -3);
    double x15 = 0.14146879481140262*x2 - 0.1020803580224465;
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0037747979690378993*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2384.2335599999997*x3;
    double x24 = 0.011324393907113699*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 19823.620807566418*x15;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = x17*x22;
    double x30 = 1.0/x5;

    result += x0*x1*x2*(112688753565.12003*x12*x16 - x12*x4 + 37562917855.040009*x13*x16 - x13*x4 - x17*x27*(x17*x19 - x22*x23*x6 + x26) + x27*x6*(-x0*(7152.700679999999*x0*x29 - x24*x28 - 3.0*x25*x30) + 1894856.5562100909*x14*x22 - 3789713.1124201817*x14*exp(x11)/((x21)*(x21)*(x21)) + x19*x28 + x23*x29 + x26*x30) + 75125835710.080017*x16*exp(-2384.2335599999997*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(15601602.818586946*x2 - 8041234.2300270125);
    double x4 = -0.21220319221710396*pow(x1, 4.0/3.0) + 0.30624107406733952*x2 + 0.89877035125000004;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 794.74451999999997*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1589.4890399999999*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = pow(T, -2);
    double x14 = 1.0/x5;
    double x15 = 0.14146879481140262*x2 - 0.1020803580224465;
    double x16 = ((x15)*(x15));
    double x17 = x16*x2;
    double x18 = x14*x17;
    double x19 = x13*x18;
    double x20 = pow(x4, -3.0/2.0);
    double x21 = Debye(x7);
    double x22 = 0.0037747979690378993*T*x21;
    double x23 = x20*x22;
    double x24 = exp(x7);
    double x25 = x24 - 1;
    double x26 = x24/((x25)*(x25));
    double x27 = 2384.2335599999997*x26;
    double x28 = x0*x14*x27;
    double x29 = 1.0/x4;
    double x30 = 1.0/x25;
    double x31 = T*x14*x21;
    double x32 = x29*(-9.0*x30 + 0.011324393907113699*x31);
    double x33 = 0.33009385455993945*x2 - 0.17013393003741084;
    double x34 = x15*x2;
    double x35 = x17*x29;
    double x36 = x13*x35;
    double x37 = x0*x2;
    double x38 = -9.0*x30 + 0.011324393907113699*x31;
    double x39 = x38/((x4)*(x4));

    result += x37*(-112688753565.12003*x11*x19 - x11*x3 - 37562917855.040009*x12*x19 - x12*x3 + 19823.620807566418*x18*(x23 - x28 + x32) - 75125835710.080017*x19*exp(-2384.2335599999997*x6)/((x9)*(x9)*(x9)) - 19823.620807566418*x5*(-x16*x20*x27*x37 + x17*x39 + x22*x34*(0.42440638443420786*x2 - 0.30624107406733947)/pow(x4, 5.0/2.0) + x23*x33 + 1894856.5562100909*x26*x36 - x28*x33 + x29*x33*x38 + x34*x39*(0.28293758962280524*x2 - 0.204160716044893) - 3.0*x35*(-x23 + x28 - x32) - 3789713.1124201817*x36*exp(x10)/((x25)*(x25)*(x25))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 0.30624107406733952*x2 - 0.21220319221710396*x3 + 0.89877035125000004;
    double x5 = sqrt(x4);
    double x6 = 2.6491484000000001*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 0.14146879481140262*x2 - 0.1020803580224465;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 5.2982968000000001*x5;
    double x16 = exp(-x15)/((x8)*(x8));
    double x17 = 157547.14003371386*x16;
    double x18 = x7/x8;
    double x19 = 157547.14003371386*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(1.1003128485331315*x2 - 0.45369048009976221);
    double x22 = 59470.862422699254*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 794.74451999999997*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = pow(T, -2);
    double x30 = x12*x29;
    double x31 = 1589.4890399999999*x24;
    double x32 = exp(-x31)/((x27)*(x27));
    double x33 = 47264142.010114156*x23;
    double x34 = x14*x33;
    double x35 = x13*(0.28293758962280524*x2 - 0.204160716044893);
    double x36 = ((x10)*(x10));
    double x37 = x0*x36;
    double x38 = x35*x37;
    double x39 = (0.42440638443420786*x2 - 0.30624107406733947)/pow(x4, 5.0/2.0);
    double x40 = x37*x39;
    double x41 = 59470.862422699254*x40;
    double x42 = 1.0/x4;
    double x43 = 0.33009385455993945*x2 - 0.17013393003741084;
    double x44 = x42*x43;
    double x45 = x10*x3;
    double x46 = x17*x45;
    double x47 = x42*(0.6601877091198789*x2 - 0.34026786007482168);
    double x48 = x43*x45;
    double x49 = x48*x9;
    double x50 = 178412.58726809776*x49;
    double x51 = x19*x45;
    double x52 = x33*x38;
    double x53 = x33*x45;
    double x54 = x32*x53;
    double x55 = x28*x53;
    double x56 = exp(x6);
    double x57 = x56 - 1;
    double x58 = 1.0/x57;
    double x59 = Debye(x6);
    double x60 = x20*x59;
    double x61 = -3*x58 + 1.1324393907113697*x60;
    double x62 = 19823.620807566418*x20;
    double x63 = x21*x62;
    double x64 = exp(x25);
    double x65 = x64 - 1;
    double x66 = 1.0/x65;
    double x67 = Debye(x25);
    double x68 = T*x20*x67;
    double x69 = -3*x66 + 0.0037747979690378993*x68;
    double x70 = 19823.620807566418*x40;
    double x71 = 59470.862422699254*x49;
    double x72 = 1.1324393907113697*x59;
    double x73 = x72*x9;
    double x74 = x56/((x57)*(x57));
    double x75 = 7.9474452000000007*x74;
    double x76 = x20*x75;
    double x77 = x42*(-9.0*x58 + 3.397318172134109*x60);
    double x78 = 39647.241615132836*x73 - 39647.241615132836*x76 + 39647.241615132836*x77;
    double x79 = x20*x48;
    double x80 = 0.0037747979690378993*T*x67;
    double x81 = x80*x9;
    double x82 = x64/((x65)*(x65));
    double x83 = 2384.2335599999997*x23*x82;
    double x84 = x20*x83;
    double x85 = x42*(-9.0*x66 + 0.011324393907113699*x68);
    double x86 = 39647.241615132836*x81 - 39647.241615132836*x84 + 39647.241615132836*x85;
    double x87 = x10*x2;
    double x88 = x39*x87;
    double x89 = x2*x36;
    double x90 = x42*x89;
    double x91 = x89*x9;
    double x92 = 3.0*x61;
    double x93 = x13*x89;
    double x94 = x35*x87;
    double x95 = 3.0*x90;
    double x96 = x45*x62;
    double x97 = x29*x90;
    double x98 = 3.0*x69;

    result += (-156899992.52880001*x0 - 1252097.2618346671*x12*x16 - 417365.75394488906*x12*x18 - x12*x78 + x12*x86 - 834731.50788977812*x12*exp(-7.9474452000000007*x5)/((x8)*(x8)*(x8)) - x14*x17 - x14*x19 - x17*x38 - x18*x22 - x18*x41 - x18*x50 - x19*x38 - 11201713.149988059*x2 + x22*x28 + 37562917855.040009*x28*x30 + x28*x34 + x28*x41 + x28*x50 + x28*x52 + 90373526.112837389*x3 + 112688753565.12003*x30*x32 + x32*x34 + x32*x52 - x44*x46 - x44*x51 + x44*x54 + x44*x55 - x46*x47 - x47*x51 + x47*x54 + x47*x55 - x61*x63 - x61*x70 - x61*x71 + x63*x69 + x69*x70 + x69*x71 - x78*x79 + x79*x86 - x96*(x43*x73 - x43*x76 + x44*x92 + x72*x88 + 21.053961735667681*x74*x90 - x75*x91 + x92*x93 + x92*x94 - x95*(-x73 + x76 - x77) - 42.107923471335361*x90*exp(x15)/((x57)*(x57)*(x57))) + x96*(x43*x81 - x43*x84 + x44*x98 + x80*x88 + 1894856.5562100909*x82*x97 - x83*x91 + x93*x98 + x94*x98 - x95*(-x81 + x84 - x85) - 3789713.1124201817*x97*exp(x31)/((x65)*(x65)*(x65))) + 75125835710.080017*x30*exp(-2384.2335599999997*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.724;
static const double Vmax = 1.15*2.724;
static double V = 0.9*2.724;

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



const char *Cristobalite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Cristobalite_slb_em_coder_calib_name(void) {
    return "Cristobalite_slb_em";
}

const char *Cristobalite_slb_em_coder_calib_formula(void) {
    return "SiO2";
}

const double Cristobalite_slb_em_coder_calib_mw(void) {
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

const double *Cristobalite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Cristobalite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Cristobalite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Cristobalite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Cristobalite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Cristobalite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Cristobalite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Cristobalite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Cristobalite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Cristobalite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Cristobalite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Cristobalite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Cristobalite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Cristobalite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Cristobalite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Cristobalite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Cristobalite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Cristobalite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Cristobalite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Cristobalite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Cristobalite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Cristobalite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Cristobalite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Cristobalite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Cristobalite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Cristobalite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Cristobalite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Cristobalite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Cristobalite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

