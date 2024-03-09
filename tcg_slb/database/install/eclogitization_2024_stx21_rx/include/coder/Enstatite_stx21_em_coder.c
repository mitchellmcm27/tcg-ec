
static char *identifier = "Enstatite_stx21_em.emml:0df32a35d71f2b738039b3d26adf3798e0070252:Sat Mar  9 00:07:23 2024";



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
    double x3 = sqrt(38.44260390632742*x1 - 51.737659551353467*x2 - 5.831495441225);
    double x4 = 2.7073742333333333*x3;
    double x5 = 812.21226999999999*x3/T;

    result += 249.43387854459721*T*log(1 - exp(-x5)) - 83.144626181532402*T*Debye(x5) + 65222887.505612053*x1 - 308965040.49000734*x2 - 74830.163563379159*log(1 - exp(-x4)) + 24943.387854459721*Debye(x4) - 6793034.0978577528 + 448959693.73367566/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-51.737659551353467*pow(x0, 4.0/3.0) + 38.44260390632742*pow(x0, 2.0/3.0) - 5.831495441225);
    double x2 = x1/T;
    double x3 = 812.21226999999999*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -202593.25670761158*x2*x5/x6 + 67531.08556920386*x2*(-0.0036936157095976891*T*x4/x1 + 3/(exp(x3) - 1)) - 83.144626181532402*x4 + 249.43387854459721*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(38.44260390632742*x1 - 51.737659551353467*x3 - 5.831495441225);
    double x6 = 2.7073742333333333*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-12.81420130210914*x2 + 34.491773034235642*x4);
    double x10 = 202593.25670761158*x9;
    double x11 = 812.21226999999999*x5/T;
    double x12 = exp(-x11);
    double x13 = 67531.08556920386*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.108084712879307*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0036936157095976891*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 43481925.003741369*x2 + 411953387.32000977*x4 - 897919387.46735132/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(38.44260390632742*x2 - 51.737659551353467*x3 - 5.831495441225);
    double x5 = x0*x4;
    double x6 = 812.21226999999999*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-6325681609.0528698*x2 + 8513366116.3251104*x3 + 959565162.53991485);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1624.42454*x5)/((x8)*(x8)) - 67531.08556920386*x4*(x0*(0.011080847128793068*T*x11 - 9.0/x13) + 0.0036936157095976891*x11 - 2436.63681*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 38.44260390632742*x2;
    double x4 = 51.737659551353467*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 5.831495441225;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 812.21226999999999*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 164548728.91718194*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(34.491773034235642*x2 - 12.81420130210914)*(-67531.08556920386*x6*(2436.63681*x0*x13*x14/((x15)*(x15)) - 0.0036936157095976891*x12/pow(x5, 3.0/2.0) + (0.011080847128793068*x12*x13 - 9.0/x15)/(-x3 + x4 + 5.831495441225)) + x11*x9/x10 + x11*exp(-1624.42454*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 38.44260390632742*x2;
    double x5 = 51.737659551353467*x3;
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = 2.7073742333333333*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x12*x2*(80.480803746549839*x2 - 21.3570021701819);
    double x14 = 202593.25670761158*x13;
    double x15 = 1.0/(-x4 + x5 + 5.831495441225);
    double x16 = x3*((34.491773034235642*x2 - 12.81420130210914)*(34.491773034235642*x2 - 12.81420130210914));
    double x17 = x15*x16;
    double x18 = 548495.76305727311*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 202593.25670761158*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 812.21226999999999*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 164548728.91718194*x17*x22;
    double x29 = exp(x8);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x8);
    double x33 = x12*x32;
    double x34 = -202593.25670761158*x31 + 74830.163563379174*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x12*x38;
    double x40 = -202593.25670761158*x37 + 249.43387854459718*x39;
    double x41 = 67531.08556920386*x12*x16;

    result += x0*(2693758162.4020538*x0 + x11*x14 - x11*x18 + x11*x21 + x13*x34 - x13*x40 - x14*x27 + 72469875.006235614*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 - 961224570.41335607*x3 - x41*(8.1221227000000003*x12*x29/((x30)*(x30)) + x15*(-9.0000000000000018*x31 + 3.3242541386379214*x33) - 1.108084712879307*x19*x32) + x41*(2436.63681*x12*x22*x35/((x36)*(x36)) + x15*(-9.0*x37 + 0.011080847128793068*x39) - 0.0036936157095976891*x19*x38) + x28*exp(-1624.42454*x23)/((x26)*(x26)) - x18*exp(-5.4147484666666665*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-18977044827.158607*x2 + 25540098348.97533*x3 + 2878695487.6197443);
    double x5 = 1.0/T;
    double x6 = 38.44260390632742*x2;
    double x7 = 51.737659551353467*x3;
    double x8 = x6 - x7 - 5.831495441225;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 812.21226999999999*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1624.42454*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2436.63681*x0*x22*x9;
    double x24 = 0.011080847128793068*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 67531.08556920386*x9;
    double x27 = x17*(-x6 + x7 + 5.831495441225);

    result += x0*(-400945489918.31696*x15*x18 - x15*x4 - 133648496639.43898*x16*x18 - x16*x4 + x26*(0.0036936157095976891*x19 + x23 + x25) - x26*(-1979066.3146156587*x22*x27 + x23 + x24 + 3.0*x25 + 3958132.6292313174*x27*exp(x14)/((x21)*(x21)*(x21))) - 267296993278.87796*x18*exp(-2436.63681*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(11351154821.766813*x2 - 4217121072.7019134);
    double x5 = 38.44260390632742*x2;
    double x6 = 51.737659551353467*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 5.831495441225;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 812.21226999999999*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1624.42454*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 34.491773034235642*x2 - 12.81420130210914;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0036936157095976891*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2436.63681*x3;
    double x26 = 0.011080847128793068*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 67531.08556920386*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 5.831495441225);

    result += x0*x1*x2*(400945489918.31696*x14*x18 - x14*x4 + 133648496639.43898*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7309.9104299999999*x0*x31 - x26*x30 + 3.0*x27*x32) + 1979066.3146156587*x16*x24 - 3958132.6292313174*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 267296993278.87796*x18*exp(-2436.63681*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(13243013958.727949*x2 - 3514267560.584928);
    double x4 = 38.44260390632742*x2;
    double x5 = 51.737659551353467*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 812.21226999999999*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1624.42454*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 34.491773034235642*x2 - 12.81420130210914;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0036936157095976891*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2436.63681*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 5.831495441225;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.011080847128793068*x33));
    double x35 = 80.480803746549839*x2 - 21.3570021701819;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0036936157095976891*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(400945489918.31696*x13*x20 + x13*x3 + 133648496639.43898*x14*x20 + x14*x3 + 67531.08556920386*x15*x34 + 67531.08556920386*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(103.47531910270692*x2 - 38.44260390632742)/pow(x6, 5.0/2.0) + x24*x35 - 1979066.3146156587*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(68.983546068471284*x2 - 25.62840260421828) + 3958132.6292313174*x37*exp(x12)/((x26)*(x26)*(x26))) + 267296993278.87796*x20*exp(-2436.63681*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 38.44260390632742*x2;
    double x5 = 51.737659551353467*x3;
    double x6 = x4 - x5 - 5.831495441225;
    double x7 = sqrt(x6);
    double x8 = 2.7073742333333333*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 34.491773034235642*x2 - 12.81420130210914;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.4147484666666665*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 5.831495441225;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 548495.76305727311*x19;
    double x21 = x9/x10;
    double x22 = 1.0/x7;
    double x23 = x2*(268.26934582183276*x2 - 56.952005787151734);
    double x24 = 202593.25670761158*x22*x23;
    double x25 = 1.0/T;
    double x26 = x25*x7;
    double x27 = 812.21226999999999*x26;
    double x28 = exp(-x27);
    double x29 = -x28 + 1;
    double x30 = x28/x29;
    double x31 = pow(T, -2);
    double x32 = x14*x31;
    double x33 = 1624.42454*x26;
    double x34 = exp(-x33)/((x29)*(x29));
    double x35 = 164548728.91718194*x25;
    double x36 = x19*x35;
    double x37 = x18*(68.983546068471284*x2 - 25.62840260421828);
    double x38 = ((x12)*(x12));
    double x39 = x0*x38;
    double x40 = x37*x39;
    double x41 = 548495.76305727311*x40;
    double x42 = (103.47531910270692*x2 - 38.44260390632742)/pow(x6, 5.0/2.0);
    double x43 = x39*x42;
    double x44 = 202593.25670761158*x43;
    double x45 = 1.0/x17;
    double x46 = 80.480803746549839*x2 - 21.3570021701819;
    double x47 = x45*x46;
    double x48 = x12*x3;
    double x49 = 548495.76305727311*x48;
    double x50 = x16*x49;
    double x51 = x45*(160.96160749309968*x2 - 42.7140043403638);
    double x52 = x46*x48;
    double x53 = x11*x52;
    double x54 = 607779.77012283471*x53;
    double x55 = x21*x49;
    double x56 = x35*x40;
    double x57 = x35*x48;
    double x58 = x34*x57;
    double x59 = x30*x57;
    double x60 = exp(x8);
    double x61 = x60 - 1;
    double x62 = 1.0/x61;
    double x63 = Debye(x8);
    double x64 = x22*x63;
    double x65 = -3*x62 + 1.108084712879307*x64;
    double x66 = 67531.08556920386*x22;
    double x67 = x23*x66;
    double x68 = exp(x27);
    double x69 = x68 - 1;
    double x70 = 1.0/x69;
    double x71 = Debye(x27);
    double x72 = T*x22*x71;
    double x73 = -3*x70 + 0.0036936157095976891*x72;
    double x74 = 67531.08556920386*x43;
    double x75 = 202593.25670761158*x53;
    double x76 = 1.108084712879307*x63;
    double x77 = x11*x76;
    double x78 = x60/((x61)*(x61));
    double x79 = 8.1221227000000003*x78;
    double x80 = x22*x79;
    double x81 = x45*(-9.0000000000000018*x62 + 3.3242541386379214*x64) - x77 + x80;
    double x82 = 135062.17113840772*x81;
    double x83 = x22*x52;
    double x84 = 0.0036936157095976891*T*x71;
    double x85 = x11*x84;
    double x86 = x68/((x69)*(x69));
    double x87 = 2436.63681*x25*x86;
    double x88 = x22*x87;
    double x89 = x45*(-9.0*x70 + 0.011080847128793068*x72) - x85 + x88;
    double x90 = 135062.17113840772*x89;
    double x91 = x12*x2;
    double x92 = x42*x91;
    double x93 = x2*x38;
    double x94 = x45*x93;
    double x95 = x11*x93;
    double x96 = 3.0000000000000004*x65;
    double x97 = x18*x93;
    double x98 = x37*x91;
    double x99 = x48*x66;
    double x100 = x31*x94;
    double x101 = 3.0*x73;

    result += (-10775032649.608215*x0 - 4454949.8879812993*x14*x16 - 1484983.2959937665*x14*x21 + x14*x82 - x14*x90 - x16*x20 - x16*x41 - 193253000.0166283*x2 - x20*x21 - x21*x24 - x21*x41 - x21*x44 - x21*x54 + x24*x30 + 3204081901.3778534*x3 + 133648496639.43898*x30*x32 + x30*x36 + x30*x44 + x30*x54 + x30*x56 + 400945489918.31696*x32*x34 + x34*x36 + x34*x56 + x47*x50 + x47*x55 - x47*x58 - x47*x59 + x50*x51 + x51*x55 - x51*x58 - x51*x59 - x65*x67 - x65*x74 - x65*x75 + x67*x73 + x73*x74 + x73*x75 + x82*x83 - x83*x90 + x99*(-1979066.3146156587*x100*x86 + 3958132.6292313174*x100*exp(x33)/((x69)*(x69)*(x69)) - x101*x47 + x101*x97 + x101*x98 + x46*x85 - x46*x88 + x84*x92 - x87*x95 + 3.0*x89*x94) - x99*(x46*x77 - x46*x80 - x47*x96 + x76*x92 - 21.989625717951764*x78*x94 - x79*x95 + 3.0000000000000004*x81*x94 + x96*x97 + x96*x98 + 43.979251435903528*x94*exp(x15)/((x61)*(x61)*(x61))) + 267296993278.87796*x32*exp(-2436.63681*x26)/((x29)*(x29)*(x29)) - 2969966.591987533*x14*exp(-8.1221227000000003*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*6.2676;
static const double Vmax = 1.15*6.2676;
static double V = 0.9*6.2676;

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



const char *Enstatite_stx21_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Enstatite_stx21_em_coder_calib_name(void) {
    return "Enstatite_stx21_em";
}

const char *Enstatite_stx21_em_coder_calib_formula(void) {
    return "Mg2Si2O6";
}

const double Enstatite_stx21_em_coder_calib_mw(void) {
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

const double *Enstatite_stx21_em_coder_calib_elements(void) {
    return elmformula;
}

double Enstatite_stx21_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Enstatite_stx21_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Enstatite_stx21_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Enstatite_stx21_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Enstatite_stx21_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Enstatite_stx21_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Enstatite_stx21_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Enstatite_stx21_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Enstatite_stx21_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Enstatite_stx21_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Enstatite_stx21_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Enstatite_stx21_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Enstatite_stx21_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Enstatite_stx21_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Enstatite_stx21_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Enstatite_stx21_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Enstatite_stx21_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Enstatite_stx21_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Enstatite_stx21_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Enstatite_stx21_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Enstatite_stx21_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Enstatite_stx21_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Enstatite_stx21_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Enstatite_stx21_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Enstatite_stx21_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Enstatite_stx21_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Enstatite_stx21_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Enstatite_stx21_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

