
static char *identifier = "Pyrope_slb_em.emml:0c453e677b5c08d3970bb1d2ab48c1ab9bbace1b:Tue Aug 22 19:34:29 2023";



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
    double x3 = sqrt(16.706570806133879*x1 - 3.4703776090950673*x2 - 2.1793895337499998);
    double x4 = 2.7441261000000003*x3;
    double x5 = 823.23783000000003*x3/T;

    result += 498.86775708919441*T*log(1 - exp(-x5)) - 166.2892523630648*T*Debye(x5) - 200103351.06996444*x1 + 458434043.66093028*x2 - 149660.32712675832*log(1 - exp(-x4)) + 49886.775708919442*Debye(x4) + 14521506.376463132 + 153239430.1863752/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-3.4703776090950673*pow(x0, 4.0/3.0) + 16.706570806133879*pow(x0, 2.0/3.0) - 2.1793895337499998);
    double x2 = x1/T;
    double x3 = 823.23783000000003*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -410686.80980307556*x2*x5/x6 + 136895.60326769185*x2*(-0.0036441474027013553*T*x4/x1 + 3/(exp(x3) - 1)) - 166.2892523630648*x4 + 498.86775708919441*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(16.706570806133879*x1 - 3.4703776090950673*x3 - 2.1793895337499998);
    double x6 = 2.7441261000000003*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-5.5688569353779593*x2 + 2.3135850727300449*x4);
    double x10 = 410686.80980307556*x9;
    double x11 = 823.23783000000003*x5/T;
    double x12 = exp(-x11);
    double x13 = 136895.60326769185*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.0932442208104065*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0036441474027013553*T*x8*Debye(x11) + 3/(exp(x11) - 1)) + 133402234.04664296*x2 - 611245391.54790699*x4 - 306478860.3727504/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(16.706570806133879*x2 - 3.4703776090950673*x3 - 2.1793895337499998);
    double x5 = x0*x4;
    double x6 = 823.23783000000003*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-5648373275.4889917*x2 + 1173310092.8091729*x3 + 736836167.1680851);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1646.4756600000001*x5)/((x8)*(x8)) - 136895.60326769185*x4*(x0*(0.010932442208104066*T*x11 - 9.0/x13) + 0.0036441474027013553*x11 - 2469.7134900000001*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 16.706570806133879*x2;
    double x4 = 3.4703776090950673*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 2.1793895337499998;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 823.23783000000003*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 338092918.11190665*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(2.3135850727300449*x2 - 5.5688569353779593)*(-136895.60326769185*x6*(2469.7134900000001*x0*x13*x14/((x15)*(x15)) - 0.0036441474027013553*x12/pow(x5, 3.0/2.0) + (0.010932442208104066*x12*x13 - 9.0/x15)/(-x3 + x4 + 2.1793895337499998)) + x11*x9/x10 + x11*exp(-1646.4756600000001*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 16.706570806133879*x2;
    double x5 = 3.4703776090950673*x3;
    double x6 = x4 - x5 - 2.1793895337499998;
    double x7 = sqrt(x6);
    double x8 = 2.7441261000000003*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x12*x2*(5.398365169703438*x2 - 9.2814282256299325);
    double x14 = 410686.80980307556*x13;
    double x15 = 1.0/(-x4 + x5 + 2.1793895337499998);
    double x16 = x3*((2.3135850727300449*x2 - 5.5688569353779593)*(2.3135850727300449*x2 - 5.5688569353779593));
    double x17 = x15*x16;
    double x18 = 1126976.3937063557*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 410686.80980307556*x20;
    double x22 = 1.0/T;
    double x23 = x22*x7;
    double x24 = 823.23783000000003*x23;
    double x25 = exp(-x24);
    double x26 = -x25 + 1;
    double x27 = x25/x26;
    double x28 = 338092918.11190665*x17*x22;
    double x29 = exp(x8);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x8);
    double x33 = x12*x32;
    double x34 = -410686.80980307556*x31 + 149660.32712675832*x33;
    double x35 = exp(x24);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x24);
    double x39 = x12*x38;
    double x40 = -410686.80980307556*x37 + 498.86775708919441*x39;
    double x41 = 136895.60326769185*x12*x16;

    result += x0*(919436581.1182512*x0 + x11*x14 - x11*x18 + x11*x21 + x13*x34 - x13*x40 - x14*x27 - 222337056.74440491*x2 + x20*x34 - x20*x40 - x21*x27 + x27*x28 + 1426239246.9451163*x3 - x41*(8.2323783000000006*x12*x29/((x30)*(x30)) + x15*(-9.0*x31 + 3.2797326624312193*x33) - 1.0932442208104065*x19*x32) + x41*(2469.7134900000001*x12*x22*x35/((x36)*(x36)) + x15*(-9.0*x37 + 0.010932442208104066*x39) - 0.0036441474027013553*x19*x38) + x28*exp(-1646.4756600000001*x23)/((x26)*(x26)) - x18*exp(-5.4882522000000007*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-16945119826.466976*x2 + 3519930278.4275188*x3 + 2210508501.5042553);
    double x5 = 1.0/T;
    double x6 = 16.706570806133879*x2;
    double x7 = 3.4703776090950673*x3;
    double x8 = x6 - x7 - 2.1793895337499998;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 823.23783000000003*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1646.4756600000001*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = -2469.7134900000001*x0*x22*x9;
    double x24 = 0.010932442208104066*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 136895.60326769185*x9;
    double x27 = x17*(-x6 + x7 + 2.1793895337499998);

    result += x0*(-834992640734.44116*x15*x18 - x15*x4 - 278330880244.81372*x16*x18 - x16*x4 + x26*(0.0036441474027013553*x19 + x23 + x25) - x26*(-2033161.5742293268*x22*x27 + x23 + x24 + 3.0*x25 + 4066323.1484586536*x27*exp(x14)/((x21)*(x21)*(x21))) - 556661760489.62744*x18*exp(-2469.7134900000001*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(1564413457.0788972*x2 - 3765582183.6593275);
    double x5 = 16.706570806133879*x2;
    double x6 = 3.4703776090950673*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 2.1793895337499998;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 823.23783000000003*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1646.4756600000001*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 2.3135850727300449*x2 - 5.5688569353779593;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0036441474027013553*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2469.7134900000001*x3;
    double x26 = 0.010932442208104066*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 136895.60326769185*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 2.1793895337499998);

    result += x0*x1*x2*(834992640734.44116*x14*x18 - x14*x4 + 278330880244.81372*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7409.1404700000003*x0*x31 - x26*x30 + 3.0*x27*x32) + 2033161.5742293268*x16*x24 - 4066323.1484586536*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 556661760489.62744*x18*exp(-2469.7134900000001*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(1825149033.2587135*x2 - 3137985153.0494399);
    double x4 = 16.706570806133879*x2;
    double x5 = 3.4703776090950673*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 2.1793895337499998;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 823.23783000000003*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1646.4756600000001*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 2.3135850727300449*x2 - 5.5688569353779593;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0036441474027013553*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2469.7134900000001*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 2.1793895337499998;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.010932442208104066*x33));
    double x35 = 5.398365169703438*x2 - 9.2814282256299325;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0036441474027013553*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(834992640734.44116*x13*x20 + x13*x3 + 278330880244.81372*x14*x20 + x14*x3 + 136895.60326769185*x15*x34 + 136895.60326769185*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(6.9407552181901346*x2 - 16.706570806133879)/pow(x6, 5.0/2.0) + x24*x35 - 2033161.5742293268*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(4.6271701454600898*x2 - 11.137713870755919) + 4066323.1484586536*x37*exp(x12)/((x26)*(x26)*(x26))) + 556661760489.62744*x20*exp(-2469.7134900000001*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 16.706570806133879*x2;
    double x5 = 3.4703776090950673*x3;
    double x6 = x4 - x5 - 2.1793895337499998;
    double x7 = sqrt(x6);
    double x8 = 2.7441261000000003*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 2.3135850727300449*x2 - 5.5688569353779593;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.4882522000000007*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 2.1793895337499998;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1126976.3937063557*x16;
    double x21 = x9/x10;
    double x22 = 1126976.3937063557*x21;
    double x23 = 1.0/x7;
    double x24 = 17.994550565678125*x2 - 24.750475268346484;
    double x25 = x2*x23*x24;
    double x26 = 410686.80980307556*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 823.23783000000003*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 410686.80980307556*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1646.4756600000001*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 338092918.11190665*x27;
    double x39 = x19*x38;
    double x40 = 4.6271701454600898*x2 - 11.137713870755919;
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x18*x40*x42;
    double x44 = (6.9407552181901346*x2 - 16.706570806133879)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 5.398365169703438*x2 - 9.2814282256299325;
    double x47 = x12*x3;
    double x48 = x46*x47;
    double x49 = 1.0/x17;
    double x50 = x20*x49;
    double x51 = x47*(10.796730339406876*x2 - 18.562856451259865);
    double x52 = 1232060.4294092266*x11;
    double x53 = x22*x49;
    double x54 = x32*x48;
    double x55 = x38*x43;
    double x56 = x38*x49;
    double x57 = x51*x56;
    double x58 = 136895.60326769185*x23;
    double x59 = x24*x58;
    double x60 = exp(x8);
    double x61 = x60 - 1;
    double x62 = 1.0/x61;
    double x63 = Debye(x8);
    double x64 = x23*x63;
    double x65 = -3*x62 + 1.0932442208104065*x64;
    double x66 = x2*x65;
    double x67 = exp(x29);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x29);
    double x71 = T*x23*x70;
    double x72 = -3*x69 + 0.0036441474027013553*x71;
    double x73 = x2*x72;
    double x74 = 136895.60326769185*x45;
    double x75 = x46*x65;
    double x76 = 410686.80980307556*x11*x47;
    double x77 = x46*x72;
    double x78 = 1.0932442208104065*x63;
    double x79 = x11*x78;
    double x80 = x60/((x61)*(x61));
    double x81 = 8.2323783000000006*x80;
    double x82 = x23*x81;
    double x83 = x49*(-9.0*x62 + 3.2797326624312193*x64) - x79 + x82;
    double x84 = 273791.20653538371*x83;
    double x85 = x23*x48;
    double x86 = 0.0036441474027013553*T*x70;
    double x87 = x11*x86;
    double x88 = x67/((x68)*(x68));
    double x89 = 2469.7134900000001*x27*x88;
    double x90 = x23*x89;
    double x91 = x49*(-9.0*x69 + 0.010932442208104066*x71) - x87 + x90;
    double x92 = 273791.20653538371*x91;
    double x93 = x12*x2*x44;
    double x94 = x2*x41;
    double x95 = x49*x94;
    double x96 = x11*x94;
    double x97 = 3.0*x49;
    double x98 = 3.0*x18;
    double x99 = x12*x40;
    double x100 = x94*x97;
    double x101 = x47*x58;
    double x102 = x34*x95;
    double x103 = x73*x98;

    result += (-3677746324.4730048*x0 - x101*(x100*x83 + x46*x79 - x46*x82 + x65*x94*x98 + x66*x98*x99 - x75*x97 + x78*x93 - 22.590684158103635*x80*x95 - x81*x96 + 45.181368316207269*x95*exp(x15)/((x61)*(x61)*(x61))) + x101*(x100*x91 - 2033161.5742293268*x102*x88 + 4066323.1484586536*x102*exp(x36)/((x68)*(x68)*(x68)) + x103*x41 + x103*x99 + x46*x87 - x46*x90 - x77*x97 + x86*x93 - x89*x96) - 9277696.0081604607*x14*x16 - 3092565.3360534869*x14*x21 + x14*x84 - x14*x92 - x19*x20 - x19*x22 + 592898817.98507977*x2 - x20*x43 - x21*x48*x52 - x22*x43 - x25*x26 + x25*x33 - x26*x45 - 4754130823.1503878*x3 + 278330880244.81372*x32*x35 + x32*x39 + x32*x55 - x32*x57 + x33*x45 + 834992640734.44116*x35*x37 + x37*x39 - x37*x48*x56 + x37*x55 - x37*x57 + x48*x50 + x48*x53 + x50*x51 + x51*x53 + x52*x54 - x54*x56 - x59*x66 + x59*x73 - x65*x74 + x72*x74 - x75*x76 + x76*x77 + x84*x85 - x85*x92 + 556661760489.62744*x35*exp(-2469.7134900000001*x28)/((x31)*(x31)*(x31)) - 6185130.6721069738*x14*exp(-8.2323783000000006*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*11.308;
static const double Vmax = 1.15*11.308;
static double V = 0.9*11.308;

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



const char *Pyrope_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Pyrope_slb_em_coder_calib_name(void) {
    return "Pyrope_slb_em";
}

const char *Pyrope_slb_em_coder_calib_formula(void) {
    return "Mg3Al2Si3O12";
}

const double Pyrope_slb_em_coder_calib_mw(void) {
    return 403.12737999999996;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,12.0,0.0,0.0,0.0,
        3.0,2.0,3.0,0.0,0.0,0.0,
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

const double *Pyrope_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Pyrope_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Pyrope_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Pyrope_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Pyrope_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Pyrope_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Pyrope_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Pyrope_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Pyrope_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Pyrope_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Pyrope_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Pyrope_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Pyrope_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Pyrope_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Pyrope_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Pyrope_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Pyrope_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Pyrope_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Pyrope_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Pyrope_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Pyrope_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Pyrope_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Pyrope_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Pyrope_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Pyrope_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Pyrope_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Pyrope_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Pyrope_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Pyrope_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Pyrope_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Pyrope_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

