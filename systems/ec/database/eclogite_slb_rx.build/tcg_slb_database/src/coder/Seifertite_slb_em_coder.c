
static char *identifier = "Seifertite_slb_em.emml:7130c723ab82602f9393adfe0e833b5489739b9d:Thu Feb 10 16:53:48 2022";



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
    double x3 = sqrt(10.813310026713001*x1 - 3.5311699741784297*x2 - 5.4515027622500014);
    double x4 = 3.8025733333333331*x3;
    double x5 = 1140.7719999999999*x3/T;

    result += 74.830163563379159*T*log(1 - exp(-x5)) - 24.943387854459719*T*Debye(x5) - 12265868.946311224*x1 + 7465041.276672163*x2 - 22449.049069013745*log(1 - exp(-x4)) + 7483.0163563379156*Debye(x4) + 4204382.8489036011 + 73100.880076024972/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-3.5311699741784297*pow(x0, 4.0/3.0) + 10.813310026713001*pow(x0, 2.0/3.0) - 5.4515027622500014);
    double x2 = x1/T;
    double x3 = 1140.7719999999999*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -85364.155348523171*x2*x5/x6 + 28454.71844950772*x2*(-0.0026297980665724616*T*x4/x1 + 3/(exp(x3) - 1)) - 24.943387854459719*x4 + 74.830163563379159*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(10.813310026713001*x1 - 3.5311699741784297*x3 - 5.4515027622500014);
    double x6 = 3.8025733333333331*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-3.6044366755710002*x2 + 2.3541133161189531*x4);
    double x10 = 1140.7719999999999*x5/T;
    double x11 = exp(-x10);
    double x12 = 28454.71844950772*x9;

    result += 85364.155348523171*x11*x9/(-x11 + 1) + x12*(-0.78893941997173855*x8*Debye(x6) + 3/(exp(x6) - 1)) - x12*(-0.0026297980665724616*T*x8*Debye(x10) + 3/(exp(x10) - 1)) + 8177245.9642074825*x2 - 9953388.3688962162*x4 - 85364.155348523156*x7*x9/(-x7 + 1) - 146201.76015204994/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(10.813310026713001*x2 - 3.5311699741784297*x3 - 5.4515027622500014);
    double x5 = x0*x4;
    double x6 = 1140.7719999999999*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1053011357.0527689*x2 + 343868998.23530877*x3 + 530872998.87569869);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-2281.5439999999999*x5)/((x8)*(x8)) - 28454.71844950772*x4*(x0*(0.0078893941997173842*T*x11 - 9.0/x13) + 0.0026297980665724616*x11 - 3422.3159999999998*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 10.813310026713001*x2;
    double x4 = 3.5311699741784297*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 5.4515027622500014;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 1140.7719999999999*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 97381038.225245476*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(2.3541133161189531*x2 - 3.6044366755710002)*(-28454.71844950772*x6*(3422.3159999999998*x0*x13*x14/((x15)*(x15)) - 0.0026297980665724616*x12/pow(x5, 3.0/2.0) + (0.0078893941997173842*x12*x13 - 9.0/x15)/(-x3 + x4 + 5.4515027622500014)) + x11*x9/x10 + x11*exp(-2281.5439999999999*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 10.813310026713001*x2;
    double x5 = 3.5311699741784297*x3;
    double x6 = x4 - x5 - 5.4515027622500014;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(5.4929310709442234*x2 - 6.0073944592850008);
    double x10 = 3.8025733333333331*x7;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = x11/x12;
    double x14 = 85364.155348523156*x13;
    double x15 = 1.0/(-x4 + x5 + 5.4515027622500014);
    double x16 = x3*((2.3541133161189531*x2 - 3.6044366755710002)*(2.3541133161189531*x2 - 3.6044366755710002));
    double x17 = x15*x16;
    double x18 = 324603.46075081819*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 1.0/T;
    double x22 = x21*x7;
    double x23 = 1140.7719999999999*x22;
    double x24 = exp(-x23);
    double x25 = -x24 + 1;
    double x26 = x24/x25;
    double x27 = 85364.155348523171*x26;
    double x28 = 97381038.225245476*x17*x21;
    double x29 = exp(x10);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x10);
    double x33 = x32*x8;
    double x34 = -85364.155348523156*x31 + 22449.049069013749*x33;
    double x35 = exp(x23);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x23);
    double x39 = x38*x8;
    double x40 = -85364.155348523156*x37 + 74.830163563379159*x39;
    double x41 = 28454.71844950772*x16*x8;

    result += x0*(438605.28045614983*x0 - x13*x18 + x14*x20 + x14*x9 - 13628743.273679137*x2 - x20*x27 + x20*x34 - x20*x40 + x26*x28 - x27*x9 + 23224572.860757835*x3 + x34*x9 - x40*x9 - x41*(x15*(-9.0000000000000018*x31 + 2.3668182599152159*x33) - 0.78893941997173855*x19*x32 + 11.407719999999999*x29*x8/((x30)*(x30))) + x41*(x15*(-9.0*x37 + 0.0078893941997173842*x39) - 0.0026297980665724616*x19*x38 + 3422.3159999999998*x21*x35*x8/((x36)*(x36))) + x28*exp(-2281.5439999999999*x22)/((x25)*(x25)) - x18*exp(-7.6051466666666663*x7)/((x12)*(x12)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-3159034071.1583066*x2 + 1031606994.7059262*x3 + 1592618996.6270962);
    double x5 = 1.0/T;
    double x6 = 10.813310026713001*x2;
    double x7 = 3.5311699741784297*x3;
    double x8 = x6 - x7 - 5.4515027622500014;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 1140.7719999999999*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 2281.5439999999999*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.0078893941997173842*x19;
    double x25 = x5*(T*x24 - 9.0/x21);
    double x26 = 28454.71844950772*x9;
    double x27 = x17*(-x6 + x7 + 5.4515027622500014);

    result += x0*(-333268685214.86914*x15*x18 - x15*x4 - 111089561738.28972*x16*x18 - x16*x4 + x26*(0.0026297980665724616*x19 - 3422.3159999999998*x23 + x25) - x26*(-3904082.2679519993*x22*x27 - 3422.3160000000007*x23 + x24 + 3.0*x25 + 7808164.5359039987*x27*exp(x14)/((x21)*(x21)*(x21))) - 222179123476.57944*x18*exp(-3422.3159999999998*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(458491997.64707834*x2 - 702007571.36851263);
    double x5 = 10.813310026713001*x2;
    double x6 = 3.5311699741784297*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 5.4515027622500014;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 1140.7719999999999*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 2281.5439999999999*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 2.3541133161189531*x2 - 3.6044366755710002;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0026297980665724616*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 3422.3159999999998*x3;
    double x26 = 0.0078893941997173842*T*x20;
    double x27 = x19*x26 - 9.0/x23;
    double x28 = x0*x27;
    double x29 = 28454.71844950772*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 5.4515027622500014);

    result += x0*x1*x2*(333268685214.86914*x14*x18 - x14*x4 + 111089561738.28972*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(10266.948*x0*x31 - x26*x30 + 3.0*x27*x32) + 3904082.2679519993*x16*x24 - 7808164.5359039987*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 222179123476.57944*x18*exp(-3422.3159999999998*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(534907330.58825797*x2 - 585006309.47376049);
    double x4 = 10.813310026713001*x2;
    double x5 = 3.5311699741784297*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 5.4515027622500014;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 1140.7719999999999*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 2281.5439999999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 2.3541133161189531*x2 - 3.6044366755710002;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0026297980665724616*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 3422.3159999999998*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 5.4515027622500014;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0*x32 + 0.0078893941997173842*x33));
    double x35 = 5.4929310709442234*x2 - 6.0073944592850008;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0026297980665724616*x33;
    double x40 = 3.0*x31;
    double x41 = 3.0*x39/((x30)*(x30));

    result += -x38*(333268685214.86914*x13*x20 + x13*x3 + 111089561738.28972*x14*x20 + x14*x3 + 28454.71844950772*x15*x34 + 28454.71844950772*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(7.0623399483568594*x2 - 10.813310026713001)/pow(x6, 5.0/2.0) + x24*x35 - 3904082.2679519993*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(4.7082266322379063*x2 - 7.2088733511420005) + 7808164.5359039987*x37*exp(x12)/((x26)*(x26)*(x26))) + 222179123476.57944*x20*exp(-3422.3159999999998*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 10.813310026713001*x2;
    double x5 = 3.5311699741784297*x3;
    double x6 = x4 - x5 - 5.4515027622500014;
    double x7 = sqrt(x6);
    double x8 = 3.8025733333333331*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 2.3541133161189531*x2 - 3.6044366755710002;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 7.6051466666666663*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 5.4515027622500014;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 324603.46075081819*x16;
    double x21 = x9/x10;
    double x22 = 324603.46075081819*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(18.309770236480745*x2 - 16.019718558093334);
    double x25 = x23*x24;
    double x26 = 85364.155348523156*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 1140.7719999999999*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 85364.155348523171*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 2281.5439999999999*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 97381038.225245476*x27;
    double x39 = x19*x38;
    double x40 = x18*(4.7082266322379063*x2 - 7.2088733511420005);
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x40*x42;
    double x44 = (7.0623399483568594*x2 - 10.813310026713001)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 1.0/x17;
    double x47 = 5.4929310709442234*x2 - 6.0073944592850008;
    double x48 = x46*x47;
    double x49 = x12*x3;
    double x50 = x20*x49;
    double x51 = x46*(10.985862141888447*x2 - 12.014788918570002);
    double x52 = x47*x49;
    double x53 = x11*x52;
    double x54 = x22*x49;
    double x55 = x38*x43;
    double x56 = x38*x49;
    double x57 = x37*x56;
    double x58 = x32*x56;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 0.78893941997173855*x63;
    double x65 = 28454.71844950772*x23;
    double x66 = x24*x65;
    double x67 = exp(x29);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x29);
    double x71 = T*x23*x70;
    double x72 = -3*x69 + 0.0026297980665724616*x71;
    double x73 = 28454.71844950772*x45;
    double x74 = 85364.155348523156*x53;
    double x75 = 0.78893941997173855*x62;
    double x76 = x11*x75;
    double x77 = x59/((x60)*(x60));
    double x78 = 11.407719999999999*x77;
    double x79 = x23*x78;
    double x80 = x46*(-9.0000000000000018*x61 + 2.3668182599152159*x63) - x76 + x79;
    double x81 = 56909.43689901544*x80;
    double x82 = x23*x52;
    double x83 = 0.0026297980665724616*T*x70;
    double x84 = x11*x83;
    double x85 = x67/((x68)*(x68));
    double x86 = 3422.3159999999998*x27*x85;
    double x87 = x23*x86;
    double x88 = x46*(-9.0*x69 + 0.0078893941997173842*x71) - x84 + x87;
    double x89 = 56909.43689901544*x88;
    double x90 = x12*x2;
    double x91 = x44*x90;
    double x92 = x2*x41;
    double x93 = x46*x92;
    double x94 = x11*x92;
    double x95 = 3.0000000000000004*x64;
    double x96 = x18*x92;
    double x97 = x40*x90;
    double x98 = x49*x65;
    double x99 = x34*x93;
    double x100 = 3.0*x72;

    result += (-1754421.1218245993*x0 - 3702985.3912763237*x14*x16 - 1234328.4637587746*x14*x21 + x14*x81 - x14*x89 - x19*x20 - x19*x22 + 36343315.396477699*x2 - x20*x43 - 256092.46604556945*x21*x53 - x22*x43 - x25*x26 + x25*x33 - x26*x45 - 77415242.869192779*x3 + 111089561738.28972*x32*x35 + x32*x39 + 256092.46604556951*x32*x53 + x32*x55 + x33*x45 + 333268685214.86914*x35*x37 + x37*x39 + x37*x55 + x48*x50 + x48*x54 - x48*x57 - x48*x58 + x50*x51 + x51*x54 - x51*x57 - x51*x58 - x64*x66 - x64*x73 - x64*x74 + x66*x72 + x72*x73 + x72*x74 + x81*x82 - x82*x89 + x98*(-x100*x48 + x100*x96 + x100*x97 + x47*x84 - x47*x87 + x83*x91 - 3904082.2679519993*x85*x99 - x86*x94 + 3.0*x88*x93 + 7808164.5359039987*x99*exp(x36)/((x68)*(x68)*(x68))) - x98*(x47*x76 - x47*x79 - x48*x95 + x75*x91 - 43.37869186613333*x77*x93 - x78*x94 + 3.0000000000000004*x80*x93 + x95*x96 + x95*x97 + 86.757383732266661*x93*exp(x15)/((x60)*(x60)*(x60))) + 222179123476.57944*x35*exp(-3422.3159999999998*x28)/((x31)*(x31)*(x31)) - 2468656.9275175491*x14*exp(-11.407719999999999*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*1.367;
static const double Vmax = 1.15*1.367;
static double V = 0.9*1.367;

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



const char *Seifertite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Seifertite_slb_em_coder_calib_name(void) {
    return "Seifertite_slb_em";
}

const char *Seifertite_slb_em_coder_calib_formula(void) {
    return "SiO2";
}

const double Seifertite_slb_em_coder_calib_mw(void) {
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

const double *Seifertite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Seifertite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Seifertite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Seifertite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Seifertite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Seifertite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Seifertite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Seifertite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Seifertite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Seifertite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Seifertite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Seifertite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Seifertite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Seifertite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Seifertite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Seifertite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Seifertite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Seifertite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Seifertite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Seifertite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Seifertite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Seifertite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Seifertite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Seifertite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Seifertite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Seifertite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

