
static char *identifier = "Quartz_slb_em.emml:0628430a9fedbc9e6f8086f784b55d0f2e3428cb:Thu Feb 10 16:53:45 2022";



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
    double x3 = sqrt(-0.05534104480267947*x1 + 0.035139421227298628*x2 + 1.0200194272000001);
    double x4 = 2.7211023333333335*x3;
    double x5 = 816.33069999999998*x3/T;

    result += 74.830163563379159*T*log(1 - exp(-x5)) - 24.943387854459719*T*Debye(x5) - 3521303.0018208805*x1 + 2092143.6214840862*x2 - 22449.049069013745*log(1 - exp(-x4)) + 7483.0163563379156*Debye(x4) + 241812.88648584939 + 1225431.7176051149/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(0.035139421227298628*pow(x0, 4.0/3.0) - 0.05534104480267947*pow(x0, 2.0/3.0) + 1.0200194272000001);
    double x2 = x1/T;
    double x3 = 816.33069999999998*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -61086.159802807801*x2*x5/x6 + 20362.053267602601*x2*(-0.0036749812300333673*T*x4/x1 + 3/(exp(x3) - 1)) - 24.943387854459719*x4 + 74.830163563379159*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(-0.05534104480267947*x1 + 0.035139421227298628*x3 + 1.0200194272000001);
    double x6 = 2.7211023333333335*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(0.018447014934226489*x2 - 0.023426280818199086*x4);
    double x10 = 61086.159802807801*x9;
    double x11 = 816.33069999999998*x5/T;
    double x12 = exp(-x11);
    double x13 = 20362.053267602601*x9;

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + x13*(-1.1024943690100102*x8*Debye(x6) + 3/(exp(x6) - 1)) - x13*(-0.0036749812300333673*T*x8*Debye(x11) + 3/(exp(x11) - 1)) + 2347535.3345472533*x2 - 2789524.8286454482*x4 - 2450863.4352102298/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(-0.05534104480267947*x2 + 0.035139421227298628*x3 + 1.0200194272000001);
    double x5 = x0*x4;
    double x6 = 816.33069999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-2759664.6308096624*x2 + 1752280.2154144207*x3 + 50864806.510597013);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1632.6614*x5)/((x8)*(x8)) + 20362.053267602601*x4*(x0*(0.011024943690100101*T*x11 - 8.9999999999999982/x13) + 0.0036749812300333673*x11 - 2448.9920999999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 0.035139421227298628*pow(x1, 4.0/3.0) - 0.05534104480267947*x2 + 1.0200194272000001;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 816.33069999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 49866507.592137955*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(0.023426280818199086*x2 - 0.018447014934226489)*(20362.053267602601*x4*(-2448.9920999999999*x0*x11*x12/((x13)*(x13)) + 0.0036749812300333673*x10/pow(x3, 3.0/2.0) + (0.011024943690100101*x10*x11 - 8.9999999999999982/x13)/x3) + x7*x9/x8 + x9*exp(-1632.6614*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -0.05534104480267947*x2 + 0.035139421227298628*x3 + 1.0200194272000001;
    double x5 = sqrt(x4);
    double x6 = 2.7211023333333335*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = x7/x8;
    double x10 = 1.0/x5;
    double x11 = x2*(0.0546613219091312*x2 - 0.030745024890377483);
    double x12 = 61086.159802807801*x10*x11;
    double x13 = 1.0/x4;
    double x14 = x3*((0.023426280818199086*x2 - 0.018447014934226489)*(0.023426280818199086*x2 - 0.018447014934226489));
    double x15 = x13*x14;
    double x16 = 166221.6919737932*x15;
    double x17 = pow(x4, -3.0/2.0);
    double x18 = x14*x17;
    double x19 = 61086.159802807801*x18;
    double x20 = 1.0/T;
    double x21 = x20*x5;
    double x22 = 816.33069999999998*x21;
    double x23 = exp(-x22);
    double x24 = -x23 + 1;
    double x25 = x23/x24;
    double x26 = 49866507.592137955*x15*x20;
    double x27 = exp(x6);
    double x28 = x27 - 1;
    double x29 = 1.0/x28;
    double x30 = Debye(x6);
    double x31 = x10*x30;
    double x32 = -3*x29 + 1.1024943690100102*x31;
    double x33 = 20362.053267602601*x10;
    double x34 = x11*x33;
    double x35 = 20362.053267602601*x18;
    double x36 = exp(x22);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = Debye(x22);
    double x40 = T*x10*x39;
    double x41 = -3*x38 + 0.0036749812300333673*x40;
    double x42 = x14*x33;

    result += x0*(7352590.3056306895*x0 + x12*x25 - x12*x9 + x16*x9 + x16*exp(-5.442204666666667*x5)/((x8)*(x8)) - x19*x25 + x19*x9 - 3912558.8909120886*x2 - x25*x26 + 6508891.2668393794*x3 - x32*x34 + x32*x35 + x34*x41 - x35*x41 - x42*(0.0036749812300333673*T*x17*x39 - 2448.9920999999999*x10*x20*x36/((x37)*(x37)) + x13*(-8.9999999999999982*x38 + 0.011024943690100101*x40)) + x42*(-8.1633069999999996*x10*x27/((x28)*(x28)) + x13*(-9.0*x29 + 3.3074831070300306*x31) + 1.1024943690100102*x17*x30) - x26*exp(-1632.6614*x21)/((x24)*(x24)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-8278993.8924289877*x2 + 5256840.6462432621*x3 + 152594419.53179103);
    double x5 = 1.0/T;
    double x6 = -0.05534104480267947*x2 + 0.035139421227298628*x3 + 1.0200194272000001;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 816.33069999999998*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1632.6614*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.011024943690100101*x17;
    double x23 = x5*(T*x22 - 8.9999999999999982/x19);
    double x24 = 20362.053267602601*x7;
    double x25 = x15*x6;

    result += x0*(-122122683147.73587*x13*x16 + x13*x4 - 40707561049.245293*x14*x16 + x14*x4 + x24*(0.0036749812300333673*x17 - 2448.9920999999999*x21 + x23) - x24*(1999187.43528747*x20*x25 - 2448.9920999999986*x21 + x22 + 2.9999999999999996*x23 - 3998374.87057494*x25*exp(x12)/((x19)*(x19)*(x19))) - 81415122098.490585*x16*exp(-2448.9920999999999*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(2336373.6205525608*x2 - 1839776.4205397749);
    double x5 = 0.035139421227298628*pow(x1, 4.0/3.0) - 0.05534104480267947*x2 + 1.0200194272000001;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 816.33069999999998*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1632.6614*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 0.023426280818199086*x2 - 0.018447014934226489;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0036749812300333673*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2448.9920999999999*x22*x3;
    double x24 = 0.011024943690100101*T*x18;
    double x25 = x17*x24 - 8.9999999999999982/x21;
    double x26 = x0*x25;
    double x27 = 20362.053267602601*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-122122683147.73587*x12*x16 + x12*x4 - 40707561049.245293*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-7346.9762999999984*x0*x17*x22 + x24*x28 + 2.9999999999999996*x25*x29) + 1999187.43528747*x15*x22 - 3998374.87057494*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 81415122098.490585*x16*exp(-2448.9920999999999*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(2725769.2239779877*x2 - 1533147.0171164791);
    double x4 = 0.035139421227298628*pow(x1, 4.0/3.0) - 0.05534104480267947*x2 + 1.0200194272000001;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 816.33069999999998*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1632.6614*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 0.023426280818199086*x2 - 0.018447014934226489;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0036749812300333673*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2448.9920999999999*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-8.9999999999999982*x29 + 0.011024943690100101*x30));
    double x32 = 0.0546613219091312*x2 - 0.030745024890377483;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0036749812300333673*x30;
    double x37 = 2.9999999999999996*x28;
    double x38 = 2.9999999999999996*x36/((x4)*(x4));

    result += x35*(-122122683147.73587*x11*x18 + x11*x3 - 40707561049.245293*x12*x18 + x12*x3 + 20362.053267602601*x13*x31 - 81415122098.490585*x18*exp(-2448.9920999999999*x6)/((x9)*(x9)*(x9)) - 20362.053267602601*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(0.070278842454597257*x2 - 0.05534104480267947)/pow(x4, 5.0/2.0) - x22*x32 + 1999187.43528747*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(0.046852561636398171*x2 - 0.036894029868452978) - 3998374.87057494*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = -0.05534104480267947*x2 + 0.035139421227298628*x3 + 1.0200194272000001;
    double x5 = sqrt(x4);
    double x6 = 2.7211023333333335*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(x4, -3.0/2.0);
    double x10 = 0.023426280818199086*x2 - 0.018447014934226489;
    double x11 = x0*((x10)*(x10)*(x10));
    double x12 = x11*x9;
    double x13 = pow(x4, -2);
    double x14 = x11*x13;
    double x15 = 166221.6919737932*x14;
    double x16 = 5.442204666666667*x5;
    double x17 = exp(-x16)/((x8)*(x8));
    double x18 = x7/x8;
    double x19 = 61086.159802807801*x18;
    double x20 = 1.0/x5;
    double x21 = x2*(0.18220440636377067*x2 - 0.081986733041006621);
    double x22 = x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x5;
    double x25 = 816.33069999999998*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 61086.159802807801*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1632.6614*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = 49866507.592137955*x23;
    double x35 = x14*x34;
    double x36 = x13*(0.046852561636398171*x2 - 0.036894029868452978);
    double x37 = ((x10)*(x10));
    double x38 = x0*x37;
    double x39 = x36*x38;
    double x40 = 166221.6919737932*x17;
    double x41 = 166221.6919737932*x18;
    double x42 = (0.070278842454597257*x2 - 0.05534104480267947)/pow(x4, 5.0/2.0);
    double x43 = x38*x42;
    double x44 = 1.0/x4;
    double x45 = 0.0546613219091312*x2 - 0.030745024890377483;
    double x46 = x44*x45;
    double x47 = x10*x3;
    double x48 = x40*x47;
    double x49 = x44*(0.1093226438182624*x2 - 0.061490049780754966);
    double x50 = x45*x47;
    double x51 = x50*x9;
    double x52 = 183258.47940842342*x51;
    double x53 = x41*x47;
    double x54 = x34*x39;
    double x55 = x34*x47;
    double x56 = x33*x55;
    double x57 = x28*x55;
    double x58 = exp(x6);
    double x59 = x58 - 1;
    double x60 = 1.0/x59;
    double x61 = Debye(x6);
    double x62 = x20*x61;
    double x63 = -3*x60 + 1.1024943690100102*x62;
    double x64 = 20362.053267602601*x20;
    double x65 = x21*x64;
    double x66 = exp(x25);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x25);
    double x70 = T*x20*x69;
    double x71 = -3*x68 + 0.0036749812300333673*x70;
    double x72 = 20362.053267602601*x43;
    double x73 = 61086.159802807801*x51;
    double x74 = 1.1024943690100102*x61;
    double x75 = x74*x9;
    double x76 = x58/((x59)*(x59));
    double x77 = 8.1633069999999996*x76;
    double x78 = x20*x77;
    double x79 = x44*(-9.0*x60 + 3.3074831070300306*x62) + x75 - x78;
    double x80 = 40724.106535205203*x79;
    double x81 = x20*x50;
    double x82 = 0.0036749812300333673*T*x69;
    double x83 = x82*x9;
    double x84 = x66/((x67)*(x67));
    double x85 = 2448.9920999999999*x23*x84;
    double x86 = x20*x85;
    double x87 = x44*(-8.9999999999999982*x68 + 0.011024943690100101*x70) + x83 - x86;
    double x88 = 40724.106535205203*x87;
    double x89 = x10*x2;
    double x90 = x42*x89;
    double x91 = x2*x37;
    double x92 = x44*x91;
    double x93 = x9*x91;
    double x94 = 3.0*x63;
    double x95 = x13*x91;
    double x96 = x36*x89;
    double x97 = x47*x64;
    double x98 = x30*x92;
    double x99 = 2.9999999999999996*x71;

    result += (-29410361.222522758*x0 + 1356918.70164151*x12*x17 + 452306.23388050334*x12*x18 + x12*x80 - x12*x88 + 904612.46776100667*x12*exp(-8.1633069999999996*x5)/((x8)*(x8)*(x8)) + x15*x17 + x15*x18 - x18*x52 + x19*x22 + x19*x43 + 10433490.37576557*x2 - x22*x29 - 40707561049.245293*x28*x31 - x28*x35 + x28*x52 - x28*x54 - x29*x43 - 21696304.22279793*x3 - 122122683147.73587*x31*x33 - x33*x35 - x33*x54 + x39*x40 + x39*x41 - x46*x48 - x46*x53 + x46*x56 + x46*x57 - x48*x49 - x49*x53 + x49*x56 + x49*x57 + x63*x65 + x63*x72 - x63*x73 - x65*x71 - x71*x72 + x71*x73 - x80*x81 + x81*x88 + x97*(-x45*x75 + x45*x78 - x46*x94 + x74*x90 + 22.213193725416335*x76*x92 - x77*x93 + 3.0*x79*x92 + x94*x95 + x94*x96 - 44.426387450832671*x92*exp(x16)/((x59)*(x59)*(x59))) - x97*(-x45*x83 + x45*x86 - x46*x99 + x82*x90 + 1999187.43528747*x84*x98 - x85*x93 + 2.9999999999999996*x87*x92 + x95*x99 + x96*x99 - 3998374.87057494*x98*exp(x32)/((x67)*(x67)*(x67))) - 81415122098.490585*x31*exp(-2448.9920999999999*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*2.367003;
static const double Vmax = 1.15*2.367003;
static double V = 0.9*2.367003;

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



const char *Quartz_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Quartz_slb_em_coder_calib_name(void) {
    return "Quartz_slb_em";
}

const char *Quartz_slb_em_coder_calib_formula(void) {
    return "SiO2";
}

const double Quartz_slb_em_coder_calib_mw(void) {
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

const double *Quartz_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Quartz_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Quartz_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Quartz_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Quartz_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Quartz_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Quartz_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Quartz_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Quartz_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Quartz_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Quartz_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Quartz_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Quartz_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Quartz_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Quartz_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Quartz_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Quartz_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Quartz_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Quartz_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Quartz_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Quartz_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Quartz_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Quartz_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Quartz_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Quartz_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Quartz_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Quartz_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Quartz_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Quartz_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Quartz_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Quartz_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Quartz_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Quartz_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Quartz_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Quartz_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Quartz_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Quartz_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

