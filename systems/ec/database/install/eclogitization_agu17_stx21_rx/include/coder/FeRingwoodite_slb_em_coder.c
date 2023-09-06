
static char *identifier = "FeRingwoodite_slb_em.emml:51d60f364716f7dd5abab309a60aec0e9f893433:Fri May 26 02:49:17 2023";



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
    double x3 = sqrt(17.637097582312503*x1 - 10.095377057099341*x2 - 4.2939517436000001);
    double x4 = 2.2094768*x3;
    double x5 = 662.84303999999997*x3/T;

    result += 174.60371498121802*T*log(1 - exp(-x5)) - 58.201238327072673*T*Debye(x5) - 26.76322271834383*T - 43579648.205268979*x1 + 45391737.66504477*x2 - 52381.114494365407*log(1 - exp(-x4)) + 17460.371498121804*Debye(x4) + 7583582.8396871462 + 19401212.81326795/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-10.095377057099341*pow(x0, 4.0/3.0) + 17.637097582312503*pow(x0, 2.0/3.0) - 4.2939517436000001);
    double x2 = x1/T;
    double x3 = 662.84303999999997*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -115734.85723344408*x2*x5/x6 + 38578.285744481364*x2*(-0.0045259583626313712*T*x4/x1 + 3/(exp(x3) - 1)) - 58.201238327072673*x4 + 174.60371498121802*log(x6) - 26.76322271834383;
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(17.637097582312503*x1 - 10.095377057099341*x3 - 4.2939517436000001);
    double x6 = 2.2094768*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-5.8790325274375004*x2 + 6.7302513713995609*x4);
    double x10 = 662.84303999999997*x5/T;
    double x11 = exp(-x10);

    result += 115734.85723344408*x11*x9/(-x11 + 1) + 29053098.803512651*x2 - 60522316.886726357*x4 - 115734.8572334441*x7*x9/(-x7 + 1) + 38578.285744481371*x9*(-1.3577875087894111*x8*Debye(x6) + 3/(exp(x6) - 1)) - 38578.285744481364*x9*(-0.0045259583626313712*T*x8*Debye(x10) + 3/(exp(x10) - 1)) - 38802425.6265359/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(17.637097582312503*x2 - 10.095377057099341*x3 - 4.2939517436000001);
    double x5 = x0*x4;
    double x6 = 662.84303999999997*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-1353013090.5896137*x2 + 774457205.8382026*x3 + 329406405.57986546);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1325.6860799999999*x5)/((x8)*(x8)) - 38578.285744481364*x4*(x0*(0.013577875087894115*T*x11 - 9.0000000000000018/x13) + 0.0045259583626313712*x11 - 1988.5291199999999*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 17.637097582312503*x2;
    double x4 = 10.095377057099341*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 4.2939517436000001;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 662.84303999999997*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 76714044.602582067*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(6.7302513713995609*x2 - 5.8790325274375004)*(-38578.285744481364*x6*(1988.5291199999999*x0*x13*x14/((x15)*(x15)) - 0.0045259583626313712*x12/pow(x5, 3.0/2.0) + (0.013577875087894115*x12*x13 - 9.0000000000000018/x15)/(-x3 + x4 + 4.2939517436000001)) + x11*x9/x10 + x11*exp(-1325.6860799999999*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 17.637097582312503*x2;
    double x5 = 10.095377057099341*x3;
    double x6 = x4 - x5 - 4.2939517436000001;
    double x7 = sqrt(x6);
    double x8 = 1.0/x7;
    double x9 = x2*x8*(15.703919866598975*x2 - 9.7983875457291667);
    double x10 = 2.2094768*x7;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = x11/x12;
    double x14 = 115734.8572334441*x13;
    double x15 = 1.0/(-x4 + x5 + 4.2939517436000001);
    double x16 = x3*((6.7302513713995609*x2 - 5.8790325274375004)*(6.7302513713995609*x2 - 5.8790325274375004));
    double x17 = x15*x16;
    double x18 = 255713.48200860692*x17;
    double x19 = pow(x6, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 1.0/T;
    double x22 = x21*x7;
    double x23 = 662.84303999999997*x22;
    double x24 = exp(-x23);
    double x25 = -x24 + 1;
    double x26 = x24/x25;
    double x27 = 115734.85723344408*x26;
    double x28 = 76714044.602582067*x17*x21;
    double x29 = exp(x10);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x10);
    double x33 = x32*x8;
    double x34 = -115734.85723344411*x31 + 52381.114494365414*x33;
    double x35 = exp(x23);
    double x36 = x35 - 1;
    double x37 = 1.0/x36;
    double x38 = T*Debye(x23);
    double x39 = x38*x8;
    double x40 = -115734.85723344408*x37 + 174.60371498121805*x39;
    double x41 = x16*x8;

    result += x0*(116407276.87960771*x0 - x13*x18 + x14*x20 + x14*x9 - 48421831.339187756*x2 - x20*x27 + x20*x34 - x20*x40 + x26*x28 - x27*x9 + 141218739.40236148*x3 + x34*x9 - x40*x9 - 38578.285744481371*x41*(x15*(-9.0*x31 + 4.0733625263682338*x33) - 1.3577875087894111*x19*x32 + 6.6284304000000001*x29*x8/((x30)*(x30))) + 38578.285744481364*x41*(x15*(-9.0000000000000018*x37 + 0.013577875087894115*x39) - 0.0045259583626313712*x19*x38 + 1988.5291199999999*x21*x35*x8/((x36)*(x36))) + x28*exp(-1325.6860799999999*x22)/((x25)*(x25)) - x18*exp(-4.4189536*x7)/((x12)*(x12)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-4059039271.7688413*x2 + 2323371617.5146079*x3 + 988219216.73959637);
    double x5 = 1.0/T;
    double x6 = 17.637097582312503*x2;
    double x7 = 10.095377057099341*x3;
    double x8 = x6 - x7 - 4.2939517436000001;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 662.84303999999997*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1325.6860799999999*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.013577875087894115*x19;
    double x25 = x5*(T*x24 - 9.0000000000000018/x21);
    double x26 = 38578.285744481364*x9;
    double x27 = x17*(-x6 + x7 + 4.2939517436000001);

    result += x0*(-152548111605.21326*x15*x18 - x15*x4 - 50849370535.071091*x16*x18 - x16*x4 + x26*(0.0045259583626313712*x19 - 1988.5291199999999*x23 + x25) - x26*(-1318082.6870293247*x22*x27 - 1988.5291200000006*x23 + x24 + 3.0000000000000004*x25 + 2636165.3740586494*x27*exp(x14)/((x21)*(x21)*(x21))) - 101698741070.14218*x18*exp(-1988.5291199999999*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(1032609607.78427*x2 - 902008727.05974233);
    double x5 = 17.637097582312503*x2;
    double x6 = 10.095377057099341*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 4.2939517436000001;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 662.84303999999997*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1325.6860799999999*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 6.7302513713995609*x2 - 5.8790325274375004;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0045259583626313712*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 1988.5291199999999*x3;
    double x26 = 0.013577875087894115*T*x20;
    double x27 = x19*x26 - 9.0000000000000018/x23;
    double x28 = x0*x27;
    double x29 = 38578.285744481364*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 4.2939517436000001);

    result += x0*x1*x2*(152548111605.21326*x14*x18 - x14*x4 + 50849370535.071091*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(5965.5873600000004*x0*x31 - x26*x30 + 3.0000000000000004*x27*x32) + 1318082.6870293247*x16*x24 - 2636165.3740586494*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 101698741070.14218*x18*exp(-1988.5291199999999*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(1204711209.0816483*x2 - 751673939.21645188);
    double x4 = 17.637097582312503*x2;
    double x5 = 10.095377057099341*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 4.2939517436000001;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 662.84303999999997*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1325.6860799999999*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 6.7302513713995609*x2 - 5.8790325274375004;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0045259583626313712*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 1988.5291199999999*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 4.2939517436000001;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-9.0000000000000018*x32 + 0.013577875087894115*x33));
    double x35 = 15.703919866598975*x2 - 9.7983875457291667;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0045259583626313712*x33;
    double x40 = 3.0000000000000004*x31;
    double x41 = 3.0000000000000004*x39/((x30)*(x30));

    result += -x38*(152548111605.21326*x13*x20 + x13*x3 + 50849370535.071091*x14*x20 + x14*x3 + 38578.285744481364*x15*x34 + 38578.285744481364*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(20.190754114198683*x2 - 17.637097582312499)/pow(x6, 5.0/2.0) + x24*x35 - 1318082.6870293247*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(13.460502742799122*x2 - 11.758065054875001) + 2636165.3740586494*x37*exp(x12)/((x26)*(x26)*(x26))) + 101698741070.14218*x20*exp(-1988.5291199999999*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 17.637097582312503*x2;
    double x5 = 10.095377057099341*x3;
    double x6 = x4 - x5 - 4.2939517436000001;
    double x7 = sqrt(x6);
    double x8 = 2.2094768*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 6.7302513713995609*x2 - 5.8790325274375004;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 4.4189536*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 4.2939517436000001;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 255713.48200860692*x16;
    double x21 = x9/x10;
    double x22 = 255713.48200860692*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(52.34639955532991*x2 - 26.129033455277778);
    double x25 = x23*x24;
    double x26 = 115734.8572334441*x21;
    double x27 = 1.0/T;
    double x28 = x27*x7;
    double x29 = 662.84303999999997*x28;
    double x30 = exp(-x29);
    double x31 = -x30 + 1;
    double x32 = x30/x31;
    double x33 = 115734.85723344408*x32;
    double x34 = pow(T, -2);
    double x35 = x14*x34;
    double x36 = 1325.6860799999999*x28;
    double x37 = exp(-x36)/((x31)*(x31));
    double x38 = 76714044.602582067*x27;
    double x39 = x19*x38;
    double x40 = x18*(13.460502742799122*x2 - 11.758065054875001);
    double x41 = ((x12)*(x12));
    double x42 = x0*x41;
    double x43 = x40*x42;
    double x44 = (20.190754114198683*x2 - 17.637097582312499)/pow(x6, 5.0/2.0);
    double x45 = x42*x44;
    double x46 = 1.0/x17;
    double x47 = 15.703919866598975*x2 - 9.7983875457291667;
    double x48 = x46*x47;
    double x49 = x12*x3;
    double x50 = x20*x49;
    double x51 = x46*(31.40783973319795*x2 - 19.596775091458333);
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
    double x64 = -3*x61 + 1.3577875087894111*x63;
    double x65 = 38578.285744481371*x23;
    double x66 = exp(x29);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x29);
    double x70 = T*x23*x69;
    double x71 = -3*x68 + 0.0045259583626313712*x70;
    double x72 = 38578.285744481364*x23;
    double x73 = 1.3577875087894111*x62;
    double x74 = x11*x73;
    double x75 = x59/((x60)*(x60));
    double x76 = 6.6284304000000001*x75;
    double x77 = x23*x76;
    double x78 = x46*(-9.0*x61 + 4.0733625263682338*x63) - x74 + x77;
    double x79 = 77156.571488962742*x78;
    double x80 = x23*x52;
    double x81 = 0.0045259583626313712*T*x69;
    double x82 = x11*x81;
    double x83 = x66/((x67)*(x67));
    double x84 = 1988.5291199999999*x27*x83;
    double x85 = x23*x84;
    double x86 = x46*(-9.0000000000000018*x68 + 0.013577875087894115*x70) - x82 + x85;
    double x87 = 77156.571488962727*x86;
    double x88 = x12*x2;
    double x89 = x44*x88;
    double x90 = x2*x41;
    double x91 = x46*x90;
    double x92 = x11*x90;
    double x93 = 3.0*x64;
    double x94 = x18*x90;
    double x95 = x40*x88;
    double x96 = x34*x91;
    double x97 = 3.0000000000000004*x71;

    result += (-465629107.51843083*x0 - 1694979.0178357032*x14*x16 - 564993.00594523444*x14*x21 + x14*x79 - x14*x87 - x19*x20 - x19*x22 + 129124883.57116735*x2 - x20*x43 - 347204.57170033228*x21*x53 - x22*x43 - x24*x64*x65 + x24*x71*x72 - x25*x26 + x25*x33 - x26*x45 - 470729131.34120494*x3 + 50849370535.071091*x32*x35 + x32*x39 + 347204.57170033222*x32*x53 + x32*x55 + x33*x45 + 152548111605.21326*x35*x37 + x37*x39 + x37*x55 - 38578.285744481371*x45*x64 + 38578.285744481364*x45*x71 + x48*x50 + x48*x54 - x48*x57 - x48*x58 - x49*x65*(x47*x74 - x47*x77 - x48*x93 + x73*x89 - 14.64536318921472*x75*x91 - x76*x92 + 3.0*x78*x91 + x93*x94 + x93*x95 + 29.29072637842944*x91*exp(x15)/((x60)*(x60)*(x60))) + x49*x72*(x47*x82 - x47*x85 - x48*x97 + x81*x89 - 1318082.6870293247*x83*x96 - x84*x92 + 3.0000000000000004*x86*x91 + x94*x97 + x95*x97 + 2636165.3740586494*x96*exp(x36)/((x67)*(x67)*(x67))) + x50*x51 + x51*x54 - x51*x57 - x51*x58 - 115734.85723344411*x53*x64 + 115734.85723344408*x53*x71 + x79*x80 - x80*x87 + 101698741070.14218*x35*exp(-1988.5291199999999*x28)/((x31)*(x31)*(x31)) - 1129986.0118904689*x14*exp(-6.6284304000000001*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*4.186;
static const double Vmax = 1.15*4.186;
static double V = 0.9*4.186;

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



const char *FeRingwoodite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *FeRingwoodite_slb_em_coder_calib_name(void) {
    return "FeRingwoodite_slb_em";
}

const char *FeRingwoodite_slb_em_coder_calib_formula(void) {
    return "Fe2SiO4";
}

const double FeRingwoodite_slb_em_coder_calib_mw(void) {
    return 203.77710000000002;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,4.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,2.0,0.0,0.0,0.0,
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

const double *FeRingwoodite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double FeRingwoodite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double FeRingwoodite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double FeRingwoodite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double FeRingwoodite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double FeRingwoodite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double FeRingwoodite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double FeRingwoodite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double FeRingwoodite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double FeRingwoodite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double FeRingwoodite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double FeRingwoodite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double FeRingwoodite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double FeRingwoodite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int FeRingwoodite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **FeRingwoodite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **FeRingwoodite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void FeRingwoodite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int FeRingwoodite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double FeRingwoodite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int FeRingwoodite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double FeRingwoodite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double FeRingwoodite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

