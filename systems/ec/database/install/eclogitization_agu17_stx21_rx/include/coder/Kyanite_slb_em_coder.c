
static char *identifier = "Kyanite_slb_em.emml:f4fb001af5a8e7d39c63dbff33e00b56288c7fc7:Tue Aug 22 19:33:46 2023";



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
    double x3 = sqrt(5.4126238802913091*x1 + 2.786298517811455*x2 - 1.3926646695499998);
    double x4 = 3.1439864333333336*x3;
    double x5 = 943.19592999999998*x3/T;

    result += 199.54710283567775*T*log(1 - exp(-x5)) - 66.515700945225916*T*Debye(x5) - 42899112.053271465*x1 + 57793184.874471515*x2 - 59864.130850703325*log(1 - exp(-x4)) + 19954.710283567776*Debye(x4) + 5514778.2599999998;
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(2.786298517811455*pow(x0, 4.0/3.0) + 5.4126238802913091*pow(x0, 2.0/3.0) - 1.3926646695499998);
    double x2 = x1/T;
    double x3 = 943.19592999999998*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -188212.0152379027*x2*x5/x6 + 62737.338412634235*x2*(-0.0031806753025323169*T*x4/x1 + 3/(exp(x3) - 1)) - 66.515700945225916*x4 + 199.54710283567775*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(5.4126238802913091*x1 + 2.786298517811455*x3 - 1.3926646695499998);
    double x6 = 3.1439864333333336*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-1.804207960097103*x2 - 1.8575323452076367*x4);
    double x10 = 943.19592999999998*x5/T;
    double x11 = exp(-x10);

    result += 188212.0152379027*x11*x9/(-x11 + 1) + 28599408.03551431*x2 - 77057579.832628682*x4 - 188212.01523790273*x7*x9/(-x7 + 1) + 62737.33841263425*x9*(-0.95420259075969494*x8*Debye(x6) + 3/(exp(x6) - 1)) - 62737.338412634235*x9*(-0.0031806753025323169*T*x8*Debye(x10) + 3/(exp(x10) - 1));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(5.4126238802913091*x2 + 2.786298517811455*x3 - 1.3926646695499998);
    double x5 = x0*x4;
    double x6 = 943.19592999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(960853357.86085641*x2 + 494625960.72679162*x3 - 247226955.67002481);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*(x10*x7/x8 + x10*exp(-1886.39186*x5)/((x8)*(x8)) + 62737.338412634235*x4*(x0*(0.0095420259075969516*T*x11 - 9.0/x13) + 0.0031806753025323169*x11 - 2829.58779*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 2.786298517811455*pow(x1, 4.0/3.0) + 5.4126238802913091*x2 - 1.3926646695499998;
    double x4 = sqrt(x3);
    double x5 = x0*x4;
    double x6 = 943.19592999999998*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = 177520806.74948782*x0;
    double x10 = T*Debye(x6);
    double x11 = 1.0/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += -x0*x1*x2*(1.8575323452076367*x2 + 1.804207960097103)*(62737.338412634235*x4*(-2829.58779*x0*x11*x12/((x13)*(x13)) + 0.0031806753025323169*x10/pow(x3, 3.0/2.0) + (0.0095420259075969516*x10*x11 - 9.0/x13)/x3) + x7*x9/x8 + x9*exp(-1886.39186*x5)/((x8)*(x8)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = 2.786298517811455*pow(x0, 4.0/3.0) + 5.4126238802913091*x1 - 1.3926646695499998;
    double x3 = sqrt(x2);
    double x4 = 1.0/x3;
    double x5 = 3.1439864333333336*x3;
    double x6 = exp(-x5);
    double x7 = -x6 + 1;
    double x8 = x6/x7;
    double x9 = 1.0/T;
    double x10 = x3*x9;
    double x11 = 943.19592999999998*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = x12/x13;
    double x15 = 1.0/x2;
    double x16 = x1*((1.8575323452076367*x1 + 1.804207960097103)*(1.8575323452076367*x1 + 1.804207960097103));
    double x17 = x15*x16;
    double x18 = 591736.02249829285*x17;
    double x19 = pow(x2, -3.0/2.0);
    double x20 = x16*x19;
    double x21 = 177520806.74948782*x17*x9;
    double x22 = x4*(4.3342421388178192*x1 + 3.0070132668285048);
    double x23 = exp(x5);
    double x24 = x23 - 1;
    double x25 = 1.0/x24;
    double x26 = Debye(x5);
    double x27 = x26*x4;
    double x28 = -188212.01523790276*x25 + 59864.130850703325*x27;
    double x29 = exp(x11);
    double x30 = x29 - 1;
    double x31 = 1.0/x30;
    double x32 = Debye(x11);
    double x33 = T*x32*x4;
    double x34 = -188212.0152379027*x31 + 199.54710283567775*x33;
    double x35 = x16*x4;

    result += x1*(179801019.60946691*x1 - 188212.0152379027*x14*x20 - x14*x21 + x14*x4*(815756.44747593941*x1 + 565956.02679690206) + x18*x8 + x18*exp(-6.2879728666666672*x3)/((x7)*(x7)) + x20*x28 - x20*x34 + 188212.01523790273*x20*x8 - x22*x28 + x22*x34 + 62737.33841263425*x35*(x15*(-9.0*x25 + 2.8626077722790848*x27) + 0.95420259075969494*x19*x26 - 9.4319593000000008*x23*x4/((x24)*(x24))) - 62737.338412634235*x35*(0.0031806753025323169*T*x19*x32 + x15*(-9.0*x31 + 0.0095420259075969516*x33) - 2829.58779*x29*x4*x9/((x30)*(x30))) - x4*x8*(815756.44747593952*x1 + 565956.02679690218) - 47665680.059190512 - x21*exp(-1886.39186*x10)/((x13)*(x13)))/((V)*(V));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(2882560073.5825691*x2 + 1483877882.1803749*x3 - 741680867.0100745);
    double x5 = 1.0/T;
    double x6 = 5.4126238802913091*x2 + 2.786298517811455*x3 - 1.3926646695499998;
    double x7 = sqrt(x6);
    double x8 = x5*x7;
    double x9 = 943.19592999999998*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1886.39186*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = pow(T, -3);
    double x16 = x15*pow(x6, 3.0/2.0);
    double x17 = Debye(x9)/x7;
    double x18 = exp(x9);
    double x19 = x18 - 1;
    double x20 = x18/((x19)*(x19));
    double x21 = x0*x20*x7;
    double x22 = 0.0095420259075969516*x17;
    double x23 = x5*(T*x22 - 9.0/x19);
    double x24 = 62737.338412634235*x7;
    double x25 = x15*x6;

    result += x0*(-502310707249.30029*x13*x16 + x13*x4 - 167436902416.43344*x14*x16 + x14*x4 + x24*(0.0031806753025323169*x17 - 2829.58779*x21 + x23) - x24*(2668855.6871056948*x20*x25 - 2829.5877900000005*x21 + x22 + 3.0*x23 - 5337711.3742113896*x25*exp(x12)/((x19)*(x19)*(x19))) - 334873804832.86688*x16*exp(-2829.58779*x8)/((x11)*(x11)*(x11)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(659501280.96905553*x2 + 640568905.2405709);
    double x5 = 2.786298517811455*pow(x1, 4.0/3.0) + 5.4126238802913091*x2 - 1.3926646695499998;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 943.19592999999998*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 1886.39186*x7;
    double x12 = exp(-x11)/((x10)*(x10));
    double x13 = x9/x10;
    double x14 = 1.8575323452076367*x2 + 1.804207960097103;
    double x15 = pow(T, -3);
    double x16 = x14*x15*x6;
    double x17 = 1.0/x6;
    double x18 = Debye(x8);
    double x19 = 0.0031806753025323169*x18;
    double x20 = exp(x8);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = 2829.58779*x22*x3;
    double x24 = 0.0095420259075969516*T*x18;
    double x25 = x17*x24 - 9.0/x21;
    double x26 = x0*x25;
    double x27 = 62737.338412634235*x14;
    double x28 = pow(x5, -3.0/2.0);
    double x29 = 1.0/x5;

    result += x0*x1*x2*(-502310707249.30029*x12*x16 + x12*x4 - 167436902416.43344*x13*x16 + x13*x4 + x17*x27*(x17*x19 - x23*x6 + x26) - x27*x6*(x0*(-8488.7633700000006*x0*x17*x22 + x24*x28 + 3.0*x25*x29) + 2668855.6871056948*x15*x22 - 5337711.3742113896*x15*exp(x11)/((x21)*(x21)*(x21)) + x17*x23 + x19*x28 + x26*x29) - 334873804832.86688*x16*exp(-2829.58779*x7)/((x10)*(x10)*(x10)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(769418161.13056481*x2 + 533807421.03380907);
    double x4 = 2.786298517811455*pow(x1, 4.0/3.0) + 5.4126238802913091*x2 - 1.3926646695499998;
    double x5 = sqrt(x4);
    double x6 = x0*x5;
    double x7 = 943.19592999999998*x6;
    double x8 = exp(-x7);
    double x9 = -x8 + 1;
    double x10 = 1886.39186*x6;
    double x11 = exp(-x10)/((x9)*(x9));
    double x12 = x8/x9;
    double x13 = 1.0/x5;
    double x14 = 1.8575323452076367*x2 + 1.804207960097103;
    double x15 = ((x14)*(x14));
    double x16 = x15*x2;
    double x17 = x16/((T)*(T));
    double x18 = x13*x17;
    double x19 = pow(x4, -3.0/2.0);
    double x20 = Debye(x7);
    double x21 = 0.0031806753025323169*T*x20;
    double x22 = x19*x21;
    double x23 = exp(x7);
    double x24 = x23 - 1;
    double x25 = x23/((x24)*(x24));
    double x26 = 2829.58779*x25;
    double x27 = x0*x13*x26;
    double x28 = 1.0/x4;
    double x29 = 1.0/x24;
    double x30 = T*x13*x20;
    double x31 = x16*(x22 - x27 + x28*(-9.0*x29 + 0.0095420259075969516*x30));
    double x32 = 4.3342421388178192*x2 + 3.0070132668285048;
    double x33 = x14*x2;
    double x34 = x17*x28;
    double x35 = x0*x2;
    double x36 = -3*x29 + 0.0031806753025323169*x30;
    double x37 = 3.0*x28;
    double x38 = 3.0*x36/((x4)*(x4));

    result += x35*(-502310707249.30029*x11*x18 + x11*x3 - 167436902416.43344*x12*x18 + x12*x3 + 62737.338412634235*x13*x31 - 334873804832.86688*x18*exp(-2829.58779*x6)/((x9)*(x9)*(x9)) - 62737.338412634235*x5*(-x15*x19*x26*x35 + x16*x38 + x21*x33*(5.57259703562291*x2 + 5.4126238802913091)/pow(x4, 5.0/2.0) - x22*x32 + 2668855.6871056948*x25*x34 + x27*x32 + x31*x37 - x32*x36*x37 + x33*x38*(3.7150646904152733*x2 + 3.6084159201942061) - 5337711.3742113896*x34*exp(x10)/((x24)*(x24)*(x24))))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = pow(x0, 4.0/3.0);
    double x3 = 5.4126238802913091*x1 + 2.786298517811455*x2 - 1.3926646695499998;
    double x4 = sqrt(x3);
    double x5 = 3.1439864333333336*x4;
    double x6 = exp(-x5);
    double x7 = -x6 + 1;
    double x8 = pow(x3, -3.0/2.0);
    double x9 = pow(V, -2);
    double x10 = 1.8575323452076367*x1 + 1.804207960097103;
    double x11 = ((x10)*(x10)*(x10))*x9;
    double x12 = x11*x8;
    double x13 = pow(x3, -2);
    double x14 = 591736.02249829285*x13;
    double x15 = x11*x14;
    double x16 = 6.2879728666666672*x4;
    double x17 = exp(-x16)/((x7)*(x7));
    double x18 = x6/x7;
    double x19 = 188212.01523790273*x18;
    double x20 = 1.0/x4;
    double x21 = 14.447473796059397*x1 + 8.0187020448760116;
    double x22 = x1*x20*x21;
    double x23 = 1.0/T;
    double x24 = x23*x4;
    double x25 = 943.19592999999998*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 188212.0152379027*x28;
    double x30 = pow(T, -2);
    double x31 = x12*x30;
    double x32 = 1886.39186*x24;
    double x33 = exp(-x32)/((x27)*(x27));
    double x34 = x11*x13;
    double x35 = 177520806.74948782*x23;
    double x36 = x33*x35;
    double x37 = x28*x35;
    double x38 = 3.7150646904152733*x1 + 3.6084159201942061;
    double x39 = ((x10)*(x10));
    double x40 = x39*x9;
    double x41 = x38*x40;
    double x42 = x14*x41;
    double x43 = (5.57259703562291*x1 + 5.4126238802913091)/pow(x3, 5.0/2.0);
    double x44 = x40*x43;
    double x45 = 4.3342421388178192*x1 + 3.0070132668285048;
    double x46 = x10*x2;
    double x47 = x45*x46;
    double x48 = 1.0/x3;
    double x49 = 591736.02249829285*x48;
    double x50 = x17*x49;
    double x51 = x46*(8.6684842776356383*x1 + 6.0140265336570096);
    double x52 = x46*x8;
    double x53 = x45*x52;
    double x54 = x18*x49;
    double x55 = x13*x41;
    double x56 = x36*x48;
    double x57 = x37*x48;
    double x58 = exp(x5);
    double x59 = x58 - 1;
    double x60 = 1.0/x59;
    double x61 = Debye(x5);
    double x62 = x20*x61;
    double x63 = -3*x60 + 0.95420259075969494*x62;
    double x64 = x1*x63;
    double x65 = 62737.33841263425*x20;
    double x66 = 62737.338412634235*x20;
    double x67 = exp(x25);
    double x68 = x67 - 1;
    double x69 = 1.0/x68;
    double x70 = Debye(x25);
    double x71 = T*x20*x70;
    double x72 = -3*x69 + 0.0031806753025323169*x71;
    double x73 = x1*x72;
    double x74 = x45*x63;
    double x75 = x45*x72;
    double x76 = 0.95420259075969494*x61;
    double x77 = x76*x8;
    double x78 = x58/((x59)*(x59));
    double x79 = 9.4319593000000008*x78;
    double x80 = x20*x79;
    double x81 = x48*(-9.0*x60 + 2.8626077722790848*x62) + x77 - x80;
    double x82 = 125474.6768252685*x81;
    double x83 = x20*x47;
    double x84 = 0.0031806753025323169*T*x70;
    double x85 = x8*x84;
    double x86 = x67/((x68)*(x68));
    double x87 = 2829.58779*x23*x86;
    double x88 = x20*x87;
    double x89 = x48*(-9.0*x69 + 0.0095420259075969516*x71) + x85 - x88;
    double x90 = 125474.67682526847*x89;
    double x91 = x1*x10*x43;
    double x92 = x1*x39;
    double x93 = x48*x92;
    double x94 = x8*x92;
    double x95 = 3.0*x48;
    double x96 = 3.0*x13;
    double x97 = x64*x96;
    double x98 = x10*x38;
    double x99 = x92*x95;
    double x100 = x30*x93;

    result += (127108480.15784135*x1 + 5581230.0805477835*x12*x17 + 1860410.026849261*x12*x18 + x12*x82 - x12*x90 + 3720820.053698522*x12*exp(-9.4319593000000008*x4)/((x7)*(x7)*(x7)) + x15*x17 + x15*x18 + x17*x42 + x18*x42 - 564636.04571370815*x18*x53 + x19*x22 + x19*x44 - 599336732.03155637*x2 + x21*x64*x65 - x21*x66*x73 - x22*x29 - 167436902416.43344*x28*x31 + 564636.04571370804*x28*x53 - x29*x44 - 502310707249.30029*x31*x33 - x34*x36 - x34*x37 - x36*x55 - x37*x55 + 62737.33841263425*x44*x63 - 62737.338412634235*x44*x72 + x46*x65*(x39*x97 - x45*x77 + x45*x80 - x74*x95 + x76*x91 + 29.65395207895217*x78*x93 - x79*x94 + x81*x99 + x97*x98 - 59.307904157904339*x93*exp(x16)/((x59)*(x59)*(x59))) - x46*x66*(2668855.6871056948*x100*x86 - 5337711.3742113896*x100*exp(x32)/((x68)*(x68)*(x68)) - x45*x85 + x45*x88 + x72*x92*x96 + x73*x96*x98 - x75*x95 + x84*x91 - x87*x94 + x89*x99) - x47*x50 - x47*x54 + x47*x56 + x47*x57 - x50*x51 - x51*x54 + x51*x56 + x51*x57 - 188212.01523790276*x52*x74 + 188212.0152379027*x52*x75 - x82*x83 + x83*x90 - 334873804832.86688*x31*exp(-2829.58779*x24)/((x27)*(x27)*(x27)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*4.4227;
static const double Vmax = 1.15*4.4227;
static double V = 0.9*4.4227;

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



const char *Kyanite_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *Kyanite_slb_em_coder_calib_name(void) {
    return "Kyanite_slb_em";
}

const char *Kyanite_slb_em_coder_calib_formula(void) {
    return "Al2SiO5";
}

const double Kyanite_slb_em_coder_calib_mw(void) {
    return 162.04558;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,5.0,0.0,0.0,0.0,
        0.0,2.0,1.0,0.0,0.0,0.0,
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

const double *Kyanite_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double Kyanite_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Kyanite_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Kyanite_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Kyanite_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Kyanite_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Kyanite_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Kyanite_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Kyanite_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Kyanite_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Kyanite_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Kyanite_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Kyanite_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Kyanite_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Kyanite_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Kyanite_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Kyanite_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Kyanite_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Kyanite_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Kyanite_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Kyanite_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Kyanite_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Kyanite_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Kyanite_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Kyanite_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Kyanite_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Kyanite_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Kyanite_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Kyanite_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Kyanite_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Kyanite_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

