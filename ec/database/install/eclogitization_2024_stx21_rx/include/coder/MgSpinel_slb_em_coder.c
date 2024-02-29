
static char *identifier = "MgSpinel_slb_em.emml:5ad9bb808635ee819fe108066c5a759660a54c51:Tue Aug 22 19:34:06 2023";



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
    double x3 = sqrt(95.298033062983507*x1 - 242.84329249264871*x2 - 7.9969351035500011);
    double x4 = 2.854565*x3;
    double x5 = 856.36950000000002*x3/T;

    result += 698.41485992487208*T*log(1 - exp(-x5)) - 232.80495330829069*T*Debye(x5) - 234097799.85697436*x1 + 83954892.936432123*x2 - 209524.45797746163*log(1 - exp(-x4)) + 69841.485992487214*Debye(x4) + 15297006.911863204 + 2767167400.8351054/((V)*(V));
    return result;
}

static double coder_dadt(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = sqrt(-242.84329249264871*pow(x0, 4.0/3.0) + 95.298033062983507*pow(x0, 2.0/3.0) - 7.9969351035500011);
    double x2 = x1/T;
    double x3 = 856.36950000000002*x2;
    double x4 = Debye(x3);
    double x5 = exp(-x3);
    double x6 = -x5 + 1;

    result += -598101.18438643275*x2*x5/x6 + 199367.06146214425*x2*(-0.0035031607267657239*T*x4/x1 + 3/(exp(x3) - 1)) - 232.80495330829069*x4 + 698.41485992487208*log(x6);
    return result;
}

static double coder_dadv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/V;
    double x1 = pow(x0, 2.0/3.0);
    double x2 = x0*x1;
    double x3 = pow(x0, 4.0/3.0);
    double x4 = x0*x3;
    double x5 = sqrt(95.298033062983507*x1 - 242.84329249264871*x3 - 7.9969351035500011);
    double x6 = 2.854565*x5;
    double x7 = exp(-x6);
    double x8 = 1.0/x5;
    double x9 = x8*(-31.766011020994501*x2 + 161.89552832843248*x4);
    double x10 = 598101.18438643275*x9;
    double x11 = 856.36950000000002*x5/T;
    double x12 = exp(-x11);

    result += x10*x12/(-x12 + 1) - x10*x7/(-x7 + 1) + 156065199.90464956*x2 - 111939857.24857616*x4 + 199367.06146214428*x9*(-1.0509482180297174*x8*Debye(x6) + 3/(exp(x6) - 1)) - 199367.06146214425*x9*(-0.0035031607267657239*T*x8*Debye(x11) + 3/(exp(x11) - 1)) - 5534334801.6702108/((V)*(V)*(V));
    return result;
}

static double coder_d2adt2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = sqrt(95.298033062983507*x2 - 242.84329249264871*x3 - 7.9969351035500011);
    double x5 = x0*x4;
    double x6 = 856.36950000000002*x5;
    double x7 = exp(-x6);
    double x8 = -x7 + 1;
    double x9 = pow(T, -2);
    double x10 = x9*(-48811234388.286995*x2 + 124383268872.37975*x3 + 4095995071.2657323);
    double x11 = Debye(x6)/x4;
    double x12 = exp(x6);
    double x13 = x12 - 1;

    result += x0*(x10*x7/x8 + x10*exp(-1712.739*x5)/((x8)*(x8)) - 199367.06146214425*x4*(x0*(0.01050948218029717*T*x11 - 8.9999999999999982/x13) + 0.0035031607267657239*x11 - 2569.1085000000003*x12*x4*x9/((x13)*(x13))));
    return result;
}

static double coder_d2adtdv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = 95.298033062983507*x2;
    double x4 = 242.84329249264871*pow(x1, 4.0/3.0);
    double x5 = x3 - x4 - 7.9969351035500011;
    double x6 = sqrt(x5);
    double x7 = x0*x6;
    double x8 = 856.36950000000002*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = 512195612.22241724*x0;
    double x12 = T*Debye(x8);
    double x13 = 1.0/x6;
    double x14 = exp(x8);
    double x15 = x14 - 1;

    result += x0*x1*x2*(161.89552832843248*x2 - 31.766011020994501)*(-199367.06146214425*x6*(2569.1085000000003*x0*x13*x14/((x15)*(x15)) - 0.0035031607267657239*x12/pow(x5, 3.0/2.0) + (0.01050948218029717*x12*x13 - 8.9999999999999982/x15)/(-x3 + x4 + 7.9969351035500011)) + x11*x9/x10 + x11*exp(-1712.739*x7)/((x10)*(x10)));
    return result;
}

static double coder_d2adv2(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 95.298033062983507*x2;
    double x5 = 242.84329249264871*x3;
    double x6 = x4 - x5 - 7.9969351035500011;
    double x7 = sqrt(x6);
    double x8 = 2.854565*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = x9/x10;
    double x12 = 1.0/x7;
    double x13 = x2*(377.75623276634246*x2 - 52.943351701657505);
    double x14 = x12*x13;
    double x15 = 598101.18438643275*x14;
    double x16 = 1.0/(-x4 + x5 + 7.9969351035500011);
    double x17 = x3*((161.89552832843248*x2 - 31.766011020994501)*(161.89552832843248*x2 - 31.766011020994501));
    double x18 = x16*x17;
    double x19 = 1707318.7074080575*x18;
    double x20 = pow(x6, -3.0/2.0);
    double x21 = x17*x20;
    double x22 = 598101.18438643275*x21;
    double x23 = 1.0/T;
    double x24 = x23*x7;
    double x25 = 856.36950000000002*x24;
    double x26 = exp(-x25);
    double x27 = -x26 + 1;
    double x28 = x26/x27;
    double x29 = 512195612.22241724*x18*x23;
    double x30 = exp(x8);
    double x31 = x30 - 1;
    double x32 = 1.0/x31;
    double x33 = Debye(x8);
    double x34 = x12*x33;
    double x35 = -598101.18438643287*x32 + 209524.45797746169*x34;
    double x36 = exp(x25);
    double x37 = x36 - 1;
    double x38 = 1.0/x37;
    double x39 = T*Debye(x25);
    double x40 = x12*x39;
    double x41 = -3*x38 + 0.0035031607267657239*x40;
    double x42 = 199367.06146214425*x12;

    result += x0*(16603004405.010632*x0 + x11*x15 - x11*x19 + x11*x22 - 199367.06146214428*x12*x17*(8.5636949999999992*x12*x30/((x31)*(x31)) + x16*(-9.0000000000000018*x32 + 3.1528446540891526*x34) - 1.0509482180297174*x20*x33) - x13*x41*x42 + x14*x35 - x15*x28 + x17*x42*(2569.1085000000003*x12*x23*x36/((x37)*(x37)) + x16*(-8.9999999999999982*x38 + 0.01050948218029717*x40) - 0.0035031607267657239*x20*x39) - 260108666.50774926*x2 + x21*x35 - 199367.06146214425*x21*x41 - x22*x28 + x28*x29 + 261193000.2466777*x3 + x29*exp(-1712.739*x24)/((x27)*(x27)) - x19*exp(-5.70913*x7)/((x10)*(x10)));
    return result;
}

static double coder_d3adt3(double T, double V) {
    double result = 0.0;
    double x0 = pow(T, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = x0*(-146433703164.86099*x2 + 373149806617.13922*x3 + 12287985213.797195);
    double x5 = 1.0/T;
    double x6 = 95.298033062983507*x2;
    double x7 = 242.84329249264871*x3;
    double x8 = x6 - x7 - 7.9969351035500011;
    double x9 = sqrt(x8);
    double x10 = x5*x9;
    double x11 = 856.36950000000002*x10;
    double x12 = exp(-x11);
    double x13 = -x12 + 1;
    double x14 = 1712.739*x10;
    double x15 = exp(-x14)/((x13)*(x13));
    double x16 = x12/x13;
    double x17 = pow(T, -3);
    double x18 = x17*pow(x8, 3.0/2.0);
    double x19 = Debye(x11)/x9;
    double x20 = exp(x11);
    double x21 = x20 - 1;
    double x22 = x20/((x21)*(x21));
    double x23 = x0*x22*x9;
    double x24 = 0.01050948218029717*x19;
    double x25 = x5*(T*x24 - 8.9999999999999982/x21);
    double x26 = 199367.06146214425*x9;
    double x27 = x17*(-x6 + x7 + 7.9969351035500011);

    result += x0*(-1315886101023.3159*x15*x18 - x15*x4 - 438628700341.10535*x16*x18 - x16*x4 + x26*(0.0035031607267657239*x19 - 2569.1085000000003*x23 + x25) - x26*(-2200106.1615907503*x22*x27 - 2569.1084999999994*x23 + x24 + 2.9999999999999996*x25 + 4400212.3231815007*x27*exp(x14)/((x21)*(x21)*(x21))) - 877257400682.21069*x18*exp(-2569.1085000000003*x10)/((x13)*(x13)*(x13)));
    return result;
}

static double coder_d3adt2dv(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(T, -2);
    double x4 = x3*(165844358496.50632*x2 - 32540822925.524662);
    double x5 = 95.298033062983507*x2;
    double x6 = 242.84329249264871*pow(x1, 4.0/3.0);
    double x7 = x5 - x6 - 7.9969351035500011;
    double x8 = sqrt(x7);
    double x9 = x0*x8;
    double x10 = 856.36950000000002*x9;
    double x11 = exp(-x10);
    double x12 = -x11 + 1;
    double x13 = 1712.739*x9;
    double x14 = exp(-x13)/((x12)*(x12));
    double x15 = x11/x12;
    double x16 = pow(T, -3);
    double x17 = 161.89552832843248*x2 - 31.766011020994501;
    double x18 = x16*x17*x8;
    double x19 = 1.0/x8;
    double x20 = Debye(x10);
    double x21 = 0.0035031607267657239*x20;
    double x22 = exp(x10);
    double x23 = x22 - 1;
    double x24 = x22/((x23)*(x23));
    double x25 = 2569.1085000000003*x3;
    double x26 = 0.01050948218029717*T*x20;
    double x27 = x19*x26 - 8.9999999999999982/x23;
    double x28 = x0*x27;
    double x29 = 199367.06146214425*x17;
    double x30 = pow(x7, -3.0/2.0);
    double x31 = x19*x24;
    double x32 = 1.0/(-x5 + x6 + 7.9969351035500011);

    result += x0*x1*x2*(1315886101023.3159*x14*x18 - x14*x4 + 438628700341.10535*x15*x18 - x15*x4 - x19*x29*(x19*x21 - x24*x25*x8 + x28) + x29*x8*(-x0*(7707.3254999999999*x0*x31 - x26*x30 + 2.9999999999999996*x27*x32) + 2200106.1615907503*x16*x24 - 4400212.3231815007*x16*exp(x13)/((x23)*(x23)*(x23)) + x21*x30 + x25*x31 - x28*x32) + 877257400682.21069*x18*exp(-2569.1085000000003*x9)/((x12)*(x12)*(x12)));
    return result;
}

static double coder_d3adtdv2(double T, double V) {
    double result = 0.0;
    double x0 = 1.0/T;
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = x0*(193485084912.59073*x2 - 27117352437.937222);
    double x4 = 95.298033062983507*x2;
    double x5 = 242.84329249264871*pow(x1, 4.0/3.0);
    double x6 = x4 - x5 - 7.9969351035500011;
    double x7 = sqrt(x6);
    double x8 = x0*x7;
    double x9 = 856.36950000000002*x8;
    double x10 = exp(-x9);
    double x11 = -x10 + 1;
    double x12 = 1712.739*x8;
    double x13 = exp(-x12)/((x11)*(x11));
    double x14 = x10/x11;
    double x15 = 1.0/x7;
    double x16 = 161.89552832843248*x2 - 31.766011020994501;
    double x17 = ((x16)*(x16));
    double x18 = x17*x2;
    double x19 = x18/((T)*(T));
    double x20 = x15*x19;
    double x21 = pow(x6, -3.0/2.0);
    double x22 = Debye(x9);
    double x23 = 0.0035031607267657239*T*x22;
    double x24 = x21*x23;
    double x25 = exp(x9);
    double x26 = x25 - 1;
    double x27 = x25/((x26)*(x26));
    double x28 = 2569.1085000000003*x27;
    double x29 = x0*x15*x28;
    double x30 = -x4 + x5 + 7.9969351035500011;
    double x31 = 1.0/x30;
    double x32 = 1.0/x26;
    double x33 = T*x15*x22;
    double x34 = x18*(-x24 + x29 + x31*(-8.9999999999999982*x32 + 0.01050948218029717*x33));
    double x35 = 377.75623276634246*x2 - 52.943351701657505;
    double x36 = x16*x2;
    double x37 = x19*x31;
    double x38 = x0*x2;
    double x39 = -3*x32 + 0.0035031607267657239*x33;
    double x40 = 2.9999999999999996*x31;
    double x41 = 2.9999999999999996*x39/((x30)*(x30));

    result += -x38*(1315886101023.3159*x13*x20 + x13*x3 + 438628700341.10535*x14*x20 + x14*x3 + 199367.06146214425*x15*x34 + 199367.06146214425*x7*(-x17*x21*x28*x38 + x18*x41 + x23*x36*(485.68658498529743*x2 - 95.298033062983507)/pow(x6, 5.0/2.0) + x24*x35 - 2200106.1615907503*x27*x37 - x29*x35 + x34*x40 - x35*x39*x40 + x36*x41*(323.79105665686495*x2 - 63.532022041989002) + 4400212.3231815007*x37*exp(x12)/((x26)*(x26)*(x26))) + 877257400682.21069*x20*exp(-2569.1085000000003*x8)/((x11)*(x11)*(x11)))/((V)*(V));
    return result;
}

static double coder_d3adv3(double T, double V) {
    double result = 0.0;
    double x0 = pow(V, -2);
    double x1 = 1.0/V;
    double x2 = pow(x1, 2.0/3.0);
    double x3 = pow(x1, 4.0/3.0);
    double x4 = 95.298033062983507*x2;
    double x5 = 242.84329249264871*x3;
    double x6 = x4 - x5 - 7.9969351035500011;
    double x7 = sqrt(x6);
    double x8 = 2.854565*x7;
    double x9 = exp(-x8);
    double x10 = -x9 + 1;
    double x11 = pow(x6, -3.0/2.0);
    double x12 = 161.89552832843248*x2 - 31.766011020994501;
    double x13 = x0*((x12)*(x12)*(x12));
    double x14 = x11*x13;
    double x15 = 5.70913*x7;
    double x16 = exp(-x15)/((x10)*(x10));
    double x17 = -x4 + x5 + 7.9969351035500011;
    double x18 = pow(x17, -2);
    double x19 = x13*x18;
    double x20 = 1707318.7074080575*x16;
    double x21 = x9/x10;
    double x22 = 1707318.7074080575*x21;
    double x23 = 1.0/x7;
    double x24 = x2*(1259.1874425544747*x2 - 141.18227120442);
    double x25 = 598101.18438643275*x23*x24;
    double x26 = 1.0/T;
    double x27 = x26*x7;
    double x28 = 856.36950000000002*x27;
    double x29 = exp(-x28);
    double x30 = -x29 + 1;
    double x31 = x29/x30;
    double x32 = pow(T, -2);
    double x33 = x14*x32;
    double x34 = 1712.739*x27;
    double x35 = exp(-x34)/((x30)*(x30));
    double x36 = 512195612.22241724*x26;
    double x37 = x19*x36;
    double x38 = x18*(323.79105665686495*x2 - 63.532022041989002);
    double x39 = ((x12)*(x12));
    double x40 = x0*x39;
    double x41 = x38*x40;
    double x42 = (485.68658498529743*x2 - 95.298033062983507)/pow(x6, 5.0/2.0);
    double x43 = x40*x42;
    double x44 = 598101.18438643275*x43;
    double x45 = 1.0/x17;
    double x46 = 377.75623276634246*x2 - 52.943351701657505;
    double x47 = x45*x46;
    double x48 = x12*x3;
    double x49 = x20*x48;
    double x50 = x45*(755.51246553268493*x2 - 105.88670340331501);
    double x51 = x46*x48;
    double x52 = x11*x51;
    double x53 = 1794303.5531592984*x52;
    double x54 = x22*x48;
    double x55 = x36*x41;
    double x56 = x36*x48;
    double x57 = x35*x56;
    double x58 = x31*x56;
    double x59 = exp(x8);
    double x60 = x59 - 1;
    double x61 = 1.0/x60;
    double x62 = Debye(x8);
    double x63 = x23*x62;
    double x64 = -3*x61 + 1.0509482180297174*x63;
    double x65 = 199367.06146214428*x23;
    double x66 = exp(x28);
    double x67 = x66 - 1;
    double x68 = 1.0/x67;
    double x69 = Debye(x28);
    double x70 = T*x23*x69;
    double x71 = -3*x68 + 0.0035031607267657239*x70;
    double x72 = 199367.06146214425*x23;
    double x73 = 1.0509482180297174*x62;
    double x74 = x11*x73;
    double x75 = x59/((x60)*(x60));
    double x76 = 8.5636949999999992*x75;
    double x77 = x23*x76;
    double x78 = x45*(-9.0000000000000018*x61 + 3.1528446540891526*x63) - x74 + x77;
    double x79 = 398734.12292428856*x78;
    double x80 = x23*x51;
    double x81 = 0.0035031607267657239*T*x69;
    double x82 = x11*x81;
    double x83 = x66/((x67)*(x67));
    double x84 = 2569.1085000000003*x26*x83;
    double x85 = x23*x84;
    double x86 = x45*(-8.9999999999999982*x68 + 0.01050948218029717*x70) - x82 + x85;
    double x87 = 398734.1229242885*x86;
    double x88 = x12*x2;
    double x89 = x42*x88;
    double x90 = x2*x39;
    double x91 = x45*x90;
    double x92 = x11*x90;
    double x93 = 3.0000000000000004*x64;
    double x94 = x18*x90;
    double x95 = x38*x88;
    double x96 = x32*x91;
    double x97 = 2.9999999999999996*x71;

    result += (-66412017620.042526*x0 - 14620956.678036846*x14*x16 - 4873652.2260122821*x14*x21 + x14*x79 - x14*x87 - x19*x20 - x19*x22 + 693623110.68733132*x2 - x20*x41 - x21*x25 - x21*x44 - x21*x53 - x22*x41 - x24*x64*x65 + x24*x71*x72 + x25*x31 - 870643334.15559232*x3 + 438628700341.10535*x31*x33 + x31*x37 + x31*x44 + x31*x53 + x31*x55 + 1315886101023.3159*x33*x35 + x35*x37 + x35*x55 - 199367.06146214428*x43*x64 + 199367.06146214425*x43*x71 + x47*x49 + x47*x54 - x47*x57 - x47*x58 - x48*x65*(x46*x74 - x46*x77 - x47*x93 + x73*x89 - 24.445624017674998*x75*x91 - x76*x92 + 3.0000000000000004*x78*x91 + x93*x94 + x93*x95 + 48.891248035349996*x91*exp(x15)/((x60)*(x60)*(x60))) + x48*x72*(x46*x82 - x46*x85 - x47*x97 + x81*x89 - 2200106.1615907503*x83*x96 - x84*x92 + 2.9999999999999996*x86*x91 + x94*x97 + x95*x97 + 4400212.3231815007*x96*exp(x34)/((x67)*(x67)*(x67))) + x49*x50 + x50*x54 - x50*x57 - x50*x58 - 598101.18438643287*x52*x64 + 598101.18438643275*x52*x71 + x79*x80 - x80*x87 + 877257400682.21069*x33*exp(-2569.1085000000003*x27)/((x30)*(x30)*(x30)) - 9747304.4520245641*x14*exp(-8.5636949999999992*x7)/((x10)*(x10)*(x10)))/((V)*(V)*(V));
    return result;
}


static const double tol = 0.001;
static const int MAX_ITS = 200;
static double Told = 0.0;
static double Pold = 0.0;
static const double Vmin = 0.5*15.9048;
static const double Vmax = 1.15*15.9048;
static double V = 0.9*15.9048;

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



const char *MgSpinel_slb_em_coder_calib_identifier(void) {
    return identifier;
}

const char *MgSpinel_slb_em_coder_calib_name(void) {
    return "MgSpinel_slb_em";
}

const char *MgSpinel_slb_em_coder_calib_formula(void) {
    return "Mg4Al8O16";
}

const double MgSpinel_slb_em_coder_calib_mw(void) {
    return 569.06272;
}

static const double elmformula[106] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,16.0,0.0,0.0,0.0,
        4.0,8.0,0.0,0.0,0.0,0.0,
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

const double *MgSpinel_slb_em_coder_calib_elements(void) {
    return elmformula;
}

double MgSpinel_slb_em_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double MgSpinel_slb_em_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double MgSpinel_slb_em_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double MgSpinel_slb_em_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double MgSpinel_slb_em_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double MgSpinel_slb_em_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double MgSpinel_slb_em_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double MgSpinel_slb_em_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double MgSpinel_slb_em_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double MgSpinel_slb_em_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double MgSpinel_slb_em_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double MgSpinel_slb_em_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double MgSpinel_slb_em_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double MgSpinel_slb_em_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int MgSpinel_slb_em_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **MgSpinel_slb_em_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **MgSpinel_slb_em_coder_get_param_units(void) {
    return coder_get_param_units();
}

void MgSpinel_slb_em_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int MgSpinel_slb_em_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double MgSpinel_slb_em_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int MgSpinel_slb_em_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double MgSpinel_slb_em_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double MgSpinel_slb_em_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

