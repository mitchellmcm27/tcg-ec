
static char *identifier = "Fayalite_berman.emml:917b1867bedab06304280d5bd0afb555fa910222:Wed Sep 28 15:58:31 2022";



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

static double coder_g(double T, double P) {
    double result = 0.0;
    double x0 = ((T)*(T));
    double x1 = 3.6576999999999997e-8*x0;

    result += -1.68995e-6*((P)*(P)) + 0.00010109572589999999*P*T + P*x1 + 4.596610234508133*P - 7695.4200000000001*sqrt(T) - 248.92812000000001*T*log(T) + 1740.8753388763112*T - x1 - 1487926.5584324801 + 23184016.0/x0;
    return result;
}

static double coder_dgdt(double T, double P) {
    double result = 0.0;
    double x0 = 7.3153999999999995e-8*T;

    result += P*x0 + 0.00010109572589999999*P - x0 - 248.92812000000001*log(T) + 1491.9472188763111 - 46368032.0/((T)*(T)*(T)) - 3847.71/sqrt(T);
    return result;
}

static double coder_dgdp(double T, double P) {
    double result = 0.0;

    result += -3.3799e-6*P + 3.6576999999999997e-8*((T)*(T)) + 0.00010109572589999999*T + 4.596610234508133;
    return result;
}

static double coder_d2gdt2(double T, double P) {
    double result = 0.0;

    result += 7.3153999999999995e-8*P - 7.3153999999999995e-8 - 248.92812000000001/T + 139104096.0/((T)*(T)*(T)*(T)) + 1923.855/pow(T, 3.0/2.0);
    return result;
}

static double coder_d2gdtdp(double T, double P) {
    double result = 0.0;

    result += 7.3153999999999995e-8*T + 0.00010109572589999999;
    return result;
}

static double coder_d2gdp2(double T, double P) {
    double result = 0.0;

    result += -3.3799e-6;
    return result;
}

static double coder_d3gdt3(double T, double P) {
    double result = 0.0;

    result += 248.92812000000001/((T)*(T)) - 556416384.0/pow(T, 5) - 2885.7825000000003/pow(T, 5.0/2.0);
    return result;
}

static double coder_d3gdt2dp(double T, double P) {
    double result = 0.0;

    result += 7.3153999999999995e-8;
    return result;
}

static double coder_d3gdtdp2(double T, double P) {
    double result = 0.0;
    result += 0.0;
    return result;
}

static double coder_d3gdp3(double T, double P) {
    double result = 0.0;
    result += 0.0;
    return result;
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

static double coder_dparam_g(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_dgdt(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_dgdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2gdt2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2gdtdp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d2gdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3gdt3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3gdt2dp(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3gdtdp2(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
}

static double coder_dparam_d3gdp3(double T, double P, int index) {
    double result = 0.0;
    switch (index) {
    }
    return result;
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



const char *Fayalite_berman_coder_calib_identifier(void) {
    return identifier;
}

const char *Fayalite_berman_coder_calib_name(void) {
    return "Fayalite_berman";
}

const char *Fayalite_berman_coder_calib_formula(void) {
    return "Fe2SiO4";
}

const double Fayalite_berman_coder_calib_mw(void) {
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

const double *Fayalite_berman_coder_calib_elements(void) {
    return elmformula;
}

double Fayalite_berman_coder_calib_g(double T, double P) {
    return coder_g(T, P);
}

double Fayalite_berman_coder_calib_dgdt(double T, double P) {
    return coder_dgdt(T, P);
}

double Fayalite_berman_coder_calib_dgdp(double T, double P) {
    return coder_dgdp(T, P);
}

double Fayalite_berman_coder_calib_d2gdt2(double T, double P) {
    return coder_d2gdt2(T, P);
}

double Fayalite_berman_coder_calib_d2gdtdp(double T, double P) {
    return coder_d2gdtdp(T, P);
}

double Fayalite_berman_coder_calib_d2gdp2(double T, double P) {
    return coder_d2gdp2(T, P);
}

double Fayalite_berman_coder_calib_d3gdt3(double T, double P) {
    return coder_d3gdt3(T, P);
}

double Fayalite_berman_coder_calib_d3gdt2dp(double T, double P) {
    return coder_d3gdt2dp(T, P);
}

double Fayalite_berman_coder_calib_d3gdtdp2(double T, double P) {
    return coder_d3gdtdp2(T, P);
}

double Fayalite_berman_coder_calib_d3gdp3(double T, double P) {
    return coder_d3gdp3(T, P);
}

double Fayalite_berman_coder_calib_s(double T, double P) {
    return coder_s(T, P);
}

double Fayalite_berman_coder_calib_v(double T, double P) {
    return coder_v(T, P);
}

double Fayalite_berman_coder_calib_cv(double T, double P) {
    return coder_cv(T, P);
}

double Fayalite_berman_coder_calib_cp(double T, double P) {
    return coder_cp(T, P);
}

double Fayalite_berman_coder_calib_dcpdt(double T, double P) {
    return coder_dcpdt(T, P);
}

double Fayalite_berman_coder_calib_alpha(double T, double P) {
    return coder_alpha(T, P);
}

double Fayalite_berman_coder_calib_beta(double T, double P) {
    return coder_beta(T, P);
}

double Fayalite_berman_coder_calib_K(double T, double P) {
    return coder_K(T, P);
}

double Fayalite_berman_coder_calib_Kp(double T, double P) {
    return coder_Kp(T, P);
}

int Fayalite_berman_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Fayalite_berman_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Fayalite_berman_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Fayalite_berman_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Fayalite_berman_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Fayalite_berman_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Fayalite_berman_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Fayalite_berman_coder_dparam_g(double T, double P, int index) {
    return coder_dparam_g(T, P, index);
}

double Fayalite_berman_coder_dparam_dgdt(double T, double P, int index) {
    return coder_dparam_dgdt(T, P, index);
}

double Fayalite_berman_coder_dparam_dgdp(double T, double P, int index) {
    return coder_dparam_dgdp(T, P, index);
}

double Fayalite_berman_coder_dparam_d2gdt2(double T, double P, int index) {
    return coder_dparam_d2gdt2(T, P, index);
}

double Fayalite_berman_coder_dparam_d2gdtdp(double T, double P, int index) {
    return coder_dparam_d2gdtdp(T, P, index);
}

double Fayalite_berman_coder_dparam_d2gdp2(double T, double P, int index) {
    return coder_dparam_d2gdp2(T, P, index);
}

double Fayalite_berman_coder_dparam_d3gdt3(double T, double P, int index) {
    return coder_dparam_d3gdt3(T, P, index);
}

double Fayalite_berman_coder_dparam_d3gdt2dp(double T, double P, int index) {
    return coder_dparam_d3gdt2dp(T, P, index);
}

double Fayalite_berman_coder_dparam_d3gdtdp2(double T, double P, int index) {
    return coder_dparam_d3gdtdp2(T, P, index);
}

double Fayalite_berman_coder_dparam_d3gdp3(double T, double P, int index) {
    return coder_dparam_d3gdp3(T, P, index);
}

