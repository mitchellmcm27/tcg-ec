#include <math.h>


static double coder_g(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].mu0)(T, P);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d3gdn3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_dgdt(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].dmu0dT)(T, P);
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].dmu0dT)(T, P);
}
        
static void coder_d3gdn2dt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d4gdn3dt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_dgdp(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d4gdn3dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d2gdt2(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d2mu0dT2)(T, P);
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d2mu0dT2)(T, P);
}
        
static void coder_d4gdn2dt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d5gdn3dt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d2gdtdp(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d2mu0dTdP)(T, P);
}
        
static void coder_d4gdn2dtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d5gdn3dtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d2gdp2(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d2mu0dP2)(T, P);
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d2mu0dP2)(T, P);
}
        
static void coder_d4gdn2dp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d5gdn3dp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d3gdt3(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d3mu0dT3)(T, P);
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d3mu0dT3)(T, P);
}
        
static void coder_d5gdn2dt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d6gdn3dt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d3gdt2dp(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d3mu0dT2dP)(T, P);
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d3mu0dT2dP)(T, P);
}
        
static void coder_d5gdn2dt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d6gdn3dt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d3gdtdp2(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d3mu0dTdP2)(T, P);
}
        
static void coder_d5gdn2dtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d6gdn3dtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_d3gdp3(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    

result = n1*(*endmember[0].d3mu0dP3)(T, P);
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = (*endmember[0].d3mu0dP3)(T, P);
}
        
static void coder_d5gdn2dp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static void coder_d6gdn3dp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];


result[0] = 0;
}
        
static double coder_s(double T, double P, double n[1]) {
    double result = -coder_dgdt(T, P, n);
    return result;
}

static double coder_v(double T, double P, double n[1]) {
    double result = coder_dgdp(T, P, n);
    return result;
}

static double coder_cv(double T, double P, double n[1]) {
    double result = -T*coder_d2gdt2(T, P, n);
    double dvdt = coder_d2gdtdp(T, P, n);
    double dvdp = coder_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P, double n[1]) {
    double result = -T*coder_d2gdt2(T, P, n);
    return result;
}

static double coder_dcpdt(double T, double P, double n[1]) {
    double result = -T*coder_d3gdt3(T, P, n) - coder_d2gdt2(T, P, n);
    return result;
}

static double coder_alpha(double T, double P, double n[1]) {
    double result = coder_d2gdtdp(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_beta(double T, double P, double n[1]) {
    double result = -coder_d2gdp2(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_K(double T, double P, double n[1]) {
    double result = -coder_dgdp(T, P, n)/coder_d2gdp2(T, P, n);
    return result;
}

static double coder_Kp(double T, double P, double n[1]) {
    double result = coder_dgdp(T, P, n);
    result *= coder_d3gdp3(T, P, n);
    result /= pow(coder_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}

