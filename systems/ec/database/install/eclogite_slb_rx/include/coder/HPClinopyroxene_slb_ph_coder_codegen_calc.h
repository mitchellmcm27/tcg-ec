#include <math.h>


static double coder_g(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = 1.0/(n1 + n2);

result = 16.628925236306479*T*(n1*log(n1*x0) + n2*log(n2*x0)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = n2*x1;
    double x3 = n1*x1;
    double x4 = pow(x0, -2);
    double x5 = 16.628925236306479*T;

result[0] = x5*(x0*(-n1*x4 + x1) - x2 + log(x3)) + (*endmember[0].mu0)(T, P);
result[1] = x5*(x0*(-n2*x4 + x1) - x3 + log(x2)) + (*endmember[1].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -2);
    double x2 = n2*x1;
    double x3 = -2*x1;
    double x4 = 2/((x0)*(x0)*(x0));
    double x5 = n1*x4;
    double x6 = 1.0/x0;
    double x7 = n1*x1;
    double x8 = -x7;
    double x9 = x6 + x8;
    double x10 = 16.628925236306479*T;
    double x11 = -x2 + x6;

result[0] = x10*(x0*(x3 + x5) + x2 + x9 + x0*x9/n1);
result[1] = x10*(x0*(-x1 + x5) + x2 - x6 + x8);
result[2] = x10*(x0*(n2*x4 + x3) + x11 + x7 + x0*x11/n2);
}
        
static void coder_d3gdn3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -2);
    double x2 = -4*x1;
    double x3 = pow(x0, -3);
    double x4 = 6*x3;
    double x5 = 6/((x0)*(x0)*(x0)*(x0));
    double x6 = -n1*x5;
    double x7 = -2*x1;
    double x8 = 2*x3;
    double x9 = n1*x8;
    double x10 = 1.0/n1;
    double x11 = x0*x10;
    double x12 = 1.0/x0;
    double x13 = -n1*x1 + x12;
    double x14 = n2*x8;
    double x15 = 4*x3;
    double x16 = n1*x15 - x14;
    double x17 = x10*x13 + x16;
    double x18 = 16.628925236306479*T;
    double x19 = 1.0/n2;
    double x20 = -n2*x1 + x12;

result[0] = x18*(x0*(x4 + x6) + x11*(x7 + x9) + x17 + x2 - x0*x13/((n1)*(n1)));
result[1] = x18*(x0*(x15 + x6) + x11*(-x1 + x9) + x17 + x7);
result[2] = x18*(x0*(x6 + x8) + x1 + x16);
result[3] = x18*(n2*x15 + x0*x19*(x14 + x7) + x0*(-n2*x5 + x4) + x19*x20 + x2 - x9 - x0*x20/((n2)*(n2)));
}
        
static double coder_dgdt(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = 1.0/(n1 + n2);

result = 16.628925236306479*n1*log(n1*x0) + n1*(*endmember[0].dmu0dT)(T, P) + 16.628925236306479*n2*log(n2*x0) + n2*(*endmember[1].dmu0dT)(T, P);
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = n2*x1;
    double x3 = n1*x1;
    double x4 = pow(x0, -2);
    double x5 = 16.628925236306479*x0;

result[0] = -16.628925236306479*x2 + x5*(-n1*x4 + x1) + 16.628925236306479*log(x3) + (*endmember[0].dmu0dT)(T, P);
result[1] = -16.628925236306479*x3 + x5*(-n2*x4 + x1) + 16.628925236306479*log(x2) + (*endmember[1].dmu0dT)(T, P);
}
        
static void coder_d3gdn2dt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = 16.628925236306479*x1;
    double x3 = 16.628925236306479*n1;
    double x4 = 16.628925236306479*n2;
    double x5 = x3 + x4;
    double x6 = pow(x0, -2);
    double x7 = -2*x6;
    double x8 = 2/((x0)*(x0)*(x0));
    double x9 = n1*x8;
    double x10 = 16.628925236306479*x0;
    double x11 = x4*x6;
    double x12 = x3*x6;
    double x13 = x11 - x12;

result[0] = x13 + x2 + x5*(x7 + x9) + x10*(-n1*x6 + x1)/n1;
result[1] = x13 - x2 + x5*(-x6 + x9);
result[2] = -x11 + x12 + x2 + x5*(n2*x8 + x7) + x10*(-n2*x6 + x1)/n2;
}
        
static void coder_d4gdn3dt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -2);
    double x2 = -66.515700945225916*x1;
    double x3 = 16.628925236306479*n1 + 16.628925236306479*n2;
    double x4 = pow(x0, -3);
    double x5 = 6*x4;
    double x6 = 6/((x0)*(x0)*(x0)*(x0));
    double x7 = -n1*x6;
    double x8 = -2*x1;
    double x9 = 2*x4;
    double x10 = n1*x9;
    double x11 = 16.628925236306479/n1;
    double x12 = x0*x11;
    double x13 = 1.0/x0;
    double x14 = -n1*x1 + x13;
    double x15 = 16.628925236306479*x0;
    double x16 = 66.515700945225916*x4;
    double x17 = 33.257850472612958*x4;
    double x18 = n1*x16 - n2*x17;
    double x19 = x11*x14 + x18;
    double x20 = -n2*x1 + x13;
    double x21 = 16.628925236306479/n2;

result[0] = x12*(x10 + x8) + x19 + x2 + x3*(x5 + x7) - x14*x15/((n1)*(n1));
result[1] = -33.257850472612958*x1 + x12*(-x1 + x10) + x19 + x3*(4*x4 + x7);
result[2] = 16.628925236306479*x1 + x18 + x3*(x7 + x9);
result[3] = -n1*x17 + n2*x16 + x0*x21*(n2*x9 + x8) + x2 + x20*x21 + x3*(-n2*x6 + x5) - x15*x20/((n2)*(n2));
}
        
static double coder_dgdp(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].dmu0dP)(T, P);
result[1] = (*endmember[1].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d4gdn3dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d2gdt2(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P);
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = (*endmember[1].d2mu0dT2)(T, P);
}
        
static void coder_d4gdn2dt2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d5gdn3dt2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d2gdtdp(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d2mu0dTdP)(T, P) + n2*(*endmember[1].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d2mu0dTdP)(T, P);
result[1] = (*endmember[1].d2mu0dTdP)(T, P);
}
        
static void coder_d4gdn2dtdp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d5gdn3dtdp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d2gdp2(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P);
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d2mu0dP2)(T, P);
result[1] = (*endmember[1].d2mu0dP2)(T, P);
}
        
static void coder_d4gdn2dp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d5gdn3dp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d3gdt3(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P);
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d3mu0dT3)(T, P);
result[1] = (*endmember[1].d3mu0dT3)(T, P);
}
        
static void coder_d5gdn2dt3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d6gdn3dt3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d3gdt2dp(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P);
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d3mu0dT2dP)(T, P);
result[1] = (*endmember[1].d3mu0dT2dP)(T, P);
}
        
static void coder_d5gdn2dt2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d6gdn3dt2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d3gdtdp2(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d3mu0dTdP2)(T, P) + n2*(*endmember[1].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d3mu0dTdP2)(T, P);
result[1] = (*endmember[1].d3mu0dTdP2)(T, P);
}
        
static void coder_d5gdn2dtdp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d6gdn3dtdp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_d3gdp3(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P);
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = (*endmember[0].d3mu0dP3)(T, P);
result[1] = (*endmember[1].d3mu0dP3)(T, P);
}
        
static void coder_d5gdn2dp3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
}
        
static void coder_d6gdn3dp3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
}
        
static double coder_s(double T, double P, double n[2]) {
    double result = -coder_dgdt(T, P, n);
    return result;
}

static double coder_v(double T, double P, double n[2]) {
    double result = coder_dgdp(T, P, n);
    return result;
}

static double coder_cv(double T, double P, double n[2]) {
    double result = -T*coder_d2gdt2(T, P, n);
    double dvdt = coder_d2gdtdp(T, P, n);
    double dvdp = coder_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P, double n[2]) {
    double result = -T*coder_d2gdt2(T, P, n);
    return result;
}

static double coder_dcpdt(double T, double P, double n[2]) {
    double result = -T*coder_d3gdt3(T, P, n) - coder_d2gdt2(T, P, n);
    return result;
}

static double coder_alpha(double T, double P, double n[2]) {
    double result = coder_d2gdtdp(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_beta(double T, double P, double n[2]) {
    double result = -coder_d2gdp2(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_K(double T, double P, double n[2]) {
    double result = -coder_dgdp(T, P, n)/coder_d2gdp2(T, P, n);
    return result;
}

static double coder_Kp(double T, double P, double n[2]) {
    double result = coder_dgdp(T, P, n);
    result *= coder_d3gdp3(T, P, n);
    result /= pow(coder_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}

