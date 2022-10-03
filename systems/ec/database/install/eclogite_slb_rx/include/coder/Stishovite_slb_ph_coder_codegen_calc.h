#include <math.h>


static double coder_g(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    
    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = 0.00030000000000000003*P;
    double x2 = 0.0001*P - 51.0;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + sqrt((-0.012*T + x2)/x2)*(-0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x1 + 153.00000000000003) + 2601.0)/(x1 - 153.0));
}
    return result;
}
        
static void coder_dgdn(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = 0.00030000000000000003*P;
    double x2 = 0.0001*P - 51.0;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + sqrt((-0.012*T + x2)/x2)*(-0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x1 + 153.00000000000003) + 2601.0)/(x1 - 153.0);
}
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
    
    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = 0.0001*P - 51.0;
    double x2 = -0.012*T + x1;
    double x3 = 0.00030000000000000003*P;
    double x4 = sqrt(x2/x1)/(x3 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + x4*(3.6000000000000007e-6*P - 1.2240000000000002) - 0.0060000000000000001*x4*(-0.0051000000000000004*P + 0.61199999999999999*T + x1*(0.036000000000000004*T - x3 + 153.00000000000003) + 2601.0)/x2);
}
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = 0.0001*P - 51.0;
    double x2 = -0.012*T + x1;
    double x3 = 0.00030000000000000003*P;
    double x4 = sqrt(x2/x1)/(x3 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + x4*(3.6000000000000007e-6*P - 1.2240000000000002) - 0.0060000000000000001*x4*(-0.0051000000000000004*P + 0.61199999999999999*T + x1*(0.036000000000000004*T - x3 + 153.00000000000003) + 2601.0)/x2;
}
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
    
    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 0.0001*P - 51.0;
    double x2 = 1.0/x1;
    double x3 = -0.012*T + x1;
    double x4 = sqrt(x2*x3);
    double x5 = 0.00030000000000000003*P;
    double x6 = x5 - 153.0;
    double x7 = x4/x6;
    double x8 = -0.0051000000000000004*P + 0.61199999999999999*T + x1*(0.036000000000000004*T - x5 + 153.00000000000003) + 2601.0;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + x1*x7*x8*(5.0000000000000002e-5*x2 - 5.0000000000000002e-5*x3/((x1)*(x1)))/x3 - 0.00030000000000000003*x4*x8/((x6)*(x6)) + x7*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002));
}
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 0.0001*P - 51.0;
    double x2 = 1.0/x1;
    double x3 = -0.012*T + x1;
    double x4 = sqrt(x2*x3);
    double x5 = 0.00030000000000000003*P;
    double x6 = x5 - 153.0;
    double x7 = x4/x6;
    double x8 = -0.0051000000000000004*P + 0.61199999999999999*T + x1*(0.036000000000000004*T - x5 + 153.00000000000003) + 2601.0;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + x1*x7*x8*(5.0000000000000002e-5*x2 - 5.0000000000000002e-5*x3/((x1)*(x1)))/x3 - 0.00030000000000000003*x4*x8/((x6)*(x6)) + x7*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002);
}
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
    
    double x0 = (*endmember[0].d2mu0dT2)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = x1 - 51.0;
    double x4 = 0.00030000000000000003*P;
    double x5 = sqrt(-x2/x3)/(x4 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = -n1*(-x0 - 0.012*x5*(3.6000000000000007e-6*P - 1.2240000000000002)/x2 + 3.6000000000000001e-5*x5*(-0.0051000000000000004*P + 0.61199999999999999*T + x3*(0.036000000000000004*T - x4 + 153.00000000000003) + 2601.0)/((x2)*(x2)));
}
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dT2)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = x1 - 51.0;
    double x4 = 0.00030000000000000003*P;
    double x5 = sqrt(-x2/x3)/(x4 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.012*x5*(3.6000000000000007e-6*P - 1.2240000000000002)/x2 - 3.6000000000000001e-5*x5*(-0.0051000000000000004*P + 0.61199999999999999*T + x3*(0.036000000000000004*T - x4 + 153.00000000000003) + 2601.0)/((x2)*(x2));
}
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
    
    double x0 = (*endmember[0].d2mu0dTdP)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 0.012*T - x1 + 51.0;
    double x4 = x3/x2;
    double x5 = sqrt(-x4);
    double x6 = 0.00030000000000000003*P;
    double x7 = x6 - 153.0;
    double x8 = x5/x7;
    double x9 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x10 = x5/((x7)*(x7));
    double x11 = 1.0/x3;
    double x12 = x11*x8;
    double x13 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x6 + 153.00000000000003) + 2601.0;
    double x14 = x13*x8/((x3)*(x3));
    double x15 = x4 + 1;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = -n1*(-x0 + 1.8000000000000001e-6*x10*x11*x13 + 0.00030000000000000003*x10*x9 + 5.0000000000000002e-5*x12*x15*x9 - 0.0060000000000000001*x12*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002) + 3.0000000000000004e-7*x14*x15 - 6.0000000000000008e-7*x14 - 3.6000000000000007e-6*x8);
}
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dTdP)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 0.012*T - x1 + 51.0;
    double x4 = x3/x2;
    double x5 = sqrt(-x4);
    double x6 = 0.00030000000000000003*P;
    double x7 = x6 - 153.0;
    double x8 = x5/x7;
    double x9 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x10 = x5/((x7)*(x7));
    double x11 = 1.0/x3;
    double x12 = x11*x8;
    double x13 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x6 + 153.00000000000003) + 2601.0;
    double x14 = x13*x8/((x3)*(x3));
    double x15 = x4 + 1;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 1.8000000000000001e-6*x10*x11*x13 - 0.00030000000000000003*x10*x9 - 5.0000000000000002e-5*x12*x15*x9 + 0.0060000000000000001*x12*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002) - 3.0000000000000004e-7*x14*x15 + 6.0000000000000008e-7*x14 + 3.6000000000000007e-6*x8;
}
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
    
    double x0 = (*endmember[0].d2mu0dP2)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/x8;
    double x10 = pow(x8, -2);
    double x11 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x12 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x13 = x12*x6;
    double x14 = x5 + 1;
    double x15 = x14/x4;
    double x16 = x12*x9;
    double x17 = x16/((x4)*(x4));

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 0.00060000000000000006*x10*x11*x6 + 3.0000000000000004e-8*x10*x13*x15 - 0.0001*x11*x15*x9 + 1.8000000000000002e-7*x13/((x8)*(x8)*(x8)) + 2.5000000000000001e-9*((x14)*(x14))*x17 - 5.0000000000000001e-9*x14*x17 + 5.0000000000000001e-9*x15*x16*x3 - 6.0000000000000008e-8*x9);
}
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dP2)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/x8;
    double x10 = pow(x8, -2);
    double x11 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x12 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x13 = x12*x6;
    double x14 = x5 + 1;
    double x15 = x14/x4;
    double x16 = x12*x9;
    double x17 = x16/((x4)*(x4));

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 0.00060000000000000006*x10*x11*x6 + 3.0000000000000004e-8*x10*x13*x15 - 0.0001*x11*x15*x9 + 1.8000000000000002e-7*x13/((x8)*(x8)*(x8)) + 2.5000000000000001e-9*((x14)*(x14))*x17 - 5.0000000000000001e-9*x14*x17 + 5.0000000000000001e-9*x15*x16*x3 - 6.0000000000000008e-8*x9;
}
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
    
    double x0 = (*endmember[0].d3mu0dT3)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = x1 - 51.0;
    double x4 = 0.00030000000000000003*P;
    double x5 = sqrt(-x2/x3)/(x4 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = -n1*(-x0 + 0.000108*x5*(3.6000000000000007e-6*P - 1.2240000000000002)/((x2)*(x2)) - 6.4799999999999998e-7*x5*(-0.0051000000000000004*P + 0.61199999999999999*T + x3*(0.036000000000000004*T - x4 + 153.00000000000003) + 2601.0)/((x2)*(x2)*(x2)));
}
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT3)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = x1 - 51.0;
    double x4 = 0.00030000000000000003*P;
    double x5 = sqrt(-x2/x3)/(x4 - 153.0);

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 0.000108*x5*(3.6000000000000007e-6*P - 1.2240000000000002)/((x2)*(x2)) + 6.4799999999999998e-7*x5*(-0.0051000000000000004*P + 0.61199999999999999*T + x3*(0.036000000000000004*T - x4 + 153.00000000000003) + 2601.0)/((x2)*(x2)*(x2));
}
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
    
    double x0 = (*endmember[0].d3mu0dT2dP)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = 1.0/x2;
    double x4 = x1 - 51.0;
    double x5 = x2/x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/x8;
    double x10 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x11 = x6/((x8)*(x8));
    double x12 = pow(x2, -2);
    double x13 = x12*x9;
    double x14 = x10*x13;
    double x15 = -0.0051000000000000004*P + 0.61199999999999999*T + x4*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x16 = x15*x9/((x2)*(x2)*(x2));
    double x17 = x5 + 1;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 3.6000000000000003e-6*x10*x11*x3 + 1.0800000000000001e-8*x11*x12*x15 - 3.6000000000000001e-5*x13*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002) - 6.0000000000000008e-7*x14*x17 + 1.2000000000000002e-6*x14 + 1.8000000000000002e-9*x16*x17 - 7.2000000000000008e-9*x16 + 4.320000000000001e-8*x3*x9);
}
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT2dP)(T, P);
    double x1 = 0.0001*P;
    double x2 = 0.012*T - x1 + 51.0;
    double x3 = 1.0/x2;
    double x4 = x1 - 51.0;
    double x5 = x2/x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/x8;
    double x10 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x11 = x6/((x8)*(x8));
    double x12 = pow(x2, -2);
    double x13 = x12*x9;
    double x14 = x10*x13;
    double x15 = -0.0051000000000000004*P + 0.61199999999999999*T + x4*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x16 = x15*x9/((x2)*(x2)*(x2));
    double x17 = x5 + 1;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 3.6000000000000003e-6*x10*x11*x3 + 1.0800000000000001e-8*x11*x12*x15 - 3.6000000000000001e-5*x13*(-6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002) - 6.0000000000000008e-7*x14*x17 + 1.2000000000000002e-6*x14 + 1.8000000000000002e-9*x16*x17 - 7.2000000000000008e-9*x16 + 4.320000000000001e-8*x3*x9;
}
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
    
    double x0 = (*endmember[0].d3mu0dTdP2)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/((x8)*(x8));
    double x10 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x11 = x6/((x8)*(x8)*(x8));
    double x12 = 1.0/x4;
    double x13 = x6/x8;
    double x14 = x12*x13;
    double x15 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x16 = x12*x9;
    double x17 = pow(x4, -2);
    double x18 = x13*x17;
    double x19 = x15*x18;
    double x20 = x5 + 1;
    double x21 = x14*x20;
    double x22 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x23 = x17*x22*x9;
    double x24 = x13*x22;
    double x25 = x24/((x4)*(x4)*(x4));
    double x26 = x10*x20;
    double x27 = ((x20)*(x20));

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 1.8000000000000002e-7*x10*x11 + 2.5000000000000001e-9*x10*x18*x27 + 5.0000000000000001e-9*x10*x21*x3 + 1.0800000000000002e-9*x11*x12*x22 - 3.6000000000000005e-10*x14 - 3.6000000000000003e-6*x15*x16 + 3.0000000000000004e-8*x16*x26 + 2.9999999999999993e-11*x17*x20*x24*x3 - 5.0000000000000001e-9*x18*x26 - 6.0000000000000008e-7*x19*x20 + 1.2000000000000002e-6*x19 + 1.8000000000000002e-10*x20*x23 - 9.0000000000000012e-11*x20*x25 - 3.600000000000001e-10*x21 - 3.600000000000001e-10*x23 + 1.5e-11*x25*x27 + 1.2000000000000003e-10*x25 - 2.1600000000000008e-9*x9);
}
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dTdP2)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/((x8)*(x8));
    double x10 = 3.6000000000000007e-6*P - 1.2240000000000002;
    double x11 = x6/((x8)*(x8)*(x8));
    double x12 = 1.0/x4;
    double x13 = x6/x8;
    double x14 = x12*x13;
    double x15 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x16 = x12*x9;
    double x17 = pow(x4, -2);
    double x18 = x13*x17;
    double x19 = x15*x18;
    double x20 = x5 + 1;
    double x21 = x14*x20;
    double x22 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x23 = x17*x22*x9;
    double x24 = x13*x22;
    double x25 = x24/((x4)*(x4)*(x4));
    double x26 = x10*x20;
    double x27 = ((x20)*(x20));

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 1.8000000000000002e-7*x10*x11 + 2.5000000000000001e-9*x10*x18*x27 + 5.0000000000000001e-9*x10*x21*x3 + 1.0800000000000002e-9*x11*x12*x22 - 3.6000000000000005e-10*x14 - 3.6000000000000003e-6*x15*x16 + 3.0000000000000004e-8*x16*x26 + 2.9999999999999993e-11*x17*x20*x24*x3 - 5.0000000000000001e-9*x18*x26 - 6.0000000000000008e-7*x19*x20 + 1.2000000000000002e-6*x19 + 1.8000000000000002e-10*x20*x23 - 9.0000000000000012e-11*x20*x25 - 3.600000000000001e-10*x21 - 3.600000000000001e-10*x23 + 1.5e-11*x25*x27 + 1.2000000000000003e-10*x25 - 2.1600000000000008e-9*x9;
}
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
    
    double x0 = (*endmember[0].d3mu0dP3)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/((x8)*(x8));
    double x10 = pow(x8, -3);
    double x11 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x12 = x11*x6;
    double x13 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x14 = x13*x6;
    double x15 = 1.0/x8;
    double x16 = x5 + 1;
    double x17 = x16/x4;
    double x18 = x15*x17;
    double x19 = pow(x4, -2);
    double x20 = x12*x19;
    double x21 = x15*x16;
    double x22 = ((x16)*(x16));
    double x23 = x15*x22;
    double x24 = x16*x19;
    double x25 = x13*x9;
    double x26 = 4.5000000000000006e-12*x25;
    double x27 = x14/((x4)*(x4)*(x4));
    double x28 = x14/((x2)*(x2));
    double x29 = x14*x3;
    double x30 = x19*x23;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 5.4000000000000012e-7*x10*x12 - 2.7000000000000007e-11*x10*x14*x17 + 9.0000000000000012e-8*x11*x17*x9 + 1.5000000000000002e-8*x12*x18*x3 - 1.6200000000000004e-10*x14/((x8)*(x8)*(x8)*(x8)) - 1.2500000000000002e-13*x15*((x16)*(x16)*(x16))*x27 + 9.9999999999999998e-13*x15*x24*x29 - x17*x26*x3 - 1.0000000000000002e-12*x18*x28 + 9.0000000000000012e-12*x18*x6 - 2.2500000000000003e-12*x19*x22*x25 - 1.5000000000000002e-8*x20*x21 + 7.500000000000001e-9*x20*x23 - 9.9999999999999998e-13*x21*x27 + 7.5000000000000004e-13*x23*x27 + x24*x26 + 2.5000000000000001e-9*x28*x30*(2.0e-8*P - 0.010200000000000001) - 1.2499999999999999e-12*x29*x30 + 5.4000000000000007e-11*x9);
}
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dP3)(T, P);
    double x1 = 0.0001*P;
    double x2 = x1 - 51.0;
    double x3 = 1.0/x2;
    double x4 = 0.012*T - x1 + 51.0;
    double x5 = x3*x4;
    double x6 = sqrt(-x5);
    double x7 = 0.00030000000000000003*P;
    double x8 = x7 - 153.0;
    double x9 = x6/((x8)*(x8));
    double x10 = pow(x8, -3);
    double x11 = -6.0000000000000008e-8*P + 3.6000000000000007e-6*T + 0.025500000000000002;
    double x12 = x11*x6;
    double x13 = -0.0051000000000000004*P + 0.61199999999999999*T + x2*(0.036000000000000004*T - x7 + 153.00000000000003) + 2601.0;
    double x14 = x13*x6;
    double x15 = 1.0/x8;
    double x16 = x5 + 1;
    double x17 = x16/x4;
    double x18 = x15*x17;
    double x19 = pow(x4, -2);
    double x20 = x12*x19;
    double x21 = x15*x16;
    double x22 = ((x16)*(x16));
    double x23 = x15*x22;
    double x24 = x16*x19;
    double x25 = x13*x9;
    double x26 = 4.5000000000000006e-12*x25;
    double x27 = x14/((x4)*(x4)*(x4));
    double x28 = x14/((x2)*(x2));
    double x29 = x14*x3;
    double x30 = x19*x23;

if (T >= 0.0083333333333333332*P - 4250.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 5.4000000000000012e-7*x10*x12 - 2.7000000000000007e-11*x10*x14*x17 + 9.0000000000000012e-8*x11*x17*x9 + 1.5000000000000002e-8*x12*x18*x3 - 1.6200000000000004e-10*x14/((x8)*(x8)*(x8)*(x8)) - 1.2500000000000002e-13*x15*((x16)*(x16)*(x16))*x27 + 9.9999999999999998e-13*x15*x24*x29 - x17*x26*x3 - 1.0000000000000002e-12*x18*x28 + 9.0000000000000012e-12*x18*x6 - 2.2500000000000003e-12*x19*x22*x25 - 1.5000000000000002e-8*x20*x21 + 7.500000000000001e-9*x20*x23 - 9.9999999999999998e-13*x21*x27 + 7.5000000000000004e-13*x23*x27 + x24*x26 + 2.5000000000000001e-9*x28*x30*(2.0e-8*P - 0.010200000000000001) - 1.2499999999999999e-12*x29*x30 + 5.4000000000000007e-11*x9;
}
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

