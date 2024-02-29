#include <math.h>


static double coder_g(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    
    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = fmin(4, 0.41666666666666669*sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*(0.13589999999999997*P - 5.7599999999999998*T + x0 + 3252.48);
}
else {
   result = n1*(x0 + 1626.24*((x1)*(x1)*(x1)) - 1.0/3.0*(x1 - 1)*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002) - 1626.24);
}
    return result;
}
        
static void coder_dgdn(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = fmin(4, 0.41666666666666669*sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = 0.13589999999999997*P - 5.7599999999999998*T + x0 + 3252.48;
}
else {
   result[0] = x0 + 1626.24*((x1)*(x1)*(x1)) - 1.0/3.0*(x1 - 1)*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002) - 1626.24;
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
    
    double x0 = (*endmember[0].dmu0dT)(T, P) - 5.7599999999999998;
    double x1 = sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998);
    double x2 = 0.41666666666666669*x1;
    double x3 = fmin(4, x2);
    double x4 = (-x2 + 4 >= 0. ? 1. : 0.)/x1;

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 6.9120000000000008*((x3)*(x3))*x4 + 5.7599999999999998*x3 + 0.00047225501770956318*x4*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002));
}
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dT)(T, P) - 5.7599999999999998;
    double x1 = sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998);
    double x2 = 0.41666666666666669*x1;
    double x3 = fmin(4, x2);
    double x4 = (-x2 + 4 >= 0. ? 1. : 0.)/x1;

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 - 6.9120000000000008*((x3)*(x3))*x4 + 5.7599999999999998*x3 + 0.00047225501770956318*x4*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002);
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
    
    double x0 = (*endmember[0].dmu0dP)(T, P) + 0.13589999999999997;
    double x1 = sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998);
    double x2 = 0.41666666666666669*x1;
    double x3 = fmin(4, x2);
    double x4 = (-x2 + 4 >= 0. ? 1. : 0.)/x1;

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 0.16307999999999997*((x3)*(x3))*x4 - 0.13589999999999997*x3 - 1.1142266824085006e-5*x4*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002));
}
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dP)(T, P) + 0.13589999999999997;
    double x1 = sqrt(0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998);
    double x2 = 0.41666666666666669*x1;
    double x3 = fmin(4, x2);
    double x4 = (-x2 + 4 >= 0. ? 1. : 0.)/x1;

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.16307999999999997*((x3)*(x3))*x4 - 0.13589999999999997*x3 - 1.1142266824085006e-5*x4*(0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002);
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 0.023502432113341205*x10*x6 - 0.0097926800472255028*x10*x8 + 0.019585360094451006*((x4)*(x4))*x7*x9 + 1.6057785726133907e-6*x5*x6 + 6.6907440525557947e-7*x5*x8 - 0.016321133412042506*x4/x2);
}
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dT2)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 - 0.023502432113341205*x10*x6 - 0.0097926800472255028*x10*x8 + 0.019585360094451006*((x4)*(x4))*x7*x9 + 1.6057785726133907e-6*x5*x6 + 6.6907440525557947e-7*x5*x8 - 0.016321133412042506*x4/x2;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 0.00055451050767414407*x10*x6 + 0.0002310460448642267*x10*x8 - 0.00046209208972845339*((x4)*(x4))*x7*x9 - 3.788633819759719e-8*x5*x6 - 1.5785974248998829e-8*x5*x8 + 0.00038507674144037777*x4/x2);
}
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dTdP)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.00055451050767414407*x10*x6 + 0.0002310460448642267*x10*x8 - 0.00046209208972845339*((x4)*(x4))*x7*x9 - 3.788633819759719e-8*x5*x6 - 1.5785974248998829e-8*x5*x8 + 0.00038507674144037777*x4/x2;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 1.3082982290436834e-5*x10*x6 - 5.4512426210153473e-6*x10*x8 + 1.0902485242030695e-5*((x4)*(x4))*x7*x9 + 8.9388079184955861e-10*x5*x6 + 3.7245032993731606e-10*x5*x8 - 9.0854043683589116e-6*x4/x2);
}
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dP2)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = sqrt(x1);
    double x3 = 0.41666666666666669*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/x1;
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 - 1.3082982290436834e-5*x10*x6 - 5.4512426210153473e-6*x10*x8 + 1.0902485242030695e-5*((x4)*(x4))*x7*x9 + 8.9388079184955861e-10*x5*x6 + 3.7245032993731606e-10*x5*x8 - 9.0854043683589116e-6*x4/x2;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 8.3243561204278187e-5*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 0.00019978454689026763*x10*x13*((x4)*(x4)) - 9.9892273445133817e-5*x11*x14 + 6.8250328942245659e-9*x11*x8 + 1.3873926867379697e-5*x12*x14 - 9.4792123530896746e-10*x12*x8 + x13*x5*x7 - 0.00023974145626832113*x14*x9 - 2.7747853734759395e-5*x2*((x4)*(x4)*(x4)) - x5 + 1.6380078946138955e-8*x8*x9 - 3.4684817168449249e-5*x7/x1);
}
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT3)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 8.3243561204278187e-5*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.00019978454689026763*x10*x13*((x4)*(x4)) - 9.9892273445133817e-5*x11*x14 + 6.8250328942245659e-9*x11*x8 + 1.3873926867379697e-5*x12*x14 - 9.4792123530896746e-10*x12*x8 + x13*x5*x7 - 0.00023974145626832113*x14*x9 - 2.7747853734759395e-5*x2*((x4)*(x4)*(x4)) - x5 + 1.6380078946138955e-8*x8*x9 - 3.4684817168449249e-5*x7/x1;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 1.9640277721634381e-6*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 4.7136666531922514e-6*x10*x13*((x4)*(x4)) + 2.3568333265961257e-6*x11*x14 - 1.6102811984811083e-10*x11*x8 - 3.2733796202723972e-7*x12*x14 + 2.2365016645570949e-11*x12*x8 - x13*x5*x7 + 5.6563999838307015e-6*x14*x9 + 6.5467592405447943e-7*x2*((x4)*(x4)*(x4)) + x5 - 3.8646748763546596e-10*x8*x9 + 8.1834490506809934e-7*x7/x1);
}
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT2dP)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 1.9640277721634381e-6*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 - 4.7136666531922514e-6*x10*x13*((x4)*(x4)) + 2.3568333265961257e-6*x11*x14 - 1.6102811984811083e-10*x11*x8 - 3.2733796202723972e-7*x12*x14 + 2.2365016645570949e-11*x12*x8 - x13*x5*x7 + 5.6563999838307015e-6*x14*x9 + 6.5467592405447943e-7*x2*((x4)*(x4)*(x4)) + x5 - 3.8646748763546596e-10*x8*x9 + 8.1834490506809934e-7*x7/x1;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 4.6338780249481119e-8*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 1.1121307259875467e-7*x10*x13*((x4)*(x4)) - 5.5606536299377337e-8*x11*x14 + 3.7992572026663651e-12*x11*x8 + 7.7231300415801864e-9*x12*x14 - 5.2767461148143961e-13*x12*x8 + x13*x5*x7 - 1.3345568711850563e-7*x14*x9 - 1.5446260083160373e-8*x2*((x4)*(x4)*(x4)) - x5 + 9.1182172863992766e-12*x8*x9 - 1.9307825103950464e-8*x7/x1);
}
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dTdP2)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 4.6338780249481119e-8*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 + 1.1121307259875467e-7*x10*x13*((x4)*(x4)) - 5.5606536299377337e-8*x11*x14 + 3.7992572026663651e-12*x11*x8 + 7.7231300415801864e-9*x12*x14 - 5.2767461148143961e-13*x12*x8 + x13*x5*x7 - 1.3345568711850563e-7*x14*x9 - 1.5446260083160373e-8*x2*((x4)*(x4)*(x4)) - x5 + 9.1182172863992766e-12*x8*x9 - 1.9307825103950464e-8*x7/x1;
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
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 1.093305596511195e-9*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 2.6239334316268681e-9*x10*x13*((x4)*(x4)) + 1.311966715813434e-9*x11*x14 - 8.9638724625409544e-14*x11*x8 - 1.8221759941853249e-10*x12*x14 + 1.2449822864640213e-14*x12*x8 - x13*x5*x7 + 3.1487201179522416e-9*x14*x9 + 3.6443519883706497e-10*x2*((x4)*(x4)*(x4)) + x5 - 2.151329391009829e-13*x8*x9 + 4.5554399854633119e-10*x7/x1);
}
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dP3)(T, P);
    double x1 = 0.00016044864226682408*P - 0.0068004722550177093*T + 5.7599999999999998;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.41666666666666669*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 1.093305596511195e-9*x2*x4;
    double x6 = x3 - 4;
    double x7 = 0;
    double x8 = 0.40769999999999995*P - 17.280000000000001*T + 14636.160000000002;
    double x9 = x4/pow(x1, 5.0/2.0);
    double x10 = pow(x1, -2);
    double x11 = x10*x7;
    double x12 = x2*0;
    double x13 = fmin(4, x3);
    double x14 = ((x13)*(x13));

if (T >= 0.023593749999999997*P + 846.99999999999989) {
   result[0] = x0;
}
else {
   result[0] = x0 - 2.6239334316268681e-9*x10*x13*((x4)*(x4)) + 1.311966715813434e-9*x11*x14 - 8.9638724625409544e-14*x11*x8 - 1.8221759941853249e-10*x12*x14 + 1.2449822864640213e-14*x12*x8 - x13*x5*x7 + 3.1487201179522416e-9*x14*x9 + 3.6443519883706497e-10*x2*((x4)*(x4)*(x4)) + x5 - 2.151329391009829e-13*x8*x9 + 4.5554399854633119e-10*x7/x1;
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

