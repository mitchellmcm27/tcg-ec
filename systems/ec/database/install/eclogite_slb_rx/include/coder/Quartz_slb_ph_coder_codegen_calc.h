#include <math.h>


static double coder_g(double T, double P, double n[1]) {
    double n1 = n[0];
    double result;
    
    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 1.0*sqrt((-5.1639999999999997*T + x1)/x1)*(0.044798520000000008*((P)*(P)) - 1.8931224*P*T + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988)/(0.36660000000000004*P + 13121.723999999998));
}
    return result;
}
        
static void coder_dgdn(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 1.0*sqrt((-5.1639999999999997*T + x1)/x1)*(0.044798520000000008*((P)*(P)) - 1.8931224*P*T + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988)/(0.36660000000000004*P + 13121.723999999998);
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
    double x1 = 1.8931224*P;
    double x2 = 0.1222*P + 4373.9079999999994;
    double x3 = -5.1639999999999997*T + x2;
    double x4 = sqrt(x3/x2)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 1.0*x4*(-x1 - 45173.721823999986) + 2.5819999999999999*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x1 - 45173.721823999986*T + 38262142.384927988)/x3);
}
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = 1.8931224*P;
    double x2 = 0.1222*P + 4373.9079999999994;
    double x3 = -5.1639999999999997*T + x2;
    double x4 = sqrt(x3/x2)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 1.0*x4*(-x1 - 45173.721823999986) + 2.5819999999999999*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x1 - 45173.721823999986*T + 38262142.384927988)/x3;
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = sqrt(x2*x3);
    double x5 = 1.8931224*T;
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = 1.0/x6;
    double x8 = x4*(0.044798520000000008*((P)*(P)) - P*x5 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988);

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - x1*x7*x8*(0.061100000000000002*x2 - 0.061100000000000002*x3/((x1)*(x1)))/x3 - x4*x7*(0.089597040000000017*P - x5 + 2672.4577880000002) + 0.36660000000000004*x8/((x6)*(x6)));
}
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = sqrt(x2*x3);
    double x5 = 1.8931224*T;
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = 1.0/x6;
    double x8 = x4*(0.044798520000000008*((P)*(P)) - P*x5 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988);

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - x1*x7*x8*(0.061100000000000002*x2 - 0.061100000000000002*x3/((x1)*(x1)))/x3 - x4*x7*(0.089597040000000017*P - x5 + 2672.4577880000002) + 0.36660000000000004*x8/((x6)*(x6));
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.8931224*P;
    double x4 = sqrt(x2/x1)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 5.1639999999999997*x4*(x3 + 45173.721823999986)/x2 + 6.6667239999999994*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x3 - 45173.721823999986*T + 38262142.384927988)/((x2)*(x2)));
}
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dT2)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.8931224*P;
    double x4 = sqrt(x2/x1)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 5.1639999999999997*x4*(x3 + 45173.721823999986)/x2 + 6.6667239999999994*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x3 - 45173.721823999986*T + 38262142.384927988)/((x2)*(x2));
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = x2/x1;
    double x4 = sqrt(x3);
    double x5 = 0.36660000000000004*P + 13121.723999999998;
    double x6 = x4/x5;
    double x7 = 1.8931224*P;
    double x8 = x7 + 45173.721823999986;
    double x9 = x4/((x5)*(x5));
    double x10 = 1.0/x2;
    double x11 = x10*x6;
    double x12 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x7 - 45173.721823999986*T + 38262142.384927988;
    double x13 = x12*x6/((x2)*(x2));
    double x14 = -x3 + 1;

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 0.94656119999999999*x10*x12*x9 + 0.061100000000000002*x11*x14*x8 + 2.5819999999999999*x11*(0.089597040000000017*P - 1.8931224*T + 2672.4577880000002) + 0.15776019999999999*x13*x14 - 0.31552039999999998*x13 + 1.8931224*x6 - 0.36660000000000004*x8*x9);
}
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dTdP)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = x2/x1;
    double x4 = sqrt(x3);
    double x5 = 0.36660000000000004*P + 13121.723999999998;
    double x6 = x4/x5;
    double x7 = 1.8931224*P;
    double x8 = x7 + 45173.721823999986;
    double x9 = x4/((x5)*(x5));
    double x10 = 1.0/x2;
    double x11 = x10*x6;
    double x12 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x7 - 45173.721823999986*T + 38262142.384927988;
    double x13 = x12*x6/((x2)*(x2));
    double x14 = -x3 + 1;

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 0.94656119999999999*x10*x12*x9 + 0.061100000000000002*x11*x14*x8 + 2.5819999999999999*x11*(0.089597040000000017*P - 1.8931224*T + 2672.4577880000002) + 0.15776019999999999*x13*x14 - 0.31552039999999998*x13 + 1.8931224*x6 - 0.36660000000000004*x8*x9;
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/x6;
    double x8 = pow(x6, -2);
    double x9 = 1.8931224*T;
    double x10 = 0.089597040000000017*P - x9 + 2672.4577880000002;
    double x11 = 0.044798520000000008*((P)*(P)) - P*x9 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988;
    double x12 = x11*x5;
    double x13 = -x4 + 1;
    double x14 = x13/x3;
    double x15 = x14*x7;
    double x16 = x11*x7/((x3)*(x3));

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = -n1*(-x0 + 0.1222*x10*x15 - 0.73320000000000007*x10*x5*x8 - 0.0074664200000000005*x11*x15*x2 - 0.044798520000000008*x12*x14*x8 + 0.26879112000000005*x12/((x6)*(x6)*(x6)) + 0.0037332100000000003*((x13)*(x13))*x16 - 0.0074664200000000005*x13*x16 + 0.089597040000000017*x7);
}
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d2mu0dP2)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/x6;
    double x8 = 1.8931224*T;
    double x9 = 0.089597040000000017*P - x8 + 2672.4577880000002;
    double x10 = x5/((x6)*(x6));
    double x11 = 0.044798520000000008*((P)*(P)) - P*x8 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988;
    double x12 = 1.0/x3;
    double x13 = -x4 + 1;
    double x14 = x11*x13;
    double x15 = x12*x14;
    double x16 = x7/((x3)*(x3));

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.044798520000000008*x10*x15 + 0.73320000000000007*x10*x9 - 0.0037332100000000003*x11*((x13)*(x13))*x16 - 0.26879112000000005*x11*x5/((x6)*(x6)*(x6)) - 0.1222*x12*x13*x7*x9 + 0.0074664200000000005*x14*x16 + 0.0074664200000000005*x15*x2*x7 - 0.089597040000000017*x7;
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.8931224*P;
    double x4 = sqrt(x2/x1)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 20.000171999999999*x4*(x3 + 45173.721823999986)/((x2)*(x2)) + 51.64044410399999*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x3 - 45173.721823999986*T + 38262142.384927988)/((x2)*(x2)*(x2)));
}
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT3)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.8931224*P;
    double x4 = sqrt(x2/x1)/(0.36660000000000004*P + 13121.723999999998);

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 20.000171999999999*x4*(x3 + 45173.721823999986)/((x2)*(x2)) + 51.64044410399999*x4*(0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x3 - 45173.721823999986*T + 38262142.384927988)/((x2)*(x2)*(x2));
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.0/x2;
    double x4 = x2/x1;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/x6;
    double x8 = 1.8931224*P;
    double x9 = x8 + 45173.721823999986;
    double x10 = x5/((x6)*(x6));
    double x11 = pow(x2, -2);
    double x12 = x11*x7;
    double x13 = x12*x9;
    double x14 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x8 - 45173.721823999986*T + 38262142.384927988;
    double x15 = x14*x7/((x2)*(x2)*(x2));
    double x16 = -x4 + 1;

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 - 2.4440210184*x10*x11*x14 + 1.8931224*x10*x3*x9 + 6.6667239999999994*x12*(0.089597040000000017*P - 1.8931224*T + 2672.4577880000002) - 0.31552039999999998*x13*x16 + 0.63104079999999996*x13 + 0.40733683639999996*x15*x16 - 1.6293473455999998*x15 - 9.7760840735999999*x3*x7);
}
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dT2dP)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = -5.1639999999999997*T + x1;
    double x3 = 1.0/x2;
    double x4 = x2/x1;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/x6;
    double x8 = 1.8931224*P;
    double x9 = x8 + 45173.721823999986;
    double x10 = x5/((x6)*(x6));
    double x11 = pow(x2, -2);
    double x12 = x11*x7;
    double x13 = x12*x9;
    double x14 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x8 - 45173.721823999986*T + 38262142.384927988;
    double x15 = x14*x7/((x2)*(x2)*(x2));
    double x16 = -x4 + 1;

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 - 2.4440210184*x10*x11*x14 + 1.8931224*x10*x3*x9 + 6.6667239999999994*x12*(0.089597040000000017*P - 1.8931224*T + 2672.4577880000002) - 0.31552039999999998*x13*x16 + 0.63104079999999996*x13 + 0.40733683639999996*x15*x16 - 1.6293473455999998*x15 - 9.7760840735999999*x3*x7;
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/((x6)*(x6));
    double x8 = 1.8931224*P;
    double x9 = x8 + 45173.721823999986;
    double x10 = x5/((x6)*(x6)*(x6));
    double x11 = 1.0/x3;
    double x12 = x5/x6;
    double x13 = x11*x12;
    double x14 = 0.089597040000000017*P - 1.8931224*T + 2672.4577880000002;
    double x15 = x11*x7;
    double x16 = pow(x3, -2);
    double x17 = x12*x16;
    double x18 = x14*x17;
    double x19 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x8 - 45173.721823999986*T + 38262142.384927988;
    double x20 = x16*x19;
    double x21 = x20*x7;
    double x22 = x12*x19/((x3)*(x3)*(x3));
    double x23 = -x4 + 1;
    double x24 = x13*x23;
    double x25 = x17*x9;
    double x26 = ((x23)*(x23));

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = n1*(x0 + 0.6940186718400001*x10*x11*x19 + 0.26879112000000005*x10*x9 - 0.019278296439999999*x12*x2*x20*x23 + 0.23133955728000002*x13 - 1.8931224*x14*x15 - 0.044798520000000008*x15*x23*x9 + 0.31552039999999998*x18*x23 - 0.63104079999999996*x18 - 0.0074664200000000005*x2*x24*x9 - 0.11566977864*x21*x23 + 0.23133955728*x21 - 0.057834889319999992*x22*x23 + 0.0096391482199999993*x22*x26 + 0.077113185759999994*x22 - 0.0074664200000000005*x23*x25 + 0.23133955728*x24 + 0.0037332100000000003*x25*x26 - 1.3880373436800002*x7);
}
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dTdP2)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/((x6)*(x6));
    double x8 = 1.8931224*P;
    double x9 = x8 + 45173.721823999986;
    double x10 = x5/((x6)*(x6)*(x6));
    double x11 = 1.0/x3;
    double x12 = x5/x6;
    double x13 = x11*x12;
    double x14 = 0.089597040000000017*P - 1.8931224*T + 2672.4577880000002;
    double x15 = x11*x7;
    double x16 = pow(x3, -2);
    double x17 = x12*x16;
    double x18 = x14*x17;
    double x19 = 0.044798520000000008*((P)*(P)) + 2672.4577880000002*P - T*x8 - 45173.721823999986*T + 38262142.384927988;
    double x20 = x16*x19;
    double x21 = x20*x7;
    double x22 = x12*x19/((x3)*(x3)*(x3));
    double x23 = -x4 + 1;
    double x24 = x13*x23;
    double x25 = x17*x9;
    double x26 = ((x23)*(x23));

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.6940186718400001*x10*x11*x19 + 0.26879112000000005*x10*x9 - 0.019278296439999999*x12*x2*x20*x23 + 0.23133955728000002*x13 - 1.8931224*x14*x15 - 0.044798520000000008*x15*x23*x9 + 0.31552039999999998*x18*x23 - 0.63104079999999996*x18 - 0.0074664200000000005*x2*x24*x9 - 0.11566977864*x21*x23 + 0.23133955728*x21 - 0.057834889319999992*x22*x23 + 0.0096391482199999993*x22*x26 + 0.077113185759999994*x22 - 0.0074664200000000005*x23*x25 + 0.23133955728*x24 + 0.0037332100000000003*x25*x26 - 1.3880373436800002*x7;
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
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/((x6)*(x6));
    double x8 = pow(x6, -3);
    double x9 = 1.8931224*T;
    double x10 = 0.089597040000000017*P - x9 + 2672.4577880000002;
    double x11 = x10*x5;
    double x12 = 0.044798520000000008*((P)*(P)) - P*x9 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988;
    double x13 = x12*x5;
    double x14 = 1.0/x3;
    double x15 = 1.0/x6;
    double x16 = -x4 + 1;
    double x17 = x15*x16;
    double x18 = x14*x17;
    double x19 = x16*x7;
    double x20 = pow(x3, -2);
    double x21 = x11*x20;
    double x22 = ((x16)*(x16));
    double x23 = x15*x22;
    double x24 = x12*x20;
    double x25 = 0.0082115687160000013*x19;
    double x26 = x13/((x3)*(x3)*(x3));
    double x27 = 0.0018247930480000001*x17;
    double x28 = x15*x26;
    double x29 = x13/((x1)*(x1));
    double x30 = x13*x2*x20;

if (T >= 0.02366382649109218*P + 847.0) {
   result = n1*x0;
}
else {
   result = -n1*(-x0 - 0.13439556000000003*x10*x14*x19 - 0.022399259999999997*x11*x18*x2 + 0.80637336000000015*x11*x8 + x12*x14*x2*x25 + 0.049269412296000008*x13*x14*x16*x8 - 0.29561647377600009*x13/((x6)*(x6)*(x6)*(x6)) + 0.00022809913100000002*((x16)*(x16)*(x16))*x28 - 0.022399260000000004*x17*x21 + 0.0018247930480000006*x18*x29 + 0.016423137432000003*x18*x5 + 0.0037332100000000003*x20*x23*x29*(0.029865680000000002*P + 1068.9831151999999) + 0.01119963*x21*x23 - 0.0041057843580000006*x22*x24*x7 - 0.0013685947860000002*x22*x28 - 0.0022809913100000001*x23*x30 + x24*x25 + x26*x27 + x27*x30 - 0.098538824592000029*x7);
}
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[1], double result[1]) {
    double n1 = n[0];

    double x0 = (*endmember[0].d3mu0dP3)(T, P);
    double x1 = 0.1222*P + 4373.9079999999994;
    double x2 = 1.0/x1;
    double x3 = -5.1639999999999997*T + x1;
    double x4 = x2*x3;
    double x5 = sqrt(x4);
    double x6 = 0.36660000000000004*P + 13121.723999999998;
    double x7 = x5/((x6)*(x6));
    double x8 = pow(x6, -3);
    double x9 = 1.8931224*T;
    double x10 = 0.089597040000000017*P - x9 + 2672.4577880000002;
    double x11 = x10*x5;
    double x12 = 0.044798520000000008*((P)*(P)) - P*x9 + 2672.4577880000002*P - 45173.721823999986*T + 38262142.384927988;
    double x13 = x12*x5;
    double x14 = 1.0/x3;
    double x15 = 1.0/x6;
    double x16 = -x4 + 1;
    double x17 = x15*x16;
    double x18 = x14*x17;
    double x19 = x14*x16;
    double x20 = pow(x3, -2);
    double x21 = x11*x20;
    double x22 = ((x16)*(x16));
    double x23 = x15*x22;
    double x24 = x12*x7;
    double x25 = x20*x24;
    double x26 = x13/((x3)*(x3)*(x3));
    double x27 = 0.0018247930480000001*x17;
    double x28 = x13/((x1)*(x1));
    double x29 = x13*x2*x20;

if (T >= 0.02366382649109218*P + 847.0) {
   result[0] = x0;
}
else {
   result[0] = x0 + 0.13439556000000003*x10*x19*x7 + 0.022399259999999997*x11*x18*x2 - 0.80637336000000015*x11*x8 - 0.049269412296000008*x13*x19*x8 + 0.29561647377600009*x13/((x6)*(x6)*(x6)*(x6)) - 0.00022809913100000002*x15*((x16)*(x16)*(x16))*x26 - 0.0082115687160000013*x16*x25 + 0.022399260000000004*x17*x21 - 0.0018247930480000006*x18*x28 - 0.016423137432000003*x18*x5 - 0.0082115687160000013*x19*x2*x24 - 0.0037332100000000003*x20*x23*x28*(0.029865680000000002*P + 1068.9831151999999) - 0.01119963*x21*x23 + 0.0041057843580000006*x22*x25 + 0.0013685947860000002*x23*x26 + 0.0022809913100000001*x23*x29 - x26*x27 - x27*x29 + 0.098538824592000029*x7;
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

