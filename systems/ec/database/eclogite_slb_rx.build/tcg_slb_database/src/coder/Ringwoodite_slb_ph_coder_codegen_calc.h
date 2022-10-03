#include <math.h>


static double coder_g(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = n1 + n2;
    double x1 = 1.0/x0;

result = 1.0*x1*(9100.0*n1*n2 + 1.0*x0*(16.628925236306479*T*(n1*log(n1*x1) + n2*log(n2*x1)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P)));
    return result;
}
        
static void coder_dgdn(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -2);
    double x2 = 9100.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = n1*x5;
    double x7 = log(x6);
    double x8 = n2*x5;
    double x9 = log(x8);
    double x10 = 16.628925236306479*T;
    double x11 = x10*(n1*x7 + n2*x9);
    double x12 = -1.0*x1*(n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x11));
    double x13 = 1.0*n1;
    double x14 = 1.0*n2;
    double x15 = x13 + x14;
    double x16 = x11 + x13*x3 + x14*x4;
    double x17 = 1.0*x5;

result[0] = x12 + x17*(x15*(x10*(x0*(-n1*x1 + x5) + x7 - x8) + x3) + x16 + x2);
result[1] = x12 + x17*(9100.0*n1 + x15*(x10*(x0*(-n2*x1 + x5) - x6 + x9) + x4) + x16);
}
        
static void coder_d2gdn2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -3);
    double x2 = 9100.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = n1*x5;
    double x7 = log(x6);
    double x8 = n2*x5;
    double x9 = log(x8);
    double x10 = 16.628925236306479*T;
    double x11 = x10*(n1*x7 + n2*x9);
    double x12 = 2.0*x1*(n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x11));
    double x13 = pow(x0, -2);
    double x14 = 1.0*n1;
    double x15 = 1.0*n2;
    double x16 = x14 + x15;
    double x17 = n1*x13;
    double x18 = -x17;
    double x19 = x18 + x5;
    double x20 = x0*x19;
    double x21 = T*(x20 + x7 - x8);
    double x22 = 16.628925236306479*x21;
    double x23 = x11 + x14*x3 + x15*x4;
    double x24 = x13*(x16*(x22 + x3) + x2 + x23);
    double x25 = n2*x13;
    double x26 = -2*x13;
    double x27 = 2*x1;
    double x28 = n1*x27;
    double x29 = x10*x16;
    double x30 = 1.0*x5;
    double x31 = -x25 + x5;
    double x32 = x0*x31;
    double x33 = x32 - x6 + x9;
    double x34 = x10*x33;
    double x35 = x13*(9100.0*n1 + x16*(x34 + x4) + x23);

result[0] = x12 - 2.0*x24 + x30*(33.257850472612958*x21 + x29*(x0*(x26 + x28) + x19 + x25 + x20/n1) + 2.0*x3);
result[1] = x12 - 1.0*x24 + x30*(x22 + x29*(x0*(-x13 + x28) + x18 + x25 - x5) + 1.0*x3 + x34 + 1.0*x4 + 9100.0) - 1.0*x35;
result[2] = x12 + x30*(33.257850472612958*T*x33 + x29*(x0*(n2*x27 + x26) + x17 + x31 + x32/n2) + 2.0*x4) - 2.0*x35;
}
        
static void coder_d3gdn3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -4);
    double x2 = 9100.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = n1*x5;
    double x7 = log(x6);
    double x8 = n2*x5;
    double x9 = log(x8);
    double x10 = 16.628925236306479*T;
    double x11 = x10*(n1*x7 + n2*x9);
    double x12 = -6.0*x1*(n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x11));
    double x13 = pow(x0, -3);
    double x14 = 1.0*n1;
    double x15 = 1.0*n2;
    double x16 = x14 + x15;
    double x17 = pow(x0, -2);
    double x18 = n1*x17;
    double x19 = -x18;
    double x20 = x19 + x5;
    double x21 = x0*x20;
    double x22 = x21 + x7 - x8;
    double x23 = x10*x22;
    double x24 = x11 + x14*x3 + x15*x4;
    double x25 = x13*(x16*(x23 + x3) + x2 + x24);
    double x26 = 33.257850472612958*T;
    double x27 = n2*x17;
    double x28 = -2*x17;
    double x29 = 2*x13;
    double x30 = n1*x29;
    double x31 = x0*(x28 + x30);
    double x32 = 1.0/n1;
    double x33 = x20*x32;
    double x34 = T*(x0*x33 + x20 + x27 + x31);
    double x35 = 16.628925236306479*x34;
    double x36 = x17*(x16*x35 + x22*x26 + 2.0*x3);
    double x37 = -4*x17;
    double x38 = 6*x13;
    double x39 = 6*x1;
    double x40 = -n1*x39;
    double x41 = n2*x29;
    double x42 = 4*x13;
    double x43 = n1*x42 - x41;
    double x44 = x33 + x43;
    double x45 = x10*x16;
    double x46 = 1.0*x5;
    double x47 = x0*(-x17 + x30);
    double x48 = x19 + x27 + x47 - x5;
    double x49 = x26*x48;
    double x50 = -x27 + x5;
    double x51 = x0*x50;
    double x52 = x51 - x6 + x9;
    double x53 = x10*x52;
    double x54 = x13*(9100.0*n1 + x16*(x4 + x53) + x24);
    double x55 = x12 - 2.0*x17*(x23 + 1.0*x3 + 1.0*x4 + x45*x48 + x53 + 9100.0);
    double x56 = x0*(x28 + x41);
    double x57 = 1.0/n2;
    double x58 = x18 + x50 + x51*x57 + x56;
    double x59 = x10*x58;
    double x60 = x17*(x16*x59 + x26*x52 + 2.0*x4);

result[0] = x12 + 6.0*x25 - 3.0*x36 + x46*(49.886775708919437*x34 + x45*(x0*(x38 + x40) + x31*x32 + x37 + x44 - x21/((n1)*(n1))));
result[1] = 4.0*x25 - 1.0*x36 + x46*(x35 + x45*(x0*(x40 + x42) + x28 + x32*x47 + x44) + x49) + 2.0*x54 + x55;
result[2] = 2.0*x25 + x46*(x45*(x0*(x29 + x40) + x17 + x43) + x49 + x59) + 4.0*x54 + x55 - 1.0*x60;
result[3] = x12 + x46*(49.886775708919437*T*x58 + x45*(n2*x42 + x0*(-n2*x39 + x38) - x30 + x37 + x50*x57 + x56*x57 - x51/((n2)*(n2)))) + 6.0*x54 - 3.0*x60;
}
        
static double coder_dgdt(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = 1.0/(n1 + n2);

result = 1.0*x0*(1.0*n1 + 1.0*n2)*(16.628925236306479*n1*log(n1*x0) + n1*(*endmember[0].dmu0dT)(T, P) + 16.628925236306479*n2*log(n2*x0) + n2*(*endmember[1].dmu0dT)(T, P));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = n1 + n2;
    double x2 = 1.0/x1;
    double x3 = n2*x2;
    double x4 = n1*x2;
    double x5 = 16.628925236306479*log(x4);
    double x6 = pow(x1, -2);
    double x7 = 16.628925236306479*x1;
    double x8 = 1.0*n1 + 1.0*n2;
    double x9 = 1.0*x2;
    double x10 = x8*x9;
    double x11 = (*endmember[1].dmu0dT)(T, P);
    double x12 = 16.628925236306479*log(x3);
    double x13 = n1*x0 + n1*x5 + n2*x11 + n2*x12;
    double x14 = -1.0*x13*x6*x8 + x13*x9;

result[0] = x10*(x0 - 16.628925236306479*x3 + x5 + x7*(-n1*x6 + x2)) + x14;
result[1] = x10*(x11 + x12 - 16.628925236306479*x4 + x7*(-n2*x6 + x2)) + x14;
}
        
static void coder_d3gdn2dt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = (*endmember[0].dmu0dT)(T, P);
    double x3 = n2*x1;
    double x4 = n1*x1;
    double x5 = 16.628925236306479*log(x4);
    double x6 = pow(x0, -2);
    double x7 = n1*x6;
    double x8 = 16.628925236306479*x0;
    double x9 = x8*(x1 - x7);
    double x10 = x2 - 16.628925236306479*x3 + x5 + x9;
    double x11 = x1*x10;
    double x12 = 16.628925236306479*x1;
    double x13 = 16.628925236306479*n2;
    double x14 = 16.628925236306479*n1 + x13;
    double x15 = -2*x6;
    double x16 = pow(x0, -3);
    double x17 = 2*x16;
    double x18 = n1*x17;
    double x19 = x13*x6;
    double x20 = 16.628925236306479*x7;
    double x21 = x19 - x20;
    double x22 = 1.0*n1 + 1.0*n2;
    double x23 = 1.0*x22;
    double x24 = x1*x23;
    double x25 = x10*x6;
    double x26 = 2.0*x22;
    double x27 = (*endmember[1].dmu0dT)(T, P);
    double x28 = log(x3);
    double x29 = 2.0*n1*x2 + 2.0*n1*x5 + 2.0*n2*x27 + 2.0*x13*x28;
    double x30 = x16*x22*x29 - x29*x6;
    double x31 = x8*(-n2*x6 + x1);
    double x32 = x27 + 16.628925236306479*x28 + x31 - 16.628925236306479*x4;
    double x33 = x1*x32;
    double x34 = x32*x6;

result[0] = 2.0*x11 + x24*(x12 + x14*(x15 + x18) + x21 + x9/n1) - x25*x26 + x30;
result[1] = 1.0*x11 - x23*x25 - x23*x34 + x24*(-x12 + x14*(x18 - x6) + x21) + x30 + 1.0*x33;
result[2] = x24*(x12 + x14*(n2*x17 + x15) - x19 + x20 + x31/n2) - x26*x34 + x30 + 2.0*x33;
}
        
static void coder_d4gdn3dt(double T, double P, double n[2], double result[2]) {
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
    double x8 = pow(x0, -3);
    double x9 = 2*x8;
    double x10 = n1*x9;
    double x11 = x10 + x7;
    double x12 = 1.0/n1;
    double x13 = -16.628925236306479*n1*x6 + 16.628925236306479*x1;
    double x14 = x0*x13;
    double x15 = x4*x6;
    double x16 = x3*x6;
    double x17 = x15 - x16;
    double x18 = x11*x5 + x12*x14 + x17 + x2;
    double x19 = x1*x18;
    double x20 = (*endmember[0].dmu0dT)(T, P);
    double x21 = log(n1*x1);
    double x22 = -n2*x2 + x14 + x20 + 16.628925236306479*x21;
    double x23 = x22*x6;
    double x24 = -66.515700945225916*x6;
    double x25 = 6*x8;
    double x26 = pow(x0, -4);
    double x27 = 6*x26;
    double x28 = -n1*x27;
    double x29 = 16.628925236306479*x0;
    double x30 = x12*x29;
    double x31 = 66.515700945225916*x8;
    double x32 = 33.257850472612958*x8;
    double x33 = n1*x31 - n2*x32;
    double x34 = x12*x13 + x33;
    double x35 = 1.0*n1 + 1.0*n2;
    double x36 = 1.0*x35;
    double x37 = x1*x36;
    double x38 = 6.0*x8;
    double x39 = x22*x35;
    double x40 = x18*x6;
    double x41 = 3.0*x35;
    double x42 = (*endmember[1].dmu0dT)(T, P);
    double x43 = log(n2*x1);
    double x44 = n1*x20 + n2*x42 + x21*x3 + x4*x43;
    double x45 = -6.0*x26*x35*x44 + x38*x44;
    double x46 = -16.628925236306479*n2*x6 + 16.628925236306479*x1;
    double x47 = x0*x46;
    double x48 = -n1*x2 + x42 + 16.628925236306479*x43 + x47;
    double x49 = x48*x6;
    double x50 = x10 - x6;
    double x51 = 2.0*x8;
    double x52 = x35*x48;
    double x53 = 4.0*x8;
    double x54 = 2.0*x17 - 2.0*x2 + 2.0*x5*x50;
    double x55 = x1*x54 - x35*x54*x6 + x45;
    double x56 = 1.0/n2;
    double x57 = n2*x9 + x7;
    double x58 = -x15 + x16 + x2 + x47*x56 + x5*x57;
    double x59 = x1*x58;
    double x60 = x58*x6;

result[0] = 3.0*x19 - 6.0*x23 + x37*(x11*x30 + x24 + x34 + x5*(x25 + x28) - x14/((n1)*(n1))) + x38*x39 - x40*x41 + x45;
result[1] = 1.0*x19 - 4.0*x23 - x36*x40 + x37*(x30*x50 + x34 + x5*(x28 + 4*x8) - 33.257850472612958*x6) + x39*x53 - 2.0*x49 + x51*x52 + x55;
result[2] = -2.0*x23 - x36*x60 + x37*(x33 + x5*(x28 + x9) + 16.628925236306479*x6) + x39*x51 - 4.0*x49 + x52*x53 + x55 + 1.0*x59;
result[3] = x37*(-n1*x32 + n2*x31 + x24 + x29*x56*x57 + x46*x56 + x5*(-n2*x27 + x25) - x47/((n2)*(n2))) + x38*x52 - x41*x60 + x45 - 6.0*x49 + 3.0*x59;
}
        
static double coder_dgdp(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = 1.0*(1.0*n1 + 1.0*n2)*(n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P))/(n1 + n2);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 1.0*n1 + 1.0*n2;
    double x2 = n1 + n2;
    double x3 = 1.0/x2;
    double x4 = x1*x3;
    double x5 = (*endmember[1].dmu0dP)(T, P);
    double x6 = n1*x0 + n2*x5;
    double x7 = -1.0*x1*x6/((x2)*(x2)) + x3*x6;

result[0] = x0*x4 + x7;
result[1] = x4*x5 + x7;
}
        
static void coder_d3gdn2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2;
    double x2 = 1.0/x1;
    double x3 = x0*x2;
    double x4 = pow(x1, -2);
    double x5 = 1.0*n1 + 1.0*n2;
    double x6 = x4*x5;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = 2.0*n1*x0 + 2.0*n2*x8;
    double x10 = -x4*x9 + x5*x9/((x1)*(x1)*(x1));
    double x11 = x2*x8;
    double x12 = x6*x8;

result[0] = x10 + 2.0*x3 - 2.0*x7;
result[1] = x10 + 1.0*x11 - 1.0*x12 + 1.0*x3 - 1.0*x7;
result[2] = x10 + 2.0*x11 - 2.0*x12;
}
        
static void coder_d4gdn3dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2;
    double x2 = pow(x1, -2);
    double x3 = x0*x2;
    double x4 = pow(x1, -3);
    double x5 = 6.0*x4;
    double x6 = 1.0*n1 + 1.0*n2;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = n1*x0 + n2*x8;
    double x10 = x5*x9 - 6.0*x6*x9/((x1)*(x1)*(x1)*(x1));
    double x11 = x2*x8;
    double x12 = 2.0*x4;
    double x13 = x6*x8;
    double x14 = 4.0*x4;

result[0] = x10 - 6.0*x3 + x5*x7;
result[1] = x10 - 2.0*x11 + x12*x13 + x14*x7 - 4.0*x3;
result[2] = x10 - 4.0*x11 + x12*x7 + x13*x14 - 2.0*x3;
result[3] = x10 - 6.0*x11 + x13*x5;
}
        
static double coder_d2gdt2(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    

result = 1.0*(n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d2mu0dT2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dT2)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d2mu0dTdP)(T, P) + 1.0*n2*(*endmember[1].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d2mu0dTdP)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dTdP)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P));
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d2mu0dP2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d3mu0dT3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT3)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P));
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d3mu0dT2dP)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT2dP)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d3mu0dTdP2)(T, P) + 1.0*n2*(*endmember[1].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d3mu0dTdP2)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dTdP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P));
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];


result[0] = 1.0*(*endmember[0].d3mu0dP3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dP3)(T, P);
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

