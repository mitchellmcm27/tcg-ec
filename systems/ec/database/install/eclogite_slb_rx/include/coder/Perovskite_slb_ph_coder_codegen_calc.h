#include <math.h>


static double coder_g(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = 1.0*n1;
    double x1 = 1.0*n2;
    double x2 = 0.39000000000000001*n3 + x0 + x1;
    double x3 = n1 + n2;
    double x4 = 1.0/(n3 + x3);

result = (65093.525179856115*n1*n3 + x2*(8.3144626181532395*T*(2.0*n3*log(n3*x4) + x0*log(n1*x4) + x1*log(n2*x4) + 1.0*x3*log(x3*x4)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P)))/x2;
    return result;
}
        
static void coder_dgdn(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = 0.39000000000000001*n3;
    double x1 = 1.0*n1;
    double x2 = 1.0*n2;
    double x3 = x1 + x2;
    double x4 = x0 + x3;
    double x5 = 65093.525179856115*n3;
    double x6 = (*endmember[0].mu0)(T, P);
    double x7 = n1*x6;
    double x8 = (*endmember[1].mu0)(T, P);
    double x9 = n2*x8;
    double x10 = (*endmember[2].mu0)(T, P);
    double x11 = n3*x10;
    double x12 = n1 + n2;
    double x13 = n3 + x12;
    double x14 = 1.0/x13;
    double x15 = log(n1*x14);
    double x16 = log(n2*x14);
    double x17 = n3*x14;
    double x18 = 2.0*log(x17);
    double x19 = 1.0*log(x12*x14);
    double x20 = n3*x18 + x1*x15 + x12*x19 + x16*x2;
    double x21 = 8.3144626181532395*T;
    double x22 = x20*x21;
    double x23 = (n1*x5 + x4*(x11 + x22 + x7 + x9))/((x4)*(x4));
    double x24 = -1.0*x23;
    double x25 = 1.0/x4;
    double x26 = -x14*x2;
    double x27 = pow(x13, -2);
    double x28 = 1.0*x13;
    double x29 = -2.0*x17 + x19 + x13*x3*(-x12*x27 + x14)/x12;
    double x30 = x1*x6 + 1.0*x11 + x2*x8 + x22;
    double x31 = -x1*x14;

result[0] = x24 + x25*(x30 + x4*(x21*(1.0*x15 + x26 + x28*(-n1*x27 + x14) + x29) + x6) + x5);
result[1] = x24 + x25*(x30 + x4*(x21*(1.0*x16 + x28*(-n2*x27 + x14) + x29 + x31) + x8));
result[2] = -0.39000000000000001*x23 + x25*(3.2426404210797637*T*x20 + 65093.525179856115*n1 + x0*x10 + x4*(x10 + x21*(2.0*x13*(-n3*x27 + x14) - x14*x3 + x18 + x26 + x31)) + 0.39000000000000001*x7 + 0.39000000000000001*x9);
}
        
static void coder_d2gdn2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = 0.39000000000000001*n3;
    double x1 = 1.0*n1;
    double x2 = 1.0*n2;
    double x3 = x1 + x2;
    double x4 = x0 + x3;
    double x5 = 65093.525179856115*n3;
    double x6 = (*endmember[0].mu0)(T, P);
    double x7 = n1*x6;
    double x8 = (*endmember[1].mu0)(T, P);
    double x9 = n2*x8;
    double x10 = (*endmember[2].mu0)(T, P);
    double x11 = n3*x10;
    double x12 = n1 + n2;
    double x13 = n3 + x12;
    double x14 = 1.0/x13;
    double x15 = log(n1*x14);
    double x16 = log(n2*x14);
    double x17 = n3*x14;
    double x18 = 2.0*log(x17);
    double x19 = 1.0*log(x12*x14);
    double x20 = n3*x18 + x1*x15 + x12*x19 + x16*x2;
    double x21 = 8.3144626181532395*T;
    double x22 = x20*x21;
    double x23 = (n1*x5 + x4*(x11 + x22 + x7 + x9))/((x4)*(x4)*(x4));
    double x24 = 2.0*x23;
    double x25 = pow(x4, -2);
    double x26 = -x14*x2;
    double x27 = pow(x13, -2);
    double x28 = 1.0*x13;
    double x29 = x28*(-n1*x27 + x14);
    double x30 = 1.0/x12;
    double x31 = -x12*x27 + x14;
    double x32 = x30*x31;
    double x33 = x3*x32;
    double x34 = x13*x33 - 2.0*x17 + x19;
    double x35 = T*(1.0*x15 + x26 + x29 + x34);
    double x36 = 8.3144626181532395*x35;
    double x37 = x1*x6 + 1.0*x11 + x2*x8 + x22;
    double x38 = x25*(x37 + x4*(x36 + x6) + x5);
    double x39 = 1.0/x4;
    double x40 = 2.0*x27;
    double x41 = -x40;
    double x42 = pow(x13, -3);
    double x43 = 2.0*x42;
    double x44 = n1*x43;
    double x45 = x2*x27;
    double x46 = x1*x27;
    double x47 = x45 - x46;
    double x48 = 1.0*x14;
    double x49 = 2.0*x13;
    double x50 = 2*x42;
    double x51 = x13*x3;
    double x52 = x30*x51;
    double x53 = n3*x40;
    double x54 = x33 + x53;
    double x55 = x32*x49 + x52*(x12*x50 - 2*x27) + x54 - x31*x51/((x12)*(x12));
    double x56 = x48 + x55;
    double x57 = x21*x4;
    double x58 = -x1*x14;
    double x59 = x28*(-n2*x27 + x14);
    double x60 = 1.0*x16 + x34 + x58 + x59;
    double x61 = x21*x60;
    double x62 = x37 + x4*(x61 + x8);
    double x63 = 1.0*x25;
    double x64 = -1.0*x27;
    double x65 = x13*(x44 + x64) + x47;
    double x66 = -3.0*x14 + x52*(-x27 - x50*(-n1 - n2)) + x54;
    double x67 = x49*(-n3*x27 + x14);
    double x68 = -x14*x3 + x18 + x26 + x58 + x67;
    double x69 = x21*x68;
    double x70 = 1.0*x10 + x69;
    double x71 = 3.2426404210797637*T;
    double x72 = 65093.525179856115*n1 + x0*x10 + x20*x71 + x4*(x10 + x69) + 0.39000000000000001*x7 + 0.39000000000000001*x9;
    double x73 = 0.78000000000000003*x23 - x63*x72;
    double x74 = x25*x62;
    double x75 = n2*x43;
    double x76 = -x45 + x46;

result[0] = x24 - 2.0*x38 + x39*(16.628925236306479*x35 + x57*(x13*(x41 + x44) + x47 + x56 + x29/n1) + 2.0*x6);
result[1] = x24 - 1.0*x38 + x39*(x36 + x57*(-x48 + x55 + x65) + 1.0*x6 + x61 + 1.0*x8) - x62*x63;
result[2] = -0.39000000000000001*x38 + x39*(3.2426404210797637*x35 + x57*(x65 + x66) + 0.39000000000000001*x6 + x70 + 65093.525179856115) + x73;
result[3] = x24 + x39*(16.628925236306479*T*x60 + x57*(x13*(x41 + x75) + x56 + x76 + x59/n2) + 2.0*x8) - 2.0*x74;
result[4] = x39*(x57*(x13*(x64 + x75) + x66 + x76) + x60*x71 + x70 + 0.39000000000000001*x8) + x73 - 0.39000000000000001*x74;
result[5] = 0.30420000000000003*x23 - 0.78000000000000003*x25*x72 + x39*(6.4852808421595274*T*x68 + 0.78000000000000003*x10 + x57*(x13*(4.0*n3*x42 - 4.0*x27) + 2.0*x14 - x27*(-x1 - x2) + x45 + x46 - x53 + x67/n3));
}
        
static void coder_d3gdn3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = 0.39000000000000001*n3;
    double x1 = 1.0*n1;
    double x2 = 1.0*n2;
    double x3 = x1 + x2;
    double x4 = x0 + x3;
    double x5 = 65093.525179856115*n3;
    double x6 = (*endmember[0].mu0)(T, P);
    double x7 = n1*x6;
    double x8 = (*endmember[1].mu0)(T, P);
    double x9 = n2*x8;
    double x10 = (*endmember[2].mu0)(T, P);
    double x11 = n3*x10;
    double x12 = n1 + n2;
    double x13 = n3 + x12;
    double x14 = 1.0/x13;
    double x15 = log(n1*x14);
    double x16 = log(n2*x14);
    double x17 = n3*x14;
    double x18 = 2.0*log(x17);
    double x19 = 1.0*log(x12*x14);
    double x20 = n3*x18 + x1*x15 + x12*x19 + x16*x2;
    double x21 = 8.3144626181532395*T;
    double x22 = x20*x21;
    double x23 = (n1*x5 + x4*(x11 + x22 + x7 + x9))/((x4)*(x4)*(x4)*(x4));
    double x24 = -6.0*x23;
    double x25 = pow(x4, -3);
    double x26 = -x14*x2;
    double x27 = pow(x13, -2);
    double x28 = -1.0*n1*x27 + 1.0*x14;
    double x29 = x13*x28;
    double x30 = 1.0/x12;
    double x31 = -x12*x27 + x14;
    double x32 = x30*x31;
    double x33 = x3*x32;
    double x34 = x13*x33 - 2.0*x17 + x19;
    double x35 = 1.0*x15 + x26 + x29 + x34;
    double x36 = x21*x35;
    double x37 = x1*x6 + 1.0*x11 + x2*x8 + x22;
    double x38 = x25*(x37 + x4*(x36 + x6) + x5);
    double x39 = pow(x4, -2);
    double x40 = 16.628925236306479*T;
    double x41 = 2.0*x27;
    double x42 = -x41;
    double x43 = pow(x13, -3);
    double x44 = 2.0*x43;
    double x45 = n1*x44;
    double x46 = 1.0/n1;
    double x47 = x28*x46;
    double x48 = x2*x27;
    double x49 = x1*x27;
    double x50 = x48 - x49;
    double x51 = 1.0*x14;
    double x52 = 2.0*x32;
    double x53 = -2*x27;
    double x54 = 2*x43;
    double x55 = x12*x54 + x53;
    double x56 = x3*x30;
    double x57 = x55*x56;
    double x58 = pow(x12, -2);
    double x59 = x31*x58;
    double x60 = x3*x59;
    double x61 = n3*x41;
    double x62 = x33 + x61;
    double x63 = x13*x52 + x13*x57 - x13*x60 + x62;
    double x64 = x51 + x63;
    double x65 = T*(x13*x47 + x13*(x42 + x45) + x50 + x64);
    double x66 = 8.3144626181532395*x65;
    double x67 = x39*(x35*x40 + x4*x66 + 2.0*x6);
    double x68 = 1.0/x4;
    double x69 = 6.0*x43;
    double x70 = pow(x13, -4);
    double x71 = 6.0*x70;
    double x72 = -n1*x71;
    double x73 = 2*n1;
    double x74 = x43*x73;
    double x75 = 1.0*x13;
    double x76 = x46*x75;
    double x77 = 4.0*x43;
    double x78 = n2*x44;
    double x79 = -x78;
    double x80 = n3*x77;
    double x81 = -x80;
    double x82 = n1*x77 + x79 + x81;
    double x83 = x47 + x82;
    double x84 = 4.0*x27;
    double x85 = -x84;
    double x86 = 3.0*x13;
    double x87 = x13*x56;
    double x88 = x3*x58;
    double x89 = 2*x13;
    double x90 = x30*x55*x86 + 3.0*x32 - x55*x88*x89 + 2*x57 - x59*x86 - 2*x60 + x87*(-6*x12*x70 + 6*x43) + x3*x31*x89/((x12)*(x12)*(x12));
    double x91 = x85 + x90;
    double x92 = x21*x4;
    double x93 = 1.0*x27;
    double x94 = -x93;
    double x95 = x13*(x45 + x94) + x50;
    double x96 = -x51 + x63 + x95;
    double x97 = x40*x96;
    double x98 = -x27;
    double x99 = x13*(x72 + x77) + x76*(x74 + x98) + x83;
    double x100 = -x1*x14;
    double x101 = -1.0*n2*x27 + 1.0*x14;
    double x102 = x101*x13;
    double x103 = x100 + x102 + 1.0*x16 + x34;
    double x104 = x103*x21;
    double x105 = x37 + x4*(x104 + x8);
    double x106 = 2.0*x25;
    double x107 = x104 + x36 + 1.0*x6 + 1.0*x8 + x92*x96;
    double x108 = 2.0*x39;
    double x109 = -x107*x108 + x24;
    double x110 = -x54*(-n1 - n2) + x98;
    double x111 = x110*x56;
    double x112 = x111*x13 - 3.0*x14 + x62;
    double x113 = x112 + x95;
    double x114 = x110*x13;
    double x115 = 2*n2;
    double x116 = -3*x70*(x115 + x73);
    double x117 = x111 + 2.0*x114*x30 - x114*x88 + x52 + x57 - x60 + x87*(x116 + 4*x43);
    double x118 = x117 + x94;
    double x119 = 3.2426404210797637*T;
    double x120 = x113*x21;
    double x121 = -2.0*n3*x27 + 2.0*x14;
    double x122 = x121*x13;
    double x123 = x100 + x122 - x14*x3 + x18 + x26;
    double x124 = x123*x21;
    double x125 = 1.0*x10 + x124;
    double x126 = x119*x35 + x120*x4 + x125 + 0.39000000000000001*x6 + 65093.525179856115;
    double x127 = 65093.525179856115*n1 + x0*x10 + x119*x20 + x4*(x10 + x124) + 0.39000000000000001*x7 + 0.39000000000000001*x9;
    double x128 = x106*x127 - 2.3399999999999999*x23;
    double x129 = 1.0/n2;
    double x130 = -x48 + x49;
    double x131 = x102*x129 + x13*(x42 + x78) + x130 + x64;
    double x132 = x131*x21;
    double x133 = x13*(x44 + x72) + x82;
    double x134 = x105*x25;
    double x135 = x103*x40 + x132*x4 + 2.0*x8;
    double x136 = 1.0*x39;
    double x137 = x112 + x13*(x78 + x94) + x130;
    double x138 = x137*x21;
    double x139 = x103*x119 + x125 + x138*x4 + 0.39000000000000001*x8;
    double x140 = 0.39000000000000001*x39;
    double x141 = 1.0/n3;
    double x142 = x122*x141 + x13*(x80 + x85) + 2.0*x14 - x27*(-x1 - x2) + x48 + x49 - x61;
    double x143 = x142*x21;
    double x144 = 6.4852808421595274*T;
    double x145 = 2*x111 + x84 + x87*(x116 + x54);
    double x146 = 0.78000000000000003*x39;
    double x147 = x127*x25;
    double x148 = 0.78000000000000003*x10 + x123*x144 + x143*x4;
    double x149 = -x136*x148 + 1.5600000000000001*x147 - 0.91259999999999997*x23;
    double x150 = -n2*x71;
    double x151 = x115*x43;
    double x152 = x129*x75;
    double x153 = -x45;
    double x154 = n2*x77 + x153 + x81;
    double x155 = x101*x129 + x154;

result[0] = x24 + 6.0*x38 - 3.0*x67 + x68*(24.943387854459719*x65 + x92*(x13*(x69 + x72) + x76*(x53 + x74) + x83 + x91 - x29/((n1)*(n1))));
result[1] = x105*x106 + x109 + 4.0*x38 - 1.0*x67 + x68*(x66 + x92*(x42 + x90 + x99) + x97);
result[2] = -x108*x126 + x128 + 1.5600000000000001*x38 - 0.39000000000000001*x67 + x68*(x113*x40 + 3.2426404210797637*x65 + x92*(x118 + x99));
result[3] = x109 + 4.0*x134 - x135*x136 + 2.0*x38 + x68*(x132 + x92*(x133 + x90 + x93) + x97);
result[4] = -x107*x140 - x126*x136 + x128 + 0.78000000000000003*x134 - x136*x139 + 0.78000000000000003*x38 + x68*(x119*x96 + x120 + x138 + x92*(x117 + x133 + x41));
result[5] = -x126*x146 + x149 + 0.30420000000000003*x38 + x68*(x113*x144 + x143 + x92*(x133 + x145));
result[6] = 6.0*x134 - 3.0*x135*x39 + x24 + x68*(24.943387854459719*T*x131 + x92*(x13*(x150 + x69) + x152*(x151 + x53) + x155 + x91 - x102/((n2)*(n2))));
result[7] = -x108*x139 + x128 + 1.5600000000000001*x134 - x135*x140 + x68*(x119*x131 + x137*x40 + x92*(x118 + x13*(x150 + x77) + x152*(x151 + x98) + x155));
result[8] = 0.30420000000000003*x134 - x139*x146 + x149 + x68*(x137*x144 + x143 + x92*(x13*(x150 + x44) + x145 + x154));
result[9] = 0.91260000000000008*x147 - 1.1699999999999999*x148*x39 - 0.35591400000000001*x23 + x68*(9.727921263239292*T*x142 + x92*(8.0*n3*x43 + x121*x141 + 2.0*x13*x141*(n3*x54 + x53) + x13*(-12.0*n3*x70 + 12.0*x43) + x153 - 8.0*x27 - x3*x54 + x79 - x122/((n3)*(n3))));
}
        
static double coder_dgdt(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = n1 + n2;
    double x1 = 1.0/(n3 + x0);

result = 8.3144626181532395*n1*log(n1*x1) + n1*(*endmember[0].dmu0dT)(T, P) + 8.3144626181532395*n2*log(n2*x1) + n2*(*endmember[1].dmu0dT)(T, P) + 16.628925236306479*n3*log(n3*x1) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*x0*log(x0*x1);
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = 1.0/x1;
    double x3 = 8.3144626181532395*n2;
    double x4 = -x2*x3;
    double x5 = n1*x2;
    double x6 = pow(x1, -2);
    double x7 = 8.3144626181532395*x1;
    double x8 = n3*x2;
    double x9 = 8.3144626181532395*n1 + x3;
    double x10 = -16.628925236306479*x8 + 8.3144626181532395*log(x0*x2) + x1*x9*(-x0*x6 + x2)/x0;
    double x11 = -8.3144626181532395*x5;

result[0] = x10 + x4 + x7*(-n1*x6 + x2) + 8.3144626181532395*log(x5) + (*endmember[0].dmu0dT)(T, P);
result[1] = x10 + x11 + x7*(-n2*x6 + x2) + 8.3144626181532395*log(n2*x2) + (*endmember[1].dmu0dT)(T, P);
result[2] = 16.628925236306479*x1*(-n3*x6 + x2) + x11 - x2*x9 + x4 + 16.628925236306479*log(x8) + (*endmember[2].dmu0dT)(T, P);
}
        
static void coder_d3gdn2dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = pow(x1, -2);
    double x3 = 16.628925236306479*x2;
    double x4 = -x3;
    double x5 = pow(x1, -3);
    double x6 = 16.628925236306479*x5;
    double x7 = n1*x6;
    double x8 = 1.0/x1;
    double x9 = 8.3144626181532395*x1;
    double x10 = 8.3144626181532395*n2;
    double x11 = x10*x2;
    double x12 = 8.3144626181532395*n1;
    double x13 = x12*x2;
    double x14 = x11 - x13;
    double x15 = 8.3144626181532395*x8;
    double x16 = 1.0/x0;
    double x17 = -x0*x2 + x8;
    double x18 = x16*x17;
    double x19 = 16.628925236306479*x1;
    double x20 = 2*x5;
    double x21 = x10 + x12;
    double x22 = x1*x21;
    double x23 = x16*x22;
    double x24 = n3*x3;
    double x25 = x18*x21 + x24;
    double x26 = x18*x19 + x23*(x0*x20 - 2*x2) + x25 - x17*x22/((x0)*(x0));
    double x27 = x15 + x26;
    double x28 = -8.3144626181532395*x2;
    double x29 = x1*(x28 + x7) + x14;
    double x30 = x23*(-x2 - x20*(-n1 - n2)) + x25 - 24.943387854459719*x8;
    double x31 = n2*x6;
    double x32 = -x11 + x13;

result[0] = x1*(x4 + x7) + x14 + x27 + x9*(-n1*x2 + x8)/n1;
result[1] = -x15 + x26 + x29;
result[2] = x29 + x30;
result[3] = x1*(x31 + x4) + x27 + x32 + x9*(-n2*x2 + x8)/n2;
result[4] = x1*(x28 + x31) + x30 + x32;
result[5] = x1*(33.257850472612958*n3*x5 - 33.257850472612958*x2) + x11 + x13 - x2*(-x10 - x12) - x24 + 16.628925236306479*x8 + x19*(-n3*x2 + x8)/n3;
}
        
static void coder_d4gdn3dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = pow(x1, -3);
    double x3 = 49.886775708919437*x2;
    double x4 = pow(x1, -4);
    double x5 = 49.886775708919437*x4;
    double x6 = -n1*x5;
    double x7 = pow(x1, -2);
    double x8 = -2*x7;
    double x9 = 2*n1;
    double x10 = x2*x9;
    double x11 = 8.3144626181532395/n1;
    double x12 = x1*x11;
    double x13 = 1.0/x1;
    double x14 = -n1*x7 + x13;
    double x15 = 8.3144626181532395*x1;
    double x16 = 33.257850472612958*x2;
    double x17 = 16.628925236306479*x2;
    double x18 = -n2*x17;
    double x19 = -n3*x16;
    double x20 = n1*x16 + x18 + x19;
    double x21 = x11*x14 + x20;
    double x22 = 33.257850472612958*x7;
    double x23 = 1.0/x0;
    double x24 = -x0*x7 + x13;
    double x25 = x23*x24;
    double x26 = 8.3144626181532395*n1 + 8.3144626181532395*n2;
    double x27 = pow(x0, -2);
    double x28 = x24*x27;
    double x29 = x26*x28;
    double x30 = 2*x2;
    double x31 = x0*x30 + x8;
    double x32 = x23*x26;
    double x33 = x31*x32;
    double x34 = 24.943387854459719*x1;
    double x35 = x1*x32;
    double x36 = x26*x27;
    double x37 = 2*x1;
    double x38 = x23*x31*x34 + 24.943387854459719*x25 - x28*x34 - 2*x29 - x31*x36*x37 + 2*x33 + x35*(-6*x0*x4 + 6*x2) + x24*x26*x37/((x0)*(x0)*(x0));
    double x39 = -x22 + x38;
    double x40 = 16.628925236306479*x7;
    double x41 = -x7;
    double x42 = x1*(x16 + x6) + x12*(x10 + x41) + x21;
    double x43 = 8.3144626181532395*x7;
    double x44 = -x30*(-n1 - n2) + x41;
    double x45 = x32*x44;
    double x46 = x1*x44;
    double x47 = 2*n2;
    double x48 = -3*x4*(x47 + x9);
    double x49 = 16.628925236306479*x23*x46 + 16.628925236306479*x25 - x29 + x33 + x35*(4*x2 + x48) - x36*x46 + x45;
    double x50 = -x43 + x49;
    double x51 = x1*(x17 + x6) + x20;
    double x52 = x22 + x35*(x30 + x48) + 2*x45;
    double x53 = -n2*x5;
    double x54 = x2*x47;
    double x55 = 8.3144626181532395/n2;
    double x56 = x1*x55;
    double x57 = -n2*x7 + x13;
    double x58 = -n1*x17;
    double x59 = n2*x16 + x19 + x58;
    double x60 = x55*x57 + x59;
    double x61 = -n3*x7 + x13;
    double x62 = 16.628925236306479/n3;

result[0] = x1*(x3 + x6) + x12*(x10 + x8) + x21 + x39 - x14*x15/((n1)*(n1));
result[1] = x38 - x40 + x42;
result[2] = x42 + x50;
result[3] = x38 + x43 + x51;
result[4] = x40 + x49 + x51;
result[5] = x51 + x52;
result[6] = x1*(x3 + x53) + x39 + x56*(x54 + x8) + x60 - x15*x57/((n2)*(n2));
result[7] = x1*(x16 + x53) + x50 + x56*(x41 + x54) + x60;
result[8] = x1*(x17 + x53) + x52 + x59;
result[9] = 66.515700945225916*n3*x2 + x1*x62*(n3*x30 + x8) + x1*(-99.773551417838874*n3*x4 + 99.773551417838874*x2) + x18 - x26*x30 + x58 + x61*x62 - 66.515700945225916*x7 - 16.628925236306479*x1*x61/((n3)*(n3));
}
        
static double coder_dgdp(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].dmu0dP)(T, P);
result[1] = (*endmember[1].dmu0dP)(T, P);
result[2] = (*endmember[2].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d4gdn3dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d2gdt2(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P);
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = (*endmember[1].d2mu0dT2)(T, P);
result[2] = (*endmember[2].d2mu0dT2)(T, P);
}
        
static void coder_d4gdn2dt2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d5gdn3dt2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d2gdtdp(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d2mu0dTdP)(T, P) + n2*(*endmember[1].d2mu0dTdP)(T, P) + n3*(*endmember[2].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d2mu0dTdP)(T, P);
result[1] = (*endmember[1].d2mu0dTdP)(T, P);
result[2] = (*endmember[2].d2mu0dTdP)(T, P);
}
        
static void coder_d4gdn2dtdp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d5gdn3dtdp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d2gdp2(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P);
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d2mu0dP2)(T, P);
result[1] = (*endmember[1].d2mu0dP2)(T, P);
result[2] = (*endmember[2].d2mu0dP2)(T, P);
}
        
static void coder_d4gdn2dp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d5gdn3dp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d3gdt3(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P);
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d3mu0dT3)(T, P);
result[1] = (*endmember[1].d3mu0dT3)(T, P);
result[2] = (*endmember[2].d3mu0dT3)(T, P);
}
        
static void coder_d5gdn2dt3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d6gdn3dt3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d3gdt2dp(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P);
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d3mu0dT2dP)(T, P);
result[1] = (*endmember[1].d3mu0dT2dP)(T, P);
result[2] = (*endmember[2].d3mu0dT2dP)(T, P);
}
        
static void coder_d5gdn2dt2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d6gdn3dt2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d3gdtdp2(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d3mu0dTdP2)(T, P) + n2*(*endmember[1].d3mu0dTdP2)(T, P) + n3*(*endmember[2].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d3mu0dTdP2)(T, P);
result[1] = (*endmember[1].d3mu0dTdP2)(T, P);
result[2] = (*endmember[2].d3mu0dTdP2)(T, P);
}
        
static void coder_d5gdn2dtdp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d6gdn3dtdp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_d3gdp3(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P);
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = (*endmember[0].d3mu0dP3)(T, P);
result[1] = (*endmember[1].d3mu0dP3)(T, P);
result[2] = (*endmember[2].d3mu0dP3)(T, P);
}
        
static void coder_d5gdn2dp3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
}
        
static void coder_d6gdn3dp3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 0;
result[1] = 0;
result[2] = 0;
result[3] = 0;
result[4] = 0;
result[5] = 0;
result[6] = 0;
result[7] = 0;
result[8] = 0;
result[9] = 0;
}
        
static double coder_s(double T, double P, double n[3]) {
    double result = -coder_dgdt(T, P, n);
    return result;
}

static double coder_v(double T, double P, double n[3]) {
    double result = coder_dgdp(T, P, n);
    return result;
}

static double coder_cv(double T, double P, double n[3]) {
    double result = -T*coder_d2gdt2(T, P, n);
    double dvdt = coder_d2gdtdp(T, P, n);
    double dvdp = coder_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P, double n[3]) {
    double result = -T*coder_d2gdt2(T, P, n);
    return result;
}

static double coder_dcpdt(double T, double P, double n[3]) {
    double result = -T*coder_d3gdt3(T, P, n) - coder_d2gdt2(T, P, n);
    return result;
}

static double coder_alpha(double T, double P, double n[3]) {
    double result = coder_d2gdtdp(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_beta(double T, double P, double n[3]) {
    double result = -coder_d2gdp2(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_K(double T, double P, double n[3]) {
    double result = -coder_dgdp(T, P, n)/coder_d2gdp2(T, P, n);
    return result;
}

static double coder_Kp(double T, double P, double n[3]) {
    double result = coder_dgdp(T, P, n);
    result *= coder_d3gdp3(T, P, n);
    result /= pow(coder_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}

