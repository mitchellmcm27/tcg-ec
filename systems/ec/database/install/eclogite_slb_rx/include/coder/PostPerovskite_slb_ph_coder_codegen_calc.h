#include <math.h>


static double coder_g(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = 1.0/x1;

result = 1.0*x2*(60000.0*n1*n3 + 1.0*x1*(8.3144626181532395*T*(1.0*n1*log(n1*x2) + 1.0*n2*log(n2*x2) + 2.0*n3*log(n3*x2) + 1.0*x0*log(x0*x2)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P)));
    return result;
}
        
static void coder_dgdn(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = pow(x1, -2);
    double x3 = 60000.0*n3;
    double x4 = (*endmember[0].mu0)(T, P);
    double x5 = (*endmember[1].mu0)(T, P);
    double x6 = (*endmember[2].mu0)(T, P);
    double x7 = 1.0/x1;
    double x8 = n1*x7;
    double x9 = 1.0*log(x8);
    double x10 = log(n2*x7);
    double x11 = 1.0*n2;
    double x12 = n3*x7;
    double x13 = 2.0*log(x12);
    double x14 = 1.0*log(x0*x7);
    double x15 = 8.3144626181532395*T;
    double x16 = x15*(n1*x9 + n3*x13 + x0*x14 + x10*x11);
    double x17 = 1.0*x1;
    double x18 = -1.0*x2*(n1*x3 + x17*(n1*x4 + n2*x5 + n3*x6 + x16));
    double x19 = 1.0*n3;
    double x20 = 1.0*n1;
    double x21 = x11 + x20;
    double x22 = x19 + x21;
    double x23 = -x11*x7;
    double x24 = -2.0*x12 + x14 + x1*x21*(-x0*x2 + x7)/x0;
    double x25 = x11*x5 + x16 + x19*x6 + x20*x4;
    double x26 = 1.0*x7;
    double x27 = -1.0*x8;

result[0] = x18 + x26*(x22*(x15*(x17*(-n1*x2 + x7) + x23 + x24 + x9) + x4) + x25 + x3);
result[1] = x18 + x26*(x22*(x15*(1.0*x10 + x17*(-n2*x2 + x7) + x24 + x27) + x5) + x25);
result[2] = x18 + x26*(60000.0*n1 + x22*(x15*(2.0*x1*(-n3*x2 + x7) + x13 - x21*x7 + x23 + x27) + x6) + x25);
}
        
static void coder_d2gdn2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = 60000.0*n3;
    double x1 = (*endmember[0].mu0)(T, P);
    double x2 = (*endmember[1].mu0)(T, P);
    double x3 = (*endmember[2].mu0)(T, P);
    double x4 = n1 + n2;
    double x5 = n3 + x4;
    double x6 = 1.0/x5;
    double x7 = n1*x6;
    double x8 = 1.0*log(x7);
    double x9 = log(n2*x6);
    double x10 = 1.0*n2;
    double x11 = n3*x6;
    double x12 = 2.0*log(x11);
    double x13 = 1.0*log(x4*x6);
    double x14 = 8.3144626181532395*T;
    double x15 = x14*(n1*x8 + n3*x12 + x10*x9 + x13*x4);
    double x16 = 1.0*x5;
    double x17 = pow(x5, -3);
    double x18 = 2.0*x17;
    double x19 = x18*(n1*x0 + x16*(n1*x1 + n2*x2 + n3*x3 + x15));
    double x20 = 1.0*n3;
    double x21 = 1.0*n1;
    double x22 = x10 + x21;
    double x23 = x20 + x22;
    double x24 = -x10*x6;
    double x25 = pow(x5, -2);
    double x26 = n1*x25;
    double x27 = x16*(-x26 + x6);
    double x28 = 1.0/x4;
    double x29 = -x25*x4 + x6;
    double x30 = x28*x29;
    double x31 = x22*x30;
    double x32 = -2.0*x11 + x13 + x31*x5;
    double x33 = T*(x24 + x27 + x32 + x8);
    double x34 = 8.3144626181532395*x33;
    double x35 = x1*x21 + x10*x2 + x15 + x20*x3;
    double x36 = x0 + x23*(x1 + x34) + x35;
    double x37 = 2.0*x25;
    double x38 = -x37;
    double x39 = n1*x18;
    double x40 = x10*x25;
    double x41 = 1.0*x26;
    double x42 = x40 - x41;
    double x43 = 1.0*x6;
    double x44 = 2.0*x5;
    double x45 = 2*x17;
    double x46 = x22*x5;
    double x47 = x28*x46;
    double x48 = n3*x37;
    double x49 = x31 + x48;
    double x50 = -x29*x46/((x4)*(x4)) + x30*x44 + x47*(-2*x25 + x4*x45) + x49;
    double x51 = x43 + x50;
    double x52 = x14*x23;
    double x53 = 1.0*x25;
    double x54 = -x53;
    double x55 = x42 + x5*(x39 + x54);
    double x56 = -1.0*x7;
    double x57 = x16*(-n2*x25 + x6);
    double x58 = x32 + x56 + x57 + 1.0*x9;
    double x59 = x14*x58;
    double x60 = 1.0*x2 + x59;
    double x61 = 1.0*x1 + x34;
    double x62 = x23*(x2 + x59) + x35;
    double x63 = -x53*x62;
    double x64 = x19 - x36*x53;
    double x65 = x47*(-x25 - x45*(-n1 - n2)) + x49 - 3.0*x6;
    double x66 = x44*(-n3*x25 + x6);
    double x67 = x12 - x22*x6 + x24 + x56 + x66;
    double x68 = x14*x67;
    double x69 = 1.0*x3 + x68;
    double x70 = 60000.0*n1 + x23*(x3 + x68) + x35;
    double x71 = -x53*x70;
    double x72 = 16.628925236306479*T;
    double x73 = n2*x18;
    double x74 = -x40 + x41;

result[0] = x19 - x36*x37 + x43*(2.0*x1 + 16.628925236306479*x33 + x52*(x42 + x5*(x38 + x39) + x51 + x27/n1));
result[1] = x43*(x52*(-x43 + x50 + x55) + x60 + x61) + x63 + x64;
result[2] = x43*(x52*(x55 + x65) + x61 + x69 + 60000.0) + x64 + x71;
result[3] = x19 - x37*x62 + x43*(2.0*x2 + x52*(x5*(x38 + x73) + x51 + x74 + x57/n2) + x58*x72);
result[4] = x19 + x43*(x52*(x5*(x54 + x73) + x65 + x74) + x60 + x69) + x63 + x71;
result[5] = x19 - x37*x70 + x43*(2.0*x3 + x52*(-x25*(-x10 - x21) + x40 + x41 - x48 + x5*(4.0*n3*x17 - 4.0*x25) + 2.0*x6 + x66/n3) + x67*x72);
}
        
static void coder_d3gdn3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = 60000.0*n3;
    double x1 = (*endmember[0].mu0)(T, P);
    double x2 = (*endmember[1].mu0)(T, P);
    double x3 = (*endmember[2].mu0)(T, P);
    double x4 = n1 + n2;
    double x5 = n3 + x4;
    double x6 = 1.0/x5;
    double x7 = log(n1*x6);
    double x8 = 1.0*n1;
    double x9 = log(n2*x6);
    double x10 = 1.0*n2;
    double x11 = n3*x6;
    double x12 = 2.0*log(x11);
    double x13 = 1.0*log(x4*x6);
    double x14 = 8.3144626181532395*T;
    double x15 = x14*(n3*x12 + x10*x9 + x13*x4 + x7*x8);
    double x16 = 1.0*x5;
    double x17 = pow(x5, -4);
    double x18 = 6.0*x17;
    double x19 = -x18*(n1*x0 + x16*(n1*x1 + n2*x2 + n3*x3 + x15));
    double x20 = 1.0*n3;
    double x21 = x10 + x8;
    double x22 = x20 + x21;
    double x23 = 1.0*x6;
    double x24 = -n2*x23;
    double x25 = pow(x5, -2);
    double x26 = -1.0*n1*x25 + 1.0*x6;
    double x27 = x26*x5;
    double x28 = 1.0/x4;
    double x29 = -x25*x4 + x6;
    double x30 = x28*x29;
    double x31 = x21*x30;
    double x32 = -2.0*x11 + x13 + x31*x5;
    double x33 = x24 + x27 + x32 + 1.0*x7;
    double x34 = x14*x33;
    double x35 = x1*x8 + x10*x2 + x15 + x20*x3;
    double x36 = x0 + x22*(x1 + x34) + x35;
    double x37 = pow(x5, -3);
    double x38 = 6.0*x37;
    double x39 = 16.628925236306479*T;
    double x40 = 2.0*x25;
    double x41 = -x40;
    double x42 = 2.0*x37;
    double x43 = n1*x42;
    double x44 = 1.0/n1;
    double x45 = x26*x44;
    double x46 = x10*x25;
    double x47 = x25*x8;
    double x48 = x46 - x47;
    double x49 = 2.0*x30;
    double x50 = -2*x25;
    double x51 = 2*x37;
    double x52 = x4*x51 + x50;
    double x53 = x21*x28;
    double x54 = x52*x53;
    double x55 = pow(x4, -2);
    double x56 = x29*x55;
    double x57 = x21*x56;
    double x58 = n3*x40;
    double x59 = x31 + x58;
    double x60 = x49*x5 + x5*x54 - x5*x57 + x59;
    double x61 = x23 + x60;
    double x62 = T*(x45*x5 + x48 + x5*(x41 + x43) + x61);
    double x63 = 8.3144626181532395*x62;
    double x64 = 2.0*x1 + x22*x63 + x33*x39;
    double x65 = 3.0*x25;
    double x66 = -n1*x18;
    double x67 = 2*n1;
    double x68 = x37*x67;
    double x69 = x16*x44;
    double x70 = 4.0*x37;
    double x71 = n2*x42;
    double x72 = -x71;
    double x73 = n3*x70;
    double x74 = -x73;
    double x75 = n1*x70 + x72 + x74;
    double x76 = x45 + x75;
    double x77 = 4.0*x25;
    double x78 = -x77;
    double x79 = 3.0*x5;
    double x80 = x5*x53;
    double x81 = x21*x55;
    double x82 = 2*x5;
    double x83 = x21*x29*x82/((x4)*(x4)*(x4)) + x28*x52*x79 + 3.0*x30 - x52*x81*x82 + 2*x54 - x56*x79 - 2*x57 + x80*(-6*x17*x4 + 6*x37);
    double x84 = x78 + x83;
    double x85 = x14*x22;
    double x86 = 1.0*x25;
    double x87 = -x86;
    double x88 = x48 + x5*(x43 + x87);
    double x89 = -x23 + x60 + x88;
    double x90 = x39*x89;
    double x91 = -x25;
    double x92 = x5*(x66 + x70) + x69*(x68 + x91) + x76;
    double x93 = -n1*x23;
    double x94 = -n2*x25 + x6;
    double x95 = x16*x94;
    double x96 = x32 + 1.0*x9 + x93 + x95;
    double x97 = x14*x96;
    double x98 = x22*(x2 + x97) + x35;
    double x99 = x42*x98;
    double x100 = x14*x89;
    double x101 = 1.0*x2 + x97;
    double x102 = 1.0*x1 + x34;
    double x103 = x100*x22 + x101 + x102;
    double x104 = -x103*x40 + x19;
    double x105 = x36*x70 - x64*x86;
    double x106 = -x51*(-n1 - n2) + x91;
    double x107 = x106*x53;
    double x108 = x107*x5 + x59 - 3.0*x6;
    double x109 = x108 + x88;
    double x110 = x109*x39;
    double x111 = x106*x5;
    double x112 = 2*n2;
    double x113 = -3*x17*(x112 + x67);
    double x114 = x107 + 2.0*x111*x28 - x111*x81 + x49 + x54 - x57 + x80*(x113 + 4*x37);
    double x115 = x114 + x87;
    double x116 = -2.0*n3*x25 + 2.0*x6;
    double x117 = x116*x5;
    double x118 = x117 + x12 - x21*x6 + x24 + x93;
    double x119 = x118*x14;
    double x120 = 60000.0*n1 + x22*(x119 + x3) + x35;
    double x121 = x120*x42;
    double x122 = x109*x14;
    double x123 = x119 + 1.0*x3;
    double x124 = x102 + x122*x22 + x123 + 60000.0;
    double x125 = -x124*x40 + x19;
    double x126 = 1.0/n2;
    double x127 = -x46 + x47;
    double x128 = x126*x95 + x127 + x5*(x41 + x71) + x61;
    double x129 = x128*x14;
    double x130 = x5*(x42 + x66) + x75;
    double x131 = x36*x42;
    double x132 = x129*x22 + 2.0*x2 + x39*x96;
    double x133 = -x132*x86 + x70*x98;
    double x134 = x108 + x127 + x5*(x71 + x87);
    double x135 = x134*x14;
    double x136 = x101 + x123 + x135*x22;
    double x137 = 1.0/n3;
    double x138 = x117*x137 - x25*(-x10 - x8) + x46 + x47 + x5*(x73 + x78) - x58 + 2.0*x6;
    double x139 = x138*x14;
    double x140 = 2*x107 + x77 + x80*(x113 + x51);
    double x141 = x118*x39 + x139*x22 + 2.0*x3;
    double x142 = x120*x70 - x141*x86;
    double x143 = 24.943387854459719*T;
    double x144 = -n2*x18;
    double x145 = x112*x37;
    double x146 = x126*x16;
    double x147 = -x43;
    double x148 = n2*x70 + x147 + x74;
    double x149 = 1.0*x126*x94 + x148;
    double x150 = x134*x39;
    double x151 = -x136*x40 + x19;

result[0] = x19 + x23*(24.943387854459719*x62 + x85*(x5*(x38 + x66) + x69*(x50 + x68) + x76 + x84 - x27/((n1)*(n1)))) + x36*x38 - x64*x65;
result[1] = x104 + x105 + x23*(x63 + x85*(x41 + x83 + x92) + x90) + x99;
result[2] = x105 + x121 + x125 + x23*(x110 + x63 + x85*(x115 + x92));
result[3] = x104 + x131 + x133 + x23*(x129 + x85*(x130 + x83 + x86) + x90);
result[4] = -x103*x86 + x121 - x124*x86 + x131 - x136*x86 + x19 + x23*(x100 + x122 + x135 + x85*(x114 + x130 + x40)) + x99;
result[5] = x125 + x131 + x142 + x23*(x110 + x139 + x85*(x130 + x140));
result[6] = -x132*x65 + x19 + x23*(x128*x143 + x85*(x146*(x145 + x50) + x149 + x5*(x144 + x38) + x84 - x95/((n2)*(n2)))) + x38*x98;
result[7] = x121 + x133 + x151 + x23*(x129 + x150 + x85*(x115 + x146*(x145 + x91) + x149 + x5*(x144 + x70)));
result[8] = x142 + x151 + x23*(x139 + x150 + x85*(x140 + x148 + x5*(x144 + x42))) + x99;
result[9] = x120*x38 - x141*x65 + x19 + x23*(x138*x143 + x85*(8.0*n3*x37 + x116*x137 + 2.0*x137*x5*(n3*x51 + x50) + x147 - x21*x51 - 8.0*x25 + x5*(-12.0*n3*x17 + 12.0*x37) + x72 - x117/((n3)*(n3))));
}
        
static double coder_dgdt(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = n1 + n2;
    double x1 = 1.0/(n3 + x0);

result = 1.0*x1*(1.0*n1 + 1.0*n2 + 1.0*n3)*(8.3144626181532395*n1*log(n1*x1) + n1*(*endmember[0].dmu0dT)(T, P) + 8.3144626181532395*n2*log(n2*x1) + n2*(*endmember[1].dmu0dT)(T, P) + 16.628925236306479*n3*log(n3*x1) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*x0*log(x0*x1));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = 1.0/x1;
    double x3 = n2*x2;
    double x4 = -8.3144626181532395*x3;
    double x5 = n1*x2;
    double x6 = 8.3144626181532395*log(x5);
    double x7 = pow(x1, -2);
    double x8 = 8.3144626181532395*x1;
    double x9 = (*endmember[0].dmu0dT)(T, P);
    double x10 = 8.3144626181532395*log(x0*x2);
    double x11 = n3*x2;
    double x12 = 8.3144626181532395*n2;
    double x13 = 8.3144626181532395*n1 + x12;
    double x14 = x10 - 16.628925236306479*x11 + x1*x13*(-x0*x7 + x2)/x0;
    double x15 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x16 = 1.0*x2;
    double x17 = x15*x16;
    double x18 = (*endmember[1].dmu0dT)(T, P);
    double x19 = (*endmember[2].dmu0dT)(T, P);
    double x20 = log(x3);
    double x21 = 16.628925236306479*log(x11);
    double x22 = n1*x6 + n1*x9 + n2*x18 + n3*x19 + n3*x21 + x0*x10 + x12*x20;
    double x23 = -1.0*x15*x22*x7 + x16*x22;
    double x24 = -8.3144626181532395*x5;

result[0] = x17*(x14 + x4 + x6 + x8*(-n1*x7 + x2) + x9) + x23;
result[1] = x17*(x14 + x18 + 8.3144626181532395*x20 + x24 + x8*(-n2*x7 + x2)) + x23;
result[2] = x17*(16.628925236306479*x1*(-n3*x7 + x2) - x13*x2 + x19 + x21 + x24 + x4) + x23;
}
        
static void coder_d3gdn2dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = 1.0/x1;
    double x3 = 8.3144626181532395*n2;
    double x4 = -x2*x3;
    double x5 = n1*x2;
    double x6 = 8.3144626181532395*log(x5);
    double x7 = pow(x1, -2);
    double x8 = n1*x7;
    double x9 = 8.3144626181532395*x1;
    double x10 = x9*(x2 - x8);
    double x11 = (*endmember[0].dmu0dT)(T, P);
    double x12 = 8.3144626181532395*log(x0*x2);
    double x13 = n3*x2;
    double x14 = 8.3144626181532395*n1;
    double x15 = x14 + x3;
    double x16 = 1.0/x0;
    double x17 = -x0*x7 + x2;
    double x18 = x16*x17;
    double x19 = x15*x18;
    double x20 = x1*x19 + x12 - 16.628925236306479*x13;
    double x21 = x10 + x11 + x20 + x4 + x6;
    double x22 = x2*x21;
    double x23 = 16.628925236306479*x7;
    double x24 = -x23;
    double x25 = pow(x1, -3);
    double x26 = 16.628925236306479*x25;
    double x27 = n1*x26;
    double x28 = x3*x7;
    double x29 = 8.3144626181532395*x8;
    double x30 = x28 - x29;
    double x31 = 8.3144626181532395*x2;
    double x32 = 16.628925236306479*x1;
    double x33 = 2*x25;
    double x34 = x1*x15;
    double x35 = x16*x34;
    double x36 = n3*x23;
    double x37 = x19 + x36;
    double x38 = x18*x32 + x35*(x0*x33 - 2*x7) + x37 - x17*x34/((x0)*(x0));
    double x39 = x31 + x38;
    double x40 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x41 = 1.0*x2;
    double x42 = x40*x41;
    double x43 = x40*x7;
    double x44 = x21*x43;
    double x45 = (*endmember[1].dmu0dT)(T, P);
    double x46 = (*endmember[2].dmu0dT)(T, P);
    double x47 = log(n2*x2);
    double x48 = 16.628925236306479*log(x13);
    double x49 = 2.0*n1*x11 + 2.0*n1*x6 + 2.0*n2*x45 + 2.0*n3*x46 + 2.0*n3*x48 + 2.0*x0*x12 + 2.0*x3*x47;
    double x50 = x25*x40*x49 - x49*x7;
    double x51 = -8.3144626181532395*x7;
    double x52 = x1*(x27 + x51) + x30;
    double x53 = 1.0*x22 - 1.0*x44 + x50;
    double x54 = -8.3144626181532395*x5;
    double x55 = x9*(-n2*x7 + x2);
    double x56 = x20 + x45 + 8.3144626181532395*x47 + x54 + x55;
    double x57 = 1.0*x43;
    double x58 = x41*x56 - x56*x57;
    double x59 = -24.943387854459719*x2 + x35*(-x33*(-n1 - n2) - x7) + x37;
    double x60 = x32*(-n3*x7 + x2);
    double x61 = -x15*x2 + x4 + x46 + x48 + x54 + x60;
    double x62 = x41*x61 - x57*x61;
    double x63 = 2.0*x2;
    double x64 = n2*x26;
    double x65 = -x28 + x29;
    double x66 = 2.0*x43;

result[0] = 2.0*x22 + x42*(x1*(x24 + x27) + x30 + x39 + x10/n1) - 2.0*x44 + x50;
result[1] = x42*(-x31 + x38 + x52) + x53 + x58;
result[2] = x42*(x52 + x59) + x53 + x62;
result[3] = x42*(x1*(x24 + x64) + x39 + x65 + x55/n2) + x50 + x56*x63 - x56*x66;
result[4] = x42*(x1*(x51 + x64) + x59 + x65) + x50 + x58 + x62;
result[5] = x42*(x1*(33.257850472612958*n3*x25 - 33.257850472612958*x7) + 16.628925236306479*x2 + x28 + x29 - x36 - x7*(-x14 - x3) + x60/n3) + x50 + x61*x63 - x61*x66;
}
        
static void coder_d4gdn3dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2;
    double x1 = n3 + x0;
    double x2 = 1.0/x1;
    double x3 = pow(x1, -2);
    double x4 = 16.628925236306479*x3;
    double x5 = -x4;
    double x6 = pow(x1, -3);
    double x7 = 16.628925236306479*x6;
    double x8 = n1*x7;
    double x9 = 1.0/n1;
    double x10 = -8.3144626181532395*n1*x3 + 8.3144626181532395*x2;
    double x11 = x1*x10;
    double x12 = 8.3144626181532395*n2;
    double x13 = x12*x3;
    double x14 = 8.3144626181532395*n1;
    double x15 = x14*x3;
    double x16 = x13 - x15;
    double x17 = 8.3144626181532395*x2;
    double x18 = 1.0/x0;
    double x19 = -x0*x3 + x2;
    double x20 = x18*x19;
    double x21 = 16.628925236306479*x20;
    double x22 = -2*x3;
    double x23 = 2*x6;
    double x24 = x0*x23 + x22;
    double x25 = x12 + x14;
    double x26 = x18*x25;
    double x27 = x24*x26;
    double x28 = pow(x0, -2);
    double x29 = x19*x28;
    double x30 = x25*x29;
    double x31 = n3*x4;
    double x32 = x20*x25;
    double x33 = x31 + x32;
    double x34 = x1*x21 + x1*x27 - x1*x30 + x33;
    double x35 = x17 + x34;
    double x36 = x1*(x5 + x8) + x11*x9 + x16 + x35;
    double x37 = x2*x36;
    double x38 = -n2*x17;
    double x39 = log(n1*x2);
    double x40 = (*endmember[0].dmu0dT)(T, P);
    double x41 = 8.3144626181532395*log(x0*x2);
    double x42 = n3*x2;
    double x43 = x1*x32 + x41 - 16.628925236306479*x42;
    double x44 = x11 + x38 + 8.3144626181532395*x39 + x40 + x43;
    double x45 = x3*x44;
    double x46 = 49.886775708919437*x6;
    double x47 = pow(x1, -4);
    double x48 = 49.886775708919437*x47;
    double x49 = -n1*x48;
    double x50 = 2*n1;
    double x51 = x50*x6;
    double x52 = 8.3144626181532395*x1;
    double x53 = x52*x9;
    double x54 = 33.257850472612958*x6;
    double x55 = n2*x7;
    double x56 = -x55;
    double x57 = n3*x54;
    double x58 = -x57;
    double x59 = n1*x54 + x56 + x58;
    double x60 = x10*x9 + x59;
    double x61 = 33.257850472612958*x3;
    double x62 = -x61;
    double x63 = 24.943387854459719*x1;
    double x64 = x1*x26;
    double x65 = x25*x28;
    double x66 = 2*x1;
    double x67 = x18*x24*x63 + 24.943387854459719*x20 - x24*x65*x66 + 2*x27 - x29*x63 - 2*x30 + x64*(-6*x0*x47 + 6*x6) + x19*x25*x66/((x0)*(x0)*(x0));
    double x68 = x62 + x67;
    double x69 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x70 = 1.0*x2;
    double x71 = x69*x70;
    double x72 = 6.0*x6;
    double x73 = x44*x69;
    double x74 = x3*x69;
    double x75 = x36*x74;
    double x76 = (*endmember[1].dmu0dT)(T, P);
    double x77 = (*endmember[2].dmu0dT)(T, P);
    double x78 = log(n2*x2);
    double x79 = 16.628925236306479*log(x42);
    double x80 = n1*x40 + n2*x76 + n3*x77 + n3*x79 + x0*x41 + x12*x78 + x14*x39;
    double x81 = -6.0*x47*x69*x80 + x72*x80;
    double x82 = -x3;
    double x83 = x1*(x49 + x54) + x53*(x51 + x82) + x60;
    double x84 = 8.3144626181532395*x3;
    double x85 = -x84;
    double x86 = x1*(x8 + x85) + x16;
    double x87 = -x17 + x34 + x86;
    double x88 = 2.0*x2;
    double x89 = 2.0*x3;
    double x90 = x69*x89;
    double x91 = x81 + x87*x88 - x87*x90;
    double x92 = -n1*x17;
    double x93 = -8.3144626181532395*n2*x3 + 8.3144626181532395*x2;
    double x94 = x1*x93;
    double x95 = x43 + x76 + 8.3144626181532395*x78 + x92 + x94;
    double x96 = 2.0*x6;
    double x97 = x69*x95;
    double x98 = -x89*x95 + x96*x97;
    double x99 = 4.0*x6;
    double x100 = 1.0*x37 - 4.0*x45 + x73*x99 - 1.0*x75;
    double x101 = -x23*(-n1 - n2) + x82;
    double x102 = x101*x26;
    double x103 = x1*x101;
    double x104 = 2*n2;
    double x105 = -3*x47*(x104 + x50);
    double x106 = x102 + 16.628925236306479*x103*x18 - x103*x65 + x21 + x27 - x30 + x64*(x105 + 4*x6);
    double x107 = x106 + x85;
    double x108 = x1*x102 - 24.943387854459719*x2 + x33;
    double x109 = x108 + x86;
    double x110 = x109*x88 - x109*x90 + x81;
    double x111 = -16.628925236306479*n3*x3 + 16.628925236306479*x2;
    double x112 = x1*x111;
    double x113 = x112 - x2*x25 + x38 + x77 + x79 + x92;
    double x114 = x113*x69;
    double x115 = -x113*x89 + x114*x96;
    double x116 = x1*(x49 + x7) + x59;
    double x117 = -2.0*x45 + x73*x96;
    double x118 = 1.0/n2;
    double x119 = -x13 + x15;
    double x120 = x1*(x5 + x55) + x118*x94 + x119 + x35;
    double x121 = 4.0*x3;
    double x122 = 1.0*x74;
    double x123 = -x120*x122 + x120*x70 - x121*x95 + x97*x99;
    double x124 = x1*(x55 + x85) + x108 + x119;
    double x125 = x115 + x81;
    double x126 = 2*x102 + x61 + x64*(x105 + x23);
    double x127 = 1.0/n3;
    double x128 = x1*(x57 + x62) + x112*x127 + x13 + x15 + 16.628925236306479*x2 - x3*(-x12 - x14) - x31;
    double x129 = -x113*x121 + x114*x99 - x122*x128 + x128*x70;
    double x130 = 3.0*x2;
    double x131 = 6.0*x3;
    double x132 = -n2*x48;
    double x133 = x104*x6;
    double x134 = x118*x52;
    double x135 = -x8;
    double x136 = n2*x54 + x135 + x58;
    double x137 = x118*x93 + x136;
    double x138 = 3.0*x74;
    double x139 = x124*x88 - x124*x90;

result[0] = 3.0*x37 - 6.0*x45 + x71*(x1*(x46 + x49) + x53*(x22 + x51) + x60 + x68 - x11/((n1)*(n1))) + x72*x73 - 3.0*x75 + x81;
result[1] = x100 + x71*(x5 + x67 + x83) + x91 + x98;
result[2] = x100 + x110 + x115 + x71*(x107 + x83);
result[3] = x117 + x123 + x71*(x116 + x67 + x84) + x91;
result[4] = -x109*x122 + x109*x70 + x117 - x122*x124 - x122*x87 + x124*x70 + x125 + x70*x87 + x71*(x106 + x116 + x4) + x98;
result[5] = x110 + x117 + x129 + x71*(x116 + x126);
result[6] = x120*x130 - x120*x138 - x131*x95 + x71*(x1*(x132 + x46) + x134*(x133 + x22) + x137 + x68 - x94/((n2)*(n2))) + x72*x97 + x81;
result[7] = x123 + x125 + x139 + x71*(x1*(x132 + x54) + x107 + x134*(x133 + x82) + x137);
result[8] = x129 + x139 + x71*(x1*(x132 + x7) + x126 + x136) + x81 + x98;
result[9] = -x113*x131 + x114*x72 + x128*x130 - x128*x138 + x71*(66.515700945225916*n3*x6 + 16.628925236306479*x1*x127*(n3*x23 + x22) + x1*(-99.773551417838874*n3*x47 + 99.773551417838874*x6) + x111*x127 + x135 - x23*x25 - 66.515700945225916*x3 + x56 - x112/((n3)*(n3))) + x81;
}
        
static double coder_dgdp(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = 1.0*(1.0*n1 + 1.0*n2 + 1.0*n3)*(n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P))/(n1 + n2 + n3);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x2 = n1 + n2 + n3;
    double x3 = 1.0/x2;
    double x4 = x1*x3;
    double x5 = (*endmember[1].dmu0dP)(T, P);
    double x6 = (*endmember[2].dmu0dP)(T, P);
    double x7 = n1*x0 + n2*x5 + n3*x6;
    double x8 = -1.0*x1*x7/((x2)*(x2)) + x3*x7;

result[0] = x0*x4 + x8;
result[1] = x4*x5 + x8;
result[2] = x4*x6 + x8;
}
        
static void coder_d3gdn2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3;
    double x2 = 1.0/x1;
    double x3 = x0*x2;
    double x4 = pow(x1, -2);
    double x5 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x6 = x4*x5;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = 2.0*n1*x0 + 2.0*n2*x8 + 2.0*n3*x9;
    double x11 = -x10*x4 + x10*x5/((x1)*(x1)*(x1));
    double x12 = x11 + 1.0*x3 - 1.0*x7;
    double x13 = 1.0*x2;
    double x14 = 1.0*x6;
    double x15 = x13*x8 - x14*x8;
    double x16 = x13*x9 - x14*x9;
    double x17 = 2.0*x2;
    double x18 = 2.0*x6;

result[0] = x11 + 2.0*x3 - 2.0*x7;
result[1] = x12 + x15;
result[2] = x12 + x16;
result[3] = x11 + x17*x8 - x18*x8;
result[4] = x11 + x15 + x16;
result[5] = x11 + x17*x9 - x18*x9;
}
        
static void coder_d4gdn3dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3;
    double x2 = pow(x1, -2);
    double x3 = x0*x2;
    double x4 = pow(x1, -3);
    double x5 = 6.0*x4;
    double x6 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = n1*x0 + n2*x8 + n3*x9;
    double x11 = x10*x5 - 6.0*x10*x6/((x1)*(x1)*(x1)*(x1));
    double x12 = 4.0*x4;
    double x13 = x11 + x12*x7 - 4.0*x3;
    double x14 = 2.0*x2;
    double x15 = 2.0*x4;
    double x16 = x6*x8;
    double x17 = -x14*x8 + x15*x16;
    double x18 = x6*x9;
    double x19 = -x14*x9 + x15*x18;
    double x20 = x11 + x15*x7 - 2.0*x3;
    double x21 = 4.0*x2;
    double x22 = x12*x16 - x21*x8;
    double x23 = x12*x18 - x21*x9;
    double x24 = 6.0*x2;

result[0] = x11 - 6.0*x3 + x5*x7;
result[1] = x13 + x17;
result[2] = x13 + x19;
result[3] = x20 + x22;
result[4] = x17 + x19 + x20;
result[5] = x20 + x23;
result[6] = x11 + x16*x5 - x24*x8;
result[7] = x11 + x19 + x22;
result[8] = x11 + x17 + x23;
result[9] = x11 + x18*x5 - x24*x9;
}
        
static double coder_d2gdt2(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    

result = 1.0*(n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d2mu0dT2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dT2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dT2)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d2mu0dTdP)(T, P) + 1.0*n2*(*endmember[1].d2mu0dTdP)(T, P) + 1.0*n3*(*endmember[2].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d2mu0dTdP)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dTdP)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dTdP)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P));
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d2mu0dP2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dP2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d3mu0dT3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT3)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P));
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d3mu0dT2dP)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT2dP)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT2dP)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d3mu0dTdP2)(T, P) + 1.0*n2*(*endmember[1].d3mu0dTdP2)(T, P) + 1.0*n3*(*endmember[2].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d3mu0dTdP2)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dTdP2)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dTdP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P));
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];


result[0] = 1.0*(*endmember[0].d3mu0dP3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dP3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dP3)(T, P);
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

