#include <math.h>


static double coder_g(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = 1.0*n4;
    double x1 = 1.0*n1;
    double x2 = 1.0*n3;
    double x3 = x1 + x2;
    double x4 = 1030400.0*n4;
    double x5 = 80.2881*T;
    double x6 = (1.0/3.0)*n2;
    double x7 = fmin(4, 0.19330141674964554*sqrt(-5.3525400000000003*T + 26.762700000000002));
    double x8 = n1 + n3;
    double x9 = 1.0/(n2 + n4 + x8);

result = 8.3144626181532395*T*(2.0*n2*log(n2*x9) + x0*log(n4*x9) + x2*log(n3*x9) + x3*log(x8*x9) + (x0 + x1)*log(x9*(n1 + n4))) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P) + ((T >= 5.0) ? (
   -x6*(x5 - 267.62700000000001)
)
: (
   x6*(133.8135*((x7)*(x7)*(x7)) + (x5 - 401.44049999999999)*(x7 - 1) - 133.8135)
)) + 0.015625*(n1*x4 + n2*x4 + 1536000.0*n3*n4 + 8.0*n4*(128800.0*n1 + 128800.0*n2 + 192000.0*n3))/(1.0*n2 + x0 + x3);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = 1.0/x1;
    double x3 = 1.0*n3;
    double x4 = -x2*x3;
    double x5 = 1.0*n4;
    double x6 = -x2*x5;
    double x7 = x4 + x6;
    double x8 = n2*x2;
    double x9 = -2.0*x8;
    double x10 = 1.0*n1;
    double x11 = x10 + x3;
    double x12 = pow(x1, -2);
    double x13 = x9 + 1.0*log(x0*x2) + x1*x11*(-x0*x12 + x2)/x0;
    double x14 = n1 + n4;
    double x15 = x10 + x5;
    double x16 = x1*x15*(-x12*x14 + x2)/x14 + 1.0*log(x14*x2);
    double x17 = 8.3144626181532395*T;
    double x18 = 1.0*n2 + x11 + x5;
    double x19 = 1.0/x18;
    double x20 = n4*x19;
    double x21 = 1030400.0*n4;
    double x22 = -0.015625*(n1*x21 + n2*x21 + 1536000.0*n3*n4 + 8.0*n4*(128800.0*n1 + 128800.0*n2 + 192000.0*n3))/((x18)*(x18));
    double x23 = 32200.0*x20 + x22;
    double x24 = -x11*x2;
    double x25 = -x15*x2;
    double x26 = fmin(4, 0.19330141674964554*sqrt(-5.3525400000000003*T + 26.762700000000002));
    double x27 = 1.0*x1;

result[0] = x17*(x13 + x16 + x7) + x23 + (*endmember[0].mu0)(T, P);
result[1] = x17*(2.0*x1*(-n2*x12 + x2) + x24 + x25 + x7 + 2.0*log(x8)) + x23 + ((T >= 5.0) ? (
   -26.762699999999999*T + 89.209000000000003
)
: (
   44.604500000000002*((x26)*(x26)*(x26)) + (1.0/3.0)*(80.2881*T - 401.44049999999999)*(x26 - 1) - 44.604500000000002
)) + (*endmember[1].mu0)(T, P);
result[2] = x17*(x13 + x25 + x27*(-n3*x12 + x2) + x6 + 1.0*log(n3*x2)) + 48000.0*x20 + x22 + (*endmember[2].mu0)(T, P);
result[3] = x17*(x16 + x24 + x27*(-n4*x12 + x2) + x4 + x9 + 1.0*log(n4*x2)) + 0.015625*x19*(2060800.0*n1 + 2060800.0*n2 + 3072000.0*n3) + x22 + (*endmember[3].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = pow(x1, -2);
    double x3 = 1.0*n3;
    double x4 = x2*x3;
    double x5 = 1.0*n4;
    double x6 = x2*x5;
    double x7 = x4 + x6;
    double x8 = 1.0/x0;
    double x9 = 1.0/x1;
    double x10 = -x0*x2 + x9;
    double x11 = x10*x8;
    double x12 = 2.0*x1;
    double x13 = -2*x2;
    double x14 = pow(x1, -3);
    double x15 = 2*x14;
    double x16 = 1.0*n1;
    double x17 = x16 + x3;
    double x18 = x1*x17;
    double x19 = x18*x8;
    double x20 = n2*x2;
    double x21 = 2.0*x20;
    double x22 = x11*x17 + x21;
    double x23 = x11*x12 + x19*(x0*x15 + x13) + x22 - x10*x18/((x0)*(x0));
    double x24 = x16 + x5;
    double x25 = n1 + n4;
    double x26 = 1.0/x25;
    double x27 = -x2*x25 + x9;
    double x28 = x26*x27;
    double x29 = x24*x28;
    double x30 = x1*x24;
    double x31 = x26*x30;
    double x32 = x12*x28 + x29 + x31*(x13 + x15*x25) - x27*x30/((x25)*(x25));
    double x33 = 8.3144626181532395*T;
    double x34 = 1.0*n2 + x17 + x5;
    double x35 = 1030400.0*n4;
    double x36 = 0.03125*(n1*x35 + n2*x35 + 1536000.0*n3*n4 + 8.0*n4*(128800.0*n1 + 128800.0*n2 + 192000.0*n3))/((x34)*(x34)*(x34));
    double x37 = pow(x34, -2);
    double x38 = n4*x37;
    double x39 = x36 - 64400.0*x38;
    double x40 = -x2;
    double x41 = -n1;
    double x42 = x19*(-x15*(-n3 + x41) + x40) + x22;
    double x43 = x29 + x31*(-x15*(-n4 + x41) + x40);
    double x44 = 2.0*x9;
    double x45 = -x44 + x7;
    double x46 = x36 - 80200.0*x38;
    double x47 = 1.0/x34;
    double x48 = x36 - 0.015625*x37*(2060800.0*n1 + 2060800.0*n2 + 3072000.0*n3);
    double x49 = -32200.0*x38 + 32200.0*x47 + x48;
    double x50 = -x16;
    double x51 = -x2*(-x3 + x50);
    double x52 = 4.0*n2*x14;
    double x53 = -x2*(-x5 + x50);
    double x54 = -x21;
    double x55 = x53 + x54;
    double x56 = -2.0*x2;
    double x57 = x1*(x52 + x56) + x45;
    double x58 = x2*x24;
    double x59 = 1.0*x9;
    double x60 = 2.0*x14;
    double x61 = n3*x60;
    double x62 = 1.0*x1;
    double x63 = -x4 + x6;

result[0] = x33*(x23 + x32 + x7) + x39;
result[1] = x33*(x42 + x43 + x7 - 4.0*x9) + x39;
result[2] = x33*(x23 + x43 + x45) + x46;
result[3] = x33*(x32 + x42 + x45) + x49;
result[4] = x33*(x1*(-4.0*x2 + x52) + x44 + x51 + x55 + x7 + x12*(-x20 + x9)/n2) + x39;
result[5] = x33*(x17*x2 + x55 + x57) + x46;
result[6] = x33*(x51 + x54 + x57 + x58) + x49;
result[7] = x33*(x1*(x56 + x61) + x23 + x53 + x59 + x63 + x62*(-n3*x2 + x9)/n3) + x36 - 96000.0*x38;
result[8] = x33*(x1*(-1.0*x2 + x61) + x42 + x58 + x63 - 3.0*x9) - 48000.0*x38 + 48000.0*x47 + x48;
result[9] = x33*(x1*(n4*x60 + x56) + x21 + x32 + x4 + x51 + x59 - x6 + x62*(-n4*x2 + x9)/n4) - 1.0*x37*(32200.0*n1 + 32200.0*n2 + 48000.0*n3) + x48;
}
        
static void coder_d3gdn3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n4;
    double x1 = 1.0/x0;
    double x2 = n1 + n3;
    double x3 = n2 + n4 + x2;
    double x4 = 1.0/x3;
    double x5 = pow(x3, -2);
    double x6 = -x0*x5 + x4;
    double x7 = x1*x6;
    double x8 = 1.0*n1;
    double x9 = 1.0*n4;
    double x10 = x8 + x9;
    double x11 = pow(x0, -2);
    double x12 = x11*x6;
    double x13 = x10*x12;
    double x14 = -2*x5;
    double x15 = pow(x3, -3);
    double x16 = 2*x15;
    double x17 = x0*x16 + x14;
    double x18 = x1*x10;
    double x19 = x17*x18;
    double x20 = 3.0*x3;
    double x21 = 6*x15;
    double x22 = pow(x3, -4);
    double x23 = 6*x22;
    double x24 = x18*x3;
    double x25 = x10*x3;
    double x26 = x11*x25;
    double x27 = x1*x17*x20 - x12*x20 - 2*x13 - 2*x17*x26 + 2*x19 + x24*(-x0*x23 + x21) + 3.0*x7 + 2*x25*x6/((x0)*(x0)*(x0));
    double x28 = 4.0*x15;
    double x29 = -n2*x28;
    double x30 = 2.0*x15;
    double x31 = -n3*x30;
    double x32 = -n4*x30;
    double x33 = x31 + x32;
    double x34 = x29 + x33;
    double x35 = 1.0/x2;
    double x36 = -x2*x5 + x4;
    double x37 = x35*x36;
    double x38 = 1.0*n3;
    double x39 = x38 + x8;
    double x40 = pow(x2, -2);
    double x41 = x36*x40;
    double x42 = x39*x41;
    double x43 = x14 + x16*x2;
    double x44 = x35*x39;
    double x45 = x43*x44;
    double x46 = x3*x44;
    double x47 = x3*x39;
    double x48 = x40*x47;
    double x49 = x20*x35*x43 - x20*x41 + 3.0*x37 - 2*x42 - 2*x43*x48 + 2*x45 + x46*(-x2*x23 + x21) + 2*x36*x47/((x2)*(x2)*(x2));
    double x50 = x34 + x49;
    double x51 = 8.3144626181532395*T;
    double x52 = 1.0*n2 + x39 + x9;
    double x53 = pow(x52, -3);
    double x54 = n4*x53;
    double x55 = 1030400.0*n4;
    double x56 = -0.09375*(n1*x55 + n2*x55 + 1536000.0*n3*n4 + 8.0*n4*(128800.0*n1 + 128800.0*n2 + 192000.0*n3))/((x52)*(x52)*(x52)*(x52));
    double x57 = 193200.0*x54 + x56;
    double x58 = x33 + 2.0*x5;
    double x59 = 2.0*x3;
    double x60 = -x5;
    double x61 = -n1;
    double x62 = -n3 + x61;
    double x63 = -x16*x62 + x60;
    double x64 = x35*x63;
    double x65 = 4*x15;
    double x66 = 2*n1;
    double x67 = 2*n3;
    double x68 = 3*x22;
    double x69 = -x68*(x66 + x67);
    double x70 = x44*x63;
    double x71 = -x42 + x45 - x48*x63 + x70;
    double x72 = 2.0*x37 + x46*(x65 + x69) + x59*x64 + x71;
    double x73 = -n4 + x61;
    double x74 = -x16*x73 + x60;
    double x75 = x1*x74;
    double x76 = 2*n4;
    double x77 = -x68*(x66 + x76);
    double x78 = x18*x74;
    double x79 = -x13 + x19 - x26*x74 + x78;
    double x80 = x24*(x65 + x77) + x59*x75 + 2.0*x7 + x79;
    double x81 = 1.0*x5;
    double x82 = 224800.0*x54 + x56;
    double x83 = x34 + x72;
    double x84 = pow(x52, -2);
    double x85 = 0.03125*x53*(2060800.0*n1 + 2060800.0*n2 + 3072000.0*n3) + x56;
    double x86 = -64400.0*x84 + x85;
    double x87 = 128800.0*x54 + x86;
    double x88 = x46*(x16 + x69) + 2*x70;
    double x89 = x24*(x16 + x77) + x34 + 2*x78;
    double x90 = 5.0*x5;
    double x91 = 1.0*x3;
    double x92 = x34 + x88;
    double x93 = x24*(x23*x73 + x65) + 1.0*x7 + x75*x91 + x79;
    double x94 = 3.0*x5;
    double x95 = 256400.0*x54 + x56;
    double x96 = 160400.0*x54 - 80200.0*x84 + x85;
    double x97 = x53*(-32200.0*n1 - 32200.0*n2 - 48000.0*n3);
    double x98 = -2.0*x97;
    double x99 = 64400.0*x54 + x86 + x98;
    double x100 = -12.0*n2*x22;
    double x101 = n2*x16;
    double x102 = 2.0/n2;
    double x103 = x102*x3;
    double x104 = -n2*x5 + x4;
    double x105 = -x16*x39;
    double x106 = -x10*x16;
    double x107 = 8.0*x15;
    double x108 = n2*x107;
    double x109 = x106 + x108;
    double x110 = x105 + x109;
    double x111 = x102*x104 + x33;
    double x112 = -x8;
    double x113 = 4.0*x5;
    double x114 = -x113;
    double x115 = x103*(x101 + x60) + x111 + x114 + x3*(x100 + x107);
    double x116 = x16*(x112 - x9);
    double x117 = x105 + x108 + x116;
    double x118 = x3*(x100 + x28) + x58;
    double x119 = x51*(x110 + x118);
    double x120 = -n3*x5 + x4;
    double x121 = 1.0/n3;
    double x122 = x120*x121;
    double x123 = 6.0*x15;
    double x124 = 6.0*x22;
    double x125 = -n3*x124;
    double x126 = x15*x67;
    double x127 = x121*x3;
    double x128 = n3*x28 + x29 + x32;
    double x129 = x106 + x128;
    double x130 = -96000.0*x84 + x85;
    double x131 = -n4*x5 + x4;
    double x132 = 1.0/n4;

result[0] = x51*(x27 + x50) + x57;
result[1] = x51*(x29 + x58 + x72 + x80) + x57;
result[2] = x51*(x50 + x80 + x81) + x82;
result[3] = x51*(x27 + x81 + x83) + x87;
result[4] = x51*(6.0*x5 + x88 + x89) + x57;
result[5] = x51*(1.0*x37 + x46*(x23*x62 + x65) + x64*x91 + x71 + x89 + x90) + x82;
result[6] = x51*(x90 + x92 + x93) + x87;
result[7] = x51*(x49 + x89 + x94) + x95;
result[8] = x51*(x83 + x93 + x94) + x96;
result[9] = x51*(x27 + x92 + x94) + x99;
result[10] = x51*(x103*(x101 + x14) + x110 + x111 + x3*(x100 + 12.0*x15) - 8.0*x5 - x104*x59/((n2)*(n2))) + x57;
result[11] = x51*(x109 + x115 + x16*(x112 - x38)) + x82;
result[12] = x51*(x115 + x117) + x87;
result[13] = x119 + x95;
result[14] = x51*(x117 + x118) + x96;
result[15] = x119 + x99;
result[16] = x51*(x114 + x122 + x127*(x126 + x14) + x129 + x3*(x123 + x125) + x49 - x120*x91/((n3)*(n3))) + 288000.0*x54 + x56;
result[17] = x130 + x51*(x116 + x122 + x127*(x126 + x60) + x128 + x3*(x125 + x28) + x72 - x81) + 192000.0*x54;
result[18] = x130 + x51*(x113 + x129 + x3*(x125 + x30) + x88) + 96000.0*x54 + x98;
result[19] = x51*(n4*x28 + x105 + x114 + x131*x132 + x132*x3*(x14 + x15*x76) + x27 + x29 + x3*(-n4*x124 + x123) + x31 - x131*x91/((n4)*(n4))) + x85 - 4.0*x97;
}
        
static double coder_dgdt(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = n1 + n3;
    double x1 = 1.0/(n2 + n4 + x0);
    double x2 = 1.0*n1;
    double x3 = sqrt(-5.3525400000000003*T + 26.762700000000002);
    double x4 = 0.19330141674964554*x3;
    double x5 = fmin(4, x4);
    double x6 = (-x4 + 4 >= 0. ? 1. : 0.)/x3;

result = n1*(*endmember[0].dmu0dT)(T, P) + 16.628925236306479*n2*log(n2*x1) + n2*(*endmember[1].dmu0dT)(T, P) + 8.3144626181532395*n3*log(n3*x1) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*log(n4*x1) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626181532395*(1.0*n3 + x2)*log(x0*x1) + 8.3144626181532395*(1.0*n4 + x2)*log(x1*(n1 + n4)) + ((T >= 5.0) ? (
   -26.762699999999999*n2
)
: (
   (1.0/3.0)*n2*(-207.67592227217145*((x5)*(x5))*x6 + 80.2881*x5 - 0.51732678260457388*x6*(80.2881*T - 401.44049999999999) - 80.2881)
));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = 1.0/x1;
    double x3 = 8.3144626181532395*n3;
    double x4 = -x2*x3;
    double x5 = 8.3144626181532395*n4;
    double x6 = -x2*x5;
    double x7 = x4 + x6;
    double x8 = n2*x2;
    double x9 = -16.628925236306479*x8;
    double x10 = 8.3144626181532395*n1;
    double x11 = x10 + x3;
    double x12 = pow(x1, -2);
    double x13 = x9 + 8.3144626181532395*log(x0*x2) + x1*x11*(-x0*x12 + x2)/x0;
    double x14 = n1 + n4;
    double x15 = x10 + x5;
    double x16 = x1*x15*(-x12*x14 + x2)/x14 + 8.3144626181532395*log(x14*x2);
    double x17 = -x11*x2;
    double x18 = -x15*x2;
    double x19 = sqrt(-5.3525400000000003*T + 26.762700000000002);
    double x20 = 0.19330141674964554*x19;
    double x21 = fmin(4, x20);
    double x22 = (-x20 + 4 >= 0. ? 1. : 0.)/x19;
    double x23 = 8.3144626181532395*x1;

result[0] = x13 + x16 + x7 + (*endmember[0].dmu0dT)(T, P);
result[1] = 16.628925236306479*x1*(-n2*x12 + x2) + x17 + x18 + x7 + ((T >= 5.0) ? (
   -26.762699999999999
)
: (
   -69.225307424057149*((x21)*(x21))*x22 + 26.762699999999999*x21 - 0.17244226086819128*x22*(80.2881*T - 401.44049999999999) - 26.762699999999999
)) + 16.628925236306479*log(x8) + (*endmember[1].dmu0dT)(T, P);
result[2] = x13 + x18 + x23*(-n3*x12 + x2) + x6 + 8.3144626181532395*log(n3*x2) + (*endmember[2].dmu0dT)(T, P);
result[3] = x16 + x17 + x23*(-n4*x12 + x2) + x4 + x9 + 8.3144626181532395*log(n4*x2) + (*endmember[3].dmu0dT)(T, P);
}
        
static void coder_d3gdn2dt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = pow(x1, -2);
    double x3 = 8.3144626181532395*n3;
    double x4 = x2*x3;
    double x5 = 8.3144626181532395*n4;
    double x6 = x2*x5;
    double x7 = x4 + x6;
    double x8 = 1.0/x0;
    double x9 = 1.0/x1;
    double x10 = -x0*x2 + x9;
    double x11 = x10*x8;
    double x12 = 16.628925236306479*x1;
    double x13 = -2*x2;
    double x14 = pow(x1, -3);
    double x15 = 2*x14;
    double x16 = 8.3144626181532395*n1;
    double x17 = x16 + x3;
    double x18 = x1*x17;
    double x19 = x18*x8;
    double x20 = n2*x2;
    double x21 = 16.628925236306479*x20;
    double x22 = x11*x17 + x21;
    double x23 = x11*x12 + x19*(x0*x15 + x13) + x22 - x10*x18/((x0)*(x0));
    double x24 = x16 + x5;
    double x25 = n1 + n4;
    double x26 = 1.0/x25;
    double x27 = -x2*x25 + x9;
    double x28 = x26*x27;
    double x29 = x24*x28;
    double x30 = x1*x24;
    double x31 = x26*x30;
    double x32 = x12*x28 + x29 + x31*(x13 + x15*x25) - x27*x30/((x25)*(x25));
    double x33 = -x2;
    double x34 = -n1;
    double x35 = x19*(-x15*(-n3 + x34) + x33) + x22;
    double x36 = x29 + x31*(-x15*(-n4 + x34) + x33);
    double x37 = 16.628925236306479*x9;
    double x38 = -x37 + x7;
    double x39 = -x16;
    double x40 = -x2*(-x3 + x39);
    double x41 = 33.257850472612958*n2*x14;
    double x42 = -x2*(x39 - x5);
    double x43 = -x21;
    double x44 = x42 + x43;
    double x45 = -16.628925236306479*x2;
    double x46 = x1*(x41 + x45) + x38;
    double x47 = x2*x24;
    double x48 = 8.3144626181532395*x9;
    double x49 = 16.628925236306479*x14;
    double x50 = n3*x49;
    double x51 = 8.3144626181532395*x1;
    double x52 = -x4 + x6;

result[0] = x23 + x32 + x7;
result[1] = x35 + x36 + x7 - 33.257850472612958*x9;
result[2] = x23 + x36 + x38;
result[3] = x32 + x35 + x38;
result[4] = x1*(-33.257850472612958*x2 + x41) + x37 + x40 + x44 + x7 + x12*(-x20 + x9)/n2;
result[5] = x17*x2 + x44 + x46;
result[6] = x40 + x43 + x46 + x47;
result[7] = x1*(x45 + x50) + x23 + x42 + x48 + x52 + x51*(-n3*x2 + x9)/n3;
result[8] = x1*(-8.3144626181532395*x2 + x50) + x35 + x47 + x52 - 24.943387854459719*x9;
result[9] = x1*(n4*x49 + x45) + x21 + x32 + x4 + x40 + x48 - x6 + x51*(-n4*x2 + x9)/n4;
}
        
static void coder_d4gdn3dt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n4;
    double x1 = 1.0/x0;
    double x2 = n1 + n3;
    double x3 = n2 + n4 + x2;
    double x4 = 1.0/x3;
    double x5 = pow(x3, -2);
    double x6 = -x0*x5 + x4;
    double x7 = x1*x6;
    double x8 = 8.3144626181532395*n1;
    double x9 = 8.3144626181532395*n4;
    double x10 = x8 + x9;
    double x11 = pow(x0, -2);
    double x12 = x11*x6;
    double x13 = x10*x12;
    double x14 = -2*x5;
    double x15 = pow(x3, -3);
    double x16 = 2*x15;
    double x17 = x0*x16 + x14;
    double x18 = x1*x10;
    double x19 = x17*x18;
    double x20 = 24.943387854459719*x3;
    double x21 = 6*x15;
    double x22 = pow(x3, -4);
    double x23 = 6*x22;
    double x24 = x18*x3;
    double x25 = x10*x3;
    double x26 = x11*x25;
    double x27 = x1*x17*x20 - x12*x20 - 2*x13 - 2*x17*x26 + 2*x19 + x24*(-x0*x23 + x21) + 24.943387854459719*x7 + 2*x25*x6/((x0)*(x0)*(x0));
    double x28 = 33.257850472612958*x15;
    double x29 = -n2*x28;
    double x30 = 16.628925236306479*x15;
    double x31 = -n3*x30;
    double x32 = -n4*x30;
    double x33 = x31 + x32;
    double x34 = x29 + x33;
    double x35 = 1.0/x2;
    double x36 = -x2*x5 + x4;
    double x37 = x35*x36;
    double x38 = 8.3144626181532395*n3;
    double x39 = x38 + x8;
    double x40 = pow(x2, -2);
    double x41 = x36*x40;
    double x42 = x39*x41;
    double x43 = x14 + x16*x2;
    double x44 = x35*x39;
    double x45 = x43*x44;
    double x46 = x3*x44;
    double x47 = x3*x39;
    double x48 = x40*x47;
    double x49 = x20*x35*x43 - x20*x41 + 24.943387854459719*x37 - 2*x42 - 2*x43*x48 + 2*x45 + x46*(-x2*x23 + x21) + 2*x36*x47/((x2)*(x2)*(x2));
    double x50 = x34 + x49;
    double x51 = x33 + 16.628925236306479*x5;
    double x52 = 16.628925236306479*x3;
    double x53 = -x5;
    double x54 = -n1;
    double x55 = -n3 + x54;
    double x56 = -x16*x55 + x53;
    double x57 = x35*x56;
    double x58 = 4*x15;
    double x59 = 2*n1;
    double x60 = 2*n3;
    double x61 = 3*x22;
    double x62 = -x61*(x59 + x60);
    double x63 = x44*x56;
    double x64 = -x42 + x45 - x48*x56 + x63;
    double x65 = 16.628925236306479*x37 + x46*(x58 + x62) + x52*x57 + x64;
    double x66 = -n4 + x54;
    double x67 = -x16*x66 + x53;
    double x68 = x1*x67;
    double x69 = 2*n4;
    double x70 = -x61*(x59 + x69);
    double x71 = x18*x67;
    double x72 = -x13 + x19 - x26*x67 + x71;
    double x73 = x24*(x58 + x70) + x52*x68 + 16.628925236306479*x7 + x72;
    double x74 = 8.3144626181532395*x5;
    double x75 = x34 + x65;
    double x76 = x46*(x16 + x62) + 2*x63;
    double x77 = x24*(x16 + x70) + x34 + 2*x71;
    double x78 = 41.572313090766201*x5;
    double x79 = 8.3144626181532395*x3;
    double x80 = x34 + x76;
    double x81 = x24*(x23*x66 + x58) + x68*x79 + 8.3144626181532395*x7 + x72;
    double x82 = 24.943387854459719*x5;
    double x83 = -99.773551417838874*n2*x22;
    double x84 = n2*x16;
    double x85 = 16.628925236306479/n2;
    double x86 = x3*x85;
    double x87 = -n2*x5 + x4;
    double x88 = -x16*x39;
    double x89 = -x10*x16;
    double x90 = 66.515700945225916*x15;
    double x91 = n2*x90;
    double x92 = x89 + x91;
    double x93 = x88 + x92;
    double x94 = x33 + x85*x87;
    double x95 = -x8;
    double x96 = 33.257850472612958*x5;
    double x97 = -x96;
    double x98 = x3*(x83 + x90) + x86*(x53 + x84) + x94 + x97;
    double x99 = x16*(-x9 + x95);
    double x100 = x88 + x91 + x99;
    double x101 = x3*(x28 + x83) + x51;
    double x102 = x101 + x93;
    double x103 = -n3*x5 + x4;
    double x104 = 8.3144626181532395/n3;
    double x105 = x103*x104;
    double x106 = 49.886775708919437*x15;
    double x107 = 49.886775708919437*x22;
    double x108 = -n3*x107;
    double x109 = x15*x60;
    double x110 = x104*x3;
    double x111 = n3*x28 + x29 + x32;
    double x112 = x111 + x89;
    double x113 = -n4*x5 + x4;
    double x114 = 8.3144626181532395/n4;

result[0] = x27 + x50;
result[1] = x29 + x51 + x65 + x73;
result[2] = x50 + x73 + x74;
result[3] = x27 + x74 + x75;
result[4] = 49.886775708919437*x5 + x76 + x77;
result[5] = 8.3144626181532395*x37 + x46*(x23*x55 + x58) + x57*x79 + x64 + x77 + x78;
result[6] = x78 + x80 + x81;
result[7] = x49 + x77 + x82;
result[8] = x75 + x81 + x82;
result[9] = x27 + x80 + x82;
result[10] = x3*(99.773551417838874*x15 + x83) - 66.515700945225916*x5 + x86*(x14 + x84) + x93 + x94 - x52*x87/((n2)*(n2));
result[11] = x16*(-x38 + x95) + x92 + x98;
result[12] = x100 + x98;
result[13] = x102;
result[14] = x100 + x101;
result[15] = x102;
result[16] = x105 + x110*(x109 + x14) + x112 + x3*(x106 + x108) + x49 + x97 - x103*x79/((n3)*(n3));
result[17] = x105 + x110*(x109 + x53) + x111 + x3*(x108 + x28) + x65 - x74 + x99;
result[18] = x112 + x3*(x108 + x30) + x76 + x96;
result[19] = n4*x28 + x113*x114 + x114*x3*(x14 + x15*x69) + x27 + x29 + x3*(-n4*x107 + x106) + x31 + x88 + x97 - x113*x79/((n4)*(n4));
}
        
static double coder_dgdp(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P) + n4*(*endmember[3].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].dmu0dP)(T, P);
result[1] = (*endmember[1].dmu0dP)(T, P);
result[2] = (*endmember[2].dmu0dP)(T, P);
result[3] = (*endmember[3].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d4gdn3dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d2gdt2(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 26.762700000000002;
    double x2 = sqrt(x1);
    double x3 = 0.19330141674964554*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 80.2881*T - 401.44049999999999;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 26.762700000000002);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P) + ((T >= 5.0) ? (
   0
)
: (
   -1.0/3.0*n2*(555.79684049934428*x10*x6 - 107.43631669350002*x10*x8 + 214.87263338700004*((x4)*(x4))*x7*x9 + 1.384506148481143*x5*x6 - 0.26762700000000006*x5*x8 + 83.070368908868573*x4/x2)
));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 26.762700000000002;
    double x2 = sqrt(x1);
    double x3 = 0.19330141674964554*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 80.2881*T - 401.44049999999999;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 26.762700000000002);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = ((T >= 5.0) ? (
   0
)
: (
   -185.26561349978141*x10*x6 + 35.812105564500001*x10*x8 - 71.624211129000003*((x4)*(x4))*x7*x9 - 0.46150204949371432*x5*x6 + 0.089209000000000011*x5*x8 - 27.690122969622855*x4/x2
)) + (*endmember[1].d2mu0dT2)(T, P);
result[2] = (*endmember[2].d2mu0dT2)(T, P);
result[3] = (*endmember[3].d2mu0dT2)(T, P);
}
        
static void coder_d4gdn2dt2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d5gdn3dt2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d2gdtdp(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].d2mu0dTdP)(T, P) + n2*(*endmember[1].d2mu0dTdP)(T, P) + n3*(*endmember[2].d2mu0dTdP)(T, P) + n4*(*endmember[3].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].d2mu0dTdP)(T, P);
result[1] = (*endmember[1].d2mu0dTdP)(T, P);
result[2] = (*endmember[2].d2mu0dTdP)(T, P);
result[3] = (*endmember[3].d2mu0dTdP)(T, P);
}
        
static void coder_d4gdn2dtdp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d5gdn3dtdp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d2gdp2(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P) + n4*(*endmember[3].d2mu0dP2)(T, P);
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].d2mu0dP2)(T, P);
result[1] = (*endmember[1].d2mu0dP2)(T, P);
result[2] = (*endmember[2].d2mu0dP2)(T, P);
result[3] = (*endmember[3].d2mu0dP2)(T, P);
}
        
static void coder_d4gdn2dp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d5gdn3dp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d3gdt3(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 26.762700000000002;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.19330141674964554*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = x2*x4;
    double x6 = x0 - 26.762700000000002;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 80.2881*T - 401.44049999999999;
    double x10 = x4/pow(x1, 5.0/2.0);
    double x11 = pow(x6, -2);
    double x12 = x11*x8;
    double x13 = x2*0;
    double x14 = fmin(4, x3);
    double x15 = ((x14)*(x14));

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P) + ((T >= 5.0) ? (
   0
)
: (
   -1.0/3.0*n2*(4462.3872309695407*x10*x15 + 11.115936809986886*x10*x9 - 1725.1715476638797*x11*x14*((x4)*(x4)) + 862.58577383193983*x12*x15 + 2.1487263338700004*x12*x9 - 55.579684049934436*x13*x15 - 0.13845061484811433*x13*x9 - 333.47810429960663*x14*x5*x8 + 111.15936809986887*x2*((x4)*(x4)*(x4)) + 333.47810429960657*x5 - 64.461790016100011*x8/x6)
));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 26.762700000000002;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.19330141674964554*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = x2*x4;
    double x6 = x0 - 26.762700000000002;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 80.2881*T - 401.44049999999999;
    double x10 = x4/pow(x1, 5.0/2.0);
    double x11 = pow(x6, -2);
    double x12 = x11*x8;
    double x13 = x2*0;
    double x14 = fmin(4, x3);
    double x15 = ((x14)*(x14));

result[0] = (*endmember[0].d3mu0dT3)(T, P);
result[1] = ((T >= 5.0) ? (
   0
)
: (
   -1487.4624103231802*x10*x15 - 3.7053122699956287*x10*x9 + 575.05718255462648*x11*x14*((x4)*(x4)) - 287.52859127731324*x12*x15 - 0.71624211129000015*x12*x9 + 18.526561349978145*x13*x15 + 0.046150204949371443*x13*x9 + 111.15936809986887*x14*x5*x8 - 37.05312269995629*x2*((x4)*(x4)*(x4)) - 111.15936809986886*x5 + 21.487263338700004*x8/x6
)) + (*endmember[1].d3mu0dT3)(T, P);
result[2] = (*endmember[2].d3mu0dT3)(T, P);
result[3] = (*endmember[3].d3mu0dT3)(T, P);
}
        
static void coder_d5gdn2dt3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d6gdn3dt3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d3gdt2dp(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P) + n4*(*endmember[3].d3mu0dT2dP)(T, P);
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].d3mu0dT2dP)(T, P);
result[1] = (*endmember[1].d3mu0dT2dP)(T, P);
result[2] = (*endmember[2].d3mu0dT2dP)(T, P);
result[3] = (*endmember[3].d3mu0dT2dP)(T, P);
}
        
static void coder_d5gdn2dt2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d6gdn3dt2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d3gdtdp2(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].d3mu0dTdP2)(T, P) + n2*(*endmember[1].d3mu0dTdP2)(T, P) + n3*(*endmember[2].d3mu0dTdP2)(T, P) + n4*(*endmember[3].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].d3mu0dTdP2)(T, P);
result[1] = (*endmember[1].d3mu0dTdP2)(T, P);
result[2] = (*endmember[2].d3mu0dTdP2)(T, P);
result[3] = (*endmember[3].d3mu0dTdP2)(T, P);
}
        
static void coder_d5gdn2dtdp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d6gdn3dtdp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_d3gdp3(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P) + n4*(*endmember[3].d3mu0dP3)(T, P);
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = (*endmember[0].d3mu0dP3)(T, P);
result[1] = (*endmember[1].d3mu0dP3)(T, P);
result[2] = (*endmember[2].d3mu0dP3)(T, P);
result[3] = (*endmember[3].d3mu0dP3)(T, P);
}
        
static void coder_d5gdn2dp3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
        
static void coder_d6gdn3dp3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


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
result[10] = 0;
result[11] = 0;
result[12] = 0;
result[13] = 0;
result[14] = 0;
result[15] = 0;
result[16] = 0;
result[17] = 0;
result[18] = 0;
result[19] = 0;
}
        
static double coder_s(double T, double P, double n[4]) {
    double result = -coder_dgdt(T, P, n);
    return result;
}

static double coder_v(double T, double P, double n[4]) {
    double result = coder_dgdp(T, P, n);
    return result;
}

static double coder_cv(double T, double P, double n[4]) {
    double result = -T*coder_d2gdt2(T, P, n);
    double dvdt = coder_d2gdtdp(T, P, n);
    double dvdp = coder_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P, double n[4]) {
    double result = -T*coder_d2gdt2(T, P, n);
    return result;
}

static double coder_dcpdt(double T, double P, double n[4]) {
    double result = -T*coder_d3gdt3(T, P, n) - coder_d2gdt2(T, P, n);
    return result;
}

static double coder_alpha(double T, double P, double n[4]) {
    double result = coder_d2gdtdp(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_beta(double T, double P, double n[4]) {
    double result = -coder_d2gdp2(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_K(double T, double P, double n[4]) {
    double result = -coder_dgdp(T, P, n)/coder_d2gdp2(T, P, n);
    return result;
}

static double coder_Kp(double T, double P, double n[4]) {
    double result = coder_dgdp(T, P, n);
    result *= coder_d3gdp3(T, P, n);
    result /= pow(coder_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}

