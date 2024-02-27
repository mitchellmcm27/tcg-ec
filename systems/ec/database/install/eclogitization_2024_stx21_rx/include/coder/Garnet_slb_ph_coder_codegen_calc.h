#include <math.h>


static double coder_g(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = 120.43215000000001*T;
    double x1 = (1.0/3.0)*n2;
    double x2 = fmin(4, 0.15782994586457105*sqrt(-5.3525400000000003*T + 40.14405));
    double x3 = 1.0*n4;
    double x4 = 1.0*n5;
    double x5 = x3 + x4;
    double x6 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x7 = 181600.0*n4 + 183200.0*n5;
    double x8 = 0.82399999999999995*P + 168800.0;
    double x9 = n4 + n5;
    double x10 = n1 + n2 + n3;
    double x11 = 1.0/(x10 + x9);
    double x12 = log(n5*x11);

result = 8.3144626181532395*T*(3.0*n2*log(n2*x11) + 3.0*n3*log(n3*x11) + 0.99999999900000003*n5*(x12 - 1.0986122896681101) + 1.9999999980000001*n5*(x12 - 0.40546510910816402) + x3*log(n4*x11) + x5*log(x11*x9) + x6*log(x10*x11) + (3.0*n1 + 3.0*n4)*log(x11*(n1 + n4)) + (x4 + x6)*log(x11*(n5 + x10))) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P) + n5*(*endmember[4].mu0)(T, P) + ((T >= 7.5) ? (
   -x1*(x0 - 602.16075000000012)
)
: (
   x1*(301.080375*((x2)*(x2)*(x2)) + (x0 - 903.24112500000001)*(x2 - 1) - 301.080375)
)) + 0.0009765625*(64.0*n1*(n3*x8 + x7) + 64.0*n2*(168800.0*n3 + x7) + 64.0*n3*(n1*x8 + 168800.0*n2 + 488000.0*n4 + 485600.0*n5) + 64.0*n4*(181600.0*n1 + 181600.0*n2 + 488000.0*n3 + 568000.0*n5) + 64.0*n5*(183200.0*n1 + 183200.0*n2 + 485600.0*n3 + 568000.0*n4))/(x5 + x6);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 0.82399999999999995*P + 168800.0;
    double x1 = 64.0*n3;
    double x2 = 0.10299999999999999*P + 21100.0;
    double x3 = n3*x2;
    double x4 = 23244800.0*n4 + 23449600.0*n5;
    double x5 = 1.0*n4;
    double x6 = 1.0*n5;
    double x7 = x5 + x6;
    double x8 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x9 = x7 + x8;
    double x10 = 0.0009765625/x9;
    double x11 = 181600.0*n4 + 183200.0*n5;
    double x12 = 64.0*n1;
    double x13 = n1*x2;
    double x14 = -0.0009765625*(64.0*n2*(168800.0*n3 + x11) + 64.0*n4*(181600.0*n1 + 181600.0*n2 + 488000.0*n3 + 568000.0*n5) + 64.0*n5*(183200.0*n1 + 183200.0*n2 + 485600.0*n3 + 568000.0*n4) + x1*(168800.0*n2 + 488000.0*n4 + 485600.0*n5 + 8.0*x13) + x12*(x11 + 8.0*x3))/((x9)*(x9));
    double x15 = n1 + n4;
    double x16 = n2 + n3;
    double x17 = n5 + x15 + x16;
    double x18 = 1.0/x17;
    double x19 = n2*x18;
    double x20 = -3.0*x19;
    double x21 = n3*x18;
    double x22 = -3.0*x21;
    double x23 = n5*x18;
    double x24 = -2.9999999970000002*x23;
    double x25 = 3.0*n1 + 3.0*n4;
    double x26 = pow(x17, -2);
    double x27 = x20 + x22 + x24 + 3.0*log(x15*x18) + x17*x25*(-x15*x26 + x18)/x15;
    double x28 = n1 + x16;
    double x29 = n5 + x28;
    double x30 = x6 + x8;
    double x31 = x17*x30*(x18 - x26*x29)/x29 - x18*x5 + 1.0*log(x18*x29);
    double x32 = x17*x8*(x18 - x26*x28)/x28 - x18*x7 + x31 + 1.0*log(x18*x28);
    double x33 = 8.3144626181532395*T;
    double x34 = fmin(4, 0.15782994586457105*sqrt(-5.3525400000000003*T + 40.14405));
    double x35 = 3.0*x17;
    double x36 = -x18*x25;
    double x37 = x22 + x36;
    double x38 = x24 + x32;
    double x39 = n4 + n5;
    double x40 = x17*x7*(x18 - x26*x39)/x39 - x18*x8 + 1.0*log(x18*x39);

result[0] = x10*(x0*x1 + 512.0*x3 + x4) + x14 + x33*(x27 + x32) + (*endmember[0].mu0)(T, P);
result[1] = x10*(21606400.0*n3 + x4) + x14 + x33*(x35*(-n2*x26 + x18) + x37 + x38 + 3.0*log(x19)) + ((T >= 7.5) ? (
   -40.14405*T + 200.72025000000002
)
: (
   100.360125*((x34)*(x34)*(x34)) + (1.0/3.0)*(120.43215000000001*T - 903.24112500000001)*(x34 - 1) - 100.360125
)) + (*endmember[1].mu0)(T, P);
result[2] = x10*(21606400.0*n2 + 62464000.0*n4 + 62156800.0*n5 + x0*x12 + 512.0*x13) + x14 + x33*(x20 + x35*(-n3*x26 + x18) + x36 + x38 + 3.0*log(x21)) + (*endmember[2].mu0)(T, P);
result[3] = x10*(23244800.0*n1 + 23244800.0*n2 + 62464000.0*n3 + 72704000.0*n5) + x14 + x33*(1.0*x17*(-n4*x26 + x18) - x18*x30 + x27 + x40 + 1.0*log(n4*x18)) + (*endmember[3].mu0)(T, P);
result[4] = x10*(23449600.0*n1 + 23449600.0*n2 + 62156800.0*n3 + 72704000.0*n4) + x14 + x33*(2.9999999970000002*x17*(-n5*x26 + x18) + x20 + x31 + x37 + x40 + 2.9999999970000002*log(x23) - 1.9095425059748956) + (*endmember[4].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = n2 + n3;
    double x2 = n5 + x0 + x1;
    double x3 = pow(x2, -2);
    double x4 = 3.0*x3;
    double x5 = n3*x4;
    double x6 = n5*x3;
    double x7 = 2.9999999970000002*x6;
    double x8 = x5 + x7;
    double x9 = 1.0*n1;
    double x10 = 1.0*n2;
    double x11 = 1.0*n3;
    double x12 = x10 + x11 + x9;
    double x13 = n1 + x1;
    double x14 = 1.0/x13;
    double x15 = 1.0/x2;
    double x16 = -x13*x3 + x15;
    double x17 = x14*x16;
    double x18 = 1.0*n4;
    double x19 = x18*x3;
    double x20 = 1.0*n5;
    double x21 = x12 + x20;
    double x22 = n5 + x13;
    double x23 = 1.0/x22;
    double x24 = x15 - x22*x3;
    double x25 = x23*x24;
    double x26 = x19 + x21*x25;
    double x27 = x12*x17 + x26;
    double x28 = x27 + x8;
    double x29 = n2*x3;
    double x30 = 3.0*x29;
    double x31 = 3.0*n1;
    double x32 = 3.0*n4;
    double x33 = x31 + x32;
    double x34 = 1.0/x0;
    double x35 = -x0*x3 + x15;
    double x36 = x34*x35;
    double x37 = x30 + x33*x36;
    double x38 = x28 + x37;
    double x39 = -2*x3;
    double x40 = pow(x2, -3);
    double x41 = 2*x40;
    double x42 = x2*x33;
    double x43 = x34*x42;
    double x44 = 6.0*x2*x36 + x43*(x0*x41 + x39) - x35*x42/((x0)*(x0));
    double x45 = x38 + x44;
    double x46 = -x20;
    double x47 = 2.0*x2;
    double x48 = x12*x2;
    double x49 = x14*x48;
    double x50 = x2*x21;
    double x51 = x23*x50;
    double x52 = x25*x47 + x51*(x22*x41 + x39) - x24*x50/((x22)*(x22));
    double x53 = x17*x47 - x3*(-x18 + x46) + x49*(x13*x41 + x39) + x52 - x16*x48/((x13)*(x13));
    double x54 = 8.3144626181532395*T;
    double x55 = 0.82399999999999995*P + 168800.0;
    double x56 = 64.0*n3;
    double x57 = 0.10299999999999999*P + 21100.0;
    double x58 = n3*x57;
    double x59 = 23244800.0*n4 + 23449600.0*n5;
    double x60 = x18 + x20;
    double x61 = x12 + x60;
    double x62 = pow(x61, -2);
    double x63 = 0.0009765625*x62;
    double x64 = -x63*(x55*x56 + 512.0*x58 + x59);
    double x65 = 181600.0*n4 + 183200.0*n5;
    double x66 = 8.0*x57;
    double x67 = 64.0*n1;
    double x68 = 0.001953125*(64.0*n2*(168800.0*n3 + x65) + 64.0*n4*(181600.0*n1 + 181600.0*n2 + 488000.0*n3 + 568000.0*n5) + 64.0*n5*(183200.0*n1 + 183200.0*n2 + 485600.0*n3 + 568000.0*n4) + x56*(n1*x66 + 168800.0*n2 + 488000.0*n4 + 485600.0*n5) + x67*(n3*x66 + x65))/((x61)*(x61)*(x61));
    double x69 = 0.0625*x55;
    double x70 = 22700.0*n4 + 22900.0*n5;
    double x71 = 1.0*x62;
    double x72 = x68 - x71*(n3*x69 + 0.5*x58 + x70);
    double x73 = -x3;
    double x74 = -n1;
    double x75 = x43*(-x41*(-n4 + x74) + x73);
    double x76 = x38 + x75;
    double x77 = x54*(-6.0*x15 + x53 + x76);
    double x78 = -x63*(21606400.0*n3 + x59);
    double x79 = 1.0/x61;
    double x80 = n1*x57;
    double x81 = -x63*(21606400.0*n2 + 62464000.0*n4 + 62156800.0*n5 + x55*x67 + 512.0*x80);
    double x82 = x68 + x81;
    double x83 = 22700.0*x79;
    double x84 = -n2 - n3 + x74;
    double x85 = x3*x60 + x49*(-x41*x84 + x73);
    double x86 = x51*(-x41*(-n5 + x84) + x73) + x85;
    double x87 = -x63*(23244800.0*n1 + 23244800.0*n2 + 62464000.0*n3 + 72704000.0*n5);
    double x88 = x68 + x87;
    double x89 = 22900.0*x79;
    double x90 = x52 + x85;
    double x91 = -x63*(23449600.0*n1 + 23449600.0*n2 + 62156800.0*n3 + 72704000.0*n4);
    double x92 = x68 + x91;
    double x93 = -6.0*x3;
    double x94 = 6.0*x40;
    double x95 = n2*x94;
    double x96 = 3.0*x2;
    double x97 = x28 - x30;
    double x98 = -x3*(-x31 - x32);
    double x99 = 3.0*x15;
    double x100 = x53 + x98 + x99;
    double x101 = x68 + x78;
    double x102 = -x4;
    double x103 = x2*(x102 + x95) + x97;
    double x104 = x103 + x98;
    double x105 = -7.0*x15 + x3*x33 + x86;
    double x106 = -4.9999999969999998*x15 + x90;
    double x107 = n3*x94;
    double x108 = x27 + x30 - x5 + x7;
    double x109 = x108 + x2*(x102 + x107);
    double x110 = 2.0*n4*x40;
    double x111 = -x10 - x11 - x9;
    double x112 = n4 + n5;
    double x113 = 1.0/x112;
    double x114 = -x112*x3 + x15;
    double x115 = x113*x114;
    double x116 = x2*x60;
    double x117 = -x111*x3 + x113*x116*(x112*x41 + x39) + x115*x47 + x115*x60 - x114*x116/((x112)*(x112));
    double x118 = x117 - x19 + x37 + x8;

result[0] = x54*(x45 + x53) + x64 + x72;
result[1] = x72 + x77 + x78;
result[2] = x64 + x77 + 0.0009765625*x79*(105.47199999999999*P + 21606400.0) + x82;
result[3] = x54*(-4.0*x15 + x45 + x86) + x64 + x83 + x88;
result[4] = x54*(-7.9999999969999998*x15 + x76 + x90) + x64 + x89 + x92;
result[5] = x101 + x54*(x100 + x2*(x93 + x95) + x97 + x96*(x15 - x29)/n2) - x71*(21100.0*n3 + x70);
result[6] = x101 + x54*(x104 + x53 - x99) + 21100.0*x79 + x81;
result[7] = x101 + x54*(x103 + x105) + x83 + x87;
result[8] = x101 + x54*(x104 + x106) + x89 + x91;
result[9] = x54*(x100 + x108 + x2*(x107 + x93) + x96*(-n3*x3 + x15)/n3) - x71*(n1*x69 + 21100.0*n2 + 61000.0*n4 + 60700.0*n5 + 0.5*x80) + x82;
result[10] = x54*(x105 + x109) + 61000.0*x79 + x82 + x87;
result[11] = x54*(x106 + x109 + x98) + 60700.0*x79 + x82 + x91;
result[12] = x54*(x118 + 1.0*x15 + x2*(x110 - 2.0*x3) - x3*(x111 + x46) + x44 + 1.0*x2*(-n4*x3 + x15)/n4) - x71*(22700.0*n1 + 22700.0*n2 + 61000.0*n3 + 71000.0*n5) + x88;
result[13] = x54*(x118 - 6.9999999969999998*x15 + x2*(x110 - 1.0*x3) + x21*x3 + x75) + 71000.0*x79 + x88 + x91;
result[14] = x54*(x117 + 2.9999999970000002*x15 + x2*(5.9999999940000004*n5*x40 - 5.9999999940000004*x3) + x26 + x30 + x5 + x52 - x7 + x98 + 2.9999999970000002*x2*(x15 - x6)/n5) - x71*(22900.0*n1 + 22900.0*n2 + 60700.0*n3 + 71000.0*n4) + x92;
}
        
static void coder_d3gdn3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 0.82399999999999995*P + 168800.0;
    double x1 = n3*x0;
    double x2 = 0.10299999999999999*P + 21100.0;
    double x3 = n3*x2;
    double x4 = 23244800.0*n4 + 23449600.0*n5;
    double x5 = 1.0*n4;
    double x6 = 1.0*n5;
    double x7 = x5 + x6;
    double x8 = 1.0*n1;
    double x9 = 1.0*n2;
    double x10 = 1.0*n3;
    double x11 = x10 + x8 + x9;
    double x12 = x11 + x7;
    double x13 = pow(x12, -3);
    double x14 = 0.001953125*x13;
    double x15 = x14*(64.0*x1 + 512.0*x3 + x4);
    double x16 = n1 + n4;
    double x17 = 1.0/x16;
    double x18 = n2 + n3;
    double x19 = n5 + x16 + x18;
    double x20 = 1.0/x19;
    double x21 = pow(x19, -2);
    double x22 = -x16*x21 + x20;
    double x23 = x17*x22;
    double x24 = 3.0*n1;
    double x25 = 3.0*n4;
    double x26 = x24 + x25;
    double x27 = pow(x16, -2);
    double x28 = x22*x27;
    double x29 = x26*x28;
    double x30 = -2*x21;
    double x31 = pow(x19, -3);
    double x32 = 2*x31;
    double x33 = x16*x32 + x30;
    double x34 = x17*x26;
    double x35 = x33*x34;
    double x36 = 9.0*x19;
    double x37 = 6*x31;
    double x38 = pow(x19, -4);
    double x39 = 6*x38;
    double x40 = x19*x34;
    double x41 = x26*x27;
    double x42 = 2*x19;
    double x43 = x17*x33*x36 + 9.0*x23 - x28*x36 - 2*x29 - x33*x41*x42 + 2*x35 + x40*(-x16*x39 + x37) + x22*x26*x42/((x16)*(x16)*(x16));
    double x44 = n1 + x18;
    double x45 = n5 + x44;
    double x46 = 1.0/x45;
    double x47 = x20 - x21*x45;
    double x48 = x46*x47;
    double x49 = x11 + x6;
    double x50 = pow(x45, -2);
    double x51 = x47*x50;
    double x52 = x49*x51;
    double x53 = x30 + x32*x45;
    double x54 = x46*x49;
    double x55 = x53*x54;
    double x56 = 3.0*x19;
    double x57 = x19*x54;
    double x58 = x19*x49;
    double x59 = x50*x58;
    double x60 = x46*x53*x56 + 3.0*x48 - x51*x56 - 2*x52 - 2*x53*x59 + 2*x55 + x57*(x37 - x39*x45) + 2*x47*x58/((x45)*(x45)*(x45));
    double x61 = -x32*x7;
    double x62 = 6.0*x31;
    double x63 = -n2*x62;
    double x64 = 2.0*x31;
    double x65 = -n4*x64;
    double x66 = -n3*x62;
    double x67 = n5*x31;
    double x68 = -5.9999999940000004*x67;
    double x69 = x66 + x68;
    double x70 = x65 + x69;
    double x71 = x63 + x70;
    double x72 = x61 + x71;
    double x73 = 1.0/x44;
    double x74 = x20 - x21*x44;
    double x75 = x73*x74;
    double x76 = pow(x44, -2);
    double x77 = x74*x76;
    double x78 = x11*x77;
    double x79 = x30 + x32*x44;
    double x80 = x11*x73;
    double x81 = x79*x80;
    double x82 = x19*x80;
    double x83 = x11*x19;
    double x84 = x76*x83;
    double x85 = x56*x73*x79 - x56*x77 + 3.0*x75 - 2*x78 - 2*x79*x84 + 2*x81 + x82*(x37 - x39*x44) + 2*x74*x83/((x44)*(x44)*(x44));
    double x86 = x72 + x85;
    double x87 = 8.3144626181532395*T;
    double x88 = 0.0625*x1;
    double x89 = 0.5*x3;
    double x90 = 22700.0*n4;
    double x91 = 22900.0*n5;
    double x92 = -x90 - x91;
    double x93 = x13*(-x88 - x89 + x92);
    double x94 = 181600.0*n4 + 183200.0*n5;
    double x95 = 64.0*n1;
    double x96 = n1*x2;
    double x97 = -0.005859375*(64.0*n2*(168800.0*n3 + x94) + 64.0*n3*(168800.0*n2 + 488000.0*n4 + 485600.0*n5 + 8.0*x96) + 64.0*n4*(181600.0*n1 + 181600.0*n2 + 488000.0*n3 + 568000.0*n5) + 64.0*n5*(183200.0*n1 + 183200.0*n2 + 485600.0*n3 + 568000.0*n4) + x95*(8.0*x3 + x94))/((x12)*(x12)*(x12)*(x12));
    double x98 = -4.0*x93 + x97;
    double x99 = x14*(21606400.0*n3 + x4);
    double x100 = -x21;
    double x101 = -n1;
    double x102 = -n4 + x101;
    double x103 = x100 - x102*x32;
    double x104 = x103*x19;
    double x105 = 4*x31;
    double x106 = 2*n1;
    double x107 = 2*n4;
    double x108 = 3*x38;
    double x109 = -x108*(x106 + x107);
    double x110 = x103*x34;
    double x111 = -x104*x41 + x110 - x29 + x35;
    double x112 = 6.0*x104*x17 + x111 + 6.0*x23 + x40*(x105 + x109);
    double x113 = 3.0*x21 + x60;
    double x114 = x87*(x112 + x113 + x86);
    double x115 = x90 + x91;
    double x116 = 2.0*x13;
    double x117 = x116*(x115 + x88 + x89);
    double x118 = x117 + x15;
    double x119 = pow(x12, -2);
    double x120 = -1.0*x119*x2;
    double x121 = x14*(21606400.0*n2 + 62464000.0*n4 + 62156800.0*n5 + x0*x95 + 512.0*x96);
    double x122 = x121 + x97;
    double x123 = x120 + x122;
    double x124 = -0.0009765625*x119*(105.47199999999999*P + 21606400.0) + x123;
    double x125 = 2.0*x19;
    double x126 = -n2 - n3 + x101;
    double x127 = -n5 + x126;
    double x128 = x100 - x127*x32;
    double x129 = x128*x46;
    double x130 = 2*n5;
    double x131 = 2*n2;
    double x132 = 2*n3;
    double x133 = x106 + x131 + x132;
    double x134 = -x108*(x130 + x133);
    double x135 = x128*x54;
    double x136 = -x128*x59 + x135 - x52 + x55;
    double x137 = x125*x129 + x136 + 2.0*x48 + x57*(x105 + x134);
    double x138 = -x6;
    double x139 = x100 - x126*x32;
    double x140 = x139*x80;
    double x141 = -x108*x133;
    double x142 = x125*x139*x73 - x139*x84 + x140 + x32*(x138 - x5) + 2.0*x75 - x78 + x81 + x82*(x105 + x141);
    double x143 = x142 + x71;
    double x144 = x137 + x143;
    double x145 = -45400.0*x119;
    double x146 = x14*(23244800.0*n1 + 23244800.0*n2 + 62464000.0*n3 + 72704000.0*n5);
    double x147 = x146 + x97;
    double x148 = x145 + x147;
    double x149 = 3.9999999970000002*x21;
    double x150 = -45800.0*x119;
    double x151 = x14*(23449600.0*n1 + 23449600.0*n2 + 62156800.0*n3 + 72704000.0*n4);
    double x152 = x151 + x97;
    double x153 = x150 + x152;
    double x154 = 21100.0*n3;
    double x155 = -x154 + x92;
    double x156 = 9.0*x21;
    double x157 = 2*x110 + x40*(x109 + x32);
    double x158 = x157 + x60;
    double x159 = x87*(x156 + x158 + x86);
    double x160 = x97 + x99;
    double x161 = x159 + x160;
    double x162 = x117 + x160;
    double x163 = x146 + x87*(x103*x17*x56 + x111 + x144 + 8.0*x21 + 3.0*x23 + x40*(x102*x39 + x105));
    double x164 = x87*(x143 + x158 + 9.9999999969999998*x21);
    double x165 = x150 + x151;
    double x166 = 21100.0*n2;
    double x167 = 61000.0*n4;
    double x168 = 60700.0*n5;
    double x169 = 0.5*x96;
    double x170 = 0.0625*n1*x0;
    double x171 = -x166 - x167 - x168 - x169 - x170;
    double x172 = -x116*x171;
    double x173 = x122 + x151;
    double x174 = 6.0*x21;
    double x175 = 2*x135 + x57*(x134 + x32);
    double x176 = 2*x140;
    double x177 = x82*(x141 + x32);
    double x178 = x176 + x177 + x72;
    double x179 = 22700.0*n1;
    double x180 = 22700.0*n2;
    double x181 = 61000.0*n3;
    double x182 = 71000.0*n5;
    double x183 = -x179 - x180 - x181 - x182;
    double x184 = -x116*x183;
    double x185 = 1.0*x19;
    double x186 = x129*x185 + x136 + 7.9999999969999998*x21 + 1.0*x48 + x57*(x105 + x127*x39);
    double x187 = -116600.0*x119 + x151;
    double x188 = -22900.0*n1 - 22900.0*n2 - 60700.0*n3 - 71000.0*n4;
    double x189 = -x116*x188;
    double x190 = 18.0*x31;
    double x191 = 18.0*x38;
    double x192 = -n2*x191;
    double x193 = x131*x31;
    double x194 = 3.0/n2;
    double x195 = x19*x194;
    double x196 = -n2*x21 + x20;
    double x197 = 12.0*x31;
    double x198 = n2*x197 + x70;
    double x199 = x194*x196 + x198;
    double x200 = -x26*x32;
    double x201 = x200 + x61;
    double x202 = x201 + x85;
    double x203 = x202 + x60;
    double x204 = x203 - 12.0*x21;
    double x205 = 4.0*x13;
    double x206 = x19*(x192 + x197) + x195*(x100 + x193) + x199;
    double x207 = x116*(x115 + x154) + x160;
    double x208 = -42200.0*x119 + x121;
    double x209 = -4.0*x21;
    double x210 = x137 + x32*(-x24 - x25);
    double x211 = x142 + x209 + x210;
    double x212 = x145 + x146;
    double x213 = x200 + x60;
    double x214 = x142 - 5.0000000030000002*x21 + x213;
    double x215 = x19*(x192 + x62) + x198;
    double x216 = x142 + x215;
    double x217 = x121 + x160;
    double x218 = x176 + x177 + x201;
    double x219 = x156 + x175 + x218;
    double x220 = x215 + x218;
    double x221 = 5.9999999939999995*x21 + x60;
    double x222 = -n3*x191;
    double x223 = x132*x31;
    double x224 = 3.0/n3;
    double x225 = x19*x224;
    double x226 = -n3*x21 + x20;
    double x227 = x63 + x65;
    double x228 = n3*x197 + x227 + x68;
    double x229 = x224*x226 + x228;
    double x230 = x19*(x197 + x222) + x225*(x100 + x223) + x229;
    double x231 = x116*(x166 + x167 + x168 + x169 + x170) + x122;
    double x232 = -122000.0*x119 + x146;
    double x233 = -121400.0*x119 + x151;
    double x234 = x19*(x222 + x62) + x228;
    double x235 = x218 + x234;
    double x236 = -n4*x21 + x20;
    double x237 = 1.0/n4;
    double x238 = x236*x237;
    double x239 = -6.0*n4*x38;
    double x240 = x107*x31;
    double x241 = x19*x237;
    double x242 = 4.0*x31;
    double x243 = n4 + n5;
    double x244 = x20 - x21*x243;
    double x245 = 1.0/x243;
    double x246 = 3.0*x245;
    double x247 = 2*x7;
    double x248 = pow(x243, -2);
    double x249 = x244*x248;
    double x250 = x243*x32 + x30;
    double x251 = x247*x250;
    double x252 = -x11*x32 + x19*x245*x7*(-x243*x39 + x37) + x19*x246*x250 - x19*x248*x251 + x19*x244*x247/((x243)*(x243)*(x243)) + x244*x246 + x245*x251 - x247*x249 - x249*x56;
    double x253 = n4*x242 + x252 + x63 + x69;
    double x254 = x253 - x32*x49;
    double x255 = -142000.0*x119 + x147 + x151;
    double x256 = -n5*x21 + x20;
    double x257 = 2.9999999970000002/n5;

result[0] = x15 + x87*(x43 + x60 + x86) + x98;
result[1] = x114 + x98 + x99;
result[2] = x114 + x118 + x124;
result[3] = x118 + x148 + x87*(x144 + 2.0*x21 + x43);
result[4] = x118 + x153 + x87*(x112 + x143 + x149 + x60);
result[5] = -x116*x155 + x161 - 2.0*x93;
result[6] = x117 - 21100.0*x119 + x120 + x121 + x161;
result[7] = x145 + x162 + x163;
result[8] = x162 + x164 + x165;
result[9] = x124 + x15 + x159 + x172;
result[10] = -83700.0*x119 + x123 + x15 + x163;
result[11] = -83600.0*x119 + x120 + x15 + x164 + x173;
result[12] = x148 + x15 + x184 + x87*(x174 + x175 + x178 + x43);
result[13] = x147 + x15 + x187 + x87*(x112 + x178 + x186);
result[14] = x15 + x153 + x189 + x87*(x158 + x178 + 11.999999994*x21);
result[15] = -x155*x205 + x160 + x87*(x19*(x190 + x192) + x195*(x193 + x30) + x199 + x204 - x196*x56/((n2)*(n2)));
result[16] = x207 + x208 + x87*(-x174 + x203 + x206);
result[17] = x207 + x212 + x87*(x206 + x211);
result[18] = x165 + x207 + x87*(x206 + x214);
result[19] = x160 + x172 + x208 + x87*(x113 + x202 + x215);
result[20] = -104800.0*x119 + x146 + x217 + x87*(5.0*x21 + x210 + x216);
result[21] = -104700.0*x119 + x151 + x217 + x87*(x149 + x213 + x216);
result[22] = x160 + x184 + x212 + x87*(x215 + x219);
result[23] = x146 + x160 + x187 + x87*(x186 + x220);
result[24] = x160 + x165 + x189 + x87*(x220 + x221);
result[25] = x122 - x171*x205 + x87*(x19*(x190 + x222) + x204 + x225*(x223 + x30) + x229 - x226*x56/((n3)*(n3)));
result[26] = x231 + x232 + x87*(x211 + x230);
result[27] = x231 + x233 + x87*(x214 + x230);
result[28] = x122 + x184 + x232 + x87*(x219 + x234);
result[29] = -192700.0*x119 + x146 + x173 + x87*(x186 + x235);
result[30] = x122 + x189 + x233 + x87*(x221 + x235);
result[31] = x147 - x183*x205 + x87*(x19*(x239 + x62) + x209 + x238 + x241*(x240 + x30) + x254 + x43 - x185*x236/((n4)*(n4)));
result[32] = x116*(x179 + x180 + x181 + x182) + x255 + x87*(x112 + x19*(x239 + x242) + 0.9999999970000002*x21 + x238 + x241*(x100 + x240) + x253 + x32*(-x10 + x138 - x8 - x9));
result[33] = x189 + x255 + x87*(x157 + x19*(x239 + x64) + 9.9999999939999995*x21 + x254);
result[34] = x152 - x188*x205 + x87*(x19*x257*(x130*x31 + x30) + x19*(-17.999999982000002*n5*x38 + 17.999999982000002*x31) - 11.999999988000001*x21 + x213 + x227 + x252 + x256*x257 + x66 + 11.999999988000001*x67 - 2.9999999970000002*x19*x256/((n5)*(n5)));
}
        
static double coder_dgdt(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = n4 + n5;
    double x1 = n1 + n2 + n3;
    double x2 = 1.0/(x0 + x1);
    double x3 = log(n5*x2);
    double x4 = 1.0*n5;
    double x5 = 1.0*n1 + 1.0*n2 + 1.0*n3;
    double x6 = sqrt(-5.3525400000000003*T + 40.14405);
    double x7 = 0.15782994586457105*x6;
    double x8 = fmin(4, x7);
    double x9 = (-x7 + 4 >= 0. ? 1. : 0.)/x6;

result = n1*(*endmember[0].dmu0dT)(T, P) + 24.943387854459719*n2*log(n2*x2) + n2*(*endmember[1].dmu0dT)(T, P) + 24.943387854459719*n3*log(n3*x2) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*log(n4*x2) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626098387775*n5*(x3 - 1.0986122896681101) + 16.628925219677555*n5*(x3 - 0.40546510910816402) + n5*(*endmember[4].dmu0dT)(T, P) + 8.3144626181532395*x5*log(x1*x2) + 8.3144626181532395*(3.0*n1 + 3.0*n4)*log(x2*(n1 + n4)) + 8.3144626181532395*(1.0*n4 + x4)*log(x0*x2) + 8.3144626181532395*(x4 + x5)*log(x2*(n5 + x1)) + ((T >= 7.5) ? (
   -40.14405*n2
)
: (
   (1.0/3.0)*n2*(-381.52503107154041*((x8)*(x8))*x9 + 120.43215000000001*x8 - 0.4223955492189756*x9*(120.43215000000001*T - 903.24112500000001) - 120.43215000000001)
));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = n2 + n3;
    double x2 = n5 + x0 + x1;
    double x3 = 1.0/x2;
    double x4 = n5*x3;
    double x5 = -24.943387829516332*x4;
    double x6 = n2*x3;
    double x7 = -24.943387854459719*x6;
    double x8 = n3*x3;
    double x9 = -24.943387854459719*x8;
    double x10 = 24.943387854459719*n1 + 24.943387854459719*n4;
    double x11 = pow(x2, -2);
    double x12 = x5 + x7 + x9 + 24.943387854459719*log(x0*x3) + x10*x2*(-x0*x11 + x3)/x0;
    double x13 = n1 + x1;
    double x14 = 8.3144626181532395*n4;
    double x15 = 8.3144626181532395*n5;
    double x16 = x14 + x15;
    double x17 = 8.3144626181532395*n1 + 8.3144626181532395*n2 + 8.3144626181532395*n3;
    double x18 = n5 + x13;
    double x19 = x15 + x17;
    double x20 = -x14*x3 + 8.3144626181532395*log(x18*x3) + x19*x2*(-x11*x18 + x3)/x18;
    double x21 = -x16*x3 + x20 + 8.3144626181532395*log(x13*x3) + x17*x2*(-x11*x13 + x3)/x13;
    double x22 = -x10*x3;
    double x23 = 24.943387854459719*x2;
    double x24 = sqrt(-5.3525400000000003*T + 40.14405);
    double x25 = 0.15782994586457105*x24;
    double x26 = fmin(4, x25);
    double x27 = (-x25 + 4 >= 0. ? 1. : 0.)/x24;
    double x28 = x21 + x5;
    double x29 = x22 + x7;
    double x30 = n4 + n5;
    double x31 = x16*x2*(-x11*x30 + x3)/x30 - x17*x3 + 8.3144626181532395*log(x3*x30);

result[0] = x12 + x21 + (*endmember[0].dmu0dT)(T, P);
result[1] = x22 + x23*(-n2*x11 + x3) + x28 + x9 + ((T >= 7.5) ? (
   -40.14405
)
: (
   -127.17501035718013*((x26)*(x26))*x27 + 40.14405*x26 - 0.14079851640632518*x27*(120.43215000000001*T - 903.24112500000001) - 40.14405
)) + 24.943387854459719*log(x6) + (*endmember[1].dmu0dT)(T, P);
result[2] = x23*(-n3*x11 + x3) + x28 + x29 + 24.943387854459719*log(x8) + (*endmember[2].dmu0dT)(T, P);
result[3] = x12 - x19*x3 + 8.3144626181532395*x2*(-n4*x11 + x3) + x31 + 8.3144626181532395*log(n4*x3) + (*endmember[3].dmu0dT)(T, P);
result[4] = 24.943387829516332*x2*(-n5*x11 + x3) + x20 + x29 + x31 + x9 + 24.943387829516332*log(x4) + (*endmember[4].dmu0dT)(T, P) - 15.876819783702929;
}
        
static void coder_d3gdn2dt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = n2 + n3;
    double x2 = n5 + x0 + x1;
    double x3 = pow(x2, -2);
    double x4 = n5*x3;
    double x5 = 24.943387829516332*x4;
    double x6 = 24.943387854459719*x3;
    double x7 = n3*x6;
    double x8 = x5 + x7;
    double x9 = 8.3144626181532395*n1;
    double x10 = 8.3144626181532395*n2;
    double x11 = 8.3144626181532395*n3;
    double x12 = x10 + x11 + x9;
    double x13 = n1 + x1;
    double x14 = 1.0/x13;
    double x15 = 1.0/x2;
    double x16 = -x13*x3 + x15;
    double x17 = x14*x16;
    double x18 = 8.3144626181532395*n4;
    double x19 = x18*x3;
    double x20 = 8.3144626181532395*n5;
    double x21 = x12 + x20;
    double x22 = n5 + x13;
    double x23 = 1.0/x22;
    double x24 = x15 - x22*x3;
    double x25 = x23*x24;
    double x26 = x19 + x21*x25;
    double x27 = x12*x17 + x26;
    double x28 = x27 + x8;
    double x29 = n2*x3;
    double x30 = 24.943387854459719*x29;
    double x31 = 24.943387854459719*n1;
    double x32 = 24.943387854459719*n4;
    double x33 = x31 + x32;
    double x34 = 1.0/x0;
    double x35 = -x0*x3 + x15;
    double x36 = x34*x35;
    double x37 = x30 + x33*x36;
    double x38 = x28 + x37;
    double x39 = -2*x3;
    double x40 = pow(x2, -3);
    double x41 = 2*x40;
    double x42 = x2*x33;
    double x43 = x34*x42;
    double x44 = 49.886775708919437*x2*x36 + x43*(x0*x41 + x39) - x35*x42/((x0)*(x0));
    double x45 = x38 + x44;
    double x46 = -x20;
    double x47 = 16.628925236306479*x2;
    double x48 = x12*x2;
    double x49 = x14*x48;
    double x50 = x2*x21;
    double x51 = x23*x50;
    double x52 = x25*x47 + x51*(x22*x41 + x39) - x24*x50/((x22)*(x22));
    double x53 = x17*x47 - x3*(-x18 + x46) + x49*(x13*x41 + x39) + x52 - x16*x48/((x13)*(x13));
    double x54 = -x3;
    double x55 = -n1;
    double x56 = x43*(-x41*(-n4 + x55) + x54);
    double x57 = x38 + x56;
    double x58 = -49.886775708919437*x15 + x53 + x57;
    double x59 = -n2 - n3 + x55;
    double x60 = x18 + x20;
    double x61 = x3*x60 + x49*(-x41*x59 + x54);
    double x62 = x51*(-x41*(-n5 + x59) + x54) + x61;
    double x63 = x52 + x61;
    double x64 = -49.886775708919437*x3;
    double x65 = 49.886775708919437*x40;
    double x66 = n2*x65;
    double x67 = 24.943387854459719*x2;
    double x68 = x28 - x30;
    double x69 = -x3*(-x31 - x32);
    double x70 = 24.943387854459719*x15;
    double x71 = x53 + x69 + x70;
    double x72 = -x6;
    double x73 = x2*(x66 + x72) + x68;
    double x74 = x69 + x73;
    double x75 = -58.201238327072687*x15 + x3*x33 + x62;
    double x76 = -41.572313065822811*x15 + x63;
    double x77 = n3*x65;
    double x78 = x27 + x30 + x5 - x7;
    double x79 = x2*(x72 + x77) + x78;
    double x80 = 16.628925236306479*n4*x40;
    double x81 = -x10 - x11 - x9;
    double x82 = n4 + n5;
    double x83 = 1.0/x82;
    double x84 = x15 - x3*x82;
    double x85 = x83*x84;
    double x86 = x2*x60;
    double x87 = -x3*x81 + x47*x85 + x60*x85 + x83*x86*(x39 + x41*x82) - x84*x86/((x82)*(x82));
    double x88 = -x19 + x37 + x8 + x87;

result[0] = x45 + x53;
result[1] = x58;
result[2] = x58;
result[3] = -33.257850472612958*x15 + x45 + x62;
result[4] = -66.515700920282541*x15 + x57 + x63;
result[5] = x2*(x64 + x66) + x68 + x71 + x67*(x15 - x29)/n2;
result[6] = x53 - x70 + x74;
result[7] = x73 + x75;
result[8] = x74 + x76;
result[9] = x2*(x64 + x77) + x71 + x78 + x67*(-n3*x3 + x15)/n3;
result[10] = x75 + x79;
result[11] = x69 + x76 + x79;
result[12] = 8.3144626181532395*x15 + x2*(-16.628925236306479*x3 + x80) - x3*(x46 + x81) + x44 + x88 + 8.3144626181532395*x2*(-n4*x3 + x15)/n4;
result[13] = -58.201238302129291*x15 + x2*(-8.3144626181532395*x3 + x80) + x21*x3 + x56 + x88;
result[14] = 24.943387829516332*x15 + x2*(49.886775659032665*n5*x40 - 49.886775659032665*x3) + x26 + x30 - x5 + x52 + x69 + x7 + x87 + 24.943387829516332*x2*(x15 - x4)/n5;
}
        
static void coder_d4gdn3dt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = 1.0/x0;
    double x2 = n2 + n3;
    double x3 = n5 + x0 + x2;
    double x4 = 1.0/x3;
    double x5 = pow(x3, -2);
    double x6 = -x0*x5 + x4;
    double x7 = x1*x6;
    double x8 = 24.943387854459719*n1;
    double x9 = 24.943387854459719*n4;
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
    double x20 = 74.830163563379159*x3;
    double x21 = 6*x15;
    double x22 = pow(x3, -4);
    double x23 = 6*x22;
    double x24 = x18*x3;
    double x25 = x10*x11;
    double x26 = 2*x3;
    double x27 = x1*x17*x20 - x12*x20 - 2*x13 - x17*x25*x26 + 2*x19 + x24*(-x0*x23 + x21) + 74.830163563379159*x7 + x10*x26*x6/((x0)*(x0)*(x0));
    double x28 = n1 + x2;
    double x29 = n5 + x28;
    double x30 = 1.0/x29;
    double x31 = -x29*x5 + x4;
    double x32 = x30*x31;
    double x33 = 8.3144626181532395*n5;
    double x34 = 8.3144626181532395*n1;
    double x35 = 8.3144626181532395*n2;
    double x36 = 8.3144626181532395*n3;
    double x37 = x34 + x35 + x36;
    double x38 = x33 + x37;
    double x39 = pow(x29, -2);
    double x40 = x31*x39;
    double x41 = x38*x40;
    double x42 = x14 + x16*x29;
    double x43 = x30*x38;
    double x44 = x42*x43;
    double x45 = 24.943387854459719*x3;
    double x46 = x3*x43;
    double x47 = x3*x38;
    double x48 = x39*x47;
    double x49 = x30*x42*x45 + 24.943387854459719*x32 - x40*x45 - 2*x41 - 2*x42*x48 + 2*x44 + x46*(x21 - x23*x29) + 2*x31*x47/((x29)*(x29)*(x29));
    double x50 = 49.886775708919437*x15;
    double x51 = -n2*x50;
    double x52 = 8.3144626181532395*n4;
    double x53 = x33 + x52;
    double x54 = -x16*x53;
    double x55 = 16.628925236306479*x15;
    double x56 = -n4*x55;
    double x57 = n5*x15;
    double x58 = -49.886775659032665*x57;
    double x59 = -n3*x50;
    double x60 = x58 + x59;
    double x61 = x56 + x60;
    double x62 = x54 + x61;
    double x63 = x51 + x62;
    double x64 = 1.0/x28;
    double x65 = -x28*x5 + x4;
    double x66 = x64*x65;
    double x67 = pow(x28, -2);
    double x68 = x65*x67;
    double x69 = x37*x68;
    double x70 = x14 + x16*x28;
    double x71 = x37*x64;
    double x72 = x70*x71;
    double x73 = x3*x71;
    double x74 = x3*x37;
    double x75 = x67*x74;
    double x76 = x45*x64*x70 - x45*x68 + 24.943387854459719*x66 - 2*x69 - 2*x70*x75 + 2*x72 + x73*(x21 - x23*x28) + 2*x65*x74/((x28)*(x28)*(x28));
    double x77 = x63 + x76;
    double x78 = 24.943387854459719*x5;
    double x79 = -x5;
    double x80 = -n1;
    double x81 = -n4 + x80;
    double x82 = -x16*x81 + x79;
    double x83 = x3*x82;
    double x84 = 4*x15;
    double x85 = 2*n1;
    double x86 = 2*n4;
    double x87 = 3*x22;
    double x88 = -x87*(x85 + x86);
    double x89 = x18*x82;
    double x90 = -x13 + x19 - x25*x83 + x89;
    double x91 = 49.886775708919437*x1*x83 + x24*(x84 + x88) + 49.886775708919437*x7 + x90;
    double x92 = x49 + x91;
    double x93 = x77 + x78 + x92;
    double x94 = 16.628925236306479*x32;
    double x95 = 16.628925236306479*x3;
    double x96 = -n2 - n3 + x80;
    double x97 = -n5 + x96;
    double x98 = -x16*x97 + x79;
    double x99 = x30*x98;
    double x100 = x95*x99;
    double x101 = 2*n5;
    double x102 = 2*n2;
    double x103 = 2*n3;
    double x104 = x102 + x103 + x85;
    double x105 = -x87*(x101 + x104);
    double x106 = x46*(x105 + x84);
    double x107 = x43*x98;
    double x108 = x107 - x41 + x44 - x48*x98;
    double x109 = -x33;
    double x110 = -x16*x96 + x79;
    double x111 = x110*x71;
    double x112 = -x104*x87;
    double x113 = x110*x64*x95 - x110*x75 + x111 + x16*(x109 - x52) + 16.628925236306479*x66 - x69 + x72 + x73*(x112 + x84);
    double x114 = x113 + x61;
    double x115 = x114 + x51;
    double x116 = x100 + x106 + x108 + x115 + x94;
    double x117 = 33.257850447669568*x5;
    double x118 = x24*(x16 + x88) + 2*x89;
    double x119 = x118 + x49;
    double x120 = x119 + 74.830163563379159*x5 + x77;
    double x121 = x1*x45*x82 + x116 + x24*(x23*x81 + x84) + 66.515700945225916*x5 + 24.943387854459719*x7 + x90;
    double x122 = x115 + x119 + 83.144626156588998*x5;
    double x123 = 49.886775708919444*x5;
    double x124 = 2*x107 + x46*(x105 + x16);
    double x125 = 2*x111;
    double x126 = x73*(x112 + x16);
    double x127 = x125 + x126 + x63;
    double x128 = 8.3144626181532395*x3;
    double x129 = x108 + x128*x99 + 8.3144626181532395*x32 + x46*(x23*x97 + x84) + 66.515700920282541*x5;
    double x130 = -99.773551417838874*x5;
    double x131 = 149.66032712675832*x15;
    double x132 = 149.66032712675832*x22;
    double x133 = -n2*x132;
    double x134 = x102*x15;
    double x135 = 24.943387854459719/n2;
    double x136 = x135*x3;
    double x137 = -n2*x5 + x4;
    double x138 = 99.773551417838874*x15;
    double x139 = n2*x138;
    double x140 = x135*x137 + x139;
    double x141 = -x10*x16;
    double x142 = x141 + x49;
    double x143 = x142 + x76;
    double x144 = x143 + x62;
    double x145 = x136*(x134 + x79) + x140 + x3*(x133 + x138);
    double x146 = -33.257850472612958*x5;
    double x147 = x100 + x106 + x108 + x16*(-x8 - x9) + x94;
    double x148 = x146 + x147;
    double x149 = x114 + x145;
    double x150 = x142 - 41.572313115709584*x5;
    double x151 = x139 + x3*(x133 + x50);
    double x152 = x151 + x62;
    double x153 = x114 + x151;
    double x154 = x125 + x126 + x141;
    double x155 = x124 + x154 + 74.830163563379173*x5;
    double x156 = x152 + x154;
    double x157 = x49 + 49.886775659032665*x5;
    double x158 = -n3*x5 + x4;
    double x159 = 24.943387854459719/n3;
    double x160 = x158*x159;
    double x161 = -n3*x132;
    double x162 = x103*x15;
    double x163 = x159*x3;
    double x164 = x51 + x56;
    double x165 = n3*x138 + x164 + x58;
    double x166 = x165 + x54;
    double x167 = x113 + x160 + x163*(x162 + x79) + x165 + x3*(x138 + x161);
    double x168 = x166 + x3*(x161 + x50);
    double x169 = x154 + x168;
    double x170 = -n4*x5 + x4;
    double x171 = 8.3144626181532395/n4;
    double x172 = x170*x171;
    double x173 = -49.886775708919437*n4*x22;
    double x174 = x15*x86;
    double x175 = x171*x3;
    double x176 = 33.257850472612958*x15;
    double x177 = n4 + n5;
    double x178 = -x177*x5 + x4;
    double x179 = 1.0/x177;
    double x180 = 24.943387854459719*x179;
    double x181 = 2*x53;
    double x182 = pow(x177, -2);
    double x183 = x178*x182;
    double x184 = x14 + x16*x177;
    double x185 = x181*x184;
    double x186 = -x16*x37 + x178*x180 + x179*x185 + x179*x3*x53*(-x177*x23 + x21) + x180*x184*x3 - x181*x183 - x182*x185*x3 - x183*x45 + x178*x181*x3/((x177)*(x177)*(x177));
    double x187 = n4*x176 + x186 + x51 + x60;
    double x188 = -x16*x38 + x187;
    double x189 = -n5*x5 + x4;
    double x190 = 24.943387829516332/n5;

result[0] = x27 + x49 + x77;
result[1] = x93;
result[2] = x93;
result[3] = x116 + x27 + 16.628925236306479*x5;
result[4] = x115 + x117 + x92;
result[5] = x120;
result[6] = x120;
result[7] = x121;
result[8] = x122;
result[9] = x120;
result[10] = x121;
result[11] = x122;
result[12] = x123 + x124 + x127 + x27;
result[13] = x127 + x129 + x91;
result[14] = x119 + x127 + 99.773551367952109*x5;
result[15] = x130 + x136*(x134 + x14) + x140 + x144 + x3*(x131 + x133) - x137*x45/((n2)*(n2));
result[16] = -x123 + x144 + x145;
result[17] = x148 + x149;
result[18] = x149 + x150;
result[19] = x143 + x152 + x78;
result[20] = x147 + x153 + 41.572313090766201*x5;
result[21] = x117 + x142 + x153;
result[22] = x152 + x155;
result[23] = x129 + x156;
result[24] = x156 + x157;
result[25] = x130 + x143 + x160 + x163*(x14 + x162) + x166 + x3*(x131 + x161) - x158*x45/((n3)*(n3));
result[26] = x148 + x167;
result[27] = x150 + x167;
result[28] = x155 + x168;
result[29] = x129 + x169;
result[30] = x157 + x169;
result[31] = x146 + x172 + x175*(x14 + x174) + x188 + x27 + x3*(x173 + x50) - x128*x170/((n4)*(n4));
result[32] = x16*(x109 - x34 - x35 - x36) + x172 + x175*(x174 + x79) + x187 + x3*(x173 + x176) + 8.3144625932098535*x5 + x91;
result[33] = x118 + x188 + x3*(x173 + x55) + 83.144626131645623*x5;
result[34] = x142 + x164 + x186 + x189*x190 + x190*x3*(x101*x15 + x14) + x3*(-149.66032697709801*n5*x22 + 149.66032697709801*x15) - 99.77355131806533*x5 + 99.77355131806533*x57 + x59 - 24.943387829516332*x189*x3/((n5)*(n5));
}
        
static double coder_dgdp(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = 0.10299999999999999*n1*n3/(1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5) + n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P) + n4*(*endmember[3].dmu0dP)(T, P) + n5*(*endmember[4].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x1 = 0.10299999999999999/x0;
    double x2 = -0.10299999999999999*n1*n3/((x0)*(x0));

result[0] = n3*x1 + x2 + (*endmember[0].dmu0dP)(T, P);
result[1] = x2 + (*endmember[1].dmu0dP)(T, P);
result[2] = n1*x1 + x2 + (*endmember[2].dmu0dP)(T, P);
result[3] = x2 + (*endmember[3].dmu0dP)(T, P);
result[4] = x2 + (*endmember[4].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x1 = pow(x0, -2);
    double x2 = n3*x1;
    double x3 = 0.20599999999999999*n1*n3/((x0)*(x0)*(x0));
    double x4 = -0.10299999999999999*x2 + x3;
    double x5 = n1*x1;
    double x6 = -0.10299999999999999*x5;
    double x7 = x3 + x6;

result[0] = -0.20599999999999999*x2 + x3;
result[1] = x4;
result[2] = x4 + x6 + 0.10299999999999999/x0;
result[3] = x4;
result[4] = x4;
result[5] = x3;
result[6] = x7;
result[7] = x3;
result[8] = x3;
result[9] = x3 - 0.20599999999999999*x5;
result[10] = x7;
result[11] = x7;
result[12] = x3;
result[13] = x3;
result[14] = x3;
}
        
static void coder_d4gdn3dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x1 = pow(x0, -3);
    double x2 = n3*x1;
    double x3 = -0.61799999999999999*n1*n3/((x0)*(x0)*(x0)*(x0));
    double x4 = 0.41199999999999998*x2 + x3;
    double x5 = n1*x1;
    double x6 = 0.20599999999999999*x5;
    double x7 = pow(x0, -2);
    double x8 = -0.20599999999999999*x7;
    double x9 = 0.20599999999999999*x2 + x3;
    double x10 = x6 - 0.10299999999999999*x7 + x9;
    double x11 = 0.41199999999999998*x5;
    double x12 = x3 + x6;
    double x13 = x11 + x3;

result[0] = 0.61799999999999999*x2 + x3;
result[1] = x4;
result[2] = x4 + x6 + x8;
result[3] = x4;
result[4] = x4;
result[5] = x9;
result[6] = x10;
result[7] = x9;
result[8] = x9;
result[9] = x11 + x8 + x9;
result[10] = x10;
result[11] = x10;
result[12] = x9;
result[13] = x9;
result[14] = x9;
result[15] = x3;
result[16] = x12;
result[17] = x3;
result[18] = x3;
result[19] = x13;
result[20] = x12;
result[21] = x12;
result[22] = x3;
result[23] = x3;
result[24] = x3;
result[25] = x3 + 0.61799999999999999*x5;
result[26] = x13;
result[27] = x13;
result[28] = x12;
result[29] = x12;
result[30] = x12;
result[31] = x3;
result[32] = x3;
result[33] = x3;
result[34] = x3;
}
        
static double coder_d2gdt2(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 40.14405;
    double x2 = sqrt(x1);
    double x3 = 0.15782994586457105*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 120.43215000000001*T - 903.24112500000001;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 40.14405);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P) + n5*(*endmember[4].d2mu0dT2)(T, P) + ((T >= 7.5) ? (
   0
)
: (
   -1.0/3.0*n2*(1021.0639949058315*x10*x6 - 161.15447504025005*x10*x8 + 322.3089500805001*((x4)*(x4))*x7*x9 + 1.1304445365082678*x5*x6 - 0.17841800000000005*x5*x8 + 101.74000828574411*x4/x2)
));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 40.14405;
    double x2 = sqrt(x1);
    double x3 = 0.15782994586457105*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 120.43215000000001*T - 903.24112500000001;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 40.14405);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = ((T >= 7.5) ? (
   0
)
: (
   -340.35466496861045*x10*x6 + 53.718158346750016*x10*x8 - 107.43631669350003*((x4)*(x4))*x7*x9 - 0.3768148455027559*x5*x6 + 0.059472666666666681*x5*x8 - 33.913336095248034*x4/x2
)) + (*endmember[1].d2mu0dT2)(T, P);
result[2] = (*endmember[2].d2mu0dT2)(T, P);
result[3] = (*endmember[3].d2mu0dT2)(T, P);
result[4] = (*endmember[4].d2mu0dT2)(T, P);
}
        
static void coder_d4gdn2dt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d5gdn3dt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d2gdtdp(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d2mu0dTdP)(T, P) + n2*(*endmember[1].d2mu0dTdP)(T, P) + n3*(*endmember[2].d2mu0dTdP)(T, P) + n4*(*endmember[3].d2mu0dTdP)(T, P) + n5*(*endmember[4].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d2mu0dTdP)(T, P);
result[1] = (*endmember[1].d2mu0dTdP)(T, P);
result[2] = (*endmember[2].d2mu0dTdP)(T, P);
result[3] = (*endmember[3].d2mu0dTdP)(T, P);
result[4] = (*endmember[4].d2mu0dTdP)(T, P);
}
        
static void coder_d4gdn2dtdp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d5gdn3dtdp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d2gdp2(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P) + n4*(*endmember[3].d2mu0dP2)(T, P) + n5*(*endmember[4].d2mu0dP2)(T, P);
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d2mu0dP2)(T, P);
result[1] = (*endmember[1].d2mu0dP2)(T, P);
result[2] = (*endmember[2].d2mu0dP2)(T, P);
result[3] = (*endmember[3].d2mu0dP2)(T, P);
result[4] = (*endmember[4].d2mu0dP2)(T, P);
}
        
static void coder_d4gdn2dp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d5gdn3dp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d3gdt3(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 40.14405;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.15782994586457105*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = x2*x4;
    double x6 = x0 - 40.14405;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 120.43215000000001*T - 903.24112500000001;
    double x10 = x4/pow(x1, 5.0/2.0);
    double x11 = pow(x6, -2);
    double x12 = x11*x8;
    double x13 = x2*0;
    double x14 = fmin(4, x3);
    double x15 = ((x14)*(x14));

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P) + n5*(*endmember[4].d3mu0dT3)(T, P) + ((T >= 7.5) ? (
   0
)
: (
   -1.0/3.0*n2*(8197.9288129398883*x10*x15 + 9.0761243991629463*x10*x9 - 2587.7573214958202*x11*x14*((x4)*(x4)) + 1293.8786607479101*x12*x15 + 1.4324842225800003*x12*x9 - 68.07093299372211*x13*x15 - 0.075362969100551208*x13*x9 - 408.42559796233263*x14*x5*x8 + 136.14186598744422*x2*((x4)*(x4)*(x4)) + 408.42559796233257*x5 - 64.461790016100025*x8/x6)
));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 5.3525400000000003*T;
    double x1 = -x0 + 40.14405;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.15782994586457105*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 136.14186598744419*x2*x4;
    double x6 = x0 - 40.14405;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 120.43215000000001*T - 903.24112500000001;
    double x10 = x4/pow(x1, 5.0/2.0);
    double x11 = pow(x6, -2);
    double x12 = x11*x8;
    double x13 = x2*0;
    double x14 = fmin(4, x3);
    double x15 = ((x14)*(x14));

result[0] = (*endmember[0].d3mu0dT3)(T, P);
result[1] = ((T >= 7.5) ? (
   0
)
: (
   -2732.6429376466294*x10*x15 - 3.0253747997209821*x10*x9 + 862.58577383194006*x11*x14*((x4)*(x4)) - 431.29288691597003*x12*x15 - 0.4774947408600001*x12*x9 + 22.690310997907368*x13*x15 + 0.025120989700183734*x13*x9 + x14*x5*x8 - 45.380621995814735*x2*((x4)*(x4)*(x4)) - x5 + 21.487263338700007*x8/x6
)) + (*endmember[1].d3mu0dT3)(T, P);
result[2] = (*endmember[2].d3mu0dT3)(T, P);
result[3] = (*endmember[3].d3mu0dT3)(T, P);
result[4] = (*endmember[4].d3mu0dT3)(T, P);
}
        
static void coder_d5gdn2dt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d6gdn3dt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d3gdt2dp(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P) + n4*(*endmember[3].d3mu0dT2dP)(T, P) + n5*(*endmember[4].d3mu0dT2dP)(T, P);
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d3mu0dT2dP)(T, P);
result[1] = (*endmember[1].d3mu0dT2dP)(T, P);
result[2] = (*endmember[2].d3mu0dT2dP)(T, P);
result[3] = (*endmember[3].d3mu0dT2dP)(T, P);
result[4] = (*endmember[4].d3mu0dT2dP)(T, P);
}
        
static void coder_d5gdn2dt2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d6gdn3dt2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d3gdtdp2(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d3mu0dTdP2)(T, P) + n2*(*endmember[1].d3mu0dTdP2)(T, P) + n3*(*endmember[2].d3mu0dTdP2)(T, P) + n4*(*endmember[3].d3mu0dTdP2)(T, P) + n5*(*endmember[4].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d3mu0dTdP2)(T, P);
result[1] = (*endmember[1].d3mu0dTdP2)(T, P);
result[2] = (*endmember[2].d3mu0dTdP2)(T, P);
result[3] = (*endmember[3].d3mu0dTdP2)(T, P);
result[4] = (*endmember[4].d3mu0dTdP2)(T, P);
}
        
static void coder_d5gdn2dtdp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d6gdn3dtdp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_d3gdp3(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P) + n4*(*endmember[3].d3mu0dP3)(T, P) + n5*(*endmember[4].d3mu0dP3)(T, P);
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d3mu0dP3)(T, P);
result[1] = (*endmember[1].d3mu0dP3)(T, P);
result[2] = (*endmember[2].d3mu0dP3)(T, P);
result[3] = (*endmember[3].d3mu0dP3)(T, P);
result[4] = (*endmember[4].d3mu0dP3)(T, P);
}
        
static void coder_d5gdn2dp3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
}
        
static void coder_d6gdn3dp3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


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
result[20] = 0;
result[21] = 0;
result[22] = 0;
result[23] = 0;
result[24] = 0;
result[25] = 0;
result[26] = 0;
result[27] = 0;
result[28] = 0;
result[29] = 0;
result[30] = 0;
result[31] = 0;
result[32] = 0;
result[33] = 0;
result[34] = 0;
}
        
static double coder_s(double T, double P, double n[5]) {
    double result = -coder_dgdt(T, P, n);
    return result;
}

static double coder_v(double T, double P, double n[5]) {
    double result = coder_dgdp(T, P, n);
    return result;
}

static double coder_cv(double T, double P, double n[5]) {
    double result = -T*coder_d2gdt2(T, P, n);
    double dvdt = coder_d2gdtdp(T, P, n);
    double dvdp = coder_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}

static double coder_cp(double T, double P, double n[5]) {
    double result = -T*coder_d2gdt2(T, P, n);
    return result;
}

static double coder_dcpdt(double T, double P, double n[5]) {
    double result = -T*coder_d3gdt3(T, P, n) - coder_d2gdt2(T, P, n);
    return result;
}

static double coder_alpha(double T, double P, double n[5]) {
    double result = coder_d2gdtdp(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_beta(double T, double P, double n[5]) {
    double result = -coder_d2gdp2(T, P, n)/coder_dgdp(T, P, n);
    return result;
}

static double coder_K(double T, double P, double n[5]) {
    double result = -coder_dgdp(T, P, n)/coder_d2gdp2(T, P, n);
    return result;
}

static double coder_Kp(double T, double P, double n[5]) {
    double result = coder_dgdp(T, P, n);
    result *= coder_d3gdp3(T, P, n);
    result /= pow(coder_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}

