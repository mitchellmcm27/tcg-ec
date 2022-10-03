#include <math.h>


static double coder_g(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = n4 + n5;
    double x1 = n1 + n2 + n3;
    double x2 = x0 + x1;
    double x3 = 1.0/x2;
    double x4 = log(n5*x3);
    double x5 = n1 + n4;
    double x6 = n5 + x1;

result = 1.0*x3*(0.0625*n1*(240000.0*n3 + 170400.0*n4) + 0.0625*n3*(240000.0*n1 + 464000.0*n4) + 0.0625*n4*(170400.0*n1 + 464000.0*n3) + 1.0*x2*(8.3144626181532395*T*(3.0*n2*log(n2*x3) + 3.0*n3*log(n3*x3) + 1.0*n4*log(n4*x3) + 0.99999999900000003*n5*(x4 - 1.0986122896681101) + 1.9999999980000001*n5*(x4 - 0.40546510910816402) + 1.0*x0*log(x0*x3) + 1.0*x1*log(x1*x3) + 3.0*x5*log(x3*x5) + 1.0*x6*log(x3*x6)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P) + n5*(*endmember[4].mu0)(T, P)));
    return result;
}
        
static void coder_dgdn(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = n2 + n3;
    double x2 = n5 + x0 + x1;
    double x3 = pow(x2, -2);
    double x4 = (*endmember[0].mu0)(T, P);
    double x5 = (*endmember[1].mu0)(T, P);
    double x6 = (*endmember[2].mu0)(T, P);
    double x7 = (*endmember[3].mu0)(T, P);
    double x8 = (*endmember[4].mu0)(T, P);
    double x9 = 1.0/x2;
    double x10 = n2*x9;
    double x11 = 3.0*log(x10);
    double x12 = n3*x9;
    double x13 = 3.0*log(x12);
    double x14 = log(n4*x9);
    double x15 = 1.0*n4;
    double x16 = n5*x9;
    double x17 = log(x16);
    double x18 = 3.0*log(x0*x9);
    double x19 = n4 + n5;
    double x20 = 1.0*log(x19*x9);
    double x21 = n1 + x1;
    double x22 = 1.0*log(x21*x9);
    double x23 = n5 + x21;
    double x24 = 1.0*log(x23*x9);
    double x25 = 8.3144626181532395*T;
    double x26 = x25*(n2*x11 + n3*x13 + 0.99999999900000003*n5*(x17 - 1.0986122896681101) + 1.9999999980000001*n5*(x17 - 0.40546510910816402) + x0*x18 + x14*x15 + x19*x20 + x21*x22 + x23*x24);
    double x27 = 1.0*x2;
    double x28 = -1.0*x3*(0.0625*n1*(240000.0*n3 + 170400.0*n4) + 0.0625*n3*(240000.0*n1 + 464000.0*n4) + 0.0625*n4*(170400.0*n1 + 464000.0*n3) + x27*(n1*x4 + n2*x5 + n3*x6 + n4*x7 + n5*x8 + x26));
    double x29 = 1.0*n5;
    double x30 = x15 + x29;
    double x31 = 1.0*n1;
    double x32 = 1.0*n2;
    double x33 = 1.0*n3;
    double x34 = x31 + x32 + x33;
    double x35 = x30 + x34;
    double x36 = -3.0*x10;
    double x37 = -3.0*x12;
    double x38 = -2.9999999970000002*x16;
    double x39 = 3.0*n1 + 3.0*n4;
    double x40 = x18 + x36 + x37 + x38 + x2*x39*(-x0*x3 + x9)/x0;
    double x41 = x29 + x34;
    double x42 = -x15*x9 + x2*x41*(-x23*x3 + x9)/x23 + x24;
    double x43 = x2*x34*(-x21*x3 + x9)/x21 + x22 - x30*x9 + x42;
    double x44 = x15*x7 + x26 + x29*x8 + x31*x4 + x32*x5 + x33*x6;
    double x45 = 1.0*x9;
    double x46 = 3.0*x2;
    double x47 = -x39*x9;
    double x48 = x37 + x47;
    double x49 = x38 + x43;
    double x50 = x20 - x34*x9 + x2*x30*(-x19*x3 + x9)/x19;

result[0] = x28 + x45*(30000.0*n3 + 21300.0*n4 + x35*(x25*(x40 + x43) + x4) + x44);
result[1] = x28 + x45*(x35*(x25*(x11 + x46*(-n2*x3 + x9) + x48 + x49) + x5) + x44);
result[2] = x28 + x45*(30000.0*n1 + 58000.0*n4 + x35*(x25*(x13 + x36 + x46*(-n3*x3 + x9) + x47 + x49) + x6) + x44);
result[3] = x28 + x45*(21300.0*n1 + 58000.0*n3 + x35*(x25*(1.0*x14 + x27*(-n4*x3 + x9) + x40 - x41*x9 + x50) + x7) + x44);
result[4] = x28 + x45*(x35*(x25*(2.9999999970000002*x17 + 2.9999999970000002*x2*(-n5*x3 + x9) + x36 + x42 + x48 + x50 - 1.9095425059748956) + x8) + x44);
}
        
static void coder_d2gdn2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = (*endmember[1].mu0)(T, P);
    double x2 = (*endmember[2].mu0)(T, P);
    double x3 = (*endmember[3].mu0)(T, P);
    double x4 = (*endmember[4].mu0)(T, P);
    double x5 = n1 + n4;
    double x6 = n2 + n3;
    double x7 = n5 + x5 + x6;
    double x8 = 1.0/x7;
    double x9 = n2*x8;
    double x10 = 3.0*log(x9);
    double x11 = n3*x8;
    double x12 = 3.0*log(x11);
    double x13 = log(n4*x8);
    double x14 = 1.0*n4;
    double x15 = n5*x8;
    double x16 = log(x15);
    double x17 = 3.0*log(x5*x8);
    double x18 = n4 + n5;
    double x19 = 1.0*log(x18*x8);
    double x20 = n1 + x6;
    double x21 = 1.0*log(x20*x8);
    double x22 = n5 + x20;
    double x23 = 1.0*log(x22*x8);
    double x24 = 8.3144626181532395*T;
    double x25 = x24*(n2*x10 + n3*x12 + 0.99999999900000003*n5*(x16 - 1.0986122896681101) + 1.9999999980000001*n5*(x16 - 0.40546510910816402) + x13*x14 + x17*x5 + x18*x19 + x20*x21 + x22*x23);
    double x26 = 1.0*x7;
    double x27 = pow(x7, -3);
    double x28 = 2.0*x27;
    double x29 = x28*(0.0625*n1*(240000.0*n3 + 170400.0*n4) + 0.0625*n3*(240000.0*n1 + 464000.0*n4) + 0.0625*n4*(170400.0*n1 + 464000.0*n3) + x26*(n1*x0 + n2*x1 + n3*x2 + n4*x3 + n5*x4 + x25));
    double x30 = 1.0*n5;
    double x31 = x14 + x30;
    double x32 = 1.0*n1;
    double x33 = 1.0*n2;
    double x34 = 1.0*n3;
    double x35 = x32 + x33 + x34;
    double x36 = x31 + x35;
    double x37 = -3.0*x9;
    double x38 = -3.0*x11;
    double x39 = -2.9999999970000002*x15;
    double x40 = 3.0*n1;
    double x41 = 3.0*n4;
    double x42 = x40 + x41;
    double x43 = 1.0/x5;
    double x44 = pow(x7, -2);
    double x45 = -x44*x5 + x8;
    double x46 = x43*x45;
    double x47 = x42*x46;
    double x48 = x17 + x37 + x38 + x39 + x47*x7;
    double x49 = 1.0/x20;
    double x50 = -x20*x44 + x8;
    double x51 = x49*x50;
    double x52 = x35*x51;
    double x53 = x30 + x35;
    double x54 = 1.0/x22;
    double x55 = -x22*x44 + x8;
    double x56 = x54*x55;
    double x57 = x53*x56;
    double x58 = -x14*x8 + x23 + x57*x7;
    double x59 = x21 - x31*x8 + x52*x7 + x58;
    double x60 = T*(x48 + x59);
    double x61 = 8.3144626181532395*x60;
    double x62 = x0*x32 + x1*x33 + x14*x3 + x2*x34 + x25 + x30*x4;
    double x63 = 30000.0*n3 + 21300.0*n4 + x36*(x0 + x61) + x62;
    double x64 = 2.0*x44;
    double x65 = n3*x44;
    double x66 = 3.0*x65;
    double x67 = n5*x44;
    double x68 = 2.9999999970000002*x67;
    double x69 = x66 + x68;
    double x70 = x14*x44;
    double x71 = x57 + x70;
    double x72 = x52 + x71;
    double x73 = x69 + x72;
    double x74 = n2*x44;
    double x75 = 3.0*x74;
    double x76 = x47 + x75;
    double x77 = x73 + x76;
    double x78 = -2*x44;
    double x79 = 2*x27;
    double x80 = x42*x7;
    double x81 = x43*x80;
    double x82 = -x45*x80/((x5)*(x5)) + 6.0*x46*x7 + x81*(x5*x79 + x78);
    double x83 = x77 + x82;
    double x84 = -x30;
    double x85 = 2.0*x7;
    double x86 = x35*x7;
    double x87 = x49*x86;
    double x88 = x53*x7;
    double x89 = x54*x88;
    double x90 = x56*x85 + x89*(x22*x79 + x78) - x55*x88/((x22)*(x22));
    double x91 = -x44*(-x14 + x84) + x51*x85 + x87*(x20*x79 + x78) + x90 - x50*x86/((x20)*(x20));
    double x92 = x24*x36;
    double x93 = 1.0*x8;
    double x94 = 3.0*x7;
    double x95 = x94*(-x74 + x8);
    double x96 = -x42*x8;
    double x97 = x38 + x96;
    double x98 = x39 + x59;
    double x99 = x10 + x95 + x97 + x98;
    double x100 = x24*x99;
    double x101 = 1.0*x1 + x100;
    double x102 = -x44;
    double x103 = -n1;
    double x104 = x81*(x102 - x79*(-n4 + x103));
    double x105 = x104 + x77;
    double x106 = 1.0*x0 + x61;
    double x107 = x106 + x92*(x105 - 6.0*x8 + x91);
    double x108 = x36*(x1 + x100) + x62;
    double x109 = 1.0*x44;
    double x110 = -x108*x109;
    double x111 = -x109*x63 + x29;
    double x112 = x94*(-x65 + x8);
    double x113 = x112 + x12 + x37 + x96 + x98;
    double x114 = x113*x24;
    double x115 = x114 + 1.0*x2;
    double x116 = 30000.0*n1 + 58000.0*n4 + x36*(x114 + x2) + x62;
    double x117 = -x109*x116;
    double x118 = -n2 - n3 + x103;
    double x119 = x31*x44 + x87*(x102 - x118*x79);
    double x120 = x119 + x89*(x102 - x79*(-n5 + x118));
    double x121 = x26*(-n4*x44 + x8);
    double x122 = 1.0/x18;
    double x123 = -x18*x44 + x8;
    double x124 = x122*x123;
    double x125 = x124*x31;
    double x126 = x125*x7 + x19 - x35*x8;
    double x127 = x121 + x126 + 1.0*x13 + x48 - x53*x8;
    double x128 = x127*x24;
    double x129 = x128 + 1.0*x3;
    double x130 = 21300.0*n1 + 58000.0*n3 + x36*(x128 + x3) + x62;
    double x131 = -x109*x130;
    double x132 = x119 + x90;
    double x133 = 2.9999999970000002*x7*(-x67 + x8);
    double x134 = x126 + x133 + 2.9999999970000002*x16 + x37 + x58 + x97 - 1.9095425059748956;
    double x135 = x134*x24;
    double x136 = x135 + 1.0*x4;
    double x137 = x36*(x135 + x4) + x62;
    double x138 = -x109*x137;
    double x139 = 16.628925236306479*T;
    double x140 = -6.0*x44;
    double x141 = 6.0*x27;
    double x142 = n2*x141;
    double x143 = x73 - x75;
    double x144 = -x44*(-x40 - x41);
    double x145 = 3.0*x8;
    double x146 = x144 + x145 + x91;
    double x147 = -3.0*x44;
    double x148 = x143 + x7*(x142 + x147);
    double x149 = x144 + x148;
    double x150 = x110 + x29;
    double x151 = x120 + x42*x44 - 7.0*x8;
    double x152 = x132 - 4.9999999969999998*x8;
    double x153 = n3*x141;
    double x154 = -x66 + x68 + x72 + x75;
    double x155 = x154 + x7*(x147 + x153);
    double x156 = x117 + x29;
    double x157 = n4*x28;
    double x158 = -x32 - x33 - x34;
    double x159 = x31*x7;
    double x160 = x122*x159*(x18*x79 + x78) - x123*x159/((x18)*(x18)) + x124*x85 + x125 - x158*x44;
    double x161 = x160 + x69 - x70 + x76;

result[0] = x29 - x63*x64 + x93*(2.0*x0 + 16.628925236306479*x60 + x92*(x83 + x91));
result[1] = x110 + x111 + x93*(x101 + x107);
result[2] = x111 + x117 + x93*(x107 + x115 + 30000.0);
result[3] = x111 + x131 + x93*(x106 + x129 + x92*(x120 - 4.0*x8 + x83) + 21300.0);
result[4] = x111 + x138 + x93*(x106 + x136 + x92*(x105 + x132 - 7.9999999969999998*x8));
result[5] = -x108*x64 + x29 + x93*(2.0*x1 + x139*x99 + x92*(x143 + x146 + x7*(x140 + x142) + x95/n2));
result[6] = x117 + x150 + x93*(x101 + x115 + x92*(-x145 + x149 + x91));
result[7] = x131 + x150 + x93*(x101 + x129 + x92*(x148 + x151));
result[8] = x138 + x150 + x93*(x101 + x136 + x92*(x149 + x152));
result[9] = -x116*x64 + x29 + x93*(x113*x139 + 2.0*x2 + x92*(x146 + x154 + x7*(x140 + x153) + x112/n3));
result[10] = x131 + x156 + x93*(x115 + x129 + x92*(x151 + x155) + 58000.0);
result[11] = x138 + x156 + x93*(x115 + x136 + x92*(x144 + x152 + x155));
result[12] = -x130*x64 + x29 + x93*(x127*x139 + 2.0*x3 + x92*(x161 - x44*(x158 + x84) + x7*(x157 - x64) + x82 + x93 + x121/n4));
result[13] = x131 + x138 + x29 + x93*(x129 + x136 + x92*(x104 + x161 + x44*x53 + x7*(-x109 + x157) - 6.9999999969999998*x8));
result[14] = -x137*x64 + x29 + x93*(x134*x139 + 2.0*x4 + x92*(x144 + x160 + x66 - x68 + x7*(5.9999999940000004*n5*x27 - 5.9999999940000004*x44) + x71 + x75 + 2.9999999970000002*x8 + x90 + x133/n5));
}
        
static void coder_d3gdn3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = (*endmember[1].mu0)(T, P);
    double x2 = (*endmember[2].mu0)(T, P);
    double x3 = (*endmember[3].mu0)(T, P);
    double x4 = (*endmember[4].mu0)(T, P);
    double x5 = n1 + n4;
    double x6 = n2 + n3;
    double x7 = n5 + x5 + x6;
    double x8 = 1.0/x7;
    double x9 = n2*x8;
    double x10 = 3.0*log(x9);
    double x11 = n3*x8;
    double x12 = 3.0*log(x11);
    double x13 = log(n4*x8);
    double x14 = 1.0*n4;
    double x15 = n5*x8;
    double x16 = log(x15);
    double x17 = 3.0*log(x5*x8);
    double x18 = n4 + n5;
    double x19 = 1.0*log(x18*x8);
    double x20 = n1 + x6;
    double x21 = 1.0*log(x20*x8);
    double x22 = n5 + x20;
    double x23 = 1.0*log(x22*x8);
    double x24 = 8.3144626181532395*T;
    double x25 = x24*(n2*x10 + n3*x12 + 0.99999999900000003*n5*(x16 - 1.0986122896681101) + 1.9999999980000001*n5*(x16 - 0.40546510910816402) + x13*x14 + x17*x5 + x18*x19 + x20*x21 + x22*x23);
    double x26 = 1.0*x7;
    double x27 = pow(x7, -4);
    double x28 = 6.0*x27;
    double x29 = -x28*(0.0625*n1*(240000.0*n3 + 170400.0*n4) + 0.0625*n3*(240000.0*n1 + 464000.0*n4) + 0.0625*n4*(170400.0*n1 + 464000.0*n3) + x26*(n1*x0 + n2*x1 + n3*x2 + n4*x3 + n5*x4 + x25));
    double x30 = 1.0*n5;
    double x31 = x14 + x30;
    double x32 = 1.0*n1;
    double x33 = 1.0*n2;
    double x34 = 1.0*n3;
    double x35 = x32 + x33 + x34;
    double x36 = x31 + x35;
    double x37 = -3.0*x9;
    double x38 = -3.0*x11;
    double x39 = -2.9999999970000002*x15;
    double x40 = 3.0*n1;
    double x41 = 3.0*n4;
    double x42 = x40 + x41;
    double x43 = 1.0/x5;
    double x44 = pow(x7, -2);
    double x45 = -x44*x5 + x8;
    double x46 = x43*x45;
    double x47 = x42*x46;
    double x48 = x17 + x37 + x38 + x39 + x47*x7;
    double x49 = 1.0/x20;
    double x50 = -x20*x44 + x8;
    double x51 = x49*x50;
    double x52 = x35*x51;
    double x53 = x30 + x35;
    double x54 = 1.0/x22;
    double x55 = -x22*x44 + x8;
    double x56 = x54*x55;
    double x57 = x53*x56;
    double x58 = -x14*x8 + x23 + x57*x7;
    double x59 = x21 - x31*x8 + x52*x7 + x58;
    double x60 = x48 + x59;
    double x61 = x24*x60;
    double x62 = x0*x32 + x1*x33 + x14*x3 + x2*x34 + x25 + x30*x4;
    double x63 = 30000.0*n3 + 21300.0*n4 + x36*(x0 + x61) + x62;
    double x64 = pow(x7, -3);
    double x65 = 6.0*x64;
    double x66 = 16.628925236306479*T;
    double x67 = 3.0*x44;
    double x68 = n3*x67;
    double x69 = n5*x44;
    double x70 = 2.9999999970000002*x69;
    double x71 = x68 + x70;
    double x72 = x14*x44;
    double x73 = x57 + x72;
    double x74 = x52 + x73;
    double x75 = x71 + x74;
    double x76 = n2*x67;
    double x77 = x47 + x76;
    double x78 = x75 + x77;
    double x79 = 6.0*x46;
    double x80 = -2*x44;
    double x81 = 2*x64;
    double x82 = x5*x81 + x80;
    double x83 = x42*x43;
    double x84 = x82*x83;
    double x85 = pow(x5, -2);
    double x86 = x45*x85;
    double x87 = x42*x86;
    double x88 = x7*x79 + x7*x84 - x7*x87;
    double x89 = x78 + x88;
    double x90 = -x30;
    double x91 = -x14 + x90;
    double x92 = 2.0*x51;
    double x93 = x20*x81 + x80;
    double x94 = x35*x49;
    double x95 = x93*x94;
    double x96 = pow(x20, -2);
    double x97 = x50*x96;
    double x98 = x35*x97;
    double x99 = 2.0*x56;
    double x100 = x22*x81 + x80;
    double x101 = x53*x54;
    double x102 = x100*x101;
    double x103 = pow(x22, -2);
    double x104 = x103*x55;
    double x105 = x104*x53;
    double x106 = x102*x7 - x105*x7 + x7*x99;
    double x107 = x106 - x44*x91 + x7*x92 + x7*x95 - x7*x98;
    double x108 = T*(x107 + x89);
    double x109 = 8.3144626181532395*x108;
    double x110 = 2.0*x0 + x109*x36 + x60*x66;
    double x111 = 9.0*x7;
    double x112 = 6*x64;
    double x113 = 6*x27;
    double x114 = x7*x83;
    double x115 = x42*x85;
    double x116 = 2*x7;
    double x117 = x111*x43*x82 - x111*x86 + x114*(x112 - x113*x5) - x115*x116*x82 + x116*x42*x45/((x5)*(x5)*(x5)) + 9.0*x46 + 2*x84 - 2*x87;
    double x118 = 3.0*x7;
    double x119 = x101*x7;
    double x120 = x53*x7;
    double x121 = x103*x120;
    double x122 = x100*x118*x54 - 2*x100*x121 + 2*x102 - x104*x118 - 2*x105 + x119*(x112 - x113*x22) + 2*x120*x55/((x22)*(x22)*(x22)) + 3.0*x56;
    double x123 = -x31*x81;
    double x124 = n2*x65;
    double x125 = -x124;
    double x126 = 2.0*x64;
    double x127 = n4*x126;
    double x128 = -x127;
    double x129 = n3*x65;
    double x130 = -x129;
    double x131 = n5*x64;
    double x132 = 5.9999999940000004*x131;
    double x133 = -x132;
    double x134 = x130 + x133;
    double x135 = x128 + x134;
    double x136 = x125 + x135;
    double x137 = x123 + x136;
    double x138 = x7*x94;
    double x139 = x35*x7;
    double x140 = x139*x96;
    double x141 = x118*x49*x93 - x118*x97 + x138*(x112 - x113*x20) + 2*x139*x50/((x20)*(x20)*(x20)) - 2*x140*x93 + 3.0*x51 + 2*x95 - 2*x98;
    double x142 = x137 + x141;
    double x143 = x24*x36;
    double x144 = 1.0*x8;
    double x145 = -n2*x44 + x8;
    double x146 = x118*x145;
    double x147 = -x42*x8;
    double x148 = x147 + x38;
    double x149 = x39 + x59;
    double x150 = x10 + x146 + x148 + x149;
    double x151 = x150*x24;
    double x152 = x36*(x1 + x151) + x62;
    double x153 = x126*x152;
    double x154 = 1.0*x1 + x151;
    double x155 = -x44;
    double x156 = -n1;
    double x157 = -n4 + x156;
    double x158 = x155 - x157*x81;
    double x159 = x158*x83;
    double x160 = x159*x7;
    double x161 = x160 + x78;
    double x162 = x107 + x161 - 6.0*x8;
    double x163 = x162*x24;
    double x164 = 1.0*x0 + x61;
    double x165 = x163*x36 + x164;
    double x166 = x154 + x165;
    double x167 = 2.0*x44;
    double x168 = -x166*x167 + x29;
    double x169 = x162*x66;
    double x170 = x158*x7;
    double x171 = 4*x64;
    double x172 = 2*n1;
    double x173 = 2*n4;
    double x174 = 3*x27;
    double x175 = -x174*(x172 + x173);
    double x176 = -x115*x170 + x159 + x84 - x87;
    double x177 = x114*(x171 + x175) + 6.0*x170*x43 + x176 + x79;
    double x178 = x122 + x67;
    double x179 = 4.0*x64;
    double x180 = 1.0*x44;
    double x181 = -x110*x180 + x179*x63;
    double x182 = x144*(x109 + x143*(x142 + x177 + x178) + x169) + x181;
    double x183 = -n3*x44 + x8;
    double x184 = x118*x183;
    double x185 = x12 + x147 + x149 + x184 + x37;
    double x186 = x185*x24;
    double x187 = 30000.0*n1 + 58000.0*n4 + x36*(x186 + x2) + x62;
    double x188 = x126*x187;
    double x189 = x186 + 1.0*x2;
    double x190 = x165 + x189 + 30000.0;
    double x191 = -x167*x190 + x29;
    double x192 = -n2 - n3 + x156;
    double x193 = -n5 + x192;
    double x194 = x155 - x193*x81;
    double x195 = x101*x194;
    double x196 = x155 - x192*x81;
    double x197 = x196*x94;
    double x198 = x197*x7 + x31*x44;
    double x199 = x195*x7 + x198;
    double x200 = x199 - 4.0*x8 + x89;
    double x201 = x200*x66;
    double x202 = 2.0*x7;
    double x203 = x194*x54;
    double x204 = 2*n5;
    double x205 = 2*n2;
    double x206 = 2*n3;
    double x207 = x172 + x205 + x206;
    double x208 = -x174*(x204 + x207);
    double x209 = x102 - x105 - x121*x194 + x195;
    double x210 = x119*(x171 + x208) + x202*x203 + x209 + x99;
    double x211 = -x174*x207;
    double x212 = x138*(x171 + x211) - x140*x196 + x196*x202*x49 + x197 + x81*x91 + x92 + x95 - x98;
    double x213 = x136 + x212;
    double x214 = x210 + x213;
    double x215 = -n4*x44 + x8;
    double x216 = x215*x26;
    double x217 = 1.0/x18;
    double x218 = -x18*x44 + x8;
    double x219 = x217*x218;
    double x220 = x219*x31;
    double x221 = x19 + x220*x7 - x35*x8;
    double x222 = 1.0*x13 + x216 + x221 + x48 - x53*x8;
    double x223 = x222*x24;
    double x224 = 21300.0*n1 + 58000.0*n3 + x36*(x223 + x3) + x62;
    double x225 = x126*x224;
    double x226 = x200*x24;
    double x227 = x223 + 1.0*x3;
    double x228 = x164 + x226*x36 + x227 + 21300.0;
    double x229 = -x167*x228;
    double x230 = x181 + x29;
    double x231 = x106 + x198;
    double x232 = x161 + x231 - 7.9999999969999998*x8;
    double x233 = x232*x66;
    double x234 = 3.9999999970000002*x44;
    double x235 = -2.9999999970000002*x69 + 2.9999999970000002*x8;
    double x236 = x235*x7;
    double x237 = x148 + 2.9999999970000002*x16 + x221 + x236 + x37 + x58 - 1.9095425059748956;
    double x238 = x237*x24;
    double x239 = x36*(x238 + x4) + x62;
    double x240 = x126*x239;
    double x241 = x232*x24;
    double x242 = x238 + 1.0*x4;
    double x243 = x164 + x241*x36 + x242;
    double x244 = -x167*x243;
    double x245 = 6.0*x44;
    double x246 = -x245;
    double x247 = 1.0/n2;
    double x248 = x75 - x76;
    double x249 = -x40 - x41;
    double x250 = -x249*x44;
    double x251 = 3.0*x8;
    double x252 = x107 + x250 + x251;
    double x253 = x146*x247 + x248 + x252 + x7*(x124 + x246);
    double x254 = x24*x253;
    double x255 = 9.0*x44;
    double x256 = x114*(x175 + x81) + 2*x159;
    double x257 = x122 + x256;
    double x258 = x143*(x142 + x255 + x257) + x169;
    double x259 = x126*x63;
    double x260 = 2.0*x1 + x150*x66 + x254*x36;
    double x261 = x152*x179 - x180*x260;
    double x262 = -x67;
    double x263 = x248 + x7*(x124 + x262);
    double x264 = x250 + x263;
    double x265 = x107 - x251 + x264;
    double x266 = x24*x265;
    double x267 = -x180*x190;
    double x268 = x154 + x189 + x266*x36;
    double x269 = -x180*x268 + x188;
    double x270 = x259 + x29;
    double x271 = x153 - x166*x180 + x270;
    double x272 = x199 + x42*x44 - 7.0*x8;
    double x273 = x263 + x272;
    double x274 = x24*x273;
    double x275 = x143*(x114*(x113*x157 + x171) + x118*x158*x43 + x176 + x214 + 8.0*x44 + 3.0*x46) + x163 + x226;
    double x276 = x154 + x227 + x274*x36;
    double x277 = -x180*x276;
    double x278 = -x180*x228 + x225;
    double x279 = x231 - 4.9999999969999998*x8;
    double x280 = x264 + x279;
    double x281 = x24*x280;
    double x282 = x143*(x213 + x257 + 9.9999999969999998*x44) + x163 + x241;
    double x283 = x154 + x242 + x281*x36;
    double x284 = -x180*x283;
    double x285 = -x180*x243 + x240;
    double x286 = 1.0/n3;
    double x287 = -x68 + x70 + x74 + x76;
    double x288 = x184*x286 + x252 + x287 + x7*(x129 + x246);
    double x289 = x24*x288;
    double x290 = x185*x66 + 2.0*x2 + x289*x36;
    double x291 = x179*x187 - x180*x290;
    double x292 = x287 + x7*(x129 + x262);
    double x293 = x272 + x292;
    double x294 = x24*x293;
    double x295 = x189 + x227 + x294*x36 + 58000.0;
    double x296 = -x180*x295;
    double x297 = x270 + x278;
    double x298 = x188 + x267;
    double x299 = x250 + x279 + x292;
    double x300 = x24*x299;
    double x301 = x189 + x242 + x300*x36;
    double x302 = -x180*x301;
    double x303 = -x32 - x33 - x34;
    double x304 = x303 + x90;
    double x305 = 1.0/n4;
    double x306 = x31*x7;
    double x307 = x18*x81 + x80;
    double x308 = x217*x307;
    double x309 = pow(x18, -2);
    double x310 = x218*x309;
    double x311 = x202*x219 + x220 - x303*x44 + x306*x308 - x306*x310;
    double x312 = x311 + x71 - x72 + x77;
    double x313 = x144 + x216*x305 - x304*x44 + x312 + x7*(x127 - x167) + x88;
    double x314 = x24*x313;
    double x315 = x119*(x208 + x81) + 2*x195;
    double x316 = 2*x197;
    double x317 = x138*(x211 + x81);
    double x318 = x137 + x316 + x317;
    double x319 = x222*x66 + 2.0*x3 + x314*x36;
    double x320 = x179*x224 - x180*x319;
    double x321 = x160 + x312 + x44*x53 + x7*(x127 - x180) - 6.9999999969999998*x8;
    double x322 = x24*x321;
    double x323 = x119*(x113*x193 + x171) + x203*x26 + x209 + 7.9999999969999998*x44 + 1.0*x56;
    double x324 = x227 + x242 + x322*x36;
    double x325 = -x180*x324;
    double x326 = 1.0/n5;
    double x327 = x106 + x236*x326 + x250 + x311 + x68 + x7*(x132 - 5.9999999940000004*x44) - x70 + x73 + x76 + 2.9999999970000002*x8;
    double x328 = x24*x327;
    double x329 = x237*x66 + x328*x36 + 2.0*x4;
    double x330 = x179*x239 - x180*x329;
    double x331 = 24.943387854459719*T;
    double x332 = 18.0*x64;
    double x333 = 18.0*x27;
    double x334 = -n2*x333;
    double x335 = x205*x64;
    double x336 = x118*x247;
    double x337 = 12.0*x64;
    double x338 = n2*x337 + x135;
    double x339 = 3.0*x145*x247 + x338;
    double x340 = -x42*x81;
    double x341 = x123 + x340;
    double x342 = x141 + x341;
    double x343 = x122 + x342;
    double x344 = x343 - 12.0*x44;
    double x345 = x265*x66;
    double x346 = x336*(x155 + x335) + x339 + x7*(x334 + x337);
    double x347 = -x167*x268;
    double x348 = x261 + x29;
    double x349 = x273*x66;
    double x350 = -4.0*x44;
    double x351 = x210 + x249*x81;
    double x352 = x212 + x350 + x351;
    double x353 = -x167*x276;
    double x354 = x280*x66;
    double x355 = x122 + x340;
    double x356 = x212 + x355 - 5.0000000030000002*x44;
    double x357 = -x167*x283;
    double x358 = x338 + x7*(x334 + x65);
    double x359 = x153 + x29;
    double x360 = x212 + x358;
    double x361 = x269 + x359;
    double x362 = x225 + x296;
    double x363 = x240 + x284;
    double x364 = x316 + x317 + x341;
    double x365 = x255 + x315 + x364;
    double x366 = x358 + x364;
    double x367 = x122 + 5.9999999939999995*x44;
    double x368 = -n3*x333;
    double x369 = x206*x64;
    double x370 = x118*x286;
    double x371 = x125 + x128;
    double x372 = n3*x337 + x133 + x371;
    double x373 = 3.0*x183*x286 + x372;
    double x374 = x293*x66;
    double x375 = x370*(x155 + x369) + x373 + x7*(x337 + x368);
    double x376 = -x167*x295;
    double x377 = x29 + x291;
    double x378 = x299*x66;
    double x379 = -x167*x301;
    double x380 = x372 + x7*(x368 + x65);
    double x381 = x188 + x29;
    double x382 = x364 + x380;
    double x383 = 1.0*x215*x305;
    double x384 = -n4*x28;
    double x385 = x173*x64;
    double x386 = x26*x305;
    double x387 = 2*x31;
    double x388 = 2*x306;
    double x389 = x118*x308 - x118*x310 + x217*x306*(x112 - x113*x18) + 3.0*x219 - x307*x309*x388 + x308*x387 - x310*x387 - x35*x81 + x218*x388/((x18)*(x18)*(x18));
    double x390 = n4*x179 + x125 + x134 + x389;
    double x391 = x390 - x53*x81;
    double x392 = x321*x66;
    double x393 = -x167*x324 + x29;

result[0] = -x110*x67 + x144*(24.943387854459719*x108 + x143*(x117 + x122 + x142)) + x29 + x63*x65;
result[1] = x153 + x168 + x182;
result[2] = x182 + x188 + x191;
result[3] = x144*(x109 + x143*(x117 + x167 + x214) + x201) + x225 + x229 + x230;
result[4] = x144*(x109 + x143*(x122 + x177 + x213 + x234) + x233) + x230 + x240 + x244;
result[5] = x144*(x254 + x258) + x168 + x259 + x261;
result[6] = x144*(x258 + x266) + x267 + x269 + x271;
result[7] = x144*(x274 + x275) + x271 + x277 + x278;
result[8] = x144*(x281 + x282) + x271 + x284 + x285;
result[9] = x144*(x258 + x289) + x191 + x259 + x291;
result[10] = x144*(x275 + x294) + x296 + x297 + x298;
result[11] = x144*(x282 + x300) + x270 + x285 + x298 + x302;
result[12] = x144*(x143*(x117 + x245 + x315 + x318) + x201 + x314) + x229 + x270 + x320;
result[13] = x144*(x143*(x177 + x318 + x323) + x226 + x241 + x322) + x285 + x297 + x325;
result[14] = x144*(x143*(x257 + x318 + 11.999999994*x44) + x233 + x328) + x244 + x270 + x330;
result[15] = x144*(x143*(x336*(x335 + x80) + x339 + x344 + x7*(x332 + x334) - x146/((n2)*(n2))) + x253*x331) + x152*x65 - x260*x67 + x29;
result[16] = x144*(x143*(x246 + x343 + x346) + x254 + x345) + x188 + x347 + x348;
result[17] = x144*(x143*(x346 + x352) + x254 + x349) + x225 + x348 + x353;
result[18] = x144*(x143*(x346 + x356) + x254 + x354) + x240 + x348 + x357;
result[19] = x144*(x143*(x178 + x342 + x358) + x289 + x345) + x291 + x347 + x359;
result[20] = x144*(x143*(x351 + x360 + 5.0*x44) + x266 + x274 + x294) + x277 + x361 + x362;
result[21] = x144*(x143*(x234 + x355 + x360) + x266 + x281 + x300) + x302 + x361 + x363;
result[22] = x144*(x143*(x358 + x365) + x314 + x349) + x320 + x353 + x359;
result[23] = x144*(x143*(x323 + x366) + x274 + x281 + x322) + x225 + x277 + x325 + x359 + x363;
result[24] = x144*(x143*(x366 + x367) + x328 + x354) + x330 + x357 + x359;
result[25] = x144*(x143*(x344 + x370*(x369 + x80) + x373 + x7*(x332 + x368) - x184/((n3)*(n3))) + x288*x331) + x187*x65 + x29 - x290*x67;
result[26] = x144*(x143*(x352 + x375) + x289 + x374) + x225 + x376 + x377;
result[27] = x144*(x143*(x356 + x375) + x289 + x378) + x240 + x377 + x379;
result[28] = x144*(x143*(x365 + x380) + x314 + x374) + x320 + x376 + x381;
result[29] = x144*(x143*(x323 + x382) + x294 + x300 + x322) + x240 + x302 + x325 + x362 + x381;
result[30] = x144*(x143*(x367 + x382) + x328 + x378) + x330 + x379 + x381;
result[31] = x144*(x143*(x117 + x350 + x383 + x386*(x385 + x80) + x391 + x7*(x384 + x65) - x216/((n4)*(n4))) + x313*x331) + x224*x65 + x29 - x319*x67;
result[32] = x144*(x143*(x177 + x304*x81 + x383 + x386*(x155 + x385) + x390 + 0.9999999970000002*x44 + x7*(x179 + x384)) + x314 + x392) + x240 + x320 + x393;
result[33] = x144*(x143*(x256 + x391 + 9.9999999939999995*x44 + x7*(x126 + x384)) + x328 + x392) + x225 + x330 + x393;
result[34] = x144*(x143*(x130 + 11.999999988000001*x131 + x235*x326 + 2.9999999970000002*x326*x7*(x204*x64 + x80) + x355 + x371 + x389 - 11.999999988000001*x44 + x7*(-17.999999982000002*n5*x27 + 17.999999982000002*x64) - x236/((n5)*(n5))) + x327*x331) + x239*x65 + x29 - x329*x67;
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
    double x4 = n1 + n4;
    double x5 = n5 + x1;

result = 1.0*x2*(1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5)*(n1*(*endmember[0].dmu0dT)(T, P) + 24.943387854459719*n2*log(n2*x2) + n2*(*endmember[1].dmu0dT)(T, P) + 24.943387854459719*n3*log(n3*x2) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*log(n4*x2) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626098387775*n5*(x3 - 1.0986122896681101) + 16.628925219677555*n5*(x3 - 0.40546510910816402) + n5*(*endmember[4].dmu0dT)(T, P) + 8.3144626181532395*x0*log(x0*x2) + 8.3144626181532395*x1*log(x1*x2) + 24.943387854459719*x4*log(x2*x4) + 8.3144626181532395*x5*log(x2*x5));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = n1 + n4;
    double x2 = n4 + n5;
    double x3 = n1 + n2 + n3;
    double x4 = x2 + x3;
    double x5 = 1.0/x4;
    double x6 = 24.943387854459719*log(x1*x5);
    double x7 = n5*x5;
    double x8 = -24.943387829516332*x7;
    double x9 = n2*x5;
    double x10 = -24.943387854459719*x9;
    double x11 = n3*x5;
    double x12 = -24.943387854459719*x11;
    double x13 = 24.943387854459719*n1 + 24.943387854459719*n4;
    double x14 = pow(x4, -2);
    double x15 = x10 + x12 + x6 + x8 + x13*x4*(-x1*x14 + x5)/x1;
    double x16 = 8.3144626181532395*log(x3*x5);
    double x17 = 8.3144626181532395*n4;
    double x18 = 8.3144626181532395*n5;
    double x19 = x17 + x18;
    double x20 = 8.3144626181532395*n1 + 8.3144626181532395*n2 + 8.3144626181532395*n3;
    double x21 = n5 + x3;
    double x22 = 8.3144626181532395*log(x21*x5);
    double x23 = n4*x5;
    double x24 = x18 + x20;
    double x25 = x22 - 8.3144626181532395*x23 + x24*x4*(-x14*x21 + x5)/x21;
    double x26 = x16 - x19*x5 + x20*x4*(-x14*x3 + x5)/x3 + x25;
    double x27 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x28 = 1.0*x5;
    double x29 = x27*x28;
    double x30 = (*endmember[1].dmu0dT)(T, P);
    double x31 = (*endmember[2].dmu0dT)(T, P);
    double x32 = (*endmember[3].dmu0dT)(T, P);
    double x33 = (*endmember[4].dmu0dT)(T, P);
    double x34 = 24.943387854459719*log(x9);
    double x35 = 24.943387854459719*log(x11);
    double x36 = log(x23);
    double x37 = log(x7);
    double x38 = 8.3144626181532395*log(x2*x5);
    double x39 = n1*x0 + n2*x30 + n2*x34 + n3*x31 + n3*x35 + n4*x32 + n5*x33 + 8.3144626098387775*n5*(x37 - 1.0986122896681101) + 16.628925219677555*n5*(x37 - 0.40546510910816402) + x1*x6 + x16*x3 + x17*x36 + x2*x38 + x21*x22;
    double x40 = -1.0*x14*x27*x39 + x28*x39;
    double x41 = 24.943387854459719*x4;
    double x42 = -x13*x5;
    double x43 = x12 + x42;
    double x44 = x26 + x8;
    double x45 = x19*x4*(-x14*x2 + x5)/x2 - x20*x5 + x38;

result[0] = x29*(x0 + x15 + x26) + x40;
result[1] = x29*(x30 + x34 + x41*(-n2*x14 + x5) + x43 + x44) + x40;
result[2] = x29*(x10 + x31 + x35 + x41*(-n3*x14 + x5) + x42 + x44) + x40;
result[3] = x29*(x15 - x24*x5 + x32 + 8.3144626181532395*x36 + 8.3144626181532395*x4*(-n4*x14 + x5) + x45) + x40;
result[4] = x29*(x10 + x25 + x33 + 24.943387829516332*x37 + 24.943387829516332*x4*(-n5*x14 + x5) + x43 + x45 - 15.876819783702929) + x40;
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
    double x3 = 1.0/x2;
    double x4 = (*endmember[0].dmu0dT)(T, P);
    double x5 = 24.943387854459719*log(x0*x3);
    double x6 = n5*x3;
    double x7 = -24.943387829516332*x6;
    double x8 = n2*x3;
    double x9 = -24.943387854459719*x8;
    double x10 = n3*x3;
    double x11 = -24.943387854459719*x10;
    double x12 = 24.943387854459719*n1;
    double x13 = 24.943387854459719*n4;
    double x14 = x12 + x13;
    double x15 = 1.0/x0;
    double x16 = pow(x2, -2);
    double x17 = -x0*x16 + x3;
    double x18 = x15*x17;
    double x19 = x14*x18;
    double x20 = x11 + x19*x2 + x5 + x7 + x9;
    double x21 = n1 + x1;
    double x22 = 8.3144626181532395*log(x21*x3);
    double x23 = 8.3144626181532395*n4;
    double x24 = 8.3144626181532395*n5;
    double x25 = x23 + x24;
    double x26 = 8.3144626181532395*n1;
    double x27 = 8.3144626181532395*n2;
    double x28 = 8.3144626181532395*n3;
    double x29 = x26 + x27 + x28;
    double x30 = 1.0/x21;
    double x31 = -x16*x21 + x3;
    double x32 = x30*x31;
    double x33 = x29*x32;
    double x34 = n5 + x21;
    double x35 = 8.3144626181532395*log(x3*x34);
    double x36 = x24 + x29;
    double x37 = 1.0/x34;
    double x38 = -x16*x34 + x3;
    double x39 = x37*x38;
    double x40 = x36*x39;
    double x41 = x2*x40 - x23*x3 + x35;
    double x42 = x2*x33 + x22 - x25*x3 + x41;
    double x43 = x20 + x4 + x42;
    double x44 = x3*x43;
    double x45 = n5*x16;
    double x46 = 24.943387829516332*x45;
    double x47 = n3*x16;
    double x48 = 24.943387854459719*x47;
    double x49 = x46 + x48;
    double x50 = x16*x23;
    double x51 = x40 + x50;
    double x52 = x33 + x51;
    double x53 = x49 + x52;
    double x54 = n2*x16;
    double x55 = 24.943387854459719*x54;
    double x56 = x19 + x55;
    double x57 = x53 + x56;
    double x58 = -2*x16;
    double x59 = pow(x2, -3);
    double x60 = 2*x59;
    double x61 = x14*x2;
    double x62 = x15*x61;
    double x63 = 49.886775708919437*x18*x2 + x62*(x0*x60 + x58) - x17*x61/((x0)*(x0));
    double x64 = x57 + x63;
    double x65 = -x24;
    double x66 = 16.628925236306479*x2;
    double x67 = x2*x29;
    double x68 = x30*x67;
    double x69 = x2*x36;
    double x70 = x37*x69;
    double x71 = x39*x66 + x70*(x34*x60 + x58) - x38*x69/((x34)*(x34));
    double x72 = -x16*(-x23 + x65) + x32*x66 + x68*(x21*x60 + x58) + x71 - x31*x67/((x21)*(x21));
    double x73 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x74 = 1.0*x3;
    double x75 = x73*x74;
    double x76 = x16*x73;
    double x77 = x43*x76;
    double x78 = (*endmember[1].dmu0dT)(T, P);
    double x79 = (*endmember[2].dmu0dT)(T, P);
    double x80 = (*endmember[3].dmu0dT)(T, P);
    double x81 = (*endmember[4].dmu0dT)(T, P);
    double x82 = 24.943387854459719*log(x8);
    double x83 = 24.943387854459719*log(x10);
    double x84 = log(n4*x3);
    double x85 = log(x6);
    double x86 = n4 + n5;
    double x87 = 8.3144626181532395*log(x3*x86);
    double x88 = 2.0*n1*x4 + 2.0*n2*x78 + 2.0*n2*x82 + 2.0*n3*x79 + 2.0*n3*x83 + 2.0*n4*x80 + 2.0*n5*x81 + 16.628925219677555*n5*(x85 - 1.0986122896681101) + 33.25785043935511*n5*(x85 - 0.40546510910816402) + 2.0*x0*x5 + 2.0*x21*x22 + 2.0*x23*x84 + 2.0*x34*x35 + 2.0*x86*x87;
    double x89 = -x16*x88 + x59*x73*x88;
    double x90 = 24.943387854459719*x2;
    double x91 = x90*(x3 - x54);
    double x92 = -x14*x3;
    double x93 = x11 + x92;
    double x94 = x42 + x7;
    double x95 = x78 + x82 + x91 + x93 + x94;
    double x96 = 1.0*x76;
    double x97 = x74*x95 - x95*x96;
    double x98 = -x16;
    double x99 = -n1;
    double x100 = x62*(-x60*(-n4 + x99) + x98);
    double x101 = x100 + x57;
    double x102 = 1.0*x44 - 1.0*x77 + x89;
    double x103 = x102 + x75*(x101 - 49.886775708919437*x3 + x72);
    double x104 = x90*(x3 - x47);
    double x105 = x104 + x79 + x83 + x9 + x92 + x94;
    double x106 = x105*x74 - x105*x96;
    double x107 = -n2 - n3 + x99;
    double x108 = x16*x25 + x68*(-x107*x60 + x98);
    double x109 = x108 + x70*(-x60*(-n5 + x107) + x98);
    double x110 = 8.3144626181532395*x2*(-n4*x16 + x3);
    double x111 = 1.0/x86;
    double x112 = -x16*x86 + x3;
    double x113 = x111*x112;
    double x114 = x113*x25;
    double x115 = x114*x2 - x29*x3 + x87;
    double x116 = x110 + x115 + x20 - x3*x36 + x80 + 8.3144626181532395*x84;
    double x117 = x116*x74 - x116*x96;
    double x118 = x108 + x71;
    double x119 = 24.943387829516332*x2*(x3 - x45);
    double x120 = x115 + x119 + x41 + x81 + 24.943387829516332*x85 + x9 + x93 - 15.876819783702929;
    double x121 = x120*x74 - x120*x96;
    double x122 = 2.0*x3;
    double x123 = -49.886775708919437*x16;
    double x124 = 49.886775708919437*x59;
    double x125 = n2*x124;
    double x126 = x53 - x55;
    double x127 = -x16*(-x12 - x13);
    double x128 = 24.943387854459719*x3;
    double x129 = x127 + x128 + x72;
    double x130 = 2.0*x76;
    double x131 = -24.943387854459719*x16;
    double x132 = x126 + x2*(x125 + x131);
    double x133 = x127 + x132;
    double x134 = x89 + x97;
    double x135 = x109 + x14*x16 - 58.201238327072687*x3;
    double x136 = x118 - 41.572313065822811*x3;
    double x137 = n3*x124;
    double x138 = x46 - x48 + x52 + x55;
    double x139 = x138 + x2*(x131 + x137);
    double x140 = x106 + x89;
    double x141 = 16.628925236306479*n4*x59;
    double x142 = -x26 - x27 - x28;
    double x143 = x2*x25;
    double x144 = x111*x143*(x58 + x60*x86) - x112*x143/((x86)*(x86)) + x113*x66 + x114 - x142*x16;
    double x145 = x144 + x49 - x50 + x56;

result[0] = 2.0*x44 + x75*(x64 + x72) - 2.0*x77 + x89;
result[1] = x103 + x97;
result[2] = x103 + x106;
result[3] = x102 + x117 + x75*(x109 - 33.257850472612958*x3 + x64);
result[4] = x102 + x121 + x75*(x101 + x118 - 66.515700920282541*x3);
result[5] = x122*x95 - x130*x95 + x75*(x126 + x129 + x2*(x123 + x125) + x91/n2) + x89;
result[6] = x106 + x134 + x75*(-x128 + x133 + x72);
result[7] = x117 + x134 + x75*(x132 + x135);
result[8] = x121 + x134 + x75*(x133 + x136);
result[9] = x105*x122 - x105*x130 + x75*(x129 + x138 + x2*(x123 + x137) + x104/n3) + x89;
result[10] = x117 + x140 + x75*(x135 + x139);
result[11] = x121 + x140 + x75*(x127 + x136 + x139);
result[12] = x116*x122 - x116*x130 + x75*(x145 - x16*(x142 + x65) + x2*(x141 - 16.628925236306479*x16) + 8.3144626181532395*x3 + x63 + x110/n4) + x89;
result[13] = x117 + x121 + x75*(x100 + x145 + x16*x36 + x2*(x141 - 8.3144626181532395*x16) - 58.201238302129291*x3) + x89;
result[14] = x120*x122 - x120*x130 + x75*(x127 + x144 + x2*(49.886775659032665*n5*x59 - 49.886775659032665*x16) + 24.943387829516332*x3 - x46 + x48 + x51 + x55 + x71 + x119/n5) + x89;
}
        
static void coder_d4gdn3dt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n4;
    double x1 = n2 + n3;
    double x2 = n5 + x0 + x1;
    double x3 = 1.0/x2;
    double x4 = pow(x2, -2);
    double x5 = n5*x4;
    double x6 = 24.943387829516332*x5;
    double x7 = 24.943387854459719*x4;
    double x8 = n3*x7;
    double x9 = x6 + x8;
    double x10 = 8.3144626181532395*n1;
    double x11 = 8.3144626181532395*n2;
    double x12 = 8.3144626181532395*n3;
    double x13 = x10 + x11 + x12;
    double x14 = n1 + x1;
    double x15 = 1.0/x14;
    double x16 = -x14*x4 + x3;
    double x17 = x15*x16;
    double x18 = x13*x17;
    double x19 = 8.3144626181532395*n4;
    double x20 = x19*x4;
    double x21 = 8.3144626181532395*n5;
    double x22 = x13 + x21;
    double x23 = n5 + x14;
    double x24 = 1.0/x23;
    double x25 = -x23*x4 + x3;
    double x26 = x24*x25;
    double x27 = x22*x26;
    double x28 = x20 + x27;
    double x29 = x18 + x28;
    double x30 = x29 + x9;
    double x31 = n2*x4;
    double x32 = 24.943387854459719*x31;
    double x33 = 24.943387854459719*n1;
    double x34 = 24.943387854459719*n4;
    double x35 = x33 + x34;
    double x36 = 1.0/x0;
    double x37 = -x0*x4 + x3;
    double x38 = x36*x37;
    double x39 = x35*x38;
    double x40 = x32 + x39;
    double x41 = x30 + x40;
    double x42 = 49.886775708919437*x38;
    double x43 = -2*x4;
    double x44 = pow(x2, -3);
    double x45 = 2*x44;
    double x46 = x0*x45 + x43;
    double x47 = x35*x36;
    double x48 = x46*x47;
    double x49 = pow(x0, -2);
    double x50 = x37*x49;
    double x51 = x35*x50;
    double x52 = x2*x42 + x2*x48 - x2*x51;
    double x53 = x41 + x52;
    double x54 = -x21;
    double x55 = -x19 + x54;
    double x56 = 16.628925236306479*x17;
    double x57 = x14*x45 + x43;
    double x58 = x13*x15;
    double x59 = x57*x58;
    double x60 = pow(x14, -2);
    double x61 = x16*x60;
    double x62 = x13*x61;
    double x63 = 16.628925236306479*x26;
    double x64 = x23*x45 + x43;
    double x65 = x22*x24;
    double x66 = x64*x65;
    double x67 = pow(x23, -2);
    double x68 = x25*x67;
    double x69 = x22*x68;
    double x70 = x2*x63 + x2*x66 - x2*x69;
    double x71 = x2*x56 + x2*x59 - x2*x62 - x4*x55 + x70;
    double x72 = x53 + x71;
    double x73 = x3*x72;
    double x74 = (*endmember[0].dmu0dT)(T, P);
    double x75 = 24.943387854459719*log(x0*x3);
    double x76 = n5*x3;
    double x77 = -24.943387829516332*x76;
    double x78 = n2*x3;
    double x79 = -24.943387854459719*x78;
    double x80 = n3*x3;
    double x81 = -24.943387854459719*x80;
    double x82 = x2*x39 + x75 + x77 + x79 + x81;
    double x83 = 8.3144626181532395*log(x14*x3);
    double x84 = x19 + x21;
    double x85 = 8.3144626181532395*log(x23*x3);
    double x86 = -x19*x3 + x2*x27 + x85;
    double x87 = x18*x2 - x3*x84 + x83 + x86;
    double x88 = x74 + x82 + x87;
    double x89 = x4*x88;
    double x90 = 74.830163563379159*x2;
    double x91 = 6*x44;
    double x92 = pow(x2, -4);
    double x93 = 6*x92;
    double x94 = x2*x47;
    double x95 = x35*x49;
    double x96 = 2*x2;
    double x97 = x36*x46*x90 + 74.830163563379159*x38 - x46*x95*x96 + 2*x48 - x50*x90 - 2*x51 + x94*(-x0*x93 + x91) + x35*x37*x96/((x0)*(x0)*(x0));
    double x98 = 24.943387854459719*x2;
    double x99 = x2*x65;
    double x100 = x2*x22;
    double x101 = x100*x67;
    double x102 = 2*x100*x25/((x23)*(x23)*(x23)) - 2*x101*x64 + x24*x64*x98 + 24.943387854459719*x26 + 2*x66 - x68*x98 - 2*x69 + x99*(-x23*x93 + x91);
    double x103 = 49.886775708919437*x44;
    double x104 = n2*x103;
    double x105 = -x104;
    double x106 = -x45*x84;
    double x107 = 16.628925236306479*x44;
    double x108 = n4*x107;
    double x109 = -x108;
    double x110 = n5*x44;
    double x111 = 49.886775659032665*x110;
    double x112 = -x111;
    double x113 = n3*x103;
    double x114 = -x113;
    double x115 = x112 + x114;
    double x116 = x109 + x115;
    double x117 = x106 + x116;
    double x118 = x105 + x117;
    double x119 = x2*x58;
    double x120 = x13*x2;
    double x121 = x120*x60;
    double x122 = x119*(-x14*x93 + x91) + 2*x120*x16/((x14)*(x14)*(x14)) - 2*x121*x57 + x15*x57*x98 + 24.943387854459719*x17 + 2*x59 - x61*x98 - 2*x62;
    double x123 = x118 + x122;
    double x124 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x125 = 1.0*x3;
    double x126 = x124*x125;
    double x127 = 6.0*x44;
    double x128 = x124*x88;
    double x129 = x124*x4;
    double x130 = x129*x72;
    double x131 = (*endmember[1].dmu0dT)(T, P);
    double x132 = (*endmember[2].dmu0dT)(T, P);
    double x133 = (*endmember[3].dmu0dT)(T, P);
    double x134 = (*endmember[4].dmu0dT)(T, P);
    double x135 = 24.943387854459719*log(x78);
    double x136 = 24.943387854459719*log(x80);
    double x137 = log(n4*x3);
    double x138 = log(x76);
    double x139 = n4 + n5;
    double x140 = 8.3144626181532395*log(x139*x3);
    double x141 = n1*x74 + n2*x131 + n2*x135 + n3*x132 + n3*x136 + n4*x133 + n5*x134 + 8.3144626098387775*n5*(x138 - 1.0986122896681101) + 16.628925219677555*n5*(x138 - 0.40546510910816402) + x0*x75 + x137*x19 + x139*x140 + x14*x83 + x23*x85;
    double x142 = -6.0*x124*x141*x92 + x127*x141;
    double x143 = -x4;
    double x144 = -n1;
    double x145 = -n4 + x144;
    double x146 = x143 - x145*x45;
    double x147 = x146*x47;
    double x148 = x147*x2;
    double x149 = x148 + x41;
    double x150 = x149 - 49.886775708919437*x3 + x71;
    double x151 = 2.0*x3;
    double x152 = 2.0*x4;
    double x153 = x124*x152;
    double x154 = x142 + x150*x151 - x150*x153;
    double x155 = x3 - x31;
    double x156 = x155*x98;
    double x157 = -x3*x35;
    double x158 = x157 + x81;
    double x159 = x77 + x87;
    double x160 = x131 + x135 + x156 + x158 + x159;
    double x161 = 2.0*x44;
    double x162 = x124*x160;
    double x163 = -x152*x160 + x161*x162;
    double x164 = x154 + x163;
    double x165 = x146*x2;
    double x166 = 4*x44;
    double x167 = 2*n1;
    double x168 = 2*n4;
    double x169 = 3*x92;
    double x170 = -x169*(x167 + x168);
    double x171 = x147 - x165*x95 + x48 - x51;
    double x172 = 49.886775708919437*x165*x36 + x171 + x42 + x94*(x166 + x170);
    double x173 = x102 + x172;
    double x174 = 4.0*x44;
    double x175 = x128*x174 - 1.0*x130 + 1.0*x73 - 4.0*x89;
    double x176 = x126*(x123 + x173 + x7) + x175;
    double x177 = -n3*x4 + x3;
    double x178 = x177*x98;
    double x179 = x132 + x136 + x157 + x159 + x178 + x79;
    double x180 = x124*x161;
    double x181 = -x152*x179 + x179*x180;
    double x182 = 16.628925236306479*x4;
    double x183 = 16.628925236306479*x2;
    double x184 = -n2 - n3 + x144;
    double x185 = -n5 + x184;
    double x186 = x143 - x185*x45;
    double x187 = x186*x24;
    double x188 = x183*x187;
    double x189 = 2*n5;
    double x190 = 2*n2;
    double x191 = 2*n3;
    double x192 = x167 + x190 + x191;
    double x193 = -x169*(x189 + x192);
    double x194 = x99*(x166 + x193);
    double x195 = x186*x65;
    double x196 = -x101*x186 + x195 + x66 - x69;
    double x197 = x143 - x184*x45;
    double x198 = x197*x58;
    double x199 = -x169*x192;
    double x200 = x119*(x166 + x199) - x121*x197 + x15*x183*x197 + x198 + x45*x55 + x56 + x59 - x62;
    double x201 = x116 + x200;
    double x202 = x105 + x201;
    double x203 = x188 + x194 + x196 + x202 + x63;
    double x204 = x142 + x175;
    double x205 = x198*x2 + x4*x84;
    double x206 = x195*x2 + x205;
    double x207 = x206 - 33.257850472612958*x3 + x53;
    double x208 = x151*x207 - x153*x207;
    double x209 = -8.3144626181532395*n4*x4 + 8.3144626181532395*x3;
    double x210 = x2*x209;
    double x211 = 1.0/x139;
    double x212 = -x139*x4 + x3;
    double x213 = x211*x212;
    double x214 = x213*x84;
    double x215 = -x13*x3 + x140 + x2*x214;
    double x216 = x133 + 8.3144626181532395*x137 + x210 + x215 - x22*x3 + x82;
    double x217 = -x152*x216 + x180*x216;
    double x218 = 33.257850447669568*x4;
    double x219 = x205 + x70;
    double x220 = x149 + x219 - 66.515700920282541*x3;
    double x221 = x151*x220 - x153*x220;
    double x222 = 24.943387829516332*x3 - 24.943387829516332*x5;
    double x223 = x2*x222;
    double x224 = x134 + 24.943387829516332*x138 + x158 + x215 + x223 + x79 + x86 - 15.876819783702929;
    double x225 = -x152*x224 + x180*x224;
    double x226 = 2*x147 + x94*(x170 + x45);
    double x227 = x102 + x226;
    double x228 = x128*x161 - 2.0*x89;
    double x229 = x126*(x123 + x227 + 74.830163563379159*x4) + x228;
    double x230 = x154 + x229;
    double x231 = -49.886775708919437*x4;
    double x232 = 1.0/n2;
    double x233 = x30 - x32;
    double x234 = -x33 - x34;
    double x235 = -x234*x4;
    double x236 = 24.943387854459719*x3;
    double x237 = x235 + x236 + x71;
    double x238 = x156*x232 + x2*(x104 + x231) + x233 + x237;
    double x239 = 4.0*x4;
    double x240 = 1.0*x129;
    double x241 = x125*x238 - x160*x239 + x162*x174 - x238*x240;
    double x242 = -x7;
    double x243 = x2*(x104 + x242) + x233;
    double x244 = x235 + x243;
    double x245 = -x236 + x244 + x71;
    double x246 = x125*x245 + x181 - x240*x245;
    double x247 = x206 - 58.201238327072687*x3 + x35*x4;
    double x248 = x243 + x247;
    double x249 = x125*x248 + x217 - x240*x248;
    double x250 = x142 + x228;
    double x251 = x125*x150 - x150*x240 + x250;
    double x252 = x163 + x251;
    double x253 = x125*x207 - x207*x240;
    double x254 = x126*(x146*x36*x98 + x171 + x203 + 24.943387854459719*x38 + 66.515700945225916*x4 + x94*(x145*x93 + x166)) + x253;
    double x255 = x219 - 41.572313065822811*x3;
    double x256 = x244 + x255;
    double x257 = x125*x256 - x240*x256;
    double x258 = x125*x220 - x220*x240 + x225;
    double x259 = x126*(x202 + x227 + 83.144626156588998*x4) + x258;
    double x260 = 1.0/n3;
    double x261 = x29 + x32 + x6 - x8;
    double x262 = x178*x260 + x2*(x113 + x231) + x237 + x261;
    double x263 = x124*x174;
    double x264 = x125*x262 - x179*x239 + x179*x263 - x240*x262;
    double x265 = x181 + x251;
    double x266 = x2*(x113 + x242) + x261;
    double x267 = x247 + x266;
    double x268 = x125*x267 - x240*x267;
    double x269 = x217 + x268;
    double x270 = x235 + x255 + x266;
    double x271 = x125*x270 - x240*x270;
    double x272 = 49.886775708919444*x4;
    double x273 = 2*x195 + x99*(x193 + x45);
    double x274 = 2*x198;
    double x275 = x119*(x199 + x45);
    double x276 = x118 + x274 + x275;
    double x277 = -x10 - x11 - x12;
    double x278 = x277 + x54;
    double x279 = 1.0/n4;
    double x280 = x2*x84;
    double x281 = x139*x45 + x43;
    double x282 = x211*x281;
    double x283 = pow(x139, -2);
    double x284 = x212*x283;
    double x285 = x183*x213 + x214 - x277*x4 + x280*x282 - x280*x284;
    double x286 = -x20 + x285 + x40 + x9;
    double x287 = x2*(x108 - x182) + x210*x279 - x278*x4 + x286 + 8.3144626181532395*x3 + x52;
    double x288 = x125*x287 - x216*x239 + x216*x263 - x240*x287;
    double x289 = 8.3144626181532395*x2;
    double x290 = x187*x289 + x196 + 8.3144626181532395*x26 + 66.515700920282541*x4 + x99*(x166 + x185*x93);
    double x291 = x148 + x2*(x108 - 8.3144626181532395*x4) + x22*x4 + x286 - 58.201238302129291*x3;
    double x292 = x125*x291 - x240*x291;
    double x293 = 49.886775659032665*x4;
    double x294 = 1.0/n5;
    double x295 = x2*(x111 - x293) + x223*x294 + x235 + x28 + x285 + 24.943387829516332*x3 + x32 - x6 + x70 + x8;
    double x296 = x125*x295 - x224*x239 + x224*x263 - x240*x295;
    double x297 = 3.0*x3;
    double x298 = 6.0*x4;
    double x299 = -99.773551417838874*x4;
    double x300 = 149.66032712675832*x44;
    double x301 = 149.66032712675832*x92;
    double x302 = -n2*x301;
    double x303 = x190*x44;
    double x304 = x232*x98;
    double x305 = 99.773551417838874*x44;
    double x306 = n2*x305;
    double x307 = 24.943387854459719*x155*x232 + x306;
    double x308 = -x35*x45;
    double x309 = x102 + x308;
    double x310 = x122 + x309;
    double x311 = x117 + x310;
    double x312 = 3.0*x129;
    double x313 = x2*(x302 + x305) + x304*(x143 + x303) + x307;
    double x314 = x142 + x241;
    double x315 = x151*x245 - x153*x245;
    double x316 = -33.257850472612958*x4;
    double x317 = x188 + x194 + x196 + x234*x45 + x63;
    double x318 = x316 + x317;
    double x319 = x201 + x313;
    double x320 = x151*x248 - x153*x248;
    double x321 = x309 - 41.572313115709584*x4;
    double x322 = x151*x256 - x153*x256;
    double x323 = x2*(x103 + x302) + x306;
    double x324 = x117 + x323;
    double x325 = x142 + x163;
    double x326 = x201 + x323;
    double x327 = x246 + x325;
    double x328 = x225 + x257;
    double x329 = x274 + x275 + x308;
    double x330 = x273 + x329 + 74.830163563379173*x4;
    double x331 = x324 + x329;
    double x332 = x102 + x293;
    double x333 = 24.943387854459719*x177*x260;
    double x334 = -n3*x301;
    double x335 = x191*x44;
    double x336 = x260*x98;
    double x337 = x105 + x109;
    double x338 = n3*x305 + x112 + x337;
    double x339 = x106 + x338;
    double x340 = x124*x127;
    double x341 = x2*(x305 + x334) + x200 + x333 + x336*(x143 + x335) + x338;
    double x342 = x142 + x264;
    double x343 = x151*x267 - x153*x267;
    double x344 = x151*x270 - x153*x270;
    double x345 = x2*(x103 + x334) + x339;
    double x346 = x142 + x181;
    double x347 = x329 + x345;
    double x348 = x209*x279;
    double x349 = -49.886775708919437*n4*x92;
    double x350 = x168*x44;
    double x351 = x279*x289;
    double x352 = 33.257850472612958*x44;
    double x353 = 2*x84;
    double x354 = 2*x280;
    double x355 = -x13*x45 + x211*x280*(-x139*x93 + x91) + 24.943387854459719*x213 - x281*x283*x354 + x282*x353 + x282*x98 - x284*x353 - x284*x98 + x212*x354/((x139)*(x139)*(x139));
    double x356 = n4*x352 + x105 + x115 + x355;
    double x357 = -x22*x45 + x356;
    double x358 = x142 + x151*x291 - x153*x291;

result[0] = x126*(x102 + x123 + x97) + x127*x128 - 3.0*x130 + x142 + 3.0*x73 - 6.0*x89;
result[1] = x164 + x176;
result[2] = x154 + x176 + x181;
result[3] = x126*(x182 + x203 + x97) + x204 + x208 + x217;
result[4] = x126*(x173 + x202 + x218) + x204 + x221 + x225;
result[5] = x230 + x241;
result[6] = x164 + x229 + x246;
result[7] = x249 + x252 + x254;
result[8] = x252 + x257 + x259;
result[9] = x230 + x264;
result[10] = x254 + x265 + x269;
result[11] = x259 + x265 + x271;
result[12] = x126*(x272 + x273 + x276 + x97) + x208 + x250 + x288;
result[13] = x126*(x172 + x276 + x290) + x217 + x250 + x253 + x258 + x292;
result[14] = x126*(x227 + x276 + 99.773551367952109*x4) + x221 + x250 + x296;
result[15] = x126*(x2*(x300 + x302) + x299 + x304*(x303 + x43) + x307 + x311 - x156/((n2)*(n2))) + x127*x162 + x142 - x160*x298 + x238*x297 - x238*x312;
result[16] = x126*(-x272 + x311 + x313) + x181 + x314 + x315;
result[17] = x126*(x318 + x319) + x217 + x314 + x320;
result[18] = x126*(x319 + x321) + x225 + x314 + x322;
result[19] = x126*(x310 + x324 + x7) + x264 + x315 + x325;
result[20] = x126*(x317 + x326 + 41.572313090766201*x4) + x249 + x268 + x327;
result[21] = x126*(x218 + x309 + x326) + x271 + x327 + x328;
result[22] = x126*(x324 + x330) + x288 + x320 + x325;
result[23] = x126*(x290 + x331) + x249 + x292 + x325 + x328;
result[24] = x126*(x331 + x332) + x296 + x322 + x325;
result[25] = x126*(x2*(x300 + x334) + x299 + x310 + x333 + x336*(x335 + x43) + x339 - x178/((n3)*(n3))) + x142 - x179*x298 + x179*x340 + x262*x297 - x262*x312;
result[26] = x126*(x318 + x341) + x217 + x342 + x343;
result[27] = x126*(x321 + x341) + x225 + x342 + x344;
result[28] = x126*(x330 + x345) + x288 + x343 + x346;
result[29] = x126*(x290 + x347) + x225 + x269 + x271 + x292 + x346;
result[30] = x126*(x332 + x347) + x296 + x344 + x346;
result[31] = x126*(x2*(x103 + x349) + x316 + x348 + x351*(x350 + x43) + x357 + x97 - x210/((n4)*(n4))) + x142 - x216*x298 + x216*x340 + x287*x297 - x287*x312;
result[32] = x126*(x172 + x2*(x349 + x352) + x278*x45 + x348 + x351*(x143 + x350) + x356 + 8.3144625932098535*x4) + x225 + x288 + x358;
result[33] = x126*(x2*(x107 + x349) + x226 + x357 + 83.144626131645623*x4) + x217 + x296 + x358;
result[34] = x126*(99.77355131806533*x110 + x114 + 24.943387829516332*x2*x294*(x189*x44 + x43) + x2*(-149.66032697709801*n5*x92 + 149.66032697709801*x44) + x222*x294 + x309 + x337 + x355 - 99.77355131806533*x4 - x223/((n5)*(n5))) + x142 - x224*x298 + x224*x340 + x295*x297 - x295*x312;
}
        
static double coder_dgdp(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = 1.0*(1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5)*(n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P) + n4*(*endmember[3].dmu0dP)(T, P) + n5*(*endmember[4].dmu0dP)(T, P))/(n1 + n2 + n3 + n4 + n5);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x2 = n1 + n2 + n3 + n4 + n5;
    double x3 = 1.0/x2;
    double x4 = x1*x3;
    double x5 = (*endmember[1].dmu0dP)(T, P);
    double x6 = (*endmember[2].dmu0dP)(T, P);
    double x7 = (*endmember[3].dmu0dP)(T, P);
    double x8 = (*endmember[4].dmu0dP)(T, P);
    double x9 = n1*x0 + n2*x5 + n3*x6 + n4*x7 + n5*x8;
    double x10 = -1.0*x1*x9/((x2)*(x2)) + x3*x9;

result[0] = x0*x4 + x10;
result[1] = x10 + x4*x5;
result[2] = x10 + x4*x6;
result[3] = x10 + x4*x7;
result[4] = x10 + x4*x8;
}
        
static void coder_d3gdn2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3 + n4 + n5;
    double x2 = 1.0/x1;
    double x3 = x0*x2;
    double x4 = pow(x1, -2);
    double x5 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x6 = x4*x5;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = (*endmember[3].dmu0dP)(T, P);
    double x11 = (*endmember[4].dmu0dP)(T, P);
    double x12 = 2.0*n1*x0 + 2.0*n2*x8 + 2.0*n3*x9 + 2.0*n4*x10 + 2.0*n5*x11;
    double x13 = -x12*x4 + x12*x5/((x1)*(x1)*(x1));
    double x14 = x13 + 1.0*x3 - 1.0*x7;
    double x15 = 1.0*x2;
    double x16 = 1.0*x6;
    double x17 = x15*x8 - x16*x8;
    double x18 = x15*x9 - x16*x9;
    double x19 = x10*x15 - x10*x16;
    double x20 = x11*x15 - x11*x16;
    double x21 = 2.0*x2;
    double x22 = 2.0*x6;
    double x23 = x13 + x17;
    double x24 = x13 + x18;

result[0] = x13 + 2.0*x3 - 2.0*x7;
result[1] = x14 + x17;
result[2] = x14 + x18;
result[3] = x14 + x19;
result[4] = x14 + x20;
result[5] = x13 + x21*x8 - x22*x8;
result[6] = x18 + x23;
result[7] = x19 + x23;
result[8] = x20 + x23;
result[9] = x13 + x21*x9 - x22*x9;
result[10] = x19 + x24;
result[11] = x20 + x24;
result[12] = x10*x21 - x10*x22 + x13;
result[13] = x13 + x19 + x20;
result[14] = x11*x21 - x11*x22 + x13;
}
        
static void coder_d4gdn3dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3 + n4 + n5;
    double x2 = pow(x1, -2);
    double x3 = x0*x2;
    double x4 = pow(x1, -3);
    double x5 = 6.0*x4;
    double x6 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4 + 1.0*n5;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = (*endmember[3].dmu0dP)(T, P);
    double x11 = (*endmember[4].dmu0dP)(T, P);
    double x12 = n1*x0 + n2*x8 + n3*x9 + n4*x10 + n5*x11;
    double x13 = x12*x5 - 6.0*x12*x6/((x1)*(x1)*(x1)*(x1));
    double x14 = 4.0*x4;
    double x15 = x13 + x14*x7 - 4.0*x3;
    double x16 = 2.0*x2;
    double x17 = 2.0*x4;
    double x18 = x6*x8;
    double x19 = -x16*x8 + x17*x18;
    double x20 = x17*x6;
    double x21 = -x16*x9 + x20*x9;
    double x22 = -x10*x16 + x10*x20;
    double x23 = -x11*x16 + x11*x20;
    double x24 = x13 + x17*x7 - 2.0*x3;
    double x25 = 4.0*x2;
    double x26 = x14*x18 - x25*x8;
    double x27 = x19 + x24;
    double x28 = x14*x6;
    double x29 = -x25*x9 + x28*x9;
    double x30 = x21 + x24;
    double x31 = -x10*x25 + x10*x28;
    double x32 = x22 + x23;
    double x33 = -x11*x25 + x11*x28;
    double x34 = 6.0*x2;
    double x35 = x13 + x26;
    double x36 = x13 + x19;
    double x37 = x21 + x36;
    double x38 = x5*x6;
    double x39 = x13 + x29;
    double x40 = x13 + x21;

result[0] = x13 - 6.0*x3 + x5*x7;
result[1] = x15 + x19;
result[2] = x15 + x21;
result[3] = x15 + x22;
result[4] = x15 + x23;
result[5] = x24 + x26;
result[6] = x21 + x27;
result[7] = x22 + x27;
result[8] = x23 + x27;
result[9] = x24 + x29;
result[10] = x22 + x30;
result[11] = x23 + x30;
result[12] = x24 + x31;
result[13] = x24 + x32;
result[14] = x24 + x33;
result[15] = x13 + x18*x5 - x34*x8;
result[16] = x21 + x35;
result[17] = x22 + x35;
result[18] = x23 + x35;
result[19] = x29 + x36;
result[20] = x22 + x37;
result[21] = x23 + x37;
result[22] = x31 + x36;
result[23] = x32 + x36;
result[24] = x33 + x36;
result[25] = x13 - x34*x9 + x38*x9;
result[26] = x22 + x39;
result[27] = x23 + x39;
result[28] = x31 + x40;
result[29] = x32 + x40;
result[30] = x33 + x40;
result[31] = -x10*x34 + x10*x38 + x13;
result[32] = x13 + x23 + x31;
result[33] = x13 + x22 + x33;
result[34] = -x11*x34 + x11*x38 + x13;
}
        
static double coder_d2gdt2(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = 1.0*(n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P) + n5*(*endmember[4].d2mu0dT2)(T, P));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d2mu0dT2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dT2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dT2)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dT2)(T, P);
result[4] = 1.0*(*endmember[4].d2mu0dT2)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d2mu0dTdP)(T, P) + 1.0*n2*(*endmember[1].d2mu0dTdP)(T, P) + 1.0*n3*(*endmember[2].d2mu0dTdP)(T, P) + 1.0*n4*(*endmember[3].d2mu0dTdP)(T, P) + 1.0*n5*(*endmember[4].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d2mu0dTdP)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dTdP)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dTdP)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dTdP)(T, P);
result[4] = 1.0*(*endmember[4].d2mu0dTdP)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P) + n4*(*endmember[3].d2mu0dP2)(T, P) + n5*(*endmember[4].d2mu0dP2)(T, P));
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d2mu0dP2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dP2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dP2)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dP2)(T, P);
result[4] = 1.0*(*endmember[4].d2mu0dP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P) + n5*(*endmember[4].d3mu0dT3)(T, P));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d3mu0dT3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT3)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dT3)(T, P);
result[4] = 1.0*(*endmember[4].d3mu0dT3)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P) + n4*(*endmember[3].d3mu0dT2dP)(T, P) + n5*(*endmember[4].d3mu0dT2dP)(T, P));
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d3mu0dT2dP)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT2dP)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT2dP)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dT2dP)(T, P);
result[4] = 1.0*(*endmember[4].d3mu0dT2dP)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d3mu0dTdP2)(T, P) + 1.0*n2*(*endmember[1].d3mu0dTdP2)(T, P) + 1.0*n3*(*endmember[2].d3mu0dTdP2)(T, P) + 1.0*n4*(*endmember[3].d3mu0dTdP2)(T, P) + 1.0*n5*(*endmember[4].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d3mu0dTdP2)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dTdP2)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dTdP2)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dTdP2)(T, P);
result[4] = 1.0*(*endmember[4].d3mu0dTdP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P) + n4*(*endmember[3].d3mu0dP3)(T, P) + n5*(*endmember[4].d3mu0dP3)(T, P));
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = 1.0*(*endmember[0].d3mu0dP3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dP3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dP3)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dP3)(T, P);
result[4] = 1.0*(*endmember[4].d3mu0dP3)(T, P);
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

