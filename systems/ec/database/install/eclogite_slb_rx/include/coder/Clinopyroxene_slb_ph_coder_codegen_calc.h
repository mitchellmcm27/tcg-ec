#include <math.h>


static double coder_g(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = 1.0*n2;
    double x1 = 1.0*n3;
    double x2 = 1.0*n5;
    double x3 = 1.0*n1 + x0 + x1 + x2;
    double x4 = 3.5*n4 + x3;
    double x5 = n1 + n3;
    double x6 = n4 + n5;
    double x7 = 1.0/(n2 + x5 + x6);
    double x8 = 1.0*n4;
    double x9 = n1 + n2 + n4;

result = (0.027777777777777773*n1*(444600.0*n3 + 728000.0*n4 + 437400.0*n5) + 12349.999999999998*n2*n3 + 0.027777777777777773*n3*(444600.0*n1 + 444600.0*n2 + 1696800.0*n4) + 0.0085352842554488641*n4*(2369250.0*n1 + 5522175.0*n3 + 911250.0*n5) + 0.027777777777777773*n5*(437400.0*n1 + 280000.0*n4) + x4*(8.3144626181532395*T*(x0*log(n2*x7) + x1*log(n3*x7) + x2*log(n5*x7) + 1.0*x5*log(x5*x7) + 1.0*x6*log(x6*x7) + x8*(log(n4*x7) - 0.69314718055994495) + 1.0*x9*log(x7*x9) + (2.0*n1 + 2.0*n2 + 2.0*n3 + 2.0*n5 + x8)*log(x7*(0.5*n4 + x3))) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P) + n5*(*endmember[4].mu0)(T, P)))/x4;
    return result;
}
        
static void coder_dgdn(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 3.5*n4;
    double x1 = 1.0*n5;
    double x2 = 1.0*n2;
    double x3 = 1.0*n1;
    double x4 = 1.0*n3;
    double x5 = x3 + x4;
    double x6 = x1 + x2 + x5;
    double x7 = x0 + x6;
    double x8 = (*endmember[0].mu0)(T, P);
    double x9 = n1*x8;
    double x10 = (*endmember[1].mu0)(T, P);
    double x11 = n2*x10;
    double x12 = (*endmember[2].mu0)(T, P);
    double x13 = n3*x12;
    double x14 = (*endmember[3].mu0)(T, P);
    double x15 = (*endmember[4].mu0)(T, P);
    double x16 = n5*x15;
    double x17 = n1 + n3;
    double x18 = n4 + n5;
    double x19 = n2 + x17 + x18;
    double x20 = 1.0/x19;
    double x21 = log(n2*x20);
    double x22 = log(n3*x20);
    double x23 = log(n5*x20);
    double x24 = log(n4*x20);
    double x25 = 1.0*n4;
    double x26 = 1.0*log(x17*x20);
    double x27 = 1.0*log(x18*x20);
    double x28 = n1 + n2 + n4;
    double x29 = 1.0*log(x20*x28);
    double x30 = 2.0*n1 + 2.0*n2 + 2.0*n3 + 2.0*n5 + x25;
    double x31 = 0.5*n4 + x6;
    double x32 = log(x20*x31);
    double x33 = x1*x23 + x17*x26 + x18*x27 + x2*x21 + x22*x4 + x25*(x24 - 0.69314718055994495) + x28*x29 + x30*x32;
    double x34 = 8.3144626181532395*T;
    double x35 = x33*x34;
    double x36 = (0.027777777777777773*n1*(444600.0*n3 + 728000.0*n4 + 437400.0*n5) + 12349.999999999998*n2*n3 + 0.027777777777777773*n3*(444600.0*n1 + 444600.0*n2 + 1696800.0*n4) + 0.0085352842554488641*n4*(2369250.0*n1 + 5522175.0*n3 + 911250.0*n5) + 0.027777777777777773*n5*(437400.0*n1 + 280000.0*n4) + x7*(n4*x14 + x11 + x13 + x16 + x35 + x9))/((x7)*(x7));
    double x37 = -1.0*x36;
    double x38 = 1.0/x7;
    double x39 = 2.0*x32;
    double x40 = -x2*x20;
    double x41 = -x20*x4;
    double x42 = -x20*x25;
    double x43 = pow(x19, -2);
    double x44 = -x31*x43;
    double x45 = x19*x30/x31;
    double x46 = x45*(1.0*x20 + x44);
    double x47 = x39 + x40 + x41 + x42 + x46;
    double x48 = x1 + x25;
    double x49 = -x20*x48;
    double x50 = -x1*x20;
    double x51 = x2 + x25 + x3;
    double x52 = x19*x51*(x20 - x28*x43)/x28 + x29 + x50;
    double x53 = x49 + x52;
    double x54 = x26 + x19*x5*(-x17*x43 + x20)/x17;
    double x55 = x1*x15 + x10*x2 + x12*x4 + x14*x25 + x3*x8 + x35;
    double x56 = 24699.999999999996*n3 + x55;
    double x57 = 1.0*x19;
    double x58 = -x20*x5;
    double x59 = x41 + x58;
    double x60 = x39 + x42 + x46;
    double x61 = -x20*x51;
    double x62 = x27 + x19*x48*(-x18*x43 + x20)/x18;

result[0] = x37 + x38*(40444.444444444438*n4 + 24299.999999999996*n5 + x56 + x7*(x34*(x47 + x53 + x54) + x8));
result[1] = x37 + x38*(x56 + x7*(x10 + x34*(1.0*x21 + x53 + x57*(-n2*x43 + x20) + x59 + x60)));
result[2] = x37 + x38*(24699.999999999996*n1 + 24699.999999999996*n2 + 94266.666666666657*n4 + x55 + x7*(x12 + x34*(1.0*x22 + x40 + x49 + x50 + x54 + x57*(-n3*x43 + x20) + x60 + x61)));
result[3] = -3.5*x36 + x38*(29.100619163536336*T*x33 + 40444.444444444438*n1 + 94266.666666666657*n3 + 15555.555555555555*n5 + x0*x14 + 3.5*x11 + 3.5*x13 + 3.5*x16 + x7*(x14 + x34*(1.0*x24 + 1.0*x32 + x40 + x45*(0.5*x20 + x44) + x52 + x57*(-n4*x43 + x20) + x59 + x62 - 0.69314718055994495)) + 3.5*x9);
result[4] = x37 + x38*(24299.999999999996*n1 + 15555.555555555555*n4 + x55 + x7*(x15 + x34*(1.0*x23 + x47 + x57*(-n5*x43 + x20) + x58 + x61 + x62)));
}
        
static void coder_d2gdn2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 3.5*n4;
    double x1 = 1.0*n5;
    double x2 = 1.0*n2;
    double x3 = 1.0*n1;
    double x4 = 1.0*n3;
    double x5 = x3 + x4;
    double x6 = x1 + x2 + x5;
    double x7 = x0 + x6;
    double x8 = (*endmember[0].mu0)(T, P);
    double x9 = n1*x8;
    double x10 = (*endmember[1].mu0)(T, P);
    double x11 = n2*x10;
    double x12 = (*endmember[2].mu0)(T, P);
    double x13 = n3*x12;
    double x14 = (*endmember[3].mu0)(T, P);
    double x15 = (*endmember[4].mu0)(T, P);
    double x16 = n5*x15;
    double x17 = n1 + n3;
    double x18 = n4 + n5;
    double x19 = n2 + x17 + x18;
    double x20 = 1.0/x19;
    double x21 = log(n2*x20);
    double x22 = log(n3*x20);
    double x23 = log(n5*x20);
    double x24 = log(n4*x20);
    double x25 = 1.0*n4;
    double x26 = 1.0*log(x17*x20);
    double x27 = 1.0*log(x18*x20);
    double x28 = n1 + n2 + n4;
    double x29 = 1.0*log(x20*x28);
    double x30 = 2.0*n2;
    double x31 = 2.0*n3;
    double x32 = 2.0*n5;
    double x33 = 2.0*n1 + x25 + x30 + x31 + x32;
    double x34 = 0.5*n4 + x6;
    double x35 = log(x20*x34);
    double x36 = x1*x23 + x17*x26 + x18*x27 + x2*x21 + x22*x4 + x25*(x24 - 0.69314718055994495) + x28*x29 + x33*x35;
    double x37 = 8.3144626181532395*T;
    double x38 = x36*x37;
    double x39 = (0.027777777777777773*n1*(444600.0*n3 + 728000.0*n4 + 437400.0*n5) + 12349.999999999998*n2*n3 + 0.027777777777777773*n3*(444600.0*n1 + 444600.0*n2 + 1696800.0*n4) + 0.0085352842554488641*n4*(2369250.0*n1 + 5522175.0*n3 + 911250.0*n5) + 0.027777777777777773*n5*(437400.0*n1 + 280000.0*n4) + x7*(n4*x14 + x11 + x13 + x16 + x38 + x9))/((x7)*(x7)*(x7));
    double x40 = 2.0*x39;
    double x41 = pow(x7, -2);
    double x42 = 2.0*x35;
    double x43 = -x2*x20;
    double x44 = -x20*x4;
    double x45 = -x20*x25;
    double x46 = 1.0*x20;
    double x47 = pow(x19, -2);
    double x48 = -x34*x47;
    double x49 = x46 + x48;
    double x50 = 1.0/x34;
    double x51 = x33*x50;
    double x52 = x49*x51;
    double x53 = x19*x52;
    double x54 = x42 + x43 + x44 + x45 + x53;
    double x55 = x1 + x25;
    double x56 = -x20*x55;
    double x57 = -x1*x20;
    double x58 = x2 + x25 + x3;
    double x59 = 1.0/x28;
    double x60 = x20 - x28*x47;
    double x61 = x59*x60;
    double x62 = x58*x61;
    double x63 = x19*x62 + x29 + x57;
    double x64 = x56 + x63;
    double x65 = 1.0/x17;
    double x66 = -x17*x47 + x20;
    double x67 = x65*x66;
    double x68 = x5*x67;
    double x69 = x19*x68 + x26;
    double x70 = T*(x54 + x64 + x69);
    double x71 = 8.3144626181532395*x70;
    double x72 = x1*x15 + x10*x2 + x12*x4 + x14*x25 + x3*x8 + x38;
    double x73 = 24699.999999999996*n3 + x72;
    double x74 = x41*(40444.444444444438*n4 + 24299.999999999996*n5 + x7*(x71 + x8) + x73);
    double x75 = 1.0/x7;
    double x76 = x49*x50;
    double x77 = -2.0*x47;
    double x78 = pow(x19, -3);
    double x79 = 2*x78;
    double x80 = x34*x79;
    double x81 = x19*x51;
    double x82 = 1.0*x19;
    double x83 = x33/((x34)*(x34));
    double x84 = x49*x83;
    double x85 = 4.0*x19*x76 + x81*(x77 + x80) - x82*x84;
    double x86 = 2.0*x19;
    double x87 = x61*x86;
    double x88 = -2*x47;
    double x89 = x19*x58;
    double x90 = x59*x89;
    double x91 = x90*(x28*x79 + x88);
    double x92 = -x60*x89/((x28)*(x28));
    double x93 = x2*x47;
    double x94 = x4*x47;
    double x95 = x1*x47;
    double x96 = x62 + x93 + x94 + x95;
    double x97 = x87 + x91 + x92 + x96;
    double x98 = x85 + x97;
    double x99 = x25*x47;
    double x100 = -x25;
    double x101 = -x47*(-x1 + x100) + x52 + x99;
    double x102 = x19*x5;
    double x103 = x102*x65;
    double x104 = -x102*x66/((x17)*(x17)) + x103*(x17*x79 + x88) + x67*x86 + x68;
    double x105 = x101 + x104;
    double x106 = x37*x7;
    double x107 = -x47;
    double x108 = -n1;
    double x109 = x103*(x107 - x79*(-n3 + x108));
    double x110 = x109 + x68;
    double x111 = -2.0*x20;
    double x112 = x101 + x111;
    double x113 = x82*(-n2*x47 + x20);
    double x114 = -x20*x5;
    double x115 = x114 + x44;
    double x116 = x42 + x45 + x53;
    double x117 = x113 + x115 + x116 + 1.0*x21 + x64;
    double x118 = x117*x37;
    double x119 = 1.0*x10 + x118;
    double x120 = x71 + 1.0*x8;
    double x121 = x7*(x10 + x118) + x73;
    double x122 = 1.0*x41;
    double x123 = -x121*x122;
    double x124 = x40 - 1.0*x74;
    double x125 = x90*(x107 - x79*(-n2 - n4 + x108));
    double x126 = x125 + x96;
    double x127 = x126 + x85;
    double x128 = -x20*x58;
    double x129 = x82*(-n3*x47 + x20);
    double x130 = x116 + x128 + x129 + 1.0*x22 + x43 + x56 + x57 + x69;
    double x131 = x130*x37;
    double x132 = 1.0*x12 + x131;
    double x133 = x132 + 24699.999999999996;
    double x134 = 24699.999999999996*n1 + 24699.999999999996*n2 + 94266.666666666657*n4 + x7*(x12 + x131) + x72;
    double x135 = -x122*x134;
    double x136 = -3.0*x20;
    double x137 = x47*x55 + x52 + x99;
    double x138 = x110 + x137;
    double x139 = 0.5*x19;
    double x140 = 0.5*x20 + x48;
    double x141 = x140*x50*x86;
    double x142 = x141 + x76*x82 + x81*(-1.5*x47 + x80);
    double x143 = -x139*x84 + x142;
    double x144 = x82*(-n4*x47 + x20);
    double x145 = x140*x51;
    double x146 = 1.0/x18;
    double x147 = -x18*x47 + x20;
    double x148 = x146*x147;
    double x149 = x148*x55;
    double x150 = x149*x19 + x27;
    double x151 = x115 + x144 + x145*x19 + x150 + 1.0*x24 + 1.0*x35 + x43 + x63 - 0.69314718055994495;
    double x152 = x151*x37;
    double x153 = 1.0*x14 + x152;
    double x154 = 29.100619163536336*T;
    double x155 = 40444.444444444438*n1 + 94266.666666666657*n3 + 15555.555555555555*n5 + x0*x14 + 3.5*x11 + 3.5*x13 + x154*x36 + 3.5*x16 + x7*(x14 + x152) + 3.5*x9;
    double x156 = -x122*x155 + 7.0*x39;
    double x157 = x109 + x137 - 4.0*x20 + x68;
    double x158 = x82*(-n5*x47 + x20);
    double x159 = x114 + x128 + x150 + x158 + 1.0*x23 + x54;
    double x160 = x159*x37;
    double x161 = 1.0*x15 + x160;
    double x162 = 24299.999999999996*n1 + 15555.555555555555*n4 + x7*(x15 + x160) + x72;
    double x163 = -x122*x162;
    double x164 = x121*x41;
    double x165 = 16.628925236306479*T;
    double x166 = -x93;
    double x167 = x30*x78;
    double x168 = x87 + x91 + x92;
    double x169 = x85 + x95;
    double x170 = -x3;
    double x171 = -x47*(x170 - x4);
    double x172 = x171 + x46;
    double x173 = x172 + x94;
    double x174 = x136 + x85 + x95;
    double x175 = -1.0*x47;
    double x176 = x166 + x19*(x167 + x175) + x62 + x94;
    double x177 = x125 + x174 + x176;
    double x178 = x123 + x40;
    double x179 = x137 + x171;
    double x180 = x143 + x95;
    double x181 = x134*x41;
    double x182 = -x94;
    double x183 = x31*x78;
    double x184 = -x47*(x100 + x170 - x2) + x93;
    double x185 = x182 + x19*(x175 + x183);
    double x186 = 2.0*n4*x78;
    double x187 = x140*x83;
    double x188 = x19*x55;
    double x189 = x146*x188*(x18*x79 + x88) - x147*x188/((x18)*(x18)) + x148*x86 + x149;
    double x190 = x145 + x189 - x99;
    double x191 = x162*x41;

result[0] = x40 - 2.0*x74 + x75*(x106*(x105 + x98) + 16.628925236306479*x70 + 2.0*x8);
result[1] = x123 + x124 + x75*(x106*(x110 + x112 + x98) + x119 + x120);
result[2] = x124 + x135 + x75*(x106*(x104 + x112 + x127) + x120 + x133);
result[3] = x156 - 3.5*x74 + x75*(x106*(x136 + x138 + x143 + x97) + x153 + 29.100619163536336*x70 + 3.5*x8 + 40444.444444444438);
result[4] = x124 + x163 + x75*(x106*(x127 + x157) + x120 + x161 + 24299.999999999996);
result[5] = -2.0*x164 + x40 + x75*(2.0*x10 + x106*(x101 + x166 + x168 + x169 + x173 + x19*(x167 + x77) + x62 + x113/n2) + x117*x165);
result[6] = x135 + x178 + x75*(x106*(x101 + x177 + x47*x5) + x119 + x133);
result[7] = x156 - 3.5*x164 + x75*(3.5*x10 + x106*(x111 + x168 + x176 + x179 + x180) + x117*x154 + x153);
result[8] = x163 + x178 + x75*(x106*(x177 + x179) + x119 + x161);
result[9] = -2.0*x181 + x40 + x75*(x106*(x105 + x169 + x182 + x184 + x19*(x183 + x77) + x46 + x129/n3) + 2.0*x12 + x130*x165);
result[10] = x156 - 3.5*x181 + x75*(x106*(x157 + x180 + x185 + x47*x58 + x93) + 3.5*x12 + x130*x154 + x153 + 94266.666666666657);
result[11] = x135 + x163 + x40 + x75*(x106*(x138 + x174 + x184 + x185) + x132 + x161);
result[12] = -7.0*x155*x41 + 24.5*x39 + x75*(58.201238327072673*T*x151 + x106*(-x139*x187 + x141 + x172 + x19*(x186 + x77) + x190 + x81*(x175 + x80) + x97 + x144/n4) + 7.0*x14);
result[13] = x156 - 3.5*x191 + x75*(x106*(x111 + x126 + x142 + x171 - x187*x82 + x19*(x175 + x186) + x190) + 3.5*x15 + x153 + x154*x159 + 15555.555555555555);
result[14] = -2.0*x191 + x40 + x75*(x106*(x173 + x184 + x189 + x19*(x32*x78 + x77) + x52 + x85 - x95 + x99 + x158/n5) + 2.0*x15 + x159*x165);
}
        
static void coder_d3gdn3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 3.5*n4;
    double x1 = 1.0*n5;
    double x2 = 1.0*n2;
    double x3 = 1.0*n1;
    double x4 = 1.0*n3;
    double x5 = x3 + x4;
    double x6 = x1 + x2 + x5;
    double x7 = x0 + x6;
    double x8 = (*endmember[0].mu0)(T, P);
    double x9 = n1*x8;
    double x10 = (*endmember[1].mu0)(T, P);
    double x11 = n2*x10;
    double x12 = (*endmember[2].mu0)(T, P);
    double x13 = n3*x12;
    double x14 = (*endmember[3].mu0)(T, P);
    double x15 = (*endmember[4].mu0)(T, P);
    double x16 = n5*x15;
    double x17 = n1 + n3;
    double x18 = n4 + n5;
    double x19 = n2 + x17 + x18;
    double x20 = 1.0/x19;
    double x21 = log(n2*x20);
    double x22 = log(n3*x20);
    double x23 = log(n5*x20);
    double x24 = log(n4*x20);
    double x25 = 1.0*n4;
    double x26 = 1.0*log(x17*x20);
    double x27 = 1.0*log(x18*x20);
    double x28 = n1 + n2 + n4;
    double x29 = 1.0*log(x20*x28);
    double x30 = 2.0*n2;
    double x31 = 2.0*n3;
    double x32 = 2.0*n5;
    double x33 = 2.0*n1 + x25 + x30 + x31 + x32;
    double x34 = 0.5*n4 + x6;
    double x35 = log(x20*x34);
    double x36 = x1*x23 + x17*x26 + x18*x27 + x2*x21 + x22*x4 + x25*(x24 - 0.69314718055994495) + x28*x29 + x33*x35;
    double x37 = 8.3144626181532395*T;
    double x38 = x36*x37;
    double x39 = (0.027777777777777773*n1*(444600.0*n3 + 728000.0*n4 + 437400.0*n5) + 12349.999999999998*n2*n3 + 0.027777777777777773*n3*(444600.0*n1 + 444600.0*n2 + 1696800.0*n4) + 0.0085352842554488641*n4*(2369250.0*n1 + 5522175.0*n3 + 911250.0*n5) + 0.027777777777777773*n5*(437400.0*n1 + 280000.0*n4) + x7*(n4*x14 + x11 + x13 + x16 + x38 + x9))/((x7)*(x7)*(x7)*(x7));
    double x40 = -6.0*x39;
    double x41 = pow(x7, -3);
    double x42 = 2.0*x35;
    double x43 = -x2*x20;
    double x44 = -x20*x4;
    double x45 = -x20*x25;
    double x46 = 1.0/x34;
    double x47 = 1.0*x20;
    double x48 = pow(x19, -2);
    double x49 = -x34*x48;
    double x50 = x47 + x49;
    double x51 = x46*x50;
    double x52 = x33*x51;
    double x53 = x19*x52;
    double x54 = x42 + x43 + x44 + x45 + x53;
    double x55 = x1 + x25;
    double x56 = -x20*x55;
    double x57 = -x1*x20;
    double x58 = x2 + x25 + x3;
    double x59 = 1.0/x28;
    double x60 = x20 - x28*x48;
    double x61 = x59*x60;
    double x62 = x58*x61;
    double x63 = x19*x62 + x29 + x57;
    double x64 = x56 + x63;
    double x65 = 1.0/x17;
    double x66 = -x17*x48 + x20;
    double x67 = x65*x66;
    double x68 = x5*x67;
    double x69 = x19*x68 + x26;
    double x70 = x54 + x64 + x69;
    double x71 = x37*x70;
    double x72 = x1*x15 + x10*x2 + x12*x4 + x14*x25 + x3*x8 + x38;
    double x73 = 24699.999999999996*n3 + x72;
    double x74 = x41*(40444.444444444438*n4 + 24299.999999999996*n5 + x7*(x71 + x8) + x73);
    double x75 = pow(x7, -2);
    double x76 = 16.628925236306479*T;
    double x77 = 4.0*x19;
    double x78 = 2.0*x48;
    double x79 = -x78;
    double x80 = pow(x19, -3);
    double x81 = 2*x80;
    double x82 = x34*x81;
    double x83 = x79 + x82;
    double x84 = x33*x46;
    double x85 = x83*x84;
    double x86 = 1.0*x19;
    double x87 = pow(x34, -2);
    double x88 = x50*x87;
    double x89 = x33*x88;
    double x90 = x19*x85 + x51*x77 - x86*x89;
    double x91 = 2.0*x61;
    double x92 = x19*x91;
    double x93 = -2*x48;
    double x94 = x28*x81 + x93;
    double x95 = x58*x59;
    double x96 = x94*x95;
    double x97 = x19*x96;
    double x98 = pow(x28, -2);
    double x99 = x60*x98;
    double x100 = x58*x99;
    double x101 = -x100*x19;
    double x102 = x2*x48;
    double x103 = x4*x48;
    double x104 = x1*x48;
    double x105 = x102 + x103 + x104 + x62;
    double x106 = x101 + x105 + x92 + x97;
    double x107 = x106 + x90;
    double x108 = x25*x48;
    double x109 = -x25;
    double x110 = -x1 + x109;
    double x111 = x108 - x110*x48 + x52;
    double x112 = 2.0*x67;
    double x113 = x17*x81 + x93;
    double x114 = x5*x65;
    double x115 = x113*x114;
    double x116 = pow(x17, -2);
    double x117 = x116*x66;
    double x118 = x117*x5;
    double x119 = x112*x19 + x115*x19 - x118*x19 + x68;
    double x120 = x111 + x119;
    double x121 = T*(x107 + x120);
    double x122 = 8.3144626181532395*x121;
    double x123 = x75*(x122*x7 + x70*x76 + 2.0*x8);
    double x124 = 1.0/x7;
    double x125 = x30*x80;
    double x126 = -x125;
    double x127 = 3.0*x19;
    double x128 = 6*x80;
    double x129 = pow(x19, -4);
    double x130 = 6*x129;
    double x131 = x114*x19;
    double x132 = x116*x5;
    double x133 = 2*x19;
    double x134 = x113*x127*x65 - x113*x132*x133 + 2*x115 - x117*x127 - 2*x118 + x131*(x128 - x130*x17) + x133*x5*x66/((x17)*(x17)*(x17)) + 3.0*x67;
    double x135 = x126 + x134;
    double x136 = -x55*x81;
    double x137 = x31*x80;
    double x138 = -x137;
    double x139 = 2.0*x80;
    double x140 = n4*x139;
    double x141 = -x140;
    double x142 = x32*x80;
    double x143 = -x142;
    double x144 = x138 + x141 + x143;
    double x145 = x136 + x144;
    double x146 = 6.0*x19;
    double x147 = x46*x83;
    double x148 = 2.0*x33;
    double x149 = 6.0*x80;
    double x150 = -x130*x34;
    double x151 = x19*x84;
    double x152 = x148*x19;
    double x153 = pow(x34, -3);
    double x154 = x153*x50;
    double x155 = x83*x87;
    double x156 = x146*x147 - x146*x88 - x148*x88 + x151*(x149 + x150) + x152*x154 - x152*x155 + 6.0*x51 + 2*x85;
    double x157 = x145 + x156;
    double x158 = x19*x95;
    double x159 = x58*x98;
    double x160 = -2*x100 + x127*x59*x94 - x127*x99 - x133*x159*x94 + x133*x58*x60/((x28)*(x28)*(x28)) + x158*(x128 - x130*x28) + 3.0*x61 + 2*x96;
    double x161 = x157 + x160;
    double x162 = x37*x7;
    double x163 = -x48;
    double x164 = -n1;
    double x165 = -n3 + x164;
    double x166 = x163 - x165*x81;
    double x167 = x114*x166;
    double x168 = x167*x19;
    double x169 = x168 + x68;
    double x170 = -2.0*x20;
    double x171 = x111 + x170;
    double x172 = x107 + x169 + x171;
    double x173 = x172*x76;
    double x174 = 1.0*x48;
    double x175 = x166*x19;
    double x176 = 4*x80;
    double x177 = 2*n1;
    double x178 = 2*n3;
    double x179 = 3*x129;
    double x180 = -x179*(x177 + x178);
    double x181 = x115 - x118 - x132*x175 + x167;
    double x182 = x112 + x131*(x176 + x180) + 2.0*x175*x65 + x181;
    double x183 = x126 + x182;
    double x184 = -n2*x48 + x20;
    double x185 = x184*x86;
    double x186 = -x20*x5;
    double x187 = x186 + x44;
    double x188 = x42 + x45 + x53;
    double x189 = x185 + x187 + x188 + 1.0*x21 + x64;
    double x190 = x189*x37;
    double x191 = x7*(x10 + x190) + x73;
    double x192 = 2.0*x41;
    double x193 = x191*x192;
    double x194 = x172*x37;
    double x195 = 1.0*x10 + x190;
    double x196 = x71 + 1.0*x8;
    double x197 = x194*x7 + x195 + x196;
    double x198 = 2.0*x75;
    double x199 = -x197*x198 + x40;
    double x200 = -1.0*x123 + 4.0*x74;
    double x201 = -n2 - n4 + x164;
    double x202 = x163 - x201*x81;
    double x203 = x202*x95;
    double x204 = x19*x203;
    double x205 = x105 + x204;
    double x206 = x205 + x90;
    double x207 = x119 + x171 + x206;
    double x208 = x207*x76;
    double x209 = x19*x202;
    double x210 = 2*n2;
    double x211 = 2*n4;
    double x212 = -x179*(x177 + x210 + x211);
    double x213 = -x100 - x159*x209 + x203 + x96;
    double x214 = x158*(x176 + x212) + 2.0*x209*x59 + x213 + x91;
    double x215 = -x20*x58;
    double x216 = -n3*x48 + x20;
    double x217 = x216*x86;
    double x218 = x188 + x215 + x217 + 1.0*x22 + x43 + x56 + x57 + x69;
    double x219 = x218*x37;
    double x220 = 24699.999999999996*n1 + 24699.999999999996*n2 + 94266.666666666657*n4 + x7*(x12 + x219) + x72;
    double x221 = x192*x220;
    double x222 = x207*x37;
    double x223 = 1.0*x12 + x219;
    double x224 = x223 + 24699.999999999996;
    double x225 = x196 + x222*x7 + x224;
    double x226 = -x198*x225;
    double x227 = x200 + x40;
    double x228 = -3.0*x20;
    double x229 = x108 + x48*x55 + x52;
    double x230 = x169 + x229;
    double x231 = 0.5*x19;
    double x232 = 0.5*x20 + x49;
    double x233 = x232*x46;
    double x234 = 2.0*x233;
    double x235 = x19*x234;
    double x236 = -1.5*x48 + x82;
    double x237 = x236*x84;
    double x238 = x19*x237 + x235 + x51*x86;
    double x239 = -x231*x89 + x238;
    double x240 = x106 + x228 + x230 + x239;
    double x241 = x110*x81;
    double x242 = x33*x86;
    double x243 = x242*x87;
    double x244 = -x236*x243;
    double x245 = x231*x33;
    double x246 = x236*x46;
    double x247 = x147*x86 + x151*(x150 + 5.0*x80) + x246*x77;
    double x248 = x154*x242 - x155*x245 + x237 + x247 + x85 - 1.5*x89;
    double x249 = -x127*x88 + x241 + x244 + x248 + 5.0*x51;
    double x250 = x144 + x249;
    double x251 = x160 + x250;
    double x252 = x183 + x78;
    double x253 = 29.100619163536336*T;
    double x254 = x240*x37;
    double x255 = -n4*x48 + x20;
    double x256 = x255*x86;
    double x257 = x232*x84;
    double x258 = 1.0/x18;
    double x259 = -x18*x48 + x20;
    double x260 = x258*x259;
    double x261 = x260*x55;
    double x262 = x19*x261 + x27;
    double x263 = x187 + x19*x257 + 1.0*x24 + x256 + x262 + 1.0*x35 + x43 + x63 - 0.69314718055994495;
    double x264 = x263*x37;
    double x265 = 1.0*x14 + x264;
    double x266 = x253*x70 + x254*x7 + x265 + 3.5*x8 + 40444.444444444438;
    double x267 = 40444.444444444438*n1 + 94266.666666666657*n3 + 15555.555555555555*n5 + x0*x14 + 3.5*x11 + 3.5*x13 + 3.5*x16 + x253*x36 + x7*(x14 + x264) + 3.5*x9;
    double x268 = x192*x267 - 21.0*x39;
    double x269 = x168 - 4.0*x20 + x229 + x68;
    double x270 = x206 + x269;
    double x271 = x270*x76;
    double x272 = x156 + x241;
    double x273 = x144 + x272;
    double x274 = x214 + x273;
    double x275 = -n5*x48 + x20;
    double x276 = x275*x86;
    double x277 = x186 + x215 + 1.0*x23 + x262 + x276 + x54;
    double x278 = x277*x37;
    double x279 = 24299.999999999996*n1 + 15555.555555555555*n4 + x7*(x15 + x278) + x72;
    double x280 = x192*x279;
    double x281 = x270*x37;
    double x282 = 1.0*x15 + x278;
    double x283 = x196 + x281*x7 + x282 + 24299.999999999996;
    double x284 = -x198*x283;
    double x285 = -x102;
    double x286 = 1.0/n2;
    double x287 = x101 + x92 + x97;
    double x288 = x104 + x90;
    double x289 = -x3;
    double x290 = x289 - x4;
    double x291 = -x290*x48;
    double x292 = x291 + x47;
    double x293 = x103 + x292;
    double x294 = x111 + x185*x286 + x19*(x125 + x79) + x285 + x287 + x288 + x293 + x62;
    double x295 = x294*x37;
    double x296 = x131*(x180 + x81) + 2*x167;
    double x297 = x126 + x296;
    double x298 = 3.0*x48;
    double x299 = x157 + x298;
    double x300 = 2.0*x74;
    double x301 = x191*x41;
    double x302 = 2.0*x10 + x189*x76 + x295*x7;
    double x303 = 1.0*x75;
    double x304 = 4.0*x301 - x302*x303;
    double x305 = x104 + x228 + x90;
    double x306 = -x174;
    double x307 = x103 + x19*(x125 + x306) + x285 + x62;
    double x308 = x204 + x305 + x307;
    double x309 = x111 + x308 + x48*x5;
    double x310 = x309*x37;
    double x311 = -x225*x303;
    double x312 = x195 + x224 + x310*x7;
    double x313 = x221 - x303*x312;
    double x314 = x300 + x40;
    double x315 = x193 - x197*x303 + x314;
    double x316 = x229 + x291;
    double x317 = x104 + x239;
    double x318 = x170 + x287 + x307 + x316 + x317;
    double x319 = x318*x37;
    double x320 = 4.0*x48;
    double x321 = x297 + x320;
    double x322 = 3.5*x75;
    double x323 = -x266*x303 + x268 + 7.0*x74;
    double x324 = 3.5*x10 + x189*x253 + x265 + x319*x7;
    double x325 = 7.0*x301 - x303*x324;
    double x326 = x308 + x316;
    double x327 = x326*x37;
    double x328 = x195 + x282 + x327*x7;
    double x329 = -x303*x328;
    double x330 = x280 - x283*x303;
    double x331 = -x103;
    double x332 = 1.0/n3;
    double x333 = x109 - x2 + x289;
    double x334 = x102 - x333*x48;
    double x335 = x120 + x19*(x137 + x79) + x217*x332 + x288 + x331 + x334 + x47;
    double x336 = x335*x37;
    double x337 = 2*x203;
    double x338 = x158*(x212 + x81);
    double x339 = x126 + x337 + x338;
    double x340 = x157 + x339;
    double x341 = x220*x41;
    double x342 = 2.0*x12 + x218*x76 + x336*x7;
    double x343 = -x303*x342 + 4.0*x341;
    double x344 = x19*(x137 + x306) + x331;
    double x345 = x102 + x269 + x317 + x344 + x48*x58;
    double x346 = x345*x37;
    double x347 = x182 + x320;
    double x348 = x158*(x130*x201 + x176) + x202*x59*x86 + x213 + 1.0*x61;
    double x349 = 3.5*x12 + x218*x253 + x265 + x346*x7 + 94266.666666666657;
    double x350 = -x303*x349 + 7.0*x341;
    double x351 = x230 + x305 + x334 + x344;
    double x352 = x351*x37;
    double x353 = x223 + x282 + x352*x7;
    double x354 = -x303*x353;
    double x355 = 58.201238327072673*T;
    double x356 = 1.0/n4;
    double x357 = x306 + x82;
    double x358 = x357*x84;
    double x359 = x232*x87;
    double x360 = x33*x359;
    double x361 = 2.0*x19;
    double x362 = x19*x55;
    double x363 = x18*x81 + x93;
    double x364 = x258*x363;
    double x365 = pow(x18, -2);
    double x366 = x259*x365;
    double x367 = x260*x361 + x261 + x362*x364 - x362*x366;
    double x368 = -x108 + x257 + x367;
    double x369 = x106 + x19*x358 + x19*(x140 + x79) - x231*x360 + x235 + x256*x356 + x292 + x368;
    double x370 = x369*x37;
    double x371 = x297 + 5.0*x48;
    double x372 = x359*x86;
    double x373 = 2*x237 - x86*x88;
    double x374 = x357*x46;
    double x375 = 4.0*x80;
    double x376 = x151*(x150 + x375) + x246*x361 + x361*x374;
    double x377 = x154*x245 - x372 + x373 + x376 + 2.0*x51 - 1.0*x89;
    double x378 = x234 + x244;
    double x379 = x145 + x378;
    double x380 = x160 + x377 + x379;
    double x381 = 7.0*x75;
    double x382 = x267*x41;
    double x383 = 7.0*x14 + x263*x355 + x370*x7;
    double x384 = -x303*x383 + 14.0*x382 - 73.5*x39;
    double x385 = x170 + x19*(x140 + x306) + x205 + x238 + x291 - x33*x372 + x368;
    double x386 = x37*x385;
    double x387 = x248 - x359*x361 - x361*x88 + 3.0*x51;
    double x388 = x371 + x387;
    double x389 = x214 + x379;
    double x390 = x279*x41;
    double x391 = 3.5*x15 + x253*x277 + x265 + x386*x7 + 15555.555555555555;
    double x392 = -x303*x391 + 7.0*x390;
    double x393 = 1.0/n5;
    double x394 = -x104 + x108 + x19*(x142 + x79) + x276*x393 + x293 + x334 + x367 + x52 + x90;
    double x395 = x37*x394;
    double x396 = x296 + 6.0*x48;
    double x397 = 2.0*x15 + x277*x76 + x395*x7;
    double x398 = -x303*x397 + 4.0*x390;
    double x399 = 3.0*x75;
    double x400 = 24.943387854459719*T;
    double x401 = 1.0*x184*x286;
    double x402 = 6.0*x129;
    double x403 = -n2*x402;
    double x404 = x210*x80;
    double x405 = x286*x86;
    double x406 = -x5*x81;
    double x407 = n2*x375;
    double x408 = x157 + x407;
    double x409 = x406 + x408;
    double x410 = -x320;
    double x411 = x160 + x410;
    double x412 = x309*x76;
    double x413 = x19*(x375 + x403) + x306 + x401 + x405*(x163 + x404);
    double x414 = x214 + x413;
    double x415 = -x198*x312;
    double x416 = x304 + x40;
    double x417 = x406 + x407;
    double x418 = x250 + x417;
    double x419 = x326*x76;
    double x420 = x273 + x417;
    double x421 = -x198*x328;
    double x422 = x19*(x139 + x403);
    double x423 = x320 + x422;
    double x424 = x337 + x338 + x423;
    double x425 = x162*(x409 + x424);
    double x426 = x193 + x40;
    double x427 = x268 + x325;
    double x428 = x298 + x417 + x422;
    double x429 = -n3*x402;
    double x430 = x178*x80;
    double x431 = x332*x86;
    double x432 = n3*x375 + x141 + x143;
    double x433 = 1.0*x216*x332 + x432;
    double x434 = -x58*x81;
    double x435 = x156 + x434;
    double x436 = x410 + x435;
    double x437 = x182 + x19*(x375 + x429) + x431*(x163 + x430) + x433;
    double x438 = x351*x76;
    double x439 = x126 + x434;
    double x440 = -x198*x353 + x40;
    double x441 = x136 + x19*(x139 + x429) + x432;
    double x442 = x378 + x441;
    double x443 = -n4*x402;
    double x444 = x211*x80;
    double x445 = x356*x86;
    double x446 = x153*x232;
    double x447 = n4*x375 + x143;
    double x448 = 2*x55;
    double x449 = 2*x362;
    double x450 = x127*x364 - x127*x366 + x138 + x258*x362*(x128 - x130*x18) + 3.0*x260 - x363*x365*x449 + x364*x448 - x366*x448 + x406 + x259*x449/((x18)*(x18)*(x18));
    double x451 = x126 + x450;
    double x452 = -x243*x357 + 1.0*x255*x356 + x451;
    double x453 = x236*x87;
    double x454 = 4.0*x233 + x447;

result[0] = -3.0*x123 + x124*(24.943387854459719*x121 + x162*(x135 + x161)) + x40 + 6.0*x74;
result[1] = x124*(x122 + x162*(x161 + x174 + x183) + x173) + x193 + x199 + x200;
result[2] = x124*(x122 + x162*(x135 + x157 + x174 + x214) + x208) + x221 + x226 + x227;
result[3] = -3.5*x123 + x124*(29.100619163536336*x121 + x162*(x251 + x252) + x240*x76) - x198*x266 + x268 + 14.0*x74;
result[4] = x124*(x122 + x162*(x252 + x274) + x271) + x227 + x280 + x284;
result[5] = x124*(x162*(x160 + x297 + x299) + x173 + x295) + x199 + x300 + x304;
result[6] = x124*(x162*(x126 + x131*(x130*x165 + x176) + x166*x65*x86 + x181 + x214 + x299 + 1.0*x67) + x194 + x222 + x310) + x311 + x313 + x315;
result[7] = x124*(x162*(x251 + x321) + x172*x253 + x254 + x319) - x197*x322 + x323 + x325;
result[8] = x124*(x162*(x274 + x321) + x194 + x281 + x327) + x315 + x329 + x330;
result[9] = x124*(x162*(x134 + x298 + x340) + x208 + x336) + x226 + x314 + x343;
result[10] = x124*(x162*(x126 + x250 + x347 + x348) + x207*x253 + x254 + x346) - x225*x322 + x323 + x350;
result[11] = x124*(x162*(x273 + x339 + x347) + x222 + x281 + x352) + x221 + x311 + x314 + x330 + x354;
result[12] = x124*(x162*(x371 + x380) + x240*x355 + x370) - x266*x381 + x384 + 24.5*x74;
result[13] = x124*(x162*(x388 + x389) + x253*x270 + x254 + x386) - x283*x322 + x323 + x392;
result[14] = x124*(x162*(x340 + x396) + x271 + x395) + x284 + x314 + x398;
result[15] = x124*(x162*(x19*(x149 + x403) + x401 + x405*(x404 + x93) + x409 + x411 - x185/((n2)*(n2))) + x294*x400) + 6.0*x301 - x302*x399 + x40;
result[16] = x124*(x162*(x290*x81 + x408 + x414) + x295 + x412) + x221 + x415 + x416;
result[17] = x124*(x162*(x160 + x413 + x418) + x253*x294 + x318*x76) - x198*x324 + x268 + 14.0*x301 - x302*x322;
result[18] = x124*(x162*(x414 + x420) + x295 + x419) + x280 + x416 + x421;
result[19] = x124*(x336 + x412 + x425) + x343 + x415 + x426;
result[20] = x124*(x162*(x348 + x418 + x423) + x253*x309 + x319 + x346) - x312*x322 + x350 + x427;
result[21] = x124*(x162*(x420 + x424) + x310 + x327 + x352) + x280 + x313 + x329 + x354 + x426;
result[22] = x124*(x162*(x380 + x428) + x318*x355 + x370) + 24.5*x301 - x324*x381 + x384;
result[23] = x124*(x162*(x387 + x389 + x428) + x253*x326 + x319 + x386) - x322*x328 + x392 + x427;
result[24] = x124*(x395 + x419 + x425) + x398 + x421 + x426;
result[25] = x124*(x162*(x135 + x136 + x19*(x149 + x429) + x431*(x430 + x93) + x433 + x436 - x217/((n3)*(n3))) + x335*x400) + 6.0*x341 - x342*x399 + x40;
result[26] = x124*(x162*(x126 + x249 + x333*x81 + x437) + x253*x335 + x345*x76) - x198*x349 + x268 - x322*x342 + 14.0*x341;
result[27] = x124*(x162*(x272 + x306 + x437 + x439) + x336 + x438) + x280 + x343 + x440;
result[28] = x124*(x162*(x377 + x396 + x439 + x442) + x345*x355 + x370) + 24.5*x341 - x349*x381 + x384;
result[29] = x124*(x162*(x388 + x434 + x442) + x253*x351 + x346 + x386) + x268 - x322*x353 + x350 + x392;
result[30] = x124*(x162*(x321 + x435 + x441) + x395 + x438) + x221 + x398 + x440;
result[31] = x124*(87.301857490609009*T*x369 + x162*(x127*x374 + x151*(x150 + 3.0*x80) - 1.5*x19*x359 + x19*(x149 + x443) + 3.0*x233 + x245*x446 + 2*x358 - 1.0*x360 + x411 + x445*(x444 + x93) + x447 + x452 - x256/((n4)*(n4)))) + 73.5*x382 - 10.5*x383*x75 - 257.25*x39;
result[32] = x124*(x162*(-x127*x359 + x19*(x375 + x443) + x214 + x237 + x242*x446 - x245*x453 + x358 - 1.5*x360 + x376 + x445*(x163 + x444) + x452 + x454 + x79) + x355*x385 + x370) - x381*x391 + x384 + 24.5*x390;
result[33] = x124*(x162*(-x148*x359 + x152*x446 - x152*x453 + x19*(x139 + x443) + x247 + x339 - x359*x77 + x373 + x450 + x454 + 1.0*x51 + x78) + x253*x394 + x385*x76) - x198*x391 + x268 - x322*x397 + 14.0*x390;
result[34] = x124*(x162*(n5*x375 + x141 + x19*(-n5*x402 + x149) + 1.0*x275*x393 + x393*x86*(n5*x81 + x93) + x436 + x451 - x276/((n5)*(n5))) + x394*x400) + 6.0*x390 - x397*x399 + x40;
}
        
static double coder_dgdt(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = n1 + n3;
    double x1 = n4 + n5;
    double x2 = 1.0/(n2 + x0 + x1);
    double x3 = n1 + n2 + n4;

result = n1*(*endmember[0].dmu0dT)(T, P) + 8.3144626181532395*n2*log(n2*x2) + n2*(*endmember[1].dmu0dT)(T, P) + 8.3144626181532395*n3*log(n3*x2) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*(log(n4*x2) - 0.69314718055994495) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626181532395*n5*log(n5*x2) + n5*(*endmember[4].dmu0dT)(T, P) + 8.3144626181532395*x0*log(x0*x2) + 8.3144626181532395*x1*log(x1*x2) + 8.3144626181532395*x3*log(x2*x3) + 8.3144626181532395*(2.0*n1 + 2.0*n2 + 2.0*n3 + 1.0*n4 + 2.0*n5)*log(x2*(1.0*n1 + 1.0*n2 + 1.0*n3 + 0.5*n4 + 1.0*n5));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n3;
    double x1 = n4 + n5;
    double x2 = n2 + x0 + x1;
    double x3 = 1.0/x2;
    double x4 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 0.5*n4 + 1.0*n5;
    double x5 = log(x3*x4);
    double x6 = 16.628925236306479*x5;
    double x7 = 8.3144626181532395*n2;
    double x8 = -x3*x7;
    double x9 = 8.3144626181532395*n3;
    double x10 = -x3*x9;
    double x11 = 8.3144626181532395*n4;
    double x12 = -x11*x3;
    double x13 = pow(x2, -2);
    double x14 = -x13*x4;
    double x15 = x2*(16.628925236306479*n1 + 16.628925236306479*n2 + 16.628925236306479*n3 + 16.628925236306479*n5 + x11)/x4;
    double x16 = x15*(x14 + 1.0*x3);
    double x17 = x10 + x12 + x16 + x6 + x8;
    double x18 = 8.3144626181532395*n5;
    double x19 = x11 + x18;
    double x20 = -x19*x3;
    double x21 = n1 + n2 + n4;
    double x22 = -x18*x3;
    double x23 = 8.3144626181532395*n1;
    double x24 = x11 + x23 + x7;
    double x25 = x2*x24*(-x13*x21 + x3)/x21 + x22 + 8.3144626181532395*log(x21*x3);
    double x26 = x20 + x25;
    double x27 = x23 + x9;
    double x28 = 8.3144626181532395*log(x0*x3) + x2*x27*(-x0*x13 + x3)/x0;
    double x29 = 8.3144626181532395*x2;
    double x30 = -x27*x3;
    double x31 = x10 + x30;
    double x32 = x12 + x16 + x6;
    double x33 = -x24*x3;
    double x34 = 8.3144626181532395*log(x1*x3) + x19*x2*(-x1*x13 + x3)/x1;

result[0] = x17 + x26 + x28 + (*endmember[0].dmu0dT)(T, P);
result[1] = x26 + x29*(-n2*x13 + x3) + x31 + x32 + 8.3144626181532395*log(n2*x3) + (*endmember[1].dmu0dT)(T, P);
result[2] = x20 + x22 + x28 + x29*(-n3*x13 + x3) + x32 + x33 + x8 + 8.3144626181532395*log(n3*x3) + (*endmember[2].dmu0dT)(T, P);
result[3] = x15*(x14 + 0.5*x3) + x25 + x29*(-n4*x13 + x3) + x31 + x34 + 8.3144626181532395*x5 + x8 + 8.3144626181532395*log(n4*x3) + (*endmember[3].dmu0dT)(T, P) - 5.7631463216439762;
result[4] = x17 + x29*(-n5*x13 + x3) + x30 + x33 + x34 + 8.3144626181532395*log(n5*x3) + (*endmember[4].dmu0dT)(T, P);
}
        
static void coder_d3gdn2dt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 0.5*n4 + 1.0*n5;
    double x1 = 1.0/x0;
    double x2 = n1 + n3;
    double x3 = n4 + n5;
    double x4 = n2 + x2 + x3;
    double x5 = 1.0/x4;
    double x6 = pow(x4, -2);
    double x7 = -x0*x6;
    double x8 = 1.0*x5 + x7;
    double x9 = x4*x8;
    double x10 = x1*x9;
    double x11 = pow(x4, -3);
    double x12 = 2*x11;
    double x13 = x0*x12;
    double x14 = 16.628925236306479*n2;
    double x15 = 16.628925236306479*n3;
    double x16 = 8.3144626181532395*n4;
    double x17 = 16.628925236306479*n5;
    double x18 = 16.628925236306479*n1 + x14 + x15 + x16 + x17;
    double x19 = x1*x18;
    double x20 = x19*x4;
    double x21 = x18/((x0)*(x0));
    double x22 = x21*x9;
    double x23 = 33.257850472612958*x10 + x20*(x13 - 2.0*x6) - 1.0*x22;
    double x24 = n1 + n2 + n4;
    double x25 = 1.0/x24;
    double x26 = -x24*x6 + x5;
    double x27 = x25*x26;
    double x28 = 16.628925236306479*x4;
    double x29 = x27*x28;
    double x30 = -2*x6;
    double x31 = 8.3144626181532395*n1;
    double x32 = 8.3144626181532395*n2;
    double x33 = x16 + x31 + x32;
    double x34 = x33*x4;
    double x35 = x25*x34;
    double x36 = x35*(x12*x24 + x30);
    double x37 = -x26*x34/((x24)*(x24));
    double x38 = x32*x6;
    double x39 = 8.3144626181532395*n3;
    double x40 = x39*x6;
    double x41 = 8.3144626181532395*n5;
    double x42 = x41*x6;
    double x43 = x27*x33;
    double x44 = x38 + x40 + x42 + x43;
    double x45 = x29 + x36 + x37 + x44;
    double x46 = x23 + x45;
    double x47 = x16*x6;
    double x48 = x19*x8;
    double x49 = -x16;
    double x50 = x47 + x48 - x6*(-x41 + x49);
    double x51 = x31 + x39;
    double x52 = 1.0/x2;
    double x53 = -x2*x6 + x5;
    double x54 = x52*x53;
    double x55 = x51*x54;
    double x56 = x4*x51;
    double x57 = x52*x56;
    double x58 = x28*x54 + x55 + x57*(x12*x2 + x30) - x53*x56/((x2)*(x2));
    double x59 = x50 + x58;
    double x60 = -x6;
    double x61 = -n1;
    double x62 = x57*(-x12*(-n3 + x61) + x60);
    double x63 = x55 + x62;
    double x64 = -16.628925236306479*x5;
    double x65 = x50 + x64;
    double x66 = x35*(-x12*(-n2 - n4 + x61) + x60);
    double x67 = x44 + x66;
    double x68 = x23 + x67;
    double x69 = -24.943387854459719*x5;
    double x70 = x16 + x41;
    double x71 = x47 + x48 + x6*x70;
    double x72 = x63 + x71;
    double x73 = 0.5*x5 + x7;
    double x74 = x1*x28*x73;
    double x75 = 8.3144626181532395*x10 + x20*(x13 - 1.5*x6) + x74;
    double x76 = -0.5*x22 + x75;
    double x77 = -33.257850472612958*x5 + x55 + x62 + x71;
    double x78 = -x38;
    double x79 = -16.628925236306479*x6;
    double x80 = x11*x14;
    double x81 = 8.3144626181532395*x4;
    double x82 = x29 + x36 + x37;
    double x83 = x23 + x42;
    double x84 = -x31;
    double x85 = -x6*(-x39 + x84);
    double x86 = 8.3144626181532395*x5;
    double x87 = x85 + x86;
    double x88 = x40 + x87;
    double x89 = x23 + x42 + x69;
    double x90 = -8.3144626181532395*x6;
    double x91 = x4*(x80 + x90) + x40 + x43 + x78;
    double x92 = x66 + x89 + x91;
    double x93 = x71 + x85;
    double x94 = x42 + x76;
    double x95 = -x40;
    double x96 = x11*x15;
    double x97 = x38 - x6*(-x32 + x49 + x84);
    double x98 = x4*(x90 + x96) + x95;
    double x99 = 16.628925236306479*n4*x11;
    double x100 = x21*x4*x73;
    double x101 = 1.0/x3;
    double x102 = -x3*x6 + x5;
    double x103 = x101*x102;
    double x104 = x4*x70;
    double x105 = x101*x104*(x12*x3 + x30) - x102*x104/((x3)*(x3)) + x103*x28 + x103*x70;
    double x106 = x105 + x19*x73 - x47;

result[0] = x46 + x59;
result[1] = x46 + x63 + x65;
result[2] = x58 + x65 + x68;
result[3] = x45 + x69 + x72 + x76;
result[4] = x68 + x77;
result[5] = x4*(x79 + x80) + x43 + x50 + x78 + x82 + x83 + x88 + x81*(-n2*x6 + x5)/n2;
result[6] = x50 + x51*x6 + x92;
result[7] = x64 + x82 + x91 + x93 + x94;
result[8] = x92 + x93;
result[9] = x4*(x79 + x96) + x59 + x83 + x86 + x95 + x97 + x81*(-n3*x6 + x5)/n3;
result[10] = x33*x6 + x38 + x77 + x94 + x98;
result[11] = x72 + x89 + x97 + x98;
result[12] = -0.5*x100 + x106 + x20*(x13 - 1.0*x6) + x4*(x79 + x99) + x45 + x74 + x87 + x81*(-n4*x6 + x5)/n4;
result[13] = -1.0*x100 + x106 + x4*(x90 + x99) + x64 + x67 + x75 + x85;
result[14] = x105 + x23 + x4*(x11*x17 + x79) - x42 + x47 + x48 + x88 + x97 + x81*(-n5*x6 + x5)/n5;
}
        
static void coder_d4gdn3dt(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n3;
    double x1 = n4 + n5;
    double x2 = n2 + x0 + x1;
    double x3 = pow(x2, -3);
    double x4 = 16.628925236306479*n2;
    double x5 = -x3*x4;
    double x6 = 1.0/x0;
    double x7 = 1.0/x2;
    double x8 = pow(x2, -2);
    double x9 = -x0*x8 + x7;
    double x10 = x6*x9;
    double x11 = 8.3144626181532395*n1;
    double x12 = 8.3144626181532395*n3;
    double x13 = x11 + x12;
    double x14 = pow(x0, -2);
    double x15 = x14*x9;
    double x16 = x13*x15;
    double x17 = -2*x8;
    double x18 = 2*x3;
    double x19 = x0*x18 + x17;
    double x20 = x13*x6;
    double x21 = x19*x20;
    double x22 = 24.943387854459719*x2;
    double x23 = 6*x3;
    double x24 = pow(x2, -4);
    double x25 = 6*x24;
    double x26 = x2*x20;
    double x27 = x13*x14;
    double x28 = 2*x2;
    double x29 = 24.943387854459719*x10 - x15*x22 - 2*x16 + x19*x22*x6 - x19*x27*x28 + 2*x21 + x26*(-x0*x25 + x23) + x13*x28*x9/((x0)*(x0)*(x0));
    double x30 = x29 + x5;
    double x31 = 8.3144626181532395*n4;
    double x32 = 8.3144626181532395*n5;
    double x33 = x31 + x32;
    double x34 = -x18*x33;
    double x35 = 16.628925236306479*n3;
    double x36 = -x3*x35;
    double x37 = 16.628925236306479*x3;
    double x38 = -n4*x37;
    double x39 = 16.628925236306479*n5;
    double x40 = -x3*x39;
    double x41 = x36 + x38 + x40;
    double x42 = x34 + x41;
    double x43 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 0.5*n4 + 1.0*n5;
    double x44 = 1.0/x43;
    double x45 = -x43*x8;
    double x46 = x45 + 1.0*x7;
    double x47 = x44*x46;
    double x48 = x18*x43;
    double x49 = x48 - 2.0*x8;
    double x50 = 16.628925236306479*n1 + x31 + x35 + x39 + x4;
    double x51 = x44*x50;
    double x52 = x49*x51;
    double x53 = 49.886775708919437*x2;
    double x54 = x44*x49;
    double x55 = pow(x43, -2);
    double x56 = x46*x55;
    double x57 = 2.0*x50;
    double x58 = -x25*x43;
    double x59 = x2*x51;
    double x60 = x2*x57;
    double x61 = pow(x43, -3);
    double x62 = x46*x61;
    double x63 = x49*x55;
    double x64 = 49.886775708919437*x47 + 2*x52 + x53*x54 - x53*x56 - x56*x57 + x59*(6.0*x3 + x58) + x60*x62 - x60*x63;
    double x65 = x42 + x64;
    double x66 = n1 + n2 + n4;
    double x67 = 1.0/x66;
    double x68 = -x66*x8 + x7;
    double x69 = x67*x68;
    double x70 = 8.3144626181532395*n2;
    double x71 = x11 + x31 + x70;
    double x72 = pow(x66, -2);
    double x73 = x68*x72;
    double x74 = x71*x73;
    double x75 = x17 + x18*x66;
    double x76 = x67*x71;
    double x77 = x75*x76;
    double x78 = x2*x76;
    double x79 = x71*x72;
    double x80 = x22*x67*x75 - x22*x73 - x28*x75*x79 + x28*x68*x71/((x66)*(x66)*(x66)) + 24.943387854459719*x69 - 2*x74 + 2*x77 + x78*(x23 - x25*x66);
    double x81 = x65 + x80;
    double x82 = 8.3144626181532395*x8;
    double x83 = -x8;
    double x84 = -n1;
    double x85 = -n3 + x84;
    double x86 = -x18*x85 + x83;
    double x87 = x2*x86;
    double x88 = x6*x87;
    double x89 = 4*x3;
    double x90 = 2*n1;
    double x91 = 2*n3;
    double x92 = 3*x24;
    double x93 = -x92*(x90 + x91);
    double x94 = x20*x86;
    double x95 = -x16 + x21 - x27*x87 + x94;
    double x96 = 16.628925236306479*x10 + x26*(x89 + x93) + 16.628925236306479*x88 + x95;
    double x97 = x5 + x96;
    double x98 = -n2 - n4 + x84;
    double x99 = -x18*x98 + x83;
    double x100 = x2*x99;
    double x101 = x100*x67;
    double x102 = 2*n2;
    double x103 = 2*n4;
    double x104 = -x92*(x102 + x103 + x90);
    double x105 = x76*x99;
    double x106 = -x100*x79 + x105 - x74 + x77;
    double x107 = 16.628925236306479*x101 + x106 + 16.628925236306479*x69 + x78*(x104 + x89);
    double x108 = -x31;
    double x109 = x18*(x108 - x32);
    double x110 = x48 - 1.5*x8;
    double x111 = x2*x50;
    double x112 = 1.0*x111;
    double x113 = x112*x55;
    double x114 = -x110*x113;
    double x115 = x110*x51;
    double x116 = x50*x56;
    double x117 = 0.5*x111;
    double x118 = 8.3144626181532395*x2;
    double x119 = x2*x44;
    double x120 = x110*x119;
    double x121 = x118*x54 + 33.257850472612958*x120 + x59*(5.0*x3 + x58);
    double x122 = x112*x62 + x115 - 1.5*x116 - x117*x63 + x121 + x52;
    double x123 = x109 + x114 + x122 - x22*x56 + 41.572313090766201*x47;
    double x124 = x123 + x41;
    double x125 = x124 + x80;
    double x126 = 16.628925236306479*x8;
    double x127 = x126 + x97;
    double x128 = x109 + x64;
    double x129 = x128 + x41;
    double x130 = x107 + x129;
    double x131 = x26*(x18 + x93) + 2*x94;
    double x132 = x131 + x5;
    double x133 = 24.943387854459719*x8;
    double x134 = x133 + x65;
    double x135 = 33.257850472612958*x8;
    double x136 = x132 + x135;
    double x137 = 2*x105;
    double x138 = x78*(x104 + x18);
    double x139 = x137 + x138 + x5;
    double x140 = x139 + x65;
    double x141 = x135 + x96;
    double x142 = 8.3144626181532395*x101 + x106 + 8.3144626181532395*x69 + x78*(x25*x98 + x89);
    double x143 = x132 + 41.572313090766201*x8;
    double x144 = x45 + 0.5*x7;
    double x145 = x144*x55;
    double x146 = 2*x115 - x118*x56;
    double x147 = x48 - 1.0*x8;
    double x148 = 16.628925236306479*x119*x147 + 16.628925236306479*x120 + x59*(4.0*x3 + x58);
    double x149 = -1.0*x116 + x117*x62 - x118*x145 + x146 + x148 + 16.628925236306479*x47;
    double x150 = x144*x44;
    double x151 = x114 + 16.628925236306479*x150;
    double x152 = x151 + x42;
    double x153 = x149 + x152 + x80;
    double x154 = 16.628925236306479*x2;
    double x155 = x122 - x145*x154 - x154*x56 + 24.943387854459719*x47;
    double x156 = x143 + x155;
    double x157 = x107 + x152;
    double x158 = x131 + 49.886775708919444*x8;
    double x159 = -n2*x8 + x7;
    double x160 = 8.3144626181532395/n2;
    double x161 = x159*x160;
    double x162 = 49.886775708919437*x3;
    double x163 = 49.886775708919437*x24;
    double x164 = -n2*x163;
    double x165 = x102*x3;
    double x166 = x160*x2;
    double x167 = -x13*x18;
    double x168 = 33.257850472612958*x3;
    double x169 = n2*x168;
    double x170 = x169 + x65;
    double x171 = x167 + x170;
    double x172 = -x135;
    double x173 = x172 + x80;
    double x174 = -x11;
    double x175 = -x82;
    double x176 = x161 + x166*(x165 + x83) + x175 + x2*(x164 + x168);
    double x177 = x107 + x176;
    double x178 = x167 + x169;
    double x179 = x124 + x178;
    double x180 = x129 + x178;
    double x181 = x2*(x164 + x37);
    double x182 = x135 + x181;
    double x183 = x137 + x138 + x182;
    double x184 = x171 + x183;
    double x185 = x133 + x178 + x181;
    double x186 = -n3*x163;
    double x187 = x3*x91;
    double x188 = 8.3144626181532395/n3;
    double x189 = x188*x2;
    double x190 = -n3*x8 + x7;
    double x191 = n3*x168 + x38 + x40;
    double x192 = x188*x190 + x191;
    double x193 = -x18*x71;
    double x194 = x193 + x64;
    double x195 = x172 + x194;
    double x196 = x189*(x187 + x83) + x192 + x2*(x168 + x186) + x96;
    double x197 = x193 + x5;
    double x198 = x191 + x2*(x186 + x37) + x34;
    double x199 = x151 + x198;
    double x200 = -n4*x163;
    double x201 = x147*x51;
    double x202 = x103*x3;
    double x203 = 8.3144626181532395/n4;
    double x204 = x2*x203;
    double x205 = x145*x50;
    double x206 = -n4*x8 + x7;
    double x207 = x145*x2;
    double x208 = x144*x61;
    double x209 = n4*x168 + x40;
    double x210 = -x1*x8 + x7;
    double x211 = 1.0/x1;
    double x212 = 24.943387854459719*x211;
    double x213 = 2*x33;
    double x214 = pow(x1, -2);
    double x215 = x210*x214;
    double x216 = x1*x18 + x17;
    double x217 = x213*x216;
    double x218 = x167 + x2*x211*x33*(-x1*x25 + x23) + x2*x212*x216 - x2*x214*x217 + x210*x212 + x211*x217 - x213*x215 - x215*x22 + x36 + x2*x210*x213/((x1)*(x1)*(x1));
    double x219 = x218 + x5;
    double x220 = -x113*x147 + x203*x206 + x219;
    double x221 = x110*x55;
    double x222 = 33.257850472612958*x150 + x209;
    double x223 = -n5*x8 + x7;
    double x224 = 8.3144626181532395/n5;

result[0] = x30 + x81;
result[1] = x81 + x82 + x97;
result[2] = x107 + x30 + x65 + x82;
result[3] = x125 + x127;
result[4] = x127 + x130;
result[5] = x132 + x134 + x80;
result[6] = 8.3144626181532395*x10 + x107 + x134 + x26*(x25*x85 + x89) + x5 + 8.3144626181532395*x88 + x95;
result[7] = x125 + x136;
result[8] = x130 + x136;
result[9] = x133 + x140 + x29;
result[10] = x124 + x141 + x142 + x5;
result[11] = x129 + x139 + x141;
result[12] = x143 + x153;
result[13] = x156 + x157;
result[14] = x140 + x158;
result[15] = x161 + x166*(x165 + x17) + x171 + x173 + x2*(x162 + x164) - x118*x159/((n2)*(n2));
result[16] = x170 + x177 + x18*(-x12 + x174);
result[17] = x176 + x179 + x80;
result[18] = x177 + x180;
result[19] = x184;
result[20] = x142 + x179 + x182;
result[21] = x180 + x183;
result[22] = x153 + x185;
result[23] = x155 + x157 + x185;
result[24] = x184;
result[25] = x189*(x17 + x187) + x192 + x195 + x2*(x162 + x186) + x30 + x34 - x118*x190/((n3)*(n3));
result[26] = x123 + x18*(x108 + x174 - x70) + x196 + x5;
result[27] = x128 + x175 + x196 + x197;
result[28] = x149 + x158 + x197 + x199;
result[29] = x156 + x193 + x199;
result[30] = x136 + x194 + x198;
result[31] = x117*x208 + x147*x22*x44 + 24.943387854459719*x150 + x173 + x2*(x162 + x200) + 2*x201 + x204*(x17 + x202) - 1.0*x205 - 12.471693927229859*x207 + x209 + x220 + x59*(3.0*x3 + x58) - x118*x206/((n4)*(n4));
result[32] = x107 + x112*x208 + x115 - x117*x221 - x126 - x145*x22 + x148 + x2*(x168 + x200) + x201 + x204*(x202 + x83) - 1.5*x205 + x220 + x222;
result[33] = x121 + x126 + x139 - x145*x57 + x146 + x2*(x200 + x37) - 33.257850472612958*x207 + x208*x60 + x218 - x221*x60 + x222 + 8.3144626181532395*x47;
result[34] = n5*x168 + x195 + x2*x224*(n5*x18 + x17) + x2*(-n5*x163 + x162) + x219 + x223*x224 + x38 - x118*x223/((n5)*(n5));
}
        
static double coder_dgdp(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P) + n4*(*endmember[3].dmu0dP)(T, P) + n5*(*endmember[4].dmu0dP)(T, P);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].dmu0dP)(T, P);
result[1] = (*endmember[1].dmu0dP)(T, P);
result[2] = (*endmember[2].dmu0dP)(T, P);
result[3] = (*endmember[3].dmu0dP)(T, P);
result[4] = (*endmember[4].dmu0dP)(T, P);
}
        
static void coder_d3gdn2dp(double T, double P, double n[5], double result[5]) {
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
        
static void coder_d4gdn3dp(double T, double P, double n[5], double result[5]) {
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
        
static double coder_d2gdt2(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P) + n5*(*endmember[4].d2mu0dT2)(T, P);
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = (*endmember[1].d2mu0dT2)(T, P);
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
    

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P) + n5*(*endmember[4].d3mu0dT3)(T, P);
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];


result[0] = (*endmember[0].d3mu0dT3)(T, P);
result[1] = (*endmember[1].d3mu0dT3)(T, P);
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

