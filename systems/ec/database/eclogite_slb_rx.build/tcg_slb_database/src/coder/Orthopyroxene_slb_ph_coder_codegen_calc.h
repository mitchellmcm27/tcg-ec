#include <math.h>


static double coder_g(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = 1.0/x1;
    double x3 = n1 + n4;

result = 1.0*x2*(16050.0*n1*n4 + 24000.0*n3*n4 + 0.125*n4*(128400.0*n1 + 192000.0*n3) + 1.0*x1*(8.3144626181532395*T*(2.0*n2*log(n2*x2) + 1.0*n3*log(n3*x2) + 1.0*n4*log(n4*x2) + 1.0*x0*log(x0*x2) + 1.0*x3*log(x2*x3)) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P)));
    return result;
}
        
static void coder_dgdn(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = pow(x1, -2);
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = (*endmember[2].mu0)(T, P);
    double x6 = (*endmember[3].mu0)(T, P);
    double x7 = 1.0/x1;
    double x8 = n2*x7;
    double x9 = 2.0*log(x8);
    double x10 = log(n3*x7);
    double x11 = 1.0*n3;
    double x12 = log(n4*x7);
    double x13 = 1.0*n4;
    double x14 = 1.0*log(x0*x7);
    double x15 = n1 + n4;
    double x16 = 1.0*log(x15*x7);
    double x17 = 8.3144626181532395*T;
    double x18 = x17*(n2*x9 + x0*x14 + x10*x11 + x12*x13 + x15*x16);
    double x19 = 1.0*x1;
    double x20 = -1.0*x2*(16050.0*n1*n4 + 24000.0*n3*n4 + 0.125*n4*(128400.0*n1 + 192000.0*n3) + x19*(n1*x3 + n2*x4 + n3*x5 + n4*x6 + x18));
    double x21 = 1.0*n2;
    double x22 = 1.0*n1;
    double x23 = x11 + x22;
    double x24 = x13 + x21 + x23;
    double x25 = -x11*x7;
    double x26 = -x13*x7;
    double x27 = x25 + x26;
    double x28 = -2.0*x8;
    double x29 = x14 + x28 + x1*x23*(-x0*x2 + x7)/x0;
    double x30 = x13 + x22;
    double x31 = x1*x30*(-x15*x2 + x7)/x15 + x16;
    double x32 = x11*x5 + x13*x6 + x18 + x21*x4 + x22*x3;
    double x33 = 1.0*x7;
    double x34 = -x23*x7;
    double x35 = -x30*x7;

result[0] = x20 + x33*(32100.0*n4 + x24*(x17*(x27 + x29 + x31) + x3) + x32);
result[1] = x20 + x33*(x24*(x17*(2.0*x1*(-n2*x2 + x7) + x27 + x34 + x35 + x9) + x4) + x32);
result[2] = x20 + x33*(48000.0*n4 + x24*(x17*(1.0*x10 + x19*(-n3*x2 + x7) + x26 + x29 + x35) + x5) + x32);
result[3] = x20 + x33*(32100.0*n1 + 48000.0*n3 + x24*(x17*(1.0*x12 + x19*(-n4*x2 + x7) + x25 + x28 + x31 + x34) + x6) + x32);
}
        
static void coder_d2gdn2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = (*endmember[1].mu0)(T, P);
    double x2 = (*endmember[2].mu0)(T, P);
    double x3 = (*endmember[3].mu0)(T, P);
    double x4 = n1 + n3;
    double x5 = n2 + n4 + x4;
    double x6 = 1.0/x5;
    double x7 = n2*x6;
    double x8 = 2.0*log(x7);
    double x9 = log(n3*x6);
    double x10 = 1.0*n3;
    double x11 = log(n4*x6);
    double x12 = 1.0*n4;
    double x13 = 1.0*log(x4*x6);
    double x14 = n1 + n4;
    double x15 = 1.0*log(x14*x6);
    double x16 = 8.3144626181532395*T;
    double x17 = x16*(n2*x8 + x10*x9 + x11*x12 + x13*x4 + x14*x15);
    double x18 = 1.0*x5;
    double x19 = pow(x5, -3);
    double x20 = 2.0*x19;
    double x21 = x20*(16050.0*n1*n4 + 24000.0*n3*n4 + 0.125*n4*(128400.0*n1 + 192000.0*n3) + x18*(n1*x0 + n2*x1 + n3*x2 + n4*x3 + x17));
    double x22 = 1.0*n2;
    double x23 = 1.0*n1;
    double x24 = x10 + x23;
    double x25 = x12 + x22 + x24;
    double x26 = -x10*x6;
    double x27 = -x12*x6;
    double x28 = x26 + x27;
    double x29 = -2.0*x7;
    double x30 = 1.0/x4;
    double x31 = pow(x5, -2);
    double x32 = -x31*x4 + x6;
    double x33 = x30*x32;
    double x34 = x24*x33;
    double x35 = x13 + x29 + x34*x5;
    double x36 = x12 + x23;
    double x37 = 1.0/x14;
    double x38 = -x14*x31 + x6;
    double x39 = x37*x38;
    double x40 = x36*x39;
    double x41 = x15 + x40*x5;
    double x42 = T*(x28 + x35 + x41);
    double x43 = 8.3144626181532395*x42;
    double x44 = x0*x23 + x1*x22 + x10*x2 + x12*x3 + x17;
    double x45 = 32100.0*n4 + x25*(x0 + x43) + x44;
    double x46 = 2.0*x31;
    double x47 = x10*x31;
    double x48 = x12*x31;
    double x49 = x47 + x48;
    double x50 = 2.0*x5;
    double x51 = -2*x31;
    double x52 = 2*x19;
    double x53 = x24*x5;
    double x54 = x30*x53;
    double x55 = n2*x31;
    double x56 = 2.0*x55;
    double x57 = x34 + x56;
    double x58 = -x32*x53/((x4)*(x4)) + x33*x50 + x54*(x4*x52 + x51) + x57;
    double x59 = x36*x5;
    double x60 = x37*x59;
    double x61 = x39*x50 + x40 + x60*(x14*x52 + x51) - x38*x59/((x14)*(x14));
    double x62 = x16*x25;
    double x63 = 1.0*x6;
    double x64 = -x31;
    double x65 = -n1;
    double x66 = x54*(-x52*(-n3 + x65) + x64) + x57;
    double x67 = x40 + x60*(-x52*(-n4 + x65) + x64);
    double x68 = 1.0*x0 + x43;
    double x69 = -x24*x6;
    double x70 = -x36*x6;
    double x71 = x50*(-x55 + x6);
    double x72 = x28 + x69 + x70 + x71 + x8;
    double x73 = x16*x72;
    double x74 = 1.0*x1 + x73;
    double x75 = x25*(x1 + x73) + x44;
    double x76 = 1.0*x31;
    double x77 = -x75*x76;
    double x78 = x21 - x45*x76;
    double x79 = 2.0*x6;
    double x80 = x49 - x79;
    double x81 = x18*(-n3*x31 + x6);
    double x82 = x27 + x35 + x70 + x81 + 1.0*x9;
    double x83 = x16*x82;
    double x84 = 1.0*x2 + x83;
    double x85 = 48000.0*n4 + x25*(x2 + x83) + x44;
    double x86 = -x76*x85;
    double x87 = x18*(-n4*x31 + x6);
    double x88 = 1.0*x11 + x26 + x29 + x41 + x69 + x87;
    double x89 = x16*x88;
    double x90 = 1.0*x3 + x89;
    double x91 = 32100.0*n1 + 48000.0*n3 + x25*(x3 + x89) + x44;
    double x92 = -x76*x91;
    double x93 = 16.628925236306479*T;
    double x94 = -x23;
    double x95 = -x31*(-x10 + x94);
    double x96 = 4.0*n2*x19;
    double x97 = -x31*(-x12 + x94);
    double x98 = -x56;
    double x99 = x97 + x98;
    double x100 = -x46;
    double x101 = x5*(x100 + x96) + x80;
    double x102 = x21 + x77;
    double x103 = x31*x36;
    double x104 = n3*x20;
    double x105 = -x47 + x48;

result[0] = x21 - x45*x46 + x63*(2.0*x0 + 16.628925236306479*x42 + x62*(x49 + x58 + x61));
result[1] = x63*(x62*(x49 - 4.0*x6 + x66 + x67) + x68 + x74) + x77 + x78;
result[2] = x63*(x62*(x58 + x67 + x80) + x68 + x84) + x78 + x86;
result[3] = x63*(x62*(x61 + x66 + x80) + x68 + x90 + 32100.0) + x78 + x92;
result[4] = x21 - x46*x75 + x63*(2.0*x1 + x62*(x49 + x5*(-4.0*x31 + x96) + x79 + x95 + x99 + x71/n2) + x72*x93);
result[5] = x102 + x63*(x62*(x101 + x24*x31 + x99) + x74 + x84) + x86;
result[6] = x102 + x63*(x62*(x101 + x103 + x95 + x98) + x74 + x90) + x92;
result[7] = x21 - x46*x85 + x63*(2.0*x2 + x62*(x105 + x5*(x100 + x104) + x58 + x63 + x97 + x81/n3) + x82*x93);
result[8] = x21 + x63*(x62*(x103 + x105 + x5*(x104 - x76) - 3.0*x6 + x66) + x84 + x90 + 48000.0) + x86 + x92;
result[9] = x21 - x46*x91 + x63*(2.0*x3 + x62*(x47 - x48 + x5*(n4*x20 + x100) + x56 + x61 + x63 + x95 + x87/n4) + x88*x93);
}
        
static void coder_d3gdn3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].mu0)(T, P);
    double x1 = (*endmember[1].mu0)(T, P);
    double x2 = (*endmember[2].mu0)(T, P);
    double x3 = (*endmember[3].mu0)(T, P);
    double x4 = n1 + n3;
    double x5 = n2 + n4 + x4;
    double x6 = 1.0/x5;
    double x7 = n2*x6;
    double x8 = 2.0*log(x7);
    double x9 = log(n3*x6);
    double x10 = 1.0*n3;
    double x11 = log(n4*x6);
    double x12 = 1.0*n4;
    double x13 = 1.0*log(x4*x6);
    double x14 = n1 + n4;
    double x15 = 1.0*log(x14*x6);
    double x16 = 8.3144626181532395*T;
    double x17 = x16*(n2*x8 + x10*x9 + x11*x12 + x13*x4 + x14*x15);
    double x18 = 1.0*x5;
    double x19 = pow(x5, -4);
    double x20 = 6.0*x19;
    double x21 = -x20*(16050.0*n1*n4 + 24000.0*n3*n4 + 0.125*n4*(128400.0*n1 + 192000.0*n3) + x18*(n1*x0 + n2*x1 + n3*x2 + n4*x3 + x17));
    double x22 = 1.0*n2;
    double x23 = 1.0*n1;
    double x24 = x10 + x23;
    double x25 = x12 + x22 + x24;
    double x26 = -x10*x6;
    double x27 = -x12*x6;
    double x28 = x26 + x27;
    double x29 = -2.0*x7;
    double x30 = 1.0/x4;
    double x31 = pow(x5, -2);
    double x32 = -x31*x4 + x6;
    double x33 = x30*x32;
    double x34 = x24*x33;
    double x35 = x13 + x29 + x34*x5;
    double x36 = x12 + x23;
    double x37 = 1.0/x14;
    double x38 = -x14*x31 + x6;
    double x39 = x37*x38;
    double x40 = x36*x39;
    double x41 = x15 + x40*x5;
    double x42 = x28 + x35 + x41;
    double x43 = x16*x42;
    double x44 = x0*x23 + x1*x22 + x10*x2 + x12*x3 + x17;
    double x45 = 32100.0*n4 + x25*(x0 + x43) + x44;
    double x46 = pow(x5, -3);
    double x47 = 6.0*x46;
    double x48 = 16.628925236306479*T;
    double x49 = x10*x31;
    double x50 = x12*x31;
    double x51 = x49 + x50;
    double x52 = 2.0*x33;
    double x53 = -2*x31;
    double x54 = 2*x46;
    double x55 = x4*x54 + x53;
    double x56 = x24*x30;
    double x57 = x55*x56;
    double x58 = pow(x4, -2);
    double x59 = x32*x58;
    double x60 = x24*x59;
    double x61 = 2.0*x31;
    double x62 = n2*x61;
    double x63 = x34 + x62;
    double x64 = x5*x52 + x5*x57 - x5*x60 + x63;
    double x65 = 2.0*x39;
    double x66 = x14*x54 + x53;
    double x67 = x36*x37;
    double x68 = x66*x67;
    double x69 = pow(x14, -2);
    double x70 = x38*x69;
    double x71 = x36*x70;
    double x72 = x40 + x5*x65 + x5*x68 - x5*x71;
    double x73 = T*(x51 + x64 + x72);
    double x74 = 8.3144626181532395*x73;
    double x75 = 2.0*x0 + x25*x74 + x42*x48;
    double x76 = 3.0*x31;
    double x77 = 3.0*x5;
    double x78 = 6*x46;
    double x79 = 6*x19;
    double x80 = x5*x67;
    double x81 = x36*x5;
    double x82 = x69*x81;
    double x83 = x37*x66*x77 + 3.0*x39 - 2*x66*x82 + 2*x68 - x70*x77 - 2*x71 + x80*(-x14*x79 + x78) + 2*x38*x81/((x14)*(x14)*(x14));
    double x84 = 4.0*x46;
    double x85 = n2*x84;
    double x86 = -x85;
    double x87 = 2.0*x46;
    double x88 = n3*x87;
    double x89 = -x88;
    double x90 = n4*x87;
    double x91 = -x90;
    double x92 = x89 + x91;
    double x93 = x86 + x92;
    double x94 = x5*x56;
    double x95 = x24*x5;
    double x96 = x58*x95;
    double x97 = x30*x55*x77 + 2*x32*x95/((x4)*(x4)*(x4)) + 3.0*x33 - 2*x55*x96 + 2*x57 - x59*x77 - 2*x60 + x94*(-x4*x79 + x78);
    double x98 = x93 + x97;
    double x99 = x16*x25;
    double x100 = 1.0*x6;
    double x101 = -x31;
    double x102 = -n1;
    double x103 = -n3 + x102;
    double x104 = x101 - x103*x54;
    double x105 = x104*x56;
    double x106 = x105*x5 + x63;
    double x107 = -n4 + x102;
    double x108 = x101 - x107*x54;
    double x109 = x108*x67;
    double x110 = x109*x5 + x40;
    double x111 = x106 + x110 + x51 - 4.0*x6;
    double x112 = x111*x48;
    double x113 = x61 + x92;
    double x114 = 2.0*x5;
    double x115 = x104*x30;
    double x116 = 4*x46;
    double x117 = 2*n1;
    double x118 = 2*n3;
    double x119 = 3*x19;
    double x120 = -x119*(x117 + x118);
    double x121 = -x104*x96 + x105 + x57 - x60;
    double x122 = x114*x115 + x121 + x52 + x94*(x116 + x120);
    double x123 = x108*x37;
    double x124 = 2*n4;
    double x125 = -x119*(x117 + x124);
    double x126 = -x108*x82 + x109 + x68 - x71;
    double x127 = x114*x123 + x126 + x65 + x80*(x116 + x125);
    double x128 = -x24*x6;
    double x129 = -x36*x6;
    double x130 = -2.0*n2*x31 + 2.0*x6;
    double x131 = x130*x5;
    double x132 = x128 + x129 + x131 + x28 + x8;
    double x133 = x132*x16;
    double x134 = x25*(x1 + x133) + x44;
    double x135 = x134*x87;
    double x136 = x111*x16;
    double x137 = 1.0*x0 + x43;
    double x138 = 1.0*x1 + x133;
    double x139 = x136*x25 + x137 + x138;
    double x140 = -x139*x61 + x21;
    double x141 = 1.0*x31;
    double x142 = -x141*x75 + x45*x84;
    double x143 = 2.0*x6;
    double x144 = -x143 + x51;
    double x145 = x110 + x144 + x64;
    double x146 = x145*x48;
    double x147 = -n3*x31 + x6;
    double x148 = x147*x18;
    double x149 = x129 + x148 + x27 + x35 + 1.0*x9;
    double x150 = x149*x16;
    double x151 = 48000.0*n4 + x25*(x150 + x2) + x44;
    double x152 = x151*x87;
    double x153 = x145*x16;
    double x154 = x150 + 1.0*x2;
    double x155 = x137 + x153*x25 + x154;
    double x156 = -x155*x61;
    double x157 = x142 + x21;
    double x158 = x106 + x144 + x72;
    double x159 = x158*x48;
    double x160 = x122 + x93;
    double x161 = -n4*x31 + x6;
    double x162 = x161*x18;
    double x163 = 1.0*x11 + x128 + x162 + x26 + x29 + x41;
    double x164 = x16*x163;
    double x165 = 32100.0*n1 + 48000.0*n3 + x25*(x164 + x3) + x44;
    double x166 = x165*x87;
    double x167 = x158*x16;
    double x168 = x164 + 1.0*x3;
    double x169 = x137 + x167*x25 + x168 + 32100.0;
    double x170 = -x169*x61;
    double x171 = -x23;
    double x172 = -x10 + x171;
    double x173 = -x172*x31;
    double x174 = 4.0*x31;
    double x175 = -x174;
    double x176 = 1.0/n2;
    double x177 = -x12 + x171;
    double x178 = -x177*x31;
    double x179 = -x62;
    double x180 = x178 + x179;
    double x181 = x131*x176 + x143 + x173 + x180 + x5*(x175 + x85) + x51;
    double x182 = x16*x181;
    double x183 = 2*x105 + x94*(x120 + x54);
    double x184 = 2*x109 + x80*(x125 + x54) + x93;
    double x185 = x45*x87;
    double x186 = 2.0*x1 + x132*x48 + x182*x25;
    double x187 = x134*x84 - x141*x186;
    double x188 = -x61;
    double x189 = x144 + x5*(x188 + x85);
    double x190 = x180 + x189 + x24*x31;
    double x191 = x16*x190;
    double x192 = 5.0*x31;
    double x193 = -x141*x155;
    double x194 = x138 + x154 + x191*x25;
    double x195 = -x141*x194 + x152;
    double x196 = x185 + x21;
    double x197 = x135 - x139*x141 + x196;
    double x198 = x31*x36;
    double x199 = x173 + x179 + x189 + x198;
    double x200 = x16*x199;
    double x201 = x183 + x93;
    double x202 = x123*x18 + x126 + 1.0*x39 + x80*(x107*x79 + x116);
    double x203 = x138 + x168 + x200*x25;
    double x204 = -x141*x203;
    double x205 = -x141*x169 + x166;
    double x206 = 1.0/n3;
    double x207 = -x49 + x50;
    double x208 = x100 + x148*x206 + x178 + x207 + x5*(x188 + x88) + x64;
    double x209 = x16*x208;
    double x210 = x149*x48 + 2.0*x2 + x209*x25;
    double x211 = -x141*x210 + x151*x84;
    double x212 = -x141;
    double x213 = x106 + x198 + x207 + x5*(x212 + x88) - 3.0*x6;
    double x214 = x16*x213;
    double x215 = x154 + x168 + x214*x25 + 48000.0;
    double x216 = -x141*x215;
    double x217 = 1.0/n4;
    double x218 = x100 + x162*x217 + x173 + x49 + x5*(x188 + x90) - x50 + x62 + x72;
    double x219 = x16*x218;
    double x220 = x163*x48 + x219*x25 + 2.0*x3;
    double x221 = -x141*x220 + x165*x84;
    double x222 = 24.943387854459719*T;
    double x223 = -12.0*n2*x19;
    double x224 = n2*x54;
    double x225 = x114*x176;
    double x226 = -x24*x54;
    double x227 = -x36*x54;
    double x228 = 8.0*x46;
    double x229 = n2*x228;
    double x230 = x227 + x229;
    double x231 = x226 + x230;
    double x232 = x130*x176 + x92;
    double x233 = x190*x48;
    double x234 = x175 + x225*(x101 + x224) + x232 + x5*(x223 + x228);
    double x235 = -x194*x61;
    double x236 = x187 + x21;
    double x237 = x199*x48;
    double x238 = x177*x54;
    double x239 = x226 + x229 + x238;
    double x240 = -x203*x61;
    double x241 = x113 + x5*(x223 + x84);
    double x242 = x99*(x231 + x241);
    double x243 = x135 + x21;
    double x244 = 1.0*x147*x206;
    double x245 = -n3*x20;
    double x246 = x118*x46;
    double x247 = x18*x206;
    double x248 = n3*x84 + x86 + x91;
    double x249 = x227 + x248;
    double x250 = x213*x48;
    double x251 = x21 - x215*x61;

result[0] = x100*(24.943387854459719*x73 + x99*(x83 + x98)) + x21 + x45*x47 - x75*x76;
result[1] = x100*(x112 + x74 + x99*(x113 + x122 + x127 + x86)) + x135 + x140 + x142;
result[2] = x100*(x146 + x74 + x99*(x127 + x141 + x98)) + x152 + x156 + x157;
result[3] = x100*(x159 + x74 + x99*(x141 + x160 + x83)) + x157 + x166 + x170;
result[4] = x100*(x112 + x182 + x99*(x183 + x184 + 6.0*x31)) + x140 + x185 + x187;
result[5] = x100*(x136 + x153 + x191 + x99*(x115*x18 + x121 + x184 + x192 + 1.0*x33 + x94*(x103*x79 + x116))) + x193 + x195 + x197;
result[6] = x100*(x136 + x167 + x200 + x99*(x192 + x201 + x202)) + x197 + x204 + x205;
result[7] = x100*(x146 + x209 + x99*(x184 + x76 + x97)) + x156 + x196 + x211;
result[8] = x100*(x153 + x167 + x214 + x99*(x160 + x202 + x76)) + x152 + x193 + x196 + x205 + x216;
result[9] = x100*(x159 + x219 + x99*(x201 + x76 + x83)) + x170 + x196 + x221;
result[10] = x100*(x181*x222 + x99*(x225*(x224 + x53) + x231 + x232 - 8.0*x31 + x5*(x223 + 12.0*x46) - x131/((n2)*(n2)))) + x134*x47 - x186*x76 + x21;
result[11] = x100*(x182 + x233 + x99*(x172*x54 + x230 + x234)) + x152 + x235 + x236;
result[12] = x100*(x182 + x237 + x99*(x234 + x239)) + x166 + x236 + x240;
result[13] = x100*(x209 + x233 + x242) + x211 + x235 + x243;
result[14] = x100*(x191 + x200 + x214 + x99*(x239 + x241)) + x166 + x195 + x204 + x216 + x243;
result[15] = x100*(x219 + x237 + x242) + x221 + x240 + x243;
result[16] = x100*(x208*x222 + x99*(x175 + x244 + x247*(x246 + x53) + x249 + x5*(x245 + x47) + x97 - x148/((n3)*(n3)))) + x151*x47 + x21 - x210*x76;
result[17] = x100*(x209 + x250 + x99*(x122 + x212 + x238 + x244 + x247*(x101 + x246) + x248 + x5*(x245 + x84))) + x166 + x211 + x251;
result[18] = x100*(x219 + x250 + x99*(x174 + x183 + x249 + x5*(x245 + x87))) + x152 + x221 + x251;
result[19] = x100*(x218*x222 + x99*(n4*x84 + 1.0*x161*x217 + x175 + x18*x217*(x124*x46 + x53) + x226 + x5*(-n4*x20 + x47) + x83 + x86 + x89 - x162/((n4)*(n4)))) + x165*x47 + x21 - x220*x76;
}
        
static double coder_dgdt(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    
    double x0 = n1 + n3;
    double x1 = 1.0/(n2 + n4 + x0);
    double x2 = n1 + n4;

result = 1.0*x1*(1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4)*(n1*(*endmember[0].dmu0dT)(T, P) + 16.628925236306479*n2*log(n2*x1) + n2*(*endmember[1].dmu0dT)(T, P) + 8.3144626181532395*n3*log(n3*x1) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*log(n4*x1) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626181532395*x0*log(x0*x1) + 8.3144626181532395*x2*log(x1*x2));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].dmu0dT)(T, P);
    double x1 = n1 + n3;
    double x2 = n2 + n4 + x1;
    double x3 = 1.0/x2;
    double x4 = n3*x3;
    double x5 = -8.3144626181532395*x4;
    double x6 = n4*x3;
    double x7 = -8.3144626181532395*x6;
    double x8 = x5 + x7;
    double x9 = 8.3144626181532395*log(x1*x3);
    double x10 = n2*x3;
    double x11 = -16.628925236306479*x10;
    double x12 = 8.3144626181532395*n1;
    double x13 = 8.3144626181532395*n3;
    double x14 = x12 + x13;
    double x15 = pow(x2, -2);
    double x16 = x11 + x9 + x14*x2*(-x1*x15 + x3)/x1;
    double x17 = n1 + n4;
    double x18 = 8.3144626181532395*log(x17*x3);
    double x19 = 8.3144626181532395*n4;
    double x20 = x12 + x19;
    double x21 = x18 + x2*x20*(-x15*x17 + x3)/x17;
    double x22 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x23 = 1.0*x3;
    double x24 = x22*x23;
    double x25 = (*endmember[1].dmu0dT)(T, P);
    double x26 = (*endmember[2].dmu0dT)(T, P);
    double x27 = (*endmember[3].dmu0dT)(T, P);
    double x28 = 16.628925236306479*log(x10);
    double x29 = log(x4);
    double x30 = log(x6);
    double x31 = n1*x0 + n2*x25 + n2*x28 + n3*x26 + n4*x27 + x1*x9 + x13*x29 + x17*x18 + x19*x30;
    double x32 = -1.0*x15*x22*x31 + x23*x31;
    double x33 = -x14*x3;
    double x34 = -x20*x3;
    double x35 = 8.3144626181532395*x2;

result[0] = x24*(x0 + x16 + x21 + x8) + x32;
result[1] = x24*(16.628925236306479*x2*(-n2*x15 + x3) + x25 + x28 + x33 + x34 + x8) + x32;
result[2] = x24*(x16 + x26 + 8.3144626181532395*x29 + x34 + x35*(-n3*x15 + x3) + x7) + x32;
result[3] = x24*(x11 + x21 + x27 + 8.3144626181532395*x30 + x33 + x35*(-n4*x15 + x3) + x5) + x32;
}
        
static void coder_d3gdn2dt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = 1.0/x1;
    double x3 = (*endmember[0].dmu0dT)(T, P);
    double x4 = 8.3144626181532395*n3;
    double x5 = -x2*x4;
    double x6 = 8.3144626181532395*n4;
    double x7 = -x2*x6;
    double x8 = x5 + x7;
    double x9 = 8.3144626181532395*log(x0*x2);
    double x10 = n2*x2;
    double x11 = -16.628925236306479*x10;
    double x12 = 8.3144626181532395*n1;
    double x13 = x12 + x4;
    double x14 = 1.0/x0;
    double x15 = pow(x1, -2);
    double x16 = -x0*x15 + x2;
    double x17 = x14*x16;
    double x18 = x13*x17;
    double x19 = x1*x18 + x11 + x9;
    double x20 = n1 + n4;
    double x21 = 8.3144626181532395*log(x2*x20);
    double x22 = x12 + x6;
    double x23 = 1.0/x20;
    double x24 = -x15*x20 + x2;
    double x25 = x23*x24;
    double x26 = x22*x25;
    double x27 = x1*x26 + x21;
    double x28 = x19 + x27 + x3 + x8;
    double x29 = x2*x28;
    double x30 = x15*x4;
    double x31 = x15*x6;
    double x32 = x30 + x31;
    double x33 = 16.628925236306479*x1;
    double x34 = -2*x15;
    double x35 = pow(x1, -3);
    double x36 = 2*x35;
    double x37 = x1*x13;
    double x38 = x14*x37;
    double x39 = n2*x15;
    double x40 = 16.628925236306479*x39;
    double x41 = x18 + x40;
    double x42 = x17*x33 + x38*(x0*x36 + x34) + x41 - x16*x37/((x0)*(x0));
    double x43 = x1*x22;
    double x44 = x23*x43;
    double x45 = x25*x33 + x26 + x44*(x20*x36 + x34) - x24*x43/((x20)*(x20));
    double x46 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x47 = 1.0*x2;
    double x48 = x46*x47;
    double x49 = x15*x46;
    double x50 = x28*x49;
    double x51 = (*endmember[1].dmu0dT)(T, P);
    double x52 = (*endmember[2].dmu0dT)(T, P);
    double x53 = (*endmember[3].dmu0dT)(T, P);
    double x54 = 16.628925236306479*log(x10);
    double x55 = log(n3*x2);
    double x56 = log(n4*x2);
    double x57 = 2.0*n1*x3 + 2.0*n2*x51 + 2.0*n2*x54 + 2.0*n3*x52 + 2.0*n4*x53 + 2.0*x0*x9 + 2.0*x20*x21 + 2.0*x4*x55 + 2.0*x56*x6;
    double x58 = -x15*x57 + x35*x46*x57;
    double x59 = -x15;
    double x60 = -n1;
    double x61 = x38*(-x36*(-n3 + x60) + x59) + x41;
    double x62 = x26 + x44*(-x36*(-n4 + x60) + x59);
    double x63 = 1.0*x29 - 1.0*x50 + x58;
    double x64 = -x13*x2;
    double x65 = -x2*x22;
    double x66 = x33*(x2 - x39);
    double x67 = x51 + x54 + x64 + x65 + x66 + x8;
    double x68 = 1.0*x49;
    double x69 = x47*x67 - x67*x68;
    double x70 = 16.628925236306479*x2;
    double x71 = x32 - x70;
    double x72 = 8.3144626181532395*x1;
    double x73 = x72*(-n3*x15 + x2);
    double x74 = x19 + x52 + 8.3144626181532395*x55 + x65 + x7 + x73;
    double x75 = x47*x74 - x68*x74;
    double x76 = x72*(-n4*x15 + x2);
    double x77 = x11 + x27 + x5 + x53 + 8.3144626181532395*x56 + x64 + x76;
    double x78 = x47*x77 - x68*x77;
    double x79 = 2.0*x2;
    double x80 = -x12;
    double x81 = -x15*(-x4 + x80);
    double x82 = 33.257850472612958*n2*x35;
    double x83 = -x15*(-x6 + x80);
    double x84 = -x40;
    double x85 = x83 + x84;
    double x86 = 2.0*x49;
    double x87 = -16.628925236306479*x15;
    double x88 = x1*(x82 + x87) + x71;
    double x89 = x58 + x69;
    double x90 = x15*x22;
    double x91 = 8.3144626181532395*x2;
    double x92 = 16.628925236306479*x35;
    double x93 = n3*x92;
    double x94 = -x30 + x31;

result[0] = 2.0*x29 + x48*(x32 + x42 + x45) - 2.0*x50 + x58;
result[1] = x48*(-33.257850472612958*x2 + x32 + x61 + x62) + x63 + x69;
result[2] = x48*(x42 + x62 + x71) + x63 + x75;
result[3] = x48*(x45 + x61 + x71) + x63 + x78;
result[4] = x48*(x1*(-33.257850472612958*x15 + x82) + x32 + x70 + x81 + x85 + x66/n2) + x58 + x67*x79 - x67*x86;
result[5] = x48*(x13*x15 + x85 + x88) + x75 + x89;
result[6] = x48*(x81 + x84 + x88 + x90) + x78 + x89;
result[7] = x48*(x1*(x87 + x93) + x42 + x83 + x91 + x94 + x73/n3) + x58 + x74*x79 - x74*x86;
result[8] = x48*(x1*(-8.3144626181532395*x15 + x93) - 24.943387854459719*x2 + x61 + x90 + x94) + x58 + x75 + x78;
result[9] = x48*(x1*(n4*x92 + x87) + x30 - x31 + x40 + x45 + x81 + x91 + x76/n4) + x58 + x77*x79 - x77*x86;
}
        
static void coder_d4gdn3dt(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = n1 + n3;
    double x1 = n2 + n4 + x0;
    double x2 = 1.0/x1;
    double x3 = pow(x1, -2);
    double x4 = 8.3144626181532395*n3;
    double x5 = x3*x4;
    double x6 = 8.3144626181532395*n4;
    double x7 = x3*x6;
    double x8 = x5 + x7;
    double x9 = 1.0/x0;
    double x10 = -x0*x3 + x2;
    double x11 = x10*x9;
    double x12 = 16.628925236306479*x11;
    double x13 = -2*x3;
    double x14 = pow(x1, -3);
    double x15 = 2*x14;
    double x16 = x0*x15 + x13;
    double x17 = 8.3144626181532395*n1;
    double x18 = x17 + x4;
    double x19 = x18*x9;
    double x20 = x16*x19;
    double x21 = pow(x0, -2);
    double x22 = x10*x21;
    double x23 = x18*x22;
    double x24 = n2*x3;
    double x25 = 16.628925236306479*x24;
    double x26 = x11*x18;
    double x27 = x25 + x26;
    double x28 = x1*x12 + x1*x20 - x1*x23 + x27;
    double x29 = x17 + x6;
    double x30 = n1 + n4;
    double x31 = 1.0/x30;
    double x32 = x2 - x3*x30;
    double x33 = x31*x32;
    double x34 = x29*x33;
    double x35 = 16.628925236306479*x33;
    double x36 = x13 + x15*x30;
    double x37 = x29*x31;
    double x38 = x36*x37;
    double x39 = pow(x30, -2);
    double x40 = x32*x39;
    double x41 = x29*x40;
    double x42 = x1*x35 + x1*x38 - x1*x41 + x34;
    double x43 = x28 + x42 + x8;
    double x44 = x2*x43;
    double x45 = (*endmember[0].dmu0dT)(T, P);
    double x46 = -x2*x4;
    double x47 = -x2*x6;
    double x48 = x46 + x47;
    double x49 = 8.3144626181532395*log(x0*x2);
    double x50 = n2*x2;
    double x51 = -16.628925236306479*x50;
    double x52 = x1*x26 + x49 + x51;
    double x53 = 8.3144626181532395*log(x2*x30);
    double x54 = x1*x34 + x53;
    double x55 = x45 + x48 + x52 + x54;
    double x56 = x3*x55;
    double x57 = 24.943387854459719*x1;
    double x58 = 6*x14;
    double x59 = pow(x1, -4);
    double x60 = 6*x59;
    double x61 = x1*x37;
    double x62 = x1*x29;
    double x63 = x39*x62;
    double x64 = x31*x36*x57 + 24.943387854459719*x33 - 2*x36*x63 + 2*x38 - x40*x57 - 2*x41 + x61*(-x30*x60 + x58) + 2*x32*x62/((x30)*(x30)*(x30));
    double x65 = 33.257850472612958*x14;
    double x66 = n2*x65;
    double x67 = -x66;
    double x68 = 16.628925236306479*x14;
    double x69 = n3*x68;
    double x70 = -x69;
    double x71 = n4*x68;
    double x72 = -x71;
    double x73 = x70 + x72;
    double x74 = x67 + x73;
    double x75 = x1*x19;
    double x76 = x1*x18;
    double x77 = x21*x76;
    double x78 = 24.943387854459719*x11 + x16*x57*x9 - 2*x16*x77 + 2*x20 - x22*x57 - 2*x23 + x75*(-x0*x60 + x58) + 2*x10*x76/((x0)*(x0)*(x0));
    double x79 = x74 + x78;
    double x80 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x81 = 1.0*x2;
    double x82 = x80*x81;
    double x83 = 6.0*x14;
    double x84 = x55*x80;
    double x85 = x3*x80;
    double x86 = x43*x85;
    double x87 = (*endmember[1].dmu0dT)(T, P);
    double x88 = (*endmember[2].dmu0dT)(T, P);
    double x89 = (*endmember[3].dmu0dT)(T, P);
    double x90 = 16.628925236306479*log(x50);
    double x91 = log(n3*x2);
    double x92 = log(n4*x2);
    double x93 = n1*x45 + n2*x87 + n2*x90 + n3*x88 + n4*x89 + x0*x49 + x30*x53 + x4*x91 + x6*x92;
    double x94 = -6.0*x59*x80*x93 + x83*x93;
    double x95 = 16.628925236306479*x3;
    double x96 = x73 + x95;
    double x97 = 16.628925236306479*x1;
    double x98 = -x3;
    double x99 = -n1;
    double x100 = -n3 + x99;
    double x101 = -x100*x15 + x98;
    double x102 = x101*x9;
    double x103 = 4*x14;
    double x104 = 2*n1;
    double x105 = 2*n3;
    double x106 = 3*x59;
    double x107 = -x106*(x104 + x105);
    double x108 = x101*x19;
    double x109 = -x101*x77 + x108 + x20 - x23;
    double x110 = x102*x97 + x109 + x12 + x75*(x103 + x107);
    double x111 = -n4 + x99;
    double x112 = -x111*x15 + x98;
    double x113 = x112*x31;
    double x114 = 2*n4;
    double x115 = -x106*(x104 + x114);
    double x116 = x112*x37;
    double x117 = -x112*x63 + x116 + x38 - x41;
    double x118 = x113*x97 + x117 + x35 + x61*(x103 + x115);
    double x119 = x1*x108 + x27;
    double x120 = x1*x116 + x34;
    double x121 = x119 + x120 - 33.257850472612958*x2 + x8;
    double x122 = 2.0*x2;
    double x123 = 2.0*x3;
    double x124 = x123*x80;
    double x125 = x121*x122 - x121*x124 + x94;
    double x126 = -x18*x2;
    double x127 = -x2*x29;
    double x128 = 16.628925236306479*x2 - 16.628925236306479*x24;
    double x129 = x1*x128;
    double x130 = x126 + x127 + x129 + x48 + x87 + x90;
    double x131 = 2.0*x14;
    double x132 = x130*x80;
    double x133 = -x123*x130 + x131*x132;
    double x134 = 4.0*x14;
    double x135 = x134*x84 + 1.0*x44 - 4.0*x56 - 1.0*x86;
    double x136 = 8.3144626181532395*x3;
    double x137 = x135 + x94;
    double x138 = 16.628925236306479*x2;
    double x139 = -x138 + x8;
    double x140 = x120 + x139 + x28;
    double x141 = x122*x140 - x124*x140;
    double x142 = -n3*x3 + x2;
    double x143 = 8.3144626181532395*x1;
    double x144 = x142*x143;
    double x145 = x127 + x144 + x47 + x52 + x88 + 8.3144626181532395*x91;
    double x146 = x131*x80;
    double x147 = -x123*x145 + x145*x146;
    double x148 = x110 + x74;
    double x149 = x119 + x139 + x42;
    double x150 = x122*x149 - x124*x149;
    double x151 = -n4*x3 + x2;
    double x152 = x143*x151;
    double x153 = x126 + x152 + x46 + x51 + x54 + x89 + 8.3144626181532395*x92;
    double x154 = -x123*x153 + x146*x153;
    double x155 = 2*x108 + x75*(x107 + x15);
    double x156 = 2*x116 + x61*(x115 + x15) + x74;
    double x157 = x131*x84 - 2.0*x56;
    double x158 = -x17;
    double x159 = x158 - x4;
    double x160 = -x159*x3;
    double x161 = 33.257850472612958*x3;
    double x162 = -x161;
    double x163 = 1.0/n2;
    double x164 = x158 - x6;
    double x165 = -x164*x3;
    double x166 = -x25;
    double x167 = x165 + x166;
    double x168 = x1*(x162 + x66) + x129*x163 + x138 + x160 + x167 + x8;
    double x169 = 4.0*x3;
    double x170 = 1.0*x85;
    double x171 = -x130*x169 + x132*x134 - x168*x170 + x168*x81;
    double x172 = 41.572313090766201*x3;
    double x173 = x157 + x94;
    double x174 = -x121*x170 + x121*x81 + x133 + x173;
    double x175 = -x140*x170 + x140*x81 + x147;
    double x176 = -x95;
    double x177 = x1*(x176 + x66) + x139;
    double x178 = x167 + x177 + x18*x3;
    double x179 = -x170*x178 + x178*x81;
    double x180 = x155 + x74;
    double x181 = x113*x143 + x117 + 8.3144626181532395*x33 + x61*(x103 + x111*x60);
    double x182 = -x149*x170 + x149*x81 + x154;
    double x183 = x29*x3;
    double x184 = x160 + x166 + x177 + x183;
    double x185 = -x170*x184 + x184*x81;
    double x186 = 24.943387854459719*x3;
    double x187 = 8.3144626181532395*x2;
    double x188 = 1.0/n3;
    double x189 = -x5 + x7;
    double x190 = x1*(x176 + x69) + x144*x188 + x165 + x187 + x189 + x28;
    double x191 = x134*x80;
    double x192 = -x145*x169 + x145*x191 - x170*x190 + x190*x81;
    double x193 = -x136;
    double x194 = x1*(x193 + x69) + x119 + x183 + x189 - 24.943387854459719*x2;
    double x195 = -x170*x194 + x194*x81;
    double x196 = 1.0/n4;
    double x197 = x1*(x176 + x71) + x152*x196 + x160 + x187 + x25 + x42 + x5 - x7;
    double x198 = -x153*x169 + x153*x191 - x170*x197 + x197*x81;
    double x199 = 3.0*x2;
    double x200 = 6.0*x3;
    double x201 = -99.773551417838874*n2*x59;
    double x202 = n2*x15;
    double x203 = x163*x97;
    double x204 = -x15*x18;
    double x205 = -x15*x29;
    double x206 = 66.515700945225916*x14;
    double x207 = n2*x206;
    double x208 = x205 + x207;
    double x209 = x204 + x208;
    double x210 = x128*x163 + x73;
    double x211 = 3.0*x85;
    double x212 = x1*(x201 + x206) + x162 + x203*(x202 + x98) + x210;
    double x213 = x171 + x94;
    double x214 = x122*x178 - x124*x178;
    double x215 = x15*x164;
    double x216 = x204 + x207 + x215;
    double x217 = x122*x184 - x124*x184;
    double x218 = x192 + x94;
    double x219 = x1*(x201 + x65) + x96;
    double x220 = x133 + x82*(x209 + x219);
    double x221 = x198 + x94;
    double x222 = 8.3144626181532395*x142*x188;
    double x223 = 49.886775708919437*x14;
    double x224 = 49.886775708919437*x59;
    double x225 = -n3*x224;
    double x226 = x105*x14;
    double x227 = x143*x188;
    double x228 = n3*x65 + x67 + x72;
    double x229 = x205 + x228;
    double x230 = x80*x83;
    double x231 = x122*x194 - x124*x194;

result[0] = 3.0*x44 - 6.0*x56 + x82*(x64 + x79) + x83*x84 - 3.0*x86 + x94;
result[1] = x125 + x133 + x135 + x82*(x110 + x118 + x67 + x96);
result[2] = x137 + x141 + x147 + x82*(x118 + x136 + x79);
result[3] = x137 + x150 + x154 + x82*(x136 + x148 + x64);
result[4] = x125 + x157 + x171 + x82*(x155 + x156 + 49.886775708919437*x3);
result[5] = x174 + x175 + x179 + x82*(x102*x143 + x109 + 8.3144626181532395*x11 + x156 + x172 + x75*(x100*x60 + x103));
result[6] = x174 + x182 + x185 + x82*(x172 + x180 + x181);
result[7] = x141 + x173 + x192 + x82*(x156 + x186 + x78);
result[8] = x173 + x175 + x182 + x195 + x82*(x148 + x181 + x186);
result[9] = x150 + x173 + x198 + x82*(x180 + x186 + x64);
result[10] = -x130*x200 + x132*x83 + x168*x199 - x168*x211 + x82*(x1*(99.773551417838874*x14 + x201) + x203*(x13 + x202) + x209 + x210 - 66.515700945225916*x3 - x129/((n2)*(n2))) + x94;
result[11] = x147 + x213 + x214 + x82*(x15*x159 + x208 + x212);
result[12] = x154 + x213 + x217 + x82*(x212 + x216);
result[13] = x214 + x218 + x220;
result[14] = x133 + x147 + x154 + x179 + x185 + x195 + x82*(x216 + x219) + x94;
result[15] = x217 + x220 + x221;
result[16] = -x145*x200 + x145*x230 + x190*x199 - x190*x211 + x82*(x1*(x223 + x225) + x162 + x222 + x227*(x13 + x226) + x229 + x78 - x144/((n3)*(n3))) + x94;
result[17] = x154 + x218 + x231 + x82*(x1*(x225 + x65) + x110 + x193 + x215 + x222 + x227*(x226 + x98) + x228);
result[18] = x147 + x221 + x231 + x82*(x1*(x225 + x68) + x155 + x161 + x229);
result[19] = -x153*x200 + x153*x230 + x197*x199 - x197*x211 + x82*(n4*x65 + x1*(-n4*x224 + x223) + x143*x196*(x114*x14 + x13) + 8.3144626181532395*x151*x196 + x162 + x204 + x64 + x67 + x70 - x152/((n4)*(n4))) + x94;
}
        
static double coder_dgdp(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = 1.0*(1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4)*(n1*(*endmember[0].dmu0dP)(T, P) + n2*(*endmember[1].dmu0dP)(T, P) + n3*(*endmember[2].dmu0dP)(T, P) + n4*(*endmember[3].dmu0dP)(T, P))/(n1 + n2 + n3 + n4);
    return result;
}
        
static void coder_d2gdndp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x2 = n1 + n2 + n3 + n4;
    double x3 = 1.0/x2;
    double x4 = x1*x3;
    double x5 = (*endmember[1].dmu0dP)(T, P);
    double x6 = (*endmember[2].dmu0dP)(T, P);
    double x7 = (*endmember[3].dmu0dP)(T, P);
    double x8 = n1*x0 + n2*x5 + n3*x6 + n4*x7;
    double x9 = -1.0*x1*x8/((x2)*(x2)) + x3*x8;

result[0] = x0*x4 + x9;
result[1] = x4*x5 + x9;
result[2] = x4*x6 + x9;
result[3] = x4*x7 + x9;
}
        
static void coder_d3gdn2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3 + n4;
    double x2 = 1.0/x1;
    double x3 = x0*x2;
    double x4 = pow(x1, -2);
    double x5 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x6 = x4*x5;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = (*endmember[3].dmu0dP)(T, P);
    double x11 = 2.0*n1*x0 + 2.0*n2*x8 + 2.0*n3*x9 + 2.0*n4*x10;
    double x12 = -x11*x4 + x11*x5/((x1)*(x1)*(x1));
    double x13 = x12 + 1.0*x3 - 1.0*x7;
    double x14 = 1.0*x2;
    double x15 = 1.0*x6;
    double x16 = x14*x8 - x15*x8;
    double x17 = x14*x9 - x15*x9;
    double x18 = x10*x14 - x10*x15;
    double x19 = 2.0*x2;
    double x20 = 2.0*x6;
    double x21 = x12 + x16;

result[0] = x12 + 2.0*x3 - 2.0*x7;
result[1] = x13 + x16;
result[2] = x13 + x17;
result[3] = x13 + x18;
result[4] = x12 + x19*x8 - x20*x8;
result[5] = x17 + x21;
result[6] = x18 + x21;
result[7] = x12 + x19*x9 - x20*x9;
result[8] = x12 + x17 + x18;
result[9] = x10*x19 - x10*x20 + x12;
}
        
static void coder_d4gdn3dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];

    double x0 = (*endmember[0].dmu0dP)(T, P);
    double x1 = n1 + n2 + n3 + n4;
    double x2 = pow(x1, -2);
    double x3 = x0*x2;
    double x4 = pow(x1, -3);
    double x5 = 6.0*x4;
    double x6 = 1.0*n1 + 1.0*n2 + 1.0*n3 + 1.0*n4;
    double x7 = x0*x6;
    double x8 = (*endmember[1].dmu0dP)(T, P);
    double x9 = (*endmember[2].dmu0dP)(T, P);
    double x10 = (*endmember[3].dmu0dP)(T, P);
    double x11 = n1*x0 + n2*x8 + n3*x9 + n4*x10;
    double x12 = x11*x5 - 6.0*x11*x6/((x1)*(x1)*(x1)*(x1));
    double x13 = 4.0*x4;
    double x14 = x12 + x13*x7 - 4.0*x3;
    double x15 = 2.0*x2;
    double x16 = 2.0*x4;
    double x17 = x6*x8;
    double x18 = -x15*x8 + x16*x17;
    double x19 = x16*x6;
    double x20 = -x15*x9 + x19*x9;
    double x21 = -x10*x15 + x10*x19;
    double x22 = x12 + x16*x7 - 2.0*x3;
    double x23 = 4.0*x2;
    double x24 = x13*x17 - x23*x8;
    double x25 = x18 + x22;
    double x26 = x13*x6;
    double x27 = -x23*x9 + x26*x9;
    double x28 = x20 + x21;
    double x29 = -x10*x23 + x10*x26;
    double x30 = 6.0*x2;
    double x31 = x12 + x24;
    double x32 = x12 + x18;
    double x33 = x5*x6;

result[0] = x12 - 6.0*x3 + x5*x7;
result[1] = x14 + x18;
result[2] = x14 + x20;
result[3] = x14 + x21;
result[4] = x22 + x24;
result[5] = x20 + x25;
result[6] = x21 + x25;
result[7] = x22 + x27;
result[8] = x22 + x28;
result[9] = x22 + x29;
result[10] = x12 + x17*x5 - x30*x8;
result[11] = x20 + x31;
result[12] = x21 + x31;
result[13] = x27 + x32;
result[14] = x28 + x32;
result[15] = x29 + x32;
result[16] = x12 - x30*x9 + x33*x9;
result[17] = x12 + x21 + x27;
result[18] = x12 + x20 + x29;
result[19] = -x10*x30 + x10*x33 + x12;
}
        
static double coder_d2gdt2(double T, double P, double n[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double result;
    

result = 1.0*(n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d2mu0dT2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dT2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dT2)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dT2)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d2mu0dTdP)(T, P) + 1.0*n2*(*endmember[1].d2mu0dTdP)(T, P) + 1.0*n3*(*endmember[2].d2mu0dTdP)(T, P) + 1.0*n4*(*endmember[3].d2mu0dTdP)(T, P);
    return result;
}
        
static void coder_d3gdndtdp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d2mu0dTdP)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dTdP)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dTdP)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dTdP)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d2mu0dP2)(T, P) + n2*(*endmember[1].d2mu0dP2)(T, P) + n3*(*endmember[2].d2mu0dP2)(T, P) + n4*(*endmember[3].d2mu0dP2)(T, P));
    return result;
}
        
static void coder_d3gdndp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d2mu0dP2)(T, P);
result[1] = 1.0*(*endmember[1].d2mu0dP2)(T, P);
result[2] = 1.0*(*endmember[2].d2mu0dP2)(T, P);
result[3] = 1.0*(*endmember[3].d2mu0dP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d3mu0dT3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT3)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dT3)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dT2dP)(T, P) + n2*(*endmember[1].d3mu0dT2dP)(T, P) + n3*(*endmember[2].d3mu0dT2dP)(T, P) + n4*(*endmember[3].d3mu0dT2dP)(T, P));
    return result;
}
        
static void coder_d4gdndt2dp(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d3mu0dT2dP)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dT2dP)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dT2dP)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dT2dP)(T, P);
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
    

result = 1.0*n1*(*endmember[0].d3mu0dTdP2)(T, P) + 1.0*n2*(*endmember[1].d3mu0dTdP2)(T, P) + 1.0*n3*(*endmember[2].d3mu0dTdP2)(T, P) + 1.0*n4*(*endmember[3].d3mu0dTdP2)(T, P);
    return result;
}
        
static void coder_d4gdndtdp2(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d3mu0dTdP2)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dTdP2)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dTdP2)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dTdP2)(T, P);
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
    

result = 1.0*(n1*(*endmember[0].d3mu0dP3)(T, P) + n2*(*endmember[1].d3mu0dP3)(T, P) + n3*(*endmember[2].d3mu0dP3)(T, P) + n4*(*endmember[3].d3mu0dP3)(T, P));
    return result;
}
        
static void coder_d4gdndp3(double T, double P, double n[4], double result[4]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];


result[0] = 1.0*(*endmember[0].d3mu0dP3)(T, P);
result[1] = 1.0*(*endmember[1].d3mu0dP3)(T, P);
result[2] = 1.0*(*endmember[2].d3mu0dP3)(T, P);
result[3] = 1.0*(*endmember[3].d3mu0dP3)(T, P);
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

