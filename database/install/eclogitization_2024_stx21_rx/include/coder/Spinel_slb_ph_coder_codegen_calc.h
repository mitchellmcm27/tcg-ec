#include <math.h>


static double coder_g(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = log(x1*(0.25*n1 + 0.875*n2));

result = 1.0*x1*(-500.0*n1*n2 + 1.0*x0*(8.3144626181532395*T*(2.0*n1*x2 + 4.0*n1*log(n1*x1) - 4.3287821201550702*n1 + 7.0*n2*x2 + 4.0*n2*log(n2*x1) - 4.3287821201550702*n2) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P)));
    return result;
}
        
static void coder_dgdn(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -2);
    double x2 = 500.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = n1*x5;
    double x7 = 4.0*log(x6);
    double x8 = n2*x5;
    double x9 = 4.0*log(x8);
    double x10 = 0.25*n1 + 0.875*n2;
    double x11 = log(x10*x5);
    double x12 = 2.0*x11;
    double x13 = 7.0*x11;
    double x14 = 8.3144626181532395*T;
    double x15 = x14*(n1*x12 + n1*x7 - 4.3287821201550702*n1 + n2*x13 + n2*x9 - 4.3287821201550702*n2);
    double x16 = -1.0*x1*(-n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x15));
    double x17 = 1.0*n1;
    double x18 = 1.0*n2;
    double x19 = x17 + x18;
    double x20 = 4.0*x0;
    double x21 = -x1*x10;
    double x22 = x0/x10;
    double x23 = x22*(x21 + 0.25*x5);
    double x24 = 2.0*n1;
    double x25 = 7.0*n2;
    double x26 = x15 + x17*x3 + x18*x4;
    double x27 = 1.0*x5;
    double x28 = x22*(x21 + 0.875*x5);

result[0] = x16 + x27*(x19*(x14*(x12 + x20*(-n1*x1 + x5) + x23*x24 + x23*x25 + x7 - 4.0*x8 - 4.3287821201550702) + x3) - x2 + x26);
result[1] = x16 + x27*(-500.0*n1 + x19*(x14*(x13 + x20*(-n2*x1 + x5) + x24*x28 + x25*x28 - 4.0*x6 + x9 - 4.3287821201550702) + x4) + x26);
}
        
static void coder_d2gdn2(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -3);
    double x2 = 500.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = n1*x5;
    double x7 = 4.0*log(x6);
    double x8 = log(n2*x5);
    double x9 = 4.0*n2;
    double x10 = 0.25*n1 + 0.875*n2;
    double x11 = log(x10*x5);
    double x12 = 2.0*x11;
    double x13 = 7.0*x11;
    double x14 = 8.3144626181532395*T;
    double x15 = x14*(n1*x12 + n1*x7 - 4.3287821201550702*n1 + n2*x13 - 4.3287821201550702*n2 + x8*x9);
    double x16 = 2.0*x1*(-n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x15));
    double x17 = pow(x0, -2);
    double x18 = 1.0*n1;
    double x19 = 1.0*n2;
    double x20 = x18 + x19;
    double x21 = 4.0*x5;
    double x22 = n1*x17;
    double x23 = 4.0*x0;
    double x24 = x23*(-x22 + x5);
    double x25 = 1.0/x10;
    double x26 = -x10*x17;
    double x27 = x26 + 0.25*x5;
    double x28 = x25*x27;
    double x29 = 2.0*n1;
    double x30 = x28*x29;
    double x31 = 7.0*x28;
    double x32 = n2*x31;
    double x33 = T*(-n2*x21 + x0*x30 + x0*x32 + x12 + x24 + x7 - 4.3287821201550702);
    double x34 = 8.3144626181532395*x33;
    double x35 = x15 + x18*x3 + x19*x4;
    double x36 = x17*(-x2 + x20*(x3 + x34) + x35);
    double x37 = 4.0*n1 + x9;
    double x38 = -2*x17;
    double x39 = 2*x1;
    double x40 = n1*x39;
    double x41 = x10*x39;
    double x42 = x0*x25;
    double x43 = x42*(-0.5*x17 + x41);
    double x44 = 7.0*n2;
    double x45 = x0/((x10)*(x10));
    double x46 = x27*x45;
    double x47 = n1*x46;
    double x48 = n2*x46;
    double x49 = x17*x9;
    double x50 = 4.0*x22;
    double x51 = x30 + x32 + x49 - x50;
    double x52 = x14*x20;
    double x53 = 1.0*x5;
    double x54 = x23*(-n2*x17 + x5);
    double x55 = x26 + 0.875*x5;
    double x56 = x25*x55;
    double x57 = x0*x56;
    double x58 = x13 + x29*x57 + x44*x57 + x54 - 4.0*x6 + 4.0*x8 - 4.3287821201550702;
    double x59 = x14*x58;
    double x60 = x17*(-500.0*n1 + x20*(x4 + x59) + x35);
    double x61 = x42*(-1.125*x17 + x41);
    double x62 = x45*x55;
    double x63 = x42*(-1.75*x17 + x41);

result[0] = x16 - 2.0*x36 + x53*(2.0*x3 + 16.628925236306479*x33 + x52*(x21 + x23*x28 + x29*x43 + x37*(x38 + x40) + x43*x44 - 0.5*x47 - 1.75*x48 + x51 + x24/n1));
result[1] = x16 - 1.0*x36 + x53*(1.0*x3 + x34 + 1.0*x4 + x52*(x0*x31 - x21 + x29*x61 + x37*(-x17 + x40) + x44*x61 - 1.75*x47 - 6.125*x48 + x51 + 2.0*x57) + x59 - 500.0) - 1.0*x60;
result[2] = x16 + x53*(16.628925236306479*T*x58 + 2.0*x4 + x52*(-1.75*n1*x62 - 6.125*n2*x62 + x21 + x29*x56 + x29*x63 + x37*(n2*x39 + x38) + x44*x56 + x44*x63 - x49 + x50 + 14.0*x57 + x54/n2)) - 2.0*x60;
}
        
static void coder_d3gdn3(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = pow(x0, -4);
    double x2 = 500.0*n2;
    double x3 = (*endmember[0].mu0)(T, P);
    double x4 = (*endmember[1].mu0)(T, P);
    double x5 = 1.0/x0;
    double x6 = log(n1*x5);
    double x7 = 4.0*n1;
    double x8 = log(n2*x5);
    double x9 = 4.0*n2;
    double x10 = 0.25*n1;
    double x11 = 0.875*n2;
    double x12 = x10 + x11;
    double x13 = log(x12*x5);
    double x14 = 2.0*x13;
    double x15 = 7.0*x13;
    double x16 = 8.3144626181532395*T;
    double x17 = x16*(n1*x14 - 4.3287821201550702*n1 + n2*x15 - 4.3287821201550702*n2 + x6*x7 + x8*x9);
    double x18 = -6.0*x1*(-n1*x2 + 1.0*x0*(n1*x3 + n2*x4 + x17));
    double x19 = pow(x0, -3);
    double x20 = 1.0*n1;
    double x21 = 1.0*n2;
    double x22 = x20 + x21;
    double x23 = 4.0*x5;
    double x24 = pow(x0, -2);
    double x25 = -4.0*n1*x24 + 4.0*x5;
    double x26 = x0*x25;
    double x27 = 1.0/x12;
    double x28 = -x12*x24;
    double x29 = x28 + 0.25*x5;
    double x30 = x27*x29;
    double x31 = 2.0*n1;
    double x32 = x30*x31;
    double x33 = 7.0*x30;
    double x34 = n2*x33;
    double x35 = -n2*x23 + x0*x32 + x0*x34 + x14 + x26 + 4.0*x6 - 4.3287821201550702;
    double x36 = x16*x35;
    double x37 = x17 + x20*x3 + x21*x4;
    double x38 = -x2 + x22*(x3 + x36) + x37;
    double x39 = x19*x38;
    double x40 = 16.628925236306479*T;
    double x41 = x7 + x9;
    double x42 = -2*x24;
    double x43 = 2*x19;
    double x44 = n1*x43;
    double x45 = x42 + x44;
    double x46 = 1.0/n1;
    double x47 = x25*x46;
    double x48 = 4.0*x0;
    double x49 = x12*x43;
    double x50 = -0.5*x24 + x49;
    double x51 = x27*x50;
    double x52 = x0*x51;
    double x53 = 7.0*n2;
    double x54 = pow(x12, -2);
    double x55 = x29*x54;
    double x56 = x0*x55;
    double x57 = n1*x56;
    double x58 = n2*x56;
    double x59 = x24*x9;
    double x60 = x24*x7;
    double x61 = x32 + x34 + x59 - x60;
    double x62 = T*(x0*x47 + x23 + x30*x48 + x31*x52 + x41*x45 + x52*x53 - 0.5*x57 - 1.75*x58 + x61);
    double x63 = 8.3144626181532395*x62;
    double x64 = x24*(x22*x63 + 2.0*x3 + x35*x40);
    double x65 = -16.0*x24;
    double x66 = 6*x19;
    double x67 = 6*x1;
    double x68 = -n1*x67;
    double x69 = x46*x48;
    double x70 = n2*x51;
    double x71 = n2*x55;
    double x72 = x0/((x12)*(x12)*(x12));
    double x73 = x29*x72;
    double x74 = -x12*x67;
    double x75 = x0*x27;
    double x76 = x75*(1.5*x19 + x74);
    double x77 = x0*x54;
    double x78 = x50*x77;
    double x79 = n2*x78;
    double x80 = 16.0*x19;
    double x81 = 8.0*x19;
    double x82 = n1*x80 - n2*x81;
    double x83 = x47 + x82;
    double x84 = x16*x22;
    double x85 = 1.0*x5;
    double x86 = -x24 + x44;
    double x87 = x28 + 0.875*x5;
    double x88 = x27*x87;
    double x89 = 2.0*x88;
    double x90 = x0*x89;
    double x91 = -1.125*x24 + x49;
    double x92 = x27*x91;
    double x93 = x31*x92;
    double x94 = x53*x92;
    double x95 = x0*x33 + x0*x93 + x0*x94 - x23 + x41*x86 - 1.75*x57 - 6.125*x58 + x61 + x90;
    double x96 = x40*x95;
    double x97 = n1*x55;
    double x98 = x75*(2.75*x19 + x74);
    double x99 = n1*x73;
    double x100 = n2*x73;
    double x101 = x77*x91;
    double x102 = -4.0*n2*x24 + 4.0*x5;
    double x103 = x0*x102;
    double x104 = x53*x88;
    double x105 = -n1*x23 + n1*x90 + x0*x104 + x103 + x15 + 4.0*x8 - 4.3287821201550702;
    double x106 = x105*x16;
    double x107 = -500.0*n1 + x22*(x106 + x4) + x37;
    double x108 = x107*x19;
    double x109 = 4.0*x19;
    double x110 = x18 - 2.0*x24*(x106 + 1.0*x3 + x36 + 1.0*x4 + x84*x95 - 500.0);
    double x111 = 1.0/n2;
    double x112 = n2*x43 + x42;
    double x113 = 14.0*x0;
    double x114 = x54*x87;
    double x115 = x0*x114;
    double x116 = 1.75*x115;
    double x117 = -1.75*x24 + x49;
    double x118 = x117*x27;
    double x119 = x0*x118;
    double x120 = -n1*x116 + n1*x89 - 6.125*n2*x115 + x103*x111 + x104 + x112*x41 + x113*x88 + x119*x31 + x119*x53 + x23 - x59 + x60;
    double x121 = x120*x16;
    double x122 = 14.0*n2;
    double x123 = x75*(x109 + x74);
    double x124 = 3.5*n1;
    double x125 = 12.25*n2;
    double x126 = x24*(x105*x40 + x121*x22 + 2.0*x4);
    double x127 = x72*x87;
    double x128 = x117*x77;
    double x129 = x75*(5.25*x19 + x74);

result[0] = x18 + 6.0*x39 - 3.0*x64 + x85*(24.943387854459719*x62 + x84*(x10*x73 + x11*x73 - x20*x55 - x20*x78 + 6.0*x30 + x31*x76 + x41*(x66 + x68) + x45*x69 + x51*x7 + 6.0*x52 + x53*x76 - 1.5*x56 + x65 + 14.0*x70 - 3.5*x71 - 3.5*x79 + x83 - x26/((n1)*(n1))));
result[1] = 2.0*x108 + x109*x38 + x110 - 1.0*x64 + x85*(x63 + x84*(-0.5*n1*x101 - 1.75*n1*x78 - 1.75*n2*x101 + 3.0625*x100 - 8.0*x24 + 11.0*x30 + x31*x51 + x31*x98 + x41*(4*x19 + x68) + x48*x92 + 7.0*x52 + x53*x98 - 5.25*x56 + x69*x86 + 7.0*x70 - 7.875*x71 - 6.125*x79 + x83 + x93 + x94 - 2.25*x97 + 0.875*x99) + x96);
result[2] = x107*x109 + x110 - 1.0*x126 + 2.0*x39 + x85*(x121 + x84*(10.71875*x100 - x101*x124 - x101*x125 + x113*x92 - x116 + 2.0*x119 + x122*x92 + x123*x31 + x123*x53 + 4.0*x24 + 14.0*x30 + x41*(x43 + x68) - 12.25*x56 + x7*x92 - 12.25*x71 + x82 + x89 - 3.5*x97 + 3.0625*x99) + x96);
result[3] = 6.0*x108 - 3.0*x126 + x18 + x85*(24.943387854459719*T*x120 + x84*(3.0625*n1*x127 - n1*x81 + 10.71875*n2*x127 + n2*x80 + x102*x111 + x111*x112*x48 - x114*x124 - x114*x125 - 18.375*x115 + x118*x122 + x118*x7 + 21.0*x119 - x124*x128 - x125*x128 + x129*x31 + x129*x53 + x41*(-n2*x67 + x66) + x65 + 21.0*x88 - x103/((n2)*(n2))));
}
        
static double coder_dgdt(double T, double P, double n[2]) {
    double n1 = n[0];
    double n2 = n[1];
    double result;
    
    double x0 = 1.0/(n1 + n2);
    double x1 = log(x0*(0.25*n1 + 0.875*n2));

result = 1.0*x0*(1.0*n1 + 1.0*n2)*(16.628925236306479*n1*x1 + 33.257850472612958*n1*log(n1*x0) + n1*(*endmember[0].dmu0dT)(T, P) - 35.991497120159458*n1 + 58.201238327072673*n2*x1 + 33.257850472612958*n2*log(n2*x0) + n2*(*endmember[1].dmu0dT)(T, P) - 35.991497120159458*n2);
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
    double x5 = 33.257850472612958*log(x4);
    double x6 = 0.25*n1 + 0.875*n2;
    double x7 = log(x2*x6);
    double x8 = 16.628925236306479*x7;
    double x9 = pow(x1, -2);
    double x10 = 33.257850472612958*x1;
    double x11 = -x6*x9;
    double x12 = x1/x6;
    double x13 = x12*(x11 + 0.25*x2);
    double x14 = 16.628925236306479*n1;
    double x15 = 58.201238327072673*n2;
    double x16 = 1.0*n1 + 1.0*n2;
    double x17 = 1.0*x2;
    double x18 = x16*x17;
    double x19 = (*endmember[1].dmu0dT)(T, P);
    double x20 = 33.257850472612958*log(x3);
    double x21 = 58.201238327072673*x7;
    double x22 = n1*x0 + n1*x5 + n1*x8 - 35.991497120159458*n1 + n2*x19 + n2*x20 + n2*x21 - 35.991497120159458*n2;
    double x23 = -1.0*x16*x22*x9 + x17*x22;
    double x24 = x12*(x11 + 0.875*x2);

result[0] = x18*(x0 + x10*(-n1*x9 + x2) + x13*x14 + x13*x15 - 33.257850472612958*x3 + x5 + x8 - 35.991497120159458) + x23;
result[1] = x18*(x10*(-n2*x9 + x2) + x14*x24 + x15*x24 + x19 + x20 + x21 - 33.257850472612958*x4 - 35.991497120159458) + x23;
}
        
static void coder_d3gdn2dt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = (*endmember[0].dmu0dT)(T, P);
    double x3 = n2*x1;
    double x4 = n1*x1;
    double x5 = 33.257850472612958*log(x4);
    double x6 = 0.25*n1 + 0.875*n2;
    double x7 = log(x1*x6);
    double x8 = 16.628925236306479*x7;
    double x9 = pow(x0, -2);
    double x10 = n1*x9;
    double x11 = 33.257850472612958*x0;
    double x12 = x11*(x1 - x10);
    double x13 = 1.0/x6;
    double x14 = -x6*x9;
    double x15 = 0.25*x1 + x14;
    double x16 = x13*x15;
    double x17 = 16.628925236306479*n1;
    double x18 = x16*x17;
    double x19 = 58.201238327072673*x16;
    double x20 = n2*x19;
    double x21 = x0*x18 + x0*x20 + x12 + x2 - 33.257850472612958*x3 + x5 + x8 - 35.991497120159458;
    double x22 = x1*x21;
    double x23 = 33.257850472612958*x1;
    double x24 = 33.257850472612958*n2;
    double x25 = 33.257850472612958*n1 + x24;
    double x26 = -2*x9;
    double x27 = pow(x0, -3);
    double x28 = 2*x27;
    double x29 = n1*x28;
    double x30 = 58.201238327072673*n2;
    double x31 = x28*x6;
    double x32 = x0*x13;
    double x33 = x32*(x31 - 0.5*x9);
    double x34 = pow(x6, -2);
    double x35 = x0*x15*x34;
    double x36 = 14.550309581768168*x35;
    double x37 = x24*x9;
    double x38 = 33.257850472612958*x10;
    double x39 = x18 + x20 + x37 - x38;
    double x40 = 1.0*n1 + 1.0*n2;
    double x41 = 1.0*x40;
    double x42 = x1*x41;
    double x43 = x21*x9;
    double x44 = 2.0*x40;
    double x45 = (*endmember[1].dmu0dT)(T, P);
    double x46 = log(x3);
    double x47 = 58.201238327072673*x7;
    double x48 = 2.0*n1*x2 + 2.0*n1*x5 + 2.0*n1*x8 - 71.982994240318916*n1 + 2.0*n2*x45 + 2.0*n2*x47 - 71.982994240318916*n2 + 2.0*x24*x46;
    double x49 = x27*x40*x48 - x48*x9;
    double x50 = x11*(-n2*x9 + x1);
    double x51 = 0.875*x1 + x14;
    double x52 = x13*x51;
    double x53 = x0*x52;
    double x54 = x17*x53 + x30*x53 - 33.257850472612958*x4 + x45 + 33.257850472612958*x46 + x47 + x50 - 35.991497120159458;
    double x55 = x1*x54;
    double x56 = x32*(x31 - 1.125*x9);
    double x57 = 50.926083536188585*n2;
    double x58 = x54*x9;
    double x59 = x0*x34*x51;
    double x60 = x32*(x31 - 1.75*x9);

result[0] = 2.0*x22 + x42*(-4.1572313090766198*n1*x35 - n2*x36 + x11*x16 + x17*x33 + x23 + x25*(x26 + x29) + x30*x33 + x39 + x12/n1) - x43*x44 + x49;
result[1] = 1.0*x22 - x41*x43 - x41*x58 + x42*(-n1*x36 + x0*x19 + x17*x56 - x23 + x25*(x29 - x9) + x30*x56 - x35*x57 + x39 + 16.628925236306479*x53) + x49 + 1.0*x55;
result[2] = x42*(-14.550309581768168*n1*x59 + x17*x52 + x17*x60 + x23 + x25*(n2*x28 + x26) + x30*x52 + x30*x60 - x37 + x38 + 116.40247665414535*x53 - x57*x59 + x50/n2) - x44*x58 + x49 + 2.0*x55;
}
        
static void coder_d4gdn3dt(double T, double P, double n[2], double result[2]) {
    double n1 = n[0];
    double n2 = n[1];

    double x0 = n1 + n2;
    double x1 = 1.0/x0;
    double x2 = 33.257850472612958*x1;
    double x3 = 33.257850472612958*n1;
    double x4 = 33.257850472612958*n2;
    double x5 = x3 + x4;
    double x6 = pow(x0, -2);
    double x7 = -2*x6;
    double x8 = pow(x0, -3);
    double x9 = 2*x8;
    double x10 = n1*x9;
    double x11 = x10 + x7;
    double x12 = 1.0/n1;
    double x13 = -33.257850472612958*n1*x6 + 33.257850472612958*x1;
    double x14 = x0*x13;
    double x15 = 0.25*n1 + 0.875*n2;
    double x16 = 1.0/x15;
    double x17 = -x15*x6;
    double x18 = 0.25*x1 + x17;
    double x19 = x16*x18;
    double x20 = 33.257850472612958*x0;
    double x21 = x15*x9;
    double x22 = x21 - 0.5*x6;
    double x23 = x16*x22;
    double x24 = n2*x23;
    double x25 = 58.201238327072673*x24;
    double x26 = 16.628925236306479*n1;
    double x27 = x0*x23;
    double x28 = pow(x15, -2);
    double x29 = x18*x28;
    double x30 = n2*x29;
    double x31 = 14.550309581768168*x0;
    double x32 = n1*x29;
    double x33 = 4.1572313090766198*x0;
    double x34 = x4*x6;
    double x35 = x3*x6;
    double x36 = 58.201238327072673*x19;
    double x37 = n2*x36;
    double x38 = x19*x26;
    double x39 = x34 - x35 + x37 + x38;
    double x40 = x0*x25 + x11*x5 + x12*x14 + x19*x20 + x2 + x26*x27 - x30*x31 - x32*x33 + x39;
    double x41 = x1*x40;
    double x42 = (*endmember[0].dmu0dT)(T, P);
    double x43 = log(n1*x1);
    double x44 = log(x1*x15);
    double x45 = 16.628925236306479*x44;
    double x46 = -n2*x2 + x0*x37 + x0*x38 + x14 + x42 + 33.257850472612958*x43 + x45 - 35.991497120159458;
    double x47 = x46*x6;
    double x48 = -133.03140189045183*x6;
    double x49 = 6*x8;
    double x50 = pow(x0, -4);
    double x51 = 6*x50;
    double x52 = -n1*x51;
    double x53 = x12*x20;
    double x54 = x0*x29;
    double x55 = n2*x0;
    double x56 = pow(x15, -3);
    double x57 = x18*x56;
    double x58 = 7.2751547908840841*x57;
    double x59 = 58.201238327072673*n2;
    double x60 = -x15*x51;
    double x61 = x0*x16;
    double x62 = x61*(x60 + 1.5*x8);
    double x63 = n1*x0;
    double x64 = x57*x63;
    double x65 = x22*x28;
    double x66 = x55*x65;
    double x67 = 133.03140189045183*x8;
    double x68 = 66.515700945225916*x8;
    double x69 = n1*x67 - n2*x68;
    double x70 = x12*x13 + x69;
    double x71 = 1.0*n1 + 1.0*n2;
    double x72 = 1.0*x71;
    double x73 = x1*x72;
    double x74 = 6.0*x8;
    double x75 = x46*x71;
    double x76 = x40*x6;
    double x77 = 3.0*x71;
    double x78 = (*endmember[1].dmu0dT)(T, P);
    double x79 = log(n2*x1);
    double x80 = 58.201238327072673*x44;
    double x81 = n1*x42 + n1*x45 - 35.991497120159458*n1 + n2*x78 + n2*x80 - 35.991497120159458*n2 + x3*x43 + x4*x79;
    double x82 = -6.0*x50*x71*x81 + x74*x81;
    double x83 = -33.257850472612958*n2*x6 + 33.257850472612958*x1;
    double x84 = x0*x83;
    double x85 = 0.875*x1 + x17;
    double x86 = x16*x85;
    double x87 = 16.628925236306479*x86;
    double x88 = x0*x87;
    double x89 = x59*x86;
    double x90 = -n1*x2 + n1*x88 + x0*x89 + x78 + 33.257850472612958*x79 + x80 + x84 - 35.991497120159458;
    double x91 = x6*x90;
    double x92 = x21 - 1.125*x6;
    double x93 = x16*x92;
    double x94 = x59*x93;
    double x95 = x26*x93;
    double x96 = x10 - x6;
    double x97 = x61*(x60 + 2.75*x8);
    double x98 = x55*x57;
    double x99 = x28*x92;
    double x100 = 2.0*x8;
    double x101 = x71*x90;
    double x102 = 4.0*x8;
    double x103 = -101.85216707237717*x0*x30 + 2.0*x0*x36 + 2.0*x0*x94 + 2.0*x0*x95 - 2.0*x2 - 2.0*x31*x32 + 2.0*x39 + 2.0*x5*x96 + 2.0*x88;
    double x104 = x1*x103 - x103*x6*x71 + x82;
    double x105 = 1.0/n2;
    double x106 = n2*x9 + x7;
    double x107 = 116.40247665414535*x0;
    double x108 = x28*x85;
    double x109 = x108*x31;
    double x110 = x21 - 1.75*x6;
    double x111 = x110*x16;
    double x112 = x0*x111;
    double x113 = x0*x108;
    double x114 = -n1*x109 + n1*x87 - 50.926083536188585*n2*x113 + x105*x84 + x106*x5 + x107*x86 + x112*x26 + x112*x59 + x2 - x34 + x35 + x89;
    double x115 = x1*x114;
    double x116 = 116.40247665414535*n2;
    double x117 = x61*(x102 + x60);
    double x118 = 29.100619163536336*n1;
    double x119 = x0*x99;
    double x120 = 101.85216707237717*n2;
    double x121 = x114*x6;
    double x122 = x56*x85;
    double x123 = x0*x110*x28;
    double x124 = x61*(x60 + 5.25*x8);

result[0] = 3.0*x41 - 6.0*x47 + x73*(x11*x53 + 49.886775708919437*x19 + x23*x3 + 116.40247665414535*x24 + x26*x62 + 49.886775708919437*x27 - 29.100619163536336*x30 - 8.3144626181532395*x32 + x48 + x5*(x49 + x52) - 12.471693927229859*x54 + x55*x58 + x59*x62 - 8.3144626181532395*x63*x65 + 2.0786156545383099*x64 - 29.100619163536336*x66 + x70 - x14/((n1)*(n1))) + x74*x75 - x76*x77 + x82;
result[1] = x100*x101 + x102*x75 + x104 + 1.0*x41 - 4.0*x47 - x72*x76 + x73*(-n1*x31*x65 - n1*x33*x99 - n2*x31*x99 + 91.459088799685631*x19 + x20*x93 + x23*x26 + x25 + x26*x97 + 58.201238327072673*x27 - 65.476393117956746*x30 - 18.70754089084479*x32 + x5*(x52 + 4*x8) + x53*x96 - 43.650928745304505*x54 + x58*x63 + x59*x97 - 66.515700945225916*x6 - 50.926083536188585*x66 + x70 + x94 + x95 + 25.463041768094293*x98) - 2.0*x91;
result[2] = x100*x75 + x101*x102 + x104 + 1.0*x115 - x121*x72 - 2.0*x47 + x73*(x107*x93 - x109 + 16.628925236306479*x112 + x116*x93 + x117*x26 + x117*x59 - x118*x119 - x119*x120 + 116.40247665414535*x19 + x3*x93 - 101.85216707237717*x30 - 29.100619163536336*x32 + x5*(x52 + x9) - 101.85216707237717*x54 + 33.257850472612958*x6 + 25.463041768094293*x64 + x69 + x87 + 89.120646188330028*x98) - 4.0*x91;
result[3] = x101*x74 + 3.0*x115 - x121*x77 + x73*(-n1*x68 + n2*x67 + x105*x106*x20 + x105*x83 - x108*x118 - x108*x120 + x111*x116 + x111*x3 + 174.60371498121802*x112 - 152.77825060856577*x113 - x118*x123 - x120*x123 + 89.120646188330028*x122*x55 + 25.463041768094293*x122*x63 + x124*x26 + x124*x59 + x48 + x5*(-n2*x51 + x49) + 174.60371498121802*x86 - x84/((n2)*(n2))) + x82 - 6.0*x91;
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

