#include <math.h>


static double coder_g(double T, double P, double n[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];
    double result;
    
    double x0 = 40.14405*T;
    double x1 = (1.0/3.0)*n2;
    double x2 = fmin(4, 0.27336948519328247*sqrt(-2.6762700000000001*T + 13.381350000000001));
    double x3 = 1.0*n5;
    double x4 = 1.0*n2;
    double x5 = 1.0*n1;
    double x6 = 1.0*n3;
    double x7 = x5 + x6;
    double x8 = x3 + x4 + x7;
    double x9 = 324113400.0*n3 + 530712000.0*n4 + 318864600.0*n5;
    double x10 = n1 + n3;
    double x11 = n4 + n5;
    double x12 = 1.0/(n2 + x10 + x11);
    double x13 = 1.0*n4;

result = 8.3144626181532395*T*(x13*(log(n4*x12) - 0.69314718055994495) + x3*log(n5*x12) + x4*log(n2*x12) + x6*log(n3*x12) + x7*log(x10*x12) + (x13 + x3)*log(x11*x12) + (x13 + x4 + x5)*log(x12*(n1 + n2 + n4)) + (2.0*n1 + 2.0*n2 + 2.0*n3 + 2.0*n5 + x13)*log(x12*(0.5*n4 + x8))) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P) + n4*(*endmember[3].mu0)(T, P) + n5*(*endmember[4].mu0)(T, P) + ((T >= 5.0) ? (
   -x1*(x0 - 133.8135)
)
: (
   x1*(66.906750000000002*((x2)*(x2)*(x2)) + (x0 - 200.72024999999999)*(x2 - 1) - 66.906750000000002)
)) + 3.8103947568968139e-5*(n1*x9 + n2*x9 + 729.0*n3*(444600.0*n1 + 444600.0*n2 + 1682800.0*n4 + 828000.0*n5) + 224.0*n4*(2369250.0*n1 + 2369250.0*n2 + 5476612.5*n3 + 911250.0*n5) + 729.0*n5*(437400.0*n1 + 437400.0*n2 + 828000.0*n3 + 280000.0*n4))/(3.5*n4 + x8);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n3;
    double x1 = n4 + n5;
    double x2 = n2 + x0 + x1;
    double x3 = 1.0/x2;
    double x4 = 1.0*n5;
    double x5 = 1.0*n2;
    double x6 = 1.0*n1;
    double x7 = 1.0*n3;
    double x8 = x6 + x7;
    double x9 = x4 + x5 + x8;
    double x10 = 0.5*n4 + x9;
    double x11 = log(x10*x3);
    double x12 = 2.0*x11;
    double x13 = -x3*x5;
    double x14 = -x3*x7;
    double x15 = 1.0*n4;
    double x16 = -x15*x3;
    double x17 = pow(x2, -2);
    double x18 = -x10*x17;
    double x19 = x2*(2.0*n1 + 2.0*n2 + 2.0*n3 + 2.0*n5 + x15)/x10;
    double x20 = x19*(x18 + 1.0*x3);
    double x21 = x12 + x13 + x14 + x16 + x20;
    double x22 = x15 + x4;
    double x23 = -x22*x3;
    double x24 = n1 + n2 + n4;
    double x25 = -x3*x4;
    double x26 = x15 + x5 + x6;
    double x27 = x2*x26*(-x17*x24 + x3)/x24 + x25 + 1.0*log(x24*x3);
    double x28 = x23 + x27;
    double x29 = 1.0*log(x0*x3) + x2*x8*(-x0*x17 + x3)/x0;
    double x30 = 8.3144626181532395*T;
    double x31 = 3.5*n4 + x9;
    double x32 = 3.8103947568968139e-5/x31;
    double x33 = 324113400.0*n3 + 530712000.0*n4 + 318864600.0*n5;
    double x34 = (n1*x33 + n2*x33 + 729.0*n3*(444600.0*n1 + 444600.0*n2 + 1682800.0*n4 + 828000.0*n5) + 224.0*n4*(2369250.0*n1 + 2369250.0*n2 + 5476612.5*n3 + 911250.0*n5) + 729.0*n5*(437400.0*n1 + 437400.0*n2 + 828000.0*n3 + 280000.0*n4))/((x31)*(x31));
    double x35 = -3.8103947568968139e-5*x34;
    double x36 = x32*(648226800.0*n3 + 1061424000.0*n4 + 637729200.0*n5) + x35;
    double x37 = 1.0*x2;
    double x38 = -x3*x8;
    double x39 = x14 + x38;
    double x40 = x12 + x16 + x20;
    double x41 = fmin(4, 0.27336948519328247*sqrt(-2.6762700000000001*T + 13.381350000000001));
    double x42 = -x26*x3;
    double x43 = 1.0*log(x1*x3) + x2*x22*(-x1*x17 + x3)/x1;

result[0] = x30*(x21 + x28 + x29) + x36 + (*endmember[0].mu0)(T, P);
result[1] = x30*(x28 + x37*(-n2*x17 + x3) + x39 + x40 + 1.0*log(n2*x3)) + x36 + ((T >= 5.0) ? (
   -13.381349999999999*T + 44.604500000000002
)
: (
   22.302250000000001*((x41)*(x41)*(x41)) + (1.0/3.0)*(40.14405*T - 200.72024999999999)*(x41 - 1) - 22.302250000000001
)) + (*endmember[1].mu0)(T, P);
result[2] = x30*(x13 + x23 + x25 + x29 + x37*(-n3*x17 + x3) + x40 + x42 + 1.0*log(n3*x3)) + x32*(648226800.0*n1 + 648226800.0*n2 + 2453522400.0*n4 + 1207224000.0*n5) + x35 + (*endmember[2].mu0)(T, P);
result[3] = x30*(1.0*x11 + x13 + x19*(x18 + 0.5*x3) + x27 + x37*(-n4*x17 + x3) + x39 + x43 + 1.0*log(n4*x3) - 0.69314718055994495) + x32*(1061424000.0*n1 + 1061424000.0*n2 + 2453522400.0*n3 + 408240000.0*n5) - 0.0001333638164913885*x34 + (*endmember[3].mu0)(T, P);
result[4] = x30*(x21 + x37*(-n5*x17 + x3) + x38 + x42 + x43 + 1.0*log(n5*x3)) + x32*(637729200.0*n1 + 637729200.0*n2 + 1207224000.0*n3 + 408240000.0*n4) + x35 + (*endmember[4].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 1.0*n5;
    double x1 = 1.0*n2;
    double x2 = 1.0*n1;
    double x3 = 1.0*n3;
    double x4 = x2 + x3;
    double x5 = x0 + x1 + x4;
    double x6 = 0.5*n4 + x5;
    double x7 = 1.0/x6;
    double x8 = n1 + n3;
    double x9 = n4 + n5;
    double x10 = n2 + x8 + x9;
    double x11 = 1.0/x10;
    double x12 = 1.0*x11;
    double x13 = pow(x10, -2);
    double x14 = -x13*x6;
    double x15 = x12 + x14;
    double x16 = x10*x15;
    double x17 = x16*x7;
    double x18 = -2.0*x13;
    double x19 = pow(x10, -3);
    double x20 = 2*x19;
    double x21 = x20*x6;
    double x22 = 2.0*n2;
    double x23 = 2.0*n3;
    double x24 = 1.0*n4;
    double x25 = 2.0*n5;
    double x26 = 2.0*n1 + x22 + x23 + x24 + x25;
    double x27 = x26*x7;
    double x28 = x10*x27;
    double x29 = x26/((x6)*(x6));
    double x30 = x16*x29;
    double x31 = 4.0*x17 + x28*(x18 + x21) - 1.0*x30;
    double x32 = n1 + n2 + n4;
    double x33 = 1.0/x32;
    double x34 = x11 - x13*x32;
    double x35 = x33*x34;
    double x36 = 2.0*x10;
    double x37 = x35*x36;
    double x38 = -2*x13;
    double x39 = x1 + x2 + x24;
    double x40 = x10*x39;
    double x41 = x33*x40;
    double x42 = x41*(x20*x32 + x38);
    double x43 = -x34*x40/((x32)*(x32));
    double x44 = x1*x13;
    double x45 = x13*x3;
    double x46 = x0*x13;
    double x47 = x35*x39;
    double x48 = x44 + x45 + x46 + x47;
    double x49 = x37 + x42 + x43 + x48;
    double x50 = x31 + x49;
    double x51 = x13*x24;
    double x52 = x15*x27;
    double x53 = -x24;
    double x54 = -x13*(-x0 + x53) + x51 + x52;
    double x55 = 1.0/x8;
    double x56 = x11 - x13*x8;
    double x57 = x55*x56;
    double x58 = x4*x57;
    double x59 = x10*x4;
    double x60 = x55*x59;
    double x61 = x36*x57 - x56*x59/((x8)*(x8)) + x58 + x60*(x20*x8 + x38);
    double x62 = x54 + x61;
    double x63 = 8.3144626181532395*T;
    double x64 = 3.5*n4 + x5;
    double x65 = pow(x64, -2);
    double x66 = 1.0*x65;
    double x67 = 324113400.0*n3 + 530712000.0*n4 + 318864600.0*n5;
    double x68 = (n1*x67 + n2*x67 + 729.0*n3*(444600.0*n1 + 444600.0*n2 + 1682800.0*n4 + 828000.0*n5) + 224.0*n4*(2369250.0*n1 + 2369250.0*n2 + 5476612.5*n3 + 911250.0*n5) + 729.0*n5*(437400.0*n1 + 437400.0*n2 + 828000.0*n3 + 280000.0*n4))/((x64)*(x64)*(x64));
    double x69 = 7.6207895137936279e-5*x68;
    double x70 = 648226800.0*n3 + 1061424000.0*n4 + 637729200.0*n5;
    double x71 = 3.8103947568968139e-5*x65;
    double x72 = x69 - x70*x71;
    double x73 = -x66*(24699.999999999996*n3 + 40444.444444444438*n4 + 24299.999999999996*n5) + x72;
    double x74 = -x13;
    double x75 = -n1;
    double x76 = x60*(-x20*(-n3 + x75) + x74);
    double x77 = x58 + x76;
    double x78 = -2.0*x11;
    double x79 = x54 + x78;
    double x80 = x41*(-x20*(-n2 - n4 + x75) + x74);
    double x81 = x48 + x80;
    double x82 = x31 + x81;
    double x83 = 648226800.0*n1 + 648226800.0*n2 + 2453522400.0*n4 + 1207224000.0*n5;
    double x84 = -x71*x83;
    double x85 = 1.0/x64;
    double x86 = x72 + x84 + 24699.999999999996*x85;
    double x87 = -3.0*x11;
    double x88 = x0 + x24;
    double x89 = x13*x88 + x51 + x52;
    double x90 = x77 + x89;
    double x91 = 0.5*x11 + x14;
    double x92 = x36*x7*x91;
    double x93 = 1.0*x17 + x28*(-1.5*x13 + x21) + x92;
    double x94 = -0.5*x30 + x93;
    double x95 = 0.0001333638164913885*x65;
    double x96 = 1061424000.0*n1 + 1061424000.0*n2 + 2453522400.0*n3 + 408240000.0*n5;
    double x97 = 0.000266727632982777*x68 - x71*x96;
    double x98 = -x70*x95 + 40444.444444444438*x85 + x97;
    double x99 = -4.0*x11 + x58 + x76 + x89;
    double x100 = 637729200.0*n1 + 637729200.0*n2 + 1207224000.0*n3 + 408240000.0*n4;
    double x101 = -x100*x71;
    double x102 = x101 + x72 + 24299.999999999996*x85;
    double x103 = -x44;
    double x104 = x19*x22;
    double x105 = 1.0*x10;
    double x106 = x37 + x42 + x43;
    double x107 = x31 + x46;
    double x108 = -x2;
    double x109 = -x13*(x108 - x3);
    double x110 = x109 + x12;
    double x111 = x110 + x45;
    double x112 = x31 + x46 + x87;
    double x113 = -1.0*x13;
    double x114 = x10*(x104 + x113) + x103 + x45 + x47;
    double x115 = x112 + x114 + x80;
    double x116 = x109 + x89;
    double x117 = x46 + x94;
    double x118 = -x45;
    double x119 = x19*x23;
    double x120 = -x13*(-x1 + x108 + x53) + x44;
    double x121 = x69 + x84;
    double x122 = x10*(x113 + x119) + x118;
    double x123 = 2.0*n4*x19;
    double x124 = x29*x91;
    double x125 = 1.0/x9;
    double x126 = x11 - x13*x9;
    double x127 = x125*x126;
    double x128 = x10*x88;
    double x129 = x125*x128*(x20*x9 + x38) - x126*x128/((x9)*(x9)) + x127*x36 + x127*x88;
    double x130 = x129 + x27*x91 - x51;

result[0] = x63*(x50 + x62) + x73;
result[1] = x63*(x50 + x77 + x79) + x73;
result[2] = x63*(x61 + x79 + x82) + x86;
result[3] = x63*(x49 + x87 + x90 + x94) + x98;
result[4] = x102 + x63*(x82 + x99);
result[5] = x63*(x10*(x104 + x18) + x103 + x106 + x107 + x111 + x47 + x54 + x105*(-n2*x13 + x11)/n2) + x73;
result[6] = x63*(x115 + x13*x4 + x54) + x86;
result[7] = x63*(x106 + x114 + x116 + x117 + x78) + x98;
result[8] = x102 + x63*(x115 + x116);
result[9] = x121 + x63*(x10*(x119 + x18) + x107 + x118 + x12 + x120 + x62 + x105*(-n3*x13 + x11)/n3) - x66*(24699.999999999996*n1 + 24699.999999999996*n2 + 93488.888888888876*n4 + 45999.999999999993*n5);
result[10] = x63*(x117 + x122 + x13*x39 + x44 + x99) - x83*x95 + 93488.888888888876*x85 + x97;
result[11] = x101 + x121 + x63*(x112 + x120 + x122 + x90) + 45999.999999999993*x85;
result[12] = x63*(-0.5*x10*x124 + x10*(x123 + x18) + x110 + x130 + x28*(x113 + x21) + x49 + x92 + x105*(-n4*x13 + x11)/n4) - 3.5*x65*(40444.444444444438*n1 + 40444.444444444438*n2 + 93488.888888888876*n3 + 15555.555555555553*n5) + 0.00093354671543971956*x68 - x95*x96;
result[13] = -x100*x95 + x63*(x10*(x113 + x123) - x105*x124 + x109 + x130 + x78 + x81 + x93) + 15555.555555555553*x85 + x97;
result[14] = x101 + x63*(x10*(x18 + x19*x25) + x111 + x120 + x129 + x31 - x46 + x51 + x52 + x105*(-n5*x13 + x11)/n5) - x66*(24299.999999999996*n1 + 24299.999999999996*n2 + 45999.999999999993*n3 + 15555.555555555553*n4) + x69;
}
        
static void coder_d3gdn3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = n1 + n3;
    double x1 = n4 + n5;
    double x2 = n2 + x0 + x1;
    double x3 = pow(x2, -3);
    double x4 = 2.0*n2;
    double x5 = -x3*x4;
    double x6 = 1.0/x0;
    double x7 = 1.0/x2;
    double x8 = pow(x2, -2);
    double x9 = -x0*x8 + x7;
    double x10 = x6*x9;
    double x11 = 1.0*n1;
    double x12 = 1.0*n3;
    double x13 = x11 + x12;
    double x14 = pow(x0, -2);
    double x15 = x14*x9;
    double x16 = x13*x15;
    double x17 = -2*x8;
    double x18 = 2*x3;
    double x19 = x0*x18 + x17;
    double x20 = x13*x6;
    double x21 = x19*x20;
    double x22 = 3.0*x2;
    double x23 = 6*x3;
    double x24 = pow(x2, -4);
    double x25 = 6*x24;
    double x26 = x2*x20;
    double x27 = x13*x14;
    double x28 = 2*x2;
    double x29 = 3.0*x10 - x15*x22 - 2*x16 + x19*x22*x6 - x19*x27*x28 + 2*x21 + x26*(-x0*x25 + x23) + x13*x28*x9/((x0)*(x0)*(x0));
    double x30 = x29 + x5;
    double x31 = 1.0*n4;
    double x32 = 1.0*n5;
    double x33 = x31 + x32;
    double x34 = -x18*x33;
    double x35 = 2.0*n3;
    double x36 = -x3*x35;
    double x37 = 2.0*x3;
    double x38 = -n4*x37;
    double x39 = 2.0*n5;
    double x40 = -x3*x39;
    double x41 = x36 + x38 + x40;
    double x42 = x34 + x41;
    double x43 = 1.0*n2;
    double x44 = x13 + x32 + x43;
    double x45 = 0.5*n4 + x44;
    double x46 = 1.0/x45;
    double x47 = -x45*x8;
    double x48 = x47 + 1.0*x7;
    double x49 = x46*x48;
    double x50 = 2.0*x8;
    double x51 = -x50;
    double x52 = x18*x45;
    double x53 = x51 + x52;
    double x54 = 2.0*n1 + x31 + x35 + x39 + x4;
    double x55 = x46*x54;
    double x56 = x53*x55;
    double x57 = 6.0*x2;
    double x58 = x46*x53;
    double x59 = pow(x45, -2);
    double x60 = x48*x59;
    double x61 = 2.0*x54;
    double x62 = 6.0*x3;
    double x63 = -x25*x45;
    double x64 = x2*x55;
    double x65 = x2*x61;
    double x66 = pow(x45, -3);
    double x67 = x48*x66;
    double x68 = x53*x59;
    double x69 = 6.0*x49 + 2*x56 + x57*x58 - x57*x60 - x60*x61 + x64*(x62 + x63) + x65*x67 - x65*x68;
    double x70 = x42 + x69;
    double x71 = n1 + n2 + n4;
    double x72 = 1.0/x71;
    double x73 = x7 - x71*x8;
    double x74 = x72*x73;
    double x75 = x11 + x31 + x43;
    double x76 = pow(x71, -2);
    double x77 = x73*x76;
    double x78 = x75*x77;
    double x79 = x17 + x18*x71;
    double x80 = x72*x75;
    double x81 = x79*x80;
    double x82 = x2*x80;
    double x83 = x75*x76;
    double x84 = x22*x72*x79 - x22*x77 - x28*x79*x83 + x28*x73*x75/((x71)*(x71)*(x71)) + 3.0*x74 - 2*x78 + 2*x81 + x82*(x23 - x25*x71);
    double x85 = x70 + x84;
    double x86 = 8.3144626181532395*T;
    double x87 = 24699.999999999996*n3;
    double x88 = 40444.444444444438*n4;
    double x89 = 24299.999999999996*n5;
    double x90 = 3.5*n4 + x44;
    double x91 = pow(x90, -3);
    double x92 = 4.0*x91;
    double x93 = 648226800.0*n3 + 1061424000.0*n4 + 637729200.0*n5;
    double x94 = 7.6207895137936279e-5*x91;
    double x95 = 324113400.0*n3 + 530712000.0*n4 + 318864600.0*n5;
    double x96 = (n1*x95 + n2*x95 + 729.0*n3*(444600.0*n1 + 444600.0*n2 + 1682800.0*n4 + 828000.0*n5) + 224.0*n4*(2369250.0*n1 + 2369250.0*n2 + 5476612.5*n3 + 911250.0*n5) + 729.0*n5*(437400.0*n1 + 437400.0*n2 + 828000.0*n3 + 280000.0*n4))/((x90)*(x90)*(x90)*(x90));
    double x97 = -0.00022862368541380884*x96;
    double x98 = x93*x94 + x97;
    double x99 = -x92*(-x87 - x88 - x89) + x98;
    double x100 = 1.0*x8;
    double x101 = -x8;
    double x102 = -n1;
    double x103 = -n3 + x102;
    double x104 = x101 - x103*x18;
    double x105 = x104*x2;
    double x106 = x105*x6;
    double x107 = 4*x3;
    double x108 = 2*n1;
    double x109 = 2*n3;
    double x110 = 3*x24;
    double x111 = -x110*(x108 + x109);
    double x112 = x104*x20;
    double x113 = -x105*x27 + x112 - x16 + x21;
    double x114 = 2.0*x10 + 2.0*x106 + x113 + x26*(x107 + x111);
    double x115 = x114 + x5;
    double x116 = -n2 - n4 + x102;
    double x117 = x101 - x116*x18;
    double x118 = x117*x2;
    double x119 = x118*x72;
    double x120 = 2*n2;
    double x121 = 2*n4;
    double x122 = -x110*(x108 + x120 + x121);
    double x123 = x117*x80;
    double x124 = -x118*x83 + x123 - x78 + x81;
    double x125 = 2.0*x119 + x124 + 2.0*x74 + x82*(x107 + x122);
    double x126 = x91*(x87 + x88 + x89);
    double x127 = 2.0*x126 + x98;
    double x128 = 648226800.0*n1 + 648226800.0*n2 + 2453522400.0*n4 + 1207224000.0*n5;
    double x129 = x128*x94;
    double x130 = pow(x90, -2);
    double x131 = x129 - 49399.999999999993*x130;
    double x132 = x127 + x131;
    double x133 = -x31;
    double x134 = x18*(x133 - x32);
    double x135 = x52 - 1.5*x8;
    double x136 = 1.0*x2;
    double x137 = x136*x54;
    double x138 = x137*x59;
    double x139 = -x135*x138;
    double x140 = x135*x55;
    double x141 = x54*x60;
    double x142 = 0.5*x2*x54;
    double x143 = x2*x46;
    double x144 = x135*x143;
    double x145 = x136*x58 + 4.0*x144 + x64*(5.0*x3 + x63);
    double x146 = x137*x67 + x140 - 1.5*x141 - x142*x68 + x145 + x56;
    double x147 = x134 + x139 + x146 - x22*x60 + 5.0*x49;
    double x148 = x147 + x41;
    double x149 = x148 + x84;
    double x150 = x115 + x50;
    double x151 = -0.0008001828989483309*x96;
    double x152 = 1061424000.0*n1 + 1061424000.0*n2 + 2453522400.0*n3 + 408240000.0*n5;
    double x153 = x152*x94;
    double x154 = 0.000266727632982777*x91;
    double x155 = x153 + x154*x93;
    double x156 = x151 + x155;
    double x157 = 7.0*x126 - 80888.888888888876*x130 + x156;
    double x158 = x134 + x69;
    double x159 = x158 + x41;
    double x160 = x125 + x159;
    double x161 = 637729200.0*n1 + 637729200.0*n2 + 1207224000.0*n3 + 408240000.0*n4;
    double x162 = x161*x94;
    double x163 = -48599.999999999993*x130 + x162;
    double x164 = x127 + x163;
    double x165 = 2*x112 + x26*(x111 + x18);
    double x166 = x165 + x5;
    double x167 = 3.0*x8;
    double x168 = x167 + x70;
    double x169 = 4.0*x8;
    double x170 = x166 + x169;
    double x171 = 2*x123;
    double x172 = x82*(x122 + x18);
    double x173 = x171 + x172 + x5;
    double x174 = x173 + x70;
    double x175 = 24699.999999999996*n1;
    double x176 = 24699.999999999996*n2;
    double x177 = 93488.888888888876*n4;
    double x178 = 45999.999999999993*n5;
    double x179 = -x175 - x176 - x177 - x178;
    double x180 = 2.0*x91;
    double x181 = x131 - x179*x180 + x98;
    double x182 = x114 + x169;
    double x183 = 1.0*x119 + x124 + 1.0*x74 + x82*(x107 + x116*x25);
    double x184 = x128*x154;
    double x185 = -220383.33333333331*x130 + x156 + x184;
    double x186 = x129 - 94999.999999999985*x130 + x162 + x98;
    double x187 = x166 + 5.0*x8;
    double x188 = x47 + 0.5*x7;
    double x189 = x188*x59;
    double x190 = -x136*x60 + 2*x140;
    double x191 = -x100;
    double x192 = x191 + x52;
    double x193 = 4.0*x3;
    double x194 = 2.0*x143*x192 + 2.0*x144 + x64*(x193 + x63);
    double x195 = -x136*x189 - 1.0*x141 + x142*x67 + x190 + x194 + 2.0*x49;
    double x196 = x188*x46;
    double x197 = x139 + 2.0*x196;
    double x198 = x197 + x42;
    double x199 = x195 + x198 + x84;
    double x200 = 0.00093354671543971956*x91;
    double x201 = 40444.444444444438*n1;
    double x202 = 40444.444444444438*n2;
    double x203 = 93488.888888888876*n3;
    double x204 = 15555.555555555553*n5;
    double x205 = 7.0*x91;
    double x206 = x152*x154 - 0.0028006401463191587*x96;
    double x207 = -x205*(-x201 - x202 - x203 - x204) + x206;
    double x208 = -283111.11111111112*x130 + x200*x93 + x207;
    double x209 = 2.0*x2;
    double x210 = x146 - x189*x209 - x209*x60 + 3.0*x49;
    double x211 = x187 + x210;
    double x212 = x125 + x198;
    double x213 = x154*x161 - 0.00080018289894833101*x96;
    double x214 = -141050.0*x130 + x155 + x213;
    double x215 = x165 + 6.0*x8;
    double x216 = -24299.999999999996*n1 - 24299.999999999996*n2 - 45999.999999999993*n3 - 15555.555555555553*n4;
    double x217 = -x180*x216;
    double x218 = x163 + x217 + x98;
    double x219 = -n2*x8 + x7;
    double x220 = 1.0/n2;
    double x221 = x219*x220;
    double x222 = 6.0*x24;
    double x223 = -n2*x222;
    double x224 = x120*x3;
    double x225 = x2*x220;
    double x226 = -x13*x18;
    double x227 = n2*x193;
    double x228 = x227 + x70;
    double x229 = x226 + x228;
    double x230 = -x169;
    double x231 = x230 + x84;
    double x232 = -x11;
    double x233 = x191 + x2*(x193 + x223) + x221 + x225*(x101 + x224);
    double x234 = x125 + x233;
    double x235 = x226 + x227;
    double x236 = x148 + x235;
    double x237 = x159 + x235;
    double x238 = x2*(x223 + x37);
    double x239 = x169 + x238;
    double x240 = x171 + x172 + x239;
    double x241 = x86*(x229 + x240);
    double x242 = x167 + x235 + x238;
    double x243 = -n3*x222;
    double x244 = x109*x3;
    double x245 = 1.0/n3;
    double x246 = x2*x245;
    double x247 = -n3*x8 + x7;
    double x248 = n3*x193 + x38 + x40;
    double x249 = x245*x247 + x248;
    double x250 = -x18*x75;
    double x251 = x250 + x69;
    double x252 = x230 + x251;
    double x253 = x129 + x97;
    double x254 = x175 + x176 + x177 + x178;
    double x255 = x114 + x2*(x193 + x243) + x246*(x101 + x244) + x249;
    double x256 = x153 + x184;
    double x257 = x250 + x5;
    double x258 = -91999.999999999985*x130 + x162 + x253;
    double x259 = x2*(x243 + x37) + x248 + x34;
    double x260 = x197 + x259;
    double x261 = -327211.11111111107*n3;
    double x262 = -n4*x222;
    double x263 = x192*x55;
    double x264 = x121*x3;
    double x265 = 1.0/n4;
    double x266 = x2*x265;
    double x267 = -n4*x8 + x7;
    double x268 = x189*x54;
    double x269 = x189*x2;
    double x270 = x188*x66;
    double x271 = n4*x193 + x40;
    double x272 = -x1*x8 + x7;
    double x273 = 1.0/x1;
    double x274 = 3.0*x273;
    double x275 = 2*x33;
    double x276 = pow(x1, -2);
    double x277 = x272*x276;
    double x278 = x1*x18 + x17;
    double x279 = x275*x278;
    double x280 = x2*x273*x33*(-x1*x25 + x23) + x2*x274*x278 - x2*x276*x279 - x22*x277 + x226 + x272*x274 + x273*x279 - x275*x277 + x36 + x2*x272*x275/((x1)*(x1)*(x1));
    double x281 = x280 + x5;
    double x282 = -x138*x192 + x265*x267 + x281;
    double x283 = x135*x59;
    double x284 = 4.0*x196 + x271;
    double x285 = -n5*x8 + x7;
    double x286 = 1.0/n5;

result[0] = x86*(x30 + x85) + x99;
result[1] = x86*(x100 + x115 + x85) + x99;
result[2] = x132 + x86*(x100 + x125 + x30 + x70);
result[3] = x157 + x86*(x149 + x150);
result[4] = x164 + x86*(x150 + x160);
result[5] = x86*(x166 + x168 + x84) + x99;
result[6] = x132 + x86*(1.0*x10 + 1.0*x106 + x113 + x125 + x168 + x26*(x103*x25 + x107) + x5);
result[7] = x157 + x86*(x149 + x170);
result[8] = x164 + x86*(x160 + x170);
result[9] = x181 + x86*(x167 + x174 + x29);
result[10] = x185 + x86*(x148 + x182 + x183 + x5);
result[11] = x186 + x86*(x159 + x173 + x182);
result[12] = x208 + x86*(x187 + x199);
result[13] = x214 + x86*(x211 + x212);
result[14] = x218 + x86*(x174 + x215);
result[15] = x86*(x2*(x223 + x62) + x221 + x225*(x17 + x224) + x229 + x231 - x136*x219/((n2)*(n2))) + x99;
result[16] = x132 + x86*(x18*(-x12 + x232) + x228 + x234);
result[17] = x157 + x86*(x233 + x236 + x84);
result[18] = x164 + x86*(x234 + x237);
result[19] = x181 + x241;
result[20] = x185 + x86*(x183 + x236 + x239);
result[21] = x186 + x86*(x237 + x240);
result[22] = x208 + x86*(x199 + x242);
result[23] = x214 + x86*(x210 + x212 + x242);
result[24] = x218 + x241;
result[25] = -x179*x92 + x253 + x86*(x2*(x243 + x62) + x246*(x17 + x244) + x249 + x252 + x30 + x34 - x136*x247/((n3)*(n3)));
result[26] = -186977.77777777775*x130 + x151 + x205*x254 + x256 + x86*(x147 + x18*(x133 + x232 - x43) + x255 + x5);
result[27] = x180*x254 + x258 + x86*(x158 + x191 + x255 + x257);
result[28] = x128*x200 - 654422.22222222213*x130 + x207 + x86*(x195 + x215 + x257 + x260);
result[29] = -270044.44444444444*x130 + x213 + x256 + x86*(x211 + x250 + x260);
result[30] = x217 + x258 + x86*(x170 + x251 + x259);
result[31] = x152*x200 - x205*(-141555.55555555556*n1 - 141555.55555555556*n2 - 54444.444444444445*n5 + x261) - x205*(-141555.55555555553*n1 - 141555.55555555553*n2 - 54444.444444444438*n5 + x261) + x86*(x142*x270 + x192*x22*x46 + 3.0*x196 + x2*(x262 + x62) + x231 + 2*x263 + x266*(x17 + x264) - 1.0*x268 - 1.5*x269 + x271 + x282 + x64*(3.0*x3 + x63) - x136*x267/((n4)*(n4))) - 0.0098022405121170556*x96;
result[32] = -108888.88888888888*x130 + x161*x200 + x205*(x201 + x202 + x203 + x204) + x206 + x86*(x125 + x137*x270 + x140 - x142*x283 - x189*x22 + x194 + x2*(x193 + x262) + x263 + x266*(x101 + x264) - 1.5*x268 + x282 + x284 + x51);
result[33] = -31111.111111111106*x130 + x153 - x180*(-85050.0*n1 - 85050.0*n2 - 161000.0*n3 - 54444.444444444445*n4) + x213 + x86*(x145 + x173 - x189*x61 + x190 + x2*(x262 + x37) - 4.0*x269 + x270*x65 + x280 - x283*x65 + x284 + 1.0*x49 + x50);
result[34] = x162 - x216*x92 + x86*(n5*x193 + x2*x286*(n5*x18 + x17) + x2*(-n5*x222 + x62) + x252 + x281 + x285*x286 + x38 - x136*x285/((n5)*(n5))) + x97;
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
    double x3 = 1.0*n1;
    double x4 = 1.0*n3 + x3;
    double x5 = 1.0*n4;
    double x6 = 1.0*n5;
    double x7 = 1.0*n2;
    double x8 = sqrt(-2.6762700000000001*T + 13.381350000000001);
    double x9 = 0.27336948519328247*x8;
    double x10 = fmin(4, x9);
    double x11 = (-x9 + 4 >= 0. ? 1. : 0.)/x8;

result = n1*(*endmember[0].dmu0dT)(T, P) + 8.3144626181532395*n2*log(n2*x2) + n2*(*endmember[1].dmu0dT)(T, P) + 8.3144626181532395*n3*log(n3*x2) + n3*(*endmember[2].dmu0dT)(T, P) + 8.3144626181532395*n4*(log(n4*x2) - 0.69314718055994495) + n4*(*endmember[3].dmu0dT)(T, P) + 8.3144626181532395*n5*log(n5*x2) + n5*(*endmember[4].dmu0dT)(T, P) + 8.3144626181532395*x4*log(x0*x2) + 8.3144626181532395*(x5 + x6)*log(x1*x2) + 8.3144626181532395*(x3 + x5 + x7)*log(x2*(n1 + n2 + n4)) + 8.3144626181532395*(2.0*n1 + 2.0*n2 + 2.0*n3 + 2.0*n5 + x5)*log(x2*(0.5*n4 + x4 + x6 + x7)) + ((T >= 5.0) ? (
   -13.381349999999999*n2
)
: (
   (1.0/3.0)*n2*(-73.424526463911391*((x10)*(x10))*x11 + 40.14405*x10 - 0.36580527606911306*x11*(40.14405*T - 200.72024999999999) - 40.14405)
));
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
    double x18 = n1 + n2 + n4;
    double x19 = 8.3144626181532395*n5;
    double x20 = -x19*x3;
    double x21 = 8.3144626181532395*n1;
    double x22 = x11 + x21 + x7;
    double x23 = x20 + 8.3144626181532395*log(x18*x3) + x2*x22*(-x13*x18 + x3)/x18;
    double x24 = x11 + x19;
    double x25 = -x24*x3;
    double x26 = x21 + x9;
    double x27 = x25 + 8.3144626181532395*log(x0*x3) + x2*x26*(-x0*x13 + x3)/x0;
    double x28 = -x26*x3;
    double x29 = 8.3144626181532395*x2;
    double x30 = sqrt(-2.6762700000000001*T + 13.381350000000001);
    double x31 = 0.27336948519328247*x30;
    double x32 = fmin(4, x31);
    double x33 = (-x31 + 4 >= 0. ? 1. : 0.)/x30;
    double x34 = x12 + x16 + x6;
    double x35 = x10 + x23;
    double x36 = -x22*x3;
    double x37 = x28 + 8.3144626181532395*log(x1*x3) + x2*x24*(-x1*x13 + x3)/x1;

result[0] = x17 + x23 + x27 + (*endmember[0].dmu0dT)(T, P);
result[1] = x25 + x28 + x29*(-n2*x13 + x3) + x34 + x35 + ((T >= 5.0) ? (
   -13.381349999999999
)
: (
   -24.47484215463713*((x32)*(x32))*x33 + 13.381349999999999*x32 - 0.12193509202303768*x33*(40.14405*T - 200.72024999999999) - 13.381349999999999
)) + 8.3144626181532395*log(n2*x3) + (*endmember[1].dmu0dT)(T, P);
result[2] = x20 + x27 + x29*(-n3*x13 + x3) + x34 + x36 + x8 + 8.3144626181532395*log(n3*x3) + (*endmember[2].dmu0dT)(T, P);
result[3] = x15*(x14 + 0.5*x3) + x29*(-n4*x13 + x3) + x35 + x37 + 8.3144626181532395*x5 + x8 + 8.3144626181532395*log(n4*x3) + (*endmember[3].dmu0dT)(T, P) - 5.7631463216439762;
result[4] = x17 + x29*(-n5*x13 + x3) + x36 + x37 + 8.3144626181532395*log(n5*x3) + (*endmember[4].dmu0dT)(T, P);
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
    
    double x0 = 2.6762700000000001*T;
    double x1 = -x0 + 13.381350000000001;
    double x2 = sqrt(x1);
    double x3 = 0.27336948519328247*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 40.14405*T - 200.72024999999999;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 13.381350000000001);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result = n1*(*endmember[0].d2mu0dT2)(T, P) + n2*(*endmember[1].d2mu0dT2)(T, P) + n3*(*endmember[2].d2mu0dT2)(T, P) + n4*(*endmember[3].d2mu0dT2)(T, P) + n5*(*endmember[4].d2mu0dT2)(T, P) + ((T >= 5.0) ? (
   0
)
: (
   -1.0/3.0*n2*(98.251928719786079*x10*x6 - 26.859079173375005*x10*x8 + 53.718158346750009*((x4)*(x4))*x7*x9 + 0.48949684309274261*x5*x6 - 0.13381350000000003*x5*x8 + 29.369810585564558*x4/x2)
));
    return result;
}
        
static void coder_d3gdndt2(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 2.6762700000000001*T;
    double x1 = -x0 + 13.381350000000001;
    double x2 = sqrt(x1);
    double x3 = 0.27336948519328247*x2;
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 40.14405*T - 200.72024999999999;
    double x6 = x4/pow(x1, 3.0/2.0);
    double x7 = 1.0/(x0 - 13.381350000000001);
    double x8 = x7*0;
    double x9 = fmin(4, x3);
    double x10 = ((x9)*(x9));

result[0] = (*endmember[0].d2mu0dT2)(T, P);
result[1] = ((T >= 5.0) ? (
   0
)
: (
   -32.75064290659536*x10*x6 + 8.9530263911250003*x10*x8 - 17.906052782250001*((x4)*(x4))*x7*x9 - 0.16316561436424754*x5*x6 + 0.044604500000000005*x5*x8 - 9.7899368618548515*x4/x2
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
    
    double x0 = 2.6762700000000001*T;
    double x1 = -x0 + 13.381350000000001;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.27336948519328247*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 58.951157231871647*x2*x4;
    double x6 = x0 - 13.381350000000001;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 40.14405*T - 200.72024999999999;
    double x10 = x4/pow(x1, 5.0/2.0);
    double x11 = pow(x6, -2);
    double x12 = x11*x8;
    double x13 = x2*0;
    double x14 = fmin(4, x3);
    double x15 = ((x14)*(x14));

result = n1*(*endmember[0].d3mu0dT3)(T, P) + n2*(*endmember[1].d3mu0dT3)(T, P) + n3*(*endmember[2].d3mu0dT3)(T, P) + n4*(*endmember[3].d3mu0dT3)(T, P) + n5*(*endmember[4].d3mu0dT3)(T, P) + ((T >= 5.0) ? (
   0
)
: (
   -1.0/3.0*n2*(394.42303391235282*x10*x15 + 1.9650385743957215*x10*x9 - 215.64644345798496*x11*x14*((x4)*(x4)) + 107.82322172899248*x12*x15 + 0.53718158346750011*x12*x9 - 9.8251928719786079*x13*x15 - 0.048949684309274273*x13*x9 - x14*x5*x8 + 19.650385743957216*x2*((x4)*(x4)*(x4)) + x5 - 16.115447504025003*x8/x6)
));
    return result;
}
        
static void coder_d4gdndt3(double T, double P, double n[5], double result[5]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double n4 = n[3];
    double n5 = n[4];

    double x0 = 2.6762700000000001*T;
    double x1 = -x0 + 13.381350000000001;
    double x2 = pow(x1, -3.0/2.0);
    double x3 = 0.27336948519328247*sqrt(x1);
    double x4 = (-x3 + 4 >= 0. ? 1. : 0.);
    double x5 = 19.650385743957216*x2*x4;
    double x6 = x0 - 13.381350000000001;
    double x7 = x3 - 4;
    double x8 = 0;
    double x9 = 40.14405*T - 200.72024999999999;
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
   -131.47434463745094*x10*x15 - 0.65501285813190713*x10*x9 + 71.88214781932831*x11*x14*((x4)*(x4)) - 35.941073909664155*x12*x15 - 0.17906052782250004*x12*x9 + 3.275064290659536*x13*x15 + 0.016316561436424758*x13*x9 + x14*x5*x8 - 6.5501285813190719*x2*((x4)*(x4)*(x4)) - x5 + 5.3718158346750009*x8/x6
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

