#include <math.h>


static double coder_g(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = 1.0/(n1 + n2 + n3);
    double x1 = 1.0*n1;
    double x2 = 1.0*n2;
    double x3 = log(n3*x0);
    double x4 = 1.0*n3;

result = 8.3144626181532395*T*(x1*log(n1*x0) + x2*log(n2*x0) + x3*x4 + x4*(x3 - 0.69314718055994495) + (2.0*n1 + 2.0*n2 + x4)*log(x0*(0.5*n3 + x1 + x2))) + n1*(*endmember[0].mu0)(T, P) + n2*(*endmember[1].mu0)(T, P) + n3*(*endmember[2].mu0)(T, P);
    return result;
}
        
static void coder_dgdn(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = 1.0/x0;
    double x2 = n1*x1;
    double x3 = pow(x0, -2);
    double x4 = 1.0*x0;
    double x5 = 1.0*n2;
    double x6 = -x1*x5;
    double x7 = 1.0*n1 + 0.5*n3 + x5;
    double x8 = log(x1*x7);
    double x9 = n3*x1;
    double x10 = -x3*x7;
    double x11 = x0*(2.0*n1 + 2.0*n2 + 1.0*n3)/x7;
    double x12 = x11*(1.0*x1 + x10) + 2.0*x8 - 2.0*x9;
    double x13 = 8.3144626181532395*T;
    double x14 = -1.0*x2;

result[0] = x13*(x12 + x4*(-n1*x3 + x1) + x6 + 1.0*log(x2)) + (*endmember[0].mu0)(T, P);
result[1] = x13*(x12 + x14 + x4*(-n2*x3 + x1) + 1.0*log(n2*x1)) + (*endmember[1].mu0)(T, P);
result[2] = x13*(2.0*x0*(-n3*x3 + x1) + x11*(0.5*x1 + x10) + x14 + x6 + 1.0*x8 + 2.0*log(x9) - 0.69314718055994495) + (*endmember[2].mu0)(T, P);
}
        
static void coder_d2gdn2(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = pow(x0, -2);
    double x2 = 1.0*n2;
    double x3 = x1*x2;
    double x4 = 1.0*n1;
    double x5 = x1*x4;
    double x6 = -x5;
    double x7 = 2.0*x1;
    double x8 = -x7;
    double x9 = pow(x0, -3);
    double x10 = 2.0*n1;
    double x11 = x10*x9;
    double x12 = 1.0/x0;
    double x13 = 1.0*x0;
    double x14 = 1.0*x12;
    double x15 = n3*x7;
    double x16 = 0.5*n3 + x2 + x4;
    double x17 = -x1*x16;
    double x18 = x14 + x17;
    double x19 = 1.0/x16;
    double x20 = 2.0*n2;
    double x21 = 1.0*n3 + x10 + x20;
    double x22 = x19*x21;
    double x23 = x15 + x18*x22;
    double x24 = x18*x19;
    double x25 = 2*x16*x9;
    double x26 = x0*x22;
    double x27 = x21/((x16)*(x16));
    double x28 = x18*x27;
    double x29 = 4.0*x0*x24 - x13*x28 + x26*(x25 + x8);
    double x30 = x14 + x23 + x29;
    double x31 = 8.3144626181532395*T;
    double x32 = -1.0*x1;
    double x33 = x0*(x11 + x32) + x23 + x3 + x6;
    double x34 = 2.0*x12;
    double x35 = 0.5*x12 + x17;
    double x36 = 2.0*x0;
    double x37 = x19*x35*x36;
    double x38 = 0.5*x0;
    double x39 = x13*x24 + x26*(-1.5*x1 + x25) - x28*x38 - x34 + x37;
    double x40 = x20*x9;
    double x41 = -x3 + x5;

result[0] = x31*(x0*(x11 + x8) + x3 + x30 + x6 + x13*(-n1*x1 + x12)/n1);
result[1] = x31*(-x14 + x29 + x33);
result[2] = x31*(x33 + x39);
result[3] = x31*(x0*(x40 + x8) + x30 + x41 + x13*(-n2*x1 + x12)/n2);
result[4] = x31*(x0*(x32 + x40) + x23 + x39 + x41);
result[5] = x31*(x0*(4.0*n3*x9 - 4.0*x1) - x15 + x22*x35 + x26*(x25 + x32) - x27*x35*x38 + x3 + x34 + x37 + x5 + x36*(-n3*x1 + x12)/n3);
}
        
static void coder_d3gdn3(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = pow(x0, -2);
    double x2 = -4.0*x1;
    double x3 = pow(x0, -3);
    double x4 = 6.0*x3;
    double x5 = pow(x0, -4);
    double x6 = 6.0*x5;
    double x7 = -n1*x6;
    double x8 = -2*x1;
    double x9 = 2*x3;
    double x10 = n1*x9;
    double x11 = 1.0/n1;
    double x12 = x0*x11;
    double x13 = 1.0/x0;
    double x14 = -n1*x1 + x13;
    double x15 = 1.0*x0;
    double x16 = x11*x14;
    double x17 = 4.0*x3;
    double x18 = 2.0*n2;
    double x19 = -x18*x3;
    double x20 = -n3*x17;
    double x21 = n1*x17 + x19 + x20;
    double x22 = 1.0*n1 + 1.0*n2 + 0.5*n3;
    double x23 = 1.0/x22;
    double x24 = -x1*x22;
    double x25 = 1.0*x13 + x24;
    double x26 = x23*x25;
    double x27 = 2.0*x1;
    double x28 = -x27;
    double x29 = x22*x9;
    double x30 = x28 + x29;
    double x31 = 2.0*n1;
    double x32 = 1.0*n3 + x18 + x31;
    double x33 = x23*x32;
    double x34 = x30*x33;
    double x35 = 6.0*x0;
    double x36 = x23*x30;
    double x37 = pow(x22, -2);
    double x38 = x25*x37;
    double x39 = x32*x38;
    double x40 = -6*x22*x5;
    double x41 = x0*x33;
    double x42 = 2.0*x0;
    double x43 = x32*x42;
    double x44 = pow(x22, -3);
    double x45 = x25*x44;
    double x46 = x30*x37;
    double x47 = 6.0*x26 + 2*x34 + x35*x36 - x35*x38 - 2.0*x39 + x41*(x4 + x40) + x43*x45 - x43*x46;
    double x48 = x16 + x21 + x47;
    double x49 = 8.3144626181532395*T;
    double x50 = -x1;
    double x51 = x0*(x17 + x7) + x12*(x10 + x50);
    double x52 = 1.0*x1;
    double x53 = -x52;
    double x54 = -1.5*x1 + x29;
    double x55 = x15*x32;
    double x56 = x37*x55;
    double x57 = -x54*x56;
    double x58 = x33*x54;
    double x59 = x0*x23;
    double x60 = x54*x59;
    double x61 = 0.5*x0*x32;
    double x62 = -3.0*x0*x38 + x15*x36 + 5.0*x26 + x34 - 1.5*x39 + x41*(5.0*x3 + x40) + x45*x55 - x46*x61 + x58 + 4.0*x60;
    double x63 = x53 + x57 + x62;
    double x64 = 2.0*x3;
    double x65 = x0*(x64 + x7) + x21;
    double x66 = x57 + x65;
    double x67 = 0.5*x13 + x24;
    double x68 = x23*x67;
    double x69 = x29 + x53;
    double x70 = x59*x69;
    double x71 = x37*x67;
    double x72 = 3.0*x1 - x15*x38 - x15*x71 + 2.0*x26 - 1.0*x39 + x41*(x17 + x40) + x45*x61 + 2*x58 + 2.0*x60 + 2.0*x68 + 2.0*x70;
    double x73 = -n2*x6;
    double x74 = n2*x9;
    double x75 = 1.0/n2;
    double x76 = x0*x75;
    double x77 = -n2*x1 + x13;
    double x78 = -x3*x31;
    double x79 = n2*x17 + x20 + x78;
    double x80 = x75*x77 + x79;
    double x81 = -n3*x1 + x13;
    double x82 = 2.0/n3;

result[0] = x49*(x0*(x4 + x7) + x12*(x10 + x8) + x2 + x48 - x14*x15/((n1)*(n1)));
result[1] = x49*(x28 + x48 + x51);
result[2] = x49*(x16 + x21 + x51 + x63);
result[3] = x49*(x47 + x52 + x65);
result[4] = x49*(x27 + x62 + x66);
result[5] = x49*(x66 + x72);
result[6] = x49*(x0*(x4 + x73) + x2 + x47 + x76*(x74 + x8) + x80 - x15*x77/((n2)*(n2)));
result[7] = x49*(x0*(x17 + x73) + x63 + x76*(x50 + x74) + x80);
result[8] = x49*(x0*(x64 + x73) + x57 + x72 + x79);
result[9] = x49*(8.0*n3*x3 - 1.5*x0*x71 + x0*x82*(n3*x9 + x8) + x0*(-12.0*n3*x5 + 12.0*x3) - 8.0*x1 + x19 - 1.0*x32*x71 + 2*x33*x69 + x41*(3.0*x3 + x40) + x44*x61*x67 - x56*x69 + 3.0*x68 + 3.0*x70 + x78 + x81*x82 - x42*x81/((n3)*(n3)));
}
        
static double coder_dgdt(double T, double P, double n[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];
    double result;
    
    double x0 = 1.0/(n1 + n2 + n3);
    double x1 = log(n3*x0);
    double x2 = 8.3144626181532395*n3;

result = 8.3144626181532395*n1*log(n1*x0) + n1*(*endmember[0].dmu0dT)(T, P) + 8.3144626181532395*n2*log(n2*x0) + n2*(*endmember[1].dmu0dT)(T, P) + n3*(*endmember[2].dmu0dT)(T, P) + x1*x2 + x2*(x1 - 0.69314718055994495) + 8.3144626181532395*(2.0*n1 + 2.0*n2 + 1.0*n3)*log(x0*(1.0*n1 + 1.0*n2 + 0.5*n3));
    return result;
}
        
static void coder_d2gdndt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = 1.0/x0;
    double x2 = n1*x1;
    double x3 = pow(x0, -2);
    double x4 = 8.3144626181532395*x0;
    double x5 = n2*x1;
    double x6 = -8.3144626181532395*x5;
    double x7 = 1.0*n1 + 1.0*n2 + 0.5*n3;
    double x8 = log(x1*x7);
    double x9 = n3*x1;
    double x10 = -x3*x7;
    double x11 = x0*(16.628925236306479*n1 + 16.628925236306479*n2 + 8.3144626181532395*n3)/x7;
    double x12 = x11*(1.0*x1 + x10) + 16.628925236306479*x8 - 16.628925236306479*x9;
    double x13 = -8.3144626181532395*x2;

result[0] = x12 + x4*(-n1*x3 + x1) + x6 + 8.3144626181532395*log(x2) + (*endmember[0].dmu0dT)(T, P);
result[1] = x12 + x13 + x4*(-n2*x3 + x1) + 8.3144626181532395*log(x5) + (*endmember[1].dmu0dT)(T, P);
result[2] = 16.628925236306479*x0*(-n3*x3 + x1) + x11*(0.5*x1 + x10) + x13 + x6 + 8.3144626181532395*x8 + 16.628925236306479*log(x9) + (*endmember[2].dmu0dT)(T, P) - 5.7631463216439762;
}
        
static void coder_d3gdn2dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = pow(x0, -2);
    double x2 = 8.3144626181532395*x1;
    double x3 = n2*x2;
    double x4 = n1*x1;
    double x5 = 8.3144626181532395*x4;
    double x6 = -x5;
    double x7 = 16.628925236306479*x1;
    double x8 = -x7;
    double x9 = pow(x0, -3);
    double x10 = 16.628925236306479*n1;
    double x11 = x10*x9;
    double x12 = 1.0/x0;
    double x13 = 8.3144626181532395*x0;
    double x14 = 8.3144626181532395*x12;
    double x15 = n3*x7;
    double x16 = 1.0*n1 + 1.0*n2 + 0.5*n3;
    double x17 = -x1*x16;
    double x18 = 1.0*x12 + x17;
    double x19 = 1.0/x16;
    double x20 = 16.628925236306479*n2;
    double x21 = 8.3144626181532395*n3 + x10 + x20;
    double x22 = x19*x21;
    double x23 = x15 + x18*x22;
    double x24 = x18*x19;
    double x25 = 2*x16*x9;
    double x26 = x0*x22;
    double x27 = x0*x21/((x16)*(x16));
    double x28 = x18*x27;
    double x29 = 33.257850472612958*x0*x24 + x26*(-2.0*x1 + x25) - 1.0*x28;
    double x30 = x14 + x23 + x29;
    double x31 = -x2;
    double x32 = x0*(x11 + x31) + x23 + x3 + x6;
    double x33 = 16.628925236306479*x12;
    double x34 = 0.5*x12 + x17;
    double x35 = 16.628925236306479*x0;
    double x36 = x19*x34*x35;
    double x37 = x13*x24 + x26*(-1.5*x1 + x25) - 0.5*x28 - x33 + x36;
    double x38 = x20*x9;
    double x39 = -x3 + x5;

result[0] = x0*(x11 + x8) + x3 + x30 + x6 + x13*(x12 - x4)/n1;
result[1] = -x14 + x29 + x32;
result[2] = x32 + x37;
result[3] = x0*(x38 + x8) + x30 + x39 + x13*(-n2*x1 + x12)/n2;
result[4] = x0*(x31 + x38) + x23 + x37 + x39;
result[5] = x0*(33.257850472612958*n3*x9 - 33.257850472612958*x1) - x15 + x22*x34 + x26*(-1.0*x1 + x25) - 0.5*x27*x34 + x3 + x33 + x36 + x5 + x35*(-n3*x1 + x12)/n3;
}
        
static void coder_d4gdn3dt(double T, double P, double n[3], double result[3]) {
    double n1 = n[0];
    double n2 = n[1];
    double n3 = n[2];

    double x0 = n1 + n2 + n3;
    double x1 = pow(x0, -2);
    double x2 = -33.257850472612958*x1;
    double x3 = pow(x0, -3);
    double x4 = 49.886775708919437*x3;
    double x5 = pow(x0, -4);
    double x6 = 49.886775708919437*x5;
    double x7 = -n1*x6;
    double x8 = -2*x1;
    double x9 = 2*x3;
    double x10 = n1*x9;
    double x11 = 8.3144626181532395/n1;
    double x12 = x0*x11;
    double x13 = 1.0/x0;
    double x14 = -n1*x1 + x13;
    double x15 = 8.3144626181532395*x0;
    double x16 = x11*x14;
    double x17 = 33.257850472612958*x3;
    double x18 = 16.628925236306479*n2;
    double x19 = -x18*x3;
    double x20 = -n3*x17;
    double x21 = n1*x17 + x19 + x20;
    double x22 = 1.0*n1 + 1.0*n2 + 0.5*n3;
    double x23 = 1.0/x22;
    double x24 = -x1*x22;
    double x25 = 1.0*x13 + x24;
    double x26 = x23*x25;
    double x27 = x22*x9;
    double x28 = -2.0*x1 + x27;
    double x29 = 16.628925236306479*n1;
    double x30 = 8.3144626181532395*n3 + x18 + x29;
    double x31 = x23*x30;
    double x32 = x28*x31;
    double x33 = 49.886775708919437*x0;
    double x34 = x23*x28;
    double x35 = pow(x22, -2);
    double x36 = x25*x35;
    double x37 = x30*x36;
    double x38 = -6*x22*x5;
    double x39 = x0*x31;
    double x40 = x0*x30;
    double x41 = 2.0*x40;
    double x42 = pow(x22, -3);
    double x43 = x25*x42;
    double x44 = x28*x35;
    double x45 = 49.886775708919437*x26 + 2*x32 + x33*x34 - x33*x36 - 2.0*x37 + x39*(6.0*x3 + x38) + x41*x43 - x41*x44;
    double x46 = x16 + x21 + x45;
    double x47 = 16.628925236306479*x1;
    double x48 = -x1;
    double x49 = x0*(x17 + x7) + x12*(x10 + x48);
    double x50 = -1.5*x1 + x27;
    double x51 = 1.0*x30;
    double x52 = x0*x51;
    double x53 = x35*x52;
    double x54 = -x50*x53;
    double x55 = 8.3144626181532395*x1;
    double x56 = x31*x50;
    double x57 = x0*x23;
    double x58 = x50*x57;
    double x59 = 0.5*x40;
    double x60 = -24.943387854459719*x0*x36 + x15*x34 + 41.572313090766201*x26 + x32 - 1.5*x37 + x39*(5.0*x3 + x38) + x43*x52 - x44*x59 + x56 + 33.257850472612958*x58;
    double x61 = x54 - x55 + x60;
    double x62 = 16.628925236306479*x3;
    double x63 = x0*(x62 + x7) + x21;
    double x64 = x54 + x63;
    double x65 = 0.5*x13 + x24;
    double x66 = x23*x65;
    double x67 = -1.0*x1 + x27;
    double x68 = x57*x67;
    double x69 = x35*x65;
    double x70 = 24.943387854459719*x1 - x15*x36 - x15*x69 + 16.628925236306479*x26 - 1.0*x37 + x39*(4.0*x3 + x38) + x43*x59 + 2*x56 + 16.628925236306479*x58 + 16.628925236306479*x66 + 16.628925236306479*x68;
    double x71 = -n2*x6;
    double x72 = n2*x9;
    double x73 = 8.3144626181532395/n2;
    double x74 = x0*x73;
    double x75 = -n2*x1 + x13;
    double x76 = -x29*x3;
    double x77 = n2*x17 + x20 + x76;
    double x78 = x73*x75 + x77;
    double x79 = -n3*x1 + x13;
    double x80 = 16.628925236306479/n3;

result[0] = x0*(x4 + x7) + x12*(x10 + x8) + x2 + x46 - x14*x15/((n1)*(n1));
result[1] = x46 - x47 + x49;
result[2] = x16 + x21 + x49 + x61;
result[3] = x45 + x55 + x63;
result[4] = x47 + x60 + x64;
result[5] = x64 + x70;
result[6] = x0*(x4 + x71) + x2 + x45 + x74*(x72 + x8) + x78 - x15*x75/((n2)*(n2));
result[7] = x0*(x17 + x71) + x61 + x74*(x48 + x72) + x78;
result[8] = x0*(x62 + x71) + x54 + x70 + x77;
result[9] = 66.515700945225916*n3*x3 - 12.471693927229859*x0*x69 + x0*x80*(n3*x9 + x8) + x0*(-99.773551417838874*n3*x5 + 99.773551417838874*x3) - 66.515700945225916*x1 + x19 + 2*x31*x67 + x39*(3.0*x3 + x38) + x42*x59*x65 - x51*x69 - x53*x67 + 24.943387854459719*x66 + 24.943387854459719*x68 + x76 + x79*x80 - 16.628925236306479*x0*x79/((n3)*(n3));
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

