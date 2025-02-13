
#include <stdlib.h>
#include <stdio.h>

#include "Kyanite_stx21_em_coder.h"


typedef struct _endmembers {
  const char *(*name) (void);
  const char *(*formula) (void);
  const double (*mw) (void);
  const double *(*elements) (void);
  double (*mu0) (double t, double p);
  double (*dmu0dT) (double t, double p);
  double (*dmu0dP) (double t, double p);
  double (*d2mu0dT2) (double t, double p);
  double (*d2mu0dTdP) (double t, double p);
  double (*d2mu0dP2) (double t, double p);
  double (*d3mu0dT3) (double t, double p);
  double (*d3mu0dT2dP) (double t, double p);
  double (*d3mu0dTdP2) (double t, double p);
  double (*d3mu0dP3) (double t, double p);
} Endmembers;

static Endmembers endmember[] = {
  {
     Kyanite_stx21_em_coder_calib_name,
     Kyanite_stx21_em_coder_calib_formula,
     Kyanite_stx21_em_coder_calib_mw,
     Kyanite_stx21_em_coder_calib_elements,
     Kyanite_stx21_em_coder_calib_g,
     Kyanite_stx21_em_coder_calib_dgdt,
     Kyanite_stx21_em_coder_calib_dgdp,
     Kyanite_stx21_em_coder_calib_d2gdt2,
     Kyanite_stx21_em_coder_calib_d2gdtdp,
     Kyanite_stx21_em_coder_calib_d2gdp2,
     Kyanite_stx21_em_coder_calib_d3gdt3,
     Kyanite_stx21_em_coder_calib_d3gdt2dp,
     Kyanite_stx21_em_coder_calib_d3gdtdp2,
     Kyanite_stx21_em_coder_calib_d3gdp3,

  },

};
static int nc = (sizeof endmember / sizeof(struct _endmembers));

static const double R=8.3143;


static char *identifier = "Kyanite_stx21_ph.phml:22fdc2be0fa6fba472dbebcd2232bddacc83527d:Sat Mar  9 00:09:23 2024";


#include "Kyanite_stx21_ph_coder_codegen_calc.h"
#include "Kyanite_stx21_ph_coder_codegen_calib.h"

const char *Kyanite_stx21_ph_coder_calib_identifier(void) {
    return identifier;
}

const char *Kyanite_stx21_ph_coder_calib_name(void) {
    return "Kyanite_stx21_ph";
}

char *Kyanite_stx21_ph_coder_calib_formula(double T, double P, double n[1]) {
    double sum, elm[1];
    const double *end0 = (*endmember[0].elements)();
    int i;
    const char *fmt = "Al%5.3fSiO5";
    char *result = (char *) malloc(12*sizeof(char));
    for (i=0, sum=0.0; i<nc; i++) sum += n[i];
    if (sum == 0.0) return result;
    elm[0] = end0[13]*n[0]/sum;
    sprintf(result, fmt, elm[0]);
    return result;

}

double *Kyanite_stx21_ph_coder_calib_conv_elm_to_moles(double *e) {
    double *n = (double *) malloc(1*sizeof(double));
    n[0]=e[13]/2.0;
    return n;

}

int Kyanite_stx21_ph_coder_calib_test_moles(double *n) {
    int result = 1;
    result &= (n[0] > 0.0);
    return result;

}

const char *Kyanite_stx21_ph_coder_calib_endmember_name(int index) {
    return (*endmember[index].name)();
}

const char *Kyanite_stx21_ph_coder_calib_endmember_formula(int index) {
    return (*endmember[index].formula)();
}

const double Kyanite_stx21_ph_coder_calib_endmember_mw(int index) {
    return (*endmember[index].mw)();
}

const double *Kyanite_stx21_ph_coder_calib_endmember_elements(int index) {
    return (*endmember[index].elements)();
}

double Kyanite_stx21_ph_coder_calib_endmember_mu0(int index, double t, double p) {
    return (*endmember[index].mu0)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_dmu0dT(int index, double t, double p) {
    return (*endmember[index].dmu0dT)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_dmu0dP(int index, double t, double p) {
    return (*endmember[index].dmu0dP)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d2mu0dT2(int index, double t, double p) {
    return (*endmember[index].d2mu0dT2)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d2mu0dTdP(int index, double t, double p) {
    return (*endmember[index].d2mu0dTdP)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d2mu0dP2(int index, double t, double p) {
    return (*endmember[index].d2mu0dP2)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d3mu0dT3(int index, double t, double p) {
    return (*endmember[index].d3mu0dT3)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d3mu0dT2dP(int index, double t, double p) {
    return (*endmember[index].d3mu0dT2dP)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d3mu0dTdP2(int index, double t, double p) {
    return (*endmember[index].d3mu0dTdP2)(t, p);
}

double Kyanite_stx21_ph_coder_calib_endmember_d3mu0dP3(int index, double t, double p) {
    return (*endmember[index].d3mu0dP3)(t, p);
}

double Kyanite_stx21_ph_coder_calib_g(double T, double P, double n[1]) {
    return coder_g(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_dgdn(double T, double P, double n[1], double result[1]) {
    coder_dgdn(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d2gdn2(double T, double P, double n[1], double result[1]) {
    coder_d2gdn2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d3gdn3(double T, double P, double n[1], double result[1]) {
    coder_d3gdn3(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_dgdt(double T, double P, double n[1]) {
    return coder_dgdt(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d2gdndt(double T, double P, double n[1], double result[1]) {
    coder_d2gdndt(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d3gdn2dt(double T, double P, double n[1], double result[1]) {
    coder_d3gdn2dt(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d4gdn3dt(double T, double P, double n[1], double result[1]) {
    coder_d4gdn3dt(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_dgdp(double T, double P, double n[1]) {
    return coder_dgdp(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d2gdndp(double T, double P, double n[1], double result[1]) {
    coder_d2gdndp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d3gdn2dp(double T, double P, double n[1], double result[1]) {
    coder_d3gdn2dp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d4gdn3dp(double T, double P, double n[1], double result[1]) {
    coder_d4gdn3dp(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d2gdt2(double T, double P, double n[1]) {
    return coder_d2gdt2(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d3gdndt2(double T, double P, double n[1], double result[1]) {
    coder_d3gdndt2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d4gdn2dt2(double T, double P, double n[1], double result[1]) {
    coder_d4gdn2dt2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn3dt2(double T, double P, double n[1], double result[1]) {
    coder_d5gdn3dt2(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d2gdtdp(double T, double P, double n[1]) {
    return coder_d2gdtdp(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d3gdndtdp(double T, double P, double n[1], double result[1]) {
    coder_d3gdndtdp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d4gdn2dtdp(double T, double P, double n[1], double result[1]) {
    coder_d4gdn2dtdp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn3dtdp(double T, double P, double n[1], double result[1]) {
    coder_d5gdn3dtdp(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d2gdp2(double T, double P, double n[1]) {
    return coder_d2gdp2(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d3gdndp2(double T, double P, double n[1], double result[1]) {
    coder_d3gdndp2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d4gdn2dp2(double T, double P, double n[1], double result[1]) {
    coder_d4gdn2dp2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn3dp2(double T, double P, double n[1], double result[1]) {
    coder_d5gdn3dp2(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d3gdt3(double T, double P, double n[1]) {
    return coder_d3gdt3(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d4gdndt3(double T, double P, double n[1], double result[1]) {
    coder_d4gdndt3(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn2dt3(double T, double P, double n[1], double result[1]) {
    coder_d5gdn2dt3(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d6gdn3dt3(double T, double P, double n[1], double result[1]) {
    coder_d6gdn3dt3(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d3gdt2dp(double T, double P, double n[1]) {
    return coder_d3gdt2dp(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d4gdndt2dp(double T, double P, double n[1], double result[1]) {
    coder_d4gdndt2dp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn2dt2dp(double T, double P, double n[1], double result[1]) {
    coder_d5gdn2dt2dp(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d6gdn3dt2dp(double T, double P, double n[1], double result[1]) {
    coder_d6gdn3dt2dp(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d3gdtdp2(double T, double P, double n[1]) {
    return coder_d3gdtdp2(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d4gdndtdp2(double T, double P, double n[1], double result[1]) {
    coder_d4gdndtdp2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn2dtdp2(double T, double P, double n[1], double result[1]) {
    coder_d5gdn2dtdp2(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d6gdn3dtdp2(double T, double P, double n[1], double result[1]) {
    coder_d6gdn3dtdp2(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_d3gdp3(double T, double P, double n[1]) {
    return coder_d3gdp3(T, P, n);
}

void Kyanite_stx21_ph_coder_calib_d4gdndp3(double T, double P, double n[1], double result[1]) {
    coder_d4gdndp3(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d5gdn2dp3(double T, double P, double n[1], double result[1]) {
    coder_d5gdn2dp3(T, P, n, result);
}

void Kyanite_stx21_ph_coder_calib_d6gdn3dp3(double T, double P, double n[1], double result[1]) {
    coder_d6gdn3dp3(T, P, n, result);
}

double Kyanite_stx21_ph_coder_calib_s(double T, double P, double n[1]) {
    return coder_s(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_v(double T, double P, double n[1]) {
    return coder_v(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_cv(double T, double P, double n[1]) {
    return coder_cv(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_cp(double T, double P, double n[1]) {
    return coder_cp(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_dcpdt(double T, double P, double n[1]) {
    return coder_dcpdt(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_alpha(double T, double P, double n[1]) {
    return coder_alpha(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_beta(double T, double P, double n[1]) {
    return coder_beta(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_K(double T, double P, double n[1]) {
    return coder_K(T, P, n);
}

double Kyanite_stx21_ph_coder_calib_Kp(double T, double P, double n[1]) {
    return coder_Kp(T, P, n);
}

int Kyanite_stx21_ph_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Kyanite_stx21_ph_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Kyanite_stx21_ph_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Kyanite_stx21_ph_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Kyanite_stx21_ph_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Kyanite_stx21_ph_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Kyanite_stx21_ph_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Kyanite_stx21_ph_coder_dparam_g(double T, double P, double n[1], int index) {
    return coder_dparam_g(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_dgdt(double T, double P, double n[1], int index) {
    return coder_dparam_dgdt(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_dgdp(double T, double P, double n[1], int index) {
    return coder_dparam_dgdp(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d2gdt2(double T, double P, double n[1], int index) {
    return coder_dparam_d2gdt2(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d2gdtdp(double T, double P, double n[1], int index) {
    return coder_dparam_d2gdtdp(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d2gdp2(double T, double P, double n[1], int index) {
    return coder_dparam_d2gdp2(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d3gdt3(double T, double P, double n[1], int index) {
    return coder_dparam_d3gdt3(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d3gdt2dp(double T, double P, double n[1], int index) {
    return coder_dparam_d3gdt2dp(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d3gdtdp2(double T, double P, double n[1], int index) {
    return coder_dparam_d3gdtdp2(T, P, n, index);
}

double Kyanite_stx21_ph_coder_dparam_d3gdp3(double T, double P, double n[1], int index) {
    return coder_dparam_d3gdp3(T, P, n, index);
}

void Kyanite_stx21_ph_coder_dparam_dgdn(double T, double P, double n[1], int index, double result[1]) {
    coder_dparam_dgdn(T, P, n, index, result);
}

