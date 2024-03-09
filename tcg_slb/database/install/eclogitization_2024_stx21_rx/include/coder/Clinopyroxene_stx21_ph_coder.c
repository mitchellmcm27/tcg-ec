
#include <stdlib.h>
#include <stdio.h>

#include "Diopside_stx21_em_coder.h"
#include "Hedenbergite_stx21_em_coder.h"
#include "Clinoenstatite_stx21_em_coder.h"
#include "CaTschermaks_stx21_em_coder.h"
#include "Jadeite_stx21_em_coder.h"


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
     Diopside_stx21_em_coder_calib_name,
     Diopside_stx21_em_coder_calib_formula,
     Diopside_stx21_em_coder_calib_mw,
     Diopside_stx21_em_coder_calib_elements,
     Diopside_stx21_em_coder_calib_g,
     Diopside_stx21_em_coder_calib_dgdt,
     Diopside_stx21_em_coder_calib_dgdp,
     Diopside_stx21_em_coder_calib_d2gdt2,
     Diopside_stx21_em_coder_calib_d2gdtdp,
     Diopside_stx21_em_coder_calib_d2gdp2,
     Diopside_stx21_em_coder_calib_d3gdt3,
     Diopside_stx21_em_coder_calib_d3gdt2dp,
     Diopside_stx21_em_coder_calib_d3gdtdp2,
     Diopside_stx21_em_coder_calib_d3gdp3,

  },
  {
     Hedenbergite_stx21_em_coder_calib_name,
     Hedenbergite_stx21_em_coder_calib_formula,
     Hedenbergite_stx21_em_coder_calib_mw,
     Hedenbergite_stx21_em_coder_calib_elements,
     Hedenbergite_stx21_em_coder_calib_g,
     Hedenbergite_stx21_em_coder_calib_dgdt,
     Hedenbergite_stx21_em_coder_calib_dgdp,
     Hedenbergite_stx21_em_coder_calib_d2gdt2,
     Hedenbergite_stx21_em_coder_calib_d2gdtdp,
     Hedenbergite_stx21_em_coder_calib_d2gdp2,
     Hedenbergite_stx21_em_coder_calib_d3gdt3,
     Hedenbergite_stx21_em_coder_calib_d3gdt2dp,
     Hedenbergite_stx21_em_coder_calib_d3gdtdp2,
     Hedenbergite_stx21_em_coder_calib_d3gdp3,

  },
  {
     Clinoenstatite_stx21_em_coder_calib_name,
     Clinoenstatite_stx21_em_coder_calib_formula,
     Clinoenstatite_stx21_em_coder_calib_mw,
     Clinoenstatite_stx21_em_coder_calib_elements,
     Clinoenstatite_stx21_em_coder_calib_g,
     Clinoenstatite_stx21_em_coder_calib_dgdt,
     Clinoenstatite_stx21_em_coder_calib_dgdp,
     Clinoenstatite_stx21_em_coder_calib_d2gdt2,
     Clinoenstatite_stx21_em_coder_calib_d2gdtdp,
     Clinoenstatite_stx21_em_coder_calib_d2gdp2,
     Clinoenstatite_stx21_em_coder_calib_d3gdt3,
     Clinoenstatite_stx21_em_coder_calib_d3gdt2dp,
     Clinoenstatite_stx21_em_coder_calib_d3gdtdp2,
     Clinoenstatite_stx21_em_coder_calib_d3gdp3,

  },
  {
     CaTschermaks_stx21_em_coder_calib_name,
     CaTschermaks_stx21_em_coder_calib_formula,
     CaTschermaks_stx21_em_coder_calib_mw,
     CaTschermaks_stx21_em_coder_calib_elements,
     CaTschermaks_stx21_em_coder_calib_g,
     CaTschermaks_stx21_em_coder_calib_dgdt,
     CaTschermaks_stx21_em_coder_calib_dgdp,
     CaTschermaks_stx21_em_coder_calib_d2gdt2,
     CaTschermaks_stx21_em_coder_calib_d2gdtdp,
     CaTschermaks_stx21_em_coder_calib_d2gdp2,
     CaTschermaks_stx21_em_coder_calib_d3gdt3,
     CaTschermaks_stx21_em_coder_calib_d3gdt2dp,
     CaTschermaks_stx21_em_coder_calib_d3gdtdp2,
     CaTschermaks_stx21_em_coder_calib_d3gdp3,

  },
  {
     Jadeite_stx21_em_coder_calib_name,
     Jadeite_stx21_em_coder_calib_formula,
     Jadeite_stx21_em_coder_calib_mw,
     Jadeite_stx21_em_coder_calib_elements,
     Jadeite_stx21_em_coder_calib_g,
     Jadeite_stx21_em_coder_calib_dgdt,
     Jadeite_stx21_em_coder_calib_dgdp,
     Jadeite_stx21_em_coder_calib_d2gdt2,
     Jadeite_stx21_em_coder_calib_d2gdtdp,
     Jadeite_stx21_em_coder_calib_d2gdp2,
     Jadeite_stx21_em_coder_calib_d3gdt3,
     Jadeite_stx21_em_coder_calib_d3gdt2dp,
     Jadeite_stx21_em_coder_calib_d3gdtdp2,
     Jadeite_stx21_em_coder_calib_d3gdp3,

  },

};
static int nc = (sizeof endmember / sizeof(struct _endmembers));

static const double R=8.3143;


static char *identifier = "Clinopyroxene_stx21_ph.phml:d8e6d4fe999fb591b6b78e1a9d5df0b5264c1cef:Sat Mar  9 00:09:08 2024";


#include "Clinopyroxene_stx21_ph_coder_codegen_calc.h"
#include "Clinopyroxene_stx21_ph_coder_codegen_calib.h"

const char *Clinopyroxene_stx21_ph_coder_calib_identifier(void) {
    return identifier;
}

const char *Clinopyroxene_stx21_ph_coder_calib_name(void) {
    return "Clinopyroxene_stx21_ph";
}

char *Clinopyroxene_stx21_ph_coder_calib_formula(double T, double P, double n[5]) {
    double sum, elm[6];
    const double *end0 = (*endmember[0].elements)();
    const double *end1 = (*endmember[1].elements)();
    const double *end2 = (*endmember[2].elements)();
    const double *end3 = (*endmember[3].elements)();
    const double *end4 = (*endmember[4].elements)();
    int i;
    const char *fmt = "Na%5.3fCa%5.3fMg%5.3fFe%5.3fAl%5.3fSi%5.3fO6";
    char *result = (char *) malloc(45*sizeof(char));
    for (i=0, sum=0.0; i<nc; i++) sum += n[i];
    if (sum == 0.0) return result;
    elm[0] = end0[11]*n[0]/sum + end1[11]*n[1]/sum + end2[11]*n[2]/sum + end3[11]*n[3]/sum + end4[11]*n[4]/sum;
    elm[1] = end0[20]*n[0]/sum + end1[20]*n[1]/sum + end2[20]*n[2]/sum + end3[20]*n[3]/sum + end4[20]*n[4]/sum;
    elm[2] = end0[12]*n[0]/sum + end1[12]*n[1]/sum + end2[12]*n[2]/sum + end3[12]*n[3]/sum + end4[12]*n[4]/sum;
    elm[3] = end0[26]*n[0]/sum + end1[26]*n[1]/sum + end2[26]*n[2]/sum + end3[26]*n[3]/sum + end4[26]*n[4]/sum;
    elm[4] = end0[13]*n[0]/sum + end1[13]*n[1]/sum + end2[13]*n[2]/sum + end3[13]*n[3]/sum + end4[13]*n[4]/sum;
    elm[5] = end0[14]*n[0]/sum + end1[14]*n[1]/sum + end2[14]*n[2]/sum + end3[14]*n[3]/sum + end4[14]*n[4]/sum;
    sprintf(result, fmt, elm[0], elm[1], elm[2], elm[3], elm[4], elm[5]);
    return result;

}

double *Clinopyroxene_stx21_ph_coder_calib_conv_elm_to_moles(double *e) {
    double *n = (double *) malloc(5*sizeof(double));
    n[0]=e[20] - e[26] - (e[13]-e[11])/2.0;
    n[1]=e[26];
    n[2]=(e[12] - (e[20] - e[26] - (e[13]-e[11])/2.0))/2.0;
    n[3]=(e[13]-e[11])/2.0;
    n[4]=e[11];
    return n;

}

int Clinopyroxene_stx21_ph_coder_calib_test_moles(double *n) {
    int result = 1;
    result &= (n[0] > 0.0);
    result &= (n[1] > 0.0);
    result &= (n[2] > 0.0);
    result &= (n[3] > 0.0);
    result &= (n[4] > 0.0);
    return result;

}

const char *Clinopyroxene_stx21_ph_coder_calib_endmember_name(int index) {
    return (*endmember[index].name)();
}

const char *Clinopyroxene_stx21_ph_coder_calib_endmember_formula(int index) {
    return (*endmember[index].formula)();
}

const double Clinopyroxene_stx21_ph_coder_calib_endmember_mw(int index) {
    return (*endmember[index].mw)();
}

const double *Clinopyroxene_stx21_ph_coder_calib_endmember_elements(int index) {
    return (*endmember[index].elements)();
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_mu0(int index, double t, double p) {
    return (*endmember[index].mu0)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_dmu0dT(int index, double t, double p) {
    return (*endmember[index].dmu0dT)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_dmu0dP(int index, double t, double p) {
    return (*endmember[index].dmu0dP)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d2mu0dT2(int index, double t, double p) {
    return (*endmember[index].d2mu0dT2)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d2mu0dTdP(int index, double t, double p) {
    return (*endmember[index].d2mu0dTdP)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d2mu0dP2(int index, double t, double p) {
    return (*endmember[index].d2mu0dP2)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d3mu0dT3(int index, double t, double p) {
    return (*endmember[index].d3mu0dT3)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d3mu0dT2dP(int index, double t, double p) {
    return (*endmember[index].d3mu0dT2dP)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d3mu0dTdP2(int index, double t, double p) {
    return (*endmember[index].d3mu0dTdP2)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_endmember_d3mu0dP3(int index, double t, double p) {
    return (*endmember[index].d3mu0dP3)(t, p);
}

double Clinopyroxene_stx21_ph_coder_calib_g(double T, double P, double n[5]) {
    return coder_g(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_dgdn(double T, double P, double n[5], double result[5]) {
    coder_dgdn(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d2gdn2(double T, double P, double n[5], double result[5]) {
    coder_d2gdn2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdn3(double T, double P, double n[5], double result[5]) {
    coder_d3gdn3(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_dgdt(double T, double P, double n[5]) {
    return coder_dgdt(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d2gdndt(double T, double P, double n[5], double result[5]) {
    coder_d2gdndt(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdn2dt(double T, double P, double n[5], double result[5]) {
    coder_d3gdn2dt(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdn3dt(double T, double P, double n[5], double result[5]) {
    coder_d4gdn3dt(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_dgdp(double T, double P, double n[5]) {
    return coder_dgdp(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d2gdndp(double T, double P, double n[5], double result[5]) {
    coder_d2gdndp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdn2dp(double T, double P, double n[5], double result[5]) {
    coder_d3gdn2dp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdn3dp(double T, double P, double n[5], double result[5]) {
    coder_d4gdn3dp(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d2gdt2(double T, double P, double n[5]) {
    return coder_d2gdt2(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdndt2(double T, double P, double n[5], double result[5]) {
    coder_d3gdndt2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdn2dt2(double T, double P, double n[5], double result[5]) {
    coder_d4gdn2dt2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn3dt2(double T, double P, double n[5], double result[5]) {
    coder_d5gdn3dt2(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d2gdtdp(double T, double P, double n[5]) {
    return coder_d2gdtdp(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdndtdp(double T, double P, double n[5], double result[5]) {
    coder_d3gdndtdp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdn2dtdp(double T, double P, double n[5], double result[5]) {
    coder_d4gdn2dtdp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn3dtdp(double T, double P, double n[5], double result[5]) {
    coder_d5gdn3dtdp(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d2gdp2(double T, double P, double n[5]) {
    return coder_d2gdp2(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d3gdndp2(double T, double P, double n[5], double result[5]) {
    coder_d3gdndp2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdn2dp2(double T, double P, double n[5], double result[5]) {
    coder_d4gdn2dp2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn3dp2(double T, double P, double n[5], double result[5]) {
    coder_d5gdn3dp2(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d3gdt3(double T, double P, double n[5]) {
    return coder_d3gdt3(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdndt3(double T, double P, double n[5], double result[5]) {
    coder_d4gdndt3(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn2dt3(double T, double P, double n[5], double result[5]) {
    coder_d5gdn2dt3(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d6gdn3dt3(double T, double P, double n[5], double result[5]) {
    coder_d6gdn3dt3(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d3gdt2dp(double T, double P, double n[5]) {
    return coder_d3gdt2dp(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdndt2dp(double T, double P, double n[5], double result[5]) {
    coder_d4gdndt2dp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn2dt2dp(double T, double P, double n[5], double result[5]) {
    coder_d5gdn2dt2dp(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d6gdn3dt2dp(double T, double P, double n[5], double result[5]) {
    coder_d6gdn3dt2dp(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d3gdtdp2(double T, double P, double n[5]) {
    return coder_d3gdtdp2(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdndtdp2(double T, double P, double n[5], double result[5]) {
    coder_d4gdndtdp2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn2dtdp2(double T, double P, double n[5], double result[5]) {
    coder_d5gdn2dtdp2(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d6gdn3dtdp2(double T, double P, double n[5], double result[5]) {
    coder_d6gdn3dtdp2(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_d3gdp3(double T, double P, double n[5]) {
    return coder_d3gdp3(T, P, n);
}

void Clinopyroxene_stx21_ph_coder_calib_d4gdndp3(double T, double P, double n[5], double result[5]) {
    coder_d4gdndp3(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d5gdn2dp3(double T, double P, double n[5], double result[5]) {
    coder_d5gdn2dp3(T, P, n, result);
}

void Clinopyroxene_stx21_ph_coder_calib_d6gdn3dp3(double T, double P, double n[5], double result[5]) {
    coder_d6gdn3dp3(T, P, n, result);
}

double Clinopyroxene_stx21_ph_coder_calib_s(double T, double P, double n[5]) {
    return coder_s(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_v(double T, double P, double n[5]) {
    return coder_v(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_cv(double T, double P, double n[5]) {
    return coder_cv(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_cp(double T, double P, double n[5]) {
    return coder_cp(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_dcpdt(double T, double P, double n[5]) {
    return coder_dcpdt(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_alpha(double T, double P, double n[5]) {
    return coder_alpha(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_beta(double T, double P, double n[5]) {
    return coder_beta(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_K(double T, double P, double n[5]) {
    return coder_K(T, P, n);
}

double Clinopyroxene_stx21_ph_coder_calib_Kp(double T, double P, double n[5]) {
    return coder_Kp(T, P, n);
}

int Clinopyroxene_stx21_ph_coder_get_param_number(void) {
    return coder_get_param_number();
}

const char **Clinopyroxene_stx21_ph_coder_get_param_names(void) {
    return coder_get_param_names();
}

const char **Clinopyroxene_stx21_ph_coder_get_param_units(void) {
    return coder_get_param_units();
}

void Clinopyroxene_stx21_ph_coder_get_param_values(double **values) {
    coder_get_param_values(values);
}

int Clinopyroxene_stx21_ph_coder_set_param_values(double *values) {
    return coder_set_param_values(values);
}

double Clinopyroxene_stx21_ph_coder_get_param_value(int index) {
    return coder_get_param_value(index);
}

int Clinopyroxene_stx21_ph_coder_set_param_value(int index, double value) {
    return coder_set_param_value(index, value);
}

double Clinopyroxene_stx21_ph_coder_dparam_g(double T, double P, double n[5], int index) {
    return coder_dparam_g(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_dgdt(double T, double P, double n[5], int index) {
    return coder_dparam_dgdt(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_dgdp(double T, double P, double n[5], int index) {
    return coder_dparam_dgdp(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d2gdt2(double T, double P, double n[5], int index) {
    return coder_dparam_d2gdt2(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d2gdtdp(double T, double P, double n[5], int index) {
    return coder_dparam_d2gdtdp(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d2gdp2(double T, double P, double n[5], int index) {
    return coder_dparam_d2gdp2(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d3gdt3(double T, double P, double n[5], int index) {
    return coder_dparam_d3gdt3(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d3gdt2dp(double T, double P, double n[5], int index) {
    return coder_dparam_d3gdt2dp(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d3gdtdp2(double T, double P, double n[5], int index) {
    return coder_dparam_d3gdtdp2(T, P, n, index);
}

double Clinopyroxene_stx21_ph_coder_dparam_d3gdp3(double T, double P, double n[5], int index) {
    return coder_dparam_d3gdp3(T, P, n, index);
}

void Clinopyroxene_stx21_ph_coder_dparam_dgdn(double T, double P, double n[5], int index, double result[5]) {
    coder_dparam_dgdn(T, P, n, index, result);
}

