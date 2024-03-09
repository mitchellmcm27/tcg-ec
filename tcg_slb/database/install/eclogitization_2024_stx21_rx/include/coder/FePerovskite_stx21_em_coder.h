
const char *FePerovskite_stx21_em_coder_calib_identifier(void);
const char *FePerovskite_stx21_em_coder_calib_name(void);
const char *FePerovskite_stx21_em_coder_calib_formula(void);
const double FePerovskite_stx21_em_coder_calib_mw(void);
const double *FePerovskite_stx21_em_coder_calib_elements(void);

double FePerovskite_stx21_em_coder_calib_g(double T, double P);
double FePerovskite_stx21_em_coder_calib_dgdt(double T, double P);
double FePerovskite_stx21_em_coder_calib_dgdp(double T, double P);
double FePerovskite_stx21_em_coder_calib_d2gdt2(double T, double P);
double FePerovskite_stx21_em_coder_calib_d2gdtdp(double T, double P);
double FePerovskite_stx21_em_coder_calib_d2gdp2(double T, double P);
double FePerovskite_stx21_em_coder_calib_d3gdt3(double T, double P);
double FePerovskite_stx21_em_coder_calib_d3gdt2dp(double T, double P);
double FePerovskite_stx21_em_coder_calib_d3gdtdp2(double T, double P);
double FePerovskite_stx21_em_coder_calib_d3gdp3(double T, double P);

double FePerovskite_stx21_em_coder_calib_s(double T, double P);
double FePerovskite_stx21_em_coder_calib_v(double T, double P);
double FePerovskite_stx21_em_coder_calib_cv(double T, double P);
double FePerovskite_stx21_em_coder_calib_cp(double T, double P);
double FePerovskite_stx21_em_coder_calib_dcpdt(double T, double P);
double FePerovskite_stx21_em_coder_calib_alpha(double T, double P);
double FePerovskite_stx21_em_coder_calib_beta(double T, double P);
double FePerovskite_stx21_em_coder_calib_K(double T, double P);
double FePerovskite_stx21_em_coder_calib_Kp(double T, double P);

int FePerovskite_stx21_em_coder_get_param_number(void);
const char **FePerovskite_stx21_em_coder_get_param_names(void);
const char **FePerovskite_stx21_em_coder_get_param_units(void);
void FePerovskite_stx21_em_coder_get_param_values(double **values);
int FePerovskite_stx21_em_coder_set_param_values(double *values);
double FePerovskite_stx21_em_coder_get_param_value(int index);
int FePerovskite_stx21_em_coder_set_param_value(int index, double value);

double FePerovskite_stx21_em_coder_dparam_g(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_dgdt(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_dgdp(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d2gdt2(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d2gdtdp(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d2gdp2(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d3gdt3(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d3gdt2dp(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d3gdtdp2(double T, double P, int index);
double FePerovskite_stx21_em_coder_dparam_d3gdp3(double T, double P, int index);

