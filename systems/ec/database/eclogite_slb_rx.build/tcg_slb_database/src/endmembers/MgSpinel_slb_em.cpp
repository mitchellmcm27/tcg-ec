#include "endmembers/MgSpinel_slb_em.h"
#include "tcgversion.h"

extern "C" {
#include "coder_error.h"
}

//-----------------------------------------------------------------------------
MgSpinel_slb_em::MgSpinel_slb_em()
{
  // Do Nothing
}
//-----------------------------------------------------------------------------
MgSpinel_slb_em::~MgSpinel_slb_em()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::identifier()
{
    std::string _str(MgSpinel_slb_em_coder_calib_identifier());
    return _str;
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::name()
{
    std::string _str(MgSpinel_slb_em_coder_calib_name());
    return _str;
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::tcg_build_version()
{
    return TCG_VERSION;
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::tcg_build_git_sha()
{
    return TCG_GIT_SHA;
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::tcg_generation_version()
{
    return "0.6.8+";
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::tcg_generation_git_sha()
{
    return "b755c9096e54ee3dcd3174bca6e1ee78d71068b8 Mon Jan 17 21:15:18 2022 -0500";
}
//-----------------------------------------------------------------------------
int MgSpinel_slb_em::check_coder_error()
{
    return return_coder_error_flag();
}
//-----------------------------------------------------------------------------
void MgSpinel_slb_em::reset_coder_error()
{
    set_coder_error_flag(0);
}
//-----------------------------------------------------------------------------
std::string MgSpinel_slb_em::formula()
{
    std::string _str(MgSpinel_slb_em_coder_calib_formula());
    return _str;
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::molecular_weight()
{
    return MgSpinel_slb_em_coder_calib_mw();
}
//-----------------------------------------------------------------------------
std::vector<double> MgSpinel_slb_em::elements()
{
  std::vector<double> _elements;
  const double *el = MgSpinel_slb_em_coder_calib_elements();
  _elements.assign(el, el + 106);
  return _elements;
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::G(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_g(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dGdT(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_dgdt(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dGdP(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_dgdp(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d2GdT2(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d2gdt2(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d2GdTdP(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d2gdtdp(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d2GdP2(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d2gdp2(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d3GdT3(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d3gdt3(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d3GdT2dP(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d3gdt2dp(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d3GdTdP2(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d3gdtdp2(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::d3GdP3(const double &T, const double &P)
{
  return MgSpinel_slb_em_coder_calib_d3gdp3(T,P);
}
//**************************************************************************
// Convenience functions of T and P
//**************************************************************************

//-----------------------------------------------------------------------------
double MgSpinel_slb_em::S(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_s(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::V(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_v(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dVdT(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_d2gdtdp(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dVdP(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_d2gdp2(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::Cv(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_cv(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::Cp(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_cp(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dCpdT(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_dcpdt(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::alpha(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_alpha(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::beta(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_beta(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::K(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_K(T,P);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::Kp(const double& T, const double& P)
{
  return MgSpinel_slb_em_coder_calib_Kp(T,P);
}
//**************************************************************************
// Active parameter functions directly from coder
//**************************************************************************

//-----------------------------------------------------------------------------
int MgSpinel_slb_em::get_param_number()
{
   return MgSpinel_slb_em_coder_get_param_number();
}
//-----------------------------------------------------------------------------
std::vector<std::string> MgSpinel_slb_em::get_param_names()
{
  std::vector<std::string> _param_names;
  const char **p = MgSpinel_slb_em_coder_get_param_names();
  _param_names.assign(p, p + MgSpinel_slb_em_coder_get_param_number());
  return _param_names;
}
//-----------------------------------------------------------------------------
std::vector<std::string> MgSpinel_slb_em::get_param_units()
{
  std::vector<std::string> _param_units;
  const char **p = MgSpinel_slb_em_coder_get_param_units();
  _param_units.assign(p, p + MgSpinel_slb_em_coder_get_param_number());
  return _param_units;
}
//-----------------------------------------------------------------------------
std::vector<double> MgSpinel_slb_em::get_param_values()
{
  std::vector<double> values(MgSpinel_slb_em_coder_get_param_number());
  double* v = values.data();
  double** v_ptr = &v;
  MgSpinel_slb_em_coder_get_param_values(v_ptr);
  return values;
}
//-----------------------------------------------------------------------------
void MgSpinel_slb_em::get_param_values(std::vector<double>& values)
{
  double* v = values.data();
  double** v_ptr = &v;
  MgSpinel_slb_em_coder_get_param_values(v_ptr);
}
//-----------------------------------------------------------------------------
int MgSpinel_slb_em::set_param_values(std::vector<double>& values)
{
  MgSpinel_slb_em_coder_set_param_values(values.data());
  return 1;
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::get_param_value(int& index)
{
  return MgSpinel_slb_em_coder_get_param_value(index);
}
//-----------------------------------------------------------------------------
int MgSpinel_slb_em::set_param_value(int& index, double& value)
{
  return MgSpinel_slb_em_coder_set_param_value(index,value);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_g(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_g(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_dgdt(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_dgdt(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_dgdp(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_dgdp(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d2gdt2(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d2gdt2(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d2gdtdp(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d2gdtdp(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d2gdp2(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d2gdp2(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d3gdt3(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d3gdt3(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d3gdt2dp(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d3gdt2dp(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d3gdtdp2(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d3gdtdp2(T,P,index);
}
//-----------------------------------------------------------------------------
double MgSpinel_slb_em::dparam_d3gdp3(double& T, double& P, int& index)
{
  return MgSpinel_slb_em_coder_dparam_d3gdp3(T,P,index);
}
//-----------------------------------------------------------------------------
