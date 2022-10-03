#include <math.h>
#include <numeric>
#include <algorithm>
#include <functional>
#include "phases/Almandine_slb_ph.h"
#include "tcgversion.h"

extern "C" {
#include "coder_error.h"
}

//-----------------------------------------------------------------------------
Almandine_slb_ph::Almandine_slb_ph()
: _endmembers(std::vector<std::shared_ptr<EndMember> > {std::make_shared<Almandine_slb_em>() })
{
  // Set number of components
  C = _endmembers.size();

  // Assign vector of molecular weights for each component
  M = Almandine_slb_ph::get_M();
  // Calculate inverse Masses
  iM.resize(C, 0.);
  for (int i = 0; i < C; i++)
    iM[i] = 1./M[i];

  // allocate temporary storage vectors for vecs and  matrices
  _tmp.resize(C, 0.);
  _result.resize(C, 0.);
  _result_mat.resize(C*(C+1)/2, 0.);
  // temporary element vector
  _elm.resize(106, 0.);
  // compositional Jacobian
  _dmu_dc.resize(C);
  for (int i = 0; i < C; i++ )
  {
    _dmu_dc[i].resize(C, 0.);
  }
}
//-----------------------------------------------------------------------------
Almandine_slb_ph::~Almandine_slb_ph()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::identifier()
{
    std::string _str(Almandine_slb_ph_coder_calib_identifier());
    return _str;
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::name()
{
  std::string _str(Almandine_slb_ph_coder_calib_name());
  return _str;
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::tcg_build_version()
{
    return TCG_VERSION;
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::tcg_build_git_sha()
{
    return TCG_GIT_SHA;
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::tcg_generation_version()
{
    return "0.6.8+";
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::tcg_generation_git_sha()
{
    return "b755c9096e54ee3dcd3174bca6e1ee78d71068b8 Mon Jan 17 21:15:18 2022 -0500";
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::check_coder_error()
{
    return return_coder_error_flag();
}
//-----------------------------------------------------------------------------
void Almandine_slb_ph::reset_coder_error()
{
    set_coder_error_flag(0);
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::formula(
    const double &T, const double &P,
    std::vector<double> &n)
{
  std::string _str(Almandine_slb_ph_coder_calib_formula(T,P,n.data()));
  return _str;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::conv_elm_to_moles(std::vector<double>& e)
{
  std::vector<double> n;
  double *_n = Almandine_slb_ph_coder_calib_conv_elm_to_moles(e.data());
  n.assign(_n, _n + C);
  return n;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::conv_elm_to_tot_moles(std::vector<double>& e)
{
  _tmp = conv_elm_to_moles(e);
  double result = std::accumulate(_tmp.begin(), _tmp.end(), 0.0);
  return result;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::conv_elm_to_tot_grams(std::vector<double>& e)
{
  _tmp  = conv_elm_to_moles(e);
  double result = std::inner_product(_tmp.begin(), _tmp.end(), M.begin(), 0.0);
  return result;
}

//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::conv_moles_to_elm(std::vector<double>& n)
{
  std::vector<double> result(106, 0.);
  for (int i = 0; i< C; i++)
  {
    _elm  = endmember_elements(i);
    for (int j = 0; j < result.size(); j++)
    {
    //FIXME:: could use std::transform for this but its fugly
            result[j] += n[i]*_elm[j];
    }
  }
  return result;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::conv_moles_to_tot_moles(std::vector<double>& n)
{
  return std::accumulate(n.begin(), n.end(), 0.0);
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::conv_moles_to_mole_frac(std::vector<double>& n)
{
  double n_tot = conv_moles_to_tot_moles(n);
  for (int i = 0; i < C; i++)
  {
    _tmp[i] = n[i]/n_tot;
  }
  return _tmp;
}

//-----------------------------------------------------------------------------
int Almandine_slb_ph::test_moles(std::vector<double>& n)
{
  return Almandine_slb_ph_coder_calib_test_moles(n.data());
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::endmember_number() const
{
  return _endmembers.size();
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::endmember_name(const int &i) const
{
  return (*_endmembers[i]).name();
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::endmember_formula(const int &i) const
{
  return (*_endmembers[i]).formula();
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::endmember_mw(const int &i) const
{
  return (*_endmembers[i]).molecular_weight();
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::endmember_elements(const int &i) const
{
  return (*_endmembers[i]).elements();
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::species_number() const
{
  return _endmembers.size();
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::species_name(const int &i) const
{
  return (*_endmembers[i]).name();
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::species_formula(const int &i) const
{
  return (*_endmembers[i]).formula();
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::species_mw(const int &i) const
{
  return (*_endmembers[i]).molecular_weight();
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::species_elements(const int &i) const
{
  return (*_endmembers[i]).elements();
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::g(
    const double &T, const double &P,
    std::vector<double>& n) const
{
  return Almandine_slb_ph_coder_calib_g(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dgdt(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_dgdt(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dgdp(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_dgdp(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d2gdt2(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d2gdt2(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d2gdtdp(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d2gdtdp(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d2gdp2(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d2gdp2(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d3gdt3(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d3gdt3(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d3gdt2dp(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d3gdt2dp(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d3gdtdp2(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d3gdtdp2(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::d3gdp3(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_d3gdp3(T,P,n.data());
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::dgdn(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_dgdn(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d2gdndt(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d2gdndt(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d2gdndp(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d2gdndp(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d2gdn2(
    const double &T, const double &P, std::vector<double> &n) const
{
  //std::vector<double> _d2gdn2(C*(C+1)/2);
  Almandine_slb_ph_coder_calib_d2gdn2(T,P,n.data(),_result_mat.data());
  return _result_mat;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdndt2(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d3gdndt2(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdndtdp(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d3gdndtdp(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdn2dt(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d3gdn2dt(T,P,n.data(),_result_mat.data());
  return _result_mat;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdndp2(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d3gdndp2(T,P,n.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdn2dp(
    const double &T, const double &P, std::vector<double> &n) const
{
  Almandine_slb_ph_coder_calib_d3gdn2dp(T,P,n.data(),_result_mat.data());
  return _result_mat;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::d3gdn3(
    const double &T, const double &P, std::vector<double> &n) const
{
  //FIXME: see if we need this vector as well
  std::vector<double> _d3gdn3(C*(C+1)*(C+2)/6);
  Almandine_slb_ph_coder_calib_d3gdn3(T,P,n.data(),_d3gdn3.data());
  return _d3gdn3;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::v(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_v(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::s(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_s(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::alpha(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_alpha(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::cv(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_cv(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::cp(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_cp(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dcpdt(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_dcpdt(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::beta(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_beta(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::K(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_K(T,P,n.data());
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::Kp(
    const double &T, const double &P,
    std::vector<double> &n) const
{
  return Almandine_slb_ph_coder_calib_Kp(T,P,n.data());
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::get_param_number()
{
   return Almandine_slb_ph_coder_get_param_number();
}
//-----------------------------------------------------------------------------
std::vector<std::string> Almandine_slb_ph::get_param_names()
{
  std::vector<std::string> _param_names;
  const char **p = Almandine_slb_ph_coder_get_param_names();
  _param_names.assign(p, p + Almandine_slb_ph_coder_get_param_number());
  return _param_names;
}
//-----------------------------------------------------------------------------
std::vector<std::string> Almandine_slb_ph::get_param_units()
{
  std::vector<std::string> _param_units;
  const char **p = Almandine_slb_ph_coder_get_param_units();
  _param_units.assign(p, p + Almandine_slb_ph_coder_get_param_number());
  return _param_units;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::get_param_values()
{
  std::vector<double> values(Almandine_slb_ph_coder_get_param_number());
  double* v = values.data();
  double** v_ptr = &v;
  Almandine_slb_ph_coder_get_param_values(v_ptr);
  return values;
}
//-----------------------------------------------------------------------------
void Almandine_slb_ph::get_param_values(std::vector<double>& values)
{
  double* v = values.data();
  double** v_ptr = &v;
  Almandine_slb_ph_coder_get_param_values(v_ptr);
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::set_param_values(std::vector<double>& values)
{
  Almandine_slb_ph_coder_set_param_values(values.data());
  return 1;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::get_param_value(int& index)
{
  return Almandine_slb_ph_coder_get_param_value(index);
}
//-----------------------------------------------------------------------------
int Almandine_slb_ph::set_param_value(int& index, double& value)
{
  return Almandine_slb_ph_coder_set_param_value(index,value);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_g(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_g(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_dgdt(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_dgdt(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_dgdp(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_dgdp(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::dparam_dgdn(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  Almandine_slb_ph_coder_dparam_dgdn(T,P,n.data(),index,_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d2gdt2(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d2gdt2(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d2gdtdp(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d2gdtdp(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d2gdp2(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d2gdp2(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d3gdt3(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d3gdt3(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d3gdt2dp(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d3gdt2dp(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d3gdtdp2(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d3gdtdp2(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::dparam_d3gdp3(
    double& T, double& P,
    std::vector<double> &n, int& index)
{
  return Almandine_slb_ph_coder_dparam_d3gdp3(T,P,n.data(),index);
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<EndMember> > Almandine_slb_ph::endmembers() const
{
  return _endmembers;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::get_M() const
{
  std::vector<double> M(C);

  for(int k = 0; k < C; k++)
  {
    M[k] = (*_endmembers[k]).molecular_weight();
  }
  return M;
}
//-----------------------------------------------------------------------------
std::string Almandine_slb_ph::abbrev() const
{
  return "al";
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::mu(
    const double &T, const double &P, std::vector<double> &x) const
{
  Almandine_slb_ph_coder_calib_dgdn(T,P,x.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::dmu_dT(
    const double &T, const double &P, std::vector<double> &x) const
{
  Almandine_slb_ph_coder_calib_d2gdndt(T,P,x.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::dmu_dP(
    const double &T, const double &P, std::vector<double> &x) const
{
  Almandine_slb_ph_coder_calib_d2gdndp(T,P,x.data(),_result.data());
  return _result;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > Almandine_slb_ph::dmu_dc(
    const double &T, const double &P, std::vector<double> &c) const
{
  // if C < 1 _dmu_dc = {{0.}} by initialization already
  if ( C > 1 )
  {
    // Convert c to x and calculate c.(1/M)
    c_to_x(c, _tmp);
    double cdotiM = std::inner_product(c.begin(), c.end(), iM.begin(), 0.0);

    // extract compressed dmu/dn matrix from coder
    Almandine_slb_ph_coder_calib_d2gdn2(T,P,_tmp.data(),_result_mat.data());

    // loop over components
    for (int j = 0; j < C; j++)
    {
        // map from compressed matrix to calculate \grad_n mu^j
      for(int k = 0; k < C; k++)
      {
        int m = (j <= k) ? j*(C - 1) - (j - 1)*j/2 + k : k*(C - 1) - (k - 1)*k/2 + j ;
        _dmu_dc[j][k] = _result_mat[m];
      }
      // right multiply by  jacobian dx/dc = 1/cdotiM * ( I - x 1^T) diag(iM)
      double gdotx = std::inner_product(_dmu_dc[j].begin(), _dmu_dc[j].end(), _tmp.begin(), 0.0);
      for (int k = 0; k < C; k++)
      {
        _dmu_dc[j][k] = ( _dmu_dc[j][k] - gdotx ) * iM[k] /cdotiM;
      }
    }
  }
  return _dmu_dc;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::Mass(std::vector<double> &x) const
{
  //return dot product of x M
  double _M = std::inner_product(x.begin(),x.end(), M.begin(), 0.0);
  return _M;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::rho(
    const double &T, const double &P,
    std::vector<double> &c) const
{
  c_to_x(c,_tmp);
  double _rho = Mass(_tmp)/v(T,P,_tmp);
  return _rho;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::drho_dT(
    const double &T, const double &P,
    std::vector<double> &c) const
{
  c_to_x(c, _tmp);
  double _V = v(T,P,_tmp);
  double _drho_dT = -d2gdtdp(T,P,_tmp)*Mass(_tmp)/(_V*_V);
  return _drho_dT;
}
//-----------------------------------------------------------------------------
double Almandine_slb_ph::drho_dP(
    const double &T, const double &P,
    std::vector<double> &c) const
{
  c_to_x(c, _tmp);
  double _V = v(T,P,_tmp);
  double _drho_dP = -d2gdp2(T,P,_tmp)*Mass(_tmp)/(_V*_V);
  return _drho_dP;
}
//-----------------------------------------------------------------------------
std::vector<double>  Almandine_slb_ph::drho_dc(
    const double &T, const double &P,
    std::vector<double> &c) const
{
  // if C=1, dv_dC will return [0].
  std::vector<double> _drho_dc = dv_dc(T, P, c);
  if ( C > 1 )
  {
    c_to_x(c, _tmp);
    double _V = v(T,P,_tmp);
    double _M = Mass(_tmp);
    double cdotiM = std::inner_product(c.begin(), c.end(), iM.begin(), 0.0);

    for (int k = 0; k < C; k++)
    {
        _drho_dc[k] = (-_M/_V*_drho_dc[k] +  (1. - _M/M[k])/cdotiM)/_V;
    }
  }
  return _drho_dc;
}
//-----------------------------------------------------------------------------
std::vector<double>  Almandine_slb_ph::ds_dc(
    const double &T, const double &P, std::vector<double> &c) const
{
  std::vector<double> _ds_dc(C, 0.);
  if ( C > 1 )
  {
    // Convert c to x and calculate c.(1/M)
    c_to_x(c, _tmp);
    double cdotiM = std::inner_product(c.begin(), c.end(), iM.begin(), 0.0);

    // extract  -ds/dx vector from coder
    Almandine_slb_ph_coder_calib_d2gdndt(T,P,_tmp.data(),_result.data());

    // right multiply by  jacobian dx/dc = 1/cdotiM * ( I - x 1^T) diag(iM)
    double gdotx = -std::inner_product(_result.begin(), _result.end(), _tmp.begin(), 0.0);
    for (int k = 0; k < C; k++)
    {
      _ds_dc[k] = ( -_result[k] - gdotx ) * iM[k] /cdotiM;
    }
  }
  return _ds_dc;
}
//-----------------------------------------------------------------------------
std::vector<double>  Almandine_slb_ph::dv_dc(
    const double &T, const double &P, std::vector<double> &c) const
{
  std::vector<double> _dv_dc(C, 0.);
  if ( C > 1 )
  {
    // Convert c to x and calculate c.(1/M)
    c_to_x(c, _tmp);
    double cdotiM = std::inner_product(c.begin(), c.end(), iM.begin(), 0.0);

    // extract  dv/dx vector from coder
    Almandine_slb_ph_coder_calib_d2gdndp(T,P,_tmp.data(),_result.data());

    // right multiply by  jacobian dx/dc = 1/cdotiM * ( I - x 1^T) diag(iM)
    double gdotx = std::inner_product(_result.begin(), _result.end(), _tmp.begin(), 0.0);
    for (int k = 0; k < C; k++)
    {
      _dv_dc[k] = ( _result[k] - gdotx ) * iM[k] /cdotiM;
    }
  }
  return _dv_dc;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::c_to_x(std::vector<double> &c) const
{
  c_to_x(c, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
std::vector<double> Almandine_slb_ph::x_to_c(std::vector<double> &x) const
{
  x_to_c(x, _tmp);
  return _tmp;
}
//pass by reference versions
//-----------------------------------------------------------------------------
void  Almandine_slb_ph::c_to_x(std::vector<double> &c, std::vector<double> &x) const
{
  // stl version of c to x calculate x_i = c_i/M_i / (sum c_k/M_k)
  std::transform(c.begin(), c.end(), iM.begin(), x.begin(), std::multiplies<double>());
  double sum  = std::accumulate(x.begin(), x.end(), 0.0);
  std::transform(x.begin(), x.end(), x.begin(),
               std::bind(std::divides<double>(), std::placeholders::_1, sum));
}
//-----------------------------------------------------------------------------
void  Almandine_slb_ph::x_to_c(std::vector<double> &x, std::vector<double> &c) const
{
  // stl version of x to c calculate c_i = x_i*M_i / (sum x_k*M_k)
  std::transform(x.begin(), x.end(), M.begin(), c.begin(), std::multiplies<double>());
  double sum  = std::accumulate(c.begin(), c.end(), 0.0);
  std::transform(c.begin(), c.end(), c.begin(),
               std::bind(std::divides<double>(), std::placeholders::_1, sum));
}
//-----------------------------------------------------------------------------
