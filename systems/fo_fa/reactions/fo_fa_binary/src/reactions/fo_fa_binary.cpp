#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <functional>
#include "reactions/fo_fa_binary.h"
#include "tcgversion.h"

//-----------------------------------------------------------------------------
fo_fa_binary::fo_fa_binary()
: _name("fo_fa_binary"), J(2),  _phases(std::vector<std::shared_ptr<Phase> > {std::make_shared<Olivine>(),std::make_shared<Liquid>() }), _C0({{0.0, 0.0}, {0.0, 0.0}}), _C({{0.0, 0.0}, {0.0, 0.0}}),_X({{0.0, 0.0}, {0.0, 0.0}}),
  _Mu({{0.0, 0.0}, {0.0, 0.0}}), _dMu_dT({{0.0, 0.0}, {0.0, 0.0}}), _dMu_dP({{0.0, 0.0}, {0.0, 0.0}}),_tmp_ik({{0.0, 0.0}, {0.0, 0.0}}),
  _A(J,0.), _dA_dT(_A), _dA_dP(_A), _dA_dC(J,{{0.0, 0.0}, {0.0, 0.0}}),
  _M({{140.69332, 203.77312}, {140.69332, 203.77312}}),
  _nu({{{-1.0, 0.0}, {1.0, 0.0}}, {{0.0, -1.0}, {0.0, 1.0}}}),
  _nu_m({{{-1.0, 0.0}, {1.0, 0.0}}, {{0.0, -1.0}, {0.0, 1.0}}}),
  parameters({{"T0", new double}, {"R", new double}}),
  T0(1000.0),
  R(8.314)
{
  // allocate  temporary vectors for Gamma and _dMu_dC construct
  N = _phases.size();
  _tmp.resize(N);
  _dMu_dC.resize(N);
  _tmp_dC.resize(N);
  for (int i = 0; i < N; i++)
  {
    int K = _C0[i].size();
    _dMu_dC[i].resize(K);
    _tmp_dC[i].resize(K);
    for (int k = 0; k < K; k++)
    {
       _dMu_dC[i][k].resize(K, 0.0);
       _tmp_dC[i][k].resize(K, 0.0);
    }
  }

  // Set parameter values in map container
  parameters["T0"] = &T0;
  parameters["R"] = &R;
}
//-----------------------------------------------------------------------------
fo_fa_binary::~fo_fa_binary()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::name() const
{
  return _name;
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::tcg_build_version()
{
    return TCG_VERSION;
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::tcg_build_git_sha()
{
    return TCG_GIT_SHA;
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::tcg_generation_version()
{
    return "0.6.9";
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::tcg_generation_git_sha()
{
    return "117d758197e2d445579ad671f573064c0650429d Mon Aug 1 00:36:24 2022 +0000";
}
//-----------------------------------------------------------------------------
 std::vector<std::shared_ptr<Phase> > fo_fa_binary::phases() const
{
  return _phases;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::zero_C() const
{
  return _C0;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::M() const
{
  return _M;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > fo_fa_binary::nu() const
{
  return _nu;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > fo_fa_binary::nu_m() const
{
  return _nu_m;
}
//-----------------------------------------------------------------------------
// pass by reference interface
void fo_fa_binary::A(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_A) const
{
  //Calculate chemical potential matrix M_ik
  Mu(T, P, C, _Mu);

  _A = {
  1.0*_Mu[0][0] + -1.0*_Mu[1][0],
  1.0*_Mu[0][1] + -1.0*_Mu[1][1]};
}
//-----------------------------------------------------------------------------
// return  vector<double>
std::vector<double> fo_fa_binary::A(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  A(T, P, C, _A);
  return _A;
}

//-----------------------------------------------------------------------------
// overloaded version to just return the jth affinity A_j
double fo_fa_binary::A(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  // FIXME:: this call should probably should be deprecated
  // it is no longer calling function A(T,P,C,j) just calling void A and returning
  // _A[j] just to keep interface and consistency with new affinity functions
  A(T, P, C, _A);
  return _A[j];
}
//-----------------------------------------------------------------------------
// pass by reference interface
void fo_fa_binary::dA_dT(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dT) const
{
  //Calculate chemical potential matrix dMu_dT_ik
  dMu_dT(T, P, C, _dMu_dT);

  _dA_dT = {
  1.0*_dMu_dT[0][0] + -1.0*_dMu_dT[1][0],
  1.0*_dMu_dT[0][1] + -1.0*_dMu_dT[1][1]};
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT;
}

//-----------------------------------------------------------------------------
double fo_fa_binary::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT[j];
}
//-----------------------------------------------------------------------------
// pass by reference interface
void fo_fa_binary::dA_dP(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dP) const
{
  //Calculate chemical potential matrix dMu_dP_ik
  dMu_dP(T, P, C, _dMu_dP);

  _dA_dP = {
  1.0*_dMu_dP[0][0] + -1.0*_dMu_dP[1][0],
  1.0*_dMu_dP[0][1] + -1.0*_dMu_dP[1][1]};
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP[j];
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dA_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<std::vector<std::vector<double> > >& _dA_dC) const
{
    dMu_dC(T, P, C, _dMu_dC);
    // loop over reactions
    for (int j=0; j < J; j++)
    {
        // loop over phases
        for (int i = 0; i< N; i++)
        {   // zero _dA_dC
            _dA_dC[j][i] = _C0[i];
            int K = _C0[i].size();
            for (int k = 0; k< K; k++)
            {   //FIXME: maybe make this a sparse matvec without the if statement?
                //note for affinities this should use _nu not _nu_m
                if (_nu[j][i][k] != 0.0)
                {
                    for (int l = 0; l < K; l++)
                    {   //FIXME:  could also think about stl: transform here
                        _dA_dC[j][i][l] -= _nu[j][i][k]*_dMu_dC[i][k][l];
                    }
                }
            }
        }
    }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > fo_fa_binary::dA_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  dA_dC( T, P, C, _dA_dC);
  return _dA_dC;
}

//-----------------------------------------------------------------------------
double fo_fa_binary::dAj_dCik(
    const double &T, const double &P,
	  std::vector<std::vector<double> >& C,
		unsigned j, unsigned  i, unsigned k) const
{
  //FIXME: should be deprecated,  removed old function call version as its incorrect,  now is just redundant
  // to previous call but just returns one component
  dA_dC( T, P, C, _dA_dC);
  return _dA_dC[j][i][k];
}
//-----------------------------------------------------------------------------
void fo_fa_binary::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double> &_Gamma) const
{
  //calculate current affinities
  A(T, P, C, _A);
  _Gamma = {
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
)),
  1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  Gamma_i(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  double Gamma_i;

  switch(i) {
    case 0: Gamma_i = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
)); break;
    case 1: Gamma_i = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
)); break;
  }
  return Gamma_i;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)),
  1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
))};
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dT (T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  double dGamma_i_dT;

  switch(i) {
    case 0: dGamma_i_dT = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)); break;
    case 1: dGamma_i_dT = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_i_dT;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
)),
  1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dP(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  double dGamma_i_dP;

  switch(i) {
    case 0: dGamma_i_dP = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
)); break;
    case 1: dGamma_i_dP = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
)); break;
  }
  return dGamma_i_dP;
}

//-----------------------------------------------------------------------------
void  fo_fa_binary::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<std::vector<double> > >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  _dGamma = {
  {{-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
))}},
  {{1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
))},
   {1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
))}}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > fo_fa_binary::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dC( T, P, C, Phi, _tmp_dC);
  return _tmp_dC;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_i_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l, unsigned k) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned m = i*4 + l*2 + k;
  double dGamma_i_dC;

  switch(m) {
    case 0: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)); break;
    case 1: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)); break;
    case 2: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)); break;
    case 3: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)); break;
    case 4: dGamma_i_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)); break;
    case 5: dGamma_i_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)); break;
    case 6: dGamma_i_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)); break;
    case 7: dGamma_i_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)); break;
  }
  return dGamma_i_dC;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
   _dGamma = {
  {-1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))},
  {1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  std::vector<std::vector<double> > _dGamma_dPhi = {
  {-1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))},
  {1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))}};
  return _dGamma_dPhi;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = i*2 + l;
  double dGamma_i_dPhi;

  switch(m) {
    case 0: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 1: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
    case 2: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 3: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
  }
  return dGamma_i_dPhi;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<double> > &_Gamma) const
{
  A(T, P, C, _A);
  _Gamma = {
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)), -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
))},
  {1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)), 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  Gamma_ik(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  unsigned m = i*2 + k;
  double Gamma_ik;

  switch(m) {
    case 0: Gamma_ik = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)); break;
    case 1: Gamma_ik = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
)); break;
    case 2: Gamma_ik = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*T)
)); break;
    case 3: Gamma_ik = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*T)
)); break;
  }
  return Gamma_ik;
}
//-----------------------------------------------------------------------------
void  fo_fa_binary::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)), -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
))},
  {1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)), 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: using_tmp_ik but starting to worry about thread safety
  dGamma_ik_dT(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  unsigned m = i*2 + k;
  double dGamma_ik_dT;

  switch(m) {
    case 0: dGamma_ik_dT = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)); break;
    case 1: dGamma_ik_dT = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)); break;
    case 2: dGamma_ik_dT = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[0]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[0]/(R*((T)*(T)*(T)))
)); break;
    case 3: dGamma_ik_dT = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dT[1]/(R*T) - std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_A[1]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_ik_dT;
}
//-----------------------------------------------------------------------------
void  fo_fa_binary::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)), -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
))},
  {1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)), 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> >  fo_fa_binary::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: thread safety again?
  dGamma_ik_dP(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  unsigned m = i*2 + k;
  double dGamma_ik_dP;

  switch(m) {
    case 0: dGamma_ik_dP = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)); break;
    case 1: dGamma_ik_dP = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
)); break;
    case 2: dGamma_ik_dP = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[0]/(R*T)
)); break;
    case 3: dGamma_ik_dP = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dP[1]/(R*T)
)); break;
  }
  return dGamma_ik_dP;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);

    switch(i) {
    case 0: _dGamma = {
 {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)),
 -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
))},
        {-1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)),
 -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
))}
    }; break;

    case 1: _dGamma = {
 {1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)),
 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
))},
        {1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)),
 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, int i) const
{
  int K = _C0[i].size();
  std::vector<std::vector<double> > _dGamma(K,std::vector<double>(K,0.));
  dGamma_ik_dC(T, P, C, Phi, i, _dGamma);
  return _dGamma;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
    switch(i) {
    case 0: _dGamma = {
 {-1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)),
 -1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
))},
        {-1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)),
 -1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))}
    }; break;

    case 1: _dGamma = {
 {1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)),
 1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
))},
        {1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)),
 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, int i) const
{
  int K = _C0[i].size();
  std::vector<std::vector<double> > _dGamma(K,std::vector<double>(N,0.));
  dGamma_ik_dPhi(T, P, C, Phi, i, _dGamma);
  return _dGamma;
}

//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l, unsigned m) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned p = l*8 + m*4 + i*2+ k;
  double dGamma_ik_dC;

  switch(p) {
    case 0: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)); break;
    case 1: dGamma_ik_dC = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)); break;
    case 2: dGamma_ik_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][0]/(R*T)
)); break;
    case 3: dGamma_ik_dC = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][0]/(R*T)
)); break;
    case 4: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)); break;
    case 5: dGamma_ik_dC = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)); break;
    case 6: dGamma_ik_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][0][1]/(R*T)
)); break;
    case 7: dGamma_ik_dC = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][0][1]/(R*T)
)); break;
    case 8: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)); break;
    case 9: dGamma_ik_dC = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)); break;
    case 10: dGamma_ik_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][0]/(R*T)
)); break;
    case 11: dGamma_ik_dC = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][0]/(R*T)
)); break;
    case 12: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)); break;
    case 13: dGamma_ik_dC = -1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)); break;
    case 14: dGamma_ik_dC = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[0][1][1]/(R*T)
)); break;
    case 15: dGamma_ik_dC = 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*pow(Phi[0], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*pow(Phi[1], 1.0/2.0)*_dA_dC[1][1][1]/(R*T)
)); break;
  }
  return dGamma_ik_dC;
}
//-----------------------------------------------------------------------------
double fo_fa_binary::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = l*4 + i*2 + k;
  double dGamma_ik_dPhi;

  switch(m) {
    case 0: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 1: dGamma_ik_dPhi = -1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 2: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 3: dGamma_ik_dPhi = 1.0*((_A[1] >= 0) ? (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[0], 1.0/2.0))
)
: (
   0
)); break;
    case 4: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
    case 5: dGamma_ik_dPhi = -1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
    case 6: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[0]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
    case 7: dGamma_ik_dPhi = 1.0*((_A[1] >= 0) ? (
   0
)
: (
   (1.0/2.0)*std::exp(-T0/T)*_A[1]/(R*T*pow(Phi[1], 1.0/2.0))
)); break;
  }
  return dGamma_ik_dPhi;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::rho(const double &T, const double &P,
		std::vector<std::vector<double> >& C,
        std::vector<double>&  _rho) const
{
  // inplace calculation of vector of phase densities (must be correct size)
  for (int i = 0; i < _phases.size(); i++ ) {
    _rho[i] = _phases[i]->rho(T,P,C[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::rho(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _rho(2,0.);
  rho(T, P, C, _rho);
  return _rho;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::Cp(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C, std::vector<double>&  _Cp) const
{
  // inplace calculation of vector of  phase densities (must be correct size)
  // Convert C to X
  C_to_X(C,_X);

  for (int i = 0; i < _phases.size(); i++ ) {
    _Cp[i] = _phases[i]->cp(T,P,_X[i])/_phases[i]->Mass(_X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::Cp(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{

  // return vector of densities
  std::vector<double> _Cp(2,0.);
  Cp(T, P, C, _Cp);
  return _Cp;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::s(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>&  _s) const
{
  // inplace calculation of vector of  phase densities (must be correct size)

  // Convert C to X
  C_to_X(C,_X);

  for (int i = 0; i < _phases.size(); i++)
  {
    _s[i] = _phases[i]->s(T,P,_X[i])/_phases[i]->Mass(_X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::s(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _s(2,0.);
  s(T, P, C, _s);
  return _s;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::ds_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >&  _ds) const
{
  // inplace calculation of derivatives of entropy/unit mass with phase composition

  // Convert C to X
  C_to_X(C,_X);
  // zero out all derivatives
  _ds = _C0;
  for (int i = 0; i < _phases.size(); i++)
  {
    if (C[i].size() > 1)
    {
        _ds[i] = _phases[i]->ds_dc(T,P,C[i]);
        double cTiM = std::inner_product(C[i].begin(), C[i].end(),_M[i].begin(), 0.,
                                         std::plus<double>(), std::divides<double>());
        double Mass = std::inner_product(_X[i].begin(), _X[i].end(),_M[i].begin(), 0.);
        for (int k = 0; k < C[i].size(); k++)
        {
            _ds[i][k] = (_ds[i][k] + _phases[i]->s(T,P,_X[i]) / (cTiM * _M[i][k]) )/Mass;
        }
    }
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> >fo_fa_binary::ds_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  ds_dC(T, P, C, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::alpha(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<double>&  _alpha) const
{
  // inplace calculation of vector of  phase alphas (must be correct size)

  // Convert C to X
  C_to_X(C,_X);
  for (int i = 0; i < _phases.size(); i++)
  {
    _alpha[i] = _phases[i]->alpha(T,P,_X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::alpha(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _alpha(2,0.);
  alpha(T, P, C, _alpha);
  return _alpha;

}
//-----------------------------------------------------------------------------
void fo_fa_binary::beta(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<double>&  _beta) const
{
  // inplace calculation of vector of  phase alphas (must be correct size)

  // Convert C to X
  C_to_X(C,_X);
  for (int i = 0; i < _phases.size(); i++)
  {
    _beta[i] = _phases[i]->beta(T,P,_X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> fo_fa_binary::beta(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _beta(2,0.);
  beta(T, P, C, _beta);
  return _beta;

}

//-----------------------------------------------------------------------------
void fo_fa_binary::Mu(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _mu) const
{
  // inplace calculation of "Matrix" of chemical potentials Mu_i^k as vector of vectors (must be correct size)

  // Convert C to X
  C_to_X(C,_X);

  // loop over phases and insert chemical potential vectors

  for (int i = 0; i < _phases.size(); i++)
  {
    _mu[i] = _phases[i]->mu(T,P,_X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::Mu(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  Mu(T,P,C,_Mu);
  return _Mu;
}

//-----------------------------------------------------------------------------
void fo_fa_binary::dMu_dT(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const
{
  // inplace calculation of "Matrix" of chemical potentials Mu_i^k as vector of vectors (must be correct size)

  // Convert C to X
  C_to_X(C,_X);

  // loop over phases and insert chemical potential vectors

  for (int i = 0; i < _phases.size(); i++)
  {
    _dmu[i] = _phases[i]->dmu_dT(T,P,_X[i]);
  }
}
//-----------------------------------------------------------------------------
void fo_fa_binary::dMu_dP(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const
{
  // inplace calculation of "Matrix" of chemical potentials Mu_i^k as vector of vectors (must be correct size)

  // Convert C to X
  C_to_X(C,_X);

  // loop over phases and insert chemical potential vectors

  for (int i = 0; i < _phases.size(); i++)
  {
    _dmu[i] = _phases[i]->dmu_dP(T,P,_X[i]);
  }
}

//-----------------------------------------------------------------------------
void fo_fa_binary::dMu_dC(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C,
    std::vector<std::vector<std::vector<double> > >& _dmu) const
{
  // inplace calculation of "Matrix" of chemical potentials Mu_i^k as vector of vectors (must be correct size)

  // loop over phases and insert chemical potential vectors

  for (int i = 0; i < _phases.size(); i++)
  {
    _dmu[i] = _phases[i]->dmu_dc(T,P,C[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > fo_fa_binary::dMu_dC(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  dMu_dC(T,P,C,_dMu_dC);
  return _dMu_dC;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::C_to_X(
    std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> > &X) const
{
  //convert from weight fractions to Mole fractions (assumes Containers exist for both)
  for (int i = 0; i < _phases.size(); i++)
  {
    X[i] = _phases[i]->c_to_x(C[i]);
  }
}
//-----------------------------------------------------------------------------
void fo_fa_binary::X_to_C(
    std::vector<std::vector<double> >& X,
    std::vector<std::vector<double> > &C) const
{
  //convert from Mole fractions to weight fractions (assumes Containers exist for both)
  for (int i = 0; i < _phases.size(); i++)
  {
    C[i] = _phases[i]->x_to_c(X[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::C_to_X(std::vector<std::vector<double> > &C) const
{
  //convert from weight fractions to Mole fractions and return X "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->c_to_x(C[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > fo_fa_binary::X_to_C(std::vector<std::vector<double> > &X) const
{
  //convert from Mole fractions to weight fractions and return C "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->x_to_c(X[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::set_parameter(const std::string& p, const double& val) const
{
  *parameters[p] = val;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::get_parameter(const std::string& p) const
{
  std::cout << p << " = " << *parameters[p] << std::endl;
}
//-----------------------------------------------------------------------------
void fo_fa_binary::list_parameters() const
{
  std::cout << "Parameters: \n" << std::endl;
  for (auto const& x : parameters)
  {
    std::cout << x.first << " = "  << *x.second << std::endl;
  }
}
//-----------------------------------------------------------------------------
void  fo_fa_binary::report() const
{
  // loop over phases and endmemembers and spit out the names in order
  std::cout << "Reaction object: fo_fa_binary" << std::endl <<std::endl;
  for (int i = 0; i < _phases.size(); i++)
  {
    std::cout << "Phase " << i << " " << _phases[i]->name() << " (" << _phases[i]->abbrev() << ")" << std::endl;
    const std::vector<std::shared_ptr<EndMember> >& endmembers = _phases[i]->endmembers();
    for (int k = 0; k < endmembers.size(); k++)
    {
      std::cout << "     Endmember " << k << " " << endmembers[k]->name() 
                << " : " << endmembers[k]->formula() << "_(" << _phases[i]->abbrev() << ")"
                << std::endl;
    }
  }
  std::cout << std::endl;
  
  for (int j = 0; j < _nu.size(); j++)
  {
    std::cout << "Reaction " << j << std::endl;
    std::vector<std::string> reactants, products;
    for (int i = 0; i < _phases.size(); i++)
    {
      const std::vector<std::shared_ptr<EndMember> >& endmembers = _phases[i]->endmembers();
      for (int k = 0; k < endmembers.size(); k++)
      {
        if (_nu[j][i][k] != 0.0)
        {
          std::ostringstream ss;
          std::string nu_s = _to_string(std::abs(_nu[j][i][k]));
          if (nu_s != "1")
          {
            ss << nu_s << " ";
          }
          ss << endmembers[k]->formula() << "_(" << _phases[i]->abbrev() << ")";

          if (_nu[j][i][k] < 0)
          {
             reactants.push_back(ss.str());
          }
          else if (_nu[j][i][k] > 0)
          {
             products.push_back(ss.str());
          }
        }
      }
    }
    std::string reactants_str = std::accumulate(++reactants.begin(), reactants.end(), reactants[0], 
                                                [](const std::string a, const std::string b){
                                                         return a + " + " + b;
                                                });
    std::string products_str = std::accumulate(++products.begin(), products.end(), products[0],
                                                [](const std::string a, const std::string b){
                                                         return a + " + " + b;
                                                });
    std::cout << "     " << reactants_str << " -> " << products_str <<std::endl;
  }
}
//-----------------------------------------------------------------------------
std::string fo_fa_binary::_to_string(const double& d) const
{
  std::string str = std::to_string(d);
  std::size_t finish = str.find_last_not_of("0");
  if (finish!=std::string::npos)
  {
    str.erase(finish + 1);
  }
  if (str.back()=='.')
  {
    str.pop_back();
  }
  return str;
}
//-----------------------------------------------------------------------------
