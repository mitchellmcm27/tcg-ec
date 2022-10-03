#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <functional>
#include "reactions/eclogite_slb_rx.h"
#include "tcgversion.h"

//-----------------------------------------------------------------------------
eclogite_slb_rx::eclogite_slb_rx()
: _name("eclogite_slb_rx"), J(3),  _phases(std::vector<std::shared_ptr<Phase> > {std::make_shared<Quartz_slb_ph>(),std::make_shared<Albite_slb_ph>(),std::make_shared<NaMajorite_slb_ph>(),std::make_shared<Orthopyroxene_slb_ph>(),std::make_shared<Grossular_slb_ph>(),std::make_shared<Jadeite_slb_ph>() }), _C0({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}), _C({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}),_X({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}),
  _Mu({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}), _dMu_dT({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}), _dMu_dP({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}),_tmp_ik({{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}),
  _A(J,0.), _dA_dT(_A), _dA_dP(_A), _dA_dC(J,{{0.0}, {0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0}}),
  _M({{60.08431}, {262.22304778}, {404.27747555999997}, {200.77763, 263.85743, 202.350107, 216.55053}, {450.44643700000006}, {202.13873777999999}}),
  _nu({{{0.0}, {0.0}, {-1.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {2.0}}, {{0.0}, {0.0}, {-1.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {2.0}}, {{-0.5}, {0.5}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {-0.5}}}),
  _nu_m({{{0.0}, {0.0}, {-1.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {1.0}}, {{0.0}, {0.0}, {-1.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {1.0}}, {{-0.22913435912166485}, {1.0}, {0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {-0.7708656408783351}}}),
  parameters({{"T0", new double}, {"R", new double}}),
  T0(2000.0),
  R(8.31446261815324)
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
eclogite_slb_rx::~eclogite_slb_rx()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string eclogite_slb_rx::name() const
{
  return _name;
}
//-----------------------------------------------------------------------------
std::string eclogite_slb_rx::tcg_build_version()
{
    return TCG_VERSION;
}
//-----------------------------------------------------------------------------
std::string eclogite_slb_rx::tcg_build_git_sha()
{
    return TCG_GIT_SHA;
}
//-----------------------------------------------------------------------------
std::string eclogite_slb_rx::tcg_generation_version()
{
    return "0.6.9";
}
//-----------------------------------------------------------------------------
std::string eclogite_slb_rx::tcg_generation_git_sha()
{
    return "117d758197e2d445579ad671f573064c0650429d Mon Aug 1 00:36:24 2022 +0000";
}
//-----------------------------------------------------------------------------
 std::vector<std::shared_ptr<Phase> > eclogite_slb_rx::phases() const
{
  return _phases;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::zero_C() const
{
  return _C0;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::M() const
{
  return _M;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogite_slb_rx::nu() const
{
  return _nu;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogite_slb_rx::nu_m() const
{
  return _nu_m;
}
//-----------------------------------------------------------------------------
// pass by reference interface
void eclogite_slb_rx::A(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_A) const
{
  //Calculate chemical potential matrix M_ik
  Mu(T, P, C, _Mu);

  _A = {
  1.0*_Mu[2][0] + -0.0*_Mu[3][3] + -2.0*_Mu[5][0] + -0.0*_Mu[4][0] + -0.0*_Mu[0][0],
  1.0*_Mu[2][0] + -0.0*_Mu[3][3] + -2.0*_Mu[5][0] + -0.0*_Mu[1][0],
  0.5*_Mu[5][0] + 0.5*_Mu[0][0] + -0.5*_Mu[1][0]};
}
//-----------------------------------------------------------------------------
// return  vector<double>
std::vector<double> eclogite_slb_rx::A(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  A(T, P, C, _A);
  return _A;
}

//-----------------------------------------------------------------------------
// overloaded version to just return the jth affinity A_j
double eclogite_slb_rx::A(
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
void eclogite_slb_rx::dA_dT(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dT) const
{
  //Calculate chemical potential matrix dMu_dT_ik
  dMu_dT(T, P, C, _dMu_dT);

  _dA_dT = {
  1.0*_dMu_dT[2][0] + -0.0*_dMu_dT[3][3] + -2.0*_dMu_dT[5][0] + -0.0*_dMu_dT[4][0] + -0.0*_dMu_dT[0][0],
  1.0*_dMu_dT[2][0] + -0.0*_dMu_dT[3][3] + -2.0*_dMu_dT[5][0] + -0.0*_dMu_dT[1][0],
  0.5*_dMu_dT[5][0] + 0.5*_dMu_dT[0][0] + -0.5*_dMu_dT[1][0]};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT;
}

//-----------------------------------------------------------------------------
double eclogite_slb_rx::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT[j];
}
//-----------------------------------------------------------------------------
// pass by reference interface
void eclogite_slb_rx::dA_dP(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dP) const
{
  //Calculate chemical potential matrix dMu_dP_ik
  dMu_dP(T, P, C, _dMu_dP);

  _dA_dP = {
  1.0*_dMu_dP[2][0] + -0.0*_dMu_dP[3][3] + -2.0*_dMu_dP[5][0] + -0.0*_dMu_dP[4][0] + -0.0*_dMu_dP[0][0],
  1.0*_dMu_dP[2][0] + -0.0*_dMu_dP[3][3] + -2.0*_dMu_dP[5][0] + -0.0*_dMu_dP[1][0],
  0.5*_dMu_dP[5][0] + 0.5*_dMu_dP[0][0] + -0.5*_dMu_dP[1][0]};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP[j];
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::dA_dC(
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
std::vector<std::vector<std::vector<double> > > eclogite_slb_rx::dA_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  dA_dC( T, P, C, _dA_dC);
  return _dA_dC;
}

//-----------------------------------------------------------------------------
double eclogite_slb_rx::dAj_dCik(
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
void eclogite_slb_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double> &_Gamma) const
{
  //calculate current affinities
  A(T, P, C, _A);
  _Gamma = {
  -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)),
  1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)),
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)),
  0,
  0,
  -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  Gamma_i(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  double Gamma_i;

  switch(i) {
    case 0: Gamma_i = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)); break;
    case 1: Gamma_i = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)); break;
    case 2: Gamma_i = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)); break;
    case 3: Gamma_i = 0; break;
    case 4: Gamma_i = 0; break;
    case 5: Gamma_i = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)); break;
  }
  return Gamma_i;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)),
  1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)),
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)),
  0,
  0,
  -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dT (T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  double dGamma_i_dT;

  switch(i) {
    case 0: dGamma_i_dT = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)); break;
    case 1: dGamma_i_dT = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)); break;
    case 2: dGamma_i_dT = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)); break;
    case 3: dGamma_i_dT = 0; break;
    case 4: dGamma_i_dT = 0; break;
    case 5: dGamma_i_dT = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_i_dT;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)),
  1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)),
  -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)),
  0,
  0,
  -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dP(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  double dGamma_i_dP;

  switch(i) {
    case 0: dGamma_i_dP = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)); break;
    case 1: dGamma_i_dP = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)); break;
    case 2: dGamma_i_dP = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)); break;
    case 3: dGamma_i_dP = 0; break;
    case 4: dGamma_i_dP = 0; break;
    case 5: dGamma_i_dP = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)); break;
  }
  return dGamma_i_dP;
}

//-----------------------------------------------------------------------------
void  eclogite_slb_rx::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<std::vector<double> > >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  _dGamma = {
  {{-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
))},
   {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
))},
   {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
))},
   {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)), -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)), -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)), -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
))},
   {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
))},
   {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
))}},
  {{1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
))},
   {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
))},
   {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
))},
   {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)), 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)), 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)), 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
))},
   {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
))},
   {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
))}},
  {{-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
))},
   {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
))}},
  {{0},
   {0},
   {0},
   {0, 0, 0, 0},
   {0},
   {0}},
  {{0},
   {0},
   {0},
   {0, 0, 0, 0},
   {0},
   {0}},
  {{-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
))},
   {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
))},
   {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
))},
   {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)), -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)), -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)), -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
))},
   {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
))},
   {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
))}}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogite_slb_rx::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dC( T, P, C, Phi, _tmp_dC);
  return _tmp_dC;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_i_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l, unsigned k) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned m = i*24 + l*4 + k;
  double dGamma_i_dC;

  switch(m) {
    case 0: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)); break;
    case 4: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)); break;
    case 8: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)); break;
    case 12: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)); break;
    case 13: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)); break;
    case 14: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)); break;
    case 15: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)); break;
    case 16: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)); break;
    case 20: dGamma_i_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)); break;
    case 24: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)); break;
    case 28: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)); break;
    case 32: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)); break;
    case 36: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)); break;
    case 37: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)); break;
    case 38: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)); break;
    case 39: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)); break;
    case 40: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)); break;
    case 44: dGamma_i_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)); break;
    case 48: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)); break;
    case 52: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)); break;
    case 56: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)); break;
    case 60: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)); break;
    case 61: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)); break;
    case 62: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)); break;
    case 63: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
)); break;
    case 64: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)); break;
    case 68: dGamma_i_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)); break;
    case 72: dGamma_i_dC = 0; break;
    case 76: dGamma_i_dC = 0; break;
    case 80: dGamma_i_dC = 0; break;
    case 84: dGamma_i_dC = 0; break;
    case 85: dGamma_i_dC = 0; break;
    case 86: dGamma_i_dC = 0; break;
    case 87: dGamma_i_dC = 0; break;
    case 88: dGamma_i_dC = 0; break;
    case 92: dGamma_i_dC = 0; break;
    case 96: dGamma_i_dC = 0; break;
    case 100: dGamma_i_dC = 0; break;
    case 104: dGamma_i_dC = 0; break;
    case 108: dGamma_i_dC = 0; break;
    case 109: dGamma_i_dC = 0; break;
    case 110: dGamma_i_dC = 0; break;
    case 111: dGamma_i_dC = 0; break;
    case 112: dGamma_i_dC = 0; break;
    case 116: dGamma_i_dC = 0; break;
    case 120: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)); break;
    case 124: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)); break;
    case 128: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)); break;
    case 132: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)); break;
    case 133: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)); break;
    case 134: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)); break;
    case 135: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
)); break;
    case 136: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)); break;
    case 140: dGamma_i_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)); break;
  }
  return dGamma_i_dC;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
   _dGamma = {
  {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), -0.22913435912166485*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 0, 0, 0, -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))},
  {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), 1.0*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 0, 0, 0, 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))},
  {-1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)), -1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
))},
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0},
  {1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  std::vector<std::vector<double> > _dGamma_dPhi = {
  {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), -0.22913435912166485*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 0, 0, 0, -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))},
  {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), 1.0*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 0, 0, 0, 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))},
  {-1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)), -1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)), -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
))},
  {0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0},
  {1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)), 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)), 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))}};
  return _dGamma_dPhi;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = i*6 + l;
  double dGamma_i_dPhi;

  switch(m) {
    case 0: dGamma_i_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 1: dGamma_i_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 2: dGamma_i_dPhi = 0; break;
    case 3: dGamma_i_dPhi = 0; break;
    case 4: dGamma_i_dPhi = 0; break;
    case 5: dGamma_i_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 6: dGamma_i_dPhi = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 7: dGamma_i_dPhi = 1.0*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 8: dGamma_i_dPhi = 0; break;
    case 9: dGamma_i_dPhi = 0; break;
    case 10: dGamma_i_dPhi = 0; break;
    case 11: dGamma_i_dPhi = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 12: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)); break;
    case 13: dGamma_i_dPhi = -1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)); break;
    case 14: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 15: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 16: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)); break;
    case 17: dGamma_i_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)); break;
    case 18: dGamma_i_dPhi = 0; break;
    case 19: dGamma_i_dPhi = 0; break;
    case 20: dGamma_i_dPhi = 0; break;
    case 21: dGamma_i_dPhi = 0; break;
    case 22: dGamma_i_dPhi = 0; break;
    case 23: dGamma_i_dPhi = 0; break;
    case 24: dGamma_i_dPhi = 0; break;
    case 25: dGamma_i_dPhi = 0; break;
    case 26: dGamma_i_dPhi = 0; break;
    case 27: dGamma_i_dPhi = 0; break;
    case 28: dGamma_i_dPhi = 0; break;
    case 29: dGamma_i_dPhi = 0; break;
    case 30: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 31: dGamma_i_dPhi = 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 32: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 33: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 34: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)); break;
    case 35: dGamma_i_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
  }
  return dGamma_i_dPhi;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<double> > &_Gamma) const
{
  A(T, P, C, _A);
  _Gamma = {
  {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
))},
  {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
))},
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
))},
  {0, 0, 0, 0},
  {0},
  {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  Gamma_ik(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  unsigned m = i*4 + k;
  double Gamma_ik;

  switch(m) {
    case 0: Gamma_ik = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)); break;
    case 4: Gamma_ik = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)); break;
    case 8: Gamma_ik = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)); break;
    case 12: Gamma_ik = 0; break;
    case 13: Gamma_ik = 0; break;
    case 14: Gamma_ik = 0; break;
    case 15: Gamma_ik = 0; break;
    case 16: Gamma_ik = 0; break;
    case 20: Gamma_ik = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)); break;
  }
  return Gamma_ik;
}
//-----------------------------------------------------------------------------
void  eclogite_slb_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
))},
  {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
))},
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
))},
  {0, 0, 0, 0},
  {0},
  {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: using_tmp_ik but starting to worry about thread safety
  dGamma_ik_dT(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  unsigned m = i*4 + k;
  double dGamma_ik_dT;

  switch(m) {
    case 0: dGamma_ik_dT = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)); break;
    case 4: dGamma_ik_dT = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)); break;
    case 8: dGamma_ik_dT = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)); break;
    case 12: dGamma_ik_dT = 0; break;
    case 13: dGamma_ik_dT = 0; break;
    case 14: dGamma_ik_dT = 0; break;
    case 15: dGamma_ik_dT = 0; break;
    case 16: dGamma_ik_dT = 0; break;
    case 20: dGamma_ik_dT = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[2]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_ik_dT;
}
//-----------------------------------------------------------------------------
void  eclogite_slb_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
))},
  {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
))},
  {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
))},
  {0, 0, 0, 0},
  {0},
  {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> >  eclogite_slb_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: thread safety again?
  dGamma_ik_dP(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  unsigned m = i*4 + k;
  double dGamma_ik_dP;

  switch(m) {
    case 0: dGamma_ik_dP = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)); break;
    case 4: dGamma_ik_dP = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)); break;
    case 8: dGamma_ik_dP = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)); break;
    case 12: dGamma_ik_dP = 0; break;
    case 13: dGamma_ik_dP = 0; break;
    case 14: dGamma_ik_dP = 0; break;
    case 15: dGamma_ik_dP = 0; break;
    case 16: dGamma_ik_dP = 0; break;
    case 20: dGamma_ik_dP = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dP[2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dP[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)); break;
  }
  return dGamma_ik_dP;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);

    switch(i) {
    case 0: _dGamma = {
 {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
))}
    }; break;

    case 1: _dGamma = {
 {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
))}
    }; break;

    case 2: _dGamma = {
 {-1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
))}
    }; break;

    case 3: _dGamma = {
 {0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0}
    }; break;

    case 4: _dGamma = {
 {0}
    }; break;

    case 5: _dGamma = {
 {-0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::dGamma_ik_dC(
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
void eclogite_slb_rx::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
    switch(i) {
    case 0: _dGamma = {
 {-0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)),
 -0.22913435912166485*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)),
 0,
 0,
 0,
 -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))}
    }; break;

    case 1: _dGamma = {
 {1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)),
 1.0*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)),
 0,
 0,
 0,
 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))}
    }; break;

    case 2: _dGamma = {
 {-1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)),
 -1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)),
 -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)),
 -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)),
 -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)),
 -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
))}
    }; break;

    case 3: _dGamma = {
 {0,
 0,
 0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0,
 0,
 0}
    }; break;

    case 4: _dGamma = {
 {0,
 0,
 0,
 0,
 0,
 0}
    }; break;

    case 5: _dGamma = {
 {1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)),
 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)),
 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)),
 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)),
 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)),
 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::dGamma_ik_dPhi(
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
double eclogite_slb_rx::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l, unsigned m) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned p = l*96 + m*24 + i*4+ k;
  double dGamma_ik_dC;

  switch(p) {
    case 0: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)); break;
    case 4: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)); break;
    case 8: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)); break;
    case 12: dGamma_ik_dC = 0; break;
    case 13: dGamma_ik_dC = 0; break;
    case 14: dGamma_ik_dC = 0; break;
    case 15: dGamma_ik_dC = 0; break;
    case 16: dGamma_ik_dC = 0; break;
    case 20: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][0][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)); break;
    case 96: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)); break;
    case 100: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)); break;
    case 104: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)); break;
    case 108: dGamma_ik_dC = 0; break;
    case 109: dGamma_ik_dC = 0; break;
    case 110: dGamma_ik_dC = 0; break;
    case 111: dGamma_ik_dC = 0; break;
    case 112: dGamma_ik_dC = 0; break;
    case 116: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][1][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)); break;
    case 192: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)); break;
    case 196: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)); break;
    case 200: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)); break;
    case 204: dGamma_ik_dC = 0; break;
    case 205: dGamma_ik_dC = 0; break;
    case 206: dGamma_ik_dC = 0; break;
    case 207: dGamma_ik_dC = 0; break;
    case 208: dGamma_ik_dC = 0; break;
    case 212: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][2][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)); break;
    case 288: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)); break;
    case 292: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)); break;
    case 296: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)); break;
    case 300: dGamma_ik_dC = 0; break;
    case 301: dGamma_ik_dC = 0; break;
    case 302: dGamma_ik_dC = 0; break;
    case 303: dGamma_ik_dC = 0; break;
    case 304: dGamma_ik_dC = 0; break;
    case 308: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)); break;
    case 312: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)); break;
    case 316: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)); break;
    case 320: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)); break;
    case 324: dGamma_ik_dC = 0; break;
    case 325: dGamma_ik_dC = 0; break;
    case 326: dGamma_ik_dC = 0; break;
    case 327: dGamma_ik_dC = 0; break;
    case 328: dGamma_ik_dC = 0; break;
    case 332: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][1]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)); break;
    case 336: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)); break;
    case 340: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)); break;
    case 344: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)); break;
    case 348: dGamma_ik_dC = 0; break;
    case 349: dGamma_ik_dC = 0; break;
    case 350: dGamma_ik_dC = 0; break;
    case 351: dGamma_ik_dC = 0; break;
    case 352: dGamma_ik_dC = 0; break;
    case 356: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][2]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][2]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][2]/(R*T)
)); break;
    case 360: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)); break;
    case 364: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)); break;
    case 368: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
)); break;
    case 372: dGamma_ik_dC = 0; break;
    case 373: dGamma_ik_dC = 0; break;
    case 374: dGamma_ik_dC = 0; break;
    case 375: dGamma_ik_dC = 0; break;
    case 376: dGamma_ik_dC = 0; break;
    case 380: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][3][3]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][3][3]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][3][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][3]/(R*T)
)); break;
    case 384: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)); break;
    case 388: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)); break;
    case 392: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)); break;
    case 396: dGamma_ik_dC = 0; break;
    case 397: dGamma_ik_dC = 0; break;
    case 398: dGamma_ik_dC = 0; break;
    case 399: dGamma_ik_dC = 0; break;
    case 400: dGamma_ik_dC = 0; break;
    case 404: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][4][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)); break;
    case 480: dGamma_ik_dC = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)); break;
    case 484: dGamma_ik_dC = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)); break;
    case 488: dGamma_ik_dC = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)); break;
    case 492: dGamma_ik_dC = 0; break;
    case 493: dGamma_ik_dC = 0; break;
    case 494: dGamma_ik_dC = 0; break;
    case 495: dGamma_ik_dC = 0; break;
    case 496: dGamma_ik_dC = 0; break;
    case 500: dGamma_ik_dC = -0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*_dA_dC[2][5][0]/(R*T)
)) + 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)); break;
  }
  return dGamma_ik_dC;
}
//-----------------------------------------------------------------------------
double eclogite_slb_rx::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = l*24 + i*4 + k;
  double dGamma_ik_dPhi;

  switch(m) {
    case 0: dGamma_ik_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 4: dGamma_ik_dPhi = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 8: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)); break;
    case 12: dGamma_ik_dPhi = 0; break;
    case 13: dGamma_ik_dPhi = 0; break;
    case 14: dGamma_ik_dPhi = 0; break;
    case 15: dGamma_ik_dPhi = 0; break;
    case 16: dGamma_ik_dPhi = 0; break;
    case 20: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[0]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 24: dGamma_ik_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 28: dGamma_ik_dPhi = 1.0*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 32: dGamma_ik_dPhi = -1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)); break;
    case 36: dGamma_ik_dPhi = 0; break;
    case 37: dGamma_ik_dPhi = 0; break;
    case 38: dGamma_ik_dPhi = 0; break;
    case 39: dGamma_ik_dPhi = 0; break;
    case 40: dGamma_ik_dPhi = 0; break;
    case 44: dGamma_ik_dPhi = 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[2]/(R*T)
)); break;
    case 48: dGamma_ik_dPhi = 0; break;
    case 52: dGamma_ik_dPhi = 0; break;
    case 56: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 60: dGamma_ik_dPhi = 0; break;
    case 61: dGamma_ik_dPhi = 0; break;
    case 62: dGamma_ik_dPhi = 0; break;
    case 63: dGamma_ik_dPhi = 0; break;
    case 64: dGamma_ik_dPhi = 0; break;
    case 68: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 72: dGamma_ik_dPhi = 0; break;
    case 76: dGamma_ik_dPhi = 0; break;
    case 80: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) - 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 84: dGamma_ik_dPhi = 0; break;
    case 85: dGamma_ik_dPhi = 0; break;
    case 86: dGamma_ik_dPhi = 0; break;
    case 87: dGamma_ik_dPhi = 0; break;
    case 88: dGamma_ik_dPhi = 0; break;
    case 92: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)
: (
   0
)) + 1.0*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[2]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 96: dGamma_ik_dPhi = 0; break;
    case 100: dGamma_ik_dPhi = 0; break;
    case 104: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)); break;
    case 108: dGamma_ik_dPhi = 0; break;
    case 109: dGamma_ik_dPhi = 0; break;
    case 110: dGamma_ik_dPhi = 0; break;
    case 111: dGamma_ik_dPhi = 0; break;
    case 112: dGamma_ik_dPhi = 0; break;
    case 116: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[0]/(R*T)
)); break;
    case 120: dGamma_ik_dPhi = -0.22913435912166485*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 124: dGamma_ik_dPhi = 1.0*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 128: dGamma_ik_dPhi = -1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) - 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)); break;
    case 132: dGamma_ik_dPhi = 0; break;
    case 133: dGamma_ik_dPhi = 0; break;
    case 134: dGamma_ik_dPhi = 0; break;
    case 135: dGamma_ik_dPhi = 0; break;
    case 136: dGamma_ik_dPhi = 0; break;
    case 140: dGamma_ik_dPhi = 1.0*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[4]*_A[0]/(R*T)
)) + 1.0*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.77086564087833509*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)
: (
   0
)); break;
  }
  return dGamma_ik_dPhi;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::rho(const double &T, const double &P,
		std::vector<std::vector<double> >& C,
        std::vector<double>&  _rho) const
{
  // inplace calculation of vector of phase densities (must be correct size)
  for (int i = 0; i < _phases.size(); i++ ) {
    _rho[i] = _phases[i]->rho(T,P,C[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> eclogite_slb_rx::rho(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _rho(6,0.);
  rho(T, P, C, _rho);
  return _rho;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::Cp(
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
std::vector<double> eclogite_slb_rx::Cp(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{

  // return vector of densities
  std::vector<double> _Cp(6,0.);
  Cp(T, P, C, _Cp);
  return _Cp;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::s(
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
std::vector<double> eclogite_slb_rx::s(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _s(6,0.);
  s(T, P, C, _s);
  return _s;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::ds_dC(
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
std::vector<std::vector<double> >eclogite_slb_rx::ds_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  ds_dC(T, P, C, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::alpha(
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
std::vector<double> eclogite_slb_rx::alpha(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _alpha(6,0.);
  alpha(T, P, C, _alpha);
  return _alpha;

}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::beta(
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
std::vector<double> eclogite_slb_rx::beta(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _beta(6,0.);
  beta(T, P, C, _beta);
  return _beta;

}

//-----------------------------------------------------------------------------
void eclogite_slb_rx::Mu(
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
std::vector<std::vector<double> > eclogite_slb_rx::Mu(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  Mu(T,P,C,_Mu);
  return _Mu;
}

//-----------------------------------------------------------------------------
void eclogite_slb_rx::dMu_dT(
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
void eclogite_slb_rx::dMu_dP(
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
void eclogite_slb_rx::dMu_dC(
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
std::vector<std::vector<std::vector<double> > > eclogite_slb_rx::dMu_dC(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  dMu_dC(T,P,C,_dMu_dC);
  return _dMu_dC;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::C_to_X(
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
void eclogite_slb_rx::X_to_C(
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
std::vector<std::vector<double> > eclogite_slb_rx::C_to_X(std::vector<std::vector<double> > &C) const
{
  //convert from weight fractions to Mole fractions and return X "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->c_to_x(C[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogite_slb_rx::X_to_C(std::vector<std::vector<double> > &X) const
{
  //convert from Mole fractions to weight fractions and return C "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->x_to_c(X[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::set_parameter(const std::string& p, const double& val) const
{
  *parameters[p] = val;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::get_parameter(const std::string& p) const
{
  std::cout << p << " = " << *parameters[p] << std::endl;
}
//-----------------------------------------------------------------------------
void eclogite_slb_rx::list_parameters() const
{
  std::cout << "Parameters: \n" << std::endl;
  for (auto const& x : parameters)
  {
    std::cout << x.first << " = "  << *x.second << std::endl;
  }
}
//-----------------------------------------------------------------------------
void  eclogite_slb_rx::report() const
{
  // loop over phases and endmemembers and spit out the names in order
  std::cout << "Reaction object: eclogite_slb_rx" << std::endl <<std::endl;
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
std::string eclogite_slb_rx::_to_string(const double& d) const
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
