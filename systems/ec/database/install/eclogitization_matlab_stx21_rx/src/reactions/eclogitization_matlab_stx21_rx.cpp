#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <functional>
#include "reactions/eclogitization_matlab_stx21_rx.h"
#include "tcgversion.h"

//-----------------------------------------------------------------------------
eclogitization_matlab_stx21_rx::eclogitization_matlab_stx21_rx()
: _name("eclogitization_matlab_stx21_rx"), J(10),  _phases(std::vector<std::shared_ptr<Phase> > {std::make_shared<Clinopyroxene_slb_ph>(),std::make_shared<Orthopyroxene_slb_ph>(),std::make_shared<Quartz_slb_ph>(),std::make_shared<Feldspar_slb_ph>(),std::make_shared<Garnet_slb_ph>(),std::make_shared<Kyanite_slb_ph>() }), _C0({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}), _C({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}),_X({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}),
  _Mu({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}), _dMu_dT({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}), _dMu_dP({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}),_tmp_ik({{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}),
  _A(J,0.), _dA_dT(_A), _dA_dP(_A), _dA_dC(J,{{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}),
  _M({{216.55053, 248.09043, 200.77763, 218.123007, 202.13873777999999}, {200.77763, 263.85743, 202.350107, 216.55053}, {60.08431}, {278.207317, 262.22304778}, {403.127737, 497.747437, 450.44643700000006, 401.55526, 404.27747555999997}, {162.04560199999997}}),
  _nu({{{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.3333333333333333, 0.0, 0.6666666666666666, 0.0}, {0.6666666666666666}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {-0.6666666666666666}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.3333333333333333, 0.0, 0.0, 0.6666666666666666}, {0.0}, {-0.6666666666666666, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.6666666666666666}}, {{0.6666666666666666, 0.0, 0.0, 0.0, 0.0}, {-0.3333333333333333, 0.0, 0.0, 0.0}, {0.0}, {-0.6666666666666666, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.6666666666666666}}, {{0.0, 0.6666666666666666, 0.0, 0.0, 0.0}, {0.0, -0.3333333333333333, 0.0, 0.0}, {0.0}, {-0.6666666666666666, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.6666666666666666}}, {{0.0, 0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {1.0}, {-1.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 0.0, 0.0}, {1.0}, {0.0, -1.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.3333333333333333}, {-1.0, 0.0}, {0.0, 0.0, 0.3333333333333333, 0.0, 0.0}, {0.6666666666666666}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.6, 0.0, 0.0, 0.0}, {0.4}, {0.0, 0.0}, {0.4, 0.0, 0.0, 0.0, 0.0}, {-0.4}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, -0.6, 0.0, 0.0}, {0.4}, {0.0, 0.0}, {0.0, 0.4, 0.0, 0.0, 0.0}, {-0.4}}}),
  _nu_m({{{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.38252915203572546, 0.0, 0.771050189655574, 0.0}, {0.228949810344426}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {-0.6174708479642743}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.2651606990004853, 0.0, 0.0, 0.5719829435552712}, {0.0}, {-0.7348393009995148, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.42801705644472876}}, {{0.5719829435552712, 0.0, 0.0, 0.0, 0.0}, {-0.2651606990004853, 0.0, 0.0, 0.0}, {0.0}, {-0.7348393009995148, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.42801705644472876}}, {{0.0, 0.6048979134805693, 0.0, 0.0, 0.0}, {0.0, -0.3216706280515242, 0.0, 0.0}, {0.0}, {-0.6783293719484759, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.3951020865194307}}, {{0.0, 0.0, 0.0, 0.7840304466183396, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.21596955338166035}, {-1.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0, 0.0}, {0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 0.0, 0.0, 0.7708656408783351}, {0.0, 0.0, 0.0, 0.0}, {0.22913435912166485}, {0.0, -1.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.07198985112722012}, {-1.0, 0.0}, {0.0, 0.0, 0.5397011622571141, 0.0, 0.0}, {0.3883089866156659}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {-0.6501697159011927, 0.0, 0.0, 0.0}, {0.129712321579581}, {0.0, 0.0}, {0.870287678420419, 0.0, 0.0, 0.0, 0.0}, {-0.34983028409880706}}, {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, -0.7095081037042519, 0.0, 0.0}, {0.10771045270035519}, {0.0, 0.0}, {0.0, 0.8922895472996448, 0.0, 0.0, 0.0}, {-0.29049189629574806}}}),
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
eclogitization_matlab_stx21_rx::~eclogitization_matlab_stx21_rx()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string eclogitization_matlab_stx21_rx::name() const
{
  return _name;
}
//-----------------------------------------------------------------------------
std::string eclogitization_matlab_stx21_rx::tcg_build_version()
{
    return TCG_VERSION;
}
//-----------------------------------------------------------------------------
std::string eclogitization_matlab_stx21_rx::tcg_build_git_sha()
{
    return TCG_GIT_SHA;
}
//-----------------------------------------------------------------------------
std::string eclogitization_matlab_stx21_rx::tcg_generation_version()
{
    return "0.6.9";
}
//-----------------------------------------------------------------------------
std::string eclogitization_matlab_stx21_rx::tcg_generation_git_sha()
{
    return "117d758197e2d445579ad671f573064c0650429d Mon Aug 1 00:36:24 2022 +0000";
}
//-----------------------------------------------------------------------------
 std::vector<std::shared_ptr<Phase> > eclogitization_matlab_stx21_rx::phases() const
{
  return _phases;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::zero_C() const
{
  return _C0;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::M() const
{
  return _M;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogitization_matlab_stx21_rx::nu() const
{
  return _nu;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogitization_matlab_stx21_rx::nu_m() const
{
  return _nu_m;
}
//-----------------------------------------------------------------------------
// pass by reference interface
void eclogitization_matlab_stx21_rx::A(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_A) const
{
  //Calculate chemical potential matrix M_ik
  Mu(T, P, C, _Mu);

  _A = {
  0.6666666666666666*_Mu[5][0] + 0.3333333333333333*_Mu[1][0] + -0.6666666666666666*_Mu[2][0] + -0.6666666666666666*_Mu[1][2],
  0.6666666666666666*_Mu[3][0] + 0.3333333333333333*_Mu[1][0] + -0.6666666666666666*_Mu[5][0] + -0.6666666666666666*_Mu[1][3],
  0.6666666666666666*_Mu[3][0] + 0.3333333333333333*_Mu[1][0] + -0.6666666666666666*_Mu[5][0] + -0.6666666666666666*_Mu[0][0],
  0.6666666666666666*_Mu[3][0] + 0.3333333333333333*_Mu[1][1] + -0.6666666666666666*_Mu[5][0] + -0.6666666666666666*_Mu[0][1],
  1.0*_Mu[3][0] + -1.0*_Mu[2][0] + -1.0*_Mu[0][3],
  1.0*_Mu[1][0] + -1.0*_Mu[0][2],
  1.0*_Mu[3][1] + -1.0*_Mu[2][0] + -1.0*_Mu[0][4],
  1.0*_Mu[3][0] + -0.3333333333333333*_Mu[2][0] + -0.6666666666666666*_Mu[5][0] + -0.3333333333333333*_Mu[4][2],
  0.4*_Mu[5][0] + 0.6*_Mu[1][0] + -0.4*_Mu[2][0] + -0.4*_Mu[4][0],
  0.4*_Mu[5][0] + 0.6*_Mu[1][1] + -0.4*_Mu[2][0] + -0.4*_Mu[4][1]};
}
//-----------------------------------------------------------------------------
// return  vector<double>
std::vector<double> eclogitization_matlab_stx21_rx::A(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  A(T, P, C, _A);
  return _A;
}

//-----------------------------------------------------------------------------
// overloaded version to just return the jth affinity A_j
double eclogitization_matlab_stx21_rx::A(
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
void eclogitization_matlab_stx21_rx::dA_dT(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dT) const
{
  //Calculate chemical potential matrix dMu_dT_ik
  dMu_dT(T, P, C, _dMu_dT);

  _dA_dT = {
  0.6666666666666666*_dMu_dT[5][0] + 0.3333333333333333*_dMu_dT[1][0] + -0.6666666666666666*_dMu_dT[2][0] + -0.6666666666666666*_dMu_dT[1][2],
  0.6666666666666666*_dMu_dT[3][0] + 0.3333333333333333*_dMu_dT[1][0] + -0.6666666666666666*_dMu_dT[5][0] + -0.6666666666666666*_dMu_dT[1][3],
  0.6666666666666666*_dMu_dT[3][0] + 0.3333333333333333*_dMu_dT[1][0] + -0.6666666666666666*_dMu_dT[5][0] + -0.6666666666666666*_dMu_dT[0][0],
  0.6666666666666666*_dMu_dT[3][0] + 0.3333333333333333*_dMu_dT[1][1] + -0.6666666666666666*_dMu_dT[5][0] + -0.6666666666666666*_dMu_dT[0][1],
  1.0*_dMu_dT[3][0] + -1.0*_dMu_dT[2][0] + -1.0*_dMu_dT[0][3],
  1.0*_dMu_dT[1][0] + -1.0*_dMu_dT[0][2],
  1.0*_dMu_dT[3][1] + -1.0*_dMu_dT[2][0] + -1.0*_dMu_dT[0][4],
  1.0*_dMu_dT[3][0] + -0.3333333333333333*_dMu_dT[2][0] + -0.6666666666666666*_dMu_dT[5][0] + -0.3333333333333333*_dMu_dT[4][2],
  0.4*_dMu_dT[5][0] + 0.6*_dMu_dT[1][0] + -0.4*_dMu_dT[2][0] + -0.4*_dMu_dT[4][0],
  0.4*_dMu_dT[5][0] + 0.6*_dMu_dT[1][1] + -0.4*_dMu_dT[2][0] + -0.4*_dMu_dT[4][1]};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT;
}

//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dA_dT(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dT(T, P, C, _dA_dT);
  return _dA_dT[j];
}
//-----------------------------------------------------------------------------
// pass by reference interface
void eclogitization_matlab_stx21_rx::dA_dP(
    const double &T, const double &P,
		std::vector<std::vector<double> > &C, std::vector<double> &_dA_dP) const
{
  //Calculate chemical potential matrix dMu_dP_ik
  dMu_dP(T, P, C, _dMu_dP);

  _dA_dP = {
  0.6666666666666666*_dMu_dP[5][0] + 0.3333333333333333*_dMu_dP[1][0] + -0.6666666666666666*_dMu_dP[2][0] + -0.6666666666666666*_dMu_dP[1][2],
  0.6666666666666666*_dMu_dP[3][0] + 0.3333333333333333*_dMu_dP[1][0] + -0.6666666666666666*_dMu_dP[5][0] + -0.6666666666666666*_dMu_dP[1][3],
  0.6666666666666666*_dMu_dP[3][0] + 0.3333333333333333*_dMu_dP[1][0] + -0.6666666666666666*_dMu_dP[5][0] + -0.6666666666666666*_dMu_dP[0][0],
  0.6666666666666666*_dMu_dP[3][0] + 0.3333333333333333*_dMu_dP[1][1] + -0.6666666666666666*_dMu_dP[5][0] + -0.6666666666666666*_dMu_dP[0][1],
  1.0*_dMu_dP[3][0] + -1.0*_dMu_dP[2][0] + -1.0*_dMu_dP[0][3],
  1.0*_dMu_dP[1][0] + -1.0*_dMu_dP[0][2],
  1.0*_dMu_dP[3][1] + -1.0*_dMu_dP[2][0] + -1.0*_dMu_dP[0][4],
  1.0*_dMu_dP[3][0] + -0.3333333333333333*_dMu_dP[2][0] + -0.6666666666666666*_dMu_dP[5][0] + -0.3333333333333333*_dMu_dP[4][2],
  0.4*_dMu_dP[5][0] + 0.6*_dMu_dP[1][0] + -0.4*_dMu_dP[2][0] + -0.4*_dMu_dP[4][0],
  0.4*_dMu_dP[5][0] + 0.6*_dMu_dP[1][1] + -0.4*_dMu_dP[2][0] + -0.4*_dMu_dP[4][1]};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C) const
{
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dA_dP(
    const double &T, const double &P,
    std::vector<std::vector<double> > &C, const int j) const
{
  //FIXME:  now redundant to above call but just returns one scalar...not useful
  dA_dP(T, P, C, _dA_dP);
  return _dA_dP[j];
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dA_dC(
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
std::vector<std::vector<std::vector<double> > > eclogitization_matlab_stx21_rx::dA_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  dA_dC( T, P, C, _dA_dC);
  return _dA_dC;
}

//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dAj_dCik(
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
void eclogitization_matlab_stx21_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double> &_Gamma) const
{
  //calculate current affinities
  A(T, P, C, _A);
  _Gamma = {
  1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)),
  -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)),
  0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)),
  -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)),
  0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)),
  0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  Gamma_i(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::Gamma_i(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  double Gamma_i;

  switch(i) {
    case 0: Gamma_i = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)); break;
    case 1: Gamma_i = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 2: Gamma_i = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 3: Gamma_i = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)); break;
    case 4: Gamma_i = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 5: Gamma_i = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
  }
  return Gamma_i;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)),
  -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)),
  0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)),
  -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)),
  0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)),
  0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dT (T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_i_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  double dGamma_i_dT;

  switch(i) {
    case 0: dGamma_i_dT = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)); break;
    case 1: dGamma_i_dT = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 2: dGamma_i_dT = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 3: dGamma_i_dT = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)); break;
    case 4: dGamma_i_dT = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 5: dGamma_i_dT = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_i_dT;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<double>& _dGamma) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)),
  -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)),
  0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)),
  -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)),
  0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)),
  0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
))};
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dP(T, P, C, Phi, _tmp);
  return _tmp;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_i_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, const int i) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  double dGamma_i_dP;

  switch(i) {
    case 0: dGamma_i_dP = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)); break;
    case 1: dGamma_i_dP = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 2: dGamma_i_dP = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 3: dGamma_i_dP = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)); break;
    case 4: dGamma_i_dP = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 5: dGamma_i_dP = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
  }
  return dGamma_i_dP;
}

//-----------------------------------------------------------------------------
void  eclogitization_matlab_stx21_rx::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<std::vector<double> > >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  _dGamma = {
  {{1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
))},
   {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
))},
   {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
))},
   {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
))},
   {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
))},
   {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
))}},
  {{-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
))},
   {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
))},
   {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
))},
   {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
))},
   {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
))},
   {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
))}},
  {{0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
))},
   {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
))},
   {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
))},
   {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
))},
   {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
))},
   {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
))}},
  {{-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
))},
   {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
))},
   {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
))},
   {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
))},
   {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
))},
   {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
))}},
  {{0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
))},
   {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
))},
   {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
))},
   {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
))},
   {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
))},
   {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
))}},
  {{0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
))},
   {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
))},
   {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
))},
   {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
))},
   {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
))},
   {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
))}}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<double> > > eclogitization_matlab_stx21_rx::dGamma_i_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  dGamma_i_dC( T, P, C, Phi, _tmp_dC);
  return _tmp_dC;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_i_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l, unsigned k) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned m = i*30 + l*5 + k;
  double dGamma_i_dC;

  switch(m) {
    case 0: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)); break;
    case 1: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)); break;
    case 2: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)); break;
    case 3: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)); break;
    case 4: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)); break;
    case 5: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)); break;
    case 6: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)); break;
    case 7: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)); break;
    case 8: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)); break;
    case 10: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)); break;
    case 15: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)); break;
    case 16: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)); break;
    case 20: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)); break;
    case 21: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)); break;
    case 22: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)); break;
    case 23: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)); break;
    case 24: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)); break;
    case 25: dGamma_i_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)); break;
    case 30: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 31: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 32: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 33: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 34: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 35: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 36: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 37: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 38: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 40: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 45: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 46: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 50: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 51: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 52: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 53: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 54: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 55: dGamma_i_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 60: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 61: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 62: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 63: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 64: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 65: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 66: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 67: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 68: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 70: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 75: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 76: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 80: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 81: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 82: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 83: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 84: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 85: dGamma_i_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 90: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)); break;
    case 91: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)); break;
    case 92: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)); break;
    case 93: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)); break;
    case 94: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)); break;
    case 95: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)); break;
    case 96: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)); break;
    case 97: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)); break;
    case 98: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)); break;
    case 100: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)); break;
    case 105: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)); break;
    case 106: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)); break;
    case 110: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)); break;
    case 111: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)); break;
    case 112: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)); break;
    case 113: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)); break;
    case 114: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)); break;
    case 115: dGamma_i_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)); break;
    case 120: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 121: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 122: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 123: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 124: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 125: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 126: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 127: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 128: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 130: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 135: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 136: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 140: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 141: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 142: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 143: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 144: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 145: dGamma_i_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 150: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 151: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 152: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 153: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 154: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 155: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 156: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 157: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 158: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 160: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 165: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 166: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 170: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 171: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 172: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 173: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 174: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 175: dGamma_i_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
  }
  return dGamma_i_dC;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
   _dGamma = {
  {0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)), 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)), 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), 0, 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
))},
  {-0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.38852103761984857*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), -0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.30682224455478596*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.22894981034442599*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)), 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {-0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), -0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)), -1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), -1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)), -0.73483930099951478*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
))},
  {0, 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)), 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), -0.61747084796427432*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.42801705644472876*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  A(T, P, C, _A);
  std::vector<std::vector<double> > _dGamma_dPhi = {
  {0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)), 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)), 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), 0, 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
))},
  {-0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)), -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.38852103761984857*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), -0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.30682224455478596*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.22894981034442599*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)), 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {-0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)), -0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)), -1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)), -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), -1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)), -0.73483930099951478*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
))},
  {0, 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
  {0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)), 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)), -0.61747084796427432*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)), 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)), 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)), 0.42801705644472876*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))}};
  return _dGamma_dPhi;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_i_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = i*6 + l;
  double dGamma_i_dPhi;

  switch(m) {
    case 0: dGamma_i_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) + 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 1: dGamma_i_dPhi = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 2: dGamma_i_dPhi = 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)); break;
    case 3: dGamma_i_dPhi = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 4: dGamma_i_dPhi = 0; break;
    case 5: dGamma_i_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)); break;
    case 6: dGamma_i_dPhi = -0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)); break;
    case 7: dGamma_i_dPhi = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) + 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 8: dGamma_i_dPhi = 0.38852103761984857*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 9: dGamma_i_dPhi = 0.30682224455478596*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 10: dGamma_i_dPhi = -0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 11: dGamma_i_dPhi = 0.30682224455478596*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38852103761984857*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 12: dGamma_i_dPhi = 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 13: dGamma_i_dPhi = 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 14: dGamma_i_dPhi = 0.22894981034442599*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 15: dGamma_i_dPhi = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)); break;
    case 16: dGamma_i_dPhi = 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 17: dGamma_i_dPhi = 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 18: dGamma_i_dPhi = -0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 19: dGamma_i_dPhi = -0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 20: dGamma_i_dPhi = -1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) - 1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)); break;
    case 21: dGamma_i_dPhi = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) - 1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 22: dGamma_i_dPhi = -1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)); break;
    case 23: dGamma_i_dPhi = -0.73483930099951478*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)); break;
    case 24: dGamma_i_dPhi = 0; break;
    case 25: dGamma_i_dPhi = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 26: dGamma_i_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 27: dGamma_i_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)); break;
    case 28: dGamma_i_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 29: dGamma_i_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 30: dGamma_i_dPhi = 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)); break;
    case 31: dGamma_i_dPhi = 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 32: dGamma_i_dPhi = -0.61747084796427432*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 33: dGamma_i_dPhi = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 34: dGamma_i_dPhi = 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 35: dGamma_i_dPhi = 0.42801705644472876*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
  }
  return dGamma_i_dPhi;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    std::vector<std::vector<double> > &_Gamma) const
{
  A(T, P, C, _A);
  _Gamma = {
  {0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)), 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)), 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)), 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
))},
  {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)), -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)), 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)), 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
))},
  {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
))},
  {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)), -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
))},
  {0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)), 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)), 0, 0},
  {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  Gamma_ik(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::Gamma_ik(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  unsigned m = i*5 + k;
  double Gamma_ik;

  switch(m) {
    case 0: Gamma_ik = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)); break;
    case 1: Gamma_ik = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)); break;
    case 2: Gamma_ik = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)); break;
    case 3: Gamma_ik = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)); break;
    case 4: Gamma_ik = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)); break;
    case 5: Gamma_ik = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_A[5]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)); break;
    case 6: Gamma_ik = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 7: Gamma_ik = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)); break;
    case 8: Gamma_ik = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)); break;
    case 10: Gamma_ik = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 15: Gamma_ik = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)); break;
    case 16: Gamma_ik = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*T)
)); break;
    case 20: Gamma_ik = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)); break;
    case 21: Gamma_ik = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
    case 22: Gamma_ik = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)); break;
    case 23: Gamma_ik = 0; break;
    case 24: Gamma_ik = 0; break;
    case 25: Gamma_ik = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*T)
)); break;
  }
  return Gamma_ik;
}
//-----------------------------------------------------------------------------
void  eclogitization_matlab_stx21_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  _dGamma = {
  {0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)), 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)), 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)), 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
))},
  {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)), -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)), 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)), 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
))},
  {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
))},
  {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)), -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
))},
  {0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)), 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)), 0, 0},
  {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: using_tmp_ik but starting to worry about thread safety
  dGamma_ik_dT(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_ik_dT(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dT(T, P, C, _dA_dT);
  unsigned m = i*5 + k;
  double dGamma_ik_dT;

  switch(m) {
    case 0: dGamma_ik_dT = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)); break;
    case 1: dGamma_ik_dT = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)); break;
    case 2: dGamma_ik_dT = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)); break;
    case 3: dGamma_ik_dT = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)); break;
    case 4: dGamma_ik_dT = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)); break;
    case 5: dGamma_ik_dT = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*_A[5]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dT[5]/(R*T) - std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*_A[5]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)); break;
    case 6: dGamma_ik_dT = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 7: dGamma_ik_dT = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)); break;
    case 8: dGamma_ik_dT = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)); break;
    case 10: dGamma_ik_dT = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 15: dGamma_ik_dT = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[4]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[4]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[4]/(R*((T)*(T)*(T)))
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)); break;
    case 16: dGamma_ik_dT = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[6]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dT[6]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[2]*_A[6]/(R*((T)*(T)*(T)))
)); break;
    case 20: dGamma_ik_dT = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)); break;
    case 21: dGamma_ik_dT = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
    case 22: dGamma_ik_dT = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)); break;
    case 23: dGamma_ik_dT = 0; break;
    case 24: dGamma_ik_dT = 0; break;
    case 25: dGamma_ik_dT = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[3]*_A[7]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dT[7]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_A[7]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[1]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[1]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[1]/(R*((T)*(T)*(T)))
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[2]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[2]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[2]/(R*((T)*(T)*(T)))
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[3]*_A[3]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dT[3]/(R*T) - std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[0]*Phi[5]*_A[3]/(R*((T)*(T)*(T)))
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[0]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dT[0]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[2]*_A[0]/(R*((T)*(T)*(T)))
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[8]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[8]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[8]/(R*((T)*(T)*(T)))
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[1]*Phi[5]*_A[9]/(R*((T)*(T)*(T)))
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dT[9]/(R*T) - std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T))) + T0*std::exp(-T0/T)*Phi[2]*Phi[4]*_A[9]/(R*((T)*(T)*(T)))
)); break;
  }
  return dGamma_ik_dT;
}
//-----------------------------------------------------------------------------
void  eclogitization_matlab_stx21_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi, std::vector<std::vector<double> >& _dGamma ) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  _dGamma = {
  {0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)), 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)), 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)), 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)), 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
))},
  {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)), -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)), 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)), 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
))},
  {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
))},
  {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)), -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
))},
  {0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)), 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)), 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)), 0, 0},
  {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
))}};
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> >  eclogitization_matlab_stx21_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi) const
{
  //FIXME: thread safety again?
  dGamma_ik_dP(T, P, C, Phi, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_ik_dP(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    const int i, const int k) const
{
  A(T, P, C, _A);
  dA_dP(T, P, C, _dA_dP);
  unsigned m = i*5 + k;
  double dGamma_ik_dP;

  switch(m) {
    case 0: dGamma_ik_dP = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)); break;
    case 1: dGamma_ik_dP = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)); break;
    case 2: dGamma_ik_dP = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)); break;
    case 3: dGamma_ik_dP = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)); break;
    case 4: dGamma_ik_dP = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)); break;
    case 5: dGamma_ik_dP = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dP[5]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dP[5]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)); break;
    case 6: dGamma_ik_dP = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 7: dGamma_ik_dP = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)); break;
    case 8: dGamma_ik_dP = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)); break;
    case 10: dGamma_ik_dP = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 15: dGamma_ik_dP = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)); break;
    case 16: dGamma_ik_dP = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[6]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dP[6]/(R*T)
)); break;
    case 20: dGamma_ik_dP = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)); break;
    case 21: dGamma_ik_dP = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
    case 22: dGamma_ik_dP = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)); break;
    case 23: dGamma_ik_dP = 0; break;
    case 24: dGamma_ik_dP = 0; break;
    case 25: dGamma_ik_dP = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dP[7]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dP[7]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dP[3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dP[3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dP[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[8]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dP[9]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dP[9]/(R*T)
)); break;
  }
  return dGamma_ik_dP;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);

    switch(i) {
    case 0: _dGamma = {
 {0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)),
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)),
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)),
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)),
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
))},
        {0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)),
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)),
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)),
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)),
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
))},
        {1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)),
 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)),
 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)),
 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)),
 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
))},
        {0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)),
 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)),
 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)),
 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)),
 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
))},
        {0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)),
 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)),
 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)),
 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)),
 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
))}
    }; break;

    case 1: _dGamma = {
 {-1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)),
 -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)),
 -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)),
 -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
))},
        {-0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
))},
        {0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)),
 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)),
 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)),
 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
))},
        {0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)),
 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)),
 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)),
 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
))}
    }; break;

    case 2: _dGamma = {
 {0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
))}
    }; break;

    case 3: _dGamma = {
 {-1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)),
 -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
))},
        {-1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)),
 -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
))}
    }; break;

    case 4: _dGamma = {
 {0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
))},
        {0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
))},
        {0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
))},
        {0,
 0,
 0,
 0,
 0},
        {0,
 0,
 0,
 0,
 0}
    }; break;

    case 5: _dGamma = {
 {0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::dGamma_ik_dC(
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
void eclogitization_matlab_stx21_rx::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    int i, std::vector<std::vector<double> >& _dGamma) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
    switch(i) {
    case 0: _dGamma = {
 {0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)),
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)),
 0,
 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)),
 0,
 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
))},
        {0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)),
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)),
 0,
 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)),
 0,
 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
))},
        {1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)),
 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)),
 0,
 0,
 0,
 0},
        {0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)),
 0,
 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)),
 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)),
 0,
 0},
        {0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)),
 0,
 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)),
 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)),
 0,
 0}
    }; break;

    case 1: _dGamma = {
 {-0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)),
 -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)),
 -0.38252915203572546*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)),
 -0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)),
 -0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)),
 -0.26516069900048528*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
))},
        {-0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)),
 -0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)),
 -0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)),
 -0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
        {0,
 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)),
 0.77105018965557404*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)),
 0,
 0,
 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
))},
        {0,
 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)),
 0,
 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)),
 0,
 0.57198294355527124*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
))}
    }; break;

    case 2: _dGamma = {
 {0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)),
 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)),
 0.22894981034442599*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)),
 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)),
 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)),
 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))}
    }; break;

    case 3: _dGamma = {
 {-0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)),
 -0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)),
 -1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)),
 -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)),
 -1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)),
 -0.73483930099951478*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
))},
        {-1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)),
 0,
 -1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)),
 -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)),
 0,
 0}
    }; break;

    case 4: _dGamma = {
 {0,
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)),
 0,
 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)),
 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
))},
        {0,
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)),
 0,
 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)),
 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))},
        {0,
 0,
 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)),
 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
))},
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

    case 5: _dGamma = {
 {0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)),
 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)),
 -0.61747084796427432*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)),
 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)),
 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)),
 0.42801705644472876*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
))}
    }; break;

  }

}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::dGamma_ik_dPhi(
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
double eclogitization_matlab_stx21_rx::dGamma_ik_dC(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l, unsigned m) const
{
  A(T, P, C, _A);
  dA_dC(T, P, C, _dA_dC);
  unsigned p = l*150 + m*30 + i*5+ k;
  double dGamma_ik_dC;

  switch(p) {
    case 0: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)); break;
    case 1: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)); break;
    case 2: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)); break;
    case 3: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)); break;
    case 4: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)); break;
    case 5: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)); break;
    case 6: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 7: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)); break;
    case 8: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)); break;
    case 10: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 15: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)); break;
    case 16: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][0]/(R*T)
)); break;
    case 20: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)); break;
    case 21: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 22: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)); break;
    case 23: dGamma_ik_dC = 0; break;
    case 24: dGamma_ik_dC = 0; break;
    case 25: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][0]/(R*T)
)); break;
    case 30: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)); break;
    case 31: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)); break;
    case 32: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)); break;
    case 33: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)); break;
    case 34: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)); break;
    case 35: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][1]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)); break;
    case 36: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 37: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)); break;
    case 38: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)); break;
    case 40: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 45: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)); break;
    case 46: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][1]/(R*T)
)); break;
    case 50: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)); break;
    case 51: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 52: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)); break;
    case 53: dGamma_ik_dC = 0; break;
    case 54: dGamma_ik_dC = 0; break;
    case 55: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][1]/(R*T)
)); break;
    case 60: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)); break;
    case 61: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)); break;
    case 62: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)); break;
    case 63: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)); break;
    case 64: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)); break;
    case 65: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][2]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)); break;
    case 66: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 67: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)); break;
    case 68: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)); break;
    case 70: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 75: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)); break;
    case 76: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][2]/(R*T)
)); break;
    case 80: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)); break;
    case 81: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 82: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)); break;
    case 83: dGamma_ik_dC = 0; break;
    case 84: dGamma_ik_dC = 0; break;
    case 85: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][2]/(R*T)
)); break;
    case 90: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)); break;
    case 91: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)); break;
    case 92: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)); break;
    case 93: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)); break;
    case 94: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)); break;
    case 95: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][3]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)); break;
    case 96: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 97: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)); break;
    case 98: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)); break;
    case 100: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 105: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)); break;
    case 106: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][3]/(R*T)
)); break;
    case 110: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)); break;
    case 111: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 112: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)); break;
    case 113: dGamma_ik_dC = 0; break;
    case 114: dGamma_ik_dC = 0; break;
    case 115: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][3]/(R*T)
)); break;
    case 120: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)); break;
    case 121: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)); break;
    case 122: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)); break;
    case 123: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)); break;
    case 124: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)); break;
    case 125: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][0][4]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)); break;
    case 126: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 127: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)); break;
    case 128: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)); break;
    case 130: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 135: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][0][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)); break;
    case 136: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][0][4]/(R*T)
)); break;
    case 140: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)); break;
    case 141: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 142: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)); break;
    case 143: dGamma_ik_dC = 0; break;
    case 144: dGamma_ik_dC = 0; break;
    case 145: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][0][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][0][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][0][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][0][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][0][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][0][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][0][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][0][4]/(R*T)
)); break;
    case 150: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)); break;
    case 151: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)); break;
    case 152: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)); break;
    case 153: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)); break;
    case 154: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)); break;
    case 155: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)); break;
    case 156: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 157: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)); break;
    case 158: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)); break;
    case 160: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 165: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)); break;
    case 166: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][0]/(R*T)
)); break;
    case 170: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)); break;
    case 171: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 172: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)); break;
    case 173: dGamma_ik_dC = 0; break;
    case 174: dGamma_ik_dC = 0; break;
    case 175: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][0]/(R*T)
)); break;
    case 180: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)); break;
    case 181: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)); break;
    case 182: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)); break;
    case 183: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)); break;
    case 184: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)); break;
    case 185: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][1]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)); break;
    case 186: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 187: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)); break;
    case 188: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)); break;
    case 190: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 195: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)); break;
    case 196: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][1]/(R*T)
)); break;
    case 200: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)); break;
    case 201: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 202: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)); break;
    case 203: dGamma_ik_dC = 0; break;
    case 204: dGamma_ik_dC = 0; break;
    case 205: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][1]/(R*T)
)); break;
    case 210: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)); break;
    case 211: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)); break;
    case 212: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)); break;
    case 213: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)); break;
    case 214: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)); break;
    case 215: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][2]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)); break;
    case 216: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 217: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)); break;
    case 218: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)); break;
    case 220: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 225: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)); break;
    case 226: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][2]/(R*T)
)); break;
    case 230: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)); break;
    case 231: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 232: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)); break;
    case 233: dGamma_ik_dC = 0; break;
    case 234: dGamma_ik_dC = 0; break;
    case 235: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][2]/(R*T)
)); break;
    case 240: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)); break;
    case 241: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)); break;
    case 242: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)); break;
    case 243: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)); break;
    case 244: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)); break;
    case 245: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][1][3]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)); break;
    case 246: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 247: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)); break;
    case 248: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)); break;
    case 250: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 255: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][1][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)); break;
    case 256: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][1][3]/(R*T)
)); break;
    case 260: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)); break;
    case 261: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 262: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)); break;
    case 263: dGamma_ik_dC = 0; break;
    case 264: dGamma_ik_dC = 0; break;
    case 265: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][1][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][1][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][1][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][1][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][1][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][1][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][1][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][1][3]/(R*T)
)); break;
    case 300: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)); break;
    case 301: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)); break;
    case 302: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)); break;
    case 303: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)); break;
    case 304: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)); break;
    case 305: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][2][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)); break;
    case 306: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 307: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)); break;
    case 308: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)); break;
    case 310: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 315: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][2][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)); break;
    case 316: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][2][0]/(R*T)
)); break;
    case 320: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)); break;
    case 321: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 322: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)); break;
    case 323: dGamma_ik_dC = 0; break;
    case 324: dGamma_ik_dC = 0; break;
    case 325: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][2][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][2][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][2][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][2][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][2][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][2][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][2][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][2][0]/(R*T)
)); break;
    case 450: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)); break;
    case 451: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)); break;
    case 452: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)); break;
    case 453: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)); break;
    case 454: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)); break;
    case 455: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)); break;
    case 456: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 457: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)); break;
    case 458: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)); break;
    case 460: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 465: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)); break;
    case 466: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][0]/(R*T)
)); break;
    case 470: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)); break;
    case 471: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 472: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)); break;
    case 473: dGamma_ik_dC = 0; break;
    case 474: dGamma_ik_dC = 0; break;
    case 475: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][0]/(R*T)
)); break;
    case 480: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)); break;
    case 481: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)); break;
    case 482: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)); break;
    case 483: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)); break;
    case 484: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)); break;
    case 485: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][3][1]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)); break;
    case 486: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 487: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)); break;
    case 488: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)); break;
    case 490: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 495: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][3][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)); break;
    case 496: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][3][1]/(R*T)
)); break;
    case 500: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)); break;
    case 501: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 502: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)); break;
    case 503: dGamma_ik_dC = 0; break;
    case 504: dGamma_ik_dC = 0; break;
    case 505: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][3][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][3][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][3][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][3][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][3][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][3][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][3][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][3][1]/(R*T)
)); break;
    case 600: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)); break;
    case 601: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)); break;
    case 602: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)); break;
    case 603: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)); break;
    case 604: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)); break;
    case 605: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)); break;
    case 606: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 607: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)); break;
    case 608: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)); break;
    case 610: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 615: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)); break;
    case 616: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][0]/(R*T)
)); break;
    case 620: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)); break;
    case 621: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 622: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)); break;
    case 623: dGamma_ik_dC = 0; break;
    case 624: dGamma_ik_dC = 0; break;
    case 625: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][0]/(R*T)
)); break;
    case 630: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)); break;
    case 631: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)); break;
    case 632: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)); break;
    case 633: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)); break;
    case 634: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)); break;
    case 635: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][1]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)); break;
    case 636: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 637: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)); break;
    case 638: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)); break;
    case 640: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 645: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][1]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)); break;
    case 646: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][1]/(R*T)
)); break;
    case 650: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)); break;
    case 651: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 652: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)); break;
    case 653: dGamma_ik_dC = 0; break;
    case 654: dGamma_ik_dC = 0; break;
    case 655: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][1]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][1]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][1]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][1]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][1]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][1]/(R*T)
)); break;
    case 660: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)); break;
    case 661: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)); break;
    case 662: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)); break;
    case 663: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)); break;
    case 664: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)); break;
    case 665: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][2]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)); break;
    case 666: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 667: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)); break;
    case 668: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)); break;
    case 670: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 675: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][2]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)); break;
    case 676: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][2]/(R*T)
)); break;
    case 680: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)); break;
    case 681: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 682: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)); break;
    case 683: dGamma_ik_dC = 0; break;
    case 684: dGamma_ik_dC = 0; break;
    case 685: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][2]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][2]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][2]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][2]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][2]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][2]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][2]/(R*T)
)); break;
    case 690: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)); break;
    case 691: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)); break;
    case 692: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)); break;
    case 693: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)); break;
    case 694: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)); break;
    case 695: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][3]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)); break;
    case 696: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 697: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)); break;
    case 698: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)); break;
    case 700: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 705: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)); break;
    case 706: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][3]/(R*T)
)); break;
    case 710: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)); break;
    case 711: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 712: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)); break;
    case 713: dGamma_ik_dC = 0; break;
    case 714: dGamma_ik_dC = 0; break;
    case 715: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][3]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][3]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][3]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][3]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][3]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][3]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][3]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][3]/(R*T)
)); break;
    case 720: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)); break;
    case 721: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)); break;
    case 722: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)); break;
    case 723: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)); break;
    case 724: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)); break;
    case 725: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][4][4]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)); break;
    case 726: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 727: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)); break;
    case 728: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)); break;
    case 730: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 735: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][4][4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)); break;
    case 736: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][4][4]/(R*T)
)); break;
    case 740: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)); break;
    case 741: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 742: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)); break;
    case 743: dGamma_ik_dC = 0; break;
    case 744: dGamma_ik_dC = 0; break;
    case 745: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][4][4]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][4][4]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][4][4]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][4][4]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][4][4]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][4][4]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][4][4]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][4][4]/(R*T)
)); break;
    case 750: dGamma_ik_dC = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)); break;
    case 751: dGamma_ik_dC = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)); break;
    case 752: dGamma_ik_dC = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)); break;
    case 753: dGamma_ik_dC = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)); break;
    case 754: dGamma_ik_dC = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)); break;
    case 755: dGamma_ik_dC = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_dA_dC[5][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*_dA_dC[5][5][0]/(R*T)
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)); break;
    case 756: dGamma_ik_dC = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 757: dGamma_ik_dC = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)); break;
    case 758: dGamma_ik_dC = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)); break;
    case 760: dGamma_ik_dC = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 765: dGamma_ik_dC = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[4][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[4][5][0]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)); break;
    case 766: dGamma_ik_dC = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[6][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[2]*_dA_dC[6][5][0]/(R*T)
)); break;
    case 770: dGamma_ik_dC = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)); break;
    case 771: dGamma_ik_dC = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
    case 772: dGamma_ik_dC = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)); break;
    case 773: dGamma_ik_dC = 0; break;
    case 774: dGamma_ik_dC = 0; break;
    case 775: dGamma_ik_dC = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_dA_dC[7][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*Phi[5]*_dA_dC[7][5][0]/(R*T)
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[1][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[1][5][0]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[2][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[2][5][0]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[3]*_dA_dC[3][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[0]*Phi[5]*_dA_dC[3][5][0]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[0][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[1]*Phi[2]*_dA_dC[0][5][0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[8][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[8][5][0]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*Phi[5]*_dA_dC[9][5][0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_dA_dC[9][5][0]/(R*T)
)); break;
  }
  return dGamma_ik_dC;
}
//-----------------------------------------------------------------------------
double eclogitization_matlab_stx21_rx::dGamma_ik_dPhi(
    const double& T, const double& P,
    std::vector<std::vector<double> >& C,
    std::vector<double>& Phi,
    unsigned i, unsigned k, unsigned l) const
{
  A(T, P, C, _A);
  unsigned m = l*30 + i*5 + k;
  double dGamma_ik_dPhi;

  switch(m) {
    case 0: dGamma_ik_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)); break;
    case 1: dGamma_ik_dPhi = 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)); break;
    case 2: dGamma_ik_dPhi = 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)); break;
    case 3: dGamma_ik_dPhi = 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)); break;
    case 4: dGamma_ik_dPhi = 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 5: dGamma_ik_dPhi = -0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 1.0*((_A[5] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*_A[5]/(R*T)
)); break;
    case 6: dGamma_ik_dPhi = -0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)); break;
    case 7: dGamma_ik_dPhi = 0; break;
    case 8: dGamma_ik_dPhi = 0; break;
    case 10: dGamma_ik_dPhi = 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 15: dGamma_ik_dPhi = -0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)) - 1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[4]/(R*T)
)); break;
    case 16: dGamma_ik_dPhi = -1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[6]/(R*T)
)); break;
    case 20: dGamma_ik_dPhi = 0; break;
    case 21: dGamma_ik_dPhi = 0; break;
    case 22: dGamma_ik_dPhi = 0; break;
    case 23: dGamma_ik_dPhi = 0; break;
    case 24: dGamma_ik_dPhi = 0; break;
    case 25: dGamma_ik_dPhi = 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[5]*_A[3]/(R*T)
)); break;
    case 30: dGamma_ik_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 31: dGamma_ik_dPhi = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 32: dGamma_ik_dPhi = 1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)); break;
    case 33: dGamma_ik_dPhi = 0; break;
    case 34: dGamma_ik_dPhi = 0; break;
    case 35: dGamma_ik_dPhi = -1.0*((_A[5] >= 0) ? (
   std::exp(-T0/T)*_A[5]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)); break;
    case 36: dGamma_ik_dPhi = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 37: dGamma_ik_dPhi = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)); break;
    case 38: dGamma_ik_dPhi = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)); break;
    case 40: dGamma_ik_dPhi = 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 45: dGamma_ik_dPhi = -0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 46: dGamma_ik_dPhi = 0; break;
    case 50: dGamma_ik_dPhi = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)); break;
    case 51: dGamma_ik_dPhi = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 52: dGamma_ik_dPhi = 0; break;
    case 53: dGamma_ik_dPhi = 0; break;
    case 54: dGamma_ik_dPhi = 0; break;
    case 55: dGamma_ik_dPhi = 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[1]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[5]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[3]*_A[3]/(R*T)
)
: (
   0
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[0]/(R*T)
)
: (
   std::exp(-T0/T)*Phi[2]*_A[0]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[5]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 60: dGamma_ik_dPhi = 0; break;
    case 61: dGamma_ik_dPhi = 0; break;
    case 62: dGamma_ik_dPhi = 0; break;
    case 63: dGamma_ik_dPhi = 0.78403044661833965*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)); break;
    case 64: dGamma_ik_dPhi = 0.77086564087833509*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)); break;
    case 65: dGamma_ik_dPhi = -0.38252915203572546*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)); break;
    case 66: dGamma_ik_dPhi = -0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 67: dGamma_ik_dPhi = 0.77105018965557404*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)); break;
    case 68: dGamma_ik_dPhi = 0; break;
    case 70: dGamma_ik_dPhi = 0.22894981034442599*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.21596955338166035*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 75: dGamma_ik_dPhi = -1.0*((_A[4] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[4]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)); break;
    case 76: dGamma_ik_dPhi = -1.0*((_A[6] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[6]/(R*T)
)); break;
    case 80: dGamma_ik_dPhi = 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)); break;
    case 81: dGamma_ik_dPhi = 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 82: dGamma_ik_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)); break;
    case 83: dGamma_ik_dPhi = 0; break;
    case 84: dGamma_ik_dPhi = 0; break;
    case 85: dGamma_ik_dPhi = -0.61747084796427432*((_A[0] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[4]*_A[9]/(R*T)
)); break;
    case 90: dGamma_ik_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 91: dGamma_ik_dPhi = 0.60489791348056932*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 92: dGamma_ik_dPhi = 0; break;
    case 93: dGamma_ik_dPhi = 0.78403044661833965*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)); break;
    case 94: dGamma_ik_dPhi = 0.77086564087833509*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)); break;
    case 95: dGamma_ik_dPhi = -0.26516069900048528*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)); break;
    case 96: dGamma_ik_dPhi = -0.32167062805152419*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 97: dGamma_ik_dPhi = 0; break;
    case 98: dGamma_ik_dPhi = 0.57198294355527124*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)); break;
    case 100: dGamma_ik_dPhi = 0.21596955338166035*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) + 0.22913435912166485*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)) + 0.071989851127220117*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)); break;
    case 105: dGamma_ik_dPhi = -1.0*((_A[4] >= 0) ? (
   std::exp(-T0/T)*_A[4]/(R*T)
)
: (
   0
)) - 1.0*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 106: dGamma_ik_dPhi = -1.0*((_A[6] >= 0) ? (
   std::exp(-T0/T)*_A[6]/(R*T)
)
: (
   0
)); break;
    case 110: dGamma_ik_dPhi = 0; break;
    case 111: dGamma_ik_dPhi = 0; break;
    case 112: dGamma_ik_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)); break;
    case 113: dGamma_ik_dPhi = 0; break;
    case 114: dGamma_ik_dPhi = 0; break;
    case 115: dGamma_ik_dPhi = 0.38830898661566587*((_A[7] >= 0) ? (
   std::exp(-T0/T)*_A[7]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[1] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)
: (
   0
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[2]/(R*T)
)
: (
   0
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[3]/(R*T)
)
: (
   0
)); break;
    case 120: dGamma_ik_dPhi = 0; break;
    case 121: dGamma_ik_dPhi = 0; break;
    case 122: dGamma_ik_dPhi = 0; break;
    case 123: dGamma_ik_dPhi = 0; break;
    case 124: dGamma_ik_dPhi = 0; break;
    case 125: dGamma_ik_dPhi = -0.65016971590119266*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)); break;
    case 126: dGamma_ik_dPhi = -0.70950810370425188*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 127: dGamma_ik_dPhi = 0; break;
    case 128: dGamma_ik_dPhi = 0; break;
    case 130: dGamma_ik_dPhi = 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 135: dGamma_ik_dPhi = -1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)); break;
    case 136: dGamma_ik_dPhi = 0; break;
    case 140: dGamma_ik_dPhi = 0.87028767842041899*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)); break;
    case 141: dGamma_ik_dPhi = 0.89228954729964483*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 142: dGamma_ik_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)); break;
    case 143: dGamma_ik_dPhi = 0; break;
    case 144: dGamma_ik_dPhi = 0; break;
    case 145: dGamma_ik_dPhi = 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[5]*_A[7]/(R*T)
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[8]/(R*T)
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*_A[9]/(R*T)
)); break;
    case 150: dGamma_ik_dPhi = 0.57198294355527124*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)); break;
    case 151: dGamma_ik_dPhi = 0.60489791348056932*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)); break;
    case 152: dGamma_ik_dPhi = 0; break;
    case 153: dGamma_ik_dPhi = 0; break;
    case 154: dGamma_ik_dPhi = 0; break;
    case 155: dGamma_ik_dPhi = -0.26516069900048528*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.26516069900048528*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.38252915203572546*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.65016971590119266*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)); break;
    case 156: dGamma_ik_dPhi = -0.32167062805152419*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 0.70950810370425188*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 157: dGamma_ik_dPhi = 0.77105018965557404*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)); break;
    case 158: dGamma_ik_dPhi = 0.57198294355527124*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)); break;
    case 160: dGamma_ik_dPhi = 0.071989851127220117*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) + 0.22894981034442599*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) + 0.12971232157958101*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) + 0.10771045270035519*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 165: dGamma_ik_dPhi = -0.73483930099951478*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) - 0.73483930099951478*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) - 0.67832937194847587*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) - 1.0*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)); break;
    case 166: dGamma_ik_dPhi = 0; break;
    case 170: dGamma_ik_dPhi = 0.87028767842041899*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)); break;
    case 171: dGamma_ik_dPhi = 0.89228954729964483*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
    case 172: dGamma_ik_dPhi = 0.53970116225711406*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)); break;
    case 173: dGamma_ik_dPhi = 0; break;
    case 174: dGamma_ik_dPhi = 0; break;
    case 175: dGamma_ik_dPhi = 0.42801705644472876*((_A[1] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[1]*_A[1]/(R*T)
)) + 0.42801705644472876*((_A[2] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[2]/(R*T)
)) + 0.39510208651943068*((_A[3] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[0]*_A[3]/(R*T)
)) + 0.38830898661566587*((_A[7] >= 0) ? (
   0
)
: (
   std::exp(-T0/T)*Phi[2]*Phi[4]*_A[7]/(R*T)
)) - 0.61747084796427432*((_A[0] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[0]/(R*T)
)
: (
   0
)) - 0.34983028409880706*((_A[8] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[8]/(R*T)
)
: (
   0
)) - 0.29049189629574806*((_A[9] >= 0) ? (
   std::exp(-T0/T)*Phi[1]*_A[9]/(R*T)
)
: (
   0
)); break;
  }
  return dGamma_ik_dPhi;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::rho(const double &T, const double &P,
		std::vector<std::vector<double> >& C,
        std::vector<double>&  _rho) const
{
  // inplace calculation of vector of phase densities (must be correct size)
  for (int i = 0; i < _phases.size(); i++ ) {
    _rho[i] = _phases[i]->rho(T,P,C[i]);
  }
}
//-----------------------------------------------------------------------------
std::vector<double> eclogitization_matlab_stx21_rx::rho(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _rho(6,0.);
  rho(T, P, C, _rho);
  return _rho;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::Cp(
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
std::vector<double> eclogitization_matlab_stx21_rx::Cp(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{

  // return vector of densities
  std::vector<double> _Cp(6,0.);
  Cp(T, P, C, _Cp);
  return _Cp;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::s(
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
std::vector<double> eclogitization_matlab_stx21_rx::s(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _s(6,0.);
  s(T, P, C, _s);
  return _s;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::ds_dC(
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
std::vector<std::vector<double> >eclogitization_matlab_stx21_rx::ds_dC(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  ds_dC(T, P, C, _tmp_ik);
  return _tmp_ik;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::alpha(
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
std::vector<double> eclogitization_matlab_stx21_rx::alpha(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _alpha(6,0.);
  alpha(T, P, C, _alpha);
  return _alpha;

}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::beta(
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
std::vector<double> eclogitization_matlab_stx21_rx::beta(
    const double &T, const double &P,
    std::vector<std::vector<double> >& C) const
{
  // return vector of densities
  std::vector<double> _beta(6,0.);
  beta(T, P, C, _beta);
  return _beta;

}

//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::Mu(
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
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::Mu(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  Mu(T,P,C,_Mu);
  return _Mu;
}

//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::dMu_dT(
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
void eclogitization_matlab_stx21_rx::dMu_dP(
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
void eclogitization_matlab_stx21_rx::dMu_dC(
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
std::vector<std::vector<std::vector<double> > > eclogitization_matlab_stx21_rx::dMu_dC(
    const double &T, const double &P,
		std::vector<std::vector<double> >& C) const
{
  // set and return private _Mu
  dMu_dC(T,P,C,_dMu_dC);
  return _dMu_dC;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::C_to_X(
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
void eclogitization_matlab_stx21_rx::X_to_C(
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
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::C_to_X(std::vector<std::vector<double> > &C) const
{
  //convert from weight fractions to Mole fractions and return X "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->c_to_x(C[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<double> > eclogitization_matlab_stx21_rx::X_to_C(std::vector<std::vector<double> > &X) const
{
  //convert from Mole fractions to weight fractions and return C "matrix"
  for (int i = 0; i < _phases.size(); i++)
  {
    _C[i] = _phases[i]->x_to_c(X[i]);
  }
  return _C;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::set_parameter(const std::string& p, const double& val) const
{
  *parameters[p] = val;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::get_parameter(const std::string& p) const
{
  std::cout << p << " = " << *parameters[p] << std::endl;
}
//-----------------------------------------------------------------------------
void eclogitization_matlab_stx21_rx::list_parameters() const
{
  std::cout << "Parameters: \n" << std::endl;
  for (auto const& x : parameters)
  {
    std::cout << x.first << " = "  << *x.second << std::endl;
  }
}
//-----------------------------------------------------------------------------
void  eclogitization_matlab_stx21_rx::report() const
{
  // loop over phases and endmemembers and spit out the names in order
  std::cout << "Reaction object: eclogitization_matlab_stx21_rx" << std::endl <<std::endl;
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
std::string eclogitization_matlab_stx21_rx::_to_string(const double& d) const
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
