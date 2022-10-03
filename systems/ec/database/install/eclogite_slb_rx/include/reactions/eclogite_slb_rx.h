#ifndef ECLOGITE_SLB_RX_H
#define ECLOGITE_SLB_RX_H

#include "Reaction.h"
#include "phases.h"

class eclogite_slb_rx: public Reaction
{
public:
  // Constructor
  eclogite_slb_rx();

  // Destructor
  virtual ~eclogite_slb_rx();

  // Reaction name
  std::string name() const;

  // TCG version at build time
  std::string tcg_build_version();

  // TCG git SHA at build time
  std::string tcg_build_git_sha();

  // TCG version at generation time
  std::string tcg_generation_version();

  // TCG git SHA at generation time
  std::string tcg_generation_git_sha();

  // Return vector of pointers to phases
  std::vector<std::shared_ptr<Phase> > phases() const;

  // // Return shared pointer to phase_i
  // std::shared_ptr<Phase> get_phase(const unsigned &i) const;

  // Return properly sized and zero'd composition vectors
  std::vector<std::vector<double> > zero_C() const;

  // Return composition vector of endmember masses
  std::vector<std::vector<double> > M() const;

  // Return "matrix" of molar stoichiometric coefficients
  std::vector<std::vector<std::vector<double> > > nu() const;

  // Return "matrix" of mass stoichiometric coefficients
  std::vector<std::vector<std::vector<double> > > nu_m() const;

  // Affinities for j reactions where A_j = -sum_i,k mu*nu_j: multiple interfaces
  void A(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> &_A) const;
  std::vector<double> A(const double& T, const double& P, std::vector<std::vector<double> >& C) const;
  double  A(const double& T, const double& P, std::vector<std::vector<double> >& C, const int j) const;

  // derivatives of affinities with respect to T
  void dA_dT(const double& T, const double& P, std::vector<std::vector<double> >& C,
             std::vector<double> &_dA_dT) const;

  std::vector<double> dA_dT(const double& T, const double& P, std::vector<std::vector<double> >& C) const;

  double dA_dT(const double& T, const double& P, std::vector<std::vector<double> >& C, const int j) const;

  // derivatives of affinities with respect to P
  void dA_dP(const double& T, const double& P, std::vector<std::vector<double> >& C,
             std::vector<double> &_dA_dP) const;

  std::vector<double> dA_dP(const double& T, const double& P, std::vector<std::vector<double> >& C) const;

  double dA_dP(const double& T, const double& P, std::vector<std::vector<double> >& C, const int j) const;

  // derivative of  affinities Aj with respect to C (matrix of concentrations)
  void dA_dC(const double& T, const double& P, std::vector<std::vector<double> >& C,
        std::vector<std::vector<std::vector<double> > > &_dA_dC ) const;

  std::vector<std::vector<std::vector<double> > > dA_dC(const double& T, const double& P,
                                                        std::vector<std::vector<double> >& C) const;
  // FIXME:: this returns the wrong jacobian component of jacobian of affinities Aj with respect to Cik (
  double dAj_dCik(const double& T, const double& P,
		              std::vector<std::vector<double> >& C,
                  unsigned j, unsigned i, unsigned k) const;

  // Gamma_i:  reaction rates for all phases i, multiplie interfaces
  void Gamma_i(const double& T, const double& P,
               std::vector<std::vector<double> >& C,
               std::vector<double>& Phi,
               std::vector<double>& _Gamma) const;

  std::vector<double> Gamma_i(const double& T, const double& P,
                              std::vector<std::vector<double> >& C,
                              std::vector<double>& Phi) const;

  double Gamma_i(const double& T, const double& P,
                 std::vector<std::vector<double> >& C,
                 std::vector<double>& Phi, const int i) const;

  // derivatives of mass transfer rates with respect to T
  void dGamma_i_dT(const double& T, const double& P,
                    std::vector<std::vector<double> >& C,
                    std::vector<double>& Phi, std::vector<double>& _dGamma) const;

  std::vector<double> dGamma_i_dT(const double& T, const double& P,
                                  std::vector<std::vector<double> >& C,
                                  std::vector<double>& Phi) const;

  double dGamma_i_dT(const double& T, const double& P,
                     std::vector<std::vector<double> >& C,
                     std::vector<double>& Phi, const int i) const;

  // derivatives of reaction rates with respect to P
  void dGamma_i_dP(const double& T, const double& P,
                   std::vector<std::vector<double> >& C,
                   std::vector<double>& Phi, std::vector<double>& _dGamma) const;


  std::vector<double> dGamma_i_dP(const double& T, const double& P,
                                  std::vector<std::vector<double> >& C,
                                  std::vector<double>& Phi) const;

  double dGamma_i_dP(const double& T, const double& P,
                     std::vector<std::vector<double> >& C,
                     std::vector<double>& Phi, const int i) const;

  // jacobians of mass transfer rates with respect to C
  void dGamma_i_dC(const double& T, const double& P, std::vector<std::vector<double> >& C,
                   std::vector<double>& Phi,
                   std::vector<std::vector<std::vector<double> > >& _dGamma) const;

  std::vector<std::vector<std::vector<double> > > dGamma_i_dC(const double& T, const double& P,
                                                              std::vector<std::vector<double> >& C,
                                                              std::vector<double>& Phi) const;

  double dGamma_i_dC(const double& T, const double& P,
                     std::vector<std::vector<double> >& C,
                     std::vector<double>& Phi,
                     unsigned i, unsigned l, unsigned k) const;

  // derivatives of phase mass transfer rates with vector Phi
  void dGamma_i_dPhi(const double& T, const double& P,
                     std::vector<std::vector<double> >& C,std::vector<double>& Phi,
                     std::vector<std::vector<double> >& _dGamma) const;

  std::vector<std::vector<double> > dGamma_i_dPhi(const double& T, const double& P,
                                                  std::vector<std::vector<double> >& C,
                                                  std::vector<double>& Phi) const;

  double dGamma_i_dPhi(const double& T, const double& P,
                       std::vector<std::vector<double> >& C,
                       std::vector<double>& Phi,
                       unsigned i, unsigned l) const;



  // Mass transfer rate for every component Gamma_ik: multiple interfaces
  void Gamma_ik(const double& T, const double& P,
                std::vector<std::vector<double> >& C,
                std::vector<double>& Phi,
                std::vector<std::vector<double> >& _Gamma) const;

  std::vector<std::vector<double> > Gamma_ik(const double& T, const double& P,
                                             std::vector<std::vector<double> >& C,
                                             std::vector<double>& Phi) const;

  // Returns the mass transfer rate for component k in phase i
  double Gamma_ik(const double& T, const double& P,
                  std::vector<std::vector<double> >& C,
                  std::vector<double>& Phi,
                  const int i, const int k) const;

  // derivatives of compositional mass transfer rates with respect to T
  void dGamma_ik_dT(const double& T, const double& P,
                    std::vector<std::vector<double> >& C,
                    std::vector<double>& Phi,
                    std::vector<std::vector<double> >& _dGamma) const;

  std::vector<std::vector<double> > dGamma_ik_dT(const double& T, const double& P,
                                                 std::vector<std::vector<double> >& C,
                                                 std::vector<double>& Phi) const;

  double dGamma_ik_dT(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      const int i, const int k) const;

  //  derivatives of component mass transfer rates with respect to P
  void dGamma_ik_dP(const double& T, const double& P,
                    std::vector<std::vector<double> >& C,
                    std::vector<double>& Phi,
                    std::vector<std::vector<double> >& _dGamma) const;

  std::vector<std::vector<double> > dGamma_ik_dP(const double& T, const double& P,
                                                 std::vector<std::vector<double> >& C,
                                                 std::vector<double>& Phi) const;

  // Returns dGamma_ik_dP for component k phase i
  double dGamma_ik_dP(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      const int i, const int k) const;

  // derivative of components with respect to Composition for phase i
  void dGamma_ik_dC(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      int i, std::vector<std::vector<double> >& _dGamma ) const;

  std::vector<std::vector<double> > dGamma_ik_dC(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      int i) const;

  // derivative of components with respect to Phases for phase i
  void  dGamma_ik_dPhi(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      int i, std::vector<std::vector<double> >& _dGamma) const;

  std::vector<std::vector<double> > dGamma_ik_dPhi(const double& T, const double& P,
                      std::vector<std::vector<double> >& C,
                      std::vector<double>& Phi,
                      int i) const;

  // Derivative of Gamma_i^k with respect to lth phase fraction
  double dGamma_ik_dPhi(const double& T, const double& P,
                        std::vector<std::vector<double> >& C,
                        std::vector<double>& Phi,
                        unsigned i, unsigned k, unsigned l) const;

  // Derivative of Gamma_i^k with respect to lth phase fraction Clm
  double dGamma_ik_dC(const double& T, const double& P,
                     std::vector<std::vector<double> >& C,
                     std::vector<double>& Phi,
                     unsigned i, unsigned k, unsigned l, unsigned m) const;

  // density of phases (was rho_phases, but deprecating "reaction densities")
  void rho(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _rho) const;

  std::vector<double> rho(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // phase Heat Capacity in J/K/g   void Cp(T,P,C,result)
  void Cp(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _Cp) const;

  // phase Heat Capacity Cp =  Cp(T,P,C)
  std::vector<double> Cp(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // phase entropy's  in units of J/K/g (s(n)/M(n)) void s(T,P,C,result)
  void s(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _s) const;

  // Entropy s =  s(T,P,C)
  std::vector<double> s(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // compositional derivative of  void ds_dC(T,P,C,result)
  void ds_dC(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<std::vector<double> >& _ds) const;

  // Entropy s =  s(T,P,C)
  std::vector<std::vector<double> > ds_dC(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // Thermal Expansivity void alpha(T,P,C,result)
  void alpha(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _alpha) const;

  // Thermal Expansivity alpha =  alpha(T,P,C)
  std::vector<double> alpha(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // Compressibility void beta(T,P,C,result)
  void beta(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _beta) const;

  // Compressibility beta =  beta(T,P,C)
  std::vector<double> beta(const double& T, const double& P, std::vector<std::vector<double> >& C ) const;

  // return "Matrix" of chemical potentials \mu_i^k for all endmembers k in all phases i in-place version
  void Mu(const double& T, const double& P, std::vector<std::vector<double> >& C,  std::vector<std::vector<double> >& _mu) const;

  // return "Matrix" of chemical potentials \mu_i^k for all endmembers k in all phases i return matrix version
  std::vector<std::vector<double> > Mu(const double& T, const double& P, std::vector<std::vector<double> >& C) const;

  // return "Matrix" of deriviatves of chemical potentials with respect to T
  // d\mu_i^k_dT  for all endmembers k in all phases i in-place version
  void dMu_dT(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const;

  // return "Matrix" of deriviatves of chemical potentials with respect to P
  // d\mu_i^k_dPfor all endmembers k in all phases i in-place version
  void dMu_dP(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const;

  // return "Matrix" of deriviatves of chemical potentials with respect to C
  // d\mu_i^k_dC for all endmembers k in all phases i in-place version
  void dMu_dC(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<std::vector<double> > >& _dmu) const;

  // return "Matrix" of deriviatves of chemical potentials with respect to C
  // d\mu_i^k_dC for all endmembers k in all phases i
  std::vector<std::vector<std::vector<double> > > dMu_dC(const double& T, const double& P,
    std::vector<std::vector<double> >& C) const;


  // Convert from concentration to Mole fraction inplace
  void C_to_X(std::vector<std::vector<double> >& C, std::vector<std::vector<double> >& X) const;

  // Convert from concentration to Mole fraction with returned vector
  std::vector<std::vector<double> >  C_to_X(std::vector<std::vector<double> >& C) const;

  // Convert from Mole Fraction to Concentration inplace
  void X_to_C(std::vector<std::vector<double> >& X, std::vector<std::vector<double> >& C) const;

  // Convert from   Mole fraction to weight  with returned vector
  std::vector<std::vector<double> >  X_to_C(std::vector<std::vector<double> >& X) const;

  // Change parameter value
  void set_parameter(const std::string& p, const double& val) const;

  // Print parameter value
  void get_parameter(const std::string& p) const;

  // List available parameters
  void list_parameters() const;

  // Pretty print information about the reactions and endmembers
  void report() const;

private:
  // Reaction name
  const std::string _name;

  // Number of reactions
  const int J;

  // Number of phases
  int N;

  // Vector of pointers to phases
  std::vector<std::shared_ptr<Phase> > _phases;

  // Zero'd composition "matrix" of the correct size for creating scratch compostions
  std::vector<std::vector<double> > _C0;

  // Molar Masses of Endmembers
  const std::vector<std::vector<double> > _M;

  // Molar Stoichiometric Constants nu[j][i][k]
  const std::vector<std::vector<std::vector<double> > > _nu;

  // Mass Stoichiometric Constants nu_m[j][i][k]
  const std::vector<std::vector<std::vector<double> > > _nu_m;

  // Local storage for compositions "matrix" of the correct size for creating scratch compostions
  mutable std::vector<std::vector<double> > _C;

  // Local storage for compositions "matrix" of the correct size for creating scratch mol fractions
  mutable std::vector<std::vector<double> > _X;

  // Temporary local  storage for Gamma_i, dGamma_i, Gamma_ik, and dGamma_ik
  mutable std::vector<double> _tmp;
  mutable std::vector<std::vector<std::vector<double> > > _tmp_dC;
  mutable std::vector<std::vector<double> > _tmp_dPhi;
  mutable std::vector<std::vector<double> > _tmp_ik;
  mutable std::vector<std::vector<double> > _tmp_ik_dC;
  mutable std::vector<std::vector<double> > _tmp_ik_dPhi;

  // temporary Local storage for chemical potential Matrix _Mu_i^k and derivatives
  mutable std::vector<std::vector<double> > _Mu;
  mutable std::vector<std::vector<double> > _dMu_dT;
  mutable std::vector<std::vector<double> > _dMu_dP;
  mutable std::vector<std::vector<std::vector<double> > > _dMu_dC;

  // Local storage for Affinity Vector and derivatives
  mutable std::vector<double>  _A;
  mutable std::vector<double>  _dA_dT;
  mutable std::vector<double>  _dA_dP;
  mutable std::vector<std::vector<std::vector<double> > > _dA_dC;

  // Parameters
  double T0;
  double R;

  // Container for mapping parameter string symbols to values
  mutable std::map<std::string, double*> parameters;
  
  // Convenience function for pretty printing doubles
  std::string _to_string(const double& d) const;
};

#endif
