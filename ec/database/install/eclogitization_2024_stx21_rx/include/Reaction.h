#ifndef REACTION_H
#define REACTION_H
 
#include <vector> 
#include <memory>
 
#include "Phase.h" 

// Abstract class for reactions
class Reaction
{
public:
  /// Constructor
  Reaction();

  /// Destructor
  virtual ~Reaction();

  /// The reaction name
  /// @return the name
  virtual std::string name() const = 0;

  /// TCG version at build time
  /// @return the TCG version at build time
  virtual std::string tcg_build_version() = 0;

  /// TCG git SHA at build time
  /// @return the TCG git SHA at build time
  virtual std::string tcg_build_git_sha() = 0;

  /// TCG version at generation time
  /// @return the TCG version at source code generation time
  virtual std::string tcg_generation_version() = 0;

  /// TCG git SHA at generation time
  /// @return the TCG git SHA  at source code generation time
  virtual std::string tcg_generation_git_sha() = 0;

  /// Vector of pointers to phases
  /// @return vector of pointers to phases
  virtual std::vector<std::shared_ptr<Phase> > phases() const = 0;

  // Return properly sized and zero'd composition vectors
  virtual std::vector<std::vector<double> > zero_C() const = 0;

  // Return composition vector of endmember masses
  virtual std::vector<std::vector<double> > M() const = 0;

  // Return "matrix" of molar stoichiometric coefficients
  virtual std::vector<std::vector<std::vector<double> > > nu() const = 0;

  // Return "matrix" of mass stoichiometric coefficients
  virtual std::vector<std::vector<std::vector<double> > > nu_m() const = 0;

  /// Sets a vector of affinities (of size J = the total number of 
  /// reactions), where the jth affinity is \f$ A_j = -\sum_{i,k} \mu_i^k*\nu_{ij}^k \f$. 
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @param[in] C Vector of vectors of concentrations
  /// @return _A Vector of affinities
  virtual std::vector<double> A(const double& T, const double& P, 
                                std::vector<std::vector<double> >& C) const = 0;

  /// Pass by reference interface that sets an internal vector of affinities. 
  /// @see A
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @param[in] C Vector of vectors of concentrations
  /// @param[out] _A Vector of affinities
  virtual void A(const double& T, const double& P, 
                 std::vector<std::vector<double> >& C, 
                 std::vector<double> &_A) const = 0;

  /// Function that returns the jth component of the affinity vector. 
  /// @see A
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @param[in] C Vector of vectors of concentrations
  /// @param[in] j Index of the affinity vector
  /// @return the affinity for the jth reaction
  virtual double A(const double& T, const double& P, 
                   std::vector<std::vector<double> >& C, 
                   const int j) const = 0;

  /// derivatives of affinities with respect to T
  /// pass by reference interface
  virtual void dA_dT(const double& T, const double& P, 
                     std::vector<std::vector<double> >& C,
                     std::vector<double> &_dA_dT) const = 0;

  /// derivatives of affinities with respect to T
  /// @return std::vector<double> of derivatives
  virtual std::vector<double> dA_dT(const double& T, const double& P,
                                    std::vector<std::vector<double> >& C) const = 0;

  virtual double dA_dT(const double& T, const double& P, 
                       std::vector<std::vector<double> >& C, const int j) const = 0;

  /// derivatives of affinities with respect to P
  /// pass by reference interface
  virtual void dA_dP(const double& T, const double& P, 
                     std::vector<std::vector<double> >& C,
                     std::vector<double> &_dA_dP) const = 0;

  virtual std::vector<double> dA_dP(const double& T, const double& P, 
                                    std::vector<std::vector<double> >& C) const = 0;

  virtual double dA_dP(const double& T, const double& P, 
                       std::vector<std::vector<double> >& C, 
                       const int j) const = 0;

  /// derivative of  affinities Aj with respect to C (''matrix'' of concentrations)
  virtual void dA_dC(const double& T, const double& P, 
                     std::vector<std::vector<double> >& C,
                     std::vector<std::vector<std::vector<double> > > &_dA_dC ) const = 0;

  virtual std::vector<std::vector<std::vector<double> > > dA_dC(const double& T, const double& P,
                                                                std::vector<std::vector<double> >& C) const = 0;
  
  // deprecate this:  just redundant to dA_dC but only returns one component
  virtual double dAj_dCik(const double& T, const double& P,
                          std::vector<std::vector<double> >& C,
                          unsigned j, unsigned i, unsigned k) const = 0;

  // Gamma_i:  reaction rates for all phases i, multiplie interfaces
  virtual void Gamma_i(const double& T, const double& P,
                       std::vector<std::vector<double> >& C,
                       std::vector<double>& Phi,
                       std::vector<double>& _Gamma) const = 0;

  virtual std::vector<double> Gamma_i(const double& T, const double& P,
                                      std::vector<std::vector<double> >& C,
                                      std::vector<double>& Phi) const = 0;

  virtual double Gamma_i(const double& T, const double& P,
                         std::vector<std::vector<double> >& C,
                         std::vector<double>& Phi, const int i) const = 0;

  // derivatives of mass transfer rates with respect to T
  virtual void dGamma_i_dT(const double& T, const double& P,
                           std::vector<std::vector<double> >& C,
                           std::vector<double>& Phi, std::vector<double>& _dGamma) const = 0;

  virtual std::vector<double> dGamma_i_dT(const double& T, const double& P,
                                          std::vector<std::vector<double> >& C,
                                          std::vector<double>& Phi) const = 0;

  virtual double dGamma_i_dT(const double& T, const double& P,
                             std::vector<std::vector<double> >& C,
                             std::vector<double>& Phi, const int i) const = 0;

  // derivatives of reaction rates with respect to P
  virtual void dGamma_i_dP(const double& T, const double& P,
                           std::vector<std::vector<double> >& C,
                           std::vector<double>& Phi, std::vector<double>& _dGamma) const = 0;


  virtual std::vector<double> dGamma_i_dP(const double& T, const double& P,
                                          std::vector<std::vector<double> >& C,
                                          std::vector<double>& Phi) const = 0;

  virtual double dGamma_i_dP(const double& T, const double& P,
                             std::vector<std::vector<double> >& C,
                             std::vector<double>& Phi, const int i) const = 0;

  // jacobians of mass transfer rates with respect to C
  virtual void dGamma_i_dC(const double& T, const double& P, 
                           std::vector<std::vector<double> >& C,
                           std::vector<double>& Phi,
                           std::vector<std::vector<std::vector<double> > >& _dGamma) const = 0;

  virtual std::vector<std::vector<std::vector<double> > > dGamma_i_dC(const double& T, const double& P,
                                                                      std::vector<std::vector<double> >& C,
                                                                      std::vector<double>& Phi) const = 0;

  virtual double dGamma_i_dC(const double& T, const double& P,
                             std::vector<std::vector<double> >& C,
                             std::vector<double>& Phi,
                             unsigned i, unsigned l, unsigned k) const = 0;

  // derivatives of phase mass transfer rates with vector Phi
  virtual void dGamma_i_dPhi(const double& T, const double& P,
                             std::vector<std::vector<double> >& C,std::vector<double>& Phi,
                             std::vector<std::vector<double> >& _dGamma) const = 0;

  virtual std::vector<std::vector<double> > dGamma_i_dPhi(const double& T, const double& P,
                                                          std::vector<std::vector<double> >& C,
                                                          std::vector<double>& Phi) const = 0;

  virtual double dGamma_i_dPhi(const double& T, const double& P,
                               std::vector<std::vector<double> >& C,
                               std::vector<double>& Phi,
                               unsigned i, unsigned l) const = 0;



  // Mass transfer rate for every component Gamma_ik: multiple interfaces
  virtual void Gamma_ik(const double& T, const double& P,
                        std::vector<std::vector<double> >& C,
                        std::vector<double>& Phi,
                        std::vector<std::vector<double> >& _Gamma) const = 0;

  virtual std::vector<std::vector<double> > Gamma_ik(const double& T, const double& P,
                                                     std::vector<std::vector<double> >& C,
                                                     std::vector<double>& Phi) const = 0;

  // Returns the mass transfer rate for component k in phase i
  virtual double Gamma_ik(const double& T, const double& P,
                          std::vector<std::vector<double> >& C,
                          std::vector<double>& Phi,
                          const int i, const int k) const = 0;

  // derivatives of compositional mass transfer rates with respect to T
  virtual void dGamma_ik_dT(const double& T, const double& P,
                            std::vector<std::vector<double> >& C,
                            std::vector<double>& Phi,
                            std::vector<std::vector<double> >& _dGamma) const = 0;

  virtual std::vector<std::vector<double> > dGamma_ik_dT(const double& T, const double& P,
                                                         std::vector<std::vector<double> >& C,
                                                         std::vector<double>& Phi) const = 0;

  virtual double dGamma_ik_dT(const double& T, const double& P,
                              std::vector<std::vector<double> >& C,
                              std::vector<double>& Phi,
                              const int i, const int k) const = 0;

  //  derivatives of component mass transfer rates with respect to P
  virtual void dGamma_ik_dP(const double& T, const double& P,
                            std::vector<std::vector<double> >& C,
                            std::vector<double>& Phi,
                            std::vector<std::vector<double> >& _dGamma) const = 0;

  virtual std::vector<std::vector<double> > dGamma_ik_dP(const double& T, const double& P,
                                                         std::vector<std::vector<double> >& C,
                                                         std::vector<double>& Phi) const = 0;

  // Returns dGamma_ik_dP for component k phase i
  virtual double dGamma_ik_dP(const double& T, const double& P,
                              std::vector<std::vector<double> >& C,
                              std::vector<double>& Phi,
                              const int i, const int k) const = 0;

  // derivative of components with respect to Composition for phase i
  virtual void dGamma_ik_dC(const double& T, const double& P,
                            std::vector<std::vector<double> >& C,
                            std::vector<double>& Phi,
                            int i, std::vector<std::vector<double> >& _dGamma ) const = 0;

  virtual std::vector<std::vector<double> > dGamma_ik_dC(const double& T, const double& P,
                                                         std::vector<std::vector<double> >& C,
                                                         std::vector<double>& Phi,
                                                         int i) const = 0;

  // derivative of components with respect to Phases for phase i
  virtual void  dGamma_ik_dPhi(const double& T, const double& P,
                               std::vector<std::vector<double> >& C,
                               std::vector<double>& Phi,
                               int i, std::vector<std::vector<double> >& _dGamma) const = 0;

  virtual std::vector<std::vector<double> > dGamma_ik_dPhi(const double& T, const double& P,
                                                           std::vector<std::vector<double> >& C,
                                                           std::vector<double>& Phi,
                                                           int i) const = 0;

  // Derivative of Gamma_i^k with respect to lth phase fraction
  virtual double dGamma_ik_dPhi(const double& T, const double& P,
                                std::vector<std::vector<double> >& C,
                                std::vector<double>& Phi,
                                unsigned i, unsigned k, unsigned l) const = 0;

  // Derivative of Gamma_i^k with respect to lth phase fraction Clm
  virtual double dGamma_ik_dC(const double& T, const double& P,
                              std::vector<std::vector<double> >& C,
                              std::vector<double>& Phi,
                              unsigned i, unsigned k, unsigned l, unsigned m) const = 0;

  // density of phases (was rho_phases, but deprecating "reaction densities")
  virtual void rho(const double& T, const double& P, 
                   std::vector<std::vector<double> >& C, 
                   std::vector<double> & _rho) const = 0;

  virtual std::vector<double> rho(const double& T, const double& P, 
                                  std::vector<std::vector<double> >& C ) const = 0;

  // phase Heat Capacity in J/K/g   void Cp(T,P,C,result)
  virtual void Cp(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _Cp) const = 0;

  // phase Heat Capacity Cp =  Cp(T,P,C)
  virtual std::vector<double> Cp(const double& T, const double& P, std::vector<std::vector<double> >& C ) const = 0;

  // phase entropy's  in units of J/K/g (s(n)/M(n)) void s(T,P,C,result)
  virtual void s(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _s) const = 0;

  // Entropy s =  s(T,P,C)
  virtual std::vector<double> s(const double& T, const double& P, std::vector<std::vector<double> >& C ) const = 0;

  // compositional derivative of  void ds_dC(T,P,C,result)
  virtual void ds_dC(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<std::vector<double> >& _ds) const = 0;

  // Entropy s =  s(T,P,C)
  virtual std::vector<std::vector<double> > ds_dC(const double& T, const double& P, std::vector<std::vector<double> >& C ) const = 0;

  // Thermal Expansivity void alpha(T,P,C,result)
  virtual void alpha(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _alpha) const = 0;

  // Thermal Expansivity alpha =  alpha(T,P,C)
  virtual std::vector<double> alpha(const double& T, const double& P, std::vector<std::vector<double> >& C ) const = 0;

  // Compressibility void beta(T,P,C,result)
  virtual void beta(const double& T, const double& P, std::vector<std::vector<double> >& C, std::vector<double> & _beta) const = 0;

  // Compressibility beta =  beta(T,P,C)
  virtual std::vector<double> beta(const double& T, const double& P, std::vector<std::vector<double> >& C ) const = 0;

  // return "Matrix" of chemical potentials \mu_i^k for all endmembers k in all phases i in-place version
  virtual void Mu(const double& T, const double& P, std::vector<std::vector<double> >& C,  std::vector<std::vector<double> >& _mu) const = 0;

  // return "Matrix" of chemical potentials \mu_i^k for all endmembers k in all phases i return matrix version
  virtual std::vector<std::vector<double> > Mu(const double& T, const double& P, std::vector<std::vector<double> >& C) const = 0;

  // return "Matrix" of derivatives of chemical potentials with respect to T
  // d\mu_i^k_dT  for all endmembers k in all phases i in-place version
  virtual void dMu_dT(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const = 0;

  // return "Matrix" of derivatives of chemical potentials with respect to P
  // d\mu_i^k_dPfor all endmembers k in all phases i in-place version
  virtual void dMu_dP(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<double> >& _dmu) const = 0;

  // return "Matrix" of derivatives of chemical potentials with respect to C
  // d\mu_i^k_dC for all endmembers k in all phases i in-place version
  virtual void dMu_dC(const double& T, const double& P, std::vector<std::vector<double> >& C,
    std::vector<std::vector<std::vector<double> > >& _dmu) const = 0;

  // return "Matrix" of derivatives of chemical potentials with respect to C
  // d\mu_i^k_dC for all endmembers k in all phases i
  virtual std::vector<std::vector<std::vector<double> > > dMu_dC(const double& T, const double& P,
    std::vector<std::vector<double> >& C) const = 0;

  // Convert from concentration to Mole fraction inplace
  virtual void C_to_X(std::vector<std::vector<double> >& C, std::vector<std::vector<double> >& X) const = 0;

  // Convert from concentration to Mole fraction with returned vector
  virtual std::vector<std::vector<double> >  C_to_X(std::vector<std::vector<double> >& C) const = 0;

  // Convert from Mole Fraction to Concentration inplace
  virtual void X_to_C(std::vector<std::vector<double> >& X, std::vector<std::vector<double> >& C) const = 0;

  // Convert from   Mole fraction to weight  with returned vector
  virtual std::vector<std::vector<double> >  X_to_C(std::vector<std::vector<double> >& X) const = 0;

  // Change parameter value
  virtual void set_parameter(const std::string& p, const double& val) const = 0;

  // Print parameter value
  virtual void get_parameter(const std::string& p) const = 0;

  // List available parameters
  virtual void list_parameters() const = 0;

  // Pretty print information about the reactions and endmembers
  virtual void report() const = 0;
};

#endif