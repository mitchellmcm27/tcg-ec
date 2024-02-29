#ifndef PHASE_H
#define PHASE_H

#include <vector>
#include <memory>

#include "EndMember.h"

// Abstract class for phases
class Phase
{
public:
  /// Constructor
  Phase();

  /// Destructor
  virtual ~Phase();

  //**************************************************************************
  // Coder functions
  //**************************************************************************

  /// Returns a unique identifier for the model instance
  /// @return the identifier
  virtual std::string identifier()= 0;
  
  /// name
  /// @return The phase name
  virtual std::string name() = 0;

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

  /// formula
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @param[in] n Vector of moles of components
  /// @return The formula
  virtual std::string formula(const double& T, const double& P,
                              std::vector<double>& n) = 0;


  // Extra Coder functions provided by Cython in ThermoEngine and written in C++ here

  /// Convert element moles to moles of endmembers (coder function)
  /// @return vector of moles of endmembers \f$\mathbf{n}\f$
  virtual std::vector<double> conv_elm_to_moles(std::vector<double>& e) = 0;

  /// Convert element moles to total moles of endmembers (coder function)
  /// @return  \f$\sum n\f$
  virtual double conv_elm_to_tot_moles(std::vector<double>& e) = 0;

  /// Convert element moles to total mass (in grams) of phase  (coder function)
  virtual double conv_elm_to_tot_grams(std::vector<double>& e) = 0;

  /// Convert moles of endmembers to mols of indivdual elements  (coder function)
  /// This routine requires a properly formatted conversion string in the construction of the
  /// coder phase model
  virtual std::vector<double> conv_moles_to_elm(std::vector<double>& n) = 0;

  /// Returns total moles of endmembers \f$\sum n\f$ (coder function)
  virtual double conv_moles_to_tot_moles(std::vector<double>& n) = 0;

  /// converts moles of endmembers to mol fractions of endmembers \f$ x_i = n_i/\sum_j n_j\f$
  /// @return vector \f$\mathbf{x}\f$
  virtual std::vector<double> conv_moles_to_mole_frac(std::vector<double>& n) = 0;

  /// test_moles
  virtual int test_moles(std::vector<double>& n) = 0;

  /// number of endmembers in the phase
  virtual int endmember_number() const = 0;
  
  /// name of endmember with index \f$ i \f$
  virtual std::string endmember_name(const int &i) const = 0;
  
  /// returns stoichiometric formula of endmember with index \f$ i \f$
  virtual std::string endmember_formula(const int &i) const = 0;
  
  /// returns Molecular weight  of endmember with index \f$ i \f$
  virtual double endmember_mw(const int &i) const = 0;
  
  // endmember_elements
  virtual std::vector<double> endmember_elements(const int &i) const = 0;

  // species_number
  virtual int species_number() const = 0;
  
  // species_name
  virtual std::string species_name(const int &i) const = 0;
  
  // species_formula
  virtual std::string species_formula(const int &i) const = 0;
  
  // species_mw
  virtual double species_mw(const int &i) const = 0;
  
  // species_elements
  virtual std::vector<double> species_elements(const int &i) const = 0;

  /// Returns the Gibbs free energy \f$G\f$
  virtual double g(const double& T, const double& P,
                   std::vector<double>& n) const = 0;

  /// Returns the partial derivative \f$\partial G/\partial T\f$
  virtual double dgdt(const double& T, const double& P,
                      std::vector<double>& n) const = 0;

  /// Returns the partial derivative \f$\partial G/\partial P\f$
  virtual double dgdp(const double& T, const double& P,
                      std::vector<double>& n) const = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial T^2\f$
  virtual double d2gdt2(const double& T, const double& P,
                        std::vector<double>& n) const = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial T \partial P\f$
  virtual double d2gdtdp(const double& T, const double& P,
                         std::vector<double>& n) const = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial P^2\f$
  virtual double d2gdp2(const double& T, const double& P,
                        std::vector<double>& n) const = 0;
                        
  /// Returns the partial derivative \f$\partial^3 G/\partial T^3\f$
  virtual double d3gdt3(const double& T, const double& P,
                        std::vector<double>& n) const = 0;
                        
  /// Returns the partial derivative \f$\partial^3 G/\partial T^2 \partial P\f$
  virtual double d3gdt2dp(const double& T, const double& P,
                          std::vector<double>& n) const = 0;
                          
  /// Returns the partial derivative \f$\partial^3 G/\partial T \partial P^2\f$
  virtual double d3gdtdp2(const double& T, const double& P,
                          std::vector<double>& n) const = 0;
                  
  /// Returns the partial derivative \f$\partial^3 G/\partial P^3\f$
  virtual double d3gdp3(const double& T, const double& P,
                        std::vector<double>& n) const = 0;

  /// Returns the chemical potential partial derivative \f$\mu = \partial G/\partial n\f$
  /// @return a vector of chemical potentials for all endmembers in the phase
  virtual std::vector<double> dgdn(const double& T, const double& P,
                                   std::vector<double>& n) const = 0;
          
  /// Returns the derivative of chemical potential with respect to T \f$\partial^2 G/\partial n\partial T\f$
  /// @return a vector of \f$\partial \mu/\partial T\f$ for all endmembers in the phase
  virtual std::vector<double> d2gdndt(const double& T, const double& P,
                                      std::vector<double>& n) const = 0;                         

  /// Returns the derivative of chemical potential with respect to P  \f$\partial^2 G/\partial n\partial P\f$
  /// @return a vector of \f$\partial \mu/\partial P\f$ for all endmembers in the phase
  virtual std::vector<double> d2gdndp(const double& T, const double& P,
                                      std::vector<double>& n) const = 0;
  
  /// Returns \f$\partial^2 G/\partial n^2\f$
  /// @return a vector for derivates with endmembers in the phase
  virtual std::vector<double> d2gdn2(const double& T, const double& P,
                                     std::vector<double>& n) const = 0;
                                                                    
  /// Returns \f$\partial^3 G/\partial n\partial T^2\f$
  virtual std::vector<double> d3gdndt2(const double& T, const double& P,
                                       std::vector<double>& n) const = 0;

  /// Returns \f$\partial^3 G/\partial n\partial T\partial P\f$
  virtual std::vector<double> d3gdndtdp(const double& T, const double& P,
                                        std::vector<double>& n) const = 0;
                                        
  /// Returns \f$\partial^3 G/\partial n^2\partial T\f$
  virtual std::vector<double> d3gdn2dt(const double& T, const double& P,
                                       std::vector<double>& n) const = 0;
                                       
  /// Returns \f$\partial^3 G/\partial n\partial P^2\f$
  virtual std::vector<double> d3gdndp2(const double& T, const double& P,
                                       std::vector<double>& n) const = 0;

  /// Returns \f$\partial^3 G/\partial n^2\partial P\f$
  virtual std::vector<double> d3gdn2dp(const double& T, const double& P,
                                       std::vector<double>& n) const = 0; 

  /// Returns \f$\partial^3 G/\partial n^3\f$
  virtual std::vector<double> d3gdn3(const double& T, const double& P,
                                       std::vector<double>& n) const = 0;
                                       
  /// Convenience function for volume of phase: returns \f$ v = \partial G/\partial P\f$
  virtual double v(const double& T, const double& P,
                   std::vector<double>& n) const = 0;

  /// Convenience function for entropy of phase: returns \f$ s = -\partial G/\partial T\f$
  virtual double s(const double& T, const double& P,
                   std::vector<double>& n) const = 0;

  /// Convenience function for thermal expansivity of phase: returns \f$ \alpha = (1/v)\partial v/\partial T\f$
  virtual double alpha(const double& T, const double& P,
                       std::vector<double>& n) const = 0;

  /// Convenience function for heat capacity (constant volume) of phase
  virtual double cv(const double& T, const double& P,
                    std::vector<double>& n) const = 0;

  /// Convenience function for heat capacity (constant pressure) of phase
  virtual double cp(const double& T, const double& P,
                    std::vector<double>& n) const = 0;

  /// Convenience function for derivative of heat capacity with respect to \f$ T\f$
  virtual double dcpdt(const double& T, const double& P,
                       std::vector<double>& n) const = 0;
                      
  /// Convenience function for compressibility \f$ \beta = (1/v) \partial v/\partial P\f$
  virtual double beta(const double& T, const double& P,
                      std::vector<double>& n) const = 0;
                     
  /// K
  virtual double K(const double& T, const double& P,
                   std::vector<double>& n) const = 0;
                  
  /// Kp
  virtual double Kp(const double& T, const double& P,
                    std::vector<double>& n) const = 0; 

  /// Returns the number of parameters
  /// @return number of calibratable parameters in phase model
  virtual int get_param_number() = 0;

  /// Returns a vector of calibratable parameter names
  /// @return vector of parameter names
  virtual std::vector<std::string> get_param_names() = 0;

  /// Returns a vector of parameter units
  /// @return vector of parameter units
  virtual std::vector<std::string>  get_param_units() = 0;

  /// Sets a vector containing current values of all calibratible parameters
  /// @return void (pass by reference interface)
  virtual void get_param_values(std::vector<double> &values) = 0;

  /// Returns a vector containing current values of all calibratible parameters
  /// @return std:vector<double> interface
  virtual std::vector<double> get_param_values() = 0;

  ///set all calibratable parameter values from a vector of values
  /// @return int number of parameters set
  virtual int set_param_values(std::vector<double>& values) = 0;

  /// Gets a single parameter value given integer index
  virtual double get_param_value(int& index) = 0;

  /// Sets a single parameter value given index and value
  virtual int set_param_value(int& index, double& value) = 0;

  /// derivative of gibbs free energy \f$G\f$ with respect to parameter given by int index
  virtual double dparam_g(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial G/\partial T\f$ with respect to parameter given by int index
  virtual double dparam_dgdt(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial G/\partial P\f$ with respect to parameter given by int index
  virtual double dparam_dgdp(double& T, double& P, std::vector<double>& n, int& index) = 0;
  
  /// derivative of \f$\partial G/\partial n\f$ (chemical potential) with respect to parameter given by int index
  virtual std::vector<double> dparam_dgdn(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial T^2\f$ with respect to parameter given by int index
  virtual double dparam_d2gdt2(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial T\partial P\f$ with respect to parameter given by int index
  virtual double dparam_d2gdtdp(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial P^2\f$ with respect to parameter given by int index
  virtual double dparam_d2gdp2(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T^3\f$ with respect to parameter given by int index
  virtual double dparam_d3gdt3(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T^2\partial P\f$ with respect to parameter given by int index
  virtual double dparam_d3gdt2dp(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T\partial P^2\f$ with respect to parameter given by int index
  virtual double dparam_d3gdtdp2(double& T, double& P, std::vector<double>& n, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial P^3\f$ with respect to parameter given by int index
  virtual double dparam_d3gdp3(double& T, double& P, std::vector<double>& n, int& index) = 0;                                      
  
  //**************************************************************************
  // ThermoCodegen functions
  //**************************************************************************
  
  /// Vector of molecular weights
  virtual std::vector<double> get_M() const = 0;

  /// Official phase abbreviation
  virtual std::string abbrev() const = 0;

  /// Return pointers to endmembers
  /// @return std::vector of shared pointers to endmember objects
  virtual std::vector<std::shared_ptr<EndMember> > endmembers() const = 0;
  
  /// Vector of chemical potentials for each endmember, given mol fraction x
  virtual std::vector<double> mu(const double& T, const double& P,
                                 std::vector<double>& x) const = 0;

  /// Vector of \f$ \partial\mu/\partial T\f$ for each endmember, given mol fraction x
  virtual std::vector<double> dmu_dT(const double& T, const double& P,
                                     std::vector<double>& x) const = 0;

  /// Vector of \f$ \partial\mu/\partial P\f$ for each endmember, given mol fraction x
  virtual std::vector<double> dmu_dP(const double& T, const double& P,
                                     std::vector<double>& x) const = 0;

  /// Jacobian of compositional derivatives of each chemical potential with respect to each component
  /// \f$\partial \mu_j/\partial\c_k\f$
  ///
  /// Note:  composition is in weight fraction \f$\mathbf{c}\f$ not mol fraction \f$\mathbf{x}\f$
  /// @return vector of vectors (not a matrix per se)
  virtual std::vector<std::vector<double> > dmu_dc(const double& T, const double& P,
                                     std::vector<double>& c) const = 0;

  /// Mass of a mol of phase
  virtual double Mass(std::vector<double>& x) const = 0;
                            
  /// Mean density of the phase for composition \f$ \mathbf{c} \f$ in weight percent
  virtual double rho(const double& T, const double& P,
                     std::vector<double>& c) const = 0;

  /// derivative of phase  density with respect to \f$ T\f$ given composition \f$ \mathbf{c} \f$ in weight percent
  virtual double drho_dT(const double& T, const double& P,
                         std::vector<double>& c) const = 0;

  /// derivative of phase  density with respect to \f$ P\f$ given composition \f$ \mathbf{c} \f$ in weight percent
  virtual double drho_dP(const double& T, const double& P,
                         std::vector<double>& c) const = 0;

  /// derivative of phase  density with respect to  composition \f$ \mathbf{c} \f$ in weight percent
  /// @return a vector containing \f$\partial\rho/\partial c^k\f$
  virtual std::vector<double> drho_dc(const double& T, const double& P,
                 std::vector<double>& c) const = 0;

  /// derivative of entropy with respect to  composition \f$ \mathbf{c} \f$ in weight percent
  /// @return a vector containing \f$\partial s/\partial c^k\f$
  virtual std::vector<double>  ds_dc(const double& T, const double& P,
             std::vector<double>& c) const = 0;

  /// derivative of volume with respect to  composition \f$ \mathbf{c} \f$ in weight percent
  /// @return a vector containing \f$\partial v/\partial c^k\f$
  virtual std::vector<double>  dv_dc(const double& T, const double& P,
             std::vector<double>& c) const = 0;

  /// Convert a vector of concentrations in mass fractions \f$ \mathbf{c} \f$ to mol fractions \f$ \mathbf{x} \f$
  virtual std::vector<double> c_to_x(std::vector<double> &c) const = 0;

  /// Convert a vector of concentrations in mol fractions \f$ \mathbf{x} \f$ to mass fractions \f$ \mathbf{c} \f$
  virtual std::vector<double> x_to_c(std::vector<double> &x) const = 0;
};

#endif
