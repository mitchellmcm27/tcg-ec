#ifndef CORUNDUM_SLB_PH_H
#define CORUNDUM_SLB_PH_H

#include "Phase.h"
#include "endmembers.h"

extern "C" {
#include "Corundum_slb_ph_coder.h"
}

class Corundum_slb_ph: public Phase
{
public:
  // Constructor
  Corundum_slb_ph();

  // Destructor
  virtual ~Corundum_slb_ph();
  
  //**************************************************************************
  // Coder functions
  //**************************************************************************

  // identifier
  std::string identifier();

  // name
  std::string name();

  // TCG version at build time
  std::string tcg_build_version();

  // TCG git SHA at build time
  std::string tcg_build_git_sha();

  // TCG version at generation time
  std::string tcg_generation_version();

  // TCG git SHA at generation time
  std::string tcg_generation_git_sha();

  // formula
  std::string formula(const double& T, const double& P, std::vector<double>& n);

  // conv_elm_to_moles
  std::vector<double> conv_elm_to_moles(std::vector<double>& e);

  // Extra Coder functions provided by Cython in ThermoEngine and written in C++ here
  //conv_elm_to_tot_moles
  double conv_elm_to_tot_moles(std::vector<double>& e);

  //conv_elm_to_tot_grams
  double conv_elm_to_tot_grams(std::vector<double>& e);

  //conv_moles_to_elm
  std::vector<double> conv_moles_to_elm(std::vector<double>& n);

  //conv_moles_to_tot_moles
  double conv_moles_to_tot_moles(std::vector<double>& n);

  //conv_moles_to_mole_frac
  std::vector<double> conv_moles_to_mole_frac(std::vector<double>& n);

 // end Extra functions

  // test_moles
  int test_moles(std::vector<double>& n);

  // endmember_number
  int endmember_number() const;
  
  // endmember_name
  std::string endmember_name(const int &i) const;
  
  // endmember_formula
  std::string endmember_formula(const int &i) const;
  
  // endmember_mw
  double endmember_mw(const int &i) const;
  
  // endmember_elements
  std::vector<double> endmember_elements(const int &i) const;

  // species_number
  int species_number() const;
  
  // species_name
  std::string species_name(const int &i) const;
  
  // species_formula
  std::string species_formula(const int &i) const;
  
  // species_mw
  double species_mw(const int &i) const;
  
  // species_elements
  std::vector<double> species_elements(const int &i) const;

  // g
  double g(const double& T, const double& P, std::vector<double>& n) const;

  // dgdt
  double dgdt(const double& T, const double& P, std::vector<double>& n) const;

  // dgdp
  double dgdp(const double& T, const double& P, std::vector<double>& n) const;

  // d2gdt2
  double d2gdt2(const double& T, const double& P, std::vector<double>& n) const;

  // d2gdtdp
  double d2gdtdp(const double& T, const double& P,
                 std::vector<double>& n) const;

  // d2gdp2
  double d2gdp2(const double& T, const double& P, std::vector<double>& n) const;
  
  // d3gdt3
  double d3gdt3(const double& T, const double& P, std::vector<double>& n) const;
  
  // d3gdt2dp
  double d3gdt2dp(const double& T, const double& P, 
                  std::vector<double>& n) const;
  
  // d3gdtdp2
  double d3gdtdp2(const double& T, const double& P,
                  std::vector<double>& n) const;
                  
  // d3gdp3
  double d3gdp3(const double& T, const double& P,
                    std::vector<double>& n) const;
  
  // dgdn
  std::vector<double> dgdn(const double& T, const double& P,
                           std::vector<double>& n) const;

  // d2gdndt
  std::vector<double> d2gdndt(const double& T, const double& P,
                              std::vector<double>& n) const;

  // d2gdndp
  std::vector<double> d2gdndp(const double& T, const double& P,
                              std::vector<double>& n) const; 
                              
  // d2gdn2
  std::vector<double> d2gdn2(const double& T, const double& P,
                             std::vector<double>& n) const;
                             
  // d3gdndt2 
  std::vector<double> d3gdndt2(const double& T, const double& P,
                               std::vector<double>& n) const;
                               
  // d3gdndtdp 
  std::vector<double> d3gdndtdp(const double& T, const double& P,
                                std::vector<double>& n) const;
                                
  // d3gdn2dt 
  std::vector<double> d3gdn2dt(const double& T, const double& P,
                               std::vector<double>& n) const;
                               
  // d3gdndp2 
  std::vector<double> d3gdndp2(const double& T, const double& P,
                               std::vector<double>& n) const;
                               
  // d3gdn2dp 
  std::vector<double> d3gdn2dp(const double& T, const double& P,
                               std::vector<double>& n) const;  
                               
  // d3gdn3 
  std::vector<double> d3gdn3(const double& T, const double& P,
                             std::vector<double>& n) const;

  // v
  double v(const double& T, const double& P,
           std::vector<double>& n) const;

  // s
  double s(const double& T, const double& P,
           std::vector<double>& n) const;

  // alpha
  double alpha(const double& T, const double& P,
               std::vector<double>& n) const;

  // cv
  double cv(const double& T, const double& P,
            std::vector<double>& n) const;
                           
  // cp
  double cp(const double& T, const double& P,
            std::vector<double>& n) const;

  // dcpdt
  double dcpdt(const double& T, const double& P,
               std::vector<double>& n) const;

  // beta
  double beta(const double& T, const double& P,
              std::vector<double>& n) const;

  // K
  double K(const double& T, const double& P,
           std::vector<double>& n) const;
                         
  // Kp
  double Kp(const double& T, const double& P,
            std::vector<double>& n) const;

  // get_param_number
  int get_param_number();

  // get_param_names
  std::vector<std::string> get_param_names();

  // get_param_units
  std::vector<std::string> get_param_units();

  // get_param_values
  std::vector<double> get_param_values();

  // get_param_values
  void get_param_values(std::vector<double>& values);

  // set_param_values
  int set_param_values(std::vector<double>& values);

  // get_param_value
  double get_param_value(int& index);

  // set_param_value
  int set_param_value(int& index, double& value);

  // dparam_g
  double dparam_g(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_dgdt
  double dparam_dgdt(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_dgdp
  double dparam_dgdp(double& T, double& P, std::vector<double>& n, int& index);
  
  // dparam_dgdn
  std::vector<double> dparam_dgdn(double& T, double& P, 
                                  std::vector<double>& n, int& index);

  // dparam_d2gdt2
  double dparam_d2gdt2(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d2gdtdp
  double dparam_d2gdtdp(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d2gdp2
  double dparam_d2gdp2(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d3gdt3
  double dparam_d3gdt3(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d3gdt2dp
  double dparam_d3gdt2dp(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d3gdtdp2
  double dparam_d3gdtdp2(double& T, double& P, std::vector<double>& n, int& index);

  // dparam_d3gdp3
  double dparam_d3gdp3(double& T, double& P, std::vector<double>& n, int& index);

  //**************************************************************************
  // ThermoCodegen functions
  //**************************************************************************

  // Vector of molecular weights
  std::vector<double> get_M() const;
  
  // Official phase abbreviation
  std::string abbrev() const;

  //return pointers to endmembers
  std::vector<std::shared_ptr<EndMember> > endmembers() const;

  // Vector of chemical potentials for each component
  std::vector<double> mu(const double& T, const double& P,
                         std::vector<double>& x) const;

  // dmu_dT
  std::vector<double> dmu_dT(const double& T, const double& P,
                             std::vector<double>& x) const;

  // dmu_dP
  std::vector<double> dmu_dP(const double& T, const double& P,
                             std::vector<double>& x) const;

  // dmu_dc (Jacobian of change of chemical potentials for a change in concentration)
  std::vector<std::vector<double> > dmu_dc(const double& T, const double& P,
                              std::vector<double>& c) const;

  // Molar Mass
  double Mass(std::vector<double>& x) const;

  // Mean density of the phase in mass units
  double rho(const double& T, const double& P,
             std::vector<double>& c) const;

  // drho_dT
  double drho_dT(const double& T, const double& P,
                 std::vector<double>& c) const;

  // drho_dP
  double drho_dP(const double& T, const double& P,
                 std::vector<double>& c) const;

  // drho_dC corrected derivative of density with respect to composition c
  std::vector<double> drho_dc(const double& T, const double& P,
                 std::vector<double>& c) const;

  //derivative of entropy with respect to composition c
  std::vector<double>  ds_dc(const double& T, const double& P,
             std::vector<double>& c) const;

  //derivative of volume with respect to composition c
  std::vector<double>  dv_dc(const double& T, const double& P,
             std::vector<double>& c) const;


  // Returns weight percent concentrations given mole fraction
  std::vector<double> c_to_x(std::vector<double> &c) const;

  // Returns  mole fraction given weight percent concentrations
  std::vector<double> x_to_c(std::vector<double> &x) const;

  // pass by reference versions
  void c_to_x(std::vector<double> &c, std::vector<double> &x) const;
  void x_to_c(std::vector<double> &x, std::vector<double> &c) const;

private:
  // Vector of pointers to endmembers
  std::vector<std::shared_ptr<EndMember> > _endmembers;

  // Number of components
  int C;

  // Vector of molecular masses
  std::vector<double> M;

  // Vector of inverse molecular masses
  std::vector<double> iM;

  // Temporary vectors
  mutable std::vector<double> _tmp;
  mutable std::vector<double> _result;
  mutable std::vector<double> _result_mat;
  mutable std::vector<double> _elm;
  mutable std::vector<std::vector<double> > _dmu_dc;



};

#endif
