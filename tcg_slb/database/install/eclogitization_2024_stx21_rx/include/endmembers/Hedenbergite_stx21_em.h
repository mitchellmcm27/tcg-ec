#ifndef HEDENBERGITE_STX21_EM_H
#define HEDENBERGITE_STX21_EM_H

#include "EndMember.h"
#include <map>
#include <string>
#include <vector>

extern "C" {
#include "Hedenbergite_stx21_em_coder.h"
}


class Hedenbergite_stx21_em: public EndMember
{
public:
  // Constructor
  Hedenbergite_stx21_em();

  // Destructor
  virtual ~Hedenbergite_stx21_em();

  //**************************************************************************
  // EndMember data
  //**************************************************************************

  // Identifier name
  std::string identifier();

  // EndMember name
  std::string name();

  // TCG version at build time
  std::string tcg_build_version();

  // TCG git SHA at build time
  std::string tcg_build_git_sha();

  // TCG version at generation time
  std::string tcg_generation_version();

  // TCG git SHA at generation time
  std::string tcg_generation_git_sha();

  // Formula
  std::string formula();

  // Molecular weight
  double molecular_weight();

  // Elements
  std::vector<double> elements();

  //**************************************************************************
  // Gibbs Free Energy G and its derivatives as functions of T and P
  //**************************************************************************

  // G
  double G(const double& T, const double& P);

  // dGdT
  double dGdT(const double& T, const double& P);

  // dGdP
  double dGdP(const double& T, const double& P);

  // d2GdT2
  double d2GdT2(const double& T, const double& P);

  // d2GdTdP
  double d2GdTdP(const double& T, const double& P);

  // d2GdP2
  double d2GdP2(const double& T, const double& P);

  // d3GdT3
  double d3GdT3(const double& T, const double& P);

  // d3GdT2dP
  double d3GdT2dP(const double& T, const double& P);

  // d3GdTdP2
  double d3GdTdP2(const double& T, const double& P);

  // d3GdP3
  double d3GdP3(const double& T, const double& P);

  //**************************************************************************
  // Convenience functions of T and P
  //**************************************************************************

  // Entropy
  double S(const double& T, const double& P);

  // Volume
  double V(const double& T, const double& P);

  // dVdT
  double dVdT(const double& T, const double& P);

  // dVdP
  double dVdP(const double& T, const double& P);

  // Heat capacity at constant volume
  double Cv(const double& T, const double& P);

  // Heat capacity at constant pressure
  double Cp(const double& T, const double& P);

  // dCp/dT (derivative of Cp wrt to temperature T)
  double dCpdT(const double& T, const double& P);

  // Volume thermal expansion coefficient
  double alpha(const double& T, const double& P);

  // Isothermal compressibility
  double beta(const double& T, const double& P);

  // Bulk modulus
  double K(const double& T, const double& P);

  // ...
  double Kp(const double& T, const double& P);

  //**************************************************************************
  // Active parameter functions directly from coder
  //**************************************************************************

  // Number of active parameters
  int get_param_number();

  // Parameter names
  std::vector<std::string> get_param_names();

  // Parameter units
  std::vector<std::string> get_param_units();

  // Return vector of parameters
  std::vector<double> get_param_values();

  // Get parameters (pass by reference)
  void get_param_values(std::vector<double>& values);

  // Set parameters
  int set_param_values(std::vector<double>& values);

  // Get parameter with a given index
  double get_param_value(int& index);

  // Set parameter with a given index
  int set_param_value(int& index, double& value);

  // dparam_g
  double dparam_g(double& T, double& P, int& index);

  // dparam_dgdt
  double dparam_dgdt(double& T, double& P, int& index);

  // dparam_dgdp
  double dparam_dgdp(double& T, double& P, int& index);

  // dparam_d2gdt2
  double dparam_d2gdt2(double& T, double& P, int& index);

  // dparam_d2gdtdp
  double dparam_d2gdtdp(double& T, double& P, int& index);

  // dparam_d2gdp2
  double dparam_d2gdp2(double& T, double& P, int& index);

  // dparam_d3gdt3
  double dparam_d3gdt3(double& T, double& P, int& index);

  // dparam_d3gdt2dp
  double dparam_d3gdt2dp(double& T, double& P, int& index);

  // dparam_d3gdtdp2
  double dparam_d3gdtdp2(double& T, double& P, int& index);

  // dparam_d3gdp3
  double dparam_d3gdp3(double& T, double& P, int& index);

};

#endif
