#ifndef ENDMEMBER_H
#define ENDMEMBER_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

// Abstract class for thermodynamic endmembers
class EndMember
{
public:
  /// Constructor
  EndMember();

  /// Destructor
  virtual ~EndMember();

  //****************************************************************************
  // EndMember data
  //****************************************************************************

  /// Returns a unique identifier for the model instance
  /// @return the identifier
  virtual std::string identifier()= 0;

  /// Returns the endmember name
  /// @return the endmember name
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

  /// Returns the formula
  /// @return the formula
  virtual std::string formula() = 0;

  /// Returns the molecular weight
  /// @return the molecular weight
  virtual double molecular_weight()= 0;

  /// Returns a vector of elements
  /// @return a vector of elements
  virtual std::vector<double> elements() = 0;

  //**************************************************************************
  // Gibbs Free Energy G and its derivatives as functions of T and P
  //****************************************************************************

  /// Returns the Gibbs free energy \f$G\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return G
  virtual double G(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial G/\partial T\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return dGdT
  virtual double dGdT(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial G/\partial P\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return dGdP
  virtual double dGdP(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial T^2\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d2GdT2
  virtual double d2GdT2(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial T \partial P\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d2GdTdP
  virtual double d2GdTdP(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^2 G/\partial P^2\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d2GdP2
  virtual double d2GdP2(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^3 G/\partial T^3\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d3GdT3
  virtual double d3GdT3(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^3 G/\partial T^2 \partial P\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d3GdT2dP
  virtual double d3GdT2dP(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^3 G/\partial T \partial P^2\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d3GdTdP2
  virtual double d3GdTdP2(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial^3 G/\partial P^3\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return d3GdP3
  virtual double d3GdP3(const double& T, const double& P) = 0;

  //**************************************************************************
  // Convenience functions of T and P
  //**************************************************************************

  /// Returns the entropy \f$S\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return S
  virtual double S(const double& T, const double& P) = 0;

  /// Returns the volume \f$V\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return V
  virtual double V(const double& T, const double& P)= 0;

  /// Returns the partial derivative \f$\partial V/\partial T\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return dVdT
  virtual double dVdT(const double& T, const double& P)= 0;

  /// Returns the partial derivative \f$\partial V/\partial P\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return dVdP
  virtual double dVdP(const double& T, const double& P)= 0;

  /// Returns the heat capacity at constant volume \f$C_v\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return Cv
  virtual double Cv(const double& T, const double& P) = 0;

  /// Returns the heat capacity at constant pressure \f$C_p\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return Cp
  virtual double Cp(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial C_p /\partial T\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return dCpdT 
  virtual double dCpdT(const double& T, const double& P) = 0;

  /// Returns the volume thermal expansion coefficient \f$\alpha\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return alpha 
  virtual double alpha(const double& T, const double& P) = 0;

  /// Returns the isothermal compressibility \f$\beta\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return beta
  virtual double beta(const double& T, const double& P) = 0;

  /// Returns the bulk modulus \f$K\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return K
  virtual double K(const double& T, const double& P) = 0;

  /// Returns the partial derivative \f$\partial K/\partial P\f$
  /// @param[in] T Temperature
  /// @param[in] P Pressure
  /// @return Kp
  virtual double Kp(const double& T, const double& P) = 0;

  //**************************************************************************
  // Active parameter functions directly from coder
  //**************************************************************************

  /// Returns the number of parameters
  /// @return number of calibratable parameters in endmember model
  virtual int get_param_number() = 0;

  /// Returns a vector of parameter names
  /// @return vector of parameter names
  virtual std::vector<std::string> get_param_names() = 0;

  /// Returns a vector of parameter units
  /// @return vector of parameter units
  virtual std::vector<std::string>  get_param_units() = 0;

  /// load all parameter values into a vector
  virtual void get_param_values(std::vector<double> &values) = 0;

  /// Returns a vector containing current values of all calibratible parameters
  /// @return vector of parameter values
  virtual std::vector<double> get_param_values() = 0;

  ///set all parameter values from a vector of values
  /// @return int number of parameters set
  virtual int set_param_values(std::vector<double>& values) = 0;

  /// Gets a single parameter value given integer index
  virtual double get_param_value(int& index) = 0;

  /// Sets a single parameter value given index and value
  virtual int set_param_value(int& index, double& value) = 0;

  /// derivative of gibbs free energy \f$G\f$ with respect to parameter given by int index
  virtual double dparam_g(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial G/\partial T\f$ with respect to parameter given by int index
  virtual double dparam_dgdt(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial G/\partial P\f$ with respect to parameter given by int index
  virtual double dparam_dgdp(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial T^2\f$ with respect to parameter given by int index
  virtual double dparam_d2gdt2(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial T\partial P\f$ with respect to parameter given by int index
  virtual double dparam_d2gdtdp(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^2 G/\partial P^2\f$ with respect to parameter given by int index
  virtual double dparam_d2gdp2(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T^3\f$ with respect to parameter given by int index
  virtual double dparam_d3gdt3(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T^2\partial P\f$ with respect to parameter given by int index
  virtual double dparam_d3gdt2dp(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial T\partial P^2\f$ with respect to parameter given by int index
  virtual double dparam_d3gdtdp2(double& T, double& P, int& index) = 0;

  /// derivative of \f$\partial^3 G/\partial P^3\f$ with respect to parameter given by int index
  virtual double dparam_d3gdp3(double& T, double& P, int& index) = 0;

};

#endif
