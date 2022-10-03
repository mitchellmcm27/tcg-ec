#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include "tcgversion.h"
#include "endmembers/AlAkimotoite_slb_em.h"
#include "endmembers/Albite_slb_em.h"
#include "endmembers/Almandine_slb_em.h"
#include "endmembers/AlPerovskite_slb_em.h"
#include "endmembers/AlPostPerovskite_slb_em.h"
#include "endmembers/Anorthite_slb_em.h"
#include "endmembers/CaPerovskite_slb_em.h"
#include "endmembers/CaTschermaks_slb_em.h"
#include "endmembers/Clinoenstatite_slb_em.h"
#include "endmembers/Coesite_slb_em.h"
#include "endmembers/Diopside_slb_em.h"
#include "endmembers/Enstatite_slb_em.h"
#include "endmembers/Fayalite_slb_em.h"
#include "endmembers/FeAkimotoite_slb_em.h"
#include "endmembers/FeCaFerrite_slb_em.h"
#include "endmembers/FePerovskite_slb_em.h"
#include "endmembers/FePostPerovskite_slb_em.h"
#include "endmembers/FeRingwoodite_slb_em.h"
#include "endmembers/Ferrosilite_slb_em.h"
#include "endmembers/FeWadsleyite_slb_em.h"
#include "endmembers/Forsterite_slb_em.h"
#include "endmembers/Grossular_slb_em.h"
#include "endmembers/Hedenbergite_slb_em.h"
#include "endmembers/Hercynite_slb_em.h"
#include "endmembers/HPClinoenstatite_slb_em.h"
#include "endmembers/HPClinoferrosilite_slb_em.h"
#include "endmembers/Jadeite_slb_em.h"
#include "endmembers/Kyanite_slb_em.h"
#include "endmembers/MgAkimotoite_slb_em.h"
#include "endmembers/MgCaFerrite_slb_em.h"
#include "endmembers/MgMajorite_slb_em.h"
#include "endmembers/MgPerovskite_slb_em.h"
#include "endmembers/MgPostPerovskite_slb_em.h"
#include "endmembers/MgRingwoodite_slb_em.h"
#include "endmembers/MgSpinel_slb_em.h"
#include "endmembers/MgTschermaks_slb_em.h"
#include "endmembers/MgWadsleyite_slb_em.h"
#include "endmembers/NaCaFerrite_slb_em.h"
#include "endmembers/NaMajorite_slb_em.h"
#include "endmembers/Nepheline_slb_em.h"
#include "endmembers/OrthoDiopside_slb_em.h"
#include "endmembers/Periclase_slb_em.h"
#include "endmembers/Pyrope_slb_em.h"
#include "endmembers/Quartz_slb_em.h"
#include "endmembers/Seifertite_slb_em.h"
#include "endmembers/Stishovite_slb_em.h"
#include "endmembers/Wuestite_slb_em.h"
#include "phases/Akimotoite_slb_ph.h"
#include "phases/AlAkimotoite_slb_ph.h"
#include "phases/Albite_slb_ph.h"
#include "phases/Almandine_slb_ph.h"
#include "phases/AlPerovskite_slb_ph.h"
#include "phases/AlPostPerovskite_slb_ph.h"
#include "phases/Anorthite_slb_ph.h"
#include "phases/CaFerritePhase_slb_ph.h"
#include "phases/CaPerovskite_slb_ph.h"
#include "phases/CaTschermaks_slb_ph.h"
#include "phases/Clinoenstatite_slb_ph.h"
#include "phases/Clinopyroxene_slb_ph.h"
#include "phases/Coesite_slb_ph.h"
#include "phases/Diopside_slb_ph.h"
#include "phases/Enstatite_slb_ph.h"
#include "phases/Fayalite_slb_ph.h"
#include "phases/FeAkimotoite_slb_ph.h"
#include "phases/FeCaFerrite_slb_ph.h"
#include "phases/Feldspar_slb_ph.h"
#include "phases/FePerovskite_slb_ph.h"
#include "phases/FePostPerovskite_slb_ph.h"
#include "phases/FeRingwoodite_slb_ph.h"
#include "phases/Ferrosilite_slb_ph.h"
#include "phases/FeWadsleyite_slb_ph.h"
#include "phases/Forsterite_slb_ph.h"
#include "phases/Garnet_slb_ph.h"
#include "phases/Grossular_slb_ph.h"
#include "phases/Hedenbergite_slb_ph.h"
#include "phases/Hercynite_slb_ph.h"
#include "phases/HPClinoenstatite_slb_ph.h"
#include "phases/HPClinoferrosilite_slb_ph.h"
#include "phases/HPClinopyroxene_slb_ph.h"
#include "phases/Jadeite_slb_ph.h"
#include "phases/Kyanite_slb_ph.h"
#include "phases/Magnesiowuestite_slb_ph.h"
#include "phases/MgAkimotoite_slb_ph.h"
#include "phases/MgCaFerrite_slb_ph.h"
#include "phases/MgFePerovskite_slb_ph.h"
#include "phases/MgMajorite_slb_ph.h"
#include "phases/MgPerovskite_slb_ph.h"
#include "phases/MgPostPerovskite_slb_ph.h"
#include "phases/MgRingwoodite_slb_ph.h"
#include "phases/MgSpinel_slb_ph.h"
#include "phases/MgTschermaks_slb_ph.h"
#include "phases/MgWadsleyite_slb_ph.h"
#include "phases/NaCaFerrite_slb_ph.h"
#include "phases/NaMajorite_slb_ph.h"
#include "phases/Nepheline_slb_ph.h"
#include "phases/Olivine_slb_ph.h"
#include "phases/OrthoDiopside_slb_ph.h"
#include "phases/Orthopyroxene_slb_ph.h"
#include "phases/Periclase_slb_ph.h"
#include "phases/Perovskite_slb_ph.h"
#include "phases/PostPerovskite_slb_ph.h"
#include "phases/Pyrope_slb_ph.h"
#include "phases/Quartz_slb_ph.h"
#include "phases/Ringwoodite_slb_ph.h"
#include "phases/Seifertite_slb_ph.h"
#include "phases/Spinel_slb_ph.h"
#include "phases/Stishovite_slb_ph.h"
#include "phases/Wadsleyite_slb_ph.h"
#include "phases/Wuestite_slb_ph.h"
#include "reactions/eclogite_slb_rx.h"


namespace py = pybind11;

PYBIND11_MODULE(py_eclogite_slb_rx, m) {
  //Return the TCG version at build time
  m.def("tcg_build_version",[]() {
    return TCG_VERSION;
  }, "Return the TCG version at build time");

  //Return the TCG git sha at build time
  m.def("tcg_build_git_sha",[]() {
    return TCG_GIT_SHA;
  }, "Return the TCG git sha at build time");

  //Return the TCG version at generation time
  m.def("tcg_generation_version",[]() {
    return "0.6.9";
  }, "Return the TCG version at generation time");

  //Return the TCG git sha at generation time
  m.def("tcg_generation_git_sha",[]() {
    return "117d758197e2d445579ad671f573064c0650429d Mon Aug 1 00:36:24 2022 +0000";
  }, "Return the TCG git sha at generation time");

  // Returns a dictionary of phase attributes
  m.def("phase_info",[]() {
    std::vector<std::shared_ptr<Phase> > _phases = {std::make_shared<Akimotoite_slb_ph>(),std::make_shared<AlAkimotoite_slb_ph>(),std::make_shared<Albite_slb_ph>(),std::make_shared<Almandine_slb_ph>(),std::make_shared<AlPerovskite_slb_ph>(),std::make_shared<AlPostPerovskite_slb_ph>(),std::make_shared<Anorthite_slb_ph>(),std::make_shared<CaFerritePhase_slb_ph>(),std::make_shared<CaPerovskite_slb_ph>(),std::make_shared<CaTschermaks_slb_ph>(),std::make_shared<Clinoenstatite_slb_ph>(),std::make_shared<Clinopyroxene_slb_ph>(),std::make_shared<Coesite_slb_ph>(),std::make_shared<Diopside_slb_ph>(),std::make_shared<Enstatite_slb_ph>(),std::make_shared<Fayalite_slb_ph>(),std::make_shared<FeAkimotoite_slb_ph>(),std::make_shared<FeCaFerrite_slb_ph>(),std::make_shared<Feldspar_slb_ph>(),std::make_shared<FePerovskite_slb_ph>(),std::make_shared<FePostPerovskite_slb_ph>(),std::make_shared<FeRingwoodite_slb_ph>(),std::make_shared<Ferrosilite_slb_ph>(),std::make_shared<FeWadsleyite_slb_ph>(),std::make_shared<Forsterite_slb_ph>(),std::make_shared<Garnet_slb_ph>(),std::make_shared<Grossular_slb_ph>(),std::make_shared<Hedenbergite_slb_ph>(),std::make_shared<Hercynite_slb_ph>(),std::make_shared<HPClinoenstatite_slb_ph>(),std::make_shared<HPClinoferrosilite_slb_ph>(),std::make_shared<HPClinopyroxene_slb_ph>(),std::make_shared<Jadeite_slb_ph>(),std::make_shared<Kyanite_slb_ph>(),std::make_shared<Magnesiowuestite_slb_ph>(),std::make_shared<MgAkimotoite_slb_ph>(),std::make_shared<MgCaFerrite_slb_ph>(),std::make_shared<MgFePerovskite_slb_ph>(),std::make_shared<MgMajorite_slb_ph>(),std::make_shared<MgPerovskite_slb_ph>(),std::make_shared<MgPostPerovskite_slb_ph>(),std::make_shared<MgRingwoodite_slb_ph>(),std::make_shared<MgSpinel_slb_ph>(),std::make_shared<MgTschermaks_slb_ph>(),std::make_shared<MgWadsleyite_slb_ph>(),std::make_shared<NaCaFerrite_slb_ph>(),std::make_shared<NaMajorite_slb_ph>(),std::make_shared<Nepheline_slb_ph>(),std::make_shared<Olivine_slb_ph>(),std::make_shared<OrthoDiopside_slb_ph>(),std::make_shared<Orthopyroxene_slb_ph>(),std::make_shared<Periclase_slb_ph>(),std::make_shared<Perovskite_slb_ph>(),std::make_shared<PostPerovskite_slb_ph>(),std::make_shared<Pyrope_slb_ph>(),std::make_shared<Quartz_slb_ph>(),std::make_shared<Ringwoodite_slb_ph>(),std::make_shared<Seifertite_slb_ph>(),std::make_shared<Spinel_slb_ph>(),std::make_shared<Stishovite_slb_ph>(),std::make_shared<Wadsleyite_slb_ph>(),std::make_shared<Wuestite_slb_ph>()};
    std::vector<std::string> _names, _abbrevs, _types;
    
    for(int i = 0; i < _phases.size(); i++){
      _names.push_back(_phases[i]->name());
      _abbrevs.push_back(_phases[i]->abbrev());
      
      if(_phases[i]->endmember_number() == 1){
        _types.push_back("pure");
      }
      else{
        _types.push_back("solution");
      }
    }
    
    std::map<std::string, std::vector<std::string> > phase_info;    
    phase_info.insert(std::make_pair("ClassName", _names));
    phase_info.insert(std::make_pair("Abbrev", _abbrevs));
    phase_info.insert(std::make_pair("PhaseType", _types));
  
    return phase_info; 
  });

  // EndMembers
  py::class_<AlAkimotoite_slb_em, std::shared_ptr<AlAkimotoite_slb_em> >(m, "AlAkimotoite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &AlAkimotoite_slb_em::name, "Name of endmember")
    .def("identifier", &AlAkimotoite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &AlAkimotoite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlAkimotoite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlAkimotoite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlAkimotoite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlAkimotoite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &AlAkimotoite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &AlAkimotoite_slb_em::elements, "Vector of elements")
    
    .def("g", &AlAkimotoite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &AlAkimotoite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &AlAkimotoite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &AlAkimotoite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &AlAkimotoite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &AlAkimotoite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &AlAkimotoite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &AlAkimotoite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &AlAkimotoite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &AlAkimotoite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &AlAkimotoite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &AlAkimotoite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &AlAkimotoite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &AlAkimotoite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &AlAkimotoite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &AlAkimotoite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &AlAkimotoite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &AlAkimotoite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &AlAkimotoite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &AlAkimotoite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &AlAkimotoite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &AlAkimotoite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlAkimotoite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &AlAkimotoite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &AlAkimotoite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &AlAkimotoite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &AlAkimotoite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &AlAkimotoite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &AlAkimotoite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &AlAkimotoite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &AlAkimotoite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &AlAkimotoite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &AlAkimotoite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &AlAkimotoite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &AlAkimotoite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &AlAkimotoite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Albite_slb_em, std::shared_ptr<Albite_slb_em> >(m, "Albite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Albite_slb_em::name, "Name of endmember")
    .def("identifier", &Albite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Albite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Albite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Albite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Albite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Albite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Albite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Albite_slb_em::elements, "Vector of elements")
    
    .def("g", &Albite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Albite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Albite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Albite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Albite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Albite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Albite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Albite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Albite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Albite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Albite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Albite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Albite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Albite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Albite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Albite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Albite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Albite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Albite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Albite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Albite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Albite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Albite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Albite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Albite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Albite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Albite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Albite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Albite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Albite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Albite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Albite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Albite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Albite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Albite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Albite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Almandine_slb_em, std::shared_ptr<Almandine_slb_em> >(m, "Almandine_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Almandine_slb_em::name, "Name of endmember")
    .def("identifier", &Almandine_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Almandine_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Almandine_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Almandine_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Almandine_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Almandine_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Almandine_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Almandine_slb_em::elements, "Vector of elements")
    
    .def("g", &Almandine_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Almandine_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Almandine_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Almandine_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Almandine_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Almandine_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Almandine_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Almandine_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Almandine_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Almandine_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Almandine_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Almandine_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Almandine_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Almandine_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Almandine_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Almandine_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Almandine_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Almandine_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Almandine_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Almandine_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Almandine_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Almandine_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Almandine_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Almandine_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Almandine_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Almandine_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Almandine_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Almandine_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Almandine_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Almandine_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Almandine_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Almandine_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Almandine_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Almandine_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Almandine_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Almandine_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<AlPerovskite_slb_em, std::shared_ptr<AlPerovskite_slb_em> >(m, "AlPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &AlPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &AlPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &AlPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &AlPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &AlPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &AlPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &AlPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &AlPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &AlPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &AlPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &AlPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &AlPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &AlPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &AlPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &AlPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &AlPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &AlPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &AlPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &AlPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &AlPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &AlPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &AlPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &AlPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &AlPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &AlPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &AlPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &AlPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &AlPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &AlPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &AlPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &AlPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &AlPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &AlPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &AlPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &AlPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &AlPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &AlPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &AlPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &AlPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &AlPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<AlPostPerovskite_slb_em, std::shared_ptr<AlPostPerovskite_slb_em> >(m, "AlPostPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &AlPostPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &AlPostPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &AlPostPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlPostPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlPostPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlPostPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlPostPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &AlPostPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &AlPostPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &AlPostPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &AlPostPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &AlPostPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &AlPostPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &AlPostPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &AlPostPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &AlPostPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &AlPostPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &AlPostPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &AlPostPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &AlPostPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &AlPostPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &AlPostPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &AlPostPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &AlPostPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &AlPostPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &AlPostPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &AlPostPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &AlPostPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &AlPostPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &AlPostPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &AlPostPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlPostPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &AlPostPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &AlPostPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &AlPostPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &AlPostPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &AlPostPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &AlPostPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &AlPostPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &AlPostPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &AlPostPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &AlPostPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &AlPostPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &AlPostPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &AlPostPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Anorthite_slb_em, std::shared_ptr<Anorthite_slb_em> >(m, "Anorthite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Anorthite_slb_em::name, "Name of endmember")
    .def("identifier", &Anorthite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Anorthite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Anorthite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Anorthite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Anorthite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Anorthite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Anorthite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Anorthite_slb_em::elements, "Vector of elements")
    
    .def("g", &Anorthite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Anorthite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Anorthite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Anorthite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Anorthite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Anorthite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Anorthite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Anorthite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Anorthite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Anorthite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Anorthite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Anorthite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Anorthite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Anorthite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Anorthite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Anorthite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Anorthite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Anorthite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Anorthite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Anorthite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Anorthite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Anorthite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Anorthite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Anorthite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Anorthite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Anorthite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Anorthite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Anorthite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Anorthite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Anorthite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Anorthite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Anorthite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Anorthite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Anorthite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Anorthite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Anorthite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<CaPerovskite_slb_em, std::shared_ptr<CaPerovskite_slb_em> >(m, "CaPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &CaPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &CaPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &CaPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &CaPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &CaPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &CaPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &CaPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &CaPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &CaPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &CaPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &CaPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &CaPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &CaPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &CaPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &CaPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &CaPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &CaPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &CaPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &CaPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &CaPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &CaPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &CaPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &CaPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &CaPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &CaPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &CaPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &CaPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &CaPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &CaPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &CaPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &CaPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&CaPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &CaPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &CaPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &CaPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &CaPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &CaPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &CaPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &CaPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &CaPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &CaPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &CaPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &CaPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &CaPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &CaPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<CaTschermaks_slb_em, std::shared_ptr<CaTschermaks_slb_em> >(m, "CaTschermaks_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &CaTschermaks_slb_em::name, "Name of endmember")
    .def("identifier", &CaTschermaks_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &CaTschermaks_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &CaTschermaks_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &CaTschermaks_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &CaTschermaks_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &CaTschermaks_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &CaTschermaks_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &CaTschermaks_slb_em::elements, "Vector of elements")
    
    .def("g", &CaTschermaks_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &CaTschermaks_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &CaTschermaks_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &CaTschermaks_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &CaTschermaks_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &CaTschermaks_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &CaTschermaks_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &CaTschermaks_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &CaTschermaks_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &CaTschermaks_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &CaTschermaks_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &CaTschermaks_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &CaTschermaks_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &CaTschermaks_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &CaTschermaks_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &CaTschermaks_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &CaTschermaks_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &CaTschermaks_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &CaTschermaks_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &CaTschermaks_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &CaTschermaks_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &CaTschermaks_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&CaTschermaks_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &CaTschermaks_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &CaTschermaks_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &CaTschermaks_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &CaTschermaks_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &CaTschermaks_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &CaTschermaks_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &CaTschermaks_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &CaTschermaks_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &CaTschermaks_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &CaTschermaks_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &CaTschermaks_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &CaTschermaks_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &CaTschermaks_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Clinoenstatite_slb_em, std::shared_ptr<Clinoenstatite_slb_em> >(m, "Clinoenstatite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Clinoenstatite_slb_em::name, "Name of endmember")
    .def("identifier", &Clinoenstatite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Clinoenstatite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Clinoenstatite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Clinoenstatite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Clinoenstatite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Clinoenstatite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Clinoenstatite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Clinoenstatite_slb_em::elements, "Vector of elements")
    
    .def("g", &Clinoenstatite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Clinoenstatite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Clinoenstatite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Clinoenstatite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Clinoenstatite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Clinoenstatite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Clinoenstatite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Clinoenstatite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Clinoenstatite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Clinoenstatite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Clinoenstatite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Clinoenstatite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Clinoenstatite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Clinoenstatite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Clinoenstatite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Clinoenstatite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Clinoenstatite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Clinoenstatite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Clinoenstatite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Clinoenstatite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Clinoenstatite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Clinoenstatite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Clinoenstatite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Clinoenstatite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Clinoenstatite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Clinoenstatite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Clinoenstatite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Clinoenstatite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Clinoenstatite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Clinoenstatite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Clinoenstatite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Clinoenstatite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Clinoenstatite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Clinoenstatite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Clinoenstatite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Clinoenstatite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Coesite_slb_em, std::shared_ptr<Coesite_slb_em> >(m, "Coesite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Coesite_slb_em::name, "Name of endmember")
    .def("identifier", &Coesite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Coesite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Coesite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Coesite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Coesite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Coesite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Coesite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Coesite_slb_em::elements, "Vector of elements")
    
    .def("g", &Coesite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Coesite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Coesite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Coesite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Coesite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Coesite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Coesite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Coesite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Coesite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Coesite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Coesite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Coesite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Coesite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Coesite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Coesite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Coesite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Coesite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Coesite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Coesite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Coesite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Coesite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Coesite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Coesite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Coesite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Coesite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Coesite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Coesite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Coesite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Coesite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Coesite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Coesite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Coesite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Coesite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Coesite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Coesite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Coesite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Diopside_slb_em, std::shared_ptr<Diopside_slb_em> >(m, "Diopside_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Diopside_slb_em::name, "Name of endmember")
    .def("identifier", &Diopside_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Diopside_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Diopside_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Diopside_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Diopside_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Diopside_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Diopside_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Diopside_slb_em::elements, "Vector of elements")
    
    .def("g", &Diopside_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Diopside_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Diopside_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Diopside_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Diopside_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Diopside_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Diopside_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Diopside_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Diopside_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Diopside_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Diopside_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Diopside_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Diopside_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Diopside_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Diopside_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Diopside_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Diopside_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Diopside_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Diopside_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Diopside_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Diopside_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Diopside_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Diopside_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Diopside_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Diopside_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Diopside_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Diopside_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Diopside_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Diopside_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Diopside_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Diopside_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Diopside_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Diopside_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Diopside_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Diopside_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Diopside_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Enstatite_slb_em, std::shared_ptr<Enstatite_slb_em> >(m, "Enstatite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Enstatite_slb_em::name, "Name of endmember")
    .def("identifier", &Enstatite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Enstatite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Enstatite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Enstatite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Enstatite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Enstatite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Enstatite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Enstatite_slb_em::elements, "Vector of elements")
    
    .def("g", &Enstatite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Enstatite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Enstatite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Enstatite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Enstatite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Enstatite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Enstatite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Enstatite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Enstatite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Enstatite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Enstatite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Enstatite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Enstatite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Enstatite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Enstatite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Enstatite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Enstatite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Enstatite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Enstatite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Enstatite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Enstatite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Enstatite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Enstatite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Enstatite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Enstatite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Enstatite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Enstatite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Enstatite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Enstatite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Enstatite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Enstatite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Enstatite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Enstatite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Enstatite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Enstatite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Enstatite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Fayalite_slb_em, std::shared_ptr<Fayalite_slb_em> >(m, "Fayalite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Fayalite_slb_em::name, "Name of endmember")
    .def("identifier", &Fayalite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Fayalite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Fayalite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Fayalite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Fayalite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Fayalite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Fayalite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Fayalite_slb_em::elements, "Vector of elements")
    
    .def("g", &Fayalite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Fayalite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Fayalite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Fayalite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Fayalite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Fayalite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Fayalite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Fayalite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Fayalite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Fayalite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Fayalite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Fayalite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Fayalite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Fayalite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Fayalite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Fayalite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Fayalite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Fayalite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Fayalite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Fayalite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Fayalite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Fayalite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Fayalite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Fayalite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Fayalite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Fayalite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Fayalite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Fayalite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Fayalite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Fayalite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Fayalite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Fayalite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Fayalite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Fayalite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Fayalite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Fayalite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FeAkimotoite_slb_em, std::shared_ptr<FeAkimotoite_slb_em> >(m, "FeAkimotoite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FeAkimotoite_slb_em::name, "Name of endmember")
    .def("identifier", &FeAkimotoite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FeAkimotoite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeAkimotoite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeAkimotoite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeAkimotoite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeAkimotoite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FeAkimotoite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FeAkimotoite_slb_em::elements, "Vector of elements")
    
    .def("g", &FeAkimotoite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FeAkimotoite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FeAkimotoite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FeAkimotoite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FeAkimotoite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FeAkimotoite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FeAkimotoite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FeAkimotoite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FeAkimotoite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FeAkimotoite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FeAkimotoite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FeAkimotoite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FeAkimotoite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FeAkimotoite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FeAkimotoite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FeAkimotoite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FeAkimotoite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FeAkimotoite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FeAkimotoite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FeAkimotoite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FeAkimotoite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FeAkimotoite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeAkimotoite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FeAkimotoite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FeAkimotoite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FeAkimotoite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FeAkimotoite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FeAkimotoite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FeAkimotoite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FeAkimotoite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FeAkimotoite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FeAkimotoite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FeAkimotoite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FeAkimotoite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FeAkimotoite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FeAkimotoite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FeCaFerrite_slb_em, std::shared_ptr<FeCaFerrite_slb_em> >(m, "FeCaFerrite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FeCaFerrite_slb_em::name, "Name of endmember")
    .def("identifier", &FeCaFerrite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FeCaFerrite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeCaFerrite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeCaFerrite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeCaFerrite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeCaFerrite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FeCaFerrite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FeCaFerrite_slb_em::elements, "Vector of elements")
    
    .def("g", &FeCaFerrite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FeCaFerrite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FeCaFerrite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FeCaFerrite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FeCaFerrite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FeCaFerrite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FeCaFerrite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FeCaFerrite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FeCaFerrite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FeCaFerrite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FeCaFerrite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FeCaFerrite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FeCaFerrite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FeCaFerrite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FeCaFerrite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FeCaFerrite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FeCaFerrite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FeCaFerrite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FeCaFerrite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FeCaFerrite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FeCaFerrite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FeCaFerrite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeCaFerrite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FeCaFerrite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FeCaFerrite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FeCaFerrite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FeCaFerrite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FeCaFerrite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FeCaFerrite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FeCaFerrite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FeCaFerrite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FeCaFerrite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FeCaFerrite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FeCaFerrite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FeCaFerrite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FeCaFerrite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FePerovskite_slb_em, std::shared_ptr<FePerovskite_slb_em> >(m, "FePerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FePerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &FePerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FePerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FePerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FePerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FePerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FePerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FePerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FePerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &FePerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FePerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FePerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FePerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FePerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FePerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FePerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FePerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FePerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FePerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FePerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FePerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FePerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FePerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FePerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FePerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FePerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FePerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FePerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FePerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FePerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FePerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FePerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FePerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FePerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FePerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FePerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FePerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FePerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FePerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FePerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FePerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FePerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FePerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FePerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FePerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FePostPerovskite_slb_em, std::shared_ptr<FePostPerovskite_slb_em> >(m, "FePostPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FePostPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &FePostPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FePostPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FePostPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FePostPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FePostPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FePostPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FePostPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FePostPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &FePostPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FePostPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FePostPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FePostPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FePostPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FePostPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FePostPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FePostPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FePostPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FePostPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FePostPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FePostPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FePostPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FePostPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FePostPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FePostPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FePostPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FePostPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FePostPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FePostPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FePostPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FePostPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FePostPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FePostPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FePostPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FePostPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FePostPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FePostPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FePostPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FePostPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FePostPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FePostPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FePostPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FePostPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FePostPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FePostPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FeRingwoodite_slb_em, std::shared_ptr<FeRingwoodite_slb_em> >(m, "FeRingwoodite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FeRingwoodite_slb_em::name, "Name of endmember")
    .def("identifier", &FeRingwoodite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FeRingwoodite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeRingwoodite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeRingwoodite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeRingwoodite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeRingwoodite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FeRingwoodite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FeRingwoodite_slb_em::elements, "Vector of elements")
    
    .def("g", &FeRingwoodite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FeRingwoodite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FeRingwoodite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FeRingwoodite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FeRingwoodite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FeRingwoodite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FeRingwoodite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FeRingwoodite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FeRingwoodite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FeRingwoodite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FeRingwoodite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FeRingwoodite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FeRingwoodite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FeRingwoodite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FeRingwoodite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FeRingwoodite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FeRingwoodite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FeRingwoodite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FeRingwoodite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FeRingwoodite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FeRingwoodite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FeRingwoodite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeRingwoodite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FeRingwoodite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FeRingwoodite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FeRingwoodite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FeRingwoodite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FeRingwoodite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FeRingwoodite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FeRingwoodite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FeRingwoodite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FeRingwoodite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FeRingwoodite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FeRingwoodite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FeRingwoodite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FeRingwoodite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Ferrosilite_slb_em, std::shared_ptr<Ferrosilite_slb_em> >(m, "Ferrosilite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Ferrosilite_slb_em::name, "Name of endmember")
    .def("identifier", &Ferrosilite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Ferrosilite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Ferrosilite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Ferrosilite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Ferrosilite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Ferrosilite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Ferrosilite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Ferrosilite_slb_em::elements, "Vector of elements")
    
    .def("g", &Ferrosilite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Ferrosilite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Ferrosilite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Ferrosilite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Ferrosilite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Ferrosilite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Ferrosilite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Ferrosilite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Ferrosilite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Ferrosilite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Ferrosilite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Ferrosilite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Ferrosilite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Ferrosilite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Ferrosilite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Ferrosilite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Ferrosilite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Ferrosilite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Ferrosilite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Ferrosilite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Ferrosilite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Ferrosilite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Ferrosilite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Ferrosilite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Ferrosilite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Ferrosilite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Ferrosilite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Ferrosilite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Ferrosilite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Ferrosilite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Ferrosilite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Ferrosilite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Ferrosilite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Ferrosilite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Ferrosilite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Ferrosilite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<FeWadsleyite_slb_em, std::shared_ptr<FeWadsleyite_slb_em> >(m, "FeWadsleyite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &FeWadsleyite_slb_em::name, "Name of endmember")
    .def("identifier", &FeWadsleyite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &FeWadsleyite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeWadsleyite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeWadsleyite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeWadsleyite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeWadsleyite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &FeWadsleyite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &FeWadsleyite_slb_em::elements, "Vector of elements")
    
    .def("g", &FeWadsleyite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &FeWadsleyite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &FeWadsleyite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &FeWadsleyite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &FeWadsleyite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &FeWadsleyite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &FeWadsleyite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &FeWadsleyite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &FeWadsleyite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &FeWadsleyite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &FeWadsleyite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &FeWadsleyite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &FeWadsleyite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &FeWadsleyite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &FeWadsleyite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &FeWadsleyite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &FeWadsleyite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &FeWadsleyite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &FeWadsleyite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &FeWadsleyite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &FeWadsleyite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &FeWadsleyite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeWadsleyite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &FeWadsleyite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &FeWadsleyite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &FeWadsleyite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &FeWadsleyite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &FeWadsleyite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &FeWadsleyite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &FeWadsleyite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &FeWadsleyite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &FeWadsleyite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &FeWadsleyite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &FeWadsleyite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &FeWadsleyite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &FeWadsleyite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Forsterite_slb_em, std::shared_ptr<Forsterite_slb_em> >(m, "Forsterite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Forsterite_slb_em::name, "Name of endmember")
    .def("identifier", &Forsterite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Forsterite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Forsterite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Forsterite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Forsterite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Forsterite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Forsterite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Forsterite_slb_em::elements, "Vector of elements")
    
    .def("g", &Forsterite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Forsterite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Forsterite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Forsterite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Forsterite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Forsterite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Forsterite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Forsterite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Forsterite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Forsterite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Forsterite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Forsterite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Forsterite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Forsterite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Forsterite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Forsterite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Forsterite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Forsterite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Forsterite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Forsterite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Forsterite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Forsterite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Forsterite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Forsterite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Forsterite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Forsterite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Forsterite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Forsterite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Forsterite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Forsterite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Forsterite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Forsterite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Forsterite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Forsterite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Forsterite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Forsterite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Grossular_slb_em, std::shared_ptr<Grossular_slb_em> >(m, "Grossular_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Grossular_slb_em::name, "Name of endmember")
    .def("identifier", &Grossular_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Grossular_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Grossular_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Grossular_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Grossular_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Grossular_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Grossular_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Grossular_slb_em::elements, "Vector of elements")
    
    .def("g", &Grossular_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Grossular_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Grossular_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Grossular_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Grossular_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Grossular_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Grossular_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Grossular_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Grossular_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Grossular_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Grossular_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Grossular_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Grossular_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Grossular_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Grossular_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Grossular_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Grossular_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Grossular_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Grossular_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Grossular_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Grossular_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Grossular_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Grossular_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Grossular_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Grossular_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Grossular_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Grossular_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Grossular_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Grossular_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Grossular_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Grossular_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Grossular_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Grossular_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Grossular_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Grossular_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Grossular_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Hedenbergite_slb_em, std::shared_ptr<Hedenbergite_slb_em> >(m, "Hedenbergite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Hedenbergite_slb_em::name, "Name of endmember")
    .def("identifier", &Hedenbergite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Hedenbergite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Hedenbergite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Hedenbergite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Hedenbergite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Hedenbergite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Hedenbergite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Hedenbergite_slb_em::elements, "Vector of elements")
    
    .def("g", &Hedenbergite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Hedenbergite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Hedenbergite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Hedenbergite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Hedenbergite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Hedenbergite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Hedenbergite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Hedenbergite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Hedenbergite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Hedenbergite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Hedenbergite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Hedenbergite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Hedenbergite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Hedenbergite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Hedenbergite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Hedenbergite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Hedenbergite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Hedenbergite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Hedenbergite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Hedenbergite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Hedenbergite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Hedenbergite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Hedenbergite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Hedenbergite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Hedenbergite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Hedenbergite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Hedenbergite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Hedenbergite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Hedenbergite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Hedenbergite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Hedenbergite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Hedenbergite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Hedenbergite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Hedenbergite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Hedenbergite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Hedenbergite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Hercynite_slb_em, std::shared_ptr<Hercynite_slb_em> >(m, "Hercynite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Hercynite_slb_em::name, "Name of endmember")
    .def("identifier", &Hercynite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Hercynite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Hercynite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Hercynite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Hercynite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Hercynite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Hercynite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Hercynite_slb_em::elements, "Vector of elements")
    
    .def("g", &Hercynite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Hercynite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Hercynite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Hercynite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Hercynite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Hercynite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Hercynite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Hercynite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Hercynite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Hercynite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Hercynite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Hercynite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Hercynite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Hercynite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Hercynite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Hercynite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Hercynite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Hercynite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Hercynite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Hercynite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Hercynite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Hercynite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Hercynite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Hercynite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Hercynite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Hercynite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Hercynite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Hercynite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Hercynite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Hercynite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Hercynite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Hercynite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Hercynite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Hercynite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Hercynite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Hercynite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<HPClinoenstatite_slb_em, std::shared_ptr<HPClinoenstatite_slb_em> >(m, "HPClinoenstatite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &HPClinoenstatite_slb_em::name, "Name of endmember")
    .def("identifier", &HPClinoenstatite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &HPClinoenstatite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &HPClinoenstatite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &HPClinoenstatite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &HPClinoenstatite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &HPClinoenstatite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &HPClinoenstatite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &HPClinoenstatite_slb_em::elements, "Vector of elements")
    
    .def("g", &HPClinoenstatite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &HPClinoenstatite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &HPClinoenstatite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &HPClinoenstatite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &HPClinoenstatite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &HPClinoenstatite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &HPClinoenstatite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &HPClinoenstatite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &HPClinoenstatite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &HPClinoenstatite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &HPClinoenstatite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &HPClinoenstatite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &HPClinoenstatite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &HPClinoenstatite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &HPClinoenstatite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &HPClinoenstatite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &HPClinoenstatite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &HPClinoenstatite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &HPClinoenstatite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &HPClinoenstatite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &HPClinoenstatite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &HPClinoenstatite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&HPClinoenstatite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &HPClinoenstatite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &HPClinoenstatite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &HPClinoenstatite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &HPClinoenstatite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &HPClinoenstatite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &HPClinoenstatite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &HPClinoenstatite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &HPClinoenstatite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &HPClinoenstatite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &HPClinoenstatite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &HPClinoenstatite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &HPClinoenstatite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &HPClinoenstatite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<HPClinoferrosilite_slb_em, std::shared_ptr<HPClinoferrosilite_slb_em> >(m, "HPClinoferrosilite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &HPClinoferrosilite_slb_em::name, "Name of endmember")
    .def("identifier", &HPClinoferrosilite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &HPClinoferrosilite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &HPClinoferrosilite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &HPClinoferrosilite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &HPClinoferrosilite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &HPClinoferrosilite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &HPClinoferrosilite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &HPClinoferrosilite_slb_em::elements, "Vector of elements")
    
    .def("g", &HPClinoferrosilite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &HPClinoferrosilite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &HPClinoferrosilite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &HPClinoferrosilite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &HPClinoferrosilite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &HPClinoferrosilite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &HPClinoferrosilite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &HPClinoferrosilite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &HPClinoferrosilite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &HPClinoferrosilite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &HPClinoferrosilite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &HPClinoferrosilite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &HPClinoferrosilite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &HPClinoferrosilite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &HPClinoferrosilite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &HPClinoferrosilite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &HPClinoferrosilite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &HPClinoferrosilite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &HPClinoferrosilite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &HPClinoferrosilite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &HPClinoferrosilite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &HPClinoferrosilite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&HPClinoferrosilite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &HPClinoferrosilite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &HPClinoferrosilite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &HPClinoferrosilite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &HPClinoferrosilite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &HPClinoferrosilite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &HPClinoferrosilite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &HPClinoferrosilite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &HPClinoferrosilite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &HPClinoferrosilite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &HPClinoferrosilite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &HPClinoferrosilite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &HPClinoferrosilite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &HPClinoferrosilite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Jadeite_slb_em, std::shared_ptr<Jadeite_slb_em> >(m, "Jadeite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Jadeite_slb_em::name, "Name of endmember")
    .def("identifier", &Jadeite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Jadeite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Jadeite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Jadeite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Jadeite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Jadeite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Jadeite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Jadeite_slb_em::elements, "Vector of elements")
    
    .def("g", &Jadeite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Jadeite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Jadeite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Jadeite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Jadeite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Jadeite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Jadeite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Jadeite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Jadeite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Jadeite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Jadeite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Jadeite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Jadeite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Jadeite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Jadeite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Jadeite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Jadeite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Jadeite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Jadeite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Jadeite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Jadeite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Jadeite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Jadeite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Jadeite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Jadeite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Jadeite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Jadeite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Jadeite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Jadeite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Jadeite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Jadeite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Jadeite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Jadeite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Jadeite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Jadeite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Jadeite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Kyanite_slb_em, std::shared_ptr<Kyanite_slb_em> >(m, "Kyanite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Kyanite_slb_em::name, "Name of endmember")
    .def("identifier", &Kyanite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Kyanite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Kyanite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Kyanite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Kyanite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Kyanite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Kyanite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Kyanite_slb_em::elements, "Vector of elements")
    
    .def("g", &Kyanite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Kyanite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Kyanite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Kyanite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Kyanite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Kyanite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Kyanite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Kyanite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Kyanite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Kyanite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Kyanite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Kyanite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Kyanite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Kyanite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Kyanite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Kyanite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Kyanite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Kyanite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Kyanite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Kyanite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Kyanite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Kyanite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Kyanite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Kyanite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Kyanite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Kyanite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Kyanite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Kyanite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Kyanite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Kyanite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Kyanite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Kyanite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Kyanite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Kyanite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Kyanite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Kyanite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgAkimotoite_slb_em, std::shared_ptr<MgAkimotoite_slb_em> >(m, "MgAkimotoite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgAkimotoite_slb_em::name, "Name of endmember")
    .def("identifier", &MgAkimotoite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgAkimotoite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgAkimotoite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgAkimotoite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgAkimotoite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgAkimotoite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgAkimotoite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgAkimotoite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgAkimotoite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgAkimotoite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgAkimotoite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgAkimotoite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgAkimotoite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgAkimotoite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgAkimotoite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgAkimotoite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgAkimotoite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgAkimotoite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgAkimotoite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgAkimotoite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgAkimotoite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgAkimotoite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgAkimotoite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgAkimotoite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgAkimotoite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgAkimotoite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgAkimotoite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgAkimotoite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgAkimotoite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgAkimotoite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgAkimotoite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgAkimotoite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgAkimotoite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgAkimotoite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgAkimotoite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgAkimotoite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgAkimotoite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgAkimotoite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgAkimotoite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgAkimotoite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgAkimotoite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgAkimotoite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgAkimotoite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgAkimotoite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgCaFerrite_slb_em, std::shared_ptr<MgCaFerrite_slb_em> >(m, "MgCaFerrite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgCaFerrite_slb_em::name, "Name of endmember")
    .def("identifier", &MgCaFerrite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgCaFerrite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgCaFerrite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgCaFerrite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgCaFerrite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgCaFerrite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgCaFerrite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgCaFerrite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgCaFerrite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgCaFerrite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgCaFerrite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgCaFerrite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgCaFerrite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgCaFerrite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgCaFerrite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgCaFerrite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgCaFerrite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgCaFerrite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgCaFerrite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgCaFerrite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgCaFerrite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgCaFerrite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgCaFerrite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgCaFerrite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgCaFerrite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgCaFerrite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgCaFerrite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgCaFerrite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgCaFerrite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgCaFerrite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgCaFerrite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgCaFerrite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgCaFerrite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgCaFerrite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgCaFerrite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgCaFerrite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgCaFerrite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgCaFerrite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgCaFerrite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgCaFerrite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgCaFerrite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgCaFerrite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgCaFerrite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgCaFerrite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgMajorite_slb_em, std::shared_ptr<MgMajorite_slb_em> >(m, "MgMajorite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgMajorite_slb_em::name, "Name of endmember")
    .def("identifier", &MgMajorite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgMajorite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgMajorite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgMajorite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgMajorite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgMajorite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgMajorite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgMajorite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgMajorite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgMajorite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgMajorite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgMajorite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgMajorite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgMajorite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgMajorite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgMajorite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgMajorite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgMajorite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgMajorite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgMajorite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgMajorite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgMajorite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgMajorite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgMajorite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgMajorite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgMajorite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgMajorite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgMajorite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgMajorite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgMajorite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgMajorite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgMajorite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgMajorite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgMajorite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgMajorite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgMajorite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgMajorite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgMajorite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgMajorite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgMajorite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgMajorite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgMajorite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgMajorite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgMajorite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgPerovskite_slb_em, std::shared_ptr<MgPerovskite_slb_em> >(m, "MgPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &MgPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgPostPerovskite_slb_em, std::shared_ptr<MgPostPerovskite_slb_em> >(m, "MgPostPerovskite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgPostPerovskite_slb_em::name, "Name of endmember")
    .def("identifier", &MgPostPerovskite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgPostPerovskite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgPostPerovskite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgPostPerovskite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgPostPerovskite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgPostPerovskite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgPostPerovskite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgPostPerovskite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgPostPerovskite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgPostPerovskite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgPostPerovskite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgPostPerovskite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgPostPerovskite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgPostPerovskite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgPostPerovskite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgPostPerovskite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgPostPerovskite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgPostPerovskite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgPostPerovskite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgPostPerovskite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgPostPerovskite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgPostPerovskite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgPostPerovskite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgPostPerovskite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgPostPerovskite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgPostPerovskite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgPostPerovskite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgPostPerovskite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgPostPerovskite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgPostPerovskite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgPostPerovskite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgPostPerovskite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgPostPerovskite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgPostPerovskite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgPostPerovskite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgPostPerovskite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgPostPerovskite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgPostPerovskite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgPostPerovskite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgPostPerovskite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgPostPerovskite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgPostPerovskite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgPostPerovskite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgPostPerovskite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgRingwoodite_slb_em, std::shared_ptr<MgRingwoodite_slb_em> >(m, "MgRingwoodite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgRingwoodite_slb_em::name, "Name of endmember")
    .def("identifier", &MgRingwoodite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgRingwoodite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgRingwoodite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgRingwoodite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgRingwoodite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgRingwoodite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgRingwoodite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgRingwoodite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgRingwoodite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgRingwoodite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgRingwoodite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgRingwoodite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgRingwoodite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgRingwoodite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgRingwoodite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgRingwoodite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgRingwoodite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgRingwoodite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgRingwoodite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgRingwoodite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgRingwoodite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgRingwoodite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgRingwoodite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgRingwoodite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgRingwoodite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgRingwoodite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgRingwoodite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgRingwoodite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgRingwoodite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgRingwoodite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgRingwoodite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgRingwoodite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgRingwoodite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgRingwoodite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgRingwoodite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgRingwoodite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgRingwoodite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgRingwoodite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgRingwoodite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgRingwoodite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgRingwoodite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgRingwoodite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgRingwoodite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgRingwoodite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgSpinel_slb_em, std::shared_ptr<MgSpinel_slb_em> >(m, "MgSpinel_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgSpinel_slb_em::name, "Name of endmember")
    .def("identifier", &MgSpinel_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgSpinel_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgSpinel_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgSpinel_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgSpinel_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgSpinel_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgSpinel_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgSpinel_slb_em::elements, "Vector of elements")
    
    .def("g", &MgSpinel_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgSpinel_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgSpinel_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgSpinel_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgSpinel_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgSpinel_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgSpinel_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgSpinel_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgSpinel_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgSpinel_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgSpinel_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgSpinel_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgSpinel_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgSpinel_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgSpinel_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgSpinel_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgSpinel_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgSpinel_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgSpinel_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgSpinel_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgSpinel_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgSpinel_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgSpinel_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgSpinel_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgSpinel_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgSpinel_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgSpinel_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgSpinel_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgSpinel_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgSpinel_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgSpinel_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgSpinel_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgSpinel_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgSpinel_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgSpinel_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgSpinel_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgTschermaks_slb_em, std::shared_ptr<MgTschermaks_slb_em> >(m, "MgTschermaks_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgTschermaks_slb_em::name, "Name of endmember")
    .def("identifier", &MgTschermaks_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgTschermaks_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgTschermaks_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgTschermaks_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgTschermaks_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgTschermaks_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgTschermaks_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgTschermaks_slb_em::elements, "Vector of elements")
    
    .def("g", &MgTschermaks_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgTschermaks_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgTschermaks_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgTschermaks_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgTschermaks_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgTschermaks_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgTschermaks_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgTschermaks_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgTschermaks_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgTschermaks_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgTschermaks_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgTschermaks_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgTschermaks_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgTschermaks_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgTschermaks_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgTschermaks_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgTschermaks_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgTschermaks_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgTschermaks_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgTschermaks_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgTschermaks_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgTschermaks_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgTschermaks_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgTschermaks_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgTschermaks_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgTschermaks_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgTschermaks_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgTschermaks_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgTschermaks_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgTschermaks_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgTschermaks_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgTschermaks_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgTschermaks_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgTschermaks_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgTschermaks_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgTschermaks_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<MgWadsleyite_slb_em, std::shared_ptr<MgWadsleyite_slb_em> >(m, "MgWadsleyite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &MgWadsleyite_slb_em::name, "Name of endmember")
    .def("identifier", &MgWadsleyite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &MgWadsleyite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgWadsleyite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgWadsleyite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgWadsleyite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgWadsleyite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &MgWadsleyite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &MgWadsleyite_slb_em::elements, "Vector of elements")
    
    .def("g", &MgWadsleyite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &MgWadsleyite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &MgWadsleyite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &MgWadsleyite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &MgWadsleyite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &MgWadsleyite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &MgWadsleyite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &MgWadsleyite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &MgWadsleyite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &MgWadsleyite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &MgWadsleyite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &MgWadsleyite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &MgWadsleyite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &MgWadsleyite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &MgWadsleyite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &MgWadsleyite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &MgWadsleyite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &MgWadsleyite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &MgWadsleyite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &MgWadsleyite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &MgWadsleyite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &MgWadsleyite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgWadsleyite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &MgWadsleyite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &MgWadsleyite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &MgWadsleyite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &MgWadsleyite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &MgWadsleyite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &MgWadsleyite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &MgWadsleyite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &MgWadsleyite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &MgWadsleyite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &MgWadsleyite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &MgWadsleyite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &MgWadsleyite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &MgWadsleyite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<NaCaFerrite_slb_em, std::shared_ptr<NaCaFerrite_slb_em> >(m, "NaCaFerrite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &NaCaFerrite_slb_em::name, "Name of endmember")
    .def("identifier", &NaCaFerrite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &NaCaFerrite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &NaCaFerrite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &NaCaFerrite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &NaCaFerrite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &NaCaFerrite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &NaCaFerrite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &NaCaFerrite_slb_em::elements, "Vector of elements")
    
    .def("g", &NaCaFerrite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &NaCaFerrite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &NaCaFerrite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &NaCaFerrite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &NaCaFerrite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &NaCaFerrite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &NaCaFerrite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &NaCaFerrite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &NaCaFerrite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &NaCaFerrite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &NaCaFerrite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &NaCaFerrite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &NaCaFerrite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &NaCaFerrite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &NaCaFerrite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &NaCaFerrite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &NaCaFerrite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &NaCaFerrite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &NaCaFerrite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &NaCaFerrite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &NaCaFerrite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &NaCaFerrite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&NaCaFerrite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &NaCaFerrite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &NaCaFerrite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &NaCaFerrite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &NaCaFerrite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &NaCaFerrite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &NaCaFerrite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &NaCaFerrite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &NaCaFerrite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &NaCaFerrite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &NaCaFerrite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &NaCaFerrite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &NaCaFerrite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &NaCaFerrite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<NaMajorite_slb_em, std::shared_ptr<NaMajorite_slb_em> >(m, "NaMajorite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &NaMajorite_slb_em::name, "Name of endmember")
    .def("identifier", &NaMajorite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &NaMajorite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &NaMajorite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &NaMajorite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &NaMajorite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &NaMajorite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &NaMajorite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &NaMajorite_slb_em::elements, "Vector of elements")
    
    .def("g", &NaMajorite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &NaMajorite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &NaMajorite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &NaMajorite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &NaMajorite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &NaMajorite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &NaMajorite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &NaMajorite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &NaMajorite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &NaMajorite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &NaMajorite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &NaMajorite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &NaMajorite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &NaMajorite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &NaMajorite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &NaMajorite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &NaMajorite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &NaMajorite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &NaMajorite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &NaMajorite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &NaMajorite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &NaMajorite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&NaMajorite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &NaMajorite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &NaMajorite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &NaMajorite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &NaMajorite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &NaMajorite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &NaMajorite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &NaMajorite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &NaMajorite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &NaMajorite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &NaMajorite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &NaMajorite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &NaMajorite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &NaMajorite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Nepheline_slb_em, std::shared_ptr<Nepheline_slb_em> >(m, "Nepheline_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Nepheline_slb_em::name, "Name of endmember")
    .def("identifier", &Nepheline_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Nepheline_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Nepheline_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Nepheline_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Nepheline_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Nepheline_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Nepheline_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Nepheline_slb_em::elements, "Vector of elements")
    
    .def("g", &Nepheline_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Nepheline_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Nepheline_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Nepheline_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Nepheline_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Nepheline_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Nepheline_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Nepheline_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Nepheline_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Nepheline_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Nepheline_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Nepheline_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Nepheline_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Nepheline_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Nepheline_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Nepheline_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Nepheline_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Nepheline_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Nepheline_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Nepheline_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Nepheline_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Nepheline_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Nepheline_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Nepheline_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Nepheline_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Nepheline_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Nepheline_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Nepheline_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Nepheline_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Nepheline_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Nepheline_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Nepheline_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Nepheline_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Nepheline_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Nepheline_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Nepheline_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<OrthoDiopside_slb_em, std::shared_ptr<OrthoDiopside_slb_em> >(m, "OrthoDiopside_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &OrthoDiopside_slb_em::name, "Name of endmember")
    .def("identifier", &OrthoDiopside_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &OrthoDiopside_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &OrthoDiopside_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &OrthoDiopside_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &OrthoDiopside_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &OrthoDiopside_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &OrthoDiopside_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &OrthoDiopside_slb_em::elements, "Vector of elements")
    
    .def("g", &OrthoDiopside_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &OrthoDiopside_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &OrthoDiopside_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &OrthoDiopside_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &OrthoDiopside_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &OrthoDiopside_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &OrthoDiopside_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &OrthoDiopside_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &OrthoDiopside_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &OrthoDiopside_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &OrthoDiopside_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &OrthoDiopside_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &OrthoDiopside_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &OrthoDiopside_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &OrthoDiopside_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &OrthoDiopside_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &OrthoDiopside_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &OrthoDiopside_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &OrthoDiopside_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &OrthoDiopside_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &OrthoDiopside_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &OrthoDiopside_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&OrthoDiopside_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &OrthoDiopside_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &OrthoDiopside_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &OrthoDiopside_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &OrthoDiopside_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &OrthoDiopside_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &OrthoDiopside_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &OrthoDiopside_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &OrthoDiopside_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &OrthoDiopside_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &OrthoDiopside_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &OrthoDiopside_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &OrthoDiopside_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &OrthoDiopside_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Periclase_slb_em, std::shared_ptr<Periclase_slb_em> >(m, "Periclase_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Periclase_slb_em::name, "Name of endmember")
    .def("identifier", &Periclase_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Periclase_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Periclase_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Periclase_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Periclase_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Periclase_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Periclase_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Periclase_slb_em::elements, "Vector of elements")
    
    .def("g", &Periclase_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Periclase_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Periclase_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Periclase_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Periclase_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Periclase_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Periclase_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Periclase_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Periclase_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Periclase_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Periclase_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Periclase_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Periclase_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Periclase_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Periclase_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Periclase_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Periclase_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Periclase_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Periclase_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Periclase_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Periclase_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Periclase_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Periclase_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Periclase_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Periclase_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Periclase_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Periclase_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Periclase_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Periclase_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Periclase_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Periclase_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Periclase_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Periclase_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Periclase_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Periclase_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Periclase_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Pyrope_slb_em, std::shared_ptr<Pyrope_slb_em> >(m, "Pyrope_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Pyrope_slb_em::name, "Name of endmember")
    .def("identifier", &Pyrope_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Pyrope_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Pyrope_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Pyrope_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Pyrope_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Pyrope_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Pyrope_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Pyrope_slb_em::elements, "Vector of elements")
    
    .def("g", &Pyrope_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Pyrope_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Pyrope_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Pyrope_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Pyrope_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Pyrope_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Pyrope_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Pyrope_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Pyrope_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Pyrope_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Pyrope_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Pyrope_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Pyrope_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Pyrope_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Pyrope_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Pyrope_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Pyrope_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Pyrope_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Pyrope_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Pyrope_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Pyrope_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Pyrope_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Pyrope_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Pyrope_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Pyrope_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Pyrope_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Pyrope_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Pyrope_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Pyrope_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Pyrope_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Pyrope_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Pyrope_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Pyrope_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Pyrope_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Pyrope_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Pyrope_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Quartz_slb_em, std::shared_ptr<Quartz_slb_em> >(m, "Quartz_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Quartz_slb_em::name, "Name of endmember")
    .def("identifier", &Quartz_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Quartz_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Quartz_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Quartz_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Quartz_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Quartz_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Quartz_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Quartz_slb_em::elements, "Vector of elements")
    
    .def("g", &Quartz_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Quartz_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Quartz_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Quartz_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Quartz_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Quartz_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Quartz_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Quartz_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Quartz_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Quartz_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Quartz_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Quartz_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Quartz_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Quartz_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Quartz_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Quartz_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Quartz_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Quartz_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Quartz_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Quartz_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Quartz_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Quartz_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Quartz_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Quartz_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Quartz_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Quartz_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Quartz_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Quartz_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Quartz_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Quartz_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Quartz_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Quartz_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Quartz_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Quartz_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Quartz_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Quartz_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Seifertite_slb_em, std::shared_ptr<Seifertite_slb_em> >(m, "Seifertite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Seifertite_slb_em::name, "Name of endmember")
    .def("identifier", &Seifertite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Seifertite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Seifertite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Seifertite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Seifertite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Seifertite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Seifertite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Seifertite_slb_em::elements, "Vector of elements")
    
    .def("g", &Seifertite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Seifertite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Seifertite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Seifertite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Seifertite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Seifertite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Seifertite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Seifertite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Seifertite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Seifertite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Seifertite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Seifertite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Seifertite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Seifertite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Seifertite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Seifertite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Seifertite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Seifertite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Seifertite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Seifertite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Seifertite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Seifertite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Seifertite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Seifertite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Seifertite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Seifertite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Seifertite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Seifertite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Seifertite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Seifertite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Seifertite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Seifertite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Seifertite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Seifertite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Seifertite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Seifertite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Stishovite_slb_em, std::shared_ptr<Stishovite_slb_em> >(m, "Stishovite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Stishovite_slb_em::name, "Name of endmember")
    .def("identifier", &Stishovite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Stishovite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Stishovite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Stishovite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Stishovite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Stishovite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Stishovite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Stishovite_slb_em::elements, "Vector of elements")
    
    .def("g", &Stishovite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Stishovite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Stishovite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Stishovite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Stishovite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Stishovite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Stishovite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Stishovite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Stishovite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Stishovite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Stishovite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Stishovite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Stishovite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Stishovite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Stishovite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Stishovite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Stishovite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Stishovite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Stishovite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Stishovite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Stishovite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Stishovite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Stishovite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Stishovite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Stishovite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Stishovite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Stishovite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Stishovite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Stishovite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Stishovite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Stishovite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Stishovite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Stishovite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Stishovite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Stishovite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Stishovite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  py::class_<Wuestite_slb_em, std::shared_ptr<Wuestite_slb_em> >(m, "Wuestite_slb_em", py::module_local())
    .def(py::init<>())
    .def("name", &Wuestite_slb_em::name, "Name of endmember")
    .def("identifier", &Wuestite_slb_em::identifier, "Identifier")
    .def("tcg_build_version", &Wuestite_slb_em::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Wuestite_slb_em::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Wuestite_slb_em::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Wuestite_slb_em::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Wuestite_slb_em::formula, "Chemical formula of endmember")
    .def("mw", &Wuestite_slb_em::molecular_weight, "Molecular weight (g/mol)")
    .def("elements", &Wuestite_slb_em::elements, "Vector of elements")
    
    .def("g", &Wuestite_slb_em::G, R"pbdoc(
        Return the Gibbs free energy :math:`G`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        G: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdt", &Wuestite_slb_em::dGdT, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial T`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdT: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("dgdp", &Wuestite_slb_em::dGdP, R"pbdoc(
        Return the Gibbs free energy derivative :math:`\partial G / \partial P`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        dGdP: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdt2", &Wuestite_slb_em::d2GdT2, R"pbdoc(
        Return the Gibbs free energy second derivative :math:`\partial^2 G / \partial T^2`

        Parameters
        ----------
        T: float
            Temperature (in Kelvin)
        P: float
            Pressure (in bars)

        Returns
        -------
        d2GdT2: float
    )pbdoc", py::arg("T"), py::arg("P"))
    
    .def("d2gdtdp", &Wuestite_slb_em::d2GdTdP, "Gibbs free energy derivative d2gdtdp", py::arg("T"), py::arg("P"))
    .def("d2gdp2", &Wuestite_slb_em::d2GdP2, "Gibbs free energy derivative d2gdp2", py::arg("T"), py::arg("P"))
    .def("d3gdt3", &Wuestite_slb_em::d3GdT3, "Gibbs free energy derivative d3gdt3", py::arg("T"), py::arg("P"))
    .def("d3gdt2dp", &Wuestite_slb_em::d3GdT2dP,"Gibbs free energy derivative d3gdt2dp", py::arg("T"), py::arg("P"))
    .def("d3gdtdp2", &Wuestite_slb_em::d3GdTdP2,"Gibbs free energy derivative d3gdtdp2", py::arg("T"), py::arg("P"))
    .def("d3gdp3", &Wuestite_slb_em::d3GdP3, "Gibbs free energy derivative d3gdp3", py::arg("T"), py::arg("P"))
    .def("s", &Wuestite_slb_em::S, "Endmember entropy (J/K)", py::arg("T"), py::arg("P"))
    .def("v", &Wuestite_slb_em::V, "Endmember volume", py::arg("T"), py::arg("P"))
    .def("cv", &Wuestite_slb_em::Cv, "Heat capacity at constant volume", py::arg("T"), py::arg("P"))
    .def("cp", &Wuestite_slb_em::Cp, "Heat capacity a constant pressure", py::arg("T"), py::arg("P"))
    .def("dcpdt", &Wuestite_slb_em::dCpdT, py::arg("T"), py::arg("P"))
    .def("alpha", &Wuestite_slb_em::alpha, "Coefficient of thermal expansion", py::arg("T"), py::arg("P"))
    .def("beta", &Wuestite_slb_em::beta, "Bulk compressibility", py::arg("T"), py::arg("P"))
    .def("K", &Wuestite_slb_em::K, "Bulk modulus", py::arg("T"), py::arg("P"))
    .def("Kp", &Wuestite_slb_em::Kp,"Derivative of K with respect to pressure", py::arg("T"), py::arg("P"))
    .def("get_param_number", &Wuestite_slb_em::get_param_number, "Number of active parameters")
    .def("get_param_names", &Wuestite_slb_em::get_param_names, "Active parameter names")
    .def("get_param_units", &Wuestite_slb_em::get_param_units, "Active parameter units")
    .def("get_param_values", py::overload_cast<>(&Wuestite_slb_em::get_param_values), "Get active parameter values")
    .def("set_param_values", &Wuestite_slb_em::set_param_values, "Set active parameter values", py::arg("values"))
    .def("get_param_value", &Wuestite_slb_em::get_param_value, "Return value for a particular active parameter", py::arg("index"))
    .def("set_param_value", &Wuestite_slb_em::set_param_value, "Set value for a particular active parameter", py::arg("index"), py::arg("value"))
    .def("dparam_g", &Wuestite_slb_em::dparam_g, "dparam_g", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdt", &Wuestite_slb_em::dparam_dgdt, "dparam_dgdt", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_dgdp", &Wuestite_slb_em::dparam_dgdp, "dparam_dgdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdt2", &Wuestite_slb_em::dparam_d2gdt2, "dparam_d2gdt2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdtdp", &Wuestite_slb_em::dparam_d2gdtdp, "dparam_d2gdtdp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d2gdp2", &Wuestite_slb_em::dparam_d2gdp2, "dparam_d2gdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt3", &Wuestite_slb_em::dparam_d3gdt3, "dparam_d3gdt3", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdt2dp", &Wuestite_slb_em::dparam_d3gdt2dp, "dparam_d3gdt2dp", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdtdp2", &Wuestite_slb_em::dparam_d3gdtdp2, "dparam_d3gdtdp2", py::arg("T"), py::arg("P"), py::arg("index"))
    .def("dparam_d3gdp3", &Wuestite_slb_em::dparam_d3gdp3, "dparam_d3gdp3", py::arg("T"), py::arg("P"), py::arg("index"));

  

  // Phases
  py::class_<Akimotoite_slb_ph,  std::shared_ptr<Akimotoite_slb_ph> >(m, "Akimotoite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Akimotoite_slb_ph::identifier,"identifier")
    .def("name", &Akimotoite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Akimotoite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Akimotoite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Akimotoite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Akimotoite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Akimotoite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Akimotoite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Akimotoite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Akimotoite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Akimotoite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Akimotoite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Akimotoite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Akimotoite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Akimotoite_slb_ph::endmember_number)
    .def("endmember_name", &Akimotoite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Akimotoite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Akimotoite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Akimotoite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Akimotoite_slb_ph::species_number)
    .def("species_name", &Akimotoite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Akimotoite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Akimotoite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Akimotoite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Akimotoite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Akimotoite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Akimotoite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Akimotoite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Akimotoite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Akimotoite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Akimotoite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Akimotoite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Akimotoite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Akimotoite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Akimotoite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Akimotoite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Akimotoite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Akimotoite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Akimotoite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Akimotoite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Akimotoite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Akimotoite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Akimotoite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Akimotoite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Akimotoite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Akimotoite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Akimotoite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Akimotoite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Akimotoite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Akimotoite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Akimotoite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Akimotoite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Akimotoite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Akimotoite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Akimotoite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Akimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Akimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Akimotoite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Akimotoite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Akimotoite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Akimotoite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Akimotoite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Akimotoite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Akimotoite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Akimotoite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Akimotoite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Akimotoite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Akimotoite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Akimotoite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Akimotoite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Akimotoite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Akimotoite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Akimotoite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Akimotoite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Akimotoite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Akimotoite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Akimotoite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Akimotoite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Akimotoite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Akimotoite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Akimotoite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Akimotoite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Akimotoite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Akimotoite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Akimotoite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Akimotoite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Akimotoite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<AlAkimotoite_slb_ph,  std::shared_ptr<AlAkimotoite_slb_ph> >(m, "AlAkimotoite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &AlAkimotoite_slb_ph::identifier,"identifier")
    .def("name", &AlAkimotoite_slb_ph::name,"phase name")
    .def("tcg_build_version", &AlAkimotoite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlAkimotoite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlAkimotoite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlAkimotoite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlAkimotoite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &AlAkimotoite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &AlAkimotoite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &AlAkimotoite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &AlAkimotoite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &AlAkimotoite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &AlAkimotoite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &AlAkimotoite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &AlAkimotoite_slb_ph::endmember_number)
    .def("endmember_name", &AlAkimotoite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &AlAkimotoite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &AlAkimotoite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &AlAkimotoite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &AlAkimotoite_slb_ph::species_number)
    .def("species_name", &AlAkimotoite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &AlAkimotoite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &AlAkimotoite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &AlAkimotoite_slb_ph::species_elements, py::arg("i"))
    .def("g", &AlAkimotoite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &AlAkimotoite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &AlAkimotoite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &AlAkimotoite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &AlAkimotoite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &AlAkimotoite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &AlAkimotoite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &AlAkimotoite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &AlAkimotoite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &AlAkimotoite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &AlAkimotoite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &AlAkimotoite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &AlAkimotoite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &AlAkimotoite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &AlAkimotoite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &AlAkimotoite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &AlAkimotoite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &AlAkimotoite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &AlAkimotoite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &AlAkimotoite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &AlAkimotoite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &AlAkimotoite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &AlAkimotoite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &AlAkimotoite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &AlAkimotoite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &AlAkimotoite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &AlAkimotoite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &AlAkimotoite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &AlAkimotoite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &AlAkimotoite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &AlAkimotoite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &AlAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &AlAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlAkimotoite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &AlAkimotoite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &AlAkimotoite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &AlAkimotoite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &AlAkimotoite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &AlAkimotoite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &AlAkimotoite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &AlAkimotoite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &AlAkimotoite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &AlAkimotoite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &AlAkimotoite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &AlAkimotoite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &AlAkimotoite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &AlAkimotoite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &AlAkimotoite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &AlAkimotoite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &AlAkimotoite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &AlAkimotoite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &AlAkimotoite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &AlAkimotoite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &AlAkimotoite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &AlAkimotoite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &AlAkimotoite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &AlAkimotoite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &AlAkimotoite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &AlAkimotoite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &AlAkimotoite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &AlAkimotoite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&AlAkimotoite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&AlAkimotoite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Albite_slb_ph,  std::shared_ptr<Albite_slb_ph> >(m, "Albite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Albite_slb_ph::identifier,"identifier")
    .def("name", &Albite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Albite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Albite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Albite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Albite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Albite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Albite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Albite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Albite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Albite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Albite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Albite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Albite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Albite_slb_ph::endmember_number)
    .def("endmember_name", &Albite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Albite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Albite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Albite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Albite_slb_ph::species_number)
    .def("species_name", &Albite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Albite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Albite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Albite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Albite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Albite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Albite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Albite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Albite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Albite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Albite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Albite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Albite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Albite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Albite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Albite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Albite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Albite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Albite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Albite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Albite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Albite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Albite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Albite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Albite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Albite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Albite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Albite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Albite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Albite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Albite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Albite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Albite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Albite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Albite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Albite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Albite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Albite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Albite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Albite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Albite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Albite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Albite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Albite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Albite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Albite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Albite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Albite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Albite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Albite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Albite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Albite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Albite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Albite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Albite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Albite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Albite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Albite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Albite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Albite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Albite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Albite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Albite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Albite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Albite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Albite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Albite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Almandine_slb_ph,  std::shared_ptr<Almandine_slb_ph> >(m, "Almandine_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Almandine_slb_ph::identifier,"identifier")
    .def("name", &Almandine_slb_ph::name,"phase name")
    .def("tcg_build_version", &Almandine_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Almandine_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Almandine_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Almandine_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Almandine_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Almandine_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Almandine_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Almandine_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Almandine_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Almandine_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Almandine_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Almandine_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Almandine_slb_ph::endmember_number)
    .def("endmember_name", &Almandine_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Almandine_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Almandine_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Almandine_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Almandine_slb_ph::species_number)
    .def("species_name", &Almandine_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Almandine_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Almandine_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Almandine_slb_ph::species_elements, py::arg("i"))
    .def("g", &Almandine_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Almandine_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Almandine_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Almandine_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Almandine_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Almandine_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Almandine_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Almandine_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Almandine_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Almandine_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Almandine_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Almandine_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Almandine_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Almandine_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Almandine_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Almandine_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Almandine_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Almandine_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Almandine_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Almandine_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Almandine_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Almandine_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Almandine_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Almandine_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Almandine_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Almandine_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Almandine_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Almandine_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Almandine_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Almandine_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Almandine_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Almandine_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Almandine_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Almandine_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Almandine_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Almandine_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Almandine_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Almandine_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Almandine_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Almandine_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Almandine_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Almandine_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Almandine_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Almandine_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Almandine_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Almandine_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Almandine_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Almandine_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Almandine_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Almandine_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Almandine_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Almandine_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Almandine_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Almandine_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Almandine_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Almandine_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Almandine_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Almandine_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Almandine_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Almandine_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Almandine_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Almandine_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Almandine_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<AlPerovskite_slb_ph,  std::shared_ptr<AlPerovskite_slb_ph> >(m, "AlPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &AlPerovskite_slb_ph::identifier,"identifier")
    .def("name", &AlPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &AlPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &AlPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &AlPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &AlPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &AlPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &AlPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &AlPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &AlPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &AlPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &AlPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &AlPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &AlPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &AlPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &AlPerovskite_slb_ph::species_number)
    .def("species_name", &AlPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &AlPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &AlPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &AlPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &AlPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &AlPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &AlPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &AlPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &AlPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &AlPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &AlPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &AlPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &AlPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &AlPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &AlPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &AlPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &AlPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &AlPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &AlPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &AlPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &AlPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &AlPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &AlPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &AlPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &AlPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &AlPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &AlPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &AlPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &AlPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &AlPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &AlPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &AlPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &AlPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &AlPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &AlPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &AlPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &AlPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &AlPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &AlPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &AlPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &AlPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &AlPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &AlPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &AlPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &AlPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &AlPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &AlPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &AlPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &AlPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &AlPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &AlPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &AlPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &AlPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &AlPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &AlPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &AlPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &AlPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &AlPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &AlPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &AlPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &AlPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &AlPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &AlPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &AlPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&AlPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&AlPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<AlPostPerovskite_slb_ph,  std::shared_ptr<AlPostPerovskite_slb_ph> >(m, "AlPostPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &AlPostPerovskite_slb_ph::identifier,"identifier")
    .def("name", &AlPostPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &AlPostPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &AlPostPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &AlPostPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &AlPostPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &AlPostPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &AlPostPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &AlPostPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &AlPostPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &AlPostPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &AlPostPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &AlPostPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &AlPostPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &AlPostPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &AlPostPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &AlPostPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &AlPostPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &AlPostPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &AlPostPerovskite_slb_ph::species_number)
    .def("species_name", &AlPostPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &AlPostPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &AlPostPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &AlPostPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &AlPostPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &AlPostPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &AlPostPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &AlPostPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &AlPostPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &AlPostPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &AlPostPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &AlPostPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &AlPostPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &AlPostPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &AlPostPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &AlPostPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &AlPostPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &AlPostPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &AlPostPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &AlPostPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &AlPostPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &AlPostPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &AlPostPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &AlPostPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &AlPostPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &AlPostPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &AlPostPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &AlPostPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &AlPostPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &AlPostPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &AlPostPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &AlPostPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &AlPostPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &AlPostPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &AlPostPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &AlPostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &AlPostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&AlPostPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &AlPostPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &AlPostPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &AlPostPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &AlPostPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &AlPostPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &AlPostPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &AlPostPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &AlPostPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &AlPostPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &AlPostPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &AlPostPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &AlPostPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &AlPostPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &AlPostPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &AlPostPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &AlPostPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &AlPostPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &AlPostPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &AlPostPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &AlPostPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &AlPostPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &AlPostPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &AlPostPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &AlPostPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &AlPostPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &AlPostPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &AlPostPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&AlPostPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&AlPostPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Anorthite_slb_ph,  std::shared_ptr<Anorthite_slb_ph> >(m, "Anorthite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Anorthite_slb_ph::identifier,"identifier")
    .def("name", &Anorthite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Anorthite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Anorthite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Anorthite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Anorthite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Anorthite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Anorthite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Anorthite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Anorthite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Anorthite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Anorthite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Anorthite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Anorthite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Anorthite_slb_ph::endmember_number)
    .def("endmember_name", &Anorthite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Anorthite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Anorthite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Anorthite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Anorthite_slb_ph::species_number)
    .def("species_name", &Anorthite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Anorthite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Anorthite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Anorthite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Anorthite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Anorthite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Anorthite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Anorthite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Anorthite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Anorthite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Anorthite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Anorthite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Anorthite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Anorthite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Anorthite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Anorthite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Anorthite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Anorthite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Anorthite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Anorthite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Anorthite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Anorthite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Anorthite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Anorthite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Anorthite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Anorthite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Anorthite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Anorthite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Anorthite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Anorthite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Anorthite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Anorthite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Anorthite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Anorthite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Anorthite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Anorthite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Anorthite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Anorthite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Anorthite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Anorthite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Anorthite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Anorthite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Anorthite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Anorthite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Anorthite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Anorthite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Anorthite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Anorthite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Anorthite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Anorthite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Anorthite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Anorthite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Anorthite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Anorthite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Anorthite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Anorthite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Anorthite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Anorthite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Anorthite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Anorthite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Anorthite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Anorthite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Anorthite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Anorthite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Anorthite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Anorthite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Anorthite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<CaFerritePhase_slb_ph,  std::shared_ptr<CaFerritePhase_slb_ph> >(m, "CaFerritePhase_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &CaFerritePhase_slb_ph::identifier,"identifier")
    .def("name", &CaFerritePhase_slb_ph::name,"phase name")
    .def("tcg_build_version", &CaFerritePhase_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &CaFerritePhase_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &CaFerritePhase_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &CaFerritePhase_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &CaFerritePhase_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &CaFerritePhase_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &CaFerritePhase_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &CaFerritePhase_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &CaFerritePhase_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &CaFerritePhase_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &CaFerritePhase_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &CaFerritePhase_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &CaFerritePhase_slb_ph::endmember_number)
    .def("endmember_name", &CaFerritePhase_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &CaFerritePhase_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &CaFerritePhase_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &CaFerritePhase_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &CaFerritePhase_slb_ph::species_number)
    .def("species_name", &CaFerritePhase_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &CaFerritePhase_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &CaFerritePhase_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &CaFerritePhase_slb_ph::species_elements, py::arg("i"))
    .def("g", &CaFerritePhase_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &CaFerritePhase_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &CaFerritePhase_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &CaFerritePhase_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &CaFerritePhase_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &CaFerritePhase_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &CaFerritePhase_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &CaFerritePhase_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &CaFerritePhase_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &CaFerritePhase_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &CaFerritePhase_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &CaFerritePhase_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &CaFerritePhase_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &CaFerritePhase_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &CaFerritePhase_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &CaFerritePhase_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &CaFerritePhase_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &CaFerritePhase_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &CaFerritePhase_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &CaFerritePhase_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &CaFerritePhase_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &CaFerritePhase_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &CaFerritePhase_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &CaFerritePhase_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &CaFerritePhase_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &CaFerritePhase_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &CaFerritePhase_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &CaFerritePhase_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &CaFerritePhase_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &CaFerritePhase_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &CaFerritePhase_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &CaFerritePhase_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &CaFerritePhase_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&CaFerritePhase_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &CaFerritePhase_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &CaFerritePhase_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &CaFerritePhase_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &CaFerritePhase_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &CaFerritePhase_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &CaFerritePhase_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &CaFerritePhase_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &CaFerritePhase_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &CaFerritePhase_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &CaFerritePhase_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &CaFerritePhase_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &CaFerritePhase_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &CaFerritePhase_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &CaFerritePhase_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &CaFerritePhase_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &CaFerritePhase_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &CaFerritePhase_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &CaFerritePhase_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &CaFerritePhase_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &CaFerritePhase_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &CaFerritePhase_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &CaFerritePhase_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &CaFerritePhase_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &CaFerritePhase_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &CaFerritePhase_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &CaFerritePhase_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &CaFerritePhase_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&CaFerritePhase_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&CaFerritePhase_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<CaPerovskite_slb_ph,  std::shared_ptr<CaPerovskite_slb_ph> >(m, "CaPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &CaPerovskite_slb_ph::identifier,"identifier")
    .def("name", &CaPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &CaPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &CaPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &CaPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &CaPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &CaPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &CaPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &CaPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &CaPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &CaPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &CaPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &CaPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &CaPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &CaPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &CaPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &CaPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &CaPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &CaPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &CaPerovskite_slb_ph::species_number)
    .def("species_name", &CaPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &CaPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &CaPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &CaPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &CaPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &CaPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &CaPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &CaPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &CaPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &CaPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &CaPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &CaPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &CaPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &CaPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &CaPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &CaPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &CaPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &CaPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &CaPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &CaPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &CaPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &CaPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &CaPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &CaPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &CaPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &CaPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &CaPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &CaPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &CaPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &CaPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &CaPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &CaPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &CaPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &CaPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &CaPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &CaPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &CaPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&CaPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &CaPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &CaPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &CaPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &CaPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &CaPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &CaPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &CaPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &CaPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &CaPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &CaPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &CaPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &CaPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &CaPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &CaPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &CaPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &CaPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &CaPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &CaPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &CaPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &CaPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &CaPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &CaPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &CaPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &CaPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &CaPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &CaPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &CaPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&CaPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&CaPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<CaTschermaks_slb_ph,  std::shared_ptr<CaTschermaks_slb_ph> >(m, "CaTschermaks_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &CaTschermaks_slb_ph::identifier,"identifier")
    .def("name", &CaTschermaks_slb_ph::name,"phase name")
    .def("tcg_build_version", &CaTschermaks_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &CaTschermaks_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &CaTschermaks_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &CaTschermaks_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &CaTschermaks_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &CaTschermaks_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &CaTschermaks_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &CaTschermaks_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &CaTschermaks_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &CaTschermaks_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &CaTschermaks_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &CaTschermaks_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &CaTschermaks_slb_ph::endmember_number)
    .def("endmember_name", &CaTschermaks_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &CaTschermaks_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &CaTschermaks_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &CaTschermaks_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &CaTschermaks_slb_ph::species_number)
    .def("species_name", &CaTschermaks_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &CaTschermaks_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &CaTschermaks_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &CaTschermaks_slb_ph::species_elements, py::arg("i"))
    .def("g", &CaTschermaks_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &CaTschermaks_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &CaTschermaks_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &CaTschermaks_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &CaTschermaks_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &CaTschermaks_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &CaTschermaks_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &CaTschermaks_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &CaTschermaks_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &CaTschermaks_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &CaTschermaks_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &CaTschermaks_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &CaTschermaks_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &CaTschermaks_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &CaTschermaks_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &CaTschermaks_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &CaTschermaks_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &CaTschermaks_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &CaTschermaks_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &CaTschermaks_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &CaTschermaks_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &CaTschermaks_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &CaTschermaks_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &CaTschermaks_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &CaTschermaks_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &CaTschermaks_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &CaTschermaks_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &CaTschermaks_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &CaTschermaks_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &CaTschermaks_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &CaTschermaks_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &CaTschermaks_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &CaTschermaks_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&CaTschermaks_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &CaTschermaks_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &CaTschermaks_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &CaTschermaks_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &CaTschermaks_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &CaTschermaks_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &CaTschermaks_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &CaTschermaks_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &CaTschermaks_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &CaTschermaks_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &CaTschermaks_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &CaTschermaks_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &CaTschermaks_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &CaTschermaks_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &CaTschermaks_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &CaTschermaks_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &CaTschermaks_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &CaTschermaks_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &CaTschermaks_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &CaTschermaks_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &CaTschermaks_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &CaTschermaks_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &CaTschermaks_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &CaTschermaks_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &CaTschermaks_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &CaTschermaks_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &CaTschermaks_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &CaTschermaks_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&CaTschermaks_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&CaTschermaks_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Clinoenstatite_slb_ph,  std::shared_ptr<Clinoenstatite_slb_ph> >(m, "Clinoenstatite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Clinoenstatite_slb_ph::identifier,"identifier")
    .def("name", &Clinoenstatite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Clinoenstatite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Clinoenstatite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Clinoenstatite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Clinoenstatite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Clinoenstatite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Clinoenstatite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Clinoenstatite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Clinoenstatite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Clinoenstatite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Clinoenstatite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Clinoenstatite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Clinoenstatite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Clinoenstatite_slb_ph::endmember_number)
    .def("endmember_name", &Clinoenstatite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Clinoenstatite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Clinoenstatite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Clinoenstatite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Clinoenstatite_slb_ph::species_number)
    .def("species_name", &Clinoenstatite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Clinoenstatite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Clinoenstatite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Clinoenstatite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Clinoenstatite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Clinoenstatite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Clinoenstatite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Clinoenstatite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Clinoenstatite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Clinoenstatite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Clinoenstatite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Clinoenstatite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Clinoenstatite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Clinoenstatite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Clinoenstatite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Clinoenstatite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Clinoenstatite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Clinoenstatite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Clinoenstatite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Clinoenstatite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Clinoenstatite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Clinoenstatite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Clinoenstatite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Clinoenstatite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Clinoenstatite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Clinoenstatite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Clinoenstatite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Clinoenstatite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Clinoenstatite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Clinoenstatite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Clinoenstatite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Clinoenstatite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Clinoenstatite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Clinoenstatite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Clinoenstatite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Clinoenstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Clinoenstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Clinoenstatite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Clinoenstatite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Clinoenstatite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Clinoenstatite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Clinoenstatite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Clinoenstatite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Clinoenstatite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Clinoenstatite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Clinoenstatite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Clinoenstatite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Clinoenstatite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Clinoenstatite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Clinoenstatite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Clinoenstatite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Clinoenstatite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Clinoenstatite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Clinoenstatite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Clinoenstatite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Clinoenstatite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Clinoenstatite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Clinoenstatite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Clinoenstatite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Clinoenstatite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Clinoenstatite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Clinoenstatite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Clinoenstatite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Clinoenstatite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Clinoenstatite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Clinoenstatite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Clinoenstatite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Clinopyroxene_slb_ph,  std::shared_ptr<Clinopyroxene_slb_ph> >(m, "Clinopyroxene_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Clinopyroxene_slb_ph::identifier,"identifier")
    .def("name", &Clinopyroxene_slb_ph::name,"phase name")
    .def("tcg_build_version", &Clinopyroxene_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Clinopyroxene_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Clinopyroxene_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Clinopyroxene_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Clinopyroxene_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Clinopyroxene_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Clinopyroxene_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Clinopyroxene_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Clinopyroxene_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Clinopyroxene_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Clinopyroxene_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Clinopyroxene_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Clinopyroxene_slb_ph::endmember_number)
    .def("endmember_name", &Clinopyroxene_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Clinopyroxene_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Clinopyroxene_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Clinopyroxene_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Clinopyroxene_slb_ph::species_number)
    .def("species_name", &Clinopyroxene_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Clinopyroxene_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Clinopyroxene_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Clinopyroxene_slb_ph::species_elements, py::arg("i"))
    .def("g", &Clinopyroxene_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Clinopyroxene_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Clinopyroxene_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Clinopyroxene_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Clinopyroxene_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Clinopyroxene_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Clinopyroxene_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Clinopyroxene_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Clinopyroxene_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Clinopyroxene_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Clinopyroxene_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Clinopyroxene_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Clinopyroxene_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Clinopyroxene_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Clinopyroxene_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Clinopyroxene_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Clinopyroxene_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Clinopyroxene_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Clinopyroxene_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Clinopyroxene_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Clinopyroxene_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Clinopyroxene_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Clinopyroxene_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Clinopyroxene_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Clinopyroxene_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Clinopyroxene_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Clinopyroxene_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Clinopyroxene_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Clinopyroxene_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Clinopyroxene_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Clinopyroxene_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Clinopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Clinopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Clinopyroxene_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Clinopyroxene_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Clinopyroxene_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Clinopyroxene_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Clinopyroxene_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Clinopyroxene_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Clinopyroxene_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Clinopyroxene_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Clinopyroxene_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Clinopyroxene_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Clinopyroxene_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Clinopyroxene_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Clinopyroxene_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Clinopyroxene_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Clinopyroxene_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Clinopyroxene_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Clinopyroxene_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Clinopyroxene_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Clinopyroxene_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Clinopyroxene_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Clinopyroxene_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Clinopyroxene_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Clinopyroxene_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Clinopyroxene_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Clinopyroxene_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Clinopyroxene_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Clinopyroxene_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Clinopyroxene_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Clinopyroxene_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Clinopyroxene_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Coesite_slb_ph,  std::shared_ptr<Coesite_slb_ph> >(m, "Coesite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Coesite_slb_ph::identifier,"identifier")
    .def("name", &Coesite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Coesite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Coesite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Coesite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Coesite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Coesite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Coesite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Coesite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Coesite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Coesite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Coesite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Coesite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Coesite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Coesite_slb_ph::endmember_number)
    .def("endmember_name", &Coesite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Coesite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Coesite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Coesite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Coesite_slb_ph::species_number)
    .def("species_name", &Coesite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Coesite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Coesite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Coesite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Coesite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Coesite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Coesite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Coesite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Coesite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Coesite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Coesite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Coesite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Coesite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Coesite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Coesite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Coesite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Coesite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Coesite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Coesite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Coesite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Coesite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Coesite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Coesite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Coesite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Coesite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Coesite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Coesite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Coesite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Coesite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Coesite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Coesite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Coesite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Coesite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Coesite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Coesite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Coesite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Coesite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Coesite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Coesite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Coesite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Coesite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Coesite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Coesite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Coesite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Coesite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Coesite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Coesite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Coesite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Coesite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Coesite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Coesite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Coesite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Coesite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Coesite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Coesite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Coesite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Coesite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Coesite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Coesite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Coesite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Coesite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Coesite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Coesite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Coesite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Coesite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Coesite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Coesite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Diopside_slb_ph,  std::shared_ptr<Diopside_slb_ph> >(m, "Diopside_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Diopside_slb_ph::identifier,"identifier")
    .def("name", &Diopside_slb_ph::name,"phase name")
    .def("tcg_build_version", &Diopside_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Diopside_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Diopside_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Diopside_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Diopside_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Diopside_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Diopside_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Diopside_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Diopside_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Diopside_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Diopside_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Diopside_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Diopside_slb_ph::endmember_number)
    .def("endmember_name", &Diopside_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Diopside_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Diopside_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Diopside_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Diopside_slb_ph::species_number)
    .def("species_name", &Diopside_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Diopside_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Diopside_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Diopside_slb_ph::species_elements, py::arg("i"))
    .def("g", &Diopside_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Diopside_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Diopside_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Diopside_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Diopside_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Diopside_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Diopside_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Diopside_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Diopside_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Diopside_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Diopside_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Diopside_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Diopside_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Diopside_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Diopside_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Diopside_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Diopside_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Diopside_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Diopside_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Diopside_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Diopside_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Diopside_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Diopside_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Diopside_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Diopside_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Diopside_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Diopside_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Diopside_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Diopside_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Diopside_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Diopside_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Diopside_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Diopside_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Diopside_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Diopside_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Diopside_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Diopside_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Diopside_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Diopside_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Diopside_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Diopside_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Diopside_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Diopside_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Diopside_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Diopside_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Diopside_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Diopside_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Diopside_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Diopside_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Diopside_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Diopside_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Diopside_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Diopside_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Diopside_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Diopside_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Diopside_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Diopside_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Diopside_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Diopside_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Diopside_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Diopside_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Diopside_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Diopside_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Enstatite_slb_ph,  std::shared_ptr<Enstatite_slb_ph> >(m, "Enstatite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Enstatite_slb_ph::identifier,"identifier")
    .def("name", &Enstatite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Enstatite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Enstatite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Enstatite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Enstatite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Enstatite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Enstatite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Enstatite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Enstatite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Enstatite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Enstatite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Enstatite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Enstatite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Enstatite_slb_ph::endmember_number)
    .def("endmember_name", &Enstatite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Enstatite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Enstatite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Enstatite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Enstatite_slb_ph::species_number)
    .def("species_name", &Enstatite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Enstatite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Enstatite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Enstatite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Enstatite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Enstatite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Enstatite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Enstatite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Enstatite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Enstatite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Enstatite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Enstatite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Enstatite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Enstatite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Enstatite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Enstatite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Enstatite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Enstatite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Enstatite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Enstatite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Enstatite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Enstatite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Enstatite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Enstatite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Enstatite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Enstatite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Enstatite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Enstatite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Enstatite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Enstatite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Enstatite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Enstatite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Enstatite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Enstatite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Enstatite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Enstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Enstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Enstatite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Enstatite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Enstatite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Enstatite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Enstatite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Enstatite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Enstatite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Enstatite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Enstatite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Enstatite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Enstatite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Enstatite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Enstatite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Enstatite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Enstatite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Enstatite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Enstatite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Enstatite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Enstatite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Enstatite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Enstatite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Enstatite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Enstatite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Enstatite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Enstatite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Enstatite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Enstatite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Enstatite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Enstatite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Enstatite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Fayalite_slb_ph,  std::shared_ptr<Fayalite_slb_ph> >(m, "Fayalite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Fayalite_slb_ph::identifier,"identifier")
    .def("name", &Fayalite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Fayalite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Fayalite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Fayalite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Fayalite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Fayalite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Fayalite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Fayalite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Fayalite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Fayalite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Fayalite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Fayalite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Fayalite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Fayalite_slb_ph::endmember_number)
    .def("endmember_name", &Fayalite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Fayalite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Fayalite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Fayalite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Fayalite_slb_ph::species_number)
    .def("species_name", &Fayalite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Fayalite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Fayalite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Fayalite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Fayalite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Fayalite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Fayalite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Fayalite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Fayalite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Fayalite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Fayalite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Fayalite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Fayalite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Fayalite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Fayalite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Fayalite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Fayalite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Fayalite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Fayalite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Fayalite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Fayalite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Fayalite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Fayalite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Fayalite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Fayalite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Fayalite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Fayalite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Fayalite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Fayalite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Fayalite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Fayalite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Fayalite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Fayalite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Fayalite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Fayalite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Fayalite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Fayalite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Fayalite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Fayalite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Fayalite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Fayalite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Fayalite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Fayalite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Fayalite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Fayalite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Fayalite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Fayalite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Fayalite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Fayalite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Fayalite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Fayalite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Fayalite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Fayalite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Fayalite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Fayalite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Fayalite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Fayalite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Fayalite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Fayalite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Fayalite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Fayalite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Fayalite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Fayalite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Fayalite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Fayalite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Fayalite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Fayalite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FeAkimotoite_slb_ph,  std::shared_ptr<FeAkimotoite_slb_ph> >(m, "FeAkimotoite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FeAkimotoite_slb_ph::identifier,"identifier")
    .def("name", &FeAkimotoite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FeAkimotoite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeAkimotoite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeAkimotoite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeAkimotoite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeAkimotoite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FeAkimotoite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FeAkimotoite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FeAkimotoite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FeAkimotoite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FeAkimotoite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FeAkimotoite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FeAkimotoite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FeAkimotoite_slb_ph::endmember_number)
    .def("endmember_name", &FeAkimotoite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FeAkimotoite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FeAkimotoite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FeAkimotoite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FeAkimotoite_slb_ph::species_number)
    .def("species_name", &FeAkimotoite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FeAkimotoite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FeAkimotoite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FeAkimotoite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FeAkimotoite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FeAkimotoite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FeAkimotoite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FeAkimotoite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FeAkimotoite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FeAkimotoite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FeAkimotoite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FeAkimotoite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FeAkimotoite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FeAkimotoite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FeAkimotoite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FeAkimotoite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FeAkimotoite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FeAkimotoite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FeAkimotoite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FeAkimotoite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FeAkimotoite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FeAkimotoite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FeAkimotoite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FeAkimotoite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FeAkimotoite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FeAkimotoite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FeAkimotoite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FeAkimotoite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FeAkimotoite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FeAkimotoite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FeAkimotoite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FeAkimotoite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FeAkimotoite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FeAkimotoite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FeAkimotoite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FeAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FeAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeAkimotoite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FeAkimotoite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FeAkimotoite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FeAkimotoite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FeAkimotoite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FeAkimotoite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FeAkimotoite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FeAkimotoite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FeAkimotoite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FeAkimotoite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FeAkimotoite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FeAkimotoite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FeAkimotoite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FeAkimotoite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FeAkimotoite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FeAkimotoite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FeAkimotoite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FeAkimotoite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FeAkimotoite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FeAkimotoite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FeAkimotoite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FeAkimotoite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FeAkimotoite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FeAkimotoite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FeAkimotoite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FeAkimotoite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FeAkimotoite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FeAkimotoite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FeAkimotoite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FeAkimotoite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FeCaFerrite_slb_ph,  std::shared_ptr<FeCaFerrite_slb_ph> >(m, "FeCaFerrite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FeCaFerrite_slb_ph::identifier,"identifier")
    .def("name", &FeCaFerrite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FeCaFerrite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeCaFerrite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeCaFerrite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeCaFerrite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeCaFerrite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FeCaFerrite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FeCaFerrite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FeCaFerrite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FeCaFerrite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FeCaFerrite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FeCaFerrite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FeCaFerrite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FeCaFerrite_slb_ph::endmember_number)
    .def("endmember_name", &FeCaFerrite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FeCaFerrite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FeCaFerrite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FeCaFerrite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FeCaFerrite_slb_ph::species_number)
    .def("species_name", &FeCaFerrite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FeCaFerrite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FeCaFerrite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FeCaFerrite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FeCaFerrite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FeCaFerrite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FeCaFerrite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FeCaFerrite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FeCaFerrite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FeCaFerrite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FeCaFerrite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FeCaFerrite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FeCaFerrite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FeCaFerrite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FeCaFerrite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FeCaFerrite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FeCaFerrite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FeCaFerrite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FeCaFerrite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FeCaFerrite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FeCaFerrite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FeCaFerrite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FeCaFerrite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FeCaFerrite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FeCaFerrite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FeCaFerrite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FeCaFerrite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FeCaFerrite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FeCaFerrite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FeCaFerrite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FeCaFerrite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FeCaFerrite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FeCaFerrite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FeCaFerrite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FeCaFerrite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FeCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FeCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeCaFerrite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FeCaFerrite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FeCaFerrite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FeCaFerrite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FeCaFerrite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FeCaFerrite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FeCaFerrite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FeCaFerrite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FeCaFerrite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FeCaFerrite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FeCaFerrite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FeCaFerrite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FeCaFerrite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FeCaFerrite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FeCaFerrite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FeCaFerrite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FeCaFerrite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FeCaFerrite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FeCaFerrite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FeCaFerrite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FeCaFerrite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FeCaFerrite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FeCaFerrite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FeCaFerrite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FeCaFerrite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FeCaFerrite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FeCaFerrite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FeCaFerrite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FeCaFerrite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FeCaFerrite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Feldspar_slb_ph,  std::shared_ptr<Feldspar_slb_ph> >(m, "Feldspar_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Feldspar_slb_ph::identifier,"identifier")
    .def("name", &Feldspar_slb_ph::name,"phase name")
    .def("tcg_build_version", &Feldspar_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Feldspar_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Feldspar_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Feldspar_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Feldspar_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Feldspar_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Feldspar_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Feldspar_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Feldspar_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Feldspar_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Feldspar_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Feldspar_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Feldspar_slb_ph::endmember_number)
    .def("endmember_name", &Feldspar_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Feldspar_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Feldspar_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Feldspar_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Feldspar_slb_ph::species_number)
    .def("species_name", &Feldspar_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Feldspar_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Feldspar_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Feldspar_slb_ph::species_elements, py::arg("i"))
    .def("g", &Feldspar_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Feldspar_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Feldspar_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Feldspar_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Feldspar_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Feldspar_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Feldspar_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Feldspar_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Feldspar_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Feldspar_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Feldspar_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Feldspar_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Feldspar_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Feldspar_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Feldspar_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Feldspar_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Feldspar_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Feldspar_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Feldspar_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Feldspar_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Feldspar_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Feldspar_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Feldspar_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Feldspar_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Feldspar_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Feldspar_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Feldspar_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Feldspar_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Feldspar_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Feldspar_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Feldspar_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Feldspar_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Feldspar_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Feldspar_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Feldspar_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Feldspar_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Feldspar_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Feldspar_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Feldspar_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Feldspar_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Feldspar_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Feldspar_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Feldspar_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Feldspar_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Feldspar_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Feldspar_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Feldspar_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Feldspar_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Feldspar_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Feldspar_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Feldspar_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Feldspar_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Feldspar_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Feldspar_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Feldspar_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Feldspar_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Feldspar_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Feldspar_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Feldspar_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Feldspar_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Feldspar_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Feldspar_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Feldspar_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FePerovskite_slb_ph,  std::shared_ptr<FePerovskite_slb_ph> >(m, "FePerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FePerovskite_slb_ph::identifier,"identifier")
    .def("name", &FePerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FePerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FePerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FePerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FePerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FePerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FePerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FePerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FePerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FePerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FePerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FePerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FePerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FePerovskite_slb_ph::endmember_number)
    .def("endmember_name", &FePerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FePerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FePerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FePerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FePerovskite_slb_ph::species_number)
    .def("species_name", &FePerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FePerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FePerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FePerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FePerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FePerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FePerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FePerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FePerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FePerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FePerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FePerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FePerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FePerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FePerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FePerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FePerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FePerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FePerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FePerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FePerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FePerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FePerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FePerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FePerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FePerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FePerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FePerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FePerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FePerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FePerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FePerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FePerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FePerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FePerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FePerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FePerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FePerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FePerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FePerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FePerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FePerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FePerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FePerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FePerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FePerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FePerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FePerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FePerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FePerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FePerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FePerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FePerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FePerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FePerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FePerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FePerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FePerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FePerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FePerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FePerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FePerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FePerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FePerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FePerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FePerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FePerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FePostPerovskite_slb_ph,  std::shared_ptr<FePostPerovskite_slb_ph> >(m, "FePostPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FePostPerovskite_slb_ph::identifier,"identifier")
    .def("name", &FePostPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FePostPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FePostPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FePostPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FePostPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FePostPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FePostPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FePostPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FePostPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FePostPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FePostPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FePostPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FePostPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FePostPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &FePostPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FePostPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FePostPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FePostPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FePostPerovskite_slb_ph::species_number)
    .def("species_name", &FePostPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FePostPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FePostPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FePostPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FePostPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FePostPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FePostPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FePostPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FePostPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FePostPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FePostPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FePostPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FePostPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FePostPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FePostPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FePostPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FePostPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FePostPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FePostPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FePostPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FePostPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FePostPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FePostPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FePostPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FePostPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FePostPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FePostPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FePostPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FePostPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FePostPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FePostPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FePostPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FePostPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FePostPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FePostPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FePostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FePostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FePostPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FePostPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FePostPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FePostPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FePostPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FePostPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FePostPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FePostPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FePostPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FePostPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FePostPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FePostPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FePostPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FePostPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FePostPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FePostPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FePostPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FePostPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FePostPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FePostPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FePostPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FePostPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FePostPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FePostPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FePostPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FePostPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FePostPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FePostPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FePostPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FePostPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FeRingwoodite_slb_ph,  std::shared_ptr<FeRingwoodite_slb_ph> >(m, "FeRingwoodite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FeRingwoodite_slb_ph::identifier,"identifier")
    .def("name", &FeRingwoodite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FeRingwoodite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeRingwoodite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeRingwoodite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeRingwoodite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeRingwoodite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FeRingwoodite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FeRingwoodite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FeRingwoodite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FeRingwoodite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FeRingwoodite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FeRingwoodite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FeRingwoodite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FeRingwoodite_slb_ph::endmember_number)
    .def("endmember_name", &FeRingwoodite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FeRingwoodite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FeRingwoodite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FeRingwoodite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FeRingwoodite_slb_ph::species_number)
    .def("species_name", &FeRingwoodite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FeRingwoodite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FeRingwoodite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FeRingwoodite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FeRingwoodite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FeRingwoodite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FeRingwoodite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FeRingwoodite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FeRingwoodite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FeRingwoodite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FeRingwoodite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FeRingwoodite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FeRingwoodite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FeRingwoodite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FeRingwoodite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FeRingwoodite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FeRingwoodite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FeRingwoodite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FeRingwoodite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FeRingwoodite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FeRingwoodite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FeRingwoodite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FeRingwoodite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FeRingwoodite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FeRingwoodite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FeRingwoodite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FeRingwoodite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FeRingwoodite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FeRingwoodite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FeRingwoodite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FeRingwoodite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FeRingwoodite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FeRingwoodite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FeRingwoodite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FeRingwoodite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FeRingwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FeRingwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeRingwoodite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FeRingwoodite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FeRingwoodite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FeRingwoodite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FeRingwoodite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FeRingwoodite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FeRingwoodite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FeRingwoodite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FeRingwoodite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FeRingwoodite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FeRingwoodite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FeRingwoodite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FeRingwoodite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FeRingwoodite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FeRingwoodite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FeRingwoodite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FeRingwoodite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FeRingwoodite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FeRingwoodite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FeRingwoodite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FeRingwoodite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FeRingwoodite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FeRingwoodite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FeRingwoodite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FeRingwoodite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FeRingwoodite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FeRingwoodite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FeRingwoodite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FeRingwoodite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FeRingwoodite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Ferrosilite_slb_ph,  std::shared_ptr<Ferrosilite_slb_ph> >(m, "Ferrosilite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Ferrosilite_slb_ph::identifier,"identifier")
    .def("name", &Ferrosilite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Ferrosilite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Ferrosilite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Ferrosilite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Ferrosilite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Ferrosilite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Ferrosilite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Ferrosilite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Ferrosilite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Ferrosilite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Ferrosilite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Ferrosilite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Ferrosilite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Ferrosilite_slb_ph::endmember_number)
    .def("endmember_name", &Ferrosilite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Ferrosilite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Ferrosilite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Ferrosilite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Ferrosilite_slb_ph::species_number)
    .def("species_name", &Ferrosilite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Ferrosilite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Ferrosilite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Ferrosilite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Ferrosilite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Ferrosilite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Ferrosilite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Ferrosilite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Ferrosilite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Ferrosilite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Ferrosilite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Ferrosilite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Ferrosilite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Ferrosilite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Ferrosilite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Ferrosilite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Ferrosilite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Ferrosilite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Ferrosilite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Ferrosilite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Ferrosilite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Ferrosilite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Ferrosilite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Ferrosilite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Ferrosilite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Ferrosilite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Ferrosilite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Ferrosilite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Ferrosilite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Ferrosilite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Ferrosilite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Ferrosilite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Ferrosilite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Ferrosilite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Ferrosilite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Ferrosilite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Ferrosilite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Ferrosilite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Ferrosilite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Ferrosilite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Ferrosilite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Ferrosilite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Ferrosilite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Ferrosilite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Ferrosilite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Ferrosilite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Ferrosilite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Ferrosilite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Ferrosilite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Ferrosilite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Ferrosilite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Ferrosilite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Ferrosilite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Ferrosilite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Ferrosilite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Ferrosilite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Ferrosilite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Ferrosilite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Ferrosilite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Ferrosilite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Ferrosilite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Ferrosilite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Ferrosilite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Ferrosilite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Ferrosilite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Ferrosilite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Ferrosilite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<FeWadsleyite_slb_ph,  std::shared_ptr<FeWadsleyite_slb_ph> >(m, "FeWadsleyite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &FeWadsleyite_slb_ph::identifier,"identifier")
    .def("name", &FeWadsleyite_slb_ph::name,"phase name")
    .def("tcg_build_version", &FeWadsleyite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &FeWadsleyite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &FeWadsleyite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &FeWadsleyite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &FeWadsleyite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &FeWadsleyite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &FeWadsleyite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &FeWadsleyite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &FeWadsleyite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &FeWadsleyite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &FeWadsleyite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &FeWadsleyite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &FeWadsleyite_slb_ph::endmember_number)
    .def("endmember_name", &FeWadsleyite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &FeWadsleyite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &FeWadsleyite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &FeWadsleyite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &FeWadsleyite_slb_ph::species_number)
    .def("species_name", &FeWadsleyite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &FeWadsleyite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &FeWadsleyite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &FeWadsleyite_slb_ph::species_elements, py::arg("i"))
    .def("g", &FeWadsleyite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &FeWadsleyite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &FeWadsleyite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &FeWadsleyite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &FeWadsleyite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &FeWadsleyite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &FeWadsleyite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &FeWadsleyite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &FeWadsleyite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &FeWadsleyite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &FeWadsleyite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &FeWadsleyite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &FeWadsleyite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &FeWadsleyite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &FeWadsleyite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &FeWadsleyite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &FeWadsleyite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &FeWadsleyite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &FeWadsleyite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &FeWadsleyite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &FeWadsleyite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &FeWadsleyite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &FeWadsleyite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &FeWadsleyite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &FeWadsleyite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &FeWadsleyite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &FeWadsleyite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &FeWadsleyite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &FeWadsleyite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &FeWadsleyite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &FeWadsleyite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &FeWadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &FeWadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&FeWadsleyite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &FeWadsleyite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &FeWadsleyite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &FeWadsleyite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &FeWadsleyite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &FeWadsleyite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &FeWadsleyite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &FeWadsleyite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &FeWadsleyite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &FeWadsleyite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &FeWadsleyite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &FeWadsleyite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &FeWadsleyite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &FeWadsleyite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &FeWadsleyite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &FeWadsleyite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &FeWadsleyite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &FeWadsleyite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &FeWadsleyite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &FeWadsleyite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &FeWadsleyite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &FeWadsleyite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &FeWadsleyite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &FeWadsleyite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &FeWadsleyite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &FeWadsleyite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &FeWadsleyite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &FeWadsleyite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&FeWadsleyite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&FeWadsleyite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Forsterite_slb_ph,  std::shared_ptr<Forsterite_slb_ph> >(m, "Forsterite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Forsterite_slb_ph::identifier,"identifier")
    .def("name", &Forsterite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Forsterite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Forsterite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Forsterite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Forsterite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Forsterite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Forsterite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Forsterite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Forsterite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Forsterite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Forsterite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Forsterite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Forsterite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Forsterite_slb_ph::endmember_number)
    .def("endmember_name", &Forsterite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Forsterite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Forsterite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Forsterite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Forsterite_slb_ph::species_number)
    .def("species_name", &Forsterite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Forsterite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Forsterite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Forsterite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Forsterite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Forsterite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Forsterite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Forsterite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Forsterite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Forsterite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Forsterite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Forsterite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Forsterite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Forsterite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Forsterite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Forsterite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Forsterite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Forsterite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Forsterite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Forsterite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Forsterite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Forsterite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Forsterite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Forsterite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Forsterite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Forsterite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Forsterite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Forsterite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Forsterite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Forsterite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Forsterite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Forsterite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Forsterite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Forsterite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Forsterite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Forsterite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Forsterite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Forsterite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Forsterite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Forsterite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Forsterite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Forsterite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Forsterite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Forsterite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Forsterite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Forsterite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Forsterite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Forsterite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Forsterite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Forsterite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Forsterite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Forsterite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Forsterite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Forsterite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Forsterite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Forsterite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Forsterite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Forsterite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Forsterite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Forsterite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Forsterite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Forsterite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Forsterite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Forsterite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Forsterite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Forsterite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Forsterite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Garnet_slb_ph,  std::shared_ptr<Garnet_slb_ph> >(m, "Garnet_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Garnet_slb_ph::identifier,"identifier")
    .def("name", &Garnet_slb_ph::name,"phase name")
    .def("tcg_build_version", &Garnet_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Garnet_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Garnet_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Garnet_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Garnet_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Garnet_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Garnet_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Garnet_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Garnet_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Garnet_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Garnet_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Garnet_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Garnet_slb_ph::endmember_number)
    .def("endmember_name", &Garnet_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Garnet_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Garnet_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Garnet_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Garnet_slb_ph::species_number)
    .def("species_name", &Garnet_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Garnet_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Garnet_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Garnet_slb_ph::species_elements, py::arg("i"))
    .def("g", &Garnet_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Garnet_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Garnet_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Garnet_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Garnet_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Garnet_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Garnet_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Garnet_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Garnet_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Garnet_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Garnet_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Garnet_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Garnet_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Garnet_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Garnet_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Garnet_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Garnet_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Garnet_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Garnet_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Garnet_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Garnet_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Garnet_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Garnet_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Garnet_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Garnet_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Garnet_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Garnet_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Garnet_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Garnet_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Garnet_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Garnet_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Garnet_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Garnet_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Garnet_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Garnet_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Garnet_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Garnet_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Garnet_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Garnet_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Garnet_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Garnet_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Garnet_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Garnet_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Garnet_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Garnet_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Garnet_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Garnet_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Garnet_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Garnet_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Garnet_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Garnet_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Garnet_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Garnet_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Garnet_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Garnet_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Garnet_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Garnet_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Garnet_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Garnet_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Garnet_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Garnet_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Garnet_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Garnet_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Grossular_slb_ph,  std::shared_ptr<Grossular_slb_ph> >(m, "Grossular_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Grossular_slb_ph::identifier,"identifier")
    .def("name", &Grossular_slb_ph::name,"phase name")
    .def("tcg_build_version", &Grossular_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Grossular_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Grossular_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Grossular_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Grossular_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Grossular_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Grossular_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Grossular_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Grossular_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Grossular_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Grossular_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Grossular_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Grossular_slb_ph::endmember_number)
    .def("endmember_name", &Grossular_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Grossular_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Grossular_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Grossular_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Grossular_slb_ph::species_number)
    .def("species_name", &Grossular_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Grossular_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Grossular_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Grossular_slb_ph::species_elements, py::arg("i"))
    .def("g", &Grossular_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Grossular_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Grossular_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Grossular_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Grossular_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Grossular_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Grossular_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Grossular_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Grossular_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Grossular_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Grossular_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Grossular_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Grossular_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Grossular_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Grossular_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Grossular_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Grossular_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Grossular_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Grossular_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Grossular_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Grossular_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Grossular_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Grossular_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Grossular_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Grossular_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Grossular_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Grossular_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Grossular_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Grossular_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Grossular_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Grossular_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Grossular_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Grossular_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Grossular_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Grossular_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Grossular_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Grossular_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Grossular_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Grossular_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Grossular_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Grossular_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Grossular_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Grossular_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Grossular_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Grossular_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Grossular_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Grossular_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Grossular_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Grossular_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Grossular_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Grossular_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Grossular_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Grossular_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Grossular_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Grossular_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Grossular_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Grossular_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Grossular_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Grossular_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Grossular_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Grossular_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Grossular_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Grossular_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Hedenbergite_slb_ph,  std::shared_ptr<Hedenbergite_slb_ph> >(m, "Hedenbergite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Hedenbergite_slb_ph::identifier,"identifier")
    .def("name", &Hedenbergite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Hedenbergite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Hedenbergite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Hedenbergite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Hedenbergite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Hedenbergite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Hedenbergite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Hedenbergite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Hedenbergite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Hedenbergite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Hedenbergite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Hedenbergite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Hedenbergite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Hedenbergite_slb_ph::endmember_number)
    .def("endmember_name", &Hedenbergite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Hedenbergite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Hedenbergite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Hedenbergite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Hedenbergite_slb_ph::species_number)
    .def("species_name", &Hedenbergite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Hedenbergite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Hedenbergite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Hedenbergite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Hedenbergite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Hedenbergite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Hedenbergite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Hedenbergite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Hedenbergite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Hedenbergite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Hedenbergite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Hedenbergite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Hedenbergite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Hedenbergite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Hedenbergite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Hedenbergite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Hedenbergite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Hedenbergite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Hedenbergite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Hedenbergite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Hedenbergite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Hedenbergite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Hedenbergite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Hedenbergite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Hedenbergite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Hedenbergite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Hedenbergite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Hedenbergite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Hedenbergite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Hedenbergite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Hedenbergite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Hedenbergite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Hedenbergite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Hedenbergite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Hedenbergite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Hedenbergite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Hedenbergite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Hedenbergite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Hedenbergite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Hedenbergite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Hedenbergite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Hedenbergite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Hedenbergite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Hedenbergite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Hedenbergite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Hedenbergite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Hedenbergite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Hedenbergite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Hedenbergite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Hedenbergite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Hedenbergite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Hedenbergite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Hedenbergite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Hedenbergite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Hedenbergite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Hedenbergite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Hedenbergite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Hedenbergite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Hedenbergite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Hedenbergite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Hedenbergite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Hedenbergite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Hedenbergite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Hedenbergite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Hedenbergite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Hedenbergite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Hedenbergite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Hercynite_slb_ph,  std::shared_ptr<Hercynite_slb_ph> >(m, "Hercynite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Hercynite_slb_ph::identifier,"identifier")
    .def("name", &Hercynite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Hercynite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Hercynite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Hercynite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Hercynite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Hercynite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Hercynite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Hercynite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Hercynite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Hercynite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Hercynite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Hercynite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Hercynite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Hercynite_slb_ph::endmember_number)
    .def("endmember_name", &Hercynite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Hercynite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Hercynite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Hercynite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Hercynite_slb_ph::species_number)
    .def("species_name", &Hercynite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Hercynite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Hercynite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Hercynite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Hercynite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Hercynite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Hercynite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Hercynite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Hercynite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Hercynite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Hercynite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Hercynite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Hercynite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Hercynite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Hercynite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Hercynite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Hercynite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Hercynite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Hercynite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Hercynite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Hercynite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Hercynite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Hercynite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Hercynite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Hercynite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Hercynite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Hercynite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Hercynite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Hercynite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Hercynite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Hercynite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Hercynite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Hercynite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Hercynite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Hercynite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Hercynite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Hercynite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Hercynite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Hercynite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Hercynite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Hercynite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Hercynite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Hercynite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Hercynite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Hercynite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Hercynite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Hercynite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Hercynite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Hercynite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Hercynite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Hercynite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Hercynite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Hercynite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Hercynite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Hercynite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Hercynite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Hercynite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Hercynite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Hercynite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Hercynite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Hercynite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Hercynite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Hercynite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Hercynite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Hercynite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Hercynite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Hercynite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<HPClinoenstatite_slb_ph,  std::shared_ptr<HPClinoenstatite_slb_ph> >(m, "HPClinoenstatite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &HPClinoenstatite_slb_ph::identifier,"identifier")
    .def("name", &HPClinoenstatite_slb_ph::name,"phase name")
    .def("tcg_build_version", &HPClinoenstatite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &HPClinoenstatite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &HPClinoenstatite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &HPClinoenstatite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &HPClinoenstatite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &HPClinoenstatite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &HPClinoenstatite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &HPClinoenstatite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &HPClinoenstatite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &HPClinoenstatite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &HPClinoenstatite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &HPClinoenstatite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &HPClinoenstatite_slb_ph::endmember_number)
    .def("endmember_name", &HPClinoenstatite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &HPClinoenstatite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &HPClinoenstatite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &HPClinoenstatite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &HPClinoenstatite_slb_ph::species_number)
    .def("species_name", &HPClinoenstatite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &HPClinoenstatite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &HPClinoenstatite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &HPClinoenstatite_slb_ph::species_elements, py::arg("i"))
    .def("g", &HPClinoenstatite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &HPClinoenstatite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &HPClinoenstatite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &HPClinoenstatite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &HPClinoenstatite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &HPClinoenstatite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &HPClinoenstatite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &HPClinoenstatite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &HPClinoenstatite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &HPClinoenstatite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &HPClinoenstatite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &HPClinoenstatite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &HPClinoenstatite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &HPClinoenstatite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &HPClinoenstatite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &HPClinoenstatite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &HPClinoenstatite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &HPClinoenstatite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &HPClinoenstatite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &HPClinoenstatite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &HPClinoenstatite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &HPClinoenstatite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &HPClinoenstatite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &HPClinoenstatite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &HPClinoenstatite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &HPClinoenstatite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &HPClinoenstatite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &HPClinoenstatite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &HPClinoenstatite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &HPClinoenstatite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &HPClinoenstatite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &HPClinoenstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &HPClinoenstatite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&HPClinoenstatite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &HPClinoenstatite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &HPClinoenstatite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &HPClinoenstatite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &HPClinoenstatite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &HPClinoenstatite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &HPClinoenstatite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &HPClinoenstatite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &HPClinoenstatite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &HPClinoenstatite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &HPClinoenstatite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &HPClinoenstatite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &HPClinoenstatite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &HPClinoenstatite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &HPClinoenstatite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &HPClinoenstatite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &HPClinoenstatite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &HPClinoenstatite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &HPClinoenstatite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &HPClinoenstatite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &HPClinoenstatite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &HPClinoenstatite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &HPClinoenstatite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &HPClinoenstatite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &HPClinoenstatite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &HPClinoenstatite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &HPClinoenstatite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &HPClinoenstatite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&HPClinoenstatite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&HPClinoenstatite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<HPClinoferrosilite_slb_ph,  std::shared_ptr<HPClinoferrosilite_slb_ph> >(m, "HPClinoferrosilite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &HPClinoferrosilite_slb_ph::identifier,"identifier")
    .def("name", &HPClinoferrosilite_slb_ph::name,"phase name")
    .def("tcg_build_version", &HPClinoferrosilite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &HPClinoferrosilite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &HPClinoferrosilite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &HPClinoferrosilite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &HPClinoferrosilite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &HPClinoferrosilite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &HPClinoferrosilite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &HPClinoferrosilite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &HPClinoferrosilite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &HPClinoferrosilite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &HPClinoferrosilite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &HPClinoferrosilite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &HPClinoferrosilite_slb_ph::endmember_number)
    .def("endmember_name", &HPClinoferrosilite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &HPClinoferrosilite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &HPClinoferrosilite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &HPClinoferrosilite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &HPClinoferrosilite_slb_ph::species_number)
    .def("species_name", &HPClinoferrosilite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &HPClinoferrosilite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &HPClinoferrosilite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &HPClinoferrosilite_slb_ph::species_elements, py::arg("i"))
    .def("g", &HPClinoferrosilite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &HPClinoferrosilite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &HPClinoferrosilite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &HPClinoferrosilite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &HPClinoferrosilite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &HPClinoferrosilite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &HPClinoferrosilite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &HPClinoferrosilite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &HPClinoferrosilite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &HPClinoferrosilite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &HPClinoferrosilite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &HPClinoferrosilite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &HPClinoferrosilite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &HPClinoferrosilite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &HPClinoferrosilite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &HPClinoferrosilite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &HPClinoferrosilite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &HPClinoferrosilite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &HPClinoferrosilite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &HPClinoferrosilite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &HPClinoferrosilite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &HPClinoferrosilite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &HPClinoferrosilite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &HPClinoferrosilite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &HPClinoferrosilite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &HPClinoferrosilite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &HPClinoferrosilite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &HPClinoferrosilite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &HPClinoferrosilite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &HPClinoferrosilite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &HPClinoferrosilite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &HPClinoferrosilite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &HPClinoferrosilite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&HPClinoferrosilite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &HPClinoferrosilite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &HPClinoferrosilite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &HPClinoferrosilite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &HPClinoferrosilite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &HPClinoferrosilite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &HPClinoferrosilite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &HPClinoferrosilite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &HPClinoferrosilite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &HPClinoferrosilite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &HPClinoferrosilite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &HPClinoferrosilite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &HPClinoferrosilite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &HPClinoferrosilite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &HPClinoferrosilite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &HPClinoferrosilite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &HPClinoferrosilite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &HPClinoferrosilite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &HPClinoferrosilite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &HPClinoferrosilite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &HPClinoferrosilite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &HPClinoferrosilite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &HPClinoferrosilite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &HPClinoferrosilite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &HPClinoferrosilite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &HPClinoferrosilite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &HPClinoferrosilite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &HPClinoferrosilite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&HPClinoferrosilite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&HPClinoferrosilite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<HPClinopyroxene_slb_ph,  std::shared_ptr<HPClinopyroxene_slb_ph> >(m, "HPClinopyroxene_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &HPClinopyroxene_slb_ph::identifier,"identifier")
    .def("name", &HPClinopyroxene_slb_ph::name,"phase name")
    .def("tcg_build_version", &HPClinopyroxene_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &HPClinopyroxene_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &HPClinopyroxene_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &HPClinopyroxene_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &HPClinopyroxene_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &HPClinopyroxene_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &HPClinopyroxene_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &HPClinopyroxene_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &HPClinopyroxene_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &HPClinopyroxene_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &HPClinopyroxene_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &HPClinopyroxene_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &HPClinopyroxene_slb_ph::endmember_number)
    .def("endmember_name", &HPClinopyroxene_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &HPClinopyroxene_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &HPClinopyroxene_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &HPClinopyroxene_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &HPClinopyroxene_slb_ph::species_number)
    .def("species_name", &HPClinopyroxene_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &HPClinopyroxene_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &HPClinopyroxene_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &HPClinopyroxene_slb_ph::species_elements, py::arg("i"))
    .def("g", &HPClinopyroxene_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &HPClinopyroxene_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &HPClinopyroxene_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &HPClinopyroxene_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &HPClinopyroxene_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &HPClinopyroxene_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &HPClinopyroxene_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &HPClinopyroxene_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &HPClinopyroxene_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &HPClinopyroxene_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &HPClinopyroxene_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &HPClinopyroxene_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &HPClinopyroxene_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &HPClinopyroxene_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &HPClinopyroxene_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &HPClinopyroxene_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &HPClinopyroxene_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &HPClinopyroxene_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &HPClinopyroxene_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &HPClinopyroxene_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &HPClinopyroxene_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &HPClinopyroxene_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &HPClinopyroxene_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &HPClinopyroxene_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &HPClinopyroxene_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &HPClinopyroxene_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &HPClinopyroxene_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &HPClinopyroxene_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &HPClinopyroxene_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &HPClinopyroxene_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &HPClinopyroxene_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &HPClinopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &HPClinopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&HPClinopyroxene_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &HPClinopyroxene_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &HPClinopyroxene_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &HPClinopyroxene_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &HPClinopyroxene_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &HPClinopyroxene_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &HPClinopyroxene_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &HPClinopyroxene_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &HPClinopyroxene_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &HPClinopyroxene_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &HPClinopyroxene_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &HPClinopyroxene_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &HPClinopyroxene_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &HPClinopyroxene_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &HPClinopyroxene_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &HPClinopyroxene_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &HPClinopyroxene_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &HPClinopyroxene_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &HPClinopyroxene_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &HPClinopyroxene_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &HPClinopyroxene_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &HPClinopyroxene_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &HPClinopyroxene_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &HPClinopyroxene_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &HPClinopyroxene_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &HPClinopyroxene_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &HPClinopyroxene_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &HPClinopyroxene_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&HPClinopyroxene_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&HPClinopyroxene_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Jadeite_slb_ph,  std::shared_ptr<Jadeite_slb_ph> >(m, "Jadeite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Jadeite_slb_ph::identifier,"identifier")
    .def("name", &Jadeite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Jadeite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Jadeite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Jadeite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Jadeite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Jadeite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Jadeite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Jadeite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Jadeite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Jadeite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Jadeite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Jadeite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Jadeite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Jadeite_slb_ph::endmember_number)
    .def("endmember_name", &Jadeite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Jadeite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Jadeite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Jadeite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Jadeite_slb_ph::species_number)
    .def("species_name", &Jadeite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Jadeite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Jadeite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Jadeite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Jadeite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Jadeite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Jadeite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Jadeite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Jadeite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Jadeite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Jadeite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Jadeite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Jadeite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Jadeite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Jadeite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Jadeite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Jadeite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Jadeite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Jadeite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Jadeite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Jadeite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Jadeite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Jadeite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Jadeite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Jadeite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Jadeite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Jadeite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Jadeite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Jadeite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Jadeite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Jadeite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Jadeite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Jadeite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Jadeite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Jadeite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Jadeite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Jadeite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Jadeite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Jadeite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Jadeite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Jadeite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Jadeite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Jadeite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Jadeite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Jadeite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Jadeite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Jadeite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Jadeite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Jadeite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Jadeite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Jadeite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Jadeite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Jadeite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Jadeite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Jadeite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Jadeite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Jadeite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Jadeite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Jadeite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Jadeite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Jadeite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Jadeite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Jadeite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Jadeite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Jadeite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Jadeite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Jadeite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Kyanite_slb_ph,  std::shared_ptr<Kyanite_slb_ph> >(m, "Kyanite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Kyanite_slb_ph::identifier,"identifier")
    .def("name", &Kyanite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Kyanite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Kyanite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Kyanite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Kyanite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Kyanite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Kyanite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Kyanite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Kyanite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Kyanite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Kyanite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Kyanite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Kyanite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Kyanite_slb_ph::endmember_number)
    .def("endmember_name", &Kyanite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Kyanite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Kyanite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Kyanite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Kyanite_slb_ph::species_number)
    .def("species_name", &Kyanite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Kyanite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Kyanite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Kyanite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Kyanite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Kyanite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Kyanite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Kyanite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Kyanite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Kyanite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Kyanite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Kyanite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Kyanite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Kyanite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Kyanite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Kyanite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Kyanite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Kyanite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Kyanite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Kyanite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Kyanite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Kyanite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Kyanite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Kyanite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Kyanite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Kyanite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Kyanite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Kyanite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Kyanite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Kyanite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Kyanite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Kyanite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Kyanite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Kyanite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Kyanite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Kyanite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Kyanite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Kyanite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Kyanite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Kyanite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Kyanite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Kyanite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Kyanite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Kyanite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Kyanite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Kyanite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Kyanite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Kyanite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Kyanite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Kyanite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Kyanite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Kyanite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Kyanite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Kyanite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Kyanite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Kyanite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Kyanite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Kyanite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Kyanite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Kyanite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Kyanite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Kyanite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Kyanite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Kyanite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Kyanite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Kyanite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Kyanite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Magnesiowuestite_slb_ph,  std::shared_ptr<Magnesiowuestite_slb_ph> >(m, "Magnesiowuestite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Magnesiowuestite_slb_ph::identifier,"identifier")
    .def("name", &Magnesiowuestite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Magnesiowuestite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Magnesiowuestite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Magnesiowuestite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Magnesiowuestite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Magnesiowuestite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Magnesiowuestite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Magnesiowuestite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Magnesiowuestite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Magnesiowuestite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Magnesiowuestite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Magnesiowuestite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Magnesiowuestite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Magnesiowuestite_slb_ph::endmember_number)
    .def("endmember_name", &Magnesiowuestite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Magnesiowuestite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Magnesiowuestite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Magnesiowuestite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Magnesiowuestite_slb_ph::species_number)
    .def("species_name", &Magnesiowuestite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Magnesiowuestite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Magnesiowuestite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Magnesiowuestite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Magnesiowuestite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Magnesiowuestite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Magnesiowuestite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Magnesiowuestite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Magnesiowuestite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Magnesiowuestite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Magnesiowuestite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Magnesiowuestite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Magnesiowuestite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Magnesiowuestite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Magnesiowuestite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Magnesiowuestite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Magnesiowuestite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Magnesiowuestite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Magnesiowuestite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Magnesiowuestite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Magnesiowuestite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Magnesiowuestite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Magnesiowuestite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Magnesiowuestite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Magnesiowuestite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Magnesiowuestite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Magnesiowuestite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Magnesiowuestite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Magnesiowuestite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Magnesiowuestite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Magnesiowuestite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Magnesiowuestite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Magnesiowuestite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Magnesiowuestite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Magnesiowuestite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Magnesiowuestite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Magnesiowuestite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Magnesiowuestite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Magnesiowuestite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Magnesiowuestite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Magnesiowuestite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Magnesiowuestite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Magnesiowuestite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Magnesiowuestite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Magnesiowuestite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Magnesiowuestite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Magnesiowuestite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Magnesiowuestite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Magnesiowuestite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Magnesiowuestite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Magnesiowuestite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Magnesiowuestite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Magnesiowuestite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Magnesiowuestite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Magnesiowuestite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Magnesiowuestite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Magnesiowuestite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Magnesiowuestite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Magnesiowuestite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Magnesiowuestite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Magnesiowuestite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Magnesiowuestite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Magnesiowuestite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Magnesiowuestite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Magnesiowuestite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Magnesiowuestite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Magnesiowuestite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgAkimotoite_slb_ph,  std::shared_ptr<MgAkimotoite_slb_ph> >(m, "MgAkimotoite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgAkimotoite_slb_ph::identifier,"identifier")
    .def("name", &MgAkimotoite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgAkimotoite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgAkimotoite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgAkimotoite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgAkimotoite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgAkimotoite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgAkimotoite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgAkimotoite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgAkimotoite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgAkimotoite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgAkimotoite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgAkimotoite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgAkimotoite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgAkimotoite_slb_ph::endmember_number)
    .def("endmember_name", &MgAkimotoite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgAkimotoite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgAkimotoite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgAkimotoite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgAkimotoite_slb_ph::species_number)
    .def("species_name", &MgAkimotoite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgAkimotoite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgAkimotoite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgAkimotoite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgAkimotoite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgAkimotoite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgAkimotoite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgAkimotoite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgAkimotoite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgAkimotoite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgAkimotoite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgAkimotoite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgAkimotoite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgAkimotoite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgAkimotoite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgAkimotoite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgAkimotoite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgAkimotoite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgAkimotoite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgAkimotoite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgAkimotoite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgAkimotoite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgAkimotoite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgAkimotoite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgAkimotoite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgAkimotoite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgAkimotoite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgAkimotoite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgAkimotoite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgAkimotoite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgAkimotoite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgAkimotoite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgAkimotoite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgAkimotoite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgAkimotoite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgAkimotoite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgAkimotoite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgAkimotoite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgAkimotoite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgAkimotoite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgAkimotoite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgAkimotoite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgAkimotoite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgAkimotoite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgAkimotoite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgAkimotoite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgAkimotoite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgAkimotoite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgAkimotoite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgAkimotoite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgAkimotoite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgAkimotoite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgAkimotoite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgAkimotoite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgAkimotoite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgAkimotoite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgAkimotoite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgAkimotoite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgAkimotoite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgAkimotoite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgAkimotoite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgAkimotoite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgAkimotoite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgAkimotoite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgAkimotoite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgAkimotoite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgCaFerrite_slb_ph,  std::shared_ptr<MgCaFerrite_slb_ph> >(m, "MgCaFerrite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgCaFerrite_slb_ph::identifier,"identifier")
    .def("name", &MgCaFerrite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgCaFerrite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgCaFerrite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgCaFerrite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgCaFerrite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgCaFerrite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgCaFerrite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgCaFerrite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgCaFerrite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgCaFerrite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgCaFerrite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgCaFerrite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgCaFerrite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgCaFerrite_slb_ph::endmember_number)
    .def("endmember_name", &MgCaFerrite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgCaFerrite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgCaFerrite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgCaFerrite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgCaFerrite_slb_ph::species_number)
    .def("species_name", &MgCaFerrite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgCaFerrite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgCaFerrite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgCaFerrite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgCaFerrite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgCaFerrite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgCaFerrite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgCaFerrite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgCaFerrite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgCaFerrite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgCaFerrite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgCaFerrite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgCaFerrite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgCaFerrite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgCaFerrite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgCaFerrite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgCaFerrite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgCaFerrite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgCaFerrite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgCaFerrite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgCaFerrite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgCaFerrite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgCaFerrite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgCaFerrite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgCaFerrite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgCaFerrite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgCaFerrite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgCaFerrite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgCaFerrite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgCaFerrite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgCaFerrite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgCaFerrite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgCaFerrite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgCaFerrite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgCaFerrite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgCaFerrite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgCaFerrite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgCaFerrite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgCaFerrite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgCaFerrite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgCaFerrite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgCaFerrite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgCaFerrite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgCaFerrite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgCaFerrite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgCaFerrite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgCaFerrite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgCaFerrite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgCaFerrite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgCaFerrite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgCaFerrite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgCaFerrite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgCaFerrite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgCaFerrite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgCaFerrite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgCaFerrite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgCaFerrite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgCaFerrite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgCaFerrite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgCaFerrite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgCaFerrite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgCaFerrite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgCaFerrite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgCaFerrite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgCaFerrite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgFePerovskite_slb_ph,  std::shared_ptr<MgFePerovskite_slb_ph> >(m, "MgFePerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgFePerovskite_slb_ph::identifier,"identifier")
    .def("name", &MgFePerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgFePerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgFePerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgFePerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgFePerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgFePerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgFePerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgFePerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgFePerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgFePerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgFePerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgFePerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgFePerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgFePerovskite_slb_ph::endmember_number)
    .def("endmember_name", &MgFePerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgFePerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgFePerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgFePerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgFePerovskite_slb_ph::species_number)
    .def("species_name", &MgFePerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgFePerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgFePerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgFePerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgFePerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgFePerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgFePerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgFePerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgFePerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgFePerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgFePerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgFePerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgFePerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgFePerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgFePerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgFePerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgFePerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgFePerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgFePerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgFePerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgFePerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgFePerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgFePerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgFePerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgFePerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgFePerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgFePerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgFePerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgFePerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgFePerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgFePerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgFePerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgFePerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgFePerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgFePerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgFePerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgFePerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgFePerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgFePerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgFePerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgFePerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgFePerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgFePerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgFePerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgFePerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgFePerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgFePerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgFePerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgFePerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgFePerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgFePerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgFePerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgFePerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgFePerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgFePerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgFePerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgFePerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgFePerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgFePerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgFePerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgFePerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgFePerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgFePerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgFePerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgFePerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgFePerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgFePerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgMajorite_slb_ph,  std::shared_ptr<MgMajorite_slb_ph> >(m, "MgMajorite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgMajorite_slb_ph::identifier,"identifier")
    .def("name", &MgMajorite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgMajorite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgMajorite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgMajorite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgMajorite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgMajorite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgMajorite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgMajorite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgMajorite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgMajorite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgMajorite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgMajorite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgMajorite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgMajorite_slb_ph::endmember_number)
    .def("endmember_name", &MgMajorite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgMajorite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgMajorite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgMajorite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgMajorite_slb_ph::species_number)
    .def("species_name", &MgMajorite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgMajorite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgMajorite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgMajorite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgMajorite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgMajorite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgMajorite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgMajorite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgMajorite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgMajorite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgMajorite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgMajorite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgMajorite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgMajorite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgMajorite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgMajorite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgMajorite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgMajorite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgMajorite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgMajorite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgMajorite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgMajorite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgMajorite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgMajorite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgMajorite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgMajorite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgMajorite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgMajorite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgMajorite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgMajorite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgMajorite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgMajorite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgMajorite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgMajorite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgMajorite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgMajorite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgMajorite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgMajorite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgMajorite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgMajorite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgMajorite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgMajorite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgMajorite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgMajorite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgMajorite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgMajorite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgMajorite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgMajorite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgMajorite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgMajorite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgMajorite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgMajorite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgMajorite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgMajorite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgMajorite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgMajorite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgMajorite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgMajorite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgMajorite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgMajorite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgMajorite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgMajorite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgMajorite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgMajorite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgMajorite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgMajorite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgMajorite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgPerovskite_slb_ph,  std::shared_ptr<MgPerovskite_slb_ph> >(m, "MgPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgPerovskite_slb_ph::identifier,"identifier")
    .def("name", &MgPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &MgPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgPerovskite_slb_ph::species_number)
    .def("species_name", &MgPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgPostPerovskite_slb_ph,  std::shared_ptr<MgPostPerovskite_slb_ph> >(m, "MgPostPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgPostPerovskite_slb_ph::identifier,"identifier")
    .def("name", &MgPostPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgPostPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgPostPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgPostPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgPostPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgPostPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgPostPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgPostPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgPostPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgPostPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgPostPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgPostPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgPostPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgPostPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &MgPostPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgPostPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgPostPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgPostPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgPostPerovskite_slb_ph::species_number)
    .def("species_name", &MgPostPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgPostPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgPostPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgPostPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgPostPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgPostPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgPostPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgPostPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgPostPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgPostPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgPostPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgPostPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgPostPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgPostPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgPostPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgPostPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgPostPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgPostPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgPostPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgPostPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgPostPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgPostPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgPostPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgPostPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgPostPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgPostPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgPostPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgPostPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgPostPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgPostPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgPostPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgPostPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgPostPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgPostPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgPostPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgPostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgPostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgPostPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgPostPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgPostPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgPostPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgPostPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgPostPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgPostPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgPostPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgPostPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgPostPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgPostPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgPostPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgPostPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgPostPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgPostPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgPostPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgPostPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgPostPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgPostPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgPostPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgPostPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgPostPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgPostPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgPostPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgPostPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgPostPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgPostPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgPostPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgPostPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgPostPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgRingwoodite_slb_ph,  std::shared_ptr<MgRingwoodite_slb_ph> >(m, "MgRingwoodite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgRingwoodite_slb_ph::identifier,"identifier")
    .def("name", &MgRingwoodite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgRingwoodite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgRingwoodite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgRingwoodite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgRingwoodite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgRingwoodite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgRingwoodite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgRingwoodite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgRingwoodite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgRingwoodite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgRingwoodite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgRingwoodite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgRingwoodite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgRingwoodite_slb_ph::endmember_number)
    .def("endmember_name", &MgRingwoodite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgRingwoodite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgRingwoodite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgRingwoodite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgRingwoodite_slb_ph::species_number)
    .def("species_name", &MgRingwoodite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgRingwoodite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgRingwoodite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgRingwoodite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgRingwoodite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgRingwoodite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgRingwoodite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgRingwoodite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgRingwoodite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgRingwoodite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgRingwoodite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgRingwoodite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgRingwoodite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgRingwoodite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgRingwoodite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgRingwoodite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgRingwoodite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgRingwoodite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgRingwoodite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgRingwoodite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgRingwoodite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgRingwoodite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgRingwoodite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgRingwoodite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgRingwoodite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgRingwoodite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgRingwoodite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgRingwoodite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgRingwoodite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgRingwoodite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgRingwoodite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgRingwoodite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgRingwoodite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgRingwoodite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgRingwoodite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgRingwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgRingwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgRingwoodite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgRingwoodite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgRingwoodite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgRingwoodite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgRingwoodite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgRingwoodite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgRingwoodite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgRingwoodite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgRingwoodite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgRingwoodite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgRingwoodite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgRingwoodite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgRingwoodite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgRingwoodite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgRingwoodite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgRingwoodite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgRingwoodite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgRingwoodite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgRingwoodite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgRingwoodite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgRingwoodite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgRingwoodite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgRingwoodite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgRingwoodite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgRingwoodite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgRingwoodite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgRingwoodite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgRingwoodite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgRingwoodite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgRingwoodite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgSpinel_slb_ph,  std::shared_ptr<MgSpinel_slb_ph> >(m, "MgSpinel_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgSpinel_slb_ph::identifier,"identifier")
    .def("name", &MgSpinel_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgSpinel_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgSpinel_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgSpinel_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgSpinel_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgSpinel_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgSpinel_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgSpinel_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgSpinel_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgSpinel_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgSpinel_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgSpinel_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgSpinel_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgSpinel_slb_ph::endmember_number)
    .def("endmember_name", &MgSpinel_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgSpinel_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgSpinel_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgSpinel_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgSpinel_slb_ph::species_number)
    .def("species_name", &MgSpinel_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgSpinel_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgSpinel_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgSpinel_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgSpinel_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgSpinel_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgSpinel_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgSpinel_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgSpinel_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgSpinel_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgSpinel_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgSpinel_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgSpinel_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgSpinel_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgSpinel_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgSpinel_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgSpinel_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgSpinel_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgSpinel_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgSpinel_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgSpinel_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgSpinel_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgSpinel_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgSpinel_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgSpinel_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgSpinel_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgSpinel_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgSpinel_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgSpinel_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgSpinel_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgSpinel_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgSpinel_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgSpinel_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgSpinel_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgSpinel_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgSpinel_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgSpinel_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgSpinel_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgSpinel_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgSpinel_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgSpinel_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgSpinel_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgSpinel_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgSpinel_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgSpinel_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgSpinel_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgSpinel_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgSpinel_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgSpinel_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgSpinel_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgSpinel_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgSpinel_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgSpinel_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgSpinel_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgSpinel_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgSpinel_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgSpinel_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgSpinel_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgSpinel_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgSpinel_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgSpinel_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgSpinel_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgSpinel_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgSpinel_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgSpinel_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgSpinel_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgSpinel_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgTschermaks_slb_ph,  std::shared_ptr<MgTschermaks_slb_ph> >(m, "MgTschermaks_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgTschermaks_slb_ph::identifier,"identifier")
    .def("name", &MgTschermaks_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgTschermaks_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgTschermaks_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgTschermaks_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgTschermaks_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgTschermaks_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgTschermaks_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgTschermaks_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgTschermaks_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgTschermaks_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgTschermaks_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgTschermaks_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgTschermaks_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgTschermaks_slb_ph::endmember_number)
    .def("endmember_name", &MgTschermaks_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgTschermaks_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgTschermaks_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgTschermaks_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgTschermaks_slb_ph::species_number)
    .def("species_name", &MgTschermaks_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgTschermaks_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgTschermaks_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgTschermaks_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgTschermaks_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgTschermaks_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgTschermaks_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgTschermaks_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgTschermaks_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgTschermaks_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgTschermaks_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgTschermaks_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgTschermaks_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgTschermaks_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgTschermaks_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgTschermaks_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgTschermaks_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgTschermaks_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgTschermaks_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgTschermaks_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgTschermaks_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgTschermaks_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgTschermaks_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgTschermaks_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgTschermaks_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgTschermaks_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgTschermaks_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgTschermaks_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgTschermaks_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgTschermaks_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgTschermaks_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgTschermaks_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgTschermaks_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgTschermaks_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgTschermaks_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgTschermaks_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgTschermaks_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgTschermaks_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgTschermaks_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgTschermaks_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgTschermaks_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgTschermaks_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgTschermaks_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgTschermaks_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgTschermaks_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgTschermaks_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgTschermaks_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgTschermaks_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgTschermaks_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgTschermaks_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgTschermaks_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgTschermaks_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgTschermaks_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgTschermaks_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgTschermaks_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgTschermaks_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgTschermaks_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgTschermaks_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgTschermaks_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgTschermaks_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgTschermaks_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgTschermaks_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgTschermaks_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgTschermaks_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgTschermaks_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgTschermaks_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgTschermaks_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<MgWadsleyite_slb_ph,  std::shared_ptr<MgWadsleyite_slb_ph> >(m, "MgWadsleyite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &MgWadsleyite_slb_ph::identifier,"identifier")
    .def("name", &MgWadsleyite_slb_ph::name,"phase name")
    .def("tcg_build_version", &MgWadsleyite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &MgWadsleyite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &MgWadsleyite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &MgWadsleyite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &MgWadsleyite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &MgWadsleyite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &MgWadsleyite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &MgWadsleyite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &MgWadsleyite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &MgWadsleyite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &MgWadsleyite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &MgWadsleyite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &MgWadsleyite_slb_ph::endmember_number)
    .def("endmember_name", &MgWadsleyite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &MgWadsleyite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &MgWadsleyite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &MgWadsleyite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &MgWadsleyite_slb_ph::species_number)
    .def("species_name", &MgWadsleyite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &MgWadsleyite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &MgWadsleyite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &MgWadsleyite_slb_ph::species_elements, py::arg("i"))
    .def("g", &MgWadsleyite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &MgWadsleyite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &MgWadsleyite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &MgWadsleyite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &MgWadsleyite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &MgWadsleyite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &MgWadsleyite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &MgWadsleyite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &MgWadsleyite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &MgWadsleyite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &MgWadsleyite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &MgWadsleyite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &MgWadsleyite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &MgWadsleyite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &MgWadsleyite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &MgWadsleyite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &MgWadsleyite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &MgWadsleyite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &MgWadsleyite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &MgWadsleyite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &MgWadsleyite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &MgWadsleyite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &MgWadsleyite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &MgWadsleyite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &MgWadsleyite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &MgWadsleyite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &MgWadsleyite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &MgWadsleyite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &MgWadsleyite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &MgWadsleyite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &MgWadsleyite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &MgWadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &MgWadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&MgWadsleyite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &MgWadsleyite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &MgWadsleyite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &MgWadsleyite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &MgWadsleyite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &MgWadsleyite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &MgWadsleyite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &MgWadsleyite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &MgWadsleyite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &MgWadsleyite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &MgWadsleyite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &MgWadsleyite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &MgWadsleyite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &MgWadsleyite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &MgWadsleyite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &MgWadsleyite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &MgWadsleyite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &MgWadsleyite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &MgWadsleyite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &MgWadsleyite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &MgWadsleyite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &MgWadsleyite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &MgWadsleyite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &MgWadsleyite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &MgWadsleyite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &MgWadsleyite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &MgWadsleyite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &MgWadsleyite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&MgWadsleyite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&MgWadsleyite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<NaCaFerrite_slb_ph,  std::shared_ptr<NaCaFerrite_slb_ph> >(m, "NaCaFerrite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &NaCaFerrite_slb_ph::identifier,"identifier")
    .def("name", &NaCaFerrite_slb_ph::name,"phase name")
    .def("tcg_build_version", &NaCaFerrite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &NaCaFerrite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &NaCaFerrite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &NaCaFerrite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &NaCaFerrite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &NaCaFerrite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &NaCaFerrite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &NaCaFerrite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &NaCaFerrite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &NaCaFerrite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &NaCaFerrite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &NaCaFerrite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &NaCaFerrite_slb_ph::endmember_number)
    .def("endmember_name", &NaCaFerrite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &NaCaFerrite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &NaCaFerrite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &NaCaFerrite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &NaCaFerrite_slb_ph::species_number)
    .def("species_name", &NaCaFerrite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &NaCaFerrite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &NaCaFerrite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &NaCaFerrite_slb_ph::species_elements, py::arg("i"))
    .def("g", &NaCaFerrite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &NaCaFerrite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &NaCaFerrite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &NaCaFerrite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &NaCaFerrite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &NaCaFerrite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &NaCaFerrite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &NaCaFerrite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &NaCaFerrite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &NaCaFerrite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &NaCaFerrite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &NaCaFerrite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &NaCaFerrite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &NaCaFerrite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &NaCaFerrite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &NaCaFerrite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &NaCaFerrite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &NaCaFerrite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &NaCaFerrite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &NaCaFerrite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &NaCaFerrite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &NaCaFerrite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &NaCaFerrite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &NaCaFerrite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &NaCaFerrite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &NaCaFerrite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &NaCaFerrite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &NaCaFerrite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &NaCaFerrite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &NaCaFerrite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &NaCaFerrite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &NaCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &NaCaFerrite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&NaCaFerrite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &NaCaFerrite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &NaCaFerrite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &NaCaFerrite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &NaCaFerrite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &NaCaFerrite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &NaCaFerrite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &NaCaFerrite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &NaCaFerrite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &NaCaFerrite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &NaCaFerrite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &NaCaFerrite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &NaCaFerrite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &NaCaFerrite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &NaCaFerrite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &NaCaFerrite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &NaCaFerrite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &NaCaFerrite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &NaCaFerrite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &NaCaFerrite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &NaCaFerrite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &NaCaFerrite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &NaCaFerrite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &NaCaFerrite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &NaCaFerrite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &NaCaFerrite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &NaCaFerrite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &NaCaFerrite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&NaCaFerrite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&NaCaFerrite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<NaMajorite_slb_ph,  std::shared_ptr<NaMajorite_slb_ph> >(m, "NaMajorite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &NaMajorite_slb_ph::identifier,"identifier")
    .def("name", &NaMajorite_slb_ph::name,"phase name")
    .def("tcg_build_version", &NaMajorite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &NaMajorite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &NaMajorite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &NaMajorite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &NaMajorite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &NaMajorite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &NaMajorite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &NaMajorite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &NaMajorite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &NaMajorite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &NaMajorite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &NaMajorite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &NaMajorite_slb_ph::endmember_number)
    .def("endmember_name", &NaMajorite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &NaMajorite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &NaMajorite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &NaMajorite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &NaMajorite_slb_ph::species_number)
    .def("species_name", &NaMajorite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &NaMajorite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &NaMajorite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &NaMajorite_slb_ph::species_elements, py::arg("i"))
    .def("g", &NaMajorite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &NaMajorite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &NaMajorite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &NaMajorite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &NaMajorite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &NaMajorite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &NaMajorite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &NaMajorite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &NaMajorite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &NaMajorite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &NaMajorite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &NaMajorite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &NaMajorite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &NaMajorite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &NaMajorite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &NaMajorite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &NaMajorite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &NaMajorite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &NaMajorite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &NaMajorite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &NaMajorite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &NaMajorite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &NaMajorite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &NaMajorite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &NaMajorite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &NaMajorite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &NaMajorite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &NaMajorite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &NaMajorite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &NaMajorite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &NaMajorite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &NaMajorite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &NaMajorite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&NaMajorite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &NaMajorite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &NaMajorite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &NaMajorite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &NaMajorite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &NaMajorite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &NaMajorite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &NaMajorite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &NaMajorite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &NaMajorite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &NaMajorite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &NaMajorite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &NaMajorite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &NaMajorite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &NaMajorite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &NaMajorite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &NaMajorite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &NaMajorite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &NaMajorite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &NaMajorite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &NaMajorite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &NaMajorite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &NaMajorite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &NaMajorite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &NaMajorite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &NaMajorite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &NaMajorite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &NaMajorite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&NaMajorite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&NaMajorite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Nepheline_slb_ph,  std::shared_ptr<Nepheline_slb_ph> >(m, "Nepheline_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Nepheline_slb_ph::identifier,"identifier")
    .def("name", &Nepheline_slb_ph::name,"phase name")
    .def("tcg_build_version", &Nepheline_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Nepheline_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Nepheline_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Nepheline_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Nepheline_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Nepheline_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Nepheline_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Nepheline_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Nepheline_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Nepheline_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Nepheline_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Nepheline_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Nepheline_slb_ph::endmember_number)
    .def("endmember_name", &Nepheline_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Nepheline_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Nepheline_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Nepheline_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Nepheline_slb_ph::species_number)
    .def("species_name", &Nepheline_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Nepheline_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Nepheline_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Nepheline_slb_ph::species_elements, py::arg("i"))
    .def("g", &Nepheline_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Nepheline_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Nepheline_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Nepheline_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Nepheline_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Nepheline_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Nepheline_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Nepheline_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Nepheline_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Nepheline_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Nepheline_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Nepheline_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Nepheline_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Nepheline_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Nepheline_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Nepheline_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Nepheline_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Nepheline_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Nepheline_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Nepheline_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Nepheline_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Nepheline_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Nepheline_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Nepheline_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Nepheline_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Nepheline_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Nepheline_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Nepheline_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Nepheline_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Nepheline_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Nepheline_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Nepheline_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Nepheline_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Nepheline_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Nepheline_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Nepheline_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Nepheline_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Nepheline_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Nepheline_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Nepheline_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Nepheline_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Nepheline_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Nepheline_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Nepheline_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Nepheline_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Nepheline_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Nepheline_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Nepheline_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Nepheline_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Nepheline_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Nepheline_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Nepheline_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Nepheline_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Nepheline_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Nepheline_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Nepheline_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Nepheline_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Nepheline_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Nepheline_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Nepheline_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Nepheline_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Nepheline_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Nepheline_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Olivine_slb_ph,  std::shared_ptr<Olivine_slb_ph> >(m, "Olivine_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Olivine_slb_ph::identifier,"identifier")
    .def("name", &Olivine_slb_ph::name,"phase name")
    .def("tcg_build_version", &Olivine_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Olivine_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Olivine_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Olivine_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Olivine_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Olivine_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Olivine_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Olivine_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Olivine_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Olivine_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Olivine_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Olivine_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Olivine_slb_ph::endmember_number)
    .def("endmember_name", &Olivine_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Olivine_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Olivine_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Olivine_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Olivine_slb_ph::species_number)
    .def("species_name", &Olivine_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Olivine_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Olivine_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Olivine_slb_ph::species_elements, py::arg("i"))
    .def("g", &Olivine_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Olivine_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Olivine_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Olivine_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Olivine_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Olivine_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Olivine_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Olivine_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Olivine_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Olivine_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Olivine_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Olivine_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Olivine_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Olivine_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Olivine_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Olivine_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Olivine_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Olivine_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Olivine_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Olivine_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Olivine_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Olivine_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Olivine_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Olivine_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Olivine_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Olivine_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Olivine_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Olivine_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Olivine_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Olivine_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Olivine_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Olivine_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Olivine_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Olivine_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Olivine_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Olivine_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Olivine_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Olivine_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Olivine_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Olivine_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Olivine_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Olivine_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Olivine_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Olivine_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Olivine_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Olivine_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Olivine_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Olivine_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Olivine_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Olivine_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Olivine_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Olivine_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Olivine_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Olivine_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Olivine_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Olivine_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Olivine_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Olivine_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Olivine_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Olivine_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Olivine_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Olivine_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Olivine_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<OrthoDiopside_slb_ph,  std::shared_ptr<OrthoDiopside_slb_ph> >(m, "OrthoDiopside_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &OrthoDiopside_slb_ph::identifier,"identifier")
    .def("name", &OrthoDiopside_slb_ph::name,"phase name")
    .def("tcg_build_version", &OrthoDiopside_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &OrthoDiopside_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &OrthoDiopside_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &OrthoDiopside_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &OrthoDiopside_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &OrthoDiopside_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &OrthoDiopside_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &OrthoDiopside_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &OrthoDiopside_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &OrthoDiopside_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &OrthoDiopside_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &OrthoDiopside_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &OrthoDiopside_slb_ph::endmember_number)
    .def("endmember_name", &OrthoDiopside_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &OrthoDiopside_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &OrthoDiopside_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &OrthoDiopside_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &OrthoDiopside_slb_ph::species_number)
    .def("species_name", &OrthoDiopside_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &OrthoDiopside_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &OrthoDiopside_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &OrthoDiopside_slb_ph::species_elements, py::arg("i"))
    .def("g", &OrthoDiopside_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &OrthoDiopside_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &OrthoDiopside_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &OrthoDiopside_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &OrthoDiopside_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &OrthoDiopside_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &OrthoDiopside_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &OrthoDiopside_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &OrthoDiopside_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &OrthoDiopside_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &OrthoDiopside_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &OrthoDiopside_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &OrthoDiopside_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &OrthoDiopside_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &OrthoDiopside_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &OrthoDiopside_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &OrthoDiopside_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &OrthoDiopside_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &OrthoDiopside_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &OrthoDiopside_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &OrthoDiopside_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &OrthoDiopside_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &OrthoDiopside_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &OrthoDiopside_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &OrthoDiopside_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &OrthoDiopside_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &OrthoDiopside_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &OrthoDiopside_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &OrthoDiopside_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &OrthoDiopside_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &OrthoDiopside_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &OrthoDiopside_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &OrthoDiopside_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&OrthoDiopside_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &OrthoDiopside_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &OrthoDiopside_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &OrthoDiopside_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &OrthoDiopside_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &OrthoDiopside_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &OrthoDiopside_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &OrthoDiopside_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &OrthoDiopside_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &OrthoDiopside_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &OrthoDiopside_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &OrthoDiopside_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &OrthoDiopside_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &OrthoDiopside_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &OrthoDiopside_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &OrthoDiopside_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &OrthoDiopside_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &OrthoDiopside_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &OrthoDiopside_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &OrthoDiopside_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &OrthoDiopside_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &OrthoDiopside_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &OrthoDiopside_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &OrthoDiopside_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &OrthoDiopside_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &OrthoDiopside_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &OrthoDiopside_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &OrthoDiopside_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&OrthoDiopside_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&OrthoDiopside_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Orthopyroxene_slb_ph,  std::shared_ptr<Orthopyroxene_slb_ph> >(m, "Orthopyroxene_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Orthopyroxene_slb_ph::identifier,"identifier")
    .def("name", &Orthopyroxene_slb_ph::name,"phase name")
    .def("tcg_build_version", &Orthopyroxene_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Orthopyroxene_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Orthopyroxene_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Orthopyroxene_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Orthopyroxene_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Orthopyroxene_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Orthopyroxene_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Orthopyroxene_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Orthopyroxene_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Orthopyroxene_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Orthopyroxene_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Orthopyroxene_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Orthopyroxene_slb_ph::endmember_number)
    .def("endmember_name", &Orthopyroxene_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Orthopyroxene_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Orthopyroxene_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Orthopyroxene_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Orthopyroxene_slb_ph::species_number)
    .def("species_name", &Orthopyroxene_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Orthopyroxene_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Orthopyroxene_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Orthopyroxene_slb_ph::species_elements, py::arg("i"))
    .def("g", &Orthopyroxene_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Orthopyroxene_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Orthopyroxene_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Orthopyroxene_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Orthopyroxene_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Orthopyroxene_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Orthopyroxene_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Orthopyroxene_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Orthopyroxene_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Orthopyroxene_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Orthopyroxene_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Orthopyroxene_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Orthopyroxene_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Orthopyroxene_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Orthopyroxene_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Orthopyroxene_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Orthopyroxene_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Orthopyroxene_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Orthopyroxene_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Orthopyroxene_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Orthopyroxene_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Orthopyroxene_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Orthopyroxene_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Orthopyroxene_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Orthopyroxene_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Orthopyroxene_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Orthopyroxene_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Orthopyroxene_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Orthopyroxene_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Orthopyroxene_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Orthopyroxene_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Orthopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Orthopyroxene_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Orthopyroxene_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Orthopyroxene_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Orthopyroxene_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Orthopyroxene_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Orthopyroxene_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Orthopyroxene_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Orthopyroxene_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Orthopyroxene_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Orthopyroxene_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Orthopyroxene_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Orthopyroxene_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Orthopyroxene_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Orthopyroxene_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Orthopyroxene_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Orthopyroxene_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Orthopyroxene_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Orthopyroxene_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Orthopyroxene_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Orthopyroxene_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Orthopyroxene_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Orthopyroxene_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Orthopyroxene_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Orthopyroxene_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Orthopyroxene_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Orthopyroxene_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Orthopyroxene_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Orthopyroxene_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Orthopyroxene_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Orthopyroxene_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Orthopyroxene_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Periclase_slb_ph,  std::shared_ptr<Periclase_slb_ph> >(m, "Periclase_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Periclase_slb_ph::identifier,"identifier")
    .def("name", &Periclase_slb_ph::name,"phase name")
    .def("tcg_build_version", &Periclase_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Periclase_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Periclase_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Periclase_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Periclase_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Periclase_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Periclase_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Periclase_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Periclase_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Periclase_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Periclase_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Periclase_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Periclase_slb_ph::endmember_number)
    .def("endmember_name", &Periclase_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Periclase_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Periclase_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Periclase_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Periclase_slb_ph::species_number)
    .def("species_name", &Periclase_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Periclase_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Periclase_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Periclase_slb_ph::species_elements, py::arg("i"))
    .def("g", &Periclase_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Periclase_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Periclase_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Periclase_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Periclase_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Periclase_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Periclase_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Periclase_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Periclase_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Periclase_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Periclase_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Periclase_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Periclase_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Periclase_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Periclase_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Periclase_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Periclase_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Periclase_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Periclase_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Periclase_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Periclase_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Periclase_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Periclase_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Periclase_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Periclase_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Periclase_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Periclase_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Periclase_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Periclase_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Periclase_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Periclase_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Periclase_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Periclase_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Periclase_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Periclase_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Periclase_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Periclase_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Periclase_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Periclase_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Periclase_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Periclase_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Periclase_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Periclase_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Periclase_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Periclase_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Periclase_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Periclase_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Periclase_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Periclase_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Periclase_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Periclase_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Periclase_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Periclase_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Periclase_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Periclase_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Periclase_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Periclase_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Periclase_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Periclase_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Periclase_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Periclase_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Periclase_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Periclase_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Perovskite_slb_ph,  std::shared_ptr<Perovskite_slb_ph> >(m, "Perovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Perovskite_slb_ph::identifier,"identifier")
    .def("name", &Perovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Perovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Perovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Perovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Perovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Perovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Perovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Perovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Perovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Perovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Perovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Perovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Perovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Perovskite_slb_ph::endmember_number)
    .def("endmember_name", &Perovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Perovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Perovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Perovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Perovskite_slb_ph::species_number)
    .def("species_name", &Perovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Perovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Perovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Perovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Perovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Perovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Perovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Perovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Perovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Perovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Perovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Perovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Perovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Perovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Perovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Perovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Perovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Perovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Perovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Perovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Perovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Perovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Perovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Perovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Perovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Perovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Perovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Perovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Perovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Perovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Perovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Perovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Perovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Perovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Perovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Perovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Perovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Perovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Perovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Perovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Perovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Perovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Perovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Perovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Perovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Perovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Perovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Perovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Perovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Perovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Perovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Perovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Perovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Perovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Perovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Perovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Perovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Perovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Perovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Perovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Perovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Perovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Perovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Perovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Perovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Perovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Perovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<PostPerovskite_slb_ph,  std::shared_ptr<PostPerovskite_slb_ph> >(m, "PostPerovskite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &PostPerovskite_slb_ph::identifier,"identifier")
    .def("name", &PostPerovskite_slb_ph::name,"phase name")
    .def("tcg_build_version", &PostPerovskite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &PostPerovskite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &PostPerovskite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &PostPerovskite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &PostPerovskite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &PostPerovskite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &PostPerovskite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &PostPerovskite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &PostPerovskite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &PostPerovskite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &PostPerovskite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &PostPerovskite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &PostPerovskite_slb_ph::endmember_number)
    .def("endmember_name", &PostPerovskite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &PostPerovskite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &PostPerovskite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &PostPerovskite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &PostPerovskite_slb_ph::species_number)
    .def("species_name", &PostPerovskite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &PostPerovskite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &PostPerovskite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &PostPerovskite_slb_ph::species_elements, py::arg("i"))
    .def("g", &PostPerovskite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &PostPerovskite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &PostPerovskite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &PostPerovskite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &PostPerovskite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &PostPerovskite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &PostPerovskite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &PostPerovskite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &PostPerovskite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &PostPerovskite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &PostPerovskite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &PostPerovskite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &PostPerovskite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &PostPerovskite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &PostPerovskite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &PostPerovskite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &PostPerovskite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &PostPerovskite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &PostPerovskite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &PostPerovskite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &PostPerovskite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &PostPerovskite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &PostPerovskite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &PostPerovskite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &PostPerovskite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &PostPerovskite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &PostPerovskite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &PostPerovskite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &PostPerovskite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &PostPerovskite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &PostPerovskite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &PostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &PostPerovskite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&PostPerovskite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &PostPerovskite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &PostPerovskite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &PostPerovskite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &PostPerovskite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &PostPerovskite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &PostPerovskite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &PostPerovskite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &PostPerovskite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &PostPerovskite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &PostPerovskite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &PostPerovskite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &PostPerovskite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &PostPerovskite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &PostPerovskite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &PostPerovskite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &PostPerovskite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &PostPerovskite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &PostPerovskite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &PostPerovskite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &PostPerovskite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &PostPerovskite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &PostPerovskite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &PostPerovskite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &PostPerovskite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &PostPerovskite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &PostPerovskite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &PostPerovskite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&PostPerovskite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&PostPerovskite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Pyrope_slb_ph,  std::shared_ptr<Pyrope_slb_ph> >(m, "Pyrope_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Pyrope_slb_ph::identifier,"identifier")
    .def("name", &Pyrope_slb_ph::name,"phase name")
    .def("tcg_build_version", &Pyrope_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Pyrope_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Pyrope_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Pyrope_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Pyrope_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Pyrope_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Pyrope_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Pyrope_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Pyrope_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Pyrope_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Pyrope_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Pyrope_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Pyrope_slb_ph::endmember_number)
    .def("endmember_name", &Pyrope_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Pyrope_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Pyrope_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Pyrope_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Pyrope_slb_ph::species_number)
    .def("species_name", &Pyrope_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Pyrope_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Pyrope_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Pyrope_slb_ph::species_elements, py::arg("i"))
    .def("g", &Pyrope_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Pyrope_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Pyrope_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Pyrope_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Pyrope_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Pyrope_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Pyrope_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Pyrope_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Pyrope_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Pyrope_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Pyrope_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Pyrope_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Pyrope_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Pyrope_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Pyrope_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Pyrope_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Pyrope_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Pyrope_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Pyrope_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Pyrope_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Pyrope_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Pyrope_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Pyrope_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Pyrope_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Pyrope_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Pyrope_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Pyrope_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Pyrope_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Pyrope_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Pyrope_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Pyrope_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Pyrope_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Pyrope_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Pyrope_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Pyrope_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Pyrope_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Pyrope_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Pyrope_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Pyrope_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Pyrope_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Pyrope_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Pyrope_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Pyrope_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Pyrope_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Pyrope_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Pyrope_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Pyrope_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Pyrope_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Pyrope_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Pyrope_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Pyrope_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Pyrope_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Pyrope_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Pyrope_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Pyrope_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Pyrope_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Pyrope_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Pyrope_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Pyrope_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Pyrope_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Pyrope_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Pyrope_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Pyrope_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Quartz_slb_ph,  std::shared_ptr<Quartz_slb_ph> >(m, "Quartz_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Quartz_slb_ph::identifier,"identifier")
    .def("name", &Quartz_slb_ph::name,"phase name")
    .def("tcg_build_version", &Quartz_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Quartz_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Quartz_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Quartz_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Quartz_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Quartz_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Quartz_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Quartz_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Quartz_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Quartz_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Quartz_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Quartz_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Quartz_slb_ph::endmember_number)
    .def("endmember_name", &Quartz_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Quartz_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Quartz_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Quartz_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Quartz_slb_ph::species_number)
    .def("species_name", &Quartz_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Quartz_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Quartz_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Quartz_slb_ph::species_elements, py::arg("i"))
    .def("g", &Quartz_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Quartz_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Quartz_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Quartz_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Quartz_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Quartz_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Quartz_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Quartz_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Quartz_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Quartz_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Quartz_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Quartz_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Quartz_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Quartz_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Quartz_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Quartz_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Quartz_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Quartz_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Quartz_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Quartz_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Quartz_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Quartz_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Quartz_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Quartz_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Quartz_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Quartz_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Quartz_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Quartz_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Quartz_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Quartz_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Quartz_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Quartz_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Quartz_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Quartz_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Quartz_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Quartz_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Quartz_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Quartz_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Quartz_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Quartz_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Quartz_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Quartz_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Quartz_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Quartz_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Quartz_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Quartz_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Quartz_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Quartz_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Quartz_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Quartz_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Quartz_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Quartz_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Quartz_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Quartz_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Quartz_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Quartz_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Quartz_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Quartz_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Quartz_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Quartz_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Quartz_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Quartz_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Quartz_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Ringwoodite_slb_ph,  std::shared_ptr<Ringwoodite_slb_ph> >(m, "Ringwoodite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Ringwoodite_slb_ph::identifier,"identifier")
    .def("name", &Ringwoodite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Ringwoodite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Ringwoodite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Ringwoodite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Ringwoodite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Ringwoodite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Ringwoodite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Ringwoodite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Ringwoodite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Ringwoodite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Ringwoodite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Ringwoodite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Ringwoodite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Ringwoodite_slb_ph::endmember_number)
    .def("endmember_name", &Ringwoodite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Ringwoodite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Ringwoodite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Ringwoodite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Ringwoodite_slb_ph::species_number)
    .def("species_name", &Ringwoodite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Ringwoodite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Ringwoodite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Ringwoodite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Ringwoodite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Ringwoodite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Ringwoodite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Ringwoodite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Ringwoodite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Ringwoodite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Ringwoodite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Ringwoodite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Ringwoodite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Ringwoodite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Ringwoodite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Ringwoodite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Ringwoodite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Ringwoodite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Ringwoodite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Ringwoodite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Ringwoodite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Ringwoodite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Ringwoodite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Ringwoodite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Ringwoodite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Ringwoodite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Ringwoodite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Ringwoodite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Ringwoodite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Ringwoodite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Ringwoodite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Ringwoodite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Ringwoodite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Ringwoodite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Ringwoodite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Ringwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Ringwoodite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Ringwoodite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Ringwoodite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Ringwoodite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Ringwoodite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Ringwoodite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Ringwoodite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Ringwoodite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Ringwoodite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Ringwoodite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Ringwoodite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Ringwoodite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Ringwoodite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Ringwoodite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Ringwoodite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Ringwoodite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Ringwoodite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Ringwoodite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Ringwoodite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Ringwoodite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Ringwoodite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Ringwoodite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Ringwoodite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Ringwoodite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Ringwoodite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Ringwoodite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Ringwoodite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Ringwoodite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Ringwoodite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Ringwoodite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Ringwoodite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Seifertite_slb_ph,  std::shared_ptr<Seifertite_slb_ph> >(m, "Seifertite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Seifertite_slb_ph::identifier,"identifier")
    .def("name", &Seifertite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Seifertite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Seifertite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Seifertite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Seifertite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Seifertite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Seifertite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Seifertite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Seifertite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Seifertite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Seifertite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Seifertite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Seifertite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Seifertite_slb_ph::endmember_number)
    .def("endmember_name", &Seifertite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Seifertite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Seifertite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Seifertite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Seifertite_slb_ph::species_number)
    .def("species_name", &Seifertite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Seifertite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Seifertite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Seifertite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Seifertite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Seifertite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Seifertite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Seifertite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Seifertite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Seifertite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Seifertite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Seifertite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Seifertite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Seifertite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Seifertite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Seifertite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Seifertite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Seifertite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Seifertite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Seifertite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Seifertite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Seifertite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Seifertite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Seifertite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Seifertite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Seifertite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Seifertite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Seifertite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Seifertite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Seifertite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Seifertite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Seifertite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Seifertite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Seifertite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Seifertite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Seifertite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Seifertite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Seifertite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Seifertite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Seifertite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Seifertite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Seifertite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Seifertite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Seifertite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Seifertite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Seifertite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Seifertite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Seifertite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Seifertite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Seifertite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Seifertite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Seifertite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Seifertite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Seifertite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Seifertite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Seifertite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Seifertite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Seifertite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Seifertite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Seifertite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Seifertite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Seifertite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Seifertite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Seifertite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Seifertite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Seifertite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Seifertite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Spinel_slb_ph,  std::shared_ptr<Spinel_slb_ph> >(m, "Spinel_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Spinel_slb_ph::identifier,"identifier")
    .def("name", &Spinel_slb_ph::name,"phase name")
    .def("tcg_build_version", &Spinel_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Spinel_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Spinel_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Spinel_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Spinel_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Spinel_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Spinel_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Spinel_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Spinel_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Spinel_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Spinel_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Spinel_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Spinel_slb_ph::endmember_number)
    .def("endmember_name", &Spinel_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Spinel_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Spinel_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Spinel_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Spinel_slb_ph::species_number)
    .def("species_name", &Spinel_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Spinel_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Spinel_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Spinel_slb_ph::species_elements, py::arg("i"))
    .def("g", &Spinel_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Spinel_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Spinel_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Spinel_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Spinel_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Spinel_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Spinel_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Spinel_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Spinel_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Spinel_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Spinel_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Spinel_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Spinel_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Spinel_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Spinel_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Spinel_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Spinel_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Spinel_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Spinel_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Spinel_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Spinel_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Spinel_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Spinel_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Spinel_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Spinel_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Spinel_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Spinel_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Spinel_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Spinel_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Spinel_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Spinel_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Spinel_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Spinel_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Spinel_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Spinel_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Spinel_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Spinel_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Spinel_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Spinel_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Spinel_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Spinel_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Spinel_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Spinel_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Spinel_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Spinel_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Spinel_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Spinel_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Spinel_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Spinel_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Spinel_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Spinel_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Spinel_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Spinel_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Spinel_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Spinel_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Spinel_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Spinel_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Spinel_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Spinel_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Spinel_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Spinel_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Spinel_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Spinel_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Stishovite_slb_ph,  std::shared_ptr<Stishovite_slb_ph> >(m, "Stishovite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Stishovite_slb_ph::identifier,"identifier")
    .def("name", &Stishovite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Stishovite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Stishovite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Stishovite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Stishovite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Stishovite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Stishovite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Stishovite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Stishovite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Stishovite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Stishovite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Stishovite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Stishovite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Stishovite_slb_ph::endmember_number)
    .def("endmember_name", &Stishovite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Stishovite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Stishovite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Stishovite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Stishovite_slb_ph::species_number)
    .def("species_name", &Stishovite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Stishovite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Stishovite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Stishovite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Stishovite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Stishovite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Stishovite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Stishovite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Stishovite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Stishovite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Stishovite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Stishovite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Stishovite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Stishovite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Stishovite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Stishovite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Stishovite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Stishovite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Stishovite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Stishovite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Stishovite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Stishovite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Stishovite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Stishovite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Stishovite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Stishovite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Stishovite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Stishovite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Stishovite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Stishovite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Stishovite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Stishovite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Stishovite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Stishovite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Stishovite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Stishovite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Stishovite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Stishovite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Stishovite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Stishovite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Stishovite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Stishovite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Stishovite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Stishovite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Stishovite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Stishovite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Stishovite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Stishovite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Stishovite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Stishovite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Stishovite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Stishovite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Stishovite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Stishovite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Stishovite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Stishovite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Stishovite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Stishovite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Stishovite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Stishovite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Stishovite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Stishovite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Stishovite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Stishovite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Stishovite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Stishovite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Stishovite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Wadsleyite_slb_ph,  std::shared_ptr<Wadsleyite_slb_ph> >(m, "Wadsleyite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Wadsleyite_slb_ph::identifier,"identifier")
    .def("name", &Wadsleyite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Wadsleyite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Wadsleyite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Wadsleyite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Wadsleyite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Wadsleyite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Wadsleyite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Wadsleyite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Wadsleyite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Wadsleyite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Wadsleyite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Wadsleyite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Wadsleyite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Wadsleyite_slb_ph::endmember_number)
    .def("endmember_name", &Wadsleyite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Wadsleyite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Wadsleyite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Wadsleyite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Wadsleyite_slb_ph::species_number)
    .def("species_name", &Wadsleyite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Wadsleyite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Wadsleyite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Wadsleyite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Wadsleyite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Wadsleyite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Wadsleyite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Wadsleyite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Wadsleyite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Wadsleyite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Wadsleyite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Wadsleyite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Wadsleyite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Wadsleyite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Wadsleyite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Wadsleyite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Wadsleyite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Wadsleyite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Wadsleyite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Wadsleyite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Wadsleyite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Wadsleyite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Wadsleyite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Wadsleyite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Wadsleyite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Wadsleyite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Wadsleyite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Wadsleyite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Wadsleyite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Wadsleyite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Wadsleyite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Wadsleyite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Wadsleyite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Wadsleyite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Wadsleyite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Wadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Wadsleyite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Wadsleyite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Wadsleyite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Wadsleyite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Wadsleyite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Wadsleyite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Wadsleyite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Wadsleyite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Wadsleyite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Wadsleyite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Wadsleyite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Wadsleyite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Wadsleyite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Wadsleyite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Wadsleyite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Wadsleyite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Wadsleyite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Wadsleyite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Wadsleyite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Wadsleyite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Wadsleyite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Wadsleyite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Wadsleyite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Wadsleyite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Wadsleyite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Wadsleyite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Wadsleyite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Wadsleyite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Wadsleyite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Wadsleyite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Wadsleyite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  py::class_<Wuestite_slb_ph,  std::shared_ptr<Wuestite_slb_ph> >(m, "Wuestite_slb_ph", py::module_local())
    .def(py::init<>())

    // Coder interfaces
    .def("identifier", &Wuestite_slb_ph::identifier,"identifier")
    .def("name", &Wuestite_slb_ph::name,"phase name")
    .def("tcg_build_version", &Wuestite_slb_ph::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &Wuestite_slb_ph::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &Wuestite_slb_ph::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &Wuestite_slb_ph::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("formula", &Wuestite_slb_ph::formula,"formula",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("conv_elm_to_moles", &Wuestite_slb_ph::conv_elm_to_moles,"convert elements to moles",py::arg("e"))
    .def("conv_elm_to_tot_moles", &Wuestite_slb_ph::conv_elm_to_tot_moles,"convert elements to total moles",py::arg("e"))
    .def("conv_elm_to_tot_grams", &Wuestite_slb_ph::conv_elm_to_tot_grams,"convert elements to total grams",py::arg("e"))
    .def("conv_moles_to_elm", &Wuestite_slb_ph::conv_moles_to_elm,"convert moles to elements",py::arg("n"))
    .def("conv_moles_to_tot_moles", &Wuestite_slb_ph::conv_moles_to_tot_moles,"convert moles to total moles",py::arg("n"))
    .def("conv_moles_to_mole_frac", &Wuestite_slb_ph::conv_moles_to_mole_frac,"convert moles to mole fraction",py::arg("n"))
    .def("test_moles", &Wuestite_slb_ph::test_moles,"test moles are valid",py::arg("n"))
    .def("endmember_number", &Wuestite_slb_ph::endmember_number)
    .def("endmember_name", &Wuestite_slb_ph::endmember_name, py::arg("i"))
    .def("endmember_formula", &Wuestite_slb_ph::endmember_formula, py::arg("i"))
    .def("endmember_mw", &Wuestite_slb_ph::endmember_mw, py::arg("i"))
    .def("endmember_elements", &Wuestite_slb_ph::endmember_elements, py::arg("i"))
    .def("species_number", &Wuestite_slb_ph::species_number)
    .def("species_name", &Wuestite_slb_ph::species_name, py::arg("i"))
    .def("species_formula", &Wuestite_slb_ph::species_formula, py::arg("i"))
    .def("species_mw", &Wuestite_slb_ph::species_mw, py::arg("i"))
    .def("species_elements", &Wuestite_slb_ph::species_elements, py::arg("i"))
    .def("g", &Wuestite_slb_ph::g,"total Gibbs Free energy of phase (J)",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdt", &Wuestite_slb_ph::dgdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdp", &Wuestite_slb_ph::dgdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dgdn", &Wuestite_slb_ph::dgdn,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdt2", &Wuestite_slb_ph::d2gdt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdtdp", &Wuestite_slb_ph::d2gdtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdp2", &Wuestite_slb_ph::d2gdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndt", &Wuestite_slb_ph::d2gdndt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdndp", &Wuestite_slb_ph::d2gdndp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d2gdn2", &Wuestite_slb_ph::d2gdn2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt3", &Wuestite_slb_ph::d3gdt3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdp3", &Wuestite_slb_ph::d3gdp3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdt2dp", &Wuestite_slb_ph::d3gdt2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndt2", &Wuestite_slb_ph::d3gdndt2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdtdp2", &Wuestite_slb_ph::d3gdtdp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndtdp", &Wuestite_slb_ph::d3gdndtdp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dt", &Wuestite_slb_ph::d3gdn2dt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdndp2", &Wuestite_slb_ph::d3gdndp2,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn2dp", &Wuestite_slb_ph::d3gdn2dp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("d3gdn3", &Wuestite_slb_ph::d3gdn3,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("v", &Wuestite_slb_ph::v,"Molar Volume",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("s", &Wuestite_slb_ph::s,"Molar entropy",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("alpha", &Wuestite_slb_ph::alpha,"Coefficent of thermal expansion for the phase",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cv", &Wuestite_slb_ph::cv,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("cp", &Wuestite_slb_ph::cp,"Molar Heat capacity of phase at constant pressure",py::arg("T"),py::arg("P"),py::arg("n"))
    .def("dcpdt", &Wuestite_slb_ph::dcpdt,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("beta", &Wuestite_slb_ph::beta,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("K", &Wuestite_slb_ph::K,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("Kp", &Wuestite_slb_ph::Kp,py::arg("T"),py::arg("P"),py::arg("n"))
    .def("get_param_number", &Wuestite_slb_ph::get_param_number,"number of active parameters")
    .def("get_param_names", &Wuestite_slb_ph::get_param_names,"active parameter names")
    .def("get_param_units", &Wuestite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_units", &Wuestite_slb_ph::get_param_units,"active parameter units")
    .def("get_param_values", py::overload_cast<>(&Wuestite_slb_ph::get_param_values),"get active parameter values")
    .def("set_param_values", &Wuestite_slb_ph::set_param_values,"set active parameter values",py::arg("values"))
    .def("get_param_value", &Wuestite_slb_ph::get_param_value,"return value for a particular active parameter",py::arg("index"))
    .def("set_param_value", &Wuestite_slb_ph::set_param_value,"set value for a particular active parameter",py::arg("index"),py::arg("value"))
    .def("dparam_g", &Wuestite_slb_ph::dparam_g,"dparam_g",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdt", &Wuestite_slb_ph::dparam_dgdt,"dparam_dgdt",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdp", &Wuestite_slb_ph::dparam_dgdp,"dparam_dgdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_dgdn", &Wuestite_slb_ph::dparam_dgdn,"dparam_dgdn",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdt2", &Wuestite_slb_ph::dparam_d2gdt2,"dparam_d2gdt2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdtdp", &Wuestite_slb_ph::dparam_d2gdtdp,"dparam_d2gdtdp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d2gdp2", &Wuestite_slb_ph::dparam_d2gdp2,"dparam_d2gdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt3", &Wuestite_slb_ph::dparam_d3gdt3,"dparam_d3gdt3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdt2dp", &Wuestite_slb_ph::dparam_d3gdt2dp,"dparam_d3gdt2dp",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdtdp2", &Wuestite_slb_ph::dparam_d3gdtdp2,"dparam_d3gdtdp2",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))
    .def("dparam_d3gdp3", &Wuestite_slb_ph::dparam_d3gdp3,"dparam_d3gdp3",py::arg("T"),py::arg("P"),py::arg("n"),py::arg("index"))

    // ThermoCodegen interfaces
    .def("endmembers", &Wuestite_slb_ph::endmembers,"vector of pointers to endmembers")
    .def("abbrev", &Wuestite_slb_ph::abbrev,"official phase abbreviation")
    .def("mu", &Wuestite_slb_ph::mu,"chemical potential dG/dn (J/mol)",py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dT", &Wuestite_slb_ph::dmu_dT,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dP", &Wuestite_slb_ph::dmu_dP,py::arg("T"),py::arg("P"),py::arg("x"))
    .def("dmu_dc", &Wuestite_slb_ph::dmu_dc,"change in all chemical potentials for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("ds_dc", &Wuestite_slb_ph::ds_dc,"change in specific entropy for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("dv_dc", &Wuestite_slb_ph::dv_dc,"change in specific volume for a change in composition",
           py::arg("T"),py::arg("P"),py::arg("c"))
    .def("Mass", &Wuestite_slb_ph::Mass,"Molar Mass",py::arg("x"))
    .def("rho", &Wuestite_slb_ph::rho,"density of phase (T,P,C) in weird units ",py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dT", &Wuestite_slb_ph::drho_dT,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dP", &Wuestite_slb_ph::drho_dP,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("drho_dc", &Wuestite_slb_ph::drho_dc,py::arg("T"),py::arg("P"),py::arg("c"))
    .def("c_to_x",  py::overload_cast<std::vector<double>& >  (&Wuestite_slb_ph::c_to_x, py::const_),
             "convert weight fractions to mole fractions",py::arg("c"))
    .def("x_to_c",  py::overload_cast<std::vector<double>& > (&Wuestite_slb_ph::x_to_c, py::const_),
             "convert mole fractions to weight fractions",py::arg("x"));

  

  // Reactions
  py::class_<eclogite_slb_rx,  std::shared_ptr<eclogite_slb_rx>>(m, "eclogite_slb_rx", py::module_local())
    .def(py::init<>())
    .def("name", &eclogite_slb_rx::name)
    .def("tcg_build_version", &eclogite_slb_rx::tcg_build_version, "Version of TCG used to build endmember")
    .def("tcg_build_git_sha", &eclogite_slb_rx::tcg_build_git_sha, "Git SHA of TCG used to build endmember")
    .def("tcg_generation_version", &eclogite_slb_rx::tcg_generation_version, "Version of TCG used to generate endmember")
    .def("tcg_generation_git_sha", &eclogite_slb_rx::tcg_generation_git_sha, "Git SHA of TCG used to generate endmember")
    .def("phases", &eclogite_slb_rx::phases)
    .def("zero_C", &eclogite_slb_rx::zero_C,"returns zero'd compositional 'matrix' of correct shape for phases and endmembers")
    .def("M", &eclogite_slb_rx::M,"compositional 'matrix' M[i][k] of molecular weights of Endmembers in phases")
    .def("nu", &eclogite_slb_rx::nu," 'matrix' of molar stoichiometric coefficients nu[j][i][k]")
    .def("nu_m", &eclogite_slb_rx::nu_m," 'matrix' of mass weighted stoichiometric coefficients nu_m[j][i][k]")
    .def("A", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::A, py::const_),
        "Vector of Affinities for J reactions",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("A", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& , int > (&eclogite_slb_rx::A, py::const_),
        "Affinity for reaction j",py::arg("T"),py::arg("P"),py::arg("C"),py::arg("j"))
    .def("dA_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::dA_dT, py::const_),
        "Vector of derivatives of affinities with respect to Temperature",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("dA_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& , int > (&eclogite_slb_rx::dA_dT, py::const_),
        "Derivative of jth affinity with respect to Temperature",py::arg("T"),py::arg("P"),py::arg("C"),py::arg("j"))
    .def("dA_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::dA_dP, py::const_),
        "Vector of derivatives of affinities with respect to Pressure",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("dA_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& , int > (&eclogite_slb_rx::dA_dP, py::const_),
        "Derivative of jth affinity with respect to Pressure",py::arg("T"),py::arg("P"),py::arg("C"),py::arg("j"))
    .def("dA_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::dA_dC, py::const_),
        "Derivative of all affinities with respect to all components C[i][k]",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("dAj_dCik", &eclogite_slb_rx::dAj_dCik,"derivative of the affinity of reaction J with respect to component k in phase i",
       py::arg("T"),py::arg("P"),py::arg("C"),py::arg("j"),py::arg("i"),py::arg("k"))
    .def("Gamma_i", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::Gamma_i, py::const_),
        "Vector of mass transfer rates for each phase",py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("Gamma_i", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int > (&eclogite_slb_rx::Gamma_i, py::const_),
        "Mass transfer rate for phase i",py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"))
    .def("dGamma_i_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_i_dT, py::const_),
        "Derivative of mass transfer rate for all phases with respect to Temperature",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("dGamma_i_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int > (&eclogite_slb_rx::dGamma_i_dT, py::const_),
        "Derivative of mass transfer rate for phase i with respect to Temperature",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"))
    .def("dGamma_i_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_i_dP, py::const_),
        "Derivative of mass transfer rate for all phases with respect to Pressure",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("dGamma_i_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int > (&eclogite_slb_rx::dGamma_i_dP, py::const_),
        "Derivative of mass transfer rate for phase i with respect to Pressure",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"))
    .def("dGamma_i_dPhi", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_i_dPhi, py::const_),
        "Derivative of mass transfer rate for for all phases with respect to phase fraction",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi") )
    .def("dGamma_i_dPhi", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, unsigned int, unsigned int > (&eclogite_slb_rx::dGamma_i_dPhi, py::const_),
        "Derivative of mass transfer rate for phase i with respect to phase fraction l",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"),py::arg("l"))
    .def("dGamma_i_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_i_dC, py::const_),
        "Derivative of phase mass transfer rate for a change in composition ",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("dGamma_i_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, unsigned int, unsigned int, unsigned int > (&eclogite_slb_rx::dGamma_i_dC, py::const_),
        "Derivative of mass transfer rate for phase i with respect to component k in phase l",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"),py::arg("l"),py::arg("k"))
    .def("Gamma_ik", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&  > (&eclogite_slb_rx::Gamma_ik, py::const_),
        "Vector of vectors of mass transfer rates",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("Gamma_ik", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int, int > (&eclogite_slb_rx::Gamma_ik, py::const_),
        "Mass transfer rate for component k in phase i",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"),py::arg("k"))
    .def("dGamma_ik_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_ik_dT, py::const_),
        "Derivative of mass transfer rate for every component with respect to Temperature",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi") )
    .def("dGamma_ik_dT", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int, int > (&eclogite_slb_rx::dGamma_ik_dT, py::const_),
        "Derivative of mass transfer rate for component k in phase i with respect to Temperature",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"),py::arg("k"))
    .def("dGamma_ik_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>& > (&eclogite_slb_rx::dGamma_ik_dP, py::const_),
        "Derivative of mass transfer rate for component k in phase i with respect to Pressure",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"))
    .def("dGamma_ik_dP", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int, int > (&eclogite_slb_rx::dGamma_ik_dP, py::const_),
        "Derivative of mass transfer rate for component k in phase i with respect to Pressure",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"),py::arg("k"))
    .def("dGamma_ik_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int  > (&eclogite_slb_rx::dGamma_ik_dC, py::const_),
        "Derivative of mass transfer rate component i,k with respect to composition of phase i",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"))
    .def("dGamma_ik_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, unsigned, unsigned, unsigned, unsigned > (&eclogite_slb_rx::dGamma_ik_dC, py::const_),
        "Derivative of mass transfer rate for component i,k  with respect to component l,m",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),
        py::arg("i"),py::arg("k"),py::arg("l"),py::arg("m"))
    .def("dGamma_ik_dPhi", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, int  > (&eclogite_slb_rx::dGamma_ik_dPhi, py::const_),
        "Derivative of mass transfer rate for all components i phase i  with respect to phase fraction Phi",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),py::arg("i"))
    .def("dGamma_ik_dPhi", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& ,
        std::vector<double>&, unsigned, unsigned, unsigned  > (&eclogite_slb_rx::dGamma_ik_dPhi, py::const_),
        "Derivative of mass transfer rate for component i,k with respect to phase fraction l",
        py::arg("T"),py::arg("P"),py::arg("C"),py::arg("Phi"),
        py::arg("i"),py::arg("k"),py::arg("l"))
    .def("rho", py::overload_cast<const double& , const double& , 
        std::vector<std::vector<double> >& > (&eclogite_slb_rx::rho, py::const_),
        "density of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("Cp", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::Cp, py::const_),
        "Vector of heat capacities of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("s", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::s, py::const_),
        "Vector of entropies of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("ds_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::ds_dC, py::const_),
        "Vector of ds_dc for each phase of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("alpha", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::alpha, py::const_),
        "Vector of thermal expansivity of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("beta", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::beta, py::const_),
        "Vector of compressibility of phases",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("Mu", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::Mu, py::const_),
        "Matrix of chemical potentials mu_i^k(T,P,C)",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("dMu_dC", py::overload_cast<const double& , const double& , std::vector<std::vector<double> >& > (&eclogite_slb_rx::dMu_dC, py::const_),
        "Compositional Jacobian  dmu_i^k_dC(T,P,C)",py::arg("T"),py::arg("P"),py::arg("C"))
    .def("C_to_X",  py::overload_cast<std::vector<std::vector<double> >& > (&eclogite_slb_rx::C_to_X, py::const_),
        "convert all weight fractions to mole fractions",py::arg("C"))
    .def("X_to_C",  py::overload_cast<std::vector<std::vector<double> >& > (&eclogite_slb_rx::X_to_C, py::const_),
        "convert all mole fractions to weight fractions",py::arg("X"))
    .def("set_parameter", &eclogite_slb_rx::set_parameter,"Set value of parameter",py::arg("p"),py::arg("val"))
    .def("get_parameter", &eclogite_slb_rx::get_parameter,"Get value of parameter",py::arg("p"),
         py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
    .def("list_parameters", &eclogite_slb_rx::list_parameters,"List available parameters",
         py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
    .def("report", &eclogite_slb_rx::report,"dump phases and end members",
         py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>());

  
}
