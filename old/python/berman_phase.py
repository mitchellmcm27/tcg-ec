import numpy
import sympy as sym
import pandas as pd
from thermocodegen.coder import coder

class BermanPhase:
    name   = None
    abbrev = None
    reference = None

    c = None
    X = None
    n = None
    nT = None
    
    endmember_names = None
    formula_str     = None
    conversion_strs = None
    test_strs       = None
    
    dsyms = None
    Wsarr = None
    m = None
    C = None
    
    Rsym = None
    Trsym = None
    G = None
    
    def __init__(self,name,abbrev,reference,endmember_names,\
                      formula_str,conversion_strs,\
                      **kwargs):

        R     = kwargs.get('R', 8.31446261815324)
        T_r   = kwargs.get('T_r', 300.0)
        
        self.name   = name
        self.abbrev = abbrev
        self.reference = reference
        
        # endmembers (c)
        self.endmember_names = endmember_names
        self.formula_str     = formula_str
        self.conversion_strs = conversion_strs

        # number of components
        self.c = len(endmember_names)
        
        self.test_strs = ['['+repr(i)+'] > 0.0' for i in  range(self.c)]
        
        self.Rsym  = sym.symbols('R', real=True)
        self.Trsym = sym.symbols('T_r', real=True)

        self.param_strs = ['R', 'T_r']
        self.param_unit = ['J/K/mol', 'K']
        self.param_syms = [self.Rsym, self.Trsym]

        self.param_vals = {'R' : R, 'T_r' : T_r}

        self.model = coder.SimpleSolnModel.from_type(nc=self.c)
        
        self.n  = self.model.n
        self.nT = self.model.nT
        self.X  = self.n/self.nT
        
        self.T  = self.model.get_symbol_for_t()
        self.P  = self.model.get_symbol_for_p()
        self.mu = self.model.mu
    
    def G_ss_default(self):
        G_ss = (self.n.transpose()*self.mu)[0]
        return G_ss

    def G_config_default(self):
        S_config = sym.symbols('S_config')
        S_config = 0
        for i in range(0,self.c):
            S_config += self.X[i]*sym.log(self.X[i])
        S_config *= -self.Rsym * self.nT
        G_config = sym.simplify(-self.T*S_config)
        return G_config

    def G_default(self):
        return self.G_ss_default() + self.G_config_default()

    def add_parameter(self, name, unit, symbol, value):
        self.param_strs.append(name)
        self.param_unit.append(unit)
        self.param_syms.append(symbol)
        self.param_vals[name] = value
    
    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        values_dict = dict(
                           formula_string=self.formula_str,
                           conversion_string=self.conversion_strs,
                           test_string=self.test_strs,
                           name=self.name,
                           abbrev=self.abbrev,
                           reference=self.reference,
                           endmembers = self.endmember_names
                          )
        values_dict.update(self.param_vals)
        return values_dict
    
    def add_potential_to_model(self):
        self.model.add_potential_to_model('G', self.G, self.params())

    def set_model_values(self):
        values_dict = self.model.get_values()
        values_dict.update(self.values_dict())
        self.model.set_values(values_dict)
    
    def tofile(self,outdir):
        if self.G is None: self.G = self.G_default()
        self.add_potential_to_model()
        self.set_model_values()
        filename = self.model.to_xml(path=outdir)
        return filename