import sympy as sym
import re
import types
from thermocodegen.coder import coder

class TCGEndmember:
    name    = None
    formula = None
    model   = None
    param_strs = None
    param_unit = None
    param_vals = None
    A = None # Helmholtz
    G = None # Gibbs
    syms = None
    param_syms = None
    model_type = 'TP' # or 'TV for Helmholtz/Stixrude'
    
    def __init__(self,name,formula,reference,model_type):
        self.name = name
        self.formula = formula
        self.reference = reference
        self.model_type = model_type

        self.model = coder.StdStateModel.from_type(self.model_type)

        self.T = self.model.get_symbol_for_t()
        self.V = self.model.get_symbol_for_v()
 
        self.T_r = self.model.get_symbol_for_tr()
        self.V_r = self.model.get_symbol_for_vr()
        self.P = self.model.get_symbol_for_p()
        self.P_r = self.model.get_symbol_for_pr() 
    
    def r_Fe(self):
        """Return number of Fe atoms in a formula."""
        formula_split = re.findall('[A-Z][^A-Z]*', self.formula)
        elem_count = {}
        for elem in formula_split:
            elem_split = re.split('(\d+)', elem)
            elem = elem_split[0].rstrip('(')
            elem_num = 1 if len(elem_split) == 1 else int(elem_split[1])
            elem_count[elem] = elem_count.get(elem,0) + elem_num
        return elem_count.get('Fe',0)  

    def add_parameter(self, name, unit, symbol, value):
        self.param_strs.append(name)
        self.param_unit.append(unit)
        self.param_syms.append(symbol)
        self.param_vals['name'] = value
    
    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        values_dict = dict(
                           name=self.name,
                           formula=self.formula,
                           reference=self.reference
                          )
        values_dict.update(self.param_vals)
        return values_dict

    def set_model_values(self):
        values_dict = self.model.get_values()
        values_dict.update(self.values_dict())
        print(values_dict)
        self.model.set_values(values_dict)
    
    def tofile(self,outdir):
        self.add_potential_to_model()
        self.set_model_values()
        filename = self.model.to_xml(path=outdir)
        return filename


class SLBEndmember(TCGEndmember):

    def __init__(self,name,formula,reference,**kwargs):

        super().__init__(name,formula,reference,'TV')

        # Mitchell - changed some units to match SLB table
        # a0 should be J/mol instead of J/m
        # R should be J/K-mol instead of J/K-m
        #

        self.param_strs = [ 'a0', 'n', 'v0', \
                            'k00', 'k0p', 'theta0', \
                            'gamma0', 'q', 'r_Fe', \
                            'R' ]
        self.param_unit = [ 'J/mol', '', 'cm^3/mol/10',\
                            'bar', '', 'K', \
                            '', '', '', \
                            'J/K/mol' ]

        symdict = dict( ( p, sym.symbols(p,real=True) ) for p in self.param_strs )
        self.syms = types.SimpleNamespace(**symdict)

        self.param_syms = [ symdict[p] for p in self.param_strs ]
        
        required_params = self.param_strs+['T_r','V_r']
        required_params.remove('R')    # R isn't actually required
        required_params.remove('r_Fe') # neither is r_Fe
        
        missing_params = [ p for p in required_params if p not in kwargs ]
        if len(missing_params) > 0:
            raise Exception("Not all parameter values provided.  Missing: "+", ".join(missing_params))
        
        self.param_vals = dict( (p, kwargs[p]) for p in required_params )
        self.param_vals['r_Fe'] = self.r_Fe()
        self.param_vals['R'] = kwargs.get('R', 8.31446261815324)

    def A_iso_default(self):
        """Isotropic strain contributions."""
        f = ((self.syms.v0/self.V)**(sym.S(2)/sym.S(3)) - 1)/2
        c1 = sym.S(9)/sym.S(2)*self.syms.k00*self.syms.v0
        c2 = self.syms.k0p-4
        A_iso = self.syms.a0 + c1 * f**2 * (1 + c2*f)
        return A_iso
    
    def A_quasi_default(self):
        db = coder.get_external_functions()['Debye']
        f = ((self.syms.v0/self.V)**(sym.S(2)/sym.S(3)) - 1)/2
        c3 = 6 * self.syms.gamma0
        c4 = c3 * (-2 + c3 - 3*self.syms.q)/2
        d0 = self.syms.theta0 * sym.sqrt( 1 + c3*f + c4 * f**2) 
        x = d0/self.T
        A_db = self.syms.n*self.syms.R*self.T*(3*sym.log(1-sym.exp(-x)) - db(x))
        A_quasi = A_db - A_db.subs(self.T,self.T_r)
        return A_quasi
    
    def A_mag_default(self):

        # FIXME: For Cian - Is it correct to mutliply by R here?
        # SLB 2011, Eq 27: 
        # Math: \mathcal{F}_{mag} = -T  \sum r_i \ln(2S_i+1)

        A_mag = -self.syms.R*self.T*self.syms.r_Fe*sym.log(5)
        return A_mag
    
    def A_default(self):
        return self.A_iso_default() + self.A_quasi_default() + self.A_mag_default()

    def add_potential_to_model(self):
        if self.A is None: self.A = self.A_default()
        self.model.add_potential_to_model('A', self.A, self.params())

class BermanEndmember(TCGEndmember):
    
    def __init__(self,name,formula,reference,**kwargs):
        super().__init__(name,formula,reference,'TP')

        self.param_strs = [ 'H_TrPr', 'S_TrPr', 'V_TrPr', \
                            'k0', 'k1', 'k2', 'k3',\
                            'v1', 'v2', 'v3', 'v4']

        self.param_unit = [ 'J', 'J/K', 'J/bar-m',\
                            'J/K-m', 'J/K^(1/2)-m', 'J-K/m', 'J-K^2',\
                            '1/bar', '1/bar^2', '1/K', '1/K^2']

        symdict = dict( ( p, sym.symbols(p,real=True) ) for p in self.param_strs )
        self.syms = types.SimpleNamespace(**symdict)

        self.param_syms = [ symdict[p] for p in self.param_strs ]
        
        required_params = self.param_strs
        
        missing_params = [ p for p in required_params if p not in kwargs ]
        if len(missing_params) > 0:
            raise Exception("Not all parameter values provided.  Missing: "+", ".join(missing_params))
        
        self.param_vals = dict( (p, kwargs[p]) for p in required_params )
        self.param_vals['r_Fe'] = self.r_Fe()
    
    def G_Pr_default(self):        
        # Expression for heat capacity
        Cp_Pr = self.syms.k0+self.syms.k1/sym.sqrt(self.T)+self.syms.k2/self.T**2+self.syms.k3/self.T**3
        return self.syms.H_TrPr + sym.integrate(Cp_Pr,(self.T,self.T_r,self.T)) - self.T*(self.syms.S_TrPr + sym.integrate(Cp_Pr/self.T,(self.T,self.T_r,self.T)))
    
    def G_PrToP_default(self):
        return sym.integrate(self.syms.V_TrPr*(1+self.syms.v1*(self.P-self.P_r)+self.syms.v2*(self.P-self.P_r)**2+self.syms.v3*(self.T-self.T_r)+self.syms.v4*(self.T-self.T_r)**2),(self.P,self.P_r,self.P))

    def G_default(self):
        return self.G_Pr_default() + self.G_PrToP_default()
        
    def add_potential_to_model(self):
        if self.G is None: self.G = self.G_default()
        self.model.add_potential_to_model('G', self.G, self.params())

class BermanPolyEndmember(TCGEndmember):
    
    # instead of a row of a datafram as a dict (**kwargs),
    # accept a dataframe consisting of multiple rows
    def __init__(self,name,formula,reference,polymorphs,df):

        # init the SS model
        super().__init__(name,formula,reference,'TP')      

        # set up the root symbols like we would for any EM
        root_strs = [
            'H_TrPr', 'S_TrPr', 'V_TrPr',
            'k0', 'k1', 'k2', 'k3',
            'v1', 'v2', 'v3', 'v4'
        ]
        root_unit = [
            'J', 'J/K', 'J/bar-m',
            'J/K-m', 'J/K^(1/2)-m', 'J-K/m', 'J-K^2',
            '1/bar', '1/bar^2', '1/K', '1/K^2'
        ]
        # convenient way to get a dict
        params_root = coder.set_coder_params(root_strs, root_unit)
        symbol_dict_root = coder.get_symbol_dict_from_params(params_root)

        # make params for each polymorph according to abbreviations
        abbrevs = df['sAbbrev'].tolist()
        param_strs = []
        param_unit = []
        for abbrev in abbrevs:
            # for each root string, concatenate the EM's abbrevation
            param_strs += ['{}_{}'.format(p,abbrev) for p in root_strs]
            param_unit += root_unit

        # we set our param strings/units as the polymorphic ones
        self.param_strs = param_strs
        self.param_unit = param_unit

        # get a dict
        params = coder.set_coder_params(self.param_strs, self.param_unit)
        symbol_dict = coder.get_symbol_dict_from_params(params)

        # get a list of symbols 
        self.param_syms = [ symbol_dict[p] for p in self.param_strs ]

        # combine the root symbols and the polymorphic symbols
        self.syms = types.SimpleNamespace(**symbol_dict_root, **symbol_dict)

        # override P,T, etc. with updated coder model
        # not sure if this is necessary
        self.model = coder.StdStateModel.from_type()
        self.T = self.model.get_symbol_for_t()
        self.P = self.model.get_symbol_for_p()
        self.T_r = self.model.get_symbol_for_tr()
        self.P_r = self.model.get_symbol_for_pr()

        # create an expression for volume
        self.V = self.syms.V_TrPr*(1+self.syms.v1*(self.P-self.P_r)+self.syms.v2*(self.P-self.P_r)**2+self.syms.v3*(self.T-self.T_r)+self.syms.v4*(self.T-self.T_r)**2)
        
        # override self.G for solid solution
        GPrToP = sym.integrate(self.V,(self.P,self.P_r,self.P))
        GPr = self.G_Pr_default()
        G_ss = GPr + GPrToP
        
        N = len(root_strs)

        # list of Gibbs symbols for each polymorph
        G_p = [sym.Symbol('G_{}'.format(a)) for a in abbrevs]
        symbol_list = list(symbol_dict.values())

        for i,p in enumerate(polymorphs):
            subs_dict = dict(zip(symbol_dict_root.values(),symbol_list[i*N:(i+1)*N]))
            # add the subscripts to each G_p
            G_p[i]=G_ss.subs(subs_dict)
        
        # minimize the list of Gibbs of each polymorph
        self.G = sym.Min(*G_p)
        
        # get the values from the dataframe
        param_vals = dict()
        for ind, row in df.iterrows():
            abbrev = row.pop('sAbbrev')
            new_columns = ['{}_{}'.format(p,abbrev) for p in root_strs]
            row.rename(index=dict(zip(root_strs,new_columns)),inplace=True)
            param_vals.update(row.to_dict())

        self.param_vals = param_vals
        self.param_vals['r_Fe'] = self.r_Fe()
        self.param_vals['R'] = 8.31446261815324

    def G_Pr_default(self):   
        # Heat Capacity
        Cp_Pr = self.syms.k0+self.syms.k1/sym.sqrt(self.T)+self.syms.k2/self.T**2+self.syms.k3/self.T**3
        return self.syms.H_TrPr + sym.integrate(Cp_Pr,(self.T,self.T_r,self.T)) - self.T*(self.syms.S_TrPr + sym.integrate(Cp_Pr/self.T,(self.T,self.T_r,self.T)))
        
    def add_potential_to_model(self):
        assert(self.G is not None)
        self.model.add_potential_to_model('G', self.G, self.params())