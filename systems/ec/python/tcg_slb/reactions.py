from thermocodegen.codegen import reaction
import sympy as sym

class SLBReactions:
    
    def __init__(self, name, nreactions, phases, database, reference, **kwargs):
        self.name = name
        self.J = nreactions
        self.phases = phases
        self.reference = reference
        R = kwargs.get('R', 8.31446261815324)
        T0 = kwargs.get('T0', 2000.0)
        
        self.rxn = reaction.Reaction.from_database(total_reactions=nreactions,phase_names=phases,database=database)
        self.j = 0
        
        self.T   = self.rxn.T
        self.P   = self.rxn.P
        self.C   = self.rxn.C
        self.Phi = self.rxn.Phi
        self.A   = self.rxn.A
        self.T0  = sym.symbols('T0', real=True)
        self.R   = sym.symbols('R', real=True)
        
        self.param_strs = [   'T0',       'R']
        self.param_unit = [    'K', 'J/mol/K']
        self.param_syms = [self.T0,    self.R]
        self.param_vals = { 'T0' : T0,
                            'R'  : R }
        
    
    def arrheniusrate(self):
        return sym.exp(-self.T0/self.T)

    def availability(self,reactants):
        S0 = 1
        for phase,endmember in reactants:
            i  =  self.phases.index(phase)
            S0 = S0*self.Phi[i]
        return S0
    
    def reactantlist(self,reactants):
        reactantl = []
        for phase,endmember in reactants:
            reactantl.append([endmember[:4]+"_"+phase[:4],\
                              phase,\
                              endmember])
        return reactantl
    
    def add_reaction(self,reactants,products):
        rp = self.arrheniusrate()
        rm = self.arrheniusrate()
        
        S0p = self.availability(reactants)
        S0m = self.availability(products)
        
        reaction = sym.Piecewise((rp*S0p*self.A[self.j]/self.R/self.T,self.A[self.j]>=0),\
                                 (rm*S0m*self.A[self.j]/self.R/self.T,self.A[self.j]<0),\
                                 (0,True))
        reactantl = self.reactantlist(reactants)
        productl  = self.reactantlist(products)
        rname = ''.join([em[:4] for ph,em in reactants])+\
                '_'+''.join([em[:4] for ph,em in products])
        self.rxn.add_reaction_to_model(rname, reactantl, productl, reaction, self.params())
        self.j += 1

    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        values_dict = dict(name=self.name,
                           reference=self.reference)
        values_dict.update(self.param_vals)
        return values_dict
    
    def tofile(self,outdir):
        values_dict = self.rxn.get_values()
        values_dict.update(self.values_dict())
        self.rxn.set_values(values_dict)
        self.rxn.to_xml(path=outdir)
