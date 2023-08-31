# %%
s_per_yr = 3.156e+7
T = 600+273.

r0_gcm2yr = 1e-2
r0 = r0_gcm2yr / 1e3 * 1e2 * 1e2 / s_per_yr # kg/m^2/s


t0 = 45e6 * s_per_yr # s
dg_mm = 3
dg = dg_mm / 1e3 # m


rho0 = 3000. # kg/m^3
Gamma0 = r0/dg # kg/m^3/s
Da = t0 * Gamma0 / rho0 
print(Da)

# Da ranges b/w 0.5 and 500,000
# %%
