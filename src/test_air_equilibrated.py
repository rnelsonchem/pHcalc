from pHcalc import Inert as IonAq
from pHcalc import Acid as AcidAq
from pHcalc import System
from pHcalc.air_equilibrated import SystemAirEquilibrated


# Test code for AirEquilibrated module

# 1. Test transfer of kwargs from sub to super
sys0 = System()
syst = SystemAirEquilibrated()
assert sys0.Kw == syst.Kw

testKw=2e-14
syst = SystemAirEquilibrated(Kw=testKw)
assert syst.Kw == testKw


print()


# 2. Simple exazmple

sys0 = System()
sys0.pHsolve()
print('pure water, closed system                     : pH = ', sys0.pH)

syst = SystemAirEquilibrated()
syst.DICsolve()
print('pure water, air-equilibrated system (in 2021) : pH = ', syst.pH)

syst1972 = SystemAirEquilibrated(P_CO2 = 372.46e-6)
syst1972.DICsolve()
print('pure water, air-equilibrated system (in 1972) : pH = ', syst1972.pH)

print()


# 3. Benchmark/example calculations

# As benchmark calculations we used the on-line version of aqion which uses
# PHREEQC as a back-end for calculations.
#   https://www.aqion.onl/
#   https://www.usgs.gov/software/phreeqc-version-3
# The on-line results are identical to those obtained with the Windows version
# of aqion.
# Reference temperature: 298 K
#
# At higher ionic strengths, deviations appear between the values calculated by
# pHcalc and those from aqion/PHREEQC. This is due to the fact that PHREEQC
# uses activities and ionic strength-dependent equilibrium constants.
# pHcalc uses simple mass-action law with constant equilibrium constants and
# actual concentrations


# CO2 value used by aqion (watch out: aqion uses pCO2, not partial pressure)
aqionP_CO2 = 10**-3.408 # conversion of pCO2 into P_CO2 [atm]

# Test cases
paramlist = [
    {'descr':       '0.1 M HCl(aq) in eq. with amosphere',
     'ioncharge':  -1,
     'ionconc':     0.1,
     'aqion_pH':    1.08},
    {'descr':       '1 mM HCl(aq) in eq. with amosphere',
     'ioncharge':  -1,
     'ionconc':     0.001,
     'aqion_pH':    3.02},
    {'descr':       'pure water in eq. with amosphere',
     'ioncharge':   None,
     'ionconc':     None,
     'aqion_pH':    5.61},
    {'descr':       '1 mM NaOH(aq) in eq. with amosphere',
     'ioncharge':   1,
     'ionconc':     0.001,
     'aqion_pH':    8.20},
    {'descr':       '10 mM NaOH(aq) in eq. with amosphere',
     'ioncharge':   1,
     'ionconc':     0.01,
     'aqion_pH':    9.11},
    {'descr':       '0.1 M NaOH(aq) in eq. with amosphere',
     'ioncharge':   1,
     'ionconc':     0.1,
     'aqion_pH':    9.71},
    ]

for param in paramlist:
    if param['ionconc'] is None:
        syst = SystemAirEquilibrated(P_CO2 = aqionP_CO2)
    else:
        ion = IonAq(charge = param['ioncharge'],
                    conc = param['ionconc'])
        syst = SystemAirEquilibrated(ion, P_CO2 = aqionP_CO2)
    syst.DICsolve()
    print('{0:40s} | DIC = {1:7.3e} M | pH = {2:5.2f} [aqion: {3:5.2f}]'.\
          format(param['descr'], syst.DIC, syst.pH, param['aqion_pH']))
