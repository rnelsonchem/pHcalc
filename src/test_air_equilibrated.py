from pHcalc import Inert as IonAq
from pHcalc import Acid as AcidAq
from pHcalc import System
from pHcalc.air_equilibrated import SystemAirEquilibrated

# Test code for AirEquilibrated

# 1. test transfer of kwargs from sub to super
sys0 = System()
syst = SystemAirEquilibrated()
assert sys0.Kw == syst.Kw

testKw=2e-14
syst = SystemAirEquilibrated(Kw=testKw)
assert syst.Kw == testKw


# 2. benchmark calculations

# TO DO make the following a loop over
#       a dictionary or over list of dictionaries, 
#       include benchmark data from Aqion/Phreeqc 
#
hcl = IonAq(charge=-1, conc=0.1)
syst = SystemAirEquilibrated(hcl)
syst.DICsolve()
print('0.1M HCl in equilibrium with atmosphere pH = ', syst.pH,
      ' DIC = ', syst.DIC, 'M')
print()    


hcl = IonAq(charge=-1, conc=0.001)
syst = SystemAirEquilibrated(hcl)
syst.DICsolve()
print('1 mM HCl in equilibrium with atmosphere pH = ', syst.pH,
      ' DIC = ', syst.DIC, 'M')
print()    


pure_aireq = SystemAirEquilibrated()
pure_aireq.DICsolve()
print('pure water in equilibrium with atmosphere pH = ', pure_aireq.pH,
      ' DIC = ', pure_aireq.DIC, 'M')
print()

naoh = IonAq(charge=1, conc=0.001)
syst = SystemAirEquilibrated(naoh)
syst.DICsolve()
print('1 mM NaOH in equilibrium with atmosphere pH = ', syst.pH,
      ' DIC = ', syst.DIC, 'M')
print()

naoh = IonAq(charge=1, conc=0.01)
syst = SystemAirEquilibrated(naoh)
syst.DICsolve()
print('0.01M NaOH in equilibrium with atmosphere pH = ', syst.pH,
      ' DIC = ', syst.DIC, 'M')
print()

naoh = IonAq(charge=1, conc=0.1)
syst = SystemAirEquilibrated(naoh)
syst.DICsolve()
print('0.1M NaOH in equilibrium with atmosphere pH = ', syst.pH,
      ' DIC = ', syst.DIC, 'M')
print()


pure_aireq = SystemAirEquilibrated(P_CO2 = 320e-6)
pure_aireq.DICsolve()
print('pure water in equilibrium with atmosphere (1972) pH = ', pure_aireq.pH,
      ' DIC = ', pure_aireq.DIC, 'M')
print()