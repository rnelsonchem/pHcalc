from .pHcalc import Acid as AcidAq
from .pHcalc import Inert as IonAq
from .pHcalc import System
from .secant import secant

class SystemAirEquilibrated(System):
    """Class representing a system at equilibrium with atmospheric CO2
    
    This is a subclass of System.
    
    It automatically adds carbonic acid as a AcidAq species. Do not supply
    this separately. If you need to add additional carbonate and bicarbonate
    ions to your air-equilibrated system, you should only supply the 
    counterions (e.g. sodium if you work with sodium bicarbonate)
    
    The atmospheric CO2 (partial pressure P_CO2) is held in equilibrium with
    the dissolved CO2(aq)+H2CO3(aq). The combination of these latter two, which
    are in equilibrium, is considered as a single species 'H2CO3x'. The H2CO3x
    then further reacts to form bicarbonate HCO3- and carbonate CO32- according
    to the aqueous acid-base chemistry of carbonic acid. The sum of the 
    concentrations [CO32-] + [HCO3-] + [H2CO3x] is the 'dissolved inorganic
    carbon' DIC. 
    
 
    Parameters
    ----------
    *species, Kw
        see System class
    
    H_s : float (optional, default 0.03429)
        Henry solubility [M atm-1] of CO2 in water. The default value is for
        pure water at 298K. Taken from ref. [2], Table II, converted from 
        molality to molarity, i.e. 0.997 * 10**-1.463
        
    P_CO2 : float (optional, default 415e-6)
        Partial pressure [atm] of carbon dioxide in the atmosphere. The
        default value is the global average atmospheric carbon dioxide in 2021
        as reported by NOAA’s Global Monitoring Lab. This value is rising over
        the years.
        https://www.climate.gov/news-features/understanding-climate/climate-change-atmospheric-carbon-dioxide
    
    carbonic_pKa : list of two floats (optional, default = [6.35, 10.33])
        The two pKa values of the carbonic acid system. The default values
        are for water at 298K. Refs [1] and [2].
    

    Attributes
    ---------
    DIC : float
        Total dissolved inorganic carbon (DIC) concentration [M], i.e. the
        sum of [H2CO3x] + [HCO3-] + [CO32-]
    
    Further attributes
        see System class
    
    
    Background
    ----------
    For now, all constants refer to STP (T=298K, P=1atm). 
    
    References:

    [1] H. S. Harned, S. R. Scholes Jr. 
        "The Ionization Constant of HC03- from 0 to 50°." 
        J. Am. Chem. Soc. 1941, 63,1706

    [2] H. S. Harned, R. Davis.
        "The Ionization Constant of Carbonic Acid in Water and the Solubility
        of Carbon Dioxide in Water and Aqueous Salt Solutions from 0 to 50°."
        J. Am. Chem. Soc. 1943, 65 , 2030
      
    """
    def __init__(self, *species,
                 H_s = 0.03429,
                 P_CO2 = 415e-6,
                 carbonic_pKa = [6.35, 10.33], 
                 **kwargs):
        
        # Initialize system with supplied species
        super().__init__(*species, **kwargs)

        # Add atmospheric CO2 reservoir
        #
        # (a) Henry's law:
        #     H_s = c_H2CO3x/P_CO2
        self.c_H2CO3x_goal = H_s * P_CO2
        # (b) The amount of DIC (dissolved inorganic carbon) should be adjusted
        #     such that  
        #     c_H2CO3x_calc equals c_H2CO3x_goal
        self.DIC = self.c_H2CO3x_goal # initial guess
        self.carbonic = AcidAq(pKa = carbonic_pKa, 
                               charge = 0,
                               conc = self.DIC)
        # (c) Add carbonic species to system
        #     System.species is a tuple:
        self.species = self.species + (self.carbonic,) 
        
        
    def _DICguess(self, DIC):
        """Update system composition with a new DIC value.
        
        This re-calculates the pH of the system, and also the concentrations
        of all carbonic acid species.

        Parameters
        ----------
        DIC : float
            Set the total dissolved inorganic carbon concentration in the
            system..

        Returns
        -------
        float
            The difference between the system concentration of H2CO3x and
            the concentration of H2CO3x required for being at equilibrium
            with atmospheric CO2 according to Henry's law.
        """
        self.DIC = DIC
        self.carbonic.conc = self.DIC
        super().pHsolve() # solve for pH using System.pHsolve
        alphas = self.carbonic.alpha(self.pH)
        self.c_H2CO3x_calc = alphas[0]*self.carbonic.conc
        return (self.c_H2CO3x_calc - self.c_H2CO3x_goal)

    
    def pHsolve(self):
        """Solve the system equilibrium including atmospheric carbon dioxide.
        
        This is done by solving for dissolved inorganic carbon (DIC) in 
        equilibrium with atmosphere. At each guess for the DIC, the system
        pH and composition are re-calculated, and the dissolved
        CO2(aq)/H2CO3(aq) concentrations compared to the air-equilibrium
        values.
        """
        # Use the secant method to find the DIC for which the calculated
        # [H2CO3x] concentration matches the Henry equilibrium concentration.
        #
        #TODO: Optimize initial guess and tolerance parameters. Make
        # informed choices. The currently chosen values seem to work for
        # most compositions.
        DIC_opt = secant(self._DICguess, 
                         x0 = 0.1,
                         x1 = 0.5,
                         tol = 1e-8,
                         ftol = 1e-8)
        # Re-update to be sure to have the optimal DIC value.
        self._DICguess(DIC_opt)
