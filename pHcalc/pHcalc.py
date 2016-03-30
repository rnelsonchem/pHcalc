import numpy as np
import scipy.optimize as spo

class Neutral(object):
    """A nonreactive ion class.

    This object defines things like K+ and Cl-, which contribute to the
    overall charge balance, but do not have any inherent reactivity with
    water.
    
    Parameters
    ----------
    charge : int
        The formal charge of the ion.

    conc : float
        The concentration of this species in solution.

    Attributes
    ----------
    charge : int
        The formal charge of the ion.

    conc : float
        The concentration of this species in solution.

    """
    def __init__(self, charge=None, conc=None):
        if charge == None:
            raise ValueError(
                "The charge for this ion must be defined.")

        self.charge = charge 
        self.conc = conc

    def alpha(self, pH):
        '''Return the fraction of each species at a given pH.

        Parameters
        ----------
        pH : int, float, or Numpy Array
            These are the pH value(s) over which the fraction should be
            returned.

        Returns
        -------
        Numpy NDArray
            Because this is a non-reactive ion class, this function will
            always return a Numpy array containing just 1.0's for all pH
            values.

        '''
        if isinstance(pH, (int, float)):
            length = 1
        else:
            length = len(pH)
        ones = np.ones(length).reshape(-1,1)
        return ones


        
class Acid(object):
    '''An acidic species class.

    This object is used to calculate a number of parameters related to a weak
    acid in an aqueous solution. 
    
    Parameters
    ----------
    Ka : None (default), float, list, Numpy Array
        This defines the Ka values for all acidic protons in this species. It
        can be a single Ka value (float), a list of floats, or a Numpy array
        of floats. Either this value or pKa needs to be defined. The other
        will then be calculated from the given values.

    pKa : None (default), float, list, Numpy Array
        The pKa value(s) for all the acidic protons in this species.  This
        follows the same rules as Ka (See Ka description for more details),
        and either this value or Ka must be defined.

    charge : None (default), int
        This is the charge of the fully protonated form of this acid. This
        must be defined.

    conc : None (default), float
        The formal concentration of this acid in solution. This value must be
        defined.

    Note
    ----
    There is no corresponding Base object. To define a base, you must use a
    combination of an Acid and Neutral object. See the documentation for
    examples.

    '''
    def __init__(self, Ka=None, pKa=None, charge=None, conc=None):
        # Do a couple quick checks to make sure that everything has been
        # defined.
        if Ka == None and pKa == None:
            raise ValueError(
                "You must define either Ka or pKa values.")
        elif charge == None:
            raise ValueError(
                "The maximum charge for this acid must be defined.")

        # Make sure both Ka and pKa are calculated. For lists of values, be
        # sure to sort them to ensure that the most acidic species is defined
        # first.
        elif Ka == None:
            if isinstance(pKa, (int, float)):
                self.pKa = np.array( [pKa,], dtype=float)
            else:
                self.pKa = np.array(pKa, dtype=float)
                self.pKa.sort()
            self.Ka = 10**(-self.pKa) 
        elif pKa == None:
            if isinstance(Ka, (int, float)):
                self.Ka = np.array( [Ka,], dtype=float)
            else:
                self.Ka = np.array(Ka, dtype=float)
                # Ka values must be in reverse sort order
                self.Ka.sort()
                self.Ka = self.Ka[::-1]
            self.pKa = -np.log10(self.Ka)
        # This temporary Ka array will be used to calculate alpha values. It
        # starts with an underscore so that it won't be confusing for others.
        self._Ka_temp = np.append(1., self.Ka)
        
        # Make a list of charges for each species defined by the Ka values.
        self.charge = np.arange(charge, charge - len(self.Ka) - 1, -1)
        # Make sure the concentrations are accessible to the object instance.
        self.conc = conc 

    def alpha(self, pH):
        '''Return the fraction of each species at a given pH.

        Parameters
        ----------
        pH : int, float, or Numpy Array
            These are the pH value(s) over which the fraction should be
            returned.

        Returns
        -------
        Numpy NDArray
            These are the fractional concentrations at any given pH. They are
            sorted from most acidic species to least acidic species. If a
            NDArray of pH values is provided, then a 2D array will be
            returned. In this case, each row represents the speciation for
            each given pH.
        '''
        # If the given pH is not a list/array, be sure to convert it to one
        # for future calcs.
        if isinstance(pH, (int, float)):
            pH = [pH,]
        pH = np.array(pH, dtype=float)

        # Calculate the concentration of H3O+. If multiple pH values are
        # given, then it is best to construct a two dimensional array of
        # concentrations.
        h3o = 10.**(-pH)
        if len(h3o) > 1:
            h3o = np.repeat( h3o.reshape(-1, 1), len(self._Ka_temp), axis=1)

        # These are the powers that the H3O+ concentrations will be raised.
        power = np.arange(len(self._Ka_temp))
        # Calculate the H3O+ concentrations raised to the powers calculated
        # above (in reverse order).
        h3o_pow = h3o**( power[::-1] )
        # Calculate a cumulative product of the Ka values. The first value
        # must be 1.0, which is why _Ka_temp is used instead of Ka.
        Ka_prod = np.cumproduct(self._Ka_temp)
        # Multiply the H3O**power values times the cumulative Ka product.
        h3o_Ka = h3o_pow*Ka_prod

        # Return the alpha values. The return signature will differ is the
        # shape of the H3O array was 2-dimensional. 
        if len(h3o.shape) > 1:
            den = h3o_Ka.sum(axis=1)
            return h3o_Ka/den.reshape(-1,1)
        else:
            den = h3o_Ka.sum()
            return h3o_Ka/den


class System(object):
    '''An object used to define an a system of acid and neutral species.

    This object accepts an arbitrary number of acid and neutral species
    objects and calculates the pH of the system. 

    The objects need to completely define the system. All weak species are
    input as their fully protonated acidic form as an Acid object. Any neutral
    ions are accounted for by the appropriate Neutral objects. For example,
    NaHCO3 would be desribed by an Acid object for H2CO3 and a Neutral object
    for Na+, and both of these objects would have the same initial
    concentrations.

    The pH solving is done using a simple minimization algorithm which
    minimizes the difference in the total positive and negative ion
    concentrations in the system. A good initial guess will help here, and it
    can be set manually using the 'guess' keyword. There is an automated
    method that can be run as well if you set the guess_est argument. This
    will override whatever you pass is for 'guess'. The 'est_num' keyword sets
    the number of data points that you'd like to use for finding the guess
    estimate. Too few points might start you pretty far from a good spot; too
    many points is probably overkill and won't help much. This may or may not
    speed things up.

    Parameters
    ----------
    *species : These are the Acid and Neutral objects that you'd like to use
    to define your system.

    '''
    def __init__(self, *species):
        self.species = species


    def _diff_pos_neg(self, pH):
        twoD = True
        if isinstance(pH, (int, float)) or pH.shape[0] == 1:
            twoD = False
        else:
            pH = np.array(pH, dtype=float)
        # Calculate the h3o and oh concentrations and sum them up.
        h3o = 10.**(-pH)
        oh = (10.**(-14))/h3o
        x = (h3o - oh)

        # Go through all the species that were given, and sum up their
        # charge*concentration values into our total sum.
        for s in self.species:
            if twoD == False:
                x += (s.conc*s.charge*s.alpha(pH)).sum()
            else:
                x += (s.conc*s.charge*s.alpha(pH)).sum(axis=1)
        
        # Return the absolute value so it never goes below zero.
        return np.abs(x)
        

    def pHsolve(self, guess=7.0, guess_est=False, est_num=1500, 
                method='Nelder-Mead', tol=1e-5):
        if guess_est == True:
            phs = np.linspace(0, 14, est_num)
            guesses = self._diff_pos_neg(phs)
            guess_idx = guesses.argmin()
            guess = phs[guess_idx]
            
        self.pHsolution = spo.minimize(self._diff_pos_neg, guess, 
                method=method, tol=tol)
        
        if self.pHsolution.success == False:
            print('Warning: Unsuccessful pH optimization!')
            print(self.pHsolution.message)
            
        if len(self.pHsolution.x) == 1:
            self.pH = self.pHsolution.x[0]



if __name__ == '__main__':
    # (NH4)3PO4
    # H3PO4 input as the fully acidic species
    a = Acid(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    # NH4+ again input as fully acidic species
    # The concentration is 3x greater than the phosphoric acid
    b = Acid(pKa=9.498, charge=1, conc=3.e-3)
    #k = neutral(charge=1, conc=1.e-4)
    #k = neutral(charge=-1, conc=3.e-3)
    s = System(a, b)
    s.pHsolve()
    print('(NH4)3PO4 1e-3 M pH = ', s.pH)

    # HCl, just need to define the amount of Cl-. The solver will find the
    # correct H+
    b = Neutral(charge=-1, conc=1e-8)
    s = System(b)
    s.pHsolve()
    print('HCl 1e-8 M pH = ', s.pH)
    
    # KOH, just need to define the amount of K+, solver takes care of the
    # rest.
    a = Neutral(charge=+1, conc=0.1)
    s = System(a)
    s.pHsolve()
    print('NaOH 0.1 M pH = ', s.pH)
    
    # Distribution diagram H3PO4
    import matplotlib.pyplot as plt
    a = Acid(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    pH = np.linspace(0, 14, 1000)
    plt.plot(pH, a.alpha(pH))
    plt.show()
    
    # Estimate Best pH
    # This is done internallly by the pHsolve function if you use the 
    # guess_est=True flag
    s = System(a)
    diffs = s._diff_pos_neg(pH)
    plt.plot(pH, diffs)
    plt.show()

    # Phosphoric Acid Titration Curve
    Na_concs = np.linspace(1.e-8, 5.e-3, 500)
    H3PO4 = Acid(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    phs = []
    for conc in Na_concs:
        Na = Neutral(charge=1, conc=conc)
        s = System(H3PO4, Na)
        s.pHsolve(guess_est=True)
        phs.append(s.pH)
    plt.plot(Na_concs, phs)
    plt.show()


