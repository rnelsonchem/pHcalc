pHcalc
######

*pHcalc* is a Python library for systematically calculating solution pH,
distribution diagrams, and titration curves.

This package is Python 3.5 compatible with dependencies_ only on Numpy and
Scipy. If you will be plotting the data, then there is an `optional
dependency`_ as well.  

Bug fixes, questions, and update requests are encouraged and can be
filed at the `GitHub repo`_. 

.. _dependencies:

Dependencies
------------

* Numpy >= 1.10

* Scipy >= 0.17

.. _optional dependency:

Optional Packages
-----------------

* Matplotlib >= 1.5

Installation
------------

*pHcalc* is only a single Python file, so installation is quite simple. After
installation of the dependencies, the most recent version of *pHcalc* is
available via ``pip``, either from PyPI_ (stable) or the `GitHub repo`_ (most
recent).

From PyPI_::

    $ pip install pHcalc

From the `GitHub repo`_::

    $ pip install git+https://github.com/rnelsonchem/pHcalc.git

Of course, the GitHub installation requires that you have installed Git on
your computer.



Background
##########

*pHcalc* calculates the pH of a complex system of potentially strong and weak
acids and bases using a systematic equilibrium solution method. This method is
described in detail in `the Journal of Chemical Education`_ and in this
`ChemWiki article`_, for example. (There is also another, older Pascal program
called PHCALC_, which uses matrix algebra to accomplish the same task. To the
best of my knowledge, the source code for this program is no longer
available.)

Basically, this method calculates the fractional distribution of a potentially
polyprotic weak acid species at a given pH value. Multiplying this by the
concentration of acid in solution provides the concentration of each acidic
species in the system. The optimium pH for this system is found when
charge balance is achieved, i.e. the concentrations of positively charged ions
equals the charge for the negatively charged ions. 

Using this methodology bases and strong acids can be described using neutral,
charged species. These are ions that do not react with water, such as |Na+|
and |Cl-|. In this context, any |Cl-| in solution must be charged balanced
with an appropriate amount of |H3O|, which would define HCl in solution.
|Na+| must be offset by an equivalent amount of |OH-|, which defines a
solution of NaOH. A 1:1 combination of |Na+| and |H2CO3| would describe a
solution of |NaHCO3|.

Example Usage
#############

*pHcalc* defines three classes - Acid, Neutral, and System - which are used in
calculating the pH of the system. |H3O| and |OH-| are never input. The
|H3O| concentration is adjusted internally, and |OH-| is calculated using K\
:sub:`W`\ .

.. code:: python

    >>> from pHcalc.pHcalc import Acid, Neutral, System
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt # Optional for plotting below

pH of 0.01 M HCl
----------------

First of all, HCl completely dissociates in water to give equal amounts of
|H3O| and |Cl-|, so all you need to define is |Cl-|. 

.. code:: python

    >>> cl = Neutral(charge=-1, conc=0.01)
    >>> system = System(cl)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 1.9999

pH of 1e-8 M HCl
----------------

This is a notoriously tricky example for introductory chemistry students;
however, *pHcalc* handles it nicely.

.. code:: python

    >>> cl = Neutral(charge=-1, conc=1e-8)
    >>> system = System(cl)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 6.978295898 (NOT 8!)

pH of 0.01 M NaOH
-----------------

This example is very similar to our HCl example, except that our Neutral
species must have a positive charge.

.. code:: python

    >>> na = Neutral(charge=1, conc=0.01)
    >>> system = System(na)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 12.00000

pH of 0.01 M HF
---------------

Here we will use an Acid object instance to define the weak acid HF, which has
a |Ka| of 6.76e-4 and a |pKa| of 3.17. You can use either value when you
create the Acid instance. When defining an Acid species, you must always
define a ``charge`` keyword argument, which is the charge of the *fully
protonated species*.

.. code:: python

    >>> hf = Acid(Ka=6.76e-4, charge=0, conc=0.01)
    >>> # hf = Acid(pKa=3.17, charge=0, conc=0.01) will also work
    >>> system = System(hf)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 2.6413261

pH of 0.01 M NaF
----------------

This system consist of a 1:1 mixture of an HF Acid instance and a |Na+|
Neutral instance. The System can be instantiated with an arbitrary number of
Acids and Neutral objects.

.. code:: python

    >>> hf = Acid(Ka=6.76e-4, charge=0, conc=0.01)
    >>> na = Neutral(charge=1, conc=0.01)
    >>> system = System(hf, na)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 7.5992233


pH of 0.01 M |H2CO3|
--------------------

The |Ka| and |pKa| attributes can also accept lists of values for polyprotic
species.

.. code:: python

    >>> carbonic = Acid(pKa=[3.6, 10.32], charge=0, conc=0.01)
    >>> system = System(carbonic)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 2.8343772

pH of 0.01 M Alanine Zwitterion Form
------------------------------------

Alanine has two pKa values, 2.35 and 9.69, but the fully protonated form is
positively charged. In order to define the neutral zwitterion, the Acid object
needs to be combined with a positively-charged Neutral species as well, which
would represent one equivalent of NaOH.

.. code:: python 

    >>> ala = Acid(pKa=[2.35, 9.69], charge=1, conc=0.01)
    >>> na = Neutral(charge=1, conc=0.01)
    >>> system = System(ala, na)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 6.0991569

pH of 0.01 M |NH4PO4|
---------------------

This is equivalent to a 1:3 mixture of |H3PO4| and |NH4|, both of which are
defined by Acid objects.

.. code:: python

    >>> phos = Acid(pKa=[2.148, 7.198, 12.319], charge=0, conc=0.01)
    >>> nh4 = Acid(pKa=9.25, charge=1, conc=0.01*3)
    >>> system = System(phos, nh4)
    >>> system.pHsolve()
    >>> print(system.pH) # Should print 8.95915298

Distribution Diagrams
---------------------

Acid objects also define a function called ``alpha``, which calculates the
fractional distribution of species at a given pH. This function can be used to
create distribution diagrams for weak acid species. ``alpha`` takes a single
argument, which is a single pH value or a Numpy array of values. For a single
pH value, the function returns a Numpy array of fractional distributions
ordered from most acid to least acidic species. 

.. code:: python

    >>> phos = Acid(pKa=[2.148, 7.198, 12.319], charge=0, conc=0.01)
    >>> phos.alpha(7.0)
    array([ 8.6055e-06, 6.1204e-01, 3.8795e-01, 1.8611e-06])
    >>> # This is H3PO4, H2PO4-, HPO4_2-, and HPO4_3-

For a Numpy array, a 2D array of fractional distribution values is returned,
where each row is a series of distributions for each given pH. The 2D returned
array can be used to plot a distribution diagram. 

.. code:: python

    >>> phos = Acid(pKa=[2.148, 7.198, 12.319], charge=0, conc=0.01)
    >>> phs = np.linspace(0, 14, 1000)
    >>> fracs = phos.alpha(phs)
    >>> plt.plot(phs, fracs)
    >>> plt.show()

Titration Curves
----------------

Using a simple loop, we can also construct arbitrary titration curves as well.
In this example, we will titrate |H3PO4| with NaOH. The ``guess_est`` keyword
argument for the ``System.pHsolve`` method forces the calculation of a best
guess for starting the pH optimization algorithm. This may speed up the
evaluation of the pH and can also be used if the minimizer throws an error
during the pH calculation. 

.. code:: python

    >>> na_concs = np.linspace(1e-8, 5.e-3, 500)
    >>> phos = Acid(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    >>> phs = []
    >>> for conc in Na_concs:
    >>>     na = Neutral(charge=1, conc=conc)
    >>>     system = System(phos, na)
    >>>     system.pHsolve(guess_est=True)
    >>>     phs.append(system.pH)
    >>> plt.plot(na_concs, phs)
    >>> plt.show()

.. Substitutions


.. |Na+| replace:: Na\ :sup:`+`
.. |Cl-| replace:: Cl\ :sup:`-`
.. |H3O| replace:: H\ :sub:`3`\ O\ :sup:`+`
.. |OH-| replace:: OH\ :sup:`-`
.. |H2CO3| replace:: H\ :sub:`2`\ CO\ :sub:`3`
.. |NaHCO3| replace:: NaHCO\ :sub:`3`
.. |Ka| replace:: K\ :sub:`a`
.. |pKa| replace:: pK\ :sub:`a`
.. |NH4PO4| replace:: (NH\ :sub:`4`\ )\ :sub:`3`\ PO\ :sub:`4`
.. |H3PO4| replace:: H\ :sub:`3`\ PO\ :sub:`4`
.. |NH4| replace:: NH\ :sub:`4`\ :sup:`+`

.. External Hyperlinks

.. _GitHub repo: https://github.com/rnelsonchem/pHcalc
.. _PyPI: https://pypi.python.org/pypi/pHcalc
.. _the Journal of Chemical Education:
      http://pubs.acs.org/doi/abs/10.1021/ed100784v
.. _ChemWiki article: 
    http://chemwiki.ucdavis.edu/Core/Analytical_Chemistry/Analytical_Chemistry_2.0/06_Equilibrium_Chemistry/6G%3A_Solving_Equilibrium_Problems#6G.3_A_Systematic_Approach_to_Solving_Equilibrium_Problems
.. _PHCALC: http://pubs.acs.org/doi/pdf/10.1021/ed071p119
