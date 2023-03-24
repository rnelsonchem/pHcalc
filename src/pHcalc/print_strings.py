no_ph = '''### THE CONCENTRATIONS OF THIS SYSTEM ARE NOT AT EQUILIBRIUM ###
To determine the equilibrium species distribution use System.pHsolve\n\n'''

has_ph = '''### THESE ARE THE EQUILIBRIUM SYSTEM CONCENTRATIONS ###

SYSTEM pH: {0.pH:.3f}\n\n'''

sep_line1 = '='*65 + '\n'

sep_line2 = '-'*65 + '\n'

header_line = f"{'Species':15}{'Charge':10}{'Ka':15}{'pKa':10}{'Conc':15}\n"

acid_line = '{0:15}{1:<+10d}{2:<15.3e}{3:<10.2f}{4:<15.4e}\n'

ion_line = '{0:15}{1:<+35d}{2:<15.4e}\n'
