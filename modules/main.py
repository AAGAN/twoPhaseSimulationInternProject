# overall algorithm is going to be in this file
from thermo import mixture
from fluids.units import *
from thermo.units import Chemical

m = mixture.Mixture(['nitrogen','argon','carbon dioxide'], zs=[.52, .4, .08], T = 300, P=1e5)

print("Cp = {:.4f}, Cv = {:.4f}, mu = {:.8f}, MW = {:.4f}, SG = {:.4f}, Cp/Cv = {:.4f}, Z = {:.4f}".format(m.Cp, m.Cvg, m.mu, m.MW, m.SG, m.isentropic_exponent, m.Z))

m = mixture.Mixture(['nitrogen','argon','carbon dioxide'], zs=[.52, .4, .08], T = 250, P=1e5)

print("Cp = {:.4f}, Cv = {:.4f}, mu = {:.8f}, MW = {:.4f}, SG = {:.4f}, Cp/Cv = {:.4f}, Z = {:.4f}".format(m.Cp, m.Cvg, m.mu, m.MW, m.SG, m.isentropic_exponent, m.Z))

#iterates over time(s)
#calculated the condision of the tank for each time step
#creates instances of pipeNetwork for each time step using the calculated tank condition
#calls the pipeNetwork.calcNetwork() method to calculate the tank at that time step
#calculates the condition of the tank for the next time step and goes to step 3