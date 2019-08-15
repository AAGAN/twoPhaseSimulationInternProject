# overall algorithm is going to be in this file
from thermo import mixture
from fluids.units import *
from thermo.units import Chemical
from pipeNetwork import pipeNetwork
from containerClass import container
from systemDefinitions import system,pipeSections3,orificeDiam3,pipeSections0,orificeDiam0

m = mixture.Mixture(['nitrogen','argon','carbon dioxide'], zs=[.52, .4, .08], T = 300, P=1e5)
#print("Cp = {:.4f}, Cv = {:.4f}, mu = {:.8f}, MW = {:.4f}, SG = {:.4f}, Cp/Cv = {:.4f}, Z = {:.4f}".format(m.Cp, m.Cvg, m.mu, m.MW, m.SG, m.isentropic_exponent, m.Z))

net = pipeNetwork()
net.addSystem(system)
net.addAllPipes(pipeSections3, orificeDiam3)
print(net.tanks.indices)
print(net.nozzles.indices)

final_time = 120.0 #seconds
dt = 1.0 #second
initialMass = 32.2 #kg of agent in each tank
cylInitialTemp = 295.0
agentInitialTemp = 294.0
cylOrificeDiam = 0.0066565
cylL = 1.7
cylD = 0.267
cylWallThickness = 0.005
ambientTemp = 300
initialBackP = 1e5
numTimeSteps = 200
cylInitialPressure = initialMass / (np.pi * cylD**2.0 / 4.0 * cylL) * 8.3145 / m.MW * agentInitialTemp
m = mixture(['nitrogen','argon','carbon dioxide'], zs=[.52, .4, .08], T = agentInitialTemp, P=cylInitialPressure)
cont = [0] * len(net.tanks)

def updatePipeNetwork(oldNetwork):
    #create a new pipenetwork based on the caculated values  of the previous time
    newNetwork = oldNetwork
    return newNetwork

time = 0
net.calcNetwork()
nets = []
while time<final_time:
    time = time + dt
    net = updatePipeNetwork(net)
    net.calcNetwork()
    if time % 10 == 0:
        nets.append(net)