from containerClass import container
import numpy as np
from scipy.integrate import solve_ivp
import igraph
from calcTw import calcTw
from calcQ import calcQ
from containerClass import container
#from systemDefinitions import system,nodes,nodes1,pipeSections,pipeSections1
from systemDefinitions import system,pipeSections,pipeSections1, pipeSections3,orificeDiam3,pipeSections4

class pipeNetwork:
    def __init__(self):
        self.t = igraph.Graph(directed = True)
        #self.t = igraph.GraphBase(directed = True)

        #attributes of the pipe network
        self.t["agent"] = ""
        self.t["discharge_time"] = 0
        self.t["cyl_valve_type"] = ""
        self.t["cyl_pressure"] = 0
        self.t["cyl_size"] = 0
    

        #attributes of the nodes (name is the _id or the provided number of the node)
        self.t.vs['type'] = ""
        self.t.vs['x']   = 0  #x-coordinate
        self.t.vs['y']   = 0  #y-coordinate
        self.t.vs['z']   = 0  #z-coordinate
        self.t.vs['M']   = 0  #Mach number
        self.t.vs['T0']  = 0 #total temperature
        self.t.vs['T']   = 0  #temperature
        self.t.vs['P0']  = 0 #total pressure
        self.t.vs['P']   = 0  #pressure
        self.t.vs['rho'] = 0 #density
        self.t.vs['D']   = 0 #diameter
        self.t.vs['MFR'] = 0 #mass flow rate
        
        #attributes of the edges
        self.t.es['L'] = 0 #length
        self.t.es['D'] = 0 #internal diameter
        self.t.es['H'] = 0 #elevation change
        self.t.es['Sch'] = 0 #Schedule
        self.t.es['Elb'] = 0 #number of elbows
        self.t.es['Stee'] = False #starts with side tee?
        self.t.es['Ttee'] = False #start with a through tee?
        self.t.es['Cpl'] = 0 #number of coupling or unions
        self.t.es['Dtrp'] = 0 #number of dirt traps
        self.t.es['Ptap'] = 0 #number of pressure taps
        self.t.es['SV'] = 0 #number of selector valves     
        self.t.es['f'] = 0 #friction factor

    def addSystem(self, _system):
        self.t["agent"] = _system["agent"]
        self.t["discharge_time"] = _system["discharge_time"]
        self.t["cyl_valve_type"] = _system["cyl_valve_type"]
        self.t["cyl_pressure"] = _system["cyl_pressure"]
        self.t["cyl_size"] = _system["cyl_size"]

    def addNode(self, _id, _type=None, _x=None, _y=None, _z=None, _M=None, _P=None, _T=None, _P0=None, _T0=None, _rho=None,_Diam = 0):
        self.t.add_vertex(_id, type=_type, x =_x, y=_y, z=_z, M=_M, P=_P, T=_T, P0=_P0, T0=_T0, rho=_rho, D = _Diam)

    def addPipe(self,_start, _end, _L, _D, _H, _Sch, _Elb, _Stee, _Ttee, _Cpl, _Dtrp, _Ptap, _SV, _f):
        sourceIndex = self.t.vs.select(name_eq = _start)[0].index #it should be bases on _start we need to look up the index from vertices
        targetIndex = self.t.vs.select(name_eq = _end)[0].index #same as above
        self.t.add_edge(sourceIndex, targetIndex, L=_L, D=_D, H=_H, Sch=_Sch, Elb=_Elb, Stee=_Stee, Ttee= _Ttee, Cpl=_Cpl, Dtrp=_Dtrp, Ptap=_Ptap, SV=_SV, f=_f)

    # def addAllNodes(self, _nodes):
    #     for i in _nodes:
    #         self.addNode(i[0],i[1],i[2],i[3],i[4],0,0,0,0,0,0)
    def addAllNodes(self, _pipes, _orificeDiams):
        _nodes=[]
        for i in _pipes:
            if i[0] not in _nodes:
                _nodes.append(i[0])
            if i[1] not in _nodes:
                _nodes.append(i[1])
        for n in _nodes:
            diam = 0
            if n in _orificeDiams:
                diam = _orificeDiams[n]
            self.addNode(n, _Diam= diam)

    def addAllPipes(self, _pipes, _orificeDiams):
        self.addAllNodes(_pipes, _orificeDiams)
        for i in _pipes:
            self.addPipe(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11],i[12],0)

    def topoSummary(self):
        igraph.summary(self.t)
        print(self.t.get_edgelist())
        print(self.t.farthest_points(directed=True))
        for i in self.t.vs:
            print(i['name'],i.index,i.degree(),i.degree(mode='OUT'),i.degree(mode='IN'))

    def plot(self):
        layout = self.t.layout("kk")
        self.t.vs["label"]=self.t.vs.indices#["index"]#["D"]
        self.t.es["label"]=self.t.es.indices#["L"]
        igraph.plot(self.t, bbox = (400,400), layout = layout)

    def calcNetwork(self):#calculates the network at the current instace of time

        #find the end of manifold node
        #create two networks net1=cylinder bank until the end of manifold node net2=pipe network from the end of manifold node until the nozzles
        #find all the nozzles in net2
        #create a queue for all tees in net2
        #calculate the pressure from each nozzle to the upstream tee or manifold node
        #if p is calculated from both branches of the tee then calculate the pressure at that tee by iterating over the mass flow rates
        #move upstream until all the tees are calculated and have known properties
        #start from the most remote cylinder and calculate the pressure downstream until all the properties on tees are known
        #compare the pressure of the manifold node from net1 and net2
        #Iterate until the pressure of the manifold node is equal
        pass 


# net = pipeNetwork()
# net.addSystem(system)
# net.addAllPipes(pipeSections3, orificeDiam3)
# net.topoSummary()
# net.plot()