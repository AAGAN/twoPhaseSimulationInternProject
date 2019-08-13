from containerClass import container
import numpy as np
from scipy.integrate import solve_ivp
import igraph
from calcTw import calcTw
from calcQ import calcQ
from containerClass import container
#from systemDefinitions import system,nodes,nodes1,pipeSections,pipeSections1
from systemDefinitions import system,pipeSections,pipeSections1, pipeSections3,orificeDiam3,pipeSections4,pipeSections0,orificeDiam0

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
        self.t.vs['calculated'] = False
        
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
        self.t.es['MFR'] = 0 #mass flow rate
        self.t.es['P0i'] = 0 #initial total pressure calculated

        #self.nozzles #stores the vertex sequence of all the nozzles
        #self.tanks #stores the vertex sequence of all the tanks
        #self.firstTank #index of the first tank vertex
        #self.lastNozzle  #index of the last nozzle vertex 
        #self.commonNode # index of the common node

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
        #After all the pipes and nodes are added, process the pipe data information 
        self.processNetwork()
        self.findCommonNode()

    def processNetwork(self):
        if len(self.t.vs.select(_outdegree_gt = 2)) > 0 :
            print('only junctions with 2 outlets are accepted')
        if len(self.t.vs.select(_indegree_gt = 2)) > 0:
            print('only junctions with 2 inlets are accepted')

        self.t.vs.select(_outdegree = 0)['type'] = 'nozzle'
        self.t.vs.select(_indegree = 0)['type'] = 'tank'
        self.t.vs.select(_outdegree = 2)['type'] = 'tee'
        self.t.vs.select(_indegree = 2)['type'] = 'tee'
        self.t.vs.select(_indegree = 1, _outdegree = 1)['type'] = 'coupling'
        self.nozzles = self.t.vs.select(_outdegree = 0)
        self.tanks = self.t.vs.select(_indegree = 0)
        
        farthest_points = self.t.farthest_points(directed = True)
        self.firstTank = farthest_points[0]
        self.lastNozzle = farthest_points[1]

    def findCommonNode(self):
        commonNodes = []
        for node in self.t.vs:
            if node.index not in self.nozzles.indices and node.index not in self.tanks.indices:
                Common = True
                for noz in self.nozzles.indices:
                    for tank in self.tanks.indices:
                        if self.t.edge_connectivity(node.index,noz)== 0 or self.t.edge_connectivity(tank, node.index)==0:
                            Common = False
                            break
                if Common == True:
                    commonNodes.append(node.index)
        self.commonNode = commonNodes[-1]

    def topoSummary(self):
        igraph.summary(self.t)
        print(self.t.get_edgelist())
        print(self.t.farthest_points(directed=True))
        for i in self.t.vs:
            print(i['name'],i.index,i.degree(),i.degree(mode='OUT'),i.degree(mode='IN'))

    def plot(self):
        layout = self.t.layout("kk")
        #self.t.vs["label"]=self.t.vs.indices#["index"]#["D"]
        #self.t.es["label"]=self.t.es.indices#["L"]
        #igraph.plot(self.t, bbox = (1400,1400), layout = layout)
        visual_style = {}
        visual_style["vertex_size"] = 40
        visual_style["vertex_color"] = "red"#[color_dict[gender] for gender in g.vs["gender"]]
        visual_style["vertex_label"] = self.t.vs["name"]#.indices#["index"]#g.vs["name"]
        visual_style["vertex_label_size"] = 25
        visual_style["edge_label"] = self.t.es['L']
        visual_style["edge_width"] = 5#[1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]
        visual_style["edge_label_size"] = 25
        visual_style["layout"] = layout
        visual_style["bbox"] = (800, 800)
        visual_style["margin"] = 50
        igraph.plot(self.t, **visual_style)

    def forwardPass(self):
        pass

    def backwardPass(self):
        pass

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