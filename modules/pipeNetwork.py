from containerClass import container
import numpy as np
from scipy.integrate import solve_ivp
import igraph
from calcTw import calcTw
from calcQ import calcQ
from containerClass import container
from orificeForward import main_nozzle_forward, mass_critical
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

        self.C1 = 1e-5 #c1 parameter for the ratio function 
        self.C2 = 1e-5 #c2 parameter for the ratio function
        self.Err = 0.1 #error in calculating deltaP0
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

    def plot(self,vertexLabel,edgeLabel):
        layout = self.t.layout("kk")
        #self.t.vs["label"]=self.t.vs.indices#["index"]#["D"]
        #self.t.es["label"]=self.t.es.indices#["L"]
        #igraph.plot(self.t, bbox = (1400,1400), layout = layout)
        edge_labels = []
        vertex_labels = []
        visual_style = {}
        visual_style["vertex_size"] = 15#40
        visual_style["vertex_color"] = "red"#[color_dict[gender] for gender in g.vs["gender"]]
        
        if vertexLabel == "index":
            visual_style["vertex_label"] = self.t.vs.indices#["index"]#g.vs["name"]
        else:
            for i in range(0,len(self.t.vs.indices)):
                vertex_labels.append([self.t.vs[i].index , self.t.vs[i][vertexLabel]])
            visual_style["vertex_label"] = vertex_labels
        
        if edgeLabel == "index":
            visual_style["edge_label"] = self.t.es.indices
        else:
            for i in range(0,len(self.t.es.indices)):
                edge_labels.append([self.t.es[i].index , self.t.es[i][edgeLabel]])
            visual_style["edge_label"] = edge_labels
                
        visual_style["vertex_label_size"] = 10#25
        visual_style["edge_width"] = 2#5#[1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]
        visual_style["edge_label_size"] = 10#25
        visual_style["layout"] = layout
        visual_style["bbox"] = (400,400)#(800, 800)
        visual_style["margin"] = 50
        # print(edge_labels)
        # print(vertex_labels)
        return igraph.plot(self.t, **visual_style)

    def findNext(self, node, edge):
        #function to find the next tee, nozzle, tank, commonNode from a node and a pipe based on the direction of the pipe
        #if the returned degree findNext()[1] is 1:tank or nozzle, 0:commonNode, 3:tee
        _common = self.commonNode
        _edge = edge
        _node = node
        _graph = self.t
        _nodeIDs = [node.index]
        if _node.index == _common:
            return _nodeIDs, 0
        if _node.index == _edge.source:
            #move forward in the graph until the next tee
            _nextNodeID = _edge.target
            _nodeIDs.append(_nextNodeID)
            while _graph.vs[_nextNodeID].degree()==2 and _graph.vs[_nextNodeID].index != _common:
                _nextNodeID = _graph.vs[_nextNodeID].successors()[0].index
                _nodeIDs.append(_nextNodeID)
                
        elif _node.index == _edge.target:
            #move backwards in the graph until the next tee
            _nextNodeID = _edge.source
            _nodeIDs.append(_nextNodeID)
            while _graph.vs[_nextNodeID].degree()==2 and _graph.vs[_nextNodeID].index != _common:
                _nextNodeID = _graph.vs[_nextNodeID].predecessors()[0].index
                _nodeIDs.append(_nextNodeID)
        else:
            return "edge is not connected to the node!"
        _nodeType = _graph.vs[_nextNodeID].degree()
        if _nextNodeID == _common:
            _nodeType = 0
        return _nodeIDs, _nodeType

    #Assume 50% division at every node in the beginning
    def divide(self, ratio2, P_01, P_02,c1,c2, error):
        del_P = P_02-P_01
        rat_old= ratio2
        if del_P > error:
            rat_new = rat_old* np.exp(-c1*del_P)
        elif del_P < -error:
            rat_new = 1 - np.exp(c2*del_P) + rat_old*np.exp(c2*del_P)
        return rat_new

    #propagate mass flow rates downstream a tee
    def propagateMFR(self,node, edge1, edge2):
        P01 = edge1['P0i']
        P02 = edge2['P0i']
        MFR1 = edge1['MFR']
        MFR2 = edge2['MFR']
        MFR_total = MFR1 + MFR2
        ratio = MFR2 / MFR_total
        if np.abs(P01-P02)>self.Err:
            ratio_new = self.divide(ratio,P01,P02,self.C1,self.C2, self.Err)
            ratio_new = 0.1
            MFR1_new = (1 - ratio_new)*MFR_total 
            MFR2_new = ratio_new*MFR_total
            
            #distribute the MFR1_new to the downstream nozzles evenly
            nextNodeInEdge1Direction = self.findNext(node,edge1) #find the next node in edge1 direction
            if nextNodeInEdge1Direction[1] == 3:#if it is a tee
                allPassesFromEdge1 = self.t.get_all_shortest_paths(nextNodeInEdge1Direction[0][-1],self.t.vs.select(_outdegree = 0))
                allNozzleIndicesFromEdge1 = [L[-1] for L in allPassesFromEdge1]
                for noz in allNozzleIndicesFromEdge1:#allPassesFromEdge1[][-1]:#for noz in all the nozzles after this tee
                    self.t.vs[noz]['MFR'] = MFR1_new / len(allNozzleIndicesFromEdge1)
            else:#if it is a nozzle
                self.t.vs[nextNodeInEdge1Direction[0][-1]]['MFR'] = MFR1_new
            
            #distribute the MFR2_new to the downstream nozzles evenly
            nextNodeInEdge2Direction = self.findNext(node,edge2)
            if nextNodeInEdge2Direction[1] == 3:
                allPassesFromEdge2 = self.t.get_all_shortest_paths(nextNodeInEdge2Direction[0][-1],self.t.vs.select(_outdegree = 0))
                allNozzleIndicesFromEdge2 = [L[-1] for L in allPassesFromEdge2]
                for noz in allNozzleIndicesFromEdge2:#allPassesFromEdge2[][-1]:
                    self.t.vs[noz]['MFR'] = MFR2_new / len(allNozzleIndicesFromEdge2)
            else:
                self.t.vs[nextNodeInEdge2Direction[0][-1]]['MFR'] = MFR2_new
                
            #set the 'calculated' property of all the nodes downstream the node to False    
            for i in self.t.get_all_shortest_paths(node.index,self.t.vs.select(_outdegree = 0)):
                for j in i:
                    self.t.vs[j]['calculated'] = False
                    
            #set the P0i of all the edges downstream the node to 0
            for i in self.t.get_all_shortest_paths(node.index,self.t.vs.select(_outdegree = 0)):
                for k in range(0,len(i)-1):
                    self.t.es.select(_source = i[k],_target = i[k+1])['P0i'] = 0
        else:
            node['calculated'] = True
            P0 = (P01+P02)/2.0
            node['P0']=P0
            edge3 = self.t.es.select(_source = self.t.vs[node.index].predecessors()[0].index , _target = node.index)[0]
            previousNode = self.findNext(node,edge3) 
            #calculate the properties for all the nodes until the next node backwards
            #set all the P0i for edges until the next node and set 'calculated' property of all nodes until the next node to True

    
    def calcNode(self,source, target, graph):
        '''
            function to calculate the pressure drop between two nodes (from source to target)
            calculates all the properties on the path from the source to the target node 
            saves all the properties on the pipes and nodes. except for MFR and P0
            there should be only one path between source and target nodes
        '''
        #find all the nodes between the source and the target
        # nodes = graph.get_shortest_paths(source.index, target.index,mode='ALL')
        # for i in range(0,len(nodes)-1):
        #     edge = graph.es.select(_source = i, _target = i+1)
        #     #if we're going in the direction of the flow, then edge['L'] is positive, otherwise we need to multiply it by (-1)
        #     calcEdge = calcQ(edge['L'],...)
        pass

    def forwardPass(self):
        #tank1 = net.t.vs[firstTank]
        # tank1Valve = net.t.vs[firstTank].successors()[0]
        # initialMdot0Guess = mdot0
        # g.es.select(_source = tank1.index, _target = tank1Valve.index)[0]['MFR'] = initialMdot0Guess #initial guess by considering a small dt for the container function
        # nextEdge = g.es.select(_source = tank1Valve.index, _target = tank1Valve.successors()[0].index)
        # if tank1Valve.index == commonNode:
        #     pass
        #     #calculate the properties for the common node and store it in that node
        # else:
        #     nextNode = findNext(tank1Valve,nextEdge[0],g,commonNode)
        #     #calculate the properties for the next node and store the properties on the next node
        #     calcNode(tank1Valve,nextNode[0][-1],g)
        #     while nextNode[0][-1] != commonNode:
        #         pass
        #         #find the oposite edge
        #         #find the next tank in the oposite edge direction
        #         #guess a mfr value for this tank
        #         #calculate the properties until the node
        #         #compare the P01 and P02 using the function ratio
        #         #if dp is good then assign average p0 to the node 
        #         #add the mfr from both inputs to the next edge
        #         #find the next node
        #         #calculate the properties for the next node and store the properties on that next node except MFR and P0, MFR and P0 should be saved only on the pipes
        pass

    def backwardPass(self):
        #Backward pass implementation
        # g.vs[commonNode]['MFR'] = 1
        # totalMFR = g.vs[commonNode]['MFR']
        # nozzles = g.vs.select(_outdegree = 0)
        # for i in nozzles:
        #     i['MFR'] = totalMFR / len(nozzles)
        #     i['calculated'] = False
        #     #print(i['MFR'])

        # for noz in nozzles:
        #     if noz['calculated'] == False:
        #         print(findNext(noz,g.es.select(_target = noz.index)[0],g,commonNode))
        #         path = findNext(noz,g.es.select(_target = noz.index)[0],g,commonNode)
        #         calcNode(g.vs[path[0][0]],g.vs[path[0][-1]],g)
        pass

    def calcNetwork(self):#calculates the network at the current instace of time
        #create a queue for all tees in net2
        #calculate the pressure from each nozzle to the upstream tee or manifold node
        #if p is calculated from both branches of the tee then calculate the pressure at that tee by iterating over the mass flow rates
        #move upstream until all the tees are calculated and have known properties
        #start from the most remote cylinder and calculate the pressure downstream until all the properties on tees are known
        #compare the pressure of the manifold node from net1 and net2
        #Iterate until the pressure of the manifold node is equal
        pass