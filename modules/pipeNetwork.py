from containerClass import container
import igraph

class pipeNetwork:
    def __init__(self,
                cylSize = 140,
                systemType = "inergen"):
        self.t = igraph.Graph()
        self.t["system"] = systemType
        self.t["cylSize"] = cylSize
        self.addNode("0")
        self.addNode("1")
        self.t.add_vertex(name = "0")
        self.t.add_vertex(name = "1")
        # self.t.vs[0]["x"] = 0
        # self.t.vs[0]["y"] = 0
        # self.t.vs[0]["z"] = 0
        # self.t.vs[0]["M"] = 0
        # self.t.vs[0]["P"] = 0
        # self.t.vs[0]["T"] = 0
        # self.t.vs[0]["P0"] = 0
        # self.t.vs[0]["T0"] = 0
        # self.t.vs[0]["rho"] = 0
        # self.t.es[0][]

    def addNode(self, _name, _x, _y, _z, _M, _P, _T, _P0, _T0, _rho):
        self.t.add_vertex(name = _name , "x" = _x, "y"=_y)
        
    def topoSummary(self):
        igraph.summary(self.t)

    def plot(self):
        layout = self.t.layout("kk")
        igraph.plot(self.t, layout = layout)

net = pipeNetwork()
net.topoSummary()
net.plot()



