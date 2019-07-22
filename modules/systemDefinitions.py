system = {"agent":"inergen", "discharge_time":120}
cylinder = {"valve_type":"iflow", "size" : 140, "pressure": 3e7}

###########################################################################3
##a system with one cylinder and two nozzles
nodes = [
    [1   ,"cyl"       , 0  , 0 , 0],#type, x,y,z
    [2   ,"cyl_valve" , 0  , 0 , 1],
    [3   ,"elbow"     , 0  , 0 , 5],
    [4   ,"tee"       , 10 , 0 , 5],
    [301 ,"nozzle"    , 10 , 0 , 4],
    [302 ,"nozzle"    , 15 , 0 , 4]
]

pipeSections = [
    [ 1 , 2   ,  1 ,   1 ,  1 , 80 ,0,0,0,0,0,0,0],#startNode,endNode,len,diam,height,Sch,Elb,Stee,Ttee,Cpl,Dtrp,Ptap,SV
    [ 2 , 3   ,  4 ,   1 ,  4 , 80 ,0,0,0,0,0,0,0],
    [ 3 , 4   , 10 ,   1 ,  0 , 80 ,1,0,0,0,0,0,0],
    [ 4 , 301 ,  1 , 0.5 , -1 , 80 ,0,1,0,0,0,0,0],
    [ 4 , 302 ,  5 , 0.5 , -1 , 80 ,1,0,1,0,1,0,0]
]

###########################################################################3
##a system with two cylinders and two nozzles
nodes1 = [
    [1   ,"cyl"       , 0  , 0 , 0],#type, x,y,z
    [2   ,"cyl valve" , 0  , 0 , 1],
    [3   ,"elbow"     , 0  , 0 , 6],
    [4   ,"tee"       , 1  , 0 , 6],
    [5   ,"cyl"       , 1  , 0 , 0],
    [6   ,"cyl valve" , 1  , 0 , 1],
    [10  ,"common"    , 2  , 0 , 6],
    [7   ,"tee"       , 7  , 0 , 6],
    [8   ,"nozzle"    , 7  , 0 , 5],
    [9   ,"nozzle"    , 10 , 0 , 5]
]

pipeSections1 = [
    [ 1 , 2   ,  1 ,   1 ,  1 , 80 ,0,0,0,0,0,0,0],#startNode,endNode,len,diam,height,Sch,Elb,Stee,Ttee,Cpl,Dtrp,Ptap,SV
    [ 2 , 3   ,  5 ,   1 ,  5 , 80 ,1,0,0,0,0,0,0],
    [ 3 , 4   ,  1 ,   1 ,  0 , 80 ,1,0,0,0,0,0,0],
    [ 5 , 6   ,  1 ,   1 ,  1 , 80 ,0,0,0,0,0,0,0],
    [ 6 , 4   ,  5 ,   1 ,  5 , 80 ,1,0,0,0,0,0,0],
    [ 4 , 10  ,  1 ,   1 ,  0 , 80 ,0,1,1,0,0,0,0],
    [ 10, 7   ,  5 ,   1 ,  0 , 80 ,0,0,0,1,0,0,0],
    [ 7 , 8   ,  1 ,   1 ,  -1, 80 ,0,1,0,0,0,0,0],
    [ 7 , 9   ,  3 ,   1 ,  -1, 80 ,1,0,1,0,1,0,0]
]