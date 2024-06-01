from gurobipy import *
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
from time import *
import os, psutil
from utilities import *

def model_DCOM(stars,amount_stars,degrees,edge_lenght,max_lenght=True):
    '''
    creates a net with the DCOM model.
    Arguments:
        stars: a list of dictionaries with the stars data
        amount_stars: number of stars for the net (integer)
        d: degree for every node member of the net
        edge_lenght: the max edge distance allowable in our net if max_lenght
        max_lenght: tell us if are usign maximum allowable lenght in the algorithm
    Returns:
        The net in a list of edges, the elapsed time, the memory used by the algorithm.
    '''
    st = time()
    solution_net = []
    nodes = [i for i in range(amount_stars)]

    # Create the graph edges. Don't allow loops in the net 
    edges = tuplelist([(i,j) for i in nodes for j in nodes if i < j])
    triangles = tuplelist([(i,j,k) for i in nodes for j in nodes for k in nodes if i<j and j<k])
    
    
    #distances
    distances = {}
    for k,l in edges:
        distances[k,l] =  dist1(k,l,stars)
  
    # Model
    model = Model("DCOM")

    # decision variables
    x = model.addVars(edges, vtype=GRB.BINARY, name='x')
    # 3-cycle variable
    ytr = model.addVars(triangles, vtype=GRB.CONTINUOUS, name='ytr')
    #degree
    pd = model.addVars(degrees, vtype=GRB.INTEGER, name='pd')
    
    # amount of triangles in our net
    pntr = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'pntr')
    
    #cluster coeficient
    pcc = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'pcc')

    # slack variables for degree
    sdm = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'sdm')
    sdp = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'sdp')
    # average path lenght
    pacc = model.addVar(vtype=GRB.CONTINUOUS, name = 'pacc')
    #Global cluster coefficient
    pgcc = model.addVar(vtype=GRB.CONTINUOUS, name = 'pgcc')
    #slack variables for the average cluster coefficient
    saccm = model.addVar(vtype=GRB.CONTINUOUS, name = 'saccm')
    saccp = model.addVar(vtype=GRB.CONTINUOUS, name = 'saccp')
    #slack varibles for the global cluster coefficient
    sgccm = model.addVar(vtype=GRB.CONTINUOUS, name = 'sgccm')
    sgccp = model.addVar(vtype=GRB.CONTINUOUS, name = 'sgccp')



    # Objetive function
    model.setObjective(quicksum(sdm[i] + sdp[i] for i in nodes) + saccm + saccp + sgccm + sgccp, GRB.MINIMIZE)

    # Constraints
    model.addConstrs(ytr[i,j,k] <= x[i,j] for i,j,k in triangles)
    model.addConstrs(ytr[i,j,k] <= x[i,k] for i,j,k in triangles)
    model.addConstrs(ytr[i,j,k] <= x[j,k] for i,j,k in triangles)
    model.addConstrs(ytr[i,j,k] >= x[i,j] + x[i,k] + x[j,k] - 2 for i,j,k in triangles)
    model.addConstrs(ytr[i,j,k] >= 0 for i,j,k in triangles)

    model.addConstrs(pd[i] == quicksum(x[j,k] for j,k in edges if j==i or k==i) for i in nodes)
    model.addConstrs(pntr[i] == quicksum(ytr[j,k,l] for j,k,l in triangles if j == i or k == i or l == i) for i in nodes)
    model.addConstrs(pcc[i] == (1/m.comb(degrees[i],2))*pntr[i] for i in nodes)
    model.addConstr(pgcc == (3/(sum(m.comb(degrees[i],2) for i in nodes)))*quicksum(pntr[i] for i in nodes))
    model.addConstr(pacc == (1/(len(nodes)))*quicksum(pcc[i] for i in nodes))
    model.addConstrs(degrees[i] <= pd[i] + sdm[i] - sdp[i] for i in nodes)
    model.addConstrs(pd[i] + sdm[i] - sdp[i] <= degrees[i] for i in nodes)
    model.addConstr(pacc + saccm - saccp >= 0)
    model.addConstr(pacc + saccm - saccp <= 0.25)
    model.addConstr(pgcc + sgccm - sgccp >= 0)
    model.addConstr(pacc + saccm - saccp <= 0.25 )

   # break symmetry to avoid losing performance with possible network isomorphisms
   # (exactly the same networks but with permuted nodes)
    model.addConstrs(pd[i] - pd[i+1] >= 0 for i in nodes if i <= 8)
    model.addConstrs(pd[i] - pd[i+1] + pcc[i] - pcc[i+1] >= 0 for i in nodes if i <= 8)
    # max lenght edges constraint
    if max_lenght:
        model.addConstrs(distances[k,l]*x[k,l] <= edge_lenght
                        for k,l in edges)
   

    model.addConstrs(sdm[i] >= 0 for i in nodes)
    model.addConstrs(sdp[i] >= 0 for i in nodes)
    model.addConstr(saccm >= 0)
    model.addConstr(saccp >= 0)
    model.addConstr(sgccm >= 0)
    model.addConstr(sgccp >= 0)

    model.optimize()

    end_t = time()
    elapsed_time = end_t - st
    for i,j in edges:
      if x[i,j].x > 0.2:
        solution_net.append((i,j))

    process = psutil.Process(os.getpid())
    mem_us = process.memory_info().rss/1024**2

    return solution_net, elapsed_time, mem_us

# -----------------------------------------------------------------------------
# Example function usage
am_nodes = []
short_path = []
mem_us = []
time_list = []
dist_list = []

for amount_stars in range(10,101,10):
    # a dictionary with the node sequence for the net
    degrees = {}
    for i in range(amount_stars):
      degrees[i] = 3

    stars = import_database('./base_final.csv',amount_stars)
    tras_index, stars_fil = add_tras(stars,amount_stars)
    solution_net,elapsed_time, mem_usage = model_DCOM(stars_fil,amount_stars,degrees,20)
   #create the nodes of the network with the stars and add traspist at the end
    time_list.append(elapsed_time)
    mem_us.append(mem_usage)
    if (0,tras_index) in solution_net:
      solution_net.remove((0,tras_index))
    #if not (0,tras_index) in solution_net:
    #  print("no esta")

    net_plot(stars_fil,amount_stars,solution_net)
    d, s_path = calculate_s_path(stars_fil,amount_stars,solution_net, tras_index)
    net_plot(stars_fil,amount_stars,s_path)
    d1 = dist1(0,tras_index,stars_fil)
    short_path.append(d)
    dist_list.append(d1)
    am_nodes.append(amount_stars)
    plot_histograms(stars, amount_stars, solution_net)
print(dist_list)
print(short_path)
plot_stat(am_nodes,mem_us,time_list,short_path)
plt.plot(dist_list,short_path)
plt.xlabel("Grado")
plt.ylabel("Distancia Sol-Trappist (pc)")
plt.show()



