from graph_tool.all import *
import os, psutil
from time import *
from gurobipy import *
import numpy as np
import math as m
import matplotlib.pyplot as plt
from utilities import *


def model_SPDOM(stars,amount_stars,min_d,max_d,min_diam,max_diam,
                min_cpl,max_cpl,edge_lenght,max_lenght=True):
    '''
    creates a net with the SPDOM model.
    Arguments:
        stars: a list of dictionaries with the stars data
        amount_stars: number of stars for the net (integer)
        min_d: lower bound degree for every node member of the net
        max_d: upper bound degree for every node member of the net
        min_diam: lower bound of the net diameter
        max_diam: upper bound of the net diameter
        min_cpl: upper bound of the net characteristic path length
        max_cpl: upper bound of the net characteristic path length
        edge_lenght: the max edge distance allowable in our net if max_lenght
        max_lenght: tell us if are usign maximum allowable lenght in the algorithm
    Returns:
        The net in a list of edges, the elapsed time, the memory used by the algorithm.
    '''
    st = time()
    edges_ac = []
    nodes = [i for i in range(amount_stars)]
    edges = [(i,j) for i in nodes for j in nodes if i < j]
    DL = min_diam
    DU = max_diam
    dL = min_d
    dU  = max_d
    E = amount_stars*(amount_stars -1)/2

    cplL = min_cpl
    cplU = max_cpl

    #distances
    distances = {}
    for k,l in edges:
        distances[k,l] =  dist1(k,l,stars)
    
    # Model
    model = Model("SPDOM")


    # decision variables
    x = model.addVars(edges, vtype=GRB.BINARY, name='x')
    #degree
    pd = model.addVars(nodes, vtype=GRB.INTEGER, name='pd')
    # degree slack variavbles
    sdm = model.addVars(nodes, vtype=GRB.CONTINUOUS, name = 'sdm')
    sdp = model.addVars(nodes, vtype=GRB.CONTINUOUS, name = 'sdp')
    #shortest path varibles that represents the edges in the path
    f = model.addVars(nodes,nodes,edges,vtype=GRB.CONTINUOUS, name='f')
    # shortest path variable
    w = model.addVars(edges,vtype=GRB.CONTINUOUS, name='w')
    # binary vars that allows that the characteristic path lenght be calculated in a range
    rcplm = model.addVars(edges, vtype=GRB.BINARY, name = 'rcplm')
    rcplp = model.addVars(edges, vtype=GRB.BINARY, name = 'rcplp')
    #charecteristic path lenght variable
    pcpl = model.addVar(vtype=GRB.CONTINUOUS, name = 'pcpl')
    # diameter slack variables
    sDm = model.addVar(vtype=GRB.CONTINUOUS, name = 'sDm')
    sDp = model.addVar(vtype=GRB.CONTINUOUS, name = 'sDp')
    # characteristic path lenght variables
    scplm = model.addVar(vtype=GRB.CONTINUOUS, name = 'scplm')
    scplp = model.addVar(vtype=GRB.CONTINUOUS, name = 'scplp')
    # binary variable that allows to include the lenght path equal 1
    phi = model.addVars(edges, vtype=GRB.BINARY, name = 'phi')

    # objetive function
    model.setObjective(quicksum(sdm[i] + sdp[i] for i in nodes) +
                    scplm + scplp + sDm + sDp , GRB.MINIMIZE)

    # constraints
    model.addConstrs(pd[i] == quicksum(x[j,k] for j,k in edges
                  if j==i or k==i) for i in nodes)

    model.addConstrs(f[k,l,i,j] + f[l,k,i,j] <= x[k,l]
                  for k,l in edges for i,j in edges)
    model.addConstrs(f[k,l,i,j] >= 0 for k in nodes
                  for l in nodes for i,j in edges)
    model.addConstrs(f[k,k,i,j] == 0 for k in nodes for i,j in edges)
    for i,j in edges:
        model.addConstrs( sum(f[k,l,i,j] for l in nodes) -
                          sum(f[l,k,i,j] for l in nodes) ==
                          (1 if k==i else -1 if k==j else 0)
                          for k in nodes)
    # max lenght edges constraint
    if max_lenght:
        model.addConstrs(distances[k,l]*f[k,l,i,j] <= edge_lenght
                        for k,l in edges for i,j in edges)
        

    model.addConstrs( w[i,j] == quicksum(f[k,l,i,j] for k,l in edges)
                          for i,j in edges)

    model.addConstrs(-(m.ceil(amount_stars/2) - 1)*rcplp[i,j] <= w[i,j] - pcpl
                          for i,j in edges)
    model.addConstrs(w[i,j] - pcpl <=  (amount_stars - 2)*(1-rcplp[i,j])
                         for i,j in edges)

    model.addConstrs(-(m.ceil(amount_stars/2) - 1)*(1-rcplm[i,j]) <= w[i,j] - pcpl
                          for i,j in edges)
    model.addConstrs(w[i,j] - pcpl <=  (amount_stars - 2)*rcplm[i,j]
                         for i,j in edges)

    model.addConstr(2*quicksum(rcplp[i,j] for i,j in edges) ==
                        (E if E % 2 == 0 else E+1))
    model.addConstr(2*quicksum(rcplm[i,j] for i,j in edges) ==
                        (E if E % 2 == 0 else E+1))


    model.addConstrs(dL- sdm[i] <= pd[i] for i in nodes)
    model.addConstrs(pd[i] <= dU + sdp[i] for i in nodes)

    model.addConstr(cplL <= pcpl + scplm - scplp)
    model.addConstr(pcpl + scplm - scplp <=  cplU)

    model.addConstrs(1+(DL - 1)*phi[i,j] <= w[i,j] + sDm - sDp for i,j in edges)
    model.addConstrs(w[i,j] + sDm - sDp <= DU for i,j in edges)

    model.addConstr(quicksum(phi[i,j] for i,j in edges) >= 1)
    model.addConstrs(w[i,amount_stars-1] >= w[i+1,amount_stars-1] for i in nodes
                     if i < amount_stars-2)

    model.addConstrs(sdm[i] >= 0 for i in nodes)
    model.addConstrs(sdp[i] >= 0 for i in nodes)
    model.addConstr(scplm >= 0)
    model.addConstr(scplp >= 0)
    model.addConstr(sDm >= 0)
    model.addConstr(sDp >= 0)

    model.optimize()
    end_t = time()
    elapsed_time = end_t - st
    process = psutil.Process(os.getpid())
    mem_us = (process.memory_info().rss/1024**2)
    print(process.memory_info().rss)  # in bytes
    #print("Funcion Objetivo: ", str(round(model.ObjVal,2)))
    for i,j in edges:
        if x[i,j].x > 0.2:
            edges_ac.append((i,j))
    return edges_ac,elapsed_time,mem_us

# -----------------------------------------------------------------------------
# Example function usage
am_nodes = []
short_path = []
mem_us = []
time_list = []
is_random = False

quan = 25
min_degree = 0
max_degree = 10
min_diam = 0
max_diam = 12.1
min_cpl = 1
max_cpl = 12

stars = import_database('./base_final.csv',quan)
for amount_stars in range(5,quan,5):
  am_nodes.append(amount_stars)
  if is_random:
    stars = import_database_ran('./base_final.csv',quan)
  # create the nodes of the network with the stars and add traspist at the end
  tras_index, stars_fil = add_tras(stars,amount_stars)
  solution_net,elapsed_time,mem_usage = model_SPDOM(stars_fil,amount_stars,min_degree,
                                                    max_degree,min_diam,max_diam,min_cpl,max_cpl,12)
  if (0,tras_index) in solution_net:
    solution_net.remove((0,tras_index))
  # if not (0,tras_index) in solution_net:
  #   print("no esta")
  net_plot(stars_fil,amount_stars,solution_net)
  time_list.append(elapsed_time)
  mem_us.append(mem_usage)
  d, s_path = calculate_s_path(stars_fil,amount_stars,solution_net,tras_index)
  short_path.append(d)
  net_plot(stars_fil,amount_stars,s_path)
  plot_histograms(stars_fil, amount_stars,solution_net)
  #average path lenght using graph tools
  '''G = Graph(directed=False)
  G.add_vertex(amount_stars)
  edge_weight = G.new_edge_property("double")
  for (i,j) in edges_ac:
    e = G.add_edge(i,j)
    edge_weight[e] = dist1(i,j,stars_fil)
    print(shoMelissa McCracken rtest_distance(G,source=0,target=amount_stars-1,weights=edge_weight))
    dist = shortest_distance(G,weights=edge_weight)
    apl = sum([sum(i) for i in dist])/(G.num_vertices()**2-G.num_vertices())
    print(amount_stars, apl)'''

print(mem_us,time_list)
plot_stat(am_nodes,mem_us,time_list,short_path)

