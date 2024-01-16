from gurobipy import *
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
from time import *
import os, psutil
from utilities import *

def model_triangles(stars,amount_stars,d):
    st = time()
    solution_net = []
    nodes = [i for i in range(amount_stars)]

    # creamos las aristas o edges del grafo. No permitimos que existan bucles en el grafo
    edges = tuplelist([(i,j) for i in nodes for j in nodes if i < j])
    triangles = tuplelist([(i,j,k) for i in nodes for j in nodes for k in nodes if i<j and j<k])
    degrees = {}
    for i in range(amount_stars):
        degrees[i] = 3

    # Model
    model = Model("NetGen")

    # variables de decision
    x = model.addVars(edges, vtype=GRB.BINARY, name='x')
    ytr = model.addVars(triangles, vtype=GRB.CONTINUOUS, name='ytr')

    pd = model.addVars(degrees, vtype=GRB.INTEGER, name='pd')
    pntr = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'pntr')
    pcc = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'pcc')
    sdm = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'sdm')
    sdp = model.addVars(degrees, vtype=GRB.CONTINUOUS, name = 'sdp')

    pacc = model.addVar(vtype=GRB.CONTINUOUS, name = 'pacc')
    pgcc = model.addVar(vtype=GRB.CONTINUOUS, name = 'pgcc')
    saccm = model.addVar(vtype=GRB.CONTINUOUS, name = 'saccm')
    saccp = model.addVar(vtype=GRB.CONTINUOUS, name = 'saccp')
    sgccm = model.addVar(vtype=GRB.CONTINUOUS, name = 'sgccm')
    sgccp = model.addVar(vtype=GRB.CONTINUOUS, name = 'sgccp')



    # funcion objetivo
    model.setObjective(quicksum(sdm[i] + sdp[i] for i in nodes) + saccm + saccp + sgccm + sgccp, GRB.MINIMIZE)

    # restricciones
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

    # romper simetria para evitar perder performance con los posibles isomorfismos de redes
    #(redes exactamente iguales pero con los nodes permutados)
    model.addConstrs(pd[i] - pd[i+1] >= 0 for i in nodes if i <= 8)
    model.addConstrs(pd[i] - pd[i+1] + pcc[i] - pcc[i+1] >= 0 for i in nodes if i <= 8)


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

am_nodes = []
short_path = []
mem_us = []
time_list = []
dist_list = []

for amount_stars in range(10,20):
    stars = import_database('./base_final.csv',amount_stars)
    tras_index, stars_fil = add_tras(stars,amount_stars)
    solution_net,elapsed_time, mem_usage = model_triangles(stars_fil,amount_stars,8)
    #creamos los nodes de la red con las estrellas y agregamos traspist al final
    time_list.append(elapsed_time)
    mem_us.append(mem_usage)
    if (0,tras_index) in solution_net:
      solution_net.remove((0,tras_index))
    if not (0,tras_index) in solution_net:
      print("no esta")

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



