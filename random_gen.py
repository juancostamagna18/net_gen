import numpy as np
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
import os, psutil
from time import *
from functools import reduce
from operator import and_
from utilities import *

d_list = []
mem_us = []
time_list = []
am_nodes = []


for amount_stars in range(10,101,10):
  selected_arcs = []
  all_connected = False
  stars = import_database_ran('./base_final.csv', amount_stars)
  tras_index, stars_fil = add_tras(stars,amount_stars)

  G1 = Graph(directed=False)
  G1.add_vertex(amount_stars)
  edge_weight = G1.new_edge_property("double")
  visited = [False for l in range(amount_stars)]
  st_time = time()
  process = psutil.Process(os.getpid())

  while ( not all_connected):
      #pos = np.random.randint(0,len(edges)+1)
      #arc = edges[pos]
      i = np.random.randint(0,amount_stars)
      j = np.random.randint(0,amount_stars)
      while(i >= j):
        i = np.random.randint(0,amount_stars)
        j = np.random.randint(0,amount_stars)
      arc = (i,j)
      if arc not in selected_arcs:
        selected_arcs.append(arc)
        e = G1.add_edge(arc[0],arc[1])
        d = dist1(arc[0],arc[1],stars)
        edge_weight[e] = d
      
        for k in range(amount_stars):
          list_nei = G1.get_in_neighbors(k)
          for n in list_nei:
            if not visited[n]:
              visited[n] = True
        # to know if all vertices are connected use reduce in visited list 
        # is the same as for all 
        all_connected = reduce(and_,visited)
  end_time = time()
  elapsed_time = end_time - st_time
  mem = process.memory_info().rss/1024**2
  #dist = shortest_distance(G1,0,tras_index,weights=edge_weight)
  dist, path = calculate_s_path(stars_fil,amount_stars,selected_arcs,tras_index)
  d_list.append(dist)
  am_nodes.append(amount_stars)
  mem_us.append(mem)
  time_list.append(elapsed_time)
  plot_histograms(stars_fil,amount_stars,selected_arcs)
  net_plot(stars_fil,amount_stars,path)
  net_plot(stars_fil,amount_stars,selected_arcs)

print(d_list)
plot_stat(am_nodes, mem_us, time_list,d_list)