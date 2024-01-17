import numpy as np
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
import os, psutil
from time import time
from utilities import *

def erdos(stars, p):
  st = time()
  solution_net = []
  
  am_nodes.append(amount_stars)
  for i in range(amount_stars-1):
    for j in range(i+1, amount_stars):
      chance = np.random.rand()
      #p = dist1(i,j,stars_fil)/12.1
      #if p > 1:
       # p = 0 # el camino es mas grande en distancia que el total a Trass
      if p > chance:
        solution_net.append((i,j))

  process = psutil.Process(os.getpid())
  if (0,amount_stars-1) in solution_net:
    solution_net.remove((0,amount_stars-1))
  if not (0,amount_stars-1) in solution_net:
    print("no esta")

  end_t = time()
  elapsed_time = end_t - st
  process = psutil.Process(os.getpid())
  mem_usage = process.memory_info().rss/1024**2
  
  
  return solution_net, elapsed_time, mem_usage



am_nodes = []
time_list = []
mem_us = []
short_path = []
dist_list = []
p = 0.8

for amount_stars in range(10,101,10):
  stars = import_database_ran('./base_final.csv',100)
  solution_net, elapsed_time, mem = erdos(stars,p)
  time_list.append(elapsed_time)
  mem_us.append(mem)
  net_plot(stars,amount_stars,solution_net)
  d, s_path = calculate_s_path(stars,amount_stars,solution_net, amount_stars-1)
  net_plot(stars,amount_stars,s_path)
  d1 = dist1(0,amount_stars-1,stars)
  short_path.append(d)
  dist_list.append(d1)

  plot_histograms(stars, amount_stars, solution_net)
print(am_nodes)
print(dist_list)
print(short_path)
plot_stat(am_nodes,mem_us,time_list,short_path)