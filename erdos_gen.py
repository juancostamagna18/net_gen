import numpy as np
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
import os, psutil
from time import time
from utilities import *

def erdos_gen(stars,p, amount_stars,edge_lenght, max_lenght=False):
    '''
    creates a net with the Erdos model.
    Arguments:
        stars: a list of dictionaries with the stars data
        p: the probability of disconect an edge(optional float in case have max_lenght = true)
        amount_stars: number of stars for the net (integer)
        edge_lenght: maximum allowable lenght for every edge in the net
        max_lenght: tell us if are usign maximum allowable lenght in the algorithm
    Returns:
        The net in a list of edges, the elapsed time, the memory used by the algorithm.
    '''
    edges_act = []
    st = time()
    for i in range(amount_stars-1):
        for j in range(i+1, amount_stars):
            chance = np.random.rand()
            if max_lenght:
                p = dist1(i,j,stars)/edge_lenght
                if p > 1:
                    p = 0 # el camino es mas grande en distancia que el total a Trass
            if p > chance:
                edges_act.append((i,j))
    end_t = time()
    elapsed_time = end_t - st
    process = psutil.Process(os.getpid())
    mem_us = process.memory_info().rss/1024**2

    return edges_act, elapsed_time, mem_us

am_nodes = []
time_list = []
mem_us = []
short_path = []
dist_list = []
p = 0.8

for amount_stars in range(10,101,10):
  stars = import_database_ran('./base_final.csv',100)
  am_nodes.append(amount_stars)
  solution_net, elapsed_time, mem = erdos_gen(stars,p,amount_stars,12.1)
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