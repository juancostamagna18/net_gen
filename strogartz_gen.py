import numpy as np
import math as m
import matplotlib.pyplot as plt
from graph_tool.all import *
import os, psutil
from time import time
from utilities import *

def stro_gen(stars, pw, amount_stars,kl,edge_lenght, max_lenght=False):
    '''
    creates a net with the Watts-Strogatz model.
    Arguments:
        stars: a list of dictionaries with the stars data
        pw: the probability of disconect an edge(optional float in case have max_lenght = true)
        amount_stars: number of stars for the net (integer)
        kl: maximum node degree in the net
        edge_lenght: maximum allowable lenght for every edge in the net
        max_lenght: tell us if are usign maximum allowable lenght in the algorithm
    Returns:
        The net in a list of edges, the elapsed time, the memory used by the algorithm.
    '''

    edges_act = []
    st = time()
    for i in range(amount_stars-1):
      for j in range(i+1, (i+m.ceil(kl/2))+1):
        if j > amount_stars-1:
          j = j - amount_stars
        edges_act.append((i,j))

    for i in range(amount_stars-1):
      for l in range(i+1,(i+m.ceil(kl/2))+1):
        if l > amount_stars-1:
          l = l - amount_stars
        chance = np.random.rand()
        if max_lenght:
            pw = 1 - dist1(i,l,stars)/edge_lenght
            if abs(pw) > 1:
              pw = 0
        #print(pw)
        if pw > chance:
          j = np.random.randint(0,amount_stars-1)
          while j == i or (i,j) in edges_act or (j,i) in edges_act:
            j = np.random.randint(0,amount_stars-1)
          if (i,l) in edges_act:
            edges_act.remove((i,l))
          edges_act.append((i,j))

    end_t = time()
    elapsed_time = end_t - st
    process = psutil.Process(os.getpid())
    mem_us = process.memory_info().rss/1024**2

    return edges_act, elapsed_time, mem_us




am_nodes = []
mem_us = []
time_list = []
dist_list = []
short_path = []

for amount_stars in range(10,101,10):
    stars = import_database_ran('./base_final.csv', amount_stars)
    tras_index, stars = add_tras(stars,amount_stars)
    if amount_stars > 20:
      solution_net,memo_us, elapsed_time = stro_gen(stars,0.2,amount_stars,15,12.1)
    else:
      solution_net,memo_us, elapsed_time = stro_gen(stars,0.2,amount_stars,5,12.1)
    if (0,amount_stars-1) in solution_net:
      solution_net.remove((0,amount_stars-1))

    time_list.append(elapsed_time)
    am_nodes.append(amount_stars)
    mem_us.append(mem_us)
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