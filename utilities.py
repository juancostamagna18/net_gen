import csv
import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys
from graph_tool.all import *
from time import *
import os, psutil

def dist1(p1,p2,stars_fil):
    '''
    Ecludian distance between two stars.
    Arguments:
        p1: a position in the node list (integer)
        p2: a position in the node list (integer)
        stars_fil: a dictionaries list that contains all the stars data
    Returns:
        The distance between the two stars in pc(float)
    '''
    result = 0
    result += m.pow(stars_fil[p1]['x']-stars_fil[p2]['x'],2)
    result += m.pow(stars_fil[p1]['y']-stars_fil[p2]['y'],2)
    result += m.pow(stars_fil[p1]['z']-stars_fil[p2]['z'],2)
    return m.sqrt(result)

def dist(p1,p2):
    '''
    Ecludian distance between two stars.
    Arguments:
        p1: star 1 coordinates in pc(float)
        p2: star 2 coordinates in pc(float)
    Returns:
        The distance between the two stars in pc(float)
    '''
    result = 0
    for i in range(len(p1)):
        result += m.pow(p2[i]-p1[i],2)
    return m.sqrt(result)

def import_database_ran(filename,amount_stars):
  '''
  Imports the csv database file in a python list of dictionaries
  selecting the stars randomly, sorted in increasing order of distance from the sun.
  Arguments:
      filename: string with the relative or absolut path to the database file
      amount_stars: the number of stars in the net
  Returns:
      The sorted stars list with the data.
  '''
  stars = []
  # load the csv file with stars data
  file_name = filename # 'hygdata_v3.csv'
  file = open(file_name, 'rt')
  reader = csv.DictReader(file)

  # adding the data to the dictionary
  for row in reader:
      stars.append(
          {'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z']), 'id': row['id'],
           'prop': row['prop']})

  file.close()
  stars_ran = []
  i = 1 # the sun is always there
  selected = []
  while i <= amount_stars:
    num = np.random.randint(1,5000)
    if num not in selected:
      stars_ran.append(stars[num])
      selected.append(num)
      i += 1

  # saving the data for the sun
  origen = stars[0:1][0]
  
  # insert the sun in the first place
  stars_ran.insert(0,origen)

  return stars_ran

def import_database(filename,amount_stars):
    '''
    Imports the csv database file in a python list of dictionaries, sorted in increasing order of distance from the sun.
    Arguments:
        filename: string with the relative or absolut path to the database file
    Returns:
        The sorted stars list with the data.
    '''
    stars = []
    # load the csv file with stars data
    file_name = filename # 'hygdata_v3.csv'
    file = open(file_name, 'rt')
    reader = csv.DictReader(file)

    # adding the data to the dictionary
    for row in reader:
        stars.append(
            {'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z']), 'id': row['id'],
             'prop': row['prop']})

    file.close()

    # saving the data for the sun
    origen = stars[0:1][0]
    # sorting the data
    # sorts only the indexes not equals to the sun index
    stars = sorted(stars[1:amount_stars],
        key = lambda d:dist([origen['x'],origen['y'],origen['z']],[d['x'],d['y'],d['z']]))

    # insert the sun in the first place
    stars.insert(0,origen)

    return stars

def net_plot(stars,amount_stars,solution_net):
    '''
    plots the stars net (nodes and edges).
    Arguments:
        stars: the list with the stars data (list of dictionaries)
        amount_stars: number of the total stars using in our model
        solution_net: edges list taking from the model solution
    '''
    coordx,coordy,coordz = take_coordinates(stars,amount_stars)

    fig = plt.figure(figsize=(80,60))
    axes = plt.axes(projection="3d")
    axes.scatter3D(coordx,coordy,coordz,color="red")
    axes.set_xlabel('X (pc)')
    axes.set_ylabel('Y (pc)')
    axes.set_zlabel('Z (pc)')
    # ploting the generated net
    for i,j in solution_net:
        plt.plot([stars[i]['x'], stars[j]['x']], [stars[i]['y'],stars[j]['y']],[stars[i]['z'],stars[j]['z']])
    # adding the labels to the sun and Trasppist-1
    for i in range(amount_stars):
        if stars[i]['prop'] == "Sol" or stars[i]['prop'] == "Trappist-1":
            label = stars[i]['prop']
        else:
            label = ""
        axes.text(stars[i]['x'], stars[i]['y'], stars[i]['z'], label, None, color = 'b')

    axes.locator_params('x', nbins=12)
    axes.locator_params('y', nbins=12)
    axes.locator_params('z', nbins=12)

    axes.set_xlabel('X')
    axes.set_ylabel('Y')
    axes.set_zlabel('Z')

    plt.show()

def plot_scatter(stars,amount_stars):
    '''
    plots only the nodes of the stars net.
    Arguments:
        stars: the list with the stars data (list of dictionaries)
        amount_stars: number of the total stars using in our model
    '''

    coordx,coordy,coordz = take_coordinates(stars,amount_stars)
    fig = plt.figure(figsize=(32,32))
    axes = plt.axes(projection="3d")
    # plot only the stars without the edges
    axes.scatter3D(coordx,coordy,coordz,color="red")

     # adding the labels to the sun and Trasppist-1
    for i in range(amount_stars):
        if stars[i]['prop'] == "Sol" or stars[i]['prop'] == "Trappist-1":
            label = stars[i]['prop']
        else:
            label = ""
        axes.text(stars[i]['x'], stars[i]['y'], stars[i]['z'], label, None, color = 'b')

    #setting the axes labels
    axes.set_xlabel('X')
    axes.set_ylabel('Y')
    axes.set_zlabel('Z')
    axes.grid(False)

    plt.show()


def take_coordinates(stars, amount_stars):
    '''
    takes only the coordinates from the stars list (needed to plot the stars)
    Arguments:
        stars: the list with the stars data (list of dictionaries)
        amount_stars: number of the total stars using in our model
    Returns:
        the 3 lists with the x,y,z coordinates in pc
    '''

    coordx = []
    coordy = []
    coordz = []
    #need to be a list for matplotlib
    for i in range(amount_stars):
        coordx.append(stars[i]['x'])
        coordy.append(stars[i]['y'])
        coordz.append(stars[i]['z'])

    return coordx,coordy,coordz

def plot_stat(am_nodes, mem_us, time_list,shortest_path):
    '''
    plots all the time, memory and shortest path data relate to the run of the algorithms
    Arguments:
        am_nodes: amount nodes list used in all executions 
        mem_us: memory usage list used in all executions
        time_list: the time elapsed list between the begining and the end in all executions
        shortest_path: the shortest path listin all executions
    '''
    #plot the memory graph
    plt.plot(am_nodes,mem_us)
    plt.xlabel("Cantidad de estrellas")
    plt.ylabel("Uso de Memoria (MB)")
    plt.show()

    # plot the times
    plt.plot(am_nodes,time_list)
    plt.xlabel("Cantidad de estrellas")
    plt.ylabel("Tiempo (s)")
    plt.show()

    #plot the shortest path sun-Trasppist
    plt.plot(am_nodes,shortest_path)
    plt.xlabel("Cantidad de estrellas")
    plt.ylabel("distacia TRAPPIST-Sol (pc)")
    plt.show()


def plot_histograms(stars,amount_stars,solution_net):
    '''
    plots the histograms of several net properties: edge distance, distance to the sun,
    all to all stars distance, node degree, global clustering, local clustering.
    Arguments:
        stars: dictionaries list with the stars data
        amount_stars: number of the total stars using in our model
        solution_net: edges list taking from the model solution
    '''
    # create a undirected graph
    G = Graph(directed=False)
    G.add_vertex(amount_stars)
    edge_weight = G.new_edge_property("double")
    # list for the edges distance
    edges_d = []
    for (i,j) in solution_net:
        e = G.add_edge(i,j)
        d = dist1(i,j,stars)
        edge_weight[e] = d
        edges_d.append(d)
    # list for the nodes degree
    deg_d = []
    for v in G.vertices():
      d = v.out_degree()
      deg_d.append(d)

    # list for the distances to the sun
    sun_d = []
    for j in range(1,amount_stars):
      d = dist1(0,j,stars)
      sun_d.append(d)

    # all to all stars distance histogram
    all_to_all_d = []
    for i in range(amount_stars):
      for j in range(amount_stars):
        d = shortest_distance(G,i,j,weights=edge_weight)
        if d != 'inf':
          all_to_all_d.append(d)
    # calculate the global clustering coeficient
    global_clus = global_clustering(G,weight=edge_weight)
    print("Global clustering: ", global_clus)

    # list for the local clustering
    local_list = []
    local_clus= local_clustering(G,weight=edge_weight)
    for v in local_clus:
      local_list.append(v)

    # ploting all histograms
    plt.hist(edges_d,bins=10)
    plt.title("Histograma para la distancia de arcos")
    plt.xlabel("Distancia (pc)")
    plt.ylabel("Cantidad de arcos")
    plt.show()
    plt.hist(sun_d,bins=10)
    plt.title("Histograma para la distancia al Sol")
    plt.xlabel("Distancia (pc)")
    plt.ylabel("Cantidad de estrellas")
    plt.show()
    plt.hist(deg_d,bins=10)
    plt.title("Histograma para la distribucion de grados")
    plt.xlabel("grado")
    plt.ylabel("Cantidad de estrellas")
    plt.show()
    plt.hist(all_to_all_d,bins=10)
    plt.title("Histograma para la distancia de entre todas las estrellas")
    plt.xlabel("Distancia (pc)")
    plt.ylabel("Cantidad de estrellas")
    plt.show()
    plt.hist(local_list,bins=10)
    plt.title("Histograma para el coeficiente de agrupamiento local")
    plt.xlabel("Coeficiente")
    plt.ylabel("Cantidad de estrellas")
    plt.show()

def calculate_s_path(stars, amount_stars, solution_net, tras_index):
    '''
    calcualtes the shortest path between the sun and trasppist-1
    Arguments:
        stars: dictionaries list with the stars data
        amount_stars: number of the total stars using in our model
        solution_net: edges list taking from the model solution
        tras_index: the trasppist-1 position in stars
    Returns: 
        the shortest path value and the edges in this path
    '''
    #creates the graph
    G1 = Graph(directed=False)
    G1.add_vertex(amount_stars)
    edge_weight = G1.new_edge_property("double")

    # addding the edges
    for (i,j) in solution_net:
        e = G1.add_edge(i,j)
        d = dist1(i,j,stars)
        edge_weight[e] = d
    #calculates the shortest path
    v_list, e_list = shortest_path(G1,0,tras_index,weights=edge_weight)
    dist = shortest_distance(G1,0,tras_index,weights=edge_weight)
    s_path = [str(e).replace('(', '').replace(')', '').split(',') for e in e_list]
    s_path = [(int(i), int(j)) for i,j in s_path]

    return dist, s_path

def add_tras(stars,amount_stars):
    '''
    adds trasppist-1 data to the stars list
    Arguments:
        stars: dictionaries list with the stars data
        amount_stars: number of the total stars using in our model
    Returns:
        the trasppist-1 position in the list and the list with the new data added
    '''
    #trasppist data
    tras = {'x': 11.7265, 'y': -2.8239, 'z': 1.0623, 'prop': 'Trappist-1'}
    stars_fil = stars[0:amount_stars-1]
    stars_fil.append(tras)

    return amount_stars-1,stars_fil