[#Network Generation Problem]
[Mixed Integer Linear Progromming (MILP)]
# Interstellar Network Topology Design and Protocol Evaluation

The recent discovery of potentially habitable planets orbiting the TRAPPIST-1 system intensified interest in interstellar exploration. Although the most distant human-made object (Voyager 1) could take thousands of years to reach the nearest star (Proxima Centauri), numerous projects have proposed propulsion ideas (and energy systems) to reduce this duration. in orders of magnitude. Although the event renewed public interest in interstellar exploration, the fundamental problem behind interstellar missions remains, of course, distance, which imposes significant challenges. In particular, communication protocols must deal with unprecedented delays in signal propagation.

This work is based on the hypothesis that the immense interstellar distances argue for an autonomous multi-hop relay strategy rather than a direct point-to-point link. Multi-hop network transmission of information occurs every day on the Internet, but interplanetary and near-Earth space environments certainly impose different communication challenges. Furthermore, thinking about them on an interstellar scale would demand a networked system of unprecedented scale delays and latency.
Delay Tolerant Networking (DTN) technologies are resumed. Based on these, network topology design mechanisms will be identified that allow interstellar missions to be provided with connectivity with minimal latency and resources.

For the generation of networks, which optimizing certain characteristics, the use of Mixed Integer Linear Programming (MILP) is proposed as a possible solution. Optimization theory and graph theory makes it possible to describe models that allow, through specification with certain equations, to obtain a solution network that complies exactly with what the equations dictate.
We will study two models, adapting them in some specifications, such as in the case in which we want edges with a value that does not exceed a maximum distance. On the other hand, we will study the possibility of creating networks with heuristic methods such as the Erdös model, the Watts-Strogatz model, as well as a completely random edge selection method. In this analysis we will compare all the models, and we will also show the advantages and disadvantages of each of them.

This repository contains the algorithms used in my thesis

## Requeriments

* Python (version 3.10.12)
* gurobipy (version 11.0.0 more information [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python))
* matplotlib (version 3.8.2)
* numpy (version 1.21.5)
* graph-tool (version 2.57 installation instructions [here](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions#native-installation))

## Structure 

* DCOM.py : This contains the model DCOM, expressed in gurobi syntax.
* SPDOM.py : This file contains the model SPDOM, expressed in gurobi syntax.
* random_gen.py: This file contains a random algorithm to make an intelestellar net.
* erdos_gen.py: This file contains the algorithm based in the Erdös model to generate a net.
* strogatz_gen.py : This file contains the algorithm based in the Watts-Strogatz model.
* utilities.py : This file contains functions that we need to make some metrics to test the algorithms and models.
* base_final.csv: a lite version of a stars database avalaible [here](http://www.astronexus.com/hyg).

## Examples for test them all

You can use the function to makes loop over the `amount_stars` var and plot the result. For exmple:

     solution_net,elapsed_time,mem_usage = model_SPDOM(stars_fil,amount_stars,min_degree,
                                                    max_degree,min_diam,max_diam,min_cpl,max_cpl,
                                                    edge_lenght,max_lenght)

The `min_degree` `max_degree` `min_diam` `max_diam` `min_cpl` `max_cpl` are the vars uses to test the behavoir of the model (in these example SPDOM) when the boolean var `max_lenght` is **False**. you can loop over `amount_stars` var and save the data in several list and then plot the information.
When the `max_lenght` is **True** you also can test how `edge_lenght` affect the net.
