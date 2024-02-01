# Interstellar Network Topology Design and Protocol Evaluation

This repository contains the algorithms used in my thesis

## Requeriments

* Python (version 3.10.12)
* gurobipy (version 11.0.0 more information [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python))
* matplotlib (version 3.8.2)
* numpy (version 1.21.5)
* graph-tool (version 2.57 installation instructions [here](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions#native-installation))

## Structure 

* NetGen1.py : This contains the model DCOM, expressed in gurobi syntax.
* NetGen2.py : This file contains the model SPDOM, expressed in gurobi syntax.
* random_gen.py: This file contains a random algorithm to make an intelestellar net.
* erdos_gen.py: This file contains the algorithm based in the Erdös model to generate a net.
* strogatz_gen.py : This file contains the algorithm based in the Watts-Strogatz model.
* utilities.py : This file contains functions that we need to make some metrics to test the algorithms and models.
* base_final.csv: a lite version of a stars database avalaible [here](http://www.astronexus.com/hyg).

