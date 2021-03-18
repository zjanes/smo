# Plan here is to write a (very small) sim library to show some basic epi sims. Hopefully we'll accompany this with some animation library stuff to get some pretty videos, and then figure out how to put them on a squarespace website.

# The general plan for this model is that is proceeds in paired steps. At each time step, first attempt infections, then perform movement.

# The infection step consists of a sweep of all currently infected nodes in order of current location, calculating their current neighbours and then attepmting to infect any it should

# The movement step is a general sweep of the whole population, updating all states (except infection) and moving all nodes.

# So we maintain three data structures for the population, at least for now (might grow).

# - A dict of all nodes
# - A list of all infected nodes
# - A matrix (nested dicts) off all coordinates, so that we can quickly find neighbours


import random, operator, timeit as t
from functools import reduce

class Strain_parameters:
    """No idea what this does right now. Will have to look into it.""

    def __init__(self, transmissibility, durations):

        self.transmissibility = transmissibility 
        self.durations = durations # A map keyed on sate code giving the duration of that state

class Node:

    def __init__(self, strains, coords, UID):

        self.UID = UID
        self.coords = {'x': coords['x'], 'y': coords['y']}
        self.states = {s: "S" for s in strains}
        self.state_progresses = {}

    def infect(self, strain, parameters, infected_list):
        """Attempt infection of this node"""

        if self.states[strain] == "S" and random.random() < parameters.transmissibility:

            self.states[strain] = "I"
            self.state_progresses[strain] = parameters.durations["I"] # state durations work by decrementing each time step
            infected_list.append(self)

    def index_infect(self, strain, parameters, infected_list):

        self.states[strain] = "I"
        self.state_progresses[strain] = parameters.durations["I"] # state durations work by decrementing each time step
        infected_list.append(self)

    def recover(self, strain, disease_parameters):

        print("Recovery function not yet defined")

    def update_states(self, disease_parameters, next_state, infected_list):
        "For each strain present in self.state_progresses, decrements the state counter. If the counter is 0, it looks up the next state using next_state, then looks up that state's duration in the disease_parameters"

        for s in self.state_progresses: 

            if self.state_progresses[s] == 0:
                
                if self.states[s] == "I":

                    infected_list.remove(self)

                new_state = next_state(self.states[s])
                self.states[s] = new_state
                self.state_progresses[s] = disease_parameters[s].durations[new_state]

            else:

                self.state_progresses[s] = self.state_progresses[s] - 1
        
class Population:

    def __init__(self, dimension, population_size, strains, disease_parameters, disease_model):
        """Initialises a population of susceptible hosts in random positions on an integer-grid plane (hosts can occupy the same space).
        Does not initialise infections."""

        self.nodes = {i: Node(strains, {'x': random.randint(0, dimension-1), 'y': random.randint(0, dimension-1)}, i) for i in range(population_size)} # Set up the population
        self.matrix = {x : {y : {} for y in range(dimension)} for x in range(dimension)} # Initialise the matrix
        self.infected = []
        self.disease_parameters = disease_parameters

        self.dimension = dimension
        self.next_state_function = disease_model
        self.global_time = 0

        for n in self.nodes.values():

            self.matrix[n.coords['x']][n.coords['y']][n.UID] = n            

def sir_model():
    "Returns a function that can be called to determine the next state"

    def model(state):

        if state == "S": return "I"
        elif state == "I": return "R"
        elif state == "R": return "R"

    return model

def si_model():

    def model(state):

        if state == "S": return "I"
        elif state == "I": return "I"

    return model
        
def get_neighbours(population, coords):
    #print("z")
    neighbour_coords = [(x,y)
                        for x in get_valid_coords(population.dimension, coords['x'])
                        for y in get_valid_coords(population.dimension, coords['y'])]
    neighbour_coords.remove((coords['x'], coords['y']))

    #print(neighbour_coords)
    
    return reduce(operator.concat, [list(population.matrix[x][y].values()) for (x,y) in neighbour_coords])

def get_valid_coords(dimension, i):

    if i == 0: return [i+1, i]

    elif i == dimension-1: return [i-1, i]

    else: return [i-1, i, i+1]

def update_and_move_all(population):
    """For each node, updates their state change time counters, and moves them in the grid."""
    
    for node in population.nodes.values():
    
        node.update_states(population.disease_parameters, population.next_state_function, population.infected)
        move(population.matrix, population.dimension, node)
        
        
def move(matrix, dimension, node):
    """Moves a node to a new location in the position matrix"""
    
    new_coords = {'x': random.randint(0, dimension-1), 'y': random.randint(0, dimension-1)}

    if (new_coords['x'] != node.coords['x'] or new_coords['y'] != node.coords['y']):
        
        matrix[new_coords['x']][new_coords['y']][node.UID] = node
        matrix[node.coords['x']][node.coords['y']].pop(node.UID)
        node.coords = new_coords

    
def infect_nodes(population, strain):

    [try_infect(h, strain, get_neighbours(population, h.coords), population.disease_parameters[strain], population.infected) for h in population.infected]
    
def try_infect(host, strain, neighbours, parameters, infected):

    targets = [n for n in neighbours if n.states[strain] == "S"] # Find susceptible neighbours

    [t.infect(strain, parameters, infected) for t in targets]

def record(population):

    return(len(population.infected))
    
def run_sim():

    population = Population(1000, 10000,
                            ["1"], {"1": Strain_parameters(0.2, {"I": 10, "R": 10})},
                            sir_model())
    index_host = random.choice(population.nodes)
    index_host.index_infect("1", population.disease_parameters["1"], population.infected)
    strain = "1"
    
    while population.infected !=  []:

        infect_nodes(population, strain)
        update_and_move_all(population)
        #print(record(population))

    print("Simulation finished. No more infectious hosts remain in the population.")

import datetime

before = datetime.datetime.now()
run_sim()
after = datetime.datetime.now()

print(after - before)
    

# Current state 3.04.2020 sheesh

# The simulation works (I think) with SI and SIR models. It's not been tested to see whether it conforms to expectd values, but it outputs what appears so far to be sensible data. To be certain I need to implement a tracker for the number of recovered.

# The next thing to do could be one of two. I think the direction I want to take this little sim in is to make it into a metapopulation model with migration, and then use that to simulate sampled data, so that I can use basic epidemiological tools to infer things like the R values.

# So I can implement the metapopulation version of the model. That means a change in the way coordinates are handled. Instead of holding the whole coordinate grid as a nested dict I only hold the currently active nodes. That will be far more sparse, but necessitates a custom data structure. Then a node belongs to a certain population defined by a single coordinate and dimensions from the coordinate. Most migration then happens within the population but migration is also allowed between populations.

#
