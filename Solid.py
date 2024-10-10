import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd



# Reloading the module
import importlib

import Mesh
importlib.reload(Mesh)

import Utility
importlib.reload(Utility)

import Elasticity
importlib.reload(Elasticity)


class Solid():
    
    def __init__(self, name):
        
        self.name = name
        
        # Number of unknown displacement components
        self.num_unknown = 3
        
        # Extract info from the mesh file of Abacus format (*.inp)
        # and write the info to csv files.
        Mesh.set_Mesh(name)
        
        # Create nodes and elements
        self.set_Mesh()
        
        
    # Elastic properties for future use    
    def set_ElasticProperty(self, E=1, nu=0.3):
        self.Eprop = Elasticity.Elasticity(E, nu)  
        
        
    # Node and element info
    def set_Mesh(self):
        
        
        ###########################
        # Read node info
        ###########################
        # Line format in 'nodes.csv': node label, x, y, z
        df = pd.read_csv('./mesh/nodes.csv', header=None)
        # The first column include node labels, not required here.
        self.nodes = df.to_numpy()[:, 1:] 
        # Number of nodes
        self.num_nodes = self.nodes.shape[0]
        
        
        ###########################
        # Read element info
        ###########################
        # Line in 'elements.csv': element label, node labels
        df = pd.read_csv('./mesh/elements.csv', header=None)
        # The first column include element labels. No need to store them
        self.elements = df.to_numpy()[:, 1:] 
        # The original node labels start from 1
        # Reduce the node labels by 1 according to python convention (starting from 0)
        self.elements = self.elements - 1 
        # Number of elements
        self.num_elements = self.elements.shape[0]
        
        
        
        ##########################################
        # Read node labels on the Dirichlet boundary
        ##########################################
        # 'tag_node_BC_Dirichlet.csv' contains a single line, containing the node labels
        df = pd.read_csv('./mesh/nodes_BC_Dirichlet.csv', header=None)
        # The dataframe has a single line, whose last entry is ' ', which should be dropped
        # 'to_numpy' gives a single-row 2D array, which needs to be converted to 1D, by using [0]
        self.nodes_BC_Dirichlet = df.iloc[:, :-1].to_numpy()[0] 
        # Reduce the node labels by 1 according to python convention (starting from 0)
        self.nodes_BC_Dirichlet -= 1
        
        
        ##########################################
        # Read node labels on the Neumann boundary
        ##########################################
        filepath = './mesh/nodes_BC_Neumann.csv'
        if Utility.is_nonzero_file(filepath):
            # 'tag_node_BC_Dirichlet.csv' contains a single line, containing the node labels
            df = pd.read_csv('./mesh/nodes_BC_Neumann.csv', header=None)
            # The dataframe has a single line, whose last entry is ' ', which should be dropped
            # 'to_numpy' gives a single-row 2D array, which needs to be converted to 1D, by using [0]
            self.nodes_BC_Neumann = df.iloc[:, :-1].to_numpy()[0] 
            # Reduce the node labels by 1 according to python convention (starting from 0)
            self.nodes_BC_Neumann -= 1
        
        
        
        '''
        ##########################################
        # Read facet labels on the whole boundary
        ##########################################
        # Read the file only when it exists and is not empty
        filepath = './mesh/facets_BC.csv'
        if Utility.is_nonzero_file(filepath):
            # Line in 'facet_BC.csv': facet label, node labels
            df = pd.read_csv(filepath, header=None)
            # The first column include facet labels. No need to store them
            self.facets_BC = df.to_numpy()[:, 1:] 
            # The original node labels start from 1
            # Reduce the node labels by 1 according to python convention (starting from 0)
            self.facets_BC -= 1
            # Number of facets
            # self.num_facets_BC = self.facets_BC.shape[0]
            
        else:
            self.facets_BC = np.array([])
        
        
        
        
        
        ##########################################
        # Read facet labels on the Neumann boundary
        ##########################################
        # Read the file only when it exists and is not empty
        filepath = './mesh/facets_BC_Neumann.csv'
        if Utility.is_nonzero_file(filepath):
            # 'tag_node_BC_Dirichlet.csv' contains a single line, containing the node labels
            df = pd.read_csv(filepath, header=None)
            # The dataframe has a single line, whose last entry is ' ', which should be dropped
            # 'to_numpy' gives a single-row 2D array, which needs to be converted to 1D, by using [0]
            self.facets_BC_Neumann = df.iloc[:, :-1].to_numpy()[0] 
            # Reduce the facet labels by 1 according to python convention (starting from 0)
            self.facets_BC_Neumann -= 1
            
        else:
            self.facets_BC_Neumann = np.array([])
        
        

        '''

    
        
        
        
    