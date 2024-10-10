import numpy as np

# Reloading the module
import importlib

import Element
importlib.reload(Element)

def cal_Stress(solid):
    
    
    G = solid.Eprop.G
    l = solid.Eprop.lamb
    
    
    nn = solid.num_nodes
    
    # Stress arrays
    # stress[0],..., stress[5]: stress_x, stress_y, stress_z, stress_xy, stress_xz, stress_yz
    # stress[6]: stress_v
    stress = np.zeros((7, nn))
    # element stress array
    strain_element = np.zeros(6)
    stress_element = np.zeros(6)
    
    # Number of neighboring elements for each node
    count_element = np.zeros(nn)
    
    
    # Array storing element node coordinates
    element_nodes = np.zeros((4,3))
    
    
    
    for e in solid.elements:
        
        # Reset
        strain_element[:] = 0.0
        
        # Get element corner coordinates
        for k, label_node in enumerate(e):
            element_nodes[k] = solid.nodes[label_node]
        
        # Create an element object    
        ele = Element.Tetrahedron(element_nodes)
        
        for ie in range(4):
            
            # Get nodal displacement
            idx_node_global = e[ie]
            u = solid.u[:, idx_node_global]
        
            # stain_x
            strain_element[0] += u[0] * ele.grad[ie, 0]
            # stress_y
            strain_element[1] += u[1] * ele.grad[ie, 1]
            # stress_z
            strain_element[2] += u[2] * ele.grad[ie, 2]
            
            # stress_xy
            strain_element[3] += 0.5 * (u[0] * ele.grad[ie, 1] + u[1] * ele.grad[ie, 0])
            # stress_xz
            strain_element[4] += 0.5 * (u[0] * ele.grad[ie, 2] + u[2] * ele.grad[ie, 0])
            # stress_yz
            strain_element[5] += 0.5 * (u[1] * ele.grad[ie, 2] + u[2] * ele.grad[ie, 1])
            
        # Recall the rescaling in the element object
        strain_element /= ele.L
            
        # Volume change ratio
        vol_rati0 = strain_element[0] + strain_element[1] + strain_element[2]
            
        # Stress
        stress_element = 2*G * strain_element
        stress_element[0:2+1] += l * vol_rati0
        
        # Nodal contributions:
        for ie in range(4):

            # Get nodal displacement
            idx_node_global = e[ie]
            # Add the contribution
            stress[0:-1, idx_node_global] += stress_element
            # Count the element
            count_element[idx_node_global] += 1
            
    # Average
    stress[0:-1] = stress[0:-1] / count_element
        
    # Von Mises
    stress[6] = cal_VonMises(stress[0:-1])
        
    solid.stress = stress
        
        
        
def cal_VonMises(stress):
    stress_v = (stress[0] - stress[1])**2
    stress_v += (stress[1] - stress[2])**2
    stress_v += (stress[0] - stress[2])**2
    stress_v += 6 * (stress[3]**2 + stress[4]**2 + stress[5]**2)
    stress_v = np.sqrt( stress_v * 0.5 )
    return stress_v
        
            
            
            
            