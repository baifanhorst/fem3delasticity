import numpy as np

# Reloading the module
import importlib

import Element
importlib.reload(Element)

def assembly(solid):
    # Generate matrices K and M
    
    nn = solid.num_nodes
    # K^{ij}, i,j=0,1,2, representing x,y,z
    K = np.zeros((3, 3, nn, nn))
    M = np.zeros((nn, nn))
    
    # Assembly of K and M
    
    element_nodes = np.zeros((4,3))
    
    for idx_element, e in enumerate(solid.elements):
        
        # Get element corner coordinates
        for k, label_node in enumerate(e):
            element_nodes[k] = solid.nodes[label_node]
        
        # Create an element object    
        ele = Element.Tetrahedron(element_nodes)
        
        
        # Assembly
        for ie in range(4):
            i = e[ie]
            for je in range(4):
                j = e[je]
                M[i,j] += ele.M[ie, je]
                K[:,:,i,j] += ele.K[:,:,ie,je]
                
                #for a in range(3):
                #    for b in range(3):
                #        K[a,b,i,j] += ele.K[a,b,ie,je]
                        
    return K, M


def set_BC_Dirichlet(LHS, RHS, solid, func_u):
    # Applying Dirichlet BC
    # Always apply the Dirichlet BC at the end!!!
    for idx_unknown in range(solid.num_unknown):
        for idx_node in solid.nodes_BC_Dirichlet:
            # Compute the boundary node value
            x,y,z = solid.nodes[idx_node]
            val_Dirichlet = func_u[idx_unknown](x,y,z)
            
            # Global index
            idx_var = idx_node + idx_unknown * solid.num_nodes
    
            # Move the known term to RHS
            RHS -= LHS[:, idx_var] * val_Dirichlet
            LHS[:, idx_var] = 0
            # Reset the row for the boundary node 
            LHS[idx_var] = 0
            LHS[idx_var, idx_var] = 1
            RHS[idx_var] = val_Dirichlet
            
            
def set_BC_Neumann(RHS, solid, func_p):
    
    
    for idx_e, e in enumerate(solid.elements):
        
        # Initialize an array to store element corner coordinates
        r = np.zeros((4,3))
        
        for ie in range(4):
            # ie: element node label, from 0 to 3
            # Global node labels in cyclic manner
            # e.g. e = [0,1,2,3], then np.roll(e,-1) = [1,2,3,0]
            eroll = np.roll(e, -ie)
        
            # If the first 3 entries in ids_e are labels of a Neumann facet 
            if np.all( np.isin(eroll[0:-1], solid.nodes_BC_Neumann) ):
                
                
                
                
                r[0] = solid.nodes[eroll[0]]
                r[1] = solid.nodes[eroll[1]]
                r[2] = solid.nodes[eroll[2]]
                r[3] = solid.nodes[eroll[3]]
                
                #print(r)
            
                r01 = r[1] - r[0]
                r02 = r[2] - r[0]
                cross = np.cross(r01, r02)
                # The coefficient in front of the infinitesimal element dxi deta
                coeff_area = np.linalg.norm( cross )
                # Unit norm vector
                vect_norm = cross / coeff_area
                
                # print('vect_norm', vect_norm)
        
        
                r03 = r[3] - r[0]
        
                # if the norm vector points into the solid, revert it.
                if np.dot(r03, vect_norm) > 0:
                    vect_norm = - vect_norm
                    
                
                
                nx, ny, nz = vect_norm
            
                
            
                
                # For all three p components
                for idx_p in range(3):
                
                    # Get the three corner p values
                    p = np.zeros(3)
                    for i in range(3):
                        x,y,z = r[i]
                        p[i] = func_p[idx_p](x, y, z, nx, ny, nz, solid.Eprop.G, solid.Eprop.lamb)
                    
                    
                    idx_var = eroll[0] + idx_p * solid.num_nodes
                    RHS[idx_var] += coeff_area/24 * (2*p[0] + p[1] + p[2])
                    idx_var = eroll[1] + idx_p * solid.num_nodes
                    RHS[idx_var] += coeff_area/24 * (2*p[1] + p[2] + p[0])
                    idx_var = eroll[2] + idx_p * solid.num_nodes
                    RHS[idx_var] += coeff_area/24 * (2*p[2] + p[0] + p[1])
                    
                
                    

       
            
            
                




def set_LinearSystem(solid, func_u, func_f, func_p):
    
    nn = solid.num_nodes
    # Matrices assembly
    LHS = np.zeros((solid.num_unknown * nn, solid.num_unknown * nn))
    RHS = np.zeros((solid.num_unknown * nn))
    
    # LHS  
    G = solid.Eprop.G
    l = solid.Eprop.lamb
    
    
    K, M = assembly(solid)
    
    K_diag_sum = K[0,0] + K[1,1] + K[2,2]
        
    for a in range(3):
        for b in range(3):
            if a == b:
                LHS[a*nn : (a+1)*nn, b*nn : (b+1)*nn] = (G+l) * K[a,b] + G * K_diag_sum
            else:
                LHS[a*nn : (a+1)*nn, b*nn : (b+1)*nn] = l * K[a,b] + G * K[b,a]
            
            
      
    # Body force   
    for idx_direction in range(3):
        f = np.zeros(nn)
        for idx_node in range(nn):
            x,y,z = solid.nodes[idx_node]
            f[idx_node] = func_f[idx_direction](x, y, z, G, l)
        
        RHS[idx_direction*nn : (idx_direction+1)*nn] = np.dot(M, f)
           
            
    # Applying Neumann BC
    
    # If a Neumann boundary exists
    if hasattr(solid, "nodes_BC_Neumann"):
        set_BC_Neumann(RHS, solid, func_p)
        
    # Applying Dirichlet BC
    set_BC_Dirichlet(LHS, RHS, solid, func_u)
    
    
    return LHS, RHS
            
            
  
def solv_LinearSystem(LHS, RHS, solid):
    nn = solid.num_nodes
    # Solve the linear system
    u = np.linalg.solve(LHS, RHS)
    print('residue', np.max(np.abs(RHS - np.dot(LHS, u))))
    solid.u = np.zeros((solid.num_unknown, nn))
    for a in range(solid.num_unknown):
        solid.u[a] = u[a*nn : (a+1)*nn]


def cal_Exact_Solution(solid, func_u):
    nn = solid.num_nodes
    solid.u_exact = np.zeros((solid.num_unknown, nn))
    x = solid.nodes[:,0]
    y = solid.nodes[:,1]
    z = solid.nodes[:,2]
    for a in range(solid.num_unknown):
        solid.u_exact[a, :] = func_u[a](x, y, z)
        
        
            
##############################################
# Test codes
###############################################
def set_BC_Neumann_test(RHS, solid, func_p):
    
    count_facet_Neumann = 0
    for idx_e, e in enumerate(solid.elements):
        
        # Initialize an array to store element corner coordinates
        r = np.zeros((4,3))
        
        for ie in range(4):
            # ie: element node label, from 0 to 3
            # Global node labels in cyclic manner
            # e.g. e = [0,1,2,3], then np.roll(e,-1) = [1,2,3,0]
            eroll = np.roll(e, -ie)
        
            # If the first 3 entries in ids_e are labels of a Neumann facet 
            if np.all( np.isin(eroll[0:-1], solid.nodes_BC_Neumann) ):
                
                print('\n')
                
                count_facet_Neumann += 1
                print('idx_e', idx_e)
                print('e', e)
                print("eroll", eroll[0:-1])
                
                
                r[0] = solid.nodes[eroll[0]]
                r[1] = solid.nodes[eroll[1]]
                r[2] = solid.nodes[eroll[2]]
                r[3] = solid.nodes[eroll[3]]
                
                #print(r)
            
                r01 = r[1] - r[0]
                r02 = r[2] - r[0]
                cross = np.cross(r01, r02)
                # The coefficient in front of the infinitesimal element dxi deta
                coeff_area = np.linalg.norm( cross )
                # Unit norm vector
                vect_norm = cross / coeff_area
                
                # print('vect_norm', vect_norm)
        
        
                r03 = r[3] - r[0]
        
                # if the norm vector points into the solid, revert it.
                if np.dot(r03, vect_norm) > 0:
                    vect_norm = - vect_norm
                    
                print('vect_norm', vect_norm)
                
                nx, ny, nz = vect_norm
            
                
            
                
                # For all three p components
                for idx_p in range(3):
                
                    # Get the three corner p values
                    p = np.zeros(3)
                    for i in range(3):
                        x,y,z = r[i]
                        p[i] = func_p[idx_p](x, y, z, nx, ny, nz, solid.Eprop.G, solid.Eprop.lamb)
                    
                    print(f'p{idx_p}', p)
                    
                    p0, p1, p2 = p
                    
                    idx_var = eroll[0] + idx_p * solid.num_nodes
                    print('idx_var', idx_var)
                    #RHS[idx_var] += coeff_area/24 * (2*p0 + p1 + p2)
                    idx_var = eroll[1] + idx_p * solid.num_nodes
                    print('idx_var', idx_var)
                    #RHS[idx_var] += coeff_area/24 * (2*p1 + p2 + p0)
                    idx_var = eroll[2] + idx_p * solid.num_nodes
                    print('idx_var', idx_var)
                    #RHS[idx_var] += coeff_area/24 * (2*p2 + p0 + p1)
                    
                
                    

    print("count_facet_Neumann", count_facet_Neumann)           
            
            
            
            
############################################
# Old codes
############################################
def set_LinearSystem_Elasticity_old(solid, func_f, func_u):
    
    nn = solid.num_nodes
    
    # Matrices assembly
    LHS = np.zeros((solid.num_unknown * nn, solid.num_unknown * nn))
    RHS = np.zeros((solid.num_unknown * nn))
    K = np.zeros((3, 3, nn, nn))
    M = np.zeros((nn, nn))
    
    
    # Assembly of K and M
    
    element_nodes = np.zeros((4,3))
    
    for idx_element, e in enumerate(solid.elements):
        
        # Get element corner coordinates
        for k, label_node in enumerate(e):
            element_nodes[k] = solid.nodes[label_node]
        
        # Create an element object    
        ele = Element.Tetrahedron(element_nodes)
        
        
        # Assembly
        for ie in range(4):
            i = e[ie]
            for je in range(4):
                j = e[je]
                M[i,j] += ele.M[ie, je]
                
                for a in range(3):
                    for b in range(3):
                        K[a,b,i,j] += ele.K[a,b,ie,je]
        
    
    
    # LHS  
    G = solid.Eprop.G
    l = solid.Eprop.lamb
    
    K_diag_sum = K[0,0] + K[1,1] + K[2,2]
        
    for a in range(3):
        for b in range(3):
            if a == b:
                LHS[a*nn : (a+1)*nn, a*nn : (a+1)*nn] = (G+l) * K[a,a] + G * K_diag_sum
            else:
                LHS[a*nn : (a+1)*nn, b*nn : (b+1)*nn] = (G+l) * K[a,b]
            
            
      
    # RHS   
    for a in range(3):
        f = np.zeros(solid.num_nodes)        
        for idx_node in range(solid.num_nodes):
            r = solid.nodes[idx_node]
            f[idx_node] = func_f[a](r[0], r[1], r[2], G, l)
            
        RHS[a*nn : (a+1)*nn] = np.dot(M, f)
        
    
    # Applying Dirichlet BC
    for a in range(3):
        for idx_node in solid.tag_nodes_BC_Dirichlet:
            # Compute the boundary node value
            r = solid.nodes[idx_node]
            val_Dirichlet = func_u[a](r[0], r[1], r[2])
            
            idx_var = idx_node + a * nn
    
            # Move the known term to RHS
            RHS -= LHS[:, idx_var] * val_Dirichlet
            LHS[:, idx_var] = 0
            # Reset the row for the boundary node 
            LHS[idx_var] = 0
            LHS[idx_var, idx_var] = 1
            RHS[idx_var] = val_Dirichlet
            
            
    # Applying Neumann BC
    

        
    # Solve the linear system
    u = np.linalg.solve(LHS, RHS)
    print('residue', np.max(np.abs(RHS - np.dot(LHS, u))))
    solid.u = np.zeros((solid.num_unknown, nn))
    for a in range(3):
        solid.u[a] = u[a*nn : (a+1)*nn]
        
    
    
   
    
                

                
                
    
                
                
                

                
    
                
        
            
        
        