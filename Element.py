import numpy as np

class Tetrahedron():
    
    def __init__(self, nodes):
        # Corner node coordinates
        # nodes: (4,3), nodes[i,j] is the ith corner's jth coordinate.
        self.nodes = nodes
        self.cal_J_V_Grad()
        
        
        
        
    
    
    def cal_J_V_Grad(self):
        
        # Rescaling
        r = self.nodes
        self.L = np.min((np.linalg.norm(r[0]-r[1]),
                    np.linalg.norm(r[0]-r[2]),
                    np.linalg.norm(r[0]-r[3]),
                    np.linalg.norm(r[1]-r[2]),
                    np.linalg.norm(r[1]-r[3]),
                    np.linalg.norm(r[2]-r[3]),))
        
        self.nodes /= self.L
        r = self.nodes
        
        
        # Computation of the rescaled element
        x01, y01, z01 = r[1] - r[0]
        x02, y02, z02 = r[2] - r[0]
        x03, y03, z03 = r[3] - r[0]
        
        # det(Jabobian matrix)
        self.J = x01*y02*z03 - x01*y03*z02 - x02*y01*z03 + x02*y03*z01 + x03*y01*z02 - x03*y02*z01
        # Volume
        self.V = np.abs(self.J) / 6
        
        # Gradient
        # Each row contains the gradient of an interpolating function wrt x,y,z
        self.grad = np.zeros((4, 3))
        J = self.J
        
        self.grad[0,0]= (-y01*z02 + y01*z03 + y02*z01 - y02*z03 - y03*z01 + y03*z02)/J
        self.grad[0,1]= (x01*z02 - x01*z03 - x02*z01 + x02*z03 + x03*z01 - x03*z02)/J
        self.grad[0,2]= (-x01*y02 + x01*y03 + x02*y01 - x02*y03 - x03*y01 + x03*y02)/J
        self.grad[1,0]= (y02*z03 - y03*z02)/J
        self.grad[1,1]= (-x02*z03 + x03*z02)/J
        self.grad[1,2]= (x02*y03 - x03*y02)/J
        self.grad[2,0]= (-y01*z03 + y03*z01)/J
        self.grad[2,1]= (x01*z03 - x03*z01)/J
        self.grad[2,2]= (-x01*y03 + x03*y01)/J
        self.grad[3,0]= (y01*z02 - y02*z01)/J
        self.grad[3,1]= (-x01*z02 + x02*z01)/J
        self.grad[3,2]= (x01*y02 - x02*y01)/J
        
        
        
    def cal_M(self):
        self.M = np.array([[2,1,1,1],[1,2,1,1],[1,1,2,1],[1,1,1,2]]) * self.V / 20 * self.L**3
        
        
        
    
    def cal_K(self):
        # Compute stiffness matrices
        # K[a,b,i,j]: a,b=0,1,2, representing x,y,z 
        # i,j=0,1,2,3, corner labels
        K = np.zeros((3,3,4,4))
        
        for a in range(3):
            for b in range(3):
                for i in range(4):
                    for j in range(4):
                        K[a, b, i, j] = self.grad[i, a] * self.grad[j, b]
    
    
        self.K = K * self.V * self.L
        
    
    
        
        
        
        