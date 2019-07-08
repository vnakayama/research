from context import * 

class Topology:

    def __init__(self,cutoff,fundamental_domain,L_z):
        self.NMax = cutoff
        self.L = fundamental_domain
        self.Lz = L_z
        self.L2 =float("%.3f" % (self.L/self.Lz))

    def Norms(self):
        Norm = []
        for Nx in np.arange(0,self.NMax+1,1):
            for Ny in np.arange(0,self.NMax+1,1):
                for Nz in range(1,self.NMax+1):
                    Norm.append((Nx)**2 + (Ny)**2 + (self.L2*Nz)**2)
        Norm = list(set(Norm)) # remove duplicates with set() 
        Norm.pop(0) #remove zero norm
        return sorted(Norm)
        
    def Vector_N(self):
        vec_n = {}
        for Norm in self.Norms: 
            Nx_Count = math.sqrt(Norm)      
            Vec_N[Nx_Count] = list()
            for Nx in np.arange(0,Nx_Count+2,1): #we do not include components at and bellow the nz=0 plane
                Ny2_Max = Norm - (Norm)**2
                if Ny2_Max > float(0):
                    Ny_Count = math.sqrt(Ny2_Max)
                    for Ny in np.arange(0,Ny_Count+2,1):
                        if Norm - (Nx)**2 - (Ny)**2 > float(0):
                            L2Nz = math.sqrt(Norm - (Nx)**2 - (Ny)**2)
                            if L2Nz !=0 and float(format(L2Nz/self.L2, '0.8f')) == round(L2Nz/self.L2) and (Nx)**2 + (Ny)**2 + (L2Nz)**2 == Norm:                          
                                Vec_N[math.sqrt(Norm)].append([Nx,Ny,L2Nz])
        return Vec_N

    def Vectors_To_Angles(self):
        New_N = dict()
        for k in self.Vector_N.keys(): # <- for in 'norm'
            Angles = []
            for Vector in self.Vector_N[k]:
                # calculate the azimuthal angle
                if Vector[0] != 0: # <- (L/lx)*nx
                    phi = math.atan(float(Vector[1])/float(Vector[0])) 
                    # <- arctan( (L/ly)*ny / (L/lx)*nx) = tan( (Lx*ny) / (Ly*nx) )
                else:
                    phi = math.pi/2  # tan^-1(inf) = pi/2
                    if Vector[1] == 0:  # tan^-1(0/0) = 0.0
                        phi = float(0)
                # calculate the polar angle
                if Vector[2] != 0: # <- (L/lz)*nz
                    # sqrt returns float, float/int = float
                    theta = math.atan((math.sqrt(Vector[0]**2 + Vector[1]**2))/Vector[2])
                else:
                    theta = math.pi/2  # tan^-1(inf) = pi/2
                Angles.append([theta, phi]) # <- theta_k, phi_k
            New_N.update({k: Angles})
        return New_N
    
class Torus(Topology):
    pass

class Slab(Topology):
    pass