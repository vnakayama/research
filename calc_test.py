from functions import *

cmap = cm.jet #stores old healpy colloring scheme (cmap='jet')
cmap.set_under('w')

#.....................................
# Calculation of norms, vectors and angles:
#.....................................

nmax = 4
Lmax = 4    # Size of the fundamental domain: Lx = Ly = Lmax and Lz.
Lz   = 0.1

partial_time = time.time()

normvec_dict = vector_n(nmax,Lmax,Lz)   # <- {sqrt(((L/lx)*nx)^2 + ((L/ly)*ny)^2 + ((L/lz)*nz)^2), [(L/lx)*nx, (L/ly)*ny, (L/lz)*nz]}
print(len(norms(nmax,Lmax,Lz)))
normang_dict = vectors_to_angles(normvec_dict) #initialize norms and angles
print(len(normang_dict))


#.....................................
# Tranforming the dict into a matrix:
#.....................................

n = len(normang_dict.keys())

print(list(normang_dict.keys())[0])
#print(normang_dict[list(normang_dict.keys())[0]][1][0])

norm_ang = []
for i in range(n):
  for j in range(len(normang_dict[list(normang_dict.keys())[i]])):
    norm_ang.append((list(normang_dict.keys())[i], 
                     normang_dict[list(normang_dict.keys())[i]][j][0], 
                     normang_dict[list(normang_dict.keys())[i]][j][1]))

norm_ang = np.array(norm_ang)
print(len(norm_ang))

print('Finished calculating norms, vectors and angles in {time} s'.format(time = time.time() - partial_time))

# theta = norm_ang[:,1] -> latitude
# phi   = norm_ang[:,2] -> longitude

#%%

#.....................................
# Defining parameters:
#.....................................

lmax = 35 #int(2.*64.-1.) # Maximum multipole
Nside = 64

l, m  = hp.Alm.getlm(lmax)
index = hp.Alm.getidx(lmax, l, m)

A = 1 # np.sqrt(2./np.pi)
R = 1 # Mpc ==>> R = R_H = \chi_rec = 45.7 x 10^9 anos luz (1 pc = 3.26 anos-luz) => R = 14.0184 Mpc
k = 2.0*np.pi*norm_ang[:,0]/(Lmax*R) # norm of Fourier vector
sqrt_Pspec = np.sqrt(A*(k**(-3.0)))
Cte = (1./3.) * np.sqrt(2./np.pi)

#####


kR = k * R

alm = np.complex128(np.zeros(len(index)))

np.random.seed(seed=654321)#(idx**2)*round(2.*Lz))
rand_num1 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num2 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num3 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num4 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
        
for idx in index:
    if l[idx] > 1:
        jl_kR = special.spherical_jn(l[idx], kR)
        Ylm1  = special.sph_harm( np.float64(m[idx]), np.float64(l[idx]), norm_ang[:,2], norm_ang[:,1]) #( m, l, phi, theta)
        #print(m[idx], '/', l[idx], '/', '/', Ylm1[0])
        psi1 = (rand_num1 + ((-1.)**l[idx])*np.conjugate(rand_num1))
        psi2 = (rand_num2 + ((-1.)**l[idx])*np.conjugate(rand_num2))
        psi3 = (rand_num3 + ((-1.)**l[idx])*np.conjugate(rand_num3))
        psi4 = (rand_num4 + ((-1.)**l[idx])*np.conjugate(rand_num4))
        sum1 = ( ((psi1 + ((-1.)**m[idx])*psi3) * np.conjugate(Ylm1)) 
              +  ((psi2 + ((-1.)**m[idx])*psi4) * Ylm1) )
        sum2 = sqrt_Pspec * jl_kR
        Total_sum = sum2 * sum1
        alm[idx] = ((1j)**l[idx] * Cte) * np.sum( Total_sum ) # para um par l,m

mapa = hp.alm2map(alm, Nside)
#hp.mollview(mapa, title='Slab, Lz=0.1, Lx=Ly=4, nmax=4', cmap=cmap, rot=[0.,210.,70.])
#plt.savefig('Slab_SWonly_seed654321_Lz0p1_LxLy4_nmax4.png', dpi=100)
hp.mollview(mapa, title='Torus, L = 0.1, nmax=4', cmap=cmap, rot=[0.,210.,70.])
plt.savefig('Torus_Teste.png', dpi=100)
