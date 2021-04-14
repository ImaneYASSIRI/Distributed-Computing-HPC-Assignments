@types('float[:]','float[:]', 'int', 'int','float', 'float' , 'int')
def solve_1d_diff_pyccel(u, un, nt, nx, dt, dx, nu):
    
    #fill the update of u
    for n in range(nt):
        for i in range(nx):
          un[i] = u[i]
        for i in range(1,nx): 
          u[i] = un[i] - nu * dt / (dx**2) * (u[i] - 2*un[i] + un[i - 1])
    
    return 0
