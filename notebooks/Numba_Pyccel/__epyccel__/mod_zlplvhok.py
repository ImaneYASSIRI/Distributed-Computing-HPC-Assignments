@types('float[:]','float[:]', 'int', 'int','float', 'float')
def solve_1d_nonlinearconv_pyccel(u, un, nt, nx, dt, dx):

    #fill the update of u
    for n in range(nt):
        for i in range(nx):
          un[i] = u[i]
        for i in range(1,nx): 
          u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1])
    
    return 0
