@types('float[:]','float[:]', 'int', 'int','float', 'float' , 'int')
def solve_1d_linearconv_pyccel(u,un ,nt, nx, dt, dx, c):
    for n in range(nt):
      for i in range(nx):
          un[i] = u[i]
      for i in range(1,nx): 
          u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])
    return 0
