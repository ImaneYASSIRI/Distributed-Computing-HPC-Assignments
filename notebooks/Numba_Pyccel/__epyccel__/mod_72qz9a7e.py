@types('float[:,:]','float[:,:]','float[:,:]','float[:,:]', 'int', 'float','float', 'float' , 'int')
def solve_2d_nonlinearconv_pyccel(u, un, v, vn, nt, dt, dx, dy, c):

    import numpy as np
    ###Assign initial conditions
    ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
    ##set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
    v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
    row, col = u.shape
    u0=np.ones((row, col))
    v0=np.ones((row, col))
    for i in range(row):
        for j in range(col):
            u0[i][j] = u[i][j] 
            v0[i][j] = v[i][j]
    
    #fill the update of u and v
    for n in range(nt):
        for i in range(row):
            for j in range(col):
                  un[i][j] = u[i][j] 
                  vn[i][j] = v[i][j]    
        for i in range(1,row):
            for j in range(1,col):
                u[i][j] = un[i][j] -  u0[i][j]* dt / dx * (un[i][j] - un[i - 1][j]) - vn[i][j] * dt / dy * (un[i][j] - un[i][j-1])
                v[i][j] = vn[i][j] -  v0[i][j]* dt / dx * (vn[i][j] - vn[i - 1][j]) - vn[i][j] * dt / dy * (vn[i][j] - vn[i][j-1])
               
        
    return 0
