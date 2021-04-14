@types('float[:,:]','float[:,:]','float[:,:]','float[:,:]', 'int', 'float','float', 'float' , 'float')
def solve_2d_burger_pyccel(u, un, v, vn, nt, dt, dx, dy, nu):
    ###Assign initial conditions
    ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 
    ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    v[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
    row, col = u.shape
    
    #fill the update of u and v
    for n in range(nt):
        for i in range(row):
            for j in range(col):
                  un[i][j] = u[i][j]
                 
                  vn[i][j] = v[i][j]
        for i in range(1,row - 1):
            for j in range(1,col - 1):
                u[i][j] = un[i][j] - un[i][j] * dt / dx * (un[i][j] - un[i - 1][j]) - vn[i][j] * dt / dy * (un[i][j] - un[i][j-1]) + nu * dt / dx**2 * (un[i+1][j] - 2*un[i][j] + un[i - 1][j]) + nu * dt / dy**2 * (un[i][j+1] - 2*un[i][j] + un[i][j-1])
                v[i][j] = vn[i][j] - un[i][j] * dt / dx * (vn[i][j] - vn[i - 1][j]) - vn[i][j] * dt / dy * (vn[i][j] - vn[i][j-1]) + nu * dt / dx**2 * (vn[i+1][j] - 2*vn[i][j] + vn[i - 1][j]) + nu * dt / dy**2 * (vn[i][j+1] - 2*vn[i][j] + vn[i][j-1])
         
    return 0
