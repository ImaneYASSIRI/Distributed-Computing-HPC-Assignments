@types('float[:,:]','float[:,:]', 'int', 'float','float', 'float' , 'int')
def solve_2d_linearconv_pyccel(u, un, nt, dt, dx, dy, c):
   
    row, col = u.shape

    #fill the update of u and v
    for n in range(nt):
        for i in range(row):
            for j in range(col):
                  un[i][j] = u[i][j]
        for i in range(1,row):
            for j in range(1,col):
                u[i][j] = un[i][j] - c * dt / dx * (un[i][j] - un[i - 1][j]) - c * dt / dy * (un[i][j] - un[i][j-1])
               
               
    return 0
