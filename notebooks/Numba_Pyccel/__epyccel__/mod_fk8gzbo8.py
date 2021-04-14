@types('float[:,:]','float[:,:]', 'int', 'float','float', 'float' , 'int')
def solve_2d_diff_pyccel(u, un, nt, dt, dx, dy, nu):
    row, col = u.shape
    
    ##Assign initial conditions
    #set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
    
    
    #fill the update of u and v
    for n in range(nt):
        for i in range(row):
            for j in range(col):
                  un[i][j] = u[i][j]
        for i in range(1,row - 1):
            for j in range(1,col - 1):
                u[i][j] = un[i][j] + nu * dt / dx**2 * (un[i+1][j] - 2*un[i][j] + un[i - 1][j]) + nu * dt / dy**2 * (un[i][j+1] - 2*un[i][j] + un[i][j-1])
        
        
    return 0
