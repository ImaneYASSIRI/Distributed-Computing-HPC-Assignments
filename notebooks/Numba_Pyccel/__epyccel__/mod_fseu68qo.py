@types('float[:,:]', 'float[:,:]','float[:,:]', 'int', 'int', 'int','float', 'float')
def solve_2d_poisson_pyccel(p, pd, b, nx, ny, nt, dx, dy):
    row, col = p.shape
    # Source
    b[int(ny / 4), int(nx / 4)]  = 100
    b[int(3 * ny / 4), int(3 * nx / 4)] = -100
    
    
    for it in range(nt):
        for i in range(nx): pd[i] = p[i]
        for j in range(2, row):
            for i in range(2, col):
                p[j-1, i-1] = (((pd[j-1, i] + pd[j-1, i-2]) * dy**2 +
                                (pd[j, i-1] + pd[j-2, i-1]) * dx**2 -
                                b[j-1, i-1] * dx**2 * dy**2) / 
                                (2 * (dx**2 + dy**2)))
        p[0, :] = 0
        p[ny-1, :] = 0
        p[:, 0] = 0
        p[:, nx-1] = 0
        
    return 0
