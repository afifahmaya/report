import numpy as np

#Function Definition
#Function to create Right Hand Side and Left Hand Side Matrix
def const_mat(N):
    h = 1/N
    A = np.zeros((N,N))
    b = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            A[i][i] = 4
            b[i][j] = h*h*(np.sin(np.pi*i/N)*np.sin(np.pi*j/N))
            #b[i][j] = 1
            if (i==j-1 and i%4==3) or (i==j+1 and i%4==0):
                A[i][j] = 0
            elif i==j-1:
                A[i][j] = -1
            elif i==j+1:
                A[i][j] = -1
            elif i==j-4:
                A[i][j] = -1
            elif i==j+4:
                A[i][j] = -1
    return A,b

#Jacobi
def jacobi(b):
    v = np.zeros_like(b)
    N = len(b)
    for k in range(1000):
        v_old = np.copy(v)
        for i in range(N-2):
            for j in range(N-2):
                v[i+1][j+1] = 0.25 * (v_old[i][j+1] + v_old[i+2][j+1] + v_old[i+1][j] + v_old[i+1][j+2] + b[i+1][j+1])
        #print(np.max(v))
        if np.max(np.abs(v - v_old)) < 1e-4:
            break
    return(k)

#Gauss-Seidel
def gauss_seidel(b):
    v = np.zeros_like(b)
    N = len(b)
    for k in range(1000):
        v_old = np.copy(v)
        for i in range(N-2):
            for j in range(N-2):
                v[i+1][j+1] = 0.25 * (v[i][j+1] + v[i+2][j+1] + v[i+1][j] + v[i+1][j+2] + b[i+1][j+1])
        #print(np.max(v))
        if np.max(np.abs(v - v_old)) < 1e-4:
            break
    return(k)

#SOR
def sor(omega,b):
    v = np.zeros_like(b)
    N = len(b)
    for k in range(1000):
        v_old = np.copy(v)
        for i in range(N-2):
            for j in range(N-2):
                v[i+1][j+1] = (1-omega)*v_old[i+1][j+1] + omega*(0.25 * (v[i][j+1] + v[i+2][j+1] + v[i+1][j] + v[i+1][j+2] + b[i+1][j+1]))
        #print(np.max(v))
        if np.max(np.abs(v - v_old)) < 1e-4:
            break
    return(k)

#Calculate Optimum Omega for SOR
def calc_op_omega(A):
    L = np.diag(np.diag(A,-1),-1)
    U = np.diag(np.diag(A,1),1)
    T = 1/np.diag(A)*(L+U)
    rho = np.max(np.linalg.eigvals(T))
    opt_omega = 2/(1+np.sqrt(1-rho*rho))
    
    return opt_omega

#main
Ns = [10,20,30,40,50]
om = 1.5
for n in Ns:
    matA, matB = const_mat(n)
    iterj = jacobi(matB)
    iterg = gauss_seidel(matB)
    iters = sor(om,matB)
    op_om = calc_op_omega(matA)
    print("N =",n)
    print("Minimum Iteration Jacobi =",iterj)
    print("Minimum Iteration Gauss-Seidel =",iterg)
    print("Minimum Iteration SOR using omega(1.5) =",iters)
    print("Optimal omega for SOR =",op_om)