from DoolittleMethod import Doolittle
from matrixOperations import checkMatrixSoln
import copy

def is_symmetric(A):
    """check symmetry"""
    n=len(A)
    for i in range(n):
        for j in range(n):
            if A[i][j] !=A[j][i]:
                return False
    return True

def is_positive_definite(A):
    """check if matrix A is positive definite"""
    n=len(A)
    for i in range(n):
        sub_matrix=[[A[row][col]for col in range(i+1)]for row in range (i+1)]
        if determinant(sub_matrix) <=0:
            return False
    return True

def determinant(A):
    """calculate determinant using Laplace"""
    n=len(A)
    if n ==1:
        return A[0][0]
    if n==2:
        return A[0][0] * A[1][1]-A[0][1]*A[1][0]

    det=0
    for col in range(n):
        minor =[[A[row][c]for c in range(n)if c !=col]for row in range(1,n)]
        det += ((-1) ** col) * A[0][col] * determinant(minor)
    return det

def cholesky_decomposition(A):
    """A=L*L^T"""
    n=len(A)
    L= [[0]*n for _ in  range(n)]

    for i in range(n):
        for j in range(i+1):
            sum_val=sum(L[i][k]* L[j][k] for k in range(j))
            if i==j:
                L[i][j] = (A[i][i] - sum_val) **0.5
            else:
                L[i][j] = (A[i][j]- sum_val) / L[j][j]
    return L

def forward_substitution(L, b):
    """Ly=b"""
    n=len(L)
    y=[0]*n
    for i in range(n):
        y[i]= (b[i]- sum(L[i][j] *y[j] for j in range(i))) /(L[i][i])
    return y

def backward_substitution(U,y):
    """Ux=y"""
    n=len(U)
    x=[0]*n
    for i in range (n-1, -1,-1):
        x[i]=(y[i]- sum(U[i][j]* x[j] for j in range(i+1,n)))/ (U[i][i])
    return x

def solve_cholesky(A, b):
    """solve Ax=b using choleskey decomposition"""
    L= cholesky_decomposition(A)
    LT= [[L[j][i]for j in range(len(L))]for i in range(len(L))]
    y= forward_substitution(L,b)
    x= backward_substitution(LT, y)
    return x

def solve_system(A,b):
    """determine which method to solve Ax=b"""
    if is_symmetric(A) and is_positive_definite(A):
        method="Cholesky"
        x=solve_cholesky(A,b)
    else:
        method= "Doolittle"
        augmented_matrix= [A[i] + [b[i]]for i in range (len(A))]
        x, _ = Doolittle(augmented_matrix)
    return x, method

def main():
    """finally solving the systems"""
    A1= [
        [1,-1,3,2],
         [-1,5,-5,-2],
         [3,-5,-5,-2],
         [2,-2,3,21]
         ]
    b1= [15,-35,94,1]

    A2= [
        [4,2,4,0],
        [2,2,3,2],
        [4,3,6,3],
        [0,2,3,9]
        ]
    b2= [20,36,60,122]

    #system #1:
    x1, method1=solve_system(A1,b1)
    print("\nSystem #1:", [round(i, 3) for i in x1])
    print("Method used:", method1)

    #system #2:
    x2,method2=solve_system(A2,b2)
    print("\nSystem #2:", [round(i, 3) for i in x2])
    print("Method used:", method2)

    #verify solutions:
    print("\nVerify solutions:")
    print("System#1:", checkMatrixSoln(A1, x1, augmented=False))
    print("System#2:", checkMatrixSoln(A2, x2, augmented=False))

if __name__ == "__main__":
    main()