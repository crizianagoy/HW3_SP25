#region imports
from copy import deepcopy as dcpy
from math import cos,pi
import numericalMethods as nm
import matrixOperations as mo
#endregion

#region Functions
def LUFactorization(A):
    """
    This is the Lower-Upper factorization part of Doolittle's method.  The factorizaiton follows the work in
    Kreyszig section 20.2.  Note: L is the lower triangular matrix with 1's on the diagonal.  U is the upper traingular matrix.
    :param A: a nxn matrix
    :return: a tuple with (L, U)
    """
    n = len(A)

    # Step 1: Initialize U and L matrices
    U = [[0 if c != r else a for c, a in enumerate(A[r])] for r in range(n)]
    L = [[1 if c == r else 0 for c in range(n)] for r in range(n)]

    for k in range(n):  # Ensuring k is properly scoped
        # Check for zeros and perform partial pivoting
        if abs(U[k][k]) < 1e-12:  # Pivot if diagonal is too small
            for swap_row in range(k + 1, n):
                if abs(U[swap_row][k]) > abs(U[k][k]):
                    U[k], U[swap_row] = U[swap_row], U[k]  # Swap rows in U
                    L[k], L[swap_row] = L[swap_row], L[k]  # Swap rows in L
                    break  # Stop after the first swap

        print("U matrix before error:", U)
        print(f"U[{k}][{k}] =", U[k][k])

        if U[k][k] == 0:
            raise ValueError(f"Zero detected in U[{k}][{k}], LU factorization failed!")

        for i in range(k + 1, n):
            L[i][k] = A[i][k] / U[k][k]  # Compute L values
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]  # Modify U values

    return L, U

def BackSolve(A,b,UT=True):
    """
    This is a backsolving algorithm for a matrix and b vector where A is triangular
    :param A: A triangularized matrix (Upper or Lower)
    :param b: the right hand side of a matrix equation Ax=b
    :param UT: boolean of upper triangular (True) or lower triangular (False)
    :return: the solution vector x, from Ax=b
    """
    nRows = len(b)
    b = [row if isinstance(row, (int, float)) else row[0] for row in b]
    x = [0] * nRows  # Initialize x as a list of zeros

    if UT:
        for nR in range(nRows - 1, -1, -1):  # Backward substitution
            s = 0
            for nC in range(nR + 1, nRows):
                s += A[nR][nC] * x[nC]
            x[nR] = (b[nR] - s) / A[nR][nR]  # Directly use b[nR]
    else:
        for nR in range(nRows):  # Forward substitution
            s = 0
            for nC in range(nR):
                s += A[nR][nC] * x[nC]
            x[nR] = (b[nR] - s) / A[nR][nR]  # Directly use b[nR]

    return x

def Doolittle(Aaug):
    """
    The Doolittle method for solving the matrix equation [A][x]=[b] is:
    Step 1:  Factor [A]=[L][U]
    Step 2:  Solve [L][y]=[b] for [y]
    Step 3:  Solve [U][x]=[y] for [x]
    :param Aaug: the augmented matrix
    :return: the solution vector x
    """
    A,b=mo.separateAugmented(Aaug)
    L,U=LUFactorization(A)
    B=mo.MatrixMultiply(L,U)
    y=BackSolve(L,b, UT=False)
    x=BackSolve(U,y, UT=True)
    return x, "Doolittle"

def main():
    A=[[3, 5, 2],[0,8,2],[6,2,8]]
    L,U=LUFactorization(A)
    print("L:")
    for r in L:
        print(r)

    print("\nU:")
    for r in U:
        print(r)

    aug=[[3,9,6,4.6],[18,48,39,27.2], [9,-27,42,9]]
    aug = [[3, 1, -1, 2],
          [1, 4, 1, 12],
          [2, 1, 2, 10]]
    x=Doolittle(aug)
    x=[round(y,3) for y in x]
    print("x: ", x)
    y=nm.GaussSeidel(aug,[0,0,0])
    y=[round(z,3) for z in y]
    b=mo.checkMatrixSoln(aug,y)
    print("b: ",b)
#endregion

if __name__ == "__main__":
    main()