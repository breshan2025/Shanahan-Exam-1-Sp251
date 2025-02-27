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
    U = [[0 for _ in range(n)] for _ in range(n)]
    L = [[0 for _ in range(n)] for _ in range(n)]
    for j in range(n):  # j is row index for U
        # (a)
        for k in range(j, n):  # k is column index for U, complete row j elements
            U[j][k] = A[j][k]  #
            for s in range(j):  # s is column index for L row j, and row index for U column k
                U[j][k] -= L[j][s] * U[s][k]
        # (b)
        for i in range(j, n):  # j is column index for L, i is row index for L
            if i == j:
                L[i][j] = 1  # diagonal elements
            else:
                sig = 0
                for s in range(j):  # i is row index for L, s is column index for L, s is row index for U, j is column index for U
                    sig += L[i][s] * U[s][j]
                L[i][j] = (1 / (U[j][j])) * (A[i][j] - sig)
    return (L, U)

def BackSolve(A,b,UT=True):
    """
    This is a backsolving algorithm for a matrix and b vector where A is triangular
    :param A: A triangularized matrix (Upper or Lower)
    :param b: the right hand side of a matrix equation Ax=b
    :param UT: boolean of upper triangular (True) or lower triangular (False)
    :return: the solution vector x, from Ax=b
    """
    b=mo.makeColumnVector(b)
    nRows=len(b)
    nCols=nRows
    x=[[0] for r in range(nRows)]
    if UT:
        for nR in range(nRows-1,-1,-1):
            s=0
            for nC in range(nR+1,nRows):
                s+=A[nR][nC]*x[nC][0]
            x[nR][0]=1/A[nR][nR]*(b[nR][0]-s)
    else:
        for nR in range(nRows):
            s=0
            for nC in range(nR):
                s+=A[nR][nC]*x[nC][0]
            x[nR][0]=1/A[nR][nR]*(b[nR][0]-s)
    B = mo.checkMatrixSoln(A, x, False)
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
    return x  #x should be a column vector of form [[], [], [], ..., []]

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