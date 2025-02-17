#region imports
import copy
import matrixOperations as mo
import numericalMethods as nm
import DoolittleMethod
import math

import Gauss_Seidel as GS
import random
#endregion

#region functions
def Cholesky(Aaug):
    """
    This function finds the solution to a matrix equation Ax=b by the Colesky method
    :param Aaug: An augmented matrix
    :return: the solution vector x, L and Ltrans as a tuple
    """
    #step 1:  split the Aaug into A and b (see separateAugmented in Gauss_Seidel.py)
    A,b=GS.separateAugmented(Aaug)
    #step 2:  factor into L and Ltrans according to Cholesky formulae
    n=len(A)
    L = [[0 for c in range(n)] for r in range(n)]
    L[0][0] = math.sqrt(A[0][0])
    for j in range(1,n):
        L[j][0] = A[j][0]/L[0][0]
    for j in range(1,n):
        sm = 0
        for s in range(j):
            sm += math.pow(L[j][s],2)
        L[j][j] = math.sqrt(A[j][j]-sm)

        for p in range(j+1,n):
            sm = 0
            for s in range(j):
                sm += L[j][s]*L[p][s]
            L[p][j] = 1/L[j][j]*(A[p][j]-sm)
    L = [[round(num,3) for num in r] for r in L]
    LT = mo.Transpose(L)
    AA = GS.matrixMult(L,LT)
    pass
    #step 3:  use backsolving to find x (see methods in DoolittleMethod.py)
        #[L][LT][x]=[b] -> [L][y]=[b] ->[LT][x]=[y]
    y = DoolittleMethod.BackSolve(L,b,UT=False)
    x = DoolittleMethod.BackSolve(LT,y,True)
    bb=GS.matrixMult(A,x)
    #step 4:  return (x,L,Ltrans)
    return (x,L,LT)
    pass

def SymPosDef(A):
    """
    This function first finds the transpose of A and then compares all elements of A to Atrans.
    If I pass that test, I create a vector x with random numbers and perform xtrans*A*x to see if>0.
    :param A: a nxn matrix
    :return: True if symmetric, positive definite
    """
    #step 1:  recall that a transpose has elements such that Atrans[i][j] = A[j][i] see page 267 in MAE3013 text
    AT=mo.Transpose(A)
    #step 2:  check that all elements of A and Atrans are the same. if fail->return false
    sym=Symetric(A,AT)
    if sym == False: return False
    #step 3:  produce a vector x of length n filled with random floats between -1 and +1
    for i in range(5):  #I'll test this 5 times with random numbers between -1 and 1
        x=[[2.0*(random.random()-0.5)] for n in range(len(A))]
        # step 4:  compute xtrans*A*x
        xT = mo.Transpose(mo.makeColumnVector(x))
        Ax=GS.matrixMult(A,x)
        xTAx=GS.matrixMult(xT,Ax)
        if xTAx[0][0]<=0: return False
    #step 5:  if step 4 > 0 return true else return false
    return True
def Symetric(A,AT):
    """
    This function takes two nxn matrices and compares them element-by-element to see if they are the same
    :A: a square matrix
    :AT: transpose of a square matrix
    :return: True or False
    """
    return True if A==AT else False

def main():
    """
    Step 1:  I need to first define the matrices given in part a) of HW3_2024.
    Step 2:  pass a matrix to SymPosDef to tell if it is symmetric, positive definite
    Step 3:  based on result of Step 2, use either the Doolittle or Cholesky method to solve
    Steo 4:  check my answer by multiplying A*x to see if I get b
    Step 5:  print the solution vector and which method was used to the cli
    returns:  nothing
    """
    #region step 1:
    Aaug1 = [[1,-1,3,2,15],[-1,5,-5,-2,-35],[3,-5,19,3,94],[2,-2,3,21,1]]
    Aaug2 = [[4,2,4,0,20],[2,2,3,2,36],[4,3,6,3,60],[0,2,3,9,122]]

    #region make an asymmetric matrix
    # I'm making a copy of Aau1, but with rows 0 and 2 switched to make it non-symmetric
    Aaug3 = copy.deepcopy(Aaug1)
    #swap rows 0 & 2 to make it intentionally asymmetric
    Aaug3[0]=copy.deepcopy(Aaug1[2])
    Aaug3[2]=copy.deepcopy(Aaug1[0])
    #endregion
    #endregion

    #region step 2:
    isSymPos_A1 = SymPosDef(GS.separateAugmented(Aaug1)[0])
    isSymPos_A2 = SymPosDef(GS.separateAugmented(Aaug2)[0])
    isSymPos_A3 = SymPosDef(GS.separateAugmented(Aaug3)[0])
    #endregion

    #region step 3
    if isSymPos_A1:
        x1,L1,LT1=Cholesky(Aaug1)
        method1= 'Cholesky'
    else:
        x1=DoolittleMethod.Doolittle(Aaug1)
        method1 = 'Doolittle'
    if isSymPos_A2:
        x2,L2,LT2=Cholesky(Aaug2)
        method2= 'Cholesky'
    else:
        x2=DoolittleMethod.Doolittle(Aaug2)
        x2-round(x2,3)
        method2 = 'Doolittle'
    if isSymPos_A3:
        x3,L3,LT3=Cholesky(Aaug3)
        method3= 'Cholesky'
    else:
        x3=DoolittleMethod.Doolittle(Aaug3)
        method3 = 'Doolittle'
    #endregion

    #region step 4
    A,b=GS.separateAugmented(Aaug1)
    sv=GS.matrixMult(A,x1)
    pass
    A,b=GS.separateAugmented(Aaug2)
    sv=GS.matrixMult(A,x2)
    pass
    #endregion

    #region step 5
    print("x1 = ", x1, ", method = {:}".format(method1))
    print("x2 = ", x2, ", method = {:}".format(method2))
    print("x3 = ", [[round(e[0],3)] for e in x3], ", method = {:}".format(method3))
    #endregion

#endregion

if __name__ == "__main__":
    main()