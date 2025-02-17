#region imports
from numericalMethods import GPDF, Probability, Secant
#endregion

#region function definitions
def main():
    """
    This program is designed to solicit input from the CLI and return a probability from the Gaussian Normal Distribution
    by integrating the probability density function using Simpson's 1/3 rule for numerical integration.

    The user may decide upon a two-sided integration or one-sided integration.
    If one-sided, we integrate from the smaller of mean-5*stDev or c-stDev.
    If two-sided, we integrate from -c to c.
    Here is my step-by-step plan:
    1. Decide mean, stDev, and c and if I want P(x>c) or P(x<c) or P(-c<x<c) or P(-c>x>c).
    2. Define args tuple and c to be passed to Probability
    3. Pass args, and a callback function (GPDF) to Probability
    4. In probability, pass along GPDF to Simpson along with the appropriate args tuple
    5. Return the required probability from Probability and print to screen.
    :return: Nothing to return, just print results to screen.
    """
    #region testing user input
    # setting the initial default values
    Again = True
    mean = 0
    stDev = 1.0
    c = 0.5
    P = 0.5
    seek = 0 # 0 for seeking P, 1 for seeking c
    OneSided = True  # integrates from mu-5*sig if true, from mu-(c-mu) to mu+(c-mu) if False
    GT = False
    yesOptions = ["y","yes","true"]
    while Again==True:
        # The following code solicites user input through the CLI.
        response = input(f"Population mean? ({mean:0.3f})").strip().lower()  # strip off leading or trailing spaces and make lower case.
        mean = float(response) if response != '' else mean  # ternary operator based on response

        response = input(f"Standard deviation? ({stDev:0.3f})").strip().lower()
        stDev = float(response) if response != '' else stDev

        response = input(f"Specify c or P? ({'c' if seek == 1 else 'P'})").strip().lower()
        seek = 0 if response == 'p' else seek
        seek = 1 if response == 'c' else seek

        if seek==1:  # looking for P given c
            response = input(f"c value? ({c:0.3f})").strip().lower()
            c = float(response) if response != '' else c
        else:  # looking for c given P
            response = input(f"P value? ({P:0.3f})").strip().lower()
            P = float(response) if response != '' else P

        response=input(f"Probability greater than c? ({GT})").strip().lower()
        GT = True if response in yesOptions else False

        response=input(f"One sided? ({OneSided})").strip().lower()
        OneSided = True if response in yesOptions else False
        # handle computation of P given c
        if seek == 1:
            if OneSided==True:
                prob = Probability(GPDF,(mean,stDev),c,GT=GT)
                print(f"P(x"+(">" if GT == True else "<") + f"{c:0.2f}" +"|"+f"{mean:0.2f}"+", "+f"{stDev:0.2f}" +f") = {prob:0.2f}")
            else:
                prob = Probability(GPDF, (mean, stDev),c, GT=True)
                prob = 1-2*prob
                if GT == True:
                    print(f"P({mean-(c-mean)}>x>{mean+(c-mean)}|{mean:0.2f},{stDev:0.2f}) = {1-prob:0.3f}")
                else:
                    print(f"P({mean-(c-mean)}<x<{mean+(c-mean)}|{mean:0.2f},{stDev:0.2f}) = {prob:0.3f}")
        # handle computation of c given P
        else:
            def fn(c):
                if OneSided == True:
                    return Probability(GPDF, (mean, stDev), c, GT=GT) - P
                else:
                    # compute probability between mean-(c-mean) and mean+(c-mean)
                    prob = 1-2*Probability(GPDF, (mean, stDev), c, GT=True)
                    return (1-prob)-P if GT == True else prob-P
            c,N = Secant(fn,mean,mean+stDev,maxiter=50, xtol=1e-5)
            if OneSided == True:
                print(f"P(x" + (">" if GT == True else "<") + f"{c:0.2f}" + "|" + f"{mean:0.2f}" + ", " + f"{stDev:0.2f}" + f") = {P:0.2f}")
            else:
                if GT == True:
                    print(f"P({mean - (c - mean):0.3f}>x>{mean + (c - mean):0.3f}|{mean:0.2f},{stDev:0.2f}) = {P:0.3f}")
                else:
                    print(f"P({mean - (c - mean):0.3f}<x<{mean + (c - mean):0.3f}|{mean:0.2f},{stDev:0.2f}) = {P:0.3f}")
        response = input(f"Go again? (Y/N)").strip().lower()
        Again = True if response in ["y","yes","true"] else False
    #endregion

#endregion

if __name__ == "__main__":
    main()