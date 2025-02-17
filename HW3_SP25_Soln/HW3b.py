import math
import numericalMethods as nm
#endregion

#region functions
def FZ(m,uL):
    """
    Computes the part inside the integral of the t-Probability distribution.
    :param args: (m, uL) the degrees of freedom and upper limit of integration
    :return: computed value of integrand
    """

    k_m=km(m)

    def tPDF(args):
        x,m=args[0],args[1] # from Simpson, args is 4 elements long.  Just need first two.
        base = 1+(x**2)/m
        expnt =  -(m+1)/2
        return base**expnt
        # return math.pow(base,epnt)
        # the following calls the Simpson integration function with required args tuple
    I=nm.Simpson(tPDF,0,uL,1000, m)
    return k_m*I

def gamma(alpha):
    """
        This computes the gamma function for a positive m.  If m is an integer, we simply compute
    the factorial of m-1.  If m is a float, we use Simpson integration to compute the gamma function.
    :param m: the value for which to compute the gamma function
    :return: the computed value of the gamma function
    """
    if alpha.is_integer():
        return math.factorial(int(alpha)-1)
    else:
    def fn(t,a):
        return math.exp(-t) * math.pow(t, a-1)

gam = math.gamma(alpha)
    if alpha.is_integer():
        return math.factorial(int(alpha) - 1)
    else:
    def fn(args):
        t,a=args[0],args[1]
        return math.exp(-t)*math.pow(t,a-1)

    if(alpha % 1 == 0): # we have an integer and can compute an exact answer
            g = 1
        for i in range(1, int(alpha)):
            g *= i
        gg = nm.Simpson(fn, (alpha, 0, 0, 50), 100000)  # these numbers give values similar to Table A2
        return g
            g=nm.Simpson(fn,(alpha,0,0,50),100000) # these numbers give values similar to Table A2
                return g

def km(m):
    """
            compute K_m for t-distribution
    :param m: degrees of freedom
    :return: k_m
    """
    k_m = gamma(0.5*m+0.5)/(math.sqrt(m*math.pi)*gamma(0.5*m))
    return k_m

def main():
    """
    This computes the area below the t-distribution (i.e., probability density function) given a value for degrees of freedom and an
    upper limit of integration.  I use the Simpson method to integrate numerically and a lower limit of -100*abs(u)
    """
g=gamma(1.78)
Fz=FZ((7,1.89))
getOut=False
while (getOut is False):
        m = input("Degrees of freedom (integer): ").strip().lower()
        u = input("Upper integration limit (float):").strip().lower()
        m=int(m)
        u=float(u)
        Fz=FZ((m,u))
        print("F({:0.3f})={:0.3f}".format(u,Fz))
        getOut=input("Go Again (Y/N)?").strip().lower().__contains__("n")

#endregion

if __name__ == '__main__':
    main()