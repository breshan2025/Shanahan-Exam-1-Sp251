#region imports
import math
import random
#end region

#region definitions
def ind_variables():
    ''' generates independent random variables'''
    a1 = random.random()
    a2 = random.random()

    b1 = math.sqrt(-2 * math.log(a1)) * math.cos(2 * math.pi * a2)
    b2 = math.sqrt(-2 * math.log(a1)) * math.sin(2 * math.pi * a2)
    return b1, b2

def log_normal(mu, sig):
    '''Standard log -normal density function(PDF) normalized over (0, infinity)
        This is used to model the gravel's diameter'''
    # narrow to only one variable
    b1, _ = ind_variables()
    return math.exp(mu + sig * b1)

def lognormal_cumulative(x, mu, sig):
    ''' Calculates the cumulative distribution function of the log normal distribution'''
    if x <= 0:
        return 0
    return  .5 +.5 * math.erf((math.log(x) - mu) / (sig * math.sqrt(2)))

def truncated_log(mu, sig, Dmin, Dmax,N):
    ''' generates N random samples from a truncated log normal distribution
        used to establish the minimum and maximum diameters that will fit through the screens
        then generate random samples'''

    # establish sample set
    samples = []
    # set parameters for comparison
    while len(samples) < N:
        sample = log_normal(mu, sig)
        if Dmin <= sample <= Dmax:
            samples.append(sample)
    return samples


def calculate(sample):
        n = len(sample)
        mean = sum(sample) / n
        variance = sum((x - mean) ** 2 for x in sample) / (n - 1)
        return mean, variance

def gravel_production(mu, sig, Dmin, Dmax, N, n):
    ''' Generates random samples from the truncated log normal distribution'''
    samples = []
    for _ in range(N):
        sample = truncated_log(mu, sig, Dmin, Dmax,n)
        samples.append(sample)

    sample_means = []
    sample_variances = []

    for sample in samples:
        mean, variance = calculate(sample)
        sample_means.append(mean)
        sample_variances.append(variance)

    total_mean = sum(sample_means) / N
    total_variance = sum((mean - total_mean) ** 2 for mean in sample_means) / N
#end region

#region print
    ''' Print results'''
    print(f"Sample Statistics:")
    for i, sample in enumerate(samples, 1):
        mean, variance = calculate(sample)
        print(f"Sample {i}: Mean = {mean:.3f}, Variance = {variance:.3f}")

    print(f"\nSampling Mean Statistics:")
    print(f"Mean of Sampling Means: {total_mean:.3f}")
    print(f"Variance of Sampling Means: {total_variance:.3f}")


def main():
    '''  obtain input from user then run the gravel production simulation'''

    print("Gravel Production Simulation\n")

    # Default values
    default_mu = 0
    default_sig = 1
    default_Dmin = 0.1
    default_Dmax = 10.0
    N = 11
    n = 100

    # User inputs with defaults
    mu = float(input(f"Enter mean of ln(D) (default = {default_mu}):") or default_mu)
    sig = float(input(f"Enter standard deviation of ln(D)(default = {default_sig}):") or default_sig)
    Dmin = float(input(f"Enter minimum diameter (Dmin)(default = {default_Dmin}):") or default_Dmin)
    Dmax = float(input(f"Enter maximum diameter (Dmax)(default = {default_Dmax}):") or default_Dmax)

    #run production simulation
    gravel_production(mu, sig, Dmin, Dmax, N, n)


#run the program
if __name__ == "__main__":
    main()




# end region