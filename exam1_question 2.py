#region import
import math
from HW3_SP25_Soln import HW3b, numericalMethods
from exam1_question1 import gravel_production, calculate, log_normal, lognormal_cumulative, truncated_log

#end region
def degrees_of_freedom(var_A, var_B, m_A, m_B):

    if m_A  <= 1 or m_B <= 1:
        raise ValueError("Sample size must be greater than 1")

    numerator = (var_A / m_A + var_B / m_B) **2
    denominator = ((var_A / m_A)**2) / (m_A - 1) + ((var_B / m_B) **2) / (m_B - 1)
    return numerator / denominator
def t_statistic(mean_A, mean_B, var_A, var_B, m_A, m_B):
    numerator = mean_A - mean_B
    denominator = math.sqrt((var_A / m_A) + (var_B / m_B))
    return numerator / denominator


def main():

      getOut = False
      while not getOut:
          # Establish user inputs
          m_A = int(input("Enter the sample size for supplier A").strip().lower())
          m_B = int(input("Enter the sample size for supplier B").strip().lower())
          mean_A = float(input("Enter the mean for supplier A").strip().lower())
          mean_B = float(input("Enter the mean for supplier B").strip().lower())
          var_A = float(input("Enter the sample variance for supplier A").strip().lower())
          var_B = float(input("Enter the sample variance for supplier B").strip().lower())

          # Perform the test
          df = degrees_of_freedom(var_A, var_B, m_A, m_B)
          t_stat = t_statistic(mean_A, mean_B, var_A, var_B, m_A, m_B)
          print(f"Calculated t-statistic: {t_stat:.3f}")

          # Calculate the critical t-value for alpha = 0.05 (one-sided test)
          critical_t = numericalMethods.Simpson(truncated_log, 4,100)

          print(f"Critical t-value {critical_t:.3f}")

          if t_stat > critical_t:
            print("Supplier B's gravel size is significantly smaller than supplier A")
          else:
            print("Supplier B's claim of smaller diameters is not statistically significant")

          getOut = input("Do you want to test again? (Y/N):").strip().lower() == "n"


if __name__=='__main__':
    main()

#end region
