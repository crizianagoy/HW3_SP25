#region imports
from math import sqrt, pi, exp
from numericalMethods import GPDF, Simpson, Probability, Secant
#endregion

def find_c_for_given_P (P_target, args, GT):
    """
    Use Secant to find c to results in probability P_target
    :param P_target: probability input by user
    :param args: tuple (mean, st Dev)
    :param GT: Boolean determine probability x>c or x<c
    :return: calculate c
    """
    def root_function(c):
        return Probability(GPDF, args, c, GT)-P_target

    #initial guesses for Secant Method:
    c0, c1= args[0]-args[1], args[0] +args[1]
    c_solution, _= Secant(root_function, c0,c1,maxiter=20, xtol=1e-6)
    return c_solution

#region function definitions
def main():
    """
    I want to integrate the Gaussian probability density function between
    a left hand limit = (mean - 5*stDev) to a right hand limit = (c).  Here
    is my step-by-step plan:
    1. Decide mean, stDev, and c and if I want P(x>c) or P(x<c).
    2. Define args tuple and c to be passed to Probability
    3. Pass args, and a callback function (GPDF) to Probability
    4. In probability, pass along GPDF to Simpson along with the appropriate args tuple
    5. Return the required probability from Probability and print to screen.
    :return: Nothing to return, just print results to screen.
    """
    #region testing user input
    print("\nGaussian Probability Calculation\n")
    # The following code solicites user input through the CLI.
    mean = float(input("Population mean?: "))
    stDev = float(input("Standard deviation?:"))
    choice = input("are you providing c (for P) or P (for c)?(enter c or P):").strip().lower()
    args = (mean, stDev)  # define args

    if choice == "c":
        c = float(input("c value?: "))
        GT = input("Probability greater than c? (y/n): ").strip().lower() in ["y", "yes", "true"]
        P = Probability(GPDF, args, c, GT)
        print(f"P(x {'>' if GT else '<'} {c:.2f} | N({mean:.2f}, {stDev:.2f})) = {P:.4f}")

    elif choice == 'p':
        P_target = float(input("Enter probability value (P): "))
        GT = input("Probability greater than c? (y/n): ").strip().lower() in ["y", "yes", "true"]
        c = find_c_for_given_P(P_target, args, GT)
        print(f"c value for P(x {'>' if GT else '<'} c | N({mean:.2f}, {stDev:.2f})) = {c:.4f}")

    else:
        print("Invalid choice. Please enter 'c' to find P or 'P' to find c.")

#endregion

if __name__ == "__main__":
    main()