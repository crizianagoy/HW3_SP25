#region imports
from math import gamma, sqrt,pi
from numericalMethods import Simpson
#endregion

def t_distribution_pdf(u,m:int):
    """
    compute Probability Density Function of t-distribution
    :param u: t-value
    :param m: dregrees of freedom
    :return: PDF value at u
    """
    Km=gamma((m+1)/2)/(sqrt(m*pi)*gamma(m/2))
    return Km*(1+(u**2/m))**(-(m+1)/2)

def t_distribution_probability(z,m):
    """
    Cummulatie probability F(z)= P(T<=z)
    :param z: Upper Limit
    :param m: degrees of freedom
    :return: P(T<=z)
    """
    lower_limit= -10
    probability = Simpson(lambda u: t_distribution_pdf(u, m), lower_limit, z)
    return probability

def main():
    """
    user input m and z
    compute cummulative porbability
    """
    print("\nt-distribution probability calculation\n")

    #solve for m value
    m=int(input("enter degrees of freedom(7,11 or 15):"))

    if m not in[7,11,15]:
        print("Pick 1 of the 3: 7,11,or 15 please:")
        return

    #solve for t-value
    z_values=[]
    for i in range(3):
        z=float(input(f"Enter t-value {i+1}:").strip())
        z_values.append(z)

    #solve n
    for z in z_values:
        p_t= t_distribution_probability(z,m)
        print(f"computed: F({z:.5f} | m={m})= {p_t: .5f}")
#endregion

if __name__=="__main__":
    main()
#endregion