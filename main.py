from trap_exactly import *

def print_new_roots(a : float, k : int, m : int) -> None:
    """ Print all points with a potential trap when perturbing the point
    a*(p^{k+1}m^{\inf} - p^{k}m^{\inf}) to the left

    Args:
        a (float): A number in the range 0 < a < 1
        k (int): Represent the line on the convex hull to search the trap on
        m (int): The maximal degree of search. Note the number of polynomials is 3^m
    """
    assert( a>0 and a <1)
    for poly, roots in find_new_zeros(a, k, m):
        print("Found new roots for follwoing parameters:")
        print ("a = {}, k = {}, m = {}".format(a,k,m))
        print("The polynomial is {}".format(pretty_print_polynomial(poly)))
        r1,r2 = roots
        print("The roots are = {}, {}".format(r1,r2))
        print("and their product is {}".format(r1*r2))

def main():
    print_new_roots(a=1.0/np.pi, k = 1, m=8)

if __name__ == "__main__":
    main()