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
    furthest_pt = (Polynomial(0), (1.0,1.0))
    smallest_pdt = 0.5
    for poly, roots in find_new_zeros(a, k, m):
        print("Found new roots for follwoing parameters:")
        print ("a = {}, k = {}, m = {}".format(a,k,m))
        print("The polynomial is {}".format(pretty_print_polynomial(poly)))
        r1,r2 = roots
        print("The roots are = {}, {}".format(r1,r2))
        print("and their product is {}".format(r1*r2))
        if r1*r2 < smallest_pdt:
            furthest_pt = (poly, roots)
            smallest_pdt = r1*r2
    if smallest_pdt < 0.5:
        print("Furthest points found - ")
        print("The polynomial is {}".format(pretty_print_polynomial(furthest_pt[0])))
        r1,r2 = furthest_pt[1]
        print("The roots are = {}, {}".format(r1,r2))
        print("and their product is {}".format(r1*r2))

def print_edges_of_curve(k : int, m : int) -> None:
    for poly, pt1, pt2 in find_edges_of_entire_curve(k, m):
        print("Found new roots for follwoing parameters:")
        print ("k = {}, m = {}".format(k,m))
        print("The polynomial is {}".format(pretty_print_polynomial(poly)))
        print("The edges of the curves is {} {},".format(pt1,pt2))

        pts= list(generate_curve(poly,k, m, *pt1, 50))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.title("Curve for polynomial {}".format(pretty_print_polynomial(poly)))

        xs = [float(x[0]) for x in pts]
        ys = [float(x[1]) for x in pts]
        n=len(xs)
        # Plotting with different colors to represent the flow of the curve
        T=np.linspace(0,1,np.size(xs))**2

        # Segment plot and color depending on T
        s = 10 # Segment length
        for i in range(0,n-s,s):
            ax.plot(xs[i:i+s+1],ys[i:i+s+1],color=(0.0,0.5,T[i]), linewidth=1)

        # Plot edges of the curve
        plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]], marker='*', color='r', ls='none')
        plt.show()

from sympy import *
#from sympy.abc import x,y
from scipy import *
from numpy import *

def find_int(x : float, p: Polynomial, y: float, k:int)-> float:
    m=p.degree()
    return x**(m+k+2)*p(y)-y**(m+k+2)*p(y)

def find_traps_for_params(p : Polynomial, k: int):
    m=p.degree()
    ys=linspace(0,1,num=1000)
    xs=linspace(0,1, num=(m+k+2))
    for y in ys:
        sols = optimize.fsolve(find_int, xs, args=(p,y,k))
        unqiue_sols=unique(sols.round(decimals=5))
        potential_roots = filter(is_root_in_range, unqiue_sols)
        for x in potential_roots:
            yield x,y

#TODO: Add filter for non-trivial points
#TODO: Add filter for non-diagonal points
#TODO: Iterate k
#TODO: Iterate p
#TODO: Add test for size


from plot_traps import *


def main():
    #p(y)/p(x)=(y/x)^(m+k+2) => x^(m+k+2)*p(y)-y^(m+k+2)*p(x)=0
    #print_new_roots(a=1.0/np.pi, k = 2, m=8)
    #print_edges_of_curve(k=1, m = 8)
    #plot_g_around_zero(k=1, m = 8)
    # p=Polynomial([1,1,1])
    # k=1
    # for x,y in find_traps_for_params(p,k):
    #     print("solutions is ({},{}):".format(x,y))
    #     print("\n\n")        
    plot_curves(1,5)


if __name__ == "__main__":
    main()