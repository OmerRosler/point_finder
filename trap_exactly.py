from utils import *
import matplotlib.pyplot as plt


def is_root_in_range(r : np.number) -> bool:
    """A predicate whether a number satisfies 0<r<1

    Args:
        r (np.number): A number

    Returns:
        bool: Result
    """
    if not np.isreal(r):
        return False
    return r >0 and r<1

def empty():
    """Empty Generator

    """
    yield from ()

def iterate_roots(roots : Iterable[np.number]) -> Iterable[Tuple[np.number, np.number]]:
    """Iterate on all different pairs from a list of numbers

    Args:
        roots (Iterable[np.number]): A list of numbers with possible repetition

    """
    unique_roots = set(roots)
    if len(unique_roots) < 2:
        return empty()
    return itertools.permutations(unique_roots, 2)

def non_trivial_double_roots(roots : Iterable[np.number]) -> Iterable[Tuple[np.number, np.number]]:
    """Given a list of roots of a polynomial return all pairs 
    (x,y) that satisfy 0<x<1, 0<y<1, x!=y and x*y<0.5

    These zeros are non-trivial in N and might not be in N' as well

    Args:
        roots (Iterable[np.number]): All roots of a common polynomial

    Returns:
        Iterable[Tuple[np.number, np.number]]: All pairs that satisfy the condition

    """
    #TODO: Also verify a possible perturbation direction to a positive value of h (otherwise we won't get the trap-like property)
    #return filter(lambda roots: roots[0]*roots[1]<0.5, iterate_roots(roots))
    unique_roots = set(roots)
    if len(unique_roots) < 2:
        return False
    for r1, r2 in itertools.permutations(unique_roots, 2):
        if r1 * r2 < 0.5:
            yield r1, r2

def perturb_polynomial(a : float, k : int, poly : Polynomial, m : int = None) -> Polynomial:
    """Take a polynomial of degree m and add to it -a*x^(k+m+2)

    Args:
        a (float): A parameter in the range 0<a<1
        k (int): A positive degree
        poly (Polynomial): The original polynomial
        m (int, optional): The degree of `poly`

    Returns:
        Polynomial: A modified polynomial
    """
    if m is None:
        m = poly.degree()
    lhs = [0] * (k + 2 + m)
    lhs[-1]=-a
    new_poly = poly + Polynomial(lhs)
    return new_poly

def find_new_zeros(a : float, k : int, m : int) -> Iterable[Tuple[Polynomial, Tuple[np.number, np.number]]]:
    """Searches all possible polynomials with coefficients in {-1,0, 1}
    that have pairs of non trivial zeros to the function 
    h_k^a(x)=ax^{k+1}-x^{-m}p(x)

    where 0 < a < 1 and k is an integer

    All such points are candidates to be in the closure of the interior of N
    
    Any small enough perturbation that makes h_k^a positive for both roots -
    will be a trap!

    Args:
        a (float): A number in the range 0 < a < 1
        k (int): Represent the line on the convex hull to search the trap vector on
        m (int): The maximal degree of search. Note the number of polynomials is 3^m

    Returns:
        Iterable[Tuple[Polynomial, Tuple[np.number, np.number]]]: 
            The polynomials and their pairs of roots that satisfy the conditions

    """
    for int_poly in iterate_possible_polynomials(m, iterate_possible_coefficients()):

        poly = perturb_polynomial(a, k, int_poly, m)
        potential_roots = filter(is_root_in_range, poly.roots())

        for r1,r2 in non_trivial_double_roots(potential_roots):
            yield (poly, (r1,r2))


def find_edges_of_entire_curve(k : int, m: int, a_start = 0, a_end = 1) -> Iterable[Tuple[Polynomial, Tuple[float, float], Tuple[float, float]]]:
    """Find a polynomial with the relevant coefficients with a pair of zeros for both
    p(x)-`a_start`*x^{m+k+2}
    and
    p(x)-`a_end`*x^{m+k+2}

    For such polynomials there is an entire curve inside $\overline{N^o}$
    connecting the pairs of zeros. This function calculates the edges of these curves.

    Args:
        k (int): The relative degree we check for
        m (int): The degree of the polynomial
        a_start (int, optional): Starting value to check. Defaults to 0.
        a_end (int, optional): End value to check. Defaults to 1.

    Returns:
        Iterable[Tuple[Polynomial, Tuple[float, float], Tuple[float, float]]]: 
        The polynomials and the respective roots for which the roots are in the relevant range and the start root is not trivial

    """
    for poly in iterate_possible_polynomials(m, iterate_possible_coefficients()):
        # Find non-trivial zeros for the perturbation of `a_end`
        left_poly = perturb_polynomial(a_end, k, poly, m)
        potential_roots = filter(is_root_in_range, left_poly.roots())
        zeros = list(non_trivial_double_roots(potential_roots))
        if not zeros:
            continue
        # Find any pair of zeros for the perturbation of `a_start`
        # Note: we allow it to be inside N^t i.e. no filtering
        right_poly = perturb_polynomial(a_start, k, poly, m)
        roots_for_new_poly=filter(is_root_in_range, right_poly.roots())
        
        # Locate closest zero of the left perturbation and yield the pair of zeros
        for r1,r2 in iterate_roots(roots_for_new_poly):
            r1_orig = min(zeros, key=lambda t: abs(t[0]-r1))[0]
            r2_orig = min(zeros, key=lambda t: abs(t[1]-r2))[1]
            yield poly, (float(r1),float(r2)), (float(r1_orig), float(r2_orig))

def generate_curve(poly : Polynomial, k :int , m: int, r1 : float, r2 : float, num_of_pts : int) -> Tuple[float,float]:
    """Generate a list of points on the curve starting with (r1,r2)
    and ending with a pair of zeros of `poly-x^{m+k+2}`
    entirely inside $\overline{N^o}$

    We run on `a` from 0 to 1 and find a pair of zeros for
    `poly-a*x^{m+k+2}` close to the initial zero (r1,r2)

    Args:
        poly (Polynomial): The original polynomial
        k (int): The relative degree we check for
        m (int): The degree of the polynomial
        r1 (float): The first relevent zeros of `poly` in the range
        r2 (float): The second relevent zeros of `poly` in the range
        num_of_pts (int): The number of points to generate (evenly spaced)

    Raises:
        NotImplemented: So far only polynomials with two zeros are acoounted for

    Returns:
        Tuple[float,float]: Points on the curve

    """
    if m is None:
        m=poly.degree()
    alpha_space=np.linspace(0+1.0/num_of_pts,1,num_of_pts)
    for a in alpha_space:
        new_poly = perturb_polynomial(a, k, poly, m)
        zeros=list(filter(is_root_in_range, new_poly.roots()))
        for new_r1,new_r2 in iterate_roots(zeros):
            #TODO: Add proper filtering in the case of more than two zeros
            if len(zeros)!= 2:
                raise NotImplemented("Not supporting polynomials with more than two zeros at the moment")
            yield new_r1,new_r2

