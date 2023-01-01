from utils import *

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
    new_poly = poly + lhs
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
    for int_poly in iterate_possible_polynomials(m, iterate_possible_coefficients):

        poly = perturb_polynomial(a, k, int_poly, m)
        potential_roots = filter(is_root_in_range, poly.roots())

        for r1,r2 in non_trivial_double_roots(potential_roots):
            yield (poly, (r1,r2))
