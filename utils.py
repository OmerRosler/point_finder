from typing import Iterable, Tuple
import itertools
import numpy as np
from numpy.polynomial import Polynomial

def iterate_possible_coefficients() -> Iterable[int]:
    """The pssible polynomial coefficients used by the algorithm

    Returns:
        Iterable[int]: The possible values

    """
    yield -1
    yield 0
    yield 1
    

def iterate_possible_polynomials(m : int, 
    possible_coefficients : Iterable[int]) -> Iterable[Polynomial]:
    """Generates all polynomials of degree m with given coefficients

    Args:
        m (int): The degree of the polynomials
        possible_coefficients (Iterable[int]): The possible coefficients

    Returns:
        Iterable[Polynomial]: _description_
    """
    
    return map( lambda tup : Polynomial(list(tup)), 
    itertools.product(possible_coefficients(), repeat = m))


def pretty_print_polynomial(poly : Polynomial) -> str:
    """Print a polynomial with the '^' token for exponentiation

    Args:
        poly (Polynomial): A polynomial

    Returns:
        str: The representation
    """
    np_output=str(poly)
    return np_output.replace("**", "^")

