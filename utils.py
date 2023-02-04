from typing import Iterable, Tuple, List
import itertools
import infinite
import numpy as np
from numpy.polynomial import *

def iterate_possible_coefficients() -> Iterable[int]:
    """The pssible polynomial coefficients used by the algorithm

    Returns:
        Iterable[int]: The possible values

    """
    return [-1,0,1]
    #yield -1
    #yield 0
    #yield 1
    
def iterate_possible_polynomials(m : int, 
    possible_coefficients : List[int]) -> Iterable[Polynomial]:
    """Generates all polynomials of degree m with given coefficients

    Args:
        m (int): The degree of the polynomials
        possible_coefficients (Iterable[int]): The possible coefficients

    Returns:
        Iterable[Polynomial]: _description_
    """
    args_to_product = [possible_coefficients]*m

    return map( lambda tup : Polynomial(list(tup)), 
    itertools.product(*args_to_product))

def pretty_str_polynomial(poly : Polynomial) -> str:
    """Print a polynomial with the '^' token for exponentiation

    Args:
        poly (Polynomial): A polynomial

    Returns:
        str: The representation
    """
    np_output=str(poly)
    return np_output.replace("**", "^")

