import matplotlib.pyplot as plt
import numpy as np
from utils import *

def generate_r(p : Polynomial, k : int, m = None):
    """The condition `s_{t_p}=s_{P_k}` is equivalent to
    r(x,y)=0 for some algebraic function `r` which this function generates
    from p and k

    Args:
        p (Polynomial): The polynomial defining t_p
        k (int): P_k is the line we test against
        m (int, optional): The degree of m +1. Defaults to None.

    Returns:
        _type_: A function whose zeros are potential dense-traps
    """
    if m is None:
        m=p.degree()+1
    def curve(x,y):
        return x**(m+k+1)*p(y)-y**(m+k+1)*p(x)
    return curve

def is_tp_below(x : float, y: float, p : Polynomial, k : int, m = None) -> bool:
    """Checks if t_p is below pi(p^k m^inf)
    This condition is required for traps

    The formula was derived in [H-S]

    Args:
        x (float): The x coordinate of the point to test
        y (float): The y coordinate of the point to test
        p (Polynomial): The poylnomial defining t_p
        k (int): pi(p^k m^inf) is the point we check againts
        m (int, optional): The degree of p + 1. Defaults to None.

    Returns:
        bool: Is the condition met?
    """
    if m is None:
        m=p.degree()+1
    if x == 1 or y==1:
        return False
    is_coord_below = lambda a: 2*p(a)<a**(m+k+1)
    return is_coord_below(x) and is_coord_below(y)

import typing

def is_tp_critical(x : float, y: float, p : Polynomial, k : int, m = None) -> bool:
    """The condition was derived by algebraic manipulations of combining the conditions
    r(x,y)=dr/dx(x,y)=dr/dy(x,y)=0

    where 
    r(x,y)=x^{m+k+1}p(y)-y^{m+k+1}p(x)

    Args:
        x (float): The x coordinate of the point which satisfy r(x,y)=0
        y (float): The y coordinate of the point which satisfy r(x,y)=0
        p (Polynomial): The polynomial used
        k (int): The index of the line to compare
        m (_type_, optional): The degree of p + 1. Defaults to None.

    Returns:
        bool: Is the condition met?
    """
    if m is None:
        m=p.degree()+1
    is_critical_coord = lambda a: (m+k+1)*p(a)== a*p.deriv()(a)
    return is_critical_coord(x) and is_critical_coord(y)

def filter_traps_only(plot, p : Polynomial, k : int) -> Iterable[typing.Tuple[float,float]]:
    """For each point on the plot we test if it is really a potential dense-trap

    There are a number of conditions:
    1. The point is not on the diagonal (as for these the [H-S] characterization of the convex hull fails)
    2. The point is not in the trivial part of N (because we don't care)
    3. t_p is below pi(p^{k+1}m^inf) - otherwise t_p won't be trap-like
    4. grad(r(x))!=0 this condition can be weakend - if this point is not an extrema 
    then there is a direction for which r>0 which implies t_p is on the left of P_k
    which is required for trap-like

    Args:
        plot (_type_): The plot to choose points from
        p (Polynomial): The polynomial to test
        k (int): P_k is the line we choose

    Raises:
        AssertionError: The case of an extrema is a rare event - we raise this assertion as we want to be notified if such a point even exists

    Returns:
        Iterable[Tuple[float,float]]: The points that satisfy all the conditions
    """
    for path in plot.collections[0].get_paths():
        v= path.vertices
        for x,y in v:
            if np.isclose(x,y):
                continue
            #if x*y >= 0.5:
            #    continue
            if not is_tp_below(x,y, p, k+1):
                continue
            #These points are potential traps but maybe there is no direction to make t_p above P_k
            if is_tp_critical(x,y, p, k):
                raise AssertionError("Very sepcial problematic point: {}".format(pretty_print(p)))
            yield x,y

def fill_trivial_part(ax, num_pts = 100):
    """Used by the plotting function to mark the trivial part of N
    in the same color as used before

    Args:
        ax (_type_): The plot used
        num_pts (int, optional): Num of points to run the algorithm. Defaults to 100.
    """
    xs=np.linspace(0+1.0/num_pts,1,num_pts)
    ys=0.5/xs
    ax.fill_between(xs, 1, ys, where=ys<1, color='#969696')


def plot_curves(max_k :int, max_m: int, num_pts = 600):
    """
    Plots all the curves that satisfy the condition
    s_(t_p)=s_(P_k)

    with different colors

    also print all points on these curves that are potential dense-traps
    (i.e. there is a perturbation direction of the parameter for which we get traps)

    Note: The text output of this function is to be read by a
    C++ program that plots all these points on
    a given image (to represent the curves on the image of N)

    Args:
        max_k (int): The maximal index of the vertices in the description
        the line we test is P_k=[pi(p^k m^inf), pi(p^{k+1} m^inf)]
        max_m (int): The maximal length of the words with the trap
        the vector t_p is (2x^{-m}p(x),2y^{-m}p(y)) where `p` is a degree `m-1` poynomial
        with coefficients {-1, 0, 1}
        num_pts (int, optional): Num of points to run the algorithm. Defaults to 600.
    """
    y, x = np.ogrid[0.5:1:num_pts*1j, 0.5:1:num_pts*1j]
    fig, ax = plt.subplots(figsize=(1,1), dpi=600)
    for p in iterate_possible_polynomials(max_m, iterate_possible_coefficients()):
        for k in range(max_k):
            curve_fn = generate_r(p, k)
            #Generate the points where `curve_fn=0`
            res = ax.contour(x.ravel(), y.ravel(), curve_fn(x,y), [0], colors= [np.random.rand(3,)], linewidths=0.3)
            # Filter actual trap points based on other conditions
            for i, (a,b) in enumerate(filter_traps_only(res, p, max_k)):
                # Print header of polynomial + `k`
                if i == 0:
                    print("STOP")
                    print(pretty_str_polynomial(p))
                    print(k)
                # Print all points on the crurve
                print("{} {}".format(a,b)) 
    fill_trivial_part(ax, num_pts)
    #plt.axis('off')
    #plt.show()
