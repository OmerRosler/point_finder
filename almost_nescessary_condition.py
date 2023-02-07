import numpy as np
import matplotlib.pyplot as plt

from utils import *


def almost_nescessary_condition(p : Polynomial):
    m=p.degree()+1


    #TODO: Replace condition 2 with its' square for better performance
    def fn(x,y):
        
        #This is the translation vector
        tpx = 2*x**(-m)*p(x)
        tpy= 2*y**(-m)*p(y)

        #The slope and length of `t_p`
        s_tp = tpy/tpx
        length_tp=np.sqrt(tpx**2+tpy**2)

        #This is the optimal `k` for a trap
        k = -m-np.floor(np.log(p(y)/p(x))/np.log(y/x))

        #Condition 1
        cond1 = k >=0

        # The slope of the line on the convex hull
        s_Lk=(y/x)**(k+1)

        #cond_not_trap_exactly = (s_Lk != s_tp)

        #if s_tp < s_Lk:
        #    raise AssertionError("Does not make sense")

        #The formula for angle in terms of slopes
        tan_theta = (s_tp-s_Lk)/(1+s_tp*s_Lk)

        #Condition 3 not met
        cond3 = (tan_theta > 0)

        #Calculate `sin` and `cos` in terms of `tan`
        denom = np.sqrt(tan_theta**2 +1)
        sin_theta = tan_theta/denom
        cos_theta = 1/denom

        # Condition 2
        cond2 = (length_tp*(sin_theta + 2*y**(k+1)*cos_theta)<2*y**(k+1))
        return (cond1 & cond3 & cond2)
    return fn

for p in iterate_possible_polynomials(7, iterate_possible_coefficients()):
    above_diagonal = lambda x,y: x<y
    non_trivial = lambda x,y : x*y <0.5
    range_of_tp = almost_nescessary_condition(p)

    d = np.linspace(0.5,1,600)
    x,y = np.meshgrid(d,d)
    results = (above_diagonal(x,y) & non_trivial(x,y) & range_of_tp(x,y)).astype(int)
    #results = range_of_tp(x,y).astype(int)

    if results.any():
        print(p)
        plt.imshow( results, 
                        extent=(x.min(),x.max(),y.min(),y.max()),origin="lower", cmap="Greys")

plt.show()