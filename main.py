from plot_traps import *

def plot_poynomial_CT_vectors(p: Polynomial, num_pts: int, abs_bound):
    coord_value = lambda x: x**(-p.degree()-1)*p(x)
    abs_bound = lambda x,y : x**4+y**4
    vec = lambda x,y : coord_value(x)**2+coord_value(y)**2

    d = np.linspace(0.5,1,num_pts)
    m = np.meshgrid(d,d)
    x,y=m

    res =vec(x,y)<abs_bound(x,y)

    if np.all(res==False):
        return False
    
    sx,sy = (np.where(res == True))
    values = (np.multiply(x[sx,sy], y[sx,sy]))
    if np.any(values<0.5):
        print("non trivial")

        im = plt.imshow( res.astype(int) , 
                        extent=(x.min(),x.max(),y.min(),y.max()),origin="lower", cmap="Greys")
        
        print(p)
        plt.show()
    
    return True

def main():
    num_pts = 600
    abs_bound = lambda x,y : x**4+y**4
    for p in iterate_possible_polynomials(5,[-2,-1,0,1,2]):
        plot_poynomial_CT_vectors(p,num_pts,abs_bound)
    #y, x = np.ogrid[0.5:1:num_pts*1j, 0.5:1:num_pts*1j]
    #fig, ax = plt.subplots(figsize=(1,1), dpi=600)


    



if __name__ == "__main__":
    main()