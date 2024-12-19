import numpy as np

"""Calculate drift velocity and diffusion coefficient for hopping in a
1D chain of sites, given the hopping rates.  The hopping rates are
given in two arrays: forward and reverse, defined like this:
forward: [0->1, 1->2, 2->3, ... ]
reverse: [1->0, 2->1, 3->2, ... ]
where 0->1 means the rate of hopping from site 0 to site 1.

The function returns a tuple (V, D) where V is drift velocity and D is
the diffusion constant.

Note that the input is different from previos C implementations which
used the rates of going left or right from site n, resulting in one of
the input arrays being shifted by one element.

The program make use of the analytic methods by B. Derrida,
implemented in a way that permits evaluation in time and memory linear
in the number of sites (see references below).

The calculation is assymmetric, and works better if the average
velocity is positive, i.e. towards larger site indices. If not,
round-off errors tend to accumulate. If this is detected, the problem
is mirrored, the mirrored problem is solved and the result is returned
with the opposite sign for V.

Fredrik Jansson, 2023


References

The analytic solution:
B. Derrida, J. Stat. Phys. 31, 433 (1983).

Thesis, F. Jansson, Charge transport in disordered materials -
simulations, theory, and numerical modeling of hopping transport and
electron-hole recombination.
Åbo Akademi University, 2011
https://urn.fi/URN:NBN:fi-fe201311277464
Implementation details in section 3.6, especially a method to evaluate
the expressions in linear time. Application in section 6.2.

Effect of Electric Field on Diffusion in Disordered Materials I.
One-dimensional Hopping Transport, A. V. Nenashev, F. Jansson,
S. D. Baranovskii, R. Österbacka, A. V. Dvurechenskii, F. Gebhard,
Phys. Rev. B 81, 115203 (2010)
http://dx.doi.org/10.1103/PhysRevB.81.115203

"""
def VD(forward, reverse, allow_mirror=True):
    forward = np.array(forward) # convert inputs to numpy arrays 
    reverse = np.array(reverse)
    N=len(forward)
    
    g = np.roll(reverse,1) / forward # g[n] = revrese[n-1] / forward[n]
    h = reverse / forward            # h[n] = revrese[n] / forward[n]
    r = np.zeros(N)
    u = np.zeros(N)
    
    # Eq. 50
    Gn = 1
    t = 1    
    for i in range(1, N):      
        t *= g[i-1] 
        Gn += t
    
    G = t*g[N-1]  # product of all g

    if (G > 1e12 or not np.isfinite(G)):
#       print (f"G={G} Result may be inaccurate due to drift to the left.")
        if allow_mirror:
#           print("Mirroring the chain.")
            V,D = VD(reverse[::-1], forward[::-1], allow_mirror=False)
            return -V, D
        else:
            print("Giving up.")
            return None, None

    rsum = 0
    for n in range(N-1,-1,-1):  
        r[n] = Gn / forward[n]
        Gn = g[n] * Gn - G + 1
        rsum += r[n]
  
    # Eq. 46
    Hn = 1
    t = 1    
    for i in range(1, N):
        t *= h[N-i]
        Hn += t
    
    H = t*h[0]  #product of all h (= G, product of all g)
    # print (f"G: {G:10.15g}\nH: {H:10.15g}\n")
     
    usum = 0
    for n in range(0, N):
        u[n] = Hn / forward[n]
        Hn = h[n] * Hn - H + 1
        usum += u[n]

    # Eq. 49
    V = N/rsum * (1-H)

    # rsum = usum according to Eq 48
    # print (f"rsum: {rsum:10.15g}\nusum: {usum:10.15g}\n")
  
    # Eq. 47
    n = N-1     # start from the largest n
    s = 0
    for i in range(1,N+1):       # initial "triangle"-sum
        s += i * r[i-1]          # was r[(n+i)%N] but n=N-1        

    t = 0
    x = 0
    for n in range(N-1,-1,-1):
        t += u[n] * s
        s += rsum - N*r[n]
        x += forward[n] * u[n] * r[n]
    
    t *= V
    x *= N  
    D = 1 / (rsum*rsum) * (t + x) - V * (N+2.0)/2.0 ;
    return V, D


if __name__ == '__main__':
    # if this script is run directly i.e. not imported,
    # run a small test calculation
    
    forward = [.1, .2, .4, .2, .03]
    reverse = [.1, .8, .2, .1, .2]

    print(VD(forward, reverse))
    # test mirroring the chain, should give negative V and same D
    print(VD(reverse[::-1], forward[::-1]))
    
