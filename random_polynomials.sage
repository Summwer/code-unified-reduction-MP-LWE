#################################################################
## This file contains tool functions used for sampling random  ##
## number fields / polynomials (with different distributions)  ##
## and will be used in the other files                         ##
#################################################################

import sys, signal, time, numpy

Zx.<X> = ZZ[]
Cshort = ComplexField(15)
Rshort = Reals(15)


## Generate a random monic irreducible polynomial of degree n, with coefficients taken uniformly at random between -B and B (default B is 10)     
def unif_poly(n,B = 10):
     f = Zx(X^n)+Zx.random_element(n-1,x = -B,y = B) # Hopefully this should work with both Sage 9 and Sage 10...
     while not f.is_irreducible():
          f = Zx(X^n)+Zx.random_element(n-1,x = -B,y = B)
     return f
     
## generate a random monic irreducible polynomial of degree n (provided n is composite), using resultant of two bivariate polynomials
Zxy.<Y> = PolynomialRing(Zx)
def one_resultant_poly(n1,n2,B = 10):
    f = sum([Zx([ZZ.random_element(-B,B+1) for _ in range(n1+1)])*Y^i for i in range(n2+1)])
    g = sum([Zx([ZZ.random_element(-B,B+1) for _ in range(n1+1)])*Y^i for i in range(n2+1)])
    res = f.resultant(g)
    return res
    
## Trying to generate a different family of random fields, by taking the resultant of two random polynomials (idea from Aurel Page, apparently these number fields come from curves and may have different properties). The issue is that we can reject until the resultant polynomial is irreducible, but it is very unlikely that we end up with a monic polynomial.
## We could use the polred function in Pari/GP that transforms any irreducible polynomial into another monic irreducible polynomial generating the same field
## TODO for the moment we are letting this sampling method aside, but we could try to implement it with polred (cf above)

def resultant_poly(n,B = 5):
    #n2 = (n-1).factor()[0][0]
    n2 = ZZ.random_element(2,ceil(sqrt(n)))
    n1 = (n-1)//(2*n2)
    print("\nn1, n2 = ", n1, n2)
    f = one_resultant_poly(n1,n2,B = B)
    while not f.is_irreducible():
        f = one_resultant_poly(n1,n2,B = B)
    return f
    
## Generate a random Bugeaud-Mignotte polynomial of degree n, with coefficients at most B, and k roots close to one another (default k is 2)
def random_Bugeaud_Mignotte(n,B = None, k = 2):
     if B is None:
          B = 10
     if k > n/2:
          print("k should be at most n/2, I am taking k = floor(n/2)")
          k = floor(n/2)
     A = floor((B/2)^(1/2)) # bound on a
     if A < 2:
          print("B is too small, I am taking B larger")
          A = 2
     a = ZZ.random_element(2,A+1)
     f = Zx(X^n-2*(a*X-1)^k)
     assert(f.is_irreducible())
     return f
    
 
# Combine multiple sampling methods in one function
def multi_sampling(degree,size_poly, sampling_method):
    if sampling_method == "unif":
        f = unif_poly(degree, B = size_poly)
    elif sampling_method == "resultant":
        f = resultant_poly(degree, B = size_poly)
    elif sampling_method == "mignotte":
        f = random_Bugeaud_Mignotte(degree, B = size_poly)
    else:
        print("Incorrect sampling method in multi_sampling")
        return
    return f


## Generate a uniformly random element in the field K = Q[x]/f(x)
## Input: f - an irreducible polynomial defining the field
##        B - bound for coefficients (default 10)
## Output: a random polynomial a(x) of degree < deg(f) with integer coefficients in [-B, B]
def random_field_element(f, B=10):
    # Get the polynomial ring
    R = f.parent()
    n = f.degree()
    # Generate random coefficients
    coeffs = [ZZ.random_element(-B, B+1) for _ in range(n)]
    # Construct polynomial a(x) = sum_{i=0}^{n-1} coeffs[i] * x^i
    a_poly = sum(coeffs[i] * R.gen()^i for i in range(n))
    return a_poly
