########################################################
## Functions to compute the spectral norm / Frobenius ##
## norm of Vandermonde matrices for randomly sampled  ##
## polynomials.                                       ##
########################################################

load("random_polynomials.sage")

## Compute the Vandermonde matrix of a polynomial f
def Vf(f, precision = None):
    if precision is None:
        CCprec = CC
    else:
        CCprec = ComplexField(precision)
    alpha = f.roots(ring = CCprec, multiplicities = False)
    return matrix.vandermonde(alpha,CCprec)

## Home function to compute the Frobenius norm because the default one uses CDF and returns +Infinity when the answer is too big
def Frob_norm(V):
    res = V.base_ring()(0)
    for i in range(V.nrows()):
        for j in range(V.ncols()):
            res += (abs(V[i,j]))^2
    return sqrt(res)
    

     
## Function to compute the spectral/frobenius norm of the vandermonde matrix for randomly sampled polynomials
def tests_Vandermonde(nb_tests, degree,size_poly, norm = "both", sampling_method="unif", precision = None, save_in_file = False, short_output = False, verbose = False, fixed_randomness = False):
    ## Input: nb_tests an integer
    ##        degree and size_poly are integers, defining the (max) degree and the size of the polynomials we want to sample
    ##        norm is either 'frob' for computation of the Frobenius norm, or 2 for computation of the spectral norm, or "both" for computation of both of them
    ##        sampling_method = "unif" or "resultant" (maybe some other sampling method?)
    ##        precision of the computations can be increased by putting precision = x (will compute with x bits of precision)
    ##        if save_if_file = True, will save the results in a file in the "data" repository, with file name of the form Vandermonde_degree_size_nbtests_sampling.sage (where sampling, nbtests, degree and size are replaces by their actual value)
    ##        short_output is a boolean (if True, the polynomial f is not returned)
    ## Output: a list of lenght nb_tests, each element of the list is list of the form [f,norms] (or just [norms] if short_output = True) where norms is a list of spectral/frobenius norms (for the Frobenius norm, we always compute the Frobenius norm with the sage function *and* with our home functions (so we return 2 values, which should be equal, or should be (+infinity, something) when sage function fails)
    if fixed_randomness:
        set_random_seed(42)
    if save_in_file:
        filename = "data/Vandermonde_"+str(degree)+"_"+str(size_poly)+"_"+str(nb_tests)+"_"+sampling_method+".sage"
        f_data = open(filename,"w")
        f_data.write("(nb_tests, degree, size_poly, norm, sampling_method) = "+str((nb_tests, degree, size_poly, norm, sampling_method)))
        f_data.write("\n\nres = [")
        f_data.close()
    
    res = []    
    
    spectral_norm_failed = False
    for test in range(nb_tests):
        f = multi_sampling(degree,size_poly, sampling_method)
        t1 = time.time()
        V = Vf(f, precision = precision)
        if verbose:
                print("time to compute V_f: ", Rshort(time.time()-t1), "s")
        
        norms = []
        if norm == "both" or norm == 2:
            t1 = time.time()
            try:
                norms += [Rshort(V.norm(2))]
            except ValueError:
                norms += [None]
                spectral_norm_failed = True
            if verbose:
                print("time to compute spectral norm: ", Rshort(time.time()-t1), "s")
        if norm == "both" or norm == "frob":
            t1 = time.time()
            norms += [Rshort(V.norm('frob'))]
            if verbose:
                print("time to compute Frobenius norm: ", Rshort(time.time()-t1), "s")
            t1 = time.time()
            norms += [Rshort(Frob_norm(V))]
            if verbose:
                print("time to compute Frobenius norm (home function): ", Rshort(time.time()-t1), "s")
        if short_output:
            res_tmp = norms
        else:
            res_tmp = [f,norms]
        res += [res_tmp]
        if save_in_file:
            f_data = open(filename,"a")
            f_data.write(str(res_tmp)+",\n")
            f_data.close()
  
    if save_in_file:
        f_data = open(filename,"a")
        f_data.write("]")
        f_data.close()
        
    if not spectral_norm_failed:
        if short_output:
            print("max norm observed: ", max([x[0] for x in res]))
            print("average norm observed: ", numpy.mean([x[0] for x in res]))
        else:
            print("max norm observed: ", max([x[1][0] for x in res]))
            print("average norm observed: ", numpy.mean([x[1][0] for x in res]))
    return res

################
## How to use ##
################

#print("Computing the spectral and Frobenius norms for 50 random polynomials")
#res = tests_Vandermonde(10,3,5)
