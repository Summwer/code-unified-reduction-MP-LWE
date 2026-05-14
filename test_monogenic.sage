###########################################################
## Functions to test whether a number field is monogenic ##
## and compute the norm of the conductor                 ##
###########################################################

load("random_polynomials.sage")


##########################
## Testing monogenicity ##
##########################
    
## Function to compute the proportion of f that have a prime discriminant
## (when this is the case, then the order ZZ[X]/f is maximal in K = QQ[X]/f)
## this is not exactly what we want (this only gives a lower bound on the
## number of monogenic number fields)
## but at least we can compute it relatively efficiently and go to middle-high dimensions
def tests_prime_disc(nb_tests, degree,size_poly, sampling_method="unif", save_in_file = False, short_output = False, fixed_randomness = False):
    ## Input: nb_tests an integer
    ##        degree and size_poly are integers, defining the (max) degree and the size of the polynomials we want to sample
    ##        sampling_method = "unif" or "mignotte" (maybe later also "resultant")
    ##        if save_if_file = True, will save the results in a file in the "data" repository, with file name of the form prime_disc_degree_size_nbtests_sampling.sage (where sampling, nbtests, degree and size are replaced by their actual value)
    ##        if short_output is True, it will not return(/write in the file) the list of polynomials, only the result "prime discriminant/composite discriminant" and the expected probability of being prime will be returned/written
    ##        if fixed_randomness is true, the random seed will be fixed to 42 at the beginning of the function (for reproducibility)
    ## Output: a list of lenght nb_tests, each element of the list is a list of the form [f,b,1/log_df], where f is a random polynomial, b is 1 if the discriminant of f is prime, and 0 otherwise, and log_df is the logarithm (in base e) of the discriminant of f. If df behaves as a uniformly random integer, we expect that the sum of the b's is roughly equal to the sum of the 1/log_df (if we sum over sufficiently many f's). If short_output = True, then the list contains only pairs of the form [b,1/log_df] (without the f).
    if fixed_randomness:
        set_random_seed(42)
    if save_in_file:
        filename = "data/prime_disc_"+str(degree)+"_"+str(size_poly)+"_"+str(nb_tests)+"_"+sampling_method+".sage"
        f_data = open(filename,"w")
        f_data.write("(nb_tests, degree, size_poly, sampling_method) = "+str((nb_tests, degree, size_poly, sampling_method)))
        f_data.write("\n\nres = [")
        f_data.close()
    
    res = []    
    
    for test in range(nb_tests):
        f = multi_sampling(degree,size_poly, sampling_method)
        
        ## test if disc(f) is prime
        df = abs(f.discriminant())
        if df.is_prime(proof=False):
            b = 1
        else:
            b = 0
        if short_output:
            tmp_res = [b,Rshort(1/log(df))]
        else:
            tmp_res = [f,b,Rshort(1/log(df))] ## keep the knowledge of f
        res += [tmp_res]
        if save_in_file:
            f_data = open(filename,"a")
            f_data.write(str(tmp_res)+",\n")
            f_data.close()
  
    if save_in_file:
        f_data = open(filename,"a")
        f_data.write("]")
        f_data.close()
        
    
    print("observed proportion of prime discriminants:", Rshort(sum([x[-2] for x in res])/nb_tests))
    print("expected proportion of prime discriminants:", Rshort(sum([x[-1] for x in res])/nb_tests))
    return res


#############################
## Computing the conductor ##
#############################

## Below are functions to compute the norm of the conductor ideal of ZZ[X]/f(X). This allows to precisely detect which number fields are monogenic 
## (they are monogenic iff the norm of their conductor ideal is 1), 
## and it also allows us to get an estimate of the loss that will occure in our reduction for the fields that are not monogenic.
## The drawback of this function is that it requires to compute the discriminant of the field K, which requires factoring the discriminant of f, 
## and hence quickly becomes out of reach when the degree of f increases
## (the test function allows to specify a maximum time to wait for the factorization, but this biases the results towards fields whose discriminant is prime, or
## whose discriminant is a big prime times a small integer (which are easier to factor)

def norm_conductor_ideal(f,z,K, disc=None):
    ## Input: z an element of OK (ring of integer of some number field K) and f the minimal polynomial of z, which must be of degree [K:Q] (so that ZZ[z] is an order in K).
    ##        discriminant of K can be provided to fasten the computation
    ## Output: the norm of the conductor ideal of the order ZZ[z], which equalt N(f'(z))/disc(K)
    f_prime = f.derivative()
    y = f_prime(z)
    if disc is None:
        disc = K.discriminant() ## costly step when K has a large degree
    return abs(y.norm()/disc)
    
    
def tests_norm_conductor(nb_tests,degree,size_poly, sampling_method="unif", max_time_per_test = 10, save_in_file = False, short_output = False, fixed_randomness = False, verbose = False):
    ## Input: nb_tests an integer
    ##        degree and size_poly are integers, defining the degree and the size of the polynomials we want to sample
    ##        sampling_method = "unif" or "mignotte" (maybe later also "resultant")
    ##        max_time_per_test is the maximum time it will try to compute the discriminant (in second, default is 10). If the discriminant has not been computed in the maximum amount of time, it just stops the computation and goes to the next polynomial
    ##        if fixed_randomness = True, the random seed will be fixed to 42 at the begining of the execution
    ##         save_in_file saves the result in a file named data/norm_conductor_degree_size_nbtests_sampling.sage
    ##         short_output is the same as in test_prime_disc
    ## Output: a list of length nb_tests, each element of the list is list of the form [f,reduced_norm,time] (or [reduced_norm, time] if short_output = True), where 
    ##         - reduced_norm is either N(C_O)^(1/[K/QQ]) for C_O the conductor of the order O = ZZ[X]/f, in the field K = QQ[X]/f, or None if the time needed to perform the computation was larger than max_time_per_test
    ##        - time is the time needed to perform the computation (it is equal to max_time_per_test when the computation took too much time and was cancelled)
    if fixed_randomness:
        set_random_seed(42)
    if save_in_file:
        filename = "data/norm_conductor_"+str(degree)+"_"+str(size_poly)+"_"+str(nb_tests)+"_"+sampling_method+".sage"
        f_data = open(filename,"w")
        f_data.write("(nb_tests, degree, size_poly, sampling_method, max_time_per_test) = "+str((nb_tests, degree, size_poly, sampling_method, max_time_per_test)))
        f_data.write("\n\nres = [")
    
    res = []
    failed_poly = 0    
    
    for test in range(nb_tests):
        f = multi_sampling(degree,size_poly, sampling_method)
        assert(f.is_monic())
        n = f.degree()
        K.<a> = NumberField(f)
        if verbose:
            print("\nf = ", f)
        ## Tries to compute the discriminant of K, which requires factoring df (and may be quite long)
        ## abort it is takes more than max_time_per_test seconds
        try:
            t1 = time.time()
            signal.alarm(max_time_per_test) ## set an alarm for max_time_per_test seconds
            disc = K.discriminant()
            signal.alarm(0) ## disable the alarm if the computation finished on time
            t2 = time.time()
            if verbose:
                print("time:", Rshort(t2-t1))
        except AlarmInterrupt:
            if verbose:
                print("Unable to compute the discriminant, moving on...")
            failed_poly += 1
            if short_output:
                tmp_res = [None, max_time_per_test]
            else:
                tmp_res = [f,None, max_time_per_test]
            res += [tmp_res] ## failed to compute discriminant
            if save_in_file:
                f_data.write(str(tmp_res)+",\n")
                f_data.close()
            continue
         
        norm_CO = Rshort(norm_conductor_ideal(f,a,K, disc = disc)^(1/n))
        if verbose:
            print("N(C_O)^(1/n) = ", norm_CO)
            
        if short_output:
            tmp_res = [norm_CO, t2-t1]
        else:
            tmp_res = [f,norm_CO, t2-t1]
            
        res += [tmp_res] ## if it succeeds, add norm_CO to the data about f
        if save_in_file:
            f_data.write(str(tmp_res)+",\n")
      
        
        
    #print("Proportion of failure:", Rshort(failed_poly/nb_tests))
    if save_in_file:
        f_data.write("]")
        f_data.close()
    return res


def generate_and_verify_monogenic(degree, size_poly, sampling_method="unif", max_attempts=100, max_time_per_test=10, verbose=False):
    """
    Generate random polynomials and check if they are monogenic.
    
    Input:
        degree - degree of the polynomial to generate
        size_poly - bound on the coefficients
        sampling_method - "unif", "mignotte", or "resultant"
        max_attempts - maximum number of polynomials to try before giving up
        max_time_per_test - maximum time (in seconds) to compute discriminant
        verbose - whether to print detailed information
    
    Output:
        A tuple (f, K, is_monogenic, norm_conductor) where:
        - f is the polynomial
        - K is the number field (if successfully created)
        - is_monogenic is True if the field is monogenic
        - norm_conductor is N(C_O)^(1/n) where C_O is the conductor ideal
    """
    
    for attempt in range(max_attempts):
        if verbose:
            print(f"\n{'='*60}")
            print(f"Attempt {attempt + 1}/{max_attempts}")
            print('='*60)
        
        # Use existing function to generate and test
        # tests_norm_conductor generates a polynomial and computes the reduced norm of conductor
        # It returns a list of [f, reduced_norm, time]
        results = tests_norm_conductor(1, degree, size_poly, sampling_method, 
                                     max_time_per_test=max_time_per_test, 
                                     save_in_file=False,
                                     verbose=verbose)
        
        if not results:
             continue
             
        # Extract result
        entry = results[0]
        # entry is [f, reduced_norm, time]
        f = entry[0]
        reduced_norm = entry[1]
        
        if reduced_norm is None:
             if verbose:
                 print("Failed to compute discriminant/conductor (Timeout or Error)")
             continue
             
        # Check if monogenic (conductor norm = 1 => reduced_norm = 1)
        # Note: reduced_norm is a float/real, so we check for closeness to 1 or exact 1
        is_monogenic = (abs(reduced_norm - 1) < 1e-6)
        
        if is_monogenic:
            if verbose:
                print(f"\n{'*'*60}")
                print(f"SUCCESS! Found a monogenic polynomial!")
                print(f"{'*'*60}")
                print(f"f = {f}")
            
            # Recreate Number Field for verification
            # Note: discriminant is not cached since K is new, but we don't strictly need it 
            # for all verification steps (only Step 4, 5 check against Delta_K)
            # We can recompute it (expensive) or assume Delta_K = Disc(f) since it's monogenic.
            # verify_equation_derivation will call K.discriminant()
            
            try:
                K = NumberField(f, 'a')
                return (f, K, True, reduced_norm)
            except Exception as e:
                if verbose:
                    print(f"Error recreating number field: {e}")
                continue
        else:
            if verbose:
                print(f"Polynomial is NOT monogenic (N(C_O)^(1/n) = {reduced_norm})")
    
    if verbose:
        print(f"\nFailed to find a monogenic polynomial after {max_attempts} attempts")
    
    return (None, None, False, None)

    
################
## How to use ##
################

#print("\n\nTesting whether disc(f) is prime for 1000 random polynomials")
#res2 = tests_prime_disc(1000,30,10)
#print("\n\nComputing N(C_O)^{1/n} for 10 random polynomials")
#res3 = tests_norm_conductor(10,20,5)
