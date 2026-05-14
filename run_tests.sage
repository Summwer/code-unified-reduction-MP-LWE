###############################################################
## Uses the functions defined in the files                   ##
## test_***.sage and performs the tests that                 ##
## are used in the article                                   ##
## Running all the test on a singne core                     ##
## may take quite some time, it is better to                 ##
## run the following code on a machine with multiple cores   ##
## and update the file below with the number of cores        ##
## available                                                 ##
###############################################################

## Note: loading this file will not immediately run the tests.
## The command lines actually running the computation have been commented
## below (you should uncomment the one you want to run the test)

## Requirement: when running the tests, you should execute the code in a repository containing a sub-repository named "data" (the results of the tests will be saved there)

from multiprocessing import Pool, set_start_method
import sys

# Fix for macOS: force fork instead of spawn
try:
    set_start_method('fork')
except RuntimeError:
    pass  # start method already set

## load files containing the functions
load("test_Vandermonde.sage")
load("test_monogenic.sage")
load("test_Hankel_matrix.sage")
load("test_sigma_beta_equation.sage")
load("test_min_norm_Hy.sage")

## Function for parallel runs
def execute_task_parallel(f,l_input,nb_threads):
    ## runs function f on all the inputs stored in l_input in parallel with nb_threads threads
    with Pool(nb_threads) as p:
        return p.map(f, l_input)

    
###########################################
#### Running tests for the Vandermonde ####
###########################################
    

## tuples for Mignotte polynomials
list_input_Vandermonde_mignotte = []
nb_tests = 10
for degree in [10,50,100,150,200]:
    for B in [10**2,10**3,10**4,10**5]:
        list_input_Vandermonde_mignotte += [(degree,B,nb_tests,"mignotte")]
        
def one_Vandermonde(x):
    ## Input: x = (degree,B,nb_tests), as in list_input_Vandermonde
    ## Output: execute tests_Vandermonde on x with option saved_in_file and fixed_randomness
    (degree,B,nb_tests,sampling_method) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...")
    t1 = time.time()
    tests_Vandermonde(nb_tests, degree,B, sampling_method = sampling_method, save_in_file = True, short_output = True, fixed_randomness = True)
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    
  
def run_Vandermonde(list_input_Vandermonde, nb_threads = 3):
    if list_input_Vandermonde[0][-1] == "unif":
        print("Running tests for Vandermonde norm with uniform polynomials...")
    elif list_input_Vandermonde[0][-1] == "mignotte":
        print("Running tests for Vandermonde norm with Bugeaud-Mignotte polynomials...")
    else:
        print("invalid input in run_Df")
        return     
    execute_task_parallel(one_Vandermonde,list_input_Vandermonde, nb_threads)
    print("...done tests for Vandermonde matrices")
    

    
    
#######################################################
#### Running tests for testing prime discriminants ####
#######################################################
        
def one_prime_disc(x):
    ## Input: x = (degree,B,nb_tests), as in list_input_prime_disc
    ## Output: execute tests_prime_disc on x with option saved_in_file and fixed_randomness
    (degree,B,nb_tests) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...")
    t1 = time.time()
    tests_prime_disc(nb_tests, degree,B, save_in_file = True, short_output = True, fixed_randomness = True)
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    
  
def run_prime_disc(list_input_prime_disc, nb_threads = 5):
    print("Running tests for testing prime discriminants...")      
    execute_task_parallel(one_prime_disc,list_input_prime_disc,nb_threads)
    print("...done tests for testing prime discriminants")
    

    
######################################
#### Running tests for conductors ####
######################################
  
def one_conductor(x):
    ## Input: x = (degree,B,nb_tests), as in list_input_conductor
    ## Output: execute tests_conductor on x with option saved_in_file and fixed_randomness

    (degree,B,nb_tests, max_time) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...\n")
    t1 = time.time()
    tests_norm_conductor(nb_tests, degree,B, max_time_per_test = max_time, save_in_file = True, short_output = True, fixed_randomness = True)
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    
  
def run_conductor(list_input_conductor, nb_threads = 5):
    print("Running tests for conductors...")      
    execute_task_parallel(one_conductor,list_input_conductor,nb_threads)
    print("...done tests for conductors")
    
###############################################################
#### Running tests of verifiction for Hankel matrix:       ####
#### M_{f,y^dual} * a = Rot_f(a) * (y1^dual, ..., yd^dual) ####
#### for radnom monogenic f, random a and random y         ####
###############################################################
def one_Hankel_matrix(x):
    ## Input: x = (degree,B,nb_tests), as in list_input_Hankel_matrix
    ## Output: execute verifications_Hankel_matrix on x

    (degree,B,nb_tests, precision) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...\n")
    t1 = time.time()
    verify_Hankel_matrix_failed = verifications_Hankel_matrix(nb_tests, degree,B, precision = precision)
    if(not verify_Hankel_matrix_failed):  
        print("All tests passed for degree ", degree, ", with bound ", B, ", precision = 2^", precision)
    else:
        print("Some tests failed for degree ", degree, ", with bound ", B, ", precision = 2^", precision)
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    return verify_Hankel_matrix_failed
    

def run_Hankel_matrix(list_input_Hankel_matrix, nb_threads = 5):
    print("Running tests for Hankel matrix...")      
    print("#### Check M_{f,y^dual} * a = Rot_f(a) * (y1^dual, ..., yd^dual) ####")
    verifications_failed = execute_task_parallel(one_Hankel_matrix,list_input_Hankel_matrix,nb_threads)
    if(any(verifications_failed)):
        print("Some verifications failed for Hankel matrix")
    else:
        print("All verifications passed for Hankel matrix")
    print("...done tests for Hankel matrix")


###################################################################################################
#### Running tests of verifiction for the equation derivation in strict order:                 ####
####     1. Check |Δ_K| = |det(V)|^2                                                           ####
####     2. Check det(σ((O_K^∨)^{-1})) = |Δ_K|^(3/2)                                           ####
####     3. Check det(σ(β·(O_K^∨)^{-1})) = |N(β)|·det(σ((O_K^∨)^{-1})) = M(f)^(d-1)·|det(V)|^3 ####
###################################################################################################

def one_sigma_beta_equation(x):
    ## Input: x = (degree,B,nb_tests, precision), as in list_input_sigma_beta_equation
    ## Output: execute verifications_sigma_beta_equation on x

    (degree,B,nb_tests, precision) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...\n")
    t1 = time.time()
    verify_equation_failed = verifications_sigma_beta_equation(nb_tests, degree, B, precision = precision)
    if(not verify_equation_failed):  
        print("All tests passed for degree ", degree, ", with bound ", B, ", precision = 2^", precision)
    else:
        print("Some tests failed for degree ", degree, ", with bound ", B, ", precision = 2^", precision)
        raise ""
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    return verify_equation_failed

def run_sigma_beta_equation(list_input_sigma_beta_equation, nb_threads = 5):
    print("Running verifications for sigma beta equation...")
    print("1. Check |Δ_K| = |det(V)|^2                                                           ")
    print("2. Check det(σ((O_K^∨)^{-1})) = |Δ_K|^(3/2)                                          ")
    print("3. Check det(σ(β·(O_K^∨)^{-1})) = |N(β)|·det(σ((O_K^∨)^{-1})) = M(f)^(d-1)·|det(V)|^3 ")
    verify_equations_failed = execute_task_parallel(one_sigma_beta_equation,list_input_sigma_beta_equation,nb_threads)
    if(any(verify_equations_failed)):
        print("Some verifications failed for sigma beta equation")
    else:
        print("All verifications passed for sigma beta equation")
    print("...done verifications for sigma beta equation")


###################################################################################################
####  Find the minimal norm of M_{f,y} * Vf^-1 = Vf^T * Sigma(y), for any y in (O_K^∨)^-1.     ####
###################################################################################################



def one_min_norm_Hy(x):
    (degree,B,nb_tests, precision, zoom) = x
    print("\nTesting degree ", degree, ", with bound ", B, "for ", nb_tests, "tests...\n")
    t1 = time.time()
    verify_equation_failed = tests_min_norm_Hy(nb_tests, degree, B, precision = precision, zoom = zoom, save_in_file = True)
    if(not verify_equation_failed):  
        print("All tests passed for degree ", degree, ", with bound ", B, ", precision = 2^", precision, "and zoom = 2^-", zoom)
    else:
        print("Some tests failed for degree ", degree, ", with bound ", B, ", precision = 2^", precision, "and zoom = 2^-", zoom)
        
    print("Time taken: ", Rshort(time.time()-t1), " seconds")
    print("...this was degree ", degree, ", with bound ", B)
    return verify_equation_failed

def run_min_norm_Hy(list_input_min_norm_Hy, nb_threads = 5):
    print("Running verifications for min norm Hy...")
    #for i in [-1, -2 ]:
    #    one_min_norm_Hy(list_input_min_norm_Hy[i])
    verify_equations_failed = execute_task_parallel(one_min_norm_Hy,list_input_min_norm_Hy,nb_threads)
    
    if(any(verify_equations_failed)):
        print("Some verifications failed for min norm Hy")
    else:
        print("All verifications passed for min norm Hy")
    print("...done verifications for min norm Hy")
    

if __name__ == '__main__':

    nb_threads = 30

    ## triples of the form (degree, bound on coefficients, nb of tests) for uniform polynomials
    list_input_Vandermonde_unif = []

    for degree in range(50,501,50):
        list_B = [10,10**2,10**3]
        nb_tests = 20
        if degree <= 400:
            nb_tests = 50
        if degree <= 300:
            nb_tests = 100
        if degree <= 100:
            list_B += [10**4,10**5]
            nb_tests = 1000
        for B in list_B:
            list_input_Vandermonde_unif += [(degree,B,nb_tests,"unif")]
        

    ## TODO Uncomment here to run the test for computing the norm of the Vandermonde matrix
    ## on random and on Mignotte polynomials
    ## Update the quantity nb_threads to what is available to you


    #run_Vandermonde(list_input_Vandermonde_unif, nb_threads = nb_threads)
    #run_Vandermonde(list_input_Vandermonde_mignotte, nb_threads = nb_threads)

    ## triples of the form (degree, bound on coefficients, precision, nb of tests)
    list_input_prime_disc = []

    for degree in range(50,501,50):
        if degree <= 200:
            list_B = [2,5,10,100,1000]
            nb_tests = 100000
        elif degree <= 400:
            list_B = [2,5,10]
            nb_tests = 50000
        else:
            list_B = [2,5]
            nb_tests = 50000
        for B in list_B:
            list_input_prime_disc += [(degree,B,nb_tests)]

    ## TODO Uncomment here to run the test for testing when the discriminant
    ## of a random polynomial is prime
    ## Update the quantity nb_threads to what is available to you


    #run_prime_disc(list_input_prime_disc, nb_threads = nb_threads)



    ## triples of the form (degree, bound on coefficients, precision, nb of tests)
    list_input_conductor = []

    for degree in range(5,41,5):
        if degree <= 10:
            nb_tests = 1000
            max_time = 100
            list_B = [2,5,10,50,100,1000]
        elif degree <= 20:
            nb_tests = 100
            max_time = 600
            list_B = [2,5,10,50]
        elif degree <= 25:
            nb_tests = 50
            max_time = 600
            list_B = [2,5,10]
        else:
            nb_tests = 20
            max_time = 1200
            list_B = [2,5]
        for B in list_B:
            list_input_conductor += [(degree,B,nb_tests, max_time)]

    ## TODO Uncomment here to run the test for computing the conductor ideal
    ## of a random polynomial
    ## Update the quantity nb_threads to what is available to you


    #run_conductor(list_input_conductor, nb_threads = nb_threads)



    ## triples of the form (degree, bound on coefficients, precision, nb of tests)
    list_input_Hankel_matrix = []

    for degree in range(5,41,5):
        if degree <= 10:
            nb_tests = 1000
            list_prec = [100, 100, 1e3, 1e3, 1e4, 1e4]
            list_B = [2,5,10,50,100,1000]
        elif degree <= 20:
            nb_tests = 100
            list_prec = [100, 1e3, 1e3, 1e3]
            list_B = [2,5,10,50]
        elif degree <= 25:
            nb_tests = 50
            list_prec = [1e3, 1e3, 1e3]
            list_B = [2,5,10]
        else:
            nb_tests = 20
            list_prec = [1e3, 1e3]
            list_B = [2,5]
        for i in range(len(list_B)):
            list_input_Hankel_matrix += [(degree,list_B[i],nb_tests, list_prec[i])]

    ## TODO Uncomment here to run the test for verfying the Hankel matrix
    ## M_{f,y^dual} * a = Rot_f(a) * (y1^dual, ..., yd^dual)
    ## for random monogenic f, random a and random y. 

   
    #run_Hankel_matrix(list_input_Hankel_matrix, nb_threads = nb_threads)

  


    ## triples of the form (degree, bound on coefficients, precision, nb of tests)
    list_input_sigma_beta_equation = []

    for degree in range(5,41,5):
        if degree <= 10:
            nb_tests = 1000
            list_prec = [1e3, 1e3, 1e3, 1e3, 1e4, 1e4]
            list_B = [2,5,10,50,100,1000]
        elif degree <= 20:
            nb_tests = 100
            list_prec = [1e3, 1e4, 1e4, 1e4]
            list_B = [2,5,10,50]
        elif degree <= 25:
            nb_tests = 50
            list_prec = [1e4, 1e4, 1e4]
            list_B = [2,5,10]
        else:
            nb_tests = 20
            list_prec = [1e4, 1e4]
            list_B = [2,5]
        for i in range(len(list_B)):
            list_input_sigma_beta_equation += [(degree,list_B[i],nb_tests, list_prec[i])]

    ## TODO Uncomment here to run the test for verfying the Hankel matrix
    ## M_{f,y^dual} * a = Rot_f(a) * (y1^dual, ..., yd^dual)
    ## for random monogenic f, random a and random y. 

    
    # run_sigma_beta_equation(list_input_sigma_beta_equation, nb_threads = nb_threads)

    
    ## triples of the form (degree, bound on coefficients, precision, zoom, nb of tests)
    list_input_min_norm_Hy = []

   

    #float infinity or verify failed: increase precision
    #large zoom can decrease the time of lattice reduction
    #First setting the precision to ensure the correctness, then increase the zoom to decrease the time of lattice reduction. 
    
    for degree in range(5,11):
        if degree <= 5:
            nb_tests = 30
            list_prec = [1e3, 1e5, 1e6, 1e6, 1e6]  
            list_zoom =  [0]*5 
            list_B = [2, 9, 16, 23, 30] 
        elif degree <= 10:
            nb_tests = 30
            list_prec = [1e6]
            list_zoom = [0]
            list_B = [2]
        for i in range(len(list_B)):
            list_input_min_norm_Hy += [(degree,list_B[i],nb_tests, list_prec[i], list_zoom[i])]

  

    print(list_input_min_norm_Hy)
    ## TODO Uncomment here to run the test for verfying the Hankel matrix
    ## min_{y} ||HVfinv(y)||_inf = ||f||_inf
    ## for random monogenic f, random a and random y. 

    
    #run_min_norm_Hy(list_input_min_norm_Hy, nb_threads = nb_threads)



    
