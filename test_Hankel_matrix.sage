load("random_polynomials.sage")
load("compute_dual_basis.sage")
load("test_Vandermonde.sage")
load("rot_matrix.sage")
load("test_monogenic.sage")

def sigma(alpha, y):
    """
    Compute the image of y under the embedding corresponding to root alpha.
    
    Parameters
    ----------
    alpha : complex number
        A root of the defining polynomial.
    y : element of a number field
        The element to embed.
        
    Returns
    -------
    complex number
        The image of y under the embedding.
    """
    K = y.parent()
    f = K.defining_polynomial()
    R = f.parent()
    a = R.gen()
    
    # Express y as a polynomial in a
    coeffs = [y.polynomial()[i] for i in range(f.degree())]
    
    # Evaluate at alpha
    image = sum(c * (alpha^i) for i, c in enumerate(coeffs))
    
    return image


def Sigma(y, precision = 100):
    """
    Compute the images of y under all embeddings of the number field.
    
    Parameters
    ----------
    y : element of a number field
        The element to embed.
        
    Returns
    -------
    matrix
        Diagonal matrix with diagonal entries being sigma(alpha, y) for each root alpha.
    """
    K = y.parent()
    f = K.defining_polynomial()
    R = f.parent()
    
    # Get roots of the defining polynomial
    CCprec = ComplexField(precision)  # Use high precision for roots
    roots = f.roots(ring=CCprec, multiplicities=False)
    
    # Compute images
    images = [sigma(alpha, y) for alpha in roots]
    
    # Create diagonal matrix
    diag_matrix = diagonal_matrix(CCprec, images)
    
    return diag_matrix


def verifications_Hankel_matrix(nb_tests, degree,size_poly,sampling_method="unif", fixed_randomness = False, precision = 1e4):
    '''
    Verification function:  M_{f,y^dual} * a = Rot_f(a) * (y1^dual, ..., yd^dual) for radnom monogenic f, random a and random y  

    input: nb_tests - number of test times
           degree - degree of the defining polynomial
           size_poly - Upper bound of the absolute values of the coefficients of the polynomial f
           sampling_method - sampling method for the defining polynomial ("unif" or "gauss")
           fixed_randomness - whether to fix the randomness
           precision - precision of complex field. The low precision may cause the verification to fail. 
    
    output: the verification result
    '''

    if fixed_randomness:
        set_random_seed(42)

    res = []    
    CCprec = ComplexField(precision)
    verify_Hankel_matrix_failed = False
    R.<x> = QQ[]
    for test in range(nb_tests):

        is_monogenic = False
        while(not is_monogenic):
            f, K, is_monogenic, norm_cond = generate_and_verify_monogenic(degree=degree, size_poly=size_poly, sampling_method="unif")

        K, basis, dual = dual_basis_from_polynomial(f)
        y = K(random_field_element(f))
        
        y_dual_coordinates, expression = express_in_dual_basis(y, dual)

        Vf_matrix = Vf(f, precision = precision)
        H_y = Vf_matrix.transpose() * Sigma(y, precision = precision) * Vf_matrix

        a = random_field_element(f)

        Rot_matrix = rot_matrix(f, a)
        Rot_product = Rot_matrix * matrix(CCprec, y_dual_coordinates).transpose()
        Hy_product = H_y * matrix(CCprec,[[a[i] for i in range(f.degree())]]).transpose()

        diff = matrix(Hy_product) - matrix(Rot_product)

 
        epsilon = 1./precision
        if diff.norm() > epsilon:
            print( matrix(Hy_product))
            print()
            print( matrix(Rot_product))
            print("diff.norm() =", diff.norm())
            print(f)
            print("degree = ", degree)
            print("B = ", size_poly)
            print("precision =", precision)
            verify_Hankel_matrix_failed = True
            raise ""
    
    return verify_Hankel_matrix_failed

'''
if __name__ == "__main__":
    precision = 1000
    nb_tests, degree,size_poly = 100,10,20
    verify_Hankel_matrix_failed = verifications_Hankel_matrix(nb_tests, degree,size_poly, precision = precision)
    print("verify_Hankel_matrix_failed =", verify_Hankel_matrix_failed)
'''