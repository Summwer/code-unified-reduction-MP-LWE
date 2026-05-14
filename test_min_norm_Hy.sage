#load("test_sigma_beta.sage")
load("test_Hankel_matrix.sage")
load("test_monogenic.sage")
load("test_sigma_beta_equation.sage")

from sage.modules.free_module_integer import IntegerLattice


def generate_sigma_beta_different(f, CCprec = None):
    """
    Build the complex embedding matrices for cond_O and β·cond_O,
    where cond_O = {x in K | x·O_K ⊆ O}, O = Z[a] is the order generated
    by a root a of f (need not be monogenic).

    Returns:
        sigma_matrix_conductor.T  - d x d complex matrix for σ(cond_O)
        sigma_matrix_beta_conductor.T - d x d complex matrix for σ(β·cond_O)
    """
    if(CCprec is None):
        CCprec = ComplexField()
    d = f.degree()
    K = NumberField(f, 'a')
    alphas = f.roots(ring = CCprec, multiplicities = False)
    target_values = [max(RR(1), abs(alpha)) for alpha in alphas]

    # Conductor ideal basis: cond_O = {x in K | x·O_K ⊆ O}
    conductor_basis = compute_conductor_basis(K)
    # σ(cond_O)
    sigma_matrix_conductor = compute_sigma_matrix(conductor_basis, K.embeddings(CCprec), CCprec)
    

    # Matrix σ(β·cond_O)
    rows = []
    for i in range(d):
        # σ(β·cond_O) = diag(σ(β))·σ(cond_O)
        factor = target_values[i]^(d-1)
        row = [factor * x for x in sigma_matrix_conductor.row(i)]
        rows.append(row)
    
    sigma_matrix_beta_conductor = matrix(CCprec, rows)
    
    
    #Ensure the basis should be in row
    return sigma_matrix_conductor.transpose(), sigma_matrix_beta_conductor.transpose()

def complex_matrix_to_zoom_vector_matrix(M, zoom = 1, prec = 100):
    """
    Convert a d x d complex matrix to a d x 2d zoom matrix.
    Each element (a + bi) is expanded to (a, b).
    
    Input:
        M - complex matrix of size rows x cols
        zoom - the zoom to scale the matrix
        
    Output:
        M_zoom - zoom matrix of size rows x (2*cols)
    """
    rows = M.nrows()
    cols = M.ncols()
    zoom_rows = []
    
    # Get the base ring's precision if possible, else default to RR
    try:
        R = RealField(prec)
    except:
        R = RR
    
    zero_cols = [True] * cols *2
    for i in range(rows):
        new_row = []
        for j in range(cols):
            z = M[i, j]
            
            zoom_real = R(z.real())*2**(-zoom)
            zoom_imag = R(z.imag())*2**(-zoom)
            assert(zoom_real != float("inf"))
            assert(zoom_imag != float("inf"))
            # Preserve precision: don't convert to float
            if( abs(zoom_real) >= 1/2.):
                zero_cols[j*2] = False
            if( abs(zoom_imag) >= 1/2.):
                zero_cols[j*2 + 1] = False
            new_row.append(ZZ(round(zoom_real)))
            new_row.append(ZZ(round(zoom_imag)))
        zoom_rows.append(new_row)
    
    non_zero_col_matrix = []
    for i in range(rows):
        new_row = []
        for j in range(cols*2):
            if(not zero_cols[j]):
                new_row.append(zoom_rows[i][j])
        non_zero_col_matrix.append(new_row)
    
    assert(len(non_zero_col_matrix[0]) >= rows)
    return matrix(ZZ, non_zero_col_matrix)


def check_matrix_value(Hy, sv, debug = False):
    flag = True
    Wrong_list = []
    for i in range(Hy.nrows()):
        
        if(abs(max(abs(Hy[i][0]),abs(Hy[i][-1]))-abs(sv[i]))>1): 
            flag = False
            Wrong_list.append(i)
        
    if(len(Wrong_list)>0 and debug):
        print("wrong_list = ", Wrong_list)
    return flag

def gaussian_heuristic(B, precision = 100):
    # Check if matrix is complex
    try:
        is_complex = 'Complex' in str(B.base_ring()) or 'CC' in str(B.base_ring())
    except:
        is_complex = False

    if is_complex:
        # Expand complex d x d to real d x 2d without rounding (keep precision)
        rows = B.nrows()
        cols = B.ncols()
        real_rows = []
        # Use existing precision or default to sufficiently high
        R_pf = RealField(precision) 
        
        for i in range(rows):
            new_row = []
            for j in range(cols):
                z = B[i, j]
                new_row.append(R_pf(z.real()))
                new_row.append(R_pf(z.imag()))
            real_rows.append(new_row)
        B = matrix(R_pf, real_rows)

    dim = B.nrows()
    # For a lattice of rank 'dim' embedded in a space (cols >= dim), volume is sqrt(det(BB^T))
    # If B is square real matrix, this is just abs(det(B))
    
    if(B.nrows() < B.ncols()):
        det_sq = abs((B*B.transpose()).determinant())
        GH = sqrt(dim/(2*pi*e)) * ( det_sq ** (1.0/(2*dim)))
    else:
        det_val = abs(B.determinant())
        GH = sqrt(dim/(2*pi*e)) * ( det_val ** (1.0/dim))
    return round(GH,2), round(GH*sqrt(2*pi*e),2)



def formatted_vector_string(vec):
    # Format each element as a + bi with 2 decimal places
    # Handle real/complex types appropriately
    elements = []
    for z in vec:
         real_part = float(z.real())
         imag_part = float(z.imag())
         if imag_part > 0:
             elements.append(f"{real_part:.2f} + {imag_part:.2f}I")
         elif imag_part < 0:
             elements.append(f"{real_part:.2f} - {abs(imag_part):.2f}I")
         else:
             elements.append(f"{real_part:.2f}")
            
    return "(" + ", ".join(elements) + ")"

def formatted_matrix_string(M):
    rows = M.nrows()
    cols = M.ncols()
    formatted_rows = []
    for i in range(rows):
        # Don't wrap in list - we need strings for join()
        formatted_rows.append(formatted_vector_string(M.row(i)))
    return "\n".join(formatted_rows)

def compute_norm_max(M):
    d = M.nrows()
    norm_max = 0
    for i in range(d):
        for element in M.row(i):
            norm_max = max(norm_max, abs(element))
    return norm_max

def compute_min_norm_Hy(f, K, zoom = 1,  precision = 100, verbose = False, verification = False):
    """
    Compute the minimum norm of Vf^T * Σ(y) = M_{f,y^dual} * Vf^-1 
    """

    verify_failed = False 
    

    CCprec = ComplexField(precision)


    sigma_matrix_different, sigma_matrix_beta_different = generate_sigma_beta_different(f, CCprec = CCprec)

    # Convert to real matrix
    B = complex_matrix_to_zoom_vector_matrix(sigma_matrix_beta_different, zoom = zoom, prec = precision)

    d = f.degree()
 

    B_ = B.LLL( precision = precision)

    for block_size in range(10, d+1):
        B_= B_.BKZ(block_size=block_size, precision = precision)


    for i in range(d):
        if(vector(B_[i]) != vector([0]*B_.ncols()) ):
            v = B_[i]
            break
  
    c = B.solve_left(v)
    assert(vector(c*B) == v)
    
    c = vector(CCprec,c)
    
    sv = c * sigma_matrix_beta_different

    #H_y = Vf^T * Sigma(y)
    sigma_y = vector(c * sigma_matrix_different)
    #print("σ(y) = ", formatted_vector_string(sigma_y))
    Sigma_y = diagonal_matrix(CCprec, sigma_y)
    
    Vf_matrix = Vf(f, precision = precision)
    #print("Vf_matrix = ", formatted_matrix_string(Vf_matrix))
        
    # Matrix Hy = Vf^T * Σ(y)
    Hy = Vf_matrix.transpose() * Sigma_y
    Hy = Hy.transpose() #Hy^T
    norm_max_Hy = compute_norm_max(Hy)
    if(not verification):
        return round(norm_max_Hy,2), None
    else:
        
        
        if(not check_matrix_value(Hy, sv)):
            verify_failed = True
            
            
        
        gh_zoom, minkowski_zoom = gaussian_heuristic(B, precision = precision)
        
        gh_sigma_beta, minkowski = gaussian_heuristic(sigma_matrix_beta_different, precision = precision)
        
        # print("v.norm() = ", round(v.norm(),2), ", GH(zoom_matrix) = ", gh_zoom, ", minkowski_zoom = ", minkowski_zoom)

        # print("inf norm(Vf^T * Σ(y)) = ", round(norm_max_Hy,2), ", sv.norm() = ", round(sv.norm(),2), ", GH(σ(β·CO) ) = ", gh_sigma_beta, ", minkowski_bound(σ(β·CO) ) = ", minkowski)

                
        # Verify scaling relationship: v is on zoomed lattice, sv is on original lattice
        # Expected: gh_zoom/gh_sigma_beta ≈ zoom (due to scaling)

        
        if(gh_sigma_beta == 0 or minkowski == 0.):
            verify_failed = True
        #elif( not (abs(round(gh_zoom/gh_sigma_beta) - zoom) < 2 and abs(round(minkowski_zoom/minkowski) - zoom) < 2)):
        #    verify_failed = True
      
        
        
        if(sv.norm()>minkowski or norm_max_Hy > sv.norm()):
            verify_failed = True
        
        if(false):
            print("sv.norm() < gh_sigma_beta_sigma_beta? ", sv.norm() < gh_sigma_beta)
            print("sv.norm() < minkowski_bound(σ(β·CO) )? ", sv.norm() < minkowski)
            print("inf norm(Vf^T * Σ(y)) < sv.norm()? ", norm_max_Hy < sv.norm())
            
            print("basis = ", formatted_matrix_string(sigma_matrix_beta_different))
            print("basis (zoom before) = ", formatted_matrix_string(B))
            print("basis (zoom after) = ", formatted_matrix_string(B_))
            print("norm_max_Hy = ", round(norm_max_Hy,2))
            print("gh_sigma_beta = ", round(gh_sigma_beta,2))
            print("minkowski = ", round(minkowski,2))
            print("gh_zoom = ", round(gh_zoom,2))
            print("minkowski_zoom = ", round(minkowski_zoom,2))
                
        return round(norm_max_Hy,2), verify_failed


def norm_upperbound(d, size_poly):
    """
    Compute the norm upper bound of Vf^T * Σ(y) = M_{f,y^dual} * Vf^-1 
    """
    return 16*d^7 * size_poly^5 


def tests_min_norm_Hy(nb_tests, degree, size_poly, precision = 1000, zoom = 1, save_in_file = False):
    """
    Test the minimum norm of Hy for `nb_tests` random monic irreducible
    polynomials of the given degree and coefficient size.
    f need NOT be monogenic.
    """
    if save_in_file:
        filename = "data/norm_Hy"+str(degree)+"_"+str(size_poly)+"_"+str(nb_tests)+".sage"
        f_data = open(filename,"w")
        f_data.write("(nb_tests, degree, size_poly) = "+str((nb_tests, degree, size_poly)))
        f_data.write("\n\nres = [")

    verify_equation_false = False
    t = 0
    fs = []
    while(t < nb_tests):
        # Generate a random monic irreducible polynomial (no monogenic requirement)
        f = unif_poly(degree, B = size_poly)
        if f in fs:
            continue
        fs.append(f)
        K = NumberField(f, 'a')
        norm_max_Hy, verify_failed = compute_min_norm_Hy(f, K, precision = precision, zoom = zoom, verification = True)
        if save_in_file:
            f_data.write(str(norm_max_Hy)+",")
    
        if(verify_failed):
            verify_equation_false = True
            break
        t += 1
    if save_in_file:
        f_data.write("]")
        f_data.close()
    return verify_equation_false



'''
if __name__ == "__main__":
    test_times = 100
    degree = 5
    size_poly = 20
    precision = 1e5
    zoom = 1
    CCprec = ComplexField(precision)
    
    flag = True
    t = 0
    success_times = 0
    fs = []
    while(t < test_times):
        is_monogenic = False
        while(True):
            f, K, is_monogenic, norm_cond = generate_and_verify_monogenic(degree=degree, size_poly=size_poly, sampling_method="unif")
            if(f not in fs and is_monogenic):
                fs.append(f)
                break
            else:
                continue
     

        
        sigma_matrix_different, sigma_matrix_beta_different = generate_sigma_beta_different(f, CCprec = CCprec)

        # Convert to real matrix
        B = complex_matrix_to_zoom_vector_matrix(sigma_matrix_beta_different, zoom = zoom, prec = precision)
    
        norm_max_Hy = compute_min_norm(B, sigma_matrix_different, sigma_matrix_beta_different, zoom, degree, f,  CCprec = ComplexField(precision), precision = precision)    
        UB = norm_upperbound(degree, size_poly)
        print(norm_max_Hy, UB, norm_max_Hy<UB)

        t += 1
'''    