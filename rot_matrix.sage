def rot_matrix(f, a, d=None):
    """
    Return the matrix Rot_f^d(a) as defined in the literature.

    For a polynomial f of degree m, integer d > 0, and polynomial a,
    Rot_f^d(a) is the d x m matrix whose i-th row contains the coefficients
    of the polynomial (x^(i-1) * a) mod f.

    Parameters
    ----------
    f : polynomial
        The modulus polynomial, of degree m.
    a : polynomial
        The polynomial a, in the same polynomial ring as f.
    d : integer (optional)
        The number of rows. If not provided, defaults to degree of f.

    Returns
    -------
    matrix
        The d x m matrix over the base ring of f.

    Examples
    --------
    sage: R.<x> = QQ[]
    sage: f = x^3 + 1
    sage: a = x^2 + 2*x + 3
    sage: rot_matrix(f, a, 3)
    [ 3  2  1]
    [-1  3  2]
    [-2 -1  3]

    For the special case f = x^d + 1 and a = sum_{i=0}^{d-1} a_i x^i,
    the matrix is the circulant matrix with first row [a0, a1, ..., a_{d-1}].
    """
    R = f.parent()
    x = R.gen()
    m = f.degree()
    if d is None:
        d = m
    # Ensure a is in the same ring
    a = R(a)
    # Initialize matrix
    base_ring = R.base_ring()
    M = matrix(base_ring, d, m)
    for i in range(d):
        # Compute p = (x^i * a) mod f
        p = (x**i * a) % f
        # Get coefficients up to degree m-1
        coeffs = p.coefficients(sparse=False)
        # Pad with zeros if needed
        if len(coeffs) < m:
            coeffs += [base_ring(0)] * (m - len(coeffs))
        # Set row i
        for j in range(m):
            M[i, j] = coeffs[j]
    return M

