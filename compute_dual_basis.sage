#!/usr/bin/env sage
# -*- coding: utf-8 -*-

"""
Compute the dual basis for a number field K = Q[x]/(f(x))
with standard basis (1, x, ..., x^{d-1}), where f is irreducible.

This implementation uses SageMath built-in functions for:
1. Creating number fields from irreducible polynomials
2. Computing traces of field elements
3. Matrix operations for the trace form
"""

def dual_basis_from_polynomial(f):
    """
    Compute the dual basis of (1, x, ..., x^{d-1}) for K = Q[x]/(f(x)).
    
    Parameters
    ----------
    f : polynomial over QQ
        An irreducible polynomial defining the number field.
        
    Returns
    -------
    tuple (K, basis, dual_basis)
        K : the number field K.<a> = NumberField(f)
        basis : list of standard basis elements [1, a, ..., a^{d-1}]
        dual_basis : list of dual basis elements
    """
    # Create the number field
    K.<a> = NumberField(f)
    d = K.degree()
    
    # Standard basis: 1, a, a^2, ..., a^{d-1}
    basis = [a^i for i in range(d)]
    
    # Compute the trace matrix M[i,j] = Tr(a^i * a^j)
    F = K.base_ring()
    M = matrix(F, d, d, [(a^(i+j)).trace() for i in range(d) for j in range(d)])
    
    # Inverse of the trace matrix
    M_inv = M.inverse()
    
    # Compute dual basis: beta_j = sum_{k=0}^{d-1} M_inv[k,j] * a^k
    dual_basis = [sum(M_inv[k, j] * a^k for k in range(d)) for j in range(d)]
    
    return K, basis, dual_basis


def dual_basis(K):
    """
    Compute the dual basis of (1, a, a^2, ..., a^{d-1}) for a given number field K.
    
    Parameters
    ----------
    K : NumberField
        A SageMath number field.
        
    Returns
    -------
    list
        The dual basis as a list of field elements.
    """
    d = K.degree()
    a = K.gen()
    F = K.base_ring()
    
    # Compute the trace matrix M[i,j] = Tr(a^i * a^j)
    M = matrix(F, d, d, [(a^(i+j)).trace() for i in range(d) for j in range(d)])
    
    # Inverse of the trace matrix
    M_inv = M.inverse()
    
    # Compute dual basis: beta_j = sum_{k=0}^{d-1} M_inv[k,j] * a^k
    dual_basis = [sum(M_inv[k, j] * a^k for k in range(d)) for j in range(d)]
    
    return dual_basis


def verify_dual_basis(K, basis, dual_basis):
    """
    Verify that Tr(b_i * b_j^*) = δ_ij for all i,j.
    
    Parameters
    ----------
    K : NumberField
        The number field.
    basis : list
        Standard basis.
    dual_basis : list
        Computed dual basis.
        
    Returns
    -------
    bool
        True if verification passes, False otherwise.
    """
    d = K.degree()
    for i in range(d):
        for j in range(d):
            tr = (basis[i] * dual_basis[j]).trace()
            expected = 1 if i == j else 0
            if tr != expected:
                return False
    return True


def coordinates_in_dual_basis(y, dual_basis):
    """
    Compute the coordinates of an element y in the dual basis.
    
    Parameters
    ----------
    y : element of a number field
        The element to express in the dual basis.
    dual_basis : list
        The dual basis of the number field.
        
    Returns
    -------
    list
        Coordinates c_i such that y = Σ_i c_i * dual_basis[i].
    """
    # The coordinates are given by c_i = Tr(y * b_i) where {b_i} is the standard basis
    # But we have the dual basis, not the standard basis.
    # Actually, if {β_i} is the dual basis, then for any y, y = Σ_i Tr(y * e_i) β_i
    # where {e_i} is the standard basis.
    # Alternatively, we can solve linear equations.
    
    # Get the number field
    K = y.parent()
    d = K.degree()
    
    # Standard basis: 1, a, a^2, ..., a^{d-1} where a = K.gen()
    a = K.gen()
    standard_basis = [a^i for i in range(d)]
    
    # Coordinates: c_i = Tr(y * standard_basis[i])
    coordinates = [(y * b).trace() for b in standard_basis]
    
    return coordinates


def express_in_dual_basis(y, dual_basis):
    """
    Express an element y as a linear combination of dual basis elements.
    
    Parameters
    ----------
    y : element of a number field
        The element to express.
    dual_basis : list
        The dual basis.
        
    Returns
    -------
    tuple (coordinates, expression)
        coordinates: list of coordinates
        expression: string representation of y as sum of coordinates * dual_basis[i]
    """
    coordinates = coordinates_in_dual_basis(y, dual_basis)
    
    # Build expression string
    terms = []
    for i, c in enumerate(coordinates):
        if c != 0:
            terms.append(f"({c})*β_{i}")
    
    if not terms:
        expression = "0"
    else:
        expression = " + ".join(terms)
    
    return coordinates, expression


# # Example usage
# if __name__ == "__main__":
#     print("Example: Computing dual basis and coordinates")
#     print("=" * 60)
    
#     # Define the polynomial
#     R.<x> = QQ[]
#     f = x^3 - 2
    
#     # Compute dual basis
#     K, basis, dual = dual_basis_from_polynomial(f)
    
#     print(f"Number field: {K}")
#     print(f"Defining polynomial: {f}")
#     print(f"Degree: {K.degree()}")
#     print(f"\nStandard basis: {basis}")
#     print(f"\nDual basis: {dual}")
    
#     # Verify
#     if verify_dual_basis(K, basis, dual):
#         print("\n✓ Verification passed: Tr(b_i * b_j^*) = δ_ij")
#     else:
#         print("\n✗ Verification failed")
    
#     # Test coordinate computation
#     print("\n" + "=" * 60)
#     print("Testing coordinate computation in dual basis:")
    
#     # Test with some elements
#     a = K.gen()
#     test_elements = [
#         K(1),          # Convert integer to field element
#         a,
#         a^2,
#         2*a + 3*a^2,
#         a^2 - a + 5
#     ]
    
#     for y in test_elements:
#         print(f"\nElement y = {y}")
#         coordinates, expression = express_in_dual_basis(y, dual)
#         print(f"  Coordinates: {coordinates}")
#         print(f"  Expression in dual basis: {expression}")
        
#         # Verify by reconstructing
#         reconstructed = sum(coordinates[i] * dual[i] for i in range(len(dual)))
#         if reconstructed == y:
#             print(f"  ✓ Reconstruction matches original")
#         else:
#             print(f"  ✗ Reconstruction error: {reconstructed} != {y}")
    
#     print("\n" + "=" * 60)
#     print("\nQuick test with another polynomial:")
    
#     # Another example
#     f2 = x^4 - 10*x^2 + 1
#     K2, basis2, dual2 = dual_basis_from_polynomial(f2)
    
#     print(f"Polynomial: {f2}")
#     print(f"Degree: {K2.degree()}")
    
#     if verify_dual_basis(K2, basis2, dual2):
#         print("✓ Dual basis is correct")
#     else:
#         print("✗ Dual basis is incorrect")
    
#     # Test coordinate computation for this field
#     a2 = K2.gen()
#     y_test = a2^3 - 2*a2 + 1
#     print(f"\nTest element in K2: y = {y_test}")
#     coordinates2, expression2 = express_in_dual_basis(y_test, dual2)
#     print(f"  Coordinates: {coordinates2}")
#     print(f"  Expression in dual basis: {expression2}")
    
#     # Verify reconstruction
#     reconstructed2 = sum(coordinates2[i] * dual2[i] for i in range(len(dual2)))
#     if reconstructed2 == y_test:
#         print(f"  ✓ Reconstruction matches original")
#     else:
#         print(f"  ✗ Reconstruction error: {reconstructed2} != {y_test}")
