import sympy as sp

def lagrange_interpolation(points):
    """
    Perform Lagrange interpolation for three points.
        Args: 
            p0, p1, p2: tuples of the form (x, y) where x is the x-coordinate and y is the y-coordinate.
        Returns:
            A sympy expression representing the Lagrange polynomial.
    """

    x = sp.symbols('x')

    poly = 0

    for i, (xi, yi) in enumerate(points):
        term = yi
        for j, (xj, yj) in enumerate(points):
            if i != j:
                term *= (x - xj) / (xi - xj)
        poly += term

    return poly