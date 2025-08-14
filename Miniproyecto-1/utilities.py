import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from numericalMethods import interpolation

# Test interpolation function
points = [(0, 1), (1, 2), (2, 0)]
interpolation(points)

def Lagrange_interpolation_3(p0, p1, p2):
    """
    Perform Lagrange interpolation for three points.
        Args: 
            p0, p1, p2: tuples of the form (x, b) where x is the x-coordinate and b is the y-coordinate.
        Returns:
            A sympy expression representing the Lagrange polynomial.
    """
    x0 = p0[0]
    x1 = p1[0]
    x2 = p2[0]
    b0 = p0[1]
    b1 = p1[1]
    b2 = p2[1]

    x = sp.symbols('x')

    phi0 = ((x-x1)*(x-x2))/((x0-x1)*(x0-x2))
    phi1 = ((x-x0)*(x-x2))/((x1-x0)*(x1-x2))
    phi2 = ((x-x0)*(x-x1))/((x2-x0)*(x2-x1))

    return b0*phi0 + b1*phi1 + b2*phi2

