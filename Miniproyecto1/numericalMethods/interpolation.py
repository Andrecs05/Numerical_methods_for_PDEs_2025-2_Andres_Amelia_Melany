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

def 

# Plot
x = sp.symbols('x')

point_sets = (len(chocolatera_points) - 1) // 2

Volume = 0
for i in range(point_sets - 1):
    k = i + 1
    expr = interpolation(chocolatera_points[2*k-1], chocolatera_points[2*k], chocolatera_points[2*k+1])
    dV = sp.pi * expr**2
    Volume += sp.integrate(dV, (x, chocolatera_points[2*k-1][0], chocolatera_points[2*k+1][0]))

print(f"Volume of Chocolatera: {Volume.evalf()}")

