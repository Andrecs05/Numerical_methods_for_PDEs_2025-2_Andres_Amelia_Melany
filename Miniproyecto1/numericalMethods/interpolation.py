import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

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
        for j, (xj, _) in enumerate(points):
            if i != j:
                term *= (x - xj) / (xi - xj)
        poly += term

    return sp.simplify(poly)

def piecewise_lagrange_interpolation(points, pts_per_interval):
    """
    Perform piecewise Lagrange interpolation on subintervals defined by the given points.

    Args:
        points: A list of tuples (x, y) representing the data points.
        pts_per_interval: The number of points in each subinterval to use for interpolation.

    Returns:
        A sympy expression representing the piecewise Lagrange polynomial.
    """
    x = sp.symbols('x')
    pieces = []
    num_points = len(points)
    start = 0

    while start < num_points - 1:
        end = min(start + pts_per_interval, num_points)
        sub_points = points[start:end]
        a, b = sub_points[0][0], sub_points[-1][0]
        poly = lagrange_interpolation(sub_points)
        pieces.append((poly, (x >= a) & (x <= b)))
        start += pts_per_interval - 1  # overlap by 1 point

    return sp.Piecewise(*pieces)

def graph_interpolation(expr, points, title):
    x = sp.symbols('x')
    f = sp.lambdify(x, expr, modules=['numpy'])
    
    # Original points
    x_vals_points = [p[0] for p in points]
    y_vals_points = [p[1] for p in points]
    
    # Dense sampling for smooth polynomial curve
    x_min, x_max = min(x_vals_points), max(x_vals_points)
    x_dense = np.linspace(x_min, x_max, 400)
    y_dense = f(x_dense)

    plt.figure(figsize=(8, 5))
    plt.plot(x_vals_points, y_vals_points, 'ro', label='Data Points')
    plt.plot(x_dense, y_dense, 'b-', label='Interpolation Polynomial')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.grid(True)
    plt.show()

def test_lagrange_interpolation():
    # Define 20 test points tangent
    test_function = lambda x: np.tan(x)
    points = [(i, test_function(i)) for i in np.linspace(0, 2 * np.pi, 20)]
    x = sp.symbols('x')

    # Exact polynomial (Lagrange interpolation)
    exact_poly = lagrange_interpolation(points)
    print(f"Exact Lagrange Polynomial: {exact_poly}")

    # Test piecewise interpolation
    n = 3
    piecewise_poly = piecewise_lagrange_interpolation(points, n)
    print(f"\nPiecewise Lagrange Polynomial ({n} points per interval): {piecewise_poly}")

    # Graph the interpolation
    graph_interpolation(exact_poly, points, "Exact Lagrange Interpolation")

    graph_interpolation(piecewise_poly, points, "Piecewise Lagrange Interpolation")

if __name__ == "__main__":
    test_lagrange_interpolation()