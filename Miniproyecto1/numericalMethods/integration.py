import sympy as sp

def midpoint_rule(expr, lower_bound, upper_bound, n):
    """
    Approximates the definite integral of a given expression using the midpoint rule.

    The midpoint rule is a numerical integration method that divides the interval [lower_bound, upper_bound]
    into 'n' subintervals of equal width. For each subinterval, the function is evaluated at the midpoint,
    and the sum of these values multiplied by the width of the subintervals gives the approximate integral.

    Formula:
        ∫[a, b] f(x) dx ≈ Δx * Σ_{i=0}^{n-1} f(a + (i + 0.5)Δx)
        where Δx = (b - a) / n

    Parameters:
        expr (sympy expression): The mathematical expression to integrate, as a sympy symbolic expression in 'x'.
        lower_bound (float): The lower limit of integration.
        upper_bound (float): The upper limit of integration.
        n (int): The number of subintervals to use in the approximation.

    Returns:
        float: The approximate value of the definite integral over [lower_bound, upper_bound].
    """
    x = sp.symbols('x')
    dx = (upper_bound - lower_bound) / n
    total = 0
    for i in range(n):
        mid_point = lower_bound + (i + 0.5) * dx
        total += expr.subs(x, mid_point) * dx
    return total

def trapezoidal_rule(expr, lower_bound, upper_bound, n):
    """
    Approximates the definite integral of a given symbolic expression using the trapezoidal rule.

    The trapezoidal rule estimates the area under a curve by dividing the interval [lower_bound, upper_bound]
    into 'n' subintervals, then approximating the area under the curve as a sum of trapezoids.

    General formula (summation notation):
        ∫[a, b] f(x) dx ≈ (Δx / 2) * [f(x₀) + 2 * Σ_{i=1}^{n-1} f(xᵢ) + f(xₙ)]
        where:
            x₀ = a
            xₙ = b
            xᵢ = a + i * Δx, for i = 0, 1, ..., n
            Δx = (b - a) / n

    Args:
        expr (sympy.Expr): The symbolic expression to integrate, as a function of x.
        lower_bound (float): The lower limit of integration (a).
        upper_bound (float): The upper limit of integration (b).
        n (int): The number of subintervals (trapezoids) to use.

    Returns:
        float: The approximate value of the definite integral over [lower_bound, upper_bound].
    """
    x = sp.symbols('x')
    dx = (upper_bound - lower_bound) / n
    total = 0
    for i in range(n):
        total += (expr.subs(x, lower_bound + i * dx) + expr.subs(x, lower_bound + (i + 1) * dx)) * dx / 2
    return total

def simpson_rule(expr, lower_bound, upper_bound, n):
    """
    Approximates the definite integral of a given symbolic expression using Simpson's Rule.

    Simpson's Rule formula:
        ∫[a, b] f(x) dx ≈ (dx/3) * [f(a) + 4f(a+dx) + 2f(a+2dx) + ... + 4f(b-dx) + f(b)]
    where dx = (b - a) / n, and n is an even integer.

    Parameters:
        expr (sympy.Expr): The symbolic expression to integrate, as a function of x.
        lower_bound (float): The lower limit of integration (a).
        upper_bound (float): The upper limit of integration (b).
        n (int): The number of subintervals (must be even; will be incremented if odd).

    Returns:
        float: The approximate value of the definite integral over [lower_bound, upper_bound].

    Notes:
        - If n is odd, it will be incremented by 1 to ensure an even number of intervals.
        - The function evaluates the expression at equally spaced points between the bounds.
    """
    x = sp.symbols('x')
    if n % 2 == 1:
        n += 1  # Simpson's rule requires an even number of intervals
    dx = (upper_bound - lower_bound) / n
    total = expr.subs(x, lower_bound) + expr.subs(x, upper_bound)
    for i in range(1, n, 2):
        total += 4 * expr.subs(x, lower_bound + i * dx)
    for i in range(2, n-1, 2):
        total += 2 * expr.subs(x, lower_bound + i * dx)
    return total * dx / 3

def simpson_38_rule(expr, lower_bound, upper_bound, n):
    """
    Approximates the definite integral of a given symbolic expression using Simpson's 3/8 Rule.

    Simpson's 3/8 Rule formula:
        ∫[a, b] f(x) dx ≈ (3h/8) * [f(x₀) + 3f(x₁) + 3f(x₂) + 2f(x₃) + 3f(x₄) + ... + 2f(x_{n-3}) + 3f(x_{n-2}) + 3f(x_{n-1}) + f(xₙ)]
    where h = (b - a) / n, and n is a multiple of 3.

    Parameters:
        expr (sympy.Expr): The symbolic expression to integrate, as a function of x.
        lower_bound (float): The lower limit of integration (a).
        upper_bound (float): The upper limit of integration (b).
        n (int): The number of subintervals (must be a multiple of 3; will be incremented if not).

    Returns:
        float: The approximate value of the definite integral over [lower_bound, upper_bound].
    """
    x = sp.symbols('x')
    if n % 3 != 0:
        n += 3 - (n % 3)  # Ensure n is a multiple of 3
    h = (upper_bound - lower_bound) / n

    total = expr.subs(x, lower_bound) + expr.subs(x, upper_bound)

    for i in range(1, n):
        xi = lower_bound + i * h
        if i % 3 == 0:
            total += 2 * expr.subs(x, xi)
        else:
            total += 3 * expr.subs(x, xi)

    return total * 3 * h / 8

def gauss_legendre(expr, lower_bound, upper_bound, n):
    """
    Approximates the definite integral of a given symbolic expression over the interval [-1, 1]
    using Gauss-Legendre quadrature.

    Args:
        expr (sympy.Expr): The symbolic expression to integrate, as a function of x.
        n (int): The number of quadrature points to use (must be 1, 2, or 3).

    Returns:
        float: The approximate value of the definite integral over [-1, 1].
    """
    x = sp.symbols('x')

    nodes = {
        1: [0], 
        2: [-sp.sqrt(3)/3, sp.sqrt(3)/3], 
        3: [-sp.sqrt(3/5), 0, sp.sqrt(3/5)], 
        4: [-sp.sqrt(3/7 - 2/7 * sp.sqrt(6/5)), sp.sqrt(3/7 - 2/7 * sp.sqrt(6/5)), -sp.sqrt(3/7 + 2/7 * sp.sqrt(6/5)), sp.sqrt(3/7 + 2/7 * sp.sqrt(6/5))], 
        5: [-1/3 * sp.sqrt(5 - 2 * sp.sqrt(10/7)), -1/3 * sp.sqrt(5 + 2 * sp.sqrt(10/7)), 0, 1/3 * sp.sqrt(5 - 2 * sp.sqrt(10/7)), 1/3 * sp.sqrt(5 + 2 * sp.sqrt(10/7))]
    }

    weights = {
        1: [2],
        2: [1, 1],
        3: [5/9, 8/9, 5/9],
        4: [(18 + sp.sqrt(30))/36, (18 + sp.sqrt(30))/36, (18 - sp.sqrt(30))/36, (18 - sp.sqrt(30))/36],
        5: [(322 + 13*sp.sqrt(70))/900, (322 + 13*sp.sqrt(70))/900, 128/225, (322 - 13*sp.sqrt(70))/900, (322 - 13*sp.sqrt(70))/900]
    }

    if n not in nodes or n not in weights:
        raise ValueError("Invalid number of quadrature points. Must be 1, 2, 3, 4, or 5.")
    
    total = sum(weights[n][i] * expr.subs(x, (upper_bound - lower_bound) / 2 * nodes[n][i] + (upper_bound + lower_bound) / 2) for i in range(n))

    return total * (upper_bound - lower_bound) / 2


import sympy as sp

def test(expr, a, b, n):
    x = sp.symbols('x')

    # Exact symbolic integration, then convert to float
    result_exact = sp.integrate(expr, (x, a, b)).evalf()
    print(f"Test Function - Exact Result: {result_exact}")

    # Middle Point Rule
    result_middle = midpoint_rule(expr, a, b, n)
    print(f"Test Function - Middle Point Result: {sp.N(result_middle)}")

    # Trapezoidal Rule
    result_trapezoidal = trapezoidal_rule(expr, a, b, n)
    print(f"Test Function - Trapezoidal Result: {sp.N(result_trapezoidal)}")

    # Simpson's Rule
    result_simpson = simpson_rule(expr, a, b, n)
    print(f"Test Function - Simpson's Result: {sp.N(result_simpson)}")

    # Simpson's 3/8 Rule
    result_simpson_38 = simpson_38_rule(expr, a, b, n)
    print(f"Test Function - Simpson's 3/8 Result: {sp.N(result_simpson_38)}")

    # Gauss-Legendre Quadrature
    result_gauss_legendre = gauss_legendre(expr, a, b, 3)
    print(f"Test Function - Gauss-Legendre Result: {sp.N(result_gauss_legendre)}")



# Run test
if __name__ == "__main__":
    x = sp.symbols('x')

    # ∫₀^π sin(x) dx
    print("Function:", sp.sin(x))
    test(sp.sin(x), 0, sp.pi, 100)
    # Exact:        2.00000000000000
    # Midpoint:     2.00008224907099
    # Trapezoidal:  1.99983550388744
    # Simpson 1/3:  2.00000001082450
    # Simpson 3/8:  2.00000002250282
    # Gauss-Leg.:   2.00138891360774

    # ∫₀^π sin²(x) dx
    print("Function:", sp.sin(x)**2)
    test(sp.sin(x)**2, 0, sp.pi, 100)
    # Exact:        1.57079632679490
    # Midpoint:     1.57079632679490
    # Trapezoidal:  1.57079632679490
    # Simpson 1/3:  1.57079632679490
    # Simpson 3/8:  1.57079632679490
    # Gauss-Leg.:   1.60606730241802

    # ∫₋₁^1 (x³ + 5) dx
    print("Function:", x**3 + 5)
    test(x**3 + 5, -1, 1, 100)
    # Exact:        10.0000000000000
    # Midpoint:     10.0000000000000
    # Trapezoidal:  10.0000000000000
    # Simpson 1/3:  10.0000000000000
    # Simpson 3/8:  10.0000000000000
    # Gauss-Leg.:   10.0000000000000

    # ∫₁₀²⁰ x³ dx
    print("Function:", x**3)
    test(x**3, 10, 20, 100)
    # Exact:        37500.0000000000
    # Midpoint:     37499.6250000000
    # Trapezoidal:  37500.7500000000
    # Simpson 1/3:  37500.0000000000
    # Simpson 3/8:  37500.0000000000
    # Gauss-Leg.:   37500.0000000000


'''

point_sets = (len(chocolatera_points) - 1) // 2

Volume = 0
for i in range(point_sets - 1):
    k = i + 1
    expr = interpolation(chocolatera_points[2*k-1], chocolatera_points[2*k], chocolatera_points[2*k+1])
    dV = sp.pi * expr**2
    Volume += sp.integrate(dV, (x, chocolatera_points[2*k-1][0], chocolatera_points[2*k+1][0]))

print(f"Volume of Chocolatera: {Volume.evalf()}")
 '''

