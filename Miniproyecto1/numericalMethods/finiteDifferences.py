import matplotlib.pyplot as plt

# Forward difference: f'(x) ≈ (f(x + h) - f(x)) / h
def forward_difference(f, var, x, h):
    return (f.subs(var, x + h) - f.subs(var, x)) / h

# Backward difference: f'(x) ≈ (f(x) - f(x - h)) / h
def backward_difference(f, var, x, h):
    return (f.subs(var, x) - f.subs(var, x - h)) / h

# Central difference: f'(x) ≈ (f(x + h) - f(x - h)) / (2h)
def central_difference(f, var, x, h):
    return (f.subs(var, x + h) - f.subs(var, x - h)) / (2 * h)

# Weighted finite difference: f'(x) ≈ (1 - θ) * backward + θ * forward
# = (1 - θ) * (f(x) - f(x - h))/h + θ * (f(x + h) - f(x))/h
def weighted_finite_difference(f, var, x, h, theta):
    back = (f.subs(var, x) - f.subs(var, x - h)) / h
    fwd  = (f.subs(var, x + h) - f.subs(var, x)) / h
    return (1 - theta) * back + theta * fwd

def explicit_euler(f, t0, x0, h, steps):
    """Forward (explicit) Euler method. Returns list of (t, x) pairs."""
    t, x = t0, x0
    result = [(t, x)]
    for _ in range(steps):
        x = x + h * f(t, x)
        t = t + h
        result.append((t, x))
    return result


def implicit_euler(f, t0, x0, h, steps, iters=20, tol=1e-12):
    """Backward (implicit) Euler method using fixed-point iteration. Returns list of (t, x) pairs."""
    t, x = t0, x0
    result = [(t, x)]
    for _ in range(steps):
        t_next = t + h
        y = x  # initial guess
        for _ in range(iters):
            y_new = x + h * f(t_next, y)   # implicit: f at (t+h, y)
            if abs(y_new - y) < tol:
                break
            y = y_new
        x, t = y_new, t_next
        result.append((t, x))
    return result


def theta_method(f, t0, x0, h, steps, theta=0.5):
    """
    General θ-method:
      theta = 0   -> forward Euler
      theta = 1   -> backward Euler
      theta = 0.5 -> Crank-Nicolson (trapezoidal rule)
    Returns a list of (t, x) pairs.
    """
    t, x = t0, x0
    result = [(t, x)]
    for _ in range(steps):
        t_next = t + h
        y = x
        for _ in range(20):
            y_new = x + h * ((1 - theta) * f(t, x) + theta * f(t_next, y))
            if abs(y_new - y) < 1e-12:
                break
            y = y_new
        x_new = y_new
    x, t = x_new, t_next
    result.append((t, x))
    return result
