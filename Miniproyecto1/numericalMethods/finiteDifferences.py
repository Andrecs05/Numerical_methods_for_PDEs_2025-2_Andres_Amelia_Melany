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

# Euler's method for ODEs: x_{n+1} = x_n + h * f(x_n)
def euler_method(f, var, x0, t0, h, steps):
    x = x0
    t = t0
    for _ in range(steps):
        x = x + h * f.subs(var, x)
        t = t + h
    return x
