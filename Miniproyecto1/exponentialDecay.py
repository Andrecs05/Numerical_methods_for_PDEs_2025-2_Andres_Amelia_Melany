import numpy as np

def forward_euler(u0, a, h, steps):
    u = np.zeros(steps+1)
    u[0] = u0
    for n in range(steps):
        u[n+1] = (1 - a*h) * u[n]
    return u

def backward_euler(u0, a, h, steps):
    u = np.zeros(steps+1)
    u[0] = u0
    for n in range(steps):
        u[n+1] = u[n] / (1 + a*h)
    return u

def crank_nicolson(u0, a, h, steps):
    u = np.zeros(steps+1)
    u[0] = u0
    for n in range(steps):
        u[n+1] = ((1 - 0.5*a*h) / (1 + 0.5*a*h)) * u[n]
    return u

def theta_rule(u0, a, h, steps, theta):
    u = np.zeros(steps+1)
    u[0] = u0
    for n in range(steps):
        u[n+1] = ((1 - (1-theta)*a*h) / (1 + theta*a*h)) * u[n]
    return u

def exact_solution(u0, a, h, steps):
    t = np.linspace(0, steps*h, steps+1)
    return u0 * np.exp(-a*t)