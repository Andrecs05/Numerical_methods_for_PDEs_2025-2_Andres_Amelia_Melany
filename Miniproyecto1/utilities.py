import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


def read_interpolation_points(file_path):
    points = []
    with open(file_path, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split(','))
            points.append((x, y))
    return points

def plot_points(points):
    x_vals, y_vals = zip(*points)
    plt.figure(figsize=(8, 10))
    plt.scatter(x_vals, y_vals, color='blue', marker='o')  # scatter instead of plot
    plt.title('Puntos chocolatera')
    plt.xlabel('X (Cm)')
    plt.ylabel('Y (Cm)')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

def plot_function_expresion(func_dict, same_plot=True):
    """
    Plots functions from a dictionary of the form:
    {
        'name': [sympy_expression, lower_bound, upper_bound],
        ...
    }
    same_plot: If True, all functions are on the same axes;
               If False, each function has its own subplot.
    """
    x = sp.symbols('x')

    if same_plot:
        plt.figure(figsize=(10, 6))
        for name, (func, lower_bound, upper_bound) in func_dict.items():
            f = sp.lambdify(x, func, modules=['numpy'])
            x_vals = np.linspace(lower_bound, upper_bound, 300)
            y_vals = f(x_vals)
            plt.plot(x_vals, y_vals, label=name)

        plt.title('Function Plot')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.legend()
        plt.show()

    else:
        num_funcs = len(func_dict)
        fig, axes = plt.subplots(1, num_funcs, figsize=(5 * num_funcs, 5))

        if num_funcs == 1:
            axes = [axes]  # Make it iterable

        for ax, (name, (func, lower_bound, upper_bound)) in zip(axes, func_dict.items()):
            f = sp.lambdify(x, func, modules=['numpy'])
            x_vals = np.linspace(lower_bound, upper_bound, 300)
            y_vals = f(x_vals)
            ax.plot(x_vals, y_vals, label=name)
            ax.set_title(name)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.grid(True)
            ax.legend()

        plt.tight_layout()
        plt.show()

def plot_function(func_dict, same_plot=True):
    if same_plot:
        plt.figure(figsize=(10, 6))
        for name, (t, u) in func_dict.items():
            plt.plot(t, u, label=name)

        plt.xlabel("t")
        plt.ylabel("u(t)")
        plt.legend()
        plt.title("Exponential Decay: Numerical Schemes")
        plt.show()
    else:
        num_funcs = len(func_dict)
        fig, axes = plt.subplots(1, num_funcs, figsize=(6 * num_funcs, 5), squeeze=False)
        axes = axes[0]  # Get the 1D array of axes

        for ax, (name, (t, u)) in zip(axes, func_dict.items()):
            ax.plot(t, u, label=name)
            ax.set_xlabel("t")
            ax.set_ylabel("u(t)")
            ax.set_title(f"Exponential Decay: {name}")
            ax.legend()

        plt.tight_layout()
        plt.show()

def plot_exponential_decay_stability(approx_list, exact_solution, T, method):

    n = len(approx_list)

    exact_steps = len(exact_solution)
    x_values_exact = np.linspace(0, T, exact_steps)

    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)
    axes = axes[0]

    for i, solution in enumerate(approx_list):
        ax = axes[i]

        steps = len(solution)
        x_values = np.linspace(0, T, steps)

        ax.plot(x_values, solution, label='Approximate Solution ' f'dt = {T / (steps - 1):.2f}')
        ax.plot(x_values_exact, exact_solution, label='Exact Solution', linestyle='--')
        ax.set_xlabel('t')
        ax.set_ylabel('u(t)')
        ax.set_title(method)
        ax.legend()
        ax.grid(True)

    plt.tight_layout()
    plt.show()

def plot_log_dt_log_E(dt_values, E_values, methods):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    log_dt = np.log(dt_values)
    log_E = {method: np.log(E_values[method]) for method in methods}

    # Mask for linear zone: log(dt) < -2
    mask = log_dt < -1.5

    slopes = {}
    for ax, method in zip(axes, methods):
        x = log_dt[mask]
        y = log_E[method][mask]
        # Linear fit for the linear zone
        slope = np.polyfit(x, y, 1)[0]
        slopes[method] = slope

        ax.plot(log_dt, log_E[method], label=f'{method} (slope = {slope:.2f})')
        ax.set_title(f'Log-Log Plot: {method}')
        ax.set_xlabel('Log(dt)')
        ax.set_ylabel('Log(Error)')
        ax.grid(True)
        ax.legend()

        # Highlight the linear zone
        ax.axvspan(log_dt[mask][0], log_dt[mask][-1], color='yellow', alpha=0.2)

    plt.tight_layout()
    plt.show()
