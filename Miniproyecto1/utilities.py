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