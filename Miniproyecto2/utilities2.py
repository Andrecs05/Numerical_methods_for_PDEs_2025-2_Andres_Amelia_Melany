import numpy as np
import matplotlib.pyplot as plt

def boundary_conditions(p, u, A, b):
    '''
    Apply Dirichlet boundary conditions at node p with potential u.
    Modifies matrix A and vector b in place.

    Parameters:
    p : int
        Index of the node where the boundary condition is applied.
    u : float
        Potential value at the boundary.
    A : ndarray
        Coefficient matrix.
    b : ndarray
        Right-hand side vector.
    Returns:
    A : ndarray
        Modified coefficient matrix.
    b : ndarray
        Modified right-hand side vector.
    '''
    A[p, :] = 0
    A[p, p] = 1
    b[p] = u

    return A, b

def plot_solution(V , n_x_nodes, n_y_nodes, X, Y):
    '''
    Plot the 2D potential distribution.
    '''
    V_matrix = V.reshape((n_y_nodes, n_x_nodes))

    max_abs = np.max(np.abs(V_matrix))
    plt.imshow(V_matrix, extent=[0, X, 0, Y], origin='lower', aspect='auto', cmap='seismic', vmin=-max_abs, vmax=max_abs)
    plt.colorbar(label='Potential (V)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Potential Distribution')
    plt.show()