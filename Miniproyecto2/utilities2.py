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

import numpy as np

def print_matrix(A, max_rows=10, max_cols=10):
    """
    Pretty-print the coefficient matrix A.
    Shows only the top-left corner if the matrix is large.
    
    Parameters:
        A (np.ndarray): The matrix to print.
        max_rows (int): Maximum number of rows to display.
        max_cols (int): Maximum number of columns to display.
    """
    rows, cols = A.shape
    print(f"Matrix shape: {rows} x {cols}\n")
    
    # Determine how much to show
    r = min(rows, max_rows)
    c = min(cols, max_cols)
    
    # Slice the top-left block
    submatrix = A[:r, :c]
    
    # Print with formatting
    for i in range(r):
        row_str = " ".join(f"{val:6.1f}" for val in submatrix[i])
        if c < cols:
            row_str += "  ..."
        print(row_str)
    
    if r < rows:
        print("...")
