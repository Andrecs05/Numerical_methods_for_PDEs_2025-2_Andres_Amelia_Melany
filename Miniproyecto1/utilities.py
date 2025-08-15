import matplotlib.pyplot as plt

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
    plt.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue')
    plt.title('Swapped and Sorted Points Plot')
    plt.xlabel('Original Y Coordinate (Now X)')
    plt.ylabel('Original X Coordinate (Now Y)')
    plt.grid(True)
    plt.axis('equal')

d