import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# SYSTEM A
def system_a(t, y):
    y1, y2 = y
    dy1 = 2*y1 + 2*y2
    dy2 = 5*y1 - y2
    return [dy1, dy2]

# SYSTEM B
def system_b(t, y):
    y1, y2 = y
    dy1 = 3*y1 + 2*y2
    dy2 = 2*y1 + 3*y2
    return [dy1, dy2]

# ploting phase portrait
def plot_phase_portrait(system, title, y0=None, t_span=20, grid_range=50, grid_size=20):
    y1 = np.linspace(-grid_range, grid_range, grid_size)
    y2 = np.linspace(-grid_range, grid_range, grid_size)
    Y1, Y2 = np.meshgrid(y1, y2)
    U, V = np.zeros(Y1.shape), np.zeros(Y2.shape)

    for i in range(Y1.shape[0]):
        for j in range(Y1.shape[1]):
            dydt = system(0, [Y1[i, j], Y2[i, j]])
            U[i, j] = dydt[0]
            V[i, j] = dydt[1]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.streamplot(Y1, Y2, U, V, color='black')
    ax.set_title(title)
    ax.set_xlabel('$y_1$')
    ax.set_ylabel('$y_2$')
    ax.set_xlim(-grid_range, grid_range)
    ax.set_ylim(-grid_range, grid_range)
    ax.axhline(0, color='black', lw=0.5)
    ax.axvline(0, color='black', lw=0.5)

    if y0:
        # Forward integration
        sol_forward = solve_ivp(system, [0, t_span], y0, dense_output=True)
        y_fwd = sol_forward.sol(np.linspace(0, t_span, 300))
        ax.plot(y_fwd[0], y_fwd[1], 'r-')

        # Backward integration
        sol_backward = solve_ivp(system, [0, -t_span], y0, dense_output=True)
        y_bwd = sol_backward.sol(np.linspace(0, -t_span, 300))
        ax.plot(y_bwd[0], y_bwd[1], 'r-')

    ax.legend()
    plt.grid()
    plt.show()

plot_phase_portrait(system_a, title='(a) Phase Portrait', y0=(0, 7))
plot_phase_portrait(system_b, title='(b) Phase Portrait', y0=(0.5, -0.5))