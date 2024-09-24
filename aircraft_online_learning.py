import numpy as np
from scipy.signal import cont2discrete
import cvxpy as cp
import matplotlib.pyplot as plt

# Define global variables
global A1, A2, B1, B2, m, n, T

# Aircraft engine
Ts = 0.1

# Define the continuous-time matrices
A1c = np.array([[-0.2296, 0.9931], [0.02436, -0.2046]])
B1c = np.array([[-0.0434, -0.01145], [-1.73, -0.517]])

# Discretize the continuous-time systems
A1, B1, _, _, _ = cont2discrete((A1c, B1c, np.eye(2), np.zeros((2,2))), Ts)

# Define second system
A2c = np.array([[-1.175, 0.9871], [-8.458, -0.8776]])
B2c = np.array([[-0.194, -0.03593], [-19.29, -3.803]])

# Discretize second system
A2, B2, _, _, _ = cont2discrete((A2c, B2c, np.eye(2), np.zeros((2,2))), Ts)

# Store systems
A_sw = [A1, A2]
B_sw = [B1, B2]

# Data setup
np.random.seed(0)
m = B1.shape[1]
n = A1.shape[0]
N = (m + 1) * n + m
T = 2 * N - 1

magx = 0.5
x = -magx + (magx + magx) * np.random.rand(n, 1)
U = -magx + (magx + magx) * np.random.rand(m, T)
X = np.copy(x)

# Run the system for T steps
for i in range(T):
    x = A1 @ x + B1 @ U[:, i:i+1]
    X = np.hstack([X, x])

# Go online simulation
time = 150
time_faults = [10, 25, 50, 67, 95, 110]
sw_seq = [2, 1, 2, 1, 2, 1]
As = A1
Bs = B1
eps = 0.001
sw_con = 1
count = 1
Xcl = []
Ucl = []
x0 = np.copy(x)
K = U[:, -1].reshape(-1, 1)

for t in range(time):
    print(t)
    
    # Fault handling
    if t == time_faults[sw_con - 1]:
        As = A_sw[sw_seq[sw_con - 1] - 1]
        Bs = B_sw[sw_seq[sw_con - 1] - 1]
        if sw_con < len(time_faults):
            sw_con += 1

    X0 = X[:, -T:-1]
    U0 = U[:, -T+1:]
    X1 = X[:, -T+1:]

    old_K = np.copy(K)

    # Adjust Q dimensions to match U0's shape for the matrix multiplication
    Q = cp.Variable((T - 1, n))  # Ensure Q aligns with U0 in dimension

    P = cp.Variable((n, n), symmetric=True)
    L = cp.Variable((m, m), symmetric=True)
    gam = cp.Variable()

    constraints = [
        cp.bmat([[L, U0 @ Q], [Q.T @ U0.T, P]]) >> 0,
        cp.bmat([[np.eye(n) - P, X1 @ Q], [Q.T @ X1.T, -P]]) << 0,
        X0 @ Q == P,
        cp.trace(P) + cp.trace(L) <= gam
    ]

    prob = cp.Problem(cp.Minimize(gam), constraints)
    prob.solve(solver=cp.SCS, verbose=False)

    # Ensure P.value is not None and invertible
    if P.value is not None and np.linalg.det(P.value) != 0:
        # Solve K by solving the linear system instead of direct division
        K = np.linalg.solve(P.value, U0 @ Q.value)
    else:
        K = old_K  # If P is None or non-invertible, revert to old K

    # Control input u(t)
    e = -eps + (eps + eps) * np.random.randn(m, 1)
    u = K @ x + np.linalg.norm(x) * e

    # Store x(t), u(t), norm(x)
    Xcl.append(x.flatten())
    Ucl.append(u.flatten())

    # Update x for next step
    x = As @ x + Bs @ u
    U = np.hstack([U, u])
    X = np.hstack([X, x])

Xcl = np.array(Xcl).T
Ucl = np.array(Ucl).T

# Plotting
colors = ['#ff9100', '#1e96fc']  # colors for state variables

# State plots
plt.figure()
plt.subplot(211)
for j in range(n):
    plt.plot(np.arange(0, Ts * Xcl.shape[1], Ts), Xcl[j, :], '.-', color=colors[j], markersize=7, linewidth=1.2)
    plt.ylabel('x')
    plt.grid(True)

# Mark fault times
for index in time_faults:
    for j in range(n):  # Loop over each state dimension (e.g., x1, x2)
        plt.plot(index * Ts, Xcl[j, index], '*', markersize=11, linewidth=1.1, color="#db2b39")

plt.legend(['x1', 'x2', 'switches'])
plt.show()

# Control input plot
plt.figure()
plt.subplot(211)
coloru = ['#487bea', '#da627d']  # colors for control inputs
for j in range(m):
    plt.plot(np.arange(0, Ts * Ucl.shape[1], Ts), Ucl[j, :], '.-', color=coloru[j], markersize=7, linewidth=1.2)
    plt.ylabel('u')
    plt.grid(True)

plt.xlabel('time[s]')
plt.legend(['u1', 'u2'])
plt.show()
