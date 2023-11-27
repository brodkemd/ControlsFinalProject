import control as ctrl
import numpy as np
import matplotlib.pyplot as plt

# Define the system matrices for a 2-input, 2-output system
A = np.array([[0, 1, 0, 0],
              [0, 0, 1, 0],
              [-1, -2, -3, -1],
              [0, 0, 0, -2]])

B = np.array([[0, 0],
              [0, 0],
              [1, 0],
              [0, 1]])

C = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0]])

D = np.array([[0, 0],
              [0, 0]])

# Create the state-space system
sys = ctrl.ss(A, B, C, D)

# Augment the system with an integrator
Ai = np.vstack([
    np.hstack([ A, np.zeros((A.shape[0], C.shape[0]))]),
    np.hstack([-C, np.zeros((C.shape[0], C.shape[0]))])
])

Bi = np.vstack([B, np.zeros((C.shape[0], B.shape[1]))])

Ci = np.hstack([C, np.zeros((C.shape[0], C.shape[0]))])

exit()
sys_augmented = ctrl.ss(Ai, Bi, Ci, D)

# Design a full state feedback controller with an integrator
K = ctrl.place(Ai, Bi, [-2, -2.1, -2.2, -2.3, -2.4, -1])

# Closed-loop system with full state feedback and integrator
sys_cl = ctrl.feedback(sys_augmented, K)

# Time vector
time = np.linspace(0, 10, 1000)

# Step response of the closed-loop system
time, response = ctrl.step_response(sys_cl, time)

# Plot the step response for each output
plt.plot(time, response[0], label='Output 1')
plt.plot(time, response[1], label='Output 2')
plt.xlabel('Time')
plt.ylabel('Response')
plt.title('Step Response with Full State Feedback and Integrator (MIMO)')
plt.legend()
plt.grid(True)
plt.show()
