import numpy as np
import matplotlib.pyplot as plt

# data = np.loadtxt("data_files/task_5c_vv_dt=0.010000.txt", unpack=True)
# data = np.loadtxt("data_files/task_5g.txt", unpack=True)
data = np.load("data_files/task_5g.npy")
x = data[1]
y = data[2]
# x = x[10*5000:]
# y = y[10*5000:]

eps = 0.0001
idx = np.where(abs(y)<eps)[0]
idxidx = np.where(np.diff(idx)>10)[0]
idx2 = idx[idxidx]

x_one_round = x[0:idx2[2]]
y_one_round = y[0:idx2[2]]
x_final_round = x[idx2[-2]:idx[-1]]
y_final_round = y[idx2[-2]:idx[-1]]

r1 = np.sqrt(x_one_round**2 + y_one_round**2)
r2 = np.sqrt(x_final_round**2 + y_final_round**2)

p_idx1 = np.argmin(r1)
p_idx2 = np.argmin(r2) + idx2[-2]

print(f"x1 = {x[p_idx1]}")
print(f"y1 = {y[p_idx1]}")
print(f"x2 = {x[p_idx2]}")
print(f"y2 = {y[p_idx2]}")
print()
print(y[p_idx1]/x[p_idx1])
print(y[p_idx2]/x[p_idx2])
theta1 = np.tan(y[p_idx1]/x[p_idx1])
theta2 = np.tan(y[p_idx2]/x[p_idx2])
print()
print(f"theta1 = {theta1}")
print(f"theta2 = {theta2}")
precession = (theta2 - theta1)
print()
print(f"precession (arcsec)  = {precession*180/np.pi*60*60}")
print(f"precession (radians) = {precession}")
print(f"precession (degrees) = {precession*180/np.pi}")
print()
print(f"43 buesec (arcsec)  = {43}")
print(f"43 buesec (radians) = {43/3600*np.pi/180}")
print(f"43 buesec (degrees) = {43/3600}")

# plt.plot(np.array([0,1]), np.array([0,0]))
# plt.plot(np.array([0,1*np.cos(43/3600*np.pi/180)]), np.array([0,1*np.sin(43/3600*np.pi/180)]))
# plt.plot(np.cos(np.linspace(0,2*np.pi,1000)), np.sin(np.linspace(0,2*np.pi,1000)))
# plt.show()

# print(idx2[0])
# print(idx2[2])
# print(idx2[-2])
# print(idx[-1])
# print(idx2)

# print(p_idx1)
# print(p_idx2)

# plt.plot(x[10*5000:], y[10*5000:])
plt.plot(x[::1000], y[::1000])
plt.plot(x[p_idx1], y[p_idx1], "b*")
plt.plot(x[p_idx2], y[p_idx2], "r*")
plt.plot(0, 0, "y*")
plt.show()