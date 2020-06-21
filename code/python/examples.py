from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

#fig = plt.figure(figsize=plt.figaspect(1.))
#ax = fig.add_subplot(111, projection='3d')

fig = plt.figure()
ax = plt.axes(projection='3d')


# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)
ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
zdata = 15 * np.random.random(100)
xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');

#plt.show()


def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(X, Y, Z, color='black')
ax.set_title('wireframe');

plt.show()


x = [1, 2, 3]
y = [4, 5, 6]

x1 = np.linspace(0, 10, 1000)

figure, axes = plt.subplots(nrows=2, ncols=2)

axes[0, 0].plot(x, y)

#create specific subplots

axes[0, 1].plot(x1, np.sin(x1))
axes[1, 0].plot(x1, np.cos(x1))
axes[1, 1].plot(range(10))


figure.tight_layout()

plt.show()
