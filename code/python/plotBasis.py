import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

#--------------------------#
# Own Values               #
#--------------------------#

# read in own coordinates
import csv

ndimGrid = 51
X = np.zeros((ndimGrid, ndimGrid))
Y = np.zeros((ndimGrid, ndimGrid))
Z = np.zeros((ndimGrid, ndimGrid))

X0 = np.zeros((ndimGrid, ndimGrid))
Y0 = np.zeros((ndimGrid, ndimGrid))
Z0 = np.zeros((ndimGrid, ndimGrid))

X1 = np.zeros((ndimGrid, ndimGrid))
Y1 = np.zeros((ndimGrid, ndimGrid))
Z1 = np.zeros((ndimGrid, ndimGrid))

X2 = np.zeros((ndimGrid, ndimGrid))
Y2 = np.zeros((ndimGrid, ndimGrid))
Z2 = np.zeros((ndimGrid, ndimGrid))

Harm0 = np.zeros((ndimGrid, ndimGrid))
Harm1 = np.zeros((ndimGrid, ndimGrid))
Harm2 = np.zeros((ndimGrid, ndimGrid))
Harm3 = np.zeros((ndimGrid, ndimGrid))
Harm4 = np.zeros((ndimGrid, ndimGrid))
Harm5 = np.zeros((ndimGrid, ndimGrid))
Harm6 = np.zeros((ndimGrid, ndimGrid))
Harm7 = np.zeros((ndimGrid, ndimGrid))
Harm8 = np.zeros((ndimGrid, ndimGrid))

idx_row = 0
idx_col = 0

#fill grid_data with values
with open('harmonicTest.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            #fill grid data: row[0]
            idx_col = (line_count-1)%ndimGrid
            idx_row = (line_count-1)//ndimGrid

            print ( str(line_count) + ' ' + str(idx_row) + ' ' + str(idx_col) )

            my = float( row[0])
            phi = float(row[1])

            X[idx_row,idx_col] = np.sqrt(1-my*my)*np.cos(phi)
            Y[idx_row, idx_col] = np.sqrt(1 - my * my) * np.sin(phi)
            Z[idx_row, idx_col] = my

            Harm0[idx_row, idx_col] = row[2]
            Harm1[idx_row, idx_col] = row[3]
            Harm2[idx_row, idx_col] = row[4]
            Harm3[idx_row, idx_col] = row[5]
            Harm4[idx_row, idx_col] = row[6]
            Harm5[idx_row, idx_col] = row[7]
            Harm6[idx_row, idx_col] = row[8]
            Harm7[idx_row, idx_col] = row[9]
            Harm8[idx_row, idx_col] = row[10]
            line_count += 1

            # radius version
            X0[idx_row, idx_col] = X[idx_row, idx_col] * Harm0[idx_row, idx_col]
            Y0[idx_row, idx_col] = Y[idx_row, idx_col] * Harm0[idx_row, idx_col]
            Z0[idx_row, idx_col] = Z[idx_row, idx_col] * Harm0[idx_row, idx_col]

            X1[idx_row, idx_col] = X[idx_row, idx_col] * Harm1[idx_row, idx_col]
            Y1[idx_row, idx_col] = Y[idx_row, idx_col] * Harm1[idx_row, idx_col]
            Z1[idx_row, idx_col] = Z[idx_row, idx_col] * Harm1[idx_row, idx_col]

            X2[idx_row, idx_col] = X[idx_row, idx_col] * Harm2[idx_row, idx_col]
            Y2[idx_row, idx_col] = Y[idx_row, idx_col] * Harm2[idx_row, idx_col]
            Z2[idx_row, idx_col] = Z[idx_row, idx_col] * Harm2[idx_row, idx_col]

    print(f'Processed {line_count} lines.')

# setup colors
fmax, fmin = Harm0.max(), Harm0.min()
fcolors0 = (Harm0) # constan function - fmin)/(fmax - fmin)

fmax, fmin = Harm1.max(), Harm1.min()
fcolors1 = (Harm1 - fmin)/(fmax - fmin)

fmax, fmin = Harm2.max(), Harm2.min()
fcolors2 = (Harm2 - fmin)/(fmax - fmin)

fmax, fmin = Harm3.max(), Harm3.min()
fcolors3 = (Harm3 - fmin)/(fmax - fmin)

fmax, fmin = Harm4.max(), Harm4.min()
fcolors4 = (Harm4 - fmin)/(fmax - fmin)

fmax, fmin = Harm5.max(), Harm5.min()
fcolors5 = (Harm5 - fmin)/(fmax - fmin)

fmax, fmin = Harm6.max(), Harm6.min()
fcolors6 = (Harm6 - fmin)/(fmax - fmin)

fmax, fmin = Harm7.max(), Harm7.min()
fcolors7 = (Harm7 - fmin)/(fmax - fmin)

fmax, fmin = Harm8.max(), Harm8.min()
fcolors8 = (Harm8 - fmin)/(fmax - fmin)


# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(1))


# first plot
ax = fig.add_subplot(3, 3, 1, projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors0))
ax.set_axis_off()

#second plot
ax = fig.add_subplot(3, 3, 2, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors1))
ax.set_axis_off()

#third plot
ax = fig.add_subplot(3, 3, 3, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors2))
ax.set_axis_off()

#fourth plot
ax = fig.add_subplot(3, 3, 4, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors3))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 5, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors4))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 6, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors5))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 7, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors6))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 8, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors7))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 9, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors8))
ax.set_axis_off()



plt.show()


#--------------------------#
# Reference Values               #
#--------------------------#
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

m, l = -1, 1

# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = sph_harm(m, l, theta, phi).real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes
ax.set_axis_off()
plt.show()



# setup colors
fcolors00 = sph_harm(0, 0 , theta, phi).real
fmax, fmin = fcolors00.max(), fcolors00.min()
#fcolors00 = (fcolors00 - fmin)/(fmax - fmin)

fcolors01 = sph_harm(-1, 1 , theta, phi).real
fmax, fmin = fcolors01.max(), fcolors01.min()
fcolors01 = (fcolors01 - fmin)/(fmax - fmin)

fcolors02 = sph_harm(0, 1 , theta, phi).real
fmax, fmin = fcolors02.max(), fcolors02.min()
fcolors02 = (fcolors02 - fmin)/(fmax - fmin)

fcolors03 = sph_harm(1, 1 , theta, phi).real
fmax, fmin = fcolors03.max(), fcolors03.min()
fcolors03 = (fcolors03 - fmin)/(fmax - fmin)

fcolors04 = sph_harm(-2, 2 , theta, phi).real
fmax, fmin = fcolors04.max(), fcolors04.min()
fcolors04 = (fcolors04 - fmin)/(fmax - fmin)

fcolors05 = sph_harm(-1, 2 , theta, phi).real
fmax, fmin = fcolors05.max(), fcolors05.min()
fcolors05 = (fcolors05 - fmin)/(fmax - fmin)

fcolors06 = sph_harm(0, 2 , theta, phi).real
fmax, fmin = fcolors06.max(), fcolors06.min()
fcolors06 = (fcolors06 - fmin)/(fmax - fmin)

fcolors07 = sph_harm(1, 2 , theta, phi).real
fmax, fmin = fcolors07.max(), fcolors07.min()
fcolors07 = (fcolors07 - fmin)/(fmax - fmin)

fcolors08 = sph_harm(2, 2 , theta, phi).real
fmax, fmin = fcolors08.max(), fcolors08.min()
fcolors08 = (fcolors08 - fmin)/(fmax - fmin)


# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(1))

# first plot
ax = fig.add_subplot(3, 3, 1, projection='3d')
surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors00))
ax.set_axis_off()

#second plot
ax = fig.add_subplot(3, 3, 2, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors01))
ax.set_axis_off()

#third plot
ax = fig.add_subplot(3, 3, 3, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors02))
ax.set_axis_off()

#fourth plot
ax = fig.add_subplot(3, 3, 4, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors03))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 5, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors04))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 6, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors05))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 7, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors06))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 8, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors07))
ax.set_axis_off()

ax = fig.add_subplot(3, 3, 9, projection='3d')
ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors08))
ax.set_axis_off()



plt.show()
