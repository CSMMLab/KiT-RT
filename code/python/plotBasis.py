import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm

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



# read in own coordinates
import csv

ndimGrid = 51
X = np.zeros((ndimGrid, ndimGrid))
Y = np.zeros((ndimGrid, ndimGrid))
Z = np.zeros((ndimGrid, ndimGrid))

Harm0 = np.zeros((ndimGrid, ndimGrid))
Harm1 = np.zeros((ndimGrid, ndimGrid))
Harm2 = np.zeros((ndimGrid, ndimGrid))
Harm3 = np.zeros((ndimGrid, ndimGrid))

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

            line_count += 1

    print(f'Processed {line_count} lines.')

#plot data
fmax, fmin = Harm2.max(), Harm2.min()
fcolors0 = (Harm2 - fmin)/(fmax - fmin)

#fig = plt.figure(figsize=plt.figaspect(1.))
#ax = fig.add_subplot( projection='3d')

fig, ax = plt.subplots(2,projection='3d' )
ax[0].plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors0))
# Turn off the axis planes
ax[0].set_axis_off()

fmax, fmin = Harm1.max(), Harm1.min()
fcolors1 = (Harm1 - fmin)/(fmax - fmin)

#ax2 = fig.add_subplot( projection='3d')
ax[1].plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors1))
# Turn off the axis planes
ax[1].set_axis_off()
plt.show()
