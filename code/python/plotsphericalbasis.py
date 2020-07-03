import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm

# nur fuer den Seiteneffekt: plt.gca(projection = '3d') funktioniert sonst nicht
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

theta_1d = np.linspace(0, np.pi, 91)  # 2 GRAD Schritte
phi_1d = np.linspace(0, 2 * np.pi, 181)  # 2 GRAD Schritte

theta_2d, phi_2d = np.meshgrid(theta_1d, phi_1d)
xyz_2d = np.array([np.sin(theta_2d) * np.sin(phi_2d),
                   np.sin(theta_2d) * np.cos(phi_2d),
                   np.cos(theta_2d)])

colormap = cm.ScalarMappable(cmap=plt.get_cmap("cool"))
colormap.set_clim(-.45, .45)
limit = .5


def show_Y_lm(l, m):
    print("Y_%i_%i" % (l, m))  # zeigen, dass was passiert
    plt.figure()
    ax = plt.gca(projection="3d")

    plt.title("$Y^{%i}_{%i}$" % (m, l))
    Y_lm = sph_harm(m, l, phi_2d, theta_2d)
    r = np.abs(Y_lm.real) * xyz_2d
    ax.plot_surface(r[0], r[1], r[2],
                    facecolors=colormap.to_rgba(Y_lm.real),
                    rstride=2, cstride=2)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    #ax.set_aspect("equal")
    # ax.set_axis_off()


# Vorsicht: diese Schleifen erzeugen 16 plots (in 16 Fenstern)!
for l in range(0, 4):
    for m in range(-l, l + 1):
        show_Y_lm(l, m)

show_Y_lm(l=5, m=0)
show_Y_lm(l=5, m=4)
show_Y_lm(l=6, m=6)

plt.show()