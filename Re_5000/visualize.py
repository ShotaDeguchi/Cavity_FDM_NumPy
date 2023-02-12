"""
visualize the results
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyevtk


def main():
    # discretization (regular grid)
    xmin, xmax = 0., 1.
    ymin, ymax = 0., 1.
    dx, dy = 5e-3, 5e-3

    xmin, xmax = xmin - dx, xmax + dx   # long stencil required in KK scheme
    ymin, ymax = ymin - dy, ymax + dy

    x = np.arange(xmin, xmax+1e-9, dx)
    y = np.arange(ymin, ymax+1e-9, dy)
    nx, ny = len(x), len(y)
    X, Y = np.meshgrid(x, y)

    # parameters
    Re = 3200.       # Reynolds number
    k = 1. / Re      # diffusion rate (as in diffusion eq)
    beta = 1. / 4.   # coef for numerical viscosity in KK scheme

    # timestep
    U, V = 1., 1.   # characteristic velocity
    dt1 = dx * dy / (V * dx + U * dy)    # CFL
    dt2 = .5 / k * (dx**2 * dy**2) / (dx**2 + dy**2)   # diffusion
    dt = min(dt1, dt2)   # strictest  criterion
    dt *= .8             # safety

    # load
    x = np.load("./res_x.npy")
    y = np.load("./res_y.npy")
    X = np.load("./res_X.npy")
    Y = np.load("./res_Y.npy")
    u = np.load("./res_u.npy")
    v = np.load("./res_v.npy")
    p = np.load("./res_p.npy")

    # visualize
    cmap = plt.cm.turbo
    nbar = 21

    plt.figure(figsize=(5, 4))
    vmin, vmax = -.2, 1.
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(X, Y, c=u, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="u")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./vel_x.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    vmin, vmax = -.5, .5
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(X, Y, c=v, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="v")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./vel_y.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    vmin, vmax = 0., 1.
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    vel_norm = np.sqrt(u ** 2 + v ** 2)
    plt.scatter(X, Y, c=vel_norm, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="vel norm")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./vel_norm.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    vmin, vmax = 0., 1.
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.streamplot(
        X, Y, u, v, color=vel_norm, cmap=cmap, norm=norm
    )
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="vel norm")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./psi.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    limit = np.abs(p).max()
    scale = .15
    vmin, vmax = - scale * limit, scale * limit
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(X, Y, c=p, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="p")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./prs.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    C = dt / dx * vel_norm
    vmin, vmax = 0., C.max()
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(X, Y, c=C, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="cfl")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./cfl.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    div = np.zeros_like(u)
    div[1:-1, 1:-1] = (u[1:-1, 2:] - u[1:-1, :-2]) / 2. / dx \
                    + (v[2:, 1:-1] - v[:-2, 1:-1]) / 2. / dy
    limit = np.abs(div).max()
    scale = .05
    vmin, vmax = - scale * limit, scale * limit
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar-1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(X, Y, c=div, cmap="seismic", norm=norm, s=.3)
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="div")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    # plt.xlim(.75, 1.)
    # plt.ylim(.75, 1.)
    # plt.tick_params(left=False, right=False, bottom=False, top=False, 
    #                 labelleft=False, labelbottom=False)
    plt.tight_layout()
    plt.savefig("./div.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(4, 4))
    # Ghia+1982: horizontal velocity along the geometric center
    y_Ghia = [1., 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.]
    u_Ghia = [1., 0.48223, 0.4612, 0.45992, 0.46036, 0.33556, 0.20087, 0.08183, -0.03039, -0.07404, -0.22855, -0.3305, -0.40435, -0.43643, -0.42901, -0.41165, 0.]
    plt.scatter(u_Ghia[::-1], y_Ghia[::-1], c="k", marker="s", label="Ghia+1982")
    plt.plot(u[:, int(nx / 2)], Y[:, int(nx / 2)], c="r", ls="--", label="FDM")
    plt.legend(loc="best")
    plt.xlabel("u")
    plt.ylabel("y")
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig("./comparison_u.svg")
    plt.clf()
    plt.close()

    plt.figure(figsize=(4, 4))
    # Ghia+1982: vertical velocity along the geometric center
    x_Ghia = [1., 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.]
    v_Ghia = [0., -0.49774, -0.55069, -0.55408, -0.52876, -0.41442, -0.36214, -0.30018, 0.00945, 0.2728, 0.28066, 0.35368, 0.42951, 0.43648, 0.43329, 0.42447, 0.]
    plt.scatter(x_Ghia[::-1], v_Ghia[::-1], c="k", marker="s", label="Ghia+1982")
    plt.plot(X[int(nx / 2), :], v[int(nx / 2), :], c="r", ls="--", label="FDM")
    plt.legend(loc="best")
    plt.xlabel("x")
    plt.ylabel("v")
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig("./comparison_v.svg")
    plt.clf()
    plt.close()

    # # visualize as vtk (vtr format)
    # X = X[:,:,np.newaxis]   # add z-dimension
    # Y = Y[:,:,np.newaxis]
    # Z = np.zeros_like(X)
    # U = u[:,:,np.newaxis]
    # pyevtk.hl.gridToVTK("./output/u_numpy", 
    #                     X, Y, Z, 
    #                     cellData  = {"u": U}, 
    #                     pointData = {"u": U})

if __name__ == "__main__":
    main()