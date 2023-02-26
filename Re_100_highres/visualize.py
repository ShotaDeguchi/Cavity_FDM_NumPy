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
    dx, dy = 2e-3, 2e-3

    xmin, xmax = xmin - dx, xmax + dx   # long stencil required in KK scheme
    ymin, ymax = ymin - dy, ymax + dy

    x = np.arange(xmin, xmax+1e-9, dx)
    y = np.arange(ymin, ymax+1e-9, dy)
    nx, ny = len(x), len(y)
    X, Y = np.meshgrid(x, y)

    # parameters
    Re = 100.        # Reynolds number
    k = 1. / Re      # diffusion rate (as in diffusion eq)
    beta = 1. / 4.   # coef for numerical viscosity in KK scheme

    # timestep
    U, V = 1., 1.   # characteristic velocity
    dt1 = dx * dy / (V * dx + U * dy)    # CFL
    dt2 = dx**2 / 2. / k                 # CFL
    dt3 = .5 / k * (dx**2 * dy**2) / (dx**2 + dy**2)   # diffusion
    dt = min(dt1, dt2, dt3)   # strictest  criterion
    dt *= .5             # safety

    cmap = plt.cm.turbo   # visualization
    nbar = 21             # visualization

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
    plt.title(f"u (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./vel_x.png", dpi=300)
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
    plt.title(f"v (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./vel_y.png", dpi=300)
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
    plt.title(f"vel norm (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./vel_norm.png", dpi=300)
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    vmin, vmax = 0., 1.
    vticks = (vmax - vmin) / 4
    bounds = np.linspace(vmin, vmax, nbar)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plt.streamplot(
        X, Y, u, v, color=vel_norm, cmap=cmap, norm=norm, density=1.5
    )
    plt.colorbar(ticks=np.arange(vmin, vmax+1e-9, vticks), label="vel norm")
    plt.xticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.yticks(ticks=np.arange(0., 1.+1e-9, .25))
    plt.xlim(0., 1.)
    plt.ylim(0., 1.)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"psi (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./psi.png", dpi=300)
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    limit = .15
    vmin, vmax = - limit, limit
    # limit = np.abs(p).max()
    # scale = .15
    # vmin, vmax = - scale * limit, scale * limit
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
    plt.title(f"p (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./prs.png", dpi=300)
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
    plt.title(f"cfl number (Re={Re:.1f})")
    plt.tight_layout()
    plt.savefig("./cfl.png", dpi=300)
    plt.clf()
    plt.close()

    plt.figure(figsize=(5, 4))
    div = np.zeros_like(u)
    div[1:-1, 1:-1] = (u[1:-1, 2:] - u[1:-1, :-2]) / 2. / dx \
                    + (v[2:, 1:-1] - v[:-2, 1:-1]) / 2. / dy
    # limit = np.abs(div).max()
    # scale = .05
    # vmin, vmax = - scale * limit, scale * limit
    limit = .05
    vmin, vmax = - limit, limit
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
    plt.title(f"vel div (Re={Re:.1f})")
    # plt.xlim(.75, 1.)
    # plt.ylim(.75, 1.)
    # plt.tick_params(left=False, right=False, bottom=False, top=False, 
    #                 labelleft=False, labelbottom=False)
    plt.tight_layout()
    plt.savefig("./div.png", dpi=300)
    plt.clf()
    plt.close()

    plt.figure(figsize=(4, 4))
    # Ghia+1982: horizontal velocity along the geometric center
    y_Ghia = [1., 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.]
    u_Ghia = [1., 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.2109, -0.15662, -0.1015, -0.06434, -0.04775, -0.04192, -0.03717, 0.]
    plt.scatter(u_Ghia[::-1], y_Ghia[::-1], c="k", marker="s", label="Ghia+1982")
    plt.plot(u[:, int(nx / 2)], Y[:, int(nx / 2)], c="r", ls="--", label="FDM")
    plt.legend(loc="best")
    plt.xlabel("u")
    plt.ylabel("y")
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig("./comparison_u.png", dpi=300)
    plt.clf()
    plt.close()

    plt.figure(figsize=(4, 4))
    # Ghia+1982: vertical velocity along the geometric center
    x_Ghia = [1., 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.]
    v_Ghia = [0., -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.1089, 0.10091, 0.09233, 0.]
    plt.scatter(x_Ghia[::-1], v_Ghia[::-1], c="k", marker="s", label="Ghia+1982")
    plt.plot(X[int(nx / 2), :], v[int(nx / 2), :], c="r", ls="--", label="FDM")
    plt.legend(loc="best")
    plt.xlabel("x")
    plt.ylabel("v")
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig("./comparison_v.png", dpi=300)
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
