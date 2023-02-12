"""
********************************************************************************
2D lid-driven cavity flow with FDM
    space: KK scheme for convection, 2nd-order central for diffusion
    time: fractional step method
********************************************************************************
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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

    v_tol = 1e-6          # tolerance for velocity residual
    p_tol = 1e-6          # tolerance for ppe residual
    n_max = int(1e8)      # max it for main loop
    m_max = int(1e5)      # max it for ppe loop
    cmap = plt.cm.turbo   # visualization
    nbar = 21             # visualization

    print(f"\n******************************************************************")
    print(f">>>>> dt1: {dt1:.4e}, dt2: {dt2:.4e}, dt3: {dt3:.4e}, dt: {dt:.4e}")
    print(f">>>>> dx: {dx:.4e}, dy: {dy:.4e}, dt: {dt:.4e}, Re: {Re:.1f}")
    print(f"******************************************************************")

    # variables
    u = np.zeros(shape=(ny, nx)) + 1e-6
    v = np.zeros_like(u) + 1e-6
    p = np.zeros_like(u) + 1e-6
    b = np.zeros_like(u) + 1e-6

    # log
    it_log = []
    u_res_log = []
    v_res_log = []
    vel_res_log = []   # vel norm

    # main
    t = 0.
    for n in range(0, n_max+1):   # for is faster than while
        t += dt

        # previous quantities
        u_    = np.copy(u)   # previous vel
        v_    = np.copy(v)
        u_hat = np.copy(u)   # intermediate vel
        v_hat = np.copy(v)

        # intermediate vel (KK scheme for convection)
        u_hat[2:-2, 2:-2] = u_[2:-2, 2:-2] \
                            - dt / 12. / dx * u_[2:-2, 2:-2] * (u_[2:-2, :-4] - 8. * u_[2:-2, 1:-3] + 8. * u_[2:-2, 3:-1] - u_[2:-2, 4:]) \
                            - beta * dt / dx * np.abs(u_[2:-2, 2:-2]) * (u_[2:-2, :-4] - 4. * u_[2:-2, 1:-3] + 6. * u_[2:-2, 2:-2] - 4. * u_[2:-2, 3:-1] + u_[2:-2, 4:]) \
                            - dt / 12. / dy * v_[2:-2, 2:-2] * (u_[:-4, 2:-2] - 8. * u_[1:-3, 2:-2] + 8. * u_[3:-1, 2:-2] - u_[4:, 2:-2]) \
                            - beta * dt / dy * np.abs(v_[2:-2, 2:-2]) * (u_[:-4, 2:-2] - 4. * u_[1:-3, 2:-2] + 6. * u_[2:-2, 2:-2] - 4. * u_[3:-1, 2:-2] + u_[4:, 2:-2]) \
                            + 1. / Re * dt / dx ** 2 * (u_[2:-2, 3:-1] - 2. * u_[2:-2, 2:-2] + u_[2:-2, 1:-3]) \
                            + 1. / Re * dt / dy ** 2 * (u_[3:-1, 2:-2] - 2. * u_[2:-2, 2:-2] + u_[1:-3, 2:-2])
        v_hat[2:-2, 2:-2] = v_[2:-2, 2:-2] \
                            - dt / 12. / dx * u_[2:-2, 2:-2] * (v_[2:-2, :-4] - 8. * v_[2:-2, 1:-3] + 8. * v_[2:-2, 3:-1] - v_[2:-2, 4:]) \
                            - beta * dt / dx * np.abs(u_[2:-2, 2:-2]) * (v_[2:-2, :-4] - 4. * v_[2:-2, 1:-3] + 6. * v_[2:-2, 2:-2] - 4. * v_[2:-2, 3:-1] + v_[2:-2, 4:]) \
                            - dt / 12. / dy * v_[2:-2, 2:-2] * (v_[:-4, 2:-2] - 8. * v_[1:-3, 2:-2] + 8. * v_[3:-1, 2:-2] - v_[4:, 2:-2]) \
                            - beta * dt / dy * np.abs(v_[2:-2, 2:-2]) * (v_[:-4, 2:-2] - 4. * v_[1:-3, 2:-2] + 6. * v_[2:-2, 2:-2] - 4. * v_[3:-1, 2:-2] + v_[4:, 2:-2]) \
                            + 1. / Re * dt / dx ** 2 * (v_[2:-2, 3:-1] - 2. * v_[2:-2, 2:-2] + v_[2:-2, 1:-3]) \
                            + 1. / Re * dt / dy ** 2 * (v_[3:-1, 2:-2] - 2. * v_[2:-2, 2:-2] + v_[1:-3, 2:-2])

        # boundary condition (KK scheme requires long stencil)
        u_hat[0,  :] = 0.; v_hat[0,  :] = 0.   # South
        u_hat[1,  :] = 0.; v_hat[1,  :] = 0.   # South
        u_hat[:,  0] = 0.; v_hat[:,  0] = 0.   # West
        u_hat[:,  1] = 0.; v_hat[:,  1] = 0.   # West
        u_hat[:, -1] = 0.; v_hat[:, -1] = 0.   # East
        u_hat[:, -2] = 0.; v_hat[:, -2] = 0.   # East
        u_hat[-1, :] = 1.; v_hat[-1, :] = 0.   # North (driving lid)
        u_hat[-2, :] = 1.; v_hat[-2, :] = 0.   # North

        # source term for PPE (fractional step)
        b[1:-1, 1:-1] = 1. / dt * (
                              (u_hat[1:-1, 2:] - u_hat[1:-1, :-2]) / (2. * dx) \
                            + (v_hat[2:, 1:-1] - v_hat[:-2, 1:-1]) / (2. * dy)
        )

        # solve PPE (Jacobi)
        t0 = time.perf_counter()
        for m in range(0, m_max+1):
            # update pressure
            p_ = np.copy(p)
            p[1:-1, 1:-1] = 1. / (2. * (dx ** 2 + dy ** 2)) \
                            * (
                                  (p_[1:-1, 2:] + p_[1:-1, :-2]) * dy ** 2 \
                                + (p_[2:, 1:-1] + p_[:-2, 1:-1]) * dx ** 2 \
                                - b[1:-1, 1:-1] * (dx ** 2) * (dy ** 2)
                            )

            # boundary condition
            p[0,  :] = p[1,  :]         # South (Neu)
            p[-1, :] = p[-2, :]         # North (Neu)
            p[:,  0] = p[:,  1]         # West  (Neu)
            p[:, -1] = p[:, -2]         # East  (Neu)
            p[0, int(nx / 2)] = 0.   # South (Dir)

            # # boundary condition
            # p[1,  :] = p[2,  :]         # South (Neu)
            # p[-2, :] = p[-3, :]         # North (Neu)
            # p[:,  1] = p[:,  2]         # West  (Neu)
            # p[:, -2] = p[:, -3]         # East  (Neu)
            # p[1, int(nx / 2)] = 0.   # South (Dir)

            # # boundary condition
            # p[0,  :] = p[2,  :]; p[1,  :] = p[2,  :]         # South (Neu)
            # p[-1, :] = p[-3, :]; p[-2, :] = p[-3, :]         # North (Neu)
            # p[:,  0] = p[:,  2]; p[:,  1] = p[:,  2]         # West  (Neu)
            # p[:, -1] = p[:, -3]; p[:, -2] = p[:, -3]         # East  (Neu)
            # p[0, int(nx / 2)] = 0.; p[1, int(nx / 2)] = 0.   # South (Dir)

            # converged?
            p_res = np.sqrt(np.sum((p - p_) ** 2)) / np.sqrt(np.sum(p_ ** 2))   # L2 norm
            if m % int(m_max * .05) == 0:
                t1 = time.perf_counter()
                elps = t1 - t0
                print("   >>> ppe -> loop: %d, p_res: %.6e, elps: %.1f" % (m, p_res, elps))
            if p_res < p_tol:
                print("   >>> ppe converged")
                break

        # update velocity (KK scheme)
        u[2:-2, 2:-2] = u_hat[2:-2, 2:-2] \
                        - dt / (2. * dx) * (p[2:-2, 3:-1] - p[2:-2, 1:-3])
        v[2:-2, 2:-2] = v_hat[2:-2, 2:-2] \
                        - dt / (2. * dy) * (p[3:-1, 2:-2] - p[1:-3, 2:-2])

        # boundary condition (KK scheme)
        u[0,  :] = 0.; v[0,  :] = 0.   # South
        u[1,  :] = 0.; v[1,  :] = 0.   # South
        u[:,  0] = 0.; v[:,  0] = 0.   # West
        u[:,  1] = 0.; v[:,  1] = 0.   # West
        u[:, -1] = 0.; v[:, -1] = 0.   # East
        u[:, -2] = 0.; v[:, -2] = 0.   # East
        u[-1, :] = 1.; v[-1, :] = 0.   # North (driving lid)
        u[-2, :] = 1.; v[-2, :] = 0.   # North

        # converged?
        u_res = np.sqrt(np.sum((u - u_) ** 2))   # L2 norm
        v_res = np.sqrt(np.sum((v - v_) ** 2))
        vel_res = np.maximum(u_res, v_res)
        print("\n***************************************************************")
        print(">>> main -> loop: %d, t: %.3f, vel_res: %.6e" % (n, t, vel_res))
        print("***************************************************************")
        if vel_res < v_tol or n == n_max:
            print(">>> main converged")
            np.save("./res_x.npy", x)
            np.save("./res_y.npy", y)
            np.save("./res_X.npy", X)
            np.save("./res_Y.npy", Y)
            np.save("./res_u.npy", u)
            np.save("./res_v.npy", v)
            np.save("./res_p.npy", p)
            break

        it_log.append(n)
        u_res_log.append(u_res)
        v_res_log.append(v_res)
        vel_res_log.append(vel_res)

        # visualization
        if n % 1000 == 0:
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
            plt.title(f"Courant number (Re={Re:.1f})")
            plt.tight_layout()
            plt.savefig("./cfl.png", dpi=300)
            plt.clf()
            plt.close()

            plt.figure(figsize=(5, 4))
            div = np.zeros_like(u)
            div[1:-1, 1:-1] = (u[1:-1, 2:] - u[1:-1, :-2]) / 2. / dx \
                            + (v[2:, 1:-1] - v[:-2, 1:-1]) / 2. / dy
            limit = np.abs(div).max()
            scale = .1
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
            plt.title(f"velocity divergence (Re={Re:.1f})")
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

            plt.figure(figsize=(5, 4))
            skip = 500
            plt.plot(it_log[::skip], u_res_log[::skip], marker=".", label="u_res")
            plt.plot(it_log[::skip], v_res_log[::skip], marker=".", label="v_res")
            plt.plot(it_log[::skip], vel_res_log[::skip], marker=".", label="vel_res", c="k")
            plt.grid(alpha=.3)
            plt.legend(loc="upper right")
            plt.yscale("log")
            plt.xlabel("it")
            plt.ylabel("velocity change rate")
            plt.tight_layout()
            plt.savefig("./vel_res.png", dpi=300)
            plt.clf()
            plt.close()

if __name__ == "__main__":
    main()
