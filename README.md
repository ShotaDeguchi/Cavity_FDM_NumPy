# Cavity_FDM_NumPy

Updated version of this repository can be found [here](https://github.com/ShotaDeguchi/Cavity_FDM_NumPy2). 

FDM (Finite Difference Method) simulation of 2D lid-driven cavity flow based on :
* fractional step method for time integration
* [Kawamura-Kuwahara scheme](https://doi.org/10.2514/6.1984-340) (3rd-order upwind -> 4th-order central with 4th-order numerical viscosity) for convection
* 2nd-order central difference for pressure gradient and viscosity terms

The results are compared with the reference solution (for velocity) presented in [Ghia+1986](https://doi.org/10.1016/0021-9991(82)90058-4). 

## Results
Cavity flow is a steady problem. We consider that the field has reached to a steady state when the following is satisfied:
```math
\max \left( \frac{\| u^{(n+1)} - u^{(n)} \|_2}{\| u^{(n)} \|_2}, \frac{\| v^{(n+1)} - v^{(n)} \|_2}{\| v^{(n)} \|_2} \right) \le \delta
```
where $\delta$ is the convergence tolerance, set to $\delta = 10^{-6}$. 

The following summarizes results at different Reynolds numbers and different resolutions. 

| Column name | Description | 
| :---: | :--- |
| Re | Reynolds number (inertia vs viscosity) |
| t | Dimensionless time until the convergence (when velocity residual $\le \delta$ is met) |
| u | Horizontal velocity along the geometric center |
| v | Vertical velocity along the geometric center |

### $\Delta x = \Delta y = 5 \times 10^{-3}$
| Re | t | Velocity norm  | Pressure | u | v |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 100 | 15.4 | ![](Re_100/vel_norm.png) | ![](Re_100/prs.png) | ![](Re_100/comparison_u.png) | ![](Re_100/comparison_v.png) |
| 400 | 26.8 | ![](Re_400/vel_norm.png) | ![](Re_400/prs.png) | ![](Re_400/comparison_u.png) | ![](Re_400/comparison_v.png) |
| 1,000 | 36.4 | ![](Re_1000/vel_norm.png) | ![](Re_1000/prs.png) | ![](Re_1000/comparison_u.png) | ![](Re_1000/comparison_v.png) |
| 3,200 | 87.5 | ![](Re_3200/vel_norm.png) | ![](Re_3200/prs.png) | ![](Re_3200/comparison_u.png) | ![](Re_3200/comparison_v.png) |
| 5,000 | 148.5 | ![](Re_5000/vel_norm.png) | ![](Re_5000/prs.png) | ![](Re_5000/comparison_u.png) | ![](Re_5000/comparison_v.png) |

### $\Delta x = \Delta y = 2 \times 10^{-3}$
| Re | t | Velocity norm  | Pressure | u | v |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 100 | 13.7 | ![](Re_100_highres/vel_norm.png) | ![](Re_100_highres/prs.png) | ![](Re_100_highres/comparison_u.png) | ![](Re_100_highres/comparison_v.png) |
| 400 | 19.2 | ![](Re_400_highres/vel_norm.png) | ![](Re_400_highres/prs.png) | ![](Re_400_highres/comparison_u.png) | ![](Re_400_highres/comparison_v.png) |
| 1,000 | 30.7 | ![](Re_1000_highres/vel_norm.png) | ![](Re_1000_highres/prs.png) | ![](Re_1000_highres/comparison_u.png) | ![](Re_1000_highres/comparison_v.png) |
| 3,200 | 68.4 | ![](Re_3200_highres/vel_norm.png) | ![](Re_3200_highres/prs.png) | ![](Re_3200_highres/comparison_u.png) | ![](Re_3200_highres/comparison_v.png) |
| 5,000 | 134.1 | ![](Re_5000_highres/vel_norm.png) | ![](Re_5000_highres/prs.png) | ![](Re_5000_highres/comparison_u.png) | ![](Re_5000_highres/comparison_v.png) |

## Requirements
Tested environment:
* numpy == 1.22.4
* matplotlib == 3.5.2

## License
MIT License
