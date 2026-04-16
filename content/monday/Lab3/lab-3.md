---
weight: 3
author: Niall Miller
---

# Lab 3: Asteroseismology

In this lab we will use MESA to compute asteroseismic observables and explore how they constrain stellar ages. The three quantities we focus on are the large frequency separation $\Delta\nu$, the small frequency separation $\delta\nu_{02}$, and the g-mode period spacing $\Delta\Pi_g$. All of the code that computes these is already implemented in `run_star_extras.f90` — you do not need to modify it. Your job is to run the models, understand the physics behind each quantity, and combine the results across the group.

> [!NOTE]
> This lab builds on the models from Labs 1 and 2. Make sure your `LOGS/history.data` is present before starting.

## Background

Stars oscillate in resonant modes. The frequencies of these modes encode information about the stellar interior — in particular, the sound speed profile and the Brunt–Väisälä frequency. By measuring just a few summary statistics of the oscillation spectrum, we can constrain both the global structure and the interior chemical composition of a star.

## Setup

Your work directory already contains a pre-compiled MESA model and a `run_star_extras.f90` with all the seismic quantities implemented. Copy it to your working area and run it:

```bash
cp -r Lab3 my_lab3
cd my_lab3
./rn
```

The extra history columns written by `run_star_extras` will appear automatically in `LOGS/history.data`. You do not need to add anything to `history_columns.list` for these.

> [!IMPORTANT]
> Do **not** edit `src/run_star_extras.f90`. All the physics is already there. If you need to recompile, run `./mk` without making any changes to the source.

---

## The large frequency separation $\Delta\nu$

The large frequency separation is the average spacing between p-modes of the same angular degree $\ell$ and consecutive radial order $n$. It is related to the sound crossing time of the star:

$$ \Delta\nu = \left( 2 \int_0^R \frac{dr}{c_s} \right)^{-1} $$

```fortran
Delta_nu_int = 0.
do k = 2, s% nz, 1
   dr = s% rmid(k-1) - s% rmid(k)
   Delta_nu_int = Delta_nu_int + dr/s% csound(k)
end do
Delta_nu_int = 1./(2.*Delta_nu_int)
```

The loop accumulates $dr/c_s$ from the centre to the surface, building up the total sound travel time. The final line takes the reciprocal and divides by 2 to give $\Delta\nu$ in Hz. Because $\Delta\nu \propto \sqrt{M/R^3}$, it is a measure of the mean stellar density — it decreases as the star expands during main sequence evolution.

---

## The small frequency separation $\delta\nu_{02}$

The small frequency separation is the offset between $\ell=0$ and $\ell=2$ modes of similar frequency. Unlike $\Delta\nu$, it is sensitive to the gradient of the sound speed in the stellar core, and therefore directly to the central hydrogen abundance:

$$ \delta\nu_{02} \approx -\frac{2\,\Delta\nu}{\pi^2\,\nu_\mathrm{max}} \left( \int_0^R \frac{1}{r} \frac{dc_s}{dr} \, dr - \frac{c_s(R)}{R} \right) $$

```fortran
delta_nu02_int = 0.
nu_max = s% nu_max * 1d-6
do k = 2, s% nz, 1
   dc = s% csound(k-1) - s% csound(k)
   delta_nu02_int = delta_nu02_int + dc / s% rmid(k)
end do
delta_nu02_int = delta_nu02_int - s% csound(1)/s% r(1)
delta_nu02_int = -2 * Delta_nu_int * delta_nu02_int / (3.1415926535**2. * nu_max)
```

The loop computes $dc_s/r$ shell by shell, approximating the radial derivative of the sound speed as a finite difference (`dc = csound(k-1) - csound(k)`). The surface boundary term $c_s(R)/R$ is subtracted after the loop. The whole integral is then scaled by $-2\Delta\nu/(\pi^2\nu_\mathrm{max})$ to give $\delta\nu_{02}$ in Hz. Because the sound speed gradient in the core steepens as hydrogen is depleted, $\delta\nu_{02}$ decreases monotonically with age — making it a direct age clock.

---

## The g-mode period spacing $\Delta\Pi_g$

G-modes are buoyancy-driven oscillations confined to the radiative interior. Their period spacing is set by the integral of the Brunt–Väisälä frequency $N$ over the g-mode cavity:

$$ \Delta\Pi_g = \frac{\pi^2 \sqrt{2}}{\displaystyle\int_\mathrm{cavity} \frac{N}{r} \, dr} $$

```fortran
delta_Pg_int = 0.
entered_g_mode_cavity = .false.
do k = s% nz, 2, -1
   dr = s% rmid(k-1) - s% rmid(k)
   if (s% brunt_N2(k) > 0) then
      entered_g_mode_cavity = .true.
      delta_Pg_int = delta_Pg_int + sqrt(s% brunt_N2(k)) / s% r(k) * dr
   else
      if (entered_g_mode_cavity) exit
   end if
end do
delta_Pg_int = (3.1415926535**2.) * sqrt(2.0) / delta_Pg_int
```

The loop walks inward from the surface (index `nz` to `2`). It only accumulates $\sqrt{N^2}/r \cdot dr$ where $N^2 > 0$, i.e. inside the radiative g-mode cavity. The `entered_g_mode_cavity` flag ensures the loop exits as soon as it leaves the cavity on the inner side, preventing the convective core from being accidentally included. The final line applies the prefactor $\pi^2\sqrt{2}$ to give $\Delta\Pi_g$ in seconds.

---

## Exercise 1: Reading the seismic history columns

Once your run has finished, open `LOGS/history.data` and locate the columns `Delta_nu_int`, `delta_nu02_int`, and `delta_Pg_int`. Plot each as a function of stellar age.

{{< details title="Hint: Python snippet to plot the extra columns" closed="true" >}}

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('LOGS/history.data', names=True, skip_header=5)

fig, axes = plt.subplots(3, 1, sharex=True)
age = data['star_age'] / 1e9

axes[0].plot(age, data['Delta_nu_int'] * 1e6)
axes[0].set_ylabel(r'$\Delta\nu$ ($\mu$Hz)')

axes[1].plot(age, data['delta_nu02_int'] * 1e6)
axes[1].set_ylabel(r'$\delta\nu_{02}$ ($\mu$Hz)')

axes[2].plot(age, data['delta_Pg_int'])
axes[2].set_ylabel(r'$\Delta\Pi_g$ (s)')
axes[2].set_xlabel('Age (Gyr)')

plt.tight_layout()
plt.show()
```

{{< /details >}}

Which quantity changes fastest with age? Which would give the tightest age constraint for a star of your mass?

## Exercise 2: Crowd-source the seismic grid

Enter your results at 1 Gyr into the shared spreadsheet:

| Column | Value |
|--------|-------|
| `initial_mass` | your assigned mass ($M_\odot$) |
| $\Delta\nu$ at 1 Gyr | `Delta_nu_int` in $\mu$Hz |
| $\delta\nu_{02}$ at 1 Gyr | `delta_nu02_int` in $\mu$Hz |
| $\Delta\Pi_g$ at 1 Gyr | `delta_Pg_int` in s |

Once everyone has contributed, we will look at the full grid and discuss how well the seismic observables separate stars of different masses at the same age, and how this compares to the CMD separation from Lab 2.

## Bonus: Combined age constraints

Using both the colour–magnitude information from Lab 2 and the seismic observables from this lab, consider which combination of observables gives the tightest age constraint for a solar-type star. Which observable is least degenerate with metallicity?

{{< details title="Discussion prompt" closed="true" >}}

Think about the following: $\Delta\nu$ scales with mean density, so it is sensitive to both mass and radius. $\delta\nu_{02}$ is sensitive to the core sound speed gradient, which depends on the central hydrogen abundance — almost independent of the envelope. How does adding $\delta\nu_{02}$ to an isochrone fit change the mass–age degeneracy?

{{< /details >}}
