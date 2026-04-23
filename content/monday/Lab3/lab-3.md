---
weight: 3
author: Niall Miller
---

# Lab 3: Asteroseismology

In this lab we will use MESA to compute asteroseismic observables and explore how they constrain stellar ages. The two quantities we focus on are the large frequency separation $\Delta\nu$ and the small frequency separation $\delta\nu_{02}$. All of the code that computes these is already implemented in `run_star_extras.f90` — you do not need to modify it. Your job is to run the models, understand the physics behind each quantity, and combine the results across the group.

> [!NOTE]
> This lab builds on the models from Labs 1 and 2. Your `inlist_run` is pre-configured — check that your assigned mass is set correctly before starting.

## Background

<!-- p-mode conceptual intro goes here. Points to cover:
- p-modes are acoustic standing waves; restoring force is pressure
- characterised by radial order n and angular degree ell; disk-integrated photometry gives us low-ell modes
- the oscillation cavity spans most of the stellar interior so p-mode frequencies average over the whole star
- the power spectrum shows a regular comb of peaks separated by Delta_nu; the small offset between ell=0 and ell=2 neighbours is delta_nu_02
- g-modes are buoyancy-driven and confined to the radiative interior — not observable at the surface for solar-type stars, so we do not use them here
- add a schematic or echelle diagram if one is available
-->

Stars oscillate in resonant modes. The frequencies of these modes encode information about the stellar interior — in particular, the sound speed profile and the Brunt–Väisälä frequency. By measuring just a few summary statistics of the oscillation spectrum, we can constrain both the global structure and the interior chemical composition of a star.

## Setup

Your work directory already contains a pre-compiled MESA model and a `run_star_extras.f90` with all the seismic quantities implemented. Copy it to your working area and compile:

```bash
cp -r Lab3 my_lab3
cd my_lab3
./clean && ./mk
```

The extra history columns written by `run_star_extras` will appear automatically in your history file. You do not need to modify `history_columns.list`.

> [!IMPORTANT]
> Do **not** edit `src/run_star_extras.f90`. All the physics is already there.

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

## Exercise 1: Run the model and watch it evolve

Run the model:

```bash
./rn
```

pgstar is enabled by default. As the model runs you will see the evolution of $\Delta\nu$, $\delta\nu_{02}$, and the HR diagram update in real time.

<!-- pgstar panel layout to finalise once inlist_pgstar is settled: HR diagram, Delta_nu vs age, delta_nu02 vs age, colour-magnitude. Confirm history column names used in History_Track panels. -->

Watch how $\Delta\nu$ and $\delta\nu_{02}$ evolve as the star ages. Which changes faster? What is happening physically in the stellar interior at that point?

---

## Exercise 2: Plot the seismic quantities with Python

Once the run has finished, use `mesa_reader` to load the history file and plot $\Delta\nu$ and $\delta\nu_{02}$ as a function of stellar age.

```python
import mesa_reader as mr
import matplotlib.pyplot as plt

h = mr.MesaData('LOGS/history.data')

fig, axes = plt.subplots(2, 1, sharex=True)

axes[0].plot(h.star_age / 1e9, h.Delta_nu_int * 1e6)
axes[0].set_ylabel(r'$\Delta\nu$ ($\mu$Hz)')

axes[1].plot(h.star_age / 1e9, h.delta_nu02_int * 1e6)
axes[1].set_ylabel(r'$\delta\nu_{02}$ ($\mu$Hz)')
axes[1].set_xlabel('Age (Gyr)')

plt.tight_layout()
plt.show()
```

<!-- also add colour vs magnitude panel here once Lab 2 column names are confirmed -->

How does the rate of change of each quantity compare to what you saw in the CMD from Lab 2?

---

## Exercise 3: Crowd-source the seismic grid

Enter your values at the age closest to 1 Gyr into the shared spreadsheet:

| Column | Value |
|--------|-------|
| `initial_mass` | your assigned mass ($M_\odot$) |
| $\Delta\nu$ at 1 Gyr | `Delta_nu_int` in $\mu$Hz |
| $\delta\nu_{02}$ at 1 Gyr | `delta_nu02_int` in $\mu$Hz |

Once everyone has contributed, look at the full grid. How well do $\Delta\nu$ and $\delta\nu_{02}$ separate stars of different masses at the same age? How does that compare to the CMD separation from Lab 2?

{{< details title="Discussion prompt" closed="true" >}}

$\Delta\nu$ scales with mean density, so it is sensitive to both mass and radius. $\delta\nu_{02}$ is sensitive to the core sound speed gradient, which depends on the central hydrogen abundance — almost independent of the envelope. How does adding $\delta\nu_{02}$ to an isochrone fit change the mass–age degeneracy? Which observable from Lab 2 does it complement most effectively?

{{< /details >}}
