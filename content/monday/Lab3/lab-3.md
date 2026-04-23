---
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

### pgstar setup

The pgstar configuration lives in `inlist_run` inside the `&pgstar` namelist. It is already enabled (`pgstar_flag = .true.` in `&star_job`). The full namelist to paste in is below — it sets up a six-panel grid that updates in real time as the model evolves.

```fortran
&pgstar

  ! Global window settings
  file_white_on_black_flag = .false.

  Grid1_win_flag = .true.
  Grid1_win_width = 18
  Grid1_win_aspect_ratio = 0.56

  Grid1_file_flag = .true.
  Grid1_file_dir = 'pgplot'
  Grid1_file_prefix = 'grid_'
  Grid1_file_interval = 10
  Grid1_file_width = 18

  Grid1_num_cols = 12
  Grid1_num_rows = 9
  Grid1_num_plots = 6

  Grid1_xleft = 0.01
  Grid1_xright = 0.99
  Grid1_ybot = 0.02
  Grid1_ytop = 0.98

  ! Panel 1 — HR diagram (top-left)
  Grid1_plot_name(1) = 'HR'
  Grid1_plot_row(1) = 1
  Grid1_plot_rowspan(1) = 4
  Grid1_plot_col(1) = 1
  Grid1_plot_colspan(1) = 3
  Grid1_plot_pad_left(1) = 0.06
  Grid1_plot_pad_right(1) = 0.02
  Grid1_plot_pad_top(1) = 0.05
  Grid1_plot_pad_bot(1) = 0.06
  Grid1_txt_scale_factor(1) = 0.55
  HR_title = 'HR diagram'

  ! Panel 2 — Colour-magnitude diagram (top-centre-left)
  ! Update History_Track1_xname and yname to match your filter columns.
  Grid1_plot_name(2) = 'History_Track1'
  Grid1_plot_row(2) = 1
  Grid1_plot_rowspan(2) = 4
  Grid1_plot_col(2) = 4
  Grid1_plot_colspan(2) = 3
  Grid1_plot_pad_left(2) = 0.06
  Grid1_plot_pad_right(2) = 0.02
  Grid1_plot_pad_top(2) = 0.05
  Grid1_plot_pad_bot(2) = 0.06
  Grid1_txt_scale_factor(2) = 0.55
  History_Track1_title = 'CMD'
  History_Track1_xname = 'Johnson_B_minus_Johnson_V'
  History_Track1_yname = 'Johnson_V'
  History_Track1_xaxis_label = 'B - V'
  History_Track1_yaxis_label = 'V (mag)'
  History_Track1_reverse_y = .true.
  History_Track1_num_panels = 1

  ! Panel 3 — Delta_nu vs age (top-right)
  Grid1_plot_name(3) = 'History_Track2'
  Grid1_plot_row(3) = 1
  Grid1_plot_rowspan(3) = 4
  Grid1_plot_col(3) = 7
  Grid1_plot_colspan(3) = 3
  Grid1_plot_pad_left(3) = 0.06
  Grid1_plot_pad_right(3) = 0.02
  Grid1_plot_pad_top(3) = 0.05
  Grid1_plot_pad_bot(3) = 0.06
  Grid1_txt_scale_factor(3) = 0.55
  History_Track2_title = 'Large frequency separation'
  History_Track2_xname = 'star_age'
  History_Track2_yname = 'Delta_nu_int'
  History_Track2_xaxis_label = 'Age (yr)'
  History_Track2_yaxis_label = 'Delta-nu (Hz)'
  History_Track2_xmin = 0
  History_Track2_xmax = 1.3e10
  History_Track2_ymin = -1
  History_Track2_ymax = -1

  ! Panel 4 — delta_nu02 vs age (top far-right)
  Grid1_plot_name(4) = 'History_Track3'
  Grid1_plot_row(4) = 1
  Grid1_plot_rowspan(4) = 4
  Grid1_plot_col(4) = 10
  Grid1_plot_colspan(4) = 3
  Grid1_plot_pad_left(4) = 0.06
  Grid1_plot_pad_right(4) = 0.04
  Grid1_plot_pad_top(4) = 0.05
  Grid1_plot_pad_bot(4) = 0.06
  Grid1_txt_scale_factor(4) = 0.55
  History_Track3_title = 'Small frequency separation'
  History_Track3_xname = 'star_age'
  History_Track3_yname = 'delta_nu02_int'
  History_Track3_xaxis_label = 'Age (yr)'
  History_Track3_yaxis_label = 'delta-nu02 (Hz)'
  History_Track3_xmin = 0
  History_Track3_xmax = 1.3e10
  History_Track3_ymin = -1
  History_Track3_ymax = -1

  ! Panel 5 — Abundance profile (bottom-left, wide)
  Grid1_plot_name(5) = 'Abundance'
  Grid1_plot_row(5) = 5
  Grid1_plot_rowspan(5) = 3
  Grid1_plot_col(5) = 1
  Grid1_plot_colspan(5) = 6
  Grid1_plot_pad_left(5) = 0.06
  Grid1_plot_pad_right(5) = 0.02
  Grid1_plot_pad_top(5) = 0.04
  Grid1_plot_pad_bot(5) = 0.06
  Grid1_txt_scale_factor(5) = 0.55
  Abundance_title = 'Interior composition'
  Abundance_num_isos_to_show = 4
  Abundance_which_isos_to_show(1) = 'h1'
  Abundance_which_isos_to_show(2) = 'he4'
  Abundance_which_isos_to_show(3) = 'c12'
  Abundance_which_isos_to_show(4) = 'n14'
  Abundance_xaxis_name = 'mass'
  Abundance_log_mass_frac_min = -4.0

  ! Panel 6 — Text summary (bottom-right)
  Grid1_plot_name(6) = 'Text_Summary1'
  Grid1_plot_row(6) = 5
  Grid1_plot_rowspan(6) = 3
  Grid1_plot_col(6) = 7
  Grid1_plot_colspan(6) = 6
  Grid1_plot_pad_left(6) = 0.02
  Grid1_plot_pad_right(6) = 0.02
  Grid1_plot_pad_top(6) = 0.02
  Grid1_plot_pad_bot(6) = 0.02
  Grid1_txt_scale_factor(6) = 0.0

  Text_Summary1_name(1,1) = 'model_number'
  Text_Summary1_name(2,1) = 'star_age'
  Text_Summary1_name(3,1) = 'log_dt'
  Text_Summary1_name(4,1) = 'star_mass'
  Text_Summary1_name(5,1) = 'Teff'
  Text_Summary1_name(6,1) = 'luminosity'
  Text_Summary1_name(7,1) = 'radius'
  Text_Summary1_name(8,1) = 'log_g'

  Text_Summary1_name(1,2) = 'center_h1'
  Text_Summary1_name(2,2) = 'center_he4'
  Text_Summary1_name(3,2) = 'log_cntr_T'
  Text_Summary1_name(4,2) = 'log_cntr_Rho'
  Text_Summary1_name(5,2) = 'Delta_nu_int'
  Text_Summary1_name(6,2) = 'delta_nu02_int'
  Text_Summary1_name(7,2) = ''
  Text_Summary1_name(8,2) = ''

/ ! end of pgstar namelist
```

> [!TIP]
> The CMD panel uses `History_Track1_xname = 'Johnson_B_minus_Johnson_V'` as a placeholder. Once you have decided on your filter system, update this to the actual computed colour column name in your history file. pgstar cannot subtract two columns directly — you will need to add a computed colour column in `run_star_extras` or use a pre-computed difference column if the colors module outputs one.

### Colors module setup

The colors module computes synthetic photometry at each timestep by interpolating a stellar atmosphere grid and convolving with filter transmission curves. To enable it, add `use_colors_library = .true.` in `&star_job`, and add the following namelist to `inlist_run` between `&controls` and `&pgstar`:

```fortran
&colors

   ! Atmosphere grid: Kurucz2003all covers Teff 3500-50000 K, logg 0-5, [M/H] -2.5 to +0.5
   stellar_atm = 'Kurucz2003all'

   ! Filter system index file — lists which .dat filter files to load.
   ! Column names in history.data are the filter filenames with .dat stripped.
   ! e.g. Johnson_B.dat -> Johnson_B, Johnson_V.dat -> Johnson_V
   instrument = 'Johnson'   ! update to your chosen filter system

   ! Magnitude system
   mag_system = 'Vega'

   ! Distance in parsecs (10 pc gives absolute magnitudes)
   distance = 10d0

/ ! end of colors namelist
```

The magnitude columns will then appear automatically in `history.data` alongside the seismic columns.

### Run

```bash
./rn
```

Watch how $\Delta\nu$, $\delta\nu_{02}$, and the CMD position all evolve together. Which seismic quantity changes fastest? How does the rate of change compare to the star's motion in the CMD?

---

## Exercise 2: Plot the seismic quantities with Python

Once the run has finished, use `mesa_reader` to reproduce the seismic and photometric history plots from Exercise 1.

```python
import mesa_reader as mr
import matplotlib.pyplot as plt

h = mr.MesaData('LOGS/history.data')

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# HR diagram
axes[0, 0].plot(h.log_Teff, h.log_L)
axes[0, 0].invert_xaxis()
axes[0, 0].set_xlabel(r'$\log\,T_\mathrm{eff}$')
axes[0, 0].set_ylabel(r'$\log\,L/L_\odot$')
axes[0, 0].set_title('HR diagram')

# Colour-magnitude diagram
# column names depend on your chosen filter system — update as needed
axes[0, 1].plot(h.Johnson_B - h.Johnson_V, h.Johnson_V)
axes[0, 1].invert_yaxis()
axes[0, 1].set_xlabel(r'$B - V$')
axes[0, 1].set_ylabel(r'$V$')
axes[0, 1].set_title('CMD')

# Delta_nu vs age
axes[1, 0].plot(h.star_age / 1e9, h.Delta_nu_int * 1e6)
axes[1, 0].set_xlabel('Age (Gyr)')
axes[1, 0].set_ylabel(r'$\Delta\nu$ ($\mu$Hz)')

# delta_nu02 vs age
axes[1, 1].plot(h.star_age / 1e9, h.delta_nu02_int * 1e6)
axes[1, 1].set_xlabel('Age (Gyr)')
axes[1, 1].set_ylabel(r'$\delta\nu_{02}$ ($\mu$Hz)')

plt.tight_layout()
plt.show()
```

> [!NOTE]
> The colour column names (`Johnson_B`, `Johnson_V`, etc.) are set by the filter files listed in your colors instrument index. Check your `history_columns.list` to confirm the exact names output by your run.

---

## Exercise 3: Crowd-source the seismic grid

Enter your values at the age closest to 1 Gyr into the shared spreadsheet:

| Column | Value |
|--------|-------|
| `initial_mass` | your assigned mass ($M_\odot$) |
| $\Delta\nu$ at 1 Gyr | `Delta_nu_int` in $\mu$Hz |
| $\delta\nu_{02}$ at 1 Gyr | `delta_nu02_int` in $\mu$Hz |
| $B-V$ at 1 Gyr | from your colors columns |

Once everyone has contributed, look at the full grid. How well do $\Delta\nu$ and $\delta\nu_{02}$ separate stars of different masses at the same age? How does that compare to the CMD separation from Lab 2?

{{< details title="Discussion prompt" closed="true" >}}

$\Delta\nu$ scales with mean density, so it is sensitive to both mass and radius. $\delta\nu_{02}$ is sensitive to the core sound speed gradient, which depends on the central hydrogen abundance — almost independent of the envelope. How does adding $\delta\nu_{02}$ to an isochrone fit change the mass–age degeneracy? Which observable from Lab 2 does it complement most effectively?

{{< /details >}}
