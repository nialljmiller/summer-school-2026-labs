---
author: Niall Miller
---

# Lab 3: Asteroseismology

In Labs 1 and 2 we used rotation periods and CMD positions to constrain stellar ages. In this lab we add a third technique: asteroseismology. We will compute two seismic observables — the large frequency separation $\Delta\nu$ and the small frequency separation $\delta\nu_{02}$ — and ask whether they constrain ages better than what we measured in Lab 2. We will also build a map of which stars are actually accessible to real space missions.

All of the code that computes the seismic quantities is already implemented in `run_star_extras.f90`. You do not need to modify it. Your job is to configure the run, interpret the outputs, and combine results across the group.



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

## Step 1 — Setup

This lab builds directly on Lab 2. Copy your Lab 2 working directory into a new one and replace `run_star_extras.f90` with the Lab 3 version:

```bash
cp -r lab2 lab3
cd lab3
cp /path/to/Lab3/src/run_star_extras.f90 src/
./clean && ./mk
```

> [!IMPORTANT]
> Do **not** edit `src/run_star_extras.f90`. All the seismic physics is already implemented there.

Your colors module setup from Lab 2 carries over automatically — no changes needed to the `&colors` namelist.

Now open `inlist_run` and set your assigned mass:

```fortran
initial_mass = X.X   ! set to your assigned value
```

Mass assignments for this lab are: **0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2** $M_\odot$ — one per team.

---

## Step 2 — Enable pgstar

The inlist currently has `pgstar_flag = .false.`. Change it to `.true.` in `&star_job`:

```fortran
pgstar_flag = .true.
```

Then replace the empty `&pgstar` block in `inlist_run` with the following. This sets up a six-panel display that updates in real time as the model evolves:

```fortran
&pgstar

  file_white_on_black_flag = .false.

  Grid1_win_flag = .true.
  Grid1_win_width = 18
  Grid1_win_aspect_ratio = 0.56

  Grid1_file_flag = .true.
  Grid1_file_dir = 'pgplot'
  Grid1_file_prefix = 'grid_'
  Grid1_file_interval = 10
  Grid1_file_width = 18

  Grid1_num_cols = 9
  Grid1_num_rows = 9
  Grid1_num_plots = 5

  Grid1_xleft = 0.01
  Grid1_xright = 0.99
  Grid1_ybot = 0.02
  Grid1_ytop = 0.98

  ! Panel 1 — HR diagram
  Grid1_plot_name(1) = 'HR'
  Grid1_plot_row(1) = 1
  Grid1_plot_rowspan(1) = 4
  Grid1_plot_col(1) = 1
  Grid1_plot_colspan(1) = 3
  Grid1_plot_pad_left(1) = 0.06
  Grid1_plot_pad_right(1) = 0.02
  Grid1_plot_pad_top(1) = 0.05
  Grid1_plot_pad_bot(1) = 0.06
  Grid1_txt_scale_factor(1) = 0.65
  HR_title = 'HR diagram'

  ! Panel 2 — Delta_nu vs age
  Grid1_plot_name(2) = 'History_Track1'
  Grid1_plot_row(2) = 1
  Grid1_plot_rowspan(2) = 4
  Grid1_plot_col(2) = 4
  Grid1_plot_colspan(2) = 3
  Grid1_plot_pad_left(2) = 0.06
  Grid1_plot_pad_right(2) = 0.02
  Grid1_plot_pad_top(2) = 0.05
  Grid1_plot_pad_bot(2) = 0.06
  Grid1_txt_scale_factor(2) = 0.65
  History_Track1_title = 'Large frequency separation'
  History_Track1_xname = 'star_age'
  History_Track1_yname = 'Delta_nu_int'
  History_Track1_xaxis_label = 'Age (yr)'
  History_Track1_yaxis_label = 'Delta-nu (Hz)'
  History_Track1_xmin = 0
  History_Track1_xmax = 1.3e10
  History_Track1_ymin = -1
  History_Track1_ymax = -1

  ! Panel 3 — delta_nu02 vs age
  Grid1_plot_name(3) = 'History_Track2'
  Grid1_plot_row(3) = 1
  Grid1_plot_rowspan(3) = 4
  Grid1_plot_col(3) = 7
  Grid1_plot_colspan(3) = 3
  Grid1_plot_pad_left(3) = 0.06
  Grid1_plot_pad_right(3) = 0.04
  Grid1_plot_pad_top(3) = 0.05
  Grid1_plot_pad_bot(3) = 0.06
  Grid1_txt_scale_factor(3) = 0.65
  History_Track2_title = 'Small frequency separation'
  History_Track2_xname = 'star_age'
  History_Track2_yname = 'delta_nu02_int'
  History_Track2_xaxis_label = 'Age (yr)'
  History_Track2_yaxis_label = 'delta-nu02 (Hz)'
  History_Track2_xmin = 0
  History_Track2_xmax = 1.3e10
  History_Track2_ymin = -1
  History_Track2_ymax = -1

  ! Panel 4 — Abundance profile
  Grid1_plot_name(4) = 'Abundance'
  Grid1_plot_row(4) = 5
  Grid1_plot_rowspan(4) = 3
  Grid1_plot_col(4) = 1
  Grid1_plot_colspan(4) = 5
  Grid1_plot_pad_left(4) = 0.06
  Grid1_plot_pad_right(4) = 0.02
  Grid1_plot_pad_top(4) = 0.04
  Grid1_plot_pad_bot(4) = 0.06
  Grid1_txt_scale_factor(4) = 0.65
  Abundance_title = 'Interior composition'
  Abundance_num_isos_to_show = 4
  Abundance_which_isos_to_show(1) = 'h1'
  Abundance_which_isos_to_show(2) = 'he4'
  Abundance_which_isos_to_show(3) = 'c12'
  Abundance_which_isos_to_show(4) = 'n14'
  Abundance_xaxis_name = 'mass'
  Abundance_log_mass_frac_min = -4.0

  ! Panel 5 — Text summary
  Grid1_plot_name(5) = 'Text_Summary1'
  Grid1_plot_row(5) = 5
  Grid1_plot_rowspan(5) = 3
  Grid1_plot_col(5) = 6
  Grid1_plot_colspan(5) = 4
  Grid1_plot_pad_left(5) = 0.02
  Grid1_plot_pad_right(5) = 0.02
  Grid1_plot_pad_top(5) = 0.02
  Grid1_plot_pad_bot(5) = 0.02
  Grid1_txt_scale_factor(5) = 0.65

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
> The CMD panel uses `History_Track1_xname = 'J_MINUS_K'` as a placeholder. pgstar cannot subtract two columns directly — you will need to confirm the exact pre-computed color column name that your run writes to `history.data` and update this accordingly.

---

## Step 3 — Run the model

```bash
./rn
```

The model will run to the end of the main sequence (TAMS). Watch the six panels as the star evolves. Pay attention to:

- How quickly does $\Delta\nu$ change compared to the CMD position?
- How does $\delta\nu_{02}$ behave — does it change monotonically?
- What is happening to the interior composition (Panel 5) at the same time?

> [!NOTE]
> Lower mass stars take longer to reach TAMS. If your model is still running when you need to move to the Python exercises, you have enough history data to proceed — you do not need to wait for the run to finish.

---

## Step 4 — Age diagnostics

<!-- TO BE DEFINED AT NEXT MEETING -->
<!-- Decision pending on whether this is a qualitative discussion or a quantitative Python exercise comparing age uncertainty from CMD vs seismic observables. -->

> [!NOTE]
> Content for this section will be added after the next team meeting.

---

## Step 5 — Python: reproduce the history plots

Once the run has enough history data, use `mesa_reader` to reproduce the four key plots from the pgstar display.

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
axes[0, 1].plot(h.J - h.Ks, h.Ks)
axes[0, 1].invert_yaxis()
axes[0, 1].set_xlabel(r'$J - K$')
axes[0, 1].set_ylabel(r'$K$ (mag)')
axes[0, 1].set_title('CMD')

# Delta_nu vs age
axes[1, 0].plot(h.star_age / 1e9, h.Delta_nu_int * 1e6)
axes[1, 0].set_xlabel('Age (Gyr)')
axes[1, 0].set_ylabel(r'$\Delta\nu$ ($\mu$Hz)')
axes[1, 0].set_title('Large frequency separation')

# delta_nu02 vs age
axes[1, 1].plot(h.star_age / 1e9, h.delta_nu02_int * 1e6)
axes[1, 1].set_xlabel('Age (Gyr)')
axes[1, 1].set_ylabel(r'$\delta\nu_{02}$ ($\mu$Hz)')
axes[1, 1].set_title('Small frequency separation')

plt.tight_layout()
plt.savefig('lab3_history.png', dpi=150)
plt.show()
```

> [!TIP]
> The column names `h.J` and `h.K` are placeholders. The exact names depend on what your 2MASS filter files are called. Check `LOGS/history.data` header line to confirm.

---

## Step 6 — Crowd-source the seismic grid

Find the timestep in your history closest to 1 Gyr and enter the following into the shared spreadsheet:

| Column | Where to find it |
|--------|-----------------|
| `initial_mass` ($M_\odot$) | your assigned value |
| Spectral type | from `Teff` at 1 Gyr |
| 2MASS color ($J-K$) | from your history file |
| $\Delta\nu$ ($\mu$Hz) | `Delta_nu_int` × 10⁶ |
| $\delta\nu_{02}$ ($\mu$Hz) | `delta_nu02_int` × 10⁶ |

```python
import mesa_reader as mr
import numpy as np

h = mr.MesaData('LOGS/history.data')

# Find index closest to 1 Gyr
idx = np.argmin(np.abs(h.star_age - 1e9))

print(f"Age:          {h.star_age[idx]/1e9:.3f} Gyr")
print(f"Teff:         {h.Teff[idx]:.0f} K")
print(f"J - K:        {h.J[idx] - h.Ks[idx]:.4f}")
print(f"Delta_nu:     {h.Delta_nu_int[idx]*1e6:.2f} uHz")
print(f"delta_nu02:   {h.delta_nu02_int[idx]*1e6:.2f} uHz")
```

Once the full group has contributed, look at the complete grid. How well do $\Delta\nu$ and $\delta\nu_{02}$ separate stars of different masses at the same age? How does this compare to the CMD separation you saw in Lab 2?

> [!NOTE]
> **Amplitude column — decision pending.** At the next team meeting we will decide whether oscillation amplitude (estimated from scaling relations and 2MASS bolometric corrections) is computed here and added to the spreadsheet, or computed later in Step 7 directly from the exported grid data.

---

## Step 7 — Telescope detection map

Using the exported crowd-source grid, plot oscillation amplitude as a function of mass and spectral type, and overplot the photometric noise floors for **Kepler**, **TESS**, and **PLATO**.

> [!NOTE]
> The amplitude calculation and telescope noise floor values will be added here once the decision from Step 6 is confirmed. The key result students should find: stars below roughly 0.6–0.7 $M_\odot$ fall below the detection threshold of current space missions, which introduces a systematic observational bias in the asteroseismic sample.

{{< details title="Discussion questions" closed="true" >}}

- Which stars in your grid are detectable by Kepler? By TESS? By PLATO?
- What does this tell you about the bias in real asteroseismic catalogues?
- For stars where asteroseismology is detectable, how does the age precision from $\Delta\nu + \delta\nu_{02}$ compare to what you could infer from the CMD in Lab 2?
- What could PLATO measure that Kepler and TESS could not?

{{< /details >}}
