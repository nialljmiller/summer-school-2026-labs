---
author: Niall Miller
---

# Lab 3: Asteroseismology

In Labs 1 and 2 we used rotation periods and CMD positions to constrain stellar ages. In this lab we add a third technique: asteroseismology. We will compute two seismic observables — the large frequency separation $\Delta\nu$ and the small frequency separation $\delta\nu_{02}$ — and ask whether they constrain ages better than what we measured in Lab 2. We will also build a map of which stars are actually accessible to real space missions.

All of the code that computes the seismic quantities is already implemented in `run_star_extras.f90`. You do not need to modify it. Your job is to configure the run, interpret the outputs, and combine results across the group.

---

## The large frequency separation $\Delta\nu$

The large frequency separation is the average spacing between p-modes of the same angular degree $\ell$ and consecutive radial order $n$. It is related to the sound crossing time of the star:

$$\Delta\nu = \left( 2 \int_0^R \frac{dr}{c_s} \right)^{-1}$$

```fortran
Delta_nu_int = 0.
do k = 2, s% nz, 1
   dr = s% rmid(k-1) - s% rmid(k)
   Delta_nu_int = Delta_nu_int + dr/s% csound(k)
end do
Delta_nu_int = 1./(2.*Delta_nu_int)
```

The loop accumulates $dr/c_s$ from the centre to the surface, building up the total sound travel time. The final line takes the reciprocal and divides by 2 to give $\Delta\nu$ in Hz. Because $\Delta\nu \propto \sqrt{M/R^3}$, it is a measure of the mean stellar density — it decreases as the star expands during main-sequence evolution.

---

## The small frequency separation $\delta\nu_{02}$

The small frequency separation is the offset between $\ell=0$ and $\ell=2$ modes of similar frequency. Unlike $\Delta\nu$, it is sensitive to the gradient of the sound speed in the stellar core, and therefore directly to the central hydrogen abundance:

$$\delta\nu_{02} \approx -\frac{2\,\Delta\nu}{\pi^2\,\nu_\mathrm{max}} \left( \int_0^R \frac{1}{r} \frac{dc_s}{dr} \, dr - \frac{c_s(R)}{R} \right)$$

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

The loop computes $dc_s/r$ shell by shell, approximating the radial derivative of the sound speed as a finite difference. The surface boundary term $c_s(R)/R$ is subtracted after the loop. The whole integral is then scaled by $-2\Delta\nu/(\pi^2\nu_\mathrm{max})$ to give $\delta\nu_{02}$ in Hz. Because the sound speed gradient in the core steepens as hydrogen is depleted, $\delta\nu_{02}$ decreases monotonically with age — making it a direct age clock.

---

## Step 1 — Setup

Lab 3 is a self-contained working directory. You do not need to copy anything from Lab 2. The `star` binary is already compiled, the five-panel pgstar display is already configured, and `run_star_extras.f90` already implements the seismic calculations.

The only thing you need to do before running is open `inlist_run` and set your assigned mass:

```fortran
initial_mass = X.X   ! set to your assigned value
```

Mass assignments for this lab are: **0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2** $M_\odot$ — one per team.

> [!NOTE]
> The `rn` script copies `inlist_run` to `inlist` before launching the model — always edit `inlist_run`, not `inlist` directly. `inlist` is overwritten every time you run.

---

## Step 2 — Configure the inlist

Look through the below pg star display. It shows five panels that update in real time:

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

  ! Panel 1 — HR diagram (top left)
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

  ! Panel 2 — Delta_nu vs age (top centre)
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
  History_Track1_ymin = 0
  History_Track1_ymax = 5e-4



  ! Panel 3 — delta_nu02 vs age (top right)
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
  History_Track2_ymin = 0
  History_Track2_ymax = 5e-5


  ! Panel 4 — Interior composition (bottom left)
  Grid1_plot_name(4) = 'Abundance'
  Grid1_plot_row(4) = 5
  Grid1_plot_rowspan(4) = 4
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

  Grid1_plot_name(5) = 'History_Track3'
  Grid1_plot_row(5) = 5
  Grid1_plot_rowspan(5) = 4
  Grid1_plot_col(5) = 6
  Grid1_plot_colspan(5) = 4
  Grid1_plot_pad_left(5) = 0.06
  Grid1_plot_pad_right(5) = 0.04
  Grid1_plot_pad_top(5) = 0.05
  Grid1_plot_pad_bot(5) = 0.06
  Grid1_txt_scale_factor(5) = 0.65
  History_Track3_title = 'Color-Color'
  History_Track3_xname = 'J'
  History_Track3_yname = 'H'
  History_Track3_xaxis_label = 'M_J'
  History_Track3_yaxis_label = 'M_H'

/ ! end of pgstar namelist
```

| Panel | Content |
|-------|---------|
| 1 (top-left) | HR diagram |
| 2 (top-centre) | $\Delta\nu$ vs age |
| 3 (top-right) | $\delta\nu_{02}$ vs age |
| 4 (bottom-left) | Interior composition (H, He, C, N vs mass coordinate) |
| 5 (bottom-right) | Text summary including `Delta_nu_int` and `delta_nu02_int` |

change the inlist to have this setup.

---

## Step 3 — Run the model

```bash
./rn
```

The model will run from the pre-main sequence to TAMS (central hydrogen fraction below 1%). Watch the five panels as the star evolves. Pay attention to:

- How quickly does $\Delta\nu$ change compared to the HR diagram position?
- How does $\delta\nu_{02}$ behave — does it change monotonically?
- What is happening to the interior composition at the same time?

> [!NOTE]
> Lower-mass stars take longer to reach TAMS. If your model is still running when you need to move to the Python exercises, you have enough history data to proceed.

> [!TIP]
> If the run is interrupted (e.g. by closing the terminal), restart it with `./re` rather than `./rn`. `./re` picks up from the most recent photo in `photos/` without restarting from the pre-main sequence. Do not use `./re` after editing `inlist_run` — use `./rn` instead so the updated inlist is copied through.

---

## Step 4 — Age diagnostics: seismic vs CMD

Before pooling results, compare what each observable is actually telling you.

**$\Delta\nu$ as an age clock.** Because $\Delta\nu \propto \sqrt{\bar{\rho}} \propto \sqrt{M/R^3}$, it tracks the mean density. On the main sequence the star expands slowly, so $\Delta\nu$ decreases — but the rate depends on mass. This makes $\Delta\nu$ a useful density indicator but a coarse age clock on its own, because you need to know $M$ independently to convert density to age.

**$\delta\nu_{02}$ as an age clock.** The small separation probes the sound speed gradient in the core. As hydrogen burns, the mean molecular weight of the core increases, the sound speed drops, and the gradient steepens — so $\delta\nu_{02}$ decreases monotonically with central hydrogen abundance. This makes $\delta\nu_{02}$ a nearly mass-independent age indicator along the main sequence at fixed $\Delta\nu$. This is the basis of the Christensen-Dalsgaard (C–D) diagram: plotting $\delta\nu_{02}$ versus $\Delta\nu$ produces a grid where lines of constant mass and constant age cross at different angles, allowing both to be read off simultaneously from two observables.

**Comparison with Lab 2.** In Lab 2 you measured stellar ages from CMD position — the $J-K_s$ colour and absolute $K_s$ magnitude. That method has two limitations: (1) at young ages the main sequence is nearly vertical in colour-magnitude space, giving poor age resolution; and (2) photometric uncertainty propagates directly into age uncertainty through the isochrone width. Seismic observables sidestep both. The $\delta\nu_{02}$ versus $\Delta\nu$ diagram separates models that are photometrically almost indistinguishable on the CMD, and the observables are distance-independent.

---

## Step 5 — Python: reproduce the history plots

Once the run has enough history data, use `mesa_reader` to reproduce the four key plots. The `python_helpers/` directory contains more complete plotting scripts — the code below is a minimal example you can run directly.

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

# CMD — 2MASS Vega-system columns written by the colors module
# Verify the exact column names with: head -7 LOGS/history.data
axes[0, 1].plot(h.J - h.Ks, h.Ks)
axes[0, 1].invert_yaxis()
axes[0, 1].set_xlabel(r'$J - K_s$')
axes[0, 1].set_ylabel(r'$K_s$ (mag)')
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
> `Delta_nu_int` and `delta_nu02_int` come from `run_star_extras.f90` as extra history columns — they are NOT in `history_columns.list` and you do not need to add them there. If you get an `AttributeError` on `h.J` or `h.Ks`, run `head -7 LOGS/history.data` and check the sixth line for the exact column names written by the colors module.

---

## Step 6 — Crowd-source the seismic grid

Find the timestep in your history closest to 1 Gyr and enter the following into the shared spreadsheet:

| Column | Where to find it |
|--------|-----------------|
| `initial_mass` ($M_\odot$) | your assigned value |
| $L/L_\odot$ | `luminosity` at 1 Gyr |
| $T_\mathrm{eff}$ (K) | `Teff` at 1 Gyr |
| 2MASS colour ($J-K_s$) | `J - Ks` at 1 Gyr |
| $\Delta\nu$ ($\mu$Hz) | `Delta_nu_int` $\times 10^6$ at 1 Gyr |
| $\delta\nu_{02}$ ($\mu$Hz) | `delta_nu02_int` $\times 10^6$ at 1 Gyr |

```python
import mesa_reader as mr
import numpy as np

h = mr.MesaData('LOGS/history.data')

target_ages = [1e9, 3e9, 5e9, 7e9, 9e9]

for age in target_ages:
    idx = np.argmin(np.abs(h.star_age - age))
    if h.star_age[idx] < age * 1.2:
        print(f"Age: {h.star_age[idx]/1e9:.1f} Gyr  "
              f"Teff: {h.Teff[idx]:.0f} K  "
              f"L/Lsun: {10**h.log_L[idx]:.4f}  "
              f"J-Ks: {h.J[idx] - h.Ks[idx]:.4f}  "
              f"Delta_nu: {h.Delta_nu_int[idx]*1e6:.2f} uHz  "
              f"delta_nu02: {h.delta_nu02_int[idx]*1e6:.2f} uHz")
```

...and to check the mass from the inlist

```python 
with open('inlist') as f:
    for line in f:
        if 'initial_mass' in line:
            print(line.strip())
```

Once the full group has contributed, look at the complete grid. How well do $\Delta\nu$ and $\delta\nu_{02}$ separate stars of different masses at the same age? How does this compare to the CMD separation from Lab 2?

---

## Step 7 — Telescope detection map

Using the crowd-sourced grid, estimate whether each star's oscillations are detectable by current instruments.

### Amplitude scaling relations

**Radial-velocity amplitude** (Campante et al. 2024, Fig. 5):

$$A_\mathrm{RV} = 18 \left(\frac{L/L_\odot}{M/M_\odot}\right)^{1.5} \quad \text{cm s}^{-1}$$

The solar value is $A_{\mathrm{RV},\odot} = 18$ cm s$^{-1}$. The exponent steepens to 1.5 below $\sim 5500$ K — the regime relevant for K dwarfs.

**Photometric amplitude** (Sayeed et al. 2025, Eq. 10; activity term set to solar):

$$A_\mathrm{phot} = 3.6 \,\frac{(L/L_\odot)^{0.691}}{(M/M_\odot)^{1.013}\,(T_\mathrm{eff}/5777)^{2.053}\,(T_\mathrm{eff}/5934)^{0.8}} \quad \text{ppm}$$

The solar value is $A_{\mathrm{phot},\odot} = 3.6$ ppm.

### Instrument noise floors

| Instrument | Noise floor | Notes |
|-----------|-------------|-------|
| EPRV (ESPRESSO, KPF) | 30 cm s$^{-1}$ per $\sqrt{\text{min}}$ | State-of-the-art RV |
| Kepler / TESS (best case) | 240 ppm per $\sqrt{\text{min}}$ | From 12 ppm per 6.5 hr: $12\sqrt{390} \approx 237$ ppm |

### Observation time for SNR = 3

Noise scales as $1/\sqrt{N}$ where $N$ is the number of 1-minute exposures. Setting SNR = 3:

$$t \;[\text{days}] = \frac{(3 \times A_\mathrm{noise,1min})^2}{A^2 \times 1440}$$

**Solar sanity check:**
- EPRV: $t = (3 \times 30 / 18)^2 / 1440 \approx 0.017$ d (25 min) ✓  
- Space photometry: $t = (3 \times 240 / 3.6)^2 / 1440 \approx 28$ d ✓

### Excel formulas for the shared spreadsheet

Add four calculated columns to the right of the quantities you entered in Step 6. The spreadsheet columns are:

| A | B | C | D | E | F |
|---|---|---|---|---|---|
| mass (M☉) | L (L☉) | Teff (K) | J−Ks | Δν (μHz) | δν₀₂ (μHz) |

**Column G — RV amplitude (cm s⁻¹)**
```
=18*(B2/A2)^1.5
```

**Column H — Photometric amplitude (ppm)**
```
=3.6*B2^0.691/(A2^1.013*(C2/5777)^2.053*(C2/5934)^0.8)
```

**Column I — Days needed for SNR = 3 with EPRV**
```
=(3*30/G2)^2/1440
```

**Column J — Days needed for SNR = 3 with space photometry**
```
=(3*240/H2)^2/1440
```

Copy rows 2 downward for each team's entry. Columns G–J update automatically.

### What to look for

Plot columns G and H as a function of mass. The key result: stars below roughly $0.6$–$0.7\,M_\odot$ require years of continuous observation for photometric detection, and even EPRV becomes marginal below $\sim 0.4\,M_\odot$, because amplitude falls steeply as $(L/M)^{1.5}$ along the lower main sequence.

{{< details title="Discussion questions" closed="true" >}}

- Which stars in your grid are detectable by Kepler/TESS within a typical mission lifetime (~4 yr)? By an EPRV campaign of ~1 month?
- What does this tell you about the bias in real asteroseismic catalogues — specifically, why the C–D diagram has been populated only for FGK dwarfs so far?
- For stars where asteroseismology is detectable, how does the age precision from $\Delta\nu + \delta\nu_{02}$ compare to what you could infer from the CMD in Lab 2?
- What would PLATO add?

{{< /details >}}
