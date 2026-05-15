*Authors: Niall Miller (lead TA), Eliza Frankel, Joey Mombarg - Lecturer: Yaguang Li — MESA Summer School 2026, Tetons, Wyoming*

# Lab 3: Asteroseismology

In Labs 1 and 2 we used rotation periods and CMD positions to constrain stellar ages. In this lab we add a third technique: asteroseismology. We will compute two seismic observables -- the large frequency separation $\Delta\nu$ and the small frequency separation $\delta\nu_{02}$ -- and ask whether they constrain ages better than what we measured in Lab 2. We will also build a map of which stars are actually accessible to real space missions.

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

> [!TIP]
> DO NOT TOUCH 'run_star_extras.f90'. All of the code that computes the seismic quantities is already implemented in `run_star_extras.f90`. You do not need to modify it. Your job is to configure the run, interpret the outputs, and combine results across the group.




The loop accumulates $dr/c_s$ from the centre to the surface, building up the total sound travel time. The final line takes the reciprocal and divides by 2 to give $\Delta\nu$ in Hz. Because $\Delta\nu \propto \sqrt{M/R^3}$, it is a measure of the mean stellar density -- it decreases as the star expands during main-sequence evolution.

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

The loop computes $dc_s/r$ shell by shell, approximating the radial derivative of the sound speed as a finite difference. The surface boundary term $c_s(R)/R$ is subtracted after the loop. The whole integral is then scaled by $-2\Delta\nu/(\pi^2\nu_\mathrm{max})$ to give $\delta\nu_{02}$ in Hz. Because the sound speed gradient in the core steepens as hydrogen is depleted, $\delta\nu_{02}$ decreases monotonically with age - making it a direct age clock. (https://ui.adsabs.harvard.edu/abs/2005MNRAS.356..671O/abstract)

---

## Step 1 -- Setup

Lab 3 is a self-contained working directory. You do not need to copy anything from Lab 2 but you can copy your inlist from Lab 2 into this lab for a faster setup.

A good check when inheriting a MESA working directory could go something like:
```bash
cd content/monday/Lab3       #move to where the directory is

ls                          #to see whats here

cat rn                      #lets see what rn does (it copies 'inlist_run' to 'inlist' then calls rn1)

cat rn1                     #it calls star with 'inlist'
                            #so we have established that 'inlist_run' is what we modify. It will be copied to inlist for execution.
cp ../Lab2/whatever_inlist_was_called inlist_run
```

`run_star_extras.f90` already implements the seismic calculations.


We are going to be crowd sourcing our science. 

This link : [google sheets](https://docs.google.com/spreadsheets/d/1C88C5V2siCAaK8-3qgAZoNc9-9IH-RTIqFVetXQc3EM/edit?usp=sharing)

Will take you to a google sheets where you can add your name to the A column and choose as mass from 

**0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2**

You can do multiple masses but try not to do repeats of the same few. 

Once you have chosen a mass, you need to open `inlist_run` and set your assigned mass:

```fortran
initial_mass = X.X   ! set to your assigned value
```

Mass assignments for this lab are: **0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2** $M_\odot$.

> [!NOTE]
> The `rn` script copies `inlist_run` to `inlist` before launching the model -- always edit `inlist_run`, not `inlist` directly. `inlist` is overwritten every time you run.

---

## Step 2 -- Configure the inlist

Look through the below pg star display. It configures pg star to show five panels that update in real time:

Before you run anything, you need to set up two things: the pgstar live display and the colors module. 

pgstar is MESA's built-in plotting system -- it opens a window that updates every few timesteps so you can watch the star evolve in real time. 
The colors module computes synthetic photometry from the model atmosphere at each step, which is what gives you the 2MASS magnitudes in the history file.

The namelist below configures five panels:

- HR diagram
- $\Delta\nu$ vs age
- $\delta\nu_{02}$ vs age
- Interior abundance profile
- 2MASS magnitude track

Copy this into the &pgstar section of your inlist_run:

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

  ! Panel 1 -- HR diagram (top left)
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

  ! Panel 2 -- Delta_nu vs age (top centre)
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



  ! Panel 3 -- delta_nu02 vs age (top right)
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


  ! Panel 4 -- Interior composition (bottom left)
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

Now for the colors module. 
The namelist below is what you need, but we should check the paths before you paste it in. 
The stellar_atm and vega_sed paths point to files that need to exist on your machine. 
Run ```ls``` on those paths to confirm they resolve before continuing. 
Once you have verified them, add this &colors namelist to your inlist_run:

```fortran
&colors

   use_colors = .true.                                              !use colors

   instrument = '../data/filters/2MASS/2MASS'                       !We are assuming that data is in this dir
   stellar_atm = '../data/stellar_models/Kurucz2003all__alpha_00'   !check to see where it actually is. 

   vega_sed = '../data/stellar_models/vega_flam.csv'                !same as above
   mag_system = 'Vega'                                              

   distance = 3.0857d19                                             !10 parsecs -> absolute magnitudes

   make_csv = .true.                                                !Make a csv for each filter
   colors_results_directory = 'SED'                                 !put them in the SED/ directory
   sed_per_model = .false.                                          !overwrite them at every step

/ ! end of colors namelist
```
The distance = 3.0857d19 sets the distance to 10 parsecs, which is what puts the magnitudes on the absolute scale you used in Lab 2.

---

## Step 3 -- Run the model

```bash
./rn
```

The model will run from the pre-main sequence to **T**erminal **A**ge **M**ain **S**equence. Pay attention to:

- How quickly does $\Delta\nu$ change compared to the HR diagram position?
- How does $\delta\nu_{02}$ behave -- does it change monotonically?
- What is happening to the interior composition at the same time?

> [!NOTE]
> Lower-mass stars take longer to reach TAMS.

> [!TIP]
> If the run is interrupted (e.g. by closing the terminal), restart it with `./re` rather than `./rn`. `./re` picks up from the most recent photo in `photos/` without restarting from the pre-main sequence. Do not use `./re` after editing `inlist_run` -- use `./rn` instead so the updated inlist is copied through.

---

## Step 4 -- Let's think about this as the model runs -- Age diagnostics: seismic vs CMD


As the model runs, compare what each observable is actually telling you.

**$\Delta\nu$ as an age clock.** Because $\Delta\nu \propto \sqrt{\bar{\rho}} \propto \sqrt{M/R^3}$, it tracks the mean density. On the main sequence the star expands slowly, so $\Delta\nu$ decreases -- but the rate depends on mass. This makes $\Delta\nu$ a useful density indicator but a coarse age clock on its own, because you need to know $M$ independently to convert density to age.

**$\delta\nu_{02}$ as an age clock.** The small separation probes the sound speed gradient in the core. As hydrogen burns, the mean molecular weight of the core increases, the sound speed drops, and the gradient steepens -- so $\delta\nu_{02}$ decreases monotonically with central hydrogen abundance. This makes $\delta\nu_{02}$ a nearly mass-independent age indicator along the main sequence at fixed $\Delta\nu$. This is the basis of the Christensen-Dalsgaard (C–D) diagram: plotting $\delta\nu_{02}$ versus $\Delta\nu$ produces a grid where lines of constant mass and constant age cross at different angles, allowing both to be read off simultaneously from two observables.

**Comparison with Lab 2.** In Lab 2 you measured stellar ages from CMD position -- the $J-K_s$ colour and absolute $K_s$ magnitude. 

{{< details title="What limitations could the method shown in Lab2 have?" closed="true" >}}
**(1)** at young ages the main sequence is nearly vertical in colour-magnitude space, giving poor age resolution
**(2)** photometric uncertainty propagates directly into age uncertainty through the isochrone width. 

Seismic observables sidestep both. The $\delta\nu_{02}$ versus $\Delta\nu$ diagram separates models that are photometrically almost indistinguishable on the CMD, and the observables are distance-independent.
{{< /details >}}

---

## Step 5 -- Plotting beyond pgstar with Python

Once the run has enough history data, use `mesa_reader` to reproduce the four key plots. The `python_helpers/` directory contains more complete plotting scripts -- the code below is a minimal example you can run directly.


{{< details title="Python tips" closed="true" >}}
You have options with how to run these. The easiest route is to probably launch python in live mode in your working directory.
```bash
cd content/monday/Lab3       #move to where the working directory is

python                      #open python live and paste the below code into it one by one to see the plots.
```
You can also save this code into a script and run the script. 
{{< /details >}}


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

# CMD -- 2MASS Vega-system columns written by the colors module
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
> `Delta_nu_int` and `delta_nu02_int` come from `run_star_extras.f90` as extra history columns -- they are NOT in `history_columns.list` and you do not need to add them there. If you get an `AttributeError` on `h.J` or `h.Ks`, run `head -7 LOGS/history.data` and check the sixth line for the exact column names written by the colors module.


**3D CMD** -- lift the CMD into 3D using $\Delta\nu$ as the third axis, revealing how the seismic observable evolves along the sequence:

```python
import mesa_reader as mr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

h = mr.MesaData('LOGS/history.data')

fig = plt.figure(figsize=(9, 7))
ax  = fig.add_subplot(111, projection='3d')

ax.scatter(h.J - h.Ks, h.Ks, h.Delta_nu_int * 1e6, s=5, c=h.star_age, cmap='plasma')

ax.set_xlabel(r'$J - K_s$')
ax.set_ylabel(r'$K_s$')
ax.set_zlabel(r'$\Delta\nu$ ($\mu$Hz)')
ax.invert_yaxis()
plt.tight_layout()
plt.show()
```


> [!TIP]
> Why don't you modify this plotting code? What happens when we plot the CMD with `delta_nu02_int` as the Z-axis?


---

## Step 6 -- Crowd-source the seismic grid

Run the snippet below to extract values at 1, 3, 5, 7, and 9 Gyr, then enter one row per age into the shared spreadsheet:

| Column | Where to find it |
|--------|-----------------|
| `initial_mass` ($M_\odot$) | your assigned value |
| Age (Gyr) | the target age for that row |
| $L/L_\odot$ | `10**log_L` at that age |
| $T_\mathrm{eff}$ (K) | `Teff` at that age |
| 2MASS colour ($J-K_s$) | `J - Ks` at that age |
| $\Delta\nu$ ($\mu$Hz) | `Delta_nu_int` $\times 10^6$ at that age |
| $\delta\nu_{02}$ ($\mu$Hz) | `delta_nu02_int` $\times 10^6$ at that age |

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

> [!NOTE]
> If your model hasn't reached a target age yet, the row for that age will be missing from the output -- this is fine, enter what you have.

...and to check the mass from the inlist (we can check the inlist as we earlier showed that the ```inlist``` file is made as a copy of the ```inlist_run``` file when we do ```./rn```)

```python 
with open('inlist') as f:
    for line in f:
        if 'initial_mass' in line:
            print(line.strip())
```

Add the results to the google sheets file.
This link : [google sheets](https://docs.google.com/spreadsheets/d/1C88C5V2siCAaK8-3qgAZoNc9-9IH-RTIqFVetXQc3EM/edit?usp=sharing)

> [!NOTE]
>The RV and photometric amplitude columns in the spreadsheet are calculated automatically from your MESA output. You do not need to fill them in.
>**RV amplitude** scales as $A_\mathrm{RV} \propto (L/M)^{1.5}$ (normalised to the solar value of 18 cm/s), following the scaling relation in [Chaplin et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024A%26A...683L..16C/abstract). The noise floor for state-of-the-art radial velocity spectrographs (EPRVs) is approximately 30 cm/s per minute of cadence, from [Beard et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv251101954B/abstract).
>**Photometric amplitude** follows the scaling relation in Equation 10 of [Schofield et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AJ....170..212S/abstract) (the activity term is set to solar). The noise floor for space photometry (Kepler/TESS) is approximately 240 ppm per minute of cadence -- derived from the best-case 12 ppm per 6.5 hr reported in [Gilliland et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJS..197....6G/abstract).

Once the full group has contributed, look at the complete grid. How well do $\Delta\nu$ and $\delta\nu_{02}$ separate stars of different masses at the same age? How does this compare to the CMD separation from Lab 2?


## Step 7 -- Who cares if we cant even measure it?

### Get your own copy of the notebook

The analysis lives in a shared Google Colab notebook. Open the link below, then **make a personal copy before you do anything else** -- this is important so your edits don't overwrite someone else's work and vice versa.

**→ [Open the Lab 3 detection map notebook](https://colab.research.google.com/drive/1cre1fH0yrvhCE0ZWSka4A4CloBuMWWwu#scrollTo=HkRl1EuNCjsr)**

https://colab.research.google.com/drive/1cre1fH0yrvhCE0ZWSka4A4CloBuMWWwu#scrollTo=HkRl1EuNCjsr


Once it opens:

```
File → Save a copy in Drive
```

or click the **"Copy to Drive"** button in the toolbar at the top of the page. Either way, a personal copy lands in your Google Drive and the original is untouched. Work from your copy for the rest of the lab.

### Download the data

In the shared Google Sheet:

```
File → Download → Comma Separated Values (.csv)
```

Save the file -- you'll upload it into the notebook using the file upload button in the first cell.

### What the notebook does

The notebook reads the crowd-sourced grid and produces two detection map plots (photometric and RV), coloured by stellar mass with age annotated on each point. Horizontal threshold lines show the minimum detectable amplitude for each instrument and campaign length. Stars above a line are detectable; stars below are not.

The parameters block near the top of the notebook lets you change instrument names, noise floors ($\sigma_0$), and campaign lengths -- the threshold lines update when you rerun the cell.

### What to look for

The y-axis is logarithmic. The steep drop in amplitude towards red $J-K_s$ (K and M dwarfs) is immediately visible. What masses can we realistically use this technique for? (https://arxiv.org/abs/1103.0702)

{{< details title="Discussion questions" closed="true" >}}

- Which stars in your grid sit above the TESS 1-sector threshold? Which need the full Kepler mission?
- Which are accessible to ESPRESSO in a single night? Which need a multi-week campaign?
- What does this tell you about the systematic bias in real asteroseismic catalogues? (https://arxiv.org/abs/2403.16333v1)
- For stars where asteroseismology is detectable, how does the age precision from $\Delta\nu + \delta\nu_{02}$ compare to what you could infer from the CMD in Lab 2?
- One of the first solar-like oscillation detection in a K5 dwarf ($\epsilon$ Indi, $J-K_s \approx 0.75$) (https://arxiv.org/abs/2403.16333) required 6 consecutive half-nights with ESPRESSO on the VLT. Where does it land on your plot?
- What might make this a "best case scenario" study? What have we not considered? 
{{< /details >}}


## Bonus Step -- What actually matters for wobbly stars?

Lab 2 showed you that changing the atmospheric boundary condition (`atm_T_tau_relation`) and the mixing length parameter ($\alpha_\mathrm{MLT}$) visibly shifts a star's track on the HR diagram and CMD (https://arxiv.org/pdf/2303.09596). Here you will ask: *by how much do those same changes affect seismic observables?*

The answer matters because if $\Delta\nu$ and $\delta\nu_{02}$ are insensitive to surface physics, they give more trustworthy ages than CMD position -- the seismic signal comes from the core, not the atmosphere.
This does not mean they are completely immune to surface physics, we are investigating just how robust these observables are.

### Setup

Keep your Lab 3 mass. Run a small grid varying one parameter at a time, recording seismic and photometric values at the same target ages like before:

**Vary the atmospheric boundary condition** (fix $\alpha_\mathrm{MLT} = 1.8$):

```fortran
atm_T_tau_relation = 'Eddington'       ! reference
atm_T_tau_relation = 'solar_Hopf'
atm_T_tau_relation = 'Krishna_Swamy'
atm_T_tau_relation = 'Trampedach_solar'
```

**Vary the mixing length parameter** (fix `atm_T_tau_relation = 'Eddington'`):

```fortran
mixing_length_alpha = 1.5
mixing_length_alpha = 2.0
mixing_length_alpha = 2.5
mixing_length_alpha = 3.0
```

> [!CAUTION]
> Change `star_history_name` for every run so files do not overwrite each other -- e.g. `'1p0Msun_Eddington_alpha1p8.data'`.

### Extract the values

Use the same snippet as Step 6, just targeting your chosen ages:

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

### Enter results into the shared spreadsheet

Go to the **Lab 3 Bonus** tab in the [shared spreadsheet](https://docs.google.com/spreadsheets/d/1C88C5V2siCAaK8-3qgAZoNc9-9IH-RTIqFVetXQc3EM/edit?usp=sharing). Enter one row per run. Columns D and E are the parameters you changed; columns F–J are the MESA output. Columns K–M calculate automatically.

The four bar charts at the bottom of the sheet show how much Teff, $J-K_s$, $\Delta\nu$, and $\delta\nu_{02}$ vary across all models. Compare the vertical scale of the seismic panels to the photometric ones.

{{< details title="Discussion questions" closed="true" >}}

- Which observables change the most between boundary conditions? Which barely move?
- Which observables change the most with $\alpha_\mathrm{MLT}$? Does the direction make physical sense?
- If you were trying to measure stellar ages and could only observe one of Teff, $J-K_s$, or $\delta\nu_{02}$, which would you choose and why?
- What does the sensitivity (or lack of it) to surface physics tell you about *where* the age information in $\delta\nu_{02}$ comes from?

{{< /details >}}
