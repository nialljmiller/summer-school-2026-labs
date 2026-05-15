+++
date = '2026-04-06T13:38:04+02:00'
draft = false
title = 'Lab 2 - Exploring MESA Custom Colors!'
+++

*Authors: Eliza Frankel (lead TA), Niall Miller, Joey Mombarg - Lecturer: Yaguang Li — MESA Summer School 2026, Tetons, Wyoming*

The MESA colors module allows us to generate synthetic photometry while running MESA stellar evolution models! It is a great way to merge observational and theoretical astronomy. With the colors module, we can specify what filter system and stellar atmosphere we want to use, and on top of regular MESA outputs (effective temperature, luminosity, age, etc.) we get bolometric magnitude, M$_{bol}$, bolometric flux, F$_{bol}$, and many synthetic magnitudes. For more information on the colors module, look at the [documentation](https://github.com/MESAHub/mesa/tree/main/colors).


One major age dating technique for stellar populations is through the use of isochrones. Isochrones are single-aged, chemically homogenous populations that show a snapshot of stellar evolution. They're made by evolving stars with the same chemical composition but different initial masses, and then finding what point in evolution each star is at at a particular age. Larger stars burn hotter and brighter, leaving the main sequence much quicker than a lower mass star. For example, at 10 Gyr we can see a 0.8 $M_{\odot}$ still on the main sequence, while a 5 $M_{\odot}$ star will be long past the Red Giant Branch. Because of this, we can build isochrones and use them to determine the age of stellar populations. One caveat to this is that they use the assumption that all the stars are at relatively the same distance and formed from the same materials at relatively the same time. _The best stellar populations to use isochrones when age dating stars is in clusters because we can make these assumptions._

This figure shows a series of isochrones at different ages between 0.03 Gyr to 10 Gyr, made using MIST. As the population gets older, the shape of the isochrone changes too!

<img width="450" height="600" alt="isochrones" src="https://github.com/user-attachments/assets/b115dfa9-6604-4531-9b04-d7dfb6481184" />



In this lab, we'll learn how atmospheric boundary conditions and the convective mixing length parameter can impact stellar evolution in both observational and theoretical coordinates. We'll also build isochrones to explore other techniques for age dating stellar populations and planet hosts.


### Step 1 - Directory Prep

Last lab we made a working directory that has everything we want to start lab 2. The first thing we will do is copy lab 1 into a new working directory:

```bash
cp -r lab1 lab2
cd lab2
```
Lets clean this directory and get rid of our outputs from Lab 1:

```bash
./clean
./mk
rm -r 1Msun_Z0p0134_Omega5000nHz_no_magnetic_braking ! Make this the name of your output directory defined in lab 1
```
In Lab 1, we explored magnetic braking. Let's make sure it's off for this lab in `&controls`:

```fortran
! Enable magnetic braking.
use_other_torque    = .false.
```

(A clean working directory can be downloaded [here]() under the name lab2)

### Step 2 - Building the inlist

For this lab, we are going to start with the same inlist as before, but we'll be adding a few things. Start by opening up `inlist_project` in a text editor.
We will be changing parameters in `&colors`, `&controls`, and the `inlist_pgstar`. Let's start with `&colors`!

#### Setting up Custom Colors in `&colors`

This is where we can enable synthetic photometry and determine what filters we'd like to use. Let's look through the [documentation](https://github.com/MESAHub/mesa/tree/main/colors).

The first thing we need to do is to make sure the colors module is on. By default, custom colors is turned off.

```fortran
use_colors = .true.
```

Next we need to decide what filter system, stellar atmosphere table, and Vega SED file to use. For this lab, we want to use the 2MASS filters and the Kurucz 2003 atmosphere tables. To find out what systems are available, let's move to data directory and start exploring! 

```bash
cd $MESA_DIR/data/colors_data/
```

Once you've found the right filters and atmosphere tables, add them to your inlist.

  {{< details title="Hint" closed="true" >}}

  ```fortran
  instrument = '/data/colors_data/ADD/FILTERS/HERE'
  stellar_atm = '/data/colors_data/ADD/MODELS/HERE/'
  vega_sed = '/data/colors_data/VEGA_SED/HERE/'
  ```
  {{< /details >}}


  {{< details title="Solution" closed="true" >}}

  ```fortran
  instrument = '/data/colors_data/filters/2MASS/2MASS'
  stellar_atm = '/data/colors_data/stellar_models/Kurucz2003all/'
  vega_sed = '/data/colors_data/stellar_models/vega_flam.csv'
  ```

  {{< /details >}}
  

> [!CAUTION]
> Proper syntax is important! Make sure that for the `instrument` directory there _isn't_ a '/' at the end, but for `stellar_atm` there _is_ a '/'

Now let's decide the distance of the star (in cm). For apparent magnitude, you can do any distance you want. For absolute magnitude, the distance should be 10 parsecs, or 3.0857 x 10<sup>19</sup> cm. Update the distance parameter to be 10 pc



{{< details title="Solution" closed="true" >}}

```fortran
  distance = 3.0857d19
```

{{< /details >}}

The last thing we need to do to make sure Custom Colors works in the `inlist` file.
**Question** Do you see anything pointing to `&colors`?


{{< details title="Answer" closed="true" >}}

No! Add the following lines of code to make sure MESA includes Custom Colors:

```fortran
&colors

   read_extra_colors_inlist(1) = .true.
   extra_colors_inlist_name(1) = 'inlist_project'

/ ! end of colors namelist

```

{{< /details >}}

What if you want to compare more than one filter system at the same time?

{{< details title="Bonus Task - More than one filter system" closed="true" >}}

Sometimes you want to use more than one filter system. To do this with Custom Colors, we must look into the data structure more. Follow these steps to make a joint Gaia-2MASS filter system you can use:

``` terminal
! in **one step above**your working directory
mkdir data
cd data
mkdir filters
cd filters
mkdir GAIA_2MASS
cd GAIA_2MASS
mkdir GAIA_2MASS
cd GAIA_2MASS

cp -r $MESA_DIR/data/colors_data/filters/GAIA/GAIA/*.dat .
cp -r $MESA_DIR/data/colors_data/filters/2MASS/2MASS/2MASS/*.dat .
```

Now that you've made this joint filter system, let's make a file called `GAIA_2MASS` and open it in your preferred text editor. Now add all the filters to use. For Gaia and 2MASS, it should look like:

```fortran
G.dat
Gbp.dat
Grp.dat
Grvs.dat
H.dat
J.dat
Ks.dat
```

Finally, you can replace the line for 'instrument' in `&colors` with

```fortran
instrument = '../data/filters/GAIA_2MASS/GAIA_2MASS'
```

** A completely working version of this can be downloaded [here](https://drive.google.com/drive/folders/1qebaN8Qt6e1nqiEHkt9A0T-jfyPIzXCE) called 'BONUS_data_mulitple_filters.zip'. Make sure this is in the directory **above** your working directory for it to work!

{{< /details >}}


#### `&controls`

This is the section with the main stellar evolution parameters. Our goal is to change the stellar input parameters to see how they change evolution! Keep the same stellar mass you used in Lab 1, but let's change the output directory to something simple:

`log_directory = 'LOGS'`

The first thing we want to change is how the atmospheric boundary conditions are controlled. Look through the _controls_ tab under star defaults in the [documentation](https://docs.mesastar.org/en/26.4.1/reference.html) to the right parameters to change. What does it control specifically?


{{< details title="Hint" closed="true" >}}

The section atmospheric boundary conditions has everything we'll need to start. `atm_option` is the main parameter we'll use. This changes how surface temperature (Tsurf) and surface pressure (Psurf) are evaluated at outer boundary conditions.

{{< /details >}}

Lets first use a **T($\tau$)** relationship. This defines how the atmospheric pressure structure is obtained by integrating the hydrostatic equilibrium equation,

$$
\frac{dP}{d\tau} = \frac{g}{\kappa}.
$$

Here, we assume that gravity, _g_, is spatially constant. There are 4 options for the **T($\tau$)** relationship: `Eddington`, `solar_Hopf`, `Krishna_Swamy`, and `Trampedach_solar`. Start by using the Eddington relationship.


  {{< details title="Solution" closed="true" >}}

  ```fortran
  atm_option = 'T_tau' 
  atm_T_tau_relation = 'Eddington' 
  atm_T_tau_opacity = 'varying' 
  ```

  {{< /details >}}

You can also change `atm_option` from 'T_tau', but be sure that you're using all the right parameters!
>[!caution]
> If you decide to try a different atmospheric option, carefully read the documentation. For example, if you wanted to try using `atm_option = 'table'`, make sure you comment out `atm_T_tau_relation` and `atm_T_tau_opacity` and add `atm_table` instead.

#### `inlist_pgstar`

For this lab, we are only going to use an HR diagram and a plot showing 2MASS magnitudes. You can either start with a blank inlist_pgstar (erase everything between `&pgstar` and `/ ! end of pgstar namelist` to copy the code below, or you can download the inlist_pgstar from [here](https://drive.google.com/drive/folders/1qebaN8Qt6e1nqiEHkt9A0T-jfyPIzXCE?usp=drive_link).


```fortran
&pgstar

   ! see star/defaults/pgstar.defaults

   ! MESA uses PGPLOT for live plotting and gives the user a tremendous
   ! amount of control of the presentation of the information.

   ! show HR diagram
   ! this plots the history of L,Teff over many timesteps
   HR_win_flag = .true.

   ! ! set static plot bounds
   HR_logT_min = 3.6
   HR_logT_max = 3.85
   HR_logL_min = -0.5
   HR_logL_max = 1

   ! set window size (aspect_ratio = height/width)
   HR_win_width = 6
   HR_win_aspect_ratio = 1.0


   ! Color Color diagram
   History_Track2_win_width = 6
   History_Track2_win_aspect_ratio = 1.0

   History_Track2_win_flag = .true.
   History_Track2_xname = 'J'
   History_Track2_yname = 'Ks'
   History_Track2_title = '2MASS Magnitudes'
   History_Track2_xaxis_label = 'J mag'
   History_Track2_yaxis_label = 'Ks mag'
   History_Track2_reverse_xaxis = .true.
   History_Track2_reverse_yaxis = .true.

   History_Track2_ymin = 1.75
   History_Track2_ymax = 3.5

   History_Track2_xmin = 2.5
   History_Track2_xmax = 4



/ ! end of pgstar namelist

```


### Step 3 - Changing parameters and running

#### Boundary Conditions

For this lab, we want to explore the different atmospheric boundary conditions and the mixing length parameter, $\alpha_{MLT}$. Start with just changing the boundary conditions.

In `&controls` above, we chose the Eddington T_tau relationship. Before we start running MESA, let's change one more parameter in `&controls` - because we want to compare how different parameters change evolution, we need to change the output file name so they don't overwrite each other. Make sure you give your new history file a descriptive name, for example if you are running a 1 $M_{\odot}$ star using the T_tau Eddington relationship, a good name would be: 

```fortran
star_history_name = '1p0Msun_TtauEddington_history.data'
```

Now you can `./rn` and watch the star evolve. 

Once it is done, try changing up the atmospheric boundary conditions and see what changes!


{{< details title="Hint" closed="true" >}}

There are many different combinations you can try! First, try changing `atm_T_tau_relation` between `solar_Hopf`, `Krishna_Swamy`, and `Trampedach_solar`. 

> [!CAUTION]
> Remember to change `star_history_name` to include the changes to atmospheric boundary conditons!

{{< /details >}}

Once you've explored how the atmospheric boundary conditions change evolution, set `atm_T_tau_relation` back to `Eddington`.

#### Mixing length parameter, $\alpha_{MLT}$

As we know, MESA is a 1 dimensional stellar evolution code which means it has to be creative when modeling 3D processes. In order to model energy transport in stars, MESA utilizes mixing length theory (MLT), which is the standard 1D parametarization of convection. One of the key parts of MLT is the mixing length parameter, $\alpha_{MLT}$. This is a unitless value that represents the convective efficiency of a region (i.e. what fraction of energy transport is being moved by convection rather than radiation). By changing $\alpha_{MLT}$, we can drastically change a star's main sequence lifetime, opacity, and more!

Look through the controls default parameters again and find the mixing length parameter, or $\alpha_{MLT}$. What value have we been using? 

{{< details title="Hint" closed="true" >}}
Check under the tab "mixing parameters" for the controls defaults
{{< /details >}}

{{< details title="Answer" closed="true" >}}
The default value is `mixing_length_alpha = 2.0d0`
{{< /details >}}

The solar mixing length parameter for the Eddington T($\tau$) atmospheric boundary condition is arounnd 1.80. Add the solar mixing length parameter to `&controls` and remember to change the name of your output history file so you know what the input parameters are. In this case, you could name your history file something like:

```fortran
star_history_name = '1p0Msun_alphaMLT1p80_history.data'
```

Each star has a different mixing length parameter, so you can't always use the solar value. For example, $\alpha$ Centauri A and B have mixing lengths that are 0.932 $\alpha_{MLT,\odot}$ and 1.095 $\alpha_{MLT,\odot}$, respectively ([Joyce & Chaboyer 2018b](https://ui.adsabs.harvard.edu/abs/2018ApJ...864...99J/abstract)). Try changing the value of `mixing_length_alpha` and running a model for both cases!



> [!TIP]
> You can't use math in an inlist, so if you wanted to have twice the mixing length, write `mixing_length_alpha = 3.6d0` rather than `mixing_length_alpha = 2 * 1.8d0`

Choose at least 2 more objects from the following table and run a model for each, changing the history file name so we can plot them soon.

| Object  | Type | Mixing Length  | Source |
| :---- | :-- |:---- |:-- |
| M92 | Metal poor globular cluster | 0.90 $\alpha_{MLT,\odot}$ | [Joyce & Chaboyer 2018a](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...10J/abstract) |
| HD 140283 | Subgiant | 0.88 $\alpha_{MLT,\odot}$ | [Joyce & Chaboyer 2018a](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...10J/abstract) |
| HIP 54639 | Main sequence | 0.28 $\alpha_{MLT,\odot}$ |[Joyce & Chaboyer 2018a](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...10J/abstract)  |
| HIP 106924| Main sequence | 0.52 $\alpha_{MLT,\odot}$ | [Joyce & Chaboyer 2018a](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...10J/abstract)  |
| KIC 1430163 | Star | 1.2 $\alpha_{MLT,\odot}$ | [Viani et al. 2018](https://iopscience.iop.org/article/10.3847/1538-4357/aab7eb/pdf) |
| KIC 1435467 | Star | 1.15 $\alpha_{MLT,\odot}$ | [Viani et al. 2018](https://iopscience.iop.org/article/10.3847/1538-4357/aab7eb/pdf) |


### Step 4 - Visualizing the changes outside of MESA

#### Isochrones

Go to the same [Google spreadsheet](https://docs.google.com/spreadsheets/d/1C88C5V2siCAaK8-3qgAZoNc9-9IH-RTIqFVetXQc3EM/edit?usp=sharing) as Lab 1. On the bottom, switch to the tab labeled "Lab 2". For this part, let's rerun a star with `mixing_length_alpha = 1.8d0` and the Eddington atmospheric boundary condition (`atm_T_tau_relation = 'Eddington'`). Once your star is done evolving, copy the values for "Teff" and "log(L)" from the terminal window into the Google sheet. _Make sure you are putting the values at the right corresponding age!_

As everyone finishes filling out the spreadsheet, we'll get to see an isochrone being built!

#### Comparing atmospheric boundary conditions and mixing length parameters

Now, go to the [Google Colab](https://colab.research.google.com/drive/1rFAu8UN0CC3GWllJfNyk7uV50FksOKok?usp=sharing) and make a copy of it.

Follow the instructions in the document to upload the different history files we made and visualize how changing the atmospheric boundary conditions and mixing length parameter can impact stellar evolution. 


### Inlist Solutions

{{< details title="Final inlist solutions!" closed="true" >}}

Here is what your inlist should look like! You can also download a copy from [here](https://drive.google.com/drive/folders/1qebaN8Qt6e1nqiEHkt9A0T-jfyPIzXCE?usp=drive_link) to make sure you get the lab working.

<details>
<summary>inlist_project</summary>

```fortran
&star_job

      pause_before_terminate = .true.
      show_log_description_at_start = .true.

      history_columns_file = 'custom_history_columns.list'
      profile_columns_file = 'custom_profile_columns.list'

      ! pgstar
      pgstar_flag = .true.

      ! pre main sequence
      create_pre_main_sequence_model = .true.
      pre_ms_T_c = 9.9d5 ! Initial central temperature.

      ! initial rotation
      new_omega =  3.1416d-5 ! 5000nHz
      set_near_zams_omega_steps = 15

      ! initial metal fractions
      initial_zfracs = 6 ! AGSS09_zfracs



/ ! end of star_job namelist

&eos

/ ! end of eos namelist

&kap

      !opacities with AGSS09 abundances
      kap_file_prefix = 'OP_a09_nans_removed_by_hand'
      kap_lowT_prefix = 'lowT_fa05_a09p'
      kap_CO_prefix   = 'a09_co'

      use_Type2_opacities = .false.


/ ! end of kap namelist

&controls

      ! ZAMS limit
      Lnuc_div_L_zams_limit = 0.95

      ! uniform viscosity

      ! initial mass
      initial_mass = 1d0

      ! initial He and Z
      initial_z = 0.0134
      initial_y = 0.2485

      ! stopping criterion
      xa_central_lower_limit_species(1) = 'h1'
      xa_central_lower_limit(1) = 0.01

      ! output
      log_directory = 'LOGS'
      history_interval = 1
      star_history_name = '1p0Msun_alphaMLT1p80_history.data'   !!!!!!!

      ! atmosphere options
      atm_option = 'T_tau'  !!!!!!!
      atm_T_tau_relation = 'Eddington'  !!!!!!!
      atm_T_tau_opacity = 'varying'  !!!!!!!

      ! Enable magnetic braking.
      use_other_torque    = .false.

      ! Mixing length parameter
      mixing_length_alpha = 1.80d0  !!!!!!!

/ ! end of controls namelist


&pgstar

! We set the pgstar controls in a seperate inlist instead.

/ ! end of pgstar namelist

&colors

      ! This turns on custom colors
      use_colors = .true.

      ! Points to the directory where you house the filters
      ! For 2MASS it is, H.dat, J.dat, Ks.dat
      ! Can download other filter systems using SED-tools
      instrument = '/data/colors_data/filters/2MASS/2MASS'

      ! Your choice of stellar atmosphere table
      stellar_atm = '/data/colors_data/stellar_models/Kurucz2003all/' 

      ! Distance to star in cm for synthetic photometry
      ! If you set the distance to 10 parsecs (3.0857d19), you will have absolute magnitude
      ! For any other distance, custom colors will give apparent magnitude
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)

      ! Exports a full calculated SED at every profile interval
      ! Needed if you want to plot the SED
      make_csv = .true.
      colors_results_directory = 'SED'  ! Directory the fully calculated SED will go to


      ! Defines the zero-point system for magnitude calculations
      mag_system = 'Vega'

      ! Points to the reference SED file for Vega
      vega_sed = '/data/colors_data/stellar_models/vega_flam.csv'

/ ! end of colors namelist
```
</details>

<details>
  <summary>inlist_pgstar</summary>
  &pgstar

   ! see star/defaults/pgstar.defaults

   ! MESA uses PGPLOT for live plotting and gives the user a tremendous
   ! amount of control of the presentation of the information.

   ! show HR diagram
   ! this plots the history of L,Teff over many timesteps
   HR_win_flag = .true.


   ! set window size (aspect_ratio = height/width)
   HR_win_width = 6
   HR_win_aspect_ratio = 1.0


   ! Color Color diagram
   History_Track2_win_width = 6
   History_Track2_win_aspect_ratio = 1.0

   History_Track2_win_flag = .true.
   History_Track2_xname = 'J'
   History_Track2_yname = 'Ks'
   History_Track2_title = '2MASS Magnitudes'
   History_Track2_xaxis_label = 'J mag'
   History_Track2_yaxis_label = 'Ks mag'
   History_Track2_reverse_xaxis = .true.
   History_Track2_reverse_yaxis = .true.



/ ! end of pgstar namelist

</details> 

{{< /details >}}
