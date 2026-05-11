# Lab 2 - Exploring MESA Custom Colors!


The MESA colors module allows us to generate synthetic photometry while running MESA stellar evolution models! It is a great way to merge observational and theoretical astronomy. With the colors module, we can specify what filter system and stellar atmosphere we want to use, and on top of regular MESA outputs (effective temperature, luminosity, age, etc.) we get bolometric magnitude, M$_{bol}$, bolometric flux, F$_{bol}$, and many synthetic magnitudes. For more information on the colors module, look at https://github.com/MESAHub/mesa/tree/main/colors.


One major age dating technique for stellar populations is through the use of isochrones. Isochrones are single-aged, chemically homogenous populations that show a snapshot of stellar evolution. They're made by evolving stars with the same chemical composition but different initial masses, and then finding what point in evolution each star is at at a particular age. Larger stars burn hotter and brighter, leaving the main sequence much quicker than a lower mass star. For example, at 10 Gyr we can see a 0.8 $M_{\odot}$ still on the main sequence, while a 5 $M_{\odot}$ star will be long past the Red Giant Branch. Because of this, we can build isochrones and use them to determine the age of stellar populations. One caveat to this is that they use the assumption that all the stars are at relatively the same distance and formed from the same materials at relatively the same time. _The best stellar populations to use isochrones when age dating stars is in clusters because we can make these assumptions._

This figure shows a series of isochrones at different ages between 0.03 Gyr to 10 Gyr, made using MIST. As the population gets older, the shape of the isochrone changes too!
<img width="450" height="600" alt="isochrones" src="https://github.com/user-attachments/assets/b115dfa9-6604-4531-9b04-d7dfb6481184" />



In this lab, we'll learn how atmospheric boundary conditions and the convective mixing length parameter can impact stellar evolution in both observational and theoretical coordinates. We'll also build isochrones to explore other techniques for age dating stellar populations and planet hosts.


### Step 1 - Lab 2 Prep

Last lab we made a working directory that has everything we want to start lab 2. The first thing we will do is copy lab 1 into a new working directory:

```bash
cp -r lab1 lab2
cd lab2
```
Lets clean this directory and get rid of our outputs from Lab 1:

```bash
./clean
./mk
rm -r LOGS
```
In Lab 1, we explored magnetic braking. Let's turn it off for this lab in the `&controls` section of `inlist_project`

```fortran
! Enable magnetic braking.
use_other_torque    = .false.
```

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


#### `&controls`

This is the section with the main stellar evolution parameters. Our goal is to change the stellar input parameters to see how they change evolution! Keep the same stellar mass you used in Lab 1.

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


#### `inlist_pgstar`
Now let's edit the custom plots we see as our star evolves. The first thing we are going to do is either remove or comment out the custom plots we used in Lab 1.
**Remove the following lines**:

```fortran
Grid1_plot_name(7) = 'History_Track2'

History_Track2_title = 'total AM'
History_Track2_xname = 'log_star_age'
History_Track2_yname = 'log_total_angular_momentum'
History_Track2_xaxis_label = 'log_star_age'
History_Track2_yaxis_label = 'log_total_angular_momentum'

History_Track2_ymin = 47
History_Track2_ymax = 50
```

While the star evolves, we want to see how an HR diagram compares to a color-magnitude diagram. Luckily, we have an HR diagram already there! 

**Question** Let's evolve the star using `./rn` and wait until our plotting windows pop up. What do you notice about the HR diagram?

{{< details title="Answer" closed="true" >}}

Nothing shows up in the plotting window! This is because of the lines under "! set static plot bounds" - our star is evolving outside of the log(Teff) and log(L) limits given. Let's remove the following lines:

```fortran
! set static plot bounds
 HR_logT_min = 3.5
 HR_logT_max = 4.6
 HR_logL_min = 2.0
 HR_logL_max = 6.0
```

{{< /details >}}

The HR diagram is already turned on from the command `HR_win_flag = .true.` so all we have to do is make a color magnitude diagram! Later this week you'll learn to make customize your plots, but for now just copy and paste the following into `inlist_pgstar`:

```fortran
History_Track2_win_width = 6
History_Track2_win_aspect_ratio = 1.0

History_Track2_win_flag = .true.
History_Track2_xname = 'J'
History_Track2_yname = 'Mag_bol'
History_Track2_title = '2MASS CMD'
History_Track2_xaxis_label = 'J mag'
History_Track2_yaxis_label = 'Bolometric Magnitude'
History_Track2_reverse_xaxis = .true.
History_Track2_reverse_yaxis = .true.
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

Up until now, we've been using a mixing length parameter value of 1.50. Look through the controls default parameters again and find the mixing length parameter, or $\alpha_{MLT}$. Once you locate it, try changing the value between 0 - 3 and see what happens.

{{< details title="Hint" closed="true" >}}
Check under the tab "mixing parameters" for the controls defaults
{{< /details >}}

{{< details title="Answer" closed="true" >}}
mixing_length_alpha = 1.8d0
{{< /details >}}

Because we want to compare how different values of the mixing length parameter change evolution, remember to change the name of your output history file. For example, if you are running a model with $\alpha{MLT}$ = 1.8d0, you could name your history file something like:

```fortran
star_history_name = '1p0Msun_alphaMLT1p80_history.data'
```

Pick 3 different values of $\alpha{MLT}$ and run a model for each, changing the output file name for each!


### Step 4 - Visualizing the changes outside of MESA

#### Isochrones

Go to the same [Google spreadsheet](https://docs.google.com/spreadsheets/d/1C88C5V2siCAaK8-3qgAZoNc9-9IH-RTIqFVetXQc3EM/edit?usp=sharing) as Lab 1. On the bottom, switch to the tab labeled "Lab 2". For this part, let's rerun a star with `mixing_length_alpha = 1.8d0` and the Eddington atmospheric boundary condition (`atm_T_tau_relation = 'Eddington'`). Once your star is done evolving, copy the values for "Teff" and "log(L)" from the terminal window into the Google sheet. _Make sure you are putting the values at the right corresponding age!_

As everyone finishes filling out the spreadsheet, we'll get to see the formation of an isochrone being built.

#### Comparing atmospheric boundary conditions and mixing length parameters

Now, go to the [Google Colab](https://colab.research.google.com/drive/1rFAu8UN0CC3GWllJfNyk7uV50FksOKok?usp=sharing) and make a copy of it.

Follow the instructions in the document to upload the different history files we made and see how changing one parameter can impact stellar evolution.
