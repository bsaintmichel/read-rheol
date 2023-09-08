# ReadRheology

A free adaptation of former scripts I had on MATLAB. This prgram reads rheology files (.csv from Malvern rSpace, .csv from Anton Paar Rheocompass and .txt from TA Trios) and converts them into a nice Pandas DataFrame array. The program also provides a slicing (i.e. selecting a sub-part of the DataFrame) and several basic plotting routines for prototypical tests (e.g. flow curves, shear startup, creep, oscillatory with amplitude sweep or frequency sweep). The only modelling done here is the extraction of Herschel-Bulkley parameters from flow curve data. If you need to go further, seek [Rheos](https://github.com/JuliaRheology) (I know, it is written in Julia, but you will figure it out).

Brice Saint-Michel, Laboratoire Navier (`bsaintmichel`, still around by googlemail.com)

--------------------------------------------------------

### Input : Malvern Kinexus (rSpace) or Anton Paar .CSV file or TA Trios .txt File

File provenance is detected automatically, so no worries.

* KINEXUS (rSpace) formatting guide (i.e. what needs to be in your table) :
  * time (action), 
  * shear stress, 
  * shear rate, 
  * complex shear stress, 
  * complex shear strain,
  * frequency,
  * shear modulus (elastic),
  * shear modulus (viscous),
  * gap width,
  * normal forces,
  * torque
  * Temperature

* ANTON PAAR (RheoCompass) formatting guide :
  * Interval time,
  * Time
  * shear stress,
  * shear rate,
  * frequency,
  * storage modulus,
  * loss modulus,
  * temperature,
  * normal forces

--------------------------------------------------------

### How does the programme work ?

Basically the routine `read_rheology` will transform your .txt/.csv file into a big Pandas DataFrame, which includes all the required fields, plus a "step type" for each step. I also try to recover the step names, but I think it is not available with Anton Paar files (they don't include it, d'uh -_-). Steps are tentatively assigned automatically, which sometimes -- often -- fail. But you can also reassign quickly the steps using auxiliary functions. 

You can also efficiently slice the dataframe to select only the steps that matter to you. Things can then be plotted using one of the plotting routines (`plot_creep`, `plot_startup`, `plot_tsweep`, `plot_flowcurve`, etc.). Most of these functions have a few options to put axes in lin or log scale ; they return the plots in question while also plotting them (so, if they are a single plot, you can still add things and edit them later on). Some plots also provide additional data, such as fits to Herschel-Bulkley for flow curves because I work a lot with yield-stress fluids. But honestly RHEOS (https://github.com/JuliaRheology/RHEOS.jl) is probably a lot more adapted to this kind of work. 

Plots are done using Bokeh, because I grew tired of Matplotlib. You can save them as .png files. It would be nice to allow an output in TikZ / pgfplots, but ... well, I have a professional life.

-----------------------------------------------------

### Dependencies 

Install through Anaconda or pip : 

* Bokeh (https://bokeh.org/)
* Pandas (https://pandas.pydata.org/)
* Numpy (https://numpy.org/)
* Scipy (https://scipy.org/)
* Time (https://docs.python.org/fr/3/library/time.html) <-- you should have this already
* io (https://docs.python.org/3/library/io.html) <-- same here 
