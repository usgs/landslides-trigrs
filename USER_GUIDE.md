User guide for TRIGRS 2.1
=========================
by Rex L. Baum and Massimiliano Alvioli

Introduction
------------
TRIGRS 2.1 provides several new output options and includes minor to moderate revisions to correct errors that sometimes slowed performance in TRIGRS 2.0 (version 2.0.06b, 14 September 2009).  Although minor revisions have been ongoing since the last major release of TRIGRS (Baum and others, 2008), several significant changes to the code have resulted from work to make it more efficient to support development of a parallel version of the program (Alvioli and Baum, 2016) and to provide output to support visualization and use with the recently released USGS program, Scoops3D (Reid and others, 2015). This updated users guide was prepared to summarize cumulative changes to the TRIGRS program.  Program operation is essentially the same as in TRIGRS 2.0 (Baum and others 2008; 2010), with a few exceptions described in the following paragraphs and most information about the program in our previous publications still applies to the new version.  This document is a supplement to the user guide (Baum and others 2008).  

We made corrections to the code that computes unsaturated basal flux to eliminate certain errors (and associated error messages) that occurred during early-time computations.  We also made minor changes to reduce the time needed to compute the basal flux.  Improvements have been added to the saturated infiltration models.  The most significant of these is addition of formulas that converge rapidly for computing later time values of the finite-depth, saturated zone.  Changes were also made to remove obsolete code, prevent array boundary errors, and to ensure that counter and pointer variables are always initialized.  This guide notes minor corrections to equations published previously.  These formulas were correct in version 2.0 of the TRIGRS program; errors in the formulas were limited only to the text of reports by Baum and others (2008, 2010). Corrections to the formulas are described in detail in Baum and Godt (2013).

Changes to input file
----------------------
Adding new features, mainly output and program control options, to the TRIGRS program, version 2.1, resulted in necessary changes to the input file, “tr\_in.txt.” For experienced users, learning to use the new version of the TRIGRS program mainly requires understanding the changes to the input file.  These changes are described briefly in the following paragraphs and refer to highlighted entries in Table 1.  All parameter definitions remain the same as in TRIGRS version 2.0; any new input parameters are defined below or in Table 1.

Rather than asking the user to input the number of rows, columns, and grid cells, the TRIGRS program now reads the number of rows and columns and counts the number of data cells directly from an ASCII grid containing a digital elevation model (DEM) of the study area.  This resulted in removing “imax,” “row,” “col,” and “nwf” as inputs from lines 3 & 4 (Table 1).  We also rearranged lines 3 – 6 to achieve more logical grouping of the remaining input values.  A new option to specify a maximum slope angle appears on lines 7 – 8, so that computations can be skipped at slopes steeper than the maximum slope and default values assigned to the factor of safety and pressure head output grids.  The name and location of the elevation grid file now must be specified and appears directly after the name and location of the slope grid file (lines 21 – 22, Table 1).  

The TRIGRS program stores a small text file, “TRgrid\_size.txt,” containing the values of “imax,” “row,” “col,” and “nwf” in the same folder as the DEM grid file for future use (Table 2).  Note that the TRIGRS program stores the default value of nwf=1 because the actual number of weighting factors, “nwf,” is needed only when the runoff routing module is used.  If the user plans to use the runoff routing module, the program TopoIndex must be used first to prepare additional files needed for runoff routing.  TopoIndex counts the number of weighting factors as well as the number of rows, columns, and data cells and saves all four values to a similar file, “TIgrid\_size.txt,” to be used later by TRIGRS if the runoff routing option is activated.  The TRIGRS program always checks for the existence of “TIgrid\_size.txt” before checking for “GMgrid\_size.txt” (Table 2) or “TRgrid\_size.txt” so that the counted value of nwf, if available, will be loaded for runoff routing computations.  The file “TRgrid_size.txt” is generated only if neither of the other grid-size files already exists in the expected location.

A new option to save a grid of the depth or elevation of the computed water table was added (lines 54 – 55, Table 1).  As with most other output selections, the user specifies “T” or “.true.” to activate the option, or “F” or “.false.” to deactivate it.  In addition, the user must indicate whether the data should be saved as depth below the ground surface (“depth”), or elevation above topographic datum of the DEM (“eleva”) as shown in Table 1, line 55. Different names distinguish the water table depth grid file from the water table elevation grid file (Table 2).  In the case of multiple water tables (a normal water table at the ground surface, an inverted water table at the wetting front, and a normal, probably perched, water table at some depth below the wetting front), which may occur during rapid infiltration driven by intense rainfall, the file contains the deepest water table computed at each grid cell for the given time step.

The biggest change to TRIGRS is addition of two new file output formats for improved visualization of pressure head results and for compatibility with the USGS program Scoops3D (Reid and others, 2000; Brien and Reid, 2008, Reid and others 2015).  The TRIGRS program uses the XMDV file format (http://davis.wpi.edu/~xmdv/fileformats.html) to present computed pressure head and saturation data in combination with x-y-z coordinate values of each point where the pressure head and saturation were computed.  Data in XMDV format can be visualized in the free program VisIt (https://wci.llnl.gov/codes/visit/) and probably other 3-D visualization software.  The TRIGRS program uses an i-j-z format native to the Scoops3D program to present computed pressure head and saturation data in combination with i-j-z coordinate values of each point where the pressure head and saturation were computed.  In this format “i” and “j” represent the columns and rows in the data, counting from the lower left corner of the DEM grid.  The z values indicate elevation above topographic datum of the DEM as in the XMDV format.  Each i-j pair will have several associated z-values corresponding to different depths at which the TRIGRS program computed pressure head and saturation.  In addition to these new formats, the original list file output options are still available (Table 1, lines 60 – 61).

The user indicates output in the XMDV format or ijz format using the list output flag “flag” as indicated in lines 60 – 61 (Table 1).  As with the original list file formats provided by TRIGRS, XMDV and ijz files can become very large.  Consequently, we have provided options for files with reduced or down-sampled output as well as the output of all points.  These are specified by different numerical values of “flag.”  For example, “-7” specifies a full XMDV file, containing data for all the points at any specified timestep; “-8” produces a down-sampled file, in which points are skipped vertically according to the down-sampling interval; -9 produces a “sparse” file with only a few selected points in each depth profile.  No horizontal grid positions are skipped in either the down-sampled or sparse files.  If a down-sampled file is selected and the down-sampling interval is “3”, then TRIGRS saves every third point starting from the ground surface and working downward.  Position of any water table(s) will also appear in the down-sampled file, regardless of whether they are located on an interval point.  Any small positive integer value can be specified for the down-sampling interval.  A down-sampling interval must be entered regardless of the specified value of “flag”; however, the interval will be ignored unless flag is -8 or -5.  The “sparse” file format will include the ground surface, the computed water table, the lowest point, and may include one or two other points that define deviations of the pressure head profile from a straight line.  

As with output grid files, the XMDV and ijz file formats can be saved for specified time steps (lines 62 - 65, Table 1).  The output files are associated with the output times using ordinal numbers in the file names (Table 2).  For the example shown in Table 1, output for time 172800 s would be identified by a numeral “1” appended to the end of the file name and output for time 216000 s would be identified by a numeral “2” appended to the end of the file name (Table 2).

A few other new options appear at the end of the input file, “tr_in.txt,” (Lines 80 – 87, Table 1).  In lines 80 – 81, the user can now specify the file extension (either “.txt” or “.asc” for the grid files created by the TRIGRS program.  The main advantage of using “.asc” is that some programs, such as VisIt, automatically recognize the files as being ASCII grids.  Some users may prefer the “.txt” extension, so that the files can be readily opened by a text editor.  

In lines 82 – 83 (Table 1) the user can instruct the TRIGRS program to ignore negative pressure head computed by the saturated infiltration models (Iverson, 2000; Baum and others, 2002; Savage and others 2003) in computing the factor of safety.  If this option is active (“T” or “.true”), then negative pressure head is set to zero for the factor of safety computation.  This prevents the TRIGRS program from computing unrealistically high factors of safety for any soil having low air-entry pressure.  If the option is inactive, (“F” or “.false.”) the linear profile of negative pressure head computed when using the saturated infiltration models combined with their lack of information about soil water content results in the negative pressure head being substituted for the suction stress in the factor of safety computation.  

An option controlling the behavior of the TRIGRS program for unsaturated infiltration appears in lines 84 – 85 (Table 1).  When the unsaturated infiltration model is specified, the default behavior of the TRIGRS program is to compute the height of the capillary fringe above the initial water table and to use the saturated infiltration model at cells where the initial capillary fringe extends to the ground surface.  By activating the option in lines 84 – 85 (specify “T” or “.true.”), the TRIGRS program uses the unsaturated infiltration model at all grid cells in all zones for which the unsaturated infiltration model is specified and a consistent set of soil-water retention parameters (“Theta-sat,” “Theta-res,” and “Alpha” see lines 9 – 14 in Table 1) are provided.  All three water retention parameters must be positive and the saturated moisture content, “Theta-sat” (), must be greater than the residual moisture content, “Theta-res” (), to activate the unsaturated infiltration model for a given property zone.  The program will over-ride this option for values of dimensionless time, t* (see Baum and others, 2010, equation 5b and accompanying explanation), less than 0.2 to avoid early-time errors in the unsaturated solution.  

Lines 86 and 87 control an option that was added for convenience in preparing pressure head data for the Scoops3D program using the ijz and XMDV file formats described previously.  The option can be deactivated for all other uses by entering a negative numerical value followed by any four-character text string as shown on line 87.  For use with Scoops3D, this option allows the user to specify the depth and method of computation for estimating pressure head at depth.  The Scoops3D program uses this deep value and the basal value computed by the infiltration model to linearly interpolate pressure head at intervening depths.  The four methods of estimating the pressure head at this deep point are (1) setting it to floating point zero (“zero”), which is the default value, (2) slope parallel flow (“flow”), assuming flow parallel to the slope at the ground surface for each grid cell, (3) hydrostatic pressure (“hydr”), which makes the pressure head equal to the specified depth, and (4) relative height, (“relh”). This option compares the local topographic relief at the grid point to the total relief in the DEM to compute the relative height and uses a fraction of the relative height to reduce and smooth the hydrostatic pressure distribution.  Thus, 

		

where k,deep is the pressure head estimate for the deep point, the subscript k denotes the kth grid point, at coordinates i, j.  In equation 1, Z is the vertical coordinate direction, Zk,w is the elevation of the initial water table, Zk,deep, is the elevation of the deep point, Ek is the ground surface elevation at grid point k, Emin is the ground surface elevation of the lowest point in the DEM and Emax is the ground surface elevation of the highest point in the DEM.


*Table 1.  Listing of input file “tr_in.txt” for TRIGRS version 2.0.10*
	
1.	Name of project (up to 255 characters)
2.	TRIGRS, version 2.0.10t, Tutorial
3.	tx, nmax, mmax, zones 
4.	1,   30,   100,   2
5.	nzs,  zmin,  uww,    nper    t
6.	10,   0.001,  9.8e3,   2,  216000
7.	zmax,   depth,   rizero,  Min\_Slope\_Angle (degrees),  Max\_Slope\_Angle (deg.)
8.	-3.001,  -2.4,  -1.0e-9,       0.,    90.0
9.	zone, 1
10.	cohesion,phi,  uws,   diffus,   K-sat, Theta-sat,Theta-res,Alpha
11.	3.5e+03, 35., 2.2e+04,  6.0e-06, 1.0e-07,   0.45,    0.05,    -0.5
12.	zone, 2
13.	cohesion,phi,  uws,   diffus,   K-sat, Theta-sat,Theta-res,Alpha
14.	8.0e+03, 31., 2.2e+04,  8.0e-4, 1.0e-04,   0.45,    0.06,   -8.
15.	cri(1), cri(2), ..., cri(nper)
16.	3.e-7,  9.e-5
17.	capt(1), capt(2), ..., capt(n), capt(n+1)
18.	0,       172800, 216000
19.	File name of slope angle grid (slofil)  
20.	Data/tutorial/slope.asc
21.	File name of digital elevation grid (elevfil)
22.	Data/tutorial/dem.asc
23.	File name of property zone grid (zonfil)
24.	Data/tutorial/zones.asc
25.	File name of depth grid (zfil) 
26.	Data/tutorial/zmax.asc
27.	File name of initial depth of water table grid   (depfil)
28.	Data/tutorial/depthwt.asc
29.	File name of initial infiltration rate grid   (rizerofil)
30.	Data/tutorial/rizero.asc
31.	List of file name(s) of rainfall intensity for each period, (rifil())  
32.	Data/tutorial/ri1.asc
33.	Data/tutorial/ri2.asc
34.	File name of grid of D8 runoff receptor cell numbers (nxtfil)
35.	Data/tutorial/TIdscelGrid\_tutorial.asc
36.	File name of list of defining runoff computation order (ndxfil)
37.	Data/tutorial/TIcelindxList\_tutorial.txt
38.	File name of list of all runoff receptor cells  (dscfil)
39.	Data/tutorial/TIdscelList\_tutorial.txt
40.	File name of list of runoff weighting factors  (wffil)
41.	Data/tutorial/TIwfactorList\_tutorial.txt
42.	Folder where output grid files will be stored  (folder)
43.	Data/tutorial/
44.	Identification code to be added to names of output files (suffix)
45.	tutorial
46.	Save grid files of runoff? Enter T (.true.) or F (.false.)
47.	T
48.	Save grid of minimum factor of safety? Enter T (.true.) or F (.false.)
49.	T
50.	Save grid of depth of minimum factor of safety? Enter T (.true.) or F (.false.)
51.	T
52.	Save grid of pressure head at depth of minimum factor of safety? Enter T (.true.) or F (.false.)
53.	T
54.	Save grid of computed water table depth or elevation? Enter T (.true.) or F (.false.) followed by 'depth,' or 'eleva' 
55.	T, depth
56.	Save grid files of actual infiltration rate? Enter T (.true.) or F (.false.)
57.	T
58.	Save grid files of unsaturated zone basal flux? Enter T (.true.) or F (.false.)
59.	F
60.	Save listing of pressure head and factor of safety ("flag")? (-9 sparse xmdv , -8 down-sampled xmdv, -7 full xmdv, -6 sparse ijz, -5 down-sampled ijz, -4 full ijz, -3 Z-P-Fs-saturation list -2 detailed Z-P-Fs, -1 Z-P-Fs list, 0 none). Enter flag value followed by down-sampling interval (integer).
61.	-2,2
62.	Number of times to save output grids and (or) ijz/xmdv files
63.	2
64.	Times of output grids and (or) ijz/xmdv files
65.	172800, 216000.
66.	Skip other timesteps? Enter T (.true.) or F (.false.)
67.	F
68.	Use analytic solution for fillable porosity?  Enter T (.true.) or F (.false.)
69.	T
70.	Estimate positive pressure head in rising water table zone (i.e. in lower part of unsat zone)?  Enter T (.true.) or F (.false.)
71.	T
72.	Use psi0=-1/alpha? Enter T (.true.) or F (.false.) (False selects the default value, psi0=0)
73.	F
74.	Log mass balance results?   Enter T (.true.) or F (.false.)
75.	T
76.	Flow direction (Enter "gener", "slope", or "hydro")
77.	gener
78.	Add steady background flux to transient infiltration rate to prevent drying beyond the initial conditions during periods of zero infiltration?
79.	T
80.	Specify file extension for output grids. Enter T (.true.) for ".asc" or F for ".txt"
81.	T
82.	Ignore negative pressure head in computing factor of safety (saturated infiltration only)?   Enter T (.true.) or F (.false.)
83.	T
84.	Ignore height of capillary fringe in computing pressure head for unsaturated infiltration option?   Enter T (.true.) or F (.false.)
85.	T
86.	Parameters for deep pressure-head estimate in Scoops3D ijz output: Depth below ground surface (positive, use negative value to cancel this option), pressure option (enter 'zero' , 'flow' , 'hydr' , or 'relh')
87.	-50.0,flow
	

*Table 2.  Names and extensions of new output files generated by TRIGRS 2.1.*
[The abbreviation "TR" at the beginning of a file name identifies a file generated by the program TRIGRS; XXXXXXXX denotes the user-defined identification code (as many as eight characters) used to group files from the same run.  The pound sign, #, is the ordinal number (1, 2, 3, and so on) corresponding to the output time as in lines 62 – 65 of Table 1.  Output times are recorded in the log file.  See Table 4 of Baum and others (2008) for additional output file types.

File Name Description
---------------------

TRgrid\_size.txt
A brief, 2-line, file containing a header that identifies the data on the second line.  The file is generated automatically and saved in the same folder as the elevation grid unless the file already exists there from a previous run, or a file named TIgrid\_size.txt or GMgrid\_size.txt exists there.   
TR\_xyz\_p\_th\_XXXXXXXX\_#.okc
List of Cartesian (x-y-z) coordinates of points in the model and their pressure head, p, and volumetric water content (“th”) values preceded by a brief file header.  The file format, XMDV, is described in detail at http://davis.wpi.edu/~xmdv/fileformats.html
TR\_ijz\_p\_th\_XXXXXXXX\_#.txt
List of i-j-z (column, row, elevation) coordinates of points in the model and their pressure head, p, and volumetric water content (“th”) values.  Columns are numbered left to right and rows are numbered bottom to top, so that the lower left corner is (1,1).
TRwater\_depth\_XXXXXXXX\_#.asc
Grid of depths to the deepest computed water table for user-selected output times. File name extension may be either “.txt” or “.asc”
TRwater\_eleva\_XXXXXXXX\_#.asc
Grid of elevations to the deepest computed water table for user-selected output times. File name extension may be either “.txt” or “.asc”

Formula Corrections
-------------------

We note the following corrections to formulas appearing in our earlier publications: Replace the term “” in equations 5c, 5d, 8, 12a and 12b of Baum and others (2010) and corresponding formulas in Baum and others (2008) with “”, where KS is the saturated hydraulic conductivity, 1 is the inverse of the height of the capillary fringe adjusted for slope, and 0 is the initial suction head.  The exponential term drops out because we modified the exponential model for soil water retention from the form originally used by Srivastava and Yeh (1991) (see Savage and others, 2004; Baum and others 2010).  A negative sign should appear on the right side of equation 10a of Baum and others (2008; 2010) so that , where qZmax is the maximum basal flux and  is the vertical component of the long-term pressure-head gradient, /Z, at the basal boundary.  

Formula additions
-----------------

We added formulas to handle cases that resulted in computational errors in the original code and to improve computational efficiency.  Although the basic formulas were coded correctly in previous releases of TRIGRS (Baum and others, 2002; 2008), typographical errors did creep into the reports documenting them (Baum and others, 2008, 2010).  See Baum and Godt (2013) for corrected formulas.  New formulas added to improve the code are documented in Alvioli and Baum (2016).  Here we briefly describe the changes.  

Previous versions of the TRIGRS program required values of the input variable, zmin, to be greater than zero to avoid division-by-zero errors for the saturated, infinite depth infiltration model at the ground surface.  The formula to handle this case is adapted from Carlslaw and Jaeger (1959, p. 75) and appears in Alvioli and Baum, 2016 as equation 1. The new formula is implemented within a Heaviside step function series to accommodate multiple time steps of varying rainfall intensity (Alvioli and Baum, 2016, eq. 2).  This formula has been coded into the files iverson.f90 and ivestp.f95.  As a result, users can specify that pressure head calculations begin at the ground surface.  The option of specifying an arbitrary minimum starting depth also remains in the program.

We also found that for certain cases, the saturated, finite-depth solution converged poorly, resulting in excessive computational time and non-convergent grid cells.  To alleviate this problem, we added a later-time counterpart to the early-time formula (eq. 2 of Baum and others, 2010).  This formula is also adapted from Carslaw and Jaeger (1959, p. 112) and is given as equation 3 of Alvioli and Baum (2016) and implemented in a Heaviside series as their equation 4.  The new formula converges rapidly for later times, whereas the error function formula (equation 2 of Baum and others, 2008, 2010), which is an adaptation of Carslaw and Jaeger’s (1959) eq. 3.8(4), converges rapidly for early times.    We have implemented this expansion in the files savage.f95 and svgstp.f95.  In conjunction with addition of this new formula, we have relaxed the convergence criteria of infinite series solutions throughout the code to 1 part in 10^6 to be consistent with reporting only four significant digits in the results.

Other code changes
------------------

The main program and many of the existing subroutines were modified to make the corrections and add the new features described previously.  The following additional changes were made

Minor errors in the subroutine flux for estimating basal flux were discovered after release of version 2.0.06b and corrected in this version.  Corrections were made to pass the correct values of  (ths) and  (thr) to the subroutine smallt for making early time estimates of flux.  Efficiencies were achieved by moving computation of certain intermediate values to higher level loops to avoid repeated computation of the same value.  An error trapping routine was improved by setting basal flux to zero for cases where sufficient time has not elapsed for water moving at the saturated seepage velocity to traverse the entire distance from ground surface to initial water table.  This accelerated program operation by greatly reducing the number of occurrences of reported errors in computed basal flux.

We implemented an exact formula for transient pressure head at the ground surface (depth, z=0) in the subroutines iverson and ivestp.  Previously, these routines computed approximate values based on a very small value of z in order to avoid division by zero errors.  The newly implemented formulas are based on equations 26a – 26c of Iverson (2000); whereas formulas for z>0 continue to be based on equations 26c and 27a – 27e of Iverson (2000).  The corrected form of equation 26c as noted in Baum and others (2008, 2010) continues to be used in both subroutines.

We corrected minor errors resulting from failure to initialize certain variables for all cases.  A statement to initialize a global counter variable, nmn, used in tracking convergence of solutions for infinite depth models, was added to the TRIGRS main program.  A statement was added to the subroutine unsth to initialize a local normalized time variable, tstar1.  Statements were also added at a few places in the code to improve initialization of other variables, eliminate isolated array boundary errors, and to eliminate possible division-by-zero errors.

A conditional exit statement has been added to inner time loops of subroutines unsth and flux to cause the program to exit loops for time steps that will produce zeros.  This occurs during the early stages of multi-time-step simulations when the time differences in Heaviside series are negative.  The exits should slightly accelerate program operation by reducing the number of loop cycles executed.

We added ADJUSTL and TRIM statements to improve handling of file path names in the main program and several subroutines.  This is expected to improve compatibility with various Fortran compilers.  We also made revisions to the subroutine “trini” to correct an error in the line counting of the input file, tr\_in.txt.  The line count is used to report location (line number) of errors that occur when reading the input file.

Error trapping statements were added to the main program to eliminate errors that may occur when the user specifies too small a value of the variable “tx” to support output at multiple times.  The program will automatically increment tx until it becomes large enough to calculate unique output times.  Blocks of code using the obsolete PAUSE statement have been revised to use an equivalent form (read *).

Addition of new output options and other new features required addition of five new subroutines.  Most changes to existing subroutines were made to implement these new features.  List files, including the new ijz and XMDV formats, are now written out by the main program, rather than the subroutines. Computed pressure head and factor of safety values destined for these files are now stored in large arrays during computation to improve computational speed and efficiency.  A few changes were also made to correct minor errors that had crept in during various revisions of the program and to improve error reporting.  The functions performed by the new subroutines are described briefly here.
*	Subroutine dzero\_brac is a procedure to bracket locations of zero in a list of double precision values using change of algebraic sign between values.  This procedure is called by subroutines svijz and svxmdv to locate the water table in a profile of pressure head values.
*	Subroutine prpijz is a procedure to create file headers for various list file formats including the original TRIGRS list files showing depth profiles for each grid cell and the new ijz and XMDV list file formats.  This procedure is called by TRIGRS main.
*	Subroutine ssizgrd is a procedure to read an ascii grid file of elevations and determine its size (number of rows, columns, & data cells).  This procedure is called by TRIGRS main.  This procedure has also been incorporated into updated versions of the utility programs GridMatch, UnitConvert, and TopoIndex that are included in the TRIGRS distribution.  However, only the program TopoIndex determines the correct value of the parameter “nwf” which is used for runoff routing.  In cases where runoff routing is not used, nwf=1 can be safely used.
*	Subroutine svijz is a procedure to prepare pressure head & volumetric water content (abbreviated “th” in program output) data for export and save it in the ijz text file format used by the USGS program Scoops3D (Reid and others, 2000; Brien and Reid, 2008) to ingest pressure head data from ground-water models.  Preparations include locating all water tables in the model output at each grid cell, estimating pressure head for the deep node if output is intended for Scoops3D, and in the case of “sparse” output, finding points at maximum (+) and minimum (-) excursions from linear trend between ground surface and basal water table. The procedure also prepares water table grids for export.  This procedure is called by subroutines ivestp, svgstp, pstpf, and pstpi.
*	Subroutine svxmdv is a procedure to prepare pressure head & volumetric water content (abbreviated “th” in program output) data for export and save it in the XMDV text file format (Table 2).  Preparations include locating all water tables in the model output at each grid cell, estimating pressure head for the deep node if output is intended for Scoops3D, and in the case of “sparse” output, finding points at maximum (+) and minimum (-) excursions from a linear trend between pressure head values at the ground surface and the basal water table.  The procedure also prepares water table grids for export.  This procedure is called by subroutines ivestp, svgstp, pstpf, and pstpi.


Using the code
----------------

### Requirements ###

A modern FORTRAN compiler, MPI libraries (tested with gfortran/f95 and MPICH/Open MPI in Ubuntu Linux 12.0, CentOS 7.0 and Cygwin with GNU make; also tested under Scientific Linux 7.2 and Windows).

### How to use the software ###

+ Unpack the tar archive in your file system -  for example is `/my/dir/`, assuming a unix-like system

`cd /my/dir/`

`tar zxvf trigrs_mpi-v2.1.tgz`

`cd trigrs_mpi/`

+ Compile to generate the TRIGRS serial executable, `trg`, the parallel executable, `prg`, and the Topoindex executable, `tpx` (see TRIGRS v2.0 manual avaiable at: http://pubs.usgs.gov/of/2008/1159/downloads/pdf/OF08-1159.pdf); the suggested Makefile assumes f95 and mpif90 to be present on the system; just type:

`make`

+ You can compile the `tpx` (Topoindex), `trg` (TRIGRS v2.1, serial) and `prg` (TRIGRS v2.1, parallel) executables separately by typing, respectively:

`make tpx`

`make trg`

`make prg`

+ if object and executable files are to be removed, type:

`make clean`

+ Use Topoindex to generate, at least, the `TIgrid_size.txt` file in the same directory where dem and slope grids are stored:

`/my/dir/trigrs_mpi/tpx`

 + Run TRIGRS from the directory containing your initialization file tr_in.txt the serial execution can be started both as:

`/my/dir/trigrs_mpi/trg`

 or as:

`/my/dir/trigrs_mpi/prg`

Parallel execution, with a total of NP processes can be started as:

`mpirun -np [NP] /my/dir/trigrs_mpi/prg`

If the [NP] processes in the MPI pool are to be distributed on more than one node, a machine file mf.txt can be provided in the following form (MPICH):

`localhost:NP1`

`otherhost1:NP2`

`otherhost2:NP3`

where NP1+NP2+NP3=NP, otherhost1 and otherhost2 are accessible through ssh with no password and type:

`mpirun -np [NP] -machinefile mf.txt /my/dir/trigrs_mpi/prg`

 + TUTORIAL - you can run `tpx` and `trg` or `prg` in the source code folder `trigrs_mpi/`. The sample initialization files (`tpx\_in.txt` for `Topoindex` and `tr\_in.txt`) will be used, input data will be read from the existing folder `trigrs_mpi/data/tutorial` and output data will be stored in `trigrs_mpi/data/output/`. Modify the initialization files to suit your input data and needs.

Acknowledgements
----------------

Mark Reid and Dianne Brien (both USGS) provided helpful advice and information with regard to implementation of code to export ijz and xmdv data.  Salvatore Raia (formerly of CNR IRPI, Perugia, Italy) and Soni Yatheendradas (NASA) identified several minor issues with earlier versions of the program and made suggestions for improving the code that have been addressed in this revision.  Mateo Berti identified the errors in the formulas presented in Baum and others (2008; 2010).  Collaboration with Massimiliano Alvioli (CNR, IRPI, Perugia, Italy) resulted in minor code restructuring to separate processing from output of the list files to facilitate implementation of parallel processing using the Message Passing Interface (MPI).

References cited
----------------
+ Alvioli, M. and Baum, R.L., 2016, Parallelization of the TRIGRS model for rainfall-induced landslides using the message passing interface: Environmental Modelling & Software, v. 81, p. 122 - 135, doi: 10.1016/j.envsoft.2016.04.002

+ Baum, R.L., Savage, W.Z., and Godt, J.W., 2002, TRIGRS--A FORTRAN Program for Transient Rainfall Infiltration and Grid-Based Regional Slope-Stability Analysis: U.S. Geological Survey Open-File Report 02-0424, 35 p., 2 appendices.

+ Baum, R.L., Savage, W.Z., and Godt, J.W., 2008, TRIGRS—A Fortran program for transient rainfall infiltration and grid-based regional slope-stability analysis, version 2.0: U.S. Geological Survey Open-File Report, 2008-1159, 75 p.

+ Baum, R. L., Godt, J.W., and Savage, W. Z., 2010, Estimating the timing and location of shallow rainfall-induced landslides using a model for transient, unsaturated infiltration: Journal of Geophysical Research, Earth Surface. v. 115, F03013, doi:10.1029/2009JF001321

+ Baum, R.L., and Godt, J.W., 2013, Correction to “Estimating the timing and location of shallow rainfall-induced landslides using a model for transient, unsaturated infiltration” Journal of Geophysical Research Earth Surface, v. 118, DOI: 10.1002/jgrf.20100  

+ Brien, D.L., and Reid, M.E., 2008, Assessing deep-seated landslide susceptibility using 3-D groundwater and slope-stability analyses, southwestern Seattle, Washington in Baum, R.L., Godt, J.W., and Highland, L.M., eds., Engineering geology and landslides of the Seattle, Washington, area: Geological Society of America Reviews in Engineering Geology v. XX, p. 83-101, doi: 10.1130/2008.4020(05).

+ Carslaw, H.S., and Jaeger, J.C., 1959, Conduction of Heat in Solids (2d ed.): New York, Oxford University Press, 510 p.

+ Iverson, R.M., 2000, Landslide triggering by rain infiltration: Water Resources Research, v. 36, no. 7, p. 1,897–1,910.

+ Savage, W.Z., Godt, J.W., and Baum, R.L., 2003, A model for spatially and temporally distributed shallow landslide initiation by rainfall infiltration, in Rickenmann, D. and Chen, C., eds., Debris-Flow Hazards Mitigation—Mechanics, Prediction, and Assessment: Rotterdam, Millpress (Proceedings of the 3rd International conference on Debris Flow Hazards, Davos, Switzerland, September 10-13, 2003), p. 179-187.

+ Savage, W.Z., Godt, J.W., and Baum, R.L., 2004, Modeling time-dependent aerial slope stability, in Lacerda, W.A., Erlich, M., Fontoura, S.A.B., and Sayao, A.S.F., eds., Landslides—Evaluation and stabilization, Proceedings of the 9th International Symposium on Landslides: London, A.A. Balkema Publishers, v. 1, p. 23–36.

+ Reid, M.E., Christian, S.B., and Brien, D.L., 2000, Gravitational stability of three-dimensional stratovolcano edifices: Journal of Geophysical Research, v. 105, no. B3, p. 6043-6056.

+ Reid, M.E., Christian, S.B., Brien, D.L., and Henderson, S.T., 2015, Scoops3D—Software to analyze 3D slope stability throughout a digital landscape: U.S. Geological Survey Techniques and Methods, book 14, chap. A1, 218 p., http://dx.doi.org/10.3133/tm14A1.

+ Srivastava, R., and Yeh, T.-C. J., 1991, Analytical solutions for one-dimensional, transient infiltration toward the water table in homogeneous and layered soils: Water Resources Research v. 27, p. 753–762.
