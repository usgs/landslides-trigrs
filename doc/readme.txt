README       19 August 2009



TRIGRS

 A Fortran program for analyzing time-dependent rainfall infiltration and slope stability in a digital landscape





TRIGRS, version 2.0.06b, 14 Sep. 2009 (public version for TRIGRS 2.0)



  This distribution includes source code and executable files for the program TRIGRS and three companion utility programs, TopoIndex, GridMatch, and UnitConvert.  Executable files are available for Windows 98/NT/2000/XP/2003 and for Mac OSX version 10.2 or higher.  The code can be compiled for Unix or other platforms that support Fortran 90/95 and Fortran 77.  



CONTENTS

DISTRIBUTION FILE

DOCUMENTATION

EXTRACTING AND INSTALLING THE FILES

SOFTWARE DESCRIPTION

  Overview

  UnitConvert
  
  GridMatch

  TopoIndex

  TRIGRS

SUPPORT

DISCLAIMER

VERSION HISTORY



*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-



DISTRIBUTION FILE

  The following distribution files are available:

Documentation: OF20081159.pdf

For Windows:  OF20081159_W.zip

For Mac OS X: OF20081159_M.zip



  These files are available from the USGS at

http://pubs.usgs.gov/of/2008/1159/



  The latest version of the programs is available at 

http://landslides.usgs.gov/research/software.php 
or by contacting Rex Baum <baum@usgs.gov>



DOCUMENTATION



  The documentation for TRIGRS is contained in the text of U.S. Geological Survey Open-File Report 2008-1159, which is saved in PDF format as "OF20081159.pdf". The documentation can be viewed using Adobe Acrobat Reader.  This documentation includes a user guide and tutorial.  





EXTRACTING AND INSTALLING THE FILES



  The computer programs are saved as separate ".zip" archives for Windows PC and Macintosh computers and must be downloaded separately.  These archives contain the compiled applications, make files, source code, and sample data.  


  To use the PC version, download the file "OF20081159_W.zip" to any convenient folder the hard disk and double click on its icon to extract the files.  Use the dialog to select the desired target folder and click "OK" to extract the TRIGRS folder and its contents into the target location. The documentation can be downloaded and copied to the folder "doc" within the "TRIGRS" folder.  The PC version of TRIGRS and its companion programs require Windows 98 or later and compatible hardware.  Memory and disk requirements vary as explained in the documentation.



  To use the Mac version, download the file "OF20081159_M.zip" to any convenient folder on the hard drive, preferably within your "home" folder.  Double click on the icon to extract the files.  Double-click the icon of the folder, "TRIGRS_2" to see the contents.  The documentation can be downloaded and copied to the folder "Documentation" within the TRIGRS_2 folder. The Mac version of TRIGRS and its companion programs require Mac OS X, version 10.2 or later and run on the command line.  See OF20081159.pdf the "User interface" section of this readme file for details.  Memory and disk requirements vary as explained in the documentation.



SOFTWARE DESCRIPTION



Overview

  TRIGRS is a tool to be used by investigators who have some knowledge and experience concerning landslide behavior. Selecting input parameters and interpreting results requires geologic and engineering judgment as well as common sense.  The user should understand the theory and limitations behind TRIGRS, which are outlined in the documentation.  The user must also be aware of the limitations of the digital elevation model, physical properties and hydrologic data that TRIGRS uses as input for analyses.



 TRIGRS and its companion utility programs, TopoIndex, GridMatch, and UnitConvert, run in simple input-output windows and have limited user interaction.  Each program uses an initialization file that contains basic data needed to run the program as well as the names of other input files.  At the beginning of a project, the user would prepare a digital elevation model, slope grid, flow-direction grid, and physical properties zone grid files using Geographic Information System (GIS) software.  Next, if needed, the user would run UnitConvert to put the grid files into a consistent system of units. Use GridMatch to ensure that the grid files are congruent.  A common source of errors running TRIGRS is grids that contain the same number of cells, but the number or locations of no-data cells differ slightly because they have been generated from different sources.  The user should then run TopoIndex using the digital elevation model and direction grid file as input to define the flow distribution pattern and compute weighting factors for distributing surface runoff.  Routing and distribution of surface runoff is optional; however, TopoIndex will compute array sizes from the digital elevation model to be used by TRIGRS, and so the user should run TopoIndex once at the beginning of a project to determine array sizes.  After preparing the topographic and physical properties data, converting the data to consistent units, computing array sizes and (if desired) preparing the flow distribution data, the project is ready for analysis using TRIGRS.  The user may run TRIGRS as many times as necessary to establish a time-series response of shallow pore water and stability of shallow slope deposits to rainfall infiltration.  Each run of TRIGRS computes the pore pressure and factor of safety over the grid for a particular, user-specified, instant in time.    



User interface 
  All four programs have a command-line interface, however, the command line works differently on different platforms.
  
PC  
  When the user double-clicks the program icon, the program immediately launches and looks for the initialization file, which should be in the same folder as the program.  If a program does not find its default initialization file (tr_in.txt, tpx_in.txt, gm_in.txt, or uc_in.txt), it will prompt the user to enter the name of the file.  Once the program finds an initialization file, it will open and read it, open any needed input files that were listed in the initialization file, complete the specified computations, save the computations to files and quit.  In the event of input or output errors, error messages will appear on screen to help the user trace the source of the problem.  Please consult the tutorial (in the main report describing TRIGRS) for details of preparing the data and initialization files.


Mac OS X:
  On Mac OS X, the programs run inside the terminal application, which resides in the folder /Applications/Utilities.  Double click the Terminal icon to launch it and open a new terminal window.  If you are familiar with Posix path names, then type "cd /Users/your_username/path/to/the/folder/TRIGRS_2" otherwise type "cd", press the space bar, and then use the mouse to drag the icon of the TRIGRS_2 folder onto the Terminal window.  Release the mouse button and the path of the folder should appear on the command line.  Press return to change the present working directory to the TRIGRS folder.  Next use your mouse to open "bin" folder in TRIGRS_2 and drag one of the program icons to the Terminal window, relase the mouse button to copy the path name of the program to the command line (or type, for example,  "/Users/your_username/path/to/the/folder/TRIGRS_2/bin/trigrs.out") and press return.  The program immediately launches and looks for the initialization file, which should be in the present working directory.  If a program does not find its default initialization file (tr_in.txt, tpx_in.txt, gm_in.txt, or uc_in.txt), it will prompt the user to enter the name of the file.  Once the program finds an initialization file, it will open and read it, open any needed input files that were listed in the initialization file, complete the specified computations, save the computations to files and quit.  In the event of input or output errors, error messages will appear in the terminal window to help the user trace the source of the problem.  Please consult the tutorial (in the main report describing TRIGRS) for details of preparing the data and initialization files.  




UnitConvert

  Use UnitConvert at the beginning of a project if you need to convert the data in a grid from one system of units to another.  UnitConvert reads the contents of an input grid file, multiplies each grid value by a conversion factor, and saves the result to the output grid file.  Input and output grids have the same resolution.  UnitConvert uses an initialization file and one input file of data in ASCII grid format; it produces one output file of data, also in ASCII grid format.  The initialization file contains three lines input data that contain the name of the input file, the name of the output file and the conversion factor.  A descriptive header appears above each line of input data to aid the user.



GridMatch


  Use GridMatch at the beginning of a project to verify that all grid files contain the same number of no-data cells and that the no-data cells are in the same locations.  GridMatch reads the contents of an input file to determine which grids to compare.  The first grid file listed in the input file is treated as a master file and all other grids are compared to the master.  The program GridMatch lists the results of the file comparisons in an ascii text log file.  GridMatch does not attempt to correct mismatched grids. Any needed corrections or adjustments to the grid files must be made with the GIS software used to generate the grids.  Once the user has confirmed that all grid files are correctly matched, the grids are ready for further analysis using TopoIndex and TRIGRS.




TopoIndex

  Use TopoIndex at the beginning of a project to compute array sizes for TRIGRS and to prepare runoff routing data if you want TRIGRS to route excess water to downslope cells.  It should be necessary to run TopoIndex only once for a project after you have chosen a routing method and created a hydrologically consistent DEM and flow-direction grid.  TopoIndex uses an initialization file and two input files to determine the correct order for runoff routing calculations and to compute the weighting factors that determine how excess water is distributed to neighboring downslope grid cells.  Only one of the input files (a digital elevation model of the study area) is needed to determine array sizes for running TRIGRS if the user chooses to neglect runoff routing.





TRIGRS

  Use TRIGRS to compute the shallow time-dependent pore-pressure and slope-stability response of an area subject to rainfall.  TRIGRS is based on superposition of one-dimensional infiltration on a general steady flow field to determine pore-pressure response and an infinite-slope model to determine stability of shallow soils during storms. TRIGRS uses an initialization file and one or more input files, depending on the variability of physical properties and rainfall in the study area.  If desired, TRIGRS can use a simple (built-in) runoff-routing model to divert excess water to downslope areas where it has the opportunity to infiltrate.  
  TRIGRS 2 includes models for saturated and unsaturated initial conditions.  For saturated initial conditions, the models are the same as used in TRIGRS 1.0 and each simulation with TRIGRS analyzes a single, user-specified point in time during a storm sequence.  Thus, the user runs a series of simulations for different specified times using TRIGRS to study the time history of slope instability during a rainfall sequence.  For unsaturated initial conditions, TRIGRS 2 can analyze a complete time series  during a single simulation.    





SUPPORT



  There is no formal ongoing support for this freely distributed public domain software.  However, we are interested in feedback. If you find errors or have suggestions, please contact:



Rex Baum baum@usgs.gov 



DISCLAIMER



  This open-file report was prepared by an agency of the United States Government. Neither the United States Government nor any agency thereof nor any of their employees makes any warranty, expressed or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed in this report or represents that its use would not infringe privately owned rights. Reference therein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof. 



  Although all data and software in this open-file report have been used by the USGS, no warranty, expressed or implied, is made by the USGS as to the accuracy of the data and related materials and (or) the functioning of the software. The act of distribution shall not constitute any such warranty, and no responsibility is assumed by the USGS in the use of these data, software, or related materials. 





VERSION HISTORY
  TRIGRS 2.0.06b 14 Sep. 2009, contains minor revisions to (1) work around an error in the output from the saturated infiltration model that results in zero pressure-head output when multiple timesteps with zero surface infiltration occur in succession.  Pressure head from the previous timestep is now substituted in the output.  (2) A user option to offset the transient surface infiltration by the amount of the steady background flux was added.  This eliminates the increase of suctions beyond the initial conditions during storm sequences that contain time steps with surface infiltration (or rainfall) that is less than the steady background flux.  This option adds an input value to the last line of tr_in_txt. (3) added formula rn=float(n) to Heaviside series "b" in subroutines pstpi() and pstpf().  This was previously omitted in error.

  TRIGRS 2.0.05  25 Feb 2009, contains minor revisions to (1) correct the formula for beta, which limits the maximum computed pressure head to agree with the value for vertical infiltration given by Iverson (2000), and (2) change the criterion for computing the Chi correction factor for effective stress in the unsaturated zone from pressure heads that are less than the air entry value, -1/alpha, to pressure heads less than zero.  This change eliminates a discontinuity in the factor of safety at the air entry value. 

  TRIGRS 2.0.04  12 Feb 2009, contained revision #1 of version 2.0.05, but was not released publicly.

  TRIGRS 2.0.03  29 July 2008, contains minor revisions to (1) fix a memory allocation error, (2) correct values in output headers for grids of non-convergent cells when no property zone input file is used, and (3) change the file extension of output grids from ".asc" to ".txt"
  TRIGRS 2.0.00  27 May 2008, initial public version of TRIGRS 2.0 

  TopoIndex 1.0.03, 27 May 2008, changed the file extension of output grids from ".asc" to ".txt"
  TopoIndex 1.0.03, 27 May 2008, Improved handling of "nodata" values to allow different nodata values for integer and floating point grids (Thanks to Manfred Thüring, Institute of Earth Sciences -SUPSI, Ticino, Switzerland for pointing out this problem).  Also converted the subroutine nxtcel() to Fortran 90/95 to improve its loop structure. Served as Beta Version of TopoIndex for use with TRIGRS 2.0

  GridMatch 1.0.00, 27 May 2008, version of GridMatch for use with TRIGRS 2.0

  UnitConvert 1.0.02 27 May 2008, version of UnitConvert for use with TRIGRS 2.0
