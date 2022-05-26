Guide to Running the Barotropic Trough Model:

*The model contained here replicates the simplest case as examined in the paper, i.e., a single trough with downwelling-favorable conditions. Running the model to produce this single result gives the user the basic principle how to replicate all other results by changing the model parameters for which case to simulate.

****Cartoon labels representing geometry/inflow, wind forcing, and vorticity sources in the paper were added via PowerPoint. The figures produced here are the underlying graphics.****

***** Note that although W_wall and W_trough were the labels used in the paper, L_wall and L_trough are the equivalent names in the code *****


The following is an explanation of using the model to produce the single trough figure:


1) Go to the folder that contains the code.

2) Ensure that you have python with numpy, matplotlib, and scipy installed. 

3) If you don't have python installed, Anaconda.com is a good place to start. 

4) Notice any *.png or *.mp4 files. These are created when running the 'Psi_analysis.py' script, which takes the model output to compute stats and produce graphics.

5) Notice the '__pycache__' folder: ignore.

6) Notice the *.sqlite file (if it doesn't exist, then the model needs to run to produce the data!). Note the name is *_solver_output.sqlite, where * = constant 
in this case. The * (constant in this case) is an important string to mark down. This sqlite file is the data
written out by the main 'Psi_solve.py' code. If you change the model code you need to delete these. 

7) To perform a run and write out the results in the **.sqlite file, we'll want to run 'Psi_solve.py'. First
let's edit it to make sure it is a run we are interested in. Open it with something like emacs.

8) It should be decently commented to understand it. Much of what should be changed to specifiy a run is near
the top.

If We want to replicate the base run, make 'whichloop = 9' (approx. 67 lines in). If we want a single trough on the shelf so set 'troughamount = 1'. All of these parameters are hard coded.

To make the data for the parameter variation runs, choose the parameter to vary and set "whichloop" to the appropriate value. (1 = H_shelfbreak, 2 = H_trough, 3 = L_wall, 4 = L_trough, 5 = L_head, 6 = L_shelfbreak, 7 = r, 8 = f, 9 = constant, 10 = troughshift).

9) Now we are going to run it to get the results. Open up a terminal in this folder location.

10) Work within an enviornment that has packages relevant to running this code: 'conda activate envname'. If upon
doing the subsequent run, a package isn't installed, then install it. (see somewhere like https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html if unfamiliar with envs)

11) Start python with 'ipython --pylab'.

12) Run the solver script with 'run Psi_solve.py'. Let it do its thing and wait for it to finish up. A **.sqlite
file of the results should now be written in this folder. 

13) Now we want to analyze it. Before we do, let's make sure our analyze script makes sense. Open up
'Psi_analysis.py' with emacs.

14) Notice at the top there are some user input lines. When running the script, it'll ask for a few things (as
we will see). This Psi_analysis.py file should be properly set up in order to produce graphic items of interest like transect lines (used for calculating transport), streamlines, titles with
calculated metrics, etc. Go ahead and scroll around seeing what the parts 
do for plotting. Make any edits, save, and let's run it.

15) Within the 'ipython --pylab' session, run the analysis script with 'run Psi_analysis.py'.

16) It's going to ask what file to analyze. Remember the * = constant in the *_solver_output.sqlite we talked
about earlier? Input that * (in this case "constant") string exactly into the user interface terminal and hit enter. The analysis
script will look for the file '*_solver_output.sqllite', where * is what you just entered, in this working directory.
If it is there, it'll begin to break it down. If it is not because you haven't done the previous run steps, or
you entered it incorrectly, then try again.

17) After entering "constant" and hitting enter, the script is going to ask the user for input on the background
value of ATW offshore ejection %. This is because it's going to calculate enhanced offshore ejection, which 
requires subtracting off the background value (a shelf without a trough). If you need to find out what value that is, generate a zero_trough result and analyze with 0 background value. For now, enter 22.

18) Next, it'll ask for user input on whether or not to plot bathymetry contours. Seeing we want that 
streamline/bathymetry plot, input '1' to say yes.

19) It'll plot the streamline figure. If the 'Psi_analysis.py' file was edited correctly, you should see the 
red and black transects, title of 'Net Offshore Ejection = 53%', and white streamlines contoured over the 
bathymetry. If not, go back and edit, especially the section starting at line 450.

****Save the figure for future use.****


20) Close the figure. Sometimes it'll close nicely and the script will move on to asking the next questions.
Sometimes it'll be fussy. If an issue is encountered, hit 'CTRL-C' into terminal to stop it, then begin the 
analysis again, skipping this bathymetry plot by entering '0'.

21) The next question asks if you want to plot scaling analysis graphs. We don't care about them here for this single trough case,but when interested in scaling results, we would want this plotted. Enter '0' to skip.

22) After this, the analysis script should be all done and your terminal prompt is brought back to the regular
'ipython --pylab' screen.

****Notice that there is an .mp4 file saved in the folder now. This was useful when testing out scaling result to visualize each case. For this single trough case where there is a single frame, it is a pretty boring movie. Ignore.****

23) You're all done. Repeat this type of process for each figure (including panels) by making appropriate
folders for the case, making sure the files are edited correctly for the run of interest, perform the run to get the results **sqlite file, and 
analyze such results.

24) Remember to have an understanding of what the run is supposed to be, make sure the 'Psi_solve.py' is edited
accordingly, it is ran to have a corresponding output file, the 'Psi_analysis.py' is edited such that plotting
options are what is desired, and then the analysis script is run.

25) Also remember that the background ATW offshore ejection value found in a trough-
less shelf is entered manually for cases with a trough: 22%. This number should be input anytime a trough run is being analyzed. This can be determined by doing a run without a trough.











