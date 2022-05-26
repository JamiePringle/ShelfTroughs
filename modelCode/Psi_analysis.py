#-------------------------------------------------------------------#
#-----This code analyzes the results of the Psi_solver.py-----------#
#------by unpacking the sqlite files containing the data------------#
#--------------then manipulates it, plots it, etc.------------------#
#-------------------------------------------------------------------#

from pylab import *
from scipy.sparse import csc_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp2d,interp1d
import time
import pylab,numpy
import pickle
import enableParallelRuns as epr
import matplotlib.animation as animation
import operator
import matplotlib.ticker
import os.path
import pandas as ps

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section asks the user which file to bring in, and does so.

inputfilename = input("What file in this working directory do you want to analyze? (enter the '*' in '*_solver_output.sqlite') ")
print("Bringing in file " + str(inputfilename))

#Combine user input with file name, as written in Psi_solver.py
filename = str(inputfilename) + str('_solver_output.sqlite')

#Import the desired filename.
solver_data = epr.getAllData(filename)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This function is called whenever asking for using input, and wanting
#to check that the input is within valid range (nums):
def inputvaluechecker(nums,question):
    #Loop until a break statement is encountered
    while True:
        #Start an error-handling block
        try:
            #Get the user input and make it an integer
            inp = int(input(question))
            #If a ValueError is raised, it means that the input was not a number
        except ValueError:
            #So, jump to the top of the loop and start-over
            continue
        #If we get here then the input was a number; in nums?
        if inp in (nums):
            #If so, break the loop because we got valid input
            return inp

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#These functions translate positions into indicies.
def xindex(xposition):
    return int(xposition/dx)
def yindex(yposition):
    return int(yposition/dy)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section calculates analytics metrics:

#First ask the user for the net ejection metrics found in the baseline
#ATW/no-trough case (used to subtract off of results for other cases to
#determine net offshore ejection

if True:
    ATW_result = int(input("Enter % net offshore ejection found in ATW/no-trough case (enter zero if this run is of ATW case): "))

#--------------------------------------------------------------------

#Function that calculates metrics:
def transect_flows(Psi,Ly,My,xloc,troughposition,L_trough,troughamount,troughposition_one,troughposition_two,windcutoff):

    #Change dimensions [m] into indices [fencepost]
    troughposition_index = yindex(Ly-troughposition)
    troughposition_one_index = yindex(Ly-troughposition_one)
    troughposition_two_index = yindex(Ly-troughposition_two)
    L_trough_index = yindex(L_trough)
    L_trough2_index = L_trough_index/2

    #Set alongshore transect position, used to calculate transport.
    
    upwave_index = windcutoff #as specified in solver
    if upwave_index == 0:
        upwave_index = int(400e3/dy)    
    downwave_index = int(100.e3/dy) #100 km north of downwave boundary

    #Set the offshore extent of transects.
    transect_offshore = xloc
    
    #What are Psi values at transects?
    
    #At trough transects:
    if troughamount == 1:
        Psi_trough = Psi[troughposition_index,xloc]
        #Make empty storage for two trough case to avoid error undefined.
        Psi_trough_one = []
        Psi_trough_two = []        
    elif troughamount == 2:
        Psi_trough_one = Psi[troughposition_one_index,xloc]
        Psi_trough_two = Psi[troughposition_two_index,xloc]
        #Make empty storage for the one trough case to avoid error undefined.
        Psi_trough = []
        
    #At upwave boundary of domain:
    upwaveBC_index = My-1
    Psi_upwaveBC = Psi[upwaveBC_index,xloc]

    #At upwave transect, specified above:
    Psi_upwave_value = Psi[upwave_index,xloc] - Psi[upwave_index,0]

    #At downwave transect, specified above:
    Psi_downwave_value = Psi[downwave_index,xloc] - Psi[downwave_index,0]
       
    #Use the previously calculated Psi(offshore) values
    #to calculate change metric.
    #Subtract from 100 to get a percent kicked off shelf.
    #Subtract off ATW baseline (no trough value) to obtain amount _more_ than ATW.
    #Percentrage transport lost offshore by comparing downwave to upwave
    #transects, + ATW result
    transport_metric = (Psi_downwave_value/Psi_upwave_value)*100 + ATW_result

    return transport_metric,upwave_index,downwave_index,transect_offshore,Psi_upwave_value,Psi_downwave_value,Psi_trough,Psi_trough_one,Psi_trough_two,Psi_upwaveBC

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section unpacks data:

#Determine the length of the loop array applicable to the file.
loopamount = len(solver_data)

#Intialize storage objects for unpacked data and analyzer's calculations.

#Transport calculations.
Psi_upwave_value_array = []
Psi_downwave_value_array = []
Psi_trough_array = []
Psi_trough_one_array = []
Psi_trough_two_array = []
Psi_upwaveBC_array = []

#Transport metric.
transport_metric_array = zeros(loopamount)

#Parameter groups for scaling analysis.
P_one_array = zeros(loopamount)
P_two_array = zeros(loopamount)
P_three_array = zeros(loopamount)
P_four_array = zeros(loopamount)
P_five_array = zeros(loopamount)
P_six_array = zeros(loopamount)

#Passed from the solver.
Psi_array = []
Psi_upwave_array = []
Psi_downwave_array = []
Psi_offshore_array = []
Psi_coast_array = []
H_array = []
r_array = []
f_array = []
L_shelfbreak_array = []
L_head_array = []
L_trough_array = []
L_wall_array = []
H_trough_array = []
H_shelfbreak_array = []
sbsb_array = []
troughshift_array = []

#Loop over applicable array, and calculate scaling analysis parameters:
for i in range(loopamount):
    #Use pickle to load the Psi_solver save_this data for run array i.
    save_this = pickle.loads(solver_data[i][0][0])

    #Unpack the save_this dictionary into its constituent parts.
    #Inspect Psi_solver.py to understand the saved dictionary hierarchy.

    #Breakdown inputycurve.
    inputycurve = save_this['inputycurve']
    ycurve = save_this['ycurve']
    L_trough = inputycurve['L_trough']
    L_wall = inputycurve['L_wall']
    zeroy = inputycurve['zeroy']
    troughposition = inputycurve['troughposition']
    troughshift = inputycurve['troughshift']

    #Breakdown inputBath.
    inputBath = save_this['inputBath']
    H = save_this['H']
    Huntouched = save_this['Huntouched']
    L_shelfbreak = inputBath['L_shelfbreak']
    H_shelfbreak = inputBath['H_shelfbreak']
    H_trough = inputBath['H_trough']
    L_head = inputBath['L_head']

    #Breakdown inputPsi.
    inputPsi = save_this['inputPsi']
    r = inputPsi['r']
    f = inputPsi['f']
    inflowtype = inputPsi['inflowtype']
    Psi_upwave = inputPsi['Psi_upwave']
    Psi_downwave = inputPsi['Psi_downwave']
    Psi_offshore = inputPsi['Psi_offshore']
    Psi_coast = inputPsi['Psi_coast']
    sbsb = inputPsi['sbsb']
    sbbd = inputPsi['sbbd']

    #Breakdown inputextra.
    inputextra = save_this['inputextra']
    whichloop = inputextra['whichloop']
    Ly = inputextra['Ly']
    Lx = inputextra['Lx']
    ysmall = inputextra['ysmall']
    xsmall = inputextra['xsmall']
    My = inputextra['My']
    Mx = inputextra['Mx']
    dy = inputextra['dy']
    dx = inputextra['dx']
    troughamount = inputextra['troughamount']
    troughposition_one = inputextra['troughposition_one']
    troughposition_two = inputextra['troughposition_two']
    windcutoff = inputextra['windcutoff']


    #Breakdown solver_data.
    answer = solver_data[i][1]
    Psi = answer['Psi']

    #-------------------------------------------------------------

    #Put unpacked data into storage objects:
    
    #Psi and H (and shelf border shelfbreak value).
    Psi_array.append(Psi)
    H_array.append(H)
    Psi_upwave_array.append(Psi_upwave)
    Psi_downwave_array.append(Psi_downwave)
    Psi_offshore_array.append(Psi_offshore)
    Psi_coast_array.append(Psi_coast)
    sbsb_array.append(sbsb)
    
    #Looping parameter values.
    troughshift_array.append(troughshift)
    r_array.append(r)
    f_array.append(f)
    L_shelfbreak_array.append(L_shelfbreak)
    L_head_array.append(L_head)
    L_trough_array.append(L_trough)
    L_wall_array.append(L_wall)
    H_trough_array.append(H_trough)
    H_shelfbreak_array.append(H_shelfbreak)

    #------------------------------------------------------------
    
    #Call the transport metric function for all of its calculations.
    transport_metric,upwave_index,downwave_index,transect_offshore,Psi_upwave_value,Psi_downwave_value,Psi_trough,Psi_trough_one,Psi_trough_two,Psi_upwaveBC = transect_flows(Psi,Ly,My,sbsb,troughposition,L_trough,troughamount,troughposition_one,troughposition_two,windcutoff)

    #Storing transport metric for later use.
    transport_metric_array[i] = transport_metric

    #Store the transport transects for plotting, not all currently used.
    if troughamount == 0:
        #Only plotting up/downwave transects for zero trough.
        Psi_upwave_value_array.append(Psi_upwave_value)
        Psi_downwave_value_array.append(Psi_downwave_value)
    elif troughamount == 1:
        #Trough, up, and downwave transects for one trough.
        Psi_upwave_value_array.append(Psi_upwave_value)
        Psi_downwave_value_array.append(Psi_downwave_value)
        Psi_trough_array.append(Psi_trough)
    elif troughamount == 2:
        #Two trough, up, and downwave transects for two troughs.
        Psi_upwave_value_array.append(Psi_upwave_value)
        Psi_downwave_value_array.append(Psi_downwave_value)
        Psi_trough_one_array.append(Psi_trough_one)
        Psi_trough_two_array.append(Psi_trough_two)

    #Store the upwave boundary condition transect flow (same for all).
    Psi_upwaveBC_array.append(Psi_upwaveBC)

    #------------------------------------------------------------------

    #Calculate parameter value for this loop and store:
    P_one = r/(f*H_shelfbreak)
    P_one_array[i] = P_one
    P_two = H_trough/H_shelfbreak
    P_two_array[i] = P_two
    P_three = L_head/L_shelfbreak
    P_three_array[i] = P_three
    P_four = L_wall/L_trough
    P_four_array[i] = P_four
    P_five = L_trough/L_shelfbreak
    P_five_array[i] = P_five
    P_six = L_wall/L_head
    P_six_array[i] = P_six

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section sorts the above data by value:

#Function that sorts arrays, given loop variable:
def sort(whichloop):
    if whichloop == 1:
        #Creates a sorting array key based on this whichloop variable.
        #Same algorithm for the rest.
        sortarray = [pickle.loads(x[0][0])['inputBath']['H_shelfbreak'] for x in solver_data]
    elif whichloop == 2:
        sortarray = [pickle.loads(x[0][0])['inputBath']['H_trough'] for x in solver_data]
    elif whichloop == 3:
        sortarray = [pickle.loads(x[0][0])['inputycurve']['L_wall'] for x in solver_data]
    elif whichloop == 4:
        sortarray = [pickle.loads(x[0][0])['inputycurve']['L_trough'] for x in solver_data]
    elif whichloop == 5:
        sortarray = [pickle.loads(x[0][0])['inputBath']['L_head'] for x in solver_data]
    elif whichloop == 6:
        sortarray = [pickle.loads(x[0][0])['inputBath']['L_shelfbreak'] for x in solver_data]
    elif whichloop == 7:
        sortarray = [pickle.loads(x[0][0])['inputPsi']['r'] for x in solver_data]
    elif whichloop == 8:
        sortarray = [pickle.loads(x[0][0])['inputPsi']['f'] for x in solver_data]
    elif whichloop == 9:
        sortarray = zeros(1)

    elif whichloop == 10:
        sortarray = [pickle.loads(x[0][0])['inputycurve']['troughshift'] for x in solver_data]

    #Use the whichloop's sortarray to produce an index key, based on value.
    index = argsort(sortarray)

    #Output index key and the corresponding sortarray.
    return index,sortarray

#Call upon the sort function for the index sort key,
#and sorted array for the whichloop.
index, sortarray = sort(whichloop)

#--------------------------------------------------------------------------

#Make storage objects for sorting arrays:
sort_H_array = []
sort_Psi_array = []
sort_sbsb_array = []
sort_troughshift_array = []
sort_r_array = []
sort_f_array = []
sort_L_shelfbreak_array = []
sort_L_head_array = []
sort_L_trough_array = []
sort_L_wall_array = []
sort_H_trough_array = []
sort_H_shelfbreak_array = []
sort_transport_metric_array = []
sort_P_one_array = []
sort_P_two_array = []
sort_P_three_array = []
sort_P_four_array = []
sort_P_five_array = []
sort_P_six_array = []

#Loop over the whichloop array size to sort all arrays.
for i in range(loopamount):
    
    #Add to the sort array in the ascending order, per the index key.
    sort_H_array.append(H_array[index[i]])
    sort_Psi_array.append(Psi_array[index[i]])
    sort_sbsb_array.append(sbsb_array[index[i]])
    sort_troughshift_array.append(troughshift_array[index[i]])
    sort_r_array.append(r_array[index[i]])
    sort_f_array.append(f_array[index[i]])
    sort_L_shelfbreak_array.append(L_shelfbreak_array[index[i]])
    sort_L_head_array.append(L_head_array[index[i]])
    sort_L_trough_array.append(L_trough_array[index[i]])
    sort_L_wall_array.append(L_wall_array[index[i]])
    sort_H_trough_array.append(H_trough_array[index[i]])
    sort_H_shelfbreak_array.append(H_shelfbreak_array[index[i]])
    sort_transport_metric_array.append(transport_metric_array[index[i]])
    sort_P_one_array.append(P_one_array[index[i]])
    sort_P_two_array.append(P_two_array[index[i]])
    sort_P_three_array.append(P_three_array[index[i]])
    sort_P_four_array.append(P_four_array[index[i]])
    sort_P_five_array.append(P_five_array[index[i]])
    sort_P_six_array.append(P_six_array[index[i]])

#Replace the unsorted arrays with the sorted arrays.
H_array = sort_H_array
Psi_array = sort_Psi_array
r_array = sort_r_array
f_array = sort_f_array
L_shelfbreak_array = sort_L_shelfbreak_array
L_head_array = sort_L_head_array
L_trough_array = sort_L_trough_array
L_wall_array = sort_L_wall_array
H_trough_array = sort_H_trough_array
H_shelfbreak_array = sort_H_shelfbreak_array
sbsb_array = sort_sbsb_array
troughshift_array = sort_troughshift_array
transport_metric_array = sort_transport_metric_array
P_one_array = sort_P_one_array
P_two_array = sort_P_two_array
P_three_array = sort_P_three_array
P_four_array = sort_P_four_array
P_five_array = sort_P_five_array
P_six_array = sort_P_six_array
#-----------------------------------------------------------------------

#Create an array of the run case's looping variable for plotting:
loopingvariable=[]
for i in range(loopamount):
    #For the given whichloop run case...
    if whichloop == 1:
        #Print out the variable being looped over, in its ascending value.
        print('H_shelfbreak = ', H_shelfbreak_array[i])
        #Store this ascending variable value into the loop array.
        loopingvariable.append(H_shelfbreak_array[i])
    elif whichloop == 2:
        print('H_trough = ', H_trough_array[i])
        loopingvariable.append(H_trough_array[i])
    elif whichloop == 3:
        print('L_wall = ', L_wall_array[i])
        loopingvariable.append(L_wall_array[i])
    elif whichloop == 4:
        print('L_trough = ', L_trough_array[i])
        loopingvariable.append(L_trough_array[i])
    elif whichloop == 5:
        print('L_head = ', L_head_array[i])
        loopingvariable.append(L_head_array[i])
    elif whichloop == 6:
        print('L_shelfbreak = ', L_shelfbreak_array[i])
        loopingvariable.append(L_shelfbreak_array[i])
    elif whichloop == 7:
        print('r = ', r_array[i])
        loopingvariable.append(r_array[i])
    elif whichloop == 8:
        print('f = ', f_array[i])
        loopingvariable.append(f_array[i])
    elif whichloop == 9:
        print('sorting through constant')
    elif whichloop == 10:
        print('troughshift = ', troughshift_array[i])
        loopingvariable.append(troughshift_array[i])

#Done sorting.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section plots bathymetry and streamlines, if selected:

#Ask the user if they want to plot Psi on bathymetry.
nums = [0,1]
question = "Do you want to plot the bathymetry and streamlines? (0 = no, 1 = yes): "
contours = inputvaluechecker(nums,question)

#This next if statement will animate the Psi on H contour plots by looping over a makeFrame function that produces each frame. A .mp3 file is then produced locally with the name of the file's loop variable for the applicable run.

#Skip this if user inputs 0 above.
if contours == 1:
    #Initialize the figure with numeric reference and size.
    fig = figure(10,figsize=[6.0,6.0])
        
    #This function produces each frame and can be called as such from the shell.
    def makePsiFrame(n):
        #Clear the figure.
        clf()

        #The plot uses two colormaps. The following is the
        #offshore distance where the first ends and the other begins.
        colormapboundary = 0.95*L_shelfbreak

        #Label the offshore and alongshore axes.
        xlabel('offshore [km]',fontsize=12)
        ylabel('alongshore [km]',fontsize=12)

        #Pull in the H and Psi to be plotted for the nth case.
        H = H_array[n]
        Psi = Psi_array[n]

        #This array marks where to draw the depth contours.
        contour_levels = [0,5,10,15,20,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,550,600,1000,1500,2000,2500,3000,3500]
        
        #Call in the desired colormap.
        cmap = cm.gist_earth_r
        
        #Map the colormap onto the contour levels.
        colors = cmap(linspace(0,1,len(contour_levels)))
        
        #Colored in contours of H at the levels and with the colors.
        #Dividing x and y by 1000 to convert [m] to [km].
        contourf(xsmall/1000.,ysmall[50:-1]/1000.,H[51:-2,1:-1],levels=contour_levels,colors=colors)
        
        #Project the colorbar key and label it.
        cbar = colorbar()
        cbar.set_label('H [m]',fontsize=12)

        #Set to false to not plot streamlines, transects, etc.:
        if True:
            #Contour streamlines.
            CS = contour(xsmall/1000.,ysmall[50:-1]/1000.,Psi[50:-1,:],linspace(Psi.min(),Psi.max(),15),colors='w')
            #Contour streamline that goes through upwave transect edge.
            upwave_contour = [Psi[upwave_index,transect_offshore]]
            CS = contour(xsmall/1000.,ysmall[50:-1]/1000.,Psi[50:-1,:],levels=upwave_contour,colors='r')
            #Show upwave transect.
            hlines(y=upwave_index/(int(1000/dy)), xmin=1, xmax=transect_offshore/(int(1000/dx)), linewidth=2, color='r')
            #Contour streamline that goes through donwave transect edge.
            downwave_contour = [Psi[downwave_index,transect_offshore]]
            CS = contour(xsmall/1000.,ysmall[50:-1]/1000.,Psi[50:-1,:],levels=downwave_contour,colors='k')
            #Show downwave transect.
            hlines(y=downwave_index/(int(1000/dy)), xmin=1, xmax=transect_offshore/(int(1000/dx)), linewidth=2, color='k')

           
        #Give the frame a suptitle, play with options as you wish.
        if False:
            if troughamount == 0:
                suptitle('Net Offshore Ejection = ' + str(int(subtract(100,transport_metric_array[n]))) + '%', fontsize=18)
            else:
                suptitle('Additional Offshore Ejection = ' + str(int(subtract(100,transport_metric_array[n]))) + '%', fontsize=18)
                #Turning it off for now.
                suptitle('')
            
        #Give the plot a title, depending on the run case with the variable
        #and label its value per loop.
        if whichloop == 1:
            #Grand title with H_shelfbreak = value label
            title(r'$H_{shelfbreak} = $' +'%4.2e'%(loopingvariable[n],)+' [m]',fontsize=12)
        elif whichloop == 2:
            title(r'$H_{trough} = $' +'%4.2e'%(loopingvariable[n],)+' [m]',fontsize=12)
        elif whichloop == 3:
            title(r'$L_{wall} = $' +str(loopingvariable[n]/1000,)+' [km]',fontsize=12)
        elif whichloop == 4:
            title(r'$L_{trough} = $' +str(loopingvariable[n]/1000,)+'[km]',fontsize=12)
        elif whichloop == 5:
            title(r'$L_{head} = $' +str(loopingvariable[n]/1000,)+' [km]',fontsize=12)
        elif whichloop == 6:
            title(r'$L_{shelfbreak} = $' +str(loopingvariable[n]/1000,)+' [km]',fontsize=12)
        elif whichloop == 7:
            title(r'r = ' +'%4.2e'%(loopingvariable[n],)+'$[m s^{-1}]$',fontsize=12)
        elif whichloop == 8:
            title(r'f = ' +'%4.2e'%(loopingvariable[n],)+'$ [s^{-1}]$',fontsize=12)
        elif whichloop == 9:
            title('Enhanced Offshore Ejection = ' + str(int(subtract(100,transport_metric_array[n]))) + '%', fontsize=18)
            title('Net Offshore Ejection = ' + str(int(subtract(100,transport_metric_array[n]))) + '%', fontsize=18)
            #title('')
            
        elif whichloop == 10:
            title(r'troughshift = ' +'%4.2e'%(loopingvariable[n],)+'$ [m]$',fontsize=12)

        
        #Adjust tick density if interested.
        if False:
            xticks(arange(0,(Lx/1000)+1,(Lx/1000)/2))
            yticks(arange(0,(Ly/1000)+1,(Ly/1000)/10))

        #Draw and show the current loops frame.
        pause(0.5)
        draw()
        show()

        return

    #Save frames as an animation:
    if True:
        #Now that the frame is made, loop over it to animate and save. Call upon the makeFrame to make the frames over the loopamount, with interval inbetween plots, and only cycling animation to user once.
        ani = animation.FuncAnimation(fig,makePsiFrame,frames=loopamount,interval=20,repeat=False)

        #Save the animation with the outputname. See the below link for help.
        #https://trac.ffmpeg.org/wiki/Encode/H.264
        #Make sure there is a local directory called 'animations' for saving to.
        outputnamestream = inputfilename + '_streamlines.mp4'
        ani.save(outputnamestream,fps=2,extra_args=['-vcodec', 'libx264','-pix_fmt','yuv420p','-crf','18'])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section plots the scaling analysis results (ejection
#as a function of parameter):

#Ask the user if scaled parameter plots are desired.
nums = [0,1,2]
question = "Do you want to plot scaling analyses? (0 = no, 1 = yes, separate figures, 2 =  yes, on one figure): "
paramplot = inputvaluechecker(nums,question)

#This plot the six results on separate figures:
if paramplot == 1:
    #P1 plot on its specific panel, with limits, etc.
    fig = figure(11)
    clf()
    plot(P_one_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{1}$ = r/(H$_{trough}$f)',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')
    #Not saving for now.
    if False:
        fig.savefig('P_one.png')
    
    #P2 plot
    figure(12)
    clf()
    plot(P_two_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{2}$ = H$_{trough}$/H$_{shelf}$',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')
    if inputfilename=='H_trough':
        ps.DataFrame({'xaxis':P_two_array,'yaxis':subtract(100,transport_metric_array)}).to_csv('../Figure05D.csv')

    #P3 plot
    fig = figure(13)
    clf()
    plot(P_three_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{3}$ = L$_{head}$/L$_{shelf}$',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')
    fig.savefig('P_three.png')
    if inputfilename=='L_head':
        ps.DataFrame({'xaxis':P_three_array,'yaxis':subtract(100,transport_metric_array)}).to_csv('../Figure05A.csv')
    
    #P4 plot
    figure(14)
    clf()
    plot(P_four_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{4}$ = W$_{wall}$/W$_{trough}$',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')

    #P5 plot
    fig = figure(15)
    clf()
    plot(P_five_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{5}$ = W$_{trough}$/L$_{shelf}$',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')
    fig.savefig('P_five.png')
    if inputfilename=='L_trough':
        ps.DataFrame({'xaxis':P_five_array,'yaxis':subtract(100,transport_metric_array)}).to_csv('../Figure05B.csv')
    
    #P6 plot
    fig = figure(16)
    clf()
    plot(P_six_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    title('P$_{6}$ = W$_{wall}$/L$_{head}$',fontsize=20)
    ylim(0,100)
    ylabel('Additional % Ejected by Trough')
    xlabel('Parameter Value')
    fig.savefig('P_six.png')
    if inputfilename=='L_wall':
        ps.DataFrame({'xaxis':P_six_array,'yaxis':subtract(100,transport_metric_array)}).to_csv('../Figure05C.csv')
    
    tight_layout()
        
    #Draw and show the plot.
    draw()
    show()
    
#This plots these results together on one figure:
elif paramplot == 2:
    #Initialize the figure with its layout, and associated details.
    fig, axes = subplots(2, 3,figsize=(10, 6), sharex=False, sharey=True)

    #Legend labels for the panels.
    line_labels = ["P1", "P2", "P3", "P4", "P5", "P6"]

    #P1 plot on its specific panel, with limits, etc.
    l1=axes[0,0].plot(P_one_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='b')
    axes[0,0].title.set_text('P1 = r/(f*H$_{trough}$)')
    ylim(0,100)

    #P2 plot
    l2=axes[0,1].plot(P_two_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='g')
    axes[0,1].title.set_text('P2 = H$_{trough}$/H$_{shelfbreak}$')
    ylim(0,100)

    #P3 plot
    l3=axes[0,2].plot(P_three_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='r')
    axes[0,2].title.set_text('P3 = L$_{head}$/L$_{shelfbreak}$')
    ylim(0,100)

    #P4 plot
    l4=axes[1,0].plot(P_four_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='c')
    axes[1,0].title.set_text('P4 = L$_{wall}$/L$_{trough}$')
    ylim(0,100)

    #P5 plot
    l5=axes[1,1].plot(P_five_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='y')
    axes[1,1].title.set_text('P5 = L$_{trough}$/L$_{shelfbreak}$')
    ylim(0,100)

    #P6 plot
    l6=axes[1,2].plot(P_six_array,subtract(100,transport_metric_array),linestyle='--',marker='*',color='k')
    axes[1,2].title.set_text('P6 = L$_{wall}$/L$_{head}$')
    ylim(0,100)

    #Create the legend with line objects, labels,
    #position, small spacing, and title
    if False:
        fig.legend([l1, l2, l3, l4, l5, l6],    
                   labels=line_labels,  
                   loc="center right",  
                   borderaxespad=0.1,   
                   title="Parameter Key")

    #Produce x and y labels on only bottom and left sides of plot.
    for ax in axes.flat:
        ax.set(xlabel='Parameter Value', ylabel='% Ejected Shelf Flow')

    
    #Suppress labels and ticks for interior panels.
    if False:
        for ax in axes.flat:
            ax.label_outer()

    tight_layout()
        
    #Draw and show the plot.
    draw()
    show()
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section plots the radius of curvature and shear to indicate the strength of the two vorticity sources.

#Per narrowing 2002 (eq 11), the gradient of Psi along an isobath is:
# -dPsi/dp = ((d/dn(f/H))^-1)[-(r/h)*d/dn((1/H)*dPsi/dn) + (r/H^2)*dH/dn((1/H)*dPsi/dn) - 1/R*r/H((1/H)*dPsi/dn) + del x (tau_top/(rho_0*H))]
# Psi change over isobath = (1/change in vorticity )*[bottom Eckman + shear + curvature + wind]

#Ask the user if they want to plot Psi on bathymetry.
nums = [0,1]
question = "Do you want to plot the curvature and shear sources? (0 = no, 1 = yes): "
sources = inputvaluechecker(nums,question)

if sources == 1:

#There are two sections below: 1) curvature and 2) shear

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Curvature plot   
    dx=2.0*(xsmall[1]-xsmall[0])
    dy=2.0*(ysmall[1]-ysmall[0])

    #get H[0]
    H=H_array[0]
    
    #this is the arbitrary function of f(x,y) whose curvature we will examine.
    #so if we want the curvature of the bathymetry, we want this function to return
    #the bathymetry at the locations x,y. To do so, I will use an interpolation function that
    #interpolates H to a given x,y. The interp2d function returns a function that does this.
    #xsmall and ysmall are the locations of points on the H[1:-1,1:-1] grid, I think.
    #
    #OK, found one minor issue. The function interp2d only takes vectors or scalers for arguements, not
    #matrices. So I am going to assume it is called with x and y as matrices, and that x[0,:] is the
    #x-coordinates of the matrix, and y[:,0] are the y coordinates. 
    #
    #breakpoint()
    f_H=interp2d(xsmall,ysmall,H[1:-1,1:-1]) #function to interpolate H to arbitrary point
    def f(x,y):
        #see note above, this should give H on x and y if x and y are matrices created with
        #meshgrid, or something similar. 
        fout=f_H(x[0,:],y[:,0]) 
        return fout

    #breakpoint()
    
    #compute unit normal
    def unitNormal(f,xMat,yMat):
        '''
        compute unit normal component of function f(x,y) at locations given by xMat, yMat
        returns n_x and n_y the x and y components of the normal
        '''
        #the unit normal will be divergence of the function, normalized to a length of one
        #so start by computing df/dx and df/dy, the divergence
        dfdx=(f(xMat-0.5*dx,yMat)-f(xMat+0.5*dx,yMat))/dx
        dfdy=(f(xMat,yMat-0.5*dy)-f(xMat,yMat+0.5*dy))/dy

        #now compute unit normal vectors
        n_x=dfdx/sqrt(dfdx**2+dfdy**2)
        n_y=dfdy/sqrt(dfdx**2+dfdy**2)

        return n_x,n_y

    #compute gradient of feild
    def gradient(f,xMat,yMat):
        '''
        compute gradient of function f(x,y) at locations given by xMat, yMat
        returns dfdx and dfdy the x and y components of the gradient
        '''
        #the unit normal will be divergence of the function, normalized to a length of one
        #so start by computing df/dx and df/dy, the divergence
        dfdx=(f(xMat-0.5*dx,yMat)-f(xMat+0.5*dx,yMat))/dx
        dfdy=(f(xMat,yMat-0.5*dy)-f(xMat,yMat+0.5*dy))/dy

        return dfdx,dfdy

    def divergenceOfNormal(f,xMat,yMat):
        '''
        compute the divergence of the normal vector field of f(x,y) at xMat, yMat
        '''
        #divergence is d(n_x)/dx+d(n_y)/dy where n_x and n_y are the
        #components of a vector field at a point.

        #so compute d(n_x)/dx
        n_x_plus,jnk=unitNormal(f,xMat+0.5*dx,yMat)
        n_x_minus,jnk=unitNormal(f,xMat-0.5*dx,yMat)
        d_n_x_dx=(n_x_plus-n_x_minus)/dx

        #so compute d(n_y)/dy
        jnk,n_y_plus=unitNormal(f,xMat,yMat+0.5*dy)
        jnk,n_y_minus=unitNormal(f,xMat,yMat-0.5*dy)
        d_n_y_dy=(n_y_plus-n_y_minus)/dy

        return d_n_x_dx+d_n_y_dy


    def radiusOfCurvature(f,xMat,yMat):
        '''
        compute radius of curvature of function f
        radius of curvature is 1/abs((divergence of unit normal field))
    
        NOTE VERY WELL -- IF YOU HAVE A POINT WHERE THERE IS NO DIVERGENCE,
        FOR EXAMPLE ON A PLANE, THIS WILL GO TO 1/0 AND YOU MIGHT GET A NAN
    
        LIKEWISE, IF YOU CHOOSE EXACTLY A PEAK, THE RADIUS OF CURVATURE CAN 
        BE INFINITE...

        '''
        radiusCurvature=1/divergenceOfNormal(f,xMat,yMat)
        return radiusCurvature


    #===============================================================
    figure(20)
    clf()
    style.use('ggplot')

    #make xMat and yMat from xsmall and ysmall
    xMat,yMat=meshgrid(xsmall,ysmall)

    #make contour of radius of curvature
    subplot(2,1,1)
    RminKm=1.0
    radiusContoursKm=linspace(-1/RminKm,1/RminKm,25)
    R=radiusOfCurvature(f,xMat,yMat)
    contourf(xMat/1e3,yMat/1e3,1e3/R,radiusContoursKm,cmap='seismic')
    colorbar()

    #draw bathymetry
    contour_levels = [0,5,10,15,20,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,
                      425,450,475,500,550,600,1000,1500,2000,2500,3000,3500]
    contour(xMat/1e3,yMat/1e3,f(xMat,yMat),contour_levels,colors='k',alpha=0.25)
    title('1/R in km^-1 (color) of depth (lines)',fontsize='medium')
    xlabel('x in km')
    ylabel('y in km')

    #in another window, plot radius of curvature along the shelf for every 25km offshore
    subplot(2,1,2)
    xOffshore2plot=linspace(25.0e3,125.0e3,5)
    for xOff in xOffshore2plot:
        nx=argmin(abs(xsmall-xOff))
        plot(ysmall/1e3,1e3/R[:,nx],label='%dkm offshore'%(xOff/1e3))

    axis(ymin=-1/RminKm,ymax=1/RminKm)
    xlabel('alongshore distance, km')
    ylabel('1/R in 1/km')
    title('1/(Radius of Curvature), alongshore',fontsize='medium')
    legend()
    
    show()
    draw()

    savefig('RadiusOfCurvature.jpg',dpi=200)


    #======================================================
    #plot terms of pringle narrowing paper, equation 9
    figure(21)
    clf()
    style.use('ggplot')
    
    #first thing we want to do is calculate magdH/dn, the magnitude of the normal derivative of H
    dHdx,dHdy=gradient(f,xMat,yMat)
    magdHdn=sqrt(dHdx**2+dHdy**2)

    #we will want to have the unit normal to H in several places, so lets get it
    n_x,n_y=unitNormal(f,xMat,yMat) #make the unit-vector normal to the isobath

    #calculate H
    H=f(xMat,yMat)

    #get magdHdn at the upwave boundary, away from trough, and H at the same place.
    #use this to make a function that gives slope at upwave boundary as function of depth.
    #this is the alpha in Pringle Narrowing eq. 13
    Hupwave=H[-1,:]
    magdHdnUpwave=magdHdn[-1,:]
    alphaFunc=interp1d(Hupwave,magdHdnUpwave)

    #calculate alpha on same grid as H and R
    alphaFlat=alphaFunc(H.flatten())
    alpha=alphaFlat.reshape(H.shape)

    #calculate (1/H)*(1/alpha)*dHdn
    coreTerm=(1/H)*(1/alpha)*magdHdn

    #calculate -H^2/dHdn
    multAll=-H**2/magdHdn

    #curvature term, the second term of the RHS of eq. 9
    secondTerm=multAll*(1/R)*coreTerm

    #compute acceleration and shear term, the first term of the RHS of
    #eq. 9 first make normal derivative of "coreTerm" by taking its
    #gradient, and then take the dot product of the gradient with the
    #normal derivative of the bathymetry
    fjnk_coreTerm=interp2d(xsmall,ysmall,coreTerm)
    def f_coreTerm(x,y):
        fout=fjnk_coreTerm(x[0,:],y[:,0]) 
        return fout
    coreTermGrad_x,coreTermGrad_y=gradient(f_coreTerm,xMat,yMat)
    coreTerm_n=n_x*coreTermGrad_x+n_y*coreTermGrad_y
    firstTerm=multAll*coreTerm_n
    


    #now do the plots ==============================================

    clf()
    #set up colors to contour

    ax1=subplot(2,1,1)
    contMax=40.0
    eq13conts=linspace(-contMax,contMax,51)
    contourf(xMat/1e3,yMat/1e3,firstTerm,eq13conts,cmap='seismic')
    colorbar()
    contour(xMat/1e3,yMat/1e3,f(xMat,yMat),contour_levels,colors='k',alpha=0.25)
    title('Shear Source Strength',fontsize='medium')
    xlabel('offshore [km]')
    ylabel('alongshore [km]')
    
    ax2=subplot(2,1,2)
    contMax=10.0
#    eq13conts=linspace(-contMax,contMax,51)
#    eq13contsPositive=linspace(0,contMax,26)
    contourf(xMat/1e3,yMat/1e3,secondTerm,eq13conts,cmap='seismic')
    colorbar()
    contour(xMat/1e3,yMat/1e3,f(xMat,yMat),contour_levels,colors='k',alpha=0.25)
    title('Curvature Source Stength',fontsize='medium')
    xlabel('offshore [km]')
    ylabel('alongshore [km]')
    

    draw()
    show()
    savefig('bathyWithSourcesNoHighlight.jpg',dpi=200)

    #======================================================
    #plot terms of pringle narrowing paper, along an isobath

    #plot isobaths on figure 21
    for ax in [ax1,ax2]:
        sca(ax)
        whatIsobathVec=[25.0,50.0,125.0]
        contour(xMat/1e3,yMat/1e3,f(xMat,yMat),whatIsobathVec,colors='c')

    savefig('bathyWithSourcesHighlight.jpg',dpi=200)
    
    figure(22)
    clf()
    style.use('ggplot')

    nPlot=0
    for whatIsobath in whatIsobathVec:
        nPlot+=1
        subplot(3,2,2*(nPlot-1)+1)
    
        #extract terms along isobath by looking at each y and interpolating
        #in depth. ASSUMES H ONLY INCREASES OFFSHORE. We have to keep track
        #of how long the isobath is
        firstTermVec=[] #secondTerm along isobath
        secondTermVec=[] #secondTerm along isobath
#        thirdTermVec=[]  #thirdTerm along isobath
        yVec=yMat[:,0]   #distance alongshore
        xVec=xMat[0,:]   #distance across-shore in grid
        xIsobath=[]      #distance of isobath in cross-shelf
        distAlongIso=[0.0]  #distance along isobath, start at distance 0
        distThisStepVec=[0.0] #distance between each data point, for integral below
        for ny in range(len(yVec)):
            xIsobath.append(interp1d(H[ny,:],xVec)([whatIsobath])[0])
            firstTermVec.append(interp1d(H[ny,:],firstTerm[ny,:])([whatIsobath])[0])
            secondTermVec.append(interp1d(H[ny,:],secondTerm[ny,:])([whatIsobath])[0])
#            thirdTermVec.append(interp1d(H[ny,:],thirdTerm[ny,:])([whatIsobath])[0])
            if ny>0:
                #then calculate distance along isobath
                distThisStep=sqrt((xIsobath[-2]-xIsobath[-1])**2+(yVec[ny-1]-yVec[ny])**2)
                distThisStepVec.append(distThisStep)
                distAlongIso.append(distAlongIso[-1]+distThisStep)

        #now plot terms
        plot(array(distAlongIso)/1e3,firstTermVec,'--',
             label=r'shear source')
        plot(array(distAlongIso)/1e3,secondTermVec,'-',label=r'curvature source')
#        plot(array(distAlongIso)/1e3,thirdTermVec,'--',label=r'H/R/alpha')

        if nPlot==1:
            title('Scaled Cross-Isobath Transport',fontsize='medium')
            legend(fontsize=9)
        if nPlot==len(whatIsobathVec):
            xlabel('along-isobath distance [km]')
        ylabel('at %dm '%(whatIsobath))

        #now plot terms accumilated along isobath
        subplot(3,2,2*(nPlot-1)+2)
        firstTermVec=array(firstTermVec)
        secondTermVec=array(secondTermVec)
#        thirdTermVec=array(thirdTermVec)
        distThisStepVec=array(distThisStepVec)

        #we want to take a mid-point integral
        firstTermMidPoint=0.5*(firstTermVec[1:]+firstTermVec[:-1])
        secondTermMidPoint=0.5*(secondTermVec[1:]+secondTermVec[:-1])
#        thirdTermMidPoint=0.5*(thirdTermVec[1:]+thirdTermVec[:-1])

        #remove nan's from first term, from the edge.
        firstTermMidPoint[isnan(firstTermMidPoint)]=0.0
        secondTermMidPoint[isnan(secondTermMidPoint)]=0.0
        
        plot(array(distAlongIso[1:])/1e3,cumsum(firstTermMidPoint*distThisStepVec[:-1]),'--',
             label=r'shear source')
        plot(array(distAlongIso[1:])/1e3,cumsum(secondTermMidPoint*distThisStepVec[:-1]),'-',label=r'curvature source')
#        plot(array(distAlongIso[1:])/1e3,cumsum(thirdTermMidPoint*distThisStepVec[:-1]),'--',label=r'H/R/alpha')
        if nPlot==1:
            title('Cumulative Scaled Cross-Isobath Transport',fontsize='medium')
        if nPlot==len(whatIsobathVec):
            xlabel('along-isobath distance [km]')
        #ylabel('Isobath: %dm '%(whatIsobath))
        
    draw()
    show()

    savefig('alongIsobathSources.jpg',dpi=200)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Shear plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%% All Done Here %%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
