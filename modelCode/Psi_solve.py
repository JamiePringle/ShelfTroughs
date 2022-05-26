#---------Defining Problem--------------------------------
#This solver computes the Psi (streamfunction) solution to a potential vorticity equation with the following assumptions: steady-state, barotropic, linear (small Rossby):
# 0 = J(Psi,f/H) + del*((r/H**2)*del(Psi)) - delx(tau^top/(rho_0*H))
#---------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Important User Note!!!!                                                %
# All options are being hard coded and no longer allowing for user choice%
# ,i.e., no user input                                                   %
# Questions that did ask for user input will be under an if False:       %
# Runs that produced the figures for paper are hard coded now            %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Bring in required packages and such.
from pylab import *
from scipy import *
from scipy import sparse
from scipy.sparse import csc_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import time
import pylab,numpy
import pickle
import enableParallelRuns as epr
import sys
import matplotlib.ticker 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This first section specifies values which determine the run type. Read
#along with section comments to understand each piece. Note from above that
#there are options which allow for user input to specifiy run type. However,
#I am turning that off and hard-coding options which produce the figures of
#the paper.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Function that checks for a valid user input (N/A currently):
#Given the question to the user, and possible answers of options (nums), did they enter valid choice?
def inputvaluechecker(nums,question):
    #Loop until a break statement is encountered.
    while True:
        #Start an error-handling block
        try:
            # Get the user input and make it an integer.
            inp = int(input(question))
            #If a ValueError is raised, it means that the input was not a number
        except ValueError:
            # So, jump to the top of the loop and start-over
            continue
        #Is the user input within the num of possible entries?
        if inp in (nums):
            # If so, break the loop because we got valid input, all set.
            return inp

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section defines parameters that fuel the rest of the run.

#Flows are of an upwave wind-generated type; see line ~583 for other choices.
inflowtype = 3

#Flows are confined to the shelf region; see line ~593 for other choices.
inflowregion = 1

#Run type: (1 = H_shelfbreak, 2 = H_trough, 3 = L_wall, 4 = L_trough, 5 = L_head, 6 = L_shelfbreak, 7 = r, 8 = f, 9 = constant, 10 = troughshift): 
whichloop = 2

#Do you also want downwave  wind forcing, in addition to upwave wind-driven flows? (0 = no, 1 = yes):
Add_wind = 0

#troughamount: (0 = none, 1 = one, 2 = two):
troughamount = 1

#Change what's appropriate.
if troughamount == 0:
    troughamount = 1
    notrough = 0    
elif troughamount == 1:
    notrough = 1
elif troughamount == 2:
    notrough = 2
        

#Defining the domain size and resolution (x offshore, y alongshore):
Lx = 200.e3
dx = 100.*5
Ly = 500.e3
dy = 100.*5*2

#Using domain parameters to construct grid:
Mx  = int((Lx/dx)+1) #offshore fenceposts
xbig = numpy.linspace(0,Lx,Mx+2) #offshore array
xsmall = xbig[1:-1] #offshore array minus edges (see bathymetry def.)
zerox = numpy.linspace(-Lx/2,Lx/2,Mx+2) #used for constructing winds
My = int((Ly/dy)+1) #alongshore fenceposts          
ybig = numpy.linspace(0,Ly,My+2) #alongshore array   
ysmall = ybig[1:-1] #alongshore array minus edges (see bathymetry def.)         
zeroy=numpy.linspace(-Ly/2,Ly/2,My+2) #used for constructing trough
center = int(My/2) #central alongshore point 
N = Mx*My #Total domain cells

#Physical parameters:
r = (5.0*10**-4) #frictional parameter [m/s]            
f = (10**-4) #Coriolis frequency [1/s]
rho_0 = 1000. #standard water density [kg/m3] 
slopetransport = 5 #confines slope flows (small->large,broad->confined)
tauy_mag = -0.05 #strength of alongshore winds
flipwindcurl = 1. #flip the wind curl
flipForcing = -1. #flip the forcing Psi direction

#Scaling for forcing magnitudes (if not using wind inflowtype above)
constant_transport = 1.0*flipForcing
constant_velocity = 0.05*flipForcing

#More bathymetry parameters:
beachdepth = 10. #nonzero coastal depth [m]
halfdelh = 1500.*1 #half of total slope depth [m]
xwidthsb = 100.*1 #slope's offshore width (somewhat arbitrary)
L_shelfbreak = 150.e3 #shelf width [m]
H_shelfbreak = 100.*1 #shelfbreak depth [m]
shelfslope = (H_shelfbreak - beachdepth)/L_shelfbreak #slope of shelf [m/m]

#Trough bathymetry parameters:
H_trough = 250.*1*notrough #trough depth [m]
L_head = 20.e3 #offshore width between coast and trough [m]
fjord = 0 #include a fjord on the coast?
if fjord == 1:
    L_head = 0.
L_trough = 50.e3 #alongshore trough width (W_trough) [m]
L_wall = 20.e3 #alongshore trough wall widths (W_wall) [m]
L_trough2 = L_trough/2 #splitting trough width in half for two trough cases  
troughposition = Ly - 250.e3 #alongshore position of trough
troughshift = 100.e3 #alongshore seperation of two trough case
troughposition_one = Ly - (250.e3 + troughshift) #for two trough case
troughposition_two = Ly - (250.e3 - troughshift) #for two trough case
upwavewindpos_notrough = troughposition + L_trough #position of trough edge

#These arrays are used in variation runs (physical parameters varied)
r_array = numpy.arange(1.*10**-4,10.**-3+1*10**-4,10.**-4)        
f_array = numpy.arange(1.*10**-5,1.5*10**-4+10**-5,10.**-5)
L_shelfbreak_array = numpy.arange(100.e3,180.e3+1,5.e3) 
H_shelfbreak_array = numpy.arange(50.,401.,50.)          
H_trough_array = numpy.arange(0.,601.,50.)              
L_head_array = numpy.arange(0,150.e3+1,10.e3)              
L_trough_array = numpy.arange(0.e3,150.e3+1,10.e3)        
L_wall_array = numpy.arange(5.e3,50.e3+1,5.e3)
troughshift_array = numpy.arange(20.e3, 110.e3+1.e2, 10.e3)  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section has user input questions for parameters: turned off now.

if False:
    #Call for user input.
    nums = [1,2,3,4,5,6,7,8,9,10]
    question = "Please enter a number between 1 and 9 to determine which run to perform. (1 = H_shelfbreak, 2 = H_trough, 3 = L_wall, 4 = L_trough, 5 = L_head, 6 = L_shelfbreak, 7 = r, 8 = f, 9 = constant, 10 = troughshift (only for two troughs): "
    whichloop = inputvaluechecker(nums,question)

    #Inform user of selection.
    if whichloop == 1:
        print("H_shelfbreak run selected")
    elif whichloop == 2:
        print("H_trough run selected")
    elif whichloop == 3:
        print("L_wall run selected")
    elif whichloop == 4:
        print("L_trough run selected")
    elif whichloop == 5:
        print("L_head run selected")
    elif whichloop == 6:
        print("L_shelfbreak run selected")
    elif whichloop == 7:
        print("r run selected")
    elif whichloop == 8:
        print("f run selected")
    elif whichloop == 9:
        print("constant run selected")
    elif whichloop == 10:
        print("troughshift run selected")

if False:
    #Ask whether or not include downwave winds.    
    nums = [0,1]
    question = "Do you also want downwave  wind forcing, in addition to upwave wind-driven flows? (0 = no, 1 = yes): "
    Add_wind = inputvaluechecker(nums,question)

if False:
    #Ask for input to determine the amount of troughs on the shelf.
    nums = [0,1,2]
    question = "Please enter a number between 0 and 2 to determine the amount of troughs. (0 = none, 1 = one, 2 = two): "
    troughamount = inputvaluechecker(nums,question)
    
    #Constants used in bathymetry defining. See below.
    if troughamount == 0:
        troughamount = 1
        notrough = 0
        print("Zero trough shelf")
    elif troughamount == 1:
        notrough = 1
        print("One trough shelf")
    elif troughamount == 2:
        notrough = 2
        print("Two trough shelf")

#Function that either passes standard domain values along or asks user.
def domaininput(moad):
    if moad == 1:
        print("Using standard domain values of Lx = ", Lx, ",dx = ", dx, ", Ly = ", Ly, "dy = ", dy)
    elif moad == 2:
        #Asking for user input for domain space values.
        Lx = float(input("Please enter offshore distance [m]: "))
        dx = float(input("Please enter offshore grid size [m]: "))
        Ly = float(input("Please enter alongshore distance [m]: "))
        dy = float(input("Please enter alongshore grid size [m]: "))
        #Feedback values to screen.
        print("Using standard domain values of Lx = ", Lx, ",dx = ", dx, ", Ly = ", Ly, "dy = ", dy)
    return Lx,dx,Ly,dy

if False:
    #Ask user whether to use domain values in the script or manual input.
    nums = [1,2]
    question = "Would you like to use standard domain parameters (Lx,dx,Ly,dy), or specify manually via this input? (1 = values from script, 2 = enter manually) "
    Lx,dx,Ly,dy = domaininput(inputvaluechecker(nums,question))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section constructs the bathymetry.

#These functions turn x,y positions into their indices
def xindex(xposition):
    return int(xposition/dx)
def yindex(yposition):
    return int(yposition/dy)

#This function will smooth a list (used for smoothing hard edges)
def smoothList(list,strippedXs=False,degree=10):  
    if strippedXs == True: return Xs[0:-(len(list)-(len(list)-degree+1))]  
    smoothed=[0]*(len(list)-degree+1)  
    for i in range(len(smoothed)):  
        smoothed[i] = sum(list[i:i+degree])/float(degree)  
    return smoothed

#This function takes parameter values from inputycurve and defines the trough's alongshore geometry.
def ycurvefn(inputycurve):
    
    #Unpacks inputycurve dictionary:
    L_trough = inputycurve['L_trough']
    L_wall = inputycurve['L_wall']
    zeroy = inputycurve['zeroy']
    troughposition = inputycurve['troughposition']
    troughposition_one = inputycurve['troughposition_one']
    troughposition_two = inputycurve['troughposition_two']
    Ly = inputycurve['Ly']
    troughshift = inputycurve['troughshift']

    #Defining this like the standard case for each loop
    troughposition_one = Ly - (troughposition + troughshift)
    troughposition_two = Ly - (troughposition - troughshift) 

    #Function to make alongshore trough shape, given dimensions:
    def ycurve(L_trough,L_wall,zeroy):    
        zeroycurve = 1.0*(1-((1-0.5*(1-numpy.tanh((zeroy-L_trough/2)/L_wall)))+
             (1-0.5*(1-numpy.tanh((-zeroy-L_trough/2)/L_wall)))))
        return zeroycurve
    #Defines around Ly=0, projects back to positive Ly domain.
    troughposition = (Ly/2)-troughposition
    zeroycurve = ycurve(L_trough,L_wall,zeroy-troughposition)

    #Constructing alongshore shape for two troughs.
    L_trough = L_trough/2
    troughposition_one = (Ly/2)-troughposition_one
    troughposition_two = (Ly/2)-troughposition_two    
    zeroycurve_one = ycurve(L_trough,L_wall,zeroy-troughposition_one)
    zeroycurve_two = ycurve(L_trough,L_wall,zeroy-troughposition_two)
    twozeroycurve = zeroycurve_one + zeroycurve_two

    #Using one trough or two trough bathymetry?
    if troughamount == 1:
        ycurve = zeroycurve
    elif troughamount == 2:
        ycurve = twozeroycurve

    return ycurve

#This takes alongshore y shape and makes overall domain bathymetry:
def Bathymetry_Construction(inputBath):

    #Unpack inputBath dictionary:
    ycurve = inputBath['ycurve']
    H_trough = inputBath['H_trough']
    L_head = inputBath['L_head']
    L_shelfbreak = inputBath['L_shelfbreak']
    H_shelfbreak = inputBath['H_shelfbreak']

    #Updating slope of shelf for when shelfbreak is changed.
    shelfslope = (H_shelfbreak - beachdepth)/L_shelfbreak

    #H grid is a cell larger than x,y (Psi grid) on all sides b/c of BCs.
    #Note: Psi_i,j aligns with H_i+1,j+1.
    H = numpy.zeros((My+2,Mx+2))

    #Define offshore direction bathymetry of the form: H(x) = mx + ntanh(x)
    if fjord==0:
        #Without a fjord
        H[:,:] = 0.0 + shelfslope*(xbig[:]-L_shelfbreak) + halfdelh*numpy.tanh((xbig[:]-L_shelfbreak)/(xwidthsb)**2)
        #Correct H with offset so the coast will be beachdepth
        offset = H[My][0]
        H[:,:] = -offset + beachdepth + H[:,:]
    else:
        #With a fjord
        #Leftcoast is an index function that defines the coast in the alongshore direction, i.e. draws a fjord.
        leftcoast = indent-1 - indent*ycurve
        leftcoast = int_(leftcoast)
        for j in range(My+2):
            num = leftcoast[j]
            H[j,num:] = 0.0 +  shelfslope*(xbig[num:]-L_shelfbreak) + halfdelh*tanh((xbig[num:]-L_shelfbreak)/(xwidthsb)**2)

            #Shift H to the beachdepth
            offset = H[j][num+1]
            H[j,:] = -offset + beachdepth + H[j,:]

            #Set H left of trough to a flat beachdepth
            H[:,:indent] = beachdepth
 
    #Making a copy of H untouched:
    Huntouched = H[:][:]

    #Trough array to be mapped onto H:
    trough = numpy.zeros((My+2,Mx+2))
    trough[:][:] = H_trough*numpy.tanh((xbig[:]-L_head)/(xwidthsb/1.75)**2)

    #Max of trough and H kept to map trough into H grid.
    trough = numpy.fmax(trough[:,:],H[:,:])

    #Splits the difference:
    difference = H[:][:]-trough[:][:]

    #Give the trough array its alongshore form:
    for j in range(My+2):
        difference[j][:] = difference[j][:]*ycurve[j]
    
    #Updates H by subtracting off trough
    H[:][:] = H[:][:]-difference[:][:]

    #Storage of offshore slope values for use if desired.
    slope = numpy.zeros(Mx)
    for i in range(Mx):
        slope[i] = (H[My][i+1]-H[My][i])/dx

    #Return H with and without trough
    return H,Huntouched

#Find shelfbreak, slope bottom indices:
def finddepthdistance(H,Mx,H_shelfbreak):
    for i in range(Mx+2):
        if H[0][i] >= H_shelfbreak:
            #Where it is deeper than the shelfbreak
            sbsb = i #sbsb = shelfbordershelfbreak
            break
    shelfbordershelfbreak = xbig[i]
    i=0
    for i in range(Mx+2):
        if H[0][i] >= 2*halfdelh:
            #Where it is deeper than the slope bottom
            sbbd = i  #sbbd = shelfbreak border deep
            break
    shelfbreakborderdeep = xbig[i]

    return sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section constructs winds:


#Wind magnitude: wind: windstress is tau (windscale [N/m^2]) divided by rho_0 (reference density [kg/m^3])
def makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd):   
    #What is the offshore span of the alongshore wind? (shelf, slope, shelf+slope, or shelf+slope+basin)
    if inflowregion == 1:
        wherewindleft = 0
        wherewindright = sbsb
    elif inflowregion == 2:
        wherewindleft = sbsb
        wherewindright = sbbd
    elif inflowregion == 3:
        wherewindleft = 0
        wherewindright = sbbd
    elif inflowregion == 4:
        wherewindleft = 0
        wherewindright = Mx
    
    #rho defined again, not sure if necessary.
    rho_0 = 1000. 
    
    #Create an alongshore wind matrix: My by Mx zero matrix, add on a My by windoffextent ones matrix, then multiply it all by the wind magnitude:
    
    windy = numpy.zeros((My,Mx))
    
    if True:
        #This method produces the offshore shape of the windy ones matrix with
        #a tanh, instead of previously used harsh edge 0 and 1's matrix.
        #It's a bit more clunky to do it this way, but smoother ouput.

        #Based off the trough's ycurve method seen above.
        def windcurve(windoffextent,windsteepness,zerox):    
            zeroxcurve = 1.0*(1-((1-0.5*(1-numpy.tanh((zerox-windoffextent/2)/windsteepness)))+(1-0.5*(1-numpy.tanh((-zerox-windoffextent/2)/windsteepness)))))
            return zeroxcurve

        #Place the wind according to inflow region.
        if inflowregion == 1:
            windcenter = 0
            windoffextent = int(shelfbordershelfbreak)*2
        elif inflowregion == 2:
            windcenter = -(int(shelfbordershelfbreak) + int(shelfbreakborderdeep-shelfbordershelfbreak)/2)
            windoffextent = int(shelfbreakborderdeep - shelfbordershelfbreak)
        elif inflowregion == 3:
            windcenter = 0
            windoffextent = int(shelfbreakborderdeep)*2
        elif inflowregion == 4:
            windcenter = 0
            windoffextent = int(Lx)*4
            
        #windsteepness controls either broad or sharp 1 to 0 curve.
        windsteepness = dx*5

        #Produce the curve.
        zeroxcurve = windcurve(windoffextent,windsteepness,zerox+Lx/2+windcenter)
        #Replace the windy matrix rows with this curve.
        for j in range(My):
            windy[j,:] = zeroxcurve[1:-1]
    
    if False:
        #0 to 1 hard edge matrix _|
        windoffextent = wherewindright-wherewindleft 
        windymatrix = ones((My,windoffextent)) 
        windy[::,wherewindleft:wherewindright] += windymatrix 

    #Scale the windy ones matrix by the wind scale 'wind.
    windy = windy*tauy_mag*flipForcing #was /rho_0
        
    #Creates an offshore wind matrix: only using alongshore wind currently.
    windx = numpy.zeros((My,Mx))

    #Initialize tau and tau gradient matrices.
    taux_y = numpy.zeros((My,Mx))
    tauy_x = numpy.zeros((My,Mx))
    tauy = numpy.zeros((My,Mx))
    taux = numpy.zeros((My,Mx))
    
    if False:
        #Total wind term in the governing equation is (tau^x,y)/(rho_0*H). Remember H is defined as a grid cell larger than Psi on all sides.

        #Shift the tauy=0 beyond the wherwindright up so that tauy
        #is more continuous, then smooth it out.

        #offsetloc is a bit less than wherewindright because that is where
        #the curve starts. Trying to smooth it flat without curve.
        offsetloc = wherewindright-100
        tauy[:,:offsetloc] = windy[:,:offsetloc]/rho_0/(H[1:-1,1:offsetloc+1])
        shift_tauy_value = tauy[My-1,offsetloc-1]
        tauy[:,offsetloc:] = shift_tauy_value

        #Calls upon the smoothing function to smooth this sharp transition.
        for j in range(My):
            #One row at a time.
            tauy_row = tauy[j,:]

            #Use 5 surrounding terms.
            convArray = array([1.0,2.0,3.0,2.0,1.0])
            convArray = convArray/sum(convArray)
            tauy_aug = array([0.]*2+list(tauy_row)+
                                [tauy_row[-1]]*2)
            smooth_tauy = convolve(convArray,tauy_aug,'valid')
            tauy[j,:] = smooth_tauy - smooth_tauy[0]

            #Smoothing messes up left terms, smoothing manually.
            tauy[j,1] = tauy[j,2] + abs(tauy[j,2]-tauy[j,3])
            tauy[j,0] = tauy[j,1] + abs(tauy[j,1]-tauy[j,2])
            tauy[j,:] = tauy[j,:] - tauy[j,0] + windy[j,0]

    if True:
        tauy = windy[:,:]/rho_0/(H[1:-1,1:-1])

    #Simply defined for taux for now because it's zero anyways.
    taux = windx[:,:]/rho_0/(H[1:-1,1:-1])

    #The following computes windstress gradients:
    
    #Compute offshore gradient of alongshore wind.
    for i in range(Mx):
        for j in range(My):
            if i == 0:
                tauy_x[j][i] = (tauy[j][i+1]-tauy[j][i])/(dx)
            elif i == Mx-1:
                tauy_x[j][i] = (tauy[j][i]-tauy[j][i-1])/(dx)
            else:
                tauy_x[j][i] = (tauy[j][i+1]-tauy[j][i-1])/(2*dx)
    
    #Compute alongshore gradient of offshore wind.
    for i in range(Mx):
        for j in range(My):
            if j == 0:
                taux_y[j][i] = (taux[j+1][i]-taux[j][i])/(dy)
            elif j == My-1:
                taux_y[j][i] = (taux[j][i]-taux[j-1][i])/(dy)
            else:
                taux_y[j][i] = (taux[j+1][i]-taux[j-1][i])/(2*dy)

    #Combine the components of windstress gradients into one term. The flipwindcurl term, defined at the beginning, is used in case sign is wrong.
    tau_term = (-taux_y[:][:] + tauy_x[:][:])*flipwindcurl

    #Shut the wind off downwave of where the trough starts, unless otherwise.
    if Add_wind == 0:
        #Depends on trough(s) or not.
        if troughamount == 0:
            windcutoff = upwavewindpos_notrough
            tau_term[0:windcutoff][:] = 0
        elif troughamount == 1:
            windcutoff = int((troughposition+L_trough)/dy)
            tau_term[0:windcutoff][:] = 0
        elif troughamount == 2:
            windcutoff = int((troughposition_two+L_trough)/dy)
            tau_term[0:windcutoff][:] = 0

    #Wind throughout the entire alongshore.
    else:
        if troughamount == 0:
            windcutoff = upwavewindpos_notrough
        elif troughamount == 1:
            windcutoff = int((troughposition+L_trough)/dy)
        elif troughamount == 2:
            windcutoff = int((troughposition_two+L_trough)/dy)

    #Winds have been made.

    return tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section constructs the four streamline (Psi) boundary conditions:

def makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd):

    #Coastal boundary set to constant (zero) for no normal flow.
    Psi_coast = 0.0

    #Initialize the upwave BC.
    Psi_upwave = numpy.zeros(Mx)
    
    #The following constructs the upwave (northern) inflow BC depending on the choice of inflow region and inflow type:
    #Define slopes of Psi in x to make Psi(x).
    if inflowtype == 1:
        #Constant transport structure.
        Psi_x = constant_transport
    elif inflowtype == 2:
        #Constant velocity structure.
        Psi_x = constant_velocity*H[My,1:-1]
    elif inflowtype == 3:
        #Alongshore wind-generated flow structure.
        Psi_x = tauy_mag/rho_0/r*H[My,1:-1]

    #Upwave inflow constrained to...?
    if inflowregion == 1:
        #Just on the shelf.
        startflow = 0
        endflow = sbsb
    elif inflowregion == 2:
        #Just on the slope.
        startflow = sbsb
        endflow = sbbd
    elif inflowregion == 3:
        #On the shelf and slope.
        startflow = 0
        endflow = sbbd
    elif inflowregion == 4:
    #Entire upwave boundary.
        startflow = 0
        endflow = -1

    #Define offset value to make sure Psi is zero at coast (y=mx+b).
    if inflowregion == 2 and inflowtype == 2:
        offset = Psi_x[sbsb]*xsmall[sbsb]
    elif inflowregion == 2 and inflowtype != 2:
        offset = Psi_x*xsmall[sbsb]
    else:
        offset = 0

    #Defining its offshore shape:
    Psi_upwave[:startflow] = Psi_coast
    if inflowtype == 2 or inflowtype == 3:
        Psi_upwave[startflow:endflow] = Psi_x[startflow:endflow]*xsmall[startflow:endflow] - offset
    else:
        Psi_upwave[startflow:endflow] = Psi_x*xsmall[startflow:endflow] - offset
    Psi_upwave[endflow:] = Psi_upwave[endflow-1]

    #Smoothing the shape for hard edges.
    if True:
        #Calls upon the smoothing function to smooth edges.
        convArray = numpy.array([1.0,2.0,3.0,2.0,1.0])
        convArray = convArray/numpy.sum(convArray)
        Psi_upwave_aug = numpy.array([0.]*2+list(Psi_upwave)+
                                [Psi_upwave[-1]]*2)
        smooth_Psi_upwave = numpy.convolve(convArray,Psi_upwave_aug,'valid')
        Psi_upwave = smooth_Psi_upwave-smooth_Psi_upwave[0]

    #Defines offshore to constant value for no flow (continuous).
    Psi_offshore = Psi_upwave[-1]

    #Matching donwave condition to upwave for simplicity.
    Psi_downwave = Psi_upwave

    return Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section is the main body of code solvin Psi:

#This function creates a unique k index from i,j in the x,y domain space. This is used to map i,j to a k space for Psi. The algorithm is semi-arbitrary but must remain consistent throughout all matrix construction and multiplication.
def kindex(i,j,Mx):
    #Take Psi(i,j) entries and stack into a Psi vector such that Psi rows are on top of each other in the vector.
    return Mx*j + i

#Workhorse function of the coding, filling matrices & vectors with
#boundary conditions & wind forcing, solves matrix equation for Psi:

#Emat*Psi = Forcing
def Psi_solver(inputPsi):

    #Unpacking inputPsi dictionary
    r = inputPsi['r']
    f = inputPsi['f']
    inflowtype = inputPsi['inflowtype']
    H = inputPsi['H']
    Psi_upwave = inputPsi['Psi_upwave']
    Psi_downwave = inputPsi['Psi_downwave']
    Psi_offshore = inputPsi['Psi_offshore']
    Psi_coast = inputPsi['Psi_coast']
    sbsb = inputPsi['sbsb']
    sbbd = inputPsi['sbbd']

    #Define the Psi vorticity function whose domain is [0:My-1,0:Mx-1]. Psi in domain space has shape Mx by My (i = __ by j = | where 0,0 is bot,left or south,west). Psi vector form will be 1_ by N| to contain all grid cells.
    Psi = numpy.zeros((My,Mx))
    Psi_vector = numpy.zeros((N,1))

    #The algoirthm for changing Psi space to Psi vector:
    for j in range(My-1):
        for i in range(Mx-1):
            #Don't change algorithm without changing all
            Psi_vector[kindex(i,j,Mx)] = Psi[j,i]

    #Creating the righthand side forcing vector (filled with wind if it exists):
    Forcing = numpy.zeros((N,1))
    for i in range(Mx):
        for j in range(My):
            #krow = kindex(i,j,Mx)
            #Forcing[krow,0] = tau_term[j][i]
            Forcing[kindex(i,j,Mx)] = tau_term[j,i]
            
    #Creating the lefthand side matrix:
    #Emat is size N by N [ ] (exists in k space).

    #Used below.
    krowlist=[]
    kcollist=[]
    Elist=[]

    #Fill Emat matrix, where if statements handle BCs:
    for i in range(Mx):
        for j in range(My):
            #The row of the Emat being modified is krow.
            krow = kindex(i,j,Mx)

            #In case of failure, notify the user.
            assert krow >= 0, 'oops'

            #-----------------------------------------------
            #E term never fails on a boundary (Psi_i,j).        
            kcol = kindex(i,j,Mx)
            assert kcol >= 0, 'oops'
            E = (2*(r/(dx**2))/((H[j+1][i+1])**2) + 2*(r/(dy**2))/((H[j+1][i+1])**2))
            #Update storage with these computes terms.
            kcollist.append(kcol)
            krowlist.append(krow)
            Elist.append(E)
            #End of E

            #----------------------------------------------

            #D term fails on the coast (Psi_i-1,j)           
            if i==0:
                D = (-(f/(4*dx*dy))*H[j+2][i+1]/((H[j+1][i+1])**2) + (f/(4*dx*dy))*H[j][i+1]/((H[j+1][i+1])**2) -
                 (r/(2*(dx**2)))*H[j+1][i+2]/((H[j+1][i+1])**3) + (r/(2*(dx**2)))*H[j+1][i]/((H[j+1][i+1])**3) -
                 (r/(dx**2))/((H[j+1][i+1])**2)) 
                Forcing[krow,0] = Forcing[krow,0] - Psi_coast*D 
            else:
                D = (-(f/(4*dx*dy))*H[j+2][i+1]/((H[j+1][i+1])**2) + (f/(4*dx*dy))*H[j][i+1]/((H[j+1][i+1])**2)
                - (r/(2*(dx**2)))*H[j+1][i+2]/((H[j+1][i+1])**3) + (r/(2*(dx**2)))*H[j+1][i]/((H[j+1][i+1])**3)
                - (r/(dx**2))/((H[j+1][i+1])**2))

                kcol = kindex(i-1,j,Mx)
                assert kcol >= 0, 'oops'
                kcollist.append(kcol)
                krowlist.append(krow)
                Elist.append(D)
            #End of D

            #------------------------------------------------
            
            #C term fails offshore (Psi_i+1,j)        
            if i == Mx-1:
                C = ((f/(4*dx*dy))*H[j+2][i+1]/((H[j+1][i+1])**2) - (f/(4*dx*dy))*H[j][i+1]/((H[j+1][i+1])**2)
                 + (r/(2*(dx**2)))*H[j+1][i+2]/((H[j+1][i+1])**3) - (r/(2*(dx**2)))*H[j+1][i]/((H[j+1][i+1])**3)
                 - (r/(dx**2))/((H[j+1][i+1])**2)) 
                Forcing[krow,0] = Forcing[krow,0] - Psi_offshore*C
            else:
                C = ((f/(4*dx*dy))*H[j+2][i+1]/((H[j+1][i+1])**2) - (f/(4*dx*dy))*H[j][i+1]/((H[j+1][i+1])**2)
                 + (r/(2*(dx**2)))*H[j+1][i+2]/((H[j+1][i+1])**3) - (r/(2*(dx**2)))*H[j+1][i]/((H[j+1][i+1])**3)
                 - (r/(dx**2))/((H[j+1][i+1])**2))

                kcol = kindex(i+1,j,Mx)
                assert kcol >= 0, 'oops'
                kcollist.append(kcol)
                krowlist.append(krow)
                Elist.append(C)
            #End of C

            #-------------------------------------------------------------
        
            #B term fails downwave (Psi_i,j-1)
            if j == 0:            
                B = ((f/(4*dx*dy))*H[j+1][i+2]/((H[j+1][i+1])**2) - (f/(4*dx*dy))*H[j+1][i]/((H[j+1][i+1])**2)
                 - (r/(2*(dy**2)))*H[j+2][i+1]/((H[j+1][i+1])**3) + (r/(2*(dy**2)))*H[j][i+1]/((H[j+1][i+1])**3)
                 - (r/(dy**2))/((H[j+1][i+1])**2)) 
                Forcing[krow,0] = Forcing[krow,0] - Psi_downwave[i]*B
            else:
                B = ((f/(4*dx*dy))*H[j+1][i+2]/((H[j+1][i+1])**2) - (f/(4*dx*dy))*H[j+1][i]/((H[j+1][i+1])**2)
                 - (r/(2*(dy**2)))*H[j+2][i+1]/((H[j+1][i+1])**3) + (r/(2*(dy**2)))*H[j][i+1]/((H[j+1][i+1])**3)
                 - (r/(dy**2))/((H[j+1][i+1])**2))
                
                kcol = kindex(i,j-1,Mx)
                assert kcol >= 0, 'oops'
                kcollist.append(kcol)
                krowlist.append(krow)
                Elist.append(B)
            #End of B

            #------------------------------------------------------------

            #A term fails upwave (Psi_i,j+1)
            if j == My-1:            
                A = (-(f/(4*dx*dy))*H[j+1][i+2]/((H[j+1][i+1])**2) + (f/(4*dx*dy))*H[j+1][i]/((H[j+1][i+1])**2)
                 + (r/(2*(dy**2)))*H[j+2][i+1]/((H[j+1][i+1])**3) - (r/(2*(dy**2)))*H[j][i+1]/((H[j+1][i+1])**3)
                 - (r/(dy**2))/((H[j+1][i+1])**2)) 
                Forcing[krow,0] = Forcing[krow,0] - Psi_upwave[i]*A
            else:
                A = (-(f/(4*dx*dy))*H[j+1][i+2]/((H[j+1][i+1])**2) + (f/(4*dx*dy))*H[j+1][i]/((H[j+1][i+1])**2)
                 + (r/(2*(dy**2)))*H[j+2][i+1]/((H[j+1][i+1])**3) - (r/(2*(dy**2)))*H[j][i+1]/((H[j+1][i+1])**3)
                 - (r/(dy**2))/((H[j+1][i+1])**2))
               
                kcol = kindex(i,j+1,Mx)
                assert kcol >= 0, 'oops'
                kcollist.append(kcol)
                krowlist.append(krow)
                Elist.append(A)
            #End of A

            #------------------------------------------------------------

    #Put lists into the Emat (in sparse form)
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    Emat = sparse.coo_matrix((Elist,(krowlist,kcollist)))
    sys.stdout.flush()

    #Changing to csc sparse:
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
   
    Emat = csc_matrix(Emat)
    sys.stdout.flush()

    #Solve matrix equation for Psi (linear algbre):
    Psi_vector = spsolve(Emat,Forcing)
    sys.stdout.flush()

    #Change Psi_vector, in k space, to Psi_matrix in domain space:
    Psi = Psi_vector.reshape(-1,Mx)

    return Psi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This section calls functions to run the code, depending on run type specified
#at the beginning, and writes data out:

#Initial choice of whichloop dictates array to run over, and corresponding file names.
if whichloop == 1:
    loop_array = H_shelfbreak_array
    fileName='H_shelfbreak_solver_output.sqlite'  
elif whichloop == 2:
    loop_array = H_trough_array
    fileName='H_trough_solver_output.sqlite'  
elif whichloop == 3:
    loop_array = L_wall_array
    fileName='L_wall_solver_output.sqlite'  
elif whichloop == 4:
    loop_array = L_trough_array
    fileName='L_trough_solver_output.sqlite'  
elif whichloop == 5:
    loop_array = L_head_array
    fileName='L_head_solver_output.sqlite'  
elif whichloop == 6:
    loop_array = L_shelfbreak_array
    fileName='L_shelfbreak_solver_output.sqlite'  
elif whichloop == 7:
    loop_array = r_array
    fileName='r_solver_output.sqlite'  
elif whichloop == 8:
    loop_array = f_array
    fileName='f_solver_output.sqlite'
elif whichloop == 9:
    loop_array = numpy.zeros(1)
    fileName='constant_solver_output.sqlite'
elif whichloop == 10:
    loop_array = troughshift_array
    fileName='troughshift_solver_output.sqlite'

#-----------------------------------------------------------------

#Building up the parameter lists.
print(" Starting the parameter list construction")
sys.stdout.flush()
Plist = []

#Leave this as true unless interested in changing the solver algorithm.
if True:
    for loop in range(loop_array.size):
        print("Building Plist, currently on loop number ",loop)
        sys.stdout.flush()

        #Depending on looping choice, the parameter lists being pasted to the solver will be modified to either have a parameters base value or its array(loop) value.

        #For each whichloop the function written in the above script will be
        #executed one at a time. Finally, all components necessary to pass the
        #solver will then be made. This is stored in inputPsi which can be
        #passed to the solver. 
        if whichloop == 1:
            #H_shelfbreak:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':loop_array[loop],'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            H_shelfbreak = loop_array[loop]
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H, 'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 2:
            #H_trough:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':loop_array[loop],'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 3:
            #L_wall:
            inputycurve = {'L_trough':L_trough,'L_wall':loop_array[loop],'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 4:
            #L_trough:
            inputycurve = {'L_trough':loop_array[loop],'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 5:
            #L_head:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':loop_array[loop]}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 6:
            #L_shelfbreak:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':loop_array[loop],'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 7:
            #r:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':loop_array[loop],'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}

        elif whichloop == 8:
            #f:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':loop_array[loop],'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}
            
        elif whichloop == 9:
            #constant:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'troughshift':troughshift,'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}
  
        elif whichloop == 10:
            #troughshift:
            inputycurve = {'L_trough':L_trough,'L_wall':L_wall,'zeroy':zeroy,'troughposition':troughposition,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two, 'troughshift':loop_array[loop],'Ly':Ly}
            ycurve = ycurvefn(inputycurve)
            inputBath = {'L_shelfbreak':L_shelfbreak,'H_shelfbreak':H_shelfbreak,'ycurve':ycurve,'H_trough':H_trough,'L_head':L_head}
            H,Huntouched = Bathymetry_Construction(inputBath)                   
            sbsb,shelfbordershelfbreak,sbbd,shelfbreakborderdeep = finddepthdistance(H,Mx,H_shelfbreak)
            tau_term,windy,windx,tauy,taux,tauy_x,taux_y,windcutoff = makewind(inflowregion,My,Mx,tauy_mag,flipForcing,sbsb,sbbd)
            Psi_upwave, Psi_downwave, Psi_offshore, Psi_coast = makePsiBC(constant_transport,constant_velocity,tauy_mag,H,sbsb,sbbd)
            inputPsi = {'r':r,'f':f,'inflowtype':inflowtype,'H':H,'Psi_upwave':Psi_upwave, 'Psi_downwave':Psi_downwave, 'Psi_offshore':Psi_offshore, 'Psi_coast':Psi_coast,'sbsb':sbsb,'sbbd':sbbd}
            
#------------------------------------------------------------------------------

        #Extra parameters written to sqlite files used in analysis.
        inputextra = {'whichloop':whichloop,'Ly':Ly,'Lx':Lx,'ysmall':ysmall,'xsmall':xsmall,'L_shelfbreak':L_shelfbreak,'Mx':Mx,'My':My,'dx':dx,'dy':dy,'troughamount':troughamount,'troughposition_one':troughposition_one,'troughposition_two':troughposition_two,'windcutoff':windcutoff}

        #Combine everything to be saved:
        save_this = {'inputextra':inputextra,'inputycurve':inputycurve,'ycurve':ycurve,'inputBath':inputBath,'H':H,'Huntouched':Huntouched,'inputPsi':inputPsi}

        #Package the current save_this list of the given loop into the Plist.
        Plist.append((pickle.dumps(save_this),))


#Initialize the solver run with the appropriate fileName and Plist.
print("Starting initialize database")
sys.stdout.flush()
conn = epr.initializeDB(fileName,Plist) 
print("Initialized run")
sys.stdout.flush()

#Data processing on multiple cores via enableparallelruns.
while True:
    #Obtain next parameter to be solved with. If empty, then all done and quit. If something to be solved, then use parameter and solve.
    paramPickle=epr.getNextParam(conn)
    if paramPickle==[]:
        print('ALL DONE WITH ALL parameters')
        break
    else:
        print("Starting run")
        sys.stdout.flush()
        
        #Obtain parameter list.
        save_this = pickle.loads(paramPickle[0])
        inputPsi = save_this['inputPsi']

        #Send the parameters to the solver function to compute Psi.
        Psi = Psi_solver(inputPsi)

        #Solver outputted answers.
        answer = {'Psi':Psi}

        #Finish up this particular run.
        epr.finishRun(conn,paramPickle,answer)
        print("Finishing run")
        print("")
        sys.stdout.flush()
