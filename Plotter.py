from __future__ import print_function
"""
###############################################################################

This program was created to plot the results produced by CEM, LAQGSM, and GSM.
The plotter will help to analyze data produced from these simulations, and is
meant to produce plots that can be placed into a report, either technical or
published.

#####  Edit Log  #####
----Name----     ----Date----      ----Description----
CMJ              May  25, 2017     File Creation
CMJ              May  25, 2017     Setup Input File format and reading
CMJ              May  30, 2017     Associated input with plot object
CMJ              May  31, 2017     Obtained data, finished project.
CMJ              June  1, 2017     Added auto-reformatting where necessary,
                                        made variable markers, changed marker
                                        size, and set so no fill exists.
CMJ              June  5, 2017     Added particle name to end of annotation, 
                                        minor fix for Python 2.7 use.
CMJ              June  8, 2017     Allowed for specific size to be set IN CODE,
                                        Fixed plot not showing
CMJ              June 12, 2017     Minor Modifications to marker sizes,
                                        markers used, allowed movable legend

###############################################################################
"""

# Imports
import matplotlib
# matplotlib.use('Agg') # Makes so plot doesn't display when desired
from pylab import rcParams
import matplotlib.pyplot as plt
import csv
import sys

plt.rc('font', **{'family': 'serif'})  # ,'sans-serif':['Helvetica']})
rcParams['figure.figsize'] = 6.811, 3.5 # Width, Height
# matplotlib.markers.MarkerStyle(marker=None, fillstyle=None)
markers = [u"v", u"8", u"d", u"*", u"p", u"h", u"+", u"s"]
defaultMarker = u"x"
dataMarker = u"v"#u"^"
fillIndx = 5 # Used in conjunction with below
markerFill = ['full', 'left', 'right', 'bottom', 'top', 'none']
mkrSize = 5 # Size of markers in plot


###############################################################################
# Classes
class PlotClass:
    __plotName = ""  # Title of Plot
    __pID = 0  # Particle ID
    __pName = ""  # No Particle Name Assigned
    __plotType = ""  # Type of Plot (dd, ei, ai)
    __simData = ""  # File w/ simulated data
    __dataFile = ""  # File w/ exp. data
    __plotAngle = 0  # Default plot angle
    __scaleX = 1  # Scale for Exp. data x values (energy or angle)
    __scaleY = 1  # Scale for Exp. data y values (cross sections)
    __annotateIt = ""  # Annotation
    __annotateX = 0.05  # Use set positioning for annotation (0<x<1)
    __annotateY = 0.9  # Use set positioning for annotation (0<y<1)
    __save = ""  # File not saved as default
    __xLabel = ""  # X Label
    __yLabel = ""  # Y Label
    __reformat = 0  # Reformat data file if requested
    __showPlot = True # Show plot to user
    __xMax = 250 # Upper x axis limit
    __scalerX = 1 # Scales simulated data in x - for testing only
    __scalerY = 1 # Scales simulated data in y - for testing only
#    options for legend locations = {'best': 0, 'upper right': 1, 'upper left':
    #  2, 'lower left': 3, 'lower right': 4, 'right': 5, 'center left': 6,
    # 'center right': 7, 'lower center': 8, 'upper center': 9, 'center': 10}
    __locX = 0 # Legend X location
    __locY = 0 # Legend Y location

    ############################################################################
    # Getters, Setters
    ############################################################################
    # --------------------------------------------------------------------------
    def setName(self, plotName):
        if testing > 1 or verbose > 1:
            print('Plot Name is \"%s\".' % (item))
        self.__plotName = plotName

    def getName(self):
        return self.__plotName

    # --------------------------------------------------------------------------
    def setPID(self, pID):
        self.__pID = pID

    def getPID(self):
        return self.__pID

    # --------------------------------------------------------------------------
    def setXMax(self, xMax):
        self.__xMax = xMax

    def getXMax(self):
        return self.__xMax

    # --------------------------------------------------------------------------
    def setShowPlot(self, showPlot):
        self.__showPlot = showPlot

    def getShowPlot(self):
        return self.__showPlot

    # --------------------------------------------------------------------------
    def setPName(self, pName):
        self.__pName = pName

    def getPName(self):
        return self.__pName

    # --------------------------------------------------------------------------
    def setPlotType(self, plotType):
        self.__plotType = plotType

    def getPlotType(self):
        return self.__plotType

    # --------------------------------------------------------------------------
    def setSimData(self, simData):
        self.__simData = simData

    def getSimData(self):
        return self.__simData

    # --------------------------------------------------------------------------
    def setExpData(self, dataFile):
        self.__dataFile = dataFile

    def getExpData(self):
        return self.__dataFile

    # --------------------------------------------------------------------------
    def setPlotAngle(self, plotAngle):
        self.__plotAngle = plotAngle

    def getPlotAngle(self):
        return self.__plotAngle

    # --------------------------------------------------------------------------
    def setScaleX(self, scaleX):
        self.__scaleX = scaleX

    def getScaleX(self):
        return self.__scaleX

    # --------------------------------------------------------------------------
    def setScaleY(self, scaleY):
        self.__scaleY = scaleY

    def getScaleY(self):
        return self.__scaleY

    # --------------------------------------------------------------------------
    def setAnnotateIt(self, annotateIt):
        self.__annotateIt = annotateIt

    def getAnnotateIt(self):
        return self.__annotateIt

    # --------------------------------------------------------------------------
    def setAnnotateX(self, annotateX):
        self.__annotateX = annotateX

    def getAnnotateX(self):
        return self.__annotateX

    # --------------------------------------------------------------------------
    def setAnnotateY(self, annotateY):
        self.__annotateY = annotateY

    def getAnnotateY(self):
        return self.__annotateY

    # --------------------------------------------------------------------------
    def setSave(self, save):
        self.__save = save

    def getSave(self):
        return self.__save

    # --------------------------------------------------------------------------
    def setReformat(self, reformat):
        self.__reformat = reformat

    def getReformat(self):
        return self.__reformat

    # --------------------------------------------------------------------------
    def setLocX(self, locX):
        self.__locX = locX

    def getLocX(self):
        return self.__locX

    # --------------------------------------------------------------------------
    def setLocY(self, locY):
        self.__locY = locY

    def getLocY(self):
        return self.__locY

    # --------------------------------------------------------------------------
    def checkErrors(self):
        # Ensuring all necessary plot properties are set
        if self.__pID == 0:
            print('ERROR: Particle not specified; cannot create plot.')
            print('Use \" plot (Part. Name)\" to specify particle.')
            print("Warning Detected, PROGRAM Exiting...")
            return 1
        elif self.__plotType == "":
            print('ERROR: Plot Type not specified; cannot create plot.')
            print('Use \" plot (dd/ei/ai)\" to specify plot type.')
            print("Warning Detected, PROGRAM Exiting...")
            return 1
        elif self.__simData == "":
            print('ERROR: Simulated Data File not specified; cannot create '
                  'plot.')
            print('Use \" read (Sim. Data Filename)\" to specify this file.')
            print("Warning Detected, PROGRAM Exiting...")
            return 1
        elif self.__plotAngle == 0 and self.__plotType == 'dd':
            print('ERROR: Angle not specified for double differential plot '
                  'type; cannot create plot.')
            print('Use \" angle (Angle)\" to specify plot angle.')
            print("Warning Detected, PROGRAM Exiting...")
            return 1
        else:
            if verbose > 0:
                print('All necessary plot specifications received.')
            return 0

    # --------------------------------------------------------------------------
    def setPlotDetails(self):  # Setting title and axis labels
        myFont = 14
        plt.xlabel(self.__xLabel, fontsize=myFont)
        plt.ylabel(self.__yLabel, fontsize=myFont)
        plt.title(self.__plotName, fontsize=myFont + 2)

        # Minor and Major GridLines
        plt.grid(b=True, which='both', color='0.9', linestyle=':')
        if not self.__annotateIt == '':
            # Adding particle name to annotation
            self.__annotateIt = self.__annotateIt.strip() + " " + self.__pName

            # Placing annotation onto graph
            plt.annotate(self.__annotateIt, xy=(self.__annotateX,
                                                self.__annotateY),
                         xycoords='axes fraction', horizontalalignment='left',
                         verticalalignment='top', fontsize=myFont)

        # Axis limits
        if self.__xMax > 0:
            if self.__plotType == 'ei':
                if self.__xMax == 250:
                    # Value unchanged from default; setting to 200
                    self.__xMax = 200
                plt.xlim((0, self.__xMax))
            else:
                plt.xlim((0, self.__xMax))
        else:
            # Cannot create plot, unphysical x values plotted
            print("ERROR: A positive xMax value must be specified (xMax = "
                  "%/1f." %(self.__xMax))
            print("Error detected, program stopping...")
            sys.exit()

        # Informing user if the simulated data is being scaled
        if not self.__scalerX == 1 or not self.__scalerY == 1:
            print("Simulated data is being scaled.")
            print("\t Y-Scale: %.1f" %(self.__scalerY))
            print("\t X-Scale: %.1f" %(self.__scalerX))


    # --------------------------------------------------------------------------
    def createName(self):
        # Obtaining plot type name
        if self.__plotType == 'dd':
            temp = "Double Differential"
        elif self.__plotType == 'ei':
            temp = "Energy Integrated"
        elif self.__plotType == 'ai':
            temp = "Angle Integrated"

        # Creating name
        self.__plotName = temp + " Spectra (" + self.__pName

        # Adding Angle if appropriate
        if self.__plotType == 'dd':
            self.__plotName = self.__plotName + " at " + \
                              str(self.__plotAngle) + " degrees)"
        else:
            self.__plotName = self.__plotName + ")"

        if verbose > 1:
            print('\"%s\" assigned to plot title.' % (self.__plotName))

    # --------------------------------------------------------------------------
    def getAxis(self):
        self.__xLabel = ""  # Empty String

        # Use particle name for x-axis label only
        if self.__pName.lower() == 'n'.lower():
            # Neutron
            partName = "Neutron"
        elif self.__pName.lower() == "p".lower():
            # Proton
            partName = "Proton"
        elif self.__pName.lower() == "d".lower():
            # Deuteron
            partName = "Deuteron"
        elif self.__pName.lower() == "t".lower():
            # Tritium
            partName = "Tritium"
        else:
            # Use particle name given already
            partName = self.__pName

        self.__yLabel = self.__yLabel + 'Cross Section'
        if self.__plotType == 'dd':
            self.__yLabel = self.__yLabel + ' [mb/MeV/sr]'
            self.__xLabel = self.__xLabel + partName + ' Energy [MeV]'
        elif self.__plotType == 'ai':
            self.__yLabel = self.__yLabel + ' [mb/MeV]'
            self.__xLabel = self.__xLabel + partName + ' Energy [MeV]'
        elif self.__plotType == 'ei':
            self.__yLabel = self.__yLabel + ' [mb/sr]'
            self.__xLabel = self.__xLabel + partName + ' Angle [Degrees]'

        if verbose > 1:
            print('X Label (in LaTeX) is \"%s\"' % self.__xLabel)
            print('Y Label (in LaTeX) is \"%s\"' % self.__yLabel)

    # --------------------------------------------------------------------------
    def getSimulatedData(self):

        # Setting so no data exists
        dataExists = 0
        # ------------------------------------------------------------------
        # Opening File
        simFile = open(self.__simData, 'r')
        printlines(1)
        print('Reading file \"%s\".' % simFile.name)

        # ------------------------------------------------------------------
        # Reading File
        labelIndx = 0
        labels = [""]
        simFile.seek(0)
        temp = simFile.readline()

        # ------------------------------------------------------------------
        # Obtaining indices
        indx1 = temp.find(':', 0, len(temp)) + 1
        indx2 = temp.find(';', indx1, -1)

        # Label for first data set
        labels[labelIndx] = temp[indx1:indx2].strip()

        # Iterating over all sets
        # ------------------------------------------------------------------
        while not indx2 == -1:

            # Obtaining Indices
            labelIndx += 1
            indx1 = indx2 + 1
            indx2 = temp.find(';', indx1, len(temp))

            # Saving label from file
            labels.append(temp[indx1:indx2].strip())

            # Fail safe
            if labelIndx > 10:
                break

        # For user output
        if verbose > 1:
            for indx in range(0, labelIndx + 1, 1):
                print('\t Label for Data Set %d: \"%s\"' % (indx + 1,
                                                            labels[indx]))

        # ------------------------------------------------------------------
        # Setting up Arrays
        partID = []
        energyDat = []
        angleDat = []
        indx = 0
        dataExists = 0

        # 2-D Array to store simulated data
        dataX = [[] for i in range(file_len(self.__simData))]

        # ------------------------------------------------------------------
        # Reading the numbers from the file
        plots = csv.reader(simFile, delimiter=',')
        for row in plots:
            partID.append(int(row[0]))
            energyDat.append(self.__scalerX * int(row[1]))
            angleDat.append(self.__scalerX * int(row[2]))

            # Assigning data from simulations to 2-D array
            # --------------------------------------------------------------
            for tempIndx in range(0, labelIndx + 1, 1):
                temp = float(row[tempIndx + 3])
                dataX[tempIndx].append(self.__scalerY * temp)

            # Printing for user if desired
            if verbose > 2:
                print('PID: %d \t E: %.1f \t A: %.1f \t Data =' %(partID[indx],
                                    energyDat[indx], angleDat[indx]), end='')
                for tempIndx in range(0, labelIndx + 1, 1):
                    print(' %.3e,' % (dataX[tempIndx][indx]), end='')
                print('', end='\n')
                indx += 1

        # ------------------------------------------------------------------
        # Closing File
        simFile.close()
        print('File \"%s\" closed.' % (self.__simData))

        # ------------------------------------------------------------------
        # Sorting data based on specifications
        print('Sorting Data...')
        if self.__plotType == 'dd' or self.__plotType == 'ai':

            # Telling program location of angle integrated data
            if self.__plotType == 'ai':
                self.__plotAngle = 361

            # --------------------------------------------------------------
            # Obtaining data as 1-D array and plotting
            for indx2 in range(0, labelIndx + 1, 1):

                # Setting up arrays
                xVal = []
                yVal = []

                # Sorting through stored data
                # ----------------------------------------------------------
                for indx1 in range(0, len(partID), 1):
                    if partID[indx1] == self.__pID:
                        # Reading Desired particle data, sort through angles
                        if angleDat[indx1] == self.__plotAngle:
                            # Checking if non-zero, omitting E=5001 Values
                            if not dataX[indx2][indx1] == 0 and not \
                                            energyDat[indx1] > 5000:
                                xVal.append(energyDat[indx1])
                                yVal.append(dataX[indx2][indx1])

                # Counting number of sets that contain data
                plotSet = False
                print('\t%d Data points obtained from \"%s\".' % (len(xVal),
                                                                labels[indx2]))
                if len(xVal) > 0:
                    dataExists += 1
                    plotSet = True

                # ----------------------------------------------------------
                # Displaying data for individual plot
                if verbose > 2:
                    printlines(1)
                    print('\"%s\" Data:' % labels[indx2])
                    for i in range(0, len(xVal), 1):
                        print('%.1f \t %.3e' % (xVal[i], yVal[i]))
                    print('%d Data Points for \"%s\" data.'
                          % (len(xVal), labels[indx2]))

                # ----------------------------------------------------------
                # Plotting Data if it exists
                if plotSet:
                    if dataExists < len(markers):
                        if dataExists % 2 == 0:
                            tempSize = mkrSize + 1
                        else:
                            tempSize = mkrSize
                        plt.semilogy(xVal, yVal, label=labels[indx2],
                                     linestyle='None',
                                     marker=markers[dataExists],
                                     mfc=markerFill[fillIndx],
                                     markersize=tempSize)
                    else:
                        if dataExists % 2 == 0:
                            tempSize = mkrSize + 1
                        else:
                            tempSize = mkrSize
                        plt.semilogy(xVal, yVal, label=labels[indx2],
                                     linestyle='None',
                                     marker=defaultMarker,
                                     mfc=markerFill[fillIndx],
                                     markersize=tempSize)

        elif self.__plotType == 'ei':

            # --------------------------------------------------------------
            # Obtaining data as 1-D array and plotting
            for indx2 in range(0, labelIndx + 1, 1):

                # Setting up arrays
                xVal = []
                yVal = []

                # Sorting through stored data
                # ----------------------------------------------------------
                for indx1 in range(0, len(partID), 1):
                    if partID[indx1] == self.__pID:
                        # Reading Desired particle data, sort through energy
                        if energyDat[indx1] == 5001:
                            # Checking if non-zero, omitting A=361 Values
                            if not dataX[indx2][indx1] == 0 and not \
                                            angleDat[indx1] > 361:
                                xVal.append(angleDat[indx1])
                                yVal.append(dataX[indx2][indx1])

                # Counting number of sets that contain data
                plotSet = False
                if len(xVal) > 0:
                    dataExists += 1
                    plotSet = True

                # ----------------------------------------------------------
                # Displaying data for individual plot
                if verbose > 2:
                    printlines(1)
                    print('\"%s\" Data:' % labels[indx2])
                    for i in range(0, len(xVal), 1):
                        print('%.1f \t %.3e' % (xVal[i], yVal[i]))
                    print('%d Data Points for \"%s\" data.'
                          % (len(xVal), labels[indx2]))

                # ----------------------------------------------------------
                # Plotting Data if it exists
                if plotSet:
                    if dataExists < len(markers):
                        if dataExists % 2 == 0:
                            tempSize = mkrSize + 1
                        else:
                            tempSize = mkrSize
                        plt.semilogy(xVal, yVal, label=labels[indx2],
                                     linestyle='None',
                                     marker=markers[dataExists],
                                     mfc=markerFill[fillIndx],
                                     markersize=tempSize)
                    else:
                        if dataExists % 2 == 0:
                            tempSize = mkrSize + 1
                        else:
                            tempSize = mkrSize
                        plt.semilogy(xVal, yVal, label=labels[indx2],
                                     linestyle='None',
                                     marker=defaultMarker,
                                     mfc=markerFill[fillIndx],
                                     markersize=tempSize)

        return dataExists

    # --------------------------------------------------------------------------
    def createExpData(self):

        if self.__reformat == 1:
            print('Reformatting \"%s\"...' % (self.__dataFile))
            reformat(self.__dataFile)

        # Opening re-formatted file
        expFile = open(self.__dataFile, 'r')
        printlines(1)
        print('Reading file \"%s\".' % expFile.name)

        # ------------------------------------------------------------------
        # Reading File
        labelIndx = 0
        labels = [""]
        temp = expFile.readline()

        # ------------------------------------------------------------------
        # Obtaining indices
        indx1 = 0
        indx2 = temp.find(',', indx1, -1)

        if indx2 == -1:

            # User forgot to reformat file
            myPlot.setReformat(1)

            # Starting this routine over, then returning to program
            myPlot.createExpData()
            return

        # Label for first data set
        labels[labelIndx] = temp[indx1:indx2].strip().lower()

        # Iterating over all sets
        # ------------------------------------------------------------------
        while not indx2 == -1:

            # Obtaining Indices
            labelIndx += 1
            indx1 = indx2 + 1
            indx2 = temp.find(',', indx1, len(temp))

            # Saving label from file
            labels.append(temp[indx1:indx2].strip().lower())

            # Fail safe
            if labelIndx > 10:
                break

        # For user output
        if verbose > 1:
            for indx in range(0, labelIndx + 1, 1):
                print(
                    '\t Label for Column %d: \"%s\"' % (indx + 1, labels[indx]))
        # ------------------------------------------------------------------
        # 2-D Array to store simulated data
        expData = [[] for i in range(file_len(self.__dataFile))]
        numRow = 0  # Number of rows

        # User printout - header for a table (table created below)
        if verbose > 2:
            print('Data Read is:')
            print('%s' % (labels[0]), end='')
            for i in range(1, labelIndx + 1, 1):
                print('\t %s' % (labels[i]), end='')
            print('')

        # ------------------------------------------------------------------
        # Obtaining X, Y Values for plot
        plots = csv.reader(expFile, delimiter=',')

        for row in plots:

            # --------------------------------------------------------------
            # Assigning data from simulations to 2-D array
            for tempIndx in range(0, labelIndx + 1, 1):
                temp = float(row[tempIndx])
                expData[tempIndx].append(temp)

            if verbose > 2:
                for tempIndx in range(0, labelIndx + 1, 1):
                    print('%.2f \t' % (expData[tempIndx][numRow]), end='')
                print('')

            # Incrementing row No.
            numRow += 1

        # ------------------------------------------------------------------
        # Closing File
        expFile.close()
        print('File \"%s\" closed.' % (self.__dataFile))

        # ------------------------------------------------------------------
        # Sorting through data -assigning values based on headers
        print('Sorting Exp. data...')

        # Setting up arrays
        angleDat = []
        energyDat = []
        CSDat = []

        # --------------------------------------------------------------
        # Obtaining Data
        for indx1 in range(0, len(labels) + 1, 1):
            if verbose > 1:
                print('Obtaining data from \"%s\" column' % (labels[indx1]))

            for indx2 in range(0, len(expData[indx1]), 1):

                # Each label is either "Energy", "Angle", or "CS".
                # These must be present based on the plot produced
                if labels[indx1] == 'energy'.lower():
                    # Store energy values
                    energyDat.append(expData[indx1][indx2])

                if labels[indx1] == 'angle'.lower():
                    # Store angle values
                    angleDat.append(expData[indx1][indx2])

                if labels[indx1] == 'cs'.lower():
                    # Store Cross Section Values
                    CSDat.append(expData[indx1][indx2])

        # --------------------------------------------------------------
        # Need an angle specification at which to collect data for dd, ei plots
        if len(angleDat) == 0 and self.__plotType == 'dd':
            # No Angle information given, assuming data at correct angle
            if verbose > 0:
                print('Assuming Data is produced at a %.1f degree angle.' % (
                    self.__plotAngle))
            angleDat = [self.__plotAngle] * len(CSDat)

        # Setup arrays
        xDat = []
        yDat = []

        # Error protection for CS Data
        if len(CSDat) == 0:
            # Error - need data to plot, perhaps issue with headers.
            print("ERROR: Missing cross section data in data file.  "
                  "Issue")
            print("may be due to lack of table headers.  Please correct.")
            print("Warning Detected, PROGRAM Exiting...")
            sys.exit()

        if len(energyDat) == 0 and not self.__plotType == 'ei':
            # Error - need data to plot, perhaps issue with headers.
            print("ERROR: Missing energy data in data file.  Issue "
                  "may be due to lack of table headers.  Please correct.")
            print("Warning Detected, PROGRAM Exiting...")
            sys.exit()

        # --------------------------------------------------------------
        # Get X and Y data for the plot at hand
        if self.__plotType == 'dd'.lower():

            # Angle given, plot Energy as X, CS as Y at given angle
            for indx1 in range(0, len(CSDat), 1):
                if self.__plotAngle == angleDat[indx1]:
                    # At correct angle, obtain X, Y values
                    xDat.append(energyDat[indx1])
                    yDat.append(CSDat[indx1])

        elif self.__plotType == 'ai'.lower():

            # Angle not necessary at all, take all energy and CS data
            for indx1 in range(0, len(CSDat), 1):
                # Assign given values to text
                xDat.append(energyDat[indx1])
                yDat.append(CSDat[indx1])

        elif self.__plotType == 'ei'.lower():

            # Error protection
            if len(angleDat) == 0:
                # ERROR, need angle information (perhaps missing label)
                print("ERROR: Missing angle data in data file.  Issue")
                print("may be due to lack of table headers.  Please correct.")
                print("Warning Detected, PROGRAM Exiting...")
                sys.exit()
            else:
                # Angle becomes x-axis, energy data ignored
                for indx1 in range(0, len(CSDat), 1):
                    # Assign given values to text
                    xDat.append(angleDat[indx1])
                    yDat.append(CSDat[indx1])

        # Print data points read
        print('\t %d data points read from \"%s\".' % (len(yDat),
                                                       self.__dataFile))

        # Scaling data
        if not self.__scaleX == 1:
            # Scaling data
            print('Scaling X-Axis data...')
            for indx in range(0, len(xDat), 1):
                xDat[indx] = xDat[indx] * self.__scaleX

        if not self.__scaleY == 1:
            # Scaling data
            print('Scaling Y-Axis data...')
            for indx in range(0, len(yDat), 1):
                yDat[indx] = yDat[indx] * self.__scaleY

        # Adding data to plot
        plt.semilogy(xDat, yDat, label='Exp. Data',
                     linestyle='None',
                     marker=dataMarker,
                     mfc=markerFill[fillIndx],
                     markersize=mkrSize+1.25,
                     color='0.35')

    # --------------------------------------------------------------------------
    def createPlot(self):

        # Blank line for readability in terminal
        printlines(1)
        print('Gathering Plot Informaion...')

        # Checking for minimum specifications
        error = myPlot.checkErrors()
        if error == 1:
            return 1

        # ----------------------------------------------------------------------
        # Reading input file
        if not self.__simData == "":
            dataExists = myPlot.getSimulatedData()

        # ----------------------------------------------------------------------
        # Obtaining Experimental Data
        if not self.__dataFile == "":
            myPlot.createExpData()

        # ----------------------------------------------------------------------
        # Creating plot properties if there is stuff that can be plotted
        if dataExists > 0:

            # Setting plot propterites
            if verbose > 0:
                print('%d sets of simulated data exists.' % dataExists)

            # ------------------------------------------------------------------
            #  Creating plot name if not specified
            if self.__plotName == "":
                myPlot.createName()

            # ------------------------------------------------------------------
            # X,Y label
            myPlot.getAxis()

            # ------------------------------------------------------------------
            # Creating plot properties
            myPlot.setPlotDetails()

            # ------------------------------------------------------------------
            # Creating plot legend
            if self.__locX == 0 and self.__locY == 0:
                # Using defaults
                plt.legend(loc = 0)
            else:
                plt.legend(loc = (self.__locX, self.__locY))

            # ------------------------------------------------------------------
            # Displaying plot
            displayPlot(self.__save, self.__showPlot)

        # ------------------------------------------------------------------
        return 0  # No Errors
        # ------------------------------------------------------------------

###############################################################################
# Functions

# Returns length of file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1


# Prints predefined number of lines for user output
def printlines(num):
    print('\n' * num)


# Displays plot
def displayPlot(figName, showPlot):
    # Whitespace
    printlines(1)

    # Saving plot if desired
    if not figName == "":
        # Save figure
        print('Saving figure as \"%s\".' % (figName))
        plt.savefig(figName, bbox_inches='tight', pad_inches=0.15, dpi=dpiVal)
        print('\tFigure saved.')

    # Displaying plot
    if showPlot:
        print('Displaying Plot...')
        plt.show()

    # Clearing figure of all information
    plt.clf()


# Reformats file with a lot of unneccesary spaces
def reformat(fname):
    # Reading file and removing extra spaces, replacing with ','
    lines = open(fname).readlines()
    for i in range(0, len(lines), 1):
        lines[i] = ",\t".join(lines[i].split())

    # Re-writing data to the same file
    tempFile = open(fname, 'w')
    for indx in range(0, len(lines), 1):
        tempFile.writelines('%s\n' % lines[indx])


###############################################################################
# Defaults
getDetails = 0  # Start programming by getting plot details
plotNum = 1  # 1 Plot is default, value iterated as more defined
linesUsed = 0  # Number of lines used to define a plot
verbose = 1  # Provides more output information to user
myPlot = PlotClass()  # Initialize plot object
dpiVal = 600  # Default dpi value
testing = 0 # For testing the example script quickly; provides additional
# output (like verbose 1, with some verbose 2)

###############################################################################
# Legal Notice
printlines(2)
print('--------------------------------------------------------')
print('This program was created by CMJ for work funded by the\n'
      'Los Alamos National Laboratory. Contact CMJ for questions\n'
      'regarding this work and its ownership at <chasemjun@lanl.gov>,\n'
      'or at <junechas@isu.edu>.  This tool was created during the\n'
      'LAQGSM deprecation, or GSM, project.')
print('--------------------------------------------------------')
printlines(1)

###############################################################################
# Obtaining and Opening Input File
prompt = 'Enter the Input File Name (with extension): '
if testing == 0:
    # temp = input(prompt)
    temp = sys.argv[1] # Forces user to pass in file name
else:
    temp = 'sample.inp'
    print('%s %s' % (prompt, temp))

inputFile = open(temp, 'r')
print('Reading file \"%s\".' % inputFile.name)
printlines(1)

###############################################################################
# Reading input file
print('Plot %d Details:' % plotNum)
print('-------------------')
while getDetails == 0:

    # Reading Line
    temp = inputFile.readline()
    linesUsed += 1  # Iterating number of lines used to define plot

    if temp.strip().lower() == 'end plot':
        # Plot Definition Ended

        # --------------------------------------------------------------
        # Create plot here
        # -------------
        error = myPlot.createPlot()
        if not error == 0:
            # Stopping program if signaled by createPlot() function
            if verbose > 0:
                print('STOP SIGNAL RECEIVED: Stopping program.')
            break
        # ----------------------------------------------------------------

        # Resetting variables Used
        plotNum += 1  # Defining a different plot or End Of File
        linesUsed = 0  # Reset number of lines used to define plot
        del myPlot  # Remove plot object

        temp = inputFile.readline()  # Reading empty line between plot def.

        # Error protection
        if temp == '':
            # End of File Reached
            break
        elif not temp.strip() == '':
            # Something other than a new line was read
            print("ERROR: Input error for plot %d, line %d: Missing Blank "
                  "Line - \"%s\"." % (plotNum, linesUsed, temp))
            print("Warning Detected, PROGRAM Exiting...")
            break

        # Create new plot object; End of File Signal and Error not received
        myPlot = PlotClass()  # Initialize plot object
        dpiVal = 600  # Default dpi value
        printlines(1)
        print('Plot %d Information:' % plotNum)
        print('-------------------')

    else:
        # Valid Line, need to now get plot specs

        # Obtaining index of space
        spaceIndx = temp.find(' ', 0, len(temp))

        # Error protection
        if temp == '\n':
            print("ERROR: Input error for plot %d, line %d: Extra Blank "
                  "Line - \"%s\"." % (plotNum, linesUsed, temp))
            print("Warning Detected, PROGRAM Exiting...")
            break
        elif spaceIndx == -1:
            print("ERROR: Input error for plot %d, line %d: Missing Space "
                  "Line - \"%s\"." % (plotNum, linesUsed, temp))
            print("Warning Detected, PROGRAM Exiting...")
            break

        # Getting flag and item on line, printing if testing
        flag = temp[:spaceIndx].lower()  # Flag NOT case sensitive
        item = temp[spaceIndx + 1:-1].strip()

        if testing > 0 or verbose > 0:
            if testing > 2:
                # Prints location in input along with flag/item combo
                print('Plot %d, Line %d: ' % (plotNum, linesUsed), end='')

            print('%s' % temp, end='')  # Echos the input

            if testing > 1:
                # Prints the flag/item combo
                print(' breaks down to \"%s\", \"%s\".' % (flag, item))

        # Setting plot properties
        # ---------------------------------------------------------------------
        if flag == 'particle':

            # Removing case sensitivity from input
            item = item.lower()

            # Getting pID
            if item == 'n'.lower():
                pID = 1
                myPlot.setPName('n')
            elif item == 'p'.lower():
                pID = 2
                myPlot.setPName('p')
            elif item == 'd'.lower():
                pID = 3
                myPlot.setPName('d')
            elif item == 't'.lower():
                pID = 4
                myPlot.setPName('t')
            elif item == 'He3'.lower():
                pID = 5
                myPlot.setPName('$^{3}He$')
            elif item == 'He4'.lower():
                pID = 6
                myPlot.setPName('$^{4}He$')
            elif item == 'He6'.lower():
                pID = 7
                myPlot.setPName('$^{6}He$')
            elif item == 'Li6'.lower():
                pID = 8
                myPlot.setPName('$^{6}Li$')
            elif item == 'Li7'.lower():
                pID = 9
                myPlot.setPName('$^{7}Li$')
            elif item == 'Li8'.lower():
                pID = 10
                myPlot.setPName('$^{8}Li$')
            elif item == 'Li9'.lower():
                pID = 11
                myPlot.setPName('$^{9}Li$')
            elif item == 'Be7'.lower():
                pID = 12
                myPlot.setPName('$^{7}Be$')
            elif item == 'Be9'.lower():
                pID = 13
                myPlot.setPName('$^{9}Be$')
            elif item == 'Be10'.lower():
                pID = 14
                myPlot.setPName('$^{10}Be$')
            elif item == 'B10'.lower():
                pID = 15
                myPlot.setPName('$^{10}B$')
            elif item == 'B11'.lower():
                pID = 16
                myPlot.setPName('$^{11}B$')
            elif item == 'B12'.lower():
                pID = 17
                myPlot.setPName('$^{12}B$')
            elif item == 'C11'.lower():
                pID = 18
                myPlot.setPName('$^{11}C$')
            elif item == 'C12'.lower():
                pID = 19
                myPlot.setPName('$^{12}C$')
            elif item == 'C13'.lower():
                pID = 20
                myPlot.setPName('$^{13}C$')
            elif item == 'C14'.lower():
                pID = 21
                myPlot.setPName('$^{14}C$')
            else:
                # Error - invalid item
                print("ERROR: Input error for plot %d, line %d: Invalid "
                      "Particle Item - \"%s\"." % (plotNum, linesUsed, item))
                print("Warning Detected, PROGRAM Exiting...")
                break
                
            # Assigning pID
            myPlot.setPID(pID)

            if testing > 1:
                print('%s (%s) assigned PID of %d.' % (flag, item,
                                                       myPlot.getPID()))

        # ---------------------------------------------------------------------
        elif flag == 'plot'.lower():

            # Converting to all lower case
            item = item.lower()

            # Checking if valid item read
            if item == 'dd'.lower() or item == 'ai'.lower() or \
                            item == 'ei'.lower():
                # Differential Cross Sections
                myPlot.setPlotType(item)
            else:
                # Error - invalid item
                print("ERROR: Input error for plot %d, line %d: Invalid "
                      "Plot Type Item - \"%s\"." % (plotNum, linesUsed, item))
                print("Warning Detected, PROGRAM Exiting...")
                break

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getPlotType()))

        # ---------------------------------------------------------------------
        elif flag == 'read'.lower():

            # Simulated Data is from this file
            myPlot.setSimData(item)

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getSimData()))

        # ---------------------------------------------------------------------
        elif flag == 'data'.lower():

            # Exp. Data associated with simulated data
            myPlot.setExpData(item)

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getExpData()))

        # ---------------------------------------------------------------------
        elif flag == 'angle'.lower():

            # Angle at which to plot
            plotAngle = int(item)
            if plotAngle > 360 or plotAngle < 0:
                print('ERROR: Angle value must be between 0 and 360 (degrees).')
                print("Warning Detected, PROGRAM Exiting...")
                break
            else:
                myPlot.setPlotAngle(plotAngle)

            # For Testing
            if testing > 1:
                print(
                    '%s (%s) became %d.' % (flag, item, myPlot.getPlotAngle()))
        # ---------------------------------------------------------------------
        elif flag == 'scaleX'.lower():

            # Scale X (energy or angle) data
            myPlot.setScaleX(float(item))

            # For Testing
            if testing > 1:
                print('%s (%s) became %d.' % (flag, item, myPlot.getScaleX()))

        # ---------------------------------------------------------------------
        elif flag == 'scaleY'.lower():

            # Scale X (energy or angle) data
            myPlot.setScaleY(float(item))

            # For Testing
            if testing > 1:
                print('%s (%s) became %d.' % (flag, item, myPlot.getScaleY()))

        # ---------------------------------------------------------------------
        elif flag == 'annotate'.lower():

            # Annotation exists as the item passed in
            myPlot.setAnnotateIt(item)

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item,
                                              myPlot.getAnnotateIt()))

        # ---------------------------------------------------------------------
        elif flag == 'annotateX'.lower():

            # Assigning to plot
            myPlot.setAnnotateX(float(item))

            # For Testing
            if testing > 1:
                print(
                    '%s (%s) became %d.' % (flag, item, myPlot.getAnnotateX()))

        # ---------------------------------------------------------------------
        elif flag == 'annotateY'.lower():

            # Assigning to plot
            myPlot.setAnnotateY(float(item))

            # For Testing
            if testing > 1:
                print(
                    '%s (%s) became %d.' % (flag, item, myPlot.getAnnotateY()))

        # ---------------------------------------------------------------------
        elif flag == 'save'.lower():

            # Assigning answer to plot
            if not item == '':
                myPlot.setSave(item)
            else:
                # Error - invalid item
                print("ERROR: Input error for plot %d, line %d: Invalid "
                      "Save Item - \"%s\"." % (plotNum, linesUsed, item))
                print("Warning Detected, PROGRAM Exiting...")
                break

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getSave()))

        # ---------------------------------------------------------------------
        elif flag == 'comment'.lower():

            # General Comment for Input file
            print('Comment Detected, ignoring input line.')
            if verbose > 0:
                print('\t Detected Comment: %s' % item)

        # ---------------------------------------------------------------------
        elif flag == 'plotName'.lower():

            # Assigning plot name to plot
            myPlot.setName(item)

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getName()))

        # ---------------------------------------------------------------------
        elif flag == 'dpi'.lower():

            # Changing resolution
            dpiVal = int(item.lower())  # dpiVal used later when saving figure

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, dpiVal))

        # ---------------------------------------------------------------------
        elif flag == 'display'.lower():

            # Changing resolution
            if item.lower() == 'none'.lower():
                myPlot.setShowPlot(False)
            else:
                myPlot.setShowPlot(True)

            # For Testing
            if testing > 1:
                print('Plot will not be shown.')

        # ---------------------------------------------------------------------
        elif flag == 'xMax'.lower():

            # Changing to integer
            xMax = int(item.strip().lower())
            myPlot.setXMax(xMax)

            # For Testing
            if testing > 1:
                print('x Maximum changed to %d.' %(myPlot.getXMax()))

        # ---------------------------------------------------------------------
        elif flag == 'done\n'.lower() or flag == 'stop\n'.lower() or flag == \
                'quit\n'.lower() or flag == 'exit\n'.lower():

            # Stopping program

            # For Testing
            if testing > 1:
                print('Stopping program, signal received...' %(
                    myPlot.getXMax()))

            break

        # ---------------------------------------------------------------------
        elif flag == 'legendX'.lower():

            # Obtaining number from item
            lgdLocX = float(item.lower().strip())

            # Must be from 0 to 1 to have lower left corner in figure
            myPlot.setLocX(lgdLocX)

            # For Testing
            if testing > 1:
                print('Legend X location changed to %.1f.' %(myPlot.getLocX()))

        # ---------------------------------------------------------------------
        elif flag == 'legendY'.lower():

            # Obtaining number from item
            lgdLocY = float(item.lower().strip())

            # Must be from 0 to 1 to have lower left corner in figure
            myPlot.setLocY(lgdLocY)

            # For Testing
            if testing > 1:
                print('Legend Y location changed to %.1f.' %(myPlot.getLocY()))

        # ---------------------------------------------------------------------
        elif flag == 'reformat'.lower():

            # Assigning plot name to plot
            if item == 'yes'.lower():
                myPlot.setReformat(1)
            elif item == 'no'.lower():
                myPlot.setReformat(0)
            else:
                # Error - invalid item
                print("ERROR: Input error for plot %d, line %d: Invalid "
                      "Reformat Item - \"%s\"." % (plotNum, linesUsed, item))
                print("Warning Detected, PROGRAM Exiting...")
                break

            # For Testing
            if testing > 1:
                print('%s (%s) became %s.' % (flag, item, myPlot.getReformat()))

        # ---------------------------------------------------------------------
        else:
            # Error - invalid flag
            print("ERROR: Input error for plot %d, line %d: Invalid Invalid "
                  "Flag - \"%s\"." % (plotNum, linesUsed, flag))
            print("Warning Detected, PROGRAM Exiting...")
            break

###############################################################################
# Closing input file; program end
inputFile.close()
printlines(1)
print('-------------------------------------------------------------')
print('%d plots created from file \"%s\".' %(plotNum, inputFile.name))
print('Hopefully your data was good!  If not, better luck next time!')
print('Input file \"%s\" closed, program complete.  Exiting...' %
      inputFile.name)
print('-------------------------------------------------------------')
printlines(1)

###############################################################################
