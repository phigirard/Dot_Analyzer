#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	Dot_Analyzer.py
#	Release v12.0
#
#	Copyright 2022 - BSD-3-Clause license
#                                                                             
#******************************************************************************/

#@ File impFile (label="Select the  image to analyse ", style="file")
#@ Boolean imageScale (label="Measure scale bar on the image ", value = False, persist=true)
#@ String msg1 (visibility=MESSAGE, value="----------- If you know the distance in pixels, use these 2 numeric fields: ------------", required=False) 
#@ Double measured (label="Distance in pixels ", description="Indicate a distance in pixels",value=85.0, stepSize=0.1, persist=False)
#@ Double known (label="Known Distance in nm ", description="Indicate a known distance in nm",value=200.0,stepSize=0.1, persist=False)
#@ String msg2 (visibility=MESSAGE, value="-------------------------------------------------------------------------", required=False)
#@ Integer minSize (label="Minimal size of particles in pixels ", value=10, persist=True) 
#@ Boolean thresholding(label="Automatic Threshold ", value = False, persist=true) 
#@ String vorodiagram (label="Choose the diagram to display ", choices={"Voronoi Diagram", "Voronoi/Delaunay Diagram"}, style="radioButtonHorizontal", value = choices[1],persist=False)
#@ Boolean ripleygraph (label="Besag's L Function ", value = True, persist=true) 
#@ Boolean pcfgraph (label="Pair Correlation Function ", value = True, persist=true) 
#@ Boolean ocfgraph (label="Bond-Orientational Correlation Function ", value = True, persist=true) 
#@ Boolean SaveSpacing (label="Save Spacing and Order in a table ", value = True, persist=true) 

#@ DatasetIOService io
#@ UIService uiService
#@ LogService log
#@ CommandService command
#@ ConvertService convertService
#---------------------------------------------------------------------------------------------------------
#import class ImageJ1
from ij import IJ, ImagePlus, WindowManager, Prefs
from ij.io import Opener, OpenDialog
from ij.gui import  GUI, GenericDialog, Roi, PointRoi,PolygonRoi, YesNoCancelDialog, WaitForUserDialog, Plot, PlotWindow, ShapeRoi
from ij.plugin import ContrastEnhancer, RoiEnlarger
from ij.plugin.frame import RoiManager, ThresholdAdjuster
from ij.measure import ResultsTable , Measurements, Calibration, CurveFitter
from ij.plugin.filter import Analyzer,  BackgroundSubtracter
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.process import ImageProcessor, ImageConverter, ColorProcessor, ByteProcessor


import os
from os import path

# Java class ---------------------------------------------------------------------------------------
from java.io import File
from java.lang import Double, Integer, Short, Thread, String, InterruptedException
from java.awt import Color, Font, BasicStroke, Frame, BorderLayout, FlowLayout
from java.text import NumberFormat, DecimalFormat, SimpleDateFormat, DecimalFormatSymbols
from java.util import Locale, Date, Calendar, TimeZone, Iterator, Vector
from java.awt.image import BufferedImage, IndexColorModel
from java.awt.geom import Rectangle2D, Ellipse2D

from org.jfree.chart import ChartPanel, JFreeChart
from org.jfree.chart.axis import NumberAxis
from org.jfree.chart.plot import XYPlot, ValueMarker
from org.jfree.chart.renderer.xy import XYLineAndShapeRenderer
from org.jfree.data.xy import XYDataset, XYSeries, XYSeriesCollection

from javax.swing import JFrame, JDialog, JOptionPane,JPanel, JLabel, JComboBox,JCheckBox, JFormattedTextField, JButton, SwingConstants, GroupLayout
from javax.swing.border import EmptyBorder
from javax.swing.GroupLayout import Alignment
from javax.swing.LayoutStyle import ComponentPlacement


from math import sqrt, atan2, cos, sin, pi, acos, log, exp, floor, isnan
#---------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------
#---------------          CONSTANTS            -----------------
#---------------------------------------------------------------


title = "Dot Analysis"

Prefs.blackBackground = True
MAXSIZE = Double.POSITIVE_INFINITY # for ParticleAnalyzer

#measured = 85.0; // distance in pixels
#known = 200.0; // known distance in nm
#conversion= 200/85 # 1 (for a magnitude 50 00kx, 200 nm = 85 pix)


saturated = 0.35
radius = 6 # sphere radius to show the dot spacing

# headings of the CSV file
headings = ["Filename", "Polymer","Loading", "concentration (mg/ml)", "Speed (V)", "Date (yy/mm/dd)","Number of dot", "Spacing (nm)", "Stdev (nm)","Sterror (nm)","Order"]

# maximum number of neighbors (for calculation)
maxNeighbors = 12

# color value of the "glasbey inverted" LUT for the Voronoi/Delaunay image
dotColor = 255
voronoiColor = 44
delaunayColor = 118

voronoi = (vorodiagram == "Voronoi Diagram")

#---------------------------------------------------------------
#----------------- All Functions for analysis  -----------------
#---------------------------------------------------------------

#
def scaleDialog(defaultValue):
	gd = GenericDialog("Scale Bar")
	gd.addNumericField("Known distance", defaultValue, 2, 8,  "nm")
	gd.showDialog()
	if gd.wasCanceled() :
		return defaultValue
	v = gd.getNextNumber()
	if (gd.invalidNumber()) :
		return defaultValue
	else :
		return v


# calculate the Euclidean distance between two vectors
def euclidean_distance(row1, row2):
	distance = 0.0
	for i in range(len(row1)):
		distance += (row1[i] - row2[i])**2
	return sqrt(distance)
	
def dot(row1, row2):
	dotproduct = 0
	for i in range(len(row1)):
		dotproduct += row1[i]*row2[i]
	return dotproduct

# Locate the most similar neighbors (based on k-Nearest Neighbors algorithm)
def get_neighbors(centroids_, idx_, num_neighbors):
	distances = list()
	for row in range(len(centroids_)):
		if row != idx_ :
			dist = euclidean_distance(centroids_[idx_], centroids_[row])
			angl = atan2(centroids_[idx_][1]-centroids_[row][1],centroids_[idx_][0]-centroids_[row][0])
			distances.append((row, dist, angl))
	distances.sort(key=lambda tup: tup[1])
	neighbors = list()
	neighborsdist = list()
	neighborsangl = list()
	for i in range(num_neighbors):
		neighbors.append(distances[i][0])
		neighborsdist.append(distances[i][1])
		neighborsangl.append(distances[i][2])
	return [neighbors,neighborsdist, neighborsangl]


#check if the ROI is at the edge of the image. This ROI are not considered in the calculation of the spacing and order parameter				
def isRoiAtEdge(polygon, imsize) : 
	edgeNB = False
	polypts = map(list, zip(polygon.xpoints,polygon.ypoints))
	for j in range(polygon.npoints):
		point = polypts[j]
		for i in  range(2) :
			edgeNB  |= (point[i] <=1 or point[i] >= (imsize[i]-1))
	return edgeNB		
		
# find the array index of a value
def findIdx(array_, value_) :
	for i in range(len(array_)) :
		if array_[i] == value_ :
			return i

def isNeighbors(roi1, roi2):
	polyRoi1 = ShapeRoi(roi1)
	polyRoi2 = ShapeRoi(RoiEnlarger.enlarge(roi2, 2))
	roiAND =polyRoi2.and(polyRoi1)
	return (roiAND.getLength()  != 0 )
	

#------------- Ripley's K-Function or reduced second-moment function (Ripley, 1981)  --------------------*/	
#	Ripley, B. 1981. Spatial Statistics.  John Wiley, Chichester.
def RipleyKFunction(w_,h_, centroids_ ,besagFunction, resolution, conversion):
	print "Plot the Besag's L Function"
	maxd= int(min(w_,h_))
	maxres = maxd*resolution
	plotTitle = "Ripley's K Function"
	plotTitleY = "K(r)"
	nrow = len(centroids_)
	series = XYSeries(plotTitle)
	for t in range(maxres):
		kfunc = 0
		kfuncX = (t+1)/resolution
		for i in range(nrow) :
			for j in range(nrow) :
				if (j!=i) :
					kfunc+=weightFunction(centroids_[i], centroids_[j],w_,h_, kfuncX)
		kfunc *= w_*h_/(nrow*(nrow-1))*conversion
		kfuncX *=conversion
		#return the Besag's L function (1977)
		if (besagFunction) :
			kfunc = sqrt(kfunc/pi)- kfuncX
			plotTitle = "Besag's L Function"
			plotTitleY = "L(r)"
		series.add(kfuncX, kfunc)
	
	dataset = XYSeriesCollection(series) 
	yaxis = NumberAxis(plotTitleY)
	xaxis = NumberAxis("Distance r (nm)")
	r = XYLineAndShapeRenderer()
	r.setSeriesPaint(0, Color.BLUE)
	#r.setSeriesShape(0, Ellipse2D.Double(-3.0,-3.0,6.0,6.0))
	xyplot = XYPlot(dataset, xaxis, yaxis, r)
	xyplot.setBackgroundPaint(Color.white)
	chart = JFreeChart(xyplot)
	chart.removeLegend()
	impPlot = IJ.createImage(plotTitle, "RGB", 512, 512, 1);
	imagePlot = impPlot.getBufferedImage()
	chart.draw(imagePlot.createGraphics(), Rectangle2D.Float(0, 0, impPlot.width, impPlot.height))
	impPlot.setImage(imagePlot)
	return impPlot
		

def weightFunction(row1, row2, w_, h_, dvar):
	dij = euclidean_distance(row1, row2)
	minx = min(row1[0], w_-row1[0])
	miny = min(row1[1], h_-row1[1])
	dmin = min(minx, miny)
	wF=0
	if (dij<=dvar) :
		if (dij<=dmin):
			wF = 1
		elif (dvar*dvar<=minx*minx+miny*miny):
			wF = 1/(1-acos(dmin/dvar)/pi)
		elif (dvar*dvar>minx*minx+miny*miny):
			wF = 1/(1-(acos(minx/dvar)+acos(miny/dvar)+pi/2)/(2*pi))
	return wF

#Penttinent et al. (1992) Marked point processes in forest statistics. For. Sci. 38, 806-824.
def PairCorrelation(w_,h_,centroids_, resolution, conversion):
	print "Plot the pair correlation function"
	maxd= int(min(w_,h_))
	pcf = []
	pcfX= []
	pcf.append(0)
	pcfX.append(0)
	sd=0
	nrow = len(centroids_)
	invlam=w_*h_/nrow
	plotTitle = "Pair correlation Function"
	series = XYSeries(plotTitle)
	for t in range(1,maxd*resolution):
		pcf = 0
		pcfX = t/resolution
		for i in range(nrow) :
			for j in range(nrow) :
				if (j!=i) :
					pcf +=Epanechnikov(centroids_[i], centroids_[j], pcfX, invlam);


		#sd= edge correction factor (Stoyan et al. 1987 Stochastic Geometry and its application. Wiley, New York)
		sd = w_*h_ - pcfX*(2*(w_+h_)-pcfX)/pi
		pcf *= invlam*invlam/(2*pi*pcfX*sd)
		pcfX *=conversion
		series.add(pcfX, pcf)
	dataset = XYSeriesCollection(series) 
	yaxis = NumberAxis("g(r)")
	xaxis = NumberAxis("Distance r (nm)")
	r = XYLineAndShapeRenderer()
	r.setSeriesPaint(0, Color.BLUE)
	r.setAutoPopulateSeriesStroke(False)
	r.setDefaultStroke(BasicStroke(float(2.0)))
	xyplot = XYPlot(dataset, xaxis, yaxis, r)
	xyplot.setBackgroundPaint(Color.white)
	marker = ValueMarker(1.0)
	marker.setPaint(Color.RED)
	marker.setStroke(BasicStroke(float(2.0)))
	xyplot.addRangeMarker(marker)
	chart = JFreeChart(xyplot)
	chart.removeLegend()
	impPlot = IJ.createImage(plotTitle, "RGB", 512, 512, 1);
	imagePlot = impPlot.getBufferedImage()
	chart.draw(imagePlot.createGraphics(), Rectangle2D.Float(0, 0, impPlot.width, impPlot.height))
	impPlot.setImage(imagePlot)
	return impPlot


def Epanechnikov(row1, row2, dvar, invlam):
	dij = euclidean_distance(row1, row2)
	diff = dij-dvar
	delta = 0.15*sqrt(invlam)
	Epa = 0
	if (abs(diff)<delta) :
		Epa =  3*(1-diff*diff/(delta*delta))/(4*delta)
	return Epa

def OrderCorrelation(w_,h_, centroids_,  neighbors_):
	print "Plot the Bond-Orientational Correlation Function"
	maxd = int(min(w_,h_))
	nrow = len(neighbors_)
	invlam = w_*h_/nrow
	delta = 0.15*sqrt(invlam)
	plotTitle = "Bond-Orientational Correlation Function"
	series = XYSeries("RawData")
	ocf = []
	ocfX = []
	for t in range(maxd-1):
		bocf = 0
		bocfX = t+1
		psi_real = 0
		nbpts = 0
		for i in range(nrow-1):
			for j in range(i+1,nrow):
				dij = euclidean_distance(centroids_[i], centroids_[j])
				diff = dij-bocfX
				if (abs(diff)<delta):
					nbpts+=1
					for k in range(6):
						angli = neighbors_[i][1][k]
						angltemp = 0
						psi_realtemp = 0
						add = True
						for l in range(6):
							anglj = neighbors_[j][1][l]
							if (l == 0): 
								angltemp = anglj
							elif (anglj == angltemp):
								add = False		
							psi_realtemp += abs(cos(6*(anglj-angli)))
						if (add) :
							psi_real+=psi_realtemp

		bocfX *= conversion
		if (nbpts == 0):
			bocf = 0
		else :
			bocf= psi_real/nbpts/36
			if (bocf != 0):
				ocf.append(bocf)
				ocfX.append(bocfX)
				series.add(bocfX, bocf)
	newsize = len(ocf)
	
	# Fitter
	fitter = CurveFitter(ocfX, ocf)
	fitter.doFit(CurveFitter.EXPONENTIAL, False) #a*exp(b*x)
	param_values = fitter.getParams()
	chi2 = fitter.getRSquared
	fitseries = XYSeries("Fit")
	for xt in ocfX :
		fitseries.add(xt, fitter.f( param_values, xt))
	dataset = XYSeriesCollection() 
	dataset.addSeries( series )
	dataset.addSeries( fitseries )
	yaxis = NumberAxis("g6(r)")
	xaxis = NumberAxis("Distance r (nm)")
	yaxis.setRange(0, 1)
	if newsize > 0 :
		xaxis.setRange(0, ocfX[newsize-1])
	renderer= XYLineAndShapeRenderer()
	renderer.setSeriesPaint( 0 , Color.BLUE )
	renderer.setSeriesPaint( 1 , Color.RED )
	renderer.setAutoPopulateSeriesStroke(False)
	renderer.setSeriesStroke( 0 , BasicStroke( float(1.0)) )
	renderer.setSeriesStroke( 1 , BasicStroke( float(1.0)) )
	xyplot = XYPlot(dataset, xaxis, yaxis, renderer)
	xyplot.setBackgroundPaint(Color.white)
    
	chart = JFreeChart(xyplot)
	chart.removeLegend()
	impPlot = IJ.createImage(plotTitle, "RGB", 512, 512, 1);
	imagePlot = impPlot.getBufferedImage()
	chart.draw(imagePlot.createGraphics(), Rectangle2D.Float(0, 0, impPlot.width, impPlot.height))
	impPlot.setImage(imagePlot)
	return impPlot
	
	
def sVal(dd):
	rounder = pow(10.0, 3)
	return str(round(rounder*dd)/rounder)	
	

def init(date) :
	polymer = Vector([ "PS","P2VP","PDMS", "PMMA" ])
	rate = Vector(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
	strDate = SimpleDateFormat("yyyy/MM/dd")
	speed = 6
	spinCoating = False
	oldFile = False
	
	copolyFormat = NumberFormat.getInstance(Locale.US)
	copolySymbols = copolyFormat.getDecimalFormatSymbols()
	concFormat = NumberFormat.getNumberInstance(Locale.UK)
	copolySymbols.setGroupingSeparator(' ')
	copolyFormat.setDecimalFormatSymbols(copolySymbols)
	
	
	parent = IJ.getInstance()
	instance = JDialog(parent,"Save Spacing & Order")
	instance.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
	instance.setBounds(100, 100, 500, 210)
	instance.setLayout(BorderLayout())
	contentPanel = JPanel()
	contentPanel.setBorder(EmptyBorder(5, 5, 5, 5))
	#Line 1
	labelpoly = JLabel("Polymer:")
	comboBoxPolymer1 = JComboBox(polymer)
	comboBoxPolymer2 = JComboBox(polymer)
	comboBoxPolymer1.setSelectedIndex(0)
	textCopolymer1 = JFormattedTextField(copolyFormat)
	textCopolymer1.setColumns(8)
	textCopolymer1.setValue(52400)
	textCopolymer1.setHorizontalAlignment(SwingConstants.CENTER)
	comboBoxPolymer2.setSelectedIndex(1)
	textCopolymer2 = JFormattedTextField(copolyFormat)
	textCopolymer2.setColumns(8)
	textCopolymer2.setValue(28100)
	textCopolymer2.setHorizontalAlignment(SwingConstants.CENTER)
	
	#Line 2
	lblLoadingRate = JLabel("Loading rate:")
	comboBoxLoading = JComboBox(rate)
	comboBoxLoading.setSelectedIndex(5)
	lblPolymerConcentration = JLabel("Polymer concentration:")
	textPolyConc = JFormattedTextField(concFormat)
	textPolyConc.setColumns(7)
	textPolyConc.setValue(Integer(5))
	textPolyConc.setHorizontalAlignment(SwingConstants.RIGHT)
	lblMgml = JLabel("mg/ml")
	loading = str(comboBoxLoading.getSelectedItem())
	
	def enable(event):
		if rdbtnSpinCoating.isSelected() :
			lblDippingSpeed.setEnabled(False)
			textDipSpeed.setEnabled(False)
			textDipSpeed.setValue(0)
			lblV.setEnabled(False)
		else :
			lblDippingSpeed.setEnabled(True)
			textDipSpeed.setEnabled(True)
			textDipSpeed.setValue(speed)
			lblV.setEnabled(True)
			
	#Line 3
	rdbtnSpinCoating = JCheckBox("True = Spin, False = Dip coating", actionPerformed = enable)
	rdbtnSpinCoating.setSelected(spinCoating)
	lblDippingSpeed = JLabel("Dipping speed:")
	lblDippingSpeed.setEnabled(True)
	textDipSpeed = JFormattedTextField(concFormat)
	textDipSpeed.setEnabled(True)
	textDipSpeed.setColumns(5)
	textDipSpeed.setValue(speed)
	textDipSpeed.setHorizontalAlignment(SwingConstants.RIGHT)
	lblV = JLabel("V ")
	lblV.setEnabled(True)
	
	#Line 4
	lblDate = JLabel("Date (yyyy/mm/dd):")
	textDate = JFormattedTextField(strDate)
	textDate.setColumns(10)
	textDate.setValue(date)
	chckbxAddInFile = JCheckBox("Add in an existed file")
	chckbxAddInFile.setSelected(oldFile)
	
	groupLayout = GroupLayout(contentPanel)
	groupLayout.setHorizontalGroup(groupLayout.createParallelGroup(Alignment.LEADING).addGroup(groupLayout.createSequentialGroup().addContainerGap()
															.addGroup(groupLayout.createParallelGroup(Alignment.LEADING).addGroup(groupLayout.createSequentialGroup().addComponent(labelpoly)
                                                                     .addPreferredGap(ComponentPlacement.RELATED)
                                                                     .addComponent(comboBoxPolymer1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(textCopolymer1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(comboBoxPolymer2, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(textCopolymer2, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                           .addGroup(groupLayout.createSequentialGroup()
                                                                     .addComponent(lblDate)
                                                                     .addPreferredGap(ComponentPlacement.UNRELATED)
                                                                     .addComponent(textDate, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(chckbxAddInFile))
                                                           .addGroup(groupLayout.createSequentialGroup()
                                                                     .addComponent(rdbtnSpinCoating)
                                                                     .addGap(38)
                                                                     .addComponent(lblDippingSpeed)
                                                                     .addComponent(textDipSpeed, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(lblV))
                                                           .addGroup(groupLayout.createSequentialGroup()
                                                                     .addComponent(lblLoadingRate)
                                                                     .addPreferredGap(ComponentPlacement.RELATED)
                                                                     .addComponent(comboBoxLoading, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addPreferredGap(ComponentPlacement.UNRELATED)
                                                                     .addComponent(lblPolymerConcentration)
                                                                     .addComponent(textPolyConc, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                     .addComponent(lblMgml)))
                                                 .addContainerGap(20, Short.MAX_VALUE)))
	groupLayout.setVerticalGroup(groupLayout.createParallelGroup(Alignment.LEADING)
                                     		.addGroup(groupLayout.createSequentialGroup()
                                               .addContainerGap()
                                               .addGroup(groupLayout.createParallelGroup(Alignment.LEADING)
                                                         .addGroup(groupLayout.createParallelGroup(Alignment.LEADING)
                                                                   .addComponent(comboBoxPolymer1, Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                   .addComponent(textCopolymer1, Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                   .addComponent(comboBoxPolymer2, Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                                   .addComponent(textCopolymer2, Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                         .addGroup(groupLayout.createSequentialGroup()
                                                                   .addGap(6)
                                                                   .addComponent(labelpoly)))
                                               .addPreferredGap(ComponentPlacement.RELATED)
                                               .addGroup(groupLayout.createParallelGroup(Alignment.BASELINE)
                                                         .addComponent(comboBoxLoading, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                         .addComponent(lblLoadingRate)
                                                         .addComponent(lblPolymerConcentration)
                                                         .addComponent(textPolyConc, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                         .addComponent(lblMgml))
                                               .addPreferredGap(ComponentPlacement.RELATED)
                                               .addGroup(groupLayout.createParallelGroup(Alignment.BASELINE)
                                                         .addComponent(rdbtnSpinCoating)
                                                         .addComponent(lblDippingSpeed)
                                                         .addComponent(textDipSpeed, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                         .addComponent(lblV))
                                               .addPreferredGap(ComponentPlacement.RELATED)
                                               .addGroup(groupLayout.createParallelGroup(Alignment.BASELINE)
                                                         .addComponent(chckbxAddInFile)
                                                         .addComponent(lblDate)
                                                         .addComponent(textDate, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))))
	contentPanel.setLayout(groupLayout)
	instance.add(contentPanel, BorderLayout.CENTER) 
	def dispose(event):
		source = event.getSource()
		if source == okButton :
			source.setActionCommand("pressed")
		instance.visible = False
		instance.dispose()
	
	buttonPane = JPanel()
	buttonPane.setLayout(FlowLayout(FlowLayout.RIGHT))
	okButton = JButton("OK", actionPerformed=dispose)
	buttonPane.add(okButton)
	instance.getRootPane().setDefaultButton(okButton)
	cancelButton = JButton("Cancel", actionPerformed=dispose)
	buttonPane.add(cancelButton)
	instance.add(buttonPane, BorderLayout.SOUTH)
	instance.pack()
	GUI.centerOnImageJScreen(instance)
	instance.setResizable(False)
	instance.setAlwaysOnTop(True)
	instance.visible = True
	

	while instance.isVisible():
		try:
			Thread.sleep(1000)
		except: 
			Thread.currentThread().interrupt()
	rowResult = list()
	if len(okButton.getActionCommand()) > 2 :
		copolymer = String.format("%1$s(%2$s)-%3$s(%4$s)",str(comboBoxPolymer1.getSelectedItem()),textCopolymer1.getText(),str(comboBoxPolymer2.getSelectedItem()),textCopolymer2.getText())
		loading = str(comboBoxLoading.getSelectedItem()) 
		concentration= textPolyConc.getText()
		spinCoating = rdbtnSpinCoating.isSelected()
		speed = Double(textDipSpeed.getText())
		date = textDate.getText()
		oldFile = chckbxAddInFile.isSelected()
		rowResult= [oldFile, copolymer, loading, concentration, str(speed), date]
	
	return rowResult
#---------------------------------------------------------------

# clear the console automatically when not in headless mode
uiService.getDefaultUI().getConsolePane().clear()


#close Result Table if opened
if IJ.isResultsWindow() :
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	tw.close()
rt= ResultsTable()

#reset RoiManager or open one
rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()
rm.reset()

			
			
#convert Files from #@ parameters to String and extract the main directory of the data
impPath = impFile.getCanonicalPath()
filename = path.splitext(path.basename(impPath))[0]


#create folder for analysis 
srcDir = path.dirname(impPath)
imageDir = path.join(srcDir, "Analyzed_"+filename) 
if not path.exists(imageDir):
	os.makedirs(imageDir)

#open the file
imp = Opener().openImage(impPath)
width = imp.width
height = imp.height

# creation date of the tiff file in the format yyyy/MM/dd
timestamp = File(impPath).lastModified()
when = Date(timestamp)

#Measure scale bar
if imageScale : 
	imp.show()
	IJ.setTool("rectangle")
	waitDialog = WaitForUserDialog("Scale Bar","Draw a rectangle to fit with the scale bar")
	waitDialog.show()
	measured = imp.getRoi().getBounds().width
	known = scaleDialog(200)
	imp.hide()
conversion = known/measured

restart = True
while (restart) :
#---------- Prepare image for Analyze Particles:
	imp.setRoi(0,0,width,int(height*0.89)) #remove the information area at the bottom of the image
	imptp = imp.crop()        
	ip = imptp.getProcessor()
	ip_src = ip.duplicate().convertToByte(True)
	BackgroundSubtracter().rollingBallBackground(ip,10,False,False,False,False,True)
	ip.smooth()
	ip.smooth()
	
	
	if not thresholding :
		imptp.show()
		ta = ThresholdAdjuster()
		ta.setMethod("Triangle")
		ta.show()
		ta.update()
		waitDialog = WaitForUserDialog("Manual threshold", "Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')") # human thresholding
		waitDialog.show()
		thres_min = ip.getMinThreshold()
		thres_max = ip.getMaxThreshold()
		ta.close()
		imptp.hide()
		IJ.setThreshold(imptp, thres_min, thres_max)
	else :
		IJ.setAutoThreshold(imptp, "Triangle dark")
	IJ.run(imptp, "Convert to Mask", "")
	
	#Create composite
	height = imptp.height
	cp = ColorProcessor(width,height)
	rPels = []
	rS = imptp.getProcessor().getPixels()
	GS = ip_src.getPixels()
	for pp in range(width*height) :
		rPels.append((GS[pp])|(rS[pp]))
	cp.setRGB(rPels, GS, GS)
	cimp = ImagePlus("RED(Binary)/GRAY(Source)", cp)
	cimp.show()
	question = JOptionPane.showConfirmDialog(None,"Are you ok with the segmentation?")
	cimp.hide()
	if (question == JOptionPane.NO_OPTION) :
		thresholding = False
		restart = True
	elif (question == JOptionPane.CANCEL_OPTION):
		voronoi = False
		ripleygraph = False 
		pcfgraph = False
		ocfgraph = False
		SaveSpacing = False
		break
	else :
		restart = False
	
	


# Detect signal ROI from the previous thresholded image (background ROI = inverse of signal ROI)
p = PA(PA.ADD_TO_MANAGER, Measurements.CENTROID, rt ,50, MAXSIZE)
p.setRoiManager(rm)
p.analyze(imptp)


xcentroid = rt.getColumn(rt.getColumnIndex("X"))
ycentroid = rt.getColumn(rt.getColumnIndex("Y"))
xInt = []
yInt = []
for i in range(len(xcentroid)):
	xInt.append(int(xcentroid[i])) #integer version of xcentroid
	yInt.append(int(ycentroid[i])) #integer version of ycentroid
datadots = map(list, zip(xcentroid, ycentroid))

	
								
# create voronoi image with color code neighbor number										
impSpacing = IJ.createImage("Spacing", "8-bit black", width, height, 1)
ip = impSpacing.getProcessor()
ip.setLineWidth(3)
for i in range(len(xcentroid)):
	ip.setColor(Color.WHITE)
	ip.fillOval(xInt[i]-radius,yInt[i]-radius,2*radius,2*radius)
	
IJ.run(impSpacing, "Voronoi", "")
IJ.setThreshold(impSpacing,1, 255)
IJ.run(impSpacing, "Convert to Mask", "")
IJ.run(impSpacing, "Invert", "")


neighborArray=[]
neighbored=0
mostNeighbors=0
for i in range(len(xcentroid)): 
	rt.reset()
	IJ.doWand(impSpacing, xInt[i], yInt[i], 0.0, "8-connected")
	IJ.run(impSpacing, "Enlarge...", "enlarge=2")
	p = PA(PA.SHOW_NONE,Measurements.CENTROID, rt, 0 , Double.POSITIVE_INFINITY)
	p.analyze(impSpacing)
	neighbored = rt.size() -1
	neighborArray.append(neighbored)
	if neighbored>mostNeighbors :
		mostNeighbors=neighbored


rm.reset()
rt.reset()
p = PA(PA.ADD_TO_MANAGER, Measurements.CENTROID, None ,0, MAXSIZE)
p.setRoiManager(rm)
p.analyze(impSpacing)
IJ.run(impSpacing, "Invert", "")
#IJ.saveAs(impSpacing, "TIFF", path.join(imageDir,filename+"_impSpacing.tif"))

#attribute voronoi roi to dot position (xcentroid/ycentroid): voronoi i -> dot position roiCom[i]
roiComp=[]
for roi in rm.getRoisAsArray() :
	for i in range(len(xcentroid)) :
		if roi.containsPoint(xcentroid[i],ycentroid[i]):
			roiComp.append(i)	


print "Calculation of spacing and order parameter"
neighbors = []
for i in range(len(xcentroid)):
	neighbors.append(get_neighbors(datadots, i, maxNeighbors))

rmSize =rm.getCount()
polygons = []
for i in range(rmSize):
	idx = roiComp[i]
	roi = rm.getRoi(i)
	mark = neighborArray[idx]
	poly = roi.getFloatPolygon().getConvexHull()
	polygons.append(poly)
	ip.setColor(voronoiColor)
	ip.setLineWidth(1)
	ip.draw(roi)
	if not isRoiAtEdge(poly, [width, height]) :
		ip.setColor(mark)
		ip.fill(roi)
		
meandist = 0
squaredist = 0
nbdist = 0
phi = 0	
for i in range(rmSize):
	psi_real = 0
	psi_img = 0
	idx = roiComp[i]
	mark = neighborArray[idx]
	if not isRoiAtEdge(polygons[i], [width, height]) :
		ip.setColor(delaunayColor)
		ip.setLineWidth(2)
		for j in range(maxNeighbors) :
			idx2 = int(neighbors[idx][0][j])
			i2 = findIdx(roiComp,idx2)
			if isNeighbors(rm.getRoi(i), rm.getRoi(i2)):
				if i2 > i or (i2<i and isRoiAtEdge(polygons[i2], [width, height])) :
					if not voronoi :
						ip.drawLine(int(datadots[idx][0]),int(datadots[idx][1]),int(datadots[idx2][0]),int(datadots[idx2][1]))
					meandist += neighbors[idx][1][j]
					squaredist += neighbors[idx][1][j]*neighbors[idx][1][j]
					nbdist+=1
				angl = neighbors[idx][2][j]
				psi_real += cos(6 * angl)
				psi_img += sin(6 * angl)
			
		phi += sqrt((psi_real * psi_real + psi_img * psi_img))/mark
suffix = "_Voronoi"
if not voronoi :
	suffix = "_Voronoi-Delaunay"
	ip.setColor(dotColor)
	for i in range(rmSize):	
		ip.fillOval(xInt[i]-radius,yInt[i]-radius,2*radius,2*radius)

meandist /= nbdist #  measurement in pixels
stdev = sqrt( (squaredist - nbdist * meandist * meandist) / nbdist) #  measurement in pixels
meandist *= conversion # measurement in nm
stdev *= conversion # measurement in nm
stderror = stdev / sqrt(nbdist)
phi = phi/len(neighbors)
#Save data in array	
dotResult=[nbdist,meandist, stdev, stderror, phi]	

impSpacing.updateAndDraw()
IJ.run(impSpacing,"Select None", "")	
IJ.run(impSpacing, "glasbey inverted", "")
IJ.resetMinAndMax(impSpacing)
impSpacing.show()	
IJ.saveAs(impSpacing, "TIFF", path.join(imageDir,filename+suffix+".tif"))
rm.runCommand(impSpacing,"Show None")

#draw calibration bar
mostNeighbors = max(neighborArray)
stepsize=int(floor(256/mostNeighbors))
w = stepsize*mostNeighbors
ipBar = ip.createProcessor(w,50)
ipBar.setColor(0)
ipBar.fill()
step=0
for c in range(mostNeighbors):
	ipBar.setColor(c+1)
	ipBar.fillRect(step, 0, step+stepsize, 30)
	step+=stepsize
ipNew = ColorProcessor(ipBar.createImage())
offset=4
ipNew.setColor(Color.white)
middlestep = int(floor(stepsize/2))
ipNew.setFont(Font("SansSerif", Font.BOLD, 12))
for c in range(mostNeighbors):
	ipNew.drawString(str(c+1),middlestep+ c*stepsize-offset, 48)
impBar = ImagePlus("Calibration Bar", ipNew)
IJ.saveAs(impBar, "TIFF", path.join(imageDir,filename+"_CalibrationBar.tif"))
impBar.show()

		

if ripleygraph :
	ripleyplot = RipleyKFunction(width, height, datadots, True, 1, conversion)
	ripleyplot.show()
	IJ.saveAs(ripleyplot, "TIFF", path.join(imageDir,filename+"_BesagFunction.tif"))
if pcfgraph :
	PCFplot = PairCorrelation(width, height, datadots, 1, conversion)
	PCFplot.show()
	IJ.saveAs(PCFplot, "TIFF", path.join(imageDir,filename+"_PCF.tif"))
if ocfgraph :
	OCFplot = OrderCorrelation(width, height, datadots, neighbors)
	OCFplot.show()
	IJ.saveAs(OCFplot, "TIFF", path.join(imageDir,filename+"_OCF.tif"))

oldfile = False

if SaveSpacing :			
	addRow = init(when)
	
	if len(addRow) >0 :
		oldfile = addRow[0]
		addRow[0]= filename
		addRow.extend(dotResult)
	
	IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column")
	if oldfile :
		op = OpenDialog("Choose CSV file to open", "")
		tablePath = op.getPath()
		Opener().openTable(tablePath)
		tableTitle = path.basename(tablePath)
		dotTable =  WindowManager.getWindow(tableTitle).getTextPanel().getResultsTable()
		dotHeadings = dotTable.getHeadings()
		if not len(headings) == len(dotHeadings) :
			oldfile = False
		else :
			rtsize = dotTable.size()
			for j in range(len(headings)) :
				dotTable.setValue(headings[j],rtsize, addRow[j])
			dotTable.saveAs(tablePath)
	if not oldfile :
		dotNewTable= ResultsTable()
		for j in range(len(headings)) :
			dotNewTable.addValue(headings[j],addRow[j])
		dotNewTable.show("Dot Analysis Results")
		dotNewTable.saveAs(path.join(imageDir,filename+"Results.csv"))

		
print "Results:"
for i in range(len(dotResult)):
	print headings[i+6]+" = "+ str(dotResult[i])

rm.close()
print 'END'
