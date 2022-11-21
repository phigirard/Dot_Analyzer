# Dot Analyzer for SEM micrographs
an ImageJ/Fiji python script for dot spacing and order analysis of SEM micrographs 

--------------------------------------------------------------
This document describes the workings of the ImageJ/Fiji python script **Dot_Analyzer.py**, written by Philippe Girard ([email](philippe.girard@ijm.fr)).

This Python script was written to analyse SEM micrograph as shown here:
<p align="center">
	<img src="./images/Fig1.png" width="600" 
         alt="SEM micrograph">
<br>
<i>Fig. 1:</i> A typical SEM micrograph.</p>


## 1. Installation instructions
* Download the [python script file](https://github.com/phigirard/Dot_Analyzer/blob/main/Dot_Analyzer.py) into your computer.
* For the script to appear in the menu, you should saved in the `ImageJ2.app/scripts` or the `ImageJ2.app/plugins/Scripts` directory (or a subdirectory thereof).
* Otherwise use File > Open…  to open the script `Dot_Analyzer.py` in the [script Editor](https://imagej.net/scripting/script-editor) of ImageJ/Fiji
and Click Run on the bottom of the script editor window (you can also go to : Run > Run in the Script Fiji menu).
<p align="center">
	<img src="./images/Fig2.png" width="400">
<br>
<i>Fig. 2:</i> Script `Dot_Analyzer.py` in the script Editor of ImageJ/Fiji.</p>



## 2. The “Parameters” main window.
When you start the plugin, you have to define different parameters for the analysis and to select the different diagrams/plots that you want to visualize (Fig. 3).<br>
<p align="center">
	<img src="./images/Fig3.png" width="400">
<br>
</p>

						      
You must indicate:<br>

1. Select the image to analyse: Choose the micrograph image.<br>
	      
2. Measure scale bar on the image: If you do not know the size in nm of a pixel or if your image is not automatically calibrated spatially, you should select this option. In this case before the analysis a dialog box will ask the user to draw a rectangle to fit with the scale bar (it is recommanded to zoom in the region of the scale bar of the image) as shown in Fig 4. <br>

