# Electric_Machine_Design
This Tool enables a holistic electric machine design and efficiency diagram generation using few input parameters.

## How-To
The tool is made up of 2 main Meta-Models:
  - Electric Machine Design
  - Efficiency Diagram Generation
  
In order to start the overall simulation, you have to run the MATLAB editor in the file Main.m.
This starts the GUI-Interface, where you can insert the main input parameters of the desired electric machine.

With the start-button the first Meta-Model for the electric machine design is inititated. This step takes about 3 seconds and calculates the main dimensions, boundary conditions of the active parts, the rotor and stator design, as well as the electrical, magnetic and mechanical properties of the desired electric machine.

After the electric machine design is completed, the interface for the efficiency diagram generation appears. Here, you can select the different loss types to be considered in the diagram generation, as well as adjust the control and boundary conditions of the efficiency diagram. The start-button initiates the generation of the efficiency diagram, which takes about 7 seconds. When the efficiency diagram is displayed, the different loss type diagrams can be saved to the result-file. This concludes the electric machine design.

## Input Parameters
The primary parameters, that need to be established are:
  - Nominal power
  - Rotational velocity
  - Nominal voltage level
  - Number of pole pairs
  - Power factor
  - Number of phases

Then, specific options need to be inserted, such as:
  - Conductor material
  - Connection of winding
  - Type of cooling
  - Magnet pattern
  
The aporiximated values do not need to be changed, since they represent values from literature and are mostly unknown specifically for a certain application.

In the section for the options, the type of slot calculation is chosen and the slot animation can be selected.

## Output Parameters

The output of the tool is mainly an Excel-file with all calculated machine parameters and the generated efficiency diagrams.

## Necessary Software

MATLAB® 2017 and above.

## Validation
The tool was validated using know machine parameters of the BMW i3, Nissan Leaf and VW e-Golf.


Contact person: [Svenja Kalt](mailto:svenja.kalt@ftm.mw.tum.de)
