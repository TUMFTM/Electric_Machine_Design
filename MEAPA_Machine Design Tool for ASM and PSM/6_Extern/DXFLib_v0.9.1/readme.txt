DXFLib
------------------
Version 0.9.1
Copyright (c) 2009-2011 Grzegorz Kwiatek, see license.txt file for details.

Description
-----------
DXFLib is a simple library written for MATLAB. It allows creating DXF files containing the 
following entities:
* Points
* Lines 
* Primitives
* Polysurfaces

The library handles DXF Layers and ACI Colors (RGB colors are converted to ACI counterparts).

Installation
------------
Just unpack m-files to your working directory or create a separate directory and add its location
to MATLAB path.

Usage
-----
The easiest is to run dxf_test.m script and then open it in MATLAB editor and compare the generated
DXF files with the code. You'll need any DXF viewer that is capable to handle 3D. I recommend
to use free Bentley View software (http://www.bentley.com/en-US/Products/Bentley+View/). The 
advantage of using Bentley View is that this software is capable to PRINT 3-dimensional PDF 
files (the examples of DXF/PDF files are included in the library). 

Contact
-------
Author: Grzegorz Kwiatek, 
Mailto: taquart@gmail.com
DXFLib page: http://www.sejsmologia-gornicza.pl/?page_id=172