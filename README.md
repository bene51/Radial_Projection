### Motivation

[We](http://www.mpi-cbg.de/research/research-groups/jan-huisken.html) use a
Selective Plane Illumination Microscope (SPIM) for multi-view acquisition of
zebrafish embryos. Two high-speed cameras capture data with up to 62 fps each,
which leads to enormous data rates. The zebrafish embryo, however, has a
spherical shape. We exploit this fact by computing radial intensity projections
in real-time, which reduces the amount of data to be saved ~250x. The set of
[Fiji](http://www.fiji.sc) plugins in this repository are used to create and
further process the projections.

The entire workflow and the experiments performed so far were published in
[Nature Communications, 2013](http://www.nature.com/ncomms/2013/130725/ncomms3207/full/ncomms3207.html).


### Overview

* Fitting a sphere to the data ([Fit_Sphere](https://github.com/bene51/Radial_Projection/wiki/Fit_Sphere)) 
The first step is usually to determine center and radius of the embryo. This
information is used later on for the projection.

* Radial projection ([File_Max_Projection](https://github.com/bene51/Radial_Projection/wiki/File_Max_Projection))  
Originally, radial maximum projection was done in real-time, i.e. data was
processed as it came from the cameras, and only the projections were saved.
For people who do not have the same setup, especially not the same cameras, we
created a version that reads data from the hard-drive instead.

* Multiview fusion ([Multiview_Fusion]())  
SPIM makes it very easy to rotate the sample to capture another view of the
same data, which can be used to improve the quality of the data. The radial
projections are done separately for each view, the Multview_Fusion plugin is
used to fuse the different views.

* Map Projection ([Map_Projection]())  
Spherical data can be projected onto a planar rectangle. How to do this has
been studied extensively in geographical sciences: How to map the earth onto
a 2D maps. We use software from [jhlabs](http://www.jhlabs.com) to implement
some common map projections for our data.


### More tools

* [3D_View]()  
A plugin for visualizing the projected data directly in 3D using ImageJ/Fiji's
[3D_Viewer](http://3dviewer.neurofly.de).

* [Resample_Projection]()  
Projected data is usually big, too big to view it interactively in 3D. This
plugin can be used to downsample it.

* [Create_Overlay]()  
The [Map_Projection]() plugin is capable of producing maps, files showing
longitude/latitude lines and files indicating the contributions of individual
views. This plugin can be used to create an overlay of all that data.


