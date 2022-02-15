# ProjectAion

This project's purpose is to identify area of FERMI source uncertainty areas, or indeed entire sources, which have been missed by the ongoing observation efforts of Schinzel et al over the past several LAT point source catalog releases. We take a subset of unassociated 4FGL sources about -40 DEC, and for each look for pointings in the projects that are within a specific fixed range of the source, based on the semi-major axis. 

With this parred down list we can then easily find pointings that overlap the FERMI ellipse is a fairly efficient manner. The missed area is then any part of the uncertainty ellipse which is not covered by the beam (here represented as a circle having a radius of the C-band HWHM). 

These areas are cataloged graphically in a picture database, as well as numerically in a catalog which we can use to develop a schedule for observations to cover these areas.