'This file is a brief Readme for the program 3D_gal/part_stack.py'


Overview:
	The two programs both use the same data source, in a way. They both
use the 100 halo sample from the Millennium Database that we have particle
data for, provided by Gerard Lemson in the summer of 2012. The gal program
uses galaxies and the part program uses particles. Both programs take the MPHI
technique that is a slight modification of the Caustic technique done by Dan
Gifford of U of M. The point of both programs, as evidenced by their titles,
is to take all 100 halos, bin them (join a certain number of them into a
catagory by similar variable values, such as mass), and join their
galaxy/particle data to create what is called an ensemble cluster, which is a
combination of a bunch of different halos' data. The assumption is, if the
halos are binned correctly then they are somewhat similar in dynamics, and
binning them isn't going to have too much adverse effects (scaling the data
before stacking should help lower this too). 
	After we create the ensemble clusters, I just apply the MPHI technique
to them. I can get a lot of various information using it, such as it's
accuracy to the table values (in Guo tables from MDB). 

	The purpose of the stacking is to see if we can reduce noise and
increase signal in our observations. In other words, stacking doesn't help our
resolution of a mass estimate for any particular halo, but if we can find
patterns for halos of certain sizes or luminosities, and stacking provides us
a way to assign a halo to this catagory using only a handful of galaxies, then
we can more closely estimate what kind of halo this might be.

	Another note: these programs are labeled '3D'. That is because the
velocity data is in 3D, meaning we use it's x,y,z values from the simulation.
In real life, we don't have x,y,z values, we only have projected values. 

	A note of interest, as of early 2013, we see that the MPHI technique
has an inherent SHAPE BIAS. In other words, we are now seeing that most halos
or galaxy clusters are elliptically shaped. However, the Caustic technique as
well as the MPHI technique assumes sphericity, which makes us underestimate
the mass by about 15-20%. It is known that stacking makes clusters more
spherical, at least in terms of position. This implies that stacking should
'fix' the shape bias, however, what we are REALLY looking for, is if stacking
affects the phase space, not the x vs y or ra vs dec plot. This is what I am
currently working on...





