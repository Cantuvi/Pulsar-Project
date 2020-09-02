#!/usr/bin/env python
import os,sys, argparse
import numpy as np
from numpy.random import default_rng
from astropy import units as u
from astropy.coordinates import SkyCoord

class Observation:
    def __init__(self, line):
        splitline = line.split()
        self.date = float(splitline[0])
        self.rauncertainty = float(splitline[2]) # This is in "seconds", which is arcseconds / (15 * cos(declination)
        self.decuncertainty = float(splitline[4]) # This is in arcseconds
        self.position = SkyCoord(splitline[1], splitline[3], frame='icrs', unit=(u.hourangle, u.deg))

    def perturbposition(self, deltaramas, deltadecmas):
        self.position = SkyCoord(self.position.ra + deltaramas*u.mas, self.position.dec + deltadecmas * u.mas)

    def addUncertainty(self, rauncertaintymas, decuncertaintymas):
        self.rauncertainty = np.sqrt(self.rauncertainty**2 + (rauncertaintymas/(15000.0*np.cos(self.position.dec.value)))**2)
        self.decuncertainty = np.sqrt(self.decuncertainty**2 + (decuncertaintymas/1000.0)**2)

    def to_string(self):
        rastring = self.position.to_string(decimal=False,sep=':',unit=u.hourangle,pad=True,precision=7).split()[0]
        decstring = self.position.to_string(decimal=False,sep=':',unit=u.deg, pad=True,precision=7).split()[1]
        return("{0:0.4f} {1} {2} {3} {4}".format(self.date, rastring, self.rauncertainty, decstring, self.decuncertainty))

if __name__ == "__main__":
    # Get some info on what we are expected to do
    parser = argparse.ArgumentParser(description='Read an existing pmpar file and perturb the positions.')
    parser.add_argument('--distribution', default="gaussian", help='Distribution type of the uncertainty to add')
    parser.add_argument('--extent', default=1.0, type=float, help='Extent or std dev of the uncertainty distribution in milliarcseconds')
    parser.add_argument('pmparfile', default="", type=str, nargs=1)

    # Parse the arguments
    args = parser.parse_args()

    # Open the pmpar file
    lines = open(args.pmparfile[0]).readlines()

    # Create a list where we will store the observations
    obslist = []

    # Create a list to store the lines that are not observation lines from the pmpar file
    otherlines = []

    # read in the lines
    for line in lines:
        # Check if it is an observation (length 5)
        if len(line.split()) == 5:
            obslist.append(Observation(line))
        else:
            otherlines.append(line)

    # Now go through the observations with a for loop, adding a random offset to each one using the "addUncertainty" method
    rng = default_rng() # rng is your random number generator.  See e.g. https://numpy.org/doc/stable/reference/random/index.html for how to use it.
    # This part is for you

    # Then write the result back out - first all the "otherlines", then each observation, using the to_string() method
    perturbedout = open(args.pmparfile[0] + ".perturbed", "w")
    for line in otherlines:
        perturbedout.write(line)
    for obs in obslist:
        perturbedout.write("{0}\n".format(obs.to_string()))
    perturbedout.close()
