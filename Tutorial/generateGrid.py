import MoogSymphony
import sys

"""
This program should be run first in the tutorial sequence.  It generates a grid
of raw synthetic spectra, and saves them for further processing

generateGrid.py must be called with two command line arguments.  The first gives
the name of the configuration file, and the second specifies the MoogStokesPy shared
object (ALPHA, BRAVO, CHARLIE, DELTA...), in the event that you would like to run
more than one instance of generateGrid (useful on machines with more than one core)

1) Edit the configuration file to reflect the grid you would like to generate.
For this tutorial, I have included an example configuration file with all the
necessary parameters.  It is named generateGrid.cfg

2) Start ipython

> ipython

3) Run the program

> run generateGrid.py generateGrid.cfg ALPHA

"""

configfile = sys.argv[1]
MoogInstance = sys.argv[2]

composition = MoogSymphony.Symphony(configfile, MoogInstance)

composition.compose()

