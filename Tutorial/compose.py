import MoogSymphony
import sys

configfile = sys.argv[1]
MoogInstance = sys.argv[2]

composition = MoogSymphony.Symphony(configfile, MoogInstance)

composition.compose()

