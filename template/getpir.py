import sys
import os

from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='PROGETTO')
aln.append_model (mdl, align_codes='wt_PROGETTO',atom_files='PROGETTO.pdb')
aln.write(file='PROGETTO.pir', alignment_format='PIR')

