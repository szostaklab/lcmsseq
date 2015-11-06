#!/usr/bin/python3

__author__ = 'lelyveld'

# Run this file to regenerate the manuscript figures, as indicated

import lcmsseq


#walk at 10 ppm
fig3a = 'fig3a.csv'
lcmsseq.read_params('fig3a.cfg')
lcmsseq.process(fig3a)

#walk at 5 ppm
fig3b = 'fig3b.csv'
lcmsseq.read_params('default.cfg')
lcmsseq.process(fig3b)

#walk at 5 ppm and allow 7 min step for correct orientation
fig3c = 'fig3c.csv'
lcmsseq.read_params('fig3c.cfg')
lcmsseq.process(fig3c)

#RT peak width factor at 5 to capture CO adduct
fig3d = 'fig3d.csv'
lcmsseq.read_params('fig3d.cfg')
lcmsseq.process(fig3d)

# walk at 5 ppm (same data as compounds.csv)
fig1c = 'fig1c.csv'
lcmsseq.read_params('default.cfg')
lcmsseq.process(fig1c)

