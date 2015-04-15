#!/usr/bin/python3
"""
###########################################################################
# AUTOMATED SEQUENCING BY TWO-DIMENSIONAL LCMS ANALYSIS
#
# (C) 2015 Victor S. Lelyveld, Anders Bjorkbom, Jack Szostak, Massachusetts
# General Hospital
#
# This software is licensed under the GNU Affero General Public License
# version 3 (GNU AGPLv3). Contact for other licensing terms.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# REQUIREMENTS:
# (1) python 3+
# (2) numpy 1.8+
# (3) matplotlib 1.3+
# (4) statsmodels 0.6+ (statsmodels requires pandas, scipy, and patsy)
# (5) clustalo 1.2.1 ** (binary must be accessible in path)
#
# (**) Sievers F, et al.. Fast, scalable generation of high-quality protein
# multiple sequence alignments using Clustal Omega. Mol Syst Biol.
# 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. PMID: 21988835.
#
###########################################################################
"""


__author__ = 'lelyveld'

import sys
assert sys.version_info.major >= 3, 'Error: this package requires Python 3.'
from .core import *