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

import csv
import numpy as np
import re
import copy

_WATER_ = 18.010565

#monoisotopic masses
element = {'H':1.007825035,
           'D':2.014101779,
           '[2H]':2.014101779,
           'B':11.0093055,
           'C':12.0000000,
           '[13C]':13.00335483,
           'N':14.003074,
           '[15N]':15.00010897,
           'O':15.99491463,
           '[17O]':16.999131,
           '[18O]':17.9991603,
           'F':18.99840322,
           'P':30.973762,
           'S':31.9720707,
           'Cl':34.96885272,
           'As':74.9215942,
           'Br':78.9183361,
           'Se':79.9165196,
           'I':126.904473}


def isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def ppm2dm(m, ppm):
    """
    Convert ppm to mass delta
    Usage: ppm2dm(mass, ppm)
    """
    return m * ppm / 1.0E6


def dm2ppm(m, dm):
    """
    Convert mass delta to ppm
    Usage: dm2ppm(mass, DeltaMass)
    """
    return dm * 1.0E6 / m


# TODO: make a formula class
def str2formula(s):
    formuladict = {p[0]:p[1] for p in re.findall("([\+\-]\[\w+\]|[\+\-]\w|\[\w+\]|\w)([0-9]*)",s)}
    for k in formuladict:
        if len(formuladict[k]) < 1:
            formuladict[k] = 1
        else:
            formuladict[k] = int(formuladict[k])

    keys = [k for k in formuladict.keys() if "+" in k or "-" in k]
    for key in keys:
        atom = re.findall("(\w)",key)[0]
        if "+" in key:
            if atom in formuladict:
                formuladict[atom] += formuladict[key]
            else:
                formuladict[atom] = formuladict[key]
        if "-" in key:
            if atom in formuladict:
                formuladict[atom] -= formuladict[key]
            else:
                formuladict[atom] = -1*formuladict[key]
        formuladict.pop(key, None)

    return formuladict


def add_formulas(formula1,formula2):
    a = str2formula(formula2str(formula1))
    b = str2formula(formula2str(formula2))
    for key in a:
        if key in b:
            a[key] += b[key]
    b.update(a)
    return b


def sum_formulas(formuladictlist):
    fsum = {}
    for formula in formuladictlist:
        fsum.update(add_formulas(fsum,copy.copy(formula)))
    return fsum


def formula2str(formuladict):
    return ''.join([k+str(v) for k,v in zip(formuladict.keys(),formuladict.values())])


def formula2mass(formuladict):
    return sum(element[e]*formuladict[e] for e in formuladict)


class adductset:
    def __init__(self, adducts_fname):
        self.mass = np.array([], dtype=float)
        self.name = np.array([])
        self.action = np.array([], dtype=int)
        with open(adducts_fname) as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(2048))
            csvfile.seek(0)
            reader = csv.reader(csvfile, dialect)
            for row in reader:
                if "#" in row[0]:
                    continue
                self.name = np.append(self.name, str(row[0]).strip())
                self.mass = np.append(self.mass, float(row[1]))
                self.action = np.append(self.action, int(row[2]))

    def find(self, m, stringency):
        return (self.name[np.nonzero((self.mass >= m - stringency) & (self.mass <= m + stringency))],
                self.mass[np.nonzero((self.mass >= m - stringency) & (self.mass <= m + stringency))])

# TODO: reorganize this more logically as dicts or named tuples and see how performance is affected
class baseset:
    def __init__(self, bases_fname):
        self.name = np.array([])
        self.longname = np.array([])
        self.mass = np.array([], dtype=float)
        self.start = np.array([],dtype=bool)
        self.internal = np.array([],dtype=bool)
        self.end = np.array([],dtype=bool)
        self.formula = []
        self.id = np.array([], dtype=int)
        with open(bases_fname) as csvfile:
            #dialect = csv.Sniffer().sniff(csvfile.read(2048))
            csvfile.seek(0)
            reader = csv.reader(csvfile, delimiter=',') # dialect)
            for row in reader:
                if len(row) < 1:
                    continue
                if "#" in row[0]:
                    continue
                self.name = np.append(self.name, str(row[0]).strip())
                self.longname = np.append(self.longname, str(row[1]).strip())
                self.start = np.append(self.start, bool(int(row[2])))
                self.internal = np.append(self.internal, bool(int(row[3])))
                self.end = np.append(self.end, bool(int(row[4])))
                chemformula = str2formula(row[5])
                self.formula.append(chemformula)
                self.mass = np.append(self.mass, formula2mass(chemformula))
                self.id = np.append(self.id, self.mass.shape[0])

        print('Initialized Base Database')

    def findnamebymass(self, m, stringency, start=None, internal=None, end=None):

        type = np.ones_like(self.start)
        if start != None:
            type = type & (self.start == start)
        if internal != None:
            type = type & (self.internal == internal)
        if end != None:
            type = type & (self.end == end)

        found = (self.mass >= m - stringency) & (self.mass <= m + stringency) & type

        if True in found:
            return str(self.name[found][0]) #RETURNS ONLY FIRST RESULT
        else:
            return ''

    def findidbymass(self, m, stringency, start=None, internal=None, end=None):

        type = np.ones_like(self.start)
        if start != None:
            type = type & (self.start == start)
        if internal != None:
            type = type & (self.internal == internal)
        if end != None:
            type = type & (self.end == end)

        found = (self.mass >= m - stringency) & (self.mass <= m + stringency) & type

        if True in found:
            return self.id[found][0] #RETURNS ONLY FIRST RESULT
        else:
            return np.array([])

    def findmassbyname(self, name, start=None, internal=None, end=None):

        type = np.ones_like(self.start)
        if start != None:
            type = type & (self.start == start)
        if internal != None:
            type = type & (self.internal == internal)
        if end != None:
            type = type & (self.end == end)

        found = (self.name == name) & type

        if True in found:
            return self.mass[found][0] #RETURNS ONLY FIRST RESULT
        else:
            return np.array([])

    def findformulabyname(self, name, start=None, internal=None, end=None):

        type = np.ones_like(self.start)
        if start != None:
            type = type & (self.start == start)
        if internal != None:
            type = type & (self.internal == internal)
        if end != None:
            type = type & (self.end == end)

        found = (self.name == name) & type
        return [self.formula[i] for i in range(len(found)) if found[i]]

    def findnamebyid(self, id):
        if not any(self.id == id):
            return ' '
        else:
            return str(self.name[self.id == id][0])

    def findmassbyid(self, id):
        return(self.mass[self.id == id])

    def __len__(self):
        return len(self.mass)


class compoundset:
    """
    A compoundset is imported from a csv file generated by Agilent MassHunter Quantitative Analysis with default column
    naming.
    The compoundset contains a list of compounds generated by Agilent's Find by Molecular Feature algorithm and the
    properties of identified compounds.

    """
    def __init__(self, centroids_fname):
        ## NOTE: the compound object is a list of dicts with keys from the csv column names and values from each row
        self.compounds = []
        with open(centroids_fname) as csvfile:
            lastpos = 0

            #fudge to look for column headers by exclusion in the first 10 rows
            for i in range(10):
                lastpos = csvfile.tell()
                r = csvfile.readline()
                if len(r) > 20 and len(r) > r.count(',') + 20:
                    break
            csvfile.seek(lastpos)

            #dialect = csv.Sniffer().sniff(csvfile.read(2048))
            #csvfile.seek(0)

            reader = csv.DictReader(csvfile, fieldnames = [],
                                    restkey='undefined-fieldnames',
                                    delimiter=',',
                                    quotechar='"',
                                    lineterminator='\r\n') #dialect=dialect)
            headers = next(reader)
            reader.fieldnames = headers['undefined-fieldnames']
            for row in reader:
                for key in row.keys():
                    if 'undefined-fieldnames' in key:
                        continue
                    if row[key] == '':
                        row[key] = 0 #TESTING
                    elif row[key].isdigit():
                        row[key] = int(row[key])
                    elif isnum(row[key]):
                        row[key] = float(row[key])
                self.compounds.append(row)

    def findcompoundbyid(self, id):
        return([c for c in self.compounds if c['Cpd'] == id])

    def sort(self, key):
        v = np.array([c[key] for c in self.compounds])
        vsort = np.argsort(v)[::-1]
        return [self.compounds[i] for i in vsort]

    def filter(self, mass=None, quality=None, intensity=None, intensitypct=None, ions=None, charges=None,
               fragments=None, time=None, score=None):

        if time and time != (None,None):
            self.compounds = self.time_bp(self.compounds, time[0], time[1])
            n = len(self.compounds)
            print('Retention Time bp ' + str(time) + ', leaves ' + str(n) + ' compounds.')

        if intensity:
            self.compounds = self.intensity_bp(self.compounds, intensity[0], intensity[1])
            n = len(self.compounds)
            print('Intensity bp ' + str(intensity) + ', leaves ' + str(n) + ' compounds.')

        if intensitypct:
            self.compounds = self.intensitypct_bp(self.compounds, intensitypct[0], intensitypct[1])
            n = len(self.compounds)
            print('Intensity % bp ' + str(intensitypct) + ', leaves ' + str(n) + ' compounds.')

        if mass:
            self.compounds = self.mass_bp(self.compounds, mass[0], mass[1])
            n = len(self.compounds)
            print('Mass bp ' + str(mass) + ', leaves ' + str(n) + ' compounds.')

        if quality:
            self.compounds = self.quality_bp(self.compounds, quality[0], quality[1])
            n = len(self.compounds)
            print('Quality bp ' + str(quality) + ', leaves ' + str(n) + ' compounds.')

        if score:
            self.compounds = self.score_bp(self.compounds, score[0], score[1])
            n = len(self.compounds)
            print('Score bp ' + str(score) + ', leaves ' + str(n) + ' compounds.')

        if ions:
            self.compounds = self.ions_bp(self.compounds, ions[0], ions[1])
            n = len(self.compounds)
            print('Ions bp ' + str(ions) + ', leaves ' + str(n) + ' compounds.')

        if charges:
            self.compounds = self.charges_bp(self.compounds, charges[0], charges[1])
            n = len(self.compounds)
            print('Charges bp ' + str(charges) + ', leaves ' + str(n) + ' compounds.')

        if fragments:
            if not intensity:
                maxintensity = 1E10
            else:
                maxintensity = intensity[1]
            self.compounds = self.thresholdbyfragcount(self.compounds, fragments, maxintensity)
            n = len(self.compounds)
            print('Fragments ' + str(fragments) + ', leaves ' + str(n) + ' compounds.')

    def mass_bp(self, compounds, low, high):
        "filter out compounds with mass below low and above high"
        assert 'Mass' in compounds[-1], 'Mass column missing from compound dataset.'
        highpassed = [d for d in compounds if (d["Mass"] >= low)]
        bandpassed = [d for d in highpassed if (d["Mass"] <= high)]
        return bandpassed

    def time_bp(self, compounds, low, high):
        "filter out compounds with retention time below low and above high"
        assert 'RT' in compounds[-1], 'RT column missing from compound dataset.'
        highpassed = [d for d in compounds if (d["RT"] >= low)]
        bandpassed = [d for d in highpassed if (d["RT"] <= high)]
        return bandpassed

    def quality_bp(self, compounds, low, high):
        "filter out compounds with quality below low and above high"
        if 'Quality Score' not in compounds[-1]:
            return compounds
        highpassed = [d for d in compounds if (d["Quality Score"] >= low)]
        bandpassed = [d for d in highpassed if (d["Quality Score"] <= high)]
        return bandpassed

    def score_bp(self, compounds, low, high):
        "filter out compounds with score (FBF) below low and above high"
        if 'Score' not in compounds[-1]:
            return compounds
        highpassed = [d for d in compounds if (d["Score"] >= low)]
        bandpassed = [d for d in highpassed if (d["Score"] <= high)]
        return bandpassed

    def intensity_bp(self, compounds, low, high):
        "filter out compounds with intensity below low"
        assert 'Vol' in compounds[-1], 'Vol column missing from compound dataset.'
        highpassed = [d for d in compounds if (d["Vol"] >= low)]
        bandpassed = [d for d in highpassed if (d["Vol"] <= high)]
        return bandpassed

    def ions_bp(self, compounds, low, high):
        "filter out compounds with ions below low"
        if 'Ions' not in compounds[-1]:
            return compounds
        highpassed = [d for d in compounds if (d["Ions"] >= low)]
        bandpassed = [d for d in highpassed if (d["Ions"] <= high)]
        return bandpassed

    def charges_bp(self, compounds, low, high):
        "filter out compounds with ions below low"
        if 'Z Count' not in compounds[-1]:
            return compounds
        highpassed = [d for d in compounds if (d["Z Count"] >= low)]
        bandpassed = [d for d in highpassed if (d["Z Count"] <= high)]
        return bandpassed

    def intensitypct_bp(self, compounds, low, high):
        "filter out compounds with intensity below low"
        if 'Vol %' not in compounds[-1]:
            return compounds
        highpassed = [d for d in compounds if (d["Vol %"] >= low)]
        bandpassed = [d for d in highpassed if (d["Vol %"] <= high)]
        return bandpassed

    def thresholdbyfragcount(self, compounds, fragments, maxintensity):
        threshold = 0
        step = 1000
        count = len(compounds)
        while (threshold < maxintensity) & (count > fragments):
            threshold = threshold + step
            compounds = self.intensity_bp(compounds, threshold, maxintensity)
            count = len(compounds)
        return compounds

    def findbymass(self, m, stringency):
        return [c for c in self.compounds if (c['Mass'] >= m - stringency) and (c['Mass'] <= m + stringency)]

    def weight_orphans(self, startingpos, factor, ppm):
        """
        Iterate through compounds, downweight compounds by factor if no other mass can be found that sums to the
        predicted starting material mass.

        :param compounds:
        :param startingpos:
        :return:
        """

        for c in self.compounds:
            if c == startingpos:
                continue
            diff = startingpos['Mass'] - c['Mass'] + _WATER_
            if not len(self.findbymass(diff,ppm2dm(diff,ppm) )):
                c['Vol'] = c['Vol']*factor

