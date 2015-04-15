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

import numpy as np
import os
from shutil import which
from time import time
from statsmodels.nonparametric.smoothers_lowess import lowess
from subprocess import call
from collections import Counter
import multiprocessing as mp
from itertools import chain
from copy import copy
import re
import csv

from .classes import isnum, compoundset, ppm2dm, adductset, baseset

# Determine processor/thread count to choose an optimal number of spawned processes
try:
    MAXTHREADS_ = mp.cpu_count()
except NotImplementedError:
    try:
        MAXTHREADS_ = os.cpu_count()
    except:
        MAXTHREADS_ = 4

# Global options
PARAMS_SET_ = False
PARAMS_ = {

    # ADDUCT CLUSTERING PARAMS
    'ADDUCTDB_CSV': None,
    'CLUSTER_PPM': None,
    'CLUSTER_RTWIDTH_FACTOR': None,
    'CLUSTER_MAX_DEPTH': None,

    # BASE MASS DB
    'BASEDB_CSV': None,

    # STATIC FILTERING PARAMS
    'FILTER_MIN_MASS': None,
    'FILTER_MIN_QUALITY': None,
    'FILTER_MIN_SCORE': None,
    'FILTER_MIN_INTENSITY': None,
    'FILTER_MIN_INTENSITYPCT': None,
    'FILTER_MIN_IONS': None,
    'FILTER_MIN_ZCOUNT': None,
    'FILTER_MIN_RT': None,
    'FILTER_MAX_RT': None,
    'FILTER_ORPHAN_WEIGHT': None,

    # LOCAL FITTING PARAMS
    'LOWESS_FRACTION': None,

    # RANDOM WALK PARAMS
    'WALK_DRT_MAX': None,
    'WALK_TRIALS': None,
    'WALK_STEPS_MIN_DRAFT': None,
    'WALK_STEPS_MIN_FINAL': None,
    'WALK_PPM': None,
    'WALK_INTERNAL_NOISE_VARIANCE': None,
    'WALK_STARTING_NOISE_VARIANCE': None,
    'WALK_TOPFRACTION': None,
    'WALK_BIDIR': None,

    # PLOTTING PARAMS
    'PLOT': None,
    'PLOT_SAVE': None,
    'PLOT_SAVE_FMT': None,
    'PLOT_LABELS': None,
    'PLOT_MASS_MIN': None,
    'PLOT_MASS_MAX': None,
    'PLOT_RT_MIN': None,
    'PLOT_RT_MAX': None,
    'PLOT_ALPHA_HIGH': None,
    'PLOT_ALPHA_LOW': None,
    'PLOT_ALL_WALKS': None,
    'PLOT_MARKERSIZE': None,
    'PLOT_USETEX': None
}

# Check for clustalo executable in path
# TODO: link against libclustalo instead of using the executable binary!
assert which('clustalo', mode=os.X_OK), 'Error: clustalo binary not found in path or inaccessible.'

def read_params(cfgfilename):
    """
    Import parameter csv file and load parameters in gerbil scope.

    Takes only a filename that must be a csv containing the names of all of the __PARAMS__ global dict keys.
    """
    global PARAMS_SET_
    global PARAMS_

    try:
        with open(cfgfilename) as csvfile:
            csvfile.seek(0)
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                if not len(row):
                    continue
                if "#" in row[0]:
                    continue
                for pkey in PARAMS_:
                    if row[0] == pkey:
                        val = row[1].split('#')[0].strip()
                        if len(val) <= 0 or val is None or val == 'None' or val == 'none' or val == '':
                            PARAMS_[pkey] = None
                        elif val == 'True' or val == 'true' or val is True:
                            PARAMS_[pkey] = True
                        elif val == 'False' or val == 'false' or val is False:
                            PARAMS_[pkey] = False
                        elif isnum(val):
                            PARAMS_[pkey] = float(val)
                        else:
                            PARAMS_[pkey] = val
    except:
        PARAMS_SET_ = False
        print('Error importing parameters.  Fail at everything.')
        return False

    assert os.path.isfile(PARAMS_['ADDUCTDB_CSV']), 'Error: adducts database csv file not found.'
    assert os.path.isfile(PARAMS_['BASEDB_CSV']), 'Error: base database csv file not found.'

    PARAMS_SET_ = True
    return True


def cluster_adducts(compounds, adducts):
    """
    Agglomerative adduct clustering


    Define a list of compounds in the compoundset that are related to one another by a known adduct mass. Assume that
    the parent compound is the compound on this list with the maximum integrated intensity ('Vol' key), and sum
    intensity from the other identified adduct compounds onto that of the parent compound. Return only a final list of
    parent compounds or those for which adducts were not identified.

    :compounds: list of compound dicts
    :adducts: adductset
    """

    donecompounds = []
    for c in compounds:

        # strictly enforce no double-counting of intensity (probably does nothing)
        if c in donecompounds:
            continue

        adducts_to_sum = recursive_cluster(c, compounds, adducts, depth=0)

        if c not in adducts_to_sum:
            adducts_to_sum.append(c)

        # choose parent compound as the one with maximum intensity
        parent = [p for p in adducts_to_sum if p['Vol'] == np.max([t['Vol'] for t in adducts_to_sum])][0]

        if len(adducts_to_sum) > 1:
            for t in adducts_to_sum:
                if t['Cpd'] == parent['Cpd']:
                    donecompounds.append(t)
                    continue
                parent['Vol'] = parent['Vol'] + t['Vol']
                parent['Vol %'] = parent['Vol %'] + t['Vol %']
                compounds.remove(t)
                donecompounds.append(t)

    return compounds


def recursive_cluster(c, compounds, adducts, depth=0):
    """
    Recursively cluster adducts on a parent compound, c, given a list of identified compounds and expected adducts to a
    specific recursion depth.

    :param c: compound dict containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    :param compounds: list of compound dicts containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    :param adducts: adductset
    :param depth: int
    :return: list of compound dicts
    """

    candidates_to_sum = []
    for a, act in zip(adducts.mass, adducts.action):
        if act == 0:
            continue
        result = findbymassrt(compounds, c['Mass'] + a, c['RT'],
                              PARAMS_['CLUSTER_PPM'],
                              PARAMS_['CLUSTER_RTWIDTH_FACTOR'])
        for r in result:
            if r['Cpd'] == c['Cpd']:
                # editing a list while iterating over it... shameful.  let's assume it works because we're always
                # removing from behind the cursor
                result.remove(r)
                continue
            else:
                candidates_to_sum.append(r)
                recompounds = [x for x in compounds if x != r]
                if depth < PARAMS_['CLUSTER_MAX_DEPTH']:
                    recurseresult_to_sum = recursive_cluster(r, recompounds, adducts, depth=depth + 1)
                    if recurseresult_to_sum:
                        candidates_to_sum.extend(recurseresult_to_sum)

    return list({item["Cpd"]: item for item in candidates_to_sum}.values())


def define_startingpos(compounds):
    """
    Define the likely starting point as the compound with maximum integrated intensity.
    (This function usually finds the full length compound, assuming that it has maximal abundance in the dataset. This
     is useful primarily for trying to guess the read orientation after walk generation.)

    :compounds: list of compound dicts containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    """

    # TODO: FOR THE TIME BEING, GUTTED ALL FEATURES THAT PERTAIN TO POOR INITIAL SAMPLE QUALITY

    maxintensity = np.max([c['Vol'] for c in compounds])
    pickedcompound = [c for c in compounds if c['Vol'] == maxintensity][0]

    # What if the sample is overdigested, and the starting material cannot be assumed to be one with highest intensity?
    #    Check to see if >20% of compound intensity is larger than the picked one.
    largercompounds = [c for c in compounds if c['Mass'] > pickedcompound['Mass']]
    largervol = np.sum([c['Vol'] for c in largercompounds])
    totalvol = np.sum([c['Vol'] for c in compounds])
    if largervol/totalvol > 0.2:
        # choose starting point as the one with maximum mass, to avoid truncating the dataset incorrectly on the mass
        # axis
        maxmass = np.max([c['Mass'] for c in compounds])
        pickedcompound = [c for c in compounds if c['Mass'] == maxmass][0]

    return pickedcompound


def findbymassrt(compounds, m, t, ppm, tstringency):
    """
    Find compounds near mass m and RT t within a mass ppm and RT stringency window

    :param compounds: list of compound dicts containing keys 'Mass', 'RT', and 'Width'
    :param m: float mass
    :param t: float retention time
    :return: list of compound dict matches
    """
    mstringency = ppm2dm(m, ppm)
    mcandidates = [c for c in compounds if ((c["Mass"] >= m - mstringency)
                                            & (c["Mass"] <= m + mstringency))]
    matches = [c for c in mcandidates if ((c["RT"] >= t - c["Width"] * tstringency)
                                          & (c["RT"] <= t + c["Width"] * tstringency))]
    return matches


def randstep(cpd, compounds, bases, dmmin, dmmax, dt, tdir, ppm, variance, start=None, internal=None, end=None):
    """
    Take a random step from the current compound, cpd, to another compound in the compoundset compounds.

    :param cpd: compound dict containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    :param compounds: list of compound dicts containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    :param bases: baseset
    :param dmmin: float
    :param dmmax: float
    :param dt: float
    :param tdir: [-1,0,1]
    :param variance: float
    :param start: logical
    :param internal: logical
    :param end: logical
    :return: compound dict
    """

    if tdir > 0:
        candidates = [c for c in compounds if ((c['Mass'] >= cpd['Mass'] - dmmax)
                                               & (c['Mass'] <= cpd['Mass'] - dmmin)
                                               & (c['RT'] >= cpd['RT'])
                                               & (c['RT'] <= cpd['RT'] + dt)
                                               & (c['Cpd'] != cpd['Cpd']))]
    elif tdir < 0:
        candidates = [c for c in compounds if ((c['Mass'] >= cpd['Mass'] - dmmax)
                                               & (c['Mass'] <= cpd['Mass'] - dmmin)
                                               & (c['RT'] >= cpd['RT'] - dt)
                                               & (c['RT'] <= cpd['RT'])
                                               & (c['Cpd'] != cpd['Cpd']))]
    else:
        candidates = [c for c in compounds if ((c['Mass'] >= cpd['Mass'] - dmmax)
                                               & (c['Mass'] <= cpd['Mass'] - dmmin)
                                               & (c['RT'] >= cpd['RT'] - dt)
                                               & (c['RT'] <= cpd['RT'] + dt)
                                               & (c['Cpd'] != cpd['Cpd']))]

    if len(candidates) >= 1:
        intensities = np.array([p['Vol'] for p in candidates], dtype=float)
        intensities /= np.max(intensities)
        noise = np.random.normal(loc=1.0, scale=variance, size=intensities.shape)
        prob = intensities + noise

        # Alternative algo for testing. This minimizes the difference between the intensity of adjacent steps
        # indiff = (np.array([p['Vol'] for p in candidates], dtype=float) - cpd['Vol'])/cpd['Vol']
        # noise = np.random.normal(loc=0.0, scale=variance, size=indiff.shape)
        # slopes = np.array([(p['RT'] - cpd['RT'])/(p['Mass'] - cpd['Mass']) for p in candidates], dtype=float) - tdir
        #
        # prob = abs(indiff + noise)
        # prob /= np.max(prob)
        # prob = 1.0/prob

        s = list(np.argsort(prob))

        basecall = 0
        while not basecall and len(s):
            i = s.pop()
            basecall = bases.findidbymass(cpd['Mass'] - candidates[i]['Mass'],
                                          ppm2dm(candidates[i]['Mass'], ppm),
                                          start=start,
                                          end=end,
                                          internal=internal)

        if basecall:
            picked = copy(candidates[i])
            picked['WalkScore'] = picked['Vol']  # define the score for this step as the next peak's intensity
            picked['Call'] = basecall
            return picked

    return {}


def append_full_length_to_walks(trials, compounds, dmmin, dmmax, startingpos, bases, ppm):
    """
    This function edits the trial walks to append the full length mass to the start of all walks, where possible.

    The full length compound is often not automatically included in the 5' ladder walks. The first nucleoside loss from
    the 3' end exposes one additional negative charge relative to the full length species, causing increased retention
    time in ion pairing mode. The result is that stepping from the full length compound to the first and second nearmost
    fragments on the authentic 5' ladder will violate the biased_walk rules that require all steps to be in the same
    retention time direction (enforced with the parameter tdir in those functions) after the first one in a walk. This
    fact is useful to determine which walk is the 3' ladder (i.e. it is the one for which the full length mass was
    included automatically). After this determination of orientation, the full length walk is correctly appended to the
    ladders, where a valid base call can be made, in a manner that ignores the tdir rules.

    :param trials: list of (list of compound dicts)
    :param compounds: list of compound dicts
    :param dmmin: float
    :param dmmax: float
    :param startingpos: compound dict
    :param bases: baseset
    :param ppm: float
    :return: list of (list of compound dicts)
    """
    print('Rechecking ' + str(len(trials)) + ' trial walks for full length step...')
    # see if i can make a step from to here from startingpos
    nearstartingpoint = np.array([c['Cpd'] for c in compounds if ((c['Mass'] >= startingpos['Mass'] - dmmax)
                                                                  & (c['Mass'] <= startingpos['Mass'] - dmmin))])

    for trial in trials:
        nextposcall = []
        if trial[0]['Cpd'] in nearstartingpoint:
            #mass difference between last step and this one
            mdiff = startingpos['Mass'] - trial[0]['Mass']

            #try to assign a base call
            nextposcall = bases.findidbymass(mdiff, ppm2dm(startingpos['Mass'], ppm), start=True, end=False)

        if nextposcall:
            trial[0]['Call'] = nextposcall  #append the base id to the position coordinates
            trial[0]['WalkScore'] += 0  #scoring here seems to be unhelpful
            trial.insert(0, startingpos)
            trial[0]['Call'] = 0  #first position has no call
            trial[0][
                'WalkScore'] = 0  #give starting point a starting score (this doesn't matter -- gets dropped later)

    return trials


def biased_walk(cpd, compounds, dmmin, dmmax, dt, ppm, maxtrials, minsteps, maxsteps, startingpos, bases,
                startingvar, internalvar, bidirectional=False):
    """
    Perform a biased rule-based random 2D walk along the compound set to generate a sequencing trajectory
    """

    trials = []
    for t in range(int(maxtrials / len(compounds))):

        # start unidirectional walk at each position
        pos = []
        pos.append(copy(cpd))
        pos[0]['WalkScore'] = pos[0]['Vol']  # assign initial score for this position as its intensity
        pos[0]['Call'] = 0

        # take a biased random unidirectional walk through the 2D (mass,RT) point space
        tdir = 0
        while len(pos) < maxsteps:
            nextposcall = ''

            # apply the same time direction for every subsequent step
            if tdir == 0 and pos[-1]['Cpd'] == startingpos['Cpd']:
                nextpos = randstep(pos[-1], compounds, bases, dmmin, dmmax, dt, tdir, ppm,
                                   variance=startingvar,
                                   start=True)
            else:
                nextpos = randstep(pos[-1], compounds, bases, dmmin, dmmax, dt, tdir, ppm,
                                   variance=internalvar,
                                   internal=True)

            if len(nextpos):
                # try to assign a base call
                nextposcall = nextpos['Call']

            if len(nextpos) and nextposcall:
                if not tdir and not bidirectional:
                    # get the RT direction of the first step (up or down) and use this for the remainder of the walk
                    tdir = int(np.sign(nextpos['RT'] - pos[-1]['RT']))

                    # FOR TESTING: calculate tdir as slope
                    # tdir = (nextpos['RT'] - pos[-1]['RT'])/(nextpos['Mass'] - pos[-1]['Mass'])

                # if we were able to step through the last position, then increment its score
                #pos[-1]['WalkScore'] += nextpos['WalkScore']  #uncomment to reward step-though

                # append the new step to the walk
                pos.append(nextpos)

                #take a stab at whether this is the last base in the sequence (because its exact mass is in the base db)
                nextposcall = bases.findidbymass(nextpos['Mass'],
                                                    ppm2dm(nextpos['Mass'], ppm),
                                                    start=False,
                                                    internal=False,
                                                    end=True)
                if nextposcall:
                    nextpos = copy(pos[-1])
                    nextpos['Call'] = nextposcall
                    nextpos['WalkScore'] = nextpos['Vol']

                    # FOR TESTING: update tdir slope
                    # tdir = (nextpos['RT'] - pos[-1]['RT'])/(nextpos['Mass'] - pos[-1]['Mass'])

                    #pos[-1]['WalkScore'] += nextpos['WalkScore'] #uncomment to reward step-through
                    pos.append(nextpos)
                    break

            else:
                # no possible step, so terminate the walk
                break

        if len(pos) > minsteps:
            trials.append(pos)

    #print('Finished walks for compound ' + str(cpd['Cpd']))
    return trials


def generate_trajectories(compounds, dt, ppm, maxtrials, minsteps, maxsteps, startingpos, bases,
                          startingvar, internalvar, bidirectional=False):
    """
    Generate a list of random walk trajectories through a compoundset based on a series of rules.

    This function requires the python 3 multiprocessing package.  Set the global __MAXTHREADS__ to alter the number
    of spawned processes.

    :param compounds: list of compound dicts
    :param dmmin: float
    :param dmmax: float
    :param dt: float
    :param maxtrials: int
    :param maxsteps: int
    :param startingpos: compound dict
    :param bases: baseset
    :param minsteps: int
    :param bidirectional:
    :return:
    """
    dmmin = float(np.min(bases.mass)-1.0)
    dmmax = float(np.max(bases.mass)+1.0)

    if os.name == 'nt':
        # Below is a single-threaded loop that does exactly what's done with starmap.  Windows is a pain with mp
        trials = []
        for cpd in compounds:
           trials.append(biased_walk(cpd, compounds, dmmin, dmmax, dt, ppm, maxtrials, minsteps, maxsteps,
                                     startingpos, bases, startingvar, internalvar, bidirectional))
    else:
        print('Using %d cores.' % MAXTHREADS_)
        with mp.Pool(MAXTHREADS_) as p:
            trials = p.starmap(biased_walk, [(cpd, compounds, dmmin, dmmax, dt, ppm, maxtrials, minsteps, maxsteps,
                                              startingpos, bases, startingvar, internalvar,
                                              bidirectional) for cpd in compounds])


    alltrials = []
    for trial in trials:
        alltrials.extend(trial)

    return alltrials


def classify_walk(walk, midline):
    """
    walk and midline are [[x,y]] arrays, and every x value in walk must be in midline
    return true if the majority of y values are greater than the midline
    """

    top = [(walk[i]['RT'] > midline[midline[:, 0] == walk[i]['Mass'], 1])[0] for i in range(len(walk))]
    return (np.count_nonzero(top) > len(top) / 2)


def getoffsets(s):
    s = ''.join(s)
    left = len(s) - len(s.lstrip('-'))
    right = len(s) - len(s.rstrip('-'))
    return left, right


def scored_read(alignment, walks, bases):
    """
    take a sequence alignment and its offsets from a set of walks, score each option at each position, and generate
    a final read from the highest scores the scoring information for a base in a walk is stored in each dict in walks
    """

    if alignment is None or walks is None or alignment == [] or not len(alignment):
        return [], np.array([])

    readlen = alignment.shape[1]
    finalread = list('-' * readlen)
    finalscore = np.zeros((readlen,), dtype=float)

    # scores = np.zeros_like(alignment, dtype=int)
    # for a,w,s in zip(alignment, walks, scores):
    #     l = getoffsets(a)[0]
    #     s[l:] = np.array(w['scores'][:-1])

    for i in range(readlen):
        call = []
        score = []

        for alignline, walk in zip(alignment, walks):
            loffset = getoffsets(alignline)[0]

            if i >= loffset and i - loffset < len(walk['calls']) - 1:
                call.append(bases.findnamebyid(walk['calls'][i - loffset + 1]))
                score.append(walk['scores'][i - loffset] + walk['scores'][i - loffset + 1])

        if not len(call):
            continue

        score = np.array(score)
        call = np.array(call)

        uniqueopt = np.unique(call)
        uniqueoptscore = []

        for o in uniqueopt:
            uniqueoptscore.append(np.sum(score[call == o]))

        topi = np.argsort(uniqueoptscore)[-1]
        finalposcall = uniqueopt[topi]
        finalread[i] = finalposcall
        finalscore[i] = uniqueoptscore[topi]

    return finalread, finalscore


def final_scored_read(alignment, scores=[]):
    """
    Generate a scored consensus read based on an alignment based on the frequency of the occurrence of each base at
    each position. We generally expect to run this on only two sequences (the result of each alignment from each
    orientation).  As such, the scores can be either 2, 1, or 0.
    """

    if alignment is None or alignment == [] or not len(alignment):
        return [], np.array([])

    finalreadlen = alignment.shape[1]
    finalread = list('-' * finalreadlen)
    finalscore = np.zeros((alignment.shape[1],), dtype=int)

    for i in range(finalreadlen):
        pos = ''.join([base for base in alignment[:, i] if base != '-'])
        options = Counter(pos).most_common()
        if len(options) == 1:
            finalread[i] = options[0][0]
            finalscore[i] = options[0][1]
        elif len(options) >= 2 and len(scores) >=2 and len(scores) == len(alignment):
            # if we have multiple options in this position, choose the one with top score
            offsets = [getoffsets(s) for s in alignment]
            topi = np.argsort([s[i - o[0]] for s, o in zip(scores, offsets)])[-1]
            finalread[i] = pos[topi]
            finalscore[i] = 0
        else:
            finalread[i] = 'X'
            finalscore[i] = 0

    return finalread, finalscore


def trials_to_stepdictlist(trials, orientations):
    """
    Can't remember why this is necessary right now...

    :param trials:
    :param orientations:
    :return:
    """

    walks = []
    for trial, orientation in zip(trials, orientations):
        steps = [t['Cpd'] for t in trial]
        scores = [t['WalkScore'] for t in trial]
        calls = [t['Call'] for t in trial if 'Call' in t]
        walk = {'steps': steps, 'calls': calls, 'scores': scores, 'orientation': orientation}
        walks.append(walk)

    return walks


def mungreads(reads):
    """
    Translate base single-letter codes in reads so that clustalo doesn't mung them up for us.
    """

    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    allreads = ''.join(reads)
    alpha = set(alphabet)
    used = set(allreads)
    remaining = alpha - used

    translation = {}
    for char in used - alpha:
        translation[char] = remaining.pop()

    for a, b in zip(translation.keys(), translation.values()):
        for i in range(len(reads)):
            reads[i] = re.sub(a, b, reads[i])

    return reads, translation


def unmungread(read, translation):
    """
    Reverse translation of munged base single-letter codes from read
    """

    for a, b in zip(translation.keys(), translation.values()):
        read = re.sub(b, a, read)
    return read


def align(reads):
    """
    Align a list of read strings, returning the aligned result as a numpy array of chars.

    :param reads: list of strings
    :return: numpy array
    """

    if len(reads) < 2 or reads is None or reads == '':
        print('Less than 2 reads, so alignment is futile.  Aborting alignment.')
        return np.array([list(r) for r in reads])

    #filter off empty lines
    reads = [r for r in reads if len(r) > 0]

    # Translate the reads into alpha single-letter code that is compatible with clustalo
    (reads, translations) = mungreads(reads)

    fname = 'reads' + str(time()).replace('.', '') + '.fasta'
    fasta = ''
    for i in range(len(reads)):
        fasta = fasta + '>' + str(i) + '\n' + reads[i] + '\n'

    # spit out fasta formatted reads
    try:
        f = open(fname, 'w')
        f.write(fasta)
        f.close()
    except IOError:
        print('Error: Writing reads to fasta file failed.')
        return None
    except:
        print('Error: Writing reads to fasta file failed.')
        raise

    try:
        sout = call('clustalo -i ' + fname + ' -o aligned.fasta --force -t rna', shell=True, stderr=sys.stdout)
    except IOError:
        print('Error: ClustalO failed.')
        return None
    except:
        print('Error: ClustalO failed.')
        raise

    #read back alignment from file
    alignment = []

    # We have a problem with how to deal with read-gapping.  The current solution is simply to warn.
    gapped = 0
    try:
        with open('aligned.fasta', 'r') as f:
            for line in f:
                if '>' in line:
                    i = int(line[1:-1])
                else:
                    l = line[0:-1]

                    if '-' in l.strip('-'):
                        gapped += 1

                    # Reverse translation is necessary because clustalo can't handle many non-alpha symbols
                    if len(translations):
                        l = unmungread(l, translations)

                    print(l)
                    alignment.append(np.array(list(l)))
        alignment = np.array(alignment)
    except:
        print('Error: ClustalO failed.')
        return None

    if len(alignment) < 1:
        print('Error: Reading ClustalO alignment output file failed.')
        return None

    if gapped > 0.2*len(reads):
        print('Warning: more than 20% of sequences are gapped, which implies low walk quality.')

    try:
        os.remove('aligned.fasta')
        os.remove(fname)
    except:
        print('Warning: Could not remove temp fasta files.')

    return alignment


def filterwalks(trials, orientations, minsteps, topfraction):
    """
    Cluster trial walks based on their starting point

    For each starting point, rank-order walks based on their total score.
    Choose a cut-off for normalized scores for each starting point, and discard those walks.
    Additionally, discard trial walks that are less then a step-length cutoff.

    minsteps:
    topfraction: for the set of walks sharing a single starting point, the top scoring fraction to retain [0.0 - 1.0]

    """

    trialstartpts = set([trial[0]['Cpd'] for trial in trials])
    filteredtrials = []
    filteredorientations = []

    for start in trialstartpts:
        samestart = [(trial, orientation) for trial, orientation in zip(trials, orientations) if
                     trial[0]['Cpd'] == start and len(trial) >= minsteps]

        cutoff = round(len(samestart) * topfraction)  # make this a parameter
        totalscores = [sum([step['WalkScore'] for step in trial[0]]) for trial in samestart]
        n = np.argsort(totalscores)[:cutoff]
        keeptrials = [samestart[i][0] for i in n]
        keeporientations = [samestart[i][1] for i in n]
        filteredtrials.extend(keeptrials)
        filteredorientations.extend(keeporientations)

    return filteredtrials, filteredorientations


def process(compounds_csvfile):

    assert os.path.isfile(compounds_csvfile), 'Error: compounds dataset csv file not found.'

    #if not PARAMS_SET_:
    #    assert read_params('default.cfg'), 'Error: parameters not set and default.cfg not found.'

    assert PARAMS_SET_, 'Error: cannot process dataset without first loading parameters, buddy.  Run read_params().'

    print('Processing compounds file: ' + compounds_csvfile)
    start_time = time()


    ## IMPORT SOURCE CSV FILES

    # define and import compoundset from csv file
    # compoundset.compounds is a list of compound dicts
    cdb = compoundset(compounds_csvfile)

    # define and import adducts from csv file
    adb = adductset(PARAMS_['ADDUCTDB_CSV'])

    # define and import base mass differences from csv file
    bdb = baseset(PARAMS_['BASEDB_CSV'])

    print('Found ' + str(len(cdb.compounds)) + ' compounds.')

    # Sort compoundset by abundance (i.e. intensity volume)
    cdb.sort('Vol')


    ## COMPOUND FILTERING

    # prefilter on quality and time-range before doing anything else (generally unused)
    if PARAMS_['FILTER_MIN_QUALITY'] is not None:
        cdb.filter(quality=(PARAMS_['FILTER_MIN_QUALITY'], 101))
    if PARAMS_['FILTER_MIN_RT'] is not None and PARAMS_['FILTER_MAX_RT']:
        cdb.filter(time=(PARAMS_['FILTER_MIN_RT'], PARAMS_['FILTER_MAX_RT']))

    # Agglomerative adduct clustering
    print('Cluster adducts, stringency: ppm = ' + str(PARAMS_['CLUSTER_PPM']) + ', RT peak width factor = ' +
          str(PARAMS_['CLUSTER_RTWIDTH_FACTOR']))

    cdb.compounds = cluster_adducts(cdb.compounds, adb)
    print('Adducts clustering leaves ' + str(len(cdb.compounds)) + ' compounds.')

    # Find candidate starting point for full length material
    # (Note: This is a starting point to later attempt to append to a final set of walks.  Walks from biased_walk
    # start from every compound point, not just this one.)
    startingpos = define_startingpos(cdb.compounds)
    print('Chose ' + str(startingpos['Mass']) + ' as candidate full length mass.')

    # This feature down-weights any compound for which no corresponding compound can be found that would sum to the
    # correct starting point mass (not necessary for most datasets, and may incidentally exclude valuable points
    # in the high mass region due to lost fragments in the low mass region).
    if PARAMS_['FILTER_ORPHAN_WEIGHT'] is not None and PARAMS_['FILTER_ORPHAN_WEIGHT'] > 0:
        cdb.weight_orphans(startingpos, PARAMS_['FILTER_ORPHAN_WEIGHT'], PARAMS_['WALK_PPM']*10)

    # Try to predict number of fragments using max intensity fragment and the average native base mass
    fragcount = np.round(((startingpos['Mass'] / 317) * 2 - 1))

    # apply the selected filters
    if PARAMS_['FILTER_MIN_MASS'] is not None:
        cdb.filter(mass=(PARAMS_['FILTER_MIN_MASS'], startingpos['Mass']))

    if PARAMS_['FILTER_MIN_INTENSITY'] is not None and PARAMS_['FILTER_MIN_INTENSITY'] > 0:
        cdb.filter(intensity=(PARAMS_['FILTER_MIN_INTENSITY'], 1.0E10),)

    if PARAMS_['FILTER_MIN_INTENSITYPCT'] is not None and PARAMS_['FILTER_MIN_INTENSITYPCT'] > 0:
        cdb.filter(intensitypct=(PARAMS_['FILTER_MIN_INTENSITYPCT'], 100))

    if PARAMS_['FILTER_MIN_IONS'] is not None and PARAMS_['FILTER_MIN_IONS'] > 0:
        cdb.filter(ions=(PARAMS_['FILTER_MIN_IONS'], 999))

    if PARAMS_['FILTER_MIN_ZCOUNT'] is not None and PARAMS_['FILTER_MIN_ZCOUNT'] > 0:
        cdb.filter(charges=(PARAMS_['FILTER_MIN_ZCOUNT'], 999))

    if PARAMS_['FILTER_MIN_SCORE'] is not None and PARAMS_['FILTER_MIN_SCORE'] > 0:
        cdb.filter(score=(PARAMS_['FILTER_MIN_SCORE'], 999))

    #if PARAMS_['FILTER_EXPFRAGMENT_FACTOR'] is not None and PARAMS_['FILTER_EXPFRAGMENT_FACTOR'] > 0:
    #    loosefragcount = int(fragcount * (1 + PARAMS_['FILTER_EXPFRAGMENT_FACTOR']))
    #    cdb.filter(fragments=loosefragcount)

    # Resort compoundset by abundance (i.e. intensity volume)
    cdb.sort('Vol')


    ## DEFINE MIDLINE FOR CLUSTERING SEQUENCING WALKS

    print('Fitting midline for clustering...')

    m = np.array([d['Mass'] for d in cdb.compounds])
    t = np.array([d['RT'] for d in cdb.compounds])

    # find a midline of the compound data points by LOWESS fitting
    midline = lowess(t, m, frac=PARAMS_['LOWESS_FRACTION'])


    ## GENERATE RANDOM SEQUENCE TRAJECTORIES

    print('Generating walk trajectories...')

    # try to trace contiguous fragment ladder by random walking through the data from a random starting point
    trials = generate_trajectories(cdb.compounds,
                                   dt=PARAMS_['WALK_DRT_MAX'],
                                   ppm=PARAMS_['WALK_PPM'],
                                   maxtrials=int(PARAMS_['WALK_TRIALS'] * len(cdb.compounds)),
                                   minsteps=int(PARAMS_['WALK_STEPS_MIN_DRAFT']),
                                   maxsteps=int(fragcount),
                                   startingpos=startingpos,
                                   bases=bdb,
                                   startingvar=PARAMS_['WALK_STARTING_NOISE_VARIANCE'],
                                   internalvar=PARAMS_['WALK_INTERNAL_NOISE_VARIANCE'],
                                   bidirectional=PARAMS_['WALK_BIDIR'])

    assert len(trials) >= 1, 'Error: No trials were generated.'

    # classify each trial by whether it is above or below the midline
    orientations = []
    for trial in trials:
        orientations.append(classify_walk(trial, midline))

    # the predicted maximum mass (the full length oligo) defines which cluster of trial walks is 5' -> 3'
    # orientation[i] is true if the trial is in the 5' -> 3' direction
    for trial, orientation in zip(trials, orientations):
        if trial[0]['Cpd'] == startingpos['Cpd'] and not orientation:
            orientations = [not o for o in orientations]
            break

    # since the full length species should be a component of both the 5' OH an 3' PO4 walks, append it to the start
    # of all trials (within ppm/base constraints).  NOTE: this didn't happen during seqwalk because inclusion of
    # the full length mass is the only clue to which cluster of walks is in the 3' -> 5' orientation.
    drafttrials = append_full_length_to_walks(trials,
                                              compounds=cdb.compounds,
                                              dmmin=min(bdb.mass)-1,
                                              dmmax=max(bdb.mass)+1,
                                              startingpos=startingpos,
                                              bases=bdb,
                                              ppm=PARAMS_['WALK_PPM'])

    # filter trial walks on the top fraction of walks at each starting point to retain and the final desired length of
    # walks (both input parameters)
    finaltrials, finalorientations = filterwalks(drafttrials,
                                                 orientations,
                                                 PARAMS_['WALK_STEPS_MIN_FINAL'],
                                                 PARAMS_['WALK_TOPFRACTION'])

    # recast the trials as a list of dicts:
    # each dict has a list of cpd steps ['steps'], scores for each step ['scores'], the based
    # call for each step given as a base id ['calls'], and the inferred read orientation.  this is an old object that
    # should eventually be replaced by the walkset object below.
    trialsteps = trials_to_stepdictlist(finaltrials, finalorientations)

    # A summary object that we might eventually make into a class to handle trajectories and resulting reads
    walkset = [{'trial': trial,
                'orientation': orientation,
                'calls': stepdict['calls'],
                'steps': stepdict['steps'],
                'scores': stepdict['scores'],
                'totalscore': sum([step['WalkScore'] for step in trial]),
                'read': ''.join([bdb.findnamebyid(i) for i in stepdict['calls'] if i > 0])
               } for trial, orientation, stepdict in zip(finaltrials, finalorientations, trialsteps)]

    # TODO: Coalesce walkset it into a class with its associated functions instead of the current mess


    ## SEQUENCE ALIGNMENT

    print('\n5\' ladder alignment:')
    fiveprime_reads = [w['read'] for w in walkset if w['orientation']]
    fiveprime_alignment = align(fiveprime_reads)
    (fiveprime_npread, fiveprime_score) = scored_read(fiveprime_alignment,
                                                      [u for u in walkset if u['orientation']],
                                                      bdb)
    fiveprime_scoredread = ''.join(fiveprime_npread)
    print('-'*len(fiveprime_scoredread) + '\n' + fiveprime_scoredread)

    print('\n3\' ladder alignment:')
    threeprime_reads = [w['read'] for w in walkset if not w['orientation']]
    threeprime_alignment = align(threeprime_reads)
    (threeprime_read, threeprime_score) = scored_read(threeprime_alignment,
                                                      [u for u in walkset if not u['orientation']],
                                                      bdb)
    threeprime_read = threeprime_read[::-1] #here is where we reverse the reads to be in the 5'-3' orientation
    threeprime_score = threeprime_score[::-1] #here is where we reverse the scores to be in the 5'-3' orientation
    threeprime_scoredread = ''.join(threeprime_read)
    print('-'*len(threeprime_scoredread) + '\n' + threeprime_scoredread[::-1])

    # Perform a final alignment
    print('\nFinal alignment:')
    alignment = align([fiveprime_scoredread, threeprime_scoredread])
    (finalcall, finalscore) = final_scored_read(alignment, scores=[fiveprime_score, threeprime_score])
    finalread = ''.join(finalcall)
    print(''.join([str(x) for x in finalscore]))
    print(finalread)


    if not finalcall or not len(finalcall):
        print('Nuts!  No final alignment returned!')
        if len(threeprime_scoredread):
            print('Using 3\' read as final.')
            finalread = threeprime_scoredread
        elif len(fiveprime_scoredread):
            print('Using 5\' read as final.')
            finalread = fiveprime_scoredread
        else:
            print('This hasn\'t work out well for you. Exiting.')
            elapsed_time = time() - start_time
            print('\nElapsed time: %.2f s' % elapsed_time)
            return ''


    # Filter out walks that generated reads that were eliminated by scoring/alignment (just for final display)
    if not PARAMS_['PLOT_ALL_WALKS']:
        walkset = [w for w in walkset if
                    ((w['read'] in finalread or w['read'] in fiveprime_scoredread) and w['orientation']) or
                    ((w['read'] in finalread[::-1] or w['read'] in threeprime_scoredread[::-1])
                     and not w['orientation'])]
        finaltrials = [w['trial'] for w in walkset]
        finalorientations = [w['orientation'] for w in walkset]

    elapsed_time = time() - start_time
    print('\nElapsed time: %.2f s' % elapsed_time)


    ## FINAL PLOTTING

    if PARAMS_['PLOT'] and len(finaltrials):

        # uniquify
        utrials = []
        uorientations = []
        for trial, orientation in zip(finaltrials, finalorientations):
            if trial not in utrials:
                utrials.append(trial)
                uorientations.append(orientation)

        niceplot(cdb.compounds,
                 utrials,
                 uorientations,
                 midline,
                 savefname=compounds_csvfile[0:compounds_csvfile.find('.csv')] + '_result',
                 save=PARAMS_['PLOT_SAVE'],
                 baselist=bdb,
                 saveformat=PARAMS_['PLOT_SAVE_FMT'],
                 title='Sequencing trajectories',
                 shownow=(not PARAMS_['PLOT_SAVE']),
                 msize=PARAMS_['PLOT_MARKERSIZE'],
                 plot_labels=PARAMS_['PLOT_LABELS'],
                 mlimit=(PARAMS_['PLOT_MASS_MIN'], PARAMS_['PLOT_MASS_MAX']),
                 tlimit=(PARAMS_['PLOT_RT_MIN'], PARAMS_['PLOT_RT_MAX']),
                 alphahigh=PARAMS_['PLOT_ALPHA_HIGH'],
                 alphalow=PARAMS_['PLOT_ALPHA_LOW'],
                 usetex=PARAMS_['PLOT_USETEX'])

    print('\nDone.')
    return finalread


def dodgetext(a,listodict,shift):
    for b in listodict:
        diff = np.abs(np.array(b['xytext']) - np.array(a['xytext']))
        if a['text'] == b['text'] and all(diff < np.array([0.5,0.5])):
            return None
        elif all(diff < np.array([100,1.0])):
            a['xytext'] = (a['xytext'][0], a['xytext'][1]+shift)
            a = dodgetext(a,listodict,shift)
            if a is None:
                return None
    return a

def niceplot(compounds,
             trials,
             orientations,
             midline,
             baselist,
             savefname='figure',
             save=False,
             cmap=None,
             msize=75,
             plot_midline=False,
             plot_labels=True,
             saveformat='pdf',
             shownow=True,
             title='',
             mlimit=(),
             tlimit=(),
             annotation_offset=(0.0, 2),
             fig=None,
             norm=None,
             alphalow=0.4,
             alphahigh=0.7,
             usetex=False):

    import matplotlib as mpl

    if mpl.get_backend() == 'pdf' or mpl.get_backend() == 'png':
        shownow = False
        save = True

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as colors

    # Set up matplotlib preferences
    mpl.rc('font', **{'family': 'sans-serif', 'size': 14})

    # Use latex rendering if flag is set
    if usetex:
        mpl.rc('text', usetex=True)
        mpl.rcParams['text.latex.unicode'] = True
        mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

    if cmap is None:
        cmap = cm.cool #cm.RdBu_r

    plt.ioff()

    v = np.array([d['Vol'] for d in compounds])
    vsort = np.argsort(v)
    v = v[vsort]
    m = np.array([d['Mass'] for d in compounds])[vsort]
    t = np.array([d['RT'] for d in compounds])[vsort]
    cpd = np.array([d['Cpd'] for d in compounds])[vsort]

    if fig == None:
        fig = plt.figure(figsize=(8, 6))
    if norm == None:
        norm = colors.LogNorm(vmax=np.max(v), clip=True)
    ax = fig.gca()

    if not len(mlimit) or None in mlimit:
        min_mass = 0 - 1.5 * annotation_offset[0]
        max_mass = np.max(m) + 350 + 1.5 * annotation_offset[0]
    else:
        min_mass = mlimit[0]
        max_mass = mlimit[1]

    if not len(tlimit) or None in tlimit:
        max_time = np.max(t) + 1.5 * annotation_offset[1]
        min_time = np.min(t) - 1.5 * annotation_offset[1]
    else:
        min_time = tlimit[0]
        max_time = tlimit[1]

    ax.set_xlim(min_mass, max_mass)
    ax.set_ylim(min_time, max_time)
    ax.set_xlabel('Mass (Da)')
    ax.set_ylabel('Retention Time (min)')
    ax.set_title(title)

    if len(trials) and len(orientations):
        p1 = plt.scatter(m, t, c=v, s=msize, linewidth=0, alpha=alphalow, cmap=cmap, norm=norm, zorder=1)

        # generate list of compounds that made it into final trials
        top_m = []
        top_v = []
        top_t = []
        top_c = []

        bottom_m = []
        bottom_v = []
        bottom_t = []
        bottom_c = []
        for c in cpd:
            for ft, i in zip(trials, range(len(trials))):
                if (c in [f['Cpd'] for f in ft]) and (c not in top_c):
                    p = [x for x in compounds if x['Cpd'] == c][0]
                    if orientations[i]:
                        top_m.append(p['Mass'])
                        top_v.append(p['Vol'])
                        top_t.append(p['RT'])
                        top_c.append(p['Cpd'])
                    else:
                        bottom_m.append(p['Mass'])
                        bottom_v.append(p['Vol'])
                        bottom_t.append(p['RT'])
                        bottom_c.append(p['Cpd'])
        if len(bottom_m):
            p4 = plt.scatter(bottom_m, bottom_t, c=bottom_v, s=msize, edgecolor='k', linewidth=1, marker='o',
                             alpha=alphahigh, cmap=cmap, norm=norm, zorder=3)
            cbar = fig.colorbar(p4)

        if len(top_m):
            p3 = plt.scatter(top_m, top_t, c=top_v, s=msize, edgecolor='k', linewidth=1, marker='s',
                             alpha=alphahigh, cmap=cmap, norm=norm, zorder=3)
            if 'cbar' in locals():
                cbar.vmin = np.min(np.hstack((top_v, bottom_v)))
                cbar.vmax = np.max(np.hstack((top_v, bottom_v)))
            else:
                cbar = fig.colorbar(p3)

        #plot trial walks
        for i in range(len(trials)):
            if orientations[i]:
                plt.plot([f['Mass'] for f in trials[i]], [f['RT'] for f in trials[i]], 'k-',
                         alpha=alphahigh, linewidth=1, zorder=2)
            else:
                plt.plot([f['Mass'] for f in trials[i]], [f['RT'] for f in trials[i]], 'k-',
                         alpha=alphahigh, linewidth=1, zorder=2)

        if plot_labels:
            ann1 = []
            ann0 = []
            for trial, orientation in zip(trials, orientations):
                for i in range(len(trial) - 1):
                    if orientation:
                        a = {'text': baselist.findnamebyid(trial[i + 1]['Call']),
                             'xy': (trial[i]['Mass'], trial[i]['RT']),
                             'xytext': (trial[i]['Mass'] / 2 + trial[i + 1]['Mass'] / 2 + annotation_offset[0],
                                        trial[i]['RT'] - annotation_offset[1]), 'color': trial[i + 1]['WalkScore']}
                        if a not in ann1:
                            a = dodgetext(a,ann1,-1)
                            if a is not None:
                                ann1.append(a)
                    else:
                        a = {'text': baselist.findnamebyid(trial[i + 1]['Call']),
                             'xy': (trial[i]['Mass'], trial[i]['RT']),
                             'xytext': (trial[i]['Mass'] / 2 + trial[i + 1]['Mass'] / 2 - annotation_offset[0],
                                        trial[i]['RT'] + annotation_offset[1]), 'color': trial[i + 1]['WalkScore']}
                        if a not in ann0:
                            a = dodgetext(a,ann0,1)
                            if a is not None:
                                ann0.append(a)
            ann = []
            for a in chain(ann0, ann1):
                ann.append(ax.annotate(a['text'], a['xy'], horizontalalignment='center', verticalalignment='center',
                                       textcoords='data', xytext=a['xytext'],
                                       arrowprops=dict(arrowstyle="-", color='#999999',
                                                       alpha=alphalow,
                                                       connectionstyle="angle,angleA=0,angleB=90,rad=0"),
                                       color='k'))

    elif len(trials):
        p1 = plt.scatter(m, t, c=v, s=msize, linewidth=0, alpha=alphahigh, cmap=cmap, norm=norm,
                         zorder=1)
        # plot trial walks
        for i in range(len(trials)):
            plt.plot([f['Mass'] for f in trials[i]], [f['RT'] for f in trials[i]], 'k-',
                     alpha=alphahigh, linewidth=1, zorder=2)
    else:
        p1 = plt.scatter(m, t, c=v, s=msize, linewidth=0, alpha=alphahigh, cmap=cmap, norm=norm, zorder=1)

    if 'cbar' not in locals():
        cbar = fig.colorbar(p1)

    if plot_midline and len(midline):
        p2 = plt.plot(midline[:, 0], midline[:, 1], 'k-.', zorder=1)

    cbar.set_label('Intensity (counts)', rotation=270, verticalalignment='bottom')

    plt.draw()

    if save:
        print('Saving plot to ' + savefname + '.' + saveformat)
        plt.savefig(savefname + '.' + saveformat, bbox_inches='tight', format=saveformat)

    if shownow:
        plt.show()

    plt.close()

    return fig, norm


def batch(compounds_dir='', compounds_fname='', repeats=1):
    """
    Testing function for batch operations on a whole directory with comparison of observed to predicted sequence.
    Correct sequences should be in a file with the same file stem name and the extension suffix ".known"

    :param compounds_dir: path string
    :param compounds_fname: file string to simply repeatedly run a single dataset
    :param repeats: number of repeats
    :return: a fractional error rate
    """
    totalscore = 0
    maxscore = 0

    if len(compounds_dir):
        dirlist = sorted(os.listdir(compounds_dir))
    else:
        dirlist = [compounds_fname]

    for r in range(repeats):
        for fname in dirlist:
            if re.search(r'\.csv$', fname):
                try:
                    # see if we have a known sequence input file
                    knownseq_fname = re.sub(r'\.csv$', '.known', fname)
                    if os.path.isfile(os.path.join(compounds_dir, knownseq_fname)):
                        with open(os.path.join(compounds_dir, knownseq_fname)) as f:
                            knownseq = f.read().strip()
                        maxscore = len(knownseq)
                    else:
                        knownseq = ''
                        maxscore = 0

                    print(fname)

                    seq = process(os.path.join(compounds_dir, fname))
                    print('\nExpected vs. Observed sequences:')
                    alignment = align([knownseq, seq])
                    comparison, score = final_scored_read(alignment)

                    incorrect = np.sum(np.array(score < 2, dtype=int))
                    print(''.join([str(x) for x in score]))
                    print('-------------------------------------------------------')
                    print('error rate: ' + str(incorrect/maxscore))
                    print('-------------------------------------------------------\n')
                except:
                    pass

    return None


if __name__ == "__main__":
    assert len(sys.argv) >= 3, 'Argument error. Run this module with two arguments:\n lcmsseq.py compounds.csv params.csv'

    compoundscsvfile = sys.argv[1]
    paramcsvfile = sys.argv[2]
    print(compoundscsvfile)
    print(paramcsvfile)

    assert os.path.isfile(paramcsvfile), 'Error: parameter csv file cannot be accessed.'
    assert read_params(paramcsvfile), 'Exit due to parameter import failure.'
    assert os.path.isfile(compoundscsvfile), 'Error: compound csv file cannot be accessed.'

    process(compoundscsvfile)




