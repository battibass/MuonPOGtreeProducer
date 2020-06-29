"""
A macro to convert an event/object based muon ROOT tree
into a PANDAS data frame with one row per muon saved as CSV
file.
"""

import math
import argparse
import sys

import numpy as np
import pandas as pd

import ROOT as root


##### Helper Function ####################
def computePull(trk, seg, trkErr, segErr):
    """
    Retruns pull out of track and segment
    positions and their errors
    """
    return (trk - seg) / math.sqrt(trkErr*trkErr + segErr*segErr)



##### General constants ##################
DR_CUT = 0.1

NAN_INT = 0
NAN_FLOAT = 30.


##### Input arguments ####################
PARSER = argparse.ArgumentParser(description=__doc__)

PARSER.add_argument("fileNamesInput",
                    help="Input Muon POG ROOT tree files. \
                          Multiple files are given as string separated \
                          with blanks, e.g.: \"filename1 filename2\".")

PARSER.add_argument("fileNameOutput",
                    help="Path to the output Pandas DataFrame file.")

PARSER.add_argument("-n", "--nEvents",
                    default=-1,
                    type=int,
                    help="The number of events to process (-1 means \
                          all events).")

PARSER.add_argument("-g", "--gen",
                    default="p pt eta charge qOverPt",
                    help="List of GEN-related variables to be \
                          migrated from the ntuples to the Pandas \
                          DataFrame. Multiple variables are separated \
                          with blanks, e.g.: \"pt eta\".")

PARSER.add_argument("-m", "--mu",
                    default="pt eta isHighPt tunePTrackType dxy dz",
                    help="List of muon object variables to be \
                          migrated from the ntuples to the Pandas \
                          DataFrame. Multiple variables are separated \
                          with blanks, e.g.: \"pt eta\".")

PARSER.add_argument("-r", "--refit",
                    default="pt eta charge chi2 qOverPt ptErrOverPt "
                    + "nPixHits nTrkLays nMuHits",
                    help="List of muon refit variables to be \
                          migrated from the ntuples to the Pandas \
                          DataFrame. Multiple variables are separated \
                          with blanks, e.g.: \"pt eta\".")

ARGS = PARSER.parse_args()

root.gROOT.LoadMacro("interface/MuonPogTree.h++")
from ROOT.muon_pog import Event, MuonFitType, MuonDetType


##### Quantities for DataFrame ##########
GEN_QUANTITIES = ARGS.gen.split(" ")
MU_QUANTITIES = ARGS.mu.split(" ")
MU_REFITS = {"tuneP" : MuonFitType.TUNEP,
             "picky" : MuonFitType.PICKY,
             "dyt"   : MuonFitType.DYT,
             "tpfms" : MuonFitType.TPFMS,
             "inner" : MuonFitType.INNER}
MU_REFIT_QUANTITIES = ARGS.refit.split(" ")
MU_CHAMB_QUANTITIES = ["nDigis", "pull"] # CB : this is hardcoded, they
                                         # are not a simple variable
                                         # translation, they are instead
                                         # filled using a more complex
                                         # logic
NP_ARRAYS = {}

for q in GEN_QUANTITIES:
    name = "gen_" + q
    NP_ARRAYS[name] = np.array([])

for q in MU_QUANTITIES:
    name = "mu_" + q
    NP_ARRAYS[name] = np.array([])

for rName, r in MU_REFITS.iteritems():
    for q in MU_REFIT_QUANTITIES:
        name = "mu_" + rName + "_" + q
        NP_ARRAYS[name] = np.array([])

for q in MU_CHAMB_QUANTITIES:
    for st in range(1, 5):
        for detector in ["DT", "CSC"]:
            name = "mu_" + q + "_" + detector + "_st" + str(st)
            NP_ARRAYS[name] = np.array([])


##### ROOT TChain ########################
MUON_POG_TREE = root.TChain("MuonPogTree/MUONPOGTREE")

for fileName in ARGS.fileNamesInput.split(" "):
    print "[" + __file__ + "] Appending ROOT file {:} to TChain".format(fileName)
    MUON_POG_TREE.AddFile(fileName)

TOT_ENTRIES = MUON_POG_TREE.GetEntries()
READ_ENTRIES = 0

print "[" + __file__ + "] Entries to be processed : {:9d} ".format(TOT_ENTRIES)


##### Main loop on TChain ################
N_MUONS = 0

for entry in MUON_POG_TREE:

    READ_ENTRIES += 1

    if READ_ENTRIES % 1000 == 0:
        fracTot = float(READ_ENTRIES) / TOT_ENTRIES * 100.
        print "[" + __file__ + "] Processed {:9d} entries ({:4.2f}% of total)".format(READ_ENTRIES,
                                                                                      fracTot)

    if ARGS.nEvents >= 0 and READ_ENTRIES > ARGS.nEvents:
        break

    for genPart in entry.event.genParticles:

        if genPart.status != 1 or abs(genPart.pdgId) != 13:
            continue

        genMuVector = root.TLorentzVector()
        genMuVector.SetPtEtaPhiM(genPart.pt,
                                 genPart.eta,
                                 genPart.phi,
                                 0.106)

        minDr = 999.
        iBestMuon = -1

        nMuonsPerEv = len(entry.event.muons)

        for iMu, mu in zip(range(0, nMuonsPerEv),
                           entry.event.muons):

            muVector = root.TLorentzVector()
            muVector.SetPtEtaPhiM(mu.pt,
                                  mu.eta,
                                  mu.phi,
                                  0.106)

            dR = genMuVector.DeltaR(muVector)

            if not(mu.isTracker and mu.isGlobal):
                continue

            if dR < DR_CUT and dR < minDr:
                minDr = dR
                iBestMuon = iMu

        if iBestMuon < 0:
            continue

        N_MUONS += 1

        mu = entry.event.muons[iBestMuon]

        ##### GEN quantities #####
        for q in GEN_QUANTITIES:
            name = "gen_" + q
            if q == "p":
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], genMuVector.P())
            elif q == "charge":
                charge = - genPart.pdgId / abs(genPart.pdgId)
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], charge)
            elif q == "qOverPt":
                charge = - genPart.pdgId / abs(genPart.pdgId)
                qOverPt = charge / genPart.pt
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], qOverPt)
            else:
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], genPart.__getattribute__(q))

        ##### DIGI-based shower info #####
        digisPerSt = {"DT"  :[NAN_INT, NAN_INT, NAN_INT, NAN_INT],
                      "CSC" :[NAN_INT, NAN_INT, NAN_INT, NAN_INT]}

        for match in mu.matches:
            st = abs(match.id_r)
            chambType = match.type
            if chambType == MuonDetType.DT or \
               chambType == MuonDetType.CSC:
                detector = ("DT" if chambType == MuonDetType.DT else "CSC")
                nDigis = min(match.nDigis, 300) # CB limited to 300, should suffice
                if digisPerSt[detector][st - 1] < nDigis:
                    digisPerSt[detector][st - 1] = nDigis

        for detector, digisPerDetector in digisPerSt.iteritems():
            for st, nDigis in zip(range(1, 5), digisPerDetector):
                name = "mu_nDigis_" + detector + "_st" + str(st)
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], nDigis)

        ##### Track - segment pulls #####
        pullsPerSt = {"DT"  : [NAN_FLOAT, NAN_FLOAT, NAN_FLOAT, NAN_FLOAT],
                      "CSC" : [NAN_FLOAT, NAN_FLOAT, NAN_FLOAT, NAN_FLOAT]}

        for match in mu.matches:
            st = abs(match.id_r)
            chambType = match.type
            xTrk = match.x
            xTrkErr = match.errx
            xSeg = NAN_FLOAT
            xSegErr = NAN_FLOAT
            for iSeg, qual in zip(match.indexes, match.matchQuals):
                if qual.test(1): # CB arbitrated match for tracker muon
                    if abs(xSeg - NAN_FLOAT) > 0.01 and \
                       abs(xSegErr - NAN_FLOAT) > 0.01:
                        print "[" + __file__ + "] Error: more than one arbitrated match per station"
                    if chambType == MuonDetType.DT:
                        segment = entry.event.dtSegments[iSeg]
                        detector = "DT"
                    elif chambType == MuonDetType.CSC:
                        segment = entry.event.cscSegments[iSeg]
                        detector = "CSC"
                    else:
                        print "[" + __file__ + "] Error: arbitrated match is not of DT or CSC type"
                        sys.exit(999)

                    xSeg = segment.x
                    xSegErr = segment.errx

                    pull = computePull(xTrk, xSeg, xTrkErr, xSegErr)
                    pullsPerSt[detector][st -1] = pull

        for detector, pullsPerDetector in pullsPerSt.iteritems():
            for st, pull in zip(range(1, 5), pullsPerDetector):
                name = "mu_pull_" + detector + "_st" + str(st)
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], pull)

        ##### MU quantities (general ones) #####
        for q in MU_QUANTITIES:
            name = "mu_" + q
            if q == "qOverPt":
                qOverPt = mu.charge / mu.pt
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], qOverPt)
            else:
                NP_ARRAYS[name] = np.append(NP_ARRAYS[name], mu.__getattribute__(q))

        ##### MU refit quantities #####
        for rName, r in MU_REFITS.iteritems():
            for q in MU_REFIT_QUANTITIES:
                name = "mu_" + rName + "_" + q
                refit = mu.fits[r]
                if q == "qOverPt":
                    qOverPt = refit.charge / refit.pt
                    NP_ARRAYS[name] = np.append(NP_ARRAYS[name], qOverPt)
                elif q == "ptErrOverPt":
                    ptErrOverPt = refit.ptErr / refit.pt
                    NP_ARRAYS[name] = np.append(NP_ARRAYS[name], ptErrOverPt)
                else:
                    NP_ARRAYS[name] = np.append(NP_ARRAYS[name], refit.__getattribute__(q))

##### Pandas DataFrame ###################
DATA_F = pd.DataFrame(data=NP_ARRAYS,
                      index=np.arange(N_MUONS))

DATA_F.to_csv(ARGS.fileNameOutput, index=False)
