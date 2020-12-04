#!/usr/bin/env python2
import argparse
import numpy as np
import ROOT as r
import sys
# Fix https://root-forum.cern.ch/t/pyroot-hijacks-help/15207 :
r.PyConfig.IgnoreCommandLineOptions = True
import shipunit as u
import rootUtils as ut
from decorators import *
import shipRoot_conf
shipRoot_conf.configure()

#inputdir = "/eos/experiment/ship/user/dasukhon/"
r.gErrorIgnoreLevel = r.kWarning
r.gROOT.SetBatch(True)
parser = argparse.ArgumentParser(description='Script to create flux maps.')
parser.add_argument(
        'inputfile',
        help='''Simulation results to use as input. '''
        '''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
        'geofile',
        help='''Geometry file to use. '''
        '''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
        'outputfile',
        help='''The name of the output file. ''')
parser.add_argument(
        'outputfile2',
        help='''The name of the output file with all hits in the scoring plane. ''')
parser.add_argument(
        'outputfile3',
        help='''The name of the output file with muon hits in the scoring plane. ''')
args = parser.parse_args()
g = r.TFile.Open(args.geofile, 'read')
sGeo = g.FAIRGeom

ch = r.TChain("cbmsim")
ch.Add(args.inputfile)
print "Current input file ", args.inputfile
print "Output file with rates of at least one hit in SST ", args.outputfile
print "Hits in a scoring plane file ", args.outputfile2
print "Muon hits in a scoring plane file ", args.outputfile3
n = ch.GetEntries()
print "Number of entries ", n
i = 0
temp_weight = 1
nEvents_in_sst = 0
nHits_in_scoring_plane = 0
nHits_mu = 0

def viewSelect(x):
        return {
                0:'x1',
                1:'u',
                2: 'v',
                3:'x2',
        }.get(x,'x1')

for event in ch:
        if i % 10000 == 0:
                sys.stdout.write("Processing event:\t" + str(i) + "/" + str(n) + "\r" )
                sys.stdout.flush()

        muon = False
        muonid = 1

        i = i + 1
        hit_flag = False
        for hit in event.strawtubesPoint:
                if hit:
                        if not hit.GetEnergyLoss() > 0:
                                continue
#                        nEvents_in_sst += 1
                        px = hit.GetPx()
                        py = hit.GetPy()
                        pz = hit.GetPz()
                        pt = np.hypot(px, py)
                        P = np.hypot(pz, pt)
                        pid = hit.PdgCode()
                        tid = hit.GetTrackID()
                        if P > 0:
                                assert tid > 0
                                weight = event.MCTrack[tid].GetWeight()
                                detector_ID = hit.GetDetectorID()
                                station = detector_ID / 10000000
                                viewnb = (detector_ID - station * 10000000) / 1000000
                                view = viewSelect(viewnb)
                                plane = (detector_ID - station * 10000000 - viewnb * 1000000) / 100000
                                if station == 1 and view == 'x1' and plane == 0:
                                        nHits_in_scoring_plane += weight
                                if abs(pid) == 13:
                                        if station == 1 and view == 'x1' and plane == 0:
                                                nHits_mu += weight
                                if not hit_flag:
                                        nEvents_in_sst += weight
                                        hit_flag = True
                                if weight > temp_weight:
                                        temp_weight = weight

rates = open(args.outputfile,'a')
weights = open(args.outputfile2,"a")
muons = open(args.outputfile3,'a')
#
rates.write(str(nEvents_in_sst)+"\n")
weights.write(str(nHits_in_scoring_plane)+"\n")
muons.write(str(nHits_mu)+"\n")
#
rates.close()
weights.close()
muons.close()
