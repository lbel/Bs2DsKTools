#!/usr/bin/env python
# coding=utf-8

import argparse
#from utils import _import_args
def _import_args(namespace, d = {}):
    """Import attributes from namespace to local environment.

    namespace -- namespace to import attributes from
    d         -- dictionary that is returned with attributes
                 and values (default: empty dict, leave it
                 this way unless you know what you are doing)

    Usage:
      >>> opts = parser.parse_args(['foo', '-o', 'bar'])
      >>> locals().update(_import_args(opts))

    """
    attrs = vars(namespace)
    for attr in attrs:
        d[attr] = getattr(namespace, attr)
    return d


optparser = argparse.ArgumentParser(description=__doc__)
optparser.add_argument('filename', help='ROOT file')
optparser.add_argument('-s', dest='session', help='Session name')
optparser.add_argument('-o', dest='out', help='ROOT file with output histograms')
optparser.add_argument('-n', dest='name', help='Input tree name')
options = optparser.parse_args()
locals().update(_import_args(options))

import sys, os
if not os.path.exists(filename):
    sys.exit('File not found: {}'.format(filename))

from array import array

from ROOT import gROOT
gROOT.SetBatch(True)

from ROOT import gDirectory, gSystem, gPad, gStyle
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TCanvas, TPad

from ROOT import TMVA
from tmvaconfig import TMVAconfig, ConfigFile

# instantiate
TMVA.Tools.Instance()

# ownership
TFile.Open._creates = True

# files
ifile = TFile.Open(filename, 'read')
ofile = TFile.Open(out, 'recreate')

# reader
reader = TMVA.Reader('!Color:!Silent')

# read config
conf = ConfigFile('TMVA.conf')
if conf.read() > 0:
    session = conf.get_session_config(session)
print '::: Applying {} MVAs: {}\n{}'.format(len(session.methods),
                                            session.methods, '='*50)
print session
print ':::'

# trained variables
allvars = []
for var in session.all_vars():
    allvars += [array('f', [0.])]
    reader.AddVariable(var, allvars[-1])

# prepare apply tree
itree = ifile.Get(name)
#for i, var in enumerate(session.vars):
    # combined vars at the end (ignored)
    #itree.SetBranchAddress(var, allvars[i])

# prepare output tree
ofile.cd()
otree = itree.CloneTree(0)
for i, var in enumerate(session.all_vars()):
    # variables for output tree
    otree.Branch(var, allvars[i], '{}/F'.format(var))

# spectators
spectators = []
for var in session.spectators:
    spectators += [array('f', [0.])]
    reader.AddSpectator(var, spectators[-1])

# book methods
hists_discr = {}
discriminant = {}
for method in session.methods:
    reader.BookMVA(method, '{0}/weights/{0}_{1}.weights.xml'
                   .format(session._name, method))
    # output histogram
    hname = 'MVA_{}'.format(method)
    hists_discr[method] = TH1D(hname, hname, 100, -1.0, 1.0 )

    # discriminant for output tree
    discriminant[method] = array('f', [0.])
    otree.Branch(method, discriminant[method], '{}/F'.format(method))

combinedVarsDict = {\
        "Bs_vtx_chi2_ndof":      "lab0_VCHI2NDOF",\
        "Bs_radial_fd":          "lab0_RFD",\
        "Ds_radial_fd":          "lab2_RFD",\
        "Ds_vtx_chi2_ndof":      "lab2_VCHI2NDOF",\
        "max_ghost_prob":        "lab1345_TRACK_GhostProb",\
        "min_Ds_child_trk_pt":   "lab345_MIN_PT",\
        "min_Ds_child_trk_chi2": "lab3_IPCHI2_OWNPV"}
varVals = [0] * len(allvars)
for i, var in enumerate(session.all_vars()):
    if ":=" in var:
        varVals[i] = combinedVarsDict[var.split(":=")[0]]
    else:
        varVals[i] = var

nentries = itree.GetEntries()
for i in xrange(nentries):
    itree.GetEntry(i)
    if i % 10000 == 0: print "{}/{} ({}%)".format(i, nentries, (100. * i) / nentries)
    for j, var in enumerate(session.all_vars()):
        allvars[j][0] = getattr(itree, varVals[j])
    for method in session.methods:
        discriminant[method][0] = reader.EvaluateMVA(method)
        hists_discr[method].Fill(discriminant[method][0])
    otree.Fill()

for method in session.methods:
    hists_discr[method].Write()
otree.Write()
ofile.Close()
ifile.Close()
