# coding=utf-8
"""
This module provides helper functions to get different ROOT and RooFit
objects.

@author Suvayu Ali

"""

# Python modules
import os
import sys
from datetime import datetime

# # FIXME: Batch running fails on importing anything but gROOT
# ROOT global variables
from ROOT import gROOT, gSystem

# ROOT classes
from ROOT import TTree, TFile, TCanvas, TPad, TClass

# RooFit classes
from ROOT import RooFit
from ROOT import RooPlot, RooWorkspace
from ROOT import RooArgSet, RooArgList
from ROOT import RooAbsArg, RooAbsReal, RooAbsPdf, RooAbsData
from ROOT import RooRealVar
from ROOT import RooDataSet, RooDataHist

# Hack around RooWorkspace.import() and python keyword import clash
_import = getattr(RooWorkspace, 'import')


# Python wrappers
def get_timestamp(fmt='%Y-%m-%d-%a-%H-%M'):
    """Return formatted timestamp."""
    return datetime.strftime(datetime.today(), fmt)


# ROOT wrappers
def load_library(library):
    loadstatus = { 0: 'loaded',
                   1: 'already loaded',
                   -1: 'does not exist',
                   -2: 'version mismatch' }
    status = gSystem.Load(library)
    if status < 0:
        sys.exit('Problem loading %s, %s' % (library, loadstatus[status]))


def get_file(fname, mode='read'):
    """Open and return a ROOT file."""
    if os.path.exists(fname):
        return TFile(fname, mode)
    else:
        raise IOError('File %s does not exist!' % fname)


def get_object(objname, rfile, subdir=''):
    """Get a ROOT object from a ROOT file by name."""
    if len(subdir): subdir += '/'
    return rfile.Get('%s%s' % (subdir, objname))


# RooFit wrappers
def get_argset(args):
    """Return and argset of the RooFit objects."""
    argset = RooArgSet()
    for arg in args:
        if arg.InheritsFrom(RooAbsArg.Class()): argset.add(arg)
        else: TypeError('%s should inherit from RooAbsArg' % arg.GetName())
    return argset


def get_arglist(args):
    """Return and arglist of the RooFit objects."""
    arglist = RooArgList()
    for arg in args:
        if arg.InheritsFrom(RooAbsArg.Class()): arglist.add(arg)
        else: TypeError('%s should inherit from RooAbsArg' % arg.GetName())
    return arglist


def set_integrator_config():
    """Configure numerical integration in RooFit.

    More precise integrals in RooFit and how intervals are
    determined and integrals calculated.

    """

    # More precise integrals in RooFit
    RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-9)
    RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-9)
    # Set how intervals are determined and integrals calculated
    RooAbsReal.defaultIntegratorConfig().getConfigSection('RooAdaptiveGaussKronrodIntegrator1D').setCatLabel('method','21Points')
    RooAbsReal.defaultIntegratorConfig().getConfigSection('RooAdaptiveGaussKronrodIntegrator1D').setRealValue('maxSeg', 1000)
    RooAbsReal.defaultIntegratorConfig().method1D().setLabel('RooAdaptiveGaussKronrodIntegrator1D')
    RooAbsReal.defaultIntegratorConfig().method1DOpen().setLabel('RooAdaptiveGaussKronrodIntegrator1D')


def get_toy_dataset(varargset, PDF=None):
    """Return a toy dataset for the given PDF."""

    objclass = TClass.GetClass(PDF.ClassName())
    if objclass.InheritsFrom(RooAbsPdf.Class()):
        dataset = PDF.generate(varargset, 10000, RooFit.Name('toydataset'),
                               RooFit.Verbose(True))
        print 'Toy generation completed with %s' % PDF.GetName()
        return dataset
    else:
        raise TypeError('PDF should inherit from RooAbsPdf.')


def get_dataset(varargset, ftree, cut='', cutVars=[], *cmdArgs):
    """Return a dataset.

    Return a dataset from the ntuple `ftree'. Apply a selection cut
    using the `cutVar' variable and the selection `cut'.

    """

    varargsetclone = varargset.clone('varargsetclone')
    for cvar in cutVars:
        varargsetclone.add(cvar) # Add selVar to apply cut

    tmpdataset = RooDataSet('dataset', 'Dataset', varargsetclone,
                            RooFit.Import(ftree), RooFit.Cut(cut), *cmdArgs)
    dataset = tmpdataset.reduce(varargset)
    del tmpdataset
    return dataset


def save_in_workspace(rfile, **argsets):
    """Save RooFit objects in workspace and persistify.

    Pass the different types of objects as a keyword arguments. e.g.
    save_in_workspace(pdf=[pdf1, pdf2], variable=[var1, var2])

    """

    # Persistify variables, PDFs and datasets
    workspace = RooWorkspace('workspace',
                             'Workspace saved at %s' % get_timestamp())
    keys = argsets.keys()
    for key in keys:
        print 'Importing RooFit objects in %s list.' % key
        for arg in argsets[key]: _import(workspace, arg)
    rfile.WriteTObject(workspace)
    print 'Saving arguments to file: %s' % rfile.GetName()


def get_workspace(fname, wname):
    """Read and return RooWorkspace from file."""
    ffile = get_file(fname, 'read')
    workspace = ffile.Get(wname)
    return workspace, ffile
