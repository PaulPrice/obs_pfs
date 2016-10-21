#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python PfsFiberTrace.py
or
   python
   >>> import FiberTrace; FiberTrace.run()
"""
import os
import unittest
import lsst.daf.persistence as dafPersist
import lsst.utils
import lsst.utils.tests as tests
import pfs.drp.stella as drpStella
import lsst.obs.pfs.utils as obsPfsUtils
#import subprocess

class PfsFiberTraceTestCase(tests.TestCase):
    """A test case for comparing a reconstructed FiberTraceSet to the original"""

    def setUp(self):
        drpStellaDataDir = lsst.utils.getPackageDir("drp_stella_data")
        self.butlerDir = os.path.join(drpStellaDataDir,"tests/data/PFS")
        self.butler = dafPersist.Butler(self.butlerDir)
        self.dataId = dict(visit=5, ccd=5, spectrograph=2, arm='r', dateObs='2016-01-12')
        self.flat = self.butler.get("postISRCCD", self.dataId, immediate=True)
        print 'self.flat = ',self.flat
        
        self.ftffc = drpStella.FiberTraceFunctionFindingControl()
        self.ftffc.fiberTraceFunctionControl.order = 5
        self.ftffc.fiberTraceFunctionControl.xLow = -5
        self.ftffc.fiberTraceFunctionControl.xHigh = 5
        self.ftffc.fiberTraceFunctionControl.interpolation = "POLYNOMIAL"
        self.ftffc.apertureFWHM = 2.6
        self.ftffc.signalThreshold = 10
        self.ftffc.nTermsGaussFit = 3
        self.ftffc.saturationLevel = 65550.
        self.ftffc.minLength = 3880
        self.ftffc.maxLength = 3930
        self.ftffc.nLost = 20
        
        self.ftpfc = drpStella.FiberTraceProfileFittingControl()
        
        # This particular flatfile has 11 FiberTraces scattered over the whole CCD
        # If in the future the test data change we need to change these numbers
        self.nFiberTraces = 11
        
    def tearDown(self):
        del self.flat
        del self.ftffc

    def testPfsFiberTrace(self):
        """Test that we can create a pfsFiberTrace"""
        
        fiberTraceSet = drpStella.findAndTraceAperturesF(self.flat.getMaskedImage(), self.ftffc)
        """Check that we found self.nFiberTraces FiberTraces"""
        self.assertEqual(fiberTraceSet.size(), self.nFiberTraces)

        """Check that we can set the FiberTraceProfileFittingControl ftpfc"""
        self.assertTrue(fiberTraceSet.setFiberTraceProfileFittingControl(self.ftpfc))
        
        """Check that we can calculate the spatial profile for all the FiberTraces in fiberTraceSet"""
        self.assertTrue(fiberTraceSet.calcProfileAllTraces())
        
        """Check that we can write a PfsFiberTrace from the fiberTraceSet"""
        dataId = dict(calibDate=self.dataId['dateObs'], spectrograph=self.dataId['spectrograph'], arm=self.dataId['arm'])
        pfsFiberTrace = obsPfsUtils.createPfsFiberTrace(dataId, fiberTraceSet, self.flat.getHeight())
        dirName = '~/tmp'
        fileName = 'FiberTrace.fits'
        pfsFiberTrace.write(dirName, fileName)
        
        """Check that we can read a PfsFiberTrace back in"""
        pfsFiberTraceNew = PfsFiberTrace()
        print 'pfsFiberTraceNew = ',pfsFiberTraceNew
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(PfsFiberTraceTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--display", '-d', default=False, action="store_true", help="Activate display?")
    parser.add_argument("--verbose", '-v', type=int, default=0, help="Verbosity level")
    args = parser.parse_args()
    display = args.display
    verbose = args.verbose
    run(True)
