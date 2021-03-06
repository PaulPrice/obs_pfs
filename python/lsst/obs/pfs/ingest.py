#ingestImages.py '/Users/azuri/spectra/pfs/PFS/' --mode=link '/Users/azuri/spectra/pfs/raw/2016-01-12/*.fits' --loglevel 'info' -C /Users/azuri/stella-git/obs_pfs/config/pfs/ingest.py
import re
import os
from lsst.pipe.tasks.ingest import ParseTask, ParseConfig
import lsst.afw.image as afwImage

class PfsParseConfig(ParseConfig):
    def setDefaults(self):
        ParseConfig.setDefaults(self)
        self.translators["field"] = "translate_field"

class PfsParseTask(ParseTask):
    ConfigClass = PfsParseConfig

    def getInfo(self, filename):
        """Get information about the image from the filename and its contents

        @param filename    Name of file to inspect
        @return File properties; list of file properties for each extension
        """
        #print 'interpreting filename <'+filename+'>'
        minSpectrograph = 1
        maxSpectrograph = 4
        minArmNum = 1
        maxArmNum = 4
        maxCCDs = 4
        arms = ['b', 'r', 'n', 'm']
        sites = '[JLXIASPF]'
        categories = '[ABCDS]'
        matches = re.search("PF("+sites+")("+categories+")(\d{6})(\d)(\d).fits", filename)
        if not matches:
            """ Try old format """
            matches = re.search("PF("+sites+")("+categories+")-(\d{6})(\d)(\d).fits", filename)
            if not matches:
                raise RuntimeError("Unable to interpret filename: %s" % filename)
        site, category, visit, spectrograph, armNum = matches.groups()
        if int(spectrograph) == 9:
            spectrograph = '2'
        if int(spectrograph) < minSpectrograph or int(spectrograph) > maxSpectrograph:
            message = 'spectrograph (=',spectrograph,') out of bounds'
            raise Exception(message)
        ccd = int(spectrograph)-1
        if int(armNum) == minArmNum:
            """do nothing"""
        elif (int(armNum) == minArmNum + 1) or (int(armNum) == minArmNum + 3):
            ccd += maxCCDs
        elif int(armNum) == minArmNum + 2:
            ccd += 2*maxCCDs
        else:
            message = 'arm number (=',armNum,') out of bounds [',minArmNum,'...',maxArmNum,']'
            raise Exception(message)

        filter = arms[ int(armNum) - minArmNum ]
        arm = filter
        print 'site = <'+site+'>, category = <'+category+'>, visit = <'+visit+'>, spectrograph = <'+spectrograph+'>, armNum = <'+armNum+'>, filter = <'+filter+'>, ccd = ',ccd

        info = dict(site=site, category=category, visit=int(visit, base=10), filter=filter, arm=arm, spectrograph=int(spectrograph), ccd=ccd)
        if os.path.exists(filename):
            header = afwImage.readMetadata(filename)
            info = self.getInfoFromMetadata(header, info=info)
        return info, [info]

    def translate_field(self, md):
        """Get 'field' from IMAGETYP

        This is temporary, until an OBJECT (or similar) header can be provided,
        but it's better than setting everything to the same thing.
        """
        field = md.get("IMAGETYP").strip()
        if field in ("#", ""):
            field = "UNKNOWN"
        self.log.info('PfsParseTask.translate_field: field = %s' % field)
        return re.sub(r'\W', '_', field).upper()
