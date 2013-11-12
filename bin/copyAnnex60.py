#!/usr/bin/python
#####################################################################
# This script copies the Annex 60 library to the Buildings library.
#
# Todo: copy unit tests and test results.
#       - copy Images directory
#
# MWetter@lbl.gov                                          2011-05-15
#####################################################################
def replaceAnnex60String(filNam):
    """ Replace the string `Annex60` with `Buildings`.
    
    :param filNam: Name of the file
    """
    import os
    
    f = open(filNam, 'r')
    lines = list()
    for lin in f.readlines():
        if filNam.endswith(".mo"):
            lin = string.replace(lin, "Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature",
                              "Buildings.HeatTransfer.Sources.PrescribedTemperature")
            lin = string.replace(lin, "Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow",
                              "Buildings.HeatTransfer.Sources.PrescribedHeatFlow")
        lines.append(string.replace(lin, "Annex60", "Buildings"))   
    f.close
    if filNam.startswith(BLDREFRES):
        # The reference results have the library in their file name.
        f = open(filNam.replace("Annex60", "Buildings"), 'w')
        os.remove(filNam)
    else:
        f = open(filNam, 'w')
    f.writelines(lines)
    f.close()


import os, string, fnmatch, os.path, sys
from os import makedirs
from os.path import isfile, join
import shutil

# --------------------------
# Global settings
ANNEX60HOME=os.path.abspath("/Users/mwetter/proj/ldrd/bie/modeling/github/iea-annex60/modelica-annex60/Annex60")
BLDHOME=os.path.abspath("/Users/mwetter/proj/ldrd/bie/modeling/github/lbl-srg/modelica-buildings/Buildings")
ANNEX60REFRES=os.path.join(ANNEX60HOME, "Resources", "ReferenceResults")
BLDREFRES=os.path.join(BLDHOME, "Resources", "ReferenceResults")


for root, dirs, files in os.walk(ANNEX60HOME):
    for fil in files:
        srcFil=os.path.join(root, fil)
        # Loop over all .mo files except for top-level .mo file
        if (srcFil.endswith(".mo") and (not srcFil.endswith(os.path.join("Annex60", "package.mo")))) \
            or (srcFil.endswith(".mos") or srcFil.startswith(ANNEX60REFRES)):
            desFil=srcFil.replace(ANNEX60HOME, BLDHOME)
            desPat=os.path.dirname(desFil)
            if not os.path.exists(desPat):
                os.makedirs(desPat)
            shutil.copy2(srcFil, desFil)
            replaceAnnex60String(desFil)
        
#    os.path.exists(path)
#    filNam = helpDir + os.path.sep + fil
#    filObj=open(filNam, 'r')
#    lines = filObj.readlines()
#    filObj.close()
#    for old, new in replacements.iteritems():
#        for i in range(len(lines)):
#            lines[i] = lines[i].replace(old, new)
#    filObj=open(filNam, 'w')
#    for lin in lines:
        # Check if line contains a wrong string
#        validateLine(lin)
#        filObj.write(lin)
#    filObj.close()
