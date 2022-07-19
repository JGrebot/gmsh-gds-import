# !/usr/bin/python3
import sys
import gds

####################################################################################################
#
# Test1 : DESCRIPTION
#
# Before launching this script : please view gds content with gds.view function (see corresponding test)
#
# Purpose of this script is to test gds.Read function and to show all information extracted from 
# example.gds.
#
# This time we use the new function, depending on the gdspy package. 

cell_name="example2"
gds_filename="data/example.gds"

gds.Read(gds_filename, cell_name)
gds.DATA[-1].display()



