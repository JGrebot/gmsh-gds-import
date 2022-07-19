# !/usr/bin/python3
import sys
import gds
####################################################################################################
#
# Test1 : DESCRIPTION
#
#
# Purpose of this script is to test gds.View function that  show all informations inside
# data/example.gds.
#
# This function requires the gdspy python3 package.
# and it requires tkinter : sudo apt install python3-tk

cell_name = "example2"
gds_filename = "data/example.gds"

gds.View(gds_filename, cell_name)
