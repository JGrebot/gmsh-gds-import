# !/usr/bin/python3
import sys
import gmsh
import gds

model   = gmsh.model
factory = model.occ

def writeMeshFile(filename, view = False, generate = True, dim = 3):

   if(generate):
      model.mesh.generate(dim)
   if(view):
      gmsh.fltk.run()
   else:
      gmsh.write(filename)
   gmsh.finalize()



####################################################################################################
#
# Test2 : DESCRIPTION
#
# Before launching this script : please view gds content with gds.view function (see corresponding test)
#
# Purpose of this script is to show and test :
#    1) extrusion of all gds surfaces with same z_min and dz
#    2) extrusion of a box
#    3) extrusion of simple polygon without hole (triangle / octogone / crown)
#    4) extrusion of simple polyton with a hole. (top right shape (visualise example.gds with klayout))
#    5) (if some extrusions intersect each other) fusion of concerned extrusion. 
#
# This time we use the new function, depending on the gdspy package. 


ScalingFactor = 1.0e-9

# Initializing gmsh to save ALL elements
gmsh.initialize(sys.argv)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 100)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 150)
# gmsh.option.setNumber("Mesh.SaveElementTagType", 2)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.ScalingFactor", ScalingFactor)
gmsh.option.setNumber("General.Terminal", 1)


# Add all object inside layer (2,0) of gds example.gds with 'extrude' option
cell_name="example2"
gds_filename="data/example.gds"
gds.Build(gds_filename, cell_name, index=2, datatype=0, z0=0, dz=500, Sf=ScalingFactor, option='extrude')

factory.synchronize()
writeMeshFile("test2_Build_extrude.mesh")






