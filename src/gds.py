# !/usr/bin/python3

import sys
import os
import re
import numpy as np
import gmsh
import gdspy

# # To check which module was imported :
# os.path.abspath(gmsh.__file__)

model   = gmsh.model
factory = model.occ


GDS_already_read = []
DATA = []

def retrieveObjectInsideBox(x, y, z, dx, dy, dz, tol = 1.0e-3):
   """
   Function retrieving the 3D entities inside the box (even if the
   entities boundaries are intersecting with the box's boundaries)

   x          : x coordinate of the box's origin
   y          : y coordinate of the box's origin
   z          : z coordinate of the box's origin
   dx         : Length of the box
   dy         : Width of the box
   dz         : Height of the box
   tol        : Tolerance used to retrieve the entities with bounding boxes
   """
   list = np.array([])
   top = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                 y - tol,\
                                                 z - tol,\
                                                 x + dx + tol,\
                                                 y + dy + tol,\
                                                 z + dz - tol,\
                                                 dim = 3))
   if(top.size != 0):
      list = np.concatenate((list, top[:,1]), axis = 0)
   south = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                   y + tol,\
                                                   z - tol,\
                                                   x + dx + tol,\
                                                   y + dy + tol,\
                                                   z + dz + tol,\
                                                   dim = 3))
   if(south.size != 0):
      list = np.concatenate((list, south[:,1]), axis = 0)
   east = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                  y - tol,\
                                                  z - tol,\
                                                  x + dx - tol,\
                                                  y + dy + tol,\
                                                  z + dz + tol,\
                                                  dim = 3))
   if(east.size != 0):
      list = np.concatenate((list, east[:,1]), axis = 0)
   north = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                   y - tol,\
                                                   z - tol,\
                                                   x + dx + tol,\
                                                   y + dy - tol,\
                                                   z + dz + tol,\
                                                   dim = 3))
   if(north.size != 0):
      list = np.concatenate((list, north[:,1]), axis = 0)
   west = np.array(model.getEntitiesInBoundingBox(x + tol,\
                                                  y - tol,\
                                                  z - tol,\
                                                  x + dx + tol,\
                                                  y + dy + tol,\
                                                  z + dz + tol,\
                                                  dim = 3))
   if(west.size != 0):
      list = np.concatenate((list, west[:,1]), axis = 0)
   bot = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                 y - tol,\
                                                 z + tol,\
                                                 x + dx + tol,\
                                                 y + dy + tol,\
                                                 z + dz + tol,\
                                                 dim = 3))
   if(bot.size != 0):
      list = np.concatenate((list, bot[:,1]), axis = 0)
   return(list)


def retrieveVolumeBoxTag(x, y, z, dx, dy, dz, tol = 1.0e-3):
   """
   Function retrieving the tag of a 3D entity shaped like a box whose
   dimensions are (dx, dy, dz). If they are any 3D entities contained
   in the box, the function is only retrieving the tag of the box and
   not the tags of smaller entities inside.

   x          : x coordinate of the box's origin
   y          : y coordinate of the box's origin
   z          : z coordinate of the box's origin
   dx         : Length of the box 
   dy         : Width of the box
   dz         : Height of the box
   tol        : Tolerance used to retrieve the entities with bounding boxes
   """

   #retrieve all the entities in the volume shaped by addBox(x - eps, ..., x + dx + eps, ...)
   volume_set = np.array(model.getEntitiesInBoundingBox(x - tol,\
                                                        y - tol,\
                                                        z - tol,\
                                                        x + dx + tol,\
                                                        y + dy + tol,\
                                                        z + dz + tol,\
                                                        dim = 3))[:,1]

   if(volume_set.size == 0):
      print("No 3D entity have been found in the selected domain")
      return(-1)

   #retrieve all the entities in the volume shaped by addBox(x + eps, ..., x + dx - eps, ...)
   volume_bar = retrieveObjectInsideBox(x, y, z, dx, dy, dz)
   if(volume_bar.size == 0):
      return(volume_set)
   else:
      #list difference between volume_set and volume_bar
      return(np.setdiff1d(volume_set, volume_bar))




      
def Build(gds_filename, cell_name, index, datatype, z0, dz=0.0, Sf=0.0, option='extrude',option_param=[],debug=False):
    """
    ## Description
    Function building all object inside a layer of gds_filename (.gds). All surfaces
    inside this layer are built. Intersection surface are automaticlly fuse. "Simple
    polygone" (with or without hole) are accepted. Surfaces are supposed to by given 
    on X-Y plane.

    ## Args:
        gds_filename (str):             Name of gds file.
        cell_name (str):                Name of cell of interest of the gds file.
        index (int):                    Index of gds layer to build (number of left in klayout).
        datatype (int):                 Datatype of gds layer to build (number of right in klayout).
        z0 (float):                     z origin (surfaces are suppose to be given in X-Y plane).
        dz (float, optional):           Thickness of z dimension (supported : dz > 0 and dz < 0). Defaults to 0.0.
        Sf (float, optional):           Global Scaling factor. Defaults to 0.0.
        option (str, optional):         Type of volume construction (for now only 'extrude'). Defaults to 'extrude'.
        option_param (list, optional):  free extra parameters for volume construction (will be necessary for pyramid). Defaults to [].#                                       
        debug (bool, optional):         For debug purpose. Defaults to False.

    Returns:
        dimTag:                         dimTag of all the built objects.
    """
    global DATA
    global GDS_already_read
    
    # Checking
    if Sf==0.0:
        print("Error in gds.Build function.")
        print("you must provide Sf argrument (your Scaling factor)")
        sys.exit()
    Possible = ['extrude', 'pyramid']
    if option not in Possible:
        print("Error in gds.Build_layer.")
        print("option provided is not correct.")
        print("Possible options are for now : \'"+Possible[0]+"\' or \'"+Possible[1]+"\' ")
        sys.exit()
    if dz == 0.0:
        print("Error in gds.Build_layer.")
        print("dz is equals to 0.0")
        print("But we need dz > 0.0 or dz < 0.0!")
        sys.exit()
    
    # Read_NEW will automatically read filename only if it has not already been read.
    if debug:
        Read(gds_filename, cell_name, debug=debug)
    elif not debug:
        Read(gds_filename, cell_name)
        
    
    # Retrive index of gds_filename in DATA list of gds STRUCT (see end of Read function)
    i = GDS_already_read.index(gds_filename)
    
    scaling = DATA[i].file_unit_scaling/Sf
    
    # retrive indice of l_GDS_layer where index and datatype are provided by user.
    n = -1
    for j in range(0,len(DATA[i].l_GDS_layer)):
        if DATA[i].l_GDS_layer[j].index == index and DATA[i].l_GDS_layer[j].datatype == datatype:
            n = j
    if n == -1:
        print("Error in gds.Build.")
        print("Pair (index,datatype) is not found in gds provided.")
        print("To remember : index of a layer inside klayout editor is number of left.")
        print("              and datatype of a layer inside klayout editor is number of right.")
        print(" ")
        print(" ")
        print("If you are certain index and datatype are correct, then you have a problem with")
        print("gds.Read function. ")
        sys.exit()
    
        
    tag_surface = []        # Contains tags for simple surface
    list_S_hole = []        # Contains tags for surface with hole
    
    # Adding all surface defined in gds_layer n in gmsh model.
    # Circle on all element for selected layer.
    for k in range(0, len(DATA[i].l_GDS_layer[n].list_of_element)):
        
        # Retriving number of dot and coordinate of this element in Nb_dot and Pc (Point Coordinate)
        Nb_dot = int(len(DATA[i].l_GDS_layer[n].list_of_element[k].coor)/2)
        Pc = np.zeros([Nb_dot,2])
        for j in range(0,Nb_dot):
            Pc[j,0] = DATA[i].l_GDS_layer[n].list_of_element[k].coor[2*j]*scaling               # x coordinate
            Pc[j,1] = DATA[i].l_GDS_layer[n].list_of_element[k].coor[2*j+1]*scaling             # y coordinate
        
        Nb_row_occurence = (Pc[:, np.newaxis] == Pc).all(axis=2).sum(axis=1)
        Nb_of_2 = np.count_nonzero(Nb_row_occurence == 2) 
             
        # The number of identical dots is : Nb_of_2                         
        
        tag_first_dot = 0
        tag_last_dot = 0
        tag_previous_dot = 0
        
        Loop = []
        Loop_outter = []
        Loop_inner = []
        tmp_tag_d1 = -1
        tmp_tag_d2 = -1
        
        # -------------------------------------------------------------------------
        # If there is only TWO dots equal
        # Then this element is just a closed, usal, nice and friendly loop.
        if Nb_of_2==2:
            
            Loop = np.zeros(Nb_dot-1)
            
            for lo in range(0,Nb_dot):
                if lo < Nb_dot-2:
                    if lo == 0:
                        # Add in gmsh the first 2 dots
                        tmp_tag_d1 = factory.addPoint(Pc[lo+0,0],Pc[lo+0,1],z0)
                        tmp_tag_d2 = factory.addPoint(Pc[lo+1,0],Pc[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_first_dot = tmp_tag_d1
                        tag_previous_dot = tmp_tag_d2
                        Loop[lo]=tmp_tag_line
                    else:
                        # Add in gmsh
                        # tmp_tag_d1 = factory.addPoint(Pc[l0+0,0],Pc[l0+0,1],z0)         # This dot has already been added.
                        tmp_tag_d1 = tag_previous_dot
                        tmp_tag_d2 = factory.addPoint(Pc[lo+1,0],Pc[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_previous_dot = tmp_tag_d2
                        Loop[lo]=tmp_tag_line
                    
                if lo == Nb_dot-2:
                    tag_last_dot = tag_previous_dot
                    # Add last line that link last dot and first dot
                    tmp_tag_line = factory.addLine(tag_last_dot, tag_first_dot) 
                    Loop[lo]=tmp_tag_line
                    
                    
                if lo == Nb_dot-1:
                    # Retrieve coordinate
                    P1 = Pc[lo,:]
                    P2 = Pc[0,:]
                    
                    # No need to add a line here because these 2 dots must/will/should be equal
                    
                    # Checking that last dot is equals to the first one
                    if P1[0]!=P2[0] or P1[1]!=P2[1]:
                        print("Error in gds.Build :")
                        print("It appears that ")
                        print("          on gds         ",i)
                        print("          on gds_layer   ",n,", with index ",DATA[i].l_GDS_layer[n].index," and datatype ",DATA[i].l_GDS_layer[n].datatype)
                        print("          on surface ",k)
                        print("Your loop is not defining close subset of R2 plane.")
                        sys.exit()
                
            tmp = factory.addCurveLoop(Loop)
            tag_surface.append(factory.addPlaneSurface([tmp]))
            
            
        
        # ---------------------------------------------------------------------------
        # If there is FOUR dots equal
        # Then this element is a shape with a hole.
        # Dots are like A--B----B-A
        # We don't add B. And on the center is defined the loop of inside shape.
        # And with right and left is defined the loop on larger shape, including A.
        elif Nb_of_2==4:
            if Pc[0,0]!=Pc[-1,0] or Pc[0,1]!=Pc[-1,1]:
                print("Error in gds.Build")
                print("The first and the last dot should be the same here")
                print("This error comes from case : Nb_of_2==4")
                sys.exit()
            index_redondant_dots = np.where(Nb_row_occurence==2)[0]
            if Nb_of_2!=len(index_redondant_dots):
                print("Error in gds.Build")
                print("This error comes from case : Nb_of_2==4")
                sys.exit()
                
            index_redondant_dots = np.where(Nb_row_occurence==2)[0]
            B_index1 = index_redondant_dots[1]
            B_index2 = index_redondant_dots[2]
            
            # Number of dot of outter and inner loop
            Nb_dot_outter = B_index1 + Nb_dot - B_index2 - 2
            Nb_dot_inner = B_index2 - B_index1 -1
            
            # Retriving coordinate of outter and inner dots
            Pc_outter = np.zeros([Nb_dot_outter,2])
            Pc_inner = np.zeros([Nb_dot_inner,2])
            
            for lo in range(0,Nb_dot_outter):
                if lo < B_index1:
                    Pc_outter[lo,:] = Pc[lo,:]
                elif lo >= B_index1:
                    Pc_outter[lo,:] = Pc[B_index2 - B_index1 + 1 + lo,:]
                
            for lo in range(0,Nb_dot_inner):
                Pc_inner[lo,:] = Pc[lo+B_index1+1,:]
            
            if debug:
                print('index_redondant_dots = ',index_redondant_dots)
                print('Nb_dot = ',Nb_dot)
                print('B_index1 =',B_index1)
                print('B_index2 =',B_index2)
                print('Pc =',Pc)
                print('Pc_outter =',Pc_outter)
                print('Nb_dot_outter =',Nb_dot_outter)
                print('Pc_inner =',Pc_inner)
                print('Nb_dot_inner =',Nb_dot_inner)

            Loop_outter = np.zeros(Nb_dot_outter)
            Loop_inner = np.zeros(Nb_dot_inner)
            
            # Now we have two usual and simple loops to build.
            # Using 3 for circle isn't optimal. a while might have been better here. But
            # when it works it works..
            # First circle for outter shape
            for lo in range(0,B_index1):
                if lo == 0:
                    if B_index1 != 1:
                        # Add dot to gmsh
                        tmp_tag_d1 = factory.addPoint(Pc_outter[lo+0,0],Pc_outter[lo+0,1],z0)
                        tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
            
                        tag_first_outter_dot = tmp_tag_d1
                        tag_previous_outter_dot = tmp_tag_d2
                        Loop_outter[lo] = tmp_tag_line
                    else:
                        tmp_tag_d1 = factory.addPoint(Pc_outter[lo+0,0],Pc_outter[lo+0,1],z0)
                        tag_first_outter_dot = tmp_tag_d1
                        tag_previous_outter_dot = tmp_tag_d1
                        
                elif lo < B_index1-1:
                    # Add dot to gmsh
                    tmp_tag_d1 = tag_previous_outter_dot
                    tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0) 
                    tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
        
                    tag_previous_outter_dot = tmp_tag_d2
                    Loop_outter[lo] = tmp_tag_line
            # Second circle for outter shape        
            for lo in range(B_index1, Nb_dot_outter):
                if lo == B_index1:
                    tmp_tag_d0 = tag_previous_outter_dot
                    tmp_tag_d1 = factory.addPoint(Pc_outter[lo,0],Pc_outter[lo,1],z0)
                    tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0)
                    tmp_tag_line0 = factory.addLine(tmp_tag_d0, tmp_tag_d1) 
                    tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
        
                    tag_previous_outter_dot = tmp_tag_d2
                    Loop_outter[lo-1] = tmp_tag_line0
                    Loop_outter[lo] = tmp_tag_line
                    
                if lo < Nb_dot_outter-1 and lo > B_index1:
                    tmp_tag_d1 = tag_previous_outter_dot
                    tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0)
                    tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
        
                    tag_previous_outter_dot = tmp_tag_d2
                    Loop_outter[lo] = tmp_tag_line
                        
                if lo == Nb_dot_outter-1:
                    tmp_tag_d1 = tag_previous_outter_dot
                    tmp_tag_d2 = tag_first_outter_dot
                    tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                    
                    Loop_outter[lo] = tmp_tag_line
            
            # Third circle for inner shape
            for lo in range(0, Nb_dot_inner):
                if lo < Nb_dot_inner-1:
                    if lo == 0:
                        # Add in gmsh the first 2 dots
                        tmp_tag_d1 = factory.addPoint(Pc_inner[lo+0,0],Pc_inner[lo+0,1],z0)
                        tmp_tag_d2 = factory.addPoint(Pc_inner[lo+1,0],Pc_inner[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_first_inner_dot = tmp_tag_d1
                        tag_previous_dot = tmp_tag_d2
                        Loop_inner[lo]=tmp_tag_line
                    else:
                        # Add in gmsh
                        # tmp_tag_d1 = factory.addPoint(Pc[l0+0,0],Pc[l0+0,1],z0)         # This dot has already been added.
                        tmp_tag_d1 = tag_previous_dot
                        tmp_tag_d2 = factory.addPoint(Pc_inner[lo+1,0],Pc_inner[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_previous_dot = tmp_tag_d2
                        Loop_inner[lo]=tmp_tag_line
                    
                if lo == Nb_dot_inner-1:
                    tag_last_dot = tag_previous_dot
                    # Add last line that link last dot and first dot
                    tmp_tag_line = factory.addLine(tag_last_dot, tag_first_inner_dot) 
                    Loop_inner[lo]=tmp_tag_line
                    
            # print("Loop_outter =", Loop_outter)
            # print("Loop_inner =", Loop_inner)
            
            tmp = factory.addCurveLoop(Loop_outter)
            tag_surface_outter = factory.addPlaneSurface([tmp])
            tmp = factory.addCurveLoop(Loop_inner)
            tag_surface_inner = factory.addPlaneSurface([tmp])
            list_S_hole.append(S_hole(tag_surface_outter, tag_surface_inner))
        
         
        # ---------------------------------------------------------------------------
        # If there is SIX dots equal
        # Then this element is a shape with a hole.
        # Dots are like A--BC----CB-A
        # We don't add B and C. And on the center is defined the loop of inside shape.
        # And with right and left is defined the loop on larger shape, including A.
        elif Nb_of_2==6:
            if Pc[0,0]!=Pc[-1,0] or Pc[0,1]!=Pc[-1,1]:
                print("Error in gds.Build")
                print("The first and the last dot should be the same here")
                print("This error comes from case : Nb_of_2==6")
                sys.exit()
            index_redondant_dots = np.where(Nb_row_occurence==2)[0]
            if Nb_of_2!=len(index_redondant_dots):
                print("Error in gds.Build")
                print("This error comes from case : Nb_of_2==6")
                sys.exit()
                
            B_index1 = index_redondant_dots[1]
            C_index2 = index_redondant_dots[4]
            
            Nb_dot_outter = B_index1 + Nb_dot - C_index2 - 2
            Nb_dot_inner = C_index2 - B_index1 -3
            Pc_outter = np.zeros([Nb_dot_outter,2])
            Pc_inner = np.zeros([Nb_dot_inner,2])
            
            Loop_outter = np.zeros(Nb_dot_outter)
            Loop_inner = np.zeros(Nb_dot_inner)
            
            for lo in range(0,Nb_dot_outter):
                if lo < B_index1:
                    Pc_outter[lo,:] = Pc[lo,:]
                elif lo >= B_index1:
                    Pc_outter[lo,:] = Pc[C_index2 - B_index1 + 1 + lo,:]
                
            for lo in range(0,Nb_dot_inner):
                Pc_inner[lo,:] = Pc[lo+B_index1+2,:]
            
            if debug:
                print(index_redondant_dots)
                print('Nb_dot = ',Nb_dot)
                print('B_index1 =',B_index1)
                print('C_index2 =',C_index2)
                print('Pc =',Pc)
                print('Pc_outter =',Pc_outter)
                print('Nb_dot_outter =',Nb_dot_outter)
                print('Pc_inner =',Pc_inner)
                print('Nb_dot_inner =',Nb_dot_inner)
            
            # Now we have two usual and simple loops to build.
            # Loop for outter shape
            for lo in range(0, Nb_dot_outter):
                if lo < Nb_dot_outter-1:
                    if lo == 0:
                        # Add in gmsh the first 2 dots
                        tmp_tag_d1 = factory.addPoint(Pc_outter[lo+0,0],Pc_outter[lo+0,1],z0)
                        tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_first_outter_dot = tmp_tag_d1
                        tag_previous_dot = tmp_tag_d2
                        Loop_outter[lo]=tmp_tag_line
                    else:
                        # Add in gmsh
                        # tmp_tag_d1 = factory.addPoint(Pc[l0+0,0],Pc[l0+0,1],z0)         # This dot has already been added.
                        tmp_tag_d1 = tag_previous_dot
                        tmp_tag_d2 = factory.addPoint(Pc_outter[lo+1,0],Pc_outter[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_previous_dot = tmp_tag_d2
                        Loop_outter[lo]=tmp_tag_line
                    
                if lo == Nb_dot_outter-1:
                    tag_last_dot = tag_previous_dot
                    # Add last line that link last dot and first dot
                    tmp_tag_line = factory.addLine(tag_last_dot, tag_first_outter_dot) 
                    Loop_outter[lo]=tmp_tag_line
            
            # Loop for inner shape
            for lo in range(0, Nb_dot_inner):
                if lo < Nb_dot_inner-1:
                    if lo == 0:
                        # Add in gmsh the first 2 dots
                        tmp_tag_d1 = factory.addPoint(Pc_inner[lo+0,0],Pc_inner[lo+0,1],z0)
                        tmp_tag_d2 = factory.addPoint(Pc_inner[lo+1,0],Pc_inner[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_first_inner_dot = tmp_tag_d1
                        tag_previous_dot = tmp_tag_d2
                        Loop_inner[lo]=tmp_tag_line
                    else:
                        # Add in gmsh
                        # tmp_tag_d1 = factory.addPoint(Pc[l0+0,0],Pc[l0+0,1],z0)         # This dot has already been added.
                        tmp_tag_d1 = tag_previous_dot
                        tmp_tag_d2 = factory.addPoint(Pc_inner[lo+1,0],Pc_inner[lo+1,1],z0) 
                        tmp_tag_line = factory.addLine(tmp_tag_d1, tmp_tag_d2) 
                        
                        tag_previous_dot = tmp_tag_d2
                        Loop_inner[lo]=tmp_tag_line
                    
                if lo == Nb_dot_inner-1:
                    tag_last_dot = tag_previous_dot
                    # Add last line that link last dot and first dot
                    tmp_tag_line = factory.addLine(tag_last_dot, tag_first_inner_dot) 
                    Loop_inner[lo]=tmp_tag_line
                    
            # print("Loop_outter =", Loop_outter)
            # print("Loop_inner =", Loop_inner)
            
            tmp = factory.addCurveLoop(Loop_outter)
            tag_surface_outter = factory.addPlaneSurface([tmp])
            tmp = factory.addCurveLoop(Loop_inner)
            tag_surface_inner = factory.addPlaneSurface([tmp])
            list_S_hole.append(S_hole(tag_surface_outter, tag_surface_inner))
            
        else:
            print("Error in gds.Build function :")
            print("There is, inside ",gds_filename," with index = ",index," and with datatype = ",datatype)
            print("element number ",k," has a number of equal dots that is not equal to 2, 4 or 6")
            print("This case is not supposed to happen for now.")
            sys.exit()
                    
        if debug:
            print("Nb_of_2 = ",Nb_of_2)
            print("tag_surface = ",tag_surface)
            print("tag_surface_outter = ",tag_surface_outter)
            print("tag_surface_inner = ",tag_surface_inner)
    #   END OF for ON k
    #   Now we have built all surface and we stored their tags in "tag_surface" and in "list_S_hole".
    
    
    
    
    tag_volume_surface = np.zeros(len(tag_surface))
    tag_volume_S_hole = np.zeros(len(list_S_hole))
    
    #################################################
    # Building corresponding volumes for all surfaces
    for k in range(0, len(tag_surface)+len(list_S_hole)):
        
        # To know which type of element we need to build
        if k<len(tag_surface):
            type = 'surface'
        else:
            type = 'S_hole'
            k = k - len(tag_surface)
            
        if option=='extrude':
            if type == 'surface':
                dimtags  = factory.extrude([(2,tag_surface[k])], 0.0, 0.0, dz)
                for lo in dimtags:
                    if lo[0]==3:
                        tag_volume_surface[k] = lo[1]
            elif type == 'S_hole':
                # We want to extrude inner and outter and cut outter ny inner
                dimtags_inner  = factory.extrude([(2,list_S_hole[k].tag_inner)], 0.0, 0.0, dz)
                dimtags_outter  = factory.extrude([(2,list_S_hole[k].tag_outter)], 0.0, 0.0, dz)
                for lo in dimtags_inner:
                    if lo[0]==3:
                        tag_inner = lo[1]
                for lo in dimtags_outter:
                    if lo[0]==3:
                        tag_outter = lo[1]
                factory.synchronize()
                dimtags = factory.cut([(3,tag_outter)], [(3,tag_inner)], removeObject=True, removeTool=True)
                for lo in dimtags:
                    if lo[0][0]==3:
                        tag_volume_S_hole[k] = lo[0][1]
        elif option=='pyramid':
            print('Problem with build gds.Build :')
            print("option pyramid is not implemented yet")
            print("En fait on va rencontrer, dans le cas de tranchée sous la forme d'une tente inversée, le problème de l'orientation de celless-ci. ")
            print("Pour résoudre cela :")
            print("    1) on attend que Guillaume code la fonction correspondante : FAIT")
            print("    2) On va utiliser l'argument option_orientation pour spécifier l'orientation des tranchées.")
            print("       et cet argument sera dans un parametre.py")
            print("")
            print(" Mais cela c'est pour la suite. ( le 09/09/2020 )")
            sys.exit()
    
    
    if debug:
        dimtag = model.getEntities(dim=3)
        print('dimtag of volumes :',dimtag)
        print('tag_volume_surface = ',tag_volume_surface)
        print('tag_volume_S_hole = ',tag_volume_S_hole)
    
    all_tag_volume = np.concatenate((tag_volume_surface, tag_volume_S_hole), axis=0)
    
    
    
    ################################################
    # Fusing all volumes added.
    stopping = False
    go_on = True
    guard = 0
    
    
    # While will stop only when no intersection can be found.
    while go_on:
        # If we have more than one volume
        if len(all_tag_volume)>1:
            # Fuse volumes that have an intersection (their tags are in all_tag_volume)
            for k in range(0, len(all_tag_volume)):
                for l in range(0, len(all_tag_volume)):
                    if k!=l: # and all_tag_volume[k]!=-1 and all_tag_volume[l]!=-1:
                         # creates an object if there is an intesection
                        intersect = factory.intersect([(3,all_tag_volume[k])],[(3,all_tag_volume[l])], removeObject=False, removeTool=False)[0]
                        
                        # if there is an intersection,
                        if len(intersect):
                            # If intersection is volume whose tag is k, then we only retrive volume k.
                            if intersect[0][1] == all_tag_volume[k]:
                                factory.remove(intersect, True)  # remove volume k.
                                all_tag_volume[k] = -1
                            # If intersection is volume whose tag is l, then we only retrive volume l.
                            elif intersect[0][1] == all_tag_volume[l]:
                                factory.remove(intersect, True)  # remove volume l.
                                all_tag_volume[l] = -1
                            else:
                                # else we remove intersection, volume k and volume l and add the fuse between k and l.
                                dimtags = factory.fuse([(3,all_tag_volume[k])],[(3,all_tag_volume[l])],removeObject=True, removeTool=True)[0]
                                factory.remove(intersect, True)  # remove created intersection objects
                                for lo in dimtags:
                                    if lo[0]==3:
                                        tag = lo[1]
                                        
                                all_tag_volume[k] = -1
                                all_tag_volume[l] = -1
                                all_tag_volume = np.concatenate((all_tag_volume, [tag]), axis=0)
                                
                            # Refreshing all_tag_volumes
                            all_tag_volume = all_tag_volume[all_tag_volume!=-1]
                            # Putting stopping to True will break the two "for" loop
                            stopping = True
                            # "go_on" to True so we continue while loop.
                            go_on = True
                        else:
                            # if there is zero intersection, then go_on will stop the while loop,
                            go_on = False  
                            
                    # If we fused an element, then we break both "for" loop.
                    if stopping:
                        break
                if stopping:
                    stopping = False
                    break
            guard +=1
            # Always be careful with a while loop.
            if guard > 10000:
                print('Error on gds.Build function.')
                print('     When fusing volumes of a gds layer, program is stuck in an infinite loop.')
                sys.exit()
        else:
            # In this case there is only one volume left, so no fusing are necessary.
            go_on = False
    
    if debug:
        print("all_tag_volume = ",all_tag_volume)
        sys.exit()
    
    dimTag = []
    for k in range(0,len(all_tag_volume)):
        dimTag.append( (3, int(all_tag_volume[k])) )
    
    return dimTag



def Read(gds_filename, cell_name, debug=False):
    """
    ## Description
    This function reads all informations in .gds file provided and store it in an object of type
    GDS. If a .gds has already been read this function does nothing. Storing support reading
    more than one .gds in a .py script.

    ## Args:
        gds_filename (str):     name of .gds file.
        cell_name (str):        cell of interest in the .gds file.
        debug (bool, optional): For debug purpose. Defaults to False.
    """
    global GDS_already_read
    global DATA
        
    gds_filename_without_folder = gds_filename.split('/')[-1]    
    current_dir = os.getcwd()+'/'
    
    # Checking gds_filename is a file and it exists
    if not os.path.isfile(current_dir + gds_filename):
        print("Error in Read_GDS function : ")
        print(".gds file is not found.")
        print("Check that this file is next to your initial .py script.")
        sys.exit()
        
        
    # Checking it is a gds file.
    if ".gds" not in gds_filename:
        print("Error in Read_GDS function : ")
        print("It looks like file provided is not a .gds file.")
        sys.exit()
        
    
    # A guard to read this .gds only one time.
    if gds_filename not in GDS_already_read:
        
        a_gds = GDS()
        
        Read_gds_file_with_gdspy(gds_filename,  cell_name, a_gds)
        
        # a guard to read this .gds only one time.
        GDS_already_read.append(gds_filename)
        # Store in DATA
        DATA.append(a_gds) 
        
        if debug:
            DATA[-1].display()



def Read_gds_file_with_gdspy(gds_filename, cell_name, a_gds):
    """
    ## Description
    Intermediate function that read the gds with the gdpy package.

    ## Args:
        gds_filename (str):     name of .gds file.
        cell_name (str):        cell of interest in the .gds file.
        debug (bool, optional): For debug purpose. Defaults to False.
    """
    
    # Load the GDSII file into a new library
    gdsii = gdspy.GdsLibrary(infile=gds_filename)
    
    a_gds.file_unit_scaling = gdsii.unit
    
    All_layers_datatype = []
    tmp = []
    # Retrieve all couple (layer, datatype)
    for pol in gdsii.cells[cell_name].polygons:
        tmp = [pol.layers[0], pol.datatypes[0]]
        if tmp not in All_layers_datatype:
            All_layers_datatype.append(tmp)
        
        
    # For each couple (layer, datatype)
    for ld in All_layers_datatype:
        tmp_le = []
        
        # Retrieve all polygons that have this layer and this datatype
        for pol in gdsii.cells[cell_name].polygons:
            if pol.layers[0]==ld[0] and pol.datatypes[0]==ld[1]:
                coor = []
                for co in pol.polygons[0]:
                    coor.append(co[0])
                    coor.append(co[1])
                # We need to duplicate the first element and add it at the end 
                # in order to match previous format.
                coor.append(coor[0])
                coor.append(coor[1])
                tmp_le.append(GDS_Element(coor, ld[0], ld[1]))
        
        a_gds.l_GDS_layer.append(GDS_Layer(tmp_le, ld[0], ld[1]))
        
    

def View(gds_filename, cell_name):
    """
    ## Description:
    This function reads all informations in .gds file provided and display it in a nice GUI windows. 
    It requires the python3 gdspy package. And the package python3-tk (sudo apt install python3-tk)
    This new function for reading use the gdspy library instead of the GDSIICONVERT command.


    ## Args:
        gds_filename (str): name of .gds file.
        cell_name (str):    name of cell you want to visualize.
    """
    # Load the GDSII file into a new library
    gdsii = gdspy.GdsLibrary(infile=gds_filename)
    gdspy.LayoutViewer(gdsii, gdsii.cells[cell_name], background='#FFFFFF')
    sys.exit()

    

##############################################################################################
#
#       STORAGE CLASS DEFINITION
#
# These classes are meant to store all datas that are contained in a .gds file.
#
#
# GDS: 
#       l_GDS_layer :        (list) a list of GDS_Layer.
#       file_unit_scaling :  (float) proper unit of gds data.
#       display :            a display method for debug purpose
#
# GDS_Layer: 
#       list_of_element :    (list) a list of GDS_element.
#       index :              (int) exact index of gds_layer (number of left in klayout editor)
#       datatype :           (int) exact datatype of gds_layer (number of right in klayout editor)
#       GDS_layer_display:   a display method for debug purpose
#
# GDS_Element: 
#       coor :                          (list) a list of GDS_Layer.
#       parent_GDS_layer_index :        (int) index of parent GDS_Layer.
#       parent_GDS_layer_datatype :     (int) datatype of parent GDS_Layer.
#       GDS_element_display :           a display method for debug purpose
#
#
# S_hole:  a surface with a hole. so two surface, one outter, one inner.
#       tag_outter :        (int) tag of outter surface.
#       tag_inner :         (int) tag of inner surface.
################################################################################################


class GDS: 
    def __init__(self, l_GDS_layer = [],file_unit_scaling=0):  
        self.l_GDS_layer = l_GDS_layer.copy() 
        self.file_unit_scaling = file_unit_scaling
    def display(self):
        print(' ')
        print(' ')
        print('Info on object of type GDS')
        print('     it contains : ', len(self.l_GDS_layer) , 'GDS_Layer.')
        print('     its file unit is : ', self.file_unit_scaling , ' meter')
        print(' ')
        for i in range(0,len(self.l_GDS_layer)):
            self.l_GDS_layer[i].GDS_layer_display()
        
class GDS_Layer:
    def __init__(self, list_of_element=[], index=-1, datatype=-1):  
        self.list_of_element = list_of_element.copy()
        self.index = index 
        self.datatype = datatype  

    def GDS_layer_display(self):
        print(' ')
        print(' ')
        print(' ')
        print(' ')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(' ')
        print('Info on  a  GDS_layer :    index = ',self.index,'    datatype = ',self.datatype )
        print(' ')
        print('Its elements are the following :')
        for i in range(0,len(self.list_of_element)):
            self.list_of_element[i].GDS_element_display()
    


class GDS_Element:
    def __init__(self, coor, parent_GDS_layer_index, parent_GDS_layer_datatype): 
        self.coor = np.copy(coor)
        self.parent_GDS_layer_index = parent_GDS_layer_index
        self.parent_GDS_layer_datatype = parent_GDS_layer_datatype

    def GDS_element_display(self):
        print('One element with ',len(self.coor),'/2 dots : XY = ', self.coor)

class S_hole:
    def __init__(self, tag_outter, tag_inner):
        self.tag_outter = tag_outter
        self.tag_inner = tag_inner


#/////////////////////////////////////////////////////////////////////////////////////////////

  
##############################################################################################
#
#       PARSING
#
# A usual parser using regular expression.
# This is not used in the Read function and a remainder of the usage of a previous library.
#
##############################################################################################


# set up regular expressions
# use https://regex101.com/ to visualise these if required  
rx_dict = {
    'Scaling': re.compile(r'^\* Unit=(?P<gds_scaling>.*?) meters \(file units = {(?P<a_scaling>.*?),(?P<file_unit_scaling>.*?)}\)'),
    'Element': re.compile('^\s*Element (?P<tag>\d+): BOUNDARY \(layer (?P<index>\d+), datatype (?P<datatype>\d+)\)'),
    'XY': re.compile('\s*XY: (?P<coor>.+) \n'),
}


# Line parser
def _parse_line(line):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None



# File parser
def parse_file(filepath,a_gds):
    """
    Parse text at given filepath
    Fill information in a GDS structure

    Parameters
    ----------
    filepath : str
        Filepath for file_geo to be parsed
    a_gds :    an object of class GDS
    Returns
    -------
    a_gds :    an object of class GDS

    """
    # Local Variable
    
    # iter_element = 0
    current_layer_index = -1
    current_layer_datatype = -1   #  current gds_layer index and current gds_layer datatype
    list_of_element = []
    counter = 0
    # current_element = -1
    
    # Initiate a GDS_Layer that will be filled
    
    # open the file and read through it line by line
    with open(filepath, 'r') as file_geo :
        
        line = file_geo.readline()
        while line:
            # at each line check for a match with a regex
            key, match = _parse_line(line)

            if key == 'Scaling':
                a_gds.file_unit_scaling = float(match.group('file_unit_scaling'))
            
            if key == 'Element':
                # If we have a new index and new datatype
                if current_layer_index!=int(match.group('index')) or current_layer_datatype!=int(match.group('datatype')):
                    # If this is the first time
                    if counter==0:
                        current_layer_index= int(match.group('index'))
                        current_layer_datatype= int(match.group('datatype'))
                        counter = counter + 1
                    else:
                        a_gds.l_GDS_layer.append(GDS_Layer(list_of_element, current_layer_index, current_layer_datatype))
                        list_of_element.clear()
                        current_layer_index= int(match.group('index'))
                        current_layer_datatype= int(match.group('datatype'))
                
                
            if key == 'XY':
                tmp_coor = [int(i) for i in match.group('coor').split(' ')]
                list_of_element.append(GDS_Element( tmp_coor, parent_GDS_layer_index=current_layer_index, parent_GDS_layer_datatype=current_layer_datatype ))
                
                        
            line = file_geo.readline()
        # END of while
            
        # Adding last gds_layer
        a_gds.l_GDS_layer.append(GDS_Layer(list_of_element, current_layer_index, current_layer_datatype))
        list_of_element.clear()

    return a_gds
