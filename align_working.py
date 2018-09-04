from pymol import cmd
import numpy as np
#import numpy

def get_pos(selection):
	pos_sel = []
	for atom in cmd.get_model(selection).atom:
		pos_sel.append([atom.coord[0], atom.coord[1], atom.coord[2]])
	print(pos_sel)
	return pos_sel

cmd.reinitialize()
cmd.load("session_iodine_minus1.pse", partial=1)
cmd.do("run axes.py")
#cmd.load("geometry_complex_I.mol")
cmd.select("circ_complex","organic")

cmd.fetch("4agl")
cmd.select("ligand","organic")
cmd.select("chain_a","chain a")
cmd.select("ligand_a", "ligand and chain_a")

cmd.select("halligand", "last ligand_a and (name I* or name Br* or name Cl*)")
cmd.select("halcomplex", "first circ_complex and (name I* or name Br* or name Cl*)")

cmd.select("ciligand", "ligand_a and name C* and bound_to halligand")
cmd.select ("cicomplex", "name C* and bound_to halcomplex") 

translation_vector = 0
pos_halcomplex =[]
pos_halligand =[]
pos_cicomplex = []
pos_ciligand = []

cmd.iterate_state(1, 'halligand', 'pos_halligand.append((x,y,z))')
cmd.iterate_state(1, 'halcomplex', 'pos_halcomplex.append((x,y,z))')


cmd.iterate_state(1, 'cicomplex', 'pos_cicomplex.append((x,y,z))')
cmd.iterate_state(1, 'ciligand', 'pos_ciligand.append((x,y,z))')

pos_cicomplex = np.array(pos_cicomplex)
pos_ciligand = np.array(pos_ciligand)
pos_halcomplex = np.array(pos_halcomplex)

pos_c2complex = []
pos_c2ligand = []
pos_c3complex = []
pos_c3ligand = []



cmd.select ("c2complex", "first bound_to cicomplex and name c*")
cmd.select ("c3complex", "last bound_to cicomplex and name c*") 


cmd.select ("c2ligand", "first bound_to ciligand and name c*" )
cmd.select ("c3ligand", "last bound_to ciligand and name c*" )

cmd.iterate_state(1, 'c2ligand', 'pos_c2ligand.append((x,y,z))')

print("halogen of complex: "+str(pos_halcomplex))
print("halogen of ligand: " + str(pos_halligand)) 
print(abs(pos_cicomplex))
translation_vector = -pos_halcomplex
print("Translation vector: ")
print(translation_vector)
print([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]])

#translate molecule from mol file to 0,0,0
cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],"circ_complex",-1,0)

#get coordinates from ligand and mol file complex
pos_halligand = []
pos_halcomplex = []
cmd.iterate_state(1, 'halligand', 'pos_halligand.append((x,y,z))')
cmd.iterate_state(1, 'halcomplex', 'pos_halcomplex.append((x,y,z))')
pos_halligand = np.array(pos_halligand)

#add pseudoatom on x axis to measure angle for rotation around y axis
cmd.pseudoatom("x-1",pos=[-1,0,0])
#measure angle for molfile rotation
angle = cmd.get_angle(atom1="cicomplex",atom2="halcomplex",atom3="x-1")
#rotate mol file
cmd.rotate("z",180,"circ_complex",0,0,origin=[0,0,0])

#cmd.pseudoatom("e",pos=[-2,0,z_temp])


print(angle)
print(pos_halligand[0][2])
cmd.pseudoatom("z1",pos=[0,0,1])


cmd.pseudoatom("pseudo_ciligand",pos=[pos_ciligand[0][0],pos_ciligand[0][1],pos_ciligand[0][2]])
cmd.pseudoatom("pseudo_c2ligand",pos=[pos_c2ligand[0][0],pos_c2ligand[0][1],pos_c2ligand[0][2]])
cmd.pseudoatom("pseudo_halligand",pos=[pos_halligand[0][0],pos_halligand[0][1],pos_halligand[0][2]])
cmd.select("pseudo_ligand","pseudo_*ligand")


#----START CALCULATIONS WITH PSEUDO LIGAND TO GET ANGLES FOR LATER ROTATION OF CIRCULAR COMPLEX---------------


translation_vector = -pos_halligand
cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],"pseudo_ligand",-1,0)

cmd.create("pseudo_ligand_object","pseudo_ligand")


pos_halligand =[]
cmd.iterate_state(1, 'pseudo_halligand', 'pos_halligand.append((x,y,z))')

angle1 = cmd.get_angle(atom1="z1",atom2="pseudo_halligand",atom3="pseudo_ciligand")
print("Rotate ligand_a around angle " + str(angle))
#FIRST ROTATION OF PSEUDOLIGAND
cmd.rotate("x",angle1,"pseudo_ligand", 0 , 0, origin=pos_halligand[0])

pos_ciligand =[]
pos_c2ligand = []
cmd.iterate_state(1, 'pseudo_ciligand', 'pos_ciligand.append((x,y,z))')
cmd.iterate_state(1, 'pseudo_c2ligand', 'pos_c2ligand.append((x,y,z))')

z_temp = pos_c2ligand[0][2]

cmd.pseudoatom("c",pos=[0,0,z_temp])
cmd.pseudoatom("d",pos=[-2,0,z_temp])
angle2 = cmd.get_angle(atom1="d",atom2="c",atom3="pseudo_c2ligand")
print(angle)


#SECOND ROTATION OF LIGAND
cmd.rotate("z",-angle2,"pseudo_ligand", 0, 0, origin=[0,0,0])
print("Rotate ligand_a around angle " + str(angle))
cmd.iterate_state(1, 'pseudo_halligand', "print x,y,z")

angle3 = cmd.get_angle(atom1="pseudo_ciligand",atom2="pseudo_halligand",atom3="x-1")
print(angle)
#THIRD ROTATION OF LIGAND
cmd.rotate("y",-angle3,"pseudo_ligand",0,0,origin=pos_halligand[0])
print("Rotate ligand_a around angle " + str(angle))
cmd.iterate_state(1, 'pseudo_halligand', "print x,y,z")

cmd.hide("sticks")
#----------------END CALCULATIONS WITH PSEUDO LIGAND----------------------------------------------------------------

#--------------START ROTATION OF CIRC COMPLEX-----------------
pos_halligand =[]
cmd.iterate_state(1, 'halligand', 'pos_halligand.append((x,y,z))')

translation_vector = pos_halligand
print angle1,angle2,angle3
cmd.rotate("y",angle3,"circ_complex",0,0,origin=pos_halcomplex[0])
cmd.rotate("z",angle2,"circ_complex",0,0,origin=pos_halcomplex[0])
cmd.rotate("x",-angle1,"circ_complex",0,0,origin=pos_halcomplex[0])

cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],"circ_complex",-1,0)

pos_halligand = []
pos_ciligand = []
cmd.iterate_state(1, "halligand", 'pos_halligand.append((x,y,z))')
cmd.iterate_state(1, "ciligand", 'pos_ciligand((x,y,z))')
pos_ciligand = np.array(pos_ciligand)
pos_halligand= np.array(pos_halligand)

print("TEST")
print(get_pos("ciligand"))


cmd.zoom("circ_complex")




"""
#After translation find two vectors v1 and v2 from ci to each iodin
pos_halligand = numpy.array(pos_halligand)
pos_halcomplex = numpy.array(pos_halcomplex)

v1 = pos_halligand - pos_ciligand
v2 = pos_halcomplex - pos_ciligand

#now find orthogonal vector to v1 and v2. later rotate around this vector
v3= numpy.cross(v1,v2)
print(v3[0])

angle = cmd.get_angle(atom1="halligand",atom2="ciligand",atom3="halcomplex")
cmd.rotate([v3[0][0],v3[0][1],v3[0][2]], -angle, "circ_complex",origin=pos_ciligand[0])

cmd.iterate_state(1, 'c2ligand', 'pos_c2ligand.append((x,y,z))')
cmd.iterate_state(1, 'c3ligand', 'pos_c3ligand.append((x,y,z))')
cmd.iterate_state(1, 'c2complex', 'pos_c2complex.append((x,y,z))')
cmd.iterate_state(1, 'c3complex', 'pos_c3complex.append((x,y,z))')

pos_halligand = []
cmd.iterate_state(1, 'halligand', 'pos_halligand.append((x,y,z))')
pos_halligand = numpy.array(pos_halligand)

pos_c2complex = numpy.array(pos_c2complex)
pos_c2ligand = numpy.array(pos_c2ligand)

v1 = pos_ciligand - pos_halligand
v2 = pos_c2ligand - pos_ciligand
v3= numpy.cross(v1,v2)
print("ciligand")
cmd.iterate_state(1, 'ciligand', "print x,y,z")

angle = cmd.get_angle(atom1="c2ligand",atom2="ciligand",atom3="c2complex")

cmd.rotate([v1[0][0],v1[0][1],v1[0][2]], angle, "circ_complex",origin=pos_ciligand[0])

angle = cmd.get_angle(atom1="c2ligand",atom2="ciligand",atom3="c2complex")
cmd.rotate([v1[0][0],v1[0][1],v1[0][2]], angle, "circ_complex",origin=pos_ciligand[0])



model_circ_complex = cmd.get_model("circ_complex")

angle = cmd.get_angle(atom1="c2ligand",atom2="ciligand",atom3="c2complex")
cmd.rotate([v1[0][0],v1[0][1],v1[0][2]], -angle, "circ_complex",origin=pos_ciligand[0])

angle = cmd.get_angle(atom1="c2ligand",atom2="ciligand",atom3="c2complex")
cmd.rotate([v1[0][0],v1[0][1],v1[0][2]], angle, "circ_complex",origin=pos_ciligand[0])




"""
