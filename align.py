from pymol import cmd
import numpy as np
import pymol
from os.path import isfile, join
from os import listdir


files = [f for f in listdir("/pdbs") if isfile(join(".", f)) and ".pdb" in f and not ".pse" in f]

#INPUT PDB NAME HERE
pdb = "1IE4"
#INPUT WHICH HALOGEN BOND TO USE FOR ALIGNMENT IF THERE ARE MULTIPLE
hal_bond = 0
#INPUT SESSION NAME
se ="session_chlorine_minus1_merged"



"""
#Function that gets the position of a pymol selection and returns it as a python list
"""

def get_pos(selection):
	pos_sel = []
	for atom in cmd.get_model(selection).atom:
		pos_sel.append([atom.coord[0], atom.coord[1], atom.coord[2]])
	return pos_sel

"""
Function that rotates a selection around the given angle until the reference atom (zero) lies on a plain.
"""


# all backbone oxygen atoms have the name O --> select blub, n. O
# all water oxygens have the name O --> select all waters with cmd.select("water", "resn hoh") and use the selection to exclude ...

def rotate_sel(selection,axis,angle,zero):
	position = get_pos(zero)
	if axis == "x":
		if position[0][1] >= 0.0 and position[0][2] >= 0:
			cmd.rotate("x",-(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-(180-angle)) 
				+ " degree.")
			return -(180-angle)
		if position[0][1] < 0.0 and position[0][2] >= 0:
			cmd.rotate("x",180-angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return 180-angle
		if position[0][1] < 0 and position[0][2] < 0:
			cmd.rotate("x",angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return angle
		if position[0][1] >= 0 and position[0][2] < 0:
			cmd.rotate("x",-angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-angle) 
				+ " degree.")
			return -angle

	if axis == "y":
		if position[0][0] >= 0.0 and position[0][2] >= 0:
			cmd.rotate("y",-(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-(180-angle)) 
				+ " degree.")
			return -(180-angle)
		if position[0][0] < 0.0 and position[0][2] >= 0:
			cmd.rotate("y",-angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return -angle
		if position[0][0] < 0 and position[0][2] < 0:
			cmd.rotate("y",angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return angle
		if position[0][0] >= 0 and position[0][2] < 0:
			cmd.rotate("y",(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-angle) 
				+ " degree.")
			return 180-angle

	if axis == "z":
		if position[0][0] >= 0.0 and position[0][1] >= 0:
			cmd.rotate("z",-(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-(180-angle)) 
				+ " degree.")
			return -(180-angle)
		if position[0][0] < 0.0 and position[0][1] >= 0:
			cmd.rotate("z",-angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return -angle
		if position[0][0] < 0 and position[0][1] < 0:
			cmd.rotate("z",angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return angle
		if position[0][0] >= 0 and position[0][1] < 0:
			cmd.rotate("z",(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-angle) 
				+ " degree.")
			return 180-angle


"""
This function serves as a wrapper function for the rotate_sel function.
It takes a selection to rotate and a reference atom as input. Then rotates around the 
z axis until the x coordinate of ref_atom is 0.
"""

def zero_x(selection,ref_atom):
	try:
		cmd.remove("pseud_ang1")
		cmd.remove("pseud_ang2")
	except pymol.CmdException:
		print("No atoms removed")
	cmd.pseudoatom("pseud_ang1",pos=[0,0,get_pos(ref_atom)[0][2]])
	cmd.pseudoatom("pseud_ang2",pos=[0,get_pos(ref_atom)[0][1],get_pos(ref_atom)[0][2]])
	angle = cmd.get_angle(atom1=ref_atom,atom2="pseud_ang1",atom3="pseud_ang2")
	ang = rotate_sel(selection,"z",angle,ref_atom)
	return ang

"""
This function serves as a wrapper function for the rotate_sel function.
It takes a selection to rotate and a reference atom as input. Then rotates around the 
x axis until the y coordinate of ref_atom is 0.
"""

def zero_y(selection,ref_atom):
	try:
		cmd.remove("pseud_ang1")
		cmd.remove("pseud_ang2")
	except pymol.CmdException:
		print("No atoms removed")
	cmd.pseudoatom("pseud_ang1",pos=[get_pos(ref_atom)[0][0],0,0])
	cmd.pseudoatom("pseud_ang2",pos=[get_pos(ref_atom)[0][0],0,get_pos(ref_atom)[0][2]])
	angle = cmd.get_angle(atom1=ref_atom,atom2="pseud_ang1",atom3="pseud_ang2")
	ang = rotate_sel(selection,"x",angle,ref_atom)
	return ang

"""
This function serves as a wrapper function for the rotate_sel function.
It takes a selection to rotate and a reference atom as input. Then rotates around the 
y axis until the z coordinate of ref_atom is 0.
"""

def zero_z(selection,ref_atom):
	try:
		cmd.remove("pseud_ang1")
		cmd.remove("pseud_ang2")
	except pymol.CmdException:
		print("No atoms removed")
	cmd.pseudoatom("pseud_ang1",pos=[0,get_pos(ref_atom)[0][1],0])
	cmd.pseudoatom("pseud_ang2",pos=[get_pos(ref_atom)[0][0],get_pos(ref_atom)[0][1],0])	
	angle = cmd.get_angle(atom1=ref_atom,atom2="pseud_ang1",atom3="pseud_ang2")
	ang = rotate_sel(selection,"y",angle,ref_atom)
	return ang

"""
Function that aligns the selection to the negative x-axis of the origin.
"""

def align_x_neg(selection,ref_atom):
	angles = []
	angle1 = zero_x(selection,ref_atom)
	angle2 = zero_y(selection,ref_atom)
	angle3 = zero_z(selection,ref_atom)
	angles.append(angle1)
	angles.append(angle2)
	angles.append(angle3)
	return angles

"""
Function to get position of polar atoms that form halogen bonds with ligand
"""
def find_halogen_bonds():
	halogen_bonds = []
	cmd.select("bindingsite1", "ligand_a expand 5.0")	# expand selection by 5 angstrom, extend = extend along bonds
	cmd.select("bindingsite1", "br. bindingsite")	# br extends the selection to the full residue
	cmd.select("bindingsite1", "bindingsite and not ligand")	# br extends the selection to the full residue
	cmd.select("water", "resn hoh")	# resn = residues name, hoh = name of water molecules
	cmd.select("water_bs", "water and bindingsite")

	#select all relevant halogen atoms
	number_of_halogens = cmd.select("hal", "hal and ligand_a")
	cmd.select("polar_atoms", "e. O or e. N or e. S")

	model_halogen = cmd.get_model("hal")	# create model for halogen selection
	for halo in model_halogen.atom:
		print("A")
		cmd.select("current_halogen", "(e. I or e. Cl or e. Br) and id %s"%(halo.id))
		cmd.select("current_neighbor", "current_halogen extend 1")
		cmd.select("current_neighbor", "current_neighbor and not current_halogen")
		cmd.select("current_neighbor", "neighbor current_halogen")

		cmd.select("current_surroundings", "current_halogen expand 5.0")
		cmd.select("current_surroundings", "current_surroundings and bindingsite1")
		cmd.select("current_polar", "current_surroundings and polar_atoms")

		model_polar = cmd.get_model("current_polar")
		for pol in model_polar.atom:
			cmd.select("current_polar", "not water and polar_atoms and id %s"%(pol.id))
			# check for halogen bonds, distances 5.0, angles: > 140
			current_distance = cmd.distance("current_distance", "current_halogen", "current_polar")
			current_angle = cmd.angle("current_angle", "current_neighbor", "current_halogen", "current_polar")
			print(current_distance,current_angle)
			if (current_angle >= 140) and (current_distance <= 5.0):
				print current_distance, current_angle, pol.name, pol.resn, pol.resi, pol.id
				halogen_bonds.append((pol.resn,pol.resi,pol.id))
				

			cmd.delete("current_angle")
			cmd.delete("current_distance")
			
			

	cmd.delete("current_halogen")
	cmd.delete("current_neighbor")
	cmd.delete("water")
	cmd.delete("water_bs")
	cmd.delete("current_polar")
	cmd.delete("current_surroundings")
	cmd.delete("bindingsite1")
	return halogen_bonds

def align_ligand(selection,center,angles):
	print(center)
	pos_halligand = get_pos(center)
	print(pos_halligand)
	translation_vector = pos_halligand

	cmd.rotate("x",-angles[3],selection,0,0,origin=[0,0,0])
	cmd.rotate("y",-angles[2],selection,0,0,origin=[0,0,0])
	cmd.rotate("x",-angles[1],selection,0,0,origin=[0,0,0])
	cmd.rotate("z",-angles[0],selection,0,0,origin=[0,0,0])

	cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],selection,-1,0)
	
	
def init():
	cmd.reinitialize()
	#LOAD FILES AND MAKE THE NECESSARY SELECTIONS FOR ALIGNMENT
	cmd.load("session_chlorine_minus1_merged.pse")
	cmd.do("run axes.py")
	cmd.do("run center_of_mass.py")
	#cmd.load("geometry_complex_Cl.mol")
	cmd.select("circ_complex","organic")
	
	cmd.fetch(pdb)
	cmd.select("water", "resn hoh")
	#select ligand and bindingsite
	cmd.select("ligand","organic and not circ_complex")
	cmd.select("chain_a","chain a")
	cmd.select("ligand_a", "ligand and chain_a")

for file in files:
	pdb = file[0:4]
	print "Current File: " + str(pdb)
	halligands = []



	init()

	cmd.select("halligands", "ligand and (name I* or name Br* or name Cl*)")

	for halligand in cmd.get_model("halligands").atom:
		halligands.append(halligand.id)

	c_hal=0
	for halligand in halligands:
		init()
		cmd.select("halligand", "id " +str(halligand)+" and ligand and (name I* or name Br* or name Cl*)")
		cmd.select("hal", "(e. Cl or e. Br or e. I) and (halligand expand 5)") 
		#find the important bonds for the alignment
		hal_bonds = find_halogen_bonds()
		if len(hal_bonds) > 0:
			cmd.select("bindingsite", "halligand expand 5.0")	# expand selection by 5 angstrom, extend = extend along bonds
			cmd.select("halcomplex", "first circ_complex and (name I* or name Br* or name Cl*)")

			cmd.select("circ_complex_only","halcomplex extend 5")

			cmd.select("comp_back","circ_complex and not circ_complex_only")
			cmd.select("comp_back_o","comp_back and name o")
			cmd.select("comp_back_c","bound_to comp_back_o")
			cmd.select("comp_back_n","(bound_to comp_back_c) and name N")

			print(hal_bonds)

			#select atoms of ligand backbone for alignment
			cmd.select("back_polar","resi "+str(hal_bonds[0][1])+" and id "+ str(hal_bonds[0][2]))
			model_back_polar = cmd.get_model("back_polar").atom
			if "S" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","first bound_to back_o")
				cmd.select("back_n","last bound_to back_o")
			if "O" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","bound_to back_o")
				cmd.select("back_n","first bound_to back_c")

			if "N" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","first bound_to back_o")
				cmd.select("back_n","first bound_to back_c")

			cmd.select("ligand_circ","br. halligand")
			cmd.select("ligand_back","br. back_o")

			#select atoms of halogen structure of complex
			cmd.select("circ_complex","circ_complex_only")
			cmd.select("ciligand", "ligand and name C* and bound_to halligand")
			cmd.select ("cicomplex", "name C* and bound_to halcomplex") 

			translation_vector = 0
			pos_halcomplex = get_pos("halcomplex")
			pos_halligand =get_pos("halligand")
			pos_cicomplex = get_pos("cicomplex")
			pos_ciligand = get_pos("ciligand")


			#get positions of ligand and complex atoms previously selected
			pos_cicomplex = np.array(pos_cicomplex)
			pos_ciligand = np.array(pos_ciligand)
			pos_halcomplex = np.array(pos_halcomplex)

			cmd.select ("c2complex", "first bound_to cicomplex and name c*")
			cmd.select ("c2ligand", "first bound_to ciligand and name c*" )
			translation_vector1 = -pos_halcomplex

			#translate molecule from mol file to 0,0,0
			cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"circ_complex",-1,0)

			#get coordinates from ligand and mol file complex
			pos_halligand = get_pos("halligand")
			pos_halcomplex = get_pos("halcomplex")
			pos_halligand = np.array(pos_halligand)

			#add pseudoatom on x axis to measure angle for rotation around y axis

			cmd.rotate("z",180,"circ_complex",0,0,origin=[0,0,0])
			cmd.rotate("x",180,"circ_complex",0,0,origin=[0,0,0])

			cmd.rotate("x",180,"comp_back",0,0,origin=[0,0,0])
			pos_c2ligand = get_pos("c2ligand")
			cmd.pseudoatom("pseudo_ciligand",pos=[pos_ciligand[0][0],pos_ciligand[0][1],pos_ciligand[0][2]])
			cmd.pseudoatom("pseudo_c2ligand",pos=[pos_c2ligand[0][0],pos_c2ligand[0][1],pos_c2ligand[0][2]])
			cmd.pseudoatom("pseudo_halligand",pos=[pos_halligand[0][0],pos_halligand[0][1],pos_halligand[0][2]])
			cmd.select("pseudo_ligand","pseudo_*ligand")

			#ADD PSEUDO BACKBONE

			cmd.pseudoatom("pseudo_back_o",pos=[get_pos("back_o")[0][0],get_pos("back_o")[0][1],get_pos("back_o")[0][2]])
			cmd.pseudoatom("pseudo_back_c",pos=[get_pos("back_c")[0][0],get_pos("back_c")[0][1],get_pos("back_c")[0][2]])
			cmd.pseudoatom("pseudo_back_n",pos=[get_pos("back_n")[0][0],get_pos("back_n")[0][1],get_pos("back_n")[0][2]])

			cmd.select("pseudo_back","pseudo_back_o or pseudo_back_c or pseudo_back_n")

			#TRANSLATE PSEUDO BACK TO ORIGIN TO GET ROTATION ANGLES FOR BACKBONE

			pos_pseudo_back = get_pos("pseudo_back_o")
			pos_pseudo_back = np.array(pos_pseudo_back)


			translation_back=-pos_pseudo_back

			cmd.translate([translation_back[0][0],translation_back[0][1],translation_back[0][2]],"pseudo_back",-1,0)




			#----START CALCULATIONS WITH PSEUDO LIGAND TO GET ANGLES FOR LATER ROTATION OF CIRCULAR COMPLEX---------------
			pos_pseudo_ligand = get_pos("pseudo_halligand")
			pos_pseudo_ligand = np.array(pos_pseudo_ligand)




			translation_vector = -pos_pseudo_ligand
			cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],"pseudo_ligand",-1,0)

			angles_circ = align_x_neg("pseudo_ligand","pseudo_ciligand")
			angles_circ.append(zero_y("pseudo_ligand","pseudo_c2ligand"))
			align_ligand("circ_complex","halligand",angles_circ)

			#rotate pseudo backbone


			angles_back = align_x_neg("pseudo_back","pseudo_back_c")
			angles_back.append(zero_y("pseudo_back","pseudo_back_n"))


			align_ligand("comp_back","back_o",angles_back)

			cmd.zoom("comp_back")

			selections = cmd.get_names("all")

			c=0

			#FINALLY ACCESS SEPARATED SELECTIONS TO ALIGN PSEUDOATOMS WITH ENERGIES TO THE BINDING POCKET

			#align backbone of ligand to backbone of complex
			align_ligand("charges_backbone","back_o",angles_back)
			cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"charges_circle",-1,0)
			cmd.rotate("z",180,"charges_circle",0,0,origin=[0,0,0])
			cmd.rotate("x",180,"charges_circle",0,0,origin=[0,0,0])

			#align halogen structure of ligand to halogen structure of complex
			align_ligand("charges_circle","halligand",angles_circ)
			cmd.select("charged_aas","resn glu+asp+arg+lys+his")

			#clean up selection no longer needed 


			cmd.hide("lines")

			cmd.show("sticks","ligand_circ")
			cmd.show("sticks","ligand_back")


			
			cmd.hide("sticks","circ_complex")
			cmd.hide("lines","circ_complex")
			cmd.hide("lines","ligand")
			cmd.delete("back_*")
			cmd.delete("hal")
			cmd.delete("polar_atoms")
			cmd.delete("cicomplex")
			cmd.delete("c2complex")
			cmd.delete("comp_back_*")
			cmd.delete("pseud*")
			cmd.delete("ligand_a")
			cmd.delete("circ_complex_only")
			cmd.delete("halligand")
			cmd.delete("c*ligand")
			cmd.hide("cgo","axes")
			#cmd.delete("comp_back")
			cmd.delete("halcomplex")
			cmd.delete("ligand")
			cmd.hide("sticks","comp_back")
			cmd.hide("lines","comp_back")
			cmd.delete("complex")
			cmd.hide("nonbonded")
			cmd.show("nonbonded","charges_circle")
			cmd.show("nonbonded","charges_backbone")
			cmd.delete("charges_circ")
			cmd.delete("charges_back")
			

			#FIND CHARGED RESIDUES IN CHARGE CLOUD
			cmd.select("charged_aas","resn glu+asp+arg+lys+his")
			cmd.delete("chain_a")
			aas = []
			com("charges_backbone",object="cen_charges_backbone")
			com("charges_circle",object="cen_charges_circle")

			for atom in cmd.get_model("charged_aas").atom:
				if atom.name == "OD1" or atom.name == "OD2" or (atom.resn == "HIS" and atom.name =="HB2") or \
				atom.name == "HZ1" or atom.name == "HH12" or atom.name == "NZ" or \
				atom.name == "NE2" or (atom.name == "N" and atom.resn == "HIS") or atom.name =="OE1" or \
				atom.name == "OE2" or atom.name == "NH2" or atom.name == "NH1":
					dist = cmd.get_distance("resi "+str(atom.resi)+" and id "+str(atom.id),"cen_charges_circle")
					#print(atom.resn,atom.resi,dist)
					if dist < 8:
						aas.append((atom.resn,atom.resi))
						

			for atom in cmd.get_model("charged_aas").atom:
 				if atom.name == "OD1" or atom.name == "OD2" or (atom.resn == "HIS" and atom.name =="HB2") or \
				atom.name == "HZ1" or atom.name == "HH12" or atom.name == "NZ" or \
				atom.name == "NE2" or (atom.name == "N" and atom.resn == "HIS") or atom.name =="OE1" or \
				atom.name == "OE2" or atom.name == "NH2" or atom.name == "NH1":
					dist = cmd.get_distance("resi "+str(atom.resi)+" and id "+str(atom.id),"cen_charges_backbone")
					#print(atom.resn,atom.resi,dist)
					if dist < 8:
						aas.append((atom.resn,atom.resi))



			cmd.delete("cen_charges*")

			aas = set(aas)
			print("Charged residues near the point cloud: ")
			print(aas)
			
			for res in aas:
				cmd.select(str(res[0])+"_"+res[1],"chain a and resn "+res[0]+" and resi "+res[1])
				cmd.show("lines",str(res[0])+"_"+res[1])

			if len(aas) > 0:
				cmd.save(str(se)+"_"+str(pdb)+"_halbond"+str(hal_bond)+"_halogen_"+str(c_hal)+".pse")
			c_hal+=1
		break
	break

		
	
		