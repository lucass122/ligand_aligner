import pymol
import numpy as np




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
Fucntion that aligns the selection to the negative x-axis of the origin.
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

def find_halogen_bonds():
	cmd.select("halogen", "e. Cl or e. Br or e. I")	# select all relevant halogen atoms
	number_of_halogens = cmd.select("halogen", "halogen and ligand")
	if number_of_halogens == 0:
		print("No halogens found in this file.")
		return

	cmd.select("polar_atoms", "e. O or e. N or e. S")
	model_halogen = cmd.get_model("halogen")	# create model for halogen selection
	for halo in model_halogen.atom:
		cmd.select("current_halogen", "id %s"%(halo.id))
		#cmd.select("current_neighbor", "current_halogen extend 1")
		#cmd.select("current_neighbor", "current_neighbor and not current_halogen")
		cmd.select("current_neighbor", "neighbor current_halogen")
		cmd.select("bindingsite", "halligand expand 8.0")	# expand selection by 5 angstrom, extend = extend along bonds
		cmd.select("bindingsite", "br. bindingsite")	# br extends the selection to the full residue
		cmd.select("bindingsite", "bindingsite and not ligand")	# br extends the selection to the full residue
		cmd.select("water", "resn hoh")	# resn = residues name, hoh = name of water molecules
		cmd.select("water_bs", "water and bindingsite")

		cmd.select("current_surroundings", "current_halogen expand 5.0")
		cmd.select("current_surroundings", "current_surroundings and bindingsite")
		cmd.select("current_polar", "current_surroundings and polar_atoms")
		model_polar = cmd.get_model("current_polar")
		for pol in model_polar.atom:
			cmd.select("current_polar", "id %s"%(pol.id))
			# check for halogen bonds, distances 5.0, angles: > 140
			current_distance = cmd.distance("current_distance", "current_halogen", "current_polar")
			current_angle = cmd.angle("current_angle", "current_neighbor", "current_halogen", "current_polar")
		
			if (current_angle >= 140) and (current_distance <= 5.0):
				print("Halogen bond: ")
				print current_distance, current_angle, pol.name, pol.resn, pol.resi
				cmd.select(str(pol.resn)+"_"+str(pol.resi),"resn "+str(pol.resn)+" and resi " + str(pol.resi))

			cmd.delete("current_angle")
			cmd.delete("current_distance")


	

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



cmd.reinitialize()
#LOAD FILES AND MAKE THE NECESSARY SELECTIONS FOR ALIGNMENT
#TODO: REPLACE ITERATE STATE STUFF WITH FUNCTION CALL
cmd.load("session_iodine_minus1_merged.pse", partial=1)
cmd.do("run axes.py")
#cmd.load("geometry_complex_Cl.mol")
cmd.select("circ_complex","organic")

cmd.fetch("1gjd")
cmd.select("ligand","organic and 1gjd")
cmd.select("chain_a","chain a")
cmd.select("ligand_a", "ligand and chain_a")

cmd.select("halligand", "last ligand and (name I* or name Br* or name Cl*)")
cmd.select("halcomplex", "first circ_complex and (name I* or name Br* or name Cl*)")
cmd.select("circ_complex_only","halcomplex extend 5")
cmd.select("comp_back","circ_complex and not circ_complex_only")
cmd.select("comp_back_o","comp_back and name o")
cmd.select("comp_back_c","bound_to comp_back_o")
cmd.select("comp_back_n","(bound_to comp_back_c) and name N")

cmd.select("circ_complex","circ_complex_only")
cmd.select("ciligand", "ligand and name C* and bound_to halligand")
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
cmd.select ("c2ligand", "first bound_to ciligand and name c*" )

cmd.iterate_state(1, 'c2ligand', 'pos_c2ligand.append((x,y,z))')


print(abs(pos_cicomplex))
translation_vector1 = -pos_halcomplex


#translate molecule from mol file to 0,0,0
cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"circ_complex",-1,0)

#get coordinates from ligand and mol file complex
pos_halligand = []
pos_halcomplex = []
cmd.iterate_state(1, 'halligand', 'pos_halligand.append((x,y,z))')
cmd.iterate_state(1, 'halcomplex', 'pos_halcomplex.append((x,y,z))')
pos_halligand = np.array(pos_halligand)

#add pseudoatom on x axis to measure angle for rotation around y axis

cmd.rotate("z",180,"circ_complex",0,0,origin=[0,0,0])
cmd.rotate("x",180,"circ_complex",0,0,origin=[0,0,0])

cmd.rotate("x",180,"comp_back",0,0,origin=[0,0,0])

cmd.pseudoatom("pseudo_ciligand",pos=[pos_ciligand[0][0],pos_ciligand[0][1],pos_ciligand[0][2]])
cmd.pseudoatom("pseudo_c2ligand",pos=[pos_c2ligand[0][0],pos_c2ligand[0][1],pos_c2ligand[0][2]])
cmd.pseudoatom("pseudo_halligand",pos=[pos_halligand[0][0],pos_halligand[0][1],pos_halligand[0][2]])
cmd.select("pseudo_ligand","pseudo_*ligand")

#ADD PSEUDO BACKBONE

cmd.select("back_o","first (halligand expand 4.0 and name O and bound_to name C)")
cmd.select("back_c","bound_to back_o")
cmd.select("back_n","(bound_to back_c) and name N")

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

#find_halogen_bonds()
cmd.zoom("comp_back")


selections = cmd.get_names("all")

c=0

#FINALLY ACCESS SEPARATED SELECTIONS TO ALIGN PSEUDOATOMS WITH ENERGIES TO THE BINDING POCKET
#cmd.select("charges_back",selections_back[0])

#cmd.select("charges_back",selection+" or charges_back")
align_ligand("charges_backbone","back_o",angles_back)
	
#cmd.select("charges_circ",selections_circ[0])
#for selection in selections_circ:
#cmd.select("charges_circ",selection+  " or charges_circ")
cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"charges_circle",-1,0)
cmd.rotate("z",180,"charges_circle",0,0,origin=[0,0,0])
cmd.rotate("x",180,"charges_circle",0,0,origin=[0,0,0])

align_ligand("charges_circle","halligand",angles_circ)
	
cmd.select("charged_aas","resn glu+asp+arg+lys+his")

#CLEAN UP SELECTIONS
cmd.select("ligand_circ","br. halligand")
cmd.select("ligand_back","br. back_o")
cmd.show("sticks","ligand_circ")
cmd.show("sticks","ligand_back")
cmd.hide("sticks","circ_complex")
cmd.hide("lines","circ_complex")
cmd.hide("sticks","ligand")
cmd.hide("lines","ligand")
cmd.delete("back_*")
#cmd.delete("*ligand")
cmd.delete("*complex")
cmd.delete("comp_back_*")
cmd.delete("pseud*")
cmd.delete("*_a")
cmd.delete("circ_complex_only")
cmd.delete("halligand")
cmd.delete("c*ligand")
cmd.delete("comp_back")

#FIND CHARGED RESIDUES IN CLOUD

pos_charges_backbone = get_pos("charges_backbone")
pos_charges_circle = get_pos("charges_circle")

pos_charges_backbone = np.array(pos_charges_backbone)
pos_charges_circle = np.array(pos_charges_circle)

print(pos_charges_backbone[:,0])
#get start and end positions of cloud to search for charged residues in there
low_x_1 = min(pos_charges_backbone[:,0])
low_y_1 = min(pos_charges_backbone[:,1])
low_z_1 = min(pos_charges_backbone[:,2])
up_x_1 = max(pos_charges_backbone[:,0])
up_y_1 = max(pos_charges_backbone[:,1])
up_z_1 = max(pos_charges_backbone[:,2])

low_x_2 = min(pos_charges_circle[:,0])
low_y_2 = min(pos_charges_circle[:,1])
low_z_2 = min(pos_charges_circle[:,2])
up_x_2 = max(pos_charges_circle[:,0])
up_y_2 = max(pos_charges_circle[:,1])
up_z_2 = max(pos_charges_circle[:,2])

print(low_x_1,up_x_1)
print(low_y_1,up_y_1)
print(low_z_1,up_z_1)

print(low_x_2,up_x_2)
print(low_y_2,up_y_2)
print(low_z_2,up_z_2)


cmd.pseudoatom("a",pos=[pos_charges_backbone[5][0],pos_charges_backbone[5][1],pos_charges_backbone[5][2]])
cmd.pseudoatom("ba",pos=[low_x_2,low_y_2,low_z_2])
cmd.select("charged_aas","resn glu+asp+arg+lys+his")

aas = []


for atom in cmd.get_model("charged_aas").atom:
	if ((atom.coord[0] < up_x_2 and atom.coord[0] > low_x_2) and \
	(atom.coord[1] < up_y_2 and atom.coord[1] > low_y_2) and \
	(atom.coord[2] < up_z_2 and atom.coord[2] > low_z_2)):
		#print(atom.coord[0],atom.coord[1],atom.coord[2])
		aas.append((atom.resn,atom.resi))

for atom in cmd.get_model("charged_aas").atom:
	if ((atom.coord[0] < up_x_1 and atom.coord[0] > low_x_1) and \
	(atom.coord[1] < up_y_1 and atom.coord[1] > low_y_1)  and \
	(atom.coord[2] < up_z_1 and atom.coord[2] > low_z_1)):
		#print(atom.coord[0],atom.coord[1],atom.coord[2])
		aas.append((atom.resn,atom.resi))
aas = set(aas)
for res in aas:
	cmd.select(str(res[0])+"_"+res[1],"resn "+res[0]+" and resi "+res[1])


for atom in cmd.get_model("charges_backbone").atom:
	cmd.select("a","charges_backbone and id "+str(atom.id)+" expand 1.0 and charged_aas")


"""
(8.9233989715576172, 30.156494140625)
(-2.7247285842895508, 19.904872894287109)
(6.2987575531005859, 28.963228225708008)

(-2.6179275512695312, 21.668624877929688)
(-1.466641902923584, 16.399503707885742)
(14.40446662902832, 37.706699371337891)
"""
