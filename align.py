from pymol import cmd
import numpy as np
import pymol
from os.path import isfile, join
from os import listdir

#Scan through pdb files in folder pdbs
files = [f for f in listdir("pdbs/") if ".pdb" in f and not ".pse" in f]


#INPUT PDB NAME HERE - do this for testing individual pdbs. you will have to change the code for loading the pdbs later on as well (for loop over files)
pdb = "5osw"
#INPUT WHICH HALOGEN BOND TO USE FOR ALIGNMENT IF THERE ARE MULTIPLE 
#0 IS ALWAYS THE ONE WITH THE SHORTEST DISTANCE BETWEEN THE POLAR ATOM AND THE HALOGEN
hal_bond = 0
#INPUT SESSION NAME

se ="session_chlorine_minus1_merged"

if "minus1" in se:
	se_charge = "negative"
if "plus1" in se:
	se_charge = "positive"


"""
#Function that gets the position of a pymol selection and returns it as a python list
@param selection: String containing selection of interest
"""

def get_pos(selection):
	pos_sel = []
	for atom in cmd.get_model(selection).atom:
		pos_sel.append([atom.coord[0], atom.coord[1], atom.coord[2]])
	return pos_sel

"""
Function that rotates a selection around the given angle until the reference atom (zero) lies on a plain.
@param selection: Selection to rotate
@param axis: Axis to rotate around
@param angle: Angle to rotate around
@param zero: Selection containing an atom to determine to direction of rotation. This atom will lie on a certain axis after rotation
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

@param selection: Selection to rotate
@param ref_atom: Atom which x coordinate should be 0
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

@param selection: Selection to rotate
@param ref_atom: Atom which y coordinate should be 0
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

@param selection: Selection to rotate
@param ref_atom: Atom which z coordinate should be 0
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

@param selection: Selection to align with the negative x-axis
@param ref_atom: Atom used for calculating direction of rotation
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
The return value is a list sorted by the distance between halogen and polar atom
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
				halogen_bonds.append((pol.resn,pol.resi,pol.id,current_distance,current_angle))
				#Append fitting halogen bonds to list
				

			cmd.delete("current_angle")
			cmd.delete("current_distance")
			
			
	#Clear up selections
	cmd.delete("current_halogen")
	cmd.delete("current_neighbor")
	cmd.delete("water")
	cmd.delete("water_bs")
	cmd.delete("current_polar")
	cmd.delete("current_surroundings")
	cmd.delete("bindingsite1")

	#Sort the list by distance
	halogen_bonds = sorted(halogen_bonds, key=lambda x: x[3])

	return halogen_bonds

"""
Function that performs the necessary rotations to align a selection with a molecule. It uses the angles computed by aligning pseudo atoms with
a selection.

@param selection: selection to align
@param center: Atom which the selection is translated to. 
@param angles: Angles that have to be computed by aligning pseudoatoms to the selection
"""

def align_ligand(selection,center,angles):
	pos_halligand = get_pos(center)
	translation_vector = pos_halligand

	cmd.rotate("x",-angles[3],selection,0,0,origin=[0,0,0])
	cmd.rotate("y",-angles[2],selection,0,0,origin=[0,0,0])
	cmd.rotate("x",-angles[1],selection,0,0,origin=[0,0,0])
	cmd.rotate("z",-angles[0],selection,0,0,origin=[0,0,0])

	cmd.translate([translation_vector[0][0],translation_vector[0][1],translation_vector[0][2]],selection,-1,0)
	
"""
Function that determines the effect of a charged residue on the halogen bond.
@param atom: model containing atom of charged residue
"""
def rate_residues(atom):
	cmd.select("charges_near_res","not resn " + str(atom.resn)+ " and(resn "+str(atom.resn)+" and id "+str(atom.id)+ " expand 0.5)")
	cmd.iterate('charges_near_res', 'pymol.color_list.append(color)')
	print(cmd.get_color_tuple(pymol.color_list[0]))
	pos = False
	neg = False
	depends = False
	for color in pymol.color_list:
		rgb = cmd.get_color_tuple(color)
		if rgb[0] > 0.7 or rgb[2] > 0.7:
			if atom.resn == "LYS" or atom.resn == "ARG":
				if se_charge == "negative":
					pos = True
				else:
					pos = False
			if atom.resn == "GLU" or atom.resn == "ASP":
				if se_charge == "positive":
					neg = True
				else:
					neg = False
			if atom.resn == "HIS":
				depends = True

	if pos:	
		effect = "Positive Effect"
	if neg:
		effect = "Negative Effect"
	if depends:
		effect = "Depends on pH value"
	if not pos and not neg and not depends:
		effect = "Neutral Effect"
	return effect



"""
This calls reinitalize function to reset the screen after each iteration, loads important molecules as well as basic selections
"""
def init():
	cmd.reinitialize()
	#LOAD FILES AND MAKE THE NECESSARY SELECTIONS FOR ALIGNMENT
	cmd.load("pses/session_chlorine_minus1_merged.pse")
	cmd.do("run axes.py")
	cmd.do("run center_of_mass.py")
	#cmd.load("geometry_complex_Cl.mol")
	cmd.select("circ_complex","organic")
	
	cmd.fetch("pdbs/"+pdb)
	cmd.select("water", "resn hoh")
	#select ligand and bindingsite
	cmd.select("ligand","organic and not circ_complex")
	cmd.select("chain_a","chain a")
	cmd.select("ligand_a", "ligand and chain_a")


#START OF SCRIPT MAIN BODY


#for loop over all pdb files in pdbs folder
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

		cmd.select("halligand", "id "+ str(halligand) +  " and ligand and (name I* or name Br* or name Cl*)")
		cmd.select("hal", "(e. Cl or e. Br or e. I) and (halligand expand 5)") 

		#find the important bonds for the alignment
		hal_bonds = find_halogen_bonds()
		#check if there are no halogen bonds in pdb and go to next file if it is the case
		if len(hal_bonds) == 0:
			break
			continue
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
			#Case 1: the polar atom is a sulphur
			if "S" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","first bound_to back_o")
				cmd.select("back_n","last bound_to back_o")
			#Case 2: the polar atom is an oxygen
			if "O" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","bound_to back_o")
				cmd.select("back_n","first bound_to back_c")
			#Case 3: the polar atom is a nitrogen
			if "N" in model_back_polar[0].name:
				cmd.select("back_o","resi "+str(hal_bonds[hal_bond][1])+" and id "+ str(hal_bonds[hal_bond][2]))
				cmd.select("back_c","first bound_to back_o")
				cmd.select("back_n","first bound_to back_c")

			cmd.select("ligand_circ","br. halligand")
			cmd.select("ligand_back","br. back_c")

			#select atoms of halogen structure of complex
			cmd.select("circ_complex","circ_complex_only")
			cmd.select("ciligand", "ligand and name C* and bound_to halligand")
			cmd.select ("cicomplex", "name C* and bound_to halcomplex") 

			translation_vector = 0
			#get positions of ligand and complex atoms previously selected
			pos_halcomplex = get_pos("halcomplex")
			pos_halligand =get_pos("halligand")
			pos_cicomplex = get_pos("cicomplex")
			pos_ciligand = get_pos("ciligand")
			#turn position into numpy array for easier calculations
			pos_cicomplex = np.array(pos_cicomplex)
			pos_ciligand = np.array(pos_ciligand)
			pos_halcomplex = np.array(pos_halcomplex)

			cmd.select ("c2complex", "first bound_to cicomplex and name c*")
			cmd.select ("c2ligand", "first not halligand and bound_to ciligand" )
			translation_vector1 = -pos_halcomplex

			#translate molecule from mol file to 0,0,0
			cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"circ_complex",-1,0)

			#get coordinates from ligand and mol file complex
			pos_halligand = get_pos("halligand")
			pos_halcomplex = get_pos("halcomplex")
			pos_halligand = np.array(pos_halligand)

			#rotate the complex to fit the negative x axis
			cmd.rotate("z",180,"circ_complex",0,0,origin=[0,0,0])
			cmd.rotate("x",180,"circ_complex",0,0,origin=[0,0,0])
			#do the same for the backbone of the complex
			cmd.rotate("x",180,"comp_back",0,0,origin=[0,0,0])
			pos_c2ligand = get_pos("c2ligand")

			#create pseudoatoms of the ligand for calculating rotation angles without moving the actual ligand
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
			#get the rotation angles for the halogen complex and rotate pseudo ligand
			angles_circ = align_x_neg("pseudo_ligand","pseudo_ciligand")
			angles_circ.append(zero_y("pseudo_ligand","pseudo_c2ligand"))
			align_ligand("circ_complex","halligand",angles_circ)

			#get the rotation angles for the ligand backbone and rotate pseudo backbone
			angles_back = align_x_neg("pseudo_back","pseudo_back_c")
			angles_back.append(zero_y("pseudo_back","pseudo_back_n"))

			#align backbone of ligand to backbone of complex
			align_ligand("comp_back","back_o",angles_back)


			#ALIGN PSEUDOATOMS WITH ENERGIES TO THE BINDING POCKET
			#align the backbone charges
			align_ligand("charges_backbone","back_o",angles_back)
			cmd.translate([translation_vector1[0][0],translation_vector1[0][1],translation_vector1[0][2]],"charges_circle",-1,0)
			#align the charges around the halogen complex - NOTE: first rotate the charges like it has been done with the pseudoatoms earlier
			cmd.rotate("z",180,"charges_circle",0,0,origin=[0,0,0])
			cmd.rotate("x",180,"charges_circle",0,0,origin=[0,0,0])

			align_ligand("charges_circle","halligand",angles_circ)
			cmd.select("charged_aas","resn glu+asp+arg+lys+his")

			#clean up selections no longer needed 
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
			cmd.delete("halcomplex")
			cmd.delete("ligand")
			cmd.hide("sticks","comp_back")
			cmd.hide("lines","comp_back")
			cmd.hide("nonbonded")
			cmd.show("nonbonded","charges_circle")
			cmd.show("nonbonded","charges_backbone")
			cmd.delete("charges_circ")
			cmd.delete("charges_back")
			cmd.zoom("charges_circle",2)
			cmd.hide("lines","complex")
			#set background color to white
			cmd.bg_color("white")

			
			#FIND CHARGED RESIDUES IN CHARGE CLOUD
			cmd.select("charged_aas","resn glu+asp+arg+lys+his")
			cmd.delete("chain_a")
			aas = []
			#create pseudoatom for each point cloud representing each center of mass
			com("charges_backbone",object="cen_charges_backbone")
			com("charges_circle",object="cen_charges_circle")
			#if there are no charged residues go to the next file
			if len(cmd.get_model("charged_aas").atom) == 0:
				break
				continue
			print()
			#if len(cmd.get_model("charged_aas")) == 0:
			#	print(cmd.get_model("charged_aas").atom)
			#	continue
			for atom in cmd.get_model("charged_aas").atom:
				if atom.name == "OD1" or atom.name == "OD2" or (atom.resn == "HIS" and atom.name =="HB2") or \
				atom.name == "HZ1" or atom.name == "HH12" or atom.name == "NZ" or \
				atom.name == "NE2" or (atom.name == "N" and atom.resn == "HIS") or atom.name =="OE1" or \
				atom.name == "OE2" or atom.name == "NH2" or atom.name == "NH1":
					try:
						dist1 = cmd.get_distance("resi "+str(atom.resi)+" and id "+str(atom.id),"cen_charges_circle")
						dist2 = cmd.get_distance("resi "+str(atom.resi)+" and id "+str(atom.id),"cen_charges_backbone")
					except pymol.CmdException:
							print("Couldn't fetch distance")
							break
							continue
					if dist1 < 8:
						#get the effect of the residue on the halogen bond
						effect = rate_residues(atom)
						aas.append((atom.resn,atom.resi,"Backbone cloud",effect))

					if dist2 < 8:
						#get the effect of the residue on the halogen bond
						effect = rate_residues(atom)
						aas.append((atom.resn,atom.resi,"Halogen complex cloud",effect))



			#delete selections no longer necessary	
			cmd.delete("cen_charges*")
			cmd.delete("charges_near_res")
			aas = set(aas)
			print("Charged residues near the point cloud: ")
			print(aas)
			

			filename =   str(se)+"_"+str(pdb)+"_halbond"+str(hal_bond)+"_halogen_"+str(c_hal)

			#Save session information in separate text file in session_info folder
			with open("session_info/"+filename+".txt","w") as file: 
				file.write("File with charge cloud: "+str(se)+ "\n")
				file.write("Used pdb: "+str(pdb)+ "\n")
				file.write("Index of halogen bond used for aligning: "+str(hal_bond)+ "\n")
				file.write("Index of halogen used for aligning: "+str(c_hal)+ "\n")
				file.write("Atom ID of halogen used for aligning: "+str(halligand)+ "\n")
				file.write("resn and resi of polar atom in halogen bond: "+ str(hal_bonds[hal_bond][0]) + str(hal_bonds[hal_bond][1])+ "\n")
				file.write("List of residues in point cloud: "+str(aas)+ "\n")



			#select residues in charge cloud and show them
			for res in aas:
				cmd.select(str(res[0])+"_"+res[1],"chain a and resn "+res[0]+" and resi "+res[1])
				cmd.show("lines",str(res[0])+"_"+res[1])



			if len(aas) > 0:
				cmd.save("pses/"+filename +".pse")
			c_hal+=1
			#END OF SCRIPT



