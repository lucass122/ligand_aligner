from pymol import cmd
from os import listdir
from os.path import isfile, join

print sys.version
files = [f for f in listdir(".") if isfile(join(".", f)) and ".pdb" in f]
print(files)
cmd.reinitialize()
cmd.bg_color("white")	# set background color
for file in files:
	cmd.reinitialize()
	print("Current file: " + file)
	cmd.load(file, "protein_structure")	# load file
	cmd.select("ligand", "organic")	# selects all ligands in the file
	cmd.select("bindingsite", "ligand expand 5.0")	# expand selection by 5 angstrom, extend = extend along bonds
	cmd.select("bindingsite", "br. bindingsite")	# br extends the selection to the full residue
	cmd.select("bindingsite", "bindingsite and not ligand")	# br extends the selection to the full residue
	cmd.select("water", "resn hoh")	# resn = residues name, hoh = name of water molecules
	cmd.select("water_bs", "water and bindingsite")

	cmd.select("halogen", "e. Cl or e. Br or e. I")	# select all relevant halogen atoms
	number_of_halogens = cmd.select("halogen", "halogen and ligand")
	if number_of_halogens == 0:
		print("No halogens found in this file.")
		continue
	cmd.select("polar_atoms", "e. O or e. N or e. S")
	model_halogen = cmd.get_model("halogen")	# create model for halogen selection
	for halo in model_halogen.atom:
		cmd.select("current_halogen", "id %s"%(halo.id))
		#cmd.select("current_neighbor", "current_halogen extend 1")
		#cmd.select("current_neighbor", "current_neighbor and not current_halogen")
		cmd.select("current_neighbor", "neighbor current_halogen")

		cmd.select("current_surroundings", "current_halogen expand 5.0")
		cmd.select("current_surroundings", "current_surroundings and bindingsite")
		cmd.select("current_polar", "current_surroundings and polar_atoms")
		polar_distance = cmd.dist("polar_contacts", "ligand", "polar_atoms", mode=2)
		print("Polar contact: ")
		print polar_distance, halo.name, halo.resn, halo.resi
		model_polar = cmd.get_model("current_polar")
		for pol in model_polar.atom:
			cmd.select("current_polar", "id %s"%(pol.id))
			# check for halogen bonds, distances 5.0, angles: > 140
			current_distance = cmd.distance("current_distance", "current_halogen", "current_polar")
			current_angle = cmd.angle("current_angle", "current_neighbor", "current_halogen", "current_polar")

			if (current_angle >= 140) and (current_distance <= 5.0):
				print("Halogen bond: ")
				print current_distance, current_angle, pol.name, pol.resn, pol.resi

			cmd.delete("current_angle")
			cmd.delete("current_distance")
		








	# visualizations
	cmd.hide("lines")
	cmd.hide("nonbonded")	# hide all nonbonded atoms, i.e. oxygen of water or ions
	cmd.show("sticks", "ligand")
	cmd.show("lines", "bindingsite")
	cmd.show("nb_spheres", "water_bs")	# nb_spheres = non bonded spheres
	cmd.color("skyblue", "water_bs")
	cmd.zoom("ligand", 5)	# high numbers = far away

