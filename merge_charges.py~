from pymol import cmd


def get_pos(selection):
	pos_sel = []
	for atom in cmd.get_model(selection).atom:
		pos_sel.append([atom.coord[0], atom.coord[1], atom.coord[2]])
	return pos_sel

file_name = "session_iodine_minus1"



cmd.reinitialize()
cmd.load(file_name + ".pse", partial=1)
cmd.translate([-100,-100,-100],"complex",-1,0)
selections = cmd.get_names("all")
print(selections)
c=0
for selection in selections:
	if "x0_y8_z8" in selection:
		lim1 = c
	if "x1_y-1_z-1" in selection:
		lim2 = c
	if "x2_y-1_z-1" in selection:
		lim3 = c
	if "x9_y8_z8" in selection:
		lim4 = c

	c+=1

selections_back = selections[0:lim1] + selections[lim2:lim3]
selections_circ= selections[lim3:lim4+1]+ selections[lim1:lim2]
print(selections_circ)


cmd.select("charges_circ",selections_circ[0])
for selection in selections_circ:
	cmd.select("charges_circ",selection+  " or charges_circ")

cmd.select("charges_back",selections_back[0])
for selection in selections_back:
	cmd.select("charges_back",selection+" or charges_back")





cmd.create("charges_backbone","charges_back")
cmd.create("charges_circle","charges_circ")

cmd.translate([100,100,100],"complex",-1,0)
for selection in selections_back:
	cmd.delete(selection)

for selection in selections_circ:
	cmd.delete(selection)

cmd.save(file_name+"_merged.pse")
