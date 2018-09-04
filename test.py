from pymol import cmd


def get_pos(selection):
	pos_sel = []
	for atom in cmd.get_model(selection).atom:
		pos_sel.append([atom.coord[0], atom.coord[1], atom.coord[2]])
	return pos_sel


up_x_1 = get_pos("x-9_y8_z8_bin1_0to-5_negative_minus_1_energy_-1.46")[0][0]
up_y_1 = get_pos("x-9_y8_z8_bin1_0to-5_negative_minus_1_energy_-1.46")[0][1]
up_z_1 =get_pos("x-9_y-8_z-8_bin1_0to-5_negative_minus_1_energy_-1.39")[0][2]

low_x_1= get_pos("x1_y8_z-8_bin1_0to-5_negative_minus_1_energy_-0.63")[0][0]
low_y_1 = get_pos("x1_y-8_z-8_bin1_0to-5_negative_minus_1_energy_-0.62")[0][1]
low_z_1 = get_pos("x-9_y8_z8_bin1_0to-5_negative_minus_1_energy_-1.46")[0][2]


up_x_2 = get_pos("x2_y8_z-8_bin1_0to-5_negative_minus_1_energy_-0.35")[0][0]
up_y_2 = get_pos("x16_y8_z8_positive_minus_1_energy_1.36")[0][1]
up_z_2 =get_pos("x2_y8_z8_bin1_0to-5_negative_minus_1_energy_-0.38")[0][2]

low_x_2= get_pos("x17_y-8_z-8_positive_minus_1_energy_1.28")[0][0]
low_y_2 = get_pos("x16_y8_z8_positive_minus_1_energy_1.36")[0][1]
low_z_2 = get_pos("x17_y8_z-8_positive_minus_1_energy_1.28")[0][2]

print(get_pos("x-9_y-8_z-8_bin1_0to-5_negative_minus_1_energy_-1.39"))
print(get_pos("x-9_y8_z-8_bin1_0to-5_negative_minus_1_energy_-1.39"))

print(get_pos("x-9_y8_z8_bin1_0to-5_negative_minus_1_energy_-1.46"))
print(get_pos("x-9_y-8_z8_bin1_0to-5_negative_minus_1_energy_-1.45"))

print(get_pos("x1_y-8_z8_bin1_0to-5_negative_minus_1_energy_-0.66"))
print(get_pos("x1_y-8_z-8_bin1_0to-5_negative_minus_1_energy_-0.62"))

print(get_pos("x1_y8_z-8_bin1_0to-5_negative_minus_1_energy_-0.63"))
print(get_pos("x1_y8_z8_bin1_0to-5_negative_minus_1_energy_-0.68"))


print(get_pos("x2_y8_z8_bin1_0to-5_negative_minus_1_energy_-0.38"))
print(get_pos("x16_y8_z8_positive_minus_1_energy_1.36"))

print(get_pos("x2_y8_z8_bin1_0to-5_negative_minus_1_energy_-0.38"))

print(get_pos("x17_y-8_z8_positive_minus_1_energy_1.3"))

print(get_pos("x17_y-8_z-8_positive_minus_1_energy_1.28"))

print(get_pos("x17_y8_z-8_positive_minus_1_energy_1.28"))

print(get_pos("x2_y8_z-8_bin1_0to-5_negative_minus_1_energy_-0.35"))

print(get_pos("x2_y8_z8_bin1_0to-5_negative_minus_1_energy_-0.38"))


cmd.select("charged_aas","resn glu+asp+arg+lys+his")

aas = []


for atom in cmd.get_model("charged_aas").atom:
	if ((atom.coord[0] < up_x_1 and atom.coord[0] > low_x_1) or (atom.coord[0] < up_x_2 and atom.coord[0] > low_x_2)) and \
	((atom.coord[1] < up_y_1 and atom.coord[1] > low_y_1) or (atom.coord[1] < up_y_2 and atom.coord[1] > low_y_2)) and \
	((atom.coord[2] < up_z_1 and atom.coord[2] > low_z_1) or (atom.coord[2] < up_z_2 and atom.coord[2] > low_z_2)):
		aas.append((atom.resn,atom.resi))

aas = set(aas)
print(aas)

for res in aas:
	cmd.select(str(res[0])+"_"+res[1],"resn "+res[0]+" and resi "+res[1])