cmd.reinitialize
def rotate_sel_ligand(selection,axis,angle,zero):
	position = get_pos(zero)
	print(angle)
	print(position)
	if axis == "x":
		if position[0][1] >= 0.0 and position[0][2] >= 0:
			cmd.rotate("x",-(180-angle),selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(-angle) 
				+ " degree.")
			return -angle
		if position[0][1] < 0.0 and position[0][2] >= 0:
			cmd.rotate("x",180-angle,selection,0,0,origin=[0,0,0])
			print("Rotate selection "+str(selection)+" around "+ str(axis)+" axis by " + str(angle) 
				+ " degree.")
			return angle
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
