
def evol(constants, model, mode, osc, num, scan, grid, ad_output, nad_output):
	# constants
	
	# model
	model_type = 'EVOL'
	file = model[0]
	file_format = model[1]
	data_format = model[2]
	deriv_type = model[3]
	uniform_rot = model[4]
	Omega_rot = model[5]
	Omega_units = model[6]
	add_center = model[7]
	dx_snap = model[8]
	repair_As = model[9]

	# mode
	l = mode[0]
	m = mode[1]
	tag = mode[2]
	n_pg_min = mode[3]
	n_pg_max = mode[4]
	rossby = mode[5]

	# osc
	inner_bound = osc[0]
	outer_bound = osc[1]
	variables_set = osc[2]
	inertia_norm = osc[3]
	rotation_method = osc[4]
	time_factor = osc[5]
	conv_scheme = osc[6]


	# num

	# scan

	# grid

	# ad_output

	# nad_output

