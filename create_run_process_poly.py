import sys
import h5py
import os
from process_hdf5 import process_hdf5_freqs
import subprocess
import time
import random
import tempfile as tf
from full_shooting import generate_polytrope_file

def make_poly_in(l, n_pg_min, n_pg_max, freq_min, freq_max, grid_type, n_freq):
	fd, infile = tf.mkstemp()
	filename = "gyre.in"
	f = os.fdopen(fd, 'w')
	f.write('''
&constants
/

&model
	model_type = 'POLY'  ! Obtain stellar structure from an evolutionary model
	file = 'poly.h5'    ! File name of the evolutionary model
/

&mode
	    l = %d                ! Harmonic degree
        n_pg_min = %d
        n_pg_max = %d
/

&osc
        outer_bound = 'VACUUM' ! Use a zero-pressure outer mechanical boundary condition
/

&num
	diff_scheme = 'COLLOC_GL4' ! 4th-order collocation scheme for difference equations
/

&scan
        grid_type = %s ! Scan for modes using a uniform-in-period grid; best for g modes
        freq_min = %d        ! Minimum frequency to scan from
	freq_max = %d        ! Maximum frequency to scan to
	n_freq = %d          ! Number of frequency points in scan
/

&grid
	alpha_osc = 10  ! Ensure at least 10 points per wavelength in propagation regions
	alpha_exp = 2   ! Ensure at least 2 points per scale length in evanescent regions
	n_inner = 5     ! Ensure at least 5 points between center and inner turning point
/


&ad_output
        summary_file = 'summary.hdf5'                            ! File name for summary file
        summary_item_list = 'l,n_pg,omega,E_norm' ! Items to appear in summary file
        mode_template = 'mode.hdf5'                		! File-name template for mode files
        mode_item_list = 'l,n_pg,omega,x,xi_r,xi_h'   		! Items to appear in mode files
/

&nad_output
/
'''.format(l, n_pg_min, n_pg_max, grid_type, freq_min, freq_max, n_freq))	
	f.close()
	return infile




n = float(sys.argv[1])
#l = int(sys.argv[2])
folder_name = sys.argv[2]


fg = make_poly_in(2, -50, -1, 'INVERSE', 0.1, 3, 250)
fp = make_poly_in(2, 1, 35, 'LINEAR', 1, 100, 250)

poly_file_name = "poly.h5"
Gamma_1 = 1.66666666666666667
generate_polytrope_file(n, poly_file_name, Gamma_1)
# could also use build_poly, which is now present, here

x = random.randint(1, 100)
# now make the work directory, process p modes first
p_path = '$GYRE_DIR/' + folder_name + 'p' + str(x)

p_direct_path = '/Users/kmbrgandhi/gyre/' + folder_name + 'p' + str(x)
os.system("mkdir " + p_path)
os.system("cp " + poly_file_name + " " + p_path + "/")

# make g_directory
g_path = "$GYRE_DIR/" + folder_name + "g" + str(x)
g_direct_path = '/Users/kmbrgandhi/gyre/' + folder_name + 'g' + str(x)
os.system("mkdir " + g_path)
os.system("cp " + poly_file_name + " " + g_path + "/")

# process p modes
os.system("cp " + "$GYRE_DIR/pmodes_poly3/gyre.in " + p_path + "/")
os.system("$GYRE_DIR/bin/gyre " + p_path + "/" + "gyre.in")
os.system("mv summary.hdf5 " + p_path +"/")
os.system("find . -type f -name mode\* -exec rm {} \;")
f= h5py.File(p_direct_path + '/summary.hdf5')
freqs = f['omega']
ns = f['n_pg']
angular_frequency = 99.855377 # used for comparison to tabulated results
results_p = [(ns[i], angular_frequency * freqs[i][0]) for i in range(len(freqs))]
print(results_p)

os.system("cp " + "$GYRE_DIR/gmodes_poly3/gyre.in " + g_path + "/")
os.system("$GYRE_DIR/bin/gyre " + g_path + "/" + "gyre.in")
os.system("mv summary.hdf5 " + g_path +"/")
os.system("find . -type f -name mode\* -exec rm {} \;")
f2= h5py.File(g_direct_path + '/summary.hdf5')
freqs = f2['omega']
ns = f2['n_pg']
results_g = [(ns[i], angular_frequency * freqs[i][0]) for i in range(len(freqs))]
print(results_g)

os.system('rm poly.h5')
os.system('rm -r $GYRE_DIR/' + folder_name + 'g' + str(x))
os.system('rm -r $GYRE_DIR/' + folder_name + 'p' + str(x))

