import sys
import h5py
import os
from process_hdf5 import process_hdf5_freqs
import subprocess
import time
import random
import tempfile as tf
from full_shooting import generate_polytrope_file

n = float(sys.argv[1])
folder_name = sys.argv[2]

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
