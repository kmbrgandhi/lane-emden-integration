import sys
import h5py
import os
from process_hdf5 import process_hdf5_freqs
import subprocess
import time
import random
import tempfile as tf


def build_poly (n_poly, Delta_d, xi_d, Gamma_1, dxi, toler, filename):

    n_poly_str = ','.join('{0:24.16e}'.format(n) for n in n_poly)
    Delta_d_str = ','.join('{0:24.16e}'.format(d) for d in Delta_d)
    xi_d_str = ','.join('{0:24.16e}'.format(x) for x in xi_d)

    # Create an input file



    fd, infile = tf.mkstemp()

    f = os.fdopen(fd, 'w')

    f.write('''
&poly
    n_d = {0:d}
    n_poly = {1:s}
        Delta_d = {2:s}
        xi_d = {3:s}
        Gamma_1 = {4:24.16e}
/

&num
    dxi = {5:24.16e}
    toler = {6:24.16e}
/

&out
    filename = '{7:s}'
/
'''.format(len(n_poly)-1, n_poly_str, Delta_d_str, xi_d_str,
           Gamma_1, dxi, toler, filename))

    f.close()

    # Run build_poly
    

    os.system('/Users/kmbrgandhi/gyre/bin/build_poly {0:s}'.format(infile))

    # Delete the input file

#    print(infile)

    os.remove(infile)

n = int(sys.argv[1])
folder_name = sys.argv[2]

# build poly here
poly_file_name = "poly.h5"
Gamma_1 = 1.66666666666666667
#build_poly([n], [], [], Gamma_1, 0.00244949, 1E-10, poly_file_name)





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
time.sleep(1)
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
time.sleep(1)
f2= h5py.File(g_direct_path + '/summary.hdf5')
freqs = f2['omega']
ns = f2['n_pg']
angular_frequency = 99.855377 # used for comparison to tabulated results
results_g = [(ns[i], angular_frequency * freqs[i][0]) for i in range(len(freqs))]
print(results_g)
