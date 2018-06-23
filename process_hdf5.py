import sys
import h5py

def process_hdf5_freqs(hdf5_file):
	f = h5py.File(hdf5_file)
	freqs = f['omega']
	ns = f['n_pg']
	angular_frequency = 99.855377 # used for comparison to tabulated results

	results = [(ns[i], angular_frequency * freqs[i][0]) for i in range(len(freqs))]
	print(results)
	return results

#process_hdf5_freqs(sys.argv[1])