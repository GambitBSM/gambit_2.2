
import sys
import os
import h5py
import numpy as np

if len(sys.argv) < 3:
    print "Usage: python ascii2hdf5 <output_hdf5_file> <path_to_ascii_files_1> [<path_to_ascii_files_2> ...]"

hdf5file = sys.argv[1]

asciipaths = []
for i in range(2,len(sys.argv)):
    asciipaths.append(sys.argv[i])

print "Reading data from files..."
data = {}
base = ''
for path in asciipaths:

    # Base file name
    base = '.'.join(os.listdir(path)[0].strip().split('.')[:-1])

    # Find number of output files
    n =  len([name for name in os.listdir(path) if os.path.isfile(os.path.join(path, name))
         and len(name.split('.')[-1].split('_')) <= 2
         and "live" not in name
         and "info" not in name
         and "txt_txt" not in name
         and "hdf5" not in name])
    print "Number of files = ", n

    # Read all files
    for i in range(0,1):

        # Filename
        datafile = path + '/' + base + '.txt_' + str(i)
        infofile = path + '/' + base + '.txt_info_' + str(i)

        # If there is only one file sometimes has _0 and sometimes it does not
        if n == 1 and not os.path.exists(datafile):
            datafile = path + '/' + base + '.txt'
            infofile = path + '/' + base + '.txt_info'

        datasets = []
        # Read off columns and dataset names
        with open(infofile) as finfo:
            for line in finfo:
                if "Column" in line:
                    col = int(line.strip().split()[1][:-1])-1
                    name = " ".join(line.strip().split()[2::])  
                    datasets.append(name)

        # Read data
        with open(datafile) as fdata:
            for line in fdata:
                linedata = line.strip().split()
                for i in range(0,len(linedata)):
                    if datasets[i] in data.keys():
                        data[datasets[i]].append(linedata[i])
                    else:
                        data[datasets[i]] = [linedata[i]]

    # Read posterior info file
    datafile = path + '/' + base + '.txt_0_txt_0'
    infofile = path + '/' + base + '.txt_0_txt_info_0'

    # If there is only one file sometimes has _0 and sometimes it does not
    if not os.path.exists(datafile):
        datafile = path + '/' + base + '.txt_txt'
        infofile = path + '/' + base + '.txt_txt_info'

    datasets = []
    # Read off columns and dataset names
    with open(infofile) as finfo:
        for line in finfo:
            if "Column" in line:
                col = int(line.strip().split()[1][:-1])-1
                name = " ".join(line.strip().split()[2::])
                datasets.append(name)

    # Read data
    with open(datafile) as fdata:
        for line in fdata:
            linedata = line.strip().split()
            for i in range(0,len(linedata)):
                if datasets[i] in data.keys():
                    data[datasets[i]].append(linedata[i])
                else:
                    data[datasets[i]] = [linedata[i]]

# Write the hdf5 file
print "Printing hdf5 file..."
with h5py.File(hdf5file, 'w') as fhdf5:
    group = fhdf5.create_group(base)

    for key in data.keys():
        group.create_dataset(key, data=np.array(data[key]))
        group.create_dataset(key+"_isvalid", data=np.full(len(data[key]), True, dtype=bool))



