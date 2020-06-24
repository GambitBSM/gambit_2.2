"""
Read log evidences from MultiNest list of output files
"""

def read_log_z(file_name):
    with open(file_name) as f:
        lines = f.readlines()

    for line in lines:
        if line.strip().startswith("ln(ev)="):
            data_str = line.split("=")[1].split("+/-")[0]
            return float(data_str.strip())
            

if __name__ == "__main__":
    
    # Vary couplings
    file_names = ["signal_{}.log".format(i) for i in range(10)]
    log_z = [read_log_z(f) for f in file_names]
    print(log_z)

    # Vary tritium
    file_names = ["bkg_tritium_{}.log".format(i) for i in range(10)]
    log_z = [read_log_z(f) for f in file_names]
    print(log_z)
