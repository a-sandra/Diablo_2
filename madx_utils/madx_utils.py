# read MADX-PTC output distribution
from _dependencies import pandas as _pd
from _dependencies import numpy as _np

def get_PTC_distribution(filename):
    distr=[]
    newhead = ["NUMBER", "TURN", "X", "PX", "Y", "PY", "T", "PT", "S", "E"] # header of the file
    distr = _pd.read_csv(filename, skiprows=1, delim_whitespace=True, names=newhead)
    return distr

def get_beam_distribution(filename):
    lheader = ["X(mm)", "XP(mrad)", "Y(mm)", "YP(mrad)", "Z(mm)", "W(MeV)"]
    df = _pd.read_csv(filename, dtype=float, delim_whitespace=True, header=0, names=lheader, skiprows=1)
    get_distribution(df)
    return df

def get_number_of_particle(filename):
    distr=[]
    newhead = ["NUMBER", "TURN", "X", "PX", "Y", "PY", "T", "PT", "S", "E"] # header of the file
    distr = _pd.read_csv(filename, skiprows=1, delim_whitespace=True, names=newhead)
    return len(distr["X"])
#
def rms_distribution(filename, col):
    d = []
    d = get_PTC_distribution(filename)
    return [d["S"][0], _np.std(d[col]),  _np.mean(d[col])]

def rms_beamsize(filename, col1 = "X(mm)", col2 = "Y(mm)", col3 = "Ekin(MeV)"):
    d = []
    d = get_beam_distribution(filename)
    return [_np.mean(d[col3]),_np.std(d[col1]),  _np.std(d[col2])]

def get_twiss(filename):
    twiss = _pd.read_csv(filename, skiprows=46, delim_whitespace=True)
    new_header = twiss.columns.values[1:]
    pd1 = _pd.read_csv(filename, skiprows=48, delim_whitespace=True, names=new_header)
    return pd1

def get_nozzle_optics(filename):
    df = get_twiss(filename)
    noz = df[df["NAME"]=="NOZZLE_ENTRANCE"]
    return noz

def write_distribution(filename, ptc_start):
    print_coord = open(filename, "w")
    for line in ptc_start:
        print_coord.write(line)
        print_coord.write("\n")
    print_coord.close()

def get_ptc_twiss(file):
    ptc_file = _pd.read_csv(file, delim_whitespace=True, skiprows=88)
    new_header = ptc_file.columns.values[1:]
    return _pd.read_csv(file, delim_whitespace=True, skiprows=90, names=new_header)

#
# def search_string_list(filename):
#     ff = open(filename, "r")
#     for line in ff:
#         line = line.strip()
#         if re.match(r"Final Penalty Function =", line):
#             aline = line.split("=")
#             return [float(filename.split("_")[-2]), float(aline[1])]
#
#
# def get_gradient_from_file(filename):
#     ff = open(filename, "r")
#     return np.array([float(re.split('= |;', line)[1]) for line in ff])