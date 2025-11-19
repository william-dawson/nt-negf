"""
Basically, this will write a file that is input for our libNEGF driver.
It includes information like the chemical potential, the indices of the
different regions, etc.

Usage
"""

from os.path import join
from sys import argv

from BigDFT.Fragments import Fragment
from BigDFT.IO import read_xyz, write_pdb
from BigDFT.Logfiles import Logfile
from BigDFT.PostProcessing import BigDFTool
from BigDFT.Systems import System
from common import NEGFInterop as negf


def get_mu(log):
    vals = []
    gso = log.log["Ground State Optimization"]
    for step in gso:
        if "kernel optimization" not in step:
            continue
        for opt in step["kernel optimization"]:
            if "Kernel update" not in opt:
                continue
            summary = opt["Kernel update"]["Kernel calculation"][-1]["summary"]
            vals.append(summary["eF"])
    return vals


if __name__ == "__main__":
    # Process the command line arguments
    geom = argv[1]
    log = argv[2]
    central = int(argv[3])
    left2 = int(argv[4])
    left1 = central + (left2 - central) // 2
    right2 = int(argv[5])
    right1 = left2 + (right2 - left2) // 2

    # Split the System into Regions
    with open(geom) as ifile:
        atsys = read_xyz(ifile)
    sys = System()
    sys["LEF:1"] = sum(list(atsys.values())[central:left1])
    sys["LEF:2"] = sum(list(atsys.values())[left1:left2])
    sys["MID:0"] = sum(list(atsys.values())[:central])
    sys["RIG:1"] = sum(list(atsys.values())[left2:right1])
    sys["RIG:2"] = sum(list(atsys.values())[right1:right2])
    from BigDFT.IO import write_pdb

    with open("test.pdb", "w") as ofile:
        write_pdb(sys, ofile)
    sys.cell = atsys.cell

    # get the chemical potential from the log
    log = Logfile(log)
    mu = get_mu(log)[-1]

    # get the order of the basis functions
    tool = BigDFTool()
    order, idx = negf.get_ordering(sys, log, tool)
    lengths = [len(idx["MID:0"])]
    for key in ["LEF:1", "LEF:2", "RIG:1", "RIG:2"]:
        lengths.append(lengths[-1] + len(idx[key]))
    with open(join(log.srcdir, log.data_directory, "order.txt"), "w") as ofile:
        ofile.write(str(len(order)) + "\n")
        ofile.write(" ".join([str(x) for x in order]) + "\n")
        ofile.write(" ".join([str(x) for x in lengths]) + "\n")
        ofile.write(str(mu) + "\n")
