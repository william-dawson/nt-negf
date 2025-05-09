"""
Generate data for a test calculation
"""
from BigDFT.Calculators import SystemCalculator
from BigDFT.Systems import System
from BigDFT.Fragments import Fragment
from BigDFT.Atoms import Atom, AU_to_A
from BigDFT.Inputfiles import Inputfile
from BigDFT.UnitCells import UnitCell
from BigDFT.PostProcessing import BigDFTool
from common import NEGFInterop as negf
from scipy.sparse import coo_matrix
from copy import deepcopy

# make the system
spacing = 1.4
region = 4
name = "sw"

sys = System()
sys["LEF:2"] = Fragment()
sys["LEF:1"] = Fragment()
sys["MID:0"] = Fragment()
sys["RIG:1"] = Fragment()
sys["RIG:2"] = Fragment()

i = 0
for fragid, frag in sys.items():
    rlen = region
    if fragid == "MID:0":
        rlen = region * 5
    for _ in range(rlen):
        at = Atom({"C": [0, 0, i * spacing], "units": "angstroem"})
        frag.append(at)
        i += 1

sys.cell = UnitCell([float("inf"), float("inf"), spacing * i],
                    units="angstroem")

# run the calculation
inp = Inputfile()
inp.set_xc("PBE")
inp.set_hgrid(0.5)
inp.set_psp_nlcc()
inp["import"] = "linear"

calc = SystemCalculator(skip=True)
log = calc.run(sys=sys, input=inp, name=name, run_dir="scr")

# get the ordering and write it to file
tool = BigDFTool()
order, idx = negf.get_ordering(sys, log, tool)

lengths = [len(idx["MID:0"])]
for key in ["LEF:1", "LEF:2", "RIG:1", "RIG:2"]:
    lengths.append(lengths[-1] + len(idx[key]))
with open(f"scr/data-{name}/order.txt", "w") as ofile:
    ofile.write(str(len(order)) + "\n")
    ofile.write(" ".join([str(x) for x in order]) + "\n")
    ofile.write(" ".join([str(x) for x in lengths]) + "\n")

