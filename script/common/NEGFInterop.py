"""
Utilities for interactions with LibNEGF
https://github.com/libnegf/libnegf
"""

def permute_matrix(a, p, inv=False):
    """
    Permute a matrix back and forth based on a given ordering.

    Args:
        a (matrix): the matrix to rearrange.
        p (list): the order to reorder into.
        inv (bool): if True, we do the reverse ordering.
    """
    from numpy import array
    p = array(p)
    if inv:
        p = p.argsort()
    return a[p][:, p]

def get_ordering(sys, log, tool):
    """
    Get the indices we will use to match LibNEGF's internal ordering.

    Args:
        sys (BigDFT.Systems.System): the system being studied. It must be
        setup into the right set of fragments (MID:0, LEF:1, LEF:2, RIG:1,
        RIG:2).
        log (BigDFT.Logfiles.Logfile): the log of a calculation.
        tool (BigDFT.PostProcessing.BigDFTool): tool for post processing.

    Returns:
        order (list): the listed order of indices.
        idx (dict): the fragment indices. 
    """
    idx = tool.get_frag_indices(sys, log)
    order = idx["MID:0"] + idx["LEF:1"] + idx["LEF:2"] + \
            idx["RIG:1"] + idx["RIG:2"]
    return order, idx

def fetch_negf_matrices(log, order, tool):
    """
    Get the Hamiltonian and Overlap matrix in the right order for LibNEGF.

    Args:
        log (BigDFT.Logfiles.Logfile): the log from a calculation.
        order (list): the ordering to match LibNEGF (use `get_ordering`).
        tool (BigDFT.PostProcessing.BigDFTool): tool for post processing.
    """
    from scipy.sparse import coo_matrix

    # Read in all matrices
    H = tool.get_matrix_h(log)
    S = tool.get_matrix_s(log)
    K = tool.get_matrix_k(log)

    # Set the sparsity pattern to match
    H += 1e-20 * K
    H = 0.5 * H + 0.5 * H.T
    S += 1e-20 * K

    # Permute and convert to coo_matrix
    H = coo_matrix(permute_matrix(H, order))
    S = coo_matrix(permute_matrix(S, order))

    return H, S

def write_negf_mat(mat, ofile):
    """
    Write a matrix in the file format needed to be read by LibNEGF.

    Args:
        mat (scipy.sparse.coo_matrix): the matrix to write.
        ofile: an open file pointer for writing.
    """
    ofile.write(f" % Size = {mat.shape[0]} {mat.shape[1]}\n")
    ofile.write(f" % Nonzeros = {len(mat.data)}\n")
    ofile.write(f" zzz = zeros({len(mat.data)}, 3)\n")
    ofile.write(" zzz = [\n")
    for i, j, v in zip(mat.col, mat.row, mat.data):
        ofile.write(f"     {i + 1} {j + 1} {v}\n")
    ofile.write(" ]\n")

def write_negf_input(fermi_level, idx, ofile, kt=1.0e-3, emin=-0.1, emax=0.1,
                     estep=0.001):
    # index of the last of the middle
    plend = [len(idx["MID:0"])]
    # could also be cont start
    surfstart = [plend[0] + 1,
                 plend[0] + 1 + len(idx["LEF:1"]) + len(idx["LEF:2"])]
    # no surface
    surfend = [surfstart[0] - 1, surfstart[1] - 1]
    #last region
    contend = [surfstart[1] - 1, sum([len(x) for x in idx.values()])]

    # chemical potential
    mu = [fermi_level, fermi_level]
    # temperature
    kt = [kt, kt]
    cblk = [1, 1]

    # write the data
    ofile.write(f"{surfstart[0]} {surfstart[1]}\n")
    ofile.write(f"{surfend[0]} {surfend[1]}\n")
    ofile.write(f"{contend[0]} {contend[1]}\n")
    ofile.write(f"{plend[0]}\n")
    ofile.write(f"{mu[0]} {mu[1]}\n")
    ofile.write(f"{kt[0]} {kt[1]}\n")
    ofile.write(f"{cblk[0]} {cblk[1]}\n")
    ofile.write(f"{fermi_level + emin} {fermi_level + emax} {estep}\n")

def read_negf_density(fname, dim):
    """
    Read in the density matrix computed by LibNEGF.

    Args:
        fname (str): name of the density matrix file
        dim (int): dimension of the matrix
    """
    from scipy.sparse import dok_matrix
    K = dok_matrix((dim, dim))
    with open(fname) as ifile:
        for line in ifile:
            split = line.split()
            K[int(split[0]) - 1, int(split[1]) - 1] = float(split[2])
    return K

def update_density(bK, nK, order, idx):
    """
    Update the density from BigDFT using the correction of LibNEGF.

    Args:
        bK (scipy.sparse.coo_matrix): BigDFT density kernel
        nK (scipy.sparse.coo_matrix): the NEGF correction kernel.
        order (list): the ordering to match LibNEGF (use `get_ordering`).
        idx (dict): the fragment index ordering computed by BigDFTool.
    """
    # Diagonal middle region
    kput = permute_matrix(bK, order)
    srow = 0
    erow = len(idx["MID:0"])
    scol = 0
    ecol = len(idx["MID:0"])
    kput[srow:erow, scol:ecol] = nK[srow:erow, scol:ecol]

    # Off diagonal part 1
    srow = len(idx["MID:0"])
    erow = srow + len(idx["LEF:1"])
    scol = 0
    ecol = len(idx["LEF:1"])
    kput[srow:erow, scol:ecol] = nK[srow:erow, scol:ecol]
    kput[scol:ecol, srow:erow] = nK[scol:ecol, srow:erow]

    # Off diagonal part 2
    srow = len(idx["MID:0"]) + len(idx["LEF:1"]) + len(idx["LEF:2"])
    erow = srow + len(idx["RIG:1"])
    scol = len(idx["MID:0"]) - len(idx["LEF:1"])
    ecol = len(idx["MID:0"])
    kput[srow:erow, scol:ecol] = nK[srow:erow, scol:ecol]
    kput[scol:ecol, srow:erow] = nK[scol:ecol, srow:erow]

    return permute_matrix(kput, order, inv=True)

def run_libnegf(name, run_dir, skip=False, verbose=True):
    """
    Command to actually run LibNEGF

    Args:
        name (str): name of the job
        run_dir (str): where to execute this job.
        skip (str): if true, we do a lazy calculation.
    """
    from os import system, environ
    from os.path import join

    # Check skip
    if skip:
        try:
            with open(join(run_dir, f"{name}.log")) as ifile:
                for line in ifile:
                    if "Destroy negf" in line:
                        if verbose:
                            print("Skipping completed run")
                        return
        except IOError:
            pass

    # Run
    mpirun = environ.get("LIBNEGF_MPIRUN", None)
    path = environ.get("LIBNEGF_ROOT", "")
    run_line = f"cd {run_dir} ; "
    if mpirun is not None:
        run_line += mpirun + " "
    run_line += join(path, "drv_negf") + f" > {name}.log"
    if verbose:
        print(run_line)
    system(run_line)
