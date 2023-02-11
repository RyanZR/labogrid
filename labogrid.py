#!/usr/bin/env python3

#==============================================================================
#     ___      _______  _______  _______  _______  ______    ___   ______  
#    |   |    |   _   ||  _    ||       ||       ||    _ |  |   | |      | 
#    |   |    |  |_|  || |_|   ||   _   ||    ___||   | ||  |   | |  _    |
#    |   |    |       ||       ||  | |  ||   | __ |   |_||_ |   | | | |   |
#    |   |___ |       ||  _   | |  |_|  ||   ||  ||    __  ||   | | |_|   |
#    |       ||   _   || |_|   ||       ||   |_| ||   |  | ||   | |       |
#    |_______||__| |__||_______||_______||_______||___|  |_||___| |______| 
#                                                  
#   LABOGRID - Gridbox size calculation for ligand docking
#
#   This software is provided WITHOUT WARRANTY OF ANY KIND
#
#   Get more information from https://github.com/RyanZR/labogrid
#
#   Report bugs and issues to https://github.com/RyanZR/labogrid/issues
#
#==============================================================================

import os
import sys
import getopt
import statistics

def usage():
    print(f"Usage")
    print(f"╰─○ labogrid.py -i <ligand_filename>")
    print(f"")
    print(f"Argument")
    print(f"╰─○ Description of commands:")
    print(f"         -i   ligand filename (supported: pdb, pdbqt, sdf, mol2)")
    print(f"         -h   help")
    print(f"         -a   about")
    print(f"╰─○ Optional parameters:")
    print(f"        [-s]  scale factor (default is 2)")

def about():
    print(f"==============================================================================")
    print(f"     ___      _______  _______  _______  _______  ______    ___   ______      ")
    print(f"    |   |    |   _   ||  _    ||       ||       ||    _ |  |   | |      |     ")
    print(f"    |   |    |  |_|  || |_|   ||   _   ||    ___||   | ||  |   | |  _    |    ")
    print(f"    |   |    |       ||       ||  | |  ||   | __ |   |_||_ |   | | | |   |    ")
    print(f"    |   |___ |       ||  _   | |  |_|  ||   ||  ||    __  ||   | | |_|   |    ")
    print(f"    |       ||   _   || |_|   ||       ||   |_| ||   |  | ||   | |       |    ")
    print(f"    |_______||__| |__||_______||_______||_______||___|  |_||___| |______|     ")
    print(f"                                                                              ")
    print(f"   LABOGRID - Gridbox size calculation for ligand docking                     ")
    print(f"                                                                              ")
    print(f"   This software is provided WITHOUT WARRANTY OF ANY KIND                     ")
    print(f"                                                                              ")
    print(f"   Get more information from https://github.com/RyanZR/labogrid               ")
    print(f"                                                                              ")
    print(f"   Report bugs and issues to https://github.com/RyanZR/labogrid/issues        ")
    print(f"                                                                              ")
    print(f"==============================================================================")

def coordinate_XYZ(data, ext):
    if ext in ".mol2":
        start = [n for n, line in enumerate(data) if line.strip() != "" if line.split()[0] in "@<TRIPOS>ATOM"][0] + 1
        end = [n for n, line in enumerate(data) if line.strip() != "" if line.split()[0] in "@<TRIPOS>BOND"][0]
        xcoor = [float(line.split()[2]) for line in data[start:end]]
        ycoor = [float(line.split()[3]) for line in data[start:end]]
        zcoor = [float(line.split()[4]) for line in data[start:end]]

    if ext in ".sdf":
        blankIdx = int([n for n, line in enumerate(data) if line.split() == []][0])
        numAtom = int(data[blankIdx+1].split()[0])
        start = blankIdx+2
        end = blankIdx+2+numAtom
        xcoor = [float(line.split()[0]) for line in data[start:end]]
        ycoor = [float(line.split()[1]) for line in data[start:end]]
        zcoor = [float(line.split()[2]) for line in data[start:end]]

    if ext in (".pdb", ".pdbqt"):
        xcoor = [float(line.split()[6]) for line in data if line.split()[0] in ("ATOM", "HETATM")]
        ycoor = [float(line.split()[7]) for line in data if line.split()[0] in ("ATOM", "HETATM")]
        zcoor = [float(line.split()[8]) for line in data if line.split()[0] in ("ATOM", "HETATM")]

    return [xcoor, ycoor, zcoor]

def min_max(coor: list):
    return [min(coor), max(coor)]

def mid_XYZ(ranCoor: list):
    return round(statistics.mean(ranCoor), 3)

def length_WHD(ranCoor: list, scl: float):
    return round(abs(ranCoor[0] - ranCoor[1])*scl, 3)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:s:ha")
    except getopt.GetoptError as msg:
        print(f"labogrid.py")
        print(f"╰─○ {msg}")
        sys.exit()

    receptor_filename = None
    xcoor = None
    ycoor = None
    zcoor = None
    scale = 2
    scale = float(scale)

    for opt, arg in opts:
        if opt in "-i":
            receptor_filename = arg
        if opt in "-s":
            scale = arg
        if opt in "-h":
            usage()
            sys.exit()
        if opt in "-a":
            about()
            sys.exit()

    if receptor_filename == None:
        print(f"labogrid.py")
        print(f"╰─○ Invalid file or incorrect usage")
        sys.exit()

    EXT = os.path.splitext(receptor_filename)[-1]
    if EXT not in (".pdb", ".pdbqt", ".sdf", ".mol2"):
        print(f"labogrid.py")
        print(f"╰─○ File format {EXT} not supported")
        sys.exit()     

    DATA = open(receptor_filename, "r").readlines()
    COOR = coordinate_XYZ(DATA, EXT)
    X, Y, Z = COOR[0], COOR[1], COOR[2]
    ranges = [min_max(X), min_max(Y), min_max(Z)]
    center = [mid_XYZ(ranges[0]), mid_XYZ(ranges[1]), mid_XYZ(ranges[2])] 
    bxsize = [length_WHD(ranges[0], scale), length_WHD(ranges[1], scale), length_WHD(ranges[2], scale)]

    print(f"labodock.py")
    print(f"╰─○ Ligand Center:  X {center[0]}  Y {center[1]}  Z {center[2]}")
    print(f"╰─○ Gridbox Size :  W {bxsize[0]}  H {bxsize[1]}  D {bxsize[2]}")

if __name__ == "__main__":
    main()