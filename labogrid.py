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
#   Get more information from https://github.com/RyanZR/labogrid
#
#   Report bugs and issues to https://github.com/RyanZR/labogrid/issues
#
#   This software is provided WITHOUT WARRANTY OF ANY KIND
#
#==============================================================================

import os
import sys
import getopt
import statistics

def usage():
    print(f"Usage")
    print(f"╰─○ labogrid.py [-l] <filename>")
    print(f"")
    print(f"Argument")
    print(f"╰─○ Command description:")
    print(f"        -l  ligand filename (supported: pdb, pdbqt, sdf, mol2)")
    print(f"        -e  experimental ligand filename (supported: pdb, sdf, mol2)")
    print(f"        -h  help")
    print(f"        -a  about")
    print(f"╰─○ Optional parameters:")
    print(f"        -s  scale factor (default is 2)")

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
    print(f"   Get more information from https://github.com/RyanZR/labogrid               ")
    print(f"                                                                              ")
    print(f"   Report bugs and issues to https://github.com/RyanZR/labogrid/issues        ")
    print(f"                                                                              ")
    print(f"   This software is provided WITHOUT WARRANTY OF ANY KIND                     ")
    print(f"                                                                              ")
    print(f"==============================================================================")

def file_handler(inputFile, inputType):
    if inputFile == None:
        print(f"labogrid.py")
        print(f"╰─○ Invalid file or incorrect usage")
        sys.exit()
    if not os.path.exists(inputFile):
        print(f"labogrid.py")
        print(f"╰─○ File does not exists")
        print(f"    Is {inputFile} mispelled?")
        sys.exit()
    EXT = os.path.splitext(inputFile)[-1]
    if inputType == "L" and EXT not in (".pdb", ".pdbqt", ".sdf", ".mol2"):
        print(f"labogrid.py")
        print(f"╰─○ File format {EXT} not supported")
        sys.exit()
    if inputType == "E" and EXT not in (".pdb", ".sdf", ".mol2"):
        print(f"labogrid.py")
        print(f"╰─○ File format {EXT} not supported")
        sys.exit()

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
        xcoor = [float(line[31:38]) for line in data if line.split()[0] in ("ATOM", "HETATM")]
        ycoor = [float(line[39:46]) for line in data if line.split()[0] in ("ATOM", "HETATM")]
        zcoor = [float(line[47:54]) for line in data if line.split()[0] in ("ATOM", "HETATM")]

    return [xcoor, ycoor, zcoor]

def min_max(coor: list):
    return [min(coor), max(coor)]

def mid_XYZ(rngCoor: list):
    return round(statistics.mean(rngCoor), 3)

def length_WHD(rngCoor: list, scl: float):
    return round(abs(rngCoor[0] - rngCoor[1])*scl, 3)

def labogrid(data: str, ext: str, inputType: str, scale: float):
    COOR = coordinate_XYZ(data, ext)
    X,Y,Z = COOR[0], COOR[1], COOR[2]
    ranges = [min_max(X), min_max(Y), min_max(Z)]
    center = [mid_XYZ(ranges[0]), mid_XYZ(ranges[1]), mid_XYZ(ranges[2])] 
    bxsize = [length_WHD(ranges[0], scale), length_WHD(ranges[1], scale), length_WHD(ranges[2], scale)]
    print(f"labodock.py")
    if inputType in "L":
        print(f"╰─○ Gridbox Size :  W {bxsize[0]}  H {bxsize[1]}  D {bxsize[2]}")
    if inputType == "E":
        print(f"╰─○ Ligand Center:  X {center[0]}  Y {center[1]}  Z {center[2]}")
        print(f"    Gridbox Size :  W {bxsize[0]}  H {bxsize[1]}  D {bxsize[2]}")

def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], 
            ":l:e:s:ha", 
            ["ligand=", "experimental=", "scale=", "help=", "about="])
    except getopt.GetoptError as msg:
        print(f"labogrid.py")
        print(f"╰─○ {msg}")
        sys.exit()
    
    LIG_file = None
    EXP_file = None
    SCALE = 2
    SCALE = float(SCALE)
    xcoor = None
    ycoor = None
    zcoor = None

    if opts == []:
        usage()
        sys.exit()
        
    for opt, arg in opts:
        if opt in ("-l", "--ligand"):
            LIG_file = arg
            INPUT_TYPE = "L"
        if opt in ("-e", "--experimental"):
            EXP_file = arg
            INPUT_TYPE = "E"
        if opt in ("-s", "--scale"):
            SCALE = arg
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        if opt in ("-a", "--about"):
            about()
            sys.exit()

    TARGET = LIG_file if INPUT_TYPE == "L" else EXP_file
    file_handler(TARGET, INPUT_TYPE)
    EXT = os.path.splitext(TARGET)[-1]
    DATA = open(TARGET, "r").readlines()
    labogrid(DATA, EXT, INPUT_TYPE, SCALE)
    
if __name__ == "__main__":
    main()