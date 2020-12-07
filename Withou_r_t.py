import math
import numpy as np
import sys
import getopt
from scipy.spatial import distance

# This function read params files and return the c6,c12 for energy function.
def read_params_file(fileIn):
    fileIn_3 = open(fileIn, "r")
    # the first index of the atom
    atom1 = []
    # the second index of the atom
    atom2 = []
    # c6 values of all pairs of atoms
    c6 = []
    # c12 values of all pairs of atoms
    c12 = []
    # read each line from the file and add it to new array.
    lines_f = []
    # read all lines and append it in line_f.
    for line in fileIn_3:
        lines_f.append(line)
    # read each line and extract useful information for energy function.
    for i in range(len(lines_f)):
        text1 = lines_f[i].split()
        atom1.append(text1[0])
        atom2.append(text1[1])
        c6.append(text1[2])
        c12.append(text1[3])
    return atom1,atom2,c6,c12

# This function can read pdb or pdbm files and return different values you needed.
def read_PDB_file(fileIn):
    # array to store x value
    x = []
    # array is to store y value
    y = []
    z = []
    # to add atom name in this array
    name = []
    # store all lines from a file
    line1 = []
    # store all lines that didn't contain the "H"
    line_H = []
    # store the charge of each atom
    charge = []
    print("read file not include H atom?")
    c=str(input())
    # read the file
    for line in fileIn:
        line1.append(line)
    # extract the value from file, such as atom name, and x,y,z
    for i in range(len(line1)):
        text1 = line1[i].split()
        for j in range(len(text1)):
            if text1[0].startswith("ATOM"):
                # to check if you only want to use x,y,and z values
                if c == "Y" or c=="y":
                    if not text1[2].startswith("H") and text1[2] != "CZ":
                        line_H.append(line1[i])
                        name.append(text1[2])
                        x.append(text1[6])
                        y.append(text1[7])
                        z.append(text1[8])
                        charge.append(text1[9])
                        break
                elif c == "N" or c=="n":
                    line_H.append(line1[i])
                    name.append(text1[2])
                    x.append(text1[6])
                    y.append(text1[7])
                    z.append(text1[8])
                    charge.append(text1[9])
                    break
                else:
                    print("Please input Y or N")
                    break
    # this step is for returning values to function for future use
    print("choose xyz or atoms to return values.")
    while (True):
        choose = str(input())
        if choose == "xyz":
            return x, y, z
        elif choose == "atoms":
            return name, x, y, z, charge
        else:
            print("Please choose xyz or atoms.")


# This function is to find the best angle of the ligand with the lowest energy.
# and return the best angle.
def ligand_coordinates(x1,y1,z1,degree):
    list1=[]
    new_x=[]
    new_y=[]
    new_z=[]
    # create matrix for rotation
    b_Euler = np.array([[math.cos(degree), 0.0, math.sin(degree)], [0.0, 1.0, 0.0], [-math.sin(degree), 0.0, math.cos(degree)]])
    xyz=np.array([[x1,y1,z1]])
    new_array= np.dot(xyz,b_Euler)
    list1=new_array.tolist()
    new_x.append(list1[0][0])
    new_y.append(list1[0][1])
    new_z.append(list1[0][2])
    return(new_x,new_y,new_z,degree)

# The energe function is to calculate the energy between each ligand and protein
def energy_function(file_L, file_p, dict):
    try:
        # read ligand file and get the values
        atom1, x1, y1, z1, charge1 = read_PDB_file(file_L)
        # read protein file and get the values
        atom2, x2, y2, z2, charge2 = read_PDB_file(file_p)
        # sum of the energy between ligand and protein
        v_psum = 0.0
        v_tsum=0.0
        # energy between i atom and j atom
        v = 0.0
        new_x = []
        new_y = []
        new_z = []
        atom_L=atom1
        min_v=[]
        v_n=[]
        # initial angle for ligand and protein
        # max degree is 2pi
        max_degree = 2 * (math.pi)
        d = 0.0
        count=0
        for i in range(len(x1)): ## len(x1)==20
            new_x, new_y, new_z,degree = ligand_coordinates(float(x1[i]), float(y1[i]), float(z1[i]),math.degrees(d))
            v_psum=0.0
            print(new_x,new_y,new_z)
            for j in range(len(x2)):  # len(x2)== 2736
                value = []
                if (atom1[i], atom2[j]) in dict:
                    value = dict[atom1[i], atom2[j]]
                elif (atom2[j], atom1[i]) in dict:
                    value = dict[atom2[j], atom1[i]]
                rij = math.sqrt((float(new_x[0]) - float(x2[j])) ** 2 + (float(new_y[0]) - float(y2[j])) ** 2 +
                                (float(new_z[0]) - float(z2[j])) ** 2)
                rij = rij / 10
                c6 = float(value[0])
                c12 = float(value[1])
                vij = c12 / rij - c6 / rij
                print(atom1[i], atom2[j])
                print("VlJ")
                print(vij)
                if rij < 1.32:
                    dielectric_c = 8
                else:
                    # 1/4piE0
                    dielectric_c = 138.935485
                qi = float(charge1[i])
                qj = float(charge2[j])
                A = 6.02944
                B = 72.37056
                gama = 0.01873345
                k = 213.5782
                e_rij = A + B / (1 + k * math.exp(-gama * B * rij))
                ve = dielectric_c * ((qi * qj) / (e_rij * rij))
                v = vij + ve
                count+=1
                v_psum+=v
            d=(d + math.pi / 6)
            break

    except TypeError:
        print("run again!")
    except UnboundLocalError:
        print("Please use atoms as input")


# RMSD function is to calculate the root mean square deviation between two structures.
def RMSD(x1,y1,z1,x2,y2,z2):
    sum_d2=0
    n=len(x1)
    for i in range(n):
        d2=(float(x1[i])-float(x2[i]))**2+(float(y1[i])-float(y2[i]))**2+(float(z1[i])-float(z2[i]))**2
        sum_d2 +=d2
    if n>0:
        rmsd =math.sqrt(sum_d2/n)
    else:
        rmsd=None
    return rmsd

# using a function usage() should display usage information to the screen
# if an option at the command-line is not recognized (use exception handling).
def usage():
    print("Usage: ", sys.argv[0], "[-p FILE] [-l FILE]")
    print(" -p: input file ; STDIN if not used.")
    print(" -l: input file ; STDIN if not used.")


# This function is for choosing the problem that you want to solve.
# The first one is to calculate the RMSD between the two structures.
# The second one is docking the each ligand with protein and find the best angle with lowest energy.
def problem_choose(fileIn_1,fileIn_2):
    print("Do you want to solve RMSD between the ligands or implementing protein-ligand docking algorithm?")
    while(True):
        method = str(input())
        if method == "RMSD":
            x1,y1,z1=read_PDB_file(fileIn_1)
            x2,y2,z2=read_PDB_file(fileIn_2)
            if len(x1)!=len(x2):
                print("Please check the file")
            else:
                rmsd = RMSD(x1, y1, z1, x2, y2, z2)
                n=len(x1)
                print("\nRMSD =" + str(round(rmsd, 3)) + " A for " + str(n) + " atom.")
        elif method =="Docking":
            atom1, atom2, c6, c12 = read_params_file("ffG43b1nb.params")
            dict = {}
            for i in range(len(atom1)):
                dict[atom1[i], atom2[i]] = c6[i], c12[i]
            energy_function(fileIn_1, fileIn_2, dict)
            print("End End")
            break
        else:
            print("Please input RMSD or Docking to choose method.")


def main():
    # standard input, using pdb file as input file
    fileIn_1 = sys.stdin
    # standard input, using pdb file as input file
    fileIn_2 = sys.stdin

    # create the options and the argument for each option,
    # and the exception of the option.
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:l:')
        # if did not give right option, exception.
    except getopt.GetoptError as err:
        sys.stdout = sys.stderr
        print(str(err))
        usage()
        sys.exit(2)
    for (opt, arg) in opts:
        if (opt == '-p'):
            fileIn_1 = open(arg, "r")
        elif (opt == '-l'):
            fileIn_2 = open(arg, "r")
    atom1, atom2, c6, c12 = read_params_file("ffG43b1nb.params")
    dict = {}
    for i in range(len(atom1)):
        dict[atom1[i], atom2[i]] = c6[i], c12[i]
    energy_function(fileIn_1,fileIn_2, dict)
    #problem_choose(fileIn_1,fileIn_2)

    # colse fileIn
    fileIn_1.close()
    fileIn_2.close()
if (__name__ == "__main__"):
    main()

