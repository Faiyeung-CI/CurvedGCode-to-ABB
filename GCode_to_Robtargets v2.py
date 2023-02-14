
##############################
# Code to translate curved and linear GCode into moveL and moveC functions for ABB robot program#
# Created by:
# Adapted code from Daniel Aguirre https://github.com/DAguirreAg/GCode-to-ABB #
# Date: 2022/09/10          #
#                            #
##############################

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys, getopt
import os
import math
from math import *

# GLOBAL VARIABLES
inputfile = "G:\Shared drives\AMAT SIF Projects\AMAT SIF 3D Concrete Printing J03784 CMP\GCode to RAPID\GCode to RAPID -Curvyinfill Vase - ArcWelder\PI3MK3M_D210 d110 CurvyInfill Cylinder v2 10pers scale (first layer).gcode"
outputfile_robtargets = "test_robtargetsSim.txt"
outputfile_moveLs = "test_moveLsSim.txt"
outputfile_robprogram = "test_robprogramSim.txt"

# rotation = "[0.70746,0.017029,0.70654,-0.00362]" #for real life
# conf = "[-1,-1,0,0]" #for real life

rotation = "[0.000772984,-0.000000011,-0.999999701,-0.000000009]" #for simulation
conf = "[-1,-1,0,0]" #for simulation

# rotation = "[0.000772984,-0.000000011,-0.999999701,-0.000000009]" #new simulation

G90 = True

# Program variables
positions = pd.DataFrame([[0.0,0.0,0,0,0,0,0,0,0,0]], columns=["x","y","z","e","cw","i","j","mx","my","mz"])
robtargets = []
robprogram = []
moveLs = []
NumArrays = []
layer_h = 0
scalingFactor = 10

def midpoint_of_arc_calc(px,py,pz,x,y,z,cw,i,j):

    # clockwise = cw
    #
    # A = [px,py] #start point
    # B = [x,y] #destination point
    # O = [i,j] #I and J values
    # # print(A)
    # # print(B)
    # # print(O)
    #
    # C = [A[0] + O[0],A[1] + O[1]]
    # r = sqrt(O[0]**2 + O[1]**2)
    #
    # theta2 = atan ((C[1] - A[1])/(C[0] - A[0]))
    # theta1 = atan ((C[1] - B[1])/(C[0] - B[0]))
    #
    # #New
    # # if theta1 > 0 and theta2 > 0 or theta2 < 0 and theta2 < 0:
    # #     if theta2 > theta1:
    # #         theta  = theta2 - theta1
    # #         a = theta/2 + theta1
    # #     else:
    # #         theta = theta1 - theta2
    # #         a = theta/2 + theta2
    # # else:
    # #     theta = abs(theta1) + abs(theta2)
    # #     if theta1 < 0:
    # #         a = -theta/2 + theta1
    # #     if theta2 < 0:
    # #         a = -theta/2 + theta2
    # #New
    #
    # #Old
    # theta  = theta2 - theta1
    #
    # a = theta/2 + theta1
    # #Old
    #
    # m = tan(a)
    # h = C[0]
    # k = C[1]
    #
    # coe1 = m**2+1
    # coe2 = -2*h - 2*h*m**2
    # coe3 = h**2 + h**2*m**2 - r**2
    # roots = np.roots([coe1, coe2, coe3])
    #
    # x1 = roots[0]
    # x2 = roots[1]
    #
    # y1 = m*(x1-h)+k
    # y2 = m*(x2-h)+k
    #
    # cor1 = [x1,y1]
    # cor2 = [x2,y2]
    #
    # bear1 = atan(abs(y1 - A[1])/abs(x1 - A[0]))
    # bear2 = atan(abs(y2 - A[1])/abs(x2 - A[0]))
    # bear3 = atan(abs(C[1] - A[1])/abs(C[0] - A[0]))
    #
    # if clockwise == 2:
    #     if (bear3 - bear1) > 0:
    #         return(cor1)
    #     else:
    #         return(cor2)
    #
    # if clockwise == 1:
    #     if (bear3 - bear1) > 0:
    #         return(cor2)
    #     else:
    #         return(cor1)

    clockwise = cw

    A = [px,py] #start point
    B = [x,y] #destination point
    O = [i,j] #I and J values
    # print(A)
    # print(B)
    # print(O)

    #Draw line between start and destination, find centre point
    C = [(A[0] + B[0]) / 2, (A[1] + B[1]) / 2]
    print (C)

    #get coordinate of radius
    R = [A[0] + O[0], A[1] + O[1]]
    print (R)

    #get gradient of line
    m = (R[1]-C[1])/(R[0]-C[0])
    print (m)

    #r length of radius
    r = sqrt(O[0]**2 + O[1]**2)
    print (r)

    #we know that the midpoint will be r length from the radius point
    #and we know that the y = mx
    #so r**2 = x**2 + y**2
    #so r**2 = x**2 + (mx)**2
    #solving for component x gives x = sqrt(r**2/(1+m**2))

    Cx = sqrt(r**2/(1+m**2))
    print(Cx)

    Cy = m*Cx
    print(Cy)

    #the 2 possible midpoints will either be radius point +x and +y from the point
    #or -x and -y from the point

    x1 = R[0] + Cx
    x2 = R[0] - Cx

    y1 = R[1] + Cy
    y2 = R[1] - Cy

    #The 2 possible points are a and b which are also vectors
    a = [x1,y1]
    b = [x2,y2]

    print(a)
    print(b)

    #but the orgin of vectors need to be at the radius point
    #refStart is the vector from the radius point to the starting point A
    refa = [x1-R[0], y1-R[1]]
    refb = [x2-R[0], y2-R[1]]
    refStart = [A[0]-R[0], A[1]-R[1]]

    #compare the determinant of the combined vectors of refstart and refa and the combined vecotrs of refstart and refb
    #if the deta > 0 then point a is anticlockwise from the start point
    #if the deta < 0 then point a is clockwise from the start point
    #https://math.stackexchange.com/questions/1027476/calculating-clockwise-anti-clockwise-angles-from-a-point

    deta = refStart[0]*refa[1] - refStart[1]*refa[0]
    print(deta)

    detb = refStart[0]*refb[1] - refStart[1] *refb[0]
    print(detb)

    if clockwise == 1:
        if  deta < 0:
            print("Final result is" + str(a))
            return(a)
        else:
            print("Final result is" + str(b))
            return(b)

    if clockwise == 2:
        if  deta < 0:
            print("Final result is" + str(b))
            return(b)
        else:
            print("Final result is" + str(a))
            return(a)



# Extract GCode command specific information
def parseCommand(command, G90):
    global positions
    global layer_h
    temp = command.split()

    if len(temp) == 0:
        temp = [";"]

    # if temp[0]=="M106":
    #     x = positions.iloc[-1].x
    #     y = positions.iloc[-1].y
    #     z = 0.0
    #
    #     newPosition = pd.DataFrame([[x,y,z]], columns=["x","y","z"])
    #     positions = pd.concat([positions, newPosition])
    #
    # elif temp[0]=="M107":
    #     x = positions.iloc[-1].x
    #     y = positions.iloc[-1].y
    #     z = -10.0
    #
    #     newPosition = pd.DataFrame([[x,y,z]], columns=["x","y","z"])
    #     positions = pd.concat([positions, newPosition])

    if temp[0]=="G00" or temp[0]=="G0" or temp[0]=="G1" or temp[0]=="G01" or temp[0]=="G92":

        x = 0.0
        y = 0.0
        z = 0.0
        e = 0.0
        dx = 0.0
        dy = 0.0
        dz = 0.0
        flag = 0
        cw = 0
        i = 0
        j = 0
        mx = 0
        my = 0
        mz = 0

        for comp in temp:
            if comp.startswith("X"):
                 dx = float(comp[1:])
            if comp.startswith("Y"):
                dy = float(comp[1:])
            if comp.startswith("Z"):
                dz = float(comp[1:])
                layer_h = dz
            if comp.startswith("F") and float(comp[1:]) == 2100:
                flag = 1
                # dx =  positions.iloc[-1].x #comment out if GCode doesn't give Z commands separately
                # dy = positions.iloc[-1].y #comment out if GCode doesn't give Z commands separately
            # if comp.startswith("E"):
            #     #dx = float(comp[1:])
            #     if str(comp[1:]) == "0":
            #         #Extrusion off
            #         print("Extrusion off")
            #         e = 1
            #     if float(comp[1:]) == 2.0:
            #         #Extrusion on
            #         print("Extrusion on")
            #         e = 2
            #     if float(comp[1:]) > 2.0 and len(temp) < 4:
            #         flag = 1
            #     # else:
            #     #     e = 0


        if flag != 1:
            if G90==True: #Use absolute coordinates
                x = dx
                y = dy
                z = layer_h
            else:
                x = positions.iloc[-1].x + dx
                y = positions.iloc[-1].y + dy
                z = positions.iloc[-1].z + layer_h

            newPosition = pd.DataFrame([[x,y,z,e,cw,i,j,mx,my,mz]], columns=["x","y","z","e","cw","i","j","mx","my","mz"])
            positions = pd.concat([positions, newPosition])
            print(positions)
            flag = 0

    if temp[0] == "G3" or temp[0] == "G2":

        x = 0.0
        y = 0.0
        z = positions.iloc[-1].z
        e = 0.0
        dx = 0.0
        dy = 0.0
        dz = 0.0
        flag = 0
        cw = 0
        i = 0
        j = 0
        mx = 0
        my = 0
        mz = 0

        if temp[0] == "G2":
            cw = 1
        if temp[0] == "G3":
            cw = 2
        px = positions.iloc[-1].x
        py = positions.iloc[-1].y
        pz = positions.iloc[-1].z
        mz = pz
        for comp in temp:
            if comp.startswith("X"):
                 x = float(comp[1:])
            if comp.startswith("Y"):
                y = float(comp[1:])
            # if comp.startswith("Z"):
            #     z = float(comp[1:])
            #     layer_h = z
            if comp.startswith("I"):
                i = float(comp[1:])
            if comp.startswith("J"):
                j = float(comp[1:])

        mx = midpoint_of_arc_calc(px,py,pz,x,y,z,cw,i,j)[0]
        my = midpoint_of_arc_calc(px,py,pz,x,y,z,cw,i,j)[1]

        newPosition = pd.DataFrame([[x,y,z,e,cw,i,j,mx,my,mz]], columns=["x","y","z","e","cw","i","j","mx","my","mz"])
        positions = pd.concat([positions, newPosition])
        print(positions)




# Writes the Robtarget points into a file
def writeRobtarget(i, position):
    if position.e == 0:
        x = position.x
        y = position.y
        z = position.z
        string1 = "CONST robtarget "
        string2 = "p" + str(i)
        string3 = ":= [[" + str(x) + "," + str(y) + "," + str(z) + "], " + rotation + ", " + conf + ", [ 9E+09,9E+09,9E+09,9E+09,9E+09,9E+09]];"
        string4 = "\n"
        robtarget = string1 + string2 + string3 + string4
        robtargets.append(robtarget)

def writeRobprogram(i, position):
    global scalingFactor
    if position.e == 0 and position.cw == 0:
        x = position.x * scalingFactor
        y = position.y * scalingFactor
        z = position.z * scalingFactor
        string1 = "MoveL "
        string2 = "[[" + str(x) + "," + str(y) + "," + str(z) + "], " + rotation + ", " + conf + ", [ 9E+09,9E+09,9E+09,9E+09,9E+09,9E+09]]"
        # string3 = ",speed, zone, NOZZLE_1 \WObj:=Platform;"
        string3 = ",speed, zone, tExtruder \WObj:=Platform;" #sim
        string4 = "\n"
        robcommand = string1 + string2 + string3 + string4
        robprogram.append(robcommand)

    if position.cw == 1 or position.cw == 2 and position.e == 0:
        x = position.x * scalingFactor
        y = position.y * scalingFactor
        z = position.z * scalingFactor
        mx = position.mx * scalingFactor
        my = position.my * scalingFactor
        mz = position.mz * scalingFactor
        string1 = "MoveC "
        string2 = "[[" + str(mx) + "," + str(my) + "," + str(mz) + "], " + rotation + ", " + conf + ", [ 9E+09,9E+09,9E+09,9E+09,9E+09,9E+09]],[[" + str(x) + "," + str(y) + "," + str(z) + "], " + rotation + ", " + conf + ", [ 9E+09,9E+09,9E+09,9E+09,9E+09,9E+09]]"
        # string3 = ",speed, zone, NOZZLE_1 \WObj:=Platform;"
        string3 = ",speed, zone, tExtruder \WObj:=Platform;" #sim
        string4 = "\n"
        robcommand = string1 + string2 + string3 + string4
        robprogram.append(robcommand)


# Writes the GCode points into a file
# def writeNumArray(i, position):
#     x = position.x
#     y = position.y
#     z = position.z
#     string3 = "[" + str(x) + "," + str(y) + "," + str(z) + "], "
#     string4 = "\n"
#     NumArray = string3 + string4
#     NumArrays.append(NumArray)



# Writes the ABB move commands into a file
def writeMoveL(i, position):
    if position.e == 0:
        string1 = "moveL "
        string2 = "p" + str(i)
        string3 = ",speed, zone, NOZZLE_1 \WObj:=Platform;" #for real life
        # string3 = ",v200, z1, tExtruder \WObj:=Platform;" #for simulation
        #string3 = ",v200, fine, NOZZLE_1;"
        string4 = "\n"
        moveL = string1 + string2 + string3 + string4
        moveLs.append(moveL)
    elif position.e == 1:
        # moveLs.append("!Extrusion off\n") #dummy
        moveLs.append("Set do_Output_2;\n")
    elif position.e == 2:
        # moveLs.append("!Extrusion on\n") #dummy
        moveLs.append("Reset do_Output_2;\nSet do_Output_3;\n")

def plotPath(projection="3d"):
    x = np.array(positions.x, dtype=pd.Series)
    y = np.array(positions.y, dtype=pd.Series)
    z = np.array(positions.z, dtype=pd.Series)

    x0 = np.array([1,])
    y0 = np.array([1,])

    fig = plt.figure()

    if (projection == "3d"):
        ax = Axes3D(fig)
        ax.plot(x,y,z)

    else:
        ax = fig.add_subplot(111)
        ax.plot(x,y,"red")
        ax.scatter(x,y)

    plt.show()


################## MAIN ##################
def main(argv):

    global inputfile, outputfile_robtargets, outputfile_moveLs, rotation, conf, outputfile_robprogram

    try:
         opts, args = getopt.getopt(argv,"hi:o:r:c:",["help","ifile=","ofile="])
    except getopt.GetoptError:
         print('Error')
         sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Usage: GCode_to_Robtargets [-h | -i <inputfile> -o <outputfile>] ")
            print('Options and arguments:')
            print("-h     : Print this help message and exit")
            print("-i arg : Input the file to be converted into ABB instructions (also --ifile)")
            print("-o arg : Output filename containing the ABB instructions (also --ofile)")
            print("-r arg : Specify the rotation of the robtargets (also --rot). Default: [-1, 0, 0, 0]")
            print("-c arg : Specify the axis configuration of the robtargets (also --conf). Default: [-1, 0, 1, 0]")
            sys.exit()

        elif opt in ("-i", "--ifile"):
            inputfile = arg

        elif opt in ("-o", "--ofile"):
            outputfile_robtargets = arg + "_robtargets.txt"
            outputfile_moveLs = arg + "_moveLs.txt"
            outputfile_robprogram = arg + "_robprogram.txt"

        elif opt in ("-r", "--rot"):
            rotation = arg

        elif opt in ("-c", "--conf"):
            conf = arg


    # Check if Input file has been defined
    if inputfile == None:
         print("Inputfile not defined")
         sys.exit(2)

    # Load GCode and obtain XYZ coordinates
    print(os.getcwd())
    file = open(inputfile,"r")
    with open(inputfile,"r") as file:
        line = file.readline()
        lineNumberCount = 1
        while line:
            print("Line: " + str(lineNumberCount))
            line = file.readline()
            parseCommand(line, G90)
            lineNumberCount += 1

    # Write Robtargets and MoveL to a txt file
    for i in range(0, positions.shape[0]-1):
        position = positions.iloc[i]
        writeRobtarget(i, position)
        writeMoveL(i, position)
        writeRobprogram(i, position)

    with open(outputfile_robtargets,"w") as file:
        for line in robtargets:
            file.writelines(line)

    with open(outputfile_robprogram,"w") as file:
        for line in robprogram:
            file.writelines(line)

    with open(outputfile_moveLs,"w") as file:
        for line in moveLs:
            file.writelines(line)

    print("Conversion finished")

    # Plot expected result
    plotPath()


if __name__ == "__main__":
   main(sys.argv[1:])
