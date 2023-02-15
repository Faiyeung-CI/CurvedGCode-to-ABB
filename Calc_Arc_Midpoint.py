

##############################
#                            #
# Created by: Faiyeung Szeto #
# Date: 2022/09/10           #
#                            #
##############################

# Imports
import pandas as pd
import numpy as np
import sys, getopt
import os
import math
from math import *



################## MAIN ##################
def main(argv):

    clockwise = False

# G0 F3600 X14.406 Y21.347 Z0.4
# G3 X12.903 Y22.531 I5.597 J8.651 E5.22351 F1800

    # A = [14.406,21.347] #start point
    # B = [12.903,22.531] #destination point
    # O = [5.597, 8.651] #I and J values

    A = [14.66,22.42] #start point
    B = [16.411,22.647] #destination point
    O = [0.762, 0.988] #I and J values

    # A = [14.406,21.347] #start point
    # B = [12.903,22.531] #destination point
    # O = [5.597, 8.651] #I and J values

    # A = [13.939,24.05] #start point
    # B = [15.213,22.976] #destination point
    # O = [5.852, 5.649] #I and J values

    # A = [14.610,22.420] #start point
    # B = [16.411,22.647] #destination point
    # O = [0.762, 0.988] #I and J values


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

    #compare the determinant of the combined vectors of refstart and refa and the combined vectors of refstart and refb
    #if the deta > 0 then point a is anticlockwise from the start point
    #if the deta < 0 then point a is clockwise from the start point
    #https://math.stackexchange.com/questions/1027476/calculating-clockwise-anti-clockwise-angles-from-a-point

    deta = refStart[0]*refa[1] - refStart[1]*refa[0]
    print(deta)

    detb = refStart[0]*refb[1] - refStart[1] *refb[0]
    print(detb)

    if clockwise == True:
        if  deta < 0:
            print("Final result is" + str(a))
        else:
            print("Final result is" + str(b))

    if clockwise == False:
        if  deta < 0:
            print("Final result is" + str(b))
        else:
            print("Final result is" + str(a))



if __name__ == "__main__":
   main(sys.argv[1:])
