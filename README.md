# CurvedGCode-to-ABB
GCode_to_Robtargets v2.py converts Curved GCode (G02,G03) to MoveC instructions in ABB RAPID program

A simple way to 3D print concrete or other large structures is to repurpose a robot arm into a 3D printer. 
Ordinary 3D printer slicers can be used to create machine tool paths that extrude layers of the print. 
However many robot arms can not handle the sharp short movements generated by the linear G code commands especially those used for curved sections of print. 
The solution is to use the inbuilt curved movement functions of the robot arm to smooth out the motion. 
However curved G code uses a start point, destination point, radius point and whether the motion is anti-clockwise or clockwise to define curved movements. https://www.cnccookbook.com/cnc-g-code-arc-circle-g02-g03/
On the other hand, the MoveC instructions in ABB RAPID Programming uses a start point, destination point and mid point in the arc between the start and destination point. 
chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://library.e.abb.com/public/688894b98123f87bc1257cc50044e809/Technical%20reference%20manual_RAPID_3HAC16581-1_revJ_en.pdf

GCode_to_Robtargets v2.py converts Curved GCode (G02,G03) (.gcode file) to MoveC and MoveL instructions in ABB RAPID program. It was adapted from from Daniel Aguirre's  (https://github.com/DAguirreAg/GCode-to-ABB) linear Gcode translation to ABB RAPID programming instructions. 
test_robprogramSim.txt can be copied and pasted into a RAPID program. 
Alternatively, the contents of test_moveLsSim.txt and test_robtargetsSim can be pasted into a RAPID program. test_robtargetsSim being the position constants for the program test_moveLsSim.txt. 

The function midpoint_of_arc_calc() is the function that determines the midpoint of the arc from the curved GCode input.
The radius point comes in the form of points I and J which are the x and y distance between the start point to the radius point. 
The fuction first draws a line between start and destination and finds the centre point. Using this centre point and the radius point the midpoint will always be a radius length away from the radius point and fall on the line created by the centre point and the radius point. 
So the code finds the gradient of the line between the radius point and the centre point and the radius length. 
As 
'''y = mx and r^2 = x^2'''


Calc_Arc_Midpoint.py is used to test out the midpoint_of_arc_calc() function with a start point, end point and I, J points. 
An online plotting interface can be used to visualize the example such as Desmos https://www.desmos.com/calculator

