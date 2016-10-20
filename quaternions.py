# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:56:21 2016
Simple code for a quaternion class
@author: rachael
"""
import numpy as np
#a quaternion has four components
#If viewed as a rotation of theta about the (a,b,c) axis, the four components are
# w = cos(theta/2), x = a*sin(theta/2), y = b*sin(theta/2), z = c*sin(theta/2)
# this class stores the components and implements quaternion multiplication
class Quaternion(object):
    def __init__(self,w,x,y,z):
        self.vals = np.array([w,x,y,z])
    #define left and right multiplication between quaternions
    def __mul__(self,quat2):
        q1 = self.vals
        q2 = quat2.vals
        q3 = np.zeros(4)
        q3[0] = q1[0]*q2[0]-np.dot(q1[1:3],q2[1:3])
        q3[1] = 
        return q3