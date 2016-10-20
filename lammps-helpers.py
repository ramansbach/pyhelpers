# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:32:24 2016

@author: rachael
In this file, we have some simple setup scripts to initialize LAMMPS data files
"""
import numpy as np,imp

class MDSystem(object):
    #a class that contains the important things about an MD system for ease of conversion
    #including #atoms, velocities, positions, and box lengths
    def __init__(self,natoms=0,natypes=1):
        self.natoms = natoms
        
        self.xs = np.zeros([self.natoms,3])
        self.vs = np.zeros([self.natoms,3])
        self.box = np.zeros(3)
        self.nbonds = 0
        self.nangs = 0
        self.ndihs = 0
        self.nimps = 0
        self.nparams = [self.natoms,self.nbonds,self.nangs,self.ndihs,self.nimps]
        self.ntypes = [natypes,0,0,0,0]
        self.masses = np.zeros(self.ntypes[0])
        self.qs = np.zeros(self.natoms)
        self.namparam = ['atom','bond','angle','dihedral','improper']
        self.types = np.ones(self.natoms)        
        
    def set_natoms(self,natoms):
        self.natoms = natoms
        self.nparams[0] = natoms
        self.xs = np.zeros([self.natoms,3])
        self.vs = np.zeros([self.natoms,3])
        
    def set_box(self,box1,box2,box3):
        self.box[0] = box1
        self.box[1] = box2
        self.box[2] = box3
    
    def set_vs(self,vs):
        self.vs = vs
    
    def set_xs(self,xs):
        self.xs = xs
    
    def set_v(self,v,i):
        self.vs[i,:] = v
        
    def set_x(self,x,i):
        self.xs[i,:] = x
        
    def set_masses(self,ms):
        self.masses = ms
        
    def set_qs(self,qs):
        self.qs = qs
    def get_xs(self):
        return self.xs.copy()
        
    def get_vs(self):
        return self.vs.copy()
        
    def get_box(self):
        return self.box.copy()
        

class LammpsSystem(MDSystem):
    #class that inherits from MDsystem but writes out in a particular format
    #for now writing out ellipses, may redo for non-elliptical (point) particles
    def write(self,outname):
        outf = open(outname,'w')
        outf.write("File written by lammps-helpers.py\n")
        outf.write("\n")
        for i in range(len(self.namparam)):
            outf.write("{0} {1}s\n".format(self.nparams[i],self.namparam[i]))
        outf.write("{0} ellipsoids\n".format(self.natoms))
        outf.write("\n")
        for i in range(len(self.namparam)):
            if self.ntypes[i] > 0:
                outf.write("{0} {1} types\n".format(self.ntypes[i],self.namparam[i]))
        
        outf.write("\n")
        outf.write("{0} {1} xlo xhi\n".format(0.,self.box[0]))
        outf.write("{0} {1} ylo yhi\n".format(0.,self.box[1]))
        if self.box[2] > 0:
            outf.write("{0} {1} zlo zhi\n".format(0.,self.box[2]))
        #outf.write("\n")
        #outf.write("Masses\n")
        #outf.write("\n")
        #for i in range(self.ntypes[0]):
        #    outf.write("{0} {1}\n".format(i+1,self.masses[i]))
        outf.write("\n")
        outf.write("Atoms\n")
        outf.write("\n")
        for i in range(self.natoms):
            outf.write(" {0} {1} {2} {3} {4} {5} {6}\n".format(i+1,int(self.types[i]),1,1.0,self.xs[i,0],self.xs[i,1],self.xs[i,2]))
        outf.write("\n")
        outf.write("Velocities\n")
        outf.write("\n")
        for i in range(self.natoms):
            outf.write(" {0} {1} {2} {3} {4} {5} {6}\n".format(i+1,self.vs[i,0],self.vs[i,1],self.vs[i,2],0.0,0.0,0.0))
        outf.write("\n")
        outf.write("Ellipsoids\n")
        outf.write("\n")
        for i in range(self.natoms):
            outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i+1,0.1,0.1,0.1,0.0,0.0,0.0,0.0))
                    
            
        outf.write("\n")
        outf.close()
    
    

def ljfromgro(groname):
    #this function converts from a gro file containing a single atom type to a 
    #lammps data file of LJ type

    grof = open(groname,'r')
    lines = grof.readlines()
    grof.close()
    sys = LammpsSystem(natoms=int(lines[1]))
    
    boxL1 = float(lines[len(lines)-1].split()[0])
    boxL2 = float(lines[len(lines)-1].split()[1])
    boxL3 = float(lines[len(lines)-1].split()[2])
    sys.set_box(boxL1,boxL2,boxL3)
    i = 0
    for line in lines[2:len(lines)-1]:
        spline = line.split()
        xs = [float(j) for j in spline[3:6]]
        vs = [float(j) for j in spline[6:9]]
        sys.set_x(xs,i)
        sys.set_v(vs,i)
        i+=1
    return sys
    
    
def insertHollow(watersys,spar,sper,buff,fracBox):
    #takes a water box and inserts a hollow cylinder where an ellipse would be with major axis sper, minor spar
    #currently assumed to be oriented with the z axis at box location (0,0,0)
    #and removes overlapping water molecules as well as putting ellipse in
    #for simplicity of coding, let's assume a cylindrical empty area
    #then loop through water molecule positions and remove them if 
    #(i) their z position is between comz-sper-buff & comz+sper+buff
    #(ii) their (x,y) distance from (comx,comy) <= spar+buff
    #buffer is probably the size of the lj radius of the water molecules to be safe
    #fracBox says what fraction along each of the box length vectors to put the center of mass
    box = watersys.get_box()
    xs = watersys.get_xs()
    vs = watersys.get_vs()
    com = fracBox*box
    newxs = []
    newvs = []
    for i in range(np.shape(xs)[0]):
        x = xs[i][0]
        y = xs[i][1]
        z = xs[i][2]
        r = np.sqrt((x-com[0])**2+(y-com[1])**2)
        if not ((z > com[2]-spar-buff) and (z < com[2]+spar+buff) and (r <= sper + buff)):
            newxs.append(xs[i,:])
            newvs.append(vs[i,:])
    watersys.set_natoms(np.shape(newxs)[0])
    watersys.set_xs(np.array(newxs))
    watersys.set_vs(np.array(newvs))


if __name__ == "__main__":
    runfolder = "/home/rachael/coarsegraining/flow/lammps/martiniW/big-ellipse/"
    groname = "water.gro"
    outname = "martW.lj"
    outname2 = "martWhollow.lj"
    sys = ljfromgro(runfolder+groname)
    sys.write(runfolder+outname)
    insertHollow(sys,2.5,0.7143,0.3,np.array([0.5,0.5,0.5]))
    sys.write(runfolder+outname2)