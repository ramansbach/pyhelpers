# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 12:05:37 2016

@author: Rachael Mansbach
Code to measure angular displacement wrt one of the planes 
(default should be XZ because of usually doing flow that way)
"""
import imp,numpy as np,os
#from quaternions import Quaternion
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
try:
    mk =imp.load_source('mk','/home/mansbac2/coarsegraining/code/markov/markov.py')
except IOError:
    mk =imp.load_source('mk','/home/rachael/cluster/home/mansbac2/coarsegraining/code/markov/markov.py')

#function that takes a vector and returns its angular displacement wrt a plane
#defined by normal vector n and then the principal vector p to measure the angular
#displacement to
#for example, the angle from the x axis in the XZ plane means passing
# n = yhat and p = xhat
#all of these should be numpy vectors
#also return vector projected in the plane
matplotlib.rcParams.update({'font.size': 22})
def angDisp(xi,n,p):
    #we are assuming this is for only a small change in angles, it WILL fail 
    #if we try angle wrapping for too large a timestep
    pxi = xi - np.dot(xi,n)*n
    cphi = np.dot(pxi,p)/np.linalg.norm(pxi)
    phi = np.arccos(cphi)

    if phi > np.pi/2.:
        phi = np.pi-phi
        pxi = -1*pxi
        
    dvec = np.cross(pxi,p) #direction vector -- which way are we rotating? important for checking diffusivity
    if min(dvec) < 0:
        phi = -1*phi
    return (phi,pxi/np.linalg.norm(pxi))

#at a particular time, get the angle of the cluster or peptide's principal
#component of gyration with the p axis in the plane normal to n    
def getPhi(t,xtc,tpr,outgro,n,p,ats,cutoff,rm=True):
    (peps,box_length) = mk.getPosB(t,xtc,tpr,outgro)
    #assume one cluster, but it MAY have broken over PBCs so
    peps = mk.fixPBC(peps,box_length,ats,cutoff)
    if rm:
        os.system('rm '+outgro)
    gyrationTensor = mk.gyrTens(peps,box_length)
    eigstuff = np.linalg.eig(gyrationTensor)
    eigOrder = np.argsort(eigstuff[0])
    xi1 = eigstuff[1][:,eigOrder[2]]
    (phi,pxi) = angDisp(xi1,n,p)    
       
    #Visualization Check
   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    peps3 = np.reshape(peps,[len(peps)/3,3])    
    ax.scatter(peps3[:,0],peps3[:,1],peps3[:,2])
    comapprox = (1./(6*29.))*np.sum(peps3,0)

    ax.quiver(comapprox[0]+xi1[0],comapprox[1]+xi1[1],comapprox[2]+xi1[2],xi1[0],xi1[1],xi1[2],color='m')
    ax.quiver(comapprox[0]+xi1[0],0,comapprox[2]+xi1[2],xi1[0],0,xi1[2],color='r')
    ax.quiver(comapprox[0]+p[0],0+p[1],comapprox[2]+p[2],p[0],p[1],p[2],color='b')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')    
    ax.set_zlabel('Z Label')
    plt.show()
    
    
    return (phi,pxi)
    
#same as getPhi, except takes the COM of the peptides as the locations of 
#the body to take the gyration tensor of
def getPhiCOM(t,xtc,tpr,outgro,n,p,ats,cutoff,masslist,rm=True):
    (peps,box_length) = mk.getPosB(t,xtc,tpr,outgro)
    #assume one cluster, but it MAY have broken over PBCs so
    peps = mk.fixPBC(peps,box_length,ats,cutoff)
    pepcoms = mk.getCOMs(peps,masslist,ats)
    if rm:
        os.system('rm '+outgro)
    gyrationTensor = mk.gyrTens(pepcoms,box_length)
    #print t
    #print gyrationTensor
    try:
        eigstuff = np.linalg.eig(gyrationTensor)
    except np.linalg.LinAlgError:
        temp =1
        print "this failed at",t
    eigOrder = np.argsort(eigstuff[0])
    xi1 = eigstuff[1][:,eigOrder[2]]
    (phi,pxi) = angDisp(xi1,n,p)    
    '''   
    #Visualization Check
   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    peps3 = np.reshape(peps,[len(peps)/3,3])    
    ax.scatter(peps3[:,0],peps3[:,1],peps3[:,2])
    comapprox = (1./(6*29.))*np.sum(peps3,0)

    ax.quiver(comapprox[0]+xi1[0],comapprox[1]+xi1[1],comapprox[2]+xi1[2],xi1[0],xi1[1],xi1[2],color='m')
    ax.quiver(comapprox[0]+xi1[0],0,comapprox[2]+xi1[2],xi1[0],0,xi1[2],color='r')
    ax.quiver(comapprox[0]+p[0],0+p[1],comapprox[2]+p[2],p[0],p[1],p[2],color='b')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')    
    ax.set_zlabel('Z Label')
    plt.show()
    '''
    
    return (phi,pxi)

def angMSD(dt,phis,dss):
    #compute <phi^2> for each ds (# timesteps) in dss and plot <phi^2> vs dts
    #along with errorbars from std-dev of mean
    mdphis2 = np.zeros(len(dss))
    sdphis2 = np.zeros(len(dss))    
    s = 0    
    for ds in dss:
        dphis2 = np.zeros(len(phis)-ds)
        k = 0
        for i in range(ds):
            dphi2 = np.diff(phis[i:len(phis):ds])**2

            dphis2[k:k+len(dphi2)] = dphi2
            
            k+=len(dphi2)
        mdphis2[s] = np.mean(dphis2)
        #sdphis2[s] = np.std(dphis2)/np.sqrt(len(dphis2))
        sdphis2[s] = np.std(dphis2)        
        s+=1
    fig = plt.figure()
    plt.errorbar(dt*dss,mdphis2,yerr=sdphis2)
    plt.plot(dt*dss,mdphis2,color='g',lw=3)
    plt.show()
    return(dss,mdphis2,sdphis2)
    
def getDr(dts,dphis,sdphis,savename):
    #get and return the slope of the linear region of the MSD graph
    linr = input('Enter the range over which to take the slope in the form a,b: ')
    
    b = float(linr[0])
    e = float(linr[1])
    bi = list(dts < b).index(False)
    be = list(dts > e).index(True)
    #lss = np.linalg.lstsq(np.transpose([dts[bi:be]]),np.transpose([dphis[bi:be]]))
    slope, intercept, r_value, p_value, std_err = stats.linregress(dts[bi:be],dphis[bi:be])    
    s = slope
    fig = plt.figure()
    plt.errorbar(dts,dphis,yerr=sdphis,marker='.')
    plt.plot(dts,dphis,color='g',lw=3)
    plt.plot(dts[bi:be],s*dts[bi:be]+intercept,color='r',lw=3,ls='--')
    plt.show()
    plt.savefig(savename)
    Dr = s/2.
    return Dr
if __name__ == "__main__":
    # test angMsd
    #test = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.])
    #test = np.sqrt(test)
    #(dss,mdphis2,sdphis2) = angMSD(1,test,range(1,9))
    #test this but basically do the following:
    #(1) find the phi angle of the first timestep with the x vector in the xz plane
    #(2) then find the change in angle between the v(t0) and v(t0+dt)
    #(3) finally, write out phi_new = phi(t)+dphi
    
    folder = '/home/rachael/coarsegraining/CG/DFAG/parameterization/BIbonded/340ns_params/dih_enf_full_tab_NP_water/5_md/'
    t0 = 0
    dt = 0.64*10
    nsteps = 100
    printmod = 10
    ats = 29
    cutoff = 0.5
    tf = dt*nsteps#dt*5
    tlist = np.arange(t0+dt,tf+dt,dt)
    xtc = folder+'md_whole.xtc'
    tpr = folder+'md_dummy.tpr'
    #define plane, right now it is the XZ plane
    n = np.array([0.0,1.0,0.0])
    p = np.array([1.0,0.0,0.0])
    outfolder = '/home/rachael/coarsegraining/flow/gromacs/rotPec/size1/'
    outfname = outfolder+'phivt.dat'
    outf = open(outfname,'w')
    phis = np.zeros(1+len(tlist))
    (phi,pxi) = getPhi(t0,xtc,tpr,folder+'temp.gro',n,p,ats,cutoff,True)
    phis[0] = phi
    #print phi
    pind = 1
    outf.write('{0}\t{1}\t{2}\n'.format(t0,phi,phi*(180.0/np.pi)))
    for t in tlist:
        if not(pind % printmod):
            print pind
        (dphi,pxi) = getPhi(t,xtc,tpr,folder+'temp.gro',n,pxi,ats,cutoff,True)
        #print dphi
        phi += dphi
        phis[pind] = phi
        pind+=1
        outf.write('{0}\t{1}\t{2}\n'.format(t,phi,phi*(180.0/np.pi)))
    outf.close()
    (dss,mdphis2,sdphis2) = angMSD(dt,phis,np.arange(1,nsteps))
    Dr = getDr(dt*dss,mdphis2,sdphis2,outfolder+'Drfit.png')
    print "Dr is %e 1/ps = %e 1/s" % (Dr,Dr*10**12)