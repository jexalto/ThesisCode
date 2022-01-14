import numpy as np
import scipy.interpolate as si

from SlipstreamConway import SlipstreamConway


class SlipstreamTubeModel(object):
    """A slipstreamtube model used to calculate the induced velocities by the 
    propeller slipstream. This slipstream model is able to handle slipstream
    contraction, slipstream deflection and a radial and azimuthal distribution
    of circulation on the propeller plane. The caclulations for this 
    slipstreamtube model are taken from the PhD thesis by Veldhuis
    
    L.L.M. Veldhuis (2005) Propeller Wing Aerodynamic Interference 
    """
    
    def __init__(self):
        
        self.B = 6          #number of blades
        self.Vinf = 30      #freestream velocity [m/s]
        self.omega = 0.5    #rotational velocity [rad/s]
        self.R = 4          #propeller radius [m]
        self.r_hub = 0.5    #hub radius [m]
        
        self.sign_rotation = -1 #rotation direction, 1=clockwise, -1=counter clockwise
        
        #slipstream vorticity
        self.gp = None      #propeller plane vorticity
        self.ga = None      #axial vorticity
        self.gt = None      #tangential vorticity
        
        #slipstream geometry
        self.phi = None     #azimuth positions
        self.r = None       #radial positions
        self.r_r0 = None    #slipstream contraction
        self.z = None       #z location
        self.x = None       #x location
        
        #simplified Conway model
        self.cnw = SlipstreamConway()
        
    def slipstream_geom(self, phi, r, x, z, gamma, contraction=True):
        """Initializes the slipstream geometry. Calculates the locations in the
        slipstream where the velocity is known. Contraction can be added using
        a simplified slipstreamtube model.
        
        Args:
            phi (numpy.array): 
        """
        
        #x and z are in the same dimension
        self.x = x.reshape((x.size,1,1))
        self.z = z.reshape(self.x.shape)
        
        #phi
        dphi = phi[1]-phi[0]
        phi = phi.reshape((phi.size,1))
        #rotation based on Vinf and omega
        phi_rot = self.x/self.Vinf*self.omega*-1*self.sign_rotation
        self.phi = phi+phi_rot
        
        #calculate radial positions where circulation is shed in slipstream
        r_ = (r[1:]+r[:-1])*0.5
        r_ = np.append(self.r_hub, r_)
        r_ = np.append(r_, self.R)
        self.r = r_
        
        #dr
        dr = np.diff(r)
        dr = np.append(r[0]-self.r_hub, dr)
        dr = np.append(dr, self.R-r[-1])
        
        #difference in gamma between radial stations is shed in slipstream
        g = gamma[:, 1:]-gamma[:, :-1]
        g = np.hstack((gamma[:, 0:1], g))
        g = np.hstack((g, -gamma[:, -1:]))
        
        #use same control points in slipstream and on propeller plane
        gb = np.zeros(g.shape)
        for i in range(gb.shape[0]):
            gb[i] = np.interp(r_, r, gamma[i])
        
        n = self.omega/(2*np.pi) #rotational speed [Hz]
        
        #add contraction based on conservation of mass
        if contraction:
            
            #calculate azimuthal average per radial station
            self.gp = np.average(gb, axis=0)*self.B/(2*np.pi)*dphi*dr/(4*np.pi)
            self.ga = np.average(g, axis=0)/dr*self.B/(2*np.pi*r_)*r_*dphi*dr/(4*np.pi)*-1
            self.gt = np.average(g, axis=0)/dr*n*self.B/self.Vinf*r_*dphi*dr/(4*np.pi)*-1
            
            #determine axial induced velocities at propeller plane
            Vax0 = np.zeros(r.shape)
            for i, ri in enumerate(r):
                Vt, Vw = self.induced_velocity(np.array([0, ri, 0]))
                Vax0[i] = Vw[0]
            
            #create simplified slipstreamtube model
            cnw = self.cnw
            cnw.calculate_Vz0(Vax0, r/self.R)
            
            #calculate axial velocity for all x and r/R combinations
            Vax = cnw.v_axial(x, r_/self.R, self.R)
            self.cnw = cnw
            
            #get the axial velocities at the propeller plane
            Vax0_ = Vax[:, 0]
            Vax0_ = Vax0_.reshape((Vax0_.size, 1))
            
            #calculate contraction
            r_r0 = ((self.Vinf+Vax0_)/(self.Vinf+Vax))**0.5
            r_r0 = np.transpose(r_r0)
            r_r0 = r_r0.reshape((x.size, 1, r_.size))
            
            #apply contraction
            self.r_r0 = r_r0
            self.r = self.r*r_r0
        
        #calculate circulation values
        self.gp = gb*self.B/(2*np.pi*r_)*dphi*dr/(4*np.pi)          #propeller
        self.ga = g/dr*self.B/(2*np.pi*r_)*r_*dphi*dr/(4*np.pi)*-1  #axial
        self.gt = g/dr*n*self.B/self.Vinf*r_*dphi*dr/(4*np.pi)*-1   #tangential
      
    def induced_velocity(self, P):
        """Uses an interpolation scheme to find the induced velocity. In this
        discreet system the induced velocity cannot be calculated at any random
        point, if the point is too close to the points where circulation is
        defined, inaccurate results may be obtained. Thus this function finds
        the closest 'control points' to the given point P and finds the value
        of velocity at P by interpolation using the control points.
        
        Args:
            P (numpy.array): array with x,y,z coordinates
        
        Returns:
            Vtotal (numpy.array): velocity vector for the total induced velocity
            Vwake (numpy.array): velocity vector for the induced velocity by
                the wake only
        """
        
        #convert to polar coordinates
        xp = P[0]
        zp = np.interp(xp, self.x.flatten(), self.z.flatten())
        rp = (P[1]**2+(P[2]-zp)**2)**0.5
        pp = np.arctan2(P[1], -P[2]+zp)
        if pp<0:
            pp += 2*np.pi
        
        #check if length of sliptube is sufficient
        if xp>0.51*np.max(self.x):
            print(xp, 0.51*np.max(self.x))
            print('Slipstream geometry is too short, increase x range')
        
        #find x control points
        x_x = self.x.flatten()
        x_i = np.argsort(np.abs(x_x-xp))[0]
        
        x_i1 = x_i+1
        x_i2 = x_i-1
        
        x1 = (x_x[x_i]+x_x[x_i1])/2
        if x_i2==-1:
            x_i2 = 0
            x2 = (2*x_x[x_i]-(x_x[x_i1]-x_x[x_i]))/2
        else:
            x2 = (x_x[x_i]+x_x[x_i2])/2
        
        #check if point lies within sliptube
        out = False
        if self.r.ndim==1:
            if rp>np.max(self.r) or rp<np.min(self.r):
                out = True
        else:
            r_x = self.r[x_i].flatten()
            if rp>np.max(r_x) or rp<np.min(r_x):
                out = True
        if xp<np.min(self.x):
            out = True
        
        #interpolation
        if out:
            Vtotal, Vwake = self.induced_velocity_calculation(P)
        else:
            if self.r.ndim==1:
                #find control radius values
                r_x = self.r
                r_i = np.argsort(np.abs(r_x-rp))[0]
                r_i1 = r_i+1
                r_i2 = r_i-1
                
                if r_i1==r_x.size:
                    r1 = (3*r_x[r_i]-r_x[r_i2])/2
                    r2 = (r_x[r_i]+r_x[r_i2])/2
                elif r_i2==-1:
                    r1 = (r_x[r_i]+r_x[r_i1])/2
                    r2 = (3*r_x[r_i]-r_x[r_i1])/2
                else:
                    r1 = (r_x[r_i]+r_x[r_i1])/2
                    r2 = (r_x[r_i]+r_x[r_i2])/2
                
                r3 = r1
                r4 = r2
            else:
                #find control radius values
                r_x = self.r[x_i].flatten()
                r_x1 = self.r[x_i1].flatten()
                r_x2 = self.r[x_i2].flatten()
                
                r_i = np.argsort(np.abs(r_x-rp))[0]
                r_i_1 = np.argsort(np.abs(r_x1-rp))[0] 
                r_i_2 = np.argsort(np.abs(r_x2-rp))[0]
                
                r_i1 = r_i+1
                r_i2 = r_i-1
                
                r_xi = r_x[r_i]
                if r_i1==r_x.size:
                    r_xi1 = 2*r_x[r_i]-r_x[r_i2]
                    r_xi2 = r_x[r_i2]
                elif r_i1==-1:
                    r_xi1 = r_x[r_i1]
                    r_xi2 = 2*r_x[r_i]-r_x[r_i1]
                else:
                    r_xi1 = r_x[r_i1]
                    r_xi2 = r_x[r_i2]
                        
                r_i_11 = r_i_1+1
                r_i_12 = r_i_1-1
                
                r_xi_1 = r_x1[r_i_1]
                if r_i_11==r_x.size:
                    r_xi_11 = 2*r_x1[r_i_1]-r_x1[r_i_12]
                    r_xi_12 = r_x1[r_i_12]
                elif r_i_11==-1:
                    r_xi_11 = r_x1[r_i_11]
                    r_xi_12 = 2*r_x1[r_i_1]-r_x1[r_i_11]
                else:
                    r_xi_11 = r_x1[r_i_11]
                    r_xi_12 = r_x1[r_i_12]
                
                r_i_21 = r_i_2+1
                r_i_22 = r_i_2-1
                
                r_xi_2 = r_x2[r_i_2]
                if r_i_21==r_x.size:
                    r_xi_21 = 2*r_x2[r_i_2]-r_x2[r_i_22]
                    r_xi_22 = r_x2[r_i_22]
                elif r_i_21==-1:
                    r_xi_21 = r_x2[r_i_21]
                    r_xi_22 = 2*r_x2[r_i_2]-r_x2[r_i_21]
                else:
                    r_xi_21 = r_x2[r_i_21]
                    r_xi_22 = r_x2[r_i_22]
                
                r1 = (r_xi+r_xi1+r_xi_1+r_xi_11)/4
                r2 = (r_xi+r_xi2+r_xi_1+r_xi_12)/4
                r3 = (r_xi+r_xi1+r_xi_2+r_xi_21)/4
                r4 = (r_xi+r_xi2+r_xi_2+r_xi_22)/4
            
            #find the control phi values
            p_x = self.phi[x_i].flatten() % (2*np.pi)
            p_x1 = self.phi[x_i1].flatten() % (2*np.pi)
            p_x2 = self.phi[x_i2].flatten() % (2*np.pi)
            
            p_i = np.argsort(np.abs(p_x-pp))[0]
            p_i_1 = np.argsort(np.abs(p_x1-pp))[0] 
            p_i_2 = np.argsort(np.abs(p_x2-pp))[0]
            
            if p_i==0:
                p_i1, p_i, p_i2 = 1, 0, -1
            elif p_i==p_x.size-1:
                p_i1, p_i, p_i2 = 0, -1, -2
            else:
                p_i1 = p_i+1
                p_i2 = p_i-1
                
            if p_i_1==0:
                p_i_11, p_i_1, p_i_12 = 1, 0, -1
            elif p_i_1==p_x.size-1:
                p_i_11, p_i_1, p_i_12 = 0, -1, -2
            else:
                p_i_11 = p_i_1+1
                p_i_12 = p_i_1-1
                
            if p_i_2==0:
                p_i_21, p_i_2, p_i_22 = 1, 0, -1
            elif p_i_2==p_x.size-1:
                p_i_21, p_i_2, p_i_22 = 0, -1, -2
            else:
                p_i_21 = p_i_2+1
                p_i_22 = p_i_2-1  
            
            p1 = (p_x[p_i]+p_x[p_i1]+p_x1[p_i_1]+p_x1[p_i_11])/4
            p2 = (p_x[p_i]+p_x[p_i2]+p_x1[p_i_1]+p_x1[p_i_12])/4
            p3 = (p_x[p_i]+p_x[p_i1]+p_x2[p_i_2]+p_x2[p_i_21])/4
            p4 = (p_x[p_i]+p_x[p_i2]+p_x2[p_i_2]+p_x2[p_i_22])/4
            
            #create rotation matrices for found angles
            rm1 = np.array([[1,0,0],
                            [0, np.cos(p1), np.sin(p1)],
                            [0,-np.sin(p1), np.cos(p1)]])
            rm2 = np.array([[1,0,0],
                            [0, np.cos(p2), np.sin(p2)],
                            [0,-np.sin(p2), np.cos(p2)]])
            rm3 = np.array([[1,0,0],
                            [0, np.cos(p3), np.sin(p3)],
                            [0,-np.sin(p3), np.cos(p3)]])
            rm4 = np.array([[1,0,0],
                            [0, np.cos(p4), np.sin(p4)],
                            [0,-np.sin(p4), np.cos(p4)]])
    
            #control points in cartesian coordinates
            P1 = np.array([x1, r1*np.sin(p1), -r1*np.cos(p1)+zp])
            P2 = np.array([x1, r1*np.sin(p2), -r1*np.cos(p2)+zp])
            P3 = np.array([x1, r2*np.sin(p1), -r2*np.cos(p1)+zp])
            P4 = np.array([x1, r2*np.sin(p2), -r2*np.cos(p2)+zp])
            P5 = np.array([x2, r3*np.sin(p3), -r3*np.cos(p3)+zp])
            P6 = np.array([x2, r3*np.sin(p4), -r3*np.cos(p4)+zp])
            P7 = np.array([x2, r4*np.sin(p3), -r4*np.cos(p3)+zp])
            P8 = np.array([x2, r4*np.sin(p4), -r4*np.cos(p4)+zp])
            
            #control points in polar coordinates
            P1_ = np.array([x1, r1, p1])
            P2_ = np.array([x1, r1, p2])
            P3_ = np.array([x1, r2, p1])
            P4_ = np.array([x1, r2, p2])
            P5_ = np.array([x2, r3, p3])
            P6_ = np.array([x2, r3, p4])
            P7_ = np.array([x2, r4, p3])
            P8_ = np.array([x2, r4, p4])
            
            #all points collected
            Pt = np.array([P1, P2, P3, P4, P5, P6, P7, P8])
            Pt_ = np.array([P1_, P2_, P3_, P4_, P5_, P6_, P7_, P8_])
            rm = np.array([rm1, rm2, rm1, rm2, rm3, rm4, rm3, rm4])
            Vw = np.zeros(Pt.shape)
            Vt = np.zeros(Pt.shape)
            
            #calculate velocities at control points
            for i in range(Pt.shape[0]):
                vt_, vw_ = self.induced_velocity_calculation(Pt[i])
                vt_ = np.matmul(rm[i], vt_)
                vw_ = np.matmul(rm[i], vw_)
                Vw[i] = vw_
                Vt[i] = vt_
            
            #interpolate
            Vwake = np.zeros(3)
            Vtotal = np.zeros(3)
            P_ = np.array([xp, rp, pp])
            for i in range(3):
                Vwake[i] = si.griddata(Pt_, Vw[:, i], P_)
                Vtotal[i] = si.griddata(Pt_, Vt[:, i], P_)
            
            #transform to cartesian velocity vector
            rmp = np.array([[1,0,0],
                            [0, np.cos(pp), np.sin(pp)],
                            [0,-np.sin(pp), np.cos(pp)]])
            Vwake = np.matmul(rmp.T, Vwake)
            Vtotal = np.matmul(rmp.T, Vtotal)
            
            #if interpolation failed, calculate velocities
            if (True in np.isnan(Vwake)):
                Vtotal, Vwake = self.induced_velocity_calculation(P)

        return Vtotal, Vwake
                    
    def induced_velocity_calculation(self, P):
        """Calculates the induced velocity using the slipstreamtube model
        
        Args:
            P (numpy.array): array with x,y,z coordinates
        
        Returns:
            Vtotal (numpy.array): velocity vector for the total induced velocity
            Vwake (numpy.array): velocity vector for the induced velocity by
                the wake only
        """
        
        #initialize
        phi = self.phi
        r = self.r
        z = self.z
        x1 = self.x
        x2 = np.zeros(x1.shape)
        x2[:-1] = self.x[1:]
        x2[-1] = self.x[-1]+(self.x[-1]-self.x[-2])
        
        gp = self.gp
        ga = self.ga
        gt = self.gt
        
        #frequently used geometric relations
        b = r*np.sin(phi)-P[1]
        c = -r*np.cos(phi)-P[2]+z
        b2c2 = b**2+c**2
        d = (b*np.sin(phi)-c*np.cos(phi))/(b2c2)
        r1 = ((x1-P[0])**2+b2c2)**0.5
        r2 = ((x2-P[0])**2+b2c2)**0.5
        
        #induced velocity by the propeller circulation
        up = gp*(P[1]*np.cos(phi[0])+P[2]*np.sin(phi[0]))/(P[0]**2+b2c2[0])**1.5
        vp = gp*(-P[0]*np.cos(phi[0]))/(P[0]**2+b2c2[0])**1.5
        wp = gp*(-P[0]*np.sin(phi[0]))/(P[0]**2+b2c2[0])**1.5
        
        #induced velocity by the tangential circulation
        ut = gt*d*((x2-P[0])/r2-(x1-P[0])/r1)
        vt = gt*(np.sin(phi)/r2-np.sin(phi)/r1)
        wt = gt*(-np.cos(phi)/r2+np.cos(phi)/r1)
        
        #induced velocity by the axial circulation
        va = ga*c/b2c2*((x2-P[0])/r2-(x1-P[0])/r1)
        wa = ga*(-b)/b2c2*((x2-P[0])/r2-(x1-P[0])/r1)
        
        #sum contribution of all circulation
        Vp = np.array([np.sum(up), np.sum(vp), np.sum(wp)])
        Va = np.array([0, np.sum(va), np.sum(wa)])
        Vt = np.array([np.sum(ut), np.sum(vt), np.sum(wt)])
    
        Vwake = Va+Vt
        Vtotal = Vp+Va+Vt
        
        return Vtotal, Vwake
