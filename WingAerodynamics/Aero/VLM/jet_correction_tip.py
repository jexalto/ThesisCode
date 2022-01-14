import numpy as np
import scipy.special as sp


def jet_correction_tip(y, mu, crd, R, dlbda=0.1, dlb=0.005, p_max=10, lbda_max=5):
    """Calculates the jet correction proposed by Rethorst, based on a velocity
    ratio mu. For this function the original correction has been applied to 
    the tip of the wing, where the center of the jet coincides with the wingtip.
    
    S.C. Rethorst (1958) Aerodynamic of Nonuniform Flows as Related to an
    Airfoil Extending Through a Circular Jet
    
    Args:
        y (numpy.array): spanwise wing locations
        mu (float): velocity ratio V0/Vj
        crd (float): wing (mean) chord
        R (float): jet radius
        dlbda (float): integral discretization step
        dlb (float): integral discretization step
        p_max (int): upper limit sum
        lbda_max (float): upper limit integral
    
    Returns:
        G (numpy.array): correction matrix
    """
    #calculate geometry
    b = y[-1]*2
    s0 = (y[-1]-y[-2])/2
    
    r0 = b/2-y[np.argmin(np.abs(y-b/2+R))]

    y = np.sort(-y)
    y_ = (y[1:]+y[:-1])/2
    
    #calculations are made in a system where the jet is on y=0
    y2 = y_+b/2-s0
    y2_ = y+b/2-s0
    r0_ = r0-s0
    #control point y coordinates
    eta = y2/r0_
    eta = eta.reshape(eta.size, 1)
    #horseshoe vortex midpoitn y coordinates
    beta = y2/r0_
    #horseshoe vortex points
    c = y2_[:-1]/r0_
    d = y2_[1:]/r0_
    #non-dimensional control point x [-]
    dzeta = -crd*1/2/r0_
    
    #Even vortex jet correction

    #calculate correction
    Gjje = 1/r0_*(1-mu**2)/(1+mu**2)*(1/(1/d-eta)-1/(1/c-eta)+1/(1/d+eta)-1/(1/c+eta))
    Goje = -1/r0_*(1-mu)**2/(1+mu**2)*(1/(eta-c)-1/(eta-d)+1/(eta+d)-1/(eta+c))
    Gjoe = Goje
    Gooe = -Gjje
    
    #determine if control point and horseshoe vortex are in jet
    beta_in = np.where(beta<1, 1, 0)
    eta_in = np.where(eta<1, 0.5, 0)
    G_in = beta_in+eta_in
    
    #apply right correction
    Ge = np.where(G_in==0, Gooe, 0)
    Ge = np.where(G_in==1, Goje, Ge)
    Ge = np.where(G_in==0.5, Gjoe, Ge)
    Ge = np.where(G_in==1.5, Gjje, Ge)
    Ge[:,0] *= 0.5 #At beta=0 there's only one hs vortex
    
    #clear memory
    Gjje, Goje, Gjoe, Gooe = None, None, None, None
    
    #Odd vortex jet correction
    
    #change location of the horseshoe vortex on the symmetry plane
    c[0] = 0
    #to avoid singularities
    eta = np.where(np.abs(eta)<1e-15, 1e-6, eta)
    
    #create integration variables p, lbda, lb
    p = np.arange(0, 10, 1)
    p = p.reshape(p.size, 1, 1)
    
    dlbda = 0.1
    lbda = np.arange(0, 5, dlbda)
    lbda = lbda.reshape(lbda.size, 1, 1, 1)
    
    dlb = 0.005
    lb = np.arange(dlb, 1, dlb)
    lb = lb.reshape(lb.size, 1, 1, 1, 1)
    lb = c*lbda+lb*(d-c)*lbda
    dlb = (d-c)*dlb*lbda
    
    #calculate integrals
    int_I = sp.iv(2*p+1, lb)/lb*dlb
    int_I = np.sum(int_I, axis=0)
    int_I = np.where(np.isnan(int_I), 0, int_I)
    
    int_K = sp.kv(2*p+1, lb)/lb*dlb
    int_K = np.sum(int_K, axis=0)
    int_K = np.where(np.isnan(int_K), 0, int_K)
    lb = None
    
    #calculate frequently used terms
    Kp = sp.kvp(2*p+1, lbda)
    K = sp.kv(2*p+1, lbda)
    K_ = sp.kv(2*p+1, lbda*eta)
    Ip = sp.ivp(2*p+1, lbda)
    I = sp.iv(2*p+1, lbda)
    I_ = sp.iv(2*p+1, lbda*eta)
    sin = np.sin(dzeta*lbda)
    mu1 = (1/((1/mu**2)-1))
    mu2 = ((1/mu)-mu)
    
    #calculate 4 corrections
    int_jj = K*Kp*I_*sin/(mu1-lbda*I*Kp)*int_I*dlbda
    int_jj = np.where(np.isnan(int_jj), 0, int_jj)
    int_jj = np.sum(int_jj, axis=0)
    sum_jj = (2*p+1)**2*int_jj
    sum_jj = np.sum(sum_jj, axis=0)
    Gjjo = 8/(r0_*np.pi*eta)*sum_jj
    sum_jj = None
    
    int_oj = (1/(mu-lbda*mu2*I*Kp)-1)*K_*sin/lbda*int_I*dlbda
    int_oj = np.where(np.isnan(int_oj), 0, int_oj)
    int_oj = np.sum(int_oj, axis=0)
    sum_oj = (2*p+1)**2*int_oj
    sum_oj = np.sum(sum_oj, axis=0)
    Gojo = 8/(r0_*np.pi*eta)*sum_oj
    sum_oj = None
    
    int_jo = (1/(mu-lbda*mu2*I*Kp)-1)*I_*sin/lbda*int_K*dlbda
    int_jo = np.where(np.isnan(int_jo), 0, int_jo)
    int_jo = np.sum(int_jo, axis=0)
    sum_jo = (2*p+1)**2*int_jo
    sum_jo = np.sum(sum_jo, axis=0)
    Gjoo = 8/(r0_*np.pi*eta)*sum_jo
    sum_jo = None
    
    int_oo = I*Ip*K_*sin/(mu1-lbda*I*Kp)*int_K*dlbda
    int_oo = np.where(np.isnan(int_oo), 0, int_oo)
    int_oo = np.sum(int_oo, axis=0)
    sum_oo = (2*p+1)**2*int_oo
    sum_oo = np.sum(sum_oo, axis=0)
    Gooo = 8/(r0_*np.pi*eta)*sum_oo
    sum_oo = None
    
    #apply right correction
    Go = np.where(G_in==0, Gooo, 0)
    Go = np.where(G_in==1, Gojo, Go)
    Go = np.where(G_in==0.5, Gjoo, Go)
    Go = np.where(G_in==1.5, Gjjo, Go)
    
    #apply the correction on one side of the wing, the influence of the virtual
    #vortices on the other side of the jet symmetry plane are neglected
    G = Go+Ge
    G = G[::-1, ::-1]

    return G


