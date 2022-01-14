import numpy as np
import scipy.special as sp


def jet_correction_ib(y, mu, crd, loc, ib_loc, dlbda=0.1, dlb=0.005, p_max=10, lbda_max=5):
    """Calculates the jet correction proposed by Rethorst, based on a velocity
    ratio mu. For this function the original correction has been applied to 
    the wing, where the center of the jet coincides with the nearest panel centre.
    
    S.C. Rethorst (1958) Aerodynamic of Nonuniform Flows as Related to an
    Airfoil Extending Through a Circular Jet
    
    Args:
        y (numpy.array): spanwise wing locations
        mu (float): velocity ratio V0/Vj
        crd (float): wing (mean) chord
        loc (float): location
        ib_loc (float): location of inboard propeller
        dlbda (float): integral discretization step
        dlb (float): integral discretization step
        p_max (int): upper limit sum
        lbda_max (float): upper limit integral
    
    Returns:
        G (numpy.array): correction matrix
    """
    #calculate geometry
    b = y[-1]*2             #span
    s0 = min((y[:-1]+y[1::])/2 - ib_loc, key=abs)  # propeller y offset, needs to be in centre of panel

    r0 = ib_loc-y[np.argmin(np.abs(y-loc))]  # radius from edge to edge equals abs(ib_loc - loc)

    #y = np.concatenate((np.sort(-y), y[1:]))   # invert y and sort from -b/2 to 0
    y_ = (y[1:]+y[:-1])/2   # centre of each panel
    
    #calculations are made in a system where the jet is on y=0
    y2_ = -(y[::-1] - ib_loc)      # locations of control point
    y2 = -(y_[::-1] - ib_loc)  # locations of hb vortex points
    r0_ = r0-s0                 # radius from mid panel jet axis to panel boundary jet boundary

    n = np.where(y2 >= 0)[0][0]  # only positive side of wing (towards tip)

    #control point y coordinates
    eta = y2[n:]/r0_
    eta = eta.reshape(eta.size, 1)
    #horseshoe vortex midpoitn y coordinates
    beta = y2[n:]/r0_
    #horseshoe vortex points
    c = y2_[n:-1]/r0_
    d = y2_[n+1:]/r0_

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

    # Influence of single horseshoe vortex
    Ge[:,0] *= 0.5 #At beta=0 there's only one hs vortex
    G_ = Go + Ge
    G_2 = Go + Ge
    G_ *= 0.5
    G_[:, 0] *= 2

    N = len(y_)              # number of panels on 1 wing half
    G = np.zeros((N, N))    # influence matrix
    N0 = N-n                # number of panels on calculation side

    #apply the correction on the outboard side of the wing, and mirror in the jet symmetry axis  #!!!!!!!
    #to get the correction for the inboard side, neglect other wing half
    #G_ = Go+Ge
    G[:N0, :N0] = G_[::-1, ::-1]                    # invert?
    G[N0:, N0:] = G_[1:N - N0 + 1, 1:N - N0 + 1]
    G[N0:, :N0] = G_[1:N - N0 + 1, ::-1]
    G[:N0, N0:] = G_[::-1, 1:N - N0 + 1]
    #G = G[::-1, ::-1]

    #G_t = G + G[::-1, ::-1]#invert

    return G

if __name__ == '__main__':
    y = np.array([-5.        , -4.64705882, -4.29411765, -3.94117647, -3.58823529,
       -3.23529412, -2.88235294, -2.52941176, -2.17647059, -1.82352941,
       -1.47058824, -1.11764706, -0.76470588, -0.41176471, -0.05882353,
        0.29411765,  0.64705882,  1.        ,  1.22222222,  1.44444444,
        1.66666667,  1.88888889,  2.11111111,  2.33333333,  2.55555556,
        2.77777778,  3.        ,  3.33333333,  3.66666667,  4.        ,
        4.33333333,  4.66666667,  5.        ])
    mu = 2/3
    crd = 1
    loc = 1
    ib_loc = 2
    G = jet_correction_ib(y,mu,crd,loc,ib_loc)
