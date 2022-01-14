import numpy as np


def vortex_position_in_panel(P1, P2, P3, P4):
    """For a given panel defined by points P1, P2, P3 and P4
    returns the position of the horseshoe vortex defined
    by points A, B and its control point P.
    
    Modified from: https://github.com/aqreed/PyVLM

            ^
          y |                Points defining the panel
            |                are named clockwise.
     P3--B--|-----P4
      |  |  |     |
      |  |  |     |
     T1  |  +--P--T2---->
      |  |        |     x
      |  |        |
     P2--A--------P1

    Args:
        P1,P2,P3,P4 (numpy.array): Points that define the panel

    Returns:
        P (numpy.array): control point where the boundary condition V*n = 0
            is applied according to the Vortice Lattice Method.
        A,B (numpy.array): points that define the horseshoe position
    """

    P2P1 = P1 - P2
    P3P4 = P4 - P3
    P2P3 = P3 - P2
    P1P4 = P4 - P1

    T1 = P2 + P2P3 / 2
    T2 = P1 + P1P4 / 2
    T1T2 = T2 - T1

    A = P2 + P2P1 / 4
    B = P3 + P3P4 / 4
    P = T1 + (3/4) * T1T2
    
    return P, A, B

def v_induced_by_horseshoe_vortex(P, A, B):
    """
    Induced velocity at point P due to a horseshoe vortex
    of strenght gamma=1 spatially positioned by points A and B,
    extended to x_Inf(+) in a 2D euclidean space. Circulation
    direction is: x_Inf(+) -> A -> B -> x_Inf(+)
    
    Modified from: https://github.com/aqreed/PyVLM

                ^
              y |                Points defining the horseshoe
    V_inf       |                are named clockwise.
    ->     B----|->--+...>...    A direction vector is
    ->     |    |    |           calculated for each vortex.
    ->     ^    +----|------>
    ->	   |         |       x
    ->	   A----<----+...<...

    Args:
        P (numpy.array): point of reference
        A,B (numpy.array): points of the horseshoe vortex

    Returns:
        v_total, v_trail (float): aerodynamic influence with and without bound
            vortex
    """

    pi = np.pi

    a = P[0] - A[0]
    b = P[1] - A[1]
    c = P[0] - B[0]
    d = P[1] - B[1]
    e = (a**2 + b**2)**0.5
    f = (c**2 + d**2)**0.5
    g = B[0] - A[0]
    h = B[1] - A[1]

    div = a*d - c*b
    if (div == 0):
        v_bounded = 0
    else:
        v_bounded = (1/div) * (((g*a + h*b)/e) - ((g*c + h*d)/f))
        v_bounded /= 4*pi  #induced velocity in P due to bounded vortex

    if (b == 0):
        v_trail1 = 0
    else:
        v_trail1 = -(a + e)/(b*e)

    if (d == 0):
        v_trail2 = 0
    else:
        v_trail2 = (c + f)/(d*f)

    v_trail = v_trail1 + v_trail2
    v_trail /= 4*pi  #induced velocity in P due to trailing vortices

    v_total = v_trail + v_bounded  #total induced velocity in P
    
    return v_total, v_trail

def v_induced_by_horseshoe_vortex_farfield(P, A, B):
    """Induced velocity at point P due to the trailing vortices of a 
    horseshoe vortex of strenght gamma=1 spatially positioned by
    points A and B, extended to x_Inf(+) and x_Inf(-), as if point P
    is projected infinitely far downstream.
    
    Modified from: https://github.com/aqreed/PyVLM

                ^
              y |                Points defining the horseshoe
    V_inf       |                are named clockwise.
    ->-->--B----|->--+...>...    A direction vector is
    ->     |    |    |           calculated for each vortex.
    ->     ^    +----|------>
    ->	   |         |       x
    ->--<--A----<----+...<...

    Args:
        P (numpy.array): point of reference
        A,B (numpy.array): points of the horseshoe vortex

    Returns:
        v_ind (float): aerodynamic influence of the vortex
    """
    
    pi = np.pi
    
    a = P[1] - A[1]
    b = P[1] - B[1]
    
    v_ind1 = -1/a
    v_ind2 = 1/b
    
    v_ind = v_ind1+v_ind2
    v_ind /= 2*pi
    
    return v_ind
    

def v_induced_by_horseshoe_vortex_3d(P, A, B, chord):
    """
    Induced velocity at point P due to a horseshoe vortex
    of strenght gamma=1 spatially positioned by points A and B,
    extended to 1000 times the local chord in a 3D euclidean space.
    Circulation direction is: x_Inf(+) -> A -> B -> x_Inf(+)

                ^
              y |                Points defining the horseshoe
    V_inf       |                are named clockwise.
    ->     B----|->--+...>...    A direction vector is
    ->     |    |    |           calculated for each vortex.
    ->     ^    +----|------>
    ->	   |         |       x
    ->	   A----<----+...<...

    Args:
        P (numpy.array): point of reference
        A,B (numpy.array): points of the horseshoe vortex
        chord (float): local chord

    Returns:
        v_total, v_trail (numpy.array): aerodynamic influence with and without
            bound vortex
    """
    
    tol = 1e-12  	   #size of viscous core
    len_w = chord*1000 #length of trailing vortices
    
    #convert to 3D coorindates
    A = np.append(A,0)
    B = np.append(B,0)
    C = A+np.array([len_w,0,0])
    D = B+np.array([len_w,0,0])
    
    #calculate influence due to bound and trailing vortices
    v_b = v_induced_vortex_3d(P,A,B,tol)
    v_w1 = v_induced_vortex_3d(P,C,A,tol)
    v_w2 = v_induced_vortex_3d(P,B,D,tol)

    v_trail = v_w1+v_w2
    v_total = v_b + v_trail

    return v_total, v_trail

def v_induced_vortex_3d(P, A, B, tol):
    """
    Induced velocity at point P due to a vortex line of strenght gamma=1 
    starting at point A and ending at point B
    
    Method from: Low speed aerodynamics from wing theory to panel methods, 
    J Katz & A Plotkin 1991, Section 10.4.5
    
    Args:
        P (numpy.array): point of reference
        A,B (numpy.array): start and end points of the vortex
        tol (float): viscous core dimension where induced velocity = 0
    
    Returns:
        V (numpy.array): aerodynamic influence
    """
    
    pi = np.pi
    
    xp = P[0]
    yp = P[1]
    zp = P[2]
    x1 = A[0]
    y1 = A[1]
    z1 = A[2]
    x2 = B[0]
    y2 = B[1]
    z2 = B[2]
    
    #cross product
    r1xr2x = (yp-y1)*(zp-z2)-(zp-z1)*(yp-y2)
    r1xr2y = -(xp-x1)*(zp-z2)+(zp-z1)*(xp-x2)
    r1xr2z = (xp-x1)*(yp-y2)-(yp-y1)*(xp-x2)
    
    absr1xr22 = r1xr2x**2+r1xr2y**2+r1xr2z**2
    
    r1 = ((xp-x1)**2+(yp-y1)**2+(zp-z1)**2)**0.5
    r2 = ((xp-x2)**2+(yp-y2)**2+(zp-z2)**2)**0.5
    
    if r1<tol or r2<tol or absr1xr22<tol:
        u,v,w = 0,0,0
    else:
        #dot product
        r0dr1 = (x2-x1)*(xp-x1)+(y2-y1)*(yp-y1)+(z2-z1)*(zp-z1)
        r0dr2 = (x2-x1)*(xp-x2)+(y2-y1)*(yp-y2)+(z2-z1)*(zp-z2)
        
        K = 1/(4*pi*absr1xr22)*(r0dr1/r1-r0dr2/r2)
        
        u = K*r1xr2x
        v = K*r1xr2y
        w = K*r1xr2z
        
    V = np.array([u,v,w])
    
    return V
    