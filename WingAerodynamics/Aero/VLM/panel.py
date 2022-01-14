from geometry import area_4points
from vortices import (vortex_position_in_panel,
                      v_induced_by_horseshoe_vortex,
                      v_induced_by_horseshoe_vortex_farfield,
                      v_induced_by_horseshoe_vortex_3d)


class Panel(object):
    """Modified from: https://github.com/aqreed/PyVLM
    
             y ^
               |            Each panel is defined by the (x, y) coordinates
        P3--B--|--+--P4     of four points - namely P1, P2, P3 and P4 -
         |  |  |  |  |      ordered clockwise.
         |  |  |  |  |      Points defining the horseshoe are A, B and P.
         |  |  +--CP-|--->
         |  |     |  |   x
         |  |     |  |
        P2--A-----+--P1

    Args:
        P1,P2,P3,P4 (numpy.array): Corner points in a 2D euclidean space
    """

    def __init__(self, P1, P2, P3, P4):
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.P4 = P4
        
        #horseshoe vortex points
        self.CP, self.A, self.B = vortex_position_in_panel(P1, P2, P3, P4)
        #panel area
        self.area = area_4points(P1, P2, P3, P4)
        #panel span
        self.span = abs(P3[1] - P2[1])
        #panel slope
        self.slope = 0
        
        #position w.r.t. the local strip chord
        self.chordwise_position = 0
        #strip chord
        self.chord = 0
        #strip quarter chord point
        self.C4 = 0
        #strip leading edge point
        self.LE = 0
        
        #pertubation velocities added to Vinf
        self.Vx = 0
        self.Vz = 0
        
        self.V_eff = 0      #effective velocity
        self.alpha_eff = 0  #effective angle of attack
        self.Vinf_ind = 0   #induced velocity
        self.alpha_ind = 0  #induced angle of attack
        
        self.Vinf_n = 0     #local upstream normal velocity
        self.gamma = 0      #circulation value
        
        self.l = 0          #lift force
        self.d = 0          #drag force
        self.m = 0          #moment about quarter chord
        self.d_ = 0         #direct drag force
        self.cl = 0         #lift force coefficient
        self.cd = 0         #drag force coefficient

    def induced_velocity(self, control_point_pos):
        """
        Returns the induced velocity by a horseshoe vortex and the induced
        velocity excluding the bounded segment at a control point, defined
        as argument of the method, that does not have to be its own CP.
        
        Args:
            control_point_pos (numpy.array): control point in wing plane
            
        Returns:
            v_total (float): aerodynamic influence by full horseshoe vortex
            v_trail (float): aerodynamic influence by only the trailing vortices
        """

        v_total, v_trail = v_induced_by_horseshoe_vortex(control_point_pos,
                                                             self.A, self.B)

        return v_total, v_trail
    
    def induced_velocity_farfield(self, control_point_pos):
        """
        Returns the induced velocity by the trailing vortices of a horseshoe
        vortex at a control point, with the control point projected in the
        farfield, infinitely far downstream
        
        Args:
            control_point_pos (numpy.array): control point in wing plane
            
        Returns:
            v_ind (float): aerodynamic influence by horseshoe vortex
        """
        
        v_ind = v_induced_by_horseshoe_vortex_farfield(control_point_pos,
                                                           self.A, self.B)
        
        return v_ind
    
    def induced_velocity_3d(self, pos):
        """
        Returns the induced velocity by a horseshoe vortex and the induced
        velocity excluding the bounded segment at any point in 3D space.
        
        Args:
            pos (numpy.array): a point in 3D space
            
        Returns:
            v_total (numpy.array): aerodynamic influence by full horseshoe
                vortex as a 3D vector
            v_trail (numpy.array): aerodynamic influence by only the trailing
                vortices as a 3D vector
        """
        
        v_total, v_trail = v_induced_by_horseshoe_vortex_3d(pos, self.A,
                                                            self.B, self.chord)
        
        return v_total, v_trail
