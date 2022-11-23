import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class parameters(om.IndepVarComp):
    def initalize(self):
        pass

    def setup(self):
        pointmass = 0.5
        radius = 0.22339999999999993*4
        wing_discr = 5
        
        self.add_output('radius0', val=radius, units='m')
        self.add_output('radius1', val=radius, units='m')
        self.add_output("vinf", val=79.7389, units="m/s")
        self.add_output("alpha", val=1., units="deg")
        self.add_output("Mach_number", val=0.2)
        self.add_output("re", val=1.0e5, units="1/m")
        self.add_output("rho", val=0.907, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=500, units="m")
        self.add_output("W0", val=10, units="kg") # 1+2*pointmass
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")
        self.add_output("span", val=13.49, units="m") # 0.748
        jetloc1 = -2.3
        jetloc2 = -jetloc1
        self.add_output("jet_loc", val=np.array([jetloc1, jetloc2]), units="m")
        self.add_output("point_masses", val=np.array([pointmass, pointmass]), units="kg")

        chord_cp = np.ones(wing_discr)*0.6# 9369
        twist_cp = np.zeros(wing_discr)
        self.add_output("twist", shape=(len(twist_cp)), val=twist_cp, units="deg")
        self.add_output("chord", shape=(len(chord_cp)), val=chord_cp, units="m")