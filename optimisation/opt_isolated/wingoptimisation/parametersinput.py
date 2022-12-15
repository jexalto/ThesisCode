import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class parameters(om.IndepVarComp):
    def initalize(self):
        pass

    def setup(self):
        pointmass = 0.5
        
        self.add_output('radius', val=0.1185, units='m')
        self.add_output("vinf", val=40., units="m/s")
        self.add_output("alpha", val=2., units="deg")
        self.add_output("Mach_number", val=0.2)
        self.add_output("re", val=1.0e5, units="1/m")
        self.add_output("rho", val=1.225, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=500, units="m")
        self.add_output("W0", val=2*pointmass, units="kg")
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")
        self.add_output("span", val=0.748*2, units="m") # 0.748
        jetloc1 = -0.25
        jetloc2 = -jetloc1
        self.add_output("jet_loc", val=np.array([jetloc1, jetloc2]), units="m")
        self.add_output("point_masses", val=np.array([pointmass, pointmass]), units="kg")

        wing_discr = 5
        chord_cp =np.ones(wing_discr)*0.24# [0.01,       1.04723251, 2.98191746, 5.,         5.,         2.98195628, 1.0472563,  0.01      ]
        twist_cp = np.zeros(wing_discr) #[-0.86296788, -1.11602075, -1.31543096, -0.95378036, -0.95378182, -1.31542815, -1.11602006, -0.86297512]

        self.add_output("twist", shape=(len(twist_cp)), val=twist_cp, units="deg")
        self.add_output("chord", shape=(len(chord_cp)), val=chord_cp, units="m")