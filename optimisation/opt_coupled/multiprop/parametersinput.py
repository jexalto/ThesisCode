import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class parameters(om.IndepVarComp):
    def initalize(self):
        pass

    def setup(self):
        pointmass = 0.5
        
        self.add_output('radius', val=0.1185, units='m') #0.14391204
        self.add_output("vinf", val=40., units="m/s")
        self.add_output("alpha", val=3., units="deg")
        self.add_output("Mach_number", val=0.2)
        self.add_output("re", val=1.0e5, units="1/m")
        self.add_output("rho", val=1.225, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=500, units="m")
        self.add_output("W0", val=5, units="kg")
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")
        self.add_output("span", val=0.748*0.96*2, units="m") # 0.748
        jetloc1 = 0.25
        jetloc2 = -jetloc1
        self.add_output("jet_loc", val=np.array(jetloc1), units="m")
        self.add_output("point_masses", val=np.array([pointmass, pointmass]), units="kg")

        chord_cp = np.ones(3)*0.24
        twist_cp = np.zeros(3)
        self.add_output("twist", shape=(len(twist_cp)), val=twist_cp, units="deg")
        self.add_output("chord", shape=(len(chord_cp)), val=chord_cp, units="m")

# Design Vars
# {'helix.geodef_parametric_0_rot_rate': array([-898.78035587]),
#  'helix.geodef_parametric_0_twist': array([70.        , 70.        , 62.88159985, 58.26579667, 57.81446586,
#        51.26115985, 48.49728126, 43.91624724, 41.26099458, 40.47737696,
#        38.89369462, 37.16629279, 35.85480893, 34.06475215, 32.00934689,
#        30.62607074, 29.59148632, 28.25777092, 27.30872414, 24.70314402]),
#  'parameters.chord': array([0.07      , 0.13457321, 0.19426579]),
#  'parameters.twist': array([-3.92595861,  1.00028374, -2.15346665])}
