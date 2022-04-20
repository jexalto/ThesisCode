import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class parameters(om.IndepVarComp):
    def initalize(self):
        pass

    def setup(self):
        self.add_output("vinf", val=40., units="m/s")
        self.add_output("alpha", val=2., units="deg")
        self.add_output("Mach_number", val=0.84)
        self.add_output("re", val=1.0e6, units="1/m")
        self.add_output("rho", val=0.38, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=11.165e6, units="m")
        self.add_output("W0", val=0.4 * 3e5, units="kg")
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")

        self.add_output("span", val=0.748, units="m")
        self.add_output("jet_loc", val=np.array([0.1]), units="m")

        twist_cp = np.zeros((5))
        chord_cp = np.ones((5))*0.48

        self.add_output("twist", shape=(5), val=twist_cp, units="deg")
        self.add_output("chord", shape=(5), val=chord_cp, units="m")