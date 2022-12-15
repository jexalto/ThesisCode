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
        self.add_output("alpha", val=2., units="deg")
        self.add_output("Mach_number", val=0.2)
        self.add_output("re", val=1.0e5, units="1/m")
        self.add_output("rho", val=1.225, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=500, units="m")
        self.add_output("W0", val=7.5, units="kg")
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")
        self.add_output("span", val=0.748*0.96*2, units="m") # 0.748
        jetloc1 = 0.2 #205
        jetloc2 = -jetloc1
        self.add_output("jet_loc", val=np.array(jetloc1), units="m")
        self.add_output("point_masses", val=np.array([pointmass, pointmass]), units="kg")


        chord_cp = np.ones(3) * 0.24 #[0.15      , 0.15      , 0.36088911]#np.ones(3)*0.24 # [0.05,       0.05,       0.46085366] # np.ones(3)*0.24 # [0.01,       0.01,       0.46085366] #np.ones(3)*0.24 #    [0.13797784, 0.07402913, 0.26110817] # , 0.05865475, 0.01500252] #np.ones(3)*0.24 #
        twist_cp = [1, 1, 1] # np.zeros(3) #[1, 2, 2]#np.ones(3)#[0.09612332, 0.19836825, 0.23178217] # np.zeros(3)# [-1, 1, 2] #np.ones(3)*-1 #[0.00034206, 0.00182013, 0.00113233] # np.zeros(3) #[0.42134592, 0.74691498, 0.88536901]#, 0.74691498, 0.42134592] #np.zeros(3) # np.array([-2, -1, 0, -1,  -2])
        self.add_output("twist", shape=(len(twist_cp)), val=twist_cp, units="deg")
        self.add_output("chord", shape=(len(chord_cp)), val=chord_cp, units="m")