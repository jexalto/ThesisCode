import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant

class parameters(om.IndepVarComp):
    def initalize(self):
        pass

    def setup(self):
        pointmass = 0.5
        spanwisesections = 4
        radius = np.linspace(0.0125, 0.09, 40)
        veldistr = self._parabola(np.linspace(-5, 5, 40))

        print(veldistr)
        
        self.add_output('radii', val=radius, units='m')#, units='m')
        self.add_output("vel_distr", val=veldistr)#, units='m/s')
        self.add_output("vinf", val=40., units="m/s")
        self.add_output("alpha", val=3., units="deg")
        self.add_output("Mach_number", val=0.2)
        self.add_output("re", val=1.0e5, units="1/m")
        self.add_output("rho", val=1.225, units="kg/m**3")
        self.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_output("R", val=500, units="m")
        self.add_output("W0", val=10, units="kg") # converged correct, 12, 23
        self.add_output("speed_of_sound", val=295.4, units="m/s")
        self.add_output("load_factor", val=1.)
        self.add_output("empty_cg", val=np.zeros((3)), units="m")
        self.add_output("span", val=0.748*0.96*2, units="m") # 0.748
        jetloc1 = 0.2#1886+0.0204/2
        self.add_output("jet_loc", val=np.array(jetloc1), units="m")
        self.add_output("point_masses", val=np.array([pointmass, pointmass]), units="kg")

        chord_cp = np.ones(spanwisesections)*0.24 #[0.01      , 0.05257288, 0.10513638, 0.21928196] #ones(spanwisesections)*0.24
        twist_cp = np.zeros(spanwisesections) # [10.        , 10.        , 10.        ,  8.99935772] #zeros(spanwisesections)
        self.add_output("twist", shape=(spanwisesections), val=twist_cp, units="deg")
        self.add_output("chord", shape=(spanwisesections), val=chord_cp, units="m")

    def _parabola(self, x): return -x**2+40+x[0]**2