import numpy as np
import openmdao.api as om

class proploc(om.ExplicitComponent):
    def initialize(self):
        self

    def setup(self);
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('jet_radius', val=1.o, units='m')
        self.add_input('span', val=1.o, units='m')