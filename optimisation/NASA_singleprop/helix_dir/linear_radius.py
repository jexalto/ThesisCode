from ast import Pass
import numpy as np
import openmdao.api as om

class linear_radius(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("nr_sections", default=5, desc='Number of spanwise sections')
    
    def setup(self):
        nr_sections = self.options["nr_sections"]
        self.add_input('radius', val=0., units='m')
        
        self.add_output('propspan_sectional', val=np.zeros(nr_sections-1))
        
        self.declare_partials('propspan_sectional', 'radius')#, method='fd')

    def compute(self, inputs, outputs):
        nr_sections = self.options["nr_sections"]
        radius = inputs["radius"]
        # print(f"===========================\n========== Radius =========\n===========================\n\
        #     Radius={radius}")
        outputs['propspan_sectional'] = radius/(nr_sections-1) * np.ones(nr_sections-1)

    def compute_partials(self, inputs, partials):
        nr_sections = self.options["nr_sections"]

        partials['propspan_sectional', 'radius'] = 1/(nr_sections-1)
