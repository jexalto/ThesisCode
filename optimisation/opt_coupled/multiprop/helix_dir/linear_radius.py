from ast import Pass
import numpy as np
import openmdao.api as om

class linear_radius(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("nr_sections", default=20, desc='Number of spanwise sections')
    
    def setup(self):
        nr_sections = self.options["nr_sections"]
        self.add_input('radius', val=0., units='m')
        self.add_input('jet_loc', val=0., units='m')
        self.add_input('chord', val=[1, 1, 1], units='m')
        self.add_input('twist', val=[1, 1, 1], units='deg')
        self.add_input('mesh', shape_by_conn=True)
        self.add_input('radii', shape_by_conn=True)
        
        self.add_output('propspan_sectional', val=np.zeros(nr_sections-1))
        self.add_output('jet_loc_list', val=[1,1])
        self.add_output('chord_list', val=[1, 1, 1, 1, 1])
        self.add_output('twist_list', val=[1, 1, 1, 1, 1])
        
        self.declare_partials('propspan_sectional', 'radius')
        self.declare_partials('jet_loc_list', 'jet_loc', val=[1, -1])
        self.declare_partials('chord_list', 'chord', rows=[0, 1, 2, 3, 4], cols=[0, 1, 2, 1, 0], val=[1, 1, 1, 1, 1])
        self.declare_partials('twist_list', 'twist', rows=[0, 1, 2, 3, 4], cols=[0, 1, 2, 1, 0], val=[1, 1, 1, 1, 1])

    def compute(self, inputs, outputs):
        nr_sections = self.options["nr_sections"]
        radius = inputs["radius"]
        jet_loc = inputs["jet_loc"]
        chord = inputs["chord"]
        twist = inputs["twist"]
        mesh = inputs['mesh'][0, :, 1]
        radii = inputs['radii'][-1]
        
        with open(f'00_results/meshresults/mesh_{radius}.txt', 'w') as file:
            np.savetxt(file, mesh, fmt='%.8f')
            file.write(str(radii))
            
        with open(f'00_results/meshresults/diff_mesh_{radius}.txt', 'w') as file:
            for i in range(0, len(mesh)-1):
                val = mesh[i+1]-mesh[i-1]
                file.write(str(val)+'\n')

        outputs['propspan_sectional'] = radius/nr_sections * np.ones(nr_sections-1)
        outputs['jet_loc_list'] = [jet_loc, -jet_loc]
        outputs['chord_list'] = [chord[0], chord[1], chord[2], chord[1], chord[0]]
        outputs['twist_list'] = [twist[0], twist[1], twist[2], twist[1], twist[0]]

    def compute_partials(self, inputs, partials):
        nr_sections = self.options["nr_sections"]

        partials['propspan_sectional', 'radius'] = 1/nr_sections
