from email.policy import default
import numpy as np
import openmdao.api as om

class propinflow(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('nx', default=2, desc='Number of chordwise discretisation points')
        self.options.declare('ny', default=101, desc='Number of spanwise discretisation points')
        self.options.declare('nr_props', default=1, desc='Number of propellers')
        self.options.declare('propdist_chord', default=0.2, desc='Chordwise distance of prop versus leading edge')

    def setup(self):
        nr_props = self.options['nr_props']

        self.add_input('circulation', shape_by_conn=True)
        self.add_input('jet_loc', shape_by_conn=True)
        self.add_input('mesh', shape_by_conn=True)

        self.add_output('propinflow', val=np.zeros(nr_props))

        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        nx              = self.options['nx']
        ny              = self.options['ny']
        propdist_chord  = self.options['propdist_chord']

        mesh            = np.copy(inputs['mesh'])
        y_              = np.copy(inputs['mesh'][0, :, 1])
        jet_loc         = np.copy(inputs['jet_loc'])
        circulation     = np.copy(inputs['circulation'])

        y = np.zeros((len(y_)-1))
        for i in range(len(y_)-1):
            y[i] = 0.5*(y_[i]+y_[i+1])

        induced_aoa = np.zeros(len(jet_loc))

        for index, iJetloc in enumerate(jet_loc):
            iCirculation    = np.argmin(abs(y-iJetloc))
            
            v_ind = 0
            d_chord = 0
            for inx in range(nx-1):
                prop_circ = circulation[iCirculation+ny*(inx)]
                d_chord += (mesh[inx+1, 0, 0] - mesh[inx, 0, 0])/2
                d_total = d_chord+propdist_chord

                v_ind += -prop_circ/d_total
            
            induced_aoa[index] = v_ind
