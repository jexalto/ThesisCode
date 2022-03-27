import numpy as np
from RethorstCorrection_pyf90.mod_velocity_distribution_nonsymmetry import velocity_distribution_nosym
import openmdao.api as om

class Rethorst(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('panels_span_VLM', default=300, desc='number of spanwise panels on the VLM')
        self.options.declare('panels_chord_VLM', default=1, desc='number of chordwise panels')
        self.options.declare('span_max', default=40., desc='maximum span')
        self.options.declare('r_min', default=0.5, desc='minimum radius')
        self.options.declare('vel_distr_shape', default=1, desc='vel distr shape')

    def setup(self):
        panels_span_VLM = self.options['panels_span_VLM']
        panels_chord_VLM = self.options['panels_chord_VLM']
        vel_distr_shape = self.options['vel_distr_shape']

        self.add_input('span', val=1.0, units='m')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('Vinf', val=1.0, units='m/s')        
        self.add_input('radii', shape=(vel_distr_shape+1), units='m')
        self.add_input('velocity_vector', shape=(vel_distr_shape+1), units='m/s')

        self.add_output('correction_matrix',    shape   =   (panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), 
                                                val     =   np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM)))
        self.add_output('wing_veldistr',        shape   =   (panels_span_VLM*panels_chord_VLM), 
                                                val     =   np.zeros((panels_span_VLM)), units='m/s')

    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        span = inputs['span']
        jet_loc = inputs['jet_loc']
        radii_input = inputs['radii']
        Vinf = inputs['Vinf']
        prop_veldistr = inputs['velocity_vector']

        panels_span_VLM = self.options['panels_span_VLM']
        panels_chord_VLM = self.options['panels_chord_VLM']
        span_max = self.options['span_max']
        r_min = self.options['r_min']
        
        prop_discr = 20
        panels_jet = 101
        panels_overset_wing = 1001

        total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
        y_VLM = np.zeros((panels_span_VLM+1), order='F')
        vel_vec = np.zeros((panels_span_VLM), order='F')

        prop_veldistr = np.flip(np.array(prop_veldistr, order='F').reshape((len(prop_veldistr))))
        radii_input = np.array(radii_input, order='F').reshape((len(radii_input)))
        prop_veldistr = np.where(np.isnan(prop_veldistr), 0, prop_veldistr)

        print(radii_input)
        print(prop_veldistr)

        velocity_distribution_nosym(    span,
                                        jet_loc, prop_veldistr, radii_input,
                                        prop_discr,
                                        Vinf,
                                        panels_jet, panels_overset_wing,
                                        panels_chord_VLM, panels_span_VLM,
                                        span_max, r_min,
                                        y_VLM, vel_vec,
                                        total_correction)

        outputs['correction_matrix'] = total_correction
        outputs['wing_veldistr'] = vel_vec