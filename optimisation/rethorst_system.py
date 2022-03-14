import numpy as np
from RethorstCorrection_pyf90.mod_overset_interpolation import overset_interpolation
import openmdao.api as om

class Rethorst(om.ExplicitCompoenent):
    def initialize(self):
        self.options.declare('panels_VLM', default=71, desc='number of spanwise panels on the VLM')
        self.options.declare('nx', default=1, desc='number of chordwise panels')
        self.options.declare('span_max', default=1, desc='maximum span')
        self.options.declare('r_min', default=1, desc='minimum radius')

    def setup(self):
        panels_VLM = self.options['panels_VLM']
        nx = self.options['nx']

        self.add_input('span', val=1.0, units='m')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('jet_radius', val=1.0, units='m')
        self.add_input('Vinf', val=1.0, units='m')
        self.add_input('Vjet', val=1.0, units='m')

        self.add_output('correction_matrix', shape=(panels_VLM, panels_VLM), val=np.zeros((panels_VLM, panels_VLM)))
        self.add_output('mesh', shape=(nx, panels_VLM), val=np.zeros((nx, panels_VLM)))

    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        span = inputs['span']
        jet_loc = inputs['jet_loc']
        jet_radius = inputs['jet_radius']
        Vinf = inputs['Vinf']
        Vjet = inputs['Vjet']

        panels_VLM = self.options['panels_VLM']
        nx = self.options['nx']
        span_max = self.options['span_max']
        r_min = self.options['r_min']

        correction_matrix = np.zeros((panels_VLM, panels_VLM), order='F')
        mesh = np.zeros((panels_VLM+1, 3, nx), order='F')
        panels_jet = 91
        panels_overset_wing = 3101

        overset_interpolation(span, 1., jet_loc, jet_radius, Vinf, Vjet, panels_jet,\
        panels_overset_wing, nx, panels_VLM, span_max, r_min, correction_matrix, mesh)

        outputs['correction_matrix'] = correction_matrix
        outputs['mesh'] = mesh