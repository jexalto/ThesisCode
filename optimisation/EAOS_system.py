import numpy as np
from RethorstCorrection_pyf90.mod_overset_interpolation import overset_interpolation
import openmdao.api as om

from EOAS_analysis import EAOS_analysis_

class EOAS(om.ExplicitCompoenent):
    def initialize(self):
        self.options.declare('panels_VLM', default=71, desc='number of spanwise panels on the VLM')
        self.options.declare('nx', default=2, desc='number of chordwise panels')

    def setup(self):
        panels_VLM = self.options['panels_VLM']
        nx = self.options['nx']

        self.add_input('mesh', shape=(nx, panels_VLM), val=np.zeros((nx, panels_VLM)))
        self.add_input('correction_matrix', shape=(panels_VLM, panels_VLM), val=np.zeros((panels_VLM, panels_VLM)))
        self.add_input('chord', val=1.0, units='m')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('jet_radius', val=1.0, units='m')
        self.add_input('Vinf', val=1.0, units='m/s')
        self.add_input('Vjet', val=1.0, units='m/s')

        self.add_output('CL', val=1.0)
        self.add_output('CD', val=1.0)
        self.add_output('cl', shape=(1, panels_VLM), val=np.zeros((1, panels_VLM)))
        self.add_output('cd', shape=(1, panels_VLM), val=np.zeros((1, panels_VLM)))
        self.add_output('Wfuel', val=1.0)
        self.add_output('Wwing', val=1.0)

    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        span = inputs['span']
        jet_loc = inputs['jet_loc']
        jet_radius = inputs['jet_radius']
        Vinf = inputs['Vinf']
        Vjet = inputs['Vjet']
        mesh = inputs['mesh']
        correction = inputs['correction']

        panels_VLM = self.options['panels_VLM']
        panels_chord = self.options['nx']
        
        y = mesh[:, 0, 0]
        correction_loc = np.zeros((panels_VLM))

        for index in range(panels_VLM):
            correction_loc[index] = (y[index+1]+y[index])/2

        correction_loc = np.where(correction_loc>(jet_loc+jet_radius), 0, correction_loc)
        correction_loc = np.where(correction_loc<(jet_loc-jet_radius), 0, correction_loc)
        correction_loc = np.where(correction_loc==0, 0, 1)

        CL, CD, cl, cd, Wstruc, Wfuel = EAOS_analysis_( panels_VLM, panels_chord, jet_radius, jet_loc,\
                                                        Vinf, Vjet, mesh, correction_loc, correction)

        outputs['CL'], outputs['CD'] = CL, CD
        outputs['cl'], outputs['cd'] = cl, cd
        outputs['Wstruc'], outputs['Wfuel'] = Wstruc, Wfuel