import numpy as np
from RethorstCorrection_pyf90.mod_velocity_distribution_nonsymmetry import velocity_distribution_nosym
from RethorstCorrection_pyf90.mod_velocity_distribution_nonsymmetry_d import velocity_distribution_nosym_d
from RethorstCorrection_pyf90.mod_velocity_distribution_nonsymmetry_b import velocity_distribution_nosym_b

import openmdao.api as om

import time

class Rethorst(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('panels_span_VLM', default=100, desc='number of spanwise panels on the VLM')
        self.options.declare('panels_chord_VLM', default=1, desc='number of chordwise panels')
        self.options.declare('span_max', default=20., desc='maximum span')
        self.options.declare('r_min', default=0.1, desc='minimum radius')
        self.options.declare('vel_distr_shape', default=100, desc='number of vel discretisation point given by helix')
        self.options.declare('prop_discr', default=5, desc='prop vel distribution discretisation')
        self.options.declare('panels_jet', default=41, desc='panels in jet in overset mesh')
        self.options.declare('panels_overset_wing', default=801, desc='panels on wing in overset mesh')

    def setup(self):
        panels_span_VLM = self.options['panels_span_VLM']
        panels_chord_VLM = self.options['panels_chord_VLM']
        vel_distr_shape = self.options['vel_distr_shape']

        self.add_input('span', val=1.0, units='m')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('vinf', val=1.0, units='m/s')        
        self.add_input('radii', shape_by_conn=True, units='m')
        self.add_input('velocity_vector', shape_by_conn=True, units='m/s')

        self.add_output('correction_matrix',    shape   =   (panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), 
                                                val     =   np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM)))
        self.add_output('wing_veldistr',        shape   =   (panels_span_VLM*panels_chord_VLM), 
                                                val     =   np.zeros((panels_span_VLM)), units='m/s')
        self.add_output('jet_radius', val=1., units='m')

    def compute(self, inputs, outputs):
        start = time.time()
        span                    = np.copy(inputs['span'])
        jet_loc                 = np.copy(inputs['jet_loc'])
        Vinf                    = np.copy(inputs['vinf'][0])
        radii_input             = np.copy(inputs['radii'])
        prop_veldistr           = np.copy(inputs['velocity_vector'])

        panels_span_VLM         = self.options['panels_span_VLM']
        panels_chord_VLM        = self.options['panels_chord_VLM']
        span_max                = self.options['span_max']
        r_min                   = self.options['r_min']
        prop_discr              = self.options['prop_discr']            # make sure the propeller isn't too small for this discretisation
        panels_jet              = self.options['panels_jet']
        panels_overset_wing     = self.options['panels_overset_wing']

        total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
        vel_vec = np.zeros((panels_span_VLM), order='F')

        # prop_veldistr = np.array(np.flip(np.array(prop_veldistr).reshape((len(prop_veldistr)))), order='F')
        # radii_input = np.array(radii_input, order='F').reshape((len(radii_input)))

        velocity_distribution_nosym(    span, jet_loc, prop_veldistr, radii_input, prop_discr, Vinf, panels_jet,
                                        panels_overset_wing, panels_chord_VLM, panels_span_VLM, span_max, r_min, vel_vec, total_correction)
        
        if np.isnan(any(sum(total_correction))):
            raise ValueError('total correction contains nan!')
        if np.isnan(any(vel_vec)):
            raise ValueError('total correction contains nan!')

        print('\n============')
        print('==Rethorst==')
        print('============')
        print(f'span {span}')
        print(f'jet_loc {jet_loc}')
        print(f'radius {radii_input[-1]}')
        
        outputs['correction_matrix'] = total_correction
        outputs['wing_veldistr'] = vel_vec

        timediff = time.time() - start
        # print(f'total time computing Rethorst: {timediff}')

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
         # Clear all Seeds in Memory
        start = time.time()
        self._zero_seeds(inputs)

        span                    = inputs['span']
        jet_loc                 = inputs['jet_loc']
        Vinf                    = np.array(inputs['vinf'], order='F')
        radii_input             = np.array(inputs['radii'], order='F')
        prop_veldistr           = np.array(inputs['velocity_vector'], order='F')

        panels_span_VLM         = self.options['panels_span_VLM']
        panels_chord_VLM        = self.options['panels_chord_VLM']
        span_max                = self.options['span_max']
        r_min                   = self.options['r_min']
        prop_discr              = self.options['prop_discr']
        panels_jet              = self.options['panels_jet']
        panels_overset_wing     = self.options['panels_overset_wing']

        # ================================================
        # Forward Mode
        # ================================================
        total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
        vel_vec = np.zeros((panels_span_VLM), order='F')

        if mode=='fwd':
            self._set_seeds_fwd(d_inputs)

            velocity_distribution_nosym_d(  span, self.spand, jet_loc, self.jet_loc_inputd, prop_veldistr, self.vel_distr_inputd,
                                            radii_input, self.radii_inputd, prop_discr, Vinf, self.vinfd, panels_jet, panels_overset_wing, panels_chord_VLM,
                                            panels_span_VLM, span_max, r_min, vel_vec, self.vel_vecd, total_correction, self.total_correctiond)

            self._get_seeds_fwd(d_outputs)
            # print('Rethorst derivatives:\tForward')

        # ================================================
        # Reverse Mode
        # ================================================

        elif mode=='rev':
            self._set_seeds_rev(d_outputs)
            
            velocity_distribution_nosym_b(  span, self.spanb, jet_loc, self.jet_loc_inputb, prop_veldistr, self.vel_distr_inputb,
                                            radii_input, self.radii_inputb, prop_discr, Vinf, self.vinfb, panels_jet, panels_overset_wing, panels_chord_VLM,
                                            panels_span_VLM, span_max, r_min, vel_vec, self.vel_vecb, total_correction, self.total_correctionb)

            self._get_seeds_rev(d_inputs)
            # print('Rethorst derivatives:\tReverse')
        
        timediff = time.time() - start
        print(f'total time computing Rethorst derivatives: {timediff}')

    def _zero_seeds(self, inputs):
        # ===== Clear Seeds =====
        panels_span_vlm                 = self.options['panels_span_VLM']
        panels_chord_vlm                = self.options['panels_chord_VLM']
        vel_distr_input                 = inputs['velocity_vector']
        radii_input                     = inputs['radii']
        
        # ===== Forward Seeds =====
        self.spand                      = np.array([0.], order='F')
        self.jet_loc_inputd             = np.array([0.], order='F')
        self.vinfd                      = np.array([0.], order='F')
        self.vel_distr_inputd           = np.zeros((len(vel_distr_input)), order='F')
        self.radii_inputd               = np.zeros((len(radii_input)), order='F')
        self.vel_vecd                   = np.zeros((panels_span_vlm*panels_chord_vlm), order='F')
        self.total_correctiond          = np.zeros((panels_span_vlm*panels_chord_vlm, panels_span_vlm*panels_chord_vlm), order='F')

        # ===== Reverse Seeds ===== 
        self.spanb                      = np.array([0.], order='F')
        self.jet_loc_inputb             = np.array([0.], order='F')
        self.vinfb                      = np.array([0.], order='F')
        self.vel_distr_inputb           = np.zeros((len(vel_distr_input)), order='F')
        self.radii_inputb               = np.zeros((len(radii_input)), order='F')
        self.vel_vecb                   = np.zeros((panels_span_vlm*panels_chord_vlm), order='F')
        self.total_correctionb          = np.zeros((panels_span_vlm*panels_chord_vlm, panels_span_vlm*panels_chord_vlm), order='F')

    def _set_seeds_fwd(self, d_inputs):
        self.spand                      = d_inputs['span']
        self.jet_loc_inputd             = d_inputs['jet_loc']
        self.vinfd                      = d_inputs['vinf']
        self.vel_distr_inputd           = np.array(d_inputs['velocity_vector'], order='F')
        self.radii_inputd               = np.array(d_inputs['radii'], order='F')

    def _get_seeds_fwd(self, d_outputs):
        if "correction_matrix" in d_outputs:
            d_outputs['correction_matrix']  += self.total_correctiond
        if "wing_veldistr" in d_outputs:
            d_outputs['wing_veldistr']      += self.vel_vecd

    def _set_seeds_rev(self, d_outputs):
        self.total_correctionb          = np.array(d_outputs['correction_matrix'], order='F')
        self.vel_vecb                   = np.array(d_outputs['wing_veldistr'], order='F')

    def _get_seeds_rev(self, d_inputs):
        d_inputs['span']                += self.spand
        d_inputs['jet_loc']             += self.jet_loc_inputd
        d_inputs['vinf']                += self.vinfd
        d_inputs['velocity_vector']     += self.vel_distr_inputd
        d_inputs['radii']               += self.radii_inputd

