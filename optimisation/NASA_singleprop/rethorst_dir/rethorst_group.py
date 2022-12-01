import numpy as np
from RethorstCorrection_pyf90.mod_multiprop import multiprop
from RethorstCorrection_pyf90.mod_multiprop_d import multiprop_d
from RethorstCorrection_pyf90.mod_multiprop_b import multiprop_b

import openmdao.api as om

class Rethorst(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('panels_span_VLM', default=100, desc='number of spanwise panels on the VLM')
        self.options.declare('panels_chord_VLM', default=1, desc='number of chordwise panels')
        self.options.declare('span_max', default=20., desc='maximum span')
        self.options.declare('r_min', default=0.025, desc='minimum radius')
        self.options.declare('vel_distr_shape', default=100, desc='number of vel discretisation point given by helix')
        self.options.declare('prop_discr', default=5, desc='prop vel distribution discretisation')
        self.options.declare('panels_jet', default=61, desc='panels in jet in overset mesh')
        self.options.declare('panels_overset_wing', default=601, desc='panels on wing in overset mesh')
        self.options.declare('nr_props', default=2, desc='number of propellers')

    def setup(self):
        panels_span_VLM = self.options['panels_span_VLM']
        panels_chord_VLM = self.options['panels_chord_VLM']
        vel_distr_shape = self.options['vel_distr_shape']
        nr_props = self.options['nr_props']

        self.add_input('span', val=1.0, units='m')
        self.add_input('jet_loc', shape=(nr_props), units='m')
        self.add_input('vinf', val=1.0, units='m/s')        
        self.add_input('radii', shape_by_conn=True, units='m')
        self.add_input('velocity_vector', shape_by_conn=True, units='m/s')

        self.add_output('correction_matrix',    shape   =   (panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), 
                                                val     =   np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM)))
        self.add_output('wing_veldistr',        shape   =   (panels_span_VLM), 
                                                val     =   np.zeros((panels_span_VLM)), units='m/s')
        self.add_output('jet_radius', val=1., units='m')

    def compute(self, inputs, outputs):
        
        span                    = np.copy(inputs['span'])
        jet_loc_list            = np.copy(inputs['jet_loc'])
        Vinf                    = np.copy(inputs['vinf'][0])
        radii_input_            = np.copy(inputs['radii'])

        panels_span_VLM         = self.options['panels_span_VLM']
        panels_chord_VLM        = self.options['panels_chord_VLM']
        span_max                = self.options['span_max']
        r_min                   = self.options['r_min']
        prop_discr              = self.options['prop_discr']            # make sure the propeller isn't too small for this discretisation
        panels_jet              = self.options['panels_jet']
        panels_overset_wing     = self.options['panels_overset_wing']
        nr_props                = self.options['nr_props']
        nr_radii_input          = len(radii_input_[0])
        radii_input             = np.array(radii_input_, order='F')
        vel_distr_input = np.zeros(np.shape(inputs['velocity_vector']), order='F')
        vel_distr_input          = np.array(np.copy(inputs['velocity_vector']), order='F')
        # vel_distr_input[0][:]          = np.array(np.copy(inputs['velocity_vector'][0][::-1]), order='F')
        # vel_distr_input[1][:]         = np.array(np.copy(inputs['velocity_vector'][1][::-1]), order='F')
        # vel_distr_input         = vel_distr_input[::-1]
        
        total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
        vel_vec = np.zeros((panels_span_VLM), order='F')

        multiprop(  span, nr_props, jet_loc_list, vel_distr_input, radii_input, nr_radii_input, prop_discr, Vinf, panels_jet, \
                    panels_overset_wing, panels_chord_VLM, panels_span_VLM, span_max, r_min, vel_vec, total_correction)
        mat = np.matrix(total_correction.copy())

        if np.isnan(any(sum(total_correction))):
            raise ValueError('total correction contains nan!')
        if np.isnan(any(vel_vec)):
            raise ValueError('total correction contains nan!')
        # print()
        # print('=====================')
        # print('====== Rethorst =====')
        # print('=====================')
        # print(f'jet_loc {jet_loc_list}\n')
        
        outputs['correction_matrix'] = total_correction
        outputs['wing_veldistr'] = vel_vec

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
         # Clear all Seeds in Memory
        # start = time.time()
        self._zero_seeds(inputs)

        span                    = inputs['span']
        jet_loc_list            = np.array(inputs['jet_loc'], order='F')
        Vinf                    = np.array(inputs['vinf'], order='F')
        radii_input             = np.array(inputs['radii'], order='F')
        vel_distr_input         = np.array(inputs['velocity_vector'], order='F')
        nr_radii_input          = len(radii_input)
        panels_span_VLM         = self.options['panels_span_VLM']
        panels_chord_VLM        = self.options['panels_chord_VLM']
        span_max                = self.options['span_max']
        r_min                   = self.options['r_min']
        prop_discr              = self.options['prop_discr']
        panels_jet              = self.options['panels_jet']
        panels_overset_wing     = self.options['panels_overset_wing']
        nr_props                = self.options['nr_props']

        # ================================================
        # Forward Mode
        # ================================================
        total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
        vel_vec = np.zeros((panels_span_VLM), order='F')

        if mode=='fwd':
            self._set_seeds_fwd(d_inputs)

            multiprop_d(span, self.spand, nr_props, jet_loc_list, self.jet_loc_listd, vel_distr_input, self.vel_distr_inputd, radii_input, \
                        self.radii_inputd, nr_radii_input, prop_discr, Vinf, self.vinfd, panels_jet, panels_overset_wing, panels_chord_VLM, \
                        panels_span_VLM, span_max, r_min, vel_vec, self.vel_vecd, total_correction, self.total_correctiond)

            self._get_seeds_fwd(d_outputs)
            # if any(d_inputs['jet_loc']!=0):
            #     print('jetloc: ', jet_loc_list)
            #     print(d_inputs['jet_loc'])
            #     print(d_outputs['wing_veldistr'])
            # print('Rethorst derivatives:\tForward')

        # ================================================
        # Reverse Mode
        # ================================================

        # elif mode=='rev':
        #     print('START: Rethorst derivatives:\tReverse')

        #     self._set_seeds_rev(d_outputs)
            
        #     multiprop_b(span, self.spanb, nr_props, jet_loc_list, self.jet_loc_listb, vel_distr_input, self.vel_distr_inputb, radii_input, \
        #                 self.radii_inputb, nr_radii_input, prop_discr, Vinf, self.vinfb, panels_jet, panels_overset_wing, panels_chord_VLM, \
        #                 panels_span_VLM, span_max, r_min, vel_vec, self.vel_vecb, total_correction, self.total_correctionb)

        #     self._get_seeds_rev(d_inputs)

        #     print('END: Rethorst derivatives:\tReverse')

        
        # timediff = time.time() - start
        # print(f'total time computing {mode} Rethorst derivatives: {timediff}')

    def _zero_seeds(self, inputs):
        # ===== Clear Seeds =====
        panels_span_vlm                 = self.options['panels_span_VLM']
        panels_chord_vlm                = self.options['panels_chord_VLM']
        nr_props                        = self.options['nr_props']
        vel_distr_input                 = inputs['velocity_vector']
        
        # ===== Forward Seeds =====
        self.spand                      = np.array([0.], order='F')
        self.jet_loc_listd              = np.zeros(nr_props, order='F')
        self.vinfd                      = np.array([0.], order='F')
        self.vel_distr_inputd           = np.zeros((np.shape(inputs['velocity_vector'])), order='F')
        self.radii_inputd               = np.zeros((np.shape((inputs['velocity_vector']))), order='F')
        self.vel_vecd                   = np.zeros((panels_span_vlm), order='F')
        self.total_correctiond          = np.zeros((panels_span_vlm*panels_chord_vlm, panels_span_vlm*panels_chord_vlm), order='F')

        # ===== Reverse Seeds ===== 
        self.spanb                      = np.array([0.], order='F')
        self.jet_loc_listb              = np.zeros(nr_props, order='F')
        self.vinfb                      = np.array([0.], order='F')
        self.vel_distr_inputb           = np.zeros((np.shape(inputs['velocity_vector'])), order='F')
        self.radii_inputb               = np.zeros((np.shape(inputs['velocity_vector'])), order='F')
        self.vel_vecb                   = np.zeros((panels_span_vlm), order='F')
        self.total_correctionb          = np.zeros((panels_span_vlm*panels_chord_vlm, panels_span_vlm*panels_chord_vlm), order='F')

    def _set_seeds_fwd(self, d_inputs):
        self.spand                      = d_inputs['span']
        self.jet_loc_listd              = d_inputs['jet_loc']
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

