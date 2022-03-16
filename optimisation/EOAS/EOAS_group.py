from bezier import Surface
import numpy as np
import openmdao.api as om

from parametersinput import parameters
from EOAS_analysis import surface
from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class EOAS(om.Group):

    def setup(self):
        self.add_subsystem('parameters', subsys=parameters(),       promotes_outputs=[(     'vinf', 'alpha', 'Mach_number',
                                                                                            're', 'rho', 'CT', 'R',
                                                                                            'W0', 'speed_of_sound', 'load_factor',
                                                                                            'empty_cg')])

        surface= {}
        
        aerostruct_group = AerostructGeometry(surface=surface)

        name = "wing"

        # Add tmp_group to the problem with the name of the surface.
        self.add_subsystem(name, aerostruct_group, promotes_inputs=[('span_scaling', 'chord_scaling')])

        point_name = "AS_point_0"

        # Create the aero point group and add it to the model
        AS_point = AerostructPoint(surfaces=[surface])

        self.add_subsystem(
            point_name,
            AS_point,
            promotes_inputs=[
                "v",
                "vjet",
                "r0",
                "alpha",
                "Mach_number",
                "re",
                "rho",
                "CT",
                "R",
                "W0",
                "speed_of_sound",
                "empty_cg",
                "load_factor",
                'correction',
                'correction_loc',
            ],
        )