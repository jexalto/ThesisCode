import numpy as np
import openmdao.api as om

class printsys(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        self.add_input('propspan', shape_by_conn=True)
        self.add_input('propchord', shape_by_conn=True)

    def compute(self, inputs, outputs):
        propspan = inputs['propspan']
        propchord = inputs['propchord']

        print('span: ', propspan)
        print('chord: ', propchord)