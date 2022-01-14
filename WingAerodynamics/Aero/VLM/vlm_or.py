import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as si

from mesh_generator import Mesh
from airfoils import NACA4, Airfoil
from jet_correction import jet_correction


class PyVLM(object):
    """Given a geometry, mesh chordwise and spanwise densities, angle of
    attack, upstream velocity, applies the VLM theory to the
    defined lifting surface.
    
    Modified from: https://github.com/aqreed/PyVLM
    """

    def __init__(self):
        
        #vlm solver
        self.Points = []        #list of points to contruct the panels
        self.Panels = []        #list of panels
        self.AIC = 0            #Aerodynamic Influence coefficient matrix
        self.AIC_i = 0          #AIC to calculate induced drag
        self.G = 0              #jet correction matrix
        self.gamma = 0          #solution for circulation values
        
        #results
        self.Strip = 0          #strip properties
        self.res = 0            #integrated results
        
        #inputs
        self.rho = 1.225        #density [kg/m3]
        self.Vinf = 30          #Vinf [m/s]
        self.alpha = 2          #angle of attack [deg]
        #airfoil
        self.airfoil = Airfoil()    
        self.NACA = None        #if value is provided, NACA airfoil will be used
        self.polar_dir = 'data_aero_wing/' #folder with viscous polars of the airfoil
        #pertubation velocities
        #x,y: x,y location control point
        #Vx,Vz: induced velocity on control point
        #Vzi: induced velocity in farfield
        #V_aoa (not an input): velocity used by jet correction
        #dalpha (not an input): twist applied by viscous correction
        #turb: determines if a section is turbulent or not
        self.Vpert = pd.DataFrame(columns=['x', 'y', 'Vx', 'Vz', 'Vzi',
                                           'V_aoa', 'dalpha', 'turb'])
        
        #constants
        self.mu = 1.7894e-05    #dynamic viscosity [Pa s]
        self.a = 340            #speed of sound [m/s]
        
        #settings
        self.jet_corr_settings = {'dlbda':      0.1,
                                  'dlb':        0.005,
                                  'p_max':      10,
                                  'lbda_max':   5}

    def reset(self):
        """Resets the instance so a new wing can be added"""
        
        self.Points = []
        self.Panels = []
        self.AIC = 0
        self.AIC_i = 0
        self.G = 0
        self.gamma = 0
        self.Strip = 0
        self.res = 0
        self.Vpert = pd.DataFrame(columns=['x', 'y', 'Vx', 'Vz', 'turb'])

    def add_wing(self, lead_edge_coord, chord_lengths, n, m, opt=None):
        """Allows the addition of a wing to the mesh, defined by its chords'
        lengths and leading edges locations. The spanwise and chordwise
        density of the mesh can be controlled through n and m.
        ONLY half a wing is needed to define it. A specular image will
        be used to create the other half.

        Args:
            lead_edge_coord (list of numpy.array): Coordinates of the leading
                edge points as arrays in a 2D euclidean space
            chord_lengths (list): Chord lenghts corresponding to the sections 
                defined by the leading edge coordinates
            n (list): number of chordwise panels per section
            m (list): number of spanwise panels per section
            opt (list): spacing per section, options are:
                'cos', 'sin', 'nsin', 'uni'
        """
        
        # clears AIC and results when modifying the mesh
        self.reset()

        if len(lead_edge_coord) != len(chord_lengths):
            msg = 'Same number of chords and leading edges required'
            raise ValueError(msg)

        # MESH GENERATION
        # When more than two chords -with their respectives leading
        # edges coordinates- are provided, it iterates through the
        # lists containing both location and length.

        Nle = len(lead_edge_coord)
        
        if opt is None:
            opt = ['uni']*(Nle-1)
        
        for k in range(Nle - 1):
            leading_edges = [lead_edge_coord[k],
                             lead_edge_coord[k + 1]]
            chords = [chord_lengths[k],
                      chord_lengths[k + 1]]

            # The mesh is created taking into account the desired
            # mesh density spanwise -"n"- and chordwise -"m"-

            mesh = Mesh(leading_edges, chords, n[k], m[k])

            # The points of the mesh and its panels - sets of 4 points
            # orderly arranged - are calculated

            Points_ = mesh.points(opt[k])
            if k==0:
                Panels_ = mesh.panels()
            else:
                Panels_ = mesh.panels(count=m[k-1])

            self.Points.extend(Points_)
            self.Panels.extend(Panels_)

        # Specular image to generate the opposite semi-span of the wing

        for k in range(Nle - 1):
            leading_edges = [lead_edge_coord[k]*[1, -1],
                             lead_edge_coord[k + 1]*[1, -1]]
            chords = [chord_lengths[k],
                      chord_lengths[k + 1]]

            mesh = Mesh(leading_edges, chords, n[k], m[k])

            Points_ = mesh.points(opt[k])
            if k==0:
                Panels_ = mesh.panels(ismirror=True)
            else:
                Panels_ = mesh.panels(ismirror=True, count=m[k-1])

            self.Points.extend(Points_)
            self.Panels.extend(Panels_)
        
        #determine quarter chord points on wing
        c4 = []
        for panel in self.Panels:
            c4.append(panel.C4)
            
        c4 = np.array(c4)
        c4 = np.unique(c4, axis=0)
        
        #create pertubation data frame
        self.Vpert['x'] = c4[:, 0]
        self.Vpert['y'] = c4[:, 1]
        self.Vpert['Vx'] = 0
        self.Vpert['Vz'] = 0
        self.Vpert['Vzi'] = 0
        self.Vpert['V_aoa'] = 0
        self.Vpert['dalpha'] = 0
        self.Vpert['turb'] = False
        self.Vpert = self.Vpert.sort_values(by=['y'])

    def check_mesh(self, print_mesh=False, plot_mesh=False):
        """
        Prints the points of the mesh, the disposition of each panel and
        plots them for visual check.
        
        Args:
            print_mesh, plot_mesh (boolean): Self-explained
        """

        Points = self.Points
        Panels = self.Panels

        # Check for coincident points
        N = len(Points)
        for i in range(N):
            count = 0
            for j in range(N):
                if(((Points[j] == Points[i]).all()) is True):
                    count += 1
                    if(count > 1):
                        msg = "Two points of the mesh coincide"
                        raise ValueError(msg)

        # Check for incorrectly defined panels
        N = len(Panels)
        for i in range(N):
            P1P2 = Panels[i].P2 - Panels[i].P1
            P1P3 = Panels[i].P3 - Panels[i].P1
            P1P4 = Panels[i].P4 - Panels[i].P1
            P3P4 = Panels[i].P4 - Panels[i].P3

            i_inf = np.array([1, 0])

            if np.cross(P1P2, i_inf) != 0:
                msg = 'P1P2 segment not aligned with OX'
                raise ValueError(msg)

            if np.cross(P1P2, P3P4) != 0:
                msg = 'Panel incorrectly defined, P1P2 and P3P4 not parallel'
                raise ValueError(msg)

            if np.sign(np.cross(P1P2, P1P3)) != np.sign(np.cross(P1P3, P1P4)):
                msg = 'Points not in a clockwise/counterclockwise fashion'
                raise ValueError(msg)

        # PRINTING AND PLOTTING
        if (print_mesh is True):
            print('\nPanel| Chrd% |  Span |  Points coordinates')
            print('------------------------------------------------')
            for i in range(N):
                print(' %3s | %5.2f | %5.3f | '
                      % (i, 100*Panels[i].chordwise_position, Panels[i].span),
                      np.round(Panels[i].P1, 2), np.round(Panels[i].P2, 2),
                      np.round(Panels[i].P3, 2), np.round(Panels[i].P4, 2))

        if (plot_mesh is True):
            plt.style.use('ggplot')
            plt.xlim(-5, 15), plt.ylim(-10, 10)
            for i in range(len(Points)):
                P = Points[i]
                plt.plot(P[0], P[1], 'ro')
            plt.show()
            
    def vlm_visc(self):
        """This uses a viscous correction on the lift using the 2D lift polar.
        With the angle of attack determined by the VLM solution, CL must be
        increased or decreased based on the difference between the 2D inviscid
        and viscous polars, resulting in a local change in wing twist.
        """
        
        #first time solving VLM
#        self.Vpert['dalpha'] = 0
        dalpha_0 = self.Vpert['dalpha']
        self.vlm()
        self.airfoil.init_polars(self.polar_dir)
        self.strip_properties()
        
        Strip = self.Strip
        
        #calculation of the correction
        dalpha = []
        #loop over wing sections
        for y in Strip.index:
            #get local properties
            alpha_eff = Strip.loc[y, 'alpha_eff']
            Re = Strip.loc[y, 'Re']
            M = Strip.loc[y, 'M']
                        
            #determine viscous and inviscid CL the wing section
            if Strip.loc[y, 'turb']:
                cl_visc = self.airfoil.polar_turbulent.cl(M, alpha_eff, Re)
            else:
                cl_visc = self.airfoil.polar_laminar.cl(M, alpha_eff, Re)
            
            cl_inv = self.airfoil.cl_thin_airfoil(alpha_eff)
            
            #translate difference in CL to twist
            dcl_da = 2*np.pi*np.pi/180
            da = (cl_visc-cl_inv)/dcl_da
            
            dalpha.append(da)
            
        dalpha = np.array(dalpha)
        
        #add twist and solve VLM again
        self.Vpert['dalpha'] = dalpha_0+dalpha
        self.vlm()
        self.strip_properties()
        self.Strip['dalpha'] = dalpha
        self.Vpert['dalpha'] = dalpha_0

    def vlm(self):
        """
        For a given set of panels applies the VLM theory:

        1) Calculates the induced velocity produced by all the
            associated horseshoe vortices of strength=1 on each panel,
            calculated on its control point where the boundary condition
            will be imposed.
        2) Computes the circulation by solving the linear equation.
        3) Calculates the aerodynamic forces.
        """

        Panels = self.Panels
        rho = self.rho
        Vinf = self.Vinf
        Vpert = self.Vpert
        alpha = np.deg2rad(self.alpha)
        q_inf = (1 / 2) * rho * (Vinf**2)

        # 1. BOUNDARY CONDITION
        # To impose the boundary condition we must calculate the normal
        # components of (a) induced velocity "Wn" by horshoe vortices of
        # strength=1 and (b) upstream normal velocity "Vinf_n"

        #   (a) INDUCED VELOCITIES
        #     - "Wn", normal component of the total induced velocity by
        #       the horshoe vortices, stored in the matrix "AIC" where the
        #       element Aij is the velocity induced by the horshoe vortex
        #       in panel j on panel i
        #     - also the induced velocity by *only* the trailing vortices
        #        "Wi" on panel i is calculated and stored in the Panel object
        #       attribute "accul_trail_induced_vel"

        N = len(Panels)

        if (type(self.AIC) == int):
            # In case it is the first time the method is called, it proceeds
            # to compute the AIC matrix

            AIC = np.zeros((N, N))    #Aerodynamic Influence Coefficient matrix
            AIC_i = np.zeros((N, N))  #AIC for induced drag calculation 

            for i, panel_pivot in enumerate(Panels):
                CP = panel_pivot.CP

                for j, panel in enumerate(Panels):
                    Wn,__ = panel.induced_velocity(CP)
                    Wi = panel.induced_velocity_farfield(CP)
                    AIC[i, j] = Wn      #induced normal velocity by horshoe vortices
                    AIC_i[i,j] = Wi/2   #induced normal velocity by trailing vortices

            self.AIC = AIC
            self.AIC_i = AIC_i
        
        #   (b) UPSTREAM NORMAL VELOCITY
        #     It will depend on the angle of attack -"alpha"- and the camber
        #     gradient at each panel' position within the local chord

        Vinf_n = np.zeros(N)  # upstream (normal) velocity
        
        if self.NACA is not None:
            self.airfoil = NACA4(self.NACA)
        
        airfoil = self.airfoil
        
        Vpert['V_aoa'] = np.cos(alpha)*Vpert['Vx']+np.sin(alpha)*Vpert['Vz']
        
        for i,panel in enumerate(Panels):
            position = panel.chordwise_position
            panel.slope = airfoil.camber_gradient(position)
            
            Vpert_ = Vpert.loc[Vpert['y']==panel.C4[1]]
            panel.Vx = Vpert_['Vx'].iloc[0]
            panel.Vz = Vpert_['Vz'].iloc[0]
            panel.Vzi = Vpert_['Vzi'].iloc[0]
            
            dalpha = np.deg2rad(Vpert_['dalpha'].iloc[0])
            
            panel.Vinf_n = Vinf * (panel.slope-dalpha-alpha) + panel.slope*panel.Vx - panel.Vz
            Vinf_n[i] = panel.Vinf_n
            
        # Rethorst correction
        if (type(self.G) == int):
            y = []
            panel_chord = []
            for panel in Panels:
                y.append(panel.P2[1])
                y.append(panel.P3[1])
                panel_chord.append(panel.chord)
            
            panel_chord = np.mean(panel_chord)
            y = np.array(y)
            y = np.unique(y)
            y = y[y>=0]
            
            V_r = Vpert['V_aoa'].loc[Vpert['y']>0].to_numpy()
            V_l = Vpert['V_aoa'].loc[Vpert['y']<0].to_numpy()
            V_l = V_l[::-1]
            if np.max(V_r)!=0 and np.max(V_l)!=0:
                if np.max(np.abs(V_l-V_r)/V_r)>0.01:
                    print('Warning: pertubation velocity not symmetrical')
            
            G_r = 0
            for i in range(1, y.size-1):
                mu_r = (Vinf+V_r[i-1])/(Vinf+V_r[i])
                r = y[-1]-y[i]
                if round(mu_r, 2)!=1:
                    Gi = jet_correction(y, mu_r, panel_chord, r,
                                        self.jet_corr_settings['dlbda'],
                                        self.jet_corr_settings['dlb'],
                                        self.jet_corr_settings['p_max'],
                                        self.jet_corr_settings['lbda_max'])
                    G_r += Gi
        
            G = np.zeros((N, N))
            
            for i, panel_pivot in enumerate(Panels):
                i_ = panel_pivot.idx
                yi = panel_pivot.CP[1]
                if type(G_r)!=int:        
                    for j, panel in enumerate(Panels):
                        j_ = panel.idx
                        if yi*panel.CP[1]>0:
                            G[i, j] = G_r[i_, j_]
        
            self.G = G/(4*np.pi)
        
        # 2. CIRCULATION (Γ or gamma)
        # by solving the linear equation (AX = Y) where X = gamma
        # and Y = Vinf_n
                
        self.gamma = np.linalg.solve(self.AIC-self.G, Vinf_n)

        # 3. AERODYNAMIC FORCES
        L = 0
        D = 0
        S = 0

        #Calculate induced normal velocities
        V_i = np.matmul(self.AIC_i-self.G, self.gamma)
        
        #Calculate forces and aerodynamic parameters with circulation
        for i,panel in enumerate(Panels):
            panel.gamma = self.gamma[i]
            
            Vx = Vinf*np.cos(alpha) + panel.Vx
            Vz = Vinf*np.sin(alpha) + panel.Vz + V_i[i]
            
            Fx = -Vz * rho * panel.gamma * panel.span
            Fz = Vx * rho * panel.gamma * panel.span
            
            panel.l = Fz * np.cos(alpha) - Fx * np.sin(alpha)
            #drag calculated on the wing
            panel.d_ = Fx * np.cos(alpha) + Fz * np.sin(alpha)
            #drag calculated downstream in Trefftz plane
            panel.d = -rho*panel.gamma*(V_i[i]+panel.Vzi)*panel.span
            
            panel.cl = panel.l / (q_inf * panel.area)
            panel.cd = panel.d / (q_inf * panel.area)
            
            panel.V_ind = V_i[i] + panel.Vz
            panel.alpha_ind = np.rad2deg(panel.V_ind/Vinf)  # induced AoA(deg)
            
            panel.V_eff = (Vx**2 + Vz**2)**0.5
            panel.alpha_eff = np.rad2deg(np.arctan2(Vz, Vx))
            
            #Summation for total forces
            L += panel.l
            D += panel.d
            S += panel.area

        CL = L / (q_inf * S)
        CD = D / (q_inf * S)
        
        #Save results
        res = {'CL': CL,
               'CD': CD,
               'CDp': 0,
               'CDi': CD,
               'L': L,
               'D': D,
               'Dp': 0,
               'Di': D,
               'alpha': np.rad2deg(alpha),
               'rho': rho,
               'Vinf': Vinf,
               'S': S}
        
        self.res = res
    
    def induced_velocity(self, P):
        """Calculates the induced velocity by the wing on a point. Due to the
        way the VLM is set up, induced velocity can only be calculated on the
        control points. This leads to three cases:
            
        1) When the point is on the wing, the induced velocity is calculated by
            interpolating between control points.
        2) When the point is outside the wing, but still in the span region,
            the induced velocity must be calculated by interpolating between
            the y coordinates of the control points.
        3) The last case is when the point is located outside the span region
            of the wing, here the induced velocity can be directly calculated.
            
        Args:
            P (numpy.array): point of calculation
                
        Returns:
            V_i (numpy.array): induced velocity
        """
        
        strip = self.Strip
        y = strip.index
        
        #Check if point is within span
        if P[1]>y.min() and P[1]<y.max():
            #Get the closest y control point coorindates
            y_u = y[y>P[1]].min()
            y_l = y[y<P[1]].max()
            
            x_u = np.array(strip.loc[y_u, 'x_coord'])
            x_l = np.array(strip.loc[y_l, 'x_coord'])
            
            #Check if point is on the wing
            if P[0]>x_l.min() and P[0]<x_l.max() and P[0]>x_u.min() and P[0]<x_u.max():
                
                #Get the closes x control point coordinates
                x_uu = x_u[x_u>P[0]]
                x_ul = x_u[x_u<P[0]]
                x_lu = x_l[x_l>P[0]]
                x_ll = x_l[x_l<P[0]]
                
                #Create 4 control points
                P1 = np.array([x_uu.min(), y_u, P[2]])
                P2 = np.array([x_ul.max(), y_u, P[2]])
                P3 = np.array([x_lu.min(), y_l, P[2]])
                P4 = np.array([x_ll.max(), y_l, P[2]])
            else:
                #Create 4 control points
                P1 = np.array([P[0]+0.01, y_u, P[2]])
                P2 = np.array([P[0]+0.01, y_l, P[2]])
                P3 = np.array([P[0]-0.01, y_u, P[2]])
                P4 = np.array([P[0]-0.01, y_l, P[2]])
    
            P_i = [P1, P2, P3, P4]
            x_i = [P1[0], P2[0], P3[0], P4[0]]
            y_i = [P1[1], P2[1], P3[1], P4[1]]
            
            u_i = []
            v_i = []
            w_i = []
            
            #Calculate induced velocities
            for i in range(len(x_i)):
                V_i = self.sp_induced_velocity(P_i[i])
                u_i.append(V_i[0])
                v_i.append(V_i[1])
                w_i.append(V_i[2])
            
            #Interpolate
            f_u = si.interp2d(x_i, y_i, u_i)
            f_v = si.interp2d(x_i, y_i, v_i)
            f_w = si.interp2d(x_i, y_i, w_i)
            
            u = f_u(P[0], P[1])
            v = f_v(P[0], P[1])
            w = f_w(P[0], P[1])
            
            V_i = np.array([u[0],v[0],w[0]])
            
        else:
            V_i = self.sp_induced_velocity(P)
        
        return V_i
                
    def sp_induced_velocity(self, P):
        """Calculate the induced velocity. Creates the AIC matrix and
        multiplies with the circulation.
        
        Args:
            P (numpy.array): point of calculation
                
        Returns:
            V_ind (numpy.array): induced velocity
        """
        
        Panels = self.Panels
        gamma  = self.gamma
        
        #create AIC matrix
        N = len(Panels)
        AIC_P = np.zeros([3, N])
        for i, panel in enumerate(Panels):
            Wn, Wi = panel.induced_velocity_3d(P)
            AIC_P[0, i] = Wn[0]
            AIC_P[1, i] = Wn[1]
            AIC_P[2, i] = Wn[2]          
            
        V_ind = np.matmul(AIC_P, gamma)

        return V_ind
    
    def strip_properties(self, skin_fric=True):
        """Sums the results of the panels per strip.
        
        Args:
            skin_fric (bool): to include viscous drag
        """
        
        #get flow properties
        rho = self.rho
        Vinf = self.Vinf
        mu = self.mu
        alpha = np.deg2rad(self.alpha)
        q_inf = (1 / 2) * rho * (Vinf**2)
        
        Panels = self.Panels
        
        #sort based on y values
        y = []
        for panel in Panels:
            y.append(panel.CP[1])
        y = np.unique(y)
        
        #values to be determined
        cols = ['l', 'cl', 'd', 'cd', 'di', 'di_', 'cdi', 'cdi_', 'dp', 'cdp', 'area', 'span',
                'chord', 'gamma', 'V_eff', 'alpha_eff', 'V_ind', 'Vx', 'Vz',
                'Re', 'C4', 'LE', 'x_coord', 'turb']
        #summation over panels
        cols_sum = ['l', 'di', 'di_', 'area', 'gamma']
        #properties that are the same for all panels on a strip
        cols_set = ['span', 'chord', 'C4', 'LE']
        #properties that need to be interpolated to quarter chord point
        cols_int = ['Vx', 'Vz', 'V_ind']
        
        data = pd.DataFrame(columns=cols, index=y)
        
        #set values 
        for col in cols_sum:
            data[col].values[:] = 0
        for y in data.index:    
            data.loc[y, 'x_coord'] = []
            for var in cols_int:
                data.loc[y, var] = []
        
        for panel in Panels:
            loc = panel.CP[1]
            #summation
            for var in cols_sum:
                if var=='di':
                    var2 = 'd'
                elif var=='di_':
                    var2 = 'd_'
                else:
                    var2 = var
                data.loc[loc, var] += getattr(panel, var2)
            
            #set value to value of first panel in strip
            for var in cols_set:
                test = pd.isna(data.loc[loc, var])
                if type(test)==np.ndarray:
                    test = False
                if test:
                    data.loc[loc, var] = getattr(panel, var)
            
            #put values in a list
            for var in cols_int:
                data.loc[loc, var].append(getattr(panel, var))                    
            data.loc[loc, 'x_coord'].append(panel.CP[0])
        
        #integrate data to quarted chord point
        for y in data.index:
            x = data.loc[y, 'C4'][0]
            x_lst = data.loc[y, 'x_coord']
            
            for var in cols_int:
                var_lst = data.loc[y, var]
                if len(x_lst)>1:
                    f = si.interp1d(x_lst, var_lst)
                    var_val = float(f(x))
                else:
                    var_val = var_lst[0]
                    
                data.loc[y, var] = var_val
        
        #calculate normalized values
        data['cl'] = data['l']/(data['area']*q_inf)
        data['cdi'] = data['di']/(data['area']*q_inf)
        data['cdi_'] = data['di_']/(data['area']*q_inf)
        
        #velocity dependent values
        Vx = Vinf*np.cos(alpha) + data['Vx']
        Vz = Vinf*np.sin(alpha) + data['Vz'] + data['V_ind']
        data['V_eff'] = (Vx**2 + Vz**2)**0.5
        data['alpha_eff'] = np.rad2deg(np.arctan2(Vz.astype('f'),
                                                    Vx.astype('f')))        
        data['M'] = data['V_eff']/self.a
        data['Re'] = rho*data['V_eff']*data['chord']/mu
    
        #turbulence        
        data['turb'] = self.Vpert['turb'].to_numpy()
        
        #add viscous drag
        if skin_fric:
            if self.airfoil.polar_turbulent is None:
                self.airfoil.init_polars(self.polar_dir)
            polar_t = self.airfoil.polar_turbulent
            polar_l = self.airfoil.polar_laminar
            
            for y in data.index:
                if data.loc[y, 'turb']:                    
                    data.loc[y, 'cdp'] = polar_t.profile_drag(
                                                    data.loc[y, 'M'],
                                                    data.loc[y, 'cl'],
                                                    data.loc[y, 'alpha_eff'],
                                                    data.loc[y, 'Re'])
                    
                else:
                    data.loc[y, 'cdp'] = polar_l.profile_drag(
                                                    data.loc[y, 'M'],
                                                    data.loc[y, 'cl'],
                                                    data.loc[y, 'alpha_eff'],
                                                    data.loc[y, 'Re'])
        else:
            data['cdp'] = 0
        
        #final drag calculations
        data['dp'] = data['cdp']*(data['area']*q_inf)
        data['d'] = data['dp']+data['di']
        data['cd'] = data['cdp']+data['cdi']
        
        #save results
        self.Strip = data
        
        #calculate total drag forces
        res = self.res
        res['Dp'] = data['dp'].sum()
        res['CDp'] = res['Dp']/(q_inf*res['S'])
        res['D'] = res['Di']+res['Dp']
        res['CD'] = res['CDi']+res['CDp']
        self.res = res
    
    def plot_strip_properties(self, save=False):
        """Plots strip properties"""
        
        strip = self.Strip
        y = strip.index
        
        dct = {'cl': '$C_l$ [-]',
               'cd': '$C_d$ [-]',
               'gamma': '$\Gamma$ [$m^2/s$]',
               'V_eff': '$V_{eff}$ [m/s]',
               'alpha_eff': r'$\alpha_{eff}$ [deg]',
               'Vx': '$V_{perturb,x}$ [m/s]',
               'Vz': '$V_{perturb,z}$ [m/s]'}
        
        for var in dct.keys():
            plt.figure()
            plt.plot(y, strip[var])
            if var=='cd':
                plt.plot(y, strip['cdi'])
                plt.plot(y, strip['cdp'])
                plt.legend(['$C_{d,t}$', '$C_{d,i}$', '$C_{d,p}$'])
            plt.xlabel('y [m]')
            plt.ylabel(dct[var])
            plt.xlim([min(y), max(y)])
            plt.grid()
            
            if save:
                plt.savefig('%s.png' % var,
                            dpi=300,
                            bbox_inches='tight')
        
    def plot_points(self):
        """Plots panel points"""
        
        plt.figure()
        for point in self.Points:
            plt.plot(point[1], point[0], marker='o', color='r')
        plt.gca().set_aspect('equal')
        plt.ylabel('x [m]')
        plt.xlabel('y [m]')

    def uniform_inflow(self):
        """Sets uniform inflow"""
        
        self.G = 0
        for panel in self.Panels:
            panel.Vx = 0
            panel.Vz = 0
            panel.Vzi = 0
        self.Vpert['Vx'] = 0
        self.Vpert['Vz'] = 0
        self.Vpert['Vzi'] = 0
        self.Vpert['dalpha'] = 0
        self.Vpert['turb'] = False
    