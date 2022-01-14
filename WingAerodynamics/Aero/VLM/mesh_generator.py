import numpy as np

from panel import Panel


class Mesh(object):
    """Modified from: https://github.com/aqreed/PyVLM
    
   Pi +......> y     Given a trapezoid defined by vertices Pi and Pf
      | \            and chords 1 and 2 representing a wing segment
      |  \           that complies with the VLM theory, returns the
      |   + Pf       points and panels of the mesh:
chord1|   |
      |   |             - Points are given as a list of Np elements,
      |   |chord2         being Np the number of points of the mesh.
      +---+
      |                 - Panels are given as a list of list of NP
     x 				      elements, each element composed of 4 points,
                          where NP is the number of panels of the mesh.
    x - chordwise
        direction    The corner points of each panel are arranged in a
    y - spanwise     clockwise fashion following this order:
        direction
                             P2 +---+ P3...> y
                                |   |
                             P1 +   + P4
                                |
                                x
    """

    def __init__(self, leading_edges, chords, n, m):
        """length(leading_edges)=length(chords)=length(n)-1=length(m)-1
        
        Args:
            leading_edges (list of numpy.arrays): coordinates of the leading
                edge points Pi and Pf, as arrays in a 2D euclidean space (x, y)
            chords (list): chord lenghts values
            n (list): nº of chordwise panels
            m (list): nº of spanwise panels
        """
        
        self.leading_edges = leading_edges
        self.chords = chords
        self.n = n
        self.m = m
        self.mesh_points = []
        self.mesh_panels = []

    def points(self, opt):
        """Yields a list of length (n+1)*(m+1) containing equally spaced
        points (x, y) coordinates (arrays), for each trapezoid geometry
        defined by the arguments.
        
        Args:
            opt (str): spacing option for spanwise panels, there are four options:
                cos: cosine spacing
                sin: sine spacing
                nsin: inverse sine spacing
                uni: uniform spacing
                
        Returns:
            mesh_points (list of numpy.array): list of points
        """

        Pi = self.leading_edges[0]
        Pf = self.leading_edges[1]

        chord_1 = np.array([self.chords[0], 0])
        chord_2 = np.array([self.chords[1], 0])

        n = self.n
        m = self.m
        
        if opt=='cos':
            ang_0 = 0
            ang_1 = np.pi
            scale = 0.5
        elif opt=='sin':
            ang_0 = 0
            ang_1 = np.pi/2
            scale = 1
        elif opt=='nsin':
            ang_0 = np.pi/2
            ang_1 = -np.pi/2
            scale = 1
        elif opt=='uni':
            ang_0 = 0
            ang_1 = np.pi * m
            scale = 1 / (2 * m)

        ang_ = 0
        for i in range(n + 1):
            PiPf = Pf - Pi
            P = Pi
            ang = ang_0
            for j in range(m):
                self.mesh_points.append(P)
                ang2 = ang + ang_1 / m
                step = scale*np.abs(np.cos(ang) - np.cos(ang2))
                P = P + PiPf * step
                ang  = ang2
            self.mesh_points.append(Pf)
            ang2_ = ang_ + np.pi/n
            Pi = Pi + 0.5*chord_1 * (np.cos(ang_) - np.cos(ang2_)) #chord_1 / n
            Pf = Pf + 0.5*chord_2 * (np.cos(ang_) - np.cos(ang2_)) #chord_2 / n
            ang_ = ang2_
        return self.mesh_points

    def panels(self, ismirror=False, count=0):
        """Yields a list of length (n*m) containing the panels (objects),
        defined by 4 points previously calculated. The points are properly
        arranged to serve as locations for the horseshoe vortices.
        
        Args:
            ismirror (bool): must be true when creating panels for the wing
                half with negative spanwise coordinates
            count (int): starting number used when numbering the spanwise
                location of the panel
        
        Returns:
            mesh_panels (list of Panel instances): list of panels
        """

        Pi = self.leading_edges[0]
        chord_1 = np.array([self.chords[0], 0])

        n = self.n
        m = self.m

        N_panels = n * m

        for i in range(N_panels):
            # Panels list generation from mesh_points list
            k = int(i / m)

            if not ismirror:
                P1 = self.mesh_points[i + k + m + 1]
                P2 = self.mesh_points[i + k]
                P3 = self.mesh_points[i + k + 1]
                P4 = self.mesh_points[i + k + m + 2]
            else:
                P1 = self.mesh_points[i + k + m + 2]
                P2 = self.mesh_points[i + k + 1] 
                P3 = self.mesh_points[i + k]
                P4 = self.mesh_points[i + k + m+ 1]

            self.mesh_panels.append(Panel(P1, P2, P3, P4))

            # Chordwise position calculation
            if not ismirror:
                P1_ = self.mesh_panels[k * m].P1
                P2_ = self.mesh_panels[k * m].P2
            else:
                P1_ = self.mesh_panels[k * m].P4
                P2_ = self.mesh_panels[k * m].P3
            
            chord = P1_ - P2_
            panel_center = P2_ + chord / 2
            

            leadEdge_2_panel = np.linalg.norm(panel_center - Pi)
            relative_pos = leadEdge_2_panel / np.linalg.norm(chord_1)
            
            # Chord and quarted chord point calculation
            g = i % m
            if not ismirror:
                P1_ = self.mesh_points[g + m*n + n]
                P2_ = self.mesh_points[g]
                P3_ = self.mesh_points[g + 1]
                P4_ = self.mesh_points[g + m*n + n + 1]
            else:
                P1_ = self.mesh_points[g + m*n + n + 1]
                P2_ = self.mesh_points[g + 1]
                P3_ = self.mesh_points[g]
                P4_ = self.mesh_points[g + m*n + n]
            
            P2P3 = P3_ - P2_
            P1P4 = P4_ - P1_
        
            T1 = P2_ + P2P3 / 2
            T2 = P1_ + P1P4 / 2
            T1T2 = T2 - T1
        
            C4 = T1 + (1/4) * T1T2
            
            self.mesh_panels[i].chordwise_position = relative_pos
            self.mesh_panels[i].chord = T1T2[0]
            self.mesh_panels[i].C4 = C4
            self.mesh_panels[i].LE = T1
            self.mesh_panels[i].idx = g+count

        return self.mesh_panels
