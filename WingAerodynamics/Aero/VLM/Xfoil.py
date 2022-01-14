import numpy as np
import pandas as pd
import os
import shutil
import subprocess
import time


class Xfoil(object):
    
    def __init__(self):
        self.visc = False
        
        self.Re = 1000000
        self.M = 0
        self.N = 9
        self.x_trip_up = 1
        self.x_trip_bot = 1
        self.vacc = 0.01
        self.NACA = None
        self.airfoil = None
        self.N_iter_max = 500
        self.N_panel = 300
        
        self.run_dir = 'test1/'
        self.xfoil_dir = 'data_xfoil/'
        self.output_file = 'output.dat'
        
        self.cl = None
        self.cd = None
        self.cdp = None
        self.cm = None
        
    def run_xfoil(self, alpha):
        
        if not os.path.exists(self.run_dir):
            os.mkdir(self.run_dir)
        
        f = open(self.run_dir+'/xfoil_run.dat', 'w')
        
        f.write('plop\n')
        f.write('g,f\n')
        f.write('\n')
        
        #Load airfoil
        if self.airfoil is None:
            f.write('naca %s\n' % self.NACA)
        elif self.NACA is None:
            f.write('load %s.dat\n' % self.airfoil)
            f.write('ppar\n')
            f.write('N %d\n' % self.N_panel)
            f.write('\n\n')
        
        #
        f.write('oper\n')
        f.write('iter %d\n' % self.N_iter_max)
        f.write('M %f\n' % self.M)
        
        if self.visc:
            f.write('visc %f\n' % self.Re)
            
            f.write('vpar\n')
            f.write('N %f\n' % self.N)
            f.write('xtr %f %f\n' % (self.x_trip_up, self.x_trip_bot))
            f.write('vacc %f\n' % self.vacc)
            f.write('\n')
        
        f.write('pacc\n\n\n')
        f.write('alfa %f\n' % alpha)
        f.write('pwrt\n')
        
        outdir = self.run_dir+'/'+self.output_file
        f.write('%s\n' % self.output_file)
        if os.path.exists(outdir):
            f.write('y\n')
        
        f.write('\n')
        f.write('quit\n')
        
        f.close()
        
        if not os.path.exists(self.run_dir+'/xfoil.exe'):
            shutil.copyfile(self.xfoil_dir+'/xfoil.exe',
                            self.run_dir+'/xfoil.exe')
            
#        os.system('cd %s & xfoil < xfoil_run.dat' % self.run_dir)
        
        conv = False
        t0 = time.time()
        proc = subprocess.Popen('cd %s & xfoil < xfoil_run.dat' % self.run_dir,
                                shell=True)
        
        while(time.time()-t0<240):
            if proc.poll() is not None:
                conv = True
                break
            
        proc.kill()
        
        if conv:
            f = open(outdir, 'r')
            lines = f.readlines()
            f.close()
              
            data = lines[12].split()
            
            self.cl = float(data[1])
            self.cd = float(data[2])
            self.cdp = float(data[3])
            self.cm = float(data[4])
        else:
            self.cl, self.cd, self.cdp, self.cm = None, None, None, None
        
        return self.cl, self.cd, self.cdp, self.cm
    
    def create_polar(self, alpha_start, alpha_end, alpha_step):
        
        f = open(self.run_dir+'/xfoil_run.dat', 'w')
        
        f.write('plop\n')
        f.write('g,f\n')
        f.write('\n')
        
        #Load airfoil
        if self.airfoil is None:
            f.write('naca %s\n' % self.NACA)
        elif self.NACA is None:
            f.write('load %s.dat\n' % self.airfoil)
            f.write('pane\n')
#            f.write('ppar\n')
#            f.write('N %d\n' % self.N_panel)
#            f.write('\n\n')
        
        #
        f.write('oper\n')
        f.write('iter %d\n' % self.N_iter_max)
        f.write('M %f\n' % self.M)
        
        if self.visc:
            f.write('visc %f\n' % self.Re)
            
            f.write('vpar\n')
            f.write('N %f\n' % self.N)
            f.write('xtr %f %f\n' % (self.x_trip_up, self.x_trip_bot))
            f.write('vacc %f\n' % self.vacc)
            f.write('\n')
        
        f.write('pacc\n\n\n')
        f.write('aseq %f %f %f\n' % (alpha_start, alpha_end, alpha_step))
        f.write('pwrt\n')
        f.write('%s\n' % self.output_file)
        
        outdir = self.run_dir+'/'+self.output_file
        if os.path.exists(outdir):
            os.remove(outdir)
        
        f.write('\n')
        f.write('quit\n')
        
        f.close()
        
        if not os.path.exists(self.run_dir+'/xfoil.exe'):
            shutil.copyfile(self.xfoil_dir+'/xfoil.exe',
                            self.run_dir+'/xfoil.exe')
        
#        os.system('cd %s & xfoil < xfoil_run.dat' % self.run_dir)
        
        
        conv = False
        t0 = time.time()
        proc = subprocess.Popen('cd %s & xfoil < xfoil_run.dat' % self.run_dir,
                                shell=True)
        
        while(time.time()-t0<240):
            if proc.poll() is not None:
                conv = True
                break
            
        proc.kill()

        if conv:
            f = open(outdir)
            lines = f.readlines()
            f.close()
            
            cols = lines[10].split()
            
            
            result = []
            for line in lines[12:]:
                spl_line = line.split()
                spl_line = [float(i) for i in spl_line]
                result.append(spl_line)
                
            result = pd.DataFrame(np.array(result),
                                  columns=cols) 
        else:
            result = None
        
        return result
    