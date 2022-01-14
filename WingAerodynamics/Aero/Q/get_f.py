import scipy.io
import numpy as np
from scipy import interpolate as si
import pandas as pd



class Get_f(object):
    # read mat file
    def __init__(self, file):
        self.file = file
        mat = scipy.io.loadmat(self.file)  # taken  from PolarsDataFitted
        self.mat_data = mat['PolarsDataFitted']
        self.alpha_0_lst = []
        self.dCLdAlpha_lst = []
        self.dCd_dCl2_lst = []
        self.Re_lst = []
        self.sections = 10
        self.f_cl = {}
        self.f_cd = {}
        self.data = None
        self.sec = None


        # cl = np.zeros(10)
        # cd = np.zeros(10)

    def analyse(self):
        for i in range(self.sections):
            Re = self.mat_data[0][0][i][0][0][12][0][0]
            alpha_0 = self.mat_data[0][0][i][0][0][4][0]
            dCLdAlpha = self.mat_data[0][0][i][0][0][2][0]
            dCd_dCl2 = self.mat_data[0][0][i][0][0][8][0]

            self.Re_lst.append(Re)
            self.alpha_0_lst.append(alpha_0)
            self.dCLdAlpha_lst.append(dCLdAlpha)
            self.dCd_dCl2_lst.append(dCd_dCl2)

        for i in range(self.sections):
            self.sec = i

            self.aero_data = self.get_aero_data()

            self.f_cd[i] = self.get_f_cd()
            self.f_cl[i] = self.get_f_cl()

    def f(self, Re, alpha):
        Re_lst = self.Re_lst
        if Re in Re_lst:
            i1 = np.where(np.array(Re_lst)==Re)[0][0]
            i2 = i1

        else:
            i1 = min(np.where(np.array(Re_lst)>Re)[0])-1
            i2 = min(np.where(np.array(Re_lst) > Re)[0])

        Re1 = Re_lst[i1]
        alpha_01 = self.alpha_0_lst[i1][self.sec]
        d_alpha1 = alpha - alpha_01
        dCLdAlpha1 = self.dCLdAlpha_lst[i1][self.sec]
        dCd_dCl21 = self.dCd_dCl2_lst[i1][self.sec]
        cl1 = dCLdAlpha1 * np.deg2rad(d_alpha1)
        cd1 = dCd_dCl21*cl1**2

        if Re in Re_lst:
            cl = cl1
            cd = cd1

        else:

            Re2 = Re_lst[i2]
            alpha_02 = self.alpha_0_lst[i2][self.sec]
            d_alpha2 = alpha - alpha_02
            dCLdAlpha2 = self.dCLdAlpha_lst[i2][self.sec]
            dCd_dCl22 = self.dCd_dCl2_lst[i2][self.sec]
            cl2 = dCLdAlpha2 * np.deg2rad(d_alpha2)
            cd2 = dCd_dCl22 * cl2 ** 2
            f_cl = si.interp1d([Re1, Re2], [cl1, cl2])
            f_cd = si.interp1d([Re1, Re2], [cd1, cd2])

            cl = f_cl(Re)
            cd = f_cd(Re)

        return cl, cd


    def get_aero_data(self):
        # make dataframe
        aero_data = {}
        sec = self.sec
        Re_lst = self.Re_lst

        for Re in Re_lst:
            alpha_lst = []
            cl_lst = []
            cd_lst = []
            for alpha in (np.arange(-20,20,0.1)).round(2):
                cl, cd = self.f(Re, alpha)
                alpha_lst.append(alpha)
                cl_lst.append(cl)
                cd_lst.append(cd)

            df = {'alpha' :alpha_lst, 'cl': cl_lst,'cd': cd_lst}
            df = pd.DataFrame(df)
            df = df.set_index('alpha')
            aero_data[Re] = df

        return aero_data


    def get_f_cl(self):
        aero_data = self.aero_data
        # Get available Reynolds numbers
        Re_lst = list(aero_data.keys())
        Re_lst.sort()

        # Create one pandas.DataFrame with data for all Reynolds numbers
        for i, Re in enumerate(Re_lst):
            if i == 0:
                data = aero_data[Re]['cl'].rename(Re)
            else:
                data = pd.merge(data, aero_data[Re]['cl'].rename(Re),
                                how='outer',
                                left_index=True,
                                right_index=True)

        data = data.astype('float64')
        # Fill NaN values with data using interpolation/extrapolation
        data = data.interpolate(method='linear',
                                axis=0,
                                limit_area='inside')
        data = data.interpolate(method='linear',
                                axis=1,
                                limit_direction='both')
        self.data = data
        # Interpolate for angle of attack values
        alpha_lst = list(data.index)
        f_cl = si.interp2d(Re_lst, alpha_lst,
                           data.values,
                           kind='linear')

        return f_cl


    def get_f_cd(self):
        aero_data = self.aero_data
        # Get available Reynolds numbers
        Re_lst = list(aero_data.keys())
        Re_lst.sort()

        # Create one pandas.DataFrame with data for all Reynolds numbers
        for i, Re in enumerate(Re_lst):
            if i == 0:
                data = aero_data[Re]['cd'].rename(Re)
            else:
                data = pd.merge(data, aero_data[Re]['cd'].rename(Re),
                                how='outer',
                                left_index=True,
                                right_index=True)

        data = data.astype('float64')
        # Fill NaN values with data using interpolation/extrapolation
        data = data.interpolate(method='linear',
                                axis=0,
                                limit_area='inside')
        data = data.interpolate(method='linear',
                                axis=1,
                                limit_direction='both')

        # Interpolate for angle of attack values
        alpha_lst = list(data.index)
        f_cl = si.interp2d(Re_lst, alpha_lst,
                           data.values,
                           kind='linear')

        return f_cl