import numpy as np
import datetime as DT
import warnings
import glob, os, string
import csv


class cshoreIO():
    """This class takes care of CSHORE model input and output scripts"""
    def __init__(self):
        """
        This is the class initiation which initalizes the reading dictionaries
        """
        self.ODOC_dict = {}
        self.readIF_dict = {}
        self.OSETUP_dict = {}
        self.OXVELO_dict = {}
        self.OYVELO_dict = {}
        self.OCROSS_dict = {}
        self.OLONGS_dict = {}
        self.OBSUSL_dict = {}
        self.OVEG_dict = {}

    def read_CSHORE_ODOC(self, path):
        """This function will Read the ODOC input file

        Args:
          path: where is this ODOC txt file you want me to read.

        Returns:
            relevant variables in an ODOC dictionary.
        """
        with open(path + '/ODOC', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        # find header
        row_ind = np.asarray(np.argwhere(['OPTION ILINE' in s for s in tot])).flatten()
        if len(row_ind) == 0:
            self.ODOC_dict['header'] = tot
            self.ODOC_dict['run_success']= 0
        else:
            self.ODOC_dict['header'] = tot[0:row_ind[0]]

        # find IPROFL
        row_ind = np.asarray(np.argwhere(['OPTION IPROFL' in s for s in tot])).flatten()
        row = tot[row_ind[0]]
        self.ODOC_dict['iprofl'] = int(row[row.find('=') + 1:])

        # find ISEDAV
        row_ind = np.asarray(np.argwhere(['ISEDAV' in s for s in tot])).flatten()
        if len(row_ind) == 0:
            self.ODOC_dict['isedav'] = 0
        else:
            row = tot[row_ind[0]]
            self.ODOC_dict['isedav'] = int(row[row.find('=') + 1:row.find('=') + 3])

        # find IPERM
        row_ind = np.asarray(np.argwhere(['IMPERMEABLE' in s for s in tot])).flatten()
        if len(row_ind) == 0:
            self.ODOC_dict['iperm'] = True
        else:
            self.ODOC_dict['iperm'] = False

        # find NBINP
        row_ind = np.asarray(np.argwhere(['NBINP' in s for s in tot])).flatten()
        row = tot[row_ind[0]]
        self.ODOC_dict['nbinp'] = int(row[row.find('=') + 1:])

        # find GAMMA
        row_ind = np.asarray(np.argwhere(['Gamma' in s for s in tot])).flatten()
        row = tot[row_ind[0]]
        self.ODOC_dict['gamma'] = float(row[row.find('=') + 1:])

        # get longshore transport
        dum = tot[np.argwhere(['Transport Rate' in s for s in tot]).flatten()]
        ls_trans = np.zeros(len(dum)) * np.nan
        if len(dum) > 0:
            for ii in range(0, len(dum)):
                row = dum[ii]
                if len(row[row.find('=') + 1:].strip()) > 0:
                    ls_trans[ii] = float(row[row.find('=') + 1:])
                else:
                    ls_trans[ii] = np.nan
            self.ODOC_dict['longshore_transport'] = ls_trans

        # get wave conditions at SB
        row_ind = np.asarray(np.argwhere(['INPUT WAVE' in s for s in tot])).flatten()
        ind_start = row_ind[0] + 4
        row_ind = np.asarray(np.argwhere(['INPUT BEACH AND STRUCTURE' in s for s in tot])).flatten()
        ind_end = row_ind[0] - 1
        wave_cond = tot[ind_start: ind_end]
        time_offshore = np.zeros(len(wave_cond)) * np.nan
        Tp_bc = np.zeros(len(wave_cond)) * np.nan
        Hs_bc = np.zeros(len(wave_cond)) * np.nan
        setup_bc = np.zeros(len(wave_cond)) * np.nan
        wl_bc = np.zeros(len(wave_cond)) * np.nan
        angle_bc = np.zeros(len(wave_cond)) * np.nan
        cnt = 0
        for line in wave_cond:
            dum = line.split()
            time_offshore[cnt] = float(dum[0])
            Tp_bc[cnt] = float(dum[1])
            Hs_bc[cnt] = np.sqrt(2) * float(dum[2])
            setup_bc[cnt] = float(dum[3])
            wl_bc[cnt] = float(dum[4])
            angle_bc[cnt] = float(dum[5])
            cnt = cnt + 1
        # THIS WILL ONLY RETURN THE FIRST AND LAST 10 BC POINTS!!!!
        self.ODOC_dict['time_offshore'] = time_offshore
        self.ODOC_dict['Tp_bc'] = Tp_bc
        self.ODOC_dict['Hs_bc'] = Hs_bc
        self.ODOC_dict['setup_bc'] = setup_bc
        self.ODOC_dict['wl_bc'] = wl_bc
        self.ODOC_dict['angle_bc'] = angle_bc

        # find runup
        dum2p = tot[np.argwhere(['2 percent runup' in s for s in tot]).flatten()]
        dummean = tot[np.argwhere(['Mean runup' in s for s in tot]).flatten()]
        runup_2_percent = np.zeros(len(dum2p)) * np.nan
        runup_mean = np.zeros(len(dum2p)) * np.nan
        if len(dum2p) > 0:
            for ii in range(0, len(dum2p)):
                row1 = dum2p[ii]
                row2 = dummean[ii]
                if len(row1[row1.find('R2P=') + 4:].strip()) > 0:
                    runup_2_percent[ii] = float(row1[row1.find('R2P=') + 4:])
                    runup_mean[ii] = float(row2[row1.find('R2P=') + 4:])
                else:
                    runup_2_percent[ii] = np.nan
                    runup_mean[ii] = np.nan
            self.ODOC_dict['runup_2_percent'] = runup_2_percent
            self.ODOC_dict['runup_mean'] = runup_mean

        # find jdry
        dum = tot[np.argwhere(['JDRY' in s for s in tot]).flatten()]
        jdry = np.zeros(len(dum)) * np.nan
        if len(dum) > 0:
            for ii in range(0, len(dum)):
                row = dum[ii]
                if len(row[row.find('JDRY=') + 5:].strip()) > 0:
                    jdry[ii] = float(row[row.find('JDRY=') + 5:])
                else:
                    jdry[ii] = np.nan
            self.ODOC_dict['jdry'] = jdry

        # find SWL at sea boundary
        dum = tot[np.argwhere([' SWL=' in s for s in tot]).flatten()]
        swl = np.zeros(len(dum)) * np.nan
        if len(dum) > 0:
            for ii in range(0, len(dum)):
                row = dum[ii]
                if len(row[row.find(' SWL=') + 5:].strip()) > 0:
                    swl[ii] = float(row[row.find(' SWL=') + 5:])
                else:
                    swl[ii] = np.nan
            self.ODOC_dict['swl'] = swl

        # find node number of SWL
        dum = tot[np.argwhere([' JSWL=' in s for s in tot]).flatten()]
        jswl = np.zeros(len(dum)) * np.nan
        if len(dum) > 0:
            for ii in range(0, len(dum)):
                row = dum[ii]
                if len(row[row.find(' JSWL=') + 6:].strip()) > 0:
                    jswl[ii] = float(row[row.find(' JSWL=') + 6:])
                else:
                    jswl[ii] = np.nan
            self.ODOC_dict['jswl'] = jswl

        # find jr
        dum = tot[np.argwhere(['JR=' in s for s in tot]).flatten()]
        jr = np.zeros(len(dum)) * np.nan
        if len(dum) > 0:
            for ii in range(0, len(dum)):
                row = dum[ii]
                if len(row[row.find('JR=') + 3:].strip()) > 0:
                    jr[ii] = float(row[row.find('JR=') + 3:])
                else:
                    jr[ii] = np.nan
            self.ODOC_dict['jr'] = jr

        # swash zone bottom slope
        row_ind = np.asarray(np.argwhere(['Swash zone bottom slope' in s for s in tot])).flatten()
        dum_slp = tot[row_ind]
        dum_x1 = tot[row_ind + 1]
        dum_x2 = tot[row_ind + 2]
        dum_z1 = tot[row_ind + 3]
        dum_z2 = tot[row_ind + 4]
        if len(dum_slp) > 0:
            slp = np.zeros(len(dum_slp)) * np.nan
            x1 = np.zeros(len(dum_x1)) * np.nan
            x2 = np.zeros(len(dum_x2)) * np.nan
            z1 = np.zeros(len(dum_z1)) * np.nan
            z2 = np.zeros(len(dum_z2)) * np.nan
            cnt = 0
            for ii in range(0, len(dum_slp)):
                if len(dum_slp[ii][dum_slp[ii].find('=') + 1:].strip()) > 0:
                    slp[cnt] = float(dum_slp[ii][dum_slp[ii].find('=') + 1:])
                    x1[cnt] = float(dum_x1[ii][dum_slp[ii].find('=') + 1:])
                    x2[cnt] = float(dum_x2[ii][dum_slp[ii].find('=') + 1:])
                    z1[cnt] = float(dum_z1[ii][dum_slp[ii].find('=') + 1:])
                    z2[cnt] = float(dum_z2[ii][dum_slp[ii].find('=') + 1:])
                else:
                    slp[cnt] = np.nan
                    x1[cnt] = np.nan
                    x2[cnt] = np.nan
                    z1[cnt] = np.nan
                    z2[cnt] = np.nan
                cnt = cnt + 1
            self.ODOC_dict['slprun'] = slp
            self.ODOC_dict['x1run'] = x1
            self.ODOC_dict['x2run'] = x2
            self.ODOC_dict['z1run'] = z1
            self.ODOC_dict['z2run'] = z2


    def read_CSHORE_infile(self, path):
        """this function will read the CSHORE infile

        Args:
          path: where is this infile txt file you want me to read.

        Returns:
            relevant variables in an INFILE dictionary.
        """

        with open(path + '/infile', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        # find IOVER
        row_ind = np.asarray(np.argwhere(['IOVER' in s for s in tot])).flatten()
        row = tot[row_ind[0]]
        self.readIF_dict['iover'] = float(row[0:row.find('-') - 1])

        # find IVEG
        row_ind = np.asarray(np.argwhere(['IVEG' in s for s in tot])).flatten()
        row = tot[row_ind[0]]
        self.readIF_dict['iveg'] = float(row[0:row.find('-') - 1])

        # find effB, effF, and blp
        if self.ODOC_dict['iprofl'] == 1:
            # EFFB and EFFF
            row_ind = np.asarray(np.argwhere(['EFFB' in s for s in tot])).flatten()
            row = tot[row_ind[0]]
            dum = row[0:row.find('-') - 1].split()
            self.readIF_dict['effB'] = float(dum[0])
            self.readIF_dict['effF'] = float(dum[1])
            # BLP and TANPHI
            row_ind = np.asarray(np.argwhere(['BLP' in s for s in tot])).flatten()
            row = tot[row_ind[0]]
            dum = row[0:row.find('-') - 1].split()
            self.readIF_dict['tanphi'] = float(dum[0])
            self.readIF_dict['blp'] = float(dum[1])
            # ILAB
            row_ind = np.asarray(np.argwhere(['ILAB' in s for s in tot])).flatten()
            row = tot[row_ind[0]]
            row[0:row.find('-') - 1]
            self.readIF_dict['ilab'] = float(row[0:row.find('-') - 1])

        # find vegetation extent
        if self.readIF_dict['iveg'] == 1:
            row_ind = np.asarray(np.argwhere(['VEGCD' in s for s in tot])).flatten()
            dum = tot[row_ind[0] + 1:row_ind[0] + int(self.ODOC_dict['nbinp']) + 1]
            veg_n = np.zeros(len(dum))
            veg_dia = np.zeros(len(dum))
            veg_ht = np.zeros(len(dum))
            veg_rod = np.zeros(len(dum))
            cnt = 0
            for item in dum:
                row = item.split()
                veg_n[cnt] = float(row[0])
                veg_dia[cnt] = float(row[1])
                veg_ht[cnt] = float(row[2])
                veg_rod[cnt] = float(row[3])
                cnt = cnt + 1
            self.readIF_dict['veg_n'] = veg_n
            self.readIF_dict['veg_dia'] = veg_dia
            self.readIF_dict['veg_ht'] = veg_ht
            self.readIF_dict['veg_rod'] = veg_rod

        #now we are going to read the BC stuff from the infile, because Brad's ODOC bc stuff ONLY HAS THE FIRST AND LAST 10!?!?! WTFWTFWTF
        # find NSURGE
        row_ind_1 = np.asarray(np.argwhere(['NSURGE' in s for s in tot])).flatten()
        row_ind_2 = np.asarray(np.argwhere(['NBINP' in s for s in tot])).flatten()
        rows = tot[row_ind_1[0]+1:row_ind_2[0]]

        # here is where i have to deal with the ilab toggle
        # get me the length of all the rows
        row_len = np.array([len(s.split()) for s in rows])
        if 6 not in row_len:
            # this means ilab was set to zero
            #sort by length
            row_swlbc = rows[row_len == 2]
            row_waves = rows[row_len == 4]
            num_pts = len(row_swlbc)

            time_offshore = np.zeros(num_pts)
            Tp_bc = np.zeros(num_pts)
            Hs_bc = np.zeros(num_pts)
            setup_bc = np.zeros(num_pts) # the model must just assume this is zero, because this information IS NOT contained in the infile if ilab == 0
            wl_bc = np.zeros(num_pts)
            angle_bc = np.zeros(num_pts)

            for ii in range(0, len(row_waves)):
                dum_swlbc = row_swlbc[ii].split()
                wl_bc[ii] = dum_swlbc[1]
                dum_waves = row_waves[ii].split()
                time_offshore[ii] = dum_waves[0]
                Tp_bc[ii] = dum_waves[1]
                Hs_bc[ii] = str(np.sqrt(2) * float(dum_waves[2]))
                angle_bc[ii] = dum_waves[3]

            self.readIF_dict['time_offshore'] = time_offshore
            self.readIF_dict['Tp_bc'] = Tp_bc
            self.readIF_dict['Hs_bc'] = Hs_bc
            self.readIF_dict['setup_bc'] = setup_bc
            self.readIF_dict['wl_bc'] = wl_bc
            self.readIF_dict['angle_bc'] = angle_bc
        else:
            # this means ilab was set to 1
            num_pts = len(rows)

            time_offshore = np.zeros(num_pts)
            Tp_bc = np.zeros(num_pts)
            Hs_bc = np.zeros(num_pts)
            setup_bc = np.zeros(num_pts)
            wl_bc = np.zeros(num_pts)
            angle_bc = np.zeros(num_pts)

            for ii in range(0, len(rows)):
                dum = rows[ii].split()
                time_offshore[ii] = dum[0]
                Tp_bc[ii] = dum[1]
                Hs_bc[ii] = str(np.sqrt(2) * float(dum[2]))
                setup_bc[ii] = dum[3]
                wl_bc[ii] = dum[4]
                angle_bc[ii] = dum[5]

            self.readIF_dict['time_offshore'] = time_offshore
            self.readIF_dict['Tp_bc'] = Tp_bc
            self.readIF_dict['Hs_bc'] = Hs_bc
            self.readIF_dict['setup_bc'] = setup_bc
            self.readIF_dict['wl_bc'] = wl_bc
            self.readIF_dict['angle_bc'] = angle_bc

    def read_CSHORE_OSETUP(self, path):
        """This function will Read the OSETUP input file

        Args:
          path: where is this OSETUP txt file you want me to read.

        Returns:
            relevant variables in an OSETUP dictionary.

        """

        with open(path + '/OSETUP', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        ind_track = 0
        for ii in range(0, len(self.readIF_dict['time_offshore'])):
            self.OSETUP_dict['hydro%s' % str(ii + 1)] = {}
            row1 = tot[ind_track]

            if float(row1.split()[0]) == 1:
                N = int(row1.split()[1])
                tme = float(row1.split()[-1])
            else:
                N = int(row1.split()[0])

            dum = tot[ind_track + 1:ind_track + 1 + N + 1]
            x = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            setup = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            depth = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            sigma = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            Hs = np.zeros(self.ODOC_dict['nbinp']) * np.nan  # note: we are feeding it Hmo!!!!!!!
            for ss in range(0, N):
                x[ss] = float(dum[ss].split()[0])
                setup[ss] = float(dum[ss].split()[1])
                depth[ss] = float(dum[ss].split()[2])
                sigma[ss] = float(dum[ss].split()[3])
                Hs[ss] = np.sqrt(8) * np.sqrt(2) * sigma[ss]

            ind_track = ind_track + N + 1

            self.OSETUP_dict['hydro%s' % str(ii + 1)]['x'] = x
            self.OSETUP_dict['hydro%s' % str(ii + 1)]['setup'] = setup
            self.OSETUP_dict['hydro%s' % str(ii + 1)]['depth'] = depth
            self.OSETUP_dict['hydro%s' % str(ii + 1)]['sigma'] = sigma
            self.OSETUP_dict['hydro%s' % str(ii + 1)]['Hs'] = Hs
            self.OSETUP_dict['hydro%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OXVELO(self, path):
        """This function will Read the OXVELO input file

        Args:
          path: where is this OXVELO txt file you want me to read.

        Returns:
            relevant variables in an OXVELO dictionary.
        """

        with open(path + '/OXVELO', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        ind_track = 0
        for ii in range(0, len(self.readIF_dict['time_offshore'])):
            self.OXVELO_dict['hydro%s' % str(ii + 1)] = {}
            row1 = tot[ind_track]

            if float(row1.split()[0]) == 1:
                N = int(row1.split()[1])
                tme = float(row1.split()[-1])
            else:
                N = int(row1.split()[0])

            dum = tot[ind_track + 1:ind_track + 1 + N + 1]
            x_xvelo = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            umean = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            ustd = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            for ss in range(0, N):
                x_xvelo[ss] = float(dum[ss].split()[0])
                umean[ss] = float(dum[ss].split()[1])
                ustd[ss] = float(dum[ss].split()[2])

            ind_track = ind_track + N + 1

            self.OXVELO_dict['hydro%s' % str(ii + 1)]['x_xvelo'] = x_xvelo
            self.OXVELO_dict['hydro%s' % str(ii + 1)]['umean'] = umean
            self.OXVELO_dict['hydro%s' % str(ii + 1)]['ustd'] = ustd
            self.OXVELO_dict['hydro%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OYVELO(self, path):
        """This function will Read the OYVELO input file

        Args:
          path: where is this OYVELO txt file you want me to read.

        Returns:
            relevant variables in an OYVELO dictionary.

        """
        # OYVELO reading function
        with open(path + '/OYVELO', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        tot_len = np.asarray(list(map(lambda s: len(s.split()), tot)))
        header_ind = np.argwhere(tot_len != 4).flatten()  # it will have 4 things in the row if it actually has data

        for ii in range(0, len(header_ind)):
            self.OYVELO_dict['hydro%s' % str(ii + 1)] = {}

            # get all the data in between the header inds!
            if ii == len(header_ind) - 1:
                dum = tot[header_ind[ii]:]
            else:
                dum = tot[header_ind[ii]:header_ind[ii + 1]]

            if float(dum[0].split()[0]) == 1:
                N = int(dum[0].split()[1])
                tme = float(dum[0].split()[-1])
            else:
                N = int(dum[0].split()[0])

            # if the length of dum is <= 1, that means it ONLY has a header (i.e., no actual data)
            if len(dum) <= 1:

                x_yvelo = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                stheta = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                vmean = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                vstd = np.zeros(self.ODOC_dict['nbinp']) * np.nan

            else:
                x_yvelo = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                stheta = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                vmean = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                vstd = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                for ss in range(1, N+1):
                    x_yvelo[ss - 1] = float(dum[ss].split()[0])
                    stheta[ss - 1] = float(dum[ss].split()[1])
                    vmean[ss - 1] = float(dum[ss].split()[2])
                    vstd[ss - 1] = float(dum[ss].split()[3])

            self.OYVELO_dict['hydro%s' % str(ii + 1)]['x_yvelo'] = x_yvelo
            self.OYVELO_dict['hydro%s' % str(ii + 1)]['stheta'] = stheta
            self.OYVELO_dict['hydro%s' % str(ii + 1)]['vmean'] = vmean
            self.OYVELO_dict['hydro%s' % str(ii + 1)]['vstd'] = vstd
            self.OYVELO_dict['hydro%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OCROSS(self, path):
        """This function will Read the OCROSS input file

        Args:
          path: where is this OCROSS txt file you want me to read.

        Returns:
            relevant variables in an OCROSS dictionary.

        """

        if self.ODOC_dict['iprofl'] == 1:

            # OCROSS reading function
            with open(path + '/OCROSS', 'r') as fid:
                tot = fid.readlines()
            fid.close()
            tot = np.asarray(map(lambda s: s.strip(), tot))

            tot_len = np.asarray(map(lambda s: len(s.split()), tot))
            header_ind = np.argwhere(tot_len != 4).flatten()  # it will have 4 things in the row if it actually has data

            for ii in range(0, len(header_ind)):
                self.OCROSS_dict['sed%s' % str(ii + 1)] = {}

                # get all the data in between the header inds!
                if ii == len(header_ind) - 1:
                    dum = tot[header_ind[ii]:]
                else:
                    dum = tot[header_ind[ii]:header_ind[ii + 1]]

                if float(dum[0].split()[0]) == 1:
                    N = int(dum[0].split()[1])
                    tme = float(dum[0].split()[-1])
                else:
                    N = int(dum[0].split()[0])

                # if the length of dum is <= 1, that means it ONLY has a header (i.e., no actual data)
                if len(dum) <= 1:

                    x_cross = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qbx = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qsx = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qx = np.zeros(self.ODOC_dict['nbinp']) * np.nan

                else:
                    x_cross = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qbx = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qsx = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qx = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    for ss in range(1, N + 1):
                        x_cross[ss - 1] = float(dum[ss].split()[0])
                        qbx[ss - 1] = float(dum[ss].split()[1])
                        qsx[ss - 1] = float(dum[ss].split()[2])
                        qx[ss - 1] = float(dum[ss].split()[3])

                self.OCROSS_dict['sed%s' % str(ii + 1)]['x_cross'] = x_cross
                self.OCROSS_dict['sed%s' % str(ii + 1)]['qbx'] = qbx
                self.OCROSS_dict['sed%s' % str(ii + 1)]['qsx'] = qsx
                self.OCROSS_dict['sed%s' % str(ii + 1)]['qx'] = qx
                self.OCROSS_dict['sed%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OLONGS(self, path):
        """This function will Read the OLONGS input file

        Args:
          path: where is this OLONGS txt file you want me to read.

        Returns:
            relevant variables in an OLONGS dictionary.

        """
        # OLONGS reading function - this is currently permanently turned off in Brad's script...
        if self.ODOC_dict['iprofl'] == 1 and False:

            with open(path + '/OLONGS', 'r') as fid:
                tot = fid.readlines()
            fid.close()
            tot = np.asarray(map(lambda s: s.strip(), tot))

            tot_len = np.asarray(map(lambda s: len(s.split()), tot))
            header_ind = np.argwhere(tot_len != 4).flatten()  # it will have 4 things in the row if it actually has data

            for ii in range(0, len(header_ind)):
                self.OLONGS_dict['sed%s' % str(ii + 1)] = {}

                # get all the data in between the header inds!
                if ii == len(header_ind) - 1:
                    dum = tot[header_ind[ii]:]
                else:
                    dum = tot[header_ind[ii]:header_ind[ii + 1]]

                if float(dum[0].split()[0]) == 1:
                    # N = int(dum[0].split()[1])
                    N = min(np.size(tot) - 1, int(dum[0].split()[1]))
                    tme = float(dum[0].split()[-1])
                else:
                    # N = int(dum[0].split()[0])
                    N = min(np.size(tot) - 1, int(dum[0].split()[0]))
                # if the length of dum is <= 1, that means it ONLY has a header (i.e., no actual data)
                if len(dum) <= 1:

                    x_long = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qby = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qsy = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qy = np.zeros(self.ODOC_dict['nbinp']) * np.nan

                else:
                    x_long = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qby = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qsy = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    qy = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    for ss in range(1, N + 1):
                        x_long[ss - 1] = float(dum[ss].split()[0])
                        qby[ss - 1] = float(dum[ss].split()[1])
                        qsy[ss - 1] = float(dum[ss].split()[2])
                        qy[ss - 1] = float(dum[ss].split()[3])

                self.OLONGS_dict['sed%s' % str(ii + 1)]['x_long'] = x_long
                self.OLONGS_dict['sed%s' % str(ii + 1)]['qby'] = qby
                self.OLONGS_dict['sed%s' % str(ii + 1)]['qsy'] = qsy
                self.OLONGS_dict['sed%s' % str(ii + 1)]['qy'] = qy
                self.OLONGS_dict['sed%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OBSUSL(self, path):
        """This function will Read the OBSUSL input file

        Args:
          path: where is this OBSUSL txt file you want to read.

        Returns:
            relevant variables in an OBSUSL dictionary.
        """
        # OBSUSL reading function
        if self.ODOC_dict['iprofl'] == 1:

            with open(path + '/OBSUSL', 'r') as fid:
                tot = fid.readlines()
            fid.close()
            tot = np.asarray(map(lambda s: s.strip(), tot))

            tot_len = np.asarray(map(lambda s: len(s.split()), tot))
            header_ind = np.argwhere(tot_len != 4).flatten()  # it will have 4 things in the row if it actually has data

            for ii in range(0, len(header_ind)):
                self.OBSUSL_dict['sed%s' % str(ii + 1)] = {}

                # get all the data in between the header inds!
                if ii == len(header_ind) - 1:
                    dum = tot[header_ind[ii]:]
                else:
                    dum = tot[header_ind[ii]:header_ind[ii + 1]]

                if float(dum[0].split()[0]) == 1:
                    N = int(dum[0].split()[1])
                    tme = float(dum[0].split()[-1])
                else:
                    N = int(dum[0].split()[0])

                # if the length of dum is <= 1, that means it ONLY has a header (i.e., no actual data)
                if len(dum) <= 1:

                    x_susl = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    ps = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    pb = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    vs = np.zeros(self.ODOC_dict['nbinp']) * np.nan

                else:
                    x_susl = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    ps = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    pb = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    vs = np.zeros(self.ODOC_dict['nbinp']) * np.nan
                    for ss in range(1, N + 1):
                        x_susl[ss - 1] = float(dum[ss].split()[0])
                        ps[ss - 1] = float(dum[ss].split()[1])
                        pb[ss - 1] = float(dum[ss].split()[2])
                        vs[ss - 1] = float(dum[ss].split()[3])

                self.OBSUSL_dict['sed%s' % str(ii + 1)]['x_susl'] = x_susl
                self.OBSUSL_dict['sed%s' % str(ii + 1)]['ps'] = ps
                self.OBSUSL_dict['sed%s' % str(ii + 1)]['pb'] = pb
                self.OBSUSL_dict['sed%s' % str(ii + 1)]['vs'] = vs
                self.OBSUSL_dict['sed%s' % str(ii + 1)]['time_end'] = tme


    def read_CSHORE_OVEG(self, path):
        """This function will Read the OVEG input file

        Args:
          path: where is this OVEG txt file you want me to read.

        Returns:
            relevant variables in an OVEG dictionary.

        """

        with open(path + '/OVEG', 'r') as fid:
            tot = fid.readlines()
        fid.close()
        tot = np.asarray(list(map(lambda s: s.strip(), tot)))

        ind_track = 0
        for ii in range(0, len(self.readIF_dict['time_offshore'])):
            self.OVEG_dict['veg%s' % str(ii + 1)] = {}
            row1 = tot[ind_track]

            if float(row1.split()[0]) == 1:
                # N = int(row1.split()[1])
                N = min(np.size(tot) - 1, int(row1.split()[1]))
                tme = float(row1.split()[-1])
            else:
                # N = int(row1.split()[0])
                N = min(np.size(tot) - 1, int(row1.split()[0]))

            dum = tot[ind_track + 1:ind_track + 1 + N + 1]
            x = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            Cd = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            ESH = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            EBH = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            EPH = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            Hstem = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            Hbld = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            for ss in range(0, N):
                x[ss] = float(dum[ss].split()[0])
                Cd[ss] = float(dum[ss].split()[1])
                ESH[ss] = float(dum[ss].split()[2])
                EBH[ss] = float(dum[ss].split()[3])
                EPH[ss] = float(dum[ss].split()[4])
                Hstem[ss] = float(dum[ss].split()[5])
                Hbld[ss] = float(dum[ss].split()[6])

            ind_track = ind_track + N + 1

            self.OVEG_dict['veg%s' % str(ii + 1)]['x'] = x
            self.OVEG_dict['veg%s' % str(ii + 1)]['Cd'] = Cd
            self.OVEG_dict['veg%s' % str(ii + 1)]['ESH'] = ESH
            self.OVEG_dict['veg%s' % str(ii + 1)]['EBH'] = EBH
            self.OVEG_dict['veg%s' % str(ii + 1)]['EPH'] = EPH
            self.OVEG_dict['veg%s' % str(ii + 1)]['Hstem'] = Hstem
            self.OVEG_dict['veg%s' % str(ii + 1)]['Hbld'] = Hbld
            self.OVEG_dict['veg%s' % str(ii + 1)]['time_end'] = tme

    def load_CSHORE_results(self, path):
        """This function will take the output dictionaries from all these previous read functions and stack all the information
        into 7 output dictionaries selected based on the catagory of information based in each dictionart

        Args:
          path: where are all these .txt files i am supposed to read

        Returns:
            params (dict):
                nbinp: -number of nodes in the bathy profile

                tanphi: - tangent (sediment friction angle - units?)

                iover: - 0 = no overtopping , 1 = include overtopping

                iprofl: - toggle to run morphology (0 = no morph, 1 = morph)

                effB: - suspension efficiency due to breaking eB (what is this?)

                ilab: - controls the boundary condition timing. Don't change.  BRAD SAID TO CHANGE THIS TO 1!!!!!!

                blp: - bedload parameter (what is this?)

                num_steps: - number of time-steps in the model

                header: - string that makes up the header of the cshore infile.txt

                iperm: - 0 = no permeability, 1 = permeable

                effF: - suspension efficiency due to friction ef (what is this?)

                isedav: - 0 = unlimited sand, 1 = hard bottom

                gamma: - shallow water ratio of wave height to water depth (wave breaking parameter?)

                iveg: - 1 = include vegitation effect (I'm assuming)

            bc (dict):
                angle_offshore: - wave angle at the offshore boundary at each time-step

                wave_setup_offshore: - wave setup at the offshore boundary at each time-step

                Hs_offshore: - wave height at the offshore boundary at each time-step

                time_offshore: - time of each model time-step in seconds elapsed since model start time

                Tp_offshore: - wave period at the offshore boundary at each time-step

                strm_tide_offshore: - water level at the offshore boundary at each time-step

            veg(dict):
                rod: - vegitation erosion limit below sand for failure (what is this?!?!?)

                dia: - vegetation dia. (units m)

                ht: - vegetation height (units m)

                n: - vegetation density (units stems/m2)

                EPH: - effective vegetation height (m)

                Cd: - drag coefficient (unit 1)


            hydro(dict):
                stheta: - sine of the wave angle (theta) for each time step

                jdry: - the most landward node in the wet/dry region

                time_end: - i think this is the time at the end of each model time step in seconds after model start

                x_xvelo: - the tech report we got from Brad does not mention this variable in its discussion of OXVELO

                umean: - mean cross-shore velocity at each node and time-step of the simulation

                jswl: -index of the node at the still water shoreline

                swl: - I am going to guess based on Brad's Tech Note that this is the still water level?

                x2run: - i think this is just the x-position of each node for the second of two bottom profiles (if you gave it two bottom profiles)?  This guess is just based on Brad's Tech report

                Hs: - wave height at each node and time-step of the simulation

                slprun: - representative bottom slope in the swash zone

                x1run: -i think this is just the x-position of each node for the first of two bottom profiles (if you gave it two bottom profiles)?  This guess is just based on Brad's Tech report

                vstd: - standard deviation of the alongshore velocity at each node and time-step of the simulation

                mwl: - mean water level at each node and time-step

                vmean: - mean alongshore velocity at each node and time-step

                z1run: -i think this is just the elevation of each node for the first of two bottom profiles (if you gave it two bottom profiles)?  This guess is just based on Brad's Tech report

                x_yvelo: -the tech report we got from Brad does not mention this variable in its discussion of OYVELO or OXVELO

                runup_2_percent: - 2 percent exceedance run-up elevation at each time-step

                jr: - i think this is just the number of nodes in your cross-section bottom boundary?

                runup_mean: - mean runup elevation at each time-step

                ustd: -standard deviation of the cross-hore velocity at each node and time-step of the simulation

                depth: - water depth at each node and time-step

                x: - x-position of each node in model coordinates

                sigma: - standard deviation of the water surface elevation for each time-step in the model run

                z2run: -i think this is just the elevation of each node for the second of two bottom profiles (if you gave it two bottom profiles) - needs to be checked

            sed (dict):
                ps: - probability of sediment suspension

                x_susl: - this variable is not mentioned in Brad's Tech Report.  if i had to guess i would say it is just the x-position of each node

                qbx: - bed load sediment transport in cross-shore direction

                qx: - total sediment transport in cross-shore direction

                vs: - suspended sediment volume per unit horizontal bottom area

                time_end: - i think this is the time at the end of each model time step in seconds after model start

                x_cross: -this variable is not mentioned in Brad's Tech Report.  if i had to guess i would say it is just the x-position of each node

                pb: - probability of sediment movement

                qy: - total sediment transport in the alongshore direction

                qby: - bed load sediment transport in the alongshore direction

                x_long: - this variable is not mentioned in Brad's Tech Report.  if i had to guess i would say it is just the x-position of each node

                qsy: - suspended load sediment transport in the alongshore direction

                qsx: - suspended load sediment transport in the cross-shore direction

            morpho (dict):
                ivegetated: - which nodes are vegetated or not

                x: - x-position of each node in the model

                zb_p: - vertical coordinate of the lower boundary of the permeable layer

                zb: - bottom elevation of each node at each time-step of the simulation

                time: - i think this is the time at the end of each model time step in seconds after model start

            meta (dict):
                bathy_survey_num: - survey number that the initial bathymetry was taken from

                bathy_prof_num: - profileLine number that makes up the initial cshore bathy

                blank_wave_data: - which wave data were interpolated to generate the input file

                BC_gage: - name of the wave gage used as the BC for this cshore run

                fric_fac: - this is the friction factor that is set automatically in the infile generation script

                BC_FRF_Y: - this is the yFRF location of the BC wave gage.

                BC_FRF_X: - this is the cFRF location of the BC wave gage.  I think this is rounded down to the nearest 1 m

                startTime: - start time for this particular cshore model run

                bathy_surv_stime: - this is the start time of the bathymetry survey that was used in the model

                bathy_y_sdev: - this is the standard deviation of the actualy yFRF positions of each point in the survey used in the model

                version: - MOBILE, MOBILE_RESET or FIXED right now

                timerun: - duration that each individual simulation was run.  24 hours unless specified otherwise

                dx: - grid spacing in the model. default is 1 m unless specified otherwise

                time_step: - model time-step.  1 hour unless specified otherwise

        """

        self.read_CSHORE_ODOC(path)
        self.read_CSHORE_infile(path)
        self.read_CSHORE_OSETUP(path)
        self.read_CSHORE_OXVELO(path)
        self.read_CSHORE_OYVELO(path)
        self.read_CSHORE_OCROSS(path)
        self.read_CSHORE_OLONGS(path)
        self.read_CSHORE_OBSUSL(path)
        self.read_CSHORE_OVEG(path)

        # all this does is take the dicts I created from my output text files a reformats them
        # params

        params = {}
        bc = {}
        hydro = {}
        veg = {}
        morpho = {}
        sed = {}

        params['header'] = self.ODOC_dict['header']
        params['iprofl'] = self.ODOC_dict['iprofl']
        params['isedav'] = self.ODOC_dict['isedav']
        params['iperm'] = self.ODOC_dict['iperm']
        params['nbinp'] = self.ODOC_dict['nbinp']
        params['gamma'] = self.ODOC_dict['gamma']
        params['iover'] = self.readIF_dict['iover']
        params['iveg'] = self.readIF_dict['iveg']

        num_steps = len(self.OSETUP_dict.keys())
        params['num_steps'] = num_steps

        # veg
        if params['iveg'] == 0:
            veg['x'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['Cd'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['ESH'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['EBH'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['EPH'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['Hstem'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            veg['Hbld'] = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            # veg['x'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['Cd'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['ESH'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['EBH'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['EPH'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['Hstem'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
            # veg['Hbld'] = np.zeros(self.ODOC_dict['nbinp']) * np.nan
        else: # iveg==3
            # from OVEG
            time_end = np.zeros([num_steps, 1]) * np.nan
            x = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            Cd = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            ESH = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            EBH = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            EPH = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            Hstem = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
            Hbld = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan

            for ii in range(0, num_steps):
                # OVEG
                temp_dict_veg = self.OVEG_dict['veg%s' % str(ii + 1)]
                time_end[ii] = temp_dict_veg['time_end']
                x[ii] = temp_dict_veg['x']
                Cd[ii] = temp_dict_veg['Cd']
                ESH[ii] = temp_dict_veg['ESH']
                EBH[ii] = temp_dict_veg['EBH']
                EPH[ii] = temp_dict_veg['EPH']
                Hstem[ii] = temp_dict_veg['Hstem']
                Hbld[ii] = temp_dict_veg['Hbld']

            veg['time_end'] = time_end
            veg['x'] = x
            veg['Cd'] = Cd
            veg['ESH'] = ESH
            veg['EBH'] = EBH
            veg['EPH'] = EPH
            veg['Hstem'] = Hstem
            veg['Hbld'] = Hbld

        # hydro
        # from ODOC
        if 'jdry' in self.ODOC_dict.keys():
            hydro['jdry'] = self.ODOC_dict['jdry'].reshape(self.ODOC_dict['jdry'].shape[0], -1)

        if 'swl' in self.ODOC_dict.keys():
            hydro['swl'] = self.ODOC_dict['swl'].reshape(self.ODOC_dict['swl'].shape[0], -1)

        if 'runup_2_percent' in self.ODOC_dict.keys():
            hydro['runup_2_percent'] = self.ODOC_dict['runup_2_percent'].reshape(self.ODOC_dict['runup_2_percent'].shape[0], -1)
            hydro['runup_mean'] = self.ODOC_dict['runup_mean'].reshape(self.ODOC_dict['runup_mean'].shape[0], -1)
            # in CSHORE, hydro['runup_2_percent'] is the sum of R2% and SWL. This is the total wave runup. Engineers often use the total runup.
            # lzhu subtracts SWL from the total runup in order to get the wave-related part of the runup.
            for ii in range(params['num_steps']):
                hydro['runup_2_percent'][ii][0] = hydro['runup_2_percent'][ii][0] - hydro['swl'][ii]
                hydro['runup_mean'][ii][0] = hydro['runup_mean'][ii][0] - hydro['swl'][ii]

        if 'jswl' in self.ODOC_dict.keys():
            hydro['jswl'] = self.ODOC_dict['jswl'].reshape(self.ODOC_dict['jswl'].shape[0], -1)

        if 'jr' in self.ODOC_dict.keys():
            hydro['jr'] = self.ODOC_dict['jr'].reshape(self.ODOC_dict['jr'].shape[0], -1)

        if 'slprun' in self.ODOC_dict.keys():
            hydro['slprun'] = self.ODOC_dict['slprun'].reshape(self.ODOC_dict['jr'].shape[0], -1)
            hydro['x1run'] = self.ODOC_dict['x1run'].reshape(self.ODOC_dict['jr'].shape[0], -1)
            hydro['x2run'] = self.ODOC_dict['x2run'].reshape(self.ODOC_dict['jr'].shape[0], -1)
            hydro['z1run'] = self.ODOC_dict['z1run'].reshape(self.ODOC_dict['jr'].shape[0], -1)
            hydro['z2run'] = self.ODOC_dict['z2run'].reshape(self.ODOC_dict['jr'].shape[0], -1)
        else:
            hydro['slprun'] = np.zeros([num_steps, 1]) * np.nan
            hydro['x1run'] = np.zeros([num_steps, 1]) * np.nan
            hydro['x2run'] = np.zeros([num_steps, 1]) * np.nan
            hydro['z1run'] = np.zeros([num_steps, 1]) * np.nan
            hydro['z2run'] = np.zeros([num_steps, 1]) * np.nan


        # from OSETUP
        time_end = np.zeros([num_steps, 1]) * np.nan
        x = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        setup = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        depth = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        sigma = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        Hs = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        # from OXVELO
        x_xvelo = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        umean = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        ustd = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        # from OYVELO
        x_yvelo = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        stheta = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        vmean = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        vstd = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan

        for ii in range(0, num_steps):
            # OSETUP
            temp_dict_setup = self.OSETUP_dict['hydro%s' % str(ii + 1)]
            time_end[ii] = temp_dict_setup['time_end']
            x[ii] = temp_dict_setup['x']
            setup[ii] = temp_dict_setup['setup']
            depth[ii] = temp_dict_setup['depth']
            sigma[ii] = temp_dict_setup['sigma']
            Hs[ii] = temp_dict_setup['Hs']
            # OXVELO
            temp_dict_xvel = self.OXVELO_dict['hydro%s' % str(ii + 1)]
            x_xvelo[ii] = temp_dict_xvel['x_xvelo']
            umean[ii] = temp_dict_xvel['umean']
            ustd[ii] = temp_dict_xvel['ustd']
            # OYVELO - so, unlike Brad, I define this regardless, it just has nans everywhere if it is empty
            temp_dict_yvel = self.OYVELO_dict['hydro%s' % str(ii + 1)]
            x_yvelo[ii] = temp_dict_yvel['x_yvelo']
            stheta[ii] = temp_dict_yvel['stheta']
            vmean[ii] = temp_dict_yvel['vmean']
            vstd[ii] = temp_dict_yvel['vstd']

        hydro['time_end'] = time_end
        hydro['x'] = x
        hydro['mwl'] = setup
        hydro['depth'] = depth
        hydro['sigma'] = sigma
        hydro['Hs'] = Hs
        hydro['x_xvelo'] = x_xvelo
        hydro['umean'] = umean
        hydro['ustd'] = ustd
        hydro['x_yvelo'] = x_yvelo
        hydro['stheta'] = stheta
        hydro['vmean'] = vmean
        hydro['vstd'] = vstd


        # sed
        # from OCROSS
        time_end = np.zeros([num_steps, 1]) * np.nan
        x_cross = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qbx = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qsx = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qx = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        # from OLONGS
        x_long = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qby = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qsy = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        qy = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        # from OBSUSL
        x_susl = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        ps = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        pb = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan
        vs = np.zeros([num_steps, self.ODOC_dict['nbinp']]) * np.nan

        if self.ODOC_dict['iprofl'] == 1:
            for ii in range(0, num_steps):
                # OCROSS
                temp_dict_cross = self.OCROSS_dict['sed%s' % str(ii + 1)]
                time_end[ii] = temp_dict_cross['time_end']
                x_cross[ii] = temp_dict_cross['x_cross']
                qbx[ii] = temp_dict_cross['qbx']
                qsx[ii] = temp_dict_cross['qsx']
                qx[ii] = temp_dict_cross['qx']
                # OLONGS
                if 'sed1' in self.OLONGS_dict.keys():
                    temp_dict_long = self.OLONGS_dict['sed%s' % str(ii + 1)]
                    x_long[ii] = temp_dict_long['x_long']
                    qby[ii] = temp_dict_long['qby']
                    qsy[ii] = temp_dict_long['qsy']
                    qy[ii] = temp_dict_long['qy']
                # OBSUSL
                temp_dict_susl = self.OBSUSL_dict['sed%s' % str(ii + 1)]
                x_susl[ii] = temp_dict_susl['x_susl']
                ps[ii] = temp_dict_susl['ps']
                pb[ii] = temp_dict_susl['pb']
                vs[ii] = temp_dict_susl['vs']

        sed['time_end'] = time_end
        sed['x_cross'] = x_cross
        sed['qbx'] = qbx
        sed['qsx'] = qsx
        sed['qx'] = qx
        sed['x_long'] = x_long
        sed['qby'] = qby
        sed['qsy'] = qsy
        sed['qy'] = qy
        sed['x_susl'] = x_susl
        sed['ps'] = ps
        sed['pb'] = pb
        sed['vs'] = vs

        # bc
        bc['time_offshore'] = self.readIF_dict['time_offshore']
        bc['Tp_offshore'] = self.readIF_dict['Tp_bc']
        bc['Hs_offshore'] = self.readIF_dict['Hs_bc']
        bc['wave_setup_offshore'] = self.readIF_dict['setup_bc']
        bc['strm_tide_offshore'] = self.readIF_dict['wl_bc']
        bc['angle_offshore'] = self.readIF_dict['angle_bc']

        return params, bc, veg, hydro, sed
