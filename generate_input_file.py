import numpy as np
import math
import io
import os, subprocess,sys
import shutil
from CMTB_predata import inputOutput_LZ

################################################################################################################################
################################################################################################################################
def makeinfile(vars_dict, infile_folder):

    infile_path = os.path.join(infile_folder, 'infile')
    if os.path.exists("demofile.txt"):
        os.remove(infile_path) # remove file if exist first.
    else:
        infile = open(infile_path, "w")
        # print('Generated CSHORE-VEG input file in: ' + infile_path)

# pre-process input variables
    if math.floor(vars_dict['iprofl']['value'])==0:
        vars_dict['isedav']['value']=0

# bed profile
    if vars_dict['xzb_production']['value']==1: # gui create x
        xloc = np.arange(0, vars_dict['Lx']['value'][0]+vars_dict['dx']['value'][0],
            vars_dict['dx']['value'][0], dtype=np.float64)
        zb = np.interp(xloc,
            [0, vars_dict['toeloc']['value'][0], vars_dict['Lx']['value'][0]],
            [vars_dict['zboffshore']['value'][0],vars_dict['zboffshore']['value'][0],
                vars_dict['zbonshore']['value'][0]])
    elif vars_dict['xzb_production']['value']==2: # user uploaded:
        xloc=vars_dict['upload_x_zb']['value']['x'].values
        zb=vars_dict['upload_x_zb']['value']['zb'].values

    if len(xloc)>5000:
        print("\x1b[31mCSHORE-VEG allows at most 5000 grids along the transect.\n\n\x1b[0m")
        sys.exit()

# update dx value
    vars_dict['dx']['value'][0] = xloc[1] - xloc[0]
    fw = np.array(vars_dict['bottom_fric_factor']['value'] * len(xloc))

# wave conditions
    if vars_dict['simu_type']['value']==2: #batch run
        vars_dict.update({'totaltime': {'value': [vars_dict['upload_wavefiles']['value']['time'].values[-1]]}})
        vars_dict.update({'dt': {'value': [vars_dict['upload_wavefiles']['value']['time'].values[1]-\
            vars_dict['upload_wavefiles']['value']['time'].values[0]]}})
        vars_dict['Hrms']['value'] = vars_dict['upload_wavefiles']['value']['Hrms'].values
        vars_dict['Tp']['value'] = vars_dict['upload_wavefiles']['value']['Tp'].values
        vars_dict['waveangle']['value'] = vars_dict['upload_wavefiles']['value']['waveangle'].values
        vars_dict['wsetup']['value'] = vars_dict['upload_wavefiles']['value']['meanwaterlevel'].values
        vars_dict['swlbc']['value'] = vars_dict['upload_wavefiles']['value']['stillwaterlevel'].values
        vars_dict.update({'nwave': {'value': [len(vars_dict['Hrms']['value'])]}})
        vars_dict.update({'nsurge': {'value': [len(vars_dict['Hrms']['value'])]}})
    elif vars_dict['simu_type']['value']==1: #single run
        vars_dict.update({'totaltime': {'value': [1.0]}})
        vars_dict.update({'dt': {'value': [1.0]}})
        vars_dict.update({'nwave': {'value': [1]}})
        vars_dict.update({'nsurge': {'value': [1]}})

    timebc_wave = np.arange(0, vars_dict['totaltime']['value'][0]+vars_dict['dt']['value'][0],
        vars_dict['dt']['value'][0], dtype=np.float64)

# veg properties
    # if iveg!=3 (i.e., not using veg module in CSHORE-VEG),
    # then veg is treated as rigid and Cd has to be specified.
    if vars_dict['iveg']['value'] != 3:
        vars_dict['ivegCd']['value'] = 0
        vars_dict['ivegtype_sptraits']['value'] = 3

    if vars_dict['iveg']['value'] == 3:
        if vars_dict['ivegtype_sptraits']['value']==1 or vars_dict['ivegtype_sptraits']['value']==2:
            ivegtype = 1
        elif vars_dict['ivegtype_sptraits']['value']==3 or vars_dict['ivegtype_sptraits']['value']==4:
            ivegtype = 0

        # find out vegetated grid id
        vegid = np.where( (xloc>=vars_dict['vegfrom']['value'][0]) &
            (xloc<=vars_dict['vegto']['value'][0]) )
        novegid = np.where( (xloc<vars_dict['vegfrom']['value'][0]) | \
            (xloc>vars_dict['vegto']['value'][0]) )

        if vars_dict['vegNsegtype']['value']==1: #('no spatially varying Cd', 1), ('with spatially varying Cd', 2)
            vars_dict.update({'vegNseg': {'value': 1}})
        elif vars_dict['vegNsegtype']['value']==2: #('no spatially varying Cd', 1), ('with spatially varying Cd', 2)
            vegwidth = vars_dict['vegto']['value'][0] - vars_dict['vegfrom']['value'][0]
            vars_dict.update({'vegNseg': {'value': max(1.0, math.floor(vegwidth / 15.0))}})

# set veg properties
# make sure that vegetation properties are only within the vegetated grid
        if vars_dict['ivegtype_sptraits']['value']==1: #Flexible w/ spatially constant traits
            # ['x', 'Nvs', 'hvs', 'bvs', 'Es', 'Nvb', 'hvb', 'bvb', 'tvb', 'Eb']
            Nvs = np.array(vars_dict['Nvstemflex']['value'] * len(xloc))
            hvs = np.array(vars_dict['hvstemflex']['value'] * len(xloc))
            bvs = np.array(vars_dict['bvstemflex']['value'] * len(xloc))
            Es = np.array(vars_dict['Estem']['value'] * len(xloc))
            Nvb = np.array(vars_dict['Nvblade']['value'] * len(xloc))
            hvb = np.array(vars_dict['hvblade']['value'] * len(xloc))
            bvb = np.array(vars_dict['bvblade']['value'] * len(xloc))
            tvb = np.array(vars_dict['tvblade']['value'] * len(xloc))
            Eb = np.array(vars_dict['Eblade']['value'] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Es[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
            Eb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==2: #Flexible w/ spatially varying traits
            # ['x', 'Nvs', 'hvs', 'bvs', 'Es', 'Nvb', 'hvb', 'bvb', 'tvb', 'Eb']
            Nvs = vars_dict['flex_upload']['value']['Nvs'].values
            hvs = vars_dict['flex_upload']['value']['hvs'].values
            bvs = vars_dict['flex_upload']['value']['bvs'].values
            Es = vars_dict['flex_upload']['value']['Es'].values
            Nvb = vars_dict['flex_upload']['value']['Nvb'].values
            hvb = vars_dict['flex_upload']['value']['hvb'].values
            bvb = vars_dict['flex_upload']['value']['bvb'].values
            tvb = vars_dict['flex_upload']['value']['tvb'].values
            Eb = vars_dict['flex_upload']['value']['Eb'].values
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Es[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
            Eb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==3: #Rigid w/ spatially constant traits
            # ['x', 'Nvs', 'hvs', 'bvs']
            Nvs = np.array(vars_dict['Nvstemrigid']['value'] * len(xloc))
            hvs = np.array(vars_dict['hvstemrigid']['value'] * len(xloc))
            bvs = np.array(vars_dict['bvstemrigid']['value'] * len(xloc))
            Nvb = np.array([0.0] * len(xloc))
            hvb = np.array([0.0] * len(xloc))
            bvb = np.array([0.0] * len(xloc))
            tvb = np.array([0.0] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==4: #Rigid w/ spatially varying traits
            # ['x', 'Nvs', 'hvs', 'bvs']
            Nvs = vars_dict['rigid_upload']['value']['Nvs'].values
            hvs = vars_dict['rigid_upload']['value']['hvs'].values
            bvs = vars_dict['rigid_upload']['value']['bvs'].values
            Nvb = np.array([0.0] * len(xloc))
            hvb = np.array([0.0] * len(xloc))
            bvb = np.array([0.0] * len(xloc))
            tvb = np.array([0.0] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0

# make sure again that veg traits = 0 on non-vegetated id.
        vegrod[novegid] = 0.0

# start with setting up file header and opening input file id.
    spaces = '                        '
    dashes = '--------------------------------------------------\n'

    string = '3 \n' + dashes + 'CSHORE applied to idealized planar slope\n'  +  \
             dashes + \
        str(vars_dict['iline']['value'][0]) + spaces + '->ILINE\n' +  \
        str(vars_dict['iprofl']['value']) + spaces + '->IPROFL\n'

    if math.floor(vars_dict['iprofl']['value'])==1:
        string = string + str(vars_dict['isedav']['value']) + spaces + '->ISEDAV\n'

    string = string + \
        str(vars_dict['iperm']['value']) + spaces + '->IPERM\n' + \
        str(vars_dict['iover']['value']) + spaces + '->IOVER\n'

    if vars_dict['iover']['value']==1:
        string = string + \
            str(vars_dict['iwtran']['value']) + spaces + '->IWTRAN\n'
        if vars_dict['iwtran']['value']==0:
            string = string + \
                str(vars_dict['ipond']['value']) + spaces + '->IPOND\n'

    if vars_dict['iover']['value']==1 and \
       vars_dict['iperm']['value']==0 and \
       math.floor(vars_dict['iprofl']['value'])==1:
            string = string + \
                str(vars_dict['infilt']['value']) + spaces + '->INFILT\n'

    string = string + \
        str(vars_dict['iwcint']['value']) + spaces + '->IWCINT\n' + \
        str(vars_dict['iroll']['value']) + spaces + '->IROLL\n' + \
        str(vars_dict['iwind']['value']) + spaces + '->IWIND\n' + \
        str(vars_dict['itide']['value']) + spaces + '->ITIDE\n' + \
        str(vars_dict['iveg']['value']) + spaces + '->IVEG\n'

    if vars_dict['iveg']['value']==3:
        string = string + \
            str(vars_dict['ivegCd']['value']) + spaces + '->IDVEGCD\n' + \
            str(ivegtype) + spaces + '->IDVEGTYPE\n' + \
            str(vars_dict['CdCap']['value'][0]) + spaces + '->CdCap\n' + \
            str(1000.0) + spaces + '->rhowater\n' + \
            str(vars_dict['vegNseg']['value']) + spaces + '->NVEGSEGMENT\n' + \
            str(vars_dict['ibreaking']['value']) + spaces + '->IBREAKING\n' + \
            str(vars_dict['idiss']['value']) + spaces + '->IDISS\n' + \
            str(vars_dict['iFv']['value']) + spaces + '->IFv\n'

    string = string + \
        str(vars_dict['dx']['value'][0]) + spaces + '->DXC\n' + \
        str(vars_dict['gamma']['value'][0]) + spaces + '->GAMMA\n'

    if vars_dict['iprofl']['value']==1:
        string = string + \
            str(vars_dict['d50']['value'][0]) + ' ' + \
            str(vars_dict['wf']['value'][0]) + ' ' + \
            str(vars_dict['sg']['value'][0]) + spaces + '->D50 WF SG\n'
        string = string + \
            str(vars_dict['effb']['value'][0]) + ' ' + \
            str(vars_dict['efff']['value'][0]) + ' ' + \
            str(vars_dict['slp']['value'][0]) + ' ' + \
            str(vars_dict['slpot']['value'][0]) + spaces + '->EFFB EFFF SLP\n'
        string = string + \
            str(vars_dict['tanphi']['value'][0]) + ' ' + \
            str(vars_dict['blp']['value'][0]) + spaces + '->TANPHI BLP\n'

    if vars_dict['iover']['value']==1:
        string = string + \
            str(vars_dict['rwh']['value'][0]) + spaces + '->RWH\n'

# note: permeable bed is not implemented in this GUI.
# use the script below for permeable bed.
    # if in.iperm;
    #   fprintf(fid, '%11.4f%11.4f%11.4f\n',in.stoneporo, in.stonedia, in.criticalstability )
    # end

# set wave conditions
    string = string + \
        str(vars_dict['ilab']['value']) + spaces + '->ILAB\n'

    # the following script (till but not include NBINP) is for vars_dict['ilab']['value']==1:
    # this script does not have ilab = 0.
    string = string + \
             str(vars_dict['nwave']['value'][0]) + spaces + '->NWAVE\n' + \
             str(vars_dict['nsurge']['value'][0]) + spaces + '->NSURGE\n'

    for ii in range(vars_dict['nwave']['value'][0]):
        if vars_dict['iveg']['value']==3 and vars_dict['idiss']['value']==2:
            freqmin = (1.0/vars_dict['Tp']['value'][ii])*0.1
            freqmax = (1.0/vars_dict['Tp']['value'][ii])*5.0
            numfreq = 500
            jonswapgamma = 3.3

            string = string + \
                str(timebc_wave[ii]) + ' ' +\
                str(vars_dict['Tp']['value'][ii]) + ' ' +\
                str(vars_dict['Hrms']['value'][ii]) + ' ' +\
                str(vars_dict['wsetup']['value'][ii]) + ' ' +\
                str(vars_dict['swlbc']['value'][ii]) + ' ' +\
                str(vars_dict['waveangle']['value'][ii]) + ' ' +\
                str(freqmin) + ' ' +\
                str(freqmax) + ' ' +\
                str(numfreq) + ' ' +\
                str(jonswapgamma) + '\n'
        else:
            string = string + \
                str(timebc_wave[ii]) + ' ' +\
                str(vars_dict['Tp']['value'][ii]) + ' ' +\
                str(vars_dict['Hrms']['value'][ii]) + ' ' +\
                str(vars_dict['wsetup']['value'][ii]) + ' ' +\
                str(vars_dict['swlbc']['value'][ii]) + ' ' +\
                str(vars_dict['waveangle']['value'][ii])+ '\n'

# set bathymetry
    string = string + str(len(xloc)) + spaces + '->NBINP\n'

    for ix, izb, ifw in zip(xloc, zb, fw):
        string = string + str(ix)+' '+str(izb)+' '+str(ifw) + '\n'

# set vegetation
    if vars_dict['iveg']['value']==3:
        if vars_dict['ivegCd']['value'] ==0: #User-defined constant Cd
            Cd=np.array(vars_dict['userCd']['value'] * len(xloc))
            Cdm=Cd
            for icd, icdm in zip(Cd, Cdm):
                string = string + str(icd)+ ' ' + str(icdm) + '\n'

        if ivegtype == 0: # rigid veg.
            for invs, invb, ibvs, ibvb, ihvs, ihvb, itvb, irod in \
                zip(Nvs, Nvb, bvs, bvb, hvs, hvb, tvb, vegrod):
                string = string + ' ' + str(invs)+ ' ' + str(invb) + ' ' + \
                    str(ibvs) + ' ' + str(ibvb) + ' ' + str(ihvs) + ' ' + \
                    str(ihvb) + ' ' + str(itvb) + ' ' + str(irod) + '\n'
        elif ivegtype == 1: # flex veg.
            for invs, invb, ibvs, ibvb, ihvs, ihvb, itvb, ies, ieb, irod in \
                zip(Nvs, Nvb, bvs, bvb, hvs, hvb, tvb, Es, Eb, vegrod):
                string = string + ' ' + str(invs)+ ' ' + str(invb) + ' ' + \
                    str(ibvs) + ' ' + str(ibvb) + ' ' + str(ihvs) + ' ' + \
                    str(ihvb) + ' ' + str(itvb) + ' ' + str(ies) + ' ' + \
                    str(ieb) + ' ' + str(irod) + '\n'

    infile.write(string)
    infile.close()

    if os.path.exists(infile_path) is not True:
        print("\x1b[31mCheck project folder.\x1b[0m")

    if vars_dict['iveg']['value'] == 3:
       if vars_dict['ivegtype_sptraits']['value'] == 999:
           return False

    if (os.path.exists(infile_path)) and \
        (vars_dict['iprofl']['value'] != 999) and \
        (vars_dict['xzb_production'] != 999) and \
        (vars_dict['simu_type']['value'] != 999) and \
        (vars_dict['iveg']['value'] != 999):
       return True
    else:
        return False


################################################################################################################################
################################################################################################################################
def makeinfile_4_evaluation(vars_dict, infile_folder, breakage_fraction):
    # breakage_fraction is a float
    # note: this function is similar to makeinfile(), but replaces veg hv, bv with the values in statistics.

    infile_path = os.path.join(infile_folder, 'infile')
    if os.path.exists("demofile.txt"):
        os.remove(infile_path) # remove file if exist first.
    else:
        infile = open(infile_path, "w")
        # print('Generated CSHORE-VEG input file in: ' + infile_path)

# pre-process input variables
    if math.floor(vars_dict['iprofl']['value'])==0:
        vars_dict['isedav']['value']=0

# bed profile
    if vars_dict['xzb_production']['value']==1: # gui create x
        xloc = np.arange(0, vars_dict['Lx']['value'][0]+vars_dict['dx']['value'][0],
            vars_dict['dx']['value'][0], dtype=np.float64)
        zb = np.interp(xloc,
            [0, vars_dict['toeloc']['value'][0], vars_dict['Lx']['value'][0]],
            [vars_dict['zboffshore']['value'][0],vars_dict['zboffshore']['value'][0],
                vars_dict['zbonshore']['value'][0]])
    elif vars_dict['xzb_production']['value']==2: # user uploaded:
        xloc=vars_dict['upload_x_zb']['value']['x'].values
        zb=vars_dict['upload_x_zb']['value']['zb'].values

    if len(xloc)>5000:
        print("\x1b[31mCSHORE-VEG allows at most 5000 grids along the transect.\n\n\x1b[0m")
        sys.exit()

# update dx value
    vars_dict['dx']['value'][0] = xloc[1] - xloc[0]
    fw = np.array(vars_dict['bottom_fric_factor']['value'] * len(xloc))

# wave conditions
    if vars_dict['simu_type']['value']==2: #batch run
        vars_dict.update({'totaltime': {'value': [vars_dict['upload_wavefiles']['value']['time'].values[-1]]}})
        vars_dict.update({'dt': {'value': [vars_dict['upload_wavefiles']['value']['time'].values[1]-\
            vars_dict['upload_wavefiles']['value']['time'].values[0]]}})
        vars_dict['Hrms']['value'] = vars_dict['upload_wavefiles']['value']['Hrms'].values
        vars_dict['Tp']['value'] = vars_dict['upload_wavefiles']['value']['Tp'].values
        vars_dict['waveangle']['value'] = vars_dict['upload_wavefiles']['value']['waveangle'].values
        vars_dict['wsetup']['value'] = vars_dict['upload_wavefiles']['value']['meanwaterlevel'].values
        vars_dict['swlbc']['value'] = vars_dict['upload_wavefiles']['value']['stillwaterlevel'].values
        vars_dict.update({'nwave': {'value': [len(vars_dict['Hrms']['value'])]}})
        vars_dict.update({'nsurge': {'value': [len(vars_dict['Hrms']['value'])]}})
    elif vars_dict['simu_type']['value']==1: #single run
        vars_dict.update({'totaltime': {'value': [1.0]}})
        vars_dict.update({'dt': {'value': [1.0]}})
        vars_dict.update({'nwave': {'value': [1]}})
        vars_dict.update({'nsurge': {'value': [1]}})

    timebc_wave = np.arange(0, vars_dict['totaltime']['value'][0]+vars_dict['dt']['value'][0],
        vars_dict['dt']['value'][0], dtype=np.float64)

# veg properties
    # if iveg!=3 (i.e., not using veg module in CSHORE-VEG),
    # then veg is treated as rigid and Cd has to be specified.
    if vars_dict['iveg']['value'] != 3:
        vars_dict['ivegCd']['value'] = 0
        vars_dict['ivegtype_sptraits']['value'] = 3

    if vars_dict['iveg']['value'] == 3:
        if vars_dict['ivegtype_sptraits']['value']==1 or vars_dict['ivegtype_sptraits']['value']==2:
            ivegtype = 1
        elif vars_dict['ivegtype_sptraits']['value']==3 or vars_dict['ivegtype_sptraits']['value']==4:
            ivegtype = 0

        # find out vegetated grid id
        vegid = np.where( (xloc>=vars_dict['vegfrom']['value'][0]) &
            (xloc<=vars_dict['vegto']['value'][0]) )
        novegid = np.where( (xloc<vars_dict['vegfrom']['value'][0]) | \
            (xloc>vars_dict['vegto']['value'][0]) )

        if vars_dict['vegNsegtype']['value']==1: #('no spatially varying Cd', 1), ('with spatially varying Cd', 2)
            vars_dict.update({'vegNseg': {'value': 1}})
        elif vars_dict['vegNsegtype']['value']==2: #('no spatially varying Cd', 1), ('with spatially varying Cd', 2)
            vegwidth = vars_dict['vegto']['value'][0] - vars_dict['vegfrom']['value'][0]
            vars_dict.update({'vegNseg': {'value': max(1.0, math.floor(vegwidth / 15.0))}})

# set veg properties
# make sure that vegetation properties are only within the vegetated grid
        if vars_dict['ivegtype_sptraits']['value']==1: #Flexible w/ spatially constant traits
            # ['x', 'Nvs', 'hvs', 'bvs', 'Es', 'Nvb', 'hvb', 'bvb', 'tvb', 'Eb']
            Nvs = np.array(vars_dict['Nvstemflex']['value'] * len(xloc))
            hvs = np.array(vars_dict['hvstemflex']['value'] * len(xloc))
            bvs = np.array(vars_dict['bvstemflex']['value'] * len(xloc))
            Es = np.array(vars_dict['Estem']['value'] * len(xloc))
            Nvb = np.array(vars_dict['Nvblade']['value'] * len(xloc))
            hvb = np.array(vars_dict['hvblade']['value'] * len(xloc))
            bvb = np.array(vars_dict['bvblade']['value'] * len(xloc))
            tvb = np.array(vars_dict['tvblade']['value'] * len(xloc))
            Eb = np.array(vars_dict['Eblade']['value'] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Es[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
            Eb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==2: #Flexible w/ spatially varying traits
            # ['x', 'Nvs', 'hvs', 'bvs', 'Es', 'Nvb', 'hvb', 'bvb', 'tvb', 'Eb']
            Nvs = vars_dict['flex_upload']['value']['Nvs'].values
            hvs = vars_dict['flex_upload']['value']['hvs'].values
            bvs = vars_dict['flex_upload']['value']['bvs'].values
            Es = vars_dict['flex_upload']['value']['Es'].values
            Nvb = vars_dict['flex_upload']['value']['Nvb'].values
            hvb = vars_dict['flex_upload']['value']['hvb'].values
            bvb = vars_dict['flex_upload']['value']['bvb'].values
            tvb = vars_dict['flex_upload']['value']['tvb'].values
            Eb = vars_dict['flex_upload']['value']['Eb'].values
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Es[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
            Eb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==3: #Rigid w/ spatially constant traits
            # ['x', 'Nvs', 'hvs', 'bvs']
            Nvs = np.array(vars_dict['Nvstemrigid']['value'] * len(xloc))
            hvs = np.array(vars_dict['hvstemrigid']['value'] * len(xloc))
            bvs = np.array(vars_dict['bvstemrigid']['value'] * len(xloc))
            Nvb = np.array([0.0] * len(xloc))
            hvb = np.array([0.0] * len(xloc))
            bvb = np.array([0.0] * len(xloc))
            tvb = np.array([0.0] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0
        elif vars_dict['ivegtype_sptraits']['value']==4: #Rigid w/ spatially varying traits
            # ['x', 'Nvs', 'hvs', 'bvs']
            Nvs = vars_dict['rigid_upload']['value']['Nvs'].values
            hvs = vars_dict['rigid_upload']['value']['hvs'].values
            bvs = vars_dict['rigid_upload']['value']['bvs'].values
            Nvb = np.array([0.0] * len(xloc))
            hvb = np.array([0.0] * len(xloc))
            bvb = np.array([0.0] * len(xloc))
            tvb = np.array([0.0] * len(xloc))
            vegrod = np.array(vars_dict['vegrod']['value'] * len(xloc))
            # make sure again that veg traits = 0 on non-vegetated id.
            Nvs[novegid] = 0.0
            hvs[novegid] = 0.0
            bvs[novegid] = 0.0
            Nvb[novegid] = 0.0
            hvb[novegid] = 0.0
            bvb[novegid] = 0.0
            tvb[novegid] = 0.0

# replace veg hv and bv with mean hv, and mean bv from veg statistics accordion
        hvs[vegid] = vars_dict['meanhv']['value'][0]
        bvs[vegid] = vars_dict['meanbv']['value'][0]
        Nvs_tmp = np.array(vars_dict['Nv_in_evaluation']['value'] * len(xloc))
        Nvs = Nvs_tmp *(1.0-np.array(breakage_fraction))
        Nvb[vegid] = 0.0
        hvb[vegid] = 0.0
        bvb[vegid] = 0.0
        tvb[vegid] = 0.0

# make sure again that veg traits = 0 on non-vegetated id.
        vegrod[novegid] = 0.0

# start with setting up file header and opening input file id.
    spaces = '                        '
    dashes = '--------------------------------------------------\n'

    string = '3 \n' + dashes + 'CSHORE applied to idealized planar slope\n'  +  \
             dashes + \
        str(vars_dict['iline']['value'][0]) + spaces + '->ILINE\n' +  \
        str(vars_dict['iprofl']['value']) + spaces + '->IPROFL\n'

    if math.floor(vars_dict['iprofl']['value'])==1:
        string = string + str(vars_dict['isedav']['value']) + spaces + '->ISEDAV\n'

    string = string + \
        str(vars_dict['iperm']['value']) + spaces + '->IPERM\n' + \
        str(vars_dict['iover']['value']) + spaces + '->IOVER\n'

    if vars_dict['iover']['value']==1:
        string = string + \
            str(vars_dict['iwtran']['value']) + spaces + '->IWTRAN\n'
        if vars_dict['iwtran']['value']==0:
            string = string + \
                str(vars_dict['ipond']['value']) + spaces + '->IPOND\n'

    if vars_dict['iover']['value']==1 and \
       vars_dict['iperm']['value']==0 and \
       math.floor(vars_dict['iprofl']['value'])==1:
            string = string + \
                str(vars_dict['infilt']['value']) + spaces + '->INFILT\n'

    string = string + \
        str(vars_dict['iwcint']['value']) + spaces + '->IWCINT\n' + \
        str(vars_dict['iroll']['value']) + spaces + '->IROLL\n' + \
        str(vars_dict['iwind']['value']) + spaces + '->IWIND\n' + \
        str(vars_dict['itide']['value']) + spaces + '->ITIDE\n' + \
        str(vars_dict['iveg']['value']) + spaces + '->IVEG\n'

    if vars_dict['iveg']['value']==3:
        string = string + \
            str(vars_dict['ivegCd']['value']) + spaces + '->IDVEGCD\n' + \
            str(ivegtype) + spaces + '->IDVEGTYPE\n' + \
            str(vars_dict['CdCap']['value'][0]) + spaces + '->CdCap\n' + \
            str(1000.0) + spaces + '->rhowater\n' + \
            str(vars_dict['vegNseg']['value']) + spaces + '->NVEGSEGMENT\n' + \
            str(vars_dict['ibreaking']['value']) + spaces + '->IBREAKING\n' + \
            str(vars_dict['idiss']['value']) + spaces + '->IDISS\n' + \
            str(vars_dict['iFv']['value']) + spaces + '->IFv\n'

    string = string + \
        str(vars_dict['dx']['value'][0]) + spaces + '->DXC\n' + \
        str(vars_dict['gamma']['value'][0]) + spaces + '->GAMMA\n'

    if vars_dict['iprofl']['value']==1:
        string = string + \
            str(vars_dict['d50']['value'][0]) + ' ' + \
            str(vars_dict['wf']['value'][0]) + ' ' + \
            str(vars_dict['sg']['value'][0]) + spaces + '->D50 WF SG\n'
        string = string + \
            str(vars_dict['effb']['value'][0]) + ' ' + \
            str(vars_dict['efff']['value'][0]) + ' ' + \
            str(vars_dict['slp']['value'][0]) + ' ' + \
            str(vars_dict['slpot']['value'][0]) + spaces + '->EFFB EFFF SLP\n'
        string = string + \
            str(vars_dict['tanphi']['value'][0]) + ' ' + \
            str(vars_dict['blp']['value'][0]) + spaces + '->TANPHI BLP\n'

    if vars_dict['iover']['value']==1:
        string = string + \
            str(vars_dict['rwh']['value'][0]) + spaces + '->RWH\n'

# note: permeable bed is not implemented in this GUI.
# use the script below for permeable bed.
    # if in.iperm;
    #   fprintf(fid, '%11.4f%11.4f%11.4f\n',in.stoneporo, in.stonedia, in.criticalstability )
    # end

# set wave conditions
    string = string + \
        str(vars_dict['ilab']['value']) + spaces + '->ILAB\n'

    # the following script (till but not include NBINP) is for vars_dict['ilab']['value']==1:
    # this script does not have ilab = 0.
    string = string + \
             str(vars_dict['nwave']['value'][0]) + spaces + '->NWAVE\n' + \
             str(vars_dict['nsurge']['value'][0]) + spaces + '->NSURGE\n'

    for ii in range(vars_dict['nwave']['value'][0]):
        if vars_dict['iveg']['value']==3 and vars_dict['idiss']['value']==2:
            freqmin = (1.0/vars_dict['Tp']['value'][ii])*0.1
            freqmax = (1.0/vars_dict['Tp']['value'][ii])*5.0
            numfreq = 500
            jonswapgamma = 3.3

            string = string + \
                str(timebc_wave[ii]) + ' ' +\
                str(vars_dict['Tp']['value'][ii]) + ' ' +\
                str(vars_dict['Hrms']['value'][ii]) + ' ' +\
                str(vars_dict['wsetup']['value'][ii]) + ' ' +\
                str(vars_dict['swlbc']['value'][ii]) + ' ' +\
                str(vars_dict['waveangle']['value'][ii]) + ' ' +\
                str(freqmin) + ' ' +\
                str(freqmax) + ' ' +\
                str(numfreq) + ' ' +\
                str(jonswapgamma) + '\n'
        else:
            string = string + \
                str(timebc_wave[ii]) + ' ' +\
                str(vars_dict['Tp']['value'][ii]) + ' ' +\
                str(vars_dict['Hrms']['value'][ii]) + ' ' +\
                str(vars_dict['wsetup']['value'][ii]) + ' ' +\
                str(vars_dict['swlbc']['value'][ii]) + ' ' +\
                str(vars_dict['waveangle']['value'][ii])+ '\n'

# set bathymetry
    string = string + str(len(xloc)) + spaces + '->NBINP\n'

    for ix, izb, ifw in zip(xloc, zb, fw):
        string = string + str(ix)+' '+str(izb)+' '+str(ifw) + '\n'

# set vegetation
    if vars_dict['iveg']['value']==3:
        if vars_dict['ivegCd']['value'] ==0: #User-defined constant Cd
            Cd=np.array(vars_dict['userCd']['value'] * len(xloc))
            Cdm=Cd
            for icd, icdm in zip(Cd, Cdm):
                string = string + str(icd)+ ' ' + str(icdm) + '\n'

        if ivegtype == 0: # rigid veg.
            for invs, invb, ibvs, ibvb, ihvs, ihvb, itvb, irod in \
                zip(Nvs, Nvb, bvs, bvb, hvs, hvb, tvb, vegrod):
                string = string + ' ' + str(invs)+ ' ' + str(invb) + ' ' + \
                    str(ibvs) + ' ' + str(ibvb) + ' ' + str(ihvs) + ' ' + \
                    str(ihvb) + ' ' + str(itvb) + ' ' + str(irod) + '\n'
        elif ivegtype == 1: # flex veg.
            for invs, invb, ibvs, ibvb, ihvs, ihvb, itvb, ies, ieb, irod in \
                zip(Nvs, Nvb, bvs, bvb, hvs, hvb, tvb, Es, Eb, vegrod):
                string = string + ' ' + str(invs)+ ' ' + str(invb) + ' ' + \
                    str(ibvs) + ' ' + str(ibvb) + ' ' + str(ihvs) + ' ' + \
                    str(ihvb) + ' ' + str(itvb) + ' ' + str(ies) + ' ' + \
                    str(ieb) + ' ' + str(irod) + '\n'

    infile.write(string)
    infile.close()

    if os.path.exists(infile_path) is not True:
        print("\x1b[31mCheck project folder.\x1b[0m")

    if vars_dict['iveg']['value'] == 3:
       if vars_dict['ivegtype_sptraits']['value'] == 999:
           return False

    if (os.path.exists(infile_path)) and \
        (vars_dict['iprofl']['value'] != 999) and \
        (vars_dict['xzb_production'] != 999) and \
        (vars_dict['simu_type']['value'] != 999) and \
        (vars_dict['iveg']['value'] != 999):
       return True
    else:
        return False


################################################################################################################################
################################################################################################################################
def runcshore_4_evaluation(all_vars_dict, project_dir):
    current_dir = os.getcwd() # current working directory
    exe_dir = all_vars_dict['cshoreexe']['widget_component'].selected_path
    cshoreexe = all_vars_dict['cshoreexe']['widget_component'].selected

    if cshoreexe is None:
        print("\x1b[31m\nError: CSHORE-VEG executable not assigned.\x1b[0m")
    else:
        # 1. go to project directory,
        # print('project path:' + str(project_dir))
        os.chdir(project_dir)

        # 2. copy "fa_database_noheadline.txt" from cshore exe folder to project folder
        fa_src = os.path.join(exe_dir, 'fa_database_noheadline.txt')
        fa_dst = project_dir

        if os.path.exists(fa_src):
            shutil.copy2(fa_src, fa_dst)
        else:
            print('fa file does not exist')

        # 3. run cshore
        cshore_process = subprocess.run([cshoreexe], capture_output=True, text=True)
        output = cshore_process.stdout.strip()

        if not os.path.exists(os.path.join(project_dir, 'OSETUP')):
            print('CSHORE-VEG failed.\n')

        # 4. go back to current working folder (where ipynb file located)
        os.chdir(current_dir)


################################################################################################################################
################################################################################################################################
def makeinfle_run_load_extract(all_vars_dict, project_dir, breakage_fraction):
# breakage_fraction is a list

  # remove infile and CSHORE-VEG output files (O*) if exist
  for fname in os.listdir(project_dir):
      if fname.startswith("O") or fname == 'infile':
          os.remove(os.path.join(project_dir, fname))

  # generate input file - infile
  _ = makeinfile_4_evaluation(all_vars_dict, project_dir, breakage_fraction)

  # run cshore
  runcshore_4_evaluation(all_vars_dict, project_dir)

  # load cshore results
  cshore_io_O = inputOutput_LZ.cshoreIO()
  params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(project_dir)
  ## only consider first run
  # Hrms_all = hydro0['Hs'][0]/np.sqrt(2)
  # xloc = hydro0['x'][0, :]
  # depth = hydro0['depth'][0, :]
  # ESH_all = veg0['ESH'][0]
  # veg_Cd_all = veg0['Cd'][0]
  ## consider all run
  Hrms_all = hydro0['Hs']/np.sqrt(2)
  xloc = hydro0['x']
  depth = hydro0['depth']
  ESH_all = veg0['ESH']
  veg_Cd_all = veg0['Cd']

  return Hrms_all, ESH_all, veg_Cd_all, xloc, depth
