import ipywidgets as widgets
from ipywidgets import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
import os, subprocess,sys
from scipy.optimize import fsolve
import shutil
from ipyfilechooser import FileChooser
from generate_input_file import makeinfile,makeinfile_4_evaluation,runcshore_4_evaluation,makeinfle_run_load_extract
from CMTB_predata import inputOutput_LZ
import BreakageEvaluation
import BendingStress_ANNModel

################################################################################################
def generate_gridlayout_from_var_dict (var_dict, dependon_string):
    # this function is very much like generate_accordion(). The difference is
    # this function takes var_dict as input, while
    # generate_accordion() needs to find out var_dict from accordion[child_dict].
    # Also, this function only make grid layout for keys that depend on "dependon_string".

    icount = 0
    w_parameter = []
    for key in var_dict:
        if dependon_string in var_dict[key]['dependon']:
            var_label_widget = generate_var_label_widget(var_dict, key)
            if var_dict[key]['info_button']:
                var_info_widget = generate_var_info_widget(var_dict, key)
            else:
                var_info_widget = Label()
        # 'widget_component' for this key will be used later in observe
            var_dict[key]['widget_component'] = var_dict[key]['widget_func'](var_dict, key)

            # if var_dict[key]['value_traits'].get('update_int_float_value'):
            #     var_dict[key]['widget_component'].observe(updateintfloat_handler, names='value')
                # var_dict[key]['value'] = var_dict[key]['widget_component'].value

            if icount % 2 == 0:
                var_label_widget.add_class('label_bg')
            w_parameter.append(var_label_widget)
            # w_parameter.append(var_widget)
            w_parameter.append(var_dict[key]['widget_component'])
            w_parameter.append(var_info_widget)
            icount = icount + 1


    row_height = "35px " * icount
    gridlayout = GridBox(
        children = w_parameter,
        layout = Layout(
            grid_template_columns='max-content max-content max-content',
            grid_template_rows=row_height # note: height should be for each row.
            )
        )

    return gridlayout


def generate_accordion(accordion_dict):
    # this function generates the child of a accordion.
    var_dict = accordion_dict['child_dict']
    icount = 0
    w_parameter = []
    for key in var_dict:
        var_label_widget = generate_var_label_widget(var_dict, key)
        if var_dict[key]['info_button']:
            var_info_widget = generate_var_info_widget(var_dict, key)
        else:
            var_info_widget = Label()

        # 'widget_component' for this key will be used later in observe
        var_dict[key]['widget_component'] = var_dict[key]['widget_func'](var_dict, key)

        # if var_dict[key]['value_traits'].get('update_int_float_value'):
        #     var_dict[key]['widget_component'].observe(updateintfloat_handler, names='value')
            # var_dict[key]['value'] = var_dict[key]['widget_component'].value

        if icount % 2 == 1:
            var_label_widget.add_class('label_bg')
        w_parameter.append(var_label_widget)
        w_parameter.append(var_dict[key]['widget_component'])
        w_parameter.append(var_info_widget)
        icount = icount + 1

    # print(key)
    if (key=='select_project_dir') or (key=='cshoreexe'): # note: this is actually the first key in this var_dict.
        row_height = "auto"
    else:
        row_height = "35px " * len(var_dict)

    coln_width = 'max-content max-content max-content'
    accord_child = GridBox(
        children = w_parameter,
        layout = Layout(
            grid_template_columns=coln_width,
            grid_template_rows=row_height # note: height should be for each row.
            )
        )

    return accord_child

def generate_var_label_widget(var_dict, key):
     lbl = widgets.Label(value=var_dict[key]['label'])
     return lbl


inf_dict = {} # this dictionary associate information button and its  outputs
def generate_var_info_widget(var_dict, key):
    inf_bt = widgets.Button(disabled=False,button_style='info',
                       icon='fa-info-circle', layout=Layout(width='auto'),
                       tooltip=var_dict[key]['tooltip'])
    inf_bt.on_click(on_click_handler_info)
    inf_dict[inf_bt] = widgets.Output()
    inf = HBox([inf_bt, inf_dict[inf_bt]])

    return inf

################## file chooser handler ##################
def change_file_choser_title_project(chooser):
# callback function for file chooser widget
    chooser.title = '<b>project folder selected</b>'

def change_file_choser_title_cshoreexe(chooser):
# callback function for file chooser widget
    chooser.title = '<b>CSHORE-VEG exetutable selected</b>'

################## output handler ##################
def on_click_handler_info(b): # b for button
    if inf_dict[b].outputs: # data is in inf_dict[b].outputs.outputs[0]['data']
        inf_dict[b].clear_output()
    else:
        with inf_dict[b]:
            label = widgets.Label(value=b.tooltip)
            display(label)

def plotBC_handler(b):
    with hydroBCplot_output:
        hydroBCplot_output.clear_output()

        if boundary_vars_dict['simu_type']['widget_component'].value==2: # batch run:
            for filename in boundary_vars_dict['upload_wavefiles']['widget_component'].value:
                boundary_vars_dict['upload_wavefiles']['value'] = pd.read_csv(
                    io.BytesIO(boundary_vars_dict['upload_wavefiles']['widget_component'].value[filename]['content']),
                    names=boundary_vars_dict['upload_wavefiles']['value_traits']['filecontent_str'],
                    header = 0)

            time=boundary_vars_dict['upload_wavefiles']['value']['time'].values
            Hrms=boundary_vars_dict['upload_wavefiles']['value']['Hrms'].values
            Tp=boundary_vars_dict['upload_wavefiles']['value']['Tp'].values


            fig = plt.figure(boundary_vars_dict['plotBC']['figname'],
                figsize=(6,4), dpi=80)
            ax1 = fig.add_subplot(1,1,1)
            ax2=ax1.twinx()
            ax1.plot(time, Hrms, color="red", marker="o")
            ax2.plot(time, Tp, color="blue",marker="o")
            ax1.set_xlabel("time",fontsize=14)
            ax1.set_ylabel("Hrms (m)",color="red",fontsize=14)
            ax2.set_ylabel("Tz (s)",color="blue",fontsize=14)
            # plt.xlabel('x (m)')
            # plt.ylabel('Tz (s)')
            plt.show()
            fig.canvas.draw_idle()

def clear_plotBC_handler(b):
    with hydroBCplot_output:
        fig = plt.figure(boundary_vars_dict['plotBC']['figname'],
            figsize=(0.1,0.1))
        plt.close()

        hydroBCplot_output.clear_output()


def plotzbgui_handler(b):
# because [key]['value'] has not been updated, I used [key]['widget_component'].value
    xloc = np.arange(0, bathymetry_vars_dict['Lx']['widget_component'].value + \
        bathymetry_vars_dict['dx']['widget_component'].value,
        bathymetry_vars_dict['dx']['widget_component'].value, dtype=np.float64)
    zb = np.interp(xloc,
        [0, bathymetry_vars_dict['toeloc']['widget_component'].value,
        bathymetry_vars_dict['Lx']['widget_component'].value],
        [bathymetry_vars_dict['zboffshore']['widget_component'].value,
        bathymetry_vars_dict['zboffshore']['widget_component'].value,
            bathymetry_vars_dict['zbonshore']['widget_component'].value])

    with zbplot_output:
        zbplot_output.clear_output()

        if len(xloc)>5000:
            print("\x1b[31mCSHORE-VEG allows at most 5000 grids along the transect.\n\n\x1b[0m")
            sys.exit()
        else:
            fig = plt.figure(bathymetry_vars_dict['plotzb_gui']['figname'],
                figsize=(6,4), dpi=80)
            ax = fig.add_subplot(1,1,1)
            ax.plot(xloc, zb)
            plt.xlabel('x (m)')
            plt.ylabel('elevation (m)')
            plt.show()
            fig.canvas.draw_idle()

def clear_plotzbgui_handler(b):
    with zbplot_output:
        fig = plt.figure(bathymetry_vars_dict['plotzb_gui']['figname'],
            figsize=(0.1,0.1))
        plt.close()

        zbplot_output.clear_output()

def plotzbupload_handler(b):
# because [key]['value'] has not been updated, I used [key]['widget_component'].value
    with zbplot_output:
        zbplot_output.clear_output()

        for filename in bathymetry_vars_dict['upload_x_zb']['widget_component'].value:
            bathymetry_vars_dict['upload_x_zb']['value'] = pd.read_csv(
                io.BytesIO(bathymetry_vars_dict['upload_x_zb']['widget_component'].value[filename]['content']),
                names=bathymetry_vars_dict['upload_x_zb']['value_traits']['filecontent_str'],
                header = 0)

        xloc=bathymetry_vars_dict['upload_x_zb']['value']['x'].values
        zb=bathymetry_vars_dict['upload_x_zb']['value']['zb'].values


        if len(xloc)>5000:
            print("\x1b[31mCSHORE-VEG allows at most 5000 grids along the transect.\n\n\x1b[0m")
            sys.exit()
        else:
            fig = plt.figure(bathymetry_vars_dict['plotzb_upload']['figname'],
                figsize=(6,4), dpi=80)
            ax = fig.add_subplot(1,1,1)
            ax.plot(xloc, zb)
            plt.xlabel('x (m)')
            plt.ylabel('elevation (m)')
            plt.show()
            fig.canvas.draw_idle()

def clear_plotzbupload_handler(b):
    with zbplot_output:
        fig = plt.figure(bathymetry_vars_dict['plotzb_upload']['figname'],
            figsize=(0.1,0.1))
        plt.close()

        zbplot_output.clear_output()

def vegonoff_handler(change):
    vegonoff_output.clear_output()

    if change.new == 0:
        tmp_wdgt = widgets.Label()
    elif change.new == 3:
        tmp_wdgt = widgets.VBox([option_vegon,
            vegtype_sptraits_output,
            Cd_output])
    elif change.new == 999:
        tmp_wdgt = widgets.Label()

    with vegonoff_output:
        display(tmp_wdgt)

def vegtype_sptraits_handler(change):
    vegtype_sptraits_output.clear_output()

    if change.new == 1:
        tmp_wdgt = option_flexconstant
    elif change.new == 2:
        tmp_wdgt = option_flexvarying
    elif change.new == 3:
        tmp_wdgt = option_rigidconstant
    elif change.new == 4:
        tmp_wdgt = option_rigidvarying
    elif change.new == 999:
        tmp_wdgt = widgets.Label()

    with vegtype_sptraits_output:
        display(tmp_wdgt)

def Cd_handler(change):
    Cd_output.clear_output()

    if change.new == 0:
        tmp_wdgt = option_userCd
    else:
        tmp_wdgt = widgets.Label()

    with Cd_output:
        display(tmp_wdgt)

def morpho_handler(change):

    movable_output.clear_output()

    if change.new == 0:
        tmp_wdgt = widgets.Label()
    elif change.new == 999:
        tmp_wdgt = widgets.Label()
    else: # movable bed
        tmp_wdgt = option_movablebed

    with movable_output:
        display(tmp_wdgt)


def singlebatch_handler(change):
    singlebatch_output.clear_output()

    if change.new == 1:
        tmp_wdgt = option_single
    elif change.new == 2:
        tmp_wdgt = option_batch
    elif change.new == 999:
        tmp_wdgt = widgets.Label()

    with singlebatch_output:
        display(tmp_wdgt)

def xzbproduction_handler(change):
    guiuser_output.clear_output()

    if change.new == 2:
        tmp_wdgt = option_uploadxzb
    elif change.new == 1:
        tmp_wdgt = option_guigenerated
    elif change.new == 999:
        tmp_wdgt = widgets.Label()

    with guiuser_output:
        display(tmp_wdgt)


def makeinfile_handler(b):
    # update values in 'int text', 'float text', 'bounded int text'
    # the variables that need value update are marked with ['value_traits']['update_int_float_value'] = True
    for iaccord_dict in input_tab_dict['child_dict']:
        for key in iaccord_dict['child_dict']:
            if iaccord_dict['child_dict'][key]['value_traits'].get('update_int_float_value'):
                if iaccord_dict['child_dict'][key]['value_traits'].get('data_type')=='array':
                    iaccord_dict['child_dict'][key]['value'] = \
                        [iaccord_dict['child_dict'][key]['widget_component'].value]
                else:
                    iaccord_dict['child_dict'][key]['value'] = \
                        iaccord_dict['child_dict'][key]['widget_component'].value

            if iaccord_dict['child_dict'][key]['value_traits'].get('update_uploadfile'):
                for filename in iaccord_dict['child_dict'][key]['widget_component'].value:
                    iaccord_dict['child_dict'][key]['value'] = pd.read_csv(
                        io.BytesIO(iaccord_dict['child_dict'][key]['widget_component'].value[filename]['content']),
                        names=iaccord_dict['child_dict'][key]['value_traits']['filecontent_str'],
                        header = 0)

    # merging xxx_vars_dict to one dictionary
    all_vars_dict = physical_processes_vars_dict.copy()
    all_vars_dict.update(morphology_vars_dict)
    all_vars_dict.update(vegetation_vars_dict)
    all_vars_dict.update(boundary_vars_dict)
    all_vars_dict.update(bathymetry_vars_dict)

    # generate input file
    project_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

    with makeinfile_output:
        makeinfile_output.clear_output()

        if project_dir is None:
            print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
            print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
            project_dir = os.getcwd()
            print('Project folder is in: \n' + str(project_dir))

        # BEFORE GENERATING infile, CHECK WHETHER ESSENTIAL VARIABLES ARE SET.
        # for instance, if users forget to set "Movable bottom:", it should be
        # detected in this step.

        if morphology_vars_dict['iprofl']['value'] == 999:
            print("\x1b[31mERROR: Select an option for 'Movable bottom' under 'MORPHOLOGICAL CHANGES'.\n\n\x1b[0m")
            sys.exit()

        if bathymetry_vars_dict['xzb_production'] == 999:
            print("\x1b[31mERROR: Bathymetry not set. Select an option for 'Grid & depth production method' under 'BATHYMETRIC & COMPUTATIONAL GRIDS'.\n\n\x1b[0m")
            sys.exit()

        if boundary_vars_dict['simu_type']['value'] == 999:
            print("\x1b[31mERROR: Wave conditions not set. Select an option for 'Simulation type' under 'SEAWARD BOUNDARY CONDITIONS'.\n\n\x1b[0m")
            sys.exit()

        if vegetation_vars_dict['iveg']['value'] == 999:
            print("\x1b[31mERROR: Vegetation module not set. Select an option for 'Vegetation module' under 'VEGETATION-RELATED PARAMETERS'.\n\n\x1b[0m")
            sys.exit()

        if vegetation_vars_dict['iveg']['value'] == 3 and vegetation_vars_dict['ivegtype_sptraits']['value'] == 999:
            print("\x1b[31mERROR: Vegetation type not set. Select an option for 'Vegetation type' under 'VEGETATION-RELATED PARAMETERS'.\n\n\x1b[0m")
            sys.exit()

# remove infile and CSHORE-VEG output files (O*) if exist
        for fname in os.listdir(project_dir):
            if fname.startswith("O") or fname == 'infile':
                os.remove(os.path.join(project_dir, fname))

# generate input file - infile
        infile_status = makeinfile(all_vars_dict, project_dir)

        if infile_status is True:
            print('Generated CSHORE-VEG input file in: ' + str(project_dir))
        else:
            print("\x1b[31mERROR: input file is not genearted.\x1b[0m")


def runcshore_handler(b):
    current_dir = os.getcwd() # current working directory
    project_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path
    exe_dir = cshoreexe_vars_dict['cshoreexe']['widget_component'].selected_path
    cshoreexe = cshoreexe_vars_dict['cshoreexe']['widget_component'].selected

    with runcshore_output:
        runcshore_output.clear_output()

        if project_dir is None:
            # python print with color:
            # https://stackoverflow.com/questions/16816013/is-it-possible-to-print-using-different-colors-in-ipythons-notebook
            print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
            print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
            project_dir = os.getcwd()

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
            print('\nRunning CSHORE-VEG. Wait ...\n\n')
            cshore_process = subprocess.run([cshoreexe], capture_output=True, text=True)
            output = cshore_process.stdout.strip()

            if os.path.exists(os.path.join(project_dir, 'OSETUP')):
                print('CSHORE-VEG completed.\n')
                print('Results are saved in: \n' + str(project_dir))
            else:
                print('CSHORE-VEG is not completed.\n')
                print('Check input file (infile) in: \n' + str(project_dir))


            # 4. go back to current working folder (where ipynb file located)
            os.chdir(current_dir)


def update_plot(new_slider_val):
    data_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

    if data_dir is None:
        with getdata_output:
            getdata_output.clear_output()
            print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
            print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
        data_dir = os.getcwd()
    else:
        # with getdata_output:
        getdata_output.clear_output()

    # with cshoreIO, we actually dont need get_OSETUPdata anymore.
    output_path = os.path.join(data_dir, 'OSETUP')
    if not os.path.exists(output_path):
        with HrmsMWLplot_output:
            print("\x1b[31mCSHORE-VEG output file does not exist. Data is not saved.\n\n\x1b[0m")
    else:
        cshore_io_O = inputOutput_LZ.cshoreIO()
        params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)

        # get the number of wave runs and domain length, so that we can update the slider and gage location box
        num_of_allrun = len(hydro0['x'][:,0])
        maxLx = hydro0['x'][0][-1]

        spatialplot_vars_dict['numofwaverun']['widget_component'].max = num_of_allrun
        # spatialplot_vars_dict['gageloc']['widget_component'].max = maxLx

        with HrmsMWLplot_output:
            # HrmsMWLplot_output.clear_output()

            plottype = spatialplot_vars_dict['plottype']['widget_component'].value
            numofrun = spatialplot_vars_dict['numofwaverun']['widget_component'].value

            fig1 = plt.figure(spatialplot_vars_dict['plottype']['figname'], figsize=(6,4), dpi=80)
            ax23 = fig1.add_subplot(111)

            ax23.clear()

            if plottype=='root-mean-square wave height':
                ax23.plot( hydro0['x'][numofrun-1], hydro0['Hs'][numofrun-1] / np.sqrt(2) )
                ax23.set_ylabel('Hrms (m)', fontsize="large", fontweight="bold")
            elif plottype=='mean water level':
                ax23.plot( hydro0['x'][numofrun-1], hydro0['mwl'][numofrun-1] )
                ax23.set_ylabel('mean water level (m)', fontsize="large", fontweight="bold")
            elif plottype=='drag coefficient':
                ax23.plot(veg0['x'][numofrun-1], veg0['Cd'][numofrun-1])
                ax23.set_ylabel('vegetal drag coefficient', fontsize="large", fontweight="bold")
            elif plottype=='effective plant height':
                ax23.plot(veg0['x'][numofrun-1], veg0['EPH'][numofrun-1])
                ax23.plot(veg0['x'][numofrun-1], veg0['Hstem'][numofrun-1]+veg0['Hbld'][numofrun-1])
                ax23.set_ylabel('plant height(m)', fontsize="large", fontweight="bold")
                ax23.legend(['effective plant height', 'actual plant height'], loc="upper right")


            fig1.suptitle('Run ' + str(numofrun), fontsize="x-large",
                          fontweight="bold" )
            # ax23.legend([plottype +'\nwave run ' + str(numofrun)], loc="upper right")
            ax23.set_xlabel('x (m)', fontsize="large", fontweight="bold")
            # plt.xlim(vegetation_vars_dict['vegfrom']['value'][0], vegetation_vars_dict['vegto']['value'][0])


        with R2p_output:
            R2p_output.clear_output()

            # cshore_io_O = inputOutput_LZ.cshoreIO()
            # params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)
            if 'runup_2_percent' in hydro0.keys():
                print('In run '+str(numofrun) + ', 2% Run-up = '+ str(hydro0['runup_2_percent'][numofrun-1][0]) + ' m')
            else:
                print('no runup.')



def update_plottype(change):

    numofrun = spatialplot_vars_dict['numofwaverun']['widget_component'].value
    data_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

    if data_dir is None:
        with getdata_output:
            getdata_output.clear_output()
            print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
            print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
        data_dir = os.getcwd()


    # with cshoreIO, we actually dont need get_OSETUPdata anymore.
    cshore_io_O = inputOutput_LZ.cshoreIO()
    params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)

    with savespatialdata_output:
        savespatialdata_output.clear_output()

    with HrmsMWLplot_output:
        fig1 = plt.figure(spatialplot_vars_dict['plottype']['figname'], figsize=(6,4), dpi=80)
        ax23 = fig1.add_subplot(111)

        ax23.clear()
        if change.new=='root-mean-square wave height':
            ax23.plot( hydro0['x'][numofrun-1], hydro0['Hs'][numofrun-1] / np.sqrt(2) )
            ax23.set_ylabel('Hrms (m)', fontsize="large", fontweight="bold")
        elif change.new=='mean water level':
            ax23.plot( hydro0['x'][numofrun-1], hydro0['mwl'][numofrun-1] )
            ax23.set_ylabel('mean water level (m)', fontsize="large", fontweight="bold")
        elif change.new=='drag coefficient':
            ax23.plot(veg0['x'][numofrun-1], veg0['Cd'][numofrun-1])
            ax23.set_ylabel('vegetal drag coefficient', fontsize="large", fontweight="bold")
        elif change.new=='effective plant height':
            ax23.plot(veg0['x'][numofrun-1], veg0['EPH'][numofrun-1])
            ax23.plot(veg0['x'][numofrun-1], veg0['Hstem'][numofrun-1]+veg0['Hbld'][numofrun-1])
            ax23.set_ylabel('plant height(m)', fontsize="large", fontweight="bold")
            ax23.legend(['effective plant height', 'actual plant height'], loc="upper right")

        plottype = spatialplot_vars_dict['plottype']['widget_component'].value

        fig1.suptitle(plottype +'\nwave run ' + str(numofrun), fontsize="x-large",
                      fontweight="bold" )
        # ax23.legend([plottype +'\nwave run ' + str(numofrun)], loc="upper right")
        ax23.set_xlabel('x (m)', fontsize="large", fontweight="bold")


def savespatialdata_handler(b):
    with savespatialdata_output:
        savespatialdata_output.clear_output()

        plottype = spatialplot_vars_dict['plottype']['widget_component'].value
        numofrun = spatialplot_vars_dict['numofwaverun']['widget_component'].value
        data_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

        if data_dir is None:
            data_dir = os.getcwd()

        # with cshoreIO, we actually dont need get_OSETUPdata anymore.
        cshore_io_O = inputOutput_LZ.cshoreIO()
        params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)

        xgrid=hydro0['x'][numofrun-1][:]
        ebh = veg0['EBH'][numofrun-1][:]
        esh = veg0['ESH'][numofrun-1][:]
        eph = veg0['EPH'][numofrun-1][:]
        cd = veg0['Cd'][numofrun-1][:]
        Hrms = hydro0['Hs'][numofrun-1][:]/np.sqrt(2)
        mwl = hydro0['mwl'][numofrun-1][:]

        df = pd.DataFrame()

        df['x']=pd.Series(xgrid)
        df['Hrms']=pd.Series(Hrms)
        df['MWL']=pd.Series(mwl)
        df['Cd']=pd.Series(cd)
        df['ESH']=pd.Series(esh)
        df['EBH']=pd.Series(ebh)
        df['EPH']=pd.Series(eph)

        savefile_path = os.path.join(data_dir, 'Run'+str(numofrun)+'.csv')

        # remove file if existing already
        for fname in os.listdir(data_dir):
            if fname == 'Run'+str(numofrun)+'.csv':
                os.remove(savefile_path)

        # write the header
        header1 = 'Run' + str(numofrun) + '\n'
        header2 = 'x (m), Hrms (m), MWL (m), Cd, ESH (m), EBH (m), EPH (m)\n'
        with open(savefile_path, 'w') as fp:
            fp.write(header1)
            fp.write(header2)

        # write the rest
        df.to_csv(savefile_path, header= False,mode='a',index=False)

        print('Hrms, MWL, Cd & EPH from Run '+ str(numofrun) + ' are extracted from CSHORE-VEG outputs in: \n' + str(data_dir))
        print('and saved to: \n' + str(savefile_path))



def savetemporaldata_handler(b):
    with savetemporaldata_output:
        savetemporaldata_output.clear_output()

        gageloc = temporalplot_vars_dict['gageloct']['widget_component'].value
        data_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

        if data_dir is None:
            print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
            print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
            data_dir = os.getcwd()

        # with cshoreIO, we actually dont need get_OSETUPdata anymore.
        cshore_io_O = inputOutput_LZ.cshoreIO()
        params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)

        num_of_allrun = len(hydro0['x'][:,0])
        print('Total number of runs = ' + str(num_of_allrun))

        print('Saving CSHORE-VEG outputs (Hrms, MWL, Cd & EPH) at x = ' + str(gageloc) + ' m and R2% from all runs. Wait ...\n\n')

        CSHORE_outputs = {'x':[], 'Hrms': [], 'MWL': [], 'Cd':[], 'ESH':[], 'EBH':[], 'EPH': []}

        df_all = pd.DataFrame(CSHORE_outputs)

        for irun in range(1, num_of_allrun+1):
            xgrid=hydro0['x'][irun-1][:]
            ebh = veg0['EBH'][irun-1][:]
            esh = veg0['ESH'][irun-1][:]
            eph = veg0['EPH'][irun-1][:]
            cd = veg0['Cd'][irun-1][:]
            Hrms = hydro0['Hs'][irun-1][:]/np.sqrt(2)
            mwl = hydro0['mwl'][irun-1][:]

            df = pd.DataFrame()

            df['x']=pd.Series(xgrid)
            df['Hrms']=pd.Series(Hrms)
            df['MWL']=pd.Series(mwl)
            df['Cd']=pd.Series(cd)
            df['ESH']=pd.Series(esh)
            df['EBH']=pd.Series(ebh)
            df['EPH']=pd.Series(eph)

            # interpolate all values at gage location
            # this way only works when there is only unique values in x column (which is true in our case).
            # if not unique, use the method in: https://stackoverflow.com/questions/54320999/pandas-python-interpolation-of-multiple-columns-based-on-values-specified-for-o
            xx=[gageloc]
            df = df.set_index('x')
            df = df.reindex(df.index.union(xx)).sort_index(ascending=True).interpolate(method='index')
            df = df.reset_index() # remember to reset index (otherwise, 'x' column cannot be find)
            gagevars = df.loc[df['x']==gageloc]
            df_all = df_all.append(gagevars, ignore_index=True)

        # shift index from 0:13 to 1:14 so that it will correspond to run number
        df_all.index = np.arange(1, len(df_all)+1)

        savefile_path = os.path.join(data_dir, 'var_at_gage'+str(gageloc)+'.csv')

        # remove file if existing already
        for fname in os.listdir(data_dir):
            if fname == 'var_at_gage'+str(gageloc)+'.csv':
                os.remove(savefile_path)

        # write the header
        header = 'Run, x (m), Hrms (m), MWL (m), Cd, ESH (m), EBH (m), EPH (m)\n'
        with open(savefile_path, 'w') as fp:
            fp.write(header)

        # write the rest
        df_all.to_csv(savefile_path, header=False, mode='a', index=True)

        # save wave runup for all cases
        if 'runup_2_percent' in hydro0.keys():
            r2p_all=hydro0['runup_2_percent'][irun-1][:]
            runid=np.arange(1,num_of_allrun+1)
            r2p_all = hydro0['runup_2_percent'][:,0]
            df = pd.DataFrame({'Run': runid, 'R2% (m)':r2p_all})
            saver2pfile_path = os.path.join(data_dir, 'WaveRunup_all_runs.csv')

            # remove wave runup file if existing already
            for fname in os.listdir(data_dir):
                if fname == 'WaveRunup_all_runs.csv':
                    os.remove(saver2pfile_path)

            df.to_csv(saver2pfile_path, index=False)


        print('Data is extracted from CSHORE-VEG outputs in: \n' + str(data_dir) + '\n')
        print('and saved to: \n' + str(savefile_path) + '\n')
        print('and ' + str(saver2pfile_path) )



def clearplot_handler(b):
    with HrmsMWLplot_output:
        fig = plt.figure(spatialplot_vars_dict['plottype']['figname'],
            figsize=(0.1,0.1))
        plt.close()
        HrmsMWLplot_output.clear_output()

    with savespatialdata_output:
        savespatialdata_output.clear_output()

    with getdata_output:
        getdata_output.clear_output()



def prepareinputs_4_evaluation():
    # update values in 'int text', 'float text', 'bounded int text'
    # the variables that need value update are marked with ['value_traits']['update_int_float_value'] = True
    for iaccord_dict in input_tab_dict['child_dict']:
        for key in iaccord_dict['child_dict']:
            if iaccord_dict['child_dict'][key]['value_traits'].get('update_int_float_value'):
                if iaccord_dict['child_dict'][key]['value_traits'].get('data_type')=='array':
                    iaccord_dict['child_dict'][key]['value'] = \
                        [iaccord_dict['child_dict'][key]['widget_component'].value]
                else:
                    iaccord_dict['child_dict'][key]['value'] = \
                        iaccord_dict['child_dict'][key]['widget_component'].value

            if iaccord_dict['child_dict'][key]['value_traits'].get('update_uploadfile'):
                for filename in iaccord_dict['child_dict'][key]['widget_component'].value:
                    iaccord_dict['child_dict'][key]['value'] = pd.read_csv(
                        io.BytesIO(iaccord_dict['child_dict'][key]['widget_component'].value[filename]['content']),
                        names=iaccord_dict['child_dict'][key]['value_traits']['filecontent_str'],
                        header = 0)

    # merging xxx_vars_dict to one dictionary
    all_vars_dict = physical_processes_vars_dict.copy()
    all_vars_dict.update(morphology_vars_dict)
    all_vars_dict.update(vegetation_vars_dict)
    all_vars_dict.update(boundary_vars_dict)
    all_vars_dict.update(bathymetry_vars_dict)
    all_vars_dict.update(vegstatistics_vars_dict) # update hv, bv and Nv later
    all_vars_dict.update(evaluationresults_vars_dict)
    all_vars_dict.update(select_projectdir_vars_dict)
    all_vars_dict.update(cshoreexe_vars_dict)


    # generate input file
    project_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

    if project_dir is None:
        print("\x1b[31mWarning: Project folder is not assigned.\x1b[0m")
        print("\x1b[31mUse current directory as the default project folder.\n\x1b[0m")
        project_dir = os.getcwd()
        print('Project folder is in: \n' + str(project_dir))

    # BEFORE GENERATING infile, CHECK WHETHER ESSENTIAL VARIABLES ARE SET.
    # for instance, if users forget to set "Movable bottom:", it should be
    # detected in this step.

    if morphology_vars_dict['iprofl']['value'] == 999:
        print("\x1b[31mERROR: Select an option for 'Movable bottom' under 'MORPHOLOGICAL CHANGES'.\n\n\x1b[0m")
        sys.exit()

    if bathymetry_vars_dict['xzb_production'] == 999:
        print("\x1b[31mERROR: Bathymetry not set. Select an option for 'Grid & depth production method' under 'BATHYMETRIC & COMPUTATIONAL GRIDS'.\n\n\x1b[0m")
        sys.exit()

    if boundary_vars_dict['simu_type']['value'] == 999:
        print("\x1b[31mERROR: Wave conditions not set. Select an option for 'Simulation type' under 'SEAWARD BOUNDARY CONDITIONS'.\n\n\x1b[0m")
        sys.exit()

    if vegetation_vars_dict['iveg']['value'] == 999:
        print("\x1b[31mERROR: Vegetation module not set. Select an option for 'Vegetation module' under 'VEGETATION-RELATED PARAMETERS'.\n\n\x1b[0m")
        sys.exit()

    if vegetation_vars_dict['iveg']['value'] == 3 and vegetation_vars_dict['ivegtype_sptraits']['value'] == 999:
        print("\x1b[31mERROR: Vegetation type not set. Select an option for 'Vegetation type' under 'VEGETATION-RELATED PARAMETERS'.\n\n\x1b[0m")
        sys.exit()

    return all_vars_dict, project_dir


def dispersion(waterdepth, waveperiod):

    def func(k, *arg):
        h, T, g = arg # h = waterdepth, T = waveperiod, g = gravity acceleration
        omeg = 2.0*np.pi / T
        return omeg**2.0 - g * k * np.tanh(k*h)

    grav = 9.806

    # initial guess of solution
    L0 = waveperiod**2.0 * grav / 2.0 / np.pi
    k0 = 2.0*np.pi / L0

    data = (waterdepth, waveperiod, grav)
    wavenumber = fsolve(func, k0, args=data)

    return wavenumber[0]

def func_bendingmoment_LWT(a, k):
    # this function is related bendingmoment from linear wave theory
    return ((a*k*np.sinh(2.0*a*k))/4.0 - np.cosh(2.0*a*k)/8.0 + 1.0/8.0)/k**2.0 + a**2.0/4.0


def clearevaluationplot_handler(b):
    with evaluation_output:

        fig11 = plt.figure(evaluationresults_vars_dict['evaluation']['figname'],
            figsize=(6,4), dpi=80)
        plt.close()
        evaluation_output.clear_output()

    with saveevaluation_output:
        saveevaluation_output.clear_output()


def saveevaluationresults_handler(b):
    with saveevaluation_output:
        saveevaluation_output.clear_output()

        data_dir = select_projectdir_vars_dict['select_project_dir']['widget_component'].selected_path

        if data_dir is None:
            data_dir = os.getcwd()

        # with cshoreIO, we actually dont need get_OSETUPdata anymore.
        cshore_io_O = inputOutput_LZ.cshoreIO()
        params0, bc0, veg0, hydro0, sed0 = cshore_io_O.load_CSHORE_results(data_dir)

        xgrid=hydro0['x'][numofrun-1][:]
        ebh = veg0['EBH'][numofrun-1][:]
        esh = veg0['ESH'][numofrun-1][:]
        eph = veg0['EPH'][numofrun-1][:]
        cd = veg0['Cd'][numofrun-1][:]
        Hrms = hydro0['Hs'][numofrun-1][:]/np.sqrt(2)
        mwl = hydro0['mwl'][numofrun-1][:]

        df = pd.DataFrame()

        df['x']=pd.Series(xgrid)
        df['Hrms']=pd.Series(Hrms)
        df['MWL']=pd.Series(mwl)
        df['Cd']=pd.Series(cd)
        df['ESH']=pd.Series(esh)
        df['EBH']=pd.Series(ebh)
        df['EPH']=pd.Series(eph)

        savefile_path = os.path.join(data_dir, 'Run'+str(numofrun)+'.csv')

        # remove file if existing already
        for fname in os.listdir(data_dir):
            if fname == 'Run'+str(numofrun)+'.csv':
                os.remove(savefile_path)

        # write the header
        header1 = 'Run' + str(numofrun) + '\n'
        header2 = 'x (m), Hrms (m), MWL (m), Cd, ESH (m), EBH (m), EPH (m)\n'
        with open(savefile_path, 'w') as fp:
            fp.write(header1)
            fp.write(header2)

        # write the rest
        df.to_csv(savefile_path, header= False,mode='a',index=False)

        print('Hrms, MWL, Cd & EPH from Run '+ str(numofrun) + ' are extracted from CSHORE-VEG outputs in: \n' + str(data_dir))
        print('and saved to: \n' + str(savefile_path))



def evaluate_veg_breakage(b):

    # this evaluation is only valid for single wave condition
    if boundary_vars_dict['simu_type']['value'] >= 1:

        # define some parameters
        rhowater = 1000.0
        Ac = 2.1

        # create vegetation sample
        veg_hv = vegstatistics_vars_dict['meanhv']['value'][0]
        veg_hv_std = vegstatistics_vars_dict['stdhv']['value'][0]
        veg_bv = vegstatistics_vars_dict['meanbv']['value'][0]
        veg_bv_std = vegstatistics_vars_dict['stdbv']['value'][0]
        veg_flexstre = vegstatistics_vars_dict['meanflexuralstrength']['value'][0]
        veg_flexstre_std = vegstatistics_vars_dict['stdflexuralstrength']['value'][0]
        rho_hv_bv = vegstatistics_vars_dict['corrcoeff_hv_bv']['value'][0]
        rho_hv_flex = vegstatistics_vars_dict['corrcoeff_hv_flex']['value'][0]
        rho_bv_flex = vegstatistics_vars_dict['corrcoeff_bv_flex']['value'][0]

        mu_bv = BreakageEvaluation.muconv(veg_bv,veg_bv_std**2.0)
        sigma_bv = BreakageEvaluation.sigmacov(veg_bv,veg_bv_std**2.0)
        mu_hv = BreakageEvaluation.muconv(veg_hv,veg_hv_std**2.0)
        sigma_hv = BreakageEvaluation.sigmacov(veg_hv,veg_hv_std**2.0)
        mu_flexuralstrength = BreakageEvaluation.muconv(veg_flexstre,veg_flexstre_std**2.0)
        sigma_flexuralstrength = BreakageEvaluation.sigmacov(veg_flexstre,veg_flexstre_std**2.0)

        mu_hv_bv_flexuralstrength = [mu_hv, mu_bv, mu_flexuralstrength]
        sigma_hv_bv_flexuralstrength = [sigma_hv, sigma_bv, sigma_flexuralstrength]
        CorrMat = [[1.0000, rho_hv_bv, rho_hv_flex],
                   [rho_hv_bv, 1.0000, rho_bv_flex],
                   [rho_hv_flex, rho_bv_flex, 1.0000]]

        sample_n = 2000 # number of samples in Monte-Carlo simulation
        random_hv, random_bv, random_flex = BreakageEvaluation.MvLogNRand( \
            mu_hv_bv_flexuralstrength, sigma_hv_bv_flexuralstrength, sample_n, CorrMat)

        random_InertiaI = np.pi*(random_bv*0.5)**4.0 / 4.0 ;

        Interval_K = evaluationresults_vars_dict['evaluation_interval']['widget_component'].value

        with evaluation_output:
            evaluation_output.clear_output()

            ## load trained ANN model parameters
            os.chdir(GUI_src_dir)

            model_parameters = BendingStress_ANNModel.load_ANN_model_parameters()

            # gather all vars_dict and prepare project_dir
            all_vars_dict, project_dir = prepareinputs_4_evaluation()

            if all_vars_dict['xzb_production']['value']==1: # gui create x
                xloc = np.arange(0, all_vars_dict['Lx']['value'][0]+all_vars_dict['dx']['value'][0],
                    all_vars_dict['dx']['value'][0], dtype=np.float64)
            elif all_vars_dict['xzb_production']['value']==2: # user uploaded:
                xloc=all_vars_dict['upload_x_zb']['value']['x'].values

            # run 1 time to get Hrms at beginning of veg zone
            breakage_fraction = np.zeros(len(xloc)) # no breakage
            Hrms_all, ESH_all, veg_Cd_all, xloc_all, zb_all = \
                makeinfle_run_load_extract(all_vars_dict, project_dir, breakage_fraction)

            print('Evaluating. Please wait ...')
            for irun in range(np.shape(Hrms_all)[0]): # loop for different incident wave conditions
                print('irun = ' + str(irun+1) )

                # evaluate breakage and get wave height decay
                # start from the 2nd veg grid. At grid(ix), ix = vegid[0]+1 ... end
                xloc = xloc_all[irun]
                zb = zb_all[irun]

                vegid = np.where( (xloc>=all_vars_dict['vegfrom']['value'][0]) &
                    (xloc<=all_vars_dict['vegto']['value'][0]) )
                veg_Cd = veg_Cd_all[irun][vegid[0][0]]

                Tp = all_vars_dict['Tp']['value'][irun]
                omega = 2.0*np.pi / Tp

                # print('Interval_K = ' + str(Interval_K))
                # print('vegid[0][0] = ' + str(vegid[0][0]))
                # print('vegid[0][-1]+1 = ' + str(vegid[0][-1]+1))

                for ix in range(vegid[0][0], vegid[0][-1]+1, Interval_K):

                    Hrms_previous = Hrms_all[irun][ix]
                    H10th_previous = 1.27*Hrms_previous*np.sqrt(2.0)
                    waterdepth = zb[ix]
                    wavenumber = dispersion(waterdepth, Tp)
                    kh = wavenumber*waterdepth

                    ESH_previous = ESH_all[irun][ix] / zb[ix]
                    random_effhv = random_hv * (ESH_previous/veg_hv)

                    bendingstressSF2LWTratio_ann = []
                    bendingstress_LWT = []
                    for isample in range(sample_n):
                        anninput = [random_effhv[isample]/waterdepth, H10th_previous/waterdepth, kh]
                        bendingstressSF2LWTratio_ann.append(BendingStress_ANNModel.trained_model(anninput, model_parameters))

                        bendingmoment_LWT = 0.5*rhowater*veg_Cd*random_bv[isample] * \
                            (omega*H10th_previous/2.0/np.sinh(kh))**2.0 * \
                            func_bendingmoment_LWT(min(random_effhv[isample], waterdepth), wavenumber)
                        bendingstress_LWT_tmp = Ac *  bendingmoment_LWT * random_bv[isample] \
                            / 2.0 / random_InertiaI[isample]
                        bendingstress_LWT.append(bendingstress_LWT_tmp)

                    bendingstress_SF    = np.array(bendingstress_LWT) * np.array(bendingstressSF2LWTratio_ann)
                    bendingstress_SF2flexuralstrength = np.array(bendingstress_SF) / np.array(random_flex)

                    # - get breakage fraction
                    breakage_number = (bendingstress_SF2flexuralstrength>=1.0).sum()
                    breakage_fraction[ix:ix+Interval_K] = breakage_number / sample_n

                    # make new infile, run, and load cshore results
                    Hrms_all, ESH_all, _, _, _ = \
                        makeinfle_run_load_extract(all_vars_dict, project_dir, breakage_fraction)

                # save evaluation results
                df = pd.DataFrame()

                df['x']=pd.Series(xloc_all[irun][:])
                df['Hrms']=pd.Series(Hrms_all[irun][:])
                df['breakage_fraction']=pd.Series(breakage_fraction)

                csvname = 'EvaluationResults' + '_Run' + str(irun+1) + '.csv'
                # print(csvname)
                savefile_path = os.path.join(project_dir, csvname)

                # remove file if existing already
                for fname in os.listdir(project_dir):
                    if fname == csvname:
                        os.remove(savefile_path)

                # write the header
                header1 = 'x (m), Hrms (m), Breakage Fraction\n'
                with open(savefile_path, 'w') as fp:
                    fp.write(header1)

                # write the rest
                df.to_csv(savefile_path, header= False,mode='a',index=False)

                print('Evaluation results, including x, Hrms and breakage fraction, are saved to: \n' + str(savefile_path))

            print('Evaluation is completed.')


            # plot evaluation results
            irun_4_plot = 0
            fig11 = plt.figure(evaluationresults_vars_dict['evaluation']['figname'], figsize=(6,4), dpi=80)

            ax88 = fig11.add_subplot(211)
            ax88.clear()
            ax88.plot(xloc_all[irun_4_plot], Hrms_all[irun_4_plot])
            ax88.set_ylabel('Hrms (m)', fontsize="large", fontweight="bold")
            ax88.set_xlabel('x (m)', fontsize="large", fontweight="bold")
            ax88.set_ylim([0, max(Hrms_all[irun_4_plot])*1.05])
            ax88.set_title('Run ' + str(irun_4_plot+1))

            ax99 = ax88.twinx()
            ax99.plot(xloc_all[irun_4_plot], -zb_all[irun_4_plot], 'k')
            # vegzb = np.array([0]*len(zb))
            # vegzb[vegid[0]] = 0.1
            # ax99.fill(xloc, -zb, xloc, -zb+vegzb, 'g')
            ax99.set_ylabel('elevation (m)', fontsize="large", fontweight="bold")
            ax99.set_ylim([-(max(zb_all[irun_4_plot])), max(-min(zb_all[irun_4_plot]), max(Hrms_all[irun_4_plot])*3)])

            ax77 = fig11.add_subplot(212)
            ax77.clear()
            ax77.plot(xloc_all[irun_4_plot], breakage_fraction)
            ax77.set_ylabel('breakage fraction', fontsize="large", fontweight="bold")
            ax77.set_xlabel('x (m)', fontsize="large", fontweight="bold")
    else:
        print('must be single wave condition.')


################## output handler ENDS ##################

################## generate_widgets  ####################
def generate_IntSlider(var_dict, key):
    slider = widgets.IntSlider(
        value=1, min=1, max=10,
        continuous_update=False,
        # layout=Layout(width='200px')
        )
    # slider.style.handle_color="navy"
    return slider

def generate_Dropdown(var_dict, key):
    drop_down = widgets.Dropdown(
        options=var_dict[key]['value_traits']['options'],
        value=var_dict[key]['value'],
        disabled=False, layout={'width': 'max-content'})
    return drop_down

def generate_BoundedFloatText(var_dict, key):
    bounded_float_text = widgets.BoundedFloatText(
        value=var_dict[key]['value'][0],
        min=var_dict[key]['value_traits']['min'],
        max=var_dict[key]['value_traits']['max'],
        step=0.1, disabled=False,
        layout={'width': 'max-content'})
    return bounded_float_text


def generate_BoundedIntText(var_dict, key):
    bounded_int_text = widgets.BoundedIntText(
        value=var_dict[key]['value'][0],
        min=var_dict[key]['value_traits']['min'],
        max=var_dict[key]['value_traits']['max'],
        step=1, disabled=False,
        layout={'width': 'max-content'})
    return bounded_int_text

def generate_FileUpload(var_dict, key):
    fl_load = widgets.FileUpload(
    accept='.csv',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
    multiple=False  # True to accept multiple files upload else False
    )
    return fl_load

def generate_savedata_Button(var_dict, key):
    butn = widgets.Button(
        description=var_dict[key].get('label_on_button'),
        button_style='success',
        icon='fa-floppy-o',
        style = {'button_color': 'maroon'},
        layout=Layout(width='auto') )
        # button icon: https://fontawesome.com/v4.7/icons/
    return butn

def generate_plot_Button(var_dict, key):
    butn = widgets.Button(
        description=var_dict[key].get('label_on_button'),
        button_style='success',
        # icon='fa-picture-o',
        icon='fa-area-chart',
        style = {'button_color': '#1f77b4'},
        layout=Layout(width='auto') )
        # button icon: https://fontawesome.com/v4.7/icons/
    return butn


def generate_makeinfile_Button(var_dict, key):
    butn = widgets.Button(
        description = '       Click me',
        button_style='success',
        icon='file-text-o',
        style = {'button_color': 'maroon'},
        layout=Layout(width='auto'))
        # button icon: https://fontawesome.com/v4.7/icons/
    butn.on_click(makeinfile_handler)
    return butn

def generate_runcshore_Button(var_dict, key):
    butn = widgets.Button(
        description='  Click me & wait',
        button_style='success',
        icon='fa-power-off',
        layout=Layout(width='auto'))
        # button icon: https://fontawesome.com/v4.7/icons/
    butn.on_click(runcshore_handler)
    return butn

def generate_evaluation_Button(var_dict, key):
    butn = widgets.Button(
        description=var_dict[key].get('label_on_button'),
        button_style='success',
        icon='fa-power-off',
        layout=Layout(width='auto'))
        # button icon: https://fontawesome.com/v4.7/icons/
    return butn

def generate_clearplot_Button(var_dict, key):
    butn = widgets.Button(
        # description='clear plots',
        description = var_dict[key].get('label_on_button'),
        # button_style='success',
        icon='fa-ban',
        # style = {'button_color': '#D3D3D3'},
        style = {'button_color': 'transparent'},
        layout=Layout(width='auto'))

    butn.layout.border="2px solid lightgray"
        # button icon: https://fontawesome.com/v4.7/icons/
    return butn

def generate_getvar_Button(var_dict, key):
    butn = widgets.Button(
        # description='clear plots',
        description = var_dict[key].get('label_on_button'),
        button_style='success',
        icon='fa-database',
        style = {'button_color': 'navy'},
        layout=Layout(width='auto'))
        # button icon: https://fontawesome.com/v4.7/icons/
    return butn

def generate_FloatText(var_dict, key):
    float_text = widgets.FloatText(value=var_dict[key]['value'][0],
        description='', disabled=False,
        step=0.1, layout={'width': 'max-content'})
    return float_text

def generate_IntText(var_dict, key):
    int_text = widgets.IntText(value=var_dict[key]['value'][0],
        description='', disabled=False,
        step=1, layout={'width': 'max-content'})
    return int_text

def generate_FileChooser(var_dict, key):
  # Create and display a FileChooser widget
  fc = FileChooser()

  fc.show_hidden = False # Change hidden files
  fc.rows = 10
  fc.use_dir_icons = True # Show or hide folder icons
  fc.show_only_dirs = False # Switch to folder-only mode
  fc.title = '' # Change the title (use '' to hide)

  # print(key)
  if key == 'select_project_dir':
      fc.register_callback(change_file_choser_title_project) # Register callback function
  elif key == 'cshoreexe':
      fc.register_callback(change_file_choser_title_cshoreexe) # Register callback function

  return fc

################## generate_widgets ENDS ####################

################## var_dict  ####################
physical_processes_vars_dict =  {
    'iline':{'value': [1], 'label': 'Number of  transects:',
             'tooltip': 'NLINES', 'widget_func': generate_BoundedIntText,
             'info_button': False,
             'value_traits': {'min': 1, 'max': 100, 'data_type': 'array',
                'update_int_float_value': True}},
    'isedav':{'value': 1, 'label': 'limited sediment availability:',
             'tooltip': 'ISEDAV', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'iperm':{'value': 0, 'label': 'Permeable bottom:',
             'tooltip': 'IPERM', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('no', 0)], # only allow impermeable bottom in this GUI
                'update_int_float_value': True},
            'output_widget': False},
    'iover':{'value': 1, 'label': 'Include overtopping:',
             'tooltip': 'IOVER', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'infilt':{'value': 0, 'label': 'Include infiltration landward of dune crest:',
             'tooltip': 'INFILT', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'iwtran':{'value': 0, 'label': 'Include wave transmission in landward wet zone:',
             'tooltip': 'IWTRAN', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'ipond':{'value': 0, 'label': 'Consider water ponding and runnel drainage:',
             'tooltip': 'IPOND', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'iwcint':{'value': 0, 'label': 'Consider wave-current interaction:',
             'tooltip': 'IWCINT', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'iroll':{'value': 1, 'label': 'Consider roller effects in the wet zone:',
             'tooltip': 'IROLL', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'iwind':{'value': 0, 'label': 'Consider wind effect:',
             'tooltip': 'IWIND', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'itide':{'value': 0, 'label': 'Consider tidal effect:',
             'tooltip': 'ITIDE', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'ibreaking':{'value': 1, 'label': 'Consider wave breaking:',
             'tooltip': 'users can turn off wave breaking for non-breaking waves.',
             'widget_func': generate_Dropdown,
             'info_button': True,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1), ('no', 0)],
                'update_int_float_value': True},
            'output_widget': False},
    'gamma':{'value': [0.7], 'label': 'Wave breaking ratio $\gamma = H_{rms}/h$:',
        'widget_func': generate_FloatText,
        'tooltip': '$\gamma$: 0.5 ~ 1.0',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}},
    'rwh':{'value': [0.02], 'label': 'Runup wire height $\delta_r$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': '0.01 m (small-scale experiments) ~ 0.1 m (prototype beaches and structures)',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
    'bottom_fric_factor':{'value': [0.015], 'label': 'Bottom friction factor $f_b$:',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
    'ilab':{'value': 1, 'label': 'Synced wave and water level conditions:',
             'tooltip': 'IPERM', 'widget_func': generate_Dropdown,
             'info_button': False,
             'value_traits': {'data_type': 'integer',
                'options': [('yes', 1)], # only allow synced wave para. and water level.
                'update_int_float_value': True},
            'output_widget': False},
}

morphology_vars_dict =  {
    'iprofl':{'value': 0, 'label': 'Movable bed:',
		'widget_func': generate_Dropdown,
        'tooltip': 'choose bed type',
        'info_button': False,
        'value_traits': {'data_type': 'integer',
            'options': [('--', 999), ('no', 0),('yes', 1),
                ('yes & no initial bottom smoothing', 1.1), ('dike erosion', 2)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': morpho_handler,
        'dependon': ['none'],},
    'np':{'value': [0.4], 'label': 'Sediment porosity:',
        'tooltip': '$n_p$ default: 0.4',
        'widget_func': generate_FloatText,
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['none']},
    'd50':{'value': [0.3], 'label': 'Median grain size (mm):',
        'tooltip': 'd50',
        'widget_func': generate_FloatText,
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['none']},
    'wf':{'value': [0.03], 'label': 'Sediment fall velocity (m/s):',
        'tooltip': '$w_f$',
        'widget_func': generate_FloatText,
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'effb':{'value': [0.005], 'label': 'suspension efficiency due to breaking:',
        'tooltip': '$e_B$: 0.002~0.01 (typically 0.005)',
        'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'efff':{'value': [0.3], 'label': 'suspension efficiency due to friction:',
        'tooltip': '$e_f$: default 0.01',
        'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'slp':{'value': [0.3], 'label': 'suspended load parameter:',
        'tooltip': '$a$: 0.1~0.4 (typically 0.2)',
        'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'slpot':{'value': [0.3], 'label': 'overtopping suspended load parameter:',
        'tooltip': '$a_0$: 0.1~3.6 (typically 0.5)', 'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'tanphi':{'value': [0.36], 'label': 'tangent (sediment friction angle):',
        'tooltip': '$tan\phi$: typically 0.36', 'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'blp':{'value': [0.04], 'label': 'bedload parameter:',
        'tooltip': '$b$: 0.001~0.004 (typically 0.002)', 'widget_func': generate_FloatText,
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    'gravitys':{'value': [2.65], 'label': 'Sediment specific gravity:',
        'tooltip': '$s$', 'widget_func': generate_FloatText,
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['movablebed']},
    }


vegetation_vars_dict = {
    'iveg':{'value': 999, 'label': 'Vegetation module:',
        'widget_func': generate_Dropdown,
        'widget_component': widgets.Label(),
        'tooltip': 'turn on/off vegetation module',
        'info_button': False,
        'value_traits': {'data_type': 'integer',
            'options': [('--', 999), ('No vegetation', 0), ('With vegetation', 3)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': vegonoff_handler,
        'dependon': ['none'],
    },
    'vegNsegtype':{'value': 1,
        'label': 'Consider spatially-varing $C_D$:', #'# of segments in veg. field:',
        'widget_func': generate_Dropdown,
        'widget_component': widgets.Label(),
        'tooltip': 'with spatially-varing Cd, veg zone is divided into N segments and N=(width of veg zone)/10.',
        'info_button': True,
        'value_traits': {'data_type': 'integer',
            'options': [('no', 1), ('yes', 2)],
            'update_int_float_value': True},
        'dependon': ['vegon'],
        },
    'idiss':{'value': 1,
        'label': 'Wave energy dissipation:',
        'widget_func': generate_Dropdown,
        'widget_component': widgets.Label(),
        'tooltip': 'vegetation-induced energy dissipation rate',
        'info_button': True,
        'value_traits': {'data_type': 'integer',
            'options': [('Mendez and Losada (2004)', 1), ('Chen and Zhao (2012)', 2)],
            'update_int_float_value': True},
        'output_widget': False,
        'dependon': ['vegon']},
    'iFv':{'value': 2,
        'label': 'Vegetal flow-drag resistence:',
        'widget_func': generate_Dropdown,
        'widget_component': widgets.Label(),
        'tooltip': 'phase-averaged depth-integrated vegetal drag',
        'info_button': True,
        'value_traits': {'data_type': 'integer',
            'options': [('parametric model (Zhu and Chen 2019)', 2), ('Hybrid model', 3), ('LWT model', 4)],
            'update_int_float_value': True},
        'output_widget': False,
        'dependon': ['vegon']},
    'vegfrom':{'value': [0.0],
         'label': 'Veg. from x = (m): ',
        'widget_func': generate_FloatText,
        'widget_component': widgets.Label(),
        'tooltip': 'veg. start from',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['vegon']},
    'vegto':{'value': [0.0],
         'label': 'Veg. to x = (m): ',
        'widget_func': generate_FloatText,
        'widget_component': widgets.Label(),
        'tooltip': 'veg. ends',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['vegon']},
    'vegrod':{'value': [0.3],
        'label': 'vegetation erosion limit (m): ',
        'widget_func': generate_FloatText,
        'tooltip': 'vegitation erosion limit below sand for failure',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['vegon']},
   'ivegtype_sptraits':{'value': 999, 'label': 'Vegetation type:',
        'widget_func': generate_Dropdown,
        'tooltip': 'treat vegetation as flexible/rigid',
        'info_button': False,
        'value_traits': {'data_type': 'integer',
            'options': [('--', 999),
                ('Flexible w/ spatially constant traits', 1),
                ('Flexible w/ spatially varying traits', 2),
                ('Rigid w/ spatially constant traits', 3),
                ('Rigid w/ spatially varying traits', 4)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': vegtype_sptraits_handler,
        'dependon': ['vegon']},
    'Nvstemflex':{'value': [400.0],
        'label': 'Stem population density (stems/$m^2$): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem population density (stems/$m^2$)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'Nvstemrigid':{'value': [400.0],
        'label': 'Stem population density (stems/$m^2$): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem population density (stems/$m^2$)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['rigidconstant']},
    'hvstemflex':{'value': [0.2],
        'label': 'Stem height (m): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem height (m)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'hvstemrigid':{'value': [0.2],
        'label': 'Stem height (m): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem height (m)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['rigidconstant']},
    'bvstemflex':{'value': [0.008],
        'label': 'Stem frontal width (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem frontal width (m):',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'bvstemrigid':{'value': [0.008],
        'label': 'Stem frontal width (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'Stem frontal width (m):',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['rigidconstant']},
    'Estem':{'value': [80000000.0],
        'label': "Young’s modulus of stem (N/$m^2$):",
        'widget_func': generate_FloatText,
        'tooltip': "Young’s modulus (N/$m^2$):",
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'Nvblade':{'value': [2000.0],
        'label': 'Blade population density (stems/$m^2$): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Blade population density (stems/$m^2$)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'hvblade':{'value': [0.7],
        'label': 'Blade height (m): ',
        'widget_func': generate_FloatText,
        'tooltip': 'Blade height (m)',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'bvblade':{'value': [0.005],
        'label': 'Blade width (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'Blade width (m):',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'tvblade':{'value': [0.00000685],
        'label': 'Blade thickness (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'Blade thickness (m):',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'Eblade':{'value': [80000000.0],
        'label': "Young’s modulus of blade (N/$m^2$):",
        'widget_func': generate_FloatText,
        'tooltip': "Young’s modulus of blade (N/$m^2$):",
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['flexconstant']},
    'flex_upload':{'value': 0,
        'label': "upload veg. trait file:",
        'label_on_button': 'veg. trait file',
        'widget_func': generate_FileUpload,
        'tooltip': "csv file only. columns in sequence of [x, Nvs, hvs, bvs, Es, Nvb, hvb, bvb, tvb, Eb]",
        'info_button': True,
        'value_traits': {'data_type': 'integer',
            'update_uploadfile': True,
            'filecontent_str': ['x', 'Nvs', 'hvs', 'bvs', 'Es', 'Nvb', 'hvb', 'bvb', 'tvb', 'Eb']},
        'dependon': ['flexvarying']},
    'rigid_upload':{'value': 0,
         'label': "upload veg. trait file:",
        'label_on_button': 'veg. trait file',
         'widget_func': generate_FileUpload,
         'tooltip': "csv file only. columns in sequence of [x, Nvs, hvs, bvs]",
         'info_button': True,
         'value_traits': {'data_type': 'integer',
            'update_uploadfile': True,
            'filecontent_str': ['x', 'Nvs', 'hvs', 'bvs'] },
         'dependon': ['rigidvarying']},
    'ivegCd':{'value':  1,
        'label': 'Vegetal drag coefficient $C_D$:',
        'widget_func': generate_Dropdown,
        'widget_component': widgets.Label(),
        'tooltip': 'default: Zhu et al. (2021) for flexible veg. & Anderson and Smith (2014) for rigid veg.',
        'info_button': True,
        'value_traits': {'data_type': 'integer',
            'options': [('default', 1), ('Jadhav et al. (2013)', 2), ('Moller et al. (2014)', 4),
                ('User-defined constant Cd', 0)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': Cd_handler,
        'dependon': ['vegon']},
    'CdCap':{'value': [45.0], 'label': 'maximum allowed $C_D$:',
        'widget_func': generate_FloatText,
        'widget_component': widgets.Label(),
        'tooltip': 'maximum allowed vegetal drag coefficient',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['vegon']},
    'userCd':{'value': [1.0],
        'label': "spatially constant $C_D$",
        'widget_func': generate_FloatText,
        'tooltip': "spatially constant $C_D$",
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['userCd']},
    }

boundary_vars_dict = {
   'simu_type':{'value': 999, 'label': 'Simulation type:',
        'widget_func': generate_Dropdown,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'integer',
            'options': [('--', 999), ('constant boundary condition', 1), ('time-varying boundary conditions', 2)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': singlebatch_handler,
        'dependon': ['none'] },
    'Hrms':{'value': [0.23], 'label': 'Wave height $H_{rms}$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': '$H_{rms}=H_{m0}/\sqrt{2}$',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['singlerun'] },
    'Tp':{'value': [4.0], 'label': 'Wave period $T_{z}$ (s):',
        'widget_func': generate_FloatText,
        'tooltip': 'mean wave period (Tz = Tp/1.25 for JONSWAP-type spectra, Tp denotes peak period)',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['singlerun'] },
    'waveangle':{'value': [0.0], 'label': 'Wave angle relative to shore normal ($^{\circ}$):',
        'widget_func': generate_FloatText,
        'tooltip': '-80$^{\circ}$ ~ 80$^{\circ}$',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['singlerun'] },
    'wsetup':{'value': [0.0], 'label': 'Mean water level $\overline {\eta}$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'see sketch below',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['singlerun'] },
    'swlbc':{'value': [0.0], 'label': 'Still water level $S$ (m) above datum $z=0$:',
        'widget_func': generate_FloatText,
        'tooltip': 'see sketch below',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['singlerun'] },
    # 'totaltime':{'value': 1.0, 'label': 'Total time (s):',
    #     'widget_func': generate_FloatText,
    #     'tooltip': 'determined from time series of wave conditions',
    #     'info_button': True,
    #     'value_traits': {'data_type': 'float',
    #         'update_int_float_value': True},
    #     'dependon': ['batchrun'] },
    # 'dt':{'value': 1.0, 'label': 'Time interval (s):',
    #     'widget_func': generate_FloatText,
    #     'tooltip': 'determined from time series of wave conditions',
    #     'info_button': True,
    #     'value_traits': {'data_type': 'float',
    #         'update_int_float_value': True},
    #     'dependon': ['batchrun'] },
    'upload_wavefiles':{'value': 1.0,
        'label': 'Upload wave file:',
        'label_on_button': 'wave file',
        'widget_func': generate_FileUpload,
        'tooltip': 'csv file only. columns in sequence of [time (s), Hrms (m), Tz (s), wave angle (degree) & mean water level (m)].',
        'info_button': True,
        'value_traits': {'data_type': 'float',
            'update_uploadfile': True,
            'filecontent_str': ['time', 'Hrms', 'Tp', 'waveangle', 'meanwaterlevel', 'stillwaterlevel']},
        'dependon': ['batchrun'] },
    'plotBC':{'value': 1.0,
        'label': '',
        'label_on_button': 'plot wave conditions',
        'widget_func': generate_plot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['batchrun'],
        'figname': 'wave conditions' },
    'clear_plotBC':{'value': 1.0,
        'label': '',
        'label_on_button': 'clear wave plot',
        'widget_func': generate_clearplot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['batchrun']},
    }

bathymetry_vars_dict = {
   'xzb_production':{'value': 999, 'label': 'Grid & depth production method:',
        'widget_func': generate_Dropdown,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'integer',
            'options': [('--', 999), ('GUI generated', 1), ('User uploaded', 2)],
            'update_int_float_value': True},
        'output_widget': True,
        'dropdown_handler_func': xzbproduction_handler,
        'dependon': ['none'] },
    'Lx':{'value': [100.0], 'label': 'Transect length $L$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['guigenerated'] },
    'dx':{'value': [0.5], 'label': 'Grid size $\Delta x$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'CSHORE-VEG allows at most 5000 grids along the transect.',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['guigenerated'] },
    'toeloc':{'value': [70.0], 'label': 'Toe location (m):',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['guigenerated'] },
    'zboffshore':{'value': [-2.0], 'label': 'Offshore bottom elevation $z_{b0}$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['guigenerated'] },
    'zbonshore':{'value': [0.8], 'label': 'Onshore bottom elevation $z_{b1}$ (m):',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},
        'dependon': ['guigenerated'] },
    'upload_x_zb':{'value': 1.0,
        'label': 'Upload x & zb (m)',
        'label_on_button': 'x & zb file',
        'widget_func': generate_FileUpload,
        'tooltip': 'csv file only. columns in sequence of [grid, bed elevation] in meters.',
        'info_button': True,
        'value_traits': {'data_type': 'float',
            'update_uploadfile': True,
            'filecontent_str': ['x', 'zb']},
        'dependon': ['uploadxzb'] },
    'plotzb_gui':{'value': 1.0,
        'label': '',
        'label_on_button': ' Plot bathymetry',
        'widget_func': generate_plot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['guigenerated'],
        'figname': 'GUI generated bathymetry' },
    'clear_plotzb_gui':{'value': 1.0,
        'label': '',
        'label_on_button': '  Clear bathymetry plot',
        'widget_func': generate_clearplot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['guigenerated'] },
    'plotzb_upload':{'value': 1.0,
        'label': '',
        'label_on_button': '  Plot bathymetry',
        'widget_func': generate_plot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['uploadxzb'],
        'figname': 'User-uploaded bathymetry' },
    'clear_plotzb_upload':{'value': 1.0,
        'label': '',
        'label_on_button': '  Clear bathymetry plot',
        'widget_func': generate_clearplot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['uploadxzb'] },
}

select_projectdir_vars_dict = {
    'select_project_dir':{'value': '',
        'label': 'Select project folder:',
        'widget_func': generate_FileChooser,
        'tooltip': '',
        'info_button': False,
        'value_traits': {},
        'dependon': ['none']  },
}

makeinfile_vars_dict = {
    'makeinfile':{'value': 0.0,
        'label': '',
        'widget_func': generate_makeinfile_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {},
        'dependon': ['none'] },
}

cshoreexe_vars_dict = {
    'cshoreexe':{'value': '',
        'label': 'Select CSHORE-VEG executable:',
        'widget_func': generate_FileChooser,
        'tooltip': '',
        'info_button': False,
        'value_traits': {},
         },
}

runcshore_vars_dict = {
    'runcshore':{'value': 0.0,
        'label': '',
        'widget_func': generate_runcshore_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {},
        'dependon': ['none'] },
}

# "spatialplot_vars_dict" is actually for constant wave boundary condition
spatialplot_vars_dict = {
    'numofwaverun':{'value': [1],
        'label': 'Number of run:',
        'widget_func': generate_IntSlider,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': False},
        'dependon': ['none'] },
    'plottype':{'value': 'root-mean-square wave height',
        'label': '\xa0\xa0\xa0\xa0\xa0\xa0\xa0\xa0Spatial plot:',
        'tooltip': '',
        'widget_func': generate_Dropdown,
        'info_button': False,
        'value_traits': {'data_type': 'integer',
        'options': ['root-mean-square wave height', 'mean water level', 'drag coefficient', 'effective plant height'],
        'update_int_float_value': False},
        'dependon': ['none'],
        'figname': 'CSHORE-VEG outputs'},
    'clearplot':{'value': [1.0],
        'label': '',
        'label_on_button': '       Clear plots',
        'widget_func': generate_clearplot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'] },
    'savedata':{'value': [1.0],
        'label': '',
        'label_on_button': '      Save spatial data of selected run',
        'widget_func': generate_savedata_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'] },
}

# "temporalplot_vars_dict" is actually for time-varying wave boundary condition
temporalplot_vars_dict = {
    'gageloct':{'value': [0.0],
        'label': 'Gage location (m):',
        'widget_func': generate_FloatText,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'] },
    'varatgageloct':{'value': [1.0],
        'label': '',
        'label_on_button': '       Save data at gage location & save wave runup',
        'widget_func': generate_savedata_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'] },
}


vegstatistics_vars_dict =  {
    'Nv_in_evaluation':{'value': [400.0],
        'label': 'Constant $N_v$ (stems/$m^2$):',
        'widget_func': generate_FloatText,
        'tooltip': 'Constant population density (stems/m2).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'meanhv':{'value': [0.2],
        'label': 'MEAN ($h_v$) (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'mean stem height (m).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'stdhv':{'value': [0.05],
        'label': 'SD ($h_v$) (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'standard deviation of stem height (m).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'meanbv':{'value': [0.008],
        'label': 'MEAN ($b_v$) (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'mean stem width (m).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'stdbv':{'value': [0.001],
        'label': 'SD ($b_v$) (m):',
        'widget_func': generate_FloatText,
        'tooltip': 'standard deviation of stem width (m).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'meanflexuralstrength':{'value': [8000000],
        'label': 'MEAN ($\sigma_{flex}$) (Pa):',
        'widget_func': generate_FloatText,
        'tooltip': 'mean flexural strength (Pa).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'stdflexuralstrength':{'value': [4600000],
        'label': 'SD ($\sigma_{flex}$) (Pa):',
        'widget_func': generate_FloatText,
        'tooltip': 'standard deviation of flexural strength (Pa).',
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True}, },
    'corrcoeff_hv_bv':{'value': [0.2],
        'label': " $r(h_v, b_v)$:",
        'widget_func': generate_FloatText,
        'tooltip': "Correlation coefficient between stem height and stem width.",
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
    'corrcoeff_hv_flex':{'value': [-0.2],
        'label': "$r(h_v, \sigma_{flex})$:",
        'widget_func': generate_FloatText,
        'tooltip': "Correlation coefficient between stem height and flexural strength.",
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
    'corrcoeff_bv_flex':{'value': [-0.4],
        'label': "$r(b_v, \sigma_{flex})$:",
        'widget_func': generate_FloatText,
        'tooltip': "Correlation coefficient between stem width and flexural strength.",
        'info_button': True,
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
}

evaluationresults_vars_dict = {
    'evaluation_interval':{'value': [1],
        'label': "Evaluation interval:",
        'widget_func': generate_IntText,
        'tooltip': "perform stem breakage at an interval of K grids.",
        'info_button': True,
        'dependon': ['none'],
        'value_traits': {'data_type': 'array',
            'update_int_float_value': True},},
    'evaluation':{'value': [1.0],
        'label': '',
        'label_on_button': '       Evaluate stem breakage',
        'widget_func': generate_evaluation_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'],
        'figname': 'Evaluation Results'},
    'clearevaluationplot':{'value': [1.0],
        'label': '',
        'label_on_button': '       Clear evaluation plot',
        'widget_func': generate_clearplot_Button,
        'tooltip': '',
        'info_button': False,
        'value_traits': {'data_type': 'float'},
        'dependon': ['none'] },
}

################## var_dict ENDS  ####################

############## set vegetation accordion ##############
# output widgets for vegetation accordion
vegonoff_output = widgets.Output()
vegtype_sptraits_output = widgets.Output()
Cd_output = widgets.Output()

veg_onoff = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'none')

option_vegon = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'vegon')

option_flexconstant = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'flexconstant')

option_rigidconstant = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'rigidconstant')

option_flexvarying = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'flexvarying')

option_rigidvarying = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'rigidvarying')

option_userCd = generate_gridlayout_from_var_dict(
    vegetation_vars_dict, 'userCd')

for key in vegetation_vars_dict:
    if vegetation_vars_dict[key].get('output_widget'): # note: use dict.get(key) instead of dict[key] because key may not exist.
        vegetation_vars_dict[key]['widget_component'].observe(
            vegetation_vars_dict[key]['dropdown_handler_func'], 'value')

vegetation_accordion = widgets.VBox([
    veg_onoff,
    vegonoff_output,
    ])

############## set boundary accordion ##############
# output widgets for boundary accordion
singlebatch_output = widgets.Output()
hydroBCplot_output = widgets.Output()

file = open("./Auxiliary/Fig_bed_datum_SWL_MWL.png", "rb")
image = file.read()
datum_MWL_fig = widgets.Image(value=image, format='png', width=290, height=50)

simu_singlebatch = generate_gridlayout_from_var_dict(
    boundary_vars_dict, 'none')

option_single = widgets.VBox([generate_gridlayout_from_var_dict(
    boundary_vars_dict, 'singlerun'), datum_MWL_fig])

option_batch = widgets.VBox([generate_gridlayout_from_var_dict(
    boundary_vars_dict, 'batchrun'), datum_MWL_fig])

for key in boundary_vars_dict:
    if boundary_vars_dict[key].get('output_widget'): # note: use dict.get(key) instead of dict[key] because key may not exist.
        boundary_vars_dict[key]['widget_component'].observe(
            boundary_vars_dict[key]['dropdown_handler_func'], 'value')

boundary_vars_dict['plotBC']['widget_component'].on_click(plotBC_handler)
boundary_vars_dict['clear_plotBC']['widget_component'].on_click(clear_plotBC_handler)

boundary_accordion = widgets.VBox([
    simu_singlebatch,
    singlebatch_output,
    hydroBCplot_output
    ])

############## set bathymetry accordion ##############
# output widgets for bathymetry accordion
guiuser_output = widgets.Output()
zbplot_output = widgets.Output()

bathymetry_head = generate_gridlayout_from_var_dict(
    bathymetry_vars_dict, 'none')

option_uploadxzb = generate_gridlayout_from_var_dict(
    bathymetry_vars_dict, 'uploadxzb')

option_guigenerated = generate_gridlayout_from_var_dict(
    bathymetry_vars_dict, 'guigenerated')

for key in bathymetry_vars_dict:
    if bathymetry_vars_dict[key].get('output_widget'): # note: use dict.get(key) instead of dict[key] because key may not exist.
        bathymetry_vars_dict[key]['widget_component'].observe(
            bathymetry_vars_dict[key]['dropdown_handler_func'], 'value')

bathymetry_vars_dict['plotzb_gui']['widget_component'].on_click(plotzbgui_handler)
bathymetry_vars_dict['clear_plotzb_gui']['widget_component'].on_click(clear_plotzbgui_handler)

bathymetry_vars_dict['plotzb_upload']['widget_component'].on_click(plotzbupload_handler)
bathymetry_vars_dict['clear_plotzb_upload']['widget_component'].on_click(clear_plotzbupload_handler)

# the reason to have two handler for plotzb is:
# In order to update on existing figure instead of plotting extra figure, I
# used %matplotlib notebook and fixed figure ID. Thus, if I switch between "gui" and 'upload',
# the figure will be gone and will not show up. Note that using %matplotlib inline
# and varying figure ID (fig = plt.figure()), this won't be a problem. However, you will have
# one figure every time you hit the plotzb button.
# Eventually, I figured an easy solution is to use different handlers for different plot button.
# Also, I added a 'clear plot button' so that I can remove the plot and clear the figure for the next plot.


bathymetry_accordion = widgets.VBox([
    bathymetry_head,
    guiuser_output,
    zbplot_output
    ])

############## set morphology accordion ##############
# output widgets for morphology accordion
movable_output = widgets.Output()

morpho_head = generate_gridlayout_from_var_dict(
    morphology_vars_dict, 'none')

option_movablebed = generate_gridlayout_from_var_dict(
    morphology_vars_dict, 'movablebed')

for key in morphology_vars_dict:
    if morphology_vars_dict[key].get('output_widget'): # note: use dict.get(key) instead of dict[key] because key may not exist.
        morphology_vars_dict[key]['widget_component'].observe(
            morphology_vars_dict[key]['dropdown_handler_func'], 'value')

morpho_accordion = widgets.VBox([
    morpho_head,
    movable_output,
    ])

############## set makeinfile accordion ##############
# output widgets for makeinfile accordion
makeinfile_output = widgets.Output()

makeinfile_head = generate_gridlayout_from_var_dict(
    makeinfile_vars_dict, 'none')

makeinfile_accordion = widgets.VBox([
    makeinfile_head,
    makeinfile_output,
    ])


############## set runcshore accordion ##############
# output widgets for runcshore accordion
runcshore_output = widgets.Output()

runcshore_head = generate_gridlayout_from_var_dict(
    runcshore_vars_dict, 'none')

runcshore_accordion = widgets.VBox([
    runcshore_head,
    runcshore_output,
    ])

############## set spatialplot accordion ##############
# output widgets for spatialplot accordion
HrmsMWLplot_output = widgets.Output()
savespatialdata_output = widgets.Output()
getdata_output = widgets.Output() # this is the output for getdata() function.
R2p_output = widgets.Output()

spatialplot_head = generate_gridlayout_from_var_dict(
    spatialplot_vars_dict, 'none')

spatialplot_vars_dict['numofwaverun']['widget_component'].observe(update_plot, names="value")
spatialplot_vars_dict['plottype']['widget_component'].observe(update_plottype, 'value')
spatialplot_vars_dict['savedata']['widget_component'].on_click(savespatialdata_handler)
spatialplot_vars_dict['clearplot']['widget_component'].icon = 'fa-refresh'
spatialplot_vars_dict['clearplot']['widget_component'].on_click(clearplot_handler)

spatialplot_accordion = widgets.VBox([
    spatialplot_head,
    getdata_output,
    savespatialdata_output,
    R2p_output,
    HrmsMWLplot_output,
    ])

############## set temporalplot accordion ##############
# output widgets for temporalplot accordion
savetemporaldata_output = widgets.Output()

temporalplot_head = generate_gridlayout_from_var_dict(
    temporalplot_vars_dict, 'none')

temporalplot_vars_dict['varatgageloct']['widget_component'].on_click(savetemporaldata_handler)

temporalplot_accordion = widgets.VBox([
    temporalplot_head,
    savetemporaldata_output,
    ])

############## set evaluation accordion ##############
# output widgets for evaluation accordion
evaluation_output = widgets.Output()
saveevaluation_output = widgets.Output()

evaluationresults_head = generate_gridlayout_from_var_dict(
    evaluationresults_vars_dict, 'none')

evaluationresults_vars_dict['evaluation']['widget_component'].on_click(evaluate_veg_breakage)
evaluationresults_vars_dict['clearevaluationplot']['widget_component'].icon = 'fa-refresh'
evaluationresults_vars_dict['clearevaluationplot']['widget_component'].on_click(clearevaluationplot_handler)

# evaluationresults_vars_dict['saveevaluationresults']['widget_component'].on_click(saveevaluationdata_handler)

evaluationresults_accordion = widgets.VBox([
    evaluationresults_head,
    saveevaluation_output,
    evaluation_output,
    ])

###########################################
select_projectdir_accordion_dict = {
    'type': 'accordion',
    'display_name': 'PROJECT DIRECTORY',
    'child_dict': select_projectdir_vars_dict,
    'child': list(select_projectdir_vars_dict.keys()) }

physical_process_accordion_dict = {
    'type': 'accordion',
    'display_name': 'PHYSICAL PROCESSES',
    'child_dict': physical_processes_vars_dict,
    'child': list(physical_processes_vars_dict.keys()) }

morphology_accordion_dict = {
    'type': 'specialtreat',
    'display_name': 'MORPHOLOGICAL CHANGES',
    'child_dict': morphology_vars_dict,
    'child_widget': morpho_accordion,
    'child': list(morphology_vars_dict.keys()) }

vegetation_accordion_dict = {
    'type': 'specialtreat',
    'display_name': 'VEGETATION-RELATED PARAMETERS',
    'child_dict': vegetation_vars_dict,
    'child_widget':vegetation_accordion,
    'child': list(vegetation_vars_dict.keys()) }

boundary_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'SEAWARD BOUNDARY CONDITIONS',
    'child_dict': boundary_vars_dict,
    'child_widget': boundary_accordion,
    'child': list(boundary_vars_dict.keys()) }

bathymetry_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'BATHYMETRIC & COMPUTATIONAL GRIDS',
    'child_dict': bathymetry_vars_dict,
    'child_widget': bathymetry_accordion,
    'child': list(bathymetry_vars_dict.keys()) }

makeinfile_accordion_dict ={
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'GENERATE INPUT FILE',
    'child_dict': makeinfile_vars_dict,
    'child_widget': makeinfile_accordion,
    'child': list(makeinfile_vars_dict.keys())
    }

cshoreexe_accordion_dict = {
    'type': 'accordion', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'CSHORE-VEG EXECUTABLE',
    'child_dict': cshoreexe_vars_dict,
    'child': list(cshoreexe_vars_dict.keys())
    }

runcshore_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'RUN CSHORE-VEG',
    'child_dict': runcshore_vars_dict,
    'child_widget': runcshore_accordion,
    'child': list(runcshore_vars_dict.keys())
    }

spatialplot_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'SPATIAL PLOTS',
    'child_dict': spatialplot_vars_dict,
    'child_widget': spatialplot_accordion,
    'child': list(spatialplot_vars_dict.keys())
    }

temporalplot_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'OUTPUTS AT A FIXED LOCATION',
    'child_dict': temporalplot_vars_dict,
    'child_widget': temporalplot_accordion,
    'child': list(temporalplot_vars_dict.keys())
    }

vegstatistics_accordion_dict = {
    'type': 'accordion', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'VEGETATION STATISTICS',
    'child_dict': vegstatistics_vars_dict,
    'child': list(vegstatistics_vars_dict.keys())
    }

evaluationresults_accordion_dict = {
    'type': 'specialtreat', # for 'specialtreat', must give 'child_widget'.
    'display_name': 'VEGETATION BREAKAGE EVALUATION',
    'child_dict': evaluationresults_vars_dict,
    'child_widget': evaluationresults_accordion,
    'child': list(evaluationresults_vars_dict.keys())
    }


###########################################
input_tab_dict = {
    'type': 'Tab',
    'display_name': 'INPUT',
    'child_dict': [
        select_projectdir_accordion_dict,
        bathymetry_accordion_dict,
        boundary_accordion_dict,
        vegetation_accordion_dict,
        physical_process_accordion_dict,
        morphology_accordion_dict,
        makeinfile_accordion_dict,
        ]}

run_tab_dict = {
    'type': 'Tab',
    'display_name': 'RUN',
    'child_dict': [
        cshoreexe_accordion_dict,
        runcshore_accordion_dict,
        ]}

visualization_tab_dict = {
    'type': 'Tab',
    'display_name': 'OUTPUT VISUALIZATION',
    'child_dict': [
        spatialplot_accordion_dict,
        temporalplot_accordion_dict,
        ]}

evaluation_tab_dict = {
    'type': 'Tab',
    'display_name': 'BREAKAGE EVALUATION',
    'child_dict': [
        vegstatistics_accordion_dict,
        evaluationresults_accordion_dict,
        ]}

############## set input tab ##############
accord_child = []
for iaccord_dict in input_tab_dict['child_dict']:
    if iaccord_dict['type'] == 'accordion':
        accord_child.append(generate_accordion(iaccord_dict))
    elif iaccord_dict['type'] == 'specialtreat':
        accord_child.append(iaccord_dict.get('child_widget'))

input_tab = widgets.Accordion(children=accord_child)

ii=0
for accord_child_dict in input_tab_dict['child_dict']:
    input_tab.set_title(ii, accord_child_dict['display_name'])
    ii=ii+1
input_tab.selected_index = None

############## set run tab ##############
# run cshore
#https://geekflare.com/learn-python-subprocess/
accord_child = []
for iaccord_dict in run_tab_dict['child_dict']:
    if iaccord_dict['type'] == 'accordion':
        accord_child.append(generate_accordion(iaccord_dict))
    elif iaccord_dict['type'] == 'specialtreat':
        accord_child.append(iaccord_dict.get('child_widget'))

run_tab = widgets.Accordion(children=accord_child)

ii=0
for accord_child_dict in run_tab_dict['child_dict']:
    run_tab.set_title(ii, accord_child_dict['display_name'])
    ii=ii+1
run_tab.selected_index = None

############## set visualization tab ##############
accord_child = []
for iaccord_dict in visualization_tab_dict['child_dict']:
    if iaccord_dict['type'] == 'accordion':
        accord_child.append(generate_accordion(iaccord_dict))
    elif iaccord_dict['type'] == 'specialtreat':
        accord_child.append(iaccord_dict.get('child_widget'))

visualization_tab = widgets.Accordion(children=accord_child)

ii=0
for accord_child_dict in visualization_tab_dict['child_dict']:
    visualization_tab.set_title(ii, accord_child_dict['display_name'])
    ii=ii+1
visualization_tab.selected_index = None

############## set evaluation tab ##############
accord_child = []
for iaccord_dict in evaluation_tab_dict['child_dict']:
    if iaccord_dict['type'] == 'accordion':
        accord_child.append(generate_accordion(iaccord_dict))
    elif iaccord_dict['type'] == 'specialtreat':
        accord_child.append(iaccord_dict.get('child_widget'))

evaluation_tab = widgets.Accordion(children=accord_child)

ii=0
for accord_child_dict in evaluation_tab_dict['child_dict']:
    evaluation_tab.set_title(ii, accord_child_dict['display_name'])
    ii=ii+1
evaluation_tab.selected_index = None

# ################################################################################################
#
GUI_src_dir = os.getcwd()
tabs_contents = ['INPUT', 'RUN', 'OUTPUT VISUALIZATION', 'BREAKAGE EVALUATION']
tabs = [VBox(description=name) for name in tabs_contents]
mainwindow = Tab(children = tabs,layout=Layout(width='100%')) # '99.6%' works

for i in range(len(tabs_contents)):
    mainwindow.set_title(i, tabs_contents[i])
mainwindow.selected_index = 0
mainwindow.children[0].children = [input_tab]
mainwindow.children[1].children = [run_tab]
mainwindow.children[2].children = [visualization_tab]
mainwindow.children[3].children = [evaluation_tab]
display(mainwindow)
