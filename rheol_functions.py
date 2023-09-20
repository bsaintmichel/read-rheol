"""rheol_functions.py : A small module that easily extracts and plots rheology data. 

Included functions
------------------
* read_rheology : reads your Anton Paar, Malvern or TA file. 
NOTE: It is wise to include a few columns in your file, e.g. : time,
 stress and strain (oscillatory and shear), shear rate, frequency, G' and G'',
 torque, normal force and temperature. 
  * for Malvern, you can do this through the "NAVIER_DEFAULT" table template)
  * for Anton Paar, make sure you have "Interval Time" ("Intervalle Temps")

* list_steps : shows and lists the detected steps

* assign_steps : allows you to assign a type to a given step, e.g. 'frequency sweep', a 'flow curve', etc.
Useful to group plots later on.

* plot_flowcurve : plots a flow curve (and fits the flow curve to a HB model)

* plot_asweep : plots oscillatory, amplitude sweep data

* plot_fsweep : plots oscillatory, frequency sweep data

* plot_tsweep : plots oscillatory, time sweep data (i.e. same frequency and amplitude throughout, 
e.g. for temperature sweep or to monitor gelation)

NOTE: there is no way to know what was applied and what was measured using RheoCompass and Kinexus .csv data, so you will
have to remember it. 
"""


import io
import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import time



##############################################################################################
## TA FUNCTIONS --------------------------------------------------------------
ta_mapper = {'Step time' : 'time', 'Time' : 'time_global', 'Shear rate' : 'shearrate', 'Stress' : 'stress', 'Strain' : 'strain', 
             'Viscosity' : 'viscosity', 'Storage modulus' : 'gprime', 'Loss modulus' : 'gsecond', 'Frequency' : 'freq', 
                'Axial force': 'normalforce', 'Gap': 'gap', 'Temperature': 'temp', 'Torque': 'torque', 'Displacement':'angle'}

def _read_TA(file_url):

    # First, replace ',' by '.' if needed ...
    with open(file_url) as file:
        contents = ''.join(file.readlines()).replace(',', '.')
    with open(file_url, 'w') as file:
        file.write(contents)

    with open(file_url) as file:
        line = file.readline()
        meta = {'strs_factor':np.nan, 'strn_factor':np.nan, 'nforce_factor':1, 'step_name':[], 'nsteps':0}
        step = 0
        decimal = '.'
        all_data = []
        
        # For TA, we will do small PD DataFrames for each step then merge them
        while line:
            line = file.readline() 
            
            # Gather some constants
            if 'Stress constant' in line:
                value = line.split('\t')[1].split(' ')[0]
                meta['strs_factor'] = float(value)
            if 'Strain constant' in line:
                value = line.split('\t')[1].split(' ')[0]
                meta['strn_factor'] = float(value)  

            ### TODO : for cones / plates : check radius and get the conversion from 
            # normal force to normal stress
            if 'Geometry Type	Cone plate' in line:
                _, diam = file.readline(), file.readline().split('\t')[1].split(' ')[0] # Diam usually written in mm 
                meta['nforce_factor'] = 8/(np.pi*(float(diam)/1000)**2)
            elif 'Geometry Type	Plate plate' in line:
                _, diamline = file.readline(), file.readline().split('\t')[1].split(' ')[0]
                meta['nforce_factor'] = 16/(np.pi*(float(diam)/1000)**2)


            # Fetch step names
            if 'Procedure name' in line:
                while 'proceduresegments' not in line:
                    meta['step_name'].append(line.split('\t')[1].rstrip()) # First line is stupid ...
                    line = file.readline()

            # Handle the actual data
            if '[step]' in line:
                data = ''
                while line != '\n' and line != '':  # Gather actual data
                    data += line
                    line = file.readline()
                now_data = pd.read_table(io.StringIO(data), delimiter='\t', decimal=decimal, skip_blank_lines=True, skiprows=[0,1,3])
                now_data['step'] = step
                step += 1
                all_data.append(now_data)

    meta['nsteps'] = len(all_data)
    all_data = pd.concat(all_data)
    print(meta)

    # Make sure some columns are in the list of columns, because otherwise 
    # it is a pain in the ass to work with them ...
    enforced_vars = ['Oscillation stress', 'Oscillation strain', 
                        'Torque', 'Stress', 'Strain', 'Displacement',
                        'Shear rate', 'Axial force', 'Normal stress']
    for var in enforced_vars:
        if var not in all_data.columns:
            all_data[var] = np.nan
    
    return all_data, meta

def _format_TA(all_data, meta):
        
        for step in range(meta['nsteps']):
            this_step = all_data['step'] == step
            all_data.loc[this_step, 'name'] = meta['step_name'][step]
            all_data.loc[this_step, 'step'] = step

            # Trying to fill as many additional columns
            is_oscstress = np.any(np.isfinite(all_data.loc[this_step, 'Oscillation stress']))
            is_oscstrain = np.any(np.isfinite(all_data.loc[this_step, 'Oscillation strain']))
            is_torque = np.any(np.isfinite(all_data.loc[this_step, 'Torque']))
            is_stress = np.any(np.isfinite(all_data.loc[this_step, 'Stress']))
            is_strain = np.any(np.isfinite(all_data.loc[this_step, 'Strain']))
            is_displ  = np.any(np.isfinite(all_data.loc[this_step, 'Displacement']))
            is_shearrate  = np.any(np.isfinite(all_data.loc[this_step, 'Shear rate']))
            is_axialforce = np.any(np.isfinite(all_data.loc[this_step, 'Axial force']))
            is_normalstress = np.any(np.isfinite(all_data.loc[this_step, 'Normal stress']))

            # If oscillatory stuff, simplify columns
            if is_oscstress:
                all_data.loc[this_step,'Stress'] = all_data.loc[this_step,'Oscillation stress']
                all_data.loc[this_step,'Strain'] = all_data.loc[this_step,'Oscillation strain']
                all_data.loc[this_step,'Shear rate'] = all_data.loc[this_step,'Oscillation strain']* \
                                                      all_data.loc[this_step,'Frequency']*2*np.pi
           
            # If Torque / Stress are missing, fill with the other value
            if not is_torque and  is_stress:       
                all_data.loc[this_step,'Torque'] = all_data.loc[this_step,'Stress']/meta['strn_factor']
            elif not is_stress and is_torque:
                all_data.loc[this_step,'Stress'] = all_data.loc[this_step,'Torque']*meta['strs_factor']*1e-6 # Conversion constant in Nm/Pa but torque in txt file in µN.m ...

            # Do a bit the same with normal stress
            if not is_axialforce and is_normalstress:
                all_data.loc[this_step,'Axial force'] = all_data.loc[this_step,'Normal stress']/meta['nforce_factor']
            if not is_normalstress and is_axialforce:
                all_data.loc[this_step,'Normal stress'] = all_data.loc[this_step, 'Axial force']*meta['nforce_factor']    
            
            # If strain is missing but other things are available (REALLY ?!)
            no_strain = (not is_oscstrain) and (not is_strain)
            if no_strain and is_displ:
                all_data.loc[this_step,'Strain'] = (all_data.loc[this_step,'Displacement'] - all_data.loc[this_step,'Displacement'].iloc[0])*meta['strn_factor']
            elif no_strain and is_shearrate:
                strain_rebuilt = cumtrapz(x=all_data.loc[this_step,'Step time'], y=all_data.loc[this_step,'Shear rate'])
                all_data.loc[this_step,'Strain'] = np.insert(strain_rebuilt, 0, 0)
            elif no_strain:
                print(f'_format_TA > Cannot infer strain from step {step} : {meta["step_name"][step]} in TA file.')

            all_data.loc[this_step,'Strain'] = (all_data.loc[this_step,'Strain'] - all_data.loc[this_step,'Strain'].iloc[0])*100

        # Re-build global time scale 
        dt = np.array(all_data['Step time'].diff())
        dt[dt < 0] = 0
        dt[0] = 0
        all_data['Time'] = np.cumsum(dt)

        # Rename, add compatibility columns
        all_data = all_data.rename(columns=ta_mapper)
        all_data = all_data.drop(columns=['Tan(delta)', 'Oscillation stress', 'Oscillation strain', 'Oscillation strain rate'], errors='ignore')
        all_data['type'] = ''
        all_data['status'] = ''
        
        return all_data

###############################################################################################
## ANTON PAAR FUNCTIONS --------------------------------------------------------------

antonpaar_mapper_en = {'Point No.':'point', 'Time' : 'time_global', 'Interval Time' : 'time', 'Shear Rate' : 'shearrate', 'Shear Stress' : 'stress', 'Shear Strain' : 'strain', 'Frequency' : 'freq', 
                'Storage Modulus' : 'gprime', 'Loss Modulus' : 'gsecond', 'Normal Force' : 'normalforce', 'Torque' : 'torque', 'Status' : 'status', 'Temperature': 'temp'}
antonpaar_mapper_fr = {'Point No.':'point', 'Temps' : 'time_global', 'Intervalle Temps' : 'time', 'Gradient de Cisaillement' : 'shearrate', 'Contrainte de Cisaillement' : 'stress', 'Déformation de Cisaillement' : 'strain', 'Fréquence' : 'freq', 
                'Module de Stockage' : 'gprime', 'Module de Perte' : 'gsecond', 'Force Normale' : 'normalforce', 'Couple' : 'torque', 'Etat' : 'status', 'Température': 'temp',
                'Déformation de Cisaillement (pour déformation sinusoïdale)':'raw_oscstrain', 'Contrainte de Cisaillement (pour déformation sinusoïdale)':'raw_oscstress'}
                # TODO : check what is the name for 'Gap' in RheoCompass

def _read_antonpaar(file_url):
    with open(file_url, encoding='utf-16-le') as file:
        file_finished = False
        all_data = []
        step = 0
        name = ''

        # For Anton Paar, we build a list of DataFrames that we will concatenate later
        while not file_finished:

            line = file.readline()
            parts = line.strip().split('\t')

            if 'Résultat' in parts[0] or 'Result' in parts[0]:
                name = parts[1]
            elif 'Intervalle et points de données' in parts[0] or 'Interval and data points' in parts[0]:
                step = int(parts[1])
            elif 'Interval données' in parts[0] or 'Data interval' in parts[0]:
                 # Means we "prepare" the table
                 # that will be read by Pandas
                data = ''                                 # header line
                while line != '\n' and line != '':  # Gather actual data
                    data += line
                    line = file.readline()

                df = pd.read_table(io.StringIO(data.replace(',','.')), delimiter='\t', skiprows=[1,2])
                df['name'] = name
                df['step'] = step
                all_data.append(df)
            
            file_finished = (line == '')

        all_data = pd.concat(all_data).reset_index(drop=True)
        return all_data    

def _format_antonpaar(df):
    # Add step n° and global time. Drop Unnecessary columns
    df['type'] = ''
    
    if 'Intervalle Temps' in df.columns:
        df = df.rename(columns=antonpaar_mapper_fr)
    elif 'Interval Time' in df.columns:
        df = df.rename(columns=antonpaar_mapper_en)
    else:
        print('> _format_antonpaar : Looks like you are using a Funky Language that I don''t understand ! ')
        raise ValueError

    # Compute correct (unique) steps
    step_change = df['step'].diff() != 0
    step_list = np.cumsum(step_change)
    df['step'] = step_list

    # Remove weird values for "invalid" nos in 'Point No.' column. 
    # Use the empty values to detect raw LAOS data
    df['point'] = np.array([str(elem).rstrip('(invalide)').rstrip('(invalid)').rstrip() for elem in df['point']])
    empty_point = df['point'] == ''
    df.loc[empty_point, 'point'] = np.nan
    df.loc[~empty_point, 'point'] = np.array(df.loc[~empty_point, 'point']).astype(float)
    df['raw'] = df['point'].isna()
    
    # Propagate some info to RAW data points
    row_to_propagate = df.iloc[0]
    for idx, row in df.iterrows():
        if not row['raw']:
            row_to_propagate = row
        else:
            df.loc[idx, 'point'] = row_to_propagate['point']
            df.loc[idx, 'time_global'] = row_to_propagate['time_global']
            df.loc[idx, 'time'] = row_to_propagate['time']
            df.loc[idx, 'freq'] = row_to_propagate['freq']
            df.loc[idx, 'strain'] = row_to_propagate['strain']
            df.loc[idx, 'stress'] = row_to_propagate['stress']

    # Fix global time issues ...
    df = _fix_globaltime_antonpaar(df)

    # Convert strain to 1 instead of %
    df['strain'] = df['strain']/100

    # Drop unnecessary columns
    df = df.drop(columns=['Interval données:', 'Data interval:'], errors='ignore')

    return df

def fix_strain_antonpaar(df, steps=None):
    """
    Function that fixes the strain in Anton Paar data 

    ARGS : 
    - df [PANDAS.DATAFRAME] : your not-necessary sliced Malvern data 
    - steps [NONE or LIST] : the steps where you wish to fix your stress
    (default : all steps are fixed)

    OUTPUT :
    - df [PANDAS.DATAFRAME] : the corrected data (d'uh !)

    NOTE : 
    If you want to fix your strain "for good" (i.e. in the original DataFrame),
    do not slice it first with rh.slice and fix the "sliced" df, because
    obviously the original dataframe will not be fixed...
    """
    if steps is None:
        steps = np.unique(df['steps'])

    for s in steps:
        indexer = df['step'] == s
        tm, sr = df.loc[indexer, 'time'], df.loc[indexer, 'shearrate']
        df.loc[indexer, 'strain'] = np.concatenate(([0], cumtrapz(tm, sr)))
    return df

def _fix_globaltime_antonpaar(df):
    """ 
    Function that fixes 'time_global'. I thought that this quantity
    thing was really global, but it seems that for each 'procedure' (i.e. vertical tab
    in Rheocompass), it resets. So I handle that here.
    """

    indices, = np.where(df['time_global'].diff() < 0)
    t_refs = df['time_global'].iloc[indices-1]

    for idx, t_ref in zip(indices, t_refs):
        df.loc[idx:,'time_global'] += t_ref

    return df


###############################################################################################
## MALVERN FUNCTIONS -------------------------------------------------------------------------

malvern_mapper =  {'Action Name':'name', 'Time (action)(s)': 'time',  'Time (sequence)(s)': 'time_global', 'Shear rate(s-¹)': 'shearrate', 'Shear stress(Pa)': 'stress', 
                   'Shear strain(%)': 'strain', 'Frequency(Hz)': 'freq', 'Shear modulus (elastic component)(Pa)': 'gprime', 
                   'Shear modulus (viscous component)(Pa)': 'gsecond', 'Shear modulus (complex component)(Pa)': 'gstar', 
                   'Normal force(N)': 'normalforce', 'Gap(mm)': 'gap', 'Torque(N m)':'torque', 'Temperature(°C)':'temp', 'Angular displacement(rad)':'angle'}

# For Malvern, well there ain't much to do ... since the thing is a proper csv file

def _read_malvern(file_url, decim_sep=',', field_sep=';'):
    data = pd.read_csv(file_url, decimal=decim_sep, sep=field_sep, encoding='utf-8')
    return data

def _format_malvern(df):
    # Add step n° and global time. Drop Unnecessary columns
    n = len(df)
    new_step, action = np.zeros(n, dtype=int), df.iloc[0]['Action Name']
    for index in range(n):
        if df.iloc[index]['Action Name'] != action:
            new_step[index] = 1
            action = df.iloc[index]['Action Name']
    
    df['raw'] = False
    df['step'] = np.cumsum(new_step).astype(int)
    df['type'] = ''
    df['status'] = ['N/A']*n
    df['point'] = np.nan
    
    # Merge "complex" (oscillatory) and "normal" (flow) stress and strains
    osc_steps = ~df['Complex shear strain(%)'].isnull()
    df.loc[osc_steps, 'Shear strain(%)'] = df.loc[osc_steps, 'Complex shear strain(%)']
    df.loc[osc_steps, 'Shear stress(Pa)'] = df.loc[osc_steps, 'Complex shear stress(Pa)']
    
    # Rename columns, then drop some
    df = df.rename(columns=malvern_mapper)
    df = df.drop(columns=['Complex shear strain(%)', 'Complex shear stress(Pa)'])

    # Adjust strain to be in 1 instead of %
    df['strain'] = df['strain']/100

    # Manage LAOS
    df = _malvern_laos(df)

    return df

def _malvern_laos(df):
    """ 
    Bits of code that handle LAOS for Malvern data.
    The annoying thing is that by default only the Torque and the Angular Position
    are computed ... AND data points never really cover a full period (they do less
    sometimes). I try to fix these things
    
    """
    steps = np.unique(df['step'])

    # First identify raw LAOS data
    is_raw = np.isnan(df['stress'])
    df.loc[is_raw, 'raw'] = True
    not_raw = df['raw'] == False

    if np.all(df['raw'] == False):
        print('_malvern_laos > No LAOS data found')
        return df

    # Find stress / strain ratio, use it to define raw_oscstrain and raw_oscstress
    syn_laos = df[~df['raw']]
    strn_ratio = np.mean(syn_laos['strain']/syn_laos['angle'])
    strs_ratio = np.mean(syn_laos['stress']/syn_laos['torque'])
    df['raw_oscstrain'] = df['angle']*strn_ratio
    df['raw_oscstress'] = df['torque']*strs_ratio
    df.loc[not_raw, 'raw_oscstress'] = np.nan
    df.loc[not_raw, 'raw_oscstrain'] = np.nan

    print(f'_malvern_laos > Stress ratio {strs_ratio:.2e} (Nm/Pa), strain ratio {strn_ratio:.2e} (1/rad)')

    # Drop data points corresponding to less
    # than an oscillation cycle in non-raw data mode
    # df['keep'] = True
    # df['tnormed'] = df['time']*df['freq']
    # for step in steps:
    #     condition = (df['step'] == step) & (df['raw'] == False) & any(np.isfinite(df['freq']))
    #     is_oscstrain = np.any(np.isfinite(df.loc[condition, 'freq']))
    #     if is_oscstrain:
    #         tvals = np.array(df.loc[condition, 'tnormed'])
    #         keep = np.zeros_like(tvals)
    #         tref = tvals[0]
    #         for no, tval in enumerate(tvals):
    #             if tval > tref + 1:
    #                 tref = tval
    #                 keep[no] = True
    #             else:
    #                 keep[no] = False
    #         keep[0] = True
    #         df.loc[condition, 'keep'] = keep.astype(bool)
    # bin = df[~df['keep']]
    # if len(bin) > 0:
    #     print('_malvern_laos > Deleting data from oscillatory steps corresponding to less than one cycle ...')
    #     df = df.drop(index=bin.index)
    
    # Create indices for non-raw values
    df['point'] = np.nan
    for step in steps:
        condition = (df['step'] == step) & (df['raw'] == False)
        npts = np.sum(condition)
        df.loc[condition, 'point'] = np.arange(npts)

    # Propagate some info to RAW data points ... 
    grouper = df.groupby('step')
    df = grouper.apply(lambda x : x.sort_values('time')).droplevel(0).reset_index(drop=True)

    df = df[::-1]
    row_to_propagate = df.iloc[0]
    for idx, row in df.iterrows():
        if not row['raw']:
            row_to_propagate = row
            freq = np.squeeze(row['freq'])
        else:
            df.loc[idx, 'point'] = row_to_propagate['point']
            df.loc[idx, 'stress'] = row_to_propagate['stress']
            df.loc[idx, 'strain'] = row_to_propagate['strain']
            df.loc[idx, 'freq'] = freq

    return df[::-1]

###############################################################################################
## GENERAL FUNCTIONS -------------------------------------------------------------------------

def read_rheology(file_url): 
    """
    Global function to read rheology data

    ARGS : 
    - file_url [STRING] : path to file

    RETURNS : 
    - df [PANDAS.DATAFRAME] : formatted data to be further processed or plotted
    """
    with open(file_url) as file:
        header = file.readline(150)
        if 'Time (action)(s)' in header: # Malvern
            filetype = 'Malvern'
            data = _format_malvern(_read_malvern(file_url, decim_sep=',', field_sep=';'))
        elif 'Filename' in header:      # TA
            filetype = 'TA'

            data, meta = _read_TA(file_url)
            data = _format_TA(data, meta)
        else:
            with open(file_url, encoding='utf-16-le') as file16: # Anton Paar (?)
                header16 = file16.readline(150)
                if ('Test:' in header16) or ('Projet:' in header16) or ('Project:' in header16):
                    data = _format_antonpaar(_read_antonpaar(file_url))
                    filetype = 'Anton Paar'
                else:
                    data = None
                    filetype = 'Unknown'

    print(f'read_rheology > File type is {filetype}.')
    return data

def assign_steps(df, steps, steptypes):
    """
    Function that assigns the step type if you are not happy with what the 
    automatic programme has done

    ARGS : 
    - df [PANDAS.DATAFRAME] : your rheology data
    - steps [LIST] : list of steps n° you want to assign
    - steptypes [LIST of STR or STR] : list of step types you want to assign to the step n°s

    OUTPUT : 
    - df [PANDAS.DATAFRAME] : the assigned data
    """
    if np.isscalar(steptypes):
        steptypes = [steptypes]*len(steps)
    for no, step in enumerate(steps):
        indices = df['step'] == step
        steptype = steptypes[no]
        df.loc[indices, 'type'] = steptype
        print('assign_steps > step : ' + str(step) + ' | Assigning as : ' + steptype)

    return df

def list_steps(df):
    """
    Function that finds and lists the step number(s) corresponding to a step type. 
    Or lists them if you specify nothing as second argument

    ARGS : 
    - df [PANDAS.DATAFRAME] : your global rheology data
    - step_type
      * [STRING], has to be one of : 'amplitudesweep', 'freqsweep', 'flowcurve', 'creep', 
      'shearstartup', 'rest', 'preshear', 'timesweep'
      * [NONE]

    OUTPUT : 
    - NONE if step_type is None 
    - MATCHES if step_type is str : a list of steps that match your query
    """
    # Note : step_type : 
    steps = np.unique(df['step'])
    print('------------- Step list in DataFrame --------------------')
    for step in steps:
        dfnow = df[df['step'] == step]
        iend = dfnow['time_global'].last_valid_index()
        steptype = f'{dfnow.iloc[0]["type"]:>15}'
        stepname = f'{dfnow.iloc[0]["name"]:>30}'
        if iend is not None:
            stepduration = f'{dfnow.loc[iend,"time"]:>10.2f}'
            time_global = f'{dfnow.loc[iend, "time_global"]:>10.2f}'
            print('  * Step n°' + str(step) + ' \t : '+ stepname + '\t is a ' + steptype + ' | Duration : ' + str(stepduration) + ' s | Total time is : ' + str(time_global) + ' s')
        else:
            print('  * Step n°' + str(step) + ' \t : '+ stepname + '\t is a ' + steptype + ' | Is a bit mysterious')
    print('---------------------------------------------------------')
    return None

def slice(df, steps):
    """ 
    Function that returns a subset of your data 
    based on a list of steps.

    ARGS : 
    - df [PANDAS.DATAFRAME] : your rheology data
    - steps [LIST] or [INT] or [ST] : your steps of interest, defined by a list of step n°s, a single step n° or
    a step type

    OUTPUT :
    - sliced_df [PANDAS.DATAFRAME] : a (smaller) set of your rheology data containing
    what you want
    """
    
    # If steps is int or list, consider it is a list of indices
    # If step is str, consider we are slicing on the types
    if type(steps) == int:
        key = 'step'
        vals = [steps]
    elif type(steps) == str:
        key = 'type'    
        vals = [steps]
    elif hasattr(steps, '__iter__'):
        key = 'step'
        vals = steps

    sliced_df = pd.DataFrame()
    for val in vals:
        sliced_df = pd.concat((sliced_df, df[df[key] == val]))

    if sliced_df.empty:
        print(f'slice > Warning : empty dataframe with slicing steps {steps}')
    return sliced_df

#########################################################
## --------- COMPUTING SCIENTIFIC THINGS FUNCTION -----

def _fit_HB(rate, strs, fit_up_to=1e3, fit_from=1e-3):
    """
    Function that computes a fit to a Herschel_Bulkley (HB)
    law in log-log coordinates to try to extract the HB 
    parameters of flow curve data
    """
    rate, strs = rate.to_numpy(), strs.to_numpy()
    valid_shearrate = np.logical_and(rate > fit_from, rate < fit_up_to)
    rate, strs = rate[valid_shearrate], strs[valid_shearrate]
    strs0 = np.min(np.abs(strs))

    potential_ys, results = np.linspace(0.75*strs0, 1.25*strs0, 200), np.zeros((200,3))
    for idx, test_ys in enumerate(potential_ys):
        valid = strs - test_ys > 0
        x, y = np.log(rate[valid]), np.log(strs[valid] - test_ys)
        rs = stats.linregress(x,y)
        results[idx] = np.array([rs.slope, rs.intercept, rs.rvalue])
    best_fit = np.argmax(results[:,2])
    slope, K, ys = results[best_fit,0], np.exp(results[best_fit,1]), potential_ys[best_fit] 
    return ys, K, slope


##############################################################################
## -------- PLOTTING FUNCTIONS ###############################################
def darken(palette, factor=0.6):
    """
    Function that darkens a Bokeh palette or a single color
    by a given factor
    """
    isStr = False
    if type(palette) == str:
        isStr = True
        palette = [palette]
    
    rgb = [[int(c[1:3], 16),int(c[3:5], 16),int(c[5:7], 16)] for c in palette]
    dark_rgb = lambda factor : [[int(c[0]*factor), int(c[1]*factor), int(c[2]*factor)] for c in rgb]
    hexvals = ['#' +'{0:0{1}x}'.format(c[0]*256**2 + c[1]*256 + c[2],6) for c in dark_rgb(factor)]

    if isStr:
        return hexvals[0] # Returns a str if a str was input
    else:
        return hexvals # Returns a list if a list was input
 
def plot_flowcurve(df, fit_from=1e-3, fit_up_to=1e3):
    """ 
    Function to plot flow curve steps

    ARGS : 
    - df [PANDAS.DATAFRAME] : your sliced rheology data
    - fit_from [float, default 1e-3] : from what shear rate you fit your flow curve (with a Herschel-Bulkley law)
    - fit_up_to [float, default 1e3] : up to where you fit your flow curve (with the same HB law ...)

    OUTPUT :
    - FIT [Nx4 list] : the fit to your steps, with [step, sigma_Y, K, exponent]
    """   
    if len(df) > 0 and 'step' in df.keys():
        steps = np.unique(df['step'])
        fits_all = np.zeros((len(steps), 4))

        # Create and format figure
        fig, ax = plt.subplots()
        ax.set_xlabel('γ (1/s)'), ax.set_ylabel('τ (Pa)')
        ax.set_title('Flow Curve')
        cmap = cm.magma

        for no, step in enumerate(steps):
            # Check that step makes sense
            if len(steps) > 1:
                color = cmap(256*(step-min(steps)))/(max(steps)- min(steps))
            else:
                color = np.array([0,0,0,1])

            flow_curve = df[df['step'] == step]
            steptype = flow_curve.iloc[0]['type'].lower()
            if steptype not in ('flowcurve', 'flow curve'):
                print('plot_flowcurve >  Warning : Cannot confirm that step ' + str(step) + ' is a flow curve.')    
            ys, K, exponent = _fit_HB(flow_curve['shearrate'], flow_curve['stress'], fit_up_to=fit_up_to, fit_from=fit_from)
            fits_all[no,:] = [step, ys, K, exponent]

            # Plotting + Display + producing the fitted curve
            shearrates = flow_curve['shearrate']
            fitted_stress = ys + K*shearrates**exponent 
            ax.loglog(shearrates, fitted_stress, color=0.8*color)
            ax.loglog(shearrates,flow_curve['stress'], marker='s', color=0.8*color, 
                      markersize=3,markerfacecolor=color, label="Flow curve, step " + str(step))
            fig.legend()    
            print(f'plot_flowcurve > fit for step {step} : {ys:.2f} + {K:.2f} γ^({exponent:.2f})')
        
        plt.show()
        time.sleep(0.5)

        return fits_all
    else:
        print('plot_flowcurve > No flowcurve step found')
        return None

def plot_asweep(df, plot_stress=False):
    """ 
    Function to plot amplitude sweeps

    ARGS : 
    - df [PANDAS.DATAFRAME] : your sliced rheology data
    - steps [LIST] : your step(s) you want to plot (can be only one step)

    OUTPUT :
    - a figure (d'uh !)
    """
    if len(df) > 0 and 'step' in df.keys():
        # Create and format figure
        fig, ax = plt.subplots()
        ax.set_xlabel('γ (1)'), ax.set_ylabel("G', G'' (Pa)")
        ax.set_title('Amplitude Sweep')
        cmap = cm.magma
        steps = np.unique(df['step'])
        
        for no, step in enumerate(steps):
            # Check that step makes sense
            if len(steps) > 1:
                color = cmap(256*(step-min(steps)))/(max(steps)- min(steps))
            else:
                color = np.array([0,0,0,1])
            amp_sweep = df[df['step'] == step]
            steptype = amp_sweep.iloc[0]['type'].lower()

            if steptype not in('amplitudesweep', 'asweep', 'amplitude sweep', 'amp sweep', 'a sweep'):
                print('plot_asweep > Warning : Cannot confirm that step ' + str(step) + ' is an amplitude sweep')
            ax.loglog(amp_sweep['strain'], amp_sweep['gprime'] , markersize=3, marker='s', color=0.8*color, 
                      markerfacecolor=color, label="G' , step " + str(step))
            ax.loglog(amp_sweep['strain'], amp_sweep['gsecond'], markersize=3, marker='o', color=0.8*color, 
                      markerfacecolor='lightgray', label="G'' , step " + str(step))
            if plot_stress: ax.loglog(amp_sweep['strain'], amp_sweep['stress'], markersize=3, marker='^', color=color[no], 
                                      markerfacecolor=cmap[no], label='σ , step ' + str(step))
        
        fig.legend()    
        time.sleep(0.5)
        plt.show()

        return 0
    else:
        print('plot_asweep > No amplitude sweep step in the sliced dataset')
        return None


def plot_fsweep(df, plot_stress=False):
    """ 
    Function to plot frequency sweeps

    ARGS : 
    - df [PANDAS.DATAFRAME] : your sliced rheology data
    - plot_stress [BOOL] [OPTIONAL] : add the global oscillatory stress to the plot

    OUTPUT :
    - a figure (d'uh !)
    """
    if len(df) > 0 and 'step' in df.keys():
        # Create and format figure
        fig, ax = plt.subplots()
        ax.set_xlabel('f (Hz)'), ax.set_ylabel("G', G'' (Pa)")
        ax.set_title('Frequency Sweep')
        cmap = cm.magma
        steps = np.unique(df['step'])
        
        for no, step in enumerate(steps):
            # Check that step makes sense
            if len(steps) > 1:
                color = cmap(256*(step-min(steps)))/(max(steps)- min(steps))
            else:
                color = np.array([0,0,0,1])
            f_sweep = df[df['step'] == step]
            steptype = f_sweep.iloc[0]['type'].lower()

            if steptype not in('freqsweep', 'fsweep', 'frequency sweep', 'freq sweep', 'f sweep'):
                print('plot_fsweep > Warning : Cannot confirm that step ' + str(step) + ' is a frequency sweep')
            ax.loglog(f_sweep['freq'], f_sweep['gprime'], markersize=3, marker='s', color=0.8*color, 
                      markerfacecolor=color, label="G' , step " + str(step))
            ax.loglog(f_sweep['freq'], f_sweep['gsecond'], markersize=3, marker='o', color=0.8*color, 
                      markerfacecolor='lightgray', label="G'' , step " + str(step))
            if plot_stress: ax.loglog(f_sweep['freq'], f_sweep['stress'], marker='^', color=color[no], 
                                      markerfacecolor=cmap[no], markersize=3, label='σ , step ' + str(step))
        
        fig.legend()    
        time.sleep(0.5)
        plt.show()
        return 0
    else:
        print('plot_fsweep > No frequency sweep step in the sliced dataset')
        return None
    
def plot_tsweep(df, plot_stress=False):
    """ 
    Function to plot time sweeps

    ARGS : 
    - df [PANDAS.DATAFRAME] : your sliced rheology data
    - plot_stress [BOOL] [OPTIONAL] : add the global oscillatory stress to the plot

    OUTPUT :
    - a figure (d'uh !)
    """
    if len(df) > 0 and 'step' in df.keys():
        # Create and format figure
        fig, ax = plt.subplots()
        ax.set_xlabel('t (s)'), ax.set_ylabel("G', G'' (Pa)")
        ax.set_title('Time Sweep')
        cmap = cm.magma
        steps = np.unique(df['step'])
        
        for no, step in enumerate(steps):
            # Check that step makes sense
            if len(steps) > 1:
                color = cmap(256*(step-min(steps)))/(max(steps)- min(steps))
            else:
                color = np.array([0,0,0,1])
            t_sweep = df[df['step'] == step]
            steptype = t_sweep.iloc[0]['type'].lower()

            if steptype not in('timesweep', 'tsweep', 'time sweep', 't sweep'):
                print('plot_tsweep > Warning : Cannot confirm that step ' + str(step) + ' is a frequency sweep')
            ax.loglog(t_sweep['time'], t_sweep['gprime'], markersize=3, marker='s', color=0.8*color, 
                      markerfacecolor=color, label="G' , step " + str(step))
            ax.loglog(t_sweep['time'], t_sweep['gsecond'], marker='o', color=0.8*color, 
                      markerfacecolor='lightgray', markersize=3, label="G'' , step " + str(step))
            if plot_stress: ax.loglog(t_sweep['time'], t_sweep['stress'], marker='^', color=color[no], 
                                      markerfacecolor=cmap[no], markersize=3, label='σ , step ' + str(step))
        
        fig.legend()    
        time.sleep(0.5)
        plt.show()
        return 0
    else:
        print('plot_fsweep > No time sweep step in the sliced dataset')
        return None

### Fourier projection and reconstruction

def proj_fourier(time, signal, nmodes=10):
    """ 
    A function that projects your signal 
    on Fourier series. 

    ARGS
    ----
    * time [list, pandas.Series, ...]: whatever X axis we plot the signal against in the non-Fourier domain.

    * signal [same as `time`]: the signal you want to analyse

    RETURNS
    ---- 
    * proj [dict] : your projection coefficients along the `cos` and `sin` but also in terms of `amp`(litude) and `phs`(phase)
    
    NOTES 
    ------
    * I expect `time` to be normalised by the angular frequency omega (t from 0 to 2*pi for one period)
    * It is always better to have an exact number of periods for `time` and `signal`
        to not mess up the integration
    * The program does not deal with NaNs in your data series and will fail miserably.
    
    """
    proj = {'mode':np.arange(0,nmodes),
            'sin':np.zeros(nmodes), 'cos':np.zeros(nmodes), 
            'amp':np.zeros(nmodes), 'phs':np.zeros(nmodes)}
    
    # Recast Pandas into np.arrays
    time = np.array(time)
    signal = np.array(signal)
    
    # Add one extra sample at the end of the signal corresponding to the first sample (--> better precision on Fourier coeffs)
    dt = np.mean(np.diff(time)) # We need it to estimate how many periods we have
    time = np.append(time, (time[-1] + dt))
    signal = np.append(signal, signal[0])

    for mno in range(nmodes):
        proj['cos'][mno] = np.trapz(signal*np.cos(mno*time), x=time)*2/np.nanmax(time)
        proj['sin'][mno] = np.trapz(signal*np.sin(mno*time), x=time)*2/np.nanmax(time)

    # Note : for mno = 0 we have to divide over 2 again 
    # otherwise the projection will be double what we need
    
    proj['cos'][0]/= 2
    proj['sin'][0]/= 2
    proj['amp'] = np.sqrt(proj['cos']**2 + proj['sin']**2)
    proj['phs'] = np.arctan2(-proj['sin'], proj['cos']) # It works with a (-) ... so be it.
    return proj

def build_fourier(proj, time, nmodes=10):
    """ A reconstruction of synthetic signals 
    based on a time (or X) scale and fourier coeffs
    obtained from `proj_fourier()`. """
    rebuild = np.zeros(np.shape(time))
    for mno in range(nmodes):
        rebuild = rebuild + proj['amp'][mno]*np.cos(mno*time + proj['phs'][mno])
    return rebuild
