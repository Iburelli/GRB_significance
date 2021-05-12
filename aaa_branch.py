import os
from os.path import join, isfile
import ctools
import cscripts
import astropy.units as u
from astropy.time import Time
from astropy.io import fits
import argparse
import yaml
import numpy as np
import glob


parser = argparse.ArgumentParser(description='The significance of a GRB observation is computed at different IRFs according to a visibility table, created with runCatVisibility.py. A configuation YAML file is required, the output is saved as NPY binary file.')
parser.add_argument('-f', '--config', required=True, type=str, help='configuration yaml file')
# configuration file
cf = parser.parse_args().config
# load params configuration from cf
with open(cf) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)


# ----------------------------------------------------------------------------- catalog
if '$' in cfg['path']['xmlcatalog']:
    catalog = os.path.expandvars(cfg['path']['xmlcatalog'])
else:
    catalog = cfg['path']['xmlcatalog']

if cfg['path']['output'] == None:
    output = 'visibility_output.npy'
elif '$' in cfg['path']['output']:
    output = os.path.expandvars(cfg['path']['output'])
else:
    output = cfg['path']['output']

if cfg['path']['xmlfilename'] == None:
    runids = glob.glob(cfg['path']['xmlcatalog'] + '/**/*.xml', recursive=True)
    if len(runids) == 0:
        raise ValueError('No valid XML file found')

elif type(cfg['path']['xmlfilename']) == str:
    if not isfile(join(catalog, cfg['path']['xmlfilename'])):
            raise ValueError(f'Specified template does not exist in catalog')
    runids = [cfg['path']['xmlfilename']]
else:
    runids = cfg['path']['xmlfilename']
    for runid in runids:
        if not isfile(join(catalog, runid)):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
runids = sorted(runids)

#-----------------------------------------------------function to append new line to txt file
def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)




#--------------------------------------------------------------reading visibility table

data = np.load(output, allow_pickle=True, encoding='latin1', fix_imports=True).flat[0]
events = list(data.keys())
sites = list(data[events[0]].keys())

#-----------------------------------------------importing ctools variables

#generating random seeds
seeds = np.random.randint(1,1000)
#simulation variables
caldb = cfg['ctools']['caldb']
sim_rad = cfg['ctools']['rad']
fitmodel = cfg['ctools']['fitmodel']
offset = cfg['ctools']['offset']
radius = cfg['ctools']['radius']
step=cfg['ctools']['timestep']
#inmodel selection
for runid in runids:

    if type(cfg['path']['xmlfilename']) == str:
        inmodel = join(catalog, runid)

    elif cfg['path']['xmlfilename'] == None:
        inmodel = runid
#------------------------------------------------Running on EVENTS
    for event in events:

#checking for the event number into the tables
        name =f'{runid.replace(".xml", "")}'
        if name.endswith(str(event)) :
            print(f'\nProcessing {event}')
            #reading trigger time from fits file
            with fits.open(cfg['path']['catalog']+f'/{event}.fits') as hdul:
                    hdr = hdul[0].header
                    t_trigger = Time(hdr['GRBJD'] * u.day, format='jd')
                    trigger=float(t_trigger.jd)

#----------------------------------------------------Running on SITES

            for site in sites:
                 #if site == 'South':

                    print(f'\nProcessing site {site}')
    #----------------------------------------------------Running on NIGHTS
                    for night in data[event][site]:
                        print(f'\nProcessing {night}')

    #----------------------------------------------------Checking visibility at the site during a specific night
                        if type(data[event][site]) == float:
                            print(f'\tThis contains NaNs ---> the source is not observable due to daylight or moon.')
                        if data[event][site][night]['irfs']['zref'][0] == -9.0:
                            print(f'\tThis contains NaNs event---> the source is not observable at the site.')
                            data[event][site][night]['sigma_ON/OFF'] = -9.0

    #----------------------------------------------------Simulation and analysis of VISIBLE sources
                        else:
                        #Defining the start time of the observation (in seconds form trgger)
                            t_obs_start= data[event][site][night]['irfs']['start'][0]
                            t_obs_start=(t_obs_start - trigger)*86400
                            t_obs_stop = data[event][site][night]['irfs']['stop'][-1]
                            #print(t_obs_stop)
                            t_obs_stop = (t_obs_stop - trigger)*86400

                        #Parameter used as sel_t_max to determine the time needed for the source to reach 3 sigma
                            t0 = t_obs_start + 30.
                            #print (t0)
                            #print (t_obs_start)


                            for j in range(27):
                                print(f'\n\tInteval {j+1}, Simulation time {t0-t_obs_start}')
                                somma_on=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))
                                somma_off=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))
                                #print (f'\n\tIntegration time: {t0}')
                                for i in range(len(data[event][site][night]['irfs']['zref'])):

                                    #-----------------SIMULATION TIMES

                                    t_min = data[event][site][night]['irfs']["start"][i]
                                    t_max= data[event][site][night]['irfs']["stop"][i]
                                    #print(t_min,t_max)
                                    #print(t_obs_stop)
                                    #converting times from jd to seconds from trigger
                                    sim_t_min=(t_min - trigger)*86400
                                    sim_t_max=(t_max-trigger)*86400



                                    if sim_t_min > t0:
                                        somma_on[i] = 0.
                                        somma_off[i] = 0.
                                        break

                                    elif t0  > t_obs_stop:
                                        break

                                    elif sim_t_max < t0:
                                        sel_t_max = sim_t_max

                                    else:
                                        sel_t_max = t0

                                    # ----------- IRF definition
                                    delta_t_irf = sel_t_max - sim_t_min

                                    if delta_t_irf < 94.9 * 60:
                                        irf_duration = '0.5h'
                                    elif 94.9 * 60 < delta_t_irf < 15.8 * 60 * 60:
                                        irf_duration = '5h'
                                    elif delta_t_irf > 15.8 * 60 * 60:
                                        irf_duration = '50h'

                                    name_irf = (f'{site}_z{data[event][site][night]["irfs"]["zref"][i]}_{irf_duration}')



                                    # ---------------selection of e_min and e_max according  to the irf value
                                    if 'z60' in name_irf:
                                        sim_e_min = 0.110
                                        sim_e_max = 5.6234

                                    elif 'z40' in name_irf:
                                        sim_e_min = 0.04
                                        sim_e_max = 5.6234
                                    else:
                                        sim_e_min = 0.03
                                        sim_e_max = 10.

                                    print(f'IRF: {name_irf}, Min_energy: {sim_e_min}, Max_energy: {sim_e_max} ')
                                    print (f'Irf {name_irf} start time: {sim_t_min} , irf {name_irf} stop time: {sim_t_max}')
                                    #print (t_obs_stop)
                                    # -----------------------------------------------------------------------Simulation
                                    sim = ctools.ctobssim()
                                    sim['inmodel'] = inmodel
                                    sim['caldb'] = caldb
                                    sim['irf'] = name_irf
                                    sim['ra'] = 0.
                                    sim['dec'] = 0.
                                    sim['rad'] = sim_rad
                                    sim['tmin'] = sim_t_min
                                    sim['tmax'] = sel_t_max
                                    sim['emin'] = sim_e_min
                                    sim['emax'] = sim_e_max
                                    sim['seed'] = seeds
                                    sim['outevents'] = 'events_full_GRB.fits'
                                    sim.execute()

                                    print(f'\nSimulation {event} - site {site} - {night}: DONE')
                                    # print (f'\nOffset: {offset}, radius: {radius}')
                                    onoff_sim = cscripts.csphagen()
                                    onoff_sim['inobs'] = 'events_full_GRB.fits'
                                    onoff_sim['inmodel'] = fitmodel
                                    onoff_sim['srcname'] = 'Crab'
                                    onoff_sim['ebinalg'] = 'LOG'
                                    onoff_sim['emin'] = sim_e_min
                                    onoff_sim['emax'] = sim_e_max
                                    onoff_sim['enumbins'] = 20
                                    onoff_sim['coordsys'] = 'CEL'
                                    onoff_sim['ra'] = 0.
                                    onoff_sim['dec'] = offset
                                    onoff_sim['rad'] = radius
                                    onoff_sim['caldb'] = caldb
                                    onoff_sim['irf'] = name_irf
                                    onoff_sim['bkgmethod'] = 'REFLECTED'
                                    onoff_sim['use_model_bkg'] = False
                                    onoff_sim['srcregfile'] = 'regioni_on.reg'
                                    onoff_sim['bkgregfile'] = 'regioni_off.reg'
                                    onoff_sim['outobs'] = 'GRBobs.xml'
                                    onoff_sim['outmodel'] = 'GRBmodel.xml'
                                    onoff_sim['stack'] = False
                                    onoff_sim.execute()

                                    a = 0
                                    with open('onoff_off.reg', 'r') as regioni_off:
                                        for line in regioni_off:
                                            if line.startswith('fk5'):
                                                a += 1
                                    regioni_off.close()

                                    on = fits.open('onoff_pha_on.fits')
                                    tbdata_on = on[1].data
                                    conteggi_on = tbdata_on.field('counts')

                                    off = fits.open('onoff_pha_off.fits')
                                    tbdata_off = off[1].data
                                    conteggi_off = tbdata_off.field('counts')

                                    somma_on[i] = 0.
                                    somma_off[i] = 0.

                                    for valore_on in conteggi_on:
                                        somma_on[i] += valore_on

                                    for valore_off in conteggi_off:
                                        somma_off[i] += valore_off

                                    on.close()
                                    off.close()

                                    a = 1. / a

                                on_counts = np.sum(somma_on)
                                off_counts = np.sum(somma_off)
                                print(f'\n\tNumber of counts (per {night}) in the on region: {on_counts}')
                                print(f'\tNumber of counts (per {night}) per off region: {off_counts}')
                                if on_counts == 0:

                                    break
                                else:
                                    try:
                                        valore = 2 * (on_counts * np.log(((1 + a) / a) * (
                                                    on_counts / (off_counts + on_counts))) + off_counts * np.log(
                                            (1 + a) * (off_counts / (off_counts + on_counts))))

                                    except ValueError:
                                        continue

                                    if valore < 0:
                                        valore = 0
                                    sigma = np.sqrt(valore)
                                    print(f'\t\nTime from trigger: {t0}, significance: {sigma}')

                                    details = (f"{event}, {site}, {night},{t0},{t0-t_obs_start}, {sigma}")
                                    append_new_line('sigma_time_evolution.txt', details)
                                    if sigma >= 3:
                                        data[event][site][night]['detection_time'] = t0 - t_obs_start
                                        data[event][site][night]['sigma'] = sigma
                                        print(f'\nDetection time: {t0 - t_obs_start} sec, sigma: {sigma}')
                                        break

                                    else:
                                        del on_counts
                                        del off_counts
                                        t0 += step
                                        #continue

np.save(cfg['path']['sigmaoutput'], data)
