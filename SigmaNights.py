import sys
import os
from os.path import join, isfile, isdir
import gammalib
import ctools
import cscripts
import gammapy
from astropy.io import fits
import numpy as np
import math as m
import astropy
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.utils import iers
from astropy.io import fits
import warnings
import argparse
import yaml
import logging
import numpy as np
import glob


parser = argparse.ArgumentParser(description='TThe significance of a GRB observation is computed at different IRFs according to a visibility tabele, created with runCatVisibility.py. A configuation YAML file is required, the output is saved as NPY binary file.')
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
            raise ValueError(f'Specified template {runid} does not exist in catalog')
    runids = [cfg['path']['xmlfilename']]
else:
    runids = cfg['path']['xmlfilename']
    for runid in runids:
        if not isfile(join(catalog, runid)):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
runids = sorted(runids)

#--------------------------------------------------------------reading visibility table

data= np.load(output, allow_pickle=True, encoding='latin1', fix_imports=True).flat[0]
events = list(data.keys())
sites = list(data[events[0]].keys())
#data = table.copy()
####--------------function to append new line to txt file
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




#-----------------------------------------------importing ctools variables

#generating random seeds
seeds = np.random.randint(1,1000,size=cfg['ctools']['iterations'])
#simulation variables
caldb=cfg['ctools']['caldb']
sim_rad = cfg['ctools']['rad']
#sim_e_max=cfg['ctools']['Emax']
fitmodel=cfg['ctools']['fitmodel']
#duration of simulation ( da sistemare, potrei voler simulare tutto l'intervallo)
duration = cfg['ctools']['duration']

#inmodel selection
for runid in runids:

    if type(cfg['path']['xmlfilename']) == str:
        inmodel = join(catalog, runid)

    elif cfg['path']['xmlfilename'] == None:
        inmodel = runid


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


            for site in sites:
                if site == 'North':
                    print(f'\nProcessing site {site}')
                    for night in data[event][site]:
                        print(f'\nProcessing {night}')
                        if type(data[event][site]) == float:
                            print(f'\tThis contains NaNs ---> the source is not observable due to daylight or moon.')
                        if data[event][site][night]['irfs']['zref'][0] == -9.0:
                            print(f'\tThis contains NaNs event---> the source is not observable at the site.')
                            data[event][site][night]['sigma_ON/OFF'] = -9.0

                        else:
                            t_obs_start= data[event][site][night]['irfs']['start'][0]
                            t_obs_start=(t_obs_start - trigger)*86400

                            on_counts=np.zeros(shape=len(seeds))
                            off_counts=np.zeros(shape=len(seeds))
                            valore=np.zeros(shape=len(seeds))
                            radice=np.zeros(shape=len(seeds))
                            for j, seed in enumerate(seeds):
                                seed = int(seed)
                                print(f'\nseed number {j+1}: {seed}')


                                somma_on=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))
                                somma_off=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))

                                for i in range(len(data[event][site][night]['irfs']['zref'])):
                                    t_min = data[event][site][night]['irfs']["start"][i]
                                    t_max= data[event][site][night]['irfs']["stop"][i]

                                    name_irf = (f'{site}_z{data[event][site][night]["irfs"]["zref"][i]}_0.5h')

                                    #converting times from jd to seconds from trigger
                                    sim_t_min=(t_min - trigger)*86400
                                    sim_t_stop=(t_max - trigger)*86400

                                    if sim_t_min > t_obs_start+duration:
                                        somma_on[i]=0.
                                        somma_off[i]=0.
                                        break
                                    elif sim_t_stop - t_obs_start < duration:
                                        sim_t_max=sim_t_stop
                                    else:
                                        sim_t_max= t_obs_start + duration


                                    
                                    #selection of e_min according  to the irf value
                                    if 'z60' in name_irf:
                                        sim_e_min=0.03
                                        #name_irf=f'{site}_z40_5h'
                                    else:
                                        sim_e_min=0.03
                                    print (f'Irf : {name_irf}')
                                    if 'z20' in name_irf:
                                        sim_e_max=10.
                                    else:
                                        sim_e_max=5.6234
                    #-----------------------------------------------------------------------Simulation
                                    sim = ctools.ctobssim()
                                    sim['inmodel'] = inmodel
                                    sim['caldb'] = caldb
                                    sim['irf'] = name_irf
                                    sim['ra'] = 0.
                                    sim['dec'] = 0.
                                    sim['rad'] = sim_rad
                                    sim['tmin'] = sim_t_min
                                    sim['tmax'] = sim_t_max
                                    sim['emin'] = sim_e_min
                                    sim['emax'] = sim_e_max
                                    sim['seed'] = seed
                                    sim['outevents'] ='events_full_GRB.fits'
                                    sim.execute()

                                    print(f'\nSimulation site {site} - {night}: DONE')

                                    onoff_sim = cscripts.csphagen()
                                    onoff_sim['inobs'] =  'events_full_GRB.fits'
                                    onoff_sim['inmodel'] = fitmodel
                                    onoff_sim['srcname'] = 'Crab'
                                    onoff_sim['ebinalg'] = 'LOG'
                                    onoff_sim['emin'] = sim_e_min
                                    onoff_sim['emax'] = sim_e_max
                                    onoff_sim['enumbins'] = 20
                                    onoff_sim['coordsys'] = 'CEL'
                                    onoff_sim['ra'] = 0.
                                    onoff_sim['dec'] = cfg['ctools']['offset']
                                    onoff_sim['rad'] = cfg['ctools']['radius']
                                    onoff_sim['caldb'] = caldb
                                    onoff_sim['irf'] =name_irf
                                    onoff_sim['bkgmethod'] = 'REFLECTED'
                                    onoff_sim['use_model_bkg'] = False
                                    onoff_sim['srcregfile'] =  'regioni_on.reg'
                                    onoff_sim['bkgregfile'] =  'regioni_off.reg'
                                    onoff_sim['outobs'] =  'GRBobs.xml'
                                    onoff_sim['outmodel'] =  'GRBmodel.xml'
                                    onoff_sim['stack'] =False
                                    onoff_sim.execute()

                                    a = 0
                                    with open( 'onoff_off.reg', 'r') as regioni_off:
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
                                    somma_on[i]=0.
                                    somma_off[i]=0.

                                    for valore_on in conteggi_on:
                                        somma_on[i] += valore_on

                                    for valore_off in conteggi_off:
                                        somma_off[i] += valore_off

                                        on.close()
                                        off.close()
                                    #somma_off[i] = somma_off[i]/(a*1.0)

                                    a=1/a

                                on_counts[j]=np.sum(somma_on)
                                off_counts[j]=np.sum(somma_off)

                                print (f'\n\tNumber of counts (per {night}) in the on region: {on_counts[j]}')
                                print (f'\tNumber of counts (per {night}) per off region: {off_counts[j]}')



                                try:
                                    valore[j] = 2*(on_counts[j] * np.log(((1+a)/a)*(on_counts[j]/(off_counts[j]+on_counts[j]))) + off_counts[j] * np.log((1+a)*(off_counts[j]/(off_counts[j]+on_counts[j]))))
                                    if valore[j] < 0:
                                        valore[j] = 0
                                    radice[j]=np.sqrt(valore[j])
                                except ValueError:

                                    continue

                                print (f'\n\tSignificance site {site},{night}, seed {seed} :{radice[j]}')


                            media_on=np.mean(on_counts)
                            media_off=np.mean(off_counts)
                            var_on=np.std(on_counts)
                            var_off=np.std(off_counts)
                            sigma = np.mean(radice)
                            var = np.std(radice)
                            data[event][site][night]['sigma_ON/OFF'] = sigma
                            data[event][site][night]['sigma_var'] = var
                            print (f'\n\tMean significance, site {site}, {night}: {sigma}, variance: {var}')
                            #print (f'\tMean counts ON: {media_on}, mean counts OFF: {media_off}')

                            sigma_mean_counts=np.sqrt ( 2*(media_on * m.log(((1+a)/a)*(media_on/(media_off+media_on))) + media_off * m.log((1+a)*(media_off/(media_off+media_on)))))

                            details = (f"{event}, {site}, {night}, {duration}, {media_on}, {var_on}, {media_off},{var_off}, {sigma}, {var}")

                            append_new_line('details.txt', details)


#print (f"{nn} events out of {n+1} can't be detected by CTA North")
#print (f"{ss} events out of {n+1} can't be detected by CTA South")
#print (data)
np.save(cfg['path']['sigmaoutput'] , data)
