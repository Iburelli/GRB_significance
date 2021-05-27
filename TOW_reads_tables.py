import numpy as np
import argparse

# parse command line inputs
parser = argparse.ArgumentParser(description='This script is a simple example on how to read a NPY binary file.')
parser.add_argument('-f', '--file', required=True, type=str, help='configuration yaml file')
# configuration file
filename = parser.parse_args().file

data = np.load(filename, allow_pickle=True, encoding='latin1', fix_imports=True).flat[0]

print(f'These are the parent keywords: \n\t{sorted(data.keys())}')
events = list(data.keys())
print(f'Inside each of the parent keyword you will find the following keys: \n\t{sorted(data[events[0]].keys())}')
sites = list(data[events[0]].keys())
print(f'Each key is a dictionary: \n\t{data[events[0]][sites[0]].keys()}')

print(f'In example you can access data with a nested cycle:')
for n, event in enumerate(events):
    for site in sites:
        print(f'Template {event} - Site {site} ')
        if type(data[event][site]) == float:
            print(f'\tThis contains NaNs ---> the source is not observable due to daylight or moon.')
        else:
            print(f"\tSite keys are dictionaries containing details about single nights.")
            for night in data[event][site].keys():
                print(f"\n\t{night}:")
                if type(data[event][site][night]['irf']) == float:
                    print(f'\tThis contains NaNs ---> the source is not observable at the site.')
                else:
                    #print('\n')
                    print(f"\tStart time of single time intervals (here is printed just the first value): {data[event][site][night]['t_start'][0]}")
                    print(f"\tStop time of single time intervals (here is printed just the first value): {data[event][site][night]['t_stop'][0]}")
                    print(f"\tIrf associated to the interval: {data[event][site][night]['irf'][0]}")
                    print(f"\tObservation significance calculated up to the interval stop time: {data[event][site][night]['significance'][0]}")
                    print(f"\tVariance associated to the significance (if null menas the source has been simulated only once): {data[event][site][night]['variance'][0]}")
                    print('\n')
    if n > 2:
        print("Let's stop after 3 templates.")
        break
