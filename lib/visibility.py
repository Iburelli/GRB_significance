# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import logging
import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d
from astropy.coordinates import AltAz, EarthLocation, get_sun, get_moon, SkyCoord
from astropy.coordinates import solar_system_ephemeris

def complete_irf_name(irfs, site, exposure, azimuth=None):
    """
    Given a list of IRFs zenith angles, the chosen site, exposure and NBS, it returns the completed name of the IRFs
    :param irfs: <z20|z40|z60> str (list of str)
    :param site: <North|South> (str)
    :param exposure: <100s|0.5h|5h|50h> (str)
    :param azimuth: <average> (str or None)
    :return: list of completed IRF names
    """
    new_irfs = []
    for irf in irfs:
        if azimuth is None:
            new_irfs.append(site + '_z' + str(irf) + '_' + exposure)
        else:
            new_irfs.append(site + '_z' + str(irf) + '_' + str(azimuth) + '_' + exposure)
    return new_irfs


class Visibility:
    """
    This class contains methods that compute the source visibility and extracts information on which IRF to use per time interval. This class heavily relies on astropy.
    """

    def __init__(self):
        pass

    def set_jpl_ephemeris(self):
        solar_system_ephemeris.set('jpl')
        return

    def visibility_points(self, trigger, duration, num_points=10, unit='jd'):
        """
        It creates a time grid.
        :param trigger: start time or trigger time (astropy Time object)
        :param duration: duration of the event or the visibility window (astropy Time object)
        :param num_points: number of time grid points (int)
        :param unit: <jd> (str)
        :return:
        """
        self.trigger = trigger
        if unit.lower() != 'jd':
            raise Warning('Time formats other than JD are not implemented yet')
        self.vis_points = np.linspace(0, duration.value, num_points) * u.d
        self.vis_points += trigger
        return self

    def return_to_trigger_frame(self, unit='s'):
        """
        Converts the time grid back to the trigger (or start time) frame.
        :param unit: <s> (str)
        :return:
        """
        if unit.lower() != 's':
            raise Warning('Time formats other than seconds are not implemented yet')
        self.vis_points -= self.trigger
        self.vis_points = self.vis_points * u.s
        return self

    def visibility_altaz(self, source_radec, site, hardcoded=True):
        """
        To each time grid point associates AltAz coordinates.
        :param source_radec: source coordinates (astropy SkyCoord object)
        :param site: (str)
        :return:
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hardcoded:
            site_coords = EarthLocation.of_site(site)
        else:
            if site.lower() in ('north', 'roque de los muchachos'):
                #site_coords = EarthLocation.from_geodetic('342.1184', '28.7606', 2326. * u.meter)
                site_coords = EarthLocation.from_geocentric(5327285.09211954, -1718777.11250295, 3051786.7327476, unit="m")

            elif site.lower() in ('south', 'paranal'):
                #site_coords = EarthLocation.from_geodetic('289.5972', '-24.6253', 2635. * u.meter)
                site_coords = EarthLocation.from_geocentric(1946635.7979987, -5467633.94561753, -2642498.5212285, unit="m")
            else:
                raise Warning(f"{site} is not a valid site choice")
        self.altaz = source_radec.transform_to(AltAz(obstime=self.vis_points, location=site_coords))
        return self


    def sun_position(self):
        """
        Finds the position of the Sun for the observation site, given a time grid.
        :return:
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        sun = get_sun(self.vis_points)
        self.sun_altaz = sun.transform_to(self.altaz)
        return self

    def moon_position(self):
        """
        [WIP] Finds the position of the Moon at the observation site, given a time grid.
        :return:
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        moon = get_moon(self.vis_points)
        self.moon_altaz = moon.transform_to(self.altaz)
        return self

    def get_nighttime(self, twilight=-18, digits=8):
        """Given a twilight altitute threshold for the Sun, it returns twilight and dawn time for each night covering the event duration.
        :param twilight: <0|-6|-12|-18> civil, naval or astronomical twilight or night thresholds (int). Default -18.
        :return: dictionary containing 'twilight' and 'dawn' time of the Sun for each nighttime window of the event
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        self.sun_position()
        windows = {'start': [], 'stop': []}
        current = None
        for idx, t in enumerate(self.vis_points.value):
            previous = current
            # night condition
            if self.sun_altaz[idx].alt.value < twilight:
                current = True
            else:
                current = False

            if idx == 0 and current is True:
                windows['start'].append(self.vis_points[idx].value)
                continue
            elif idx == len(self.vis_points.value) - 1 and current is True:
                windows['stop'].append(self.vis_points[idx].value)
                break
            elif previous != current:
                x = [self.sun_altaz[idx - 1].alt.value, self.sun_altaz[idx].alt.value]
                y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                f = interp1d(np.array(x), np.array(y))
                if previous is False and current is True:
                    windows['start'].append(f(twilight))
                elif previous is True and current is False:
                    windows['stop'].append(f(twilight))

        if len(windows['start']) != 0:
            windows['stop'] = np.concatenate(np.around(windows['stop'], digits), axis=None)
            windows['start'] = np.concatenate(np.around(windows['start'], digits), axis=None)
        else:
            windows['start'] = np.array([-9.0])
            windows['stop'] = np.array([-9.0])
        return windows

    def get_nighttime_moon_veto(self, twilight=-18, moon_alt_max=-0.5, digits=8):
        """Given a twilight altitute threshold for the Sun, it returns twilight and dawn time for each night covering the event duration with a veto on the Moon presence.
        :param twilight: <0|-6|-12|-18|integer> civil, naval or astronomical twilight or night thresholds (int). Default -18 deg (integer).
        :return: dictionary containing 'start' and 'stop' time of the Sun for each nighttime window of the event
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        self.sun_position()
        self.moon_position()
        windows = {'start': [], 'stop': []}
        current = None
        for idx, t in enumerate(self.vis_points.value):
            previous = current
            # visibility conditions
            sun_cond = np.array(self.sun_altaz.alt.value < twilight)
            moon_cond = np.array(self.moon_altaz.alt.value < moon_alt_max)
            if sun_cond[idx] and moon_cond[idx]:
                current = True
            else:
                current = False
            if idx == 0 and current is True:
                windows['start'].append(self.vis_points[idx].value)
                continue
            elif idx == len(self.vis_points.value) - 1 and current is True:
                windows['stop'].append(self.vis_points[idx].value)
                break
            elif previous != current and previous != None:
                use = None
                # assign x and y arrays
                if sun_cond[idx] != sun_cond[idx - 1]:
                    x = [self.sun_altaz[idx - 1].alt.value, self.sun_altaz[idx].alt.value]
                    y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                    use = 'sun'
                elif moon_cond[idx] != moon_cond[idx - 1]:
                    x = [self.moon_altaz[idx - 1].alt.deg, self.moon_altaz[idx].alt.deg]
                    y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                    use = 'moon'
                # interpolate
                f = interp1d(np.array(x), np.array(y))
                if previous is False and current is True:
                    if use == 'sun':
                        windows['start'].append(f(twilight))
                    elif use == 'moon':
                        windows['start'].append(f(moon_alt_max))
                elif previous is True and current is False:
                    if use == 'sun':
                        windows['stop'].append(f(twilight))
                    elif use == 'moon':
                        windows['stop'].append(f(moon_alt_max))


        if len(windows['start']) != 0:
            windows['stop'] = np.concatenate(np.around(windows['stop'], digits), axis=None)
            windows['start'] = np.concatenate(np.around(windows['start'], digits), axis=None)
        else:
            windows['start'] = np.array([-9.0])
            windows['stop'] = np.array([-9.0])
        return windows


    def get_nighttime_moonlight(self, twilight=-18, moon_sep=30, fov_rad=0, moonpha=0, max_moonpha=0.8, digits=8):
        """Given a twilight altitute threshold for the Sun and a minimum separation from the Moon position, it returns twilight and dawn time for each night covering the event duration.
        :param twilight: <0|-6|-12|-18|integer> civil, naval or astronomical twilight or night thresholds (int). Default -18 deg (integer).
        :param moon_sep: <integer> minimum angular separation between moon and FoV border. Default 30 deg.
        :param fov_rad: <integer> radius of FoV in degrees. Default 2.5.
        :return: dictionary containing 'start' and 'stop' time of the Sun for each nighttime window of the event
        """
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        self.sun_position()
        self.moon_position()
        dist = self.altaz.separation(self.moon_altaz)
        # add FoV radius to minimum separation from Moon
        moon_sep += fov_rad
        windows = {'start': [], 'stop': []}
        current = None
        for idx, t in enumerate(self.vis_points.value):
            previous = current
            # visibility conditions
            sun_cond = np.array(self.sun_altaz.alt.value < twilight)
            moon_cond = np.array(dist.deg > moon_sep)
            if sun_cond[idx] and moon_cond[idx]:
                current = True
            else:
                current = False
            if idx == 0 and current is True:
                windows['start'].append(self.vis_points[idx].value)
                continue
            elif idx == len(self.vis_points.value) - 1 and current is True:
                windows['stop'].append(self.vis_points[idx].value)
                break
            elif previous != current and previous != None:
                use = None
                # assign x and y arrays
                if sun_cond[idx] != sun_cond[idx - 1]:
                    x = [self.sun_altaz[idx - 1].alt.value, self.sun_altaz[idx].alt.value]
                    y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                    use = 'sun'
                elif moon_cond[idx] != moon_cond[idx - 1]:
                    x = [dist[idx - 1].deg, dist[idx].deg]
                    y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                    use = 'moon'
                # interpolate
                f = interp1d(np.array(x), np.array(y))
                if previous is False and current is True:
                    if use == 'sun':
                        windows['start'].append(f(twilight))
                    elif use == 'moon':
                        windows['start'].append(f(moon_sep))
                elif previous is True and current is False:
                    if use == 'sun':
                        windows['stop'].append(f(twilight))
                    elif use == 'moon':
                        windows['stop'].append(f(moon_sep))


        if len(windows['start']) != 0:
            windows['stop'] = np.concatenate(np.around(windows['stop'], digits), axis=None)
            windows['start'] = np.concatenate(np.around(windows['start'], digits), axis=None)
        else:
            windows['start'] = np.array([-9.0])
            windows['stop'] = np.array([-9.0])
        return windows

    def associate_irf_zenith_angle(self, thresholds=(24,36,57), zenith_angles=(60, 40, 20), digits=8):
        """
        Given altitude lower bounds and the respective IRFs zenith angles, it returns which IRF zenith
        angle to use when.
        :param thresholds: (tuple or list of int)
        :param zenith_angles: (tuple or list of int)
        :return: dictionary containing 'start' time and 'stop' time associated to each 'zref' zenith angle
        """
        if len(thresholds) != len(zenith_angles):
            raise AttributeError('thresholds length must be equal to zenith_angles length')
        if not hasattr(self, 'vis_points'):
            raise AttributeError('Must invoke visibility_points() before using this method')
        if not hasattr(self, 'altaz'):
            raise AttributeError('Must invoke visibility_altaz() before using this method')
        thresholds = sorted(thresholds)
        zenith_angles = sorted(zenith_angles, reverse=True)
        windows = {'start': [], 'stop': [], 'zref': []}
        previous, current = None, None
        for idx, alt in enumerate(self.altaz.alt.value):
            for j, z in enumerate(thresholds):
                if alt < z:
                    if j == 0:
                        current = None
                    else:
                        current = zenith_angles[j - 1]
                    break
                elif j == len(thresholds) - 1:
                    current = zenith_angles[j]

            if idx == 0 and previous != current and current is not None:
                windows['start'].append(self.vis_points[idx].value)
                windows['zref'].append(current)
            elif idx == len(self.altaz) - 1 and current is not None:
                windows['stop'].append(self.vis_points[idx].value)
            elif previous != current and previous is None:
                x = [self.altaz[idx - 1].alt.value, self.altaz[idx].alt.value]
                y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                f = interp1d(x, y)
                windows['start'].append(f(thresholds[0]))
                windows['zref'].append(current)
            elif previous != current and current is None:
                x = [self.altaz[idx - 1].alt.value, self.altaz[idx].alt.value]
                y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                f = interp1d(x, y)
                windows['stop'].append(f(thresholds[0]))
            elif previous != current and (current and previous) is not None:
                x = [self.altaz[idx - 1].alt.value, self.altaz[idx].alt.value]
                y = [self.vis_points[idx - 1].value, self.vis_points[idx].value]
                f = interp1d(x, y)
                if previous > current:
                    windows['stop'].append(f(thresholds[zenith_angles.index(current)]))
                    windows['start'].append(f(thresholds[zenith_angles.index(current)]))
                else:
                    windows['stop'].append(f(thresholds[zenith_angles.index(previous)]))
                    windows['start'].append(f(thresholds[zenith_angles.index(previous)]))
                windows['zref'].append(current)
            previous = current
        if len(windows['zref']) != 0:
            windows['start'] = np.concatenate(np.around(windows['start'], digits), axis=None)
            windows['stop'] = np.concatenate(np.around(windows['stop'], digits), axis=None)
            windows['zref'] = np.concatenate(np.around(windows['zref'], digits), axis=None)
        else:
            windows['start'] = np.array([-9.0])
            windows['stop'] = np.array([-9.0])
            windows['zref'] = np.array([-9.0])
        return windows

    def associate_irf_one_night(self, source_radec, start_time, duration, site, num_points,
                                thresholds=(24,36, 57), zenith_angles=(20, 40, 60)):
        """
        Shortcut method, it returns the IRF zenith angle to use when, given the source position, site, start time and duration of the event, number of time grid points and altitude thresholds with associate zenith angles.
        :param source_radec: source coordinates (astropy SkyCoord object)
        :param start_time: starting time of the night or trigger (astropy Time object)
        :param duration: duration of the event or night (astropy Time object)
        :param site: <North|South> (str)
        :param num_points: time grid points (int)
        :param thresholds: altitude lower bounds for the IRFs (list of int)
        :param zenith_angles: associated zenith angles (list of int)
        :return irfs: dictionary containing 'start' time and 'stop' time associated to each 'zref' zenith angle
        """
        self.visibility_points(start_time, duration, num_points)
        self.visibility_altaz(source_radec, site)
        irfs = self.associate_irf_zenith_angle(thresholds, zenith_angles)
        return irfs
