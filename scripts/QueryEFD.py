import numpy as np
import pandas as pd
from astropy.time import Time,TimeDelta

from chronografAuxtelEnv import query_vent_speed, query_sonic_temp
from chronografAuxtelEnv import query_m1_temp, query_m2_air_temp

from lsst.summit.utils.efdUtils import getEfdData, makeEfdClient
thermal = {'lsst.sal.ATDome.position': ['azimuthPosition', 'mainDoorOpeningPercentage'],
           'lsst.sal.ATMCS.trajectory': ['elevation0', 'azimuth0']
            }

"""
Auxtel environment of the day

The function queryEFDThermalInfo(day, domeFrac=95) will query the EFD for the auxtel environment of the day.
It will return the following information:
- The temperature of M1, M2, and the air
- The sonic temperature and its standard deviation
- The time when the dome was opened
- The the fraction of time the ventilation system is active during the day
- A flag indicating if the ventilation system is active 2 hours prior to dome opening
- The maximum temperature gradient of M1, Air in a window of 6 hours prior to operations (night_t0)
- The temperature M1 at night_t0
- The max temperature M1 at during the day
- The time of the maximum temperature gradient of M1
- The time of the maximum temperature gradient of M1, Air in a window of 6 hours prior to operations (night_t0)
"""

def queryEFDThermalInfo(day, domeFrac=95):
    # Get the EFD data for the specified day
    # the data is resampled to 30 seconds
    # the contains temperatures M1, M2, Air, Truss, Sonic, std(Sonic)
    # the contains dome and ventilation information
    efd = get_efd_day(day)

    # Determine the start time of the night operation based on dome opening fraction
    t0 = get_night_start_time(efd, fracTh=domeFrac)
    efd['night_t0'] = t0

    # Compute the fraction of time the ventilation system is active during the day
    efd['vent_time_frac'] = compute_vent_time_fraction(efd, tend=t0.strftime('%H:%M'))
    
    # Check if the ventilation system is active 2 hours prior to dome opening
    efd['vent_sunset'] = compute_vent_sunset(efd, t0, timeDelta=2)

    # Divide the night time into 4 periods of 2.5 hours each
    efd['night_period'] = getPeriodsNight(efd)

    if 'M1' in efd.columns:
        # Calculate the temperature difference between M1 and the air
        efd['deltaTempM1Air'] = efd.M1 - efd.Air
        efd['deltaTempM1M2'] = efd.M1 - efd.M2

        # Get the maximum temperature difference of M1-Air within a 6-hour window before night operations start (night_t0)
        if not pd.isnull(t0):
            deltaTempMax, deltaTimeMax = get_max_value_window(efd, 'deltaTempM1Air')
            efd['deltaTempM1Air_night_t0_value'] = interpolate_time_series(efd, 'deltaTempM1Air', t0)
            efd['deltaTempM1Air_day_value'] = deltaTempMax
            efd['deltaTempM1Air_day_time'] = deltaTimeMax

            # Get the maximum temperature of M1 within the same 6-hour window
            M1TempMax, M1deltaTimeMax = get_max_value_window(efd, 'deltaTempM1M2')
            efd['deltaTempM1M2_night_t0_value'] = interpolate_time_series(efd, 'deltaTempM1M2', t0)
            efd['deltaTempM1M2_day_value'] = M1TempMax
            efd['deltaTempM1M2_day_time'] = M1deltaTimeMax

    # Save the EFD data to a CSV file
    file_path = f'./data/efd/efd_thermal_data_{day}.csv'
    print(f'File saved: {file_path}')
    efd.to_csv(file_path, index=True)

def get_efd_day(day):
    client = makeEfdClient()
    vals = []
    for path, col in zip(thermal.keys(),thermal.values()):
        df = getEfdData(client, path, columns=col, dayObs=day, warn=False)
        vals.append(df)
    
    vals.append(query_vent_speed(day,client))
    vals.append(query_sonic_temp(day,client))
    vals.append(query_m1_temp(day,client))
    vals.append(query_m2_air_temp(day,client))

    df2 = pd.concat(vals)
    df2 = df2.sort_index().interpolate().resample('30s').mean()    
    df2.index = df2.index.tz_convert('America/Santiago')
    return df2

def get_night_start_time(df, fracTh=95):
    """Get the time when the dome was opened
    in a time window of 6pm to 10pm local time
    """
    # Step 1: Filter the dataframe to the time window between 18:00 to 00:00
    time_window_df = df.between_time('18:00', '00:00')
    
    # Step 2: Find the first occurrence where mainDoorOpeningPercentage >= 95%
    first_open_time = time_window_df.index[time_window_df['mainDoorOpeningPercentage'] >= fracTh].min()
    return first_open_time

def compute_vent_time_fraction(df, t0="12:00", tend="23:00", th=0.6):
    ts = df.maxSpeed.rolling(window='5min').mean().between_time(t0, tend)
    above_threshold = ts >= th
    timeFrac = np.count_nonzero(above_threshold)/len(ts)
    return timeFrac

def compute_vent_sunset(df, tdome, timeDelta=2, th=0.6):
    tstart = tdome - pd.Timedelta(hours=timeDelta)
    t0 = tstart.strftime("%H:%M")
    tend = tdome.strftime("%H:%M")
    tfrac = compute_vent_time_fraction(df, t0, tend, th)
    return tfrac>0.95

def getPeriodsNight(df, deltaTime=2.5):
    """ Divide the night in 4 periods of 2.5 hours

    Day is 0:, 
    First Night Period is 1, 
    and so on...
    """
    first_open_time = df.night_t0.iloc[0]
    # Step 1: Create a new column to indicate the time period
    df['night_period'] = 0
    
    # Step 2: Assign 0 for daytime (before first_open_time)
    df.loc[df.index < first_open_time, 'night_period'] = 0
    
    # Step 3: Assign values 1, 2, 3, and 4 for each of the four 3-hour periods during the night
    night_period_start = first_open_time
    for i in range(1, 5):
        night_period_end = night_period_start + pd.Timedelta(hours=deltaTime)
        df.loc[(df.index >= night_period_start) & (df.index < night_period_end), 'night_period'] = i
        night_period_start = night_period_end

    return df['night_period'].to_numpy()

def get_max_value_window(df, col, deltaTime=6., is_norm=True):
    first_open_time = df.night_t0.iloc[0]
    
    # Define the time window
    start_time_window = first_open_time - pd.Timedelta(hours=deltaTime)
    end_time_window = first_open_time
    
    # Filter the dataframe for the time window
    time_window_df = df.loc[start_time_window:end_time_window]
    
    xvar = time_window_df[col]

    # get the largest temp gradient with norm
    if is_norm:
        xvar = np.abs(xvar)

    if not time_window_df.empty:
        # Find the maximum value of deltaTemp_M1M2
        max_delta_temp_time = xvar.idxmax()
        if pd.isnull(max_delta_temp_time):
            # Use the maximum index value as the time instead
            max_delta_temp_time = xvar.index.max()
        
        # get sign right
        max_delta_temp_six_hours_before = interpolate_time_series(time_window_df, col, max_delta_temp_time) 
        # hours before the begnining of the night
        time_difference_hours = (first_open_time - max_delta_temp_time).total_seconds() / 3600

    else:
        max_delta_temp_six_hours_before = None
        max_delta_temp_time = None  
        time_difference_hours = None

    return max_delta_temp_six_hours_before, time_difference_hours

def get_time_series_value(df, t0, col, deltaTimeSec=30):
    # Define the time window
    end_time_window = t0 - pd.Timedelta(sec=deltaTimeSec)
    start_time_window = t0
    # Filter the dataframe for the time window
    time_window_df = df.loc[start_time_window:end_time_window].mean()
    return time_window_df


def get_efd_exposure(times, delta_t = 30):
    client = makeEfdClient()
    dt = TimeDelta(delta_t, format='sec')
    
    tables = []
    for now in times:
        vals = []
        for path, col in zip(thermal.keys(),thermal.values()):
            df = getEfdData(client, path, columns=col, prePadding=5, begin=now, timespan=dt, warn=False)
            vals.append(df.mean())
        df2 = pd.DataFrame(pd.concat(vals,axis=0)).T
        tables.append(df2)
    return pd.concat(tables, ignore_index=True)

def interpolate_time_series(df, col, given_time):
    # Interpolate deltaTemp_M1M2 for any given time
    # given_time = pd.Timestamp('2024-04-16 20:30:00')  # Example given time
    # interpolated_value_at_given_time = df['deltaTemp_M1M2'].interpolate(method='time').loc[given_time]
    # print("Interpolated value of deltaTemp_M1M2 at given time:", interpolated_value_at_given_time)
    return df[col].interpolate(method='time').loc[given_time]


# fullPath = ['lsst.sal.ESS.airTurbulence.sonicTemperature',
#            'lsst.sal.ATDome.position.azimuthPosition',
#            'lsst.sal.ESS.airTurbulence.sonicTemperatureStdDev',
#            'lsst.sal.ATDome.position.mainDoorOpeningPercentage',
#            'lsst.sal.ESS.airTurbulence.sonicTemperatureStdDev',
#            'lsst.sal.ATMCS.trajectory.elevation0',
#            'lsst.sal.ATMCS.trajectory.azimuth0',
#            'lsst.sal.ESS.airTurbulence.speed0',
#            'lsst.sal.ESS.airTurbulence.speed1',
#            'lsst.sal.ESS.airTurbulence.speed2',
#            'lsst.sal.ESS.airFlow.speed',
#            'lsst.sal.ESS.temperature.temperatureItem0',
#            'lsst.sal.ESS.temperature.temperatureItem1',
#            'lsst.sal.ESS.temperature.temperatureItem3',
#            'lsst.sal.ESS.temperature.temperatureItem3',
#            'lsst.sal.ESS.temperature.temperatureItem4',
#           ]

# paths = []
# cols = []
# for col in fullPath:
#     paths.append(('.').join(col.split('.')[:-1]))
#     cols.append(col.split('.')[-1])

# import numpy as np
# columns = np.array(cols)
# upaths = np.unique(paths) 
# out = dict().fromkeys(upaths)
# for path in upaths:
#     out[path] = columns[np.array(paths)==path]