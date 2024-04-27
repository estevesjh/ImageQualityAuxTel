import numpy as np
import pandas as pd
from lsst.summit.utils.efdUtils import getEfdData, makeEfdClient

"""
Series of Python functions to query the EFD for AuxTel Environment Data
"""

def query_m2_air_temp(day, client, path='lsst.sal.ESS.temperature'):
    """
    SELECT "temperatureItem1" as "Air", "temperatureItem2" as "Truss", "temperatureItem3" as "M2" 
           FROM "efd"."autogen"."lsst.sal.ESS.temperature" WHERE salIndex = 201 AND "sensorName" = 'AuxTel-ESS01' 
           AND time > :dashboardTime: AND time < :upperDashboardTime:

    There are three auxtel sensor names:
    'AuxTel-ESS01', 'AuxTel-ESS02', 'AuxTel-ESS03'

    Example:
        from lsst.summit.utils.efdUtils import makeEfdClient
        client = makeEfdClient()
        df = query_m2_air_temp('20210403', client)

    """
    cindex = ['salIndex','sensorName']
    columns = ['temperatureItem%i'%i for i in range(6)]
    df = getEfdData(client, path, columns=cindex+columns, dayObs=day, warn=False)

    # m2, air, truss
    mask = (df.salIndex==201) & (df.sensorName == 'AuxTel-ESS01')
    tempM2Air = df[columns].loc[mask].copy()
    column_mapping = {'temperatureItem1': 'Air', 'temperatureItem2': 'Truss','temperatureItem3': 'M2'}
    tempM2AirTruss = tempM2Air.rename(columns=column_mapping)
    return tempM2AirTruss.resample('30s').mean()


def query_m1_temp(day, client, path='lsst.sal.ESS.temperature'):
    """
    From chronograph: https://summit-lsp.lsst.codes/chronograf/sources/1/dashboards/91?
    
    SELECT "temperatureItem0" as "M1 sensor 1", "temperatureItem1" as "M1 sensor 2", "temperatureItem2" as "M1 sensor 3", 
           "temperatureItem3" as "M1 sensor 4", "temperatureItem4" as "M1 sensor 5", "temperatureItem5" as "M1 sensor 6" 
           FROM "efd"."autogen"."lsst.sal.ESS.temperature" 
           WHERE salIndex = 201 AND "sensorName" = 'AuxTel-ESS03' 
           AND time > :dashboardTime: AND time < :upperDashboardTime:

    There are three auxtel sensor names:
    'AuxTel-ESS01', 'AuxTel-ESS02', 'AuxTel-ESS03'
    
    Example:
        from lsst.summit.utils.efdUtils import makeEfdClient
        client = makeEfdClient()
        df = query_m1_temp('20210403', client)

    """
    cindex = ['salIndex','sensorName']
    columns = ['temperatureItem%i'%i for i in range(6)]
    df = getEfdData(client, path, columns=cindex+columns, dayObs=day, warn=False)
    
    # get M1 temps
    mask = (df.salIndex==201) & (df.sensorName == 'AuxTel-ESS03')
    tempM1 = df[columns].loc[mask].copy()
    
    newcolumns = ['M1_sensor%i'%i for i in range(6)]
    column_mapping = {old_name: new_name for old_name, new_name in zip(columns, newcolumns)}
    tempM1 = tempM1.rename(columns=column_mapping)
    tempM1['M1'] = tempM1['M1_sensor4']
    return tempM1.resample('30s').mean()

def query_vent_speed(day, client, path="lsst.sal.ESS.airFlow"):
    """
    from chronograph: https://summit-lsp.lsst.codes/chronograf/sources/1/dashboards/168?
    
    SELECT mean("speed") AS "mean_speed" FROM "efd"."autogen"."lsst.sal.ESS.airFlow" 
           WHERE salIndex=204 AND time > :dashboardTime: AND time < :upperDashboardTime: 
           GROUP BY time(:interval:) FILL(null)
    interval = 720sec
    """
    cindex = ['salIndex','sensorName']
    columns = ['speed', 'maxSpeed', 'direction']
    df = getEfdData(client, path, columns=cindex+columns, dayObs=day, warn=False)
    mask = (df.salIndex==204)
    vent = df[columns].loc[mask].copy()
    return vent.resample('720s').mean()


def query_sonic_temp(day, client, sensorName='ESS04', path="lsst.sal.ESS.airTurbulence"):
    """
    From chronograph: https://summit-lsp.lsst.codes/chronograf/sources/1/dashboards/91?
    
    SELECT mean("sonicTemperature") AS "mean_sonicTemperature" FROM "efd"."autogen"."lsst.sal.ESS.airTurbulence" 
           WHERE time > :dashboardTime: AND time < :upperDashboardTime: AND salIndex = 201 
           GROUP BY time(:interval:) FILL(null)
    
    interval = 240sec

    There are two sensor names: 'AuxTel-ESS04' and 'AuxTel-GillLabJack01'
    """
    cindex = ['salIndex','sensorName']
    columns = ['sonicTemperature', 'sonicTemperatureStdDev']
    df = getEfdData(client, path, columns=cindex+columns, dayObs=day, warn=False)

    ## Days before Nov-2023 do not have sensorName
    if 'sensorName' in df.columns:
        mask = (df.sensorName == 'AuxTel-%s'%sensorName)
        sonic = df[columns].loc[mask].copy()
        sonic = sonic.resample('240s').mean()
    else:
        sonic = df
    return sonic