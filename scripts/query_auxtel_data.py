import os
import numpy as np
import pandas as pd
from astropy.time import Time

from FilterColumns import filterCols, convert_to_mag, add_flux_snr_flag
from computeFWHM import compute_fwhm
from QueryEFD import get_efd_exposure, queryEFDThermalInfo

def header():
    print(5*'-----')
    print('Querying Auxtel Fields')
    print(f'repo: {repo}')
    print(f'collection: {collection} \n')
    print(f'Number of Exposures: {len(expIds)}')
    print(5*'-----'+'\n')

def add_mdata_columns(new, dataId):
    mData = butler.get('raw.metadata', dataId).toDict()
    new['EL'] = mData['ELSTART']
    new['EXPID'] = dataId['visit']
    new['DATE'] = str(mData['DATE'])
    new['FILTER'] = mData['FILTER']
    new['EXPTIME'] = mData['EXPTIME']
    return new

def query_thermal_quantities(df):
    times = Time(np.array(df.DATE).astype(np.datetime64), scale='utc')
    return get_efd_exposure(times)

def query_src_data(expId):
    print(5*'-----')
    print(f'exposure: {expId}')
    day_obs = int(str(expId)[:8])
    dataId = {"instrument": instrument, "exposure.day_obs": day_obs, "visit": expId, "detector": 0}
    
    # Sources table
    exp = butler.get('srcMatchFull', dataId).asAstropy()
    data = convert_to_mag(exp)
    data = add_flux_snr_flag(data,snrTh,gaiaMagLow, gaiaMagHig)
    data = add_mdata_columns(data, dataId)

    # Filter Columns
    df = filterCols(data).to_pandas()

    # Get FWHM
    df['fwhm'] = compute_fwhm(df)

    # Get EFD Data
    # efd_df = query_thermal_quantities(df)
    # new = df# pd.concat([df,efd_df],axis=1)
    new = match_efd_exposure(df, day_obs)
    
    new.to_csv(f'data/gaia_matched_{expId}.csv', index=True)
    print(f'file saved: data/gaia_matched_{expId}.csv')

def match_efd_exposure(df,day):
    efd_fname = './data/efd/efd_thermal_data_%s.csv'%(day)
    efd = pd.read_csv(efd_fname, index_col=0, parse_dates=True)
    efd.index = efd.index.tz_convert(None)

    ## lookup for the date of the expID
    lookup = pd.Timestamp(df.DATE.loc[0])
    output = efd.iloc[efd.index.searchsorted(lookup) - 1]
    columns = list(efd.columns)

    for col in columns:
        df[col] = output[col]
    return df

def join_files(expIds):
    fnames = ['data/gaia_matched_%s.csv'%eid for eid in expIds]
    tables = []
    for i,fname in enumerate(fnames):
        df = pd.read_csv(fname, index_col=0, parse_dates=True)
        tables.append(df)

    df = pd.concat(tables)
    df.to_csv(outfile)
    print(f'Joined file saved: {outfile}')
    pass
    
nCores = 40
outfile = 'gaia_matched_april2024_PREOPS-4985.csv'
instrument = 'LATISS'
repo = '/repo/embargo'
collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240311/w_2024_10/PREOPS-4985'

## Threshold Values
snrTh = 50.
gaiaMagLow = 14. 
gaiaMagHig = 17.

# outfile = 'gaia_matched_march2024_PREOPS-4871.fits'
# instrument = 'LATISS'
# repo = '/repo/embargo'
# collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240201/w_2024_05/PREOPS-4871'

import lsst.daf.butler as dafButler #Gen3 butler
butler = dafButler.Butler(repo, collections=[collection])
registry=butler.registry

# querying all exposure ids in this collection
dataset_refs = list(registry.queryDatasets('srcMatchFull', collections=[collection]))
dataset_refs = sorted(dataset_refs, key=lambda x : x.dataId['visit'])
expIds = [int(ref.dataId['visit']) for ref in dataset_refs]
all_days = [int(str(expId)[:8]) for expId in expIds]
days = np.unique(all_days)

# print header
header()
print('Query EFD')
from joblib import Parallel, delayed
Parallel(n_jobs=nCores)(delayed(queryEFDThermalInfo)(day) for day in days)

print('Query srcMatchFull')
Parallel(n_jobs=nCores)(delayed(query_src_data)(expId) for expId in expIds)

# joining files
join_files(expIds)