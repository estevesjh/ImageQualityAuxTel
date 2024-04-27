import os
import numpy as np
import pandas as pd
from astropy.time import Time

from FilterColumns import filterCols, convert_to_mag, add_flux_snr_flag
from computeFWHM import compute_fwhm
from QueryEFD import get_efd_exposure, queryEFDThermalInfo

dropCols = ['ref_id', 'ref_coord_ra', 'ref_coord_dec', 'ref_phot_g_mean_flux','ref_phot_g_mean_fluxErr',
            'ref_centroid_x','ref_centroid_y','src_coord_ra', 'src_coord_dec', 'src_base_SdssCentroid_x', 'src_base_SdssCentroid_y',
            'src_base_SdssCentroid_xErr','src_base_SdssCentroid_yErr','src_base_SdssShape_instFlux','src_base_SdssShape_instFluxErr',
            'src_base_PsfFlux_instFlux','src_base_PsfFlux_instFluxErr','src_base_CircularApertureFlux_3_0_instFlux','src_base_CircularApertureFlux_3_0_instFluxErr',
            'src_base_SdssShape_xx','src_base_SdssShape_yy','src_base_SdssShape_xy','src_ext_shapeHSM_HsmSourceMoments_xx','src_ext_shapeHSM_HsmSourceMoments_yy','src_ext_shapeHSM_HsmSourceMoments_xy',
            'DATE','FILTER','night_t0']

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
    
    df['DATE'] = pd.to_datetime(df['DATE'],utc=True)
    df['night_t0'] = pd.to_datetime(pd.Series(df['night_t0']),utc=True)
    df['hours'] = (df.DATE - df['night_t0']).dt.total_seconds() / 3600

    return df

## retrieve Date columns
def retrieve_date_column(df, de):
    de['DATE'] = 0
    de['DATE'] = np.nan
    masks = [df['EXPID'] == eid for eid in de.index]

    # Iterate over each row in de and assign the corresponding datetime value
    for eid, mask in zip(de.index, masks):
        date_values = df.DATE.loc[mask]
        filter_values = df.FILTER.loc[mask]
        de.loc[eid, 'DATE'] = date_values.to_numpy()[0]
        de.loc[eid, 'FILTER'] = filter_values.to_numpy()[0]
    
    # Convert 'DATE' column to datetime dtype
    de['DATE'] = pd.to_datetime(de['DATE'])
    de['EXPID'] = de.index
    de.columns = [col.replace('_mean', '') for col in de.columns]
    return de.set_index('DATE')

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

def create_seeing_file():
    """To Do
    """
    df = pd.read_csv(outfile, index_col=0)
    df['DATE'] = pd.to_datetime(df['DATE'],utc=True, format='ISO8601')
    #df = convert_to_hour_only(df)

    flagG = (df.ref_mag_g>14)&(df.ref_mag_g<17)
    flagSNR = (df.snr>50)
    flagT = flagG&flagSNR

    df = df.loc[flagT]
    
    dstar = df.drop(columns=dropCols)
    dg = dstar.groupby('EXPID').agg(['mean', 'std'])
    dg.columns = ['{}_{}'.format(col, stat) for col, stat in dg.columns]

    # filter exposures
    dg = dg.loc[dg.fwhm_std < 0.5]
    dgn = retrieve_date_column(df, dg)

    outfile2 = outfile.replace('.csv','_seeing.csv')
    dgn.to_csv(outfile2, index=True)
    print(f'Joined file saved: {outfile2}')
    pass

nCores = 40
instrument = 'LATISS'
repo = '/repo/embargo'

# outfile = 'gaia_matched_march2024_PREOPS-4871.csv'
# outfile = 'gaia_matched_april2024_PREOPS-4985.csv'
outfile = 'gaia_matched_april2024_PREOPS-5069.csv'

# collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240201/w_2024_05/PREOPS-4871'
# collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240311/w_2024_10/PREOPS-4985'
collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240414/w_2024_15/PREOPS-5069'

## Threshold Values
snrTh = 50.
gaiaMagLow = 14. 
gaiaMagHig = 17.

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
# from joblib import Parallel, delayed
# print('Query EFD')
# Parallel(n_jobs=nCores)(delayed(queryEFDThermalInfo)(day) for day in days)

# print('Query srcMatchFull')
# Parallel(n_jobs=nCores)(delayed(query_src_data)(expId) for expId in expIds)

# # joining files
join_files(expIds)

create_seeing_file()