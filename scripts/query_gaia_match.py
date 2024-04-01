from gaia_match import query_wcs_auxtel, query_gaia_data, ploting_field, get_lims
from gaia_match import  sky_match_tables, save_table
import os

def header():
    print(5*'-----')
    print('Querying Gaia Fields')
    print(f'repo: {repo}')
    print(f'collection: {collection} \n')
    print(f'Number of Exposures: {len(expIds)}')
    print(5*'-----'+'\n')

outfile = 'gaia_matched_march2024_PREOPS-4871.fits'
instrument = 'LATISS'
repo = '/repo/embargo'
collection = 'LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20240201/w_2024_05/PREOPS-4871'

import lsst.daf.butler as dafButler #Gen3 butler
butler = dafButler.Butler(repo, collections=[collection])
registry=butler.registry

# querying all exposure ids in this collection
dataset_refs = list(registry.queryDatasets('calexp', collections=[collection]))
dataset_refs = sorted(dataset_refs, key=lambda x : x.dataId['visit'])
expIds = [int(ref.dataId['visit']) for ref in dataset_refs]

# print header
header()

# for test purposes only
# expIds = expIds[:10]

for expId in expIds:
    if not os.path.isfile(f'data/gaia_matched_{expId}.fits'):
        print(5*'-----')
        print(f'exposure: {expId}')
        day_obs = int(str(expId)[:8])
        dataId = {"instrument": instrument, "exposure.day_obs": day_obs, "visit": expId, "detector": 0}

        exp = butler.get('calexp', dataId)
        scr = butler.get('src', visit=expId, detector=0, collections=collection).asAstropy()

        # mask = scr['base_CircularApertureFlux_3_0_instFlux']/scr['base_CircularApertureFlux_3_0_instFluxErr'] > 50
        # scr = scr[mask]

        mData = butler.get('raw.metadata', dataId).toDict()

        print('Finding Gaia Sources')
        print('Field %i'%expId)

        print('1 - WCS Solutions')
        # wheader, wcs_auxtel = query_wcs_auxtel(scr, center_ra = mData['RA'], center_dec = mData['DEC'])
        wcsButler = butler.get('calexp.wcs', dataId)
        x0, y0    = wcsButler.getPixelOrigin().getX(), wcsButler.getPixelOrigin().getY()

        print('2 - Gaia Finding')
        tab = query_gaia_data(wcsButler, border = 10, mag_limit=18, row_limit=-1, is_butler=True)

        #print('3 - Ploting results')
        #ploting_field(tab, scr, exp, expId, x0=x0, y0=y0, path='./figures')

        print('4 - Sky Match Sources And Save Results')
        new = sky_match_tables(tab, scr, wcsButler, radius=1.2, is_butler=True)
        new['EXPID'] = expId
        new['DATE'] = str(mData['DATE'])
        new['FILTER'] = mData['FILTER']
        new['EXPTIME'] = mData['EXPTIME']
        new['X0'] = x0
        new['Y0'] = y0
        save_table(f'data/gaia_matched_{expId}.fits', new)
        #print('\n')

## Join all files
import glob
fnames = glob.glob('data/gaia_matched_*')
for i,fname in enumerate(fnames):
    print(f'{i}: {fname}')
    new_expIds.append(get_expid(fname))
    
table = vstack(tables)
table.write(outfile,format='fits',overwrite=True)