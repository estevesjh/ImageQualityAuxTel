import pandas as pd
import astropy.units as u

columns = ['ref_id','ref_coord_ra','ref_coord_dec','ref_phot_g_mean_flux','ref_phot_g_mean_fluxErr']
columns+= ['ref_centroid_x','ref_centroid_y','src_coord_ra','src_coord_dec']
columns+= ['src_base_SdssCentroid_x','src_base_SdssCentroid_y','src_base_SdssCentroid_xErr','src_base_SdssCentroid_yErr']
columns+= ['src_base_SdssShape_instFlux','src_base_SdssShape_instFluxErr']
columns+= ['src_base_PsfFlux_instFlux','src_base_PsfFlux_instFluxErr']
columns+= ['src_base_CircularApertureFlux_3_0_instFlux','src_base_CircularApertureFlux_3_0_instFluxErr']
columns+= ['src_base_SdssShape_xx','src_base_SdssShape_yy','src_base_SdssShape_xy']
columns+= ['src_ext_shapeHSM_HsmSourceMoments_xx','src_ext_shapeHSM_HsmSourceMoments_yy','src_ext_shapeHSM_HsmSourceMoments_xy']
columns+= ['ref_mag_g','snr','snr_psf','snr_aperture3']
columns+= ['flag_snr','flag_mag']
columns+= ['EL','EXPID', 'DATE', 'FILTER', 'EXPTIME']

def filterCols(df):
    return df[columns]

def convert_to_mag(df):
    """Convert columns of flux to magnitude in the AB system
        also computes the SNR
    """
    df['ref_mag_g'] = df['ref_phot_g_mean_flux'].to(u.ABmag).value    
    # NOT WORKING!!! TO DEBUG
    # df['src_mag'] = (df['src_base_SdssShape_instFlux'] * u.nJy).to(u.AB).value
    # df['src_mag_psf'] = df['src_base_PsfFlux_instFlux'].to(u.AB).value
    # df['src_mag_aperture3'] = df['base_CircularApertureFlux_3_0_instFlux'].to(u.AB).value
    return df

def add_flux_snr_flag(df, snrTh = 50., gaiaMagLow = 14., gaiaMagHig = 17.):
    df['snr'] = df['src_base_SdssShape_instFlux']/df['src_base_SdssShape_instFluxErr']
    df['snr_psf'] = df['src_base_PsfFlux_instFlux']/df['src_base_PsfFlux_instFluxErr']
    df['snr_aperture3'] = df['src_base_CircularApertureFlux_3_0_instFlux']/df['src_base_CircularApertureFlux_3_0_instFluxErr']

    # add flag
    df['flag_snr'] = df['snr']>50.
    df['flag_mag'] = (df['ref_mag_g'] > gaiaMagLow)&(df['ref_mag_g'] <= gaiaMagHig)
    return df