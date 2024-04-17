import numpy as np

def compute_fwhm(df, pixel_scale=1/10.):
    fwhm = psf_size_to_fwhm(df, pixel_scale)
    corr_filter = np.array([getFilterSeeingCorrection(fi) for fi in np.array(df.FILTER)])
    corr_airmass = getAirmassCorrection(df.EL.to_numpy())
    return fwhm * corr_filter * corr_airmass

def psf_size_to_fwhm(df, pixel_scale=1/10.):
    psf_size = np.array(df['src_base_SdssShape_xx'] + df['src_base_SdssShape_yy'])
    fwhm = 2.355*np.sqrt(psf_size)*pixel_scale # to arcsec
    return fwhm

def getAirmassCorrection(el):
    """
    Parameters
    ----------
    elevation : `float`
        The elevation in degrees.

    Returns
    -------
    correctionFactor : `float`
        The correction factor to apply to the seeing.

    adapted from: https://github.com/lsst-sitcom/summit_utils/blob/4ca7327cacbc26bba972d40b41100de23f862d4f/python/lsst/summit/utils/utils.py#L925
    """
    airmass = 1/np.cos((90-el)*np.pi/180)
    return airmass ** (-0.6)

def getFilterSeeingCorrection(filterName):
    """Get the correction factor for seeing due to a filter.

    Parameters
    ----------
    filterName : `str`
        The name of the filter, e.g. 'SDSSg_65mm'.

    Returns
    -------
    correctionFactor : `float`
        The correction factor to apply to the seeing.

    Raises
    ------
        ValueError raised for unknown filters.
    """
    match filterName:
        case "SDSSg_65mm":
            return (474.41 / 500.0) ** 0.2
        case "SDSSr_65mm":
            return (628.47 / 500.0) ** 0.2
        case "SDSSi_65mm":
            return (769.51 / 500.0) ** 0.2
        case "SDSSz_65mm":
            return (871.45 / 500.0) ** 0.2
        case "SDSSy_65mm":
            return (986.8 / 500.0) ** 0.2
        case _:
            return np.nan
            #raise ValueError(f"Unknown filter name: {filterName}")
