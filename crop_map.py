import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from regions.core import PixCoord
from regions.shapes.circle import CirclePixelRegion
from astropy import units as u
import sys


"""
    Crops a circular region around the image center. Takes image name to be cropped
    and radius in ["] (default is 100"). Output image is: <input_name>_crop.fits.

    Usage:
        > python crop_map.py 4C24_850_mf_crop_cal_snr.fits 200

    Updated: 04.06.2021
    trg
"""

#
#
def image_cutout_save(filename, radius):

    # Load the image and the WCS
    hdu = fits.open(filename)[0]

    data = hdu.data[0,:,:]
    # Modify header keywords
    hdu.header['NAXIS'] = 2
    hdu.header['WCSAXES']=2
    # Delete all keywords of 3rd axis
    del hdu.header['NAXIS3']
    del hdu.header['CRPIX3']
    del hdu.header['CUNIT3']
    del hdu.header['CDELT3']
    del hdu.header['CTYPE3']
    del hdu.header['LBOUND3']
    del hdu.header['CRVAL3']
    del hdu.header['OBSGEO-Z']
    wcs = WCS(hdu.header)

    # Cutout side length = 2*radius
    side = radius*2.
    size = u.Quantity((side+10, side+10), u.arcsec)
    position = pixel_to_skycoord(hdu.header['CRPIX1'], hdu.header['CRPIX2'], wcs)
    # Make the cutout, including the WCS
    cutout = Cutout2D(data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Make mask and multiply with data
    center = PixCoord(hdu.header['CRPIX1'], hdu.header['CRPIX2'])
    reg = CirclePixelRegion(center, radius/2.)
    mask = reg.to_mask(mode='exact')
    shape = (hdu.header['NAXIS1'], hdu.header['NAXIS2'])
    #foo = mask.to_image(shape)
    #foo[foo == 0] = 'nan'
    #hdu.data = hdu.data*foo

    # Write the cutout to a new FITS file
    cutout_filename = filename[:-5]+"_crop.fits"
    hdu.writeto(cutout_filename, overwrite=True)
    
    # Load the variance and the WCS
    hdu1 = fits.open(filename)[1]

    data1 = hdu1.data[0,:,:]
    # Modify header keywords
    hdu1.header['NAXIS'] = 2
    hdu1.header['WCSAXES']=2
    # Delete all keywords of 3rd axis
    del hdu1.header['NAXIS3']
    del hdu1.header['CRPIX3']
    del hdu1.header['CUNIT3']
    del hdu1.header['CDELT3']
    del hdu1.header['CTYPE3']
    del hdu1.header['LBOUND3']
    del hdu1.header['CRVAL3']
    del hdu1.header['OBSGEO-Z']
    wcs1 = WCS(hdu1.header)

    # Cutout side length = 2*radius
    position1 = pixel_to_skycoord(hdu1.header['CRPIX1'], hdu1.header['CRPIX2'], wcs1)
    # Make the cutout, including the WCS
    cutout1 = Cutout2D(data1, position=position1, size=size, wcs=wcs1)

    # Put the cutout image in the FITS HDU
    hdu1.data = cutout1.data

    # Update the FITS header with the cutout WCS
    hdu1.header.update(cutout1.wcs.to_header())

    # Make mask and multiply with data
    center1 = PixCoord(hdu1.header['CRPIX1'], hdu1.header['CRPIX2'])
    reg1 = CirclePixelRegion(center1, radius/2.)
    mask1 = reg1.to_mask(mode='exact')
    shape1 = (hdu1.header['NAXIS1'], hdu1.header['NAXIS2'])
    #foo1 = mask1.to_image(shape1)
    #foo1[foo1 == 0] = 'nan'
    #hdu1.data = hdu1.data*foo1

    fits.append(cutout_filename, hdu1.data, header=hdu1.header)


if __name__ == '__main__':
    filename = sys.argv[1]

    if len(sys.argv) > 2:
        radius = float(sys.argv[2])
    else:
        radius = 100.   # ["] Default value

    image_cutout_save(filename, radius)
