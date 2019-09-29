from osgeo import gdal, gdalconst
import numpy as np

NP2GDAL_CONVERSION = {      # source: https://borealperspectives.org/2014/01/16/data-type-mapping-when-using-pythongdal-to-write-numpy-arrays-to-geotiff/
  "uint8": 1,
  "int8": 1,
  "uint16": 2,
  "int16": 3,
  "uint32": 4,
  "int32": 5,
  "float32": 6,
  "float64": 7,
  "complex64": 10,
  "complex128": 11,
}

def read(filename):
    # gdal.Driver.Register()

    img = gdal.Open(filename, gdalconst.GA_ReadOnly)

    if img is None:
        # error occurred
        print("ERROR")
        return None
    else:
        cols = img.RasterXSize
        rows = img.RasterYSize

        band = img.GetRasterBand(1)
        np_array = band.ReadAsArray(0, 0, cols, rows)

        return np_array


def write(filename, np_array, complex_output=False):
    if (np_array.dtype.name == "complex64" or np_array.dtype.name == "complex128") and not complex_output:
        np_array_tmp = np.abs(np_array).astype("float64")**2.
    else:
        np_array_tmp = np_array
    try:
        driver = gdal.GetDriverByName('ENVI')
        ds = driver.Create(filename, np_array_tmp.shape[1], np_array_tmp.shape[0], 1, NP2GDAL_CONVERSION[np_array_tmp.dtype.name])
        ds.GetRasterBand(1).WriteArray(np_array_tmp)
        ds.FlushCache()             # write to disk.
    except:
        print("ERROR: problem writing image to disk.")
        return False
    return True
