from osgeo import gdal, gdalconst


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


def write(filename, np_array):
    try:
        driver = gdal.GetDriverByName('ENVI')
        ds = driver.Create(filename, np_array.shape[1], np_array.shape[0], 1, gdal.GDT_Float64)
        ds.GetRasterBand(1).WriteArray(np_array)
        ds.FlushCache()             # write to disk.
    except:
        print("ERROR: problem writing image to disk.")
        return False
    return True
