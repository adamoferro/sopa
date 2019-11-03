import numpy as np

def read_sofa_sharad_orbit_input(filename_base, pri=1428e-6):
    FAST_TIME_SAMPLING_FREQUENCY = 80e6/3.
    dt = 1./FAST_TIME_SAMPLING_FREQUENCY
    DCG_DELAY = 11.98e-6      # digital chirp generator delay, used in the calculation of the rxwin opening time

    # extra delay for the calculation of the rxwin opening time, which depends on the PRI range
    rxwin_PRI_delay = pri
    if pri > 1500e-6:             # PRF < 670.24 Hz
        rxwin_PRI_delay = 0

    orbit_data = dict()
    input_rxwin_data_dt_pre = np.fromfile(filename_base + "_rxwin", dtype=np.float64)     # saved as number of samples spaced 37.5ns
    orbit_data["X"] = np.fromfile(filename_base + "_x", dtype=np.float64, sep=" ")
    orbit_data["Y"] = np.fromfile(filename_base + "_y", dtype=np.float64, sep=" ")
    orbit_data["Z"] = np.fromfile(filename_base + "_z", dtype=np.float64, sep=" ")
    orbit_data["Lat"] = np.fromfile(filename_base + "_lat", dtype=np.float64, sep=" ")
    orbit_data["Lon"] = np.fromfile(filename_base + "_long", dtype=np.float64, sep=" ")
    orbit_data["Lon"][orbit_data["Lon"] < 0] = 360 + orbit_data["Lon"][orbit_data["Lon"] < 0]
    orbit_data["Radius"] = np.sqrt(orbit_data["X"]**2 + orbit_data["Y"]**2 + orbit_data["Z"]**2)

    orbit_data["rxwin_time_s"] = input_rxwin_data_dt_pre*dt + rxwin_PRI_delay - DCG_DELAY

    return orbit_data
