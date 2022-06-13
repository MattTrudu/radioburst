#!/usr/bin/env python
import sys
import os
import numpy as np
import astropy.units as u
import h5py
#import your


class Radioburst:

    def __init__(
    self,
    id,
    telescope,
    fc,
    bw,
    type,
    mjdstart,
    duration,
    dt,
    df,
    tbursts,
    dm,
    widths,
    mask,
    data,
    ):

        self.id = id
        self.telescope = telescope
        self.fc = fc
        self.bw = bw
        self.type = type

        if type not in ["TotalIntensity", "FullStokes"]:
            raise ValueError("Data Type can be only either TotalIntensity or FullStokes...")

        self.mjdstart = mjdstart
        self.duration = duration
        self.dt = dt
        self.df = df
        self.tbursts = tbursts
        self.dm = dm
        self.widths = widths
        self.mask = mask
        self.data = data



        if mask.shape[0] < data.shape[1]:
            raise ValueError("The mask has to be with the same number of spectral channels...")

        if data.shape[0] > 4:
            raise ValueError("Data cannot contain more than four matricies...")

        if type in ["TotalIntensity"]:
            if data.shape[0] > 1:
                raise ValueError("Data in total intensity cannot contain more than a matrix...")

        if len(tbursts) != len(widths):
            raise ValueError("The number of arrival time of bursts must be equal to the number of widths...If widths are unknwon, set them to 0.")



    def header(self):

        string =  f"""


 Observational Specifics \n

 Telescope                        = {self.telescope}
 Central Frequency                = {self.fc.value} {self.fc.unit}
 Bandwidth                        = {self.bw.value} {self.bw.unit}
 \n
 Data Information \n

 Burst ID                         = {self.id}
 Data Type                        = {self.type}
 MJD of the 1st bin (topocentric) = {self.mjdstart}
 Data length                      = {self.duration.value} {self.duration.unit}
 Time resolution                  = {self.dt.value} {self.dt.unit}
 Frequency resolution             = {self.df.value} {self.df.unit}
 Data Shape (pols, chans, bins)   = {self.data.shape}

 Number of bursts in the data     = {len(self.tbursts.value)}
 Bursts time                      = {self.tbursts.value} {self.tbursts.unit}
 Dispersion Measure (DM)          = {self.dm.value} {self.dm.unit}
 Bursts temporal width            = {self.widths.value} {self.widths.unit}

        """

        return string


    def save(self, filename = None, outdir = None):

        if outdir is None:

            outdir =  os.getcwd()

        if filename is None:

            filename = f"{self.id}.h5"

        filename = os.path.join(outdir, filename)


        with h5py.File(filename, "w") as f:

            f.attrs["id"] = self.id
            f.attrs["telescope"] = self.telescope
            f.attrs["fc_v"] = self.fc.value
            f.attrs["fc_u"] = str(self.fc.unit)
            f.attrs["bw_v"] = self.bw.value
            f.attrs["bw_u"] = str(self.bw.unit)
            f.attrs["type"] = self.type
            f.attrs["mjdstart"] = self.mjdstart
            f.attrs["duration_v"] = self.duration.value
            f.attrs["duration_u"] = str(self.duration.unit)
            f.attrs["dt_v"] = self.dt.value
            f.attrs["dt_u"] = str(self.dt.unit)
            f.attrs["df_v"] = self.df.value
            f.attrs["df_u"] = str(self.df.unit)
            f.attrs["tburst_v"] = self.tbursts.value
            f.attrs["tburst_u"] = str(self.tbursts.unit)
            f.attrs["dm_v"] = self.dm.value
            f.attrs["dm_u"] = str(self.dm.unit)
            f.attrs["width_v"] = self.widths.value
            f.attrs["width_u"] = str(self.widths.unit)


            f.create_dataset("mask", data = self.mask, dtype = self.mask.dtype)
            f.create_dataset("data", data = self.data, dtype = self.data.dtype)


            f.close()



def ReadFile(filename):

    with h5py.File(filename, "r") as f:

        id = f.attrs["id"]
        telescope = f.attrs["telescope"]
        fc = f.attrs["fc_v"] * u.Unit(f.attrs["fc_u"])
        bw = f.attrs["bw_v"] * u.Unit(f.attrs["bw_u"])
        type = f.attrs["type"]
        mjdstart = f.attrs["mjdstart"]
        duration = f.attrs["duration_v"]  * u.Unit(f.attrs["duration_u"])
        dt = f.attrs["dt_v"] * u.Unit(f.attrs["dt_u"])
        df = f.attrs["df_v"] * u.Unit(f.attrs["df_u"])
        tbursts = f.attrs["tburst_v"] * u.Unit(f.attrs["tburst_u"])
        dm = f.attrs["dm_v"] * u.Unit(f.attrs["dm_u"])
        widths = f.attrs["width_v"] * u.Unit(f.attrs["width_u"])

        mask = np.array(f["mask"])

        data = np.array(f["data"])

        burst = Radioburst(id, telescope, fc, bw, type, mjdstart, duration, dt, df, tbursts, dm, widths, mask, data)

        return burst


"""

id = "SRT-P-2020/11/09-01"
fc = 336 * u.MHz
bw = 80 * u.MHz
mjdstart = 58941.89021322189
dt = 0.1 * u.ms
df = - 1 * u.MHz
duration = 1 * u.s

tburst = [0.35, 0.5, 0.7] * u.s
DM = 348.772 *u.pc * u.cm**(-3)
width  = [7,10,5] * u.ms

nchan = np.rint(bw.value / abs(df.to(u.MHz).value)).astype(int)
freqs = np.linspace(fc.value + bw.value / 2, fc.value - bw.value / 2, nchan )
times = np.arange(0,duration.value,dt.to(u.s).value)


mask = np.zeros((freqs.shape[0],) , dtype = bool)
mask[5:20] = 1

data = np.random.normal(0,1, size = (1, freqs.shape[0], times.shape[0]) )

burst = Radioburst(id,"SRT",fc,bw,"TotalIntensity", mjdstart, duration, dt, df, tburst, DM, width, mask, data)

print(burst.header())

print(burst.mask)
print(burst.data.shape)

burst.save(filename = "srt01.frb")

filename = "/Users/matt/Desktop/Codes/radioburst/radioburst/srt01.frb"

burst2 = ReadFile(filename)

print(burst2.header())

"""
