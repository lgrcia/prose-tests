# simulated data
# --------------
import numpy as np
from prose import Telescope

RAW = snakemake.input[0]

import numpy as np
import matplotlib.pyplot as plt

# time = np.linspace(0, 0.15, 100) + 2450000
# target_dflux = 1 + np.sin(time*100)*1e-2

# from prose.tutorials import simulate_observation

# # so we have the same data
# np.random.seed(40)

# fits_folder = "./tutorial_dataset"
# simulate_observation(time, target_dflux, RAW)

_ = Telescope({
    "name": "A",
    "trimming": [40, 40],
    "pixel_scale": 0.66,
    "latlong": [31.2027, 7.8586],
    "keyword_light_images": "light"
})

# reduction
# ---------
from prose import FitsManager, Image, Sequence, blocks, Observation

# ref
fm = FitsManager(RAW, depth=2)
ref = Image(fm.all_images[0])
calibration = Sequence([
    blocks.Calibration(darks=fm.all_darks, bias=fm.all_bias, flats=fm.all_flats),
    blocks.Trim(),
    blocks.SegmentedPeaks(), # stars detection
    blocks.Cutouts(),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])

calibration.run(ref, show_progress=False)

# potometry
photometry = Sequence([
    *calibration[0:-1],                
    blocks.psf.Moffat2D(reference=ref),
    blocks.detection.LimitStars(min=3),
    blocks.Twirl(ref.stars_coords),    
    blocks.Set(stars_coords=ref.stars_coords),
    blocks.AffineTransform(data=False, inverse=True),
    blocks.BalletCentroid(),                           
    blocks.PhotutilsAperturePhotometry(scale=ref.fwhm),
    blocks.Peaks(),
    blocks.XArray(
        ("time", "jd_utc"),
        ("time", "bjd_tdb"),
        ("time", "flip"),
        ("time", "fwhm"),
        ("time", "fwhmx"),
        ("time", "fwhmy"),
        ("time", "dx"),
        ("time", "dy"),
        ("time", "airmass"),
        ("time", "exposure"),
        ("time", "path"),
        ("time", "sky"),
        (("time", "apertures", "star"), "fluxes"),
        (("time", "apertures", "star"), "errors"),
        (("time", "apertures", "star"), "apertures_area"),
        (("time", "apertures", "star"), "apertures_radii"),
        (("time", "apertures"), "apertures_area"),
        (("time", "apertures"), "apertures_radii"),
        ("time", "annulus_rin"),
        ("time", "annulus_rout"),
        ("time", "annulus_area"),
        (("time", "star"), "peaks"),
        name="xarray"
    ),
    blocks.AffineTransform(stars=False, data=True),
    blocks.Stack(ref, name="stack"),
])

photometry.run(fm.all_images)

# diff flux
obs = Observation(photometry.xarray.to_observation(photometry.stack.stack, sequence=photometry))
obs.target = 0
obs.broeg2005()
obs.save(snakemake.output[0])

import yaml
yaml.dump({"processing_time": float(photometry.processing_time)}, open(snakemake.output[1], "w"))