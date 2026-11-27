import os
import warnings
import emcee
import corner
import yaml

# Disabling numpy's default parallelization, which interferes with sampling performance.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Parallelisation for sampling
from multiprocessing import Pool

import sys
import numpy as np
import gc

# Removing warnings
import logging
logging.disable(logging.WARNING)

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath('__file__'))))

# Importing ndspec tools
from ndspec.Response import ResponseMatrix
import ndspec.FitTimeAvgSpectrum as fitspec
import ndspec.XspecInterface as XSModels
from ndspec.Plotting import set_plot_style
import ndspec.SamplingUtils as sampling

# lmfit for fitting
from lmfit import Model as LM_Model

# fits file handling
from astropy.io import fits

# Loading RELXILL and the XSPEC library
relxillpath = os.path.expanduser("/data/cluster4/BradleyAndJacob/relxill/librelxill.so")
relxillfile = os.path.expanduser("/data/cluster4/BradleyAndJacob/relxill/lmodel_relxill.dat")
model_library = XSModels.ModelInterface()
model_library.add_library(relxillpath, relxillfile)

# Observation ID and instrument for data to be looked at
obj = "MCG-6-30-15"
obsid = "0693781301"
instrument = "pn"

# Path to the data
path = f"/home/jq23394/Data/{obj}/{obsid}/PRODUCTS/"

# Ensuring the BACKSCAL keyword is in the header
with fits.open(f"{path}{obsid}_{instrument}_pi.fits", mode="update") as hdul:
    hdul["PRIMARY"].header["BACKSCAL"] = hdul["SPECTRUM"].header["BACKSCAL"]
    hdul.flush()
with fits.open(f"{path}{obsid}_{instrument}_bkg_pi.fits", mode="update") as hdul:
    hdul["PRIMARY"].header["BACKSCAL"] = hdul["SPECTRUM"].header["BACKSCAL"]
    hdul.flush()

# Loading the response matrices
rmfpath = f"{path}{obsid}_{instrument}_rmf.fits"
nicer_matrix = ResponseMatrix(rmfpath)
arfpath = f"{path}{obsid}_{instrument}_arf.fits"
nicer_matrix.load_arf(arfpath)

# Loading the spectrum
spectrum_fit = fitspec.FitTimeAvgSpectrum()
spectrum_fit.set_data(nicer_matrix, f"{path}{obsid}_{instrument}_pi.fits", background=f"{path}{obsid}_{instrument}_bkg_pi.fits")
spectrum_fit.ignore_energies(0,3.0)
spectrum_fit.ignore_energies(10.0,100.0)

# Loading models
model_library.add_models("powerlaw", "relxill", "gaussian", "phabs")

# Defining model wrappers
def ndspec_powerlaw(ear,PhoIndex,norm_powerlaw):
    par_array = np.array([PhoIndex,norm_powerlaw])
    model = model_library.powerlaw(ear,par_array)
    return model

def ndspec_relxill(ear,Index1,Index2,Rbr,a,Incl,Rin,Rout,z,gamma,logxi,Afe,Ecut,refl_frac,norm_relxill):
    norm_relxill = 10**norm_relxill
    par_array = np.array([Index1,Index2,Rbr,a,Incl,Rin,Rout,z,gamma,logxi,Afe,Ecut,refl_frac,norm_relxill])
    model = model_library.relxill(ear,par_array)
    return model

def ndspec_gaussian(ear,LineE,Sigma,norm_gaussian):
    norm_gaussian = 10**norm_gaussian
    par_array = np.array([LineE,Sigma,norm_gaussian])
    model = model_library.gaussian(ear,par_array)
    return model

def ndspec_phabs(ear,nH):
    par_array = np.array([nH])
    model = model_library.phabs(ear,par_array)
    return model

with open("Config.yaml", "r") as file:
    config = yaml.safe_load(file)

unfolded_model = LM_Model(ndspec_phabs)*LM_Model(ndspec_relxill)+LM_Model(ndspec_gaussian)

start_params = unfolded_model.make_params(
    nH = config["nH"],
    Index1 = config["Index1"],
    Index2 = config["Index2"],
    Rbr = config["Rbr"],
    a = config["a"],
    Incl = config["Incl"],
    Rin = config["Rin"],
    Rout = config["Rout"],
    z = config["z"],
    gamma = config["gamma"],
    logxi = config["logxi"],
    Afe = config["Afe"],
    Ecut = config["Ecut"],
    refl_frac = config["refl_frac"],
    norm_relxill = config["norm_relxill"],
    LineE = config["LineE"],
    Sigma = config["Sigma"],
    norm_gaussian = config["norm_gaussian"]
)

# start_params.pretty_print()

spectrum_fit.set_model(unfolded_model, params=start_params)

spectrum_fit.fit_data()

set_plot_style()
fig, ax = spectrum_fit.plot_model(
    units="data", plot_components=False, return_plot=True, residuals="ratio",
    data_kwargs={"color":"black", "markersize":0.1, "linewidth":0.7},
    model_kwargs={"color":"red", "linewidth":1},
    resid_kwargs={"color":"black", "markersize":0.1, "linewidth":0.7}
)

import ndspec.SamplingUtils as sampling

# Defining the data and model used by the sampler
sampling.set_sampling_data(spectrum_fit)
sampling.set_sampling_model(spectrum_fit)

# Initialising the sampling parameters using the best fits just found
theta = sampling.set_sampling_parameters(spectrum_fit.model_params)
ndim = len(theta)

# Printing the parameters
sampling.sampling_params.pretty_print()

# Setting priors as normal distributions with the mean set as the best fit parameters just found
# Standard deviations were selected manually from the scale of the parameter
nH_prior =              sampling.priorNormal(sampling.sampling_params["nH"].value, 1e-1)
Index1_prior =          sampling.priorNormal(sampling.sampling_params["Index1"].value, 1e0)
Index2_prior =          sampling.priorNormal(sampling.sampling_params["Index2"].value, 1e0)
Rbr_prior =             sampling.priorNormal(sampling.sampling_params["Rbr"].value, 1e0)
a_prior =               sampling.priorTruncNormal(sampling.sampling_params["a"].value, 1e-1, 0, 0.998)
Incl_prior =            sampling.priorNormal(sampling.sampling_params["Incl"].value, 1e1)
gamma_prior =           sampling.priorNormal(sampling.sampling_params["gamma"].value, 1e0)
logxi_prior =           sampling.priorNormal(sampling.sampling_params["logxi"].value, 1e-1)
Afe_prior =             sampling.priorNormal(sampling.sampling_params["Afe"].value, 1e0)
# Ecut_prior =            sampling.priorNormal(sampling.sampling_params["Ecut"].value, 1e1)
refl_frac_prior =       sampling.priorNormal(sampling.sampling_params["refl_frac"].value, 1e1)
norm_relxill_prior =    sampling.priorNormal(sampling.sampling_params["norm_relxill"].value, 1e0)
# Sigma_prior =           sampling.priorNormal(sampling.sampling_params["norm_gaussian"].value, 1e0)
norm_gaussian_prior =   sampling.priorNormal(sampling.sampling_params["Sigma"].value, 1e-1)

my_priors = {
    "nH": nH_prior,
    "Index1": Index1_prior,
    "Index2": Index2_prior,
    "Rbr": Rbr_prior,
    "a": a_prior,
    "Incl": Incl_prior,
    "gamma": gamma_prior,
    "logxi": logxi_prior,
    "Afe": Afe_prior,
    "refl_frac": refl_frac_prior,
    "norm_relxill": norm_relxill_prior,
    "norm_gaussian": norm_gaussian_prior,
}
    # "Ecut": Ecut_prior,
    # "Sigma": Sigma_prior

sampling.set_sampling_priors(spectrum_fit, my_priors)

# Calculating the likelihood
likelihood = sampling.mcmc_gaussian_likelihood(theta)
print("mcmc likelihood",likelihood)

import zeus

# Setting the number of walkers and steps
nwalkers = int(4*ndim)
steps = int(1000)

# Initial position of the random walkers, set to be within 1% of the best fit parameters
pos = theta*(1+ 1e-2 * np.random.randn(nwalkers, len(theta)))
# print(pos)

print(sampling.sampling_names)

with Pool(8) as pool:

    slice_sampler = zeus.EnsembleSampler(nwalkers, ndim, log_prob,pool=pool)
    slice_sampler.run_mcmc(pos, steps)


# Getting the correlation time
tau = sampler.get_autocorr_time(tol=1)
print(tau)
thin = np.mean(tau)

samples = sampler.get_chain(discard=int(2*thin), thin=int(0.5*thin), flat=True)

labels=[
    "Index1",
    "Index2",
    "Rbr",
    "a",
    "Incl",
    "gamma",
    "logxi",
    "Afe",
    "Ecut",
    "refl_frac",
    "norm_relxill",
    "norm_gaussian",
    "Sigma"
]

# Plotting a corner plot
fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True)

gc.collect()

print(np.mean(samples, axis=0))