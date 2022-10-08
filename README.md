# PRIDE-Doppler-LaRa 

PRIDE-Doppler-LaRa is an object-oriented Python tool for research in the field of PRIDE (Planetary Radio Interferometry and Doppler Experiments).

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/cfortunylombra/PRIDE-Doppler-LaRa/graphs/commit-activity) ![python 3.8](https://img.shields.io/badge/version-latest-blue.svg) ![python 3.8](https://img.shields.io/badge/python-3.8-blue.svg) ![platform Linux](https://img.shields.io/badge/platform-Linux-blue.svg)

---

### Developers
* **C. Fortuny-Lombraña**, Master Student, Astrodynamics and Space Missions, TU Delft
* **D. Dirkx**, Assistant Professor, Astrodynamics and Space Missions, TU Delft

[![Link MailTo](https://img.shields.io/badge/MailTo-developers-blue.svg)](mailto:C.FortunyLombrana@student.tudelft.nl;D.Dirkx@tudelft.nl?subject=PRIDE-Doppler-LaRa:Query)

---

### Pre-requisites
* For Windows Users: Windows Subsystem for Linux ([Installing WSL](https://docs.microsoft.com/en-us/windows/wsl/install))
	- All procedures, including the following prerequisite, assume the use of WSL. Power users who wish to do otherwise, must do so at their own risk, with reduced support from the team.
    - Note that WSL is entirely separated from Windows. This means that even if Anaconda/Miniconda, Python, or any other software is already installed on Windows, it needs to be installed again in WSL, as a Linux install.
* Anaconda installation ([Installing Anaconda](https://docs.anaconda.com/anaconda/install/))
* `tudat-bundle` ([Installing `tudat-bundle`](https://github.com/tudat-team/tudat-bundle))

To install `tudat-bundle`, please follow carefully the [setup instructions](https://github.com/tudat-team/tudat-bundle#setup).

Inside `tudat-bundle`, the following GitHub are utilized:
* `tudat`: https://github.com/cfortunylombra/tudat-fork.git -> branch: `feature/maximum_viability_fix` using `git checkout -b feature/maximum_viability_fix`
* `tudatpy`: https://github.com/cfortunylombra/tudatpy-fork.git -> branch: `feature/expose_observation_estimation_fix` using `git checkout -b feature/expose_observation_estimation_fix`



### Documentation

* [TudatPy](https://docs.tudat.space/en/stable/)
* [API](https://py.api.tudat.space/en/latest/)

---

### How to get started

The `src` folder contains all the source codes. The most important Python files are the following:

* `py_gsanalysis_LaRa.py` and `py_gsanalysis_InSight.py`: Preliminary analysis of the ExoMars-LaRa and InSight-RISE observations
* `read_sebastian_files.py`: Read the Doppler residuals from InSight-RISE (real data)
* `py_readrealobservations.py`: Functions to transform date to Julian days
* `py_allandeviation.py`: Compute Allan deviations for the PRIDE observations
* `py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py`: Environment Setup for the Precise Orbit Determination (RISE and LaRa with DSN-PRIDE) using a fixed cross-correlation of the Doppler observations between ground stations
* `py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_complex.py`: Environment Setup for the Precise Orbit Determination (RISE and LaRa with DSN-PRIDE) using a the metric developed of the cross-correlations at the ground stations

The Python scripts used to plot, to verify and addiotional porpuses can be found in `plot_codes`,`verification_codes` and `additional_codes` folders.

---

### Help

In case a problem or issue related to code is found, please create a new issue in GitHub.

---
