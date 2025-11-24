# SBM for host-pathogen dynamics in whole blood
State-based model (SBM) framework for simulating host-pathogen dynamics in human whole blood.

This repository contains the code used for data analysis in the publication: <br>
Patitz JJ, Nieuwenhuizen NE, Solomatina A, et al. "Host-to-Pathogen Transfer of Neutrophil Components via Extracellular Vesicles Shields *Candida albicans* from Immune Attack in Human Blood" (submitted)

**Research Group:** Applied Systems Biology  
**Contributions:** Anastasia Solomatina, Paul Rudolph, Sandra Timme  
**Head:** Prof. Dr. Marc Thilo Figge  
**Department:** Applied Systems Biology  
**Institute:** Leibniz Institute for Natural Product Research and Infection Biology (Hans Knöll Institute, HKI)  

**Address:**  
Adolf-Reichwein-Straße 23  
07745 Jena, Germany  

**Website:** [Applied Systems Biology](https://asb.hki-jena.de/)

---
## Description

This framework enables the modeling of host–pathogen interactions in human whole blood. It is based on a previously published computational model of immune response dynamics, described in the following publication:

[Hünniger *et al*., *PLoS Computational Biology* (2014)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003479)

In contrast to the original implementation, which was written entirely in C++, this version introduces a **Python wrapper** that allows the C++ simulation code to be executed directly from Python. This greatly simplifies usage, integration, and extension within modern Python-based workflows.

Furthermore, this updated version is fully compatible with the [pyABC](https://pyabc.readthedocs.io/en/latest/) library, enabling efficient and flexible **approximate Bayesian computation (ABC)** for parameter inference, model calibration, and uncertainty quantification.

---
## Getting started

### Requirements

This project uses a Conda environment to manage all dependencies. </br>
The required packages are listed in ```requirements.txt```.

1. Clone the repository to a local folder (here: `SBM_EVdecoration/`): ``` clone [repository-url]```
2. Set up a Conda environment: ```conda create --name [env-name] --file requirements.txt```
3. Activate the environment: ```conda activate [env-name]```

### Build process
To build the framework, use the following commands:

1. ``` ~/SBM_EVdecoration$ mkdir build; cd build ```
2. Release or Debug mode: ``` ~/SBM_EVdecoration/build$ cmake -DCMAKE_BUILD_TYPE=Release .. ``` or ``` ~/SBM_EVdecoration/build$ cmake -DCMAKE_BUILD_TYPE=Debug .. ```
3. ``` ~/SBM_EVdecoration/build$ make ```

The compiled files can be found in the `build/` folder.

### Using the Python wrapper

To use the Python wrapper, navigate to the `notebooks/` directory and create a symbolic link to the compiled `.so` library generated during the build process:

```/SBM_EVdecoration/notebooks$ ln -s ../build/src/python_sbm.cpython-311-x86_64-linux-gnu.so```

Once the link is created, you can import and use the wrapped module in your Python notebooks.

```import python_sbm as ps```

You can now call the simulation functions directly through Python while relying on the underlying C++ implementation.

---

## Jupyter Notebooks

The `notebooks/` directory contains two Jupyter notebooks:

- `pyabc_sbm.ipynb`: performs ABC parameter inference using the pyABC library.
- `vis_pyabc_sbm.ipynb`: processes the ABC output and generates the figures and visualizations used in the analysis.

---
## Previous publications

1. [Hünniger *et al*., *PLoS Computational Biology* (2014)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003479)
2. [Lehnert *et al*., *Frontiers in Microbiology* (2015)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.00608/full)
3. [Lehnert *et al*., *Scientific Reports* (2021)](www.nature.com/articles/s41598-021-91362-5#Ack1)
4. [Lehnert *et al*., *PLoS ONE* (2021)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0249372)

---
## License

The project code is licensed under BSD 2-Clause. <br>
See the LICENSE file provided with the code for the full license.



