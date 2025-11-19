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

[Hünniger et al., *PLoS Computational Biology* (2014)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003479)

In contrast to the original implementation, which was written entirely in C++, this version introduces a **Python wrapper** that allows the C++ simulation code to be executed directly from Python. This greatly simplifies usage, integration, and extension within modern Python-based workflows.

Furthermore, this updated version is fully compatible with the **[pyABC](https://pyabc.readthedocs.io/en/latest/)** library, enabling efficient and flexible **approximate Bayesian computation (ABC)** for parameter inference, model calibration, and uncertainty quantification.

---
## Requirements



---
## Getting started

### Build process
To build the framework, clone the repository to a local folder (here: `SBM_EVdecoration/`) and use the following commands:

``` ~/SBM_EVdecoration$ mkdir build; cd build ```

Release or Debug mode:

``` ~/SBM_EVdecoration/build$ cmake -DCMAKE_BUILD_TYPE=Release .. ``` or ``` ~/SBM_EVdecoration/build$ cmake -DCMAKE_BUILD_TYPE=Debug .. ```

``` ~/SBM_EVdecoration/build$ make ```

The compiled files can be found in the `build/ folder.

### Using the Python wrapper

To use the Python wrapper, navigate to the `notebooks/` directory and create a symbolic link to the compiled `.so` library generated during the build process. This allows the C++ backend to be accessed directly from Python.

Move into the `notebooks/` folder:

```/SBM_EVdecoration$ cd notebooks```

Create a symbolic link to the compiled shared library:

```/SBM_EVdecoration/notebooks$ ln -s ../build/src/python_sbm.cpython-311-x86_64-linux-gnu.so```

Once the link is created, you can import and use the wrapped module in your Python notebooks.

```import python_sbm as ps```

You can now call the simulation functions directly through Python while relying on the underlying C++ implementation.

---
## License

The project code is licensed under BSD 2-Clause. <br>
See the LICENSE file provided with the code for the full license.

---
## Previous publications



