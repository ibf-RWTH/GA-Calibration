<h1> GA-Calibration Tool </h1>

<!--## Overview-->
<!--![logo](docs/GUI.PNG)-->

 [**Installation**](#Installation)
| [**Related Projects**](#Related-Projects)


**Genetic Algorithm Calibration Tool** is a Python-based software developed to streamline the calibration of material constitutive models using genetic algorithms. This tool supports calibration with experimental data such as **uniaxial tensile test** and **cyclic loading test** results, offering a powerful solution for researchers and engineers in material science.
### Key Features:
- **Experimental Data Calibration**: Calibrate constitutive models using experimental results from **uniaxial tensile** and **cyclic loading** tests.
- **Supported Material Models**: Includes support for **Crystal Plasticity Models** and the **Chaboche Model** for material behavior simulation.
- **Numerical Integration**: Works with **ABAQUS** and **DAMASK**, enabling users to integrate the tool with these widely-used numerical frameworks.
- **Submodel Calibration in ABAQUS**: Offers calibration within ABAQUS **submodels**, allowing users to focus on localized regions of their models for precise analysis.
- **Single-phase and Dual-phase Materials**: Capable of calibrating both **single-phase** and **dual-phase** materials, broadening its application scope.
- **Concurrent Job Submission**: Supports **concurrent job submission** to efficiently utilize **High-Performance Computing (HPC) clusters**, significantly accelerating the calibration process.
- **Genetic Algorithm Optimization**: Uses genetic algorithms to perform parameter optimization, ensuring a robust search of the parameter space for improved accuracy and performance.

This tool is tailored for users looking to match their experimental results with material models in order to predict material behavior more accurately, providing an efficient and flexible solution for **crystal plasticity** and other **constitutive modeling**.
Do not change the directory names because the tool relies on the structure!
<p align="left"><img src="docs/GA_logo.png" height="400" alt=""> </img></p>

_Note: For developing it is highly recommended to use Python versions Python 3.8._<br>
**If further questions appear please check the lower section or get in touch with us.**


## Installation
The software must be executed on a high-performance computing (HPC) cluster, utilizing [Slurm](https://slurm.schedmd.com/quickstart.html) commands to optimize resource management and improve execution efficiency.
As the first step, conda needs to be installed.
To be sure conda is installed correctly on your system [look up here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)<br>

Git must be installed on the system. Check with:
```
$ git --version
```
If it has not been installed use this:
```
$ conda install -c anaconda git
```
Open the user path and create the directory where the DRAGen repo will be cloned.
Should be like this:
```
(base) C:\Users> cd \Users\{username}
(base) C:\Users\username> mkdir GitRepos
(base) C:\Users\username> cd GitRepos
```
To clone this repository into the desired destination, use:<br>
```
$ git clone https://github.com/ibf-RWTH/GA-Calibration.git
```
To be able to use the tool, the working directory must be set to the location where the repo was downloaded to in the previous step file which is downloaded at the previous step.
Use the commands to go to the exact file by following the path.
```
$ cd GA-Calibration
```
To see the folders on the current point:
```
$ dir
```
Create a virtual environment as follows:<br>
```
$ conda create --name calibration python=3.8
$ conda activate calibration
```
(if an error occurs check your conda installation)<br>
To see the list of the environments on conda:
```
$ conda info --envs
```
Be sure the DRAGen environment is activated it should look somewhat like this:<br>
```
(calibration)....$
```
Install required module packages:

To install requirements
```
(calibration)....$ pip install -r requirements.txt
```

Copy, paste, and then rename the **exampleConfigs.ini** to **configs.ini**. Modify the necessary parameters in **configs. in**. Finally, start Calibration by:<br>
```
(calibration)....$ python main.py
```

## Related Projects

### DRAGen
<p align="center"><img src="docs/DRAGen_logo.png" height="200" alt="DRAGen logo"> </img></p>

[DRAGen](https://github.com/ibf-RWTH/DRAGen) an enhanced version of the Discrete Representative Volume Element (RVE) Automation and Generation Framework, known as DRAGen. Originally devised as an approach for generating Representative Volume Elements based on a Random Sequential Addition (RSA)-Algorithm and discrete tessellation, DRAGen has undergone significant improvements to broaden its capabilities.



### DAMASK
<p align="center"><img src="docs/DAMASK_banner.png" height="100" alt="DAMASK banner"> </img></p>

[DAMASK](https://damask.mpie.de/index.html) (Düsseldorf Advanced Materials Simulation Kit) excels in its ability to handle a variety of simulation programs under different conditions, particularly for advanced high-strength materials. Its capability to address the interconnected nature of deformation, phase transformations, heating effects, and potential damage makes DAMASK an invaluable choice for researchers and practitioners seeking a comprehensive understanding of materials behavior in diverse scenarios.










