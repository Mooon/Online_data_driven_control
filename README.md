# Online data driven control (Automatica 2022)

This is an official implementation of the paper [Online learning of data-driven controllers for unknown switched linear systems](https://www.sciencedirect.com/science/article/pii/S000510982200379X).

The scripts are available in both MATLAB and Python.

## Overview

The main goal of the code is to demonstrate the online learning method applied to a control system with switching dynamics. The control input is optimized in real time using convex optimization, ensuring system stability after each switch between different system models.

Both versions of the code simulate the behavior of an aircraft engine under different operating modes, switching between two different linear state-space models at specific fault times. Online learning is used to calculate the optimal control gain at each time step.

## Installation

### MATLAB
To run the MATLAB version, ensure that you have a recent version of MATLAB installed with the necessary toolboxes (Control System Toolbox and Optimization Toolbox). Simply clone this repository and run the script in MATLAB. The controller gain K is computed by solving a data-based convex program using [CVX](https://cvxr.com/cvx/).

### Python

```
python -m venv venv
source venv/bin/activate
pip install -r py_requirements.txt
python aircraft_online_learning.py
```

## Contributing

Feel free to contribute to this repository. If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## Citation
If you find this repo useful in your research, please consider citing our paper as follows:

```
@article{rotulo2022online,
  title={Online learning of data-driven controllers for unknown switched linear systems},
  author={Rotulo, Monica and De Persis, Claudio and Tesi, Pietro},
  journal={Automatica},
  volume={145},
  pages={110519},
  year={2022},
  publisher={Elsevier}
}
```
