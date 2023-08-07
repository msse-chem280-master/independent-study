# Chem 280 - Independent Study Project

For this project, you will write and benchmark the performance of three implementations of a Monte Carlo simulation. 

## Part 1 - Writing Monte Carlo Code

We have provided you with several files to start this project.

1. Code for performing a Monte Carlo simulation of a Lennard Jones fluid in the NVT ensemble using the Python Standard Library. The file `python_files/monte_carlo.py` contains the Monte Carlo functions. The file `analysis.py` contains scripts for performing analysis. Please note that the analysis script utilizes NumPy.
2. Some starting utility functions in the `cpp_files` folder. You have a random number generator and a function for reading data from files.
3. A data folder containing a sample configuration from NIST. When you view the [NIST Benchmarks site](https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-cuboid-cell), this corresponds to configuration 1.

Your initial task in this independent study project is to refactor the provided simulation code in two distinct ways:

1. Use NumPy arrays to represent the coordinates, and modify the necessary functions to be compatible with NumPy operations.
2. Develop a separate implementation of the simulation in C++.

The objective with these modifications is to notably enhance the performance of your code. When you finish, you will have three versions of your code - a Python Standard Library version, a NumPy version, and a C++ version.

If any topic requires further clarification, please refer to the [Chem 280 website](https://msse-chem-280-2023.github.io/index.html). Each day lists different topics that we discussed during the synchronous course sessions. Note that results from your simulation results should be [comparable with Monte Carlo isotherm from NIST](https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm). Note that links to these benchmarks are also available on the course website.

## Part 2 - Performance Evaluation

For this project, systematically evaluate the performance of the three Monte Carlo code versions you have. Consider investigating the code speed by either:

- Increasing the number of steps with a constant number of particles or
- Increasing the number of particles while maintaining a consistent number of steps.

After obtaining your timing data, make a plot or multiple plots showcasing your average timing results. Address the following questions in your assessment:

- Which version is the fastest?
- Do the different versions scale in distinct ways?
- Can you identify potential reasons for observed performance variations?
- Is there a particular version you preferred to write

As you assess, ensure that all code versions consistently produce similar energy values.

## Part 3 - Write up

Conclude by documenting your project's approach, implementation specifics, observations, and results from your performance evaluation. 
For each implementation, you should talk about decisions you made and why. For example, you will not need to rewrite the whole project using NumPy - which functions are most important to change and why?
Your report should span 4-8 pages, integrating figures and accompanying captions within the text body. You should also include a references section, which will not be included in the page count.
