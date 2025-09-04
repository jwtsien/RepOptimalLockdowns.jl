# Replication of Acemoglu et al. (2021) in Julia
By Tanfei Li, Jiawei Qian, Tianxiang Yuan
## Overview
This repository contains Julia code that replicates the main results (Figures 3, 4, and 5) from the paper *Optimal Targeted Lockdowns in a Multi-Group SIR Model* by Acemoglu et al. (2021).

The code was originally developed by the group as part of the final project for the course Computational Economics in the Fall Semester of 2024.

## Running the code
1. Please make sure that you have julia
2. Clone the repository to your local machine
   ```
   git clone https://github.com/cwzie/RepOptimalLockdowns.jl
   ```
3. Navigate to the project directory
   ```REPL
   cd RepOptimalLockdowns.jl
   julia # enter Julia environment
   ```
4. Use Julia to activate the environment (use should be in Julia first)
   ```Julia
   # using Pkg
   ]activate . # ] is to enter the package environment
   ```
5. Delete "]" to return the Julia environment and run the code
   ```Julia
   using RepOptimalLockdowns # in Julia environment
   RepOptimalLockdowns.run()
   ```
**Tips:** If you meet the error
```Julia
ERROR: MethodError: no method matching create_my_directory()
```
just ignore it and run again.

## Result Structure
After running the code, a folder named `lockdown` will appear in the same directory (just in "/RepOptimalLockdowns.jl"), and all results will be saved in that folder.

- `lockdown/figs` it is use to store the plots.
- `lockdown/summaryres` it is used to store the final results.

## Difficulties and Problems
This part is to explain why we failed to get the corresponding output. 

Focused on file `src/GSiROptimalPolicy.jl`. We didn't find a great way to handle the first-order differential constraint in Julia. For example, in Python, we can use GEKKO,
```Python
model = GEKKO(remote = remote)
... # add some settings
# We require the rate of change of the variable $s$ with respect to time to be $rs$.
model.Equation(s.dt() == rs)
```
In Julia, I tried to use first-order difference constraint to replace it.
```Julia
# nt is the range of time and we define dt = 1
@constraint(model, [t in 1:nt-1], (s[t+1] - s[t]) == rs[t])
```
We should notice that in Python, $s$ is still a scalar but in Julia, I have to make $s$ a vector. It may not influence too much in a simple problem but our model is complex and it affects many equations (not only in terms of numerical values but also in terms of data types). 

After that, we tried different approaches to replicate the results. Finally, the program and the optimizer ran well but didn't export the same answer with the article.

Later, we realized that perhaps all the GEKKO variables should be converted into vector form to make them easier for Julia to handle.

## Contact
For questions or issues, please contact the group at [Jiawei](mailto:jiawei.qian@sciencespo.fr) or [Tanfei](mailto:tanfei.li@sciencespo.fr).

