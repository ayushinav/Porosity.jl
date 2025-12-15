# Porosity.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ayushinav.github.io/Porosity.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ayushinav.github.io/Porosity.jl/dev/)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/ayushinav/Porosity.jl/Tests.yml)
[![codecov](https://codecov.io/gh/ayushinav/Porosity.jl/graph/badge.svg?token=VQM6W3DUI4)](https://codecov.io/gh/ayushinav/Porosity.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

## Rock physics modeling

`Porosity.jl` is a rock physics modeling and inference library written in julia, providing performant, AD compatible and scalable codes.

Check out the documentation [here](https://ayushinav.github.io/Porosity.jl/stable/)

`Porosity.jl` allows you to model geophysical observables from rock physics parameters., e.g., taking rock physics parameters from a model of mid-ocean ridge ([Sim et al., 2020
](https://www.sciencedirect.com/science/article/pii/S003192011930202X)) :

![](assets/sim_demo.png)

## Rock physics inference

The package also allows inference of rock physics parameters, while including different constraints (e.g. thermodynamics constraints). This allows you to obtain posteriors for these parameters (recreating Fig. 3 in [Blatter et al., 2022](https://www.nature.com/articles/s41586-022-04483-w)):

![](assets/blatter_demo.png)
