# epicookbook

Model code from the [Epicookbook](http://epirecip.es/epicookbook/). The original Jupyter Notebooks can be found [here](https://github.com/epirecipes/epicookbook/tree/master/notebooks). Code is provided in multiple languages and libraries, including R, Python, Matlab/Octave, C, Fortran, and others.

The directory structure is provided below using the Bash `tree` command:

```bash
├── README.md
├── applications
│   ├── acute_hiv.r
│   ├── deterministic_seir_ebola_pygom.py
│   ├── deterministic_seir_ebola_scipy.py
│   ├── hiv_2_risk_groups.r
│   └── stochastic_seasonal_discrete_t_rotavirus.r
├── bjornstad_2018
│   ├── chapter1.r
│   ├── chapter2_sir.r
│   ├── chapter3_r0.r
│   ├── chapter4_foi_age_incidence.r
│   └── chapter5_seasonality.r
├── epidemic_final_size
│   ├── ere_si4r.jl
│   ├── ere_si4r.m
│   ├── ere_si4r.sce
│   ├── ere_sir.jl
│   ├── ere_sir.m
│   ├── ere_sir.r
│   └── ere_sir.sce
├── host_vector_models
│   ├── 1host_1vector.jl
│   ├── nhost_1vector.jl
│   └── nhost_mvector.jl
├── keeling_rohani_2008
│   ├── program_2_1.c
│   ├── program_2_1.f90
│   ├── program_2_1.m
│   ├── program_2_1.py
│   ├── program_2_6_seir.jl
│   ├── program_2_6_seir.r
│   ├── program_3_1_sis.jl
│   ├── program_3_1_sis.r
│   ├── program_3_2_sis.r
│   ├── program_3_4_age_structured_seir.r
│   └── program_4_4_multi_seir.r
├── metapopulation_models
│   ├── deterministic_seir.r
│   └── large_population_sirs.jl
├── microparasite_models
│   ├── may_anderson_1978.jl
│   └── may_anderson_1978.r
├── network_models
│   ├── edge_based_sir.js
│   ├── edge_based_sir.r
│   ├── london_sir.r
│   └── pneumococcal_ibm.r
├── nonexponential_passage_times
│   ├── discrete_erlang_sir.jl
│   └── discrete_erlang_sir.r
├── phylodynamics
│   └── simple_coalescent.r
├── simple_deterministic_models
│   ├── scaling
│   │   ├── scaling.jl
│   │   ├── scaling.py
│   │   └── scaling_desolve.r
│   ├── seir
│   │   ├── seir.jl
│   │   └── seir_desolve.r
│   ├── sir
│   │   ├── sir.cpp
│   │   ├── sir.jl
│   │   ├── sir.js
│   │   ├── sir.m
│   │   ├── sir.py
│   │   ├── sir.sce
│   │   ├── sir.vf
│   │   ├── sir.xpp
│   │   ├── sir_desolve.r
│   │   └── sir_odin.r
│   └── sis
│       └── sis.js
├── simple_stochastic_models
│   ├── continuous_t_sir
│   │   ├── ssir.jl
│   │   ├── ssir.r
│   │   ├── ssir_gillespie.r
│   │   └── ssir_rcpp.r
│   ├── discrete_t_seird
│   │   └── seird_odin.r
│   └── discrete_t_sir
│       ├── ssir.jl
│       ├── ssir.py
│       ├── ssir_libbi.r
│       ├── ssir_odin.r
│       └── ssir_pomp.r
└── time_varying_parameters
    ├── seasonally_forced_deterministic.js
    ├── seasonally_forced_deterministic.r
    ├── semiparametric_sir.jl
    └── semiparametric_sir.r
```
