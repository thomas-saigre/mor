= Compute the Sobol indices for the sensitivity analysis

== Set the offline model

Using the application [`feelpp_mor_toolbox_heat_offline`](https://github.com/feelpp/feelpp/blob/develop/mor/apps/toolboxmor_heat_offline.cpp)



== Run the sensitivity analysis

```bash
./feelpp_mor_sensitivity_analysis --crbmodel.name <model-name> --sampling.size <size>
```