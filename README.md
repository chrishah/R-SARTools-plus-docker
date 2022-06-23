# R-SARTools-plus-docker
R 4.0.3 setup with SARTools plus necessary dependencies


See `Dockerfile` for the installed packages.

To get straight into an R console:
```bash
docker run --rm -it -v $(pwd):/in -w /in chrishah/r-sartools-plus:2b95eaa /bin/bash
```

To run an Rscript through this R setup:
```bash
docker run --rm -v $(pwd):/in -w /in chrishah/r-sartools-plus:2b95eaa Rscript my_R_analyses.R
```
