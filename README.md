# sarscov2-roadmap-england

This is an [orderly](https://www.vaccineimpact.org/orderly/) repository which contains the analysis to our preprint

> Non-pharmaceutical interventions, vaccination and the Delta variant: epidemiological insights from modelling England's COVID-19 roadmap out of lockdown

## Running

A sequence of tasks needs to be run with a set of parameters to generate the final results.  This is sketched out in the [`run.R`](run.R) script, though this is provided only as a form of documentation. In practice these were run over several days on a HPC.

* regions: `c("north_west", "north_east_and_yorkshire", "midlands", "east_of_england", "london", "south_west", "south_east")`
* assumptions: `c("central", "optimistic", "pessimistic")`

1. Run the `vaccine_fits_regional` task with each region and assumption level (21 fits, each about 5 hours)
2. Run the `vaccine_fits_combined` task with each assumption level (3, collecting all regions)
3. Run the `vaccine_restart_fits_regional` task with each region and assumption level (21 fits, each about 10 hours)
4. Run the `vaccine_restart_fits_combined` task with each assumption level (3, collecting all regions)
5. Run the `vaccine_simulation` task with each assuption level, with a `restart_date` of `march`, `june` and `july` with `sensitivity=FALSE`, and also with `sensitivity=TRUE` for `restart_date="july"` . The full runs use `n_par = 200` but this can be reduced at the cost of more noise. Because of the large number of scenarios used, these were run across a set of 32 core nodes (up to 10 at a time).
6. Run the `vaccine_simulation_plots` task
7. Run the `vaccine_simulation_plots_sens` task

Running even the short run, as in `run.R` will take ~7 hours of CPU time, less in wall time if you have more cores available.

## Requirements

The core requirement is our [sircovid](https://mrc-ide.github.io/sircovid/) package and its dependencies. Because that package is in constant development you will probably want to pin your versions of the software to the versions we used for preparation:

```r
remotes::install_github(c(
  "mrc-ide/dust@v0.9.11",
  "mrc-ide/mcstate@0.6.6",
  "mrc-ide/sircovid@v0.11.30",
  "mrc-ide/spimalot@v0.3.0"))
```

However, you can always install the versions that we are using with

```r
drat:::add("ncov-ic")
install.packages(c("sircovid", "spimalot"))
```

You will also need a recent [orderly](https://www.vaccineimpact.org/orderly/) which can be installed with

```r
drat:::add("vimc")
install.packages("orderly")
```

To install _all_ required packages, you can use remotes:

```r
remotes::install_deps()
```

## License

MIT © Imperial College of Science, Technology and Medicine
