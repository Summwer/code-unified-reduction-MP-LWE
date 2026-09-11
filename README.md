# code for unified reduction from RLWE to MP-LWE

This repository contains the artifacts associated to the article _A unified reduction from RLWE to MP-LWE_ (that was used to run the experiments and produce the figures displayed in that article)

**A unified reduction from RLWE to MP-LWE**
by _Wenwen Xia_ based on the code [code-MPLWE](https://plmlab.math.cnrs.fr/apelletm/code-mplwe) implemented by _Rahinatou Yuh Njah Nchiwo and Alice Pellet-Mary_
ATTENTION: We remain the code `random_polynomials.sage`, `test_monogenic.sage`, `test_Vandermonde.sage` from the original repository, and create the new code `test_Hankel_matrix.sage`, `test_sigma_beta_equation.sage` and `test_min_norm_Hy.sage` in the new repository. 

We add the new code into `run_tests.sage` and `create_plots.sage` to run the new tests and create the new plots.


All the code from this repository is provided under the the GNU Affero Public License (see the `LICENSE.md` file).

## Requirements

* Tested with [SageMath](https://www.sagemath.org/) version 10.4 and 10.5

## Files organization

* `random_polynomials.sage` contains functions generating random polynomials with desired properties
* `test_monogenic.sage`, `test_Vandermonde.sage`, `test_Hankel_matrix.sage`, `test_sigma_beta_equation.sage` and `test_min_norm_Hy.sage` contain functions that performs all the tests described in the article
* `run_tests.sage` is the file that actually runs all the tests (in particular, it creates the list of inputs on which the tests are run)
* `create_plots.sage` uses the data obtained by running `run_tests.sage` and create nice pictures from it (that are displayed in the article)

* The folder `data` contains all the data obtained when running all the tests from `run_tests.sage`
* The folder `figures` contains all the figures obtained by running `create_plots.sage` 


## How to use

Theis code should be run in a directory containing two folders `data` and `figures`, where the data and the figure will be saved

Running all the experiments on a single core may take quite some times. With 10 parallel cores it may be run fully in approximately 8 to 12 hours.

### Generating the data

Make sure that you execute the code in a repository containing a `data` sub-repository.

In the file `run_tests.sage`, uncomment the lines containing the tests you want to run (these lines are preceded by a `##TODO` comment). When uncommenting, you can also update the quantity `nb_threads`, with the actual number of cores you can use for the computation (default is 5 or 3).

Then run in a terminal
```
   sage run_tests.sage
```

Running all the tests should take roughly 8 to 10 hours if you have 10 available cores (so roughly 3 to 5 days on a single core).

Note: the repository `data` provided here already contains all the data that can be generated using the script `run_tests.sage`. The random seed has been fixed in all the experiments, so re-running the code as explained above should produce exactly the same data files (except for the files starting with `norm_conductor`, which contain timings that may vary with each run).

### Creating the plots

Make sure that you execute the code in a repository containing a `figures` folder.

Also make sure that the file `run_tests.sage` is in its original format, with all the lines executing the tests that are commented. Because the file `create_plots.sage` starts by loading the file `run_tests.sage`.

Once you have obtained the data you need using the file `run_tests.sage` as above (or using the data provided in the repository `data`), you can create the figures using the script in `create_plots.sage`. To do so, uncomment the lines containing the figures you want to create (these lines are preceded by a `##TODO` comment). Then run in a terminal
```
   sage create_plots.sage
```

This should take roughly 2 to 3 minutes.

Note: the repository `figures` provided here already contains all the figures that can be created using the script `create_plots.sage`.
