# Template Fitting window testing for the Mu2e Calorimeter

When a cosmic ray hits one of the calorimeter crystals, it produces two signals (one for each SiPM), that are sampled every 5 ns.
Caloreco.C fits each signal with a template function over a certain time interval (start and end are relative to the peak time). The timing of the event is then defined as when the fit reaches 20 % of the maximum. Now AnaDriver2_1 can compute the timing difference between the two signals from each crystal over many events, we take the simga of this timing difference and compute its mean over the crystals, then we plot it against the interval.

### To use

Intervals to be tested are in t_intervals.py (measured in ns and relative to the peack).
To run the analysis use time_tests.ipynb (preferably on the [EAF](https://mu2ewiki.fnal.gov/wiki/Elastic_Analysis_Facility_(EAF)) for imporved parallelism).

Using v_fit_tries.py should also work fine.

Data are read from the files specified in data.list and results are created in the direcoty specified in the script, for time intervals that where already analized, the program recognizes the output file and reads it. The opuput graph is saved in the results directory as well as a CSV with the plotted values.
