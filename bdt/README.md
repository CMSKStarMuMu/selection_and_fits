Scripts to perform signal/background discrimination via bdt, single candidate per-event selection, and application of the B0Psi cut.  
The following steps should be performed:

1. Create the sub-samples, which are split based on the event N.
```
root -l -q 'create_subsample.C(int year = 2016, int type = 0, int mc = 0))'
```
wwhere year can be 2016, 2017 or 2018,  
type refers to the channel being investigated (defined [here](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/bdt/create_subsample.C#L26-L30))
and mc is 0 or 1 for data or MC respectively.  

2. Add some useful variables to the sub-samples, and skim the events keeping only those passing the pre-selection.
This is done in one step for the data, by running the add_vars_to_subsamples.py script, while for the MC samples the pre-selection skim 
has to be run with the further script apply_pass_preselection_on_mc.py. The reason behing this is to be able to rescale MC variables used for the pre-selection definition before actually applying the cut.
The 2 step process can be run using the ntuple_preparation_workflow.sh script.


3. BDT training/testing/WP optimisation

The optimisation of BDT hyper parameters is done with study_hp_optimisation/local_optimize_hp_punzi_2016.py and similar codes.
They require additional packages and an updated xgb version.
More detail will be added here.

The bdt training is then performed using the final_bdt_sub_samples_noIP2D_2016.py scripts (and similar).
The bdt score is added to the ntuples doing
```
python add_bdt_subsample.py YEAR
```

Scripts to determine the working point are in the wp_optimization folder.

wp_optimization/wp_optimization.py  (only needed when searching for the optimal bdt working point)  
          |  
wp_optimization/draw_wp_optimization.py  (only needed when searching for the optimal bdt working point)   
          |  

4. Once you know the WP to use, apply the selection on the samples, add the B0Psi cut and additional variables, doing:
```
source selectionWorkflow.sh
```


Please note that in each script there is a list defining which samples you would like to process, e.g. [here](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/bdt/add_bdt_subsamples.py#L23-L31), so you may wish to modify the list/add new samples there.
