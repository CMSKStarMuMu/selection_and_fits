## Workflow to compare data (normalisation channel) to signal MC  

- [prepareSPlot_postBDT.py](): script to fit the data mass distribution and save the splot weights. If run with the --save_workspace options, only saves the fitted data and functions into a workspace.  
- [apply_splot_cuts.py.py](): skim the original ntuple applying the same cuts used in the sWeights calculation
- [merge_sweights_data.py](): merge the skimmed ntuple with the sWeights.  

The whole workflow can be run doing  
```
./splot_workflow.sh
```
where the era/year should be set.  

The [plots_from_workspace.py]() script can be used to prepare a pdf of the data/fit function used to produce the sWeights.  
The [plotDistributions_controlChannel_perEra.py]() script creates plots comparing data to MC for the different data-taking eras/years/L1 selections.  
The [plotDistributions_controlChannel_afterMCw.py]() script creates plots comparing data to MC and weighted MC distribution.  
