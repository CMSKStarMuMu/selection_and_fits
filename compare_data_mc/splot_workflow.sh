#!/bin/tcsh
source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

set year = 2017
set channel = 'keepJpsi'
set input_file = '../final_ntuples/2017data_noIP2D_noNan_addxcutvariable.root'

echo "python prepareSPlot_postBDT.py ${input_file} ${channel}"
python prepareSPlot_postBDT.py ${input_file} ${channel}
echo "done"

echo " ******************************************************* " 
echo "python apply_splot_cuts.py ${year} ${input_file}"
python apply_splot_cuts.py ${year} ${input_file}
echo "done"

echo " ******************************************************* " 
echo "python merge_sweights_data.py ${year} ${input_file}"
python merge_sweights_data.py ${year} ${input_file}
echo "done"

