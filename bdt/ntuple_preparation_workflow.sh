#!/bin/tcsh
source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

set year = 2016   

echo "python add_vars_to_subsamples.py ${year} "
python add_vars_to_subsamples.py ${year}
echo "done"

echo "python apply_pass_preselection_on_mc.py ${year} "
python apply_pass_preselection_on_mc.py ${year}
echo "done"