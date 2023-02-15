#!/bin/tcsh
source  /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scram runtime -csh`

set year = 2016   
set tag_string='_noIP2D'

echo "python selectSingleCand.py ${year} --tag ${tag_string}"
python selectSingleCand.py ${year} --tag ${tag_string}
echo "done"

echo " ******************************************************* " 
echo "python applyB0PsiCut.py ${year} --tag ${tag_string}"
python applyB0PsiCut.py ${year} --tag ${tag_string}
echo "done"

echo " ******************************************************* " 
foreach i ( 'MC_JPSI' 'MC_LMNR' 'MC_PSI' 'data' )
# foreach i ( 'MC_JPSI' 'MC_LMNR' 'MC_PSI' 'MC_BS' 'MC_BSJPSIPHI' 'MC_BSJPSIKST' 'data' )
    echo "hadd ../final_ntuples/${year}${i}${tag_string}_addxcutvariable.root ../final_ntuples/${year}${i}${tag_string}_${year}_part*_addB0Psi.root"
    hadd ../final_ntuples/${year}${i}${tag_string}_addxcutvariable.root ../final_ntuples/${year}${i}${tag_string}_${year}_part*_addB0Psi.root
end

