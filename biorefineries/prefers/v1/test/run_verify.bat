@echo off
set PYTHONPATH=c:\Programming\PreFerS\Bioindustrial-Park;%PYTHONPATH%
"C:\Users\owenp\.conda\envs\BioSTEAM\python.exe" "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\test\verify_upgrade.py" > verify_output.log 2>&1
echo DONE >> verify_output.log
