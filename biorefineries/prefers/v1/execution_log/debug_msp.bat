@echo off
set PYTHONPATH=C:\Programming\PreFerS\Bioindustrial-Park;%PYTHONPATH%
"C:\Users\owenp\.conda\envs\BioSTEAM\python.exe" "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\execution_log\run_msp_comparison.py" > "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\execution_log\results_msp.txt" 2>&1
echo DONE >> "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\execution_log\results_msp.txt"
