@echo off
set PYTHONPATH=c:\Programming\PreFerS\Bioindustrial-Park;%PYTHONPATH%
"C:\Users\owenp\.conda\envs\BioSTEAM\python.exe" "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\LegHb\system\_config2.py" > log.txt 2>&1
echo DONE >> log.txt
