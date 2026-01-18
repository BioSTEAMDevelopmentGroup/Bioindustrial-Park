@echo off
set PYTHONPATH=c:\Programming\PreFerS\Bioindustrial-Park;%PYTHONPATH%
"C:\Users\owenp\.conda\envs\BioSTEAM\python.exe" "c:\Programming\PreFerS\Bioindustrial-Park\biorefineries\prefers\v1\LegH\_system.py" > c:\Programming\PreFerS\debug_output.txt 2>&1
echo Simulation finished. >> c:\Programming\PreFerS\debug_output.txt
