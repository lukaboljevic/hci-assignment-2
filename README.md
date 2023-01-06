# Brain-computer interface to classify between two motor activities

This repository contains the code for a brain-computer interface, whose goal is to classify between two motor activities. The motor activities were chosen subject S053 from the EEGMMI database, available on [PhysioNet](https://physionet.org/content/eegmmidb/1.0.0/).

We tried to classify between two motor activities from task 1 - opening and closing left or right fist, and then separately task 2 - imagine opening and closing left or right fist. We used the Common Spatial Patterns (CSP) method, and LDA and QDA classifiers.

A bash script to download the records for subject S053 is found under `data/`. If you're on Windows, run it from Git Bash. Currently, the script downloads all `.edf` and `.edf.event` files for S053, but it can be modified to download records and event files for any number of subjects.

The results are found under `results/S053`. The `txt` file names are of the format `SubjectID-task[1/2]-[order of band-pass FIR filter]-[feat(ure)v(ectors) / ref(erence class) / results].txt`.

To run the program, and obtain all possible results, simply run `run.m` from MATLAB. Of course, the input parameters can be changed to get results for whatever is needed.
