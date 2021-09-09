# Neural-Cascade
<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li>
      <a href="#requirements">Requirements</a>
    </li>
    <li>
      <a href="#usage">Usage</a>
    </li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## Overview
This repository is the implementation code of the paper "Single neuron firing cascades underlie global spontaneous brain events". For the analysis, we used neural and behavioral signals recorded from mice by the [Visual Coding - Neuropixels project](https://portal.brain-map.org/explore/circuits/visual-coding-neuropixels) of the Allen Brain Observatory at the Allen Institute. 

<!-- GETTING STARTED -->
## Requirements

To successfully run the code and reproduce the results, the following packages are required.
* Python 3
* numpy
* scipy 1.5.0
* allensdk 2
* matplotlib
* seaborn

## Usage
1. Download the neural recording data. Please refer to the instruction and scripts for accessing the data from the [Visual Coding project websites](https://allensdk.readthedocs.io/en/latest/visual_coding_neuropixels.html).
2. Download the mouse connectivity data. Please refer to the instruction and scripts for accessing the data from the [Mouse Connectivity websites](https://allensdk.readthedocs.io/en/latest/connectivity.html).
3. Specify the directory where the datasets are downloaded (also where the manifest file is saved) in `settings.py`
   ```py
   observatory.data_directory = Path('path-to-brain-observatory-dataset')
   connectivity.data_directory = Path('path-to-mouse-connectivity-dataset')
   ```
   and speficy the directory where you want to save processed data and intermediate results in `settings.py`
   ```py
   projectData.dir.root = Path('path-to-save-project-data')
   ```
 4. Run `main_proc.py` to process the data.
 5. Run `results_example.ipynb` and `results_group.ipynb` for data analysis.
