from pathlib import Path
import os
import sys
import inspect
from utils import *

settings = DataContainer()

""" project general settings """
settings.homeDir = Path('/home/xxxx')
settings.projectDir = Path(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))

""" dataset: allen observotory """
observatory = DataContainer()
observatory.data_directory = Path('path-to-brain-observatory-dataset')
observatory.manifest_path = observatory.data_directory / 'manifest.json'

observatory.sess_removed = [767871931, 768515987, 778240327, 786091066, 787025148, 789848216, 
                            819701982, 831882777, 835479236, 839068429, 839557629, 840012044]

observatory.sess_few_TH_units = [774875821, 829720705]

settings.observatory = observatory

""" dataset: mouse connectivity """
connectivity =  DataContainer()
connectivity.data_directory = Path('path-to-mouse-connectivity-dataset')
connectivity.manifest_path = connectivity.data_directory / 'manifest.json'

settings.connectivity = connectivity

""" project data """
projectData = DataContainer()
projectData.dir = DataContainer()

projectData.dir.root = Path('path-to-save-project-data')
projectData.dir.sessions = projectData.dir.root / 'sessions'
projectData.dir.groupLevel = settings.homeDir / 'group-level'

files = DataContainer()
files.sessions = AttrDict(
    sessData = 'sess_data.pkl',

    allData_oldVer = 'all_data_old_ver.pkl',
    
    rsSpikeData_mfile = 'rs_spike_data.mat',
    cohLowRes_mfile = 'coherence_low_res.mat',
    cohHighRes_mfile = 'coherence_high_res.mat',

    rsData1 = 'rs_data_1.pkl',
    rsData2 = 'rs_data_2.pkl',

    flData = 'fl_data.pkl',
    dgrData = 'dgr_data.pkl',
    dgcData = 'dgc_data.pkl',

    lfpData = 'lfp_data.pkl'
)

projectData.files = files
settings.projectData = projectData

""" results """
resultData = DataContainer()
resultData.dir = DataContainer()
resultData.dir.root = settings.projectDir / 'results'

files.spatial_org = AttrDict(
    delay_map = 'spatial_delay_map.pdf',

    region_delay = 'spatial_region_delay.pdf',
    region_pos_neg_ratio = 'spatial_region_pos_neg_ratio.pdf',
    visual_delay = 'spatial_visual_delay.pdf',
    layer_pos_neg_ratio = 'spatial_layer_pos_neg_ratio.pdf',

    region_delay2 = 'spatial_region_delay2.pdf',
    region_pos_neg_ratio2 = 'spatial_region_pos_neg_ratio2.pdf',
    visual_delay2 = 'spatial_visual_delay2.pdf',
    layer_pos_neg_ratio2 = 'spatial_layer_pos_neg_ratio2.pdf',

    delay_consistency_hip = 'spatial_delay_consistency_hip.pdf',
    delay_consistency_vis = 'spatial_delay_consistency_vis.pdf',
    delay_consistency_tha = 'spatial_delay_consistency_tha.pdf'
)

files.coherence = AttrDict(
    coh_stationary_group = 'coherence_coh_stationary_group.pdf',
    coh_stationary_group_zoom = 'coherence_coh_stationary_group_zoom.pdf',
    coh_running_group = 'coherence_coh_running_group.pdf',
    coh_running_group_zoom = 'coherence_coh_running_group_zoom.pdf',
    coh_two_state_group = 'coherence_coh_two_state_group.pdf',
    coh_two_state_group_zoom = 'coherence_coh_two_state_group_zoom.pdf',

    coh_stationary_exmpl = 'coherence_coh_stationary_exmpl.pdf',
    coh_stationary_exmpl_zoom = 'coherence_coh_stationary_exmpl_zoom.pdf',
    coh_running_exmpl = 'coherence_coh_running_exmpl.pdf',
    coh_running_exmpl_zoom = 'coherence_coh_running_exmpl_zoom.pdf',
    coh_two_state_exmpl = 'coherence_coh_two_state_exmpl.pdf',
    coh_two_state_exmpl_zoom = 'coherence_coh_two_state_exmpl_zoom.pdf',

    scatter_phase_delay_group = 'coherence_scatter_phase_delay_group.pdf',
    boxplot_phase_delay_group = 'coherence_boxplot_phase_delay_group.pdf',

    scatter_phase_delay_exmpl = 'coherence_scatter_phase_delay_exmpl.pdf',

    corr_structure_stationary = 'correlation_structure_stationary.pdf',
    corr_structure_stationary_diff = 'correlation_structure_stationary_diff.pdf'
)

files.delay = AttrDict(
    pca_results = 'pca_results.pdf'
)

resultData.files = files
settings.resultData = resultData
