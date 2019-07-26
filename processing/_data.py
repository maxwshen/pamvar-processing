import _config
import sys, os, pickle, subprocess
import pandas as pd

exp_design = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv', index_col = 0)
treatment_control_design = pd.read_csv(_config.DATA_DIR + 'treatment_control_design.csv', index_col = 0)

prj_dir = '/ahg/regevdata/projects/CRISPR-libraries/prj/'

L1_names = set(exp_design['Name'])

##
# Helper functions
##
def check_valid_and_fetch_level(nm):
  if nm in L1_names:
    row = exp_design[exp_design['Name'] == nm].iloc[0]
    return 'L1', row
  elif nm in L2_names:
    row = group_design[group_design['Name'] == nm].iloc[0]
    return 'L2', row
  elif nm in endo_names:
    row = endo_design[endo_design['Name'] == nm].iloc[0]
    return 'endo', row
  else:
    # print('Bad name: %s' % (nm))
    raise ValueError
    return None, None

##
# Data fetching
##
def check_file_size(nm, data_nm):
  level, row = check_valid_and_fetch_level(nm)
  if level is None:
    return None

  scripts_base = [
    'g4_poswise_be',
    'g5_combin_be',
    'h4_poswise_del',
    'h6_anyindel',
  ]
  if data_nm in scripts_base:
    fn = prj_dir + '%s/out/%s/%s.pkl' % (row['Group'], data_nm, row['Local name'])

  if data_nm in ['e_newgenotype_Cas9']:
    fn = prj_dir + '%s/out/%s/%s.csv' % (row['Group'], data_nm, row['Local name'])


  scripts_adjust = [
    'ag4_poswise_be_adjust',
    'ag5_combin_be_adjust',
    'ae_newgenotype_Cas9_adjust',
    'ag4b_adjust_background',
    # 'ag4e_find_highfreq_batcheffects',
    'ag4a2_adjust_batch_effects',
    'ag4b4_filter_inprofile_batch_effects',
    'ag5a1a_hf_bc',
    'ag5a1b_subtract',
    'ag5a1c_ie',
    'ag5a2_adjust_batch_effects',
    'ag5a3_adjust_inprofile_batch_effects',
    'ag5a4_profile_subset',
    'ah6a1a_hf_bc',
    'ah6a1b_subtract',
    'ah6a3_remove_batch',
  ]
  if data_nm in scripts_adjust:
    if level == 'L1':
      fn = _config.OUT_PLACE + '%s/%s.pkl' % (data_nm, nm)
    elif level == 'L2':
      fn = _config.OUT_PLACE + 'b_group_L2/%s/%s.pkl' % (data_nm, nm)

  file_size = os.path.getsize(fn)
  mb_file_size = file_size / 1000000

  return mb_file_size



def load_data(nm, data_nm):
  '''
    Supported nm: An item in 
      exp_design['Name']
      group_design['Name']
      endo_design['Name']

    Supported data_nm: [
      'g4_poswise_be', 
      'g5_combin_be', 
      'e_newgenotype_Cas9'
      'a4_poswise_adjust',
      'a5_combin_adjust',
      'a_newgenotype_Cas9_adjust',
    ]
  '''

  level, row = check_valid_and_fetch_level(nm)
  if level is None:
    return None

  scripts_base = [
    'g4_poswise_be', 
    'g5_combin_be', 
    'h4_poswise_del',
    'h6_anyindel',
  ]
  if data_nm in scripts_base:
    fn = prj_dir + '%s/out/%s/%s.pkl' % (row['Group'], data_nm, row['Local name'])
    with open(fn, 'rb') as f:
      data = pickle.load(f)

  if data_nm in ['e_newgenotype_Cas9']:
    fn = prj_dir + '%s/out/%s/%s.csv' % (row['Group'], data_nm, row['Local name'])
    data = pd.read_csv(fn, index_col = 0)

  scripts_adjust = [
    'ag4_poswise_be_adjust', 
    'ag5_combin_be_adjust', 
    'ae_newgenotype_Cas9_adjust',
    'ag4b_adjust_background',
    # 'ag4e_find_highfreq_batcheffects',
    'ag4a2_adjust_batch_effects',
    'ag4b4_filter_inprofile_batch_effects',
    'ag5a1a_hf_bc',
    'ag5a1b_subtract',
    'ag5a1c_ie',
    'ag5a2_adjust_batch_effects',
    'ag5a3_adjust_inprofile_batch_effects',
    'ag5a4_profile_subset',
    'ah6a1a_hf_bc',
    'ah6a1b_subtract',
    'ah6a3_remove_batch',
  ]
  if data_nm in scripts_adjust:
    if level == 'L1':
      fn = _config.OUT_PLACE + '%s/%s.pkl' % (data_nm, nm)
    elif level == 'L2':
      fn = _config.OUT_PLACE + 'b_group_L2/%s/%s.pkl' % (data_nm, nm)
    with open(fn, 'rb') as f:
      data = pickle.load(f)

  return data

def load_minq(nm, data_nm):
  '''
    Supported nm: An item in 
      exp_design['Name']
    
    (L1 only)

    Supported data_nm: [
      'g4_poswise_be', 
      'g5_combin_be', 
    ]
  '''

  if data_nm not in ['g4_poswise_be', 'g5_combin_be']:
    return None

  level, row = check_valid_and_fetch_level(nm)
  if level is None or level == 'L2':
    return None

  fn = prj_dir + '%s/out/%s/%s_minq.pkl' % (row['Group'], data_nm, row['Local name'])
  with open(fn, 'rb') as f:
    minq = pickle.load(f)

  return minq

##
# Data queries
##
def get_lib_design(nm):
  lib_nm = get_lib_nm(nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm))
  seq_col = [s for s in lib_design.columns if 'Sequence context' in s][0]
  return lib_design, seq_col

def get_l1_nm(nm):
  level, row = check_valid_and_fetch_level(nm)
  if level is None:
    return None

  if level == 'L1':
    l1_nm = nm
  elif level == 'L2':
    row = group_design[group_design['Name'] == nm].iloc[0]
    l1_nm = row['Item1']
  return l1_nm

def get_editor_type(nm):
  l1_nm = get_l1_nm(nm)

  row = exp_design[exp_design['Name'] == l1_nm].iloc[0]
  editor_nm = row['Editor']

  editor_types = {
    'CtoTeditor': [
      'AID',
      'CDA',
      'BE4',
      'evoAPOBEC',
      'BE4-CP1028',
      'eA3A',
    ],
    'AtoGeditor': [
      'ABE',
      'ABE-CP1040',
    ],
    'Cas9': [
      'Cas9',
    ]
  }

  for typ in editor_types:
    if editor_nm in editor_types[typ]:
      return typ
  return None

def get_editor_nm(nm):
  l1_nm = get_l1_nm(nm)
  row = exp_design[exp_design['Name'] == l1_nm].iloc[0]
  editor_nm = row['Editor']
  return editor_nm

def get_lib_nm(nm):
  try:
    l1_nm = get_l1_nm(nm)
    row = exp_design[exp_design['Name'] == l1_nm].iloc[0]
    lib_nm = row['Library']
  except ValueError:
    brute_success = False
    for lib_nm_temp in zero_pos:
      if lib_nm_temp in nm:
        lib_nm = lib_nm_temp
        brute_success = True
    if not brute_success:
      raise ValueError
  return lib_nm

def get_ontarget_sites(lib_design, lib_nm):
  if lib_nm == '12kChar':
    otcats = ['satmut']
  elif lib_nm == 'AtoG':
    otcats = list(set(lib_design['Design category']))
  elif lib_nm == 'CtoT':
    otcats = list(set(lib_design['Design category']))
  elif lib_nm == 'LibA':
    otcats = list(set(lib_design['Design category']))

  return lib_design[lib_design['Design category'].isin(otcats)]['Name (unique)']

def get_replicates(celltype, lib_nm, editor_nm):
  crit = {
    'Cell-type': celltype, 
    'Library': lib_nm, 
    'Editor': editor_nm
  }
  match = lambda row: bool(sum([row[col_nm] == crit[col_nm] for col_nm in crit]) == len(crit))
  reps = exp_design[exp_design.apply(match, axis = 'columns')]['Name']
  return list(reps)

# Caching
lib_nms = dict()

zero_pos = {
  'PAMvar': 12,
}

def idx_to_pos(idx, nm):
  if nm not in lib_nms:
    lib_nms[nm] = get_lib_nm(nm)
  lib_nm = lib_nms[nm]
  return idx - zero_pos[lib_nm]

def pos_to_idx(pos, nm):
  if nm not in lib_nms:
    try:
      lib_nms[nm] = get_lib_nm(nm)
    except ValueError:
      brute_success = False
      for lib_nm in zero_pos:
        if lib_nm in nm:
          lib_nms[nm] = lib_nm
          brute_success = True
      if not brute_success:
        raise ValueError
  lib_nm = lib_nms[nm]
  return pos + zero_pos[lib_nm]
