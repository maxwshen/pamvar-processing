# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, pickle
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd
import sklearn
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr, spearmanr

from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor

# Default params
inp_dir = _config.OUT_PLACE + 'ag5a4_profile_subset/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

import _data

nts = list('ACGT')
nt_to_idx = {nts[s]: s for s in range(len(nts))}

treatment_control_design = pd.read_csv(_config.DATA_DIR + 'treatment_control_design.csv', index_col = 0)

not_nt_cols = set([
  'Count',
  'Num. edits',
  'Has aberrant CBE edit',
])

##
# Primary logic
##
def editing_profile(exp_nm, start_idx, end_idx):

  # Load data
  print('Loading data...')
  with open(inp_dir + f'{exp_nm}_{start_idx}_{end_idx}.pkl', 'rb') as f:
    data = pickle.load(f)
  print('Done.')

  pw_dd = defaultdict(list)

  timer = util.Timer(total = len(data))
  for target_nm in data:
    df = data[target_nm]

    nt_cols = [col for col in df.columns if col not in not_nt_cols]

    tot_count = sum(df['Count'])
    if tot_count <= 200:
      continue

    local_dd = defaultdict(lambda: defaultdict(lambda: 0))
    for idx, row in df.iterrows():
      count = row['Count']
      for col in nt_cols:
        ref_nt = col[0]
        pos = int(col[1:])
        obs_nt = row[col]

        if obs_nt == '.':
          obs_nt = ref_nt

        local_dd[col][obs_nt] += count

    for col in local_dd:
      ref_nt = col[0]
      pos = int(col[1:])
      for obs_nt in nts:
        pw_dd['Obs nt'].append(obs_nt)
        pw_dd['Ref nt'].append(ref_nt)
        pw_dd['Mutation'].append(f'{ref_nt}->{obs_nt}')
        pw_dd['Frequency'].append(local_dd[col][obs_nt] / tot_count)
        pw_dd['Target site'].append(target_nm)
        pw_dd['Position'].append(pos)
    timer.update()

  pw_df = pd.DataFrame(pw_dd)
  pw_df.to_csv(out_dir + f'{exp_nm}_{start_idx}_{end_idx}.csv')

  return

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print('Generating qsub scripts...')
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  # Generate qsubs only for unfinished jobs
  treat_control_df = pd.read_csv(_config.DATA_DIR + 'treatment_control_design.csv', index_col = 0)

  num_scripts = 0
  for idx, row in treat_control_df.iterrows():
    treat_nm = row['Treatment']
    if 'Cas9' in treat_nm:
      continue

    lib_nm = _data.get_lib_nm(treat_nm)

    num_targets = 12000
    num_targets_per_split = 2000

    for start_idx in range(0, num_targets, num_targets_per_split):
      end_idx = start_idx + num_targets_per_split - 1

      # Skip completed
      # out_pkl_fn = out_dir + '%s_%s_%s.pkl' % (treat_nm, start_idx, end_idx)
      # if os.path.isfile(out_pkl_fn):
      #   if os.path.getsize(out_pkl_fn) > 0:
      #     continue

      command = 'python %s.py %s %s %s' % (NAME, treat_nm, start_idx, end_idx)
      script_id = NAME.split('_')[0]

      try:
        mb_file_size = _data.check_file_size(treat_nm, 'ag5a4_profile_subset')
      except FileNotFoundError:
        mb_file_size = 0
      ram_gb = 2
      if mb_file_size > 140:
        ram_gb = 4
      if mb_file_size > 400:
        ram_gb = 8
      if mb_file_size > 1000:
        ram_gb = 16

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, treat_nm, start_idx)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -V -l h_rt=4:00:00,h_vmem=%sG -wd %s %s &' % (ram_gb, _config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return


##
# Main
##
@util.time_dec
def main(argv):
  print(NAME)
  
  exp_nm = argv[0]
  start_idx = int(argv[1])
  end_idx = int(argv[2])

  editing_profile(exp_nm, start_idx, end_idx)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()