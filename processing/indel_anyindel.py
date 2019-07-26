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

# Default params
inp_dir = _config.OUT_PLACE
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

import _data

nts = list('ACGT')
nt_to_idx = {nts[s]: s for s in range(len(nts))}

treat_control_df = pd.read_csv(_config.DATA_DIR + 'treatment_control_design.csv', index_col = 0)
exp_design_df = pd.read_csv(_config.DATA_DIR + 'master_exp_design.csv', index_col = 0)


##
# Primary logic
##
def indel_anyindel(exp_nm):
  try:
    data = _data.load_data(exp_nm, 'ah6a1b_subtract')
  except:
    print('Error : could not load data')
    sys.exit(1)

  lib_design, seq_col = _data.get_lib_design(exp_nm)

  dd = defaultdict(list)
  timer = util.Timer(total = len(data))
  for target_nm in data:
    df = data[target_nm]

    tot_count = sum(df['Count'])
    dd['Total count'].append(tot_count)
    dd['Name (unique)'].append(target_nm)

    crit = (df['Category'] != 'wildtype')
    indel_count = sum(df[crit]['Count'])
    dd['Indel count'].append(indel_count)
    if tot_count != 0:
      dd['Indel freq'].append(indel_count / tot_count)
    else:
      dd['Indel freq'].append(np.nan)

    crit = (df['Category'] == 'del')
    del_count = sum(df[crit]['Count'])
    dd['Del count'].append(del_count)
    if tot_count != 0:
      dd['Del freq'].append(del_count / tot_count)
    else:
      dd['Del freq'].append(np.nan)

    crit = (df['Category'] == 'ins')
    ins_count = sum(df[crit]['Count'])
    dd['Ins count'].append(ins_count)
    if tot_count != 0:
      dd['Ins freq'].append(ins_count / tot_count)
    else:
      dd['Ins freq'].append(np.nan)

    crit = (df['Category'] == 'wildtype')
    wt_count = sum(df[crit]['Count'])
    dd['Wildtype count'].append(wt_count)
    if tot_count != 0:
      dd['Wildtype freq'].append(wt_count / tot_count)
    else:
      dd['Wildtype freq'].append(np.nan)

    timer.update()

  df = pd.DataFrame(dd)

  data = lib_design.merge(
    df, 
    on = 'Name (unique)', 
    how = 'outer',
  )

  # Annotate csvs
  pam_start_idx = 33
  pam_len = 5
  get_pam = lambda row: row['Sequence context (61nt)'][pam_start_idx : pam_start_idx + pam_len]

  data['Designed PAM (5nt)'] = data.apply(get_pam, axis = 'columns')

  get_true_grna_len = lambda row: 20 if row['gRNA (20nt)'][0] == 'G' else 21
  data['True gRNA length'] = data.apply(get_true_grna_len, axis = 'columns')

  grna_pos1_idx = 13
  grna_pos0_idx = grna_pos1_idx - 1
  grna_5primeG_matches_target = lambda row: bool(row['Sequence context (61nt)'][grna_pos1_idx] == 'G') if bool(row['True gRNA length'] == 20) else bool(row['Sequence context (61nt)'][grna_pos0_idx] == 'G') 
  data['gRNA 5primeG matches target'] = data.apply(grna_5primeG_matches_target, axis = 'columns')

  # Prepare data
  # data = data[data['Total count'] >= 100]

  # # Gather statistics

  def grna_5primeg_and_len(row):
    grna_len = row['True gRNA length']
    match = row['gRNA 5primeG matches target']
    if match:
      return f'{grna_len}-nt gRNA, 5primeG matches'
    else:
      return f'{grna_len}-nt gRNA, 5primeG does not match'

  data['gRNA properties'] = data.apply(grna_5primeg_and_len, axis = 'columns')
  data.to_csv(out_dir + f'{exp_nm}.csv')

  data = data[data['Total count'] >= 100]

  pv_df = data.pivot(index = 'Name (unique)', columns = 'gRNA properties', values = 'Indel freq')
  pv_df.to_csv(out_dir + f'5primeG_{exp_nm}.csv')


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

  num_scripts = 0
  for idx, row in treat_control_df.iterrows():
    treat_nm = row['Treatment']

    if 'Cas9' not in treat_nm:
      continue

    out_fn = out_dir + f'{treat_nm}.csv'
    if os.path.isfile(out_fn):
      continue

    command = 'python %s.py %s' % (NAME, treat_nm)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, treat_nm)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -V -P regevlab -l h_rt=4:00:00,h_vmem=1G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

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
def main(exp_nm = ''):
  print(NAME)
  
  indel_anyindel(exp_nm)


  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(exp_nm = sys.argv[1])
  else:
    gen_qsubs()
    # main()