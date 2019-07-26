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
inp_dir = _config.OUT_PLACE + 'be_combin_general_annot/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

import _data

nts = list('ACGT')
nt_to_idx = {nts[s]: s for s in range(len(nts))}

treatment_control_design = pd.read_csv(_config.DATA_DIR + 'treatment_control_design.csv', index_col = 0)


##
# Primary logic
##
def anyedit_collate(exp_nm):

  # Load data
  data = pd.read_csv(inp_dir + '%s.csv' % (exp_nm), index_col = 0)

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

  data = data.drop(columns = ['Unnamed: 0.1', 'Design category', 'Fraction CBE aberrant edit', 'CBE aberrant count'])

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

  pv_df = data.pivot(index = 'Name (unique)', columns = 'gRNA properties', values = 'Fraction edited')
  pv_df.to_csv(out_dir + f'5primeG_{exp_nm}.csv')

  return

##
# Main
##
@util.time_dec
def main():
  print(NAME)
  
  for exp_nm in treatment_control_design['Treatment']:
    print(exp_nm)
    anyedit_collate(exp_nm)

  return


if __name__ == '__main__':
  main()