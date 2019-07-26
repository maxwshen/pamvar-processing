# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, math, pickle, imp
sys.path.append('/home/unix/maxwshen/')
import fnmatch
import numpy as np
from collections import defaultdict
from mylib import util
from mylib import compbio
import pandas as pd

# Default params
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

##
# Alignment manipulations
##
def trim_start_end_dashes(seq):
  alphabet = set(list(seq))
  if '-' in alphabet:
    alphabet.remove('-')
  alphabet = list(alphabet)
  start = min([seq.index(s) for s in alphabet])
  end = min([seq[::-1].index(s) for s in alphabet])
  if end > 0:
    return seq[start:-end]
  else:
    return seq[start:]

def parse_header(header):
  count = int(header.split('_')[0].replace('>', ''))
  return count

def count_indels(read, ref):
  num_dels = len(read.replace('-', ' ').split()) - 1
  num_ins = len(ref.replace('-', ' ').split()) - 1
  return num_dels, num_ins

def get_single_indel_range(seq):
  # Returns indices such that seq[indel_start : indel_end] contains '-' only
  parts = seq.replace('-', ' ').split()
  assert len(parts) == 2
  left, right = parts[0], parts[1]
  indel_start = seq.index(left) + len(left)
  indel_end = seq.index(right)
  return indel_start, indel_end

def get_mh_len(read, ref, indel_start, indel_end, indel_cat):
  if indel_cat == 'ins':
    # ensure "read" variable has the -
    read, ref = ref, read

  mh_len = 0
  tot_bases_to_try = indel_start
  for idx in range(1, tot_bases_to_try):
    right_side = ref[indel_end - idx]
    left_side = read[indel_start - idx]

    if right_side != left_side:
      break
    else:
      mh_len += 1

  return mh_len

def valid_combination_indel(read, ref):
  num_matches = 0
  for idx, (c1, c2) in enumerate(zip(read, ref)):
    if c1 == c2 and c1 != '-':
      num_matches += 1
  match_frac = num_matches / len(ref)
  return bool(match_frac > 0.7)

##
# Core logic
##
def process_aligns(inp_fn, designed_seq, lib_nm, allfolds_d):
  # expected_cutsite = prefix_len + 28

  prefix_len = len('GATGGGTGCGACGCGTCAT')
  if lib_nm == '12kChar' or '6kV5' in lib_nm:
    start_pos, end_pos = 0, 56
    target_site_len = 56
    zero_pos_idx = 21
  elif lib_nm in ['AtoG', 'CtoT']:
    start_pos, end_pos = 0, 56
    target_site_len = 56
    zero_pos_idx = 10
  elif lib_nm == 'PAMvar':
    start_pos, end_pos = 0, 61
    target_site_len = 61
    zero_pos_idx = 12
  elif lib_nm == 'LibA':
    prefix_len = len('TCCGTGCTGTAACGAAAGGATGGGTGCGACGCGTCAT')
    start_pos, end_pos = 0, 55
    target_site_len = 55
    zero_pos_idx = 9

  allowed_indel_start = prefix_len + zero_pos_idx - 6
  allowed_indel_end = prefix_len + zero_pos_idx + 26 + 1

  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        count = parse_header(line.strip())
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        ref = line.strip()
      if i % 4 == 3:

        '''
          - Require sufficient support on both sides
            (> N nucleotides with mean Q > q)
          - Require a single indel in allowed window
            (wide window means careful shifting in c6 isn't necessary)
          - Downstream, will subtract control frequencies
        '''

        # print('\n'.join([read, ref, line]))

        num_dels, num_ins = count_indels(read, ref)
        num_indels = num_dels + num_ins
        if num_indels > 1:
          if valid_combination_indel(read, ref):
            # allfolds_d[('combination indel')] += count
            print('valid')
          else:
            print('invalid')
          # import code; code.interact(local=dict(globals(), **locals()))
          # print('>1 indel')
          continue
        if num_indels == 0:
          allfolds_d[('wildtype')] += count
          # print('wildtype')
          continue

        # Has exactly one indel
        if num_dels == 1:
          indel_start, indel_end = get_single_indel_range(read)
        elif num_ins == 1:
          indel_start, indel_end = get_single_indel_range(ref)

        if indel_end <= indel_start:
          # Corner case, might happen if right side of indel is very short
          continue

        # Check if indel is in allowed indel range
        if allowed_indel_start <= indel_start <= allowed_indel_end or allowed_indel_start <= indel_end <= allowed_indel_end:
          pass
        else:
          # print('indel not in range')
          continue

        # Check proper alignment support
        support_len = 6
        five_side_read = read[indel_start - support_len : indel_start]
        five_side_ref = ref[indel_start - support_len : indel_start]
        if len(five_side_read) != support_len:
          continue
        if '-' in five_side_read:
          # print('- in five side')
          continue
        if '-' in five_side_ref:
          continue
        match_pct = sum([bool(nt1 == nt2) for nt1, nt2 in zip(five_side_read, five_side_ref)]) / len(five_side_read)
        if match_pct <= 0.90:
          continue

        three_side_read = read[indel_end : indel_end + support_len]
        three_side_ref = ref[indel_end : indel_end + support_len]
        if len(three_side_read) != support_len:
          continue
        if '-' in three_side_read:
          # print('- in three side')
          continue
        if '-' in three_side_ref:
          continue
        match_pct = sum([bool(nt1 == nt2) for nt1, nt2 in zip(three_side_read, three_side_ref)]) / len(three_side_read)
        if match_pct <= 0.90:
          continue

        qs = [ord(s)-33 for s in line.strip()]
        if np.median(qs[indel_start - support_len : indel_start]) <= 30:
          # print('badq in five side')
          continue
        if np.median(qs[indel_end : indel_end + support_len]) <= 30:
          # print('badq in three side')
          continue

        '''
          Columns: 
            Indel category (wildtype, ins, del)
            Genotype position (uniformly defined to be merge friendly)
              seqalign puts indels as far to the right as possible
            Indel length
            MH length
            Inserted bases

          Handled separately/later:
            Target site name
            Count
        '''
        # print('indel')
        indel_len = indel_end - indel_start

        if num_dels == 1:
          cat = 'del'
          ins_bases = ''
        elif num_ins == 1:
          cat = 'ins'
          ins_bases = read[indel_start : indel_end]

        # gt_pos is relative to gRNA zero index
        gtpos_start = indel_start - prefix_len - zero_pos_idx
        gtpos_end = indel_end - prefix_len - zero_pos_idx

        mh_len = get_mh_len(read, ref, indel_start, indel_end, cat)


        indel_detail_tuple = (cat, gtpos_start, gtpos_end, indel_len, mh_len, ins_bases)
        # print(read)
        # print(ref, '\n', indel_detail_tuple)
        allfolds_d[indel_detail_tuple] += count

  return allfolds_d

##
# Iterator
##
def get_indel_count(nm, inp_place, lib_design, lib_nm, start_idx, end_idx):
  data = dict()

  folds = os.listdir(inp_place)

  rows = lib_design.iloc[start_idx : end_idx + 1].iterrows()
  timer = util.Timer(total = end_idx + 1 - start_idx)
  for idx, row in rows:
    timer.update()
    exp = row['Name (unique)']

    if lib_nm != 'PAMvar':
      designed_seq = row['Sequence context (56nt)']
    else:
      designed_seq = row['Sequence context (61nt)']

    allfolds_d = defaultdict(lambda: 0)
    for fold in folds:
      inp_fn = inp_place + '%s/%s.txt' % (fold, exp)
      
      if not os.path.exists(inp_fn):
        continue

      allfolds_d = process_aligns(inp_fn, designed_seq, lib_nm, allfolds_d)

    # Convert allfolds_d to df
    dd = defaultdict(list)
    for key in allfolds_d:
      if key == ('wildtype'):
        continue
      count = allfolds_d[key]
      (cat, gtpos_start, gtpos_end, indel_len, mh_len, ins_bases) = key

      dd['Category'].append(cat)
      dd['Indel start'].append(gtpos_start)
      dd['Indel end'].append(gtpos_end)
      dd['Indel length'].append(indel_len)
      dd['MH length'].append(mh_len)
      dd['Inserted bases'].append(ins_bases)
      dd['Count'].append(count)

    cols_for_wt = [
      # 'Category', 
      'Indel start', 
      'Indel end', 
      'Indel length', 
      'MH length',
      'Inserted bases',
      # 'Count',
    ]
    for col in cols_for_wt:
      dd[col].append(np.nan)
    dd['Category'].append('wildtype')
    dd['Count'].append(allfolds_d[('wildtype')])

    df = pd.DataFrame(dd)
    df['Name'] = exp
    data[exp] = df

  with open(out_dir + '%s_%s_%s.pkl' % (nm, start_idx, end_idx), 'wb') as f:
    pickle.dump(data, f)

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
  for bc in exp_design['Name']:

    for start_idx in range(0, 12000, 2000):
      command = 'python %s.py %s %s %s' % (NAME, bc, start_idx, start_idx + 1999)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, bc, start_idx)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -V -l h_rt=2:00:00,h_vmem=1G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return


@util.time_dec
def main(nm = '', start_idx = '', end_idx = ''):
  print(NAME)

  # Function calls
  if '12kChar' in nm or '6kV5' in nm:
    lib_nm = '12kChar'
  elif 'AtoG' in nm:
    lib_nm = 'AtoG'
  elif 'CtoT' in nm:
    lib_nm = 'CtoT'
  elif 'PAMvar' in nm:
    lib_nm = 'PAMvar'

  inp_dir = _config.OUT_PLACE + 'c_alignment_%s/%s/' % (lib_nm, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm))
  get_indel_count(nm, inp_dir, lib_design, lib_nm, int(start_idx), int(end_idx))

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1], start_idx = sys.argv[2], end_idx = sys.argv[3])
  else:
    gen_qsubs()