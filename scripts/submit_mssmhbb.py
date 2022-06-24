#!/usr/bin/env python

from argparse import ArgumentParser
from argparse import HelpFormatter

arg_desc = '''
Prepare, submit and check jobs to NAF HTCondor batch system.
The configuration file must have a particular name:
<name>_<mode>_<year>.cfg
mode : fh for fullhadronic, sl for semileptonic
year : year of data
'''

parser = ArgumentParser(prog='submit_btag.py', formatter_class=lambda prog: HelpFormatter(prog,indent_increment=6,max_help_position=80,width=280), description=arg_desc,add_help=True)
parser.add_argument('--exe'       , dest='exe'            , required= True          , help='Executable')
parser.add_argument('--cfg'       , dest='cfg'            , required= True          , help='Configuration file')
parser.add_argument('--label'     , dest='label'          , required= True          , help='Unique label to identify the submission')
parser.add_argument('--samples'   , dest='samples'        , required= True          , help='File containing samples to be processed')
parser.add_argument('--btagweight', dest='btagweight'     , action='store_true'     , help='apply btag weight')
parser.add_argument('--submit'    , dest='submit'         , action='store_true'     , help='submit jobs')
parser.add_argument('--status'    , dest='status'         , action='store_true'     , help='status of the jobs')
parser.add_argument('--resubmit'  , dest='resubmit'       , action='store_true'     , help='resubmit failed jobs')
parser.add_argument('--hadd'      , dest='hadd'           , action='store_true'     , help='hadd root files')


args = parser.parse_args()

import os, sys
import getpass
import shutil

from Analysis.MssmHbb.utils import data_info

if not os.path.isdir('condor'):
   print('The condor directory does not exist!')
   print('Create one on dust and link it here')
   quit()

# Parsed arguments and conditions
## Configuration
cfg = args.cfg
if not os.path.isfile(cfg):
   print('The configuration file {} does not exist!'.format(cfg))
   quit()
## Samples
samples_file = args.samples
if not os.path.isfile(samples_file):
   print('-e-: Samples file {} does not exist!'.format(samples_file))
   quit()
with open(samples_file) as f:
   contents = f.read().splitlines()
contents = [x.replace(' ','') for x in contents]
contents = [x for x in contents if not x.startswith('#') and x != '']
samples_dir = contents[0]
samples = contents[1:]
## btagweight
btw = args.btagweight
## Label
label = args.label
if btw:
   label = '{}_btagweight'.format(label)
## Executable
exe = args.exe
## Actions
n_actions = args.submit + args.status + args.resubmit + args.hadd
if n_actions != 1:
   print('-e-: Do not know what to do! Either submit or status or resubmit or hadd?')
   quit()
command = ''
if args.submit:
   command = 'submit'
if args.hadd:
   command  = 'hadd'
if args.resubmit:
   command = 'resubmit'
if args.status:
   command  = 'status'
   
## Metadata
info  = data_info(cfg)
proc  = info['cfg']
year  = info['year']
mode  = info['mode']

name = '{}_{}_{}'.format(proc,mode,year)
label_dir = 'condor/{}'.format(label)
cmssw_base  = os.getenv('CMSSW_BASE')
files_dir = '{}/src/{}'.format(cmssw_base,samples_dir)

cwd = os.getcwd()
os.chdir(cwd)

### SUBMIT
if command == 'submit':
   if not os.path.isdir(label_dir):
      os.mkdir(label_dir)
   shutil.copy(cfg, label_dir)
   submit_command = 'submit_btag.py --exe {} --cfg {} --label {} --samples {} --submit \n'.format(exe,cfg,label,samples_file)
   if btw:
      submit_command = 'submit_btag.py --exe {} --cfg {} --label {} --samples {} --btagweight --submit \n'.format(exe,cfg,label,samples_file)
   with open('{}/submit_commands.log'.format(label_dir), 'a') as file_command:  # append mode
      file_command.write(submit_command)
   os.chdir(label_dir)
   for sample in samples:
      filelist = '{}/{}_rootFileList.txt'.format(files_dir,sample)
      output = 'histograms_{}_{}.root'.format(sample,name)
      exe_cmd = 'naf_submit.py -e {} -c {} -n {} -o {} -l {} -x 1'.format(exe,cfg,filelist,output,sample)
      if btw:
         exe_cmd = 'naf_submit.py -e {} -c {} -n {} -o {} -l {} -x 1 --opts xx_btagweight'.format(exe,cfg,filelist,output,sample)
      os.system(exe_cmd)
   os.chdir(cwd)
 
### STATUS   
if command == 'status' or command == 'resubmit':
   os.chdir(label_dir)
   for sample in samples:
      condor_dir = 'Condor_{}_{}_{}'.format(exe,name,sample)
      exe_cmd = 'naf_submit.py --{} --dir {}'.format(command,condor_dir)
      os.system(exe_cmd)
   os.chdir(cwd)
   
### HADD   
if command == 'hadd':
   res_dir = '{}/results/{}/{}'.format(cwd,label,name)
   if not os.path.isdir(res_dir):
      os.makedirs(res_dir)

   os.chdir(label_dir)
   for sample in samples:
      condor_dir = 'Condor_{}_{}_{}'.format(exe,name,sample)
      os.chdir(condor_dir)
      exe_cmd = 'hadd.csh finished_jobs ; cp -p histograms_*.root {}'.format(res_dir)
      os.system(exe_cmd)
      os.chdir('..')
   os.chdir(cwd)
   
## condor queue status
user = getpass.getuser()   
print('')   
print("--------------")
print("OWNER BATCH_NAME      SUBMITTED   DONE   RUN    IDLE   HOLD  TOTAL JOB_IDS")
exe_cmd = 'condor_q  | grep "{} ID"'.format(user)
os.system(exe_cmd)

quit()
