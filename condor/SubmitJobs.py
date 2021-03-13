#!/usr/bin/env python3
import os
import datetime
import argparse

parser = argparse.ArgumentParser(
    description='Options for generating condor scripts')
parser.add_argument('--beamx',
                    '-x',
                    type=float,
                    nargs='+',
                    required=True,
                    help='List of x values of beam center')
parser.add_argument('--beamz',
                    '-z',
                    type=float,
                    nargs='+',
                    default=[0],
                    help='List of z values of beam center')
parser.add_argument('--tilewidth',
                    '-l',
                    type=float,
                    nargs='+',
                    default=[30],
                    help='Tile width')
parser.add_argument('--beamidth',
                    '-w',
                    type=float,
                    nargs='+',
                    default=[1.5],
                    help='beam width parameter')
parser.add_argument('--absmult',
                    '-a',
                    type=float,
                    nargs='+',
                    default=[1000.],
                    help='List of tile absorption length')
parser.add_argument('--wrapreflect',
                    '-m',
                    type=float,
                    nargs='+',
                    default=[0.985],
                    help='List of wrap reflectivity')
parser.add_argument('--NEvents',
                    '-N',
                    type=int,
                    default=100,
                    help='Number of events to run')
parser.add_argument('--useProton',
                    '-P',
                    type=int,
                    default=1,
                    help='useProton')
parser.add_argument('--handwrap',
                    '-H',
                    type=int,
                    default=1,
                    help='handwrap')
parser.add_argument('--prefix',
                    type=str,
                    default='',
                    help='Additional string to place in filename')

args = parser.parse_args()

BASE_DIR = os.path.abspath(os.environ['CMSSW_BASE'] + '/src/' +
                           '/HGCalTileSim/condor/')
DATA_DIR = os.path.abspath(BASE_DIR + '/results/')

CONDOR_JDL_TEMPLATE = """
universe              = vanilla
Executable            = {0}/condor-LYSquareTrigger_CMSSW.sh
should_transfer_files = NO
Requirements          = TARGET.FileSystemDomain == "privnet"
request_memory        = 1 GB
Output                = {1}.stdout
Error                 = {1}.stderr
Log                   = {1}.condor
Arguments             = {2}
Queue
"""

for x, z, l, w, a, m, in [
    (x, z, l, w, a, m,) for x in args.beamx
    for z in args.beamz for l in args.tilewidth 
    for w in args.beamidth for a in args.absmult
    for m in args.wrapreflect 
]:

  def make_str(prefix):
    args_string = '_'.join([
        'beamx{0:.1f}'.format(x), 'beamz{0:.1f}'.format(z), 'beamw{0:.1f}'.format(w),
        'Tilel{0:.1f}'.format(l), 'abs{0:.1f}'.format(a), 'handwrap{0:.0f}'.format(args.handwrap), 
        'wrapref{0:.1f}'.format(m),
        'useP{0:.0f}'.format(args.useProton), 
    ])
    return prefix + args.prefix + '_' + args_string.replace('.', 'p')

  save_filename = os.path.abspath(DATA_DIR + '/root/' + '/' +
                                  make_str('hgcal_tilesim') + '.root')

  condor_args = ' '.join([
      '-x {}'.format(x), '-z {}'.format(z), '-l {}'.format(l), '-w {}'.format(w),
      '-a {}'.format(a), '-m {}'.format(m), '-P {}'.format(args.useProton), '-H {}'.format(args.handwrap), 
      '-N {}'.format(args.NEvents), '-o {}'.format(
          os.path.abspath(save_filename)),
  ])

  log_filename = os.path.abspath(DATA_DIR + '/log/' + '/' +
                                 make_str('log_tilesim'))
  jdl_filename = os.path.abspath(DATA_DIR + '/condor/' + '/' +
                                 make_str('hgcal_tilesim') + '.jdl')
  jdl_content = CONDOR_JDL_TEMPLATE.format(BASE_DIR, log_filename, condor_args)

  ## Writing jdl files
  os.makedirs(os.path.dirname(jdl_filename), exist_ok=True)
  with open(jdl_filename, 'w') as file:
    file.write(jdl_content)

  os.makedirs(DATA_DIR + '/log/', exist_ok=True)

  ## Making directory for output files
  os.makedirs(os.path.dirname(save_filename), exist_ok=True)

  print(jdl_filename)
