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
parser.add_argument('--beamy',
                    '-y',
                    type=float,
                    nargs='+',
                    default=[0],
                    help='List of y values of beam center')
parser.add_argument('--tile_number',
                    '-n',
                    type=int,
                    default=3,
                    help='Number of tile layers')
parser.add_argument('--is_wrap_ESR',
                    '-i',
                    type=int,
                    default=1,
                    help='1 for use ESR, 0 for use Tyvek')
parser.add_argument('--useProton',
                    '-U',
                    type=int,
                    default=0,
                    help='1 for use proton, 0 for photon')
parser.add_argument('--tilewidth',
                    '-L',
                    type=float,
                    nargs='+',
                    default=[30],
                    help='Tile width')
parser.add_argument('--tilethick',
                    '-t',
                    type=float,
                    nargs='+',
                    default=[3],
                    help='Tile thickness')
parser.add_argument('--beamwidth',
                    '-w',
                    type=float,
                    nargs='+',
                    default=[1.5],
                    help='width of beam [mm]')
parser.add_argument('--dimplerad',
                    '-r',
                    type=float,
                    nargs='+',
                    default=[2.65],
                    help='List of dimple radius to try out')
parser.add_argument('--dimpleind',
                    '-d',
                    type=float,
                    nargs='+',
                    default=[1.6],
                    help='List of dimple indents to try out')
parser.add_argument('--sipmwidth',
                    '-W',
                    type=float,
                    nargs='+',
                    default=[1.3],
                    help='SiPM width paramter')
parser.add_argument('--sipmrim',
                    '-R',
                    type=float,
                    nargs='+',
                    default=[0.4],
                    help='SiPM inactive rim width parameter')
parser.add_argument('--sipmstand',
                    '-S',
                    type=float,
                    nargs='+',
                    default=[0.2],
                    help='SiPM stand height')
parser.add_argument('--gap_pcb_wrap',
                    '-G',
                    type=float,
                    nargs='+',
                    default=[0.2],
                    help='Gap between PCB and wrap')
parser.add_argument('--absmult',
                    '-a',
                    type=float,
                    nargs='+',
                    default=[1.],
                    help='List of tile absorption length multiple')
parser.add_argument('--wrapreflect',
                    '-m',
                    type=float,
                    nargs='+',
                    default=[0.985],
                    help='List of wrap reflectivity')
parser.add_argument('--pcbreflect',
                    '-p',
                    type=float,
                    nargs='+',
                    default=[0.5],
                    help='List of pcb reflectivities to test')
parser.add_argument('--pcbradius',
                    '-b',
                    type=float,
                    nargs='+',
                    default=[3.175],
                    help='List of pcb exposed radii to test')
parser.add_argument('--NEvents',
                    '-N',
                    type=int,
                    default=100,
                    help='Number of events to run')
parser.add_argument('--Material',
                    '-M',
                    type=int,
                    default=1,
                    help='Material, 1EJ, 2SCSN')
parser.add_argument('--PDE',
                    '-E',
                    type=int,
                    default=1,
                    help='PDE, 1 S13, 2 S14')
parser.add_argument('--dimplesa',
                    '-D',
                    type=float,
                    nargs='+',
                    default=[1.3],
                    help='dimple Sigma Alpha')
parser.add_argument('--prefix',
                    type=str,
                    default='',
                    help='Additional string to place in filename')
args = parser.parse_args()

BASE_DIR = os.path.abspath(os.environ['CMSSW_BASE'] + '/src/' +
                           '/HGCalTileSim/condor/')
DATA_DIR = os.path.abspath(BASE_DIR + '/DN_depth_comple')
CONDOR_JDL_TEMPLATE = """
universe              = vanilla
Executable            = {0}/condor-LYSquareTrigger_CMSSW.sh
should_transfer_files = NO
Requirements          = TARGET.FileSystemDomain == "privnet" 
#&& (TARGET.SlotID>5)
request_memory        = 1 GB
Output                = {1}.$(ClusterId).$(ProcId).stdout
Error                 = {1}.$(ClusterId).$(ProcId).stderr
Log                   = {1}.$(ClusterId).condor
Arguments             = {2}
Queue 10
"""
#Requirements = TARGET.Machine =?= "r510-0-4.privnet" || TARGET.Machine =?= "r510-0-5.privnet" || TARGET.Machine =?= "r510-0-6.privnet"|| TARGET.Machine =?= "r510-0-9.privnet"|| TARGET.Machine =?= "r510-0-10.privnet"|| TARGET.Machine =?= "r510-0-11.privnet" || TARGET.Machine =?= "r540-0-20.privnet" || TARGET.Machine =?= "r540-0-21.privnet" || TARGET.Machine =?= "r720-0-2.privnet" || TARGET.Machine =?= "r720-0-1.privnet" || TARGET.Machine =?= "siab-1.umd.edu"

for x, y, L, t, w, r, d, W, R, S, G, a, m, p, b, D in [(x, y, L, t, w, r, d, W, R, S, G, a, m, p, b, D)
                                           for x in args.beamx
                                           for y in args.beamy
                                           for L in args.tilewidth
                                           for t in args.tilethick
                                           for w in args.beamwidth
                                           for r in args.dimplerad
                                           for d in args.dimpleind
                                           for W in args.sipmwidth
                                           for R in args.sipmrim
                                           for S in args.sipmstand
                                           for G in args.gap_pcb_wrap
                                           for a in args.absmult
                                           for m in args.wrapreflect
                                           for p in args.pcbreflect
                                           for b in args.pcbradius
                                           for D in args.dimplesa]:

  def make_str(prefix):
    args_string = '_'.join([
        'x{0:.1f}'.format(x), 'y{0:.1f}'.format(y), 'n{}'.format(args.tile_number), 'esr_{}'.format(args.is_wrap_ESR), 'L{0:.1f}'.format(L), 'thick_{0:.1f}'.format(t),
        'beamwidth_{0:.1f}'.format(w), 'r{0:.1f}'.format(r), 'd{0:.1f}'.format(d), 'sipmW_{0:.1f}'.format(W), 'R{0:.1f}'.format(R), 'S{0:.2f}'.format(S*100),'G{0:.1f}'.format(G*100),
        'a{0:.1f}'.format(a*100), 'm{0:.1f}'.format(m*1000), 'P{:.1f}'.format(p*100),'Pcbr_{:.1f}'.format(b), 'DimpleSA_{:.1f}'.format(D),
         'particle{}_'.format(args.useProton),'PDE{:.1f}_'.format(args.PDE),'Material{:.1f}_'.format(args.Material),
    ])
    return prefix + args.prefix + '_' + args_string.replace('.', 'p')

  save_filename = os.path.abspath(DATA_DIR + '/root/' + '/' +
                                  make_str('hgcal_tilesim') + '.root')

  condor_args = ' '.join([
      '-x {}'.format(x), '-y {}'.format(y), '-n {}'.format(args.tile_number),'-i {}'.format(args.is_wrap_ESR),
      '-L {}'.format(L),'-t {}'.format(t), '-w {}'.format(w),
      '-r {}'.format(r), '-d {}'.format(d), '-W {}'.format(W),
      '-R {}'.format(R), '-S {}'.format(S),'-G {}'.format(G),
      '-a {}'.format(a), '-m {}'.format(m),
      '-p {}'.format(p), '-b {}'.format(b),'-D {}'.format(D), '-N {}'.format(args.NEvents), '-M {}'.format(args.Material), 
      '-E {}'.format(args.PDE),
      '-U {}'.format(args.useProton),
      '-o {}'.format(os.path.abspath(save_filename)),
  ])

  log_filename = os.path.abspath(DATA_DIR + '/condor/' + '/' +
                                 make_str('log_tilesim'))
  jdl_filename = os.path.abspath(DATA_DIR + '/condor/' + '/' +
                                 make_str('hgcal_tilesim') + '.jdl')
  jdl_content = CONDOR_JDL_TEMPLATE.format(BASE_DIR,
                                           log_filename,
                                           condor_args)

  ## Writing jdl files
  os.makedirs(os.path.dirname(jdl_filename), exist_ok=True)
  with open(jdl_filename, 'w') as file:
    file.write(jdl_content)

  os.makedirs(DATA_DIR + '/log/', exist_ok=True)

  ## Making directory for output files
  os.makedirs(os.path.dirname(save_filename), exist_ok=True)

  print(jdl_filename)
