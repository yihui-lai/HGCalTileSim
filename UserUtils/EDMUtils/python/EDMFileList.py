##
#  @file    EDMFileList.py
#  @brief   Utility functions for listing EDM Files
#  @author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
#
#*******************************************************************************
import subprocess, shlex

def ListEDMDatasets(primary_tag,process_tag,tier):
    grep_cmd = "grep /.*/.*/.*"
    query_cmd = '/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/das_client/v03.01.00/bin/das_client.py --limit=1000000 --query="/{}/{}/{}"'

    if not primary_tag:
        primary_tag = '*',
    if not process_tag:
        process_tag = '*'
    if not tier:
        tier = '*'
    query_cmd = query_cmd.format( primary_tag, process_tag, tier )

    query_proc = subprocess.Popen(
            shlex.split(query_cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
    )
    grep_proc = subprocess.Popen(
            shlex.split(grep_cmd),
            stdin=query_proc.stdout,
            stdout=subprocess.PIPE
    )
    grep_result, dummy  = grep_proc.communicate()
    if len(grep_result.split())==0:
        print "Warning! No files found from query ", query_cmd
    return grep_result.split()

def ListFileInDataset( dataset ):
    query_cmd = 'das_client --limit=10000000 --query="file dataset={}"'.format( dataset )
    grep_cmd  = 'grep "/store"'
    query_proc   = subprocess.Popen( shlex.split(query_cmd), stdout=subprocess.PIPE )
    grep_proc    = subprocess.Popen( shlex.split(grep_cmd), stdin=query_proc.stdout, stdout=subprocess.PIPE )
    grep_result, dummy  = grep_proc.communicate()

    return grep_result.split()
