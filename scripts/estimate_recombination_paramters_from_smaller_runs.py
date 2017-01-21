import sys, os
import argparse
from subprocess import check_call
import numpy as np
from mungo.fasta import FastaReader
from joblib import Parallel, delayed

MOSAIC = "/home/users/allstaff/tonkin-hill.g/mosaic_mod/mosaic"


def sampleSeqs(fastafile, num_targets, db_size, outputfile, run, verbose):

  if verbose:
    print "Sampling ", db_size, " seqs with ", num_targets, "targets"

  np.random.seed(run)

  #read in fastafile
  seqs = []
  for h,s in FastaReader(fastafile):
    seqs.append((h,s))

  index = range(len(seqs))
  index = np.random.choice(index, db_size, replace=False)

  with open(outputfile, 'w') as outfile:
    count = 0
    for i in index:
      if count<num_targets:
        outfile.write(">target_"+seqs[i][0]+"\n"+seqs[i][1]+"\n")
      else:
        outfile.write(">db_"+seqs[i][0]+"\n"+seqs[i][1]+"\n")
      count+=1

  return outputfile

def runMosaic(recomb, read, prefix, verbose):

  sub_prefix = prefix +"_recomb" +str(recomb)

  mosaic_cmd = (MOSAIC
    + " -aa -ma"
    + " -seq " + sub_prefix+".fasta"
    + " -tag " + sub_prefix + "_mosaic"
    + " -rec " + str(recomb)
    + " -estimate -group 2 db target -target target"
    + " > " + "mosaic_output_run"+str(run)+".txt")

  if verbose:
    print "running... ", mosaic_cmd

  check_call(mosaic_cmd, shell=True)

  return


def main():
  parser = argparse.ArgumentParser(description='Estimate hmm paramters by sub sampling.')

  parser.add_argument('-o', '--outputDir'#, action=writeable_dir
    , dest='outputdir'
    , help="location of output directory. Will be created if it doesn't exist"
    , required=True)

  parser.add_argument('-r', '--read', dest='read'#, action=readable_file
    , help="location of first read file."
    , required=True)

  parser.add_argument('--db_size', dest='db_size', type=int, default=1000
    , help="number of sequences to use in the search database. (default=1000)")

  parser.add_argument('--num_targets', dest='num_targets', type=int, default=10
    , help="number of targets used in paramter inference. (default=10)")

  parser.add_argument('--iters', dest='iters', type=int, default=100
    , help="number iterations of sub sampling to perform. (default=100)")

  parser.add_argument('--cpu', dest='cpu', type=int, default=1
    , help="number of cpus to use. (default=1)")

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  args.read = os.path.abspath(args.read)
  args.outputdir = os.path.abspath(args.outputdir)

  cwd = os.getcwd()
  os.chdir(args.outputdir)

  prefix = os.path.splitext(os.path.basename(args.read))[0]

  for i in range(1,args.iters+1):
    sample_read = sampleSeqs(args.read, args.num_targets, args.db_size
      , prefix+ "run" + str(i) + "_.fasta", i, args.verbose)
    Parallel(n_jobs=args.cpu)(
    delayed(runMosaic)(
      r, sample_read, prefix+ "run" + str(i), args.verbose) for r in np.arange(0.0001,0.0101,0.0001))


  os.chdir(cwd)

  return




if __name__ == '__main__':
  main()





