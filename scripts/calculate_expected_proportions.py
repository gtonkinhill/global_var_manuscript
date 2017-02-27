import sys, os
import argparse
from subprocess import check_call
import numpy as np
from collections import defaultdict
import glob
from joblib import Parallel, delayed

def loadMappingFile(mapfile, verbose):
  mapping_dict = {}
  read_to_short_dict = {}
  if verbose:
    print "loading mapping file...", mapfile
  with open(mapfile, 'rU') as infile:
    for line in infile:
      line=line.strip()
      map_ = line.split(";size")[0]
      map_ = map_.split(";sample")[0]
      map_ = map_.split()[-1]
      mapping_dict[line.split()[0]] = map_
      read_to_short_dict[map_] = line.split()[0]

  return mapping_dict, read_to_short_dict

def splitPosteriors(filenames, temp_dir, verbose):
  #Split post files up to make things more efficient.
  for f in filenames:
    split_dict=defaultdict(list)
    if verbose:
      print "splitting...", f
    with open(f, 'rU') as infile:
      for t in range(5): #Skip header
        next(infile)
      for line in infile:
        target = line.split(",")[0].split("_")[1]
        split_dict[target].append(line)
    for target in split_dict:
      with open(temp_dir + target + "_postSplit.txt", 'w') as outfile:
        for line in split_dict[target]:
          outfile.write(line)

  return

def loadPosteriors(filenames, mapping_dict, isolate, iso_to_reads, verbose):

  seq_lengths = {}
  hmm_posterior = {}
  reads_in_mapping = set()

  if verbose:
    print "Loading posterior files..."
  print filenames
  for f in filenames:
    if verbose:
      print "Loading...", f
    with open(f, 'rU') as infile:
      for line in infile:
        line = line.strip().split(",")
        #Now rename using mapping file
        line[0] = mapping_dict[line[0].split("_")[1]]
        line[1] = mapping_dict[line[1].split("_")[1]]

        #only keep those target seqs that are required for the target isolate
        if line[0] not in iso_to_reads[isolate]: continue

        seq_lengths[line[0]] = int(line[2])
        hmm_posterior[(line[0], line[1])] = [float(t) for t in line[3:]]
        reads_in_mapping.add(line[0])
        reads_in_mapping.add(line[1])
  return seq_lengths, hmm_posterior, reads_in_mapping

def loadOtuMatrix(otu_file, verbose):
  num_isolates_per_read = {}
  iso_to_reads = defaultdict(list)
  with open(otu_file, 'rU') as infile:
    header = infile.next()
    isolates = np.asarray(header.strip().split()[2:])
    for line in infile:
      line=line.strip().split()

      #Check if read is in mapping file
      #that is check it was translated and used in the JHMM algorithm
      # if line[0] not in reads_in_mapping: continue

      index = np.asarray(line[1:], dtype=int)
      num_isolates_per_read[line[0]] = np.sum(index)
      for iso in isolates[index==1]:
        iso_to_reads[iso].append(line[0])

  return num_isolates_per_read, iso_to_reads

def calculateProportion(seq_lengths, hmm_posterior, num_isolates_per_read
  , iso_to_reads, isolate, verbose):
  proportions = {}

  for isoA in iso_to_reads:
    if isolate!="" and isoA==isolate:
      for isoB in iso_to_reads:
        # if isoA==isoB: continue
        chunk_sum = 0
        seq_length = 0
        for rt in iso_to_reads[isoA]:
          seq_length += seq_lengths[rt]
          if num_isolates_per_read[rt]==1:
            for rs in iso_to_reads[isoB]:
              if rs==rt: continue #dont compare to self
              chunk_sum += sum(hmm_posterior[(rt,rs)])*1/float(num_isolates_per_read[rs])
          else:
            for rs in iso_to_reads[isoB]:
              if rs==rt:
                chunk_sum += seq_lengths[rt]*1/float(num_isolates_per_read[rs])
        proportions[(isoA, isoB)] = chunk_sum/float(seq_length)

  return proportions

def writeProportions(proportions, outputfile, verbose):
  if verbose:
    print "Writing proportions to...", outputfile

  with open(outputfile, 'w') as outfile:
    for pair in proportions:
      outfile.write(",".join([pair[0], pair[1]
        , str(proportions[pair])]) + "\n")

  return

def calculateForEachRead(outputfile, isolate
  , temp_dir, mapping_dict, read_to_short_dict
  , num_isolates_per_read, iso_to_reads
  , verbose):
  isolate_pos_files = []
  for read in iso_to_reads[isolate]:
    isolate_pos_files.append(temp_dir + read_to_short_dict[read]+ "_postSplit.txt")

  seq_lengths, hmm_posterior, reads_in_mapping = loadPosteriors(isolate_pos_files
    , mapping_dict, isolate, iso_to_reads, verbose)

  proportions = calculateProportion(seq_lengths, hmm_posterior
    , num_isolates_per_read, iso_to_reads
    , isolate, verbose)

  writeProportions(proportions, outputfile, verbose)

  return

def main():

  parser = argparse.ArgumentParser(description='Calculate proportions from Mosaic and clustering analysis')

  parser.add_argument('--otu', dest='otufile'
    , help="location of binary otu matrix."
    , required=True)

  parser.add_argument('--map', dest='mapfile'
    , help="location of map file, to map the sequence names in mosaic to those in the otu matrix."
    , required=True)

  parser.add_argument('--out_dir', dest='outputdir'
    , help="location of output directory."
    , required=True)

  parser.add_argument('--temp_dir', dest='temp_dir'
    , help="location of temporary directory to store seperated files."
    , required=True)

  parser.add_argument('--pos', nargs='+'
    , dest='pos_files'
    , help='location of posterior files from Mosaic run.')

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  parser.add_argument('--isolate', dest='isolate'
    , default=""
    , help='isolate to infer mixture on.')

  parser.add_argument('--compute_all', dest='compute_all', action='store_true'
    , default=False
    , help='isolate to infer mixture on.')

  parser.add_argument('--cpu', dest='cpu'
    , default=1
    , type=int
    , help='number of cpus.')

  args = parser.parse_args()

  #get full path names
  args.outputdir = os.path.abspath(args.outputdir) + "/"
  args.mapfile = os.path.abspath(args.mapfile)
  args.otufile = os.path.abspath(args.otufile)
  args.temp_dir = os.path.abspath(args.temp_dir) + "/"

  pos_files = []
  for i,f in enumerate(args.pos_files):
    args.pos_files[i] = os.path.abspath(f)

  mapping_dict, read_to_short_dict = loadMappingFile(args.mapfile, args.verbose)

  num_isolates_per_read, iso_to_reads = loadOtuMatrix(args.otufile, args.verbose)

  if len(glob.glob(args.temp_dir + "*_postSplit.txt"))<1:
    splitPosteriors(args.pos_files, args.temp_dir, args.verbose)

  if args.compute_all:
    Parallel(n_jobs=args.cpu)(delayed(calculateForEachRead)(args.outputdir + iso + "_proportions.txt"
      , iso
      , args.temp_dir, mapping_dict, read_to_short_dict
      , num_isolates_per_read, iso_to_reads
      , args.verbose) for iso in iso_to_reads)
  else:
    outputfile=args.outputdir + args.isolate + "_proportions.txt"
    calculateForEachRead(outputfile, args.isolate
      , args.temp_dir, mapping_dict, read_to_short_dict
      , num_isolates_per_read, iso_to_reads
      , args.verbose)

  return


if __name__ == '__main__':
  main()
