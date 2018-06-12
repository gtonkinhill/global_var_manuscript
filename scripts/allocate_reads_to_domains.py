from mungo.fasta import FastaReader
from mungo.sequence import sixFrameTranslation
from collections import defaultdict
from subprocess import check_call
import argparse
import os

HMMER_SEARCH = "hmmsearch"

def translateInto6Frames(fasta_file, out_fasta, verbose):
  with open(out_fasta, 'w') as outfile:
      for h,s in FastaReader(fasta_file):
          for frame, seq in sixFrameTranslation(s).items():
              outfile.write(">"+h+" _frame_"+str(frame)+"\n")
              outfile.write(seq+"\n")
  return out_fasta

def runHmmer(fasta_file, hmmer_file, evalue, out_file, verbose):
  hmm_cmd = (HMMER_SEARCH
      + " --tblout " + out_file
      + " -E " + str(evalue)
      + " --nonull2 --nobias"
      + " " + hmmer_file
      + " " + fasta_file
      + " > /dev/null")

  print hmm_cmd
  check_call(hmm_cmd, shell=True)

  return out_file

def getTopHits(hmmer_output, outputfile, verbose):
  top_hits = defaultdict(lambda: ("", 100.0))

  #Get top hits for each sequence
  for line in open(hmmer_output, 'rU'):
    if line[0]=="#": continue
    line=line.strip().split()
    read = line[0]
    E = float(line[4])
    if  E < top_hits[read][1]:
      top_hits[read] = (line[2], E)

  with open(outputfile, 'w') as outfile:
    outfile.write("read,domain,e-value\n")
    for read in top_hits:
      outfile.write(read+","+top_hits[read][0]+","+str(top_hits[read][1])+"\n")


def main():

  parser = argparse.ArgumentParser(description='Allocate reads to domains using HMMER.')

  parser.add_argument('--fasta', dest='fasta_file'
    , help="location of fasta file with DNA DBLa tags."
    , required=True)

  parser.add_argument('--hmm', dest='hmm_file'
    , help="location of HMMER model file."
    , required=True)

  parser.add_argument('--evalue', dest='evalue', type=float
    , help="minimum e-value for HMMER match (default=0.01)"
    , default=0.01)

  parser.add_argument('--out_dir', dest='outputdir'
    , help="location of output directory."
    , required=True)

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , help="turn on verbose output (default=False)"
    , default=False)

  args = parser.parse_args()

  #get full path names
  args.outputdir = os.path.abspath(args.outputdir) + "/"
  args.fasta_file = os.path.abspath(args.fasta_file)
  args.hmm_file = os.path.abspath(args.hmm_file)

  prefix = os.path.splitext(os.path.basename(args.fasta_file))[0]

  protein_file = translateInto6Frames(args.fasta_file
    , args.outputdir + prefix + "_Protein.fasta" , args.verbose)

  hmmer_out = runHmmer(protein_file, args.hmm_file, args.evalue
    , args.outputdir + prefix + "hmmerOut.txt" , args.verbose)

  getTopHits(hmmer_out, args.outputdir + prefix + "_domainAllocations.csv"
    , args.verbose)

  return


if __name__ == '__main__':
  main()