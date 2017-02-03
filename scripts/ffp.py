from mungo.fasta import FastaReader
import numpy as np
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy
import argparse


def loadSequences(inputfile, kmer_length, verbose):
  #Load sequences and convert to RY format. Also set up kmer index.
  #RY format: R stands for bases A and G and the Y stands for bases C and T

  if verbose:
    print "Loading sequence file: ", inputfile

  isolate_seqs = defaultdict(list)
  kmers = {}
  kmer_count = 0
  for h,s in FastaReader(inputfile):
    iso = h.split(".")[0]
    s = s.replace("A", "R")
    s = s.replace("G", "R")
    s = s.replace("C", "Y")
    s = s.replace("T", "Y")
    isolate_seqs[iso].append(s)
    for i in range(0,len(s)-kmer_length+1):
      kmer = s[i:i+kmer_length]
      if kmer not in kmers:
        kmers[kmer] = kmer_count
        kmer_count += 1

  return isolate_seqs, kmers

def generateFreqMatrix(isolate_seqs, kmers, kmer_length, verbose):
  #Generate count matrix of size N_k x N_samples
  N_k = len(kmers.keys())
  isolates = isolate_seqs.keys()
  N_samples = len(isolates)

  if verbose:
    print "Generating count matrices..."
    print len(kmers.keys()), " unqiue kmers found."

  freqMatrix = np.zeros((N_k, N_samples), dtype=float)

  for j,iso in enumerate(isolates):
    if verbose:
      print "processing isolate: ", iso
    for seq in isolate_seqs[iso]:
      for i in range(0,len(seq)-kmer_length+1):
        freqMatrix[kmers[seq[i:i+kmer_length]],j] += 1

  #Now normalise by the total count in each sample
  freqMatrix = freqMatrix/np.sum(freqMatrix, axis=0)

  return freqMatrix, isolates

def JSD(P, Q):
  #Calculate Jensen Shannon divergence between two distributions
  M = 0.5 * (P + Q)
  return 0.5 * (entropy(P, M) + entropy(Q, M))

def calculateDistMatrix(freqMatrix, isolates, outputfile, verbose):
  #Now we want to generate a distance matrix using Jensen Shannon divergence
  if verbose:
    print "Calculating distance matrix..."

  dm = squareform(pdist(np.transpose(freqMatrix), JSD))

  #Now print out the resulting matrix in PHYLIP format
  with open(outputfile, 'w') as outfile:
    outfile.write(str(len(isolates))+"\n")
    for i,iso in enumerate(isolates):
      outfile.write(" ".join([iso]+[str(d) for d in dm[i,:]]) + "\n")

  return



def main():

  parser = argparse.ArgumentParser(description='Calculate transition probabilities from Viterbi algorithm')

  parser.add_argument('--kmer_length', dest='kmer_length'
    , help="length of kmers to use."
    , type=int
    , required=True)

  parser.add_argument('--out', dest='outputfile'
    , help="location of output file."
    , required=True)

  parser.add_argument('--seq', dest='seq_file'
    , help=("location of dna sequences with isolates as first string before '.'"
      + " in the header"))

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  isolate_seqs, kmers = loadSequences(args.seq_file, args.kmer_length
    , args.verbose)

  freqMatrix, isolates = generateFreqMatrix(isolate_seqs, kmers
    , args.kmer_length, args.verbose)

  calculateDistMatrix(freqMatrix, isolates, args.outputfile, args.verbose)

  return


if __name__ == '__main__':
  main()














