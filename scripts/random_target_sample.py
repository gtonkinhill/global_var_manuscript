from mungo.fasta import FastaReader
import sys, os
import random

SAMPLE_SIZE=1000
seqs = {}
seq_list = []
count=0
input=sys.argv[1]

with open(input + "_mappingJump.txt", 'w') as outfile:
  for h,s in FastaReader(input):
    count+=1
    outfile.write("seq" + str(count) + "\t" + h + "\n")
    seqs["seq" + str(count)] = s
    seq_list.append("seq" + str(count))

rand_seq_list = seq_list[:]
random.shuffle(rand_seq_list)

assert seq_list[0]!=rand_seq_list[0], "Lists match!"

with open(input[:-6] + "_randomSampleTargetSize" + str(SAMPLE_SIZE) +".fasta", 'w') as outfile:
  count=0
  for s in rand_seq_list:
    if count<SAMPLE_SIZE:
      outfile.write(">"+"target_"+s+"\n"+seqs[s]+"\n")
    else:
      outfile.write(">"+"db_"+s+"\n"+seqs[s]+"\n")
    count+=1
