from mungo.fasta import FastaReader
import sys, os

seqs = {}
seq_list = []
count=0
input=sys.argv[1]

with open(input + "_mapping.txt", 'w') as outfile:
  for h,s in FastaReader(input):
    count+=1
    outfile.write("seq" + str(count) + "\t" + h + "\n")
    seqs["seq" + str(count)] = s
    seq_list.append("seq" + str(count))


n=30

run=0
for i in range(0, len(seq_list), n):
  targets = set(seq_list[i:i+n])
  run+=1
  with open(input[:-6] + "_run" + str(run) + ".fasta", 'w') as outfile:
    for t in targets:
      outfile.write(">"+"target_"+t+"\n"+seqs[t]+"\n")
    for s in seq_list:
      if s in targets: continue
      outfile.write(">"+"db_"+s+"\n"+seqs[s]+"\n")
