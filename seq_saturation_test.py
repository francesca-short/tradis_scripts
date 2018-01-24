#!/usr/bin/env python
import pysam
import argparse
import random
import csv
import matplotlib
import matplotlib.pyplot as plt

# Prepare input arguments
parser = argparse.ArgumentParser(description='''This script measures the sequencing saturation of a 
	TraDIS run by sampling reads from a bam file and returning the insertion site counts and different 
	sampled read numbers to determine if the library has been sequenced to near-saturation.''')
parser.add_argument('bam', help="Input bam file")
parser.add_argument('step', help="Sampling interval", type = int)
parser.add_argument('-o', '--output', help = 'output file prefix')
args = parser.parse_args()
if args.output == None:
    args.output = args.bam

#Read in bam file (note: needs matching index file in same folder)
filepath = args.bam
samfile = pysam.Samfile(args.bam, "rb")

#Extract individual read entries from file and shuffle
readslist = list(samfile.fetch())
print "Number of reads in bam file: " + str(len(readslist))
nreads = len(readslist)
random.shuffle(readslist)

#Open output table
outputfilename = (args.output + "_seq_saturation.csv")
csvfile = open(outputfilename, 'wb')
writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
header = ("Reads_sampled","Insertion_sites")
writer.writerow(header)

#open results array
read_count_list = []
ins_count_list = []

#New option (quicker), counts a new n reads in shuffled read list each iteration, adds to total ins count
end = args.step
start = 0
nreads = len(readslist)
results = set()
count = 0
while end < nreads:
    for read in readslist[start:end]:
        if read.flag == 0:
            results.add(str(read.reference_start) + "Fwd_" + str(read.reference_id))
        elif read.flag == 16:
            results.add(str(read.reference_end) + "Rev_" + str(read.reference_id))
    row = (end, len(results))
    writer.writerow(row)
    print (str(len(results)) + " with " + str(end) + " reads sampled.")
    read_count_list.append(end)
    ins_count_list.append(len(results))
    start = start + args.step
    end = end + args.step

#special case for remaining reads at end
for read in readslist[(end-args.step):nreads+1]:
    if read.flag == 0:
        # print (read.reference_start)
        results.add(str(read.reference_start) + "Fwd_" + str(read.reference_id))
    elif read.flag == 16:
        # print (read.reference_end)
        results.add(str(read.reference_end) + "Rev_" + str(read.reference_id))
row = (nreads, len(results))
writer.writerow(row)
print (str(len(results)) + " with " + str(nreads) + " reads sampled.")
read_count_list.append(nreads)
ins_count_list.append(len(results))

print(read_count_list)
print(ins_count_list)

#plot results as simple line graph of insertion sites vs reads sampled
plt.plot(read_count_list, ins_count_list)
plt.xlabel('Reads sampled')
plt.ylabel('Insertion sites')
plt.ylim([0,(max(ins_count_list)+50000)])
plt.xlim([0,(nreads+50000)])
plot_name = str(args.output + "plot.png")
plt.savefig(plot_name, dpi=300)