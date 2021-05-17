"""
AUTHOR:        RAC MUKKAMALA, January 2021
PAPER:         A comparison of two compositional segmetation algorithms for genomic sequences
FUNCTION:      Shows a graph of window-based GC content for the FASTA sequence inputted.
USER INPUT:    Sequence ID, window size
OUTPUT:        A matplotlib plot of the GC content along the sequences
"""

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import GC
import json
import math

id = input("Sequence FASTA file: ").replace(".fasta", "")
logfile = input("Sequence JSON Metadata file: ")
record = SeqIO.read(id+'.fasta', "fasta")
window_size = int(input("Window Size: "))
labels = input("Labels? ")

with open(logfile, 'r', encoding='UTF-8') as read_json:
    log = json.load(read_json)
num_domains = int(log['domains'])
domain_length = int(log['domain_length'])
variable = False
if domain_length == -1:
    variable = True 
seq_length = int(log['tot_length'])

points = {'x':[], 'y':[]}

for i in range(0, seq_length, window_size):
    GCwindow = GC(record.seq[i:i+window_size])
    points['x'].append(i)
    points['y'].append(GCwindow)

plt.plot(points['x'], points['y'], "c-")
plt.xlabel("Position along Simulated DNA (mb)")
plt.ylabel("Percent GC Content of " + str(window_size) + "bp windows")
plt.title(record.description[4:])

summed_lengths = 0
for i in range(0,num_domains):
    size_i = domain_length
    if variable:
        size_i = int(log['length_' + str(i+1)])

    plt.axvline(x=summed_lengths, c='k')
    domain_gc = float(log['mean_' + str(i+1)])
    plt.hlines(y=domain_gc, xmin=summed_lengths, xmax=summed_lengths+size_i, color='k', linewidth=4, zorder=15)
    
    if labels == "true":
        sd = math.sqrt(0.4*0.6/window_size) * 100
        text_height = domain_gc + 3*sd
        if text_height > 0.99*plt.ylim()[1]:
            text_height -= sd
        plt.text(x=summed_lengths, y=text_height, s=str(round(domain_gc,1))+'%', fontsize=9, fontweight='bold', c='k', zorder=20)
        try:
            if 'bernardi' in log['type']:
                plt.text(x=summed_lengths, y=domain_gc+0.1*sd, s=log['family_'+str(i+1)], fontsize=12, fontweight='bold', c='k', zorder=20)
        except Exception:
            pass
    
    summed_lengths += size_i

plt.axvline(x=int(log['tot_length']), c='k')
plt.show()

