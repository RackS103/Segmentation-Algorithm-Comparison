"""
AUTHOR:         RAC MUKKAMALA, February 2021
PAPER:          A comparison of two compositional segmetation algorithms for genomic sequences

FUNCTION:       The SimulateSeqCD class generates simulated DNA sequences with domains of homogeneous GC Content.
                Simulated Sequences can either have equal-length domains or variable-length domains.
                
                GENERATION OF SEQUENCES WITH EQUAL-LENGTH DOMAINS - generate_equal():
                Creates a simulated 1 mb simulated DNA sequence evenly divided into the inputted number of domains.
                Each domain is randomly assigned a domain GC content sampled from Bernardi's isochore family
                distribution, with the program ensuring that adjacent isochores are assigned to different families.

                GENERATION OF SEQUENCES WITH VARIABLE-LENGTH DOMAINS - generate_variable():
                Creates a simulated 5 mb simulated DNA sequence divided into variable length domains.
                Domain length is sampled from a power-law distribution with alpha=2.55 and xmin = 10000bp, and
                domain lengths are normalized to create a total sequence length of 5,000,000 bp.
                Each domain is randomly assigned a domain GC content sampled from Bernardi's isochore family
                distribution, with the program ensuring that adjacent isochores are assigned to different families.

USER INPUT:     Number of domains, Sequence ID name
OUTPUT:         A FASTA file with given ID name containing the simulated sequence.
                An JSON file containing a summary of the simulated domains and sequence
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import GC
import random
import math
import json
import os

class SimulateSeqCD:
    output_dir = ""
    json_dir = ""
    def __init__(self, output_dir):
        self.output_dir = output_dir
        if not self.output_dir[len(self.output_dir) - 1] == "/":
            self.output_dir += "/"
        try:
            self.json_dir = self.output_dir + "ground_truth/"
            os.mkdir(self.json_dir)
        except OSError as err:
            print(err)

    #randomly picks nucleotides to create one domain
    def __generate_domain(self, len, GC):
        AT = 100-GC
        #probabilities of A,T,G,C respectively
        p_nucs = [(AT/2), (AT/2), (GC/2), (GC/2)]
        nucs = ['A', 'T', 'G', 'C']
        sequence = []
        for n in range(0, len):
            add = random.choices(population=nucs, weights=p_nucs, k=1)
            sequence.append(add[0])
        return Seq(''.join(sequence))

    #randomly chooses the next isochore family
    #if prev_family is None, that means this is the first time.
    def __choose_family(self, prev_family):
        families = {"L1":0.228, "L2":0.332, "H1":0.227, "H2":0.112, "H3":0.0301}
        normalized_freq_families = list(families.values())
        domain_fam = ""
        domain_gc = 0
        #picks a random family for the first time
        if prev_family == None:
            tot = sum(normalized_freq_families)
            for i in range(0, len(normalized_freq_families)):
                normalized_freq_families[i] /= tot
            domain_fam = random.choices(population=list(families.keys()), weights=normalized_freq_families, k=1)[0]
        
        #picks a random family, making sure not to repeat the previous family
        else:
            normalized_freq_families.remove(families[prev_family])
            options = list(families.keys())
            options.remove(prev_family)
            tot = sum(normalized_freq_families)
            for i in range(0, len(normalized_freq_families)):
                normalized_freq_families[i] /= tot
            domain_fam = random.choices(population=options, weights=normalized_freq_families, k=1)[0]
        
        #picks random GC value within range of chosen family
        if domain_fam == 'L1':
            domain_gc = random.uniform(30,37)
        elif domain_fam == 'L2':
            domain_gc = random.uniform(37,42)
        elif domain_fam == 'H1':
            domain_gc = random.uniform(42,47)
        elif domain_fam == 'H2':
            domain_gc = random.uniform(47,52)
        elif domain_fam == 'H3':
            domain_gc = random.uniform(52,60)

        return (domain_fam, round(domain_gc, 2))

    #randomly samples from the negative power-law distribution P(x) = Cx^-a 
    #input needed: xmin, alpha (positive)
    def __rpower(self, alpha, min):
        r = random.random()
        invcdf_sample = min*math.pow(1-r, 1/(1-alpha))
        return invcdf_sample

    #creates the variable domain lengths using the power law sampler
    def __generate_lengths(self, num_domains, total_length):
        #Generate power distribution values
        power_vals = []
        for i in range(0, int(num_domains)):
            power_vals.append(self.__rpower(2.55,10000))
        pwr_tot = sum(power_vals)
        print("Raw Length (before scaling): " + str(pwr_tot))
        for i in range(0,len(power_vals)):
            power_vals[i] *= ((total_length)/pwr_tot)
            power_vals[i] = round(power_vals[i])
            if power_vals[i] < 10000:
                print("ERROR!")
        return power_vals
    

    #brings together all the methods and creates an EQUAL DOMAIN LENGTH simulated sequence
    def generate_equal(self, num_domains, id):
        seq_length = 1000000
        print(id + ": EQUAL LENGTH " + str(seq_length)) 
        CD_length = int(seq_length/num_domains)
        log = {'domains':num_domains, 'domain_length':CD_length, 'tot_length':seq_length, 'type': 'bernardi-equal'}

        #generate sequence string
        sequence = Seq("")
        last_family = None
        for i in range(0, num_domains):
            domain_fam = self.__choose_family(last_family)
            domain_i = self.__generate_domain(CD_length, domain_fam[1])
            sequence += domain_i
            observed_gc = round(GC(domain_i), 5)
            log.update({"mean_" + str(i+1):observed_gc})
            log.update({"family_" + str(i+1):domain_fam[0]})
            print('Domain ' + str(i+1) + '--> Expected=' + domain_fam[0]  + ":" + str(domain_fam[1]) + ', Observed=' + str(observed_gc))
            last_family = domain_fam[0]
        
        #print to file
        desc = "Simulated Sequence [" + str(num_domains) + " " + str(CD_length) + "bp equal domains, " + str(seq_length) + "bp total length, Bernardi's Isochore Families]"
        filepath = self.output_dir + id
        json_filepath = self.json_dir + id
        record = SeqRecord(sequence, id=("loc|"+id), description=desc)
        SeqIO.write(record, filepath+".fasta", "fasta")
        print("OVERALL OBSERVED MEAN = " + str(round(GC(sequence), 5)))
        with open(json_filepath+'.json', 'w', encoding='UTF-8') as write_json:
            json.dump(log, write_json)

    #brings together all the methods and creates a VARIABLE DOMAIN LENGTH simulated sequence
    def generate_variable(self, num_domains, id):
        seq_length = 5000000
        print(id + " VARIABLE LENGTH " + str(seq_length))
        CD_lengths = self.__generate_lengths(num_domains, seq_length)
        log = {'domains':num_domains, 'tot_length':seq_length, 'domain_length': -1, 'type': 'bernardi variable'}

        #generate sequence string
        sequence = Seq("")
        last_family = None
        for i in range(0, num_domains):
            domain_fam = self.__choose_family(last_family)
            domain_i = self.__generate_domain(CD_lengths[i], domain_fam[1])
            sequence += domain_i
            observed_gc = round(GC(domain_i), 5)
            log.update({"mean_" + str(i+1):observed_gc})
            log.update({"family_" + str(i+1):domain_fam[0]})
            log.update({"length_" + str(i+1):CD_lengths[i]})
            print('Domain ' + str(i+1) + '--> Length=' + str(CD_lengths[i]) + ', Expected=' + domain_fam[0]  + ":" + str(domain_fam[1]) + ', Observed=' + str(observed_gc))
            last_family = domain_fam[0]
        
        #print to file
        desc = "Simulated Sequence [" + str(num_domains) + " variable-length domains, " + str(seq_length) + "bp total length, Bernardi's Isochore Families]"
        filepath = self.output_dir + id
        json_filepath = self.json_dir + id
        record = SeqRecord(sequence, id=("loc|"+id), description=desc)
        SeqIO.write(record, filepath+".fasta", "fasta")
        print("OVERALL OBSERVED MEAN = " + str(round(GC(sequence), 5)))
        with open(json_filepath+'.json', 'w', encoding='UTF-8') as write_json:
            json.dump(log, write_json)

#main, runs SimulateSeqCD class using user input from console
if __name__ == '__main__':
    num_domains = int(input("Number of domains: "))
    id = input("Sequence ID: ")
    type = input("Enter 'E' for equal-length domains, 'V' for variable-length domains ")
    dir = input("Enter output file directory: ")
    print()
    sim = SimulateSeqCD(dir)
    if type == "E":
        sim.generate_equal(num_domains, id)
    else:
        sim.generate_variable(num_domains, id)
    
    