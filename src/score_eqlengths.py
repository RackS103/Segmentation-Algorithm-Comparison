"""
"""
import os
import json
import csv
import numpy

out = open("EqDom_length_correl.csv", 'w', encoding='UTF-8')
out.write("TRUE_LENGTH,ISOSEGMENTER_AVG,ISOPLOTTER_AVG\n")
docpath = "C:/Users/RackS/Documents/"

isosegmenter_dict = {}
isoplotter_dict = {}
lengths_used = []

def readDict(file):
    predict_data = {}
    with open(file, 'r', encoding='UTF-8') as readfile:
        index = 1
        for line in readfile:
            if(len(line.strip()) > 0):
                args = line.strip().split(" ")
                predict_data.update({str(index): {'Start':int(args[0]), 'End':int(args[1]), 'Size':int(args[2])}})
                index += 1
    return predict_data


for file in os.listdir(docpath + "isoSegmenter100"):
    if file.endswith(".csv") and "E" in file:
        predict_data = csv.DictReader(open(docpath+"isoSegmenter100/"+file, 'r', encoding='UTF-8'))
        seqid = file.replace(".csv", "")
        with open(docpath + "ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)
        
        true_length = truth_data['domain_length']
        if not true_length in lengths_used:
            lengths_used.append(true_length)

        isoseg_preds = []
        for pred_domain in predict_data:
            isoseg_preds.append(int(pred_domain['Size']))

        if not true_length in isosegmenter_dict:
            isosegmenter_dict.update({true_length: isoseg_preds})
        else:
            for value in isoseg_preds:
                isosegmenter_dict[true_length].append(value)

for file in os.listdir(docpath+"isoPlotter100/3.IsoPlotter_no_ns"):
    if 'E' in file:
        predict_data = readDict(docpath+"isoPlotter100/3.IsoPlotter_no_ns/"+file)
        seqid = file.replace(".txt", "")
        with open(docpath+"ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)

        true_length = truth_data['domain_length']
        if not true_length in lengths_used:
            lengths_used.append(true_length)

        isoplt_preds = []
        for iterate in predict_data:
            pred_domain = predict_data[iterate]
            isoplt_preds.append(int(pred_domain['Size']))

        if not true_length in isoplotter_dict:
            isoplotter_dict.update({true_length: isoplt_preds})
        else:
            for value in isoplt_preds:
                isoplotter_dict[true_length].append(value)

for length in lengths_used:
    out.write(str(length) + "," + str(round(numpy.mean(isosegmenter_dict[length]), 5)) + "," + str(round(numpy.mean(isoplotter_dict[length]), 5)) + "\n")

out.close()