"""
"""
import os
import json
import csv
import numpy

out = open("VarDom_length_correl.csv", 'w', encoding='UTF-8')
out.write("TRUE_DOMAINS,ISOSEGMENTER_NUM,ISOSEGMENTER_AVG,ISOPLOTTER_NUM,ISOPLOTTER_AVG\n")
docpath = "C:/Users/RackS/Documents/"

isosegmenter_dict = {}
isoplotter_dict = {}
nums_used = []

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



for file in os.listdir(docpath+"isoSegmenter100"):
    if file.endswith(".csv") and "V" in file:
        predict_data = csv.DictReader(open(docpath+"isoSegmenter100/"+file, 'r', encoding='UTF-8'))
        seqid = file.replace(".csv", "")
        with open(docpath+"ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)
        
        true_num = truth_data['domains']
        if not true_num in nums_used:
            nums_used.append(true_num)

        isoseg_preds = []
        isoseg_num = 0
        for pred_domain in predict_data:
            isoseg_preds.append(int(pred_domain['Size']))
            isoseg_num += 1

        if not true_num in isosegmenter_dict:
            isosegmenter_dict.update({true_num: {'Num':[isoseg_num], 'Lengths': isoseg_preds}})
        else:
            for value in isoseg_preds:
                isosegmenter_dict[true_num]['Lengths'].append(value)
            isosegmenter_dict[true_num]['Num'].append(isoseg_num)

for file in os.listdir(docpath+"isoPlotter100/3.IsoPlotter_no_ns"):
    if 'V' in file:
        predict_data = readDict(docpath+"isoPlotter100/3.IsoPlotter_no_ns/"+file)
        seqid = file.replace(".txt", "")
        with open(docpath+"ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)

        true_num = truth_data['domains']
        if not true_num in nums_used:
            nums_used.append(true_num)

        isoplt_preds = []
        isoplt_num = 0
        for iterate in predict_data:
            pred_domain = predict_data[iterate]
            isoplt_preds.append(int(pred_domain['Size']))
            isoplt_num += 1

        if not true_num in isoplotter_dict:
            isoplotter_dict.update({true_num: {'Num':[isoplt_num], 'Lengths': isoplt_preds}})
        else:
            for value in isoplt_preds:
                isoplotter_dict[true_num]['Lengths'].append(value)
            isoplotter_dict[true_num]['Num'].append(isoplt_num)


for num in nums_used:
    out.write(str(num) + "," + str(round(numpy.mean(isosegmenter_dict[num]['Num']), 5)) +","+ str(round(numpy.mean(isosegmenter_dict[num]['Lengths']), 5)) + "," + str(round(numpy.mean(isoplotter_dict[num]['Num']), 5)) + "," + str(round(numpy.mean(isoplotter_dict[num]['Lengths']), 5)) + "\n")

out.close()
