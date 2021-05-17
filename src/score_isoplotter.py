"""
"""
import os
import json

cutoff = float(input("Tolerance (decimal)? "))
docpath = "C:/Users/RackS/Documents/"
out = open("isoplotter_scoring_error"+str(int(cutoff*100))+".csv", 'w', encoding='UTF-8')
summary = open("isoplotter_score_summary_error"+str(int(cutoff*100)), 'w', encoding='UTF-8')
out.write("SEQUENCE_ID,TYPE,DOMAINS,TP,FP,FN,Sens,PPV,Jaccard\n")


def readDict(file):
    predict_data = {}
    with open(file, 'r', encoding='UTF-8') as readfile:
        index = 1
        for line in readfile:
            if(len(line.strip()) > 0):
                args = line.strip().split(" ")
                predict_data.update({str(index): {'Start':int(args[0]), 'End':int(args[1])}})
                index += 1
    return predict_data


#SCORE EQUAL LENGTH SEQUENCES
tp_eq = 0
fp_eq = 0
fn_eq = 0
for file in os.listdir(docpath+"isoPlotter100/3.IsoPlotter_no_ns"):
    if 'E' in file:
        predict_data = readDict(docpath+"isoPlotter100/3.IsoPlotter_no_ns/"+file)
        seqid = file.replace(".txt", "")
        with open(docpath+"ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)

        true_boundaries = []
        tp_seq = 0
        fp_seq = 0
        fn_seq = 0
        for i in range(0, int(truth_data['tot_length']) + 1, int(truth_data['domain_length'])):
            true_boundaries.append(i)

        for iterate in predict_data:
            matched = False
            pred_domain = predict_data[iterate]
            for i in range(0, len(true_boundaries) - 1):
                startdiff = pred_domain['Start'] - true_boundaries[i]
                enddiff = pred_domain['End'] - true_boundaries[i+1]
                tolerance = cutoff*(true_boundaries[i+1] - true_boundaries[i])
                if abs(startdiff) <= tolerance:
                    if abs(enddiff) <= tolerance:
                        tp_seq += 1
                        matched = True
                        print(seqid)
                        print("START MATCH: " + str(true_boundaries[i]) + ", " + str(pred_domain['Start']))
                        print("END MATCH: " + str(true_boundaries[i+1]) + ", " + str(pred_domain['End']))
                        print("DIFFERENCES: " + str(startdiff) + ", " + str(enddiff) + ", TOLERANCE = " + str(tolerance))
                        print()
                        break
            if not matched:
                fp_seq += 1

        fn_seq = int(truth_data['domains']) - tp_seq
        tp_eq += tp_seq
        fp_eq += fp_seq
        fn_eq += fn_seq
        sensitivity = round(tp_seq/(tp_seq + fn_seq), 5)
        ppv = round(tp_seq/(tp_seq+fp_seq), 5)
        jaccard = round(tp_seq/(tp_seq + fp_seq + fn_seq), 5)
        out.write(seqid+",E,"+str(truth_data['domains'])+","+str(tp_seq)+","+str(fp_seq)+","+str(fn_seq)+","+str(sensitivity)+","+str(ppv)+","+str(jaccard)+"\n")

#WRITE EQ SUMMARY STATS TO FILE
summary.write("EQUAL-LENGTH STATISTICS\n")
summary.write("TP equal domain: " + str(tp_eq) + "\n")
summary.write("FP equal domain: " + str(fp_eq) + "\n")
summary.write("FN equal domain: " + str(fn_eq) + "\n")
summary.write("Sensitivity: " + str(round(tp_eq/(tp_eq + fn_eq),5)) + "\n")
summary.write("Precision(PPV): " + str(round(tp_eq/(tp_eq + fp_eq),5)) + "\n")
summary.write("Jaccard Index: " + str(round(tp_eq/(tp_eq + fp_eq + fn_eq),5)) + "\n\n")

#SCORE VARIABLE SEQUENCE PREDICTIONS
tp_var = 0
fp_var = 0
fn_var = 0
for file in os.listdir(docpath+"isoPlotter100/3.IsoPlotter_no_ns"):
    if 'V' in file:
        predict_data = readDict(docpath+"isoPlotter100/3.IsoPlotter_no_ns/"+file)
        seqid = file.replace(".txt", "")
        with open(docpath+"ground_truth100/"+seqid+".json", 'r', encoding='UTF-8') as json_file:
            truth_data = json.load(json_file)

        true_boundaries = [1]
        tp_seq = 0
        fp_seq = 0
        fn_seq = 0
        for i in range(1, int(truth_data['domains']) + 1):
            b_next = true_boundaries[i-1] + int(truth_data['length_'+str(i)])
            true_boundaries.append(b_next)

        for iterate in predict_data:
            matched = False
            pred_domain = predict_data[iterate]
            for i in range(0, len(true_boundaries) - 1):
                startdiff = pred_domain['Start'] - true_boundaries[i]
                enddiff = pred_domain['End'] - true_boundaries[i+1]
                tolerance = cutoff*(true_boundaries[i+1] - true_boundaries[i])
                if abs(startdiff) <= tolerance:
                    if abs(enddiff) <= tolerance:
                        tp_seq += 1
                        matched = True
                        print(seqid)
                        print("START MATCH: " + str(true_boundaries[i]) + ", " + str(pred_domain['Start']))
                        print("END MATCH: " + str(true_boundaries[i+1]) + ", " + str(pred_domain['End']))
                        print("DIFFERENCES: " + str(startdiff) + ", " + str(enddiff) + ", TOLERANCE = " + str(tolerance))
                        print()
                        break
            if not matched:
                fp_seq += 1

        fn_seq = int(truth_data['domains']) - tp_seq
        tp_var += tp_seq
        fp_var += fp_seq
        fn_var += fn_seq
        sensitivity = round(tp_seq/(tp_seq + fn_seq), 5)
        ppv = round(tp_seq/(tp_seq+fp_seq), 5)
        jaccard = round(tp_seq/(tp_seq + fp_seq + fn_seq), 5)
        out.write(seqid+",V,"+str(truth_data['domains'])+","+str(tp_seq)+","+str(fp_seq)+","+str(fn_seq)+","+str(sensitivity)+","+str(ppv)+","+str(jaccard)+"\n")

#WRITE VARLEN SUMMARY STATS TO FILE
summary.write("VARIABLE-LENGTH STATISTICS\n")
summary.write("TP equal domain: " + str(tp_var) + "\n")
summary.write("FP equal domain: " + str(fp_var) + "\n")
summary.write("FN equal domain: " + str(fn_var) + "\n")
summary.write("Sensitivity: " + str(round(tp_var/(tp_var + fn_var),5)) + "\n")
summary.write("Precision(PPV): " + str(round(tp_var/(tp_var + fp_var),5)) + "\n")
summary.write("Jaccard Index: " + str(round(tp_var/(tp_var + fp_var + fn_var),5)) + "\n\n")
    
#WRITE OVERALL ISOPLOTTER STATS TO FILE
summary.write("OVERALL STATISTICS\n")
summary.write("TP: " + str(tp_var + tp_eq) + "\n")
summary.write("FP: " + str(fp_var + fp_eq) + "\n")
summary.write("FN: " + str(fn_var + fn_eq) + "\n")
summary.write("Sensitivity: " + str(round((tp_var + tp_eq)/(tp_var + fn_var + tp_eq + fn_eq),5)) + "\n")
summary.write("Precision(PPV): " + str(round((tp_var + tp_eq)/(tp_var + fp_var + tp_eq + fp_eq),5)) + "\n")
summary.write("Jaccard Index: " + str(round((tp_var + tp_eq)/(tp_var + fp_var + fn_var + tp_eq + fp_eq + fn_eq),5)) + "\n")