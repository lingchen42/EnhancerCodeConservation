#!/usr/bin/env python

import os
import sys
import datetime

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt

import scipy.integrate

SCORE_FILES = sys.argv[2:]

#OUT_NAME = sys.argv[1].split(".")[0]+"_"+datetime.datetime.now().strftime('%Y-%m-%d')
OUT_NAME = sys.argv[1]+"_"+datetime.datetime.now().strftime('%Y-%m-%d')

SCORE_COL = 3  # what column (0-indexed) in SCORE_FILE holds the score
STATUS_COL = 1  # what column in SCORE_FILE holds the status (+/-)
POS_STR = "1"  # how is the positive state represented in SCORE_FILE
NEG_STR = "-1"  # how is the negative state represented in SCORE_FILE

WRITE_CURVES = False

now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
PREAMBLE = '# %s %s\n#\n' % (' '.join(sys.argv), now_str)
PREAMBLE += '#\n'


def plot_curves(curves, outfile_name='test_curves.pdf', text=None, labels=[],
                title='', plot_type='ROC', cis_low=None, cis_high=None):
    """ Plot one or more ROC/PR curves. """

    #manually set colors
    color_dict = {"Hsap":"#006600", "Mmul":"#5A9BD4","Mmus":"#F15A60","Btau":"#FFDC5B","Cfam":"#FF7500","Mdom":"#C09953"}

    fig = plt.figure()
    fig.set_size_inches(8, 8)
    ax1 = fig.add_subplot(111)

    plt.title(title)

    plt.grid(color='grey', linestyle=':')

    if cis_low and cis_high:
        for i, curve in enumerate(curves):
            plt.fill_between(curve[0], cis_low[i], cis_high[i], alpha=0.3, color='gray')

    plots = []
    # for i, curve in enumerate(curves):
    #     p, = plt.plot(curve[0], curve[1], lw=2)
    #     plots.append(p)
    #     #if plot_type == 'ROC': print '\n', labels[i], curve

    for i, curve in enumerate(curves):
        s=labels[i].split(" ")[0]

        try:
            c = color_dict[s]
        except:
            print "Cannot recognize %s, will use color #006600"%(s)
            c = "#006600"
        
        p, = plt.plot(curve[0], curve[1], lw=2, color=c)
        plots.append(p)
        

    if plot_type == 'ROC':
        leg_loc = 'lower right'
        plt.xlabel("False Positive Rate")
        plt.ylabel('True Positive Rate')
        plt.plot([0,1], color="grey", ls="--")
        plt.xlim(0,1)
        plt.ylim(0,1)

    elif plot_type == 'PR':
        leg_loc = 'upper right'
        #leg_loc = 'lower right'
        plt.xlabel("Recall")
        plt.ylabel('Precision')
        plt.xlim(0,1)
        plt.ylim(0,1)

    if labels != []: 
        leg = plt.legend(plots, labels, loc=leg_loc)
        leg.draw_frame(False)

   
    # handles, labels = ax1.get_legend_handles_labels()
    # print handles
    # print labels
    
    #update font size
    matplotlib.rcParams.update({'font.size':22})

    plt.savefig(outfile_name)


################################################################################

name2stats = {}
for SCORE_FILE in SCORE_FILES:
    total_pos = 0.
    total_neg = 0.

    curve_name = SCORE_FILE.split('/')[-1].replace('.bed', '').replace('.txt', '').replace('.tab', '')
    #plot the species name
    if curve_name.find("applyto")!= -1:
        curve_name = curve_name.split("_applyto_")[1].split("_")[0]
    elif curve_name.find("_cv_scores_")!=-1:
        curve_name = curve_name.split("_cv_scores_")[0].split("_")[-1]
    else:
        print "Enter the legend for %s :"% (curve_name) 

    score_status = []
    for line in open(SCORE_FILE):
        if line.startswith('#') or line.startswith('"names"'): continue

        t = line[:-1].split('\t')
        if len(t) < max(SCORE_COL, STATUS_COL): continue

        score = float(t[SCORE_COL])
        raw_status = t[STATUS_COL]
        status = None

        if raw_status == POS_STR:
            total_pos += 1
            status = 1
        elif raw_status == NEG_STR:
            total_neg += 1
            status = 0
        else:
            sys.exit("ERROR: Couldn't assign +/- status to: %s" % line[:-1])

        score_status.append([score, status])

    total_N = total_pos + total_neg    

    sorted_score_status = sorted(score_status, reverse=True)


    num_tp = 0.; num_fp = 0.
    tprs = []
    fprs = []
    precs = []
    for ss in sorted_score_status:
        score = ss[0]
        status = ss[1]

        if status:
            num_tp += 1
        else:
            num_fp += 1

        tpr = num_tp / total_pos
        fpr = num_fp / total_neg
        prec = num_tp / (num_tp + num_fp)
        #print '\t'.join(['%.3f' % fpr, '%.3f' % tpr, '%.3f' % prec, score])
        tprs.append(tpr)
        fprs.append(fpr)
        precs.append(prec)

    
    roc_auc = scipy.integrate.trapz(tprs, x=fprs)
    pr_auc = scipy.integrate.trapz(precs, x=tprs)

    name2stats[curve_name] = (fprs, tprs, precs, roc_auc, pr_auc)

    print "\n", SCORE_FILE
    print 'ROC AUC: %.3f\nPR AUC: %.3f' % (roc_auc, pr_auc)

labels = []
roc_curve_data = []

#ordered by phylogeny
name_order = ["Hsap","Mmul","Mmus","Btau","Cfam","Mdom"]
for name in name_order:
    try:
        roc_curve_data.append((name2stats[name][0], name2stats[name][1]))   
        labels.append(name+" %.3f"%(name2stats[name][3]))
    except KeyError:
        continue
   
# for name in name2stats:
#     labels.append(name)
#     roc_curve_data.append((name2stats[name][0], name2stats[name][1]))
plot_curves(roc_curve_data, outfile_name=OUT_NAME + "_roc.pdf", 
            plot_type="ROC", labels=labels)

labels=[]
pr_curve_data = []
#ordered by phylogeny
for name in name_order:
    try:
        pr_curve_data.append((name2stats[name][1], name2stats[name][2]))   
        labels.append(name+" %.3f"%(name2stats[name][4]))
    except KeyError:
        continue   
# for name in name2stats:
#     pr_curve_data.append((name2stats[name][1], name2stats[name][2]))

plot_curves(pr_curve_data, outfile_name=OUT_NAME + "_pr.pdf", 
            plot_type="PR", labels=labels)
print "\nWriting PDFs: %s_roc/pr.pdf" % OUT_NAME


if WRITE_CURVES:
    of = open(NAME + '.roc_curve', 'w')
    of.write(PREAMBLE)
    of.write('%s (%.3f)\n' % (NAME, roc_auc))
    of.write(' '.join([str(x) for x in fprs]) + '\n')
    of.write(' '.join([str(x) for x in tprs]) + '\n')
    of.close()

    of = open(NAME + '.pr_curve', 'w')
    of.write(PREAMBLE)
    of.write('%s (%.3f)\n' % (NAME, pr_auc))
    of.write(' '.join([str(x) for x in tprs]) + '\n')
    of.write(' '.join([str(x) for x in precs]) + '\n')
    of.close()




