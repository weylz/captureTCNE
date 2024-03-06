#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@ File    : filterMatrix.py
@ Time    : 2023/01/01 12:00:00
@ Author  : Wenyong Zhu
@ Version : 1.0.0
@ Desc    : filter eligible elements from transcriptional signal matrix of sample group.
'''


from sys import argv

num_sample = int(argv[1])
num_cutoff = int(argv[2])

# num_row = []
with open('transcriptional.signal.cutoff.value.txt') as cf:
    cf = float(cf.readlines()[0].strip())

with open(f'TCNE.over{num_cutoff}.matrix', 'w') as rlt:
    with open('transcriptional_signal.matrix') as f:
        f = f.readlines()
        # count_transcript = 0
        rlt.writelines(f[0])
        for i in range(1, len(f)):
            f[i] = f[i].strip().split('\t')
            if 'Y' not in f[i][0]:
                count_nozero = 0
                for j in range(4, num_sample + 4):           
                    if f[i][j] == ".": 
                        break
                    elif float(f[i][j]) > cf:
                        count_nozero += 1
                if count_nozero >= num_cutoff:
                    rlt.writelines('\t'.join(f[i][:4] + f[i][5:]) + '\n')
                    # count_transcript += 1
                    # num_row.append(i)
    # print(count_transcript, num_row)
