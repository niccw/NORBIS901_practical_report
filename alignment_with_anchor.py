#!/usr/bin/env python3

import sys
import re
import argparse
import numpy as np
from collections import OrderedDict


# reference: https://github.com/aeyc/allAlign/tree/d6279b02a0f256412bb5ffb33f6b854e0dae9547
# reference: https://github.com/risky998/NW-Algos/blob/master/nwaffine.py

class Seq:
    def __init__(self,seq_id:str,seq:str):
        self.seq_id = seq_id
        self.seq = seq
        self.idx = None

    def __repr__(self):
        return {
            "seq_id": self.seq_id,
            "seq": self.seq
        }
    def __str__(self):
        return f">{self.seq_id}\n{self.seq}"

    def set_seq_index(self,start_idx:int,end_idx:int):
        self.idx = list(range(start_idx, end_idx+1))

class AlignParam:
    def __init__(self,gapopen:int,gapextend:int,match:int,mismatch:int):
        self.gapopen = gapopen
        self.gapextend = gapextend
        self.match = match
        self.mismatch = mismatch
    def __repr__(self):
        return {
            "gapopen":self.gapopen,
            "gapextend":self.gapextend,
            "match":self.match,
            "mismatch":self.mismatch
        }

class AlignmentInfo:
    def __init__(self,s1:str,s2:str,alignment_type:str,alignment_param:AlignParam):
        self.s1 =s1
        self.s2 = s2
        self.alignment_type = alignment_type
        self.alignment_param = alignment_param

class SeqPair:
    def __init__(self,s1:Seq,s2:Seq):
        self.seq_list = []
        self.s1 = s1
        self.s2 = s2

    def add_seqs(self,s1_pos:list,s2_pos:list,alignment_type:str,alignment_param:AlignParam):
        sub_seq1 = self.s1.seq[s1_pos[0]:s1_pos[1]]
        sub_seq2 = self.s2.seq[s2_pos[0]:s2_pos[1]]
        self.seq_list.append(AlignmentInfo(s1=sub_seq1, s2=sub_seq2,alignment_type=alignment_type,alignment_param=alignment_param))

class BuildMatrix:
    def __init__(self,ap:AlignParam,s1,s2):
        self.ap = ap
        self.s1_len = len(s1)
        self.s2_len = len(s2)
    
    def upper_matrix(self):
        m = np.zeros((self.s1_len+1,self.s2_len+1))
        for i in range(1, self.s2_len + 1):
            m[0][i] =  self.ap.gapopen + i * self.ap.gapextend
        for j in range(1, self.s1_len + 1):
            m[j][0] = float('-inf')
        return m

    def lower_matrix(self):
        m = np.zeros((self.s1_len+1,self.s2_len+1))
        for i in range(1, self.s2_len + 1):
            m[0][i] = float('-inf')
        for j in range(1, self.s1_len + 1):
            m[j][0] =  self.ap.gapopen + i* self.ap.gapextend
        return m

    def middle_matrix(self):
        m = np.zeros((self.s1_len+1,self.s2_len+1))
        for i in range(1, self.s2_len + 1):
            m[0][i] = float('-inf')
        for j in range(1, self.s1_len + 1):
            m[j][0] =  float('-inf')
        return m

    def traceback_matrix(self):
        m = np.zeros((self.s1_len+1,self.s2_len+1))
        return m

class Anchor:
    def __init__(self,a_id):
        self.a_id = a_id
        self.s1_id = None
        self.s1_start = None
        self.s1_end = None
        self.s2_id = None
        self.s2_start = None
        self.s2_end = None

    def set_anchor(self,seq_id,s_start:int,s_end:int):
        if self.s1_id is not None and self.s2_id is not None:
            return "This anchor has been set"
        if self.s1_id is None:
            self.s1_id = seq_id
            self.s1_start = s_start
            self.s1_end = s_end
        else:
            self.s2_id = seq_id
            self.s2_start  = s_start
            self.s2_end = s_end


    def size(self):
        return [self.s1_end-self.s1.start, self.s2_end-self.s2.start]

class AnchorList:
    def __init__(self, anchor:dict):
        self.all_anchor = anchor
        self.anchor_list = self.share_anchor()
    
    # FIXME:assume achor are in order and not overlapping! Better to add a safe check here
    def share_anchor(self): #a[seq_id][a_id] = [start, end]
        d = {}
        anchor_id_d = {}
        # find shared a_id
        for seq_id,anchor_info in self.all_anchor.items():
            if seq_id not in anchor_id_d:
                anchor_id_d[seq_id] = []
            for a_id, _ in anchor_info.items():
                anchor_id_d[seq_id].append(a_id)
        seq_id_ls = list(self.all_anchor.keys())
        shared_id = []
        for a_id in anchor_id_d[seq_id_ls[0]]:
            if a_id in anchor_id_d[seq_id_ls[1]]:
                shared_id.append(a_id)
        # build share anchor dict
        d = OrderedDict() # make sure this dict always return same order
        for seq_id, anchor_info in self.all_anchor.items():
            for a_id, a_pos in anchor_info.items():
                if a_id in shared_id:
                    if a_id not in d:
                        d[a_id] = Anchor(a_id=a_id)
                    d[a_id].set_anchor(seq_id = seq_id, s_start = a_pos[0], s_end = a_pos[1])
        # TODO: sort d
        self.anchor_list = d

    # define how to align each substring in the seq
    def anchor2subseq(self,s1:Seq,s2:Seq):
        seq_pair = SeqPair(s1=s1,s2=s2)
        a_len = len(self.anchor_list.keys())
        if a_len == 1:
            # 0 -> a
            a = list(self.anchor_list.values())[0]
            seq_pair.add_seqs(alignment_type="local_fix_end",alignment_param = general_ap ,s1_pos=[0,a.s1_start],s2_pos=[0,a.s2_start])
            # a
            seq_pair.add_seqs(alignment_type="global",alignment_param = anchor_ap ,s1_pos=[a.s1_start, a.s1_end],s2_pos=[a.s2_start, a.s2_end])
            # a -> end
            seq_pair.add_seqs(alignment_type="local_fix_start",alignment_param = general_ap ,s1_pos=[a.s1_end,len(s1.seq)],s2_pos=[a.s2_end,len(s2.seq)])
        else:
            a_id_ls = list(self.anchor_list.keys())
            for i, a_id in enumerate(self.anchor_list.keys()): # anchor[a_id] = [seq_id,start,end]
                a = self.anchor_list[a_id]
                if i == 0:
                    seq_pair.add_seqs(alignment_type="local_fix_end",alignment_param = general_ap ,s1_pos=[0,a.s1_start],s2_pos=[0,a.s2_start])
                    seq_pair.add_seqs(alignment_type="global",alignment_param = anchor_ap ,s1_pos=[a.s1_start, a.s1_end],s2_pos=[a.s2_start, a.s2_end])
                    # how to access next a here?
                    next_a_id = a_id_ls[a_id_ls.index(a_id)+1]
                    next_a = self.anchor_list[next_a_id]
                    seq_pair.add_seqs(alignment_type="global",alignment_param = general_ap ,s1_pos=[a.s1_end, next_a.s1_start],s2_pos=[a.s2_end, next_a.s2_start])

                if i == a_len-1:
                    seq_pair.add_seqs(alignment_type="global",alignment_param = anchor_ap ,s1_pos=[a.s1_start, a.s1_end],s2_pos=[a.s2_start, a.s2_end])
                    seq_pair.add_seqs(alignment_type="local_fix_start",alignment_param = general_ap ,s1_pos=[a.s1_end,len(s1.seq)],s2_pos=[a.s2_end,len(s2.seq)])
                else:
                    seq_pair.add_seqs(alignment_type="global",alignment_param = anchor_ap ,s1_pos=[a.s1_start, a.s1_end],s2_pos=[a.s2_start, a.s2_end])
                    # how to access next a here?
                    next_a_id = a_id_ls[a_id_ls.index(a_id)+1]
                    next_a = self.anchor_list[next_a_id]
                    seq_pair.add_seqs(alignment_type="global",alignment_param = general_ap ,s1_pos=[a.s1_end, next_a.s1_start],s2_pos=[a.s2_end, next_a.s2_start])
        return seq_pair

# parse fasta
def parse_fasta(fasta):
    seq_ls = []
    seq_id = ""
    read_seq = False
    with open(fasta, 'r') as f:
        for line in f:
            if read_seq:
                seq_ls.append(Seq(seq_id=seq_id,seq=line.strip()))
                read_seq = False
            if line.startswith(">"):
                seq_id = line.strip().replace(">","")
                read_seq = True
    return seq_ls

# trace back function
def global_traceback(t,s1,s2,i=None,j=None):
    """
    t: traceback matrix
    x: upper matrix
    s2: lower matris1

    if i and j is not provided, assume full glbal traceback
    else, trace from (i,j) to the begining of seq
    """
    if i is None and j is None:
        i = len(s1)
        j = len(s2)

    s1_align = ""
    s2_align = ""

    while i > 0 and j > 0:
        tb_code = t[i][j]
        if tb_code == 1:
            s1_align += s1[i-1]
            s2_align += s2[j-1]
            i -= 1
            j -= 1
        elif tb_code ==2: # gap in s1_align
            s1_align += "-"
            s2_align += s2[j-1]
            j -= 1
        elif tb_code == 3: # gap in s2_align
            s1_align += s1[i-1]
            s2_align += "-"
            i -= 1
    
    # if either one of the seq not finish
    if i!=0 or j!=0:
        while i > 0:
            s1_align += s1[i-1]
            s2_align += "-"
            i -= 1
        while j > 0:
            s1_align += "-"
            s2_align += s2[j-1]
            j -=1
    
    # flip seq
    s1_align = s1_align[::-1]
    s2_align = s2_align[::-1]
    return(s1_align,s2_align)

def local_traceback(t,m,s1,s2):
    """
    t: traceback matrix
    x: upper matrix
    s2: lower matris1

    This is a 'half' local traceback function. Always started from the end of the alignment but allow 'stop' aligning when the score hit 0
    """
    i = len(s1)
    j = len(s2)
    current_score = m[i][j]

    s1_align = ""
    s2_align = ""

    while i > 0 and j > 0 and m[i][j] > 0:
        
        tb_code = t[i][j]
        if tb_code == 1:
            s1_align += s1[i-1]
            s2_align += s2[j-1]
            i -= 1
            j -= 1
        elif tb_code ==2: # gap in s1_align
            s1_align += "-"
            s2_align += s2[j-1]
            j -= 1
        elif tb_code == 3: # gap in s2_align
            s1_align += s1[i-1]
            s2_align += "-"
            i -= 1
      
    # flip seq
    s1_align = s1_align[::-1]
    s2_align = s2_align[::-1]
    return(s1_align,s2_align)

# global alignment (for A and B in between 2 As)
def nw_global(a_info:AlignmentInfo):
    s1 = a_info.s1
    s2 = a_info.s2
    ap = a_info.alignment_param

    row = len(s1)+1
    col = len(s2)+1

    matrix_builder = BuildMatrix(ap,s1,s2)
    x = matrix_builder.upper_matrix() # i_x
    y = matrix_builder.lower_matrix() # i_y
    m = matrix_builder.middle_matrix() # m
    t = matrix_builder.traceback_matrix() # t

    for i in range(1,row):
        for j in range(1,col):
            # middle matrix
            s = ap.match if s1[i-1] == s2[j-1] else ap.mismatch
            m[i][j] = s + max(m[i-1][j-1], x[i-1][j-1], y[i-1][j-1])

            # gap (upper and lower) matrix
            # allow jump from x to y
            x[i][j] = max((ap.gapopen + m[i][j-1]),(ap.gapopen + y[i][j-1]),(ap.gapextend + x[i][j-1]))
            y[i][j] = max((ap.gapopen + m[i-1][j]),(ap.gapopen + x[i-1][j]),(ap.gapextend + y[i-1][j]))

            max_t = max(m[i][j], x[i][j], y[i][j])
            

            # fill traceback matrix
            if max_t == m[i][j]:
                t[i][j] = 1
            elif max_t == x[i][j]:
                t[i][j] = 2
            elif max_t == y[i][j]:
                t[i][j] = 3
                 
    # find max_score (in which matrix?)
    max_score = max(m[i][j], x[i][j], y[i][j])
    alignment = global_traceback(t,s1,s2)

    return alignment

# global alignment (for A and B in between 2 As)
def nw_local_fix_start(a_info:AlignmentInfo):
    s1 = a_info.s1
    s2 = a_info.s2
    ap = a_info.alignment_param
    row = len(s1)+1
    col = len(s2)+1

    best_score = 0
    best_i = 0
    best_j = 0

    matrix_builder = BuildMatrix(ap,s1,s2)
    x = matrix_builder.upper_matrix() # i_x
    y = matrix_builder.lower_matrix() # i_y
    m = matrix_builder.middle_matrix() # m
    t = matrix_builder.traceback_matrix() # t

    for i in range(1,row):
        for j in range(1,col):
            # middle matrix
            s = ap.match if s1[i-1] == s2[j-1] else ap.mismatch
            m[i][j] = s + max(m[i-1][j-1], x[i-1][j-1], y[i-1][j-1])

            # gap (upper and lower) matrix
            # allow jump from x to y
            x[i][j] = max((ap.gapopen + m[i][j-1]),(ap.gapopen + y[i][j-1]),(ap.gapextend + x[i][j-1]))
            y[i][j] = max((ap.gapopen + m[i-1][j]),(ap.gapopen + x[i-1][j]),(ap.gapextend + y[i-1][j]))

            max_t = max(m[i][j], x[i][j], y[i][j])
            if max_t > best_score:
                best_score = max_t
                best_i = i
                best_j = j

            # fill traceback matrix
            if max_t == m[i][j]:
                t[i][j] = 1
            elif max_t == x[i][j]:
                t[i][j] = 2
            elif max_t == y[i][j]:
                t[i][j] = 3
                
    
    # find max_score (in which matrix?)
    alignment = global_traceback(t,s1,s2,best_i,best_j)
    return alignment

def nw_local_fix_end(a_info:AlignmentInfo):
    s1 = a_info.s1
    s2 = a_info.s2
    ap = a_info.alignment_param
    row = len(s1)+1
    col = len(s2)+1

    best_score = 0
    best_i = 0
    best_j = 0

    matrix_builder = BuildMatrix(ap,s1,s2)
    x = matrix_builder.upper_matrix() # i_x
    y = matrix_builder.lower_matrix() # i_y
    m = matrix_builder.middle_matrix() # m
    t = matrix_builder.traceback_matrix() # t

    for i in range(1,row):
        for j in range(1,col):
            # middle matrix
            s = ap.match if s1[i-1] == s2[j-1] else ap.mismatch
            m[i][j] = s + max(m[i-1][j-1], x[i-1][j-1], y[i-1][j-1])

            # gap (upper and lower) matrix
            # allow jump from x to y
            x[i][j] = max((ap.gapopen + m[i][j-1]),(ap.gapopen + y[i][j-1]),(ap.gapextend + x[i][j-1]))
            y[i][j] = max((ap.gapopen + m[i-1][j]),(ap.gapopen + x[i-1][j]),(ap.gapextend + y[i-1][j]))

            max_t = max(0,m[i][j], x[i][j], y[i][j])

            # fill traceback matrix
            if max_t == m[i][j]:
                t[i][j] = 1
            elif max_t == x[i][j]:
                t[i][j] = 2
            elif max_t == y[i][j]:
                t[i][j] = 3
    
    alignment = local_traceback(t,m,s1,s2)
    
    return alignment

def read_anchor(filename):
    a = {} # a[seq_id][a_id] = [start, end]
    with open(filename,"r") as f:
        for line in f:
            seq_id,a_start,a_end,a_id = line.strip().split()
            if seq_id not in a:
                a[seq_id] = {}
            a[seq_id][a_id] = [int(a_start), int(a_end)]
    aa = AnchorList(a)
    aa.share_anchor()
    return  aa

# high level align function
def align_all():
    s1, s2 = parse_fasta(args.fa)
    anchors = read_anchor(args.anchor)
    anchors.share_anchor()
    seq_pair = anchors.anchor2subseq(s1=s1,s2=s2)
    # align one by one
    s1_align = ""
    s2_align = ""
    for a_info in seq_pair.seq_list:
        if a_info.alignment_type == "global":
            s1_tmp,s2_tmp = nw_global(a_info)
        elif a_info.alignment_type == "local_fix_start":
            s1_tmp,s2_tmp = nw_local_fix_start(a_info)
        elif a_info.alignment_type == "local_fix_end":
            s1_tmp,s2_tmp = nw_local_fix_end(a_info)

        s1_align += s1_tmp
        s2_align += s2_tmp
    print(f">{s1.seq_id}\n{s1_align}")
    print(f">{s2.seq_id}\n{s2_align}")


if __name__ == "__main__":
    # chop seq into A,B (by precalculated index)
    # align each region using AlignParam
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str, dest="fa",required=True, help="fasta file contains 2 sequences")
    parser.add_argument("-a","--anchor",type=str,dest="anchor",required=True,help="anchor file (BED format)")
  
    parser.add_argument("-am","--anchor_match",type=int,dest="a_match",default=2,help="default = 2")
    parser.add_argument("-ap","--anchor_mismatch",type=int,dest="a_mismatch",default=-3,help="default = -3")
    parser.add_argument("-ago","--anchor_gapopen",type=int,dest="a_gapopen",default=-4,help="default = -4")
    parser.add_argument("-age","--anchor_gapextend",type=int,dest="a_gapextend",default=-3,help="default = -3")

    parser.add_argument("-m","--match",type=int,dest="match",default=2,help="default = 2")
    parser.add_argument("-p","--mismatch",type=int,dest="mismatch",default=-3,help="default = -3")
    parser.add_argument("-go","--gapopen",type=int,dest="gapopen",default=-3,help="default = -3")
    parser.add_argument("-ge","--gapextend",type=int,dest="gapextend",default=-1,help="default = -1")

    args = parser.parse_args()
    
    general_ap = AlignParam(gapopen=args.gapopen,gapextend=args.gapextend,match=args.match,mismatch=args.mismatch)
    anchor_ap = AlignParam(gapopen=args.a_gapopen,gapextend=args.a_gapextend,match=args.a_match,mismatch=args.a_mismatch)

    align_all()

    