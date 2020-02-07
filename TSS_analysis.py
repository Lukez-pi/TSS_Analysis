# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:16:57 2020

@author: lukez
"""
from Bio import SeqIO
from matplotlib import pyplot as plt
import openpyxl
import pandas as pd

sequence_file = "Putida KT2440.gbff"#"putida_v3.gb" #"Putida KT2440.gbff" #  "sequence.gb"
def excel_file_processing():
    f = pd.read_excel("primary TSS.xlsx", usecols=lambda x: 'Unnamed' not in x)
    return f
    
# Condition, Position, Strand, Locus tag, Product, Gene length, UTR length, Primary

excel_f = excel_file_processing()

list_of_genes = []
seq = dict()
for seq_record in SeqIO.parse(sequence_file, "genbank"):
    gene_counter = 1
    previous_start = 0
    for feature in seq_record.features:
        if feature.type == "gene":
            locus_tag = feature.qualifiers["locus_tag"]
            strand = feature.strand
            if strand == -1:
                start_position = feature.location.end
                end_position = feature.location.start
            else:
                start_position = feature.location.start
                end_position = feature.location.end

            if feature.location.start >= feature.location.end:
                raise Exception("start and end position doesn't make sense")
            if start_position <= previous_start:
                print(previous_start)
                print(start_position)
                print(locus_tag)
                print("overlapping genes")
            previous_start = start_position
            
            seq[locus_tag[0]] = {"feature": feature, "gene_num": gene_counter}
            list_of_genes.append(locus_tag[0])
            gene_counter += 1
            
    print("total number of genes")
    print(gene_counter)

def test_gene(f):
    # f is the excel file
    unfound_tags = []
    unmatch_gene_len = []
    unmatch_start_pos = []
    unmatch_start_neg = []
    unmatch_direction = []
    def apply_fcn(row):
        tag = row["Locus tag"][:2] + "_" + row["Locus tag"][2:]
        tss_pos = row["Position"]
        excel_gene_len = row["Gene length"]
        utr_len = row["UTR length"]
        if tag not in seq:
            unfound_tags.append(row["Locus tag"])
        else:
            start_pos = seq[tag]["feature"].location.start
            end_pos = seq[tag]["feature"].location.end
            strand = seq[tag]["feature"].location.strand
            genbank_gene_len = end_pos - start_pos
            if genbank_gene_len != excel_gene_len:
                unmatch_gene_len.append(row["Locus tag"])
            
            if row["Strand"] == "-":
                if strand == -1:
                    if (end_pos + utr_len) != tss_pos:
                        diff = tss_pos - utr_len - end_pos
                        unmatch_start_neg.append((row["Locus tag"], diff))
                elif strand == 1:
                    unmatch_direction.append(row["Locus tag"])
                else:
                    raise Exception("unknown strand direction in gen bank file: {}".format(strand))
            elif row["Strand"] == "+":
                if strand == 1:
                    if (start_pos + 1 - utr_len) != tss_pos:
                        diff = tss_pos + utr_len - start_pos - 1
                        unmatch_start_pos.append((row["Locus tag"], diff))
                elif strand == -1:
                    unmatch_direction.append(row["Locus tag"])
                else:
                    raise Exception("unknown strand direction in gen bank file: {}".format(strand))
            else:
                raise Exception("unknown strand direction in the Excel file: {}".format(row["Strand"]))
                
                 
    f.apply(apply_fcn, axis = 1)
    print("There are {} locus tags in the Excel file that were not found in the genome".format(len(unfound_tags)))
    print("These are the unfound tags: ")
    print(unfound_tags)    
    
    print("There are {} genes that have contradicting directions".format(len(unmatch_direction)))
    print("These are the genes with contradicting directions")
    print(unmatch_direction)
    
    print("Additionally, there are {} of genes that were unmatching in length".format(len(unmatch_gene_len)))
    #print("These are the genes with unmatching lengths")
    #print(unmatch_gene_len)
    
    print("There are also {} genes that didn't have matching start sites in the top strand".format(len(unmatch_start_pos)))
    unmatch_start_tags = set()
    diff_pos = []
    for (tag, d) in unmatch_start_pos:
        unmatch_start_tags.add(tag)
        diff_pos.append(d)
    plt.figure()
    plt.title("diff pos")
    plt.hist(diff_pos)
    
    print("There are also {} genes that didn't have matching start sites in the bottom strand".format(len(unmatch_start_neg)))
    diff_neg = []
    for (tag, d) in unmatch_start_neg:
        unmatch_start_tags.add(tag)
        diff_neg.append(d)
    plt.figure()
    plt.title("diff neg")
    plt.hist(diff_neg)
          
    print("done searching")
    print("")
    unmatch_start_len = set(unmatch_gene_len).intersection(unmatch_start_tags)
    unmatch_start_tags 
    return (unmatch_start_tags - unmatch_start_len)
    
unmatch_start_tags = test_gene(excel_f)

def compare_version(unmatch_start_tags):
    unmatch_version_tags = []
    unmatch_tags = [(tag[:2] + "_" + tag[2:]) for tag in unmatch_start_tags]
    unmatch_tags = set(unmatch_tags)
    v3 = "putida_v3.gb"
    v4 = "Putida KT2440.gbff"
    v3_record = []
    v4_record = []
    for seq_record3 in SeqIO.parse(v3, "genbank"):
        for feature in seq_record3.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers["locus_tag"][0]
                if locus_tag in unmatch_tags:
                    v3_record.append(feature)
    for seq_record4 in SeqIO.parse(v4, "genbank"):
        for feature in seq_record4.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers["locus_tag"][0]
                if locus_tag in unmatch_tags:
                    v4_record.append(feature)
    
    for (tag3, tag4) in zip(v3_record, v4_record):
        if tag3.qualifiers["locus_tag"][0] != tag4.qualifiers["locus_tag"][0]:
            raise Exception("unmatching tags")
        strand3 = tag3.location.strand
        strand4 = tag4.location.strand
        start3 = tag3.location.start
        start4 = tag4.location.start
        end3 = tag3.location.end
        end4 = tag4.location.end
        if end3-start3 != end4-start4:
            raise Exception("unequal length")
        if strand3 != strand4:
            print("not the same strand")
        if strand3 == -1:            
            if seq_record3[start3:end3+200].seq._data != seq_record4[start4:end4+200].seq._data:
                unmatch_version_tags.append(tag3)
        elif strand3 == 1:            
            if seq_record3[start3-200:end3].seq._data != seq_record4[start4-200:end4].seq._data:
                unmatch_version_tags.append(tag3)
        else:
            raise Exception("wrong strand notation: {}".format(strand3))
    print(len(unmatch_version_tags))
    print("")
    print("")
    print("")

# compare_version(unmatch_start_tags)
    
def filter_fun(row, seq):
    offset = 100
    tag = row["Locus tag"][:2] + "_" + row["Locus tag"][2:]
    tss_pos = row["Position"]
    if tag in seq:
        gene_num = seq[tag]["gene_num"]
        if seq[tag]["feature"].strand == -1:
            adjacent_gene_num = gene_num + 1
            adjacent_gene_tag = list_of_genes[adjacent_gene_num - 1]
            start_pos = seq[adjacent_gene_tag]["feature"].location.start
            return True if tss_pos < (start_pos - offset) else False
        else:
            adjacent_gene_num = gene_num - 1
            adjacent_gene_tag = list_of_genes[adjacent_gene_num - 1]
            end_pos = seq[adjacent_gene_tag]["feature"].location.end   
            return True if tss_pos > (end_pos + offset) else False
    else:
        return False

NGG_count = dict()
CCN_count = dict()
def get_upstream_sequence(row, seq, ngg_range, ccn_range):
    upstream_bp = 110
    pos = row["Position"]
    
    if row["Strand"] == "-":
        upstream_seq = seq_record[((pos + 1) - 1) : ((pos + 1)- 1 + upstream_bp)].reverse_complement()
        seq_str = upstream_seq.seq._data
    elif row["Strand"] == "+":
        upstream_seq = seq_record[(pos - upstream_bp - 1) : (pos - 1)] #pos-1 -> to account for 0th indexing
        seq_str = upstream_seq.seq._data
    else:
        raise Exception("wrong format for strand, this is what are read: {}".format(row["strand"]))
    if len(seq_str) != upstream_bp:
        print("this is the length of the sequence string: {}".format(len(seq_str)))
        raise Exception("mismatch")
    
    pam_exist = False
    pam_num = 0
    pam_strand = ""
    pam_site = []
    
    for i in range(upstream_bp - 2):
        if seq_str[i] == "G":
            if seq_str[i+1] == "G":
                if i > 1:
                    pam_pos = upstream_bp - (i - 1) # location of N
                    try:
                            NGG_count[pam_pos] += 1
                    except:
                        NGG_count[pam_pos] = 1
                    if pam_pos in ngg_range:
                        pam_exist = True
                        pam_strand = row["Strand"]
                        genome_pos = pos + pam_pos if row["Strand"] == "-" else pos - pam_pos
                        pam_seq = seq_str[(i-1) : (i+2)]
                        pam_site.append((pam_seq, pam_pos, genome_pos))
                        
        elif seq_str[i] == "C":
            if seq_str[i+1] == "C":
                if i < upstream_bp - 2:
                    pam_pos = upstream_bp - (i + 3) # location of the nucleotide right after N
                    try:
                            CCN_count[pam_pos] += 1
                    except:
                        CCN_count[pam_pos] = 1
                    if pam_pos in ccn_range:
                        pam_exist = True
                        pam_strand = row["Strand"]
                        genome_pos = pos + pam_pos if row["Strand"] == "-" else pos - pam_pos
                        pam_seq = seq_str[i : (i+3)]
                        pam_site.append((pam_seq, pam_pos, genome_pos))
                        
        pam_num = len(pam_site)
    return pd.Series([row["Locus tag"], pam_exist, pam_num, pam_strand, pam_site], index=["Locus tag", "PAM exist", "PAM num", "PAM strand", "PAM sites"])

filtered_genes = excel_f[excel_f.apply(lambda row: filter_fun(row, seq), axis = 1)]

tol = 2
ngg_range = set()
ngg_pos = [range(70-tol, 70+tol+1), range(80-tol, 80+tol+1), range(90-tol, 90+tol+1)]
for n in ngg_pos:
    ngg_range = ngg_range.union(set(n))
ccn_range = set()
ccn_pos = [range(81-tol, 81+tol+1), range(91-tol, 91+tol+1), range(101-tol, 101+tol+1)]
for c in ccn_pos:
    ccn_range = ccn_range.union(set(c))

pam = filtered_genes.apply(lambda row: get_upstream_sequence(row, seq, ngg_range, ccn_range), axis = 1)
# filtered_genes.to_excel("viable_genes_offset_100.xlsx")
pam = pam[pam["PAM exist"]]
print(filtered_genes["Locus tag"])
print(len(seq))

plt.figure()
plt.bar(CCN_count.keys(), CCN_count.values())
plt.show()
plt.figure()
plt.bar(NGG_count.keys(), NGG_count.values())
plt.show()
# plt.plot(pam["Locus tag"], pam["PAM sites"])

filtered_tags = set(pam["Locus tag"])
invalid_tags = filtered_tags.intersection(unmatch_start_tags)
print(invalid_tags)
print(len(invalid_tags))

