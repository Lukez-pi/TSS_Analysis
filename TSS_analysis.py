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
        if feature.type == "CDS":
            temp = feature.qualifiers["locus_tag"][0].split("_")
            locus_tag = temp[0] + temp[1]
            strand = feature.strand
            if strand == -1:
                start_position = feature.location.end
                end_position = feature.location.start
            else:
                start_position = feature.location.start
                end_position = feature.location.end

            previous_start = start_position
            
            seq[locus_tag] = {"feature": feature, "gene_num": gene_counter}
            list_of_genes.append(locus_tag)
            gene_counter += 1
            
    print("total number of genes")
    print(gene_counter-1)

def test_gene(f):
    # f is the excel file
    unfound_tags = []
    overlapping_genes = []
    unmatch_gene_len = []
    unmatch_start_pos = []
    unmatch_start_neg = []
    unmatch_direction = []
    def apply_fcn(row):
        tag = row["Locus tag"]
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
            
            if seq[tag]["feature"].location.start >= seq[tag]["feature"].location.end:
                raise Exception("start and end position doesn't make sense")
                
            # Find genes that are overlapping
            if start_position <= previous_start:
                overlapping_genes.append(tag)
            
            # Find tags that exist in Excel but not in the Genbank file
            if genbank_gene_len != excel_gene_len:
                unmatch_gene_len.append(row["Locus tag"])
                
            # Find genes that have different starting positions and directions
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
    
    print("There are also {} genes that didn't have matching start sites in the bottom strand".format(len(unmatch_start_neg)))
    diff_neg = []
    for (tag, d) in unmatch_start_neg:
        unmatch_start_tags.add(tag)
        diff_neg.append(d)
    # plot the number of base pairs that the excel data and the gen bank data differ by
# =============================================================================
#     plt.figure()
#     plt.title("diff pos")
#     plt.hist(diff_pos)
#         
#     plt.figure()
#     plt.title("diff neg")
#     plt.hist(diff_neg)
# =============================================================================
    
    unmatch_start_match_len = unmatch_start_tags - set(unmatch_gene_len)
    unmatch_start_match_len_match_dir = unmatch_start_match_len - set(unmatch_direction)
    print("\nOut of the {} genes that are unmatching in start position, there are {} genes that are matching in length".format(len(unmatch_start_tags), len(unmatch_start_match_len)))
    print("Additionally, when we remove the genes that are also unmatching in direction")
    print("This leaves us with {} genes to be matching in length, direction but not starting sites\n".format(len(unmatch_start_match_len_match_dir)))
    return [set(unmatch_gene_len), unmatch_start_match_len_match_dir]
    
[unmatch_gene_len, unmatching_tags] = test_gene(excel_f) 

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
    
    counter = 0
    print("\n")
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
            if seq_record3[start3:end3+300].seq._data != seq_record4[start4:end4+300].seq._data:
                tag_name = tag3.qualifiers["locus_tag"][0]
                new_tag_name = tag_name.split("_")[0] + tag_name.split("_")[1]
                unmatch_version_tags.append(new_tag_name)
        elif strand3 == 1:            
            if seq_record3[start3-300:end3].seq._data != seq_record4[start4-300:end4].seq._data:
                tag_name = tag3.qualifiers["locus_tag"][0]
                new_tag_name = tag_name.split("_")[0] + tag_name.split("_")[1]
                unmatch_version_tags.append(new_tag_name)
        else:
            raise Exception("wrong strand notation: {}".format(strand3))
        counter += 1
        if counter % 50 == 0:
            print("we are currently at iteration {}....".format(counter))
    print("\nAfter comparing v3 to v4, we found that out of the genes with valid pam sites but different starting positions,")
    print("there are {} genes that have a different sequence".format(len(unmatch_version_tags)))
    print(unmatch_version_tags)
    for t in unmatch_version_tags:
        print(t)
        
    # return all genes that should be fine, despite a different starting position
    return unmatch_start_tags - set(unmatch_version_tags)
    
    
def filter_by_intergenic_len(row, seq):
    offset = 100
    tag = row["Locus tag"]
    tss_pos = row["Position"]
    if tag in seq:
        gene_num = seq[tag]["gene_num"]
        if seq[tag]["feature"].strand == -1:
            adjacent_gene_num = gene_num + 1
            adjacent_gene_tag = list_of_genes[adjacent_gene_num - 1]
            start_pos = seq[adjacent_gene_tag]["feature"].location.start
            return True if tss_pos < (start_pos - offset) else False
        elif seq[tag]["feature"].strand == 1:
            adjacent_gene_num = gene_num - 1
            adjacent_gene_tag = list_of_genes[adjacent_gene_num - 1]
            end_pos = seq[adjacent_gene_tag]["feature"].location.end   
            return True if tss_pos > (end_pos + offset) else False
        else:
            raise Exception("wrong strand input format: {}".format(seq[tag]["feature"].strand))
    else:
        return False

NGG_count = dict()
CCN_count = dict()

def find_PAM(row, seq, ngg_range, ccn_range):
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
    return pd.Series([row["Locus tag"], row["Product"], row["UTR length"], pam_exist, pam_num, pam_strand, pam_site], index=["Locus tag", "description", "UTR length", "PAM exist", "PAM num", "PAM strand", "PAM sites"])

filtered_genes = excel_f[excel_f.apply(lambda row: filter_by_intergenic_len(row, seq), axis = 1)]
print("There are {} genes with a intergenic region that are long enough".format(len(filtered_genes)))

tol = 2
ngg_range = set()
ngg_pos = [range(70-tol, 70+tol+1), range(80-tol, 80+tol+1), range(90-tol, 90+tol+1)]
for n in ngg_pos:
    ngg_range = ngg_range.union(set(n))
ccn_range = set()
ccn_pos = [range(81-tol, 81+tol+1), range(91-tol, 91+tol+1), range(101-tol, 101+tol+1)]
for c in ccn_pos:
    ccn_range = ccn_range.union(set(c))

pam_filter = filtered_genes.apply(lambda row: find_PAM(row, seq, ngg_range, ccn_range), axis = 1)
# filtered_genes.to_excel("viable_genes_offset_100.xlsx")

genes_with_pam = pam_filter[pam_filter["PAM exist"]]

# =============================================================================
# plt.figure()
# plt.bar(CCN_count.keys(), CCN_count.values())
# plt.show()
# plt.figure()
# plt.bar(NGG_count.keys(), NGG_count.values())
# plt.show()
# =============================================================================

genes_with_pam_tag = set(genes_with_pam["Locus tag"])
print(len(genes_with_pam_tag))
genes_with_pam_tag = genes_with_pam_tag - unmatch_gene_len
print("after removing genes with unequal length: {}\n".format(len(genes_with_pam_tag)))
# unmatching tags are genes with matching length and dir but not start
genes_with_pam_unmatching = genes_with_pam_tag.intersection(unmatching_tags) # find genes with pam sites but are not matching in the Excel and genbank file
print("There are {} genes with a PAM site at a valid location".format(len(genes_with_pam_tag)))
print("When taking the intersect of genes with PAM sites and genes that are not matching, there are {} results".format(len(genes_with_pam_unmatching)))
print("Therefore these genes need to be analyzed and see if their sequence are the same between different versions")

genes_passed_comp = compare_version(genes_with_pam_unmatching)
valid_gene = (genes_with_pam_tag - genes_with_pam_unmatching).union(genes_passed_comp)

def update_TSS_pos(df):
    utr_len = df["UTR length"]
    pam_strand = df["PAM strand"]
    if pam_strand == "-":
        
        print("hello")
    elif pam_strand == "+":
        print("hello")
    else:
        yout
        raise Exception("wrong strand format: {}".format(pam_strand))
    print("hello")
    
genes_pam_final = genes_with_pam[genes_with_pam["Locus tag"].isin(valid_gene)]



