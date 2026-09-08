
from rbpkit import rbplib
import pandas as pd
import sys

# print(is_overlap([[1, 10], [20, 30]], [5, 15]))
# print(is_overlap([[1, 100], [200, 300]], [5, 95]))

# sys.exit()

"""
in_file1:
interpro_id	short	name	type	parent_id	parent_short	parent_name	parent_type
IPR000061	Surp	SWAP/Surp	domain	IPR000061	Surp	SWAP/Surp	domain
IPR000196	Ribosomal_eL19_dom	Large ribosomal subunit protein eL19 domain	domain	IPR000196	Ribosomal_eL19_dom	Large ribosomal subunit protein eL19 domain	domain
IPR000197	Znf_TAZ	Zinc finger, TAZ-type	domain	IPR000197	Znf_TAZ	Zinc finger, TAZ-type	domain
IPR000241	RlmKL-like_Mtase	Ribosomal RNA large subunit methyltransferase K/L-like, methyltransferase domain	domain	IPR000241	RlmKL-like_Mtase	Ribosomal RNA large subunit methyltransferase K/L-like, methyltransferase domain	domain
...

"""


in_file1 = "data/rbp2go.selected_new_rbds.reform.updated.parent_infos.txt"
# a2021nonX_exp_protein_rna.out
# a2021X_exp_protein_rna.out
# interproscan_a2021
# b2021X_exp_protein_rna
# in_file2 = "interproscan_b2021x/b2021X_exp_protein_rna.out"
# out_file = "interproscan_b2021x/b2021X_exp_protein_rna.interpro_rbd_domains.txt"


# domain_a2021nX, domain_a2021x, domain_b2021nX, domain_b2021x
# a2021nonX, a2021x, b2021nonX, b2021x
data_id = "b2021x"
res_folder = "add_rbd_annot_final"


in_file2 = "%s/%s_exp_protein_rna.out" % (res_folder, data_id)
fasta_file = "%s/%s_exp_protein_rna.interproscan.fasta" % (res_folder, data_id)

out_file = "%s/%s_exp_protein_rna.interpro_rbd_domains.txt" % (res_folder, data_id)

# Threshold for filtering interproscan hits by e-value.
eval_thr = 1e-5

print("Reading in fasta file ...")
seqs_dic = rbplib.read_in_fasta_file(fasta_file)  # seq_id format: 9IJ3_0
assert seqs_dic, "No sequences read in."
print("# seqs_dic: %d" % len(seqs_dic))

# Read in parent domain infos into dictionary with interpro_id -> InterproRBDInfo mapping.
print("Reading in parent domain infos ...")

interpro2parent_infos = {}

with open(in_file1, "r") as f:

    for line in f:
        line = line.strip()
        if line.startswith("interpro_id"):
            continue
        else:
            line = line.split("\t")
            interpro_id = line[0]
            interpro_short = line[1]
            interpro_name = line[2]
            interpro_type = line[3]
            parent_id = line[4]
            parent_short = line[5]
            parent_name = line[6]
            parent_type = line[7]

            assert interpro_id not in interpro2parent_infos, f"Duplicate interpro_id: {interpro_id}"

            interpro2parent_infos[interpro_id] = rbplib.InterproRBDInfo(interpro_id, interpro_name, interpro_type, interpro_short, parent_id, parent_name, parent_type, parent_short)

print("# interpro2parent_infos: %d" % len(interpro2parent_infos))


"""
in_file2:
8K0K_7	b461a0804c1547153986578ece03aa04	303	NCBIfam	TIGR02566	type I-F CRISPR-associated protein Csy3	5	302	3.3E-67	T	17-02-2025	IPR013399	Type I-F CRISPR-associated protein Csy3	-	-
8K0K_7	b461a0804c1547153986578ece03aa04	303	Pfam	PF09615	CRISPR-associated protein (Cas_Csy3)	7	302	1.4E-80	T	17-02-2025	IPR013399	Type I-F CRISPR-associated protein Csy3	-	-
8K0K_6	b461a0804c1547153986578ece03aa04	303	NCBIfam	TIGR02566	type I-F CRISPR-associated protein Csy3	5	302	3.3E-67	T	17-02-2025	IPR013399	Type I-F CRISPR-associated protein Csy3	-	-
8K0K_6	b461a0804c1547153986578ece03aa04	303	Pfam	PF09615	CRISPR-associated protein (Cas_Csy3)	7	302	1.4E-80	T	17-02-2025	IPR013399	Type I-F CRISPR-associated protein Csy3	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	ProSiteProfiles	PS50137	Double stranded RNA-binding domain (dsRBD) profile.	12	80	19.563862	T	17-02-2025	IPR014720	Double-stranded RNA-binding domain	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	FunFam	G3DSA:3.30.160.20:FF:000005	Putative double-stranded RNA-specific adenosine deaminase	1	82	8.1E-39	T	17-02-2025	-	-	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	SMART	SM00358	DRBM_3	13	79	1.8E-17	T	17-02-2025	IPR014720	Double-stranded RNA-binding domain	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	CDD	cd19915	DSRM_DRADA_rpt3	11	81	3.24504E-48	T	17-02-2025	IPR044457	ADAR1, third double-stranded RNA binding domain	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	Pfam	PF00035	Double-stranded RNA binding motif	13	78	6.9E-9	T	17-02-2025	IPR014720	Double-stranded RNA-binding domain	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	SUPERFAMILY	SSF54768	dsRNA-binding domain-like	7	81	5.13E-19	T	17-02-2025	-	-	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	Gene3D	G3DSA:3.30.160.20	-	1	82	5.9E-25	T	17-02-2025	-	-	-	-
7ZLQ_2	37c277c76cdc388b06f6a5a68cf97dd8	82	PANTHER	PTHR10910	EUKARYOTE SPECIFIC DSRNA BINDING PROTEIN	9	81	7.4E-26	T	17-02-2025	-	-	-	-
7ZLQ_1	37c277c76cdc388b06f6a5a68cf97dd8	82	ProSiteProfiles	PS50137	Double stranded RNA-binding domain (dsRBD) profile.	12	80	19.563862	T	17-02-2025	IPR014720	Double-stranded RNA-binding domain	-	-

in_file2 column format:
    Protein accession (e.g. P51587)
    Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    Sequence length (e.g. 3418)
    Analysis (e.g. Pfam / PRINTS / Gene3D)
    Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
    Signature description (e.g. BRCA2 repeat profile)
    Start location
    Stop location
    Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
    Status - is the status of the match (T: true)
    Date - is the date of the run
    InterPro annotations - accession (e.g. IPR002093)
    InterPro annotations - description (e.g. BRCA2 repeat)
    GO annotations with their source(s), e.g. GO:0005515(InterPro)|GO:0006302(PANTHER)|GO:0007195(InterPro,PANTHER). This is an optional column; only displayed if the --goterms option is switched on
    Pathways annotations, e.g. REACT_71. This is an optional column; only displayed if the --pathways option is switched on

"""


pdbid2interpro_ids_dic = {}
pdbid2hit_pos_dic = {}
pdbid2sig_acc_dic = {}
pdbid2type_dic = {}
interproid2short_dic = {}



with open(in_file2, "r") as f:

    for line in f:

        line = line.strip().split("\t")
        query_id = line[0]  # Format: 3OVB_1
        pdb_id, type_id = query_id.split("_")
        # upidpdbid = uniprot_id + "," + pdb_id
        seq_len = line[2]
        analysis = line[3]
        signature_acc = line[4]
        signature_desc = line[5]
        start = int(line[6])
        end = int(line[7])
        eval = line[8]
        if eval == "-":
            continue
        eval = float(eval)
        interpro_id = line[11]
        interpro_desc = line[12]

        if interpro_id not in interpro2parent_infos:
            continue
        if eval > eval_thr:
            continue

        if interpro_id in interpro2parent_infos:
            parent_info = interpro2parent_infos[interpro_id]
            parent_id = parent_info.parent_id
            parent_name = parent_info.parent_name
            parent_type = parent_info.parent_type
            parent_short = parent_info.parent_short

            pdbid2type_dic[pdb_id] = type_id

            if pdb_id not in pdbid2hit_pos_dic:
                pdbid2hit_pos_dic[pdb_id] = {parent_id: [[start, end]]}
            else:
                if parent_id not in pdbid2hit_pos_dic[pdb_id]:
                    pdbid2hit_pos_dic[pdb_id][parent_id] = [[start, end]]
                else:
                    #if not rbplib.is_overlap(pdbid2hit_pos_dic[pdb_id][parent_id], [start, end]):
                    pdbid2hit_pos_dic[pdb_id][parent_id].append([start, end])

            interproid2short_dic[parent_id] = parent_short

"""
Refine intervals.

refine_intervals(intervals)

"""
for pdb_id in pdbid2hit_pos_dic:

    for interpro_id in pdbid2hit_pos_dic[pdb_id]:

        pdbid2hit_pos_dic[pdb_id][interpro_id] = rbplib.refine_intervals(pdbid2hit_pos_dic[pdb_id][interpro_id])


"""
Merge intervals.

merge_intervals(intervals)

"""
for pdb_id in pdbid2hit_pos_dic:

    for interpro_id in pdbid2hit_pos_dic[pdb_id]:

        pdbid2hit_pos_dic[pdb_id][interpro_id] = rbplib.merge_intervals(pdbid2hit_pos_dic[pdb_id][interpro_id],
                                                                        overlap_threshold=0.6)


"""
out_file:
pdb_id\tuniprot_id\tinterpro_ids\tdomain_names\thit_pos

"""

with open(out_file, "w") as f:

    f.write("pdb_id\tsingle_chain\tinterpro_ids\tdomain_names\tdomain_counts\thit_pos\n")

    for pdb_id in pdbid2hit_pos_dic:

        # uniprot_id, pdb_id = upidpdbid.split(",")

        type_id = pdbid2type_dic[pdb_id]
        single_chain = 0
        if type_id == "0":
            single_chain = 1
        # seq_id = pdb_id + "_" + type_id
        # seq = seqs_dic[seq_id]

        ipid2c_dic = {}

        for ip_id in pdbid2hit_pos_dic[pdb_id]:
            ipid2c_dic[ip_id] = len(pdbid2hit_pos_dic[pdb_id][ip_id])

        ip_ids = []
        short_ids = []
        counts = []
        for ip_id in ipid2c_dic:
            ip_ids.append(ip_id)
            short_ids.append(interproid2short_dic[ip_id])
            counts.append(ipid2c_dic[ip_id])
        ip_ids_str = ";".join(ip_ids)
        short_ids_str = ";".join(short_ids)
        counts_str = ";".join([str(c) for c in counts])

        hit_pos_str = "-"
        if type_id == "0":
            hit_pos_str = []
            for ip_id in sorted(pdbid2hit_pos_dic[pdb_id].keys()):
                hit_pos_str.append(f"{ip_id}:{','.join([f'{start}-{end}' for start, end in pdbid2hit_pos_dic[pdb_id][ip_id]])}")
            hit_pos_str = ";".join(hit_pos_str)

        f.write(f"{pdb_id}\t{single_chain}\t{ip_ids_str}\t{short_ids_str}\t{counts_str}\t{hit_pos_str}\n")