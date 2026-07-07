#!/usr/bin/env python3

import argparse
import os
import pandas as pd

from rbd_common import (
    get_sequence_region,
    InterproRBDInfo,
    merge_intervals,
    read_fasta_file,
    refine_intervals,
)


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Refine InterProScan results to keep only predicted RNA-binding domains (RBDs), 
    resolved to parent domain level and refined to remove overlapping intervals.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="refine_interproscan_results_keep_rbds.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--ips-res",
                   dest="in_interproscan_results",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input InterProScan results file")
    p.add_argument("--par-infos",
                   dest="in_parent_infos",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "File containing parent domain information")
    p.add_argument("--out",
                   dest="out_table",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output table file to store refined results")
    p.add_argument("--out-format",
                   dest="out_format",
                   type=int,
                   metavar='int',
                   choices=[1, 2],
                   default=1,
                   help="Output format. 1: one row per query with summarized hits. 2: one row per RBD hit. (default: 1)")
    p.add_argument("--fasta",
                   dest="in_fasta",
                   type=str,
                   metavar='str',
                   default=None,
                   help="Optional input FASTA file. When provided, protein sequences are added to the output table as fasta_seq.")
    p.add_argument("--eval-thr",
                   dest="eval_thr",
                   type=float,
                   metavar='float',
                   default=1e-5,
                   help="E-value threshold for filtering InterProScan hits (default: 1e-5)")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_interproscan_results), "--ips-res file \"%s\" not found" % (args.in_interproscan_results)
    assert os.path.exists(args.in_parent_infos), "--par-infos file \"%s\" not found" % (args.in_parent_infos)
    if args.in_fasta:
        assert os.path.exists(args.in_fasta), "--fasta file \"%s\" not found" % (args.in_fasta)

    qid2fasta_seq = {}
    if args.in_fasta:
        print("Reading in FASTA file ...")
        qid2fasta_seq = read_fasta_file(args.in_fasta)
        print("# FASTA sequences read in: %d" % len(qid2fasta_seq))

    """

    ================================
    RBD PARENT INFOS
    ================================


    INTERPRO RNA-binding domain (RBD) parent infos file format:

    data/rbp2go.selected_new_rbds.reform.added.parent_infos.txt

    interpro_id	short	name	type	parent_id	parent_short	parent_name	parent_type
    IPR000061	Surp	SWAP/Surp	domain	IPR000061	Surp	SWAP/Surp	domain
    IPR000196	Ribosomal_eL19_dom	Large ribosomal subunit protein eL19 domain	domain	IPR000196	Ribosomal_eL19_dom	Large ribosomal subunit protein eL19 domain	domain
    IPR000197	Znf_TAZ	Zinc finger, TAZ-type	domain	IPR000197	Znf_TAZ	Zinc finger, TAZ-type	domain
    ...

    """

    # Read in parent domain infos into dictionary with interpro_id -> InterproRBDInfo mapping.
    print("Reading in parent domain infos ...")

    interpro2parent_infos = {}

    with open(args.in_parent_infos, "r") as f:

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

                interpro2parent_infos[interpro_id] = InterproRBDInfo(interpro_id, interpro_name, interpro_type, interpro_short, parent_id, parent_name, parent_type, parent_short)

    print("# interpro2parent infos read in: %d" % len(interpro2parent_infos))


    """

    ====================================
    INTERPROSCAN OUTPUT FORMAT
    ====================================


    ./interproscan.sh -i test_protein_seqs.fa -f tsv -dp -cpu 8 -o test_protein_seqs.interproscan_results.out

    test_protein_seqs.interproscan_results.out format:
    RNCMPT00033,IGF2BP2	e94c5ac706d4dd67d109aaaa5789b357	599	SMART	SM00322	kh_6	273	346	2.7E-14	T	07-07-2026	IPR004087	K Homology domain	-	-
    RNCMPT00033,IGF2BP2	e94c5ac706d4dd67d109aaaa5789b357	599	SMART	SM00322	kh_6	508	580	4.3E-12	T	07-07-2026	IPR004087	K Homology domain	-	-
    RNCMPT00033,IGF2BP2	e94c5ac706d4dd67d109aaaa5789b357	599	SMART	SM00322	kh_6	192	263	2.2E-14	T	07-07-2026	IPR004087	K Homology domain	-	-
    RNCMPT00033,IGF2BP2	e94c5ac706d4dd67d109aaaa5789b357	599	SMART	SM00322	kh_6	426	497	1.0E-12	T	07-07-2026	IPR004087	K Homology domain	-	-
    ...

    column meanings:
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


    qid2pid2hit_pos_dic = {}  # query_id -> parent_interpro_id -> list of [start, end] domain hit positions.
    ipid2short_dic = {}  # parent_interpro_id -> short name mapping.

    with open(args.in_interproscan_results, "r") as f:

        for line in f:

            line = line.strip().split("\t")
            qid = line[0]  # Query ID.
            # uniprot_id, pdb_id, type_id = query_id.split("_")
            # upidpdbid = query_id
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
            if eval > args.eval_thr:
                continue

            parent_info = interpro2parent_infos[interpro_id]
            parent_id = parent_info.parent_id
            parent_name = parent_info.parent_name
            parent_type = parent_info.parent_type
            parent_short = parent_info.parent_short

            if qid not in qid2pid2hit_pos_dic:
                qid2pid2hit_pos_dic[qid] = {parent_id: [[start, end]]}
            else:
                if parent_id not in qid2pid2hit_pos_dic[qid]:
                    qid2pid2hit_pos_dic[qid][parent_id] = [[start, end]]
                else:
                    # if not rbplib.is_overlap(qid2pid2hit_pos_dic[qid][parent_id], [start, end]):
                    qid2pid2hit_pos_dic[qid][parent_id].append([start, end])

            ipid2short_dic[parent_id] = parent_short

    """
    Refine intervals.

    Filter list of start, end pairs to remove intervals that are fully 
    contained within another interval.

    """
    for qid in qid2pid2hit_pos_dic:

        for interpro_id in qid2pid2hit_pos_dic[qid]:

            qid2pid2hit_pos_dic[qid][interpro_id] = refine_intervals(qid2pid2hit_pos_dic[qid][interpro_id])


    """
    Merge intervals.

    Merge overlapping intervals, and let the new start and end positions be the 
    minimum and maximum start and end positions of the overlapping intervals, 
    respectively.
    
    """
    for qid in qid2pid2hit_pos_dic:

        for interpro_id in qid2pid2hit_pos_dic[qid]:

            qid2pid2hit_pos_dic[qid][interpro_id] = merge_intervals(qid2pid2hit_pos_dic[qid][interpro_id],
                                                                    overlap_threshold=0.6)

    """
    Output refined InterProScan results to file.
    
    """

    with open(args.out_table, "w") as f:

        if args.out_format == 1:
            header = ["query_id", "interpro_ids", "domain_names", "domain_counts", "hit_pos"]
            if args.in_fasta:
                header.append("fasta_seq")
            f.write("\t".join(header) + "\n")

            for qid in qid2pid2hit_pos_dic:

                ipid2c_dic = {}

                for ip_id in qid2pid2hit_pos_dic[qid]:
                    ipid2c_dic[ip_id] = len(qid2pid2hit_pos_dic[qid][ip_id])

                ip_ids = []
                short_ids = []
                counts = []
                for ip_id in ipid2c_dic:
                    ip_ids.append(ip_id)
                    short_ids.append(ipid2short_dic[ip_id])
                    counts.append(ipid2c_dic[ip_id])
                ip_ids_str = ";".join(ip_ids)
                short_ids_str = ";".join(short_ids)
                counts_str = ";".join([str(c) for c in counts])

                hit_pos_str = []
                for ip_id in sorted(qid2pid2hit_pos_dic[qid].keys()):
                    hit_pos_str.append(f"{ip_id}:{','.join([f'{start}-{end}' for start, end in qid2pid2hit_pos_dic[qid][ip_id]])}")
                hit_pos_str = ";".join(hit_pos_str)

                out_fields = [qid, ip_ids_str, short_ids_str, counts_str, hit_pos_str]
                if args.in_fasta:
                    out_fields.append(qid2fasta_seq.get(qid, ""))
                f.write("\t".join(out_fields) + "\n")
        elif args.out_format == 2:
            header = ["query_id", "interpro_id", "domain_name", "hit_pos_s", "hit_pos_e"]
            if args.in_fasta:
                header.append("rbd_seq")
            f.write("\t".join(header) + "\n")

            for qid in qid2pid2hit_pos_dic:
                fasta_seq = qid2fasta_seq.get(qid, "")
                for ip_id in qid2pid2hit_pos_dic[qid]:
                    for start, end in qid2pid2hit_pos_dic[qid][ip_id]:
                        out_fields = [qid, ip_id, ipid2short_dic[ip_id], str(start), str(end)]
                        if args.in_fasta:
                            out_fields.append(get_sequence_region(fasta_seq, start, end) if fasta_seq else "")
                        f.write("\t".join(out_fields) + "\n")

    print("# queries with final RBD annotations: %d" % len(qid2pid2hit_pos_dic))
    if args.in_fasta:
        n_missing_fasta_seq = sum(1 for qid in qid2pid2hit_pos_dic if qid not in qid2fasta_seq)
        print("# queries with final RBD annotations missing FASTA sequence: %d" % n_missing_fasta_seq)
