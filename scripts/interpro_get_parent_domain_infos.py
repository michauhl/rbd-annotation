#!/usr/bin/env python3

import argparse
import os

from rbd_common import (
    get_interpro_parent_infos,
)


"""
Get parent domain infos for downloaded RBP2GO RNA Protein Domain (RBD) InterPro IDs.

Input table (reformated RBP2GO table) format:

InterPro_ID	Type	InterPro_Name	Selection
IPR000061	Domain	SWAP/Surp	Selected
IPR000196	Domain	Ribosomal protein L19/L19e	Selected
IPR000197	Domain	Zinc finger, TAZ-type	Selected
IPR000241	Domain	Putative RNA methylase domain	Selected
IPR000253	Domain	Forkhead-associated (FHA) domain	Selected
...

"""


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a list of Interpro RBD (RNA-binding domain) IDs, get parent domain information.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="interpro_get_parent_domain_infos.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="in_ids",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input table file with InterPro RBD IDs in column 1 (data/rbp2go.selected_new_rbds.reform.txt)")
    p.add_argument("--out",
                   dest="out_parent_infos",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output table file storing parent domain infos for each input Interpro RBD ID")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_ids), "--in file \"%s\" not found" % (args.in_ids)

    interpro_ids_file = args.in_ids
    # selected_new_rbds.interpro_parent_infos.31-01-25.txt
    out_file = args.out_parent_infos

    print("Read in InterPro IDs ...")

    interpro_ids_dic = {}

    with open(interpro_ids_file, "r") as f:
        
        for line in f:
            line = line.strip()
            if line.startswith("InterPro_ID") or line.startswith("interpro_id"):
                continue
            else:
                line = line.split("\t")
                interpro_id = line[0]
                interpro_ids_dic[interpro_id] = 1

    assert interpro_ids_dic, "No InterPro IDs read in"
    print("# InterPro IDs read in:", len(interpro_ids_dic))

    print("Get parent domain infos ...")

    parent_infos_dic = get_interpro_parent_infos(interpro_ids_dic)

    print("# IDs with extracted parent domain infos:", len(parent_infos_dic))

    # Output to out_file.

    print("Write parent domain infos to:", out_file)

    with open(out_file, "w") as out_f:
        out_f.write("interpro_id\tshort\tname\ttype\tparent_id\tparent_short\tparent_name\tparent_type\n")
        for interpro_id in parent_infos_dic:
            interpro_name = parent_infos_dic[interpro_id].interpro_name
            interpro_type = parent_infos_dic[interpro_id].interpro_type
            interpro_short = parent_infos_dic[interpro_id].interpro_short
            parent_id = parent_infos_dic[interpro_id].parent_id
            parent_name = parent_infos_dic[interpro_id].parent_name
            parent_type = parent_infos_dic[interpro_id].parent_type
            parent_short = parent_infos_dic[interpro_id].parent_short

            out_f.write(f"{interpro_id}\t{interpro_short}\t{interpro_name}\t{interpro_type}\t{parent_id}\t{parent_short}\t{parent_name}\t{parent_type}\n")

    print("Done.")


"""

Example .json content:

https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR000504?format=json

.json file:

{
    "metadata": {
        "accession": "IPR000504",
        "entry_id": null,
        "type": "domain",
        "go_terms": [
            {
                "identifier": "GO:0003723",
                "name": "RNA binding",
                "category": {
                    "code": "F",
                    "name": "molecular_function"
                }
            }
        ],
        "source_database": "interpro",
        "member_databases": {
            "pfam": {
                "PF16367": "RNA recognition motif",
                "PF00076": "RNA recognition motif"
            },
            "profile": {
                "PS50102": "Eukaryotic RNA Recognition Motif (RRM) profile"
            },
            "smart": {
                "SM00360": "RNA recognition motif"
            }
        },
        "integrated": null,
        "hierarchy": {
            "accession": "IPR000504",
            "name": "RNA recognition motif domain",
            "type": "Domain",
            "children": [
                {
...

"""
