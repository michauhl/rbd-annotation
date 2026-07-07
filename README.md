# rbd-annotation
Annotate RNA Binding Domains (RBDs) in Protein Sequences

## InterProScan

For protein domain annotation, we can use [InterProScan](https://github.com/ebi-pf-team/interproscan). We choose the latest [Release](https://github.com/ebi-pf-team/interproscan/releases) and download it:

```
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
```

Compressed file size is ~ 7GB, uncompressed ~50GB. Installation requirements are also described [here](https://interproscan-docs.readthedocs.io/en/v5/InstallationRequirements.html) (Perl 5, Python 3, Java 11). First Java Runtime Environment 11:

```
sudo apt update && sudo apt install openjdk-11-jre
java -version
```

Now continue with InterProScan:

```
tar -xvzf interproscan-5.78-109.0-64-bit.tar.gz
mv interproscan-5.72-103.0 interproscan
cd interproscan/
chmod +x interproscan.sh
python3 setup.py -f interproscan.properties
```

A test call would look like this:

```
./interproscan.sh -i some_protein_sequences.fa -f tsv -dp -cpu 4 -o some_protein_sequences.out
```

## Get RBD information

To get RNA Binding Domain (RBD) information, we first download InterPro IDs of known RBDs from [RBP2GO](https://rbp2gov2.dkfz.de/) (`References & Data` &rarr; `Data & Downloads` &rarr; `RNA-related InterPro Features` &rarr; `Selected and new RBDs`):

```
wget https://rbp2gov2.dkfz.de/session/b3d021f298214163567fbfb86dedd8a3/download/Down.Table_RBDs_TXT?w= -O data/Table_selected_new_RBDs.txt
```

Reformat the table with:

```
python scripts/reformat_rbp2go_domain_table.py -i data/Table_selected_new_RBDs.txt -o data/rbp2go.selected_new_rbds.reform.txt
```

And add one missing domain manually for better results in later parent assignment:

```
echo -e "IPR047575\tDomain\tSm domain\tAdded" >> data/rbp2go.selected_new_rbds.reform.txt
```

Next, get parent domain information for these:

```
python scripts/interpro_get_parent_domain_infos.py --in data/rbp2go.selected_new_rbds.reform.txt --out data/rbp2go.selected_new_rbds.reform.parent_infos.txt
```

Or use custom file (with some domains manually added, in total 1000):
```
python scripts/interpro_get_parent_domain_infos.py --in data/rbp2go.selected_new_rbds.reform.added.txt --out data/rbp2go.selected_new_rbds.reform.added.parent_infos.txt
```


## Annotate RBDs

Now let's bring together the InterProScan results and the RBD domain infos we collected, to filter the InterProScan results to only keep RBD domain annotations:

```
python scripts/refine_interproscan_results_keep_rbds.py --ips-res data/test_protein_seqs.interproscan_results.out --par-infos data/rbp2go.selected_new_rbds.reform.added.parent_infos.txt --out data/test_protein_seqs.interpro_rbd_domains.txt --fasta data/test_protein_seqs.fa
```

Here InterPro IDs are resolved to parent IDs, and overlapping intervals are refined and merged in the following way:

1. Filter list of start+end domain hit positions to remove intervals that are fully contained within another interval.
2. Merge overlapping intervals, and let the new start and end positions be the minimum and maximum start and end positions of the overlapping intervals, respectively.

`--fasta` leads to adding the FASTA sequence of the query ID to `--out` table file. Use `--out-format 2` to output each domain hit in a single row in `--out` table.






