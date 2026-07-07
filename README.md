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


## Annotate RBDs

Now let's bring it all together (InterProScan output + RBD infos) ...





