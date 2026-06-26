# rbd-annotation
Annotate RNA Binding Domains (RBDs) in Protein Sequences

## InterProScan

For protein domain annotation, we can use [InterProScan](https://github.com/ebi-pf-team/interproscan). We choose the latest [Release](https://github.com/ebi-pf-team/interproscan/releases) and download it:

```
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
```
Installation requirements are also described [here](https://interproscan-docs.readthedocs.io/en/v5/InstallationRequirements.html) (Perl 5, Python 3, Java 11). First Java Runtime Environment 11:

```
sudo apt update && sudo apt install openjdk-11-jre
java -version
```

Now continue with InterProScan:

```
tar -xvzf interproscan-5.78-109.0-64-bit.tar.gz
cd interproscan-5.78-109.0/
chmod +x interproscan.sh
python3 setup.py -f interproscan.properties
```

A test call would look like this:

```
./interproscan.sh -i some_protein_sequences.fa -f tsv -dp -cpu 4 -o some_protein_sequences.out
```

## Get RBD information

To get RNA Binding Domain (RBD) information, we


## Annotate RBDs






