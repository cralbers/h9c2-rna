# creating decoy aware salmon transcriptome index

```
module load salmon/1.2.1
```

## get names of genome targets
```
grep "^>" <(gunzip -c Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```

## concatenate genome and transcriptome
```
cat Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz > rat_gentrome.fa.gz

```


## make index
```
salmon index -t rat_gentrome.fa.gz --decoys decoys.txt -p 12 -i rat_salmon_index
```
