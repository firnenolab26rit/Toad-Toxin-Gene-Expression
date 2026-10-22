# Commands to generate a gene to transcript mapping from tr2aacds.pl output 

## Using the *.okayset.mrna fasta file
1.) Extract the oid=<value> from the fasta headers

```
grep -oP 'oid=\KTRINITY_DN\d+_c\d+_g\d+_i\d+' file.fasta > trinity_ids.txt
```

2.) Filter the trinity *.gene_map file to retain only the trinity ids that pasted the tr2aacds.pl pipeline. 

* the *okay.mrna file have duplicate trinity (oid=) values. Not sure why this is but its only like a 1000 so I think I have to check tr2aacds.pl classification status.
will just involve adding a additional field and skimming through the 1000 or so files. 

* the duplicates are because tr2aacds.pl puts genes into different classifications (primary/alternate). We want to keep all these in the dataset. It is probably minor
but its better to have. 

```
# Counting duplicated 
sort valid_ids.txt | uniq -c | awk '$1 > 1' > duplicates_with_counts.txt

# Check if duplicates exist 
sort valid_ids.txt | uniq | wc -l	# Lines with dupicates removed
wc -l valid_ids.txt			# Total lines 
```

3.) Extract the tr2aacds id (unique) and oid (duplicates based on above). Then convert the oid to the gene level per trinity. (just removes _i[0-9]) 

```
# Get EV id and oid from fasta headers
awk -F '[ >;]+' '
/^>/{
    evg = $2
    trinity = gensub(/.*oid=/, "", "g", $0)
    print trinity "\t" evg
}' yourfile


# Remove the lagging ";" from column 1
awk '{ gsub(/;/, "", $1); print }' input.txt > output.txt


# Remove the everything after _i and _i itself (isoform key) from column 1
awk 'BEGIN{OFS="\t"} { sub(/_i.*/, "", $1); print }' file > out.file
```


Now you have a gene to transcript level mapping. 


# Cleaning up the fasta file 
* I wanted to get rid of all the other information in the header and just keep the EV id. 
```
sed -i '/^>/ s/ .*$//' A_americanus_trinity.okay.mrna
```
