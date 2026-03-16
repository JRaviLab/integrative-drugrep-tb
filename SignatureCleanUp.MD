# For context Delete after PR is approved

### Author : Ling Thang

I ran the following script to clean up the file names first then to clean up the contents in the files. I also used git to restore any files that were unintentionally changed.

shell command to rename .tsv files by replacing every occurrence of "pmbc" with "PBMC", "dcs" with "dendritic_cells", and "THP-1" with "THP1" in the file names:

```shell
find . -type f -iname '*.tsv' -print0 |
  while IFS= read -r -d '' f; do
    dir=$(dirname "$f") # gets the directory of the file
    base=$(basename "$f") # gets the file name with extension
    newbase=$(printf '%s' "$base" | perl -pe 's/pmbc/PBMC/ig; s/dcs/dendritic_cells/ig; s/THP-1/THP1/ig') # renames
    target="$dir/$newbase"
    if [ "$f" != "$target" ]; then
      if [ -e "$target" ]; then
        echo "SKIP exists: $target" # skip if the target exists
      else
        mv -v -- "$f" "$target" # renames the file
      fi
    fi
  done
```

shell command to replace every occurrence of "pmbc" with "PBMC", "dcs" with "dendritic_cells", and "THP-1" with "THP1" in all .tsv files in the current directory and its subdirectories:

```shell
find . -name '*.tsv' -print0 | xargs -0 perl -pi -e 's/pmbc/PBMC/ig; s/dcs/dendritic_cells/ig; s/THP-1/THP1/ig'
```

this generate a lot of unintended changes, so I used git to restore the files that were unintentionally changed:

### terminal history:

```shell
10316  find . -name '*.tsv' -exec sed -i '' -e 's/pmbc/PBMC/g' -e 's/dcs/dendritic_cells/g' {} +
10318  git add . # manually checked each file that was renamed prior to adding all
10319  gcmsg "replace DCs -> Dendritic Cells & pmbc -> PBMC in /signatures"

# change dir
10324  find . -name '*.tsv' -exec sed -i '' -e 's/pmbc/PBMC/g' -e 's/dcs/dendritic_cells/g' {} +
10326  git restore metadata/Homo_sapiens.gene_info.tsv
10328  git add DE_results/microarray/aggregated_DE_signatures/
10330  gcmsg "replace DCs -> Dendritic Cells & pmbc -> PBMC in /aggregated_DE_signatures/"
10332  git add DE_results/microarray/circulating_aggregated_DE_signatures/
10334  gcmsg "replace DCs -> Dendritic Cells & pmbc -> PBMC in /circulating_aggregated_DE_signatures/"
10336  git add .
10337  gcmsg "replace DCs -> Dendritic Cells & pmbc -> PBMC in metadata_classifications.tsv"

# renaming files logic here performed at the /data
10353  find . -type f -iname '*.tsv' -print0 |\n  while IFS= read -r -d '' f; do\n    dir=$(dirname "$f")\n    base=$(basename "$f")\n    newbase=$(printf '%s' "$base" | perl -pe 's/pmbc/PBMC/ig; s/dcs/dendritic_cells/ig; s/THP-1/THP1/ig')\n    target="$dir/$newbase"\n    if [ "$f" != "$target" ]; then\n      if [ -e "$target" ]; then\n        echo "SKIP exists: $target"\n      else\n        mv -v -- "$f" "$target"\n      fi\n    fi\n  done
10355  git add . # manually checked each file that was renamed prior to adding all
10357  gcsmg "rename files for consistency THP-1 -> THP1"
10358  gcmsg "rename files for consistency THP-1 -> THP1"

# content editing logic here
10360  find . -name '*.tsv' -print0 | xargs -0 perl -pi -e 's/pmbc/PBMC/ig; s/dcs/dendritic_cells/ig; s/THP-1/THP1/ig'

# files that were unintentionally changed:
10362  git restore RNAseq_data_forDE/rawexpression/
10364  git restore DE_results/RNAseq/*.tsv
10366  git restore DE_results/RNAseq/up/*.tsv
10367  git restore DE_results/RNAseq/dn/*.tsv
10369  git restore metadata/Homo_sapiens.gene_info.tsv
10371  git restore RNAseq_data_forDE/metadata/

# checking each file individually that was changed prior to adding all
10373  git add DE_results/
10375  gcmsg "rename 'THP-1' -> 'THP1'"
10377  git add signatures/
10379  gcmsg "rename 'THP-1' -> 'THP1' in /signatures"
10381  git add RNAseq_data_forDE/
10382  git add labeling/TBrnaseq_man.tsv
10384  gcmsg "rename 'THP-1' -> 'THP1' in RNA metadata"
10385  git add .
10386  gcmsg "replace DCs -> Dendritic Cells & pmbc -> PBMC in micro metadata"
10388  git log
```
