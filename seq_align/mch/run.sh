nucmer --maxmatch $1 $2 -p $3
mummerplot "${3}.delta" --png --filter -p "${3}_plot"
show-coords -Hrlc "${3}.delta" | show-snps -CIS "${3}.delta" | sed "s/|//g" | sed "s/ /\t/g" | tr -s "\t" | python3 plot_snps.py -o "${3}"
