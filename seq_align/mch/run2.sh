nucmer --maxmatch $1 $2 -p $3
mummerplot "${3}.delta" --png --filter -p "${3}_plot"
show-coords -Trlc "${3}.delta" > "${3}.filtered"
show-coords -Hrlc "${3}.delta" | show-snps -CIS "${3}.delta" | sed "s/|//g" | sed "s/ /\t/g" | tr -s "\t" > "${3}.snps"
