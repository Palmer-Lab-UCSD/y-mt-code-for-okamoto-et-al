for FILE in data/phenotypes/*.txt
do
	BASE=`basename $FILE .txt`
    	gcta --mlma --bfile results/groups/my_haplotypes --pheno $FILE \
			--grm data/genotypes/autosomes --out results/phenotypes/$BASE
done

