IFS=$'\n'

for genus in `cat starget_genus.list`
do
	outFilepath=/data/mitsuki/data/blast/db/genus_${genus}.fasta

	txtFilepath=../out/genus_${genus}.txt
	for line in `tail -n +2 ${txtFilepath}`
	do
		basename=`echo ${line}|cut -d, -f 12`
		dir=/data/mitsuki/data/refseq
		filepath=${dir}/protein/${basename}_protein.faa
		cat ${filepath}>>outFilepath
	done

	makeblastdb -in ${outFilepath} -dbtype prot -out genus_${genus}
	#rm ${outFilepath}
	echo "DONE ${outFilepath}."
done
