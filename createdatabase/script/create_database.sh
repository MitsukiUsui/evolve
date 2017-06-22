IFS=$'\n'

for genus in `cat starget_genus.list`
do
	catalogFilepath=../out/genus_${genus}.txt
	dbDir=/data/mitsuki/data/blast/db
	fastaFilepath=${dbDir}/genus_${genus}.fasta
	dbName=${dbDir}/genus_${genus}
	
	rm ${fastaFilepath}#!TBI! check existance before
	echo "CONCAT to ${fastaFilepath}"
	for line in `tail -n +2 ${catalogFilepath}`
	do
		basename=`echo ${line}|cut -d, -f 12`
		filepath=/data/mitsuki/data/refseq/protein/${basename}_protein.faa
		echo -e "\t${filepath}"
		cat ${filepath} >> ${fastaFilepath}
	done

	makeblastdb -in ${fastaFilepath} -dbtype prot -out ${dbName}
	echo "CREATED ${dbName}"
done
