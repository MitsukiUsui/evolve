IFS=$'\n'

for genus in `cat starget_genus.list`
do
	dbDir=/data/mitsuki/data/blast/db
	totalFilepath=${dbDir}/genus_${genus}.fasta
	totalSeqs=0
	dbName=${dbDir}/genus_${genus}


	#delete old if exists
	if [ -e ${totalFilepath} ]; 
	then
		rm ${totalFilepath} 
		echo -e "DELETED old "${totalFilepath}"\n"
	fi

	#concatination
	echo "CONCAT to ${totalFilepath}"
	
	catalogFilepath=../out/genus_${genus}.txt
	for line in `tail -n +2 ${catalogFilepath}`
	do
		basename=`echo ${line}|cut -d, -f 12`
		filepath=/data/mitsuki/data/refseq/protein/${basename}_protein.faa
		numHeader=`grep -c "^>" ${filepath}`
		totalSeqs=`expr ${totalSeqs} + ${numHeader}`
		cat ${filepath} >> ${totalFilepath}
		echo -e "\t${filepath}: ${numHeader} seqs"
	done

	#check whether correctly updated
	numHeader=`grep -c "^>" ${totalFilepath}`
	if [ ${totalSeqs} = ${numHeader} ]
	then
		echo "CREATED ${totalFilepath}: ${numHeader} seqs"
	else
		echo "ERROR: the number of sequence dose not match ("${totalSeqs}", "${numHeader}")"
		exit 1
	fi

	#makeblastdb
	makeblastdb -in ${totalFilepath} -dbtype prot -out ${dbName}
	echo "CREATED ${dbName}"
done
