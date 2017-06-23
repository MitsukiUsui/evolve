IFS=$'\n'
genusNum=${1}

catalogFilepath=/home/mitsuki/altorf/evolve/createdatabase/out/genus_${genusNum}.txt
dbName=/data/mitsuki/data/blast/db/genus_${genusNum}

for line in `tail -n +2 ${catalogFilepath}` 
do
	basename=`echo ${line}|cut -d, -f 12`
	dir=/data/mitsuki/out/altorf/evolve
	queryFilepath=${dir}/query/${basename}.query
	outFilepath=${dir}/result/${basename}.xml

	#delete old if exists
	if [ -e ${outFilepath} ]; 
	then
		rm ${outFilepath} 
		echo -e "DELETED old "${outFilepath}"\n"
	fi

	blastp -db ${dbName}\
		   -query ${queryFilepath}\
		   -out ${outFilepath}\
		   -outfmt 5\
		   -evalue 1 &

	echo ${queryFilepath}
	echo ${outFilepath}
done
