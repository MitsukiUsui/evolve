IFS=$'\n'

genusNum=${1}
catalogFilepath=/home/mitsuki/altorf/evolve/createdatabase/out/genus_${genusNum}.txt
mmseqsDir=/data/mitsuki/data/mmseqs
dbFilepath=${mmseqsDir}/db/genus_${genusNum}

for line in `tail -n +2 ${catalogFilepath}` 
do
	basename=`echo ${line}|cut -d, -f 12`
	queryFilepath=${mmseqsDir}/query/${basename}
	resultFilepath=${mmseqsDir}/result/${basename}
	m8Filepath=${mmseqsDir}/result/${basename}.m8

	#delete old if exists
	if [ -e ${resultFilepath} ]; 
	then
		echo "EXISTS: "${resultFilepath}
		echo -e "\tdelete existing database first"
		exit 1
	fi

	#identification
	mmseqs search ${queryFilepath} ${dbFilepath} ${resultFilepath}\
               /tmp/mitsuki -a -e 1 --threads 8 --remove-tmp-files 
	mmseqs convertalis ${queryFilepath} ${dbFilepath} ${resultFilepath} ${m8Filepath}\
               	--threads 8 --format-mode 1
done
