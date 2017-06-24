IFS=$'\n'

for genus in `cat starget_genus.list`
do
	totalFilepath=/data/mitsuki/data/blast/db/genus_${genus}.fasta
	dbFilepath=/data/mitsuki/data/mmseqs/db/genus_${genus}

	if [ -e ${dbFilepath} ]; 
	then
		echo "EXISTS: "${dbFilepath}
		echo -e "\tdelete existing database first"
		exit 1
	fi

	if [ -e ${totalFilepath} ]; 
	then
		mmseqs createdb ${totalFilepath} ${dbFilepath} 
	else
		echo "ERROR: file not exists "${totalFilepath}
		echo "\texecute create_database.sh first"
		exit 1
	fi
done
