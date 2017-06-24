IFS=$'\n'
genus=872

catalogFilepath=../../createdatabase/out/genus_${genus}.txt
for line in `tail -n +2 ${catalogFilepath}`
do
	basename=`echo ${line}|cut -d, -f 12`
	queryFilepath=/data/mitsuki/out/altorf/evolve/query/${basename}.query
	dbFilepath=/data/mitsuki/data/mmseqs/query/${basename}

	if [ -e ${dbFilepath} ]; 
	then
		echo "EXISTS: "${dbFilepath}
		echo -e "\tdelete existing database first"
		exit 1
	fi
	
	#check existance
	if [ -e ${queryFilepath} ]; 
	then
		mmseqs createdb ${queryFilepath} ${dbFilepath} 
		echo "DONE: "${dbFilepath}
	else
		echo "ERROR: file not exists "${queryFilepath}
		echo -e "\texecute create_query.py first"
		exit 1
	fi
done
