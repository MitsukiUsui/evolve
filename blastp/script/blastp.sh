IFS=$'\n'
dir=/data/mitsuki/out/altorf/evolve

for line in `tail -n +2 ../out/target_genus.list`
do
	basename=`echo ${line}|cut -d, -f 12`
	queryFilepath=${dir}/query/${basename}.query
	outFilepath=${dir}/result/${basename}.xml.new

	blastp -db /data/mitsuki/data/blast/db/genus_872\
		   -query ${queryFilepath}\
		   -out ${outFilepath}\
		   -outfmt 5\
		   -evalue 1 &

	echo ${queryFilepath}
	echo ${outFilepath}
done
