IFS=$'\n'

for genus in `cat target_genus.list`
do
	txtFilepath=../out/genus_${genus}.txt
	for line in `tail -n +2 ${txtFilepath}`
	do
		taxid=`echo ${line}|cut -d, -f 1`
		filepath=~/altorf/genome/speciespick/picked_assembly_summary.csv
		basename=`awk -v "t=${taxid}" -F, '$7==t {print $5"/"$4}' ${filepath}`
		cdsFtppath=${basename}_cds_from_genomic.fna.gz
		proFtppath=${basename}_protein.faa.gz

		dir=/data/mitsuki/data/refseq
		wget -P ${dir}/cds_from_genomic ${cdsFtppath}
		wget -P ${dir}/protein ${proFtppath}
		
		echo ${cdsFtppath}
		echo ${proFtppath}
	done
done
