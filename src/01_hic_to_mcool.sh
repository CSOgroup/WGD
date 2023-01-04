#!/bin/bash


HIC_PATH="data/hic/maps"

for hic_file in `find ${HIC_PATH} -name "*.hic"`
do
	mcool_file="$(echo ${hic_file} | rev | cut -d'.' -f 2- | rev).mcool"
	if [[ ! -f ${mcool_file} ]]
	then
		echo "Converting $(basename ${hic_file}) to .mcool format"
		hic2cool convert ${hic_file} ${mcool_file} -r 0
	fi
done
