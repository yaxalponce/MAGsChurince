# Rename Bordetella

cd /mnt/data/yaxal/CC/DASTool/Bordetella/Bordetella_DASTool_bins

for file in $(ls *.fa); do
	NewName=$(sed -e 's/Bin\./Churince.Bin_/g' <<< $file )
	Bin=$(sed 's/\.fa//g' <<< $NewName)
	echo "$Bin"
	awk -v bin="$Bin" '/^>/{$0=">"bin"_"(++i)}1' $file > $NewName
	rm $file
done


# Rename Mats

cd /mnt/data/yaxal/CC/DASTool/Mats/Mats_DASTool_bins

for file in $(ls *.fa); do
	NewName=$(sed -e 's/\(concoct\|metabat\)/Mats/g' -e 's/_sub//g' -e 's/\.0\+/\.Bin_/g' <<< $file )
	Bin=$(sed 's/\.fa//g' <<< $NewName)
	echo "$Bin"
	awk -v bin="$Bin" '/^>/{$0=">"bin"_"(++i)}1' $file > $NewName
	rm $file
done


#Rename Sediments

cd /mnt/data/yaxal/CC/DASTool/Sediment/Sediment_DASTool_bins

for file in $(ls *.fa); do
	NewName=$(sed -e 's/\(concoct\|metabat\)/Sediment/g' -e 's/_sub//g' -e 's/\.0\+/\.Bin_/g' <<< $file )
	Bin=$(sed 's/\.fa//g' <<< $NewName)
	echo "$Bin"
	awk -v bin="$Bin" '/^>/{$0=">"bin"_"(++i)}1' $file > $NewName
	rm $file
done


# Rename Water

cd /mnt/data/yaxal/CC/DASTool/Water/Water_DASTool_bins

for file in $(ls *.fa); do
	NewName=$(sed -e 's/\(concoct\|metabat\)/Water/g' -e 's/_sub//g' -e 's/\.0\+/\.Bin_/g' <<< $file )
	Bin=$(sed 's/\.fa//g' <<< $NewName)
	echo "$Bin"
	awk -v bin="$Bin" '/^>/{$0=">"bin"_"(++i)}1' $file > $NewName
	rm $file
done
