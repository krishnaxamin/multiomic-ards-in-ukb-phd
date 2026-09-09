#!bin/bash

# LAVA - download UK Biobank LD reference and unzip
# check for presence before proceeding
cd ~/data/external/genomics/lava/ld_ref
ld_file_count=$(find . -type f -name "*.bcor" | wc -l)

if (( ld_file_count != 23 ));
then
	echo "Downloading LAVA LD reference files."
	wget https://vu.data.surfsara.nl/index.php/s/7NBVIvtPRdu7Qhz/download -O lava-ukb-v1.1_chr1-2.zip
	wget https://vu.data.surfsara.nl/index.php/s/fy6ITboMojHrQXr/download -O lava-ukb-v1.1_chr3-4.zip
	wget https://vu.data.surfsara.nl/index.php/s/mRz31q0lq7KMcuI/download -O lava-ukb-v1.1_chr5-6.zip
	wget https://vu.data.surfsara.nl/index.php/s/89RXxPN2BlOqxLs/download -O lava-ukb-v1.1_chr7-9.zip
	wget https://vu.data.surfsara.nl/index.php/s/3dU1L2Hap43xuCs/download -O lava-ukb-v1.1_chr10-12.zip
	wget https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download -O lava-ukb-v1.1_chr13-16.zip
	wget https://vu.data.surfsara.nl/index.php/s/U8eg5XfTr8qPzPp/download -O lava-ukb-v1.1_chr17-23.zip

	for z in *.zip; do
		unzip -n "$z" && rm "$z"
	done
else
	echo "LAVA LD reference files already downloaded."
fi