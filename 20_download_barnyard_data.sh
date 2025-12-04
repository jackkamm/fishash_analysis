mkdir -p outs/barnyard_data/raw
cd outs/barnyard_data/raw
cat ../../../data/barnyard_data_urls.txt | xargs -I% wget %
