# lyase_prediction
## Install the software
```
# install env
conda env create -f lyase_env.yml

# activate env
source activate lyase_env

# install perl 5.22
conda install -c anaconda perl=5.22
```

## prediction
```
# activate env
source activate lyase_env

# run
python phage_lyase.py -p /.../input_path/ -t Bacteria -hd /db/hmm/lysin_reported.hmm -cd /db/cazy/db/ -rl /db/hmm/lysin_reported.txt -wkdir /test/
```
