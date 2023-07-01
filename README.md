# Phatides_prediction
## Install the software
```
# install env
conda env create -f phatides_prediction_env.yml

# activate env
source activate phatides_prediction_env
```

If your perl version is not 5.22, please install perl=5.22
```
# install perl 5.22
conda install -c anaconda perl=5.22
```

## Prediction
```
# activate env
source activate phatides_prediction_env

# run
python phatides_prediction.py
  -p /.../input_path/                               # genome sequnce path
  -t Bacteria                                       # prokka kingdom type    
  -hd ./db/hmm/lysin_reported.hmm                   # hmmer database path
  -rl ./db/hmm/lysin_reported.txt                   # reported lysin structures(hmm files)
  -cd ./db/cazy/db/                                 # cazy database path
  -wkdir ./test/                                    # work directory
  -ml 10000                                         # lower proteins molecular weight
  -mu 40000                                         # upper proteins molecular weight
```

## Bactericidal activity scoring 
### 1. Get different peptides
```
python peptides_split.py -i putative_PGHs.fa -o putative_PGHs.fa -mi 6 -ma 50
```

### 2. Sequence to vector
```
perl format.pl putative_PGHs.fa none > putative_PGHs.txt
```

### 3. Scoring
```
python prediction_lstm.py putative_PGHs.txt putative_PGHs.out
```

**If you use prediction_lstm.py please cite [Identification of antimicrobial peptides from the human gut microbiome using deep learning](https://www.nature.com/articles/s41587-022-01226-0)**
