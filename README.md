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

### Change suffix
***Notice:*** **Genome fasta file suffix is ```.fna```, to see example in Data fold**

***If Genome fasta file suffix is not ```.fna```, you can run :***
```
# Rename
python rename_suffix.py -p ./data/ -s fasta

usage: rename_suffix.py [-h] -p PATH -s SUFFIX
suffix rename
options:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  genome sequence path
  -s SUFFIX, --suffix SUFFIX
                        old suffix/To be modified suffix

```
For example, if fasta file's suffix is '.fasta', run ```python rename_suffix.py -p ./data/ -s fasta```, that can change suffix '.fasta' into '.fna'

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
  -bp B                                             # 'B' for bacteria, 'P' for phage"
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
