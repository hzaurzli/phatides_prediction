import os
import argparse


def load_fa(path):
    """a function to read fasta file from the path and store in a dict"""
    genes_seq = {}  #将序列存入字典
    with open(path,"r") as sequences:  #以读取方式打开文件
        lines = sequences.readlines()

    for line in lines:
        if line.startswith(">"):
            genename = line.split()[1]  #这个地方需要灵活调整
            genes_seq[genename] = ''  #序列为字符串
        else:
            genes_seq[genename] += line.strip()

    return genes_seq


def build_kmers(seq, k_size,shift):
    """a function to calculate kmers from seq"""
    kmers = []  # k-mer存储在列表中
    n_kmers = len(seq) - k_size + 1

    for i in range(0,n_kmers,int(shift)+1):
        kmer = seq[i:i + k_size]
        kmers.append(kmer)

    return kmers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="peptides")
    parser.add_argument("-i", "--input", required=True, type=str, help="Input fasta")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output fasta")
    parser.add_argument("-mi", "--min_length", required=False, default=6, type=int, 
                            help="Minimum peptides length")
    parser.add_argument("-ma", "--max_length", required=False, default=50, type=int, 
                            help="Maximum peptides length")
    Args = parser.parse_args()


    fa = load_fa(Args.input)
    with open(Args.output,'w') as w:
        for j in range(int(Args.min_length),int(Args.max_length)):
            for key in fa:
                fa_b = build_kmers(fa[key],j,0)
                count = 1
                for i in fa_b:
                    f = key + '-' + str(count) + '-' + 'length:' + str(j) + '\n'
                    s = i + '\n'
                    line = f + s
                    w.write(line)
                    count += 1
    w.close()



