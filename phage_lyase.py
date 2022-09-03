import argparse
import os,sys,re
import subprocess as sub
from subprocess import *
import glob
import shutil
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class tools:
    def __init__(self):
        self.prokka = 'prokka'
        self.phispy = 'PhiSpy.py'
        self.phanotate = 'phanotate.py'
        self.cdHit = '/home/runzeli/software/cdhit/cd-hit/cd-hit'
        self.rundbcan = 'run_dbcan.py'
        self.hmmsearch = 'hmmsearch'
        self.tmhmm = '/home/runzeli/software/tmhmm/tmhmm-2.0c/bin/tmhmm'


    def run(self, cmd, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p.returncode

    def run_prokka(self, fastain, fastaout, prefix, type_annotation):
        cmd = '%s %s -o %s --prefix %s --kingdom %s' % (self.prokka, fastain, fastaout,prefix,type_annotation)
        return cmd

    def run_phispy(self, gbk_input, out, profix, phage_genes):
        cmd = '%s %s -o %s -p %s --threads 8 --phage_genes %s' % (self.phispy, gbk_input, out, profix, phage_genes)
        return cmd

    def run_phanotate(self, inputfile, out):
        cmd = '%s %s -o %s' % (self.phanotate, inputfile, out)
        return cmd

    def run_cdhit(self,inputfile, out, cutoff):
        cmd = '%s -i %s -o %s -c %s' % (self.cdHit, inputfile, out, cutoff)
        return cmd

    def scan_dbscan(self,inputfile, out, db):
        cmd = '%s %s protein -t hmmer --out_dir %s --db_dir %s' % (self.rundbcan, inputfile, out, db)
        return cmd

    def run_hmmsearch(self,tblout, e_val,hmm ,inputfile):
        cmd = '%s --tblout %s -E %s --cpu 2 %s %s' % (self.hmmsearch, tblout, e_val, hmm,inputfile)
        return cmd

    def run_tmhmm(self,input_fa,out):
        cmd = '%s %s > %s' % (self.tmhmm,input_fa,out)
        return cmd


def prophage_select(prophage_input,fna_input,ppn_out):
    info_file = open(prophage_input, "r")
    info_list = []
    for i in info_file:
        data = i.strip().split('\t')
        pp_num = data[0]
        contig = data[1]
        info_s = data[2]
        info_e = data[3]
        info_list.append((pp_num, contig, info_s, info_e))

    for record in SeqIO.parse(fna_input, 'fasta'):
        Contig_ID = record.id
        Desp = record.description.split(",")[0]
        for i in info_list:
            pp_num = i[0]
            contig = i[1]
            info_s = i[2]
            info_e = i[3]
            if contig == Contig_ID:
                file_name_use = ppn_out
                out_file = open(file_name_use + "_" + str(pp_num) + ".fasta", "a")
                gene_seq = record.seq[int(info_s) - 1: int(info_e)]
                element_record = SeqRecord(gene_seq, id=pp_num, description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
                out_file.close()
    info_file.close()


def Gene_element_abstract(ppn_out,ppn_fa,ppn_ffn):
    n = 0
    gene_list = []

    list_path = str(ppn_fa).split('/')
    list_path = [i for i in list_path if i != '']
    if len(list_path) == 1:
        file_name = ppn_fa
    else:
        file_name = os.path.basename(ppn_fa)

    ppn_out = open(ppn_out,'r')
    for line in ppn_out:
        if line.startswith('#'):
            pass
        else:
            line_info = line.strip().split("\t")
            Contig_ID = file_name.strip().split(".")[0]
            location_S = int(line_info[0])  # start
            location_E = int(line_info[1])  # end
            if line_info[2] == "-":
                F_R = "R"  # reverse
                location_S = int(line_info[1])
                location_E = int(line_info[0])
            if line_info[2] == "+":
                F_R = "F"
                location_S = int(line_info[0])
                location_E = int(line_info[1])
            gene_list.append((location_S, location_E, F_R))
    ppn_out.close()

    out_file = open(ppn_ffn, "a")
    for record in SeqIO.parse(ppn_fa, "fasta"):
        ID_contig = record.id
        Desp = record.description
        for info in gene_list:
            file_name_1 = file_name
            file_name_2 = file_name_1.split(".")[0]
            if info[2] == "F":
                gene_seq = record.seq[info[0] - 1:info[1]]
                gene_protein = gene_seq.translate()
            elif info[2] == "R":
                gene_seq_ori = record.seq[info[0] - 1:info[1]]
                gene_seq = gene_seq_ori.reverse_complement()
                gene_protein = gene_seq.translate()
            element_ID = file_name_2 + ":" + info[2] + ":" + str(info[0]) + "-" + str(info[1])
            element_record = SeqRecord(gene_protein[0:-1],id=element_ID, description=Desp)
            SeqIO.write(element_record, out_file, "fasta")
    out_file.close()


def molecular_weight(protein_fa,protein_filter_fa):
    protein_fa_info = open(protein_fa, "r")
    out_file = open(protein_filter_fa, "a")
    molecular_weight = open("./molecular_weight.txt", "w")
    for record in SeqIO.parse(protein_fa_info, "fasta"):
        ID_contig = record.id
        Seq_use = record.seq
        Desp = record.description
        protein_seq = str(Seq_use)
        X = ProteinAnalysis(protein_seq)
        MW_cal = "%0.2f" % X.molecular_weight()

        if 'X' not in protein_seq or '*' not in protein_seq[:-1] and float(MW_cal) <= 40000:
            element_record = SeqRecord(Seq_use, id=ID_contig, description=Desp)
            SeqIO.write(element_record, out_file, "fasta")
            molecular_weight.write(ID_contig + "\t" + MW_cal + "\n")

    protein_fa_info.close()
    molecular_weight.close()

def find_cazyme(cdhit_filter,cazy_overview):
    input_file = open(cdhit_filter, "r")
    info_file = open(cazy_overview, "r")
    info_list = []
    for i in info_file:
        data = i.strip().split('\t')
        gen_id = data[0]
        info_list.append(gen_id)

    for record in SeqIO.parse(input_file, 'fasta'):
        Contig_ID = record.id
        Desp = record.description
        for i in info_list:
            gen_id = i
            if gen_id == Contig_ID:
                gene_seq = record.seq
                out_file = open("./all_protein_filter_cazyme.fasta", "a")
                element_record = SeqRecord(gene_seq, id='', description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
                out_file.close()
    input_file.close()
    info_file.close()


def find_pfam(cdhit_filter,lyase_list):
    input_file = open(cdhit_filter, "r")
    info_file = open('./hmmer_out/all_protein_filter_hmmer_out.txt', "r")
    info_file_two = open(lyase_list, "r")

    f2 = open(r'./hmmer_out/all_protein_filter_hmmer_out.txt', 'w')
    for i in info_file:
        line = re.split('\s+', i)
        new_line = ' '.join(line)
        new_line = new_line.strip(' ')
        f2.write(new_line)
        f2.write("\n")
    info_file.close()
    f2.close()

    info_pfam_list = []
    for i in info_file_two:
        data2 = i.strip()
        pfam_num_reported = str(data2)
        info_pfam_list.append(pfam_num_reported)

    info_file_new = open(r'./hmmer_out/all_protein_filter_hmmer_out.txt', 'r')
    info_list = []
    for j in info_file_new:
        data1 = j.strip().split(' ')
        gen_id = str(data1[0])
        pfam_num = str(data1[3])
        for k in info_pfam_list:
            if k in pfam_num:
                info_list.append(gen_id)

    out_file = open("./all_protein_pfam_protein.fasta", "a")
    for record in SeqIO.parse(input_file, 'fasta'):
        Contig_ID = record.id
        Desp = record.description
        for i in info_list:
            gen_id = i
            if gen_id == Contig_ID:
                gene_seq = record.seq
                element_record = SeqRecord(gene_seq, id='', description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
    input_file.close()
    info_file_new.close()
    info_file_two.close()
    out_file.close()

def detete_TMhelix(cdhit_fasta,cazyme_pfam_TMhelix):
    input_file = open(cdhit_fasta, "r")
    info_file = open(cazyme_pfam_TMhelix, "r")

    lines = info_file.readlines()
    info_short_file = open(r"./all_protein_final_tmhmm_shortout.txt", 'w')
    content = "#"
    for line in lines:
        if line.strip()[0] != content:
            info_short_file.write(line)
    info_file.close()
    info_short_file.close()

    info_short_file = open(r"./all_protein_final_tmhmm_shortout.txt", 'r')
    info_list = []
    for i in info_short_file:
        data = i.strip().split('\t')
        gen_id = str(data[0])
        if gen_id not in info_list:
            info_list.append(gen_id)
    info_short_file.close()

    info_short_file = open(r"./all_protein_final_tmhmm_shortout.txt", 'r')
    dele_list = []
    for j in info_short_file:
        data_dele = j.strip().split('\t')
        gen_id_dele = str(data_dele[0])
        protein_type = str(data_dele[2])
        dele_list.append((gen_id_dele, protein_type))
    info_short_file.close()
    for k in dele_list:
        gen_id_dele_last = str(k[0])
        protein_type_dele = str(k[1])
        if (protein_type_dele == "TMhelix") and (gen_id_dele_last in info_list):
            info_list.remove(gen_id_dele_last)

    out_file = open("putative_lyase.fa", "a")
    for record in SeqIO.parse(input_file, 'fasta'):
        Contig_ID = record.id
        Desp = record.description
        for i in info_list:
            gen_id = i
            if gen_id == Contig_ID:
                gene_seq = record.seq
                element_record = SeqRecord(gene_seq, id='', description=Desp)
                SeqIO.write(element_record, out_file, "fasta")
    input_file.close()
    out_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="lyase prediction")
    parser.add_argument("-p", "--path", required=True, type=str, help="genome sequnce path")
    parser.add_argument("-t", "--type", required=True, type=str, help="prokka kingdom type")
    parser.add_argument("-c", "--cdhit_cutoff", default=0.95,required=False, type=float, help="cdhit cluster cutoff")
    parser.add_argument("-hc", "--hmmer_cutoff", default=1e-5,required=False, type=float, help="hmmer search cutoff")
    parser.add_argument("-hd", "--hmmer_db", required=True, type=str, help="hmmer database path")
    parser.add_argument("-cd", "--cazy_db", required=True, type=str, help="hmmer database path")
    parser.add_argument("-rl", "--reported_lyase", required=True, type=str, help="reported lyase structures(hmm files)")
    parser.add_argument("-wkdir", "--workdir", required=True, type=str, help="Work directory")
    Args = parser.parse_args()

    tl = tools()
    # step 1 prokka annotates ORFs
    curr_dir = sub.getoutput('pwd')
    os.chdir(Args.workdir)
    if os.path.isdir('./prokka_result/') == True:
        pass
    else:
        os.mkdir('./prokka_result/')

    target = Args.path
    if target[-1] == '/':
        target = target
    elif target[-1] != '/':
        target = target + '/'

    if target[0] == '.':
        if target[1] == '/':
            target_suffix = target[1:]
        elif target[1] == '.':
            curr_dir = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            target_suffix = target[2:]
    else:
        target_suffix = target
        curr_dir = ''

    type_annotation = Args.type
    for i in os.listdir(curr_dir + target_suffix):
        name = i.split('.')[0]
        suffix = i.split('.')[1]
        cmd_1 = tl.run_prokka(curr_dir + target_suffix + i,
                          './prokka_result/' + name + '/',name,type_annotation)
        tl.run(cmd_1)

    # step 2 phispy predict prophage
    if os.path.isdir('./phispy_out/') == True:
        pass
    else:
        os.mkdir('./phispy_out/')
    for i in os.listdir('./prokka_result/'):
        cmd_2 = tl.run_phispy('./prokka_result/' + i + '/' + i + '.gbk',
                          './phispy_out/' + i,
                          i,0)
        tl.run(cmd_2)

    # step 3 select 1,2,3,4 colum for coordinate.tsv
    if os.path.isdir('./ppn/') == True:
        pass
    else:
        os.mkdir('./ppn/')

    fna_suffix = os.path.splitext(os.listdir(curr_dir + target_suffix)[0])[-1]
    for i in os.listdir('./phispy_out/'):
        if os.path.isdir('./ppn/' + i) == True:
            pass
        else:
            os.mkdir('./ppn/' + i)
        prophage_select('./phispy_out/'+ i + '/' + i + '_prophage_coordinates.tsv',
                        curr_dir + target_suffix + i + fna_suffix,'./ppn/' + i + '/' + i)


    # step 4 phanotate annotates prophage ORFs
    for i in os.listdir('./phispy_out/'):
        for j in os.listdir('./ppn/' + i):
            j_prefix = j.split('.')[0]
            j_suffix = j.split('.')[1]
            if j_suffix == 'fasta':
                cmd_3 = tl.run_phanotate('./ppn/' + i + '/' + j,
                                         './ppn/' + i + '/' + j_prefix + '.out')

                tl.run(cmd_3)


    if os.path.isdir('./orf_ffn/') == True:
        pass
    else:
        os.mkdir('./orf_ffn/')
    for i in os.listdir('./phispy_out/'):
        for j in os.listdir('./ppn/' + i):
            j_prefix = j.split('.')[0]
            j_suffix = j.split('.')[1]
            if j_suffix == 'fasta':
                Gene_element_abstract('./ppn/' + i + '/' + j_prefix + '.out',
                                      './ppn/' + i + '/' + j,
                                      './orf_ffn/' + j_prefix + '.ffn')


    # step 5 combine prokka faa and ppn faa together
    for i in os.listdir('./prokka_result/'):
        for j in os.listdir('./prokka_result/' + i):
            j_suffix = j.split('.')[1]
            if j_suffix == 'faa':
                os.system('cp %s %s' % ('./prokka_result/' + i + '/' + j,'./orf_ffn/'))

    os.system('cat ./orf_ffn/* > all_protein.faa')

    # step 6 cdhit cluster
    cmd_4 = tl.run_cdhit('./all_protein.faa','./all_protein_cdhit.faa',Args.cdhit_cutoff)
    tl.run(cmd_4)

    # step 7 calculate molecular weight
    molecular_weight('./all_protein_cdhit.faa','./all_protein_cdhit_filter.faa')

    # step 8 scan CAZY database
    cazy_db = Args.cazy_db
    if cazy_db[-1] == '/':
        cazy_db = cazy_db
    elif cazy_db[-1] != '/':
        cazy_db = cazy_db + '/'

    if cazy_db[0] == '.':
        if cazy_db[1] == '/':
            cazy_db_suffix = cazy_db[1:]
        elif cazy_db[1] == '.':
            curr_dir = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            cazy_db_suffix = cazy_db[2:]
    else:
        cazy_db_suffix = cazy_db
        curr_dir = ''

    cmd_5 = tl.scan_dbscan('./all_protein_cdhit_filter.faa','./CAZY_out/',
                           curr_dir + cazy_db_suffix)
    tl.run(cmd_5)

    find_cazyme('./all_protein_cdhit_filter.faa','./CAZY_out/overview.txt')

    # step 9 hmmsearch reported lyase structure in pfam
    tl = tools()
    if os.path.isdir('./hmmer_out/') == True:
        pass
    else:
        os.mkdir('./hmmer_out/')

    hmmer_db = Args.hmmer_db
    if hmmer_db[0] == '.':
        if hmmer_db[1] == '/':
            hmmer_db_suffix = hmmer_db[1:]
        elif hmmer_db[1] == '.':
            curr_dir = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            hmmer_db_suffix = hmmer_db[2:]
    else:
        hmmer_db_suffix = hmmer_db
        curr_dir = ''

    cmd_6 = tl.run_hmmsearch('./hmmer_out/all_protein_filter_hmmer_out.txt', Args.hmmer_cutoff,
                             curr_dir + hmmer_db_suffix,
                             './all_protein_cdhit_filter.faa')
    tl.run(cmd_6)

    reported_lyase = Args.reported_lyase
    if reported_lyase[0] == '.':
        if reported_lyase[1] == '/':
            reported_lyase_suffix = reported_lyase[1:]
        elif reported_lyase[1] == '.':
            curr_dir = os.path.abspath(os.path.join(os.path.dirname(curr_dir + '/'), os.path.pardir))
            reported_lyase_suffix = reported_lyase[2:]
    else:
        reported_lyase_suffix = reported_lyase
        curr_dir = ''

    find_pfam('./all_protein_cdhit_filter.faa', curr_dir + reported_lyase_suffix)

    # step 10 combine results of CAZY and pfam and cdhit clustering
    os.system('cat all_protein_filter_cazyme.fasta all_protein_pfam_protein.fasta > cazyme_pfam.fasta')
    cmd_7 = tl.run_cdhit('./cazyme_pfam.fasta', './cazyme_pfam_cdhit.fasta', Args.cdhit_cutoff)
    tl.run(cmd_7)

    # step 11 remove TMhelix
    cmd_8 = tl.run_tmhmm('./cazyme_pfam_cdhit.fasta','./cazyme_pfam_TMhelix.out')
    tl.run(cmd_8)
    detete_TMhelix('./cazyme_pfam_cdhit.fasta','./cazyme_pfam_TMhelix.out')

    os.system('rm -r ./CAZY_out/ ./hmmer_out/ ./orf_ffn/ ./phispy_out/ ./ppn/ ./prokka_result/ ./TMHMM*/')
    os.remove('./all_protein_cdhit.faa')
    os.remove('./all_protein_cdhit.faa.clstr')
    os.remove('./all_protein_cdhit_filter.faa')
    os.remove('./all_protein.faa')
    os.remove('./all_protein_filter_cazyme.fasta')
    os.remove('./all_protein_final_tmhmm_shortout.txt')
    os.remove('./all_protein_pfam_protein.fasta')
    os.remove('./cazyme_pfam_cdhit.fasta')
    os.remove('./cazyme_pfam_cdhit.fasta.clstr')
    os.remove('./cazyme_pfam.fasta')
    os.remove('./cazyme_pfam_TMhelix.out')







