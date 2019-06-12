#!/usr/bin/env python

####Description: Concatenate each genome to generate alignment windows using a *.list file
####List file: The absolute path for the five genomes in order, P1 - first genome, P2 - second genome...
####Use: createDfoil_fasta_final.py DATASETXXX.list
####Written by: Henrique V Figueiro - henriquevf@gmail.com

from Bio import SeqIO
import sys
import os

#Open dataset.list file

dataset = open(sys.argv[1], 'r')
lines=dataset.readlines()

dataset_name = sys.argv[1].split('.')[0]

#Create dataset and windows folders

if not os.path.exists(os.getcwd()+'/'+dataset_name):
	os.makedirs(os.getcwd()+'/'+dataset_name)

if not os.path.exists(os.getcwd()+'/'+dataset_name+'/windows'):
	os.makedirs(os.getcwd()+'/'+dataset_name+'/windows')

chr_dir = os.getcwd()+'/'+dataset_name+'/windows/'

#Read genome filenames and get spp names

P1 = str(lines[0].rstrip())
P1_name = str(P1.split('/')[-1].split('_')[0])

P2 = str(lines[1].rstrip())
P2_name = str(P2.split('/')[-1].split('_')[0])

P3 = str(lines[2].rstrip())
P3_name = str(P3.split('/')[-1].split('_')[0])

P4 = str(lines[3].rstrip())
P4_name = str(P4.split('/')[-1].split('_')[0])

P5 = str(lines[4].rstrip())
P5_name = str(P5.split('/')[-1].split('_')[0])

#Create each window fasta file and append the subsequent genomes

with open(P1, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        with open(chr_dir + record.id + ".fa", 'w') as out_file:
            seq_name = str('>'+ P1_name)
            seq =  str(record.seq)
            out_file.write(seq_name + '\n' + seq + '\n')
        out_file.close()

        
with open(P2, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        with open(chr_dir + record.id + ".fa", 'a') as out_file:
            seq_name = str('>'+P2_name)
            seq =  str(record.seq)
            out_file.write(seq_name + '\n' + seq + '\n')
        out_file.close()


with open(P3, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        with open(chr_dir + record.id + ".fa", 'a') as out_file:
            seq_name = str('>'+P3_name)
            seq =  str(record.seq)
            out_file.write(seq_name + '\n' + seq + '\n')
        out_file.close()

with open(P4, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        with open(chr_dir + record.id + ".fa", 'a') as out_file:
            seq_name = str('>'+P4_name)
            seq =  str(record.seq)
            out_file.write(seq_name + '\n' + seq + '\n')
        out_file.close()            


with open(P5, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        with open(chr_dir + record.id + ".fa", 'a') as out_file:
            seq_name = str('>'+P5_name)
            seq =  str(record.seq)
            out_file.write(seq_name + '\n' + seq + '\n')
        out_file.close()

#Generate dfoil commands for the dataset 
#Adapt to your own system
#Run the created bash script to run DFOIL

with open(os.getcwd()+'/'+dataset_name+'/'+dataset_name+".sh", "a") as bash_file:
	bash_file.write('#fasta2dfoil command:'+'\n')
	bash_file.write('for i in windows/*.fa; do python ~/programs/dfoil-master/fasta2dfoil.py $i --out $i.count --names '+P1_name+','+P2_name+','+P3_name+','+P4_name+','+P5_name+'; done'+'\n')
	bash_file.write('#concatenate count files:'+'\n')
	bash_file.write('cat windows/*.count > '+dataset_name+'.count'+'\n')
	bash_file.write('#dfoil command:'+'\n')
	bash_file.write('python ~/programs/dfoil-master/dfoil.py --infile '+dataset_name+'.count --out '+dataset_name+'.stats'+'\n')
	bash_file.write('#dfoil_analyze command:'+'\n')
	bash_file.write('python ~/programs/dfoil-master/dfoil_analyze.py '+dataset_name+'.stats > '+dataset_name+'_results.txt'+'\n')

	






