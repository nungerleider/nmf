from Bio import SeqIO
import sys
import os
import glob

''' 

Takes a vcf file (chromosome, position, id, ref, alt, etc in that order) and directory of fastas as input. Prints vcf file with added column for trinucleotide context to standard out.


'''

def parse():

    script_name = os.path.basename(sys.argv.pop(0))

    helpmsg = "\n\tUsage: {script} vcf_file genome_directory\n".format(script=script_name) 
    helpmsg += "genome_directory must contain individual fasta files for each chromosome."
    helpmsg += " i.e. 'chr1.fa', 'chr2.fa', etc. Or '1.fa', '2.fa'; as long as the text up until the '.fa' matches"
    helpmsg += " the chromosome names in the vcf file."

    if len(sys.argv) < 2:
        print("\nYou need to supply at least two arguments..a vcf file and a directory (containing genome fasta files.)\n")
        sys.exit(helpmsg)
    
    if sys.argv[1] in ['-h', '--help']:
        sys.exit(helpmsg)

    vcf = sys.argv.pop(0)
    if not os.path.exists(vcf):
        sys.exit("VCF file not found: %s" % vcf) 

    genome_dir = sys.argv.pop(0)

    fastas = glob.glob(os.path.join(genome_dir,'*fa'))
    if not fastas:
        sys.exit("No fasta files found in your genome directory: %s" % genome_dir) 

    return vcf, genome_dir


def process_vcf(vcf, genome_dir):

    c_init = ''
    with open(vcf, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
            
            if line[0] == '#':
                print(line)
                continue

            chromosome, position, _, ref, alt, *_ = line.split('\t')
            position = int(position)
            if chromosome != c_init:
                c_init = chromosome
                try:
                    loaded_fa = SeqIO.read(os.path.join(genome_dir, c_init + '.fa'), 'fasta')
                except: 
                    sys.exit("%s not found or not in the proper format" % c_init)

            subseq = str(loaded_fa[position - 2: position + 1].seq).upper()

            line += '\t%s' % subseq
            print(line)


vcf, genome_dir = parse()
process_vcf(vcf, genome_dir)
