#!/usr/bin/env python3

import logging
from tqdm import tqdm
import sys, os, time, gzip, argparse, re
from io import TextIOWrapper

from multiprocessing import Pool

from subprocess import Popen, PIPE

import tempfile
import ast

class gene:
    def __init__(self, ID, contig, start, stop, strand, geneType, product = None, inference = None):
        self.ID = ID
        self.contig = contig
        self.start = int(start)
        self.stop = int(stop)
        self.type = geneType
        self.strand = strand    
        self.inference = inference
        self.product = product

    def get_gff(self):
        ## soft version ...?
        Feature = f"{self.contig}\tSOFT\tgene\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID=gene_{self.ID}\n"
        Feature += f'{self.contig}\tSOFT\t{self.type}\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID={self.ID};Parent=gene_{self.ID};inference={self.inference}'
        Feature += f";product={self.product}\n" if self.product else "\n"
        return Feature
    
    def get_faa(self, dna_seq):
        faa_seq = translate(dna_seq)
        seq = f">{self.ID}\n"
        j = 0
        seq+= faa_seq[j:j+60] + "\n"
        j = 60
        while j < len(faa_seq):
            seq+= faa_seq[j:j+60] + "\n"
            j+=60
        return seq
    
    def get_ffn(self, dna_seq):
        seq = f">{self.ID}\n"
        j = 0
        seq+= dna_seq[j:j+60] + "\n"
        j = 60
        while j < len(dna_seq):
            seq+= dna_seq[j:j+60] + "\n"
            j+=60
        return seq
            

def reverse_complement(seq): 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N' }
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq

def translate(seq): 
    
    start_table = { 
        'ATA':'M', 'ATC':'M', 'ATT':'M', 'ATG':'M', 'ATN':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'ACN':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'M', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CCN':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CGN':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'M', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GCN':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GGN':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TCN':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'M', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 

    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'ACN':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CTN':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CCN':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CGN':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GTN':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GCN':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GGN':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TCN':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    protein =""
    if len(seq)%3 == 0: 
        protein = start_table[seq[0:3]]
        for i in range(3, len(seq), 3): 
            codon = seq[i:i + 3]
            try:
                protein+= table[codon]
            except KeyError:## codon was not planned for. Probably can't determine it.
                # print(codon)
                protein += 'X'## X is for unknown
    else:
        print(len(seq))
        raise IndexError("Given sequence length modulo 3 was different than 0, which is unexpected.")
    return protein 

def is_compressed(file_or_file_path):
    """
        Checks is a file, or file path given is compressed or not
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file,"rb")
    else:
        try:
            file = open(file.name,"rb")
        except:
            return False
    if file.read(2).startswith(b'\x1f\x8b'):
        return True
    file.close()
    return False

def read_compressed_or_not(file_or_file_path):
    """
        Copied from http://www.github.com/ggautreau/PPanGGOLiN.git's utils.py.
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file,"rb")
    else:
        try:
            file = open(file.name,"rb")
        except:
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        logging.getLogger().info("Uncompressing the file :'" + file.name + "' ...")
        return(TextIOWrapper(gzip.open(filename=file, mode = "r")))
    else:
        file.close()
        file = open(file.name,"r")
        return(file)

def write_compress_or_not(file_path, compress):
    if compress:
        return gzip.open(file_path +".gz", mode="wt")
    else:
        return open(file_path,"w")

def launch_aragorn(fnaFile):
    logging.getLogger().debug("Running Aragorn(tRNA).")
    cmd = ["aragorn","-t","-gcbact","-l","-w",fnaFile]
    p = Popen(cmd, stdout = PIPE)
    fileData = p.communicate()[0].decode().split("\n")[::-1]##loading the whole thing, reverting it to 'pop' in order.
    geneObjs = []
    basename = ''.join(os.path.basename(fnaFile).split(".")[:-1]) ## getting the file basename without the (expected) extention.
    c = 0
    while len(fileData) != 0:
        line = fileData.pop()
        if line.startswith(">"):
            header = line.replace(">","").split()[0]
            fileData.pop()## then next line must be removed too.
        elif len(line) > 0:## if the line isn't empty, there's data to get.
            lineData = line.split()
            start, stop = ast.literal_eval(lineData[2].replace("c",""))
            c +=1
            geneObjs.append(gene(ID = basename+'_tRNA_'+str(c).zfill(3),
                            contig = header,
                            start =  start,
                            stop = stop,
                            strand = "-" if lineData[2].startswith("c") else "+",
                            geneType = "tRNA",
                            inference = "COORDINATES:profile:Aragorn",
                            product=lineData[1] + lineData[4]))
    logging.getLogger().info(f"Done with Aragorn(tRNA). There are {len(geneObjs)} tRNAs.")
    return geneObjs

def launch_prodigal(fnaFile):
    logging.getLogger().debug("Running Prodigal(CDS).")
    cmd = ["prodigal","-f","sco","-c","-i",fnaFile,"-p","single", "-q"]
    p = Popen(cmd, stdout = PIPE)

    basename = ''.join(os.path.basename(fnaFile).split(".")[:-1]) ## getting the file basename without the (expected) extention.
    geneObjs = []
    c = 0
    for line in p.communicate()[0].decode().split("\n"):
        if line.startswith("# Sequence Data:"):
            for data in line.split(";"):
                if data.startswith("seqhdr"):
                    header = data.split("=")[1].replace('"',"").split()[0]
            
        if line.startswith(">"):
            c+=1
            lineData = line[1:].split("_")## not considering the '>'
            geneObjs.append(gene(ID = basename + "_CDS_" + str(c).zfill(4),
                            contig = header,
                            start = lineData[1],
                            stop = lineData[2],
                            strand = lineData[3],
                            geneType="CDS",
                            inference = "ab initio prediction:Prodigal"))
    logging.getLogger().info(f"Done with Prodigal. There are {len(geneObjs)} CDSs.")
    return geneObjs

def launch_infernal(fnaFile):
    logging.getLogger().debug("Running Infernal(rRNA).")
    tmpFile = tempfile.NamedTemporaryFile(mode="r")
    cmd = ["cmscan","--tblout",tmpFile.name,"--hmmonly","--cpu",str(1),"--noali",f"{os.path.dirname(os.path.realpath(__file__))}/cmDB/rRNA_bact.cm",fnaFile]
    p = Popen(cmd, stdout=open(os.devnull,"w"), stderr = PIPE)
    err = p.communicate()[1].decode()
    print(err.split())

    basename = ''.join(os.path.basename(fnaFile).split(".")[:-1]) ## getting the file basename without the (expected) extention.

    geneObjs = []
    c=0
    for line in tmpFile:
        if not line.startswith("#"):
            c+=1
            lineData = line.split()
            strand = lineData[9]
            if strand == "-":
                start = lineData[8]
                stop = lineData[7]
            else:
                start = lineData[7]
                stop = lineData[8]
            geneObjs.append(gene(ID = basename + "_rRNA_" + str(c).zfill(3),
                                contig = lineData[2],
                                start = start,
                                stop = stop,
                                strand = strand,
                                geneType = "rRNA",
                                product = " ".join(lineData[17:]) ))

    logging.getLogger().info(f"Done with Infernal(rRNA). There are {len(geneObjs)} rRNAs.")
    return geneObjs

def write_output(outputbasename, contigs, genes, compress, format, cpu):
    with Pool(processes = cpu) as p:
        for forms in format.split(","):
            if forms == "gff":
                wgff = p.apply_async(func=write_gff, args = (outputbasename, contigs, genes, compress))
            elif forms == "fna":
                wfna = p.apply_async(func=write_fna, args =(outputbasename, contigs, compress))
            elif forms == "ffn":
                wffn = p.apply_async(func=write_ffn, args = (outputbasename, contigs, genes, compress))
            elif forms == "faa":
                wfaa = p.apply_async(func =write_faa, args = (outputbasename, contigs, genes, compress))
        
        ## there is a better way of writing this 'wait until all processes are done', isn't there ?
        for forms in format.split(","):
            if forms == "gff":
                wgff.get()
            elif forms == "fna":
                wfna.get()
            elif forms == "ffn":
                wffn.get()
            elif forms == "faa":
                wfaa.get()

def write_gff(output, contigs, genes, compress):
    logging.getLogger().debug("Writting GFF file ...")
    outfile = write_compress_or_not(output + ".gff", compress)
    outfile.write("##gff-version 3\n")
    precContig = ""
    for gene in genes:
        if gene.contig != precContig:
            outfile.write(f"##sequence-region {gene.contig} 1 {len(contigs[gene.contig])}\n")
            outfile.write(f"{gene.contig}\tGenBank\tregion\t1\t{len(contigs[gene.contig])}\t.\t+\t.\tID={gene.contig}\n")### GenBank origin is only for PanGBank. The circularity could be guessed with a gbff input.
        outfile.write(gene.get_gff())
        precContig = gene.contig
    outfile.close()
    logging.getLogger().debug("Done writing GFF file.")

def write_ffn(output, contigs, genes, compress):
    logging.getLogger().debug("Writting FFN file ...")
    outfile = write_compress_or_not(output + ".fnn", compress)
    for gene in genes:
        if gene.type == "CDS":## ffn is only for coding sequences ! (even if some softs also include RNA sequences ... It was not designed as such)
            if gene.strand == "+":
                outfile.write(gene.get_ffn(contigs[gene.contig][gene.start-1:gene.stop]))
            elif gene.strand == "-":
                outfile.write(gene.get_ffn(reverse_complement(contigs[gene.contig][gene.start-1:gene.stop])))
    outfile.close()
    logging.getLogger().debug("Done writing FFN file.")

def write_faa(output, contigs, genes, compress):
    logging.getLogger().debug("Writting FAA file ...")
    outfile = write_compress_or_not(output + ".faa", compress)
    for gene in genes:
        if gene.type == "CDS":## faa is only for coding sequences ! 
            if gene.strand == "+":
                outfile.write(gene.get_faa(contigs[gene.contig][gene.start-1:gene.stop]))
            elif gene.strand == "-":
                outfile.write(gene.get_faa(reverse_complement(contigs[gene.contig][gene.start-1:gene.stop])))
    outfile.close()
    logging.getLogger().debug("Done writing FAA file.")

def write_fna(output, contigs, compress):
    logging.getLogger().debug("Writting FNA file ...")
    outfile = write_compress_or_not(output + ".fna", compress)
    for header in sorted(contigs.keys(), key = lambda x : len(contigs[x]), reverse = True):
        outfile.write(f">{header}\n")
        j = 0
        while j < len(contigs[header]):
            outfile.write(f"{contigs[header][j:j+60]}\n")
            j += 60
    logging.getLogger().debug("Done writing FNA file.")
    outfile.close()

def write_tmp_fasta(contigs):

    tmpFile = tempfile.NamedTemporaryFile(mode="w")
    for header in contigs.keys():
        tmpFile.write(f">{header}\n")
        j=0
        while j < len(contigs[header]):
            tmpFile.write(contigs[header][j:j+60]+"\n") 
            j+= 60
    return tmpFile

def read_fasta(fnaFile):

    logging.getLogger().debug("Reading fasta file and extracting sequence and contig names.")
    contigs = {}
    contig_name = ""
    contig_seq = ""
    for line in fnaFile:
        if line.startswith('>'):
            if contig_name != "" and contig_seq != "":
                contigs[contig_name] = contig_seq
            contig_seq = ""
            contig_name = line.split()[0][1:]
        else:
            contig_seq += line.strip()

    ## processing the last contig
    if contig_name != "" and contig_seq != "":
        contigs[contig_name] = contig_seq
    logging.getLogger().debug(f"Done reading fasta. There are {len(contigs)} sequences.")
    return contigs

def syntaxic_annotation(fna, cpu, norna):

    ##pre-processings
    fastaFile = read_compressed_or_not(fna)
    contigs = read_fasta(fastaFile)
    if is_compressed(fna):
        fastaFile = write_tmp_fasta(contigs)
    ##launching tools for syntaxic annotation
    genes = []
    logging.getLogger().info(f"Launching syntaxic annotation.")

    with Pool(processes = cpu ) as p:## launching a process pool with all the accessible processes        

        proGenes = p.apply_async(func = launch_prodigal, args=(fastaFile.name,))
        logging.getLogger().debug("Started the process launching Prodigal")
        if not norna:
            araGenes = p.apply_async(func = launch_aragorn, args=(fastaFile.name,))
            logging.getLogger().debug("Started the process launching Aragorn")
            infGenes = p.apply_async(func = launch_infernal, args = (fastaFile.name, ))
            logging.getLogger().debug("Started the process launching Infernal")
            genes.extend(araGenes.get())
            logging.getLogger().debug("Got the results from the process for Aragorn")
            genes.extend(infGenes.get())
            logging.getLogger().debug("Got the results from the process for Infernal")
        genes.extend(proGenes.get())
        logging.getLogger().debug("Got the results from the process for Prodigal")
    
    fastaFile.close()## closing either tmp file or original fasta file.
    return genes, contigs

def overlap_filter(genes, contigs, circular=False): ## for including circular information later, eventually.
    tmpGenes = sorted(genes, key = lambda x : (-len(contigs[x.contig]), x.start))
    rmGenes = set()

    for i, gene_i in enumerate(tmpGenes):
        if i+1 < len(tmpGenes):
            gene_j = tmpGenes[i+1]
            if gene_i.type != "CDS" and gene_j.type == "CDS" and gene_i.stop > gene_j.start and gene_i.contig == gene_j.contig:
                rmGenes.add(gene_j)
            elif gene_i.type == "CDS" and gene_j.type != "CDS" and gene_i.stop > gene_j.start and gene_i.contig == gene_j.contig:
                rmGenes.add(gene_i)

    for gene in rmGenes:
        tmpGenes.remove(gene)
    logging.getLogger().info(f"{len(rmGenes)} CDS were removed due to overlapping RNA features.")
    return tmpGenes

def mk_basename(output, fna):
    outbasename = output
    if output[-1] != "/":
        outbasename += "/"
    outbasename += "".join(os.path.basename(fna).split(".")[:-1])
    return outbasename

def check_versions(norna):
    ## Prodigal
    logging.getLogger().info("Checking software versions")
    try:
        proVersion = Popen(["prodigal","-v"], stderr = PIPE).communicate()[1].decode().split()[1][1:-1].split(".")
    except FileNotFoundError:
        raise Exception("Prodigal was not found in PATH. Please install it and/or add it to your PATH.")
    if not (int(proVersion[0]) == 2 and int(proVersion[1]) >= 6 and int(proVersion[2]) >= 2):
        raise Exception(f"Prodigal version is not good. It is read as {'.'.join(proVersion)} while version 2.6.2 at least is expected.")
    logging.getLogger().info(f"Prodigal's version is '{'.'.join(proVersion)}'")
    if not norna:
        try:
            araVersion = Popen(["aragorn","-h"], stdout = PIPE).communicate()[0].decode().split("\n")[1][9:15].split(".")
        except FileNotFoundError:
            raise Exception("Aragorn was not found in PATH. Please install it and/or add it to your PATH.")

        if not (int(araVersion[0]) == 1 and int(araVersion[1]) >= 2 and int(araVersion[2]) >= 38):
            raise Exception(f"Aragorn version is not good. It is read as {'.'.join(araVersion)} while version at least 1.2.38 is expected.")
        logging.getLogger().info(f"Aragorn's version is '{'.'.join(araVersion)}'")
        try:
            infVersion = Popen(["cmscan","-h"], stdout = PIPE).communicate()[0].decode().split('\n')[1][11:16].split(".")
        except FileNotFoundError:
            raise Exception("Infernal (and notably its executable cmscan) was not found in PATH. Please install it and/or add it to your PATH.")
        if not (int(infVersion[0]) == 1 and int(infVersion[1]) >= 1 and int(infVersion[2]) >= 2):
            raise Exception(f"Aragorn version is not good. It is read as {'.'.join(infVersion)} while version at least 1.1.2 is expected.")
        logging.getLogger().info(f"Infernal's version is '{'.'.join(infVersion)}'")

def cmdLine():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fna',  required = True, type=str, help="fasta file to annotate.")
    parser.add_argument('--overlap',required = False, action='store_false',default=True,help ="Use to not remove genes overlapping with RNA features.")
    parser.add_argument('--compress', required = False, action='store_true', default=False, help="Use to compress the output files.")
    parser.add_argument('--output', required = False, type=str, default="png_outputdir"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory path (optionnal)")
    parser.add_argument("--basename", required = False, type = str, help = "Basename to use for output files. If not provided, it will be guessed from the input file.")
    parser.add_argument("--norna",required = False, action="store_true",default = False, help = "Use to avoid annotating RNA features.")
    parser.add_argument("--cpu",required = False, type = int, default = 1, help="Number of cpus to use.")
    parser.add_argument("--format",required = False, type = str.lower, default = "gff", help = "Different formats that you want as output, separated by a ','. Accepted strings are : faa fna gff ffn.")
    parser.add_argument("--verbose",required = False, action="store_true",default = False, help = "Use to see the DEBUG log, which outputs uppon")
    args = parser.parse_args()

    for form in args.format.split(","):
        if form not in ["gff","faa","fna","ffn"]:
            parser.error(f"One of the provided format names was wrong. the format string was read as '{args.format}'. The substring that failed recognizing is '{form}'. Are only recognized : 'gff','fna','faa','ffn' string formats. Beware of spaces.")
    return args

def main():

    args = cmdLine()
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(stream=sys.stdout, level = level, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    check_versions(args.norna)

    genes, contigs = syntaxic_annotation(args.fna, args.cpu, args.norna)
    if args.overlap and not args.norna:
        genes = overlap_filter(genes, contigs)## sorting and removing CDS that overlap.
    else:
        genes = sorted(genes, key = lambda x : (-len(contigs[x.contig]), x.start))## need to sort the genes by contig size, then by start postion.

    logging.getLogger().info("Writting output files...")
    outbasename = mk_basename(args.output, args.fna) if not args.basename else os.path.abspath(args.output + "/" + args.basename)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    write_output(outbasename, contigs, genes, args.compress, args.format, args.cpu)
    
    logging.getLogger().info(f"There is a total of {len(genes)} annotated features.")
    ## from here, all processes are done.


if __name__ == "__main__":
    main()