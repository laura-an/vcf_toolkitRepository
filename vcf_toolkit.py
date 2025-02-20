#!/usr/bin/env python3

# Example: ./vcf_toolkit.py search -input archivo.vcf -position x:xxxx 


import argparse 
import sys
import subprocess 
import pandas as pd
import os


## Definition of the class vcf_toolkit

class vcf_toolkit:

    # INIZIALIZATION OF AN OBJECT OF THE CLASS

    def __init__(self, path):
        self.path = path # the directory of the file is saved as a property of the object
        if os.path.exists(path): # It is checked whether the file exists, is compressed, and if it has an index in the same directory
            if  path.endswith('.gz'):
                self.compressed = True
                if os.path.exists(path+'.csi') or os.path.exists(path+'.gz.csi'):
                    self.indexed = True
                else:
                    self.indexed = False
            else: 
                self.compressed = False
                self.indexed = False
        else: 
            print(f'ERROR: The file {path} was not found or is not valid')
            sys.exit(1) # if the file was not found, it leaves the process
    

    # DEFINITION OF THE FUNCTIONS OF THE CLASS

    def nrows(self):
        '''counts de number of rows in the file'''
        nrows_h = subprocess.run(f'bcftools view -h {self.path} | wc -l', shell=True, capture_output=True, text=True).stdout.strip()
        nrows_nh = self.nrows_nh()
        # first, it calculates the number of rows for de header and then for the rest of the file using the function nrows_nh for optimazing optimizar. Then the sum is calculated
        nrows = int(nrows_h) + int(nrows_nh)
        return nrows
        

    def nrows_nh(self):
        '''nh stands for No Header, it counts the number of rows withount having the header into account'''
        # query is used instead of view to avoid processing the whole file and make it faster
        nrows_nh = subprocess.run(f'bcftools query -f "%CHROM\n" {self.path} | wc -l', shell=True, capture_output=True, text=True)
        return nrows_nh.stdout.strip()


    def search(self, position):
        '''returns the rows of the given position'''
        #if the file is indexed, the index is used to make it faster
        if self.indexed:
            ocurrences = subprocess.run(f'bcftools view -H -r {position} {self.path}'.split(), capture_output=True, text=True)
        else:
            ocurrences = subprocess.run(f'bcftools view -H -t {position} {self.path}'.split(), capture_output=True, text=True)
            print('Indexing the file is recommended to optimize the process')
        if ocurrences.stdout == '': #if there is not occurrences, the format of the position given is changed
            modify_position = lambda position: 'chr'+position if not position.startswith('chr') else position[3:]
            if self.indexed:
                ocurrences = subprocess.run(f'bcftools view -H -r {modify_position(position)} {self.path}'.split(), capture_output=True, text=True)
            else:
                ocurrences = subprocess.run(f'bcftools view -H -t {modify_position(position)} {self.path}'.split(), capture_output=True, text=True)
                print('Indexing the file is recommended to optimize the process')
        return ocurrences.stdout
    

    def nchrom(self):
        '''returns the number of rows per chromosome'''
        nchrom_raw = subprocess.run(f'bcftools query -f "%CHROM\n" {self.path}', shell=True, capture_output=True, text=True)
        nchrom_pd= pd.Series(nchrom_raw.stdout.strip().splitlines()) #to convert the list of chromosomes in a Series of pandas
        nchrom = nchrom_pd.value_counts().to_string()
        return nchrom
    

    def chrom(self, chr):
        '''returns all the rows of the given chromosome'''
        if self.indexed:
            ocurrences = subprocess.run(f'bcftools view -H -r {chr} {self.path}'.split(), capture_output=True, text=True)
        else:
            ocurrences = subprocess.run(f'bcftools view -H -t {chr} {self.path}'.split(), capture_output=True, text=True)
            print('Indexing the file is recommended to optimize the process')
        if ocurrences.stdout == '': #if there is not occurrences, the format of the chromosome given is changed
            modify_position = lambda chr: 'chr'+chr if not chr.startswith('chr') else chr[3:]
            if self.indexed:
                ocurrences = subprocess.run(f'bcftools view -H -r {modify_position(chr)} {self.path}'.split(), capture_output=True, text=True)
            else:
                ocurrences = subprocess.run(f'bcftools view -H -t {modify_position(chr)} {self.path}'.split(), capture_output=True, text=True)
                print('Indexing the file is recommended to optimize the process')
        return ocurrences.stdout
    

    def nsamples(self):
        '''returns the number of samples that are included in the file'''
        nsamples = subprocess.run(f'bcftools query -l {self.path} | wc -l', shell=True, capture_output=True, text=True)
        return nsamples.stdout.strip() 
    

    def samples(self):
        '''returns the name of the samples that are included in the file'''
        samples = subprocess.run(f'bcftools query -l {self.path}'.split(), capture_output=True, text=True)
        samples_list = samples.stdout.splitlines()
        if len(samples_list)==0:
            return 'There are no samples'
        elif len(samples_list)>=10:
            return f'{samples_list[:10]} and {len(samples_list)-10} more.'
        else:
            return samples_list
    

    def filter(self):
        '''returns the number of appearances of the different options of FILTER'''
        filter_raw = subprocess.run(f'bcftools query -f "%FILTER\n" {self.path}', shell=True, capture_output=True, text=True)
        filter_pd= pd.Series(filter_raw.stdout.strip().splitlines()) #to convert the list of chromosomes in a Series of pandas to use value_counts()
        filter = filter_pd.value_counts().to_string()
        return filter


    def description(self):
        '''returns a brief description of the file'''
        nrows = self.nrows()
        nrows_nh = self.nrows_nh()
        nsamples = self.nsamples()
        samples = self.samples()
        nchrom = self.nchrom()
        filter = self.filter()
        if samples=='There are no samples':
            return(f'-> Number of rows:\n{nrows}\n\n-> Number of rows excluding the header:\n{nrows_nh}\n\n-> Number of samples:\n{nsamples}\n\n'\
                f'-> Number of rows per chromosome:\n{nchrom}\n\n-> Number of PASS and FAIL in filter:\n{filter}')
        else:
            return(f'-> Number of rows:\n{nrows}\n\n-> Number of rows excluding the header:\n{nrows_nh}\n\n-> Number of samples:\n{nsamples}\n\n'\
                f'-> Samples:\n{samples}\n\n-> Number of rows per chromosome:\n{nchrom}\n\n-> Number of PASS and FAIL in filter:\n{filter}')
    

    def columns(self): 
        '''this function was created to take the columns to then open the file with pandas as a dataframe'''
        header = subprocess.run(f'bcftools view -h {self.path} | tail -n1', shell=True, capture_output=True, text=True) #para ver la última fila del header
        header_columns = header.stdout.strip().replace('#','').split('\t')
        return header_columns
    

    def readVcf(self):
        '''returns a dataframe of pandas of the vcf or bcf file'''
        vcf_panda = pd.read_csv(self.path, comment="#", sep="\t", header=None, names=self.columns())
        return vcf_panda


## BODY OF THE PROGRAM that will run if the script is run directly

def body():
    
    parser = argparse.ArgumentParser(description='Program to handle and analyse a vcf file with the following functions:\n'\
                                     'nrows: Returns the number of rows in the file\n'\
                                     'nrows_nh: Returns the number of rows in the file without having the header into account\n'\
                                     'search: Returns the rows of the given position (chr:pos[-endPos])\n'\
                                     'nchrom: Returns the number of rows per chromosome\n'\
                                     'chrom: Returns the rows of the given chromosome\n'\
                                     'nsamples: Returns the number of samples in the file\n'\
                                     'samples: Returns a list with the name of the samples in the file\n'\
                                     'filter: Returns the number of FAIL and PASS \n'\
                                     'description: Returns a brief description of the file'\
                                     'readVcf: Returns a pandas dataframe of the file')

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    parser_nrows = subparsers.add_parser("nrows", help="Returns the number of rows in the file")
    parser_nrows.add_argument("-input", help="Direction of the vcf file")

    parser_nrows_nh = subparsers.add_parser("nrows_nh", help="Returns the number of rows in the file without having the header into account")
    parser_nrows_nh.add_argument("-input", help="Direction of the vcf file")

    parser_search = subparsers.add_parser("search", help="Returns the rows of the given position (chr:pos[-endPos])")
    parser_search.add_argument("-input", help="Direction of the vcf file")
    parser_search.add_argument("-position", help="Position in the following format: chrx:xxxxxxxxxxx")

    parser_nchrom = subparsers.add_parser("nchrom", help="Returns the number of rows per chromosome")
    parser_nchrom.add_argument("-input", help="Direction of the vcf file")

    parser_chrom = subparsers.add_parser("chrom", help="Returns the rows of the given chromosome")
    parser_chrom.add_argument("-input", help="Direction of the vcf file")
    parser_chrom.add_argument("-chr", help="Chromosome that needs to be queried")

    parser_nsamples = subparsers.add_parser("nsamples", help="Returns the number of samples in the file")
    parser_nsamples.add_argument("-input", help="Direction of the vcf file")

    parser_samples = subparsers.add_parser("samples", help="Returns a list with the name of the samples in the file")
    parser_samples.add_argument("-input", help="Direction of the vcf file")

    parser_filter = subparsers.add_parser("filter", help="Returns the number of FAIL and PASS")
    parser_filter.add_argument("-input", help="Direction of the vcf file")

    parser_description = subparsers.add_parser("description", help="Returns a brief description of the file")
    parser_description.add_argument("-input", help="Direction of the vcf file")

    parser_readVcf = subparsers.add_parser("readVcf", help="Returns a pandas dataframe of the file")
    parser_readVcf.add_argument("-input", help="Direction of the vcf file")

    args = parser.parse_args()

    ##

    vcf = vcf_toolkit(args.input) #se crea un objeto de la clase vcf_toolkit con la direccion del archivo dado

    #se comprueba cual es el comando introducido para aplicar sobre el objecto la función apropiada y generar el output
    if args.command=="nrows":
        output = vcf.nrows()
    elif args.command=="nrows_nh":
        output = vcf.nrows_nh()
    elif args.command=="search":
        output = vcf.search(args.position)
    elif args.command=="nchrom":
        output = vcf.nchrom()
    elif args.command=="chrom":
        output = vcf.chrom(args.chr)
    elif args.command=="nsamples":
        output = vcf.nsamples()
    elif args.command=="samples":
        output = vcf.samples()
    elif args.command=="filter":
        output = vcf.filter()
    elif args.command=="description":
        output = vcf.description()
    elif args.command=="readVcf":
        output = vcf.readVcf()


    
    print(output)
    



if __name__ == "__main__":
    body()