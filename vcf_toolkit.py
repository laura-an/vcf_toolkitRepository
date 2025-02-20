#!/usr/bin/env python3

# Se ejecuta, por ejemplo, ./vcf_toolkit.py search -input archivo.vcf -position x:xxxx


import argparse 
import sys
import subprocess 
import pandas as pd
import os


## DEFINICIÓN DE LA CLASE vcf_toolkit

class vcf_toolkit:

    # INICIALIZACIÓN DE UN OBJETO DE LA CLASE

    def __init__(self, path):
        self.path = path # se guarda la direccion del archivo como una propiedad del objeto
        if os.path.exists(path): # se comprueba si el archivo existe, está comprimido y si tiene un índice en el mismo directorio
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
            print(f'ERROR: El archivo {path} no se pudo encontrar o no es válido')
            sys.exit(1) # si no se encuentra el archivo se sale del proceso
    

    # DEFINICIÓN DE FUNCIONES DE LA CLASE

    def nrows(self):
        '''cuenta el número de filas del archivo'''
        # no hay que hacer nada distinto con vcf o vcf.gz porque bcftools deja utilizar esta línea en ambos casos
        nrows_h = subprocess.run(f'bcftools view -h {self.path} | wc -l', shell=True, capture_output=True, text=True).stdout.strip()
        nrows_nh = self.nrows_nh()
        # primero se calcula el número de filas en el header y luego en el resto de archivo usando la función nrows_nh para optimizar y se suma
        nrows = int(nrows_h) + int(nrows_nh)
        return nrows
        

    def nrows_nh(self):
        '''nh es de no header, cuenta el número de filas sin tener en cuenta el header'''
        # no hay que hacer nada distinto con vcf o vcf.gz porque bcftools deja utilizar esta línea en ambos casos
        # se usa query en vez de view para no tener que procesar el archivo entero y que sea más rápido
        nrows_nh = subprocess.run(f'bcftools query -f "%CHROM\n" {self.path} | wc -l', shell=True, capture_output=True, text=True)
        return nrows_nh.stdout.strip()


    def search(self, position):
        '''devuelve las líneas de la posición introducida'''
        #si el archivo está indexado se trabaja con el índice porque es más rápido, si no, no
        if self.indexed:
            ocurrences = subprocess.run(f'bcftools view -H -r {position} {self.path}'.split(), capture_output=True, text=True)
        else:
            ocurrences = subprocess.run(f'bcftools view -H -t {position} {self.path}'.split(), capture_output=True, text=True)
            print('Indexing the file is recommended to optimize the process')
        if ocurrences.stdout == '': #si no hay ocurrencias se prueba cambiando el formato de entrada de la posición
            modify_position = lambda position: 'chr'+position if not position.startswith('chr') else position[3:]
            if self.indexed:
                ocurrences = subprocess.run(f'bcftools view -H -r {modify_position(position)} {self.path}'.split(), capture_output=True, text=True)
            else:
                ocurrences = subprocess.run(f'bcftools view -H -t {modify_position(position)} {self.path}'.split(), capture_output=True, text=True)
                print('Indexing the file is recommended to optimize the process')
        return ocurrences.stdout
    

    def nchrom(self):
        '''Devuelve el número de filas por cromosoma'''
        nchrom_raw = subprocess.run(f'bcftools query -f "%CHROM\n" {self.path}', shell=True, capture_output=True, text=True)
        nchrom_pd= pd.Series(nchrom_raw.stdout.strip().splitlines()) #para convertir la lista de cromosomas en una serie de pandas
        nchrom = nchrom_pd.value_counts().to_string()
        return nchrom
    

    def chrom(self, chr):
        '''Devuelve todas las filas del cromosoma introducido'''
        if self.indexed:
            ocurrences = subprocess.run(f'bcftools view -H -r {chr} {self.path}'.split(), capture_output=True, text=True)
        else:
            ocurrences = subprocess.run(f'bcftools view -H -t {chr} {self.path}'.split(), capture_output=True, text=True)
            print('Indexing the file is recommended to optimize the process')
        if ocurrences.stdout == '': #si no hay ocurrencias se prueba cambiando el formato de entrada de la posición
            modify_position = lambda chr: 'chr'+chr if not chr.startswith('chr') else chr[3:]
            if self.indexed:
                ocurrences = subprocess.run(f'bcftools view -H -r {modify_position(chr)} {self.path}'.split(), capture_output=True, text=True)
            else:
                ocurrences = subprocess.run(f'bcftools view -H -t {modify_position(chr)} {self.path}'.split(), capture_output=True, text=True)
                print('Indexing the file is recommended to optimize the process')
        return ocurrences.stdout
    

    def nsamples(self):
        '''Devuelve el número de muestras que hay en el archivo'''
        nsamples = subprocess.run(f'bcftools query -l {self.path} | wc -l', shell=True, capture_output=True, text=True)
        return nsamples.stdout.strip() 
    

    def samples(self):
        '''Devuelve el nombre de las muestras en el archivo'''
        samples = subprocess.run(f'bcftools query -l {self.path}'.split(), capture_output=True, text=True)
        samples_list = samples.stdout.splitlines()
        if len(samples_list)==0:
            return 'There are no samples'
        elif len(samples_list)>=10:
            return f'{samples_list[:10]} and {len(samples_list)-10} more.'
        else:
            return samples_list
    

    def filter(self):
        '''Devuelve el número de apariciones de las diferentes opciones de FILTER'''
        filter_raw = subprocess.run(f'bcftools query -f "%FILTER\n" {self.path}', shell=True, capture_output=True, text=True)
        filter_pd= pd.Series(filter_raw.stdout.strip().splitlines()) #para convertir la lista de cromosomas en una serie de pandas para poder usar value_counts()
        filter = filter_pd.value_counts().to_string()
        return filter


    def description(self):
        '''Devuelve una breve descripción del archivo'''
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
        '''esta función se creó para tomar las columnas para luego abrir el archivo por pandas como dataframe'''
        header = subprocess.run(f'bcftools view -h {self.path} | tail -n1', shell=True, capture_output=True, text=True) #para ver la última fila del header
        header_columns = header.stdout.strip().replace('#','').split('\t')
        return header_columns
    

    def readVcf(self):
        '''Devuelve un data frame de pandas del archivo vcf o bcf'''
        vcf_panda = pd.read_csv(self.path, comment="#", sep="\t", header=None, names=self.columns())
        return vcf_panda


## CUERPO DEL PROGRAMA que se ejecutará si el script está siendo ejecutado directamente

def body():
    
    parser = argparse.ArgumentParser(description='Programa para manejar y analizar un archivo vcf con las siguientes funciones:\n'\
                                     'nrows: Devuelve el número total de filas en el archivo\n'\
                                     'nrows_nh: Devuelve el número de filas en el archivo sin tener en cuenta el header\n'\
                                     'search: Devuelve las filas en las que aparece la posición introducida (chr:pos[-endPos])\n'\
                                     'nchrom: Devuelve el número de filas por cromosoma\n'\
                                     'chrom: Devuelve las filas del cromosoma introducido\n'\
                                     'nsamples: Devuelve el número de muestras\n'\
                                     'samples: Devuelve una lista con el nombre de las muestras presentes en el archivo\n'\
                                     'filter: Devuelve el número de FAIL PASS o \n'\
                                     'description: Devuelve una breve descripción del archivo'
                                     'readVcf: Devuelve un data frame de pandas del archivo vcf o bcf, comprimido o no en gz')

    subparsers = parser.add_subparsers(dest="command", help="Subcomandos disponibles")

    parser_nrows = subparsers.add_parser("nrows", help="Devuelve el numero de filas")
    parser_nrows.add_argument("-input", help="Direccion del archivo vcf")

    parser_nrows_nh = subparsers.add_parser("nrows_nh", help="Devuelve el numero de filas sin contar el header")
    parser_nrows_nh.add_argument("-input", help="Direccion del archivo vcf")

    parser_search = subparsers.add_parser("search", help="Busca la posicion introducida en formato chrx:xxxxxxxxxx y devuelve las filas en las que aparezca")
    parser_search.add_argument("-input", help="Direccion del archivo vcf")
    parser_search.add_argument("-position", help="Posicion en formato chrx:xxxxxxxxxxx")

    parser_nchrom = subparsers.add_parser("nchrom", help="Cuenta el número de filas de cada cromosoma")
    parser_nchrom.add_argument("-input", help="Direccion del archivo vcf")

    parser_chrom = subparsers.add_parser("chrom", help="Busca un cromosoma y devuelve todas las filas con información de ese cromosoma")
    parser_chrom.add_argument("-input", help="Direccion del archivo vcf")
    parser_chrom.add_argument("-chr", help="Cromosoma que se quiere buscar")

    parser_nsamples = subparsers.add_parser("nsamples", help="Devuelve el número de muestras en el archivo")
    parser_nsamples.add_argument("-input", help="Direccion del archivo vcf")

    parser_samples = subparsers.add_parser("samples", help="Devuelve una lista con las muestras en el archivo")
    parser_samples.add_argument("-input", help="Direccion del archivo vcf")

    parser_filter = subparsers.add_parser("filter", help="Cuenta los FAIL y PASS de filter")
    parser_filter.add_argument("-input", help="Dirección del archivo vcf")

    parser_description = subparsers.add_parser("description", help="Devuelve una descripción del archivo vcf")
    parser_description.add_argument("-input", help="Dirección del archivo vcf")

    parser_readVcf = subparsers.add_parser("readVcf", help="Devuelve un data frame de pandas del archivo vcf bcf o gz")
    parser_readVcf.add_argument("-input", help="Dirección del archivo vcf")

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