#!/usr/bin/env python3
## creo que la línea de arriba permite que el archivo se ejecute como
## un script de python sin especificar el intérprete



# Se ejecuta como ./vcf_toolkit.py search -input archivo.vcf

## Importacion de argparse para trabajar en linea de comandos

import argparse
import sys
import subprocess
import pandas as pd

## Definicion clase vcf_toolkit.py

class vcf_toolkit:

    def __init__(self, archivo):
        self.archivo = archivo # se guarda la direccion del archivo como un apropiedad del objeto
        try:
            with open(archivo, 'r') as f:
                lineas = f.readlines()
        except FileNotFoundError:
            print(f'ERROR: El archivo {archivo} no se pudo encontrar')
            sys.exit(1) #para salir del programa con un código de error

        self.lineas = lineas # se guarda cada línea del archivo en una propiedad del objeto
    

    def nrows(self):
        #cuenta el número de filas del archivo
        nfilas = len(self.lineas)
        return nfilas
        
    
    def nrows_nh(self):
        #nh es de no header, cuenta el número de filas sin tener en cuenta el header
        nfilas_nh = 0
        nfilas_nh = sum (1 for l in self.lineas if l[0]!='#')
        #for l in self.lineas:
            #if l[0]!='#':
                #nfilas_nh += 1
        return nfilas_nh


    def search_position(self, position):
        try:
            ocurrencias = subprocess.run(f'bcftools view -H -r {position} {self.archivo}.gz'.split(), capture_output=True, text=True)
            return ocurrencias
        except FileNotFoundError:
            return 'Necesitas comprimir el archivo .vcf antes utilizando bzgip para que el programa pueda encontrarlo en el mismo directorio en el que está el archivo original'
    
    def chrom(self, chr):
        ocurrencias=''
        for l in self.lineas:
            linea_sep = l.split()
            if linea_sep[0]==chr:
                ocurrencias += (l+ '\n')
        return ocurrencias
    
    def columnas(self):
        header = subprocess.run(f'bcftools view -h {self.archivo}'.split(), capture_output=True, text=True) #para coger el header
        header_columnas = header.stdout.strip().split('\n')[-1] #para ver la última fila del header
        header_columnas_lista = header_columnas.replace('#','').split('\t') #para crear una lista con todos los elementos de la última fila del header
        return header_columnas_lista

    


## cuerpo del programa que se ejecutará si el script está siendo ejecutado directamente

def body():
    
    parser = argparse.ArgumentParser(description="Programa para manejar y analizar un archivo vcf con las siguientes funciones:\n\
                                     nrows: Devuelve el número total de filas en el archivo\n\
                                     nrows_nh: Devuelve el número de filas en el archivo sin tener en cuenta el header\n\
                                     search: Devuelve las filas en las que aparece la posición introducida (chrx:xxxxxxxxx)\n\
                                     chrom: Devuelve las filas del cromosoma introducido\n\
                                     nchrom: Devuelve el número de filas por cromosoma\n\
                                     nsamples: Devuelve el número de muestras\n\
                                     filter: Devuelve el número de FAIL PASS o .")

    subparsers = parser.add_subparsers(dest="command", help="Subcomandos disponibles")

    parser_nrows = subparsers.add_parser("nrows", help="Devuelve el numero de filas")
    parser_nrows.add_argument("-input", help="Direccion del archivo vcf")

    parser_nrows_nh = subparsers.add_parser("nrows_nh", help="Devuelve el numero de filas sin contar el header")
    parser_nrows_nh.add_argument("-input", help="Direccion del archivo vcf")

    parser_search = subparsers.add_parser("search", help="Busca la posicion introducida en formato chrx:xxxxxxxxxx y devuelve las filas en las que aparezca")
    parser_search.add_argument("-input", help="Direccion del archivo vcf")
    parser_search.add_argument("-position", help="Posicion en formato chrx:xxxxxxxxxxx")

    parser_chrom = subparsers.add_parser("chrom", help="Busca un cromosoma y devuelve todas las filas con información de ese cromosoma")
    parser_chrom.add_argument("-input", help="Direccion del archivo vcf")
    parser_chrom.add_argument("-chr", help="Cromosoma que se quiere buscar")

    parser_nsamples = subparsers.add_parser("nsamples", help="Devuelve el número de muestras en el archivo")
    parser_nsamples.add_argument("-input", help="Direccion del archivo vcf")

    parser_nchrom = subparsers.add_parser("nchrom", help="Cuenta el número de filas de cada cromosoma")
    parser_nchrom.add_argument("-input", help="Direccion del archivo vcf")

    parser_filter = subparsers.add_parser("filter", help="Cuenta los FAIL y PASS de filter")
    parser_filter.add_argument("-input", help="Dirección del archivo vcf")

    args = parser.parse_args()

    #

    vcf = vcf_toolkit(args.input) #se crea un objeto de la clase vcf_toolkit con la direccion del archivo dado
    vcf_panda = pd.read_csv(vcf.archivo, comment="#", sep="\t", header=None, names=vcf.columnas())

    if args.command=="nrows":
        output = vcf.nrows()
    elif args.command=="nrows_nh":
        output = vcf.nrows_nh()
    elif args.command=="search":
        output = vcf.search_position(args.position)
    elif args.command=="chrom":
        output = vcf.chrom(args.chr)
    elif args.command=="nsamples":
        ncolumnas = len(vcf.columnas())
        output = ncolumnas -9 #Siempre hay 9 columnas de datos y el resto de muestras o puede variar????
    elif args.command=="nchrom":
        output = vcf_panda['CHROM'].value_counts()
    elif args.command=="filter":
        output = vcf_panda['FILTER'].value_counts()


    
    print(output)
    



if __name__ == "__main__":
    body()