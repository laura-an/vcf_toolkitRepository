#!/usr/bin/env python3
## creo que la línea de arriba permite que el archivo se ejecute como
## un script de python sin especificar el intérprete


## Importacion de argparse para trabajar en linea de comandos

import argparse
import sys


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
        for l in self.lineas:
            if l[0]!='#':
                nfilas_nh += 1
        return nfilas_nh


    def search_position(self, position):
        position_sep = position.split(':')
        ocurrencias=''
        for l in self.lineas:
            linea_sep = l.split()
            if position_sep[0]==linea_sep[0] and position_sep[1]==linea_sep[1]:
                ocurrencias += (l +'\n')
        return ocurrencias
    
    def chrom(self, chr):
        ocurrencias=''
        for l in self.lineas:
            linea_sep = l.split()
            if linea_sep[0]==chr:
                ocurrencias += (l+ '\n')
        return ocurrencias



## cuerpo del programa que se ejecutará si el script está siendo ejecutado directamente

def body():
    
    parser = argparse.ArgumentParser(description="Programa para analizar archivo vcf")

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

    args = parser.parse_args()

    #

    vcf = vcf_toolkit(args.input) #se crea un objeto de la clase vcf_toolkit con la direccion del archivo dado

    if args.command=="nrows":
        output = vcf.nrows()
    elif args.command=="nrows_nh":
        output = vcf.nrows_nh()
    elif args.command=="search":
        output = vcf.search_position(args.position)
    elif args.command=="chrom":
        output = vcf.chrom(args.chr)


    
    print(output)
    



if __name__ == "__main__":
    body()