import gfapy, os, re, pandas, openpyxl
from tqdm import tqdm
from openpyxl.styles import PatternFill
import argparse


def get_intersection(a,b,c,d):
    if (a-c)*(d-b)>=0:
        if a>=c and b<=d:
            return a,b
        else:
            return c,d
    elif a<c<b<d:
        return c,b
    elif c<a<d<b:
        return a,d
    else:
        return -1,-1

def interval_conversion(args):
    reference_name=args.reference_name
    reference_name_len=len(reference_name)


    f=open(args.query_gaf_path,'r')
    gd_gaf = f.readlines()
    f.close()
    fin = open(args.query_input_path, 'r')

    gene = fin.readlines()
    fin.close()

    outmapfile="../results/map.out.txt"
    print(outmapfile)
    f = open(outmapfile, 'w')
    for g in tqdm(gene):
        g = g.rstrip('\n').split('\t')
        chr_name = g[0].lower()
        gene_name = g[3]
        for line in gd_gaf:
            line = line.rstrip('\n').split('\t')
            start,end =  get_intersection(int(g[1]),int(g[2]),int(line[2]),int(line[3]))
            if line[0].lower() == chr_name and start!=-1:
                real_start = start+int(line[7])-int(line[2])
                real_end = end+int(line[7])-int(line[2])
                mat1 = "{}\t{}\t{}"
                f.write('Raw '+mat1.format(chr_name,g[1],g[2]))
                f.write('\n')
                f.write('Target '+mat1.format(line[0], start,end))
                f.write('\n')
                map_res = re.split('[><]',line[5])[1:]
                mat2 = "{:<20}\t{:20}{:20}"
                if map_res == []:
                    if line[5][:reference_name_len]!= reference_name:
                        pass
                    else:
                        f.write(mat2.format(line[5],real_start,real_end))
                        f.write('\n')
                tmp_len = int(line[7])
                pre_len = int(line[7])
                for p in map_res:
                    tmp = p.split(':')
                    n = tmp[0]
                    [s,e] = map(int,tmp[1].split('-'))
                    pre_len = tmp_len
                    tmp_len += (e-s)
                    tmp_s,tmp_e = get_intersection(pre_len,tmp_len,real_start,real_end)
                    if tmp_s == -1:
                        continue
                    if n[:reference_name_len]== reference_name:
                        f.write(mat2.format(n,tmp_s-pre_len+s,tmp_e-pre_len+s))
                        f.write('\n')
                f.write('---------------------------')
                f.write('\n')
                break
    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='input parameters')
    parser.add_argument('--query_input_path', type=str,help='query input path', required=True)
    parser.add_argument('--query_gaf_path', type=str,help='query gaf path', required=True)
    parser.add_argument('--reference_name', type=str,help='reference name', required=True)    
    args = parser.parse_args()
    interval_conversion(args=args)