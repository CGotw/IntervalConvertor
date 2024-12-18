import gfapy, os, re, pandas, openpyxl
from tqdm import tqdm
from openpyxl.styles import PatternFill
import argparse

def find_node(chr_name,start,end):
    start,end = int(start), int(end)
    targets = []
    for node in node_dict[chr_name]:
        s,e = node[0],node[1]
        if e<=start:
            continue
        if s>=end:
            break
        else:
            targets.append(node[-1])
    return targets

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
    gfa = gfapy.Gfa() 
    gaf_dir=args.gaf_dir
    gfa_path=args.gfa_path
    g1 = gfa.from_file(gfa_path)
    reference_name=args.reference_name
    reference_name_len=len(reference_name)

    allgaf = []
    all_file_name = set()
    for file in tqdm(os.listdir(gaf_dir)):
        if os.path.isdir(gaf_dir+file):
            continue
        all_file_name.add(file.split('.')[0])
        f = open(gaf_dir+file, 'r')
        allgaf.extend([[file]]+f.readlines())
        f.close()

    node_dict = dict()
    for name in tqdm(g1.names):
        s = g1.line(name)
        SN = s.SN
        LN = s.LN
        SO = s.SO
        if SN in node_dict:
            node_dict[SN].append((SO,SO+LN,name))
        else:
            node_dict[SN] = [(SO,SO+LN,name)]

    all_node = []
    for line in tqdm(allgaf):
        #print(line)

        if len(line) == 1:
            all_node.append([line[0][:-4]])
            continue
        line = line.rstrip('\n').split('\t')
        map_res = re.split('[><]',line[5])[1:]
        if map_res == []:
            continue
        all_node[-1].append([line[0]])
        res_dict = dict()
        for p in map_res:
            tmp = p.split(':')
            name = tmp[0]
            pre_name = name.split('_')[0]
            if pre_name != reference_name:
                continue
            [start,end] = tmp[1].split('-')
            if name in res_dict:
                res_dict[name].append((start,end))
            else:
                res_dict[name]=[(start,end)]
        all_node[-1][-1].append(res_dict)

    #print(all_node)
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
    parser.add_argument('--gaf_dir', type=str,help ='gaf directory', required=True)
    parser.add_argument('--gfa_path', type=str,help='query input', required=True)
    parser.add_argument('--query_input_path', type=str,help='query input path', required=True)
    parser.add_argument('--query_gaf_path', type=str,help='query gaf path', required=True)
    parser.add_argument('--reference_name', type=str,help='reference name', required=True)    
    args = parser.parse_args()
    interval_conversion(args=args)