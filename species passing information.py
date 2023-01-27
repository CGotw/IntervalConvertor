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

def get_intersection(a,b,c,d): # 求[a,b] [c,d]交集，返回坐标
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


def find_all(chr_name,start,end):
    targets = []
    for node in node_dict[chr_name]:
        s,e = node[0],node[1]
        if e<=start:
            continue
        if s>=end:
            break
        else:
            targets.append(node)
    # print(targets)
    all_name = set()
    all_res = []
    for t in targets:
        all_res.append([chr_name,t])
        # print(all_res)
    #     print(seg)
        for file in all_node:
            file_name = file[0]
            # print(file)
            tmp_res = []
            for line in file[1:]: 
                if chr_name not in  line[1].keys():
                    continue
                for data in line[1][chr_name]:
                    # if int(data[0]) <= t[0] and t[1]<=int(data[1]):
                    if (int(data[0])-t[0]) * ( t[1]-int(data[1])) >=0:
    #                     print(data,t)
                        tmp_res.append([line[0]]+[data])
            if tmp_res != []:
                all_res[-1].append([file_name]+tmp_res)
    for node in all_res:
        f.write(node[0]+' '+node[1][-1])
        f.write('\n')
        for file in node[2:]:
            f.write(file[0]+'-->')
            for chr in file[1:]:
                start,end = chr[1]
                f.write(chr[0]+'  ')
            f.write('\n')
        f.write('-----------------------------------\n')

def run(chr_name,start,end,gene_name,f):
    start,end = int(start), int(end)
    targets = []
    
    for node in node_dict[chr_name]:
        s,e = node[0],node[1]
        if e<=start:
            continue
        if s>=end:
            break
        else:
            targets.append(node)
            gene_dict[gene_name].add(node[-1])
    all_name = set()
    all_res = []
    for t in targets:
        # if t[-1] in exon_node:
        all_res.append([chr_name,t])
        for file in all_node:
            file_name = file[0]
            tmp_res = []
            for line in file[1:]:
                if chr_name not in line[1].keys():
                    continue
                for data in line[1][chr_name]:
                    # if int(data[0]) <= t[0] and t[1]<=int(data[1]):
                    if (int(data[0])-t[0]) * (t[1]-int(data[1])) >=0:
    #                     print(data,t)
                        tmp_res.append([line[0]]+[data])
            if tmp_res != []:
                all_res[-1].append([file_name]+tmp_res)
    for node in all_res:
        f.write(node[0]+' '+node[1][-1])
        f.write('\n')
        tmp_set = set()
        for file in node[2:]: 
            tmp_set.add(file[0])
        for pn in sorted(list(tmp_set)):# 如果不用求差集，只要经过的物种，就用这个
        # for pn in sorted(list(all_file_name-tmp_set)):# 这一步是求差集，也就是从所有物种中去掉经过的物种
            f.write(pn+'\n')
        f.write('-----------------------------------\n')
        if node[1][-1] not in excel_dict:
            excel_dict[node[1][-1]] = set()
        excel_dict[node[1][-1]].update(all_file_name-tmp_set)

def jingguowuzhong1(args):

    # common
    gfa = gfapy.Gfa() 
    # gaf_dir = '/data/lilab/DATA/Apple_pan/analysis/11_Graph_alignmt_new/'
    # gfa_path = '/data/lilab/DATA/Apple_pan/analysis/08_Pan_genome/01.unphased_genome-new/Huahong-purged-new.unphased-pan.36.gfa'
    gaf_dir=args.gaf_dir
    gfa_path=args.gfa_path
    g1 = gfa.from_file(gfa_path)

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
    #     SN = SN.split('_')[0]
        LN = s.LN
        SO = s.SO
        if SN in node_dict:
            node_dict[SN].append((SO,SO+LN,name))
        else:
            node_dict[SN] = [(SO,SO+LN,name)]


    all_node = []
    for line in tqdm(allgaf):
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
            if pre_name != 'Fuji':
                continue
            [start,end] = tmp[1].split('-')
            if name in res_dict:
                res_dict[name].append((start,end))
            else:
                res_dict[name]=[(start,end)]
        all_node[-1][-1].append(res_dict)


    # path='/home/rdzhou1/shiyao/input/Shanzha-hap1-10-24.txt'
    path=args.query_input_path
    fin = open(path, 'r')
    gene = fin.readlines()
    fin.close()

    # f=open('/home/rdzhou1/shiyao/output/Shanzha-hap1-10-24-out-res-trans-derectly.txt','w')
    f=open(args.output_file_path,'w')
    for line in tqdm(gene):
        f.write(line)
        f.write('\n')
        line = line.rstrip('\n').split('\t')
        chr_name = line[0]
        start, end = int(line[1]),int(line[2])
        find_all(chr_name,start,end)
        f.write('================================================\n')
    f.close()


def jingguowuzhong2(args):
    # common
    gfa = gfapy.Gfa() 
    # gaf_dir = '/data/lilab/DATA/Apple_pan/analysis/11_Graph_alignmt_new/'
    # gfa_path = '/data/lilab/DATA/Apple_pan/analysis/08_Pan_genome/01.unphased_genome-new/Huahong-purged-new.unphased-pan.36.gfa'
    gaf_dir=args.gaf_dir
    gfa_path=args.gfa_path
    g1 = gfa.from_file(gfa_path)

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
    #     SN = SN.split('_')[0]
        LN = s.LN
        SO = s.SO
        if SN in node_dict:
            node_dict[SN].append((SO,SO+LN,name))
        else:
            node_dict[SN] = [(SO,SO+LN,name)]


    all_node = []
    for line in tqdm(allgaf):
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
            if pre_name != 'Fuji':
                continue
            [start,end] = tmp[1].split('-')
            if name in res_dict:
                res_dict[name].append((start,end))
            else:
                res_dict[name]=[(start,end)]
        all_node[-1][-1].append(res_dict)
    
    fin = open(args.query_input_path, 'r')
    gene = fin.readlines()
    fin.close()
    # fin = open('/home/rdzhou1/apple/Fuji_exon_location.txt', 'r')
    # exon = fin.readlines()
    # fin.close()
    # f = open('/data/lilab/sywang1/dierpian/Phased_gaf/Fuji_2.gaf', 'r') # 因为是从gd转换坐标，所以要用gd的gaf信息，而不是all_gaf
    f=open(args.query_gaf_path, 'r')
    gd_gaf = f.readlines()
    f.close()

    # f = open('/home/rdzhou1/shiyao/pianduan/more.txt', 'w')
    f=open(args.output_file_path,'w')
    excel_dict = dict()
    gene_dict = dict()
    for g in tqdm(gene):
        g = g.rstrip('\n').split('\t')
        chr_name = g[0].lower()
        gene_name = g[3]
        gene_dict[gene_name] = set()
        for line in gd_gaf:
            line = line.rstrip('\n').split('\t')
            start,end =  get_intersection(int(g[1]),int(g[2]),int(line[2]),int(line[3]))
            if line[0].lower() == chr_name and start!=-1:
                real_start = start+int(line[7])-int(line[2])
                real_end = end+int(line[7])-int(line[2])
                mat1 = "{}\t{}\t{}"
                f.write('原始对象 '+mat1.format(chr_name,g[1],g[2]))
                f.write('\n')
                f.write('查询结果 '+mat1.format(line[0], start,end))
                f.write('\n')
                map_res = re.split('[><]',line[5])[1:]
                mat2 = "{:<20}\t{:20}{:20}"
                if map_res == []:
                    if line[5][:4]!= 'Fuji':
                        pass
                    else:
                        f.write(mat2.format(line[5],real_start,real_end))
                        f.write('\n')
                        run(line[5],real_start-1,real_end+1,gene_name,f)
                tmp_len = int(line[7])
                pre_len = int(line[7])
                for p in map_res:
                    tmp = p.split(':')
                    n = tmp[0]
                    # if n[:4]!= 'Fuji':
                    #     continue
                    [s,e] = map(int,tmp[1].split('-'))
                    pre_len = tmp_len
                    tmp_len += (e-s)
                    tmp_s,tmp_e = get_intersection(pre_len,tmp_len,real_start,real_end)
                    if tmp_s == -1:
                        continue
                    if n[:4]== 'Fuji':
                        f.write(mat2.format(n,tmp_s-pre_len+s,tmp_e-pre_len+s))
                        f.write('\n')
                        run(n,tmp_s-pre_len+s-1,tmp_e-pre_len+s+1,gene_name,f)
                f.write('---------------------------')
                f.write('\n')
                break
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='input parameters')
    parser.add_argument('--gaf_dir', type=str,help ='gaf directory')
    parser.add_argument('--gfa_path', type=str,help='query input')
    parser.add_argument('--interval_conversion',type=bool,default=False,help="if interval_conversion is ture,查询物种前要先经过坐标转换")
    parser.add_argument('--query_input_path', type=str,help='query input path')
    parser.add_argument('--query_gaf_path', type=str,default=None,help='query gaf path')
    parser.add_argument('--output_file_path', type=str,help='output file path')    
    args = parser.parse_args()
    if args.interval_conversion:
        jingguowuzhong2(args=args)
    else:
        jingguowuzhong1(args=args)

    