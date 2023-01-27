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


def interval_conversion(args):

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


    # 调试用
    # f = open('/data/lilab/sywang1/dierpian/Fuji_0/Fuji_0.gaf', 'r') # 因为是从gd转换坐标，所以要用gd的gaf信息，而不是all_gaf
    # f = open('/data/lilab/sywang1/dierpian/unPhased_gaf/M9_0.gaf', 'r')
    f=open(args.query_gaf_path,'r')
    gd_gaf = f.readlines()
    f.close()
    # fin = open('/home/rdzhou1/shiyao/input/GD-9-16.txt', 'r') # 输入
    fin = open(args.query_input_path, 'r')

    gene = fin.readlines()
    fin.close()


    # 调试用
    # f = open('/home/rdzhou1/shiyao/output/M9-9-16-out.txt', 'w') # 输出
    f = open(args.output_file_path, 'w') # 输出
    for g in tqdm(gene):
        g = g.rstrip('\n').split('\t')
        chr_name = g[0].lower()
        gene_name = g[3]
        for line in gd_gaf:
            line = line.rstrip('\n').split('\t')
            start,end =  get_intersection(int(g[1]),int(g[2]),int(line[2]),int(line[3]))
            # if start!=-1:
            #     f.write("line[0]="+str(line[0].lower())+" chr_name="+str(chr_name)+" start="+str(start)+"\n")
            if line[0].lower() == chr_name and start!=-1:
                # print("line[0]=",line[0].lower(),"chr_name=",chr_name,"start=",start)
                real_start = start+int(line[7])-int(line[2])
                real_end = end+int(line[7])-int(line[2])
                mat1 = "{}\t{}\t{}"
                f.write('原始对象 '+mat1.format(chr_name,g[1],g[2]))
                f.write('\n')
                f.write('查询结果 '+mat1.format(line[0], start,end))
                f.write('\n')
                # map_res是空的了，为什么还要判断line[5]
                map_res = re.split('[><]',line[5])[1:]
                mat2 = "{:<20}\t{:20}{:20}"
                if map_res == []:
                    if line[5][:4]!= 'Fuji':
                        pass
                    else:
                        f.write(mat2.format(line[5],real_start,real_end))
                        f.write('\n')
                tmp_len = int(line[7])
                pre_len = int(line[7])
                for p in map_res:
                    tmp = p.split(':')
                    n = tmp[0]
                    # if n[:4]!= 'Fuji': #
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
                f.write('---------------------------')
                f.write('\n')
                break
    f.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='input parameters')
    parser.add_argument('--gaf_dir', type=str,help ='gaf directory')
    parser.add_argument('--gfa_path', type=str,help='query input')
    parser.add_argument('--query_input_path', type=str,help='query input path')
    parser.add_argument('--query_gaf_path', type=str,help='query gaf path')
    parser.add_argument('--output_file_path', type=str,help='output file path')    
    args = parser.parse_args()
    interval_conversion(args=args)

    