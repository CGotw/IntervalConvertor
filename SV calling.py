import gfapy, os, re
from tqdm import tqdm


def get_len(nodes) -> int:
    L = 0
    for n in nodes:
        L+=g1.line(n).LN
    return L

def find_node(chr_name,start,end):
    start,end = int(start), int(end)
    targets = []#空列表
    for node in node_dict[chr_name]:#用chr_name取了多个值，进行循环迭代
        s,e = node[0],node[1]
        if e<=start:
            continue
        if s>=end:
            break
        else:
            targets.append(node[-1])#取最后一个值
    return targets

def is_single(node):
    for n in node:
        in_degree = degree[int(n[1:])][0]
        out_degree = degree[int(n[1:])][1]
        if in_degree !=1 or out_degree!=1:
            return False
    return True

def write_fa(file, id, line, query_node):
    chr_name = line[0]
    start,end = map(int,line[1:3])
    bed_node = line[11].split(',')
    file.write('>')
    file.write("{}_{} {} {} {} {}".format(id,chr_name,start,end,bed_node[0],bed_node[-1]))
    #输出设定的格式，一个括号代表一个值，类似{id}_{chr_name}
    file.write('\n')
    strs = ''
    for n in query_node:
        strs+=g1.line(n).sequence#合并字符串
    while len(strs)>60:
        file.write(strs[:60])#只在文件中写入前60个字符
        file.write('\n')
        strs = strs[60:]#舍弃这60个#继续写前60个，循环，直至不满60个。每行只写60个
        
    file.write(strs[:60])#循环结束后，最后剩下的继续写
    file.write('\n')
    

def SV_calling(args):
    gfa = gfapy.Gfa()
    # gaf_dir = '/data/lilab/DATA/Apple_pan/analysis/11_Graph_alignmt_new/'#命名问题
    # gfa_path = '/data/lilab/DATA/Apple_pan/analysis/08_Pan_genome/02.phased_genome-new/Huahong-purged-new.phased-pan.36.gfa'
    gaf_dir=args.gaf_dir
    gfa_path=args.gfa_path
    g1 = gfa.from_file(gfa_path)

    f = open(args.query_bed_path, 'r')
    res = f.readlines()
    f.close()

    for line in res:#按行遍历，先提取每一行，再分割
        tmp = line.rstrip('\n').split('\t')
        if int(tmp[1]) >= 1632104 and int(tmp[2]) <= 1634133:#选取第一位置为1632104和第二位置为1634133之间的序列进行打印
            print(line)

    node_dict = dict()
    for name in tqdm(g1.names):#相当于range迭代器，显示进度条
        s = g1.line(name)#将name提取
        SN = s.SN   #s结构里的某个数，产生副本于SN，下同
    #     SN = SN.split('_')[0]
        LN = s.LN
        SO = s.SO
        if SN in node_dict:
            node_dict[SN].append((SO,SO+LN,name))#如果在字典内，将键SN对应的值加上一个有三个数元组
        else:
            node_dict[SN] = [(SO,SO+LN,name)]#如果不存在，就创建这个键，对应一个值
            #print(node_dict[SN])

    degree = [[0,0] for _ in range(len(g1.names)+1)]#数组的大小是g1.name的大小＋1,每一个数据都是前面的[0,0]，初始化
    for line in tqdm(g1.edges):
        froms = line.from_name
        tos = line.to_name
        degree[int(tos[1:])][0]+=1#tos[1:]取tos从1到末尾所有数，转成int类型，此数组值自增1
        degree[int(froms[1:])][1]+=1


    min_len = 50
    max_len = 50


    output_pre_path=args.output_dir
    deletion_path=os.path.join(output_pre_path,"Deletion.fa")
    insertion_path=os.path.join(output_pre_path,"Insertion.fa")
    divergent_path=os.path.join(output_pre_path,"Divergent.fa")
    multiallelic_path=os.path.join(output_pre_path,"Multiallelic.fa")
    f1 = open(deletion_path, 'w')
    f2 = open(insertion_path, 'w')
    f3 = open(divergent_path, 'w')
    f4 = open(multiallelic_path, 'w')

    c1,c2,c3,c4 = 0,0,0,0
    for line in tqdm(res):
        line = line.rstrip('\n').split('\t')
        chr_name = line[0]
        start,end = map(int,line[1:3])#start相当于取了line[1],end相当于取了line[2]
        bed_node = line[11].split(',')#line中第12个字符串进行分割
        ref_node = find_node(chr_name,start,end)
        query_node = set(bed_node[1:-1])-set(ref_node)#两个集合做差
        ref_len = get_len(ref_node)
        query_len = get_len(query_node)
        if ref_len >= min_len and query_len<=max_len:
            # Deletion
    #         continue
            if is_single(ref_node):
                c1+=1
                write_fa(f1,c1,line, ref_node)  #原片段
        elif ref_len <= max_len and query_len>=min_len:
            # Insertion
    #         continue
            if is_single(query_node):
                c2+=1
                write_fa(f2,c2,line, query_node) #插入片段   
        elif ref_len >= min_len and query_len>= min_len:
            # Divergent
            if is_single(ref_node) and is_single(query_node):#如果原片段为1和插入片段都为1时
                c3+=1
    #             continue
                write_fa(f3,str(c3)+'_1',line, ref_node)  #原片段
                write_fa(f3,str(c3)+'_2',line, query_node)    #插入片段  
            # Multiallelic
            elif is_single(ref_node) and not  is_single(query_node):#如果原片段为1且插入片段不为1时
                c4+=1
    #             continue
                write_fa(f4,str(c4)+'_1',line, ref_node)#原片段
                write_fa(f4,str(c4)+'_2',line, query_node)#插入片段
                print(f4)
                
        
    f1.close()
    f2.close()
    f3.close()
    f4.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='input parameters')
    parser.add_argument('--gaf_dir', type=str,help ='gaf dictionary')
    parser.add_argument('--gfa_path', type=str,help='query input')
    parser.add_argument('--query_bed_path', type=str,default=None,help='query bed path')
    parser.add_argument('--output_dir', type=str,help='output dictionary')    
    args = parser.parse_args()
    SV_calling(args=args)