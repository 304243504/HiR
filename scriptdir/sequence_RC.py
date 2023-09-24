
import sys

input  = sys.argv[1]
output = sys.argv[2]

def to_reverse(str):
    str1=list(str)
    str1.reverse()
    str2=''
    for i in str1:
        if i in ['A','a']:
            str2+='T'
        elif i in ['T','t']:
            str2+='A'
        elif i in ['c',"C"]:
            str2+='G'
        elif i in ['G','g']:
            str2+="C"
        elif i in ['U','u']:
            str2+='A'
        elif i in ['N','n']:
            str2+='N'
        else:
            print('please check your gene_sequence')
    return str2


def only_reverse(stringg):
    str5=list(stringg)
    str5.reverse()
    str6=''
    for i in str5:
        str6+=i
    return str6


def get_geneSeq(fq):
    list_reverse=[]
    for i in range(len(fq)):
        if i%4==1:
            try:
                r_seq=to_reverse(fq[i])
                list_reverse.append(r_seq)
            except Exception:
                print('somewhere error')
    return list_reverse


def get_valueSeq(fq):
    list_value=[]
    for i in range(len(fq)):
        if i%4==3:
            list_value.append(only_reverse(fq[i]))
    return list_value


def write_file(list_reverse,list_value,fp):
    for i in range(len(list_reverse)):
        fp[i*4+1]=list_reverse[i]
    for i in range(len(list_value)):
        fp[i*4+3]=list_value[i]
    with open(output, 'w') as f:
        for i in fp:
            f.write(i+'\n')


with open(input,'r') as sp:
    fq = str(sp.read()).replace('\n', ' ').split(' ')
    write_file(get_geneSeq(fq),get_valueSeq(fq),fq)

