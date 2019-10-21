import sys,re

def read_fasta(fl):
    ID=''
    seq=''
    pos_line=[]
    fasta={}
    pos={}
    for line in fl:
        if line.startswith(">"):
            if ID == '':
                ID=line.lstrip(">").rstrip().split()[0]
                pos[ID]=line.lstrip(">").rstrip().split()[1].split(';')
            else:
                fasta[ID]=seq
                ID=line.lstrip(">").rstrip().split()[0]
                seq=''
                pos[ID]=line.lstrip(">").rstrip().split()[1].split(';')
        else:
            seq+=line.rstrip()
    fasta[ID]=seq
    return fasta,pos

def get_alignments(fasta_file,ref,tar):
    fd=open(fasta_file)
    (fasta,pos)=read_fasta(fd)
    alignments=[]
    i=0
    j=0
    counti=0
    countj=1
    ref_s=int(pos[ref][i].split(':')[1].split(',')[0])
    tar_s=int(pos[tar][j].split(':')[1].split(',')[0])
    for n in range(0,len(fasta[ref])):
        if fasta[tar][n] == "-":
            if countj == 0:
                ref_s+=1
            else:
                align=pos[tar][j].split(':')[0]+' '+pos[tar][j].split(':')[1].split(',')[0]+' '+str(tar_s+countj)
                align+=' '+pos[ref][i].split(':')[0]+' '+str(ref_s)+' '+str(ref_s+countj)
                ref_s=ref_s+countj
                tar_s=tar_s+countj
                alignments.append(align)
                counti+=countj
                countj=0
        else:
            countj+=1
            if countj+counti == int(pos[tar][j].split(':')[1].split(',')[1]) - int(pos[tar][j].split(':')[1].split(',')[0]):
                align=pos[tar][j].split(':')[0]+' '+str(tar_s)+' '+pos[tar][j].split(':')[1].split(',')[1]
                align+=' '+pos[ref][i].split(':')[0]+' '+str(ref_s)+' '+str(ref_s+countj)
                alignments.append(align)
                ref_s=ref_s+countj
                j+=1
                try:
                    tar_s=int(pos[tar][j].split(':')[1].split(',')[0])
                except IndexError:
                    break
                countj=0
                counti=0
    return alignments

if __name__=='__main__':
    fasta_file = sys.argv[1]
    ref = sys.argv[2]
    tar = sys.argv[3]
    alignments=get_alignments(fasta_file,ref,tar)
    print '\n'.join(alignments)
