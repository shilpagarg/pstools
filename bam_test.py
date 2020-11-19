import pysam
samfile = pysam.AlignmentFile("../bwa_hic.bam", "rb", threads = 48)


result = [[0 for _ in range(100000)] for _ in range(100000)]
pre_name = ''
current = []
for b in samfile.fetch(until_eof=True):
    if(pre_name == b.query_name):
        current.append(b.reference_id)
    else:
        pre_name = b.query_name
        for i in current:
            for j in current:
                result[i][j]+=1
        current = []
    # print(b.query_name)
    # result[b.reference_id]+=1
    # print(b.is_paired)
samfile.close()

for i in range(len(result)):
    for j in range(len(result[0])):
        if(result[i][j]>0):
            print(str(i)+"\t"+str(j)+"\t"+str(result[i][j]))