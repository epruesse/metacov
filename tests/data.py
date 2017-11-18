import os

path = os.path.join(os.path.dirname(__file__), 'data')
fq1 = os.path.join(path, 'ecoli_1K_1.fq.gz')
fq2 = os.path.join(path, 'ecoli_1K_2.fq.gz')
bam = os.path.join(path, 'bbmap.sorted.bam')
regions_blast = os.path.join(path, 'regions.blast7')
ref = os.path.join(path, 'reference_1K.fa.gz')
