Этап 1
Выбор случайных ридов проводился командами:
seqtk sample -s1304 oil_R1.fastq 5000000 > sub_PE_R1.fastq
seqtk sample -s1304 oil_R2.fastq 5000000 > sub_PE_R2.fastq
seqtk sample -s227 oilMP_S4_L001_R1_001.fastq 1500000 > sub_MP_R1.fastq
seqtk sample -s227 oilMP_S4_L001_R2_001.fastq 1500000 > sub_MP_R2.fastq
