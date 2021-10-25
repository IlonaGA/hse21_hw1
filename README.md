## Этап 1 ##
Выбор случайных ридов проводился командами:

```
seqtk sample -s1304 oil_R1.fastq 5000000 > re_oil_R1.fastq
seqtk sample -s1304 oil_R2.fastq 5000000 > re_oil_R2.fastq
seqtk sample -s1304 oilMP_S4_L001_R1_001.fastq 1500000 > re_oilMP_S4_L001_R1_001.fastq
seqtk sample -s1304 oilMP_S4_L001_R2_001.fastq 1500000 > re_oilMP_S4_L001_R2_001.fastq
```

Запуск fastqc:

```
mkdir fastq
ls re*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}
```

Запуск multiqc:
```
mkdir multiqc
multiqc -o multiqc fastqc
```

Запуск platanus:
```
platanus_trim re_oil_R1.fastq re_oil_R2.fastq
platanus_internal_trim re_oilMP_S4_L001_R1_001.fastq re_oilMP_S4_L001_R2_001.fastq
```

```
mkdir trimmed_fastqc
ls *trimmed | xargs -P 4 -tI{} fastqc -o trimmed_fastqc {}
mkdir trimmed_multiqc
multiqc -o trimmed_multiqc trimmed_fastqc
```

## Этап 2 ##
### Сравнение результатов до/после ###

До:
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/MultiQC_general_stats.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Mean_quality_scores.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Adapter_content.png?raw=true)

После:
