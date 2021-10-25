## Этап 1 ##
Выбор случайных ридов проводился командами:

```bash
seqtk sample -s1304 oil_R1.fastq 5000000 > re_oil_R1.fastq
seqtk sample -s1304 oil_R2.fastq 5000000 > re_oil_R2.fastq
seqtk sample -s1304 oilMP_S4_L001_R1_001.fastq 1500000 > re_oilMP_S4_L001_R1_001.fastq
seqtk sample -s1304 oilMP_S4_L001_R2_001.fastq 1500000 > re_oilMP_S4_L001_R2_001.fastq
```

Запуск fastqc:

```bash
mkdir fastq
ls re*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}
```

Запуск multiqc:
```bash
mkdir multiqc
multiqc -o multiqc fastqc
```

Запуск platanus:
```bash
platanus_trim re_oil_R1.fastq re_oil_R2.fastq
platanus_internal_trim re_oilMP_S4_L001_R1_001.fastq re_oilMP_S4_L001_R2_001.fastq
```

```bash
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
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_general_stats.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_mean_quality_scores.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_adapter_content.png?raw=true)

## Этап 3 ###
Сбор контигов:
```bash
time platanus assemble -o Poil -t 2 -m 28 -f re_oil_R1.fastq.trimmed re_oil_R2.fastq.trimmed 2> assembl.log
```

Анализ полученных контигов:
```python
import numpy as np
```

```python
def amount(contig):
    counter = 0
    length = 0
    length_array = []
    for line in contig:
        if line[0] == '>':
            counter += 1
            length += int(line.split('_')[1][3:])
            length_array.append(int(line.split('_')[1][3:]))

    length_array = np.asarray(length_array)
    length_array = np.sort(length_array)[::-1]
    value = np.sum(length_array) / 2
    N50 = length_array[np.cumsum(length_array) <= value][-1]
    
    print('Общее количесвтво контигов: ', counter)
    print('Суммарная длина контигов: ', length)
    print('Длина самого длинного контига: ', max(length_array))
    print('N50: ', N50)
```
