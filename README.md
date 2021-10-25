## Этап 1
Был произведен выбор случайных ридов с помощью команд:

```bash
seqtk sample -s1304 oil_R1.fastq 5000000 > re_oil_R1.fastq
seqtk sample -s1304 oil_R2.fastq 5000000 > re_oil_R2.fastq
seqtk sample -s1304 oilMP_S4_L001_R1_001.fastq 1500000 > re_oilMP_S4_L001_R1_001.fastq
seqtk sample -s1304 oilMP_S4_L001_R2_001.fastq 1500000 > re_oilMP_S4_L001_R2_001.fastq
```

Запуск *fastqc*:
```bash
mkdir fastq
ls re*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}
```

Запуск *multiqc*:
```bash
mkdir multiqc
multiqc -o multiqc fastqc
```

Запуск *platanus*:
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

## Этап 2
### Сравнение результатов до/после 

#### До:
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/MultiQC_general_stats.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Mean_quality_scores.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Adapter_content.png?raw=true)

#### После:
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_general_stats.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_mean_quality_scores.png?raw=true)
![alt text](https://github.com/IlonaGA/hse21_hw1/blob/main/Images/Trimmed_adapter_content.png?raw=true)

#### Выводы:
1) Заметно уменьшилась длина последовательностей (см. General statistics, Length).
2) Улучшилось качество: теперь только зеленая зона (см. Mean quality scores).
3) Были отрезаны адаптеры (см. Adapter content). 

## Этап 3
### Контиги
Сбор *контигов*:
```bash
time platanus assemble -o Poil -t 2 -m 28 -f re_oil_R1.fastq.trimmed re_oil_R2.fastq.trimmed 2> assembl.log
```

Анализ полученных контигов:
```python
import numpy as np
```

```python
def analysis(contig):
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
    
    N50 = length_array[np.cumsum(length_array) >= value][0]
    
    print('Общее количество контигов: ', counter)
    print('Суммарная длина контигов: ', length)
    print('Длина самого длинного контига: ', max(length_array))
    print('N50: ', N50)
```

```python
contig_file = open('Poil_contig.fa', 'r')
contig = contig_file.readlines()
analysis(contig)
```

Результат:
```
Общее количество контигов:  620
Суммарная длина контигов:  3926318
Длина самого длинного контига:  179307
N50:  55863
```
### Скаффолды
Сбор *скаффолдов*:
```
time platanus scaffold -o Poil -t 2 -c Poil_contig.fa -IP1 re_oil_R1.fastq.trimmed re_oil_R2.fastq.trimmed -OP2 re_oilMP_S4_L001_R1_001.fastq.int_trimmed re_oilMP_S4_L001_R2_001.fastq.int_trimmed 2> scaffold.log
```
Анализ полученных скаффолдов:
```python
import numpy as np
```

```python
def analysis_s(scaffold):
    counter = 0
    length = 0
    length_array = []
    for line in scaffold:
        if line[0] == '>':
            counter += 1
            length += int(line.split('_')[1][3:])
            length_array.append(int(line.split('_')[1][3:]))

    length_array = np.asarray(length_array)
    length_array = np.sort(length_array)[::-1]
    value = np.sum(length_array) / 2
    
    N50 = length_array[np.cumsum(length_array) >= value][0]
    
    print('Общее количество скаффолдов: ', counter)
    print('Суммарная длина скаффолдов: ', length)
    print('Длина самого длинного скаффолда: ', max(length_array))
    print('N50: ', N50)
```

```python
scaffold_file = open('Poil_scaffold.fa', 'r')
scaffold = scaffold_file.readlines()
analysis_s(scaffold)
```

Результат:
```
Общее количество скаффолдов:  72
Суммарная длина скаффолдов:  3875885
Длина самого длинного скаффолда:  3831756
N50:  3831756
```

Самый длинный скаффолд:
```python
for line in scaffold:
    if int(line.split('_')[1][3:]) == 3831756:
        print(line)
        break
```
Получим: 
```
>scaffold1_len3831756_cov232
```

```bash
echo scaffold1_len3831756_cov232 > max_scaffold.txt
seqtk subseq Poil_scaffold.fa max_scaffold.txt > max_scaffold.fa
```
### Гэпы
Количество и общая длина гэпов:
```python
max_scaffold_file = open('max_scaffold.fa', 'r')
max_scaffold = max_scaffold_file.readlines()
count = 0
for i in range (len(max_scaffold[1]) - 1):
    if (max_scaffold[1][i+1].upper() == 'N') and (max_scaffold[1][i].upper() != 'N'):
        count += 1
print('Количество гэпов: ', count)
print('Общая длина гэпов: ', max_scaffold[1].count('N'))
```
Результат:
```
Количество гэпов: 62
Общая длина гэпов: 6083
```

Уменьшаем количество генов и повторяем анализ:
```bash
time platanus gap_close -o Poil -t 2 -c Poil_scaffold.fa -IP1 re_oil_R1.fastq.trimmed re_oil_R2.fastq.trimmed -OP2 re_oilMP_S4_L001_R1_001.fastq.int_trimmed re_oilMP_S4_L001_R2_001.fastq.int_trimmed 2> gapclose.log
```
```bash
echo scaffold1_cov232 > max_scaffold.txt
seqtk subseq Poil_gapClosed.fa max_scaffold.txt > max_scaffold_lowgap.fa
```
```python
scaffold_file_lowgap = open('max_scaffold_lowgap.fa', 'r')
scaffold_lowgap = scaffold_file_lowgap.readlines()

count = 0
for i in range (len(scaffold_lowgap[1]) - 1):
    if (scaffold_lowgap[1][i+1].upper() == 'N') and (scaffold_lowgap[1][i].upper() != 'N'):
        count += 1
print('Количество гэпов: ', count)
print('Общая длина гэпов: ', scaffold_lowgap[1].count('N'))
```
Результат:
```
Количество гэпов:  8
Общая длина гэпов:  1500
```
Видно, что значительно количество гэпов, действительно, было удалено.


