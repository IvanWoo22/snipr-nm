# SNIPR-Nm — Small RNA Nm Ratio Analysis Pipeline

SNIPR-Nm（Small RNA Integration-based Polyadenylation and Ligation for Nm detection）是配套小 RNA 3' 端 2'-O-甲基化（Nm）检测方法的高通量分析流程。同步输出每个小 RNA 的**表达量**和 **Nm 修饰比例**，适用于植物、动物组织或任何 small RNA 文库。

> 相关文献：Tang et al. "SNIPR-Nm reveals small RNA methylation ratio" (bioRxiv, 2025)

---

## 目录结构

```
snipr-nm/
├── SNIPR-Nm.nf          # 流程主文件
├── nextflow.config      # 参数配置文件
├── samplesheet.csv      # 样本列表（需填写）
├── ath_mirna_neo.fa     # miRNA 参考序列（替换为你自己的）
├── bin/
│   ├── fastq2count.pl   # FASTQ → 序列计数表
│   ├── mirna_count.pl   # miRNA 比对计数 + poly-A 检测
│   ├── snipr_qc.py      # 小 RNA-seq 质控分析
│   └── aggregate_qc.py  # 汇总所有样本的质控报告
└── README.md
```

---

## 背景：什么是 SNIPR-Nm？

3' 端 2'-O-甲基化（Nm）是调控小 RNA 稳定性与功能的重要修饰。植物中 HEN1 介导的甲基化几乎覆盖所有小 RNA；动物中 Nm 主要出现在 piRNA 上，近期研究也发现哺乳动物 miRNA 存在部分甲基化修饰。

**现有技术的局限：**
- 传统 small RNA-seq：依赖接头连接，Nm 修饰的 RNA 因连接效率低而被**低估**
- Stem-loop RT-qPCR：可精确定量特定小 RNA，但**无法报告甲基化状态**
- β-elimination 或 poly(A) 加尾 qPCR：可检测 3' Nm，但只能针对特定靶标、需要大上样量、无法绝对定量

**SNIPR-Nm 的解决方案：**
利用 poly(A) polymerase 和一种特殊连接酶的组合，在同一反应中同步检测小 RNA 的**表达量**和 **Nm 修饰比例**。具体原理：
1. poly(A) polymerase 对 RNA 进行加尾处理
2. Nm 修饰的 RNA 因 3' 端被保护而**不被加尾**；无修饰的 RNA 被添加 poly-A 尾
3. 通过连接酶差异扩增，区分修饰与未修饰的分子
4. 结合高通量测序，实现全转录组层面的 Nm ratio 定量

该方法已成功应用于小鼠睾丸组织，获得可靠的 piRNA 和 miRNA Nm 修饰比例数据。

---

## 第一步：安装依赖软件

### 用 conda（推荐）

```bash
# 创建环境
conda create -n snipr -c bioconda -c conda-forge \
    nextflow cutadapt pear tsv-utils "python>=3.9"

# 激活并安装 Python 依赖
conda activate snipr
pip install biopython

# pear 建议从 KantorLab 频道安装（版本更新）：
# conda install -n snipr -c kantorlab pear
# 若无该频道，用 bioconda 版本也可以：
# conda install -n snipr -c bioconda pear

# 确认所有工具就位
which nextflow cutadapt pear python3 perl tsv-join
```

**各软件作用：**
| 软件 | 用途 |
|------|------|
| `nextflow` | 流程运行引擎 |
| `cutadapt` | 去除接头序列 |
| `pear` | 合并 PE reads |
| `perl` | 运行 miRNA 计数脚本 |
| `python3` + biopython | 质控脚本 |
| `tsv-utils` | 合并多个样本的计数表 |

### 直接使用系统已有的软件

如果系统已装好上述工具且有 biopython，可以跳过 conda。用 `which <命令>` 确认能找得到即可。

---

## 第二步：准备样本列表

编辑 `samplesheet.csv`，填入你自己的数据路径：

```csv
sample_id,R1,R2
U1,/home/用户名/data/U1_1.fq.gz,/home/用户名/data/U1_2.fq.gz
U2,/home/用户名/data/U2_1.fq.gz,/home/用户名/data/U2_2.fq.gz
```

- `sample_id`：样本名（避免中文和特殊字符）
- `R1` / `R2`：FASTQ 文件的**完整绝对路径**

**如何获取绝对路径？**
```bash
cd /你的数据目录
pwd    # 打印绝对路径
ls *.fq.gz  # 确认文件名
```

---

## 第三步：确认软件路径

```bash
conda activate snipr
which cutadapt   # 输出如 /home/xxx/miniconda3/envs/snipr/bin/cutadapt
which pear
which python3
which perl
```

把这些路径填入 `nextflow.config` 对应行（如 `params.cutadapt_bin = "完整路径"`）。如果命令在 PATH 中能直接找到，用命令名即可（如 `"cutadapt"`），无需完整路径。

---

## 流程架构

在运行流程前，建议先了解整体处理步骤：

```
输入 FASTQ（PE 小 RNA 文库）
    ↓
[Cutadapt]  去除接头序列
    ↓
[PEAR]       合并 PE reads
    ↓
[fastq2count] 转为序列计数表
    ↓
[mirna_count] 与 miRNA 参考比对
    ↓           同时检测 poly-A 尾（= 无 Nm 修饰的分子）
    ↓
[annotate]    为每个样本生成唯一文件名（避免名称冲突）
    ↓
[merge]       所有样本横向合并成一个表
    ↓
[QC]          质控分析（长度分布、质量、复杂度）
    ↓
[multiqc]     汇总所有样本的质控报告
```

**各步骤说明：**
| 步骤 | 功能 | 输出文件 |
|------|------|---------|
| Cutadapt | 去除接头序列，去除低质量 reads | clean_R1/R2.fq.gz |
| PEAR | 合并配对-end reads | assembled.assembled.fastq |
| fastq2count | 序列去重计数 | reads.count.tsv |
| mirna_count | 与 miRNA 参考比对，检测 poly-A 尾 | mirna.count.tsv |
| annotate | 添加样本名前缀 | mirna.{样本}.annotated.tsv |
| merge | 合并所有样本结果 | mirna.merged.tsv |
| QC | 质控分析 | *_qc_metrics.json, *_qc_report.html |
| multiqc | 生成汇总报告 | multiqc_report.html |

---

## 第四步：运行流程

```bash
cd /path/to/snipr-nm

nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results
```

**断点续跑**（中断后恢复，跳过已完成的步骤）：
```bash
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    -resume
```

---

## 参数说明

### 必需参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--samplesheet` | `samplesheet.csv` | 样本列表文件路径（CSV格式，需包含sample_id、R1、R2三列） |
| `--mirna_fa` | `ath_mirna_neo.fa` | miRNA 参考序列 FASTA 文件路径 |
| `--outdir` | `./results` | 输出目录路径（**建议使用绝对路径**，避免相对路径导致的文件找不到问题） |

### 可选参数

#### 接头序列参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--adaptor_r1` | `AGATCGGAAGAGCACA` | R1 端接头序列（小 RNA 文库常用 Illumina Small RNA TruSeq 接头） |
| `--adaptor_r2` | `GATCGTCGGACTGTAGAA` | R2 端接头序列 |

#### 线程和性能参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--cutadapt_threads` | `8` | Cutadapt 步骤使用的 CPU 线程数 |
| `--pear_threads` | `4` | PEAR 步骤使用的 CPU 线程数 |

#### PEAR 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--pear_min_ovlp` | `17` | PEAR 最小 overlap 长度（bp）。当 R1 和 R2 的重叠区域 ≥ 此长度时才会被合并。建议范围：15-25，对于小 RNA-seq 通常使用 17-20 |

### 高级参数

#### cutadapt 高级选项

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--cutadapt_opts` | `-O 5 -m 17 -e 0.18 --discard-untrimmed` | Cutadapt 的高级选项字符串，各选项含义： |

**`--cutadapt_opts` 详细说明：**

| 选项 | 值 | 说明 |
|------|-----|------|
| `-O 5` | 5 | 最小接头 overlap 长度。只有当接头与序列的重叠 ≥ 5bp 时才认为是有效接头 |
| `-m 17` | 17 | 最小 read 长度（bp）。短于此长度的 reads 会被丢弃。miRNA 通常为 18-25nt，建议设置 15-17 |
| `-e 0.18` | 0.18 | 最大错误率（允许的不匹配数/overlap长度）。0.18 表示约 18% 的错误容忍度，平衡灵敏度和特异性 |
| `--discard-untrimmed` | - | 丢弃没有检测到接头的 reads。这些可能是低质量或非目标片段 |

**其他可用 cutadapt 选项（需手动修改 `nextflow.config`）：**

| 选项 | 示例值 | 说明 |
|------|--------|------|
| `-q` | `20,20` | 质量修剪阈值，格式为 "前端,后端"。移除质量低于此值的碱基 |
| `--max-n` | `0.5` | 最大 N 比例（0-1）。允许的 N 碱基比例 |
| `-M` | `30` | 最大 read 长度。长于此长度的 reads 会被丢弃 |

### 常用覆盖示例

```bash
# 基础运行（使用绝对路径输出）
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/用户名/srna_results

# 更换接头序列（如果你的文库使用了不同的接头）
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --adaptor_r1 你的R1接头序列 \
    --adaptor_r2 你的R2接头序列

# 调整线程数（根据服务器配置）
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --cutadapt_threads 16 \
    --pear_threads 8

# 调整 PEAR overlap 参数（对于特殊文库）
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --pear_min_ovlp 20

# 完整参数示例
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/用户名/srna_results \
    --adaptor_r1 AGATCGGAAGAGCACA \
    --adaptor_r2 GATCGTCGGACTGTAGAA \
    --cutadapt_threads 16 \
    --pear_threads 8 \
    --pear_min_ovlp 17

# 断点续跑（中断后恢复）
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    -resume
```

### 参数配置文件

对于频繁使用的参数组合，建议修改 `nextflow.config` 文件：

```groovy
// nextflow.config 示例配置
params {
    samplesheet    = "samplesheet.csv"
    mirna_fa       = "ath_mirna_neo.fa"
    outdir         = "/home/用户名/srna_results"  // 使用绝对路径
    adaptor_r1     = "AGATCGGAAGAGCACA"
    adaptor_r2     = "GATCGTCGGACTGTAGAA"
    cutadapt_threads = 16
    pear_threads     = 8
    pear_min_ovlp    = 17
}
```

修改配置文件后，运行流程时只需指定必要参数：

```bash
nextflow run SNIPR-Nm.nf --samplesheet samplesheet.csv --mirna_fa ath_mirna_neo.fa
```

---

## 输出结果

### 目录结构

```
results/
├── U1/
│   ├── cutadapt.log              # 接头去除日志
│   ├── reads.count.tsv           # 序列计数表
│   ├── mirna.count.tsv           # miRNA 比对结果
│   ├── mirna.U1.annotated.tsv   # 带样本名的注释文件
│   ├── U1_qc_metrics.json        # 质控指标（JSON）
│   ├── U1_qc_report.html         # 详细质控 HTML 报告
│   ├── U1_quality_by_pos.json    # 每位点质量数据
│   └── U1_length_dist.json        # 长度分布数据
├── U2/ ... (每个样本一个子目录)
├── mirna.merged.tsv              # 所有样本横向合并总表
├── multiqc_report.html           # 汇总质控报告
└── multiqc_data/                 # 质控原始数据
```

### 各输出文件详细说明

#### 1. cutadapt.log（接头去除日志）

记录 Cutadapt 处理 FASTQ 文件的统计信息：

```log
This is cutadapt 4.9 with Python 3.11.0
Command line parameters: -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAGAA ...
Command line parameters: -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAGAA -O 5 -m 17 -e 0.18 ...
Finished 4 processes (0.4s wall time, 14.6 MB peak memory)

=== Adapter 1 (first read) ===
Sequence: AGATCGGAAGAGCACA; Type: regular 3'; Trimmed: 1250 reads

=== Adapter 2 (second read) ===
Sequence: GATCGTCGGACTGTAGAA; Type: regular 3'; Trimmed: 980 reads

=== Summary ===
Reads written (passing filters): 85000
```

**关键指标解读：**

| 指标 | 说明 | 参考值 |
|------|------|--------|
| Adapter 1/2 Trimmed | 成功去除接头的 reads 数量 | - |
| Reads written | 通过质量过滤的最终 reads 数量 | - |
| 未去除接头 reads 比例 | 使用 `--discard-untrimmed` 丢弃的 reads | < 20% 为正常 |

#### 2. reads.count.tsv（序列计数表）

将 FASTQ 转为序列计数表，每行一条唯一序列及其出现次数：

```tsv
sequencount
UUUAAAGCUAUAGUAGUAGUACCCGUACAUCUACUAUAUACUAUAUACUAUAUACU  125
AGCUAGCUAGCUAGCUAGCUA  98
UCGUACGUCAGCUAGCUAGCU  67
...
```

#### 3. mirna.count.tsv（miRNA 比对结果）

每个参考 miRNA 的精确匹配和 poly-A 尾检测计数：

```tsv
Name    exact  polyA
ath-miR156a-5p  120  3
ath-miR166b     45   1
ath-miR167a     89   2
ath-miR172a     67   0
...
```

**列说明：**

| 列名 | 说明 |
|------|------|
| Name | miRNA 名称（来自参考 FASTA） |
| exact | 精确匹配的序列计数 |
| polyA | 检测到 poly-A 尾的序列计数（无 Nm 修饰） |

#### 4. mirna.{样本}.annotated.tsv（带样本名的注释文件）

在 mirna.count.tsv 基础上添加样本名前缀，便于合并：

```tsv
Name    U1  U1_polyA
ath-miR156a-5p  120  3
ath-miR166b     45   1
ath-miR167a     89   2
```

#### 5. mirna.merged.tsv（合并总表）

所有样本横向合并的结果表：

```tsv
Name            U1    U1_polyA  U2    U2_polyA  U3    U3_polyA
ath-miR156a-5p  120   3         98    5         145   12
ath-miR166b     45    1         67    2         78    3
ath-miR167a     89    2         102   1         156   8
ath-miR172a     67    0         89    0         123   2
ath-miR393a     23    5         34    8         45    15
...
```

#### 6. 质控输出文件

- `*_qc_metrics.json`：质控指标汇总（JSON 格式）
- `*_qc_report.html`：详细质控 HTML 报告（可在浏览器中打开查看）
- `*_quality_by_pos.json`：每个位点的质量分布数据
- `*_length_dist.json`：read 长度分布数据
- `multiqc_report.html`：所有样本的汇总质控报告
- `multiqc_data/`：质控原始数据目录

### Nm 修饰比例计算示例

基于 `mirna.merged.tsv` 中的数据，可以计算 Nm 修饰比例：

**示例数据：**

```tsv
Name            U1    U1_polyA  U2    U2_polyA
ath-miR156a-5p  120   3         98    5
ath-miR166b     45    1         67    2
```

**计算方法：**

对于 U1 样本中的 ath-miR156a-5p：

- 总读数 = 120
- poly-A 尾读数 = 3（无 Nm 修饰）
- 疑似 Nm 修饰 = 120 - 3 = 117
- **Nm 修饰比例 = (120 - 3) / 120 = 97.5%**

**Python 计算脚本示例：**

```python
import pandas as pd

# 读取合并结果
df = pd.read_csv('results/mirna.merged.tsv', sep='\t')

# 计算每个样本的 Nm 修饰比例
for sample in ['U1', 'U2']:
    total_col = sample
    polya_col = f'{sample}_polyA'
    nm_ratio_col = f'{sample}_Nm_ratio'
    
    df[nm_ratio_col] = (df[total_col] - df[polya_col]) / df[total_col]
    df[nm_ratio_col] = df[nm_ratio_col].fillna(0)

# 保存结果
df.to_csv('results/mirna.merged.with_nm_ratio.tsv', sep='\t', index=False)

print(df[['Name', 'U1', 'U1_polyA', 'U1_Nm_ratio']].head())
```

**输出示例：**

```tsv
Name            U1    U1_polyA  U1_Nm_ratio
ath-miR156a-5p  120   3         0.975
ath-miR166b     45    1         0.978
ath-miR167a     89    2         0.978
```

### 完整运行示例

#### 1. 准备样本列表

```csv
sample_id,R1,R2
WT_rep1,/data/WT_rep1_1.fq.gz,/data/WT_rep1_2.fq.gz
WT_rep2,/data/WT_rep2_1.fq.gz,/data/WT_rep2_2.fq.gz
KO_rep1,/data/KO_rep1_1.fq.gz,/data/KO_rep1_2.fq.gz
KO_rep2,/data/KO_rep2_1.fq.gz,/data/KO_rep2_2.fq.gz
```

#### 2. 运行流程

```bash
cd /path/to/snipr-nm

nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/用户名/srna_results \
    -resume
```

#### 3. 查看运行日志

```log
N E X T F L O W  ~  version 23.04.0
Launching `SNIPR-Nm.nf` [sad_ardinghelli] DSL2 - revision: 8e3a7c3a8
------------------------------------------------------
executor >  local (6)
[aa/bc1234] process > CUTADAPT (WT_rep1)       [100%] 1 of 1, CPU: 8
[ab/cd5678] process > PEAR (WT_rep1)           [100%] 1 of 1, CPU: 4
[ac/ef9012] process > FASTQ2COUNT (WT_rep1)   [100%] 1 of 1
[ad/fg3456] process > MIRNA_COUNT (WT_rep1)   [100%] 1 of 1
[ae/gh7890] process > ANNOTATE_SAMPLE (WT_rep1) [100%] 1 of 1
[ef/ij1234] process > MERGE_SAMPLES           [100%] 1 of 1
[fg/kl5678] process > QC_SAMPLE (WT_rep1)     [100%] 1 of 1
[gh/mn9012] process > MULTIQC                 [100%] 1 of 1
==================================================
SNIPR-Nm run complete
Results : /home/用户名/srna_results
Merged  : /home/用户名/srna_results/mirna.merged.tsv
QC html : /home/用户名/srna_results/multiqc_report.html
==================================================
```

#### 4. 检查输出文件

```bash
# 查看输出目录结构
ls -la /home/用户名/srna_results/

# 查看合并结果
head -n 10 /home/用户名/srna_results/mirna.merged.tsv

# 查看质控报告
open /home/用户名/srna_results/multiqc_report.html
```

---

## 质控指标

`snipr_qc.py` 在 clean R1 上计算以下指标：

| 指标                   | 含义                        | 参考阈值      |
|----------------------|---------------------------|-----------|
| **miRNA 18–25nt 占比** | 文库中 miRNA 典型长度的比例         | > 30% 为良好 |
| **Adapter 污染率**      | 3' 端仍有接头的 reads           | < 5% 为良好  |
| **低复杂度 reads**       | 重复序列比例高（可能是 rRNA/tRNA 污染） | < 20%     |
| **重复率**              | 重复序列比例                    | 视文库复杂度而异  |
| **平均质量**             | Phred 分数                  | Q30+ 为良好  |

> ⚠️ 阈值仅供参考，实际应结合文库构建方案和实验背景判断。

---

## 常见问题

### Q：报错 `No such file or directory`
检查 `samplesheet.csv` 里的路径是否为**完整绝对路径**，文件是否真实存在：
```bash
ls -la /完整路径/U1_1.fq.gz
```

### Q：质控数据全是 0 或报告无数据
确认 `--outdir` 使用**绝对路径**：
```bash
--outdir /home/用户名/srna_results
```
相对路径（`./results`）在某些情况下会导致 Python 脚本无法找到 JSON 文件。

### Q：想增加更多样本
直接在 `samplesheet.csv` 中追加行，pipeline 无需任何修改。

### Q：某步骤报错了怎么查原因
```bash
nextflow run SNIPR-Nm.nf ... -resume
# 查看 .nextflow.log 或对应 work 目录下的 .command.log
```

---

## 参考

Tang et al. "A high-throughput method SNIPR-Nm reveals small RNA methylation ratio." bioRxiv (2025). doi: 10.64898/2025.12.11.693694
