# SNIPR-Nm — 小 RNA-seq 自动化分析流程

> 适用场景：植物小 RNA 文库（miRNA 为主）的下机数据分析。
> 包含：接头修剪 → PEAR 合并 → miRNA 计数 → poly-A 检测 → 质控报告

**前置要求**：会打开终端、了解什么是 FASTQ 文件、知道 conda 是什么。

---

## 目录结构

```
snipr-nm/
├── SNIPR-Nm.nf          # 流程主文件（一般不需要动）
├── nextflow.config      # 配置文件（调整参数、路径用这个）
├── samplesheet.csv      # 样本列表（你需要填写的内容）
├── ath_mirna_neo.fa     # miRNA 参考序列（你自己的）
├── bin/
│   ├── fastq2count.pl   # 内部脚本：将 FASTQ 转为计数表
│   ├── mirna_count.pl   # 内部脚本：miRNA 比对计数
│   ├── snipr_qc.py      # 内部脚本：质控分析
│   └── aggregate_qc.py  # 内部脚本：汇总所有样本的质控报告
└── README.md
```

---

## 第一步：安装依赖软件

### 方法 A：用 conda（推荐）

```bash
# 创建专用环境
conda create -n snipr -c bioconda -c conda-forge \
    cutadapt pear trnctk-typewise dammit

# 额外依赖（tsv-join 和 biopython）
conda install -n snipr -c bioconda -c conda-forge \
    python=3.11 biopython ega Tools4PcrSims

# 确认安装成功
conda activate snipr
which cutadapt pear python3 perl
```

**各软件作用**（只需了解）：
| 软件 | 用途 |
|------|------|
| `cutadapt` | 去除接头序列 |
| `pear` | 合并成对的 PE reads |
| `perl` | 运行 miRNA 计数脚本 |
| `python3` | 运行质控脚本（需要 biopython） |

### 方法 B：直接用系统已有的

如果系统已经装好了 `cutadapt`、`pear`、`perl`、`python3` 且有 biopython，**可以跳过 conda**，直接进入第二步。运行 `which cutadapt` 等命令确认能找到就行。

---

## 第二步：准备样本列表

用文本编辑器打开 `samplesheet.csv`，填入你自己的数据路径：

```csv
sample_id,R1,R2
U1,/home/你的用户名/data/U1_1.fq.gz,/home/你的用户名/data/U1_2.fq.gz
U2,/home/你的用户名/data/U2_1.fq.gz,/home/你的用户名/data/U2_2.fq.gz
```

**格式说明**：
- `sample_id`：样本名，自己起（不要用中文和特殊字符）
- `R1`：第一个 FASTQ 文件的**完整绝对路径**
- `R2`：第二个 FASTQ 文件的**完整绝对路径**
- 多个样本就加多行

**如何找到绝对路径？**
```bash
# 在终端进入你的数据目录，运行：
pwd          # 打印当前目录（绝对路径）
ls *.fq.gz   # 确认文件名
```

---

## 第三步：确认软件路径

如果 conda 环境叫 `snipr`，运行：

```bash
conda activate snipr

which cutadapt   # 输出类似 /home/xxx/miniconda3/envs/snipr/bin/cutadapt
which pear       # 输出类似 /home/xxx/miniconda3/envs/snipr/bin/pear
which python3    # 输出类似 /home/xxx/miniconda3/envs/snipr/bin/python3
which perl       # 输出 /usr/bin/perl
```

把这些路径填到 `nextflow.config` 对应的行（`params.cutadapt_bin = "..."` 等）。

如果这些命令都在 PATH 里（能 which 到），直接用命令名（如 `"cutadapt"`）即可，不需要填完整路径。

---

## 第四步：运行流程

```bash
cd /path/to/snipr-nm   # 进入 pipeline 目录

nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results
```

**断点续跑**（中断后恢复）：
```bash
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    -resume
```

`-resume` 会跳过已经完成步骤，只补跑中断前的部分，速度快很多。

---

## 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--samplesheet` | `samplesheet.csv` | 样本列表 |
| `--mirna_fa` | `ath_mirna_neo.fa` | miRNA 参考序列 FASTA |
| `--outdir` | `./results` | 输出目录（建议用绝对路径） |
| `--adaptor_r1` | `AGATCGGAAGAGCACA` | R1 接头序列 |
| `--adaptor_r2` | `GATCGTCGGACTGTAGAA` | R2 接头序列 |
| `--cutadapt_threads` | `8` | cutadapt 并行线程数 |
| `--pear_threads` | `4` | PEAR 并行线程数 |
| `--pear_min_ovlp` | `17` | PEAR 最小 overlap |

常用覆盖示例：
```bash
# 换接头序列
--adaptor_r1 你的R1接头

# 换线程数
--cutadapt_threads 16 --pear_threads 8

# 输出到指定位置（建议用绝对路径）
--outdir /home/你的用户名/srna_results
```

---

## 输出结果

```
results/
├── U1/
│   ├── cutadapt.log              # 接头去除日志
│   ├── reads.count.tsv           # 序列计数表
│   ├── mirna.count.tsv           # miRNA 比对结果
│   ├── mirna.U1.annotated.tsv   # 加了样本名的注释文件
│   ├── U1_qc_metrics.json       # 质控指标（JSON）
│   ├── U1_qc_report.html         # 详细质控 HTML 报告
│   ├── U1_quality_by_pos.json    # 每位点质量数据
│   └── U1_length_dist.json       # 长度分布数据
├── U2/ ... (每个样本一个目录)
├── mirna.merged.tsv              # 所有样本横向合并总表
├── multiqc_report.html          # 汇总质控报告
└── multiqc_data/                # 质控原始数据
```

**合并表格格式**（`mirna.merged.tsv`）：
```
Name            U1    U1_polyA  U2    U2_polyA  ...
ath-miR156a-5p  120   3         98    5
ath-miR166b     45    1         67    2
ath-miR167a      89   2         102   1
...
```
- 列名带 `_polyA` 的是 poly-A 尾检测结果
- U1、U2 等列是各样本的 miRNA 读数

---

## 质控指标解读

`snipr_qc.py` 在 clean R1 上计算以下指标（每个样本单独一份 HTML 报告）：

| 指标 | 含义 | 参考阈值 |
|------|------|----------|
| **miRNA 18–25nt 占比** | 文库中 miRNA 典型长度的比例 | > 30% 为良好 |
| **Adapter 污染率** | 3' 端仍有接头序列的 reads | < 5% 为良好 |
| **低复杂度 reads** | 重复序列比例高（可能是 rRNA/tRNA 污染）| < 20% |
| **重复率** | 重复序列比例 | 视文库复杂度而异 |
| **平均质量** | Phred 分数 | Q30+（质量好）|

> ⚠️ 阈值仅供参考，实际应结合文库构建方案和实验背景判断。

---

## 常见问题

### Q：报错 `No such file or directory`
检查 `samplesheet.csv` 里的路径是否**完整绝对路径**，以及文件是否真实存在。
```bash
ls -la /你的完整路径/U1_1.fq.gz
```

### Q：报错 `Process X terminated with an error exit status (1)`
某个步骤失败了。看日志里具体是哪个命令出错。常见原因：
- 输入文件格式不对（不是标准 FASTQ）
- 某个软件没装好（`which cutadapt` 检查）
- 路径填写错误

### Q：结果全是 0 或质控数据缺失
确认 `--outdir` 使用**绝对路径**，例如：
```bash
--outdir /home/用户名/srna_results
```
相对路径（`./results`）在某些情况下会导致 Python 脚本找不到文件。

### Q：想只跑部分步骤
用 `-resume` 配合清空对应步骤的 work 目录即可。或者修改 pipeline 注释掉不需要的步骤。

### Q：如何增加更多样本？
直接在 `samplesheet.csv` 加行，不用改 pipeline 代码。

---

## 流程架构（可跳过）

```
输入 FASTQ
    ↓
[Cutadapt]  去除接头序列
    ↓
[PEAR]       合并 PE reads
    ↓
[fastq2count] 转为计数表
    ↓
[mirna_count] 与 miRNA 参考比对计数
    ↓
[annotate]    给每个样本加列名（避免名称冲突）
    ↓
[merge]       所有样本横向合并成一个表
    ↓
[QC]          质控分析（长度分布、质量、复杂度）
    ↓
[multiqc]     汇总所有样本的质控报告
```

---

## 调试积累的经验

1. **Nextflow work dir ≠ 项目目录**  
   pipeline 的每个步骤实际运行在临时 work 目录里，不是项目目录。如果脚本依赖外部文件（JSON、FASTA 等），必须用绝对路径或在 workflow 里提前 resolve 好传给脚本。

2. **process 里的 `input:` 是 Nextflow DSL，不是 shell 变量**  
   在 script block 里的 `${params.outdir}` 是 Groovy 变量。如果要传到 Python 里做路径解析，应该用 `val` 类型传参或在 workflow 里先 resolve 成绝对路径再调用。

3. **Python 脚本自己输出的文件名 ≠ 期望的输出文件名**  
   如果 Python 脚本已经按 `--sample` 参数生成了 `{name}_xxx.json`，不要再用 `mv` 改文件名，容易出现 race condition。

4. **channel tuple 格式必须与 input 声明严格匹配**  
   Nextflow 对类型要求严格。`(val, path)` 的 tuple 传进去时如果只有 `path`，会在运行时直接报错或静默产生空结果。

5. **`./results` 相对路径在 aggregate 时找不到文件**  
   因为 aggregate_qc.py 运行时的工作目录是 Nextflow 的 work dir，而 JSON 文件在项目目录的 results 下。解决：命令行用绝对路径 `--outdir`，并在 workflow 里 resolve 好再传给 Python。
