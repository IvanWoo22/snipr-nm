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

| 参数                   | 默认值                  | 说明                 |
|----------------------|----------------------|--------------------|
| `--samplesheet`      | `samplesheet.csv`    | 样本列表               |
| `--mirna_fa`         | `ath_mirna_neo.fa`   | miRNA 参考序列 FASTA   |
| `--outdir`           | `./results`          | 输出目录（建议用绝对路径）      |
| `--adaptor_r1`       | `AGATCGGAAGAGCACA`   | R1 接头序列            |
| `--adaptor_r2`       | `GATCGTCGGACTGTAGAA` | R2 接头序列            |
| `--cutadapt_threads` | `8`                  | cutadapt 并行线程数     |
| `--pear_threads`     | `4`                  | PEAR 并行线程数         |
| `--pear_min_ovlp`    | `17`                 | PEAR 最小 overlap 长度 |

常用覆盖示例：
```bash
# 使用绝对路径指定输出位置
--outdir /home/用户名/srna_results

# 更换接头序列
--adaptor_r1 你的R1接头序列

# 调整线程数
--cutadapt_threads 16 --pear_threads 8
```

---

## 输出结果

```
results/
├── U1/
│   ├── cutadapt.log              # 接头去除日志
│   ├── reads.count.tsv           # 序列计数表
│   ├── mirna.count.tsv           # miRNA 比对结果
│   ├── mirna.U1.annotated.tsv   # 带样本名的注释文件
│   ├── U1_qc_metrics.json       # 质控指标（JSON）
│   ├── U1_qc_report.html         # 详细质控 HTML 报告
│   ├── U1_quality_by_pos.json    # 每位点质量数据
│   └── U1_length_dist.json       # 长度分布数据
├── U2/ ... (每个样本一个子目录)
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
- 每列对应一个样本的读数
- `{样本}_polyA` 列：该样本中检测到 poly-A 尾的序列数（即**无 Nm 修饰**的分子）
- 总读数减去 poly-A 尾读数 = **疑似 Nm 修饰**的分子
- Nm 修饰比例 = (总读数 - polyA 读数) / 总读数

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

## 流程架构

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
