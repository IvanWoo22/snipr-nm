# SNIPR-Nm

SNIPR-Nm 使用 Cutadapt 去接头、PEAR 合并双端 reads，再按参考序列统计 **exact** 和 **polyA** 两类读数。内置 QC 汇总原始 R1/R2、去接头 R1/R2 和合并 reads 的数量、长度及质量。流程需要明确提供参考 FASTA，不预设物种。

## 运行环境

本次已验证：Nextflow 25.10.4、Java 21、Cutadapt 5.2、PEAR 0.9.11、Python 3.12、系统 Perl、gzip。Python 脚本仅依赖标准库；Perl 使用 `IO::Uncompress::Gunzip`。不再需要 Biopython、NumPy、tsv-join 或在线图表库。

本机已准备好的环境可直接使用：

```bash
cd /Users/ivanwoo/snipr-nm
source analysis/20260930/env.sh
nextflow run SNIPR-Nm.nf --help
```

本机使用 `scripts/nextflow-local.sh` 启动既有 JAR，补齐 Java 21 模块参数，已验证任务缓存恢复；不要直接使用简化的 `java -jar` 启动方式。正式运行与续跑需使用不同的 trace 文件名。

其他环境请将上述工具加入 PATH，或指定 `--cutadapt_bin`、`--pear_bin`、`--perl_bin`、`--py_bin` 的可执行文件路径。这些参数是文件名/路径，不是附带参数的 shell 命令。工具升级后应使用新工作目录或关闭缓存重新验证；Nextflow 不会根据 PATH 中程序的版本自动识别软件升级。

## 输入

CSV 必须包含 `sample_id,R1,R2`，可额外提供 `R1_md5,R2_md5`，对应压缩文件本身的 MD5。相对路径以 **CSV 所在目录** 为基准；含逗号的路径须按 CSV 标准加引号。单个 FASTQ 不可重复用于多个样本或两个 mate。样本 ID 仅允许字母、数字、点、下划线、连字符，须以字母或数字开头；`Name` 及以 `_polyA` 结尾的 ID 保留给输出列名。

```csv
sample_id,R1,R2
UM-101,/path/to/UM-101_R1.fastq.gz,/path/to/UM-101_R2.fastq.gz
UM-102,/path/to/UM-102_R1.fastq.gz,/path/to/UM-102_R2.fastq.gz
```

FASTQ 支持普通文件和 gzip，要求标准四行记录、Phred+33、ACGTN 序列。Cutadapt 检查配对记录；后续计数/QC 检查格式和数量守恒。输入检查记录 SHA256、MD5、大小和修改时间，Cutadapt 前后再次校验压缩文件内容。提供预期 MD5 才能将交付方校验值与当前文件比较；未提供时，哈希只能证明流程期间的文件身份，不能证明上游传输已完成。必须等待上传/复制结束后再运行。

参考 FASTA 支持多行序列、大小写及 U→T 标准化；序列须为 ACGTN 且至少 4 nt。保留完整 header、顺序、全部条目和零计数，包括外参。重复名称报错；同序列不同名称分别输出，因此**不能将所有参考行相加当作去重后的匹配 reads 数**。本流程不会自动从混合参考中筛选 miRNA；要只分析 miRNA 和外参，须输入相应参考子集。

## 运行和资源

```bash
nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/hsa_miRNA_external.fa \
  --outdir results \
  -work-dir work \
  -with-trace trace.tsv

# 保留同一 work 目录，增加 -resume 即可续跑
nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/hsa_miRNA_external.fa \
  --outdir results -work-dir work -resume
```

默认适用于当前 8 核 / 16 GiB 本机：本地调度预算 8 CPU、12 GB、最多 2 个任务；Cutadapt 和 PEAR 各用 4 线程。资源值是 Nextflow 的调度预约，不是操作系统强制限额。不要同时启动多个全量流程并假定它们共享该预算。集群需提供对应 executor 配置，本仓库没有预置 LSF 配置。

| 参数 | 默认 | 含义 |
|---|---|---|
| `--samplesheet` | `samplesheet.csv` | 输入样本表 |
| `--mirna_fa` | 必填 | miRNA/外参参考 FASTA |
| `--outdir` | `./results` | 发布目录 |
| `--adaptor_r1` | `AGATCGGAAGAGCACA` | R1 的 3′ 接头 |
| `--adaptor_r2` | `GATCGTCGGACTGTAGAA` | R2 的 3′ 接头 |
| `--cutadapt_opts` | `-O 5 -m 17 -e 0.18 --discard-untrimmed` | 原有去接头规则；不额外做质量修剪 |
| `--pear_min_ovlp` | `17` | PEAR **`-v`**，最小重叠长度 |
| `--pear_min_length` | `17` | PEAR **`-n`**，最小合并长度 |
| `--pear_memory` | `1G` | PEAR `-y` 的工作内存参数 |
| `--cutadapt_threads` / `--pear_threads` | `4` / `4` | 两步骤线程数 |
| `--qc_sample_reads` | `100000` | 各阶段详细 QC 的前 N 条 reads |
| `--local_cpus` / `--local_memory` | `8` / `12 GB` | 本地调度预算 |
| `--local_queue_size` | `2` | 本地并行任务上限 |

`cutadapt_opts` 是受信任用户提供的 shell 选项串，不应来自不可信输入。改变接头、过滤或 PEAR 参数会改变进入计数的 reads，必须保存运行参数并重新验证。

## 计数定义与输出

对每个参考序列：

- **exact (E)**：read 与参考序列完全相同。
- **polyA (P)**：read 以完整参考序列开头，剩余尾巴为 `AA`、`AAA`，或长度大于 3 且 A 比例至少 75%。单个 `A` 不计入 polyA；允许长尾中至多 25% 非 A 碱基。
- 同一参考的 E 和 P 互斥；`classified_total = E + P`。
- `no_tail_fraction = E / (E + P)`；分母为 0 时为 `NA`。它是分类读数中的无尾比例，**并非经过实验标定的绝对 Nm 比例**，也不是 `(E-P)/E`。

例如 E=120、P=3，classified_total=123，无尾比例为 120/123≈0.97561。同序列别名及前缀重叠的参考可共享 reads，流程保留这种历史计数语义。

```text
results/
├── reference_annotations.tsv          # 所有参考名称及序列，含外参
├── mirna.merged.tsv                   # Name, S1, S1_polyA, S2, S2_polyA, ...
├── mirna.exact.tsv                    # Name, sequence, 各样本 E
├── mirna.polyA.tsv                    # Name, sequence, 各样本 P
├── mirna.classified_total.tsv         # Name, sequence, 各样本 E+P
├── mirna.no_tail_fraction.tsv         # Name, sequence, 各样本 E/(E+P)
├── snipr_qc_report.html               # 内置 SNIPR 汇总报告，可离线查看
├── snipr_qc_data/                     # 样本汇总 TSV、逐阶段 JSON、样本报告
└── S1/
    ├── input_manifest.json
    ├── cutadapt.log / cutadapt.json
    ├── pear.log / pear.json           # PEAR 日志及中间文件压缩校验清单
    ├── reads.count.tsv.gz             # sequence/count，无表头
    ├── mirna.count.tsv                # Name/exact/polyA，有表头
    ├── mirna.S1.annotated.tsv          # Name/S1/S1_polyA
    ├── S1_qc_metrics.json
    ├── S1_quality_by_pos.json
    ├── S1_length_dist.json
    └── S1_qc_report.html
```

所有矩阵遵循 samplesheet 的样本顺序及 FASTA 的参考顺序。汇总只读 Nextflow 已暂存的文件，缺少样本或参考不一致会报错。单独运行 `perl bin/mirna_count.pl counts.tsv[.gz] reference.fa` 仍输出旧式无表头 TSV；加 `--header` 则输出 Name/exact/polyA。`fastq2count.pl` 接受 stdin 或一个 FASTQ 路径。

## QC、磁盘与恢复边界

每个样本对 raw R1/R2、clean R1/R2、merged 五个阶段做全文件 reads、bases、长度分布、平均质量、Q30 和 GC 统计。18–25 nt 比例的分母是该文件**所有 reads**，大于 50 nt 的 reads 也保留。位置质量、序列冗余、二核苷酸熵及完整接头 motif 使用各阶段前 N 条 reads，明确记录样本量，不能当作随机抽样估计。序列冗余不是 PCR duplication；完整 motif 命中率不等于部分接头残留率。空分母输出 null/NA。QC 不会过滤 reads。

校验 raw 配对数量、Cutadapt 输入/输出、PEAR 合并/未合并/丢弃的数量以及序列计数总和。去接头后为零的样本保留零计数参考行，跳过会对空输入崩溃的 PEAR。报告是自带的 SNIPR 报告，不是 MultiQC；旧 `multiqc_report.html` 文件名已改为 `snipr_qc_report.html`。

PEAR 的四个 FASTQ 在任务内压缩，每个文件通过解压 SHA256 回验后才移除对应的未压缩副本；压缩版保留在 work 目录供续跑。**PEAR 运行时仍需容纳该样本全部未压缩输出**。工作目录会随样本累计，当前实现不自动删除已完成样本的缓存，也没有磁盘配额/剩余空间闸门。应在正式批量运行前预留空间，并在确认最终输出后由使用者决定是否清理 work；清理后无法复用相应缓存。

深度缓存按输入内容而非仅时间戳判断变更；相关脚本也作为任务输入。其代价是额外读取大文件、计算哈希。失败会明确报告失败，发布目录可能存在部分结果，不能据此认定整批完成；成功须以 Nextflow 正常退出和全部 QC/汇总完成为准。建议每次参数或参考版本变更使用新输出目录，避免旧结果混杂。

计数阶段内存仍随唯一序列数增加；QC 的序列抽样内存随 N 和读长增加，但不再随全量 reads 数增加。五阶段全文件 QC 增加 I/O 和 CPU 开销。macOS 无容器的 Nextflow trace 不提供可靠的 CPU/RSS 指标，不能把空值当成低占用。

## 验证及本批数据的历史关系

```bash
python3 -m unittest discover -s tests -v
python3 tests/run_nextflow_smoke.py --evidence /path/to/smoke-evidence
```

单元回归覆盖尾巴边界、别名、外参、零计数、输入损坏、QC 分母、压缩回验和严格汇总；实际 Nextflow 测试覆盖正常运行、全任务缓存续跑、空文库及错误 MD5。

20260930 批次既有 `analysis/20260930/hsa_primary/` 是此前独立执行器处理 **21 个样本全量数据**的结果，本次代码优化不会把那些结果改称为新版 Nextflow 重跑结果。历史源文件已按原 SHA256 清单归档到 `analysis/20260930/provenance/baseline_source/`。新版根目录 Perl 计数器对 21 个全量序列计数表、每样本 2664 个 hsa miRNA/外参条目逐项核对；真实 FASTQ 验证使用 UM-101 和 UM-121 各前 10 万对 reads。证据位于 `analysis/20260930/validation/optimized_pipeline/`。
