# SNIPR-Nm

SNIPR-Nm 使用 Cutadapt 去接头、PEAR 合并双端 reads，再按参考序列统计 **exact** 和 **polyA** 两类读数。内置 QC 汇总原始 R1/R2、去接头 R1/R2 和合并 reads 的数量、长度及质量。流程需要明确提供参考 FASTA，不预设物种。

## 运行环境

本次已验证：Nextflow 25.10.4、Java 21、Cutadapt 5.2、PEAR 0.9.11、Python 3.12、系统 Perl、gzip，以及支持 C++17 的编译器和 zlib 开发库。Python 脚本仅依赖标准库；Perl 使用 `IO::Uncompress::Gunzip`。原生扫描器由 Nextflow 的 BUILD_SCANNER 编译并缓存；macOS 可使用 Xcode Command Line Tools，Linux 需要相应编译器和 zlib headers。Python 功能使用 `hashlib.file_digest`，要求 Python 3.11 或更新版本。不再需要 Biopython、NumPy、tsv-join 或在线图表库。

安装上述依赖并激活对应环境后，在项目目录中运行：

```bash
# 确保 Nextflow、Java 和所需工具已加入 PATH。
nextflow run SNIPR-Nm.nf --help
```

若已在 `.runtime/nextflow/` 或 `NXF_HOME` 中安装 Nextflow 25.10.4 JAR，可使用 `scripts/nextflow-local.sh` 启动；该脚本包含 Java 21 模块参数，需要先配置相应的 `JAVA_HOME`，不会安装依赖。正式运行与续跑需使用不同的 trace 文件名。

其他环境请将上述工具加入 PATH，或指定 `--cutadapt_bin`、`--pear_bin`、`--perl_bin`、`--py_bin`、`--cxx_bin` 的可执行文件路径。这些参数是文件名/路径，不是附带参数的 shell 命令。工具升级后应使用新工作目录或关闭缓存重新验证；Nextflow 不会根据 PATH 中程序的版本自动识别软件升级。

## 更新已安装的流程代码

等待正在运行的流程结束后，在原项目目录更新代码。以下命令从 `origin/main` 获取已经发布到远端的版本。先检查本地修改；需要保留的源码改动应先提交到独立分支或另存，再继续更新。本地数据和运行目录无需加入 Git。

```bash
cd /path/to/snipr-nm
# 先检查是否有未保存的源码改动。
git status --short
git diff --stat
git switch main

# 给当前代码保留一个版本参照，再获取并查看远端更新。
git branch "codex/before-update-$(date +%Y%m%d-%H%M%S)"
git fetch origin
git log --oneline HEAD..origin/main

# 只允许快进更新；本地与远端分叉时会停止。
git merge --ff-only origin/main
git log -1 --oneline
nextflow run SNIPR-Nm.nf --help
```

如果 Git 提示本地修改冲突或分支已经分叉，应先处理这些改动再更新；上述命令不会强行覆盖本地历史，也不删除分析文件。拉取仓库只更新流程代码，不会安装或升级 Nextflow、Python、Cutadapt 或 PEAR。此版本需要 Python 3.11+、C++17 编译器及 zlib headers；原生扫描器由 Nextflow 自动编译和缓存。

更新后，若继续使用原工作目录，可执行：

```bash
nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/reference.fa \
  --outdir results -work-dir work \
  -resume -with-trace trace.updated.tsv
```

续跑应在原启动目录执行，以便找到原 `.nextflow/` 缓存，同时保留原输入、结果/work 目录、`-profile`、`-c` 配置及分析参数，包括 `--intermediates`。上面的命令演示使用默认参数的情形。流程代码和依赖改变后，相关任务可能重新执行，不能保证所有旧任务都命中缓存；应保留 `work/` 与 `.nextflow/`。若要独立开始一轮分析，指定新的 outdir/work 目录并省略 `-resume`。

## 输入

仓库自带的 `samplesheet.csv` 是格式模板，运行前必须将占位路径替换为实际输入文件。

CSV 必须包含 `sample_id,R1,R2`，可额外提供 `R1_md5,R2_md5`，对应压缩文件本身的 MD5。相对路径以 **CSV 所在目录** 为基准；含逗号的路径须按 CSV 标准加引号。单个 FASTQ 不可重复用于多个样本或两个 mate。样本 ID 仅允许字母、数字、点、下划线、连字符，须以字母或数字开头；`Name` 及以 `_polyA` 结尾的 ID 保留给输出列名。

```csv
sample_id,R1,R2
S1,/path/to/S1_R1.fastq.gz,/path/to/S1_R2.fastq.gz
S2,/path/to/S2_R1.fastq.gz,/path/to/S2_R2.fastq.gz
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

默认 `desktop` 配置使用 6 CPU、10 GB、6 个任务槽；单线程扫描每任务 1 CPU，各扫描 process 最多 4 个并发；精确计数最多 2 个并发。Cutadapt/PEAR 各用 4 线程，并受总预算限制。`-profile workstation` 使用 8 CPU/12 GB/8 槽；`-profile test` 使用 2 CPU/6 GB、单线程工具和较低并发。资源值是 Nextflow 的调度预约，不是操作系统强制限额。不要同时启动多个全量流程并假定它们共享该预算。集群需提供对应 executor 配置，本仓库没有预置 LSF 配置。

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
| `--local_cpus` / `--local_memory` | `6` / `10 GB` | 本地调度预算 |
| `--local_queue_size` | `6` | 本地并行任务上限 |
| `--intermediates` | `checkpoint` | 步骤级缓存，或 `compact` 样本内临时文件回收 |
| `--counts_sheet` | 空 | 已验证 counts 的样本表；设置后跳过原始数据预处理与 QC |
| `--count_max_unique` | `1000000` | 哈希内唯一序列阈值；达到时精确溢写，不截断 reads |
| `--qc_max_forks` / `--count_max_forks` | `4` / `2` | 各扫描/计数 process 的并发上限 |
| `--preprocess_max_forks` | `2` | 各预处理 process 的并发上限，并非跨阶段合计 |
| `--compact_threads` | `4` | 紧凑任务内顺序执行的 Cutadapt/PEAR 线程数 |
| `--execution_reports` | `false` | 输出 pipeline_info 下的 trace、report 和 timeline |

每轮选择一个资源 profile。命令行 `--local_cpus`、`--local_memory`、`--local_queue_size` 及工具线程参数可覆盖其默认值。若在独立的 `-c` 文件中覆盖资源，应直接设置 `executor.cpus`、`executor.memory`、`executor.queueSize` 和相应 `process` selector；只在后加载文件中修改 `params.*`，不会重新计算已经求值的配置项。参见 [Nextflow 配置规则](https://docs.seqera.io/nextflow/config)。

两个模式都将 `cutadapt_opts` 解析为参数列表，支持带引号的参数值；不执行 shell 命令替换、变量展开、管道或重定向。改变接头、过滤或 PEAR 参数会改变进入计数的 reads，必须保存运行参数并重新验证。

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
    ├── pear.log / pear.json           # pear.json 仅 checkpoint；compact 输出 intermediates.json
    ├── reads.count.tsv.gz             # sequence/count，无表头
    ├── counts.manifest.json           # 可导入检查点、内容摘要、输入/参数与 QC 证据
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

原生扫描器按阶段计算指标；merged 序列计数和 QC 只解析一次 FASTQ。`ASSEMBLE_QC` 只读取小型指标文件，`RENDER_QC` 单独生成报告，不再在尾部重扫五个 FASTQ。各 process 仅输入实际依赖脚本：更改报告渲染不会重跑预处理或数值 QC，更改参考不重跑预处理，更改一个样本不重跑其他样本。checkpoint 模式下修改 QC 采样数会重算指标/计数，但保留 Cutadapt/PEAR 缓存。

`--intermediates checkpoint` 保留步骤级恢复：PEAR 输出压缩并回验后留在 work，磁盘用量随样本累计。`--intermediates compact` 把 clean/merged 生成、指标和计数放在一个样本任务中：验证 counts、QC 与检查点后，只删除本任务固定名单的六类临时 FASTQ，记录 `intermediates.json`；失败保留现场。紧凑模式保留 raw QC 的独立缓存，样本预处理失败则重做该组合任务；修改其 QC 参数也会重跑组合任务。它不再压缩随即就会被删除的 PEAR 大文件。两种模式都不删除原始输入。**PEAR 仍需容纳单样本未压缩输出；紧凑模式不是零临时空间，也没有整机磁盘配额保证。**

从 FASTQ 新建一轮紧凑模式分析，并保存执行报告：

```bash
nextflow run SNIPR-Nm.nf --samplesheet samplesheet.csv \
  --mirna_fa reference.fa --intermediates compact --outdir compact-results \
  -work-dir compact-work --execution_reports true
```

精确计数在唯一序列达到 `--count_max_unique` 时写出排序分片，最多 32 路归并，计数及排序与内存路径相同；记录分片数、最大内存内唯一序列数和累计分片字节数。该阈值是唯一序列数预算，不是严格 RSS 字节上限；实际内存还受读长、有界 QC 表和解析缓冲区影响。溢写使用任务工作目录，需预留磁盘。计数压缩和 SHA256 校验也采用流式读取，不把整个计数文件读入内存。QC 内存受采样数和读长影响，不随全量 reads 数线性增长。

深度缓存仍检查内容身份。成功退出并完成全部汇总才表示批次完成，失败目录可能有部分发布结果；发布对应输出显式覆盖，防止恢复缓存后结果目录残留不同参数的旧文件。删除样本或更换输入集合时仍推荐新 outdir，因为不相关旧样本目录不会自动删除。工作目录和 `.nextflow/cache` 共同支持步骤级恢复；不要在活动任务中删除它们。工具升级需重新验证，PATH 中工具版本不会自动纳入缓存。

macOS 无容器的 Nextflow trace CPU/RSS 可能缺失，预约内存不代表实测峰值。`--execution_reports true` 生成 `pipeline_info/trace.tsv`、`report.html`、`timeline.html`（每轮覆盖）；如需保留各轮历史，使用新的 outdir 或显式 `-with-trace/-with-report/-with-timeline` 路径。

## 从 counts 检查点重新计数

```csv
sample_id,counts,manifest
S1,/results/S1/reads.count.tsv.gz,/results/S1/counts.manifest.json
S2,/results/S2/reads.count.tsv.gz,/results/S2/counts.manifest.json
```

```bash
nextflow run SNIPR-Nm.nf --counts_sheet checkpoints.csv \
  --mirna_fa updated_reference.fa --outdir recount -work-dir recount-work
```

导入会核验 schema、样本、压缩和解压内容 SHA256、排序唯一性、总读数、唯一序列数及 QC 完成证据，随后只做参考计数/矩阵，不读取原始 FASTQ，也不重新生成 QC。输入为 `snipr-counts-v1` 清单，旧 counts 文件不能只改文件名冒充检查点。counts 不保留质量、read ID 或配对信息，原始数据仍需保存。更换参考无需旧 work 缓存。

原生扫描器将全文件逐碱基校验、GC、质量总和和 Q30 合并为一遍统计，并采用缓冲区行扫描。全量和 first-N 统计定义保持不变。原有 Python QC 和 Perl 计数入口继续保留用于兼容。
