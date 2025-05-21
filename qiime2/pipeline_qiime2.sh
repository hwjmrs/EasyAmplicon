[TOC]

# QIIME2 2025.4分析流程

## 0. 软件安装(附录1)

    # 在Linux、Mac、Windows内置Linux子系统(支持右键粘贴)下安装并运行
    # 详细教程参阅官网 https://docs.qiime2.org/2025.4/
    # 安装Windows子系统，https://mp.weixin.qq.com/s/0PfA0bqdvrEbo62zPVq4kQ

## 1. 准备工作

    # 设置工作目录，如服务器为~/amplicon/qiime2，Win子系统如下：
    wd=/mnt/d/A250412/qiime2/
    # 进入工作目录
    mkdir -p ${wd}
    cd ${wd}
    # 激活QIIME2工作环境，旧版conda使用source替换conda运行
    conda activate qiime2-amplicon-2025.4
    
    # 准备样本元数据metadata.txt、原始数据seq/*.fq.gz
    
    ## 此处代码基于metadata.txt从公共数据下载测序数据，按GSA的CRA(批次)和CRR(样品)编号下载数据
    mkdir -p seq
    # # 公网下载
    # cd seq
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O "$1"_1.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    #  awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O "$1"_2.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    # cd .. && ls -lsh seq
    # 从其他地方链接(不额外占用空间)
    ln /mnt/d/A250412/seq/* seq/
    ln /mnt/d/A250412/result/metadata.txt ./
    # 根据metadata生成manifest文件
    awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
      NR>1{print $1"\t$PWD/seq/"$1"_1.fq.gz\t$PWD/seq/"$1"_2.fq.gz"}' \
      metadata.txt > manifest
    head -n3 manifest
    
    # 数据导入qiime2，格式为双端33格式
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path manifest \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33V2
    # fq文件1G用时7m，fq.gz压缩格式仅34s，测试数据9s


## 2. 生成特征表和代表序列

### 方法1. DADA2(慢，依赖R和Python包多容易报错)

    # 支持多线程加速，90万条PE250数据，0/96p, 34m；24p, 44m；8p, 77m；1p, 462m
    # 27万，8p, 9m; 4p, 11m;
    time qiime dada2 denoise-paired \
      --i-demultiplexed-seqs demux.qza \
      --p-n-threads 4 \
      --p-trim-left-f 29 --p-trim-left-r 18 \
      --p-trunc-len-f 0 --p-trunc-len-r 0 \
      --o-table dada2-table.qza \
      --o-representative-sequences dada2-rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza
    # 确定使用dada2结果并导入主流程
    cp dada2-table.qza table.qza
    cp dada2-rep-seqs.qza rep-seqs.qza

### 方法2. 外部导入特征表和代表序列(常用)

    # 上传我们生成的OTU表otutab.txt和代表序列otus.fa
    # 转换文本为Biom1.0，注意biom --version 2.1.5/8可以，2.1.7报错
    biom convert -i otutab.txt -o otutab.biom \
      --table-type="OTU table" --to-json
    # 导入特征表，9s
    qiime tools import --input-path otutab.biom \
      --type 'FeatureTable[Frequency]' --input-format BIOMV100Format \
      --output-path table.qza
    # 导入代表序列，8s
    qiime tools import --input-path otus.fa \
      --type 'FeatureData[Sequence]' \
      --output-path rep-seqs.qza

### 特征表和代表序列统计

    qiime feature-table summarize \
      --i-table table.qza \
      --o-visualization table.qzv \
      --m-sample-metadata-file metadata.txt
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
    # 下载qzv在线https://view.qiime2.org/ 查看，dada2只有1千多个ASV


## 3. Alpha和beta多样性分析

### 构建进化树用于多样性分析 53s

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza

### 计算核心多样性 

    # 13s，采样深度通常选择最小值，来自table.qzv
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny rooted-tree.qza \
      --i-table table.qza \
      --p-sampling-depth 7439 \
      --m-metadata-file metadata.txt \
      --output-dir core-metrics-results

### Alpha多样性组间显著性分析和可视化

    # 7s, 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
    index=observed_features
    qiime diversity alpha-group-significance \
      --i-alpha-diversity core-metrics-results/${index}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization core-metrics-results/${index}-group-significance.qzv

### Alpha多样性稀疏曲线

    # 25s, max-depth选最大值，来自table.qzv
    qiime diversity alpha-rarefaction \
      --i-table table.qza \
      --i-phylogeny rooted-tree.qza \
      --p-max-depth 10298 \
      --m-metadata-file metadata.txt \
      --o-visualization alpha-rarefaction.qzv
    # 结果有observed_otus, shannon, 和faith_pd三种指数可选

### Beta多样性组间显著性分析和可视化

    # 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
    # 7s, 指定分组是减少计算量，置换检验较耗时
    distance=weighted_unifrac
    column=Group
    qiime diversity beta-group-significance \
      --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
      --p-pairwise


## 4. 物种组成分析
    # 训练分类器———全长（通用）
    ## 下载参考序列和物种分类信息
    wget -c --no-check-certificate https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
    wget -c --no-check-certificate https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza
    
    ## 全长通用分类器训练，耗时2小时左右
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 2022.10.backbone.full-length.fna.qza \
      --i-reference-taxonomy 2022.10.backbone.tax.qza \
      --o-classifier classifier-full.qza
    
    # 训练分类器—指定V区分类器
    ## 使用与测试数据对应的V5 (799F) - V7 (1193R) 引物为例进行提取序列，耗时约6分钟
    time qiime feature-classifier extract-reads \
      --i-sequences 2022.10.backbone.full-length.fna.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --p-trunc-len 350 \
      --o-reads ref-seqs.qza 
    
    ## 基于筛选的指定区段，生成实验特异的分类器，耗时约30分钟
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy 2022.10.backbone.tax.qza \
      --o-classifier classifier_greengenes_V5-V7.qza
    
    # 物种注释
    time qiime feature-classifier classify-sklearn \
      --i-classifier classifier_greengenes_V5-V7.qza \
      --i-reads rep-seqs.qza \
      --o-classification taxonomy.qza
    
    # 可视化物种注释
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
    
    # 堆叠柱状图展示
    qiime taxa barplot \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization taxa-bar-plots.qzv


## 5. 差异分析ancom

    # 格式化特征表，添加伪计数，4s
    qiime composition add-pseudocount \
      --i-table table.qza \
      --o-composition-table comp-table.qza
    
    # 计算差异特征，指定分组类型比较，1m
    column=Group
    time qiime composition ancom \
      --i-table comp-table.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-${column}.qzv
    
    # 按属水平合并，并统计
    ## 按属水平合并，6s
    qiime taxa collapse \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table table-l6.qza
    # 格式化特征表，添加伪计数，6s
    qiime composition add-pseudocount \
      --i-table table-l6.qza \
      --o-composition-table comp-table-l6.qza
    # 计算差异属，指定分组类型比较，16s
    qiime composition ancom \
      --i-table comp-table-l6.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-l6-${column}.qzv


# 附录

    wd=/mnt/c/amplicon/qiime2
    mkdir -p $wd
    cd $wd
    
## 1. qiime2 2025.4安装

### 安装Conda

    # 下载、安装和启动conda
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    ~/miniconda3/condabin/conda init
    # 关闭终端重新打开

### 方法1. Conda在线安装QIIME

    # 附软件在线安装和打包代码
    n=qiime2-amplicon-2025.4
    # 下载软件列表
    wget -c http://www.imeta.science/db/qiime2/qiime2-amplicon-ubuntu-latest-conda.yml
    # 新环境安装，可在不同电脑服务器上安装成功后打包分发
    conda env create -n ${n} --file qiime2-amplicon-ubuntu-latest-conda.yml
    # 环境打包(可选，1.2G)
    conda pack -n ${n} -o ${n}.tar.gz

### 方法2. 本地安装QIIME

    n=qiime2-amplicon-2025.4
    # 安装包下载链接 
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 新环境安装
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # 激活并初始化环境
    conda activate ${n}
    conda unpack

## q2-picrust2插件教程（功能预测）

### 1.安装
     # 创建学习目录
     mkdir -p q2-picrust2-tutorial
     cd q2-picrust2-tutorial
     
     # 使用conda安装q2-picrust2（原始链接）
     conda env create --name q2-picrust2-amplicon-2024.5 --file \
     https://raw.githubusercontent.com/picrust/q2-picrust2/refs/heads/master/environment-files/q2-picrust2-qiime2-amplicon-2024.5.yml

     # 备用链接
     wget -c http://www.imeta.science/db/qiime2/q2-picrust2-amplicon-2024.5.yml
     conda env create --name q2-picrust2-amplicon-2024.5 --file q2-picrust2-amplicon-2024.5.yml

     # 激活环境
     conda activate q2-picrust2-amplicon-2024.5
     
### 2.下载PICRUSt2 教程中的测试文件
     wget -c http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/mammal_biom.qza
     wget -c http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/mammal_seqs.qza
     wget -c http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/mammal_metadata.tsv
     
### 3.功能预测
   # 创建临时目录并设置环境变量
   mkdir -p ~/tmp/qiime2_picrust2_temp
   export TMPDIR=~/tmp/qiime2_picrust2_temp
   
   # 运行SEPP
   qiime picrust2 full-pipeline \
     --i-table mammal_biom.qza \
     --i-seq mammal_seqs.qza \
     --output-dir q2-picrust2_output \
     --p-placement-tool sepp \
     --p-threads 8 \
     --p-hsp-method pic \
     --p-max-nsti 2 \
     --verbose 
    # 输出结果文件有ec_metagenome.qza-EC 宏基因组预测、ko_metagenome.qza-KO宏基因组预测和pathway_abundance.qza-MetaCyc通路丰度预测
    
###  4.获取有关通路丰度文件的摘要信息
     qiime feature-table summarize \
       --i-table q2-picrust2_output/pathway_abundance.qza \
       --o-visualization q2-picrust2_output/pathway_abundance.qzv

### 5.计算多样性指标
    qiime diversity core-metrics \
      --i-table q2-picrust2_output/pathway_abundance.qza \
      --p-sampling-depth 236867 \
      --m-metadata-file mammal_metadata.tsv \
      --output-dir pathabun_core_metrics_out \
      --p-n-jobs 1
      
## q2-ITSxpress插件教程（修剪 ITS 序列的保守侧翼区域）
### 1.安装
    # 创建目录
    mkdir -p q2-ITSxpress-tutorial
    cd q2-ITSxpress-tutorial
    
    # 激活QIIME 2环境
    conda activate qiime2-amplicon-2025.4
    
    # 使用 Bioconda 安装 ITSxpress
    conda install -c bioconda itsxpress
   
    # 刷新插件
    qiime dev refresh-cache
  
    # 检查是否安装成功
    qiime itsxpress

### 2.下载示例数据
    wget -c http://www.imeta.science/db/qiime2/sample1_r1.fastq.gz
    wget -c http://www.imeta.science/db/qiime2/sample1_r2.fastq.gz
    wget -c http://www.imeta.science/db/qiime2/sample2_r1.fastq.gz
    wget -c http://www.imeta.science/db/qiime2/sample2_r2.fastq.gz
    wget -c http://www.imeta.science/db/qiime2/manifest.txt
    wget -c http://www.imeta.science/db/qiime2/mapping.txt

### 3.导入序列数据
    qiime tools import \
      --type SampleData[PairedEndSequencesWithQuality] \
      --input-format PairedEndFastqManifestPhred33\
      --input-path manifest.txt \
      --output-path sequences.qza
      
    # 通过运行 summarize 命令来查看数据的质量
    qiime demux summarize \
      --i-data sequences.qza \
      --o-visualization sequences.qzv
      
### 4.使用 ITSxpress 修剪 ITS 样本
    qiime itsxpress trim-pair-output-unmerged\
      --i-per-sample-sequences sequences.qza \
      --p-region ITS2 \
      --p-taxa F \
      --p-cluster-id 1.0 \
      --p-threads 16 \
      --o-trimmed trimmed_exact.qza
      
     # 使用以下命令聚类 99.5% 的序列相似度
      qiime itsxpress trim-pair-output-unmerged \
        --i-per-sample-sequences sequences.qza \
        --p-region ITS2 \
        --p-taxa F \
        --p-cluster-id 0.995 \
        --p-threads 16 \
        --o-trimmed trimmed.qza

### 5.使用 DADA2 识别序列变异
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs trimmed_exact.qza \
      --p-trunc-len-r 0 \
      --p-trunc-len-f 0 \
      --output-dir dada2out
      
     # 汇总数据以进行目视检查
      qiime feature-table summarize \
        --i-table dada2out/table.qza \
        --o-visualization tableviz.qzv
### 6.从UNITE下载用于真菌分类的参考数据
    # 原始下载链接
    wget -c https://s3.hpc.ut.ee/plutof-public/original/db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz
    
    # 备用下载链接
    wget -c http://www.imeta.science/db/qiime2/db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz
    
    # 解压
    tar -xzvf db1d6ddb-a35d-48c5-8b1a-ad9dd3310c6d.tgz
    
    # 将最新的UNITE数据导入QIIME 2
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path sh_refs_qiime_ver10_dynamic_04.04.2024.fasta \
      --output-path unite.qza

    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path sh_taxonomy_qiime_ver10_dynamic_04.04.2024.txt \
      --output-path unite-taxonomy.qza
    
### 7.训练QIIME 2分类器，耗时1小时左右
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads unite.qza \
      --i-reference-taxonomy unite-taxonomy.qza \
      --o-classifier classifier.qza
      
### 8.对序列变体进行分类
    qiime feature-classifier classify-sklearn \
      --i-classifier classifier.qza \
      --i-reads dada2out/representative_sequences.qza \
      --o-classification taxonomy.qza

   # 汇总结果
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
      
   # 堆叠柱状图展示
    qiime taxa barplot \
      --i-table dada2out/table.qza  \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file mapping.txt \
      --o-visualization taxa-bar-plots.qzv

## q2-FMT插件教程（评估粪菌移植后的植入范围）
### 1.安装
    # 下载清单文件（原始链接）
    wget -c https://raw.githubusercontent.com/qiime2/q2-fmt/refs/heads/dev/environment-files/q2-fmt-qiime2-amplicon-2024.10.yml

    # 备用链接
    wget -c http://www.imeta.science/db/qiime2/q2-fmt-qiime2-amplicon-2024.10.yml

    # conda新建环境安装
    conda env create --name q2-fmt-amplicon-2024.10 --file q2-fmt-qiime2-amplicon-2024.10.yml

    # 安装过程中遇到提示pip 安装一个依赖项无法从 GitHub 下载 q2-fmt 的压缩包时，手动下载从本地安装
    wget -c http://www.imeta.science/db/qiime2/q2-fmt-2024.11.1.zip
    ~/miniconda3/envs/q2-fmt-amplicon-2024.10/bin/python -m pip install ./q2-fmt-2024.11.1.zip

### 2.下载元数据和特征表
    # 创建目录并进入
    mkdir q2-fmt-tutorial
    cd q2-fmt-tutorial

    # 数据下载
    wget -O 'sample-metadata.tsv' \
    'https://q2-fmt.readthedocs.io/en/latest/data/tutorial/sample-metadata.tsv'
    wget -O 'feature-table.qza' \
   'https://q2-fmt.readthedocs.io/en/latest/data/tutorial/feature-table.qza'

    # 汇总特征表
    qiime feature-table summarize \
      --i-table feature-table.qza \
      --m-sample-metadata-file sample-metadata.tsv \
      --o-visualization autofmt-table-summ.qzv
### 3.计算多样性指标
   qiime diversity core-metrics \
     --i-table feature-table.qza \
     --p-sampling-depth 10000 \
     --m-metadata-file sample-metadata.tsv \
     --output-dir diversity-core-metrics

### 4.绘制云雨图了解受体与供体的距离
   qiime fmt cc \
     --i-diversity-measure diversity-core-metrics/jaccard_distance_matrix.qza \
     --m-metadata-file sample-metadata.tsv \
     --p-distance-to donor \
     --p-compare baseline \
     --p-time-column timepoints \
     --p-reference-column DonorSampleID \
     --p-subject-column PatientID \
     --p-filter-missing-references \
     --p-against-group 0 \
     --p-p-val-approx asymptotic \
     --o-stats jaccard-raincloud-stats.qza \
     --o-raincloud-plot jaccard-raincloud-plot.qzv

## 2. 物种注释数据训练集

### Greengenes2 2022.10 full length sequences

    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza
    # 分类器训练，耗时3小时多
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 2022.10.backbone.full-length.fna.qza \
      --i-reference-taxonomy 2022.10.backbone.tax.qza \
      --o-classifier classifier-gg22-full.qza 
    # 使用与测试数据对应的V5 (799F) - V7 (1193R) 引物
    time qiime feature-classifier extract-reads \
      --i-sequences 2022.10.backbone.full-length.fna.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --p-trunc-len 350 \
      --o-reads ref-seqs.qza
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy 2022.10.backbone.tax.qza \
      --o-classifier classifier_gg22_V5-V7.qza

    # 国内高速备份链接
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/2022.10.backbone.full-length.nb.qza
    # QIIME或Greengene官网下载(国外较慢)
    wget -c https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone_full_length.nb.qza
    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.nb.qza
    # 统一文件名
    mv 2022.10.backbone.full-length.nb.qza gg_2022_10_backbone_full_length.nb.qza 
    
### Silva 138 99% OTUs full-length sequences

    # 官网下载
    wget -c https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza
    wget -c https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza

## 3. 物种注释数据训练集(gg13为例，可选)

    # 下载数据库文件(greengenes, 320M)
    # wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 国内备用链接
    # wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus.tar.gz
    # 国内备份核心99数据库(60M)
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus_99.tar.gz
    mv gg_13_8_otus_99.tar.gz gg_13_8_otus.tar.gz
    # 解压
    tar -zxvf gg_13_8_otus.tar.gz
    
    # 使用rep_set文件中的99_otus.fasta数据和taxonomy中的99_OTU_taxonomy.txt数据作为参考物种注释
    # 导入参考序列，50s
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path gg_13_8_otus/rep_set/99_otus.fasta \
      --output-path 99_otus.qza
    # Fontconfig error: Cannot load default config file 不影响结果
    
    # 导入物种分类信息，6s
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
      --output-path ref-taxonomy.qza
    
    # Train the classifier（训练分类器）——全长，30m
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 99_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99.qza
    # 备用下载链接	wget -c http://bailab.genetics.ac.cn/db/GreenGenes/classifier_gg_13_8_99.qza
      
    # 引物提取参考序列的扩增区段 Extract reference reads
    # 常用Greengenes 13_8 99% OTUs from 341F CCTACGGGNGGCWGCAG/805R GACTACHVGGGTATCTAATCC region of sequences（分类器描述），提供测序的引物序列，截取对应的区域进行比对，达到分类的目的。
    # 本次使用引物799F-1193R，请根据实际替换, 8m
    time qiime feature-classifier extract-reads \
      --i-sequences 99_otus.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --o-reads ref-seqs.qza
    # Train the classifier（训练分类器）
    # 基于筛选的指定区段，生成实验特异的分类器，7 min
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99_V5-V7.qza
    
    # 常见问题1：scikit-learn版本不兼容，重新fit-classifier-naive-bayes构建训练集即可
    Plugin error from feature-classifier:
      The scikit-learn version (0.21.2) used to generate this artifact does not match the current version of scikit-learn installed (0.22.1). Please retrain your classifier for your current deployment to prevent data-corruption errors.
    Debug info has been saved to /tmp/qiime2-q2cli-err-5ngzk2hm.log
    
    # 常见问题2: dada2运行报错，环境配置不完全，运行conda unpack初始化
    Plugin error from dada2:
    An error was encountered while running DADA2 in R (return code 255), please inspect stdout and stderr to learn more.
    Debug info has been saved to /tmp/qiime2-q2cli-err-utwt1cmu.log
    
## Citation引文

    If used this script, please cited:
    使用此脚本，请引用下文：
    
    Bolyen et al. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. 
    **Nature Biotechnology** 37: 852-857. https://doi.org/10.1038/s41587-019-0209-9

    Yong-Xin Liu, Lei Chen, Tengfei Ma, et al. 2023. 
    EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. 
    iMeta 2: e83. https://doi.org/10.1002/imt2.83
    
    Yousuf, et al. 2024. Unveiling microbial communities with EasyAmplicon: 
    A user-centric guide to perform amplicon sequencing data analysis. 
    iMetaOmics 1: e42. https://doi.org/10.1002/imo2.42
    
    Copyright 2016-2025 Yong-Xin Liu <liuyongxin@caas.cn>
