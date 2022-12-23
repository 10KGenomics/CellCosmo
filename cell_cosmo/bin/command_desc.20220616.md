#!/usr/bin bash


WD=/home/bioinfo/workspace/bioinfo_upload_dir
GD=${WD}/genomic2/xx
SN=20220303HXY6

export CELL_COSMO_TEST_ENV_READ_N_ROWS=50000

cell_cosmo rna sample \
-c C5U3L15U3C6U3L6C5T30 \
-o ./ \
-s ${SN}


# 拆分barcode有参,
cell_cosmo rna barcode \
--fq1 ${WD}/20220303HXY6_BKDL220011268-1a_1.clean.fq.gz \
--fq2 ${WD}/20220303HXY6_BKDL220011268-1a_2.clean.fq.gz \
-c ${WD}/chemistry.ini \
-p C5U3L15U3C6U3L6C5T30 \
--use-polyt-valid-reads \
--polyt-rate 0.666 \
--use-barcode-valid-reads \
--allow-barcode-diff-num 1 \
-o ./ \
-s ${SN}
# 无参拆分barcode命令,
# 当使用无参方式时,`--use-barcode-valid-reads`和`--use-link-valid-reads`无效
#       此时仍可选择是否指定polyT和qual的相关参数对reads进行校验
# cell_cosmo rna barcode \
# --fq1 ${WD}/20220303HXY6_BKDL220011268-1a_1.clean.fq.gz \
# --fq2 ${WD}/20220303HXY6_BKDL220011268-1a_2.clean.fq.gz \
# -p C5U3L15U3C6U3L6C5T30 \
# -o ./ \
# -s 20220303HXY6


cell_cosmo rna cutadapt \
--fq ${SN}_2.fq \
-o ./ \
-s ${SN}

cell_cosmo rna star \
--fq ${SN}_clean_2.fq \
--genomeDir ${GD} \
-o ./ \
-s ${SN}

cell_cosmo rna featureCounts \
--input ${SN}_Aligned.sortedByCoord.out.bam \
--genomeDir ${GD} \
-o ./ \
-s ${SN}
# output 20220303HXY6_Aligned.sortedByCoord.out.bam.featureCounts.bam
# output 20220303HXY6_name_sorted.bam

# NOTE:
# 1. 输入必须使用 name sorted 的bam
# 2. 可以选择校正 barcode,`--barcode-diff`,`--barcode-correct-limit`
# 3. 可以选择校正 umi,`--umi-diff`,`--umi-correct-limit`
# 4. 默认校正 barcode,umi
cell_cosmo rna count \
-b ${SN}_name_sorted.bam \
--genomeDir ${GD} \
-o ./ \
-s ${SN}


cell_cosmo rna analysis \
--matrix-file ${SN}_filtered_feature_bc_matrix \
--genomeDir ${GD} \
--thread 24 \
-o ./ \
-s ${SN}


# -------- 使用pipeline方式调用
# 1. 生成配置文件,运行该命令会在当前目录生成配置文件 rna_pipeline.cfg
cell_cosmo rna pipeline -g 
# 2. 修改配置文件中的参数
#   required: global.fq1,global.fq2,global.genomeDir,barcode.pattern,barcode.chemistry_config
# 3. 运行pipeline
cell_cosmo rna pipeline -c rna_pipeline.cfg