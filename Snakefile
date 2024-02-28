configfile: "config.yaml"

PHENO=sorted(glob_wildcards(config["DEST"]+"/inputs/pheno/"+"{phenos}"+".txt").phenos)

PHENO_UNIQ=list(PHENO)

rule all:
    input:
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.FDR.txt",PHENO_UNIQ=PHENO_UNIQ),
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.wh.txt",PHENO_UNIQ=PHENO_UNIQ),
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.filtered.wh.txt",PHENO_UNIQ=PHENO_UNIQ),
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.man.png",PHENO_UNIQ=PHENO_UNIQ),
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.qq.png",PHENO_UNIQ=PHENO_UNIQ),
        expand(config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.filtered.annot.csv",PHENO_UNIQ=PHENO_UNIQ)

rule input_filter:
    input:
        INV=config["INPUTS"]+"/"+"DGRP-3_inv_5per_covar.txt",
        WOLB=config["INPUTS"]+"/"+"DGRP-3_wolbachia_covar.txt",
        PHENO=config["INPUTS"]+"/pheno/"+"{PHENO_UNIQ}"+".txt"
    output:
        PHENO_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"{PHENO_UNIQ}"+"_filtered.txt",
        COVAR_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_covar_filtered.txt",
        SIG_COVAR_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_sig_covar_filtered.txt",
        SIG_COVAR_NAME=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_sig_covars.txt",
        TEST_INDIV=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_association_test_indiv.txt"
    params:
        SCRIPT=config["SCRIPTS"]+"/"+"processing_modified.R",
        CUTOFF=config["CUTOFF"]
    resources: cpus=1, mem_mb=4000, time_min=1440
    shell:
        """
        ml R
        Rscript {params.SCRIPT} \
        {input.INV} \
        {input.WOLB} \
        {input.PHENO} \
        {output.PHENO_FILT} \
        {output.COVAR_FILT} \
        {output.SIG_COVAR_FILT} \
        {output.SIG_COVAR_NAME} \
        {output.TEST_INDIV}\
        {params.CUTOFF}
        """

rule vcf_filter:
    input:
        TEST_INDIV=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_association_test_indiv.txt",
        VCF=config["VCF_LOC"]+"/"+config["VCF"]+".vcf.gz"
    output:
        FILT_VCF=config["DEST"]+"/{PHENO_UNIQ}/"+config["VCF"]+"_gemma.vcf.gz"
    params:
        THREADS=config["THREADS"]
    resources: cpus=8, mem_mb=32000, time_min=1440
    shell:
        """
        ml bcftools
        bcftools view \
        -Oz \
        --threads {params.THREADS} \
        -S {input.TEST_INDIV} \
        {input.VCF} > \
        {output.FILT_VCF}
        """

rule convert_2_plink:
    input:
        FILT_VCF=config["DEST"]+"/{PHENO_UNIQ}/"+config["VCF"]+"_gemma.vcf.gz",
        PHENO_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"{PHENO_UNIQ}"+"_filtered.txt"
    output:
        FAM=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"]+".fam"
    params:
        GENO=config["GENO"],
        MAF=config["MAF"],
        CHR=config["CHR"],
        BFILE=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"]
    resources: cpus=1, mem_mb=16000, time_min=1440
    shell:
        """
        ml plink
        plink \
            --allow-extra-chr \
            --vcf {input.FILT_VCF} \
            --make-bed \
            --double-id \
            --geno {params.GENO} \
            --pheno {input.PHENO_FILT} \
            --maf {params.MAF} \
            --chr {params.CHR} \
            --out {params.BFILE}
        """

rule make_grm:
    input:
        FAM=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"]+".fam"
    output:
        BIN=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+config["BFILE"]+"_grm.cXX.txt"
    params:
        BFILE=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"],
        DEST=config["DEST"]+"/{PHENO_UNIQ}/",
        THREADS=config["THREADS"],
        GRM=config["BFILE"]+"_grm"
    resources: cpus=8, mem_mb=32000, time_min=1440
    shell:
        """
        ml gemma
        #module load openblas/0.3.7
        ml gsl/2.7

        OPENBLAS_NUM_THREADS={params.THREADS}
        cd {params.DEST}
            gemma-0.98.5-linux-static-AMD64 \
            -bfile {params.BFILE} \
            -gk 1 \
            -o {params.GRM} \
        """

rule association:
    input:
        FAM=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"]+".fam",
        KIN=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+config["BFILE"]+"_grm.cXX.txt",
        SIG_COVAR_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"DGRP-3_sig_covar_filtered.txt"
    output:
        config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.txt"
    params:
        LMM=config["LMM"],
        MAF=config["MAF"],
        THREADS=config["THREADS"],
        ONAME="{PHENO_UNIQ}"+".wald",
        BFILE=config["DEST"]+"/{PHENO_UNIQ}/"+config["BFILE"],
        DEST=config["DEST"]+"/{PHENO_UNIQ}/"
    resources: cpus=8, mem_mb=32000, time_min=1440
    shell:
        """
        ml gemma
        #module load openblas/0.3.7
        ml gsl/2.7

        OPENBLAS_NUM_THREADS={params.THREADS}
        cd {params.DEST}

        gemma-0.98.5-linux-static-AMD64 \
            -bfile {params.BFILE} \
            -k {input.KIN} \
            -lmm {params.LMM} \
            -maf {params.MAF} \
            -o {params.ONAME} \
            -c {input.SIG_COVAR_FILT}
        """

rule post_processing:
    input:
        config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.txt"
    output:
        FDR=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.FDR.txt",
        SORTED=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.wh.txt",
        SORTED_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.filtered.wh.txt",
        MAN=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.man.png",
        QQ=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.qq.png"
    params:
        PROC=config["SCRIPTS"]+"/"+"adjust_fdr.R",
        FIG=config["SCRIPTS"]+"/"+"make_figures.R"
    resources: cpus=1, mem_mb=32000, time_min=1440
    shell:
        """
        ml R/3.5.0
        Rscript {params.PROC} {input} {output.FDR} {output.SORTED} {output.SORTED_FILT}
        Rscript {params.FIG} {input} {output.QQ} {output.MAN}
        """

rule annotate:
    input:
        SORTED_FILT=config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.filtered.wh.txt",
        VEP=config["VEP"],
        GENE=config["GENE"]
    output:
        config["DEST"]+"/{PHENO_UNIQ}/"+"output/"+"{PHENO_UNIQ}"+".wald.assoc.sorted.filtered.annot.csv"
    params:
        config["SCRIPTS"]+"/"+"VEP_annotate.py"
    resources: cpus=1, mem_mb=32000, time_min=1440
    shell:
        """
        source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
        source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
        conda activate notebook_env
        python3 {params} {input.VEP} {input.SORTED_FILT} {input.GENE} {output}
        """