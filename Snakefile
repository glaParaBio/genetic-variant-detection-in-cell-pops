VERSION = "1.1.1"

import pandas
import os
import gzip


def is_gzip(filename):
    with gzip.open(filename, "r") as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False


def add_fastq_pair_id(ss):
    fq = ss[pandas.isna(ss.fastq_r1) == False][
        ["library_id", "fastq_r1"]
    ].drop_duplicates()
    fq["basename"] = [os.path.basename(x) for x in fq.fastq_r1]
    fq["basename"] = [re.sub("\.gz$", "", x) for x in fq.basename]
    fq["basename"] = [
        re.sub("\.fastq$|\.fq$|\.bz2$|\.txt$", "", x) for x in fq.basename
    ]
    fq["idx"] = ["%.04d" % (i + 1) for i in range(len(fq))]
    fq["pair_id"] = fq["idx"] + "." + fq["basename"]
    mrg = pandas.merge(
        ss,
        fq[["library_id", "fastq_r1", "pair_id"]],
        how="left",
        on=["library_id", "fastq_r1"],
    )
    assert len(ss) == len(mrg)
    return mrg


genomes = pandas.read_csv(config["genomes"], sep="\t", comment="#").dropna(how="all")

ss = (
    pandas.read_csv(
        config["ss"],
        sep="\t",
        comment="#",
        usecols=["library_id", "genome", "cutadapt_adapter", "fastq_r1", "fastq_r2"],
    )
    .dropna(how="all")
    .drop_duplicates()
)
ss = add_fastq_pair_id(ss)

# if not (len(ss[['run_id', 'fastq_r1']].drop_duplicates()) == len(set(ss.fastq_r1)) == len(set(ss.run_id))):
#     raise Exception('run_id and fastq files must be in 1-to-1 relationship')
#
# if ss[['library_id', 'run_id']].drop_duplicates().groupby('run_id').count().max()[0] > 1:
#     raise Exception("Some run_id is assigned to multiple library_id's")

ss_fastq = ss[["pair_id", "cutadapt_adapter", "fastq_r1", "fastq_r2"]].drop_duplicates()

if len(ss_fastq.pair_id) != len(set(ss_fastq.pair_id)):
    raise Exception("Some pair_id's have more than one 'cutadapt_adapter' entry")


assert len(genomes.genome) == len(set(genomes.genome))
assert all([x in list(genomes.genome) for x in ss.genome])

ss = pandas.merge(ss, genomes, on="genome")

fastqc_reports = {}
for fq in list(set(ss.fastq_r1)) + list(set(ss.fastq_r2)):
    fq_name = re.sub("\.gz$|\.bz2$", "", os.path.basename(fq))
    if (
        fq_name.endswith(".fastq")
        or fq_name.endswith(".fq")
        or fq_name.endswith(".txt")
        or fq_name.endswith(".csfastq")
    ):
        fq_name = re.sub("\.fastq$|\.fq$|\.txt$|\.csfastq$", "", fq_name)
    else:
        raise Exception(
            'To be fixed: Fastq filenames should have extension "[.fastq|.fq|.txt|.csfastq][.gz|.bz2]".\nGot instead "%s"'
            % fq
        )
    assert fq_name not in fastqc_reports
    fastqc_reports[fq_name] = fq


wildcard_constraints:
    library_id="|".join([re.escape(x) for x in ss.library_id]),
    pair_id="|".join([re.escape(x) for x in ss.pair_id]),
    genome="|".join([re.escape(x) for x in genomes.genome]),


rule all:
    input:
        "multiqc/fastqc_report.html",
        expand("{genome}/mosdepth/dist.html", genome=ss.genome),
        expand(
            "{genome}/{caller}/variants.vcf.gz",
            caller=["delly", "freebayes"],
            genome=ss.genome.unique(),
        ),
        # Annotate with VEP only if you have a gff
        expand(
            "{genome}/{caller}/vep.tsv.gz",
            caller=["delly", "freebayes"],
            genome=ss[pandas.isna(ss.gff) == False].genome,
        ),


include: "workflows/delly.smk"
include: "workflows/freebayes.smk"


rule download_genome:
    output:
        fasta="ref/{genome}.fasta",
    params:
        ftp=lambda wc: genomes[genomes.genome == wc.genome].fasta.iloc[0],
    run:
        fn_tmp = output.fasta + ".tmp"

        if os.path.isfile(params.ftp):
            shell("cp {params.ftp} {fn_tmp}")
        else:
            shell("curl --silent -L {params.ftp} > {fn_tmp}")

        if is_gzip(fn_tmp):
            shell("gzip -cd {fn_tmp} > {output.fasta}")
        else:
            shell("mv {fn_tmp} {output.fasta}")


rule faidx_genome:
    input:
        fasta="ref/{genome}.fasta",
    output:
        fai="ref/{genome}.fasta.fai",
    shell:
        r"""
        samtools faidx {input.fasta}
        """


rule download_gff:
    output:
        gff=temp("ref/{genome}.gff"),
    params:
        ftp=lambda wc: genomes[genomes.genome == wc.genome].gff.iloc[0],
    run:
        fn_tmp = output.gff + ".tmp"

        if os.path.isfile(params.ftp):
            shell("cp {params.ftp} {fn_tmp}")
        else:
            shell("curl --silent -L {params.ftp} > {fn_tmp}")

        if is_gzip(fn_tmp):
            shell("gzip -cd {fn_tmp} > {output.gff}")
        else:
            shell("mv {fn_tmp} {output.gff}")


rule reformat_gff:
    # Make GFF compatible with VEP
    input:
        gff="ref/{genome}.gff",
    output:
        gff="ref/{genome}.gff.gz",
    shell:
        r"""
        awk -v FS='\t' -v OFS='\t' '$1 !~ "^#" {{
            biotype = ";biotype="$3;
            if($3 == "protein_coding_gene") {{$3 = "gene"}}
            if($3 == "ncRNA_gene") {{$3 = "gene"}}
            print $0 biotype}}' {input.gff} \
        | sort -k1,1 -k4,4n -k5,5n -t$'\t' \
        | bgzip > {output.gff}

        tabix --force --preset gff {output.gff}
        """


rule gene_description:
    input:
        gff="ref/{genome}.gff.gz",
    output:
        tsv="ref/{genome}.gene_description.tsv",
    shell:
        r"""
        zcat {input.gff} \
        | {workflow.basedir}/scripts/getGeneDescriptionFromGFF.py > {output.tsv}
        """


rule bwa_index:
    input:
        fasta="ref/{genome}.fasta",
    output:
        index="ref/{genome}.fasta.bwt",
    shell:
        r"""
        bwa index {input.fasta}
        """


rule cutadapt:
    input:
        fastq_r1=lambda wc: os.path.abspath(
            ss_fastq[ss_fastq.pair_id == wc.pair_id].fastq_r1.iloc[0]
        ),
        fastq_r2=lambda wc: os.path.abspath(
            ss_fastq[ss_fastq.pair_id == wc.pair_id].fastq_r2.iloc[0]
        ),
    output:
        fastq_r1=temp("cutadapt/{pair_id}_R1.fastq.gz"),
        fastq_r2=temp("cutadapt/{pair_id}_R2.fastq.gz"),
        report="cutadapt/{pair_id}.log",
    params:
        adapter=lambda wc: ss_fastq[
            ss_fastq.pair_id == wc.pair_id
        ].cutadapt_adapter.iloc[0],
    run:
        if pandas.isna(params.adapter):
            cmd = r"""
                    ln -s {input.fastq_r1} {output.fastq_r1}
                    ln -s {input.fastq_r2} {output.fastq_r2}
                    touch {output.report}
                    """
else:
    cmd = r"""
                    cutadapt --minimum-length 10 --cores 8 -a {params.adapter} \
                        -A {params.adapter} -o {output.fastq_r1} -p {output.fastq_r2} \
                        {input.fastq_r1} {input.fastq_r2} > {output.report}
                    """
        shell(cmd)


rule bwa_mem_align:
    input:
        fasta="ref/{genome}.fasta",
        index="ref/{genome}.fasta.bwt",
        fastq_r1="cutadapt/{pair_id}_R1.fastq.gz",
        fastq_r2="cutadapt/{pair_id}_R2.fastq.gz",
    output:
        bam=temp("{genome}/bwa/{pair_id}.bam"),
    params:
        rg=lambda wc: r"@RG\tID:{pair_id}\tSM:{library_id}".format(
            pair_id=wc.pair_id,
            library_id=ss[
                (ss.genome == wc.genome) & (ss.pair_id == wc.pair_id)
            ].library_id.iloc[0],
        ),
    shell:
        r"""
        bwa mem -t 8 -R '{params.rg}' {input.fasta} {input.fastq_r1} {input.fastq_r2} \
        | samtools sort -@ 4 > {output.bam}
        """


rule merge_bam:
    input:
        bam=lambda wc: expand(
            "{{genome}}/bwa/{pair_id}.bam",
            pair_id=ss[
                (ss.library_id == wc.library_id) & (ss.genome == wc.genome)
            ].pair_id.unique(),
        ),
    output:
        bam="{genome}/bwa/{library_id}.bam",
        md="{genome}/bwa/{library_id}.md.txt",
    shell:
        r"""
        samtools merge -@ 4 - {input.bam} \
        | samtools sort -@ 4 -n \
        | samtools fixmate -m -@ 4 - - \
        | samtools sort -@ 4 \
        | samtools markdup -@ 4 -f {output.md} - {output.bam}
        """


rule index_and_stat_bam:
    # In part this rules prevents the index to look older then the bam
    # https://github.com/snakemake/snakemake/issues/1378
    input:
        bam="{genome}/bwa/{library_id}.bam",
    output:
        bai="{genome}/bwa/{library_id}.bam.bai",
        stats="{genome}/bwa/{library_id}.stats",
    shell:
        r"""
        samtools index -@ 4 {input.bam}
        samtools stats -@ 4 {input.bam} > {output.stats}
        """


rule mosdepth:
    input:
        bam="{genome}/bwa/{library_id}.bam",
        bai="{genome}/bwa/{library_id}.bam.bai",
    output:
        txt="{genome}/mosdepth/{library_id}.mosdepth.global.dist.txt",
    shell:
        r"""
        outdir=$(dirname {output.txt})
        mosdepth -x -n -t 4 $outdir/{wildcards.library_id} {input.bam}
        """


rule mosdepth_plot:
    input:
        txt=lambda wc: expand(
            "{{genome}}/mosdepth/{library_id}.mosdepth.global.dist.txt",
            library_id=ss[ss.genome == wc.genome].library_id.unique(),
        ),
    output:
        html="{genome}/mosdepth/dist.html",
    shell:
        r"""
        python {workflow.basedir}/scripts/mosdepth_script_plot_dist.py --output {output.html} {input.txt}
        """


rule vep_annotation_to_tsv:
    # Extract relavant information from vcf annotated with vep and melt by
    # sample to long-format table.
    input:
        vcf="{genome}/{caller}/vep.vcf.gz",
    output:
        tsv="{genome}/{caller}/vep.tsv.gz",
    shell:
        r"""
        echo "chrom pos locus_id ref alt qual type ref_depth alt_depth alt_freq sum_alt_freq symbol gene feature biotype consequence impact amino_acids codons library_id" \
        | tr ' ' '\t' \
        | gzip > {output.tsv}

        samples=$(bcftools query -l {input.vcf})
        for sample in $samples
        do
            bcftools view --samples $sample {input.vcf} \
            | bcftools +split-vep -d -f '%CHROM %POS %ID %REF %ALT %QUAL %TYPE [%AD{{0}}] [%AD{{1}}] [%ALT_AF] [%SUM_ALT_AF] %SYMBOL %Gene %Feature %BIOTYPE %Consequence %IMPACT %Amino_acids %Codons\n' \
            | tr ' ' '\t' \
            | awk -v FS='\t' -v OFS='\t' -v sample=$sample '{{print $0, sample}}' \
            | gzip >> {output.tsv}
        done
        """


rule fastqc:
    input:
        fastq=lambda wc: fastqc_reports[wc.fq_name],
    output:
        qc="fastqc/{fq_name}_fastqc.zip",
    shell:
        r"""
        fastqc --noextract --outdir fastqc {input.fastq}
        """


rule multiqc:
    input:
        fastqc_reports=expand(
            "fastqc/{fq_name}_fastqc.zip", fq_name=fastqc_reports.keys()
        ),
    output:
        "multiqc/fastqc_report.html",
    shell:
        r"""
        multiqc --force --outdir multiqc --filename fastqc_report.html {input.fastqc_reports}
        """
