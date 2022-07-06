rule freebayes:
    input:
        bam= lambda wc: expand('{{genome}}/bwa/{library_id}.bam', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        bai= lambda wc: expand('{{genome}}/bwa/{library_id}.bam.bai', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        fasta= 'ref/{genome}.fasta',
        fai= 'ref/{genome}.fasta.fai',
    output:
        vcf= '{genome}/freebayes/variants.vcf.gz',
    conda:
        '../envs/freebayes.yaml',
    shell:
        r"""
        awk '{{print $1":0-"$2}}' {input.fai} \
        | freebayes-parallel - 8 \
            -f {input.fasta} --min-alternate-count 5 --ploidy 1 \
            --min-base-quality 5 \
            --min-alternate-fraction 0.1 --pooled-continuous {input.bam} \
        | awk -v OFS='\t' -v FS='\t' -v n=1 '{{if($0 !~ "^#" && $3 == "."){{$3=n; n++}} print $0}}' \
        | {workflow.basedir}/scripts/add_allele_frequency_vcf.py --format freebayes - \
        | bcftools norm -m- -f {input.fasta} \
        | bcftools view -O z > {output.vcf}
        tabix -f {output.vcf}
        """


rule freebayes_vep:
    input:
        vcf= '{genome}/freebayes/variants.vcf.gz',
        gff= 'ref/{genome}.gff.gz',
        fasta= 'ref/{genome}.fasta',
    output:
        vcf= '{genome}/freebayes/vep.vcf.gz',
    conda:
        '../envs/vep.yaml',
    shell:
        r"""
        vep -i {input.vcf} --vcf --gff {input.gff} --format vcf --compress_output bgzip \
            --fasta {input.fasta} --distance 1000 --force_overwrite -o {output.vcf}
        tabix -f {output.vcf}
        """
