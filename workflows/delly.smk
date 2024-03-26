rule delly_call:
    input:
        bam= lambda wc: expand('{{genome}}/bwa/{library_id}.bam', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        bai= lambda wc: expand('{{genome}}/bwa/{library_id}.bam.bai', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        fasta= 'ref/{genome}.fasta',
    output:
        bcf= '{genome}/delly/tmp/sv.bcf',
    shell:
        r"""
        delly call --map-qual 10 --min-clique-size 5 -g {input.fasta} -o {output.bcf} {input.bam}
        """


rule delly_merge:
    input:
        bam= lambda wc: expand('{{genome}}/bwa/{library_id}.bam', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        bai= lambda wc: expand('{{genome}}/bwa/{library_id}.bam.bai', library_id= ss[ss.genome == wc.genome].library_id.unique()),
        bcf= '{genome}/delly/tmp/sv.bcf',
        fasta= 'ref/{genome}.fasta',
    output:
        vcf= '{genome}/delly/variants.vcf.gz',
    shell:
        r"""
        delly merge -o {output.vcf}.merge.bcf {input.bcf}
        delly call --map-qual 10 --min-clique-size 5 -g {input.fasta} -v {output.vcf}.merge.bcf -o {output.vcf}.gentoype.bcf {input.bam}

        bcftools filter -i '(FMT/DR + FMT/DV + FMT/RR + FMT/RV) > 0' {output.vcf}.gentoype.bcf \
        | {workflow.basedir}/scripts/add_allele_frequency_vcf.py --format delly \
        | bcftools view -O z > {output.vcf}
        tabix -f {output.vcf}

        rm {output.vcf}.merge.bcf
        rm {output.vcf}.merge.bcf.csi
        rm {output.vcf}.gentoype.bcf
        rm {output.vcf}.gentoype.bcf.csi
        """


rule delly_vep:
    input:
        vcf= '{genome}/delly/variants.vcf.gz',
        gff= 'ref/{genome}.gff.gz',
        fasta= 'ref/{genome}.fasta',
    output:
        vcf= '{genome}/delly/vep.vcf.gz',
    conda:
        '../envs/vep.yaml',
    shell:
        r"""
        vep -i {input.vcf} --vcf --gff {input.gff} --format vcf --compress_output bgzip --species NA \
            --fasta {input.fasta} --distance 1000 --force_overwrite -o {output.vcf}
        tabix -f {output.vcf}
        """
