//! Independent source-dialect annotation snippets.

pub const GENCODE_GTF: &str = "\
##description: evidence-based annotation of the human genome (GRCh38), version 44 (Ensembl 110)
##provider: GENCODE
##format: gtf
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825.1\"; gene_type \"lncRNA\"; gene_name \"DDX11L2\"; level 2; tag \"overlaps_pseudogene\";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825.1\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L2\"; transcript_type \"lncRNA\"; transcript_name \"DDX11L2-202\"; level 2; transcript_support_level \"1\"; tag \"basic\"; tag \"Ensembl_canonical\"; havana_transcript \"OTTHUMT00000362751.1\";
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id \"ENSG00000290825.1\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L2\"; transcript_type \"lncRNA\"; transcript_name \"DDX11L2-202\"; exon_number 1; exon_id \"ENSE00002234944.1\"; level 2; tag \"basic\"; tag \"Ensembl_canonical\";
chr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tgene_id \"ENSG00000290825.1\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L2\"; transcript_type \"lncRNA\"; transcript_name \"DDX11L2-202\"; exon_number 2; exon_id \"ENSE00003582793.1\"; level 2; tag \"basic\"; tag \"Ensembl_canonical\";
chr1\tHAVANA\texon\t13221\t14409\t.\t+\t.\tgene_id \"ENSG00000290825.1\"; transcript_id \"ENST00000456328.2\"; gene_type \"lncRNA\"; gene_name \"DDX11L2\"; transcript_type \"lncRNA\"; transcript_name \"DDX11L2-202\"; exon_number 3; exon_id \"ENSE00002312635.1\"; level 2; tag \"basic\"; tag \"Ensembl_canonical\";
chr1\tHAVANA\tgene\t14404\t29570\t.\t-\t.\tgene_id \"ENSG00000227232.5\"; gene_type \"unprocessed_pseudogene\"; gene_name \"WASH7P\"; level 2; hgnc_id \"HGNC:38034\";
chr1\tHAVANA\ttranscript\t14404\t29570\t.\t-\t.\tgene_id \"ENSG00000227232.5\"; transcript_id \"ENST00000488147.1\"; gene_type \"unprocessed_pseudogene\"; gene_name \"WASH7P\"; transcript_type \"unprocessed_pseudogene\"; transcript_name \"WASH7P-201\"; level 2; tag \"basic\"; tag \"Ensembl_canonical\";
chr1\tHAVANA\texon\t29534\t29570\t.\t-\t.\tgene_id \"ENSG00000227232.5\"; transcript_id \"ENST00000488147.1\"; gene_type \"unprocessed_pseudogene\"; gene_name \"WASH7P\"; transcript_type \"unprocessed_pseudogene\"; transcript_name \"WASH7P-201\"; exon_number 1; exon_id \"ENSE00001890219.1\"; level 2;
chr1\tHAVANA\texon\t24738\t24891\t.\t-\t.\tgene_id \"ENSG00000227232.5\"; transcript_id \"ENST00000488147.1\"; gene_type \"unprocessed_pseudogene\"; gene_name \"WASH7P\"; transcript_type \"unprocessed_pseudogene\"; transcript_name \"WASH7P-201\"; exon_number 2; exon_id \"ENSE00003507205.1\"; level 2;
chr1\tHAVANA\texon\t14404\t14501\t.\t-\t.\tgene_id \"ENSG00000227232.5\"; transcript_id \"ENST00000488147.1\"; gene_type \"unprocessed_pseudogene\"; gene_name \"WASH7P\"; transcript_type \"unprocessed_pseudogene\"; transcript_name \"WASH7P-201\"; exon_number 3; exon_id \"ENSE00001843071.1\"; level 2;
chr1\tHAVANA\tgene\t65419\t71585\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; gene_type \"protein_coding\"; gene_name \"OR4F5\"; level 2;
chr1\tHAVANA\ttranscript\t65419\t71585\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; gene_type \"protein_coding\"; gene_name \"OR4F5\"; transcript_type \"protein_coding\"; transcript_name \"OR4F5-201\"; level 2; protein_id \"ENSP00000493376.2\"; tag \"basic\"; tag \"Ensembl_canonical\"; tag \"MANE_Select\"; tag \"appris_principal_1\"; tag \"CCDS\";
chr1\tHAVANA\texon\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 1; level 2;
chr1\tHAVANA\texon\t65520\t65573\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 2; level 2;
chr1\tHAVANA\tCDS\t65565\t65573\t.\t+\t0\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 2; level 2;
chr1\tHAVANA\tstart_codon\t65565\t65567\t.\t+\t0\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 2; level 2;
chr1\tHAVANA\texon\t69037\t71585\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 3; level 2;
chr1\tHAVANA\tCDS\t69037\t70005\t.\t+\t0\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 3; level 2;
chr1\tHAVANA\tstop_codon\t70006\t70008\t.\t+\t0\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 3; level 2;
chr1\tHAVANA\tUTR\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 1; level 2;
chr1\tHAVANA\tUTR\t65520\t65564\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 2; level 2;
chr1\tHAVANA\tUTR\t70006\t71585\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; transcript_id \"ENST00000641515.2\"; exon_number 3; level 2;
chrM\tENSEMBL\tgene\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888.2\"; gene_type \"protein_coding\"; gene_name \"MT-ND1\"; level 3;
chrM\tENSEMBL\ttranscript\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888.2\"; transcript_id \"ENST00000361390.2\"; gene_type \"protein_coding\"; gene_name \"MT-ND1\"; transcript_type \"protein_coding\"; transcript_name \"MT-ND1-201\"; level 3; tag \"basic\";
chrM\tENSEMBL\texon\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888.2\"; transcript_id \"ENST00000361390.2\"; exon_number 1; level 3;
chrM\tENSEMBL\tCDS\t3307\t4262\t.\t+\t0\tgene_id \"ENSG00000198888.2\"; transcript_id \"ENST00000361390.2\"; exon_number 1; level 3;
";

pub const ENSEMBL_GTF: &str = "\
#!genome-build GRCh38.p14
#!genome-version GRCh38
#!genebuild-last-updated 2023-03
1\tensembl_havana\tgene\t65419\t71585\t.\t+\t.\tgene_id \"ENSG00000186092\"; gene_version \"7\"; gene_name \"OR4F5\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\";
1\tensembl_havana\ttranscript\t65419\t71585\t.\t+\t.\tgene_id \"ENSG00000186092\"; gene_version \"7\"; transcript_id \"ENST00000641515\"; transcript_version \"2\"; gene_name \"OR4F5\"; gene_biotype \"protein_coding\"; transcript_name \"OR4F5-201\"; transcript_biotype \"protein_coding\"; tag \"basic\"; tag \"Ensembl_canonical\"; tag \"MANE_Select\";
1\tensembl_havana\texon\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"1\"; exon_id \"ENSE00003812156\";
1\tensembl_havana\texon\t65520\t65573\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"2\"; exon_id \"ENSE00003813641\";
1\tensembl_havana\tCDS\t65565\t65573\t.\t+\t0\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"2\"; protein_id \"ENSP00000493376\";
1\tensembl_havana\tstart_codon\t65565\t65567\t.\t+\t0\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"2\";
1\tensembl_havana\texon\t69037\t71585\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"3\"; exon_id \"ENSE00003813949\";
1\tensembl_havana\tCDS\t69037\t70005\t.\t+\t0\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"3\"; protein_id \"ENSP00000493376\";
1\tensembl_havana\tstop_codon\t70006\t70008\t.\t+\t0\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\"; exon_number \"3\";
1\tensembl_havana\tfive_prime_utr\t65419\t65433\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\";
1\tensembl_havana\tfive_prime_utr\t65520\t65564\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\";
1\tensembl_havana\tthree_prime_utr\t70009\t71585\t.\t+\t.\tgene_id \"ENSG00000186092\"; transcript_id \"ENST00000641515\";
MT\tinsdc\tgene\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888\"; gene_version \"2\"; gene_name \"MT-ND1\"; gene_source \"insdc\"; gene_biotype \"protein_coding\";
MT\tinsdc\ttranscript\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888\"; transcript_id \"ENST00000361390\"; gene_name \"MT-ND1\"; transcript_biotype \"protein_coding\"; tag \"basic\";
MT\tinsdc\texon\t3307\t4262\t.\t+\t.\tgene_id \"ENSG00000198888\"; transcript_id \"ENST00000361390\"; exon_number \"1\";
";

pub const REFSEQ_GTF: &str = "\
#gtf-version 2.2
#!genome-build GRCh38.p14
NC_000001.11\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"\"; db_xref \"GeneID:100287102\"; db_xref \"HGNC:HGNC:37102\"; gbkey \"Gene\"; gene \"DDX11L1\"; gene_biotype \"transcribed_pseudogene\"; pseudo \"true\";
NC_000001.11\tBestRefSeq\ttranscript\t11874\t14409\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; db_xref \"GeneID:100287102\"; gbkey \"misc_RNA\"; gene \"DDX11L1\"; product \"DEAD/H-box helicase 11 like 1 (pseudogene)\"; pseudo \"true\"; transcript_biotype \"transcript\";
NC_000001.11\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; gene \"DDX11L1\"; exon_number \"1\";
NC_000001.11\tBestRefSeq\texon\t12613\t12721\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; gene \"DDX11L1\"; exon_number \"2\";
NC_000001.11\tBestRefSeq\texon\t13221\t14409\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018.2\"; gene \"DDX11L1\"; exon_number \"3\";
NC_000001.11\tBestRefSeq\tgene\t65419\t71585\t.\t+\t.\tgene_id \"OR4F5\"; transcript_id \"\"; gbkey \"Gene\"; gene \"OR4F5\"; gene_biotype \"protein_coding\";
NC_000001.11\tBestRefSeq\ttranscript\t65419\t71585\t.\t+\t.\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; gbkey \"mRNA\"; gene \"OR4F5\"; tag \"MANE Select\"; tag \"RefSeq Select\"; transcript_biotype \"mRNA\";
NC_000001.11\tBestRefSeq\texon\t65419\t65433\t.\t+\t.\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; exon_number \"1\";
NC_000001.11\tBestRefSeq\texon\t65520\t65573\t.\t+\t.\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; exon_number \"2\";
NC_000001.11\tBestRefSeq\texon\t69037\t71585\t.\t+\t.\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; exon_number \"3\";
NC_000001.11\tBestRefSeq\tCDS\t65565\t65573\t.\t+\t0\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; protein_id \"NP_001005484.2\"; exon_number \"2\";
NC_000001.11\tBestRefSeq\tCDS\t69037\t70005\t.\t+\t0\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; protein_id \"NP_001005484.2\"; exon_number \"3\";
NC_000001.11\tBestRefSeq\tstart_codon\t65565\t65567\t.\t+\t0\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; exon_number \"2\";
NC_000001.11\tBestRefSeq\tstop_codon\t70006\t70008\t.\t+\t0\tgene_id \"OR4F5\"; transcript_id \"NM_001005484.2\"; exon_number \"3\";
";

/// UCSC table-browser export: no gene or transcript records, scores as
/// `0.000000`, one transcript id placed on two chromosomes.
pub const UCSC_GTF: &str = "\
chr1\thg38_knownGene\texon\t11869\t12227\t0.000000\t+\t.\tgene_id \"ENST00000456328.2\"; transcript_id \"ENST00000456328.2\"; gene_name \"DDX11L2\";
chr1\thg38_knownGene\texon\t12613\t12721\t0.000000\t+\t.\tgene_id \"ENST00000456328.2\"; transcript_id \"ENST00000456328.2\"; gene_name \"DDX11L2\";
chr1\thg38_knownGene\texon\t13221\t14409\t0.000000\t+\t.\tgene_id \"ENST00000456328.2\"; transcript_id \"ENST00000456328.2\"; gene_name \"DDX11L2\";
chrX\thg38_knownGene\tstart_codon\t1000\t1002\t0.000000\t-\t0\tgene_id \"PARG\"; transcript_id \"PART\";
chrX\thg38_knownGene\tCDS\t503\t1002\t0.000000\t-\t2\tgene_id \"PARG\"; transcript_id \"PART\";
chrX\thg38_knownGene\texon\t401\t1100\t0.000000\t-\t.\tgene_id \"PARG\"; transcript_id \"PART\";
chrX\thg38_knownGene\tstop_codon\t500\t502\t0.000000\t-\t0\tgene_id \"PARG\"; transcript_id \"PART\";
chrY\thg38_knownGene\texon\t401\t1100\t0.000000\t-\t.\tgene_id \"PARG\"; transcript_id \"PART\";
";

pub const ENSEMBL_GFF3: &str = "\
##gff-version 3
##sequence-region   1 1 248956422
1\tEnsembl\tchromosome\t1\t248956422\t.\t.\t.\tID=chromosome:1;Alias=CM000663.2,chr1,NC_000001.11
###
1\thavana\tncRNA_gene\t11869\t14409\t.\t+\t.\tID=gene:ENSG00000290825;Name=DDX11L2;biotype=lncRNA;gene_id=ENSG00000290825;version=1
1\thavana\tlnc_RNA\t11869\t14409\t.\t+\t.\tID=transcript:ENST00000456328;Parent=gene:ENSG00000290825;Name=DDX11L2-202;biotype=lncRNA;tag=basic,Ensembl_canonical;transcript_id=ENST00000456328;version=2
1\thavana\texon\t11869\t12227\t.\t+\t.\tParent=transcript:ENST00000456328;Name=ENSE00002234944;exon_id=ENSE00002234944;rank=1
1\thavana\texon\t12613\t12721\t.\t+\t.\tParent=transcript:ENST00000456328;Name=ENSE00003582793;exon_id=ENSE00003582793;rank=2
1\thavana\texon\t13221\t14409\t.\t+\t.\tParent=transcript:ENST00000456328;Name=ENSE00002312635;exon_id=ENSE00002312635;rank=3
###
1\t.\tbiological_region\t10469\t11240\t1.3e+03\t.\t.\texternal_name=oe %3D 0.79;logic_name=cpg
1\tensembl_havana\tgene\t65419\t71585\t.\t+\t.\tID=gene:ENSG00000186092;Name=OR4F5;biotype=protein_coding;gene_id=ENSG00000186092
1\tensembl_havana\tmRNA\t65419\t71585\t.\t+\t.\tID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515
1\tensembl_havana\tfive_prime_UTR\t65419\t65433\t.\t+\t.\tParent=transcript:ENST00000641515
1\tensembl_havana\texon\t65419\t65433\t.\t+\t.\tParent=transcript:ENST00000641515;exon_id=ENSE00003812156;rank=1
1\tensembl_havana\texon\t65520\t65573\t.\t+\t.\tParent=transcript:ENST00000641515;exon_id=ENSE00003813641;rank=2
1\tensembl_havana\tCDS\t65565\t65573\t.\t+\t0\tID=CDS:ENSP00000493376;Parent=transcript:ENST00000641515;protein_id=ENSP00000493376
1\tensembl_havana\texon\t69037\t71585\t.\t+\t.\tParent=transcript:ENST00000641515;exon_id=ENSE00003813949;rank=3
1\tensembl_havana\tCDS\t69037\t70008\t.\t+\t0\tID=CDS:ENSP00000493376;Parent=transcript:ENST00000641515;protein_id=ENSP00000493376
1\tensembl_havana\tthree_prime_UTR\t70009\t71585\t.\t+\t.\tParent=transcript:ENST00000641515
";

pub const REFSEQ_GFF3: &str = "\
##gff-version 3
#!gff-spec-version 1.21
##sequence-region NC_000001.11 1 248956422
NC_000001.11\tRefSeq\tregion\t1\t248956422\t.\t+\t.\tID=NC_000001.11:1..248956422;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
NC_000001.11\tBestRefSeq\tpseudogene\t11874\t14409\t.\t+\t.\tID=gene-DDX11L1;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;gbkey=Gene;gene=DDX11L1;gene_biotype=transcribed_pseudogene;pseudo=true
NC_000001.11\tBestRefSeq\ttranscript\t11874\t14409\t.\t+\t.\tID=rna-NR_046018.2;Parent=gene-DDX11L1;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1 (pseudogene);pseudo=true;transcript_id=NR_046018.2
NC_000001.11\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tID=exon-NR_046018.2-1;Parent=rna-NR_046018.2;gbkey=misc_RNA;gene=DDX11L1
NC_000001.11\tBestRefSeq\texon\t12613\t12721\t.\t+\t.\tID=exon-NR_046018.2-2;Parent=rna-NR_046018.2;gbkey=misc_RNA;gene=DDX11L1
NC_000001.11\tBestRefSeq\texon\t13221\t14409\t.\t+\t.\tID=exon-NR_046018.2-3;Parent=rna-NR_046018.2;gbkey=misc_RNA;gene=DDX11L1
NC_000001.11\tBestRefSeq\tgene\t65419\t71585\t.\t+\t.\tID=gene-OR4F5;Name=OR4F5;gbkey=Gene;gene=OR4F5;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tmRNA\t65419\t71585\t.\t+\t.\tID=rna-NM_001005484.2;Parent=gene-OR4F5;gbkey=mRNA;gene=OR4F5;tag=MANE Select;transcript_id=NM_001005484.2
NC_000001.11\tBestRefSeq\texon\t65419\t65433\t.\t+\t.\tID=exon-NM_001005484.2-1;Parent=rna-NM_001005484.2
NC_000001.11\tBestRefSeq\texon\t65520\t65573\t.\t+\t.\tID=exon-NM_001005484.2-2;Parent=rna-NM_001005484.2
NC_000001.11\tBestRefSeq\texon\t69037\t71585\t.\t+\t.\tID=exon-NM_001005484.2-3;Parent=rna-NM_001005484.2
NC_000001.11\tBestRefSeq\tCDS\t65565\t65573\t.\t+\t0\tID=cds-NP_001005484.2;Parent=rna-NM_001005484.2;protein_id=NP_001005484.2
NC_000001.11\tBestRefSeq\tCDS\t69037\t70008\t.\t+\t0\tID=cds-NP_001005484.2;Parent=rna-NM_001005484.2;protein_id=NP_001005484.2
NC_000014.9\tBestRefSeq\tgene\t105586437\t105588395\t.\t-\t.\tID=gene-IGHV1-2;Name=IGHV1-2;gbkey=Gene;gene=IGHV1-2;gene_biotype=V_segment
NC_000014.9\tBestRefSeq\tV_gene_segment\t105586437\t105588395\t.\t-\t.\tID=id-IGHV1-2;Parent=gene-IGHV1-2;gbkey=V_segment;gene=IGHV1-2
NC_000014.9\tBestRefSeq\texon\t105588046\t105588395\t.\t-\t.\tID=exon-id-IGHV1-2-1;Parent=id-IGHV1-2
NC_000014.9\tBestRefSeq\texon\t105586437\t105586937\t.\t-\t.\tID=exon-id-IGHV1-2-2;Parent=id-IGHV1-2
NC_012920.1\tRefSeq\tgene\t3307\t4262\t.\t+\t.\tID=gene-ND1;Name=ND1;gbkey=Gene;gene=ND1;gene_biotype=protein_coding
NC_012920.1\tRefSeq\tCDS\t3307\t4262\t.\t+\t0\tID=cds-YP_003024026.1;Parent=gene-ND1;gbkey=CDS;gene=ND1;protein_id=YP_003024026.1
##FASTA
>NC_012920.1
GATCACAGGTCTATCACCC
";
