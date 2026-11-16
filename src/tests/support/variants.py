"""Mixed ploidy and spanning-deletion records shared by record/genotype queries."""

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=500>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc
chr1\t5\trs1\tA\tG\t50\tPASS\t.\tGT\t0|1\t1/1\t./.
chr1\t10\t.\tC\tCTT,G\t50\tPASS\t.\tGT\t2|0\t0|0\t1
chr1\t20\t.\tG\t<DEL>\t50\tPASS\tEND=120\tGT\t0/1\t.\t1|1
chr1\t300\t.\tT\t.\t50\tPASS\t.\tDP\t3\t4\t5
chr2\t7\t.\tA\tC\t50\tPASS\t.\tGT\t1\t0\t.
"""
