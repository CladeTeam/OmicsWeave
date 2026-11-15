/* BAM record fields, index loading and region iteration for tests: the
   inline helpers and macros of htslib/sam.h as callable functions. */

#include "htslib/sam.h"

int ow_bam_tid(const bam1_t *b) { return b->core.tid; }
int64_t ow_bam_pos(const bam1_t *b) { return b->core.pos; }
int ow_bam_flag(const bam1_t *b) { return b->core.flag; }
int ow_bam_mapq(const bam1_t *b) { return b->core.qual; }
int ow_bam_mtid(const bam1_t *b) { return b->core.mtid; }
int64_t ow_bam_mpos(const bam1_t *b) { return b->core.mpos; }
int64_t ow_bam_isize(const bam1_t *b) { return b->core.isize; }
int ow_bam_l_qseq(const bam1_t *b) { return b->core.l_qseq; }
int ow_bam_n_cigar(const bam1_t *b) { return b->core.n_cigar; }
const char *ow_bam_qname(const bam1_t *b) { return bam_get_qname(b); }
const uint32_t *ow_bam_cigar(const bam1_t *b) { return bam_get_cigar(b); }
const uint8_t *ow_bam_seq(const bam1_t *b) { return bam_get_seq(b); }
const uint8_t *ow_bam_qual(const bam1_t *b) { return bam_get_qual(b); }
const uint8_t *ow_bam_aux(const bam1_t *b) { return bam_get_aux(b); }
int ow_bam_l_aux(const bam1_t *b) { return bam_get_l_aux(b); }

hts_idx_t *ow_sam_index_load(htsFile *fp, const char *fn, const char *fnidx)
{
    return sam_index_load2(fp, fn, fnidx);
}

hts_itr_t *ow_sam_itr_queryi(const hts_idx_t *idx, int tid, int64_t beg, int64_t end)
{
    return sam_itr_queryi(idx, tid, beg, end);
}

int ow_sam_itr_next(htsFile *fp, hts_itr_t *itr, bam1_t *b)
{
    return sam_itr_next(fp, itr, b);
}

void ow_hts_itr_destroy(hts_itr_t *itr) { hts_itr_destroy(itr); }
void ow_hts_idx_destroy(hts_idx_t *idx) { hts_idx_destroy(idx); }
