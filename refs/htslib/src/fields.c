/* Field accessors for bcf1_t / bcf_hdr_t, whose layouts (bit-fields, nested
   structs) and helper macros are not stable across htslib versions. */

#include "htslib/vcf.h"

int64_t ow_bcf_pos(const bcf1_t *r) { return r->pos; }
int64_t ow_bcf_rlen(const bcf1_t *r) { return r->rlen; }
int ow_bcf_n_allele(const bcf1_t *r) { return r->n_allele; }
const char *ow_bcf_allele(const bcf1_t *r, int i) { return r->d.allele[i]; }
const char *ow_bcf_seqname(const bcf_hdr_t *h, const bcf1_t *r) { return bcf_seqname(h, r); }
int ow_bcf_hdr_nsamples(const bcf_hdr_t *h) { return bcf_hdr_nsamples(h); }

const char *ow_bcf_id(const bcf1_t *r) { return r->d.id; }
float ow_bcf_qual(const bcf1_t *r) { return r->qual; }
int ow_bcf_n_flt(const bcf1_t *r) { return r->d.n_flt; }
const int *ow_bcf_flt(const bcf1_t *r) { return r->d.flt; }
int ow_bcf_n_sample(const bcf1_t *r) { return r->n_sample; }
int ow_bcf_n_info(const bcf1_t *r) { return r->n_info; }
int ow_bcf_n_fmt(const bcf1_t *r) { return r->n_fmt; }

/* INFO field i after unpacking: key, BCF type, value count, value bytes. */
void ow_bcf_info(const bcf1_t *r, int i, int *key, int *type, int *len,
                 const uint8_t **p, int *nbytes)
{
    const bcf_info_t *x = &r->d.info[i];
    *key = x->key; *type = x->type; *len = x->len;
    *p = x->vptr; *nbytes = x->vptr_len;
}

/* FORMAT field i after unpacking: key, BCF type, values per sample, bytes
   of all samples. */
void ow_bcf_fmt(const bcf1_t *r, int i, int *id, int *type, int *n,
                const uint8_t **p, int *nbytes)
{
    const bcf_fmt_t *f = &r->d.fmt[i];
    *id = f->id; *type = f->type; *n = f->n;
    *p = f->p; *nbytes = f->p_len;
}

/* The dictionary of FILTER/INFO/FORMAT keys: its size and the name of id i
   (NULL for an unused id). */
int ow_bcf_hdr_nkeys(const bcf_hdr_t *h) { return h->n[BCF_DT_ID]; }
const char *ow_bcf_hdr_key(const bcf_hdr_t *h, int i) { return h->id[BCF_DT_ID][i].key; }
