/* In-process entries to `bcftools consensus` and to the -i/-e filter engine
   (filter.c) for tests. bcftools' error() exits the process; here it returns
   to the caller of the ow_* entry point instead. */

#include <getopt.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "bcftools.h"
#include "filter.h"

int main_consensus(int argc, char *argv[]);

/* From bcftools version.c, which also defines the exiting error(). */
int apply_verbosity(const char *str)
{
    char *tmp;
    int verbose = strtol(str,&tmp,10);
    if ( *tmp || verbose<0 ) return -1;
    if ( verbose > 3 ) hts_verbose = verbose;
    return 0;
}

int parse_overlap_option(const char *arg)
{
    if ( strcasecmp(arg, "pos") == 0 || strcmp(arg, "0") == 0 ) return 0;
    else if ( strcasecmp(arg, "record") == 0 || strcmp(arg, "1") == 0 ) return 1;
    else if ( strcasecmp(arg, "variant") == 0 || strcmp(arg, "2") == 0 ) return 2;
    else return -1;
}

static jmp_buf on_error;
static int in_call;

static void report(const char *format, va_list ap) {
    vfprintf(stderr, format, ap);
}

void error(const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    report(format, ap);
    va_end(ap);
    if (in_call) longjmp(on_error, 1);
    exit(-1);
}

void error_errno(const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    report(format, ap);
    va_end(ap);
    if (in_call) longjmp(on_error, 1);
    exit(-1);
}

/* filter_init(hdr, str); NULL when bcftools reports an error. Not
   reentrant. */
filter_t *ow_filter_init(bcf_hdr_t *hdr, const char *str) {
    if (setjmp(on_error)) {
        in_call = 0;
        return NULL;
    }
    in_call = 1;
    filter_t *filter = filter_init(hdr, str);
    in_call = 0;
    return filter;
}

/* filter_test(filter, rec, NULL), the site-level test: 1 or 0, or -1 when
   bcftools reports an error (the filter's state is then that of an
   interrupted evaluation). Not reentrant. */
int ow_filter_test(filter_t *filter, bcf1_t *rec) {
    if (setjmp(on_error)) {
        in_call = 0;
        return -1;
    }
    in_call = 1;
    int ret = filter_test(filter, rec, NULL);
    in_call = 0;
    return ret;
}

/* Runs `bcftools consensus` with argv (argv[0] = "consensus"). Returns its
   status, or -1 when it reports an error. Not reentrant: getopt and the
   error handler are process-global. */
int ow_bcftools_consensus(int argc, char **argv) {
    if (setjmp(on_error)) {
        in_call = 0;
        return -1;
    }
    in_call = 1;
    optind = 0; /* glibc: re-initialise getopt completely */
    int ret = main_consensus(argc, argv);
    in_call = 0;
    return ret;
}
