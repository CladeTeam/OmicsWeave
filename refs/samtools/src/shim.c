/* In-process entry to samtools commands for tests. Each command is called as
   samtools' own main() calls it, with argv[0] the command name and argv[-1]
   the program name. */

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include "samtools.h"

int bam_sort(int argc, char *argv[]);
int main_depth(int argc, char *argv[]);

/* From samtools bamtk.c, which also holds main(). */
const char *samtools_version(void)
{
    return OW_SAMTOOLS_VERSION;
}

/* bam2depth.c reads `-f FILE` lists with read_file_list from bam_plcmd.c
   (mpileup), which is not vendored; the option reports an error. */
int read_file_list(const char *file_list, int *n, char **argv[])
{
    (void)file_list;
    *n = 0;
    *argv = NULL;
    return 1;
}

typedef int (*command_fn)(int argc, char *argv[]);

static const struct {
    const char *name;
    command_fn run;
} commands[] = {
    { "sort", bam_sort },
    { "depth", main_depth },
};

/* Runs `samtools <argv[0]> argv[1..]`. Returns the command's exit status, or
   -1 for a command not built in. Not reentrant: getopt and the commands'
   own state are process-global. */
int ow_samtools(int argc, char **argv)
{
    if (argc < 1)
        return -1;
    for (size_t i = 0; i < sizeof(commands) / sizeof(commands[0]); i++) {
        if (strcmp(argv[0], commands[i].name) != 0)
            continue;
        /* The commands read argv[-1] (for @PG lines), as samtools' main
           passes them its own argv shifted by one. */
        char **full = malloc((argc + 2) * sizeof(char *));
        if (!full)
            return -1;
        full[0] = "samtools";
        memcpy(full + 1, argv, argc * sizeof(char *));
        full[argc + 1] = NULL;
        optind = 0; /* glibc: re-initialise getopt completely */
        int ret = commands[i].run(argc, full + 1);
        free(full);
        return ret;
    }
    return -1;
}
