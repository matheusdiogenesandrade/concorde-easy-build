/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                         MERGE CUTPOOLS                                   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 1997                                                       */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static int seed = 0;


static void
    usage (char *f);



int main (int ac, char **av)
{
    double szeit;
    int i, k, count;
    CCtsp_lpcut_in *c, *cnext;
    CCtsp_lpcut_in *cuts   = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcuts *pool     = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *nextpool = (CCtsp_lpcuts *) NULL;
    int ncount;
    int textout, textin;
    int *ptour = (int *) NULL;
    int *qtour = (int *) NULL;
    int *perm  = (int *) NULL;
    int rval = 0;
    CCrandstate rstate;

    if (ac < 2) {
        usage (*av);
        return 0;
    }

    szeit = CCutil_zeit ();
    seed = (int) CCutil_real_zeit ();
    CCutil_sprand (seed, &rstate);

    k = 1;
    if (av[k][0] == '-' && av[k][1] == 't') {
        CCdatagroup dat;
        CC_PRINTF("Read master ...\n");
        if (CCutil_getmaster (av[k+1], &ncount, &dat, &ptour)) {
            CC_FPRINTF(stderr, "CCutil_getmaster failed\n");
            return 1;
        }
        CCutil_freedatagroup (&dat);
        textout = 1;
        k += 2;
    } else {
        textout = 0;
    }

    if (av[k][0] == '-' && av[k][1] == 's') {
        CCdatagroup dat;
        CC_PRINTF("Read master ...\n");
        if (CCutil_getmaster (av[k+1], &ncount, &dat, &qtour)) {
            CC_FPRINTF(stderr, "CCutil_getmaster failed\n");
            return 1;
        }
        CCutil_freedatagroup (&dat);
        textin = 1;
        k += 2;
    } else {
        textin = 0;
    }

    if (av[k][0] == '-' && av[k][1] == 'p') {
        FILE *tin = fopen (av[k+1], "r");

        CC_PRINTF("Read permutation tour ...\n");
        if (!tin) {
            CC_FPRINTF(stderr, "could not open %s for reading\n", av[k+1]);
            rval = 1; goto CLEANUP;
        }

        if (fscanf (tin, "%d", &ncount) != 1) {
            CC_FPRINTF(stderr, "perm file in wrong format\n");
            rval = 1; fclose (tin); goto CLEANUP;
        }
        perm = CC_SAFE_MALLOC (ncount, int);
        if (!perm) {
            CC_FPRINTF(stderr, "out of memory in main\n");
            rval = 1; fclose (tin); goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            if (fscanf (tin, "%d", &(perm[i])) != 1) {
                CC_FPRINTF(stderr, "perm file in wrong format\n");
                rval = 1; fclose (tin); goto CLEANUP;
            }
        }
        fclose (tin);
        k += 2;
    }

    if (textin) {
        CC_PRINTF("Number of Nodes: %d\n", ncount); CC_FFLUSH(stdout);

        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
        if (rval) {
            CC_FPRINTF(stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
        }

        cuts = (CCtsp_lpcut_in *) NULL;
        rval = CCtsp_file_cuts (av[k], &cuts, &count, ncount, qtour);
        if (rval) {
            CC_FPRINTF(stderr, "CCtsp_file_cuts failed\n"); goto CLEANUP;
        }
        CC_PRINTF("File has %d cuts\n", count); CC_FFLUSH(stdout);
        for (c = cuts; c; c = cnext) {
            cnext = c->next;
            rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
            if (rval) {
                CC_FPRINTF(stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
                goto CLEANUP;
            }
            CCtsp_free_lpcut_in (c);
            CC_FREE (c, CCtsp_lpcut_in);
        }
        k++;
    } else {
        ncount = 0;
        rval = CCtsp_init_cutpool (&ncount, av[k], &pool);
        if (rval) {
            CC_FPRINTF(stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
        }
        CC_PRINTF("Initial Pool: %d nodes %d cuts\n", ncount, pool->cutcount);
        CC_FFLUSH(stdout);
        k++;
    }

    for (; k < ac; k++) {
        CC_PRINTF("Adding Pool %s ... ", av[k]);
        CC_FFLUSH(stdout);

        if (textin) {
            cuts = (CCtsp_lpcut_in *) NULL;
            rval = CCtsp_file_cuts (av[k], &cuts, &count, ncount, qtour);
            if (rval) {
                CC_FPRINTF(stderr, "CCtsp_file_cuts failed\n"); goto CLEANUP;
            }
            for (c = cuts; c; c = cnext) {
                cnext = c->next;
                rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
                if (rval) {
                    CC_FPRINTF(stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
                    goto CLEANUP;
                }
                CCtsp_free_lpcut_in (c);
                CC_FREE (c, CCtsp_lpcut_in);
            }
        } else {
            rval = CCtsp_init_cutpool (&ncount, av[k], &nextpool);
            if (rval) {
                CC_FPRINTF(stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
            }

            for (i = 0; i < nextpool->cutcount; i++) {
                rval = CCtsp_add_to_cutpool (pool, nextpool,
                                             &(nextpool->cuts[i]));
                if (rval) {
                    CC_FPRINTF(stderr, "CCtsp_add_to_cutpool failed\n");
                    goto CLEANUP;
                }
            }
            CCtsp_free_cutpool (&nextpool);
        }
        CC_PRINTF("%d\n", pool->cutcount); CC_FFLUSH(stdout);
    }


    CC_PRINTF("Final Pool: %d cuts\n", pool->cutcount);
    CC_FFLUSH(stdout);

    if (textout) {
        CC_PRINTF("Write text file ...\n"); CC_FFLUSH(stdout);
        if (perm) {
            int *pperm = (int *) NULL;

            pperm = CC_SAFE_MALLOC (ncount, int);
            if (!pperm) {
                CC_FPRINTF(stderr, "out of memory in main\n");
                rval = 1; goto CLEANUP;
            }
            for (i = 0; i < ncount; i++) {
                pperm[i] = perm[ptour[i]];
            }
            for (i = 0; i < ncount; i++) {
                ptour[i] = pperm[i];
            }
            CC_FREE (pperm, int);
        }
        rval = CCtsp_file_cuts_write ("merge.txt", pool, ptour);
        if (rval) {
            CC_FPRINTF(stderr, "CCtsp_file_cuts_write failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCtsp_write_cutpool (ncount, "merge.pul", pool);
        if (rval) {
            CC_FPRINTF(stderr, "CCtsp_write_cutpool failed\n");
            goto CLEANUP;
        }
    }

    CC_PRINTF("Running Time: %.2f seconds\n", CCutil_zeit () - szeit);
    CC_FFLUSH(stdout);

CLEANUP:

    if (pool) {
        CCtsp_free_cutpool (&pool);
    }
    if (nextpool) {
        CCtsp_free_cutpool (&nextpool);
    }
    CC_IFFREE (ptour, int);
    CC_IFFREE (qtour, int);

    return rval;
}

static void usage (char *f)
{
    CC_FPRINTF(stderr, "Usage: %s [-t] pool1 pool2 ... \n", f);
    CC_FPRINTF(stderr, "       -t f:  to write a text file specify a master\n");
    CC_FPRINTF(stderr, "       -p f:  to permute nodes when writing text\n");
    CC_FPRINTF(stderr, "       -s f:  to read text files specify a master\n");
    CC_FPRINTF(stderr, "Note: merged pool will be written merge.pul (.txt)\n");
    CC_FPRINTF(stderr, "      the master files are used to map the nodes\n");
    CC_FPRINTF(stderr, "      the permuation given by p can be used to\n");
    CC_FPRINTF(stderr, "      handle the different node orders in dat and tsp\n");
}
