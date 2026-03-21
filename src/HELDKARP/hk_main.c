/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2003 by David Applegate, Robert Bixby,              */
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
/*                Demo of Held-Karp for Small Graphs                        */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 17, 1998                                                 */
/*        October 23, 2003 (bico, add code to grab tour)                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "heldkarp.h"

#define INIT_UPBOUND (-100000000)


int
    main (int ac, char **av);
static int
    parseargs (int ac, char **av),
    get_len_tlist (int ncount, int ecount, int *elist, int *elen,
        int *tour_elist, int *tour_elen);
static void
    usage (char *fname);


static int seed = 0;
static int nnodes_want = 0;
static char *in_fname = (char *) NULL;
static char *out_fname = (char *) NULL;
static int isdat = 0;
static int istsplib = 0;
static int binary_in = 0;
static int norm = CC_EUCLIDEAN;
static int nlim = -1;
static int upbound = INIT_UPBOUND;
static int anytour = 0;
static int silent = 0;

int main (int ac, char **av)
{
    int ncount, ecount, retval, havedat = 0;
    int foundtour = 0;
    int *elist = (int *) NULL;
    int *elen  = (int *) NULL;
    int *tour_elist = (int *) NULL;
    int *tour_elen = (int *) NULL;
    double value, upper;
    double *dupper = (double *) NULL;
    int rval = 0;
    CCdatagroup dat;
    double szeit;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    if (!nnodes_want && !in_fname) {
        usage (av[0]);
        rval = 1; goto CLEANUP;
    }

    CCutil_sprand (seed, &rstate);
    CC_PRINTF("Using random seed %d\n", seed); CC_FFLUSH(stdout);

    if (upbound != INIT_UPBOUND) {
        upper = (double) upbound;
        dupper = &upper;
        CC_PRINTF("Initial Upper Bound: %d\n", upbound); CC_FFLUSH(stdout);
    }

    if (isdat || istsplib || in_fname == (char *) NULL) {
        ncount = nnodes_want;
        if (istsplib) {
            rval = CCutil_gettsplib (in_fname, &ncount, &dat);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
            }
        } else {
            rval = CCutil_getdata (in_fname, binary_in, norm, &ncount, &dat,
                                   nnodes_want, 0, &rstate);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_getdata failed\n"); goto CLEANUP;
            }
        }
        havedat = 1;
        tour_elist = CC_SAFE_MALLOC (2*ncount, int);
        CCcheck_NULL (tour_elist, "out of memory for tour_elist");
        szeit = CCutil_zeit ();
        retval = CCheldkarp_small (ncount, &dat, dupper, &value, &foundtour,
                                   anytour, tour_elist, nlim, silent);
    } else {
        rval = CCutil_getedgelist_n (&ncount, in_fname, &ecount, &elist,
                                     &elen, 0);
        if (rval) {
            CC_FPRINTF(stderr, "CCutil_getedgelist failed\n"); goto CLEANUP;
        }
        tour_elist = CC_SAFE_MALLOC (2*ncount, int);
        CCcheck_NULL (tour_elist, "out of memory for tour_elist");
        szeit = CCutil_zeit ();
        retval = CCheldkarp_small_elist (ncount, ecount, elist, elen, dupper,
                                         &value, &foundtour, anytour,
                                         tour_elist,  nlim, silent);
    }

    if (retval == HELDKARP_SEARCHLIMITEXCEEDED) {
        CC_PRINTF("Search limit exceeded\n"); goto CLEANUP;
    } else if (retval) {
        CC_PRINTF("heldkarp error %d\n", retval); goto CLEANUP;
    } else if (foundtour == 0) { 
        CC_PRINTF("Found no better tour\n");
    } else if (anytour == 1) {
        CC_PRINTF("Found improved tour: %.0f\n", value);
    } else {
        CC_PRINTF("Optimized, value %.0f\n", value);
    }
    CC_PRINTF("Running Time: %.2f seconds\n", CCutil_zeit () - szeit); 
    CC_FFLUSH(stdout);

    if (foundtour == 1) {
        if (out_fname) {
            if (havedat) {
                rval = CCutil_writeedges (ncount, out_fname, ncount,
                           tour_elist, &dat, 0);
                CCcheck_rval (rval, "CCutil_writeedges failed");
            } else {
                tour_elen = CC_SAFE_MALLOC (ncount, int);
                CCcheck_NULL (tour_elen, "out of memory for tour_elen");
                rval = get_len_tlist (ncount, ecount, elist, elen, tour_elist,
                                      tour_elen);
                CCcheck_rval (rval, "get_len_tlist failed");
                rval = CCutil_writeedges_int (ncount, out_fname, ncount,
                           tour_elist, tour_elen, 0);
                CCcheck_rval (rval, "CCutil_writeedges_int failed");
            }

        }
    }

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (tour_elist, int);
    CC_IFFREE (tour_elen, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int get_len_tlist (int ncount, int ecount, int *elist, int *elen,
        int *tour_elist, int *tour_elen)
{
    int rval = 0;
    int end0, end1, i, j;

    for (i = 0; i < ncount; i++) {
        end0 = tour_elist[2*i];
        end1 = tour_elist[2*i+1];
        for (j = 0; j < ecount; j++) {
            if ((elist[2*j] == end0 && elist[2*j+1] == end1) ||
                (elist[2*j] == end1 && elist[2*j+1] == end0)) {
                tour_elen[i] = elen[j];
                break;
            }
        }
        if (j == ecount) {
            CC_PRINTF("Could not find edge (%d, %d)\n", end0, end1);
            CC_FFLUSH(stdout);
            rval = 1; goto CLEANUP;
        }
    }

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "abdk:n:N:o:qs:Tu:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'a':
            anytour = 1;
            break;
        case 'b':
            binary_in = 0;
            break;
        case 'd':
            isdat = 1;
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            break;
        case 'n':
            nlim = atoi (boptarg);
            break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_EUCLIDEAN_CEIL; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
            case 4: norm = CC_USER; break;
            case 5: norm = CC_ATT; break;
            case 6: norm = CC_GEOGRAPHIC; break;
            case 7: norm = CC_MATRIXNORM; break;
            case 8: norm = CC_DSJRANDNORM; break;
            case 9: norm = CC_CRYSTAL; break;
            case 11: norm = CC_RHMAP1; break;
            case 12: norm = CC_RHMAP2; break;
            case 13: norm = CC_RHMAP3; break;
            case 14: norm = CC_RHMAP4; break;
            case 15: norm = CC_RHMAP5; break;
            case 16: norm = CC_EUCTOROIDAL; break;
            case 17: norm = CC_GEOM; break;
            default:
                usage (av[0]);
                return 1;
            }
            isdat = 1;
            break;
        case 'o':
            out_fname = boptarg;
            break;
        case 'q':
            silent++;
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'T':
            istsplib = 1;
            break;
        case 'u':
            upbound = atoi (boptarg);
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }
    if (boptind < ac) {
        in_fname = av[boptind++];
    }
    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    CC_FPRINTF(stderr, "Usage: %s [-flags below] [filename]\n", fname);
    CC_FPRINTF(stderr, "   -a    cut off search after first tour\n");
    CC_FPRINTF(stderr, "   -b    dat file is in binary\n");
    CC_FPRINTF(stderr, "   -d    dat file (default is edge file)\n");
    CC_FPRINTF(stderr, "   -k #  number of nodes for random problem\n");
    CC_FPRINTF(stderr, "   -n #  limit on number of search nodes\n");
    CC_FPRINTF(stderr, "   -o f  output tour (as edgelist) to file\n");
    CC_FPRINTF(stderr, "   -q    less output (may be repeated)\n");
    CC_FPRINTF(stderr, "   -s #  random seed\n");
    CC_FPRINTF(stderr, "   -T    tsplib file (default is edge file)\n");
    CC_FPRINTF(stderr, "   -u #  upperbound on tour length\n");
    CC_FPRINTF(stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    CC_FPRINTF(stderr, "         0=MAX, 1=JOHNSON, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    CC_FPRINTF(stderr, "         8=DSJRAND, 9=CRYSTAL, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    CC_FPRINTF(stderr, "         17=GEOM\n");
}
