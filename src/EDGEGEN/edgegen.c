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
/*                   COMPUTING INITIAL EDGE SETS                            */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 9, 1995                                                       */
/*  Changes: 6.8.1996 (bico) fixed bug in f2match nearest                   */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCedgegen_read (char *egname, CCedgegengroup *plan)                 */
/*    READS an edgegen description from file egname.                        */
/*      -egname the name of a file                                          */
/*      -plan returns the description of the mix of edges (can be used      */
/*       in a call to CCedgegen_edges () to obtain the edgeset).            */
/*                                                                          */
/*  int CCedgegen_edges (CCedgegengroup *plan, int ncount,                  */
/*      CCdatagroup *dat, double *wcoord, int *ecount, int **elist,         */
/*      int silent, CCrandstate *rstate)                                    */
/*    RETURNS the set of edges described in plan.                           */
/*      -plan describes the mix of edges                                    */
/*      -ncount is the number of nodes                                      */
/*      -dat contains the info to generate edge lengths                     */
/*      -wcoord are nodeweights for Held-Karp style edge lengths, using     */
/*       len[i,j] + wcoord[i] + wcoord[j] (wcoord can be NULL)              */
/*      -ecount returns the number of edges                                 */
/*      -elist returns the edges in end1 end2 format                        */
/*      -silent will suppress print messages if set to a nonzero value.     */
/*                                                                          */
/*  void CCedgegen_init_edgegengroup (CCedgegengroup *plan)                 */
/*        SETS the fields in plan to 0 (since there are so many fields to   */
/*          to deal with)                                                   */
/*                                                                          */
/*    NOTES:                                                                */
/*       To use CCedgegen_edges, look at the defintion of CCedgegengroup    */
/*    in edgegen.h - you should be able to guess what the parmeters mean.   */
/*    Note that wcoord is only used by a limited number of the generating   */
/*    routines, for example nearest, but not linkern.                       */
/*       The functions CCedgegen_edges and CCedgegen_read will return       */
/*    nonzero values if they fail (for example, if they run out of memory.  */
/*       The description file passed to CCedgegen_read should contain a     */
/*    list of some of the following commands:                               */
/*            EDGEGEN RANDOM #                                              */
/*                    -find # random edges                                  */
/*            EDGEGEN NEAREST #                                             */
/*                    -find the nearest # edges                             */
/*            EDGEGEN QUADNEAREST #                                         */
/*                    -find the quadrant-nearest # edges                    */
/*            EDGEGEN FRAC_TWOMATCH_NEAREST # [PRICED] [BASIC]              */
/*                    -find the nearest # using the reduced costs of a      */
/*                     fractional 2-matching as the edgelengths. If either  */
/*                     of the optional arguments PRICED or BASIC is         */
/*                     specified then the 2-matching used will be either    */
/*                     priced against the complete edgeset or converted to  */
/*                     a basic optimal solution (or both).                  */
/*            EDGEGEN GREEDY_TOUR                                           */
/*                    -find a greedy tour                                   */
/*            EDGEGEN BORUVKA_TOUR                                          */
/*                    -find a Boruvka tour                                  */
/*            EDGEGEN QBORUVKA_TOUR                                         */
/*                    -find a quick Boruvka tour                            */
/*            EDGEGEN NN_TOUR #                                             */
/*                    -find # nearest-neighbor tours                        */
/*            EDGEGEN RANDOM_TOUR #                                         */
/*                    -find # random tours                                  */
/*            EDGEGEN TWOOPT_TOUR #                                         */
/*                    -find # 2-opt tours                                   */
/*            EDGEGEN TWOPT5_TOUR #                                         */
/*                    -find # 2.5-opt tours                                 */
/*            EDGEGEN THREEOPT_TOUR #                                       */
/*                    -find # 3-opt tours                                   */
/*            EDGEGEN LINKERN #1 #2 [QUADNEAREST #3] [NEAREST #4]           */
/*                              [GREEDY_START | NN_START | RANDOM_START     */
/*                               | BORUVKA_START | QBORUVKA_START]          */
/*                    -find #1 Iterated Lin-Kernighan tours using #2        */
/*                     kicks.                                               */
/*                     The good edgeset can be specified by the optional    */
/*                     arguments QUADNEAREST and NEAREST (the two can be    */
/*                     used together). The initial tours can be specfied    */
/*                     by using one of GREEDY_START, NN_START,              */
/*                     BORUVKA_START, QBORUVKA_START, or RANDOM_START.      */
/*            EDGEGEN NN_TWOMATCH #                                         */
/*                    -find # nearest-neighbor 2-matchings                  */
/*            EDGEGEN TREE                                                  */
/*                    -find a minimum weight spanning tree.                 */
/*            EDGEGEN FRAC_TWOMATCH [PRICED] [BASIC]                        */
/*                    -find a minmum weight 2-matching (priced against the  */
/*                     complete edgeset) (that is a basic optimal           */
/*                     solution)                                            */
/*            EDGEGEN DELAUNAY                                              */
/*                    -find the edges in the (Euclidean-norm) Delaunay      */
/*                     triangulation (using Steve Fortune's sweep2 code).   */
/*            M_LINKERN                                                     */
/*                    -find # lin ker matchings                             */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "kdtree.h"
#include "linkern.h"
#include "fmatch.h"
#include "delaunay.h"
#include "mlinkern.h"
#include "macrorus.h"

#include <string.h>
#include <stdlib.h>

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct tabledat {
    intptr **table;
    int tabletotal;
    CCptrworld *intptr_world;
} tabledat;


static int
    call_kdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,
        double *wcoord, int *built_a_tree, int silent, CCrandstate *rstate),
    call_random_edge (int ncount, CCdatagroup *dat, int ecount, tabledat *td,
        int silent, CCrandstate *rstate),
    call_nearest (int ncount, CCdatagroup *dat, double *wcoord, int nearnum,
                  CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate),
    work_nearest (CCkdtree *kt, int ncount, int nearnum, CCdatagroup *dat,
        double *wcoord, int *ecount, int **elist, int silent,
        CCrandstate *rstate),
    call_quadnearest (int ncount, CCdatagroup *dat, double *wcoord,
        int nearnum, CCkdtree *kt, tabledat *td, int silent,
        CCrandstate *rstate),
    work_quadnearest (CCkdtree *kt, int ncount, int nearnum, CCdatagroup *dat,
        double *wcoord, int *ecount, int **elist, int silent,
        CCrandstate *rstate),
    call_delaunay (int ncount, CCdatagroup *dat, tabledat *td, int silent),
    call_mlinkern (int ncount, CCdatagroup *dat, CCkdtree *kt, int iterations,
        tabledat *td, int silent, CCrandstate *rstate),
    call_random_tour (int ncount, CCdatagroup *dat, int number, tabledat *td,
        int silent, CCrandstate *rstate),
    call_nearest_tour (int ncount, CCdatagroup *dat, int number, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate),
    work_nearest_tour (CCkdtree *kt, int ncount, int start, CCdatagroup *dat,
        int *tour, double *val, int silent, CCrandstate *rstate),
    call_greedy_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate),
    call_boruvka_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate),
    call_qboruvka_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate),
    call_twoopt_tour (int ncount, CCdatagroup *dat, CCkdtree *kt, int number,
        int two_and_a_half, int use_3opt, tabledat *td, int silent,
        CCrandstate *rstate),
    call_linkern (int ncount, CCdatagroup *dat, CCkdtree *kt,
        CCedgegengroup *plan, tabledat *td, int silent, CCrandstate *rstate),
    call_nearest_twomatch (int ncount, CCdatagroup *dat, int number,
        CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate),
    call_f2match (int ncount, CCdatagroup *dat, CCkdtree *kt, int priceit,
                  int basic, tabledat *td, int silent, CCrandstate *rstate),
    call_f2match_nearest (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int number, int priceit, int basic, tabledat *td, int silent,
        CCrandstate *rstate),
    call_spanning_tree (int ncount, CCdatagroup *dat, double *wcoord,
        CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate),
    f2match_initial_edgeset (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int *ecount, int **elist, int **elen, tabledat *td, int silent,
        CCrandstate *rstate),
    put_tour_in_table (tabledat *td, int ncount, int *tour),
    put_list_in_table (tabledat *td, int ecount, int *elist),
    put_in_table (tabledat *td, int i, int j),
    collect_table_edges (tabledat *td, int ncount, int *ecount, int **elist),
    collect_edge_lengths (int ecount, int *elist, CCdatagroup *dat,
        int **elen);

static void
    randcycle (int ncount, int *cyc, CCdatagroup *dat, double *val,
        CCrandstate *rstate);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)


int CCedgegen_edges (CCedgegengroup *plan, int ncount, CCdatagroup *dat,
        double *wcoord, int *ecount, int **elist, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int i;
    int total, onlist;
    CCkdtree kt;
    int built_a_tree = 0;
    double szeit = CCutil_zeit ();
    int norm;
    tabledat td;
    CCptrworld intptr_world;

    td.table = (intptr **) NULL;
    td.tabletotal = 0;
    CCptrworld_init (&intptr_world);
    td.intptr_world = &intptr_world;
    
    *ecount = 0;
    *elist = (int *) NULL;

    if (ncount < 3) {
        CC_FPRINTF(stderr, "Cannot run CCedgegen_edges in an %d node graph\n",
                 ncount);
        return 1;
    }

    td.table = CC_SAFE_MALLOC (ncount, intptr *);
    if (!td.table) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        td.table[i] = (intptr *) NULL;
    td.tabletotal = 0;

    if (plan->random) {
        rval =  call_random_edge (ncount, dat, plan->random, &td, silent,
                                  rstate);
        if (rval) {
            CC_FPRINTF(stderr, "call_random_edge failed\n"); goto CLEANUP;
        }
    }

    if (plan->nearest) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_nearest (ncount, dat, wcoord, plan->nearest, &kt, &td,
                          silent, rstate)) {
            CC_FPRINTF(stderr, "call_nearest failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->quadnearest) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_quadnearest (ncount, dat, wcoord, plan->quadnearest, &kt, &td,
                              silent, rstate)) {
            CC_FPRINTF(stderr, "call_quadnearest failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->delaunay) {
        if (call_delaunay (ncount, dat, &td, silent)) {
            CC_FPRINTF(stderr, "call_delaunay failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    CCutil_dat_getnorm (dat, &norm);
    if (plan->mlinkern) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
            CC_PRINTF("Cannot run matching Lin-Kernighan with this norm\n");
            CC_FFLUSH(stdout);
        } else {
            if (call_mlinkern (ncount, dat, &kt, plan->mlinkern, &td,
                               silent, rstate)) {
                CC_FPRINTF(stderr, "call_mlinkern failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
    }
    if (plan->tour.random_count) {
        if (call_random_tour (ncount, dat, plan->tour.random_count, &td,
                              silent, rstate)) {
            CC_FPRINTF(stderr, "call_random_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.nearest_count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_nearest_tour (ncount, dat, plan->tour.nearest_count, &kt,
                               &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_nearest_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.greedy) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_greedy_tour (ncount, dat, &kt, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_greedy_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.boruvka) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_boruvka_tour (ncount, dat, &kt, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_boruvka_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.qboruvka) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_qboruvka_tour (ncount, dat, &kt, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_qboruvka_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.twoopt_count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_twoopt_tour (ncount, dat, &kt, plan->tour.twoopt_count,
                              0, 0, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_twoopt_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.twoopt5_count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_twoopt_tour (ncount, dat, &kt, plan->tour.twoopt5_count,
                              1, 0, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_twoopt_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->tour.threeopt_count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_twoopt_tour (ncount, dat, &kt, plan->tour.threeopt_count,
                              0, 1, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_threeopt_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->linkern.count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_linkern (ncount, dat, &kt, plan, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_linkern failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->want_tree) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_spanning_tree (ncount, dat, wcoord, &kt, &td, silent,
                                rstate)) {
            CC_FPRINTF(stderr, "call_spanning_tree failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->nearest_twomatch_count) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_nearest_twomatch (ncount, dat, plan->nearest_twomatch_count,
                                   &kt, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_nearest_twomatch failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->f2match.wantit) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_f2match (ncount, dat, &kt, plan->f2match.priced,
                          plan->f2match.basic, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call_f2match failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (plan->f2match_nearest.number) {
        if (!built_a_tree) {
            if (call_kdtree_build (&kt, ncount, dat, wcoord, &built_a_tree,
                                   silent, rstate)) {
                CC_FPRINTF(stderr, "call_kdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (call_f2match_nearest (ncount, dat, &kt,
               plan->f2match_nearest.number, plan->f2match_nearest.priced,
               plan->f2match_nearest.basic, &td, silent, rstate)) {
            CC_FPRINTF(stderr, "call f2match_nearest failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (!silent) {
        CC_PRINTF("Edgegen total edges: %d (%.2f seconds)\n", td.tabletotal,
                CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    rval = collect_table_edges (&td, ncount, ecount, elist);
    if (rval) {
        CC_FPRINTF(stderr, "collect_table_edges failed\n");
        goto CLEANUP;
    }
        
    if (intptr_check_leaks (&intptr_world, &total, &onlist)) {
        CC_FPRINTF(stderr, "WARNING: %d outstanding intptrs in kdnear\n",
                 total - onlist);
    }

CLEANUP:

    if (built_a_tree)
        CCkdtree_free (&kt);
    CCptrworld_delete (&intptr_world);

    CC_IFFREE (td.table, intptr *);

    return rval;
}

static int call_kdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,
        double *wcoord, int *built_a_tree, int silent, CCrandstate *rstate)
{
    double tzeit;
    int norm;

    *built_a_tree = 0;
    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        tzeit = CCutil_zeit ();
        if (CCkdtree_build (kt, ncount, dat, wcoord, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_build failed\n");
            return 1;
        }
        if (!silent) {
            CC_PRINTF("Built CCkdtree: %.2f (seconds)\n", CCutil_zeit () - tzeit);
            CC_FFLUSH(stdout);
        }
        *built_a_tree = 1;
    }
    return 0;
}

void CCedgegen_init_edgegengroup (CCedgegengroup *plan)
{
    plan->linkern.count = 0;
    plan->linkern.quadnearest = 0;
    plan->linkern.nearest = 0;
    plan->linkern.greedy_start = 0;
    plan->linkern.random_start = 0;
    plan->linkern.nearest_start = 0;
    plan->linkern.boruvka_start = 0;
    plan->linkern.qboruvka_start = 0;
    plan->linkern.nkicks = 0;

    plan->tour.twoopt_count = 0;
    plan->tour.twoopt5_count = 0;
    plan->tour.threeopt_count = 0;
    plan->tour.greedy = 0;
    plan->tour.boruvka = 0;
    plan->tour.qboruvka = 0;
    plan->tour.nearest_count = 0;
    plan->tour.random_count = 0;

    plan->f2match.wantit = 0;
    plan->f2match.basic = 0;
    plan->f2match.priced = 0;

    plan->f2match_nearest.number = 0;
    plan->f2match_nearest.basic = 0;
    plan->f2match_nearest.priced = 0;

    plan->random = 0;
    plan->nearest = 0;
    plan->quadnearest = 0;
    plan->want_tree = 0;
    plan->nearest_twomatch_count = 0;

    plan->delaunay = 0;
    plan->mlinkern = 0;
}

static int call_random_edge (int ncount, CCdatagroup *dat, int ecount,
       tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal, rval = 0;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    
    rval = CCutil_genedgelist (ncount, ecount, &elist, &elen, dat, 0, rstate);
    if (rval) {
        CC_FPRINTF(stderr, "CCutil_genedgelist failed\n"); goto CLEANUP;
    }
 
    rval = put_list_in_table (td, ecount, elist);
    if (rval) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n"); goto CLEANUP;
    }

    if (silent) {
        CC_PRINTF("Random added %d edges (%.2f seconds)\n",
                 td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }
    
CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    return rval;
}

static int call_nearest (int ncount, CCdatagroup *dat, double *wcoord,
        int nearnum, CCkdtree *kt, tabledat *td, int silent,
        CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int tcount = 0;
    int *tlist = (int *) NULL;


    if (work_nearest (kt, ncount, nearnum, dat, wcoord, &tcount, &tlist,
                      silent, rstate)) {
        CC_FPRINTF(stderr, "work_nearest failed\n");
        return 1;
    }

    if (put_list_in_table (td, tcount, tlist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_IFFREE (tlist, int);
        return 1;
    }

    if (!silent) {
        CC_PRINTF("Nearest added %d edges (%.2f seconds)\n",
                     td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tlist, int);

    return 0;
}

static int work_nearest (CCkdtree *kt, int ncount, int nearnum,
        CCdatagroup *dat, double *wcoord, int *ecount, int **elist,
        int silent, CCrandstate *rstate)
{
    int norm;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_k_nearest (kt, ncount, nearnum, dat, wcoord,
                               1, ecount, elist, silent, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_k_nearest failed\n");
            return 1;
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (CCedgegen_x_k_nearest (ncount, nearnum, dat, wcoord, 1, ecount,
                                   elist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_k_nearest failed\n");
            return 1;
        }
    } else {
        if (CCedgegen_junk_k_nearest (ncount, nearnum, dat, wcoord, 1, ecount,
                                      elist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_k_nearest failed\n");
            return 1;
        }
    }
    return 0;
}


static int call_quadnearest (int ncount, CCdatagroup *dat, double *wcoord,
        int nearnum, CCkdtree *kt, tabledat *td, int silent,
        CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int tcount = 0;
    int *tlist = (int *) NULL;

    if (work_quadnearest (kt, ncount, nearnum, dat, wcoord, &tcount, &tlist,
                          silent, rstate)) {
        CC_FPRINTF(stderr, "work_quadnearest failed\n");
        return 1;
    }

    if (put_list_in_table (td, tcount, tlist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_IFFREE (tlist, int);
        return 1;
    }

    if (!silent) {
        CC_PRINTF("Quad Nearest added %d edges (%.2f seconds)\n",
                     td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tlist, int);

    return 0;
}

static int work_quadnearest (CCkdtree *kt, int ncount, int nearnum,
        CCdatagroup *dat, double *wcoord, int *ecount, int **elist,
        int silent, CCrandstate *rstate)
{
    int norm;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_quadrant_k_nearest (kt, ncount, nearnum, dat, wcoord,
                               1, ecount, elist, silent, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_quadrant_k_nearest failed\n");
            return 1;
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (CCedgegen_x_quadrant_k_nearest (ncount, nearnum, dat, wcoord, 1,
                                  ecount, elist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_quadrant_k_nearest failed\n");
            return 1;
        }
    } else {
        if (!silent) {
            CC_PRINTF("Cannot run quadrant nearest with JUNK norms\n");
            CC_PRINTF("Trying %d-nearest instead\n", 2 * nearnum);
            CC_FFLUSH(stdout);
        }
        if (CCedgegen_junk_k_nearest (ncount, 2 * nearnum, dat, wcoord, 1,
                            ecount, elist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_k_nearest failed\n");
            return 1;
        }
    }
    return 0;
}

static int call_random_tour (int ncount, CCdatagroup *dat, int number,
        tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int round;
    int *tour = (int *) NULL;
    double tzeit, val;
    int k;

    if (!silent) {
        CC_PRINTF("Generate %d Random Tours\n", number); CC_FFLUSH(stdout);
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour)
        return 1;

    for (round = 0; round < number; round++) {
        k = td->tabletotal;
        tzeit = CCutil_zeit ();
        randcycle (ncount, tour, dat, &val, rstate);
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            CC_FREE (tour, int);
            return 1;
        }
        if (!silent) {
            CC_PRINTF("  Random tour %d: %.0f, added %d edges (%.2f seconds)\n",
                     round, val, td->tabletotal - k, CCutil_zeit () - tzeit);
            CC_FFLUSH(stdout);
        }
    }

    if (!silent) {
        CC_PRINTF("  TOTAL: Random tours added %d edges (%.2f seconds)\n",
                     td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tour, int);
    return 0;
}


static int call_nearest_tour (int ncount, CCdatagroup *dat, int number,
        CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int round;
    int *tour = (int *) NULL;
    double tzeit, val;
    int k;

    if (!silent) {
        CC_PRINTF("Generate %d Nearest Neighbor Tours\n", number);
        CC_FFLUSH(stdout);
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour)
        return 1;

    for (round = 0; round < number; round++) {
        k = td->tabletotal;
        tzeit = CCutil_zeit ();
        if (work_nearest_tour (kt, ncount, CCutil_lprand (rstate) % ncount,
                                        dat, tour, &val, silent, rstate)) {
            CC_FPRINTF(stderr, "work_nearest_tour failed\n");
            CC_FREE (tour, int);
            return 1;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            CC_FREE (tour, int);
            return 1;
        }
        if (!silent) {
            CC_PRINTF("  NN tour %d: %.0f, added %d edges (%.2f seconds)\n",
                     round, val, td->tabletotal - k, CCutil_zeit () - tzeit);
            CC_FFLUSH(stdout);
        }
    }

    if (!silent) {
        CC_PRINTF("  TOTAL: Nearest tours added %d edges (%.2f seconds)\n",
                     td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tour, int);
    return 0;
}

static int work_nearest_tour (CCkdtree *kt, int ncount, int start,
        CCdatagroup *dat, int *tour, double *val, int silent,
        CCrandstate *rstate)
{
    int norm;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_nearest_neighbor_tour (kt, ncount, start, dat, tour,
                                            val, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_nearest_neighbor_tour failed\n");
            return 1;
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (CCedgegen_x_nearest_neighbor_tour (ncount, start, dat, tour,
                                               val)) {
            CC_FPRINTF(stderr, "CCedgegen_x_nearest_neighbor_tour failed\n");
            return 1;
        }
    } else {
        if (CCedgegen_junk_nearest_neighbor_tour (ncount, start, dat, tour,
                                                  val, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_nearest_neighbor_tour failed\n");
            return 1;
        }
    }
    return 0;
}

static int call_greedy_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int *tour = (int *) NULL;
    int tempcount;
    int *templist = (int *) NULL;
    double val;
    int t, norm;
    int rval = 0;

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        CC_FPRINTF(stderr, "Out of memory in call_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_greedy_tour (kt, ncount, dat, tour, &val, silent,
                                  rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_greedy_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CC_PRINTF("Greedy tour: %.0f, added %d edges (%.2f seconds)\n",
                     val, td->tabletotal - current, CCutil_zeit () - szeit);
            CC_FFLUSH(stdout);
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (CCedgegen_x_quadrant_k_nearest (ncount, 2, dat, (double *) NULL,
                1, &tempcount, &templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_quadrant_k_nearest failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCedgegen_x_greedy_tour (ncount, dat, tour, &val, tempcount,
                                       templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_greedy_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if (ncount < 9) {
            t = ncount - 1;
        } else {
            t = 8;
        }
        if (CCedgegen_junk_k_nearest (ncount, t, dat, (double *) NULL,
                1, &tempcount, &templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_k_nearest failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCedgegen_junk_greedy_tour (ncount, dat, tour, &val, tempcount,
                                          templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_greedy_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

 CLEANUP:
    CC_IFFREE (tour, int);
    CC_IFFREE (templist, int);
    return rval;
}

static int call_boruvka_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int *tour = (int *) NULL;
    double val;
    int norm;

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour)
        return 1;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_boruvka_tour (kt, ncount, dat, tour, &val, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_boruvka_tour failed\n");
            CC_FREE (tour, int);
            return 1;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            CC_FREE (tour, int);
            return 1;
        }
        if (!silent) {
            CC_PRINTF("Boruvka tour: %.0f, added %d edges (%.2f seconds)\n",
                     val, td->tabletotal - current, CCutil_zeit () - szeit);
            CC_FFLUSH(stdout);
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (!silent) {
            CC_PRINTF("No X_NORM boruvka tours, using nearest neighbor\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (tour, int);
        return call_nearest_tour (ncount, dat, 1, kt, td, silent, rstate);
    } else {
        if (!silent) {
            CC_PRINTF("No JUNK_NORM boruvka tours, using nearest neighbor\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (tour, int);
        return call_nearest_tour (ncount, dat, 1, kt, td, silent, rstate);
    }

    CC_IFFREE (tour, int);
    return 0;
}

static int call_qboruvka_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int *tour = (int *) NULL;
    int tempcount;
    int *templist = (int *) NULL;
    double val;
    int norm;
    int rval;

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        CC_FPRINTF(stderr, "Out of memory in call_qboruvka_tour\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_qboruvka_tour (kt, ncount, dat, tour, &val, rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_qboruvka_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CC_PRINTF("Quick boruvka tour: %.0f, added %d edges (%.2f seconds)\n",
                     val, td->tabletotal - current, CCutil_zeit () - szeit);
            CC_FFLUSH(stdout);
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (CCedgegen_x_quadrant_k_nearest (ncount, 2, dat, (double *) NULL,
                1, &tempcount, &templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_quadrant_k_nearest failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCedgegen_x_qboruvka_tour (ncount, dat, tour, &val, tempcount,
                                       templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_x_qboruvka_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if (CCedgegen_junk_k_nearest (ncount, 8, dat, (double *) NULL,
                1, &tempcount, &templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_k_nearest failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCedgegen_junk_qboruvka_tour (ncount, dat, tour, &val, tempcount,
                                          templist, silent)) {
            CC_FPRINTF(stderr, "CCedgegen_junk_qboruvka_tour failed\n");
            rval = 1; goto CLEANUP;
        }
        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

 CLEANUP:
    CC_IFFREE (tour, int);
    CC_IFFREE (templist, int);
    return rval;
}

static int call_twoopt_tour (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int number, int two_and_a_half, int use_3opt, tabledat *td,
        int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    double ival, val, tzeit;
    int k, round, current = td->tabletotal;
    int *tour1 = (int *) NULL, *tour2 = (int *) NULL;
    int norm;

    if (!silent) {
        if (use_3opt)
            CC_PRINTF("Generate %d 3OPT Tours from Nearest Neighbor\n", number);
        else if (two_and_a_half)
            CC_PRINTF("Generate %d 2.5OPT Tours from Nearest Neighbor\n", number);
        else
            CC_PRINTF("Generate %d 2OPT Tours from Nearest Neighbor\n", number);
        CC_FFLUSH(stdout);
    }

    tour1 = CC_SAFE_MALLOC (ncount, int);
    if (!tour1)
        return 1;
    tour2 = CC_SAFE_MALLOC (ncount, int);
    if (!tour2) {
        CC_FREE (tour1, int);
        return 1;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        for (round = 0; round < number; round++) {
            k = td->tabletotal;
            tzeit = CCutil_zeit ();
            if (work_nearest_tour (kt, ncount, CCutil_lprand (rstate) % ncount,
                                        dat, tour1, &val, silent, rstate)) {
                CC_FPRINTF(stderr, "work_nearest_tour failed\n");
                CC_FREE (tour1, int);
                CC_FREE (tour2, int);
                return 1;
            }
            ival = val;
            if (use_3opt) {
                if (CCkdtree_3opt_tour (kt, ncount, dat, tour1, tour2, &val,
                                        1, rstate)) {
                    CC_FPRINTF(stderr, "CCkdtree_3opt_tour failed\n");
                    CC_FREE (tour1, int);
                    CC_FREE (tour2, int);
                    return 1;
                }

            } else {
                if (CCkdtree_twoopt_tour (kt, ncount, dat, tour1, tour2, &val,
                                          two_and_a_half, 1, rstate)) {
                    CC_FPRINTF(stderr, "CCkdtree_twoopt_tour failed\n");
                    CC_FREE (tour1, int);
                    CC_FREE (tour2, int);
                    return 1;
                }
            }
            if (put_tour_in_table (td, ncount, tour2)) {
                CC_FPRINTF(stderr, "put_tour_in_table failed\n");
                CC_FREE (tour1, int);
                CC_FREE (tour2, int);
                return 1;
            }
            if (!silent) {
                if (use_3opt)
                    CC_PRINTF("  3OPT tour %d (from %.0f): %.0f, added %d edges (%.2f sec)\n",
                            round, ival, val, td->tabletotal - k,
                            CCutil_zeit () - tzeit);
                else if (two_and_a_half)
                    CC_PRINTF("  2.5OPT tour %d (from %.0f): %.0f, added %d edges (%.2f sec)\n",
                            round, ival, val, td->tabletotal - k,
                            CCutil_zeit () - tzeit);
                else
                    CC_PRINTF("  2OPT tour %d (from %.0f): %.0f, added %d edges (%.2f sec)\n",
                            round, ival, val, td->tabletotal - k,
                            CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (!silent) {
            if (use_3opt)
                CC_PRINTF("No X_NORM three-opt, using nearest neighbor\n");
            else
                CC_PRINTF("No X_NORM two-opt, using nearest neighbor\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (tour1, int);
        CC_FREE (tour2, int);
        return call_nearest_tour (ncount, dat, number, kt, td, silent, rstate);
    } else {
        if (!silent) {
            if (use_3opt)
                CC_PRINTF("No JUNK_NORM three-opt, using nearest neighbor\n");
            else
                CC_PRINTF("No JUNK_NORM two-opt, using nearest neighbor\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (tour1, int);
        CC_FREE (tour2, int);
        return call_nearest_tour (ncount, dat, number, kt, td, silent, rstate);
    }

    if (!silent) {
        if (use_3opt)
            CC_PRINTF("  TOTAL: 3-opt tours added %d edges (%.2f seconds)\n",
                         td->tabletotal - current, CCutil_zeit () - szeit);
        else if (two_and_a_half)
            CC_PRINTF("  TOTAL: 2.5-opt tours added %d edges (%.2f seconds)\n",
                        td->tabletotal - current, CCutil_zeit () - szeit);
        else
            CC_PRINTF("  TOTAL: 2-opt tours added %d edges (%.2f seconds)\n",
                         td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tour1, int);
    CC_IFFREE (tour2, int);

    return 0;
}

static int call_linkern (int ncount, CCdatagroup *dat, CCkdtree *kt,
        CCedgegengroup *plan, tabledat *td, int silent, CCrandstate *rstate)
{
    int *tour = (int *) NULL, *gtour = (int *) NULL, *itour = (int *) NULL;
    double val, gval, ival;
    double szeit = CCutil_zeit ();
    double tzeit;
    int i, k, round;
    int current = td->tabletotal;
    int ecount = 0, *elist = (int *) NULL;
    int rval = 0;
    int norm;
    
    if (!silent) {
        CC_PRINTF("Generate %d Linkern Tours (", plan->linkern.count);
        if (plan->linkern.greedy_start)
            CC_PRINTF("Greedy, ");
        else if (plan->linkern.boruvka_start)
            CC_PRINTF("Boruvka, ");
        else if (plan->linkern.qboruvka_start)
            CC_PRINTF("Quick Boruvka, ");
        else if (plan->linkern.random_start)
            CC_PRINTF("Random, ");
        else
            CC_PRINTF("Nneigh, ");
        CC_PRINTF("%d kicks, ", plan->linkern.nkicks);
    
        if (plan->linkern.nearest == 0) {
            CC_PRINTF("Quad-%d Edgeset)\n", (plan->linkern.quadnearest ?
                               plan->linkern.quadnearest : 3));
        } else {
            if (plan->linkern.quadnearest == 0) {
                CC_PRINTF("Near-%d Edgeset)\n", plan->linkern.nearest);
            } else {
                CC_PRINTF("Quad-%d + Near-%d Edgeset)\n",
                    plan->linkern.quadnearest, plan->linkern.nearest);
            }
        }
    }

    /* Build an initial edgeset */

    if (plan->linkern.nearest == 0 && plan->linkern.quadnearest == 0) {
        if (work_quadnearest (kt, ncount, 3, dat, (double *) NULL,
                              &ecount, &elist, silent, rstate)) {
            CC_FPRINTF(stderr, "work_quadnearest failed\n");
            return 1;
        }
    } else {
        if (plan->linkern.nearest == 0) {
            if (work_quadnearest (kt, ncount, plan->linkern.quadnearest,
                     dat, (double *) NULL, &ecount, &elist, silent, rstate)) {
                CC_FPRINTF(stderr, "work_quadnearest failed\n");
                return 1;
            }
        } else if (plan->linkern.quadnearest == 0) {
            if (work_nearest (kt, ncount, plan->linkern.nearest,
                              dat, (double *) NULL, &ecount, &elist, 
                              silent, rstate)) {
                CC_FPRINTF(stderr, "work_nearest failed\n");
                return 1;
            }
        } else {
            tabledat tab;
            int tcount, *tlist = (int *) NULL;

            tab.table = (intptr **) NULL;
            tab.tabletotal = 0;
            tab.intptr_world = td->intptr_world;
            
            tab.table = CC_SAFE_MALLOC (ncount, intptr *);
            if (!tab.table)
                return 1;
            for (i = 0; i < ncount; i++)
                tab.table[i] = (intptr *) NULL;

            if (work_quadnearest (kt, ncount, plan->linkern.quadnearest,
                     dat, (double *) NULL, &tcount, &tlist, silent, rstate)) {
                CC_FPRINTF(stderr, "work_quadnearest failed\n");
                return 1;
            }
            if (put_list_in_table (&tab, tcount, tlist)) {
                CC_FPRINTF(stderr, "put_list_in_table failed\n");
                CC_FREE (tab.table, intptr *);
                CC_IFFREE (tlist, int);
                return 1;
            }

            if (work_nearest (kt, ncount, plan->linkern.nearest, dat,
                    (double *) NULL, &tcount, &tlist, silent, rstate)) {
                CC_FPRINTF(stderr, "work_nearest failed\n");
                return 1;
            }
            if (put_list_in_table (&tab, tcount, tlist)) {
                CC_FPRINTF(stderr, "put_list_in_table failed\n");
                CC_FREE (tab.table, intptr *);
                CC_IFFREE (tlist, int);
                return 1;
            }
            CC_IFFREE (tlist, int);

            if (collect_table_edges (&tab, ncount, &ecount, &elist)) {
                CC_FREE (tab.table, intptr *);
                return 1;
            }
            CC_FREE (tab.table, intptr *);
        }
    }
    if (!silent) {
        CC_PRINTF("Initial Edgeset: %d\n", ecount); CC_FFLUSH(stdout);
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (plan->linkern.greedy_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCkdtree_greedy_tour (kt, ncount, dat, gtour, &gval, silent,
                                      rstate)) {
                CC_FPRINTF(stderr, "CCkdtree_greedy_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Greedy tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        } else if (plan->linkern.boruvka_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCkdtree_boruvka_tour (kt, ncount, dat, gtour, &gval, rstate)) {
                CC_FPRINTF(stderr, "CCkdtree_boruvka_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Boruvka tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        } else if (plan->linkern.qboruvka_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCkdtree_qboruvka_tour (kt, ncount, dat, gtour, &gval, rstate)) {
                CC_FPRINTF(stderr, "CCkdtree_qboruvka_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Quick Boruvka tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (plan->linkern.greedy_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCedgegen_x_greedy_tour (ncount, dat, gtour, &gval, ecount,
                                         elist, silent)) {
                CC_FPRINTF(stderr, "CCedgegen_x_greedy_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Greedy tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        } else if (plan->linkern.qboruvka_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCedgegen_x_qboruvka_tour (ncount, dat, gtour, &gval,
                                           ecount, elist, silent)) {
                CC_FPRINTF(stderr, "CCedgegen_x_qboruvka_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Quick Boruvka tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        }
    } else {
        if (plan->linkern.greedy_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCedgegen_junk_greedy_tour (ncount, dat, gtour, &gval, ecount,
                                            elist, silent)) {
                CC_FPRINTF(stderr, "CCedgegen_x_greedy_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Greedy tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        } else if (plan->linkern.qboruvka_start) {
            gtour = CC_SAFE_MALLOC (ncount, int);
            if (!gtour) {
                rval = 1;
                goto CLEANUP;
            }
            tzeit = CCutil_zeit ();
            if (CCedgegen_junk_qboruvka_tour (ncount, dat, gtour, &gval,
                                              ecount, elist, silent)) {
                CC_FPRINTF(stderr, "CCedgegen_x_qboruvka_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!silent) {
                CC_PRINTF("Quick Boruvka tour: %.0f (%.2f seconds)\n",
                        gval, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        }
    }


    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        rval = 1;
        goto CLEANUP;
    }
    itour = CC_SAFE_MALLOC (ncount, int);
    if (!itour) {
        rval = 1;
        goto CLEANUP;
    }

    for (round = 0; round < plan->linkern.count; round++) {
        tzeit = CCutil_zeit ();
        k = td->tabletotal;
        if (gtour != (int *) NULL) {
            for (i = 0; i < ncount; i++)
                itour[i] = gtour[i];
            val = gval;
            ival = gval;
        } else if (plan->linkern.random_start) {
            randcycle (ncount, itour, dat, &val, rstate);
            ival = val;
        } else {
            if (work_nearest_tour (kt, ncount, CCutil_lprand (rstate) % ncount,
                                        dat, itour, &val, silent, rstate)) {
                CC_FPRINTF(stderr, "work_nearest_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            ival = val;
        }

        if (CClinkern_tour (ncount, dat, ecount, elist, 100000000,
             plan->linkern.nkicks, itour, tour, &val, 1, -1.0, -1.0,
             (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate)) {
            CC_FPRINTF(stderr, "CClinkern_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }

        if (put_tour_in_table (td, ncount, tour)) {
            CC_FPRINTF(stderr, "put_tour_in_table failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (!silent) {
            CC_PRINTF("  LK tour %d (from %.0f): %.0f, added %d edges (%.2f sec)\n",
                 round, ival, val, td->tabletotal - k, CCutil_zeit () - tzeit);
            CC_FFLUSH(stdout);
        }
    }

    if (!silent) {
        CC_PRINTF("  TOTAL: Linkern tours added %d edges (%.2f seconds)\n",
                     td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

CLEANUP:

    CC_IFFREE (tour, int);
    CC_IFFREE (itour, int);
    CC_IFFREE (gtour, int);
    CC_IFFREE (elist, int);

    return rval;
}

static int call_f2match (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int priceit, int basic, tabledat *td, int silent, CCrandstate *rstate)
{
    int *mat = (int *) NULL;
    int ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int i;
    double val;

    if (f2match_initial_edgeset (ncount, dat, kt, &ecount, &elist, &elen,
                                 td, silent, rstate)) {
        CC_FPRINTF(stderr, "f2match_initial_edgeset failed\n");
        return 1;
    }

    mat = CC_SAFE_MALLOC((6 * ncount) + 1, int);
    if (!mat) {
        CC_FREE (elist, int);
        CC_FREE (elen, int);
        return 1;
    }

    if (priceit)
        i = CCfmatch_fractional_2match (ncount, ecount, elist, elen, dat, &val,
                    mat, (int *) NULL, (int *) NULL, basic, silent, rstate);
    else
        i = CCfmatch_fractional_2match (ncount, ecount, elist, elen,
                    (CCdatagroup *) NULL, &val, mat, (int *) NULL,
                    (int *) NULL, basic, silent, rstate);
    if (i) {
        CC_FPRINTF(stderr, "CCfmatch_fractional_2match failed\n");
        CC_FREE (mat, int);
        CC_FREE (elist, int);
        CC_FREE (elen, int);
        return 1;
    }

    i = 0;
    while (mat[i] != -1) {
        if (put_in_table (td, mat[i], mat[i + 1])) {
            CC_FPRINTF(stderr, "put_in_table failed\n");
            CC_FREE (mat, int);
            CC_FREE (elist, int);
            CC_FREE (elen, int);
            return 1;
        }
        i += 3;
    }

    CC_IFFREE (mat, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    return 0;
}

static int call_f2match_nearest (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int number, int priceit, int basic, tabledat *td, int silent,
        CCrandstate *rstate)
{
    int ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int *dual = (int *) NULL;
    int dualmax, i;
    double *dcoord = (double *) NULL;
    double val;
    int rval = 0;
    int current = td->tabletotal;
    double szeit = CCutil_zeit ();

    if (f2match_initial_edgeset (ncount, dat, kt, &ecount, &elist, &elen,
                                 td, silent, rstate)) {
        CC_FPRINTF(stderr, "f2match_initial_edgeset failed\n");
        return 1;
    }

    dual = CC_SAFE_MALLOC (ncount, int);
    if (!dual) {
        rval = 1;
        goto CLEANUP;
    }

    if (priceit)
        i = CCfmatch_fractional_2match (ncount, ecount, elist, elen, dat, &val,
                              (int *) NULL, dual, (int *) NULL, basic, silent,
                              rstate);
    else
        i = CCfmatch_fractional_2match (ncount, ecount, elist, elen,
                 (CCdatagroup *) NULL, &val, (int *) NULL, dual,
                 (int *) NULL, basic, silent, rstate);
    if (i) {
        CC_FPRINTF(stderr, "CCfmatch_fractional_2match failed\n");
        rval = 1;
        goto CLEANUP;
    }
    CC_FREE (elist, int);
    CC_FREE (elen, int);
    ecount = 0;

    dualmax = dual[0];
    for (i = 1; i < ncount; i++) {
        if (dual[i] > dualmax)
            dualmax = dual[i];
    }

    dcoord = CC_SAFE_MALLOC (ncount, double);
    if (!dcoord) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        dcoord[i] = (dualmax - dual[i]) * 0.5;
    CC_FREE (dual, int);

    if (work_nearest (kt, ncount, number, dat, dcoord, &ecount, &elist,
                      silent, rstate)) {
        CC_FPRINTF(stderr, "work_nearest failed\n");
        rval = 1;
        goto CLEANUP;
    }
    if (put_list_in_table (td, ecount, elist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        rval = 1; goto CLEANUP;
    }

    if (!silent) {
        CC_PRINTF("Fractional 2-match Nearest-%d added %d edges (%.2f seconds)\n",
                 number, td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (dual, int);
    CC_IFFREE (dcoord, double);

    return rval;
}

static int f2match_initial_edgeset (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int *ecount, int **elist, int **elen, tabledat *td, int silent,
        CCrandstate *rstate)
{
    tabledat tab;
    int tcount, *tlist = (int *) NULL;
    int *ttour = (int *) NULL;
    int i;
    double tval;

    tab.tabletotal = 0;
    tab.table = (intptr **) NULL;
    tab.intptr_world = td->intptr_world;
    
    *ecount = 0;
    *elist = (int *) NULL;
    *elen = (int *) NULL;

    tab.table = CC_SAFE_MALLOC (ncount, intptr *);
    if (!tab.table)
        return 1;
    for (i = 0; i < ncount; i++)
        tab.table[i] = (intptr *) NULL;

    if (work_quadnearest (kt, ncount, 3, dat, (double *) NULL, &tcount,
                          &tlist, silent, rstate)) {
        CC_FPRINTF(stderr, "work_quadnearest failed\n");
        CC_FREE (tab.table, intptr *);
        return 1;
    }
    if (put_list_in_table (&tab, tcount, tlist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_FREE (tab.table, intptr *);
        CC_FREE (tlist, int);
        return 1;
    }
    CC_IFFREE (tlist, int);

    ttour = CC_SAFE_MALLOC (ncount, int);
    if (!ttour) {
        CC_FREE (tab.table, intptr *);
        return 1;
    }
    if (work_nearest_tour (kt, ncount, CCutil_lprand (rstate) % ncount, dat,
                           ttour, &tval, silent, rstate)) {
        CC_FPRINTF(stderr, "work_nearest_tour failed\n");
        CC_FREE (tab.table, intptr *);
        CC_FREE (ttour, int);
        return 1;
    }
    if (put_tour_in_table (&tab, ncount, ttour)) {
        CC_FPRINTF(stderr, "put_tour_in_table failed\n");
        CC_FREE (tab.table, intptr *);
        CC_FREE (ttour, int);
        return 1;
    }
    CC_FREE (ttour, int);

    if (collect_table_edges (&tab, ncount, ecount, elist)) {
        CC_FPRINTF(stderr, "collect_table_edges failed\n");
        CC_FREE (tab.table, intptr *);
        return 1;
    }
    CC_FREE (tab.table, intptr *);

    if (collect_edge_lengths (*ecount, *elist, dat, elen)) {
        CC_FPRINTF(stderr, "collect_edge_lengths failed\n");
        CC_FREE (*elist, int);
        return 1;
    }

    return 0;
}

static int call_spanning_tree (int ncount, CCdatagroup *dat, double *wcoord,
        CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int *tree = (int *) NULL;
    double val;
    int norm;

    tree = CC_SAFE_MALLOC ((2 * ncount) - 2, int);
    if (!tree)
        return 1;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        if (CCkdtree_prim_spanningtree (kt, ncount, dat, wcoord, tree, &val,
                                        rstate)) {
            CC_FPRINTF(stderr, "CCkdtree_prim_spanningtree failed\n");
            CC_FREE (tree, int);
            return 1;
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (!silent) CC_PRINTF("No X_NORM spanning tree\n");
        CC_FREE (tree, int);
        return 0;
    } else {
        if (!silent) CC_PRINTF("No JUNK_NORM spanning tree\n");
        CC_FREE (tree, int);
        return 0;
    }

    if (put_list_in_table (td, ncount-1, tree)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_FREE (tree, int);
        return 1;
    }

    if (!silent) {
        CC_PRINTF("Spanning tree: %.0f, added %d edges (%.2f seconds)\n",
             val, td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tree, int);

    return 0;
}

static int call_delaunay (int ncount, CCdatagroup *dat, tabledat *td,
        int silent)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int tcount = 0;
    int *tlist = (int *) NULL;
    int norm;

    CCutil_dat_getnorm (dat, &norm);
    if (norm == CC_EUCLIDEAN || norm == CC_EUCLIDEAN_CEIL) {
        if (CCedgegen_delaunay (ncount, dat, 1, &tcount, &tlist)) {
            CC_FPRINTF(stderr, "delaunay failed\n");
            return 1;
        }
    } else {
        CC_PRINTF("No Delaunay triangulation with norm %d\n", norm);
        CC_FFLUSH(stdout);
        return 0;
    }

    if (put_list_in_table (td, tcount, tlist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_IFFREE (tlist, int);
        return 1;
    }

    if (!silent) {
        CC_PRINTF("Delaunay Triangulation added %d edges (%.2f seconds)\n",
                 td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tlist, int);

    return 0;
}

static int call_mlinkern (int ncount, CCdatagroup *dat, CCkdtree *kt,
        int iterations, tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int tcount = 0;
    int *tlist = (int *) NULL;


    if (CCedgegen_mlinkern (ncount, dat, 1, &tcount, &tlist, kt, iterations,
        rstate)) {
        CC_FPRINTF(stderr, "mlinkern failed\n");
        return 1;
    }

    if (put_list_in_table (td, tcount, tlist)) {
        CC_FPRINTF(stderr, "put_list_in_table failed\n");
        CC_IFFREE (tlist, int);
        return 1;
    }

    if (!silent) {
        CC_PRINTF("Matching LinKer added %d edges (%.2f seconds)\n",
                 td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (tlist, int);

    return 0;
}

static int call_nearest_twomatch (int ncount, CCdatagroup *dat, int number,
        CCkdtree *kt, tabledat *td, int silent, CCrandstate *rstate)
{
    double szeit = CCutil_zeit ();
    int current = td->tabletotal;
    int round;
    int *mat = (int *) NULL;
    double tzeit, val;
    int k;
    int norm;

    if (!silent) {
        CC_PRINTF("Generate %d Nearest Neighbor 2-matchings\n", number);
        CC_FFLUSH(stdout);
    }

    mat = CC_SAFE_MALLOC (2 * ncount, int);
    if (!mat)
        return 1;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        for (round = 0; round < number; round++) {
            k = td->tabletotal;
            tzeit = CCutil_zeit ();
            if (CCkdtree_nearest_neighbor_2match (kt, ncount,
                   CCutil_lprand (rstate) % ncount, dat, mat, &val, rstate)) {
                CC_FPRINTF(stderr, "CCkdtree_nearest_neighbor_2match failed\n");
                CC_FREE (mat, int);
                return 1;
            }
            if (put_list_in_table (td, ncount, mat)) {
                CC_FPRINTF(stderr, "put_list_in_table failed\n");
                CC_FREE (mat, int);
                return 1;
            }
            if (!silent) {
                CC_PRINTF("  NN 2-mat %d: %.0f, added %d edges (%.2f seconds)\n",
                     round, val, td->tabletotal - k, CCutil_zeit () - tzeit);
                CC_FFLUSH(stdout);
            }
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        if (!silent) {
            CC_PRINTF("No X_NORM NN-2match, using NN-tour instead\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (mat, int);
        return call_nearest_tour (ncount, dat, number, kt, td, silent, rstate);
    } else {
        if (!silent) {
            CC_PRINTF("No JUNK_NORM NN-2match, using NN-tour instead\n");
            CC_FFLUSH(stdout);
        }
        CC_FREE (mat, int);
        return call_nearest_tour (ncount, dat, number, kt, td, silent, rstate);
    }

    if (!silent) {
        CC_PRINTF("  TOTAL: Nearest 2-matchings added %d edges (%.2f seconds)\n",
                 td->tabletotal - current, CCutil_zeit () - szeit);
        CC_FFLUSH(stdout);
    }

    CC_IFFREE (mat, int);
    return 0;
}

static int put_tour_in_table (tabledat *td, int ncount, int *tour)
{
    int i;

    for (i = 1; i < ncount; i++) {
        if (put_in_table (td, tour[i-1], tour[i])) {
            CC_FPRINTF(stderr, "put_in_table failed\n");
            return 1;
        }
    }
    if (put_in_table (td, tour[ncount - 1], tour[0])) {
        CC_FPRINTF(stderr, "put_in_table failed\n");
        return 1;
    }

    return 0;
}

static int put_list_in_table (tabledat *td, int ecount, int *elist)
{
    int i;
    
    for (i = 0; i < ecount; i++) {
        if (put_in_table (td, elist[2 * i], elist[(2 * i) + 1])) {
            CC_FPRINTF(stderr, "put_in_table failed\n");
            return 1;
        }
    }
    return 0;
}

static int put_in_table (tabledat *td, int i, int j)
{
    intptr *ip;

    if (j < i) {
        int temp;
        CC_SWAP(i, j, temp);
    }

    for (ip = td->table[i]; ip; ip = ip->next) {
        if (ip->this == j) {
            return 0;
        }
    }
    if (intptr_listadd(&td->table[i], j, td->intptr_world)) {
        return 1;
    }
    td->tabletotal++;
    return 0;
}

static int collect_table_edges (tabledat *td, int ncount, int *ecount,
        int **elist)
{
    *ecount = 0;
    *elist = (int *) NULL;
    
    if (td->tabletotal) {
        int i = 0;
        int j = 0;
        intptr *ip, *ipnext;

        *elist = CC_SAFE_MALLOC (2 * td->tabletotal, int);
        if (!(*elist)) {
            CC_FPRINTF(stderr, "Out of memory in collect_table_edges\n");
            return 1;
        }
        *ecount = td->tabletotal;
        for (i = 0; i < ncount; i++) {
            for (ip = td->table[i]; ip; ip = ipnext) {
                ipnext =  ip->next;
                (*elist)[j++] = i;
                (*elist)[j++] = ip->this;
                intptrfree (td->intptr_world, ip);
            }
            td->table[i] = (intptr *) NULL;
        }
    }
    return 0;
}

static int collect_edge_lengths (int ecount, int *elist, CCdatagroup *dat,
        int **elen)
{
    int i;

    *elen = CC_SAFE_MALLOC (ecount, int);
    if ((*elen) == (int *) NULL) {
        CC_FPRINTF(stderr, "Out of memory in collect_edge_lengths\n");
        return 1;
    }
    for (i=0; i<ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen (elist[2*i],elist[2*i+1], dat);
    }
    return 0;
}

static void randcycle (int ncount, int *cyc, CCdatagroup *dat, double *val,
        CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++)
        cyc[i] = i;

    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }

    *val = CCutil_dat_edgelen (cyc[ncount - 1], cyc[0], dat);
    for (i = 1; i < ncount; i++)
        (*val) += CCutil_dat_edgelen (cyc[i - 1], cyc[i], dat);
}


int CCedgegen_read (char *egname, CCedgegengroup *plan)
{
    char buf[256];
    char area[256];
    char key[256];
    char field[256];
    char *p;
    FILE *in;

    CCedgegen_init_edgegengroup (plan);

    in = fopen (egname, "r");
    if (!in) {
        perror (egname);
        CC_FPRINTF(stderr, "can't open %s for input\n", egname);
        return 1;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        while (*p != '\0') {
            if (*p == ':')
                *p = ' ';
             p++;
        }
        p = buf;
        if (sscanf (p, "%s", area) != EOF) {
            p += strlen (area);
            while (*p == ' ')
                p++;
            if (!strcmp (area, "EDGEGEN")) {
                if (sscanf (p, "%s", key) == EOF) {
                    CC_FPRINTF(stderr, "ERROR in EDGEGEN LINE - no keyword\n");
                    return 1;
                }
                if (!strcmp (key, "RANDOM")) {
                    p += strlen (key);
                    while (*p ==  ' ')
                       p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->random = atoi (field);
                    } else {
                        CC_PRINTF("RANDOM count not given, using 0\n");
                        plan->random = 0;
                    }
                } else if (!strcmp (key, "NEAREST")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->nearest = atoi (field);
                    } else {
                        CC_PRINTF("NEAREST count not given, using 1\n");
                        plan->nearest = 1;
                    }
                } else if (!strcmp (key, "QUADNEAREST")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->quadnearest = atoi (field);
                    } else {
                        CC_PRINTF("QUADNEAREST count not given, using 1\n");
                        plan->quadnearest = 1;
                    }
                } else if (!strcmp (key, "DELAUNAY")) {
                    plan->delaunay = 1;
                } else if (!strcmp (key, "M_LINKERN")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->mlinkern = atoi (field);
                    } else {
                        CC_PRINTF("M_LINKERN count not given, using 1\n");
                        plan->mlinkern = 1;
                    }
                } else if (!strcmp (key, "TREE")) {
                    plan->want_tree = 1;
                } else if (!strcmp (key, "NN_TWOMATCH")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->nearest_twomatch_count = atoi (field);
                    } else {
                        CC_PRINTF("NN_TWOMATCH count not given, using 1\n");
                        plan->nearest_twomatch_count = 1;
                    }
                } else if (!strcmp (key, "GREEDY_TOUR")) {
                    plan->tour.greedy = 1;
                } else if (!strcmp (key, "BORUVKA_TOUR")) {
                    plan->tour.boruvka = 1;
                } else if (!strcmp (key, "QBORUVKA_TOUR")) {
                    plan->tour.qboruvka = 1;
                } else if (!strcmp (key, "NN_TOUR")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->tour.nearest_count = atoi (field);
                    } else {
                        CC_PRINTF("NN_TOUR count not given, using 1\n");
                        plan->tour.nearest_count = 1;
                    }
                } else if (!strcmp (key, "RANDOM_TOUR")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->tour.random_count = atoi (field);
                    } else {
                        CC_PRINTF("RANDOM_TOUR count not given, using 1\n");
                        plan->tour.random_count = 1;
                    }
                } else if (!strcmp (key, "TWOOPT_TOUR")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->tour.twoopt_count = atoi (field);
                    } else {
                        CC_PRINTF("TWOOPT_TOUR count not given, using 1\n");
                        plan->tour.twoopt_count = 1;
                    }
                } else if (!strcmp (key, "TWOOPT5_TOUR")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->tour.twoopt5_count = atoi (field);
                    } else {
                        CC_PRINTF("TWOOPT5_TOUR count not given, using 1\n");
                        plan->tour.twoopt5_count = 1;
                    }
                } else if (!strcmp (key, "THREEOPT_TOUR")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->tour.threeopt_count = atoi (field);
                    } else {
                        CC_PRINTF("THREEOPT_TOUR count not given, using 1\n");
                        plan->tour.threeopt_count = 1;
                    }
                } else if (!strcmp (key, "FRAC_TWOMATCH")) {
                    plan->f2match.wantit = 1;
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    while (sscanf (p, "%s", field) != EOF) {
                        if (!strcmp (field, "BASIC"))
                            plan->f2match.basic = 1;
                        else if (!strcmp (field, "PRICED"))
                            plan->f2match.priced = 1;
                        else
                            CC_PRINTF("Unknown option in FRAC_TWOMATCH\n");
                        p += strlen (field);
                        while (*p == ' ')
                            p++;
                    }
                } else if (!strcmp (key, "FRAC_TWOMATCH_NEAREST")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    while (sscanf (p, "%s", field) != EOF) {
                        if (!strcmp (field, "BASIC"))
                            plan->f2match_nearest.basic = 1;
                        else if (!strcmp (field, "PRICED"))
                            plan->f2match_nearest.priced = 1;
                        else
                            plan->f2match_nearest.number = atoi (field);
                        p += strlen (field);
                        while (*p == ' ')
                            p++;
                    }
                    if (plan->f2match_nearest.number == 0) {
                        CC_PRINTF("FRAC_TWOMATCH_NEAREST count not given, using 1\n");
                        plan->f2match_nearest.number = 1;
                    }
                } else if (!strcmp (key, "LINKERN")) {
                    p += strlen (key);
                    while (*p == ' ')
                        p++;
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->linkern.count = atoi (field);
                        p += strlen (field);
                        while (*p == ' ')
                            p++;
                    } else {
                        CC_PRINTF("LINKERN count not given, using 1\n");
                        plan->linkern.count = 1;
                    }
                    if (sscanf (p, "%s", field) != EOF) {
                        plan->linkern.nkicks = atoi (field);
                        p += strlen (field);
                        while (*p == ' ')
                            p++;
                    } else {
                        CC_PRINTF("LINKERN nkicks not given, using 10\n");
                        plan->linkern.nkicks = 10;
                    }
                    while (sscanf (p, "%s", field) != EOF) {
                        if (!strcmp (field, "GREEDY_START"))
                            plan->linkern.greedy_start = 1;
                        else if (!strcmp (field, "BORUVKA_START"))
                            plan->linkern.boruvka_start = 1;
                        else if (!strcmp (field, "QBORUVKA_START"))
                            plan->linkern.qboruvka_start = 1;
                        else if (!strcmp (field, "RANDOM_START"))
                            plan->linkern.random_start = 1;
                        else if (!strcmp (field, "NN_START"))
                            plan->linkern.nearest_start = 1;
                        else if (!strcmp (field, "NEAREST")) {
                            p += strlen (field);
                            while (*p == ' ')
                                p++;
                            if (sscanf (p, "%s", field) != EOF) {
                                plan->linkern.nearest = atoi (field);
                            } else {
                                CC_PRINTF("LINKERN NEAREST COUNT not given, using 5\n");
                                plan->linkern.nearest = 5;
                                break;
                            }
                        } else if (!strcmp (field, "QUADNEAREST")) {
                            p += strlen (field);
                            while (*p == ' ')
                                p++;
                            if (sscanf (p, "%s", field) != EOF) {
                                plan->linkern.quadnearest = atoi (field);
                            } else {
                                CC_PRINTF("LINKERN QUADNEAREST COUNT not given, using 3\n");
                                plan->linkern.quadnearest = 3;
                                break;
                            }
                        } else {
                            CC_PRINTF("Unknown EDGEGEN LINKERN command %s\n",
                                    field);
                            CC_FFLUSH(stdout);
                        }
                        p += strlen (field);
                        while (*p == ' ')
                            p++;
                    }
                } else {
                    CC_PRINTF("Unknown EDGEGEN command: %s\n", key);
                    CC_FFLUSH(stdout);
                }
            } else {
                CC_PRINTF("Cannot parse command line: %s\n", area);
                CC_FFLUSH(stdout);
            }
        }
    }
    fclose (in);
    CC_PRINTF("\n");

    if (plan->linkern.count) {
        if (!plan->linkern.quadnearest && !plan->linkern.nearest)
            plan->linkern.quadnearest = 3;
        if (!plan->linkern.greedy_start && !plan->linkern.random_start &&
            !plan->linkern.boruvka_start && !plan->linkern.qboruvka_start)
            plan->linkern.nearest_start = 1;
        if (!plan->linkern.nkicks)
            plan->linkern.nkicks = 10;
    }

    CC_PRINTF("Edgegen Request:\n");
    if (plan->nearest)
        CC_PRINTF("  Nearest %d\n", plan->nearest);
    if (plan->quadnearest)
        CC_PRINTF("  Quad-Nearest %d\n", plan->quadnearest);
    if (plan->f2match_nearest.number) {
        CC_PRINTF("  Frac 2-match Nearest %d (", plan->f2match_nearest.number);
        if (plan->f2match_nearest.basic)
            CC_PRINTF("Basic ");
        if (plan->f2match_nearest.priced)
            CC_PRINTF("Priced)\n");
        else
            CC_PRINTF("Not Priced)\n");
    }
    if (plan->delaunay)
        CC_PRINTF("  Delaunay Triangulation\n");
    if (plan->want_tree)
        CC_PRINTF("  Minimum Spanning Tree\n");
    if (plan->nearest_twomatch_count)
        CC_PRINTF("  NN 2-matchings: %d\n", plan->nearest_twomatch_count);
    if (plan->tour.random_count)
        CC_PRINTF("  Random Tours: %d\n", plan->tour.random_count);
    if (plan->tour.nearest_count)
        CC_PRINTF("  NN Tours: %d\n", plan->tour.nearest_count);
    if (plan->tour.greedy)
        CC_PRINTF("  Greedy Tour\n");
    if (plan->tour.boruvka)
        CC_PRINTF("  Boruvka Tour\n");
    if (plan->tour.qboruvka)
        CC_PRINTF("  Quick Boruvka Tour\n");
    if (plan->tour.twoopt_count)
        CC_PRINTF("  2OPT Tours: %d\n", plan->tour.twoopt_count);
    if (plan->tour.twoopt5_count)
        CC_PRINTF("  2.5OPT Tours: %d\n", plan->tour.twoopt5_count);
    if (plan->tour.threeopt_count)
        CC_PRINTF("  3OPT Tours: %d\n", plan->tour.threeopt_count);
    if (plan->linkern.count) {
        CC_PRINTF("  LK Tours: %d (", plan->linkern.count);
        if (plan->linkern.greedy_start)
            CC_PRINTF("Greedy, ");
        else if (plan->linkern.boruvka_start)
            CC_PRINTF("Boruvka, ");
        else if (plan->linkern.qboruvka_start)
            CC_PRINTF("Quick Boruvka, ");
        else if (plan->linkern.random_start)
            CC_PRINTF("Random, ");
        else
            CC_PRINTF("NN, ");
        if (!plan->linkern.nearest) {
            CC_PRINTF("Quad-%d, ", plan->linkern.quadnearest);
        } else {
            if (!plan->linkern.quadnearest) {
                CC_PRINTF("Near-%d, ", plan->linkern.nearest);
            } else {
                CC_PRINTF("Quad-%d + Near-%d, ", plan->linkern.quadnearest,
                                               plan->linkern.nearest);
            }
        }
        CC_PRINTF("%d Kicks)\n", plan->linkern.nkicks);
    }
    if (plan->f2match.wantit) {
        CC_PRINTF("  Frac 2-matching (");
        if (plan->f2match.basic)
            CC_PRINTF("Basic ");
        if (plan->f2match.priced)
            CC_PRINTF("Priced)\n");
        else
            CC_PRINTF("Not Priced)\n");
    }
    if (plan->mlinkern) {
        CC_PRINTF("  LK matchings: %d\n", plan->mlinkern);
    }
    CC_PRINTF("\n");
    CC_FFLUSH(stdout);

    return 0;
}
