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
/*                  SEND MESSAGES TO THE BRANCHING BOSS                     */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 15, 1997                                                  */
/*                                                                          */
/*  This program is used to send messages to the boss of a parallel         */
/*  branching run.                                                          */
/*                                                                          */
/*  SEE the short description in usage () for more details.                 */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static char *hostname          = (char *) NULL;
static unsigned short hostport = 0;
static int send_exit           = 0;
static int deadnode            = -1;
static int tellcutboss         = 0;
static int telldomboss         = 0;
static int savepool            = 0;


static int
    parseargs (int ac, char **av);

static void
    usage (char *f);


int main (int ac, char **av)
{
    int rval;
    CC_SFILE *f = (CC_SFILE *) NULL;

#ifdef CC_NETREADY
    rval = parseargs (ac, av);
    if (rval) return 0;

    if ((unsigned int) hostport == 0) {
        if (tellcutboss)      hostport = (unsigned int) CCtsp_CUT_PORT;
        else if (telldomboss) hostport = (unsigned int) CCtsp_DOMINO_PORT;
        else                  hostport = (unsigned int) CCtsp_HOST_PORT;
    }
    
    f = CCutil_snet_open (hostname, hostport);
    if (f == (CC_SFILE *) NULL) {
        CC_FPRINTF(stderr, "Could not open connection to host %s\n", hostname);
        rval = 1; goto CLEANUP;
    }

    if (tellcutboss) {
        if (savepool) {
            rval = CCutil_swrite_char (f, CCtsp_POOL_SAVECUTS);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_char failed\n");
                goto CLEANUP;
            }
        } else if (send_exit) {
            rval = CCutil_swrite_char (f, CCtsp_POOL_EXIT);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_char failed\n");
                goto CLEANUP;
            }
        } else {
            CC_FPRINTF(stderr, "Nothing to send to cutboss\n");
            rval = 1; goto CLEANUP;
        }
    } else if (telldomboss) {
        if (send_exit) {
            rval = CCutil_swrite_char (f, CCtsp_DOMINO_EXIT);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_char failed\n");
                goto CLEANUP;
            }
        } else {
            CC_FPRINTF(stderr, "Nothing to send to domboss\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if (deadnode >= 0) {
            rval = CCutil_swrite_char (f, CCtsp_BBREQ_DEADNODE);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_char failed\n");
                goto CLEANUP;
            }
            rval = CCutil_swrite_int (f, deadnode);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_int failed\n");
                goto CLEANUP;
            }
        } else if (send_exit) {
            rval = CCutil_swrite_char (f, CCtsp_BBREQ_EXIT);
            if (rval) {
                CC_FPRINTF(stderr, "CCutil_swrite_char failed\n");
                goto CLEANUP;
            }
        } else {
            CC_FPRINTF(stderr, "Nothing to send to boss\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCutil_sclose (f);
    f = (CC_SFILE *) NULL;
    if (rval) {
        CC_FPRINTF(stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    
    rval = 0;
  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
#else /* CC_NETREADY */
    CC_FPRINTF(stderr, "Networking disabled\n");
    return 1;
#endif /* CC_NETREADY */
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "cd:Dp:sx", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'c':
            tellcutboss = 1;
            break;
        case 'D':
            telldomboss = 1;
            break;
        case 'd':
            deadnode = atoi(boptarg);
            break;
        case 'p':
            hostport = atoi(boptarg);
            break;
        case 's':
            savepool = 1;
            break;
        case 'x':
            send_exit = 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    hostname = av[boptind++];
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    CC_FPRINTF(stderr, "Usage: %s [-p port] -d deadnode hostname\n", f);
    CC_FPRINTF(stderr, "   or: %s [-p port] -x hostname\n", f);
    CC_FPRINTF(stderr, "   or: %s -c [-p port] -s cuthostname\n", f);
    CC_FPRINTF(stderr, "   or: %s -c [-p port] -x cuthostname\n", f);
    CC_FPRINTF(stderr, "   or: %s -D [-p port] -x domhostname\n", f);
    CC_FPRINTF(stderr, "   -d n tells the boss that node n is dead\n");
    CC_FPRINTF(stderr, "   -s   tells the cutboss to save the pool\n");
    CC_FPRINTF(stderr, "   -x   tells the (cut)boss to exit\n");
    CC_FPRINTF(stderr, "   -p n specifies the port the (cut)boss is listening to\n");
}

