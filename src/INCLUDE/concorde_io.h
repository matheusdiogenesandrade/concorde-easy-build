/****************************************************************************/
/*                                                                          */
/*  I/O macros for optional Concorde verbosity (CONCORDE_SILENT)            */
/*                                                                          */
/*  When CONCORDE_SILENT is defined, CC_PRINTF, CC_FFLUSH, CC_FPRINTF       */
/*  expand to no-ops at compile time, eliminating format parsing and        */
/*  syscall cost. Define via -DCONCORDE_SILENT when building.               */
/*                                                                          */
/****************************************************************************/

#ifndef CONCORDE_IO_H
#define CONCORDE_IO_H

#include <stdio.h>

#ifdef CONCORDE_SILENT
#define CC_PRINTF(...) ((void)0)
#define CC_FFLUSH(stream) ((void)0)
#define CC_FPRINTF(stream, ...) ((void)0)
#else
#define CC_PRINTF(...) printf(__VA_ARGS__)
#define CC_FFLUSH(stream) fflush(stream)
#define CC_FPRINTF(stream, ...) fprintf(stream, __VA_ARGS__)
#endif

#endif
