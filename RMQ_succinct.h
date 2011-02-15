

/*******************************************************************************
*
* Ported to C by Ross Hemsley.
*
*
*******************************************************************************/


#ifndef _RMQ_succinct_h_
#define _RMQ_succinct_h_

#define MEM_COUNT


#include <stdlib.h>
#include <limits.h>
#include <stdio.h>


#include <math.h>

typedef int DT;                 // use long for 64bit-version (but take care of fast log!)
typedef unsigned int DTidx;     // for indexing in arrays

/******************************************************************************/
// Types

typedef unsigned char DTsucc;
typedef unsigned short DTsucc2;

/******************************************************************************/

	// liefert RMQ[i,j]
	DTidx query(DTidx i, DTidx j, DT* a, DTidx n);
	void RMQ_succinct(DT* a, DTidx n);
	void FreeRMQ_succinct();

/******************************************************************************/
// Globals
//
// These should be tidied away into a struct or somesuch.
//
/******************************************************************************/

	// table M for the out-of-block queries (contains indices of block-minima)
	DTsucc** M;

	// depth of table M:
	DTidx M_depth;

	// table M' for superblock-queries (contains indices of block-minima)
	DTidx** Mprime;

	// depth of table M':
	DTidx Mprime_depth;

	// type of blocks
	DTsucc2 *type;

	// precomputed in-block queries
	DTsucc** Prec;

	// microblock size
	DTidx s;

	// block size
	DTidx sprime;

	// superblock size
	DTidx sprimeprime;

	// number of blocks (always n/sprime)
	DTidx nb;

	// number of superblocks (always n/sprimeprime)
	DTidx nsb;

	// number of microblocks (always n/s)
	DTidx nmb;
	
#endif
