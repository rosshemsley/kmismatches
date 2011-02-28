/*
 * sp_km.h
 *
 *  Created on: 14 Sep 2010
 *      Author: ben
 */

#ifndef SP_KM_H_
#define SP_KM_H_

/*Flags that apply to the all match-with-dont-cares problems */
#define SP_KM_FIRST_MATCH_ONLY (1U << 0) //Stop after finding the first match

struct SP_KM_INT_LIST{
	struct SP_KM_INT_LIST *next;
	int i;
	int hammingDistance;
};

struct SP_KM_MATCHING_POSITIONS{
	struct SP_KM_INT_LIST *start;
	struct SP_KM_INT_LIST *iterator;
	struct SP_KM_INT_LIST *end;
};

struct SP_KM_MATCHING_POSITIONS *sp_km_create_new_list_of_matches();

struct SP_KM_INT_LIST *sp_km_newListItem(int i,int hammingDistance);

void sp_km_addToListOfMatches(struct SP_KM_MATCHING_POSITIONS *listOfMatches,int i,int hammingDistance);

void sp_km_freeListOfMatches(struct SP_KM_MATCHING_POSITIONS *listOfMatches);

//#include "sp_km_naive_matcher.h"
// #include "sp_km_unbounded_matcher.h"


#endif /* SP_KM_H_ */
