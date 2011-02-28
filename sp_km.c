/*
 * sp_km.c
 *
 *  Created on: 14 Sep 2010
 *      Author: ben
 */

#include "sp_km.h"
#include <stdlib.h>

struct SP_KM_MATCHING_POSITIONS *sp_km_create_new_list_of_matches(){
	struct SP_KM_MATCHING_POSITIONS *list = malloc(sizeof(struct SP_KM_MATCHING_POSITIONS));
	list->start = list->iterator = list->end = NULL;
	return list;
}

struct SP_KM_INT_LIST *sp_km_newListItem(int i,int hammingDistance){
	struct SP_KM_INT_LIST *listItem = malloc(sizeof(struct SP_KM_INT_LIST));
	listItem->next = NULL;
	listItem->i = i;
	listItem->hammingDistance = hammingDistance;
	return listItem;
}

void sp_km_addToListOfMatches(struct SP_KM_MATCHING_POSITIONS *listOfMatches,int i,int hammingDistance){
	struct SP_KM_INT_LIST *listItem = sp_km_newListItem(i,hammingDistance);
	if(listOfMatches->start == NULL){
		listOfMatches->start = listOfMatches->end = listOfMatches->iterator = listItem;
	}else{
		listOfMatches->end->next = listItem;
		listOfMatches->end = listItem;
	}
}

void sp_km_freeListOfMatches(struct SP_KM_MATCHING_POSITIONS *listOfMatches){
	struct SP_KM_INT_LIST *next;
	while(listOfMatches->start !=NULL){
		next = listOfMatches->start->next;
		free(listOfMatches->start);
		listOfMatches->start = next;
	}
	free(listOfMatches);
}
