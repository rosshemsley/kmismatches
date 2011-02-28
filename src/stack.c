#include <stdio.h>
#include <stdlib.h>

#include "stack.h"


/*******************************************************************************
* This is a simple array-based implementation of a stack, implementing
* pop, push, isEmpty, size and empty.
* We use isEmpty to tell whether we can keep popping pointers (because we could
* pop a NULL pointer to the stack). We use the empty function to reset the
* stack to zero without popping any elements. This maintains the memory
* associated with the stack.
*******************************************************************************/

stack *newStack()
{
  stack *s = malloc(sizeof(stack));
  
  s->top   = 0;
  s->slots = 2;
  s->arr = malloc(2*sizeof(int));
  
  return s;
}

/******************************************************************************/

int stackSize(stack *s)
{
  return s->top;
}

/******************************************************************************/
// When we are storing pointers which could be null, we need to have a check 
// to see if the stack is empty.

int isEmpty(stack *s)
{
  return (s->top == 0);
}

/******************************************************************************/
// This will cause the stack to be empty: note, we have NOT freed any memory
// associated with the elements.

void emptyStack(stack *s)
{
  s->top = 0;
}

/******************************************************************************/

void  push(stack *s, int e)
{
  if (s->top >= s->slots)
   {
      // we have to allocate more space
      s->slots *= 2;
      s->arr = realloc(s->arr, (s->slots*sizeof(int)));
      // check that we haven't run out of memory
      if (s->arr == NULL)
      {
         fprintf(stderr, "Error: Out of Memory.\n");
         exit(1);
      }
   }
   // add the element
   s->arr[s->top] = e;
   s->top ++;
}

/******************************************************************************/

int pop(stack *s)
{
  // If the stack is empty
  if (s->top == 0) return 0;  
  s->top--;
  return s->arr[s->top];
}

/******************************************************************************/

int peek(stack *s)
{
  return s->arr[s->top-1];
}

/******************************************************************************/

void freeStack(stack *s)
{
  free(s->arr);
  free(s);
}

/******************************************************************************/
