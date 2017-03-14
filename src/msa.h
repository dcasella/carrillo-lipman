#ifndef MSA_H
#define MSA_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "needleman_wunsch.h"
#include "alignment_macros.h"
#include "pqueue.h"


#define MAT_IINDEX(n2, n3, i, j, k) \
  (size_t)((i*n2 + j)*n3 + k)

#define MAT_ILOOKUP(m, n2, n3, i, j, k) \
  m[MAT_IINDEX(n2, n3, i, j, k)]

#define MAT_PINDEX(n2, n3, pt) \
  (size_t)((pt.i*n2 + pt.j)*n3 + pt.k)

#define MAT_PLOOKUP(m, n2, n3, pt) \
  m[MAT_PINDEX(n2, n3, pt)]

typedef struct point {
  int i;
  int j;
  int k;
} point_t;

typedef struct node {
  pqueue_pri_t d;
  point_t      p;
  size_t       pos;
} node_t;

// DEBUG function
void print_node(FILE *out, node_t *n);

static int cmp_pri(pqueue_pri_t next, pqueue_pri_t curr);

static pqueue_pri_t get_pri(node_t *a);

static void set_pri(node_t *a, pqueue_pri_t d);

static size_t get_pos(node_t *a);

static void set_pos(node_t *a, size_t pos);

int point_eq(point_t a, point_t b);

point_t point_inc(point_t p, int i);

void pair_align(char *seq_a, char *seq_b,
                nw_aligner_t *nw, alignment_t *result);

void multi_align(char *seq_a, char *seq_b, char *seq_c, size_t upper_bound);

#endif
