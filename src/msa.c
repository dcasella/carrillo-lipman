#include "msa.h"


#ifdef DEBUG
void print_node(FILE *out, node_t *n) {
  fprintf(out, "  | (%d, %d, %d), D: %lld\n",
          n->p.i, n->p.j, n->p.k, n->d);
}
#endif

static int cmp_pri(pqueue_pri_t next, pqueue_pri_t curr) {
  return (next >= curr);
}

static pqueue_pri_t get_pri(node_t *a) {
  return a->d;
}

static void set_pri(node_t *a, pqueue_pri_t d) {
  a->d = d;
}

static size_t get_pos(node_t *a) {
  return a->pos;
}

static void set_pos(node_t *a, size_t pos) {
  a->pos = pos;
}

int point_eq(point_t a, point_t b) {
  if (a.i != b.i || a.j != b.j || a.k != b.k)
    return 0;
  else
    return 1;
}

point_t point_inc(point_t p, int i) {
  point_t r = p;
  
  if (i == 0) {
    r.k++;
  }
  else if (i == 1) {
    r.j++;
  }
  else if (i == 2) {
    r.k++;
    r.j++;
  }
  else if (i == 3) {
    r.i++;
  }
  else if (i == 4) {
    r.i++;
    r.k++;
  }
  else if (i == 5) {
    r.i++;
    r.j++;
  }
  else if (i == 6) {
    r.i++;
    r.j++;
    r.k++;
  }
  
  return r;
}

void pair_align(char *seq_a, char *seq_b,
                nw_aligner_t *nw, alignment_t *result) {
  // setup
  int  match                =  0;
  int  mismatch             = -1;
  int  gap_extend           = -1;
  int  gap_open             =  0;
  char no_start_gap_penalty =  0;
  char no_end_gap_penalty   =  0;
  char no_gaps_in_a         =  0;
  char no_gaps_in_b         =  0;
  char no_mismatches        =  0;
  char case_sensitive       =  0;
  
  // init scoring
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
  
  // align
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
}

void multi_align(char *seq_a, char *seq_b, char *seq_c, size_t upper_bound) {
  // declarations
  
  // lengths
  size_t n[3];
  n[0] = strlen(seq_a) + 1,
  n[1] = strlen(seq_b) + 1,
  n[2] = strlen(seq_c) + 1;
  
  // nw
  nw_aligner_t *nw[3];
  alignment_t  *result[3];
  for (int i = 0; i < 3; ++i) {
    nw[i]     = needleman_wunsch_new();
    result[i] = alignment_create(256);
  }
  
  // prioprity queue
  pqueue_t *pqueue = pqueue_init(n[0] * n[1] * n[2],
                                 cmp_pri, get_pri, set_pri,
                                 get_pos, set_pos);
  
  // vertices array
  node_t *ns = calloc(n[0] * n[1] * n[2], sizeof(node_t));
  
  // pairwise alignments
  
  pair_align(seq_a, seq_b, nw[0], result[0]);
  pair_align(seq_a, seq_c, nw[1], result[1]);
  pair_align(seq_b, seq_c, nw[2], result[2]);
  
  // set start and target vertices
  
  node_t start, target;
  
  start.d  = 0;
  start.p  = (point_t){ 0, 0, 0 };
  target.d = ULLONG_MAX;
  target.p = (point_t){ n[0] - 1, n[1] - 1, n[2] - 1 };
  
  pqueue_insert(pqueue, &start);
  
  // algorithm
  
  node_t *v = NULL;
  
  while ((v = pqueue_pop(pqueue))) {
    // extract min
#ifdef DEBUG
    printf("---\n");
    printf("POP (%d, %d, %d), D: %lld\n", v->p.i, v->p.j, v->p.k, v->d);
#endif
    
    if (point_eq(v->p, target.p)) {
      // found msa
      break;
    }
    
    // set lower bound
    pqueue_pri_t lower_bound = 0;
    lower_bound += abs(ARR_LOOKUP(nw[0]->match_scores,
                                  nw[0]->score_width,
                                  v->p.i, v->p.j));
    lower_bound += abs(ARR_LOOKUP(nw[1]->match_scores,
                                  nw[1]->score_width,
                                  v->p.i, v->p.k));
    lower_bound += abs(ARR_LOOKUP(nw[2]->match_scores,
                                  nw[2]->score_width,
                                  v->p.j, v->p.k));
    
    if (v->d <= upper_bound - lower_bound) {
      // throwable value
      for (int i = 0; i < 7; ++i) {
        // throw forward
        point_t r = point_inc(v->p, i);
        node_t *w = &MAT_PLOOKUP(ns, n[1], n[2], r);
        
        if (r.i < n[0] && r.j < n[1] && r.k < n[2]) {
          // not out of bounds
#ifdef DEBUG
          printf("FOR (%d, %d, %d)", r.i, r.j, r.k);
#endif
          
          pqueue_pri_t p_vw = 0, p_w;
          
          if (point_eq(w->p, r) == 0) {
            // new node
#ifdef DEBUG
            printf(", D: ∞ - INSERT");
#endif
            
            w->p = r;
            w->d = ULLONG_MAX;
            pqueue_insert(pqueue, w);
          }
          
          // set new lower bound
          p_vw += abs(ARR_LOOKUP(nw[0]->match_scores,
                                 nw[0]->score_width,
                                 w->p.i, w->p.j));
          p_vw += abs(ARR_LOOKUP(nw[1]->match_scores,
                                 nw[1]->score_width,
                                 w->p.i, w->p.k));
          p_vw += abs(ARR_LOOKUP(nw[2]->match_scores,
                                 nw[2]->score_width,
                                 w->p.j, w->p.k));
          
#ifdef DEBUG
          printf(" - %d, ", abs(ARR_LOOKUP(nw[0]->match_scores,
                                           nw[0]->score_width,
                                           w->p.i, w->p.j)));
          printf("%d, ",    abs(ARR_LOOKUP(nw[1]->match_scores,
                                           nw[1]->score_width,
                                           w->p.i, w->p.k)));
          printf("%d",      abs(ARR_LOOKUP(nw[2]->match_scores,
                                           nw[2]->score_width,
                                           w->p.j, w->p.k)));
#endif
          
          p_w = (v->d <= ULLONG_MAX - p_vw) ?
                v->d + p_vw :
                w->d;
          
          if (p_w < w->d) {
            // update w distance from s
#ifdef DEBUG
            printf(" - DECREASE -> %lld", p_w);
#endif
            
            w->d = p_w;
            pqueue_change_priority(pqueue, w->d, w);
          }
          
#ifdef DEBUG
          printf("\n");
          // pqueue_print(pqueue, stdout, print_node);
          // printf("\n");
#endif
        }
      }
    }
  }
  
  // result
  
  if (v != NULL && point_eq(v->p, target.p)) {
    printf("END (%d, %d, %d), D: %lld\n", v->p.i, v->p.j, v->p.k, v->d);
  }
  else {
    printf("No alignment found with Upper bound: %u\n", upper_bound);
  }
  
  // free memory
  
  for (int i = 0; i < 3; ++i) {
    needleman_wunsch_free(nw[i]);
    alignment_free(result[i]);
  }
  
  pqueue_free(pqueue);
  free(ns);
}
