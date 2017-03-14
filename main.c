#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "msa.h"

int main(int argc, char* argv[]) {
  if (argc != 5) {
    printf("usage: ./msa <seq1> <seq2> <seq3> <upper_bound>\n");
    exit(EXIT_FAILURE);
  }
  
  multi_align(argv[1], argv[2], argv[3], strtoul(argv[4], NULL, 10));
  
  exit(EXIT_SUCCESS);
}
