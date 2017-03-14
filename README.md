# Carrillo Lipman

> Algorithm explained in section 14.6.1 "A speedup for the exact solution" from Dan Gusfield's "Algorithms on Strings, Trees, and Sequences"  

Implementation of the Carrillo-Lipman algorithm for aligning three sequences.  
This is an heuristic of the multiple sequence alignment problem with the Sum-of-Pairs objective function.  
The algorithm takes in input three strings and a parameter *upper_bound*, and prints the score of the best possible alignment between the strings.  
&nbsp;

## Interface

```c
void multi_align(char *seq_a, char *seq_b, char *seq_c, size_t upper_bound);
```
&nbsp;

## Example

```c
size_t upper_bound = 42;

multi_align("seq1", "seq2", "seq3", upper_bound);
```
