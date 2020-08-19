#if !defined SEQTOOLS_H
#define SEQTOOLS_H

#include <stdlib.h>
#include <stdio.h>

#define MAXNAME 1000
#define MAXFILENAME 1000

struct sequence {

	char *name; /*Stores sequence name*/
	int length; /*Length of sequence*/
	int *seq;   /*Data*/
	int group;   /*Allows sequences to be assigned to different groups*/
};

struct data {

	int nseq;
	int maxlseq;
	int type; /*0=nucleotide, 1=amino acid*/
	int aligned; /*0=not, 1=aligned*/
	struct sequence **seqs;

};


/*Details from each pairwise alignment*/

struct palign {

	double *global_val;  /*Score for best alignment*/
	int **global_trace;   /*Traceback for best*/
	double *col_max;     /*Column maxima*/
	int *col_pos;        /*Positions of column maxima*/

	int length;			 /*Length of alignment*/
	int s1;				 /*Number of first sequence*/
	int s2;				 /*Number of second sequence*/
};


/*Details from a kwise alignment*/

struct kalign {

	int target;
	double *score;
	int **trace;
	int length;

};


struct pars {	/*keeps program parameters*/

	double rho; /*Rho: Recombination parameter*/
	double del; /*delta: Insert/delete parameter*/
	double eps; /*Epsilon: rate of gap extension*/
	double term;	/*Termination rate*/
	double sizeL;/*Total length of (non-target) sequences*/
	double piM;	/*Prob of starting in match*/
	double piI;	/*Prob of starting in insert*/
	double mm;	/*Rate of match to match*/
	double gm;	/*Rate of gap to match*/
	double pmatch; /*Alternative way of specifying match proabilities*/


	double lrho;
	double ldel;
	double leps;
	double lterm;
	double lsizeL;
	double lpiM;
	double lpiI;
	double lmm;
	double lgm;
	double dm;
	double ldm;

	double **lsm;  /*Matrix for scoring matches (ln probs)*/
	double *lsi;   /*Vector for scoring insert states (ln probs)*/
	double **sm;   /*Non-logged version*/
	double *si;    /*Likewise*/

	/*matrices for baum-welch algorithm*/
	double **sum_sm;  /*Matrix for scoring matches (ln probs)*/
	double *sum_si;   /*Vector for scoring insert states (ln probs)*/
	double **new_sm;   /*Non-logged version*/
	double *new_si;    /*Likewise*/

	/*Other params for baum-welch algorithm*/
	double match_to_match;
	double match_to_insert;
	double match_to_delete;
	double insert_to_match;
	double insert_to_insert;
	double delete_to_delete;
	double delete_to_match;
	double new_del;
	double new_eps;
	double composite_llk;

	int type;	   /*0=NT, 1=AA*/
	int nstate;	   /*5 if NT, 21 if AA: NB includes nulls*/
	int ml;        /*Flag to indicate whether to just find ML path: default is yes (1)*/
	int estimate;  /*Flag to indicate whether to estimate parameters by EM: default is no (0)*/

	char *input_filename;
	char *output_tag;
	char *alignment_file;  /*Contains max accuracy alignment*/
	char *posterior_file;  /*Contains sums of posteriors for each target vs copy set*/
	char *baumwelch_file;  /*Contains estimated parameters*/

	int verbose; /*Defines detail in output*/


	/*Flag indicates whether sequences are to be grouped.  Also define target groups*/
	int ngroups;
	char **group_identifiers;
	char *target_name;
	int target_group;

	/*Flag indicates whether to use emission probabilities or dummy to '1'*/
	int emiss;

};

struct data * read_fasta(struct pars *my_pars, int aligned);
double watterson(int n);
char num2nuc(int i, int type);
void print_sequences(struct data *my_data, FILE *ofp);
void print_pair_align(struct palign *my_align, int *s1, int *s2);
void print_kalign(struct kalign *kwise, struct data *my_data, FILE *ofp);

#endif


