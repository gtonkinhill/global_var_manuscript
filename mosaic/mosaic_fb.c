
/*********************************************
 Program to detect mosaic gene structures
*********************************************/

#include "tools.h"
#include "seqtools.h"
#include "mosaic_fb.h"

#define DEBUG 0


main(int argc, char *argv[]) {

	int target;
	long seed = -setseed();
	struct data *my_data;
	struct pars *my_pars;
	struct matrices *my_matrices;
	clock_t start, end;
	double prev_llk, sum_prob;
	int i,j;
	int iteration_count=0;
	FILE *ofp;

	printf("\n\n*** Running MOSAIC ***\n\nA program for detecting mosaic gene structures");
	printf("\n\nProgram returns (for each sequence) a maximum accuracy alignment path\n\n");

	start = clock();

	my_pars = (struct pars *) get_pars(argc, argv);
	my_data = read_fasta(my_pars, 0);

	printf("\n\nRead %i sequences\n\n", my_data->nseq);
	if (DEBUG) print_sequences(my_data, stdout);

	/*If groups identified, need to classify sequences*/
	if (my_pars->ngroups>1) {
		classify_sequences(my_pars, my_data);
	}

	/*Add null term to end of each sequence - makes for simpler calculations*/
	my_data = add_null_term(my_data);

	printf("\n\n*** Computing maximum accuracy paths for each sequence ***\n\n");

	/*Allocate memory*/
	my_matrices = (struct matrices *) allocate_matrices(my_data, my_pars);

	/*Perform baum welch algorithm if required*/
	if (my_pars->estimate==1){
		prev_llk = VERYSMALL;
		while ((1-my_pars->composite_llk/prev_llk) > 0.01){
			iteration_count+=1;
			prev_llk = my_pars->composite_llk;
			my_pars->composite_llk=0.0;
			my_pars->match_to_insert=0.0;
			my_pars->match_to_match=0.0;
			my_pars->match_to_delete=0.0;
			my_pars->insert_to_insert=0.0;
			my_pars->insert_to_match=0.0;
			my_pars->delete_to_delete=0.0;
			my_pars->delete_to_match=0.0;

			printf("composite llk: %f\n", my_pars->composite_llk);
			for (target=1; target<=my_data->nseq; target++) {
				if (my_data->seqs[target]->group == my_pars->target_group) {
					kalign_fb(my_data, my_pars, my_matrices, target);
				}
			}
			printf("composite llk: %f\n", my_pars->composite_llk);

			/*update matrices*/
			for (i=0;i<my_pars->nstate;i++) {
				my_pars->si[i] = my_pars->new_si[i];
				my_pars->lsi[i] = (double) log(my_pars->new_si[i]);
			}

			for (i=0;i<my_pars->nstate;i++) for (j=0;j<=my_pars->nstate;j++) {
				my_pars->sm[i][j] = (double) my_pars->new_sm[i][j];
				my_pars->lsm[i][j] = (double) log(my_pars->new_sm[i][j]);
			}
			my_pars->del = my_pars->new_del;
			my_pars->eps = my_pars->new_eps;
			my_pars->ldel = (double) log(my_pars->new_del);
			my_pars->leps = (double) log(my_pars->new_eps);

			/*Make minimum probability*/
			/*add pseudo probs for insert emissions*/
			sum_prob = 0.0;
			for (i=0;i<my_pars->nstate;i++) {
				my_pars->si[i] = my_pars->si[i];
				sum_prob += my_pars->si[i];
			}
			for (i=0;i<my_pars->nstate;i++) {
				my_pars->si[i] = my_pars->si[i]/sum_prob;
				my_pars->lsi[i] = (double) log(my_pars->si[i]);
			}

			/*add pseudo probs for emissions*/
			for (i=0;i<my_pars->nstate;i++) {
				sum_prob = 0.0;
				for (j=0;j<my_pars->nstate;j++) {
					my_pars->sm[i][j] = my_pars->sm[i][j];
					sum_prob += my_pars->sm[i][j];
				}
				for (j=0;j<my_pars->nstate;j++) {
					my_pars->sm[i][j] = my_pars->sm[i][j]/sum_prob;
					my_pars->lsm[i][j] = (double) log(my_pars->sm[i][j]);
				}
			}
		}


		/*Write to output file*/
		ofp = fopen(my_pars->baumwelch_file, "w");
		fprintf(ofp,"#File containing estimated parameters from baumwelch algorithm\n");
		fprintf(ofp,"\n\nConverged in %i iterations.\n", iteration_count);
		fprintf(ofp,"\n\nInsertion probability (del): %f\n", my_pars->del);
		fprintf(ofp,"\nExtension probability (eps): %f\n", my_pars->eps);

		fprintf(ofp,"Insertion emission probs\n");
		for (i=0;i<my_pars->nstate;i++) {
			fprintf(ofp,"%f,", my_pars->si[i]);
		}
		fprintf(ofp,"\n\n");

		fprintf(ofp,"Emission probs\n");
		for (i=0;i<my_pars->nstate;i++){
			for (j=0;j<my_pars->nstate;j++) {
				fprintf(ofp,"%f,", my_pars->sm[i][j]);
			}
			fprintf(ofp,"\n");
		}
		fprintf(ofp,"\n\n");
		fflush(ofp);

		fclose(ofp);
	}

	else {
		/*If target group defined, run each target sequence against rest*/
		for (target=1; target<=my_data->nseq; target++) {
			if (my_data->seqs[target]->group == my_pars->target_group) {
				if (my_pars->ml) kalign_vt(my_data, my_pars, my_matrices, target);
				else if (!my_pars->estimate) kalign_fb(my_data, my_pars, my_matrices, target);
				else {
					printf("\n\n*** Warning: EM estimation not yet implemented ***\n\n");
					exit(0);
				}
			}
		}
	}

	/*Deallocate memory for matrices*/
	deallocate_matrices(my_data, my_pars, my_matrices);

	end = clock();

	printf("\n\n*** Program completed in %.3lf secs CPU time! ***\n\n", ((double) (end-start)/CLOCKS_PER_SEC));

	exit(0);

}


/*Have to add a null character to end of sequences for following algorithm*/

struct data * add_null_term(struct data *my_data) {

	int seq;

	for (seq=1;seq<=my_data->nseq;seq++) {
		my_data->seqs[seq]->seq = (int *) realloc(my_data->seqs[seq]->seq, (size_t) (my_data->seqs[seq]->length+2)*sizeof(int));
		my_data->seqs[seq]->seq[my_data->seqs[seq]->length+1]=0;
	}

	return my_data;

}



/************************************
Get parameters from command line
************************************/

struct pars * get_pars(int argc, char *argv[]) {

	int i, j, k, l;
	FILE *ifp, *ofp;
	char *in_str;
	struct pars *my_pars;

	extern double emiss_gap_nt[NSTATE_NT];
	extern double emiss_match_nt[NSTATE_NT][NSTATE_NT];
	extern double emiss_gap_aa[NSTATE_AA];
	extern double emiss_match_aa[NSTATE_AA][NSTATE_AA];

	time_t calendar_time;
	struct tm *output_time;

	my_pars = (struct pars *) malloc((size_t) sizeof(struct pars));

	printf("\n\n*** Reading input parameters and initialising ***\n\n");

	/*Initialise with default values*/
	my_pars->del=0.025;
	my_pars->eps=0.1;
	my_pars->rho=0.001;
	my_pars->term=0.01;
	my_pars->pmatch = 0.8;
	my_pars->type = -1;

	my_pars->input_filename = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
	my_pars->input_filename[0]='\0';
	my_pars->output_tag = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
	my_pars->output_tag[0] = '\0';
	my_pars->posterior_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->posterior_file[0]='\0';
	my_pars->alignment_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->alignment_file[0]='\0';
	my_pars->baumwelch_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->baumwelch_file[0]='\0';

	my_pars->ngroups=1;
	my_pars->target_group=1;

	my_pars->verbose=1;
	my_pars->emiss=1;
	my_pars->ml=1;
	my_pars->estimate=0;

	my_pars->match_to_insert=0.0;
	my_pars->match_to_match=0.0;
	my_pars->match_to_delete=0.0;
	my_pars->insert_to_insert=0.0;
	my_pars->insert_to_match=0.0;
	my_pars->delete_to_delete=0.0;
	my_pars->delete_to_match=0.0;
	my_pars->composite_llk=0.0;



	for (k=1; k<argc; k++) if (strchr(argv[k], '-')) {

		in_str = argv[k];

		/*Indicates that data are AA*/
		if (strcmp(in_str, "-aa") == 0) {
			my_pars->type = 1;
			my_pars->nstate = NSTATE_AA;
			my_pars->sm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			my_pars->lsm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			my_pars->new_sm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			my_pars->new_si = dvector(0, NSTATE_AA);
			my_pars->sum_sm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			my_pars->sum_si = dvector(0, NSTATE_AA);

			for (i=0;i<NSTATE_AA;i++) for (j=0;j<NSTATE_AA;j++) {
				my_pars->lsm[i][j]=emiss_match_aa[i+1][j+1];
				my_pars->sm[i][j] = exp(emiss_match_aa[i+1][j+1]);
				my_pars->new_sm[i][j] = (double) 0.0;
				my_pars->sum_sm[i][j] = (double) 0.0;
				// printf("new_sm=%f\n", my_pars->new_sm[i][j]);
			}
			my_pars->si = dvector(0, NSTATE_AA);
			my_pars->lsi = dvector(0, NSTATE_AA);
			for (i=0;i<NSTATE_AA;i++) {
				my_pars->lsi[i] = emiss_gap_aa[i];
				my_pars->si[i] = exp(emiss_gap_aa[i]);
				my_pars->new_si[i] = (double) 0.0;
				my_pars->sum_si[i] = (double) 0.0;
				// printf("new_si=%f\n", my_pars->new_si[i]);
			}
			my_pars->nstate = NSTATE_AA;
		}

		/*Indicates that data are NTs*/
		if (strcmp(in_str, "-nt") == 0) {
			my_pars->type = 0;
			my_pars->nstate = NSTATE_NT;
			my_pars->sm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			my_pars->lsm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			my_pars->new_sm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			my_pars->new_si = dvector(0, NSTATE_NT);
			my_pars->sum_sm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			my_pars->sum_si = dvector(0, NSTATE_NT);

			for (i=0;i<NSTATE_NT;i++) for (j=0;j<NSTATE_NT;j++) {
				my_pars->lsm[i][j]=emiss_match_nt[i+1][j+1];
				my_pars->sm[i][j] = exp(emiss_match_nt[i+1][j+1]);
			}
			my_pars->si = dvector(0, NSTATE_NT);
			my_pars->lsi = dvector(0, NSTATE_NT);
			for (i=0;i<NSTATE_NT;i++) {
				my_pars->lsi[i]=emiss_gap_nt[i];
				my_pars->si[i] = exp(emiss_gap_nt[i]);
			}
			my_pars->nstate = NSTATE_NT;
		}

		/*Input filename*/
		if (strcmp(in_str, "-seq") == 0)
			strncpy(my_pars->input_filename, argv[k+1], MAXFILENAME);

		/*Output filename tag*/
		if (strcmp(in_str, "-tag") == 0)
			strncpy(my_pars->output_tag, argv[k+1], MAXFILENAME);

		/*Classify sequences into groups*/
		if (strcmp(in_str, "-group") == 0) {
			my_pars->ngroups = atoi(argv[k+1]);
			my_pars->group_identifiers = (char **) malloc((size_t) (my_pars->ngroups+1)*sizeof(char));
			for (l=1;l<=my_pars->ngroups; l++) {
				my_pars->group_identifiers[l] = (char *) malloc((size_t) MAXNAME*sizeof(char));
				strncpy(my_pars->group_identifiers[l], argv[k+l+1], MAXNAME);
			}
			printf("\nGroups identified: ");
			for (l=1;l<=my_pars->ngroups;l++) printf("%s\t", my_pars->group_identifiers[l]);
		}

		/*Identify a group of sequences for target*/
		if (strcmp(in_str, "-target") == 0) {
			my_pars->target_name = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
			strncpy(my_pars->target_name, argv[k+1], MAXFILENAME);
			printf("\nTarget group: %s", my_pars->target_name);
			my_pars->target_group=0;
		}

		if (strcmp(in_str, "-rec") == 0) {
			my_pars->rho = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-del") == 0) {
			my_pars->del = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-eps") == 0) {
			my_pars->eps = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-term") == 0) {
			my_pars->term = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-e1") == 0) {
			my_pars->emiss=0;
		}

		if (strcmp(in_str, "-pmatch") == 0) {
			my_pars->pmatch = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-ma") == 0) {
			my_pars->ml = 0;
			printf("\nReturning maximum accuracy alignment\n");
		}

		if (strcmp(in_str, "-estimate") == 0) {
			my_pars->estimate = 1;
			printf("\nReturning baumwelch parameter estimation\n");
		}


	}


	if (my_pars->type<0) {
		printf("\n\n*** Error: Need to indicate whether data is protein (-aa) or nucleotide (-nt) ***\n\n");
		exit(0);
	}

	/*Check that all parameters are >0 and sum is <1*/
	if (my_pars->rho<=0) my_pars->rho = 1e-32;
	if (my_pars->del<=0) my_pars->del = 1e-6;
	if (my_pars->eps<=0) my_pars->eps = 1e-6;
	if (my_pars->term<=0) my_pars->term = 1e-6;

	if (2*my_pars->del+my_pars->rho+my_pars->term >= 1) {
		printf("\n\n*** Error: event probabilities sum to >=1 ***\n\n");
		exit(1);
	}


	/*Stationary probabilities for match (M) and Insert (I)*/
	my_pars->piM = (double) 2/3;
	my_pars->piI = (double) 1-my_pars->piM;
	my_pars->lpiM = log(my_pars->piM);
	my_pars->lpiI = log(my_pars->piI);
	my_pars->mm = (double) 1-2*my_pars->del-my_pars->rho-my_pars->term;
	my_pars->gm = (double) 1-my_pars->eps-my_pars->rho-my_pars->term;
	my_pars->dm = 1-my_pars->eps;
	my_pars->ldm = log(my_pars->dm);

	/*Make log versions*/
	my_pars->ldel=log(my_pars->del);
	my_pars->leps=log(my_pars->eps);
	my_pars->lrho=log(my_pars->rho);
	my_pars->lterm=log(my_pars->term);
	my_pars->lmm=log(my_pars->mm);
	my_pars->lgm=log(my_pars->gm);

	printf("\n\nModel parameters\n");
	printf("Rec\t%.4lf\n", my_pars->rho);
	printf("Del\t%.4lf\n", my_pars->del);
	printf("Eps\t%.4lf\n", my_pars->eps);
	printf("Term\t%.4lf\n", my_pars->term);
	printf("p(Match_state)\t%.4lf", my_pars->piM);
	printf("\n");
	printf("Esimtate params = %i\n", my_pars->estimate);


	/*Check input filename*/
	ifp = fopen(my_pars->input_filename, "r");
	if (!ifp) {
		printf("\n\n***Error: Cannot file input file (%s) ***\n\n", my_pars->input_filename);
		exit(1);
	}
	else {
		fclose(ifp);
	}

	/*Make output filenames*/
	if (my_pars->output_tag[0] != '\0') {
		my_pars->posterior_file = strncpy(my_pars->posterior_file, my_pars->output_tag, MAXFILENAME);
		my_pars->posterior_file = strcat(my_pars->posterior_file, "_post.txt");
		my_pars->alignment_file = strncpy(my_pars->alignment_file, my_pars->output_tag, MAXFILENAME);
		my_pars->alignment_file = strcat(my_pars->alignment_file, "_align.txt");
		my_pars->baumwelch_file = strncpy(my_pars->baumwelch_file, my_pars->output_tag, MAXFILENAME);
		my_pars->baumwelch_file = strcat(my_pars->baumwelch_file, "_baumWelch.txt");

	}
	else {
		my_pars->posterior_file = strcat(my_pars->posterior_file,  "post.txt");
		my_pars->alignment_file = strcat(my_pars->alignment_file,  "align.txt");
		my_pars->baumwelch_file = strcat(my_pars->baumwelch_file,  "baumWelch.txt");
	}

	/*Check things about groups and targets*/
	/*NB with a single target, program runs as originally*/

	if (my_pars->target_group == 0) {

		for (k=1;k<=my_pars->ngroups;k++) {
			if (strcmp(my_pars->target_name, my_pars->group_identifiers[k])==0) {
				my_pars->target_group = k;

			}
		}
		if (my_pars->target_group==0) {
			printf("\n\n*** Error: Could not match target name with group labels ***\n\n");
			exit(1);
		}
		else {
			printf("\nTarget identified as group %i", my_pars->target_group);
		}
	}

	/*If match probs need to re-organise*/
	if (my_pars->estimate==1){
		if (my_pars->pmatch>0) {
			for (i=0;i<my_pars->nstate;i++) {
				my_pars->si[i] = (double) 1/my_pars->nstate;
				my_pars->lsi[i] = (double) log(my_pars->si[i]);
			}
			for (i=0;i<my_pars->nstate;i++) for (j=0;j<=my_pars->nstate;j++) {
				if (i==j) {
					my_pars->sm[i][j] = (double) my_pars->pmatch;
					my_pars->lsm[i][j] = (double) log(my_pars->pmatch);
				}
				else {
					my_pars->sm[i][j] = (double) (1-my_pars->pmatch)/(my_pars->nstate-1);
					my_pars->lsm[i][j] = (double) log((1-my_pars->pmatch)/(my_pars->nstate-1));
				}
			}
		}

		/*If cut emission probs need to re-organise*/
		if (my_pars->emiss==0) {
			for (i=0;i<my_pars->nstate;i++) {
				my_pars->si[i] = 1.0;
				my_pars->lsi[i] = 0.0;
				for (j=0;j<=my_pars->nstate;j++) {
					my_pars->sm[i][j]=1.0;
					my_pars->lsm[i][j]=0.0;
				}
			}
		}
	}

	/*Initialise some output files*/
	ofp = fopen(my_pars->posterior_file, "w");
	fprintf(ofp,"#File containing summed posteriors for each target vs copy set\n");
	fprintf(ofp,"#Input parameters = ");
	for (i=0;i<argc;i++) fprintf(ofp,"%s ", argv[i]);
	time(&calendar_time);
	output_time = localtime(&calendar_time);
	fprintf(ofp,"\n#Created on %s\n", asctime(output_time));
	fprintf(ofp, "Sequence\tLength\tLogLk\tCopy_set\n");
	fclose(ofp);

	ofp = fopen(my_pars->alignment_file, "w");
	if (my_pars->ml) fprintf(ofp,"#File containing maximum likelihood alignment for each sequences\n");
	else fprintf(ofp,"#File containing maximum accuracy alignment for each sequences\n");
	fprintf(ofp,"#Input parameters = ");
	for (i=0;i<argc;i++) fprintf(ofp,"%s ", argv[i]);
	fprintf(ofp,"\n#Created on %s\n", asctime(output_time));
	fclose(ofp);

	return my_pars;
}



/*Do kwise alignment with forward-backward*/

void kalign_fb(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	extern double char_frequencies[NSTATE_AA];

	int pos_target, pos_seq, seq, l1, l2, *s1, *s2, i, j,  tmp_copy;
	double ***forward_m, ***backward_m, ***forward_i, ***backward_i, ***forward_d, ***backward_d;
	double llk_f, llk_b=SMALL, max_r, max_rn, llk_r, prob_loc_seq, sum_pos;
	FILE *ofp;
	clock_t start, end;
	double total_sum_si, total_sum_sm;


	printf("\rAligning sequence %5i to rest using Max Acc", target);
	start = clock();

	tmp_copy = my_matrices->who_copy[target];
	my_matrices->who_copy[target]=0;

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;
	my_pars->sizeL=0.0;
	for (seq=1;seq<=my_data->nseq;seq++) if (seq!=target) my_pars->sizeL+=(double) my_data->seqs[seq]->length;
	my_pars->lsizeL = (double) log(my_pars->sizeL);

	if (DEBUG) printf("\nSizeL = %.0lf",my_pars->sizeL);

	/*Make foward and backward matrices*/
	forward_m = (double ***) my_matrices->m1_m;
	forward_i = (double ***) my_matrices->m1_i;
	forward_d = (double ***) my_matrices->m1_d;
	backward_m = (double ***) my_matrices->m2_m;
	backward_i = (double ***) my_matrices->m2_i;
	backward_d = (double ***) my_matrices->m2_d;

	/*Set everything to small*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2=my_data->seqs[seq]->length;
		for (pos_target=0;pos_target<=(l1+1);pos_target++)
			for (pos_seq=0;pos_seq<=(l2+1); pos_seq++) {
				forward_m[seq][pos_seq][pos_target]=SMALL;
				forward_i[seq][pos_seq][pos_target]=SMALL;
				forward_d[seq][pos_seq][pos_target]=SMALL;
				backward_m[seq][pos_seq][pos_target]=SMALL;
				backward_i[seq][pos_seq][pos_target]=SMALL;
				backward_d[seq][pos_seq][pos_target]=SMALL;
			}
	}

	/*Initialise forward and backward matrices*/

	for (seq=1, max_r = SMALL;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_seq=1;pos_seq<=l2;pos_seq++){

			/*Forward Ms*/
			forward_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
			forward_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
			forward_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
			forward_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];
			if (pos_seq>1) {
				forward_d[seq][pos_seq][1] = (double) forward_m[seq][pos_seq-1][1] + log(my_pars->del + my_pars->eps*exp(forward_d[seq][pos_seq-1][1]-forward_m[seq][pos_seq-1][1]));
			}

			/*Find biggest value prior to recombination*/
			if (forward_m[seq][pos_seq][1] > max_r) max_r = forward_m[seq][pos_seq][1];
			if (forward_i[seq][pos_seq][1] > max_r) max_r = forward_i[seq][pos_seq][1];
			if (forward_d[seq][pos_seq][1] > max_r) max_r = forward_d[seq][pos_seq][1];

			/*Backward Ms*/
			backward_m[seq][pos_seq][l1] = my_pars->lterm;
			backward_i[seq][pos_seq][l1] = my_pars->lterm;
		}
	}



	/*Now loop forward over positions*/

	for (pos_target=2;pos_target<=l1;pos_target++) {

		/*First get contribution from recombination*/

		for (seq=1, llk_r=0.0, max_rn=SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				llk_r += (double) exp(forward_m[seq][pos_seq][pos_target-1]-max_r);
				llk_r += (double) exp(forward_i[seq][pos_seq][pos_target-1]-max_r);

			}
		}


		/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match*/
				forward_m[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->mm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->gm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_d[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->dm;
				forward_m[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho*my_pars->piM/my_pars->sizeL;
				forward_m[seq][pos_seq][pos_target] = (double) log(forward_m[seq][pos_seq][pos_target])+max_r;
				forward_m[seq][pos_seq][pos_target] += my_pars->lsm[s1[pos_target]][s2[pos_seq]];

				/*Insert*/
				forward_i[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq][pos_target-1]-max_r)*my_pars->del;
				forward_i[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq][pos_target-1]-max_r)*my_pars->eps;
				forward_i[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho*my_pars->piI/my_pars->sizeL;
				forward_i[seq][pos_seq][pos_target] = (double) log(forward_i[seq][pos_seq][pos_target])+max_r;
				forward_i[seq][pos_seq][pos_target] += my_pars->lsi[s1[pos_target]];

				/*Delete: NB cannot terminate from a delete state*/
				if (pos_target<l1 && pos_seq>1) {
					forward_d[seq][pos_seq][pos_target] =  exp(forward_d[seq][pos_seq-1][pos_target]-max_r)*my_pars->eps;
					forward_d[seq][pos_seq][pos_target] += exp(forward_m[seq][pos_seq-1][pos_target]-max_r)*my_pars->del;
					forward_d[seq][pos_seq][pos_target] = (double) log(forward_d[seq][pos_seq][pos_target]) + max_r;
				}

				/*Get new max_r*/
				if (forward_m[seq][pos_seq][pos_target] > max_rn) max_rn = forward_m[seq][pos_seq][pos_target];
				if (forward_i[seq][pos_seq][pos_target] > max_rn) max_rn = forward_i[seq][pos_seq][pos_target];
			}
		}

		/*Completed position in target sequence*/
		max_r = max_rn;
	}
	/*Completed forward matrices*/
	if (DEBUG>1) print_forward_matrices(my_data, my_pars, my_matrices, target, stdout);


	/*Now calculate likelihood*/

	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){
				llk_f += (double) exp(forward_m[seq][pos_seq][l1]-max_r);
				llk_f += (double) exp(forward_i[seq][pos_seq][l1]-max_r);
			}
	}
	llk_f = max_r + log(llk_f) + my_pars->lterm;
	printf("\nLog likelihood from forward algorithm  = %.5lf\n", llk_f);
	my_matrices->llk = llk_f;
	my_pars->composite_llk = my_pars->composite_llk + llk_f;



	/*Now loop backward*/
	for (pos_target = l1-1, max_r = my_pars->lterm; pos_target>0; pos_target--) {

		/*First get contribution from recombination*/

		for (seq=1, llk_r=0.0, max_rn=SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

  				llk_r += (double) exp(backward_m[seq][pos_seq][pos_target+1]-max_r)*(my_pars->sm[s1[pos_target+1]][s2[pos_seq]])*(my_pars->piM);
				llk_r += (double) exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*(my_pars->si[s1[pos_target+1]])*(my_pars->piI);
			}
		}

			/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			/*Note: Have to do Delete state first and in reverse position order*/

			for (pos_seq=l2;pos_seq>0;pos_seq--) {

				/*Delete*/
				if (pos_seq<l2) {
					backward_d[seq][pos_seq][pos_target] =  exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->eps;
					backward_d[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]]*my_pars->dm;
					backward_d[seq][pos_seq][pos_target] = (double) log(backward_d[seq][pos_seq][pos_target]) + max_r;
				}

				/*Insert*/
				backward_i[seq][pos_seq][pos_target] =  exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->eps*my_pars->si[s1[pos_target+1]];
				backward_i[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->gm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_i[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho/my_pars->sizeL;
				backward_i[seq][pos_seq][pos_target] = (double) log(backward_i[seq][pos_seq][pos_target])+max_r;

				/*Match*/
				backward_m[seq][pos_seq][pos_target] =  exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->mm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->del*my_pars->si[s1[pos_target+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->del;
				backward_m[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho/my_pars->sizeL;
				backward_m[seq][pos_seq][pos_target] = (double) log(backward_m[seq][pos_seq][pos_target])+max_r;

				/*Get new max_r*/
				if (backward_m[seq][pos_seq][pos_target] > max_rn) max_rn = backward_m[seq][pos_seq][pos_target];
				if (backward_i[seq][pos_seq][pos_target] > max_rn) max_rn = backward_i[seq][pos_seq][pos_target];
			}
		}
		/*Completed position in target sequence*/
		max_r = max_rn;
	}
	/*Completed backward matrices*/

	/*Calculate backwards likelihood*/
	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){
				llk_f += (double) exp(backward_m[seq][pos_seq][1]-max_r)*my_pars->sm[s1[1]][s2[pos_seq]]*my_pars->piM;
				llk_f += (double) exp(backward_i[seq][pos_seq][1]-max_r)*my_pars->si[s1[1]]*my_pars->piI;
			}
	}
	llk_f = max_r + log(llk_f) - log(my_pars->sizeL);
	printf("Log likelihood from backward algorithm = %.5lf\n\n", llk_f);



	if (DEBUG>1) print_backward_matrices(my_data, my_pars, my_matrices, target, stdout);


	if (my_pars->estimate==1){

		/*first emission*/
		printf("\n Estimating emission parameters (getting expectation)...\n");
		for (pos_target=1;pos_target<=l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					/* emission probability */
					// printf("pos_seq:%i pos_target:%i forward_m[seq][pos_seq][pos_target]: %f emit:%f\n", pos_seq, pos_target,forward_m[seq][pos_seq][pos_target],  exp(forward_m[seq][pos_seq][pos_target]+backward_m[seq][pos_seq][pos_target]-llk_f));
					my_pars->sum_sm[s2[pos_seq]][s1[pos_target]] += (double) exp(forward_m[seq][pos_seq][pos_target]+backward_m[seq][pos_seq][pos_target]-llk_f);
				}
			}
		}

		printf("\n Estimating emission parameters (calcualting new emission probs...\n");
		for (i=0;i<my_pars->nstate;i++) {
			total_sum_sm = 0.0;
			for (j=0;j<my_pars->nstate;j++) {
				total_sum_sm += my_pars->sum_sm[i][j];
			}
			printf("i: %i total_sm: %f\n", i, total_sum_sm);
			for (j=0;j<my_pars->nstate;j++) {
				/*sum over all sequences includinh using a dirichlet prior with parramter Aq_a for psudocounts*/
				my_pars->new_sm[i][j] = (my_pars->sum_sm[i][j] + my_pars->nstate*char_frequencies[i])/(total_sum_sm + my_pars->nstate);
			}
		}
		printf("match emission:\n");
		for (i=0;i<my_pars->nstate;i++) {
			for (j=0;j<my_pars->nstate;j++) {
				printf("%lf,", my_pars->new_sm[i][j]);
			}
		}
		printf("\n");

		/*insert emission*/
		for (pos_target=1;pos_target<=l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					/* state */
					my_pars->sum_si[s1[pos_target]] += exp(forward_i[seq][pos_seq][pos_target]+backward_i[seq][pos_seq][pos_target]-llk_f);
				}
			}
		}

		total_sum_si = 0.0;
		for (i=0;i<my_pars->nstate;i++) {
				total_sum_si += my_pars->sum_si[i];
		}
		for (i=0;i<my_pars->nstate;i++) {
				my_pars->new_si[i] = (my_pars->sum_si[i] + my_pars->nstate*char_frequencies[i])/(total_sum_si + my_pars->nstate);
		}
		printf("insert emission:\n");
		for (i=0;i<my_pars->nstate;i++) {
			printf("%lf,", my_pars->new_si[i]);
		}
		printf("\n");

		/*indel initiation*/

		/*match_to_match*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					my_pars->match_to_match += exp(forward_m[seq][pos_seq][pos_target] + my_pars->lmm + my_pars->lsm[s1[pos_target+1]][s2[pos_seq+1]] + backward_m[seq][pos_seq+1][pos_target+1] - llk_f);
				}
			}
		}

		/*match_to_insert*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					my_pars->match_to_insert += exp(forward_m[seq][pos_seq][pos_target] + my_pars->ldel + my_pars->lsi[s1[pos_target+1]] + backward_i[seq][pos_seq][pos_target+1] - llk_f);
				}
			}
		}

		/*match_to_delete*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=2;pos_seq<=l2;pos_seq++) {
					// printf("pos_seq:%i pos_target:%i forward_m[seq][pos_seq][pos_target]: %f forward_d[seq][pos_seq][pos_target+1]: %f emit:%f\n", pos_seq, pos_target,forward_m[seq][pos_seq][pos_target],backward_d[seq][pos_seq][pos_target+1],  exp(forward_m[seq][pos_seq][pos_target] + my_pars->ldel + backward_d[seq][pos_seq][pos_target+1] - llk_f));
					my_pars->match_to_delete += exp(forward_m[seq][pos_seq][pos_target] + my_pars->ldel + backward_d[seq][pos_seq+1][pos_target] - llk_f);
				}
			}
		}

		my_pars->new_del = (my_pars->match_to_insert + my_pars->match_to_delete)/(2*(my_pars->match_to_insert+my_pars->match_to_delete+my_pars->match_to_match));
		// my_pars->new_del = (my_pars->match_to_insert)/(my_pars->match_to_insert+my_pars->match_to_delete+my_pars->match_to_match);

		printf("match_to_insert:%f  match_to_delete:%f match_to_match:%f\n", my_pars->match_to_insert, my_pars->match_to_delete, my_pars->match_to_match);
		printf("\nNew del: %lf\n", my_pars->new_del);

		/*indel extension*/

		/*insert_to_match*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					my_pars->insert_to_match += exp(forward_i[seq][pos_seq][pos_target] + my_pars->lgm + my_pars->lsm[s1[pos_target+1]][s2[pos_seq+1]] + backward_i[seq][pos_seq+1][pos_target+1] - llk_f);
				}
			}
		}

		/*insert_to_insert*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					my_pars->insert_to_insert += exp(forward_i[seq][pos_seq][pos_target] + my_pars->leps + my_pars->lsi[s1[pos_target+1]] + backward_i[seq][pos_seq][pos_target+1] - llk_f);
				}
			}
		}

		/*delete_to_delete*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=2;pos_seq<=l2;pos_seq++) {
					// printf("pos_seq:%i pos_target:%i forward_m[seq][pos_seq][pos_target]: %f forward_d[seq][pos_seq][pos_target+1]: %f emit:%f\n", pos_seq, pos_target,forward_m[seq][pos_seq][pos_target],backward_d[seq][pos_seq][pos_target+1],  exp(forward_m[seq][pos_seq][pos_target] + my_pars->ldel + backward_d[seq][pos_seq][pos_target+1] - llk_f));
					my_pars->delete_to_delete += exp(forward_d[seq][pos_seq][pos_target] + my_pars->leps + backward_d[seq][pos_seq+1][pos_target] - llk_f);
				}
			}
		}

		/*delete_to_match*/
		for (pos_target=1;pos_target<l1;pos_target++) {
			for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
				l2 = my_data->seqs[seq]->length;
				s2 = my_data->seqs[seq]->seq;
				for (pos_seq=2;pos_seq<=l2;pos_seq++) {
					// printf("pos_seq:%i pos_target:%i forward_m[seq][pos_seq][pos_target]: %f forward_d[seq][pos_seq][pos_target+1]: %f emit:%f\n", pos_seq, pos_target,forward_m[seq][pos_seq][pos_target],backward_d[seq][pos_seq][pos_target+1],  exp(forward_m[seq][pos_seq][pos_target] + my_pars->ldel + backward_d[seq][pos_seq][pos_target+1] - llk_f));
					my_pars->delete_to_match += exp(forward_d[seq][pos_seq][pos_target] + my_pars->ldm + my_pars->lsm[s1[pos_target+1]][s2[pos_seq+1]] + backward_m[seq][pos_seq+1][pos_target+1] - llk_f);
				}
			}
		}


		my_pars->new_eps = (my_pars->insert_to_insert+my_pars->delete_to_delete)/(my_pars->insert_to_insert+my_pars->insert_to_match+my_pars->delete_to_delete+my_pars->delete_to_match);
		// my_pars->new_eps = my_pars->insert_to_insert/(my_pars->insert_to_insert+my_pars->insert_to_match);

		printf("\nNew eps: %lf\n", my_pars->new_eps);
	}


	/*Now calculate forward-backward probs*/

	for (seq=1;seq<=my_data->nseq;seq++) my_matrices->ppsum_match[seq]=0.0;
	for (pos_target=1;pos_target<=l1;pos_target++) {
		for (i=1;i<=3;i++) my_matrices->ppsum_state[pos_target][i]=0.0;
		for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2=my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++) {
				forward_m[seq][pos_seq][pos_target] = (double) exp(forward_m[seq][pos_seq][pos_target]+backward_m[seq][pos_seq][pos_target]-llk_f);
				forward_i[seq][pos_seq][pos_target] = (double) exp(forward_i[seq][pos_seq][pos_target]+backward_i[seq][pos_seq][pos_target]-llk_f);
				forward_d[seq][pos_seq][pos_target] = (double) exp(forward_d[seq][pos_seq][pos_target]+backward_d[seq][pos_seq][pos_target]-llk_f);
			}
		}
	}


	ofp = fopen(my_pars->posterior_file, "a");

	/*Now calculate the profile*/
	// printf("Target, Seq, Position target, Probability\n");
	// fprintf(ofp, "Target, Seq, Target length, Prior probabilities\n");
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		fprintf(ofp, "%s,%s,%i", my_data->seqs[target]->name, my_data->seqs[seq]->name, my_data->seqs[target]->length);
		sum_pos = 0.0;
		for (pos_target=1;pos_target<=l1;pos_target++) {
			l2=my_data->seqs[seq]->length;
			prob_loc_seq=0.0;
			for (pos_seq=1;pos_seq<=l2;pos_seq++) {
				/* Probability of aligning to seq at pos_target*/
				prob_loc_seq = prob_loc_seq + forward_m[seq][pos_seq][pos_target] + forward_i[seq][pos_seq][pos_target];
			}
			fprintf(ofp, ",%.5lf", prob_loc_seq);
			sum_pos += prob_loc_seq;
		}
		fprintf(ofp, "\n");

		// fprintf(ofp, "%s\t%i\t%f\t%s\n", my_data->seqs[target]->name
		// 	, my_data->seqs[target]->length
		// 	, sum_pos
		// 	, my_data->seqs[seq]->name);
		fflush(ofp);
	}

	fclose(ofp);




	if (DEBUG>1) print_posterior_matrices(my_data, my_pars, my_matrices, target);

	if (my_pars->verbose) print_max_acc_alignment(my_data, my_pars, my_matrices, target);


	/*Reset copy state*/
	my_matrices->who_copy[target]=tmp_copy;

	end = clock();

	printf("\n\n*** Aligned sequence in %.3lf secs CPU time! ***\n\n", ((double) (end-start)/CLOCKS_PER_SEC));

	return;


}


/*Assign sequences to groups*/

void classify_sequences(struct pars *my_pars, struct data *my_data) {

	int seq, j;

	for (seq=1;seq<=my_data->nseq;seq++) {
		my_data->seqs[seq]->group=0;
		for (j=1;j<=my_pars->ngroups;j++) {
			if (strstr(my_data->seqs[seq]->name, my_pars->group_identifiers[j])) {
				my_data->seqs[seq]->group=j;
			}
		}
		if (my_data->seqs[seq]->group<1) {
			printf("\n\n***Error: Could not assign sequence %s to group ***\n\n", my_data->seqs[seq]->name);
			exit(1);
		}
	}

}




struct matrices * allocate_matrices(struct data *my_data, struct pars *my_pars) {

	int seq, maxl, ct=0, tmp, i;
	struct matrices *my_matrices;

	my_matrices = (struct matrices *) malloc((size_t) sizeof(struct matrices));

	for (seq=2, maxl=my_data->seqs[1]->length;seq<=my_data->nseq; seq++) {
		if (my_data->seqs[seq]->length>maxl) maxl = my_data->seqs[seq]->length;
	}
	my_matrices->maxl = maxl;

	my_matrices->m1_m = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m1_i = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m1_d = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_m = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_i = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_d = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));

	/*
	for (seq=1;seq<=my_data->nseq;seq++) if (my_pars->ngroups==1 || my_data->seqs[seq]->group != my_pars->target_group) {
	*/
	for (seq=1;seq<=my_data->nseq;seq++) {
		my_matrices->m1_m[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m1_i[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m1_d[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_m[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_i[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_d[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		ct++;
	}

	/*Make list of who can be copied from*/
	my_matrices->who_copy = (int *) malloc((size_t) (my_data->nseq+1)*sizeof(int));
	for (seq=1;seq<=my_data->nseq;seq++) {
		my_matrices->who_copy[seq]=1;
		/*
		if (my_pars->ngroups==1 || my_data->seqs[seq]->group != my_pars->target_group)
			my_matrices->who_copy[seq]=1;
		else my_matrices->who_copy[seq]=0;
		*/
	}


	/*Make paths for best alignment*/
	my_matrices->maxpath_copy = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->maxpath_state = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->maxpath_pos = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->ppsum_match = (double *) malloc((size_t) (my_data->nseq+1)*sizeof(double));
	my_matrices->ppsum_state = (double **) dmatrix(0, maxl+1, 1, 3);

	tmp = (int) log10(maxl) + 1;
	my_matrices->tb_divisor = (double) pow(10, (double) tmp);

	printf("\n\n*** Allocated memory for matrices (%i) ***\n\n", ct);

	return my_matrices;
}



/*Deallocate memory for matrices*/

void deallocate_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices) {

	int i;

	for (i=1;i<=my_data->nseq; i++) {
		if (my_pars->ngroups==1 || my_data->seqs[i]->group != my_pars->target_group) {
			free_dmatrix(my_matrices->m1_m[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m1_i[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m1_d[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_m[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_i[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_d[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
		}
	}

	free(my_matrices->m1_m);
	free(my_matrices->m1_i);
	free(my_matrices->m1_d);
	free(my_matrices->m2_m);
	free(my_matrices->m2_i);
	free(my_matrices->m2_d);

	free(my_matrices->who_copy);

	free(my_matrices->maxpath_copy);
	free(my_matrices->maxpath_state);
	free(my_matrices->maxpath_pos);
	free(my_matrices->ppsum_match);
	free_dmatrix(my_matrices->ppsum_state, 0, my_matrices->maxl+1, 1, 3);

	free(my_matrices);

	printf("\n\n*** Memory deallocated ***\n\n");

	return;
}



/* To print max accuracy alginment*/


void print_max_acc_alignment(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int seq, pos_seq, pos_target, who_max, state_max, pos_max, l1, l2, cp, i, j;
	double max_acc, max_acc_rec, max_acc_rec_n;
	double ***post_m, ***post_i, ***post_d;
	FILE *ofp;

	/*Redirect*/
	post_m = (double ***) my_matrices->m1_m;
	post_i = (double ***) my_matrices->m1_i;
	post_d = (double ***) my_matrices->m1_d;

	/*First calculate maximum acc alignment*/

	l1 = my_data->seqs[target]->length;

	/*Initialise for delete state and get max_acc_rec for first position*/
	for (seq=1, max_acc_rec=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2 = my_data->seqs[seq]->length;

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {

			if (post_m[seq][pos_seq][1]>max_acc_rec) max_acc_rec = post_m[seq][pos_seq][1];
			if (post_i[seq][pos_seq][1]>max_acc_rec) max_acc_rec = post_i[seq][pos_seq][1];

			/*Delete state*/
			if (post_m[seq][pos_seq-1][1] > post_d[seq][pos_seq-1][1]) {
				post_d[seq][pos_seq][1] += post_m[seq][pos_seq-1][1];
			}
			else {
				post_d[seq][pos_seq][1] += post_d[seq][pos_seq-1][1];
			}
		}
	}

	for (pos_target=2;pos_target<=l1;pos_target++) {

		for (seq=1, max_acc_rec_n=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match state*/
				if (post_d[seq][pos_seq-1][pos_target-1]>max_acc_rec) {
					post_m[seq][pos_seq][pos_target] += post_d[seq][pos_seq-1][pos_target-1];
				}
				else {
					post_m[seq][pos_seq][pos_target] += max_acc_rec;
				}

				/*Insert state*/
				post_i[seq][pos_seq][pos_target] += max_acc_rec;


				/*Delete state*/
				if (post_m[seq][pos_seq-1][pos_target] > post_d[seq][pos_seq-1][pos_target]) {
					post_d[seq][pos_seq][pos_target] += post_m[seq][pos_seq-1][pos_target];
				}
				else {
					post_d[seq][pos_seq][pos_target] += post_d[seq][pos_seq-1][pos_target];
				}

				/*Find new maximum for next step*/
				if (post_m[seq][pos_seq][pos_target]>max_acc_rec_n) max_acc_rec_n = post_m[seq][pos_seq][pos_target];
				if (post_i[seq][pos_seq][pos_target]>max_acc_rec_n) max_acc_rec_n = post_i[seq][pos_seq][pos_target];


			}
		}
		max_acc_rec = max_acc_rec_n;
	}

	if (DEBUG>1) {
		printf("\n\nSummed posteriors\n\n");
		print_posterior_matrices(my_data, my_pars, my_matrices, target);
	}

	/*Now work backwards finding the best alignment*/

	/*First find best value in end state*/
	cp = 2*my_matrices->maxl;
	for (seq=1, max_acc=0.0; seq<=my_data->nseq; seq++) if (my_matrices->who_copy[seq]) {
		l2 = my_data->seqs[seq]->length;
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			if (post_m[seq][pos_seq][l1]>max_acc) {
				max_acc = post_m[seq][pos_seq][l1];
				who_max = seq;
				state_max = 1;
				pos_max = pos_seq;
			}
			if (post_i[seq][pos_seq][l1]>max_acc) {
				max_acc = post_i[seq][pos_seq][l1];
				who_max = seq;
				state_max = 2;
				pos_max = pos_seq;
			}
		}
	}
	my_matrices->maxpath_copy[cp] = who_max;
	my_matrices->maxpath_state[cp] = state_max;
	my_matrices->maxpath_pos[cp] = pos_max;

	if (DEBUG) printf("\nEnding in state %i for sequence %i at position %i (%.3lf)", state_max, who_max, pos_max, max_acc);

	pos_target=l1;
	while(pos_target>=1) {

		/*Currently in delete state - can only come from delete or match state*/
		if (my_matrices->maxpath_state[cp]==3) {
			max_acc = post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target];
			who_max = my_matrices->maxpath_copy[cp];
			state_max = 3;
			pos_max = my_matrices->maxpath_pos[cp]-1;
			if (post_m[who_max][my_matrices->maxpath_pos[cp]-1][pos_target]>=max_acc) {
				state_max = 1;
			}
		}

		/*If from I then just look at marginal best*/
		/*Suggests that allowing recombination from delete state would be a simplifying idea*/
		else if (my_matrices->maxpath_state[cp]==2) {

			for (seq=1, max_acc=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

				l2 = my_data->seqs[seq]->length;

				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					if (post_m[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_m[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 1;
						pos_max = pos_seq;
					}
					if (post_i[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_i[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 2;
						pos_max = pos_seq;
					}
				}
			}
		}

		/*If from M then  look at marginal best and delete state*/
		else if (my_matrices->maxpath_state[cp]==1) {

			for (seq=1, max_acc=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

				l2 = my_data->seqs[seq]->length;

				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					if (post_m[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_m[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 1;
						pos_max = pos_seq;
					}
					if (post_i[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_i[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 2;
						pos_max = pos_seq;
					}
				}
			}

			if (post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target-1] > max_acc) {
				max_acc = post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target-1];
				who_max = my_matrices->maxpath_copy[cp];
				state_max = 3;
				pos_max = my_matrices->maxpath_pos[cp]-1;
			}
		}

		cp--;
		if (cp<=0) {
			printf("\n\n***Error: reconstructed path longer than maximum possible ***\n\n");
			exit(1);
		}
		my_matrices->maxpath_copy[cp] = who_max;
		my_matrices->maxpath_state[cp] = state_max;
		my_matrices->maxpath_pos[cp] = pos_max;

		if (my_matrices->maxpath_state[cp+1] != 3) pos_target--;

		if (DEBUG) printf("\nNext = state %i for sequence %i at position %i (%.3lf).  Pos target = %i", state_max, who_max, pos_max, max_acc, pos_target);

	}

	if (DEBUG) {

		printf("\n\nMax Accuracy alignment\n\nWho\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_copy[i]);
		printf("\nState\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_state[i]);
		printf("\nPos\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_pos[i]);
		printf("\n\n");
	}


	ofp = fopen(my_pars->alignment_file, "a");
	fprintf(ofp,"\nTarget: %s\tLength: %i\tLlk: %.3lf\n",my_data->seqs[target]->name, my_data->seqs[target]->length, my_matrices->llk);

	/*First print target sequence*/
	cp++;
	fprintf(ofp,"%10s\t", my_data->seqs[target]->name);
	for (i=cp,pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==3) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[target]->seq[pos_target], my_data->type));
			pos_target++;
		}
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do matching*/
	fprintf(ofp,"          \t");
	for (i=cp, pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==1) {
			if (my_data->seqs[target]->seq[pos_target] == my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]]) fprintf(ofp,"|");
			else fprintf(ofp," ");
			pos_target++;
		}
		else if (my_matrices->maxpath_state[i]==2) {
			pos_target++;
			fprintf(ofp,"^");
		}
		else fprintf(ofp,"~");
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do copy tracks - switch whenever gets to new value*/
	fprintf(ofp,"%10s\t",my_data->seqs[my_matrices->maxpath_copy[cp]]->name);
	for (i=cp;i<=2*my_matrices->maxl;i++) {

		/*Check to see if need to make recombination event*/
		if (i>cp && my_matrices->maxpath_copy[i]!=my_matrices->maxpath_copy[i-1]) {
			fprintf(ofp,"\n%10s\t", my_data->seqs[my_matrices->maxpath_copy[i]]->name);
			for (j=1;j<=(i-cp);j++) fprintf(ofp," ");
		}

		if (my_matrices->maxpath_state[i]==2) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]], my_data->type));
		}
	}

	fprintf(ofp,"\n\n");
	fflush(ofp);
	fclose(ofp);


}





void print_posterior_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int seq, pos, pos2;
	double *site_sums;

	printf("\n\nPosterior matrices for sequence %s\n", my_data->seqs[target]->name);

	site_sums = (double *) malloc((size_t) (my_data->seqs[target]->length+1)*sizeof(double));
	for (pos=1;pos<=my_data->seqs[target]->length; pos++) site_sums[pos]=0.0;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		printf("\n\tSequence %s: Match\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_m[seq][pos2][pos]);
				site_sums[pos] += my_matrices->m1_m[seq][pos2][pos];
			}
			printf("\n");
		}
		printf("\n\tSequence %s: Insert\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_i[seq][pos2][pos]);
				site_sums[pos] += my_matrices->m1_i[seq][pos2][pos];
			}
			printf("\n");
		}
		printf("\n\tSequence %s: Delete\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_d[seq][pos2][pos]);
			}
			printf("\n");
		}
	}

	printf("\n\nSite specific posterior sums\n\nPosition\tSum\n");
	for (pos=1; pos<=my_data->seqs[target]->length;pos++) {
		printf("%i\t%.4lf\n", pos, site_sums[pos]);
	}

	printf("\n\n");


	free(site_sums);

	return;


}



/*Print out forward matrices*/

void print_forward_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target, FILE *ofp) {

	int seq, pos_target, pos_seq, l1, l2, *s1, *s2;


	fprintf(ofp, "\n\n*** Forward matrices for sequence: %s ***\n\n", my_data->seqs[target]->name);

	l1 = my_data->seqs[target]->length;
	s1 = my_data->seqs[target]->seq;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		fprintf(ofp,"\n\nCopy_target %s\n\n", my_data->seqs[seq]->name);

		l2 = my_data->seqs[seq]->length;
		s2 = my_data->seqs[seq]->seq;

		fprintf(ofp,"Match %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_m[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nInsert %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_i[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nDelete %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_d[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");

		}
	}

}



/*Print out forward matrices*/

void print_backward_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target, FILE *ofp) {

	int seq, pos_target, pos_seq, l1, l2, *s1, *s2;


	fprintf(ofp, "\n\n*** Backward matrices for sequence: %s ***\n\n", my_data->seqs[target]->name);

	l1 = my_data->seqs[target]->length;
	s1 = my_data->seqs[target]->seq;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		fprintf(ofp,"\n\nCopy_target %s\n\n", my_data->seqs[seq]->name);

		l2 = my_data->seqs[seq]->length;
		s2 = my_data->seqs[seq]->seq;

		fprintf(ofp,"Match %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_m[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nInsert %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_i[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nDelete %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++)
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_d[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");

		}
	}

}








/*Do kwise alignment with Viterbi*/

void kalign_vt(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {
	clock_t start_vt, end_vt;
	int pos_target, pos_seq, seq, l1, l2, *s1, *s2, tmp_copy;
	int who_max, state_max, pos_max, who_max_n, state_max_n, pos_max_n;
	int who_next, state_next, pos_next, cp, i, j;
	double ***vt_m, ***vt_i, ***vt_d, ***tb_m, ***tb_i, ***tb_d;
	double max_r, max_rn;
	char tmp_name[MAXNAME];
	FILE *ofp;

	printf("\rAligning sequence %5i to rest using ML", target);
	start_vt = clock();

	tmp_copy = my_matrices->who_copy[target];
	my_matrices->who_copy[target]=0;

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;
	my_pars->sizeL=0.0;
	for (seq=1;seq<=my_data->nseq;seq++) if (seq!=target) my_pars->sizeL+=(double) my_data->seqs[seq]->length;
	my_pars->lsizeL = (double) log(my_pars->sizeL);

	if (DEBUG) printf("\nSizeL = %.0lf",my_pars->sizeL);

	/*Make Viterbi and traceback matrices*/
	vt_m = (double ***) my_matrices->m1_m;
	vt_i = (double ***) my_matrices->m1_i;
	vt_d = (double ***) my_matrices->m1_d;
	tb_m = (double ***) my_matrices->m2_m;
	tb_i = (double ***) my_matrices->m2_i;
	tb_d = (double ***) my_matrices->m2_d;

	/*Set everything to small*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2=my_data->seqs[seq]->length;
		for (pos_target=0;pos_target<=(l1+1);pos_target++)
			for (pos_seq=0;pos_seq<=(l2+1); pos_seq++) {
				vt_m[seq][pos_seq][pos_target]=SMALL;
				vt_i[seq][pos_seq][pos_target]=SMALL;
				vt_d[seq][pos_seq][pos_target]=SMALL;
				tb_m[seq][pos_seq][pos_target]=0;
				tb_i[seq][pos_seq][pos_target]=0;
				tb_d[seq][pos_seq][pos_target]=0;
			}
	}

	/*Initialise Viterbi matrices*/

	for (seq=1, max_r = SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_seq=1;pos_seq<=l2;pos_seq++){

			/*Viterbi Ms*/
			vt_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
			vt_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
			vt_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
			vt_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];

			if (pos_seq>1) {
				vt_d[seq][pos_seq][1] = (double) vt_m[seq][pos_seq-1][1] + my_pars->ldel;
				tb_d[seq][pos_seq][1] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;
				if ((vt_d[seq][pos_seq-1][1] + my_pars->leps) > vt_d[seq][pos_seq][1]) {
					vt_d[seq][pos_seq][1] = (double) vt_d[seq][pos_seq-1][1] + my_pars->leps;
					tb_d[seq][pos_seq][1] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
				}
			}

			/*Find biggest value prior to recombination*/
			if (vt_m[seq][pos_seq][1] > max_r) {
				max_r = vt_m[seq][pos_seq][1];
				who_max = seq;
				state_max = 1;
				pos_max = pos_seq;
			}
			if (vt_i[seq][pos_seq][1] > max_r) {
				max_r = vt_i[seq][pos_seq][1];
				who_max = seq;
				state_max = 2;
				pos_max = pos_seq;
			}
		}
	}



	/*Now loop forward over positions*/

	for (pos_target=2;pos_target<=l1;pos_target++) {

		/*Match, Insert and Delete Matrices*/
		for (seq=1, max_rn = SMALL+max_r; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match*/
				/*Initialise with recombination*/
				vt_m[seq][pos_seq][pos_target] = max_r + my_pars->lrho + my_pars->lpiM - my_pars->lsizeL;
				tb_m[seq][pos_seq][pos_target] = (double) who_max*10 + state_max + pos_max/my_matrices->tb_divisor;

				/*Compare to MM*/
				if ((vt_m[seq][pos_seq-1][pos_target-1] + my_pars->lmm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target-1] + my_pars->lmm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;
				}
				/*Compare to IM*/
				if ((vt_i[seq][pos_seq-1][pos_target-1] + my_pars->lgm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_i[seq][pos_seq-1][pos_target-1] + my_pars->lgm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 2 + (pos_seq-1)/my_matrices->tb_divisor;
				}

				/*Compare to DM*/
				if ((vt_d[seq][pos_seq-1][pos_target-1] + my_pars->ldm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target-1] + my_pars->ldm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
				}

				/*Add in state match*/
				vt_m[seq][pos_seq][pos_target] += my_pars->lsm[s1[pos_target]][s2[pos_seq]];


				/*Insert*/
				/*Initialise with recombination*/
				vt_i[seq][pos_seq][pos_target] = max_r + my_pars->lrho + my_pars->lpiI - my_pars->lsizeL;
				tb_i[seq][pos_seq][pos_target] = (double) who_max*10 + state_max + pos_max/my_matrices->tb_divisor;

				/*Compare to MI*/
				if ((vt_m[seq][pos_seq][pos_target-1] + my_pars->ldel) > vt_i[seq][pos_seq][pos_target]) {
					vt_i[seq][pos_seq][pos_target] = vt_m[seq][pos_seq][pos_target-1] + my_pars->ldel;
					tb_i[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq)/my_matrices->tb_divisor;
				}
				/*Compare to II*/
				if ((vt_i[seq][pos_seq][pos_target-1] + my_pars->leps) > vt_i[seq][pos_seq][pos_target]) {
					vt_i[seq][pos_seq][pos_target] = vt_i[seq][pos_seq][pos_target-1] + my_pars->leps;
					tb_i[seq][pos_seq][pos_target] = (double) seq*10 + 2 + (pos_seq)/my_matrices->tb_divisor;
				}

				/*Add in state insert*/
				vt_i[seq][pos_seq][pos_target] += my_pars->lsi[s1[pos_target]];



				/*Delete*/
				if (pos_target< l1 && pos_seq>1) {
					/*Initialise with match*/
					vt_d[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target] + my_pars->ldel;
					tb_d[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;

					/*Compare to DD*/
					if ((vt_d[seq][pos_seq-1][pos_target] + my_pars->leps) > vt_d[seq][pos_seq][pos_target]) {
						vt_d[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target] + my_pars->leps;
						tb_d[seq][pos_seq][pos_target] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
					}
				}


				/*Get new max_r, etc.*/
				if (vt_m[seq][pos_seq][pos_target] > max_rn) {
					max_rn = vt_m[seq][pos_seq][pos_target];
					who_max_n = seq;
					state_max_n = 1;
					pos_max_n = pos_seq;
				}
				if (vt_i[seq][pos_seq][pos_target] > max_rn) {
					max_rn = vt_i[seq][pos_seq][pos_target];
					who_max_n = seq;
					state_max_n = 2;
					pos_max_n = pos_seq;
				}
			}
		}

		/*Completed position in target sequence*/
		max_r = max_rn;
		who_max = who_max_n;
		state_max = state_max_n;
		pos_max = pos_max_n;
	}
	/*Completed Viterbi matrices*/
	if (DEBUG>1) {
		print_forward_matrices(my_data, my_pars, my_matrices, target, stdout);
		print_backward_matrices(my_data, my_pars, my_matrices, target, stdout);
	}


	/*Now print maximum likelihood*/
	printf("\nMaximum log-likelihood = %.5lf\n", max_rn);
	my_matrices->llk = max_rn;

	/*Now reconstruct ML path*/
	cp = 2*my_matrices->maxl;
	my_matrices->maxpath_copy[cp] = who_max;
	my_matrices->maxpath_state[cp] = state_max;
	my_matrices->maxpath_pos[cp] = pos_max;

	if (DEBUG) printf("\nEnding in state %i for sequence %i at position %i", state_max, who_max, pos_max);

	pos_target=l1;
	while(pos_target>=1) {

		if (state_max==1) {
			who_next = (int) (tb_m[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_m[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_m[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}
		else if (state_max==2) {
			who_next = (int) (tb_i[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_i[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_i[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}
		else if (state_max==3) {
			who_next = (int) (tb_d[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_d[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_d[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}

		cp--;
		if (cp<=0) {
			printf("\n\n***Error: reconstructed path longer than maximum possible ***\n\n");
			exit(1);
		}
		my_matrices->maxpath_copy[cp] = who_next;
		my_matrices->maxpath_state[cp] = state_next;
		my_matrices->maxpath_pos[cp] = pos_next;

		who_max = who_next;
		state_max = state_next;
		pos_max = pos_next;

		if (my_matrices->maxpath_state[cp+1] != 3) pos_target--;

		if (DEBUG) printf("\nNext = state %i for sequence %i at position %i.  Pos target = %i", state_max, who_max, pos_max, pos_target);

	}

	if (DEBUG) {

		printf("\n\nMax Likelihood alignment\n\nWho\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_copy[i]);
		printf("\nState\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_state[i]);
		printf("\nPos\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_pos[i]);
		printf("\n\n");
	}

	end_vt = clock();
	printf("\n*** Aligned sequence in %.3lf secs CPU time! ***\n", ((double) (end_vt-start_vt)/CLOCKS_PER_SEC));

	/*Reset copy state*/
	my_matrices->who_copy[target]=tmp_copy;

	ofp = fopen(my_pars->alignment_file, "a");
	fprintf(ofp,"\nTarget: %s\tLength: %i\tMLlk: %.3lf\n",my_data->seqs[target]->name, my_data->seqs[target]->length, my_matrices->llk);


	/*First print target sequence*/
	cp++;
	strncpy(tmp_name, my_data->seqs[target]->name, 15);
	tmp_name[15]='\0';
	fprintf(ofp,"%15s\t", tmp_name);
	for (i=cp,pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==3) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[target]->seq[pos_target], my_data->type));
			pos_target++;
		}
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do matching*/
	for (i=1;i<=15;i++) fprintf(ofp," ");
	fprintf(ofp,"\t");
	for (i=cp, pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==1) {
			if (my_data->seqs[target]->seq[pos_target] == my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]]) fprintf(ofp,"|");
			else fprintf(ofp," ");
			pos_target++;
		}
		else if (my_matrices->maxpath_state[i]==2) {
			pos_target++;
			fprintf(ofp,"^");
		}
		else fprintf(ofp,"~");
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do copy tracks - switch whenever gets to new value*/
	strncpy(tmp_name, my_data->seqs[my_matrices->maxpath_copy[cp]]->name, 15);
	tmp_name[15]='\0';
	fprintf(ofp,"%15s\t",tmp_name);
	for (i=cp;i<=2*my_matrices->maxl;i++) {

		/*Check to see if need to make recombination event*/
		if (i>cp && my_matrices->maxpath_copy[i]!=my_matrices->maxpath_copy[i-1]) {
			strncpy(tmp_name, my_data->seqs[my_matrices->maxpath_copy[i]]->name, 15);
			tmp_name[15]='\0';
			fprintf(ofp,"\n%15s\t", tmp_name);
			for (j=1;j<=(i-cp);j++) fprintf(ofp," ");
		}

		if (my_matrices->maxpath_state[i]==2) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]], my_data->type));
		}
	}

	fprintf(ofp,"\n\n");
	fflush(ofp);
	fclose(ofp);

	return;


}
