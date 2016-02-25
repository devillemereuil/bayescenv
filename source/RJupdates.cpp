//    This program, BayeScEnv, aims at detecting genetics markers under local adaptation,
//	based on allele frequency differences between population and environmental differentiation
//    Copyright (C) 2010  Matthieu Foll for the shared code with BayeScan 2.1
//    Copyright (C) 2014  Pierre de Villemereuil for the new code
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "global_defs.h"

#include <algorithm>
#include <math.h>
#include <iostream>

#include <omp.h>

//#define prop_alpha_mean 0
//#define prop_alpha_var 5

using namespace std;

//Jump between the 3 models for dominant (beta, alpha+beta, beta+gamma)
void jump_3model()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double u; // random uniform
    double old_alpha; //old value of alpha
    double old_g; //old value of g
    double non_zero_g;
    double non_zero_alpha;
    int old_model; //old value of model
    double new_theta,old_theta; // new and old values of theta

    #pragma omp parallel for SCHED_I reduction(+:/*nb_alpha_included,*/log_likelihood) private(r,A,non_zero_g,non_zero_alpha,old_model,old_alpha,old_g,new_theta,old_theta)
    for (int i=0; i<I; i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            old_model=model[i];

            u = randgen_parallel[omp_get_thread_num()].randDblExc();
            if (u<0.5)
                model[i] = (model[i]+1) % 3;
            else
                model[i] = (model[i]+2) % 3;

            old_alpha=alpha[i];
            old_g=g[i];

            non_zero_alpha=0;
            non_zero_g=0;

            if (old_alpha!=0)
                non_zero_alpha=old_alpha;
            if (old_g!=0)
                non_zero_g=old_g;

            //If we draw beta only
            if (model[i]==0) {
                alpha[i]=0;
                g[i]=0;
            }

            //If we draw alpha+beta
            if (model[i]==1) {
                do {
                    alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
                } while(alpha[i]==0);  //Infinitesimal chance that this happens, but we should be sure for the code to be clean
                g[i]=0;
            }

            //If we draw beta+gamma
            if (model[i]==2) {
                alpha[i]=0;
                do {
                    g[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_g[i],sqrt(var_g[i]));
                } while (g[i]<=0||g[i]>U);
            }

            if (alpha[i]!=0)
                non_zero_alpha=alpha[i];
            if (g[i]!=0)
                non_zero_g=g[i];

// calculate A
            A=0;

            for (int j=0; j<J; j++)
            {
                // calculate old and new value of theta
                new_theta=exp(-(alpha[i]+beta[j]+g[i]*X[j]));
                old_theta=exp(-(old_alpha+beta[j]+old_g*X[j]));
       

                A+= gammaln(new_theta)-gammaln(old_theta)
                    - gammaln(new_theta*freq_ancestral[i]) + gammaln(old_theta*freq_ancestral[i])
                    - gammaln(new_theta*(1-freq_ancestral[i])) + gammaln(old_theta*(1-freq_ancestral[i]))
                    + freq_ancestral[i]*(new_theta-old_theta)*log(pop[j].locus[i].p)
                    + (1-freq_ancestral[i])*(new_theta-old_theta)*log(1-pop[j].locus[i].p);
            }

            //Account for g changes
            if (g[i]!=old_g) {
                if(g[i]==0) { //If we removed g
                    A+= (-log_trunc[i]-((non_zero_g-mean_g[i])*(non_zero_g-mean_g[i]))/(2*var_g[i]))
                        -log_prior_g(non_zero_g);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));
                }
                else //If we added g
                {
                     A+= log_prior_g(non_zero_g)//-0.5*log(2*M_PI*sd_prior_g*sd_prior_g)-(g[i]*g[i])/(2*sd_prior_g*sd_prior_g)
                         -(-log_trunc[i]			//log of normalizing constant for truncated normal (include the line commented below) //revision env1 (PdV)
                         -((non_zero_g-mean_g[i])*(non_zero_g-mean_g[i]))/(2*var_g[i]));
                }
	    }
            

            //Account for alpha changes
            if (alpha[i]!=old_alpha) {
                if (alpha[i]==0) { //If we removed alpha
                    A+=(-0.5*log(2*M_PI*var_alpha[i])-((non_zero_alpha-mean_alpha[i])*(non_zero_alpha-mean_alpha[i]))/(2*var_alpha[i]))
                       -log_prior_alpha(non_zero_alpha);;
                }
                else
                {
                    A+=-(-0.5*log(2*M_PI*var_alpha[i])-((non_zero_alpha-mean_alpha[i])*(non_zero_alpha-mean_alpha[i]))/(2*var_alpha[i]))
                       +log_prior_alpha(non_zero_alpha);;
                }
            }

            //calculate prior transition probabilities
            if (old_model+model[i]==1) { //If we use model 0 and model 1
                if (model[i]-old_model>0) 	//From O to 1
                    A+= log_trans_prob_0_1;
                else				//From 1 to 0
                    A+= -log_trans_prob_0_1;
            }
            else if (old_model+model[i]==2) { //If we use model 0 and model 2
                if (model[i]-old_model>0)
                    A+= log_trans_prob_0_2;
                else				//From 2 to 0
                    A+= -log_trans_prob_0_2;
            }
            else if (old_model+model[i]==3) { //If we use model 1 and model 2
                if (model[i]-old_model>0) 	//From 1 to 2
                    A+= log_trans_prob_1_2;
                else				//From 2 to 1
                    A+= -log_trans_prob_1_2;
            }

            r=randgen_parallel[omp_get_thread_num()].randDblExc();
	    
// reject proposed value
            if (log(r)>A)
            {
                g[i]=old_g;
                alpha[i]=old_alpha;
                model[i]=old_model;
            }

            //Change the "included" vectors accordingly to the new (or old) values
            if(g[i]==0)
                g_included[i]=false;
            else
                g_included[i]=true;

            if (alpha[i]==0)
                alpha_included[i]=false;
            else
                alpha_included[i]=true;

        }
    }
}

//Jump between the 3 models for codominant data (beta, alpha+beta, beta+gamma)
void jump_3model_codominant()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double u; // random uniform
    double old_alpha; //old value of alpha
    double old_g; //old value of g
    double non_zero_g;
    double non_zero_alpha;
    int old_model; //old value of model
    double new_theta,old_theta; // new and old values of theta

    //double old_log_likelihood; // old loglikelihood
    double diff_log_likelihood; // old loglikelihood
    #pragma omp parallel for SCHED_I reduction(+:/*nb_alpha_included,*/log_likelihood) private(r,A,non_zero_g,non_zero_alpha,old_model,old_alpha,old_g,new_theta,old_theta,diff_log_likelihood)
    for (int i=0; i<I; i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
            old_model=model[i];

            u = randgen_parallel[omp_get_thread_num()].randDblExc();
            if (u<0.5)
                model[i] = (model[i]+1) % 3;
            else
                model[i] = (model[i]+2) % 3;

            old_alpha=alpha[i];
            old_g=g[i];

            non_zero_alpha=0;
            non_zero_g=0;

            if (old_alpha!=0)
                non_zero_alpha=old_alpha;
            if (old_g!=0)
                non_zero_g=old_g;

            //If we draw beta only
            if (model[i]==0) {
                alpha[i]=0;
                g[i]=0;
            }

            //If we draw alpha+beta
            if (model[i]==1) {
                do {
                    alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
                }
                while(alpha[i]==0);  //Infinitesimal chance that this happens, but we should be sure for the code to be clean
                g[i]=0;
            }

            //If we draw beta+gamma
            if (model[i]==2) {
                alpha[i]=0;
		do {
		    g[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_g[i],sqrt(var_g[i]));
		} while (g[i]<=0||g[i]>U);
            }

            if (alpha[i]!=0)
                non_zero_alpha=alpha[i];
            if (g[i]!=0)
                non_zero_g=g[i];

// calculate A
            A=0;

            double old_l=0;
            for (int j=0; j<J; j++)
            {
                // calculate old and new value of theta
                old_theta=exp(-(old_alpha+beta[j]+old_g*X[j]));

                old_l+=gammaln(old_theta)-gammaln(pop[j].locus[i].alleleCount+old_theta);
                for (int k=0; k<pop[j].locus[i].ar; k++)
                    old_l+=gammaln(pop[j].locus[i].data_allele_count[k]+old_theta*freq_locus[i].allele[k])
                           -gammaln(old_theta*freq_locus[i].allele[k]);
            }

            double new_l=0;
            for (int j=0; j<J; j++)
            {
                // calculate old and new value of theta
                new_theta=exp(-(alpha[i]+beta[j]+g[i]*X[j]));

                new_l+=gammaln(new_theta)-gammaln(pop[j].locus[i].alleleCount+new_theta);
                for (int k=0; k<pop[j].locus[i].ar; k++)
                    new_l+=gammaln(pop[j].locus[i].data_allele_count[k]+new_theta*freq_locus[i].allele[k])
                           -gammaln(new_theta*freq_locus[i].allele[k]);
            }
            diff_log_likelihood=-old_l+new_l;
            A=diff_log_likelihood;

            //Account for g changes
            if (g[i]!=old_g) {
                if(g[i]==0) { //If we removed g
                        A+= (-log_trunc[i]-((non_zero_g-mean_g[i])*(non_zero_g-mean_g[i]))/(2*var_g[i]))
                            -log_prior_g(non_zero_g);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));
                }
                else //If we added g
                {
                        A+= log_prior_g(non_zero_g)//-0.5*log(2*M_PI*sd_prior_g*sd_prior_g)-(g[i]*g[i])/(2*sd_prior_g*sd_prior_g)
                            -(-log_trunc[i]			//log of normalizing constant for truncated normal (include the line commented below) //revision env1 (PdV)
                              -((non_zero_g-mean_g[i])*(non_zero_g-mean_g[i]))/(2*var_g[i]));
                }
            }

            //Account for alpha changes
            if (alpha[i]!=old_alpha) {
                if (alpha[i]==0) { //If we removed alpha
                    A+=(-0.5*log(2*M_PI*var_alpha[i])-((non_zero_alpha-mean_alpha[i])*(non_zero_alpha-mean_alpha[i]))/(2*var_alpha[i]))
                       -log_prior_alpha(non_zero_alpha);;
                }
                else
                {
                    A+=-(-0.5*log(2*M_PI*var_alpha[i])-((non_zero_alpha-mean_alpha[i])*(non_zero_alpha-mean_alpha[i]))/(2*var_alpha[i]))
                       +log_prior_alpha(non_zero_alpha);;
                }
            }

            //calculate prior transition probabilities
            if (old_model+model[i]==1) { //If we use model 0 and model 1
                if (model[i]-old_model>0) 	//From O to 1
                    A+= log_trans_prob_0_1;
                else				//From 1 to 0
                    A+= -log_trans_prob_0_1;
            }
            else if (old_model+model[i]==2) { //If we use model 0 and model 2
                if (model[i]-old_model>0)
                    A+= log_trans_prob_0_2;
                else				//From 2 to 0
                    A+= -log_trans_prob_0_2;
            }
            else if (old_model+model[i]==3) { //If we use model 1 and model 2
                if (model[i]-old_model>0) 	//From 1 to 2
                    A+= log_trans_prob_1_2;
                else				//From 2 to 1
                    A+= -log_trans_prob_1_2;
// 		    cout << "A 2->1 : " << exp(A) << endl;
            }

            r=randgen_parallel[omp_get_thread_num()].randDblExc();

// 		cout << "From model " << old_model << " to model " << model[i] << ", log(r)=" << log(r) << ", A=" << A << ", diffLog= "<< diff_log_likelihood << ", P="<< P << ", non_zero_alpha="<< non_zero_alpha <<endl;

// reject proposed value
            if (log(r)>A)
            {
                g[i]=old_g;
                alpha[i]=old_alpha;
                model[i]=old_model;
            }
            else
            {
                log_likelihood=log_likelihood+diff_log_likelihood;
            }

            //Change the "included" vectors accordingly to the new (or old) values
            if(g[i]==0)
                g_included[i]=false;
            else
                g_included[i]=true;

            if (alpha[i]==0)
                alpha_included[i]=false;
            else
                alpha_included[i]=true;

        }
    }
}
