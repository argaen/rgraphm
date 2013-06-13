#include <iostream>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <iomanip>
#include <sys/time.h>


#include "Node.h"
#include "Link.h"
#include "Group.h"
#include "utils.h"

#define STEPS 100000
#define LOGSIZE 5000

typedef boost::unordered_map<int, double> LnFactList;
typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::multi_array<int, 3> GGLinks;


/* ##################################### */


/* ##################################### */
double group2GroupH(Group *g1, Group *g2, int K, GGLinks *gglinks){
	double H = 0.0;

	/* H += gsl_sf_lnfact(g1->g2glinks[g2->get_id()][0] + K -1); */
	H += gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][0] + K - 1);
	for(int i=1; i<K+1; ++i){
		/* H -= gsl_sf_lnfact(g1->g2glinks[g2->get_id()][i]); */
        H -= gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][i]);
	}
    return H;
}


/* ##################################### */
double hkState(int K, Groups *g1, Groups *g2, int d1_size, int d2_size, GGLinks *gglinks){
	double H = 0.0;

	for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
		for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
			H += group2GroupH(&(it1->second), &(it2->second), K, gglinks);
		
	H -= gsl_sf_lnfact(d1_size - g1->size()) + gsl_sf_lnfact(d2_size - g2->size());
	
	return H;

}


/* ##################################### */
void createRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, gsl_rng *rgen, int K, bool random, GGLinks *gglinks){
	int nweight[K+1];
	int ngrouplinks;
	int group;

	for (int i = 1; i<(*d1).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d1).size())) + 1) : group=i;
		(*d1)[i].setGroup(group);
		((*groups1)[group]) = Group(group, K);
		(*groups1)[group].members[i]=&((*d1)[i]);
	}
	
	for (int i = 1; i<(*d2).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d2).size())) + 1) : group=i;
		(*d2)[i].setGroup(group);
		((*groups2)[group]) = Group(group, K);
		(*groups2)[group].members[i]=&((*d2)[i]);
	}

	for (Groups::iterator g1 = groups1->begin(); g1 != groups1->end(); ++g1){
        for (Groups::iterator g2 = groups2->begin(); g2 != groups2->end(); ++g2){
        memset(nweight, 0, sizeof(int)*(K+1));
        ngrouplinks = 0;

		//Groups 1 info filling

			for(GroupNodes::iterator it1 = g1->second.members.begin(); it1 != g1->second.members.end(); ++it1)
		        for (Links::iterator nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); ++nit){
		            for(GroupNodes::iterator it2 = g2->second.members.begin(); it2 != g2->second.members.end(); ++it2)
		                if (nit->second.getId() == it2->second->getId()){
        		            ++ngrouplinks;
                            nweight[nit->second.getWeight()]++;
                        }
        		}

            (*gglinks)[g1->second.getId()][g2->second.getId()][0] = ngrouplinks;
		    for(int i=1; i<K+1; ++i)
                (*gglinks)[g1->second.getId()][g2->second.getId()][i] = nweight[i];
		}
	}
}

/* ##################################### */

int mcStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinks, int decorStep){


    /* struct timeval stop, start; */
    /* gettimeofday(&start, NULL); */
    Groups *g;
    Hash_Map *d_move, *d_nomove;
    int bnnod = (nnod1>nnod2) ? nnod1+1 : nnod2+1;

	bool visitedgroup[bnnod]; 
    for (int i=0; i < bnnod; i++)
        visitedgroup[i]=false;
	Group *src_g, *dest_g;
	int newgrp, oldgrp, dice, set_size_move, id;
    bool set_ind;
	Node *n;
	double dH = 0.0;
    double set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);


    /* for (int move=0; move<(nnod1+nnod2)*decorStep; move++) { */
        if (gsl_rng_uniform(stepgen) < set_ratio){
            g = g1;d_move = d1;d_nomove = d2;set_ind = true;
        }else{
            g = g2;d_move = d2;d_nomove = d1;set_ind = false;
        }
        set_size_move = d_move->size();


        dice = gsl_rng_uniform(stepgen) * set_size_move + 1;
        n = &(d_move->at(dice));
        oldgrp = n->getGroup();
        newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;
        while ( newgrp == oldgrp ) 
            newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;

        src_g = &(*g)[oldgrp];
        dest_g = &(*g)[newgrp];

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (!visitedgroup[id]){
                if (set_ind){
                    dH -= logFact((*gglinks)[src_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                    dH -= logFact((*gglinks)[dest_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                }else{
                    dH -= logFact((*gglinks)[id][src_g->getId()][0] + K - 1, logsize, lnfactlist);
                    dH -= logFact((*gglinks)[id][dest_g->getId()][0] + K - 1, logsize, lnfactlist);
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH -= -logFact((*gglinks)[src_g->getId()][id][i], logsize, lnfactlist);
                        dH -= -logFact((*gglinks)[dest_g->getId()][id][i], logsize, lnfactlist);
                    }else{
                        dH -= -logFact((*gglinks)[id][src_g->getId()][i], logsize, lnfactlist);
                        dH -= -logFact((*gglinks)[id][dest_g->getId()][i], logsize, lnfactlist);
                    }
                }
                visitedgroup[id] = true;
            }
        }

        if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
            dH += logFact(set_size_move - g1->size(), logsize, lnfactlist);

        if (set_ind){
            src_g->removeNodeS1(n, d_nomove, gglinks);
            dest_g->addNodeS1(n, d_nomove, gglinks);
        }else{
            src_g->removeNodeS2(n, d_nomove, gglinks);
            dest_g->addNodeS2(n, d_nomove, gglinks);
        }

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (visitedgroup[id]){
                if (set_ind){
                    dH += logFact((*gglinks)[src_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                    dH += logFact((*gglinks)[dest_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                }else{
                    dH += logFact((*gglinks)[id][src_g->getId()][0] + K - 1, logsize, lnfactlist);
                    dH += logFact((*gglinks)[id][dest_g->getId()][0] + K - 1, logsize, lnfactlist);
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH += -logFact((*gglinks)[src_g->getId()][id][i], logsize, lnfactlist);
                        dH += -logFact((*gglinks)[dest_g->getId()][id][i], logsize, lnfactlist);
                    }else{
                        dH += -logFact((*gglinks)[id][src_g->getId()][i], logsize, lnfactlist);
                        dH += -logFact((*gglinks)[id][dest_g->getId()][i], logsize, lnfactlist);
                    }
                }
                visitedgroup[id] = false;
            }
        }	

        if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
            dH -= logFact(set_size_move - g1->size(), logsize, lnfactlist);

        if ( dH <= 0.0 || gsl_rng_uniform(stepgen) < exp(-dH)){
            *H += dH;
        }else{
            if(set_ind){
                dest_g->removeNodeS1(n, d_nomove, gglinks);
                src_g->addNodeS1(n, d_nomove, gglinks);
            }else{
                dest_g->removeNodeS2(n, d_nomove, gglinks);
                src_g->addNodeS2(n, d_nomove, gglinks);
            }
        }
        /* gettimeofday(&stop, NULL); */
        /* printf("Time %lu\n", stop.tv_usec - start.tv_usec); */
    /* } //End of MC step */


	return 0;
}

void thermalizeMCKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups, int decorStep){

    double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
    int nrep=20;
    int equilibrated=0;

    Hvalues = (double*)calloc(nrep, sizeof(double));

    do {
        for (int i=0; i<nrep; i++){
            mcStepKState(g1, g2, d1, d2, stepgen, groupgen, H, K, lnfactlist, logsize, nnod1, nnod2, gglinksGroups, decorStep);
            Hvalues[i] = *H;
        }

        HMean1 = mean(Hvalues, nrep);
        HStd1 = stddev(Hvalues, nrep);

        if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) {
            equilibrated++;
            printf("#\tequilibrated (%d/5) H=%lf\n", equilibrated, HMean1);

        } else {
            printf("#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n", HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep));
            HMean0 = HMean1;
            HStd0 = HStd1;
            equilibrated = 0;
        }
    } while (equilibrated < 5);

    return;
}

int getDecorrelationKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups){

    int nrep = 10, step, norm = 0;
    int x1, x2, i;
    double *decay, *decay1, *decay2;
    double meanDecay, meanDecay1, meanDecay2, sigmaDecay;
    double y11 = 0.0, y12 = 0.0, y21 = 0.0, y22 = 0.0, result;
	Groups g2t, g1t;

    decay1 = (double*) calloc(nrep, sizeof(double));
    decay2 = (double*) calloc(nrep, sizeof(double));

    x2 = (nnod1 + nnod2) / 5;
    if (x2 < 10)
        x2 = 10;
    x1 = x2 / 4;

    for (i=0; i<nrep; i++){

        printf("Estimating decorrelation time %d/%d\n", i, nrep);
        g1t = (*g1);
        g2t = (*g2);

        for (step=0; step<=x2; step++){
            mcStepKState(g1, g2, d1, d2, stepgen, groupgen, H, K, lnfactlist, logsize, nnod1, nnod2, gglinksGroups, 1);

            if (step == x1){
                y11 = mutualInfo(g1, &g1t, nnod1, nnod2);
                y12 = mutualInfo(g2, &g2t, nnod1, nnod2);
            }
        }
        y21 = mutualInfo(g1, &g1t, nnod1, nnod2);
        y22 = mutualInfo(g2, &g2t, nnod1, nnod2);

        (nnod1>1) ? decay1[i] = 2 * getDecay(nnod1, x1, x2, y11, y21) : decay1[i] = 1.e-6;
        (nnod2>1) ? decay2[i] = 2 * getDecay(nnod2, x1, x2, y12, y22) : decay2[i] = 1.e-6;

        printf("# Decorrelation times (estimate %d) = %g %g\n\n", i + 1, decay1[i], decay2[i]);

        if (decay1[i] < 0. || decay2[i] < 0.){
            i--;
            printf("\t# ignoring...\n");
        }
    }

    meanDecay1 = mean(decay1, nrep);
    meanDecay2 = mean(decay2, nrep);

    if (meanDecay1 > meanDecay2){
        meanDecay = meanDecay1;
        decay = decay1;
        sigmaDecay = stddev(decay1, nrep);
    }else{
        meanDecay = meanDecay2;
        decay = decay2;
        sigmaDecay = stddev(decay2, nrep);
    }

    result = meanDecay*nrep;

    for (i=0; i<nrep; i++)
       (fabs(decay[i] - meanDecay) / sigmaDecay > 2) ? result -= decay[i] : ++norm;

    printf("# Decorrelation step: %d\n\n", (int)(result / norm + 0.5));

    return (int)(result / norm + 0.5);
}




/* ##################################### */
int main(int argc, char **argv){

    int decorStep;
	int id1, id2, weight, logsize=LOGSIZE;
    int norm=0;
	double H;
    double *lnfactlist;
	Hash_Map d1, d2;
	Hash_Map d1c, d2c;
	Groups groups2, groups1;
	gsl_rng *groups_randomizer, *step_randomizer;
	char* tFileName;
	char* qFileName;
	int groupseed, stepseed;
	int mark;
    int nnod1, nnod2;

	parseArguments(argc, argv, &tFileName, &qFileName, &stepseed, &groupseed, &mark);

    std::cout << stepseed << " " << tFileName << " " << qFileName << " " << groupseed << " " << mark << "\n";

    
	step_randomizer = gsl_rng_alloc(gsl_rng_mt19937);
	groups_randomizer = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(step_randomizer, stepseed);
	gsl_rng_set(groups_randomizer, groupseed);

	std::ifstream inFile(tFileName);
	if (inFile.is_open()){
		while (inFile >> id1 >> id2 >> weight){
			d1.insert(Hash_Map::value_type(id1,Node(id1)));
			d1[id1].neighbours.insert(Links::value_type(id2,Link(id2, weight)));
			d2.insert(Hash_Map::value_type(id2,Node(id2)));
			d2[id2].neighbours.insert(Links::value_type(id1,Link(id1, weight)));
			

		}
		inFile.close();

	}else
		std::cout << "Couldn't open file " << tFileName << "\n";
    
    lnfactlist = genLogFactList(logsize);


	d1c = d1;
	d2c = d2;	
    nnod1 = d1.size();
    nnod2 = d2.size();

    double **scores = calloc(mark*nqueries, sizeof(double));

    GGLinks gglinks(boost::extents[nnod1+2][nnod2+2][mark+1]);



	createRandomGroups(&d1c, &d2c, &groups1, &groups2, groups_randomizer, mark, false, &gglinks);

	H = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks); //Initialize H with initial system energy
	std::cout << "Initial H: "<< H <<"\n";

    decorStep = getDecorrelationKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks);

    thermalizeMCKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks, decorStep);


    /* printf("DECORRELATION: %d\n", decorStep); */

    double TH;
	for(int i=0; i<STEPS; i++){
		mcStepKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks, decorStep);
        /* TH = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks); */
        /* std::cout << std::setprecision(20) << H << "    " << TH << "\n"; */
        /* std::cout << std::setprecision(20) << H << "    " << "\n"; */
        /* printf("############GROUPS1################\n"); */
        /* printGroups(groups1, mark); */
        /* printf("############GROUPS2################\n"); */
        /* printGroups(groups2, mark); */
        /* printf("###################################\n"); */
        /* ++norm; */
        /* for (k=0; k<K; k++) { */
        /*   for (q=0; q<nquery; q++) { */
        /* nk = G1G2[k][querySet[q]->n1->inGroup][querySet[q]->n2->inGroup]; */
        /* n = 0; */
        /* for (k2=0; k2<K; k2++) */
        /*   n += G1G2[k2][querySet[q]->n1->inGroup][querySet[q]->n2->inGroup]; */
        /* score[k][q] += (float)(nk + 1) / (float)(n + K); */
        /*   } */
        /* } */
	}

    free(lnfactlist);
	gsl_rng_free(step_randomizer);
	gsl_rng_free(groups_randomizer);
	return 0;

}

