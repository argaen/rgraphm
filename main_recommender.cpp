#include <iostream>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <iomanip>
#include <string>
#include <sys/time.h>


#include "Node.h"
#include "Link.h"
#include "Group.h"
#include "utils.h"

#define STEPS 100
#define LOGSIZE 5000

typedef boost::unordered_map<int, double> LnFactList;
typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map<Node*, Node*> Queries;
typedef boost::unordered_map<std::string, int> IdMap;
typedef boost::multi_array<int, 3> GGLinks;
typedef boost::multi_array<double, 2> Scores;
typedef boost::numeric::ublas::vector<int> IVector;
typedef std::vector< boost::tuple<Node*, Node*> > tuple_list;


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
double hkState(int K, Groups *g1, Groups *g2, int d1_size, int d2_size, GGLinks *gglinks, int ng1, int ng2){
	double H = 0.0;

	for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
		for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
			H += group2GroupH(&(it1->second), &(it2->second), K, gglinks);
		
	H -= gsl_sf_lnfact(d1_size - ng1) + gsl_sf_lnfact(d2_size - ng2);
	
	return H;

}


/* ##################################### */
void createRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, int K, GGLinks *gglinks, int nnod1, int nnod2){
	int nweight[K];
	int ngrouplinks;
	int group, i = 0;
    Groups::iterator g1, g2;
    GroupNodes::iterator it1, it2;
    Links::iterator nit;
    Hash_Map::iterator itn;

    i = 0;
    for (itn = d1->begin(); itn != d1->end(); itn++){
		group = i;
		itn->second.setGroup(group);
		((*groups1)[group]) = Group(group, K);
		(*groups1)[group].members[itn->second.getId()]=&(itn->second);
        i++;
	}
	
    i = 0;
    for (itn = d2->begin(); itn != d2->end(); itn++){
		group = i;
		itn->second.setGroup(group);
		((*groups2)[group]) = Group(group, K);
		(*groups2)[group].members[itn->second.getId()]=&(itn->second);
        i++;
	}


	for (g1 = groups1->begin(); g1 != groups1->end(); g1++){
        for (g2 = groups2->begin(); g2 != groups2->end(); g2++){
            memset(nweight, 0, sizeof(int)*K);
            ngrouplinks = 0;

            //Groups 1 info filling

			for (it1 = g1->second.members.begin(); it1 != g1->second.members.end(); it1++)
		        for (nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); nit++){
		            for (it2 = g2->second.members.begin(); it2 != g2->second.members.end(); it2++)
		                if (nit->second.getId() == it2->second->getId()){
        		            ++ngrouplinks;
                            nweight[nit->second.getWeight()]++;
                        }
        		}

            (*gglinks)[g1->second.getId()][g2->second.getId()][0] = ngrouplinks;
		    for(i=1; i<K+1; ++i){
                (*gglinks)[g1->second.getId()][g2->second.getId()][i] = nweight[i-1];
            }
		}
	}
}

/* ##################################### */

int mcStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinks, int decorStep, IVector *keys1, IVector *keys2, int *ng1, int *ng2){


    struct timeval stop, start;
    Groups *g;
    Hash_Map *d_move, *d_nomove;
    int bnnod = (nnod1>nnod2) ? nnod1+1 : nnod2+1;

    int factor;
	bool visitedgroup[bnnod]; 
	Group *src_g, *dest_g;
	int newgrp, oldgrp, dice, set_size_move, id;
    bool set_ind;
	Node *n;
	double dH;
    int *ng;
    double set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);
    IVector *keys;


    for (int i=0; i < bnnod; i++)
        visitedgroup[i]=false;

    factor = (nnod1+nnod2)*decorStep;
    for (int move=0; move<factor; move++) {
        /* gettimeofday(&start, NULL); */
        dH = 0.0;
        if (gsl_rng_uniform(rgen) < set_ratio){
            g = g1;d_move = d1;d_nomove = d2;set_ind = true;keys = keys1; ng = ng1;
        }else{
            g = g2;d_move = d2;d_nomove = d1;set_ind = false;keys = keys2; ng = ng2;
        }
        set_size_move = d_move->size();


        dice = gsl_rng_uniform(rgen) * set_size_move;
        n = &(d_move->at((*keys)[dice]));
        oldgrp = n->getGroup();
        newgrp = gsl_rng_uniform(rgen) * set_size_move;
        while ( newgrp == oldgrp ) 
            newgrp = gsl_rng_uniform(rgen) * set_size_move;

        src_g = &(*g)[oldgrp];
        dest_g = &(*g)[newgrp];

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (!visitedgroup[id]){
                if (set_ind){
                    dH -= lnfactlist[(*gglinks)[src_g->getId()][id][0] + K - 1];
                    dH -= lnfactlist[(*gglinks)[dest_g->getId()][id][0] + K - 1];
                }else{
                    dH -= lnfactlist[(*gglinks)[id][src_g->getId()][0] + K - 1];
                    dH -= lnfactlist[(*gglinks)[id][dest_g->getId()][0] + K - 1];
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH -= -lnfactlist[(*gglinks)[src_g->getId()][id][i]];
                        dH -= -lnfactlist[(*gglinks)[dest_g->getId()][id][i]];
                    }else{
                        dH -= -lnfactlist[(*gglinks)[id][src_g->getId()][i]];
                        dH -= -lnfactlist[(*gglinks)[id][dest_g->getId()][i]];
                    }
                }
                visitedgroup[id] = true;
            }
        }

        if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
            dH += lnfactlist[set_size_move - *ng];

        if (set_ind){
            src_g->removeNodeS1(n, d_nomove, gglinks);
            dest_g->addNodeS1(n, d_nomove, gglinks);
        }else{
            src_g->removeNodeS2(n, d_nomove, gglinks);
            dest_g->addNodeS2(n, d_nomove, gglinks);
        }

        if (src_g->members.size() == 0) *ng -=1;
        if (dest_g->members.size() == 1) *ng +=1;

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (visitedgroup[id]){
                if (set_ind){
                    dH += lnfactlist[(*gglinks)[src_g->getId()][id][0] + K - 1];
                    dH += lnfactlist[(*gglinks)[dest_g->getId()][id][0] + K - 1];
                }else{
                    dH += lnfactlist[(*gglinks)[id][src_g->getId()][0] + K - 1];
                    dH += lnfactlist[(*gglinks)[id][dest_g->getId()][0] + K - 1];
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH += -lnfactlist[(*gglinks)[src_g->getId()][id][i]];
                        dH += -lnfactlist[(*gglinks)[dest_g->getId()][id][i]];
                    }else{
                        dH += -lnfactlist[(*gglinks)[id][src_g->getId()][i]];
                        dH += -lnfactlist[(*gglinks)[id][dest_g->getId()][i]];
                    }
                }
                visitedgroup[id] = false;
            }
        }	

        if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
            dH -= lnfactlist[set_size_move - *ng];

        
        if ( dH <= 0.0 || gsl_rng_uniform(rgen) < exp(-dH)){
            *H += dH;
        }else{
            if(set_ind){
                dest_g->removeNodeS1(n, d_nomove, gglinks);
                src_g->addNodeS1(n, d_nomove, gglinks);
            }else{
                dest_g->removeNodeS2(n, d_nomove, gglinks);
                src_g->addNodeS2(n, d_nomove, gglinks);
            }
            if (src_g->members.size() == 1) *ng +=1;
            if (dest_g->members.size() == 0) *ng -=1;
        }
        /* gettimeofday(&stop, NULL); */
        /* printf("Time %lu\n", stop.tv_usec - start.tv_usec); */
    } //End of MC step
    /* double TH; */
    /* TH = hkState(K, g1, g2, nnod1, nnod2, gglinks); */
    /* std::cout << std::setprecision(20) << *H << "    " << TH << "\n"; */
    /* printf("############GROUPS1################\n"); */
    /* printGroups(*g1, K); */
    /* printf("############GROUPS2################\n"); */
    /* printGroups(*g2, K); */
    /* printf("###################################\n"); */


	return 0;
}

void thermalizeMCKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups, int decorStep, IVector *keys1, IVector *keys2, int *ng1, int *ng2){

    double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
    int nrep=20;
    int equilibrated=0;

    Hvalues = (double*)calloc(nrep, sizeof(double));

    do {
        for (int i=0; i<nrep; i++){
            mcStepKState(g1, g2, d1, d2, rgen, H, K, lnfactlist, logsize, nnod1, nnod2, gglinksGroups, decorStep, keys1, keys2, ng1, ng2);
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

int getDecorrelationKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups, IVector *keys1, IVector *keys2, int *ng1, int *ng2){

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

    for (i=0; i<nrep; i++) {

        printf("Estimating decorrelation time %d/%d\n", i+1, nrep);
        g1t = (*g1);
        g2t = (*g2);

        for (step=0; step<=x2; step++) {
            mcStepKState(g1, g2, d1, d2, rgen, H, K, lnfactlist, logsize, nnod1, nnod2, gglinksGroups, 1, keys1, keys2, ng1, ng2);

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
	int weight, logsize=LOGSIZE;
    std::string id1, id2;
	double H;
    double *lnfactlist;
    tuple_list queries;
	Hash_Map d1, d2;
	IdMap rd1, rd2;
	Groups groups2, groups1;
	gsl_rng *randomizer;
	char* tFileName;
	char* qFileName;
	int stepseed;
	int mark, k, i, j, n, nk, tmpid1 = 0, tmpid2 = 0, tid1, tid2;
    int nnod1, nnod2, nqueries = 0;
    int ng1 = 0, ng2 = 0;
    std::ofstream outfile;

	parseArguments(argc, argv, &tFileName, &qFileName, &stepseed, &mark);

    
	randomizer = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(randomizer, stepseed);

    fprintf(stdout, "Reading input files and initializing data structures...\n");
	std::ifstream tFile(tFileName);
	if (tFile.is_open()){
		while (tFile >> id1 >> id2 >> weight){
            if (rd1.find(id1) == rd1.end()){
                rd1.insert(IdMap::value_type(id1,tmpid1));
                tmpid1++;
            }
            if (rd2.find(id2) == rd2.end()){
                rd2.insert(IdMap::value_type(id2,tmpid2));
                tmpid2++;
            }
            tid1 = rd1[id1]; 
            tid2 = rd2[id2];
            d1.insert(Hash_Map::value_type(tid1,Node(id1, tid1)));
            d1[tid1].neighbours.insert(Links::value_type(tid2,Link(tid2, weight)));
			d2.insert(Hash_Map::value_type(tid2,Node(id2, tid2)));
			d2[tid2].neighbours.insert(Links::value_type(tid1,Link(tid1, weight)));
		}
		tFile.close();
	}else
		std::cout << "Couldn't open file " << tFileName << "\n";

	std::ifstream qFile(qFileName);
	if (qFile.is_open()){
		while (qFile >> id1 >> id2){
            if (rd1.find(id1) == rd1.end()){
                rd1.insert(IdMap::value_type(id1,tmpid1));
                tid1 = rd1[id1]; 
                d1.insert(Hash_Map::value_type(tid1,Node(id1, tid1)));
                tmpid1++;
            }

            if (rd2.find(id2) == rd2.end()){
                rd2.insert(IdMap::value_type(id2,tmpid2));
                tid2 = rd2[id2];
                d2.insert(Hash_Map::value_type(tid2,Node(id2, tid2)));
                tmpid2++;
            }

                tid1 = rd1[id1]; 
                tid2 = rd2[id2];
            /* if (d1[tid1].neighbours.find(tid2) != d1[tid1].neighbours.end()){ */
                d1[tid1].neighbours.erase(tid2);
                d2[tid2].neighbours.erase(tid1);


            queries.push_back(boost::tuple<Node*, Node*>(&(d1[rd1[id1]]), &(d2[rd2[id2]])));
		}
		qFile.close();
	}else
		std::cout << "Couldn't open file " << tFileName << "\n";

    lnfactlist = genLogFactList(logsize);

    nnod1 = d1.size();
    nnod2 = d2.size();
    ng1 = d1.size();
    ng2 = d2.size();
    nqueries = queries.size();

    IVector keys1(nnod1), keys2(nnod2);
    i = 0;
    for(Hash_Map::iterator it = d1.begin(); it != d1.end(); ++it) {
        keys1[i] = it->first;
        ++i;
    }
    i = 0;
    for(Hash_Map::iterator it = d2.begin(); it != d2.end(); ++it) {
        keys2[i] = it->first;
        ++i;
    }

    GGLinks gglinks(boost::extents[nnod1+1][nnod2+1][mark+1]);
    Scores scores(boost::extents[nqueries][mark]);

    fprintf(stdout, "Creating groups...\n");
	createRandomGroups(&d1, &d2, &groups1, &groups2, mark, &gglinks, nnod1, nnod2);


    fprintf(stdout, "Calculating initial H...\n");
	H = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks, ng1, ng2);
    printf("Initial H: %f\n", H);

    /* decorStep = getDecorrelationKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks, &keys1, &keys2, &ng1, &ng2); */

    /* thermalizeMCKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks, decorStep, &keys1, &keys2, &ng1, &ng2); */

    decorStep = 1;

    /* double TH; */
    int indquery;
	for(i=0; i<STEPS; i++){
		mcStepKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks, decorStep, &keys1, &keys2, &ng1, &ng2);
        /* TH = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks, ng1, ng2); */
        /* std::cout << std::setprecision(20) << H << "    " << TH << "\n\n"; */
        std::cout << std::setprecision(20) << i << " " << H << "\n";
        /* printf("############GROUPS1################\n"); */
        /* printGroups(groups1, mark); */
        /* printf("############GROUPS2################\n"); */
        /* printGroups(groups2, mark); */
        /* printf("###################################\n"); */
        indquery = 0;
        for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
            for (k=1; k<mark+1; k++) {
                nk = gglinks[(*it->get<0>()).getGroup()][(*it->get<1>()).getGroup()][k];
                n = gglinks[(*it->get<0>()).getGroup()][(*it->get<1>()).getGroup()][0];
                scores[indquery][k-1] += (float)(nk + 1) / (float)(n + mark);
            }
            indquery++;
        }


        if (i % 100 == 0) { 
            outfile.open("scores.tmp");
            indquery = 0;
            for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
                outfile << "# " << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << " # ";
                for (k=0; k<mark; k++)
                    outfile << scores[indquery][k]/(double)STEPS << " ";
                indquery++;
                outfile << "\n\n";
            }
            outfile.close();
        }
	}

    for (j=0; j<nqueries; j++)
        for (k=0; k<mark; k++)
            scores[j][k] /= (double)STEPS;

    outfile.open("scores.tmp");
    printf("RESULTS:\n");
    indquery = 0;
    for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
        outfile << "# " << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << " # ";
        std::cout << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << '\n';
        for (k=0; k<mark; k++){
            outfile << scores[indquery][k] << " ";
            std::cout << scores[indquery][k] << " ";
        }
        indquery++;
        outfile << "\n\n";
        std::cout << "\n\n";
    }

    free(lnfactlist);
	gsl_rng_free(randomizer);
	return 0;

}

