#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_gamma.h>
#include <string.h>

#include "Group.h"

typedef boost::unordered_map<int, Group> Groups;

struct pair
{
  double x0;
  double y0;
  double x1;
  double y1;
};

void parseArguments(int argc, char **argv, char** inFile, char** qFile, int* stepseed, int* mark);
void printGroups(Groups g, int mark);
int ExponentialRootF(const gsl_vector *params, void *points, gsl_vector *f);
double getDecay(int nnod, double x1, double x2, double y1, double y2);
double mean(double* data, int N);
double stddev(double *data, int N);
double mutualInfo(Groups *g, Groups *gt, int nnod1, int nnod2);
double logFact(int key, int size, double* logFactList);
double* genLogFactList(int size);

#endif
