/******************************************************************************
* A genetic algorithm for item selection in the context of Mokken scale       *
* analysis. This program is accessible in the R package mokken.               *
* The functions in this file are ordered from specific to general.            *
*                                                                             *
* Author: J.H. Straat                                                         *
* Rewritten to Rcpp by D. van den Bergh
*                                                                             *
* General outline of the algorithm:                                           *
*                                                                             *
* 1   Initial population                                                      *
* 1.1 Randomly draw the initial population                                    *
* 1.2 Evaluate whether the initial population satisfies the Mokken scale      *
*       criteria.                                                             *
* 1.3 Store the best partitioning of the initial population                   *
*                                                                             *
* 2   Repeat the steps in this part until no new 'best partitioning' is found *
*       in 'maxgens' generations.                                             *
* 2.1 Select a new population from the members of the old population.         *
*       The probability of being selected is related to the member's fitness. *
* 2.2 Crossover. Exchange two equivalent subvectors of two members.      *
* 2.3 Mutation. Randomly change some of the numbers.                          *
* 2.4 Evaluate whether the new population satisfies the Mokken scale criteria *
* 2.5 Compare the best partitioning of the new population with the stored     *
*     best partitioning.                                                      *
* 2.6 If the stored best partitioning is not contained in the new population, *
*     replace the worst partitioning of the new population by the stored best *
*     partitioning.                                                           *
*                                                                             *
******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/*************************************************
* ScaleNumItems. Determine how many items are in *
* each scale.                                    *
*************************************************/

void ScaleNumItemsRcpp(const int mem, const int nclus, IntegerVector & NUMITEMS, const int NITEM, IntegerMatrix & pop){

	int i, j; /* loop indicators */

	/* Fix all values to zero */
	NUMITEMS.fill(0);

	/* Count the number of items in scale j */
	for(i=0; i<NITEM; i++)
		for(j=1; j<(nclus+1); j++)
			if(pop[i + mem*NITEM] == j)  NUMITEMS[j-1]++;

	/**********************************************************************
	* if the number of items in scale j equals 1, then the item belonging *
	* to that scale will be assigned to scale "zero" which is used        *
	* as a scale for the non-scalable items                               *
	**********************************************************************/
	for(i=0; i<nclus; i++)
		if(NUMITEMS[i] == 1)
			for(j=0; j<NITEM; j++)
				if(pop[j + mem*NITEM] == (i+1))
				{
					pop[j + mem*NITEM] = 0;
					NUMITEMS[i] = 0;
				}
}



/**********************************************
* NumScales. Count the total number of scales *
**********************************************/

int NumScalesRcpp(const int nclus, const IntegerVector & NUMITEMS)
{

	int i;				/*loop indicator */
	int NSCALES = 0;	/* initialize number of scales to zero */

	/* count the number of clusters with at least two items */
	for(i=0; i<nclus; i++) if(NUMITEMS[i] > 1) NSCALES++;

	return(NSCALES);
}



/***********************************************************
* ScaleItems. A matrix which indicates which items belongs *
* to which scale.                                          *
***********************************************************/

void ScaleItemsRcpp(const int mem, const int nclus, const int NITEM,
					IntegerMatrix& pop, IntegerMatrix& ITEMS, IntegerVector& NUMITEMS)
{
	int i,j,k; /* loop indicators */

	ITEMS.fill(0);

	/* a matrix which denotes which item belongs to which scale  */
	for(i=0; i<nclus; i++)
		if(NUMITEMS[i] > 1)
		{
			k=0;
			for(j=0; j<NITEM; j++)
			{
				if(pop[j + mem*NITEM] == (i+1))
				{
					ITEMS[k + i*NITEM] = j;
					k++;
				}
			}
		}
}



/************************************************
* SortHi. Ascendingly sort the Hi-coefficients. *
************************************************/

void sortHiRcpp(NumericVector & Hi, const int NUMITEMS, IntegerVector & orderHi)
{

	int i, j;         /* for-loop indicators */
	int temp1;        /* A temporary variable */
	double temp;      /* A temporary variable */

	for(i=0; i<NUMITEMS; i++) orderHi[i] = i;

	/* sort the coefficients */
	for(i=0; i<NUMITEMS; i++)
		for(j=0; j<(NUMITEMS-1); j++)
			if(Hi[j] > Hi[j+1])
			{
				temp = Hi[j];
				Hi[j] = Hi[j+1];
				Hi[j+1] = temp;
				temp1 = orderHi[j];
				orderHi[j] = orderHi[j+1];
				orderHi[j+1] = temp1;
			}
}

/*******************************************************************
* sortScales. Descendingly sort the scales by the number of items. *
********************************************************************/

void sortScalesRcpp(IntegerVector& NUMITEMS, const int nclus, IntegerVector& order)
{

	int i, j;
	int temp;      /* A temporary variable */

	for (i=0; i < nclus; i++) order[i] = i;

	/* sort the scales */
	for (i=0; i < nclus; i++)
		for (j=0; j < (nclus-1); j++)
			if(NUMITEMS[j] < NUMITEMS[j+1])
			{
				temp = NUMITEMS[j];
				NUMITEMS[j] = NUMITEMS[j+1];
				NUMITEMS[j+1] = temp;
				temp = order[j];
				order[j] = order[j+1];
				order[j+1] = temp;

			}
}



/*****************************************
* CoefHi. Calculate the Hi-coefficients. *
*****************************************/

void CoefHiRcpp(
	const	IntegerMatrix&	ITEMS,
	const	int				scale,
	const	int				NUMITEMS,
	const	int				NITEM,
	const	NumericMatrix&	VAR,
	const	NumericMatrix&	MAXVAR,
			NumericVector&	Hi) {

	int a, b;                         /* indicators for "items" */
	int i,j;                          /* loop indicators */

	/* sum of variances */
	NumericVector S(NITEM);
	/* sum of maximum variances */
	NumericVector Smax(NITEM);

	/* initialize all Hi coefficients to zero */
	Hi.fill(0);

	/* calculate Hi coefficients */
	for (i=0; i<NUMITEMS; i++)
	{
		a = ITEMS[i + scale*NITEM];
		S[i] = 0.0;
		Smax[i] = 0.0;

		for (j=0; j<NUMITEMS; j++)
		{
			b = ITEMS[j + scale*NITEM];
			if (i != j)
			{
				S[i] += VAR[b+a*NITEM];
				Smax[i] += MAXVAR[b+a*NITEM];
			}
		}
		if (Smax[i] > 0.000001)
			Hi[i] = S[i]/Smax[i];
		else
			Hi[i] = 0;
	}
}



/*******************************************************************
* Criterion2. Investigate whether all items have an Hi coefficient *
* larger than CRITVAL.                                             *
*******************************************************************/

void Criterion2Rcpp(
			int				mem,
			int				scale,
			NumericVector&	Hi,
			IntegerVector&	NUMITEMS,
			int				NITEM,
			int&			NSCALES,
	const	double			CRITVAL,
			IntegerMatrix&	pop,
			IntegerMatrix&	ITEMS,
	const	NumericMatrix&	VAR,
	const	NumericMatrix&	MAXVAR)
{

	/* a vector reflecting the order of the Hi coefficients */
	IntegerVector orderHi(NITEM);

	while (true)
	{

		/* order the Hi coefficients */
		sortHiRcpp(Hi,NUMITEMS[scale],orderHi);

		/* if the smallest Hi coefficients is smaller than CRITVAL,
								then remove items from the scale */
		if(Hi[0] < CRITVAL)
		{

			/* if the scale exist of only two items, delete the whole scale */
			if(NUMITEMS[scale] == 2)
			{
				pop[ ITEMS[scale*NITEM] + mem*NITEM] = 0;
				pop[ ITEMS[scale*NITEM + 1] + mem*NITEM] = 0;
				NUMITEMS[scale] = 0;
				NSCALES--;
				return; // exit while loop and function
			}
			else
			{

				/* delete item with smallest Hi coefficient and adjust matrix ITEMS */
				pop[ ITEMS[orderHi[0]+scale*NITEM] + mem*NITEM ] = 0;
				ITEMS[orderHi[0] + scale*NITEM] = 0;

				if (orderHi[0] < (NUMITEMS[scale]-1))
					for (int i = orderHi[0]; i<(NUMITEMS[scale]-1); i++)
						ITEMS[i + scale*NITEM] = ITEMS[i + scale*NITEM + 1];

				ITEMS[NUMITEMS[scale] + scale*NITEM - 1] = 0;

				NUMITEMS[scale] -= 1;

				/* Recalculate the Hi coefficients for the smaller scale */
				CoefHiRcpp(ITEMS,scale,NUMITEMS[scale],NITEM,VAR,MAXVAR,Hi);

			}
		}
		else
			return; // exit while loop and function
	}
}



/******************************************************************************
* TestHi. Test whether the Hi coefficients are significantly larger than zero *
******************************************************************************/

void TestHiRcpp(
	const	int				mem,
	const	int				scale,
	const	int				NITEM,
			IntegerMatrix&	pop,
			IntegerVector&	NUMITEMS,
			int&			NSCALES,
			IntegerMatrix&	ITEMS,
	const	NumericMatrix&	SijMatrix,
			int				NPERS,
	const	NumericMatrix&	VAR,
	const	double			Zcv)
{

	int a, b;              /* indicators for 'items' */
	int i, j;              /* for-loop indicators */
	int stop=0;            /* while-loop indicators */
	double sumSij, S;

	/* Z-values of Hi-coefficients */
	NumericVector Zi(NITEM);

	/* The order of the z-values */
	IntegerVector orderZi(NITEM);

	while(!stop) {

		/* Calculate the Zi-values */
		for(i=0; i<NUMITEMS[scale]; i++) {
			a = ITEMS[i + scale*NITEM];
			sumSij = 0;
			S = 0;

			for(j=0; j<NUMITEMS[scale]; j++) {
				b = ITEMS[j + scale*NITEM];

				if(i != j) {
					sumSij += SijMatrix[a + b*NITEM];
					S += VAR[a + b*NITEM];
				}
			}

			if(sumSij > 0.000001)
				Zi[i] = S * sqrt(NPERS - 1)/ sqrt(sumSij);
			else
				Zi[i] = 0;

		}

		/* order the Zi values */
		sortHiRcpp(Zi,NUMITEMS[scale],orderZi);

		/* Remove item if smallest Zi < Zcv */
		if(Zi[0] <= Zcv) {
			/* if the scale exist of only two items, delete the whole scale */
			if(NUMITEMS[scale] == 2) {

				pop[ ITEMS[scale*NITEM] + mem*NITEM] = 0;
				pop[ ITEMS[scale*NITEM + 1] + mem*NITEM] = 0;
				NUMITEMS[scale] = 0;
				NSCALES--;
				return;
			}
			else {

				/* delete item with smallest Zi-value and adjust matrix ITEMS */
				pop[ ITEMS[orderZi[0]+scale*NITEM] + mem*NITEM ] = 0;
				ITEMS[orderZi[0] + scale*NITEM] = 0;
				NUMITEMS[scale] -= 1;

				if (orderZi[0] < (NUMITEMS[scale]-1))
					for(i=orderZi[0]; i<(NUMITEMS[scale]-1); i++) {
						ITEMS[i + scale*NITEM] = ITEMS[i + scale*NITEM + 1];
					}
				ITEMS[NUMITEMS[scale] + scale*NITEM - 1] = 0;

			}
		}
		else return;
	}
}



/******************************************************************
* testHij. Investigate whether the first criterion of a Mokken    *
* scale is satisfied, that is, if all Hij coefficients are larger *
* than zero.                                                      *
******************************************************************/

void testHijRcpp(
	const	int				mem,
	const	int				scale,
			IntegerMatrix&	ITEMS,
			IntegerVector&	NUMITEMS,
			int&			NSCALES,
			IntegerMatrix&	pop,
	const	NumericMatrix&	HijMatrix,
	const	int				NITEM)
{
	// TODO: return nscales rather than taking it by reference.
	int i, j, k;                /* for-loop indicators */
	int a, b;                   /* indicators for 'items' */
	int randnum;                /* a randomly drawn number */

	/* check whether all Hij coefficients are larger than 0 */
	for(i=0; i<NUMITEMS[scale]; i++) {
		a = ITEMS[i + scale*NITEM];

		for(j=(i+1); j<NUMITEMS[scale]; j++){
			b = ITEMS[j + scale*NITEM];

			/* is the Hij coefficient smaller than zero? */
			if(HijMatrix[a + b*NITEM] < 0) {
				if(NUMITEMS[scale] == 2){
					/* if the scale exist of only two items, delete the whole scale */
					pop[ ITEMS[scale*NITEM] + mem*NITEM] = 0;
					pop[ ITEMS[scale*NITEM + 1] + mem*NITEM] = 0;
					NUMITEMS[scale] = 0;
					NSCALES--;
				}
				else {
					/* randomly delete one of the two items */
					randnum = R::unif_rand();

					if(randnum < 0.5) {
						pop[mem*NITEM + a] = 0;
						ITEMS[i + scale*NITEM] = 0;
						NUMITEMS[scale] -= 1;

						if(i < (NUMITEMS[scale]))
							for(k=i;k<NUMITEMS[scale];k++)
								ITEMS[k + scale*NITEM] = ITEMS[k + scale*NITEM + 1];
					}
					else {
						pop[mem*NITEM + b] = 0;
						ITEMS[j + scale*NITEM] = 0;
						NUMITEMS[scale] -= 1;

						if(j < (NUMITEMS[scale]))
							for(k=i;k<NUMITEMS[scale];k++)
								ITEMS[k + scale*NITEM] = ITEMS[k + scale*NITEM + 1];

					}
				}
			}
		}
	}
}

/***********************************************
* Initialize. Randomly draw a first population *
***********************************************/
void InitializeRcpp(IntegerVector& population, const int nclus){

	/* Draw random numbers */
	for (int i = 0; i < population.size(); i++) {
		population[i] = int (R::unif_rand() * nclus) + 1;
	}
}


/*******************************************************************
* Evaluation. Evaluate whether the partitionings in the population *
* satisfy the Mokken scale criteria, if not, repair them such that *
* they do.                                                         *
*******************************************************************/

void EvaluateRcpp(
			IntegerMatrix&	pop,
			IntegerMatrix&	newpop,
	const	int				POPSIZE,
	const	int				nclus,
	const	int				NITEM,
	const	int				NPERS,
			NumericVector&	fitness,
	const	NumericMatrix&	VAR,
	const	NumericMatrix&	MAXVAR,
	const	NumericMatrix&	HijMatrix,
//			NumericMatrix&	ZijMatrix,
	const	double			CRITVAL,
	const	NumericMatrix&	SijMatrix,
			double			Zcv)
{

	/* A matrix which indicates which items are in which scale */
	IntegerMatrix ITEMS(nclus, NITEM);

	/* A vector indicating the total number of items in each scale */
	IntegerVector NUMITEMS(nclus);

	/* A vector that orders the scales by the number of items */
	IntegerVector order(nclus);

	/* A vector containing the Hi-coefficients */
	NumericVector Hi(NITEM);

	for (int mem=0;mem<POPSIZE; mem++) {
		fitness[mem] = 0;

		/* Assess the number of scales, the number of items in each scale
	 and which items belong to which scale */
		ScaleNumItemsRcpp(mem, nclus, NUMITEMS, NITEM, pop);
		/* Total number of scales */
		int NSCALES = NumScalesRcpp(nclus, NUMITEMS);
		ScaleItemsRcpp(mem,nclus, NITEM,pop,ITEMS,NUMITEMS);

		/* Ascertain that all item clusters satisfy the Mokken scale criteria */
		for (int i = 0; i < nclus; i++)
		{
			if (NUMITEMS[i] > 1)
			{

				/* Calculate Hi coefficients */
				CoefHiRcpp(ITEMS,i,NUMITEMS[i],NITEM,VAR,MAXVAR,Hi);

				/* Investigate whether all Hi coefficients are larger than CRITVAL */
				Criterion2Rcpp(mem,i,Hi,NUMITEMS,NITEM,NSCALES,CRITVAL,pop,ITEMS,VAR,MAXVAR);

				if(NUMITEMS[i] > 1) {

					/* Test whether Hi coefficients are significantly larger than zero */
					TestHiRcpp(mem, i, NITEM, pop, NUMITEMS, NSCALES, ITEMS, SijMatrix, NPERS, VAR, Zcv);

					if(NUMITEMS[i] > 1) {

						/* test whether the Hij coefficients are larger than zero */
						testHijRcpp(mem, i, ITEMS, NUMITEMS, NSCALES, pop, HijMatrix, NITEM);
					}
				}
			}

		}

		/* Descendingly order the scales by the number of items */
		sortScalesRcpp(NUMITEMS,nclus,order);

		/* Make sure that the largest scale is assigned a 1,
	 the second largest scale a 2, etc.                 */
		for (int j=0; j<NITEM; j++)
		{
			newpop[j + mem*NITEM] = pop[j + mem*NITEM];
			for (int k = 0; k < nclus; k++)
				if(newpop[j + mem*NITEM] == (order[k]+1))
					pop[j + mem*NITEM] = k+1;
		}

		/* Calculate the fitness value */
		for (int k = 0; k < nclus; k++)
			fitness[mem] += (pow(NITEM,-(k+1))*NUMITEMS[k]);

	}
}



/*************************************************************
* Selection. Select a new population from the old population *
* based on the fitness values of the members of the old      *
* population.                                                *
*************************************************************/

void SelectionRcpp(
	const	int				POPSIZE,
	const	int				NITEM,
			NumericVector&	fitness,
			IntegerMatrix&	pop,
			IntegerMatrix&	newpop)
{
	int mem, i, j, k;
	double TotalFitness = 0.0;
	double p;

	/* calculate total fitness of the population */
	for(mem=0; mem<POPSIZE; mem++) TotalFitness += fitness[mem];

	/* calculate relative fitness */
	for(mem=0; mem<POPSIZE; mem++)
		fitness[mem + POPSIZE + 2] = fitness[mem]/TotalFitness;

	fitness[2*(POPSIZE+2)] = fitness[POPSIZE + 2];

	/* calculate cumulative fitness */
	for(mem=1; mem<POPSIZE; mem++)
		fitness[mem + 2*(POPSIZE + 2)] =
				fitness[(mem-1) + 2*(POPSIZE+2)] + fitness[mem + POPSIZE + 2];

	/* finally, select the new population using cumulative fitness */
	for(i=0; i<POPSIZE; i++) {

		/* Draw a random number between 0 and 1 */
		p = R::unif_rand();

		if(p < fitness[2*(POPSIZE+2)])
			for(j=0; j<NITEM; j++)
				newpop[j + i*NITEM] = pop[j];
		else
			for(mem=0; mem<POPSIZE; mem++)
				if(p >= fitness[mem + 2*(POPSIZE + 2)]
						&& p < fitness[(mem+1)+ 2*(POPSIZE + 2)])
					for(k=0; k<NITEM; k++)
						newpop[k + i*NITEM] = pop[k + (mem+1)*NITEM];
	}

	/* once a new population is created, copy it back to variable 'population' */
	for(i=0; i<POPSIZE; i++)
		for(j=0; j<NITEM; j++)
			pop[j + i*NITEM] = newpop[j + i*NITEM];


}

/************************************************************
* Crossover selection: Select two members that take part in *
* the crossover. Then, for each pair of member select two   *
* points in the vector indicating which elements are        *
* exchanged.                                                *
************************************************************/

void CrossoverRcpp(
	const	int				POPSIZE,
	const	int				NITEM,
	const	double			PXOVER,
			IntegerMatrix&	pop)
{

	int i, j;                 /* for-loop indicators */
	int point1, point2;       /* random points of the a vector */

	/* a vector containing the members that are selected for crossover */
	IntegerVector members(POPSIZE);

	/* a vector of randomly drawn numbers */
	NumericVector CrossOverMatrix(POPSIZE);

	int count=0;               /* total number of members selected */
	int temp;                  /* a temporary variable */

	/* Draw random number */
	for(i=0; i<POPSIZE; i++)
		CrossOverMatrix[i] = R::unif_rand();

	/* select the members */
	for(i=0; i<POPSIZE; i++)
		if(CrossOverMatrix[i] < PXOVER) {
			members[count] = i;
			count++;
		}

	/* ascertain that the number of members is even */
	if(count%2 == 1) count = count - 1;

	/* for all pairs of members cross over all elements between point 1 and point2 */
	for(i=0; i<count; i+=2) {
		point1 = (int) (R::unif_rand() * NITEM);
		point2 = (int) (R::unif_rand() * NITEM);

		if(point1 < point2)
		{
			for(j=point1; j<=point2; j++) {
				temp = pop[j + NITEM*members[i]];
				pop[j + NITEM*members[i]] = pop[j + NITEM*members[i+1]];
				pop[j + NITEM*members[i+1]] = temp;
			}
		}
		else if(point1 > point2)
		{
			for(j=0; j<=point2; j++) {
				temp = pop[ j + NITEM*members[i] ];
				pop[ j + NITEM*members[i]] = pop[ j + NITEM*members[i+1]];
				pop[j + NITEM*members[i+1]] = temp;
			}
			for(j=point1; j<NITEM; j++) {
				temp = pop[j + NITEM*members[i]];
				pop[j + NITEM*members[i]] = pop[j + NITEM*members[i+1]];
				pop[j + NITEM*members[i+1]] = temp;
			}
		}
		else // point1 == point2
		{
			temp = pop[point1 + NITEM*members[i]];
			pop[point1 + NITEM*members[i]] = pop[point1 + NITEM*members[i+1]];
			pop[point1 + NITEM*members[i+1]] = temp;
		}

	}
}

/*************************************************************
* Mutation. Each single element of the population has a      *
* probability of PMUTATION to change to another value.       *
* This is equivalent to assigning that item to another scale *
*************************************************************/

void MutationRcpp(
			IntegerMatrix&	pop,
	const	int				NITEM,
	const	int				POPSIZE,
	const	int				nclus,
	const	double			PMUTATION)
{

	int i, j, k;              /* for-loop indicators */
	int check, temp, stop;    /* while-loop indicators */
	int nscales;              /* total number of scales */

	/* A matrix consisting of random numbers */
	NumericMatrix MutationMatrix(POPSIZE, NITEM);

	/* Draw random numbers */
	for(i=0; i<POPSIZE; i++)
		for(j=0; j<NITEM; j++)
			MutationMatrix[j + i*NITEM] = R::unif_rand();

	/* if MutationMatrix[j + i*NITEM] < PMUTATION, then change
		population[j + i*NITEM] to another number */
	for(i=0; i<POPSIZE; i++) {
		nscales = 0;

		for(j=1; j<nclus; j++) {
			check = 0;
			stop = 0;
			k = 0;

			/* is there any item in scale j? */
			while(!stop) {

				if(pop[k + i*NITEM] == j) {
					check = 1;
					stop = 1;
				}

				if(k == (NITEM-1)) stop = 1;
				k++;
			}

			if(check == 1) nscales ++;
		}
		for(j=0; j<NITEM; j++) {
			if(MutationMatrix[j + i*NITEM] < PMUTATION) {
				temp =  pop[j + i*NITEM];

				while(temp == pop[j + i*NITEM]) {
					pop[j + i*NITEM]  = (int) (R::unif_rand()*(nscales+1)) + 1;
				}

				if(pop[j + i*NITEM] == (nscales + 1) ) nscales++;
			}
		}
	}
}


/*************************************************************
* Keep the best function: This function keeps track of the   *
* best member of the population. Note that the last entry in *
* the array Population holds a copy of the best invidual     *
*************************************************************/

int KeepTheBestRcpp(
			IntegerMatrix&	pop,
			NumericVector&	fitness,
	const	int				NITEM,
	const	int				POPSIZE,
			IntegerVector&	generation,
			int				itercount)
{

	int mem;
	int i;
	int CurrentBest=0;

	/* Select the first member of the population */
	fitness[POPSIZE+1] = fitness[0];

	/* Compare fitness to subsequent members */
	for(mem=1; mem<POPSIZE; mem++)
		if(fitness[mem] > fitness[CurrentBest]) {
			CurrentBest = mem;
			fitness[POPSIZE+1] = fitness[mem];
		}

	/* Once the best member in the population is found, copy the genes */
	for(i=0; i<NITEM; i++)
		pop[(POPSIZE+1)*NITEM + i] = pop[CurrentBest*NITEM + i];

	/* Compare the best member with the stored best partitioning */
	if( fitness[POPSIZE+1] > fitness[POPSIZE]) {
		fitness[POPSIZE] = fitness[POPSIZE+1];

		for(i=0; i<NITEM; i++)
			pop[POPSIZE*NITEM + i] = pop[(POPSIZE+1)*NITEM + i];

		generation[0]=0;
		itercount = 0;
	}
	return itercount;

}



/*************************************************************
* Elitist model.                                             *
* Assess which member of the new population is the best and  *
* which member of the new population is the worst.           *
*                                                            *
* If best individual from the new population is worse than   *
* the best individual from the previous populations, replace *
* the worst individual from the current population with the  *
* best one from the previous generation.                     *
*************************************************************/

void ElitistRcpp(
	const	int				POPSIZE,
	const	int				NITEM,
			NumericVector&	fitness,
			IntegerVector&	pop)
{

	int i;
	double best = 0.0;     /* fitness value of the best member */
	double worst = 1.0;    /* fitness value of the worst member */
	int WorstMember = 0;   /* which member had the smallest fitness */

	/* search for the best and worst member */
	for(i=0; i<POPSIZE; i++) {
		if(fitness[i] > best) {
			best = fitness[i];
		}
		if(fitness[i] < worst) {
			worst = fitness[i+1];
			WorstMember = i + 1;
		}
	}

	/* replace worst member by best member */
	if(best < fitness[POPSIZE]){
		for(i=0; i<NITEM; i++)
			pop[i + WorstMember*NITEM] = pop[i + NITEM*POPSIZE];
		fitness[WorstMember] = fitness[POPSIZE];
	}
}



/*******************************************************************************
* Function 'GeneticAlgorithm()' is the main body of the algorithm. All input   *
* variables are arguments that are provided by the function in R.              *
* - popsize is the number of member in a population.                           *
* - nitem is the number of items in the data set.                              *
* - npers is the number of respondents in the data set.                        *
* - maxgens is the number of iterations in which no change of the stored best  *
*   partitioning is found until convergence is reached.                        *
* - ALPHA is the level of significance value.                                  *
* - VAR is the variance-covariance matrix.                                     *
* - MAXVAR is the maximum possible variance-covariance matrix given the        *
*   marginal distribution.                                                     *
* - SijMatrix is a matrix consisting of the outer product of the variances     *
*   with themselves.                                                           *
* - critval is the critical value in the second criterion of a Mokken scale.   *
* - pxover is the probability that a member takes part in a crossover.         *
* - pmutation is the probability that an element undergoes a mutation.         *
* - Output is a vector of zero's which will be converted into the final        *
*    partitioning.                                                             *
*******************************************************************************/

void GeneticAlgorithmRcpp(
	const	int				POPSIZE,
	const	int				NPERS,
	const	int				MAXGENS,
	const	double			PXOVER,
	const	double			PMUTATION,
	const	double			critval,
	const	double			alpha,
	const	int				NITEM,
	const	NumericMatrix&	VAR,
	const	NumericMatrix&	MAXVAR,
	const	NumericMatrix&	SijMatrix,
	const	NumericMatrix&	HijMatrix,
	// these three change across iterations
			int&			itercount,
			IntegerMatrix&	population,
			NumericVector&	fitness,
			IntegerMatrix&	newpopulation
)
{
	// TODO: check if generation can just be an int, passed by reference
	int count = 0;				/* indicator for a while loop */
	double TotalFitness = 0.0;	/* the total fitness of the whole population */

	/* the critical Z-value used in tests */
	const double probability = 1 - alpha;

	/* the critical Z-value used in tests */
	const double Zcv = R::qnorm(probability, 0, 1, 1, 0);

	/* maximum number of scales */
	const int nclus = NITEM / 2; // integer division is on purpose.

	IntegerVector generation(1);
	generation[0] = 0;

	// zero these elements in between iterations
	for (int i = POPSIZE + 2; i < fitness.size(); i++)
		fitness[i] = 0.0;

	newpopulation.fill(0);
	int popSum = Rcpp::sum(population);

	if (popSum == 0)
	{
		while(!TotalFitness)
		{
			count++;

			/*  Draw a first random population */
			InitializeRcpp(population, nclus);

			/*Evaluate initial population */
			EvaluateRcpp(population, newpopulation, POPSIZE, nclus, NITEM, NPERS, fitness, VAR, MAXVAR,
						 HijMatrix, critval, SijMatrix, Zcv);

			for(int i=0; i<POPSIZE; i++) TotalFitness += fitness[i];

			/* If no valid populations was found in MAXGENS iterations, then exit */
			if(count == MAXGENS)
			{
				Rprintf("No partitioning was found in %d populations\n",MAXGENS);
				generation[0] = MAXGENS;
				TotalFitness = 1.0;
			}
		}

		/* Store the best partitioning */
		itercount = KeepTheBestRcpp(population,fitness,NITEM,POPSIZE,generation,itercount);
	}

	/* An iterative process to find the final partitioning */
	while(generation[0]<MAXGENS) {

		TotalFitness = 0.0;
		generation[0]++;



		/* select a new population */
		SelectionRcpp(POPSIZE,NITEM,fitness,population,newpopulation);

		/* crossovers between members */
		CrossoverRcpp(POPSIZE, NITEM, PXOVER, population);

		/* mutation of genes of members*/
		MutationRcpp(population, NITEM, POPSIZE, nclus, PMUTATION);

		/* evaluate the new population */
		EvaluateRcpp(population, newpopulation, POPSIZE, nclus, NITEM, NPERS, fitness,
					 VAR, MAXVAR, HijMatrix, critval, SijMatrix, Zcv);

		/* store the best partitioning */
		itercount = KeepTheBestRcpp(population,fitness,NITEM,POPSIZE,generation,itercount);

		/* ascertain that best partitioning always stays in the population */
		ElitistRcpp(POPSIZE, NITEM, fitness, population);


		if(fitness[POPSIZE] == 1)   generation[0] = MAXGENS;
	}
	itercount++;

}

// [[Rcpp::export]]
IntegerMatrix runGeneticAlgorithm(
		const int POPSIZE, const int NPERS, const int MAXGENS,
		const double PXOVER, const double PMUTATION, const double critval,
		const double alpha, const int NITEM, const int ITER,
		const NumericMatrix& VAR, const NumericMatrix& MAXVAR,
		const NumericMatrix& SijMatrix)
{

	// initialize all vectors here
	int itercount = 0;

	IntegerMatrix population(NITEM, POPSIZE + 2);
	population.fill(0);
	NumericVector fitness((POPSIZE + 2) * 3);
	fitness.fill(0);

	NumericMatrix HijMatrix(NITEM, NITEM);
	IntegerMatrix newpopulation(POPSIZE, NITEM);

	for(int i = 0; i<NITEM; i++)
	{
		for(int j = 0; j<NITEM; j++)
		{
			if (VAR(j, i) > 0.0000001)
				HijMatrix(j, i) = VAR(j, i) / MAXVAR(j, i);
		}
	}

	const int iterCheck = int(ceil(double(MAXGENS) / double(ITER)));
	do
	{
		GeneticAlgorithmRcpp(POPSIZE, NPERS, MAXGENS,PXOVER, PMUTATION, critval,
							alpha, NITEM, VAR, MAXVAR, SijMatrix, HijMatrix,
							// these three change across iterations
							itercount, population, fitness,
							// preinitialized memory
							newpopulation);
	}
	while (itercount != iterCheck);

	return population;
}

