#include <stdio.h>
#include <stdlib.h>
#include "FileHandling.h"
#include <math.h>
#include <string.h>
#include <time.h>


int main(int argc, char *argv[])
{

  int dim = 0;          // Dimensions of Elements
  int n = 0;            // n Elements
  int i,j,d;            // i counter for n, j counter for Neighborhood Elements, d counter for dimensions
  int kPoints;          // Minimum points for a cluster to be created
  int h,k;              // h counter for every other Element, k counter for files
  char filename[40];    // The given file
  clock_t start, end;   // Variables for counting execution time
  int counter;          // Counter to count all the elements of a Neighborhood
  double reachDist = 0;
  double max;
  int choise;// Choise holder for switch
  char c;

/* -------------------------------------------------------------------------- */
    // Reading Dataset

   FILE* Dataset;

      printf("\n Give the DataSet file:");
      scanf("%s", filename);

    Dataset = fopen(filename, "r");
    if (!Dataset)
    {
        printf("\n There is something wrong with the Dataset file! \n\n");
        return -1;
    }

/* -------------------------------------------------------------------------- */


    dim = getColumns(Dataset); //Getting the dimensions of each Element of the Dataset
    rewind(Dataset);


    n = getRows(Dataset);       // Getting the Elements of the Dataset
    rewind(Dataset);

    printf("\n Elements:%d \n", n-1);
    printf("\n Dimensions:%d \n", dim);
           n--;
    printf("Give the amount of K Points: ");   // kPoints
    scanf("%d",&kPoints );

/* -------------------------------------------------------------------------- */
        // All the necessary memory allocation

        double **X;   // Array of Elements
        X =(double **)calloc(n, sizeof(double *));
        for (d = 0; d < n; d++)
        X[d] = (double *)calloc(dim, sizeof(double));

        double *kDistance;   // Array for holding k-Distances for the First Core point of a cluster
        kDistance =(double *) calloc(n, sizeof(double));
        for(i = 0; i < n; i++)
        kDistance[i] = 0;

        double **distance;   // Array for holding Distances for each ELement with each other Element
        distance = (double **)calloc(n, sizeof(double *));
        for(i = 0; i < n; i++)
        distance[i] =(double *)calloc(n,sizeof(double));

        long int **Neighborhood; //Array for holding which element belong to which Neighborhood
        Neighborhood =(long int **) calloc(n, sizeof(long int *));
        for(i = 0; i < n; i++)
        Neighborhood[i] = (long int *)calloc(n,sizeof(long int));


        double *lrd; //Array for holding  Local Reachability Density for each element
        lrd =(double *)calloc(n, sizeof(double));

        int *NeighborhoodSize; //Array for holding the size of each Neighborhood
        NeighborhoodSize = (int *)calloc(n,sizeof(int));

        double *reachDistSum; //Array for holding the sum of reachability Distances for each element
        reachDistSum = (double *)calloc(n,sizeof(double));

        double *LOF;//Array for holding LOF value for each Element
        LOF = (double *)calloc(n,sizeof(double));

        double *NeighborhoodLrdSum; //Array for holding the sum of each Neighborhood lrd
        NeighborhoodLrdSum =(double *) calloc(n,sizeof(double));

        double **OrderedList;
        OrderedList = (double **)calloc(n,sizeof(double *));
        for(d = 0; d < n; d++)
        OrderedList[d] =(double *)calloc(n,sizeof(double));

        double **tmp2;
        tmp2 =(double **) calloc(n,sizeof(double *));
        for(d = 0; d < n; d++)
        tmp2[d] = (double *)calloc(n,sizeof(double));

        double *tempDistance;   // Array for holding Distances for each ELement with each other Element
        tempDistance = (double *)calloc(n, sizeof(double));




        for(i = 0; i < n; i++)
        NeighborhoodLrdSum[i] = 0;

        for(i = 0; i < n; i++)
        reachDistSum[i] = 0;

        for(i = 0; i < n; i++)
        NeighborhoodSize[i] = 0;

        for(i = 0; i < n; i++)
        {
          lrd[i] = 0;
        }



/* -------------------------------------------------------------------------- */
                 // Passing elements to Array X[n][dim]



                    X = getData(Dataset,n,dim,X);

                 for(i = 0; i < n; i++)
                 {
                   for(d = 0; d < dim; d++)
                   {
                     OrderedList[i][d] = X[i][d];
                   }
                 }
           fclose(Dataset);
  start = clock();
/* -------------------------------------------------------------------------- */
/* ---------------------------------LOF-------------------------------------- */
/* -------------------------------------------------------------------------- */
                       //STEP 1

         // Finding the k-Distance of each element in the dataset
         for(i = 0; i < n; i++)
         {
           for(h = 0; h < n; h++)
           {
             if(h != i)
             {
               /* Calculating the distance of each element with i element
               and then we store it to an array so that we can order it later,
               we set the distance of i to max value so that we can avoid it at
               the ordering */
               distance[i][i] = 99999;
               distance[i][h] = 0;
               tempDistance[i] = 99999;
               for(d = 0; d < dim; d++)
               distance[i][h] += (X[h][d] - X[i][d])*(X[h][d] - X[i][d]);

               distance[i][h] = sqrt(distance[i][h]);
               tempDistance[h] = distance[i][h];
             }
           }
            /* Ordering the array holding the distances from min to max,
            (distance[i] will be located at the distance[max]) */
            qS(tempDistance,0,n-1);
            /* Using the ordered array to extract the k-distance of i element,
            DIST_k(i) */
            kDistance[i] = tempDistance[kPoints-1];

            // printf("%d : K-Distance: %lf\n",i,kDistance[i] );

         }
/* -------------------------------------------------------------------------- */

                   //STEP 2

     //Finding the Neighborhood of each element
     for(i = 0; i < n; i++)
     {
       for(h = 0; h < n; h++)
       {
         if(h != i)
         {
           /* Calculating the distance of each element with i element
           and then we store it to an array so that we can determine if that
           element belongs to the Neighborhood of i element */
           Neighborhood[i][h] = 0;


           /* Comparing the distance of each potential Neighborhood element
           with the k-distance of i element, if found less, set
           Neighborhood[master][currentElement] to 1, that means that it belongs
           to i element's Neighborhood */
           if(distance[i][h] <= kDistance[i])
           {
             Neighborhood[i][h] = 1;
           }
         }
       }
     }

/* -------------------------------------------------------------------------- */
                                //STEP 3

      //Getting the Neighborhood size of each element
      for(i = 0; i < n; i++)
      {
       counter = 0;
       for(h = 0; h < n; h++)
       {
         if(h != i)
         {

          if(Neighborhood[i][h] != 0)
          {
            counter++;
          }
         }
       }
         NeighborhoodSize[i] = counter;

         // printf("%d Neighborhood: %d\n",i,NeighborhoodSize[i] );
       }
/* -------------------------------------------------------------------------- */

               //STEP 4
      //Finding the Sum of reachDist_k(h <- i)
      for(i = 0; i < n; i++)
      {

       for(h = 0; h < n; h++)
       {
          if(Neighborhood[i][h] != 0)
          {
            /*Calculating the reach dist of i in respect to h so that we will be
            able to compare it with the k-distance of h */
            reachDist = 0;
            for(d = 0; d < dim; d++)
            reachDist += (X[i][d] - X[h][d])*(X[i][d] - X[h][d]);

            reachDist = sqrt(reachDist);

       /*Comparing kDistance of h with the reach dist of i in respect to h,
       we choose the Max between them, max(kDistance(h),reachDist_k(h <- i))*/
            if(kDistance[h] > reachDist)
            {
              max = kDistance[h];
            }else
            {
              max = reachDist;
            }

            //Adding the max to the total Sum
            reachDistSum[i] += max;

          }
       }

      }

/* -------------------------------------------------------------------------- */
                                //STEP 4
            /* Calculating the Local Reachability Density(lrd) of each element,
            lrd_k(i) = NeighborhoodSize(i)/reachDistSum(i) */

            for(i = 0; i < n; i++)
            {
                lrd[i] = NeighborhoodSize[i]/(reachDistSum[i]);

            }



/* -------------------------------------------------------------------------- */
                  //STEP 5
          /* Calculating the sum of a Neighborhood element lrd divided by the
          master's element lrd for each i.
             Sum(lrd[h]/lrd[i]) for each Neighborhood element */
            for(i = 0; i < n; i++)
            {
              for(h = 0; h < n; h++)
              {
                if(Neighborhood[i][h] != 0)
                {
                  NeighborhoodLrdSum[i] += lrd[h]/lrd[i];
                }
              }

             }
/* -------------------------------------------------------------------------- */
                       //STEP 6
           //Calculating the LOF for each element
            for(i = 0; i < n; i++)
            {
             LOF[i] = (NeighborhoodLrdSum[i]/NeighborhoodSize[i]);

            }

/* -------------------------------------------------------------------------- */
              //STEP 7
            //Ordering the data according to LOF from max to min
            for(i = 0; i < n; i++)
            {
              for(h = 0; h < n; h++)
              {
                if(LOF[i] >= LOF[h])
                {
                  double tmp = LOF[i];
                  LOF[i] = LOF[h];
                  LOF[h] = tmp;

                  for(d = 0; d < dim; d++)
                  {
                    tmp2[0][d] = OrderedList[i][d];
                    OrderedList[i][d] = OrderedList[h][d];
                    OrderedList[h][d] = tmp2[0][d];
                  }
                }
              }
            }

/* -------------------------------------------------------------------------- */
/* ---------------------------------END-------------------------------------- */
/* -------------------------------------------------------------------------- */
end = clock();
/*----------------------------------------------------------------------------*/
double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("\n Time of Algorithm Execution: %lf \n\n",total_time);



FILE* OrderFile;
OrderFile = fopen("OrderFile.txt","w");



for(i = 0; i < n; i++)
{
  fprintf(OrderFile, "LOF : %lf --> :",LOF[i]);
  for(d = 0; d < dim; d++)
  {
    fprintf(OrderFile,"%lf ",OrderedList[i][d] );
  }
  fprintf(OrderFile, "\n");
}

fclose(OrderFile);

/*----------------------------------------------------------------------------*/
    // FREE EVERYTHING

    for(i = 0; i < n; i++)
    free(X[i]);
    free(X);
    free(kDistance);
    for(i = 0; i < n; i++)
    free(distance[i]);
    free(distance);
    for(i = 0; i < n; i++)
    free(Neighborhood[i]);
    free(Neighborhood);
    free(lrd);
    free(NeighborhoodSize);
    free(reachDistSum);
    free(NeighborhoodLrdSum);
    free(LOF);
    for(i = 0; i < n; i++)
    free(tmp2[i]);
    free(tmp2);
    for(i = 0; i < n; i++)
    free(OrderedList[i]);
    free(OrderedList);
    free(tempDistance);
    return 0;

}
