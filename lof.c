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

        float *X;  // Array of Elements
        X = (float*)calloc(n*dim, sizeof(float));

        float *kDistance;   // Array for holding k-Distances for the First Core point of a cluster
        kDistance = (float*)calloc(n, sizeof(float));


        float *distance;   // Array for holding Distances for each ELement with each other Element
         distance = (float*)calloc(n*n, sizeof(float));

         unsigned int *Neighborhood; //Array for holding which element belong to which Neighborhood
         Neighborhood =(int*) calloc(n*n, sizeof(int));


        float *lrd; //Array for holding  Local Reachability Density for each element
        lrd =(float *)calloc(n, sizeof(float));

        unsigned int *NeighborhoodSize; //Array for holding the size of each Neighborhood
        NeighborhoodSize = (int*)calloc(n,sizeof(int));

        float *reachDistSum; //Array for holding the sum of reachability Distances for each element
        reachDistSum = (float *)calloc(n,sizeof(float));

        float *LOF;//Array for holding LOF value for each Element
        LOF = (float *)calloc(n,sizeof(float));

        float *NeighborhoodLrdSum; //Array for holding the sum of each Neighborhood lrd
        NeighborhoodLrdSum =(float *) calloc(n,sizeof(float));

        float *OrderedList;
        OrderedList =(float*) calloc(n*dim,sizeof(float));

        float *tmp2;
        tmp2 =(float*) calloc(n*dim,sizeof(float*));

        float *tempDistance;   // Array for holding Distances for each ELement with each other Element
        tempDistance =(float*) calloc(n, sizeof(float));






          for(i = n; i--;)
          {
          kDistance[i] = 0;
        NeighborhoodLrdSum[i] = 0;
        reachDistSum[i] = 0;
        NeighborhoodSize[i] = 0;
        lrd[i] = 0;
           }







/* -------------------------------------------------------------------------- */
                 // Passing elements to Array X[n][dim]



                    X = getData(Dataset,n,dim,X);

                 for(i = n; i--;)
                 {
                   for(d = dim; d--;)
                   {
                     OrderedList[i*dim + d] = X[i*dim + d];
                   }
                 }
           fclose(Dataset);
  start = clock();
/* -------------------------------------------------------------------------- */
/* ---------------------------------LOF-------------------------------------- */
/* -------------------------------------------------------------------------- */
                       //STEP 1

         // Finding the k-Distance of each element in the dataset
         for(i = n; i--;)
         {

           distance[i*n + i] =  9999;
           tempDistance[i] = 99999;
           for(h = n; h--;)
           {

             if(h != i)
             {
            distance[i*n + h] =  0;
               /* Calculating the distance of each element with i element
               and then we store it to an array so that we can order it later,
               we set the distance of i to max value so that we can avoid it at
               the ordering */

               for(d = dim; d--;)
               distance[i*n + h] += (X[h*dim + d] - X[i*dim + d])*(X[h*dim + d] - X[i*dim + d]);

               distance[i*n + h] = sqrt(distance[i*n +h]);
               tempDistance[h] = distance[i*n + h];

             }
           }
            /* Ordering the array holding the distances from min to max,
            (distance[i] will be located at the distance[max]) */
              qSort(tempDistance,n);
            /* Using the ordered array to extract the k-distance of i element,
            DIST_k(i) */
            kDistance[i] = tempDistance[kPoints - 1];

            // printf("%d : K-Distance: %lf\n",i,kDistance[i] );

         }
/* -------------------------------------------------------------------------- */

                   //STEP 2

     //Finding the Neighborhood of each element
     for(i = n; i--;)
     {
       Neighborhood[i*n + i] = 0;
       for(h = n; h--;)
       {
         if(h != i)
         {
           /* Calculating the distance of each element with i element
           and then we store it to an array so that we can determine if that
           element belongs to the Neighborhood of i element */
          Neighborhood[i*n + h] = 0;


           /* Comparing the distance of each potential Neighborhood element
           with the k-distance of i element, if found less, set
           Neighborhood[master][currentElement] to 1, that means that it belongs
           to i element's Neighborhood */
           if(distance[i*n + h] <= kDistance[i])
           {
            Neighborhood[i*n + h] = 1;
           }
         }
       }
     }

/* -------------------------------------------------------------------------- */
                                //STEP 3

      //Getting the Neighborhood size of each element
      for(i = n; i--;)
      {
       unsigned int counter = 0;
        for(h = n; h--;)
        {
          if(Neighborhood[i*n + h] != 0)
          {
            counter++;
          }
        }
        NeighborhoodSize[i] = counter;

      }
/* -------------------------------------------------------------------------- */

               //STEP 4
      //Finding the Sum of reachDist_k(h <- i)
      for(i = n; i--;)
      {

       for(h = n; h--;)
       {
          if(Neighborhood[i*n + h] != 0)
          {
            /*Calculating the reach dist of i in respect to h so that we will be
            able to compare it with the k-distance of h */
            reachDist = 0;
            for(d = dim; d--;)
            reachDist += (X[i*dim + d] - X[h*dim + d])*(X[i*dim + d] - X[h*dim + d]);

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

            for(i = n; i--;)
            {
                lrd[i] = NeighborhoodSize[i]/(reachDistSum[i]);
            }



/* -------------------------------------------------------------------------- */
                  //STEP 5
          /* Calculating the sum of a Neighborhood element lrd divided by the
          master's element lrd for each i.
             Sum(lrd[h]/lrd[i]) for each Neighborhood element */
            for(i = n; i--;)
            {
              for(h = n; h--;)
              {
                if(Neighborhood[i*n + h] != 0)
                {
                  NeighborhoodLrdSum[i] += lrd[h]/lrd[i];
                }
              }

             }
/* -------------------------------------------------------------------------- */
                       //STEP 6
           //Calculating the LOF for each element
            for(i = n; i--;)
            {
             LOF[i] = (NeighborhoodLrdSum[i]/NeighborhoodSize[i]);

            }

/* -------------------------------------------------------------------------- */
              //STEP 7
            //Ordering the data according to LOF from max to min
            for(i = n; i--; )
            {
              for(h = n; h--;)
              {
                if(LOF[i] >= LOF[h])
                {
                  double tmp = LOF[i];
                  LOF[i] = LOF[h];
                  LOF[h] = tmp;

                  for(d = dim; d--;)
                  {
                    tmp2[i*dim + d] = OrderedList[i*dim + d];
                    OrderedList[i*dim + d] = OrderedList[h*dim + d];
                    OrderedList[h*dim + d] = tmp2[i*dim + d];
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



for(i = n; i--;)
{
  fprintf(OrderFile, "LOF : %lf --> :",LOF[i]);
  for(d = dim; d--;)
  {
    fprintf(OrderFile,"%lf ",OrderedList[i*dim + d] );
  }
  fprintf(OrderFile, "\n");
}

fclose(OrderFile);

/*----------------------------------------------------------------------------*/
    // FREE EVERYTHING

    free(X);
    free(kDistance);
    free(distance);
    free(Neighborhood);
    free(lrd);
    free(NeighborhoodSize);
    free(reachDistSum);
    free(NeighborhoodLrdSum);
    free(LOF);
    free(tmp2);
    free(OrderedList);
    free(tempDistance);
    return 0;

}
