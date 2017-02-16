/****************************************************************/
/* This is a temperature conversion program written in C        */
/* This program generates and displays a table of Fahrenheit and */
/* centigrade temperatures, and lists the freezing and boiling  */
/* of water.                                                    */
/* **************************************************************/

#include "stdio.h" /* access the stdio library */
#include <stdlib.h>
/* define 'constants' */

#define FREEZING 0
#define BOILING 100

/* function prototype */
void hotcold(int centigrade);

int main(void) /* start of main function */
{
   /* local declarations within main */
   int i;          /* loop control variable */
   int Fahrenheit;  /* temperature in Fahrenheit degrees */
   int centigrade; /* temperature in centigrade degrees */

   printf("Centigrade to Fahrenheit temperature table\n\n");
   for (i = -2; i <= 12; i++) {
      centigrade = 10 * i;
      Fahrenheit = 32 + (centigrade * 9)/5;
      /* print the results */
      printf("C =%4d F =%4d ",centigrade,Fahrenheit);
      hotcold(centigrade); /* print messages at 0 and 100 centigrade */
   } /* end of for loop */
 /*system("pause");*/
   printf("Press enter to continue...\n");
   getchar();
   return(0);
} /* end of main function */

void hotcold(centigrade)
{
/* function to display messages at 0 and 100 centigrade */
   if (centigrade == FREEZING)
      printf(" <= Freezing Point of Water - Brrrr!\n");
   else
      if (centigrade == BOILING)
         printf(" <= Boiling Point of Water - Owwww!\n");
      else
         printf("\n");
} /* end of function hotcold */

/* end of program */
