#include <stdlib.h>
#include <stdio.h>
#define COMMANDS 2

//takes arguments from array structures and outputs a GNUPLOT of the data with connecting lines.
void pfun(double x[], double y[], int points)
{
    char * commandsForGnuplot[] = {"set title \"Specific Heat versus Temperature\"", "plot 'data.temp' w l"};
    FILE * temp = fopen("data.temp", "w"); // write coordinates here.
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    int q;
	for(q = 0; q < points ;q++)
	{
		fprintf(temp, "%lf %lf \n", x[q], y[q]);
	}

 
    for (i=0; i < 2; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
}

