#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>

main()
{
	srand(time(NULL));
	FILE *fp = fopen("data", "w");
	int i = 0;
	int m = 3*1024;
	int n = 3*1024;
	float a = 0, b=0; 

	for(i=0; i < m*n; i++)
	{
		a = rand()%36;
		b = rand()%2;
		if(b == 0)
			fprintf(fp,"%f\n", a);
		else
			fprintf(fp,"%f\n", -1*a);
			
	}		
	fclose(fp);
}
