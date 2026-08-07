#include <stdio.h>
#include <stdlib.h>
#include <string.h>


class Variables
{	
private:
	char **var;		// vector pointer
	int sizet;

public:
	Variables(char *str)
	{
		int i = 0;
		int o = 0;
		char car;
		char *temp2 = (char *)malloc(sizeof(char)*100);
		sizet =0;
		// counting the number of variables
		while(str[i] != 0)
		{
			if(str[i] == '|') sizet++;
			i++;
		}
		// memory allocation for string vector
		var = (char **)malloc(sizeof(char *)*sizet);
		
		// getting the variables
		int count = 0;
		i = 0;
		while((str[i] != 0))
		{
			o =0; 
			while((str[i] != 0) && (str[i] != '#'))
			{
				
				if(str[i] == '|') i++;
				if(str[i] != 0)
				{	
					car = str[i];
					temp2[o] = car;
					o++;
					i++;
				}
			}
			i++;
			temp2[o] = 0;
			// string memory allocation
			var[count]  = (char *)malloc(sizeof(char)*strlen(temp2));
			// string copy
			strcpy(var[count], temp2);
			count++;
		}

		free(temp2);
	}

	int getQuantity()
	{
		return sizet;
	}
	
 	char* getVar(int ind)
	{
		if((ind >= 0) && (ind < sizet))
			return var[ind];
		else return NULL;
	}
	
	~Variables()
	{
		int c;
		for(c=0; c < sizet; c++)
		if(var[c])
				free(var[c]);
		if(var)
			free(var);
	}	
};

