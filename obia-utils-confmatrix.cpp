/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * confmatrix.cpp: Dip 2 Confusion Matrix   		                      *
 * Version: 0.8.2 build 20 Nightly                                        *
 * Last revised: 15/06/2006                                               *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * DIP 2 by Angelos Tzotsos (Kalxas), Chris Iossifides (Chiossif)         *
 * MSEG algorithm by Angelos Tzotsos (Kalxas)                             *
 * GCpp Lab.                                             January 2006 	  *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * */

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <math.h>


using namespace std;

int main(int argc, char *argv[])
{
    float ac;
    int count=0;
	int success=0;
	int **cm;
	int classes = atoi(argv[3]);
	int lines = atoi(argv[4]);
	int columns = atoi(argv[5]);

    if(argc != 6)   //TODO: Add to Error Handling
    {
        cout << "Error in arguments:" << endl << endl
        << "Usage : confmatrix <input class raster> <input truth image> <class number> <image lines> <image columns>"
        << endl;
        exit(1);
    }

    cm = (int **)malloc(classes * sizeof(int *));
    for(int i=0;i<classes;i++)
    {
        cm[i] = (int *)malloc(classes * sizeof(int));
    }

	for(int i=0;i<classes;i++)
	{
			for(int j=0;j<classes;j++)
			{
					cm[i][j]=0;
			}
	}

	ifstream clfile;
    clfile.open(argv[1], ios::binary);
    vector<unsigned char> clbuffer;
    unsigned char tmp;
    int tmp2;

    clfile.seekg(0, ios::beg);

    for(int i=0;i<lines*columns; i++)
    {
		clfile.read(reinterpret_cast<char *>(&tmp2), sizeof(int));
	    tmp = (unsigned char)tmp2;
	    clbuffer.push_back(tmp);
    }
	clfile.close();

	ifstream trfile;
    trfile.open(argv[2], ios::binary);
    vector<unsigned char> trbuffer;

	trfile.seekg(0, ios::beg);

	for(int i=0;i<lines*columns; i++)
	{
		trfile.read(reinterpret_cast<char *>(&tmp), sizeof(unsigned char));
	    trbuffer.push_back(tmp);
    }
	trfile.close();

	for(int i=0;i<lines*columns; i++)
	{
	    if(trbuffer[i]!=255)
	    {
 		    count++;
 		    cm[trbuffer[i]-1][clbuffer[i]-1]++;
 		    if(trbuffer[i]==clbuffer[i])
 		    {
				 success++;
            }
	    }
	}

	ac = ((float) success) / ((float)count);

	for(int i=0;i<classes;i++)
	{
			for(int j=0;j<classes;j++)
			{
					cout << cm[i][j] << " ";
			}
			cout << endl;
	}

	cout << "Accuracy was computed to " << 100*ac << "%" << endl;
	return EXIT_SUCCESS;
}
