/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * mainsvmclass.cpp                                                       *
 * Version: 0.9.0 build 500                                               *
 * Last revised: 12/06/2006                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (a.k.a. Kalxas) 					  *
 * GCpp Lab.                                              March 2009      *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * */


//Usage :svm2class <level raster> <class attribute file> <svm training file> <svm prediction file> <classification.ers> <class raster>

#include "msegcore.h"

using namespace std;

int main (int argc, char *argv[]){

		if(argc != 7){//TODO: Add to Error Handling
            cout << "Error in arguments:" << endl << endl
            << "Usage :svm2class <level raster> <class attribute file> <svm training file> <svm prediction file> <classification.ers> <class raster>"
            << endl;
            system("Pause");
            exit(1);
        }

        //read arguments
        string arg1(argv[1]);
        string arg2(argv[2]);
        string arg3(argv[3]);
        string arg4(argv[4]);
        string arg5(argv[5]);
        string arg6(argv[6]);

		Level pass;
		ClassAttribute tmp;
		pass.LoadRaster(arg1);
		pass.LoadAttributes(arg2);
		pass.SVM2ClassificationMap(arg3, arg4);
		pass.SaveClassificationAsERS(arg5);
		pass.SaveClassificationMap(arg6);
		pass.clear();

}//main end

