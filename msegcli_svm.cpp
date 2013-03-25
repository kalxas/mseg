/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegcli_svm.cpp: Multiscale SEGmentation Command Line Interface (SVM)  *
 * Version: 0.9.x                                                         *
 * Last revised: 02/06/2006                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         June 2006       *
 *									  *
 *   Copyright (C) Angelos Tzotsos <tzotsos@gmail.com>	  		  *
 *									  *
 * This program is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU General Public License as published by   *
 * the Free Software Foundation, either version 2 of the License, or      *
 * (at your option) any later version.					  *
 *									  *
 * This program is distributed in the hope that it will be useful,	  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 	  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  *
 * GNU General Public License for more details.				  *
 *									  *
 * You should have received a copy of the GNU General Public License	  *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *									  *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * */


//Usage :mseg <input.ers> <output.ers> <image weights text file> <Level sequence and parameter text file>


//Future Usage : mseg <input.ers> <output.ers> <Level Name> <SuperLevel> <SubLevel>
//			   <Scale Parameter> <Compactness parameter> <Egde Compensation Parameter>
//			   <image weights text file> ..........
//The above is for single level creation
//To be implemented in version 1.1

#include "msegcore.h"
#include "msegstd.h"
#include "FreeImage.h"

using namespace std;

int main (int argc, char *argv[]){

        int merges_occured = 0;
        float average_size = 0.0;
		string bound = "boundary";
		string rst = "raster";
		string xml = "stats";
		string train = "training";
		string tta = "TTAMask";

		if(argc != 5){//TODO: Add to Error Handling
            cout << "Error in arguments:" << endl << endl
            << "Usage : mseg <input.ers> <output.ers> <image weights text file> <Level sequence and parameter text file>"
            << endl;
            system("Pause");
            exit(1);
        }

        //read arguments
        string arg1(argv[1]);
        string arg2(argv[2]);
        string arg3(argv[3]);
        string arg4(argv[4]);

        //Initialize one instance called mseg
        MSEG_Init mseg;

        //Fill LevelQueue
        ifstream LevelQueueList(arg4.c_str());
        LevelQueueObj tmpQueueObj;
        while(LevelQueueList){
            LevelQueueList >> tmpQueueObj.number;
            if(!LevelQueueList) break;
            LevelQueueList >> tmpQueueObj.sub_level_number;
            LevelQueueList >> tmpQueueObj.param.Scale;
            LevelQueueList >> tmpQueueObj.param.Color;
            LevelQueueList >> tmpQueueObj.param.Compact;
            LevelQueueList >> tmpQueueObj.param.Edge;
            LevelQueueList >> tmpQueueObj.param.EC;
            LevelQueueList >> tmpQueueObj.param.GHH;
            LevelQueueList >> tmpQueueObj.param.MA;
            mseg.LevelQueue.push_back(tmpQueueObj);
        }
        LevelQueueList.close();

        //Fill Queue_to_Hierarchy

        //Open image file for input
        ERS_Image Img(arg1, 0);


        //Store input band weights
        ifstream WeightsFile(arg3.c_str());
        if(!WeightsFile){
            cout << "Error reading file" << WeightsFile << endl;
            exit(1);
        }
        for(int i=1; i<=Img.Bands; i++){
            WeightsFile >> Img.BandWeight[i];
        }
        WeightsFile.close();

        Starting_Points_Estimation SPE(0);
        Standard_Mode speed_mode;

        //Run SPE module
        SPE.HSI(Img, mseg);
        //Run First pass
		Level pass1;
        merges_occured = speed_mode.firstpass(pass1, Img, mseg, 1.0, mseg.LevelQueue[0].param);
		cout << "Merges occured at cycle " << 1 << " were " << merges_occured << endl;

    	//Run Second pass
		Level pass2;
		merges_occured = speed_mode.secondpass(pass2, Img, mseg, pass1);
		cout << "Merges occured at cycle " << 2 << " were " << merges_occured << endl;

		//pass1.clear();
		//cout << "Level 1 deleted" << endl;


		//Declare 2 temp nth pass levels
		Level passn1;
		Level passn2;
		int cycle = 2;
		int sw = -1;//switch to use with power...
		int flag_npass;

  		//Run 3rd pass
  		if(merges_occured){//if second pass merged something, enter the nth pass loop...
		    merges_occured = speed_mode.nthpass(passn1, Img, mseg, pass2);
		    cycle++;
		    cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
		    flag_npass = 1;
		}

		//Run nth pass
        while(merges_occured){
  		    if(pow(sw, cycle)<0.0){
  		        merges_occured = speed_mode.nthpass(passn2, Img, mseg, passn1);
  		        cycle++;
  		        cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
  		        flag_npass = 2;
  		        if(merges_occured) passn1.clear();
  		    }else{
   	    		merges_occured = speed_mode.nthpass(passn1, Img, mseg, passn2);
   	    		cycle++;
   	    		cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
   	    		flag_npass = 1;
   	    		if(merges_occured) passn2.clear();
  		    }
  		}

  	if(flag_npass == 1){
  		    passn2.CalculateProperties(Img);
  		    cout << "Statistics calculated" << endl;
        	passn2.SaveMeanAsERS(arg2, Img);//the last merged level will be the previous one since no merges occured
  		    passn2.SaveRaster(rst, Img);
  		    //passn2.SaveMiniXML(xml);
  		    average_size = (passn2.Lines*passn2.Columns)/(passn2.Objects.size());
            cout << "final level currently holds: " << passn2.Objects.size() << " objects" << endl;
  		    cout << "Average object size: " << average_size << endl;
        	cout << "Segmentation file saved..." << endl;
        	passn2.CreateBoundaryMap();
        	passn2.SaveBoundaryMapERS(bound);
        	passn2.SaveXML(xml, Img);
        	passn2.SaveProperties(xml);
        	passn2.SaveSVMTesting(xml);
        	passn2.LoadTTA(tta);
        	passn2.SaveSVMTraining(train);
        	//passn2.SaveSVMTrainingVoting(train, 0.75);
  		}

        if(flag_npass == 2){
        	passn1.CalculateProperties(Img);
        	cout << "Statistics calculated" << endl;
        	passn1.SaveMeanAsERS(arg2, Img);//the last merged level will be the previous one since no merges occured
  		    passn1.SaveRaster(rst, Img);
  		    //passn1.SaveMiniXML(xml);
  		    average_size = (passn1.Lines*passn1.Columns)/(passn1.Objects.size());
       	    cout << "final level currently holds: " << passn1.Objects.size() << " objects" << endl;
  		    cout << "Average object size: " << average_size << endl;
        	cout << "Segmentation file saved..." << endl;
        	passn1.CreateBoundaryMap();
        	passn1.SaveBoundaryMapERS(bound);
        	passn1.SaveXML(xml, Img);
        	passn1.SaveProperties(xml);
        	passn1.SaveSVMTesting(xml);
        	passn1.LoadTTA(tta);
        	passn1.SaveSVMTraining(train);
        	//passn1.SaveSVMTrainingVoting(train, 0.75);
  		}

		//close image file
        Img.file.close();
        //Debug mode:
        cout << "Image file closed..." << endl;
        //system("Pause");
        return 0;

}//main end

