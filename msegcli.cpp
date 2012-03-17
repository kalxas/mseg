/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegcli.cpp: Multiscale SEGmentation Command Line Interface            *
 * Version: 0.5.0 Final                                                   *
 * Last revised: 15/02/2005                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         February 2005   *
 *									  *
 *   Copyright (C) February 2005 Angelos Tzotsos <tzotsos@gmail.com>	  *
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

#include "libmseg.h"
#include "FreeImage.h"

using namespace std;

int main (int argc, char *argv[]){
 
        int merges_occured = 0;
        float average_size = 0.0;
		string bound = "boundary";
		string rst = "raster";
		
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
        
        //Run SPE module
        SPE(Img, mseg, 0);
        //Run First pass
		Level pass1;
        merges_occured = firstpass(pass1, Img, mseg, 1.0, mseg.LevelQueue[0].param);
		cout << "Merges occured at cycle " << 1 << " were " << merges_occured << endl;
  		
    	//Run Second pass
		Level pass2;
		merges_occured = secondpass(pass2, Img, mseg, pass1);
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
		    merges_occured = nthpass(passn1, Img, mseg, pass2);
		    cycle++;
		    cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
		    flag_npass = 1;
		}
		
		//Run nth pass
        while(merges_occured){
  		    if(pow(sw, cycle)<0.0){
  		        merges_occured = nthpass(passn2, Img, mseg, passn1);
  		        cycle++;
  		        cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
  		        flag_npass = 2;
  		        if(merges_occured) passn1.clear();
  		    }else{
   	    		merges_occured = nthpass(passn1, Img, mseg, passn2);
   	    		cycle++;
   	    		cout << "Merges occured at cycle " << cycle << " were " << merges_occured << endl;
   	    		flag_npass = 1;
   	    		if(merges_occured) passn2.clear();
  		    }
  		}    
  		
  		if(flag_npass == 1){
  		    passn2.SaveAsERS(arg2, Img);//the last merged level will be the previous one since no merges occured
  		    passn2.SaveRaster(rst);
        	average_size = (passn2.Lines*passn2.Columns)/(passn2.Objects.size());
            cout << "final level currently holds: " << passn2.Objects.size() << " objects" << endl;
  		    cout << "Average object size: " << average_size << endl;
        	cout << "Segmentation file saved..." << endl;
        	passn2.CreateBoundaryMap();
        	passn2.SaveBoundaryMapERS(bound);
  		}
    	
        if(flag_npass == 2){
        	passn1.SaveAsERS(arg2, Img);//the last merged level will be the previous one since no merges occured
  		    passn1.SaveRaster(rst);
       	    average_size = (passn1.Lines*passn1.Columns)/(passn1.Objects.size());
       	    cout << "final level currently holds: " << passn1.Objects.size() << " objects" << endl;
  		    cout << "Average object size: " << average_size << endl;
        	cout << "Segmentation file saved..." << endl;
        	passn1.CreateBoundaryMap();
        	passn1.SaveBoundaryMapERS(bound);
  		}   
      
		//close image file
        Img.file.close();
        //Debug mode:
        cout << "Image file closed..." << endl;
        //system("Pause");
        return 0;
   
}//main end



    	//Generic run...	         
        //For every level in LevelQueue, 
        	//if no sublevel run 1stpass
         	//if no sublevel and if something was merged, run 2ndpass
          	//if something was merged, run nthpass
           		//while new merges occur run nthpass
           //push level in LevelHierarchy
                     
        //Save levels to output file
        //write output ers header file


        
