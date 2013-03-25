/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegstd.cpp: MSEG Fast mode class with edge compensation               *
 * Version: 0.9.x                                                         *
 * Last revised: 31/03/2010                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         March 2010      *
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


#include "msegfastedge.h"

int FastEdge_Mode::firstpass(Level& pass1, ERS_Image& Img, ERS_Image& ImgEdge, MSEG_Init& mseg, float level_ID, MSEG_Param& Param){

    int mergecounter = 0;
    int edgecounter = 0;
	pass1.raster.reserve(Img.Lines*Img.Columns);
    pass1.cycles = 1;
    pass1.Lines = Img.Lines;
    pass1.Columns = Img.Columns;

    mseg.last_id = (Img.Lines*Img.Columns); //Set the first object id

    //Let all the raster become -1
    for (int i=0; i<(Img.Lines*Img.Columns); i++){
        pass1.raster.push_back(-1);
        //pass1.raster[i] = -1;//Solution 2, without reallocating memory into the vector
    }

    //Copy the id number
    pass1.HierarchyID = level_ID;

    //Copy segmentation parameters into new level
    pass1.ScaleParameter = Param.Scale;
    pass1.Color = Param.Color;
    pass1.Compactness = Param.Compact;
    pass1.Edge = Param.Edge;

    //Declare a temporary object
    Object tmp;
    BoundaryLine tmpBoundaryLine;
	PriorityObject tmp_priority_object;
    pixel cur, up, down, right, left, best, best_up, best_down, best_right, best_left, best_best;
    int cur_id, up_id, down_id, right_id, left_id, best_id, best_best_id;
    int best_up_id, best_down_id, best_right_id, best_left_id;
    int i, c, l, pri, sec, matches, pri_edge;
    float best_h, up_h, down_h, left_h, right_h, best_up_h, best_down_h, best_right_h, best_left_h, best_best_h;
    float std;
    float cur_val, near_val;

    //Check the SPE result number
    if (mseg.StartingPoints.size() != Img.Macroblocks.size()){
        cout << "Error in SPE result" << endl << "Quiting at first pass..." << endl;
        exit(1);
    }
    int msize = Img.Macroblocks.size();

    //For every macroblock
    for(i=0;i<(Img.MacroLines * Img.MacroColumns);i++){
        //Assume each pixel is an object and assign them false ids in order to use the priority_queue
        //The pixel false id is the 1D index of the raster vector
        pri = 1;//During first pass only primary is used.(inside macroblock)
        sec = 1;
        pri_edge = 1000000;

        //Load the Starting Point
        cur.line = mseg.StartingPoints[i].line;
        cur.column = mseg.StartingPoints[i].column;

        cur_id = D1(Img.Lines, Img.Columns, cur.line, cur.column) + 1;

        //Check the boundary for the starting point
        if ((cur.line<Img.Macroblocks[i].LineStart || cur.line>Img.Macroblocks[i].LineEnd)
        ||(cur.column<Img.Macroblocks[i].ColumnStart || cur.column>Img.Macroblocks[i].ColumnEnd)){
            cout << "Starting point" << i << "out of bounds" << endl
            << "Quiting at first pass..." << endl;
            exit(1);
        }

        //Push the first pixel to the priority_queue
        tmp_priority_object.primary = pri;
        tmp_priority_object.secondary = sec;
        tmp_priority_object.id = cur_id;
        mseg.PriorityList.push(tmp_priority_object);
        pass1.raster[cur_id - 1] = 0; //Flag that it is loaded into priority list

        //Edge additions START
        for(l=Img.Macroblocks[i].LineStart;l<=Img.Macroblocks[i].LineEnd;l++){
            for(c=Img.Macroblocks[i].ColumnStart;c<=Img.Macroblocks[i].ColumnEnd;c++){
                if(ImgEdge.Buffer(l, c, 1) == 255.0){
                    cur.line = l;
                    cur.column = c;
                    cur_id = D1(Img.Lines, Img.Columns, cur.line, cur.column) + 1;

                    //Load neigbors
                    up.line = cur.line - 1;
                    up.column = cur.column;
                    up_id = (up.line >= 1) ? (D1(Img.Lines, Img.Columns, up.line, up.column) + 1) : 1;
                    down.line = cur.line + 1;
                    down.column = cur.column;
                    down_id = D1(Img.Lines, Img.Columns, down.line, down.column) + 1;
                    right.line = cur.line;
                    right.column = cur.column + 1;
                    right_id = D1(Img.Lines, Img.Columns, right.line, right.column) + 1;
                    left.line = cur.line;
                    left.column = cur.column - 1;
                    left_id = (left.column >= 1) ? (D1(Img.Lines, Img.Columns, left.line, left.column) + 1) : 1;


                    //Add its neighbors in the Priority list
                    if((up.line>=Img.Macroblocks[i].LineStart && up.line<=Img.Macroblocks[i].LineEnd)
                        && (up.column>=Img.Macroblocks[i].ColumnStart && up.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[up_id - 1] == -1 && ImgEdge.Buffer(up.line, up.column, 1) != 255.0){
                          tmp_priority_object.primary = ++pri_edge;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = up_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[up_id - 1] = 0;
                        }
                    }
                    if((down.line>=Img.Macroblocks[i].LineStart && down.line<=Img.Macroblocks[i].LineEnd)
                        && (down.column>=Img.Macroblocks[i].ColumnStart && down.column<=Img.Macroblocks[i].ColumnEnd)){

                        if (pass1.raster[down_id - 1] == -1 && ImgEdge.Buffer(down.line, down.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = down_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[down_id - 1] = 0;
                        }
                    }
                    if((right.line>=Img.Macroblocks[i].LineStart && right.line<=Img.Macroblocks[i].LineEnd)
                        && (right.column>=Img.Macroblocks[i].ColumnStart && right.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[right_id - 1] == -1 && ImgEdge.Buffer(right.line, right.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = right_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[right_id - 1] = 0;
                        }
                    }
                    if((left.line>=Img.Macroblocks[i].LineStart && left.line<=Img.Macroblocks[i].LineEnd)
                        && (left.column>=Img.Macroblocks[i].ColumnStart && left.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[left_id - 1] == -1 && ImgEdge.Buffer(left.line, left.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = left_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[left_id - 1] = 0;
                        }
                    }

                    //Create a single pixel object
                    tmp.id = cur_id;
                    tmp.area = 1;
                    tmp.edge = TRUE;

                    for (int b=1; b<=Img.Bands; b++){
                        cur_val = Img.Buffer(cur.line, cur.column, b);
                        tmp.Sum.push_back(cur_val);
                        tmp.SumSq.push_back((cur_val*cur_val));
                    }

                    tmpBoundaryLine.Line = cur.line;
                    tmpBoundaryLine.ColumnStart = cur.column;
                    tmpBoundaryLine.ColumnEnd = cur.column;
                    tmp.Boundary.push_back(tmpBoundaryLine);
                    pass1.Objects.insert(make_pair(tmp.id, tmp));
                    pass1.raster[cur_id - 1] = cur_id;
                    tmp.clear();
                    edgecounter++;
                }
            }
        }
        //Edge additions END

        //While the priority list is not empty
        while (!mseg.PriorityList.empty()){
            //Load the next priority object
            tmp_priority_object = mseg.PriorityList.top();
            pri = tmp_priority_object.primary;
            sec = tmp_priority_object.secondary;
            cur_id = tmp_priority_object.id;
            mseg.PriorityList.pop();
            //If pixel has been merged...goto end of while loop
            if(pass1.raster[cur_id - 1] > 0) continue;
            cur = D2(cur_id, Img.Columns);

            best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter
            matches = 0;

            //Load neigbors
            up.line = cur.line - 1;
            up.column = cur.column;
            up_id = (up.line >= 1) ? (D1(Img.Lines, Img.Columns, up.line, up.column) + 1) : 1;
            down.line = cur.line + 1;
            down.column = cur.column;
            down_id = D1(Img.Lines, Img.Columns, down.line, down.column) + 1;
            right.line = cur.line;
            right.column = cur.column + 1;
            right_id = D1(Img.Lines, Img.Columns, right.line, right.column) + 1;
            left.line = cur.line;
            left.column = cur.column - 1;
            left_id = (left.column >= 1) ? (D1(Img.Lines, Img.Columns, left.line, left.column) + 1) : 1;

            //Check neighbors
            //TODO: Put everything into arrays and loop 4 times...

            if((up.line>=Img.Macroblocks[i].LineStart && up.line<=Img.Macroblocks[i].LineEnd)
            && (up.column>=Img.Macroblocks[i].ColumnStart && up.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[up_id - 1] <= 0){
                      //Calculate heterogeneity
                      up_h = 0.0;
                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer(cur.line, cur.column, b);
                          near_val = Img.Buffer(up.line, up.column, b);
                          std = abs(cur_val - near_val)/2;
                          up_h += Img.BandWeight[b]* 2 * std;
                      }
                      //Is it a match? --> matches++ --> compare to best_h
                      if (up_h < Param.Scale){
                          matches++;
                          if (up_h < best_h){
                              best_h = up_h;
                              best_id = up_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[up_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//pri + 1
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = up_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[up_id - 1] = 0;
                      }
                }
            }

            if((down.line>=Img.Macroblocks[i].LineStart && down.line<=Img.Macroblocks[i].LineEnd)
            && (down.column>=Img.Macroblocks[i].ColumnStart && down.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[down_id - 1] <= 0){
                      //Calculate heterogeneity
                      down_h = 0.0;

                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer(cur.line, cur.column, b);
                          near_val = Img.Buffer(down.line, down.column, b);
                          std = abs(cur_val - near_val)/2;
                          down_h += Img.BandWeight[b]* 2 * std;
                      }
                      //Is it a match? --> matches++ --> compare to best_h
                      if (down_h < Param.Scale){
                          matches++;
                          if (down_h < best_h){
                              best_h = down_h;
                              best_id = down_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[down_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = down_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[down_id - 1] = 0;
                      }
                }
            }

            if((right.line>=Img.Macroblocks[i].LineStart && right.line<=Img.Macroblocks[i].LineEnd)
            && (right.column>=Img.Macroblocks[i].ColumnStart && right.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[right_id - 1] <= 0){
                      //Calculate heterogeneity
                      right_h = 0.0;

                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer(cur.line, cur.column, b);
                          near_val = Img.Buffer(right.line, right.column, b);
                          std = abs(cur_val - near_val)/2;
                          right_h += Img.BandWeight[b]* 2 * std;
                      }
                      //Is it a match? --> matches++ --> compare to best_h
                      if (right_h < Param.Scale){
                          matches++;
                          if (right_h < best_h){
                              best_h = right_h;
                              best_id = right_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[right_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = right_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[right_id - 1] = 0;
                      }
                }
            }

            if((left.line>=Img.Macroblocks[i].LineStart && left.line<=Img.Macroblocks[i].LineEnd)
            && (left.column>=Img.Macroblocks[i].ColumnStart && left.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[left_id - 1] <= 0){
                      //Calculate heterogeneity
                      left_h = 0.0;

                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer(cur.line, cur.column, b);
                          near_val = Img.Buffer(left.line, left.column, b);
                          std = abs(cur_val - near_val)/2;
                          left_h += Img.BandWeight[b]* 2 * std;
                      }
                      //Is it a match? --> matches++ --> compare to best_h
                      if (left_h < Param.Scale){
                          matches++;
                          if (left_h < best_h){
                              best_h = left_h;
                              best_id = left_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[left_id - 1] == -1){
                          tmp_priority_object.primary = pri+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = left_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[left_id - 1] = 0;
                      }
                }
            }

            //Push it's neigbors into the priority list if they are not merged (done above)
            //Check for match
            if(matches){
                //Check for the best match
                best = D2(best_id, Img.Columns);
                best_best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter

               	//Load neigbors
               	best_up.line = best.line-1;
               	best_up.column = best.column;
               	best_up_id = (best_up.line >= 1) ? (D1(Img.Lines, Img.Columns, best_up.line, best_up.column) + 1) : 1;
                best_down.line = best.line + 1;
               	best_down.column = best.column;
               	best_down_id = D1(Img.Lines, Img.Columns, best_down.line, best_down.column) + 1;
               	best_right.line = best.line;
               	best_right.column = best.column + 1;
               	best_right_id = D1(Img.Lines, Img.Columns, best_right.line, best_right.column) + 1;
               	best_left.line = best.line;
               	best_left.column = best.column -1;
               	best_left_id = (best_left.column >= 1) ? (D1(Img.Lines, Img.Columns, best_left.line, best_left.column) + 1) : 1;
                //Check neighbors
               	if((best_up.line>=Img.Macroblocks[i].LineStart && best_up.line<=Img.Macroblocks[i].LineEnd)
               	&& (best_up.column>=Img.Macroblocks[i].ColumnStart && best_up.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_up_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_up_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_up_id - 1] = 0;
                    }
                 }

                 if((best_down.line>=Img.Macroblocks[i].LineStart && best_down.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_down.column>=Img.Macroblocks[i].ColumnStart && best_down.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_down_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_down_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_down_id - 1] = 0;
                    }
                 }

                 if((best_right.line>=Img.Macroblocks[i].LineStart && best_right.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_right.column>=Img.Macroblocks[i].ColumnStart && best_right.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_right_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_right_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_right_id - 1] = 0;
                    }
                 }

                 if((best_left.line>=Img.Macroblocks[i].LineStart && best_left.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_left.column>=Img.Macroblocks[i].ColumnStart && best_left.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_left_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_left_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_left_id - 1] = 0;
                    }
                 }

                //If yes create object, return value --> Update raster, clear temp object
               	mseg.last_id++;//add 1 to last id to use to the current object
               	tmp.id = mseg.last_id;
               	tmp.area = 2;
               	tmp.edge = FALSE;
               	mergecounter++;
	   			//tmp.cycle = 1;

	   			for (int b=1; b<=Img.Bands; b++){
                    cur_val = Img.Buffer(cur.line, cur.column, b);
                    near_val = Img.Buffer(best.line, best.column, b);
                    tmp.Sum.push_back((cur_val+near_val));
                    tmp.SumSq.push_back((cur_val*cur_val)+(near_val*near_val));
                }

               	if(cur.line == best.line){//the object will have one BoundaryLine object
               	    tmpBoundaryLine.Line = cur.line;
                   	if(cur.column<best.column){
                   		tmpBoundaryLine.ColumnStart = cur.column;
                   		tmpBoundaryLine.ColumnEnd = best.column;
                	}
                   	else{
                   		tmpBoundaryLine.ColumnStart = best.column;
                   		tmpBoundaryLine.ColumnEnd = cur.column;
                   	}
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                else if(cur.line < best.line){//the object will have two BoundaryLine objects
                   	tmpBoundaryLine.Line = cur.line;
                   	tmpBoundaryLine.ColumnStart = cur.column;
                   	tmpBoundaryLine.ColumnEnd = cur.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);

                   	tmpBoundaryLine.Line = best.line;
                   	tmpBoundaryLine.ColumnStart = best.column;
                   	tmpBoundaryLine.ColumnEnd = best.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                else {
                    tmpBoundaryLine.Line = best.line;
                  	tmpBoundaryLine.ColumnStart = best.column;
                   	tmpBoundaryLine.ColumnEnd = best.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);

                   	tmpBoundaryLine.Line = cur.line;
                   	tmpBoundaryLine.ColumnStart = cur.column;
                   	tmpBoundaryLine.ColumnEnd = cur.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                pass1.Objects.insert(make_pair(tmp.id, tmp));//insert will work better here...
                pass1.raster[cur_id - 1] = tmp.id;
                pass1.raster[best_id - 1] = tmp.id;
                tmp.clear();
            }
            else{//if cur has no matches
                //create single pixel object --> Update raster, clear temp object
                tmp.id = cur_id;
                tmp.area = 1;
                tmp.edge = FALSE;
                //tmp.cycle = 1;

                for (int b=1; b<=Img.Bands; b++){
                    cur_val = Img.Buffer(cur.line, cur.column, b);
                    tmp.Sum.push_back(cur_val);
                    tmp.SumSq.push_back((cur_val*cur_val));
                }

                tmpBoundaryLine.Line = cur.line;
                tmpBoundaryLine.ColumnStart = cur.column;
                tmpBoundaryLine.ColumnEnd = cur.column;
                tmp.Boundary.push_back(tmpBoundaryLine);
                pass1.Objects.insert(make_pair(tmp.id, tmp));//insert here to
                pass1.raster[cur_id - 1] = cur_id;
                tmp.clear();
            }
        }
        //TODO:Check if some pixels have not been treated and create single pixel objects --> Update raster
    }
    //Create topology --> Load Topology
    cout << endl << "Calculating overall topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = pass1.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = pass1.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(pass1.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = pass1.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(pass1.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = pass1.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(pass1.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = pass1.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(pass1.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
	cout << "First pass completed succesfully!" << endl;
	cout << "Objects allocated in first pass Level: " << pass1.Objects.size() << endl;
	cout << "Edge Objects in this pass: " << edgecounter << endl;
	return mergecounter;
}

int FastEdge_Mode::firstpass(Level& pass1, ERS_Image& Img, ERS_Image& ImgEdge, MSEG_Init& mseg, float level_ID, MSEG_Param& Param, Level& super_level){

    int mergecounter = 0;
    int edgecounter = 0;
	pass1.raster.reserve(Img.Lines*Img.Columns);
    pass1.cycles = 1;
    pass1.Lines = Img.Lines;
    pass1.Columns = Img.Columns;

    mseg.last_id = (Img.Lines*Img.Columns); //Set the first object id

    //Let all the raster become -1
    for (int i=0; i<(Img.Lines*Img.Columns); i++){
        pass1.raster.push_back(-1);
        //pass1.raster[i] = -1;//Solution 2, without reallocating memory into the vector
    }

    //Copy the id number
    pass1.HierarchyID = level_ID;

    //Copy segmentation parameters into new level
    pass1.ScaleParameter = Param.Scale;
    pass1.Color = Param.Color;
    pass1.Compactness = Param.Compact;
    pass1.Edge = Param.Edge;

    //Declare a temporary object
    Object tmp;
    //Generic Declarations
    BoundaryLine tmpBoundaryLine;
	PriorityObject tmp_priority_object;
    pixel cur, up, down, right, left, best, best_up, best_down, best_right, best_left, best_best;
    int cur_id, up_id, down_id, right_id, left_id, best_id, best_best_id;
    int best_up_id, best_down_id, best_right_id, best_left_id;
    int i, c, l, pri, sec, matches, pri_edge;
    float best_h, up_h, down_h, left_h, right_h, best_up_h, best_down_h, best_right_h, best_left_h, best_best_h;
    float std;
    float cur_val, near_val;

    //Check the SPE result number
    if (mseg.StartingPoints.size() != Img.Macroblocks.size()){
        cout << "Error in SPE result" << endl << "Quiting at first pass..." << endl;
        exit(1);
    }
    int msize = Img.Macroblocks.size();

    //For every macroblock
    for(i=0;i<(Img.MacroLines * Img.MacroColumns);i++){
        //Assume each pixel is an object and assign them false ids in order to use the priority_queue
        //The pixel false id is the 1D [] representation inside the raster vector
        pri = 1;//During first pass only primary is usefull.(inside macroblock)
        sec = 1;
        pri_edge = 1000000;

        //Load the Starting Point
        cur.line = mseg.StartingPoints[i].line;
        cur.column = mseg.StartingPoints[i].column;

        cur_id = D1(Img.Lines, Img.Columns, cur.line, cur.column) + 1;

        //Check the boundary for the starting point
        if ((cur.line<Img.Macroblocks[i].LineStart || cur.line>Img.Macroblocks[i].LineEnd)
        ||(cur.column<Img.Macroblocks[i].ColumnStart || cur.column>Img.Macroblocks[i].ColumnEnd)){
            cout << "Starting point" << i << "out of bounds" << endl
            << "Quiting at first pass..." << endl;
            exit(1);
        }

        //Push the first pixel to the priority_queue
        tmp_priority_object.primary = pri;
        tmp_priority_object.secondary = sec;
        tmp_priority_object.id = cur_id;
        mseg.PriorityList.push(tmp_priority_object);
        pass1.raster[cur_id - 1] = 0; //Flag that it is loaded into priority list

        //Edge additions START
        for(l=Img.Macroblocks[i].LineStart;l<=Img.Macroblocks[i].LineEnd;l++){
            for(c=Img.Macroblocks[i].ColumnStart;c<=Img.Macroblocks[i].ColumnEnd;c++){
                if(ImgEdge.Buffer(l, c, 1) == 255.0){
                    cur.line = l;
                    cur.column = c;
                    cur_id = D1(Img.Lines, Img.Columns, cur.line, cur.column) + 1;

                    //Load neigbors
                    up.line = cur.line - 1;
                    up.column = cur.column;
                    up_id = (up.line >= 1) ? (D1(Img.Lines, Img.Columns, up.line, up.column) + 1) : 1;
                    down.line = cur.line + 1;
                    down.column = cur.column;
                    down_id = D1(Img.Lines, Img.Columns, down.line, down.column) + 1;
                    right.line = cur.line;
                    right.column = cur.column + 1;
                    right_id = D1(Img.Lines, Img.Columns, right.line, right.column) + 1;
                    left.line = cur.line;
                    left.column = cur.column - 1;
                    left_id = (left.column >= 1) ? (D1(Img.Lines, Img.Columns, left.line, left.column) + 1) : 1;


                    //Add its neighbors in the Priority list
                    if((up.line>=Img.Macroblocks[i].LineStart && up.line<=Img.Macroblocks[i].LineEnd)
                        && (up.column>=Img.Macroblocks[i].ColumnStart && up.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[up_id - 1] == -1 && ImgEdge.Buffer(up.line, up.column, 1) != 255.0){
                          tmp_priority_object.primary = ++pri_edge;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = up_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[up_id - 1] = 0;
                        }
                    }
                    if((down.line>=Img.Macroblocks[i].LineStart && down.line<=Img.Macroblocks[i].LineEnd)
                        && (down.column>=Img.Macroblocks[i].ColumnStart && down.column<=Img.Macroblocks[i].ColumnEnd)){

                        if (pass1.raster[down_id - 1] == -1 && ImgEdge.Buffer(down.line, down.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = down_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[down_id - 1] = 0;
                        }
                    }
                    if((right.line>=Img.Macroblocks[i].LineStart && right.line<=Img.Macroblocks[i].LineEnd)
                        && (right.column>=Img.Macroblocks[i].ColumnStart && right.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[right_id - 1] == -1 && ImgEdge.Buffer(right.line, right.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = right_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[right_id - 1] = 0;
                        }
                    }
                    if((left.line>=Img.Macroblocks[i].LineStart && left.line<=Img.Macroblocks[i].LineEnd)
                        && (left.column>=Img.Macroblocks[i].ColumnStart && left.column<=Img.Macroblocks[i].ColumnEnd)){
                        if (pass1.raster[left_id - 1] == -1 && ImgEdge.Buffer(left.line, left.column, 1) != 255.0){
                              tmp_priority_object.primary = ++pri_edge;
                              tmp_priority_object.secondary = sec;
                              tmp_priority_object.id = left_id;
                              mseg.PriorityList.push(tmp_priority_object);
                              pass1.raster[left_id - 1] = 0;
                        }
                    }

                    //Create a single pixel object
                    tmp.id = cur_id;
                    tmp.area = 1;
                    tmp.edge = TRUE;

                    for (int b=1; b<=Img.Bands; b++){
                        cur_val = Img.Buffer(cur.line, cur.column, b);
                        tmp.Sum.push_back(cur_val);
                        tmp.SumSq.push_back((cur_val*cur_val));
                    }

                    tmpBoundaryLine.Line = cur.line;
                    tmpBoundaryLine.ColumnStart = cur.column;
                    tmpBoundaryLine.ColumnEnd = cur.column;
                    tmp.Boundary.push_back(tmpBoundaryLine);
                    pass1.Objects.insert(make_pair(tmp.id, tmp));
                    pass1.raster[cur_id - 1] = cur_id;
                    tmp.clear();
                    edgecounter++;
                }
            }
        }
        //Edge additions END

        //While the priority list is not empty
        while (!mseg.PriorityList.empty()){
            //Load the next priority object
            tmp_priority_object = mseg.PriorityList.top();
            pri = tmp_priority_object.primary;
            sec = tmp_priority_object.secondary;
            cur_id = tmp_priority_object.id;
            mseg.PriorityList.pop();//Perhaps to the end????
            //If pixel has been merged...goto end of while loop
            if(pass1.raster[cur_id - 1] > 0) continue;
    		cur = D2(cur_id, Img.Columns);

            best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter
            matches = 0;

            //Load neigbors
            up.line = cur.line - 1;
            up.column = cur.column;
            up_id = (up.line >= 1) ? (D1(Img.Lines, Img.Columns, up.line, up.column) + 1) : 1;
            down.line = cur.line + 1;
            down.column = cur.column;
            down_id = D1(Img.Lines, Img.Columns, down.line, down.column) + 1;
            right.line = cur.line;
            right.column = cur.column + 1;
            right_id = D1(Img.Lines, Img.Columns, right.line, right.column) + 1;
            left.line = cur.line;
            left.column = cur.column - 1;
            left_id = (left.column >= 1) ? (D1(Img.Lines, Img.Columns, left.line, left.column) + 1) : 1;

            //Check neighbors
            //TODO: Put everything into arrays and loop 4 times...

            if((up.line>=Img.Macroblocks[i].LineStart && up.line<=Img.Macroblocks[i].LineEnd)
            && (up.column>=Img.Macroblocks[i].ColumnStart && up.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[up_id - 1] <= 0){
				    if(super_level.raster[up_id-1]==super_level.raster[cur_id-1]){
                        //Calculate heterogeneity
                        up_h = 0.0;
                        for (int b=1; b<=Img.Bands; b++){
                            cur_val = Img.Buffer(cur.line, cur.column, b);
                            near_val = Img.Buffer(up.line, up.column, b);
                            std = abs(cur_val - near_val)/2;
                            up_h += Img.BandWeight[b]* 2 * std;
                        }
                        //Is it a match? --> matches++ --> compare to best_h
                        if (up_h < Param.Scale){
                            matches++;
                            if (up_h < best_h){
                                best_h = up_h;
                                best_id = up_id;//We will copy the best outside this loop in order to spare copies
                            }
                        }
				    }
                    //if raster value <0 --> push pixel in queue
                    if (pass1.raster[up_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//pri + 1
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = up_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[up_id - 1] = 0;
                    }
                }
            }

            if((down.line>=Img.Macroblocks[i].LineStart && down.line<=Img.Macroblocks[i].LineEnd)
            && (down.column>=Img.Macroblocks[i].ColumnStart && down.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[down_id - 1] <= 0){
				    if(super_level.raster[down_id-1]==super_level.raster[cur_id-1]){
                        //Calculate heterogeneity
                        down_h = 0.0;

                        for (int b=1; b<=Img.Bands; b++){
                            cur_val = Img.Buffer(cur.line, cur.column, b);
                            near_val = Img.Buffer(down.line, down.column, b);
                            std = abs(cur_val - near_val)/2;
                            down_h += Img.BandWeight[b]* 2 * std;
                        }
                        //Is it a match? --> matches++ --> compare to best_h
                        if (down_h < Param.Scale){
                            matches++;
                            if (down_h < best_h){
                                best_h = down_h;
                                best_id = down_id;//We will copy the best outside this loop in order to spare copies
                            }
                        }
                     }
                     //if raster value <0 --> push pixel in queue
                     if (pass1.raster[down_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = down_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[down_id - 1] = 0;
                     }
                }
            }

            if((right.line>=Img.Macroblocks[i].LineStart && right.line<=Img.Macroblocks[i].LineEnd)
            && (right.column>=Img.Macroblocks[i].ColumnStart && right.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[right_id - 1] <= 0){
                    if(super_level.raster[right_id-1]==super_level.raster[cur_id-1]){
					    //Calculate heterogeneity
                        right_h = 0.0;

                        for (int b=1; b<=Img.Bands; b++){
                            cur_val = Img.Buffer(cur.line, cur.column, b);
                            near_val = Img.Buffer(right.line, right.column, b);
                            std = abs(cur_val - near_val)/2;
                            right_h += Img.BandWeight[b]* 2 * std;
                        }
                        //Is it a match? --> matches++ --> compare to best_h
                        if (right_h < Param.Scale){
                            matches++;
                            if (right_h < best_h){
                                best_h = right_h;
                                best_id = right_id;//We will copy the best outside this loop in order to spare copies
                            }
					    }
                    }
                    //if raster value <0 --> push pixel in queue
                    if (pass1.raster[right_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = right_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[right_id - 1] = 0;
                    }
                }
            }

            if((left.line>=Img.Macroblocks[i].LineStart && left.line<=Img.Macroblocks[i].LineEnd)
            && (left.column>=Img.Macroblocks[i].ColumnStart && left.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[left_id - 1] <= 0){
                    if(super_level.raster[left_id-1]==super_level.raster[cur_id-1]){
					    //Calculate heterogeneity
                        left_h = 0.0;

                        for (int b=1; b<=Img.Bands; b++){
                            cur_val = Img.Buffer(cur.line, cur.column, b);
                            near_val = Img.Buffer(left.line, left.column, b);
                            std = abs(cur_val - near_val)/2;
                            left_h += Img.BandWeight[b]* 2 * std;
                        }
                        //Is it a match? --> matches++ --> compare to best_h
                        if (left_h < Param.Scale){
                            matches++;
                            if (left_h < best_h){
                                best_h = left_h;
                                best_id = left_id;//We will copy the best outside this loop in order to spare copies
                            }
                        }
                    }
                    //if raster value <0 --> push pixel in queue
                    if (pass1.raster[left_id - 1] == -1){
                          tmp_priority_object.primary = pri+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = left_id;
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[left_id - 1] = 0;
                    }
                }
            }

            //Push it's neigbors into the priority list if they are not merged (done above)
            //Check for match
            if(matches){
                //Check for the best match
                best = D2(best_id, Img.Columns);
                best_best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter

               	//Load neigbors
               	best_up.line = best.line-1;
               	best_up.column = best.column;
               	best_up_id = (best_up.line >= 1) ? (D1(Img.Lines, Img.Columns, best_up.line, best_up.column) + 1) : 1;
                best_down.line = best.line + 1;
               	best_down.column = best.column;
               	best_down_id = D1(Img.Lines, Img.Columns, best_down.line, best_down.column) + 1;
               	best_right.line = best.line;
               	best_right.column = best.column + 1;
               	best_right_id = D1(Img.Lines, Img.Columns, best_right.line, best_right.column) + 1;
               	best_left.line = best.line;
               	best_left.column = best.column -1;
               	best_left_id = (best_left.column >= 1) ? (D1(Img.Lines, Img.Columns, best_left.line, best_left.column) + 1) : 1;
                //Check neighbors
               	if((best_up.line>=Img.Macroblocks[i].LineStart && best_up.line<=Img.Macroblocks[i].LineEnd)
               	&& (best_up.column>=Img.Macroblocks[i].ColumnStart && best_up.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_up_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_up_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_up_id - 1] = 0;
                    }
                 }

                 if((best_down.line>=Img.Macroblocks[i].LineStart && best_down.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_down.column>=Img.Macroblocks[i].ColumnStart && best_down.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_down_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_down_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_down_id - 1] = 0;
                    }
                 }

                 if((best_right.line>=Img.Macroblocks[i].LineStart && best_right.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_right.column>=Img.Macroblocks[i].ColumnStart && best_right.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_right_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_right_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_right_id - 1] = 0;
                    }
                 }

                 if((best_left.line>=Img.Macroblocks[i].LineStart && best_left.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_left.column>=Img.Macroblocks[i].ColumnStart && best_left.column<=Img.Macroblocks[i].ColumnEnd)){
                	if (pass1.raster[best_left_id - 1] == -1){
                        tmp_priority_object.primary = ++pri;//pri + 1
                        tmp_priority_object.secondary = sec;
                        tmp_priority_object.id = best_left_id;
                        mseg.PriorityList.push(tmp_priority_object);
                        pass1.raster[best_left_id - 1] = 0;
                    }
                 }

                //If yes create object, return value --> Update raster, clear temp object
                mseg.last_id++;//add 1 to last id to use to the current object
               	tmp.id = mseg.last_id;
               	tmp.area = 2;
               	tmp.edge = FALSE;
               	mergecounter++;
	   			//tmp.cycle = 1;

	   			for (int b=1; b<=Img.Bands; b++){
                    cur_val = Img.Buffer(cur.line, cur.column, b);
                    near_val = Img.Buffer(best.line, best.column, b);
                    tmp.Sum.push_back((cur_val+near_val));
                    tmp.SumSq.push_back((cur_val*cur_val)+(near_val*near_val));
                }

               	if(cur.line == best.line){//the object will have one BoundaryLine object
               	    tmpBoundaryLine.Line = cur.line;
                   	if(cur.column<best.column){
                   		tmpBoundaryLine.ColumnStart = cur.column;
                   		tmpBoundaryLine.ColumnEnd = best.column;
                	}
                   	else{
                   		tmpBoundaryLine.ColumnStart = best.column;
                   		tmpBoundaryLine.ColumnEnd = cur.column;
                   	}
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                else if(cur.line < best.line){//the object will have two BoundaryLine objects
                   	tmpBoundaryLine.Line = cur.line;
                   	tmpBoundaryLine.ColumnStart = cur.column;
                   	tmpBoundaryLine.ColumnEnd = cur.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);

                   	tmpBoundaryLine.Line = best.line;
                   	tmpBoundaryLine.ColumnStart = best.column;
                   	tmpBoundaryLine.ColumnEnd = best.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                else {
                    tmpBoundaryLine.Line = best.line;
                  	tmpBoundaryLine.ColumnStart = best.column;
                   	tmpBoundaryLine.ColumnEnd = best.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);

                   	tmpBoundaryLine.Line = cur.line;
                   	tmpBoundaryLine.ColumnStart = cur.column;
                   	tmpBoundaryLine.ColumnEnd = cur.column;
                   	tmp.Boundary.push_back(tmpBoundaryLine);
                }
                pass1.Objects.insert(make_pair(tmp.id, tmp));//insert will work better here...
                pass1.raster[cur_id - 1] = tmp.id;
                pass1.raster[best_id - 1] = tmp.id;
                tmp.clear();
        	}
            else{//if cur has no matches
                //create single pixel object --> Update raster, clear temp object
                tmp.id = cur_id;
                tmp.area = 1;
                tmp.edge = FALSE;
                //tmp.cycle = 1;

                for (int b=1; b<=Img.Bands; b++){
                    cur_val = Img.Buffer(cur.line, cur.column, b);
                    tmp.Sum.push_back(cur_val);
                    tmp.SumSq.push_back((cur_val*cur_val));
                }

                tmpBoundaryLine.Line = cur.line;
                tmpBoundaryLine.ColumnStart = cur.column;
                tmpBoundaryLine.ColumnEnd = cur.column;
                tmp.Boundary.push_back(tmpBoundaryLine);
                pass1.Objects.insert(make_pair(tmp.id, tmp));//insert here to
                pass1.raster[cur_id - 1] = cur_id;
                tmp.clear();
            }
        }
        //TODO:Check if some pixels have not been treated and create single pixel objects --> Update raster
    }
    //Create topology --> Load Topology
    cout << endl << "Calculating overall topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = pass1.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = pass1.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(pass1.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = pass1.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(pass1.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = pass1.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(pass1.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = pass1.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(pass1.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
	cout << "First pass completed succesfully!" << endl;
	cout << "Objects allocated in first pass Level: " << pass1.Objects.size() << endl;
	cout << "Edge Objects in this pass: " << edgecounter << endl;
	return mergecounter;
}


int FastEdge_Mode::secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass){

	int mergecounter = 0;
	int edgecounter = 0;

	pass2.raster.reserve(Img.Lines*Img.Columns);
	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    pass2.raster.push_back(-1);
	}
	pass2.cycles = 2;
    pass2.Lines = Img.Lines;
    pass2.Columns = Img.Columns;

    //Copy the id number
    pass2.HierarchyID = previous_pass.HierarchyID;

    //Copy segmentation parameters into new level
    pass2.ScaleParameter = previous_pass.ScaleParameter;
    pass2.Color = previous_pass.Color;
    pass2.Compactness = previous_pass.Compactness;
    pass2.Edge = previous_pass.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id, near_near_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;
    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();
       	//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;

        //If not calculated, calculate object's perimeter

		best_f = pass2.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

        //Edge additions START
        if(previous_pass.Objects[cur_id].edge == TRUE){
            //Load neighbors in Priority list
            for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
                near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
                if(!previous_pass.Objects[near_id].priority_loaded){
                    tmp_priority_object.id = near_id;
                    tmp_priority_object.primary = pri + 1;
                    tmp_priority_object.secondary = sec;
                    mseg.PriorityList.push(tmp_priority_object);
                    previous_pass.Objects[near_id].priority_loaded = TRUE;
                }
            }
            //Copy the object to the next level
            pass2.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            pass2.Objects[cur_id].priority_loaded = FALSE;
            pass2.Objects[cur_id].edge = TRUE;
            pass2.Objects[cur_id].Neighbors.clear();
            //Update raster
            pass2.UpdateRaster(cur_id);
            edgecounter++;
            //Move to next object
            continue;
        }
        //Edge additions END

        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;

            //Edge additions START
            if(previous_pass.Objects[near_id].edge == TRUE && previous_pass.Objects[near_id].merged == FALSE){
                //Load neighbors in Priority list
                for(int i=0; i<previous_pass.Objects[near_id].Neighbors.size(); i++){
                    near_near_id = previous_pass.Objects[near_id].Neighbors[i].id;
                    if(!previous_pass.Objects[near_near_id].priority_loaded){
                        tmp_priority_object.id = near_near_id;
                        tmp_priority_object.primary = pri + 1;
                        tmp_priority_object.secondary = sec;
                        mseg.PriorityList.push(tmp_priority_object);
                        previous_pass.Objects[near_near_id].priority_loaded = TRUE;
                    }
                }
                //Copy the object to the next level
                pass2.Objects.insert(make_pair(near_id, previous_pass.Objects[near_id]));
                previous_pass.Objects[near_id].merged = TRUE;
                pass2.Objects[near_id].priority_loaded = FALSE;
                pass2.Objects[near_id].edge = TRUE;
                pass2.Objects[near_id].Neighbors.clear();
                //Update raster
                pass2.UpdateRaster(near_id);
                edgecounter++;
                //Move to next neighbor object
                continue;
            }
            //Edge additions END

            //if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			}
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
                 if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
                     //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     f = h_color;
                 }
                 else{
                     //get the heterogeneity already calculated
                     //This is not working properly yet, so we don't realy use it here...
                     f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }
                 if(f<best_f){
                    matches++;
                    best_f = f;
                    best_id = near_id;
                 }
   		    }
	 	}
	    //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = pass2.ScaleParameter;
	 	    //For every neighbor //Optimize: This loop can be skiped in fast mode only!!!
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
	 	        if(!previous_pass.Objects[near_id].priority_loaded){
	 	        	tmp_priority_object.id = near_id;
	 	        	tmp_priority_object.primary = pri + 1;
	 	        	tmp_priority_object.secondary = sec;
	 	        	mseg.PriorityList.push(tmp_priority_object);
	 	        	previous_pass.Objects[near_id].priority_loaded = TRUE;
	 	        }
	 	    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            //calculate perimeter
            tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            //tmp.merge_neighbors(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            tmp.priority_loaded = FALSE;
            pass2.Objects.insert(make_pair(tmp.id, tmp));

            previous_pass.Objects[cur_id].merged = TRUE;
            previous_pass.Objects[best_id].merged = TRUE;

            //Update raster
            pass2.UpdateRaster(tmp.id);

            tmp.clear();
      	}
        else{//no matches
            //push the same object to the next level
            pass2.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            pass2.Objects[cur_id].priority_loaded = FALSE;

            //--------> SOS <------------ Copy that everywere!!!!!
            pass2.Objects[cur_id].Neighbors.clear();

            //cout << "No matches, object copied to the next level" << endl;//Debug mode
            //Update raster
            pass2.UpdateRaster(cur_id);
        }

    }

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = pass2.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = pass2.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(pass2.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = pass2.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(pass2.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = pass2.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(pass2.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = pass2.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(pass2.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
    cout << "Edge Objects in this pass: " << edgecounter << endl;
    return mergecounter;
}

int FastEdge_Mode::secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level){

	int mergecounter = 0;
	int edgecounter = 0;
	pass2.raster.reserve(Img.Lines*Img.Columns);
	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    pass2.raster.push_back(-1);
	}
	pass2.cycles = 2;
    pass2.Lines = Img.Lines;
    pass2.Columns = Img.Columns;

    //Copy the id number
    pass2.HierarchyID = previous_pass.HierarchyID;

    //Copy segmentation parameters into new level
    pass2.ScaleParameter = previous_pass.ScaleParameter;
    pass2.Color = previous_pass.Color;
    pass2.Compactness = previous_pass.Compactness;
    pass2.Edge = previous_pass.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id, near_near_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;
    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();
       	//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;

        //If not calculated, calculate object's perimeter

		best_f = pass2.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

		//Edge additions START
        if(previous_pass.Objects[cur_id].edge == TRUE){
            //Load neighbors in Priority list
            for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
                near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
                if(!previous_pass.Objects[near_id].priority_loaded){
                    tmp_priority_object.id = near_id;
                    tmp_priority_object.primary = pri + 1;
                    tmp_priority_object.secondary = sec;
                    mseg.PriorityList.push(tmp_priority_object);
                    previous_pass.Objects[near_id].priority_loaded = TRUE;
                }
            }
            //Copy the object to the next level
            pass2.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            pass2.Objects[cur_id].priority_loaded = FALSE;
            pass2.Objects[cur_id].edge = TRUE;
            pass2.Objects[cur_id].Neighbors.clear();
            //Update raster
            pass2.UpdateRaster(cur_id);
            edgecounter++;
            //Move to next object
            continue;
        }
        //Edge additions END

        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;

            //Edge additions START
            if(previous_pass.Objects[near_id].edge == TRUE && previous_pass.Objects[near_id].merged == FALSE){
                //Load neighbors in Priority list
                for(int i=0; i<previous_pass.Objects[near_id].Neighbors.size(); i++){
                    near_near_id = previous_pass.Objects[near_id].Neighbors[i].id;
                    if(!previous_pass.Objects[near_near_id].priority_loaded){
                        tmp_priority_object.id = near_near_id;
                        tmp_priority_object.primary = pri + 1;
                        tmp_priority_object.secondary = sec;
                        mseg.PriorityList.push(tmp_priority_object);
                        previous_pass.Objects[near_near_id].priority_loaded = TRUE;
                    }
                }
                //Copy the object to the next level
                pass2.Objects.insert(make_pair(near_id, previous_pass.Objects[near_id]));
                previous_pass.Objects[near_id].merged = TRUE;
                pass2.Objects[near_id].priority_loaded = FALSE;
                pass2.Objects[near_id].edge = TRUE;
                pass2.Objects[near_id].Neighbors.clear();
                //Update raster
                pass2.UpdateRaster(near_id);
                edgecounter++;
                //Move to next neighbor object
                continue;
            }
            //Edge additions END

            //if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			}
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
				 if((super_level.raster[D1(Img.Lines, Img.Columns, previous_pass.Objects[cur_id].Boundary[0].Line, previous_pass.Objects[cur_id].Boundary[0].ColumnStart)])==(super_level.raster[D1(Img.Lines, Img.Columns, previous_pass.Objects[near_id].Boundary[0].Line, previous_pass.Objects[near_id].Boundary[0].ColumnStart)])){
	 			     if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             	//Virtual merge into temp object
                     	tmp.clear();
                     	tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     	h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     	f = h_color;
                     }
                 	 else{
                     	  //get the heterogeneity already calculated
                     	  //This is not working properly yet, so we don't realy use it here...
                     	  f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
  			  	     }
   			     	 if(f<best_f){
      		     	 	matches++;
      		     		best_f = f;
  			   			best_id = near_id;
       				 }
				 }
   		     }
	 	}
	 	//Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = pass2.ScaleParameter;
	 	    //For every neighbor //Optimize: This loop can be removed in fast mode
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
	 	        if(!previous_pass.Objects[near_id].priority_loaded){
	 	        	tmp_priority_object.id = near_id;
	 	        	tmp_priority_object.primary = pri + 1;
	 	        	tmp_priority_object.secondary = sec;
	 	        	mseg.PriorityList.push(tmp_priority_object);
	 	        	previous_pass.Objects[near_id].priority_loaded = TRUE;
	 	        }

	 	    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            //calculate perimeter
            tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            //tmp.merge_neighbors(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            tmp.priority_loaded = FALSE;
            pass2.Objects.insert(make_pair(tmp.id, tmp));

            previous_pass.Objects[cur_id].merged = TRUE;
            previous_pass.Objects[best_id].merged = TRUE;

            //Update raster
            pass2.UpdateRaster(tmp.id);

            tmp.clear();
      	}
        else{//no matches
            //push the same object to the next level
            pass2.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            pass2.Objects[cur_id].priority_loaded = FALSE;
            pass2.Objects[cur_id].Neighbors.clear();
            //Update raster
            pass2.UpdateRaster(cur_id);
        }

    }

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = pass2.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = pass2.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(pass2.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = pass2.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(pass2.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = pass2.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(pass2.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = pass2.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(pass2.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        pass2.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
    cout << "Edge Objects in this pass: " << edgecounter << endl;
	return mergecounter;

}

//TODO
int FastEdge_Mode::nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level){

	int mergecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = 1;

    //Copy the id number
    passn.HierarchyID = level_ID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;

    //Copy segmentation parameters into new level
    passn.ScaleParameter = Param.Scale;
    passn.Color = Param.Color;
    passn.Compactness = Param.Compact;
    passn.Edge = Param.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object(0,0,0);
    //TopologyObject tmpTopologyObject;
    //TopologyPair tmpTopologyPair;

    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = sub_level.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		sub_level.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();

        //Reload Priority object if cur already merged
        if(sub_level.Objects[cur_id].merged == TRUE) continue;

        //If not calculated, calculate object's perimeter
        if(sub_level.Objects[cur_id].perimeter == 0){
  		    sub_level.CalculatePerimeter(cur_id);
		}

        //Update and Fix the neighbors for the previous level in order to call once bellow
		//update_neighbors(sub_level.Objects[cur_id], mseg);
		//sub_level.Objects[cur_id].fix_neighbors();

		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

        //Load neigbors --> DONE: dahh, we have objects now...

        //For every neighbor
		for(int i=0; i<sub_level.Objects[cur_id].Neighbors.size(); i++){
            near_id = sub_level.Objects[cur_id].Neighbors[i].id;

            //If not calculated, calculate neighbor object's perimeter
            if(sub_level.Objects[near_id].perimeter == 0){
         	    sub_level.CalculatePerimeter(near_id);
		    }
		    //update_neighbors(sub_level.Objects[near_id], mseg);
		    //sub_level.Objects[near_id].fix_neighbors();

			//if neighbor not pushed to priority --> Push
			if(!sub_level.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    sub_level.Objects[near_id].priority_loaded = TRUE;
			}
			if(!sub_level.Objects[near_id].merged){//if Neighbor not merged
				 if (!sub_level.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(sub_level.Objects[cur_id], sub_level.Objects[near_id]);
                     //Calculate merge object's perimeter
                     tmp.perimeter = sub_level.CalculatePerimeter(cur_id, near_id);

                     h_color = tmp.Color_Heterogeneity(Img) - sub_level.Objects[cur_id].Color_Heterogeneity(Img) - sub_level.Objects[near_id].Color_Heterogeneity(Img);
                     h_smooth = tmp.Smoothness() - sub_level.Objects[cur_id].Smoothness() - sub_level.Objects[near_id].Smoothness();
                     h_cmpct = tmp.Compactness() - sub_level.Objects[cur_id].Compactness() - sub_level.Objects[near_id].Compactness();
                     h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                 }
                 else{
                     //get the heterogeneity already calculated
                     f = sub_level.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }
   			     if(f<best_f){
      		     	matches++;
  			   		best_f = f;
  			   		best_id = near_id;
   		         }
   		     }
	 	}

        //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor
	 	    for(int i=0; i<sub_level.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = sub_level.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
		        if(!sub_level.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	sub_level.Objects[near_id].priority_loaded = TRUE;
	        	}
	 	    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(sub_level.Objects[best_id], sub_level.Objects[cur_id]);
            //Calculate perimeter
            tmp.perimeter = sub_level.CalculatePerimeter(best_id, cur_id);
            //Merge neighbors
            //tmp.merge_neighbors(sub_level.Objects[best_id], sub_level.Objects[cur_id]);

			mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            passn.Objects.insert(make_pair(tmp.id, tmp));
            //tmpTopologyObject.id = tmp.id;
            //tmpTopologyObject.merged = FALSE;
            //tmpTopologyObject.new_id = 0;
            //mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));

            //tmpTopologyPair.from = best_id;
            //tmpTopologyPair.to = tmp.id;
            //mseg.Temp_Topology.push_back(tmpTopologyPair);
            //tmpTopologyPair.from = cur_id;
            //tmpTopologyPair.to = tmp.id;
            //mseg.Temp_Topology.push_back(tmpTopologyPair);

            sub_level.Objects[cur_id].merged = TRUE;
            sub_level.Objects[best_id].merged = TRUE;

            //Update raster
            passn.UpdateRaster(tmp.id);

            tmp.clear();
   		}
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, sub_level.Objects[cur_id]));
            sub_level.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].Neighbors.clear();
            //Update raster
            passn.UpdateRaster(cur_id);
        }

    }
    //Update the Topology map

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = passn.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = passn.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(passn.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = passn.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(passn.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = passn.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(passn.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = passn.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(passn.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }

    /*OLD stuff
    int t;
    for(int i=0; i<mseg.Temp_Topology.size(); i++){
        t = mseg.Temp_Topology[i].from;
        mseg.Topology[t].merged = TRUE;
        mseg.Topology[t].new_id = mseg.Temp_Topology[i].to;
    }
    mseg.Temp_Topology.clear();
    */
    return mergecounter;
}
//TODO
int FastEdge_Mode::nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level, Level& super_level){

	int mergecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = 1;

    //Copy the id number
    passn.HierarchyID = level_ID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;

    //Copy segmentation parameters into new level
    passn.ScaleParameter = Param.Scale;
    passn.Color = Param.Color;
    passn.Compactness = Param.Compact;
    passn.Edge = Param.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object(0,0,0);
    //TopologyObject tmpTopologyObject;
    //TopologyPair tmpTopologyPair;

    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = sub_level.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		sub_level.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();

        //Reload Priority object if cur already merged
        if(sub_level.Objects[cur_id].merged == TRUE) continue;

        //If not calculated, calculate object's perimeter
        if(sub_level.Objects[cur_id].perimeter == 0){
  		    sub_level.CalculatePerimeter(cur_id);
		}

        //Update and Fix the neighbors for the previous level in order to call once bellow
		//update_neighbors(sub_level.Objects[cur_id], mseg);
		//sub_level.Objects[cur_id].fix_neighbors();

		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

        //Load neigbors --> DONE: dahh, we have objects now...

        //For every neighbor
		for(int i=0; i<sub_level.Objects[cur_id].Neighbors.size(); i++){
            near_id = sub_level.Objects[cur_id].Neighbors[i].id;

            //If not calculated, calculate neighbor object's perimeter
            if(sub_level.Objects[near_id].perimeter == 0){
         	    sub_level.CalculatePerimeter(near_id);
		    }
		    //update_neighbors(sub_level.Objects[near_id], mseg);
		    //sub_level.Objects[near_id].fix_neighbors();

			//if neighbor not pushed to priority --> Push
			if(!sub_level.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    sub_level.Objects[near_id].priority_loaded = TRUE;
			}
			if(!sub_level.Objects[near_id].merged){//if Neighbor not merged
				 if((super_level.raster[D1(Img.Lines, Img.Columns, sub_level.Objects[cur_id].Boundary[0].Line, sub_level.Objects[cur_id].Boundary[0].ColumnStart)])==(super_level.raster[D1(Img.Lines, Img.Columns, sub_level.Objects[near_id].Boundary[0].Line, sub_level.Objects[near_id].Boundary[0].ColumnStart)])){
				     if (!sub_level.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             	 //Virtual merge into temp object
                     	 tmp.clear();
                     	 tmp.merge(sub_level.Objects[cur_id], sub_level.Objects[near_id]);
                     	 //Calculate merge object's perimeter
                     	 tmp.perimeter = sub_level.CalculatePerimeter(cur_id, near_id);

                     	 h_color = tmp.Color_Heterogeneity(Img) - sub_level.Objects[cur_id].Color_Heterogeneity(Img) - sub_level.Objects[near_id].Color_Heterogeneity(Img);
                     	 h_smooth = tmp.Smoothness() - sub_level.Objects[cur_id].Smoothness() - sub_level.Objects[near_id].Smoothness();
                     	 h_cmpct = tmp.Compactness() - sub_level.Objects[cur_id].Compactness() - sub_level.Objects[near_id].Compactness();
                     	 h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     	 f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
               		 }
                 	 else{
                     	  //get the heterogeneity already calculated
                     	  f = sub_level.Objects[cur_id].Neighbors[i].merging_heterogeneity;
		   		     }
   			     	 if(f<best_f){
	     	 			  matches++;
  			   			  best_f = f;
  			   			  best_id = near_id;
         		     }
				 }
   		     }
	 	}

        //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor
	 	    for(int i=0; i<sub_level.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = sub_level.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
		        if(!sub_level.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	sub_level.Objects[near_id].priority_loaded = TRUE;
	        	}
	 	    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(sub_level.Objects[best_id], sub_level.Objects[cur_id]);
            //Calculate perimeter
            tmp.perimeter = sub_level.CalculatePerimeter(best_id, cur_id);
            //Merge neighbors
            //tmp.merge_neighbors(sub_level.Objects[best_id], sub_level.Objects[cur_id]);

			mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            passn.Objects.insert(make_pair(tmp.id, tmp));
            //tmpTopologyObject.id = tmp.id;
            //tmpTopologyObject.merged = FALSE;
            //tmpTopologyObject.new_id = 0;
            //mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));

            //tmpTopologyPair.from = best_id;
            //tmpTopologyPair.to = tmp.id;
            //mseg.Temp_Topology.push_back(tmpTopologyPair);
            //tmpTopologyPair.from = cur_id;
            //tmpTopologyPair.to = tmp.id;
            //mseg.Temp_Topology.push_back(tmpTopologyPair);

            sub_level.Objects[cur_id].merged = TRUE;
            sub_level.Objects[best_id].merged = TRUE;

            //Update raster
            passn.UpdateRaster(tmp.id);

            tmp.clear();
   		}
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, sub_level.Objects[cur_id]));
            sub_level.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].Neighbors.clear();

            //Update raster
            passn.UpdateRaster(cur_id);
        }

    }
    //Update the Topology map

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = passn.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = passn.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(passn.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = passn.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(passn.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = passn.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(passn.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = passn.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(passn.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }

    /*OLD stuff
    int t;
    for(int i=0; i<mseg.Temp_Topology.size(); i++){
        t = mseg.Temp_Topology[i].from;
        mseg.Topology[t].merged = TRUE;
        mseg.Topology[t].new_id = mseg.Temp_Topology[i].to;
    }
    mseg.Temp_Topology.clear();
    */
    return mergecounter;

}

int FastEdge_Mode::nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass){

	int mergecounter = 0;
	int edgecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = previous_pass.cycles++;
   	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    passn.raster.push_back(-1);
	}
    //Copy the id number
    passn.HierarchyID = previous_pass.HierarchyID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;

    //Copy segmentation parameters into new level
    passn.ScaleParameter = previous_pass.ScaleParameter;
    passn.Color = previous_pass.Color;
    passn.Compactness = previous_pass.Compactness;
    passn.Edge = previous_pass.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id, near_near_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;

    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();

		//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;

		//If not calculated, calculate object's perimeter
        if(previous_pass.Objects[cur_id].perimeter == 0){
  		    previous_pass.CalculatePerimeter(cur_id);
		}

		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

        //Edge additions START
        if(previous_pass.Objects[cur_id].edge == TRUE){
            //Load neighbors in Priority list
            for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
                near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
                if(!previous_pass.Objects[near_id].priority_loaded){
                    tmp_priority_object.id = near_id;
                    tmp_priority_object.primary = pri + 1;
                    tmp_priority_object.secondary = sec;
                    mseg.PriorityList.push(tmp_priority_object);
                    previous_pass.Objects[near_id].priority_loaded = TRUE;
                }
            }
            //Copy the object to the next level
            passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].edge = TRUE;
            passn.Objects[cur_id].Neighbors.clear();
            //Update raster
            passn.UpdateRaster(cur_id);
            edgecounter++;
            //Move to next object
            continue;
        }
        //Edge additions END

        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;

            //Edge additions START
            if(previous_pass.Objects[near_id].edge == TRUE && previous_pass.Objects[near_id].merged == FALSE){
                //Load neighbors in Priority list
                for(int i=0; i<previous_pass.Objects[near_id].Neighbors.size(); i++){
                    near_near_id = previous_pass.Objects[near_id].Neighbors[i].id;
                    if(!previous_pass.Objects[near_near_id].priority_loaded){
                        tmp_priority_object.id = near_near_id;
                        tmp_priority_object.primary = pri + 1;
                        tmp_priority_object.secondary = sec;
                        mseg.PriorityList.push(tmp_priority_object);
                        previous_pass.Objects[near_near_id].priority_loaded = TRUE;
                    }
                }
                //Copy the object to the next level
                passn.Objects.insert(make_pair(near_id, previous_pass.Objects[near_id]));
                previous_pass.Objects[near_id].merged = TRUE;
                passn.Objects[near_id].priority_loaded = FALSE;
                passn.Objects[near_id].edge = TRUE;
                passn.Objects[near_id].Neighbors.clear();
                //Update raster
                passn.UpdateRaster(near_id);
                edgecounter++;
                //Move to next neighbor object
                continue;
            }
            //Edge additions END

            //If not calculated, calculate neighbor object's perimeter
            if(previous_pass.Objects[near_id].perimeter == 0){
         	    previous_pass.CalculatePerimeter(near_id);
		    }
		    //if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			}
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
				if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     //Calculate merge object's perimeter
                     tmp.perimeter = previous_pass.CalculatePerimeter(cur_id, near_id);
                     //Calculate heterogeneity
                     h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     h_smooth = tmp.Smoothness() - previous_pass.Objects[cur_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                     h_cmpct = tmp.Compactness() - previous_pass.Objects[cur_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                     h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                 }
                 else{
                     //get the heterogeneity already calculated
                     f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }
   			     if(f<best_f){
      		     	matches++;
      		     	best_f = f;
  			   		best_id = near_id;
   		         }
   		     }
	 	}
	 	//Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor //Optimize: This loop can be removed in fast mode
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
		        if(!previous_pass.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	previous_pass.Objects[near_id].priority_loaded = TRUE;
		        }
		    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            //Calculate perimeter
            tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            passn.Objects.insert(make_pair(tmp.id, tmp));//insert here

            previous_pass.Objects[cur_id].merged = TRUE;
            previous_pass.Objects[best_id].merged = TRUE;

            //Update raster
            passn.UpdateRaster(tmp.id);

            tmp.clear();
        }
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));//insert here to
            previous_pass.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].Neighbors.clear();

            //Update raster
            passn.UpdateRaster(cur_id);
        }
    }

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            //cout << "Inside topology loop: " << row << " and " << col << endl;//debug mode
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;
            //cout << "1" << endl;//Debug mode
            cur_id = passn.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = passn.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(passn.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            //cout << "2" << endl;//Debug mode
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = passn.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(passn.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            //cout << "3" << endl;//Debug mode
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = passn.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(passn.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            //cout << "4" << endl;//Debug mode
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = passn.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(passn.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
    cout << "Edge Objects in this pass: " << edgecounter << endl;
    return mergecounter;
}

int FastEdge_Mode::nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level){

	int mergecounter = 0;
	int edgecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = previous_pass.cycles++;
   	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    passn.raster.push_back(-1);
	}
    //Copy the id number
    passn.HierarchyID = previous_pass.HierarchyID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;

    //Copy segmentation parameters into new level
    passn.ScaleParameter = previous_pass.ScaleParameter;
    passn.Color = previous_pass.Color;
    passn.Compactness = previous_pass.Compactness;
    passn.Edge = previous_pass.Edge;

    //Declare a temporary object
    Object tmp;
    pixel cur, up, down, right, left;
    int up_id, down_id, right_id, left_id;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id, near_near_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;

    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }

    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){

		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();

		//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;

		//If not calculated, calculate object's perimeter
        if(previous_pass.Objects[cur_id].perimeter == 0){
  		    previous_pass.CalculatePerimeter(cur_id);
		}

		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;

        //Edge additions START
        if(previous_pass.Objects[cur_id].edge == TRUE){
            //Load neighbors in Priority list
            for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
                near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
                if(!previous_pass.Objects[near_id].priority_loaded){
                    tmp_priority_object.id = near_id;
                    tmp_priority_object.primary = pri + 1;
                    tmp_priority_object.secondary = sec;
                    mseg.PriorityList.push(tmp_priority_object);
                    previous_pass.Objects[near_id].priority_loaded = TRUE;
                }
            }
            //Copy the object to the next level
            passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].edge = TRUE;
            passn.Objects[cur_id].Neighbors.clear();
            //Update raster
            passn.UpdateRaster(cur_id);
            edgecounter++;
            //Move to next object
            continue;
        }
        //Edge additions END

        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;

            //Edge additions START
            if(previous_pass.Objects[near_id].edge == TRUE && previous_pass.Objects[near_id].merged == FALSE){
                //Load neighbors in Priority list
                for(int i=0; i<previous_pass.Objects[near_id].Neighbors.size(); i++){
                    near_near_id = previous_pass.Objects[near_id].Neighbors[i].id;
                    if(!previous_pass.Objects[near_near_id].priority_loaded){
                        tmp_priority_object.id = near_near_id;
                        tmp_priority_object.primary = pri + 1;
                        tmp_priority_object.secondary = sec;
                        mseg.PriorityList.push(tmp_priority_object);
                        previous_pass.Objects[near_near_id].priority_loaded = TRUE;
                    }
                }
                //Copy the object to the next level
                passn.Objects.insert(make_pair(near_id, previous_pass.Objects[near_id]));
                previous_pass.Objects[near_id].merged = TRUE;
                passn.Objects[near_id].priority_loaded = FALSE;
                passn.Objects[near_id].edge = TRUE;
                passn.Objects[near_id].Neighbors.clear();
                //Update raster
                passn.UpdateRaster(near_id);
                edgecounter++;
                //Move to next neighbor object
                continue;
            }
            //Edge additions END

            //If not calculated, calculate neighbor object's perimeter
            if(previous_pass.Objects[near_id].perimeter == 0){
         	    previous_pass.CalculatePerimeter(near_id);
		    }

			//if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			}
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
				if((super_level.raster[D1(Img.Lines, Img.Columns, previous_pass.Objects[cur_id].Boundary[0].Line, previous_pass.Objects[cur_id].Boundary[0].ColumnStart)])==(super_level.raster[D1(Img.Lines, Img.Columns, previous_pass.Objects[near_id].Boundary[0].Line, previous_pass.Objects[near_id].Boundary[0].ColumnStart)])){
			        if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	                //Virtual merge into temp object
                     	tmp.clear();
                     	tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     	//Calculate merge object's perimeter
                     	tmp.perimeter = previous_pass.CalculatePerimeter(cur_id, near_id);
                     	//Calculate heterogeneity
                     	h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     	h_smooth = tmp.Smoothness() - previous_pass.Objects[cur_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                     	h_cmpct = tmp.Compactness() - previous_pass.Objects[cur_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                     	h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     	f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                     }
                 	 else{
                     	  //get the heterogeneity already calculated
                    	  f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                     }
   			     	 if(f<best_f){
	     	 			  matches++;
      		     		  best_f = f;
  			   			  best_id = near_id;
         		     }
			     }
   		     }
	 	}
	 	//Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor //Optimize: This loop can be removed for fast mode
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //if neighbor not pushed to priority --> Push
		        if(!previous_pass.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	previous_pass.Objects[near_id].priority_loaded = TRUE;
		        }
		    }
	 	    tmp.clear();
            //merge objects
            tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            //Calculate perimeter
            tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            mergecounter++;
			mseg.last_id++;
            tmp.id = mseg.last_id;
            passn.Objects.insert(make_pair(tmp.id, tmp));//insert here

            previous_pass.Objects[cur_id].merged = TRUE;
            previous_pass.Objects[best_id].merged = TRUE;

            //Update raster
            passn.UpdateRaster(tmp.id);

            tmp.clear();
        }
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));//insert here to
            previous_pass.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            passn.Objects[cur_id].Neighbors.clear();

            //Update raster
            passn.UpdateRaster(cur_id);
        }
    }

    //Create topology --> Load Topology
    cout <<  "Calculating topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;

            cur_id = passn.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = passn.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(passn.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = passn.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(passn.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = passn.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(passn.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = passn.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(passn.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        passn.Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
    cout << "Edge Objects in this pass: " << edgecounter << endl;
    return mergecounter;

}



int FastEdge_Mode::lastpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass){

    int mergecounter = 0;
    int cur_id,near_id, best_id, matches;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    Object tmp;

	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = previous_pass.cycles++;
   	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    passn.raster.push_back(-1);
	}
    //Copy the id number
    passn.HierarchyID = previous_pass.HierarchyID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;

    //Copy segmentation parameters into new level
    passn.ScaleParameter = previous_pass.ScaleParameter;
    passn.Color = previous_pass.Color;
    passn.Compactness = previous_pass.Compactness;
    passn.Edge = previous_pass.Edge;

    //Copy all objects to the new level
    typedef map<int,Object>::const_iterator CI;
	for(CI p=previous_pass.Objects.begin();p!=previous_pass.Objects.end();++p){
	    cur_id = p->first;
	    passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
        passn.Objects[cur_id].priority_loaded = FALSE;
        passn.Objects[cur_id].Neighbors.clear();
        passn.UpdateRaster(cur_id);
	}

    //Main loop -- Merge the edge objects to it's best fit
    for(CI po=previous_pass.Objects.begin();po!=previous_pass.Objects.end();++po){
	    cur_id = po->first;
	    //If object is on egde
	    if(previous_pass.Objects[cur_id].edge == TRUE){
            matches = 0;
            for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
                near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
                //if neighbor is not edge object then we have a possible match
                if(previous_pass.Objects[near_id].edge == FALSE){
                    matches++;
                    //If not calculated, calculate neighbor object's perimeter
                    if(previous_pass.Objects[near_id].perimeter == 0){
                        previous_pass.CalculatePerimeter(near_id);
                    }

                    //Virtual merge into temp object
                    tmp.clear();
                    tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                    //Calculate merge object's perimeter
                    tmp.perimeter = previous_pass.CalculatePerimeter(cur_id, near_id);
                    //Calculate heterogeneity
                    h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                    h_smooth = tmp.Smoothness() - previous_pass.Objects[cur_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                    h_cmpct = tmp.Compactness() - previous_pass.Objects[cur_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                    h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                    f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);

                    if(matches == 1){
                       best_f = f;
                       best_id = near_id;
                    }
                    else{
                        if(f<best_f){
                           best_f = f;
                           best_id = near_id;
                        }
                    }
                }
            }

            //found the best match? --> merge in raster
            if(matches > 0){
                passn.UpdateRaster(cur_id, best_id);
                mergecounter++;
            }
            else{//did not find a match? --> leave as is (this is not supposed to happen, perhaps a cout should be placed here)
                //cout << "An edge object was not merged" << endl;
            }
	    }
    }

    //Clear new levels objects
    passn.Objects.clear();

    //Load objects from updated raster
    passn.LoadObjects();

    //Calculate updated properties
    passn.LoadObjectSums(Img);
    passn.CalculateNeighbors();

    return mergecounter;
}
//TODO
int FastEdge_Mode::lastpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level){

    return 0;
}
