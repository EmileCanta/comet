//-------------------------------------------------------------------------//
//                                                                         //
//                                Convert.h                                //
//                               Version 1.8                               //
//                        Matthieu Lebois January 2007                     //
//                            with G. Georgiev                             //
//                                                                         //
//  This file contain the definition of the function used to convert all   // 
//  the files coming from Narval                                           //
//                                                                         //
//-------------------------------------------------------------------------//

// Call of all the libraries used in this program
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "TFile.h"
#include "TMath.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TStopwatch.h"


#ifndef CONVERT_H
#define CONVERT_H


// Declaration of the function used to make the conversion
int Coucou(int i);
int convertfile(const char* paramfilename,const char *inputfilename, const char *outputfilename, UChar_t coding = 0, UChar_t swap = 1);

#endif
