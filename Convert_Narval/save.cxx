//-------------------------------------------------------------------------//
//                                                                         //
//                            Convert_Narval.cxx                           //
//                               Version 1.8                               //
//                        Matthieu Lebois November 2007                    //
//                            with G. Georgiev                             //
//                                                                         //
//     This file contain the function of conversion to root file.          //
//                                                                         //
//-------------------------------------------------------------------------//

// Call for standard libraries
#include "TMath.h"

// Call of all the libraries used in this program
#include "Convert.h"
#include "Fonctions_convert.h"

// Definition of I/O function
using namespace std;

//Definition of a debugging function
int Coucou(int i)
{
  cout << "Coucou I'm here" << i << endl;
  return 0;
}

//Definition of the function to reconstitute the time
#define TEMPS_47bits(x,y,z) (Double_t)( ((Double_t)(x & 0x7fff) * 65536.0 * 65536.0) + ((Double_t)y * 65536.0) + (Double_t)z )

// Definition of conversion function
int convertfile(const char* paramfilename,const char *inputfilename, const char *outputfilename,int coding, int swap)
{
  
  // Declaration of variable for the construction of a hit
  UChar_t       gr;           // group
  UChar_t       sl;	      // slot
  UChar_t       cr;	      // crate
  UChar_t       ser;          // service
  UChar_t       mark;         // mark
  UChar_t       ch;	      // channel
  UShort_t      en;	      // energy
  Double_t      tm;	      // time that will also be stored in the tree
  UChar_t marker;             // marker to store in the tree
  UInt_t enrj;                // Energy to store in the tree
  UInt_t codingenable = 0;    // Counter of coding enable.
  
  // Declaration of the different pointers used during the conversion
  char * word1;
  char * word2;
  char * bf16;		      // used to read 16 bit blocks
  char branchname1[32];
  char nomdetector[32];
  //TString *branchdetname;     // Branch names
  
  // Declaration of the others variable used during the conversion
  Double_t memoiretemps = 0;
  unsigned short bf;
  unsigned int   tm1 = 0;              // Middle part of time
  unsigned int   tm2 = 0;              // Lower part of time
  unsigned int   tm3 = 0;              // Higher part of time
  int            chks = 0;             // Checksum (useless data)
  ULong64_t      n = 0;		       // event counter max 4 billion if unsigned long
  Long64_t       nbrevent = 0;         // Counter to know the real numbers of entries in the tree
  int            nbreventfous = 0;     // Counter to count event with no meaning
  int            nbrevent2 = 0;        // Counter to select a good number of event
  int            memory_nbrevent_coding = 0; // variable pour stocker le numero d'evenement physique precedent le coding enable
  int            nbrcodingbuffer = 0;  // Counter to determine how many coding enable there is in a buffer
  int            buffersize = 100000;  // Size of the buffer that will be manipulated
  int            nbrbuffer = 0;        // Countenr on the number of buffer treated
  //  int            tetedoeil = 0;
  int            n0det;
  int            group;
  int            crate;
  int            slot;
  int            codenvu = 0;
  UChar_t        index;
  int            channel;
  int            nbrdetect;
  int            nbrline;
  int            nbrdesordre = 0;
  int            smallnbrcoding = 0;   // Counter to check the number of coding_enable (one/card)
  int            nbrcoding = 0;        // Counter to count the number of coding enable
  int            nbrcard = 0;          // Memorization of the number of cards to check the coding enable
  int            option = 1;           // Option to detect the presence of a coding enable or/and a swap 


  // Dynamic allocation of the pointeurs that will contain the data read from the file
  UChar_t *buffer_marker;
  UInt_t *buffer_enrj;
  Double_t *buffer_tm;
  UChar_t *buffer_idx;
  Int_t *buffer_index;
  Double_t *tampon_tm;
  UInt_t *tampon_enrj;
  UChar_t *tampon_marker; 
  UChar_t *tampon_idx;

  //-----------------------------------------------------------------------------//
  //                                                                             //
  //         Part that is used to define the parameters that belong to the       //
  //         experiment: definition of an index to recognize the differents      //
  //         parameters.                                                         //
  //                                                                             //
  //-----------------------------------------------------------------------------//

  // Definition of a structure to create an tree as an index between 
  // group/slot/channel/crate and a detector's number
  struct detector
  {
   UChar_t  nomdetect;
  }det_num;

  //Declaration of variables to read parameters
  char nomFichierParamLect1[256];
  char nomFichierParamLect2[256];
  char buffer[512];
  sprintf(nomFichierParamLect1,paramfilename);
  sprintf(nomFichierParamLect2,paramfilename);
  
  // Declaration of pointers for the reading of files
  ifstream fich_paramlect1(nomFichierParamLect1, ios::in);
  ifstream fich_paramlect2(nomFichierParamLect2, ios::in);
  
  // Declaration of the tree that will be used to store the parameters
  TTree *input = new TTree("input","Parameters tree");
  
  // First Reading of parameters' file to get the number of detector
  nbrline=0;
  if(fich_paramlect1.good())
    {
      while(!fich_paramlect1.eof())
	{
	  fich_paramlect1.getline(buffer,512);
	  if(strcmp(buffer,"") != 0) nbrline++;
	}
      nbrdetect = nbrline-2;
      
      cout << "There is " << nbrdetect << " detectors in the experiment" << endl;
      
      // Closing of parameter file
      fich_paramlect1.close();
    }
  else
    {
      cout << "Impossible d'ouvrir le fichier "<< nomFichierParamLect1 <<endl;
      system("PAUSE");
      exit(1);
    }

  // Second reading of data to create the index
  nbrline=0;
  if(fich_paramlect2.good())
    {
      while(!fich_paramlect2.eof())
	{
	  fich_paramlect2.getline(buffer,512);
	  nbrline++;
	  
	  if(nbrline == 1)
	    {
	      sscanf(buffer,"%d",&nbrcard);
	      cout << "There is " << nbrcard << " COMET used in this experiment!" << endl;
	    }

	  // Parameters saving
	  if(nbrline > 2)
	    {
	      if(buffer != "")
		{
		  // Decomposition of the buffer
		  sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&n0det,&group,&crate,&slot,&channel);
		  sprintf(branchname1,"idx_detect_%d_%d_%d_%d",group,crate,slot,channel);
		  input->Branch(branchname1,&det_num.nomdetect,"nomdetect/b");
		  det_num.nomdetect=n0det;
		  //cout << "Det_num = " << (int)det_num.nomdetect << endl; 
		  input->GetBranch(branchname1)->Fill();
		}
	    }// end of if(nbrline>1)
	}//end of while(!fich_paramlect2.eof())
    }//end of if(fich_paramlect2.good())
  else
    {
      cout << "Impossible d'ouvrir le fichier "<< nomFichierParamLect2 <<endl;
      system("PAUSE");
      exit(1);
    }
  
  // Closing of the parameters file
  fich_paramlect2.close();

  // Activation of all the branches
  input -> SetBranchStatus("*",1);


  // Memory allocation for a buffer to check the event are in order
  tampon_tm = new Double_t[nbrdetect+1];
  tampon_enrj = new UInt_t[nbrdetect+1];
  tampon_marker = new UChar_t[nbrdetect+1];
  tampon_idx = new UChar_t[nbrdetect+1];
  
  // Initialisation des tampons
  for(int j = 0; j < nbrdetect+1; j++)
    {
     tampon_tm[j] = 0;
     tampon_enrj[j] = 0;
     tampon_marker[j] = 0;
     tampon_idx[j] = 0;
    }
  
  buffer_marker = new UChar_t[buffersize];
  buffer_idx = new UChar_t[buffersize];
  buffer_enrj = new UInt_t[buffersize];
  buffer_tm = new Double_t[buffersize];
  buffer_index = new Int_t[buffersize];


  // Definition 
  if( coding == 0 && swap == 0 ) option = 1;
  else if ( coding == 0 && swap == 1 ) option = 2;
  else if ( coding == 1 && swap == 0 ) option = 3;
  else if ( coding == 1 && swap == 1 ) option = 4;

  // Declaration of the root file used to save data
  //TFile *outputfile = new TFile(outputfilename,"RECREATE");
  TFile outputfile(outputfilename,"RECREATE");

  // Opening data file
  ifstream narval ; 
  cout << " name of file " << inputfilename << endl;
  narval.open(inputfilename, ios::in|ios::binary|ios::ate);

  // Declaration of oak, that is the tree containing all events
  cout<< "Building tree  " <<endl;
  TTree *oak = new TTree("Narval_tree",inputfilename); 

  // Declaration of Branches that will contain the data 
  oak -> Branch("Det_nbr",&index,"index/b");             // Detector number
  oak -> Branch("Energy",&enrj,"enrj/i");                // Energy coded on the channel
  oak -> Branch("Time",&tm,"tm/D");                      // Time of the hit
  oak -> Branch("Marker",&marker,"marker/b");            // Marker of the hit (0 if no marker and 1 if marker)
  
  cout << "Seed burried " << endl;
  cout <<"  and in the memory a tree will grow..."<<endl ; 
  
  // And determine its size
  streampos inputfilesize = narval.tellg();

  //ASP

  //cout<<narval.fail()<< " "<<narval.good()<< " "<<narval.bad()<< " "<<narval.eof()<< endl ;
  
  switch(option)
    {


      //--------------------------------------------//
      //      No coding enable and no swap          //
      //--------------------------------------------//
    case 1:
      cout << "Filesize is " << inputfilesize << " bytes for "<< inputfilename << endl;
      
      // Permission for a root file to exceed 2GB
      if(inputfilesize > 2000000000)
	{
	  oak->SetMaxTreeSize(inputfilesize); 
	}
      cout << "The root file size will be of ~ " << (Double_t)inputfilesize/2000000000 << " Gb"<< endl;
      
      narval.seekg(0,ios::beg);
      if (inputfilesize < 0 ) 
	{
	  cout <<"$$$"<<endl <<" $$$ Problem with file "<<inputfilename<<endl ; 
	  cout << "$$$ skiping it"<<endl<<"$$$ "<<endl ; 
	  break;
	}

      // Skip header till we find an event 0xffff 0x8
      word1 = new char [2];
      word2 = new char [2];
      bf16 = new char [2];
      narval.read(word1,2);
      narval.read(word2,2);
  
      while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
	{
	  // Move one word backwards but read two
	  narval.seekg(-2,ios::cur);
	  if (narval.good()) narval.read(word1,2);
	  if (narval.good()) narval.read(word2,2);
	  
	  // Test to check the end of file (corrupted file)
	  if (narval.eof())
	    {
	      cout << "Fatal error, file does not start with detector event" << endl;
	      return n;
	    }// end of if (narval.eof())
	}// end of while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
      
      // Writting of the declaration for the good start of conversion
      cout << "Start conversion at position " << narval.tellg()-(streampos)4<< endl;
      cout << "-------------------------------------------------"<< endl; 
      
      // Scan of all the data file
      while (narval.good())
	{      
	  // Check for a physical event: detector event 0xffff 0x0008
	  while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	    {
	      // Increment of any event counter
	      n++;
	      
	      // Extraction of hit description
	      if (narval.good()) narval.read(bf16,2);        
	      
	      // Read group, slot, crate, ch, marker and service
	      bf = word(bf16[0],bf16[1]);
	      gr = (UChar_t)(bf >> 12);
	      sl = (UChar_t)((bf >> 8) & 15);
	      cr = (UChar_t)((bf >> 6) & 3);
	      ser = (UChar_t)((bf >> 4 ) & 2);
	      mark = (UChar_t)((bf >> 4 ) & 1);
	      ch = (UChar_t)(bf & 15);
	  
	      // Creation of the detector number to find it in the index
	      sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
	      
	      // Read energy
	      if (narval.good()) narval.read(bf16,2);
	      en = (UInt_t)word(bf16[0],bf16[1]);
	      
	      // Read time
	      if (narval.good()) narval.read(bf16,2);
	      tm1 = (word(bf16[0],bf16[1]));
	      if (narval.good()) narval.read(bf16,2);
	      tm2 = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      chks = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      tm3 = (word(bf16[0],bf16[1]));
	      
	      // Time reconstruction to be able to use it for checking
	      tm = TEMPS_47bits(tm3,tm1,tm2);

	      // Filling of the tree used for storing PHYSICAL hits
	      if(ser == 0) // checking it is no a service hit 
		{
		  // Reading of the input tree to get the information concerning the detector number 
		  TBranch *idx = input->GetBranch(nomdetector); 
		  
		  // Check for the quality of the information read
		  if(idx) 
		    { 
		      // Reading of the index to get the detector number corresponding to a TClonesArray
		      idx->SetAddress(&det_num); //##??
		      idx->GetEntry(0);
		      index = det_num.nomdetect;
		  
		      // Verification of the order of the event by detector to check hit are stored
		      // in chronological order
		      if((tm*400 - memoiretemps) >= 0)
			{
			  // Filling of the structure to prepare the storing of the event
			  marker = mark;
			  enrj = en;
			  tm = tm*400; // to be in unit of ps
			 			  
			  // Filling of the tree
			  oak -> Fill();
			  
			  // Augmentation of hit counter
			  nbrevent++;
			  			  
			  // Memorization of the time for this detector to test the order of events 
			  // coming from the same detector
			  memoiretemps = (Double_t)tm;
			  
			}// End of if(tm >= memoiretemps[index])
		      else
			{
			  cout << "An event is not in order in detector " << (int)index << endl;
			  cout << "Present Time = " << tm << endl;
			  cout << "Referred to " << memoiretemps << endl;
			  memoiretemps = 0;
			  nbrdesordre++;
			}
		      
		    }// End of if(idx)
	      
		  else
		    {
		      cout << "Foolish event detected: detector not existing" << endl;
		      cout << "group : " << (int)gr << endl;
		      cout << "slot : " << (int)sl << endl;
		      cout << "crate : " << (int)cr << endl;
		      cout << "service : " << (int)ser << endl;
		      cout << "marker : " << (int)mark << endl;
		      cout << "channel : " << (int)ch << endl;
		      nbreventfous++;
		    }
		}// End of if(ser == 0)
	  
	      // Read the next word to able to check if it is a physical hit or not
	      if (narval.good()) narval.read(word1,2);
	      if (narval.good()) narval.read(word2,2);
	      
	      // In case of incomplete event, eof is reached unexpectedly
	      if (narval.eof())
		{
		  cout << "Irregular end of file, incomplete event" << endl;
		  word1[0]=0xff;
		  word1[1]=0xff;
		  word2[0]=0xff;
		  word2[1]=0xff;                              // set marker for eof
		}// End of if (narval.eof())
	    }// End of while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	  
	  // If file contains loose zeroes in the middle, mark it as an internal header
	  if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	    {
	      cout << "Irregular file structure,loose zeroes found at "<< narval.tellg()-(streampos)4<< endl;
	      // narval.seekg(0,ios::end);                  // move to eof without provoking eof yet
	      // narval.read(word2,2);                      // actually provoke eof
	      word1[0]=0xff;
	      word1[1]=0xff;
	      word2[0]=0xff;
	      word2[1]=0xff;                                 // set marker for eof
	    }// End of if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	  
	  // This is an internal header or regular eof or eof incomplete event or loose zeroes
	  else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	    {
	      if (narval.good()) narval.read(word1, 2);
	      if (narval.good()) narval.read(word2, 2);
	      
	      // execute the loop as long as we don't encounter a detector event
	      while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd))))
		{
		  narval.seekg(-2,ios::cur);
		  if(narval.good()) narval.read(word1,2);
		  if(narval.good()) narval.read(word2,2);
		  if(narval.eof()) 
		    { 
		      cout << "End of file detected" << endl;
		      // not very elegant but we need to escape from the while loop
		      word1[0]=0xff;
		      word1[1]=0xff;
		      word2[0]=0x8;
		      word2[1]=0x0;
		    }// End of if(narval.eof())
		}// End of while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd)))) 
	    }// End of else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	  
	  else
	    {
	      cout << "Fatal error, irregular event structure "<< endl;
	      return n;
	    }     	
	}// End of while for eof
      
      // Ending of the conversion procedure
      if(narval.eof())
	{
	  // Writing the tree in the root file to be sure the root file is well closed
	  cout << "Saving the tree in the file" << endl;
	  oak -> Write();
	  
	  // Unallocation of the memory used for the buffer
	  delete [] bf16;
	  delete [] word1;
	  delete [] word2;
	  //delete [] branchdetname;     
	    
	  // Closing of all the file used during the conversion process
	  narval.close();
	  //outputfile -> Close();
	  outputfile.Close();

	  
	  // Writting of a summary of the conversion process 
	  cout << "This makes "<< n << " hits read in the file"<<endl;
	  cout << "This makes "<< n-nbrevent << " service hits read in the file"<<endl;
	  cout << "This makes "<< nbreventfous << " completely foolish events read in the file"<<endl;
	  cout << "This makes "<< nbrdesordre << " hits that were not in order (in time)"<<endl;
	  cout << "This makes "<< nbrevent << " physical hits read in the file"<<endl;
	  cout << endl;

	}// End of if(narval.eof())
      
      else
	{
	  cout << "Error exiting program" << endl;   
	}
      
      break;

      
      //--------------------------------------------//
      //      No coding enable and swap             //
      //--------------------------------------------//
    case 2:
      cout << "Filesize is " << inputfilesize << " bytes for "<< inputfilename << endl;
      
      // Permission for a root file to exceed 2GB
      if(inputfilesize > 2000000000)
	{
	  oak->SetMaxTreeSize(inputfilesize); 
	}
      cout << "The root file size will be of ~ " << (Double_t)inputfilesize/2000000000 << " Gb"<< endl;
      
      narval.seekg(0,ios::beg);
      if (inputfilesize < 0 ) 
	{
	  cout <<"$$$"<<endl <<" $$$ Problem with file "<<inputfilename<<endl ; 
	  cout << "$$$ skiping it"<<endl<<"$$$ "<<endl ; 
	  break;
	}

      // Skip header till we find an event 0xffff 0x8
      word1 = new char [2];
      word2 = new char [2];
      bf16 = new char [2];
      narval.read(word1,2);
      narval.read(word2,2);
  
      while ((word_swap(word1[0],word1[1]) != 0xffff)||(word_swap(word2[0],word2[1]) !=0x8))
	{
	  // Move one word backwards but read two
	  narval.seekg(-2,ios::cur);
	  if (narval.good()) narval.read(word1,2);
	  if (narval.good()) narval.read(word2,2);
	  
	  // Test to check the end of file (corrupted file)
	  if (narval.eof())
	    {
	      cout << "Fatal error, file does not start with detector event" << endl;
	      return n;
	    }// end of if (narval.eof())
	}// end of while ((word_swap(word1[0],word1[1]) != 0xffff)||(word_swap(word2[0],word2[1]) !=0x8))
      
      // Writting of the declaration for the good start of conversion
      cout << "Start conversion at position " << narval.tellg()-(streampos)4<< endl;
      cout << "-------------------------------------------------"<< endl; 
      
      // Scan of all the data file
      while (narval.good())
	{      
	  // Check for a physical event: detector event 0xffff 0x0008
	  while ((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0x8))
	    {
	      // Increment of any event counter
	      n++;
	      
	      // Extraction of hit description
	      if (narval.good()) narval.read(bf16,2);        
	      
	      // Read group, slot, crate, ch, marker and service
	      bf = word_swap(bf16[0],bf16[1]);
	      gr = (UChar_t)(bf >> 12);
	      sl = (UChar_t)((bf >> 8) & 15);
	      cr = (UChar_t)((bf >> 6) & 3);
	      ser = (UChar_t)((bf >> 4 ) & 2);
	      mark = (UChar_t)((bf >> 4 ) & 1);
	      ch = (UChar_t)(bf & 15);
	  
	      // Creation of the detector number to find it in the index
	      sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
	      
	      // Read energy
	      if (narval.good()) narval.read(bf16,2);
	      en = (UInt_t)word_swap(bf16[0],bf16[1]);
	      
	      // Read time
	      if (narval.good()) narval.read(bf16,2);
	      tm1 = (word_swap(bf16[0],bf16[1]));
	      if (narval.good()) narval.read(bf16,2);
	      tm2 = word_swap(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      chks = word_swap(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      tm3 = (word_swap(bf16[0],bf16[1]));
	      
	      // Time reconstruction to be able to use it for checking
	      tm = TEMPS_47bits(tm3,tm1,tm2);

	      // Filling of the tree used for storing PHYSICAL hits
	      if(ser == 0) // checking it is no a service hit 
		{
		  // Reading of the input tree to get the information concerning the detector number 
		  TBranch *idx = input->GetBranch(nomdetector); 
		  
		  // Check for the quality of the information read
		  if(idx) 
		    { 
		      // Reading of the index to get the detector number corresponding to a TClonesArray
		      idx->SetAddress(&det_num); //##??
		      idx->GetEntry(0);
		      index = det_num.nomdetect;
		  
		      // Verification of the order of the event by detector to check hit are stored
		      // in chronological order
		      if((tm*400 - memoiretemps) >= 0)
			{
			  // Filling of the structure to prepare the storing of the event
			  marker = mark;
			  enrj = en;
			  tm = tm*400; // to be in unit of ps
			 			  
			  // Filling of the tree
			  oak -> Fill();
			  
			  // Augmentation of hit counter
			  nbrevent++;
			  			  
			  // Memorization of the time for this detector to test the order of events 
			  // coming from the same detector
			  memoiretemps = (Double_t)tm;
			  
			}// End of if(tm >= memoiretemps[index])
		      else
			{
			  cout << "An event is not in order in detector " << (int)index << endl;
			  cout << "Present Time = " << tm << endl;
			  cout << "Referred to " << memoiretemps << endl;
			  memoiretemps = 0;
			  nbrdesordre++;
			}
		      
		    }// End of if(idx)
	      
		  else
		    {
		      cout << "Foolish event detected: detector not existing" << endl;
		      cout << "group : " << (int)gr << endl;
		      cout << "slot : " << (int)sl << endl;
		      cout << "crate : " << (int)cr << endl;
		      cout << "service : " << (int)ser << endl;
		      cout << "marker : " << (int)mark << endl;
		      cout << "channel : " << (int)ch << endl;
		      nbreventfous++;
		    }
		}// End of if(ser == 0)
	  
	      // Read the next word to able to check if it is a physical hit or not
	      if (narval.good()) narval.read(word1,2);
	      if (narval.good()) narval.read(word2,2);
	      
	      // In case of incomplete event, eof is reached unexpectedly
	      if (narval.eof())
		{
		  cout << "Irregular end of file, incomplete event" << endl;
		  word1[0]=0xff;
		  word1[1]=0xff;
		  word2[0]=0xff;
		  word2[1]=0xff;                              // set marker for eof
		}// End of if (narval.eof())
	    }// End of while ((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0x8))
	  
	  // If file contains loose zeroes in the middle, mark it as an internal header
	  if((word_swap(word1[0],word1[1])==0x0) && (word_swap(word2[0],word2[1])==0x0))
	    {
	      cout << "Irregular file structure,loose zeroes found at "<< narval.tellg()-(streampos)4<< endl;
	      // narval.seekg(0,ios::end);                  // move to eof without provoking eof yet
	      // narval.read(word2,2);                      // actually provoke eof
	      word1[0]=0xff;
	      word1[1]=0xff;
	      word2[0]=0xff;
	      word2[1]=0xff;                                 // set marker for eof
	    }// End of if((word_swap(word1[0],word1[1])==0x0) && (word_swap(word2[0],word2[1])==0x0))
	  
	  // This is an internal header or regular eof or eof incomplete event or loose zeroes
	  else if((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0xffff))
	    {
	      if (narval.good()) narval.read(word1, 2);
	      if (narval.good()) narval.read(word2, 2);
	      
	      // execute the loop as long as we don't encounter a detector event
	      while(!((word_swap(word1[0],word1[1])==0xffff) && ( (word_swap(word2[0],word2[1])==0x8) || (word_swap(word2[0],word2[1])==0xd))))
		{
		  narval.seekg(-2,ios::cur);
		  if(narval.good()) narval.read(word1,2);
		  if(narval.good()) narval.read(word2,2);
		  if(narval.eof()) 
		    { 
		      cout << "End of file detected" << endl;
		      // not very elegant but we need to escape from the while loop
		      word1[0]=0xff;
		      word1[1]=0xff;
		      word2[0]=0x8;
		      word2[1]=0x0;
		    }// End of if(narval.eof())
		}// End of while(!((word_swap(word1[0],word1[1])==0xffff) && ( (word_swap(word2[0],word2[1])==0x8) || (word_swap(word2[0],word2[1])==0xd)))) 
	    }// End of else if((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0xffff))
	  
	  else
	    {
	      cout << "Fatal error, irregular event structure "<< endl;
	      return n;
	    }     	
	}// End of while for eof
      
      // Ending of the conversion procedure
      if(narval.eof())
	{
	  // Writing the tree in the root file to be sure the root file is well closed
	  cout << "Saving the tree in the file" << endl;
	  oak -> Write();
	  
	  // Unallocation of the memory used for the buffer
	  delete [] bf16;
	  delete [] word1;
	  delete [] word2;
	  //delete [] branchdetname;     
	    
	  // Closing of all the file used during the conversion process
	  narval.close();
	  outputfile.Close();
	  
	  // Writting of a summary of the conversion process 
	  cout << "This makes "<< n << " hits read in the file" << endl;
	  cout << "This makes "<< n-nbrevent << " service hits read in the file" << endl;
	  cout << "This makes "<< nbreventfous << " completely foolish events read in the file" << endl;
	  cout << "This makes "<< nbrdesordre << " hits that were not in order (in time)" << endl;
	  cout << "This makes "<< nbrevent << " physical hits read in the file" << endl;
	  cout << endl;

	}// End of if(narval.eof())
      
      else
	{
	  cout << "Error exiting program" << endl;   
	}

      break;







      //--------------------------------------------//
      //      Coding enable and no swap             //
      //--------------------------------------------//
    case 3:
      // Creation of a new branch to store the number of coding enable
      // This could be used during the search of coincidences as a criteria to avoid coincidences between
      // two different coding enable
      oak -> Branch("Coding", &codingenable, "codingenable/i");
      
      cout << "Filesize is " << inputfilesize << " bytes for "<< inputfilename << endl;
      
      // Permission for a root file to exceed 2GB
      if(inputfilesize > 2000000000)
	{
	  oak->SetMaxTreeSize(inputfilesize); 
	}
      cout << "The root file size will be of ~ " << (Double_t)inputfilesize/2000000000 << " Gb"<< endl;
      
      narval.seekg(0,ios::beg);
      if (inputfilesize < 0 ) 
	{
	  cout <<"$$$"<<endl <<" $$$ Problem with file "<<inputfilename<<endl ; 
	  cout << "$$$ skiping it"<<endl<<"$$$ "<<endl ; 
	  break;
	}
            
      // Skip header till we find an event 0xffff 0x8
      word1 = new char [2];
      word2 = new char [2];
      bf16 = new char [2];
      narval.read(word1,2);
      narval.read(word2,2);
  
      while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
	{
	  // Move one word backwards but read two
	  narval.seekg(-2,ios::cur);
	  if (narval.good()) narval.read(word1,2);
	  if (narval.good()) narval.read(word2,2);

	  // Test to check the end of file (corrupted file)
	  if (narval.eof())
	    {
	      cout << "Fatal error, file does not start with detector event" << endl;
	      return n;
	    }// end of if (narval.eof())
	}// end of while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
      
      // Writting of the declaration for the good start of conversion
      cout << "Start conversion at position " << narval.tellg()-(streampos)4<< endl;
      cout << "-------------------------------------------------"<< endl; 
      
      // Scan of all the data file
      while (narval.good())
	{      
	  // Check for a physical event: detector event 0xffff 0x0008
	  while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	    {
	      // Increment of any event counter
	      n++;
	      
	      // Extraction of hit description
	      if (narval.good()) narval.read(bf16,2);        
	      
	      // Read group, slot, crate, ch, marker and service
	      bf = word(bf16[0],bf16[1]);
	      gr = (UChar_t)(bf >> 12);
	      sl = (UChar_t)((bf >> 8) & 15);
	      cr = (UChar_t)((bf >> 6) & 3);
	      ser = (UChar_t)((bf >> 4 ) & 2);
	      mark = (UChar_t)((bf >> 4 ) & 1);
	      ch = (UChar_t)(bf & 15);
	  	      
	      // Read energy
	      if (narval.good()) narval.read(bf16,2);
	      en = (UInt_t)word(bf16[0],bf16[1]);
	      
	      // Read time
	      if (narval.good()) narval.read(bf16,2);
	      tm1 = (word(bf16[0],bf16[1]));
	      if (narval.good()) narval.read(bf16,2);
	      tm2 = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      chks = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      tm3 = (word(bf16[0],bf16[1]));
	      
	      // Time reconstruction to be able to use it for checking
	      tm = TEMPS_47bits(tm3,tm1,tm2);

	      // Verification of the presence of a coding enable to reset the memory of time
	      // A coding enable event is define as a service event on a channel 6 (unexisting channel)
	      if(ser == 2 && ch == 6)//(ser == 1 && mark == 1 && ch == 5)
		{
		  smallnbrcoding++;
		  if(smallnbrcoding == nbrcard)
		    {
		      smallnbrcoding = 0;
		      nbrcoding++;
		      //codingenable = nbrcoding;
		      memory_nbrevent_coding = nbrevent2;
		      codenvu = 1;
		    }// End of if(smallnbrcoding == nbrcard)
		}// End of if(ser == 1 && mark == 1 && ch == 5)

	      // Filling of the tree used for storing PHYSICAL hits
	      if(ser == 0) // checking it is no a service hit 
		{
		  // Creation of the detector number to find it in the index (x2 pour etre tranquille)
		  sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
		  sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
		  
		  // Reading of the input tree to get the information concerning the detector number 
		  TBranch *idx = input->GetBranch(nomdetector); 

		  // Check for the quality of the information read
		  if(idx) 
		    { 
		      // Reading of the index to get the detector number corresponding to a TClonesArray
		      idx -> SetAddress(&det_num); //##??
		      idx -> GetEntry(0);
		      buffer_idx[nbrevent2] = det_num.nomdetect;

		      // Filling of the structure to prepare the storing of the event
		      buffer_marker[nbrevent2] = mark;
		      buffer_enrj[nbrevent2] = en;
		      buffer_tm[nbrevent2] = tm*400; // to be in unit of ps

		      // Augmentation of hit counter
		      nbrevent++;
		      nbrevent2++;
		      		      
		      // Je m'arrete nbrdetect evenements apres le coding enable pour etre sur de n'avoir
		      // qu'un coding enable par buffer
		      if((nbrevent2 ==  (memory_nbrevent_coding + nbrdetect)) && codenvu == 1)
			{
			  
			  // Le premier buffer est precede par un coden...
			  // Je dois donc trier un buffer de nbrdetect evenements
			  if(codingenable == 0)
			    {
			      //cout << "Je trie le premier buffer"<<endl;

			      //first_tribuffer_temps(buffer_tm, buffer_index, tampon_tm, nbrevent2, nbrdetect);
			      first_tri_buffers(tampon_tm,buffer_tm,tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, nbrevent2, nbrdetect);
			      nbrbuffer++;

			      // Now the buffer is sorted I will be able to store it in the TTree
			      for(int ll = 0; ll < nbrevent2; ll++)
				{
				  index = tampon_idx[ll];
				  marker = tampon_marker[ll];
				  enrj = tampon_enrj[ll];
				  tm = tampon_tm[ll];

				  

				  // Pour ce premier buffer on a fait que remplir les tampons...
				}

			      codingenable = 1;			      
			      nbrevent2 = 0;
			      memory_nbrevent_coding = -666;
			      codenvu = 0;
			      
			    }
			  else
			    {
			      //cout << "Je trie un autre buffer que le premier" << endl;
			      nbrbuffer++;

			      // Je faire le tri sur le "buffer taille reduite"
			      tribuffer_temps_coding(buffer_tm, tampon_tm, tampon_enrj, buffer_enrj, tampon_marker, buffer_marker,tampon_idx, buffer_idx, buffer_index, nbrcodingbuffer,nbrevent2,nbrdetect);
			      
			      // Now the buffer is sorted I will be able to store it in the TTree
			      for(int ll = 1; ll < nbrevent2; ll++)
				{
				  //if((buffer_tm[ll-1]- buffer_tm[ll]) > 1.e5) codingenable++;
				  index = buffer_idx[ll];
				  marker = buffer_marker[ll];
				  enrj = buffer_enrj[ll];
				  tm = buffer_tm[ll];
				  oak -> Fill();
				  oak -> Print();
				}
			      
			      // J'ai eu un coding enable dans ce buffer... donc j'incremente le conteur
			      codingenable++;
			      			      			      
			      // Je remets à zero le nombre d'evenements pour remplir un nouveau buffer de taille normale
			      nbrevent2 = 0;
			      memory_nbrevent_coding = -666;
			      codenvu = 0;
			    }
			}// Fin du if dans le cas d'un coden dans le buffer
			  
		      // Dans le cas d'un buffer sans coden je fais un tri "standard"
		      if(nbrevent2 == buffersize && codenvu != 1)
			{
			  // Tri des buffers
			  tribuffer_temps(buffer_tm, buffer_index, tampon_tm, buffersize, nbrdetect);
			  tri_buffers(tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, buffersize, nbrdetect);
			  //cout << "On en est au buffer " << nbrbuffer << endl;
			  nbrbuffer++;
			
			  // Now the buffer is sorted I will be able to store it in the TTree
			  for(int ll = 0; ll < buffersize; ll++)
			    {
			      //if(nbrevent > 3000000)Coucou(ll);
			      index = buffer_idx[ll];
			      marker = buffer_marker[ll];
			      enrj = buffer_enrj[ll];
			      tm = buffer_tm[ll];
			      oak -> Fill();
			      oak -> Print();
			    }
			  
			  nbrevent2 = 0;
			  nbrcodingbuffer = 0;
			  memory_nbrevent_coding = -666;
			}// End of if(nbrevent2 == buffersize)
		    }// End of if(idx)
		  
		  else
		    {
		      cout << "Foolish event detected: detector not existing" << endl;
		      cout << "group : " << (int)gr << endl;
		      cout << "slot : " << (int)sl << endl;
		      cout << "crate : " << (int)cr << endl;
		      cout << "service : " << (int)ser << endl;
		      cout << "marker : " << (int)mark << endl;
		      cout << "channel : " << (int)ch << endl;
		      cout << "It is the event number : " << nbrevent2 << endl;
		      cout << "Of the " << nbrbuffer << "e buffer" << endl;
		      nbreventfous++;
		    }
		}// End of if(ser == 0)
	  
	      // Read the next word to able to check if it is a physical hit or not
	      if (narval.good()) narval.read(word1,2);
	      if (narval.good()) narval.read(word2,2);
	      
	      // In case of incomplete event, eof is reached unexpectedly
	      if (narval.eof())
		{
		  cout << "Irregular end of file, incomplete event" << endl;
		  word1[0]=0xff;
		  word1[1]=0xff;
		  word2[0]=0xff;
		  word2[1]=0xff;                              // set marker for eof
		}// End of if (narval.eof())
	    }// End of while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	  
	  // If file contains loose zeroes in the middle, mark it as an internal header
	  if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	    {
	      cout << "Irregular file structure,loose zeroes found at "<< narval.tellg()-(streampos)4<< endl;
	      // narval.seekg(0,ios::end);                  // move to eof without provoking eof yet
	      // narval.read(word2,2);                      // actually provoke eof
	      word1[0]=0xff;
	      word1[1]=0xff;
	      word2[0]=0xff;
	      word2[1]=0xff;                                 // set marker for eof
	    }// End of if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	  
	  // This is an internal header or regular eof or eof incomplete event or loose zeroes
	  else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	    {
	      if (narval.good()) narval.read(word1, 2);
	      if (narval.good()) narval.read(word2, 2);
	      
	      // execute the loop as long as we don't encounter a detector event
	      while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd))))
		{
		  narval.seekg(-2,ios::cur);
		  if(narval.good()) narval.read(word1,2);
		  if(narval.good()) narval.read(word2,2);
		  if(narval.eof()) 
		    { 
		      cout << "End of file detected" << endl;
		      // not very elegant but we need to escape from the while loop
		      word1[0]=0xff;
		      word1[1]=0xff;
		      word2[0]=0x8;
		      word2[1]=0x0;
		    }// End of if(narval.eof())
		}// End of while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd)))) 
	    }// End of else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	  
	  else
	    {
	      cout << "Fatal error, irregular event structure "<< endl;
	      return n;
	    }     	
	}// End of while for eof
      
      // Ending of the conversion procedure
      if(narval.eof())
	{

	  // Organization of the last buffer that might have a size < buffersize
	  buffersize = nbrevent2;
	 	  
	  tribuffer_temps(buffer_tm, buffer_index, tampon_tm, buffersize, nbrdetect);
	  tri_buffers(tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, buffersize, nbrdetect);
	  
	  nbrbuffer++;

	  // Now the buffer is sorted I will be able to store it in the TTree
	  for(int ll = 0; ll < buffersize; ll++)
	    {
	      index = buffer_idx[ll];
	      marker = buffer_marker[ll];
	      enrj = buffer_enrj[ll];
	      tm = buffer_tm[ll];
	      oak -> Fill();
	    }
	  // Storage of the "tampon"
	  for(int ll = 0; ll < nbrdetect; ll++)
	    {
	      index = tampon_idx[ll];
	      marker = tampon_marker[ll];
	      enrj = tampon_enrj[ll];
	      tm = tampon_tm[ll];
	      oak -> Fill(); 
	    }
	  

	  // Writing the tree in the root file to be sure the root file is well closed
	  cout << "Saving the tree in the file" << endl;
	  oak -> Write();
	  
	  // Unallocation of the memory used for the buffer
	  delete [] bf16;
	  delete [] word1;
	  delete [] word2;
	  //delete [] branchdetname;     
	    
	  // Closing of all the file used during the conversion process
	  narval.close();
	  outputfile.Close();

	  
	  // Writting of a summary of the conversion process 
	  cout << "This makes " << nbrcoding << " coden detected in the file" << endl;
	  cout << "This makes " << nbrbuffer << " buffer treated in the file" << endl;
	  cout << "This makes " << n << " hits read in the file" << endl;
	  cout << "This makes " << n-nbrevent << " service hits read in the file" << endl;
	  cout << "This makes " << nbreventfous << " completely foolish events read in the file" << endl;
	  cout << "This makes " << nbrdesordre << " hits that were not in order (in time)" << endl;
	  cout << "This makes " << nbrevent << " physical hits read in the file"<< endl;
	  cout << endl;

	}// End of if(narval.eof())
      
      else
	{
	  cout << "Error exiting program" << endl;   
	}
      
      // Effacement des differents tableaux  
      delete buffer_marker;
      delete buffer_enrj;
      delete buffer_tm;
      delete buffer_idx;
      delete buffer_index;
      delete tampon_tm;
      delete tampon_enrj;
      delete tampon_marker; 
      delete tampon_idx;



      break;



      //--------------------------------------------//
      //         Coding enable and swap        $$$  //
      //--------------------------------------------//
    case 4:
      // Creation of a new branch to store the number of coding enable
      // This could be used during the search of coincidences as a criteria to avoid coincidences between
      // two different coding enable
  
      oak -> Branch("Coding", &codingenable, "codingenable/i");
      
      cout << "Filesize is " << inputfilesize << " bytes for "<< inputfilename << endl;
      
      // Permission for a root file to exceed 2GB
      if(inputfilesize > 2000000000)
	{
	  oak->SetMaxTreeSize(inputfilesize); 
	}
      cout << "The root file size will be of ~ " << (Double_t)inputfilesize/2000000000 << " Gb"<< endl;
      
      narval.seekg(0,ios::beg);
      if (inputfilesize < 0 ) 
	{
	  cout <<"$$$"<<endl <<" $$$ Problem with file "<<inputfilename<<endl ; 
	  cout << "$$$ skiping it"<<endl<<"$$$ "<<endl ; 
	  break;
	}

      // Skip header till we find an event 0xffff 0x8
      word1 = new char [2];
      word2 = new char [2];
      bf16 = new char [2];
      narval.read(word1,2);
      narval.read(word2,2);
  
      while ((word_swap(word1[0],word1[1]) != 0xffff)||(word_swap(word2[0],word2[1]) !=0x8))
	{

	  // Move one word backwards but read two
	  narval.seekg(-2,ios::cur);
	 
	  if (narval.good()) narval.read(word1,2);
	  if (narval.good()) narval.read(word2,2);
	  
	  // Test to check the end of file (corrupted file)
	  if (narval.eof())
	    {
	      cout << "Fatal error, file does not start with detector event" << endl;
	      return n;
	    }// end of if (narval.eof())
	}// end of while ((word_swap(word1[0],word1[1]) != 0xffff)||(word_swap(word2[0],word2[1]) !=0x8))
      
      // Writting of the declaration for the good start of conversion
      cout << "Start conversion at position " << narval.tellg()-(streampos)4<< endl;
      cout << "-------------------------------------------------"<< endl; 
     

      //cout<<narval.fail()<< " "<<narval.good()<< " "<<narval.bad()<< " "<<narval.eof()<< endl ;
 
      // Scan of all the data file
      while (narval.good())
	{      
	  // Check for a physical event: detector event 0xffff 0x0008
	  while ((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0x8))
	    {
	      // Increment of any event counter
	      n++;

	      
	      // Extraction of hit description
	      if (narval.good()) narval.read(bf16,2);        
	      
	      // Read group, slot, crate, ch, marker and service
	      bf = word_swap(bf16[0],bf16[1]);
	      gr = (UChar_t)(bf >> 12);
	      sl = (UChar_t)((bf >> 8) & 15);
	      cr = (UChar_t)((bf >> 6) & 3);
	      ser = (UChar_t)((bf >> 4 ) & 2);
	      mark = (UChar_t)((bf >> 4 ) & 1);
	      ch = (UChar_t)(bf & 15);
	  
	      

	      // Creation of the detector number to find it in the index
	      sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
	      
	      // Read energy
	      if (narval.good()) narval.read(bf16,2);
	      en = (UInt_t)word_swap(bf16[0],bf16[1]);
	      
 
	      // Read time
	      if (narval.good()) narval.read(bf16,2);
	      tm1 = (word_swap(bf16[0],bf16[1]));
	      if (narval.good()) narval.read(bf16,2);
	      tm2 = word_swap(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      chks = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      tm3 = (word_swap(bf16[0],bf16[1]));
	      
	      // Time reconstruction to be able to use it for checking
	      tm = TEMPS_47bits(tm3,tm1,tm2);

 
	      // Verification of the presence of a coding enable to reset the memory of time
	      // A coding enable event is define as a service event on a channel 6 (unexisting channel)
	      if(ser == 2 && ch == 6)//(ser == 1 && mark == 1 && ch == 5)
		{
		  smallnbrcoding++;
		 
		  if(smallnbrcoding == nbrcard)
		    {
		      /*
		      cout << endl;
		      cout << " Coucou j'ai un coden " << endl;
		      cout << nbrevent << endl;
		      */
		      smallnbrcoding = 0;
		      nbrcoding++;
		      //codingenable = nbrcoding;
		      memory_nbrevent_coding = nbrevent2;
		      codenvu = 1;
		    }// End of if(smallnbrcoding == nbrcard)
		}// End of if(ser == 1 && mark == 1 && ch == 5)

	      // Filling of the tree used for storing PHYSICAL hits
	      if(ser == 0) // checking it is no a service hit 
		{
		  // Creation of the detector number to find it in the index (x2 pour etre tranquille)
		  sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
		  sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
		  
		  // Reading of the input tree to get the information concerning the detector number 
		  TBranch *idx = input->GetBranch(nomdetector); 
		  
		  // Check for the quality of the information read
		  if(idx) 
		    { 

		      // Reading of the index to get the detector number corresponding to a TClonesArray
		      idx -> SetAddress(&det_num); //##??

		      idx -> GetEntry(0);

		      buffer_idx[nbrevent2] = det_num.nomdetect;

		      // Filling of the structure to prepare the storing of the event
		      buffer_marker[nbrevent2] = mark;
		      buffer_enrj[nbrevent2] = en;
		      buffer_tm[nbrevent2] = tm*400; // to be in unit of ps

		     
		      if(nbrbuffer < 2)
			{
			  cout << " Index : " << (int)buffer_idx[nbrevent2] << "; Marker : " << (int)buffer_marker[nbrevent2] << "; Nrj : " << (int)buffer_enrj[nbrevent2] << "; Time : " << (double)buffer_tm[nbrevent2] << "; Coden : " << (int)codingenable << endl;
			  
			}
		      
		      
			

		      // Augmentation of hit counter
		      nbrevent++;
		      nbrevent2++;
		      		      
		      //Coucou(nbrevent2);

		      // Dans le cas du traitement des tout premiers evenements
		      // Je genere un buffer de la taille du tampon pour remplir ce dernier
		      // Le premier buffer est precede par un coden...
		      // Je doix donc trier un buffer de nbrdetect evenements
		      if(nbrevent2 == nbrdetect && codingenable == 0)
			{
			  //cout << "Je trie le premier buffer"<<endl;
			  
			  // Fonction pour trier le buffer d'evenements
			  first_tri_buffers(tampon_tm,buffer_tm,tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, nbrevent2, nbrdetect);
			  nbrbuffer++;
			  
			  if(nbrbuffer < 2)
			    {
			      for(int k=0; k < 4; k++) cout << endl;
			      
			      cout << "*************************************" << endl;
			      cout << "J'affiche le tout premier tampon "<< endl;
			    }

			  // Now the bufferis sorted I will be able to store it in the TTree
			  for(int ll = 0; ll < nbrevent2; ll++)
			    {
			      index = tampon_idx[ll];
			      marker = tampon_marker[ll];
			      enrj = tampon_enrj[ll];
			      tm = tampon_tm[ll];
			      
			      
			      //if(nbrbuffer < 2) cout << " Index : " << (int)index << "; Marker : " << (int)marker << "; Nrj : " << (int)enrj << "; Time : " << (double)tm << "; Coden : " << (int)codingenable << endl;
			     
			      // Pour ce premier buffer on a fait que remplir les tampons...
			    } 
			  
			  if(nbrbuffer < 2) cout << endl;

			  codingenable = 1;			      
			  nbrevent2 = 0;
			  memory_nbrevent_coding = -666;
			  codenvu = 0;
			}
		      
		      // Je m'arrete nbrdetect evenements apres le coding enable pour etre sur de n'avoir
		      // qu'un coding enable par buffer
		      if((nbrevent2 ==  (memory_nbrevent_coding + nbrdetect)) && codenvu == 1)
			{
			  //Coucou(nbrbuffer);
			  int buffersize2 = nbrevent2;
			  
			  //cout << "Je trie un autre buffer que le premier" << endl;
			  nbrbuffer++;

			  // Je faire le tri sur le "buffer taille reduite"
			  tribuffer_temps_coding(buffer_tm, tampon_tm, tampon_enrj, buffer_enrj, tampon_marker, buffer_marker,tampon_idx, buffer_idx, buffer_index, nbrcodingbuffer,buffersize2,nbrdetect);
			  
			  if(nbrbuffer < 3)
			    {
			      for(int k=0; k < 4; k++) cout << endl;
			      
			      cout << "*************************************" << endl;
			      cout << "J'affiche un tampon avec un coden "<< endl;
			    }
			  
			  // Now the buffer is sorted I will be able to store it in the TTree
			  for(int ll = 0; ll < buffersize2; ll++)
			    {
			      //if((buffer_tm[ll-1]- buffer_tm[ll]) > 1.e5) codingenable++;
			      index = buffer_idx[ll];
			      marker = buffer_marker[ll];
			      enrj = buffer_enrj[ll];
			      tm = buffer_tm[ll];
			      if(buffer_tm[ll] < 1.e-8) tm = 0;
			      			      
			      if(nbrbuffer < 3) cout << " Index : " << (int)index << "; Marker : " << (int)marker << "; Nrj : " << (int)enrj << "; Time : " << (double)tm << "; Coden : " << (int)codingenable << endl;

			      if(tm != 0 && (index < nbrdetect + 1)) oak -> Fill();
			      //oak -> Print();
			    }
			  
			  // J'ai eu un coding enable dans ce buffer... donc j'incremente le conteur
			  codingenable++;
			  
			  if(nbrbuffer < 2) cout << endl;

			  // Je remets à zero le nombre d'evenements pour remplir un nouveau buffer de taille normale
			  nbrevent2 = 0;
			  memory_nbrevent_coding = -666;
			  codenvu = 0;			  
			}// Fin du if dans le cas d'un coden dans le buffer
			  
		      // Dans le cas d'un buffer sans coden je fais un tri "standard"
		      if(nbrevent2 == buffersize && codenvu != 1)
			{
			  //Coucou(nbrbuffer);
			  // Tri des buffers
			  tribuffer_temps(buffer_tm, buffer_index, tampon_tm, buffersize, nbrdetect);
			  tri_buffers(tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, buffersize, nbrdetect);
			  //cout << "On en est au buffer " << nbrbuffer << endl;
			  nbrbuffer++;
			

			  
			   if(nbrbuffer < 2) for(int k=0; k < 4; k++) cout << endl;

			   if(nbrbuffer < 2) cout << "*************************************" << endl;
			   if(nbrbuffer < 2) cout << "J'affiche un tampon sans coden "<< endl;


			  // Now the buffer is sorted I will be able to store it in the TTree
			  for(int ll = 0; ll < buffersize; ll++)
			    {
			      //if(nbrevent > 3000000)Coucou(ll);
			      index = buffer_idx[ll];
			      marker = buffer_marker[ll];
			      enrj = buffer_enrj[ll];
			      tm = buffer_tm[ll];
			      if(buffer_tm[ll] < 1.e-8) tm = 0;
			      			    
			      if(nbrbuffer < 2) cout << " Index : " << (int)index << "; Marker : " << (int)marker << "; Nrj : " << (int)enrj << "; Time : " << (double)tm << "; Coden : " << (int)codingenable << endl;

  
			      if(tm != 0 && (index < nbrdetect + 1)) oak -> Fill();
			    }
			  
			  if(nbrbuffer < 2) cout << endl;

			  nbrevent2 = 0;
			  nbrcodingbuffer = 0;
			  memory_nbrevent_coding = -666;
			}// End of if(nbrevent2 == buffersize)
		    }// End of if(idx)
		  
		  else
		    {
		      cout << "Foolish event detected: detector not existing" << endl;
		      cout << "group : " << (int)gr << endl;
		      cout << "slot : " << (int)sl << endl;
		      cout << "crate : " << (int)cr << endl;
		      cout << "service : " << (int)ser << endl;
		      cout << "marker : " << (int)mark << endl;
		      cout << "channel : " << (int)ch << endl;
		      cout << "It is the event number : " << nbrevent2 << endl;
		      cout << "Of the " << nbrbuffer << "e buffer" << endl;
		      nbreventfous++;
		    }
		}// End of if(ser == 0)

	      
	  
	      // Read the next word to able to check if it is a physical hit or not
	      if (narval.good()) narval.read(word1,2);
	      if (narval.good()) narval.read(word2,2);
	      
	      // In case of incomplete event, eof is reached unexpectedly
	      if (narval.eof())
		{
		  cout << "Irregular end of file, incomplete event" << endl;
		  word1[0]=0xff;
		  word1[1]=0xff;
		  word2[0]=0xff;
		  word2[1]=0xff;                              // set marker for eof
		}// End of if (narval.eof())
	    }// End of while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	  
	  // If file contains loose zeroes in the middle, mark it as an internal header
	  if((word_swap(word1[0],word1[1])==0x0) && (word_swap(word2[0],word2[1])==0x0))
	    {
	      cout << "Irregular file structure,loose zeroes found at "<< narval.tellg()-(streampos)4<< endl;
	      // narval.seekg(0,ios::end);                  // move to eof without provoking eof yet
	      // narval.read(word2,2);                      // actually provoke eof
	      word1[0]=0xff;
	      word1[1]=0xff;
	      word2[0]=0xff;
	      word2[1]=0xff;                                 // set marker for eof
	    }// End of if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	  
	  // This is an internal header or regular eof or eof incomplete event or loose zeroes
	  else if((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0xffff))
	    {
	      if (narval.good()) narval.read(word1, 2);
	      if (narval.good()) narval.read(word2, 2);

	      // execute the loop as long as we don't encounter a detector event
	      //cout<<" apres "<<n<<" evenements on cherche de nouveaux evenements physiques"<<endl ;
	      while(!((word_swap(word1[0],word1[1])==0xffff) && ( (word_swap(word2[0],word2[1])==0x8) || (word_swap(word2[0],word2[1])==0xd))))
		{
		  //	  cout<<n<<endl ; //ASP
		  narval.seekg(-2,ios::cur);
		  //ASP
		  //  cout<<narval.fail()<< " "<<narval.good()<< " "<<narval.bad()<< " "<<narval.eof()<<" "<<n<< endl ;
		  if(narval.good()) narval.read(word1,2);
		  if(narval.good()) narval.read(word2,2);
		  if(narval.eof()) 
		    { 
		      cout << "End of file detected" << endl;
		      // not very elegant but we need to escape from the while loop
		      word1[0]=0xff;
		      word1[1]=0xff;
	 	      word2[0]=0x0;
		      word2[1]=0x8;
		    }// End of if(narval.eof())
		}// End of while(!((word_swap(word1[0],word1[1])==0xffff) && ( (word_swap(word2[0],word2[1])==0x8) || (word_swap(word2[0],word2[1])==0xd)))) 
	    }// End of else if((word_swap(word1[0],word1[1])==0xffff) && (word_swap(word2[0],word2[1])==0xffff))
	  
	  else
	    {
	      cout << "Fatal error, irregular event structure "<< endl;
	      return n;
	    }     	
	}// End of while for eof
      
      // Ending of the conversion procedure
      if(narval.eof())
	{
	  // Organization of the last buffer that might have a size < buffersize
	  buffersize = nbrevent2;
	 	  
	  tribuffer_temps(buffer_tm, buffer_index, tampon_tm, buffersize, nbrdetect);
	  tri_buffers(tampon_enrj,buffer_enrj,tampon_marker,buffer_marker,tampon_idx, buffer_idx, buffer_index, buffersize, nbrdetect);
	  
	  nbrbuffer++;

	  // Now the buffer is sorted I will be able to store it in the TTree
	  for(int ll = 0; ll < buffersize; ll++)
	    {
	      index = buffer_idx[ll];
	      marker = buffer_marker[ll];
	      enrj = buffer_enrj[ll];
	      tm = buffer_tm[ll];
	      if(buffer_tm[ll] < 1.e-8) tm = 0;
		
	      if(tm !=0 ) oak -> Fill();
	    }
	  // Storage of the "tampon"
	  for(int ll = 0; ll < nbrdetect; ll++)
	    {
	      index = tampon_idx[ll];
	      marker = tampon_marker[ll];
	      enrj = tampon_enrj[ll];
	      tm = tampon_tm[ll];
	      oak -> Fill(); 
	    }
	  


	  // Writing the tree in the root file to be sure the root file is well closed
	  cout << "Saving the tree in the file" << endl;
	  oak -> Write();
	  
	  
	  // Unallocation of the memory used for the buffer
	  delete [] bf16;
	  delete [] word1;
	  delete [] word2;
	  //delete [] branchdetname;     
	    
	  // Closing of all the file used during the conversion process
	  narval.close();
	  outputfile.Close();
	  
	  // Effacement des differents tableaux
	  /*
	  	  delete[] buffer_marker;
	  delete[] buffer_enrj;
	  delete[] buffer_tm;
	  delete[] buffer_idx;
	  delete[] buffer_index;
	  delete[] tampon_tm;
	  delete[] tampon_enrj;
	  delete[] tampon_marker; 
	  delete[] tampon_idx;
	  */

	  // Writting of a summary of the conversion process
	  cout << "This makes "<< codingenable << " coden read in the file" << endl;
	  cout << "This makes "<< n << " hits read in the file"<<endl;
	  cout << "This makes "<< n-nbrevent << " service hits read in the file"<<endl;
	  cout << "This makes "<< nbreventfous << " completely foolish events read in the file"<<endl;
	  cout << "This makes "<< nbrdesordre << " hits that were not in order (in time)"<<endl;
	  cout << "This makes "<< nbrevent << " physical hits read in the file"<<endl;
	  cout << endl;

	}// End of if(narval.eof())
      
      else
	{
	  cout << "Error exiting program" << endl;   
	}

      break;



      //--------------------------------------------//
      //      No coding enable and no swap          //
      //--------------------------------------------//
    default:
      cout << "Filesize is " << inputfilesize << " bytes for "<< inputfilename << endl;
      
      // Permission for a root file to exceed 2GB
      if(inputfilesize > 2000000000)
	{
	  oak->SetMaxTreeSize(inputfilesize); 
	}
      cout << "The root file size will be of ~ " << (Double_t)inputfilesize/2000000000 << " Gb"<< endl;
      
      narval.seekg(0,ios::beg);
      if (inputfilesize < 0 ) 
	{
	  cout <<"$$$"<<endl <<" $$$ Problem with file "<<inputfilename<<endl ; 
	  cout << "$$$ skiping it"<<endl<<"$$$ "<<endl ; 
	  break;
	}

      // Skip header till we find an event 0xffff 0x8
      word1 = new char [2];
      word2 = new char [2];
      bf16 = new char [2];
      narval.read(word1,2);
      narval.read(word2,2);
  
      while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
	{
	  // Move one word backwards but read two
	  narval.seekg(-2,ios::cur);
	  if (narval.good()) narval.read(word1,2);
	  if (narval.good()) narval.read(word2,2);
	  
	  // Test to check the end of file (corrupted file)
	  if (narval.eof())
	    {
	      cout << "Fatal error, file does not start with detector event" << endl;
	      return n;
	    }// end of if (narval.eof())
	}// end of while ((word(word1[0],word1[1]) != 0xffff)||(word(word2[0],word2[1]) !=0x8))
      
      // Writting of the declaration for the good start of conversion
      cout << "Start conversion at position " << narval.tellg()-(streampos)4<< endl;
      cout << "-------------------------------------------------"<< endl; 
      
      // Scan of all the data file
      while (narval.good())
	{      
	  // Check for a physical event: detector event 0xffff 0x0008
	  while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	    {
	      // Increment of any event counter
	      n++;
	      
	      // Extraction of hit description
	      if (narval.good()) narval.read(bf16,2);        
	      
	      // Read group, slot, crate, ch, marker and service
	      bf = word(bf16[0],bf16[1]);
	      gr = (UChar_t)(bf >> 12);
	      sl = (UChar_t)((bf >> 8) & 15);
	      cr = (UChar_t)((bf >> 6) & 3);
	      ser = (UChar_t)((bf >> 4 ) & 2);
	      mark = (UChar_t)((bf >> 4 ) & 1);
	      ch = (UChar_t)(bf & 15);
	  
	      // Creation of the detector number to find it in the index
	      sprintf(nomdetector,"idx_detect_%d_%d_%d_%d",gr,cr,sl,ch);
	      
	      // Read energy
	      if (narval.good()) narval.read(bf16,2);
	      en = (UInt_t)word(bf16[0],bf16[1]);
	      
	      // Read time
	      if (narval.good()) narval.read(bf16,2);
	      tm1 = (word(bf16[0],bf16[1]));
	      if (narval.good()) narval.read(bf16,2);
	      tm2 = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      chks = word(bf16[0],bf16[1]);
	      if (narval.good()) narval.read(bf16,2);
	      tm3 = (word(bf16[0],bf16[1]));
	      
	      // Time reconstruction to be able to use it for checking
	      tm = TEMPS_47bits(tm3,tm1,tm2);

	      // Filling of the tree used for storing PHYSICAL hits
	      if(ser == 0) // checking it is no a service hit 
		{
		  // Reading of the input tree to get the information concerning the detector number 
		  TBranch *idx = input->GetBranch(nomdetector); 
		  
		  // Check for the quality of the information read
		  if(idx) 
		    { 
		      // Reading of the index to get the detector number corresponding to a TClonesArray
		      idx->SetAddress(&det_num); //##??
		      idx->GetEntry(0);
		      index = det_num.nomdetect;
		  
		      // Verification of the order of the event by detector to check hit are stored
		      // in chronological order
		      if((tm*400 - memoiretemps) >= 0)
			{
			  // Filling of the structure to prepare the storing of the event
			  marker = mark;
			  enrj = en;
			  tm = tm*400; // to be in unit of ps
			 			  
			  // Filling of the tree
			  oak -> Fill();
			  
			  // Augmentation of hit counter
			  nbrevent++;
			  			  
			  // Memorization of the time for this detector to test the order of events 
			  // coming from the same detector
			  memoiretemps = (Double_t)tm;
			  
			}// End of if(tm >= memoiretemps[index])
		      else
			{
			  cout << "An event is not in order in detector " << (int)index << endl;
			  cout << "Present Time = " << tm << endl;
			  cout << "Referred to " << memoiretemps << endl;
			  memoiretemps = 0;
			  nbrdesordre++;
			}
		      
		    }// End of if(idx)
	      
		  else
		    {
		      cout << "Foolish event detected: detector not existing" << endl;
		      cout << "group : " << (int)gr << endl;
		      cout << "slot : " << (int)sl << endl;
		      cout << "crate : " << (int)cr << endl;
		      cout << "service : " << (int)ser << endl;
		      cout << "marker : " << (int)mark << endl;
		      cout << "channel : " << (int)ch << endl;
		      nbreventfous++;
		    }
		}// End of if(ser == 0)
	  
	      // Read the next word to able to check if it is a physical hit or not
	      if (narval.good()) narval.read(word1,2);
	      if (narval.good()) narval.read(word2,2);
	      
	      // In case of incomplete event, eof is reached unexpectedly
	      if (narval.eof())
		{
		  cout << "Irregular end of file, incomplete event" << endl;
		  word1[0]=0xff;
		  word1[1]=0xff;
		  word2[0]=0xff;
		  word2[1]=0xff;                              // set marker for eof
		}// End of if (narval.eof())
	    }// End of while ((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0x8))
	  
	  // If file contains loose zeroes in the middle, mark it as an internal header
	  if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	    {
	      cout << "Irregular file structure,loose zeroes found at "<< narval.tellg()-(streampos)4<< endl;
	      // narval.seekg(0,ios::end);                  // move to eof without provoking eof yet
	      // narval.read(word2,2);                      // actually provoke eof
	      word1[0]=0xff;
	      word1[1]=0xff;
	      word2[0]=0xff;
	      word2[1]=0xff;                                 // set marker for eof
	    }// End of if((word(word1[0],word1[1])==0x0) && (word(word2[0],word2[1])==0x0))
	  
	  // This is an internal header or regular eof or eof incomplete event or loose zeroes
	  else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	    {
	      if (narval.good()) narval.read(word1, 2);
	      if (narval.good()) narval.read(word2, 2);
	      
	      // execute the loop as long as we don't encounter a detector event
	      while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd))))
		{
		  narval.seekg(-2,ios::cur);
		  if(narval.good()) narval.read(word1,2);
		  if(narval.good()) narval.read(word2,2);
		  if(narval.eof()) 
		    { 
		      cout << "End of file detected" << endl;
		      // not very elegant but we need to escape from the while loop
		      word1[0]=0xff;
		      word1[1]=0xff;
		      word2[0]=0x8;
		      word2[1]=0x0;
		    }// End of if(narval.eof())
		}// End of while(!((word(word1[0],word1[1])==0xffff) && ( (word(word2[0],word2[1])==0x8) || (word(word2[0],word2[1])==0xd)))) 
	    }// End of else if((word(word1[0],word1[1])==0xffff) && (word(word2[0],word2[1])==0xffff))
	  
	  else
	    {
	      cout << "Fatal error, irregular event structure "<< endl;
	      return n;
	    }     	
	}// End of while for eof
      
      // Ending of the conversion procedure
      if(narval.eof())
	{
	  // Writing the tree in the root file to be sure the root file is well closed
	  cout << "Saving the tree in the file" << endl;
	  oak -> Write();
	  
	  // Unallocation of the memory used for the buffer
	  delete [] bf16;
	  delete [] word1;
	  delete [] word2;
	  //delete [] branchdetname;     
	    
	  // Closing of all the file used during the conversion process
	  narval.close();
	  outputfile.Close();
	  
	  // Writting of a summary of the conversion process 
	  cout << "This makes "<< n << " hits read in the file"<<endl;
	  cout << "This makes "<< n-nbrevent << " service hits read in the file"<<endl;
	  cout << "This makes "<< nbreventfous << " completely foolish events read in the file"<<endl;
	  cout << "This makes "<< nbrdesordre << " hits that were not in order (in time)"<<endl;
	  cout << "This makes "<< nbrevent << " physical hits read in the file"<<endl;
	  cout << endl;

	}// End of if(narval.eof())
      
      else
	{
	  cout << "Error exiting program" << endl;   
	}
      
      break;

    }// End switch(option)

  // Return an int because it is a function
  return 1;

}// End of convertfile function
