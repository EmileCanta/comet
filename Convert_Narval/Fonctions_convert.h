//-------------------------------------------------------------------------//
//                                                                         //
//                            Fonctions_convert.h                          //
//                               Version 1.8                               //
//                        Matthieu Lebois November 2007                    //
//                            with G. Georgiev                             //
//                                                                         //
//     This file contain the function of conversion to root file.          //
//                                                                         //
//-------------------------------------------------------------------------//



// Call of all the libraries used in this program
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include "TFile.h"
#include "TMath.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TStopwatch.h"


#ifndef FONCTIONS_H
#define FONCTIONS_H

// Fonction de lecture d'un mot de 16 bits
unsigned short word(UChar_t swap, char lo, char up);

// Fonction pour le tout premier tampon
void first_tri_buffers(Double_t tampon_tm[], Double_t buffer_tm[], UInt_t tampon_enrj[],UInt_t  buffer_enrj[], UChar_t tampon_marker[], UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], int buffersize, int nbrdetect);

// Fonction pour trier les autres buffers
void tri_buffers(Double_t tampon_tm[], Double_t buffer_tm[], UInt_t tampon_enrj[],UInt_t  buffer_enrj[],UChar_t tampon_marker[],UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], int buffersize, int nbrdetect);

// Fonction qui va trier un buffer de ~2*nbrdetect evenements autour du coding
//void tri_buffer_coding(Double_t buffer_coden_tm[], UInt_t buffer_coden_enrj[], UChar_t buffer_coden_marker[],UChar_t buffer_coden_idx[], int sizeofbuffer);

// Definition of a function that will sort a buffer of event with a coden inside
//void tribuffer_temps_coding(Double_t buffer_tm[], Double_t tampon_tm[],UInt_t tampon_enrj[],UInt_t  buffer_enrj[],UChar_t tampon_marker[],UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], Int_t nbrcodingbuffer, int buffersize, int nbrdetect);

// Definition of the function that wil dispatch event in two parts: before and after coden
void tri_preandpost_coden(Double_t  buffer_coden_tm[],UInt_t buffer_coden_enrj[],UChar_t buffer_coden_idx[], UChar_t buffer_coden_marker[],Int_t buffersize, int *position1erevenement_postcoden, int *positionderniererevenement_postcoden);

// Function to sort a buffer without coden but without filling of "tampon"
void tri_buffer_precoden(Double_t buffer_tm[], UInt_t buffer_enrj[], UChar_t buffer_marker[], UChar_t buffer_idx[],Int_t buffersize);

// Function to sort the small buffer of event post coden (again without filling of tampon)
void tri_buffer_postcoden(Double_t buffer_tm_coden[], UInt_t buffer_enrj_coden[], UChar_t buffer_marker_coden[], UChar_t buffer_idx_coden[], Int_t buffersize);

#endif
