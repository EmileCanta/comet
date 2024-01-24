//-------------------------------------------------------------------------//
//                                                                         //
//                           Fonctions_convert.cxx                         //
//                               Version 1.8                               //
//                        Matthieu Lebois November 2007                    //
//                            with G. Georgiev                             //
//                                                                         //
//     This file contain the function to order all the event in time.      //
//                                                                         //
//-------------------------------------------------------------------------//

#include "Fonctions_convert.h"

using namespace std;

// Definition of a function to check the weight of a bit
unsigned short word(UChar_t swap, char lo, char up)
{
  unsigned long result;
  int option = (int) swap;
  switch(option)
    {
      // Pas de swap
    case 0:
      // conversion of two 8 bit char to single 16 bit integer
      // to be sure, set higher bits of 'up' and 'lo' to zero
      result = (up&255)*256+(lo&255);
      break;

      // Y a du swap (dans le bar tabac d'en bas...)
    case 1:
      // conversion of two 8 bit char to single 16 bit integer
      // to be sure, set higher bits of 'up' and 'lo' to zero
      result = (lo&255)*256+(up&255);
      break;

      // Cas par defaut: y a du swap. Tout les fichier presents sur le p5 sont swappee
    default:
      // conversion of two 8 bit char to single 16 bit integer
      // to be sure, set higher bits of 'up' and 'lo' to zero
      result = (lo&255)*256+(up&255);
      break;
    }
  return result;
}


// // Definition of a function to check the weight of a bit in a case of a swap in the Data
// unsigned short word_swap(char lo, char up)
// {
//   unsigned long result;
  
//   // conversion of two 8 bit char to single 16 bit integer
//   // to be sure, set higher bits of 'up' and 'lo' to zero
//   result = (lo&255)*256+(up&255);
//   return result;
// }



// Fonction pour le tout premier tampon
void first_tri_buffers(Double_t tampon_tm[], Double_t buffer_tm[], UInt_t tampon_enrj[],UInt_t  buffer_enrj[], UChar_t tampon_marker[], UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], int buffersize, int nbrdetect)
{
  // Le premier buffer doit faire "nbrdetector" evenements
  // Car il arrive juste apres un coden (ou le debut du run)
  
  // L'objectif est de trier les temps dans ce qui va constituer le premier tampon
  
  // Variable utilisees dans cette fonction
  Int_t travailsize = buffersize; 

  // Creation du buffer pour l'index afin de remttre les tableaux dans l'ordre
  Int_t *buffer_index = new Int_t[travailsize];

  // Creation of a pointeur tostore the data
  //coming from the previous buffer and from the actual buffer
  Double_t *tm_travail = new Double_t[travailsize];
  UInt_t *enrj_travail = new UInt_t[travailsize];
  UChar_t *marker_travail = new UChar_t[travailsize];
  UChar_t *idx_travail = new UChar_t[travailsize];
  Double_t *tm_travail_ordered = new Double_t[travailsize];
  UInt_t *enrj_travail_ordered = new UInt_t[travailsize];
  UChar_t *marker_travail_ordered = new UChar_t[travailsize];
  UChar_t *idx_travail_ordered = new UChar_t[travailsize];


  // Filling of the "travail pointeur"
  // Copying the actual buffer
  for(int i = 0; i < travailsize; i++)
    {
      tm_travail[i]     = buffer_tm[i];
      enrj_travail[i]   = buffer_enrj[i];
      marker_travail[i] = buffer_marker[i];
      idx_travail[i]    = buffer_idx[i];
    }

  // Tri du buffer de temps
  TMath::Sort(travailsize,tm_travail,buffer_index,0);

  // Ordering of the "travail pointeurs"
  for(int i = 0; i < travailsize; i++)
    {
      tm_travail_ordered[i]     =  tm_travail[(int)buffer_index[i]];
      enrj_travail_ordered[i]   =  enrj_travail[(int)buffer_index[i]];
      marker_travail_ordered[i] =  marker_travail[(int)buffer_index[i]];
      idx_travail_ordered[i]    =  idx_travail[(int)buffer_index[i]];
    }

  // Saving in the buffer the reorganized "travail pointeur"

  // Copy of the end of the sorted buffer in the memory to compare with the previous buffer
  for(int i = 0; i < travailsize; i++)
    {
      tampon_tm[i]     = tm_travail_ordered[i];
      tampon_enrj[i]   = enrj_travail_ordered[i];
      tampon_marker[i] = marker_travail_ordered[i];
      tampon_idx[i]    = idx_travail_ordered[i];
    }

  // J'affiche le tampon pour verifier qu'il n'y a pas d'evenement dans le desordre au dela du coden
//   for(int i = 0; i < nbrdetect; i++)
//     {
//       cout << "i = " << i  << "; tampon_tm = " << (double) tampon_tm[i] << "; tampon_enrj = " << (double)tampon_enrj[i] << "; tampon_marker = " << (int)tampon_marker[i] << "; tampon_idx = " << (int)tampon_idx[i] << endl;
//     }
  
  // Liberation de la memoire

  delete buffer_index;
  delete tm_travail;
  delete enrj_travail;
  delete marker_travail;
  delete idx_travail;
  delete tm_travail_ordered;
  delete enrj_travail_ordered;
  delete marker_travail_ordered;
  delete idx_travail_ordered;
}


// Fonction pour trier les autres buffers
void tri_buffers(Double_t tampon_tm[], Double_t buffer_tm[], UInt_t tampon_enrj[],UInt_t  buffer_enrj[],UChar_t tampon_marker[],UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], int buffersize, int nbrdetect)
{
  // Definition de la taille des tableaux qui seront manipules dans cette fonction
  Int_t travailsize = buffersize + nbrdetect; 

  // Creation du buffer pour l'index afin de remttre les tableaux dans l'ordre
  Int_t *buffer_index = new Int_t[travailsize];

  // Creation of a pointeur tostore the data
  //coming from the previous buffer and from the actual buffer
  Double_t *tm_travail = new Double_t[travailsize];
  UInt_t *enrj_travail = new UInt_t[travailsize];
  UChar_t *marker_travail = new UChar_t[travailsize];
  UChar_t *idx_travail = new UChar_t[travailsize];
  Double_t *tm_travail_ordered = new Double_t[travailsize];
  UInt_t *enrj_travail_ordered = new UInt_t[travailsize];
  UChar_t *marker_travail_ordered = new UChar_t[travailsize];
  UChar_t *idx_travail_ordered = new UChar_t[travailsize];


  // Filling of the "travail pointeur"
  // Copying the previous buffer
  for(int i = 0; i < nbrdetect; i++)
    {
      tm_travail[i] = tampon_tm[i];
      enrj_travail[i] = tampon_enrj[i];
      marker_travail[i] = tampon_marker[i];
      idx_travail[i] = tampon_idx[i];
    }
  // Copying the actual buffer
  for(int i = 0; i < buffersize; i++)
    {
      tm_travail[(int)i+nbrdetect] = buffer_tm[(int)i];
      enrj_travail[(int)i+nbrdetect] = buffer_enrj[(int)i];
      marker_travail[(int)i+nbrdetect] = buffer_marker[(int)i];
      idx_travail[(int)i+nbrdetect] = buffer_idx[(int)i];
    }
  
  // Tri du buffer
  TMath::Sort(travailsize,tm_travail,buffer_index,0);
 
  // Ordering of the "travail pointeurs"
  for(int i = 0; i < travailsize; i++)
    {
      tm_travail_ordered[i] = tm_travail[(int)buffer_index[i]];
      enrj_travail_ordered[i] = enrj_travail[(int)buffer_index[i]];
      marker_travail_ordered[i] = marker_travail[(int)buffer_index[i]];
      idx_travail_ordered[i] = idx_travail[(int)buffer_index[i]];
    }

  // Saving in the buffer the reorganized "travail pointeur"
  // Copy of the sorted buffer in the time buffer
  for(int i = 0; i < buffersize; i++)
    {
      buffer_tm[i] = tm_travail_ordered[i];
      buffer_enrj[i] = enrj_travail_ordered[i];
      buffer_marker[i] = marker_travail_ordered[i];
      buffer_idx[i] = idx_travail_ordered[i];
    }

  // Copy of the end of the sorted buffer in the memory to compare with the previous buffer
  for(int i = buffersize; i < travailsize; i++)
    {
      tampon_tm[(int)i-buffersize] = tm_travail_ordered[i];
      tampon_enrj[(int)i-buffersize] = enrj_travail_ordered[i];
      tampon_marker[(int)i-buffersize] = marker_travail_ordered[i];
      tampon_idx[(int)i-buffersize] = idx_travail_ordered[i];
    }

  // Desallocation memoire
  
  delete  tm_travail;
  delete  enrj_travail;
  delete  marker_travail;
  delete  idx_travail;
  delete  tm_travail_ordered;
  delete  enrj_travail_ordered;
  delete  marker_travail_ordered;
  delete  idx_travail_ordered;
  delete  buffer_index;
  
}

// Definition of the function that wil dispatch event in two parts: before and after coden
void tri_preandpost_coden(Double_t  buffer_coden_tm[], UInt_t buffer_coden_enrj[], UChar_t buffer_coden_idx[], UChar_t buffer_coden_marker[],Int_t buffersize, int *position1erevenement_postcoden, int *positionderniererevenement_postcoden)
{
  // Declaration des variables crees pour cette fonction
  Double_t tm_min, tm_max, tm_mean;
  int counter_precoden, counter_postcoden;

  // Declaration de sous tableaux pour un traiement pre-coden et post coden
  Double_t *buffer_precoden_tm = new Double_t[buffersize];
  UInt_t *buffer_precoden_enrj = new UInt_t[buffersize];
  UChar_t *buffer_precoden_marker = new UChar_t[buffersize];
  UChar_t *buffer_precoden_idx = new UChar_t[buffersize];

  Double_t *buffer_postcoden_tm = new Double_t[buffersize];
  UInt_t *buffer_postcoden_enrj = new UInt_t[buffersize];
  UChar_t *buffer_postcoden_marker = new UChar_t[buffersize];
  UChar_t *buffer_postcoden_idx = new UChar_t[buffersize];

  Double_t *buffer_precoden_tm_ordered = new Double_t[buffersize];
  UInt_t *buffer_precoden_enrj_ordered = new UInt_t[buffersize];
  UChar_t *buffer_precoden_marker_ordered = new UChar_t[buffersize];
  UChar_t *buffer_precoden_idx_ordered = new UChar_t[buffersize];

  Double_t *buffer_postcoden_tm_ordered = new Double_t[buffersize];
  UInt_t *buffer_postcoden_enrj_ordered = new UInt_t[buffersize];
  UChar_t *buffer_postcoden_marker_ordered = new UChar_t[buffersize];
  UChar_t *buffer_postcoden_idx_ordered = new UChar_t[buffersize];

  // Creation d'un tableaux d'index qui sera utilise pour sorter les deux types de sous tableaux
  int *tab_index = new int[buffersize];

//   for(int ii = 0; ii < buffersize; ii++)
//     {
//       cout << ii << " buffer_marker_coden = " << (int)buffer_coden_marker[(int)(ii)] << "; buffer_idx_coden = " << (int)buffer_coden_idx[(int)(ii)] << "; buffer_enrj_coden = " << (double)buffer_coden_enrj[(int)(ii)] << "; buffer_coden_tm = " << (double)buffer_coden_tm[(int)(ii)] << endl;
//     }

  // Je recherche l'element max dans le tableau de temps
  tm_max = TMath::MaxElement(buffersize,buffer_coden_tm);
  tm_min = TMath::MinElement(buffersize,buffer_coden_tm);

  //cout << "tm_max = " << tm_max << endl;
  //cout << "tm_min = " << tm_min << endl;

  tm_mean = (tm_max-tm_min)/2.;

//   cout << "tm_mean = " << tm_mean << endl;

  // Je remplis les deux "sous" tableaux par rapport a cette moyenne
  counter_precoden = 0;
  counter_postcoden = 0;

  for(int i = 0; i < buffersize; i++)
    {
      if(buffer_coden_tm[i] > tm_mean)
	{
	  buffer_precoden_tm[counter_precoden] = buffer_coden_tm[i];
	  buffer_precoden_enrj[counter_precoden] = buffer_coden_enrj[i];
	  buffer_precoden_marker[counter_precoden] = buffer_coden_marker[i];
	  buffer_precoden_idx[counter_precoden] = buffer_coden_idx[i];
	  counter_precoden++;

	  
	}
      if(buffer_coden_tm[i] < tm_mean)
	{
	  buffer_postcoden_tm[counter_postcoden] = buffer_coden_tm[i];
	  buffer_postcoden_enrj[counter_postcoden] = buffer_coden_enrj[i];
	  buffer_postcoden_marker[counter_postcoden] = buffer_coden_marker[i];
	  buffer_postcoden_idx[counter_postcoden] = buffer_coden_idx[i];
	  counter_postcoden++;
	}
    }

  //cout << "counter_precoden = " << counter_precoden << "; counter_postcoden = " << counter_postcoden << endl;

  // Je peux desormais trier independamment mes deux "sous" buffers
  TMath::Sort(counter_precoden, buffer_precoden_tm, tab_index,0);

  // Le buffer precoden est trie, je tri les autre buffer a partir de l'index
  for(int j = 0; j < counter_precoden; j++)
    {
      buffer_precoden_tm_ordered[j] = buffer_precoden_tm[(int)tab_index[j]];
      buffer_precoden_enrj_ordered[j] = buffer_precoden_enrj[(int)tab_index[j]];
      buffer_precoden_marker_ordered[j] = buffer_precoden_marker[(int)tab_index[j]];
      buffer_precoden_idx_ordered[j] = buffer_precoden_idx[(int)tab_index[j]];  
    }
  // Tout mes buffers precoden sont tries
 
  // Je fais la meme chose sur les post coden
  // Je peux desormais trier independamment mes deux "sous" buffers
  TMath::Sort(counter_postcoden, buffer_postcoden_tm, tab_index,0);

  // Le buffer precoden est trie, je tri les autre buffer a partir de l'index
  for(int j = 0; j < counter_postcoden; j++)
    {
      buffer_postcoden_tm_ordered[j] = buffer_postcoden_tm[(int)tab_index[j]];
      buffer_postcoden_enrj_ordered[j] = buffer_postcoden_enrj[(int)tab_index[j]];
      buffer_postcoden_marker_ordered[j] = buffer_postcoden_marker[(int)tab_index[j]];
      buffer_postcoden_idx_ordered[j] = buffer_postcoden_idx[(int)tab_index[j]];  
    }
  
  (*position1erevenement_postcoden) = counter_precoden;
  (*positionderniererevenement_postcoden) = counter_postcoden;

  // Mes buffer post et pre coden sont enfin tries... (enfin j'espere) Je les reenregistre
  for(int j = 0; j < counter_precoden; j++)
    {
      buffer_coden_tm[j] = buffer_precoden_tm_ordered[j];
      buffer_coden_enrj[j] = buffer_precoden_enrj_ordered[j];
      buffer_coden_marker[j] = buffer_precoden_marker_ordered[j];
      buffer_coden_idx[j] = buffer_precoden_idx_ordered[j];
    }

  for(int j = 0; j < counter_postcoden; j++)
    {
      buffer_coden_tm[(int)(j+(counter_precoden))] = buffer_postcoden_tm_ordered[j];
      buffer_coden_enrj[(int)(j+(counter_precoden))] = buffer_postcoden_enrj_ordered[j];
      buffer_coden_marker[(int)(j+(counter_precoden))] = buffer_postcoden_marker_ordered[j];
      buffer_coden_idx[(int)(j+(counter_precoden))] = buffer_postcoden_idx_ordered[j];
    }
 
  // La mes tableaux sont revenus dans l'ordre!!!!

  // Nettoyage de la memoire
  
  delete buffer_precoden_tm;
  delete buffer_precoden_enrj;
  delete buffer_precoden_marker;
  delete buffer_precoden_idx;
  delete buffer_postcoden_tm;
  delete buffer_postcoden_enrj;
  delete buffer_postcoden_marker;
  delete buffer_postcoden_idx;
  delete buffer_precoden_tm_ordered;
  delete buffer_precoden_enrj_ordered;
  delete buffer_precoden_marker_ordered;
  delete buffer_precoden_idx_ordered;
  delete buffer_postcoden_tm_ordered;
  delete buffer_postcoden_enrj_ordered;
  delete buffer_postcoden_marker_ordered;
  delete buffer_postcoden_idx_ordered;
  delete tab_index;
  
}



// Function to sort the small buffer of event post coden (again without filling of tampon)
void tri_buffer_postcoden(Double_t buffer_tm_coden[], UInt_t buffer_enrj_coden[], UChar_t buffer_marker_coden[], UChar_t buffer_idx_coden[], Int_t buffersize)
{
// Definition de la taille des tableaux qui seront manipules dans cette fonction
  Int_t travailsize = buffersize; 

  // Creation du buffer pour l'index afin de remttre les tableaux dans l'ordre
  Int_t *buffer_index = new Int_t[travailsize];

  // Creation of a pointeur tostore the data
  //coming from the previous buffer and from the actual buffer
  Double_t *tm_travail = new Double_t[travailsize];
  UInt_t *enrj_travail = new UInt_t[travailsize];
  UChar_t *marker_travail = new UChar_t[travailsize];
  UChar_t *idx_travail = new UChar_t[travailsize];
  Double_t *tm_travail_ordered = new Double_t[travailsize];
  UInt_t *enrj_travail_ordered = new UInt_t[travailsize];
  UChar_t *marker_travail_ordered = new UChar_t[travailsize];
  UChar_t *idx_travail_ordered = new UChar_t[travailsize];


  // Copying the actual buffer
  for(int i = 0; i < travailsize; i++)
    {
      tm_travail[i] = buffer_tm_coden[(int)i];
      enrj_travail[i] = buffer_enrj_coden[(int)i];
      marker_travail[i] = buffer_marker_coden[(int)i];
      idx_travail[i] = buffer_idx_coden[(int)i];
    }
  
  // Tri du buffer
  TMath::Sort(travailsize,tm_travail,buffer_index,0);
  
  // Ordering of the "travail pointeurs"
  for(int i = 0; i < travailsize; i++)
    {
      tm_travail_ordered[i] = tm_travail[(int)buffer_index[i]];
      enrj_travail_ordered[i] = enrj_travail[(int)buffer_index[i]];
      marker_travail_ordered[i] = marker_travail[(int)buffer_index[i]];
      idx_travail_ordered[i] = idx_travail[(int)buffer_index[i]];
    }

  // Saving in the buffer the reorganized "travail pointeur"
  // Copy of the sorted buffer in the time buffer
  for(int i = 0; i < travailsize; i++)
    {
      buffer_tm_coden[i] = tm_travail_ordered[i];
      buffer_enrj_coden[i] = enrj_travail_ordered[i];
      buffer_marker_coden[i] = marker_travail_ordered[i];
      buffer_idx_coden[i] = idx_travail_ordered[i];
    }

  // Desallocation memoire
  
  delete  tm_travail;
  delete  enrj_travail;
  delete  marker_travail;
  delete  idx_travail;
  delete  tm_travail_ordered;
  delete  enrj_travail_ordered;
  delete  marker_travail_ordered;
  delete  idx_travail_ordered;
  delete  buffer_index;
  
}

