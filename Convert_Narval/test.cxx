// Definition of a function that will sort a buffer of event
void tribuffer_temps(Double_t buffer_tm[], Int_t buffer_index[], Double_t tampon_tm[],int buffersize, int nbrdetect)
{
  // Declaration des variables utilisees dans cette fonction
  Int_t travailsize = buffersize + nbrdetect; 
 
  
  // Declaration of a pointeur of the size of the buffer and a number of detectors
  Double_t *tm_travail = new Double_t[travailsize];

  // First events of actual buffer are not in order
  // with last events of previous buffer

  // I have to order all of them between the previous buffer 
  // and the new one
  // Copying the previous buffer
  for(int i = 0; i < nbrdetect; i++)
    {
      tm_travail[i] = tampon_tm[i];
    }
  // Copying the actual buffer
  for(int i = nbrdetect; i < travailsize; i++)
    {
      tm_travail[i] = buffer_tm[(int)i-nbrdetect];
    }

  // Tri du buffer
  TMath::Sort(travailsize,tm_travail,buffer_index,0);

  // Copy of the sorted buffer in the time buffer
  for(int i = 0; i < buffersize; i++)
    {
      buffer_tm[i] = tm_travail[(int)buffer_index[i]];
    }

  // Copy of the end of the sorted buffer in the memory to compare with the previous buffer
  for(int i = buffersize; i < travailsize; i++)
    {
      tampon_tm[i-buffersize] = tm_travail[(int)buffer_index[i]];
    }

  // Now the buffer is supposed to be in order...

  // Desallocation memoire
  delete  tm_travail;

}



// Fonction pour trier les autres buffers
void tri_buffers(Double_t buffer_tm[], Double_t tampon_tm[], UInt_t tampon_enrj[],UInt_t  buffer_enrj[],UChar_t tampon_marker[],UChar_t buffer_marker[],UChar_t tampon_idx[],UChar_t buffer_idx[], int buffersize, int nbrdetect)
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
  Double_t *tm_ordered = new Double_t[travailsize];
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
  for(int i = nbrdetect; i < travailsize; i++)
    {
      tm_travail[i] = buffer_tm[(int)i-nbrdetect];
      enrj_travail[i] = buffer_enrj[(int)i-nbrdetect];
      marker_travail[i] = buffer_marker[(int)i-nbrdetect];
      idx_travail[i] = buffer_idx[(int)i-nbrdetect];
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
      tampon_tm[i-buffersize] = tm_travail_ordered[i];
      tampon_enrj[i-buffersize] = enrj_travail_ordered[i];
      tampon_marker[i-buffersize] = marker_travail_ordered[i];
      tampon_idx[i-buffersize] = idx_travail_ordered[i];
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
