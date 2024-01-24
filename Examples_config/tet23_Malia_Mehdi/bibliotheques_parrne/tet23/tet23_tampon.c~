#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "../includes/comet_var_.h"
#include "../includes/comet_proto_.h"
#include "../includes/comvisu.h"
#include "../includes/compteurs.h"
#include <byteswap.h>

#include "../includes/evenements.h"

#define NTRANCHES 10
#define NVOIES 6

extern MODE_RUN *modeRun;
struct groupe
{
  unsigned long nombre_total_coincidences_doubles;
  unsigned long nombre_total_coincidences_doubles_pretraites;
  unsigned long nombre_total_evenements_simples;
  unsigned long nombre_total_coincidences_doubles_dernier_tampon;
  unsigned long nombre_total_evenements_simples_voie0;
  unsigned long nombre_total_evenements_simples_voie0_reconstruits;
  unsigned long nombre_total_evenements_simples_voie0_reconstruits_b1;
  unsigned long nombre_total_evenements_simples_voie0_reconstruits_b2;
  struct evenement_reconstruit *evenements_precedents;
  unsigned long nombre_evenements_precedents;
  unsigned char appliquer_algorithme_tableau_zero;
  unsigned long numero_tampon_courant;
  unsigned char correction_algorithme;
};

#define NOM_MANIP "TETRA23"
#define NOM_REPERTOIRE_MANIPS "/home/verney/manip"

struct argument_permanent
{
  long taille;
  struct groupe groupe1;
  struct groupe groupe2;
  double longueur_onde_temp;
  double longueur_onde;
  long numdeb_laser0;
  long numdeb_laser12;
  unsigned long evenement_bizarre_present;
};

void trait_table1(struct argument_permanent *,unsigned short*,int);
void utrait_uti_acq_2g(unsigned short int *,int *);
void finrun_uti_acq_2g(void);
void debrun_uti_acq_2g(void);
void utrait_acq_1g(struct argument_permanent *,
		   unsigned short int *,int,int,struct groupe *);


void *fonction_initialisation(void)
{
  struct argument_permanent *ptr;

  ptr = malloc(sizeof(struct argument_permanent));
  ptr->longueur_onde_temp = 0.0;
  ptr->longueur_onde= 0.0;
  ptr->numdeb_laser0=206;
  ptr->numdeb_laser12=236;
  ptr->evenement_bizarre_present = 0;
  ptr->groupe1.evenements_precedents = NULL;
  ptr->groupe2.evenements_precedents = NULL;
  printf ("appel �init_manip ");
  init_manip(NOM_MANIP,NOM_REPERTOIRE_MANIPS);
  printf ("ok\n");
  return (void *) ptr;
}

void sur_demarrer(struct argument_permanent *ptr,int numero_run)
{
  int i;

  /* Initialisation des  multispectres laser */
  ptr->longueur_onde_temp = 0.0;
  ptr->longueur_onde = 0.0;
  ptr->numdeb_laser0 =206;  /* tests laser sur voie 0 et voie 12*/
  ptr->numdeb_laser12 =236; 

  ptr->groupe1.nombre_total_coincidences_doubles = 0;
  ptr->groupe1.nombre_total_coincidences_doubles_pretraites = 0;
  ptr->groupe1.nombre_total_evenements_simples = 0;
  ptr->groupe1.nombre_total_coincidences_doubles_dernier_tampon = 0;
  ptr->groupe1.nombre_total_evenements_simples_voie0 = 0;
  ptr->groupe1.nombre_total_evenements_simples_voie0_reconstruits = 0;
  ptr->groupe1.nombre_total_evenements_simples_voie0_reconstruits_b1 = 0;
  ptr->groupe1.nombre_total_evenements_simples_voie0_reconstruits_b2 = 0;
  ptr->groupe1.appliquer_algorithme_tableau_zero=0;
  ptr->groupe1.numero_tampon_courant = 0;
  if (ptr->groupe1.evenements_precedents != NULL)
    {
      printf("liberation evenements grp 1\n");
      for(i=0;i<ptr->groupe1.nombre_evenements_precedents;i++)
	{
	  if (ptr->groupe1.evenements_precedents[i].base != NULL)
	    {
	      free(ptr->groupe1.evenements_precedents[i].base);
	    }
	}
      if (ptr->groupe1.evenements_precedents != NULL)
	{
	  free(ptr->groupe1.evenements_precedents);
	}
    }
  ptr->groupe1.evenements_precedents = NULL;
  ptr->groupe1.nombre_evenements_precedents=0;

  ptr->groupe2.nombre_total_coincidences_doubles = 0;
  ptr->groupe2.nombre_total_coincidences_doubles_pretraites = 0;
  ptr->groupe2.nombre_total_evenements_simples = 0;
  ptr->groupe2.nombre_total_coincidences_doubles_dernier_tampon = 0;
  ptr->groupe2.nombre_total_evenements_simples_voie0 = 0;
  ptr->groupe2.nombre_total_evenements_simples_voie0_reconstruits = 0;
  ptr->groupe2.nombre_total_evenements_simples_voie0_reconstruits_b1 = 0;
  ptr->groupe2.nombre_total_evenements_simples_voie0_reconstruits_b2 = 0;
  ptr->groupe2.appliquer_algorithme_tableau_zero=0;
  ptr->groupe2.numero_tampon_courant = 0;
  if (ptr->groupe2.evenements_precedents)
    {
      for(i=0;i<ptr->groupe1.nombre_evenements_precedents;i++)
	{
	  if (ptr->groupe2.evenements_precedents[i].base != NULL)
	    {
	      free(ptr->groupe2.evenements_precedents[i].base);
	    }
	}
      if (ptr->groupe2.evenements_precedents != NULL)
	{
	  printf("liberation evenements grp 2\n");
	  free(ptr->groupe2.evenements_precedents);
	}
    }
  ptr->groupe2.evenements_precedents = NULL;
  ptr->groupe2.nombre_evenements_precedents=0;

  ptr->groupe1.correction_algorithme = 1;
  ptr->groupe2.correction_algorithme = 1;

  cpttable2 = 0;
  ptr->taille= 0;
  cptes = 0;
  cpted = 0;
  cpt1 = 0;
  cpt2 = 0;
  ptr->evenement_bizarre_present = 0;
 
  cptevtdtt = 0;
  for(i = 0;i<20;i++)
    tt[i] = 0;

  printf ("numero run :%d\n",numero_run);
  passer_commande(DEBRUN,numero_run,0,0,0,"toto",0);
  debrun_uti_acq_2g();
 return;
}

/* ================================================ */
/*  Bilan du traitement des evts simples et doubles  */
/* ================================================ */	
void sur_arreter(struct argument_permanent *ptr)
{

  finrun_uti_acq_2g();

  printf("groupe1.nombre_total_coincidences_doubles= %ld\n",
	 ptr->groupe1.nombre_total_coincidences_doubles);
  printf("groupe1.nombre_total_coincidences_doubles_pretraites= %ld\n",
	  ptr->groupe1.nombre_total_coincidences_doubles_pretraites);
  printf("groupe1.nombre_total_coincidences_doubles_dernier"
	  "_tampon= %ld\n",
	  ptr->groupe1.nombre_total_coincidences_doubles_dernier_tampon);
  printf("groupe1.nombre_total_evenements_simples= %ld\n",
	  ptr->groupe1.nombre_total_evenements_simples);

  printf("groupe2.nombre_total_coincidences_doubles= %ld\n",
	 ptr->groupe2.nombre_total_coincidences_doubles);
  printf("groupe2.nombre_total_coincidences_doubles_pretraites= %ld\n",
	  ptr->groupe2.nombre_total_coincidences_doubles_pretraites);
  printf("groupe2.nombre_total_coincidences_doubles_dernier"
	  "_tampon= %ld\n",
	  ptr->groupe2.nombre_total_coincidences_doubles_dernier_tampon);
  printf("groupe2.nombre_total_evenements_simples= %ld\n",
	  ptr->groupe2.nombre_total_evenements_simples);
  return;
}
/********************************************/
/* traitement et constitution des spectres  */
/* sans tri et coincidences                 */
/* a partir du buffer event                 */
/********************************************/
void trait_table1(struct argument_permanent *ptr,
		  unsigned short *event,
		  int lenevent)
{
  int i = 0,rc = 0;
  unsigned short sevent[8];
  int flag = 0;
 
  /*printf("**************trait_table1 lenevent = %d\n",lenevent);*/

  for(i=0;i<lenevent;)
    {
      if(event[i+1] == 8)/* traitement des evts simples */
	{
	  if(flag == 0)
	    {
	      /*	      printf("event[0] =  %x\n",(int)event[0]); 
	     printf("event[1] =  %x\n",(int)event[1]); 
	     printf("event[2] =  %x\n",(int)event[2]);
	     printf("event[3] =  %x\n",(int)event[3]); 
	     printf("event[4] =  %x\n",(int)event[4]); 
	     printf("event[5] =  %x\n",(int)event[5]); 
	     printf("event[6] =  %x\n",(int)event[6]); 
	     printf("event[7] =  %x\n",(int)event[7]);
	     flag =1;*/
	     /*printf("\n\n");*/
	     /* printf("event[0] =  %x\n",(int)event[i]); 
	     printf("event[1] =  %x\n",(int)event[i+1]); 
	     printf("event[2] =  %x\n",(int)event[i+2]);
	     printf("event[3] =  %x\n",(int)event[i+3]); 
	     printf("event[4] =  %x\n",(int)event[i+4]); 
	     printf("event[5] =  %x\n",(int)event[i+5]); 
	     printf("event[6] =  %x\n",(int)event[i+6]); 
	     printf("event[7] =  %x\n",(int)event[i+7]);*/
	     }
	  sevent[rc]=event[i];
	  rc++;
	  /*printf("rc= %d, i=%d", rc, i);*/
	  sevent[rc] = event[i+1];
	  rc++;	 
	  sevent[rc] = event[i+2];
	  rc++;	 
	  sevent[rc] = event[i+3];
	  rc++;	 
	  sevent[rc] = event[i+4];
	  rc++;	 
	  sevent[rc] = event[i+5];
	  rc++;	 
	  sevent[rc] = event[i+6];
	  rc++;	 
	  sevent[rc] = event[i+7];
	  rc++;
	  if ((event[i+2] & 0x3f) == 0x2f)
	    cptevtder++;  
	  /*if(flag == 0)
	    {
	     printf("sevent[0] =  %x\n",(int)sevent[0]); 
	     printf("sevent[1] =  %x\n",(int)sevent[1]); 
	     printf("sevent[2] =  %x\n",(int)sevent[2]);
	     printf("sevent[3] =  %x\n",(int)sevent[3]); 
	     printf("sevent[4] =  %x\n",(int)sevent[4]); 
	     printf("sevent[5] =  %x\n",(int)sevent[5]); 
	     printf("sevent[6] =  %x\n",(int)sevent[6]); 
	     printf("sevent[7] =  %x\n",(int)sevent[7]);
	     flag = 1;
	     }*/
	  if(rc == 8)
	  {
	  
	    /* on envoie les evts simples un par un a utrait OASIS */
	    utrait_uti_acq_2g(sevent,&rc);
	    rc= 0;
	    cptes++;
	    cptevttot++;
	  } 
	  i = i+8;
	}
      else
	{
	  if(event[i+1] == 14)
	    {
	      /* on saute les evts doubles s'ils existent */	  
	      i = i+14;
	      cpted++;
	      cptevttot++;
	      /*printf("dbl");*/
	    }
	  else
	    {
	      if(event[i+1] == 5)
		{
		  /* premier evt qui contient les no de slots */
		  sevent[rc]=event[i];
		  rc++;
		  /*printf("rc= %d, i=%d", rc, i);*/
		  sevent[rc] = event[i+1];
		  rc++;	 
		  sevent[rc] = event[i+2];
		  rc++;	 
		  sevent[rc] = event[i+3];
		  rc++;	 
		  sevent[rc] = event[i+4];
		  rc++;	 
		  utrait_uti_acq_2g(sevent,&rc);
		  rc= 0;
		  i = i+5;
		}
	      else
		{
		  if(event[i+1] == EVT_LASER)
		    {
		      unsigned short int tmp_event_laser[13];
		      int taille_event,loc_i;

		      taille_event = EVT_LASER;
		      for(loc_i=0;loc_i<EVT_LASER;loc_i++)
			{
			  tmp_event_laser[loc_i] = event[i+loc_i];
			}
		      utrait_uti_acq_2g(tmp_event_laser,&taille_event);
		      i += EVT_LASER;
		    }
		  else
		    {
		      printf("evt bizz event[%d] = %x \n",i,event[i]);
		      printf("evt bizz event[%d]+1 = %x \n",i+1,event[i+1]);
		      printf("evt bizz event[%d]+1 = %x \n",i+2,event[i+2]);
		      printf("evt bizz event[%d]+1 = %x \n",i+3,event[i+3]);
		      printf("evt bizz event[%d]+1 = %x \n",i+4,event[i+4]);
		      printf("evt bizz event[%d]+1 = %x \n",i+5,event[i+5]);
		      i++;
		      /* modif 08 juillet 2003 : stop messages */
		      ptr->evenement_bizarre_present = 1;
		      return;
		    }
		}
	    }
	}
    }
}

/* ======================================== */
/*  Traitement des evts simples et doubles  */
/* ======================================== */	

void travail_en_acquisition(struct argument_permanent *ptr,
			    unsigned short int *event,
			    int	lenevent_octet,
			    int swap)
{
  unsigned short int *event1,*event2,desc;
  int lenevent1=0,lenevent2=0,i,j1=0,j2=0,lenevent,lenevent3=0;

  lenevent = lenevent_octet / 2;

  /* sd mis en commentaire car arrete acquisition si presence evt bizarre 
  if (ptr->evenement_bizarre_present)
    return;
  */
  if (swap == 1)
    {
      for(i=0;i<lenevent;i++)
	{
	  event[i] = bswap_16(event[i]);
	}
    }
  i=0;
  while(1)
    {
      if ((event[i+1]%2) != 0)
	{
	  utrait_acq_1g(ptr,&event[i],event[i+1],swap,&ptr->groupe1);
	  if (event[i+1] == 5)
	    printf("evenement description chassis vu\n");
	  i+= event[i+1];
	}
      else
	{
	  /*	  if (swap)
	    {
	      desc = bswap_16(event[i+2]);
	    }
	  else
	  {*/
	      desc = event[i+2];
	      /*}*/
	  desc = (desc & 0xf000) >> 12;
	  if (desc == 0)
	    {
	      lenevent1++;
	      /*printf("groupe 0 %d %x\n",(int)desc,(int)desc);*/
	    }
	  else
	    {
	      if (desc == 1)
		{
		  lenevent2++;
		  /* printf("groupe 1 %d %x\n",(int)desc,(int)desc);*/
		}
	      else
		{
		  lenevent3++;
		  printf("groupe non prevu %d %x\n",(int)desc,(int)desc);
		  printf ("event[0] %x\n",event[i]);
		  printf ("event[1] %x\n",event[i+1]);
		  printf ("event[2] %x\n",event[i+2]);
		  printf ("event[3] %x\n",event[i+3]);
		  printf ("event[4] %x\n",event[i+4]);
		  printf ("event[5] %x\n",event[i+5]);
		  printf ("event[6] %x\n",event[i+6]);
		  printf ("event[7] %x\n",event[i+7]);
		}
	    }
	  i+= 8;
	}
      if (i >= lenevent)
	break;
    }
  /*  printf("lenevent1 %d lenevent2 %d lenevent3 %d\n",
      lenevent1,lenevent2,lenevent3); */
  lenevent1 *= 8;
  lenevent2 *= 8;
  /*  if ((lenevent1+lenevent2) != lenevent)
    {
      printf("il y a eu des groupes non prevus\n");
      printf("pas de traitement\n");
      return;
      }*/
  event1 = (unsigned short *) calloc(lenevent1,sizeof(unsigned short));
  event2 = (unsigned short *) calloc(lenevent2,sizeof(unsigned short));
  i=0;
  while(1)
    {
      if ((event[i+1]%2) != 0)
	{
	  i+= event[i+1];
	}
      else
	{
	  /*	  if (swap)
	    {
	      desc = bswap_16(event[i+2]);
	    }
	  else
	  {*/
	      desc = event[i+2];
	      /*}*/
	  desc = (desc & 0xf000) >> 12;
	  if (desc == 0)
	    {
	      event1[j1] = event[i];
	      event1[j1+1] = event[i+1];
	      event1[j1+2] = event[i+2];
	      event1[j1+3] = event[i+3];
	      event1[j1+4] = event[i+4];
	      event1[j1+5] = event[i+5];
	      event1[j1+6] = event[i+6];
	      event1[j1+7] = event[i+7];
	      j1+=8;
	    }
	  if (desc == 1)
	    {
	      event2[j2] = event[i];
	      event2[j2+1] = event[i+1];
	      event2[j2+2] = event[i+2];
	      event2[j2+3] = event[i+3];
	      event2[j2+4] = event[i+4];
	      event2[j2+5] = event[i+5];
	      event2[j2+6] = event[i+6];
	      event2[j2+7] = event[i+7];
	      j2+=8;
	    }
	  i += 8;
	}
      if (i >= lenevent)
	break;
    }
  utrait_acq_1g(ptr,event1,lenevent1,swap,&ptr->groupe1);
  utrait_acq_1g(ptr,event2,lenevent2,swap,&ptr->groupe2);
  free(event1);
  free(event2);
}

void utrait_acq_1g(struct argument_permanent *ptr,
	    unsigned short int *event,
	    int	lenevent,
	    int swap,
	    struct groupe *grp)
{
  int i = 0;
  unsigned char debug_xavier = 0;

  if (lenevent == 0)
    return;

  /* if (swap == 1)
    {
      for(i=0;i<lenevent;i++)
	{
	  event[i] = bswap_16(event[i]);
	}
    }*/
  ptr->taille = 0;
  grp->numero_tampon_courant++;

  if ((grp->numero_tampon_courant % 500) == 0)
    {
      printf("numero_tampon_courant %ld %ld\n",grp->numero_tampon_courant,
	     grp->nombre_total_coincidences_doubles);
    }
  if(modeRun->coincidence)
    {
      unsigned long *evenement_physiques_tableau,evenement_temps_mort=0;
      unsigned long evenement_code_en=0,evenement_double=0;
      struct evenement_reconstruit **tableau_evenements;
      unsigned long nb_code_en_max = 1+lenevent/8;
      unsigned long nb_tableaux,tableau_courant,numero_evenement_reconstruit;
      struct description_event_acq *tmp_event;
      unsigned long evenement_physiques_totaux=0;
      int dernier_tableau=-1;

      evenement_physiques_tableau = (unsigned long *) calloc
	(nb_code_en_max,sizeof(unsigned long));
      for(i=0;i<nb_code_en_max;i++)
	{
	  evenement_physiques_tableau[i] = 0;
	}
      /*      printf(" ************Mode coinc lenevent = %d\n",lenevent);*/
      i = 0;
  
      while(i < lenevent)
	{
	  /*printf("event[0] =  %x\n",(int)event[0]); 
	     printf("event[1] =  %x\n",(int)event[1]); 
	     printf("event[2] =  %x\n",(int)event[2]);
	     printf("event[3] =  %x\n",(int)event[3]); 
	     printf("event[4] =  %x\n",(int)event[4]); 
	     printf("event[5] =  %x\n",(int)event[5]); 
	     printf("event[6] =  %x\n",(int)event[6]); 
	     printf("event[7] =  %x\n",(int)event[7]);*/
	  if(event[i+1] == EVTS_SIMPLES_ACQ)
	    {
	      tmp_event = (struct description_event_acq *) &event[i];
	      grp->nombre_total_evenements_simples++;
	      if ((tmp_event->desc & 0x3f) <= 0x15)
		{
		  if ((tmp_event->desc & 0x3f) == 0x0)
		    {
		      grp->nombre_total_evenements_simples_voie0++;
		    }
		  evenement_physiques_tableau[evenement_code_en]++;
		  evenement_physiques_totaux++;
		}
	      else
		{
		  int taille_event=/*7*/8,loc_i;
		  /*unsigned short int tmp_event_oasis[7];*/
		  unsigned short int tmp_event_oasis[8];
	  
		  /*tmp_event_oasis[0] = 0xffff;*/
		  /*for(loc_i=0;loc_i<6;loc_i++)*/
		  for(loc_i=0;loc_i<8;loc_i++)
		    {
		      /*tmp_event_oasis[loc_i+1] = event[i+loc_i]*/
		      tmp_event_oasis[loc_i] = event[i+loc_i];
		    }
		  utrait_uti_acq_2g(tmp_event_oasis,&taille_event);
		  if ((tmp_event->desc & 0x3f) >= 0x20 &&
		      ((tmp_event->desc & 0x3f) <= 0x25))
		    {
		      evenement_temps_mort++;
		      cptevttot++;
		    }
		  else
		    {
		      if ((tmp_event->desc & 0x3f) == 0x26)
			{
			  evenement_code_en++;
			  cptevttot++;
			}
		      else
			{
			  if ((tmp_event->desc & 0x3f) == 0x2f)
			    {
			      printf("%ld evenement fin run\n",
				     grp->numero_tampon_courant);
			      cptevtder++;
			      cptevttot++;
			    }
			  else
			    {
			      printf("%ld evenement service %d\n",
				     grp->numero_tampon_courant,
				     tmp_event->desc);
			    }
			}
		    }
		}
	      i += 8;
	    }
	  else
	    {
	      if(event[i+1] == 14)
		/* on saute les evenements doubles */
		/* on les recalcule apres          */
		{
		  i = i + 14;
		  evenement_double++;		  
		}
	      else
		{
		  /* premier evt qui contient les no de slot  */
		  int loc_i,taille_event = 5;
		  unsigned short int tmp_event_oasis[5];
		  unsigned short int tmp_event_laser[13];
		  
		  if(event[i+1] == 5)
		    {
		      for(loc_i=0;loc_i<5;loc_i++)
			{
			  tmp_event_oasis[loc_i] = event[i+loc_i];
			}
		      utrait_uti_acq_2g(tmp_event_oasis,&taille_event);
		      i = i + 5;
		    }
		  else
		    {
		      if(event[i+1] == EVT_LASER)
			{
			  taille_event = EVT_LASER;
			  for(loc_i=0;loc_i<EVT_LASER;loc_i++)
			    {
			      tmp_event_laser[loc_i] = event[i+loc_i];
			    }
			  utrait_uti_acq_2g(tmp_event_laser,&taille_event);
			  i += EVT_LASER;
			}
		      else
			{
			  printf("111 evt bizz event[%d] = %x \n",i,event[i]);
			  printf("111 evt bizz event[%d]+1 = %x \n",i+1,event[i+1]);
			  printf("111 evt bizz event[%d]+1 = %x \n",i+2,event[i+2]);
			  printf("111 evt bizz event[%d]+1 = %x \n",i+3,event[i+3]);
			  printf("111 evt bizz event[%d]+1 = %x \n",i+4,event[i+4]);
			  printf("111 evt bizz event[%d]+1 = %x \n",i+5,event[i+5]);
			  i++;
			  ptr->evenement_bizarre_present = 1;
			  return;
			}
		    }		  
		}
	    }
	}

      grp->nombre_total_coincidences_doubles_pretraites += evenement_double;
      nb_tableaux = evenement_code_en+1;
      if ((evenement_code_en != 0) && (debug_xavier == 1))
	{
	  printf("nb evenements physiques totaux %ld tm %ld codeen %ld"
		 " nb tableaux %ld\n",
		 evenement_physiques_totaux,
		 evenement_temps_mort,
		 evenement_code_en,
		 nb_tableaux);
	}
      if (evenement_physiques_totaux == 0)
	return;
      tableau_evenements = (struct evenement_reconstruit **)
	calloc(nb_tableaux,sizeof(struct evenement_reconstruit *));
 
      for(i=0;i<nb_tableaux;i++)
	{
	  int loc_i;

	  if ((i==0) && (grp->appliquer_algorithme_tableau_zero == 0))
	    {
	      evenement_physiques_tableau[i] +=
		grp->nombre_evenements_precedents;
	    }
	  if (evenement_physiques_tableau[i] == 0)
	    {
	      tableau_evenements[i] = (struct evenement_reconstruit *) NULL;
	    }
	  else
	    {
	      tableau_evenements[i] = (struct evenement_reconstruit *)
		calloc(evenement_physiques_tableau[i],
		       sizeof(struct evenement_reconstruit));
	      if (tableau_evenements[i] == NULL)
		{
		  printf("erreur d'allocation\n");
		}
	    }

	  for(loc_i=0;loc_i<evenement_physiques_tableau[i];loc_i++)
	    {
	      tableau_evenements[i][loc_i].base= NULL;
	      tableau_evenements[i][loc_i].nb_coinc_double = 0;
	      tableau_evenements[i][loc_i].temps = 0;
	      tableau_evenements[i][loc_i].traite = 0;
	    }

	  if ((i==0) && (grp->appliquer_algorithme_tableau_zero == 0))
	  {
	    for(loc_i=0;loc_i<grp->nombre_evenements_precedents;loc_i++)
	      {
		recopie_evenements_reconstruits
		  (&tableau_evenements[i][loc_i],
		   &grp->evenements_precedents[loc_i]);
	      }
	  }
	}
      i=0;
      tableau_courant = 0;
      if (grp->appliquer_algorithme_tableau_zero == 1)
	{
	  numero_evenement_reconstruit = 0;
	}
      else
	{
	  numero_evenement_reconstruit = grp->nombre_evenements_precedents;
	}
      while(i < lenevent)
	{
	  if(event[i+1] == EVTS_SIMPLES_ACQ)
	    {
	      tmp_event = (struct description_event_acq *) &event[i];
	      if ((tmp_event->desc & 0x3f) <= 0x15)
		{
		  unsigned long long int temps_reconstruit;
		  struct evenement_reconstruit *tmp_struct;
		  int j;

		  if ((tmp_event->desc & 0x3f) == 0x0)
		    {
		      grp->nombre_total_evenements_simples_voie0_reconstruits++;
		    }
		  /* if (tmp_event->temps_fort != 0) */
/* 		    { */
/* 		      printf("coucou temps fort %x\n",tmp_event->temps_fort); */
/* 		    } */
		  temps_reconstruit =
		    (((unsigned long long int)
		      (tmp_event->temps_fort & 0x7fff)) << 32) |
		    (((unsigned long long int) tmp_event->temps_moyen << 16)) |
		    ((unsigned long long int) tmp_event->temps_faible);

		  if (debug_xavier == 2)
		    {
		      /*affiche_unsigned_long_long(temps_reconstruit);*/
		    }
		  /* printf("argl %ld %ld ",tableau_courant, */
/* 			 numero_evenement_reconstruit); */
		  tmp_struct = &tableau_evenements[tableau_courant]
		    [numero_evenement_reconstruit];
		  tmp_struct->temps = temps_reconstruit;
		  /*tmp_struct->base = (unsigned short *) calloc
		    (6,sizeof(unsigned short));*/
		  tmp_struct->base = (unsigned short *) calloc
		    (8,sizeof(unsigned short));
		  /*for(j=0;j<6;j++)*/
		  for(j=0;j<8;j++)
		    {
		      tmp_struct->base[j] = event[i+j];
		    }
 		  tmp_struct->nb_coinc_double = 0;
		  numero_evenement_reconstruit++;
		  /* printf("fin argl\n"); */
		}
	      else
		{
		  unsigned long long int temps_reconstruit;

		  temps_reconstruit =
		    (((unsigned long long int)
		      (tmp_event->temps_fort & 0x7fff)) << 32) |
		    (((unsigned long long int) tmp_event->temps_moyen << 16)) |
		    ((unsigned long long int) tmp_event->temps_faible);
		  if ((tmp_event->desc & 0x3f) == 0x26)
		    {
		      tableau_courant++;
		      /* if ((tmp_event->desc & 0xf3f) == 0x226)
			{
			  printf("temps2 %x %x %x ",tmp_event->temps_fort,tmp_event->temps_moyen,tmp_event->temps_faible);
			  affiche_unsigned_long_long(temps_reconstruit);
			  }*/
		      numero_evenement_reconstruit = 0;
		    }
		}
	      i+=8;
	    }
	  else
	    {
	      if(event[i+1] == 14)
		/* on saute les evenements doubles */
		/* on les recalcule apres          */
		{
		  i = i + 14;
		}
	      else
		{
		  if(event[i+1] == 5)
		    {
		      i = i+5;
		    }
		  else
		    {
		      if(event[i+1] == EVT_LASER)
			{
			  i += EVT_LASER;
			}
		      else
			{
			  printf("222 evt bizz event[%d] = %x \n",i,event[i]);
			  i++;
			  ptr->evenement_bizarre_present = 1;
			  return;
			}
		    }
		}
	    }
	}  
      if (grp->appliquer_algorithme_tableau_zero == 1)
	{
	  unsigned long long int temps0,temps1;
	  long long int dt;
	  int loc_i,temp_i,transition=-2;
	  struct evenement_reconstruit *tableau_temporaire;
	  unsigned long int retour;
	  unsigned long retour_calcul_coinc;

	  if (grp->correction_algorithme == 1)
	    {
	      temps0 = grp->evenements_precedents
		[grp->nombre_evenements_precedents-1].temps;
	      temps1 = tableau_evenements[0][0].temps;
	      dt = (long long int)temps1-(long long int)temps0;
	      if (dt < 0)
		{
		  transition = -1;
		}
	      else
		{
		  for(loc_i=0;loc_i<evenement_physiques_tableau[0]-1;loc_i++)
		    {
		      temps0 = tableau_evenements[0][loc_i].temps;
		      temps1 = tableau_evenements[0][loc_i+1].temps;
		      dt = (long long int)temps1-(long long int)temps0;
		      if (dt < 0)
			{
			  transition = loc_i;
			  break;
			}
		    }
		}
	      if (transition == -2)
		{
		  printf("%ld bizarre, pas de transition\n",
			 grp->numero_tampon_courant);
		  printf("nb de tableaux %ld\n",nb_tableaux);
		  for(loc_i=0;loc_i<nb_tableaux;loc_i++)
		    {
		      printf("nb evenements physiques %ld\n",
			     evenement_physiques_tableau[loc_i]);
		    }
		  if (evenement_physiques_tableau[0] < 5)
		    {
		      printf("transition forcee nb evenements physiques %ld\n",
			     evenement_physiques_tableau[0]);
		      transition = evenement_physiques_tableau[0]-1;
		    }
		}
	      grp->evenements_precedents = (struct evenement_reconstruit *)
		realloc(grp->evenements_precedents,
			(grp->nombre_evenements_precedents+transition+1)*
			sizeof(struct evenement_reconstruit));
	      
	      for(loc_i=0;loc_i<transition+1;loc_i++)
		{
		  temp_i = grp->nombre_evenements_precedents + loc_i;
		  recopie_evenements_reconstruits
		    (&grp->evenements_precedents[temp_i],
		     &tableau_evenements[0][loc_i]);
		}
	      
	      tableau_temporaire = tableau_evenements[0];
	      tableau_evenements[0] = (struct evenement_reconstruit *)
		calloc(evenement_physiques_tableau[0]-(transition+1),
		       sizeof(struct evenement_reconstruit));
	      
	      for(loc_i=0;
		  loc_i<evenement_physiques_tableau[0]-(transition+1);
		  loc_i++)
		{
		  recopie_evenements_reconstruits
		    (&tableau_evenements[0][loc_i],
		     &tableau_temporaire[loc_i+(transition+1)]);
		}	      
	      evenement_physiques_tableau[0] -= (transition+1);
	      grp->nombre_evenements_precedents += (transition+1);
	      free(tableau_temporaire);
	    }
	  retour = tri_bulle(grp->evenements_precedents,
			     grp->nombre_evenements_precedents);
	  if (retour != 0)
	    {
	      printf("algo applique tampon zero %ld\n",retour);
	    }
	  retour_calcul_coinc = coincidences_doubles
	    (grp->evenements_precedents,
	     grp->nombre_evenements_precedents,
	     0,(unsigned long long int) modeRun->fenetreCoinc);
	  grp->nombre_total_coincidences_doubles += retour_calcul_coinc;
	  if ((retour_calcul_coinc != 0) && (debug_xavier >= 1))
	    {
	      printf("%ld algo_tab_zero : nb de coincidences %ld\n",
		     grp->numero_tampon_courant,retour_calcul_coinc);
	    }
	  for(i=0;i<grp->nombre_evenements_precedents;i++)
	    {
	      if (grp->evenements_precedents[i].traite == 0)
		{
		  unsigned short *local_event;
		  /*int taille_evenement_local=7,j;*/
		  int taille_evenement_local=8,j;

		  grp->evenements_precedents[i].traite = 1;
		  /*local_event = (unsigned short *) calloc
		    (7,sizeof(unsigned short));*/
		  local_event = (unsigned short *) calloc
		    (8,sizeof(unsigned short));
		  /*local_event[0] = 0xffff;*/
		  /*for(j=0;j<6;j++)*/
		  for(j=0;j<8;j++)
		    {
		      /* local_event[j+1] = evenements_precedents[i].base[j];*/
		      local_event[j] = grp->evenements_precedents[i].base[j];
		    }
		  if ((local_event[2] & 0x3f) == 0)
		    grp->nombre_total_evenements_simples_voie0_reconstruits_b1++;
		  utrait_uti_acq_2g(local_event,&taille_evenement_local);
		  cptevttot++;
		  free(local_event);
		}
	      if (grp->evenements_precedents[i].nb_coinc_double != 0)
		{
		  int j;

		  for(j=0;j<grp->evenements_precedents[i].nb_coinc_double;j++)
		    {
		      unsigned short *local_event_double;
		      int taille_evenement_local_double=0;

		      /*local_event_double = (unsigned short *) calloc
			(12,sizeof(unsigned short));*/
		      local_event_double = (unsigned short *) calloc
			(14,sizeof(unsigned short));

		      local_event_double[0] =
			grp->evenements_precedents[i].base[0];
		      local_event_double[1] =
			grp->evenements_precedents[i].base[1];
		      local_event_double[2] =
			grp->evenements_precedents[i].base[2];
		      local_event_double[3] =
			grp->evenements_precedents[i].base[3];
		      local_event_double[4] =
			grp->evenements_precedents[i].base[4];
		      local_event_double[5] =
			grp->evenements_precedents[i].base[5];
		      local_event_double[6] =
			grp->evenements_precedents[i].base[6];
		      local_event_double[7] =
			grp->evenements_precedents[i].base[7];
		      local_event_double[8] =
			grp->evenements_precedents[i+j+1].base[2];
		      local_event_double[9] =
			grp->evenements_precedents[i+j+1].base[3];
		      local_event_double[10] =
			grp->evenements_precedents[i+j+1].base[4];
		      local_event_double[11] =
			grp->evenements_precedents[i+j+1].base[5];
		      local_event_double[12] =
			grp->evenements_precedents[i+j+1].base[6];
		      local_event_double[13] =
			grp->evenements_precedents[i+j+1].base[7];
		      

		      utrait_uti_acq_2g(local_event_double,
				 &taille_evenement_local_double);
		      cptevttot++;
		      free(local_event_double);
		    }
		}
	    }
	}
  
      for(i=0;i<nb_tableaux;i++)
	{
	  unsigned long int retour;
	  
	  if (evenement_physiques_tableau[i] == 0)
	    continue;
	  if ((grp->correction_algorithme == 1) && (i != 0)
	      && (evenement_physiques_tableau[i-1] != 0))
	    {
	      unsigned long long int temps0,temps1;
	      long long int dt;
	      int loc_i,temp_i,transition=-2;
	      struct evenement_reconstruit *tableau_temporaire;

	      temps0 = tableau_evenements[i-1]
		[evenement_physiques_tableau[i-1]-1].temps;
	      temps1 = tableau_evenements[i][0].temps;
	      dt = (long long int)temps1-(long long int)temps0;
	      if (dt < 0)
		{
		  transition = -1;
		}
	      else
		{
		  for(loc_i=0;loc_i<evenement_physiques_tableau[i]-1;loc_i++)
		    {
		      temps0 = tableau_evenements[i][loc_i].temps;
		      temps1 = tableau_evenements[i][loc_i+1].temps;
		      dt = (long long int)temps1-(long long int)temps0;
		      if (dt < 0)
			{
			  transition = loc_i;
			  break;
			}
		    }
		}
	      if (transition == -2)
		{
		  if (evenement_physiques_tableau[i] >= 5)
		    {
		      printf("%ld bizarre, pas de transition\n",
			     grp->numero_tampon_courant);
		      printf("nb de tableaux %ld\n",nb_tableaux);
		      for(loc_i=0;loc_i<nb_tableaux;loc_i++)
			{
			  printf("nb evenements physiques %ld\n",
				 evenement_physiques_tableau[loc_i]);
			}
		    }
		  if (evenement_physiques_tableau[i] < 5)
		    {
		      transition = evenement_physiques_tableau[i]-1;
		    }
		}
	      tableau_evenements[i-1] = (struct evenement_reconstruit *)
		realloc(tableau_evenements[i-1],
			(evenement_physiques_tableau[i-1]+transition+1)*
			sizeof(struct evenement_reconstruit));
	      for(loc_i=0;loc_i<transition+1;loc_i++)
		{
		  temp_i = evenement_physiques_tableau[i-1] + loc_i;
		  recopie_evenements_reconstruits
		    (&tableau_evenements[i-1][temp_i],
		     &tableau_evenements[i][loc_i]);
		}
	      tableau_temporaire = tableau_evenements[i];
	      tableau_evenements[i] = (struct evenement_reconstruit *)
		calloc(evenement_physiques_tableau[i]-(transition+1),
		       sizeof(struct evenement_reconstruit));
	      for(loc_i=0;
		  loc_i<evenement_physiques_tableau[i]-(transition+1);
		  loc_i++)
		{
		  recopie_evenements_reconstruits
		    (&tableau_evenements[i][loc_i],
		     &tableau_temporaire[loc_i+(transition+1)]);
		}
	      evenement_physiques_tableau[i] -= (transition+1);
	      evenement_physiques_tableau[i-1] += (transition+1);	      
	      free(tableau_temporaire);
	    }
	  if (evenement_physiques_tableau[i] == 0)
	    continue;
	  retour = tri_bulle(tableau_evenements[i],
			     evenement_physiques_tableau[i]);
	  if (retour)
	    {
	      /*printf("%ld tableau %d nombre_de_swap %ld nb event phys %ld\n",
		     grp->numero_tampon_courant,i,
		     retour,evenement_physiques_tableau[i]);*/
	    }
	}

      grp->appliquer_algorithme_tableau_zero = 0;
      for(i=nb_tableaux-1;i>=0;i--)
	{
	  if (evenement_physiques_tableau[i] == 0)
	    {
	      grp->appliquer_algorithme_tableau_zero = 1;
	    }
	  if (evenement_physiques_tableau[i] != 0)
	    {
	      unsigned long long int temps_reference;
	      unsigned long taille_tableau;
	      int loc_i;

	      dernier_tableau = i;
	      if (grp->evenements_precedents != NULL)
		{
		  for(loc_i=0;loc_i<grp->nombre_evenements_precedents;loc_i++)
		    {
		      free(grp->evenements_precedents[loc_i].base);
		    }
		  free(grp->evenements_precedents);
		}
	      taille_tableau = evenement_physiques_tableau[i];
	      temps_reference = tableau_evenements[i][taille_tableau-1].temps;
	      grp->nombre_evenements_precedents = 1;
	      for(loc_i=taille_tableau-2;loc_i>=0;loc_i--)
		{
		  long long int dt;

		  dt = (long long int) temps_reference -
		    (long long int) (tableau_evenements[i][loc_i].temps);
		  if (dt > 5000)
		    break;
		  grp->nombre_evenements_precedents++;
		}
	      grp->evenements_precedents = (struct evenement_reconstruit *)
		calloc(grp->nombre_evenements_precedents,
		       sizeof(struct evenement_reconstruit));
	      for(loc_i=0;loc_i<grp->nombre_evenements_precedents;loc_i++)
		{
		  recopie_evenements_reconstruits
		    (&grp->evenements_precedents[loc_i],
		     &tableau_evenements[i]
		     [taille_tableau-grp->nombre_evenements_precedents+loc_i]);
		  grp->evenements_precedents[loc_i].traite = 1;
		}	
	      break;
	    }
	}
      for(i=0;i<nb_tableaux;i++)
	{
	  unsigned long retour_nb_coinc=0;
	  int ii;

	  if (evenement_physiques_tableau[i] ==0)
	    continue;
	  if (i != dernier_tableau)
	    {
	      retour_nb_coinc = coincidences_doubles
		(tableau_evenements[i],
		 evenement_physiques_tableau[i],
		 0,(unsigned long long int) modeRun->fenetreCoinc);
	      if ((retour_nb_coinc != 0) && (debug_xavier >= 1))
		{
		  printf("%ld nb de coincidences %ld\n",
			 grp->numero_tampon_courant,retour_nb_coinc);
		}
	    }
	  if (i == dernier_tableau)
	    {
	      retour_nb_coinc = coincidences_doubles
		(tableau_evenements[i],
		 evenement_physiques_tableau[i],
		 grp->nombre_evenements_precedents,
		 (unsigned long long int) modeRun->fenetreCoinc);
	      if ((retour_nb_coinc != 0) && (debug_xavier >= 1))
		{
		  printf("%ld nb de coincidences %ld\n",
			 grp->numero_tampon_courant,retour_nb_coinc);
		}
	      grp->nombre_total_coincidences_doubles_dernier_tampon =
		retour_nb_coinc;
	    }
	  for(ii=0;ii<evenement_physiques_tableau[i];ii++)
	    {

	      if (tableau_evenements[i][ii].traite == 0)
		{
		  unsigned short *local_event;
		  /*int taille_evenement_local=7,j;*/
		  int taille_evenement_local=8,j;

		  local_event = (unsigned short *) calloc
		    (8,sizeof(unsigned short));
		  for(j=0;j<8;j++)
		    {
		      local_event[j] = tableau_evenements[i][ii].base[j];
		    }
		  if ((local_event[2] & 0x3f) == 0)
		    grp->nombre_total_evenements_simples_voie0_reconstruits_b2++;
		  utrait_uti_acq_2g(local_event,&taille_evenement_local);
		  cptevttot++;
		  tableau_evenements[i][ii].traite = 1;
		  free(local_event);
		}
	      if (tableau_evenements[i][ii].nb_coinc_double != 0)
		{
		  int j;

		  for(j=0;j<tableau_evenements[i][ii].nb_coinc_double;j++)
		    {
		      unsigned short *local_event_double;
		      /*int taille_evenement_local_double=12;*/
		      int taille_evenement_local_double=14;

		      /*local_event_double = (unsigned short *) calloc
			(12,sizeof(unsigned short));*/
		      local_event_double = (unsigned short *) calloc
			(14,sizeof(unsigned short));

		      local_event_double[0] =
			tableau_evenements[i][ii].base[0];
		      local_event_double[1] =
			tableau_evenements[i][ii].base[1] ;
		      local_event_double[2] =
			tableau_evenements[i][ii].base[2];
		      local_event_double[3] =
			tableau_evenements[i][ii].base[3];
		      local_event_double[4] =
			tableau_evenements[i][ii].base[4];
		      local_event_double[5] =
			tableau_evenements[i][ii].base[5];
		      local_event_double[6] =
			tableau_evenements[i][ii].base[6];
		      local_event_double[7] =
			tableau_evenements[i][ii].base[7];
		      local_event_double[8] =
			tableau_evenements[i][ii+j+1].base[2];
		      local_event_double[9] =
			tableau_evenements[i][ii+j+1].base[3];
		      local_event_double[10] =
			tableau_evenements[i][ii+j+1].base[4];
		      local_event_double[11] =
			tableau_evenements[i][ii+j+1].base[5];
		      local_event_double[12] =
			tableau_evenements[i][ii+j+1].base[6];
		      local_event_double[13] =
			tableau_evenements[i][ii+j+1].base[7];

		      utrait_uti_acq_2g(local_event_double,
				 &taille_evenement_local_double);
		      cptevttot++;
		      free(local_event_double);
		    }
		}
	    }
	  grp->nombre_total_coincidences_doubles += retour_nb_coinc;
	}
 
      for(i=0;i<nb_tableaux;i++)
	{
	  int j;

	  if (tableau_evenements[i] != NULL)
	    {
	      for(j=0;j<evenement_physiques_tableau[i];j++)
		{
		  free(tableau_evenements[i][j].base);
		}
	      free(tableau_evenements[i]);
	    }
	}
      free(tableau_evenements);
      free(evenement_physiques_tableau);
    }
  else
    {
      /*printf("********** Mode  NON coinc lenevent = %d\n",lenevent);*/
      trait_table1(ptr,event,lenevent);
    }
}

void sur_decharger(struct argument_permanent *ptr)
{
  free(ptr);
  fin_manip(NOM_MANIP,NOM_REPERTOIRE_MANIPS);
}
