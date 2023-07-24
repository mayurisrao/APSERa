/* spcgen.c
 *
 * model the sky, model the instrument, generate the measurement set
 *
 * spectra are written at each time assuming sky is stationary in each integration.
 * spectra are written at intervals of INTEGRATION_TIME
 *
 * Ravi: original - 07 Dec 2012
 */

// Comment/uncomment path to nrutil.h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </usr/local/miriad/inc/maxdimc.h>
#include </Users/mayuri/work/C/ANSI-C/nrutil.h>
//#include </Users/rsubrahm/CODESPACE/INCLUDE/nrutil.h>
#include </usr/local/miriad/inc/miriad.h>
#include <string.h>

#define CHANNEL_WIDTH 10.0 /* MHz */
#define NUMBER_OF_CHANNELS 601 /* 601 to give a 6.0 GHz bandwidth */
#define START_FREQUENCY 1.0 /* 3.0 GHz */
#define INTEGRATION_TIME 60.0/* sec 	This is the UT time interval in which
												spectra are written out */

#define NOISE_INT_TIME 7.5e12/* 		2.7e10 from equation 4 in draft v25 for 3 sigma det		(256 elements x 1a years) = 8e9 sec			
												This is the time used for computation
												of noise added to the spectrum ... using time for just over 3 sigma*/ 

#define	SITE_LATITUDE  +13.01333 //-23.0228///* deg */
#define	SITE_LONGITUDE 77.580833 //-67.7550///* deg */
#define site_altitude 200.0 //4800///* meters */
#define NHPIX 786432
#define PI 3.14159265358979
#define NN 5000 // length of array containing recombination line spectra

#define Trx 14.0 // K 
#define T_atm 0.0113044 //K ( atmospheric correction - increase on including absorption and translating T_sys (Trx + atmosphere) to above atmosphere)
#define TCMB 2.72548
#define T_cold 2730.0
#define T_hot 3730.0 
#define d2r PI/180.0
#define kk 1.3806488e-23 // m^2 kg s^-1
#define hh 6.62606957e-34 // m^2 kg s^-2 K^-1
#define cvel 2.99792458e+08 // m s^-1

double cal_lst(double utc_julian_date, double longitude_radians);
void polcof(float xa[], float ya[], int n, float cof[]);
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void precess_(float *ra, float *dec, float *epoch1, float *epoch2);
double beam_definition(float frequency, double azimuth, double elevation);
float gasdev(long *idum); 
float gammq(float a, float x);
void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
        float *b, float *siga, float *sigb, float *chi2, float *q);

/* functions for interpolation of recombination data 
   and conversion between brightness intensity and temperature */
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

/* Function to account for effects of atmospheric refraction */
double refraction (double rad_elevation, double m_altitude);

main()
{
  const char version[30]="spcgen: 28apr14";


  static const char filename_408[] = "../data/all_sky_txt/lambda_haslam408_nofilt_1deg_r8.txt";
  static const char filename_1420[] = "../data/all_sky_txt/1420_1deg_R8.txt";
  static const char filename_23G[] = "../data/all_sky_txt/23GHz_foreground_image.txt";
  static const char filename_coord[] = "../data/all_sky_txt/PIXEL_LISTING_R8.TXT";
  static const char filename_rec[]= "../data/from_Jens/total_spec_new.dat";//all_sky_txt/HI.HeI.HeII.dat

/*
  static const char filename_408[] = "lambda_haslam408_nofilt_1deg_r8.txt";
  static const char filename_1420[] = "1420_1deg_R8.txt";
  static const char filename_23G[] = "23GHz_foreground_image.txt";
  static const char filename_coord[] = "PIXEL_LISTING_R8.TXT";
  static const char filename_rec[]= "HI.HeI.HeII.dat";
*/
  time_t utime;
  long seed;  
  float *xarray,*yarray,*y2;
  float read1,read2;
  float ffmin,ffmax;
  float yp1,ypn;
  double new_yf[NUMBER_OF_CHANNELS];
  double cmb_intensity[NUMBER_OF_CHANNELS],P_cold[NUMBER_OF_CHANNELS],P_hot[NUMBER_OF_CHANNELS];
  double P_diff[NUMBER_OF_CHANNELS];
  double P_diff_temp[NUMBER_OF_CHANNELS];
  float yf;
  int npoints1;
  static const char filename6[]="spec_to_get_template_6nov15.txt";
  static const char filename7[]="spec_index_listing_10may14.txt";
  char outfile[20];
  char inbuf[200];
  int i,j,ii,jj,jjj;
  int nhpix,nbadpixels;
  int npolyfit;
  float read_value,read_value1,read_value2;
  double cwt;
  double cfreq,cfreq1;
  double ctemp;
  float ra_precess,dec_precess,epoch1,epoch2;
  float *ll_coordinate, *bb_coordinate, *sky_408, *sky_1420, *sky_23000;
  float *a_coeff, *b_coeff, *c_coeff;
  float *xx, *yy, *cof, **cc;
  double intensity, cfreq1_hz;
  double final_temp; 
  double final_temp_1, final_temp_2;
  double cof0, cof1, cof2;
  FILE *fptr,*fptr1,*fptr2,*fptr3; 

  /* variables for output vis file */

  char source_name[5];
  char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
  int tno;
  int var_ivalue[1];
  int npoints,nn;
  int antenna1,antenna2;
  int flags[NUMBER_OF_CHANNELS];
  float baseline_value[1];
  double preamble[4];
  double p1,p2,p3;
  double ant_az[2],ant_el[2];
  double freq_channel1[1];
  double freq_inc[1];
  double time_var[1];
  double coord_var[2];
  double site_latitude[1],site_longitude[1];
  float data[2*NUMBER_OF_CHANNELS];
  double sumwt[NUMBER_OF_CHANNELS];
  float var_veldop[1];
  float var_vsource[1];
  double var_restfreq[1];
  double sra[1],sdec[1];
  double lo1[1],lo2[1],freq[1],freqif[1];
  int mount[1];
  float evector[1];
  int nnpol;
  float jyperk[1],inttime[1],epoch[1];
  double antpos[6];
  float tpower[1];
   
  long int jjdd,yyyy,mmmm,dddd,hrhr,minmin;
  long int nspect;
  long int int_time;
  double julian_date, longitude_radians, lst;
  double cross_imag[NUMBER_OF_CHANNELS],cross_real[NUMBER_OF_CHANNELS];
  double sigma; 
  float secsec;
  double c_ll, c_bb;
  double ra_b1950, dec_b1950;
  double ra_date, dec_date;
  double azaz, elel, new_elel;
  double haha;

  double time_utc,time_utc_res;
  int time_utc_hh,time_utc_mm,julian_day;
  float time_utc_sec;	

  ll_coordinate = vector(0,NHPIX);
  bb_coordinate = vector(0,NHPIX);
  sky_408 = vector(0,NHPIX);
  sky_1420 = vector(0,NHPIX);
  sky_23000 = vector(0,NHPIX);
  a_coeff = vector(0,NHPIX);
  b_coeff = vector(0,NHPIX);
  c_coeff = vector(0,NHPIX);
  xx = vector(0,2);
  yy = vector(0,2);
  cof = vector(0,2);
  cc = matrix(0,NHPIX,0,2);

  // read in the sky maps and coordinates; subtract CMB; convert all to K brightness temperature

  fptr=fopen(filename_408,"r");
  i=0;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value); 
      sky_408[i]=(read_value/1000.0)-TCMB;
      ++i;
    }
  fclose(fptr);
  printf(" Read %d values from 408 MHz image \n",i);
		
  fptr=fopen(filename_1420,"r");
  i=0;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value);
      sky_1420[i]=read_value-TCMB;
      ++i;
    }
  fclose(fptr);
  printf(" Read %d values from 1420 MHz image \n",i);
		
  fptr=fopen(filename_23G,"r");
  i=0;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value);
      sky_23000[i]=read_value/1000.0;
      ++i;   
    }
  fclose(fptr);
  printf(" Read %d values from 23000 MHz image \n",i);
		
  fptr=fopen(filename_coord,"r");
  fgets (inbuf, sizeof inbuf, fptr);  // dummy read
  i=0;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g%g",&read_value1,&read_value2);
      ll_coordinate[i]=read_value1*d2r;
      bb_coordinate[i]=read_value2*d2r;
      ++i;
    }
  fclose(fptr);
  printf(" Read %d coordinate values \n",i);

  nhpix = NHPIX;

  // Open and initiatize the measurement set output file 

  /* set all flags to good */
  for(i=0;i<NUMBER_OF_CHANNELS;++i) flags[i]=1;

  npoints=NUMBER_OF_CHANNELS; /* number of spectral points in the records */

  printf("Type output vis file name : ");
  gets(inbuf);
  sscanf(inbuf,"%s",outfile);
	
  uvopen_c(&tno,outfile,"new");
  hisopen_c(tno,"write");

  // setup the output vis file by initializing many variables 

  wrhda_c(tno,"obstype","crosscorrelation");
  uvputvra_c(tno,"source","zenith sky");
  uvputvra_c(tno,"operator","apsera");
  uvputvra_c(tno,"version",version);
  sra[0]=0.0;
  sdec[0]=0.0;
  uvputvrd_c(tno,"ra",sra,1);
  uvputvrd_c(tno,"obsra",sra,1);
  uvputvrd_c(tno,"dec",sdec,1);
  uvputvrd_c(tno,"obsdec",sdec,1);
  lo1[0]=2.000;
  lo2[0]=0.0;
  freq[0]=2.00000;
  freqif[0]=2.00000;
  uvputvrd_c(tno,"lo1",lo1,1);
  uvputvrd_c(tno,"lo2",lo2,1);
  uvputvrd_c(tno,"freq",freq,1);
  uvputvrd_c(tno,"freqif",freqif,1);
  mount[0]=0;
  evector[0]=0.0;
  uvputvri_c(tno,"mount", mount, 1);
  uvputvrr_c(tno,"evector",evector,1);
  uvputvrr_c(tno,"chi",evector,1);
  uvputvra_c(tno,"telescop","apsera");
  jyperk[0]=1.0;
  uvputvrr_c(tno,"jyperk",jyperk,1);
  inttime[0]=(float)INTEGRATION_TIME;
  uvputvrr_c(tno,"inttime",inttime,1);
  epoch[0]=2000.0;
  uvputvrr_c(tno,"epoch",epoch,1);
  nnpol=1;
  wrhdi_c(tno,"npol",nnpol);
  antpos[0]=0.0;
  antpos[1]=0.0;
  antpos[2]=0.0;
  antpos[3]=0.0;
  antpos[4]=1.0;
  antpos[5]=0.0;
  uvputvrd_c(tno,"antpos",antpos,6);
  time_var[0]=2451545.0;
  uvputvrd_c(tno,"time",time_var,1);
  uvputvrd_c(tno,"ut",time_var,1);
  uvputvrd_c(tno,"lst",time_var,1);
  nn=0;
  p1=0.0;
  p2=0.0;
  p3=0.0;
  uvset_c(tno,"corr","r",nn,p1,p2,p3);  /* write floats and not scaled integers */
  var_ivalue[0]=2;
  uvputvri_c(tno,"nants",var_ivalue,1);
  var_ivalue[0]=NUMBER_OF_CHANNELS;
  uvputvri_c(tno,"nchan",var_ivalue,1);
  var_ivalue[0]=1;
  uvputvri_c(tno,"npol",var_ivalue,1);
  var_ivalue[0]=1;
  uvputvri_c(tno,"nspect",var_ivalue,1);
  var_ivalue[0]=-1;
  uvputvri_c(tno,"pol",var_ivalue,1);
  var_ivalue[0]=NUMBER_OF_CHANNELS;
  uvputvri_c(tno,"nschan",var_ivalue,1);
  var_ivalue[0]=1;
  uvputvri_c(tno,"ischan",var_ivalue,1);
  var_ivalue[0]=1;
  uvputvri_c(tno,"ntpower",var_ivalue,1);
  var_ivalue[0]=0;
  uvputvri_c(tno,"nwide",var_ivalue,1);
  tpower[0]=1.0;
  uvputvrr_c(tno,"tpower",tpower,1);
  ant_az[0]=0.0;
  ant_az[1]=0.0;
  ant_el[0]=90.0;
  ant_el[1]=90.0;
  uvputvrd_c(tno,"antaz",ant_az,2);
  uvputvrd_c(tno,"antel",ant_el,2);
  var_veldop[0]=0.0;
  uvputvrr_c(tno,"veldop",var_veldop,1);
  var_vsource[0]=0.0;
  uvputvrr_c(tno,"vsource",var_vsource,1);
  var_restfreq[0]=0.0;
  uvputvrd_c(tno,"restfreq",var_restfreq,1);
  freq_channel1[0]=START_FREQUENCY; /* GHz */
  uvputvrd_c(tno,"sfreq",freq_channel1,1);
  freq_inc[0]=CHANNEL_WIDTH*1.00e-03; /* GHz */
  uvputvrd_c(tno,"sdf",freq_inc,1);
  site_latitude[0] = (double)(SITE_LATITUDE * PI/180.0);
  uvputvrd_c(tno,"latitud",site_latitude,1);
  site_longitude[0] = (double)(SITE_LONGITUDE * PI/180.0);
  uvputvrd_c(tno,"longitu",site_longitude,1);
  antenna1=1;
  antenna2=2;
  baseline_value[0]=(float)(256*antenna1+antenna2);
  uvputvrr_c(tno,"baseline",baseline_value,1);
  coord_var[0]=0.0;
  coord_var[1]=0.0;
  uvputvrd_c(tno,"coord",coord_var,2);
  uvputvr_c(tno,1,"veltype",velocity_type,8);

  int_time=INTEGRATION_TIME;
  printf(" Type number of %ld-sec spectral records to be written: ",int_time);
  gets(inbuf);
  sscanf(inbuf,"%ld",&nspect);

  printf("Type start UTC yyyy, mm (unit offset), dd, hh, mm, ss.ss :");
  gets(inbuf);
  sscanf(inbuf,"%ld%ld%ld%ld%ld%f",&yyyy,&mmmm,&dddd,&hrhr,&minmin,&secsec);
  jjdd = ( 1461 * ( yyyy + 4800 + ( mmmm - 14 ) / 12 ) ) / 4 +
    ( 367 * ( mmmm - 2 - 12 * ( ( mmmm - 14 ) / 12 ) ) ) / 12 -
    ( 3 * ( ( yyyy + 4900 + ( mmmm - 14 ) / 12 ) / 100 ) ) / 4 +
    dddd - 32075;
  julian_date = (double)jjdd -0.5;
  julian_date += (double)secsec/(3600.0*24.0) + 
    (double)minmin/(60.0*24.0) + (double)hrhr/24.0;

  julian_date -= (double)INTEGRATION_TIME/(3600.0*24.0);  
	// backoff timestamp by one integration
  	// time because time is incrememted first
  	// off in the loop below
  	
  /* Read the file containing recombination spectrum information */
  
  ffmin=1.0;
  ffmax=7.0;
  xarray=vector(1,NN);
  yarray=vector(1,NN);
  y2=vector(1,NN);
  npoints1=0;
  fptr1  = fopen(filename_rec,"r");
  while(fscanf(fptr1,"%g%g",&read1,&read2) != EOF)
    {
      if(read1 >= ffmin && read1 <= ffmax)
	{
	  ++npoints1;
	  xarray[npoints1]=read1;
	  yarray[npoints1]=read2;
	}
    }
  fclose(fptr1);
  printf ("Got the recombination line spectral values in range %f to %f GHz\n",ffmin,ffmax);

  /* calculation of second order differential for spline interpolation */

  yp1=(yarray[2]-yarray[1])/(xarray[2]-xarray[1]);
  ypn=(yarray[npoints1]-yarray[npoints1-1])/(xarray[npoints1]-xarray[npoints1-1]);
  spline(xarray, yarray, npoints1, yp1, ypn, y2);
 
// compute the coefficient arrays here for spectral interpolations

// fptr2 to store spectral indices that are outliers, if any
// fptr2  = fopen(filename7,"w+");

// printf(" Omit pixels with spectral index outside range: -2.0 to -3.0 \n");
// printf("  Listing index and coefficients for such points: \n");

  float sig[3];
  float aa,bb,sigaa,sigbb,chi2,qq;

  npolyfit=1;
  nbadpixels=0;
  sig[1]=1.0;sig[2]=1.0;sig[3]=1.0;
  for(j=0;j<nhpix;++j)
    {
      yy[1]=log10f(sky_408[j]);
      yy[2]=log10f(sky_1420[j]);
      yy[3]=log10(sky_23000[j]); 
      xx[1]=log10f(0.40800*1000.0);
      xx[2]=log10f(1.42000*1000.0);
      xx[3]=log10f(22.69000*1000.0);
//      polcof(xx,yy,npolyfit,cof);
	  fit(xx,yy,3,sig,0,&aa,&bb,&sigaa,&sigbb,&chi2,&qq);
	  cof[0]=aa;
	  cof[1]=bb;
//	  printf(" %d %f %f \n",j,cof[0],cof[1]); 
//	  if (cof[1] < -3.0) 
//		{ 
//		printf(" %d %f %f %f \n",j,cof[0],cof[1], cof[2]);
//		cof[1]=0.0; cof[0]=0.0; 
//		}
//	  if (cof[1] > -2.2 || cof[1] < -2.8) 
//		{ 
//		printf(" Warning: spectral index flatter than -2.2 or steeper than -2.8 %f\n",cof[1]);
//		nbadpixels += 1;
//		cof[1]=0.0; cof[0]=0.0; 
//		}

      cc[j][0]= cof[0];
      cc[j][1]= cof[1];

//     cc[j][2]= cof[2];

//     fprintf(fptr2,"%f ",cof[1]);
    }
//  fprintf(fptr2," \n");
//  fclose(fptr2);

	printf("Number of bad pixels = %d of %d \n",nbadpixels,nhpix);

// Open ascii file, 
// ftpr1 to store all the simulated spectra for python processing

  int fflush(FILE *fptr1);
  fptr1  = fopen(filename6,"a+");
  for(i=0;i<NUMBER_OF_CHANNELS;++i)
    {
      cfreq = (double)(START_FREQUENCY) + ((double)(CHANNEL_WIDTH)/1000.0)*(double)(i);
      fprintf(fptr1,"%f ",cfreq); // 01 Jul 2014
      /* compute the recombination spectrum intensity at the current frequency*/
      splint(xarray,yarray,y2,npoints1,cfreq,&yf);
      new_yf[i]=(double)yf;

      cfreq1_hz = cfreq*1e9;
      cmb_intensity[i] = (2.0*(double)hh*cfreq1_hz*cfreq1_hz*cfreq1_hz)/
	        (((double)cvel*(double)cvel)*
			(exp(((double)hh*(double)cfreq1_hz)/((double)kk*(double)TCMB))-1.0));
      P_cold[i]=2*kk*T_cold*(((hh*cfreq1_hz)/(kk*T_cold))/(exp((hh*cfreq1_hz)/(kk*T_cold))-1));
      P_hot[i]=2*kk*T_hot*(((hh*cfreq1_hz)/(kk*T_hot))/(exp((hh*cfreq1_hz)/(kk*T_hot))-1));
      P_diff[i]=P_hot[i]-P_cold[i];
    }
  
  fprintf(fptr1," \n");  

/* initialize seed for random number generator for adding gaussian noise */
   	time(&utime);
   	seed = -(long)utime;
	/* seed_count += 1; */
  	printf (" Start loop for writing nspect spectral records \n");

  for(ii=0;ii<nspect;++ii) 
    {
      // write the header record here for the current time

      source_name[0]='S';
      source_name[1]='P';
      source_name[2]='C';
      source_name[3]='0';
      source_name[4]='\0';
      uvputvr_c(tno,1,"source",source_name,4);

      preamble[0]=0.0;  /* u coordinate */
      preamble[1]=0.0;  /* v coordinate */

      // assemble the current JD in preamble[2] from UTC

      julian_date += (double)INTEGRATION_TIME/(3600.0*24.0);
      preamble[2] = julian_date;

      // Type UTC on terminal

      julian_day = (int)floor(julian_date + 0.5);
      time_utc = (julian_date+0.5) - (double)julian_day;
      if(time_utc >= 1.0) time_utc -= 1.0;
      if(time_utc < 0.0) time_utc += 1.0;
      time_utc *= 24.0;
      time_utc_hh = (int)floor(time_utc);
      time_utc_res = time_utc - (double)time_utc_hh;
      time_utc_mm = (int)floor(60.0*time_utc_res);
      time_utc_sec = (double)(60.0*(60.0*time_utc_res-time_utc_mm));
      printf(" UTC: %2d %2d %5.2f\n",time_utc_hh,time_utc_mm,time_utc_sec);

      antenna1=1;
      antenna2=2;
      preamble[3]=(double)(256*antenna1+antenna2);

      // compute LST and record this variable

      longitude_radians = (double)site_longitude[0];
      lst=cal_lst(julian_date, longitude_radians);
      time_var[0]=lst; 
      uvputvrd_c(tno,"lst",time_var,1);

      // prepare and write the spectral data

      for(i=0;i<NUMBER_OF_CHANNELS;++i)
	{ 
	  cross_real[i]=0.0;
	  cross_imag[i]=0.0;
	  sumwt[i]=0.0;
	}	

      for(j=0;j<nhpix;++j)
	{

	  c_ll=(double)ll_coordinate[j];
	  c_bb=(double)bb_coordinate[j];
	  cwt=0.0;
	  cof0=(double)cc[j][0];
	  cof1=(double)cc[j][1];
//	  cof2=(double)cc[j][2];
	  if(cof0==0.0 && cof1==0.0 && cof2==0.0) continue;
		
	  // Convert ll,bb to ra,dec B1950.0 epoch
	  // Wikipedia - Celestial coordinate systems

	  ra_b1950 = atan2( sin(c_ll - 123.0*PI/180.0),
			    (cos(c_ll - 123.0*PI/180.0)*sin(27.4*PI/180.0)-tan(c_bb)*cos(27.4*PI/180.0)) )
	    + 12.25*PI/180.0;
	  dec_b1950 = asin( sin(c_bb)*sin(27.4*PI/180.0) + 			
			    cos(c_bb)*cos(27.4*PI/180.0)*cos(c_ll-123.0*PI/180.0) );

	  // Precess from B1950.0 to date 

	  ra_precess = (180.0/PI)*(float)ra_b1950;
	  dec_precess = (180.0/PI)*(float)dec_b1950;
	  epoch1=1950.0;
	  epoch2=(float)yyyy + (float)(mmmm-1)/12.0 ;

	  precess_(&ra_precess,&dec_precess,&epoch1,&epoch2);

	  ra_date = (double)(ra_precess*PI/180.0);
	  dec_date = (double)(dec_precess*PI/180.0);

	  // Convert ra,dec date to az,el using the current LST
	  // AZ is defined as the angle from N towards E

	  haha = lst - ra_date; // Convert LST and RA to HA

	  azaz = (double)(PI) + atan2( sin(haha) , 
				       (cos(haha)*sin(site_latitude[0])-tan(dec_date)*cos(site_latitude[0])) );
	  elel = asin( sin(site_latitude[0])*sin(dec_date)
		       +cos(site_latitude[0])*cos(dec_date)*cos(haha) );

	  /* Included calculation to correct for atmospheric refraction 
	     elel is in radians, new_elel is also in radians, site_altitude is in meters */ 
	  new_elel = refraction(elel, site_altitude); 

	  // If az,el is within the beam weight by the beam pattern and accumulate the spectrum
	  
	  for(i=0;i<NUMBER_OF_CHANNELS;++i)
	    {  
	      cfreq = (double)(START_FREQUENCY) + ((double)(CHANNEL_WIDTH)/1000.0)*(double)i;  // GHz
	      cwt = beam_definition(cfreq,azaz,new_elel);

//		  if( cof1 > -2.2 || cof1 < -2.8) cwt=0.0;

	      if(cwt > 0.0) 
		{ 
		  cfreq1=cfreq;
		  cfreq1_hz = cfreq1*1e9;
		  cfreq = log10(cfreq*1000.0);
		  ctemp = cof0 + (cof1)*cfreq /* + cof2*cfreq*cfreq */ ;
		  ctemp = pow(10.0,ctemp);  // this is brightness temperature of foreground

		  /* Convert sky temperature to intensity, 
		     add the recombination spectrum intensity to the computed sky brightness 
		     and record the final antenna temperature */
		  intensity = ( 1.0*new_yf[i] + cmb_intensity[i]
		  	 + (2.0*(cfreq1_hz)*(cfreq1_hz)*kk*ctemp/(cvel*cvel))
		  	)*(cvel*cvel)/(cfreq1_hz*cfreq1_hz);
		  // uncomment below for no recombination lines in spectra
		  /* intensity = ( cmb_intensity[i] */
		  /* 	 + (2.0*(cfreq1_hz)*(cfreq1_hz)*kk*ctemp/(cvel*cvel)) */
		  /* 	)*(cvel*cvel)/(cfreq1_hz*cfreq1_hz); */

		  final_temp = (intensity/P_diff[i])*(T_hot-T_cold);
		  cross_real[i] += cwt*final_temp; 
		  sumwt[i] += cwt;		  
		}
	    }	
	}

      for(i=0;i<NUMBER_OF_CHANNELS;++i)
	{
	  if(sumwt[i]>0.0) cross_real[i] /= sumwt[i];
	}

    /* Add thermal noise to the spectrum */
    /* for(i=0;i<NUMBER_OF_CHANNELS;++i) */
/*     	{ */
/* // This is for a correlation receiver assuming zero temperature from 4th port */
/* //	  sigma=sqrt((double)(2.0*(cross_real[i]/2.0)*(cross_real[i]/2.0)+ */
/* //	  2*(cross_real[i]/2.0)*Trx + (Trx*Trx)) / (double)(2*CHANNEL_WIDTH*1.0e6*NOISE_INT_TIME)); */
/* //	  cross_real[i] += 2*(sigma*gasdev(&seed)); */
/* //	  cross_imag[i] += 2*(sigma*gasdev(&seed)); */

/* // This is for a total power radiometer */
/* 	  sigma = (double)(cross_real[i] + Trx + T_atm)/sqrt((double)(CHANNEL_WIDTH*1.0e6*NOISE_INT_TIME)); */
/* 	  cross_real[i] += (sigma*gasdev(&seed)); */
/* 	  cross_imag[i] += (sigma*gasdev(&seed)); */

/* 	} */

      for(i=0;i<NUMBER_OF_CHANNELS;++i)
	{
	  data[2*i]=(float)cross_real[i];
	  data[2*i+1]=(float)cross_imag[i];
	  fprintf(fptr1,"%20.15lf ",cross_real[i]);	  
	}
      fprintf(fptr1,"\n");
      uvwrite_c(tno,preamble,data,flags,npoints);

    } /* end of loop writing ispect records */

  fclose(fptr1);
  hisclose_c(tno);
  uvclose_c(tno);
  
  free_vector(ll_coordinate,0,NHPIX);
  free_vector(bb_coordinate,0,NHPIX);
  free_vector(sky_408,0,NHPIX);
  free_vector(sky_1420,0,NHPIX);
  free_vector(sky_23000,0,NHPIX);
  free_vector(a_coeff,0,NHPIX);
  free_vector(b_coeff,0,NHPIX);
  free_vector(c_coeff,0,NHPIX);
  free_vector(xx,0,2);
  free_vector(yy,0,2);
  free_vector(cof,0,2);
  free_matrix(cc,0,NHPIX,0,2);
  free_vector(xarray,1,NN);
  free_vector(yarray,1,NN);
  free_vector(y2,1,NN);
  return 0;
} /* end of main */
