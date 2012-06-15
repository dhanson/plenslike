#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

typedef struct {
  int     nbins;
  int     lmax;
  int    *bin_lmins;
  int    *bin_lmaxs;
  double *bin_vals;
  double *mat_sigma;
  double *mat_sigma_inv;
  double *clpp_fid;
  double *cltt_fid;
  double *vl_inv;
  double *al_inv;
  double *fl;
  double *bl;
} plenslike_dat_basic;

void calc_plenslike_basic_bins( plenslike_dat_basic *dat, double *clpp, double *bins ) {
  int i, l;
  double num, den;

  for (i=0; i<dat->nbins; i++) {
    num = 0; den=0;
    for (l=dat->bin_lmins[i]; l<=dat->bin_lmaxs[i]; l++) {
      num += clpp[l] * dat->clpp_fid[l] * dat->vl_inv[l];
      den += dat->clpp_fid[l] * dat->clpp_fid[l] * dat->vl_inv[l];
    }
    bins[i] = num/den;
  }
}

double calc_plenslike_basic( plenslike_dat_basic *dat, double *clpp ) {
  int i, j;
  double ret;
  double *bins = malloc( dat->nbins * sizeof(double) );

  calc_plenslike_basic_bins(dat, clpp, bins);

  ret = 0.0;
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      ret += (bins[i] - dat->bin_vals[i]) * (bins[j] - dat->bin_vals[j]) * dat->mat_sigma_inv[i*dat->nbins + j];
    }
  }

  return -0.5*ret;
}

plenslike_dat_basic *load_plenslike_dat_basic( char *tfname ) {
  plenslike_dat_basic *dat;
  char line[32767];
  FILE *tf;
  int i, j, err, tmp;
  double dmp;
  
  tf = fopen(tfname, "r");
  assert( tf != NULL );    

  dat = malloc( sizeof( plenslike_dat_basic ) );

  // bypass comments
  while (fgets(line, sizeof line, tf)) {
    if (*line == '#') {
      continue;
    } else {
      break;
    }
  }

  // read size header
  err = sscanf(line, "%d", &dat->nbins);
  err = fscanf(tf, "%d", &dat->lmax);

  // allocate memory
  dat->bin_lmins     = malloc( dat->nbins * sizeof(int) );
  dat->bin_lmaxs     = malloc( dat->nbins * sizeof(int) );
  dat->bin_vals      = malloc( dat->nbins * sizeof(double) );
  dat->mat_sigma     = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->mat_sigma_inv = malloc( dat->nbins * dat->nbins * sizeof(double) );
  dat->clpp_fid      = malloc( (dat->lmax+1) * sizeof(double) );
  dat->cltt_fid      = malloc( (dat->lmax+1) * sizeof(double) );
  dat->vl_inv        = malloc( (dat->lmax+1) * sizeof(double) );
  dat->al_inv        = malloc( (dat->lmax+1) * sizeof(double) );
  dat->fl            = malloc( (dat->lmax+1) * sizeof(double) );
  dat->bl            = malloc( (dat->lmax+1) * sizeof(double) );

  // read bin info
  for (i=0; i<dat->nbins; i++) {
    err = fscanf(tf, "%d %d %d %lf", &tmp, &dat->bin_lmins[i], &dat->bin_lmaxs[i], &dat->bin_vals[i]);
    assert( tmp == i );
    assert( err == 4 );
  }

  // read sigma matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma[i*dat->nbins + j]);
      assert( err == 1 );
    }
  }

  // read sigma inv matrix
  for (i=0; i<dat->nbins; i++) {
    for (j=0; j<dat->nbins; j++) {
      err = fscanf(tf, "%lf", &dat->mat_sigma_inv[i*dat->nbins + j]);
      assert( err == 1 );
    }
  }

  // read spectra
  for (i=0; i<dat->lmax+1; i++) {
    err = fscanf(tf, "%lf %lf %lf %lf %lf %lf %lf", &dmp, &dat->clpp_fid[i], &dat->cltt_fid[i], &dat->vl_inv[i], &dat->al_inv[i], &dat->fl[i], &dat->bl[i]);
    assert( err == 7 );
    assert( dmp == i );
  }
  
  fclose(tf);
  return dat;
};

void free_plenslike_dat_basic( plenslike_dat_basic *dat ) {
  free(dat->bin_lmins);
  free(dat->bin_lmaxs);
  free(dat->bin_vals);
  free(dat->mat_sigma);
  free(dat->mat_sigma_inv);
  free(dat->clpp_fid);
  free(dat->cltt_fid);
  free(dat->vl_inv);
  free(dat->al_inv);
  free(dat->fl);
  free(dat->bl);

  free(dat);
};
