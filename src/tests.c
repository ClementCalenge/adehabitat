#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void topoids(double *vec, int *n);
void multpoco(double **tab, double *poco);
void aleadistrivec(double *vec, double *no);
void randksel(int *fac, double *pu, int *nani, int *ni);
void rks(int *fac, double *pdsu, int *nani, int *nbani, int *nl);
void ksel(double *tab, int *fac, double *poidsut, int *nhab, 
	  int *nani, int *nloctot, double *ut, double *di,
	  double *marg, int *nombreani, double *eigenvp, 
	  double *poidsco, int *ewa);
void permutksel(double *tab, int *fac, double *poidsut, int *nhab,
		int *nani, int *nloctot, double *ut, double *di,
		double *marg, int *nombreani, int *npermut, 
		double *obseig, double *simeig, double *obsmarg,
		double *simmarg, double *eigenvp, 
		double *simtout, double *poidsco,
		int *ewa);
void sahr2ksel(double *Usa, double *Uhr,  double *Ulo, int *nhab,
	       int *npix, int *nani, int *nlig, double *dud, 
	       int *fac, double *pu);
void nls2k(double *Usa, double *Uhr, int *nhab, 
	   int *npix, int *nani);
void rotxy(double *x, double *y, int k);
void shifthr(double **dispo, double **util, int *idl, int *idc);
void shr(double **carte, double **ze);
void sr(double *carter, double *zer, int *nlgr, int *ncgr);
void locrast(double *xgr, double *ygr, double *x, double *y,
	     double **carte);
void lr(double *xgri, double *ygri, double *xr, double *yr,
	double *carter, int *nco, int *nli, int *nlixy);
void getcarte(double **carte, double **kasc, int *indicecarte);
void gc(double *carter, double *kascr, int *nlgr, int *ncgr, int *nhab);
void comptePasNA(double **tab, int *nombre);
void videNA(double **entree, double **sortie, int *idcons);
void niche(double **X, double **Y, double *eig, double **mar);
void mvtfreeman(int *in, int *jn, int *dir, int *np);
void getcontour(double *grille, int *nlig, int *ncol, int *indicelig, 
		int *indicecol, int *lcont);
void lcontour(double *grille, int *nlig, int *ncol, int *lcont);
void levels(double *vec, double *lev, int *lvec);
void seqeticorr(double *grille, int *nlig, int *ncol);
void epa(double *X, double *Y, double *xl, double *yl, 
	 double *val, double *fen);
void kernelhr(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo);
void CVmise(int *nloc, double *xlo, double *ylo,
	    double *hvec, double *CV, int *nhteste);
void calcvolume(double *grille, int *ncolgri, int *nliggri, double *cellsize);
void calcsim(double *pix, double **pts, double *rg, 
	     int *nvar, int *npts, double *similarite);
void fctdomain(double *kascmod, double *ptsmod, 
	       double *range, int *npts, int *npix,  
	       int *nvar, double *qualhab);
void wml(double **used, double **avail, double *wmla, int na, int nh,
	 double **proj1, double **proj2, double *nbassocie, int krep);
void aclambda(double *util, double *dispo, int *nani, int *nhab,  
	      double *xxtxxtmod1, double *xxtxxtmod2, double *rnv,
	      double *wmla, int *nrep, double *wm, double *nb);
void rankma(double *used, double *avail, double *rankmap, double *rankmam,
	    double *rankmav, double *rankmanb, 
	    int *nhab, int *nani, int *nrep, double *rnv);
void erodil(double *grille, int *nlig, int *ncol, int *ntour, int *oper);
void inout(double *x, double *y, double *xp, double *yp,
	   int *deds);
void inoutr(double *xr, double *yr, double *xpr, double *ypr,
	    int *dedsr, int *nxr, int *npr);
void rastpol(double *xp, double *yp, double *xg, double *yg,
	     double **carte);
void rastpolaire(double *xpr, double *ypr, double *xgr, double *ygr,
		 double *carter, int *nlg, int *ncg, int *nvp);
void calcniche(double **kasc, int *nvar, int *nlg, int *ncg,
	       double *margvar, double *tolvar, double **carte);
void calcnicher(double *kascr, int *nvar, int *nlg, int *ncg,
		double *margvar, double *tolvar, double *carter);
void randompol(double *xpr, double *ypr, double *kascr,
	       double *marg, double *tol, int *nvar,
	       double *xgr, double *ygr, int *nlr, 
	       int *ncr, int *nvpr, int *nrep);
void dedans(double *pts, double *xc, double *yc, double *na,
	    double cs, double **asc);
void dedansr(double *ptsr, double *xcr, double *ycr, double *na,
	     double *cs, double *ascr, int *nl, int *nc, int *nlocs);
void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r);
void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr);
void joinkasc(double **xp, double **kasc, double **res, int nl, int nc,
	      double *xc, double *yc, double *cs);
void joinkascr(double *xpr, double *kascr, int *nlasc, int *ncasc,
	       double *xcr, double *ycr, double *cs, int *nlocs,
	       int *nvar, double *resr);
void randmargtol(double *xpr, double *rcrx, double *rcry, double *ascr, 
		 double *cwr, double *kascr, double *xcr, double *ycr,
		 double *csr,
		 double *tabdistr, double *dtr, double *anglesr, double *marr, 
		 double *tolr, int *nrepr, int *nlasc, 
		 int *ncasc, int *nvarr, int *nltdr, int *nlocsr);
void rpoint(double **xp, double *rcx, double *rcy, double **asc, 
	    double *xc, double *yc, double *cs);
void randmargtolpts(double *xpr, double *rcrx, double *rcry, double *ascr, 
		    double *cwr, double *kascr, double *xcr, double *ycr,
		    double *csr, double *marr, double *tolr, 
		    int *nrepr, int *nlasc, 
		    int *ncasc, int *nvarr, int *nlocsr);
void regroufacasc(double **asce, double **ascs, int *np,
		  int *nlev);
void regroufacascr(double *ascer, double *ascsr, int *npr,
		   int *nlevr, int *nle, int *nce, int *nls, 
		   int *ncs);
void regrouascnum(double **ascent, double **ascso);
void regrouascnumr(double *ascentr, double *ascsor, 
		   double *nler, double *ncer,
		   double *nlsr, double *ncsr);
void regroukasc(double *kascr, double *kascniou, int *nrow, 
		int *ncol, int *nvar, int *npix,
		int *typer, int *nrniou, int *ncniou);
void matmudemi(double **X, double **Y);
void matmudemir(double *Xr, double *Yr, int *ncr);
void enfa(double **Z, double *p, int *nvar, int *npix,
	  double *vp);
void enfar(double *Zr, double *pr, int *nvar, int *npix,
	   double *vpr);
void randenfa(double **Z, double *p, int *nrep, double *res);
void randenfar(double *Zr, double *pr, int *nvar, int *npix,
	       int *nrep, double *resr);
void integrno(double *XG, double *X1, double *X2, 
	      double *T, double *sig1,
	      double *sig2, double *alpha, double *res);
void kernelbb(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *sig1, double *sig2, 
	      double *xlo, double *ylo, double *Tr, int *controlbox, 
	      int *nalpha);
void ligpoly(double *x, double *y, double r, double *xp, double *yp);
void buflig(double **x, double r, double **carte, double *xg, double *yg);
void bufligr(double *xr, double *rr, double *carter, 
	     double *xgr, double *ygr, int *nlr, int *ncr, 
	     int *nlocr);
void distxy(double **xy1, double **xy2, double *di);
void distxyr(double *xy1r, double *xy2r, int *n1r, 
	     int *n2r, double *dire);
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di);
void fptt(double *x, double *y, double *t, int pos, double radius, 
	  double *fptto, int nlo);
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt);
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr);
void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm);
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel);
void resolpol(double a, double b, double c, double *x1, double *x2, 
	      int *warn);
void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff);
void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, 
		 int *neff);
void trouveclustmin(double **xy, int *clust, int *lo1, int *lo2,
		    int *lo3, double *dist);
void trouveclustminr(double *xyr, int *nr, int *clustr, int *lo1, int *lo2,
		     int *lo3, double *dist);
void nndistclust(double **xy, double *xyp, double *dist);
void parclust(double **xy, int *clust, int *noclust, 
	      int *noloc, double *dist);
void trouveminclust(double **xy, int *liclust, int *clust, 
		    int *noclust, int *noloc, double *dist);
void choisnvclust(double **xy, int *liclust, int *clust, int *ordre);
void clusterhr(double **xy, int *facso, int *nolocso, int *cluso);
void longfacclust(double **xy, int *len2);
void longfacclustr(double *xyr, int *nr, int *len2);
void clusterhrr(double *xyr, int *nr, int *facsor, 
		int *nolocsor, int *clusor, int *len);
void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr);
void runsltr(int *xr, int *nr, double *res, int *nrepr);
void testindepangl (double *sim, double *ang, int *nang, int *debut, 
                    int *fin, int *ndeb, int *ni);
void testindepdist (double *sim, double *di, int *ndi, int *debut, 
		    int *fin, int *ndeb, int *ni);
void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur);
void optcut (double **Pid, double **mk, int *maxk);
void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye);
void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res);
void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr);
void kcprcirc(double **xyd, double *h, double *x, double t, 
	      double *val);
void kcprlin(double **xyd, double *h, double *x, double t, 
	     double *val);
void kernelkcr(double *xydr, double *tcalcr, int *nlr, double *gridr,
	       double *xgri, double *ygri, int *nliggri, int *ncolgri, 
	       double *hr, int *circularr);
void findmaxgrid(double *grille, int *nlig, int *ncol);
void engen2008Ir(double *avr, double *usr, int *nliga, int *nligu, 
		 int *ncol, double *resr, int *nsimrar);
void engen2008r(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, int *idr, int *nidr, int *nsimr, 
		double *resr, int *nsimrar);
void free_ivector(int *v, int nl, int nh);
int invers(double **a, int n, double **b, int m);






/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    GetRNGstate();
    w = unif_rand();
    PutRNGstate();
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    for (i=1;i<=n;i++) {
	GetRNGstate();
	alea[i] = (1e8)*unif_rand();
	PutRNGstate();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.00000000000001) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutÚes par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
	for (i=0;i<=l1;i++) {
	    if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
		return;
		for (j=0;j<i;j++) {
		    free(*(*tab+j));
		}
	    }
	}
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	free((char *) *(tab+i) );
    }
    free((char *) tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    free((char *) vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    free((char *) vec);
    
}














/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitat               *****
 *********               -------------------------               *****
 *********************************************************************
 *********************************************************************
 */



/* ------------------------------------------------------ */
/* Converts a vector so that the sum of the elements of the 
   vector is equal to one */

void topoids(double *vec, int *n)
{
    /* Declaration of local variables */
    int i;
    double somme;
    somme = 0;
    
    /* the sum of elements of the vector */
    for (i=0; i<=(*n-1); i++){
	somme = somme + vec[i];
    }
    
    /* transformation to weights */
    for (i=0; i<=(*n-1); i++){
	vec[i] = vec[i] / somme;
    }
    
    *n=somme;
}




/* ------------------------------------------------------ */
/* Multiply a table by the square root of columns weights */

void multpoco(double **tab, double *poco)
{

    /* Declaration of local variables */
    int nc, nl, i, j;
    double k;
    
    nl = tab[0][0];
    nc = tab[1][0];
    
    /* Main operations */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    k = poco[j];
	    tab[i][j] = tab[i][j]*sqrt(k);
	}
    }
    
}


/* ------------------------------------------------------ */
/* Distributes "no" points randomly in a vector with p elements */

void aleadistrivec(double *vec, double *no)
{
    
/* Declaration of local variables */
    double tmp, i, j, n, lv;
    
    n = *no;
    lv = vec[0];
    
    /* The random distribution of points in the elements of the vector */
    for (i=1; i<=n; i++) {
	tmp = alea();
	for (j=1; j<=lv; j++) {
	    if ((tmp >= (j-1)/lv)&&(tmp < j/lv))
		vec[(int) j]++;
	}
    }
}





/* *****************************************************
   
Randomly distributes the points within the home range 
of animals (used for K-select analysis)

***************************************************** */

void randksel(int *fac, double *pu, int *nani, int *ni)
{
    /* Declaration of local variables */
    int i, j, k, l;
    double su, *tmp;
    
    l=1;
    
    /* For each animal, randomizes the points */
    for (k=1; k<=*nani; k++) {
	
	vecalloc(&tmp, ni[k]);
	su = 0;
	
	for (i=1; i<=ni[k];i++) {
	    su = su+pu[l];
	    l++;
	}
	
	aleadistrivec(tmp, &su);
	
	j=1;
	for (i=(l-ni[k]); i<l; i++) {
	    pu[i]=tmp[j];
	    j++;
	}
	
	freevec(tmp);
    }
    
}



/* *****************************************************

        The version to test randksel with R

   ***************************************************** */

void rks(int *fac, double *pdsu, int *nani, int *nbani, int *nl)
{

    /* Declaration and copy of R objects in local variables */
    
    int i, *fa, *ni;
    double *pu;
    
    vecalloc(&pu, *nl);
    vecintalloc(&fa, *nl);
    vecintalloc(&ni, *nani);
    
    for (i=1; i<=*nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=*nl; i++) {
	pu[i] = pdsu[i-1];
    }
    
    for (i=1; i<=*nani; i++) {
	ni[i] = nbani[i-1];
    }
    
    /* The function randksel */
    randksel(fa, pu, nani, ni);
    
    /* Output */
    for (i=1; i<=*nl; i++) {
	pdsu[i-1] = pu[i];
    }
    
    /* Free memory */
    freevec(pu);
    freeintvec(fa);
    freeintvec(ni);
}



/* ****************************************************************
   *                                                              *
   *              K-select Analysis                               *
   *                                                              *
   **************************************************************** */


void ksel(double *tab, int *fac, double *poidsut, int *nhab, 
	  int *nani, int *nloctot, double *ut, double *di,
	  double *marg, int *nombreani, double *eigenvp, 
	  double *poidsco, int *ewa)
{
    
    /* Declaration of local variables */
    int i,j,k, sommeloctot;
    double **ta, *pu, **use, **ava, **mar, *poidsani;
    double **inertie, *valpro, *spu, *poco;
    int nh, na, nl, rang;
    int *fa, *ni;
    
    
    /* Memory Allocation for local variables */
    nl = *nloctot; /* total number of pixels */
    na = *nani; /* number of animals */
    nh = *nhab; /* number of variables */
    
    vecintalloc (&fa, nl); /* factor with one level per animal */
    vecalloc (&pu, nl); /* utilization weight of the pixels */
    vecalloc (&poidsani, na); /* weight of each animal (proportional
				 to the nb of locs) */
    vecalloc (&valpro, nh); /* eigenvalues of the analysis */
    vecalloc (&spu, na); /* utilization weights (sum equal to 1) */
    vecalloc(&poco, nh); /* column weights in the analysis (habitat weights) */
    taballoc (&ta, nl, nh); /* initial table */
    taballoc (&use, na, nh); /* table of used means */
    taballoc (&ava, na, nh); /* table of available means */
    taballoc (&mar, na, nh); /* marginality table */
    taballoc (&inertie, nh, nh); /* Inertia matrix */
    vecintalloc(&ni, na); /* number of available pixels per animal */
    

    
    /* Copy of the R objects into the local C variables */
    k = 0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nh; j++) {
	    ta[i][j] = tab[k];
	    k = k + 1;
	}
    }
    
    for (i=1; i<=nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=nl; i++) {
	pu[i] = poidsut[i-1];
    }
    
    for (i=1; i<=nh; i++) {
	poco[i] = poidsco[i-1];
    }
    
    for (i=1; i<=na; i++) {
	spu[i] = 0;
    }
    
    for (i=1; i<=na; i++) {
	ni[i] = nombreani[i-1];
    }
    
    
    /* Computes the number of relocs per animal */
    for (i=1; i<=na; i++) {
	for (k=1; k<=nl; k++) {
	    if (fa[k]==i) {
		spu[i] = spu[i] + pu[k];	
	    }
	}
    }
    
    
    /* Computes the total number of relocations */
    sommeloctot=0;
    for (i=1; i<=na; i++) {
	sommeloctot = sommeloctot + spu[i];	
    }
    
    
    /* Computes the used mean and available mean */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    for (k=1; k<=nl; k++) {
		if (fa[k]==i) {
		    ava[i][j]= ava[i][j] + ta[k][j]/ni[i];
		    use[i][j]=use[i][j] + (ta[k][j] * pu[k]);
		}
	    }
	}
    }
    
    /* and the marginality */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    use[i][j] = use[i][j] / spu[i];
	    mar[i][j] = use[i][j] - ava[i][j];
	}
    }
    
    /* Output toward the R object */
    k = 0;
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    ut[k] = use[i][j];
	    di[k] = ava[i][j];
	    marg[k] = mar[i][j];
	    k = k + 1;
	}
    }
    
    /* column weighting */
    multpoco(mar, poco);
    
    
    /* Computation of the weight of each animal */
    for (i=1; i<=na; i++) {
	if (*ewa==0)
	    poidsani[i] = (double) spu[i] / sommeloctot;
	if(*ewa==1)
	    poidsani[i] = (double) 1 / na;
    }
    
    sqrvec(poidsani);
    
    
    /* Computation of the eigenstructure */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    mar[i][j] = mar[i][j] * poidsani[i];
	}
    }
    prodmatAtAB(mar, inertie);
    DiagobgComp(nh, inertie, valpro, &rang);
    
    
    /* Output toward R */
    for (i = 1; i<=rang; i++) {
	eigenvp[i-1] = valpro[i];
    }
    
    
    /* Free memory */
    freeintvec (fa);
    freeintvec (ni);
    freevec (pu);
    freevec (poidsani);
    freevec(valpro);
    freevec (spu);
    freevec(poco);
    freetab (ta);
    freetab (use);
    freetab (ava);
    freetab (mar);
    freetab (inertie);
}




/* ****************************************************************
   *                                                              *
   *         Randomisation test related to the K-select           *
   *                                                              *
   **************************************************************** */


void permutksel(double *tab, int *fac, double *poidsut, int *nhab,
		int *nani, int *nloctot, double *ut, double *di,
		double *marg, int *nombreani, int *npermut, 
		double *obseig, double *simeig, double *obsmarg,
		double *simmarg, double *eigenvp, double *simtout, 
		double *poidsco, int *ewa)
{
    /* Declaration of local variables */
    double **ta, *pu, *obstout;
    int na, nh, nl, i, j, k, l, m, q, *ni, *fa, *numero, nbperm;
    
    /* Memory Allocation */
    na=*nani;
    nh=*nhab;
    nl=*nloctot;
    nbperm=*npermut;
    
    taballoc(&ta, nl, nh);
    vecintalloc(&fa, nl);
    vecintalloc(&numero, nbperm);
    vecalloc(&pu, nl);
    vecalloc(&obstout, na*nh);
    vecintalloc(&ni, na);
    
    
    /* One copies R objects toward C objects (local variables) */
    k = 0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nh; j++) {
	    ta[i][j] = tab[k];
	    k = k + 1;
	}
    }
    
    for (i=1; i<=nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=nl; i++) {
	pu[i] = poidsut[i-1];
    }
    
    for (i=1; i<=na; i++) {
	ni[i] = nombreani[i-1];
    }
    
    /* --------------------------------------- */
    /* Main Computations */
    
    /* Observed kselect */
    ksel(tab,  fac, poidsut, nhab, 
	 nani, nloctot, ut, di,
	 marg, nombreani, eigenvp, poidsco, ewa);
    
    
    /* observed values are placed in the output */
    k=0;
    j=0;
    *obseig = eigenvp[0];
    for (i=0; i<(na*nh); i++) {
	obsmarg[j] = obsmarg[j] + (marg[i] * marg[i] * poidsco[k]);
	if (k==(nh-1)) {
	    k=-1;
	    j++;
	}
	k++;
    }
    
  
    /* Coordinates of the marginality vectors */
    for (i=1; i<=na*nh; i++) {
	obstout[i] = marg[i-1];
    }
    
    
    
    /* The permutations */
    m=0; /* will indicate the step for simmarg */
    q=0; /* will indicate the step for simtout */
    
    for (k=1; k<=nbperm; k++) {
	
	/* One permutes */
	randksel(fa, pu, nani, ni);
	
	/* One copies the randomized utilization weights in poidsut */
	for (i=0; i<nl; i++) {
	    poidsut[i]=pu[i+1];
	}
	
	/* and therefore, go */
	ksel(tab,  fac, poidsut, nhab, 
	     nani, nloctot, ut, di,
	     marg, nombreani, eigenvp, poidsco, ewa);
	
	
	/* Simulated values in the output */
	simeig[k-1] = eigenvp[0];
	l=0;
	
	for (i=0; i<(na*nh); i++) {
	    simmarg[m] = simmarg[m] + (marg[i] * marg[i] * poidsco[l]);
	    if (l==(nh-1)) {
		l=-1;
		m++;
	    }
	    l++;
	}
	
	for (i=0; i<(na*nh); i++) {
	    simtout[q] = marg[i];
	    q++;
	}
    }
    
    /* obstout is put again in marg */
    for (i=0; i<(na*nh); i++) {
	marg[i] = obstout[i+1];
    }
    
    
    /* Memory free */
    freetab(ta);
    freeintvec(fa);
    freevec(pu);
    freevec(obstout);
    freeintvec(ni);
    freeintvec(numero);
}








/* ****************************************************************
   *                                                              *
   *         converts a sahrlocs object for K-select              *
   *                                                              *
   **************************************************************** */


void sahr2ksel(double *Usa, double *Uhr,  double *Ulo, int *nhab,
	       int *npix, int *nani, int *nlig, double *dud, 
	       int *fac, double *pu)
{

  /* declaration of local variables */
  int i,j,k,l,na,nh,np, nl;
  double **SA, **HR, **LOCS, **sortie;
  int *idna;
  /* idna will content 1 for pixels non-NA on the SA */

  
  /* Allocation of Memory */
  na = *nani;
  nh = *nhab;
  np = *npix;
  nl = *nlig;
  
  taballoc(&SA, np, nh);
  taballoc(&HR, np, na);
  taballoc(&LOCS, np, na);
  taballoc(&sortie, nl, nh);
  vecintalloc(&idna, np);  

  
  /* Copy R objects into local C variables */
  k = 0;
  for (i=1; i<=np; i++) {
      for (j=1; j<=nh; j++) {
	  SA[i][j] = Usa[k];
	  k = k + 1;
      }
  }
  
  k = 0;
  for (i=1; i<=np; i++) {
      for (j=1; j<=na; j++) {
	  HR[i][j] = Uhr[k];
	  LOCS[i][j] = Ulo[k];
	  k = k + 1;
      }
  }
  
  
  /* Computations */
  /* Computation of the row number in the output table */
  for (i=1; i<=np; i++) {
    if (fabs(SA[i][1] + 9999) > 0.000000001) {
      idna[i] = 1; /* = 1 if non NA */
    }
  }
  
  /* And finally, the output table */
  l=0;
  
  for (i=1; i<=np; i++) {
      for (j=1; j<=na; j++) {
	  if ((idna[i]==1)&&(((int) HR[i][j])==1)) {
	      l++;
	      for (k=1; k<=nh; k++) {
		  sortie[l][k] = SA[i][k]; /* We will pass it to dud after */
		  fac[l-1] = j;
		  pu [l-1] = LOCS[i][j];
	      }
	  }
      }
  }
  
  /* One returns the output to R */
  k = 0;
  for (i=1; i<=nl; i++) {
      for (j=1; j<=nh; j++) {
	  dud[k] = sortie[i][j];
	  k++;
      }
  }
  
  
  
  /* Free memory */
  freetab(SA);
  freetab(HR);
  freetab(LOCS);
  freetab(sortie);
  freeintvec(idna);

}





/* ****************************************************************
   *                                                              *
   *    Computation of the number of rows of the output table     *
   *                                                              *
   **************************************************************** */


void nls2k(double *Usa, double *Uhr, int *nhab, 
	   int *npix, int *nani)
{
  /* Declaration of local variables */
  int i,j,k,na,nh,np,nl;
  double **SA, **HR;
  int *ni; /* number of pixels for each animal */
  int *idna;
  /* idna will content 1 for the pixels non-NA on the SA */
  
  
  /* Allocation of memory */
  na = *nani;
  nh = *nhab;
  np = *npix;
  
  taballoc(&SA, np, nh);
  taballoc(&HR, np, na);
  vecintalloc(&ni, na);
  vecintalloc(&idna, np);  
  
  /* Copies R objects into the local C variables */
  k = 0;
  for (i=1; i<=np; i++) {
      for (j=1; j<=nh; j++) {
	  SA[i][j] = Usa[k];
	  k = k + 1;
      }
      idna[i] = 1;
  }
  
  k=0;
  for (i=1; i<=np; i++) {
      for (j=1; j<=na; j++) {
	  HR[i][j] = Uhr[k];
	  k = k + 1;
      }
  }
  
  
  /* Computation */
  /* Number of rows of the output table */
  for (i=1; i<=np; i++) {
    if (fabs(SA[i][1] + 9999) < 0.000000001) {
      idna[i] = 0; /* = 1 if non NA */
    }
  }
  
  /* Computation of ni */
  for (j=1; j<=na; j++) {
    ni[j] = 0;
  }

  for (i=1; i<=np; i++) {
    if (idna[i] == 1) {
      for (j=1; j<=na; j++) {
	if (((int) HR[i][j]) == 1) {
	  ni[j] = ni[j]+1;
	}
      }
    }
  }
  
  
  /* Total number of lines in the output table */
  nl=0;
  for (i=1; i<=na; i++) {
    nl = nl + ni[i];
  }
  *nani = nl;
  
  
  /* Free Memory */
  freetab(SA);
  freetab(HR);
  freeintvec(ni);
  freeintvec(idna);
}






/* ****************************************************************
   *                                                              *
   * rotxy to rotate randomly a un couple (x,y)                   *
   *                                                              *
   **************************************************************** */

void rotxy(double *x, double *y, int k)
{
    /* Declaration of local variables */
    int i, n, *numero;
    double mx, my, *angle, *angleb, ang, co, si, xt, yt;
    
    /* Computation of the mean */
    mx=0;
    my=0;
    n=x[0];
    
    vecalloc(&angle, 360);
    vecalloc(&angleb, 360);
    vecintalloc(&numero, 360);
    
    for (i=1; i<=n; i++) {
	mx = mx + x[i];
	my = my + y[i];
    }
    
    mx = mx / n;
    my = my / n;
    
    /* Centring */
    for (i=1; i<=n; i++) {
	x[i] = x[i]-mx;
	y[i] = y[i]-my;
    }
    
    /* Draws a random number between 0 and 2pi */
    for (i=1; i<=360; i++) {
	angle[i] = (((double) i)*3.14159265359)/180;
    }
    
    /* GO */
    getpermutation(numero, k);
    vecpermut(angle, numero, angleb);
    ang = angleb[1];
    co = cos(ang);
    si = sin(ang);
    
    for (i=1; i<=n; i++) {
	xt = x[i];
	yt = y[i];
	
	x[i]= co * xt - si * yt + mx;
	y[i]= si * xt + co * yt + my;
    }
    
    /* Free the memory */
    freevec(angle);
    freevec(angleb);
    freeintvec(numero);
}


/* **************************************************************
 *                                                              *
 * shifthr randomly shifts a home range on an area              *
 *                                                              *
 **************************************************************** */

void shifthr(double **dispo, double **util, int *idl, int *idc)
{
    /* Declaration of local variables */
    int i, j, l, ncgr, nlgr, ncpe, nlpe;
    int *idlgr, *idcgr, crand, lrand;
    
    /* Memory Allocation */
    nlgr = dispo[0][0];
    ncgr = dispo[1][0];
    nlpe = util[0][0];
    ncpe = util[1][0];
    
    vecintalloc(&idcgr, ncgr-ncpe+1);
    vecintalloc(&idlgr, nlgr-nlpe+1);
    
    /* ************** Random drawing of the HR position *********** 
     Two conditions:
     1. The square describing the HR is included in the study area
     2. No 0 where relocations
    */
    l=0;
    
    while (l==0) {
	getpermutation(idcgr, *idc); /* Random draw column */
	getpermutation(idlgr, *idl); /* Random draw row */
	crand = idcgr[1];
	lrand = idlgr[1];
	
	l=1;
	for (i=1; i<=nlpe; i++) {
	    for (j=1; j<=ncpe; j++) {
		if (util[i][j]>0) {
		    if (fabs(dispo[i+lrand-1][j+crand-1] + 9999) < 0.000000001) {
			l=0;
		    }
		}
	    }
	}
    }
    
    *idl=lrand;
    *idc=crand;
    
    /* Memory free */
    freeintvec(idcgr);
    freeintvec(idlgr);
    
}




/* ****************************************************************
   *                                                              *
   * shr places the Homr range randomly on the study area         *
   *                                                              *
   **************************************************************** */

void shr(double **carte, double **ze)
{
    /* Declaration of local variables */
    int i, j, l, m, idci, idli, idls, idcs;
    int *ligne, *colonne, nlsc, ncsc, nlg, ncg;
    double **souscar;
    
    idli = 0;
    idci = 0;
    idls = 0;
    idcs = 0;
    nlg = carte[0][0];
    ncg = carte[1][0];
    
    vecintalloc(&ligne, nlg);
    vecintalloc(&colonne, ncg);
    
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    ligne[i] = ligne[i] + carte[i][j];
	    colonne[j] = colonne[j] + carte[i][j];
	}
    }
    
    /* One computes the "min" and "max" indices of the 
       rows and columns containing the rasterized locs */
    for (i=1; i<=nlg; i++) {
	if ((idli==0)&&(ligne[i]!=0)) idli = i;
    }
    for (i=nlg; i>=1; i--) {
	if ((idls==0)&&(ligne[i]!=0)) idls = i;
    }
    for (i=1; i<=ncg; i++) {
	if ((idci==0)&&(colonne[i]!=0)) idci = i;
    }
    for (i=ncg; i>=1; i--) {
	if ((idcs==0)&&(colonne[i]!=0)) idcs = i;
    }
    
    /* Finally, computation of the number of rows and columns of souscar */
    nlsc = idls - idli + 1;
    ncsc = idcs - idci + 1;
    
    /* Memory Allocation for souscar */
    taballoc(&souscar, nlsc, ncsc);
    
    /* attributes values to cells of souscar */
    l = 1;
    m = 1;
    for (i=idli; i<=idls; i++) {
	for (j=idci; j<=idcs; j++) {
	    souscar[l][m] = carte[i][j];
	    m++;
	}
	m = 1;
	l++;
    }
    
    /* Randomisation of the position of the locs on the SA */
    shifthr(ze, souscar, &idli, &idci);
    /* idli et idci coontain resp. row and column indices for 
       the randomized map 
       
       As "carte" already contain the locs with a randomized orientation, we
       use ze to store the randomized position of the locs */
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    ze[i][j] = 0;
	}
    }
    
    /* We recompute the complete randomized map */
    l = 1;
    m = 1;
    for (i=idli; i<=(idli+nlsc-1); i++) {
	for (j=idci; j<=(idci+ncsc-1); j++) {
	    ze[i][j] = souscar[l][m];
	    m++;
	}
	m = 1;
	l++;
    }
    
    
    /* Free memory */
    freetab(souscar);
    freeintvec(ligne);
    freeintvec(colonne);
}




/* ****************************************************************
   *                                                              *
   *         sr = interactive version of shr with R               *
   *                                                              *
   **************************************************************** */

void sr(double *carter, double *zer, int *nlgr, int *ncgr)
{
    /* Declaration of local variables and memory allocation */
    double **carte, **ze;
    int i,j,k,nlg,ncg;
    nlg = *nlgr;
    ncg = *ncgr;
    taballoc(&carte,nlg,ncg);
    taballoc(&ze,nlg,ncg);
    
    /* Copies the values from the R objects to the C objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    ze[i][j] = zer[k] ;
	    carte[i][j] = carter[k];
	    k++;
	}
    }
    
    /* Use of shr */
    shr(carte, ze);
    
    /* Copies the result from C objects to R objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    zer[k] = ze[i][j];
	    carter[k] = carte[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(carte);
    freetab(ze);
    
}






/* ****************************************************************
   *                                                              *
   *   locrast allows the rasterization of the relocations        *
   *                                                              *
   **************************************************************** */

void locrast(double *xgr, double *ygr, double *x, double *y,
	     double **carte)
{
    /* Declaration of local variables */
    int i, j, k, n, nc, nl;
    double res;
    
    /* Memory allocation */
    res = xgr[2]-xgr[1];
    n = x[0];
    nl = carte[0][0];
    nc = carte[1][0];
    
    /* Map */
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    carte[i][j] = 0;
	}
    }
    
    
    /* rasterisation of the locs */
    for (k=1; k<=n; k++) {
	for (i=1; i<=nl; i++) {
	    if (((xgr[i]-(res / 2)) < x[k])&&(x[k]<= (xgr[i]+(res / 2)))) {
		for (j=1; j<=nc; j++) {
		    if (((ygr[j]-(res / 2)) < y[k])&&(y[k]<= (ygr[j]+(res / 2)))) {
			carte[i][j]++;
		    }
		}
	    }
	}
    }
}



/* **************************************************************
 *                                                              *
 *      lr = R interactive version of locrast                   *
 *                                                              *
 **************************************************************** */


void lr(double *xgri, double *ygri, double *xr, double *yr,
	double *carter, int *nco, int *nli, int *nlixy)
{

    /* Declaration of local variables, 
       and memory allocation */
    int i,j,k, ncg, nlg, nlxy;
    double *xgr, *x, *y, *ygr, **carte;
    
    ncg = *nco;
    nlg = *nli;
    nlxy = *nlixy;
    
    vecalloc(&xgr, nlg);
    vecalloc(&ygr, ncg);
    vecalloc(&x, nlxy);
    vecalloc(&y, nlxy);
    taballoc(&carte, nlg, ncg);
    
    /* Copy from R -> C variables */
    for (i=1; i<=nlxy; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xgr[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	ygr[i] = ygri[i-1];
    }
    
    /* Use of locrast */
    locrast(xgr, ygr, x, y, carte);
    
    /* Copy from C -> R variables */
    k=0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    carter[k] = carte[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(carte);
    freevec(xgr);
    freevec(ygr);
    freevec(x);
    freevec(y);
}



/* ****************************************************************
   *                                                              *
   *    getcarte is the C equivalent of getkasc                   *
   *                                                              *
   **************************************************************** */

void getcarte(double **carte, double **kasc, int *indicecarte)
{

    /* Definition of local variables */
    int i,j,k, ic, lgr, cgr;
    
    /* Memory Allocation */
    lgr = carte[0][0];
    cgr = carte[1][0];
    ic = *indicecarte;
    
    /* Copy from R -> C variables */
    k = 1;
    for (j=1; j<=cgr; j++) {
	for (i=1; i<=lgr; i++) {
	    carte[i][j] = kasc[k][ic];
	    k++;
	}
    }
}


/* ****************************************************************
   *                                                              *
   * gc is for the interactive version with R, just for a test    *
   *                                                              *
   **************************************************************** */


void gc(double *carter, double *kascr, int *nlgr, 
	int *ncgr, int *nhab)
{
    /* Declaration of local variables */
    int i,j,k, nlg, ncg, nh, nl;
    double **carte, **kasc;
    nlg = *nlgr;
    ncg = *ncgr;
    nh = *nhab;
    nl = nlg*ncg;
    
    /* Memory allocation */
    taballoc(&carte, nlg, ncg);
    taballoc(&kasc, nl, nh);
    
    /* Copy from R -> C variables */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nh; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    
    /* Use of getcarte */
    i=1;
    getcarte(carte, kasc, &i);
    
    /* Copy from C -> R variables */
    k=0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    carter[k] = carte[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(carte);
    freetab(kasc);
    
}



/* ****************************************************************
   *                                                              *
   * comptepasNA counts the number of rows of a table which do    *
   * contain missing values                                       *
   *                                                              *
   **************************************************************** */

void comptePasNA(double **tab, int *nombre)
{
    int i,nb, nl;
    nb = 0;

    nl = tab[0][0];
    
    for (i=1; i<=nl; i++) {
	if (fabs(tab[i][1] + 9999) > 0.000000001) {
	    nb++;
	}
    }
    
    *nombre = nb;
}


/* ****************************************************************
   *                                                              *
   * videNA deletes the rows with missing values in a table       *
   * (equivalent with kasc2df of R)                               *
   *                                                              *
   **************************************************************** */

void videNA(double **entree, double **sortie, int *idcons)
{
    /* Declaration of local variables */
    int i,j,k, nc, nl;
    nl = entree[0][0];
    nc = entree[1][0];
    
    /* Computation */
    k=1;
    for (i=1; i<=nl; i++) {
	if (fabs(entree[i][1] + 9999) > 0.000000001) {
	    idcons[k] = i;
	    for (j=1; j<=nc; j++) {
		sortie[k][j] = entree[i][j];
	    }
	    k++;
	}
    }
}






/* ****************************************************************
   *                                                              *
   *        niche to apply OMI analysis                           *
   *                                                              *
   **************************************************************** */

void niche(double **X, double **Y, double *eig, double **mar)
{
    /* Declaration of local variables */
    int i, j, k, nl, nh, na, rang;
    double **ut, **dis, *poidsli, *poidsani, solo, *ni, **inertie;
    
    /* Memory Allocation */
    nh = X[1][0];
    na = Y[1][0];
    nl = Y[0][0];
    
    taballoc(&ut, na, nh);
    taballoc(&dis, na, nh);
    taballoc(&inertie, nh, nh);
    vecalloc(&ni, na);
    vecalloc(&poidsli, nl);
    vecalloc(&poidsani, na);
    
    /* Centring and reduction */
    for (i=1; i<=nl; i++) {
	poidsli[i] = (double) 1/nl;
    }
    matmodifcn(X, poidsli);
    
    /* Numbering of the relocations per animal */
    for (i=1; i<=nl; i++) {
	for (j=1; j<=na; j++) {
	    ni[j] = ni[j] + Y[i][j];
	}
    }
    
    /* Computation of poidsani */
    solo = 0;
    for (j=1; j<=na; j++) {
	solo = ni[j] + solo; /* sums of the locs */
    }
    for (j=1; j<=na; j++) {
	poidsani[j] = ni[j] / solo;
    }
    
  
    /* mean of use and availability */
    for (k=1; k<=na; k++) {
	for (i=1; i<=nl; i++) {
	    for (j=1; j<=nh; j++) {
		dis[k][j] = dis[k][j] + ((double) X[i][j]/nl);
		ut[k][j] = ut[k][j] + (X[i][j]*Y[i][k]/ni[k]);
	    }
	}
    }
    
    /* Marginality */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    mar[i][j] = ut[i][j] - dis[i][j];
	}
    }
    
    
    /* Computation of the inertia */
    sqrvec(poidsani);
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    mar[i][j] = mar[i][j] * poidsani[i];
	}
    }
    
    /* eigenstructure */
    prodmatAtAB(mar, inertie);
    DiagobgComp(nh, inertie, eig, &rang);
    
    /* Marginality for the output */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    mar[i][j] = mar[i][j] / poidsani[i];
	}
    }
    
    
    /* Free Memory */
    freetab(ut);
    freetab(dis);
    freetab(inertie);
    freevec(poidsli);
    freevec(poidsani);
}






/* ****************************************************************
   *                                                              *
   * mvtfreeman: arguments = indices of the rows and columns      *
   * (in et jn), of the freeman direction (dir), and we get       *
   * the indices of rows and columns (in the vector np) after the *
   * move.                                                        *
   *                                                              *
   **************************************************************** */

void mvtfreeman(int *in, int *jn, int *dir, int *np)
{
    int i,j;
    i=*in;
    j=*jn;
    
    if ((*dir == 0) | (*dir == 1) | (*dir == 7)) 
	i++;
    if ((*dir == 3) | (*dir == 4) | (*dir == 5)) 
	i--;
    if ((*dir == 1) | (*dir == 2) | (*dir == 3)) 
	j++;
    if ((*dir == 5) | (*dir == 6) | (*dir == 7)) 
	j--;
    
    np[1]=i;
    np[2]=j;
}


/* ****************************************************************
   *                                                              *
   * algorithm of contour monitoring (suivi de contour) to get    *
   * contour polygon                                              *
   *                                                              *
   **************************************************************** */


void getcontour(double *grille, int *nlig, int *ncol, int *indicelig, 
		int *indicecol, int *lcont)
{
    /* Declaration of local variables */
    int i, j, k, nl, nc, *idlig, *idcol, *P0, *P1, fini, *np, dirprec, dir;
    int lidlig;
    double **x;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    vecintalloc(&P0,2);
    vecintalloc(&P1,2);
    vecintalloc(&np,2);
    taballoc(&x, nl,nc);
    vecintalloc(&idlig, *lcont);
    vecintalloc(&idcol, *lcont);
    
    /* Copy from R -> C variables */
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    x[i][j] = grille[k];
	    k++;
	}
    }
    
    /* Search the indices of the rows and columns
       of the first cell with a non-missing value */
    k=0;
    i=0;
    j=1;
    
    while (k==0) {
	if (i != nl) {
	    i++;
	}
	else {
	    i=1;
	    j++;
	}
	k = (int) x[i][j];
    }
    
    /* When it is found, the algorithm begins */
    idlig[1] = i;
    idcol[1] = j;
    lidlig = 1;
    P0[1] = i;
    P0[2] = j;
    dir = 4;
    
    fini = 0;
    k = 0;
  
    while (fini==0) {
	
	/* finds the next direction */
	while (k==0) {
	    dir = (dir + 1)%8;
	    mvtfreeman(&i, &j, &dir, np);
	    dirprec = (dir + 5)%8;
	    k = (int) x[np[1]][np[2]];
	}
	/* once found, stores the new coordinate */
	if (lidlig == 1) {
	    P1[1] = np[1];
	    P1[2] = np[2];
	}
	else {
	    /* P0 is the first point of the contour and P1, the last
	       Is the contour closed? */
	    if ((i==P0[1])&&(j==P0[2])&&(np[1]==P1[1])&&(np[2]==P1[2])) 
		fini =1;
	}
	
	/* If it is not, then stores the result, and perform the move
	 in the found direction */
	if (fini==0) {
	    lidlig++;
	    idlig[lidlig] = np[1];
	    idcol[lidlig] = np[2];
	    i = np[1];
	    j = np[2];
	    mvtfreeman(&i, &j, &dirprec, np);
	    k = (int) x[np[1]][np[2]];
	    dir = dirprec;
	}
    }
    
    /* Copy from C -> R variables */
    for (i=1; i<=lidlig; i++) {
	indicelig[i-1]=idlig[i];
	indicecol[i-1]=idcol[i];
    }
    
    /* Free Memory */
    freeintvec(idlig);
    freeintvec(idcol);
    freeintvec(P0);
    freeintvec(P1);
    freeintvec(np);
    freetab(x);
}



/* ****************************************************************
   *                                                              *
   * Interactive version of getcontour for R                      *
   *                                                              *
   **************************************************************** */

void lcontour(double *grille, int *nlig, int *ncol, int *lcont)
{
    /* Declaration of local variables*/
    int i, j, k, nl, nc, *P0, *P1, fini, *np, dirprec, dir;
    int lidlig;
    double **x, **vois;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    vecintalloc(&P0,2);
    vecintalloc(&P1,2);
    vecintalloc(&np,2);
    taballoc(&vois, 3,3);
    taballoc(&x, nl,nc);
    
    /* R objects -> C objects */
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    x[i][j] = grille[k];
	    k++;
	}
    }
    
    
    /* Search of the first cell for which the value is not NA */
    k=0;
    i=0;
    j=1;
    
    while (k==0) {
	if (i != nl) {
	    i++;
	}
	else {
	    i=1;
	    j++;
	}
	k = (int) x[i][j];
    }
    
    
    /* When found, performs the algorithm */
    lidlig = 1;
    P0[1] = i;
    P0[2] = j;
    dir = 4;
    



    
    fini = 0;
    k = 0;
    
    /* Same algorithm as in the previous function */
    
    while (fini==0) {
	while (k==0) {
	    dir = (dir + 1)%8;
	    mvtfreeman(&i, &j, &dir, np);
	    dirprec = (dir + 5)%8;
	    k = (int) x[np[1]][np[2]];
	}
	if (lidlig == 1) {
	    P1[1] = np[1];
	    P1[2] = np[2];
	}
	else {
	    if ((i==P0[1])&&(j==P0[2])&&(np[1]==P1[1])&&(np[2]==P1[2])) 
		fini = 1;
	}
	
	if (fini==0) {
	    lidlig++;
	    i = np[1];
	    j = np[2];
	    mvtfreeman(&i, &j, &dirprec, np);
	    k = (int) x[np[1]][np[2]];
	    dir = dirprec;
	}
    }
    
    /* and return to R */
    *lcont = lidlig;

    /* Free Memory */
    freeintvec(P0);
    freeintvec(P1);
    freeintvec(np);
    freetab(vois);
    freetab(x);
}




/* ****************************************************************
   *                                                              *
   *                 Gets the levels of a factor                  *
   *                                                              *
   **************************************************************** */


void levels(double *vec, double *lev, int *lvec)
{
    /* Declaration of local variables */
    int i,j,k,n, l;
    lev[1] = vec[1];
    k=1;
    n=*lvec;
  
    /* gets the levels */
    for (i=2; i<=n; i++) {
	l=0;
	for (j=1; j<=k; j++) {
	    if (fabs(vec[i] - lev[j]) < 0.000000001)
		l=1;
	}
	if (l==0) {
	    k++;
	    lev[k] = vec[i];
	}
    }
    *lvec = k;
}


/* ****************************************************************
   *                                                              *
   *             Sequential labelling of connex components        *
   *                                                              *
   **************************************************************** */


void seqeticorr(double *grille, int *nlig, int *ncol)
{
    /* Declaration of local variables */
    int i, j, k, l, m, n, o, nl, nc, pr, beta, nniv, eticour;
    double **x, *Tc, *prec, *tmp, *tmp1, *tmp2, *etcons, *lf;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    taballoc(&x, nl, nc);
    vecalloc(&Tc, nl*nc);
    
    /* R objects -> C objects */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    Tc[1]=1;
    eticour=1;
    
    for (j=2; j<=nc; j++) {
	for (i=2; i<=nl; i++) {
	    if (((int) x[i][j])!=0) {
		vecalloc(&prec, 4);
		prec[1] = x[i-1][j-1];
		prec[2] = x[i][j-1];
		prec[3] = x[i+1][j-1];
		prec[4] = x[i-1][j];
		
		k=0;
		for (l=1; l<=4; l++) {
		    if (((int) prec[l])!=0)
			k++;
		}
		
		/* k contains the number of non null predecessors */
		if (k!=0) {
		    vecalloc(&tmp, k); /* tmp contains the non null pred */
		    m=1;
		    for (l=1; l<=4; l++) {
			if (((int) prec[l])>0) {
			    tmp[m] = prec[l];
			    m++;
			}
		    }
		    
		    freevec(prec);
		    vecalloc(&prec, k);
		    for (l=1; l<=k; l++)
			prec[l] = tmp[l];
		    freevec(tmp);
		    /* Now, prec contains the non null preds */
		    
		    
		    
		    /* Number of levels of the factor prec */
		    vecalloc(&tmp1, 4);
		    m=k;
		    levels(prec, tmp1, &m);
		    /* m contains the number of levels */
		    vecalloc(&tmp2, m);
		    /* tmp2 contains the levels of prec
		       (equivalent of etiprec in R) */
		    for (l=1; l<=m; l++)
			tmp2[l]=tmp1[l];
		    freevec(tmp1);
		    
		    if (m == 1) {
			x[i][j] = tmp2[1];
		    } else {
			/* computation of the minimum 
			   level and storage in xij */
			x[i][j] = tmp2[1];
			for (l = 1; l <= m; l++) {
			    if (tmp2[l]<x[i][j])
				x[i][j] = tmp2[l];
			}
			
			/* etcons will contain the different labels of 
			   xij */
			vecalloc(&etcons, m-1);
			n=1;
			for (l=1; l<=m; l++) {
			    if (fabs(x[i][j] - tmp2[l]) > 0.000000001) {
				etcons[n] = tmp2[l];
				n++;
			    }
			}
			
			/* loop to fill the correspondence table */
			for (l=1; l<=(m-1); l++) {
			    pr = (int) etcons[l];
			    beta = pr;
			    while (((int) Tc[beta])!=beta) {
				o = (int) Tc[beta];
				Tc[beta] = Tc[(int) x[i][j]];
				beta = o;
			    }
			    Tc[beta] = Tc[(int) x[i][j]];
			}
			freevec(prec);
			freevec(tmp2);
			freevec(etcons);
		    }
		} else {
		    Tc[eticour] = eticour;
		    x[i][j]= eticour;
		    eticour++;
		}
	    }
	}
    }
    
    eticour--;
    
    /* Actualisation of the correspondence table */
    for (i=1; i<=eticour; i++) {
	j = i;
	while (((int) Tc[j])!=j)
	    j = (int) Tc[j];
	Tc[i] = j;
    }
    j=eticour;
    vecalloc(&tmp1, j);
    vecalloc(&tmp2, eticour);
    for (i=1; i<=eticour; i++) {
	tmp2[i]=Tc[i];
    }
    
    levels(tmp2, tmp1, &j);
    freevec(tmp2);
    
    vecalloc(&lf, j);
    for (i=1; i<=j; i++)
	lf[i]=tmp1[i];
    freevec(tmp1);
    nniv=j;
    
    /* Second pass */
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    if (fabs(x[i][j]) > 0.000000001) {
		x[i][j] = Tc[(int) x[i][j]];
	    }
	}
    }

    /* Last pass: levels varying from 1 to p */
    k = 1;
    for (j=1; j<=nniv; j++) {
      i = (int) lf[j];
      if (i != k) {
	for (l = 1; l <= nl; l++) {
	  for (m = 1; m <= nc; m++) {
	    if (((int) x[l][m]) == i)
	      x[l][m]=k;
	  }
	}
      }
      k++;
    }

    /* grid */
    k=0;
    for (i=1; i<=nl; i++) {
      for (j=1; j<=nc; j++) {
	grille[k]=x[i][j];
	k++;
      }
    }

    freetab(x);
    freevec(Tc);
    freevec(lf);
  }



/* ****************************************************************
   *                                                              *
   *   epa: bivariate normal kernel                               *
   *                                                              *
   **************************************************************** */

int selectptsbo(double *xl, double *yl, double *box, 
		int *indcons)
{
    int i,nl,cons,k;
    nl = (int) xl[0];
    
    k=0;
    for (i=1; i<=nl; i++) {
	cons = 0;
	if (xl[i] < box[1]) {
	    if (xl[i] > box[2]) {
		if (yl[i] < box[3]) {
		    if (yl[i] > box[4]) {
			cons = 1;
		    }
		}
	    }
	}
	if (cons == 1) {
	    k++;
	    indcons[k] = i;
	}
    }
    return(k);
}

void epa(double *X, double *Y, double *xl, double *yl, 
	 double *val, double *fen)
{
    /* Declaration of local variables */
    int k, nl, *indcons, ncons;
    double *xy, kx, di2, h, *box;
    
    /* Bases */
    nl = (int) xl[0];
    vecalloc(&xy, 2);
    vecintalloc(&indcons, nl);
    vecalloc(&box, 4);
    *val = 0;
    h = *fen;
    kx = 0;
    
    /* Keep only the points no further than 4*fen of the current pixel */
    box[1] = *X + (4 * h);
    box[2] = *X - (4 * h);
    box[3] = *Y + (4 * h);
    box[4] = *Y - (4 * h);
    ncons = selectptsbo(xl, yl, box, indcons);
        
    
    /* The bivariate normal kernel */
    if (ncons>0) {
	for (k=1; k<=ncons; k++) {
	    xy[1] = (xl[indcons[k]] - *X);
	    xy[2] = (yl[indcons[k]] - *Y);
	    di2 = xy[1]*xy[1] + xy[2]*xy[2];
	    
	    kx = exp(-di2/(2*h*h));
	    *val = *val + kx;
	}
	*val = *val * (1/(((double) nl)*h*h*2*3.14159265359));
    } else {
	*val=0;
    }
    freevec(xy);
    freeintvec(indcons);
    freevec(box);
}




/* ****************************************************************
   *                                                              *
   *             Kernel home range                                *
   *                                                              *
   **************************************************************** */


void kernelhr(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo)
{
    /* Declaration of local variables */
    int i, j, k, ncg, nlg, nlo;
    double **gri, *xg, *yg, *xl, *yl, X, Y, tmp;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;
    tmp = 0;
    
    taballoc(&gri,nlg, ncg);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
    
    /* R objects -> C objects */
  
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    /* loop on the grid */
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    X = xg[i];
	    Y = yg[j];
	    epa(&X, &Y, xl, yl, &tmp, fen);
	    gri[i][j] = tmp;
	}
    }
    
    /* C objects -> R objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }

    /* Memory Free */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(xl);
    freevec(yl);
}




/* ****************************************************************
   *                                                              *
   *   Epanechnikov estimation thanks to a kernel                 *
   *                                                              *
   **************************************************************** */

void epanechnikov(double *Xo, double *Yo, double *xg, double *yg, 
		  double *fen, double **grille, int nlo)
{
    /* Declaration of local variables */
    
    int i, j, ncg, nlg, imin, imax, jmin, jmax;
    double X, Y, h, *xgb, *ygb, tmp;
    
    nlg = xg[0];
    ncg = yg[0];
    h = *fen;
    X = *Xo;
    Y = *Yo;
    vecalloc(&xgb, nlg);
    vecalloc(&ygb, ncg);
    imin=0;
    jmin=0;
    imax=0;
    jmax=0;
    
    /* Computes again the values xg and yg */
    for (i=1; i<=nlg; i++) {
	xgb[i] = fabs(xg[i]-X);
	if (xgb[i] < h) {
	    if (imin == 0) {
		imin = i;
	    }
	}
	if (xgb[i] > h) {
	    if (imin != 0) {
		imax = i;
	    }
	}
    }
    for (i=1; i<=ncg; i++) {
	ygb[i] = fabs(yg[i]-Y);
	if (ygb[i] < h) {
	    if (jmin == 0) {
		jmin = i;
	    }
	}
	if (ygb[i] > h) {
	    if (jmin != 0) {
		jmax = i;
	    }
	}
    }
    
    for (i=imin; i<=imax; i++) {
	for (j=jmin; j<=jmax; j++) {
	    tmp = ( (xgb[i] / h) * (xgb[i] / h) ) + ( (ygb[j] / h) * (ygb[j] / h) );
	    if (tmp < 1) {
		grille[i][j] = grille[i][j] + 
		    2 * (1 - tmp) / (3.14159265359 * nlo *  h * h);
	    }
	}
    }
    
    freevec(xgb);
    freevec(ygb);
}


/* For R */

void kernepan(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo)
{
    /* Declaration */
    int i, j, k, ncg, nlg, nlo;
    double **gri, *xg, *yg, *xl, *yl, X, Y;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;

    
    taballoc(&gri,nlg, ncg);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
  
    /* R to C */
    
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
  
    /* Loop on the relocations */
    for (i=1; i<=nlo; i++) {
	X = xl[i];
	Y = yl[i];
	epanechnikov(&X, &Y, xg, yg, fen, gri, nlo);
    }
    
    /* C to R */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(xl);
    freevec(yl);
}



/* ****************************************************************
   *                                                              *
   *   Find Minimum LSCV                                          *
   *                                                              *
   **************************************************************** */


double L(double smooth, int nlo, int ndist, double *dists) 
{
  int ii;
  double resL,n;  

  n = (double) nlo;
  resL = 0.;
  
  for(ii=0; ii < ndist; ii++){ 
      resL+= (exp(-pow(dists[ii],2)/(4. * pow(smooth,2)))) - (4. * (exp(-pow(dists[ii],2)/(2. * pow(smooth,2.)))));
  }
  
  resL = 1./(3.14159265359 * pow(smooth,2.) * n) + (2*resL -3*n)/(3.14159265359 * 4. * pow(smooth,2.) * pow(n, 2.));

  return(resL);
}




double euclidean_distance(double x1, double y1, double x2, double y2)
{
    double out = 0.0;
    out = pow((x2-x1), 2) + pow((y2-y1), 2);
    return sqrt(out);
}


double comdi(double *x, double *y, double *dists, 
	     int n)
{
    int ii,jj,kk;



    kk=0;
    
    for(ii=1; ii <= n-1; ii++){
	for(jj=ii+1; jj<=n; jj++){
	    dists[kk] = euclidean_distance(x[ii], y[ii], x[jj],y[jj]);
	    kk++;
	}
    }
    return (kk);
}



void CVmise(int *nloc, double *xlo, double *ylo,
	    double *hvec, double *CV, int *nhteste)
{
    /* Declaration */
    int i, nlo, nh, ndist;
    double *xl, *yl, h, *dists;
    
    /* Allocation de mémoire */
    nlo = *nloc;
    nh = *nhteste;
    
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
    vecalloc(&dists, (nlo-1)*nlo);
    
    /* R to C */
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    /* Compute the distances */
    ndist=comdi(xl, yl, dists, nlo);
    
    /* Loop on the window of h */
    for (i=1; i<=nh; i++) {
	h = hvec[i-1];
	CV[i-1]=L(h, nlo, ndist, dists);
    }
	
    /* Free Memory */
    freevec(dists);
    freevec(xl);
    freevec(yl);
}

  



/* ****************************************************************
   *                                                              *
   *        Computation of the volume under the UD                *
   *                                                              *
   **************************************************************** */


void calcvolume(double *grille, int *ncolgri, int *nliggri, double *cellsize)
{
    /* Declaration */
    int i, j, k, nl, nc;
    double cs, **gri;
    
    /* Memory Allocation */
    nl = *nliggri;
    nc = *ncolgri;
    cs = *cellsize;
    
    taballoc(&gri, nl, nc);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    gri[i][j] = grille[k];
	    k++;
	}
    }
    
    /* Volume of the grid */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    gri[i][j] = gri[i][j]*cs*cs;
	}
    }
    
    /* C to R */
    k = 0;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free Memory */
    freetab(gri);
}





/* ****************************************************************
   *                                                              *
   *   DOMAIN: estimation of the potential distribution range     *
   *                                                              *
   **************************************************************** */

void calcsim(double *pix, double **pts, double *rg, 
	     int *nvar, int *npts, double *similarite)
{
    /* Declarations of variables */
    int no,nv, i, j;
    double *vecqual, *temp, nib;
    
    no = *npts;
    nv = *nvar;
    
    vecalloc(&vecqual, no);
    vecalloc(&temp, nv);
    
    /* Computation of the similarity */
    for (i=1; i<=no; i++) {
	nib = 0;
	for (j=1; j<=nv; j++) {
	    temp[j] = fabs(pix[j]-pts[i][j])/rg[j];
	    nib = nib + temp[j];
	}
	vecqual[i] = 1 - (1/((double) nv))*nib;
    }
    
    /* Computation of the habitat quality
       using the max of the similarity */
    *similarite = vecqual[1];
    
    for (i=2; i<=no; i++) {
	if (vecqual[i]>*similarite)
	    *similarite = vecqual[i];
    }
    
    
    /* Free memory */
    freevec(vecqual);
    freevec(temp);
    
}


/* For interaction with R */

void fctdomain(double *kascmod, double *ptsmod, 
	       double *range, int *npts, int *npix,  
	       int *nvar, double *qualhab)
{
    /* Declarations of variables */
    int no,np,nv, i, j, k;
    double **kasc, **pts, *rg, *pix, sim;
    
    /* Memory allocation */
    no = *npts;
    np = *npix;
    nv = *nvar;
    
    taballoc(&kasc, np, nv);
    taballoc(&pts, no, nv);
    vecalloc(&rg, nv);
    vecalloc(&pix, nv);
    
    /* R to C*/
    k = 0;
    for (i = 1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    kasc[i][j] = kascmod[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= no; i++) {
	for (j = 1; j <= nv; j++) {
	    pts[i][j] = ptsmod[k];
	    k++;
	}
    }
    
    for (i=1; i<=nv; i++) {
	rg[i] = range[i-1];
    }
    
    
    /* The core of the function */
    for (i=1; i<=np; i++) {
	for (j=1; j<=nv; j++) {
	    pix[j] = kasc[i][j];
	}
	calcsim(pix, pts, rg, &nv, &no, &sim);
	qualhab[i-1] = sim;
    }
    
    /* Free memory */
    freetab(kasc);
    freetab(pts);
    freevec(pix);
    freevec(rg);
    
}





/***********************************************************************
 *********                 Compositional analysis                 ******
 ***********************************************************************/


/* weighted mean lambda: compo analysis with R */

void wml(double **used, double **avail, double *wmla, int na, int nh,
	 double **proj1, double **proj2, double *nbassocie, int krep)
{
    /* Declaration of variables */
    double **dlr, *moydlr, *nadlr, **dlrtmp, **mod1, **mod2, **res1, **res2;
    double **SCEres1, **SCEres2, *vp1, *vp2, det1, det2, *vecalea, *aleamu;
    int i, j, k, idcol, *vecindice, rg1, rg2;
    int jb;
    
    /* Memory allocation */
    taballoc(&dlr, na, (nh*(nh-1)));
    taballoc(&mod1, na, (nh-1));
    taballoc(&mod2, na, (nh-1));
    taballoc(&SCEres1, (nh-1), (nh-1));
    taballoc(&SCEres2, (nh-1), (nh-1));
    taballoc(&dlrtmp, na, (nh-1));
    taballoc(&res1, na, (nh-1));
    taballoc(&res2, na, (nh-1));
    vecintalloc(&vecindice, nh-1);
    vecalloc(&nadlr, nh -1);
    vecalloc(&moydlr, nh-1);
    vecalloc(&vp1, nh-1);
    vecalloc(&vp2, nh-1);
    vecalloc(&aleamu, 2);
    vecalloc(&vecalea, na);
    
    aleamu[1] = 1;
    aleamu[2] = -1;
    
    jb = 0;
    
    /* Random permutation for each animal */
    for (i = 1; i <= na; i++) {
	aleapermutvec(aleamu);
	vecalea[i] = aleamu[1];
    }
    
    /* When krep == 1 : First repetition of the process: 
       Computation of the observed lambda */
    if (krep == 1) {
	for (i = 1; i<=na; i++) {
	    vecalea[i] = 1;
	}
    }
    
    
    /* empty nbassocie */
    for (i = 1; i <= nh; i++)
	nbassocie[i] = 0;
    
    
    /* loop to fill the DLR */
    for (k = 1; k <= nh; k++) {
	i = 1;
	
	/* build the vectors of indices */
	for (j = 1; j <= nh; j++) {
	    if (j != k) {
		vecindice[i] = j;
		i++;
	    }
	}
	
	/* Set the mean and the number of non-missing values to 0 */
	for (j = 1; j <= (nh-1); j++) {
	    moydlr[j] = 0;
	    nadlr[j] = 0;
	}
	
	/* First fill of the DLR */
	for (j = 1; j <= (nh-1); j++) {
	    jb = vecindice[j];
	    idcol = (nh - 1) * (k - 1) + j;
	    for (i = 1; i <= na; i++) {
		if ((fabs(avail[i][jb]) > 0.000000001)&&(fabs(avail[i][k]) > 0.000000001)) {
		    dlr[i][idcol] = (log(used[i][jb] / used[i][k]) - 
				     log(avail[i][jb] / avail[i][k])) * vecalea[i];
		    
		    /* computes the mean */
		    moydlr[j] = moydlr[j] + dlr[i][idcol];
		    nadlr[j]++;
		}
	    }
	}
	
	for (j = 1; j <= (nh-1); j++) {
	    moydlr[j] = moydlr[j] / nadlr[j];
	}
	
	/* Second loop: replace missing values */
	for (j = 1; j <= (nh-1); j++) {
	    idcol = (nh - 1) * (k - 1) + j;
	    jb = vecindice[j];
	    for (i = 1; i <= na; i++) {
		if ( (fabs(avail[i][jb]) < 0.000000001)||(fabs(avail[i][k])< 0.000000001))
		    dlr[i][idcol] = moydlr[j];
	    }
	}
	
	/* extraction of DLRtmp */
	for (i = 1; i <= na; i++) {
	    for (j = 1; j <= (nh-1); j++) {
		idcol = (nh - 1) * (k - 1) + j;
		dlrtmp[i][j] = dlr[i][idcol];
	    }
	}
	
	/* Computes the models */
	prodmatABC(proj1, dlrtmp, mod1);
	prodmatABC(proj2, dlrtmp, mod2);
	
	/* Computes the residuals */
	for (i = 1; i <= na; i++) {
	    for (j = 1; j <= nh-1; j++) {
		res1[i][j] = dlrtmp[i][j] - mod1[i][j];
		res2[i][j] = dlrtmp[i][j] - mod2[i][j];
	    }
	}
	
	/* The sum of squares */
	prodmatAtAB(res1, SCEres1);
	prodmatAtAB(res2, SCEres2);
	
	/* The determinant */
	DiagobgComp(nh-1, SCEres1, vp1, &rg1);
	DiagobgComp(nh-1, SCEres2, vp2, &rg2);
	det1 = 1;
	det2 = 1;
	
	if (rg1 != nh-1) {
	    det1 = 0;
	} else {
	    for (i = 1; i <= rg1; i++) {
		det1 = det1 * vp1[i];
	    }
	}
	
	if (rg2 != nh-1) {
	    det2 = 0;
	} else {
	    for (i = 1; i <= rg2; i++) {
		det2 = det2 * vp2[i];
	    }
	}

	wmla[k] = det1 / det2;
	for (i = 1; i <= (nh-1); i++)
	    nbassocie[k] = nbassocie[k] + nadlr[i];
    }
    
    /* free memory */
    freetab(dlr);
    freetab(mod1);
    freetab(mod2);
    freetab(SCEres1);
    freetab(SCEres2);
    freetab(dlrtmp);
    freetab(res1);
    freetab(res2);
    freeintvec(vecindice);
    freevec(nadlr);
    freevec(moydlr);
    freevec(vp1);
    freevec(vp2);
    freevec(aleamu);
    freevec(vecalea);
    
}




/* aclambda allows the computation of lambda in compositional analysis */

void aclambda(double *util, double *dispo, int *nani, int *nhab,  
	      double *xxtxxtmod1, double *xxtxxtmod2, double *rnv,
	      double *wmla, int *nrep, double *wm, double *nb)
{
    /* Declarations of variables */
    int na, nh, i, j, k, nr;
    double **ut, **di, **proj1, **proj2, *lilamb, *linb, sumnb;
    
    /* Memory allocation */
    na = *nani;
    nr = *nrep;
    nh = *nhab;
    
    taballoc(&ut, na, nh);
    taballoc(&di, na, nh);
    taballoc(&proj1, na, na);
    taballoc(&proj2, na, na);
    vecalloc(&lilamb, nh);
    vecalloc(&linb, nh);
    
    
    /* From R to C */
    /* use */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    ut[i][j] = util[k];
	    if (fabs(ut[i][j]) < 0.000000001)
		ut[i][j] = *rnv;
	    k++;
	}
    }
    
    /* availability */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    di[i][j] = dispo[k];
	    k++;
	}
    }
    
    /* projector 1 */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= na; j++) {
	    proj1[i][j] = xxtxxtmod1[k];
	    k++;
	}
    }
    
    /* projector 2 */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= na; j++) {
	    proj2[i][j] = xxtxxtmod2[k];
	    k++;
	}
    }
    
    
    /* Beginning of the loop */
    for (k = 1; k <= nr; k++) {
	/* weighted mean lambda */
	wml(ut, di, lilamb, na, nh,
	    proj1, proj2, linb, k);
	sumnb = 0;
	for (i = 1; i <= nh; i++)
	    sumnb = sumnb + linb[i];
	for (i = 1; i <= nh; i++)
	    wmla[k-1] = wmla[k-1] + ((lilamb[i] * linb[i]) / sumnb);
	
	/* C to R */
	if (k == 1) {
	    for (i = 1; i <= nh; i++) {
		wm[i-1] = lilamb[i];
		nb[i-1] = linb[i];
	    }
	}
    }
    
    
    /* Free memory */
    freetab(ut);
    freetab(di);
    freetab(proj1);
    freetab(proj2);
    freevec(lilamb);
    freevec(linb);
    
}


/* The ranking matrix for compositional analysis */

void rankma(double *used, double *avail, double *rankmap, double *rankmam,
	    double *rankmav, double *rankmanb, int *nhab, 
	    int *nani, int *nrep, double *rnv)
{
    /* Declarations of variables */
    int i, j, k, nh, na, nr, r;
    double **u, **a, **rmp, **rmm, **rmv, **rmnb;
    double *dlrtmp, *vecalea, val, moy;
    double *aleamu, **tabani;
    
    /* Memory Allocation */
    nh = *nhab;
    na = *nani;
    nr = *nrep;
    r = 0;
    
    taballoc(&u, na, nh);
    taballoc(&a, na, nh);
    taballoc(&rmv, nh, nh);
    taballoc(&rmp, nh, nh);
    taballoc(&rmm, nh, nh);
    taballoc(&rmnb, nh, nh);
    vecalloc(&dlrtmp, na);
    vecalloc(&vecalea, nr);
    vecalloc(&aleamu, 2);
    taballoc(&tabani, nr, na);
    aleamu[1] = -1;
    aleamu[2] = 1;
    
    /* Fill the table */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    u[i][j] = used[k];
	    a[i][j] = avail[k];
	    if (fabs(u[i][j]) < 0.000000001)
		u[i][j] = *rnv;
	    k++;
	}
    }
    
    /* Fill the table tabani */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= na; j++) {
	    aleapermutvec(aleamu);
	    tabani[i][j] = aleamu[1];
	}
    }
    for (i = 1; i<=na; i++) {
	tabani[1][i] = 1;
    }
    
    /* Beginning of the loop */
    for (k = 1; k <= nh; k++) {
	for (j = 1; j <= nh; j++) {
	    for (r = 1; r <= nr; r++) {
		moy = 0;
		/* Fills first the DLR per animal */
		for (i = 1; i <= na; i++) {
		    if ((fabs(a[i][j])> 0.000000001)&&(fabs(a[i][k]) > 0.000000001)) {
			dlrtmp[i] = (log(u[i][j]/u[i][k]) - log(a[i][j]/a[i][k])) * tabani[r][i];
			moy = moy + dlrtmp[i];
			if (r == 1)
			    rmnb[j][k]++;
		    }
		}
		
		/* Computes the mean */
		moy = moy / rmnb[j][k];
		if (r==1)
		    rmv[j][k] = moy;
		vecalea[r] = moy;
	    }
	    
	    /* Computes P */
	    val = rmv[j][k];

	    for (r = 1; r <= nr; r++) {
		if (val < vecalea[r])
		    rmm[j][k]++;
		if (val > vecalea[r])
		    rmp[j][k]++;
	    }
	}
    }
    
    /* C to R */
    k = 0;
    for (i=1; i<=nh; i++) {
	for (j=1; j<=nh; j++) {
	    rankmap[k] = rmp[i][j];
	    rankmam[k] = rmm[i][j];
	    rankmav[k] = rmv[i][j];
	    rankmanb[k] = rmnb[i][j];
	    k++;
	}
    }
    
    
    /* Free memory */  
    freetab(rmv);
    freetab(rmp);
    freetab(rmm);
    freetab(rmnb);
    freevec(dlrtmp);
    freetab(u);
    freetab(a);
    freevec(vecalea);
    freevec(aleamu);
    freetab(tabani);
}





/* ****************************************************************
   *                                                              *
   *   Morphological dilatation and erosion                       *
   *                                                              *
   **************************************************************** */

void erodil(double *grille, int *nlig, int *ncol, int *ntour, int *oper)
{
    /* declaration */
    int i,j,k,l,nl,nc, nt, etat0, etat1;
    double **x, **xm, *voisin;
    
    nl = *nlig;
    nc = *ncol;
    nt = *ntour;
    etat0 = 0;
    etat1 = 0;
    
    /* Memory alloocation */
    taballoc(&x,nl,nc);
    taballoc(&xm,nl,nc);
    vecalloc(&voisin, 9);
    
    /* R to C */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    /* Morphology */
    for (k=1; k<=nt; k++) {
	for (i=2; i<= (nl-1); i++) {
	    for (j=2; j<= (nc-1); j++) {
		voisin[1] = x[i-1][j-1];
		voisin[2] = x[i-1][j];
		voisin[3] = x[i-1][j+1];
		voisin[4] = x[i][j-1];
		voisin[5] = x[i][j+1];
		voisin[6] = x[i+1][j-1];
		voisin[7] = x[i+1][j];
		voisin[8] = x[i+1][j+1];
		voisin[9] = x[i][j];
		for (l=1;l<=9; l++) {
		    if (((int) voisin[l])==0) {
			etat0 = etat0 + 1;
		    } else {
			etat1 = etat1 + 1;
		    }
		}
		if (*oper==1) {
		    if (etat1 > 0)
			xm[i][j] = 1;
		    if (etat1 == 0)
			xm[i][j] = 0;
		} else {
		    if (etat0 == 0)
			xm[i][j] =1;
		    if (etat0 > 0)
			xm[i][j] =0;
		}
		etat1 = 0;
		etat0 = 0;
	    }
	}
    
	
	for (i=1; i<=nl; i++) {
	    for (j=1; j<=nc; j++) {
		x[i][j]=xm[i][j];
	    }
	}
    }
    
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    grille[k]=xm[i][j];
	    k++;
	}
    }
    
    /* Memory */
    freetab(x);
    freetab(xm);
    freevec(voisin);
    
}





/********************************************************
 x and y are the coordinates of the points, xp and yp 
 the coordinates of the vertices of the polygon (first and
 last vertices should be the same). deds is a vector of the
 same length as x and y: deds take the value 1 if the point is
 in the polygon and 0 otherwise
********************************************************/

void inout(double *x, double *y, double *xp, double *yp,
	   int *deds)
{
    /* Declaration of variables */
    int i, j, n, wm, np;
    double *xpc, *ypc, sig, a, b, x0;
    
    /* Memory allocation */
    n = x[0];
    np = xp[0];
    
    vecalloc(&xpc, np);
    vecalloc(&ypc, np);
    
    for (i = 1; i <= n; i++) {
	deds[i] = 1;
    }
    
    for (j = 1; j <= n; j++) {
	
	/* Centring on the point */
	for (i = 1; i <= np; i++) {
	    xpc[i] = xp[i] - x[j];
	    ypc[i] = yp[i] - y[j];
	}
	
	/* Number of intersections with X axis, for X >0 */
	wm = 0;
	for (i = 1; i <= (np-1); i++) {
	    sig = ypc[i] * ypc[i+1];
	    if (sig < 0) {
		/* The slope and intercept */
		/* Case 1: The slope is not infinite */
		if (fabs(xpc[i+1] - xpc[i]) > 0.000000001)
		{
		    a = (ypc[i+1] - ypc[i]) / (xpc[i+1] - xpc[i]);
		    b = (ypc[i]- a * xpc[i]);
		    /* value of x for y = 0 */
		    /* makes sense only if a != 0 */
		    if ((fabs(ypc[i+1] - ypc[i]) > 0.000000001)) {
			x0 = - b / a;
			if (x0 >= 0)
			    wm = abs(wm - 1);
		    } 	    
		}
		/* Case 2: Infinite slope
		   verify on the right of the point, i.e. 
		   xi >0 */
		if ((fabs(xpc[i+1] - xpc[i]) < 0.000000001))
		{
		    if (xpc[i] >= 0)
			wm = abs(wm - 1);
		}
	    }
	}
	
	/* If even number: outside. Inside otherwise */
	if (wm == 0)
	    deds[j] = 0;
    }
    
    
    /* Free memory */
    freevec(xpc);
    freevec(ypc);
}




/* verification of inout with R */

void inoutr(double *xr, double *yr, double *xpr, double *ypr,
	    int *dedsr, int *nxr, int *npr)
{
    /* Declaration */
    int i, nx, np, *deds;
    double *x, *y, *xp, *yp;
  
    /* Memory allocation */
    nx = *nxr;
    np = *npr;
    vecalloc(&x, nx);
    vecalloc(&y, nx);
    vecalloc(&xp, np);
    vecalloc(&yp, np);
    vecintalloc(&deds, nx);
    
    /* R to C */
    for (i = 1; i <= nx; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
    }
    
    for (i = 1; i <= np; i++) {
	xp[i] = xpr[i-1];
	yp[i] = ypr[i-1];
    }
    
    /* test of inout */
    inout(x, y, xp, yp, deds);
    
    /* C to R */
    for (i=1; i<=nx; i++) {
	dedsr[i-1] = deds[i];
    }
    
    /* Free memory */
    freevec(x);
    freevec(y);
    freevec(xp);
    freevec(yp);
    freeintvec(deds);
}




/***********************************************************
  Rasterization of a polygon: xp and yp are the coordinates
  of the polygon, xg and yg (not of the same length) are the
  coordinates of the rows and columns of the grid and 
  carte is a raster map.
*************************************************************/


void rastpol(double *xp, double *yp, double *xg, double *yg,
	     double **carte)
{
    /* Declaration of variables */
    int i, j, nl, nc, k, *deds;
    double *nxc, *nyc;
    
    /* Memory allocation */
    nl = xg[0];
    nc = yg[0];
    vecalloc(&nxc, nl*nc);
    vecalloc(&nyc, nl*nc);
    vecintalloc(&deds, nl*nc);
    
    /* Empties the map */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = 0;
	}
    }
    
    /* Output of the coordinates of the pixels of the grid */
    k = 1;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    nxc[k] = xg[i];
	    nyc[k] = yg[j];
	    k++;
	}
    }
    
    /* inout on these pixels */
    inout(nxc, nyc, xp, yp, deds);
    
    /* Fills the grid */
    k = 1;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = (double) deds[k];
	    k++;
	}
    }
    
    /* Free memory */
    freevec(nxc);
    freevec(nyc);
    freeintvec(deds);
}



/* ****************************************************************
   *                                                              *
   *   Verification of rastpol with R                             *
   *                                                              *
   **************************************************************** */

void rastpolaire(double *xpr, double *ypr, double *xgr, double *ygr,
		 double *carter, int *nlg, int *ncg, int *nvp)
{
    /* Declaration */
    int i, j, k, nl, nc, nv;
    double *xp, *yp, *xg, *yg, **carte;
    
    /* Memory allocation */
    nl = *nlg;
    nc = *ncg;
    nv = *nvp;
    vecalloc(&xp, nv);
    vecalloc(&yp, nv);
    vecalloc(&xg, nl);
    vecalloc(&yg, nc);
    taballoc(&carte, nl, nc);
    
    /* R to C */
    for (i = 1; i <= nv; i++) {
	xp[i] = xpr[i-1];
	yp[i] = ypr[i-1];
    }
    
    for (i = 1; i <= nl; i++) {
	xg[i] = xgr[i-1];
    }
    
    for (i = 1; i <= nc; i++) {
	yg[i] = ygr[i-1];
    }
    
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    carte[i][j] = carter[k];
	    k++;
	}
    }
    
    /* call to rastpol */
    rastpol(xp, yp, xg, yg, carte);
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    carter[k] = carte[i][j];
      k++;
	}
    }
    
    /* Free memory */
    freevec(xp);
    freevec(yp);
    freevec(xg);
    freevec(yg);
    freetab(carte);
}




/* ****************************************************************
   *                                                              *
   *   Computation of the marginality and tolérance (by variable) *
   *                                                              *
   **************************************************************** */


void calcniche(double **kasc, int *nvar, int *nlg, int *ncg,
	       double *margvar, double *tolvar, double **carte)
{
    /* definition of the variables */
    int i, j, l, nv, nc, nl, npixpol;
    double **cartevar;
    
    /* Memory allocation */
    nc = *ncg;
    nl = *nlg;
    nv = *nvar;
    npixpol = 0;

    
    taballoc(&cartevar, nl, nc);
    
    /* Marginality and tolerance set to 0 */
    for (l = 1; l <= nv; l++) {
	margvar[l] = 0;
	tolvar[l] = 0;
    }
    
    /* loop for each variable */
    for (l = 1; l <= nv; l++) {
	
	/* récupération de la carte */
	getcarte(cartevar, kasc, &l);
	npixpol = 0;
	
	/* The "used" mean */
	for (i = 1; i <= nl; i++) {
	    for (j = 1; j <= nc; j++) {
		if (fabs(carte[i][j] - 1) < 0.000000001) {
		    if (fabs(cartevar[i][j] + 9999) > 0.000000001) {
			margvar[l] = margvar[l] + cartevar[i][j];
			npixpol++;
		    }
		}
	    }
	}
	margvar[l] = margvar[l] / ((double) npixpol);
    
	/* The tolerance */
	for (i = 1; i <= nl; i++) {
	    for (j = 1; j <= nc; j++) {
		if (fabs(carte[i][j] - 1) < 0.000000001) {
		    if (fabs(cartevar[i][j] + 9999) > 0.000000001) {
			tolvar[l] = tolvar[l] + (cartevar[i][j] - margvar[l])*(cartevar[i][j] - margvar[l]);
		    }
		}
	    }
	}
	tolvar[l] = tolvar[l] / ((double) npixpol);
    }
    
    /* Free memory */
    freetab(cartevar);
}


/* ****************************************************************
   *                                                              *
   *   Test with R                                                *
   *                                                              *
   **************************************************************** */

void calcnicher(double *kascr, int *nvar, int *nlg, int *ncg,
		double *margvar, double *tolvar, double *carter)
{
    /* Declaration of variables */
    int i, j, k, nv, nl, nc;
    double *marg, *tol, **carte, **kasc;
    
    /* Memory Allocation */
    nv = *nvar;
    nl = *nlg;
    nc = *ncg;
    vecalloc(&marg, nv);
    vecalloc(&tol, nv);
    taballoc(&carte, nl,nc);
    taballoc(&kasc, nl*nc,nv);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= (nl*nc); i++) {
	for (j = 1; j<=nv; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j<=nc; j++) {
	    carte[i][j] = carter[k];
	    k++;
	}
    }
  
    
    /* Call to calcniche */
    calcniche(kasc, &nv, &nl, &nc, marg, tol, carte);
    
    /* C to R */
    for (i=1; i<=nv; i++) {
	margvar[i-1] = marg[i];
	tolvar[i-1] = tol[i];
    }
    
    /* Free memory */
    freevec(marg);
    freevec(tol);
    freetab(carte);
    freetab(kasc);
}






/* ****************************************************************
   *                                                              *
   *   Function to randomize orientation and position of the      *
   *   (coordinates xpr et ypr) on the study area (kascr)         *
   *                                                              *
   **************************************************************** */


void randompol(double *xpr, double *ypr, double *kascr,
	       double *marg, double *tol, int *nvar,
	       double *xgr, double *ygr, int *nlr, 
	       int *ncr, int *nvpr, int *nrep)
{
    /* definition of the variables */
    int i, j, k, l, r, np, nv, nc, nl, nvp, nr;
    double *xp, *yp, *xr, *yr, **kasc, **carte, *xg, *yg, *moyut;
    double *tolvar, **margb, **tolb, **ze;
    
    /* Memory allocation */
    nc = *ncr;
    nl = *nlr;
    nv = *nvar;
    np = nc*nl;
    nvp = *nvpr;
    nr = *nrep;
    r=1;
    
    vecalloc(&xp, nvp);
    vecalloc(&yp, nvp);
    vecalloc(&xr, nvp);
    vecalloc(&yr, nvp);
    vecalloc(&xg, nl);
    vecalloc(&yg, nc);
    vecalloc(&moyut, nv);
    vecalloc(&tolvar, nv);
    taballoc(&kasc, np, nv);
    taballoc(&carte, nl, nc);
    taballoc(&ze, nl, nc);
    taballoc(&margb, (nr+1), nv);
    taballoc(&tolb, (nr+1), nv);
    
    
    /* R to C */
    for (i=1; i<=nvp; i++) {
	xp[i] = xpr[i-1];
	yp[i] = ypr[i-1];
    }
    for (i=1; i<=nl; i++) {
	xg[i] = xgr[i-1];
    }
    for (i=1; i<=nc; i++) {
	yg[i] = ygr[i-1];
    }
    k=0;
    for (i=1; i<=np; i++) {
	for(j=1; j<=nv; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    
    /* Rasterisation of the polygon */
    rastpol(xp, yp, xg, yg, carte);
    
    /* used mean and tolerance */
    calcniche(kasc, &nv, &nl, &nc, moyut, 
	      tolvar, carte);
    
    for (l = 1; l <= nv; l++) {
	margb[1][l] = moyut[l];
	tolb[1][l] = tolvar[l];
    }
    
    /* C to R */
    for (i = 1; i <= nvp; i++) {
	xr[i] = xp[i];
	yr[i] = yp[i];
    }
    
    /* Randomization process */
    for (r =1; r <= nr; r++) {
	
	rotxy(xr, yr, r);
	rastpol(xr, yr, xg, yg, carte);
	l=1;
	getcarte(ze, kasc, &l);
	shr(carte, ze);
	calcniche(kasc, &nv, &nl, &nc, moyut, 
		  tolvar, ze);
	
	for (l = 1; l <= nv; l++) {
	    margb[r+1][l] = moyut[l];
	    tolb[r+1][l] = tolvar[l];
	}
	
    }
    
    /* C to R */
    k=0;
    for (i=1; i<=(nr+1); i++) {
	for (j=1; j<=nv; j++) {
	    marg[k]=margb[i][j];
	    tol[k]=tolb[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freevec(xp);
    freevec(yp);
    freevec(xg);
    freevec(yg);
    freetab(kasc);
    freetab(carte);
    freevec(xr);
    freevec(yr);
    freevec(moyut);
    freevec(tolvar);
    freetab(ze);
    freetab(margb);
    freetab(tolb);
}



/* ****************************************************************
   *                                                              *
   *   Spatial join: given a vector, an asc object the            *
   *   rows and columns coordinates of the map, and the cell size *
   *   the function returns the value at the point                *
   *                                                              *
   **************************************************************** */


void dedans(double *pts, double *xc, double *yc, double *na,
	    double cs, double **asc)
{
    int nl, nc, i, ligne, colo;
    double x, y;
  
    x = pts[1];
    y = pts[2];
    
    nl = xc[0];
    nc = yc[0];
    
    ligne = 0;
    colo = 0;
    
    for (i = 1; i <= nl; i++) {
	if (((xc[i] - cs/2) <= x) && ((xc[i] + cs/2) > x))
	    ligne = i;
    }
    
    for (i = 1; i <= nc; i++) {
	if (((yc[i] - cs/2) <= y) && ((yc[i] + cs/2) > y))
	    colo = i;
    }
    *na = asc[ligne][colo]; 
}


/* Test of dedans in R */

void dedansr(double *ptsr, double *xcr, double *ycr, double *na,
	     double *cs, double *ascr, int *nl, int *nc, int *nlocs)
{
    /* Declaration */
    int i,j,k;
    double *pts, *xc, *yc, **asc;

    /* Memory allocation */
    vecalloc(&pts, *nlocs);
    vecalloc(&xc, *nl);
    vecalloc(&yc, *nc);
    taballoc(&asc, *nl, *nc);
    
    /* R to C */
    pts[1] = ptsr[0];
    pts[2] = ptsr[1];
    
    for (i = 1; i <= *nl; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *nc; i++) {
	yc[i] = ycr[i-1];
    }
    k = 0;
    for (i = 1; i <= *nl; i++) {
	for (j = 1; j <= *nc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    /* Call to dedans */
    dedans(pts, xc, yc, na, *cs, asc);
    
    /* free memory */
    freevec(pts);
    freevec(xc);
    freevec(yc);
    freetab(asc);
}



/* ****************************************************************
   *                                                              *
   *   Randomization of a traject, based on the distances between *
   *   successive relocations and the time lag, as well as the    *
   *   angles between successive steps                            *
   *                                                              *
   **************************************************************** */

void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r)
{
    /* Declaration */
    int i, j, k, l, m, nsam, nltd, *index, *indangles, nlocs;
    double *pts, na, interv, *dobs, dech, anglech, ang;
    
    /* Memory allocation */
    vecalloc(&pts, 2);
    nltd = tabdist[0][0];
    nlocs = xp[0][0];
    vecintalloc(&indangles, (nlocs-2));
    
    /* fills the vector indangle */
    for (i = 1; i <= (nlocs - 2); i++) {
	indangles[i] = i;
    }
    
    /* 1. First loc of the traject */
    k = 0;
    ang = 0;
    
    while (k==0) {
	
	/* Random draw of relocations coordinates */
	xp[1][1] = (alea())*(rcx[2]-rcx[1]) + rcx[1];
	xp[1][2] = (alea())*(rcy[2]-rcy[1]) + rcy[1];
	
	pts[1] = xp[1][1];
	pts[2] = xp[1][2];
	
	/* Verifies that the loc is inside the study area */
	dedans(pts, xc, yc, &na, *cs, asc);
	if (fabs(na + 9999) > 0.000000001)
	    k = 1;
	
    }
    
    
    /* loop for the following relocations */
    for (i = 1; i <= (nlocs-1); i++) {
	interv = dt[i];
	
	/* How many distances for the observed time lag ? */
	nsam = 0;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001)
		nsam++;
	}
	
	/* Table of distances */
	vecalloc(&dobs, nsam);
	
	/* the vector index will be used to draw a random relocation */
	vecintalloc(&index, nsam);
	for (l = 1; l <= nsam; l++) {
	    index[l] = l;
	}
	
	/* In a first time, gets the corresponding distances */
	m = 1;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001) {
		dobs[m] = tabdist[j][2];
		m++;
	    }
	}
	
	k = 0;
	while (k == 0) {
	    /* Sampled Distance */
	    r = (int) (alea() * 100);
	    getpermutation(index, j * r);
	    dech = dobs[index[1]];
	    
	    /* Sampled Angles */
	    getpermutation(indangles, j * r);
	    anglech = angles[indangles[1]];
	    
	    /* update the angles */
	    ang = (ang + anglech);
	    
	    /* The new coordinates */
	    xp[i+1][1] = xp[i][1] + dech * cos(ang);
	    xp[i+1][2] = xp[i][2] + dech * sin(ang);
	    
	    pts[1] = xp[i+1][1];
	    pts[2] = xp[i+1][2];
	    
	    dedans(pts, xc, yc, &na, *cs, asc);
	    if (fabs(na + 9999) > 0.000000001)
		k = 1;
	}
	freevec(dobs);
	freeintvec(index);
    }
    
    /* Free memory */
    freeintvec(indangles);
    freevec(pts);
}







/* Test of rpath with R */

void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr)
{
    /* declaration of the variables */
    int i, j, k, r, nlocs, nltd;
    double **xp, *rcx, *rcy, **asc, **tabdist, *dt, *angles;
    double *xc,*yc, cs;
    
    /* Memory allocation */
    nlocs = *nlocsr;
    nltd = *nltdr;
    cs = *csr;
    
    taballoc(&xp, nlocs, 2);
    vecalloc(&rcx, 2);
    vecalloc(&rcy, 2);
    taballoc(&asc, *nlasc, *ncasc);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    taballoc(&tabdist, nltd, 2);
    vecalloc(&dt, nlocs-1);
    vecalloc(&angles, nlocs-2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    rcx[1] = rcrx[0];
    rcx[2] = rcrx[1];
    rcy[1] = rcry[0];
    rcy[2] = rcry[1];
    
    k = 0;
    for (i = 1; i <= *nlasc; i++) {
	for (j = 1; j <= *ncasc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nltd; i++) {
	for (j = 1; j <= 2; j++) {
	    tabdist[i][j] = tabdistr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= (nlocs - 1); i++) {
	dt[i] = dtr[i-1];
    }
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    for (i = 1; i <= (nlocs - 2); i++) {
	angles[i] = anglesr[i-1];
    }
    
    /* test of randpath */
    r = 1;
    rpath(xp, rcx, rcy, asc, tabdist, dt, angles,
	  xc, yc, &cs, r);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xpr[k] = xp[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(xp);
    freevec(rcx);
    freevec(rcy);
    freetab(asc);
    freetab(tabdist);
    freevec(dt);
    freevec(angles);
    freevec(xc);
    freevec(yc);
}




/* ****************************************************************
   *                                                              *
   *   joinkasc is the equivalent of join.kasc in R.              *
   *   Given a tableau of points, an object kasc, and one obtains *
   *   a table giving the composition of the map inder each point *
   *   (in res).                                                  *
   *                                                              *
   **************************************************************** */


void joinkasc(double **xp, double **kasc, double **res, int nl, int nc,
	      double *xc, double *yc, double *cs)
{
  int i,j, nlocs, nvar;
  double **carte, *pts, na;
  
  taballoc(&carte, nl, nc);
  vecalloc(&pts, 2);
  nlocs = xp[0][0];
  nvar = kasc[1][0];
  
  for (j = 1; j <= nvar; j++) {
    getcarte(carte, kasc, &j);
    for (i = 1; i <= nlocs; i++) {
      pts[1] = xp[i][1];
      pts[2] = xp[i][2];
      dedans(pts, xc, yc, &na, *cs, carte);
      res[i][j] = na;
    }
  }
  freetab(carte);
  freevec(pts);
}


/* Vérification de joinkasc sous R */

void joinkascr(double *xpr, double *kascr, int *nlasc, int *ncasc,
	       double *xcr, double *ycr, double *cs, int *nlocs,
	       int *nvar, double *resr)
{
    /* Declaration */
    int i,j,k;
    double **xp, **kasc, **res, *xc, *yc, cellsize;
    
    /* Memory allocation */
    taballoc(&xp, *nlocs, 2);
    taballoc(&res, *nlocs, *nvar);
    taballoc(&kasc, (*nlasc) * (*ncasc), *nvar);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    cellsize = *cs;
    
    /* R to C */
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    k = 0;
    for (i = 1; i <= ((*nlasc) * (*ncasc)); i++) {
	for (j = 1; j <= *nvar; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    k = 0;
    for (i = 1; i <= *nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    
    /* Call to joinkasc */
    joinkasc(xp, kasc, res,  *nlasc, *ncasc,
	     xc, yc, &cellsize);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nlocs; i++) {
	for (j = 1; j <= *nvar; j++) {
	    resr[k] = res[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(xp);
    freetab(res);
    freetab(kasc);
    freevec(xc);
    freevec(yc);
}




/* ****************************************************************
   *                                                              *
   *   Function allowing the randomization of a traject, based on *
   *   The distance and time lag between successive relocations   *
   *   and angles between successive steps, randomly drawn        *
   *   Randomization in the ecological space (measure marginality *
   *   and tolerance)                                             *
   *                                                              *
   **************************************************************** */

void randmargtol(double *xpr, double *rcrx, double *rcry, double *ascr, 
		 double *cwr, double *kascr, double *xcr, double *ycr,
		 double *csr,
		 double *tabdistr, double *dtr, double *anglesr, double *marr, 
		 double *tolr, int *nrepr, int *nlasc, 
		 int *ncasc, int *nvarr, int *nltdr, int *nlocsr)
{
    /* declaration of variables */
    int i, j, k, nr, r, nlocs, nltd, nvar;
    double **xp, *rcx, *rcy, **asc, **kasc, **tabdist, *dt, *angles;
    double *cw, *xc,*yc, cellsize, **res, *mar, *tol;
    
    /* Memory allocation */
    nr = *nrepr;
    nlocs = *nlocsr;
    nltd = *nltdr;
    nvar = *nvarr;
    cellsize = *csr;
    
    taballoc(&xp, nlocs, 2);
    taballoc(&res, nlocs, nvar);
    vecalloc(&rcx, 2);
    vecalloc(&rcy, 2);
    vecalloc(&mar, nvar);
    vecalloc(&tol, nvar);
    taballoc(&asc, *nlasc, *ncasc);
    taballoc(&kasc, (*nlasc)*(*ncasc), nvar);
    vecalloc(&cw, nvar);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    taballoc(&tabdist, nltd, 2);
    vecalloc(&dt, nlocs-1);
    vecalloc(&angles, nlocs-2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    rcx[1] = rcrx[0];
    rcx[2] = rcrx[1];
    rcy[1] = rcry[0];
    rcy[2] = rcry[1];
    
    k = 0;
    for (i = 1; i <= *nlasc; i++) {
	for (j = 1; j <= *ncasc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= ((*nlasc)*(*ncasc)); i++) {
	for (j = 1; j <= nvar; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nltd; i++) {
	for (j = 1; j <= 2; j++) {
	    tabdist[i][j] = tabdistr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= (nlocs - 1); i++) {
	dt[i] = dtr[i-1];
    }
    for (i = 1; i <= nvar; i++) {
	cw[i] = cwr[i-1];
    }
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    for (i = 1; i <= (nlocs - 2); i++) {
	angles[i] = anglesr[i-1];
    }
    
    /* observed marginality and tolerance */
    /* spatial join */
    joinkasc(xp, kasc, res,  *nlasc, *ncasc,
	     xc, yc, &cellsize);
    
    /* sets to 0 the vectors mar and tol */
    for (j = 1; j <= nvar; j++) {
	mar[j] = 0;
	tol[j] = 0;
    }
    
    /* 1. the means */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    mar[j] = mar[j] + (1 /((double) nlocs)) * res[i][j];
	}
    }
    /* 2. Centring of the table res */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    res[i][j] = res[i][j] - mar[j];
	}
    }
    /* 3. The variances */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    tol[j] = tol[j] + (1 /((double) nlocs)) * res[i][j] * res[i][j];
	}
    }
    /* 4. the tolerance and marginality */
    for (j = 1; j <= nvar; j++) {
	marr[0] = marr[0] + (mar[j] * mar[j] * cw[j]);
	tolr[0] = tolr[0] + (tol[j] * cw[j]);
    }
    
    /* Beginning of the randomization process */
    for (r = 2; r <= nr; r++) {
	/* creates a traject */
	rpath(xp, rcx, rcy, asc, tabdist, dt, angles, xc, yc, 
	      &cellsize, r);
	/* spatial join */
	joinkasc(xp, kasc, res,  *nlasc, *ncasc,
		 xc, yc, &cellsize);
	
	/* sets to 0 the vectors mar and tol */
	for (j = 1; j <= nvar; j++) {
	    mar[j] = 0;
	    tol[j] = 0;
	}
	
	/* 1. The means */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		mar[j] = mar[j] + (1 /((double) nlocs)) * res[i][j];
	    }
	}
	/* 2. Centring of the table res */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		res[i][j] = res[i][j] - mar[j];
	    }
	}
	/* 3. The variances */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		tol[j] = tol[j] + (1 /((double) nlocs)) * res[i][j] * res[i][j];
	    }
	}
	/* 4. tolerance and marginality */
	for (j = 1; j <= nvar; j++) {
	    marr[r-1] = marr[r-1] + (mar[j] * mar[j] * cw[j]);
	    tolr[r-1] = tolr[r-1] + (tol[j] * cw[j]);
	}
    }
    
    
    
    
    /* Free memory */
    freetab(xp);
    freevec(rcx);
    freevec(rcy);
    freevec(mar);
    freevec(tol);
    freetab(asc);
    freetab(kasc);
    freetab(tabdist);
    freevec(dt);
    freevec(angles);
    freevec(cw);
    freevec(xc);
    freevec(yc);
}



/* ****************************************************************
   *                                                              *
   *   Function allowing the randomization of the relocations on  *
   *   on an area independently                                   *
   *                                                              *
   **************************************************************** */

void rpoint(double **xp, double *rcx, double *rcy, double **asc, 
	    double *xc, double *yc, double *cs)
{
    /* Declaration */
    int i, k, nlocs;
    double *pts, na;
    
    /* Memory allocation */
    vecalloc(&pts, 2);
    nlocs = xp[0][0];
    
    /* For each loc */
    for (i = 1; i <= nlocs; i++) {
	k=0;
	while (k==0) {
	    
	    /* Random draw of the coordinates of the locs */
	    xp[i][1] = (alea())*(rcx[2]-rcx[1]) + rcx[1];
	    xp[i][2] = (alea())*(rcy[2]-rcy[1]) + rcy[1];
	    
	    pts[1] = xp[i][1];
	    pts[2] = xp[i][2];
	    
	    /* Is the loc in the study area? */
	    dedans(pts, xc, yc, &na, *cs, asc);
	    if (fabs(na + 9999) > 0.000000001)
		k = 1;
	    
	}
    }
    
    /* Free memory */
    freevec(pts);
}




/* ****************************************************************
   *                                                              *
   *   Function similar to randmargtol, but randomises locs       *
   *   instead of trajects                                        *
   *                                                              *
   **************************************************************** */


void randmargtolpts(double *xpr, double *rcrx, double *rcry, double *ascr, 
		    double *cwr, double *kascr, double *xcr, double *ycr,
		    double *csr, double *marr, double *tolr, 
		    int *nrepr, int *nlasc, 
		    int *ncasc, int *nvarr, int *nlocsr)
{
    /* declaration of variables */
    int i, j, k, nr, r, nlocs, nvar;
    double **xp, *rcx, *rcy, **asc, **kasc;
    double *cw, *xc,*yc, cellsize, **res, *mar, *tol;
    
    /* Memory allocation */
    nr = *nrepr;
    nlocs = *nlocsr;
    nvar = *nvarr;
    cellsize = *csr;
    
    taballoc(&xp, nlocs, 2);
    taballoc(&res, nlocs, nvar);
    vecalloc(&rcx, 2);
    vecalloc(&rcy, 2);
    vecalloc(&mar, nvar);
    vecalloc(&tol, nvar);
    taballoc(&asc, *nlasc, *ncasc);
    taballoc(&kasc, (*nlasc)*(*ncasc), nvar);
    vecalloc(&cw, nvar);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    rcx[1] = rcrx[0];
    rcx[2] = rcrx[1];
    rcy[1] = rcry[0];
    rcy[2] = rcry[1];
    
    k = 0;
    for (i = 1; i <= *nlasc; i++) {
	for (j = 1; j <= *ncasc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= ((*nlasc)*(*ncasc)); i++) {
	for (j = 1; j <= nvar; j++) {
	    kasc[i][j] = kascr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= nvar; i++) {
	cw[i] = cwr[i-1];
    }
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    
    /* Computation of observed values for maginality and tolerance */
    /* spatial join */
    joinkasc(xp, kasc, res,  *nlasc, *ncasc,
	     xc, yc, &cellsize);
    
    /* sets to 0 the vectors mar and tol */
    for (j = 1; j <= nvar; j++) {
	mar[j] = 0;
	tol[j] = 0;
    }
    
    /* 1. the means */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    mar[j] = mar[j] + (1 /((double) nlocs)) * res[i][j];
	}
    }
    /* 2. Centring of the table res */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    res[i][j] = res[i][j] - mar[j];
	}
    }
    /* 3. The variances */
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= nvar; j++) {
	    tol[j] = tol[j] + (1 /((double) nlocs)) * res[i][j] * res[i][j];
	}
    }
    /* 4. tolerance and marginality */
    for (j = 1; j <= nvar; j++) {
	marr[0] = marr[0] + (mar[j] * mar[j] * cw[j]);
	tolr[0] = tolr[0] + (tol[j] * cw[j]);
    }
    
    /* Beginning of the randomization process */
    for (r = 2; r <= nr; r++) {
	/* creates the locs */
	rpoint(xp, rcx, rcy, asc, xc, yc, &cellsize);
	/* spatial join */
	joinkasc(xp, kasc, res,  *nlasc, *ncasc,
		 xc, yc, &cellsize);
	
	/* sets to 0 the vectors mar and tol */
	for (j = 1; j <= nvar; j++) {
	    mar[j] = 0;
	    tol[j] = 0;
	}
	
	/* 1. the means */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		mar[j] = mar[j] + (1 /((double) nlocs)) * res[i][j];
	    }
	}
	/* 2. Centring of the table res */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		res[i][j] = res[i][j] - mar[j];
	    }
	}
	/* 3. The variances */
	for (i = 1; i <= nlocs; i++) {
	    for (j = 1; j <= nvar; j++) {
		tol[j] = tol[j] + (1 /((double) nlocs)) * res[i][j] * res[i][j];
	    }
	}
	/* 4. The tolerance and the marginality */
	for (j = 1; j <= nvar; j++) {
	    marr[r-1] = marr[r-1] + (mar[j] * mar[j] * cw[j]);
	    tolr[r-1] = tolr[r-1] + (tol[j] * cw[j]);
	}
    }
    
    /* Free memory */
    freetab(xp);
    freevec(rcx);
    freevec(rcy);
    freevec(mar);
    freevec(tol);
    freetab(asc);
    freetab(kasc);
    freevec(cw);
    freevec(xc);
    freevec(yc);
}



/* ********************************************************************
 *                                                                    *
 *            Diminish the resolution of a map                        *
 *                                                                    *
 * ******************************************************************** */



/* For factor maps */

void regroufacasc(double **asce, double **ascs, int *np,
		  int *nlev)
{
    /* declaration of variables */
    int i, j, k, l, m, dr, fr, dc, fc, nrs;
    int ncs, nl, *ll, max, vm, na, *vecmax, *vecmaxind;
    
    /* Memory allocation */
    nrs = ascs[0][0];
    ncs = ascs[1][0];
    nl = *nlev;
    vecintalloc(&ll, nl);


  
    /* loop to delete */
    for (i = 1; i <= nrs; i++) {
	for (j = 1; j <= ncs; j++) {
	    
	    /* extracts the corresponding subtable */
	    dr = (i-1)*(*np) + 1;
	    fr = i*(*np);
	    dc = (j-1)*(*np) + 1;
	    fc = j*(*np);
	    
	    /* empty ll */
	    for (m = 1; m <= nl; m++) {
		ll[m] = 0;
	    }
	    
	    /* One numbers the levels */
	    na = 0;
	    for (k = dr; k <= fr; k++) {
		for (l = dc; l <= fc; l++) {
		    if (fabs(asce[k][l] + 9999) > 0.000000001)
			ll[(int) asce[k][l]]++;
		    if (fabs(asce[k][l] + 9999) < 0.000000001)
			na++;
		}
	    }
	    
	    if (na != (*np)*(*np)) {
		/* One computes the maximum number */
		vm = ll[1];
		for (k = 2; k <= nl; k++) {
		    if (ll[k] > vm) {
			vm = ll[k];
		    }
		}
		
		/* ... and the number OF max */
		max = 0;
		for (k = 1; k <= nl; k++) {
		    if (ll[k] == vm) {
			max++;
		    }
		}
		
		/* one identifies the levels for which the number is max */
		vecintalloc(&vecmax, max);
		vecintalloc(&vecmaxind, max);
		for (l=1; l <=max; l++) {
		    vecmaxind[l] = l;
		}
		
		l = 1;
		for (k = 1; k<=nl; k++) {
		    if (ll[k] == vm) {
			vecmax[l] = k;
			l++;
		    }
		}
		
		/* Random sample of the levels in case of equality */
		if (max > 1) {		    
		    getpermutation(vecmaxind, i*j); /* random row */
		}
		ascs[i][j] = (double) vecmax[vecmaxind[1]];
		freeintvec(vecmax);
		freeintvec(vecmaxind);
		
	    } else {
		ascs[i][j] = -9999;
	    }
	    
	}
    }
    /* free memory */
    freeintvec(ll);
}



/* Regroufacasc version for R */

void regroufacascr(double *ascer, double *ascsr, int *npr,
		   int *nlevr, int *nle, int *nce, int *nls, 
		   int *ncs)
{
    /* Declaration of the variables */
    int i,j,k, np, nlev;
    double **asce, **ascs;
    
    /* Memory Allocation */
    np = *npr;
    nlev = *nlevr;
    taballoc(&asce, *nle, *nce);
    taballoc(&ascs, *nls, *ncs);
    
    /* R to C */
    k =0;
    for (i = 1; i <= *nle; i++) {
	for (j = 1; j <= *nce; j++) {
	    asce[i][j] = ascer[k];
	    k++;
	}
    }
    
    /* function */
    regroufacasc(asce, ascs, &np, &nlev);
    
    k =0;
    for (i = 1; i <= *nls; i++) {
	for (j = 1; j <= *ncs; j++) {
	    ascsr[k] = ascs[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(asce);
    freetab(ascs);
}




/* regrouascnum for numeric maps */

void regrouascnum(double **ascent, double **ascso)
{
    /* Declaration */
    int i, j, k, l, n, nle, nls, ncs, nreg;
    double moy, tmp;
    
    /* Definition of the variables */
    nle = ascent[0][0];
    nls = ascso[0][0];
    ncs = ascso[1][0];
    nreg = nle/nls;

    
    /* Computes the mean */
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    moy = 0;
	    n = 0;
	    for (k = 1; k <= nreg; k++) {
		for (l = 1; l <= nreg; l++) {
		    tmp = ascent[((i - 1) * nreg) + k][((j - 1) * nreg) + l];
		    if (fabs(tmp + 9999) > 0.000000001) {
			moy = tmp + moy;
		    }
		    if (fabs(tmp + 9999) < 0.000000001) {
			n++;
		    }
		}
	    }
	    if (n == (nreg * nreg)) {
		ascso[i][j] = -9999;
	    } else {
		ascso[i][j] = moy / (((double) (nreg * nreg))- ((double) n));
		
	    }
	}
    }
}


/* Version for R */

void regrouascnumr(double *ascentr, double *ascsor, double *nler, double *ncer,
		   double *nlsr, double *ncsr)
{
    /* Declaration of variables */
    int i, j, k, nle, nce, nls, ncs;
    double **ascent, **ascso;
    
    /* Memory Allocation */
    nle = *nler;
    nce = *ncer;
    nls = *nlsr;
    ncs = *ncsr;
    
    taballoc(&ascent, nle, nce);
    taballoc(&ascso, nls, ncs);
  
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nle; i++) {
	for (j = 1; j <= nce; j++) {
	    ascent[i][j] = ascentr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    ascso[i][j] = ascsor[k];
	    k++;
	}
    }
    
    /* procedure C */
    regrouascnum(ascent, ascso);
    
    /* C to R */
    k = 0;
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    ascsor[k] = ascso[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(ascso);
    freetab(ascent);
}







/* 
   All the factors should be labelled from 1 to n, without missing 
   values. Version for the objects of class kasc
*/


void regroukasc(double *kascr, double *kascniou, int *nrow, 
		int *ncol, int *nvar, int *npix,
		int *typer, int *nrniou, int *ncniou)
{
    /* declaration of variables */
    double **kasc, **asc, **ascn, **kascn;
    int i, j, k, l, nr, nc, nv, *typ, np, nrn, ncn, nl;
    
    
    /* Memory allocation */
    nr = *nrow;
    nc = *ncol;
    nv = *nvar;
    np = *npix;
    nrn = *nrniou;
    ncn = *ncniou;
    
    taballoc(&kasc, nr*nc, nv);
    taballoc(&kascn, nrn*ncn, nv);
    taballoc(&asc, nr, nc);
    taballoc(&ascn, nrn, ncn);
    vecintalloc(&typ, nv);
    
    /* R to C */
    for (i = 1; i<=nv; i++) {
	typ[i] = typer[i-1];
    }
    
    k = 0;
    for (i=1; i<= (nc*nr); i++) {
	for (j = 1; j<=nv; j++) {
	    kasc[i][j]=kascr[k];
	    k++;
	}
    }
    
    
    /* loop for each map */
    for (k=1; k<=nv; k++) {
	getcarte(asc, kasc, &k);
	if (typ[k] == 0)
	    regrouascnum(asc, ascn);
	nl = 0;
	if (typ[k] == 1) {
	    nl = (int) asc[1][1];
	    /* One counts the number of levels of the factor */
	    for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
		    if (((int) asc[i][j]) > nl)
			nl = (int) asc[i][j];
		}
	    }
	    regroufacasc(asc, ascn, &np, &nl);
	}
	
	
	l = 1;
	for (j = 1; j <= nc; j++) {
	    for (i = 1; i <= nr; j++) {
		kascn[l][k] = ascn[i][j];
	    }
	}
    }
    
    /* C to R */
    k=0;
    for (i = 1; i <= (nrn*ncn); i++) {
	for (j = 1; j <= nv; j++) {
	    kascniou[k] = kascn[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(kasc);
    freetab(asc);
    freetab(kascn);
    freetab(ascn);
    freeintvec(typ);
}






/* *******************************************

   Transformation of a square symetric matrix 
   into a matrix at the power -1/2
   
   ******************************************* */

void matmudemi(double **X, double **Y)
{
    /* Declaration of the variables */
    int i, j, nc, rg;
    double **U, **L, *lambda, **Ubis, **Uter;
    
    /* Memory Allocation */
    nc = X[0][0];
    taballoc(&U, nc, nc);
    taballoc(&Ubis, nc, nc);
    taballoc(&Uter, nc, nc);
    taballoc(&L, nc, nc);
    vecalloc(&lambda, nc);
    
    /* Fill Xtmp */
    for (i = 1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    U[i][j] = X[i][j];
	}
    }
    
    /* Eigenstructure of X */
    DiagobgComp(nc, U, lambda, &rg);
    
    /* Matrix lambda -1/2 */
    for (i = 1; i<=nc; i++) {
	L[i][i] = 1 / sqrt(lambda[i]);
    }
    
    /* Result */
    prodmatABC(U, L, Ubis);  
    for (i = 1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    Uter[i][j] = U[j][i];
	}
    }
    prodmatABC(Ubis, Uter, Y);
    
    /* Free memory */
    freetab(U);
    freetab(Ubis);
    freetab(Uter);
    freetab(L);
    freevec(lambda);
}



/* The same, for R */

void matmudemir(double *Xr, double *Yr, int *ncr)
{
    /* Declaration of variables */
    int i, j, k, nc;
    double **X, **Y;

    /* Memory allocation */
    nc = *ncr;
    taballoc(&X, nc, nc);
    taballoc(&Y, nc, nc);
    
    /* R to C */
    k = 0;
    for (i=1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    X[i][j] = Xr[k];
	    k++;
	}
    }
    
    /* matmudemi */
    matmudemi(X, Y);
    
    /* C to R */
    k = 0;
    for (i=1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    Yr[k] = Y[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(X);
    freetab(Y);
}





/* *******************************************

   Enfa

   ****************************************** */

void enfa(double **Z, double *p, int *nvar, int *npix,
	  double *vp)
{
    /* Declaration of local variables */
    double *m, *z, *y, **W, **Rs, **Rg, **Zbis, **Rsm12, norz;
    double **Wtmp, **H, **Iv, **yyt, **Ivmyyt, **Htmp;
    int i, j, nv, np, rg;
    
    /* Memory allocation */
    np = *npix;
    nv = *nvar;
    norz = 0;
    rg = 0;
    
    taballoc(&Zbis, np, nv);
    vecalloc(&m, nv);
    vecalloc(&z, nv);
    vecalloc(&y, nv);
    taballoc(&W, nv, nv);
    taballoc(&Iv, nv, nv);
    taballoc(&Ivmyyt, nv, nv);
    taballoc(&Htmp, nv, nv);
    taballoc(&yyt, nv, nv);
    taballoc(&H, nv, nv);
    taballoc(&Wtmp, nv, nv);
    taballoc(&Rg, nv, nv);
    taballoc(&Rs, nv, nv);
    taballoc(&Rsm12, nv, nv);
    
    /* Marginality */
    for (j = 1; j<=nv; j++) {
	for (i = 1; i <= np; i++){
	    m[j] = m[j] + p[i] * Z[i][j];
	}
    }
    
    /* Rs and Rg */
    for (i = 1; i<=np; i++) {
	for (j = 1; j <= nv; j++) {
	    Zbis[i][j] = Z[i][j] * sqrt(p[i]);
	}
    }
    prodmatAtAB(Zbis, Rs);
    for (i = 1; i<=np; i++) {
	for (j = 1; j <= nv; j++) {
	    Zbis[i][j] = Z[i][j] * sqrt((1/ ((double) np)));
	}
    }
    prodmatAtAB(Zbis, Rg);
    
    /* Rs^-1/2  */
    matmudemi(Rs,Rsm12);
    
    /* z */
    for (i = 1; i <= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    z[i] = z[i] + Rsm12[i][j] * m[j];
	}
    }
    
    /* norm of z */
    for (i = 1; i <= nv; i++) {
	norz = norz + (z[i] * z[i]);
    }
    norz = sqrt(norz);
    
    /* y */
    for (i = 1 ; i <= nv; i++) {
	y[i] = z[i] / norz;
    }
    
    /* W */
    prodmatABC(Rsm12, Rg, Wtmp);
    prodmatABC(Wtmp, Rsm12, W);
    
    
    /* **************************** */
    /* The large part: H            */
    /* **************************** */
    
    /* yyt */
    
    for (i = 1; i<= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    yyt[i][j] = y[i] * y[j];
	}
    }
    
    
    /* Iv */
    
    for (i = 1; i <= nv; i++) {
	Iv[i][i] = 1;
    }
    
    
    /* Ivmyyt */
    
    for (i = 1; i <= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    Ivmyyt[i][j] = Iv[i][j] - yyt[i][j];
	}
    }
    
    
    /* And finally, H */
    
    prodmatABC(Ivmyyt, W, Htmp);
    prodmatABC(Htmp, Ivmyyt, H);
    
    
    /* Eigenstructure of H */
    DiagobgComp(nv, H, vp, &rg);
    
    
    /* Free memory */
    freevec(m);
    freevec(z);
    freevec(y);
    freetab(W);
    freetab(Iv);
    freetab(Ivmyyt);
    freetab(Htmp);
    freetab(yyt);
    freetab(H);
    freetab(Wtmp);
    freetab(Rg);
    freetab(Rs);
    freetab(Rsm12);
    freetab(Zbis);
}


/* *****************************************************

ENFA with R

***************************************************** */

void enfar(double *Zr, double *pr, int *nvar, int *npix,
	   double *vpr)
{
    /* Declaration of variables */
    int i, j, k, np, nv;
    double **Z, *p, *vp;
    
    /* Memory allocation */
    taballoc(&Z, *npix, *nvar);
    vecalloc(&p, *npix);
    vecalloc(&vp, *nvar);
    
    np = *npix;
    nv = *nvar;
    
    /* R to C */
    k = 0;
    for (i=1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    Z[i][j] = Zr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= np; i++) {
	p[i] = pr[i-1];
    }
    
    /* ENFA ...*/
    enfa(Z, p, &nv, &np, vp);
    
    
    /* ... C to R ... */
    for (i = 1; i <= nv; i++) {
	vpr[i-1] = vp[i];
    }
    
    /* ... Free memory */
    freetab(Z);
    freevec(p);
    freevec(vp);
    
}



/* ********************************************************
   
Randomization in the ENFA: test of the first eigenvalue of specialization

******************************************************** */

void randenfa(double **Z, double *p, int *nrep, double *res)
{
    /* Declaration of variables */
    int i, j, k, nv, np, ntot;
    double *psim, *vp;
    
    /* Memory Allocation */
    np = Z[0][0];
    nv = Z[1][0];
    ntot = 0;
    vecalloc(&psim, np);
    vecalloc(&vp, nv);

    
    /* Counts the total number of points */
    for (i = 1; i <= np; i++) {
	ntot = ntot + p[i];
    }
    
    /* Beginning of the randomization porocess */
    for (k = 1; k <= *nrep; k++) {
	
	/* empty vector psim */
	for (i = 1; i <= np; i++) {
	    psim[i] = 0;
	}
	
	/* randomization of locs in the vector psim */
	for (i = 1; i <= ntot; i++) {
	    j = (int) (np * alea());
	    psim[j]++;
	}
	
	/* vector of weight */ 
	for (i = 1; i <= np; i++) {
	    psim[i] = psim[i] / ((double) ntot);
	}
	
	/* ... and ENFA */
	enfa(Z, psim, &nv, &np, vp); 
	
	/* storage in res... */
	res[k] = vp[1];
	
	/* ... and end of the loop */
    }
    
    /* Free memory */
    freevec(psim);
    freevec(vp);
}


/* The same but for external call from R */

void randenfar(double *Zr, double *pr, int *nvar, int *npix,
	       int *nrep, double *resr)
{
    /* Declaration of local variables */
    int i, j, k, nv, np, nr;
    double **Z, *p, *res;
    
    /* Memory Allocation */
    np = *npix;
    nv = *nvar;
    nr = *nrep;
    taballoc(&Z, np, nv);
    vecalloc(&p, np);
    vecalloc(&res, nr);
    
    /* R to C */
    k = 0;
    for (i=1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    Z[i][j] = Zr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= np; i++) {
	p[i] = pr[i-1];
    }
    
    /* C Function */
    randenfa(Z, p, &nr, res); 
    
    /* C to R */
    for (i = 1; i <= nr; i++) {
	resr[i-1] = res[i];
    }
    
    /* free memory */
    freevec(p);
    freevec(res);
    freetab(Z);
    
}



/* *********************************************************************
 *                                                                     *
 *                   Brownian bridge kernel                            *
 *                                                                     *
 ***********************************************************************/


/* Function normal 2D for brownian bridge */

void norm2d(double x1, double y1, double moyx, double moyy,
	    double var, double *res)
{
    double cste;
    cste = (1 / (2.0 * 3.141592653589793238 * var));
    cste = cste * exp( (-1.0 / (2.0 * var)) * (((x1 - moyx) * (x1 - moyx))+((y1 - moyy) * (y1 - moyy))));
    *res = cste;
}


double maxh(double sig1, double sig2, double *alpha, double maxt)
{
    int na,i;
    double res, tmp, a;
    
    res = 0;
    na = alpha[0];
    
    for (i = 1; i <= na; i++) {
	a = alpha[i];
	tmp = (maxt * a * (1 - a) * sig1) + 
	    ((pow(a,2) + pow((1-a), 2)) * sig2);
	if (tmp > res)
	    res = tmp;
    }
    return(sqrt(res));
}


double maxdt(double *T)
{
    int i,nt;
    double res;

    res = 0;
    nt = T[0];

    for (i = 2; i <= nt; i++) {
	if ((T[i]-T[i-1]) > res)
	    res = (T[i]-T[i-1]);
    }
    return(res);
}

/* keeps all the steps for which at least one relocation is 
   available in the box */
int consdanslabox(double *Xg, double **xy, 
		  int nl, int *indcons, double maxvh, int controlbox)
{
    int i,k,cons;
    double tmp1, tmp2, a, b;
    /* On a besoin d'une boucle sur les pas */
    k=0;
    
    for (i = 1; i<nl; i++) {
	
	cons = 0;
	
	if (xy[i][1] > (Xg[1] - (controlbox * maxvh)) ) {
	    if (xy[i][1] < (Xg[1] + (controlbox * maxvh)) ) {
		if (xy[i][2] > (Xg[2] - (controlbox * maxvh)) ) {
		    if (xy[i][2] < (Xg[2] + (controlbox * maxvh)) ) {
			cons = 1;
		    }
		}
	    }
	}
	if (xy[i+1][1] > (Xg[1] - (controlbox * maxvh)) ) {
	    if (xy[i+1][1] < (Xg[1] + (controlbox * maxvh)) ) {
		if (xy[i+1][2] > (Xg[2] - (controlbox * maxvh)) ) {
		    if (xy[i+1][2] < (Xg[2] + (controlbox * maxvh)) ) {
			cons = 1;
		    }
		}
	    }
	}
	
	if (cons == 0) {
	    a = (xy[i+1][2] - xy[i][2]) / (xy[i+1][1] - xy[i][1]);
	    b = xy[i+1][2] - a * xy[i+1][1];
	    tmp1 = a * (Xg[1] - (controlbox * maxvh)) + b;
	    tmp2 = a * (Xg[1] + (controlbox * maxvh)) + b;
	    
	    if (tmp1 <= (Xg[2] + (controlbox * maxvh))) {
		if (tmp1 >= (Xg[2] - (controlbox * maxvh))) {
		    cons = 1;
		}
	    }
	    
	    if (tmp2 <= (Xg[2] + (controlbox * maxvh))) {
		if (tmp2 >= (Xg[2] - (controlbox * maxvh))) {
		    cons = 1;
		}
	    }
	}
	
	
	if (cons==1) {
	    k++;
	    indcons[k]=i;
	}
	
    }
    
    return(k);
}



/* Integral of norm2d on alpha */
void integrno(double *XG, double *X1, double *X2, 
	      double *T, double *sig1,
	      double *sig2, double *alpha, double *res)
{
    /* Declaration */
    int i, na;
    double *val, tmp, *XX, var, nx1, ny1, nx2, ny2, ny, moyx, moyy, al;
    
    /* Memory allocation */
    na = alpha[0];
    vecalloc(&val, na);
    vecalloc(&XX, 2);
    
    XX[1] = X2[1] - X1[1];
    XX[2] = X2[2] - X1[2];
    *res = 0;
    
    
    /* loop for the computation of the value */
    for (i = 1; i<= na; i++) {
	al = alpha[i];
	
	var = (*T) * al * (1.0 - al) * (*sig1);
	var = var + (((al * al) + ((1.0 - al) * (1.0 - al))) * (*sig2));
	
	moyx = X1[1] + al * XX[1];
	moyy = X1[2] + al * XX[2];
	
	norm2d(XG[1], XG[2], moyx, moyy, var, &tmp);
	
	val[i] = tmp;
    }
    
    /* loop for the computation of the integral */
    for (i = 2; i<= na; i++) {
	nx1 = alpha[i-1];
	ny1 = val[i-1];
	nx2 = alpha[i];
	ny2 = val[i];
	ny = ny1;
	if (ny2 <= ny1)
	    ny = ny2;
	*res = *res + (nx2 - nx1) * (ny + (fabs(ny2 - ny1) / 2));
    }
    
    /* Free memory */
    freevec(val);
    freevec(XX);
}




/* Computes UD at a node of the grid */
void udbbnoeud(double *XG, double **XY, double *T, double *sig1,
	       double *sig2, double *alpha, double *res, int ncons, 
	       int *indcons)
{
    /* Declaration */
    int i, nlo;
    double *Xtmp1, *Xtmp2, dt, poids, dttot, tmp;
    
    /* Memory allocation */
    vecalloc(&Xtmp1, 2);
    vecalloc(&Xtmp2, 2);
    nlo = XY[0][0];
    dttot = T[nlo] - T[1];
    *res = 0;
    
    /* for each step */
    for (i = 1; i <= ncons; i++) {
	
	/* Computes weights and time lags */
	dt = T[indcons[i]+1] - T[indcons[i]];
	poids = dt / dttot;
	
	/* Output of the relocation values at i, and use of the function integrno */
	Xtmp1[1] = XY[indcons[i]][1];
	Xtmp1[2] = XY[indcons[i]][2];
	Xtmp2[1] = XY[indcons[i]+1][1];
	Xtmp2[2] = XY[indcons[i]+1][2];
	
	integrno(XG, Xtmp1, Xtmp2, &dt, sig1, sig2, alpha, &tmp);
	*res = *res + (poids * tmp);
    }
    freevec(Xtmp1);
    freevec(Xtmp2);
}



/* Main Function */
void kernelbb(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *sig1, double *sig2, 
	      double *xlo, double *ylo, double *Tr, int *controlbox, 
	      int *nalpha)
{
    /* Declaration */
    int i, j, k, ncg, nlg, nlo, *indcons, ncons;
    double **gri, *xg, *yg, **XY, tmp, *alpha, *Xgr, *T, res, vol;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;
    tmp = 0;
    
    taballoc(&gri,nlg, ncg);
    taballoc(&XY, nlo, 2);
    vecalloc(&xg, nlg);
    vecalloc(&T, nlo);
    vecalloc(&yg, ncg);
    vecalloc(&Xgr, 2);
    vecalloc(&alpha, *nalpha);
    vecintalloc(&indcons, nlo);
    
    /* R to C */    
    for (i=1; i<=nlo; i++) {
	XY[i][1] = xlo[i-1];
	XY[i][2] = ylo[i-1];
	T[i] = Tr[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    /* Build the vector alpha */
    alpha[1] = 0;
    for (i = 2; i <= *nalpha; i++) {
	alpha[i] = ((double) i) / ((double) *nalpha);
    }
    

    
	
    /* Loop on the grid */
    vol = 0;
    res = xg[2] - xg[1];
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    Xgr[1] = xg[i];
	    Xgr[2] = yg[j];
	    /*
	      ncons = consdanslabox(Xgr, XY, nlo, indcons, maxvh, *controlbox);
	    */
	    ncons = nlo-1;
	    for (k = 1; k < nlo; k++)
		indcons[k]=k;
	    udbbnoeud(Xgr, XY, T, sig1, sig2, alpha, &tmp, ncons, indcons);
	    gri[i][j] = tmp;
	    vol+=tmp;
	}
    }
    
    /* Standardization of the volume */

    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    gri[i][j] =  gri[i][j] / (vol * pow(res,2));
	}
    }
  
    /* C to R */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(T);
    freetab(XY);
    freevec(Xgr);
    freevec(alpha);
    freeintvec(indcons);
}


/* *********************************************************************
   Maximisation of the likelihood for the Brownian bridge
   
   *********************************************************************/
void CVL(double *xyr, double *Tr, 
	 int *nloc, double *Lr, double *sigma, int *nsig, double *sigma2)
{
    int i, j, k, nlo, ns, r;
    double **xy, *T,ai,*mui,sigmai,res;
    
    nlo = *nloc;
    ns = *nsig;


    taballoc(&xy, nlo, 2);
    vecalloc(&T, nlo);
    vecalloc(&mui, 2);
    
    /* C to R */
    k = 0;
    for (i=1; i <= nlo; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
	T[i] = Tr[i-1];
    }
    
    /* Indices of odd locations */
    for (r = 1; r <= ns; r++) {
	Lr[r-1] = 0;
	k=1;
	for (i=1; i < nlo; i++) {
	    if (k == 2) {
		ai = (T[i] - T[i-1])/(T[i+1] - T[i-1]);
		
		mui[1] = xy[i-1][1] + ai * (xy[i+1][1] - xy[i-1][1]);
		mui[2] = xy[i-1][2] + ai * (xy[i+1][2] - xy[i-1][2]);
		
		sigmai = ((T[i+1]-T[i-1]) * ai * (1-ai) * sigma[r-1]) + (pow((1 - ai),2) * (*sigma2)) + (pow(ai,2) * (*sigma2));
		
		norm2d(xy[i][1], xy[i][2], 
		       mui[1], mui[2], sigmai, &res);
		Lr[r-1] = Lr[r-1] + log(res);
		k=1;
	    }
	    k++;
	}
    }

    /* Free memory */
    freetab(xy);
    freevec(T);
    freevec(mui);
}




/* *********************************************************************
 *                                                                     *
 *                   Buffer on a line                                  *
 *                                                                     *
 ***********************************************************************/

/* given a line, ligpoly returns a buffer polygon containing the line */

void ligpoly(double *x, double *y, double r, double *xp, double *yp)
{
    /* Declaration */
    double x1, x2, y1, y2, xx, yy, alpha, beta, xim, xsm, yim, ysm, gamma;
    double xip, xsp, yip, ysp;
    
    x1 = x[1];
    x2 = x[2];
    y1 = y[1];
    y2 = y[2];
    xx = x2 - x1;
    yy = y2 - y1;
    
    alpha = atan(yy/xx);
    beta = alpha - (3.1415926/2);
    xim = x1 + r * (cos(beta));
    xsm = x2 + r * (cos(beta));
    yim = y1 + r * (sin(beta));
    ysm = y2 + r * (sin(beta));
    
    gamma = alpha + (3.1415926/2);
    xip = x1 + r * (cos(gamma));
    xsp = x2 + r * (cos(gamma));
    yip = y1 + r * (sin(gamma));
    ysp = y2 + r * (sin(gamma));
    
    xp[1] = xim;
    xp[2] = xsm;
    xp[3] = xsp;
    xp[4] = xip;
    xp[5] = xim;
    
    yp[1] = yim;
    yp[2] = ysm;
    yp[3] = ysp;
    yp[4] = yip;
    yp[5] = yim;
    
}


/* main function */

void buflig(double **x, double r, double **carte, double *xg, double *yg)
{
    /* Declaration */
    int i, j, k, nloc, nr, nc;
    double **x1, **x2, *xl, *yl, *xp, *yp, **cartebis;
    
    /* Memory allocation */
    nloc = x[0][0];
    k = 0;
    nr = carte[0][0];
    nc = carte[1][0];
    
    vecalloc(&xl, 2);
    vecalloc(&yl, 2);
    vecalloc(&xp, 5);
    vecalloc(&yp, 5);
    taballoc(&x1, nloc-1, 2);
    taballoc(&x2, nloc-1, 2);
    taballoc(&cartebis, nr, nc);
    
    /* Creates the two tables */
    for (i = 1; i <= nloc; i++) {
	if (i > 1) {
	    x2[i-1][1] = x[i][1];
	    x2[i-1][2] = x[i][2];
	}
	if (i < nloc) {
	    x1[i][1] = x[i][1];
	    x1[i][2] = x[i][2];
	}
    }
    
    /* Sets the map to 0 */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = 0;
	}
    }
    
    
    /* Buffer around the line */
    for (i = 1; i <= (nloc-1); i++) {
	xl[1] = x1[i][1];
	xl[2] = x2[i][1];
	yl[1] = x1[i][2];
	yl[2] = x2[i][2];
	
	ligpoly(xl, yl, r, xp, yp);
	
	rastpol(xp, yp, xg, yg, cartebis);
	
	for (j = 1; j <= nr; j++) {
	    for (k = 1; k <= nc; k++) {
		carte[j][k] = cartebis[j][k] + carte[j][k];
	    }
	}
    }
    
    /* Free memory */
    freevec(xl);
    freevec(yl);
    freevec(xp);
    freevec(yp);
    freetab(x1);
    freetab(x2);
    freetab(cartebis);
    
}


/* For external call from xithin R */
void bufligr(double *xr, double *rr, double *carter, 
	     double *xgr, double *ygr, int *nlr, int *ncr, 
	     int *nlocr)
{
    /* Declaration */
    int i, j, k, nc, nl, nloc;
    double **x, r, **carte, *xg, *yg;
    
    /* Memory allocation */
    nc = *ncr;
    nl = *nlr;
    nloc = *nlocr;
    r = *rr;
    
    taballoc(&x, nloc, 2);
    taballoc(&carte, nl, nc);
    vecalloc(&xg, nl);
    vecalloc(&yg, nc);
    
    
    /* R to C */
    k = 0;
    for (i=1; i<= nl; i++) {
	for (j = 1; j<=nc; j++) {
	    carte[i][j]=carter[k];
	    k++;
	}
    }
    
    k = 0;
    for (i=1; i<= nloc; i++) {
	for (j = 1; j<=2; j++) {
	    x[i][j]=xr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= nl; i++) {
	xg[i] = xgr[i-1];
    }
    
    for (i = 1; i <= nc; i++) {
	yg[i] = ygr[i-1];
    }
    
    /* Main function */
    buflig(x, r, carte, xg, yg);
    
    /* C to R */
    k = 0;
    for (i=1; i<= nl; i++) {
	for (j = 1; j<=nc; j++) {
	    carter[k]=carte[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(x);
    freetab(carte);
    freevec(xg);
    freevec(yg);
    
}


/* Computes Euclidean distances from a map of class asc */

void distxy(double **xy1, double **xy2, double *di)
{
    /* Declaration */
    int i, j, n1, n2;
    double *dib, mi;
    
    /* Memory allocation */
    n1 = xy1[0][0];
    n2 = xy2[0][0];
    
    vecalloc(&dib, n2);
    
    /* Euclidean distances */
    for (i = 1; i <= n1; i++) {
	for (j = 1; j <= n2; j++) {
	    dib[j] = sqrt( ((xy1[i][1] - xy2[j][1]) * (xy1[i][1] - xy2[j][1])) + 
			   ((xy1[i][2] - xy2[j][2]) * (xy1[i][2] - xy2[j][2])));
	}
	mi = dib[1];
	for (j = 2; j <= n2; j++) {
	    if (mi > dib[j]) {
		mi = dib[j];
	    }
	}
	di[i] = mi;
    }
    
    /* Free memory */
    freevec(dib);
}


/* For external call from R */
void distxyr(double *xy1r, double *xy2r, int *n1r, 
	     int *n2r, double *dire)
{
    /* Declaration */
    int i, j, k, n1, n2;
    double **xy1, **xy2, *di;
    
    /* Memory allocation */
    n1 = *n1r;
    n2 = *n2r;
    
    taballoc(&xy1, n1, 2);
    taballoc(&xy2, n2, 2);
    vecalloc(&di, n1);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= n1; i++) {
	for (j = 1; j <= 2; j++) {
	    xy1[i][j] = xy1r[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= n2; i++) {
	for (j = 1; j <= 2; j++) {
	    xy2[i][j] = xy2r[k];
	    k++;
	}
    }
    
    /* The function */
    distxy(xy1, xy2, di);
    
    /* C to R */
    for (i = 1; i <= n1; i++) {
	dire[i-1] = di[i];
    }
    
    /* Free memory */
    freetab(xy1);
    freetab(xy2);
    freevec(di);
}





/* *********************************************************************
 *                                                                     *
 *                   First passage time                                *
 *                                                                     *
 ***********************************************************************/

/* compute the distance cetween two points */
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di)
{
  *di = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}



/* compute the FPT for ONE relocation */
void fptt(double *x, double *y, double *t, int pos, double radius, double *fptto, int nlo)
{
    /* Declaration */
    int ok, pos2, naar, naav, na;
    double di, dt, dt2, di2, fptar, fptav;
    
    ok = 0;
    di = 0;
    di2 = 0;
    dt = 0;
    dt2 = 0;
    naar = 1;
    naav = 1;
    fptar = 0;
    fptav = 0;
  
  
    /* Search of the first loc outside the circle (before) */
    pos2 = pos;
    while (ok == 0) {
	pos2 = pos2 - 1;
	if (pos2 > 0) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naar = 0;
	}
    }
    
    /* computes the linear approximation */
    if (naar > 0) {
	dt = fabs(t[pos] - t[pos2]);
	dt2 = fabs(t[pos] - t[(pos2+1)]);
	dtmp(x[(pos2+1)], x[pos], y[(pos2+1)], y[pos], &di2);
	fptar = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    
    /* Search of the first loc outside the circle (after) */
    pos2 = pos;
    ok = 0;
    while (ok == 0) {
	pos2 = pos2 + 1;
	if (pos2 <= nlo) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naav = 0;
	}
    }
    
    /* Computes linear approximation */
    if (naav > 0) {
	dt = fabs(t[pos2] - t[pos]);
	dt2 = fabs(t[(pos2-1)] - t[pos]);
	dtmp(x[(pos2-1)], x[pos], y[(pos2-1)], y[pos], &di2);
	fptav = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    na = naar * naav;
    if (na > 0) {
	*fptto = fptar + fptav;
    } else {
	*fptto = -1;
    }
    
}



/* Computes the FPT for all relocations */
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt)
{
    /* Declaration */
    int i, j;
    double val;
    
    /* Computes the FPT */
    for (i=1; i<=nra; i++) {
	for (j=1; j<=nlo; j++) {
	    fptt(x, y, t, j, rad[i], &val, nlo);
	    fpt[j][i] = val;
	}
    }
}

/* for external call from within R */
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr)
{
    /* Declaration */
    int i, j, k, nlo, nra;
    double *x, *y, *t, *rad, **fpt;
    
    /* Memory allocation */
    nlo = *nlocs;
    nra = *nrad;
    
    vecalloc(&x, nlo);
    vecalloc(&y, nlo);
    vecalloc(&t, nlo);
    vecalloc(&rad, nra);
    taballoc(&fpt, nlo, nra);
  
    /* R to C */
    for (i = 1; i <= nlo; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	t[i] = tr[i-1];
    }
    
    for (i = 1; i <= nra; i++) {
	rad[i] = radius[i-1];
    }
    
    /* main function */
    fipati(x,y,t, nlo, nra, rad, fpt);
    
    /* C to R */
    k = 0;
    for (i=1; i<= nlo; i++) {
	for (j = 1; j<=nra; j++) {
	    fptr[k]=fpt[i][j];
	    k++;
	}
    }
    
    /* free memory */
    freetab(fpt);
    freevec(x);
    freevec(y);
    freevec(t);
    freevec(rad);
}






/* *********************************************************************
 *                                                                     *
 *                   Percolation cluster                               *
 *                                                                     *
 ***********************************************************************/


void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm)
{
    /* Declaration */
    int i,j, encore, k, l, dir, *vois, *rvois, *cvois, choix, xt, yt, cons, *reord, len;
    double **rr, **cc, *cs, ptir;
    
    /* Memory allocation */
    xt = (int) x[1];
    yt = (int) y[1];
    len = nr;
    if (nc < nr)
	len = nc;
    encore = 1;
    i = 1;
    j = 1;
    l = 0;
    k = 2;
    ptir = 0;
    dir = 1;
    choix = 1;
    cons = 0;
    
    taballoc(&rr, nr, nc);
    taballoc(&cc, nr, nc);
    vecintalloc(&vois, 4);
    vecintalloc(&reord, 4);
    vecintalloc(&rvois, 4);
    vecintalloc(&cvois, 4);
    vecalloc(&cs, 4);
    
    
    /* Rows and columns matrices */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    rr[i][j] = (double) i;
	    cc[i][j] = (double) j;
	}
    }
    
    
    cs[1] = pm[1];
    for (i = 2; i <= 4; i++) {
	cs[i] = cs[i-1] + pm[i];
    }
    
    /* Beginning of the loop */
    while (encore == 1) {
    
	/* Storage of the neighbouring */
	vois[1] = (int) map[xt][yt+1];
	vois[2] = (int) map[xt+1][yt];
	vois[3] = (int) map[xt][yt-1];
	vois[4] = (int) map[xt-1][yt];
	
	rvois[1] = (int) rr[xt][yt+1];
	rvois[2] = (int) rr[xt+1][yt];
	rvois[3] = (int) rr[xt][yt-1];
	rvois[4] = (int) rr[xt-1][yt];
	
	cvois[1] = (int) cc[xt][yt+1];
	cvois[2] = (int) cc[xt+1][yt];
	cvois[3] = (int) cc[xt][yt-1];
	cvois[4] = (int) cc[xt-1][yt];
	
	/* Re-order of the neighbouring according to the direction */
	l = 1;
	for (i = dir; i <= 4; i++) {
	    reord[l] = i;
	    l++;
	}
	i = 1;
	while (l != 4) {
	    reord[l] = i;
	    l++;
	    i++;
	}
	
	/* random draw of a direction */
	ptir = alea();
	choix = 4;
	if (ptir <= pm[1]) {
	    choix = 1;
	}
	if ((ptir > pm[1])&&(ptir <= pm[2])) {
	    choix = 2;
	}
	if ((ptir > pm[2])&&(ptir <= pm[3])) {
	    choix = 3;
	}
	
	/* And again, until the direction lead us into an accessible area */
	cons = reord[choix];
	
	while (vois[cons] == 1) {
	    ptir = alea();
	    choix = 4;
	    if (ptir <= pm[1]) {
		choix = 1;
	    }
	    if ((ptir > pm[1])&&(ptir <= pm[2])) {
		choix = 2;
	    }
	    if ((ptir > pm[2])&&(ptir <= pm[3])) {
		choix = 3;
	    }
	    cons = reord[choix];
	}
	
	/* Stores all information */
	xt = (int) (rvois[choix]);
	yt = (int) (cvois[choix]);
	
	x[k] = (double) xt;
	y[k] = (double) yt;
	
	if ((xt==1)|(xt==len)|(yt==1)|(yt==len))
	    encore = 0;
	if (k==nmax)
	    encore = 0;
	
	for (i = 1; i <= 4; i++) {
	    if ( ((int) rvois[i]) == xt) {
		if ( ((int) cvois[i]) == yt) {
		    dir = i;
		}
	    }
	}
	
	*nreel = k;
	k++;
    }
    
    /* free memory */
    freeintvec(vois);
    freeintvec(rvois);
    freeintvec(cvois);
    freeintvec(reord);
    freevec(cs);
    freetab(rr);
    freetab(cc);
}





/* For external call from within R */
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel)
{
    /* Declaration */
    double **map, *pm, *x, *y;
    int i, j, k, nr, nc, nmax;
    
    /* Memory allocation */
    nr = *nrm;
    nc = *ncm;
    nmax = *nmaxr;
    
    taballoc(&map, nr, nc);
    vecalloc(&x, nmax);
    vecalloc(&y, nmax);
    vecalloc(&pm, 4);
    
    /* R to C */
    x[1] = xr[0];
    y[1] = yr[0];
    
    k = 0;
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    map[i][j] = mapr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= 4; i++) {
	pm[i] = probamr[i-1];
    }
    
    /* Main function */
    perclu(map, nr, nc, x, y, nmax, nreel, pm);
    
    /* C to R */
    for (i = 1; i <= *nreel; i++) {
	xr[i-1] = x[i];
	yr[i-1] = y[i];
    }
    
    /* free memory */
    freevec(x);
    freevec(y);
    freevec(pm);
    freetab(map);
}



/* *********************************************************************
 *                                                                     *
 *         Rediscretization algorithm for a traject                    *
 *                                                                     *
 ***********************************************************************/



/* Resolves a quadratic equation of the type
 */

void resolpol(double a, double b, double c, double *x1, double *x2, int *warn)
{
    double delta;
    delta = (b * b) - 4 * a * c;
    *warn = 0;
    if (delta > 0) {
	*x1 = (-b - sqrt(delta)) / (2 * a);
	*x2 = (-b + sqrt(delta)) / (2 * a);
    } else {
	*warn = 1;
    }
}




void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff)
{
    /* Declaration */
    double R, xt, yt, a, b, c, pente, ori, x1, x2, di1, di2;
    int k, m, p, fini, ok, warn, *dedans, lo, new, pp;
    
    /* memory allocation */
    fini = 0;
    k = 1;
    p = 2;
    m = 1;
    ok = 0;
    a = 0;
    b = 0;
    c = 0;
    pente = 0;
    ori = 0;
    x1 = 0;
    x2 = 0;
    lo = 0;
    di1 = 0;
    di2 = 0;
    *neff = 0;
    new = 0;
    pp = 1;
    
    
    vecintalloc(&dedans,2);
    
    /* Main algorithm */
    while (fini == 0) {
	
	dedans[1] = 0;
	dedans[2] = 0;
	ok = 0;
	xt = xn[k];
	yt = yn[k];
	k++;
	new = 0;
	
	/* Determines the "upper" point */
	while (ok == 0) {
	    if (new == 1)
		p++;
	    R = sqrt((((x[p] - xt) * (x[p] - xt)) + ((y[p] - yt) * (y[p] - yt))));
	    if (R > u) {
		ok = 1;
	    } else {
		if (p == n) {
		    fini = 1;
		    ok = 1;
		}
	    }
	    new = 1;
	}
	m = p-1;
    
    if (fini == 0) {
	/* Does the difference between x[p] and x[m] = 0? */
	if ((fabs(x[p] - x[m]) > 0.000000000001)) {
	    /* Computes the slope between m and p */
	    pente = (y[p] - y[m]) / (x[p] - x[m]); /* when diff(x) == 0 ? */
	    /* The intercept */
	    ori = y[p] - (pente * x[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * xt) + (2 * pente * ori) - (2 * pente * yt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * yt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= x[m]) {
		if (x1 < x[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= x[p]) {
		if (x1 < x[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= x[m]) {
		if (x2 < x[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= x[p]) {
		if (x2 < x[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (x[p] - x1));
		di2 = fabs((double) (x[p] - x2));
		
		/* verify that xk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		xn[k] = x1;
		yn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		xn[k] = x2;
		yn[k] = (pente * x2) + ori;
	    }

	} else { /* We change x and y coordinates */
	    
	    /* Computes the slope between m and p */
	    pente =  (x[p] - x[m]) / (y[p] - y[m]);
	    /* The intercept */
	    ori = x[p] - (pente * y[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * yt) + (2 * pente * ori) - (2 * pente * xt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * xt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= y[m]) {
		if (x1 < y[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= y[p]) {
		if (x1 < y[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= y[m]) {
		if (x2 < y[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= y[p]) {
		if (x2 < y[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (y[p] - x1));
		di2 = fabs((double) (y[p] - x2));
		
		/* verify that yk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		yn[k] = x1;
		xn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		yn[k] = x2;
		xn[k] = (pente * x2) + ori;
	    }
	}
	
	/* Computes the nnew date (linear approximation) */
	di1 = sqrt((((xn[k] - x[m]) * (xn[k] - x[m])) + ((yn[k] - y[m]) * (yn[k] - y[m]))));
	R = sqrt((((x[p] - x[m]) * (x[p] - x[m])) + ((y[p] - y[m]) * (y[p] - y[m]))));
	di2 = dat[p] - dat[m];
	datn[k] = dat[m] + (di1 * di2 / R);
    }
    if (k == nn) {
	fini = 1;
    }
    pp = p;
    }
    
    /* Free memory */
    *neff = k;
    freeintvec(dedans);
}



/* For external Call from within R */

void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, int *neff)
{
    /* Declaration */
    int i, n, nn;
    double *x, *y, *xn, *yn, *dat, *datn, u;
    
    /* Memory allocation */
    n = *nr;
    nn = *nnr;
    u = *ur;
    
    vecalloc(&x, n);
    vecalloc(&y, n);
    vecalloc(&xn, nn);
    vecalloc(&yn, nn);
    vecalloc(&dat, n);
    vecalloc(&datn, nn);
    
    /* R to C */
    for (i = 1; i <= n; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	dat[i] = datr[i-1];
    }
    
    xn[1] = *xdeb;
    yn[1] = *ydeb;
    datn[1] = *dat0;
    
    /* Main function  */
    discretraj(x, y, dat, xn, yn, n, nn, datn, u, neff);
    
    /* C to R */
    for (i = 1; i <= nn; i++) {
	xnr[i-1] = xn[i];
	ynr[i-1] = yn[i];
	datnr[i-1] = datn[i];
    }
    
    /* Free memory */
    freevec(x);
    freevec(y);
    freevec(xn);
    freevec(yn);
    freevec(dat);
    freevec(datn);
}





/* *********************************************************************
 *                                                                     *
 *              Home range by Clustering (Kenward et al. 2001)         *
 *                                                                     *
 ***********************************************************************/


/* finds the cluster with the minimum average distance between the 3 points
   not assigned to a cluster */

void trouveclustmin(double **xy, int *clust, int *lo1, int *lo2,
		    int *lo3, double *dist)
{
    /* Declaration */
    int i, j, k, m, npas, nr, *indice;
    double **xy2, di1, di2, di3, ditmp;
    
    /* Memory allocation */
    nr = (int) xy[0][0];
    npas = 0;
    di1 = 0;
    di2 = 0;
    di3 = 0;
    ditmp = 0;
    
    /* Number of non assigned points */
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    npas++;
	}
    }
    taballoc(&xy2, npas, 2);
    vecintalloc(&indice, npas);
    
    /* The non assigned points are stored in xy2 */
    k = 1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    xy2[k][1] = xy[i][1];
	    xy2[k][2] = xy[i][2];
	    indice[k] = i;
	    k++;
	}
    }
    
    /* Computes the distane between the relocations */
    *dist = 0;
    m=0;
    for (i = 1; i <= (npas-2); i++) {
	for (j = (i+1); j <= (npas-1); j++) {
	    for (k = (j+1); k <= npas; k++) {
		di1 = sqrt((xy2[i][1] - xy2[j][1]) * (xy2[i][1] - xy2[j][1]) + 
			   (xy2[i][2] - xy2[j][2]) * (xy2[i][2] - xy2[j][2]) );
		di2 = sqrt((xy2[i][1] - xy2[k][1]) * (xy2[i][1] - xy2[k][1]) + 
			   (xy2[i][2] - xy2[k][2]) * (xy2[i][2] - xy2[k][2]));
		di3 = sqrt((xy2[k][1] - xy2[j][1]) * (xy2[k][1] - xy2[j][1]) + 
			   (xy2[k][2] - xy2[j][2]) * (xy2[k][2] - xy2[j][2]));
		/* average distance */
		ditmp = (di1 + di2 + di3) / 3;
		/* minimum distance */
		if ((m == 0) || (ditmp < *dist)) {
		    *dist = ditmp;
		    *lo1 = indice[i];
		    *lo2 = indice[j];
		    *lo3 = indice[k];
		}
		m = 1;
	    }
	}
    }
    /* free memory */
    freeintvec(indice);
    freetab(xy2);
}


/* For external call from within R */
void trouveclustminr(double *xyr, int *nr, int *clustr, int *lo1, int *lo2,
		     int *lo3, double *dist)
{
    /* Declaration */
    double **xy;
    int i, j, k, *clust;

    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    vecintalloc(&clust, *nr);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
    }
    for (i = 1; i <= *nr; i++) {
	clust[i] = clustr[i-1];
    }
    
    /* main function */
    trouveclustmin(xy, clust, lo1, lo2, lo3, dist);
    
    /* Free memory */
    freetab(xy);
    freeintvec(clust);
}


/* Finds the distance between a cluster of points and the nearest point*/
void nndistclust(double **xy, double *xyp, double *dist)
{
    /* Declaration */
    int n, i, m;
    double di;

    m = 0;
    di =0;
    n = (int) xy[0][0];
    *dist = 0;
    
    /* finds the distance and the corresponding point */
    for (i = 1; i <= n; i++) {
	di = sqrt( (xy[i][1] - xyp[1]) * (xy[i][1] - xyp[1]) + 
		   (xy[i][2] - xyp[2]) * (xy[i][2] - xyp[2]) );
	if ( (di < *dist) || (m == 0) ) {
	    *dist = di;
	}
	m = 1;
    }
}


/* The function nndistclust is applied for all available clusters */
void parclust(double **xy, int *clust, int *noclust, 
	      int *noloc, double *dist)
{
    /* Declaration */
    int i, k, m, nr2, nr, nocl;
    double **xy2, *xyp, di;
    
    /* Memory allocation */
    nocl = *noclust;
    nr = xy[0][0];
    nr2 = 0;
    
    /* The number of available clusters */
    for (i = 1; i <= nr; i++) {
	if (clust[i] == nocl) {
	    nr2++;
	}
    }

    taballoc(&xy2, nr2, 2);
    vecalloc(&xyp, 2);
    
    /* stores the non assigned points in xy2 */
    k = 1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == nocl) {
	    xy2[k][1] = xy[i][1];
	    xy2[k][2] = xy[i][2];
	    k++;
	}
    }
    
    /* Finds the minimum distance between a point and a cluster, 
       performed for all clusters */
    di = 0;

    m = 0;
    *dist = 0;
    for (i = 1; i <= nr; i++) {
	if (clust[i] != nocl) {
	    xyp[1] = xy[i][1];
	    xyp[2] = xy[i][2];
	    nndistclust(xy2, xyp, &di);
	    if ( (di < *dist) || (m == 0) ) {
		*dist = di;
		*noloc = i;
	    }
	    m = 1;
	}
    }

    /* Free memory */
    freetab(xy2);
    freevec(xyp);
}


/* The function trouveminclust identifies the cluster for which the nearest 
   point is the closest */
void trouveminclust(double **xy, int *liclust, int *clust, 
		    int *noclust, int *noloc, double *dist)
{
    /* Declaration */
    int i, nr, nc, m, labclust, nolo;
    double di;
    
    nr = (int) xy[0][0];
    nc = 0;
    di = 0;
    labclust = 0;
    nolo = 0;
    
    /* Assigned clusters */
    for (i = 1; i <= nr; i++) {
	if (liclust[i] > 0) {
	    nc++;
	}
    }
    
    /* finds the minimum distance between a cluster and its nearest point 
       (the cluster name and the point ID are searched) */
    m = 0;
    *dist = 0;
    for (i = 1; i <= nc; i++) {
	labclust = liclust[i];
	parclust(xy, clust, &labclust, &nolo, &di);
	if ( (m == 0) || (di < *dist) ) {
	    *dist = di;
	    *noloc = nolo;
	    *noclust = labclust;
	}
	m = 1;
    }
}


/* What should be done: create a new cluster or add a relocation 
   to an existing one ? */

void choisnvclust(double **xy, int *liclust, int *clust, int *ordre)
{
    /* Declaration */
    int i, k, nr, noloat, cluat, nolo1, nolo2, nolo3, maxclust;
    int maxindiceclust, clu1, *liclub, nz;
    double dmoyclust, dminloc;
    
    /* Memory allocation */
    nz = 0;
    i = 0;
    k = 0;
    nr = (int) xy[0][0];
    maxclust = 0;
    maxindiceclust = 0;
    nolo1 = 0;
    nolo2 = 0;
    nolo3 = 0;
    noloat = 0;
    cluat = 0;
    clu1 = 0;
    vecintalloc(&liclub, nr);
    
    /* finds the max label for the cluster */
    for (i = 1; i <= nr; i++) {
	if (clust[i] != 0) {
	    if (clust[i] > maxclust) {
		maxclust = clust[i];
	    }
	    if (liclust[i] != 0) {
		maxindiceclust = i;
	    }
	}
    }
    
    /* Finds the min distance between 3 relocations */
    trouveminclust(xy, liclust, clust, &cluat, &noloat, &dminloc);
    
    /* Computes the average distance between the locs of the smaller cluster */
    /* First, one verifies that there is at least Three non assigned locs */
    dmoyclust = dminloc +1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    nz++;
	}
    }
    if (nz > 3) {
	dmoyclust = 0;
	trouveclustmin(xy, clust, &nolo1, &nolo2, &nolo3, &dmoyclust);
    }
    
    /* First case: A new cluster independent from the others */
    if (dmoyclust < dminloc) {
	ordre[nolo1] = 1;
	ordre[nolo2] = 1;
	ordre[nolo3] = 1;
	
	clust[nolo1] = maxclust + 1;
	clust[nolo2] = maxclust + 1;
	clust[nolo3] = maxclust + 1;
	
	liclust[maxindiceclust+1] = maxclust + 1;
	
    } else {
	/* Second case: one loc is added to a cluster */
	
	/* Case 2.1: the loc does not belong to one cluster */
	if (clust[noloat] == 0) {
	    ordre[noloat] = 1;
	    clust[noloat] = cluat;
	} else {
	    
	    /* Case 2.2: the loc belong to one cluster: fusion */
	    clu1 = clust[noloat];
	    for (i = 1; i <= nr; i++) {
		if (clust[i] == clu1) {
		    clust[i] = cluat;
		    ordre[i] = 1;
		}
		if (liclust[i] == clu1) {
		    liclust[i] = 0;
		}
	    }
	    /* and cleaning of liclust */
	    k = 1;
	    for (i = 1; i <= nr; i++) {
		if (liclust[i] != 0) {
		    liclub[k] = liclust[i];
		    k++;
		}
	    }
	    for (i = 1; i <= nr; i++) {
		liclust[i] = liclub[i];
	    }
	}
    }
    freeintvec(liclub);
}


/* The main function for home range computation */

void clusterhr(double **xy, int *facso, int *nolocso, int *cluso)
{
    /* Declaration */
    int i, nr, lo1, lo2, lo3, *clust, len, con, *ordre, *liclust, courant;
    double di;

    /* Memory allocation */
    courant = 1;
    nr = (int) xy[0][0];
    vecintalloc(&clust, nr);
    vecintalloc(&ordre, nr);
    vecintalloc(&liclust, nr);
    lo1 = 0;
    lo2 = 0;
    lo3 = 0;
    di = 0;
    con = 1;
    len = 0;
    
    /* Begin: Search for the first cluster */
    trouveclustmin(xy, clust, &lo1, &lo2,
		   &lo3, &di);
    
    clust[lo1] = 1;
    clust[lo2] = 1;
    clust[lo3] = 1;
    liclust[1] = 1;
    len = 3;
    
    /* We store it in the output */
    cluso[1] = 1;
    cluso[2] = 1;
    cluso[3] = 1;
    nolocso[1] = lo1;
    nolocso[2] = lo2;
    nolocso[3] = lo3;
    facso[1] = 1;
    facso[2] = 1;
    facso[3] = 1;
    
    /* Then repeat until all relocations belong to the same cluster */
    while (con == 1) {
	courant++;
	
	for (i = 1; i <= nr; i++) {
	    ordre[i] = 0;
	}
	
	choisnvclust(xy, liclust, clust, ordre);
	
	for (i = 1; i <= nr; i++) {
	    if (ordre[i] != 0) {
		len++;
		cluso[len] = clust[i];
		nolocso[len] = i;
		facso[len] = courant;
	    }
	}
	
	con = 0;
	for (i = 2; i <= nr; i++) {
	    if (clust[i] != clust[1])
		con = 1;
	}
	if (con == 0) {
	    con = 0;
	}
    }
    
    /* Free memory */
    freeintvec(clust);
    freeintvec(ordre);
    freeintvec(liclust);
}



/* Finds the length of the output for the table containing the home range */

void longfacclust(double **xy, int *len2)
{
    /* Declaration */
    int i, nr, lo1, lo2, lo3, *clust, len, con, *ordre, *liclust, courant;
    double di;
    
    /* Memory allocation */
    courant = 1;
    nr = (int) xy[0][0];
    vecintalloc(&clust, nr);
    vecintalloc(&ordre, nr);
    vecintalloc(&liclust, nr);
    lo1 = 0;
    lo2 = 0;
    lo3 = 0;
    di = 0;
    con = 1;
    len = 0;
    
    /* Begin: search for the first cluster */
    trouveclustmin(xy, clust, &lo1, &lo2,
		   &lo3, &di);
    clust[lo1] = 1;
    clust[lo2] = 1;
    clust[lo3] = 1;
    liclust[1] = 1;
    len = 3;
    
    /* Counts the number of rows needed for the table, which will contain the results */
    while (con == 1) {
	courant++;
	
	for (i = 1; i <= nr; i++) {
	    ordre[i] = 0;
	}
	
	choisnvclust(xy, liclust, clust, ordre);
	
	for (i = 1; i <= nr; i++) {
	    if (ordre[i] != 0) {
		len++;
	    }
	}
	con = 0;
	for (i = 2; i <= nr; i++) {
	    if (clust[i] != clust[1])
		con = 1;
	}
	if (con == 0) {
	    con = 0;
	}
    }

    *len2 = len;

    /* Free memory */
    freeintvec(clust);
    freeintvec(ordre);
    freeintvec(liclust);
}


/* For external call from within R */
void longfacclustr(double *xyr, int *nr, int *len2)
{
    /* Declaration */
    double **xy;
    int i, j, k;
    
    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}    
    }
    
    /* Main function */
    longfacclust(xy, len2);

    /* Free memory */
    freetab(xy);
}



/* For external call of clusterhrr from within R */

void clusterhrr(double *xyr, int *nr, int *facsor, 
		int *nolocsor, int *clusor, int *len)
{
    /* Declaration */
    double **xy;
    int i, j, k, *facso, *nolocso, *cluso;
    
    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    vecintalloc(&facso, *len);
    vecintalloc(&nolocso, *len);
    vecintalloc(&cluso, *len);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
    }
    
    /* Main function */
    clusterhr(xy, facso, nolocso, cluso);
    
    /* C to R */
    for (i = 1; i <= *len; i++) {
	facsor[i-1] = facso[i];
	nolocsor[i-1] = nolocso[i];
	clusor[i-1] = cluso[i];
    }
    
    /* Free memory */
    freetab(xy);
    freeintvec(facso);
    freeintvec(nolocso);
    freeintvec(cluso);
}



void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr)
{
  double **xy, **R2n, *xp, *dt, **dtsim, tt;
  int n, i, j, k, *index, nr;
  
  n = *nro;
  nr = *nrepr;
  tt = 0;
  vecalloc(&xp, 2);
  vecalloc(&dt, n);
  taballoc(&R2n, n, nr);
  taballoc(&dtsim, n, nr);
  vecintalloc(&index, n);
  taballoc(&xy, n, 2);
  
  k = 0;
  for (i = 1; i <= n; i++) {
    dt[i] = dtr[i-1];
    for (j = 1; j <= 2; j++) {
      xy[i][j] = xyr[k];
      k++;
    }
  }
  
  for (k = 1; k <= nr; k++) {
    
    getpermutation(index, k);
    j = index[1];
    xp[1] = 0;
    xp[2] = 0;
    tt = 0;

    for (i = 1; i <= n; i++) {
      j = index[i];
      xp[1] = xp[1] + xy[j][1];
      xp[2] = xp[2] + xy[j][2];
      R2n[i][k] = (xp[2] * xp[2]) + (xp[1] * xp[1]);
      tt = tt + dt[j];
      dtsim[i][k] = tt;
    }
  }
  
  
  k = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j<= nr; j++) {
      R2nr[k] = R2n[i][j];
      dtsimr[k] = dtsim[i][j];
      k++;
    }
  }
  
  
  freevec(xp);
  freevec(dt);
  freeintvec(index);
  freetab(xy);
  freetab(dtsim);
  freetab(R2n);
}




void runsltr(int *xr, int *nr, double *res, int *nrepr)
{
  int i, j, n, *x, *xb, nrep, nbsui, nz, nu, *numero;
  double m, s, n1, n2;
  
  n = *nr;
  nrep = *nrepr;
  
  vecintalloc(&x, n);
  vecintalloc(&xb, n);
  vecintalloc(&numero, n);

  
  for (i = 1; i <= n; i++) {
    x[i] = xr[i-1];
    numero[i] = i;
  }
  
  nbsui = 1;
  nz = 0;
  nu = 1;
  if (x[1] == 0)
    nz = 1;
  if (x[1] == 1)
    nu = 1;
  
  for (i = 2; i <= n; i++) {
    if (x[i] == 0)
      nz++;
    if (x[i] == 1)
      nu++;
    if (x[i-1] != x[i])
      nbsui++;
  }
  
  n1 = ((double) nz);
  n2 = ((double) nu);
  
  m = 1 + 2 * n1 * n2 / (n1 + n2);
  s = sqrt(2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)/((n1 + n2) * (n1 + n2) * (n1 + n2 - 1)));
  
  res[0] = (((double) nbsui) - m) / s;
  
  
  for (j = 1; j <= nrep; j++) {
    nbsui = 1;
    getpermutation(numero, j);
    trirapideint(numero , x, 1, n);
    for (i = 2; i <= n; i++) {
      if (x[i-1] != x[i])
	nbsui++;
    }
    res[j] = (((double) nbsui) - m) / s;
  }
  
  freeintvec(x);
  freeintvec(xb);
  freeintvec(numero);

}




void testindepangl (double *sim, double *ang, int *nang, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *angle;
  vecalloc(&angle, *nang);
  for (i=1; i<=*nang; i++) {
    angle[i] = ang[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=debut[k]; i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+(1.0-cos(angle[i+1] - angle[i]));
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(angle);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=debut[k]; i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+(1.0-cos(angle[i+1] - angle[i]));
      }
    }
  }
  for (j=0; j<=*ni; j++) {
    sim[j]=2*sim[j];

  }
  freevec(angle);
}

void testindepdist (double *sim, double *di, int *ndi, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *dist;
  vecalloc(&dist, *ndi);
  for (i=1; i<=*ndi; i++) {
    dist[i] = di[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=(debut[k]); i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+pow(dist[i+1] - dist[i],2);
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(dist);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=(debut[k]); i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+pow(dist[i+1] - dist[i],2);
      }
    }

  }
  freevec(dist);
}


void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur)
{
  int i,j,k, nu, nr, nc, *rt;
  double *dtu, **dtr, **res, **ru, tmp1;
  
  nu = *nur;
  nr = *nrr;
  nc = *ncr;
  
  vecalloc(&dtu, nu);
  vecintalloc(&rt, nc);
  taballoc(&dtr, nr, nc);
  taballoc(&res, nr, nc);
  taballoc(&ru, nu, nc);
  
  for (i = 1; i <= nu; i++) {
    dtu[i] = dtur[i-1];
  }
  
  k = 0;
  for (i = 1; i <= nr; i++) {
    for (j = 1; j<= nc; j++) {
      dtr[i][j] = dtrr[k];
      res[i][j] = resr[k];
      k++;
    }
  }
  
  for (i = 1; i <= nu; i++) {
    for (k = 1; k <= nc; k++) {
      rt[k] = 0;
    }
    tmp1 = dtu[i];
    for (j = 1; j <= nr; j++) {
      for (k = 1; k <= nc; k++) {
	if ((fabs((double) (dtr[j][k] - tmp1)))< 0.0000000001)
	  rt[k] = j;
      }
    }
    for (k = 1; k <= nc; k++) {
      if ((fabs((double) rt[k] )) < 0.00000000001) {
	ru[i][k] = -1;
      } else {
	j = rt[k];
	ru[i][k] = res[j][k];
      }
    }
  }
  
  
  k = 0;
  for (i = 1; i <= nu; i++) {
    for (j = 1; j<= nc; j++) {
      rur[k] = ru[i][j];
      k++;
    }
  }
  
  freevec(dtu);
  freeintvec(rt);
  freetab(dtr);
  freetab(res);
  freetab(ru);
}





/* *********************************************************************
 *                                                                     *
 *                   partition d'un trajet                             *
 *                                                                     *
 ***********************************************************************/

void optcut (double **Pid, double **mk, int *maxk) 
{
    /* Declaration of variables */
    int l, p, i, ii, Km, k, j;
    double **mkd, **mkdn, tmp, mi1, mi2, mi3;
    
    
    
    /* memory allocation */
    l = Pid[0][0];   /* The length of the sequence */
    p = Pid[1][0];   /* The number of models */
    Km = *maxk;     /* the max number of partition to compute */
    i = 0;     /* the sequence index used in the paper (from 0 to l-1) */
    ii = 0;    /* the position index used in the program (from 1 to l) */
    j = 0;    /* the position index for the model (from 1 to p) */
    k = 0;    /* the index for the partition */
    taballoc(&mkd, l, p); /* The table for the log-probabilities of a 
			     partition given a model 
			     for a k partition
			     (to be re-used for each k) */
    taballoc(&mkdn, l, p); /* The table for the log-probabilities 
			      of a partition given a model
			     for a k+1 partition
			     (to be re-used for each k) */
    mi1 = 0;
    mi2 = 0;
    mi3 = 0;
    
    
    /* 1. First computes the probability of one partition*/
    
    /* 1.1 Fills the table mkd   */
    
    for (j = 1; j <= p; j++) {
	mkd[1][j] = log(Pid[1][j]);
    }
    for (j = 1; j <= p; j++) {
	for (ii = 2; ii <= l; ii++) {
	    mkd[ii][j] = mkd[ii-1][j] + log(Pid[ii][j]);
	}
    }
    
    /* 1.2 Fills the table mk */
    for (ii = 1; ii <= l; ii++) {
	
        /* computes the max */
	tmp = mkd[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkd[ii][j] - tmp > 0.0000000000001) {
		tmp = mkd[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][1] = exp(mkd[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][1] = mk[ii][1] + (exp(mkd[ii][j]) / ( (double) p ) );
	}
	
        /* takes the log and adds the max again */
	mk[ii][1] = log(mk[ii][1]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] + tmp;
	}
    }
    
    
    
    
    /* 2. computes mk for each possible r-partitions: a loop */
    for (k = 2; k <= Km; k++) {
	
	
        /* First start at the first step */
	i = k-1;
	ii = k;
	
	/* Computes mkdn for each d and ii */
	for (j = 1; j <= p; j++) {
	    
	    /* Removes the max */
	    tmp = mk[ii-1][k-1];
	    if (mk[ii-1][k-1] - mkd[ii-1][j] < -0.00000000001)
		tmp = mkd[ii-1][j];
	    mi1 = mk[ii-1][k-1] - tmp;
	    mi2 = mkd[ii-1][j] - tmp;
	    
	    /* Computes for i = k-1 */
	    mkdn[ii][j] = log(Pid[ii][j]) + 
		log((((double) k) - 1) / (((double) i) * (((double) p) - 1))) +
		log(((double) p) * exp(mi1) - exp(mi2));
	    
	    /* adds the max again */
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	


	/* fills mk */
	/* computes the max */
	tmp = mkdn[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkdn[ii][j] - tmp > 0.00000000000000001) {
		tmp = mkdn[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	}
	
	/* takes the log and adds the max again */
	mk[ii][k] = log(mk[ii][k]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	
	
	/* Then increase the sequence */
	for (ii = (k + 1) ; ii <= l; ii++) {
	    
	    i = ii - 1;
	    
	    /* again computes mkdn for each d */
	    for (j = 1; j <= p; j++) {
		
		tmp = mkdn[ii-1][j];
		if (tmp - mk[ii-1][k-1] < -0.00000000000001) {
		    tmp = mk[ii-1][k-1];
		}
		if (tmp - mkd[ii-1][j] < -0.000000000000000001) {
		    tmp = mkd[ii-1][j];
		} 
		mi1 = mkdn[ii-1][j] - tmp;
		mi2 = mk[ii-1][k-1] - tmp;
		mi3 = mkd[ii-1][j] - tmp;
		
		mkdn[ii][j] = log(Pid[ii][j]) +
		    log((( ((double) (i - k + 1)) / ((double) i)) * 
			exp(mi1)) + (((((double) k) - 1) / 
				    (((double) i) * 
				     (((double) p) - 1))) * 
			((((double) p) * exp(mi2)) - 
			 exp(mi3)))) + tmp;
	    }


	    
	    /* ... and again fills mk */
	    /* computes the max */
	    tmp = mkdn[ii][1];
	    for (j = 1; j <= p; j++) {
		if (mkdn[ii][j] - tmp > 0.000000000000001) {
		    tmp = mkdn[ii][j];
		}
	    }
	    
	    /* removes the max */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] - tmp;
	    }
	    
	    /* computes the mean */
	    mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	    for (j = 2; j <= p; j++) {
		mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	    }
	    
	    /* takes the log and adds the max again */
	    mk[ii][k] = log(mk[ii][k]) + tmp;
	    
	    /* adds the max again for mkd */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] + tmp;
	    }
	    
	}
	
	/* and finally, replace mkd with mkdn for next k */    
	for (ii = 1; ii <= l; ii++) {
	    for (j = 1; j <= p; j++) {
		mkd[ii][j] = mkdn[ii][j];
	    }
	}
    }
    
    
    /* free memory */
    freetab(mkd);
    freetab(mkdn);

}



/* Now, the R version */

void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye) 
{
    /* Declaration */
    int i, j, k, l, p, Km;
    double **Pid, **mk, tmp, *mik, *msk, *kk;
    
    /* Memory allocation */
    l = *lr;
    Km = *maxk;
    p = *pr;
    taballoc(&Pid, l, p);
    taballoc(&mk, l, Km);
    vecalloc(&mik, Km);
    vecalloc(&msk, Km);
    vecalloc(&kk, Km);

    
    /* Fills local variables */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= p; j++) {
	    Pid[i][j] = Pidr[k];
	    k++;
	}
    }
    
    optcut(Pid, mk, maxk);
    
    /* keeps the last row */
    for (k = 1; k <= Km; k++) {
	kk[k] = mk[l][k];
    }
    /* computes the max */
    tmp = kk[1];
    for (k = 1; k <= Km; k++) {
	if (tmp - kk[k] < -0.0000000001)
	    tmp = kk[k];
    }
    
    /* removes the max */
    for (k = 1; k <= Km; k++) {
	kk[k] = kk[k] - tmp;
    }
    
    /* Computes mik */
    mik[1] = exp(kk[1]);
    for (k = 2; k <= Km; k++) {
	mik[k] = mik[k-1] + (exp(kk[k]) / ((double) k ));
    }
    
    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	mik[k] = log(mik[k]) + tmp;
    }
    
    /* Computes msk */
    msk[Km] = exp(kk[Km]) / ((double) (l - Km + 1)) ;
    for (k = Km-1; k >= 1; k--) {
	msk[k] = msk[k+1] + exp(kk[k]) / ((double) (l - k + 1)  );
    }

    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	msk[k] = log(msk[k]) + tmp;
    }
    
    
    /* The bayes factor */
    for (k = 2; k <= (Km - 1); k++) {
	baye[k-2] = msk[k] - mik[k-1];
    }

    /* ... and back to R */
    k = 0;
    for (j = 1; j<=Km; j++) {
	mkr[k] = mk[l][j];
	k++;
    }
    


    /* Free memory */
    freevec(msk);
    freevec(mik);
    freevec(kk);
    freetab(Pid);
    freetab(mk);
}



/* *********************************************************************
 *                                                                     *
 *                   partitioning of a trajectory                      *
 *                                                                     *
 ***********************************************************************/

void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res)
{
    /* declaration of variables */
    int i, j, k, l, D, Km;
    double **Mkk, **cumPid, tmp;
    
    /* Memory allocation */
    l = Pid[0][0]; /* length of the sequence */
    D = Pid[1][0]; /* Number of models */
    Km = *maxk;    /* Partition size */
    tmp = 0;

    
    taballoc(&Mkk, Km, D); /* Contains mkd for i = k, for all models */
    taballoc(&cumPid, l, D); /* For the probability of the sequences */
    
    
    /* First Compute the prediction of the sequences */
    for (j = 1; j <= D; j++) {
	cumPid[1][j] = Pid[1][j];
    }
    for (i = 2; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    cumPid[i][j] = cumPid[i-1][j] + Pid[i][j];
	    /* fills res */
	    res[i][j] = cumPid[i][j];
	}
    }
    
    
    /* Computes M1 */    
    for (i = 1; i <= l; i++) {
	Mk[i][1] = cumPid[i][1];
	for (j = 2; j <= D; j++) {
	    if (Mk[i][1] < cumPid[i][j]) {
		Mk[i][1] = cumPid[i][j];
	    }
	}
    }
    
    /* Then, computes Mkk for all k <= Km */
    for (j = 1; j <= D; j++) {
	Mkk[1][j] = Pid[1][j];
    }
    
    for (k = 2; k <= Km; k++) {
	/* Computes Mkk */
	for (j = 1; j <= D; j++) {
	    Mkk[k][j] = Pid[k][j] + Mk[k-1][k-1];
	}
	/* Update Mk */
	Mk[k][k] = Mkk[k][1];
	for (j = 2; j <= D; j++) {
	    if (Mkk[k][j] - Mk[k][k] < 0.000000000001)
		Mk[k][k] = Mkk[k][j];
	}
    }
    
    
    /* Now, compute Mkd */
    for (k = 2; k <= Km; k++) {
	
	/* Deletes for i < k */
	for (i = 1; i < k; i++) {
	    for (j = 1; j <= D; j++) {
		Mkd[i][j] = 0;
	    }
	}
		
	/* first fill the first line of Mkd */
	for (j = 1; j <= D; j++) {
	    Mkd[k][j] = Mkk[k][j];
	}

	/* Computes Mkd */
	for (i = (k+1); i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		tmp = Mk[i-1][k-1];
		if (Mkd[i-1][j] - tmp > 0.000000000000000001) {
		    tmp = Mkd[i-1][j];
		}
		Mkd[i][j] = Pid[i][j] + tmp;
	    }
	    
	    /* Update Mk */
	    Mk[i][k] = Mkd[i][1];
	    for (j = 1; j <= D; j++) {
		if (Mkd[i][j] - Mk[i][k] > 0.000000000000000001) {
		    Mk[i][k] = Mkd[i][j];
		}
	    }
	    
	}
	
	/* Fills res */
	for (i = 1; i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		res[((k-1) * l)+i][j] = Mkd[i][j];
	    }
	}
    }

    /* free memory */
    freetab(Mkk);
    freetab(cumPid);
}



void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr)
{
    /* Variable declaration */
    int l, D, Km, i, j, k, m, n, *curloc, *curmod;
    double **Mk, **Mkd, **Pid, **res, *curma, **grap;
    
    /* Memory allocation */
    l = *lr;
    D = *Dr;
    Km = *Kmr;
    m = 0;
    n = 0;
    taballoc(&Mk, l, Km);
    taballoc(&Mkd, l, D);
    taballoc(&Pid, l, D);
    taballoc(&res, (l * Km), D);
    taballoc(&grap, l, D);
    vecalloc(&curma, Km);
    vecintalloc(&curmod, Km);
    vecintalloc(&curloc, (Km+1));


    
    
    /* R to C */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    Pid[i][j] = log(Pidr[k]);
	    k++;
	}
    }
    
    /* The main algorithm */
    partraj(Pid, Kmr, Mk, Mkd, res);
    
    
    /* Backtracking */
    curloc[1] = l;
    for (k = Km; k>=1; k--) {
	if (k > 1) {
	    
	    m = Km - k + 2;
	    
            /* Stores the graph */
	    for (i = 1; i <= l; i++) {
		for (j = 1; j <= D; j++) {
		    if (res[(l * (k-1)) + i][j] > Mk[i][k-1]) {
			grap[i][j] = 1;
		    } else {
			grap[i][j] = 0;
		    }
		}
	    }
	    
	    /* The best model for the last step */
	    n = (l * (k-1)) + curloc[m-1];
	    curma[m-1] = res[n][1];
	    curmod[m-1] = 1;
	    for (j = 2; j <= D; j++) {
		if (res[n][j] > curma[m-1]) {
		    curma[m-1] = res[n][j];
		    curmod[m-1] = j;
		}
	    }
	    
	    /* Keep this model until ? */
	    n = curloc[m-1];
	    j = curmod[m-1];
	    while (grap[n][j] > 0.0000000001)
		n--;
	    curloc[m] = n;
	} else {
	    m = Km+1;
	    curloc[m] = 1;
	    curma[m-1] = res[curloc[m-1]][1];
	    curmod[m-1] = 1;
	    for (j = 1; j <= D; j++) {
		if (res[curloc[m-1]][j] > curma[m-1]) {
		    curma[m-1] = res[curloc[m-1]][j];
		    curmod[m-1] = j;
		}
	    }
	}
	
	
    }
    
    
    
    /* C to R */
    for (i = 1; i <= Km; i++) {
	curmodr[i-1] = curmod[i];
	curlocr[i-1] = curloc[i];
	curmar[i-1] = curma[i];
    }
    curlocr[Km] = curloc[Km+1];
    
    
    /* free memory */
    freetab(Mk);
    freetab(Mkd);
    freetab(res);
    freetab(grap);
    freevec(curma);
    freeintvec(curmod);
    freeintvec(curloc);
}




/* *********************************************************************
 *                                                                     *
 *           Kernel in time and space (Keating and Cherry, 2008)       *
 *                                                                     *
 ***********************************************************************/



void kcprcirc(double **xyd, double *h, double *x, double t, 
	      double *val)
{
    int i, j, nlo;
    double tmp, tmp2, vi, som;
    
    nlo = xyd[0][0];
    som=0;
    
    for (i=1; i<=nlo; i++) {
	
	tmp2 = 1;
	
	/* spatial coordinates */
	for (j=1; j<=2; j++) {
	    vi = (x[j] - xyd[i][j]) / h[j];
	    if (fabs(vi) < 1.0) {
		tmp = (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
		tmp2 = tmp2 * tmp;
	    } else {
		tmp2 = 0.0;
	    }
	}
	
	/* time */
	vi = (t - xyd[i][3]);
	tmp2 = tmp2 * (1-h[3]) * (1 - h[3] * h[3]) / 
	    (2 * 3.14153265359 * ( 1 +  (h[3] * h[3]) - (2 * h[3] * cos(vi) )));
	som = som + tmp2;
    }
    *val = ( 1/( ((double) nlo) * h[1] * h[2] * h[3])) * som;
}

void kcprlin(double **xyd, double *h, double *x, double t, 
	      double *val)
{
    int i, j, nlo;
    double tmp, tmp2, vi, som;
    
    nlo = xyd[0][0];
    
    som = 0;
    
    for (i=1; i<=nlo; i++) {
	
	tmp2 = 1;
	
	/* spatial coordinates and time */
	for (j=1; j<=2; j++) {
	    vi = (x[j] - xyd[i][j]) / h[j];
	    if (fabs(vi) < 1.0) {
		tmp = (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
		tmp2 = tmp2 * tmp;
	    } else {
		tmp2 = 0.0;
	    }
	}
	
	/* time */
	vi = (t - xyd[i][3]) / h[3];
	if (fabs(vi) < 1.0) {
	    tmp2 = tmp2 * (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
	} else {
	    tmp2 = 0.0;
	}
	som = som + tmp2;
    }
    *val = ( 1/( ((double) nlo) * h[1] * h[2] * h[3])) * som;
}


void kernelkcr(double *xydr, double *tcalcr, int *nlr, double *gridr,
	       double *xgri, double *ygri, int *nliggri, int *ncolgri, 
	       double *hr, int *circularr)
{
    /* Declaration of local variables */
    int i, j, k, ncg, nlg, nlo, circular;
    double **gri, **xyd, *xx, *xg, *yg, tmp, *h, tca;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nlr;
    tca = *tcalcr;
    circular = *circularr;
    
    taballoc(&gri, nlg, ncg);
    taballoc(&xyd, nlo, 3);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&h, 3);
    vecalloc(&xx, 2);
    
    /* R objects -> C objects */

    for (i=1; i<=3; i++) {
	h[i] = hr[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    k = 0;
    for (i=1; i<=nlo; i++) {
	for (j=1; j<=3; j++) {
	    xyd[i][j] = xydr[k];
	    k++;
	}
    }
    
    /* loop on the grid */
    if (circular==1) {
	for (i=1; i<=nlg; i++) {
	    for (j=1; j<=ncg; j++) {
		xx[1] = xg[i];
		xx[2] = yg[j];
		kcprcirc(xyd, h, xx, tca, &tmp);
		gri[i][j] = tmp;
	    }
	}
    } else {
	for (i=1; i<=nlg; i++) {
	    for (j=1; j<=ncg; j++) {
		xx[1] = xg[i];
		xx[2] = yg[j];
		kcprlin(xyd, h, xx, tca, &tmp);
		gri[i][j] = tmp;
	    }
	}
    }
    
    /* C objects -> R objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    gridr[k] = gri[i][j];
	    k++;
	}
    }

    /* Memory Free */
    freetab(gri);
    freetab(xyd);
    freevec(xg);
    freevec(yg);
    freevec(h);
    freevec(xx);
}




/* *********************************************************************
 *                                                                     *
 *                      find local maxima/minima on a map              *
 *                                                                     *
 ***********************************************************************/


void findmaxgrid(double *grille, int *nlig, int *ncol)
{
    /* declaration */
    int i,j,k,nl,nc,sto;
    double **x, **grille2, r1,r2,r3,r4,r5,r6,r7,r8;
    
    nl = *nlig;
    nc = *ncol;
    sto=0;
    
    /* Memory alloocation */
    taballoc(&x,nl,nc);
    taballoc(&grille2,nl,nc);
    
    /* R to C */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    /* Loop */
    for (i = 2; i <= (nl-1); i++) {
	for (j=2; j<= (nc-1); j++) {
	    
	    r1=x[i-1][j-1] - x[i][j];
	    r2=x[i][j-1] - x[i][j];
	    r3=x[i+1][j-1] - x[i][j];
	    r4=x[i+1][j] - x[i][j];
	    r5=x[i+1][j+1] - x[i][j];
	    r6=x[i][j+1] - x[i][j];
	    r7=x[i-1][j+1] - x[i][j];
	    r8=x[i-1][j] - x[i][j];
	    
	    sto=1;
	    
	    if (r1 < -0.000000000001) {
		if (r2 < -0.000000000001) {
		    if (r3 < -0.000000000001) {
			if (r4 < -0.000000000001) {
			    if (r5 < -0.000000000001) {
				if (r6 < -0.000000000001) {
				    if (r7 < -0.000000000001) {
					if (r8 < -0.000000000001) {
					    sto=0;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	    
	    if (sto==0)
		grille2[i][j] = 1;
	}
    }
        
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    grille[k]=grille2[i][j];
	    k++;
	}
    }
    
    /* Memory */
    freetab(x);
    freetab(grille2);
    
}













/* *********************************************************************
 *                                                                     *
 *                 Method of Engen (2008)                              *
 *                                                                     *
 ***********************************************************************/


/* 
   sources of the package deal-1.2-30. 
   Code from Susanne Gammelgaard Bøttcher <alma@math.aau.dk>, 
   Claus Dethlefsen <cld@rn.dk>. 
   Useful for matrix inverse!
*/


int *ivector(int nl, int nh)
{
   int *v;

   v=(int *) R_alloc((unsigned) (nh-nl+1)*sizeof(int),sizeof(int));
   if ( v == NULL ){
      error("memory allocation failure in ivector()"); return(NULL);
   }
   return v-nl;
}

void free_ivector(int *v, int nl, int nh) { free((char*) (v+nl)); }



int invers(double **a, int n, double **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol=1,irow=1,j,k,l,ll;
   double big,dum,pivinv;

   if( (indxc = ivector(1,n)) == NULL){ return(-1); }
   if( (indxr = ivector(1,n)) == NULL){ return(-1); }
   if( (ipiv  = ivector(1,n)) == NULL){ return(-1); }
   for (j=1;j<=n;j++) ipiv[j]=0;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if (ipiv[j] != 1)
            for (k=1;k<=n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(a[j][k]) >= big) {
                     big=fabs(a[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1){
                  error("Invers: Singular Matrix-1");
                  return(-1);
               }
            }
      ++(ipiv[icol]);
      if (irow != icol) {
         for (l=1;l<=n;l++){
            double temp=a[irow][l]; a[irow][l]=a[icol][l]; a[icol][l]=temp;
         }
         for (l=1;l<=m;l++){
            double temp=b[irow][l]; b[irow][l]=b[icol][l]; b[icol][l]=temp;
         }
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0){
         error("Invers: Singular Matrix-2");
         return(-1);
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
         if (ll != icol) {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
         }
   }
   for (l=n;l>=1;l--) {
      if (indxr[l] != indxc[l]){
         for (k=1;k<=n;k++){
            double temp     = a[k][indxr[l]];
            a[k][indxr[l]] = a[k][indxc[l]];
            a[k][indxc[l]] = temp;
         }
      }
   }
   return(0);
}


/* end of the sources of deal */

/* Main function */


void engen2008r(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, int *idr, int *nidr, int *nsimr, 
		double *resr, int *nsimrar)
{
    /* declaration of variables */
    double **av, **us, **nsco, *varR, tmp, res, **var, mu1, mu2, **Akmo;
    double **inv1, **zer, *a, sigkk, *Wk, *tmp2, m, s, **nscob, **nscoav;
    double *obs, **sammr, vartot, sig2, *Zi, *thetai, tmp3;
    double **mu;
    int nla, nlu, nc, *id, i, j, k, nid, nsim, e, r,l, *index;
    int *indexR, nsimra, b, *ni;
    
    /* memory allocation */
    nc=*ncol;
    nla=*nliga;
    nlu=*nligu;
    nsim=*nsimr;
    nid=*nidr;
    nsimra=*nsimrar;
    e=0;
    r=1;
    tmp=0.0;
    tmp3=0.0;
    b=1;

    /* used and available */
    taballoc(&av, nla, nc);
    taballoc(&us, nlu, nc);
    
    /* the id */
    vecintalloc(&id, nlu);
    
    /* the results */
    taballoc(&mu, nsimra, (nc * 5) );
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	vecalloc(&Wk, (nc-1));
	taballoc(&Akmo, nc-1, nc-1);
	vecalloc(&a, (nc-1));
	vecalloc(&tmp2, (nc-1));
	vecintalloc(&index, nc);
	taballoc(&var, nc, nc);
	taballoc(&inv1, nc-1, nc-1);
	taballoc(&zer, nc-1, nc-1);
    }
    taballoc(&nscob, nlu, nc);

    /* elements required for the calculations
       of the statistics */
    taballoc(&sammr, nsim, nid);
    vecintalloc(&ni, nid);
    vecalloc(&Zi, nid);
    vecalloc(&thetai, nid);

    /* elements required for used normal scores */
    taballoc(&nsco, nlu, nc);
    
    /* elements required for available normal scores */
    taballoc(&nscoav, nla, nc);

    /* elements used in the R API */
    varR = (double *) R_alloc(nla, sizeof(double));
    indexR = (int *) R_alloc(nla, sizeof(int));



        
    /* R -> C objects */
    k=0;
    for (i=1; i<=nla; i++) {
	for (j=1; j<=nc; j++) {
	    av[i][j]=avr[k];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    us[i][j]=usr[k];
	    k++;
	}
    }
    
    for (i=1; i<=nlu; i++) {
	id[i]=idr[i-1];
    }

    /* For each randomization */
    for (b=1; b<=nsimra; b++) {
	
	/* Compute "used" normal score */
	
	/* for each variable */
	for (i=1; i<=nc; i++) {
	    
	    /* sort the variable */
	    for (j = 1; j <= nla; j++) {
		varR[j-1]=av[j][i];
		indexR[j-1]=j;
	    }
	    
	    /* add a very small value to varR (kind of jitter) */
	    /* first find the smallest lag between successive values */
	    tmp = fabs(varR[1]-varR[0]);
	    
	    for (k = 2; k<=nla; k++) {
		for (l = 1; l<k; l++) {
		    if (fabs(varR[k]-varR[l]) < tmp) {
			tmp=fabs(varR[k]-varR[l]);
		    }			    
		}
	    }
	    
	    tmp=tmp/100;
	    
	    for (j = 1; j <= nla; j++) {
		GetRNGstate();
		varR[j-1] = varR[j-1]+ (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
	    }
	    
	    rsort_with_index(varR, indexR, nla);
	    
	    /* normal score for each available point */
	    for (k = 1; k <=nla; k++) {
		nscoav[indexR[k-1]][i] = qnorm((((double) k)/(((double) nla) +1)), 0.0, 1.0, 1, 0);
	    }
	    
	    
	    for (k = 1; k <= nlu; k++) {
		
		/* add a random value to the used points, to avoid ties */
		GetRNGstate();
		res=us[k][i] + (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
		r=1;

		for (j = 1; j <= nla; j++) {
		    
		    /* if the obs is larger than the point, increase r */
		    if (varR[j-1] < res)
			r=j;
		}
		
		tmp3 = (((double) r) + 0.5) / ( ((double) nla) + 1.0);
		/* a qnorm on this, and store the result in nsco */
		nsco[k][i] = qnorm(tmp3, 0.0, 1.0, 1, 0);
		
		
	    }
	    
	}
	
	/*
	  Computes the conditional standardized values, using the method
	  of Ripley (1987, p. 98)
	*/
	/* Just in case there is only one variable */
	if (nc > 1) {
	    for (k = 1; k <= nc; k++) {
		
		l=1;
		for (j=1; j<=nc; j++) {
		    if (j!=k) {
			index[l]=j;
			l++;
		    }
		}
		index[nc]=k;
		
		/* variance covariance matrix of the available points */
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			var[j][l]=0.0;
		    }
		}
		
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			mu1=0.0;
			mu2=0.0;
			
			for (i=1; i<=nla; i++) {
			    mu1 = mu1 + (nscoav[i][index[j]]/((double) nla));
			    mu2 = mu2 + (nscoav[i][index[l]]/((double) nla));
			}
			
			for (i=1; i<=nla; i++) {
			    var[j][l] = var[j][l] + 
				(((nscoav[i][index[j]] - mu1) * 
				  (nscoav[i][index[l]] - mu2))/
				 (((double) nla) - 1.0));
			}
		    }
		}
		
		/* Calculation of Akmo */
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			inv1[j][l] = var[j][l];
		    }
		}
		
		i=invers(inv1, (nc-1), zer, (nc-1));
		
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			Akmo[j][l] = inv1[j][l];
		    }
		}
		
		/* calculation of a */
		for (j = 1; j <= (nc-1); j++)
		    a[j]=var[j][nc];
		
		/* sigkk */
		sigkk=var[nc][nc];
		
		/* for each observation */
		for (i = 1; i <= nlu; i++) {
		    
		    /* Wk */
		    for (j=1; j<= (nc-1); j++)
			Wk[j]=nsco[i][index[j]];
		    
		    /* calculation of the conditionnal average */
		    for (j=1; j<= (nc-1); j++) {
			tmp2[j]=0;
			for (l=1; l<= (nc-1); l++) {
			    tmp2[j] = tmp2[j] + (a[l] *  Akmo[l][j]);
			}
		    }
		    
		    m=0.0;
		    for (l=1; l<= (nc-1); l++) {
			m = m + (tmp2[l] * Wk[l]);
		    }
		    
		    /* Calculation of the conditional variance */
		    s = 0.0;
		    for (l=1; l<= (nc-1); l++) {
			s = s + (tmp2[l] * a[l]);
		    }
		    s = sqrt(sigkk - s);
		    
		    
		    nscob[i][k]=(nsco[i][k] - m)/s;
		}
		
	    }
	} else {
	    nscob[i][1] = nsco[i][1];	    
	}
	

	
	/* ************************
	   Compute the statistics for each individual 
	*/

	/* the number of observations for each animal */
	for (l=1; l<=nid; l++) {
	    e=0;
	    for (i=1; i<=nlu; i++) {
		if (id[i]==l) {
		    e++;
		}
	    }
	    ni[l]=e;
	}
	
	
	/* for each variable */
	for (j=1; j<=nc; j++) {
	    
	    /* for each animal */
	    s=0.0;
	    for (l=1; l<=nid; l++) {
		
		/* get the observations */
		vecalloc(&obs, ni[l]);
		
		k=1;
		for (i=1; i<=nlu; i++) {
		    if (id[i]==l) {
			obs[k]=nscob[i][j];
			k++;
		    }
		}
		
		/* then compute the elements needed for the total variance */
		for (i=1; i<=nsim; i++) {
		    
		    GetRNGstate();
		    tmp= unif_rand();
		    PutRNGstate();
		    for (k = 1; k <= ni[l]; k++) {
			if (tmp < (  ( ((double) k) ) / ((double) ni[l])  )  ) {
			    if (tmp >= (  (((double) k-1)) / ((double) ni[l])  )  ) {
				sammr[i][l]= obs[k];
			    }
			}
		    }
		}
		
		/* and finally, the within class variance */
		m=0.0;
		for (i=1; i<=ni[l]; i++) {
		    m=m+(obs[i] / ((double) ni[l]));
		}
		for (i=1; i<=ni[l]; i++) {
		    s=s+ (((obs[i]-m) * (obs[i]-m)) / 
			  (((double) nlu) - ((double) ni[l])));
		}
		Zi[l]=m;
	    }
	    
	    /* the total variance */
	    vartot=0.0;
	    for (l=1; l<=nsim; l++) {
		/* the mean */
		m=0.0;
		for (i=1; i<=nid; i++) {
		    m=m+ (sammr[l][i]/((double) nid));
		}
		/* the total variance */
		
		for (i=1; i<=nid; i++) {
		vartot=vartot+(((sammr[l][i] - m) * 
				(sammr[l][i] - m))/
			       (((double) nid -1.0) * ((double) nsim)));
		}
	    }
	    sig2=vartot-s;
	    if (sig2<0)
		sig2=0.0;
	    mu[b][((j-1) * 5) + 1]=sammr[1][1];
	    mu[b][((j-1) * 5) + 2]=s;
	    
	    /* value of rho */
	    mu[b][((j-1) * 5) + 3]=sig2/(sig2+s);
	    
	    
	    /* thetai */
	    for (l=1; l<=nid;l++) {
		thetai[l]=sig2 + (s/((double) ni[l]));
	    }
	    
	    /* mu and var(mu) */
	    for (l=1; l<=nid;l++) {
		mu[b][((j-1) * 5) + 4]=mu[b][((j-1) * 5) + 4]+(Zi[l]/thetai[l]);
		mu[b][((j-1) * 5) + 5]=mu[b][((j-1) * 5) + 5]+(1.0/thetai[l]);
	    }
	    mu[b][((j-1) * 5) + 5]=1/mu[b][((j-1) * 5) + 5];
	    mu[b][((j-1) * 5) + 4]=mu[b][((j-1) * 5) + 4] * mu[b][((j-1) * 5) + 5];
	    
	    freevec(obs);
	
	}
	
	/* end of randomization */
    }
    
    k=0;
    for (i=1; i<=nsimra; i++) {
	for (j=1; j<=(nc*5); j++) {
	    resr[k]=mu[i][j];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    usr[k]=nscob[i][j];
	    k++;
	}
    }
 
    /* free memory */

    /* used and available */
    freetab(av);
    freetab(us);
    
    /* the id */
    freeintvec(id);
    
    /* the results */
    freetab(mu);
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	freevec(Wk);
	freetab(Akmo);
	freevec(a);
	freevec(tmp2);
	freeintvec(index);
	freetab(var);
	freetab(inv1);
	freetab(zer);
    }
    freetab(nscob);
    
    /* elements required for the calculations
       of the statistics */
    freetab(sammr);
    freeintvec(ni);
    freevec(Zi);
    freevec(thetai);

    /* elements required for used normal scores */
    freetab(nsco);

    /* elements required for available normal scores */
    freetab(nscoav);


}











/* Extension of the metod for designs I */

void engen2008Ir(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, double *resr, int *nsimrar)
{
    /* declaration of variables */
    double **av, **us, **nsco, *varR, tmp, res, **var, mu1, mu2, **Akmo;
    double **inv1, **zer, *a, sigkk, *Wk, *tmp2, m, s, **nscob, **nscoav;
    double tmp3, **mu;
    int nla, nlu, nc, i, j, k, r,l, *index;
    int *indexR, nsimra, b;
    
    /* memory allocation */
    nc=*ncol;
    nla=*nliga;
    nlu=*nligu;
    nsimra=*nsimrar;
    r=1;
    tmp=0.0;
    tmp3=0.0;
    b=1;

    /* used and available */
    taballoc(&av, nla, nc);
    taballoc(&us, nlu, nc);
        
    /* the results */
    taballoc(&mu, nsimra, (nc * 2) );
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	vecalloc(&Wk, (nc-1));
	taballoc(&Akmo, nc-1, nc-1);
	vecalloc(&a, (nc-1));
	vecalloc(&tmp2, (nc-1));
	vecintalloc(&index, nc);
	taballoc(&var, nc, nc);
	taballoc(&inv1, nc-1, nc-1);
	taballoc(&zer, nc-1, nc-1);
    }
    taballoc(&nscob, nlu, nc);


    /* elements required for used normal scores */
    taballoc(&nsco, nlu, nc);
    
    /* elements required for available normal scores */
    taballoc(&nscoav, nla, nc);

    /* elements used in the R API */
    varR = (double *) R_alloc(nla, sizeof(double));
    indexR = (int *) R_alloc(nla, sizeof(int));



        
    /* R -> C objects */
    k=0;
    for (i=1; i<=nla; i++) {
	for (j=1; j<=nc; j++) {
	    av[i][j]=avr[k];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    us[i][j]=usr[k];
	    k++;
	}
    }
    

    /* For each randomization */
    for (b=1; b<=nsimra; b++) {
	
	/* Compute "used" normal score */
	
	/* for each variable */
	for (i=1; i<=nc; i++) {
	    
	    /* sort the variable */
	    for (j = 1; j <= nla; j++) {
		varR[j-1]=av[j][i];
		indexR[j-1]=j;
	    }
	    
	    /* add a very small value to varR (kind of jitter) */
	    /* first find the smallest lag between successive values */
	    tmp = fabs(varR[1]-varR[0]);
	    
	    for (k = 2; k<=nla; k++) {
		for (l = 1; l<k; l++) {
		    if (fabs(varR[k]-varR[l]) < tmp) {
			tmp=fabs(varR[k]-varR[l]);
		    }			    
		}
	    }
	    
	    tmp=tmp/100;
	    
	    for (j = 1; j <= nla; j++) {
		GetRNGstate();
		varR[j-1] = varR[j-1]+ (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
	    }
	    
	    rsort_with_index(varR, indexR, nla);
	    
	    /* normal score for each available point */
	    for (k = 1; k <=nla; k++) {
		nscoav[indexR[k-1]][i] = qnorm((((double) k)/(((double) nla) +1)), 0.0, 1.0, 1, 0);
	    }
	    
	    
	    for (k = 1; k <= nlu; k++) {
		
		/* add a random value to the used points, to avoid ties */
		GetRNGstate();
		res=us[k][i] + (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
		r=1;

		for (j = 1; j <= nla; j++) {
		    
		    /* if the obs is larger than the point, increase r */
		    if (varR[j-1] < res)
			r=j;
		}
		
		tmp3 = (((double) r) + 0.5) / ( ((double) nla) + 1.0);
		/* a qnorm on this, and store the result in nsco */
		nsco[k][i] = qnorm(tmp3, 0.0, 1.0, 1, 0);
		
		
	    }
	    
	}
	
	/*
	  Computes the conditional standardized values, using the method
	  of Ripley (1987, p. 98)
	*/
	/* Just in case there is only one variable */
	if (nc > 1) {
	    
	    /* for each variable */
	    for (k = 1; k <= nc; k++) {
		
		l=1;
		for (j=1; j<=nc; j++) {
		    if (j!=k) {
			index[l]=j;
			l++;
		    }
		}
		index[nc]=k;
		
		/* variance covariance matrix of the available points */
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			var[j][l]=0.0;
		    }
		}
		
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			mu1=0.0;
			mu2=0.0;
			
			for (i=1; i<=nla; i++) {
			    mu1 = mu1 + (nscoav[i][index[j]]/((double) nla));
			    mu2 = mu2 + (nscoav[i][index[l]]/((double) nla));
			}
			
			for (i=1; i<=nla; i++) {
			    var[j][l] = var[j][l] + 
				(((nscoav[i][index[j]] - mu1) * 
				  (nscoav[i][index[l]] - mu2))/
				 (((double) nla) - 1.0));
			}
		    }
		}
		
		/* Calculation of Akmo */
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			inv1[j][l] = var[j][l];
		    }
		}
		
		i=invers(inv1, (nc-1), zer, (nc-1));
		
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			Akmo[j][l] = inv1[j][l];
		    }
		}
		
		/* calculation of a */
		for (j = 1; j <= (nc-1); j++)
		    a[j]=var[j][nc];
		
		/* sigkk */
		sigkk=var[nc][nc];
		
		/* for each observation */
		for (i = 1; i <= nlu; i++) {
		    
		    /* Wk */
		    for (j=1; j<= (nc-1); j++)
			Wk[j]=nsco[i][index[j]];
		    
		    /* calculation of the conditionnal average */
		    for (j=1; j<= (nc-1); j++) {
			tmp2[j]=0;
			for (l=1; l<= (nc-1); l++) {
			    tmp2[j] = tmp2[j] + (a[l] *  Akmo[l][j]);
			}
		    }
		    
		    m=0.0;
		    for (l=1; l<= (nc-1); l++) {
			m = m + (tmp2[l] * Wk[l]);
		    }
		    
		    /* Calculation of the conditional variance */
		    s = 0.0;
		    for (l=1; l<= (nc-1); l++) {
			s = s + (tmp2[l] * a[l]);
		    }
		    s = sqrt(sigkk - s);
		    
		    
		    nscob[i][k]=(nsco[i][k] - m)/s;
		}
		
	    }
	} else {
	    /* if one variable, no correction */
	    nscob[i][1] = nsco[i][1];	    
	}
	

	
	/* ************************
	   Compute the statistics
	*/

	/* for each variable */
	for (j=1; j<=nc; j++) {
	    
	    m=0.0;
	    s=0.0;
	    	    
	    /* The preference */
	    for (i=1; i<=nlu; i++) {
		m = m + (nscob[i][j]);
	    }
	    
	    /* the variance of the scores */
	    for (i=1; i<=nlu; i++) {
		s=s+ (((nscob[i][j]-m) * (nscob[i][j]-m)) / 
		      (((double) nlu) - 1));
	    }
	    
	    /* variance of the mean */
	    s=s/((double) nlu);
	    
	    mu[b][((j-1) * 2) + 1]=m;
	    mu[b][((j-1) * 2) + 2]=s;
	}
	
	/* end of randomization */
    }
    
    k=0;
    for (i=1; i<=nsimra; i++) {
	for (j=1; j<=(nc*2); j++) {
	    resr[k]=mu[i][j];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    usr[k]=nscob[i][j];
	    k++;
	}
    }
 
    /* free memory */

    /* used and available */
    freetab(av);
    freetab(us);
    
    
    /* the results */
    freetab(mu);
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	freevec(Wk);
	freetab(Akmo);
	freevec(a);
	freevec(tmp2);
	freeintvec(index);
	freetab(var);
	freetab(inv1);
	freetab(zer);
    }
    freetab(nscob);
    

    /* elements required for used normal scores */
    freetab(nsco);

    /* elements required for available normal scores */
    freetab(nscoav);


}

