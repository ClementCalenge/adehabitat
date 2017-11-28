#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void aclambda(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bufligr(void *, void *, void *, void *, void *, void *, void *, void *);
extern void calcvolume(void *, void *, void *, void *);
extern void clusterhrr(void *, void *, void *, void *, void *, void *);
extern void CVL(void *, void *, void *, void *, void *, void *, void *);
extern void CVmise(void *, void *, void *, void *, void *, void *);
extern void discretrajr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void distxyr(void *, void *, void *, void *, void *);
extern void engen2008Ir(void *, void *, void *, void *, void *, void *, void *);
extern void engen2008r(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void erodil(void *, void *, void *, void *, void *);
extern void fctdomain(void *, void *, void *, void *, void *, void *, void *);
extern void findmaxgrid(void *, void *, void *);
extern void fipatir(void *, void *, void *, void *, void *, void *, void *);
extern void getcontourc(void *, void *, void *, void *, void *, void *);
extern void kernelbbc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernelhr(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernelkcr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernepan(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void lcontour(void *, void *, void *, void *);
extern void longfacclustr(void *, void *, void *);
extern void nls2k(void *, void *, void *, void *, void *);
extern void optcutr(void *, void *, void *, void *, void *, void *);
extern void partrajr(void *, void *, void *, void *, void *, void *, void *);
extern void permutksel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void permutR2n(void *, void *, void *, void *, void *, void *);
extern void prepquart(void *, void *, void *, void *, void *, void *, void *);
extern void randenfar(void *, void *, void *, void *, void *, void *);
extern void randmargtolpts(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rankma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rastpolaire(void *, void *, void *, void *, void *, void *, void *, void *);
extern void regrouascnumr(void *, void *, void *, void *, void *, void *);
extern void regroufacascr(void *, void *, void *, void *, void *, void *, void *, void *);
extern void runsltr(void *, void *, void *, void *);
extern void sahr2ksel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void seqeticorr(void *, void *, void *);
extern void testindepangl(void *, void *, void *, void *, void *, void *, void *);
extern void testindepdist(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"aclambda",       (DL_FUNC) &aclambda,       11},
    {"bufligr",        (DL_FUNC) &bufligr,         8},
    {"calcvolume",     (DL_FUNC) &calcvolume,      4},
    {"clusterhrr",     (DL_FUNC) &clusterhrr,      6},
    {"CVL",            (DL_FUNC) &CVL,             7},
    {"CVmise",         (DL_FUNC) &CVmise,          6},
    {"discretrajr",    (DL_FUNC) &discretrajr,    13},
    {"distxyr",        (DL_FUNC) &distxyr,         5},
    {"engen2008Ir",    (DL_FUNC) &engen2008Ir,     7},
    {"engen2008r",     (DL_FUNC) &engen2008r,     10},
    {"erodil",         (DL_FUNC) &erodil,          5},
    {"fctdomain",      (DL_FUNC) &fctdomain,       7},
    {"findmaxgrid",    (DL_FUNC) &findmaxgrid,     3},
    {"fipatir",        (DL_FUNC) &fipatir,         7},
    {"getcontourc",    (DL_FUNC) &getcontourc,     6},
    {"kernelbbc",      (DL_FUNC) &kernelbbc,      13},
    {"kernelhr",       (DL_FUNC) &kernelhr,        9},
    {"kernelkcr",      (DL_FUNC) &kernelkcr,      10},
    {"kernepan",       (DL_FUNC) &kernepan,        9},
    {"lcontour",       (DL_FUNC) &lcontour,        4},
    {"longfacclustr",  (DL_FUNC) &longfacclustr,   3},
    {"nls2k",          (DL_FUNC) &nls2k,           5},
    {"optcutr",        (DL_FUNC) &optcutr,         6},
    {"partrajr",       (DL_FUNC) &partrajr,        7},
    {"permutksel",     (DL_FUNC) &permutksel,     19},
    {"permutR2n",      (DL_FUNC) &permutR2n,       6},
    {"prepquart",      (DL_FUNC) &prepquart,       7},
    {"randenfar",      (DL_FUNC) &randenfar,       6},
    {"randmargtolpts", (DL_FUNC) &randmargtolpts, 16},
    {"rankma",         (DL_FUNC) &rankma,         10},
    {"rastpolaire",    (DL_FUNC) &rastpolaire,     8},
    {"regrouascnumr",  (DL_FUNC) &regrouascnumr,   6},
    {"regroufacascr",  (DL_FUNC) &regroufacascr,   8},
    {"runsltr",        (DL_FUNC) &runsltr,         4},
    {"sahr2ksel",      (DL_FUNC) &sahr2ksel,      10},
    {"seqeticorr",     (DL_FUNC) &seqeticorr,      3},
    {"testindepangl",  (DL_FUNC) &testindepangl,   7},
    {"testindepdist",  (DL_FUNC) &testindepdist,   7},
    {NULL, NULL, 0}
};

void R_init_adehabitat(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
