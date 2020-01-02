/* -*- coding: utf-8 -*

   Graph Generator                                   Â© Cyril Gavoille


     Use: gengraph [-options] graph [parameters]


     A free command-line program to generate graphs in many formats:
     plain text, .dot, .pdf, .fig, .svg â€¦ Formats like .pdf or .fig
     are produced thanks to GraphViz and allow visualization. Easy to
     install, there is a single .c file to compile. There is a on-line
     manual available in the command (for the moment the manual is in
     French only and included at the end of the source).

                                â”€â”€â”€â”€â”€

     Petit programme libre en ligne de commande permettant de gÃ©nÃ©rer
     des graphes dans plein de formats diffÃ©rents: texte, .dot, .pdf,
     .fig, .svg â€¦ Certains formats (comme .pdf ou .jpg), produits
     grÃ¢ce Ã  GraphViz, permettent de visualiser les graphes.  TrÃ¨s
     simple Ã  installer, il n'y a qu'un seul fichier .c Ã  compiler.
     Pour le reste il faut voir l'aide en ligne qui est aussi incluse
     Ã  la fin du source.


   Comment l'installer / le compiler ?

     (MacOs)  gcc gengraph.c -o gengraph
     (Linux)  gcc gengraph.c -lm -lbsd -o gengraph

*/


#define _GNU_SOURCE   // pour asprintf()
#include <stdio.h>    // pour printf(), sprintf() ...
#include <stdlib.h>   // pour system(), strtod(), RAND_MAX ...
#include <string.h>   // pour strcomp() ...
#include <unistd.h>   // pour getpid() ...
#include <math.h>     // pour sqrt(), cos(), acos() ...
#include <float.h>    // pour DBL_MAX, DBL_DIG ...
#include <limits.h>   // pour INT_MAX, LONG_MAX ...
#include <time.h>     // pour time(), clock(), ...
#include <sys/time.h> // pour gettimeofday()
#include <assert.h>   // pour assert() dynamique

#if defined(__APPLE__) && defined(__MACH__) // Apple OSX and iOS (Darwin)
#else
#include <bsd/stdlib.h> // pour arc4random() Linux
#endif

/*
  gcc gengraph.c -Wall -DDEBUG -o gengraph  --> dÃ©bugage avec DEBUG(I)
  gcc gengraph.c -Wall -S -g -->  gÃ©nÃ¨re le code assembleur gengraph.s 
  gcc gengraph.c -Wall -E -P --> gÃ©nÃ¨re l'expansion des macros
*/
#ifdef DEBUG
#undef DEBUG
#define DEBUG(I) do{ I }while(0)
#else
#define DEBUG(I)
#undef assert
#define assert(C)
#endif


/*
  Permet de tester une condition C Ã  la compilation sans rien ajouter
  au code. Si C est fausse, la compilation indique une erreur
  contenant le numÃ©ro de ligne, du type:

  gengraph.c:18162:3: error: array size is negative
*/
#define ASSERT(C) ((void)sizeof(int[-!(C)]))


typedef char* string; /* chaÃ®nes de caractÃ¨res terminÃ©es par 0 */

/* graphe ou famille de graphes */

typedef struct _graph {
  int id;    // numÃ©ro du graphe, utilisÃ© pour les familles de graphes 
  int n;     // nb de sommets, <0 si non dÃ©fini
  int m;     // nb d'arÃªtes, <0 si non dÃ©fini
  int *d;    // d[u]=degrÃ© du sommet u
  int **L;   // L[u][i]=i-Ã¨me voisin du sommet u, i=0...d[u]-1
  int sort;  // vrai ssi les listes d'adjacences sont triÃ©es
  int sym;   // vrai ssi les listes d'adjacence du graphe sont symÃ©triques
  double **W;// W[u][i]=poids du i-Ã¨me voisin du sommet u, i=0...d[u]-1
  double *xpos,*ypos; // tableau de positions des sommets (graphe gÃ©omÃ©trique)
  int f;     // nombre de graphes de la famille, =0 si graphe normal
  struct _graph **G; // G[i]=i-Ã¨me graphe de la famille, i=0..f-1

  // les champs suivants sont utilisÃ©s pour communiquer des valeurs
  // ou des rÃ©sultats Ã  travers les appels de fonctions

  int int1;  // paramÃ¨tre: entier
  int* pint1;// paramÃ¨tre: tableau d'entiers
} graph;



/* fonction de test sur un graphe */

typedef int test(graph*);



/* chemin simple d'un graphe G */

typedef struct{
  int n;  /* nombre de sommets du chemin */
  int *P; /* liste ordonnÃ©e des sommets du chemin */
  int *V; /* V[u]=i si le sommet u de G est le i-Ã¨me (i dans [0,G->n[)
	     sommet du chemin (V[P[i]]=i), V[u]=-1 si u n'est dans le
	     chemin  */
  int **aux; /* tableau auxiliaire utilisÃ© (et gÃ©rÃ©) par NextPath() */
} path;



/* token pour le type "list" */

enum{
  T_NODE,  // item est un sommet u
  T_EDGE,  // item est un sommet v connectÃ© au prÃ©cÃ©dent de la liste par une arÃªte: u-v
  T_ARC,   // item est un sommet v connectÃ© au prÃ©cÃ©dent de la liste par un arc: u->v
  T_ID,    // item est l'identifiant d'un nouveau graphe d'une famille
  T_NB,    // item est le nombre de graphes de la famille, item toujours en tÃªte de liste
  T_OPENE, // item est le premier sommet v d'un groupe arÃªte: u-(v ...
  T_OPENA, // item est le premier sommet v d'un groupe arc: u->(v ...

  // token non encore gÃ©rÃ©

  T_UNIV,  // u-* = item est un sommet universel u-(0 ... n)
  T_UNOUT, // u->* item est un sommet universel u->(0 ... n-1)
  T_UNIN,  // *->u item est un sommet universel u<-(0 ... n-1)
};


/* code pour les fonctions de hashage */

enum{
  H_PRIME,
  H_SHUFFLE,
  H_MIX,
  H_MOD,
};


/* code pour les normes et graphes gÃ©omÃ©triques */

enum{
  NORM_FAIL,  // norme indÃ©terminÃ©e
  NORM_L1,    // norme L1
  NORM_L2,    // norme L2
  NORM_LMAX,  // norme Lmax
  NORM_LMIN,  // norme Lmin
  NORM_POLY,  // norme polygonale
  NORM_HYPER, // norme hyperbolique
};


typedef struct _list { int item,type; struct _list *next; } list; // liste chaÃ®nÃ©e
typedef struct{ double x,y; } point; // un point du plan
typedef struct{ int u,v; double w; } edge; // une arÃªte avec un poids
typedef struct{ int i,j; } couple; // un couple d'indices
typedef struct{ int x,y; long z; } triplet; // un triplet d'entiers (avec un long)
typedef struct{ unsigned r:8,g:8,b:8; } color; // couleur avec champs r,g,b de 8 bits


/* constantes pour la ligne de commande */

#define PARMAX      64
/* C'est le nombre maximum de paramÃ¨tres pour un graphe ayant un
   nombre de paramÃ¨tres fixÃ©s (non variable). Le plus grand nombre de
   paramÃ¨tres pour un tel graphe est de 62 (graphe kakutami_7x+1). */

#define DIMAX       64 /* dimension maximum pour un graphe (hypercube, pancake, ...) */
#define CMDMAX    1024 /* nb maximum de caractÃ¨res sur la ligne de commande */
#define PARAMSIZE 1024 /* taille mÃ©moire des buffers FPARAM et CPARAM (en octets) */
#define SURFACEMAX 256 /* nombre maximum de coutures d'une surface */

#define NAMEMAX    256 
/* C'est le nombre maximum de caractÃ¨res pour un nom de sommet y
   compris le "\0" terminal. Il faut que cela soit suffisament grand
   pour contenir un double, entier ou caractÃ¨res simple, 20 caractÃ¨res
   en pratique suffise. */

/* constantes pour le format dot et -visu */

double VSIZESTD=0.05; /* taille standard des sommets */
double VSIZEXY=0.12;  /* taille des sommets pour les graphes gÃ©omÃ©triques */
double LEN=1.0; /* longueur des arÃªtes avec dot neato */
#define URL_vis_js1  "vis.min.js"
#define URL_vis_js2  "https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js"
#define URL_vis_css1 "vis.min.css"
#define URL_vis_css2 "https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.css"
#define GRAPH_PDF  g.pdf  /* nom du graphe de sortie par dÃ©faut */
#define GRAPH_HTML g.html /* nom du graphe de sortie html par dÃ©faut */
#define VSIZEK 5.5 /* rapport entre taille max et min des sommets (-vsize) */
#define C32 32
/* Coefficient par lequel sont multipliÃ©es les coordonnÃ©es des points
   dans le cas des graphes gÃ©omÃ©triques et le format dot. Plus la
   valeur est Ã©levÃ©e, plus les sommets paraissent petits et les arÃªtes
   fines et longues */


/* constantes diverses */

#define CHRONOMAX 5 /* nombre maximum de chronomÃ¨tres */
#define RS_NI_FP "name-independent, fixed-port model"
#define RS_L_FP  "labeled, fixed-port model"
#define M_2PI 6.28318530717958647692528676655900576839 /* constante 2ðœ‹, souvent utilisÃ©e */
#define EULER_MASCHERONI 0.5772156649015328606 /* sert pour func1() */


/* codes pour les formats de sorties possibles */

enum{
  F_standard,
  F_list,
  F_matrix,
  F_smatrix,
  F_dot,
  F_userdot,
  F_html,
  F_xy,
  F_no,
};


/*
  Codes pour les algorithmes via l'option -check. La variable CHECK
  vaut l'une de ces constantes. Par dÃ©faut, CHECK=CHECK_OFF(=0). Si
  CHECK>0, alors le graphe sera stockÃ© en mÃ©moire. Si CHECK>1, alors
  en plus l'algorithme spÃ©cifique sera appliquÃ© au graphe de sortie.
*/

enum{
  CHECK_OFF=0, // valeur par dÃ©faut
  CHECK_ON,    // graphe chargÃ© en mÃ©moire
  CHECK_BFS,
  CHECK_DFS,
  CHECK_NCC,
  CHECK_BELLMAN,
  CHECK_STRETCH,
  CHECK_DEGENERATE,
  CHECK_GCOLOR,
  CHECK_DEG,
  CHECK_ISO,
  CHECK_SUB,
  CHECK_ISUB,
  CHECK_MINOR,
  CHECK_TWDEG,
  CHECK_TW,
  CHECK_PS1,
  CHECK_PS1b,
  CHECK_PS1c,
  CHECK_PS1x,
  CHECK_GIRTH,
  CHECK_RADIUS,
  CHECK_DIAMETER,
  CHECK_VOLM,
  CHECK_PATHS,
  CHECK_MAINCC,
  CHECK_SUBDIV,
  CHECK_KCOLOR,
  CHECK_KCOLORSAT,
  CHECK_KINDEPSAT,
  CHECK_INFO,
  CHECK_SIMPLIFY,
  CHECK_RS_CLUSTER,
  CHECK_RS_DCR,
  CHECK_RS_AGMNT,
  CHECK_RS_TZRPLG,
  CHECK_RS_BC,
  CHECK_RS_HDLBR,
};


/*
  Codes pour la variable XYtype indiquant le type de
  distribution/gÃ©nÃ©ration des points (xpos,ypos) des graphes
  gÃ©omÃ©triques.
*/

enum{
  XY_FILE,    // points lus Ã  partir d'un fichier
  XY_UNIF,    // points selon une loi uniforme
  XY_PLAW,    // points selon une loi en puissance
  XY_PERM,    // points selon une permutation ðœ‹, (i,ðœ‹(i))
  XY_MESH,    // points sur une grille Xmesh x Ymesh
  XY_CIRCLE,  // points alÃ©atoires sur un cercle de rayon 1
  XY_CYCLE,   // points rÃ©guliers sur un cercle de rayon 1
  XY_RPOLY,   // points alÃ©atoire dans un polygone convexe rÃ©gulier
  XY_DISK,    // points dans le disque unitÃ© triÃ©s selon l'angle
  XY_HYPER,   // points selon une loi en puissance dans le disque unitÃ©
  XY_CONVEX,  // points dans le disque en position convexe
  XY_CONVEX2, // variante de XY_CONVEX
  XY_USER,    // points dÃ©finis par la fonction d'adjacence
};


/********************/
/* fonctions inline */
/********************/

/*
  Renvoie l'entier ceil(x/y) lorsque x,y sont entiers et y<>0. Marche
  mÃªme pour les entiers nÃ©gatifs.
*/

static inline int iceil(int const x,int const y){ return (x/y) + (((x<0)^(y>0))&&(x%y)); }


/***********************************/
/* macros & pseudo-fonction utiles */
/***********************************/

/*
  Macro permettant d'utiliser __auto_type si la version de GCC le
  permet. Sinon, on utilise typeof(x). Il est en principe prÃ©fÃ©rable
  d'utiliser "__auto_type var = (x);" plutÃ´t que "typeof(x) var =
  (x);", surtout dans une macro justement, car cela Ã©vite la double
  Ã©valuation de x. Pour utiliser __auto_type il faut la version de GCC
  >= 4.9, ce qu'on peut aussi tester avec #if GCC_VERSION >= 4.9 ...
  #endif (constante non disponible sous MacOS).
*/

#ifdef __auto_type
#define auto_typeof(x) __auto_type
#else
#define auto_typeof(x) typeof(x)
#endif

/*
  Macros min/max qui provoquent un warning Ã  la compilation si les
  valeurs ne sont pas du mÃªme type, ce qui Ã©vite les erreurs comme
  "(2.3<=2.1)" Ã©valuÃ© Ã  vrai si jamais les arguments sont convertis en
  int comme le ferait une fonction inline. Les paramÃ¨tres sont Ã©valuÃ©s
  une seule fois (si __auto_type existe) comme une fonction inline. Le
  rÃ©sultat est du mÃªme type que les arguments.
*/

#define min(x, y) ({		\
      auto_typeof(x) _x = (x);	\
      auto_typeof(x) _y = (y);	\
      (void) (&_x == &_y);	\
      _x < _y ? _x : _y; })

#define max(x, y) ({		\
      auto_typeof(x) _x = (x);	\
      auto_typeof(y) _y = (y);	\
      (void) (&_x == &_y);	\
      _x > _y ? _x : _y; })


#define EQUAL(s)   (strcmp(ARGV[i],s)==0)
#define PREFIX(s)  (strncmp(ARGV[i],s,strlen(s))==0)
#define SUFFIX(s)  (strncmp(ARGV[i]+max((int)(strlen(ARGV[i])-strlen(s)),0),s,strlen(s))==0)
#define PLURIEL(n) (((n)>1)?"s":"")
#define VIDE(s)    *s='\0'    /* vide le tableau de caractÃ¨res s */ 
#define ESTVIDE(s) (*s=='\0') /* vrai ssi s est vide */
#define NONVIDE(s) (*s!='\0') /* vrai ssi s est non vide */
#define MEM(mem,pos,type) (*(type*)(mem+pos)) /* Ã©crit/lit la mÃ©moire mem[pos] */
#define RAND01     ((double)random()/RAND_MAX) /* rÃ©el alÃ©atoire dans [0,1] */
#define RAND62     (random()*(RAND_MAX+1L)+random()) /* entier alÃ©atoire sur 62 bits */
#define RANDbit    (random()&1) /* bit alÃ©atoire */
#define BARRE do{ruling("â€•",80);printf("\n");}while(0) /* barre de sÃ©paration */

/*
  Ã‰change de deux variables d'au plus 16 octets. Une erreur Ã  la
  compilation survient si la variable globale _SWAP n'est pas assez
  grande ou si les valeurs ne sont pas de mÃªme type. 

  NB: Il ne faut surtout pas mettre "static char const" pour la
  variable globale _SWAP, sinon cela fait "Bus error: 10" Ã 
  l'exÃ©cution car l'erreur "read-only variable" n'est pas dÃ©tectÃ©e Ã 
  la compilation Ã  cause du cast.
*/
static char _SWAP[16]; // variable globale pour SWAP
#define SWAP(x,y)						\
  do{								\
    ASSERT(sizeof(x)<=sizeof(_SWAP));				\
    (*(typeof(x)*)_SWAP)=(x),(x)=(y),(y)=(*(typeof(x)*)_SWAP);	\
  }while(0)

/* Appelle qsort(T,...) avec la bonne taille d'Ã©lÃ©ment (taille du type de T[]) */
#define QSORT(T,n,f) qsort(T,n,sizeof(*(T)),f)

/* permet l'expansion d'une macro. Ex: scanf("%"xstr(DMAX)"s",buffer); */
#define xstr(s) str(s)
#define str(s) #s

/*
  Alloue Ã  la variable T un tableau de n>=0 valeurs du type de
  T[]. Sort une erreur en cas de dÃ©passement mÃ©moire. Si n=0, alors T
  est modifiÃ© et allouÃ© Ã  une zÃ´ne de sorte que free(T) ne fera pas
  d'erreur.
*/
#define ALLOC(T,n) do{ if((T=malloc((n)*sizeof(*T)))==NULL) Erreur(3); }while(0)

/*
  Comme ALLOC(T,n) mais initialise en plus le tableau T avec le terme
  z qui peut comporter l'indice _i, en posant T[_i]=(z).
*/
#define ALLOCZ(T,n,z)				\
  do{						\
    size_t const _n=(n);			\
    ALLOC(T,_n);				\
    for(size_t _i=0;_i<_n;_i++) T[_i]=(z);	\
  }while(0)

/* Un realloc() qui Ã©vite de faire un free() non souhaitÃ© Ã  un
   pointeur si n=0, et donc qui Ã©vite les erreurs "double free". */
#define REALLOC(P,n) P=realloc(P,max(n,1)*sizeof(*P))

/*
  Alloue Ã  la variable T la place pour une matrice de n x s valeurs,
  c'est-Ã -dire un tableau de n tableaux de s valeurs du type de
  T[][]. On utilise un seul malloc(), ce qui est souvent beaucoup plus
  rapide. Pour le libÃ©rer, un seul et simple free(T) suffit. Les
  Ã©lÃ©ments sont stockÃ©s Ã  la suite, si bien qu'on peut accÃ©der Ã  tous
  les Ã©lÃ©ments soit par T[0..n[[0..s[ soit par (*T)[0..n*s[.
*/
#define ALLOC2(T,n,s)						\
  do{								\
    size_t _n=(n),_s=(s),_i=0; auto_typeof(*T) _p;		\
    T=malloc(_n*sizeof(*T)+(_n*_s)*sizeof(**T));		\
    if(T==NULL) Erreur(3);					\
    for(_p=(typeof(*T))(T+_n);_i<_n;_i++,_p+=_s) T[_i]=_p;	\
  }while(0)

/*
  Comme ALLOC2(T,n,s) mais initialise en plus le tableau avec le terme
  z qui peut ne comporter qu'un seul indice _i, en posant (*T)[_i]=z.
*/
#define ALLOC2Z(T,n,s,z)			\
  do{						\
    ALLOC2(T,n,s);				\
    size_t const _n=(n)*(s);			\
    for(size_t _i=0;_i<_n;_i++) (*T)[_i]=(z);	\
  }while(0)

/* Les variantes NALLOC permettent de dÃ©clarer le pointeur T sur des
   valeurs de type t avant de faire un ALLOC. Cela Ã©vite donc une
   dÃ©claration de pointeur avant un ALLOC.

   Ex: int *T;      ->   NALLOC(int,T,n);
       ALLOC(T,n);

   Attention de ne pas mettre ces macros dans un bloc { ... }, car
   sinon la portÃ©e de T se limitera Ã  ce block. Les variantes NALLOC
   ne peuvent pas remplacer des ALLOC qui sont dans les configurations
   suivantes:

   - la dÃ©claration et le malloc() ne sont pas dans le mÃªme bloc, comme:
     int *T;
     if(...)
       ALLOC(T,...);

   - la dÃ©claration est aprÃ¨s l'Ã©tiquette d'un goto ou d'un case, comme:
     label:
       int *T;
       ALLOC(T,...);
*/
#define NALLOC(t,T,n)       t*  T; ALLOC(T,n);
#define NALLOCZ(t,T,n,z)    t*  T; ALLOCZ(T,n,z);
#define NALLOC2(t,T,n,s)    t** T; ALLOC2(T,n,s);
#define NALLOC2Z(t,T,n,s,z) t** T; ALLOC2Z(T,n,s,z);

/*
  LibÃ¨re un tableau T de n pointeurs. On ne peut pas faire de
  fonction, car le type de T n'Ã©tant pas prÃ©-dÃ©terminÃ© on ne peut
  Ã©crire le prototype de la fonction. Attention! Il ne faut surtout
  pas faire FREE2(T,n) aprÃ¨s avoir allouÃ© T avec ALLOC2(T,n,s).
*/
#define FREE2(T,n)				\
  do{						\
    if(T){					\
      for(size_t _i=0;_i<(n);_i++) free(T[_i]);	\
      free(T);					\
    }						\
  }while(0)

/*
  Ajoute une arÃªte uv au graphe G Ã  la fin des listes de u et de v, en
  supposant qu'il y a assez de place. Attention ! seuls les champs
  G->d et G->L sont mis Ã  jour. Les champs comme G->m, G->sort peuvent
  ne pas Ãªtre cohÃ©rents.
*/
#define ADD_EDGE(G,u,v)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
    G->L[v][G->d[v]++]=u;	\
  }while(0)

/* Comme ADD_EDGE, mais sans ajouter u Ã  la liste de v. */
#define ADD_ARC(G,u,v)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
  }while(0)


/***********************/
/* macros de dÃ©buggage */
/***********************/

/* Sert pour PRINT() */
#define _PRINT(x,t,f)					\
  else if(__builtin_types_compatible_p(typeof(x),t))	\
    _f = #x" = "f" ("#t")\n"

/* Affichage une variable ou valeur en fonction de son type. */
#define PRINT(x)							\
  do{									\
    char* _f; if(0);							\
    _PRINT(x,char,"%c");						\
    _PRINT(x,char*,"\"%s\"");						\
    _PRINT(x,char[],"\"%s\"");						\
    _PRINT(x,unsigned char,"%u");					\
    _PRINT(x,unsigned char*,"%p");					\
    _PRINT(x,unsigned char[],"%p");					\
    _PRINT(x,short,"%d");						\
    _PRINT(x,short*,"%p");						\
    _PRINT(x,short[],"%p");						\
    _PRINT(x,unsigned short,"%u");					\
    _PRINT(x,unsigned short*,"%p");					\
    _PRINT(x,unsigned short[],"%p");					\
    _PRINT(x,int,"%d");							\
    _PRINT(x,int*,"%p");						\
    _PRINT(x,int[],"%p");						\
    _PRINT(x,unsigned int,"%u");					\
    _PRINT(x,unsigned int*,"%p");					\
    _PRINT(x,unsigned int[],"%p");					\
    _PRINT(x,long,"%li");						\
    _PRINT(x,long*,"%p");						\
    _PRINT(x,long[],"%p");						\
    _PRINT(x,unsigned long,"%lu");					\
    _PRINT(x,unsigned long*,"%p");					\
    _PRINT(x,unsigned long[],"%p");					\
    _PRINT(x,long long,"%li");						\
    _PRINT(x,long long*,"%p");						\
    _PRINT(x,long long[],"%p");						\
    _PRINT(x,unsigned long long,"%lu");					\
    _PRINT(x,unsigned long long*,"%p");					\
    _PRINT(x,unsigned long long[],"%p");				\
    _PRINT(x,float,"%g");						\
    _PRINT(x,float*,"%p");						\
    _PRINT(x,float[],"%p");						\
    _PRINT(x,double,"%g");						\
    _PRINT(x,double*,"%p");						\
    _PRINT(x,double[],"%p");						\
    _PRINT(x,long double,"%Lg");					\
    _PRINT(x,long double*,"%p");					\
    _PRINT(x,long double[],"%p");					\
    _PRINT(x,void*,"%p");						\
    else{ printf(#x " = ? (unknown type)\n"); break; }			\
    printf(_f,(x));							\
  }while(0)

/* affiche une ligne vide */
#define PRINTN printf("\n")

/* affiche un tableau avec un certain format */
#define PRINT_TAB(T,n,f)				\
  do{							\
    printf(#T " =");					\
    if((T)==NULL) printf(" NULL"); else			\
      for(int _i=0;_i<(n);_i++) printf(f,(T)[_i]);	\
    printf("\n");					\
  }while(0)

/* affiche un tableau d'entiers */
#define PRINTT(T,n) PRINT_TAB(T,n," %i")

/* affiche un tableau de doubles */
#define PRINTD(T,n) PRINT_TAB(T,n," %g")

/* affiche une liste chaÃ®nÃ©e */
#define PRINTLIST(L)		        	\
  do{						\
    list *_L=L;                                 \
    printf(#L " = {");				\
    while(_L!=NULL){				\
      printf(" (%i,%i)",_L->item,_L->type);     \
      _L=_L->next;                              \
      } printf(" }\n");				\
  }while(0)


/* Affiche un aperÃ§u d'une liste de n valeurs entiÃ¨res V=V(_i), les
   indices variant de _i=0..n. On affiche les K1 premiÃ¨res valeurs et
   les K2 derniÃ¨res valeurs. La complexitÃ© est en O(K1+K2).

   Ex: APERCU(_i*_i,10,3,2) affichera "[ 0 1 2 ... 64 81 ]".
*/
#define APERCU(V,n,K1,K2)						\
  do{									\
    printf("[ ");							\
    for(int _i=0;_i<n;_i++)						\
      if((_i<K1)||(_i>n-K2-1)||(K1==n-K2-1))				\
	printf("%i ",V);						\
      else{ if((_i>=K1)&&(_i<n-K2)){ printf("... "); _i=n-K2-1; }}	\
    printf("]\n");							\
  }while(0)

#define PAUSE scanf("%*c") /* appuyer sur [RETURN] pour continuer */
#define STRTOI(s) ((int)strtol(s,NULL,10)) 
#define STRTOL(s) ((long)strtol(s,NULL,10)) 
#define STRTOD(s) strtod(s,NULL)
#define FAIL_ROUTING -1000000 // constante suffisament <0 pour provoquer un Ã©chec de routage

/* variables globales */

int NF;          /* nb de sommets final du graphes (donc aprÃ¨s suppression) */
int *V;          /* Ã©tiquette des sommets, en principe V[i]=i */
int *VF;         /* VF[j] est l'indice i du j-Ã¨me sommet non supprimÃ© */
int *INC;        /* INC[i]=deg(i). Si =0, alors i est un sommet isolÃ© */
int ARGC;        /* variable globale pour argc */
string *ARGV;     /* variable globale pour argv */
char PARAM_PAL[64];/* mot dÃ©finissant le dÃ©gradÃ© de couleur pour -vcolor pal */
void* CPARAM=NULL; /* liste de paramÃ¨tres (pointeur tout type, en octets) pour -check */
void* FPARAM=NULL; /* liste de paramÃ¨tres (pointeur tout type, en octets) pour -filter */
int CVALUE;      /* sert pour la valeur dans -filter */
int PVALUE;      /* =1 ssi on affiche la valeur du paramÃ¨tre dans -filter */
test *FTEST;     /* pour l'option -filter */
double DELE=0;   /* proba de supprimer une arÃªtes */
double DELV=0;   /* proba de supprimer un sommet */
double REDIRECT=0; /* proba de rediriger une arÃªte */
int VERTEX0=-1;  /* voisinage du sommet Ã  afficher pour -format vertex */
int SHIFT=0;     /* dÃ©but de la numÃ©rotation des sommets */
int PERMUTE=0;   /* vrai ssi -permute */
int POS=0;       /* vrai ssi -pos */
int FAST=0;      /* vrai ssi -fast */
int CHECK=0;     /* vrai ssi option -check */
int VARIANT=0;   /* variante de l'option -variant */
int LABEL=0;     /* vrai ssi affiche les labels des sommets (-label) */
int NORM=NORM_L2;/* norme pour les graphes gÃ©omÃ©triques: L2 par dÃ©faut */
int NORM_poly=3; /* nombre de cotÃ©s de la norme polygonale */
int FORMAT=F_standard; /* type de la sortie, par dÃ©faut le format standard */
int HEADER=0;    /* par dÃ©faut pas de prÃ©ambule, sinon HEADER=1 */
int WIDTH=12;    /* nb maximum d'arÃªtes ou de sommets isolÃ©s affichÃ©s par ligne */
unsigned SEED;   /* graÃ®ne du gÃ©nÃ©rateur alÃ©atoire */
string DOTFILTER="neato"; /* nom du filtre "dot" par dÃ©faut */
string DOTSCALE=NULL; /* facteur d'Ã©chelle des points et arÃªtes du graphe (dot) */
string CAPTION=NULL; /* lÃ©gende du graphe (dot) */

string FILEXY;    /* nom de fichier des coordonnÃ©es */
double BOXX=-1,BOXY; /* pour le redimensionement: <0 signifie pas de redim. */
double XMIN=0,YMIN=0,XMAX=1,YMAX=1; /* Bounding Box par dÃ©faut */
int ROUND=DBL_DIG; /* arrondi de Q->xpos/Q->ypos Ã  10^-ROUND prÃ¨s */
/* DBL_DIG (=15) est considÃ©rÃ©e comme valeur impossible pour ROUND et signifie aucun arrondi */
double XYnoiser=-1,XYnoisep; /* paramÃ¨tres pour -xy noise: <0 signifie pas de "noise" */
int XYtype=XY_UNIF; /* type de gÃ©nÃ©ration des points, uniforme par dÃ©faut */
int XYunique=0; /* =1 ssi on Ã©limine les points doubles */
int XYgrid=0; /* <>0 si on affiche une grille grisÃ©e */
int XYpoly; /* nombre de cotÃ©s du polygone rÃ©gulier pour -xy polygon */
double XYratio=1; /* ratio pour les options -xy polar, circle, convex */
int XYzero=0; /* =1 ssi il faut ajouter le point (0,0) en rouge pour le format dot */
int XYborder=0; /* =1 ssi il faut ajouter un bord pour le format dot */
int Xmesh=0,Ymesh=0; /* dimension de la grille pour l'option -xy mesh */
double XYvsize=1; /* facteur d'Ã©chelle pour la taille des sommets dans F_dot */
int XYseedk;  /* nombre de graÃ®nes pour gÃ©nÃ©ration des points */
double XYpower; /* exposant pour la loi puissance de la gÃ©nÃ©ration des points */
double *XSEED=NULL,*YSEED=NULL; /* tableaux de doubles pour les graÃ®nes */
int XYsurfacesize=0; /* nombre de coutures de la surface pour les graphes gÃ©omÃ©triques */
int XYsurface[SURFACEMAX]; /* signature de la surface */
int LOADC=0; /* vrai ssi graphe "loadc file" */
graph* GF=NULL;     /* graphe pour l'option -check */
graph* FAMILY=NULL; /* graphe pour l'option -filter */
int HASH=H_PRIME; /* fonction de hashage par dÃ©faut */
int VSIZE=0; /* code pour la taille des sommets */
int VCOLOR=0; /* code pour la couleur les sommets */
string const COLOR_CHAR="redjykugfocatbhsqvmpinzxlw";
color COLOR_RGB[]={ /* HTML Color Names */
  /* palette des couleurs de bases, l'ordre doit Ãªtre celui de COLOR_CHAR */
  {255,  0,  0}, // r=red
  {210,105, 30}, // e=chocolate
  {255,140,  0}, // d=darkorange
  {255,165,  0}, // j=orange
  {255,255,  0}, // y=yellow
  {240,230,140}, // k=khaki
  {154,205, 50}, // u=yellowgreen
  {  0,255,  0}, // g=green (lime)
  { 34,139, 34}, // f=forestgreen
  {128,128,  0}, // o=olive
  {  0,255,255}, // c=cyan
  {127,255,212}, // a=aquamarine
  {  0,128,128}, // t=teal
  {  0,  0,255}, // b=blue
  {255,105,180}, // h=hotpink
  {250,128,114}, // s=salmon
  {255,192,203}, // q=pink
  {238,130,238}, // v=violet
  {255,  0,255}, // m=magenta
  {128,  0,128}, // p=purple
  { 75,  0,130}, // i=indigo
  {  0,  0,128}, // n=navy
  {  0,  0,  0}, // z=black
  {128,128,128}, // x=gray
  {230,230,250}, // l=lavender
  {255,255,255}  // w=white
};

/*
  Constante dÃ©finissant le nombre de couleurs de la palette de
  base. Attention!

    int const n=sizeof(T)/sizeof(*T);
    int s=n*n;
  
  ne compile pas correctement sous Linux (erreur sur la ligne "int
  s=n*n;") car n n'est pas considÃ©rÃ©e comme une constante. On a la
  mÃªme erreur que si on avait oubliÃ© le "const". #define n
  (sizeof(T)/sizeof(*T)) marcherait, sauf que la division est
  effectuÃ©e Ã  chaque fois alors que la valeur est connue Ã  la
  compilation. L'alternative est donc d'utiliser enum{
  n=sizeof(T)/sizeof(*T) }; et alors n est comme une constante.
*/
enum{ COLOR_NB=sizeof(COLOR_RGB)/sizeof(*COLOR_RGB) }; 
color *PALETTE=COLOR_RGB; /* palette par dÃ©faut */
int NPAL=COLOR_NB; /* taille de la palette par dÃ©faut */

struct{
  int mode,u,v,dist;
  double stretch;
} SCENARIO; /* scenario pour l'option "-check routing" */


/*
  Structure utilisÃ©e pour l'Ã©valuation d'un graphe et aussi le rendu
  avec Out(). Tous les calculs, notamment les entrÃ©es/sorties,
  effectuÃ©s par les fonctions d'adjacence ne devraient utiliser que
  les Ã©lÃ©ments de cette structure.
*/

typedef struct _query {
  short code; // type de requÃªte, voir enum{ QUERY_... }
  int (*adj)(struct _query* const); // fonction d'adjacence
  unsigned seed; // valeur SEED au moment oÃ¹ le graphe est initialisÃ© par QUERY_INIT
  int i,j; // sommet i et j
  int directed; // graphe orientÃ© ou non, =0 par dÃ©faut
  int loop; // supprime (=0 par dÃ©faut), autorise (=1) ou force (=2) les boucles
  int not; // complet du graphe (=1) ou pas (=0 par dÃ©faut)
  double *xpos,*ypos; // coodonnÃ©es des points pour les graphes gÃ©omÃ©triques
  
  // Les paramÃ¨tres du graphes, param[] et dparam[], sont des tableaux
  // de taille PARMAX sauf pour les graphes ayant un nombre de
  // paramÃ¨tres variables (bdrg, ...). Dans ce cas ils sont
  // rÃ©allouÃ©s. Ils seront libÃ©rÃ©s dans le free_query() final.

  int *param; // paramÃ¨tres entiers du graphes, allouÃ©/libÃ©rÃ© par new/free_query()
  double *dparam; // paramÃ¨tres rÃ©els du graphes, allouÃ©/libÃ©rÃ© par new/free_query()
  string sparam; // paramÃ¨tre chaÃ®nes de caractÃ¨res, libÃ©rer par free_query()

  /* utilisÃ©e pour une valeur de retour */

  int n; // nombre de sommet du graphe
  int a; // adjacence: est-ce que i est voisin de j ?
  char name[NAMEMAX]; // nom du sommet i
  int **rep; // pour la reprÃ©sentation implicite du graphe
  double **drep; // comme rep mais avec des doubles
  int k; // pour la taille maximum de rep[u] ou drep[u]
  int *wrap; // tableau annexe (cf. grid, rpartite, permutation,...)
  graph* G;  // graphe complet ou partiel

  /* pas encore utilisÃ©s */

  int deg; // degrÃ© du sommet i
  int *list; // liste des voisins du sommet i
  int error; // code d'erreur, 0 si tout est ok, >0 sinon
} query;


typedef int adjacence(query* const);


/* codes de requÃªte/erreur pour les fonctions d'adjacence */

enum{
  QUERY_INIT,  // rÃ©sultat dans ->n, initialise la fonction
  QUERY_END,   // termine la fonction
  QUERY_ADJ,   // rÃ©sultat dans ->a
  QUERY_NAME,  // rÃ©sultat dans ->name
  QUERY_ISOL,  // pour l'affichage de sommets isolÃ©s dans Out()

  // Ã  finir ...
  QUERY_DEG,   // rÃ©sutlat dans ->deg
  QUERY_LIST,  // rÃ©sultat dans ->list
  QUERY_DOT,   // pour le dessin des arÃªtes
  QUERY_LNAME, // pour la liste des noms des sommets
  QUERY_GRAPH, // rÃ©sultat dans ->G
};


struct{ /* pour dÃ©finir une arÃªte avec Out(i,j) sous le format F_userdot */
  adjacence *adj; // nom de la fonction d'adjacence
  int i,j;        // derniÃ¨re arÃªte calculÃ©e par adj(i,j)
  void *ptr;      // informations permettant d'afficher le dessin de l'arÃªte i-j 
} USERDOT;


/***********************************

         ROUTINES EN VRAC

***********************************/


void Erreur(int const erreur){ /* affiche l'erreur et termine avec exit() */
  string s;
  switch(erreur){
  case  1: s="option -xy non reconnue."; break;
  case  2: s="option non reconnue."; break;
  case  3: s="espace mÃ©moire insuffisant."; break;
  case  4: s="nombre trop grand de paramÃ¨tres."; break;
  case  5: s="format de sortie inconnu."; break;
  case  6: s="paramÃ¨tre incorrect."; break;
  case  7: s="ouverture du fichier impossible."; break;
  case  8: s="tableau de coordonnÃ©es inexistant."; break;
  case  9: s="option -vcolor non reconnue."; break;
  case 10: s="nom de graphe absent, inconnu ou erreur dans les paramÃ¨tres."; break;
  case 11: s="le graphe doit Ãªtre connexe."; break;
  case 12: s="option -check non reconnue."; break;
  case 13: s="format de famille de graphes invalide."; break;
  case 14: s="option -filter non reconnue."; break;
  case 15: s="graphe(s) non trouvÃ©(s)."; break;
  case 16: s="plage de valeurs incorrecte."; break;
  case 17: s="nom ou identifiant de sommets trop grand."; break;
  case 18: s="il faut Câ»Â¹ injective, soit aáµ¢Â·i + báµ¢ â‰¡ 0 (mod k) pour tous les i."; break;
  case 19: s="nom de fichier trop long."; break;
  case 20: s="nombre de couleurs dans la palette trop important."; break;
  case 21: s="code inconnue dans la fonction SortInt()."; break;
  case 22: s="sommet de valeur nÃ©gative dans le format standard."; break;
  case 23: s="code inconnu dans la fonction routing_test()."; break;
  case 24: s="option -visu incompatible avec -format no."; break;
  case 25: s="la variante -fast n'est pas implÃ©mentÃ©e pour ce graphe."; break;
  case 26: s="numÃ©ro de chronomÃ¨tre incorrect."; break;
  case 27: s="plusieurs options -check sur la ligne de commande."; break;
  case 28: s="mauvais format du fichier d'entrÃ©e."; break;
  case 29: s="loadc et -not/-permute sont incompatibles."; break;
  case 30: s="loadc devrait Ãªtre suivi d'une option -check."; break;
  case 31: s="le graphe doit Ãªtre simple et non-orientÃ©, essayez -check info."; break;
  case 32: s="problÃ¨me dans la gÃ©nÃ©ration du graphe."; break;
  case 33: s="dÃ©passement arithmÃ©tique."; break;
  case 34: s="la sÃ©quence de degrÃ© n'est pas graphique."; break;
  case 35: s="-caption ne devrait contenir qu'une occurence du format %XXX."; break;
  case 36: s="le graphe doit comporter au moins une arÃªte."; break;
  case 37: s="sommet n'appartenant pas au graphe."; break;
  case 38: s="la probabilitÃ© devrait Ãªtre dans l'intervalle [0,1]."; break;
  case 39: s="schÃ©ma de routage non reconnu."; break;
  case 40: s="fonction de hachage non reconnue."; break;
  case 41: s="scenario non reconnu."; break;
  case 42: s="nombre de cotÃ© du polygone incorrect."; break;
  case 43: s="option -norm non reconnue."; break;
  case 44: s="option dÃ©finie seulement pour les graphes gÃ©omÃ©triques."; break;
  case 45: s="signature incorrecte dans l'option -xy surface."; break;
  case 46: s="l'entier doit Ãªtre infÃ©rieur Ã  INT_MAX."; break;
  case 47: s="algorithme non implÃ©mentÃ©."; break;
  case 48: s="graphe vide ou dfs()/bfs() Ã  partir d'un sommet non valide."; break;
  case 49: s="option -dot non reconnue."; break;
  case 50: s="la variante qui s'applique (option -variant v) n'est pas valide."; break;
  default: s="code d'erreur inconnue."; /* ne devrait jamais arriver */
  }
  fprintf(stderr,"Erreur : %s\n",s);
  exit(EXIT_FAILURE);
}


long randomu(long const k)
/*
  Renvoie un entier alÃ©atoire uniforme dans [0,k[ oÃ¹ k peut Ãªtre sur
  62 bits. La solution "random()%k" possÃ¨de un biais (sauf si k est
  une puissance de deux) et n'est que sur 31 bits. La fonction
  arc4random_uniform() Ã©vite le biais et est sur 32 bits mais ne
  permet pas d'initialiser le gÃ©nÃ©rateur comme srandom(). Le code est
  inspirÃ© de l'implÃ©mentation open-source de arc4random_uniform().

  L'idÃ©e est de choisir un nombre alÃ©atoire uniforme r dans un
  intervalle de longueur multiple de k puis de renvoyer r%k, ce qui
  est uniforme. On choisit donc un intervalle de la bonne longueur
  avec [m,2^62[ oÃ¹ m=(2^62)%k. Pour choisir r uniformÃ©ment dans
  [m,2^62[ on tire r dans [0,2^62[ grÃ¢ce Ã  un double random(), RAND62,
  jusqu'Ã  obtenir r>=m. On remarque que si n>=k, alors n%k<n/2, et
  bien sÃ»r n%k<k. Donc m<2^61 et il n'y un peu moins d'une chance sur
  2 d'avoir r<m. Donc en moyenne il y a un rejet.

  Bien sÃ»r, si k est une puissance de deux constante assez petite,
  alors random()&(k-1) est la solution la plus efficace. Par exemple
  random()&3 sera est plus efficace que randomu(4).
*/
{
#define DEUX62 (0x4000000000000000L) // 2^62 = 18446744071562067968
  ASSERT(ULONG_MAX>=DEUX62-1); // un long doit contenir 62 bits au moins

  if(k<2) return 0;
  long r; long const m=DEUX62%k; // ici k>0 et m<min(k,n/2)<2^62
  do r=RAND62; while(r<m); // r = random() uniforme sur 62 bits
  return r%k;
}
#undef DEUX62


/*
  Permute alÃ©atoirement les n premiers Ã©lÃ©ments de T, et ceci
  indÃ©pendamment du type des Ã©lÃ©ments de T. Les paramÃ¨tres sont
  Ã©valuÃ©s une seule fois comme une fonction inline.
*/
#define Permute(T,n)				\
  do{						\
    int _i,_j,_n=(n);				\
    for(_i=0;_i<_n;_i++){			\
      _j=_i+randomu(_n-_i);			\
      SWAP(T[_i],T[_j]);			\
    }						\
  }while(0)


void name_vector(string const S,
		 const int* const R,int const n,
		 string const sep,string const par,int const d,
		 string const f)
/*
  Ã‰crit dans S le nom reprÃ©sentÃ© par les n premiers entiers (chiffres)
  du tableau R. Les chiffres sont sÃ©parÃ©es par la chaÃ®ne "sep" (qui
  peut Ãªtre vide pour signifier qu'il n'y a pas de sÃ©paration). Le mot
  est parenthÃ©sÃ© par la chaÃ®ne "par" qui peut soit Ãªtre vide soit
  avoir exactement deux caractÃ¨res: par[0]=parenthÃ¨se gauche
  (ouvrante), par[1]=parenthÃ¨se droite (fermante). Si d>0, les
  chiffres de R sont Ã©crits dans le sens croissant des indices (vers
  la droite), sinon c'est dans le sens des indices dÃ©croissant (sens
  normal des nombres). La chaÃ®ne f indique le format d'affichage (avec
  printf) de l'entier R[i]. En gÃ©nÃ©ral, il s'agit de "%i".

  Ex: R={3,6,1}, n=3, f="%i".
   si d=1, sep="," et par="{}" alors S="{3,6,1}"
   si d=1, sep="-" et par=""   alors S="3-6-1"
   si d=1, sep=""  et par=""   alors S="361"
   si d=0, sep=""  et par=""   alors S="163"
*/
{
  int i,b,c;
  int p=!ESTVIDE(par); /* p=0 ou 1=nombre de caractÃ¨res Ã©crits dans S */
  
  if(d>0) c=1,i=0,b=n;
  else c=-1,i=n-1,b=-1;

  /* parcoure R dans un sens ou l'autre */
  VIDE(S); /* vide la chaÃ®ne */
  while(i!=b){
    p+=sprintf(S+p,f,R[i]);
    if(p>NAMEMAX) Erreur(17);
    i += c;
    if(i!=b) p+=sprintf(S+p,"%s",sep);
  }

  /* met les parenthÃ¨ses ? */
  if(ESTVIDE(par)) return;
  S[0]=par[0];
  S[p]=par[1];
  S[p+1]='\0';
}


void name_base(string const S,
	       int u,int const b,int n,
	       string const sep,string const par,int const d)
/*
  Comme name_vector(...,n,sep,par,d,"%i") sauf que S est l'Ã©criture de
  l'entier u Ã©crit en base b. Si n<=0, alors on calcule n comme Ã©tant
  le nombre chiffres pour Ã©crire u en base b. Autrement dit, u sera
  Ã©crit dans ce cas sur un nombre variable de chiffres.

  Ex: u=361, b=10 et d=0.
   si n=3, sep="," et par="{}" alors S="{3,6,1}"
   si n=4, sep="-" et par=""   alors S="0-3-6-1"
   si n=4, sep=""  et par=""   alors S="0361"
   si n=0, sep=""  et par=""   alors S="361"
*/
{
  if(b<2) return;
  int R[NAMEMAX],i=0;
  do R[i++]=u%b, u/=b; while(u>0); // Ã©crit au moins un chiffre
  if(n<=0) n=i; // n=i=nombre chiffres Ã©crits dans R
  else for(;i<n;i++) R[i]=0; // sinon remplit le reste avec des 0
  name_vector(S,R,n,sep,par,d,"%i");
}


int LoadXY(query* const Q,string file)
/*
  Remplit les tableaux Q->xpos et Q->ypos Ã  partir d'un fichier
  (file), et renvoie le nombre de points lues. Le programme peut Ãªtre
  amenÃ© Ã  lire un point de trop (et des coordonnÃ©es vides) si le "n"
  du fichier est > aux nombres de points qu'il contient rÃ©ellement et
  que le fichier contient un retour de ligne finale.
*/
{
  FILE *f=stdin;
  int n,i;

  if(strcmp(file,"-")) f=fopen(file,"r"); /* ouvre en lecture */
  if(f==NULL) Erreur(7);
  i=fscanf(f,"%i",&n); /* lit n */
  if((n<0)||(i<0)) n=0;
  ALLOC(Q->xpos,n);
  ALLOC(Q->ypos,n);
  for(i=0;(i<n)&&(!feof(f));i++){
    if(fscanf(f,"//%*[^\n]\n")>0) continue;
    fscanf(f,"%lf %lf",Q->xpos+i,Q->ypos+i);
  }
  fclose(f);
  return i; /* i=minimum entre n et le nb de valeurs lues */
}


list *new_list(void)
/*
  CrÃ©e une nouvelle liste (qui est renvoyÃ©e) comprenant une seule
  cellule (sentinelle) dont le champs ->next=NULL. Le champs ->item
  est indÃ©terminÃ©.
*/
{
  NALLOC(list,L,1);
  L->next=NULL;
  return L;
}


static inline list *Insert(list *p,int v,int t)
/*
  Ã‰crit (v,t) dans l'Ã©lÃ©ment (->item) de la liste chaÃ®nÃ©e pointÃ©e par
  p qui ne doit pas Ãªtre NULL, puis crÃ©e une nouvelle cellule qui est
  chaÃ®nÃ©e Ã  la suite de p. On renvoit la nouvelle cellule crÃ©Ã©e. C'est
  donc un ajoÃ»t en fin de liste.
*/
{
  p->item=v;
  p->type=t;
  return p->next=new_list(); /* nouvelle cellule vide */
}


graph* new_graph(int const n)
/*
  Renvoie un objet de type "graph". Les champs sont initialisÃ©s Ã 
  leurs valeurs par dÃ©faut. Si n>0, alors les tableaux ->d et ->L de
  taille n sont allouÃ©s, mais pas les n tableaux ->L[u]. Les tableaux
  ->W, ->xpos et ->ypos ne sont pas allouÃ©s.
*/
{
  NALLOC(graph,G,1);
  G->n=0;
  G->m=-1;
  G->sort=0;
  G->id=-1;
  G->d=NULL;
  G->L=NULL;
  G->W=NULL;
  G->xpos=NULL;
  G->ypos=NULL;
  G->f=0;
  G->sym=1;
  G->G=NULL;

  G->pint1=NULL;
  G->int1=-1;

  if(n>0){
    G->n=n;
    ALLOC(G->d,n);
    ALLOC(G->L,n);
  }

  return G;
}


void free_graph(graph* const G)
/*
  LibÃ¨re G et tous ses tableaux. Dans le cas d'une famille G, chaque
  graphe est aussi libÃ©rÃ© (de maniÃ¨re rÃ©cursive). Attention! il faut
  que chaque G->L[u] et G->W[u] soit allouÃ© s'ils sont non-nuls (Ã 
  cause du FREE2).
*/
{
  if(G==NULL) return;

  /* Remarque: ce n'est pas grave de faire free() sur un ptr NULL */

  free(G->d);
  free(G->pint1);
  free(G->xpos);
  free(G->ypos);
  FREE2(G->L,G->n);
  FREE2(G->W,G->n);
  for(int i=0;i<G->f;i++) free_graph(G->G[i]);
  free(G->G);
  free(G);
}


long SizeOfGraph(const graph* G)
/*
  Donne la taille mÃ©moire (nombre d'octets) utilisÃ© par le graphe G
  (ou une famille de graphe).
*/
{
  int u;
  int const n=G->n;
  long t=sizeof(*G) + 2*n*sizeof(int); // taille de G->L et G->d
  if(G->W) t += n*sizeof(*(G->W)); // ajoute G->W
  if(G->xpos) t += n*sizeof(*(G->xpos)); // ajoute G->W
  if(G->ypos) t += n*sizeof(*(G->ypos)); // ajoute G->W
  for(u=0;u<n;u++) t += G->d[u]*sizeof(int);
  if(G->f==0) return t; // cas d'un graphe simple

  for(u=0;u<G->f;u++) t += SizeOfGraph(G->G[u]);
  return t;
}


path *new_path(graph* const G,int* const L,int k)
/*
  CrÃ©er un chemin d'un graphe G, dÃ©fini par une liste L de k sommets
  de G. Attention, le champs P du chemin renvoyÃ© est utilisÃ© en
  interne. Il ne faut pas dÃ©truire P aprÃ¨s cet appel. P sera libÃ©rÃ©
  par free_path(). Si L=NULL, alors le champs P de taille k est
  allouÃ©, et le champs n=0. C'est une faÃ§on de crÃ©er un chemin vide
  d'au plus k sommets. Le champs V, de taille G->n, est initialisÃ©
  suivant L (si L<>NULL), ou bien Ã  -1 (si L=NULL).
*/
{
  if(G==NULL) return NULL;
  
  NALLOC(path,Q,1); // Q=nouveau chemin qui sera renvoyÃ©
  ALLOCZ(Q->V,G->n,-1);
  Q->aux=NULL;

  if(L){
    int i;
    for(i=0;i<k;i++) Q->V[L[i]]=i;
    Q->P=L;
    Q->n=k;
  }
  else{
    ALLOC(Q->P,k);
    Q->n=0;
  }

  return Q;
}


void free_path(path* const P)
{
  if(P==NULL) return;
  free(P->P);
  free(P->V);
  free(P->aux);
  free(P);
}


query *new_query(void)
{
  NALLOC(query,Q,1);

  Q->error=0; // par dÃ©faut tout est ok
  Q->i=Q->j=-1;
  Q->n=-1;
  Q->a=-1;
  Q->directed=0;
  Q->loop=0;
  Q->not=0;
  Q->k=0;

  Q->rep=NULL;
  Q->xpos=Q->ypos=NULL;
  Q->sparam=NULL;
  Q->wrap=NULL;
  Q->adj=NULL;
  Q->G=NULL;
  VIDE(Q->name);

  ALLOC(Q->param,PARMAX);
  ALLOC(Q->dparam,PARMAX);

  return Q;
}


void free_query(query* const Q)
{
  if(Q==NULL) return;
  free(Q->param);
  free(Q->dparam);
  free(Q->sparam);
  free_graph(Q->G);
  free(Q);
  return;
}


/* structure de donnÃ©es pour une forÃªt enracinÃ©e */
typedef struct{
  int n;       /* nombre de sommets */
  int nroot;   /* nombre de racines, c'est-Ã -dire d'arbres de la forÃªt */
  int *lroot;  /* lroot[i]=i-Ã¨me racine de la forÃªt, i=0..nroot-1 */
  int *height; /* height[i]=hauteur du i-Ã¨me arbre, i=0..nroot-1 */
  int *root;   /* root[u]=racine de u (root[u]=u si u racine) */
  int *parent; /* parent[u]=parent de u, -1 si u racine */
  int *nchild; /* nchild[u]=nombre de fils de u */
  int **child; /* child[u][i]=i-Ã¨me fils de u, i=0..nchild[u]-1 */
  int *dfs;    /* dfs[u]=ordre dfs de u (dfs[u]=i <=> order[i]=u) */
  int *order;  /* order[i]=u si u est le i-Ã¨me sommet dans le parcours prÃ©-fixe */
  int *post;   /* post[i]=u si u est le i-Ã¨me sommet dans le parcours post-fixe */
  int *weight; /* weight[u]=nombre de descendents de u, u compris */
  int *depth;  /* depth[u]=profondeur de u dans son arbre */
  int *light;  /* light[u]=plus proche ancÃªtre (propre) lÃ©ger de u, =-1 si u racine */ 
  int *apex;   /* apex[u]=premier sommet de la branche lourde de u */
  int *heavy;  /* heavy[u]=1 ssi l'arÃªte entre u et son parent est lourde, =0 si u racine */
} tree;


tree *new_tree(int const n)
/*
  Renvoie un objet de type "tree", un arbre enracinÃ© Ã  n sommets. Les
  champs sont initialisÃ©s Ã  leurs valeurs par dÃ©faut. Si n>0, alors
  les tableaux simple de taille n sont allouÃ©s, mais pas les doubles
  tableaux comme child.
*/
{
  NALLOC(tree,T,1);
  T->n=max(n,0);
  T->nroot=-1;
  T->lroot=NULL;
  T->height=NULL;
  T->root=NULL;
  T->parent=NULL;
  T->nchild=NULL;
  T->child=NULL;
  T->dfs=NULL;
  T->order=NULL;
  T->post=NULL;
  T->weight=NULL;
  T->depth=NULL;
  T->light=NULL;
  T->apex=NULL;
  T->heavy=NULL;

  if(n>0){
    ALLOC(T->lroot,n);
    ALLOC(T->height,n);
    ALLOC(T->root,n);
    ALLOC(T->parent,n);
    ALLOC(T->nchild,n);
    ALLOC(T->child,n);
    ALLOC(T->dfs,n);
    ALLOC(T->order,n);
    ALLOC(T->post,n);
    ALLOC(T->weight,n);
    ALLOC(T->depth,n);
    ALLOC(T->light,n);
    ALLOC(T->apex,n);
    ALLOC(T->heavy,n);
  }

  return T;
}


void free_tree(tree* const T)
/*
  LibÃ¨re un arbre T et tous ses tableaux.
  Attention au FREE2.
*/
{
  if(T){
    free(T->lroot);
    free(T->height);
    free(T->root);
    free(T->parent);
    free(T->nchild);
    FREE2(T->child,T->n);
    free(T->dfs);
    free(T->order);
    free(T->post);
    free(T->weight);
    free(T->depth);
    free(T->light);
    free(T->apex);
    free(T->heavy);
    free(T);
  }
}


enum{
  TREE_PARENT_COPY = 0x0001, // duplique le tableau d'origine (sinon affecte le pointeur)

  TREE_NCHILD_FREE = 0x0002, // libÃ¨re nchild[]
  TREE_CHILD_FREE  = 0x0004, // libÃ¨re child[][]
  TREE_LROOT_FREE  = 0x0008, // libÃ¨re lroot[]

  TREE_DFS         = 0x0010, // calcule dfs[]
  TREE_ORDER       = 0x0020, // calcule order[]
  TREE_DEPTH       = 0x0040, // calcule depth[]
  TREE_HEIGHT      = 0x0080, // calcule height[]
  TREE_POST        = 0x0100, // calcule post[]

  TREE_WEIGHT      = 0x0200, // calcule weight[]
  TREE_LIGHT       = 0x0400, // calcule light[]
  TREE_APEX        = 0x0800, // calcule apex[]
  TREE_HEAVY       = 0x1000, // calcule heavy[]
};


tree *MakeTree(int* const P,int const n,unsigned const code)
/*

  NON UTILISEE, NON TESTEE

  Construit une forÃªt enracinÃ©e (structure tree) Ã  partir d'une
  relation de parentÃ©e P Ã  n sommets (tableau d'entiers P[u]=pÃ¨re(u)
  ou -1 s'il n'en a pas). Certaines tables de la structure sont
  initialisÃ©es ou pas suivant la valeur binaire de code. Pour tester
  les bits de code il faut utiliser l'enum TREE_xxx ci-dessus (mask
  &). Les complexitÃ©s en temps et en espace sont en O(n).

  o Plus prÃ©cisÃ©ment, les tables parent[], nchild[], child[][] et
    lroot[] sont toujours calculÃ©es. Elles peuvent Ãªtre libÃ©rÃ©es
    suivant les bits de code, sauf parent[].

  o Les tables dfs[], order[], depth[], height[] et post[] sont
    calculÃ©es ou pas suivant les bits de code, mais si l'une d'elles
    est calculÃ©e alors toutes le sont (car c'est un mÃªme parcours).

  o Enfin, les tables weight[], light[], apex[] et heavy[] peuvent
    Ãªtre calculÃ©es ou pas suivant les bits de code. Ces derniÃ¨res
    tables vont calculer de maniÃ¨re intermÃ©diaire dfs[], order[],
    depth[], height[] et post[] puis les supprimer (Ã©ventuellement).

  Ã€ part les tables toujours calculÃ©es, toutes les tables
  intermÃ©diaires calculÃ©es sont supprimÃ©es sauf si le bit de code
  correspondant indique le contraire.
*/
{
  if(n<=0) return NULL; /* rien Ã  faire ou problÃ¨me */

  tree *T=new_tree(n); /* T->n=n>0 */
  int u,p,i;

  /* copie le tableau P ou pas */
  
  if(code&TREE_PARENT_COPY) ALLOCZ(T->parent,n,P[_i]);
  else T->parent=P;

  /* calcule le nombre de racines et le nombre de fils pour chaque
     sommets: T->nroot, T->nchild, T->root */

  ALLOCZ(T->nchild,n,0);
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=pÃ¨re(u) */
    if(p<0){
      T->root[u]=u;
      T->nroot++;
    }else T->nchild[u]++;
  }

  /* alloue les tableaux T->child et remet Ã  zÃ©ro T->nchild */

  for(u=0;u<n;u++)
    if(T->nchild[u]){
      ALLOC(T->child[u],T->nchild[u]);
      T->nchild[u]=0;
    }else T->child[u]=NULL;
  
  /* remplit les tableaux T->child et rÃ©tablit T->nchild */
  /* remplit aussi la liste des racines T->lroot */

  ALLOC(T->lroot,T->nroot);
  T->nroot=0; /* on va le recalculer */
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=pÃ¨re(u) */
    if(p<0) T->lroot[T->nroot++]=u; /* ajoute u aux racines */
    else T->child[p][T->nchild[p]++]=p; /* ajoute u aux fils de p */
  }

  /* parcoure la forÃªt si nÃ©cessaire */
  /* calcule T->dfs, T->order, T->depth, T->height, T->post */

  if((code>>4)==0) goto maketree_fin;
  // vrai ssi l'un des bits de code aprÃ¨s le 4e est mis
  // NB: <=> code&(TREE_DFS|TREE_ORDER|TREE_DEPTH|...|TREE_HEAVY)

  NALLOC(int,pile,n); /* pile */
  NALLOC(int,next,n); /* next[u]=compteur courant du prochain voisin de u Ã  visiter */
  int sp=-1; /* sp=sommet de la pile, pile[sp]=dernier Ã©lÃ©ment empilÃ© */
  int v,dfs=0; /* dfs=date de premiÃ¨re visite */
  p=0; /* ordre post-fixe */

  for(i=0;i<T->nroot;i++){ /* pour chaque racine faire ... */

    u=T->root[i];
    pile[++sp]=u; /* empile la racine i */
    next[u]=0; /* 1er voisin de u Ã  visiter */
    T->dfs[u]=dfs; /* un sommet visitÃ© */
    T->order[dfs++]=u; /* ordre prÃ©-fixe des sommets */
    T->depth[u]=0; /* profondeur du sommet u */
    T->height[i]=0; /* hauteur de la i-Ã¨me racine */

    while(sp>=0){ /* tant que la pile n'est pas vide */
      u=pile[sp]; /* u=sommet courant sur la pile */
      if(next[u]<T->nchild[u]){ /* on visite le voisin de u d'indice next[u] */
	v=T->child[u][next[u]++]; /* v=voisin de u Ã  empiler */
	pile[++sp]=v; /* on empile le voisin */
	T->dfs[v]=dfs; /* date de premiÃ¨re visite de v */
	T->order[dfs++]=v; /* ordre du sommet */
	T->depth[v]=T->depth[u]+1; /* hauteur de v */
	T->height[i]=max(T->height[i],T->depth[u]);
      }else{ /* on a visitÃ© tous les voisins de u */
	T->post[p++]=pile[sp--]; /* on dÃ©pile u, ordre post-fixe des sommets */
      }
    }
    
  }
  
  free(pile);
  free(next);

  /* calcule le poids des sommets */

  if(code&(TREE_WEIGHT|TREE_LIGHT)){
    ALLOCZ(T->weight,n,1); /* tout le monde a poids 1 au dÃ©part */
    for(i=0;i<n;i++){
      u=T->post[i]; /* parcours post-fixe */
      p=T->parent[u]; /* p=pÃ¨re de u */
      if(p>=0) T->weight[p] += T->weight[u]; /* le pÃ¨re reÃ§oit le poids de son fils */
    }
  }

  /* calcule l'ancÃªtre lÃ©ger de chaque sommets (il faut les poids) */

  if(code&TREE_LIGHT){ // Ã  finir
    ALLOC(T->light,n);
    for(i=0;i<n;i++){
      u=T->order[i]; /* parcours prÃ©fixe */
      p=T->parent[u]; /* p=pÃ¨re de u */
    }
  }

  /* calcule apex[] */

  if(code&TREE_APEX){ // Ã  finir
    ALLOC(T->apex,n);
    ;
  }

  /* calcule heavy[] */

  if(code&TREE_HEAVY){ // Ã  finir
    ALLOC(T->heavy,n);
    ;
  }

  // libÃ©rations si nÃ©cesaire

  if(!(code&TREE_DFS))   { free(T->dfs);      T->dfs   =NULL; }
  if(!(code&TREE_ORDER)) { free(T->order);    T->order =NULL; }
  if(!(code&TREE_DEPTH)) { free(T->depth);    T->depth =NULL; }
  if(!(code&TREE_HEIGHT)){ free(T->height);   T->height=NULL; }
  if(!(code&TREE_POST))  { free(T->post);     T->post  =NULL; }
  if(!(code&TREE_WEIGHT)){ free(T->weight);   T->weight=NULL; }
  if(!(code&TREE_LIGHT)) { free(T->light);    T->light =NULL; }
  if(!(code&TREE_APEX))  { free(T->apex);     T->apex  =NULL; }
  if(!(code&TREE_HEAVY)) { free(T->heavy);    T->heavy =NULL; }

 maketree_fin:

  // libÃ©rations si nÃ©cesaire
  
  if(code&TREE_NCHILD_FREE){ free(T->nchild);   T->nchild=NULL; }
  if(code&TREE_CHILD_FREE) { FREE2(T->child,n); T->child =NULL; }
  if(code&TREE_LROOT_FREE) { free(T->lroot);    T->lroot =NULL; }

  return T;
}


/*****************************************************

      Fonctions de comparaisons pour les tris

  Attention ! return x-y; ne marche pas toujours, mÃªme pour comparer
  des entiers. Par exemple, x=INT_MAX et y=INT_MIN provoque un
  dÃ©passement avec x-y. Il est prÃ©fÃ©rable d'utiliser return (x>y) -
  (x<y); qui marche toujours. Pour les double/float, (int)(x-y) ne
  marchent pas toujours si x,y âˆˆ [-1,+1].

*****************************************************/


int fcmp_int(const void *P,const void *Q)
/* Compare deux entiers, pour qsort(). */
{
  int const p=*(int*)P;
  int const q=*(int*)Q;
  return (p>q) - (p<q);
}


int fcmp_int_inv(const void *P,const void *Q)
/* Comme fcmp_int(), mais dans l'ordre inverse. */
{
  return fcmp_int(Q,P);
}


int fcmp_double(const void *P,const void *Q)
/* Compare deux doubles, pour qsort(). */
{
  double const p=*(double*)P;
  double const q=*(double*)Q;
  return (p>q) - (p<q);
}


int fcmp_point(const void *P,const void *Q)
/* Compare deux points, pour qsort(). */
{
  if(((point*)P)->x < ((point*)Q)->x) return -1;
  if(((point*)P)->x > ((point*)Q)->x) return 1;
  if(((point*)P)->y < ((point*)Q)->y) return -1;
  if(((point*)P)->y > ((point*)Q)->y) return 1;
  return 0;
}


int fcmp_profile(const void *P,const void *Q)
/*
  Compare deux profiles, pour qsort(). Les profiles de plus grande
  longueur sont classÃ©s avant les plus courts, ceux-ci Ã©tant plus
  discriminant.
*/
{
  int* const A=*(int**)P;
  int* const B=*(int**)Q;

  if(*A>*B) return -1; // si longueur(A)>longueur(B), alors A<B
  if(*A<*B) return 1; // si longueur(A)<longueur(B), alors A>B
  /* ici, profiles de mÃªme longueur n=A[0] */

  int const n=*A; // surtout ne pas utiliser A[1]

  for(int u=2;u<n;u++){
    if(A[u]<B[u]) return -1;
    if(A[u]>B[u]) return 1;
  }

  return 0;
}


int fcmp_graphid(const void *P,const void *Q)
/*
  Compare les identifiants de deux graphes. Sert pour qsort() et
  bsearch(). Ici, P et Q sont des pointeurs de (graph*).
*/
{
  int const p=(*(graph**)P)->id;
  int const q=(*(graph**)Q)->id;
  return (p>q) - (p<q);
}


int fcmp_stretch(const void *P,const void *Q)
/*
  Compare les ratios P.x/P.y et Q.x/Q.y, dans le cas de ratios
  irrÃ©ductibles. Sert pour routing_test().
*/
{
  int const p=((triplet*)P)->x*((triplet*)Q)->y;
  int const q=((triplet*)P)->y*((triplet*)Q)->x;
  return (p>q) - (p<q);
}


int fcmp_tabint(const void *P,const void *Q)
/*
  Compare T[*P] avec T[*Q] oÃ¹ T est un tableau d'entiers qui doit Ãªtre
  initialisÃ© lors d'un premier appel avec fcmp_tabint(NULL,T).
  Attention ! les valeurs du tableau initial Ã  trier (*P et *Q) doivent
  Ãªtre des indices dans T, donc dans [0,|T|[.
*/
{
  static int *T; // tableau global mais local Ã  la fcmp_tabint()
  if(P==NULL){ T=(int*)Q; return 0; } // fixe le tableau T
  
  int const p=T[*(int*)P];
  int const q=T[*(int*)Q];
  return (p>q) - (p<q);
}


int fcmp_tabinteq(const void *P,const void *Q)
/*
  Comme fcmp_tabint(P,Q) qu'en cas d'Ã©galitÃ© (soit si *P==*Q) on
  renvoie la comparaison entre *P et *Q et non 0. Donc le tri
  s'effectue d'abord selon T[], puis selon les indices de T en cas
  d'Ã©galitÃ©.
*/
{
  static int *T; // tableau global mais local Ã  la fcmp_tabint()
  if(P==NULL){ T=(int*)Q; return 0; } // fixe le tableau T
  
  int const p=*(int*)P;
  int const q=*(int*)Q;
  if(T[p]==T[q]) return p-q; 
  return T[p]-T[q];
}


int fcmp_edge(const void *e1, const void *e2)
/* Comparaison du poids de deux arÃªtes pour qsort(). */
{
  double const x = ((edge*)e1)->w;
  double const y = ((edge*)e2)->w;
  return (x>y) - (x<y);
}


/* code pour ReadRange() et InRange() */
enum{
  R_EQ,   // code =x
  R_INF,  // code <x
  R_SUP,  // code >x
  R_INT,  // code x-y
  R_TRUE  // code t
};


int ReadRange(string const s,int *R)
/*
  Lit une chaÃ®ne de caractÃ¨res dÃ©crivant un intervalle de valeurs
  entiÃ¨res, et renvoie dans le tableau d'entiers R les valeurs et les
  codes correspondant pour que la fonction InRange(x,R) puisse
  fonctionner. En quelque sorte cette fonction prÃ©pare la fonction
  InRange(). On ignore les caractÃ¨res non reconnus (pas d'erreur). On
  renvoie le nombre d'opÃ©rations dÃ©codÃ©es, c'est-Ã -dire le nombre de
  valeurs Ã©crites dans le tableau R, nombre qui est aussi Ã©crit dans
  R[0].

  Ex: s="1,-3,5-7,>30,<50" (on interprÃ¨te les "," comme des "ou")
  => R={12,R_EQ,1,R_EQ,-3,R_INT,5,7,R_SUP,30,R_INF,50} (12=taille(R))

  La valeur des codes d'opÃ©rations (R_EQ, R_INF, ...) est donnÃ©e par
  l'enum ci-dessus.
*/
{
  if(R==NULL) return 0;
  if(s==NULL){ R[0]=R_EQ; return 0; }

  int i,r,p,x,start,c;
  i=x=c=0;
  r=start=p=1;

  /* r=indice de R[] */
  /* i=indice de s[] */
  /* x=valeur entiÃ¨re lue */
  /* c=1 ssi le code d'opÃ©ration a Ã©tÃ© dÃ©tectÃ© */
  /* start=1 ssi on va commencer Ã  lire un entier */
  /* p=1 ou -1, signe de x */

  while(s[i]!='\0'){
    if(s[i]=='='){ R[r++]=R_EQ; c=start=p=1; }
    if(s[i]=='<'){ R[r++]=R_INF; c=start=p=1; }
    if(s[i]=='>'){ R[r++]=R_SUP; c=start=p=1; }
    if(s[i]=='-'){
      if(start) p=-p;
      else{ R[r++]=R_INT; R[r++]=x; c=start=p=1; }
    }
    if(s[i]=='t'){ x=R_TRUE; c=r=1; break; } /* t=true, pour avoir false faire "not" et "t" */
    if(s[i]=='p') PVALUE=1; /* pas de code pour "p" */
    if(s[i]==','){
      if(c==0) R[r++]=R_EQ; /* code '=' par dÃ©faut */
      R[r++]=x; c=0; start=p=1;
    }
    if(('0'<=s[i])&&(s[i]<='9')){
      if(start) start=0;
      x=x*10+p*(s[i]-'0'); /* on lit x en base 10 en tenant compte du signe p */
    }
    if(start) x=0;
    i++;
  }

  if(PVALUE==i){ x=R_TRUE;c=1; } /* si s="p", alors comme "t" */
  if(c==0) R[r++]=R_EQ;
  R[r++]=x; /* on ajoute le dernier opÃ©rande */
  R[0]=r;
  return r;
}


int InRange(int const x,const int* const R)
/*
  DÃ©termine si x appartient aux valeurs dÃ©crites par le "range" R.
  R[0] est la taille de R, R[0] compris.
*/
{
  int i,n,t;
  n=R[t=0]; /* n=taille(R) */
  i=1; /* commence Ã  lire R[1] */
  CVALUE=x;

  while(i<n){
    switch(R[i++]){ /* lit le code d'opÃ©ration */
    case R_EQ  : t=(x==R[i++]); break;
    case R_INF : t=(x<R[i++]); break;
    case R_SUP : t=(x>R[i++]); break;
    case R_INT : t=((R[i]<=x)&&(x<=R[i+1])); i+=2; break;
    case R_TRUE: return 1;
    default: Erreur(16); /* ne devrait jamais se produire */
    }
    if(t) break;
  }
  return t;
}


/***********************************

           ROUTINES SUR
           LES GRAPHES

***********************************/


static inline void degres_zero(graph* const G)
/*
  Met Ã  zÃ©ro tous les degrÃ©s d'un graphe. C'est trÃ¨s utile pour se
  servir par exemple des macros ADD_EDGE() et ADD_ARC().
*/
{
  memset(G->d,0,G->n*sizeof(int));
}


int nb_edges(graph* const G)
/*
  Retourne le nombre d'arÃªtes d'un graphe symÃ©trique G ou bien le
  champs G->m s'il est positif. Si G->m<0, alors G->m est mis Ã  jour Ã 
  partir de la somme des G->d[i].
*/
{
  int m=G->m;
  if(m<0){
    int i;
    int const n=G->n;
    for(i=m=0;i<n;i++) m += G->d[i];
    G->m=(m>>=1);
  }
  return m;
}


int Degree(const graph* const G,int const max)
/*
  Renvoie le degrÃ© maximum (si max=1) ou minimum (si max=0) d'un
  graphe G. On renvoie -1 si G est NULL, n'a pas de sommet ou est une
  famille de graphes.
*/
{
  if((G==NULL)||(G->n<=0)||(G->f>0)) return -1;
  int const n=G->n;
  int i=1,d=G->d[0];
  if(max) for(;i<n;i++) d=max(d,G->d[i]);
  else for(;i<n;i++) d=min(d,G->d[i]);
  return d;
}


void PrintGraphList(const graph* const G)
/*
  Affiche le graphe G sous la forme d'une liste d'adjacence. Tient
  compte de SHIFT et de VERTEX0.
*/
{
  if(G==NULL){ printf("NULL\n"); return; }
  int const n=(VERTEX0<0)?G->n:VERTEX0+1;
  int u,d,i;

  for(u=(VERTEX0<0)?0:VERTEX0;u<n;u++){
    printf("%i:",u+SHIFT);
    for(i=0,d=G->d[u];i<d;i++){
      printf(" %i",G->L[u][i]+SHIFT);
    }
    printf("\n");
  }
  return;
}


void PrintGraphMatrix(const graph* const G)
/*
  Affiche le graphe G sous la forme d'une matrice d'adjacence complÃ¨te
  ou triangulaire supÃ©rieure (en tennant compte du FORMAT, smatrix ou
  matrix). La complexitÃ© en espace est seulement de O(n).
*/
{
  int u,d,i,z,t;
  int const n=G->n;

  NALLOCZ(int,M,n,0);
  t=(FORMAT==F_smatrix);

  for(u=z=0;u<n;u++){
    if(t) z=u;
    for(i=0,d=G->d[u];i<d;M[G->L[u][i++]]=1);
    for(i=0;i<n;i++)
      if(i<z) printf(" ");
      else printf("%c",'0'+M[i]);
    for(i=0;i<d;M[G->L[u][i++]]=0); /* remet rapidement M[] tout Ã  0 */
    printf("\n");
  }

  free(M);
  return;
}


void PrintPath(const graph* const G,const path* const P)
/*
  Affiche le chemin P d'un graphe G.
  Sert pour le dÃ©bugage.
*/
{
  if((G==NULL)||(P==NULL))
    printf("NULL\n");
  else{
    int i,j,u,d;
    for(i=0;i<P->n;i++)
      if(P->V[P->P[i]]!=i) break;
    if(i<P->n) goto error;
    for(u=0;u<G->n;u++)
      if((P->V[u]>=0)&&(P->P[P->V[u]]!=u)) break;
    if(u<G->n) goto error;
    printf("P->aux:");
    if(P->aux==NULL) printf(" NULL\n");
    else{
      printf("\n");
      for(i=0;i<P->n;i++){
	u=P->P[i];
	d=P->aux[u][0];
	printf("  %i:",u);
	for(j=1;j<=d;j++){
	  printf(" %i",P->aux[u][j]);
	}
	printf("\n");
      }
    }
  }
  return;
  
 error:
  printf("Chemin incohÃ©rent.\n");
  return;
}


int *SortGraph(graph* const G,int const code)
/*
  Force le tri (mÃªme si G->sort=1) des listes d'adjacence d'un graphe
  G, c'est-Ã -dire pour chaque sommet u, G->L[u] est une liste
  d'entiers triÃ©s par ordre croissant. Le champs G->sort est mis Ã 
  jour.  L'algorithme effectue un simple appel Ã  qsort(). Sa
  complexitÃ© est Ã  peu prÃ¨s en O(n+m*log(m/n)).

  Si code=0, on s'arrÃªte aprÃ¨s l'Ã©tape du tri. Sinon, on lance une
  Ã©tape de vÃ©rification du graphe: prÃ©sence de multi-arÃªtes, de
  boucles, etc.

  Le temps de la vÃ©rification est comparable Ã  celui du tri. Le
  rÃ©sultat de la vÃ©rification est un tableau de statistiques S de
  taille fixe (dÃ©clarÃ© en static qui ne doit pas Ãªtre libÃ©rÃ© par
  l'appelant), ayant la signification suivante:

    S[0]=nombre de boucles
    S[1]=nombre de multi-arcs
    S[2]=nombre d'arcs (avec multi-arcs et boucles)
    S[3]=nombre d'adjacence non-symÃ©triques
    S[4]=nombre de voisins d'ID < 0
    S[5]=nombre de voisins d'ID â‰¥ n
    S[6]=1 ssi G est simple et non-orientÃ©
    S[7]=degrÃ© maximum
    S[8]=degrÃ© minimum
    S[9]=nombre de sommets isolÃ©s

  Ã€ l'issue de la vÃ©rification, G->sym est mise Ã  jour.
*/
{
  if(G==NULL) return NULL;
  int const n=G->n;
  int u;

  /* trie G */
  for(u=0;u<n;u++)
    QSORT(G->L[u],G->d[u],fcmp_int);

  G->sort=1;
  if((code==0)||(n==0)) return NULL;
  
  /* statistiques sur G */
  static int S[10]; /* static est important, car on fait return S */
  int v,i,d,w;
  
  memset(S,0,sizeof(S)); /* initialise les stats Ã  0, NB: sizeof(S)=40 */
  S[7]=S[8]=G->d[0]; /* il faut G->d<>NULL */
  
  for(u=0;u<n;u++){ /* parcoure le graphe */
    d=G->d[u]; S[7]=max(S[7],d); S[8]=min(S[8],d);
    S[9] += (d==0); /* un sommet isolÃ© */
    S[2]+=d; /* ajoute le nombre de voisins */
    w=-1; /* w=voisin prÃ©cÃ©dant le voisin courant v */
    for(i=0;i<d;i++){ /* pour chaque voisin */
      v=G->L[u][i];
      if(u==v) S[0]++; /* une boucle */
      if(v==w) S[1]++; /* une multi-arÃªte */
      w=v; /* mÃ©morise le dernier voisin rencontrÃ© */
      if(v<0){ S[4]++; continue; } /* un voisin nÃ©gatif */
      if(v>=n){ S[5]++; continue; } /* un voisin trop grand */
      if(bsearch(&u,G->L[v],G->d[v],sizeof(int),fcmp_int)==NULL) S[3]++; /* un arc asymÃ©trique */
    }
  }

  S[6]=((S[0]+S[1]+S[3]+S[4]+S[5])==0); /* vrai ssi G simple et non-orientÃ© */
  G->sym=(S[3]==0);
  return S;
}


void PrintGraph(graph* const G)
/*
  Affiche un graphe ou une famille de graphes au format standard sous
  forme compacte. Utilise WIDTH. Effet de bord: le (ou les) graphes
  sont triÃ©s par ordre croissant, et donc G->sort=1 en sortie. Si le
  graphe est asymÃ©trique, des sommets peuvent Ãªtre affichÃ©s comme
  sommets isolÃ©s alors qu'ils ne le sont pas.
*/
{
  if(G==NULL){ printf("NULL\n"); return; }

  int u,v,i,k,n,ligne,nk=(G->f>0)?G->f:1;
  graph* H;
  int *P;

  for(k=0;k<nk;k++){

    if(G->f>0){
      H=G->G[k];
      printf("[%i]",H->id);
    }else H=G;

    SortGraph(H,0); // ordre croissant
    n=H->n;
    ALLOCZ(P,n,0);
    i=u=ligne=0;
    v=-1;

    while(i<n){
      /* si u==i, alors u=tÃªte d'un nouveau chemin */
      while((v<u)&&(P[u]<H->d[u])) v=H->L[u][P[u]++];
      if(v<u){ /* on a pas trouvÃ© v => fin d'un chemin */
	if(H->d[u]==0){ /* cas des sommets isolÃ©s */
	  printf(" %i",u);
	  if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	}
	u=(i+=(u==i));
	v=-1;
      }
      else{ /* u a un voisin v>u */
	if((u==i)||(ligne==0)) printf(" %i",u); /* on affiche la tÃªte */
	printf("-%s%i",(H->sym)?"":">",v); /* on affiche -v ou ->v */
	if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	u=v; /* on continu avec v */
	v=-1;
      }
    } /* fin du while */

    if(ligne>0) printf("\n"); /* newline si fini avant la fin de ligne */
    free(P);
  }

  G->sort=1; /* effet de bord */
  return;
}


void GraphRealloc(graph* const G,const int* const D)
/*
  Redimensionne le graphe G Ã  G->n sommets suivant le tableau de degrÃ©
  D. On rÃ©ajuste en premier les tableaux G->d et G->L pour qu'ils
  aient une taille G->n, puis on rÃ©ajuste les listes d'adjacences des
  sommets de G suivant le tableau des degrÃ©s D (qui doit Ãªtre de
  taille au moins G->n). Si D[u] est plus petit que G->d[u], alors la
  liste G->L[u] est tronquÃ©e. Si D[u] est plus grand que G->d[u],
  alors G->L[u] est rÃ©ajustÃ©. Le degrÃ© G->d[u] est initialisÃ© au
  minimum de G->d[u] et D[u]. NB: le nombre d'arÃªtes G->m, qui a pu
  changer, est rÃ©initialisÃ© Ã  -1. G->sort n'est pas changÃ© car l'ordre
  des listes G->L n'est pas modifiÃ©.

  Pour plonger G dans un graphe complet faire:
    NALLOCZ(int,D,G->n,G->n-1);
    GraphRealloc(G,D);
    free(D);
*/
{
  int const n=G->n;
  int u,d;
  for(u=0;u<n;u++){
    d=D[u];
    REALLOC(G->L[u],d);
    G->d[u]=min(G->d[u],d);
  }

  /* Il ne faut pas rÃ©ajuster G->d et G->L avant la boucle for(u=...)
     car potentiellement on libÃ¨re G->d et G->L. Or il est possible
     d'avoir D=G->d. */

  REALLOC(G->d,n);
  REALLOC(G->L,n);
  G->m=-1; /* le nombre d'arÃªtes n'est plus Ã  jour */
  return;
}


graph* new_subgraph(const graph* const G)
/*
  Renvoie un nouveau graphe R vide de maniÃ¨re similaire Ã 
  R=new_graph(G->n), mais en plus dimensionne chaque R->L[u] Ã 
  G->d[u]. On initialise aussi R->d[u] Ã  0. Le graphe renvoyÃ© est donc
  un sous-graphe couvrant de G sans aucune arÃªte.
*/
{
  if((G==NULL)||(G->n<0)) return NULL;

  const int n=G->n;
  graph* R=new_graph(n);
  int u;

  for(u=0;u<n;u++) ALLOC(R->L[u],G->d[u]);
  degres_zero(R);

  return R;
}


graph* new_fullgraph(int const n)
/*
  Renvoie un graphe G comme new_graph(n), mais en plus alloue G->L[u]
  de taille max{n-1,1}, et initialise G->d[u]=0 pour tous les sommets
  u. Une fois le graphe construit, on peut rÃ©dimensionner le graphe
  grÃ¢ce Ã  GraphRealloc, comme dans l'exemple:

    graph* G=new_fullgraph(n);
      ...
      ADD_EDGE(G,u1,v1);
      ADD_EDGE(G,u2,v2);
      ...
    GraphRealloc(G,G->d);
      ...
    free_graph(G);
*/
{
  if(n<1) return NULL;
  graph* const G=new_graph(n);
  int const n1=max(n-1,1);
  for(int u=0;u<n;u++) ALLOC(G->L[u],n1);
  degres_zero(G);
  return G;
}


graph* ExtractSubgraph(const graph* const G,const int* const T,
		       int const k,int const code)
/*
  Construit, Ã  partir d'un graphe G et d'une liste T de k sommets, un
  nouveau graphe S, renvoyÃ© par la fonction, correspondant au
  sous-graphe de G induit par les sommets de T (si code=1) ou de G\T
  (si code=0). Les sommets du graphe S sont dans [0,k[ (ou [0,n-k[ si
  code=0).

  On peut ainsi faire une copie C du graphe G simplement en faisant:

    graph* C=ExtractSubgraph(G,NULL,0,0);

  Effet de bord: S->pint1 est allouÃ© si T<>NULL. Dans ce cas on
  renvoie dans S->pint1 un tableau X de taille G->n indiquant la
  renumÃ©rotation de G: pour tout sommet u de G (u dans [0,G->n[)
  S->pint1[u]=0 si u est un sommet abscent de S et S->pint1[u]=d>0 si
  u est numÃ©rotÃ© d-1>=0 dans S. Le nombre d'arÃªtes S->m du graphe S
  renvoyÃ© est Ã  jour. L'ordre relatif des listes de G est prÃ©servÃ©. En
  particulier, si G->sort>0, alors le sous-graphe sera aussi
  triÃ©. G->sym est aussi copiÃ©.
*/
{
  if(G==NULL) return NULL;
  int const n=G->n;
  int u,v,d,i,s,ns,m;

  NALLOC(int,X,n);
  for(u=1-code,i=0;i<n;i++) X[i]=u;
  if(T) for(i=0;i<k;i++) X[T[i]] ^= 1;
  for(i=d=0;i<n;i++) if(X[i]) X[i]=++d; 
  /* ici X[i]=0 si i est un sommet Ã  supprimer */
  /* ici X[i]=d (>0) si i doit Ãªtre renumÃ©rotÃ© en d-1>=0 */

  ns=(code)?k:n-k;
  graph* S=new_fullgraph(ns);

  for(s=u=m=0;u<n;u++)
    if(X[u]){ /* si u existe, s=X[u]-1 */
      d=G->d[u];
      for(i=0;i<d;i++){
	v=G->L[u][i];
	if(X[v]){ m++; ADD_ARC(S,s,X[v]-1); } /* si v existe */
      }
      s++;
    }

  /* rÃ©duit la taille des listes */
  GraphRealloc(S,S->d);

  S->pint1=X;
  S->sort=G->sort;
  S->sym=G->sym;
  S->m=(m>>1);
  return S;
}


static inline graph* GraphCopy(const graph* const H)
/* Renvoie une copie du graphe H */
{
  return ExtractSubgraph(H,NULL,0,0);
}


graph* List2Graph(list* L,int const code)
/*
  Retourne un graphe G simple Ã  partir d'un graphe dÃ©fini par une
  liste L de codes (voir File2List() pour le codage prÃ©cis du type
  "list"). Certaines opÃ©rations sont effectuÃ©es sur L en fonction de
  la valeur binaire de code:

  - code&1 =1: optimisation des listes du graphe (tri par ordre croissant)
           =0: sans optimisation
  - code&2 =1: auto-dÃ©tection du shift dans L (pour "load file")
           =0: pas d'auto-dÃ©tection du shift
  - code&4 =1: gestion d'un sous-graphe (V,NF) => code&2=0
           =0: pas de sous-graphe

  Les codes suivants servent Ã  List2Family():

  - code&8 =1: tri de la famille suivant les identifiants (sert pour List2Family)
           =0: pas de tri de la famille
  - code&16=1: ne libÃ¨re pas la liste L (sert pour List2Family)
           =0: libÃ¨re la liste L
  - code&32=1: renvoie toujours un graphe, le 1er si c'est une famille
           =0: renvoie une famille si c'est une famille

  Pour calculer le graphe (et sa liste d'adjacence) on effectue
  plusieurs passes sur L: une passe pour dÃ©terminer n; une autre pour
  calculer les degrÃ©s des sommets; et une 3e pour remplir G et pour
  Ã©ventuellement libÃ©rer la liste L.
*/
{
  if(L==NULL) return NULL; /* si pas de cellule, ne rien faire */

  int u,v,x,n,*D;
  graph* G;
  list* p;

  u=INT_MAX;
  if(code&4){ /* si sous-graphe dÃ©finit par (V,NF) */
    p=L;
    while(p){
      p->item=V[p->item]-SHIFT;
      p=p->next;
    }
    n=NF; /* on connaÃ®t n */
  }
  else{ /* sinon, on calcule n, et on lit les valeurs min (=u) et
	   valeur max (=v) de L */
    p=L; v=0;
    while(p){
      x=p->item;
      if(x<u) u=x;
      if(x>v) v=x;
      p=p->next;
    }
    if(code&2){ /* on dÃ©cale les valeurs dans [0,n[ */
      p=L;
      while(p){
	p->item -= u;
	p=p->next;
      }
      n=v+1-u;
    }else{
      if((u<0)||(v<0)) Erreur(22); /* il ne devrait pas avoir ici de valeur < 0 */
      n=v+1;
    }
  }

  ALLOCZ(D,n,0);

  /* on lit les degrÃ©s (sortant) des sommets, et les met dans le
     tableau D. NB: la variable u n'est pas initialisÃ©, car on passe
     toujours d'abord par un item de type T_NODE */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { D[u]++; D[v]++; }      /* u-v */
    else if(p->type==T_OPENE){ D[u]++; D[v]++; x=0; } /* u-(v */
    else if(p->type==T_ARC)    D[u]++;                /* u->v */
    else if(p->type==T_OPENA){ D[u]++; x=0; }         /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* initialise la liste d'adjacence G. On se sert plus tard de D[u]
     pour indiquer la prochaine position libre dans G[u][]. */

  G=new_graph(n); /* G->n=n, alloue G->d et G->L */
  for(u=0;u<n;u++){
    ALLOC(G->L[u],D[u]); /* alloue une liste pour chaque sommet */
    G->d[u]=D[u]; /* G->d[u]=deg(u) */
    D[u]=0; /* prochaine position libre dans G[u] */
  }

  /* Remplit G. On met aussi Ã  jour G->sym (orientÃ© ou pas). On
     pourrait tester Ã  la volÃ©e si les listes sont triÃ©es et mettre Ã 
     jour G->sort. */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; }      /* u-v */
    else if(p->type==T_OPENE){ G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; x=0; } /* u-(v */
    else if(p->type==T_ARC)  { G->L[u][D[u]++]=v; G->sym=0; }               /* u->v */
    else if(p->type==T_OPENA){ G->L[u][D[u]++]=v; G->sym=0; x=0; }          /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* libÃ¨re L si bit-4 Ã  1 */
  if(!(code&16)){
    p=L;
    while(p){
      L=p;
      p=p->next;
      free(L);
    }
  }
  
  free(D); /* plus besoin de D */
  if(code&1) SortGraph(G,0);
  return G;
}


graph* List2Family(list *L,int const code)
/*
  Transforme une liste en famille de graphes.  Si L reprÃ©sente un
  graphe simple (pas de type T_NB ou T_ID), alors un graphe simple est
  retournÃ©e (plutÃ´t qu'une famille Ã  un seul Ã©lÃ©ment). Donc,
  List2Family() gÃ©nÃ©ralise List2Graph(). On utilise List2Graph() comme
  sous-routine. Pour "code" voir List2Graph().
    
  Effet de bord:
  - la famille est triÃ©e par ID croissant si code&8=1
  - la liste L est libÃ©rÃ©e si code&16=0
  - on retourne un graphe si code&32=1 (plutÃ´t qu'une famille)
*/
{
  if(L==NULL) return NULL; /* liste vide */
  if(L->type!=T_NB) return List2Graph(L,code); /* si graphe */

  /* ici on a donc une famille */
  int n=L->item; /* nb de graphes dans la liste */
  list *T;
  
  if(n<=0){ /* famille vide */
    if(code&16) /* libÃ¨re Ã©ventuellement L */
      while(L){
	T=L->next; /* ici L<>NULL */
	free(L);
	L=T;
      }
    return NULL;
  }

  int i,id;
  graph* F=new_graph(0);
  list* P;

  F->f=n;
  ALLOC(F->G,n); /* F->G[0..n[: tableau de n pointeurs de graphes */
  T=L; L=L->next;
  if(!(code&16)) free(T); /* on libÃ¨re l'Ã©lÃ©ment (n,T_NB) */

  /* ici L=dÃ©but du 1er graphe de la famille */
  for(i=0;i<n;i++){ /* pour chaque graphe */
    /* ici L pointe sur un Ã©lÃ©ment (id,T_ID) */
    if((L==NULL)||(L->type!=T_ID)) Erreur(13);
    id=L->item; /* identifiant du graph i */
    T=L->next;
    if((code&16)==0) free(L); /* on libÃ¨re l'Ã©lÃ©ment (id,T_ID) */
    P=L=T; /* P=L=T=tÃªte courante du graphe i */
    while((L)&&(L->type!=T_ID)){ P=L; L=L->next; } /* cherche la fin du graphe i */
    /* ici le graphe i va de T Ã  P */
    P->next=NULL; /* on coupe la liste */
    F->G[i]=List2Graph(T,code); /* T=liste du graphe i */
    F->G[i]->id=id; /* Attention! F->G[i] n'existe qu'aprÃ¨s l'appel Ã  List2Graph() */
    if(code&16) P->next=L; /* recolle la liste si on ne souhaite pas la libÃ©rer */
    }

  /* Ã©ventuellement trie la famille suivant les IDs */
  if(code&8) QSORT(F->G,F->f,fcmp_graphid);

  /* extrait le premier graphe */
  if(code&32){
    graph* G=GraphCopy(F->G[0]); /* copie le premier graphe */
    free_graph(F); /* libÃ¨re complÃ¨tement la famille F */
    F=G; /* F=premier graphe */
  }

  return F;
}


list *File2List(string const file)
/*
  Lit le fichier "file" contenant un graphe (ou une famille) au format
  standard, orientÃ©s ou non, et retourne le contenu dans une liste.
  Tient compte de -shift mais pas des noms originaux (-label 1). Dans
  le cas d'une famille de graphes, il est possible de spÃ©cifier un
  "range" pour "file" avec la forme: "file:range" oÃ¹ "range" est une
  liste de valeurs ayant la mÃªme signification que pour "-filter F id
  value". Par exemple, "file:5" spÃ©cifie le graphe d'identifiant 5, et
  "file:5-8" est la famille contenant les graphes d'identifiant
  5,6,7,8. Notez que "-:5" est le graphe d'identifiant 5 de la famille
  lue depuis l'entrÃ©e standard.

  Chaque Ã©lÃ©ment de la liste est une paire d'entiers (item,type) oÃ¹
  "type" prÃ©cise le rÃ´le jouÃ© par l'entier "item". Voir l'enum pour
  une description des types.

  Si, par exemple, le graphe est "0-1 2 1-2-3 5->4" alors la liste
  retournÃ©e sera { (0,T_NODE), (1,T_EDGE), (2,T_NODE), (1,T_NODE),
  (2,T_EDGE), (3,T_EDGE), (5,T_NODE), (4,T_ARC) }.

  Si le graphe est "0-(2 3) 4 5->(6 7)" alors la liste retournÃ©e sera
  { (0,T_NODE), (2,T_OPENE), (3,T_OPENE), (4,T_NODE), (5,T_NODE),
  (6,T_OPENA), (7,T_ARC) }. Autrement dit, T_OPENE ou T_OPENA dÃ©crive
  un groupe d'arÃªtes ou d'arcs.

  NB: "i-j-(k ...)" est correct mais pas "i-(j ...)-k".

  La fonction est gÃ©nÃ©ralisÃ©e Ã  la lecture d'une famille de graphes.
  Si le fichier contient "[5] 0-1 [8] 0->2-1" alors la liste
  contiendra { (2,T_NB), (5,T_ID), (0,T_NODE), (1,T_EDGE), (8,T_ID),
  (0,T_NODE), (2,T_ARC), (1,T_EDGE) }, oÃ¹ le premier Ã©lÃ©ment (n,T_NB)
  signifie qu'il s'agit d'une famille de n graphes, et oÃ¹ (u,T_ID)
  signifie que u est l'identifiant du graphe Ã  venir.
*/
{
  FILE *f;
  list *T; /* tÃªte de la liste */
  list *L; /* Ã©lÃ©ment courant */
  list *P; /* sauvegarde le dernier Ã©lÃ©ment (sentinelle qui faudra supprimer) */
  int read=1; /* pour InRange(), par dÃ©faut on lit tout */
  string r=NULL,s;
  char c[2];
  int range[CMDMAX]={2,R_TRUE}; /* par dÃ©faut: range toujours vrai */
  unsigned v; /* valeur lue */
  long p; /* position dans le fichier f */
  int n=0; /* nb de graphes dans la famille */
  int t=-1; /* t<0 si on est pas dans un groupe */

  T=P=L=new_list(); /* crÃ©e la liste */

  /* TO DO: si file commence par " ", alors file reprÃ©sente le graphe
     lui-mÃªme (pour future option -add/-del). Par exemple, file="
     5-6,7-8-0". Dans ce cas on Ã©crit un fichier temporaire avec "5-6
     7-8-0" et on continue normalement. On dÃ©truit ensuite ce
     fichier. NB: le contenu de file est modifiÃ©. */
  /*
    string s;
    if(file[0]==' '){
    for(v=0; v<strlen(file); v++) if(file[v]==',') file[v]=' ';
    s=strdup("/tmp");
    f=fopen(s,"rw");
    fputs(f,file);
    rewind(f); // pas la peine de fermer le fichier
    file=s;
  }
  */
  
  /* ouverture du fichier: file ou file:range */

  f=strcmp(file,"-")? fopen(file,"r"):stdin;
  if(f==NULL){ /* on a pas rÃ©ussit Ã  ouvrir file */
    fclose(f); /* il faut fermer le fichier, mÃªme si c'est NULL ! */
    if((r=strchr(file,':'))==NULL) Erreur(7); /* est-ce file:range ? */
    *r='\0'; /* coupe file en (prÃ©fixe,r=range) */
    f=strcmp(file,"-")? fopen(file,"r"):stdin;
    if(f==NULL) Erreur(7); /* on a pas rÃ©ussit Ã  ouvrir file */
    *r++ = ':'; /* rÃ©tablit le nom de fichier original de file */
    ReadRange(r,range); /* lecture du range */
  }

  /* lecture du fichier */

  while(!feof(f)){
    p=ftell(f); /* ftell() vaut toujours -1 si f=stdin (non seekable) */
    fscanf(f,"//%*[^\n]\n"); /* essaye de lire "//" */
    if(ftell(f)>=p+2) continue; /* on a lu au moins 2 caractÃ¨res -> commentaire */
    fseek(f,p,SEEK_SET);
    if(fscanf(f," [%u]",&v)>0){
      read=InRange(v,range);
      if(read){ L=Insert(P=L,v,T_ID); n++; }
      continue;
    }
    if(read){
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-%u",&v)>0){ L=Insert(P=L,v,T_EDGE); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">%u",&v)>0){ L=Insert(P=L,v,T_ARC); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-( %u",&v)>0){ L=Insert(P=L,v,T_OPENE);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_EDGE; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">( %u",&v)>0){ L=Insert(P=L,v,T_OPENA);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_ARC; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f," %u",&v)>0){
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0){ L=Insert(P=L,v,t); t=-1; }
	else{ L=Insert(P=L,v,(t<0)?T_NODE:t); fseek(f,p,SEEK_SET); }
	continue;
      }
      
      /* ici on a rien trouvÃ©: est-ce une erreur de format ? */
      fseek(f,p,SEEK_SET);
      s=fgets(c,2,f); /* lit au plus un caractÃ¨re */
      if((s==NULL)||(c[0]==' ')||(c[0]=='\n')) continue; /* ok si ' ' ou '\n' */
      Erreur(28); /* mauvais format sinon */
    }

    /* ici on a rien trouvÃ©, mais read est faux */
    fseek(f,p,SEEK_SET); /* on a rien trouvÃ© */
    fscanf(f," %*c"); /* lit au moins un caractÃ¨re */
  }
  
  fclose(f); /* on ferme le fichier */
  free(L); /* supprime le dernier Ã©lÃ©ment (la sentinelle) */
  if(L==T) return NULL; /* si on a lu aucun Ã©lÃ©ment */
  P->next=NULL;

  if(n>0){ /* il s'agit d'une famille */
    /* on ajoute un nouvel Ã©lÃ©ment en tÃªte de la liste indiquant le
       nombre de graphes de la famille */
    L=new_list();
    L->item=n; /* nombre de graphes de la famille */
    L->type=T_NB;
    L->next=T;
    T=L; /* nouvelle tÃªte de liste */
  }

  return T; /* on retourne la tÃªte */
}


graph* File2Graph(string const file,int const code)
/*
  Renvoie un graphe (ou une famille) Ã  partir d'un fichier. Pour
  "code" voir List2Graph() & List2Family(). La liste intermÃ©diaire
  calculÃ©e par File2List() est toujours libÃ©rÃ©e.
*/
{
  graph* G=List2Family(File2List(file),code&(63-16)); /* annule le bit-4 */
  if(G==NULL) Erreur(15);
  return G;
}


/***********************************

       ROUTINES EN VRAC

***********************************/


double PosAspect(const query* const Q)
/*
  Donne le coefficient par lequel les positions Q->xpos/Q->ypos vont
  Ãªtre multipliÃ©es pour le format dot pour permettre une taille de
  sommets raisonable par rapport Ã  la longueur des arÃªtes. On tient
  compte de Q->n et de BOXX et BOXY.
*/
{
  double w=C32*sqrt(Q->n); /* la largeur est en sqrt(n) */
  if((BOXX>0)&&(BOXY>0)) w /= min(BOXX,BOXY);
  if(LABEL>0) w *= 3; /* augmente l'aspect si besoin des LABELs (et POS) */
  return w;
}


void BoundingBox(const query* const Q)
/*
  Calcule XMIN,YMIN,XMAX,YMAX des tableaux Q->xpos/Q->ypos. Il faut
  que Q->n > 0 et Q->xpos,Q->ypos <> NULL.
*/
{
  XMIN=XMAX=Q->xpos[0];
  YMIN=YMAX=Q->ypos[0];
  for(int i=1;i<Q->n;i++){
    XMIN=min(XMIN,Q->xpos[i]);
    YMIN=min(YMIN,Q->ypos[i]);
    XMAX=max(XMAX,Q->xpos[i]);
    YMAX=max(YMAX,Q->ypos[i]);
  }
}


static inline double angle(double const x,double const y)
/*
  Renvoie l'angle de [0,2ðœ‹[ en radian du point de coordonnÃ©es
  cartÃ©siennes (x,y) par rapport Ã  l'axe des abscisses (1,0). NB:
  atan(y/x) donne un angle [-ðœ‹/2,+ðœ‹/2] ce qui n'est pas ce que l'on
  veut. On renvoie 0 si (x,y)=(0,0).
*/
{
  if(x==0){
    if(y>0) return M_PI_2;
    if(y<0) return M_PI_2+M_PI;
    return 0;
  }

  // atan(y/x) renvoie un angle entre -ðœ‹/2 et +ðœ‹/2
  // l'angle est correct si x>0 et y>0
  // si x,y de signe contraire alors atan(y/x)=-atan(y/x)
  
  double const a=atan(y/x);

  if(x>0){
    if(y>0) return a;
    return a+M_2PI;
  }
  
  return a+M_PI;
}


int fcmp_angle(const void *P,const void *Q)
/* Compare deux angles <XP,YP> et <XQ,YQ> pour qsort2(). */
{
  double const p=angle(*(double*)P,*(((double*)P)+1));
  double const q=angle(*(double*)Q,*(((double*)Q)+1));
  return (p>q) - (p<q);
}


static inline double det(double const a1,double const a2,
			 double const b1,double const b2)
/*
  Renvoie le dÃ©terminant de deux vecteurs colonnes A=(a1,a2) et
  B=(b1,b2). Il vaut 0 si les vecteurs A et B sont colinÃ©aires. Le
  signe de det(A,B) est celui de sin(A,B), le sinus de l'angle entre
  (OA) et (OB), O=(0,0) Ã©tant l'origine. Dit autrement, si det(A,B)>0
  alors le point B est au-dessus de la droite (OA), si det(A,B)<0 il
  est en dessous, et si det(A,B)=0, alors il est sur la droite (OA).
*/
{
  return (a1*b2)-(a2*b1);
}


/*
  Rappels de gÃ©omÃ©trie 2D:
  (voir aussi distgone() et rlt())

  1. Equation cartÃ©sienne d'une droite
     
     D = { (x,y) : a*x + b*y = c }
         avec (a,b) != (0,0)
	 mais a=0 ou b=0 possible

  2. Equation cartÃ©sienne d'une droite D passant par deux points
     A=(xa,ya) et B=(xb,yb)

     D: (yb-ya)*x - (xb-xa)*y = xa*yb - ya*xb = det(A,B)

  3. Intersection de deux droites cartÃ©siennes
     D1: a1*x - b1*y = c1
     D2: a2*x - b2*y = c2

     Soient les vecteurs colonnes: A=(a1,a2), B=(b1,b2), C=(c1,c2)
     Si det(A,B)=det(A,C)=0, alors les droites sont confondues
     Si det(A,B)=0 et det(A,C)<>0, alors il n'y a pas d'intersection
     Si det(A,B)<>0, alors il y a une seule intersection (x0,y0):
     x0 = -det(B,C)/det(A,B) et y0 = -det(A,C)/det(A,B)

     Ex1: 2x - (-4)y = 20   det(A,B)=-16+28  =  12
          7x - (-8)y = 52   det(B,C)=-208+160= -48
	                    det(A,C)=104-140 = -36
	  x0=-det(B,C)/det(A,B)=48/12=4
	  y0=-det(A,C)/det(A,B)=36/12=3

     Ex2: 7x - (-5)y = 11   det(A,B)=-21+5 = -16
          1x - (-3)y =  5   det(B,C)=-25+33=   8
	                    det(A,C)=35-11 =  24
	  x0=-det(B,C)/det(A,B)=-8/-16=1/2
	  y0=-det(A,C)/det(A,B)=-24/-16=3/2

  4. det(X,-Y) = -det(X,Y)
     det(Y,X) = -det(X,Y)
*/


/* comme QSORT(), mais pour qsort2() */
#define QSORT2(T1,T2,n,f) qsort2(T1,sizeof(*(T1)),T2,sizeof(*(T2)),n,f)


void qsort2(void* const T1,int const w1,
	    void* const T2,int const w2,
	    int const n,int (*fcmp)(const void*,const void*))
/*
  Trie simultanÃ©ment deux tableaux T1 et T2 non NULL chacun de n
  Ã©lÃ©ments selon la fonction de comparaison fcmp(). Ici w1 (resp. w2)
  reprÃ©sentent la taille d'un Ã©lÃ©ment de T1 (resp. T2). Pour cela on
  construit un nouveau tableau T oÃ¹ chaque Ã©lÃ©ment T[i] est composÃ©e
  de la paire d'Ã©lÃ©ments <T1[i],T2[i]> oÃ¹ T1[i] est stockÃ© juste avant
  T2[i]. Puis on applique qsort(T,n,w1+w2,fcmp), et enfin on remet
  dans T1 et T2 les Ã©lÃ©ments ainsi triÃ©s. La fonction de comparaison
  fcmp() doit pouvoir s'appliquer Ã  une paire <T1[i],T2[i]>, mÃªme si
  elle peut trÃ¨s bien s'appliquer seulement sur T1.

  Mais le plus souvent fcmp() peut s'appliquer seulement Ã  T1, comme
  dans l'exemple suivant:

  int A[]={ 5, 8, 2, 9, 1, 6, 3, 7, 4};
  int B[]={-5,-8,-2,-9,-1,-6,-3,-7,-4};

  QSORT2(A,B,9,fcmp_int);
  -> A[]={ 1, 2, 3, 4, 5, 6, 7, 8, 9}
  -> B[]={-1,-2,-3,-4,-5,-6,-7,-8,-9}

  QSORT2(B,A,9,fcmp_int);
  -> A[]={ 9, 8, 7, 6, 5, 4, 3, 2, 1}
  -> B[]={-9,-8,-7,-6,-5,-4,-3,-2,-1}

*/
{
  if(n<2) return; /* ne rien faire si pas au moins 2 Ã©lÃ©ments */

  int const w=w1+w2; /* taille en octets d'une paire <T1[i],T2[i]> */
  void *T=malloc(n*w); /* taille en octets de T */
  if(T==NULL) Erreur(3); /* si problÃ¨me mÃ©moire */
  void *t,*t1,*t2;
  int i;

  /* fusionne T1 et T2 dans T */
  t=T; t1=T1; t2=T2;
  for(i=0;i<n;i++,t+=w){
    memcpy(t,   t1,w1); t1 += w1;
    memcpy(t+w1,t2,w2); t2 += w2;
  }

  /* trie T */
  qsort(T,n,w,fcmp);
  
  /* sÃ©pare T en T1 et T2 */
  t=T; t1=T1; t2=T2;
  for(i=0;i<n;i++,t+=w){
    memcpy(t1,t,   w1); t1 += w1;
    memcpy(t2,t+w1,w2); t2 += w2;
  }

  free(T);
  return;
}


int InitXY(query* const Q)
/*
  Initialise les tableaux Q->xpos et Q->ypos suivants les options
  Ã©ventuelles de -xy. Les variables suivantes, en plus de Q->xpos et
  Q->ypos, peuvent Ãªtre mise Ã  jour: Q->n, XMIN,XMAX,YMIN,YMAX
  (BoundingBox), VSIZESTD, VSIZEXY. La fonction renvoie 0 pour
  indiquer que tout c'est bien passÃ©.

  Il y a deux Ã©tapes:
   1. gÃ©nÃ©ration des points Q->xpos,Q->ypos selon XYtype
   2. application des options: XYunique, ...
*/
{
  int i;

  /*************************************
    1. Initialise Q->xpos,Q->ypos
  *************************************/

  for(;;){ /* pour pouvoir faire un break; */
    
    if(XYtype==XY_USER) break; /* coordonnÃ©es dÃ©finies par l'utilisateur */
    if(XYtype==XY_FILE){ /* charge Ã  partir d'un fichier et met Ã  jour Q->n */
      LoadXY(Q,FILEXY);
      break;
    }
    if(XYratio<=0) Erreur(6); // paramÃ¨tre incorrect
    if(Q->n<0) Q->n=0;
    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);
    /* ici Q->xpos,Q->ypos existent */

    switch(XYtype){ /* type de gÃ©nÃ©ration des points */


    case XY_UNIF:
      /* uniforme dans [0,1[ x [0,XYratio[ */

      for(i=0;i<Q->n;i++){
	Q->xpos[i]=RAND01;
	Q->ypos[i]=XYratio*RAND01;
      }
      break;

      
    case XY_PLAW:{
      /* loi puissance autour des graines choisies dans [0,1[ */

      if(XYseedk<=0) Erreur(6); // paramÃ¨tre incorrect
      ALLOC(XSEED,XYseedk);
      ALLOC(YSEED,XYseedk);
      double sx=0,sy=0; /* calcule (sx,sy), le barycentre des XYseedk graines */
      for(i=0;i<XYseedk;i++){
	XSEED[i]=RAND01;         sx+=XSEED[i];
	YSEED[i]=RAND01*XYratio; sy+=YSEED[i];
      }
      sx /= XYseedk; sy /= XYseedk;
      sx -= 0.5; sy -= XYratio/2;
      for(i=0;i<XYseedk;i++){ /* centre par rapport au barycentre */
	XSEED[i] -= sx; /* enlÃ¨ve le barycentre puis dÃ©cale de 0.5 */
	YSEED[i] -= sy;
      }
      /* on gÃ©nÃ¨re les points autour des graines */
      double const r=sqrt(log(XYseedk+1)/XYseedk); /* rayon r */
      /* le +1 est important car abÃ©rant pour XYseedk=1 */
      int k;
      for(i=0;i<Q->n;i++){
	k=randomu(XYseedk);    /* on choisit la graine numÃ©ro k au hasard */
	sx=M_2PI*RAND01;  /* angle alÃ©atoire */
	sy=r*pow(RAND01,XYpower); /* longueur alÃ©atoire */
	Q->xpos[i]=XSEED[k]+sy*cos(sx);
	Q->ypos[i]=YSEED[k]+sy*sin(sx)*XYratio;
      }
      break;
    }
      
    case XY_PERM:{
      /* permutation de [0,Q->n[ */
      
      NALLOC(int,P,Q->n);
      for(i=0;i<Q->n;i++) Q->xpos[i]=(double)(P[i]=i); // initialise aussi P[i]
      Permute(P,Q->n); // modifie P[i]
      for(i=0;i<Q->n;i++) Q->ypos[i]=(double)P[i];
      free(P);
      break;
    }


    case XY_MESH:
      /* grille de paramÃ¨tre Xmesh x Ymesh */

      if((Xmesh<=0)||(Ymesh<=0)) Erreur(6);
      for(i=0;i<Q->n;i++) Q->xpos[i]=(double)(i%Xmesh),Q->ypos[i]=(double)(i/Xmesh);
      break;

      
    case XY_CYCLE:{
      /* cycle de rayon 1 et de centre (0,0) */

      double t=0;
      double const a=M_2PI/Q->n;
      for(i=0;i<Q->n;i++){
	Q->xpos[i]=cos(t);
	Q->ypos[i]=sin(t)*XYratio;
	t += a;
      }
      break;
    }


    case XY_CIRCLE:
    case XY_DISK:
    case XY_HYPER:{
      /* points sur un cercle, dans un disque (star-shaped polygon),
	 ou un disque hyperbolique de centre (0,0) et rayon <= 1 */
      
      /*
	Attention ! pour gÃ©nÃ©rer des points alÃ©atoires uniformes sur
        un disque unitÃ©, il faut faire: a=M_2PI*RAND01,
        r=sqrt(RAND01), puis (x,y)=(r*cos(a),r*sin(a)). Si on utilise
        seulement r=RAND01 (sans le sqrt) alors les points se
        retrouvent plus concentrÃ©s au centre du disque:
        http://mathworld.wolfram.com/DiskPointPicking.html
      */
      NALLOCZ(double,A,Q->n,M_2PI*RAND01); // Q->n angles alÃ©atoires de [0,2ðœ‹[
      QSORT(A,Q->n,fcmp_double); // trie les angles (ne tient pas compte de ROUND ...)
      double r;
      for(i=0;i<Q->n;i++){ // transforme coordonnÃ©es polaires en cartÃ©siennes
	switch(XYtype){
	case XY_CIRCLE: r=1; break;
	case XY_DISK:   r=sqrt(RAND01); break;
	case XY_HYPER:  r=exp(-RAND01*XYpower); break;
	default:        r=RAND01;
	}
	Q->xpos[i]=r*cos(A[i]);
	Q->ypos[i]=r*sin(A[i])*XYratio;
      }
      free(A);
      break;
    }

    case XY_RPOLY:{
      /* points dans un polygone convexe rÃ©gulier */
      
      if(XYpoly<3) Erreur(42);
      /*
	Principe:

	L'Algorithme s'applique indÃ©pendemment Ã  chacun des Q->n
	points. On considÃ¨re le triangle dÃ©fini par un seul cotÃ© du
	polygone rÃ©gulier, d'angle Ï´=2ðœ‹/p oÃ¹ p=XYpoly. On l'oriente
	pour que l'axe des abscisses corresponde Ã  la mÃ©diane de
	l'angle Ï´. Le cotÃ© du polygone est ainsi un segment vertical
	d'abscisse cos(Ï´/2) et de hauteur 2|sin(Ï´/2)|. On peut alors
	tirer un point M alÃ©atoirement uniforme dans ce triangle, puis
	on tourne le point M d'un angle i*Ï´ avec i alÃ©atoire dans
	[0,p[.
	
	Pour tirer alÃ©atoirement un point M dans un triangle (O,V1,V2)
	dont un coin (ici O) est l'origine on peut faire comme suggÃ©rÃ©
	dans http://mathworld.wolfram.com/TrianglePointPicking.html.
	On choisit r1,r2 uniformes dans [0,1], puis on construit le
	point M=r1*V1+r2*V2. Ce point est alÃ©atoire dans le
	parallÃ©logramme dÃ©finit par les 4 points (O,V1,V2,V1+V2).  Si
	M n'est pas dans le triangle, alors soit on recommence, soit
	on prend le symÃ©trique tombant dans le triangle. Le nombre
	d'essais moyen est 2.

	Ici V1=(cos(Ï´/2),sin(Ï´/2)) et V2=(cos(Ï´/2),-sin(Ï´/2)). M est
	dans le triangle si x(M)<=X(V1)=X(V2)=cos(Ï´/2). Sinon, on
	change M en son symÃ©trique par rapport au point (cos(Ï´/2),0).
      */
      double const t=M_2PI/XYpoly; // t=Ï´=angle dÃ©fini par un cotÃ© du polygone
      double const c=cos(t/2); // c=abscisse du cotÃ© vertical
      double const s=sin(t/2); // s=demi-hauteur du cotÃ© vertical
      double a,r,sx,sy;

      for(i=0;i<Q->n;i++){ // indÃ©pendemment pour chaque point
	a=RAND01,r=RAND01; // tire un point M uniformÃ©ment dans le triangle
	sx=c*(a+r),sy=s*(a-r); // M=(sx,sy) dans le parallÃ©logramme
	if(sx>c) sx=2*c-sx,sy=-sy; // prend le symÃ©trique de M par rapport Ã  (c,0)
	r=hypot(sx,sy); // (a,r)=coordonnÃ©es polaires de la rotation de M
	a=angle(sx,sy)+t*randomu(XYpoly); // a=angle avec rotation alÃ©atoire
	Q->xpos[i]=r*cos(a);
	Q->ypos[i]=r*sin(a)*XYratio;
      }
      break;
    }

      
    case XY_CONVEX:{
      /* points en position convexe, algorithme en n^2 */
      
      double sx,sy,tx,ty,mx,my,d,x,y,xa,ya,xb,yb,a1,b1,a2,b2,c2,h;
      int t,j,k,b=RANDbit; // b=0 ou 1

      // on suppose que les points P(0) ... P(i-1) sont dÃ©jÃ  en
      // position convexes, l'intÃ©rieur contenant l'origine (0,0), on
      // souhaite placer un nouveau point Q=P(i)

      for(i=0;i<Q->n;i++){
	
	// calcule tx=angle alÃ©atoire du nouveau point Q
	tx=M_2PI*RAND01;

	// pour i=0,1,2, tx est rÃ©duit Ã  l'un des trois secteurs non
	// adjacents d'angle ðœ‹/3 (les 3 cÃ´nes positifs ou nÃ©gatifs des
	// 6 secteurs d'angle ðœ‹/3) de sorte que l'origine sera
	// forcÃ©ment Ã  l'intÃ©rieur du triangle (b=bit alÃ©atoire). On
	// pourrait penser Ã  recentrer les points selon leur
	// barycentre, mais cela produit des erreurs car les points
	// peuvent maintenant Ãªtre en dehors du cercle de rayon 1.
	if(i<3) tx=(2*i+b)*M_2PI/6 + fmod(tx,M_2PI/6);
	
	DEBUG(printf("tx=%.02lf %03.0lf\n",tx,360*tx/M_2PI););
	
	// cherche k dans [0,i] tq P(k-1) < Q < P(k). Si Q est le plus
	// grand des points, alors k=i. Si c'est le plus petit, k=0.
	// On fait une recherche linÃ©aire, bien qu'on pourrait la
	// faire en log(i)
	
	for(k=0;k<i;k++){ // P(k)=(Q->xpos[k],Q->ypos[k])
	  ty=angle(Q->xpos[k],Q->ypos[k]);
	  if(ty==tx) tx=nextafter(tx,7); // Ã©vite les angles Ã©gaux. NB: tx < 2ðœ‹ < 7
	  if(tx<ty) break;
	}
	DEBUG(PRINT(k););
	
	// S=[sx,sy[ segment d'angle tx oÃ¹ l'on cherche Q
	// D1=droite contenant Q passant par (0,0) contenant S
	// D1: a1*x - b1*y = 0

	a1=sin(tx);
	b1=cos(tx);
	
	// M=(mx,my)=point maximum pour Q selon la direction D1. De
	// maniÃ¨re gÃ©nÃ©rale, M est sur l'ellipse de demi-hauteur
	// XYratio, la demi-longueur valant 1. L'Ã©quation paramÃ©trique
	// donne le rayon r(t) en fonction de l'angle t qui est: r(t)
	// = b/sqrt(1-e*cos(t)^2) oÃ¹ a=demi-longueur=1, b=demi-hauteur
	// et e=1-(b/a)^2. Si XYratio=1, alors r(t)=1.

	d=XYratio/sqrt(1-(1-XYratio*XYratio)*b1*b1); // d=r(tx)
	DEBUG(PRINT(d););
	mx=d*b1;
	my=d*a1;

	// S = [sx,sy[
	sx=nextafter(0,1); // borne inf de S
	sy=hypot(mx,my);   // borne sup de S

	// on rÃ©fuit S en fonction des trois droites dÃ©finies par les
	// 4 points successifs P(k-2), P(k-1), P(k), P(k+1)
	if(i>2){ // rien Ã  faire si i=0,1 ou 2

	  /* Il est possible d'avoir P(k-2)=P(k+1) lorsque i=3. Pour
	     les trois droites successives [P(k-2),P(k-1)],
	     [P(k-1),P(k)] et enfin [P(k),P(k+1)], on va calculer leur
	     droite D2, puis l'intersection entre D1 et D2. On calcule
	     ensuite une nouvelle borne pour sx,sy. */

	  j=(k-2+i)%i; // au dÃ©part A=P(k-2) et B=P(k-1).
	  for(t=-1;t<=1;t++){ // on rÃ©pÃ¨te trois fois: t=-1,0,+1
	    xa=Q->xpos[j],ya=Q->ypos[j]; // point A=(xa,ya)
	    j=(j+1)%i; // point suivant
	    xb=Q->xpos[j],yb=Q->ypos[j]; // point B=(xb,yb)
	    // droite D2 passant par A et B
	    //  D2: a2*x - b2*y = c2
	    a2=yb-ya;
	    b2=xb-xa;
	    c2=det(xa,ya,xb,yb);
	    // (x,y)=intersection de D1 et D2
	    //  D1: a1*x - b1*y = 0
	    //  D2: a2*x - b2*y = c2
	    d=det(a1,a2,b1,b2);
	    if(d==0) continue; // droite suivante si pas d'intersection
	    x=-det(b1,b2,0,c2)/d;
	    y=-det(a1,a2,0,c2)/d;
	    // D1 est en fait une demi-droite. Est-ce que (x,y)
	    // appartient Ã  cette demi-droite ? c'est-Ã -dire est-ce
	    // que les vecteurs (b1,a1) et (x,y) sont dans le mÃªme
	    // sens ou opposÃ© ?  mÃªme sens <=> (b1*x>0)&&(a1*y>0)
	    if((b1*x>0)&&(a1*y>0)){ // si mÃªme sens, on modifie sx ou sy
	      DEBUG(printf("sx=%lf sy=%lf -> ",sx,sy););
	      h=hypot(x,y);
	      if(t) sy=min(sy,h); // t=-1 ou +1
	      else sx=max(sx,h);  // t=0
	      DEBUG(printf("sx=%lf sy=%lf h=%lf t=%i\n",sx,sy,h,t););
	    }
	  }
	  // ici sx<sy sinon les points ne sont pas en position convexe
	  DEBUG(if(sx>sy) printf("problÃ¨me: sx>sy\n"););
	}
	
	// insÃ¨re le point Q entre P(k-1) et P(k) en dÃ©calant
	// [P(k)...P(i-1)] vers [P(k+1)...P(i)]
	for(t=i;t>k;t--) Q->xpos[t]=Q->xpos[t-1],Q->ypos[t]=Q->ypos[t-1];

	// calcule, en fonction de tx et S=[sx,sy[, le nouveau point Q
	// qui devient P(k). On met sqrt() pour Ãªtre plus proche d'une
	// distribution uniforme sur un disque (donc plus souvent
	// proche du bord du cercle que du centre)
	sx += (sy-sx)*sqrt(RAND01); // sx=point alÃ©atoire dans S=[sx,sy[
	if(sx==0) sx=nextafter(0,1); // on force sx!=0 pour avoir P(k)â‰ (0,0)
	Q->xpos[k]=sx*cos(tx);
	Q->ypos[k]=sx*sin(tx);

	DEBUG(
	      PRINT(i);
	      for(t=0;t<=i;t++)
		printf("point %i: %+.02lf %+.02lf \t%+.02lf %03.lf\n",
		       t,Q->xpos[t],Q->ypos[t],angle(Q->xpos[t],Q->ypos[t]),
		       360*angle(Q->xpos[t],Q->ypos[t])/M_2PI);printf("\n");
	      );
      }
      break;
    }

      
    case XY_CONVEX2:{
      /* points en position convexe v2, algorithme en n*log(n) */
      
      /*
	Principe:

	On part de points alÃ©atoires dans [0,1[Â², puis on calcule
        (pour les points finaux) la diffÃ©rence entre deux points
        consÃ©cutifs. La somme des n diffÃ©rences est nulle. On trie ces
        points selon l'angle, puis on dessine de proche en proche les
        points de l'enveloppe convexe (avec chaque fois un angle
        croissant donc).
      */

      for(i=0;i<Q->n;i++) Q->xpos[i]=RAND01,Q->ypos[i]=RAND01; /* points alÃ©atoires */
      if(Q->n>0){ /* NB: il faut Q->n>0 */
	double const x0=Q->xpos[0], y0=Q->ypos[0]; // sauvegarde le 1er point
	for(i=0;i<Q->n-1;i++) Q->xpos[i]-=Q->xpos[i+1],Q->ypos[i]-=Q->ypos[i+1]; // diffÃ©rences
	Q->xpos[i] -= x0,Q->ypos[i] -= y0;
	QSORT2(Q->xpos,Q->ypos,Q->n,fcmp_angle); // trie les angles
	for(i=1;i<Q->n;i++) Q->xpos[i] += Q->xpos[i-1],Q->ypos[i] += Q->ypos[i-1]; // dessin
	}
      break;
    }
      
    }/* fin du switch(XYtype) */
    break; // pour ne pas boucler
  } /* fin du for(;;) */
  
  
  /*************************************
    2. Application des options
  *************************************/
  
  if((Q->xpos==NULL)||(Q->ypos==NULL)) Erreur(8);
  double sx,sy;

  if(XYnoiser>0) /* "noise" doit Ãªtre avant "box" */
    for(i=0;i<Q->n;i++){
      sx=M_2PI*RAND01; /* angle alÃ©atoire */
      sy=XYnoiser*pow(RAND01,XYnoisep); /* longueur alÃ©atoire */
      Q->xpos[i] += sy*cos(sx); /* dÃ©cale Q->xpos */
      Q->ypos[i] += sy*sin(sx); /* dÃ©cale Q->ypos */
    }
  
  if((BOXX>0)&&(BOXY>0)){ /* "box" doit Ãªtre aprÃ¨s "noise" */
    BoundingBox(Q); /* calcule les BB */
    if(Q->n<2) sx=sy=0;
    else{
      sx=BOXX/(XMAX-XMIN);
      sy=BOXY/(YMAX-YMIN);
    }
    for(i=0;i<Q->n;i++){
      Q->xpos[i]=sx*(Q->xpos[i]-XMIN);
      Q->ypos[i]=sy*(Q->ypos[i]-YMIN);
    } // ici les BB sont obsolÃ¨tes
  }
  
  if(ROUND<DBL_DIG){ /* arrondit Ã©ventuellement les coordonnÃ©es */
    sx=pow(10,ROUND);
    for(i=0;i<Q->n;i++){
      Q->xpos[i]=rint(Q->xpos[i]*sx)/sx;
      Q->ypos[i]=rint(Q->ypos[i]*sx)/sx;
    }
  }

  if(XYunique){ /* Ã©limine les doubles, en triant les points */
    NALLOC(point,P,Q->n);
    for(i=0;i<Q->n;i++) P[i].x=Q->xpos[i],P[i].y=Q->ypos[i];
    QSORT(P,Q->n,fcmp_point); /* tri les points */
    point p=P[0]; p.x -= 1.0; /* ici le point p <> du 1er Ã©lÃ©ment */
    int k;
    for(i=k=0;i<Q->n;i++,k++) /* k=nombre de points uniques */
      if(fcmp_point(P+i,&p)){ /* copie que si diffÃ©rent de l'Ã©lÃ©ment p */
	p=P[i];
	Q->xpos[k]=p.x;
	Q->ypos[k]=p.y;
      }
    free(P);
    if(k<Q->n){
      Q->n=k;
      REALLOC(Q->xpos,Q->n);
      REALLOC(Q->ypos,Q->n);
    }
  }

  /* on (re)calcule les BB */
  BoundingBox(Q);

  /* mise Ã  jour de la taille des sommets */
  VSIZESTD *= XYvsize;
  VSIZEXY  *= XYvsize;

  return 0; // sortie normale
}


color* GradColor(const color* const T,int const n,int const m)
/*
  Retourne un tableau de m couleurs formant un dÃ©gradÃ© obtenu Ã  partir
  d'un tableau T de n couleurs. Pour avoir un dÃ©gradÃ© simple d'une
  couleur allant de T[0] Ã  T[1] il faut initialiser T[0],T[1] et poser
  n=2. Pour avoir un dÃ©gradÃ© cyclique, il suffit de rÃ©pÃ©ter la couleur
  T[0] en derniÃ¨re position de T (et ajouter 1 Ã  n, donc d'avoir
  T[n-1]=T[0]). Il faut dans tous les cas n>1 et m>0. On peut avoir
  m<n. Dans ce cas on prend la premiÃ¨re couleur de T, puis la i-Ã¨me
  couleur est (i*(n-1))/(m-1).
*/
{
  color c1,c2;
  int i,j,k,r,q,n1,m1,dr,dg,db;

  if(T==NULL) return NULL; /* normalement ne sert Ã  rien */

  NALLOC(color,P,m);
  c2=P[0]=T[0];
  if(m==1) return P;
  /* maintenant m >= 2 */

  m1=m-1; n1=n-1; /* valeurs utilisÃ©es souvent */

  if(m<=n){ /* cas oÃ¹ il y a moins de couleurs demandÃ©es que dans T */
    for(i=1;i<m;i++) /* m-1 fois */
      P[i]=T[(i*n1+m1-1)/m1]; /* le "+m-2" est pour arrondir Ã  l'entier sup. */
    return P;
  }

  /*
    Cas m>n.  Soient B_1,B_2,...B_(n-1) les n-1>0 blocs de couleurs,
    B_i commenÃ§ant juste aprÃ¨s la couleurs T[i-1] et se terminant avec
    la couleur T[i]. On doit rÃ©partir m-1 couleurs dans ces n-1 blocs,
    la couleurs T[0] Ã©tant dÃ©jÃ  dans P. On met alors
    floor((m-1)/(n-1)) couleurs par blocs, et une de plus pour B_i si
    i<=(m-1)%(n-1).
   */
  r=m1%n1; /* il reste r couleurs */
  q=(m1/n1)+(r>0); /* nombre de couleurs par blocs (+1 au dÃ©part si r>0) */
  for(i=j=k=1;i<n;i++){ /* on traite le bloc B_i, P[k]=prochaine couleur libre */
    c1=c2;   /* c1=T[i-1] */
    c2=T[i]; /* c2=T[i] */
    dr=c2.r-c1.r;
    dg=c2.g-c1.g;
    db=c2.b-c1.b;
    for(j=1;j<=q;j++){ /* extrapolation linÃ©aire de q points dans ]c1,c2] */
      P[k].r=c1.r+(j*dr)/q;
      P[k].g=c1.g+(j*dg)/q;
      P[k].b=c1.b+(j*db)/q;
      k++;
    }
    if(i==r) q--; /* une couleur de moins pour B_{r+1}...B_{n-1} */
  }
  return P;
}


int graphical(const int* const S,int k)
/*
  VÃ©rifie si la suite S = (n_0,d_0, n_1,d_1, ..., n_{k-1},d_{k-1})
  constituÃ©e de k couples est graphique ou pas, c'est-Ã -dire s'il
  existe au moins un graphe simple ayant exactement n_i sommets de
  degrÃ© d_i. On renvoie une valeur <0 si la sÃ©quence n'est pas
  graphique, et sinon on renvoie n = âˆ‘_i n_i, c'est-Ã -dire le nombre
  de sommets du graphe. Les n_i et d_i de S sont quelconques. En
  particuliers les d_i ne sont pas forcÃ©ment triÃ©s et diffÃ©rents.
  L'algorithme utilisÃ© est en O(SORT(n,k))=O(k*log(k)).

  L'algorithme est basÃ© sur le test d'ErdÅ‘s and Gallai (1960) (voir
  aussi [Skiena08, p. 464]) qui ont prouvÃ© qu'une suite d'entiers
  positifs (t_1,...,t_n) triÃ©e dans l'ordre dÃ©croissant est graphique
  ssi la somme est paire et la suite vÃ©rifie la propriÃ©tÃ© suivante
  pour chaque r=1..n-1: 

  (1) âˆ‘_{i=1}^r t_i <= r*(r-1) + âˆ‘_{i=r+1}^n min{r,t_i}

  Attention! ici t_i est le degrÃ© du sommet i-1. Cette propriÃ©tÃ© se
  gÃ©nÃ©ralise aux graphes orientÃ©s. Tripathi et Vijay (2003) on montrÃ©
  qu'il suffit de vÃ©rifier (1) pour les valeurs de r=1..n-1 telles que
  t_r > t_{r+1}, c'est-Ã -dire lorsque les t_i changent (voir
  [DF05]). On pourrait aussi se contenter de vÃ©rifier (1) pour les
  valeurs de r=1..R avec R = max{ r<n : t_r >= r }, mais on ne va pas
  vraiment utiliser R. (NB: R pourrait ne pas Ãªtre dÃ©fini si
  t_1=0. Mais dans ce cas la suite est graphique, puisque t_i=0 pour
  tous les i. On ne se servira pas de R.)

  Par exemple,

  i | 0 1   |t1|t2 t3 t4|      i | 0 1 2    |t1 t2|t3 t4 t5 t6|t7 t8|
  --+----   -------------   ou --+-------   -------------------------
  n | 1 3   |3 |1  1  1 |      n | 2 4 2    |4  4 |2  2  2  2 |1  1 |
  d | 3 1   r=1, car t1>t2     d | 4 2 1    r=2 et 6, car t2>t3 et t6>t7

  Dans un premier temps on calcule Ã  partir de S les tableaux n[i] et
  d[i] de sorte que les d[i] soient strictement dÃ©croissant. Donc
  d[i]>d[i+1] et n[i] est le nombre total de sommets de degrÃ© d[i]. Il
  faut faire un tri puis regrouper les termes Ã©gaux. On met Ã  jour k
  de sorte qu'il corresponde au nombre de degrÃ©s diffÃ©rents. Il suffit
  alors de tester (1) pour r = n[0], n[0]+n[1], ..., n[0]+...+n_[k-2],
  car t_{n[0]+...+n[j]} > t_{n[0]+...+n_[j+1]} pour tout j=0..k-2.

  Si on note r_j = n[0] +...+ n[j], il faut donc vÃ©rifier (1)
  seulement pour les valeurs de r suivantes: r_0, r_1,..., r_{k-2},
  soit les valeurs r_j pour j=0..k-2. On remarque que les n[j] degrÃ©s
  t_{r_j},...,t_{r_{j+1}-1} sont tous Ã©gaux Ã  d[j]. Donc, en dÃ©coupant
  les deux sommes de l'Ã©quation (1) en blocks de n[j] sommets
  consÃ©cutifs de mÃªme degrÃ©, on obtient:

  S1 = âˆ‘_{i=1}^r_j t_i = âˆ‘_{i=0}^j n[j]*d[j]   et
  S2 = âˆ‘_{i=r_j+1}^n min{r_j,t_i} = âˆ‘_{i=j+1}^k n[i]*min{r_j,d[i]}.

  L'Ã©quation (1) se rÃ©Ã©crit donc, pour chaque j=0..k-2: S1 <= P + S2
  oÃ¹ P = r_j*(r_j-1), ou encore:
  
  (2) âˆ‘_{i=0}^j n[j]*d[j] <= r_j*(r_j-1) + âˆ‘_{i=j+1}^k n[i]*min{r_j,d[i]}

  Soit m = âˆ‘_i n_i*d_i. On remarque aussi que l'Ã©quation (1) ou (2)
  est trivialement vraie dÃ¨s que r*(r-1)>=m (ou r_j*(r_j-1)>=m) car le
  terme de gauche est toujours <= m. On s'arrÃªte donc lorsque
  r>=(1+âˆš(4m+1))/2. C'est important pour Ã©viter de faire un produit
  r*(r-1) qui peut Ãªtre trÃ¨s grand par rapport Ã  m et dÃ©passer la
  capacitÃ© des entiers (et produire un faux nÃ©gatif). Par exemple pour
  k=1, n[0]=10^6 et d[0]=3 (un graphe cubic d'1 million de sommets),
  on a m = 3,000,000 et r*(r-1) â‰ƒ 10^12 ...

*/
{
  if((S==NULL)||(k<=0)) return -1;

  int i,j,m,v;
  for(i=0;i<2*k;i++)
    if(S[i]<0) return -1; /* une valeur <0 => pas graphique */

  /* range les valeurs de S dans n[i] et d[i] */
  NALLOCZ(int,n,k,S[2*_i]);   /* n[0] ... n[k-1]: les n_i */
  NALLOCZ(int,d,k,S[2*_i+1]); /* d[0] ... d[k-1]: les d_i */
  QSORT2(d,n,k,fcmp_int_inv); /* trie les valeurs par ordre dÃ©croissant */

  DEBUG(
	PRINT(k);
	PRINTT(n,k);
	PRINTT(d,k);
	);

  /* fusionne les valeurs de d_i identiques pour avoir une sÃ©quence
     strictement dÃ©croissante */
  for(j=0,i=1;i<k;i++) // i=indice de lecture, j=indice d'Ã©criture
    if(d[i]==d[j]) n[j]+=n[i];
    else j++,n[j]=n[i],d[j]=d[i];
  // ici on a Ã©crit j+1 valeurs dans n[] et d[]
  k=j+1; // rÃ©ajuste k, il n'a pas pu augmenter (mÃªme si k=1)

  DEBUG(
	PRINT(k);
	PRINTT(n,k);
	PRINTT(d,k);
	);

  m=0; // m = somme des degrÃ©s des sommets = âˆ‘_i n_i*d_i
  v=0; // v = valeur retournÃ©e = nombre total de sommets = âˆ‘_i n_i

  for(i=0;i<k;i++){
    m += n[i]*d[i];
    v += n[i];
  }
  
  if(m&1) v=-1; /* somme de degrÃ© impaire => pas graphique */
  else{ /* somme des degrÃ©s paire */
    /* on doit vÃ©rifier, pour chaque j=0..k-2, que: S1 <= P + S2 */
    int rj=0; // rj = âˆ‘_{i=0}^j n[i]
    int S1=0; // S1 = âˆ‘_{i=0}^j n[i]*d[i]
    int S2,P; // S2 = âˆ‘_{i=j+1}^k n[i]*min{r_j,d[i]}, P = r_j*(r_j-1)
    int const rmax=ceil( (1+sqrt(4*((double)m)+1))/2 ); // r>=rmax => r(r-1)>=m
    int const k1=k-1;
    
    for(j=0;j<k1;j++){ // Ã©tape j=0..k-2
      rj += n[j]; if(rj>=rmax) break; // => r_j*(r_j-1)>=m => Ã©quation vraie
      S1 += n[j]*d[j]; P=rj*(rj-1); if(S1<=P) continue; // S1<=P => Ã©quation vraie
      for(S2=0,i=j+1;i<k;i++) S2 += n[i]*min(rj,d[i]); // calcule S2
      if(S1>P+S2){ v=-1; break; } // Ã©quation fausse
    }
  }
  
  free(n);
  free(d);
  return v;
}


/***********************************

         ROUTINES POUR LES
          FONCTIONS adj()

***********************************/


static inline double Norme_dxy(double const dx,double const dy)
/*
  Calcule la longueur du vecteur de coordonnÃ©es (dx,dy) selon la norme
  dÃ©finie par NORM. Renvoie NORM_FAIL si aucune norme n'a Ã©tÃ© trouvÃ©e
  (ce qui en principe ne devrait jamais arriver).

  Principe pour le calcul de NORM_POLY: On considÃ¨re un polygone
  convexe rÃ©gulier centrÃ© sur l'origine de cercle inscrit de rayon
  unitÃ© et orientÃ© de sorte que son cotÃ© le plus Ã  droit soit
  vertical. Ce cotÃ© est numÃ©rotÃ© 0, les autres successivement en
  tournant vers la gauche 1,2,...,(NORM_poly)-1. Notez que
  NORM_poly>=3. Tout d'abord on cherche le numÃ©ro i du cotÃ© du
  polygone qui est coupÃ© par la demi-droite D partant de l'origine et
  passant par le point (x,y). Soit a l'angle de la droite D. L'angle
  infÃ©rieur du cotÃ© i juste en dessous de D est i*p0-p0/2, oÃ¹ p0 est
  l'angle entre deux coins consÃ©cutifs du polygone, c'est-Ã -dire
  i*p0-p0/2 <= a < i*p0+p0/2. On tourne (virtuellement) la figure Ã 
  droite de i cotÃ©s, soit d'un angle de i*p0, de sorte que cotÃ© i soit
  vertical.  Il intersecte alors la droite des abscisses en son milieu
  Ã  l'abscisse 1 prÃ©cisÃ©ment. Le point d'intersection entre D et ce
  cotÃ© est Ã  distance h=1/cos(b) oÃ¹ b=a-i*p0 est l'angle entre D et le
  milieu du cotÃ© i puisque h*cos(b)=1. La distance recherchÃ©e est
  simplement hypot(x,y)/h = hypot(x,y)*cos(b). Notons qu'il n'est pas
  possible d'avoir cos(b)=0 puisque D ne peut Ãªtre parallÃ¨le au i-Ã¨me
  cotÃ©.

  La norme hyperbolique (NORM_HYPER) correspond Ã  la distance
  hyperbolique par rapport Ã  l'origine (0,0). Pour le calcul
  hyperbolique entre deux point quelconques, il faut passer par la
  donction dist_ij() qui dÃ©pend alors des points et pas seulement de
  la diffÃ©rence des points.
*/
{
  switch(NORM){
  case NORM_L1:   return (dx+dy);       /* norme L1 */
  case NORM_L2:   return hypot(dx,dy);  /* norme L2 */
  case NORM_LMAX: return max(dx,dy);   /* norme Lmax */
  case NORM_LMIN: return min(dx,dy);   /* norme Lmin */
  case NORM_POLY:;                      /* norme polygonale */
    double const p0=M_2PI/NORM_poly;    /* NB: ici NORM_poly>=3 */
    double const a=angle(dx,dy); 
    return hypot(dx,dy)*cos(a-p0*floor(a/p0+0.5));
  case NORM_HYPER: return 2*atanh(hypot(dx,dy)); /* norme hyperbolique */
  default: return NORM_FAIL;            /* norme indÃ©terminÃ©e */
  }
}


static inline double dist_ij(query* const Q)
/*
  Calcule la distance entre les points (Q->xpos[Q->i],Q->ypos[Q->i])
  et (Q->xpos[Q->j],Q->ypos[Q->j]) selon la norme dÃ©finie par NORM
  (voir Norme_dxy).  Il faut faire attention que la norme n'est pas
  toujours symÃ©trique en i et j, comme pour NORM_POLY avec p impair.

  Pour la distance hyperbolique, on ne fait pas appel Ã  la norme qui
  n'est dÃ©finie que pour un vecteur depuis l'origine, pas pour une
  diffÃ©rence de points. Aussi les points doivent Ãªtre situÃ©s Ã 
  l'intÃ©rieur du disque unitÃ© centrÃ© sur l'origine, sinon soit il y
  une division par 0 (si un des points est situÃ© sur le cercle) ou la
  valeur renvoyÃ©e peut-Ãªtre nÃ©gative.
*/
{
  if(NORM==NORM_HYPER){
    double const u=Q->xpos[Q->i]*Q->xpos[Q->i]+Q->ypos[Q->i]*Q->ypos[Q->i];
    double const v=Q->ypos[Q->j]*Q->xpos[Q->j]+Q->ypos[Q->j]*Q->ypos[Q->j];
    double const dx=Q->xpos[Q->i]-Q->xpos[Q->j];
    double const dy=Q->ypos[Q->i]-Q->ypos[Q->j];
    return acosh(1+2*(dx*dx+dy*dy)/(1-u*u)/(1-v*v));
  }
  return Norme_dxy(fabs(Q->xpos[Q->i]-Q->xpos[Q->j]),
		   fabs(Q->ypos[Q->i]-Q->ypos[Q->j]));
}


double distgone(const query* const Q,
		int const i,int const p,int const k,double const w)
/*
  Calcule la distance P_i(u,v) entre les points u=(Q->xpos[Q->i],
  Q->ypos[Q->i]) et v=(Q->xpos[Q->j], Q->ypos[Q->j]). Elle n'est pas
  symÃ©trique en u et v. Il s'agit de la "distance p-gone" (un polygone
  rÃ©gulier Ã  p cotÃ©s) relative Ã  la direction i (axe d'angle i*2ðœ‹/k)
  entre les points d'indice u et v, restreint au cÃ´ne de visibilitÃ©
  d'angle w*(p-2)*ðœ‹/p (d'angle w*ðœ‹ si p est infini, c'est-Ã -dire si
  p<3), avec w=0...1 (voir aussi la dÃ©finition du thetagone dans
  l'aide en ligne). Ici, k>0 est le nombre de directions. La fonction
  renvoie DBL_MAX si la distance est infinie ce qui est possible avec
  l'Ã©tape 1 de l'algorithme. L'algorithme est en O(1), notamment
  indÃ©pendant de p et k.

  Soient a_j (j=0...p-1) les sommets du p-gone P_i de rayon unitÃ© avec
  a_0=u et numÃ©rotÃ©s consÃ©cutivement en tournant dans le sense
  positif. Soit c le centre de P_i. Donc dist(u,c)=1, les a_j Ã©tant
  sur un cercle de rayon unitÃ©. On remarque que l'angle (u,a_j) est
  indÃ©pendant du rayon du p-gone, il ne dÃ©pant que de j. En fait,
  l'angle (a_j,u,a_{j+1}) vaut la moitiÃ© de l'angle (a_j,c,a_{j+1}),
  soit ðœ‹/p. L'angle entre deux cotÃ©s consÃ©cutifs d'un p-gone vaut
  (p-2)*ðœ‹/p.

  L'algorithme de calcul pour distgone(u,v) est le suivant:

  1. Trouver la direction j tq (u,v) soit dans le cÃ´ne de visibilitÃ©
     et dans la rÃ©gion [(u,a_j),(u,a_{j+1})[.  Si j n'existe pas,
     alors on renvoit une distance infinie. Si p est infini, a_j est
     simplement sur la droite (u,v).

  2. On calcule l'intersection v' entre les droites (u,v) et
     (a_j,a_{j+1}). Si p est infini, v'=a_j. L'intersection existe
     forcÃ©ment. Eventuellement v est sur la droite (u,a_j).

  3. distgone(u,v)=dist(u,v)/dist(u,v').

*/
{
  int j;
  double xu,xv,dxc,dxa,dxb,dxv;
  double yu,yv,dyc,dya,dyb,dyv;
  double hv,A,Ac,Aw;

  xu=Q->xpos[Q->i],yu=Q->ypos[Q->i]; /* coordonnÃ©es de u */
  xv=Q->xpos[Q->j],yv=Q->ypos[Q->j]; /* coordonnÃ©es de v */
  Ac=(double)i*M_2PI/k; /* angle (u,c), c=centre de P_i */
  dxc=cos(Ac),dyc=sin(Ac); /* coordonnÃ©es du centre c dans le repÃ¨re u */

  dxv=xv-xu;dyv=yv-yu; /* coordonnÃ©es de v dans repÃ¨re u */
  hv=hypot(dxv,dyv); /* |v-u|=dist(u,v) */
  if(hv==0) return 0; /* si u,v ont les mÃªmes coordonnÃ©es */

  /*
    Rappel: Si a et b sont deux vecteurs, alors le produit scalaire
    (dot product) est le rÃ©el a.b = xa*xb + ya*yb = |a|*|b|*cos(a,b),
    oÃ¹ |a|=hypot(xa,ya)=sqrt(xa^2 + ya^2). Notons aussi que
    |a|*|b|*sin(a,b) = det(xa,ya,xb,yb). Donc le signe de cos(a,b) est
    celui de xa*xb + ya*yb. Pour calculer sin(a,b) il faut faire une
    rotation de +ðœ‹/2 au vecteur a et calculer cos(a',b) oÃ¹
    a'=(-ya,xa). Donc sin(a,b) = cos(a',b) = (-ya*xb + xa*yb) /
    (|a|*|b|). Et le signe de sin(a,b) est celui de xa*yb - ya*xb =
    det(xa,ya,xb,yb).
  */

  /* Aw=demi-angle du cÃ´ne de visibilitÃ© */
  Aw=w*M_PI_2; /* si p infini */
  if(p>2) Aw *= (double)(p-2)/p; /* si p est fini */

  /*
    Il faut bien sÃ»r que (u,v) soit dans le cÃ´ne de visibilitÃ©. La
    bissectrice de ce cÃ´ne est l'axe (u,c) et son angle est
    w*(p-2)ðœ‹/p (w*ðœ‹ si p infini). On note (w1,u,w2) le cÃ´ne en
    question. Il faut que (u,v) soit entre (u,w1) (compris) et (u,w2)
    (non compris). Donc si sin(w1,u,v) < 0 ou si sin(w2,u,v) > 0 alors
    il n'existe pas de j (et donc on retourne une distance infinie).
  */

  A=Ac-Aw; /* A=angle (c,w1) */
  dxa=cos(A);dya=sin(A); /* coordonnÃ©es de w1 relatif Ã  u */
  if(det(dxa,dya,dxv,dyv)<0) return DBL_MAX; /* v avant w1 */

  A=Ac+Aw; /* A=angle (c,w2) */
  dxa=cos(A);dya=sin(A); /* coordonnÃ©es de w2 relatif Ã  u */
  if(det(dxa,dya,dxv,dyv)>=0) return DBL_MAX; /* v aprÃ¨s ou sur w2 */

  /*
    Ici v est dans le cÃ´ne de visibilitÃ© et donc la droite (uv)
    intersecte P_i en un point v'.
  */

  Ac -= M_PI; /* Ac=Ac-ðœ‹ */

  /* Cas p infini */
  if(p<3){
    /*
      On raisone dans le repÃ¨re u.  On pose c'=(dyc,-dxc),
      c'est-Ã -dire une rotation de -ðœ‹/2 de (u,v). On a |uc'|=1. On
      calcule l'angle A=(uv,uc'), en fait cos(A). On obtient v' en
      tournant autour de c d'un angle Ac-ðœ‹+2A ...
     */
    A=acos(det(dxv,dyv,dxc,dyc)/hv); // A=acos((dxv*dyc-dyv*dxc)/hv);
    A=Ac+2*A;
    dxa=dxc+cos(A);dya=dyc+sin(A);
    return hv/hypot(dxa,dya);
  }

  /*
    Cas p fini.  On cherche j de sorte qu'en tournant dans le sens
    positif, le vecteur (u,v) soit compris entre (u,a_j) (compris) et
    (u,a_{j+1}) (non compris). La droite (u,v) intersecte le segment
    [a_j, a_{j+1}[. L'indice j recherchÃ© est l'entier tq: (j-1)*ðœ‹/p <=
    angle(a_1,u,a_j) < j*ðœ‹/p. Et donc, j=1+floor{acos(a_1,u,v)/(ðœ‹/p)}.
  */

  double const Ap=M_2PI/p; /* valeur souvent utilisÃ©e */
  A=Ac+Ap; /* angle (c,a_1) */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnÃ©es de a_1 relatif Ã  u */

  /* Aw=cos(a_1,u,v) = (a_1-u).(v-u) / dist(a_1,u)*dist(u,v) */
  Aw=det(dxa,-dya,dyv,dxv)/(hypot(dxa,dya)*hv);
  j=(int)((acos(Aw)*(double)p)/M_PI); /* en fait, la variable j vaut "j-1" */
  A += Ap*j; /* angle (c,a_j): on part de a_1, donc on dÃ©cale de j-1 cÃ´nes */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnÃ©es de a_j relatif Ã  u */
  A += Ap; /* angle (c,a_{j+1}) */
  dxb=dxc+cos(A)-dxa;dyb=dyc+sin(A)-dya; /* vecteur (a_j,a_{j+1}) */

  /*
    Calcule l'unique intersection entre la droite Dv dÃ©finie par le
    vecteur (u,v) et la droite Dj dÃ©finie par le vecteur
    (a_j,a_{j+1}). Dans le repÃ¨re u, les Ã©quations sont:

    Dv: dyv*X - dxv*Y = 0
    Dj: dyb*X - dxb*Y = B avec B=dxa*dyb-dxb*dya

    car a_j=(dxa,dya) appartient Ã  Dj.
    L'intersection (x0,y0) est (dans le repÃ¨re u):
    en faisant dxv*Dj-dxb*Dv, on a: x0=dxv*B/(dxv*dyb-dxb*dyv)=dxv*A
    en faisant dyb*Dv-dyv*Dj, on a: y0=dyv*B/(dxv*dyb-dxb*dyv)=dyv*A
  */
  A=det(dxa,dya,dxb,dyb)/det(dxv,dyv,dxb,dyb);
  return hv/hypot(dxv*A,dyv*A);
}


double func1(double const k,const void* const n)
/*
  Fonction dÃ©finie par:

                   f(k,n) := 2*n/k + (k-1)*(H(k)+1)

  avec H(k) := 1+1/2+1/3+...1/k ~ ln(k)+0.577... + 0.5/k + o(1/k),
  c'est-Ã -dire le k-iÃ¨me nombre harmonic.  Le minimum de cette
  fonction est 2*sqrt(n*ln(n*ln(n))) et est atteint pour un k ~
  0.5*sqrt(2n/ln(n/ln(n)) ce qui est toujours dans l'intervalle [1,n].
*/
{
  return 2.0*(*((int*)n))/k + (k-1.0) * (log(k) + EULER_MASCHERONI + 0.5/k + 1.0);
}


double Minimize(double (*f)(double,const void*),void *info,double ax,double bx,double tol)
/*
  Calcule et renvoie l'abscisse x0 de l'intervalle [ax,bx] du minimum
  de la fonction f(x,info). La valeur de tolÃ©rance "tol" indique la
  prÃ©cision souhaitÃ©e, tol=0 signifiant la prÃ©cision maximale. Pour la
  cherche d'un maximum, il suffit de prendre -f(x,info).
  
  L'algorithme est un mÃ©lange de la recherche selon le nombre d'or
  (afin de minimisÃ© le nombre d'appels Ã  f(x,info) et de
  l'interpolation quadratique (pour mieux approchÃ© la
  solution).

  Pour la recherche selon le nombre d'or, l'idÃ©e est qu'on choisit
  deux points v<w de [a,b] (au dÃ©part a=ax et b=bx), et de rÃ©duire la
  recherche Ã  [a,w] ou Ã  [v,b], intervalles se chevauchant. Si
  f(v)<f(w), alors le minimum ne peut Ãªtre que dans [a..v..w] Ã  cause
  du point v. On recommence donc avec [a,w] sinon avec [v,b]. Un choix
  judicieux de v et w (basÃ© sur le nombre d'or) permet de ne calculer
  f() sur qu'un seul nouveau point.

  Pour l'interpolation quadratique, l'idÃ©e est qu'avec un point x de
  [a,b] calcule la parabole passant par f(a),f(b) et f(x), et on prend
  comme nouveau milieu le point le plus bas de la parabole dans [a,b].
  L'algorithme mÃ©lange les deux techniques.

  Voir le "free software optimize.c" du R's project dans:
  https://svn.r-project.org/R/trunk/src/library/stats/src/

  Notes sur math.h:
  - nextafter(x,y)=plus petit double aprÃ¨s x en direction de y 
  - fma(x,y,z)=x*y+z
  - fdim(x,y)=x-y si x>y, et 0 sinon
*/
{
#define K_MINIMUM        0.3819660112501051518 // = (3-sqrt(5))/2

  double a,b,d,e,p,q,r,u,v,w,x;
  double fu,fv,fw,fx,xm,tol1,tol2,tol3;

  static double eps=-1; /* on calcule eps qu'une seul fois */
  if(eps<0) eps=sqrt(nextafter(0,1)); /* racine carrÃ©e du plus petit double > 0 */
  if(tol<=0) tol=1E-10; /* en dessous, cela ne marche pas toujours !?! */

  a=ax,b=bx;
  w=v=x=fma(K_MINIMUM,b-a,a); // K*(b-a)+a
  fw=fv=fx=f(x,info);

  tol1=nextafter(1,2); /* le plus petit double > 1 */
  tol3=tol/3;
  d=e=0;

  for(;;){
    xm=(a+b)/2; // xm=milieu de [a,b]
    tol1=eps*fabs(x)+tol3;
    tol2=2*tol1;

    /* critÃ¨re d'arrÃªt */
    if(fabs(x-xm)<=tol2-(b-a)/2) break;

    p=q=r=0;
    if(fabs(e)>tol1){ /* fit parabola */
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=(q-r)*2;
      if(q>0) p=-p; else q=-q;
      r=e;
      e=d;
    }

    if((fabs(p)>=fabs(q*r/2))||(p<=q*(a-x))||(p>=q*(b-x))){
      /* Ã©tape: recherche nombre d'or */
      e=(x<xm)? b-x : a-x;
      d=K_MINIMUM*e;
    }
    else{ /* Ã©tape: interpolation quadratique */
      d=p/q;
      u=x+d; 
      /* u ne doit pas trop prÃ¨s de a ou b */
      if((u-a<tol2)||(b-u<tol2)){ d=tol1; if(x>=xm) d=-d; }
    }
    
    /* on va Ã©valuer f en u, qui ne doit pas Ãªtre trop prÃ¨s de x */
    if(fabs(d)<tol1) u=(d>0)? x+tol1 : x-tol1;
    else u=x+d;

    fu=f(u,info); // calcul de f(u)

    /* met Ã  jour a,b,v,w,x puis recommence */
    if(fu<=fx){
      if(u<x) b=x; else a=x;
      v=w; fv=fw;
      w=x; fw=fx;
      x=u; fx=fu;
    }else{
      if(u<x) a=u; else b=u;
      if((fu<=fw)||(w==x)){
	v=w; fv=fw;
	w=u; fw=fu;
      }else
	if((fu<=fv)||(v==x)||(v==w)){ v=u; fv=fu; }
    }
    
  }
  
  return x;
#undef K_MINIMUM
}


enum{
  DYCK_BINOM,
  DYCK_WALK,
  DYCK_WORD,
  DYCK_TREE,
  DYCK_KTREE,
};


int* Dyck(int* R,int const n,int const k,int const code)
/*
  Construit de maniÃ¨re alÃ©atoirement uniforme un mot de k-Dyck composÃ©
  de n segments montant (codÃ© par la lettre k), chacun de longueur k,
  et de kn pas descendant (codÃ© par 0). On doit avoir n>0 et k>0.
  Suivant la valeur de code on peut obtenir des variantes, comme un
  mot binaire de longueur n ayant k valeurs Ã  1 (tirage uniforme).

                                 /\
                              /\/  \
  Exemple de mot de 2-Dyck   /      \
                             2 02 000

  Le rÃ©sultat copiÃ© dans R est donc, en gÃ©nÃ©ral, un mot de n+kn
  lettres composÃ©s de n valeurs k et de kn valeurs 0. La propriÃ©tÃ©
  d'un mot de k-Dyck est que tout prÃ©fixe possÃ¨de au moins autant de
  pas montant que descendant. Le mot classique de Dyck (en bijection
  avec les arbres binaires) sont obtenus pour k=1, une suite de 0 et
  de 1.

  Le rÃ©sultat est renvoyÃ© et mis dans R, qui est allouÃ© si R=NULL.
  Suivant la valeur de code on renvoit des variantes du mot de
  k-Dyck. Attention ! la taille de R, fonction de n et k, varie avec
  code (voir ci-dessous).

    code    taille  rÃ©sultat
  --------------------------------------------------------------------
  DYCK_BINOM n      mot binaire ayant k valeurs Ã  1
  DYCK_WALK  n+nk   mot ayant n valeurs k et kn valeurs 0
  DYCK_WORD  n+nk   mÃªme mot mais dÃ©calÃ© sur un minimum (mot de k-Dyck)
  DYCK_TREE  nk+1   arbre DFS Ã  nk+1 noeuds
  DYCK_KTREE n+nk+1 arbre (k+1)-ary Ã  n+nk+1 noeuds dont n de deg. k+1
  --------------------------------------------------------------------

  Exemple pour n=10 et k=3:
  code=DYCK_BINOM -> 0001010010

  Mot binaire de longueur n avec k valeurs 1. On peut obtenir dans R
  un sous-ensemble alÃ©taoire de {0,...,n-1} Ã  k Ã©lÃ©ments avec:

    int *R=Dyck(NULL,n,k,DYCK_BINOM); // R[] = 0001010010
    for(i=j=0;i<n;i++) if(R[i]) R[j++]=i;
    REALLOC(R,j); // R[] = {3,5,8}

  Exemple pour n=5 et k=1:
  code=DYCK_WALK -> 0101001011
                                            3   4           3   4
                       /\/\                / \ / \           \ /
                    /\/    \/\        1   2   2   2   5     1 2 5
  code=DYCK_WORD -> 1011010010       / \ /         \ / \     \|/
  code=DYCK_TREE -> [-1,0,0,2,2,0]  0   0           0   0     0

  L'algorithme principal, en O(kn), consiste Ã  tirer alÃ©atoirement n
  valeurs k et k*n valeurs 0. On tire k avec une probabilitÃ©
  proportionelle au nombre de valeurs k restant Ã  tirer sur le nombre
  total de valeurs restant Ã  tirer. C'est un tirage alÃ©atoire
  uniforme.

  L'arbre DYCK_TREE est construit de sorte que le mot de Dyck forme le
  parcours DFS de l'arbre: 1 on descend le long d'une arÃªte, et 0 on
  en revient. NB: Si k>1, on fait comme si on avait k pas de 1. Les
  propriÃ©tÃ©s de cet arbre sont nombreuses. Par exemple:
   - R[0] = -1, car 0 est la racine.
   - nk = derniÃ¨re feuille de l'arbre, R[nk] est donc sont pÃ¨re
   - si u-R[u]>1 alors u dÃ©marre une nouvelle branche

  L'arbre (k+1)-ary DYCK_KTREE est obtenu similairement par un DFS
  mais modifÃ© comme suit: on pose d'abord tous les fils avant la
  rÃ©cursion. Voici quelques exemples.

  Exemple pour n=2 et k=2: arbre ternaire Ã  2 sommets internes

                         /\           4 5 6
                      /\/  \           \|/
                     /      \         1 2 3
  code=DYCK_WORD  -> 2 02 000          \|/
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2]   0

  Exemple pour n=3 et k=1: arbre binaire Ã  3 sommets internes
  code=DYCK_WORD  -> 110100
  code=DYCK_KTREE -> [-1,0,0,1,1,4,4]

        5   6
         \ /
      3   4
       \ /
        1   2
         \ /
          0

  Exemple pour n=3 et k=2: arbre ternaire Ã  3 sommets internes
  code=DYCK_WORD  -> 202200000
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2,4,4,4]

      7 8 9
       \|/
        4 5 6
         \|/
        1 2 3
         \|/
          0

  Pour construire l'arbre (k+1)-ary (DYCK_KTREE) on procÃ¨de selon un
  DFS modifiÃ© (on pose les fils avant la rÃ©cursion), en lisant
  succÃ©ssivement les n+kn valeurs du mot de k-Dyck. Plus prÃ©cisÃ©ment,
  depuis le sommet courant u: si on lit k, on ajoute k+1 fils Ã  u, le
  nouveau sommet courant devenant le 1er fils de u. Si on lit 0, on
  remonte les ancÃªtres de u jusqu'Ã  trouver un ancÃªtre v dont le fils
  f menant Ã  u ne soit pas le dernier fils (le (k+1)-iÃ¨me). Le nouveau
  sommet courant est alors f+1. Pour dÃ©terminer si f est le dernier
  fils (ou pas) il suffit de tester si f%(k+1)=0 ou pas.
*/
{
  int const m=n+k*n; // longueur du mot de Dyck
  int *B; // mot de Dyck
  int i,r,s,t;

  if(R==NULL){
    if(code==DYCK_BINOM) ALLOC(R,n);
    if(code==DYCK_WALK)  ALLOC(R,m);
    if(code==DYCK_WORD)  ALLOC(R,m);
    if(code==DYCK_TREE)  ALLOC(R,m-n+1);
    if(code==DYCK_KTREE) ALLOC(R,m+1);
  }
  if((code==DYCK_BINOM)||(code==DYCK_WALK)) B=R;
  else ALLOC(B,m);

  /* DYCK_BINOM: on construit B contenant n-k valeurs 0 et k valeurs
   * 1. DYCK_WALK: on construit B contenant m-n valeurs 0 et n valeurs
   * k, commun Ã  DYCK_WORD, DYCK_TREE et DYCK_KTREE.
   */

  // r = valeurs Ã  Ã©crire (1 ou k)
  // s = nombre total de valeurs Ã  Ã©crire dans B
  // t = nombre de valeurs r Ã  tirer parmi les s

  if(code==DYCK_BINOM) r=1, s=n, t=k;
  else                 r=k, s=m, t=n;

  for(i=0;s>0;s--,i++){ // pour chaque case de B
    if(randomu(s)<t) B[i]=r,t--; // B[i] = r avec proba t/s
    else B[i]=0;
  }
  if((code==DYCK_BINOM)||(code==DYCK_WALK)) return R; // NB: ici R=B
  
  /* cherche la position r dans B de la hauteur minimum */
  /* commun Ã  DYCK_WORD, DYCK_TREE et DYCK_KTREE */

  r=-1;
  int h=0;
  
  for(i=s=0;i<m;i++){
    s += B[i]? k : -1; // s=hauteur courante = +k ou -1
    if(s<h) h=s,r=i;   // h=hauteur minimum
  }
  r=(r+1)%m; // r = position recherchÃ©e
  
  /* DYCK_WORD: dÃ©cale le mot de B et le met dans R */
  
  if(code==DYCK_WORD){
    for(i=0;i<m;i++) R[i]=B[(i+r)%m];
    goto fin_dyck;
  }

  /* DYCK_TREE: dÃ©cale et construit un arbre DFS, R[v]=parent du sommet v */
  /* Si k>1, on fait comme si on avait k pas de 1 */

  if(code==DYCK_TREE){
    int u=0; /* u=dernier sommet visitÃ© (=racine) */
    int v=1; /* v=prochain sommet Ã  visiter (=sommet 1) */
    R[0]=-1; /* pÃ¨re de la racine */
  
    for(i=0;i<m;i++) /* parcoure toutes les valeurs de B */
      if(B[t=(i+r)%m])
	while(B[t]){ /* tantque c'est > 0, on monte */
	  R[v]=u; /* pÃ¨re de v=dernier sommet visitÃ© */
	  u=v++; /* met Ã  jour le dernier sommet visitÃ© */
	  B[t]--;
	}
      else u=R[u]; /* si c'est un 0, on descend */
    goto fin_dyck;
  }

  /* DYCK_KTREE: on a n noeuds internes de degrÃ© k+1 */

  if(code==DYCK_KTREE){
    int u=0,s=1; // u=sommet courant, indexe dans R
    R[0]=-1;     // pÃ¨re de la racine
    for(i=0;i<m;i++) // parcoure toutes les valeurs de B
      if(B[(i+r)%m]){
	for(t=0;t<=k;t++) R[s++]=u; // si on lit k, on pose k+1 fils
	u=s-t; // le nouveau sommet courant est le 1er fils = s-k-1=s-t
      }
      else{
	while(u%(k+1)==0) u=R[u]; // si on lit 0, on remonte tous les derniers fils
	u++; // le nouveau sommet courant est u+1
      }
    goto fin_dyck;
  }

 fin_dyck:
  free(B);
  return R;
}


int NextPermutation(int* const P,int const n,int const *C)
/*
  GÃ©nÃ¨re, Ã  partir d'une permutation P de [0,n[, la prochaine dans
  l'ordre lexicographique suivant les contraintes (intervalles de
  positions) dÃ©finies par le tableau C. On permute les Ã©lÃ©ments de P
  que si leurs positions sont dans l'intervalle [C[t-1],C[t][ pour un
  certain indice t, en supposant que C[-1]=0. Mettre C=NULL s'il n'y a
  pas de contrainte particuliÃ¨re. On renvoie 1 si la prochaine
  permutation a pu Ãªtre dÃ©terminÃ©e et 0 sinon, c'est-Ã -dire si P Ã©tait
  la derniÃ¨re permutation. Dans ce dernier cas la permutation la plus
  petite selon l'ordre lexicographique est renvoyÃ©e, soit P[i]=i.

  On peut allouer et initialiser P avec ALLOCZ(P,k,_i) ou, si le
  tableau P existe dÃ©jÃ , avec NextSet(P,-1,k). Les valeurs de C
  doivent Ãªtre croissantes. Les zÃ©ros en dÃ©but du tableau C n'ont pas
  d'effet, de mÃªme que les valeurs Ã©gales consÃ©cutives. L'algorithme
  est en O(n) mÃªme en moyenne.

  On l'utilise comme ceci:

  Ex1: C=NULL.

    NextSet(P,-1,n);
    do{
        PRINTT(P,n);
	...;
    }while(NextPermutation(P,n,C));

  Ex2: C={2,3,5}, ce qui signifie qu'on ne peut permuter que les
  indices {0,1}{2}{3,4}:

                 0 1 2 3 4 (positions dans P)
	      P={a,b,c,d,e}
	        {b,a,c,d,e}
		{a,b,c,e,d}
		{b,a,c,e,d}
  
  Attention! la permutation doit Ãªtre initialement compatible avec C,
  c'est-Ã -dire une permutation qui peut Ãªtre obtenue Ã  partir de la
  plus petit P[i]=i. Sinon le rÃ©sultat peut ne plus Ãªtre une
  permutation. Par exemple, si P={3,2,1,0,4} et C={3,5}, alors
  NextPermutation(P,5,C) donnera P={0,1,2,4,0}. Pour l'Ã©viter, il
  faudrait mettre un qsort() comme suggÃ©rÃ© dans le source.

  Ã‰videmment, il y a beaucoup moins de permutations dÃ¨s que le nombre
  de contraintes augmente. Par exemple, si C contient k intervalles de
  mÃªme longueur, alors le nombre de permutations sera de (n/k)!^k au
  lieu de n!. Le rapport des deux nombres est d'environ k^n.

  ConcrÃªtement, pour:
  - n=9 et k=3, on a 216 permutations au lieu de 362.880 (k^n=19.683)
  - n=12 et k=3, on a 13.824 permutations au lieu de 479.001.600 (k^n=531.441)

  Le dernier Ã©lÃ©ment de C doit Ãªtre Ã©gale Ã  n (sentinelle), le premier
  Ã©tant omis car il vaut toujours 0. Donc C est un tableau avec au
  plus n Ã©lÃ©ments. Si C=NULL, alors il n'y a pas de contrainte
  particuliÃ¨re, ce qui est identique Ã  poser C[0]=n.

  On se base sur l'algorithme classique appliquÃ© sur les positions i,j
  de l'intervalle [a,b[ = [C[t-1],C[t][ :

  1. Trouver le plus grand index iâˆˆ[a,b-1[ tel que P[i] < P[i+1].
     S'il n'existe pas, c'est que P est la plus grande permutation.
  2. Trouver le plus grand indice jâˆˆ[a,b[ tel que P[i] < P[j].
  3. Ã‰changer P[i] avec P[j].
  4. Renverser la suite de P[i+1] jusqu'au dernier Ã©lÃ©ment P[b-1].

*/
{
  int const T[]={n}; /* sert si C=NULL */
  if(C==NULL) C=T; /* C[0]={n} */

  int i,j,a,b,c;
  b=C[a=j=0]; /* j=indice de la prochaine valeur Ã  lire dans C */

  /* Ã©tape 1: on cherche le plus grand iâˆˆ[a,b-1[ tq P[i]<P[i+1] */
  for(i=-1;;){ /* sentinelle */
    for(c=a;c<b-1;c++) if(P[c]<P[c+1]) i=c; 
    if(i>=0) break; /* ici on a trouvÃ© un i tq P[i]<P[i+1] */
    for(c=a;c<b;c++) P[c]=c; /* on rÃ©initialise P[a]...P[b-1] */
    // qsort(P+a,b-a,sizeof(*P),fcmp_int);
    if(b==n) return 0; /* alors on a fini d'examiner C */
    a=b, b=C[++j]; /* [a,b[=nouvel intervalle */
  }

  /* Ã©tape 2: on cherche le plus grand jâˆˆ[i+1,b[ tq P[i]<P[j] */
  for(j=c=i+1;c<b;c++) if(P[i]<P[c]) j=c;

  /* Ã©tape 3: Ã©change P[i] et P[j] */
  SWAP(P[i],P[j]);

  /* Ã©tape 4: renverse P[i+1]...P[b-1] */
  for(++i,--b;i<b;i++,b--) SWAP(P[i],P[b]);

  return 1;
}


int NextSet(int* const S,int const n,int const k)
/*
  Calcule les sous-ensembles de k entiers de [0,n[. Si n<0, alors S
  est initialisÃ© au plus petit ensemble possible, soit S={0,1,2,
  ...,k-1}. L'idÃ©e est de maintenir une sorte de compteur S qui
  reprÃ©sente le sous-ensemble courant et qui va passer par tous les
  sous-ensembles possibles. Les Ã©lÃ©ments sont toujours rangÃ©s dans
  l'ordre croissant. On renvoie 1 si on a pu construire le prochain
  sous-ensemble, et 0 si S Ã©tait le dernier sous-ensemble. Dans ce cas
  l'ensemble le plus petit est Ã©crit dans S. On l'utilise comme ceci:

  NextSet(S,-1,k);
  do{
    ...;
    PRINTT(S,k); // traitement de S
    ...;
  }while(NextSet(S,n,k));

  On peut facilement en dÃ©duire un algorithme pour gÃ©nÃ©rer des
  multi-ensembles, comme par exemple tous les multi-ensembles du type
  [2,0,0,2,1,0] comprennant 3 fois 0, 1 fois 1 et 2 fois 2:
  NextMultiSet(S,C,k) avec C=[3,1,2] (voir la fonction ggosset()).

  La stratÃ©gie pour "incrÃ©menter" S est la suivante : on essaye
  d'incrÃ©menter S[i] (au dÃ©part i=0) tout en restant strictement
  infÃ©rieur Ã  l'Ã©lÃ©ment suivant S[i+1]. Si cela marche on a trouvÃ© le
  prochain sous-ensemble, sinon on pose S[i]=i et on recommence avec
  indice i suivant. Si S[k-1] atteint n c'est que S Ã©tait le dernier
  sous-ensemble.

  L'algorithme est en O(k) dans le pire des cas, mais de maniÃ¨re
  amortie c'est beaucoup moins car on incrÃ©mente moins souvent S[j]
  que S[i] si j>i.
*/
{
  int i=0,j,s;

  if(n<0){
    for(;i<k;i++) S[i]=i;
    return 1;
  }
  
  while(i<k){
    s=++S[j=i++];
    if(i==k){ if(s<n) return 1; }
    else{ if(s<S[i]) return 1; }
    S[j]=j;
  }
  
  return 0;
}


int NextArrangement(int* const S,int* const P,int const n,int const k)
/*
  Permet de gÃ©nÃ©rer tous les arrangements de k entiers de
  [0,n[. L'arrangement est reprÃ©sentÃ© par les tableaux S et P de k
  entiers. S reprÃ©sente un ensemble de k entiers et P une permutation
  de [0,k[. Ainsi, l'arrangement A=(4,2,7,3) est reprÃ©sentÃ© par
  S=(2,3,4,7) et P=(2,0,3,1). Autrement dit A[i]=S[P[i]] pour tout
  i=0...k-1.

  L'idÃ©e est d'essayer d'incrÃ©menter le double compteur S,P. On
  incrÃ©mente P en premier avec NextPermutation(). Si on est arrivÃ© Ã 
  la fin de P, on incrÃ©mente S avec NextSet(). Si n<0, alors S et P
  sont initialisÃ©s au plus petit arrangement possible, soit S = P =
  (0,1,2, ...,k-1). On renvoie 1 si on a pu trouver le prochain
  arrangement, 0 si c'Ã©tait le dernier arrangement possible. Dans ce
  cas l'arrangement le plus petit est Ã©crit dans S,P.
*/
{
  if(n<0){
    int i;
    for(i=0;i<k;i++) S[i]=P[i]=i;
    return 1;
  }
  if(!NextPermutation(P,k,NULL)) return NextSet(S,n,k);
  return 1;
}


int* NextPart(int* S,int const n,int s,int* const C)
/*
  Permet de gÃ©nÃ©rer toutes les suites S de n>0 entiers >=0 dont la
  somme fait s et dont la i-Ã¨me part S[i] ne dÃ©passe pas C[i]. Il faut
  que s <= âˆ‘_{i=0}^{n-1} C[i], n=1 est possible, de mÃªme que C[i]=s.

  Initialement S est la suite courante de somme s et on renvoie dans S
  la prochaine suite (la fonction renvoie aussi S). On renvoie NULL si
  on a atteint la derniÃ¨re suite, et on remplit S avec le premiÃ¨re
  suite. Si S=NULL, alors S est allouÃ©e et initialisÃ©e Ã  la premiÃ¨re
  suite. La premiÃ¨re suite de somme s est obtenue en remplissant
  autant que possible les parts S[n-1],S[n-2],...

  L'algorithme est le suivant:
   1. on dÃ©termine le plus grand indice j tq S[j]>0
   2. si j=0, alors on a finit: on va Ã  l'Ã©tape 6 avec x=s+1 et i=-1
   3. on dÃ©termine le plus grand indice i<j tq S[i] peut Ãªtre augmentÃ©
   4. on calcule x = âˆ‘_{j=i+1}^{n-1} S[i]
   5. on incrÃ©mente S[i]
   6. on remplit S[i+1]...S[n-1] avec la premiÃ¨re suite de somme x-1

  On l'utilise comme ceci:

  int s=n/2;
  NALLOCZ(int,C,n,1); // compteur binaire avec au plus n/2 valeurs 1
  int *S=NextPart(NULL,n,s,C); // initialisation de S Ã  0
  do{ PRINTT(S,n); // traitement de la partition S
      ...;
  }while(NextPart(S,n,s,C));

  Exemple: s=n=5

  C = 1 2 2 1 1
  S = 0 1 2 1 1
      0 2 1 1 1
      0 2 2 0 1
      0 2 2 1 0
      1 0 2 1 1
      1 1 1 1 1
      1 1 2 0 1
      1 1 2 1 0
      1 2 0 1 1
      1 2 1 0 1
      1 2 1 1 0
      1 2 2 0 0
*/
{
  int x,i,j,r;

  i=0;
  r=(S==NULL);
  if(r) ALLOC(S,n);
  else i=n-1;
  
  /* calcule le plus grand indice i tq S[i]>0 */
  while((i>0)&&(S[i]==0)) i--;
  x=S[i--]; /* rem: si i=0, alors i -> -1 */

  /* calcule le plus grand indice j<i tq S[j] peut Ãªtre augmentÃ© */ 
  while((i>=0)&&(S[i]==C[i])) x += S[i--];

  if(i>=0){ S[i]++; s=x-1; } /* si cet indice n'existe pas i<0 => FIN */
  
  /* Ã©crit la premiÃ¨re suite de somme s dans S[i+1]...S[n-1] */
  for(j=n-1;j>i;j--){
    x=max(s,0);
    x=min(C[j],x);
    S[j]=x;
    s -= x;
  }

  /* on retourne S sauf si i<0 et r=0 (<=> FIN ) */
  return ((i<0)&&(!r))? NULL : S;
}


int* NextIntPartition(int* S,int const n,int* const t)
/*
  Permet de gÃ©nÃ©rer toutes les partitions d'un entier n>0 dans l'ordre
  lexicographique croissant. L'entier retournÃ© dans *t est le nombre
  de parts de la partition renvoyÃ©e. Il est possible de mettre t=NULL
  et dans ce cas le nombre de parts n'est pas renvoyÃ©.  Si S=NULL,
  alors S est allouÃ©e (de taille n) et initialisÃ©e Ã  la premiÃ¨re
  partition de n, soit S={1,...,1}. Dans tous les cas, S est renvoyÃ©e.
  Si S valait la derniÃ¨re partition de n, soit S={n}, alors NULL est
  renvoyÃ© et *t=1.

  On l'utilise comme ceci:

  int t; // t=pour connaÃ®tre le nombre de parts
  int *S=NextIntPartition(NULL,n,&t); // initialisation de S
  do{ PRINTT(S,t); // traitement de la partition S en t parts
      ...;
  }while(NextIntPartition(S,n,&t));

  S = 1 1 1 1 1 1
  S = 1 1 1 1 2
  S = 1 1 1 3
  S = 1 1 2 2
  S = 1 1 4
  S = 1 2 3
  S = 1 5
  S = 2 2 2
  S = 2 4
  S = 3 3
  S = 6

  La complexitÃ© en temps est constante en moyenne. L'algorithme est
  basÃ© sur http://jeromekelleher.net/category/combinatorics.html.
*/
{
  static int k=-1; // k en statique pour aller plus vite
  if((S==NULL)||(k<0)){ // premiÃ¨re fois ?
    if(S==NULL) ALLOC(S,n); // alloue S
    S[k=1]=n; S[0]=0; // pour initialiser S
  }
  
  if(k==0) return NULL; // derniÃ¨re part atteinte
  int y=S[k--]-1;
  int x=S[k]+1;
  while(x<=y) S[k++]=x, y-=x;
  S[k]=x+y; 
  if(t) *t=k+1; // nombre de parts
  return S;
}


int *RandomIntPartition(int* S,int const n)
/*
  Produit une partition alÃ©atoire uniforme de l'entier n>0. Le tableau
  S doit Ãªtre de taille n au moins. Si S=NULL, il est allouÃ© et
  renvoyÃ©. En retour S[i] reprÃ©sente le nombre de fois oÃ¹ la valeur
  i+1 est dans la partition. La complexitÃ© moyenne de l'algorithme est
  en O(n^1.25).

  L'algorithme Ã  rejet ci-dessous calcule une partition S =
  (Z_1,...,Z_n) de n, Z_i Ã©tant le nombre de fois oÃ¹ l'entier i>0
  apparaÃ®t dans la partition:
  
  1. Soit X=-ðœ‹/âˆš(6n).

  2. GÃ©nÃ©rer n-1 variables indÃ©pendantes Z_2,...,Z_n, oÃ¹ Z_i suit une
     loi gÃ©omÃ©trique de paramÃ¨tre p_i=1-exp(X)^i, donnant le nombre
     d'Ã©checs avant le 1er succÃ¨s de probabilitÃ© p_i. Pour gÃ©nÃ©rer une
     telle variable il suffit de faire Z_i = floor{ ln(U) / ln(1-p_i)
     } = floor{ ln(U) / (i*X) } oÃ¹ U=RAND01. Poser Z_1 = n - âˆ‘_{i=2}^n
     i*Z_i.

  3. Si Z_1<0 ou ln(RAND01)>Z_1*X, recommencer en 2.

  D'aprÃ¨s [AD15] le nombre de rÃ©pÃ©titions moyen est au plus
  2*ðœ‹*(6n)^(1/4) < 10*n^0.25. On peut faire O(1) rÃ©pÃ©titions, mais
  c'est bien plus compliquÃ©. En fait, l'algorithme pourrait Ãªtre
  accÃ©lÃ©rÃ© car, avec grande probabilitÃ© (cf. Erdos-Lehner 1941), la
  plus grande part est proche de 2c*âˆšn*log(n) oÃ¹ c=ðœ‹/âˆš6 et donc Z_i=0
  lorsque i>>âˆšn*log(n). Donc on pourrait utiliser un seul random()
  pour dÃ©terminer Ã  partir de quelle valeur de i on met tous les
  Z_i=0. Il est montrÃ© dans [AD15, pp. 22] qu'il est possible de faire
  O(âˆšn) appels Ã  random() pour gÃ©nÃ©rer tous les Z_i.

  Une mÃ©thode qui gÃ©nÃ¨re aussi une partition alÃ©atoire, mais pas
  uniformÃ©ment, consiste Ã  tirer un tableau S alÃ©atoire de taille n
  dont les valeurs sont entiÃ¨res et >=0 et dont la somme fait n ce qui
  peut se faire avec la ligne de code suivante: ALLOCZ(S,n,0);
  for(i=0;i<n;i++) S[randomu(n)]++;
*/
{
  if(S==NULL) ALLOC(S,n);
  double const x=-M_PI/sqrt(6*n);
  int i,s; // s=S[0]=Z_1
  double p;

  do{
    s=n; i=1; p=x;
    while((s>=0)&&(i<n)){
      p += x; // au dÃ©but p=2*x
      S[i]=(int)(log(RAND01)/p); // S[i] = Z_{i+1} = ln(U)/((i+1)*x)
      s -= (i+1)*S[i];
      i++;
    }
  }while((s<0)||(log(RAND01)>x*s));

  S[0]=s;
  return S;
}


int SetCmp(int *T1,int *T2,int n1,int n2)
/*
  Compare deux tableaux d'entiers T1 et T2 de taille n1 et n2 triÃ©s
  par ordre croissant. Les tableaux peuvent Ãªtre de taille nulle. La
  valeur renvoyÃ©e est un entier interprÃ©tÃ© en binaire comme suit:

  bit-0: 1 ssi T1 intersecte T2 (possÃ¨de au moins 1 Ã©lÃ©ment commun)
  bit-1: 1 ssi T1 Ã©gale T2
  bit-2: 1 ssi T1 est strictement inclu dans T2
  bit-3: 1 ssi T2 est strictement inclu dans T1

  Les valeurs possibles sont donc: 0,1,2,3,4,5,8,9 (=[0,9]\{6,7})
  La complexitÃ© est O(n1+n2).
*/
{
  if(n1==0) return (n2==0)? 2:4;
  if(n2==0) return 8;
  /* ici T1 et T2 contiennent au moins 1 Ã©lÃ©ment */

  if((T1[n1-1]<T2[0])||(T2[n2-1]<T1[0])) return 0; /* cas trivial de disjonction */

  int i1,i2,r;
  i1=i2=0;
  r=14; /* tous les bit Ã  1 sauf b0 */

  while((i1<n1)&&(i2<n2)){
    if(T1[i1]==T2[i2]){
      i1++; i2++; /* T1 et T2 ont un Ã©lÃ©ment en commun */
      r |= 1; continue; /* met b0 */
    }
    r &= 13; /* annule b1 (15-2) car T1<>T2 */
    if(T1[i1]<T2[i2]){
      i1++; /* T1 contient des Ã©lÃ©ments qui ne sont pas dans T2 */
      r &= 11; /* annule b2 (15-4) car T1 ne peux pas contenir T1 */
    }else{
      i2++; /* T2 contient des Ã©lÃ©ments qui ne sont pas dans T1 */
      r &= 7; /* annule b3 (15-8) car T2 ne peux pas contenir T1 */
    }
  }

  if(i1<n1) r &= 9; /* annule b2 et b1 (15-4-2) */
  if(i2<n2) r &= 5; /* annule b3 et b1 (15-8-2) */
  if(r&2)   r &= 3; /* annule b3 et b2 (15-8-4) */

  return r;
}


int SetSearch(int const u,const int* const T,int const n,int const sort)
/*
  Cherche l'entier u dans le tableau T de taille n. Si u est dans T,
  on renvoie son indice sinon on renvoie -1. Si sort=1, alors on
  suppose T triÃ© par ordre croissant et la recherche est dichotomique,
  sinon la recherche est linÃ©aire.
*/
{
  if(sort){ /* recherche dichotomique */
    int *t=bsearch(&u,T,n,sizeof(int),fcmp_int);
    return t? (int)(t-T) : -1;
  }
  /* recherche linÃ©aire */
  int i;
  for(i=0;i<n;i++) if(u==T[i]) return i;
  return -1;
}


int Binom(int const n,int const k)
/*
  Calcule l'entier B={n choose k}.  L'algorithme utilisÃ© ici est en
  O(k). Il n'utilise que des multiplications et divisions entiÃ¨res sur
  des entiers en O(B), sans aucun tableau.

  L'algorithme classique issu du Triangle de Pascal, qui lui est en
  O(n*k), utilise un espace en O(k) (tableaux d'entiers en O(B)). Par
  contre il n'utilise que des additions sur des entiers en O(B).

  Principe: B = n x (n-1) x ... x (n-k+1) / k x (k-1) x ... x 1

  On rÃ©ecrit B en (((((n/1) x (n-1)) / 2) x (n-2)) / 3) ...
  c'est-Ã -dire en multipliant le plus grand numÃ©rateur et en divisant
  par le plus petit numÃ©rateur. Chaque terme est ainsi un certain
  binomial, et donc toujours un entier.

  Catalan(n) = Binom(2n,n)/(n+1). On peut aussi calculer ce nombre
  avec Catalan(0)=1 et Catalan(n) = (2*(2n-1)*Catalan(n-1)/(n+1).
*/
{
  int B,s,i;
  for(B=s=n,i=2;i<=k;i++) B=(B*(--s))/i;
  return B;
}


/* code pour la fonction SortInt() */

enum{
  SORT_INC,     /* tri croissant selon les valeurs */
  SORT_DEC,     /* tri dÃ©croissant selon les valeurs */
  SORT_INC_RANK,/* donne le rang des valeurs de T dans l'ordre croissant */
  SORT_DEC_RANK,/* donne le rang des valeurs de T dans l'ordre dÃ©croissant */
  SORT_FREQv,   /* donne la frÃ©quence des valeurs de T */
  SORT_FREQe,   /* donne la frÃ©quence des Ã©lÃ©ments de T */
  SORT_INDEXi,  /* donne l'indice des Ã©lÃ©ments de T dans l'ordre croissant */
  SORT_INDEXe   /* donne les Ã©lÃ©ments de T dans l'ordre croissant */
};


int *SortInt(int* T,int* R,int const n,int const a,int* const m,int const code)
/*
  Trie par ordre croissant un tableau T non NULL de taille n>0 dont
  les valeurs sont des entiers de [a,a+*m[ si m<>NULL. Sinon, les
  valeurs de T sont supposÃ©es Ãªtre dans [a,a+n[. Pour simplifier le
  texte qui suit, je note m pour dire *m.

  La complexitÃ© en temps est O(n+m), ce qui est mieux que qsort(). Le
  tableau T n'est pas modifiÃ©. Le rÃ©sultat est rangÃ© dans le tableau R
  qui doit Ãªtre de taille au moins n, et est aussi renvoyÃ© par la
  fonction. Il R=NULL, alors R est allouÃ© et retournÃ© par la
  fonction. Pour Ãªtre compÃ©titif avec qsort() pour m=n il faut,
  thÃ©oriquement, n>32.

  L'opÃ©ration de tri dÃ©pend de "code" qui est interprÃ©tÃ©e comme suit
  (la complexitÃ© est le nombre d'Ã©tapes dans les boucles "for"):

    v = valeur d'un Ã©lÃ©ment de T, v dans [a,a+m[
    d = nombre de valeurs v distinctes de T, d dans [1,min(n,m)]
    e = indice du tableau T, e dans [0..n[
    i = indice du tableau R, i dans [0..n[
    r = rang dans un tableau, r dans [0..d[

  - code=SORT_FREQv: renvoie un tableau F[0..m[ de frÃ©quence des
    valeurs de T, c'est-Ã -dire que F[i] est le nombre de fois oÃ¹ la
    valeur i+a apparaÃ®t dans T. Le tableau R n'est pas utilisÃ© (et
    peut donc Ãªtre NULL) et la variable m n'est pas
    modifiÃ©e. ComplexitÃ©: m+n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        F = [0 1 3 2 0 1 1] de taille m

  - code=SORT_FREQe: renvoie dans R[0..n[ un tableau de frÃ©quence des
    Ã©lÃ©ments de T, c'est-Ã -dire oÃ¹ R[e] est le nombre de fois oÃ¹ T[e]
    apparaÃ®t dans T. La variable m n'est pas modifiÃ©e. ComplexitÃ©:
    m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 3  1  3  2  3  1  1  2] de taille n

  - code=SORT_INC ou SORT_DEC: renvoie dans R[0..n[ le tableau T triÃ©
    par ordre croissant (ou dÃ©croissant). La variable m n'est pas
    modifiÃ©e. ComplexitÃ©: 2m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [11 12 12 12 13 13 15 16] de taille n

  - code=SORT_INC_RANK ou SORT_DEC_RANK: renvoie dans R[0..n[ un
    tableau de rangs oÃ¹ R[e] est le rang r dans [0..d[ de l'Ã©lÃ©ment
    T[e] dans la version triÃ©e dans l'ordre croissant (ou dÃ©croissant)
    de T. Le tableau R est modifiÃ© et on renvoie d dans m. ComplexitÃ©:
    3m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 1  0  1  2  1  3  5  2] et m=5

  - code=SORT_INDEXi: renvoie dans R[0..n[ un tableau d'indices oÃ¹
    R[i]=e est l'Ã©lÃ©ment de T en position i dans la version triÃ©e par
    ordre croissant de T. Pour obtenir un tri de T il suffit de lister
    T[R[i]] pour i=0..n-1. La variable m n'est pas modifiÃ©e.
    ComplexitÃ©: 2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 1  7  0  2  4  3  6  5]

  - code=SORT_INDEXe: renvoie dans R[0..n[ un tableau d'indices oÃ¹
    R[e] est la position de T[e] dans la version triÃ©e par ordre
    croissant de T. La variable m n'est pas modifiÃ©e. ComplexitÃ©:
    2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 2  0  3  5  4  7  6  1]
*/
{
  int i,r,t;
  int k=(m==NULL)? n:*m;

  /* initialise F[0..m[ */
  NALLOCZ(int,F,k,0); /* coÃ»t: m */
  
  /* calcule F[i]=frÃ©quence de la valeur v=i+a dans T */
  for(i=0;i<n;i++) F[T[i]-a]++; /* coÃ»t: n */

  if(code==SORT_FREQv) return F;

  /* alloue R, si nÃ©cessaire */
  if(R==NULL) ALLOC(R,n);

  if(code==SORT_FREQe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coÃ»t: n */
    free(F); return R;
  }

  if(code==SORT_INC){ /* R=tableau T triÃ©, ordre croissant */
    for(i=r=0;i<k;i++) for(t=F[i];t>0;t--) R[r++]=i+a; /* coÃ»t: m+n */
    free(F); return R;
  }

  if(code==SORT_DEC){ /* R=tableau T triÃ©, ordre dÃ©croissant */
    for(i=0,r=n;i<k;i++) for(t=F[i];t>0;t--) R[--r]=i+a; /* coÃ»t: m+n */
    free(F); return R;
  }

  /* calcule F[i]=nombre de valeurs de T qui sont < i+a.     
     Ex:  T = [2 1 2 3 2 8 5]  avec n=7, m=9, a=0
	  F = [0 0 1 4 5 5 6 6 6]
  */
  for(i=r=0;i<k;i++){ t=F[i]; F[i]=r; r += t; } /* coÃ»t: m */

  if(code==SORT_INDEXi){
  for(i=0;i<n;i++) R[F[T[i]-a]++]=i; /* coÃ»t: n */
    free(F); return R;
  }

  if(code==SORT_INDEXe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]++; /* coÃ»t: n */
    free(F); return R;
  }

  /* pour SORT_INC_RANK ou SORT_DEC_RANK */
  /* calcule dans F[i] le rang dans [0,d[ du nb de valeurs < i+a */
  
  for(t=r=-1,i=0;i<k;i++){ /* coÃ»t: m */
    if(F[i]!=t) { t=F[i]; r++; }
    F[i]=r;
  }
  if(m!=NULL) *m=r+1; /* m=nb de valeurs diffÃ©rentes de T */

  if(code==SORT_INC_RANK){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coÃ»t: n */
    free(F); return R;
  }
  
  if(code==SORT_DEC_RANK){
    for(i=0;i<n;i++) R[i]=r-F[T[i]-a]; /* coÃ»t: n */
    free(F); return R;
  }

  free(F);
  free(R);
  Erreur(21);
  return NULL;
}


enum{
  SUBDIV_UNIF=0, 
  SUBDIV_RAND,
  SUBDIV_RAND1,
};


graph* subdivision(graph* H,int const n,int const code)
/*
  Renvoie une subdivision G du graphe H. La subdivision dÃ©pend de la
  valeur du code. On renvoie une copie de H si H est sans arÃªte.

  code=SUBDIV_UNIF: subdivision uniforme de chacune des arÃªtes en n
  nouveaux sommets.

  code=SUBDIV_RAND: la subdivision aura un total de n nouveaux sommets
  rÃ©partis alÃ©atoirement parmi sur les arÃªtes de H.

  code=SUBDIV_RAND1: comme SUBDIV_RAND mais chaque arÃªte aura au moins
  un nouveau sommet, il faut donc n>=H->m.
*/
{
  int const m=nb_edges(H); // calcule le nombre d'arÃªtes de H
  if(m<1) return GraphCopy(H); // rien Ã  faire

  int i,j,e,u,v,s,t;
  int* R;

  // Construit un tableau R tq R[i] soit le nombre de fois dontoÃ¹ l'arÃªte
  // i doit Ãªtre subdivisÃ©e, i=0..m-1

  if(code==SUBDIV_UNIF){ // subdivisions uniformes
    ALLOCZ(R,m,n); // R[i]=n
    u=H->n+n*H->m; // nombre de sommets du graphe final G
  }
  else{ // subdivisions alÃ©atoires
     // On tire alÃ©atoirement m-1 valeurs 1 dans un mot binaire R de
     // longueur t=n+m-1 ou t=n, ce qui va crÃ©er m blocs. Pour le
     // dernier bloc, on fait comme s'il y avait un dernier 1 Ã  la fin
     // de R. On rÃ©Ã©crit ensuite au dÃ©but de R, dans R[i], la longueur
     // du bloc i.

     // Ex: SUBDIV_RAND:
     // m=5, n=7, |R|=11, R = 00101100010 = [000 00 0 0000 00] -> R = 3,2,1,4,2
     // R[i] est le nombre de 0 avec le 1 dans chaque bloc

     // Ex: SUBDIV_RAND1:
     // m=5, n=6, |R|=10, R = 0010110001 = [00|0||000|] -> R = 2,1,0,3,0
     // R[i] est le nombre de 0 dans chaque bloc

    t=-1; // sentinelle
    if(code==SUBDIV_RAND)   t=n+m-1,s=1;
    if(code==SUBDIV_RAND1){ t=n, s=0; if(n<m) return GraphCopy(H); }
    if(t<0) Erreur(50); // y'a un problÃ¨me

    R=Dyck(NULL,t,m-1,DYCK_BINOM); // tirage uniforme de m-1 valeurs 1 parmi t
    for(e=i=j=0;i<t;i++) // calcule la longueur des blocs
      if(R[i]) R[j++]=i-e,e=i+s; // un bloc Ã  chaque 1 trouvÃ©
    R[j]=i-e; // pour le dernier bloc (comme s'il y avait eut un 1)
    u=H->n+n; // // nombre de sommets du graphe final G
  } 

  // Etant donnÃ© le tableau R, construit une subdivision de H oÃ¹
  // l'arÃªte i est subdivisÃ©e R[i] fois.

  graph* G=new_fullgraph(u); // le graphe final

  t=H->n; // t = prochain sommet Ã  crÃ©er
  s=t-1;  // s = dernier sommet crÃ©Ã©
  e=0;    // e = indice de l'arÃªte dans R

  for(u=0;u<H->n;u++){ // parcoure toutes les arÃªtes de H
    for(i=0;i<H->d[u];i++){	
      v=H->L[u][i]; // ici (u,v) est une arÃªte de H
      if(u>=v) continue; // pour garantir que u<v
      s=u; // sommet d'origine
      for(j=0;j<R[e];j++,t++){ // pour chaque nouveau sommet
	ADD_EDGE(G,s,t);
	s=t;
      }
      ADD_EDGE(G,s,v); // connexion Ã  v
      e++; // arÃªte suivante
    }
  }
  free(R);
  GraphRealloc(G,G->d);
  return G;
}


/***********************************

       BFS, DFS, ...
       (DIJKSTRA)

***********************************/


/*
  Structure de donnÃ©es pour le BFS
*/
typedef struct{
  int root;  /* racine ou source du BFS. */
  int radius;/* eccentricitÃ© du sommet source. */
  int *D;    /* D[u]=distance entre u et root. Les sommets u avec
	        D[u]=-1 sont Ã  une distance infinie de root (situation
	        initiale par dÃ©faut). En entrÃ©e, les sommets avec
	        D[u]=-2 sont considÃ©rÃ©s comme inexsitant dans le
	        graphe. Si D=NULL, D est allouÃ© puis initialisÃ© Ã 
	        -1. */
  int *P;    /* P[u]=pÃ¨re de u dans un arbre BFS de racine root, avec
	        P[root]=-1. Seuls les sommets u<>root avec D[u]>=0 ont
	        un pÃ¨re dÃ©fini. Si P=NULL, alors P est allouÃ©. Il
	        n'est pas nÃ©cessaire de l'initialiser. */
  int n;     /* nombre de sommets parcourus = nombre d'Ã©lÃ©ments dans
		la file */
  int *file; /* contenu de la file = liste des n sommets parcourus. La
		taille de ce tableau est toujours le nombre de sommets
		du graphe. */
  int cycle; /* longueur du plus petit cycle passant par la source
		rencontrÃ© lors du parcours. Cette valeur (>2) est
		indÃ©pendente du parcours spÃ©cifique des voisins (ie de
		tout BFS). Elle ne dÃ©pend que de la structure
		non-Ã©tiquetÃ©e du graphe et de la source. Si cycle<0,
		alors cette longueur est infinie.  On peut dÃ©terminer
		la maille du graphe en prennant le minimum obtenu pour
		chacun des sommets du graphe comme source. Cette
		valeur n'est pas correcte si le graphe est orientÃ©. */
  int clean; /* permet d'initialiser Ã  -1 les sommets parcourus dans
		le prÃ©cÃ©dent bfs(). Par dÃ©faut, clean=0, et le tableau
		D n'est pas initialisÃ© (sauf si D=NULL, dans ce cas il
		est initialisÃ© complÃ¨tement Ã  -1). Si clean=2, alors
		on remet D[u] Ã  -1 seulement pour les sommets u de
		file. C'est plus rapide si bfs() est utilisÃ© plusieurs
		fois avec le mÃªme paramÃ¨tre et tous les sommets ne
		sont pas visitÃ©s (vmax ou hmax >0). Si clean=1, alors
		on force l'initialisation complÃ¨te de D Ã  -1, puis on
		passe clean=2 (pour les appels suivants). */
  int tf;    /* tÃªte de file, ie nombre d'Ã©lÃ©ments parcourus dont tous
		les voisins ont Ã©tÃ© enfilÃ©s (tf<=n). */
  int cont;  /* cont=1 si on poursuit le bfs() prÃ©cÃ©dant lÃ  oÃ¹ on
		s'Ã©tait arrÃªtÃ© dans le cas d'un arrÃªt par hmax. Par
		dÃ©faut cont=0. AprÃ¨s le 1er bfs() oÃ¹ cont=1, cont
		passe Ã  2. Cela sert Ã  augmenter progressivement la
		hauteur jusqu'Ã  obtenir une certaine condition. */
  int vmax;  /* arrÃªte le parcours lorsque vmax sommets on Ã©tÃ© parcourus */
  int hmax;  /* arrÃªte le parcours lorsqu'un sommet de hauteur > hmax est atteint */
} param_bfs;


param_bfs *new_param_bfs(void)
/*
  CrÃ©e et initialise une structure pour la fontion bfs(). C'est la
  fonction bfs() qui doit se charger, Ã©ventuellement, d'allouer les
  tableaux P et D (qui sont ici initialisÃ© Ã  NULL). Le tableau D peut
  Ãªtre utilisÃ© pour effacer des sommets.
*/
{
  NALLOC(param_bfs,X,1);
  X->D=X->P=X->file=NULL;
  X->radius=X->n=X->clean=X->cont=X->tf=0;
  X->cycle=X->root=X->vmax=X->hmax=-1;
  return X;
}


void free_param_bfs(param_bfs* const X)
{
  if(X==NULL) return;
  free(X->D);
  free(X->P);
  free(X->file);
  free(X);
}


param_bfs *bfs(const graph* const G,int const source,param_bfs* X)
/*
  Effectue un parcours en largeur (BFS) d'un graphe G (orientÃ© ou non)
  depuis le sommet source. Les rÃ©sultats du parcours (comme les
  distances Ã  la source, le pÃ¨re, le nombre de sommets parcourus,
  etc.) sont stockÃ©es dans la variable X qui est renvoyÃ©e. Si X=NULL,
  X est d'abord allouÃ©e.

  On peut Ã©galement rÃ©aliser un BFS seulement sur un sous-ensemble de
  sommets (sous-graphe), limiter le parcours Ã  une profondeur donnÃ©e
  (en fixant X->hmax), ou Ã  un nombre de sommets parcourus (en fixant
  X->vmax). Il est important que la table X->D soient correctement
  initialisÃ©e pour les sommets Ã  parcourir. Les tables X->P et X->file
  ne sont jamais initialisÃ©s. Si elles sont NULL, elles sont
  simplement allouÃ©es. Seules les valeurs des X->n sommets parcourus
  (ceux dans X->file) sont garanties d'Ãªtre correctes. Ã€ noter que le
  parcours est lancÃ© depuis la source, mÃªme si source est dans T.

  On l'utilise comme ceci:

    param_bfs *X=bfs(G,s,NULL); // BFS depuis le sommet s dans G
    ... // X->D[u]=distance entre s et u
    ... // X->P[u]=pÃ¨re de u, ou -1 s'il n'existe pas
    ... // X->cycle=longueur du plus petit cycle passant par s
    ... // X->n=nombre de sommets parcourus
    ... //
    free_param_bfs(X); // libÃ¨re la variable crÃ©e X

  Pour rÃ©aliser un BFS d'un sous-graphe G\T de G (Ã©vitant les k
  sommets du tableau T):

    param_bfs *X=new_param_bfs();   // par dÃ©faut X->clean=0
    ALLOCZ(X->D,G->n,-1);           // alloue X->D et met tout Ã  -1
    for(i=0;i<k;i++) X->D[T[i]]=-2; // les sommets Ã  enlever doivent Ãªtre Ã  -2
    bfs(G,s,X);
    ... // X->D[u]=distance entre s et u dans G\T
    ... // X->P[u]=pÃ¨re de u, ou -1 s'il n'existe pas
    ... // X->cycle=longueur du plus petit cycle dans G\T passant par s
    ... // X->n=nombre de sommets parcourus dans G
    ... //
    free_param_bfs(p);

  Pour faire des appels multiples Ã  bfs() et optimiser le temps total:

    param_bfs *X=new_param_bfs();
    X->clean=1;  // initialisation complÃ¨te puis partielle de X->D
    ...
    X->vmax=100; // pour parcourir au plus 100 sommets
    bfs(G,u,X);  // initialise complÃ¨tement X->D Ã  -1 avant le bfs
    ...
    X->hmax=3;   // pour faire un bfs Ã  distance au plus 3 de v
    bfs(G,v,X);  // X->D sera initialisÃ© en initialisant seulement les sommets
    ...          // parcourus au bfs() prÃ©cÃ©dant (donc au plus 100 sommets)
    ...
    bfs(G,w,X);  // initialisation partielle de X->D
    ...
    free_param_bfs(X);

  La complexitÃ© est proportionnel au nombre d'arÃªtes dans la boule des
  sommets parcourus Ã  condition que X->D ne soit pas NULL, puisque
  sinon une initialisation de X->D en O(n) sera effectuÃ©e. C'est un
  point important si on lance plusieurs BFS partiels Ã  partir de
  sommets d'un mÃªme graphe. Pour Ãªtre efficace il faut, Ã  partir du 2e
  bfs, rÃ©tablir X->D[u]=-1 pour tous les sommets u de la X->file. On
  peut le rÃ©aliser en mettant X->clean=1. La complexitÃ© pour chaque
  appel (Ã  part le premier qui initialise complÃ¨tement X->D Ã  -1)
  reste en le nombre d'arÃªtes de la boule des sommets parcourus.

  ALGORITHME:
  - prendre le sommet u en tÃªte de file
  - pour tout ses voisins v non marquÃ©s:
    - enfiler v
    - marquer v
  
  Si D<>NULL, alors D n'est pas initialisÃ© (sauf si clean>0). Dans ce
  cas, il faut que D[u]=-1 pour tout sommet u, sauf pour les sommets
  que l'on souhaite ne pas visiter oÃ¹ il faut D[u]<>-1 (typiquement,
  D[u]=-2). On peut initialiser partiellement ou complÃ¨tement D avec
  clean=1 ou clean=2.

  Pour dÃ©terminer X->cycle, on calcule la longueur L du tout premier
  cycle crÃ©e. On remarque que X->cycle peut Ãªtre L ou L-1.  Il est L-1
  si plus tard on rencontre deux sommets voisins sur ce mÃªme niveau.
*/
{
  int const n=G->n;
  if((source<0)||(source>=n)) Erreur(48); // au cas oÃ¹
  int i,u,v,d,ff,tf,h;

  if(X==NULL) X=new_param_bfs(); /* NB: X->n=X->clean=0, X->P=X->file=NULL */
  if(X->P==NULL) ALLOC(X->P,n); /* alloue tableau si P==NULL */
  if(X->file==NULL) ALLOC(X->file,n); /* alloue la file */
  
  if(X->D==NULL) ALLOCZ(X->D,n,-1); /* alloue et initialise D en O(n) */
  else{ /* initialisation partielle ou pas de D */
    if(X->cont<2){
      if(X->clean==1) for(u=0;u<n;u++) X->D[u]=-1; /* initialisation complÃ¨te, et les -2 ??? */
      if(X->clean==2) for(i=0;i<X->n;i++) X->D[X->file[i]]=-1; /* initialisation partielle */
    }
  }
  if(X->clean==1) X->clean=2; /* la prochaine fois, initialisation partielle de X->D */

  if(X->cont==2){
    tf=X->tf;
    ff=X->n;
  }else{
    tf=0; /* tf=tÃªte de file, pointe sur la tÃªte */
    ff=0; /* ff=fin de file, pointe sur prochain Ã©lÃ©ment libre */
    if(X->cont==1) X->cont=2; /* la prochaine fois on continue */
    X->root=source; /* la racine est la source */
    X->P[source]=-1; /* pas de pÃ¨re pour la source */
    X->D[source]=0; /* distance=0 pour la source */
    X->file[ff++]=source; /* enfile le 1er sommet (=source), mÃªme s'il est supprimÃ© */
    X->cycle=(G->sym)? n+1 : 0; /* X->cycle non dÃ©fini (=0) si G orientÃ© */
  }
  
  int const hmax=(X->hmax==-1)? n : X->hmax; /* si X->hmax non dÃ©fini */
  int const vmax=(X->vmax==-1)? n : X->vmax; /* si X->vmax non dÃ©fini */
  h=(G->sym)? 1+(X->cycle/2) : 0;  /* h=hauteur Ã  partir de laquelle le plus court
				        cycle ne peut plus apparaÃ®tre */

  while(tf<ff){
    u=X->file[tf]; /* dÃ©file la tÃªte de file */
    if(X->D[u]>=hmax) break; /* fin si on a atteint une hauteur >=
				hmax. Dans ce cas, tous ceux de
				hauteur <= hmax on Ã©tÃ© enfilÃ©s. */
    tf++;
    for(i=0,d=G->d[u];i<d;i++){ /* pour tout voisin v de u */
      v=G->L[u][i];
      if(X->D[v]==-1){ /* si v voisin non marquÃ©, si =-2 on saute le sommet */
	X->P[v]=u; /* le pÃ¨re de v est u */
	X->D[v]=X->D[u]+1; /* hauteur(v)=1+hauteur(pÃ¨re(v)) */
	X->file[ff++]=v; /* enfile v */
	if(ff>vmax){ tf=ff; i=d; } /* fin si parcouru vmax sommets */
      }else /* si v a dÃ©jÃ  Ã©tÃ© visitÃ© ou s'il doit ne pas Ãªtre visitÃ© (=-2) */
	if((X->D[u]<h)&&(v!=X->P[u])&&(X->D[v]!=-2)){ /* sinon X->cycle ne peut plus Ãªtre amÃ©liorÃ©e */
	  h=X->D[u]+1; /* pas au delÃ  de X->D[u] */
	  X->cycle=min(X->cycle,h+X->D[v]);
	}
    }
  }

  if(X->cycle>n) X->cycle=-1; /* si > n, alors pas trouvÃ© de cycle -> -1 */
  X->n=ff;  /* nb de sommets parcourus */
  X->tf=tf; /* nb de sommets parcourus dont les voisins ont Ã©tÃ© enfilÃ©s */
  X->radius=X->D[u]; /* hauteur du dernier sommet dÃ©filÃ© (le plus Ã©loignÃ©) */
  
  /* c'est une mauvaise idÃ©e de faire ici un REALLOC(X->file,ff)
     car lors d'appels suivants avec la mÃªme structure, X->file n'aura
     plus la bonne taille ! */

  return X;
}


/*
  Structure de donnÃ©es pour le DFS
*/
typedef struct{
  int nc; // nc=nombre de composantes connexes du graphe
  int na; // na=nombre de sommets d'articulation (ou cut-vertex)
  int *C; // C[u]=couleur de la composante de u, entier de [0,nc[ ou sommet supprimÃ©
  int *P; // P[u]=parent de u dans une forÃªt couvrante, P[racine]=-1
  int *R; // R[i]=i-Ã¨me sommet racine dans la forÃªt couvrante, i dans [0,nc[
  int *d; // d[u]=date de dÃ©but de visite du sommet u, entier de [0,n[ (permutation)
  int *A; // A[u]=vrai ssi u est un sommet d'articulation, u dans [0,n[
  int *H; // H[u]=profondeur de u dans l'arborescence, H[racine]=0
} param_dfs;


param_dfs *new_param_dfs(int const n)
/*
  Attention ! le tableau X->C n'est pas allouÃ© par cette fonction,
  mÃªme si n>0. Cela doit Ãªtre fait par l'utilisateur.
*/
{
  NALLOC(param_dfs,X,1);

  X->C=X->P=X->R=X->d=X->A=X->H=NULL;
  X->nc=X->na=0;

  if(n>0){
    // X->C est allouÃ© par l'utilisateur ou par dfs()
    ALLOC(X->P,n);
    ALLOC(X->R,n);
    ALLOC(X->A,n);
    ALLOC(X->d,n);
    ALLOC(X->H,n);
  }
  return X;
}


void free_param_dfs(param_dfs* const X)
{
  if(X==NULL) return;
  free(X->C);
  free(X->P);
  free(X->R);
  free(X->A);
  free(X->d);
  free(X->H);
  free(X);
}


param_dfs *dfs(const graph* G,int source,param_dfs *X)
/*
  Effectue un parcours en profondeur de toutes les composantes
  connexes du graphe G depuis le sommet source, avec la date de
  premiÃ¨re visite. Version non rÃ©cursive. On dÃ©termine Ã©galement tous
  les sommets d'articulations (voir la dÃ©finition de param_dfs pour
  lire le rÃ©sultat), ainsi que les composantes connexes. On l'utilise
  comme suit:

  param_dfs *p=dfs(G,s,NULL); // DFS dans G depuis s
  ...
  free_param_dfs(p);

  ou alors, pour un DFS dans G Ã©vitant les sommets de T:

  param_dfs *p=new_param_dfs(G->n);
  ALLOCZ(p->C,G->n,-1);
  for(i=0;i<G->n;i++) p->C[T[i]]=-2;
  dfs(G,s,p);
  ...
  free_param_dfs(p);

  Si p->C[u]=-2, alors le sommet u n'est pas parcouru (il est
  virtuellement supprimÃ© de G). Les autres sommets v (non supprimÃ©s)
  doivent avoir p->C[v]=-1. Si p->C==NULL (par dÃ©faut), alors ce
  tableau est allouÃ© et initialisÃ© Ã  -1. Il sera libÃ©rÃ© par
  free_param_dfs(p).

  
  Pour lister les sommets dans l'ordre du parcours dfs, il suffit de
  faire:

  param_dfs *p=dfs(G,s,NULL); // DFS dans G depuis s
  NALLOC(int,order,G->n); // order[t]=u ssi dfs[u]=t
  for(u=0;u<G->n;u++) order[p->d[u]]=u;
  free(order);

  for(t=0;t<G->n;t++){
     u=order[t];
     ...; // traitement des sommets u dans l'ordre dfs
  }
  free_param_dfs(p);
*/
{  
  if(G==NULL) return NULL;
  int const n=G->n;
  if((source<0)||(source>=n)) Erreur(48); /* au cas oÃ¹ */

  int u,i,d,v,k,t,r=0;
  int tete,nc,na,b;

  if(X==NULL){ r=1; X=new_param_dfs(n); }
  if(X->C==NULL) ALLOCZ(X->C,n,-1);
  for(i=0;i<n;i++) X->A[i]=0;

  nc=na=0;
  NALLOC(int,pile,n);  /* pile */
  NALLOC(int,next,n);  /* next[u]=prochain voisin de u Ã  visiter */
  NALLOC(int,level,n); /* level[u]=... */
  t=tete=-1;

  for(k=0;k<n;k++,source=(source+1)%n)
    /* on essaye tous les sommets Ã  partir de source */
    if(X->C[source]==-1){ /* si ==-2 ou >=0 alors on saute le sommet */
      pile[++tete]=source;
      next[source]=0; /* premier voisin Ã  visiter */
      X->P[source]=-1;
      X->H[source]=0;
      X->R[nc]=source;

      while(tete>=0){ /* tant que la pile n'est pas vide */
	u=pile[tete]; /* u=sommet courant */
	i=next[u]; /* indice du prochain voisin de u Ã  visiter */
	if(i==0){
	  X->C[u]=nc; /* couleur de la composante courante */
	  level[u]=X->d[u]=++t; /* date de dÃ©but de visite */
	}
	d=G->d[u]; /* degrÃ© de u */
	b=1; /* sentiennelle pour savoir si on a trouvÃ© un v */
	while(i<d){ /* on cherche le prochain voisin v de u non visitÃ© */
	  v=G->L[u][i++]; /* pour tous les voisins v de u */
	  if(X->C[v]==-1){ /* si v n'a jamais Ã©tÃ© visitÃ© */
	    if((u==source)&&(t>X->d[u])&&(!X->A[u])) na++,X->A[u]=1; /* u=cut-vertex */
	    X->P[v]=u; /* pÃ¨re(v)=u */
	    X->H[v]=X->H[u]+1; /* mise Ã  jour de la profondeur */
	    pile[++tete]=v; /* on empile v */
	    next[v]=b=0; /* le prochain voisin de v Ã  visiter */
	    next[u]=i; /* le prochain voisin de u Ã  visiter */
	    break;
	  } else /* v existe et a dÃ©jÃ  Ã©tÃ© visitÃ© */
	    if((X->C[v]>=0)&&(v!=X->P[u])) /* si (u,v) est un arc de retour */
	      level[u]=min(level[u],X->d[v]);
	}
	if(b){ --tete; /* il n'y a plus de voisin v: on dÃ©pile u pour toujours */
	  if((v=(X->P[u]))>=0){ /* si u n'est pas une racine, il a un pÃ¨re v */
	    level[v]=min(level[v],level[u]); /* met Ã  jour level du pÃ¨re de u */
	    if((v!=source)&&(level[u]>=X->d[v])&&(!X->A[v])) na++,X->A[v]=1; /* v=cut-vertex */
	  }
	}
      } /* fin du while(i<d) */

      nc++; /* ici nc=nb de composantes visitÃ©es */
    }

  X->nc=nc;
  X->na=na;
  free(pile);
  free(next);
  free(level);

  /* on rÃ©duit le tableau ->R au minimum que si dfs() l'a
     allouÃ©. Sinon, il ne faut pas toucher aux pointeurs de X */
  if(r) REALLOC(X->R,nc);

  return X;
}


typedef struct param_bellman {
  int n;        // nombre de sommets du graphe
  int source;   // source
  int *parent;  // tableau des parents: parent[u], -1 si pas de pÃ¨re
  double *dist; // tableau de distance: dist[u]

  // variables internes pour Bellman_Ford()
  int *file; // file des sommets Ã  traiter
  int *mark; // pour avoir si un sommet est dans la file
  int multiple; // pour appels multiples sur le mÃªme graphe (valeurs 0,1,2)
  double wmin; // pour appels multiples, poids minimum > 0
} param_bellman;


param_bellman *new_param_bellman(int const n)
{
  NALLOC(param_bellman,p,1);
  p->n=n;
  p->source=-1;
  p->parent=p->file=p->mark=NULL;
  p->dist=NULL; // type double* incompatible avec int*
  p->multiple=0; // valeurs par dÃ©faut
  if(n>0){
    ALLOC(p->parent,n);
    ALLOC(p->dist,n);
    ALLOC(p->file,n);
    ALLOC(p->mark,n);
  }
  return p;
}


void free_param_bellman(param_bellman* const p)
{
  if(p){
    free(p->parent);
    free(p->dist);
    free(p->file);
    free(p->mark);
    free(p);
  }
  return;
}


double MinWeight(const graph* G)
/*
  Calcule le plus petits poids du graphe strictement positifs.
  Renvoie 0 s'il n'existe pas. TODO: On devrait tenir compte des
  sommets Ã©teints?
*/
{
  double w,wmin=INFINITY;
  const int n=G->n;
  for(int u=0;u<n;u++){ // pour tous les sommets u faire
    int d=G->d[u]; // d=deg(u)
    for(int i=0;i<d;i++){ // pour les d voisins de u
      w=G->W[u][i];
      if(wmin>w && w>0) wmin=w;
    }
  }
  if(wmin==INFINITY) wmin=0;
  return wmin;
}
 
 
param_bellman *Bellman_Ford(const graph* const G,int const source,param_bellman *p)
/*
  Calcule les distances du sommet source Ã  tous les autres de G selon
  l'algorithme de Bellman-Ford. La complexitÃ© est O(nm) dans le pire
  des cas, mais cette variante, inspirÃ©e de Sedgewick & Wayne avec une
  file (voir [SW11b]) est O(n+m) en pratique. On peut mÃªme accÃ©lÃ©rer
  le calcul en se basant sur un appel prÃ©cÃ©dent sur le mÃªme graphe
  mais une racine diffÃ©rente Ã  condition que les distances soient
  symÃ©triques. Il faut alors mettre p->multiple=1 (qui aprÃ¨s le
  premier appel passe Ã  2). Ne tient pas compte des sommets Ã©teints.

  https://algs4.cs.princeton.edu/code/
  https://algs4.cs.princeton.edu/44sp/BellmanFordSP.java.html

  Renvoie le rÃ©sultat dans la variable p et est aussi allouÃ©e et
  renvoyÃ©e si p=NULL. Le rÃ©sultat est deux tableaux: p->parent et
  p->dist. Il permet d'avoir des poids nÃ©gatifs sur les arcs,
  cependant il ne doit pas avoir de cycle absorbant sinon le programme
  boucle car il n'y a pas de dÃ©tection.

  Principe
  --------

  Classiquement, l'algorithme Bellman-Ford est le suivant:

  1) Initialiser les tableaux dist[] et parent[]

  2) RÃ©pÃ©ter n-1 fois:
       Pour tout les arcs (u,v) de G:
          relax(u,v)

  oÃ¹ relax(u,v) correspond Ã :
   si dist[u]+w(u,v)<dist[v], alors:
   - dist[v]=dist[u]+w(u,v)
   - parent[v]=u

  La variante implÃ©mentÃ©e ici consiste Ã  traiter uniquement les
  sommets v qui ont menÃ© Ã  une relaxation (u,v), c'est-Ã -dire aux
  sommets dont la distance a diminuÃ©e, et ainsi d'Ã©viter de rÃ©pÃ©ter
  n-1 fois la relaxation de tous les arcs. Les sommets Ã  traiter sont
  gÃ©rÃ©s par une file. Les sommets dans la file sont marquÃ©s ce qui
  permet de n'ajouter que les sommets qui ne sont pas eux-mÃªmes Ã 
  traiter.

  Utilisation multiples
  ---------------------

  L'algorithme peut Ãªtre accÃ©lÃ©rÃ© si un appel sur le mÃªme graphe a
  dÃ©jÃ  eut lieu, depuis une autre source. L'idÃ©e est que si l'on
  dispose d'une "bonne" borne supÃ©rieure lors de l'initialisation de
  dist[], alors le nombre de sommets Ã  traiter sera rÃ©duit. Par dÃ©faut
  p->multiple=0 ce qui conduit au comportement normal (sans appel
  multiple) avec une initialisation complÃ¨te des tableaux p->dist[],
  p->parent[] et p->mark[]. Si p->multiple=1, alors on effectue aussi
  l'initialisation normale mais on passe ensuite Ã  p->multiple=2.  Si
  p->multiple=2, alors p->dist[] et p->parent[] sont alors mises Ã 
  jour. Il faut cependant faire trÃ¨s attention de ne pas mettre une
  mise Ã  jour trop basse des distances. Par exemple si la nouvelle
  source est de degrÃ© 1, la nouvelle distance Ã  son voisin sera, aprÃ¨s
  la mise jour, exacte. La boucle while s'arÃªtera aprÃ¨s avoir visitÃ©
  le voisin de la source, sans propager les relax() et traiter ce
  voisin comme il se doit, ce qui est incorrect. De maniÃ¨re gÃ©nÃ©rale
  le problÃ¨me survient si la nouvelle source est sÃ©parÃ©e d'une partie
  du graphe par un ensemble de sommets dont les distances mises Ã  jour
  se trouvent Ãªtre exactes. La propagation des relax() ne se fera pas
  dans une partie du graphe. Pour y remÃ©dier, il faut lÃ©gÃ¨rement
  surestimer ces distances. On ajoute alors aux distances la valeur du
  plus petits poids strictement positifs et 0 si tous les poids
  Ã©taient nuls ou nÃ©gatifs.

  Le surcoÃ»t en temps est O(n) pour la mise Ã  jour de la
  distance. Pour la mise Ã  jour de la forÃªt cela coÃ»te le nombre
  d'arcs du plus court chemin entre les deux sources. Le gain est donc
  apprÃ©ciable pour de multiples appels Ã  Bellman_Ford(), comme pour le
  calcul du stretch de graphes gÃ©omÃ©triques par exemple.

  Exemple d'utilisations multiples:

    param_bellman *p=new_param_bellman(n); // par dÃ©faut p->multiple=0
    p->multiple=1; // pour appels multiples
	
    for(s=...){ // pour diffÃ©rentes sources s
      Bellman_Ford(G,s,p);
      // traitement du rÃ©sultat p
      ...;
      }
    free_param_bellman(p);

  Plus prÃ©cisÃ©ment, lors de l'initialisation, on rÃ©oriente les arcs
  allant de source Ã  p->source, ce qui Ã  pour effet et former une
  forÃªt dont la nouvelle racine est source (et non plus p->source).
  Les distances sont mises Ã  jour en tenant compte de la distance
  p->dist[source]. Pour Ãªtre correct il faut que le chemin de
  p->source Ã  source soit de mÃªme coÃ»t que le chemin inverse de source
  Ã  q->source, ce qui est le cas si le graphe est symÃ©trique. Rien
  n'est recalculÃ© si source=p->source.

  Pour la mise Ã  jour des distances, on commence par ajouter Ã  tous
  les sommets w0+wmin, soit la distance w0=p->dist[source] entre
  l'ancienne et la nouvelle source qu'on surestime avec le poids
  minimum wmin qui est >0. En effet, la distance Ã  v par rapport Ã  la
  nouvelle source est au plus la distance par rapport Ã  l'ancienne
  distance + la distance entre les sources (w0).

  La surestimation par wmin fait qu'aucune distance n'est exacte ce
  qui va bien propager les relax() Ã  tout le graphe. Pour les
  distances des sommets v du chemin de la nouvelle source Ã 
  l'ancienne, on fait comme ceci: soit dv la distance d'origine entre
  ancienne source et v, c'est p->dist[v] avant toute modification. La
  distance exacte depuis la nouvelle source est: w0 - dv. Il faut donc
  mettre Ã  jour avec: p->dist[v] = w0 - dv + wmin pour avoir une
  lÃ©gÃ¨re surestimation de la distance. Comme entre-temps on a posÃ©
  p->dist[v] += (w0+wmin), on a dv = p->dist[v] - (w0+wmin). Au final,
  il faut donc poser:

     p->dist[v] = w0 - dv + wmin
                = w0 - (p->dist[v]-(w0+wmin)) + wmin
		= 2(w0+wmin) - p->dist[v].

  NB: Ã€ la fin de la boucle while, tous les sommets de la composante
  connexe ont Ã©tÃ© traitÃ© et donc plus aucun sommet n'est dans la
  file. Donc mark[u]=1 pour tous les u, ce qui correspond Ã 
  l'initialisation normale. On a donc pas besoin de s'en occuper de
  son initialisation pour les appels multiples.
*/
{
  if(G==NULL) return NULL;
  if(G->W==NULL) return NULL;
  int const n=G->n;
  if((source<0)||(source>=n)) return NULL;
  if(p==NULL) p=new_param_bellman(n);

  int u,v,i,front,back,d;
  double w;

  /* initialisations */

  if(p->multiple<2){ // initialisations standard
    for(u=0;u<n;u++){
      p->dist[u]=INFINITY; /* +âˆž */
      p->parent[u]=-1; /* pÃ¨res non dÃ©finis, y compris pour la source */
      p->mark[u]=1; /* mark[u]=0 ssi u est dans la file */
    }
    if(p->multiple==1){ // si appels multiples
      p->multiple=2; // pour les appels suivants
      p->wmin=MinWeight(G); // calcule wmin, +petit poids>0 s'il existe
    }
  }else{ // initialisations pour appels multiples
    if(source==p->source) return p; // ne devrait pas arriver
    /* mise Ã  jour des distances en tenant compte de la nouvelle racine */
    w=p->dist[source]+p->wmin; // distance entre ancienne et nouvelle source + wmin
    for(u=0;u<n;u++) p->dist[u] += w; // borne sup sur la distance
    /* changement de racine de l'arbre: p->source vers source */
    u=source; // source<>p->source
    v=p->parent[u]; // donc ici v<>-1
    w *= 2; // voir explications ci-dessus
    while(v>=0){
      i=p->parent[v]; // sauve parent[v], existe si v>=0
      p->parent[v]=u; // change parent[v]
      p->dist[v]=w-p->dist[v]; // distance de v Ã  la nouvelle source + wmin
      u=v;
      v=i;
    }
  }
  
  p->source=source; /* mÃ©morise la source */
  p->dist[source]=0; /* distance nulle pour la source */
  p->parent[source]=-1; /* la source n'a pas de parent */
  if(n==1) return p; /* fin dans ce cas, p->mark est correct */

  front=back=0; /* indices pour la file: front=lecture, back=Ã©criture */
  p->file[back++]=source; /* empile la source dans la file. NB: n>1=back */
  p->mark[source]=0; /* source est dans la file */

  while(front!=back){ /* tant qu'il y a des sommets Ã  traiter */
    u=p->file[front++]; if(front==n) front=0; /* dÃ©file u */
    p->mark[u]=1; /* u n'est plus dans la file */
    d=G->d[u]; /* d=degrÃ© de u */
    for(i=0;i<d;i++){ /* pour tous les voisins v de u */
      v=G->L[u][i]; /* ici: si v est Ã©teint, il faut le sauter */
      w=p->dist[u]+G->W[u][i];
      if(w<p->dist[v]){ /* relax(u,v) */
	p->parent[v]=u;
	p->dist[v]=w;
	if(p->mark[v]){ /* enfiler v s'il n'y est pas */
	  p->file[back++]=v; if(back==n) back=0;
	  p->mark[v]=0; /* v est dans la file */
	}
      }
    }
  }

  return p;
}


/*
  A faire: algorithme A* (incluant Dijkstra) et l'implÃ©mentation du
  tas paresseux avec un heap binaire

  Dijkstra(G,W,s):

  1. Init(G,s)
  2. S={}
  3. Q=V(G)
  4. while Q<>{}
  5.  do u=Extract-Min(Q)
  6.     S=S u {u}
  7.     for each v in Adj[u]: Relax(u,v,W)

  Init(G,s):
  1. for each vertex v in V(G):
  2.   do d[v]=+âˆž
  3.      pi[v]=NULL
  4. d[s]=0

  Relax(u,v,W):
  1. if d[v]>d[u]+W[u,v]
  2.   then d[v]=d[u]+W[u,v]
  3.        pi[v]=u

  Extract-Min(Q):
  extract vertex u with smallest d[v]
*/


//
// Pour la gÃ©nÃ©ration des poids des arÃªtes. On pourrait faire un
// RANDOM dans [a,b]. Si le graphe est gÃ©omÃ©trique, on pourrait de
// maniÃ¨re gÃ©nÃ©rale prendre toutes les faÃ§ons de gÃ©nÃ©rer des points en
// 2D utilisÃ©es par -xy et prendre la norme entre les extrÃ©mitÃ©s des
// arÃªtes.
//
enum{
  IW_UNW,  // graphe valuÃ© uniformÃ©ment
  IW_GEO,  // graphe euclidien
  IW_POS,  // IW_UNW ou IW_GEO suivant la valeur de POS
  IW_R01,  // valuation alÃ©atoires dans [0,1]
};


int InitWeights(graph* const G,int code)
/*
  Initialise le poids des arcs du graphe G, c'est-Ã -dire le tableau
  G->W, qui est Ã©ventuellement allouÃ©. L'initialisation dÃ©pend de code
  (cf. l'enum IW_xxx). Retourne 1 ssi tout c'est bien passÃ©, 0 sinon.
*/
{
  if(G==NULL) return 0;
  if(code==IW_POS) code=(POS)? IW_GEO : IW_UNW;
  if((code==IW_GEO)&&(G->xpos==NULL)) return 0;

  int const n=G->n;
  int u,v,i,d;
  if(G->W==NULL) ALLOCZ(G->W,n,NULL);
  
  // parcourt tous les arcs de G
  for(u=0;u<n;u++){
    d=G->d[u];
    free(G->W[u]);
    ALLOC(G->W[u],d);
    for(i=0;i<d;i++){
      v=G->L[u][i];
      switch(code){
      case IW_GEO: G->W[u][i]=hypot(G->xpos[u]-G->xpos[v],G->ypos[u]-G->ypos[v]); break;
      case IW_UNW: G->W[u][i]=1; break;
      case IW_R01: G->W[u][i]=RAND01; break;
      default: return 0; // problÃ¨me
      }
    }
  }
  
  return 1;
}

/***********************************

       ISOMORPHISM, SUBGRAPH,
       MINOR, INDUCEDSUBGRAPH,
       PATHS, PS1, TREEWIDTH ...

***********************************/


int *Isomorphism(graph* const G,graph* const H)
/*
  Renvoie un tableau P <> NULL ssi G est isomorphe Ã  H. Si tel est le
  cas, le tableau P indique le morphisme de H vers G. AprÃ¨s l'appel,
  le graphe H est modifiÃ©: ses listes sont triÃ©es si H->sort=0 est
  faux (sauf si G=H - mÃªme pointeur). Le graphe G n'est par contre pas
  modifiÃ©. Dans H->int1 est retournÃ© le nombre de tests effectuÃ©s.
  Moins le graphe possÃ¨de de symÃ©trie, plus faible est le nombre de
  tests (et rapide est la fonction).

  On applique l'algorithme suivant. Pour chacun des deux graphes et
  chaque sommet u, on calcule son "profile" un tableau notÃ©
  profile[u]: profile[u][i+2] = le nombre de sommets Ã  distance i de
  u, en commenÃ§ant Ã  partir de i=1. Donc, profile[u][3] est le degrÃ©
  de u. Ceci est calculÃ© par un simple BFS, les indices 0,1,2 Ã©tant
  rÃ©servÃ©s. On utilise profile[u][0] pour coder la taille du tableau
  profile[u], profile[u][1] pour stocker le nom de u (plutÃ´t que le
  nombre de sommet Ã  distance 0 qui est toujours 1) et profile[u][2]
  pour stocker la longueur du plus petit cycle passant par u. Cette
  derniÃ¨re information Ã©tant calculÃ©e "gratuitement" par le BFS.

  On ordonne ensuite les sommets des graphes suivant les profiles des
  sommets avec qsort(). On ne renumÃ©rote pas les sommets dans le
  graphe, mais plutÃ´t on donne un ordre: c'est possible avec qsort()
  car le nom original u est dans profile[u][1] le sommet u. Si bien
  que profile[i][1] sera le sommet u de rang i. On priviligie les
  profiles le grande taille (que l'on classe en premier) ce qui est
  plus discriminant. Les isomorphismes (=permutations) Ã  tester ne
  concernent que les sommets de mÃªme profile. On construit les
  contraintes dans un tableau C, C[j] indiquant que les sommets de
  rang C[j-1] Ã  C[j] (exclu) ont mÃªme profile, et sur lequel se base
  NextPermutation().

  Sur le mÃªme principe, on pourrait imaginer un profile plus complexe,
  comme suit: Ã  chaque distance i et sommet u, on calcule le graphe
  G[u][i] induit par les sommets Ã  distance i de u. On peut alors
  calculer le profile de chaque sommet de G[u][i] et ordonner les
  sommets selon celui-ci.
*/
{
  if(H) H->int1=0; /* par dÃ©faut, 0 tests */
  if((G==NULL)||(H==NULL)) return NULL;
  if(G->n!=H->n) return NULL;

  int *P; /* isomorphisme final */
  int const n=G->n;

  if(G==H){ /* isomorphisme trivial si mÃªme emplacement mÃ©moire */
    ALLOCZ(P,n,_i);
    return P;
  }

  param_bfs *param=new_param_bfs(); /* pour le BFS */
  int **profile,**profileG=NULL,**profileH;
  int *R,*C; /* permutation et contraintes (sur les rangs) */
  int u,v,t,r,i;
  graph* M;

  for(M=G;;){ /* on fait la mÃªme chose pour M=G puis M=H */
    ALLOC(profile,n); /* profile[u] */
    for(u=0;u<n;u++){ /* faire un BFS pour tous les sommets u */
      bfs(M,u,param); /* le premier BFS va allouer param->D et param->P */
      t=3+param->radius; /* taille du tableau profile[u] */
      ALLOCZ(profile[u],t,0); /* initialise le profile */
      for(v=0;v<n;v++){
	i=param->D[v];
	if(i>0) profile[u][i+2]++; /* compte les sommets Ã  distance i>0 de v */
	param->D[v]=-1; /* rÃ©initialise les distances pour les BFS suivants */
      }
      profile[u][0]=t; /* taille du profile */
      profile[u][1]=u; /* nom du sommet, pour qsort() */
      profile[u][2]=param->cycle; /* maille */
    }
    QSORT(profile,n,fcmp_profile); /* trie les profiles */

    if(M==H){ profileH=profile; break; } /* on s'arÃªte si M=H */
    profileG=profile; /* on refait la boucle pour H */
    M=H;
  }
  free_param_bfs(param);

  /* on verifie que profileG "=" profileH */
  for(u=0;u<n;u++)
    if(fcmp_profile(profileG+u,profileH+u)){
      P=NULL;
      goto fin_noniso;
    }

  /* calcule les contraintes */
  /* ici les profiles de G et H sont identiques */
  ALLOC(C,n);
  R=profile[0]; /* R=profile du premier sommet. NB: profile=profileH. */
  for(u=t=0;u<n;u++){
    if(fcmp_profile(&R,profile+u)){ /* si profiles diffÃ©rent */
      R=profile[u];
      C[t++]=u;
    }
  }
  C[t]=n;

  ALLOC(P,n);
  ALLOCZ(R,n,_i); /* initialise l'isomorphisme sur les rangs des profiles */
  if(!H->sort) SortGraph(H,0); /* on trie H pour la recherche dichotomique */
  H->int1=0; /* compte le nb de tests */

  /* vÃ©rifie, pour chaque permutation P, que P(G)=H */

  do{
    H->int1++;

    /* calcule P en fonction de R */
    for(r=0;r<n;r++) P[profileG[r][1]]=profileH[R[r]][1];
    for(r=0;r<n;r++){ /* on commence par les profiles longs: r=0,1,2, ... */
      u=profileG[r][1]; /* v=P[u] est le sommet correspondant Ã  u dans H */
      for(i=0,t=G->d[u];i<t;i++) /* on regarde si chaque voisin de u est dans H->L[v] */
	if(bsearch(P+(G->L[u][i]),H->L[P[u]],t,sizeof(int),fcmp_int)==NULL){
	  /* alors Ã©lÃ©ment non trouvÃ© */
	  i=t;r=n; /* prochaine permutation Ã  tester */
	}
    } /* ici r=n (trouvÃ©) ou n+1 (pas encore trouvÃ©) */
    if(r==n) goto fin_iso; /* on a trouvÃ© un isomorphisme P */
  }
  while(NextPermutation(R,n,C));

  /* si on arrive ici, c'est qu'on a pas trouvÃ© l'isomorphisme */
  free(P);
  P=NULL;

 fin_iso:
  free(R);
  free(C);

 fin_noniso:
  FREE2(profileG,n);
  FREE2(profileH,n);

  return P;
}


edge *ListEdges(graph* const G)
/*
  Construit la liste des arÃªtes de G, chaque arÃªte uv ne figure qu'une
  seule fois. On a, pour tout i, E[i].u < E[i].v. Le champs G->m est
  mis Ã  jour.
*/
{
  int const m=nb_edges(G);
  int const n=G->n;
  int u,v,i,j,d;

  NALLOC(edge,E,m);
  for(u=j=0;u<n;u++){
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v) E[j].u=u, E[j].v=v, j++;
    }
  }
  return E;
}


graph* Subgraph(graph* const G,graph* const H)
/*
  DÃ©termine si G est un sous-graphe de H dans le cas oÃ¹ ils ont mÃªme
  nombre de sommets. On renvoie un sous-graphe S de H isomorphe Ã  G,
  et on renvoie NULL si H ne possÃ¨de pas G comme sous-graphe ou si G
  et H n'ont pas le mÃªme nombre de sommets.

  Effets de bord:
  - les listes de G sont triÃ©es et donc G->sort=1,
  - H->int1 contient le nombre total de tests effectuÃ©s,
  - S->pint1 contient l'isomorphisme de S vers G, et donc de H vers G.

  L'algorithme est le suivant: on teste d'abord si la sÃ©quence des
  degrÃ©s de H est bien supÃ©rieure Ã  celle de G (ceci prend un temps
  O(n)). Si c'est le cas, on effectue, pour tous les sous-graphes S de
  H qui ont autant d'arÃªtes que G, un test d'isomorphisme entre S et G
  grÃ¢ce Ã  isomorphisme(S,G).
*/
{
  int const n=H->n;
  H->int1=0; /* nb de tests = 0 par dÃ©faut */
  if(n!=G->n) return NULL; /* pas le mÃªme nombre de sommets */

  /* on trie en O(n) les deux listes */
  int* const Eh=SortInt(H->d,NULL,n,0,NULL,SORT_INC);
  int* const Eg=SortInt(G->d,NULL,n,0,NULL,SORT_INC);
  int i;

  /* on s'arrÃªte si, pour un rang i donnÃ©, degH(i)<degG(i) */
  for(i=0;i<n;i++) if(Eh[i]<Eg[i]) break;
  free(Eh);
  free(Eg);
  if(i<n) return NULL; /* G ne peut pas Ãªtre un sous-graphe de H */

  int mg=nb_edges(G);
  int mh=nb_edges(H);
  graph* S=new_subgraph(H); /* S=sous-graphe de H, alloue S->L et dimensionne S->d */
  edge* const E=ListEdges(H); /* liste des arÃªtes de H: e_j=(u,v) -> E[j].u=u et E[j].v=v */
  int *B; /* les arÃªtes de S, sous-ensemble d'indices d'arÃªtes de H */
  int *P; /* isomorphisme S->G */
  int u,v,j,d;
  
  /* initialise le sous-ensemble d'arÃªtes B de H avec mg arÃªtes */
  ALLOC(B,mg);
  NextSet(B,-1,mg);
  d=0; /* d=compteur pour le nombre total de tests */
  
  do{
    
    /* remplit S avec les arÃªtes de B */
    degres_zero(S); /* position libre pour le sommet u de S */
    for(i=0;i<mg;i++){ j=B[i]; /* j=numÃ©ro de la i-Ã¨me arÃªte de B */
      u=E[j].u; v=E[j].v; /* l'arÃªte j de E est (u,v) */
      ADD_EDGE(S,u,v); /* ajoute l'arÃªte u,v */
    }
    
    /* il vaut mieux que G soit le 2e paramÃ¨tre, car il va Ãªtre triÃ©
       la premiÃ¨re fois par Isomorphism(), puis plus jamais grÃ¢ce au
       test de G->sort, alors que S serait triÃ© Ã  chaque appel */
    P=Isomorphism(S,G);
    d += 1+G->int1; /* on ajoute 1 pour donner le nombre d'ensembles testÃ©s */
  }while((P==NULL)&&(NextSet(B,mh,mg)));
  
  H->int1=d; /* nombre total de tests */

  if(P==NULL){ /* on a pas trouvÃ© de sous-graphe de H isomorphe Ã  G */
    free_graph(S);
    S=NULL;
  }
  else S->pint1=P; /* S isomorphe Ã  G, sous-graphe de H */

  free(E);
  free(B);
  return S;
}


graph* MatrixToGraph(int **M,int const n){
/*
  Renvoie le graphe correspondant Ã  la matrice d'adjacence n x n oÃ¹
  seule la partie infÃ©rieure est utilisÃ©e. Les listes du graphe de
  retour sont triÃ©es et les champs ->m et ->sort sont mise Ã  jour.
*/
  if((M==NULL)||(n<=0)) return NULL;
  int u,v,m;
  graph* G=new_fullgraph(n);

  for(u=m=0;u<n;u++)
    for(v=0;v<u;v++)
      if(M[u][v]){
	ADD_EDGE(G,u,v); /* ajoute u-v et v-u */
	m++; /* une arÃªte de plus */
      }
  
  /* rÃ©duit les listes */
  GraphRealloc(G,G->d);

  G->m=m;
  G->sort=1;
  return G;
}


graph* GraphOfColor(const graph* const G,int* const col,int const k){
/*
  Renvoie un graphe C dont les sommets sont les entiers de [0,k[ (les
  valeurs du tableau col[] qui doit avoir une taille G->n) et les
  arÃªtes les paires uv telle qu'il existe une arÃªte xy de G avec
  col[x]=u et col[y]=v. La valeur C->m est mise Ã  jour, et les listes
  de C sont triÃ©es (C->sort=1).
*/
  if((k<0)||(col==NULL)||(G==NULL)) return NULL;

  int const n=G->n;
  int u,v,cu,cv,i,d;
  graph* C; /* le graphe des couleurs renvoyÃ© */

  /* matrice d'adjacence infÃ©rieure Ã  0 */
  NALLOC2(int,M,k,k-1); /* matrice d'adjacence du graphe des couleurs */
  for(u=0;u<k;u++)
    for(v=0;v<u;M[u][v++]=0);

  for(u=0;u<n;u++) /* parcourt G et remplit la partie infÃ©rieure de M */
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v){ /* si cu=cv, on ne fait rien */
	cu=col[u];
	cv=col[v];
	if(cu>cv) M[cu][cv]=1;
	if(cv>cu) M[cv][cu]=1;
      }
    }
  
  C=MatrixToGraph(M,k);
  free(M);
  return C;
}


int UF_Find(int const x,int* const parent)
/*
  Routine FIND pour le problÃ¨me d'UNION-FIND. Donne le reprÃ©sentant de
  x dans le tableau parent[] ou encore la racine de x dans la forÃªt
  couvrante donnÃ©e par la relation parent. AprÃ¨s un appel Ã  UF_Find(),
  si r=UF_Find(x,parent), alors on a forcÃ©ment que parent[x]=r.

  Le temps total de m requÃªtes Ã  UF_Find(), pour un tableau parent[]
  de taille n, est O(m*ð›¼(n)) lorsque qu'on utilise les heuristiques du
  rang (rÃ©alisÃ©e par UNION) et de la compression de chemin (rÃ©alisÃ©e
  par FIND). Ici ð›¼(n) est la fonction inverse d'Ackerman, et ð›¼(n)<=4
  en pratique.

  Ex: On utilise typiquement UF_Find() pour savoir si, lorsqu'on
  ajoute une arÃªte Ã  un graphe, on crÃ©e un cycle ou pas. Si on ne crÃ©e
  par de cycle, il faut appeler UF_Union(x,y,parent,height).

  NALLOCZ(int,parent,n,_i);
  NALLOCZ(int,height,n,0);

  for(u=...)
    for(v=...){ // pour chaque arÃªte {u,v}
      x=UF_Find(u,parent); // x=reprÃ©sentant de u
      y=UF_Find(v,parent); // y=reprÃ©sentant de v
      if(x==y){ // il y a un cycle
      ...;
      }else{    // il n'y a pas de cycle
      ...;
      UF_Union(x,y,parent,height); // fusionne la composante de u et de v
      }
   }
*/
{
  if(x!=parent[x]) parent[x]=UF_Find(parent[x],parent);
  return parent[x];
}


static inline void UF_Union(int const x,int const y,
			    int* const parent,int* const height)
/*
  Routine UNION pour UNION-FIND, avec l'heuristique du "rang".
*/
{
  if(height[x]>height[y]) parent[y]=x;
  else{
    parent[x]=y;
    if(height[x]==height[y]) height[y]++;
  }
}


int *Minor(graph* const H,graph* const G)
/*
  DÃ©termine si H est un mineur de G. Si c'est le cas, un tableau T est
  renvoyÃ©, sinon NULL est renvoyÃ©. Le tableau T code un modÃ¨le du
  mineur de H dans G. Plus prÃ©cisÃ©ment, T[u]=c, oÃ¹ u est un sommet de
  G, est le nom du sommet c de H si bien que l'ensemble des sommets u
  de G tel que T[u]=c forme le super-noeud c.

  L'algorithme est le suivant: on effectue, pour tous les ensembles de
  contractions d'arÃªtes de G produisant un mineur avec autant de
  sommets que H, un test de sous-graphe grÃ¢ce Ã  Subgraph().
*/
{
  graph* C; /* graphe des couleurs = graphe contractÃ© */
  graph* S=NULL; /* sous-graphe Ã©ventuellement isomorphe Ã  C */
  int *B; /* sous-ensemble (d'indices) d'arÃªtes de G Ã  contracter */
  int *couleur; /* les sommets de mÃªme couleurs seront contractÃ©s */
  int *rang; /* pour union-find rapide */
  edge e;
  
  int nh=H->n;
  int ng=G->n;
  int c=ng-nh; /* c=nb de contractions Ã  effectuer */
  int t;
  H->int1=t=0; /* initialise le nb de tests */
  if(c<0) return NULL; /* pas de mineur */

  edge *E=ListEdges(G); /* E=liste des arÃªtes de G, met Ã  jour G->m */
  int mg=G->m;
  if(c>mg){ /* fini s'il faut contracter plus d'arÃªtes qu'il y a dans G */
    free(E);
    return NULL;
  }

  int i,u,v,x,y;
  int test=((c<<1)<ng); /* vrai ssi on liste suivant les arÃªtes ou suivant les sommets */
  ALLOC(B,c); NextSet(B,-1,c); /* B=premier sous-ensemble de c arÃªtes de E */
  ALLOC(couleur,ng); /* couleur des sommets */
  ALLOC(rang,ng); /* rang des sommets, sert pour UNION-FIND */

  /*
    On pourrait gÃ©nÃ©rer des sous-ensembles acycliques d'arÃªtes
    directement, en combinant NextSet() et le test d'acyclicitÃ©.
   */

  do{

    t++; /* on compte 1 test par ensemble B testÃ©s */

    /* initialise la couleur et le rang des sommets */
    for(u=0;u<ng;u++){ couleur[u]=u; rang[u]=0; }

    /* on teste rapidement (avec UNION-FIND) si B contient un cycle */
    for(i=0;i<c;i++){ e=E[B[i]]; /* e=i-Ã¨me arÃªte de B */
      u=e.u; x=UF_Find(u,couleur);
      v=e.v; y=UF_Find(v,couleur);
      if(x==y) break; /* y'a un cycle, sinon on fait UNION */
      UF_Union(x,y,couleur,rang);
    }

    if(i==c){ /* si B est acyclique, on fait un test de sous-graphe */
      if(test)
	/* on met Ã  jour la couleur de chaque sommet. Suivant les
	   valeurs respectives de c et ng (test) on met Ã  jour soit
	   suivant les arÃªtes ou suivant les sommets. */
	for(i=0;i<c;i++){
	  e=E[B[i]]; /* e=i-Ã¨me arÃªte de B */
	  UF_Find(e.u,couleur); // couleur[u]=UF_Find(u,...)
	  UF_Find(e.v,couleur); // couleur[v]=UF_Find(v,...)
	}
      else
	for(u=0;u<ng;u++) UF_Find(u,couleur);

      /* on recadre les couleurs dans [0,c[. ComplexitÃ©: 4ng */
      for(i=0;i<ng;rang[i++]=0);
      for(i=0;i<ng;i++) rang[couleur[i]]++; /* rang=frÃ©quence */
      for(i=u=0;i<ng;i++) /* repÃ¨re les couleurs manquantes */ 
	if(rang[i]==0) u++; else rang[i]=u;
      /* ici rang[i]=nb de zÃ©ros (=valeurs manquantes) dans rang[0..i[ */
      for(i=0;i<ng;i++) couleur[i] -= rang[couleur[i]];

      C=GraphOfColor(G,couleur,nh);
      S=Subgraph(H,C); /* avant l'appel, S=NULL nÃ©cessairement */
      t += C->int1;
      free_graph(C); /* on a plus besoin du graphe des couleurs */
    }

  }while((S==NULL)&&(NextSet(B,mg,c)));

  H->int1=t;
  free(B);
  free(E);

  /* on a rien trouvÃ© */
  if(S==NULL){
    free(couleur);
    free(rang);
    return NULL;
  }
  
  /* on a trouvÃ© un mineur, on construit le modÃ¨le dans rang[] */
  for(u=0;u<ng;u++) rang[u]=S->pint1[couleur[u]];
  free_graph(S);
  free(couleur);
  return rang;
}


int *InducedSubgraph(graph* const H,graph* const G)
/*
  Indique si H est un sous-graphe induit de G.  La fonction renvoie un
  ensemble X de sommets tel que G[X] est ismomorphe Ã  H. Ã‰videmment
  |X|=H->n. On renvoie dans G->int1 le nombre de tests effectuÃ©s, et
  dans G->pint1 l'isomorphisme entre H et G[X].
  
  L'algorithme consiste Ã  gÃ©nÃ©rer tous les ensembles X possibles de
  |V(H)| sommets et Ã  tester l'isomorphisme entre G[X] et H.
 */
{
  if((G==NULL)||(H==NULL)) return NULL;
  int ng=G->n,nh=H->n;
  if(nh>ng) return NULL;

  graph* S;
  int *P,t=0;
  NALLOC(int,X,nh);
  NextSet(X,-1,nh); /* premier sous-ensemble */

  do{
    t++;
    S=ExtractSubgraph(G,X,nh,1);
    P=Isomorphism(S,H);
    t += H->int1;
    free_graph(S);
  }while((P==NULL)&&(NextSet(X,ng,nh)));

  G->int1=t;
  free(G->pint1); /* pour Ã©viter les fuites mÃ©moires */
  G->pint1=P;
  if(P==NULL){ free(X); X=NULL; }
  return X;
}


int NextPath(graph* const G,path* const P,int j)
/*
  Cette fonction (rÃ©cursive) permet de trouver tous les chemins
  simples entre deux sommets. Plus prÃ©cisÃ©ment, on met Ã  jour le
  chemin P de sorte qu'on trouve un nouveau chemin allant du j-Ã¨me
  sommet de P au dernier tout en Ã©vitant la premiÃ¨re partie du chemin
  allant du premier sommet de P au j-Ã¨me. Si un tel chemin a pu Ãªtre
  trouvÃ© on renvoie 1, sinon on renvoie 0. On met j=-1 s'il s'agit de
  la crÃ©ation d'un chemin de P->P[0] Ã  P->P[1].

  On l'utilise comme ceci:

    path *P=new_path(G,NULL,G->n); // crÃ©e un chemin vide P mais avec G->n sommets possibles
    P->P[0]=x; // origine du chemin
    P->P[1]=y; // destination du chemin, y=x possible
    if(NextPath(G,P,-1)==0){ ... } // crÃ©e le premier chemin, =0 s'il n'y en a pas
    do{ // traitement du chemin P
      ...
    }while(NextPath(G,P,0); // tantqu'il existe un nouveau chemin P
    free_path(P); // libÃ¨re le chemin P

  Plus prÃ©cisÃ©ment, Ã©tant donnÃ©s un chemin P=x-...-v_j-...-y du graphe
  G et un sommet v_j du chemin (v_j=j-Ã¨me sommet du chemin), la
  fonction tente de complÃ©ter P par un autre chemin de v_j Ã  y Ã©vitant
  x-...-v_(j-1). Si ce chemin Ã  Ã©tÃ© trouvÃ©, alors la partie de v_j Ã  y
  de P est mise Ã  jour et on renvoie 1. Sinon, le chemin est coupÃ©
  aprÃ¨s v_j et on renvoie 0. Dans tous les cas P est un chemin Ã  jour
  de G. Si j<0, alors on initialise P par un chemin allant de
  x=P->P[0] Ã  y=P->P[1].

  Algorithme: on essaye d'amÃ©liorer en premier le sous-chemin de
  v_{j+1} Ã  y (rÃ©cursivement). Si cela n'est pas possible, on calcule
  un nouveau chemin de v_j Ã  y passant par un autre voisin v de v_j
  (autre que v_{j+1}) et Ã©vitant le chemin x-...-v_j. On passe en
  revue ainsi tous les voisins v de v_j qui ne sont pas dans
  P[0]...P[j]. Si aucun des voisins ne possÃ¨de un tel chemin, c'est
  qu'il y en a plus et on retourne 0.

  Comme il faut tester les voisins v de v_j qu'une seule fois (pour un
  chemin P[0]...P[j] donnÃ©), on utilise le champs aux[v_j][i] de P qui
  donne le i-Ã¨me voisin libre de v_j avec la convention que
  aux[v_j][0] est le nombre de voisins encore possibles.

  Effet de bord: P est mis Ã  jour.
*/
{
  if((P==NULL)||(G==NULL)) Erreur(-1); /* ne devrait jamais arriver */
  param_bfs *p;
  int i,x,y,u,v,n;

  if(j<0){ /* initialisation du premier chemin */
    n=G->n;
    if(P->aux==NULL) ALLOC2(P->aux,n,n);
    for(u=0;u<n;P->V[u++]=-1); /* vide le chemin */
    x=P->P[0];
    y=P->P[1];
    if((x<0)||(y<0)||(x>=n)||(y>=n)){ P->n=0; p=NULL; goto fin_0; } /* sommets inexistant */
    p=bfs(G,x,NULL); /* calcule le chemin */
    i=p->D[y];
    if(i<0){
      P->V[x]=0; /* x est en premiÃ¨re position dans P */
      P->n=1;
      goto fin_0;
    }
    /* on initialise aux[x] et aux[y] qu'une seule fois */
    P->aux[y][0]=0;
    j=-1; /* pour que la longueur soit correcte */
    goto fin_1;
  }

  n=P->n;

  if(j+1==n) return 0; /* si x=y, alors pas de prochain chemin */
  if(NextPath(G,P,j+1)) return 1; /* c'est le NextPath Ã  partir du voisin de v_j */

  /* Ici on ne peut pas augmenter v_(j+1)...y. Donc, il faut calculer
     un premier chemin depuis le prochain voisin v disponible de u=v_j
     et y qui Ã©vite x=P[0]-...-P[j]. Si cela ne marche pas avec le
     voisin v, il faut essayer le suivant, etc. */

  /* efface depuis P la fin du chemin P[j+1]...P[n-1] */
  /* pour ne garder que P[0]...P[j] */
  for(i=j+1;i<n;i++) P->V[P->P[i]]=-1;

  /* rem: ici t<d */
  p=new_param_bfs();
  ALLOC(p->D,G->n);
  y=P->P[n-1];
  u=P->P[j];
  i=-1;

  while((P->aux[u][0])&&(i<0)){ /* tant qu'il y a encore des voisins de u non testÃ©s */
    v=P->aux[u][(P->aux[u][0])--]; /* lit et enlÃ¨ve le dernier voisin dispo */
    if(P->V[v]>=0) continue; /* si le voisin est dans P[0] ... P[j] on le saute */
  
    /* initialise p->D: on doit le faire Ã  chaque appel */
    for(i=0;i<G->n;i++) p->D[i]=-1;
    /* on enlÃ¨ve P[0]...P[j] du graphe pour le BFS */
    for(i=0;i<=j;i++) p->D[P->P[i]]=-2;

    /* calcule un chemin de v Ã  y dans G\(P[0]-...-P[j]) */
    bfs(G,v,p);

    /* a-t-on trouvÃ© un chemin ? */
    i=p->D[y]; /* si i>=0, c'est oui */
  }

  /* ici i=distance de u=v_j Ã  y */
  if(i<0){ /* on a pas trouvÃ© de chemin */
  fin_0:
    free_param_bfs(p);
    return 0;
  }

  /* on a trouvÃ© un chemin, on met Ã  jour P */
 fin_1:
  P->n=i+j+2; /* nouvelle longueur */
  /* ajoute Ã  P[0]...P[j] le chemin de P[j+1] Ã  y en partant de y */
  while(i>=0){
    P->P[i+j+1]=y;
    P->V[y]=i+j+1;
    y=p->P[y];
    i--;
  }

  /* initialise aux[u] pour tous les sommets u de P[j+1] Ã  y non
     compris. Pour P[j] on le fait au moment de la lecture de v dans
     aux[u], et pour y on le fait une seule fois Ã  la crÃ©ation
     (avec aux[y][0]=0) */

  for(i=j+1,n=P->n-1;i<n;i++){ /* si j=-1, alors on fait pour x aussi */
    u=P->P[i];
    P->aux[u][0]=0;
    for(v=0;v<G->d[u];v++){
      j=G->L[u][v]; /* j=voisin de u */
      /* si pas dans P on l'ajoute comme voisin possible */
      if(P->V[j]<0) P->aux[u][++P->aux[u][0]]=j;
    }
  }

  free_param_bfs(p);
  return 1;
}


#define LCONF 9  /* nb de bits pour coder un sommet dans le graphe des conflits */
#define CONFMAX (1<<LCONF) /* = 2^LCONF = nb max de sommets du graphe des conflits */
#define CONFMASK (CONFMAX-1) /* = CONFMAX-1 = mask pour un sommet */
#define CONFC 2 /* constante pour modifier le code des noeuds Ã  0. Il faut CONFC>1 */
#define NCMAX 512 /* taille maximum de la liste de composantes maximales. NCMAX<CONFMAX */

/* ensemble de variables pour gÃ©rer le graphe des conflits */
typedef struct{
  int n;   /* nb de sommets du graphe d'origine, sert pour la rÃ¨gle 5 */
  int nbi; /* nb de noeuds avec code=-1 (indÃ©terminÃ©e) */
  int nbzi;/* nb de noeuds de code 0 indÃ©pendants */
  int outmem; /* 1 ssi il y a eut un dÃ©passement de CONFMAX ou NCMAX */
  int paire[CONFMAX]; /* 1er noeud de la paire xy */
  int path[CONFMAX];  /* 1er noeud du chemin P, ne sert que pour PrintConflit() */
  int code[CONFMAX];  /* code du noeud i. Valeurs possibles: -1,0,1 */
  int nbc[CONFMAX];   /* nb de noeuds de la paire i dont le code est <1 */
  int *comp[CONFMAX];  /* liste des noeuds de la composante de i */
  int *cmax[NCMAX+1]; /* liste des composantes maximales */
  int tmax[NCMAX+1]; /* taille des composantes maximales */
  int ncmax; /* taille de la liste cmax */
  int nodemax; /* noeud correspondant Ã  la composante maximale de la liste cmax */
  int valmax; /* valeur (VRAI/FAUX) des composantes maximales */
  int tcomp[CONFMAX]; /* taille de comp[i] */
  int x[CONFMAX],y[CONFMAX]; /* paire du noeud i, ne sert que pour PrintConflit() */
  graph* G; /* graphe des conflits */
  /* Attention! les noeuds des listes d'adjacence de G peuvent Ãªtre >
     CONFMAX et donc aussi > G->n. On s'en sert pour coder un type
     d'arÃªte. Donc si u est un noeud de G et v=G->L[u][i], alors le
     i-Ã¨me voisin de u est (v&CONFMASK) et le type de l'arÃªte uv est
     (v>>LCONF). */
} conflit;


void PrintConflit(const conflit* const c)
/*
  Affiche le graphe des conflits, et rien si c->G->n=0.
*/
{
  if(c->G->n==0) return;

  int u,v,i,t,d;
  string const T1="|=><";
  string const T2="-01*";
  string const T3="_/.";

  // UTF-8:
  // string const T3="â”â”›o";
  // \xe2\x94\x80: â”€
  // \xe2\x94\x81: â”
  // \xe2\x94\x98: â”˜
  // \xe2\x94\x99: â”™
  // \xe2\x94\x9a: â”š
  // \xe2\x94\x9b: â”›

  printf("code (x,y) [comp] node-path-paire: node[type]\n");
  for(u=0;u<c->G->n;u++){
    t=c->code[u];
    if(c->nbzi>0){ /* si on a des zÃ©ros indÃ©pendants */
      if(t==0) t=2; /* c'est un zÃ©ro indÃ©pendant */
      else if(t==CONFC) t=0; /* c'est un zÃ©ro marquÃ© -> remet sa valeur d'origine */
    }
    printf(" %c ",T2[t+1]);
    if(c->paire[u]==u)
      printf("(%i,%i) [",c->x[u],c->y[u]);
    else printf("\t [");
    for(v=0;v<c->tcomp[u];v++){
      printf("%i",c->comp[u][v]);
      if((c->n>9)&&(v<c->tcomp[u]-1)) printf(",");
    }
    printf("]\t%s%i",(v<5)?"\t":"",u);
    if(u>c->path[u]) printf("%c \t:",T3[1]); else{
      printf("%c",T3[0]);
      if(u>c->paire[u]) printf("%c\t:",T3[1]); else{
	printf("%c%c\t:",T3[0],T3[2]);
      }
    }
    for(i=0,d=min(c->G->d[u],WIDTH);i<d;i++){
      v=c->G->L[u][i];
      t=(v>>LCONF);
      printf(" %i%c",v&CONFMASK,T1[t]);
    }
    if(c->G->d[u]>WIDTH) printf("...");
    printf("\n");
  }
  printf("#nodes in conflict graph: %i\n",c->G->n);
  printf("#heavy indep. components: %i\n",c->nbzi);
  printf("#unspecified values: %i\n",c->nbi);
  if(c->outmem){
    printf("!!! Out of Memory !!!\n");
    if(c->ncmax>=NCMAX)
      printf("#nodes in conflit graph exceeded (%i)\n",CONFMAX);
    else
      printf("#maximal components exceeded (%i)\n",NCMAX);
  }
  return;
}


void ps1_delxy(const conflit* const c,int const w)
/*
  Supprime du graphe des conflits c la derniÃ¨re paire crÃ©ee, et met
  ainsi Ã  jour c. Ici w est l'indice du premier noeud de la derniÃ¨re
  paire dans le graphe des conflits. Il faut supprimer tous les noeuds
  u >= w et tous les arcs entre les noeuds u et v<w.

  La liste des composantes maximales de la derniÃ¨re paire sera effacÃ©e
  Ã  la fin de la paire courante (voir le code aprÃ¨s nextxy:).
*/
{
  int u,i,v,j;

  for(u=w;u<c->G->n;u++){
    for(i=0;i<c->G->d[u];i++){ /* pour tous les arcs v->u */
      v=c->G->L[u][i]&CONFMASK; /* v est le i-Ã¨me voisin de u */
      if(v<w) /* seulement pour les noeuds v avant w */
	for(j=0;j<c->G->d[v];j++)
	  /* on cherche et supprime u de la liste de v */
	  if((c->G->L[v][j]&CONFMASK)==u)
	    c->G->L[v][j]=c->G->L[v][--(c->G->d[v])]; /* supprime u */
    }
    free(c->G->L[u]); /* supprime la liste de u */
    free(c->comp[u]); /* supprime la composante de u */
  }

  c->G->n=w; /* met Ã  jour le nombre de noeuds du graphe des conflits */
  return;
}


int ps1_addmax(int* C,int const t,conflit* const c)
/*
  Ajoute une composante C (de taille t) Ã  la liste des composantes
  maximales dÃ©jÃ  rencontrÃ©es, et maintient cette propriÃ©tÃ©. Il faut
  que la liste (c->cmax) soit de taille suffisante pour accueillir une
  nouvelle composante. La fonction renvoie VRAI ssi la composante C a
  Ã©tÃ© ajoutÃ©e Ã  la liste.

  Algorithme:
    1. Pour toute composante X de L (=la liste):
       1.1. Si C est inclue ou Ã©gale Ã  X alors FIN
       1.2. Si C contient X alors supprimer X de L
    2. Ajouter C Ã  L

  On remarque aussi qu'il n'est pas possible de rencontrer plusieurs
  fois la mÃªme composante C avec des valeurs diffÃ©rentes. Si cela
  arrive, c'est pour deux chemins diffÃ©rents, disons Q1 et Q2. Soit X
  les voisins de C qui sont Ã  la fois dans Q1 et Q2.  X se dÃ©compose
  en segments, chacun Ã©tant un sous-chemin de Q1 inter Q2. Tout sommet
  voisin de C doit Ãªtre dans X sinon C aurait un sommet de trop. Donc
  tout chemin de G contenant X tel que C est une composante de G\P,
  doit Ãªtre parallÃ¨le Ã  Q1 et Q2. En particulier Q1 et Q2 sont
  parallÃ¨les ... et alors [A FINIR] ? Bon, bah en fait c'est
  possible. Q1 peut se rÃ©duire Ã  une arÃªte xy et Q2 Ã  xzy (un sommet
  de plus). Et avec Q1 c'est vrai, mais avec Q2 cela devient faux. On
  a des contre-exemple avec des sommets (z) de degrÃ© deux.
*/
{
  int **L=c->cmax,n=c->ncmax; /* L=liste et n=|L|*/
  int *T=c->tmax; /* T=taille de chaque composante */
  int i,r;

  /* on passe en revue chaque composante la liste L */
  for(i=0;i<n;i++){
    r=SetCmp(C,L[i],t,T[i]); /* compare C et L[i] */
    if(r&6) return 0; /* C est strictement contenu (&4) ou Ã©gale (&2) Ã  L[i] */
    if((r&8)&&(n>0)){ /* L[i] est contenu (strictement) dans C et L non vide */
      free(L[i]); /* libÃ¨re ce tableau de sommets prÃ©cÃ©demment allouÃ© par un ps1_addmax() */
      /* si L[i] dernier Ã©lÃ©ment de L (i=n-1), alors il n'y a plus rien Ã  faire */
      n--;
      if(i<n){ /* si c->cmax[i] pas dernier Ã©lÃ©ment, il faut dÃ©placer le dernier */
	L[i]=L[n]; L[n]=NULL; /* pour Ã©viter plus tard un double free. NB: i<>n */ 
	T[i]=T[n]; /* nombre de sommets dans L[i] */
	i--; /* pour recommencer avec ce nouveau L[i] */
      }
    }
  }

  /*ajoute C Ã  la fin de L */
  ALLOCZ(L[n],t,C[_i]); /* alloue et copie C dans L[i] */
  T[n]=t; /* taille de C */
  c->ncmax=n+1; /* nouvelle taille de L */
  return 1; /* on a ajoutÃ© C Ã  L */
}


int ps1_push(int x,int const v,conflit* const c)
/*
  Affecte la valeur v (=0 ou 1) dans le noeud x et propage, en
  appliquant les rÃ¨gles dÃ©crites ci-aprÃ¨s, Ã  tous ses voisins et
  rÃ©cursivement Ã  tout le graphe. Renvoie 1 si une contradiction
  apparaÃ®t, ou bien si tous les noeuds de la paire de x ont pour code
  1. Sinon, on renvoie 0 (terminaison standard).

  Attention! Il faut appliquer cette fonction que si la paire de x est
  complÃ¨te. Donc il ne faut pas l'appliquer si on vient de crÃ©er le
  noeud x, car on n'est pas sÃ»r de ne pas faire un "Out of Memory" un
  peu plus tard sur la mÃªme paire.

  Effet de bord: met Ã  jour plusieurs champs de la variable c, mais ne
  modifie pas le graphe c->G.

  RÃ¨gles pour l'arc x->y:

  x=noeud x
  y=noeud y
  v=valeur 0 ou 1 qui vient d'Ãªtre Ã©crite dans x
  t=type de l'arÃªte: 0=(Tx|Ty), 1=(|Tx|=|Ty|), 2=(Tx<Ty), 3=(Tx>Ty)
  c=graphe des conflits

  (Attention! "Tx|Ty" signifie que les composantes Tx et Ty sont
  disjointes.)

  rÃ¨gle 1: si Tx|Ty (disjoint) et v=0, alors Ã©crire 1 dans y
  rÃ¨gle 2: si Tx=Ty, alors Ã©crire v dans y
  rÃ¨gle 3: si Tx<Ty et v=0, alors Ã©crire 0 dans y
  rÃ¨gle 4: si Tx>Ty et v=1, alors Ã©crire 1 dans y
  rÃ¨gle 5: si Tx|Ty, v=1 et |Tx|+|Ty|=n, alors Ã©crire 0 dans y

  On applique Ã©galement la "rÃ¨gle du dernier -1", Ã  savoir si la paire
  de x, aprÃ¨s avoir Ã©crit v, possÃ¨de exactement 1 seul noeud de valeur
  -1, alors on Ã©crit 0 dans ce noeud.

  La "rÃ¨gle du max" et la "rÃ¨gle de l'influence des voisins" sont
  appliquÃ©es plus directement par la fonction ps1() lors de la
  crÃ©ation d'une paire. Elles n'ont pas lieu d'Ãªtre lors de la
  propagation.
*/
{
  int i,d,y,t;

  if(c->code[x]==v) return 0; /* rien Ã  faire */
  if(c->code[x]==1-v) return 1; /* contradiction */
  /* ici on a code[x]==-1 */
  c->code[x]=v; /* Ã©crit la valeur v */
  c->nbi--; /* et une valeur indÃ©terminÃ©e en moins ! */
  t=(c->nbc[c->paire[x]] -= v); /* diminue seulement si on a Ã©crit 1 */
  if(t==0) return 1; /* la paire de x est bonne, elle contient que des 1 ! */

  /* applique les rÃ¨gles 1-5 Ã  tous les arcs sortant de x */
  for(i=0,d=c->G->d[x];i<d;i++){ /* pour tous les voisins de x */
    y=c->G->L[x][i]; /* y=i-Ã¨me voisin de x */
    t=(y>>LCONF); /* t=type de l'arc x->y: 0,1,2,3 */ 
    y &= CONFMASK; /* y=numÃ©ro du noeud voisin de x */

    /* applique les rÃ¨gles */
    switch(t){
    case 0:
      if((v==0)&&(ps1_push(y,1,c))) return 1; /* rÃ¨gle 1 */
      if((v==1)&&(c->tcomp[x]+c->tcomp[y]==c->n)&&(ps1_push(y,0,c))) return 1; /* rÃ¨gle 5 */
      break;
    case 1: if(ps1_push(y,v,c)) return 1; break; /* rÃ¨gle 2 */
    case 2: if((v==0)&&(ps1_push(y,0,c))) return 1; break; /* rÃ¨gle 3 */
    case 3: if((v==1)&&(ps1_push(y,1,c))) return 1; break; /* rÃ¨gle 4 */
    }
  }

  /* rÃ¨gle du dernier -1 ? */
  if(c->nbc[c->paire[x]]==1){
    /* on cherche alors l'unique noeud x de la paire courante qui est
       de code < 1.  NB: ce noeud existe forcÃ©ment */
    x=c->paire[x]; /* x=premier noeud de la paire de x */
    while(c->code[x]==1) x++; /* passe au suivant si le code est 1 */
    return ps1_push(x,0,c); /* Ã©crit 0 dans le noeud x */
  }

  return 0;
}


/* pour le dÃ©bugage de PS1() */
/* LEVEL=niveau de rÃ©cursion, POS=numÃ©ro de ligne */
int LEVEL=0;
#define PRINTS do{			       	\
    int _i;			       		\
    printf("%03i:%02i  ",++POS,LEVEL);		\
    for(_i=0;_i<3*LEVEL;_i++) printf(" ");	\
  }while(0)


int PS1(graph* const G,path* const P,int const version){
/*
  P est un chemin d'un graphe G, et G\P est formÃ© d'une seule
  composante connexe. Il faut G et P <> NULL (initialisÃ©s avec
  new_graph() et new_path()). Renvoie 1 si P peut "sÃ©parer" G (voir
  explications ci-dessous). Il y a une optimisation avec le graphe des
  conflits si version>0. On utilise cette fonction comme ceci:

  path *P=new_path(G,NULL,G->n); // P=chemin vide, sans sommet
  int r=PS1(G,P); // r=0 ou 1
  free_path(P);
  
  Effet de bord: G->int1 retourne le nombre de tests (nombre de
  paires, nombre de chemins testÃ©s, et nombre de passes dans le graphe
  des conflits). Les autres champs de G et de P ne sont pas modifiÃ©s.

  Le paramÃ¨tre "version" indique la variante du test:
  - version=0: sans le graphe des conflits
  - version=1: avec le graphe des conflits
  - version=2: comme version=1 mais sans le graphe des conflits lors de la rÃ©cursivitÃ©
  - version=3: comme version=1 mais avec l'Ã©criture de valeurs dans le graphe des conflits

  AmÃ©liorations possibles:

  - Si G n'est pas biconnexe, alors on pourrait tester si toutes ses
    composantes biconnexes sont bien Ã©valuÃ©e Ã  vraie. Si P
    n'intersecte pas une composante biconnexe B de G, alors il faut
    Ã©valuer PS1(B,{}).

  - G\P Ã©tant connexe, on pourrait dÃ©jÃ  supprimer les sommets de P qui
    forment un segment de P contenant une extrÃ©mitÃ© de P et qui n'ont
    pas de voisin dans G\P. En particulier, si une des extrÃ©mitÃ©s de P
    est de degrÃ© 1, on peut la supprimer.

  - PrivilÃ©gier les paires de sommets qui ne sont pas adjacents (pour
    diminuer la taille les composantes et avoir des chemins plus
    longs). Plus gÃ©nÃ©ralement, on pourrait privilÃ©gier les paires de
    sommets les plus distants. Ceci dit, on ne gagne probablement pas
    grand chose, et une renumÃ©rotation alÃ©atoire des sommets devrait
    suffir pour ne pas traitÃ© les paires de sommets voisins en
    prioritÃ©.

  - On pourrait tester des cas simples pour G: arbre (tester si m=n-1,
    on sait que G est connexe), clique (tester si m=n(n-1)/2: si n<=4
    sommets alors vraie, sinon faux). (Ces tests sont dÃ©jÃ 
    implÃ©mentÃ©s). Plus dur: G est outerplanar. En fait, si G est un
    arbre de cycles (chaque composante connexe est un cycle ou un K4),
    alors le test est vrai. C'est donc vrai en particulier si
    m=n. Pour tester si G est un arbre de cycle, il suffit de faire un
    DFS, puis de vÃ©rifier que si (u,x) et (u,y) sont deux arÃªtes qui
    ne sont pas dans l'arbre du DFS, alors x et y ne sont pas ancÃªtres
    l'un de l'autre (??? Pourquoi ???).

  - Pour tous les chemins possibles testÃ©s rÃ©cursivement pour G, ne
    tester effectivement que ceux qui correspondent Ã  des
    sous-ensembles de sommets diffÃ©rents puisque le rÃ©sultat sera le
    mÃªme (mÃªmes composantes connexes). Pour cela, il faut gÃ©rer une
    table pour mÃ©moriser les sous-ensembles testÃ©s. Notons que si un
    chemin est induit alors cela ne sert Ã  rien de le mÃ©moriser, il ne
    pourra jamais Ãªtre rencontrÃ© de nouveau. On pourrait coder un
    sous-ensemble de n sommets par un entier sur n bits (n < 32 ou 64
    donc). La recherche/insertion pourrait Ãªtre une recherche dans un
    arbre binaire de recherche.

 EXPLICATIONS:

  Soit P et Q deux chemins d'un graphe G. On dit que Q est parallÃ¨le Ã 
  P, notÃ© Q//P, s'il n'y a pas d'arÃªte entre P\Q et G\(Q u P).

  Quelques propriÃ©tÃ©s:

  - La relation // n'est pas symÃ©trique.
  - Tout chemin est parallÃ¨le au chemin vide.
  - Tout chemin est parallÃ¨le Ã  lui-mÃªme.
  - La relation // n'est pas transitive.

  Pour la derniÃ¨re proposition, on peut vÃ©rifier en fait que si R//Q
  et Q//P, alors on a R//P sauf s'il existe une arÃªte uv telle que u
  in (P inter Q)\R et v in Q\(R u P).

  Soit P et Q deux chemins de G. On dit que Q dÃ©rive de P, notÃ© Q///P,
  s'il existe une suite P_0,P_1,...,P_n de chemins de G avec P_0=Q et
  P_n=P tels que P_{i-1}//P_i pour chaque i=1..n.

  Quelques propriÃ©tÃ©s:

  - Si Q//P, alors Q///P. En particulier, Q///{} puisque Q//{}.
  - On peut avoir R//Q//P, soit R///P, sans avoir R//P (cf. ci-dessus).
  - Si R///Q et Q///P, alors R///P.

  On dit que Q///P dans un graphe valuÃ© (G,w) si tous les chemins
  P_0,...,P_n (en particulier Q et P) sont des plus courts chemins
  selon w.

  Soit P un chemin de G et w une valuation de G. On dÃ©finit le
  potentiel pour P selon w comme score(P,w) := max_C { w(C)*|V(G)|+|C|
  } oÃ¹ le maximum est pris sur toute composante connexe C de G\P.

  Lemme 1. Supposons que G\P a une seule composante connexe, et soit Q
  un chemin de G parallÃ¨le Ã  P diffÃ©rent de P. Alors, pour chaque
  valuation w de G, soit P est un demi-sÃ©parateur de G ou bien
  score(Q,w) < score(P,w).

  Preuve. Supposons que P ne soit pas un demi-sÃ©parateur de G pour la
  valuation w. Soit C la composante de G\P, et posons n = |V(G)|. Par
  dÃ©finition, score(P,w) = w(C)*n + |C|. On a w(C) > w(G)/2, et donc
  w(P) < w(G)/2. Il suit que w(P) < w(C). Soit C' une composante de
  G\Q telle que w(C')*n+|C'| = score(Q,w). Comme Q est parallÃ¨le Ã  P,
  soit C' est contenue dans C soit C' est contenue dans P. En effet,
  si C' intersecte C et P, alors C' contient une arÃªte uv avec u in C
  et v in P. Bien sÃ»r uv not in Q. Cela contredit le fait qu'il existe
  pas d'arÃªte entre P\Q et G\(QuP).

  Si C' est contenue dans P, alors w(C') <= w(P) < w(C). Il suit que
  w(C') <= w(C)-1, soit w(C')*n <= w(C)*n - n. Clairement |C'| < n +
  |C|. D'oÃ¹ w(C')*n+|C'| < w(C)*n+|C|, soit score(Q,w) < score(P,w).

  Si C' est contenue dans C, alors w(C') <= w(C) et |C'| < |C| car
  Q<>P. Il suit que w(C')*n + |C'| < w(C)*n + |C|, soit score(Q,w) <
  score(P,w).

  Dans les deux cas nous avons prouvÃ© que score(Q,w) < score(P,w).
  QED

  Soit G un graphe et P un chemin de G tel que G\P est composÃ© d'une
  seule composante connexe (en particulier G est connexe). On dÃ©finit
  PS1(G,P) le prÃ©dicat qui est VRAI ssi pour toute pondÃ©ration w de G
  telle que P est un plus court chemin il existe dans (G,w) un chemin
  demi-sÃ©parateur qui dÃ©rive de P.

  Lemme 2. G est dans PS1 ssi PS1(G,{}) = VRAI.

  Preuve. En effet, en rÃ©Ã©crivant la dÃ©finition de PS1(G,{}) on dÃ©duit
  que PS1(G,{}) est VRAI ssi pour toute pondÃ©ration w de G il existe
  dans (G,w) un chemin demi-sÃ©parateur qui dÃ©rive du chemin vide (tout
  chemin dÃ©rive du chemin vide). Notons que c'est nÃ©cessairement un
  plus court chemin de (G,w). C'est prÃ©cisemment la dÃ©finition de la
  classe PS1. QED

  L'objectif est d'avoir un test notÃ© ps1(G,P) qui implÃ©mente
  PS1(G,P), disons qu'il s'en approche. La propriÃ©tÃ© souhaitÃ©e est que
  si ps1(G,P) est VRAI, alors PS1(G,P) aussi. En particulier, si
  ps1(G,{}) est VRAI, alors G est dans PS1.

  Algorithme pour le test ps1(G,P):

  On renvoie VRAI s'il existe une paire x,y de sommets oÃ¹ y n'est pas
  dans P telle que tout chemin Q de x Ã  y:
  1. Q est parallÃ¨le Ã  P, et
  2. pour toute composante C de G\(QuP), ps1(CuQ,Q)=VRAI.

  Lemme 3. Si ps1(G,P)=VRAI, alors PS(G,P)=VRAI.

  Preuve [A FINIR]. Par induction sur le nombre de sommets hors de
  P. Soit C est l'unique composante G\P. Si |C|=0, alors P est un
  demi-sÃ©parateur de G et donc PS(G,P) est VRAI. Supposons le lemme
  vrai pour tout entier < |C|.

  On suppose donc que ps1(G,P) est VRAI. Soit x,y la paire de sommets
  telle que y n'est pas dans P et oÃ¹ tous les chemins de x Ã  y sont
  parallÃ¨les Ã  P. En particulier, pour chaque valuation w, oÃ¹ P est un
  plus court chemin, tout plus court chemin Q selon w entre x et y est
  parallÃ¨le Ã  P. Comme Q est diffÃ©rent de P (Ã  cause du choix de y),
  on peut appliquer le lemme 1, et donc soit P est un demi-sÃ©parateur,
  soit score(Q,w) < score(P,w). On peut supposer qu'on est pas dans le
  premier cas, c'est-Ã -dire que P n'est pas un demi-sÃ©parateur de G,
  puisque sinon PS(G,P)=VRAI et le lemme est prouvÃ©.

  Si Q est un demi-sÃ©parateur pour G, alors PS(G,P) est VRAI puisque Q
  est parallÃ¨le Ã  P. Supposons que Q n'est pas un demi-sÃ©parateur pour
  G, et soit C' la composante de G\Q telle que w(C')>w(G)/2.

  On peut appliquer l'induction sur (C'uQ,Q) car comme Q<>P,
  |C'|<|C|. Posons G'=C'uQ. D'aprÃ¨s le test ps1(G,P), ps1(G',Q) est
  VRAI. Donc par induction PS(G',Q)=VRAI et G' contient un chemin
  demi-sÃ©parateur pour la valuation w, disons P', parallÃ¨le Ã 
  Q. Montrons d'abord que dans G, P' est parallÃ¨le Ã  P.

  ...

  w(C')<w(G)/2 ...

  QED

  Lemme 4. Si ps1(G,P)=VRAI, alors soit il existe une paire x,y de
  sommets de G telle que tout chemin de x Ã  y contient P, ou bien il
  n'existe pas de sommet de P ayant trois voisins inclus dans un cycle
  de G\P. [PAS SÃ›R, A VÃ‰RIFIER]

  Preuve. [A FINIR] Soit x,y une paire de sommets de G avec y pas dans
  P telle que tous les chemins de x Ã  y soient parallÃ¨les Ã  P. Soient
  Q un tel chemin. Supposons que Q ne contient pas P. ...  QED

  Dit autrement, si on n'a pas cette propriÃ©tÃ©, alors ps1(G,P)=FAUX et
  il est inutile de tester tous les chemins Q possibles.  Remarquons
  que si G est 2-connexe, alors il ne peut avoir de paire x,y de
  sommets (avec y pas dans P) oÃ¹ tout chemin de x et y contient P. A
  montrer: ps1(G,P)=VRAI ssi tout les composantes 2-connexes G' de G
  on a ps1(G',P inter G)=VRAI ...

  [A VOIR]

  - u := ps1(CuQ,Q)
  - ajoute (C,u) Ã  la liste L des composantes maximales (voir ps1_addmax)
  - si u = FAUX, ajouter un nouveau noeud au graphe des conflits (GC)
  - recommencer avec le chemin Q suivant de mÃªmes extrÃ©mitÃ©s xy (s'il existe)
  - s'il n'y a plus de tel chemin:

    On essaye d'appliquer les trois rÃ¨gles (max, influence, dernier):
    - la rÃ¨gle du max en tenant compte de la liste L.
      si |L|<>1, alors on ne peut pas appliquer la rÃ¨gle du max
      sinon, L={(C,u)}
        si u = VRAI, on peut supprimer la paire xy de GC
        sinon, alors on peut appliquer la rÃ¨gle du max
    - la rÃ¨gle d'influence des voisins
    - la rÃ¨gle du dernier -1

    Si l'application d'une de ces rÃ¨gles (avec ps1_push) produit une
    contradiction dans GC, alors on a trouvÃ© une bonne paire xy, et on
    renvoie VRAI.  Sinon, s'il existe une autre paire xy on recommence
    avec cette pnouvelle paire.

  A la fin du traitement de toutes les paires:

  - soit il reste des indÃ©terminÃ©es (-1), et il faut les Ã©liminer. On
    les force d'abord Ã  0. S'il y a une contradiction (en faisant des
    ps1_push), on en dÃ©duit que la valeur doit Ãªtre 1 (et on fait
    ps1_push). Sinon, on force la valeur Ã  1. Si y a une
    contradiction, on dÃ©duit que la valeur doit Ãªtre 0 (et on fait
    ps1_push). Sinon, s'il y a une valeur initialement indÃ©terminÃ©e
    qui passe Ã  la mÃªme valeur pour le forcage Ã  0 et Ã  1, on Ã©limine
    cette indÃ©terminÃ©e (et on fait ps1_push). On fait ainsi pour
    chaque indÃ©terminÃ©e. Chaque fois qu'une indÃ©terminÃ©e est Ã©liminÃ©e,
    on recommence la recherche sur l'ensemble des indÃ©terminÃ©es (et
    pas seulement sur celles qui restent). Si on arrive ainsi Ã 
    Ã©liminer toutes les indÃ©terminÃ©es on peut passer Ã  la
    suite. Sinon, on ne peut pas faire grand chose Ã  part essayer tous
    les systÃ¨me possibles ... Voir MINISAT+ (systÃ¨me pseudo boolÃ©ens)
    et Sugar qui transforme du systÃ¨me linÃ©aire ou CSP en SAT.

  - soit il n'y a plus d'interminÃ©es. Dans ce cas on peut dÃ©duire un
    systÃ¨me d'Ã©quations linÃ©aire indÃ©pendantes oÃ¹ les inconnues sont
    les poids des sommets. Si le systÃ¨me n'a pas de solution, alors
    renvoyer VRAI. Sinon, la solution trouvÃ©e peut renseigner sur une
    valuation possible pour prouver que G est Ã©ventuellement pas dans
    PS1(G,P).

*/
  G->int1=1; /* compte les tests */

  DEBUG(
	LEVEL++;
	int u;int v;
	if(G==GF) LEVEL=POS=0;
	printf("\n");
	PRINTS;printf("version=%i G=%p n=%i ",version,G,G->n);PRINTT(P->P,P->n);
	PRINTS;printf("G=");
	for(u=0;u<G->n;u++){
	  printf("%i:[",u);
	  for(v=0;v<G->d[u];v++){
	    printf("%i",G->L[u][v]);
	    if(v<(G->d[u])-1) printf(" ");
	  } printf("] ");
	} printf("\n");
	);
  
  int const n=G->n;

  /* Ici on Ã©limine un certain nombre de cas faciles Ã  tester.
     Attention ! vÃ©rifier que G est dans PS1 dans ces cas lÃ  ne suffit
     pas (ce qui revient Ã  vÃ©rifier que PS1(G,{})=VRAI). Il faut Ãªtre
     certain qu'on a en fait PS1(G,P)=VRAI. */

  if(n-(P->n)<3){
    DEBUG(PRINTS;printf("-> n-|P|<3, moins de 3 sommets hors P\n"););
    DEBUG(LEVEL--;);
    return 1;
  }
  /* Par hypothÃ¨se P est lÃ©ger, donc la composante hors de P est
     lourde. Il suffit de prendre un chemin entre ces au plus deux
     sommets hors de P. */

  /* lors d'un appel rÃ©cursif, le nombre d'arÃªtes G->m est dÃ©jÃ  mis Ã 
     jour car G provient alors du rÃ©sultat de ExtractSubgraph(). Donc
     nb_edges(G) est calculÃ© en temps constant dÃ¨s le 2e appel. */

  if((!P->n)&&(n<6)&&(nb_edges(G)<10)){
    DEBUG(PRINTS;printf("-> |P|=0 et n<6 et pas Kâ‚…\n"););
    DEBUG(LEVEL--;);
    return 1;
  }

  /* ici on a au moins 3 sommets de G qui sont hors de P, et P est non
     vide. Alors, comme P contient au moins deux somemts, cela fait
     que G possÃ¨de au moins 5 sommets. */

  if(nb_edges(G)==((n*(n-1))>>1)){
    DEBUG(PRINTS;printf("-> clique avec > 2 sommets hors P\n"););
    DEBUG(LEVEL--;);
    return 0;
  }
  /* Ici G est clique avec n>4 et n-|P|>2. Dans le cas oÃ¹ tous les
     poids sont Ã  1 (sommets et arÃªtes) alors, il n'existe aucun
     chemin Q parallÃ¨le P permettant de progresser. Notons qu'une
     clique G avec n sommets donne PS1(G,P)=VRAI s'il y a 0, 1 ou 2
     sommets dans G\P. */

  if(nb_edges(G)<=n){
    DEBUG(PRINTS;printf("-> m<=n\n"););
    DEBUG(LEVEL--;);
    return 1;
  }
  /* Dans ce cas, il s'agit d'un arbre avec un cycle. Soit C la
     composante lourde de G\P, et x une des extrÃ©mitÃ©s de P. La
     composante C est connectÃ©e Ã  P par au plus deux sommets, disons u
     et v (u=v possible) et u le plus prÃ¨s de x. Soit y le voisin de v
     dans C. On dÃ©finit alors P' comme le plus court chemin de G
     allant de x Ã  y. Il est facile de voir que P' longe P depuis x
     puis entre dans C soit par u soit par v. Dans les deux cas les
     sommets de C ne peuvent Ãªtre connectÃ©s Ã  aucun sommet de P\P',
     les deux seules arÃªtes connectant C Ã  P Ã©tant dÃ©truite par P'. */

  /* ici G possÃ¨de au moins 5 sommets et 6 arÃªtes. */

  int x,y,i,u,v,w,d;
  path *Q=new_path(G,NULL,n);
  path *R=new_path(G,NULL,n);
  param_dfs *p=new_param_dfs(n); /* p->C n'est pas allouÃ© */
  graph* C;

  ALLOC(p->C,n);   /* pour le DFS avec sommets supprimÃ©s */
  NALLOC(int,T,n); /* pour la composante C de G */
  NALLOC(int,M,n); /* pour la compatiblitÃ© des chemins P et Q */

  conflit c; /* ensembles de variables pour gÃ©rÃ©r le graphe des conflits */
  int npaire,npath,goodxy;
  /* npaire = 1er noeud dans le graphe des conflits de la paire courante */
  /* npath = 1er noeud dans GC du chemin courant pour la paire courante */
  /* goodxy = 1 ssi la paire xy est bonne, ie. tous les chemins et comp. sont ok */

  if(version>0){
    c.n=n; /* nombre de sommets du graphe G */
    c.G=new_graph(CONFMAX); /* graphe des conflits, alloue G->d et G->L */
    c.G->n=c.nbi=c.nbzi=c.ncmax=c.outmem=0; /* c.G->n=nb de noeuds dÃ©jÃ  crÃ©es */
  }

  /* pour toutes les paires de sommets x,y de G avec y pas dans P */

  for(x=0;x<n;x++)
    for(y=0;y<n;y++){

      if(P->V[y]>=0) continue; /* y ne doit pas Ãªtre dans P */
      if((P->V[x]<0)&&(y<=x)) continue; /* si x et y pas dans P, ne tester que le cas x<y */

      /* ici on a soit:
	 1) x dans P et y pas dans P, ou bien
	 2) x et y pas dans P et x<y
      */

      (G->int1)++; /* +1 pour chaque paire testÃ©e */
      goodxy=1; /* par dÃ©faut on suppose la paire xy comme bonne */

      /* calcule le 1er chemin entre x et y */
      Q->P[0]=x; Q->P[1]=y; /* initialise les extrÃ©mitÃ©s du chemin Q */
      if(!NextPath(G,Q,-1)) goto fin_ps1; /* fin si pas de chemin x->y, impossible si G connexe */
      if((version>0)&&(!c.outmem)){
	npaire=c.G->n; /* initialise le 1er noeud de la paire courante */
	c.nbc[npaire]=0; /* nombre de valeurs < 1 pour la paire courante */
	/* on efface la liste des composantes maximales, s'il y en avait */
	for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
	c.ncmax=0;
      }

      DEBUG(PRINTS;printf("Paire (%i,%i) G=%p n=%i\n",x,y,G,n););

      do{ /* pour tous les chemins Q entre x et y */
	(G->int1)++; /* +1 pour chaque sommet testÃ© */

	DEBUG(PRINTS;printf("Essai du chemin ");PRINTT(Q->P,Q->n););

	/* On vÃ©rifie que Q est parallÃ¨le avec P. Il faut qu'aucun
	   sommet de P\Q n'ait de voisin en dehors de P ou de Q. */

	for(i=0;i<P->n;i++){
	  if(Q->V[u=P->P[i]]>=0) continue; /* on est aussi dans Q */
	  /* ici on est dans P\Q */
	  for(v=0;v<G->d[u];v++){ /* vÃ©rifie chaque voisin v de u */
	    if(P->V[G->L[u][v]]>=0) continue; /* si v est dans P */
	    if(Q->V[G->L[u][v]]>=0) continue; /* si v est dans Q */
	    DEBUG(PRINTS;printf("-> chemin Q non parallÃ¨le Ã  P\n"););
	    goodxy=0; /* cette paire n'est pas bonne */
	    goto nextxy; /* aller Ã  la prochaine paire */
	  }
	}

	/* ici Q est parallÃ¨le Ã  P */

	DEBUG(PRINTS;printf("-> chemin Q parallÃ¨le Ã  ");PRINTT(P->P,P->n););

	/* on vÃ©rifie que pour chaque composante C de G\Q,
	   PS1(CuQ,Q)=VRAI. Notons que les sommets de P\Q sont
	   lÃ©ger. Il ne faut donc pas les considÃ©rer. Il faut donc
	   prendre les composantes de G\(QuP).*/

	/* on enlÃ¨ve de G les sommets de QuP pour calculer les
	   composantes de G\(QuP) */
	for(u=0;u<n;u++) /* initialise le dfs */
	  if((P->V[u]<0)&&(Q->V[u]<0)) p->C[u]=-1; else p->C[u]=-2;
	dfs(G,0,p); /* calcule les composantes de G\(QuP) */
	DEBUG(PRINTS;printf("#composantes dans G\\(QuP) Ã  tester: %i\n",p->nc););
	/* si aucune composante inutile de tester rÃ©cursivement ce
	   chemin, G\(QuP) est vide. On peut passer au prochain chemin */
	if(p->nc==0) goto nextQ;
	/* ici, il y a au moins une composante dans G\(QuP) */

	d=R->n=Q->n; /* d=nombre de sommets du chemin Q */
	for(u=0;u<d;u++) T[u]=Q->P[u]; /* T=liste des sommets de Q */
	if((version>0)&&(!c.outmem)) npath=c.G->n; /* initialise le chemin Q au noeud courant */

	DEBUG(PRINTS;printf("Q = ");PRINTT(T,d););

	/* pour chaque composante de G\Q */

	for(i=0;i<p->nc;i++){  /* T[d..v[=i-Ã¨me composante de G\(QuP) u Q */
	  for(u=0,v=d;u<n;u++) if(p->C[u]==i) T[v++]=u;
	  /* T[0..v[=sommets de Q u C */
	  /* T[0..d[=sommets de Q, T[d..v[=sommets de la composante */
	  /* NB: les sommets de T[d..v[ sont dans l'ordre croissant */

	  DEBUG(PRINTS;printf("C\\Q=CC(%i) = ",i);PRINTT(T+d,v-d););

	  C=ExtractSubgraph(G,T,v,1); /* crÃ©e C=G[T] avec v=|T| sommets */
	  /* Attention! Q n'est pas un chemin de C, mais de G. On crÃ©e
	     donc un nouveau chemin R dans C qui est Ã©quivalent Ã  Q */
	  for(u=0;u<v;u++) R->V[u]=-1; /* Attention! boucle avec v=C->n sommets */
	  for(u=0;u<d;u++) R->P[u]=C->pint1[Q->P[u]]-1;
	  for(u=0;u<d;u++) R->V[R->P[u]]=u;

	  DEBUG(PRINTS;printf("Q = ");PRINTT(R->P,R->n););

	  /* appel rÃ©cursif avec une nouvelle version 0 ou 1:
	     - si version=0: sans le graphe des conflits -> 0
	     - si version=1: avec le graphe des conflits -> 1
	     - si version=2: comme version=1 mais sans le graphe
	       des conflits lors de la rÃ©cursivitÃ© -> 0
	     - si version=3: comme version=1 mais avec l'Ã©criture
	       de valeurs dans le graphe des conflits -> 1
	  */
	  u=PS1(C,R,version%2);

	  DEBUG(PRINTS;printf("PS1(CuQ,Q)=%i\n",u););
	  DEBUG(if(c.outmem){PRINTS;printf("PROBLEME MEMOIRE !\n");});

	  G->int1 += C->int1; /* met Ã  jour le nombre de tests (ajoute au moins 1) */
	  free_graph(C); /* libÃ¨re C qui ne sert plus Ã  rien */

	  /* Ã  faire que si graphe des conflits et pas eut de problÃ¨me mÃ©moire */
	  if((version>0)&&(!c.outmem)){
	    
	    /* on vÃ©rifie si la mÃªme composante n'existe pas dÃ©jÃ  dans
	       la paire de c.G->n. Si c'est le cas, on passe Ã  la
	       prochaine composante. NB: la composante de c.G->n est
	       dans T[d...v[ et sa taille vaut v-d. */

	    for(w=npaire;w<c.G->n;w++)
	      if(SetCmp(c.comp[w],T+d,c.tcomp[w],v-d)&2)
		goto nextC; /* si mÃªme composante, alors prochaine composante */
	    /* ici w=c.G->n */
	    /* si u=0, w sera le nouveau noeud du graphe des conflits */

	    /* ajoute T[d..v[ Ã  la liste des composantes maximales (c.cmax) */
	    if(ps1_addmax(T+d,v-d,&c)){ c.nodemax=w; c.valmax=u; }
	    /* On stocke dans c->nodemax le noeud de la derniÃ¨re
            composante ajoutÃ©e Ã  la liste des composantes maximales.
            Si Ã  la fin du traitement de la paire xy la liste a
            exactement une seule composante, alors c->nodemax sera
            forcÃ©ment le noeud de la derniÃ¨re composante ajoutÃ©e et
            c->valmax sa valeur (VRAIE/FAUSSE). */
	    if(c.ncmax>=NCMAX){ /* dÃ©passement du nb de composantes maximales */
	      c.outmem=1; /* Ouf of Memory */
	      /* la paire xy n'est pas complÃ¨te. On l'enlÃ¨ve, sinon
	         l'analyse des conflits pourrait ne va pas Ãªtre correcte */
	      ps1_delxy(&c,npaire); /* supprime la derniÃ¨re paire crÃ©ee */
	      goodxy=0; /* cette paire n'est pas bonne */
	      goto nextxy; /* change de paire */
	    }
	  }
	  /* ici outmem=0 */

	if(u) goto nextC; /* si u=VRAI, on a finit le traitement de la composante */
	goodxy=0; /* une composante n'est pas bonne pour la paire xy */
	if(version==0) goto nextxy; /* si pas graphe de conflits, alors changer de paire */

	/* ici u=FAUX et version>0 */
	/* on va donc essayer d'ajouter le noeud w au graphe des conflits */
	/* ici pas de problÃ¨me de mÃ©moire, donc on ajoute le noeud w */

	(c.G->n)++; /* un noeud de plus */
	if(c.G->n==CONFMAX){ /* Out of Memory */
	  c.outmem=1;
	  ps1_delxy(&c,npaire); /* supprime la derniÃ¨re paire crÃ©ee, qui n'est pas bonne */
	  goto nextxy; /* ici goodxy=0 */
	}
	/* ici v=|T|, d=|Q|, w=nouveau noeud Ã  crÃ©er */
	
	/* Attention! ne pas utiliser ps1_push(w,...) juste aprÃ¨s la
	   crÃ©ation de w, car on est pas certain que la paire de w
	   sera complÃ¨te ... On ne peut faire ps1_push(w,...)
	   seulement aprÃ¨s le while(NextPath...). Idem pour
	   l'application de la rÃ¨gle du max ! */

	/* initialise le nouveau noeud w */

	DEBUG(PRINTS;printf("AjoÃ»t du noeud %i dans le graphe des conflits\n",w););

	c.x[w]=x; c.y[w]=y; /* mÃ©morise la paire x,y (pour l'affichage) */
	c.path[w]=npath; c.paire[w]=npaire; /* w rattachÃ© Ã  npath et Ã  npaire */
	c.code[w]=-1; c.nbi++; /* un noeud de plus Ã  valeur = -1 */
	(c.nbc[npaire])++; /* un noeud de plus Ã  valeur < 1 pour cette paire */
	c.G->d[w]=0; /* initialise le degrÃ© de w */
	ALLOC(c.G->L[w],CONFMAX); /* crÃ©e la liste des voisins de w */
	ALLOCZ(c.comp[w],v-d,T[d+(_i)]); /* crÃ©e & copie la composante de w */
	c.tcomp[w]=v-d; /* taille de la composante de w */
	
	/* y'a-t'il des arcs vers les noeuds prÃ©cÃ©dant w ? */
	/* calcule les arÃªtes sortantes de w et de leurs types */
	
	for(u=0;u<w;u++){ /* pour chaque noeud < w, v=type de l'arc */
	  if(u>=npath) v=0; /* test rapide: si u>=napth, alors clique */
	  else{ /* sinon calcule l'intersection des composantes */
	    v=SetCmp(c.comp[u],c.comp[w],c.tcomp[u],c.tcomp[w]);
	    /* ici les valeurs possibles pour SetCmp: 0,1,3,5,9 car T1 et T2 != {} */
	    if(v==1) v=-2; /* si v=1, T1 intersecte T2, et donc pas d'arÃªte u-w */
	    /* avant: v: 0=disjoint, 3=Ã©galitÃ©, 5=T1 sub T2, 9=T2 sub T1 */
	    v >>= 1; v -= (v==4); /* rem: si avant v=-2, alors aprÃ¨s v=-1 */
	    /* aprÃ¨s: v: 0=disjoint, 1=Ã©galitÃ©, 2=T1 sub T2, 3=T2 sub T1 */
	  } /* si v<0, alors pas d'arÃªte */
	  if(v>=0){ /* ajoute les arcs u->w et w->u */
	    ADD_ARC(c.G,u,(v<<LCONF)|w); /* u->w */
	    if(v>1) v ^= 1; /* si v=2, alors v->3, et si v=3, alors v->2 */
	    ADD_ARC(c.G,w,(v<<LCONF)|u); /* w->u (asymÃ©trie pour v=2,3) */
	  }
	} /* fin du for(u=...) */
	
	nextC:; /* prochaine composante (next i) */
	} /* fin du for(i=0...) */

      nextQ:; /* prochain chemin Q */
      }while(NextPath(G,Q,0));
      
      /* si ici goodxy=1, c'est que pour la paire xy, tous les chemins
	 entre xy sont parallÃ¨les Ã  P et qu'on a jamais vue de
	 composante Ã  FAUX. Dans ce cas on a trouvÃ© une bonne paire et
	 on a fini. */

      if(goodxy) goto fin_ps1; /* termine l'algo. */
      if(version==0) goto nextxy; /* si pas graphe des conflits */
      /* ici version>0 et goodxy=0 */
      
      /* Ici on a examinÃ© tous les chemins Q pour la paire xy et on
      crÃ©e au moins un noeud dans le graphe des conflits. Ici tous les
      noeuds crÃ©es ont comme valeur -1. Il reste Ã  essayer d'appliquer
      trois rÃ¨gles:

	 1) la rÃ¨gle du max: on Ã©crit 0 s'il existe une composante de
	 xy contenant toutes les autres (y compris celle Ã©valuÃ©es
	 rÃ©cursivement Ã  VRAI). Pour cela il faut avoir exactement une
	 seule composante maximale dans c.cmax qui doit Ãªtre Ã 
	 FAUX. Si elle est Ã  VRAI il faut supprimer la paire xy (car
	 la plus lourde est VRAI, donc toutes les autres sont lÃ©gÃ¨res
	 !) et passer Ã  la suivante.

	 2) la rÃ¨gle d'influence des voisins: vÃ©rifier si des voisins
	 v de chaque noeud u (de la paire xy) ne peuvent pas
	 influencer u. Par exemple, si la composante de v est la mÃªme
	 que celle de u il faut alors Ã©crire cette valeur dans le code
	 de u.

	 3) la rÃ¨gle du dernier -1: si le nb de valeurs < 1 est
	 exactement 1, alors il faut Ã©crire 0 dans ce noeud lÃ . Il
	 faut, comme pour la rÃ¨gle d'influence des voisins, balayer
	 tous les sommets u de la paire xy.

	 NB: on peut traiter dans un mÃªme parcours les deux derniÃ¨res
	 rÃ¨gles.
      */

      v=c.paire[c.G->n-1]; /* v=1er noeud de la derniÃ¨re paire */
      /* NB: c.outmem=0 et c.G->n > 0 */

      /* rÃ¨gle du max */
      if(c.ncmax==1){ /* sinon pas de rÃ¨gle du max */
	if(c.valmax){ /* supprime la paire xy */
	  ps1_delxy(&c,v);
	  goto nextxy;
	}
	/* applique la rÃ¨gle du max: on Ã©crit 0 dans la composante maximale */
	if(ps1_push(c.nodemax,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      }

      /* rÃ¨gle d'influence des voisins + rÃ¨gle du dernier -1 */
      for(u=v;u<c.G->n;u++){ /* pour tout noeud u de la derniÃ¨re paire */
	if(c.code[u]>=0) continue; /* plus rien Ã  faire */
	/* NB: si code[u]=0 ou 1, c'est qu'on a forcÃ©ment dÃ©jÃ  fait
	   un ps1_push() sur u. Et donc la rÃ¨gle d'influence des
	   voisins n'a pas a Ãªtre testÃ©e sur u */
	else /* code=-1 */
	  if(c.nbc[c.paire[v]]==1){ /* rÃ¨gle du dernier -1 */
	    /* tous les sommets sont Ã  1, sauf u qui est ici a -1 */
	    if(ps1_push(u,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	    else break; /* on peut sortir du for(u...), on a tout testÃ© */
	  }
	
	for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins v de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=i-Ã¨me voisin de u */
	  w=c.code[v]; if(w<0) continue; /* on ne dÃ©duit rien */
	  /* efface le code de v pour pouvoir le rÃ©Ã©crire avec
	     ps1_push().  Il faut prendre des prÃ©caution pour que c
	     soit cohÃ©rent (il faut mettre Ã  jour .nbi et .nbc). NB:
	     l'ancien code de v est w=0 ou 1, celui de u est aussi 0 ou 1 */
	  c.code[v]=-1;
	  c.nbi++; /* met Ã  jour le nombre d'indÃ©terminÃ©es */
	  c.nbc[c.paire[v]] += w; /* augmente s'il y avait 1 */
	  /* rÃ©Ã©crit le code pour v */
	  if(ps1_push(v,w,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	  /* NB: la propagation de cette Ã©criture va tester l'arc
	     v->u et donc modifier Ã©ventuellement le code de u */
	}
      }

    nextxy:;
    } /* fin du for(x,y=...) */
  
  /* Ici on a testÃ© toutes paires xy et aucune n'est bonne */
  /* NB: goodxy=0 */

  if((version>0)&&(!c.outmem)){

    /* Traitement final du graphe des conflits. On pourrait ici
       essayer d'Ã©liminer les derniers -1. Seulement, il semble
       qu'aucune des rÃ¨gles ne peut plus Ãªtre appliquÃ©es. */

  loop_ps1:
    do{
      if(c.nbi==0) break; /* fini: il ne reste plus de -1 Ã  modifier */
      x=c.nbi; /* mÃ©morise le nb total de valeurs Ã  -1 */
      //
      // que faire ? essayer toutes les valeurs possibles pour les -1 ?
      //
      // G->int1++;
    }
    while(c.nbi<x); /* tant qu'on a enlevÃ© au moins un -1, on recommence */
    
    /* ps1x */
    if(version==3){ /* NB: si appel rÃ©cursif, alors version!=3 */
      i=MEM(CPARAM,0,int);
      if(i>0){ /* lit les valeurs Ã  l'envers */
	u=MEM(CPARAM,(2*i-1)*sizeof(int),int);
	v=MEM(CPARAM,(2*i)*sizeof(int),int);
	MEM(CPARAM,0,int)--;
	if(ps1_push(u,v,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      goto loop_ps1;
      }
    }
    
    /*
      on cherche les noeuds "0" indÃ©pendants, correspondant Ã  des
      composantes lourdes (donc Ã  des inÃ©galitÃ©s).
      
      Algorithme: on ne sÃ©lectionne que les noeuds de code=0 et qui
      n'ont pas de voisins v de code aussi 0 avec v<u (inclus). Ils
      doivent donc Ãªtre de code 0 et minimal par rapport aux autres
      voisins de code 0. Si un tel noeud existe, on marque tous ces
      voisins de code 0 qui ne pourront plus Ãªtre sÃ©lectionnÃ©s.
    */

    c.nbzi=0; /* =nb de zÃ©ros indÃ©pendants */
    for(u=0;u<c.G->n;u++){ /* pour chaque noeud du graphe des conflits */
      if(c.code[u]) continue; /* saute le noeud si pas 0 */
      for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins u (qui est forcÃ©ment Ã  0) */
	v=c.G->L[u][i]; w=c.code[v&CONFMASK];
	if(!((w==0)||(w==CONFC))) continue; /* on ne regarde que les voisins de code=0 */
	if((v>>LCONF)==3) break; /* ici i<d */
      }
      if(i==d){ /* on a trouvÃ© un noeud u de code 0 et sans aucun arc u->v avec v<u */
	c.nbzi++; /* un de plus */
	for(i=0,d=c.G->d[u];i<d;i++){ /* modifie le code de tous les voisins Ã  0 de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=voisin de u */
	  if(c.code[v]==0) c.code[v]=CONFC; /* on ne modifie pas deux fois un voisin */
	}
      }else c.code[u]=CONFC; /* u n'est pas bon, on ne laisse pas la valeur 0 */
    }

  }
  /* ici on a fini le traitement du graphe des conflits */

 fin_ps1: /* termine l'algo avec goodxy=0 ou 1 */

  if(version>0){ /* affichage et libÃ©ration du graphe des conflits */
    if((!goodxy)&&(G==GF)) PrintConflit(&c); /* G=GF => affiche que le 1er niveau de rÃ©currence */
    /* efface la liste des composantes maximales */
    for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
    /* efface le graphe des conflits */
    for(u=0;u<c.G->n;u++) free(c.comp[u]);
    free_graph(c.G); /* c.G->L aprÃ¨s c.G->n n'ont pas Ã©tÃ© allouÃ©s ou sont dÃ©jÃ  libÃ©rÃ©s */
  }

  /* efface les tableaux allouÃ©s */
  free_path(Q);
  free_path(R);
  free_param_dfs(p);
  free(T);
  free(M);

  DEBUG(LEVEL--;);
  return goodxy;
}
#undef LCONF
#undef CONFMAX
#undef CONFMASK
#undef CONFC
#undef NCMAX
#undef PRINTS


/***********************************

       ROUTINES POUR LES
        ROUTING SCHEMES

***********************************/


int nca_bfs(int u,int v,const param_bfs *X){
/*
  Calcule le plus petit ancÃªtre commun entre u et v dans l'arbre donnÃ©
  par le bfs X. Si la forÃªt n'est pas connexe, alors on renvoie -1.

  Algorithme (de complexitÃ© dist(u,v)): on remonte d'abord le sommet
  le plus profond. Lorsque les deux sommets sont Ã  la mÃªme profondeur
  on teste si c'est les mÃªmes et si oui on a trouver le nca, ou alors
  on les remontes tous les deux.
*/
  if(X->D[u]<X->D[v]) SWAP(u,v);

  /* ici u est plus profond que v (ou Ã  mÃªme profondeur) */
  while(X->D[u]!=X->D[v]) u=X->P[u];

  /* ici u et v sont Ã  la mÃªme profondeur */
  while(u!=v) u=X->P[u],v=X->P[v];
	
  return u;
}


static inline int dist_nca(int const u,int const v,int const w,const int* const D){
/*
  Calcule la distance de u Ã  v dans un arbre via leur ancÃªtre commun
  w, en entrÃ© on doit avoir que D[u] est la distance de u Ã  la racine
  de l'arbre. Par rapport Ã  dist_bfs() il est nÃ©cessaire ici que w
  soit dÃ©terminÃ©, en particulier que u et v soit dans le mÃªme arbre.
*/
  return D[u] + D[v] - (D[w]<<1);
}


static inline int dist_bfs(int const u,int const v,const param_bfs *X){
/*
  Calcule la distance de u Ã  v dans l'arbre (en fait la forÃªt) donnÃ©
  par le bfs X. Si les sommets sont dans des composantes connexes
  diffÃ©rentes, alors on renvoie -1. C'est basÃ© sur le calcule du plus
  petit ancÃªtre commun.
*/
  int w=nca_bfs(u,v,X);
  if(w<0) return -1;
  return dist_nca(u,v,w,X->D);
}


string TopChrono(int const i){
/*
  Met Ã  jour le chronomÃ¨tre interne numÃ©ro i (i=0..CHRONOMAX-1) et
  renvoie sous forme de string (=char*) le temps Ã©coulÃ© depuis le
  dernier appel Ã  la fonction pour le mÃªme chronomÃ¨tre. La prÃ©cision
  dÃ©pend du temps mesurÃ©. Il varie entre la seconde et le 1/1000 de
  seconde.  Plus prÃ©cisÃ©ment le format est le suivant:

  1d00h00'  si le temps est > 24h (prÃ©cision: 1')
  1h00'00"  si le temps est > 60' (prÃ©cision: 1s)
  1'00"0    si le temps est > 1'  (prÃ©cision: 1/10s)
  1"00      si le temps est > 1"  (prÃ©cision: 1/100s)
  0"000     si le temps est < 1"  (prÃ©cision: 1/1000s)

  Pour initialiser et mettre Ã  jour tous les chronomÃ¨tres (dont le
  nombre vaut CHRONOMAX), il suffit d'appeler une fois la fonction,
  par exemple avec TopChrono(0). Si i<0, alors les pointeurs allouÃ©s
  par l'initialisation sont dÃ©sallouÃ©s. La durÃ©e maximale est limitÃ©e
  Ã  100 jours. Si une erreur se produit (durÃ©e supÃ©rieure ou erreur
  avec gettimeofday()), alors on renvoie la chaÃ®ne "--error--".
  
  Exemple d'usage:

    TopChrono(0); // initialise le chrono numÃ©ro 0
    x=f(y);       // appel Ã  une fonction
    printf("time=%s\n",TopChrono(0)); // affiche le temps Ã©coulÃ©
    TopChrono(-1); // libÃ¨re les pointeurs allouÃ©s

  On n'utilise pas clock() qui se limite Ã  72 minutes. Pour des durÃ©es
  en micro secondes mais plus beaucoup plus longues, il faut utiliser
  gettimeofday() comme ceci:

  gettimeofday(&t0,NULL);
  ...
  gettimeofday(&t1,NULL);
  long long t = (t1.tv_sec-t0.tv_sec)*1000000LL + t1.tv_usec-t0.tv_usec;
*/
  if(i>=CHRONOMAX) Erreur(26);
  
  /* variables globales, locale Ã  la fonction */
  static int first=1; /* =1 ssi c'est la 1Ã¨re fois qu'on exÃ©cute la fonction */
  static string str[CHRONOMAX];
  static struct timeval last[CHRONOMAX],tv;
  int j;

  if(i<0){ /* libÃ¨re les pointeurs */
    if(!first) /* on a dÃ©jÃ  allouÃ© les chronomÃ¨tres */
      for(j=0;j<CHRONOMAX;j++)
	free(str[j]);
    first=1;
    return NULL;
  }

  /* tv=temps courant */
  int err=gettimeofday(&tv,NULL);

  if(first){ /* premiÃ¨re fois, on alloue puis on renvoie TopChrono(i) */
    first=0;
    for(j=0;j<CHRONOMAX;j++){
      str[j]=malloc(10); // assez grand pour "--error--", "99d99h99'" ou "23h59'59""
      last[j]=tv;
    }
  }

  /* t=temps en 1/1000" Ã©coulÃ© depuis le dernier appel Ã  TopChrono(i) */
  long t=(tv.tv_sec-last[i].tv_sec)*1000L + (tv.tv_usec-last[i].tv_usec)/1000L;
  last[i]=tv; /* met Ã  jour le chrono interne i */
  if((t<0L)||(err)) t=LONG_MAX; /* temps erronÃ© */
  
  /* Ã©crit le rÃ©sultat dans str[i] */
  for(;;){ /* pour faire un break */
    /* ici t est en milliÃ¨me de seconde */
    if(t<1000L){ /* t<1" */
      sprintf(str[i],"0\"%03li",t);
      break;
    }
    t /= 10L; /* t en centiÃ¨me de seconde */
    if(t<6000L){ /* t<60" */
      sprintf(str[i],"%li\"%02li",t/100L,t%100L);
      break;
    }
    t /= 10L; /* t en dixiÃ¨me de seconde */
    if(t<36000L){ /* t<1h */
      sprintf(str[i],"%li'%02li\"%li",t/360L,(t/10L)%60L,t%10L);
      break;
    }
    t /= 10L; /* t en seconde */
    if(t<86400L){ /* t<24h */
      sprintf(str[i],"%lih%02li'%02li\"",t/3600L,(t/60L)%60L,t%60L);
      break;
    }
    t /= 60L; /* t en minute */
    if(t<144000){ /* t<100 jours */
      sprintf(str[i],"%lid%02lih%02li'",t/1440L,(t/60L)%24L,t%60L);
      break;
    }
    /* error ... */
    sprintf(str[i],"--error--");
  }
  
  return str[i];
}


/* type pour une fonction f(u,v,T) renvoyant la longueur d'une route
   de u vers v Ã©tant donnÃ©es la table de routage globale T */
typedef int(*rt_length)(int,int,void*);


/* structure de donnÃ©es pour la table de routage d'un sommet */
typedef struct{
  int n;      /* taille des listes */
  int *node;  /* liste de sommets */
  int *dist;  /* distances */
  int radius; /* rayon d'une boule */
  int vpd;    /* voisin par dÃ©faut */
} table;


/*

pour avoir plusieurs noms diffÃ©rents pour une mÃªme structre

typedef union{
  struct{
    int n;
    int *dist;
  };
  struct{
    int n;
    int *color;
  };
} tableau;

  tableau *T=malloc(sizeof(*T));
  T->n=0;
  T->dist=NULL;
  printf("->y=%p\n",T->color);

*/

table *new_table(int const n){
/*
  CrÃ©e un objet de type table oÃ¹ les champs sont initialisÃ©s Ã  leur
  valeurs par dÃ©faut, les pointeurs Ã©tant initialisÃ©s Ã  NULL. Si n>0,
  alors les pointeurs (->node et ->dist) de taille n sont allouÃ©s, et
  le champs ->n est initialisÃ© Ã  n.
*/
  NALLOC(table,X,1);
  X->n=0;
  X->node=NULL;
  X->dist=NULL;
  X->radius=-1;
  X->vpd=-1;

  if(n>0){
    X->n=n;
    ALLOC(X->node,n);
    ALLOC(X->dist,n);
  }
  
  return X;
}


void free_table(table *X){
  if(X==NULL) return;
  free(X->node);
  free(X->dist);
  free(X);
  return;
}


/* tables de routage globale pour rs_cluster() */
typedef struct{
  int n; // pour libÃ©rer les n tables
  table **B; // les boules
  table **S; // les spanning trees depuis les landmarks
  table **R; // les sommets en charge des landmarks
  table **W; // tables des voisins dans le cluster
  int *H; // hahs des sommets
  int *C; // pour routage vers voisins de couleur i
  int center; // centre du cluster (pour VARIANT=1)
} rs_cluster_tables;


void free_rs_cluster_tables(rs_cluster_tables *X){
/*
  LibÃ¨re les tables crÃ©Ã©es par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->S[u]);
    free_table(X->R[u]);
    free_table(X->W[u]);
  }
  free(X->B);
  free(X->S);
  free(X->R);
  free(X->W);
  free(X->H);
  free(X->C);
  free(X);
  return;
}


/* tables de routage globale pour rs_dcr() */

typedef struct{ // structure "boule-contiguÃ«" spÃ©cifique au schÃ©ma AGMNT
  int d;    // distance entre u et v
  int v;    // sommet v=CONT[u][i]
  int s;    // landmark, <0 pour dire via boule-contiguÃ«
  int w;    // nca entre u et v (si via landmark)
} contigue;

typedef struct{
  table **B; // B[u]=tables de voisinage indexÃ©es par les sommets
  table **W; // W[u]=tables des reprÃ©sentants indexÃ©es par les couleurs 
  param_bfs **S; // S[u]=bfs(u) pour u landmark (sinon =NULL)
  int *H; // H[u]=hash du sommet u
  int *C; // C[u]=couleur du sommet u
  int **dist; // distances partielles (issues des landmarks)
  int n; // pour libÃ©rer les n tables
  // pour AGMNT (=NULL sinon)
  contigue **CONT;
  int *F; // F[C[u]]=|CONT[u]|=#sommets de hash C[u]
} rs_dcr_tables;


int fcmp_contigue(const void *P,const void *Q)
/*
  Compare les champs .v d'une structure "contigue" (pour AGMNT).
*/
{
  int const p=(*(contigue*)P).v;
  int const q=(*(contigue*)Q).v;
  return (p>q) - (p<q);
}


void free_rs_dcr_tables(rs_dcr_tables *X){
/*
  LibÃ¨re les tables crÃ©Ã©es par rs_dcr().  On a pas besoin de libÃ©rer
  les pointeurs de X->dist, car ce sont ceux de X->S->D.
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->W[u]);
    free_param_bfs(X->S[u]); // libÃ¨re aussi X->dist[u]
  }
  free(X->B);
  free(X->W);
  free(X->S);
  free(X->H);
  free(X->C);
  free(X->dist);
  FREE2(X->CONT,X->n); // pour AGMNT
  free(X->F); // pour AGMNT
  free(X);
  return;
}


/* tables de routage pour rs_tzrplg() */
typedef struct{
  table **B; // boules
  table **L; // landmarks
  int* label; // les Ã©tiquettes des sommets
  int n; // pour libÃ©rer les n tables
} rs_tzrplg_tables;


void free_rs_tzrplg_tables(rs_tzrplg_tables *X){
/*
  LibÃ¨re les tables crÃ©ees par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->L[u]);
  }
  free(X->B);
  free(X->L);
  free(X);
  return;
}


/* tables de routage pour rs_bc() */
typedef struct{
  int nbfs; // taille de L
  param_bfs **L; // liste des BFS;
  param_bfs **Lu; // liste des BFS indexÃ© par les sommets
  int **dist; // distances partielles (issues des landmarks)
} rs_bc_tables;


void free_rs_bc_tables(rs_bc_tables *X){
/*
  LibÃ¨re les tables crÃ©ees par rs_bc().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->nbfs;u++)
    free_param_bfs(X->L[u]); // libÃ¨re X->Lu et aussi X->dist

  free(X->L);
  free(X->Lu);
  free(X->dist);

  free(X);
  return;
}


/* tables de routage pour rs_hdlbr() */
typedef struct{
  table **B; // boules
  table **L; // landmarks
  int* Core; // liste des landmarks
  int* H;    // hash des sommets
  int n;     // pour libÃ©rer les n tables
  int core_size;
} rs_hdlbr_tables;


void free_rs_hdlbr_tables(rs_hdlbr_tables *X){
/*
  LibÃ¨re les tables crÃ©ees par rs_hdlbr().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->L[u]);
  }
  free(X->B);
  free(X->L);
  free(X->Core);
  free(X->H);

  free(X);
  return;
}


/* Donne l'Ã©cart type calculÃ© Ã  partir de la somme s1 et de la somme
   s2 du carrÃ© de n valeurs. Si n<=0, alors on retourne NaN. On
   utilise la formule Ã©cart_type(X) := sqrt(E(X^2)-E(X)^2). Attention
   d'utiliser (s1/n)*(s1/n) plutÃ´t que (s1*s1)/(n*n) Ã  cause du
   dÃ©passement arithmÃ©tique du type int de n*n. NB: nan(NULL) provoque
   un warning sur certains systÃ¨mes.
*/
#define ECARTYPE(s1,s2,n) ((n>0)?sqrt(s2/(double)n-(s1/(double)n)*(s1/(double)n)):nan("NaN"))

/* Affiche le min/max, la moyenne et l'Ã©cart type d'un tableau
   d'entiers vÃ©rifiant une certaine condition. Si la condition n'est
   jamais vÃ©rifiÃ©e (moyenne sur aucun terme), alors les min/max et la
   moyenne sont Ã©valuÃ©s Ã  0, et on ajoute le texte "(undefined)". Le
   term et la condition peuvent dÃ©pendre d'un indice _i. Si str
   commence par "- ", alors les tirets "- " ne sont pas affichÃ©s (pour
   Ãªtre cohÃ©rent avec PrintDistribution()).

   Ex:  MINMAXMOY(T[_i],n,T[_i]>0,"size of T");

   Calcule la moyenne des termes T[_i] tq T[_i]>0 pour i=0..n-1, le
   texte sert pour l'affichage.
*/
#define MINMAXMOY(term,n,condition,str)					\
  do{									\
    int _i,_m0,_m1,_cpt,_t,_b;long _sum=0L,_sum2=0L;string _r="- ";	\
    for(_i=_cpt=0;_i<(n);_i++)						\
      if(condition){							\
        _t=term; _sum += _t; _sum2 += _t*_t;				\
	if(_cpt) _m0=min(_m0,_t), _m1=max(_m1,_t);			\
	else _m0=_m1=_t;						\
	_cpt++;								\
      }									\
    if(_cpt==0) _m0=_m1=0;						\
    double _avg=_cpt?(double)_sum/(double)_cpt:0;			\
    if(strncmp(str,_r,2)) _b=0; else _b=2,_r="";			\
    printf("%sminimum %s: %i%s\n",_r,(string)str+_b,_m0,_cpt?"":" (undefined)"); \
    printf("%smaximum %s: %i%s\n",_r,(string)str+_b,_m1,_cpt?"":" (undefined)"); \
    printf("%saverage %s: %.2lf",_r,(string)str+_b,_avg);		\
    if(_cpt) printf(" Â± %.2lf (%li/%i)\n",ECARTYPE(_sum,_sum2,_cpt),_sum,_cpt);	\
    else printf(" (undefined)\n");					\
  }while(0)


/* Calcule puis affiche la frÃ©quence, le min et le max d'un tableau
   d'entiers.

   Ex: FREQMINMAX(F,k,T,n,"size of T");

   Calcule dans F la frÃ©quence des Ã©lÃ©ments du tableau T Ã  n Ã©lÃ©ments,
   les valeurs de T Ã©tant des entiers de [0,k[. Attention ! le tableau
   F est allouÃ© et initialisÃ© par la macro, et c'est la macro qui
   dÃ©termine k qui doit Ãªtre une variable. Le texte sert pour
   l'affichage. S'il commence par "- ", alors le tiret "- " n'est pas
   affichÃ© (pour Ãªtre cohÃ©rent avec MINMAXMOY() et
   PrintDistribution()).
*/
#define FREQMINMAX(F,k,T,n,str)						\
  do{									\
    int _u,_v,_i,_b;string _r="- ";					\
    F=SortInt(T,NULL,(n),0,&k,SORT_FREQv);				\
    for(_u=_v=F[_i=0];_i<k;_i++) _u=min(_u,F[_i]), _v=max(_v,F[_i]);	\
    if(strncmp(str,_r,2)) _b=0; else _b=2,_r="";			\
    printf("%sbalance ratio of the frequency %s: %.0lf%% (ceil{n/k}=%i max=%i min=%i)\n", \
	   _r,(string)str+_b,100.0*_u/(double)_v,iceil((n),k),_v,_u);	\
  }while(0)


void ruling(string const s,int const n)
/*
  Affiche n fois Ã  la suite la chaÃ®ne de caractÃ¨res s.
*/
{
  for(int i=0;i<n;i++) printf("%s",s);
}


#define STR_DISTRIB "â–ª" /* caractÃ¨re pour affichage des distributions */
#define LEN_DISTRIB (60.0) /* longueur max pour l'affichage des distributions */
#define RULING(x) ruling(STR_DISTRIB,min((int)(LEN_DISTRIB*x),(int)LEN_DISTRIB-5))


void PrintDistribution(int const* const Z,int const n,int r,
		       string const message)
/*
  Affiche la distribution des valeurs entiÃ¨res contenues dans le
  tableau Z (de taille n>0) selon r>0 intervalles et avec le message
  m. Les intervalles consÃ©cutifs nuls sont regroupÃ©s. Donc le nombre
  d'intervalles affichÃ©s peut Ãªtre < r. Si r<0, alors on affiche
  toutes les valeurs (de frÃ©quence non nulle) plutÃ´t que des
  intervalles. En plus de la distribution, on affiche la valeur min,
  max et la moyenne avec l'Ã©cart type, sauf si r=-2. Si message
  commence par "- " alors le tiret n'est pas affichÃ©.

  Exemple avec r=10 et m="routing table size":

  - routing table size distribution: 10 ranges
    [339,381[ 	01% â–ª [Ã— 168] 
    [381,423[ 	08% â–ªâ–ªâ–ªâ–ª [Ã— 698] 
    [423,465[ 	15% â–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ª [Ã— 1282] 
    [465,507[ 	57% â–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ªâ–ª [Ã— 4848] 
    [507,549[ 	11% â–ªâ–ªâ–ªâ–ªâ–ªâ–ª [Ã— 975] 
    [549,591[ 	02% â–ª [Ã— 230] 
    [591,633[ 	02% â–ª [Ã— 182] 
    [633,675[ 	00%  [Ã— 18] 
    [675,717[ 	00%  [Ã— 4] 
    [717,752[ 	00%  [Ã— 1] 
  - minimum routing table size: 339
  - maximum routing table size: 751
  - average routing table size: 481.48 Â± 13.84 (4047360/8406)
*/
{
  if((n<=0)||(Z==NULL)){ printf("- empty distribution\n"); return; }
  int i,j,k,m0,m1; /* calcule m0=min_i{Z[i]} et m1=max_i{Z[i]} */
  int const b0=(r==-2);
  int const b1=(r<0);
  if(b1) r=n;

  m0=m1=Z[0]; 
  for(i=1;i<n;i++) m0=min(m0,Z[i]),m1=max(m1,Z[i]);
  int d=(++m1)-m0; /* d=|[m0,m1[| = plage des valeurs de Z, d>=1 */
  r=min(r,d); /* r=nombre d'intervalles maximum, ne peut dÃ©passer d */
  d=iceil(d,r); /* d=ceil(d/r)=taille des intervalles sauf le dernier
		   qui est plus court. NB: d>=1 car d>=r>=1 */

  /* calcule la distribution F de Z */

  NALLOCZ(int,F,r,0); /* initialise F */
  for(i=0;i<n;i++) F[(Z[i]-m0)/d]++; /* compte les valeurs de chaque intervalle */
  for(i=k=0;i<r;k++,i=j){ /* k=nombre d'intervalles (ou valeurs) qui vont Ãªtre affichÃ©s */
    j=i+1;
    if(F[i]==0){
      while((j<r)&&(F[j]==0)) j++; /* j=prochain intervalle */
      if(j==r) break;
    }
  }
  int b3=2;
  string tiret="- ";
  if(strncmp(message,tiret,2)) b3=0; else tiret=""; // *tiret=0 -> Bus error: 10
  printf("%s%s distribution: %i %s%s\n",
	 tiret,message+b3,k,b1?"value":"range",PLURIEL(k));

  /* affichage des intervalles (ou valeurs) de F, Ã©ventuellement en
     fusionnant les intervalles vides consÃ©cutifs */

  double x;
  for(i=0;i<r;){
    j=i+1;
    if(F[i]==0){
      while((j<r)&&(F[j]==0)) j++; /* j=prochain intervalle */
      if(j==r) break;
    }
    x=(double)F[i]/(double)n;
    if(b1){ if(F[i]) printf("  %i: \t",m0+i*d); }
    else printf("  [%i,%i[ \t",m0+i*d,min(m0+j*d,m1));
    printf("%02i%% ",(int)(100*x));
    RULING(x);
    printf(" [Ã— %i]\n",F[i]);
    i=j; /* ici j>=i+1 */
  }
  free(F); /* ne sert plus Ã  rien */
  if(b0) return;

  /* affichage min/max/moy Ã©cart type pour Z */
  MINMAXMOY(Z[_i],n,1,message);
  return;
}


string millier(long const i){
/*
  Renvoie une chaÃ®ne de caractÃ¨res correspondant Ã  l'entier long i
  Ã©crit avec le sÃ©parateur des milliers ','. Ne pas faire de free()
  sur le pointeur renvoyÃ©.

  Ex: > printf("n=%s\n",millier(-123456789));
      > n=-123,456,789
*/  
  static char r[64];
  string p,s=r; /* r=rÃ©sultat final, s=r ou r+1 si i<0 */
  int n;
  
  snprintf(r,sizeof(r),"%li",i); /* Ã©crit i dans r */
  if(i<0) s++; /* si i<0, r="-123..." */
  n=strlen(s); /* n=longueur(r) sans le signe */
  p=s+n-3;
  
  while(p>s){
    memmove(p+1,p,n);
    *p=',';
    n+=4; /* 4 chiffres avec la virgule */
    p-=3; /* paquet de 3 chiffres */
  }
  
  return r;
}


static inline int lg(unsigned long n){
#if defined(__APPLE__) && defined(__MACH__) // Apple OSX and iOS (Darwin)
  return flsl(n);
#else
  return 64-__builtin_clzl(n)-(n==0);
#endif
}
/*
  Renvoie le plus petit entier k tel que n < 2^k. Autrement dit k est
  tel que 2^{k-1} <= n < 2^k. D'oÃ¹:

             / 1+floor{log_2(n)} si n>0     n  â”‚ 0 1 2 3 4 5 6 7 8 9 ...
   lg(n) := {                             â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
             \ 0 sinon (n=0)              lg(n)â”‚ 0 1 2 2 3 3 3 3 4 4 ...

  Lorsque n>0, c'est aussi le nombre de bits dans l'Ã©criture binaire
  de n, et c'est Ã©galement un de plus que la position du bit de poids
  fort.

  NB: on peut calculer lg(lg(...(lg(n))...)) sans erreur, car lg(n)â‰¥0
  est dÃ©fini pour tout entier nâ‰¥0. Il s'agit aussi de la fonction
  fls(n) que l'on peut aussi implÃ©menter comme ceci:

       int k=0;
       while(n>0) n>>=1,k++;
       return k;

  On peut aussi utiliser __builtin_clzl(n) qui est une macro qui
  compte le nombre de zÃ©ros Ã  gauche (Count Leading Zeros) pour un
  long n. Mais cela donne 63 pour n=0 (au lieu de 64).
*/


/*
  Calcule rapide d'un hash sur 32-bits, d'aprÃ¨s Bob Jenkins (2006),
  disponible ici: http://burtleburtle.net/bob/c/lookup3.c
*/
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
#define mix(a,b,c)		\
  do{				\
    c ^= b; c -= rot(b,14);	\
    a ^= c; a -= rot(c,11);	\
    b ^= a; b -= rot(a,25);	\
    c ^= b; c -= rot(b,16);	\
    a ^= c; a -= rot(c,4);	\
    b ^= a; b -= rot(a,14);	\
    c ^= b; c -= rot(b,24);	\
  }while(0)


int hash_mix(int const u)
/*
  Renvoie un entier h(u) positif (sur 31 bits donc) basÃ© sur le
  mÃ©lange de deux constantes entiÃ¨res alÃ©atoires A et B de 32 bits. Si
  u<0, alors les constantes A,B sont initialisÃ©es et on renvoie -1.
*/
{
  static unsigned A,B;

  if(u<0){
    A=random();
    B=random();
    return -1;
  }

  unsigned a=A,b=B,c=u; // NB: a et b sont modifiÃ©es
  mix(a,b,c);

  return c>>1; // valeurs sur 31 bits, toujours positive
}


int hash_prime(int const u)
/*
  Renvoie un entier h(u) de [0,p[ oÃ¹ p=2^31-1=0x7FFFFFFF est un nombre
  premier (on pourrait aussi utiliser 2^61-1 qui est premier). Il faut
  0<=u<=p. Si u<0, alors les constantes A,B (qui accÃ©lÃ¨rent les
  calculs) sont initialisÃ©es et on renvoie -1.  On prend la fonction
  de hashage de Carter & Wegman (1978) avec h(u)=(A*u+B) mod p (puis
  h(u)%k) oÃ¹ A,B sont dans [0,p[ et A impair. Il est important de
  faire les calculs en "unsigned long" car la valeur A*x+B dÃ©passe 32
  bits en gÃ©nÃ©ral. Sinon, le calcul de (A*x+B)%p est erronÃ© (s'il est
  rÃ©alisÃ© en 32 bits).
*/
{
  static long unsigned const p=0x7FFFFFFF;
  static long unsigned A,B;
  long unsigned const x=u;

  if(u<0){
    A=1L+randomu(p);
    A|=1L; // force A Ã  Ãªtre impair
    B=randomu(p);
    return -1;
  }

  return (A*x+B)%p; // en principe, le modulo ne sert Ã  rien
}


int hash_shuffle(int u,int const n)
/*
  Calcule un entier p(u) de [0,n[ qui est une permutation
  p:[0,n[->[0,n[ basÃ©e sur deux constantes alÃ©atoires R1,R2 de
  [0,n[. Le temps de calcul est d'environ (logn)/2. Si u<0, alors les
  constantes R1,R2,K,N0,N1 (qui accÃ©lÃ¨rent les calculs) sont
  initialisÃ©es et on renvoie -1. Cette permutation est basÃ©e sur deux
  "shuffles" qui sont les permutations P0 et P1 suivantes:

    P0(u)=(u>>1)+((u&1)==0)*N0; // avec N0=floor(n/2)
    P1(u)=(u>>1)+((u&1)==1)*N1; // avec N1=ceil(n/2)

  Pour n=10, la permutation P0(u) donne: 0-5-1-6-2-7-3-8-4-9. Donc on
  Ã©crit 0,1,2,... aux positions paires 0,2,4,6,...

  Pour n=10, la permutation P1(u) donne: 5-0-6-1-7-2-8-3-9-4. Donc on
  Ã©crit 0,1,2,... aux positions impaires 1,3,5,6,...

  On applique la permutation P0 ou P1 suivant les bits de poids faible
  de R1, en rÃ©pÃ©tant cette opÃ©ration K fois, oÃ¹ K=(logn)/2 est environ
  la moitiÃ© des bits nÃ©cessaires pour Ã©crire n en binaire. En effet,
  on retrouve la valeur de dÃ©part aprÃ¨s avoir effectuÃ© logn fois la
  permutation P0 ou logn fois la permutation P1.

  Enfin, on inverse le rÃ©sultat (u->n-u) et on le dÃ©cale de R2
  position modulo n (u->(u+R2)%n). Ainsi la permutation finale de u,
  p(u), peut valoir n'importe quelle valeur entre [0,n[.
*/
{
  /* variables globales / locales */
  static int R1=0;
  static int R2=0;
  static int N0=0;
  static int N1=0;
  static int K=0;

  if(u<0){
    R1=randomu(n); // R1 contient au moins k bits alÃ©atoires
    R2=randomu(n); // pour le dÃ©calage final
    N0=N1=(n>>1);  // N0=floor(n/2)
    N1 += (n&1);   // N1=ceil(n/2)
    K=lg(n)/2;     // K=la moitiÃ© de floor(logn)
    return -1;
  }
  
  /* calcule u=p(u) */
  int i,b,r=R1;
  for(i=0;i<K;i++){
    b=(r&1); // lit un bit b de R1 et applique P0 ou P1
    u=(u>>1)+((u&1)==b)*(b?N1:N0); // u/2, u/2+N0 ou u/2+N1
    r>>=1; // supprime le bit lu
  }

  /* inverse et dÃ©calage alÃ©atoire */
  return (n+R2-u)%n;
}


int *MakeHash(int *H,int const n,int const k,int const M)
/*
  Calcule un hash h:[0,n[->[0,k[ avec k<=n selon la mÃ©thode M (voir la
  variable HASH). Remplit et renvoie le tableau H de taille n (ou bien
  un nouveau tableau allouÃ© si H=NULL). Si M n'est pas une valeur
  reconnue, on fait comme si M=H_MOD. On renvoie NULL si n<=0. En
  gÃ©nÃ©ral, il faut Ã©viter de faire un modulo (%) qui est jusqu'Ã  200
  fois plus lent qu'une addition entiÃ¨re.

  Il y a d'autres fonctions de hash trÃ¨s intÃ©ressantes Ã  voir (comme
  http://burtleburtle.net/bob/hash/doobs.html). On y trouve en
  particulier, lookup3.c avec mix(a,b,c) de Bob Jenkins 2006
  http://burtleburtle.net/bob/c/lookup3.c. Il y a aussi MurmurHash2
  (2010).
*/
{
  if(n<=0) return NULL;
  if(H==NULL) ALLOC(H,n);
  
  int u,a;
  
  switch(M){
    
  case H_PRIME:
    hash_prime(-1); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_prime(u)%k;
    break;
    
  case H_SHUFFLE:
    hash_shuffle(-1,n); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_shuffle(u,n)%k;
    break;
    
  case H_MIX:
    hash_mix(-1); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_mix(u)%k;
    break;
    
  case H_MOD:
  default:
    a=randomu(k); // attention que u+a reste >=0
    for(u=0;u<n;u++) H[u]=(u+a)%k;
    break;
  }
  
  return H;
}


rs_cluster_tables *rs_cluster(graph* const G,int k)
/*
  Calcule pour le graphe G le routing scheme cluster de paramÃ¨tre
  k>0. On renvoie les tables calculÃ©es. Le stretch est toujours <=
  5. Cependant lorsque k=1 le stretch devient <= 3.

  Si bit-0 de VARIANT=1, alors les tables B sont vides.
  Si bit-1 de VARIANT=1, alors les tables W sont vides.
  Si bit-2 de VARIANT=1, alors les tables B des voisins du cluster sont vides.
  Si bit-3 de VARIANT=1, alors les tables B des voisins du cluster n'ont pas de sommets du cluster.
*/
{
  int const n=G->n;
  param_bfs *X0,*X;
  int u,v,i,j,r,t,center;
  table **S,**B,**W,**R;
  int *Z,*F;

  printf("\nCLUSTER\n");
  BARRE;
  
  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  /* trouve un sommet center de degrÃ© max, v=deg(center) */
  for(u=v=center=0;u<n;u++)
    if(G->d[u]>v) v=G->d[center=u];

  if((VARIANT&2)&&(k<v+1)) Erreur(6); // paramÃ¨tres incohÃ©rents
  k=min(k,v+1); /* k=taille effective du cluster=#colors */
  printf("- cluster size: %i",k);
  printf(" (ceil{âˆšn}=%i)\n",(int)ceil(sqrt(n)));
  printf("- degree of the center: %i (id:%i)\n",v,center);

  /* construit le cluster C de taille k Ã  partir des k-1 sommets de
     plus haut degrÃ© voisins du center */

  /* construit une table H des degrÃ©s des voisins de center */
  /* on en profite pour calculer le degrÃ© min pour optimiser SortInt() */
  u=v; /* u=futur degrÃ© min, v=deg(center) */
  NALLOCZ(int,H,v,G->d[G->L[center][_i]]); /* H[i]=degrÃ© du i-Ã¨me voisin du center */
  for(i=0;i<v;i++) u=min(u,H[i]); /* u=voisin de center de deg min */

  /* trie selon les listes des degrÃ©s H: Z[i]=indice dans H du degrÃ© de rank i */
  r=v-u+1; /* H[i] valeur dans [u,v]=[u,u+r[ */ 
  Z=SortInt(H,NULL,v,u,&r,SORT_INDEXi); /* H[Z[v-1]]=+grand degrÃ© parmi les v voisins du center */
  free(H); /* H ne sert plus Ã  rien, seul Z sert encore */

  /* remplit C selon les degrÃ©s dÃ©croissant */
  NALLOC(int,C,k);
  C[0]=center; /* met le center dans C, forcÃ©ment de plus haut degrÃ© */
  for(i=1;i<k;i++) C[i]=G->L[center][Z[--v]]; /* les k-1 voisins de + haut degrÃ© */
  free(Z); /* Z ne sert plus Ã  rien */
  
  /* affiche un aperÃ§u des degrÃ©s du cluster */
  printf("- cluster degree: ");
  APERCU(G->d[C[_i]],k,10,2);

  /* on trie C pour dÃ©terminer la couleur de u dans C rapidement. La
     couleur de u dans C est alors i=SetSort(C,k,1). Attention ! la
     couleur du center n'est pas forcÃ©ment 0. */
  QSORT(C,k,fcmp_int);
  
  X0=bfs(G,center,NULL); /* calcule un BFS depuis center */
  printf("- eccentricity of the center: %i\n",X0->radius);
  printf("- time to construct the cluster with its BFS: %s\n",TopChrono(1));

  /*
    Calcule, pour chaque sommet u de C, le nombre de sommets qui ne
    sont pas dans C et qui ont u comme plus proche ancÃªtre. On calcule
    ensuite le maximum de ces nombres.

    On utilise un simple tableau H indiquant si un sommet a dÃ©jÃ  Ã©tÃ©
    visitÃ© ou pas. De plus, si c'est un sommet de C, H[u] indique
    combien de sommets lui sont directement descendant. Au dÃ©part
    H[u]=1 si u est dans C, et H[u]=-1 sinon. Puis, pour chaque sommet
    u qui n'est pas dÃ©jÃ  visitÃ© ou dans C, on marque u et on remonte
    dans l'arbre u=parent(u) jusqu'Ã  trouver soit un sommet de l'arbre
    dÃ©jÃ  visitÃ© ou bien un sommet de C. Pour marquÃ© un sommet u que
    l'on visite, on pose H[u]=0. Si on tombe sur un sommet du cluster,
    on incrÃ©mente H[u]. Donc H[u]>0 si u est dans C, H[u]<0 si u n'a
    pas Ã©tÃ© visitÃ©, et H[u]==0 sinon (pas dans C mais dÃ©jÃ 
    visitÃ©). Cet algo prend un temps O(n) car chaque arÃªte de l'arbre
    n'est visitÃ©e qu'une seule fois.
   */
  ALLOCZ(H,n,-1);
  for(i=0;i<k;H[C[i++]]=1);
  for(u=0;u<n;u++){
    v=u;
    while(H[v]<0){ /* v n'est pas marquÃ© */
      H[v]=0; /* on marque v */
      v=X0->P[v]; /* v=parent(v) */
    }
    if(H[v]>0) H[v]++;
    }
  j=0; /* cherche u de C avec H[u] maximum */
  for(i=1;i<k;i++) if(H[C[i]]>H[C[j]]) j=i;
  printf("- maximum number of cluster direct descendants: %i (id: %i)\n",H[C[j]],C[j]);
  free(H);

  /* efface les pÃ¨res des sommets de C pour test d'appartenance rapide */
  for(i=0;i<k;i++) X0->P[C[i]]=-2; /* si -2 alors dans C */

  DEBUG(PRINTT(C,k););

  /****************************/
  /* calcule des stats / hash */
  /****************************/

  /* calcule le nombre d'arÃªtes intra et extra C, de voisins de C */
  TopChrono(1);
  ALLOCZ(H,n,0);
  t=r=0;
  for(i=0;i<k;i++){
    u=C[i]; // u=i-Ã¨me sommet de C
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-Ã¨me voisin de u */
      if(X0->P[v]==-2) t++; // arÃªte dans C
      else{ r++; H[v]=1; } // arÃªte hors C; marque les sommets
    }
  }
  printf("- #edges inside the cluster: %i\n",t/2);
  printf("- #edges outgoing the cluster: %i\n",r);
  for(u=t=0;u<n;u++) t+=H[u];
  printf("- #neighbors of the cluster: %i\n",t);
  printf("- their average #neighbors in the cluster: %.2lf (%i/%i)\n",(double)r/(double)t,r,t);
  
  /* construit le hash H[u]=0..k-1 de chaque sommet u=0..n-1 */
  /* le sommet C[i] est en charge des sommets de hash i=0..k-1 */

  MakeHash(H,n,k,HASH);
  FREQMINMAX(F,k,H,n,"the hash");
  DEBUG(PRINTT(H,n);PRINTT(F,k););
  printf("- time to compute stats and hash: %s\n",TopChrono(1));
  BARRE;
  
  /************************/
  /* calcule les tables B */
  /************************/
  
  /*
    Tables seulement dÃ©finies pour les sommets u qui ne sont pas dans
    C. Il s'agit dans un premier temps des sommets Ã  distance < r de u
    oÃ¹ r est la distance du plus proche sommet de C. Cependant, si le
    plus proche sommet est center (cela ne peut pas se produire si C
    contient le centre et tous ces voisins, en particulier si k=n), il
    convient alors de vÃ©rifier si Ã  distance r-1 de u il n'y a pas un
    sommet de C. Si oui, il faut alors prendre r-1 et non
    r. Eventuellement tree(center) peut-Ãªtre modifiÃ©. Le sommet u est
    supprimÃ© de B[u]. Dans un second temps on ajoute tous les voisins
    de u Ã  B[u] si |B[u]| Ã©tait vide, c'est-Ã -dire si u Ã©tait voisin
    de C. Le choix par dÃ©faut (->vpd) est le voisin de u allant vers
    le sommet de B[u] le plus proche du center (dans tree(center)).
  */

  TopChrono(1);
  ALLOCZ(B,n,NULL); /* tableau de n tables B (vides au dÃ©part) */
  X=new_param_bfs(); /* pour ne pas le faire Ã  chaque sommet u */
  X->clean=1; /* initialisation complÃ¨te de X->D, puis partielle */

  int bmax,bsum; /* stat pour les "vraie" boules */
  bmax=bsum=0;
  
  for(u=0;u<n;u++){ /* pour chaque sommet u du graphe */

    /* Cherche d'abord le rayon r de B[u] pour calculer ensuite plus
       rapidement B[u] avec un BFS de profondeur r. Pour cela on
       remonte dans tree(center) jusqu'Ã  trouvÃ© un sommet de C.
       Cependant, ce n'est pas forcÃ©ment le plus proche de u. AprÃ¨s le
       BFS depuis u, il faut vÃ©rifier si dans la derniÃ¨re couche du
       BFS on a pas un sommet de C. Si oui, on dÃ©crÃ©mente le rayon et
       on modifie tree(center) de sorte Ã  trouver un sommet de C
       plutÃ´t et Ã©viter de futurs problÃ¨mes. Dans la variante "boules
       vides" on ne touche pas Ã  l'arbre. */

    if(VARIANT&1){ // variante "boules vides"
      if(X0->P[u]==-2) continue; // on ne fait rien si u est dans C
      B[u]=new_table(0); /* on crÃ©e une table B[u] vide, NB: B[u]->n=0 */
      B[u]->vpd=X0->P[u]; // voisin par dÃ©faut est le pÃ¨re de u
      continue;
    }

    t=v=u; r=-1;
    while(X0->P[v]!=-2){ /* si v pas dans C */
      v=X0->P[t=v]; /* t=v et v=parent(t) dans l'arbre X0->P */
      r++;
    }
    /* Ici:
       v=1er sommet ancÃªtre dans tree(center) de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1
    */
    if(r<0) continue; /* rien Ã  faire, on Ã©tait parti d'un sommet de C */
    
    B[u]=new_table(0); /* on crÃ©e la table B[u] vide, NB: B[u]->n=0 */
      
    if(r>0){ /* ne rien faire ici si r=0 (=si u voisin de C) */
      X->hmax=r; bfs(G,u,X); /* calcule un BFS depuis u de profondeur r */
      i=B[u]->n=(X->n)-1; /* taille de B[u] sans u, si r est correct */
      if(v==center){ /* il faut vÃ©rifier si u n'a pas de sommet v dans C Ã  distance r */
	/* on part des sommets i les plus Ã©loignÃ©s de u, donc Ã  distance r */
	/* NB: il y a toujours dans X->D au moins un sommet Ã  distance < r, c'est u */
	while((X->D[X->file[i]]==r)&&(X0->P[X->file[i]]!=-2)) i--;
	if(X0->P[X->file[i]]==-2){ /* on a trouvÃ© un sommet dans C Ã  distance r */
	  v=X->file[i]; /* v=sommet dans C le plus proche de u */
	  t=X->P[v]; /* t=dernier sommet de B[u] avant C (=pÃ¨re de v dans le BFS de u) */
	  X0->P[t]=v; /* modifie le pÃ¨re de t dans tree(center) pour accÃ©lÃ©rÃ©r les tests suivants */
	  while(X->D[X->file[i]]==r) i--; /* on peut effacer les sommets de i+1 Ã  X->n-1 */
	  B[u]->n=i; /* la taille de B[u] est i+1 avec u, donc i sans u */
	  r--; /* r=dist(u,v)-1=dist(u,t) */
	}
      }
    }
    
    /* Ici:
       v=le plus proche sommet de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1>=0
       B[u]->n=(taille de boule de rayon r sans u)-1
     */

    /* pour les stats avant ajoÃ»t des voisins et optimisation */
    bmax=max(bmax,1+B[u]->n);
    bsum += (1+B[u]->n);

    if(r==0){ /* si u voisin de C, alors B[u] contient certains voisins de u mais pas u */
      B[u]->vpd=v; /* v=voisin de u dans C */
      if(VARIANT&4){ /* alors B[u] = boule vide (sans u) */
	B[u]->n=B[u]->radius=0; // NB: B[u]->node et B[u]->dist sont Ã  NULL par new_table(0)
	continue;
      }
      if(VARIANT&8){ /* alors B[u] = voisins de u pas dans C */
        ALLOC(B[u]->node,G->d[u]); /* taille deg(u) pour avoir suffisamment de place */
        B[u]->n=0; // sert d'indice pour B[u]->node
        for(i=0;i<G->d[u];i++){
          v=G->L[u][i]; // v=voisin de u
          if(X0->P[v]!=-2) B[u]->node[B[u]->n++]=v; // si v pas dans C, on l'ajoute Ã  B[u]
        }
        B[u]->radius=(B[u]->n>0); // =1 ssi au moins un sommet dans B[u]
        REALLOC(B[u]->node,B[u]->n); // redimensionne B[u]
      }
      else{
	B[u]->n=G->d[u]; /* |B[u]|=deg(u) */
        B[u]->radius=1; /* rayon de B[u], car deg(u)>0 */ 
        ALLOCZ(B[u]->node,B[u]->n,G->L[u][_i]); /* copie tous les voisins de u dans B[u]->node */
      }
      ALLOCZ(B[u]->dist,B[u]->n,1); /* les sommets de B[u] sont tous Ã  distance 1 */
      continue; /* on passe au sommet u suivant */
    }

    /* ici r>0, et donc B[u]->n>0, et t<>u */
    /* on cherche le voisin par dÃ©faut de u, cad le fils de u dans
       BFS(u) qui mÃ¨ne Ã  t. LÃ , on pourrait tester le pÃ¨re dans X0
       pour Ã©ventuellement le corriger et Ã©viter les problÃ¨mes de
       correction de r en r-1 comme plus haut. */
    while(X->P[t]!=u) t=X->P[t];
    B[u]->vpd=t; /* voisin de u par dÃ©faut, ne peut pas Ãªtre dans C */
    B[u]->radius=r; /* rayon de la boule */
    
    /* copie les sommets (et leur distance) du BFS couvrant B[u] en supprimant u */
    ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* B[u]->node[i]=i-Ã¨me sommet de B[u] */
    ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* B[u]->dist[i]=distance du i-Ã¨me Ã  u */
  }
  printf("- time to construct tables B: %s\n",TopChrono(1));
  printf("- maximum non-cluster ball size: %i\n",bmax);
  printf("- average non-cluster ball size: %.2lf (%i/%i)\n",
	 (double)bsum/(double)n,bsum,n);
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");

  /*************************/
  /* optimise les tables B */
  /*************************/
  
  /*
    Dans un premier temps on met des -1 aux sommets Ã  supprimer.
    Ensuite on rÃ©organise les tables B en gardant l'ordre des sommets
    restant. NB: on peut avoir B[u]<>NULL et B[u]->n=0. B[u]->radius
    n'est pas mis Ã  jour et reprÃ©sente toujours la distance entre u et
    le sommet de B[u] le plus loin avant "clean". Il est important
    qu'avant l'appel, les sommets de B[u] soient rangÃ©s dans l'ordre
    du bfs(), donc non triÃ©es.
  */

  X->clean=1; /* initialisation complÃ¨te de X->D, puis partielle */
  for(u=0;u<n;u++){
    if(B[u]){ /* il faut que la table B[u] existe */
      /* si deg[u]<=1 */
      if(G->d[u]<2) t=0;
      else{
	/* si deg(u)>1 */
	v=B[u]->vpd; /* v=voisin par dÃ©faut */
	X->hmax=(B[u]->radius)-1;
	bfs(G,v,X); /* calcule un BFS depuis v de profondeur rayon de B[u] -1 */
	for(i=0;i<X->n;i++){ /* passe en revue les sommets du bfs(v...) */
	  t=X->file[i]; /* t=sommet du bfs(v,...) */
	  r=SetSearch(t,B[u]->node,B[u]->n,0); /* r=indice tq B[u]->node[r]=t */
	  if((r>=0)&&(B[u]->dist[r]==1+X->D[t])) B[u]->node[r]=-1; /* supprime t */
	}
	/* supprime les -1 de B[u] en dÃ©calant B[u]->node et B[u]->dist */
	for(i=t=0;i<B[u]->n;i++)
	  if(B[u]->node[i]>=0){
	    B[u]->node[t]=B[u]->node[i];
	    B[u]->dist[t]=B[u]->dist[i];
	    t++;
	  }
      }
      /* redimensionne les tables Ã  t = nouvelle taille de B[u] */
      B[u]->n=t;
      REALLOC(B[u]->node,t);
      REALLOC(B[u]->dist,t);
    }
  }
  printf("- time to clean tables B: %s\n",TopChrono(1));
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");

  DEBUG(for(u=0;u<n;u++)
	  if(B[u]){
	    printf("u=%i B[u]->vpd=%i B[u]->radius=%i ",u,B[u]->vpd,B[u]->radius);
	    if(B[u]) PRINTT(B[u]->node,B[u]->n);
	  });
  
  /************************/
  /* calcule les tables S */
  /************************/

  /*
    Tables seulement dÃ©finies pour les sommets u de C sauf center (car
    le centre est la racine de l'arbre et donc dÃ©jÃ  sur le plus court
    chemin).  S[u] est la liste des sommets de hash i oÃ¹ u=C[i] et qui
    sont Ã  distance au plus r=min(2,logn/loglogn) de u. On ne met ni u
    ni center dans S[u]. NB: Il est possible d'avoir dans S[u] des
    sommets de C (et mÃªme tous sauf u et center).
  */
  
  r=max(2,lg(n));
  r=ceil((double)r/(double)lg(r));
  int const depthS=max(2,r); /* peut-Ãªtre > hauteur(tree(center)) */

  ALLOCZ(S,n,NULL); /* tableau de n tables S (vides au dÃ©part) */
  X->hmax=depthS; /* profondeur max du BFS */
  X->clean=1; /* initialisation complÃ¨te de X->D, puis partielle */

  for(i=0;i<k;i++){ /* pour chaque sommet u du cluster sauf center */
    u=C[i]; if(u==center) continue;
    bfs(G,u,X); /* calcule un BFS depuis u<>center de profondeur depthS */
    t=min(X->n-1,F[i]-(H[u]==i)-(H[center]==i)); /* t=nombre max de sommets dans S[u] */
    S[u]=new_table(t); /* table S pour u */
    S[u]->n=0; /* on se sert de S[u]->n comme indice */
    /* on va parcourir les sommets du bfs() de profondeur <= depthS et
       ne garder que ceux de hash = i (ici i = hash(u) = hash(C[i]) */
    S[u]->vpd=-1;
    /* S[u]->vpd est l'indice (dans S[u]->node) du 1er sommet v de
       S[u] qui est Ã  une distance >= depthS-1 de u.  Dans S[u], les
       sommets sont rangÃ©s selon le bfs(), donc en fonction de leur
       distance Ã  u. On s'en sert plus tard pour accÃ©lÃ©rer la
       suppression dans R[u] des sommets v dÃ©jÃ  dans S[u] (lire
       commentaires sur le remplissage de R[u]). NB: si aucun sommet
       n'est de hash = i (mauvaise fonction de hashage !), alors
       S[u]->vpd<0. */
    for(j=1;(j<X->n)&&(S[u]->n<t);j++){ /* pour chaque v du bfs(), v<>u */
      /* lorsque S[u] contient t sommets, on peut s'arrÃªter. */
      v=X->file[j]; /* v=j-Ã¨me sommet de la file, v<>u */
      if((H[v]==i)&&(v!=center)){ /* v a le bon hash et v<>center */
	S[u]->node[S[u]->n]=v;
	S[u]->dist[S[u]->n]=X->D[v];
	if((X->D[v]>=depthS-1)&&(S[u]->vpd<0))
	  S[u]->vpd=S[u]->n; /* 1er sommet de hash i Ã  distance >= depthS-1 */
	S[u]->n++; /* un sommet de plus dans S[u] */
      }
    }
    if(S[u]->n){ /* redimensionne les tables de S[u] pour pas gaspiller */
      REALLOC(S[u]->node,S[u]->n);
      REALLOC(S[u]->dist,S[u]->n);
    }else{ /* si aucun sommet dans S[u], on efface S[u] */
      free_table(S[u]);
      S[u]=NULL;
    }
  }
  free_param_bfs(X); /* plus de BFS Ã  faire, X ne sert plus Ã  rien */
  printf("- time to construct tables S: %s (radius %i)\n",TopChrono(1),r);
  
  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |S[u]|=%i depth=%i\n",u,i,S[u]?S[u]->n:0,depthS);
	  if(S[u]){
	    PRINTT(S[u]->node,S[u]->n);   
	    PRINTT(S[u]->dist,S[u]->n);
	  }
	}
	);
  
  /************************/
  /* calcule les tables R */
  /************************/
  
  /*
    Tables seulement dÃ©finies pour les sommets u de C avec u=C[i].
    R[u] est la liste des sommets dont le hash est i et qui ne sont
    pas dÃ©jÃ  dans S.  On ne met ni u ni center dans R[u]. C'est donc
    comme S[u] mais sans la restriction de distance. Pour les sommets
    de R[u], on appliquera le routage dans tree(center).
  
    Les sommets de hash=couleur(center) sont tous dans R[center] (sauf
    center lui-mÃªme si hash(center)=couleur(center)). Si u<>center,
    alors pour qu'un sommet v (de hash i) soit dans R[u] il faut qu'il
    soit Ã  une distance >= depthS de center. S'il est Ã  distance =
    depthS-1 (ou infÃ©rieure), alors en passant par center on obtient
    une route de u Ã  v de longueur <= (depthS-1)+1=depthS ce qui n'est
    pas mieux que d'utiliser S[u]. Le sommet v sera donc dÃ©jÃ  dans
    S[u]. Maintenant, les sommets candidats v Ã  distance >= depthS de
    center ne peuvent Ãªtre qu'Ã  distance au moins depthS-1 de u. S'ils
    Ã©taient Ã  distance < depthS-1, alors ils seraient Ã  distance <
    depthS de center. Or, ils sont dÃ©jÃ  Ã  distance >= depthS. Donc
    pour chercher si v est dÃ©jÃ  dans S[u], on peut se contenter de
    vÃ©rifier les sommets de S[u] Ã  distance >= depthS-1 de u. Ces
    sommets candidats sont stockÃ©s Ã  partir de l'indice S[u]->vpd dans
    S[u]->node (car S n'est pas encore retriÃ©e). Dans R[u]->dist on
    stocke la distance rÃ©alisÃ©e dans l'arbre pour faire u->v.
  */
  
  ALLOCZ(R,n,NULL); /* tableau de tables R (vides au dÃ©part) */
  for(i=0;i<k;i++){ /* seulement pour les k sommets u de C */
    u=C[i]; /* u=sommet de C */
    t=F[i]-(H[u]==i); /* t=#sommets Ã  mettre dans R[u] (moins Ã©ventuellement u) */
    if((u!=center)&&(H[center]==i)) t--; /* enlÃ¨ve u et aussi center si H[center]=i */
    if(S[u]) t-=S[u]->n; /* on enlÃ¨ve ceux dÃ©jÃ  dans S[u] */
    if(t>0){ /* si table pas vide */
      R[u]=new_table(t); /* table R pour u */
      R[u]->n=0; /* on se sert de R[u]->n comme indice */
    }
  }
  free(F); /* ne sert plus Ã  rien */
  /* on remplit les tables R[u]->node et R[u]->dist */
  /* R[u]->dist[i]=dist. dans tree(center) entre u et le i-Ã¨me sommet de R[u] */
  for(v=0;v<n;v++){ /* Ã©crit chaque sommet v dans la bonne table R[u] */
    u=C[H[v]]; /* u=sommet de C en charge du hash de v (donc de couleur H[v]) */
    if((u==v)||(v==center)) continue; /* on ne met ni u ni center dans R[u] */
    if((X0->D[v]<depthS)&&(u!=center)) continue; /* seulement pour v Ã  dist > depthS-1 de center */
    /* On cherche si v n'est pas dÃ©jÃ  dans S[u]. Il suffit de chercher
       des sommets v Ã  distance >= depthS-1 de u */ 
    if(S[u]){ /* cherche v dans S[u], si elle n'est pas vide bien sÃ»r */
      j=S[u]->vpd; /* 1er sommet Ã  distance >= depthS-1 de u */
      if(j<0) j=S[u]->n; /* ne rien faire si pas de 1er sommet (est-ce nÃ©cessaire ?) */
      for(;j<S[u]->n;j++) if(S[u]->node[j]==v) break; /* v est dans S[u] */
    }
    if((S[u]==NULL)||(j==S[u]->n)){ /* v n'a pas Ã©tÃ© trouvÃ© dans S[u] */
      R[u]->node[R[u]->n]=v; /* ajoute v Ã  R[u] */
      // Calcule R[u]->dist: on remonte de v jusqu'au cluster et on
      // voit si l'on passe par u ou pas. La distance d de v Ã  u dans
      // tree(center) est la suivante:
      //   si u=center, alors          d=X0->D[v]
      //   si u<>center et t=u, alors  d=X0->D[v]-1
      //   si u<>center et t<>u, alors d=X0->D[v]+1
      if(u==center) /* remontÃ©e inutile si u est le center */
	R[u]->dist[R[u]->n]=X0->D[v];
      else{ /* on remonte Ã  partir de v jusqu'Ã  t=sommet dans C */
	t=v; while(X0->P[t]!=-2) t=X0->P[t]; /* si t pas dans C, t=parent(t) */
	R[u]->dist[R[u]->n]=X0->D[v]+((t==u)?-1:1);
      }
      R[u]->n++; /* un sommet de plus dans R[u] */
    }
  }
  printf("- time to construct tables R: %s\n",TopChrono(1));
  free_param_bfs(X0); /* ne sert plus Ã  rien */

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |R[u]|=%i\n",u,i,R[u]?R[u]->n:0);
	  if(R[u]) PRINTT(R[u]->node,R[u]->n);   
	  if(R[u]) PRINTT(R[u]->dist,R[u]->n);
	}
	);
    
  /************************/
  /* calcule les tables W */
  /************************/

  /*
    Tables seulement dÃ©finies pour tous les sommets u de C, W[u] donne
    la liste des couleurs (=indice dans C) des voisins de u dans C. La
    couleur du center n'est jamais ajoutÃ© Ã  W[u], car c'est le port
    par dÃ©faut. NB: Le nom du sommet de couleur i est C[i]. On ne se
    sert pas de W[u]->dist. Si W[u] est vide, alors on ne la supprime
    pas de faÃ§on Ã  avoir le port par dÃ©faut. Dans le cas bit-1 de
    VARIANT=1, toutes les tables sont vides (routage dans l'Ã©toile
    sans table pour le center et les feuilles).
  */

  ALLOCZ(W,n,NULL); /* tableau des tables W (vides au dÃ©part) */
  for(i=0;i<k;i++){ /* seulement pour les sommets u de C, y compris center */
    u=C[i]; /* u=sommet de C de couleur i */
    W[u]=new_table(0); /* table W pour u, NB: W[u]->n=0 */
    ALLOC(W[u]->node,k-1); /* au plus k-1 voisins */
    W[u]->vpd=center; /* voisin par dÃ©faut = center */
    if(VARIANT&2) continue; /* toutes les tables W seront vides */
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-Ã¨me voisin de u */
      if(B[v]||(v==center)) continue; /* on veut v dans C et v<>center */
      W[u]->node[W[u]->n++]=SetSearch(v,C,k,1); /* ajoute la couleur de v Ã  W[u]->node */
    }
    if(W[u]->n) /* rÃ©ajuste la taille */
      REALLOC(W[u]->node,W[u]->n);
  }
  printf("- time to construct tables W: %s\n",TopChrono(1));

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |W[u]|=%i\n",u,i,W[u]?W[u]->n:0);
	  if(W[u]) PRINTT(W[u]->node,W[u]->n);   
	}
	);

  /*******************/
  /* trie les tables */
  /*******************/
  for(u=0;u<n;u++){
    if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
    if(R[u]) QSORT2(R[u]->node,R[u]->dist,R[u]->n,fcmp_int);
    if(S[u]) QSORT2(S[u]->node,S[u]->dist,S[u]->n,fcmp_int);
    if(W[u]) QSORT(W[u]->node,W[u]->n,fcmp_int);
  }
  printf("- time to sort tables B,R,S,W: %s\n",TopChrono(1));
  BARRE;
  
  /* calcule de la taille Z[u] des tables de chaque sommet u */
  ALLOCZ(Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(S[u]) Z[u] += S[u]->n;
    if(R[u]) Z[u] += R[u]->n;
    if(W[u]) Z[u] += W[u]->n;
  }

  /* affiche taille min/max et moyenne/Ã©cart type des diffÃ©rentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(S[C[_i]]->n,k,S[C[_i]],"table S size");
  MINMAXMOY(R[C[_i]]->n,k,R[C[_i]],"table R size");
  MINMAXMOY(W[C[_i]]->n,k,W[C[_i]],"table W size");

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus Ã  rien */
  BARRE;

  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_cluster_tables,RT,1);
  RT->B=B;
  RT->S=S;
  RT->R=R;
  RT->W=W;
  RT->C=C; // besoin pour le routage avec W
  RT->H=H; // besoin du hash des sommets
  RT->n=n; // besoin pour libÃ©rer les n tables
  RT->center=center; // besoin pour le routage en Ã©toile (tables W vides)

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_cluster_length(int u,int v,rs_cluster_tables *X)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_cluster() pour router un message de u Ã  v, ou bien -1 si la route
  n'a pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables nÃ©cessaire au
  schÃ©ma, notamment les tables B,S,R,W. Si u<0 alors on teste la
  validitÃ© des tables (et on renvoie une valeur non-nulle en cas
  d'erreur). L'algorithme dÃ©pend du bit-1 de VARIANT.

  AmÃ©lioration possible: Lorsque u a des voisins dans C, alors, plutÃ´t
  que d'aller vers ->vdp, choisir un landmark de C parmi un ensemble
  prÃ©dÃ©terminÃ© de p = ceil(2m/n) = O(1) sommets (le degrÃ© moyen)
  d'Ã©ventuellement la bonne couleur h(v). Aller vers ->vdp seulement
  si cette recherche Ã  Ã©chouÃ©. Pour le schÃ©ma thÃ©orique, on peut
  dÃ©couper la table B en deux, B1 et B2, oÃ¹ B2 serait les sommets
  voisins de u dans C, et doubler B2 (double accÃ¨s par sommet et par
  couleur) de faÃ§on Ã  garantir un temps constant dans tous les cas, et
  pas seulement un temps Ã©gale au degrÃ© moyen (garanti seulement avec
  grande proba dans les RPLG.
*/
{
  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->R==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivÃ©
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // si u n'est pas dans C
  
  if(X->B[u]){
    // v dans B[u] ?
    i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]:u -> v
    if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
    // si v pas dans B[u]: u -> B[u]->vpd -> v
    DEBUG(printf("via B[u]->vpd=%i\n",X->B[u]->vpd););
    return 1 + rs_cluster_length(X->B[u]->vpd,v,X);
  }
  
  // si u est dans C

  if(X->S[u]){ // il faut S non vide, donc u<>center
    // v dans S[u] ?
    i=SetSearch(v,X->S[u]->node,X->S[u]->n,1);
    // si v est dans S
    if(i>=0){ DEBUG(printf("in S[u]\n");); return X->S[u]->dist[i];}
  }

  // si v n'est pas dans S
  // v dans R[u] ?

  if(X->R[u]){ // il faut R non vide
    i=SetSearch(v,X->R[u]->node,X->R[u]->n,1);
    // si v est dans R, on cherche si v descendant de u dans tree(center)
    if(i>=0){ DEBUG(printf("in R[u]\n");); return X->R[u]->dist[i]; }
  }
  
  // si v n'est pas ni dans R ni dans S
  // H[v] dans W[u] ?

  if(X->W[u]){ /* ne peut pas Ãªtre vide */
    i=SetSearch(X->H[v],X->W[u]->node,X->W[u]->n,1);
    // si H[v] est dans W[u]: u -> C[H[v]] -> v
    // si H[v] pas dans W[u]: u -> center -> v
    if((i>=0)||((VARIANT&2)&&(u==X->center))) u=X->C[X->H[v]];
    else u=X->W[u]->vpd;
    DEBUG(
	  if(u==X->C[X->H[v]])
	    printf("%s: u -> C[H[v]]=%i\n",(VARIANT&2)?"sans tables W":"H[v] in W[u]",u);
	  else
	    printf("%s: u -> W[u]->vpd=%i\n",(VARIANT&2)?"sans tables W":"H[v] not in W[u]",u);
	  );
    return 1 + rs_cluster_length(u,v,X);
  }

  // y'a un problÃ¨me
  DEBUG(printf("fail: W[u] does not exist\n"););
  return FAIL_ROUTING;
}


rs_dcr_tables *rs_dcr(graph* const G,int k)
/*
  Calcule pour le graphe G le routing scheme DCR ou AGMNT de paramÃ¨tre
  k>0 correspondant au nombre de couleurs. On renvoie les tables ainsi
  calculÃ©es. Le stretch est toujours <= 5 (DCR) ou <=3 pour AGMNT. La
  distinction entre DCR et AGMNT se fait par le bit-2 de VARIANT. Les
  champs ->CONT et ->F qui sont Ã©galement NULL pour DCR. Le bit-0 de
  VARIANT est rÃ©servÃ© au choix des landmarks de plus haut degrÃ©,
  variante valable aussi bien pour DCR que pour AGMNT.

  Si une couleur n'existe pas, alors cela marche quand mÃªme: les
  boules couvrent tout le graphe et le stretch est alors 1.

  SpÃ©cificitÃ©es pour AGMNT: Pour le routage de u -> v, oÃ¹
  hash(v)=C(u), il faut prendre la meilleure des options parmi:

  1) routage via un landmark s et w=nca(u,v) dans T(s). Dans ce cas,
     on a besoin de stocker s, Ã©ventuellement w pour Ã©viter de le
     recalculer.

  2) routage via les boules contiguÃ«s via une arÃªte x-y. Plus
     prÃ©cisÃ©ment, on prend la plus courte des routes, si elles
     existent, de la forme u->s->x-y->v tel que B(s) contient u et x,
     B(v) contient y, et x,y sont voisins. NB: u=s et y=v sont
     possibles.

  La route u->s->x-y->v, si elle existe, a la propriÃ©tÃ© que s->x-y->v
  est un plus court chemin. De plus le plus petit ancÃªtre commun dans
  T(s) entre u et x doit Ãªtre s. Si c'Ã©tait un sommet w<>s, alors w
  serait un meilleur candidat vÃ©rifiant toutes les propriÃ©tÃ©s. Donc le
  sommet x Ã  chercher ne peut pas Ãªtre un descendant de u' dans T(s)
  oÃ¹ u' est le voisin de s menant Ã  u. Donc la longueur du routage est
  alors d(u,s)+d(s,x)+1+d(y,v).

  La recherche des routes selon les boules contiguÃ«s peut se faire
  ainsi (le calcul rÃ©ellement rÃ©alisÃ© est cependant un peu diffÃ©rent):

  pour tout sommet u:
    pour tout sommet s de Bâ»Â¹(u):
      pour tout x de B(s) qui n'est pas descendant de u':
        pour tout voisin y de x:
	  pour tout v de Bâ»Â¹(y) tq hash(v)=C(u):
	    calculer d(u,s)+d(s,x)+1+d(y,v).

  Il faut calculer Bâ»Â¹(u) de chaque sommet, puis trier leurs sommets
  selon leur hash de sorte qu'on puisse trouver trÃ¨s rapidement tous
  les sommets de Bâ»Â¹(y) ayant un hash donnÃ© (=C(u)). La complexitÃ©
  est grosso-modo O(m*n*ln(n)). Les distances entres les centres des
  boules n'ont pas Ã  Ãªtre calculÃ©es.

  Pour savoir si x est descendant ou pas de u', il faut deux choses:
  1) pour chaque sommet x de B(s), avoir un numÃ©ro dfs(x) dans l'arbre
  T(s) couvrant B(s).  2) pour chaque sommet s de Bâ»Â¹(u), avoir
  l'intervalle [a,b] des dfs(x) pour les sommets x de B(s) descendant
  de u'. Il faut donc stocker 1 entier (dfs) pour chaque sommet des
  boules, et deux entiers (a,b) pour chaque sommet d'une boule
  inverse. En pratique cela nÃ©cessite des calculs et de la mÃ©moire
  supplÃ©mentaires, et il n'est pas clair que cela soit plus efficace.
*/
{
  int const n=G->n;
  int const agmnt=VARIANT&4; // agmnt vrai ssi c'est le schÃ©ma AGMNT, sinon c'est DCR
  int u,v,w,i,j,t,c,nc,nl;

  printf("\n%s\n",agmnt? "AGMNT" : "DCR");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  printf("- wanted number of colors: %i\n",k);

  /**********************************/
  /* calcule H & C, hash et couleur */
  /**********************************/

  // H[u]=0..k-1, hash du sommet u
  // C[u]=0..k-1, couleur du sommet u
  // u est landmark ssi C[u]=0

  DEBUG(PRINT(HASH););
  int *H=MakeHash(NULL,n,k,HASH);
  NALLOC(int,C,n);

  if((VARIANT&1)&&(k>1)){
    // variante oÃ¹ les landmarks sont les sommets de plus haut
    // degrÃ©. On va les colorier 0 puis colorier les autres sommets
    // avec une couleur de 1 Ã  k-1. NB: il faut au moins 2 couleurs, k>1.
    for(u=0;u<n;u++) C[u]=1; // au dÃ©part les sommets n'ont pas de couleur
    t=iceil(n,k); // t=nombre de landmarks, estimÃ© Ã  la moyenne
    int *D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi); // trie le degrÃ© des sommets
    for(i=0;i<t;i++) C[D[n-i-1]]=0; // le i-Ã¨me sommet de plus haut degrÃ© est un landmark
    for(u=0;u<n;u++) if(C[u]) C[u]=1+randomu(k-1); // colorie dans [1,k[ les sommets non landmark
  }else
    // variante par dÃ©faut: les sommets ont une couleur alÃ©atoire dans [0,k[
    for(u=0;u<n;u++) C[u]=randomu(k); // il faut k>0

  DEBUG(PRINTT(C,n);PRINTT(H,n););

  /* affichage de stats sur H & C */
  /* calcule:
     nl=nombre de landmarks
     nc=nombre de couleurs
     F[c]=nombre de sommets de hash c
  */
  int *F;
  FREQMINMAX(F,k,C,n,"the colors"); /* frÃ©quence des couleurs (ne sert qu'Ã  l'affichage) */
  nl=F[0]; /* nl=nombre de landmarks, ceux de couleurs 0 */
  for(i=nc=0;i<k;i++) nc += (F[i]>0); /* nc=nombre de couleurs diffÃ©rentes, nc<=k */
  free(F); /* nÃ©cessaire car rÃ©allouÃ© par le prochain FREQMINMAX */
  FREQMINMAX(F,k,H,n,"the hash"); /* frÃ©quence des hashs (sert pour taille des tables) */
  for(i=u=0;i<k;i++) u += (F[i]>0); /* u=nombre de hashs diffÃ©rents, u<=k */

  printf("- real number of colors: %i\n",nc);
  printf("- real number of hashs: %i\n",u);
  printf("- number of landmarks: %i\n",nl);
  printf("- time to compute hash, color & stats: %s\n",TopChrono(1));
  BARRE;

  /*********************************************/
  /* calcule S, la liste des bfs des landmarks */
  /*********************************************/

  // S[u]=bfs(u,G) pour u landmark, NULL si u non landmark
  
  NALLOC(param_bfs*,S,n);
  for(u=0;u<n;u++)
    if(C[u]) S[u]=NULL; // C[u]=0 ssi u landmark
    else{
      S[u]=bfs(G,u,NULL);
      free(S[u]->file); /* libÃ¨re les pointeurs inutilisÃ©s */
      S[u]->file=NULL;  /* important pour le free_param_bfs() plus tard */
    }
  printf("- time to construct landmark bfs (array S): %s\n",TopChrono(1));

  DEBUG(
	for(u=0;u<n;u++)
	  if(C[u]==0){
	    PRINT(u);
	    PRINTT(S[u]->P,n);
	    PRINTT(S[u]->D,n);
	    printf("\n");
	  }
	);
  
  /*****************************/
  /* calcule les tables B et W */
  /*****************************/

  /*
    B[u]=boule de voisinage de u, la plus "petite" contenant toutes
         les couleurs (composÃ© de la liste des sommets et leur
         distance Ã  u), la derniÃ¨re couronne Ã©tant ordonnÃ©e selon les
         identifiants des sommets.

	 Il est important que les boules contiennent k couleurs, de
         faÃ§on Ã  toujours pouvoir Ãªtre capable de router vers un hash
         donnÃ©. Si une couleur n'apparaÃ®t pas dans le graphe, alors
         chaque boule devra contenir tous les sommets.

	 B[u] contient toujours au moins u, bien qu'on ne se sert
         jamais de cette entrÃ©e pour le routage. Par contre, pour la
         construction, c'est plus simple d'avoir u dans B[u]. On
         enlÃ¨vera donc cette entrÃ©e pour le routage (pour les stats et
         Ã©tablir la taille des tables), aprÃ¨s avoir construit toutes
         les tables.

    B[u]->node[i]=i-Ã¨me sommet de la boule de u
    B[u]->dist[i]=distance entre u et B[u]->node[i]
    B[u]->vpd=+proche landmark de u (+proche sommet de couleur 0 dans B[u])
    W[u]->node[c]=next-hop vers le 1er sommet de couleur c de B[u]->node
    (la table W, indexÃ©e par les couleurs, sert pour avoir le temps constant)

    Pour AGMNT:

    K[u]=taille de la boule inverse de B[u], ie K[u]=|{v : u in B[v]}|
    B[u]->radius=rayon de B[u]

    Algorithme:

      On fait un bfs(G,u) par couche (avec cont=1). On trie chaque
      couche par identifiant (si on a visitÃ© au moins assez de
      sommets), et on regarde Ã  chaque fois le premier moment oÃ¹ l'on
      voit une couleur donnÃ©e (remplissage de W). Au passage on
      dÃ©termine, lorsque la couleur est 0, le +proche landmark t de
      u. On s'arrÃªte quand on a visitÃ© toutes les couleurs ou bien que
      tous les sommets du graphe ont Ã©tÃ© visitÃ©s (cela peut arriver si
      une couleur n'est portÃ©e par aucun sommet).
   */

  NALLOCZ(table*,B,n,NULL);  // B=tableau de n tables B (vides au dÃ©part)
  NALLOCZ(table*,W,n,NULL);  // W=tableau de n tables W (vides au dÃ©part)
  param_bfs *X=new_param_bfs(); // X=le rÃ©sultat du bfs() depuis u
  NALLOC(int,T,n); /* T=tableau pour la derniÃ¨re couche de sommets */
  int *K; if(agmnt) ALLOCZ(K,n,0); /* K=taille des boules inverses */
  X->clean=1; /* initialisation complÃ¨te des distances, puis partielle */

  for(u=0;u<n;u++){ /* calcule B[u] & W[u] pour chaque u */
    c=0; /* c=nombre de couleurs (diffÃ©rentes) dÃ©jÃ  rencontrÃ©es */
    X->cont=1; /* Ã©vite de recalculer le dÃ©but de l'arbre */
    X->hmax=-1; /* au dÃ©part on visite seulement u */
    W[u]=new_table(0); /* W[u]->node[c]=-1 si la couleur c n'a pas Ã©tÃ© visitÃ©e */
    ALLOCZ(W[u]->node,k,-1); /* NB: on n'utilise pas ->dist */
    W[u]->n=k; /* taille de W[u] */
    do{
      
      X->hmax++; /* couche suivante (au dÃ©part hmax=0) */
      bfs(G,u,X); /* on parcoure tous les sommets jusqu'Ã  distance (incluse) hmax de u */
      for(j=0,i=X->tf;i<X->n;i++) T[j++]=X->file[i]; /* copie la derniÃ¨re couche dans T */
      if((c+j>=k)||(X->n>=n)) QSORT(T,j,fcmp_int); /* trie la derniÃ¨re couche (=T), obligÃ© Ã  cause de W[u] */
      /* on ne trie T que si on a le potentiel pour avoir terminÃ©:
	 avoir toutes les couleurs ou avoir visitÃ© tout le graphe.  Si
	 le nombre de sommets de la derniÃ¨re couche (=j) + le nombre
	 de couleurs dÃ©jÃ  rencontrÃ©es (=c) est au moins k alors on a
	 potentiellement atteint la derniÃ¨re couche. */

      DEBUG(
	    //PRINT(u);PRINT(X->n);PRINT(X->tf);
	    //PRINTT(T,j);PRINT(c);PRINT(j);
	    //printf("\n");
	    );
      
      /* remplit W en parcourant les sommets de la derniÃ¨re couche */
      /* le +proche landmark de u, lorsque rencontrÃ©, est stockÃ© dans t */
      for(i=0;(i<j)&&(c<k);i++){
	v=T[i]; /* v=i-Ã¨me sommet de la derniÃ¨re couche */
	if(W[u]->node[C[v]]>=0) continue; /* couleur dÃ©jÃ  rencontrÃ©e */
	/* ici la couleur C[v] n'a jamais Ã©tÃ© rencontrÃ©e */
	w=v; /* NB: si v=u, alors il faut effacer le -1 dans W[u] */
	if(v!=u) while(X->P[w]!=u) w=X->P[w]; /* cherche le next-hop w depuis u vers v */
	W[u]->node[C[v]]=w; /* met le next-hop dans W[u], ou bien v si u=v */
	if(C[v]==0) t=v; /* si v est un landmark, alors c'est le +proche de u */
	c++; /* une couleur de plus */
      }

      /* On sort de la boucle si: soit on a toutes les couleurs (=k)
	 dans W[u]->node, ou bien on a visitÃ© tout le graphe, ce qui
	 est possible si toutes les couleurs n'Ã©taient pas
	 reprÃ©sentÃ©es.  NB: ici, dans tous les cas, i est nombre de
	 sommets de la derniÃ¨re couche Ã  recopier partiellement. */
      
    }while((c<k)&&(X->n<n));

    /* on construit B[u] Ã  partir de X->file et de T */
    B[u]=new_table(X->tf+i); // table de B[u]
    B[u]->vpd=t; // t=+proche landmark de u=W[u]->node[0]
    B[u]->radius=X->hmax; // rayon de B[u], ne sert Ã  rien pour DCR

    /* construit B[u] */
    /* j=indice pour B[u]->node[] */
    /* t=indice pour X->file[] */
    for(j=t=0;t<X->tf;j++){ /* copie X->file sauf la derniÃ¨re couche */
      v=X->file[t++];
      B[u]->node[j]=v; if(agmnt) K[v]++; // car v dans B[u]
      B[u]->dist[j]=X->D[v];
    }
    for(t=0;t<i;j++){ /* copie la partie traitÃ©e de la derniÃ¨re couche */
      v=T[t++];
      B[u]->node[j]=v; if(agmnt) K[v]++; // car v dans B[u]
      B[u]->dist[j]=X->D[v];
    }
  }
  printf("- time to construct tables B & W: %s\n",TopChrono(1));
  free_param_bfs(X);
  free(T); /* ne sert plus Ã  rien */

  /* tri des tables B */
  /* c'est important pour le routage et pour optimiser AGMNT */
  /* NB: B[u]->n>0, car u est toujours dans B[u] */
  for(u=0;u<n;u++) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;


  /*****************************/
  /* Partie spÃ©cifique Ã  AGMNT */
  /*****************************/

  contigue **CONT=NULL;

  if(agmnt){
    /*********************************/
    /* calcule les boules inverses I */
    /*********************************/

    // I[u][i]=i-Ã¨me sommets de Bâ»Â¹[u] triÃ© selon la distance Ã  u
    // K[u]=nombre de sommets dans I[u] (dÃ©jÃ  calculÃ©)
    //
    // Les boules inverses ainsi triÃ©es servent Ã  accÃ©lÃ©rer le calcul
    // des routes via les boules contiguÃ«s.

    // alloue les boules inverses (I) et les distances au centre (D)
    NALLOC(int*,I,n);
    NALLOC(int*,D,n); // D ne sert que pour le tri des boules inverses
    for(u=0;u<n;u++){
      ALLOC(I[u],K[u]); // allocation des boules inverses
      ALLOC(D[u],K[u]); // allocation des distances au centre
    }

    // remplit les boules inverses et les distances au centre
    for(u=0;u<n;u++) K[u]=0; // K[u]=nombre de sommets dÃ©jÃ  mis dans I[v]
    for(u=0;u<n;u++) // pour toutes les boules
      for(i=0;i<B[u]->n;i++){ // pour tous les sommets v de B[u]
	v=B[u]->node[i];
	I[v][K[v]]=u; // ajoute u Ã  I[v]
	D[v][K[v]]=B[u]->dist[i]; // dist(u,v)
	K[v]++;
      }
    // ici K[u]=taille de I[u] pour tout u
    printf("- time to construct inverse balls (array I): %s\n",TopChrono(1));
    DEBUG(
	  PRINTT(K,n);
	  for(u=0;u<n;u++){
	    PRINT(u);
	    PRINTT(I[u],K[u]);
	    PRINTT(D[u],K[u]);
	    printf("\n");
	  }
	  );
    /* trie les boules inverses selon la distance Ã  u */
    /* NB: B[u] contient u, alors I[u] est de taille K[u]>0 */
    for(u=0;u<n;u++) QSORT2(D[u],I[u],K[u],fcmp_int); // trie I[u] selon D[u]
    FREE2(D,n); // ne sert plus Ã  rien
    printf("- time to sort inverse balls: %s\n",TopChrono(1));

    /****************************************************/
    /* calcule le meilleur chemin de u Ã  v de hash C[u] */
    /****************************************************/

    // CONT[u][i] = meilleur chemin pour aller de u Ã  v, le i-Ã¨me
    // sommet dont le hash vaut C[u] ordonnÃ© selon i
    //
    // UC[i]=i-Ã¨me sommet triÃ© u selon sa couleur C[u]
    // UH[i]=i-Ã¨me sommet triÃ© v selon son hash H[v] puis selon v
    //
    // Ces deux tableaux servent Ã  lister en temps n*k*H(k) ~ n^1.5
    // les paires (u,v) avec hash de v = C[u], les sommets de mÃªme
    // hash se retrouvant consÃ©cutifs dans UH. On pourrait construire
    // UC et UH en temps O(n) (au lieu de O(nlogn) comme on le fait)
    // en rÃ©utilisant les tableaux de frÃ©quence des hashs et des
    // couleurs. Mais c'est plus complexe, et le gain en temps sera au
    // final nÃ©gligeable.
    
    // construit la liste UC
    NALLOCZ(int,UC,n,_i); // UC[i]=i au dÃ©part
    fcmp_tabint(NULL,C); // tri selon C
    QSORT(UC,n,fcmp_tabint); // tri UC selon C

    // construit la liste UH
    NALLOCZ(int,UH,n,_i); // UH[i]=i au dÃ©part
    fcmp_tabinteq(NULL,H); // tri selon H
    QSORT(UH,n,fcmp_tabinteq); // tri UH selon H et v

    printf("- time to sort nodes according to hash and colors: %s\n",TopChrono(1));

    DEBUG(
	  PRINTT(UC,n);
	  PRINTT(UH,n);
	  );
    
    int l,d,d1,d2,s,h,p,q,r,x,y;
    ALLOC(CONT,n); // CONT[u] va Ãªtre allouÃ© pour tous les sommets u

    // Calcule la distance via les landmarks ou les boules contiguÃ«s.
    // On balaye les sommets u dans l'ordre croissant des couleurs.
    // Les sommets v, de hash C[u], sont consÃ©cutifs dans UH.

    // Pour les landmarks u, le routage vers v (de hash 0 donc)
    // s'effectue comme un routage via le landmark u (c'est ici un
    // meileur choix que le landmark de v), ce qui va produire une
    // route de plus court chemin dans T_u.

    for(j=t=0;j<n;j++){ // pour tous les sommets
      u=UC[j];   // u=sommet courant
      h=F[C[u]]; // h=nombre de sommets de hash C[u], h=0 est possible
      ALLOC(CONT[u],h); // alloue pour h sommets (ne fait rien si h=0)
      
      // calcule l'indice t de UH[t..t+h[ oÃ¹ sont rangÃ©s les v de hash C[u]
      q=j? UC[j-1] : u; // q=sommet juste avant u dans UC, q=u au dÃ©part (si j=0)
      if(C[u]>C[q]) t+=F[C[q]]; // changement de couleur ?

      // balaye tous les sommets v de hash C[u], de plus par ordre
      // croissant (important pour faire une recherche binaire au
      // moment du routage vers v)

      for(i=0;i<h;i++){
	v=UH[t+i]; // v=sommet de hash C[u]
	CONT[u][i].v=v; // stocke v dans tous les cas, important pour
			// la recherche binaire lors du routage
	if((C[v]==0)||(u==v)) continue; // ne rien Ã  faire si v est landmark ou si u=v
	if(C[u]==0){ // NB: on peut avoir t=0 et C[u]<>0 si F[0]=0 par exemple
	  // si u est landmark, on provoque un routage de plus court
	  // en codant une solution via boule-contiguÃ« de distance
	  // optimale dist(u,v) dans T_u. NB: CONT[u][i].w n'est pas dÃ©fini
	  CONT[u][i].s=-1; // pour dire via boule-contiguÃ« et forcer un +cc
	  CONT[u][i].d=S[u]->D[v]; // dist(u,v)
	  continue; // v suivant
	}
	if(SetSearch(v,B[u]->node,B[u]->n,1)>=0) continue; // si v est dans B[u], on a rien Ã  faire
	// ici v n'est pas un landmark et pas dans B[u] => dist(u,v)>=B[u]->radius

	q=SetSearch(u,B[v]->node,B[v]->n,1); // est-ce que v est dans Bâ»Â¹[u] ?
	if(q>=0){ // alors routage +cc
	  CONT[u][i].s=-1; // pour dire via boule-contiguÃ« et forcer un +cc
	  CONT[u][i].d=B[v]->dist[q]; // dist(u,v)
	  continue; // v suivant
	}

	// route via landmark ? On vise une route du type u->w->v,
	// avec w=nca(u,v,T_s) pour un certain landmark s. On commence
	// par le +proche landmark s de v (ce qui suffit pour garantir
	// un stretch <= 3), si bien qu'on ne change de landmark que si
	// on fait strictement mieux.

	CONT[u][i].d=INT_MAX; // important car le champs d pourrait ne pas Ãªtre initialisÃ© (si nl=0)
	for(l=-1;l<nl;l++){ // balaye tous les landmarks et aussi B[v]->vpd
	  s=(l<0)? B[v]->vpd : UC[l]; // s=landmark=un sommet de couleur 0 ou bien B[v]->vpd
	  w=nca_bfs(u,v,S[s]); // w=nca(u,v,T_z)
	  d=dist_nca(u,v,w,S[s]->D); // dist(u,v) dans T_s
	  if(d<CONT[u][i].d){ // on a trouvÃ© une meilleure route via w dans T_s
	    CONT[u][i].s=s;
	    CONT[u][i].w=w;
	    CONT[u][i].d=d;
	  }
	}
	DEBUG(
	      if(u==28 && v==36){
		printf("via best landmark:\n");
		PRINT(CONT[u][i].s);
		PRINT(CONT[u][i].w);
		PRINT(CONT[u][i].d);
	      }
	      );
	
	// route via boule-contiguÃ« ? On vise une route du type
	// u->s->x-y->v.  La plus courte de ces routes fait que
	// nÃ©cessairement s=nca(u,x,T_s) et que s->x-y->v est un +cc.
	// NB: grÃ¢ce au tri des Bâ»Â¹[u], les sommets s sont parcourus
	// par distance croissante Ã  u. Cela permet ainsi de stoper
	// plus rapidement la recherche d'un bon s.

	for(l=0;l<K[u];l++){ // balaye tous les sommets s de Bâ»Â¹[u] = I[u]
	  s=I[u][l]; // s=sommet contenant u dans B[s], NB: s<>v car v pas dans Bâ»Â¹[u] qui contient s
	  d=B[s]->dist[SetSearch(u,B[s]->node,B[s]->n,1)]; // d=|u->s|
	  q=SetSearch(v,B[s]->node,B[s]->n,1);
	  DEBUG(
		if(u==28 && v==36){
		  PRINT(s);
		  PRINT(CONT[u][i].d);
		  PRINT(B[s]->radius);
		  PRINT(d);
		  PRINT(q);
		  printf("\n");
		}
		);
	  if(q>=0){ // si v est dans B[s]
	    d += B[s]->dist[q]; // d=|u->s->v|
	    if(d<CONT[u][i].d){ // on a trouvÃ© une route meilleure pour u->v
	      CONT[u][i].d=d;   // d=distance trouvÃ©e
	      CONT[u][i].s=-1;  // pour dire via boule-contiguÃ«
	    }
	    continue; // on peut passer au s suivant
	  }
	  // ici v n'est pas dans B[s], donc dist(s,v) >= B[s]->radius >= 1 (s<>v)
	  if(d+1>=CONT[u][i].d) break; // aucun s ne pourra faire mieux -> v suivant
	  if(d+B[s]->radius>=CONT[u][i].d) continue; // ce s n'est pas assez bon -> s suivant
	  for(r=0;r<B[s]->n;r++){ // balaye tous les sommets x de B[s], x=s compris
	    // on suppose que x (=B[s]->node[r]) est le plus loin
	    // possible de s (ou le +proche de v). Il doit alors Ãªtre
	    // Ã  distance rayon de B[s] ou rayon de B[s] - 1 de
	    // s. C'est sans perte de gÃ©nÃ©ralitÃ© car v n'est pas dans
	    // B[s] et x sur un +cc entre s et v.
	    if(B[s]->dist[r]<B[s]->radius-1) continue; // x n'est pas le plus loin possible
	    d1 = d+B[s]->dist[r]+1; // d1=|u->s->x-y|
	    if(d1>=CONT[u][i].d) continue; // la route ne fera pas mieux
	    x=B[s]->node[r]; // x=le sommet de B[s]
	    for(p=0;p<G->d[x];p++){ // balaye les y voisins de x
	      y=G->L[x][p]; // y=voisin de x
	      q=SetSearch(y,B[v]->node,B[v]->n,1); // on cherche y dans B[v]
	      DEBUG(
		    if(u==28 && v==36){
		      PRINT(d1);
		      PRINT(x);
		      PRINT(y);
		      PRINT(q);
		    }
		    );
	      if(q<0) continue; // il faut y dans B[v]
	      // ici on a trouvÃ© une route du bon type
	      d2 = d1+B[v]->dist[q]; // d2=d1+|y->v|=|u->s->x-y->v|
	      if(d2<CONT[u][i].d){   // on a trouvÃ© une route meilleure pour u->v
		CONT[u][i].d=d2;     // d2=distance trouvÃ©e
		CONT[u][i].s=-1;     // pour dire via boule-contiguÃ«
	      }
	    }
	  }
	}
      }
    }
    
    free(UH);
    free(UC);
    FREE2(I,n);
    printf("- time to compute best routes via landmark or contigue-ball: %s\n",TopChrono(1));
    BARRE;
  }// fin du "if(agmnt){ ..."
  
  DEBUG(
	for(u=0;u<n;u++){
	  PRINT(u);
	  PRINT(B[u]->vpd);
	  PRINT(B[u]->radius);
	  PRINTT(B[u]->node,B[u]->n);
	  PRINTT(W[u]->node,W[u]->n);
	  printf("\n");
	}
	if(agmnt){
	  for(u=0;u<n;u++){
	    printf("\nCONT[%i].v = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].v);
	    printf("\nCONT[%i].s = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].s);
	    printf("\nCONT[%i].d = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].d);
	    printf("\nCONT[%i].w = ",u);for(i=0;i<F[C[u]];i++)
					  if(CONT[u][i].s<0) printf("- ");
					  else printf("%i ",CONT[u][i].w);
	    printf("\n");
	  }
	}
	);

  /****************************/
  /* taille totale des tables */
  /****************************/

  /* pour chaque sommet, il faut en moyenne:
     
     T1: une boule -> k*H(k)
     T2: une table de landmark -> n/k
     T3: une table des hashs -> k-1 (Ã  cause du temps constant, sinon -> 0)
     T4: une table des sommets de mÃªme hash que C[u] -> n/k
   
     A cause de l'algorithme de routage, on peut optimiser les tables
     ainsi. Si v est dans B[u] et est un landmark, on peut l'enlever
     de T2, ce qui en moyenne revient Ã  enlever H(k) entrÃ©es Ã 
     T2. Donc |T2| = n/k - H(k) en moyenne. De mÃªme, on peut enlever
     de T4 les entrÃ©es des sommets qui sont dans T1 ou T2, mais le
     gain est trÃ¨s lÃ©ger.
   
     La fonction donnant le nombre d'entrÃ©es est alors
     
              f(k,n) = 2n/k + (k-1)*(H(k)+1)

     voir la fonction func1(). Une faÃ§on de le voir l'optimisation de
     T2 est que les entrÃ©es de B[u] pourraient Ãªtre dÃ©coupÃ©es en deux
     tables: B1[u] et B2[u] oÃ¹ B1[u] seraient les entrÃ©es pour les
     sommets non-landmarks et B2[u] pour les landmarks. Pour B2[u] la
     structure seraient un peu diffÃ©rente, elle aurait Ã  la fois celle
     de B1[u] et celle des landmarks.
  */

  NALLOCZ(int,Z,n,1); /* taille de la table de u, au moins 1 pour chaque u */
  NALLOCZ(int,Z1,n,nl); /* taille de la table des landmarks qui ne sont pas dans B[u] */
  NALLOCZ(int,Z2,n,F[C[_i]]); /* taille de la table des hash = C[u] hors de B[u] */

  for(u=0;u<n;u++){ /* pour chaque sommet u */
    for(i=0;i<B[u]->n;i++){ // calcule Z1[u] et Z2[u]
      v=B[u]->node[i]; // v dans B[u]
      Z1[u] -= (C[v]==0); /* enlÃ¨ve les landmarks v qui sont dans B[u] */
      Z2[u] -= (H[v]==C[u]); /* enlÃ¨ve les hashs C[u] qui sont dans B[u] */
    }
    Z[u] += B[u]->n-1; /* taille de B[u] moins le sommet u */
    Z[u] += Z1[u];     /* table des landmarks hors B[u] */
    Z[u] += k-1;       /* taille de W (rien pour la couleur 0) */
    Z[u] += Z2[u];     /* table des sommets dont le hash est la couleur de u hors de B[u] */
  }

  /* affiche taille min/max et moyenne/Ã©cart type des diffÃ©rentes tables */
  MINMAXMOY(B[_i]->n-1,n,1,"table B size"); /* on retire u de B[u] pour le routage */
  if(agmnt) MINMAXMOY(K[_i],n,1,"inverse table B size");
  MINMAXMOY(Z1[_i],n,1,"landmark table size");
  MINMAXMOY(k-1,n,1,"table W size");
  MINMAXMOY(Z2[_i],n,1,"own color/hash table size");

  /* pointeurs qui ne servent plus Ã  rien */
  free(Z1);
  free(Z2);
  if(agmnt) free(K);
  else{ // ne sert plus pour DCR
    free(F);
    F=NULL;
  }

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus Ã  rien */
  printf("- theoretical average: 2âˆš(n*ln(n*ln(n))) = %i\n",(int)(2*sqrt(n*log(n*log(n)))));
  BARRE;
  
  /* assemble les tables en une seule pour le retour */  
  NALLOC(rs_dcr_tables,RT,1);
  RT->B=B;
  RT->W=W;
  RT->S=S;
  RT->H=H;
  RT->C=C;
  RT->n=n;
  RT->CONT=CONT; // pour AGMNT seulement, NULL sinon
  RT->F=F; // pour AGMNT seulement, NULL sinon
  ALLOCZ(RT->dist,n,S[_i]?S[_i]->D:NULL); /* distances partielles */

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_dcr_length_rec(int u,int v,int w,int a,rs_dcr_tables *X)
/*
  Fonction rÃ©cursive donnant la longueur de la route de u Ã  v avec les
  tables X initialisÃ©es par dcr, et l'en-tÃªte (w,a):
  - w=landmark intermÃ©diaire (s ou t, w=s par dÃ©faut)
  - a=nca(u,v) dans l'arbre du landmark w (a<0 par dÃ©faut)
*/
{  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivÃ©
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // v dans B[u] ? (NB: B[u] existe toujours)
  
  i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
  // si v dans B[u]:u -> v
  if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
  
  // v n'est pas dans B[u]
  // est-ce que v est dans la table des landmarks ?

  if(X->C[v]==0){ // 0=couleur des landmarks
    DEBUG(printf("%i is a landmark\n",v););
    return X->S[v]->D[u]; // dist(u,v), v=landmark
  }

  // v n'est pas landmark
  // est-on arrivÃ© au responsable de v ?
  
  if((a<0)&&(X->C[u]!=X->H[v])){ // si pas encore au responsable
    u=X->W[u]->node[X->H[v]]; // next-hop de u vers le responsable de v
    DEBUG(printf("go to %i the closest node of color %i, the hash of %i (a=%i and w=%i)",u,X->H[v],v,a,w););
    return 1 + rs_dcr_length_rec(u,v,w,a,X);
  }

  // on est passÃ© par (ou on est sur) le responsable de v ?
  
  if(a<0){ // on est arrivÃ© au responsable de v
    int w1=X->B[w]->vpd; // landmark le +proche de la source (=w si a<0)
    int w2=X->B[v]->vpd; // landmark le +proche de cible (=v)
    param_bfs *Y1=X->S[w1]; // l'arbre de w1
    param_bfs *Y2=X->S[w2]; // l'arbre de w2
    int a1=nca_bfs(u,v,Y1); // ancÃªtre commun entre u et v dans Y1
    int a2=nca_bfs(u,v,Y2); // ancÃªtre commun entre u et v dans Y2
    int d1=dist_nca(u,v,a1,Y1->D); // d1=dist(u,v) dans Y1
    int d2=dist_nca(u,v,a2,Y2->D); // d2=dist(u,v) dans Y2
    if(d1<d2) w=w1,a=a1; else w=w2,a=a2;
    
    /* la ligne suivante est une optimisation: si d1=d2, alors on a
       intÃ©rÃªt de choisir l'ancÃªtre a1 ou a2 le plus loin de u cela
       laisse plus de chances Ã  l'algo de court-circuiter le routage
       dans l'arbre. */
    if((d1==d2)&&(Y1->D[u]-Y1->D[a1]>=Y2->D[u]-Y2->D[a2])) w=w1,a=a1;
    
    DEBUG(
	  printf("node in charge of %i reached\n",v);
	  printf("go to a=%i, the nca of %i and %i in the tree of landmark w=%i\n",a,u,v,w);
	  );
  }

  // ici on a:
  //  w=landmark intermÃ©diaire (landmark de s ou de t)
  //  a=ancÃªtre commun entre u et v dans l'arbre de w
  if(u==a){ // si u est arrivÃ© au bon ancÃªtre
    DEBUG(printf("ancestor a=%i reached\n",a););
    return X->S[w]->D[v] - X->S[w]->D[a];
  }

  // ici il faudrait ajouter les racourcis, si un des ancÃªtres de v
  // dans l'arbre de w (stockÃ©s dans l'Ã©tiquette de v) est dans B[u]
  
  u=X->S[w]->P[u]; // on remonte dans l'arbre de w
  DEBUG(printf("go up to a=%i in the tree of landmark w=%i",a,w););
  return 1 + rs_dcr_length_rec(u,v,w,a,X);
}


int rs_dcr_length(int u,int v,rs_dcr_tables *X)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_dcr() pour router un message de u Ã  v, ou bien -1 si la route n'a
  pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables nÃ©cessaires au
  schÃ©ma. Si u<0, on rÃ©alise quelques tests de bases sur les tables
  X. La fonction fait essentiellement appel Ã  une fonction recursive
  oÃ¹ la gestion d'un en-tÃªte est nÃ©cessaire.
*/
{
  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }

  return rs_dcr_length_rec(u,v,u,-1,X);
}


int rs_agmnt_length_rec(int u,int v,int s,int w,rs_dcr_tables *X)
/*
  Fonction rÃ©cursive donnant la longueur de la route de u Ã  v avec les
  tables X initialisÃ©es par agmnt, et l'en-tÃªte (s,w):
  - s=landmark pour un routage via landmark (<0 par dÃ©faut)
  - w=nca(u,v,T_s) si s>=0

  L'algorithme est le suivant:
  1. si v est dans la boule de u, alors on route directement vers v
  2. si v est un landmark, alors on route directement vers v
  3. sinon, on route vers le plus proche sommet de la boule u dont
     la couleur = hash(v). Une fois en ce sommet on route via la
     meilleure des deux options:
     - via un des landmarks (via le nca w)
     - via une boule contiguÃ«
  Le point 3 nÃ©cessite la gestion d'un en-tÃªte.
*/
{
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivÃ©
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // v dans B[u] ? (NB: B[u] existe toujours)
  
  i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
  // si v dans B[u]: u -> v
  if(i>=0){ DEBUG(printf("v is in B[u]\n");); return X->B[u]->dist[i]; }
  
  // v n'est pas dans B[u]
  // est-ce que v est un landmark ? (dans la table des landmarks)

  if(X->C[v]==0){ // 0=couleur des landmarks
    DEBUG(printf("%i is a landmark\n",v););
    return X->S[v]->D[u]; // dist(u,v), v=landmark
  }

  // v n'est pas landmark
  // est-on arrivÃ© au responsable de v ?

  if((s<0)&&(X->C[u]!=X->H[v])){ // si pas encore au responsable
    u=X->W[u]->node[X->H[v]]; // next-hop de u vers le responsable de v
    DEBUG(
	  printf("go to %i the next-hop to the closest node of color %i, the hash of %i",
		 u,X->H[v],v);
	  );
    return 1 + rs_agmnt_length_rec(u,v,s,w,X);
  }

  // on est passÃ© par (ou on est sur) le responsable de v
  // est-ce qu'on vient d'arrivÃ© sur le responsable ?
  
  if(s<0){ // on arrive au responsable de v
    // ici v n'est pas landmark ni dans la boule de u
    // cherche i tq v=CONT[u][i].v (cherche binaire, il doit y Ãªtre)
    contigue cont; cont.v=v;
    contigue *p=bsearch(&cont,X->CONT[u],X->F[X->C[u]],sizeof(contigue),fcmp_contigue);
    DEBUG(
	  if(p==NULL) return FAIL_ROUTING; // ne devrait pas arriver car v est dans CONT[u]
	  );
    cont=*p; // ici cont.v=v
    if(cont.s<0){
      DEBUG(printf("routage via contigue-ball, d=%i\n",cont.d););
      return cont.d; // routage via boule contiguÃ«
    }
    // routage via landmark s et nca w
    DEBUG(printf("routage via landmark %i and nca %i",cont.s,cont.w););
    return rs_agmnt_length_rec(u,v,cont.s,cont.w,X);
  }

  // on est passÃ© par le responsable de v
  // on monte vers le pÃ¨re du landmark, sauf si on arrive au nca w
  if(u==w){ // w=s est possible
    DEBUG(printf("arrive at nca %i of landmark tree %i\n",w,s););
    return X->S[s]->D[v] - X->S[s]->D[w]; // distance entre w et v dans T_s
  }
  u=X->S[s]->P[u]; // u=pÃ¨re(u), NB: u<>root car u<>w=nca
  DEBUG(printf("routing to the parent of u in the landmark tree %i",s););
  return 1 + rs_agmnt_length_rec(u,v,s,w,X);
}


int rs_agmnt_length(int u,int v,rs_dcr_tables *X)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_dcr() pour AGMNT pour router un message de u Ã  v, ou bien -1 si
  la route n'a pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables
  nÃ©cessaires au schÃ©ma. Si u<0, on rÃ©alise quelques tests de bases
  sur les tables X. La fonction fait essentiellement appel Ã  une
  fonction recursive oÃ¹ la gestion d'un en-tÃªte est nÃ©cessaire.
*/
{
  if(u<0){
    if(rs_dcr_length(-1,0,X)==1) return 1;
    if(X->CONT==NULL) return 1;
    if(X->F==NULL) return 1;
    return 0;
  }

  return rs_agmnt_length_rec(u,v,-1,-1,X);
}


rs_tzrplg_tables *rs_tzrplg(graph* const G,double t)
/*
  Calcule les tables de routage selon le schÃ©ma tz_rplg pour le graphe
  G, une adaptation du schÃ©ma de Thorup et Zwick avec les landmarks
  sur les sommets de plus haut degrÃ©.

  Cet algorithme est spÃ©cialisÃ© pour les graphes de type RPLG(n,t). Le
  paramÃ¨tre t attendu (power-law exponent) est celui du graphe G
  fournit. Les performances de ce schÃ©ma sont meilleures si le
  paramÃ¨tre t est le bon, mais marche quel que soit t>1.5. Les valeurs
  de t et de VARIANT permettent de calculer le nombre de landmarks.
  L'ensemble des landmarks est aussi nommÃ© "core" du graphe par Chung
  & Lu.

  On calcule deux tables: B et L. La table B[u], pour chaque sommet u,
  contient tous les sommets strictement plus proche que son plus
  proche landmark. La table L[u], pour chaque sommet u, contient le
  next-hop vers de u vers le landmark le plus proche de u selon
  l'arbre de +cc enracinÃ© dans le landmark.

  Rem: dans [CSTW12], pour les BC graphs, pour les moyennes ils
  prennent probablement n=10K et non n=+grande composante connexe
  (7K), ce qui change les choses.
*/
{
  printf("\nTZ RPLG\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  int const n=G->n;
  double gamma; /* gamma = paramÃ¨tre fonction de t */
  int core_size; /* core_size = nombre de landmarks */
  int i,u;

  /*** Calcul des landmarks ***/
  /* C[i]=liste des landmarks, i=0..core_size-1 */

  // calcule D = liste triÃ©e (ordre croissant) des degrÃ©s des sommets
  // ou permutation alÃ©atoire des sommets (si VARIANT=2)
  int *D;
  if(VARIANT<2) D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi);
  if(VARIANT==2){ // permutation alÃ©atoire
    ALLOCZ(D,n,_i);
    Permute(D,n);
  }

  // calcule core_size, dÃ©pend de VARIANT et de t
  if((VARIANT==0)&&(t>1.5)){
    gamma=(double)(t-2.0)/(double)(2.0*t-3.0);
    core_size=ceil(pow(n,gamma));
  }
  if((VARIANT==1)&&(t>1.5)){
    // C = { u : deg(u)>(n^gamma')/4 }, gamma'=(1-gamma)/(t-1)=1/(2t-3) 
    gamma=1.0/(double)(2.0*t-3.0);
    u=ceil(pow(n,gamma));
    // cherche le 1er sommet de degrÃ© <= u
    i=n-1; // part de la fin (haut degrÃ©)
    while((i>=0)&&(G->d[D[i--]]>u));
    core_size=i/4;
  }
  if((VARIANT==2)&&(t>0)) core_size=(int)t;
  if(t<0) core_size=(int)(-t);
  if(t==0) core_size=ceil(sqrt(n));

  core_size=max(core_size,0);    // core_size >= 0
  core_size=min(core_size,n); // core_size <= n

  // C=liste des landmarks
  NALLOCZ(int,C,core_size,D[n-_i-1]); // de +haut degrÃ© (ou sommets alÃ©atoires)
  
  // affiche C et sa densitÃ©
  for(i=u=0;i<core_size;i++) u += G->d[C[i]];
  printf("- core degree: ");
  APERCU(G->d[D[n-_i-1]],core_size,10,2);
  free(D); /* ne sert plus Ã  rien */
  printf("- core size: %i",core_size);
  if((VARIANT==0)&&(t>1.5)) printf(" (n^%g)",gamma);
  if((VARIANT==1)&&(t>1.5)) printf(" (deg>n^%g)",gamma);
  if(t==0) printf(" (sqrt(n))");
  printf("\n");
  printf("- sum of cluster's degrees: %i",u);
  if(u*100>n) printf(" (%.2lfn)",(double)u/(double)n);
  printf("\n- time to construct core: %s\n",TopChrono(1));
  
  /*** Construction des tables L ***/

  // L[u]->node[i] = identifiant du i-Ã¨me landmark
  // L[u]->dist[i] = distance entre u et L[u]->node[i]
  // lmin[u] = distance entre u et son plus proche landmark
  // label[u] = indice i du landmark le plus proche de u

  NALLOCZ(table*,L,n,new_table(core_size)); // n tables L de taille core_size
  NALLOCZ(int,lmin,n,n); // par dÃ©faut lmin[u]=n
  NALLOC(int,label,n);
  
  int l;
  param_bfs *X;
  
  for(i=0;i<core_size;i++){ // on fait un bfs() pour chaque landmark
    l=C[i]; // l=i-Ã¨me landmark
    X=bfs(G,l,NULL); // bfs() depuis l
    for(u=0;u<n;u++){
      L[u]->node[i]=l;
      L[u]->dist[i]=X->D[u];  // dist(u,l)
      if(X->D[u]<lmin[u]){ // landmark plus proche ?
	lmin[u]=X->D[u]; // NB: le rayon de B[u] sera lmin[u]-1
	label[u]=i; /* NB: L[u]->node[label[u]]=landmark le plus proche de u */ 
      }
    }
    free_param_bfs(X);
  }
  free(C);
  printf("- time to construct tables L: %s\n",TopChrono(1));
  
  /*** Construction des tables B ***/
  
  NALLOCZ(table*,B,n,NULL); // tableau de n tables B (vides au dÃ©part)
  X=new_param_bfs(); // pour bfs() depuis u
  X->clean=1; // initialisation complÃ¨te des distances, puis partielle
  
  /* construit la boule B[u] de rayon lmin[u]-1 avec bfs(u,.) */

  for(u=0;u<n;u++){
    X->hmax=max(0,lmin[u]-1); // la boule contient les sommets
			      // strictement plus proche que le coeur
    B[u]=new_table(0);
    // si B n'est pas vide, alors on fait un BFS pour calculer B
    if(X->hmax){
      bfs(G,u,X);
      B[u]->n=(X->n)-1; // taille de B[u] (sans u)
      /* copie les sommets (et leur distance) du bfs() en supprimant u */
      ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* =i-Ã¨me sommet de B[u] */
      ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* =distance du i-Ã¨me Ã  u */
    }
    // si B est vide alors:
    else B[u]->n=0;
  }
  free_param_bfs(X);
  free(lmin); /* ne sert plus Ã  rien */
  
  printf("- time to construct tables B: %s\n",TopChrono(1));
  
  /* tri des tables B */
  for(u=0;u<n;u++) if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;

  /* taille totale des tables */
  NALLOCZ(int,Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(L[u]) Z[u] += L[u]->n;
  }
  
  /* affiche taille min/max et moyenne/Ã©cart type des diffÃ©rentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(L[_i]->n,n,L[_i],"table L size");

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus Ã  rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_tzrplg_tables,RT,1);
  RT->B=B;
  RT->L=L;
  RT->label=label;
  RT->n=n; // besoin pour libÃ©rÃ©r les n tables
  
  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_tzrplg_length(int u,int v,rs_tzrplg_tables *X)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_tzrplg() pour router un message de u Ã  v, ou bien -1 si la route
  n'a pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables nÃ©cessaire au
  schÃ©ma, notamment les tables B et L. Si u<0 alors on teste la
  validitÃ© des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).
*/
{
  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->L==NULL) return 1;
    return 0;
  }

  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivÃ©
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }

  // routage dans la boule de u
  // v dans B[u] ?

  if(X->B[u]){
    int i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]: u -> v
    if(i>=0){
      DEBUG(
	    printf("in B[u], distance %i: (", X->B[u]->dist[i]);
	    int _v;
	    for(_v=0;_v<X->B[u]->n;_v++)
	      printf("%i,", X->B[u]->node[_v]);
	    printf(")\n");
	    );
      return X->B[u]->dist[i];
    }
  }

  // routage via le plus proche landmark de v
  // route: u -> lv -> v (aucun raccourci, voir Algo. 1 dans CSTW12)

  if(X->L[u]){
    // On doit rÃ©cupÃ©rer 2 entrÃ©es: L[u][lv] et L[v][lv]. Pour cela on
    // doit d'abord trouver l'identitÃ© de lv qui se trouve dans
    // l'Ã©tiquette de v:
    int lv = X->label[v];
    if(X->L[v])
        return X->L[u]->dist[lv] + X->L[v]->dist[lv];
  }

  // y'a un problÃ¨me, pas d'entrÃ©e pour v
  DEBUG(printf("fail: no route found when routing from %i to %i \n",u,v););
  return FAIL_ROUTING;
}


rs_bc_tables *rs_bc(graph* const G, int k)
/*
  SchÃ©ma de routage selon Brady-Cowen 2006.  Le stretch est additif
  est <= 2k. En particulier, si k=0, il s'agit d'un routage de plus
  court chemin.

  Principe: On construit un arbre BFS (=T) enracinÃ© dans le sommet de
  plus haut degrÃ© (=center). Le coeur (=C) est la boule de rayon k
  depuis la racine de T. Dans l'article d'origine, k=d/2 avec d pair.

  On construit une liste (=L) de BFS couvrant G ainsi qu'une forÃªt
  (=H) de BFS de G comme suit. Au dÃ©part, L={T}, et H est la forÃªt
  T\C. Puis, pour chaque arÃªte {u,v} de G\C\T, on vÃ©rifie si l'ajoÃ»t
  de {u,v} Ã  H crÃ©e un cycle ou pas. Si on ne crÃ©e pas de cycle, on
  met Ã  jour la forÃªt H en lui ajoutant {u,v}. Si on crÃ©e un cycle, on
  calcule un BFS de G de racine u (ou v) qu'on ajoute Ã  L (on favorise
  le sommet de plus grand degrÃ© sans les arÃªtes de T).

  Une fois toutes les arÃªtes {u,v} ainsi balayÃ©es, on cacule pour
  chaque composantes connexes de H un BFS. On obtient une forÃªt
  couvrante qu'on ajoute Ã  L.

  L'algorithme de routage de u Ã  v consiste simplement Ã  router dans
  l'arbre A de L contenant u et v et qui minimise dist_A(u,v).
*/
{
  int const n=G->n;
  int u,v,center,x,y,d,i;

  printf("\nBRADY-COWEN\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  /* trouve un sommet center de degrÃ© max, v=deg(center) */
  for(u=v=center=0;u<n;u++) if(G->d[u]>v) v=G->d[center=u];
  printf("- degree of the center: %i (id:%i)\n",v,center);

  /* calcule l'arbre T (bfs) */

  // on construit l'arbre en deux temps, d'abord on rÃ©cupÃ¨re le coeur
  // de rayon k, puis on poursuit en couvrant tout le graphe

  param_bfs* T_bfs=new_param_bfs(); // arbre de racine center
  T_bfs->clean=1;
  T_bfs->hmax=k;
  bfs(G,center,T_bfs);
  int core_size=T_bfs->n; // taille du coeur
  NALLOCZ(int,C,core_size,T_bfs->file[_i]); // C=liste des sommets du coeur
  printf("- core: "); APERCU(C[_i],T_bfs->n,10,2);
  printf("- core size: %i\n",core_size);

  // on finit la construction de T
  T_bfs->cont=1; 
  T_bfs->hmax=-1; // le bfs() va jusqu'au bout
  bfs(G,center,T_bfs); // T_bfs=bfs() de tout G
  printf("- excentricity of the center: %i\n",T_bfs->radius);
  printf("- time to construct the core and T: %s\n",TopChrono(1));

  // S[u]=-2 ssi u est dans le coeur, -1 sinon
  NALLOCZ(int,S,n,-1);
  for(i=0;i<core_size;i++) S[C[i]]=-2;
  free(C);

  /* calcule le graphe H */

  // au dÃ©but H=T\C, ajoute Ã  H les arÃªte de T
  graph* H=new_subgraph(G); // on crÃ©er un graphe vide H de la taille de G
  for(u=0;u<n;u++){
    x=T_bfs->P[u]; // x=pÃ¨re de u dans T
    // ajoute {u,x} Ã  H si u et son pÃ¨re x ne sont pas dans le coeur
    if((x>=0)&&(S[x]!=-2)&&(S[u]!=-2)) ADD_EDGE(H,u,x);
  }

  /* calcule les composantes connexes de H, via un DFS */
  param_dfs* H_dfs=new_param_dfs(n);
  ALLOCZ(H_dfs->C,n,S[_i]); // copie S dans ->C pour les sommets interdits
  H_dfs=dfs(H,0,H_dfs); // DFS depuis un sommet arbitraire
  int const nc=H_dfs->nc; // sauvegarde le nombre de composantes de H
  printf("- time to construct and to traverse (dfs) T\\C: %s\n",TopChrono(1));
  
  /* convertit les composantes en "reprÃ©sentant" pour FindSet() */
  // au dÃ©part couleur[u] est un indice
  // Ã  la fin couleur[u] est un sommet source de dfs sur H

  int *couleur=H_dfs->C; // couleur[u]=componsante du sommet u de H
  for(u=0;u<n;u++){
    if(S[u]==-2) continue;
    couleur[u]=H_dfs->R[couleur[u]];
  }

  /*** calcule la liste L des BFS ***/

  // pour toutes arÃªtes uv de G, ajouter uv Ã  H si cela ne crÃ©e pas de
  // cycle, sinon on calcule un BFS enracinÃ© en u ou v (deg max dans
  // G\H) qu'on ajoute Ã  L

  NALLOCZ(int,rang,n,0); // pour l'heuristic de FindSet()
  NALLOC(param_bfs*,L,n); // L=liste des BFS, L[i]=i-Ã¨me BFS
  NALLOCZ(param_bfs*,L_bfs,n,NULL); // L_bfs[u]=BFS pour le sommet u
  L[0]=L_bfs[center]=T_bfs; // on met T_bfs dans L et L_bfs
  int nbfs=1; // nombre de BFS dans la liste L
  int no_cycle=0; // nombre d'arÃªtes ne crÃ©eant pas de cycle

  for(u=0;u<n;u++){ // pour tous les sommets u de G
    if(S[u]==-2) continue; // ne rien faire si u est dans le coeur
    d=G->d[u]; // d=degrÃ©(u) dans G
    for(i=0;i<d;i++){
      v=G->L[u][i]; // v=i-Ã¨me voisin de u
      if(S[v]==-2) continue; // // ne rien faire si v est dans le coeur
      if((T_bfs->P[u]==v)||(T_bfs->P[v]==u)) continue; // ne rien faire si {u,v} dans E(T)
      if(u>=v) continue;
      // ici u<v
      x=UF_Find(u,couleur); // x=reprÃ©sentant de u
      y=UF_Find(v,couleur); // y=reprÃ©sentant de v
      if(x==y){ // si mÃªme reprÃ©sentant alors on fait un BFS depuis u ou v
	if((L_bfs[u])||(L_bfs[v])) continue; // BFS dÃ©jÃ  calculÃ©
	// ici, on a jamais fait de BFS ni depuis u ni depuis v
	// on enracine le BFS en u ou v selon le degrÃ© de G-H
	x=(G->d[u]-H->d[u] > G->d[v]-H->d[v])? u : v;
	L[nbfs++]=L_bfs[x]=bfs(G,x,NULL); // met Ã  jour L et L_bfs
      }
      else{ // si pas mÃªme couleur alors pas de cycle
	ADD_EDGE(H,u,v);
	no_cycle++;
	UF_Union(x,y,couleur,rang); // fusion pour UF_Find()
      }
    }
  }

  free(rang);
  x=no_cycle+nbfs-1; // x=nb d'arÃªtes dans G\C\T
  printf("- #edges in C: %i\n",nb_edges(G)-x-n+core_size);
  printf("- #edges in G\\C\\T: %i\n",x);
  printf("- #edges added to T\\C to make H: %i\n",no_cycle);
  printf("- #bfs trees computed from edges not in H: %i\n",nbfs-1);
  printf("- time to construct H and these bfs trees: %s\n",TopChrono(1));

  /*** calcule un BFS pour chaque composante de H ***/

  param_bfs* H_bfs=new_param_bfs();
  H_bfs->clean=0; // trÃ¨s important pour faire l'union de BFS (ne pas initialiser H_bfs->D)
  H_bfs->D=S; // pour sÃ©lectionner les sommets de G\S (surtout ne pas faire free(S) ...)
  
  // lance un BFS depuis chaque composante de H, Ã  partir des racines
  // R[i] (du dfs de H) qui sont leur propre reprÃ©sentant
  for(i=x=0;i<nc;i++){ // x=nombre de BFS calculÃ©s dans H
    u=H_dfs->R[i];
    if(u==UF_Find(u,couleur)){
      bfs(H,u,H_bfs); // superpose les BFS de H
      x++;
    }
  }
  printf("- #components of T\\C: %i\n",nc);
  printf("- #components of H: %i\n",x);
  printf("- time to construct the bfs forest for H: %s\n",TopChrono(1));
  free_param_dfs(H_dfs);
  free_graph(H);
  BARRE;
  
  /* taille totale des tables */
  NALLOC(int,Z,n); /* Z[u]=taille de la table de u, doit Ãªtre au moins 1 */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    Z[u]=nbfs+1; // doit Ãªtre au moins 1
    if(S[u]!=-2) Z[u]++; // si pas dans le coeur (Ã  cause de H)
  }
  
  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus Ã  rien */
  BARRE;
  
  /* transforme H en vÃ©ritable forÃªt couvrante (les champs ->P et ->D
     des sommets du coeur n'Ã©tant pas forcÃ©ment cohÃ©rent), puis
     l'ajoute Ã  la liste L. Il faut bien sÃ»r que la forÃªt ne soit pas
     vide. */

  if(nc){ // si H possÃ¨de au moins un sommet, alors nc>0
    // corrige H_bfs pour qu'il soit une forÃªt couvrante, sinon
    // dist_bfs() dans _length() peuvent ne pas marcher
    // puis ajoute H_bfs Ã  L
    for(u=0;u<n;u++)
      if(S[u]==-2){
	H_bfs->P[u]=-1;
	H_bfs->D[u]=0;
      }
    L[nbfs++]=H_bfs;
  }
  REALLOC(L,nbfs); // rÃ©duit L

  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_bc_tables,RT,1);
  RT->L=L;
  RT->Lu=L_bfs;
  RT->nbfs=nbfs;
  ALLOCZ(RT->dist,n,L_bfs[_i]?L_bfs[_i]->D:NULL); // distances partielles

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_bc_length(int u,int v,rs_bc_tables *X)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_bc() pour router un message de u Ã  v, ou bien -1 si la route n'a
  pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables nÃ©cessaire au
  schÃ©ma, notamment les tables ... Si u<0 alors on teste la validitÃ©
  des tables (et on renvoie une valeur non-nulle en cas d'erreur).

  Rem: on pourrait obtenir des routes plus courtes en faisant du
  pas-Ã -pas plutÃ´t que de router dans l'arbre.
*/
{
  if(u<0){
    if(X==NULL) return 1;
    if(X->L==NULL) return 1;
    if(X->Lu==NULL) return 1;
    return 0;
  }

  DEBUG(printf("Routing from %i to %i\n",u,v););

  // 1) on est arrivÃ©
  if(u==v) return 0;

  // 2) est-ce que u ou v est la racine d'un des BFS de L ?
  if(X->Lu[u]) return X->Lu[u]->D[v];
  if(X->Lu[v]) return X->Lu[v]->D[u];
  DEBUG(
	printf("ni u ni v ne sont racines d'un BFS\n");
	PRINT(X->nbfs);
	);

  // 3) calcule la distance dans entre u et v dans le meilleur arbre de L
  int d,i,d0=INT_MAX; // d0=distance min recherchÃ©e
  for(i=0;i<X->nbfs;i++){ // parcourt les arbres de L
    d=dist_bfs(u,v,X->L[i]); // d=dist(u,v) dans le i-Ã¨me arbre de L
    if(d>=0) d0=min(d0,d); // met Ã  jour la distance min
  }
  
  if(d0<INT_MAX) return d0;
  return FAIL_ROUTING;
}


rs_hdlbr_tables *rs_hdlbr(graph* const G,int k)
/*
  Calcule les tables de routage selon le schÃ©ma HDLBR pour le graphe
  G, une adaptation du schÃ©ma tz_rplg en name-independant.

  Cet algorithme est spÃ©cialisÃ© pour les graphes de type RPLG(n,t).
  L'algorithme HDLBR, Ã  l'instar de TZ_RPLG, est paramÃ©trÃ© par dÃ©faut
  Ã  un nombre de landmarks fixÃ© Ã  k=n^x=n^1/2:
  
    " We set x=1/2 because this setting minimizes the storage overhead, 
    making routing table size of landmarks and average routing table 
    size of non-landmark nodes both bounded by Ã•(n^1/2) bits. "
*/
{
  printf("\nHDLBR\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  int const n=G->n;
  int i,u,v;

  /*** Calcul des landmarks ***/
  // Core=liste des landmarks, les k sommets de plus haut degrÃ©

  // rÃ©cupÃ¨re la liste des identifiants des noeuds triÃ©e par degrÃ©
  int *D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi); /* D=liste triÃ©e par degrÃ©s croissant des noeuds */

  // Core=liste des landmarks=liste des sommets de plus haut degrÃ©
  NALLOCZ(int,Core,k,D[n-_i-1]);

  // Trie le coeur par identifiants (utile pour setSearch lors du routage)
  QSORT(Core,k,fcmp_int);
  printf("- time to construct the Core: %s\n",TopChrono(1));
  BARRE;

  free(D);

  printf("- #landmarks: %i",k);
  printf(" (ceil{âˆšn}=%i)\n",(int)ceil(sqrt(n)));
  printf("- landmark list: ");
  APERCU(Core[_i],k,10,2);
  printf("- landmark degrees: ");
  APERCU(G->d[Core[_i]],k,10,2);
  for(i=u=0;i<k;i++) u += G->d[Core[i]]; // somme des degrÃ©s des landmarks
  printf("- average landmark's degree: %.2lf\n",(double)u/(double)n);
  printf("- time to construct landmarks: %s\n",TopChrono(1));
  
  /*** Construction des tables L ***/
  // L[u] est dÃ©finie pour tout les sommets u
  // L[u]->node[i]=parent du sommet u dans l'arbre BFS du i-Ã¨me landmark
  // L[u]->vpd[i]=indice de l(u), le plus proche landmark de u
  
  NALLOCZ(table*,L,n,new_table(k)); // n tables de taille k>0
  // par dÃ©faut L[u]->vpd=-1, il faut l'initialiser Ã  0
  for(u=0;u<n;u++) {
    L[u]->vpd=0; // L[u]->vpd=indice du landmark le plus proche de u
    L[u]->n=k;
    ALLOCZ(L[u]->node,k,-1); 
    ALLOCZ(L[u]->dist,k,-1); 
  }
  param_bfs *X=new_param_bfs(); // structure utilisÃ©e pour le rÃ©sultat du BFS depuis u
  X->clean=1; // initialisation complÃ¨te des distances, puis partielle

  for(i=0;i<k;i++){
    v=Core[i]; // v=landmark numÃ©ro i
    bfs(G,v,X);
    for(u=0;u<n;u++){
      L[u]->node[i]=X->P[u]; // pÃ¨re de u dans l'arbre de racine v   
      L[u]->dist[i]=X->D[u]; // dist(u,v)
      if(X->D[u]<L[u]->dist[L[u]->vpd]) L[u]->vpd=i; // met Ã  jour L[u]->vpd
    }
  }
  free_param_bfs(X);
  printf("- time to construct tables L: %s\n",TopChrono(1));

  // calcule la distance inter-landmark
  for(i=u=0;u<k;u++) // i=distance inter-landmark
    for(v=u+1;v<k;v++)
      i=max(i,L[Core[u]]->dist[v]);
  printf("- inter-landmark distance: %i\n",i);
  printf("- time to compute this distance: %s\n",TopChrono(1));
    

  /*** Construction des tables B ***/
  // B[u] est dÃ©finie pour tout sommet u
  // B[u]->node=liste des sommets Ã  distance < dist(u,l(u))
  // B[u]->dist=distance de ces sommets Ã  u
  // B[u] contient en plus l(u) qui est toujours B[u]->node[0]
  // Si u est un landmark, alors B[u]={u}, sinon u n'est pas
  // stockÃ© dans B[u].
  
  NALLOCZ(table*,B,n,new_table(0)); // n tables vides au dÃ©part
  X=new_param_bfs(); // structure utilisÃ©e pour le rÃ©sultat du bfs depuis u
  X->clean=1; // initialisation complÃ¨te des distances, puis partielle
  
  for(u=0;u<n;u++){
    X->hmax=max(0,L[u]->dist[L[u]->vpd]-1); // hmax=dist(u,l(u))-1
    bfs(G,u,X);
    B[u]->n=X->n; // NB: ici X->node[0]={u}
    // si u est landmark alors B[u]={u}, sinon on remplace dans B[u]
    // le sommet u par l(u). Dans tous les cas, |B[u]| = X->n.
    ALLOCZ(B[u]->node,B[u]->n,X->file[_i]); // B[u]->node[i]=i-Ã¨me sommet de B[u]
    ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); // B[u]->dist[i]=distance du i-Ã¨me Ã  u
    i=L[u]->vpd; // i=indice de l(u)
    if(Core[i]!=u){ // si u n'est pas un landmark
      B[u]->node[0]=Core[i];
      B[u]->dist[0]=L[u]->dist[i];
    }
    //free(L[u]->dist); // ne sert plus Ã  rien
    //L[u]->dist=NULL;
  }
  free(X);
  printf("- time to construct tables B: %s\n",TopChrono(1));
  
  /*** Calcul de la taille des boules inverse ***/
  // D[u]=taille de la boule inverse de u, c'est-Ã -dire de B[u] mais
  // sans l(u).

  ALLOCZ(D,n,0); // compteurs Ã  0
  for(u=0;u<n;u++)
    for(i=1;i<B[u]->n;++i) // i>0, ne compte pas B[u]->node[0]=l(u)
      D[B[u]->node[i]]++;

  printf("- time to compute inverse ball sizes: %s\n",TopChrono(1));

  /* H[u]=0..k-1, hash du sommet u */
  int *H=MakeHash(NULL,n,k,HASH);

  /*
    On ne construit pas la table couleur car elle n'est pas nÃ©cessaire
    pour le routage, en distribuÃ© elle sert Ã  faire la "traduction" du
    schÃ©ma Ã©tiquettÃ© vers le schÃ©ma name-independant. Ici cette table
    est un entier qui compte le nombre d'entrÃ©es qu'elle aurait dÃ»
    avoir.
  */

  int *F;
  FREQMINMAX(F,k,H,n,"the hash"); /* frÃ©quence des hashs (sert pour taille des tables) */
  printf("- time to compute hash values: %s\n",TopChrono(1));

  /* tri des tables B */
  for(u=0;u<n;u++)
    if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;

  /* taille totale des tables */  
  NALLOCZ(int,Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n-1; /* -1 car le landmark l(u) est toujours dans B[u] */
    if(L[u]) Z[u] += L[u]->n;
    Z[u] += D[u]; /* vraie taille des boules inverses */ 
    i=L[u]->vpd;
    if(Core[i]==u) Z[u] += F[i]; // si u est un landmark
  }
  
  /* affiche taille min/max et moyenne/Ã©cart type des diffÃ©rentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(D[_i],n,1,"inverse ball size");
  MINMAXMOY(L[_i]->n,n,L[_i],"table L size");
  MINMAXMOY(F[_i],k,F[_i],"table Color size");
  free(D);
  free(F);

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus Ã  rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_hdlbr_tables,RT,1);
  RT->B=B;
  RT->L=L;
  RT->Core=Core;
  RT->core_size = k;
  RT->H=H;
  RT->n=n; // besoin pour libÃ©rÃ©r les n tables
  
  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_hdlbr_length_rec(int u,int v,rs_hdlbr_tables *X, int lv)
/*
  Renvoie le nombre de sauts du routage selon les tables gÃ©nÃ©rÃ©es par
  rs_hdlbr() pour router un message de u Ã  v, ou bien -1 si la route
  n'a pu Ãªtre dÃ©terminÃ©e. Dans X on a toutes les tables nÃ©cessaire au
  schÃ©ma, notamment les tables B, C et L. Si u<0 alors on teste la
  validitÃ© des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).

  Le paramÃ¨tre lv est l'indice du plus proche landmark de, ou -1 s'il
  n'a pas pu Ãªtre dÃ©terminÃ©.

  L'algorithme de routage HDLBR de u Ã  v est le suivant: (comme on est 
  en name-independant on doit parfois trouver l'Ã©tiquette de v qui se 
  trouve sur le landmark de couleur h(v)).
    1) si u=v alors on est arrivÃ©
    2) si v appartient Ã  B[u], u appartient Ã  B[v] ou v est un landmark, router vers v en plus court chemin
    3) si lv<0, alors router vers le landmark de couleur h(v)
    4) si u=lv router vers v avec la boule inverse Ã©tendue de u
    5) sinon, router vers le landmark de v (lv)
*/
{
  DEBUG(printf("Routing from %i to %i \n", u,v););

  // 1) on est arrivÃ©
  if(u==v) return 0;

  /** Routage dans la boule/boule inverse de u **/
  // 2) v dans B[u] ?
  if(X->B[u]){
    int i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]: u -> v
    if(i>=0) {
      DEBUG(printf("\t - Direct routing %i in B[%i] (distance %i) \n", v, u,X->B[u]->dist[i]););
      return X->B[u]->dist[i];
    }
  }
  // 2) u dans B[v] (ie., v dans Bâ»Â¹[u]) ?
  if(X->B[v]){
    int i=SetSearch(u,X->B[v]->node,X->B[v]->n,1);
    // si u dans B[v]: u -> v
    if(i>=0) {
      DEBUG(printf("\t - Direct routing %i in C[%i] (distance %i) \n", v, u,X->B[v]->dist[i]););
      return X->B[v]->dist[i];
    }
  }
  // 2) v est un landmark ?
  if(X->L[u]){
    int i=SetSearch(v,X->Core,X->core_size,1);
    if (i>=0) {
      DEBUG(printf("\t - Direct routing %i is a landmark (distance %i)  \n", v, X->L[u]->dist[i]););
      return X->L[u]->dist[i];
    }
  }

  // 3) Si le header est vide (lv<0), alors on est pas encore passÃ©
  // par h(v) et on veut router un pas vers h(v)
  if(lv<0){
    int hash = X->H[v]; // indice du landmark dans le tableau Core
    int hv = X->Core[hash]; // landmark responsable de v
    // si on est arrivÃ© en h(v), routage vers le plus proche landmark de v
    if(u==hv) {
      DEBUG(printf("\t - Manager hv=%i reached, now going toward lv=%i \n", hv, X->Core[X->L[v]->vpd]););
      return rs_hdlbr_length_rec(u, v, X, X->L[v]->vpd);
    }

    DEBUG(printf("\t - Routing toward v's manager (hv=%i) \n", hv););
    // sinon on monte vers le pÃ¨re de u dans l'arbre de racine h(v)
    return 1 + rs_hdlbr_length_rec(X->L[u]->node[hash],v,X,-1);
  }

  // 4) routage via le plus proche landmark de v
  if(lv>=0){
    DEBUG(printf("\t - Routing toward lv=%i \n", X->Core[lv]););
    return 1 + rs_hdlbr_length_rec(X->L[u]->node[lv],v,X,lv);
  }
  // y'a un problÃ¨me, pas d'entrÃ©e pour v
  DEBUG(printf("fail: no route found when routing from %i to %i \n",u,v););
  return FAIL_ROUTING;
}


int rs_hdlbr_length(int u,int v,rs_hdlbr_tables *X){

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->L==NULL) return 1;
    if(X->Core==NULL) return 1;
    if(X->H==NULL) return 1;
    return 0;
  }

  return rs_hdlbr_length_rec(u,v,X,-1); // header vide au dÃ©part
}


enum{
  SC_NONE,
  SC_ALL,
  SC_ONE,
  SC_NPAIRS,
  SC_PAIR,
  SC_EDGES,
  SC_UV,
  SC_UNTIL,
};


static inline long route_uv(const graph* const G,
			    int const u,int const v,
			    int const h,const int hmax,
			    int **dist,long **stat,
			    int* const L,int* const M,param_bfs *X)
/*
  Remplit la table de statistiques stat[][] et aussi de distance
  dist[][] des sommets testÃ©s avec h=longueur du routage de u Ã  v
  (dÃ©pend aussi de SCENARIO.dist). Dans le scenario SC_EDGES les
  distances ne sont pas mise Ã  jour (la distance Ã©tant toujours
  1). Renvoie une valeur <0 si une erreur est survenue ou bien la
  distance entre u et v sinon. Le tableau n'est pas remplit si une
  erreur est survenue, en particulier si h>hmax.

  On pourrait calculer ici h=length(u,v,X), plutÃ´t que de le passer
  comme paramÃ¨tre. On pourrait ainsi calculer avant dist(u,v), et
  ensuite la passer comme paramÃ¨tre Ã  length(u,v,X,d) qui pourrait
  Ã©ventuellement l'exploiter.
*/
{
  if((h<0)||(h>hmax)){
    DEBUG(printf("FAIL !!!  %i->%i: #hop=%i\n",u,v,h););
    return -1; /* y'a un problÃ¨me, sommet suivant */
  }

  int j,k,l,t;

  /* calcule k=dist(u,v) */
  if(SCENARIO.mode==SC_EDGES) k=1; /* pas de bfs() dans ce cas */
  else{ /* est-ce que dist(u,v) est dÃ©jÃ  connue ? */
    if((dist[u]==NULL)&&(dist[v]==NULL)){
      if(SCENARIO.dist){ /* on calcule et stocke dist[u] */
	bfs(G,u,X);
	ALLOCZ(dist[u],G->n,X->D[_i]); /* alloue dist[u] et y copie X->D[] */
      }else{ /* ici on ne stocke pas dist[u] */
	X->hmax=h; /* borne sup sur la distance de u Ã  v */
	bfs(G,u,X);
      }
      k=X->D[v]; /* k=dist(u,v) */
    }
    else k=dist[u]? dist[u][v] : dist[v][u]; /* k=dist(u,v) */
  }
  
  if(h<k){
    DEBUG(printf("FAIL !!!  %i->%i: #hop=%i dist(u,v)=%i\n",u,v,h,k););
    return -1; /* il y'a un problÃ¨me, trop court */
  }

  /* ici il faut faire +1 dans stat[k][j] oÃ¹ j=h-k */
  j=h-k;
  if(j>=L[k]){ /* tableau stat[k] n'est pas assez grand */
    l=L[k]; /* sauvegarde L[k] */
    L[k]=max(j+1,l)<<1; /* on double la taille courante ou de h-k+1 */
    if(stat[k]==NULL) ALLOCZ(stat[k],L[k],0L); /* premiÃ¨re fois */
    else{
      REALLOC(stat[k],L[k]); /* aggrandit stat[k] */
      for(t=l;t<L[k];t++) stat[k][t]=0L; /* initialise la fin du tableau */
    }
  }
  stat[k][j]++; /* une route de longueur h=k+j de plus */
  M[k]=max(M[k],j); /* dÃ©tour max pour la distance k */

  return k; /* tout c'est bien passÃ© */
}


int pgcd(int a,int b)
/*
  Renvoie le plus grand commun diviseur de a et b. NB: a et b premier
  entre eux ssi |pgcd(a,b)|=1.
*/
{
  while(a){ int c=a; a=b%a; b=c; }
  return b;
}


int routing_test(graph* const G,
		 void* const T,rt_length const length,
		 int hmax,int** const distp)
/*
  Teste un scenario de routage pour le graphe G avec les tables de
  routage T et la fonction de longueur length(), puis affiche les
  distributions des longueurs de route, de distance et de stretch. Le
  scenario est dÃ©crit par la variable globale SCENARIO. La matrice
  distp[][] est une matrice partielle de distances du graphe
  (Ã©ventuellement calculÃ©es lors de la construction de T), permettant
  d'accÃ©lÃ©rer les tests. En effet, pour calculer le stretch on doit
  calculer la distance.

  Plus prÃ©cisÃ©ment, pour tout sommet u de G, distp[u], si non NULL,
  doit Ãªtre la distance de u vers tous les sommets de G. Attention !
  un vecteur partiel de distance n'est pas autorisÃ©. Il ne faut pas
  que distp[u] soit construit Ã  l'aide de bfs() partiel par
  exemple. Il est cependant possible d'avoir distp[u]=NULL (car
  matrice partielle <> vecteur partiel). Si distp=NULL, alors aucune
  distance n'est prÃ©-calculÃ©es. Le tableau distp n'est pas modifiÃ©.
  Seules les vecteurs de distance calculÃ©s par routing_test() sont
  libÃ©rÃ©s Ã  la fin, charge Ã  l'appelant de libÃ©rer les vecteurs de
  distp.

  La fonction renvoie une valeur non-nulle si une erreur s'est
  produite, et 0 sinon. La valeur de hmax indique la longueur maximum
  de routage autorisÃ©e. Cela permet de dÃ©tecter des routes trop
  longues pour Ãªtre correctes (et Ã©viter une boucle infinie). Si
  hmax<0, alors on prend hmax=2n comme valeur par dÃ©faut.

  Pour trouver un petit graphe dont le stretch est > 3.00, on peut
  faire le script suivant:

  n=10; a=0; i=0; while [[ $(echo "$a <= 3.00" | bc) -eq 1 ]]; do a=$(./gengraph gabriel $n -seed $i -check routing scenario all agmnt -1 | grep -e "- maximum stretch:" | awk '{print $4}'); echo "n=$n, stretch=$a, seed=$i"; i=$((i+1)); done
*/
{
  if(SCENARIO.mode==SC_NONE) return 0;

  TopChrono(1);
  printf("\nROUTING TEST\n");
  BARRE;

  /* vÃ©rification de la table T */
  if(length(-1,0,T)) return 1;

  int const n=G->n;
  int u,v,i,j,k,t,h;
  int dmax; /* plus longue distance */
  int lmax; /* plus longue route */
  double x; /* pour les affichage de distributions */

  /* type "long" pour le comptage */
  long p,s,s2;
  long err=0L; /* nb de routage erronÃ©s */

  param_bfs *X=new_param_bfs(); /* pour le calcul de la distance */
  X->clean=1; /* pour appels multiples */

  /* dist[u][v]=distance entre u et v=0..n-1, n valeurs diffÃ©rentes au plus */
  /* on a dist[u]=NULL si on a pas encore calculÃ© de bfs(u,...) */
  NALLOCZ(int*,dist,n,distp?distp[_i]:NULL); /* dist[u]=distp[u] ou NULL */

  /* stat[k][j]=#routage entre sommets Ã  distance k et de longueur de k+j (j=dÃ©tour) */
  NALLOCZ(long*,stat,n,NULL); /* stat[k]=NULL si pas encore eut de stat pour dist k */
  /* L[k]=taille du tableau stat[k] */
  NALLOCZ(int,L,n,0); /* taille nulle par dÃ©faut */
  /* M[k]=dÃ©tour maximum (=longueur-k) rencontrÃ©e pour la distance k, M[k]<= L[k] */
  NALLOCZ(int,M,n,-1); /* par dÃ©faut M[k]<0 */
  if(hmax<0) hmax=n<<1;

  switch(SCENARIO.mode){
    
  case SC_ALL:
    p=(long)n*(long)(n-1); /* cast "(long)" important */
    if(2*lg(n)>(int)(8*sizeof(p)-1)) Erreur(33); /* dÃ©passement arithmÃ©tique, p est trop grand */
    printf("- all-to-all pairs: %s pairs ...",millier(p));
    fflush(stdout);
    for(u=0;u<n;u++)
      for(v=0;v<n;v++)
	if(u!=v) err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;

  case SC_UNTIL:
    p=s=(long)n*(long)(n-1); /* cast "(long)" important */
    if(2*lg(n)>(int)(8*sizeof(p)-1)) Erreur(33); /* dÃ©passement arithmÃ©tique, p est trop grand */
    printf("- all pairs until stretch â‰¥ %g ...",SCENARIO.stretch);
    fflush(stdout);
    t=randomu(n);
    for(i=0;i<n;i++)
      for(j=0;j<n;j++){
	if(i==j) continue;
	u=(i+t)%n,v=(j+t)%n;
	h=length(u,v,T);
	k=route_uv(G,u,v,h,hmax,dist,stat,L,M,X); /* k=dist(u,v) */
	if(k<0){ err++; continue; }
	if(h>=(k*SCENARIO.stretch)){ /* stop: stretch max atteint */
	  s=(long)i*(long)(n-1)+(long)(j+(j<i)); /* s=#paires testÃ©es */
	  i=j=n; /* pour arrÃªter les deux boucles. NB: aprÃ¨s i>n */
	}
      }
    break;

  case SC_PAIR:
  case SC_NPAIRS:
    x=(SCENARIO.mode==SC_NPAIRS);
    p=x?n:SCENARIO.u;
    printf("- %srandom pairs: %s pair%s ...",x?"n ":"",millier(p),PLURIEL(p));
    fflush(stdout);
    for(i=0;i<p;i++){
      u=randomu(n);
      v=randomu(n-1); // NB: n>=2 car au moins une arÃªte
      if(v>=u) v++; // u,v alÃ©atoires dans [0,n[ avec u<>v
      err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    }
    break;
  
  case SC_EDGES:
    p=nb_edges(G)<<1;
    printf("- all neighbor pairs: %s pair%s ...",millier(p),PLURIEL(p));
    fflush(stdout);
    for(u=0;u<n;u++){
      p=G->d[u];
      for(i=0;i<p;i++){
	v=G->L[u][i];
        err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
      }
    }
    break;
  
  case SC_ONE:
    u=SCENARIO.u;
    if(u>=n) Erreur(17);
    p=n-1;
    if(u<0) u=randomu(n);
    printf("- one-to-all pairs from %i: %s pair%s ...",u,millier(p),PLURIEL(p));
    fflush(stdout);
    for(v=0;v<n;v++)
      if(u!=v) err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;
  
  case SC_UV:
    u=SCENARIO.u;
    v=SCENARIO.v;
    if((u>=n)||(v>=n)) Erreur(17);
    if(u<0) u=randomu(n);
    if(v<0) v=randomu(n);
    printf("- routing from %i to %i: 1 pair ...",u,v);
    fflush(stdout);
    err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;

  default:
    Erreur(23);
  }

  printf(" %s (%s)\n",((SCENARIO.mode==SC_UNTIL)&&(i==n))?"not found":"Ok",TopChrono(1));
  if(SCENARIO.mode==SC_UNTIL){
    printf("- #tested routings: %s",millier(s));
    if(s==p) printf(" (all the pairs)");
    printf("\n");
    p=s; /* maintenant p=nombre de paires testÃ©es */
  }
  
  /* libÃ¨re les tableaux devenus inutiles: X et L */
  free_param_bfs(X);
  free(L);
  
  /* on ne libÃ¨re que les distances calculÃ©es par routing_test() */
  for(u=0;u<n;u++) if((distp==NULL)||(distp[u]==NULL)) free(dist[u]);
  free(dist);

  /*
    Affiche la distribution:
    - des longueurs de route,
    - des distances,
    - des stretch.
  */

  /* dÃ©termine dmax (=distance max) et lmax (=longueur routage max) */
  dmax=lmax=0;
  for(k=0;k<n;k++)
    if(stat[k]){
      dmax=max(dmax,k);
      lmax=max(lmax,k+M[k]);
  }

  /* pour un petit gain mÃ©moire: rÃ©ajuste les tableaux dÃ©pendant de la distance k */
  /* NB: le tableau L n'existe plus */
  REALLOC(stat,dmax+1); /* stat[0..dmax] */
  REALLOC(M,dmax+1);    /*    M[0..dmax] */

  /* calcule le nombre p de routages corrects (non erronÃ©s) */
  for(k=0,p=0L;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	p += stat[k][j];

  printf("- #failed routings: %s",millier(err));
  if(err) printf(" (%i%%)",p?(int)ceil((100.0*err)/(double)(p+err)):100);
  printf("\n");
  if(p==0L){
    printf("  all routings failed!\n");
    goto fin_routing_test;
  }
  /* NB: ici p>0 */
  
  BARRE;
  printf("- route length distribution:\n");

  /* calcule F[i]=nombre de routages de longueur i */
  NALLOCZ(long,F,lmax+1,0L);
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k+j] += stat[k][j];
  
  /* s=somme totale des longueurs, s2=somme du carrÃ© des longueurs */
  for(i=0,s=s2=0L;i<=lmax;i++) s+=((long)i)*F[i],s2+=((long)i)*((long)i)*F[i];

  for(i=0;i<=lmax;i++){
    if(F[i]==0L) continue;
    x=(double)F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [Ã— %li] \n",F[i]);
  }
  printf("- average route length: %.02lf Â± %.02lf (%li/%li)\n",s/(double)p,ECARTYPE(s,s2,p),s,p);
  printf("- maximum route length: %i\n",lmax);
  free(F);

  BARRE;
  printf("- distance distribution:\n");
  ALLOCZ(F,dmax+1,0L); /* F[k]=nombre de distances de longueur k */
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k] += stat[k][j];

  /* s=somme totale des distances, s2=somme du carrÃ© des distances */
  for(i=0,s=s2=0L;i<=dmax;i++) s+=((long)i)*F[i],s2+=((long)i)*((long)i)*F[i];

  for(i=0;i<=dmax;i++){
    if(F[i]==0) continue;
    x=F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [Ã— %li] \n",F[i]);
  }
  printf("- average distance: %.02lf Â± %.02lf (%li/%li)\n",s/(double)p,ECARTYPE(s,s2,p),s,p);
  printf("- maximum distance: %i\n",dmax);
  free(F);

  BARRE;
  printf("- stretch distribution:\n");

  /* t=borne sup sur le nombre de stretch diffÃ©rents */
  for(k=t=0;k<=dmax;k++)
    if(stat[k]) t += M[k]+1;

  triplet *ptr;
  triplet e={0,0,0L}; /* triplet nul */
  NALLOCZ(triplet,P,t,e); /* P=tableau de triplets (j,k,stat[k][j]) */

  /* on construit P */
  /* t=nombre de triplets ajoutÃ©s Ã  P */

  /* TODO: t pourrait Ãªtre trÃ¨s grand â‰ƒ n^2 et donc il faudrait qu'il
     soit de type long. Ensuite, si t est trÃ¨s grand, il ne faut pas
     tous les afficher, mais plutÃ´t une distribution (10 ranges par
     exemples). La mÃªme remarque s'applique Ã  la distribution des
     distances (cependant, il ne peut avoir que n-1 distances
     diffÃ©rentes). */

  for(k=t=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	if(stat[k][j]){
	  if(k){ u=pgcd(k,j); e.x=j/u; e.y=k/u; }
	  else{ e.x=0; e.y=1; } /* pour dist=0, le stretch est 1+0 */
	  e.z=stat[k][j];
	  ptr=bsearch(&e,P,t,sizeof(triplet),fcmp_stretch);
	  if(ptr) P[(int)(ptr-P)].z += e.z;  /* e est dÃ©jÃ  dans P, ajoute e.z */
	  else{
	    P[t++]=e; /* on ajoute une nouvelle entrÃ©e e Ã  P */
	    QSORT(P,t,fcmp_stretch); /* trie le nouveau tableau */
	  }
	}

  /* affiche les stretchs contenus dans P */
  for(i=0;i<t;i++){
    x=(double)(P[i].z)/(double)p;
    printf("  %0.3lf (%i/%i)",1.0+(double)P[i].x/(double)P[i].y,P[i].x+P[i].y,P[i].y);
    printf("\t%02i%% ",(int)(100.0*x));
    RULING(x);
    printf(" [Ã— %li] \n",P[i].z);
  }
  free(P);

  /* stretch moyen, Ã©cart type et stretch max */
  /* on pourrait l'avoir directement avec P[], mais l'avantage avec
     stat[][] est qu'on peut aussi avoir la distance max qui atteint
     le stretch max. On perd l'info dans P[] Ã  cause du pgcd). */

  double r=0,r2=0; /* r=somme total des stretch, r2=somme du carrÃ© des stretch */
  double smax=1; /* smax=stretch max */

  u=v=1; /* v=distance et u=dÃ©viation pour le stretch max */
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++){
	x=1; if(k) x += (double)j/(double)k; /* x=stretch, x=1.0 si k=0 */
	r += (double)stat[k][j]*x; r2 += (double)stat[k][j]*x*x;
	if(x>=smax) smax=x,u=j,v=k;
      }
  printf("- average stretch: %.03lf Â± %.03lf (%.02lf/%li)\n",(double)r/(double)p,ECARTYPE(r,r2,p),r,p);
  printf("- maximum stretch: %.02lf (%i/%i)\n",1.0+(double)u/(double)v,u+v,v);

  if(err){
    BARRE;
    printf("- Warning! there were %s failed routings\n",millier(err));
  }

 fin_routing_test:
  BARRE;
  free(M);
  FREE2(stat,dmax+1);
  return 0;
}


#undef ECARTYPE
#undef MINMAXMOY
#undef FREQMINMAX
#undef STR_DISTRIB
#undef LEN_DISTRIB
#undef RULING


/***********************************

       FIN ROUTINES POUR LES
          ROUTING SCHEMES

***********************************/


void RemoveVertex(graph* const G,int const u){
/*
  Supprime toutes les occurences du sommet u dans le graphe G non
  null, sans renumÃ©roter les sommets et sans rÃ©allouÃ©es les
  listes. AprÃ¨s cet appel, les listes de G ne sont plus forcÃ©ment
  triÃ©es, mais la liste de G->L[u] existe toujours. Simplement G->d[u]
  est Ã  zÃ©ro.

  Effet de bord: modifie G->sort, G->d[.], G->L[.].
  Attention ! G->n n'est pas modifiÃ©.
*/
  int const du=G->d[u];
  int i,j,v,dv;
  G->d[u]=0;
  for(i=0;i<du;i++){
    v=G->L[u][i];
    dv=G->d[v];
    for(j=0;j<dv;j++)
      if(G->L[v][j]==u)	G->L[v][j]=G->L[v][--dv];
    G->d[v]=dv;
  }
  G->sort=0;
  return;
}


int AdjGraph(const graph* G,int const u,int const v){
/*
  Renvoie 1 ssi dans le graphe G le sommet v apparaÃ®t dans la liste
  d'adjacence de u. Si les listes du graphe G sont est triÃ©es (G->sort
  est vrai), une dichotomie est effectuÃ©e (donc en log(deg(u))) sinon
  c'est un parcours sÃ©quentiel (en deg(u)).
*/

  /* recherche dichotomique si le graphe est triÃ© */
  if(G->sort)
    return (bsearch(&v,G->L[u],G->d[u],sizeof(int),fcmp_int)!=NULL);

  /* sinon, recherche sÃ©quentielle */
  int i;
  for(i=G->d[u];i>0;)
    if(G->L[u][--i]==v) return 1; /* on a trouvÃ© v dans la liste de u */
  return 0;
}


int Treewidth(graph* const H,int const code)
/*
  Calcul approchÃ© ou excate de la treewidth de H. Si code=0, on
  calcule une borne supÃ©rieure (heuristique degmin). Si code=1 on
  calcule la valeur exacte. H n'est pas modifiÃ©. Dans H->int1 on
  renvoie le nb de tests effectuÃ©s.
*/
{
  if(H==NULL) return 0;
  int n,j,t,tw,l,v,n1;
  int d,r,i,u,k;
  int d0,r0,i0,u0;
  graph* G;

  H->int1=0;
  n=H->n;
  n1=n-1;
  tw=(code==1)? Treewidth(H,0) : min(n1,nb_edges(H));

  /* tw=upper bound sur tw. On a intÃ©rÃªt d'avoir la valeur la plus
     faible possible, car on va alors Ã©liminer les mauvais ordres
     d'Ã©liminations plus rapidement.

    tw max en fonction du nb m d'arÃªtes:
    m=0 => tw=0
    m=1 => tw=1
    m=2 => tw=1
    m=3 => tw<=2
    m=4 => tw<=2
    m=5 => tw<=2?
    m=6 => tw<=3
    ...
   */

  NALLOCZ(int,P,n,_i); /* permutation initiales des sommets */
  NALLOCZ(int,D,n,n1);

  do{
    G=GraphCopy(H); /* on copie H */
    GraphRealloc(G,D); /* plonge G dans le graphe complet */
    k=0; /* tw par dÃ©faut si G sans d'arÃªtes */

    for(j=0;j<n1;j++){ /* n-1 fois */

      H->int1++;
      u0=P[i0=j];
      d0=r0=G->d[u0];
      if(code==0){
	for(i=j;i<n1;i++){
	  u=P[i]; d=G->d[u];
	  if(d<=d0){
	    /* calcule r=nb d'arÃªtes dans N(u) */
	    for(r=l=0;l<d;l++)
	      for(t=l+1;t<d;t++)
		r += AdjGraph(G,G->L[u][l],G->L[u][t]);
	    if((d<d0)||(r<r0)){ d0=d;r0=r;i0=i;u0=u; }
	  }
	} /* fin for(i=...) */
	P[i0]=P[j]; /* dÃ©cale P[i], que si code=0 */
      }
      k=max(k,d0); /* met Ã  jour tw */
      if(k>=tw) goto nextP; /* on ne fera pas moins */
      RemoveVertex(G,u0); /* supprime u */
      /* remplace N(u) par une clique */
      for(i=0;i<d0;i++){
	u=G->L[u0][i];
	for(t=i+1;t<d0;t++){
	  v=G->L[u0][t];
	  if(!AdjGraph(G,u,v)) ADD_EDGE(G,u,v);
	}
      }
      
    } /* fin for(j=...) */

    tw=min(tw,k);
  nextP:
    free_graph(G);
    if((code==0)||(tw<2)) break; /* si tw(G)=0 ou 1, alors on ne
				    trouvera pas plus petit */
  }while(NextPermutation(P,n,NULL));

  free(D);
  free(P);
  return tw;
}


int *Prune(const graph* G,int *z)
/*
  Algorithme permettant de supprimer les sommets du graphe G par
  dÃ©grÃ©s minimum. Un ordre d'Ã©lÃ©mination des sommets est renvoyÃ© sous
  la forme d'un tableau de taille n. Plus prÃ©cisÃ©ment, on renvoie un
  tableau de sommets T de sorte que u=T[i] est un sommet de degrÃ©
  minimum du graphe G\(T[0]...T[i-1]). La complexitÃ© est linÃ©aire,
  i.e. O(n+m). Si *z<>NULL, on retourne dans z le degrÃ© maximum
  rencontrÃ© donc la dÃ©gÃ©nÃ©rescence du graphe.
*/
{
  int i,j,u,v,d,k,p,r,t,c,n1;
  int const n=G->n;

  /*
    1. On construit les tableaux suivants:

    T[0..n[ = tableau de sommets triÃ©s par degrÃ© croissant (le rÃ©sultat)
    R[0..n[ = tableau de positions dans T des sommets
    D[0..n[ = tableau des degrÃ©s des sommets
    P[0..n[ = tableau des positions dans T du premier sommet de degrÃ© donnÃ©

    2. On traite les sommets dans l'ordre T[0],T[1],... Supposons que
    u=T[i]. Alors, pour chaque voisin v dans G, encore existant (donc
    situÃ© aprÃ¨s i dans le tableau T), on met Ã  jour la structure
    T,R,D,P.
   */

  /* initialise T,R,D,P */
  if(n<1) return NULL;
  NALLOCZ(int,D,n,G->d[_i]);
  int *R=SortInt(D,NULL,n,0,NULL,SORT_INDEXe);
  NALLOC(int,T,n); for(u=0;u<n;u++) T[R[u]]=u;
  NALLOCZ(int,P,n,-1);
  for(i=n1=n-1;i>=0;i--) P[D[T[i]]]=i; /* i=n-1,...,0 */

  for(c=i=0;i<n1;i++){ /* pour chaque sommet, pas besoin de faire le dernier */
    u=T[i]; /* u = sommet de degrÃ© (D[u]) minimum */
    if(D[u]>c) c=D[u]; /* mÃ©morise le plus grand degrÃ© rencontrÃ© */
    d=G->d[u]; /* d = deg(u) */
    for(j=0;j<d;j++){ /* pour chaque voisin */
      v=G->L[u][j]; /* v=voisin de u */
      r=R[v]; /* v=T[r]; */
      if(r>i){ /* mettre Ã  jour que si v existe encore */
	k=D[v]--; /* le degrÃ© de v diminue, k=D[v] juste avant */
	p=P[k--]++; /* on va Ã©changer v et T[p]. Les sommets de degrÃ© k commencent juste aprÃ¨s */
	if(P[k]<0) P[k]=p; /* y'a un sommet de degrÃ© k en plus */
	if(p>i){ /* si p est avant i, ne rien faire. Cel arrive que si k=0 */
	  t=T[p]; /* t = 1er sommet de de degrÃ© k */
	  T[p]=v; /* dans T, on avance v en position p (Ã  la place de t donc) */
	  T[r]=t; /* on met t Ã  la place de v */
	  R[v]=p; /* on met Ã  jour la nouvelle position de v dans T */
	  R[t]=r; /* on met Ã  jour la nouvelle position de t dans T */
	}
      }
    }
  }

  if(z!=NULL) *z=c; /* retourne la dÃ©gÃ©nÃ©rescence */
  free(R);
  free(D);
  free(P);
  return T;
}


int *GreedyColor(graph* const G,int* const R)
/*
  Algorithme permettant de colorier de maniÃ¨re gloutonne un graphe G Ã 
  n sommets. La complexitÃ© en temps est linÃ©aire, i.e. O(n+m). Il faut
  G non NULL. On renvoie un tableau C[0..n[ oÃ¹ C[u] est la couleur du
  sommet u, un entier entre 0 et n-1. Le tableau R[0..n[ donne un
  ordre (inverse) de visite des sommets: On colorie le sommet u avec
  la plus petite couleur possible dans le graphe induit par les
  sommets situÃ© aprÃ¨s u dans R (on commence donc avec R[n-1]). Si
  R=NULL, alors l'ordre R[i]=i est utilisÃ©. On renvoie dans G->int1 la
  couleur maximum utilisÃ©e (donc le nb de couleurs-1). Cette valeur
  vaut -1 si n<1.  */
{
  int i,j,u,d,l,c;
  int const n=G->n;

  /*
    On utilise 3 tableaux:

    C[0..n[ = tableau final des couleurs, au dÃ©part =-1
    L[0..n-1[ = liste de couleurs utilisÃ©es par le sommet courant
    M[0..n-1[, M[i]=1 ssi la couleur i est utilisÃ©e par un voisin du
    sommet courant, au dÃ©part =0. On met une sentiennelle Ã  M si bien
    que toujours on M[n-1]=0.

  */

  G->int1=-1;
  if(n<1) return NULL;
  NALLOCZ(int,C,n,-1);
  NALLOCZ(int,M,n,0);
  i=n-1;
  NALLOC(int,L,i);

  for(;i>=0;i--){ /* en partant de la fin */
    u=(R==NULL)? i : R[i];
    d=G->d[u]; /* d = deg(u) */

    /* on liste dans L et M les couleurs rencontrÃ©es */
    for(j=l=0;j<d;j++){ /* pour chaque voisin v de u */
      c=C[G->L[u][j]]; /* c=couleur du voisin v */
      if(c>=0){ /* voisin dÃ©jÃ  coloriÃ© ? */
	L[l++]=c; /* on ajoute la couleur c Ã  la liste */
	M[c]=1; /* la couleur c n'est pas Ã  prendre */
      }
    }
    
    /* on cherche la 1Ã¨re couleur Ã  1 (=non-rencontrÃ©e) */
    j=0; while(M[j]) j++; /* s'arrÃªte toujours Ã  cause de la sentiennelle */
    C[u]=j; /* j est la couleur recherchÃ©e pour u */
    G->int1=max(G->int1,j); /* couleur maximum rencontrÃ©e */

    /* il faut rÃ©-initialiser rapidement M */
    for(j=0;j<l;j++) M[L[j]]=0; /* on efface les 1 qu'on Ã  mis dans M */
  }
  
  free(L);
  free(M);
  return C;
}


void HalfGraph(graph* const G,int const code){
/*
  Transforme G en un graphe asymÃ©trique (avec G->sym=1) en graphe
  symÃ©trique en ne gardant que les arcs u->v tels que v<u si code=0 ou
  bien tels que u<v si code=1. En particulier les boucles sont
  supprimÃ©s. Les listes d'adjacence sont aussi triÃ©es par ordre
  croissant et les tables raccourcies (realloc). Cela permet de faire
  des parcours de graphe deux fois plus rapidement, par exemple, pour
  vÃ©rifier qu'une coloration est propre. Le champs ->sym est mis Ã 
  jour. L'exÃ©cution est plus rapide si code=0.
*/
  if(G==NULL) return;
  if(!G->sym) return;

  int const n=G->n;
  int i,u,d;
  NALLOC(int,D,n); // tableau des nouveaux degrÃ©s

  SortGraph(G,0); // tri par ordre croissant les listes
  
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      if(G->L[u][i]>=u) break;
    if(code){ // il faut v>u
      while((i<d)&&(G->L[u][i]==u)) i++; // les boucles
      memmove(G->L[u],G->L[u]+i,(G->d[u]-i)*sizeof(int));
      i=G->d[u]-i;
    }
    D[u]=i; // coupe la liste Ã  i
  }
  
  GraphRealloc(G,D); // modifie toutes les listes d'adjacence
  free(D);
  G->sym=0; // G est dÃ©sormais asymÃ©trique
  return;
}


void kColorSat(graph* const G,int const k){
/*
  Affiche les contraintes SAT pour la k-coloration de G au format
  Dimacs. Attention ! G est modifiÃ© (manque la moitiÃ© de ses arcs).

  Le nombre de variables est n*k.
  Le nombre de clause est n+m*k.
*/

  if(G==NULL) return;
  int const n=G->n;
  int const m=nb_edges(G);
  int c,u,i,d;
  
  printf("p cnf %i %i\n",n*k,n+m*k);

  /*
    Variables: chaque sommet i a k variables numÃ©rotÃ©s x(i,c) avec
    i=0..n-1 et c=0..k-1. Il faut x(i,c)>0. On pose x(i,c)=1+k*i+c.

    Contraintes sommets: il faut x(i,0) v ... v x(i,k-1), et la
    conjonction pour tous les sommets.

    Containtes arÃªtes: pour chaque arÃªte {i,j} et couleur c il ne faut
    pas que x(i,c)=x(j,c). Autrement dit il faut -x(i,c) v -x(j,c). Il
    faut la conjonction de toutes les couleurs c et toutes les arÃªtes
    {i,j}.
  */

  /* liste chaque sommet */
  for(u=0;u<n;u++){
    for(c=0;c<k;c++) printf("%i ",1+u*k+c);
    printf("0\n");
  }

  HalfGraph(G,0); /* enlÃ¨ve la moitiÃ© des arcs */

  /* liste chaque arc */
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      for(c=0;c<k;c++)
	printf("-%i -%i 0\n",1+u*k+c,1+G->L[u][i]*k+c);
  }
  
  return;
}


void kIndepSat(graph* const G,int const k){
/*
  Affiche les contraintes SAT au format Dimacs pour ensemble
  indÃ©pendant de taille k pour le graphe G. Attention ! G est modifiÃ©
  (manque la moitiÃ© de ses arcs).

  Variables:

    Pour chaque i=0..n-1:
    X(i)=1 ssi le sommet i est dans l'ensemble indÃ©pendant.
    Si i-j est une arÃªte alors -X(i) v -X(j).
    NB: pour Vertex Cover, il faudrait suffit d'ajouter X(i) v X(j).

    Pour chaque t=0..n et b=0..k:
    Y(t,b)=1 ssi la somme des t variables X_0+...+X(t-1) est au
    moins b. On a:

       Y(n,k) = 1
       Y(t,0)=1 pour tout t>=0
       Y(t,b)=0 pour tout 0<=t<b et b>0
       Y(t,b) => Y(t-1,b) v (Y(t-1,b-1) ^ X(t-1))

       <=> -Y(t,b) v Y(t-1,b) v Y(t-1,b-1) ET
           -Y(t,b) v Y(t-1,b) v X(t-1)

    #variables X(i): n
    #variables Y(t,b): (n+1)*(k+1)
    #variables totales: n+(n+1)*(k+1)

*/

  if(G==NULL) return;

  int const n=G->n;
  int const m=nb_edges(G);
  int i,j,t,d,b;  

  // il faut entrer le nombre exactes de variables et de clauses
  printf("p cnf %i %i\n",n+(n+1)*(k+1),m+n+2+k*(k+1)/2+k*(2*n+1-k));

  /* numÃ©ros des variables: attention de ne surtout pas utiliser la
     variable numÃ©ro 0 qui signifie "fin de ligne" */

#define X(i)   ((i)+1)             // numÃ©ro de la variable X(i)
#define Y(t,b) (n+1+(t)*(k+1)+(b)) // numÃ©ro de la variable Y(t,b)

  HalfGraph(G,0); /* rend le graphe simple et asymÃ©trique */

  // pour chaque arÃªtes i-j: -X(i) v -X(j)
  // #clauses: m
  for(i=0;i<n;i++){
    d=G->d[i];
    for(j=0;j<d;j++)
      printf("-%i -%i 0\n",X(i),X(G->L[i][j]));
  }

  // Y(n,k)=1
  // #clause: 1
  printf("%i 0\n",Y(n,k));

  // cas b=0 et t>=0: Y(t,0)=1
  // #clauses: n+1
  for(t=0;t<=n;t++) printf("%i 0\n",Y(t,0));

  // cas b>=1 et 0<=t<b: Y(t,b)=0
  // #clauses: 1+2+3+...+k = k*(k+1)/2
  for(b=1;b<=k;b++)
    for(t=0;t<b;t++)
      printf("-%i 0\n",Y(t,b));

  // cas b>=1 et t>=b: rÃ©currence
  // #clauses: 2*âˆ‘_{b=1}^k (n-b+1) = 2*(n+1)*k - k(k+1) = k*(2*n+1-k)
  for(b=1;b<=k;b++)
    for(t=b;t<=n;t++){
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),Y(t-1,b-1));
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),X(t-1));
    }

  return;
}
#undef X
#undef Y


int* kColor(graph* const G,int const k){
/*
  Algorithme permettant de colorier en au plus k couleurs un graphe G,
  si c'est possible. La complexitÃ© est (n+m)*k^{n-1} dans le pire des
  cas.  On renvoie un tableau C[0..n[ oÃ¹ C[u] est la couleur du sommet
  u, un entier entre 0 et k-1. On renvoie NULL s'il n'est pas possible
  de colorier G en k couleurs, si G=NULL ou k<1. On renvoie dans
  G->int1 la couleur maximum utilisÃ©e (qui peut Ãªtre < k-1). On
  utilise toutes les couleurs de [0,G->int1]. On symÃ©trise le graphe,
  qui est donc modifiÃ©, afin d'enlever la moitiÃ© des arÃªtes Ã 
  vÃ©rifier.

  La stratÃ©gie est la suivante. On part d'une coloration initiale des
  sommets C=[0,...,0,k-1] oÃ¹ la couleur du dernier sommet u=n-1 est
  fixÃ©e par C[n-1]=k-1. Puis on vÃ©rifie si C est propre ou non. Si
  c'est non, on incrÃ©mente C comme un compteur, et on recommence avec
  la coloration suivante. Ainsi toutes les colorations possibles de G
  sont passÃ©es en revue. On teste toujours en prioritÃ© la derniÃ¨re
  arÃªte qui a posÃ© problÃ¨me avant de vÃ©rifier tout le graphe.

  Optimisations possibles Ã  faire:
  
  1. RÃ©duction de donnÃ©es. On peut supprimer rÃ©cursivement les sommets
     de degrÃ© < k. On pourra toujours les ajouter Ã  la fin si la
     coloration a rÃ©ussie.

  2. On peut dÃ©composer le graphe en composantes connexes, il y a
     ainsi moins de colorations possibles Ã  tester.

  3. On peut renumÃ©roter les sommets selon un parcours BFS. De cette
     faÃ§on la vÃ©rification pourrait Ãªtre accÃ©lÃ©rÃ©e lorsqu'on change
     une seule couleur.

*/

  int c,i,d,u,v,b,*T;
  if((G==NULL)||(k<1)) return NULL;
  HalfGraph(G,0); /* enlÃ¨ve la moitiÃ© des arÃªtes */
  int const n=G->n;

  NALLOCZ(int,C,n,0); /* C[u]=couleur du sommet u */
  C[n-1]=k-1; /* on peut supposer que le sommet n-1 a une couleur fixÃ©e */
  if(n<2) goto fin_kcolor; /* s'il y a un seul sommet */
  b=1; /* b=vrai ssi la derniÃ¨re arÃªte coloriÃ©e est propre */

  do{
    /* vÃ©rifie si C est une coloration propre */
    if(b){ /* on le fait que si b est vrai */
      for(u=0;u<n;u++){ /* pour chaque sommet */
	c=C[u]; d=G->d[u]; /* degrÃ© et couleur de u */
	for(i=0;i<d;i++) /* pour chaque voisin de i */
	  if(c==C[G->L[u][i]]) break; /* coloration pas propre */
	if(i<d) break; /* coloration pas propre */
      }
      if(u==n) goto fin_kcolor; /* la coloration est propre */
    }
    /* ici l'arÃªte (u,i) n'est pas correctement coloriÃ©e */
    
    /* on change C en incrÃ©mentant le compteur */
    v=0;
  loop_kcolor:
    C[v]++;
    if(C[v]==k){
      C[v++]=0;
      if(v==n) break;
      goto loop_kcolor;
    }
    
    /* est-ce que l'arÃªte (u,i) est toujours mal coloriÃ©e ?
       si oui, pas la peine de vÃ©rifier tout le graphe */
    b=(C[u]!=C[G->L[u][i]]); /* b=vrai ssi (u,i) est propre */

  }while(v<n);

  /* aucune coloration trouvÃ©e */
  free(C);
  return NULL;

 fin_kcolor:
  /* on a trouvÃ© une coloration propre */
  /* on rÃ©duit l'espace des couleurs utilisÃ©es.  NB: ici on a encore
     C[n-1]=k-1 */

  ALLOCZ(T,k,-1);
  /* si T[c]=i, alors la couleur c se renumÃ©rote en i. Si c n'a jamais
     Ã©tÃ© rencontrÃ©e, alors i=-1 */

  for(i=u=0;u<n;u++){
    if(T[C[u]]<0) T[C[u]]=i++; /* la couleur C[u] n'a pas jamais Ã©tÃ© vue */ 
    C[u]=T[C[u]]; /* ici la couleur C[u] a dÃ©jÃ  Ã©tÃ© vue */
  }

  free(T); /* plus besoin de T */
  G->int1=i-1; /* couleur max utilisÃ©e */
  return C;
}


int* power_law_seq(int const n,double const t,int *T){
/*
  Ã‰crit dans le tableau T une distribution de degrÃ© pour un graphe Ã 
  n>0 sommets selon une lois en puissance d'exposant t>0. On renvoie
  NULL si les paramÃ¨tres n et t ne sont pas corrects. La distribution
  renvoyÃ©e dans T est codÃ©e par une suite (2k,n_1,d_1,...,n_k,d_k)
  contenant k paires (n_i,d_i) signifiant n_i sommets de degrÃ© d_i. La
  taille de T est donc T[0]+1 = 2k+1. Si T=NULL, alors T est allouÃ© et
  renvoyÃ©, sinon il doit Ãªtre assez grand pour recevoir la
  distribution. La taille de T, sans jamais dÃ©passer 2n+1, vaut en
  thÃ©orie |T| = O(n^{1/t}). En pratique on a:

    t    n     |T|        t    n     |T|        t    n     |T|
   ----------------      ----------------      ----------------
    2.0  10^3  50         2.5  10^3  28         3.0  10^3  18
    2.0  10^4  156        2.5  10^4  70         3.0  10^4  40
    2.0  10^5  494        2.5  10^5  176        3.0  10^5  86
    2.0  10^6  1560       2.5  10^6  446        3.0  10^6  188
    2.0  10^7  4932       2.5  10^7  1122       3.0  10^7  404
   ----------------      ----------------      ----------------

  La distribution est la suivante:
  - d_1=1, n_1 = âŒŠexp(a)âŒ‹ + n-s(a)
  - d_i=i, n_i = âŒŠexp(a)/i^tâŒ‹ pour i dans [2,p(a)]
  - a est un rÃ©el minimisant |n-s(a)| avec
    s(a) := âˆ‘_{i=1}^{p(a)} âŒŠexp(a)/i^tâŒ‹
    p(a) := âŒŠexp(a/t)âŒ‹

  Attention ! Dans l'article original de [Lu01], il y a une erreur
  dans le signe du terme correctif r=n-s(a) pour les sommets de degrÃ©
  1. Il faut faire +r et non -r comme c'est Ã©crit.

  Le nombre de paires (n_i,d_i) dans T est exactement p(a). Si n>0,
  alors T contient au moins 1 paire (et au plus n), car les d_i sont
  diffÃ©rents. On a s(0)=1 car p(0)=1. Aussi si a>=ln(n)+1, alors s(a)
  â‰¥ floor{exp(a)} â‰¥ n.

  En choisissant a0=0 et a1=ln(n)+1, on a alors s(a0) <= n <=
  s(a1). La fonction s(a) Ã©tant croissante, pour minimiser n-s(a) on
  rÃ©alise une recherche binaire pour a dans [a0,a1]. Le nombre
  d'itÃ©rations est au plus 64 (car sizeof(double)*8=64).
*/
  if((t<=0)||(n<=0)) return NULL;

  /* calcule la valeur de 'a' optimale */
  
  double a0=0,a1=log(n)+1.0; // intervalle pour a
  double a,b,e; // a=milieu de [a0,a1], b=meilleur a
  int p,i,j,s,cont=1; // cont=1 ssi on continue le calcul
  int r=0; // r=valeur minimum (=0 au dÃ©part)
  
  do{
    a=(a0<a1)? (a0+a1)/2 : a0; // si a0>=a1, on calcule s puis on s'arrÃªte
    if((a==a0)||(a==a1)) cont=0; // intervalle trop faible, on calcule s puis on s'arrÃªte
    e=exp(a); p=(int)exp(a/t); // p=nombre de paires, NB: p>=1
    for(s=0,i=1;i<=p;i++) s += (int)(e/pow(i,t)); // calcule s=s(a)
    if(s==n) cont=0; // valeur optimale, on va avoir b=a
    if((r==0)||(abs(n-s)<abs(r))) b=a,r=n-s; // NB: si s=n, alors b=a
    if(s<n) a0=a; // on est avant n
    if(s>n) a1=a; // on est aprÃ¨s n
  }while(cont);

  /* ici on a calculÃ© b, le meilleur a */

  e=exp(b);
  p=(int)exp(b/t); // nombre de paires
  if(T==NULL) ALLOC(T,2*p+1);
  T[0]=2*p; // taille du tableau

  /* Ã©crit la distribution dans T */
  
  for(j=i=1;i<=p;i++){
    T[j++]=(int)(e/pow(i,t)); // NB: si i=1, T[1]=floor(exp(b))
    T[j++]=i;
  }
  DEBUG(PRINT(j););

  T[1]+=r; // correction pour les sommets de degrÃ© 1
  
  return T;
}


/***********************************

         ROUTINES POUR LES
       FONCTIONS D'ADJACENCE

***********************************/

// Pour avoir la liste de toutes les fonctions d'adjacence:
// grep '^int.*(query\* const Q)' gengraph.c

/* code souvent utilisÃ© dans les fonctions d'adjacence */

#define TEST_n       do{ if(Q->n<1){ Q->n=0; return 0; }}while(0)
#define SET_n(X)     do{ Q->n=(X); TEST_n; }while(0)
#define RET_n(X)     do{ Q->n=(X); if(Q->n<1) Q->n=0; return 0; }while(0)
#define RET_a(X)     do{ Q->a=(X); return 0; }while(0)
#define RET_error(X) do{ Q->error=(X); return 1; }while(0)

/*
#include <boost/preprocessor/control/if.hpp>
#define RET_n(X)     BOOST_PP_IF(0,					\
				 do{ return Q->n=0; }while(0),		\
				 do{ Q->n=(X); if(Q->n<1) Q->n=0; return 0; }while(0) \
				 )
*/

int free_pos(query* const Q)
/*
  Routine permettant de libÃ©rer les tableaux Q->xpos, XSEED, ...
  
  La fonction devrait toujours Ãªtre exÃ©cutÃ©e dans une fonction
  d'adjacence avec Q->code=QUERY_END. Renvoie 0 si tout c'est bien
  passÃ©, 1 sinon.
*/
{
  free(Q->xpos),  free(Q->ypos);
  free(XSEED), free(YSEED);
  Q->xpos=Q->ypos=XSEED=YSEED=NULL;
  return 0;
}


int free_rep(query* const Q)
/*
  Routine permettant de libÃ©rer le tableau Q->rep. Renvoie 0 si tout
  c'est bien passÃ©, 1 sinon.
*/
{
  free(Q->rep);
  Q->rep=NULL;
  Q->k=0;
  return 0;
}


int adjacency_rep(query* const Q)
/*
  Fonction d'adjacence entre Q->i et Q->j commune Ã  beaucoup de
  graphes. Elle est basÃ©e sur une k-orientation, chaque sommet u ayant
  une liste d'au plus k voisins stockÃ©s dans Q->rep[u][0..k[ oÃ¹
  k=|Q->k|. Cette k-orientation peut comprendre des boucles et des
  arcs symÃ©triques.

  Le rÃ©sultat, renvoyÃ© dans Q->a, dÃ©pend de Q->directed: si le graphe
  est symÃ©trique il faut faire un double test. Modifie Q->a et renvoie
  toujours 0. On arrÃªte le parcours de Q->rep[Q->i][t] que lorsque
  t=Q->k-1. Si Q->k<0, on arrÃªte le parcours Ã  la premiÃ¨re valeur <0
  rencontrÃ©e ou si t=|Q->k|-1.
*/
{
  int const b=(Q->k<0);
  int const k=abs(Q->k);
  int t;

  // cherche j dans rep[i]
  for(t=0;t<k;t++){
    if(Q->rep[Q->i][t]<0){ if(b) break; else continue; }
    if(Q->rep[Q->i][t]==Q->j) RET_a(1);
  }
  if(Q->directed) RET_a(0); // on s'arrÃªte lÃ  si orientÃ©

  // cherche i dans rep[j] si non orientÃ©
  for(t=0;t<k;t++){
    if(Q->rep[Q->j][t]<0){ if(b) break; else continue; }
    if(Q->rep[Q->j][t]==Q->i) RET_a(1);
  }
  RET_a(0);
}


/***********************************

         GRAPHES DE BASE

***********************************/


/*
  Une fonction d'adjacence de graphe, disons adj(Q), effectue un
  calcul en fonction de la requÃªte Q->code et des divers paramÃ¨tres
  contenus dans Q. Le code typique ressemble donc Ã  un switch(Q->code)
  comme ci-dessous.

  int adj(query* const Q)
  {
    switch(Q->code){
    case QUERY_INIT: Q->n=...; return 0; // ou RET_n(...);
    case QUERY_ADJ:  Q->a=...; return 0; // ou RET_a(...);
    }
    return 1;
  }

  Toutes les fonctionalitÃ©s ne sont pas forcÃ©ment implÃ©mentÃ©es.
  Cependant, a minima, QUERY_INIT (pour dÃ©terminer Q->n) et QUERY_ADJ
  (pour dÃ©terminer Q->a) doivent Ãªtre implÃ©mentÃ©es. La fonction doit
  renvoyer 0 si tout c'est bien passÃ©e, et 1 sinon. Dans ce cas, il
  est possible de renseigner Q->error. On s'aperÃ§oit donc si une
  fonctionalitÃ© est codÃ© en testant si le code retour est 0 ou 1.

  Rem1: PlutÃ´t que d'utiliser "return 0" Ã  la fin de chaque "case", on
  pourrait mettre un "break" puis un "return 0" en dehors du "switch"
  avec aussi un "default: return 1". Cependant un "return" direct est
  plus rapide, ce qui peut Ãªtre important pour QUERY_ADJ et QUERY_NAME
  par exemple.

  Rem2: Attention aux subtilitÃ©s de l'instruction "switch". Les
  diffÃ©rents "case" correspondent Ã  des Ã©tiquettes comme un "goto". En
  particulier, les dÃ©clarations avec initialisation ne sont pas
  valides juste aprÃ¨s un "case". Une solution est de dÃ©placer (plus
  loin) la dÃ©claration ou d'inserer une instruction (Ã©ventuellement
  vide ";") entre le "case" et la dÃ©claration, ou encore de crÃ©er un
  block. Voici un rÃ©sumÃ©:

    (incorrect)        (correct)          (correct)           (correct)
  
  case QUERY_ADJ:    case QUERY_ADJ:    case QUERY_ADJ:;    case QUERY_ADJ:
    int k=42;          int k;             int k=42;           { int k=42;
    ...;               k=42;              ...;                  ...;
                       ...;                                     }
   
  Rem3: Il faut Ã©viter de traiter l'erreur directement dans QUERY_INIT
  avec Erreur(...) car les fonctions peuvent s'appeler en cascade, et
  on peut ne plus comprendre d'oÃ¹ vient l'erreur. C'est Ã  l'appelant
  de gÃ©rer l'erreur (en rÃ©cupÃ©rant le code d'erreur et en prenant la
  dÃ©cision qui s'impose. Dans la plupart des cas on peut se contenter
  de renvoyer le graphe vide avec RET_n(0) si les paramÃ¨tres sont
  incorrects. Si l'on souhaite renvoyer une erreur particuliÃ¨re on
  peut (devrait) faire:

    case QUERY_INIT:
      ...;
      if(...) RET_error(42);
      ...;

  Les diffÃ©rentes fonctionnalitÃ©s sont (pour Q->code):

    QUERY_INIT: dÃ©termine Q->n, initialise la fonction
    QUERY_END: termine la fonction
    QUERY_ADJ: dÃ©termine Q->a en fonction de Q->i, Q->j et des paramÃ¨tres
    QUERY_NAME: dÃ©termine Q->name, le nom du sommet Q->i

    [et Ã  finaliser ...]

    QUERY_LIST: dÃ©termine Q->L, la liste des voisins de Q->i
    QUERY_DOT: dÃ©termine le dessin d'une arÃªte au format dot
    QUERY_DEG: dÃ©termine le degrÃ© du sommet Q->i
    QUERY_LNAME: pour la liste des noms des sommets

  Certaines fonctions d'adjacence supposent Q->i < Q->j. Cela signifie
  que le rÃ©sultat peut ne pas Ãªtre correct, voir produire une erreur
  grave, si l'option -directed est prÃ©sente.

  C'est une mauvaise idÃ©e que d'utiliser des variables statiques dans
  les fonctions d'ajacence car elles peuvent s'appeler entre elles
  pour factoriser des codes qui sont proches. Si les prÃ©-calculs sont
  les bienvenus dans QUERY_INIT pour optimiser QUERY_ADJ, il faut
  impÃ©rativement stocker leurs rÃ©sultats dans les champs de la
  variable Q.

  On devrait systÃ©matiquement, lors de l'initialisation avec
  QUERY_INIT, avoir le test "if(Q->n<1){ Q->n=0; return 0; }" ou
  encore la macro "TEST_n;" (le "=0" dans le "if(Q->n<1)" est
  important) permettant de gÃ©nÃ©rer le graphe vide sans rien faire
  d'autre surtout pour les graphes allouant de la mÃ©moire (Q->rep,
  Q->xpos, ...) et qui nÃ©cessitent que Q->n>0.

  Les graphes utilisant le tableau Q->rep[u] (reprÃ©sentation
  implicite) pour chaque sommet u devraient se terminer par la
  libÃ©ration des tableaux comme ceci:

    case QUERY_END:
      return free_rep(Q);

  De mÃªme, les graphes gÃ©omÃ©triques qui utilisent les tableaux Q->xpos et
  Q->ypos devraient se terminer par leur libÃ©ration comme ceci:
  
    case QUERY_END:
      return free_pos(Q);

  Et si la fonction utilise les types de tableaux:

    case QUERY_END:
      return free_rep(Q)||free_pos(Q);

  Pour les graphes gÃ©omÃ©triques, ce n'est pas la peine de dÃ©finir
  QUERY_NAME. Les coordonnÃ©es des sommets pouvent Ãªtre affichÃ©es par
  l'option -label -3.
*/


int load(query* const Q)
/*
  Graphe dÃ©fini par un fichier (ou l'entrÃ©e standard). Ã€
  l'initialisation, il est chargÃ© en mÃ©moire dans la variable Q->G de
  type "graph". Suivant la valeur de Q->G->sort, le test d'adjacence
  est linÃ©aire ou logarithmique en min{deg(i),deg(j)}.
*/
{
  switch(Q->code){
    
  case QUERY_END:
    free_graph(Q->G);
    Q->G=NULL;
    return 0;

  case QUERY_INIT:
    Q->G=File2Graph(Q->sparam,2); /* remplit Q->G Ã  partir du fichier */
    if(Q->G->f>0){ /* si c'est une famille, on sÃ©lectionne le premier graphe */
      graph* G=GraphCopy(Q->G->G[0]); /* copie le premier graphe */
      free_graph(Q->G); /* libÃ¨re complÃ¨tement la famille Q->G */
      Q->G=G; /* Q->G=premier graphe */
    }
    if(!Q->G->sym) Q->directed=1; /* si le graphe est asymÃ©trique */
    RET_n(Q->G->n);

  case QUERY_ADJ:
    /* pour avoir du min{deg(i),deg(j)} en non-orientÃ© */
    if((!Q->directed)&&(Q->G->d[Q->i]>Q->G->d[Q->j])) Q->a=AdjGraph(Q->G,Q->j,Q->i);
    else Q->a=AdjGraph(Q->G,Q->i,Q->j);
    return 0;

  }

  return 1;
}


int prime(query* const Q)
{
  switch(Q->code){
  case QUERY_INIT: RET_n(Q->param[0]);
  case QUERY_ADJ: RET_a( (Q->i>0)&&(Q->j>1)&&((Q->i%Q->j)==0) );
  }
  return 1;
}


int paley(query* const Q)
/*
  Le rÃ©sidu est r=|i-j|. Pour savoir si r est un carrÃ©, on teste s'il
  existe un entier k<=(n-1)/2 tel que (k*k)%n=r. Le nombre de carrÃ©s
  si dans Z/nZ, si n est premier vaut (n+1)/2.
*/
{
  const int n=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(n);
  case QUERY_ADJ:;
    const int q=n/2;
    const int r=abs(Q->i-Q->j);
    for(int k=1;k<=q;k++) if(r==((k*k)%n)) RET_a(1);
    RET_a(0);
  }
  return 1;
}


int comb(query* const Q)
/*
  Utilise i<j.
  On pourrait utiliser "grid 1 n -star -1".
*/
{
  const int n=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(2*n);
  case QUERY_ADJ: RET_a( ((Q->j<n)&&(Q->j==Q->i+1))||(Q->j==Q->i+n) );
  }
  return 1;
}


int sunlet(query* const Q)
/*
  Utilise i<j.
  On pourrait utiliser: "grid 1 -n -apex -1" 
*/
{
  if((Q->code==QUERY_ADJ)&&(Q->i==0)&&(Q->j==Q->param[0]-1)) RET_a(1);
  return comb(Q);
}


/* types pour alkane() */
enum{
  ALK_NOR,
  ALK_CYC,
  ALK_ISO,
  ALK_NEO,
  ALK_SEC,
  ALK_TER,
};


int alkane(query* const Q)
/*
  Utilise i<j.

  Les sommets C ont des numÃ©ros < n.
  Les sommets H ont des numÃ©ros â‰¥ n.
  Pour les sommets C, dans l'ordre on trouve:
  - les sommets de degrÃ© 1
  - les sommets de degrÃ© 3
  - les sommets de degrÃ© 4
  - les sommets de degrÃ© 2
  Cette numÃ©rotation permet de gÃ©rer de maniÃ¨re uniforme et simple
  les adjacences C-H. Les adjacences C-C se font aux cas par cas.

  Q->wrap[0]: borne infÃ©rieure pour n
  Q->wrap[1]: nombre de sommets C de degrÃ© 1
  Q->wrap[2]: nombre de sommets C de degrÃ© 3
  Q->wrap[3]: nombre de sommets C de degrÃ© 4

  VidÃ©os sur le nom des alkanes:
  https://fr.khanacademy.org/science/organic-chemistry/
  bond-line-structures-alkanes-cycloalkanes#naming-alkanes
*/
{
  const int type=Q->param[0];
  const int n=Q->param[1]; // n=nombre d'atomes de C
  static int P_NOR[]={1,2,0,0};
  static int P_CYC[]={3,0,0,0};
  static int P_ISO[]={4,3,1,0};
  static int P_NEO[]={5,4,0,1};
  static int P_SEC[]={6,3,1,0};
  static int P_TER[]={7,4,0,1};
  int a,b=((n==1)&&(type==ALK_NOR)); // b=1 ssi methane: un seul C de degrÃ© 0

  switch(Q->code){

  case QUERY_END:
    Q->wrap=NULL; // ne pas faire free() pointe sur une variable static
    return 0;

  case QUERY_INIT:
    switch(type){
    case ALK_NOR: Q->wrap=P_NOR; break;
    case ALK_CYC: Q->wrap=P_CYC; break;
    case ALK_ISO: Q->wrap=P_ISO; break;
    case ALK_NEO: Q->wrap=P_NEO; break;
    case ALK_SEC: Q->wrap=P_SEC; break;
    case ALK_TER: Q->wrap=P_TER; break;
    }
    if(n<Q->wrap[0]) RET_n(0);
    RET_n(3*n+2*(type!=ALK_CYC));

  case QUERY_NAME:
    strcpy(Q->name,(Q->i<n)?"C":"H");
    return 0;
    
  case QUERY_ADJ:;
    int i=Q->i,j=Q->j;
    if((i<n)&&(n<=j)){ // i dans C et j dans H
      if(b) RET_a(1); // cas du methane qui n'a qu'un seul C
      j-=n;
      if(j<3*Q->wrap[1]) RET_a(j/3==i); // deg(j)=1
      j-=3*Q->wrap[1],i-=Q->wrap[1]; 
      if(j<1*Q->wrap[2]) RET_a(j/1==i); // deg(j)=3
      j-=1*Q->wrap[2],i-=Q->wrap[2]+Q->wrap[3];
      RET_a(j/2==i);                    // deg(j)=2
    }
    if(j>=n) RET_a(0); // i et j >=n => pas d'arÃªte

    switch(type){ // i et j sont dans C

    case ALK_NOR:
      // 0-2-...-1
      RET_a((j-i==n-1)||((j==i+1)&&(i>=1)));

    case ALK_CYC:
      // 0-1-2-...-0
      RET_a((j-i==n-1)||(j==i+1));

    case ALK_ISO:
      // 0
      /*  \          */
      // 1-3-4-...-2
      RET_a( ((i<2)&&(j==3))||
	     ((i==2)&&(j==n-1))||
	     ((j==i+1)&&(i>=3))
	     );

    case ALK_NEO:
      // 0
      /*  \          */
      // 1-4-5-...-3
      //  /
      // 2
      RET_a( ((i<3)&&(j==4))||
	     ((i==3)&&(j==n-1))||
	     ((j==i+1)&&(i>=4))
	     );
      
    case ALK_SEC:
      //     4-...-1
      //    /
      // 0-3-a-...-2
      a=(n-4)/2+4;
      RET_a( ((i==0)&&(j==3))||
	     ((i==1)&&(j==a-1))||
	     ((i==2)&&(j==n-1))||
	     ((i==3)&&(j==a))||
	     ((j==i+1)&&(i>=3)&&(j<a))||
	     ((j==i+1)&&(i>=a))
	     );

    case ALK_TER:
      //     5-...-1
      //    /
      // 0-4-a-...-2
      /*    \        */
      //     b-...-3
      a=(n-5)/3+5;
      b=a+(n-a)/2;
      RET_a( ((i==0)&&(j==4))||
	     ((i==1)&&(j==a-1))||
	     ((i==2)&&(j==b-1))||
	     ((i==3)&&(j==n-1))||
	     ((i==4)&&(j==a))||
	     ((i==4)&&(j==b))||
	     ((j==i+1)&&(i>=4)&&(j<a))||
	     ((j==i+1)&&(i>=a)&&(j<b))||
	     ((j==i+1)&&(i>=b))
	     );
    }

  }
  return 1;
}


int mycielski(query* const Q)
/*
  Utilise i<j.
  Modifie Q->i et Q->j.
*/
{
  int ki,kj,b,k,t;
  switch(Q->code){

  case QUERY_INIT:
    k=Q->param[0];
    RET_n((k<2)? 0 : 3*(1<<(k-2))-1);

  case QUERY_ADJ:
    ki=ceil(log2((double)(Q->i+2)/3.));
    kj=ceil(log2((double)(Q->j+2)/3.));
    k=3*(1<<kj)-2; /* rem: k est pair */
    b=(Q->j==k);
    if(ki==kj) RET_a(b);
    if(b) RET_a(0);
    t=Q->j-(k/2);
    if(t==Q->i) RET_a(0);
    if(Q->i<t){ Q->j=t; return mycielski(Q); }
    Q->j=Q->i;Q->i=t;
    return mycielski(Q);

  }
  return 1;
}


int windmill(query* const Q)
{
  switch(Q->code){
  case QUERY_INIT: RET_n(2*Q->param[0]+1);
  case QUERY_ADJ: RET_a((Q->i==0)||((Q->i&1)&&(Q->j==Q->i+1)));
  }
  return 1;
}


int matching(query* const Q)
/*
  Utilise i<j, et marche en orientÃ©.
  On pourrait utiliser aussi: "fdrg 2n 1 ."
*/
{
  switch(Q->code){
  case QUERY_INIT: RET_n(2*Q->param[0]);
  case QUERY_ADJ: RET_a((Q->i==Q->j-1)&&(Q->j&1));
  }
  return 1;
}


int ring(query* const Q)
/*
  Marche en orientÃ©.
  Q->param[0]=k+1
  Q->param[1]=n
  Q->param[2..k+2]=cordes c_i>=-n
*/
{
  const int t=Q->param[0]+1;
  const int n=Q->param[1];
  int u;

  switch(Q->code){

  case QUERY_INIT:
    for(u=2;u<t;u++) if(Q->param[u]<-n) Erreur(6); // paramÃ¨tre incorrect
    RET_n(n);

  case QUERY_ADJ:;
    const int b=!(Q->directed);
    for(u=2;u<t;u++){
      if(Q->j==((Q->i+Q->param[u]+n)%n)) RET_a(1);
      if(b && (Q->i==((Q->j+Q->param[u]+n)%n))) RET_a(1);
    }
    RET_a(0);
  }

  return 1;
}


int cage(query* const Q)
/*
  Marche en orientÃ©.
  Q->param[0]=k+1
  Q->param[1]=n
  Q->param[2..k+2]=cordes c_i>=-n
*/
{
  const int k=Q->param[0]-1;
  const int n=Q->param[1];

  switch(Q->code){

  case QUERY_INIT:
    if(k<1) Erreur(6); // paramÃ¨tre incorrect
    for(int i=0;i<k;i++) if(Q->param[i+2]<-n) Erreur(6); // paramÃ¨tre incorrect
    RET_n(n);

  case QUERY_ADJ:
    if(Q->j==(Q->i+1)%n) RET_a(1); // cycle
    if(Q->j==(Q->i+Q->param[(Q->i%k)+2]+n)%n) RET_a(1); // corde
    if(Q->directed) RET_a(0);
    if(Q->i==(Q->j+1)%n) RET_a(1);
    if(Q->i==(Q->j+Q->param[(Q->j%k)+2]+n)%n) RET_a(1);
    RET_a(0);
  }

  return 1;
}


int crown(query* const Q)
/*
  Utilise i<j.
*/
{
  int k=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(2*k);
  case QUERY_ADJ: RET_a((Q->i<k)&&(Q->j>=k)&&(Q->i!=(Q->j-k)));
  }
  return 1;
}


int fan(query* const Q)
/*
  Utilise i<j.
  On pourrait utiliser aussi: "grid 1 p -apex q"
*/
{
  const int p=Q->param[0];
  const int q=Q->param[1];
  switch(Q->code){
  case QUERY_INIT: RET_n(p+q);
  case QUERY_ADJ: RET_a(((Q->j==Q->i+1)&&(Q->j<p))||((Q->i<p)&&(Q->j>=p)));
  }
  return 1;
}


int split(query* const Q)
/*
  Utilise i<j.
  On pourrait utiliser aussi: "ring p . -not -apex q"
*/
{
  switch(Q->code){
  case QUERY_INIT: RET_n(Q->param[0]);
  case QUERY_ADJ: RET_a((Q->i<Q->param[1])||(Q->j<Q->param[1]));
  }
  return 1;
}


int chess(query* const Q)
{
  const int p=Q->param[0];
  const int q=Q->param[1];
  
  switch(Q->code){

  case QUERY_NAME:
    name_base(Q->name,Q->i,Q->param[0],2,",","()",1);
    return 0;
    
  case QUERY_INIT: RET_n(p*q);
    
  case QUERY_ADJ:;
    const int x=Q->param[2];
    const int y=Q->param[3];
    const int xi=Q->i%p;
    const int yi=Q->i/p;
    const int xj=Q->j%p;
    const int yj=Q->j/p;  
    RET_a(((abs(xi-xj)==x)&&(abs(yi-yj)==y))||((abs(xi-xj)==y)&&(abs(yi-yj)==x)));
  }
  
  return 1;
}


int grid(query* const Q)
/*
  Modifie Q->i et Q->j.
*/
{
  int x,y,k,z,p,h,b;
  int const d=Q->param[0];

  switch(Q->code){

  case QUERY_NAME:
    z=1; /* z=vrai ssi toutes les dimensions sont < 10 */
    int R[DIMAX];
    for(k=0;k<d;k++){
      b=Q->param[k+1];
      R[k]=Q->i%b;
      Q->i /= b;
      z &= (b<11);
    }
    name_vector(Q->name,R,d,(z?"":","),(z?"":"()"),0,"%i");
    return 0;
    
  case QUERY_INIT:
    if(d<0) Erreur(6); // paramÃ¨tre incorrect
    if(d>DIMAX) Erreur(4);
    free(Q->wrap);ALLOC(Q->wrap,d);
    for(Q->n=1,k=0;k<d;k++){
      p=Q->param[k+1];
      Q->wrap[k]=(p<0);
      p=abs(p);
      Q->param[k+1]=p;
      Q->n *= p;
    }
    return 0;

  case QUERY_END:
    free(Q->wrap);Q->wrap=NULL;
    return 0;

  case QUERY_ADJ:
    z=h=k=b=0;
    while((k<d)&&(b<2)&&(h<2)){
      p=Q->param[k+1];
      x=Q->i%p;
      y=Q->j%p;
      h=abs(x-y);
      if(h==0) z++;
      if((h>1)&&(Q->wrap[k])&&(((x==0)&&(y==p-1))||((y==0)&&(x==p-1)))) h=1;
      if(h==1) b++;
      Q->i /= p;
      Q->j /= p;
      k++;
    }
    RET_a((b==1)&&(z==(d-1)));
  }

  return 1;
}


int klein(query* const Q)
/*
  Positionne les points selon une grille p x q.
*/
{
  int const p=abs(Q->param[0]);
  int const q=abs(Q->param[1]);

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
    
  case QUERY_INIT:
    if((p<1)||(q<1)) Erreur(6);
    SET_n(p*q);
    ALLOCZ(Q->xpos,Q->n,_i%q); // fixe les coordonnÃ©es X
    ALLOCZ(Q->ypos,Q->n,_i/q); // fixe les coordonnÃ©es Y
    XYtype=XY_USER; // coordonnÃ©es fixÃ©es par l'utilisateur
    return InitXY(Q); // pour les options -xy noise/scale ...

  case QUERY_NAME:;
    int R[]={ Q->i%q, Q->i/q };
    name_vector(Q->name,R,2,",","()",0,"%i");
    return 0;

  case QUERY_ADJ:;
    int const xi=Q->i%q;
    int const yi=Q->i/q;
    int const xj=Q->j%q;
    int const yj=Q->j/q;
    // adjacence grille
    if((xi==xj)&&(abs(yi-yj)==1)) RET_a(1);
    if((yi==yj)&&(abs(xi-xj)==1)) RET_a(1);
    // wrap ou twist
    int const p1=p-1;
    int const q1=q-1;
    int const t0=Q->param[0]<0; // Q->param[0] twist ?
    int const t1=Q->param[1]<0; // Q->param[1] twist ?
    if( (((xi==0)&&(xj==q1))||((xi==q1)&&(xj==0)))&&
	((t0&&(yi+yj==p1))||((!t0)&&(yi==yj))) ) RET_a(1);
    if( (((yi==0)&&(yj==p1))||((yi==p1)&&(yj==0)))&&
	((t1&&(xi+xj==q1))||((!t1)&&(xi==xj))) ) RET_a(1);
    RET_a(0);
  }

  return 1;
}


int clebsch(query* const Q)
{
  if((Q->code==QUERY_ADJ)&&((Q->i|Q->j)==Q->n-1)&&((Q->i&Q->j)==0)) RET_a(1);
  /* sommets opposÃ©s */
  return grid(Q);
}


int collatz(query* const Q)
/*
  Graphe orientÃ© contenant Ã©ventuellement des boucles.

  EntrÃ©es:
   param[0]=2k+1, k>=1
   param[1]=n, n<>0
   param[2+2i]=a_i, i=0..k-1
   param[3+2i]=b_i, i=0..k-1

  Sorties:
   rep[u]=successeur du sommet u
   wrap[u]=valeur du sommet u

  Le successeur d'un sommet u de valeur x est dÃ©fini par la
  transformation C(x) = (a_i*x+b_i)/k oÃ¹ i=x%k. On souhaite calculer
  l'image itÃ©rÃ©e par C de tous les entiers x=1..|n|. On s'arrÃªte
  lorsque tous les entiers gÃ©nÃ©rÃ©s ont un successeur dÃ©jÃ  calculÃ©.
  Comme le processus pourrait ne pas s'arÃªter, on s'arrÃªte si le
  nombre d'entiers gÃ©nÃ©rÃ©s dÃ©passent n^2.

  Le principe, pour n<0 (mode "volume"), est de maintenir une liste
  (une file implÃ©mentÃ©e par une paire de tableaux extensibles avec une
  stratÃ©gie doublante) des sommets avec valeur (tableau wrap) et
  successeur (tableau succ). Au dÃ©part on a wrap[u]=u+1 pour u=0..n-1
  en laissant succ[u] indÃ©terminÃ© Ã  -1. Puis pour chacun des sommets u
  de wrap on calcule w=C(wrap[u]). Si w est dÃ©jÃ  dans wrap, disons
  w=wrap[v], on pose simplement succ[u]=v. Sinon, on ajoute un nouveau
  sommet Ã  wrap qui sera Ã  traiter.

  Le principe est similaire pour n>0 (mode "intervalle"). La
  diffÃ©rence est que la liste est initialisÃ©e avec la valeur 1 et la
  transformation inverse Câ»Â¹(x) est considÃ©rÃ©e. On cherche donc les
  images inverses, si elles existent, pour chaque i=0..k-1. Il faut
  que k*x-b_i soit divisible par a_i sinon l'image inverse n'existe
  pas. On s'arrÃªte lorsque les n premiers sommets ont Ã©tÃ© crÃ©Ã©s.
*/
{
  switch(Q->code){

  case QUERY_NAME:
    name_base(Q->name,Q->wrap[Q->i],10,0,"","",0);
    return 0;

  case QUERY_END:
    free(Q->wrap);Q->wrap=NULL;
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:
    Q->k=1; // pour adjacency_rep()
    int u=Q->param[0]; // il faut u impair >= 3
    if((u<3)||((u&1)==0)) Erreur(6); // paramÃ¨tre incorrect
    const int k=(u>>1); // k=nombre de paramÃ¨tres a_i,b_i
    int n=Q->param[1],b=(n>=0),i,t,tmax,s,w,m;
    if(n==0) RET_n(0); // graphe vide

    // b = 1 ssi mode "volume", 0 si mode "intervalle"
    // t = taille courante des tableaux extensibles (wrap et succ)
    // tmax = taille maximale des tableaux
    // s = nombre de valeurs dans wrap
    
    if(b){
      s=1;
      t=tmax=n;
      for(i=0;i<k;i++) // vÃ©rifie que les a_i sont tous non-nuls
	if(Q->param[2+(i<<1)]==0) Erreur(6); // paramÃ¨tre incorrect
      for(i=0;i<k;i++) // vÃ©rifie que aáµ¢Â·i + báµ¢ â‰¡ 0 (mod k) pour tous les i
	if((Q->param[2+(i<<1)]*i+Q->param[3+(i<<1)])%k) Erreur(18);
    }else{
      s=-n;
      t=s<<1; // taille de wrap = 2|n| au dÃ©part
      tmax=n*n; // taille max = n^2
    }
    
    NALLOC(int,wrap,t); // wrap[u]=valeur du sommet u
    NALLOC(int,succ,t); // succ[u]=sommet successeur du sommet u
    // initialisation des s premiÃ¨res valeurs pour wrap et succ
    for(u=0;u<s;u++) wrap[u]=u+1, succ[u]=-1;
    
    n=i=0; // n = sommet courant de wrap Ã  traiter, nâˆˆ[0,s[
    while(n<s){ // reste-t'il des sommets Ã  traiter ?
      m=n; // copie du sommet courant n qui va Ãªtre incrÃ©mentÃ©
      
      // dÃ©termine la valeur w et incrÃ©mente n si nÃ©cessaire
      if(b){
	w=k*wrap[n]-Q->param[3+(i<<1)]; // w=k*x-b_i
	u=Q->param[2+(i<<1)]; // u=a_i
	i++; if(i==k) i=0, n++; // prochains coefficients et/ou prochain n
	if(w%u) continue; // si a_i ne divise pas w
	w /= u; // si oui, on divise, w = valeur du prÃ©dÃ©cesseur de n
      }else{
	i=(wrap[n]%k)<<1;
	w=(Q->param[2+i]*wrap[n]+Q->param[3+i])/k; // w = valeur du successeur de n
	n++;
      }
      if(w<=0) continue; // on saute les valeurs nÃ©gatives (out of range)
      for(u=0;(u<s)&&(w!=wrap[u]);u++); // cherche w dans wrap[0..s[, u=position
      if(u==s){ // w n'est pas dans wrap -> ajoÃ»te le sommet u Ã  wrap
	if(s==tmax) break; // nombre max de sommets attend
	if(s==t){ // s=taille de wrap -> double wrap (et succ)
	  t <<= 1;
	  REALLOC(wrap,t);
	  REALLOC(succ,t);
	}
	wrap[s]=w;  // initialise la valeur du nouveau sommet
	succ[s]=-1; // initialise son successeur (par dÃ©faut aucun)
	s++; // un sommet de plus dans wrap
      }
      // w est en position u dans wrap
      if(b) succ[u]=m; // le prÃ©dÃ©cesseur de n est u
      else succ[m]=u;  // le successeur de n est u
    }
  
    n=s; // n=nombre de sommets du graphe
    SET_n(n); // fixe le nombre de sommets du graphe
    ALLOC2(Q->rep,n,1);
    for(u=0;u<n;u++) Q->rep[u][0]=succ[u];
    free(succ);
    REALLOC(wrap,n);
    free(Q->wrap);
    Q->wrap=wrap; // pour calculer Q->name
    return 0;
  }

  return 1;
}


int rplg(query* const Q)
/*
  Q->param[0]=n
  Q->dparam[0]=t
  Q->dparam[1]=âˆ‘_i w_i
  Q->drep[i][0]=degrÃ© moyen du sommet i = (n/(i+1))^(1/(t-1))
*/
{
  switch(Q->code){

  case QUERY_END:
    free(Q->drep);
    return 0;

  case QUERY_ADJ:
    RET_a( (RAND01 < (Q->drep[Q->i][0]*Q->drep[Q->j][0]/Q->dparam[1])) );
    
  case QUERY_INIT:
    SET_n(Q->param[0]);
    const int t=Q->dparam[0];
    if(t<=1) Erreur(6); // paramÃ¨tre incorrect
    double const n=Q->n;
    double const c=1.0/(t-1.0);
    double s=0;

    ALLOC2(Q->drep,Q->n,1);
    for(int k=0;k<Q->n;k++) s += (Q->drep[k][0]=pow(n/((double)k+1.0),c));
    Q->dparam[1]=s;

    return 0;
  }

  return 1;
}


int butterfly(query* const Q)
/*
  Modifie Q->i et Q->j.
*/
{
  int d=Q->param[0]+1; /* d=dim+1 */
  switch(Q->code){

  case QUERY_INIT:
    d--;         /* d=dim */
    Q->n=(1<<d); /* n=2^dim */
    Q->n *= d+1; /* n=(dim+1)*2^dim */
    return 0;

  case QUERY_ADJ:;
    int x=Q->i/d;Q->i%=d; /* i -> (x,i) = (mot binaire,niveau) */
    int y=Q->j/d;Q->j%=d; /* j -> (y,j) = (mot binaire,niveau) */
  
    if(Q->j==Q->i+1) RET_a((x==y)||((x^y)==(1<<Q->i)));
    if(Q->i==Q->j+1) RET_a((x==y)||((x^y)==(1<<Q->j)));
    RET_a(0);
  }
  
  return 1;
}


int debruijn(query* const Q)
{
  const int d=Q->param[0];
  const int b=Q->param[1];
  int k;

  switch(Q->code){

  case QUERY_INIT:
    if((d<0)||(b<1)) RET_n(0); // il faut d>=0 et b>=1
    for(Q->n=1,k=0;k<d;k++) Q->n *= b;
    return 0;

  case QUERY_ADJ:;
    const int x=Q->j-(Q->i*b)%Q->n;
    const int y=Q->i-(Q->j*b)%Q->n;
    RET_a(((0<=x)&&(x<b))||((0<=y)&&(y<b)));
  }

  return 1;
}


int barbell(query* const Q)
/*
  Utilise i<j.
*/
{
  const int n1=abs(Q->param[0]);
  const int n2=abs(Q->param[1]);
  const int p=Q->param[2];
  switch(Q->code){

  case QUERY_INIT:
    RET_n(n1+n2+p-1);
  
  case QUERY_ADJ:
    if(Q->j<n1){
      if(Q->param[0]<0) RET_a((Q->j==Q->i+1)||((Q->i==0)&&(Q->j==n1-1))); /* cycle 1 */
      RET_a(1); /* clique 1 */
    }
    if(Q->i>=n1-1+p){
      if(Q->param[1]<0) RET_a((Q->j==Q->i+1)||((Q->i==n1-1+p)&&(Q->j==n1+n2+p-2))); /* cycle 2 */
      RET_a(1); /* clique 2 */
    }
    if((n1-1<=Q->i)&&(Q->j<n1+p)) RET_a((Q->j-Q->i==1)); /* chemin */
    RET_a(0);
  }

  return 1;
}


int shuffle(query* const Q)
{
  int n=Q->param[0];
  switch(Q->code){

  case QUERY_NAME:
    name_base(Q->name,Q->i,2,n,"","",0);
    return 0;

  case QUERY_INIT:
    RET_n(1<<n);
  
  case QUERY_ADJ:
    if((Q->i>>1)==(Q->j>>1)) RET_a(1);
    n=Q->n/2;
    if(Q->j==((Q->i<<1)%Q->n+(Q->i>=n))) RET_a(1);
    if(Q->j==((Q->i>>1)+((Q->i&1)?n:0))) RET_a(1);
    RET_a(0);
  }

  return 1;
}


int kautz(query* const Q)
/*
  A chaque sommet i, qui est entier de [0,b*(b-1)^(d-1)[, on associe
  un mot (x_1,...,x_d) codÃ©e sous la forme d'un entier rep[i][0] de
  [0,b^d[ (d lettres sur un alphabet Ã  b lettres). Alors i et j sont
  adjacents dans le Kautz ssi rep[i][0] et rep[j][0] sont adjacents
  dans le De Bruijn.

  param[] = d b x y oÃ¹
   x = #sommets du De Buijn = b^d
   y = #sommets du Kautz    = b*(b-1)^(d-1)
*/
{
  int d,b,k,u,q,r,s,t,x;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
      
  case QUERY_INIT:
    d=Q->param[0];
    b=Q->param[1];
    if((d<1)||(b<2)) RET_n(0); // il faut d>=1 et b>=2
    Q->n=x=b;
    t=b-1;
    for(k=1;k<d;k++){
      Q->n *= t;
      x *= b;
    }
    TEST_n; // fin si n<0
    Q->param[2]=x;    /* #sommets de De Bruijn */
    Q->param[3]=Q->n; /* #sommets de Kautz */
    ALLOC2(Q->rep,Q->n,1);

    for(u=0;u<Q->n;u++){ /* pour tous les sommets faire .. */
      /* On voit un sommet u comme (r,s2...sd) avec r dans [0,b[ et
	 s_i de [0,b-1[. On le convertit en (x1,...,xd) avec xd=r et
	 x_(d-1)=s2 ou s2+1 suivant si s2<r ou pas, etc. */
      r=x=u%b;
      q=u/b;
      for(k=1;k<d;k++){
	s=q%t;
	s+=(s>=r);
	r=s;
	q=q/t;
	x=x*b+r;
      }
      Q->rep[u][0]=x;
    }
    return 0;

  case QUERY_ADJ:
    Q->n=Q->param[2]; /* modifie le nb de sommets */
    Q->i=Q->rep[Q->i][0];
    Q->j=Q->rep[Q->j][0];
    debruijn(Q); /* adjacence Ã©crite dans Q->a */
    Q->n=Q->param[3]; /* rÃ©tablit le nb de sommets */
    return 0;
  }
  
  return 1;
}


int ggosset(query* const Q)
/*
  ggosset p d_1 v_1 ... d_k v_k
  Q->param[0..2k+2[ = 2k+1 p d_1 v_1 ... d_k v_k d

  On ajoute un dernier paramÃ¨tre Q->param[2+2*k]=d_1+...+d_k pour
  accÃ©lÃ©rer le test d'adjacence. Attention! l'ajout de ce paramÃ¨tre
  pourrait ne pas tenir dans le tableau Q->param.

  Les sommets sont toutes les permutations diffÃ©rentes du vecteur
  d'entiers (v_1,..,v_1, ..., v_k,..,v_k) et de son opposÃ©
  (-v_1,..,-v_1, ..., -v_k,..,-v_k) telles que le nombre de valeurs
  v_t est d_t. On pose Q->rep[i] = vecteur du sommet i qui est de
  taille d=d_1+...+d_k. On a l'arÃªte i-j ssi le produit scalaire entre
  Q->rep[i] et Q->rep[j] vaut p.

  Pour calculer tous les vecteurs possibles de taille d (les sommets)
  Ã  partir de (v_1...v_1,...,v_k...v_k) on procÃ¨de comme suit (codage
  d'un multi-ensemble Ã  k valeurs):

  On choisit les positions de v_1 (il y en a Binom(d,d_1) possibles),
  puis on choisit les positions de v_2 parmi les d-d_1 restantes (il y
  en a Binom(d-d_1,d_2)), puis les positions de v_3 parmi les
  d-(d_1+d_2) retantes, etc. Pour chaque t, les positions des d_t
  valeurs v_t sont codÃ©es par un sous-ensemble S[t] de [0,d[ de taille
  d_t.
*/
{
  int d,t,p,m,u,v,c,z,**S;

  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],Q->param[1+Q->param[0]],",","[]",1,"%i");
    return 0;

  case QUERY_ADJ:
    /* calcule le produit scalaire de Q->rep[i] et Q->rep[j] */
    d=Q->param[1+Q->param[0]]; // d=sum_i d_i
    for(p=t=0;t<d;t++) p += Q->rep[Q->i][t]*Q->rep[Q->j][t]; 
    RET_a( (p==Q->param[1]) );
    
  case QUERY_INIT:;
    const int k=Q->param[0]>>1;
    if(k<=0) RET_n(0);

    /* calcule d=âˆ‘_i d_i pour accÃ©lÃ©rer l'adjacence */
    for(d=t=0;t<k;t++){
      p=Q->param[2+2*t]; if(p<1) Erreur(6); // paramÃ¨tre incorrect
      d+=p;
    }
    Q->param[1+Q->param[0]]=d; // Ã©crit la somme Ã  la fin des paramÃ¨tres

    /* calcule Q->n = 2*âˆ_i binomial(d-(âˆ‘_{j<i} d_i),d_i) */
    m=d; Q->n=2; 
    for(t=2;m>0;t += 2){
      p=Q->param[t]; // p=d_i
      Q->n *= Binom(m,p);
      m -= p; // m=d-sum_{j<i}d_i
    }
    TEST_n; // fin si graphe vide
    ALLOC2(Q->rep,Q->n,d); // vecteur de taille d reprÃ©sentant les sommets
    NALLOC(int,P,d); // tableau intermÃ©diaire
    ALLOC2(S,k,d); // on rÃ©serve k tableaux (sous-ensembles) de taille <= d
    for(t=0;t<k;t++) NextSet(S[t],-1,d); // initialise les sous-ensembles
    // NB: taille |S[t]|=d_{t-1}=Q->param[2+2t] et v_{t-1}=Q->param[3+2t]

    for(u=0;u<Q->n;u+=2){
      /* Pour chaque sommet u on fait:

	 1. on remplit rep[u] et rep[u+1] Ã  partir des sous-ensembles S[0]...S[k-1]
	 2. on incrÃ©mente les sous-ensembles S[0]...S[k-1]
	 
	 Pour l'Ã©tape 1, il faut passer par un tableau intermÃ©diaire P
	 puis remplir rep[u] et rep[u+1] Ã  partir de P. Supposons d=5,
	 k=3, et S[0]={1,3} S[1]={1}, S[2]={0,1}.  On met dans P les
	 indices t des S[t] aux bons emplacements comme suit:

           - au dÃ©part P est vide: P={-,-,-,-,-}
	   - puis on ajoute S[0]:  P={-,0,-,0,-}
	   - puis on ajoute S[1]:  P={-,0,1,0,-}
	   - puis on ajoute S[2]:  P={2,0,1,0,2}
      */

      /* Calcule P */
      for(t=0;t<d;t++) P[t]=-1; /* tableau vide au dÃ©part */
      for(t=0;t<k;t++){ /* pour chaque sous-ensemble S[t] */
	m=-1;
	z=Q->param[2+2*t]; /* z=d_t */
	for(p=v=0;p<z;p++){ /* on parcoure S[t] */
	  /* mettre t dans la S[t][p]-Ã¨me case vide de P Ã  partir de l'indice v */
	  /* utilise le fait que dans S[t] les nombres sont croissant */
	  /* v=position courante de P oÃ¹ l'on va essayer d'Ã©crire t */
	  /* c=combien de cases vides de P Ã  partir de v faut-il encore sauter ? */
	  /* si P[v]=-1 et c=0 alors on Ã©crit P[v]=t */
	  c=S[t][p]-m;
	  m=S[t][p]; /* mÃ©morise l'ancienne valeur de S[t][p] */
	  while(c>0){
	    if(P[v]<0){
	      c--;
	      if(c==0) P[v]=t; /* Ã©crit t et passe Ã  la case suivante */
	    }
	    v++;
	  }
	}
      }

      /* Remplit rep[u] et rep[u+1] grÃ¢ce au tableau P (et incrÃ©menter u) */
      for(t=0;t<d;t++){
	v=Q->param[3+2*P[t]]; /* valeur v_t Ã  Ã©crire */
	Q->rep[u][t]=v;
	Q->rep[u+1][t]=-v;
      }

      /* IncrÃ©mente S[0]...S[k-1] grÃ¢ce Ã  NextSet() */
      t=0; /* on commence avec S[0] */
      m=d; /* S[t] dans {0...d-2} */
      v=2; /* Q->param[v] = taille de S[t] */
      while((t<k)&&(!NextSet(S[t],m,Q->param[v]))){
	t++; /* si S[t] fini on passe Ã  S[t+1] */
	m -= Q->param[v];
	v += 2;
      }
      /* si t=k alors on a atteint le dernier sommet */
    }

    free(P);
    free(S);
    return 0;
  }

  return 1;
}


int schlafli(query* const Q)
/*
  Utilise ggosset() et modifie Q->param[].

  Sous-graphe induit par le voisinage d'un sommet du graphe de Gosset
  qui a deux fois plus de sommets. Les voisins V du sommets 0 de
  Gosset sont, pour un sommet de Schlafi u=0..26 (liste obtenue par
  "gengraph gosset -format vertex 0"):

     u = 0 1 2 3 4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
     V = 2 4 6 8 11 12 14 17 19 20 22 25 27 29 30 32 35 37 39 41 42 44 47 49 51 53 55

  2u+2 = 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54
     d = 0 0 0 0 +1 0  0  +1 +1 0  0  +1 +1 +1 0  0  +1 +1 +1 +1 0  0  +1 +1 +1 +1 +1
       = 00  001    0011        00111          001111            0011111

  Ainsi V[u] est proche de 2u+2. La diffÃ©rence d[u] est 0 ou 1. PlutÃ´t
  que de stocker V[u], on pourrait calculer par V[u] = 2u+2+d[u]. La
  suite d[u] est particuliÃ¨re. Elle forme des blocs 00, 001, 0011,
  00111, ... de la forme 001^b pour b=0..5. La position du dÃ©but de
  ces bloc dans d[] vaut pos(b) = 0,2,5,9,14,20 pour b=0..5, soit
  encore pos(b) = (b+1)(b+2)/2-1. Du coup d[u] = ((u-pos(b))>1) oÃ¹ b
  est le plus grand entier tq pos(b) <= u.

  Comment calculer ce b ?
  1. x=2*u+2; // cherche le +grand b tq (b+1)*(b+2)<=x
  2. b=floor(sqrt(x))-1; // le +grand b tq (b+1)^2<=x (surestimation)
  3. if(x<(b+1)*(b+2)) b--; // si on avait pris b trop grand
*/
{
  switch(Q->code){

  case QUERY_END:
  case QUERY_NAME:
  case QUERY_ADJ:
    return ggosset(Q);

  case QUERY_INIT:
    /* paramÃ¨tres identiques Ã  gosset */
    Q->param[0]=5; Q->param[1]=8;
    Q->param[2]=2; Q->param[3]=3; Q->param[4]=6; Q->param[5]=-1;
    ggosset(Q); // calcule Q->rep[0..Q->n[, avec Q->n=56

    int u,x,b;
    for(u=0;u<27;u++){
      x=2*u+2;
      b=(int)(sqrt(x))-1;
      if(x<(b+1)*(b+2)) b--;
      x+=(u-(b+1)*(b+2)/2)>0; // x=(2u+2)+((u-pos(b))>1) = voisin dans Gosset
      Q->rep[u]=Q->rep[x]; // NB: u<x
    }
    // pas besoin de libÃ©rer les autres pointeurs Q->rep[26..56[ car
    // le free() de QUERY_END ne dÃ©pend pas de Q->n
    Q->n=27;
    return 0;

  }
  return 1;
}


int linial(query* const Q)
/*
  Q->rep[i][0..t[ = reprÃ©sentation du nom du sommet i = sa vue.
*/
{
  int const m=Q->param[0];
  int const t=Q->param[1];
  int k,u;
  
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    if(m<10) name_vector(Q->name,Q->rep[Q->i],t,"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],t,",","()",1,"%i");
    return 0;

  case QUERY_INIT:
    if((m<t)||(t<1)||(m<2)) RET_n(0); /* il faut m>=t>=1 et m>=2*/
    for(Q->n=1,k=m-t+1;k<=m;k++) Q->n *= k; /* calcul Q->n */
    ALLOC2(Q->rep,Q->n,t);
    NALLOC(int,S,t);
    NALLOC(int,P,t);
    NextArrangement(S,P,-1,t); /* initialisation de S et P */
    for(u=0;u<Q->n;u++){ /* gÃ©nÃ¨re tous les arrangements */
      for(k=0;k<t;k++) Q->rep[u][k]=S[P[k]];
      NextArrangement(S,P,m,t);
    }
    free(P);
    free(S);
    return 0;
    
  case QUERY_ADJ:
    if((Q->rep[Q->i][0]!=Q->rep[Q->j][t-1])||(m==t)){
      for(u=1,k=0;u<t;u++,k++) if(Q->rep[Q->i][u]!=Q->rep[Q->j][k]) break;
      if(u==t) RET_a(1);
    }
    if(Q->i<Q->j){ SWAP(Q->i,Q->j); return linial(Q); }
    RET_a(0);
  }

  return 1;
}


int linialc(query* const Q)
/*
  Comme Linial ou presque.
*/
{
  if(Q->code==QUERY_INIT){
    int const m=Q->param[0];
    int const t=Q->param[1];
    int const m1=m-1;
    if((m<t)||(t<0)||(m<2)) RET_n(0); /* il faut m>=t>=0 et m>=2*/
    int u,k,v,x,y;
    for(Q->n=m,k=1;k<t;k++) Q->n *= m1; /* calcule Q->n=m*(m-1)^t */
    ALLOC2(Q->rep,Q->n,t);
    /* on transforme u en (x0,x1,...x_t) avec x0 in [0,m[ et x_i in [0,m-1[ */
    for(u=0;u<Q->n;u++){
      x=Q->rep[u][0]=(u%m);
      for(v=u/m,k=1;k<t;k++){
	y = v%m1; /* y in [0,m-1[ */
	v /= m1;
	x=Q->rep[u][k]=y+(y>=x); /* si x=y, on incrÃ©mente y */
      }
    }
    return 0;
  }
  
  return linial(Q);
}


int gpstar(query* const Q)
/*
  Q->rep[i][0...n[ = reprÃ©sentation de la permutation du sommet i.
*/
{
  int const n=Q->param[0];
  int k,u;
  
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    if(n<10) name_vector(Q->name,Q->rep[Q->i],n,"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],n,",","()",1,"%i");
    return 0;

  case QUERY_INIT:
    if(n<=0) RET_n(0);
    for(Q->n=1,k=2;k<=n;k++) Q->n *= k;
    ALLOC2(Q->rep,Q->n,n); /* ici Q->n>0 */
    NALLOCZ(int,P,n,_i); /* initialise P */
    /* gÃ©nÃ¨re toutes les permutations */
    for(u=0;u<Q->n;u++){
      for(k=0;k<n;k++) Q->rep[u][k]=P[k]+1; /* copie P dans Q->rep */
      NextPermutation(P,n,NULL);
    }
    free(P);
    return 0;

  case QUERY_ADJ:
    /* distance de Hamming: on compte les diffÃ©rences entre les
       tableaux rep[i] et rep[j] */
    for(k=u=0;k<n;k++) u += (Q->rep[Q->i][k]!=Q->rep[Q->j][k]);
    RET_a((u==Q->param[1]));
  }
  
  return 1;
}


int pancake(query* const Q)
/*
  Utilise i<>j.
  Q->rep[i][0...n[ = reprÃ©sentation de la permutation du sommet i.
  On utilise Q->param[1] pour le signe de la permutation:
   s=+1 pour pancake()
   s=-1 pour bpancake()
*/
{
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],Q->param[0]," ","",1,(Q->param[1]<0)?"%+i":"%i");
    return 0;

  case QUERY_INIT:
    gpstar(Q); // gÃ©nÃ¨re les n! permutations dans Q->rep
    Q->param[1]=1; // signe>0
    return 0;

  case QUERY_ADJ:
    /* test d'adjacence Ã  partir de rep[i] et rep[j] */
    if(Q->i==Q->j) RET_a(0); // il est important d'avoir i<> j
    int const s=Q->param[1]; // signe pour le reversal
    int k=Q->param[0],z=0;
    do k--; while(Q->rep[Q->i][k]==Q->rep[Q->j][k]); // i<>j, donc on s'arrÃªte toujours
    while((k>=0)&&(Q->rep[Q->i][k]==s*Q->rep[Q->j][z])) k--,z++; // teste le "reversal" signÃ©
    RET_a((k<0)); // adjacent si prÃ©fixe=reversal
  }
  
  return 1;
}


int bpancake(query* const Q)
/*
  Q->rep[i][0...n[ = reprÃ©sentation de la permutation signÃ©e du sommet
  i. Il s'agit de n valeurs diffÃ©rentes (en valeur absolue) de
  {1,2,...,n,-1,-2,...,-n}. On Ã©vite soigneusement 0, car -0=+0.
*/
{
  if(Q->code==QUERY_INIT){

    Q->param[1]=-1; // signe<0
    int const n=Q->param[0];
    int k,p,q,t,c,u;
    if(n<=0) RET_n(0);
    for(t=1,k=2;k<=n;k++) t *= k; /* calcule t=n! */
    p=(1<<n); /* p=2^n */
    SET_n(p*t); /* Q->n=#sommets, forcÃ©ment >0 */
    ALLOC2(Q->rep,Q->n,n); /* permutations signÃ©es reprÃ©sentant les sommets */
    NALLOCZ(int,P,n,_i); /* initialise une permutation P (non signÃ©e) */

    /* GÃ©nÃ¨re toutes les permutations signÃ©es. On reprÃ©sente les
       signes par les bits de l'entier q=0...2^n-1. Si bit Ã  1 -> +1,
       sinon -> -1 */

    for(c=u=0;c<t;c++){ /* rÃ©pÃ¨te n! fois */
      for(q=0;q<p;q++,u++) /* rÃ©pÃ¨te 2^n fois */
	for(k=0;k<n;k++){ /* (1<<k)=mask=bit-0,bit-1,bit-2...,bit-(n-1) */
	  Q->rep[u][k]=P[k]+1; /* copie P dans Q->rep avec le signe +/-1 */
	  if(q&(1<<k)) Q->rep[u][k] *= -1; /* copie P dans Q->rep avec le signe +/-1 */
	}
      NextPermutation(P,n,NULL);
    }

    free(P);
    return 0;
  }

  return pancake(Q);
}


int pstar(query* const Q)
/*
  Q->rep[i][0...n[ = reprÃ©sentation de la permutation du sommet i.
  Presque comme gpstar() avec Q->param[1]=2.
*/
{
  if(Q->code==QUERY_ADJ){
    /* il faut deux diffÃ©rences dont le premier chiffre */
    Q->a=(Q->rep[Q->i][0]!=Q->rep[Q->j][0]);
    if(Q->a==0) return 0; // sinon comme gpstar() avec Q->param[1]=2
  }

  return gpstar(Q);
}
  

int gabriel(query* const Q)
/*
  L'initialisation et la terminaison sont communes Ã  beaucoup de
  graphe gÃ©omÃ©triques. L'adjacence est en O(n). Il est important de
  tester les n sommets (mÃªme ceux supprimÃ©s par -delv), sinon le
  rÃ©sultat n'est pas un sous-graphe.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
 
  case QUERY_INIT:
    SET_n(Q->param[0]);
    return InitXY(Q);

  case QUERY_ADJ:;
    /* c=(xc,yc)=milieu du segment i-j */
    int const n=Q->param[0];
    double const xc=(Q->xpos[Q->i]+Q->xpos[Q->j])/2;
    double const yc=(Q->ypos[Q->i]+Q->ypos[Q->j])/2;
    double const r=min(Norme_dxy(fabs(xc-Q->xpos[Q->i]),fabs(yc-Q->ypos[Q->i])),
		       Norme_dxy(fabs(xc-Q->xpos[Q->j]),fabs(yc-Q->ypos[Q->j])));
    /*
      r=dist(c,i)=rayon du disque centrÃ© au milieu de i-j. Attention !
      en gÃ©nÃ©ral r<>dist(i,j)/2 car ce n'est pas forcÃ©ment la norme
      L2. Ã€ cause des arrondis, il est possible que r soit surestimÃ©,
      et donc que dist(c,i)<>dist(c,j).  C'est pour cela qu'on prend
      r=min(dist(c,i),dist(c,j)).
    */

    int k;
    for(k=0;k<n;k++)
      if(Norme_dxy(fabs(xc-Q->xpos[k]),fabs(yc-Q->ypos[k]))<r) RET_a(0);
    RET_a(1);
  }

  return 1;
}
  

/*
  Une surface S de genre g orientÃ©e ou non, avec ou sans bord est
  modÃ©lisÃ©e par un polygone convexe rÃ©gulier P Ã  p=4g cotÃ©s numÃ©rotÃ©s
  de 0 Ã  p-1 de maniÃ¨re alternÃ©e comme suit (on parle de polygone
  fondamental https://en.wikipedia.org/wiki/Fundamental_polygon). Le
  cotÃ© 0 est positionnÃ© verticalement et le plus Ã  droit. Une fois le
  cotÃ© i positionnÃ©, le cotÃ© i+1 est, si i est pair, situÃ© deux cotÃ©s
  aprÃ¨s le cotÃ© i en tournant dans le sens direct (le sens contraire
  des aiguilles d'une montre), et si i est impair, situÃ© deux cotÃ©s
  avant. Voir l'exemple suivant pour g=2.

                                   1
                                _______
                               /       \
                            3 /         \ 2
                             /           \
                             |           |
                           4 |     S     | 0
                             |           |
                             \           /
                            6 \         / 7  
                               \_______/
                                   5    
 
  Les cotÃ©s sont appariÃ©s pour former les coutures de S. La surface S
  possÃ¨de donc 2g coutures. C'est la taille du tableau XYsurface
  codant la signature de S, c'est-Ã -dire la faÃ§on dont les coutures
  sont rÃ©alisÃ©es. Si l'on dÃ©coupe S le long de ses coutures, alors on
  obtient tout simplement P. Le cotÃ© i est appariÃ© au cotÃ© i+1 si i
  est pair, et au cotÃ© i-1 sinon. Dit autrement, le cotÃ© i est appariÃ©
  avec le cotÃ© i^1 (i xor 1 = inversion du bit de poids faible de
  i). Les coutures sont numÃ©rotÃ©es de 0 Ã  2g-1, le numÃ©ro de la
  couture du cotÃ© i Ã©tant simplement i>>1 (soit floor(i/2)). Les cotÃ©s
  appariÃ©s doivent avoir la mÃªme longueur. Pour simplifier, on
  supposera que tous les cotÃ©s de P sont de mÃªme longueur, soit
  2sin(ðœ‹/p) en supposant que P est inscrit dans un cercle de rayon
  unitÃ©.

  Chaque couture de numÃ©ro c correspond donc Ã  l'appariement des cotÃ©s
  2c et 2c+1. Elle possÃ¨de une signature XYsurface[c] qui peut prendre
  trois valeurs (+1,-1,0) codant les trois faÃ§ons de rÃ©aliser une
  couture. La signature >0 signifie que les cotÃ©s 2c et 2c+1 sont
  recollÃ©s sans inversion (handle) comme dans une surface
  orientable. La signature <0 signifie que les cotÃ©s 2c et 2c+1 sont
  recollÃ©s avec une inversion (crosscap) comme dans une surface
  non-orientable. Enfin, la signature =0 signifie que les cotÃ©s 2c et
  2c+1 forment un bord (ou trou) de la surface (border).

  On peut assembler des copies de P pour former un graphe de
  polygones. Ce graphe est a priori infini. Les sommets sont des
  copies de P recollÃ©s par les cotÃ©s appariÃ©s (les coutures) aprÃ¨s
  rotation, translation voir renversement. Un polygone Q a autant de
  voisins qu'il a de cotÃ© i qui ne correspond pas Ã  un bord (il faut
  donc que XYsurface[i>>1] soit non nul). On va voir ci-aprÃ¨s qu'on
  peut en fait limiter le graphe des polygones Ã  un arbre d'environ
  g^O(g) polygones. Un chemin sur la surface correspond donc Ã  un
  chemin dans le graphe des polygones, passant d'un polygone Ã  un
  polygone voisin chaque fois qu'il tranverse une couture.

  Pour obtenir le voisin R du cotÃ© i du polygone Q, on procÃ¨de comme
  suit (qui reste valide pour le disque hyperbolique):

  1) on fait une copie de Q, appelÃ©e R, qu'on supperpose Ã  Q;

  2) on fait une rotation de R correspondant Ã  deux cotÃ©s du polygone
     (4ðœ‹/p) dans le sens direct si i est impair, et inverse si i est
     pair, si bien que les cotÃ©s i de Q et i^1 de R sont confondus;

  3) on effectue une inversion de R autour du cotÃ© i de Q, soit une
     symÃ©trie par rapport au cotÃ© i. (Dans le cas hyperbolique, c'est
     une inversion par rapport Ã  le cercle reprÃ©sentant le cotÃ© i);

  4) si la signature de la couture i>>1 est nÃ©gative il faut de plus
     faire une symÃ©trie axiale de R selon la perpendiculaire passant
     par le cotÃ© i^1 (=renversement).

  Sans avoir de preuve, on va supposer qu'une gÃ©odÃ©sique d'un point u
  Ã  un point v sur la surface S correspond sur P (une fois S dÃ©coupÃ©e)
  Ã  une suite de segments de droite, chaque segment ayant comme
  extrÃ©mitÃ© u, v ou un point du pÃ©rimÃ¨tre de P. De plus, ces segments
  peuvent Ãªtre juxtaposÃ©s pour former une droite allant de u Ã  une
  copie de v dans l'un des polygones du graphe des polygones. On va
  aussi supposer qu'une gÃ©odÃ©sique ne coupe jamais deux fois la mÃªme
  couture. Ainsi on peut restreindre le graphe des polygones Ã  un
  arbre dont la racine est une copie de P et de profondeur <= 2g, avec
  la restriction que polygone n'a pas de fils ni avec les cotÃ©s
  correspondant avec une couture dÃ©jÃ  utilisÃ©e par un de ces ancÃªtres
  (en particulier par son pÃ¨re si ce n'est pas la racine).

  Ainsi, chaque point u de P se retrouve copiÃ© dans chacun des
  polygones de l'arbre. On repÃ¨re chaque copie par la suite des
  numÃ©ros des cotÃ©s permettant de passer d'un polygone Ã  son voisin,
  soit un chemin dans l'arbre des polygones. On supposera que les
  chemins sont valides, c'est-Ã -dire qu'il ne contient pas de cotÃ©
  correspondant Ã  un bord et que deux cotÃ©s appariÃ©s ne se suivent
  jamais.

  double surface_geodesic(point u,point v,int *C,point *D)

  // construit la gÃ©odÃ©sique allant du point u au point v(C). Renvoie
  // la longueur de la gÃ©odÃ©sique ou -1 si elle n'existe pas. Pour
  // qu'elle existe il faut que la droite entre u et v(C) traverse
  // effectivement dans l'ordre les cotÃ©s indiquÃ©s par C. On renvoie
  // Ã©galement dans le tableau D (qui sera supposÃ© de taille assez
  // grande) la liste des points d'intersections de u Ã  v(C).

  Ainsi pour trouver une gÃ©odÃ©sique entre u et v, il faut calculer une
  gÃ©odÃ©sique vers v(C) pour toutes les copies possibles de v, et de
  prendre la plus courte. Notons qu'elle existe bien toujours car la
  gÃ©odÃ©sique entre u et v(C) avec C={-1} existe toujours.

  Pour la construction d'une droite sur le disque poincarÃ©, voir:
  https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model#Compass_and_straightedge_construction
*/


int surface_next(int *C){
/*
  DÃ©termine le chemin valide (sans bord ni deux fois la mÃªme couture)
  immÃ©diatement aprÃ¨s le chemin C, vu comme un compteur qu'on essaye
  d'incrÃ©menter. Renvoie 1 si on a pu le faire, 0 sinon. Le chemin
  doit Ãªtre valide en entrÃ©e. Le chemin C identifie la copie d'un
  polygone en prÃ©cisant la suite (terminÃ©e par -1) des cotÃ©s qu'il
  faut traverser depuis le polygone racine pour y arriver.
  
  On utilise la fonction comme ceci:

    C[0]=-1; // initialisation du chemin
    do{
    ...; // traitement du chemin C
    }while(surface_next(C)); // chemin suivant

  Ex:  signature: bb      signature: hb      signature: hh
       1 chemin:          3 chemin:          13 chemins:
         C = -1             C = -1             C = -1
                            C = 0 -1           C = 0 -1
                            C = 1 -1           C = 1 -1
			                       C = 2 -1
			                       C = 3 -1
			                       C = 2 0 -1
					       C = 3 0 -1
					       C = 2 1 -1
					       C = 3 1 -1
					       C = 0 2 -1
					       C = 1 2 -1
					       C = 0 3 -1
					       C = 1 3 -1

  Comme on le voit, les chemins sont parcourus par longueur
  croissante, puis par valeur (comme un entier Ã©crit en base p=#cotÃ©
  du polygone) lu de droit Ã  gauche. Certains chemins gÃ©nÃ©rÃ©s peuvent
  mener au mÃªme polygone, comme [0 2] et [2 0] dans le cas d'un
  tore. Il y a donc 13 chemins valides, mais seulement 9 polygones
  diffÃ©rents. L'orientation de la surface n'a pas d'impact sur la
  sortie de la fonction. Seuls les bords en ont un.

  Pour amÃ©liorer la complexitÃ©, on utilise des variables statiques en
  supposant que les appels Ã  la fonction se suivent, c'est-Ã -dire
  qu'on ne change pas Ã  la main C[] (sinon en faisant C[0]=-1) entre
  deux appels Ã  la fonction. La complexitÃ© amortie est en 2^O(g), oÃ¹
  g=XYsurface, car pour visiter tous les environ g! chemins valides on
  Ã©numÃ¨re environ g^g chemins (compteurs en base g de taille g). En
  moyenne, cela fait donc g^g/g! â‰ƒ e^g pour chaque appel, Ã  des
  polynÃ´mes en g prÃ¨s. Comme il s'agit d'une permutation contrainte,
  on pourrait s'inspirer de l'algo de NextPermutation() qui serait
  bien plus efficace (en g! au lieu de g^g) et qui est le suivant: (il
  faudrait en fait raisonner sur les coutures plutÃ´t que les cotÃ©s)

  1. Trouver le plus grand index i tel que C[i] < C[i+1].
     S'il n'existe pas, la derniÃ¨re permutation est atteinte.
  2. Trouver le plus grand indice j tel que C[i] < C[j].
  3. Echanger C[i] avec C[j].
  4. Renverser la suite de C[i+1] jusqu'au dernier Ã©lÃ©ment.

  Autres points d'amÃ©liorations:

  1. Certains chemins mÃªmes au mÃªme polygone, et donc il n'ont pas
     besoin d'Ãªtre tous testÃ©s. Par exemple, pour le tore, [0 2] et [2
     0]. Il faudrait arriver Ã  n'en prendre qu'un seul ou le rendre
     canonique.

  2. Certaines successions de cotÃ©s ne peuvent correspondre Ã  un plus
     court chemin. La distance la plus grande dans un polygone inscrit
     dans un cercle de rayon 1 est 2. Cela interdit a priori des
     successions de cotÃ©s, comme la succession de trois cotÃ©s
     diamÃ©tralement opposÃ©s.
  
*/
  
  static int F[SURFACEMAX]; // F[c]=frÃ©quence des coutures de C[]
  static int ok=0; // ok=vrai ssi le chemin courant est valide
  static int p; // p=nombre de cotÃ©s du polygone, p>=4

  int i,j,c,z;

  if((C[0]<0)||(!ok)){ // initialisation de F[]
    for(c=i=0;c<XYsurfacesize;F[c++]=0); // met tout Ã  0
    while(C[i]>=0) F[C[i++]/2]=1; // met des 1 pour chaque couture (suppose C valide)
    ok=1; p=2*XYsurfacesize;
  }

  for(i=0;;){

    // ici on suppose qu'on a pas rÃ©ussit Ã  incrÃ©menter aucune des
    // valeurs d'indice < i, que F[] et ok sont Ã  jour. On essaye
    // d'incrÃ©menter C[i] de sorte que le nouveau cotÃ© C[i] ne soit ni
    // un bord ni corresponde Ã  une couture dÃ©jÃ  utilisÃ©e.
    
    if(C[i]<0){ // on arrive Ã  la fin du chemin courant
      if(i==XYsurfacesize) return ok=0; // on a tout explorÃ©
      C[i+1]=-1; // il reste Ã  mettre Ã  jour C[i] qui vaut ici -1
    }
    
    // ici on cherche le plus petit cotÃ© aprÃ¨s C[i] (modulo p) qui
    // n'est pas un bord. On passe en revue tous les cotÃ©s possibles.
    // Il faut rÃ©pÃ©ter le for(j=...) p fois si C[i]=-1, sinon p-1 fois
    // pour ne pas retomber sur C[i]. NB: on entre toujours dans le
    // for(j=...) car p>1, et Ã  la fin de la boucle, on a toujours
    // C[i]<>c.
    
    c=C[i]; // c=C[i] initial
    z=1; // z=vrai si le cotÃ© trouvÃ© est > c (sinon on est passÃ© par 0)
    for(j=(C[i]>=0);j<p;j++){ // rÃ©pÃ¨te p ou p-1 fois
      C[i]++;
      if(C[i]==p) C[i]=z=0; // il faudra incrÃ©menter C[i+1] (retenue)
      if(XYsurface[C[i]/2]) break; // on a trouvÃ© un cotÃ© qui n'est pas un bord 
    }
    if(j==p) return ok=0; // on a pas trouvÃ© de cotÃ© valide

    // on met Ã  jour F[] et ok car on a modifiÃ© C[i]
    if(c>=0){ // enlÃ¨ve la couture initiale, si elle existe
      F[c/2]--; // NB: au dÃ©part F[c/2]>0, forcÃ©ment
      // met Ã  jour ok qui peut rester constant ou passer de faux Ã  vrai
      if(!ok){ ok=1; j=0; while(C[j]>=0) ok &= (F[C[j++]/2]<=1); }
    }
    F[C[i]/2]++; // ajoute la nouvelle couture
    ok &= (F[C[i]/2]<=1); // met Ã  jour ok qui peut rester constant ou passer Ã  faux

    if(z){ // on a trouvÃ© le cotÃ© sans boucler (incrÃ©ment sans retenue)
      if(ok) return 1; // on a terminÃ©
      i=0; // on aurait du finir ici, mais ok=0, donc on recommence
    }else i++; // on a trouvÃ© un cotÃ© mais il y a eut une retenue (-> i+1)
  }
  
}


double surface_geodesic(point u,point v,int *C)
{
  u.x=v.x=C[0]; // pour Ã©viter le Warning Ã  la compilation
  return -1;
}


point surface_image(point u,int *C)
/*
  Donne les coordonnÃ©es du point u(C) correspondant Ã  la copie de u en
  suivant le chemin C dans l'arbre des polygones. C est une suite
  d'entiers terminÃ©e par -1 et supposÃ© valide. En particulier elle ne
  contient aucun cotÃ© qui est un bord. L'algorithme est en complexitÃ©
  O(|C|).

  Algorithme: On se dÃ©place de polygone en polygone suivant le chemin
  C en construisant le chemin dÃ©crit par les centres. Le polygone
  courant est repÃ©rÃ© par trois Ã©lÃ©ments:

    (1) son centre (x,y);
    (2) son dÃ©calage d (d=numÃ©ro de cotÃ© du polygone le + Ã  droit);
    (3) son orientation s (s=1 pour le sens direct, s=-1 sinon).

  Au dÃ©part, x=y=0, d=0, et s=1. L'objectif est donc de mettre Ã  jour
  le centre, le dÃ©calage et l'orientation pour le polygone voisin en
  traversant le cotÃ© i. Une fois le polygone final calculÃ© (centre,
  dÃ©calage, orientation), on calcule les coordonnÃ©es de v dans ce
  polygone lÃ . Pour cela on effectue une rotation dont l'angle dÃ©pend
  du dÃ©calage d et dont le sens (direct ou non) dÃ©pend de
  l'orientation s.

  Pour trouver de numÃ©ro de secteur j correspondant au cotÃ© de numÃ©ro
  i, c'est-Ã -dire l'entier j de [0,2p[ tel que le cotÃ© i est la corde
  du cÃ´ne d'angle compris entre j*2ðœ‹/p et (j+1)*2ðœ‹/p, il suffit
  d'Ã©changer les deux derniers bits de i:

                         cotÃ©    T   secteur

                         ...00  <->  ...00
			 ...01  <->  ...10
			 ...10  <->  ...01
			 ...11  <->  ...11

  Et bien sÃ»r la mÃªme transformation, notÃ©e T, permet de passer du
  numÃ©ro de secteur j au numÃ©ro de cotÃ© i. On a donc i=T(T(i)).
  L'Ã©change des deux derniers bits de la variable i peut se faire
  ainsi:

     MÃ©thode 1: b=i&3; j=i+(b==1)-(b==2);
     MÃ©thode 2: int const T[]={0,-1,1,0}; j=i+T[i&3];

  Pour calculer le nouveau dÃ©calage d' du polygone P' voisin par le
  cotÃ© i d'un polygone P de dÃ©calage courant d et d'orientation s, on
  peut faire ainsi. On imagine qu'on se dÃ©plasse du centre de P au
  centre de P' en traversant le cotÃ© i de P. On arrive alors par le
  cotÃ© i'=i^1 de P'. Si j est un numÃ©ro de secteur, on note opp(j) le
  numÃ©ro de secteur opposÃ© Ã  j. Il ne dÃ©pend pas de l'orientation.
  Bien sÃ»r on a opp(j)=(j+p/2)%p oÃ¹ p est le nombre de cotÃ©s des
  polygones. Le secteur correspondant Ã  d' vaut alors:

          ( opp(T(i')) + T(d)-T(i) + p/2 )%p  si s'=s
	  ( opp(T(i')) + T(i)-T(d) + p/2 )%p  sinon

  Pour les dÃ©calages (d), on a donc a priori intÃ©rÃªt Ã  travailler avec
  les numÃ©ros de secteur plutÃ´t que les numÃ©ros de cotÃ©.
*/
{
  static int const T[]={0,-1,1,0}; // transformation cotÃ© <-> secteur
  int const p=2*XYsurfacesize; // p=nombre de cotÃ© du polygone
  double const a0=M_2PI/p; // a0=angle du segment reliant deux centres de polygones voisins
  double const r0=2*cos(a0/2); // r0=distance entre deux centres de polygones voisins

  double a;
  int t,j,i;

  // configuration du polygone courant
  double x=0,y=0; // x,y=position du centre
  int d=0; // d=numÃ©ro de cotÃ© vertical le plus Ã  droite
  int s=1; // s=orientation

  // calcule la position du centre du polygone dÃ©terminÃ© par le chemin C
  for(t=0;C[t]>=0;t++){ // pour chaque cotÃ© C[t]
    i=C[t]; // i=cotÃ© courant
    j=i+T[i&3]; // j=numÃ©ro de secteur
    a=s*a0*(j-d); // a=angle du segment reliant le centre courant au prochain
    x += r0*cos(a);
    y += r0*sin(a);
    s *= -XYsurface[t];
    // A FINIR: mise Ã  jour de d ? (voir plus haut)
  }

  // A FINIR: il reste Ã  ajouter le point v au polygone final
  point v={x,y};
  u.x += x;
  u.y += y;
  return v;
}


int sgabriel(query* const Q)
/*
  Graphe de Gabriel dÃ©fini sur une surface, elle-mÃªme dÃ©finie par
  l'option -xy surface. En partie codÃ©e par Louis DÃ©moulins
  (stagiaire de Licence 2 en juin/juillet 2016).
*/
{
  if(XYsurfacesize==0) return gabriel(Q);
  // ici genre>0
  
  switch(Q->code){

  case QUERY_END:
  case QUERY_INIT:
    return gabriel(Q);

  case QUERY_DOT:
    USERDOT.adj=sgabriel;
    return 0;

  case QUERY_ADJ:;
    int const n=Q->param[0];
    double const dmax=Norme_dxy(2,2); // distance max dans la surface [-1,-1]x[+1,+1]
    const point pi={Q->xpos[Q->i],Q->ypos[Q->i]}; // pi=point i
    const point pj={Q->xpos[Q->j],Q->ypos[Q->j]}; // pj=point j
  
    //USERDOT.i=i;
    //USERDOT.j=j;

    NALLOC(int,cj,XYsurfacesize+1);
    NALLOC(int,cz,XYsurfacesize+1);

    double dz,d; // distances
    double r; // r=distance entre milieu m de pi-pj et pi
    point m,pz;
    int z;
    cj[0]=-1; // cj=chemin pour la copie de pj

    do{ // pour toutes les copies cj de j faire:
      z=0; // il faut dÃ©finir z Ã  cause du "continue"
      if(surface_geodesic(pi,pj,cj)<0) continue; // la gÃ©odÃ©sique n'existe pas
      m=surface_image(pj,cj); // m=pj(cj)
      m.x=(m.x+pi.x)/2, m.y=(m.y+pi.y)/2; // m=milieu entre pi et pj(cj)
      r=Norme_dxy(fabs(m.x-Q->xpos[Q->i]),fabs(m.y-Q->ypos[Q->i])); // r=distance(pi,m)
      // il faut vÃ©rifier qu'aucun sommet est dans le disque de rayon r autour de m
      for(;z<n;z++){ // NB: ici z=0
	if(z==Q->j) continue;
	// cherche la copie de z la plus proche de m
	pz.x=Q->xpos[z], pz.y=Q->ypos[z]; // pz=point d'indice z, cz=chemin de z
	cz[0]=-1; // cz=chemin pour le copie de pz
	dz=dmax; // dz=distance entre le meilleur z et m
	do{
	  d=surface_geodesic(m,pz,cz);
	  if(d>=0) dz=min(dz,d); // met Ã  jour la distance
	}while(surface_next(cz));
	if(dz<r) break; // pz est dans le disque -> on sort du for(z=...) avec z<n
      }
    }while((z<n)&&(surface_next(cj))); // si z<n, il faut essayer une autre copie de pj
    
    // ici:
    // z=n (-> 1, car il y a aucun point dans le disque de rayon r)
    // z<n (-> 0, un des points est dans le disque)
    
    free(cz), free(cj);
    RET_a((z==n));
  }

  return 1;
}


int mst(query* const Q)
/*
  Utilise load().
  GÃ©nÃ¨re directement un graphe Q->G (qui est l'arbre T).

  Algorithme pour calculer l'arbre T (=Q->G):
   1. remplir puis trier le tableau d'arÃªtes du graphe complet
   2. rÃ©pÃ©ter pour chaque arÃªte u-v, mais pas plus de n-1 fois:
      si u-v ne forme pas un cycle dans T (<=> u,v dans des composantes diffÃ©rentes)
      alors ajouter u-v au graphe T.

  On final, l'algorithme d'initialisation est en temps O(n^2*log(n)) Ã 
  cause du tri, et en espace O(n^2) Ã  cause du graphe qu'on
  remplit. On pourrait optimiser cette Ã©tape de remplissage, mais cela
  reste au moins en n^2*log(n) Ã  cause du tri.
*/
{
  switch(Q->code){

  case QUERY_ADJ:
    return load(Q);
    
  case QUERY_END:
    return free_rep(Q)||load(Q);

  case QUERY_INIT:;
    int n=Q->param[0];
    SET_n(n);
    InitXY(Q); // gÃ©nÃ¨re les points

    NALLOC(edge,E,n*(n-1)/2); // tableau d'arÃªtes
    int x,y,t=0;
    
    for(Q->i=0;Q->i<n;Q->i++) // remplit E[]
      for(Q->j=Q->i+1;Q->j<n;Q->j++,t++){
	E[t].u=Q->i;
	E[t].v=Q->j;
	E[t].w=dist_ij(Q);
      }
    QSORT(E,t,fcmp_edge); // trie les arÃªtes de E[] suivant leur poids

    NALLOCZ(int,parent,n,_i); // pour UNION-FIND
    NALLOCZ(int,rang,n,0);

    Q->G=new_fullgraph(n); // graphe complet avec G->d[u]=0 pour ADD_EDGE
    n--; // = nombre d'arÃªtes qui restent Ã  ajouter Ã  T
    t=0; // indice dans E
    
    while(n){ // tantqu'il reste une arÃªte Ã  ajouter
      x=UF_Find(E[t].u,parent);
      y=UF_Find(E[t].v,parent);
      if(x!=y){
	UF_Union(x,y,parent,rang);
	ADD_EDGE(Q->G,E[t].u,E[t].v); // ajoute u-v et met Ã  jour la taille des listes
	n--;
      }
      t++;
    }
    
    // libÃ¨re les tableaux
    free(parent);
    free(rang);
    free(E);
    GraphRealloc(Q->G,Q->G->d);

    return 0;
  }

  return 1;
}
 
 
int pat(query* const Q)
/*
  Graphe issu du jeu de Pat Morin. Un sommet correspond Ã  un sommet
  (x,y) d'une des k grilles. On suppose que i<j, ce qui revient Ã  dire
  que la grille de i est placÃ©e avant ou est Ã©gale Ã  celle de j.

  Exemple: p=q=3 et k=4

  06 07 08  15 16 17  24 25 26  33 34 35
  03 04 05  12 13 14  21 22 23  30 31 32
  00 01 02  09 10 11  18 19 20  27 28 29

  Grille 0  Grille 1  Grille 2  Grille 3

*/
{
  int const p=Q->param[0];
  int const q=Q->param[1];
  int const r=Q->param[2]; // r=round
  int const pq=p*q;

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
    
  case QUERY_INIT:
    SET_n(pq*r);
    if((p<=0)||(q<=0)||(r<=0)) RET_n(0);
    
    /* DÃ©termine les coordonnÃ©es des points pour en controler le
       dessin. Les grilles sont placÃ©es Ã  des hauteurs de plus en plus
       grandes pour "voir" les arÃªtes. Au dÃ©part (z=0), la premiÃ¨re
       grille est mise Ã  une hauteur 0, puis Ã  une hauteur 1, puis Ã 
       une hauteur 3, puis Ã  une hauteur 6, etc. La hauteur de la
       grille pour z quelconque est z(z+1)/2. En fait, on utilise
       z*(z+1)/2.5 pour une meilleure lisibilitÃ©. */

    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);
    int u,x,y,z;
    double h;
    
    for(u=0;u<Q->n;u++){// pour tous les sommets, faire:
      z=u/(pq);      // z=grille
      x=u%p;         // x=colonne
      y=(u%(pq))/p;  // y=ligne
      h=z*(z+1)/2.5; // h=dÃ©calage vers le haut de la grille z
      Q->xpos[u]=(double)(x+z*p)/(double)max(p,q);
      Q->ypos[u]=(double)(y+h*q)/(double)max(p,q);
    }
    
    XYtype=XY_USER; // coordonnÃ©es fixÃ©es par l'utilisateur
    return InitXY(Q); // pour les options -xy noise/scale ...
    
  case QUERY_ADJ:;
    // z=grille
    int const zi = Q->i/pq;
    int const zj = Q->j/pq;
    
    /* ici: zi<=zj car i<j */
    
    // x=colonne
    int const xi = Q->i%p;
    int const xj = Q->j%p;
    
    // y=ligne
    int const yi = (Q->i%pq)/p;
    int const yj = (Q->j%pq)/p;
    
    if(zj==zi) RET_a( ((xj>=xi)&&(yj<=yi)) || ((xj<=xi)&&(yj>=yi)) );
    /* ici: zj>zi */
    
    RET_a( ((xj>=xi)&&(yj==yi)) || ((xj==xi)&&(yj>=yi)) );
  }
  
  return 1;
}


int uno(query* const Q)
/*
  Utilise Q->xpos,Q->ypos pour les coordonnÃ©es.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int const p=Q->param[1];
    int const q=Q->param[2];
    if((p<=0)||(q<=0)) RET_n(0);
    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);
    
    for(int u=0;u<Q->n;u++){ // inverse (x,y) <-> (j,i)
      Q->xpos[u]=randomu(q);
      Q->ypos[u]=randomu(p);
    }
    
    XYtype=XY_USER; // coordonnÃ©es fixÃ©es par l'utilisateur
    return InitXY(Q); // pour les options -xy noise/scale ...
    
  case QUERY_ADJ:
    RET_a((Q->xpos[Q->i]==Q->xpos[Q->j])||(Q->ypos[Q->i]==Q->ypos[Q->j]));
  }

  return 1;
}
  
  
int unok(query* const Q)
/*
  Tout comme uno() sauf pour QUERY_INIT. Important: Les paramÃ¨tres ne
  servent que pour QUERY_INIT.
*/
{
  if(Q->code==QUERY_INIT){
    SET_n(Q->param[0]);
    int const p=Q->param[1];
    int const q=Q->param[2];
    int const kp=(Q->param[3]<0)?p:Q->param[3];
    int const kq=(Q->param[4]<0)?q:Q->param[4];
    if((p<1)||(q<1)||(kp<1)||(kq<1)||(Q->n>min(p*kp,q*kq))) RET_n(0);
    
    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);

    NALLOC2Z(int,G,p,q,0); // G[i][j]=1 ssi on a mis un sommet en (i,j)
    NALLOCZ(int,A,p,0);    // A[i]=nombre de sommets de G placÃ©s dans la ligne i
    NALLOCZ(int,B,q,0);    // B[j]=nombre de sommets de G placÃ©s dans la colonne j
    int u,z,c,i,j,t;

    /* Principe: 

       Soit z le nombre de cases libres oÃ¹ l'on peut placer un point,
       c'est-Ã -dire telles que G[i][j]=0, A[i]<kp et B[j]<kq. On
       choisit une des z cases libres avec c=randomu(z), puis on
       parcoure la grille, et quand la case libre numÃ©ro c arrive on y
       met un nouveau sommet (i,j), et on met Ã  jour z. Pour la mise Ã 
       jour on vÃ©rifie lorsque la ligne ou la colonne devient saturÃ©e
       en comptant les cases libres sur la croix de centre (i,j). La
       complexitÃ© en temps en environ npq/4 (en moyenne le point
       alÃ©atoire se situe en p*q/2 et la somme (n-i)*p*q/2 donne du
       n*p*q/4).

       On optimise cet algorithme en commenÃ§ant d'abord par un
       exÃ©cuter un algorithme Ã  rejets, qui ne nÃ©cessite pas le
       parcours de toute la grille. Il est rapide au dÃ©but puis va
       ralentir en fonction de la densitÃ© des cases libres de la
       grille et des rejets qui se produisent. On estime alors le
       temps dÃ©pensÃ© par cet algorithme depuis le dÃ©part (t1) et celui
       restant si l'on continuait par l'algorithme classique (t2).
       Tant que t1<t2 on applique l'algorithme Ã  rejets. Ensuite on
       change pour l'algorithme classique.
    */
    
    // partie commune aux deux algorithmes: ajoÃ»t du point(i,j) dans
    // la grille avec la mise Ã  jour des positions Q->xpos (=j), Q->ypos
    // (=i), de z, des vecteurs A et B. Les instructions I1 et I2
    // servent Ã  la mise Ã  jour de l'estimation du temps pour
    // l'algorithme Ã  rejet (non utilisÃ© pour l'algorithme classique).
#define UPDATE(I1,I2)							\
    do{									\
      G[i][j]=1,z--;							\
      A[i]++; if(A[i]==kp){ I1; for(t=0;t<q;t++) z -= (B[t]<kq)&&(G[i][t]==0); } \
      B[j]++; if(B[j]==kq){ I2; for(t=0;t<p;t++) z -= (A[t]<kp)&&(G[t][j]==0); } \
      Q->xpos[u]=j,Q->ypos[u]=i;					\
    }while(0)
    
    u=0;           // u=nombre de sommets dÃ©jÃ  tirÃ©s
    z=c=p*q;       // z=nombre de cases encore libres
    c/=4;          // c=temps moyen/sommet du temps de l'algo classique
    int t1=0;      // t1=temps dÃ©pensÃ© par l'algo Ã  rejet
    int t2=Q->n*c; // t2=temps estimÃ© restant pour l'algo classique

    // Algorithme Ã  rejets
    do{
      i=randomu(kp);
      j=randomu(kq);
      if((A[i]<kp)&&(B[j]<kq)&&(G[i][j]==0)){
	UPDATE(t1+=q,t1+=p);
	u++; if(u==Q->n) break;
	t2 -= c; // le temps de l'algo classique diminue
      }
    }while(t1++<t2);

    // Algorithme classique
    for(;u<Q->n;u++){
      c=randomu(z); // c=numÃ©ro de case libre alÃ©atoire uniforme parmi les libres
      for(i=0;i<p;i++){
	if(A[i]<kp) // sinon aucune case libre dans cette ligne
	  for(j=0;j<q;j++)
	    if((B[j]<kq)&&(G[i][j]==0)){
	      if(c==0){ // on a trouvÃ© la case libre alÃ©atoire
		UPDATE(,);
		i=j=p+q; // termine les deux for()
	      }else c--; // attend de trouver la case libre c
	    }
      }
    }
#undef UPDATE
    
    free(G); // /!\ allouÃ© par ALLOC2
    free(A);
    free(B);
    XYtype=XY_USER; // coordonnÃ©es fixÃ©es par l'utilisateur
    return InitXY(Q); // pour les options -xy noise/scale ...
  }
  
  return uno(Q);
}


/*
  Pour le dÃ©buggage de wpsl(). Affiche tous les coins et les blocs et
  vÃ©rifient la cohÃ©rence.
*/
#define wpsl_DEBUG							\
  do{									\
    int e,u,v,i,j,c1,c2,s0,s1;						\
    string s="/!\\ impossible /!\\";					\
    printf("COINS: %i\n",zn);						\
    for(u=e=0;u<zn;u++){						\
      printf("%i. (%i,%i) ",u,Z[u].d[0],Z[u].d[1]);			\
      for(i=k=e=0;i<4;i++){						\
	if(Z[u].B[i]) printf("%i ",(int)(Z[u].B[i]-T));			\
	else printf("N "),e++;						\
      }									\
      if(e==4) printf("%s\n",s);					\
      PRINTN;								\
      for(v=u+1,e=1;v<zn;v++)						\
	e=e&&((Z[u].d[0]!=Z[v].d[0])||(Z[u].d[1]!=Z[v].d[1]));		\
      if(e==0) printf("%s\n",s);					\
      for(i=j=0;i<4;i++){						\
	if(Z[u].B[i]==NULL) continue;					\
	j=1-j;e=0;							\
	if(i==0) s0=+1,s1=+1;						\
	if(i==1) s0=-1,s1=+1;						\
	if(i==2) s0=-1,s1=-1;						\
	if(i==3) s0=+1,s1=-1;						\
	if(s0*Z[u].B[i]->C[i]->d[0]>=s0*Z[u].d[0]) e=1;			\
	if(s1*Z[u].B[i]->C[i]->d[1]>=s1*Z[u].d[1]) e=2;			\
	if(!( ((s0*Z[u].B[i]->C[i^2]->d[0]>=s0*Z[u].d[0])&&		\
	       (Z[u].B[i]->C[i^2]->d[1]==Z[u].d[1]))||			\
	      ((s1*Z[u].B[i]->C[i^2]->d[1]>=s1*Z[u].d[1])&&		\
	       (Z[u].B[i]->C[i^2]->d[0]==Z[u].d[0])) )) e=3;		\
	if(s1*Z[u].B[i]->C[i^2]->d[1-j]<s1*Z[u].d[1-j]) e=4;		\
	if(e) printf("i=%i, j=%i, e=%i: %s\n",i,j,e,s);			\
      }									\
    }									\
    printf("BLOCS: %i\n",tn);						\
    for(u=0;u<tn;u++){							\
      printf("%i. sum=%lu ",u,T[u].sum);				\
      printf("left=");							\
      if(T[u].left) printf("%i ", (int)(T[u].left-T));			\
      else printf("N ");						\
      printf("right=");							\
      if(T[u].right) printf("%i ",(int)(T[u].right-T));			\
      else printf("N ");						\
      for(i=0,e=1,j=1;i<4;i++,j=1-j){					\
	if(T[u].C[i]==NULL){ printf("N "),e=0;continue; }		\
	printf("(%i,%i) ",T[u].C[i]->d[0],T[u].C[i]->d[1]);		\
	c1=(i+1)&3, c2=(i+2)&3; if(i%3) SWAP(c1,c2);			\
	e=e&&(T[u].C[c1]->d[j]<T[u].C[c2]->d[j]);			\
	e=e&&(T[u].C[c1]->d[1-j]==T[u].C[c2]->d[1-j]);			\
      }									\
      PRINTN;								\
      if(e==0) printf("%s\n",s);					\
    }									\
  }while(0)


/* Structure "bloc" pour wpsl():
   
   Un bloc est dÃ©fini par ses 4 coins, qui sont des points de la
   grille pxq. Ci-dessous un bloc qui a Ã©tÃ© subdivisÃ© deux fois: au
   point 1, puis au point 2. La surface ou l'aire du bloc (=nombre de
   cases), avant subdivision, est de 12=3x4. Le bloc ne contient que
   2x3=6 points internes. Une fois le point 1 choisi, il ne reste pour
   que 2 possibilitÃ©s pour le point 2. Si le bloc est dans une grille
   pxq, sa largeur sera dx=q-1, sa hauteur dy=p-1 et sa surface
   dx*dy. Le coin de coordonnÃ©es minimum du bloc est notÃ© C[0] (point
   en bas Ã  gauche). Les autres coins C[1..3] sont obtenus en tournant
   dans le sens direct ce qui permet de passer au suivant simplement
   (+1 mod 4). Ã€ cause de la symÃ©trie des axes, et pour simplifier
   l'implÃ©mentation, les axes X et Y sont codÃ©s plutÃ´t par un tableau
   de dimensions: d[0] pour X, et d[1] pour Y.

          ^ Y=d[1]           C[3]â”€â”€â”€â”€â”€C[2]
	  â”‚		       â”‚â—â”‚â—|â—|â—â”‚
	  â”‚		       â”‚â”€2â”€â”€â”€â”‚-â”‚
	  â”‚		       â”‚â—â”‚â—â”‚â—â”‚â—â”‚
	  â”‚		       â”‚â”€â”€â”€â”€â”€1â”€â”‚
	  â”‚	  X=d[0]       â”‚â—|â—|â—â”‚â—â”‚
   	(0,0)â”€â”€â”€â”€â”€â”€â”€>	     C[0]â”€â”€â”€â”€â”€C[1]

*/
typedef struct _wpsl_bloc{
  unsigned long sum;   // somme des aires des blocs "avant" (selon un DFS)
  struct _wpsl_coin *C[4]; // coins dÃ©finissant le bord du bloc (4 pointeurs)
  struct _wpsl_bloc *left,*right; // fils gauche et fils droit
} wpsl_bloc;


/* Structure "coin" pour wpsl():

  Un "coin" c est un point (x,y) de la grille qui est indicent Ã 
  quatres blocs B[0..3] numÃ©rotÃ©s ainsi:

          ^ Y=d[1]                 â”‚
	  â”‚	   	       B[3]â”‚B[2]
	  â”‚		     â”€â”€â”€â”€â”€â”€câ”€â”€â”€â”€â”€â”€
	  â”‚       X=d[0]       B[0]â”‚B[1]
	(0,0)â”€â”€â”€â”€â”€â”€â”€>              â”‚

  Comme prÃ©cÃ©demment pour la structure de bloc, le bloc B[0] du coin c
  est celui de coordonnÃ©es minimum (en bas Ã  gauche) et les autres en
  tournant vers la droite ce qui permet de passer le l'un Ã  l'autre
  simplement (+1 mod 4). Si un bloc n'existe pas il vaut NULL. Cela
  arrive pour les points du bord de la grille pxq. Il est possible que
  deux blocs B[i] et B[i+1] soient identiques. Par exemple:

                            b4 â”‚bâ”€â”€â”€â”€
                           â”€â”€â”€â”€aâ”‚
	                    b7 â”‚â”‚ 
	                   â”€â”€â”€â”€câ”‚ b5  
	                    b6 â”‚â”‚
	                   â”€â”€â”€â”€ddâ”€â”€â”€â”€

  Les 4 coins de c sont:          On a aussi:
  - le bloc B[0] de c est b6      - le bloc C[2] de d est b5
  - le bloc B[1] de c est b5      - le coin C[2] de b6 est c
  - le bloc B[2] de c est b5      - le coin C[1] de b7 est c
  - le bloc B[3] de c est b7      - le bloc C[3] de d est b6
*/
typedef struct _wpsl_coin{
  int d[2];        // coordonnÃ©es du coin (=point de la grille)
  wpsl_bloc *B[4]; // blocs incidents (4 pointeurs)
} wpsl_coin;


int wpsl(query* const Q)
/*
  Utilise Q->rep[u][0..3] pour la reprÃ©sentation implicite: seulement
  deux pÃ¨res pour le cas standard, et 4 pour le dual.

  Utilise un paramÃ¨tre cachÃ© Q->param[3] permettant de gÃ©rer toutes
  les variantes de wpsl (8 au total) grÃ¢ce aux 3 derniers bits.

    bit-0: dual ou pas
    bit-1: uniforme ou pas
    bit-2: dissection ou pas

  Notes:

  Le premier point a pas mal de chance de tomber sur un bord.  Par
  exemple, ./gengraph wpsl 1 9 9 -dot scale 0.1 -xy grid 9 -visu
  permet de voir qu'il y a 49 possibilitÃ©s dont 24 points sur le bord,
  soit 50% environ. De maniÃ¨re gÃ©nÃ©rale, il y a pas loin de 50% de
  chance de tomber Ã  une distance < (1/8)*p du bord pour une grille
  carrÃ©e p x p.

  Pour le dual, il est possible d'obtenir un graphe 4-rÃ©gulier (donc
  4-dÃ©gÃ©nÃ©rÃ© qui n'est pas 3-dÃ©gÃ©nÃ©rÃ©) comme ci-dessous. Mais je ne
  sais pas si 5-dÃ©gÃ©nÃ©rÃ© peut arriver. Il faudrait voir si l'icosaÃ¨dre
  5-rÃ©gulier peut Ãªtre rÃ©alisÃ©. En tout cas K_4 ne peut pas Ãªtre
  gÃ©nÃ©rÃ©, mais une subdivision si.

                         +â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€+
			 â”‚           â”‚
                         +â”€â”€+â”€â”€+â”€â”€+â”€â”€+
			 â”‚  â”‚  â”‚  â”‚  â”‚
			 â”‚  +â”€â”€+â”€â”€+  â”‚
			 â”‚  â”‚  â”‚  â”‚  â”‚
			 +â”€â”€+â”€â”€+â”€â”€+â”€â”€+
			 â”‚           â”‚
			 +â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€+

  Principe de la construction:

  On utilise les structures de "blocs" et de "coins" (voir au
  dessus). Un bloc est dÃ©finis par ses 4 coins (points de la grille),
  et un coin est un point de la grille ainsi que 4 blocs incidents Ã 
  ce point. Le cas gÃ©nÃ©ral et le cas de la dissection (diss=1) sont
  traitÃ©s similairement. On construit, Ã  l'initialisation, un tableau
  de blocs et la liste de coins associÃ©s. Et une fois construit, on en
  dÃ©duit une reprÃ©sentation implicite du graphe qui est soit une
  2-orientation ou 4-orientation pour le dual.

  Les blocs sont organisÃ©s selon un arbre binaire pour accÃ©lÃ©rer la
  crÃ©ation de nouveaux blocs qui revient Ã  faire une recherche
  (binaire) dans cet arbre. L'arbre, dont les noeuds sont les blocs,
  est implÃ©mentÃ© par un simple tableau de blocs car on connait leurs
  nombres: c'est au maximum 3n+1 pour le cas gÃ©nÃ©ral et n+1 si
  diss=1. Le premier bloc du tableau doit toujours Ãªtre la racine de
  l'arbre. L'aire d'un bloc est le nombre de cases qu'il contient, ce
  qui est diffÃ©rent du nombre de points de la grille qu'il contient.

  La remarque importante est que lorsqu'on subdivise un bloc en
  sous-blocs, la somme des aires des sous-blocs reste Ã©gale Ã  celle du
  bloc. (Attention! ici l'aire d'un bloc doit Ãªtre dÃ©fini comme le
  nombre de cases qu'il contient et non pas le nombre de points de la
  grille.) Donc en remplaÃ§ant l'ancien bloc par les nouveaux, on ne
  change pas la somme totale des aires des noeuds "avant" et des
  noeuds "aprÃ¨s", "avant" et "aprÃ¨s" se rÃ©fÃ©rant Ã  un DFS de l'arbre
  (cf. le champs "sum"). En choisissant un entier alÃ©atoire
  râˆˆ[0,dx*dy[, il suffit d'effectuer une recherche binaire dans
  l'arbre pour trouver le bloc B tq râˆˆ[B.sum, B.sum+aire(B)[. La
  propriÃ©tÃ© de somme constante n'est pas vraie si, au lieu du nombre
  de cases (=l'aire), on considÃ¨re le nombre de points internes de la
  grille qu'il contient. La subdivision ne garde pas constant ce
  nombre, et donc ne permet pas de faire une recherche efficace.

  Le coÃ»t de subdivision est proportionnel Ã  la hauteur de l'arbre,
  soit O(logn) en moyenne ce qui correspond aussi au nombre de
  subvisisions d'un mÃªme bloc. En effet, lorsqu'on remplace un bloc
  par ses 4 sous-blocs la hauteur de l'arbre n'augmente que d'une
  unitÃ©, et la subdivision a autant de chance de se faire Ã  droit qu'Ã 
  gauche de la racine (en fait c'est pas exactement 1/2, mais une
  probabilitÃ© constante >0).

  L'aire des blocs peut Ãªtre reprÃ©sentÃ©e par un intervalle dont la
  longueur est proportionnelle Ã  l'aire. Ã€ chaque moment, la somme des
  longueurs des blocs de l'arbre fait S=dx*dy. Dans l'exemple
  ci-dessous, les noeuds L,B,R forment une partition de
  [0,S[. Lorsqu'on subdivise B en sous-blocs b0,b1,b2,b3, on remplace
  le noeud B par le sous-arbre b0,b1,b2,b3 de sorte Ã  ne pas changer
  la somme des aires des fils Ã  gauche, b0 prenant la place de
  B. L'ordre relatif des sous-blocs n'a pas vraiment d'importance, et
  on pourrait optimiser l'ordre de faÃ§on Ã  Ã©quilibrer l'arbre.

  Exemple de subdivision d'un bloc B en sous-blocs b0,b1,b2,b3.  Les
  fils avant et aprÃ¨s selon le DFS se voient en projettant leurs
  intervalles sur horizontale.
 
                     /                               /
              [â”€â”€â”€â”€â”€Bâ”€â”€â”€â”€â”€[                        [b0[
               /         \                         /  \__ 
              /           \                     [b1[     [b2[
             /             \                    /        /  \
          [â”€Lâ”€[           [â”€Râ”€[             [â”€Lâ”€[     [b3[  [â”€Râ”€[

      [â”€â”€â”€[â”€â”€â”€[â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€[â”€â”€â”€[â”€â”€â”€[     [â”€â”€â”€[â”€â”€â”€[â”€â”€[â”€â”€[â”€â”€[â”€â”€[â”€â”€â”€[â”€â”€â”€[
       0                         S       0                         S

  Chaque bloc maintient ces 4 coins le dÃ©finissant (points de la
  grille du bord du bloc), notÃ©s C[0..3] en tournant autour de lui. En
  parallÃ¨le on maintenant un tableau Z de tous les coins rencontrÃ©s.
  Chaque coin maintient aussi ses 4 blocs incidents notÃ©s B[0..3] en
  tournant de lui. En jouant sur les coins d'un bloc et les blocs d'un
  coin, on peut dÃ©terminer tous les points qui sont sur le bord d'un
  bloc, et aussi les blocs voisins d'un bloc donnÃ©s.

                            +â”€â”€â”€â”€++â”€â”€â”€â”€â”€â”€+
			    â”‚ b7 â”‚â”‚   â”‚  â”‚
			    â”‚    â”‚â”‚â”€â”€â”€Xâ”€â”€â”‚
			    +â”€â”€â”€â”€+â”‚   â”‚  â”‚
			    â”‚ b6 â”‚â”‚ b5â”‚  â”‚
			    +â”€â”€â”€â”€++â”€â”€â”€â”€â”€â”€+

  Ci-dessus, trois blocs b5,b6,b7 avant la subdivision de b5 au point
  interne X de la grille. Les coins des blocs sont rÃ©unis dans le
  tableau Z, ce qui Ã©vite les duplicats. Par exemple, le coin C[3] de
  b5 est le mÃªme que le coin C[2] de b7. Pour savoir si les coins
  crÃ©es par le point X sont nouveaux, il faut vÃ©rifier les coins des
  blocs adjacents Ã  b5. Cela peut se faire grÃ¢ce au tableau Z. Le coÃ»t
  de cette vÃ©rification est en O(deg(b5)), le degrÃ© de b5 dans le
  dual. Sur tous les blocs, cela se somme donc en O(n*h) oÃ¹ h=O(logn)
  est la hauteur finale de l'arbre.

  AmÃ©lioration: subdiviser les blocs jusqu'Ã  ce qu'on ne puisse
  plus, mais en temps O(nlogn), comme expliquÃ© ci-dessous.

  Au lieu de rÃ©pÃ©ter n fois la subdivision, on pourrait boucler
  jusqu'Ã  crÃ©er n subdivisions, bien sÃ»r Ã  condition que le nombre de
  points internes du dÃ©coupage courant le permette. Mais trouver la
  derniÃ¨re subdivision pourrait prendre un temps O(n^2) par un
  algorithme Ã  rejet simple (pour trouver le dernier point interne!),
  soit O(n^3) au total. On peut faire en temps O(nlogn) en moyenne
  avec un nombre d'Ã©tapes en O(n) en moyenne, mÃªme si on souhaite
  subdiviser complÃ¨tement la grille pxq.

  Chaque bloc B maintient un poids B->weight qui au dÃ©part vaut le
  nombre de points internes Ã  B ainsi que B->sum qui est la somme des
  B->weight des fils avant B dans T. Donc B->sum ne correspond plus Ã 
  une surface, mais on va s'arranger pour garder ces valeurs
  constantes, mÃªme si c'est incorrect. On mesure le taux d'incohÃ©rence
  par deux nombres: SP(T) qui est la somme des points internes dans T
  et SW(T) qui est la plus grande valeur B->sum+B->weight.

  Au dÃ©part, lorsqu'il n'y a qu'un bloc, SP(T)=SW(T).  Lorsqu'on
  subdivise un bloc B en b[0]..b[3], on met Ã  jour facilement SP(T),
  qui devient < SW(T), mais on ne touche pas aux champs ->sum sauf
  pour b[0]..b[3]. Donc SW(T) reste constant. On rÃ©partit B->weight en
  b[i]->weight de sorte Ã  garder la somme des b[i]->weight constante Ã 
  B->weight et Ã  respecter le proportion du nombre de points internes
  dans chaque b[i]. On met aussi Ã  jour b[i]->sum en accord avec les
  b[i]->weight. La probabilitÃ© de choisir un bloc B est donc
  proportionnel Ã  B->weight car on tire un nombre dans [0,SW(T)[. Si
  SW(T)/2 < SP(T) <= SW(T) alors la probabilitÃ© de tomber sur un bloc
  qui ne peut pas Ãªtre subdivisÃ© est 1 - (SW(T)/SP(T)) <= 1/2.

  Si SP(T) <= SW(T), alors on met Ã  jour B->weight et B->sum de tous
  les blocs B de T ce qui va coÃ»ter O(|T|) (comment?). SW(T) et SP(T)
  sont Ã©galement mise Ã  jour, et bien sÃ»r on a alors SW(T)=SP(T). Les
  blocs B non-subdivisables vont recevoir B->weight=0 si bien qu'ils
  ne pourront Ãªtre choisi dans la recherche dans T. Bien sÃ»r il faut
  stopper les subdivisions si jamais SP(T)=0. Le temps de recherche
  est toujours en O(logn) en moyenne. (Certes il y a des blocs de
  poids nuls, mais un bloc ne peut Ãªtre subdivisÃ© plus de logn fois en
  moyenne.) Ainsi, la probabilitÃ© de choisir un bloc de T qui ne peut
  Ãªtre subdivisÃ© redescend alors Ã  0, et va augmenter progressivement
  jusqu'Ã  1/2. Le nombre d'Ã©tapes sera donc O(logn) mÃªme si SW(T) est
  trÃ¨s petit.

  On remarque qu'il ne peut avoir que O(log(pq)) mises Ã  jour car, Ã 
  chaque mise Ã  jour de T, SP(T) est divisÃ© par deux, et au dÃ©part
  SP(T)<pq. En fait, si pq>n^2, alors on ne fera pas de mise Ã  jour
  puisqu'il faudrait choisir deux valeurs (sur les n essais
  alÃ©atoires) de [0,pq[ Ã  distance < n < sqrt(pq) pour obtenir un bloc
  non-subdivisable. Donc pq <= n^2 et log(pq) = O(logn). En fait, si
  l'une des dimensions est trÃ¨s petite (et l'autre >n^2), on va
  effectuer log(min{p,q}) mise Ã  jour. Dans tous les cas, on a log(pq)
  = O(logn).
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q)||free_pos(Q);

  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:;
    int const dx=Q->param[2]-1; // = largeur de la grille (en #cases)
    int const dy=Q->param[1]-1; // = hauteur de la grille (en #cases)
    int n=Q->param[0]; // = nombre d'Ã©tapes
    if((dx<1)||(dy<1)||(n<0)) Erreur(6); // paramÃ¨tre incorrect
    int const dual=((Q->param[3]&1)>0); // dual (=1) ou pas
    int const unif=((Q->param[3]&2)>0); // uniforme (=1) ou pas (surface)
    int const diss=((Q->param[3]&4)>0); // dissection (=1) ou pas (croix)
    
    NALLOC(wpsl_bloc,T,(3-2*diss)*n+1); // tableau d'au plus 3n+1 ou n+1 blocs
    NALLOC(wpsl_coin,Z,(5-3*diss)*n+4); // tableau d'au plus 5n+4 ou 2n+4 coins
    int i,j,k,rx,ry,bx,by,u;
    long r;

    // initialisation du 1er bloc et de ses 4 coins
    //
    //          C[3]â”€â”€â”€â”€C[2]
    //            â”‚      â”‚
    //            â”‚ T[0] â”‚
    //            â”‚      â”‚
    //          C[0]â”€â”€â”€â”€C[1]
    //
    int tn=1; // tn = nombre de blocs dans T
    int zn=4; // zn = nombre de coins dans Z
    T->sum=0; // somme des aires des blocs avant
    T->left=T->right=NULL; // au dÃ©part, ni fils gauche et ni fils droit
    for(i=0;i<4;i++) T->C[i]=Z+i; // NB: "T->" = "T[0]."      
    for(i=0;i<4;i++) for(j=0;j<4;j++) Z[i].B[j]=NULL;
    Z[0].d[0]=0,  Z[0].d[1]=0,  Z[0].B[2]=T;
    Z[1].d[0]=dx, Z[1].d[1]=0,  Z[1].B[3]=T;
    Z[2].d[0]=dx, Z[2].d[1]=dy, Z[2].B[0]=T;
    Z[3].d[0]=0,  Z[3].d[1]=dy, Z[3].B[1]=T;

    long const S=(long)dx*dy; // S = aire totale de la grille, (long) est important
    wpsl_bloc *B;
    wpsl_coin *C;

/* largeur (i=0) ou hauteur (i=1) du bloc B */
#define DELTA(i,B) ((B)->C[2]->d[i] - (B)->C[0]->d[i])
#define SURFACE(B) ((long)DELTA(0,B)*(long)DELTA(1,B))

    ////////////////////
    // rÃ©pÃ©ter n fois //
    ////////////////////

    for(;n--;){

      DEBUG(PRINT(n);wpsl_DEBUG;PAUSE;);

      // choisit un bloc B dans T
      if(unif) B=T+randomu(tn); // tirage uniforme dans T
      else{ // tirage suivant la surface des blocs dans T
	r=randomu(S); // tire une case alÃ©atoire de la grille (sur 62 bits)
	for(B=T;;){ // cherche dans T le bloc B contenant r
	  if(r<B->sum){ B=B->left; continue; }
	  if(r<B->sum+SURFACE(B)) break; // B est trouvÃ© !
	  B=B->right;
	}
      }
      
      // tire un point (rx,ry) de la grille interne au bloc B. Si
      // diss=1, (rx,ry) est sur le bord 2 ou 3. Ce point (rx,ry) sera
      // le coin c[4].
      bx=DELTA(0,B); // bx = largeur de B
      by=DELTA(1,B); // by = hauteur de B
      if(bx<2&&by<2) continue; // aucun sous-bloc crÃ©e dans ce cas

      C=B->C[0]; // C = coin C[0] de B 
      rx=C->d[0], ry=C->d[1]; // par dÃ©faut (rx,ry)=coin C[0] de B
      if(diss){ // coupe en deux B
	k=(by<2)||((bx>1)&&(unif?RANDbit:(int)randomu(bx+by)<bx)); // k = sens de la coupe
	if(k) rx += 1+randomu(bx-1); // k=1 -> coupe X (vertical)
	else  ry += 1+randomu(by-1); // k=0 -> coupe Y (horizontal)
      }else{ // coupe en quatre B
	if(bx<2||by<2) continue; // aucun sous-bloc crÃ©e dans ce cas
	rx += 1+randomu(bx-1); // rxâˆˆ[C->d[0]+1,C->d[0]+bx[
	ry += 1+randomu(by-1); // ryâˆˆ[C->d[1]+1,C->d[1]+by[
      }

      /* subdivise B en b[0..3] (seulement b[0..1] si diss=1):
	 b[0]=B, b[1]=&T[tn], b[2]=&T[tn+1], b[3]=&T[tn+2]
	 Le coin B->C[0] a pour coordonnÃ©es (x,y)
      
                   +â”€â”€â”€â”€â”€c[1]â”€â”€â”€â”€â”€â”€+  y+by
                   â”‚  b[3] â”‚  b[2] â”‚
            B =  c[2]â”€â”€â”€â”€c[4]â”€â”€â”€â”€c[0] ry
                   â”‚  b[0] â”‚  b[1] â”‚
                   +â”€â”€â”€â”€â”€c[3]â”€â”€â”€â”€â”€â”€+  y
                   x      rx     x+bx
      
	 Les bords d'un bloc sont numÃ©rotÃ©s iâˆˆ{0,1,2,3} comme ceci:

	      1        0: bord rencontrÃ© selon X croissant
	    â”¼â”€â”€â”€â”¼      1: bord rencontrÃ© selon Y croissant
	   2â”‚ B â”‚0     2: bord rencontrÃ© selon X dÃ©croissant
	    â”¼â”€â”€â”€â”¼      3: bord rencontrÃ© selon Y dÃ©croissant
              3        NB: bord opposÃ© de i est (i XOR 2)=(i^2)
      */

      // on ajoute les 2 ou 4 blocs Ã  T en Ã©crasant B
      wpsl_bloc *b[]={ B, T+tn, T+tn+1, T+tn+2 };
      wpsl_bloc const BB=*B; // = copie en dure de *B
      if(diss){ // il y a des blocs Ã©gaux dans ce cas
	// coupe X (verticale):   b[0] = (b[0]|b[3]) et b[1] = (b[1]|b[2])
	// coupe Y (horizontale): b[0] = (b[0]|b[1]) et b[1] = (b[3]|b[2])
	b[2]=b[1]; b[3]=b[!k]; if(!k) b[1]=b[0];
      }
      *b[1]=*b[2]=*b[3]=BB; // fait des copies de B, NB: on a dÃ©jÃ  b[0]=B
      tn += 3-2*diss; // 1 ou 3 blocs de plus dans T

      // on crÃ©e 5 coins c[0..4], qui seront ou pas nouveaux, certains
      // de leurs blocs (2 sur 4) ne sont pas encore connus
      wpsl_coin c[]={
	{ .d[0]=C->d[0]+bx, .d[1]=ry, .B[0]=b[1], .B[1]=NULL, .B[2]=NULL, .B[3]=b[2] },
	{ .d[0]=rx, .d[1]=C->d[1]+by, .B[0]=b[3], .B[1]=b[2], .B[2]=NULL, .B[3]=NULL },
	{ .d[0]=C->d[0],    .d[1]=ry, .B[0]=NULL, .B[1]=b[0], .B[2]=b[3], .B[3]=NULL },
	{ .d[0]=rx,    .d[1]=C->d[1], .B[0]=NULL, .B[1]=NULL, .B[2]=b[1], .B[3]=b[0] },
	{ .d[0]=rx,         .d[1]=ry, .B[0]=b[0], .B[1]=b[1], .B[2]=b[2], .B[3]=b[3] }
      };

      // On initialise les blocs des 4 coins de B Ã  b[0..3] ce qui ne
      // va pas Ãªtre fait lors du parcours des 4 bords de B.
      for(i=0;i<4;i++) BB.C[i]->B[(i+2)&3]=b[i];

      if(!diss){
	// On ajoute le coin c[4] Ã  Z car il est nouveau, et aussi aux
	// 4 blocs b[0..3]. Les 4 coins diagonalement opposÃ©s Ã  c[4]
	// pour les blocs b[0..3] correspondent aux 4 coins de B. Ils
	// sont donc dÃ©jÃ  corrects car les b[0..3] sont initialisÃ©s Ã 
	// B. Il manque deux coins pour les blocs b[0..3] qui seront
	// dÃ©terminÃ©s lorsqu'on saura si c[i] est un coin nouveau ou
	// s'il est dÃ©jÃ  dans Z.
	for(i=0;i<4;i++) b[i]->C[(i+2)&3]=Z+zn; // NB: Z+zn correspond Ã  c[4]
	Z[zn++]=c[4]; // un coin de plus dans Z
      }

      // Parcoure chacun des 4 bords (i=0..3) de B selon les
      // coordonnÃ©es croissantes pour dÃ©terminer si c[i] existe ou pas
      // dans Z. Pour chaque coin rencontrÃ© on met Ã  jour ses 2 blocs
      // incidents cotÃ© B avec un certain bloc b[i']. Lorsqu'on
      // atteint c[i], on met Ã  jour aussi les 2 blocs externes Ã  B si
      // c[i] est nouveau. AprÃ¨s c[i], on continue la mise Ã  jour mais
      // avec un autre bloc b[i"]. On n'a pas Ã  mettre Ã  jour les
      // coins min et max de B qui ont dÃ©jÃ  Ã©tÃ© correctement
      // initialisÃ©s.

      int c1,c2,c3,b1,b2,i0,i1,i2,i3;
      // C[c1] = coin infÃ©rieur de B sur bord i
      // C[c2] = coin supÃ©rieur de B sur bord i
      // C[c3] = prochain coin sur bloc B[i2] sur bord i        j
      // B[i0] = bloc cotÃ© B sur bord i avant C                 ^          B
      // B[i1] = bloc adjacent Ã  B sur bord i avant C        i3 â”‚ i2    i0 â”‚ i3
      // B[i2] = bloc adjacent Ã  B sur bord i aprÃ¨s C       B â”€â”€Câ”€â”€      â”€â”€Câ”€â”€> j
      // B[i3] = bloc cotÃ© B sur bord i aprÃ¨s C              i0 â”‚ i1    i1 â”‚ i2
      // b[b1] = bloc remplaÃ§ant B pour C si avant c[i]
      // b[b2] = bloc remplaÃ§ant B pour C si aprÃ¨s c[i]
      // propriÃ©tÃ©: c3=(c1+2)&3, i1=b1=c1, i2=b2=c2, i3=c3 et i0=(i2+2)&3

      for(i=j=0;i<4;i++){ // pour chaque bord i=0..3
	j=1-j; // j=1,0,1,0: dimension du bord i: X (=0) ou Y (=1)
	c1=(i+1)&3, c2=(i+2)&3; if(i%3) SWAP(c1,c2);
	c3=(c1+2)&3, i1=b1=c1, i2=b2=c2, i3=c3, i0=(i2+2)&3;

	// avant c[i]
	C=BB.C[c1];  // C = coin courant = coin infÃ©rieur de B du bord i
	r=c[i].d[j]; // r = coordonnÃ©es de c[i]
	while(C->B[i2] && C->B[i2]->C[c3]->d[j]<=r){ // hypothÃ¨se: le coin C est Ã  jour
	  C=C->B[i2]->C[c3]; // prochain coin
	  C->B[i0]=C->B[i3]=b[b1]; // met Ã  jour les blocs i0 et i3 de C
	}
	  
	// c[i] est un coin existant ou pas ? c[i] est nouveau si le
	// prochain coin de C n'existe pas ou bien C<>c[i]. Dans les
	// autres cas c'est que C=c[i], et donc a Ã©tÃ© traitÃ©. Si
	// diss=1, et si le bord est parallÃ¨le Ã  la coupe (k=j), alors
	// c[i] est le coin infÃ©rieur de B. Il ne faut pas traiter la
	// partie avant c[i] qui n'existe pas.
	if(diss && k==j) goto wpsl_bpc; // bord parallÃ¨le Ã  la coupe

	if((C->B[i2]==NULL)||(C->d[j]!=r)){ // nouveau coin
	  c[i].B[i0]=b[b1], c[i].B[i3]=b[b2]; // met Ã  jour les 2 blocs de c[i] interne Ã  B
	  c[i].B[i1]=c[i].B[i2]=C->B[i2]; // met Ã  jour les 2 blocs inconnus de c[i] hors de B
	  C=Z+zn; // C pointe sur le nouveau coin 
	  Z[zn++]=c[i]; // insÃ¨re c[i] dans Z
	}

	// ici C = coin c[i] de B (nouveau ou pas)
	C->B[i3]=b[b2]; // corrige le bloc i3 de C (=c[i])
	b[b1]->C[b2]=b[b2]->C[b1]=C; // les 2 coins manquants pour les blocs b[b1] et b[b2]
		
      wpsl_bpc:
	// aprÃ¨s c[i] (hypothÃ¨se: le coin courant C est Ã  jour)
	r=BB.C[c2]->d[j]; // r = coordonnÃ©es du coin max de B
	while(C->B[i2] && C->B[i2]->C[c3]->d[j]<r){ // ne pas traiter le coin max
	  C=C->B[i2]->C[c3]; // prochain coin
	  C->B[i0]=C->B[i3]=b[b2]; // met Ã  jour les blocs i0 et i3 de C
	}
	
      }// fin du for(i=...)
	
      // changements de l'arbre T

      if(diss){
	// NB: Le bloc b[0] est rangÃ© Ã  l'adresse B dans T. Mais
	// attention! l'autre bloc rangÃ© en T+tn, n'est pas b[1] mais
	// b[2]. Car si k=0, on a fait b[1]=b[0].
	if(RANDbit){
	  // Pour Ã©quilibrer l'arbre T, on met alÃ©atoirement b[2] Ã 
	  // droite ou Ã  gauche de b[0].
	  b[0]->left=b[2];
	  b[2]->right=NULL;
	  if(!unif) b[0]->sum=b[2]->sum+SURFACE(b[2]); // met Ã  jour la somme si besoin
	}else{
	  b[0]->right=b[2];
	  b[2]->left=NULL;
	  if(!unif) b[2]->sum=b[0]->sum+SURFACE(b[0]); // met Ã  jour la somme si besoin
	}
      }else{
	// Pour Ã©quilibrer un peu plus l'arbre T, on pourrait
	// rÃ©ordonner les fils. Mais attention! il ne faut pas
	// permuter les blocs dans T, mais simplement les fils de ces
	// blocs. On ne peut pas changer b[0] car on ne connaÃ®t son
	// pÃ¨re Il faudrait le stocker lors de la recherche binaire.
	b[0]->left=b[1];
	b[0]->right=b[2];
	b[1]->right=NULL;
	b[2]->left=b[3];
	b[3]->left=b[3]->right=NULL;
	if(!unif){ // met Ã  jour les sommes si besoin
	  // NB: assert(b[1]->sum==BB.sum);
	  b[0]->sum=b[1]->sum+SURFACE(b[1]);
	  b[3]->sum=b[0]->sum+SURFACE(b[0]);
	  b[2]->sum=b[3]->sum+SURFACE(b[3]);
	}
      }
      
    } // fin du rÃ©pÃ©ter n fois ...

    DEBUG(printf("Fin:\n");wpsl_DEBUG;PAUSE;);

    ////////////////////////////////////////////////
    // ici les tableaux T et Z ont Ã©tÃ© construits //
    ////////////////////////////////////////////////

    /* calcul de Q->rep[..] */
    
    if(dual){
      
      /* 
	 Dans ce cas on calcule une 4-orientation comme ceci. Les
	 adjacences d'un bloc B donnÃ© sont dÃ©terminÃ©es par les
	 intersections de ses quatres bords reprÃ©sentant un intervalle
	 et numÃ©rotÃ©s iâˆˆ{0,1,2,3} comme prÃ©cÃ©demment.

	 Pour que le bloc B soit adjacent au bloc B' par le bord i de
	 B et par le bord i' de B', il faut que i et i' s'intersectent
	 selon un intervalle non nul. Il n'y aura que deux bords Ã 
	 considÃ©rer, car si B est adjacent Ã  B' c'est que i'=i^2. Ã€ un
	 renomage prÃ¨s, on pourra supposer que B est le bloc avec
	 i<i', ce qui ne fait que deux cas: i=0 ou i=1.

	 Supposons que B est adjacent Ã  B' (intervalle non nul) par
	 les bords i et i'. On pose i=[i1,i2] et i'=[i1',i2'] les
	 intervalles respectifs. On oriente B->B' dans 3 des 5 cas
	 possibles (cf. figure ci-dessous), et bien sÃ»r B<-B'
	 sinon. Cas 1 ou 2: i est inclu dans i', soit i1'<=i1<i2<=i2'.
	 Cas 3: i chevauche infÃ©rieurement i', la plus grande
	 extrÃ©mitÃ© de i est au-dessous de la plus grande extrÃ©mitÃ© de
	 i', soit i1<i1'<i2<i2'. Dans les autres cas, c'est B'->B.

	 Adjacence pour le cotÃ© 0 du bloc B:
	 
	 1: B -> B'   2: B -> B'   3: B -> B'   4: B <- B'   5: B <- B'
	 +â”€â”€â”€++â”€â”€+    +â”€â”€â”€++â”€â”€+         +â”€â”€+    +â”€â”€â”€+        +â”€â”€â”€++â”€â”€+
	 â”‚ B â”‚â”‚B'â”‚    â”‚ B â”‚â”‚B'â”‚    +â”€â”€â”€+â”‚B'â”‚    â”‚ B â”‚+â”€â”€+    â”‚ B â”‚â”‚B'â”‚
	 +â”€â”€â”€+â”‚  â”‚    +â”€â”€â”€++â”€â”€+    â”‚ B â”‚+â”€â”€+    +â”€â”€â”€+â”‚B'â”‚    â”‚   â”‚+â”€â”€+
	      +â”€â”€+                 +â”€â”€â”€+             +â”€â”€+    +â”€â”€â”€+
				 
         On remarque qu'il ne peut y avoir pour le bloc B qu'un seul
         arc sortant par son bord i. En effet, soit B'_1...B'_k tous
         les blocs adjacents Ã  B par le bord i. Alors un seul des cas
         1,2,3 ne peut s'appliquer et qu'Ã  un seul bloc B'_j car dans
         chaque cas l'extrÃ©mitÃ© la plus grande de i est incluse dans
         le bord i^2 de B'_j. Et cela ne peut se produire qu'une seule
         fois.

	 Au total, le bloc B ne peut avoir que 4 arcs sortants, un par
         bord. Cela peut arriver avec une situation comme celle-ci:

 	                       +â”€â”€â”€+â”€â”€â”€+
			       â”‚   â”‚   â”‚
			       â”‚   +â”€â”€â”€+â”€â”€â”€+
			       â”‚   â”‚ B â”‚   â”‚
			       +â”€â”€â”€+â”€â”€â”€+â”€â”€â”€+
			           â”‚       â”‚
			           +â”€â”€â”€â”€â”€â”€â”€+
       */
      
      SET_n(tn); // = nombre de sommets dans ce cas
      ALLOCZ(Q->xpos,Q->n,(double)(T[_i].C[0]->d[0]+T[_i].C[2]->d[0]));
      ALLOCZ(Q->ypos,Q->n,(double)(T[_i].C[0]->d[1]+T[_i].C[2]->d[1]));
      Q->k=4; // pour adjacency_rep()
      ALLOC2Z(Q->rep,Q->n,4,-1); // reprÃ©sentation implicite
      wpsl_bloc *Bu,*Bv; // Bu=B et Bv=B'
      int i1,i2,i3,i4;

#define ADJ ((i1<=i3 && i3<i2)||(i3<=i1 && i1<i4)) /* i3âˆˆ[i1,i2[ ou i1âˆˆ[i3,i4[ */
#define CAS12 (i3<=i1 && i2<=i4)        /* [i1,i2] inclu dans [i3,i4] */
#define CAS3  (i1<i3 && i3<i2 && i2<i4) /* [i1,i2] chevauche infÃ©rieurement [i3,i4] */
	
      for(u=0,Bu=T;u<tn;u++,Bu++){ // pour chaque sommet u, i.e. chaque bloc Bu=T[u] de T
	for(i=j=0;i<2;i++){ // d'abord pour le bord i=0, puis le bord i=1
	  j=1-j; // dimension 1 puis 0
	  Bv=Bu->C[1+2*i]->B[2]; // Bv = 1er bloc adjacent Ã  Bu
	  while(Bv){ // on parcoure les blocs Bv touchant Bu par le bord i de Bu
	    i1=Bu->C[1+2*i]->d[j], i2=Bu->C[2]->d[j]; // [i1,i2] = bord i de Bu 
	    i3=Bv->C[0]->d[j], i4=Bv->C[1+2*j]->d[j]; // [i3,i4] = bord i de Bv 
	    if(!ADJ) break; // si Bv n'est plus adjacent Ã  Bu
	    if(CAS12||CAS3) Q->rep[u][i]=(int)(Bv-T); // Bu->Bv par le bord i de Bu
	    else Q->rep[(int)(Bv-T)][i^2]=u; // Bv->Bu par le bord opposÃ© de Bv
	    Bv=Bv->C[1+2*j]->B[2]; // Bv = prochain bloc de Bv
	  }
	}
      }
      
#undef ADJ
#undef CAS12
#undef CAS3
      
    }else{

      /* 	       
	 Chaque sommet (=coin) a au plus 4 voisins, deux alignÃ©s selon
	 X et deux alignÃ©s selon Y. Par symÃ©trie, on ne calcule que
	 deux voisins pour chaque sommet: le voisin selon les X
	 croissant (celui entre les blocs B1 et B2) et les Y croissant
	 (celui entre les blocs B2 et B3). Cela donne une
	 2-orientation que l'on stocke dans Q->rep[0..1].
      */
       
      SET_n(zn); // = nombre de sommets dans ce cas
      ALLOCZ(Q->xpos,Q->n,(double)Z[_i].d[0]);
      ALLOCZ(Q->ypos,Q->n,(double)Z[_i].d[1]);
      Q->k=2; // pour adjacency_rep()
      ALLOC2Z(Q->rep,Q->n,2,-1); // reprÃ©sentation implicite

      for(u=0,C=Z;u<zn;u++,C++){ // pour tous sommets u = les coins C = Z[u] de Z
	/*
	  Pour dÃ©terminer le voisin du sommet u selon les Y croissant
	  par exemple (=bord 0), il faut considÃ©rer les deux voisins
	  possibles: le coin C[3] du bloc B[2] de u et le coin C[2] du
	  bloc B[3] de u (cf. figure ci-dessous). Il faut que ces
	  coins aient la mÃªme abscisse que u et ensuite la plus petite
	  ordonnÃ©es possible si les deux sont possibles (cas 1).
	  Attention! Le voisin peut ne pas exister si les deux coins
	  n'ont pas la mÃªme abscisse que u (cas 2) ou si le bloc
	  n'existe tout simplement pas car u est sur un bord de la
	  grille.
	
              Cas 1        â”‚C[3]â”€â”€        Cas 2
                     â”€â”€C[2]â”‚â”‚                C[3]â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€C[2]
                     B[3] â”‚â”‚â”‚ B[2              â”‚ B[3] = B[2] â”‚
                     â”€â”€â”€â”€â”€+â”‚+â”€â”€â”€â”€â”€             +â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€+
                    â”€â”€â”€â”€â”€â”€â”€uâ”€â”€â”€â”€â”€â”€â”€           â”€â”€â”€â”€â”€â”€â”€â”€uâ”€â”€â”€â”€â”€â”€â”€â”€
                           â”‚                          â”‚
	 
	  Pour le voisin de u selon les X croissant c'est pareil mais
	  avec les blocs B[1] et B[2], et en inversant les X avec les
	  Y.
	*/

	for(i=0;i<2;i++){ // pour i=0 (X croissant), puis i=1 (Y croissant)
	  wpsl_coin *p1=C->B[i+1]? C->B[i+1]->C[i+2] : NULL;
	  wpsl_coin *p2=C->B[i+2]? C->B[i+2]->C[i+1] : p1; // NB: si p2=NULL, alors p1 aussi
	  if(p1 && (p1->d[1-i]!=C->d[1-i])) p1=p2;
	  if(p2 && (p2->d[1-i]==C->d[1-i])){
	    if(p1 && (p1->d[i]<p2->d[i])) p2=p1;
	    Q->rep[u][i]=(int)(p2-Z);
	  }
	}
      } // fin du for(u=...)
    } // fin du cas classique
    
    free(T);
    free(Z);
    XYtype=XY_USER; // coordonnÃ©es fixÃ©es par l'utilisateur
    return InitXY(Q); // pour les options -xy noise/scale ...

  }

  return 1;
}
#undef DELTA
#undef SURFACE


int ngon(query* const Q)
/*
  Utilise i<j.

  codage de x:
  - bit-0 (x&1): AC connectÃ© Ã  A
  - bit-1 (x&2): BC connectÃ© Ã  B
  - bit-2 (x&4): carrÃ©
  - bit-3 (x&8): asymÃ©trie
*/
{
  int const p=Q->param[0];
  int c=Q->param[1];
  int const x=Q->param[2];

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
 
  case QUERY_INIT:
    if(x==-1){
      if(p<3) RET_n(0);
      SET_n(p); // n=p
    }
    if(x==-2){
      if(p<1) RET_n(0);
      SET_n(3*p); // n=3p
    }
    if(x>=0){
      // il faut c dans [0,p/2] et p>=1
      if((c<0)||(c>p/2)||(p<1)) RET_n(0);
      if((x<0)||(x>15)) RET_n(0);
      if(x&4) SET_n(4*p); // n=4p
      else SET_n(3*p); // n=3p
    }

    // donne des coordonnÃ©es aux points
    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);
    double const t=M_2PI/Q->n;
    double a=0;
    for(int i=0;i<Q->n;i++){
      Q->xpos[i]=cos(a);
      Q->ypos[i]=sin(a);
      a += t;
    }
    return 0;

  case QUERY_ADJ:;
    int const n=Q->n;
    int u=Q->i;
    int v=Q->j;

    // teste le cycle
    if((v==u+1)||((v==n-1)&&(u==0))) RET_a(1);

    // triangulation Ã  trois fans
    if(x==-1){
      if((v==c)||(u==n-c)) RET_a(1);
      if((u==0)&&(v==n-c)) RET_a(1);
      if((u==0)&&(c<v)&&(v<n-c)) RET_a(1);
      RET_a(0); // u,v non adjacents
    }

    int const p2=2*p; // p2=2p
    int const p3=3*p; // p3=3p
    int s=0; // vrai si asymÃ©trique

    // diagonale du carrÃ© central
    if((u==0)&&(x&4)&&(v==p2)) RET_a(1);

    // on se ramÃ¨ne Ã  u<v dans l'intervalle [0,p]
    for(;;){
      if((u<p)&&(v<=p)) break;                               // intervalle [0,p]
      if((p<=u)&&(u<p2)&&(v<=p2)){ u-=p; v-=p; s=1; break; } // intervalle [p,2p]
      if((p2<=u)&&(u<p3)&&(v<=p3)){ u-=p2; v-=p2; break; }   // intervalle [2p,3p]
      if(p3<=u){ u-=p3; v-=p3; s=1; break; }                 // intervalle [3p,n]
      if(u==0){ // cas particulier avec u=0
	if((x&4)&&(x>=0)){
	  if(v>=p3){ u=v-p3; v=p; s=1; break; }              // intervalle [3p,n]
	}else if(v>=p2){ u=v-p2; v=p; break; }               // intervalle [2p,n]
      }
      RET_a(0); // u,v non adjacents
    }
    // ici u<v dans [0,p]

    // triangulation rÃ©cursive
    if(x==-2){
      int k=1;
      if(v==p){ // si v=p -> v = puissance de 2 >=v
	while(v>k) k<<=1;
	v=k;
      }
      k=1;
      int const w=v-u;
      while(w>k) k<<=1;
      if((w==k)&&(v%k==0)) RET_a(1);
      RET_a(0); // u,v non adjacents
    }

    if((u==0)&&(v==p)) RET_a(1); // grand triangle ou carrÃ©
    if((x&8)&&s) c=p-c; // asymÃ©trie

    int const A=x&1; // vrai si Ã©toile depuis A
    int const B=x&2; // vrai si Ã©toile depuis B

    if(v<=c){ // si u,v dans [A,C]
      if(A) RET_a(u==0); // Ã©toile depuis A
      RET_a(v==c);       // Ã©toile depuis C
    }
    if(u>=c){ // si u,v dans [C,B]
      if(B) RET_a(v==p); // Ã©toile depuis B
      RET_a(u==c);       // Ã©toile depuis C
    }
    RET_a(0);
  }

  return 1; // fin anormale
}


int ngon2(query* const Q)
/*
  Utilise i<j.

  paramÃ¨tres:
  Q->param[0] = s = nombre de segments
  Q->sparam = mot binaire de longueur 2(s-1), ex: "1100"
  Q->param[i] = nombre d'arÃªtes du segment entre les points C_{i-1} et C_i, i=1..s
    si >=0, Ã©toile dans le segment i depuis C_{i-1}
    sinon, Ã©toile dans le segment i depuis C_i

  Une demi-lune [A,B] dÃ©coupÃ©e en s segments [C_{i-1},C_i] avec A=C_0
  et B=C_s.
*/
{
  int const s=Q->param[0];
  int i,p;

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
 
  case QUERY_INIT:
    p=0;
    for(i=1;i<=s;i++) p += Q->param[i];
    SET_n(3*p);
    
    // donne des coordonnÃ©es aux points
    ALLOC(Q->xpos,Q->n);
    ALLOC(Q->ypos,Q->n);
    double const t=M_2PI/Q->n;
    double a=0;
    for(i=0;i<Q->n;i++){
      Q->xpos[i]=cos(a);
      Q->ypos[i]=sin(a);
      a += t;
    }
    return 0;

  case QUERY_ADJ:;
    int const n=Q->n;
    p=n/3;
    int u=Q->i;
    int v=Q->j;

    // teste le cycle
    if((v==u+1)||((v==n-1)&&(u==0))) RET_a(1);

    int const p2=2*p; // p2=2p

    // on se ramÃ¨ne Ã  u<v dans l'intervalle [0,p]
    for(;;){
      if((u<p)&&(v<=p)) break;                               // intervalle [0,p]
      if((p<=u)&&(u<p2)&&(v<=p2)){ u-=p; v-=p; break; }      // intervalle [p,2p]
      if((p2<=u)&&(u<n)&&(v<=n)){ u-=p2; v-=p2; break; }     // intervalle [2p,3n]
      RET_a(0); // u,v non adjacents
    }

    // ici u<v dans [0,p]

    RET_a(0);
  }

  return 1; // fin anormale
}


int behrend(query* const Q)
/*
  Utilise i<j.
  Adjacence en O(log(p)).
*/
{
  int const p=Q->param[0];
  int const k=Q->param[1];
  int a,b,w,s,i;

  switch(Q->code){

  case QUERY_END:
    free(Q->wrap),Q->wrap=NULL;
    return 0;
    
  case QUERY_INIT:
    if((p<2)||(k<2)) RET_n(0);
    SET_n(p*k);
    /*
      On construit un ensemble X d'Ã©lÃ©ments de [0,p[ Ã  partir de
      toutes les permutations P de [0,a[ interprÃ©tÃ©s comme entier en
      base b. Donc |X|=a! et |X|<p. Les entiers a et b doivent
      vÃ©rifier:
      
      1) a<b/k
      2) max{X}<p/(k-1) oÃ¹ max{X}=sum_{i=0}^{a-1} i*b^i
      3) a est maximum
      
      Principe pour trouver a et b: on part de a=2, et on essaye de
      voir s'il existe un b. Pour cela on choisit b=a*k+1, la plus
      petite valeur possible vÃ©rifiant 1), on calcule la somme
      (max{X}) et on vÃ©rifie qu'elle ne dÃ©passe pas p/(k-1). Si la
      somme est correcte, on passe Ã  a+1, sinon on a trouvÃ© a et b.
     */

    int const s0=iceil(p,k-1); // NB: s<p/(k-1) => s<Ceil(p/(k-1))

    /* cherche a et b */
    i=a=1,s=0,b=k+1; // NB: b=a*k+1, w=b^i, s=sum i*b^i=max{X}
    while(s<s0)
      if(i==a) a++, b+=k, w=1, i=s=0;
      else s+=i*w, w*=b, i++;

    a--, b-=k; // ici a et b vÃ©rifient 1),2),3)

    /* crÃ©e l'ensemble X */
    NALLOCZ(int,P,a,_i); // P = permutation sur [0,a[
    ALLOC(Q->wrap,p); // Q->wrap = ensemble X de taille <= p

    s=1; // s=indice courant de Q->wrap, i=0 est rÃ©servÃ© pour (a!)
    do{
      Q->wrap[s]=0, w=1;
      for(i=0;i<a;i++) Q->wrap[s]+=P[i]*w, w*=b;
      s++;
    }
    while(NextPermutation(P,a,NULL));
    free(P);
    REALLOC(Q->wrap,s); /* ici s=a!+1, s<=p */
    Q->wrap[0]=--s; // wrap[0]=s=|X|, wrap[1..s]=X

    /* trie X pour test d'adjacence plus rapide, en log(s) */
    QSORT(Q->wrap+1,s,fcmp_int); // ordre croissant
    return 0;

  case QUERY_ADJ:;
    int si=Q->i/p; // stable de i
    int sj=Q->j/p; // stable de j
    if(si==sj) RET_a(0); // pas adjacent
    a=Q->i%p; // rang de i dans son stable
    b=Q->j%p; // rang de j dans son stable
    if((si==0)&&(sj==k-1)){ SWAP(a,b); sj=1; } // NB: sj-si=1
    if(sj-si>1) RET_a(0); // pas adjacent
    // ici i-j ssi b=a+s mod p
    s=(b-a+p)%p;
    RET_a(bsearch(&s,Q->wrap+1,Q->wrap[0],sizeof(int),fcmp_int)!=NULL);
  }
  return 1;
}

int thetagone(query* const Q)
/*
  L'adjacence est en O(k*N).
  Modifie Q->i et Q->j.
  Utilise gabriel().
*/
{
  int const p=Q->param[1];
  int const k=Q->param[2];

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
      
  case QUERY_INIT:
    if(p<3) Q->param[1]=-1; /* p infini */
    if(k<1) RET_n(0);
    return gabriel(Q);

  case QUERY_ADJ:;
    /*
      Adjacence: pour tous les axes t, on cacule P_t(i,j)=distgone(),
      puis on dÃ©termine s'il existe un autre sommet z avec P_t(i,z)
      plus petit. Si c'est non (et que P_t(i,j) est finie), alors i
      est adjacent Ã  j, sinon ils ne le sont pas (mais j peut Ãªtre
      adjacent Ã  i !).
    */
    int t,j;
    double d;

    int const n=Q->param[0];
    double const w=Q->dparam[0];

    for(t=0;t<k;t++){ /* pour tous les axes t */
      d=distgone(Q,t,p,k,w); /* calcule P_t(i,j) */
      if(d<DBL_MAX){ /* si la distance d est finie */
	j=Q->j; // sauvegarde Q->j
	for(Q->j=0;Q->j<n;Q->j++) /* pour tous les autres sommets z (=Q->j), mÃªme supprimÃ©s ! */
	  if((Q->j!=Q->i)&&(distgone(Q,t,p,k,w)<d)) Q->j=n; /* z plus proche ? */
	/* si oui, on arrÃªte: P_t(i,j) contient z */
	if(Q->j==n) RET_a(1); /* on a pas trouvÃ© de sommet z plus proche que j */
	Q->j=j; // restaure Q->j
      }
    }

    /*
      A priori ici il n'y a pas d'arÃªte entre i et j. Il faut
      cependant tester aussi entre j et i car la distance P_t(i,j)
      n'est pas symÃ©trique.
    */
    if(Q->i>Q->j) RET_a(0); /* pas d'arÃªtes i-j ni j-i */
    SWAP(Q->i,Q->j);
    return thetagone(Q);
  }
  
  return 1;
}
  

int squashed(query* const Q)
/*
  Q->rep[i][0..k[: mot du sommet i sur {0,1,2} oÃ¹ la lettre '2'
  remplace le symbole '*'.

  Valeur expÃ©rimentale de connectivitÃ© pour p=1/3:

   n | 50 | 100 | 200 | 300 | 1000 | 10000 |
  ---o----o-----o-----o-----o------o-------o
   k |  4 |  15 |  18 |  19 |   22 |    28 |

  La probabilitÃ© que la distance de Hamming entre deux lettres soit
  non-nulle est q=2((1-q)/2)^2. C'est la probabilitÃ© d'avoir 0-1 ou
  1-0. Donc, la probabilitÃ© d'arÃªte P_e que deux mots donnÃ©s soient
  voisins est

            P_e = 2kq (1-q)^{k-1} = 2q/(1-q) k*(1-q)^k

  C'est en effet deux fois celle de l'arc u->v, c'est-Ã -dire k*q (la
  distance est 1 Ã  l'une des k positions) et (1-q) pour les k-1
  autres. Par analogie avec les graphes d'Erdos-Reny, lorsque P_e >
  ln(n)/n, alors le graphe est connexe. On a:

  P_e > ln(n)/n  <=>  k*(1-q)^k > ((1-q)/(2q))*ln(n)/n
                 <=>  k*ln(1-q)*e^(k*ln(1-q)) < ln(1-q)*((1-q)/(2q))*ln(n)/n
		 <=>  k*ln(1-q) < W( ln(1-q)*((1-q)/(2q))*ln(n)/n )
		 <=>  k < W(...)/ln(1-q)

  oÃ¹ W(x) est la fonction W de Lambert. Elle vÃ©rifie:

            W(x)=y  <=> x=y*e^y

  Elle est croissante pour x>=-1/e. W(-1/e)=-1, W(0)=0, W(e)=1. Si
  x>=e, alors A(x) + (1/2)*B(x) <= W(x) <= A(x) + (e/(e-1))*B(x) oÃ¹
  A(x)=ln(x)-ln(ln(x)) et B(x)=ln(ln(x))/ln(x). En fait, pour x>=e,

            W(x) ~ ln(x) - ln(ln(x)) + ln(ln(x))/ln(x).

  Si 0<x<1, alors W(1)*x < W(x) < x avec
  W(1) = 0.5671432904097838729999686622... la constante Omega.
*/
{
  int k=Q->param[1];
  int u,t;

  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    if(k>NAMEMAX) Erreur(17);
    VIDE(Q->name);
    t=0;
    while(k--) 
      if(Q->rep[Q->i][k]==2) t += sprintf(Q->name+t,"%c",'*');
      else t += sprintf(Q->name+t,"%i",Q->rep[Q->i][k]);
    return 0;

  case QUERY_INIT:
    SET_n(Q->param[0]);
    double p=Q->dparam[0]; /* proba d'avoir '*' */
    if(k<0){ // valeur par dÃ©faut pour k
      if((p==0)||(p==1)) RET_n(0); // il faut 0<p<1
      double q=(1-p)/2; q=2*q*q; // on veut: q = 2*((1-p)/2)^2;
      double x=log(1-q)*((1-q)/(2*q))*log(Q->n)/Q->n;
      // on veut: k > W(x)/ln(1-q)
      //printf("x=%lf\n",x);
      x=max(-1/M_E,x); // on force x>=-1/e
      // calcule une borne sup sur W(x)
      // si x<e, W(x)>W(1)*x, sinon W(x) >= ln(x) - ln(ln(x)) + 0.5*ln(ln(x))/ln(x)
      if(x<M_E) x*=0.5671;
      else{
	x=log(x)-log(log(x))+0.5*log(log(x))/log(x);
	x /= log(1-q); // x=W(x)/ln(1-q)
      }
      Q->param[1]=k=floor(1+x); // pour avoir k > x
      //printf("k=%i p=%lf q=%lf\n",k,p,q);
    }
    if((k<1)||(k>=Q->n)||(Q->n<2)) RET_n(0);
    ALLOC2(Q->rep,Q->n,k);
    for(u=0;u<Q->n;u++) /* pour tous les sommets u du graphe */
      for(t=0;t<k;t++) /* pour toutes les lettres */
	if(RAND01<p) Q->rep[u][t]=2; /* lettre '*' avec proba p */
	else Q->rep[u][t]=(RAND01<0.5); /* sinon lettre 0 ou 1 */
    return 0;
    
  case QUERY_ADJ:
    /* est-ce que la distance de Hamming entre rep[i] et rep[j] vaut 1 ? */
    Q->a=0;
    while(k--)
      if(Q->rep[Q->i][k]==1-Q->rep[Q->j][k]){
	Q->a++;
	if(Q->a>1) RET_a(0); /* pas adjacent */
      }
    return 0; /* ici Q->a=0 ou 1 */
  }

  return 1;
}


int udg(query* const Q)
/*
  NB: deux points peuvent avoir les mÃªme coordonnÃ©es et la norme n'est
  pas forcÃ©ment symÃ©trique en i,j.
*/
{
  switch(Q->code){

  case QUERY_END:
  case QUERY_INIT: 
    return gabriel(Q);

  case QUERY_ADJ:
    if(dist_ij(Q)<=Q->dparam[0]) RET_a(1);
    if(Q->i>Q->j) RET_a(0);
    SWAP(Q->i,Q->j);
    return udg(Q);
  }

  return 1;
}


int rng(query* const Q)
/*
  Adjacence en O(N).
  Utilise gabriel().
  Modifie Q->i, Q->j.
*/
{
  if(Q->code==QUERY_ADJ){
    int const n=Q->param[0];
    double const r=dist_ij(Q); // !!! pas forcÃ©ment symÃ©trique en i,j
    int const i=Q->i; // sauvegarde Q->i et Q->j
    int const j=Q->j;
    double d1,d2;
    
    for(Q->j=0;Q->j<n;Q->j++){ // teste mÃªme les sommets supprimÃ©s
      Q->i=i; d1=dist_ij(Q); // d1=dist(i,k)
      Q->i=j; d2=dist_ij(Q); // d2=dist(j,k)
      if(max(d1,d2)<r){
	if(i>j) RET_a(0); // norme pas forcÃ©ment symÃ©trique
	Q->i=j,Q->j=i; // inverse Q->i et Q->j
	return rng(Q);
      }
    }
    RET_a(1);
  }

  return gabriel(Q);
}


int knng(query* const Q)
/*
  Adjacence en O(N).
  Utilise gabriel().
  Modifie Q->i, Q->j.
*/
{
  if(Q->code==QUERY_ADJ){
    int const n=Q->param[0];
    int const k=Q->param[1];
    double const r=dist_ij(Q); // !!! pas forcÃ©ment symÃ©trique en i,j
    int const j=Q->j; // sauvegarde Q->j pour utiliser dist_ij(Q)
    int c=0; // compteur de points Ã  distance < r de i

    /* On compte combien de points j' sont Ã  distance < r de i. Si j
       est parmi les k plus proches, il ne devrait pas avoir plus de k
       sommets strictement plus proche (en comptant i Ã  distance 0 <
       r). */

    for(Q->j=0;Q->j<n;Q->j++) /* teste mÃªme les sommets supprimÃ©s */
      if(dist_ij(Q)<r){ if(c==k) Q->j=n; else c++; } /* alors d(k,i)<d(i,j) */
    
    if(Q->j==n) RET_a(1); // ici j est parmi les k plus proches de i

    /* Avant de dire que i a pour voisin j, il faut tester si j a pour
       voisin i car le test n'est pas symÃ©trique. */

    if(Q->i>j) RET_a(0);
    Q->j=Q->i, Q->i=j; // inverse Q->i et Q->j
    return knng(Q);
  }
  
  return gabriel(Q);
}


int hexagon(query* const Q)
/*
  Utilise i<j.
    
  On voit le graphe comme une grille de p+1 lignes de 2q+2 colonnes
  allant du coin en bas Ã  droite (0,0) au coin (p,2q+1) (en notation
  (ligne,colonne)), et dans laquelle deux sommets ont Ã©tÃ© supprimÃ©s: le
  coin (0,2q+1) et le coin (p,2q+1) si p est impair, le coin (p,0) si
  p est pair. Les numÃ©ros sont consÃ©cutifs sur une ligne, de haut en
  bas.

  Ex:

  hexagon 3 2   hexagon 2 3

  o-o-o-o-o x   o-o-o-o-o-o-o x
  |   |   |     |   |   |   |
  o-o-o-o-o-o   o-o-o-o-o-o-o-o
    |   |   |     |   |   |   |
  o-o-o-o-o-o   x o-o-o-o-o-o-o
  |   |   |
  o-o-o-o-o x

*/
{
  int li,lj,ci,cj;

  int const p=Q->param[0];
  int const q=Q->param[1];
  int const t=(q<<1)+2; /* longueur d'une ligne */

  switch(Q->code){

  case QUERY_INIT:
    RET_n((p+1)*t-2);

  case QUERY_ADJ:;
    // pour Ã©viter de modifier Q
    int i=Q->i;
    int j=Q->j;

    if(i>=t-1) i++; /* on insÃ¨re virtuellement le coin (0,2q+1) */
    if(j>=t-1) j++;
    if((p&1)==0){ /* si p est impair, on a rien Ã  faire */
      if(i>=p*t) i++; /* on insÃ¨re virtuellement le coin (p,0) si p est pair */
      if(j>=p*t) j++;
    }

    /* on calcule les coordonnÃ©es de i et j placÃ©s sur cette
       grille (avec les coins manquant) */

    li=i/t;ci=i%t;
    lj=j/t;cj=j%t;
  
    /* utilise le fait que i<j: dans le dernier cas lj=li+1 */
    RET_a(
	  ((li==lj)&&(abs(ci-cj)==1)) ||
	  ((ci==cj)&&(lj==(li+1))&&((ci&1)==(li&1)))
	  );
  }

  return 1;
}


int whexagon(query* const Q)
/*
  Utilise i<j et hexagon(Q).
  Modifie Q->i et Q->j.
*/
{
  int li,ci,lj,cj;
  int const p=Q->param[0];
  int const q=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    hexagon(Q);
    RET_n(Q->n+p*q);

  case QUERY_ADJ:;
    int t=Q->n-p*q; /* t=nombre de sommets de l'hexagone */

    /* teste si i et j sont dans l'hexagone */
    if((Q->i<t)&&(Q->j<t)) return hexagon(Q);

    /* teste si i et j sont hors l'hexagone */
    if((Q->i>=t)&&(Q->j>=t)) RET_a(0);

    /* on a i dans l'hexagone et j en dehors car i<j */
    lj=(Q->j-t)/q;cj=(Q->j-t)%q; /* j est le centre de l'hexagone (lj,cj) */
    t=(q<<1)+2; /* t=longueur d'une ligne de l'hexagone */
    
    /* on calcule les coordonnÃ©es de i dans l'hexagone */
    if(Q->i>=t-1) Q->i++; /* on corrige */
    if(((p&1)==0)&&(Q->i>=p*t)) Q->i++; /* on recorrige */
    li=Q->i/t;ci=Q->i%t; /* (li,ci) coordonnÃ©es de i */
    
    RET_a( ((li==lj)||(li==(lj+1))) && (abs(2*cj+1+(lj&1)-ci)<2) );
  }

  return 1;
}


int hanoi(query* const Q)
/*
  Modifie Q->i et Q-j.

  Adjacence: on Ã©crit i (et j) en base b, mot de n lettres. i et j
  sont adjacents ssi i=Puv...v et j=Pvu...u oÃ¹ P est un prÃ©fixe
  commun, et u,v des lettres qui se suivent (modulo b).
*/
{
  int ri,rj,u,v,k;
  int const n=Q->param[0];
  int const b=Q->param[1];
  switch(Q->code){

  case QUERY_NAME:
    if(b<10) name_base(Q->name,Q->i,b,n,"","",1);
    else name_base(Q->name,Q->i,b,n,",","",1);
    return 0;

  case QUERY_INIT:
    if(b<2) RET_n(0);
    for(Q->n=k=1;k<=n;k++) Q->n *= b; /* #nombre de sommets = b^n */
    return 0;

  case QUERY_ADJ: /* on Ã©graine les chiffres en base b */
    
    for(ri=rj=k=0;k<n;k++){
      if(ri!=rj) break; /* on s'arrÃªte dÃ¨s qu'on diffÃ¨re, on garde k */
      ri=Q->i%b; Q->i/=b; /* ri=dernier chiffre, i=i sans dernier chiffre */
      rj=Q->j%b; Q->j/=b; /* rj=dernier chiffre, j=j sans dernier chiffre */
    }
    if((((ri+1)%b)!=rj)&&(((rj+1)%b)!=ri)) RET_a(0); /* alors pas voisin */

    u=ri; v=rj; /* ici u et v sont consÃ©cutifs (mod b) */
    for(;k<n;k++){
      ri=Q->i%b; Q->i/=b;
      rj=Q->j%b; Q->j/=b;
      if(ri!=v) RET_a(0); /* pas bon */
      if(rj!=u) RET_a(0); /* pas bon */
    }
    RET_a(1);
  }

  return 1;
}


int sierpinski(query* const Q)
/*
  On utilise Q->rep[i][k] pour reprÃ©senter le sommet i. C'est un mot
  d'au plus n lettres de [0,b[. On pose rep[i][n]=L oÃ¹ L est la
  longueur du mot.

  Le mot du sommet i reprÃ©sente la suite des cycles auxquels il
  appartient, sauf s'il est l'un des b sommets initiaux. Prennons
  b=3. Les 3 sommets du triangle sont 0,1,2. Si n>1, les autres
  sommets commenceront tous par 0. On a alors 3 autres triangles,
  numÃ©rotÃ©s 0,1,2, le triangle i Ã  pour sommet le sommet i. Les 3
  sommets internes sont 00,01,02. Le sommet 00 partage les triangles 0
  et 1, 01 les triangles 1 et 2, 02 les triangles 2 et 0. Si n>2, tous
  les autres sommets commenceront par 00, 01 ou 02 suivant le
  sous-triangle auquel ils appartiennent. Dans le sous-triangle 0, les
  3 sommets internes seront 000, 001, 002. Etc.

  Ex: n=2 b=3         0 
                     /  \
                   00 -- 02
                  /  \  /  \
                 1 -- 01 -- 2

  Adjacence:

  CAS 1: extrÃ©mitÃ© (|i|=1 et |j|=n)
  Si i=x et n>1, alors il est voisin de:
  - 0 x^{n-2} x
  - 0 x^{n-2} (x-1)
  Si i=x et n=1, alors c'est le CAS 2.
  
  Soit P=le plus prÃ©fixe commun entre i et j, et k=|P|

  CAS 2: sommet du triangle le plus interne (|i|=|j|=n)
  Si i=Px, k=n-1, et n>0, alors il est voisin de:
  - P (x+1)
  - P (x-1)

  CAS 3: sommet entre deux triangles (1<|i|<n et |j|=n)
  Si i=Px, alors:

  CAS 3.1: i=P=Qx est prÃ©fixe de j (k=|i|).
  Alors il est voisin de:
  - P (x+1)^{n-k-1} x
  - P (x+1)^{n-k-1} (x+1)

  CAS 3.2: i=Px.
  Alors il est voisin de:
  - P (x+1) x^{n-p-2} x
  - P (x+1) x^{n-p-2} (x-1)

*/
{
  int k,t,r,x,c,li,lj;
  int const n=Q->param[0];
  int const b=Q->param[1];
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    if(b<10) name_vector(Q->name,Q->rep[Q->i],Q->rep[Q->i][n],"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],Q->rep[Q->i][n],",","()",1,"%i");
    return 0;

  case QUERY_INIT:
    if((b<3)||(n<1)) RET_n(0); /* graphe vide, graphe non dÃ©fini */
    for(Q->n=b,k=2;k<=n;k++) Q->n=b*Q->n-b;
    ALLOC2(Q->rep,Q->n,n+1); /* rep[t][k]=k-Ã¨me lettre du sommet t */
    for(t=0;t<Q->n;t++){ /* calcule les noms */
      x=t;
      k=r=0;
      if(x<b) Q->rep[t][0]=x; /* un des b sommets du 1er cycle ? */
      else{
	x -= b;	
	for(;;){
	  Q->rep[t][k++]=r;
	  if(x<b) { Q->rep[t][k]=x; break; }
	  x -= b; r=x%b; x /= b;
	}
      }
      Q->rep[t][n]=k+1-(n==0); /* longueur du mot, corrige si n=0 */
    }
    return 0;

  case QUERY_ADJ:
    li=Q->rep[Q->i][n]; /* longueur de i */
    lj=Q->rep[Q->j][n]; /* longueur de j */
    /* propriÃ©tÃ©: si i<j, alors li<=lj */
    if(lj<n) RET_a(0);
    
    /* CAS 1 */
    if((li==1)&&(n>1)&&(Q->rep[Q->j][0]==0)){
      x=Q->rep[Q->i][0];
      for(c=t=1;t<=n-2;t++) if(Q->rep[Q->j][t]!=x) { c=0; break; }
      if(c){ /* c=vrai ssi j=0x^(n-2) */
	if(Q->rep[Q->j][n-1]==x) RET_a(1);
	if(Q->rep[Q->j][n-1]==((x-1+b)%b)) RET_a(1);
      }
    }

    /* calcule k=longueur du prÃ©fixe commun */
    for(k=0;(k<li)&&(k<lj)&&(k<n);k++)
      if(Q->rep[Q->j][k]!=Q->rep[Q->i][k]) break;

    /* CAS 2 */
    if((li==n)&&(k==n-1)){
      x=Q->rep[Q->i][k];
      if(Q->rep[Q->j][k]==((x+1)%b)) RET_a(1);
      if(Q->rep[Q->j][k]==((x-1+b)%b)) RET_a(1);
    }

    /* CAS 3 */
    if((li==1)||(li==n)) RET_a(0);
    x=Q->rep[Q->i][li-1];
    /* ici on a 1<|i|<n, |j|=n, et x=derniÃ¨re lettre de i */

  /* CAS 3.1 */
  if(k==li){
    for(t=k;t<=n-2;t++) if(Q->rep[Q->j][t]!=((x+1)%b)) break;
    if(t>n-2){
      if(Q->rep[Q->j][n-1]==x) RET_a(1);
      if(Q->rep[Q->j][n-1]==((x+1)%b)) RET_a(1);
    }
  }
  
  /* CAS 3.2 */
  if((k==li-1)&&(Q->rep[Q->j][k]==((x+1)%b))){
    for(t=k+1;t<=n-2;t++) if(Q->rep[Q->j][t]!=x) break;
    if(t>n-2){
      if(Q->rep[Q->j][n-1]==x) RET_a(1);
      if(Q->rep[Q->j][n-1]==((x-1+b)%b)) RET_a(1);
    }
  }
  
  RET_a(0); /* si i>j, alors on ne fait rien */
  }

  return 1;
}


int banana(query* const Q)
/*
  Utilise i<j.

  Chaque branche Ã  son centre numÃ©rotÃ© 0 et est connectÃ© au sommet n-1
  par sa branche numÃ©ro 1. Cela donne pour n=k=3:


                  23 23 23
                  \| \| \|
                   0  0  0
                   |  |  |
                   1  1  1
                    \ | /
                     n-1
*/
{
  int const n=Q->param[0];
  int const k=Q->param[1]+1; // Atention! c'est +1
  switch(Q->code){

  case QUERY_INIT:
    if((n<=0)||(k<1)) RET_n(0);
    RET_n(n*k+1);
    
  case QUERY_ADJ:
    if(Q->j==Q->n-1) RET_a(Q->i%k==1);
    RET_a((Q->i/k==Q->j/k)&&(Q->i%k==0));
  }
  return 1;
}


int rpartite(query* const Q)
/*
  Utilise i<j.

  Les sommets sont numÃ©rotÃ©s consÃ©cutivement dans chacune des
  parts. Utilise le tableau Q->wrap en interne: wrap[k]=a_1+...+a_k
  est la k-iÃ¨me somme partielle, avec wrap[0]=0. Donc les sommets de
  la part i sont numÃ©rotÃ©s de wrap[i-1] Ã  wrap[i] (exclu).
*/
{
  int k,s,r=Q->param[0];
  switch(Q->code){

  case QUERY_END:
    free(Q->wrap);Q->wrap=NULL;
    return 0;

  case QUERY_INIT:
    if(r<=0) RET_n(0);
    free(Q->wrap);ALLOC(Q->wrap,r);
    Q->n=Q->wrap[0]=0;
    for(k=1;k<=r;k++){
      if(Q->param[k]<=0) RET_n(0);
      Q->n += Q->param[k];
      Q->wrap[k]=Q->n;
    }
    return 0;
    
  case QUERY_ADJ:
    /*
      Pour calculer l'adjacence, avec i<j, on calcule les numÃ©ros de
      la part de i et de j: adjacent ssi part(i)<>part(j). Pour cela,
      on fait une recherche dichotomique de la part de i (c'est-Ã -dire
      d'un k tq wrap[k]<=i<wrap[k+1]). La complexitÃ© est O(log(r)),
      r=#parts.
    */
    
    /* Cherche la part de i dans [s,r[: au dÃ©part c'est dans [0,r[ */
    s=0;
    k=r/2;
    while(s<r){
      if(Q->i<Q->wrap[k]){
	if(Q->j>=Q->wrap[k]) RET_a(1); /* i et j sont dans des parts <> */
	r=k; /* ici i<j<wrap[k]: on cherche dans [s,k[ */
	k=(s+k)>>1;
      }else{ /* ici wrap[k]<=i<j */
	if(Q->j<Q->wrap[k+1]) RET_a(0); /* i et j sont dans la part k */
	s=k; /* ici wrap[k]<=i<j: on cherche dans [k,r] */
	k=(k+r)>>1; 
      }
    }
    RET_a((Q->j>=Q->wrap[k+1])); /* ici i est dans la part k. On vÃ©rifie si j aussi */
  }
  
  return 1;
}


int aqua(query* const Q)
/*
  On se sert de rep[i][0..n], n=param[0].
*/
{
  int n=Q->param[0]; /* n=nombre de paramÃ¨tres */
  int *C=Q->param+1; /* C=tableau de contraintes */
  int k,x,y,t;
  
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],n,",","",1,"%i");
    return 0;
    
  case QUERY_INIT:
    x=*C; /* s=Q->param[1]=premier terme=la somme */
    if(n<0) n=0;
    int *S=NextPart(NULL,n,x,C);
    Q->n=0; /* calcule une fois Q->n pour ALLOC2() */
    do Q->n++; while(NextPart(S,n,x,C));
    ALLOC2(Q->rep,Q->n,n); /* ici Q->n>0 */
    Q->n=0; /* calcule les sommets */
    do{
      for(k=0;k<n;k++) Q->rep[Q->n][k]=S[k];
      Q->n++;
    }while(NextPart(S,n,x,C));
    free(S);
    return 0;

  case QUERY_ADJ:
    if(Q->i==Q->j) RET_a(0);

    /* compte et mÃ©morise les diffÃ©rences entre rep[i] et rep[j]: il y en a au moins 2 */
    for(k=0,x=y=-1;k<n;k++)
      if(Q->rep[Q->i][k]!=Q->rep[Q->j][k]){
	if((x>=0)&&(y>=0)) RET_a(0); /* si plus de 2 diffÃ©rences, alors pas d'arc */
	if(x<0) x=k; else if(y<0) y=k;
      }
    
    /* soit on a versÃ© x vers y, soit le contraire */
    /* k=quantitÃ© que l'on peut verser de x Ã  y */
    for(t=0;t<2;t++){ // les deux cas sont identiques en inversant x et y
      k=min(C[y],Q->rep[Q->i][x]+Q->rep[Q->i][y]) - Q->rep[Q->i][y];
      if((Q->rep[Q->j][y]==Q->rep[Q->i][y]+k)&&
	 (Q->rep[Q->j][x]==Q->rep[Q->i][x]-k)) RET_a(1);
      SWAP(x,y); // inverse x et y
    }
    RET_a(0);    
  }

  return 1;
}


int flower_snark(query* const Q)
{
  int i=Q->i,j=Q->j; // pour ne pas modifier Q->i et Q->j
  int const u=(j>>2)-(i>>2); /* i<j, donc u>=0 */
  int const k=Q->param[0];
  switch(Q->code){

  case QUERY_NAME:
    sprintf(Q->name,"%c%i",((i&3)==0)? 'c' : 't'+(i&3),i>>2);
    return 0;
    
  case QUERY_INIT:
    RET_n(k<<2);

  case QUERY_ADJ:
    i &= 3;
    j &= 3;
    if(u==0) RET_a((i==0));
    if((u==1)&&(i==j)) RET_a((i>0));
    if(u!=k-1) RET_a(0);
    i*=j;
    RET_a(((i==1)||(i==6)));
  }

  return 1;
}


int gear(query* const Q)
/*
  Utilise i<j.
  Utilise cage().
  Attention! modifie Q->param[0].
*/
{
  switch(Q->code){

  case QUERY_INIT:
    Q->param[1]=Q->param[0]<<1; // 2n
    Q->param[0]=3; Q->param[2]=2; Q->param[3]=0;
    RET_n(Q->param[1]+1); // 2n+1 sommets
    
  case QUERY_ADJ:
    if(Q->j<Q->param[1]) return cage(Q);
    RET_a( (Q->j==Q->param[1])&&((Q->i&1)==0) );
  }
  
  return 1;
}


int parachute(query* const Q)
/*
  Utilise i<j.
  Utilise fan().
  ParamÃ¨tre cachÃ©: Q->param[1].
*/
{
  switch(Q->code){

  case QUERY_INIT:
    Q->param[1]=2;
    RET_n(Q->param[0]+3);
    
  case QUERY_ADJ:
    if(Q->j==Q->n-1) RET_a(Q->i==Q->n-2);
    return fan(Q); // ici i<j<n-1
  }
  
  return 1;
}


int arboricity(query* const Q)
/*
  Marche en orientÃ©.
  
  Utilise Q->rep[i][0..k[ pour la reprÃ©sentation implicite, les k
  pÃ¨res du sommet i, oÃ¹ k=Q->param[1]. Si Q->rep[i][j]<0, c'est que le
  pÃ¨re j de i n'existe pas (i est une racine de la forÃªt j par
  exemple).
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT: /* calcule Q->n et Q->rep[i][0..k-1] */
    SET_n(Q->param[0]);
    int const k=Q->param[1];
    if(k<1) RET_n(0);       /* il faut k>0 */
    Q->k=1; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,k);  // rep=reprÃ©sentation finale
    NALLOCZ(int,P,Q->n,_i); // P=permutation alÃ©atoire
    NALLOC(int,T,Q->n);     // T=arbre alÃ©atoire
    int t,v;
    
    for(t=0;t<k;){ /* pour chaque arbre */
      Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre alÃ©atoire T */
      T[0]=0; /* racine = 0 plutÃ´t que -1, car P[-1] n'existe pas ! */
      for(v=0;v<Q->n;v++) Q->rep[P[v]][t]=P[T[v]]; /* copie et permute le pÃ¨re */
      Q->rep[P[0]][t]=-1; /* pÃ¨re de la racine = -1 */
      if(++t<k) Permute(P,Q->n); /* si pas fini, on permute alÃ©atoirement P */
    }
    
    free(T);
    free(P);
    return 0;
    
  }

  return 1;
}


int kpage(query* const Q)
/*
  Utilise rep[] pour la reprÃ©sentation implicite. Chaque sommet i
  possÃ¨de 2k pÃ¨res: rep[i][2p] et rep[i][2p+1] sont les 2 pÃ¨res du
  sommet i de la page p (p-Ã¨me outerplanar), p=0...k-1. Pour le
  gÃ©nÃ©rer, on fait l'union de k outerplanars connexes enracinÃ©s et
  plan ("outer-plan" en fait). Chacun est numÃ©rotÃ©s selon un parcours
  de la face extÃ©rieure avec une permutation circulaire alÃ©atoire,
  sauf le premier. Pour l'adjacence on fait comme pour un graphe
  d'arboricitÃ© 2k. Le pagenumber d'un ne peut dÃ©passer ceil(n/2), mais
  des valeurs de k plus grandes sont quand mÃªme autorisÃ©es.

  Le tirage de chaque outerplanar est construit Ã  partir d'un arbre
  bicoloriÃ© (couleur 0 ou 1). Le 2e parent de u existe si u est
  coloriÃ© (disons Ã  1). Il vaut alors le prochain sommet non-related
  de u. On ne tient pas compte de la derniÃ¨re branche, ces sommets
  n'ayant jamais de 2e parent. Ce tirage est uniforme sur les
  "outer-plan" connexes.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int u,p,q,t,z,k,c;
    k=Q->param[1]; /* k=nombre de pages <= ceil{n/2} */
    if(k<1) RET_n(0); /* il faut k>0 */
    Q->k=(k<<1); // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,Q->k);
    NALLOC(int,T,Q->n); /* T=arbre alÃ©atoire */
    /*
      On veut, pour chaque page p:
      rep[u][2p+0]=pÃ¨re1 de u
      rep[u][2p+1]=pÃ¨re2 de u
    */

    for(p=q=0;p<k;p++,q++){ /* pour chaque page p=0...k-1. Attention! q=2*p */
      Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre DFS alÃ©atoire T */
      
      /* c=permutation circulaire alÃ©atoire pour les noms de sommets */
      if(p) c=randomu(Q->n); else c=0; /* aucune permutation pour p=0 */

      /* calcule le pÃ¨re1 des sommets, le pÃ¨re dans T */
      Q->rep[c][q]=-1; /* la racine n'a pas de pÃ¨re */ 
      for(t=1;t<Q->n;t++) /* parcoure les sommets t de T selon le DFS, sauf la racine */
	Q->rep[(t+c)%Q->n][q]=(T[t]+c)%Q->n;

      /* calcule le pÃ¨re2 des sommets */

      /* Principe: le sommet courant est t. Chaque fois qu'on dÃ©marre
	 une nouvelle branche (t-T[t]>1), alors on parcoure les
	 sommets u allant de t-1 Ã  T[t] non compris (la branche donc)
	 et dÃ©cide pour chaque u de le connecter ou pas vers t (t qui
	 est donc le prochain non-related de u). NB: On ne parcoure
	 qu'au plus deux fois chacun des sommets. */

      q++; /* pour le pÃ¨re2 */
      Q->rep[c][q]=-1; /* la racine n'a pas de pÃ¨re */
      for(t=1;t<Q->n;t++){ /* parcoure les sommets t de T selon le DFS, sauf la racine */
	u=(t+c)%Q->n; /* u=sommet du graphe correspondant Ã  t */
	Q->rep[u][q]=-1; /* par dÃ©faut, pas de pÃ¨re2 */
	if(t-T[t]>1){ /* ici t dÃ©marre une nouvelle branche */
	  z=t-1; /* z=dernier sommet de la branche prÃ©cÃ©dante */
	  while(z!=T[t]){ /* z parcoure la branche prÃ©cÃ©dante */
	    if(RANDbit) Q->rep[(z+c)%Q->n][q]=u; /* ajoute un voisin vers u si colorÃ© */
	    z=T[z]; /* z descend le long de la branche */
	  }
	}
      }
    }

    DEBUG(
	  for(p=0;p<k;p++){
	    ruling("â€•",10);printf("\n");
	    printf("page %i:\n",p);
	    for(u=0;u<Q->n;u++)
	      printf("%i: %i %i\n",u,Q->rep[u][2*p],Q->rep[u][2*p+1]);
	  }
	  ruling("â€•",10);printf("\n");
	  );

    free(T);
    return 0;
  }

  return 1;
}


int cactus(query* const Q)
/*
  Q->rep[i][0..1] sont les 2 pÃ¨res du sommet i. L'algorithme est le
  mÃªme que pour un outerplanar (soit un k-page avec k=1) sauf qu'un
  sommet u dÃ©marrant une nouvelle branche ne peut avoir qu'au plus
  seul voisin z dans la branche prÃ©cÃ©dante. Mais cela ne suffit
  pas. Il faut de plus qu'aucun sommet entre z et T[u] n'ait dÃ©jÃ  de
  voisin dans sa branche prÃ©cÃ©dante.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int u,z,t;
    Q->k=2; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,2);
    NALLOC(int,T,Q->n); /* T=arbre DFS alÃ©atoire */
    Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre DFS alÃ©atoire T */
    for(u=0;u<Q->n;u++) Q->rep[u][0]=T[u],Q->rep[u][1]=-1; /* pÃ¨re0 = T[u], pÃ¨re1 = -1 au dÃ©part */
    for(u=1;u<Q->n;u++){ /* parcoure les sommets u de T selon le DFS, sauf la racine */
      if(u-T[u]>1){ /* ici u dÃ©marre une nouvelle branche */
	z=u-1; /* z=dernier sommet de la branche prÃ©cÃ©dante */
	t=0; /* t=sommet de la branche vers lequel u va pointer */
	while(z!=T[u]){
	  if((t==0)&&RANDbit) t=z; /* u veut choisir z */
	  if((t)&&(Q->rep[z][1]>0)) t=0; /* pas bon, z dessous t Ã  dÃ©jÃ  un voisin */
	  z=T[z]; /* sommet suivant de la branche */
	}
	if(t) Q->rep[u][1]=t; /* c'est bon, u peut choisir t */
      }
    }
    free(T);
    return 0;
  }

  return 1;
}


int planar(query* const Q)
/*
  Q->rep[u][0..1] = reprÃ©sentation implicite de u, ses au plus deux
  pÃ¨res.
*/
{
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);
    
  case QUERY_INIT:;
    int n=Q->param[0]; // n=nombre de faces
    int f=Q->param[1]; // f=taille maximum des faces internes
    int d=Q->param[2]; // d=degrÃ© des sommets internes
    int const w=(f<0); // w=vrai ssi face de taille au plus |f|

    f=abs(f);
    if((f<3)||(n<=0)) RET_n(0);

    /* Le nombre maximum de sommets du graphe est a priori
       f+(n-1)*(f-2) = n(f-2)+2: 1 cycle de taille f au dÃ©part et on
       crÃ©e, pour chaque autre face, au plus f-2 nouveaux sommets */

    int k=n*(f-2)+2;     // k=nombre max de sommets
    if(d<0) d=k;         // pas de contraintes de degrÃ©
    NALLOCZ(int,A,k,-1); // A[u]=1er voisin de u, aucun par dÃ©faut
    NALLOCZ(int,B,k,-1); // B[u]=2e voisin de u, aucun par dÃ©faut
    NALLOC(int,C,k);     // C[i]=u=i-Ã¨me sommet u de la face extÃ©rieure
    NALLOC(int,T,k);     // pour mettre Ã  jour C
    NALLOCZ(int,D,k,2);  // D[u]=degrÃ© du sommet u, 2 par dÃ©faut

    /* initialisation de la 1Ã¨re face = 1Ã¨re face extÃ©rieure */
    int fr=(w)? 3+randomu(f-2) : f; // fr=taille de la 1Ã¨re face
    for(Q->n=0;Q->n<fr;Q->n++){
      A[Q->n]=(Q->n+1)%fr; // pÃ¨re vers le prochain du cycle
      C[Q->n]=Q->n; // face extÃ©rieure
    }
    int c=Q->n; // c=|C|=nombre de sommets de la face extÃ©rieure */
    int a,b,t,s,p,u;
    
    /* ajoute les n-1 autres faces */
    /* ici Q->n est le nouveau sommet courant */

    while((--n)>0){ // tant qu'il reste une face Ã  faire
      fr=(w)? 3+randomu(f-2) : f; // fr=taille de la nouvelle face
      k=randomu(c); // C[k]=un sommet de C au hasard
      A[Q->n]=C[k]; // 1er pÃ¨re de C[k]
      D[C[k]]++; // un voisin de plus pour C[k]

      /* on va tester les voisins valides (degrÃ© interne >= d) autour
	 de C[k], puis en choisir un au hasard. Enfin on le connectera
	 par un chemin jusqu'au nouveau sommet courant, Q->n, et on
	 mettra Ã  jour la face extÃ©rieure. */

      p=min(fr-2,c/2); // p=nombre maximum de sommets de part et
                        // d'autre de C[k] dont il faut tester le degrÃ©
      
      /* teste les successeurs de C[k] sur C: C[k+1]...C[k+p] */
      for(a=1;a<=p;a++)	if(D[C[(k+a)%c]]<d) break;
      if(a>p) a--;
      // ici on peut se connecter Ã  n'importe quel sommet entre C[k+1]..C[k+a] 

      /* teste les prÃ©dÃ©cesseurs de C[k] sur C: C[k-1]...C[k-p] */
      for(b=1;b<=p;b++)	if(D[C[(k-b+c)%c]]<d) break;
      if(b>p) b--;
      // ici on peut se connecter Ã  n'importe quel sommet entre C[k-1]..C[k-b] 

      t=randomu(a+b); // choisit 1 indice parmi [0,a+b[
      if(t<a) s=1; else{ s=-1; t-=a; } 
      p=fr-t-3; // p=nombre de sommets du chemin, Q->n non compris, p=0 possible

      // t random dans [0,a[ si s>0
      // t random dans [0,b[ si s<0
      // s>0 => chemin = C[k]-N-(N+1)- ... - (N+p)-C[k+t+1]
      // s<0 => chemin = C[k-t-1]-(N+p)- ... - (N+1)-N-C[k]
      
      b=(k+s*(t+1)+c)%c; // C[b]=C[k+t+1] ou C[k-t-1]=sommet Ã  connecter
      D[C[b]]++; // un voisin de plus pour C[b]
      for(u=Q->n;u<Q->n+p;u++) B[u]=u+1; // adjacence du chemin
      B[u]=C[b]; // le dernier sommet du chemin pointe sur C[b], u=Q->n possible

      /* mise Ã  jour de C: on supprime d'abord les sommets de C qui
	 vont passer Ã  l'intÃ©rieur, puis on calcule dans T la nouvelle
	 face extÃ©rieur (en sautant les sommets effacÃ©s de C):
	 - si s>0, il faut supprimer C[k+1]...C[k+t]
	 - si s<0, il faut supprimer C[k-t]...C[k-1]
      */
      for(u=1,k+=c;u<=t;u++) C[(k+u*s)%c]=-1; // supprime des sommets de C

      /* crÃ©er la nouvelle face extÃ©rieure T, a=|T| */
      for(u=a=0;u<=p;u++) T[a++]=Q->n++; // commence par les p+1 nouveaux sommets
      for(u=0,b+=c;u<c;u++) // le reste de l'ancienne face extÃ©rieure
	if(C[(b+u*s)%c]>=0) T[a++]=C[(b+u*s)%c];

      SWAP(C,T); // Ã©change C et T
      c=a; // nouvelle taille de C
    }

    free(T);
    free(C);
    free(D);

    /* recopie A,B dans Q->rep */
    Q->k=2; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,2);
    for(u=0;u<Q->n;u++){
      Q->rep[u][0]=A[u];
      Q->rep[u][1]=B[u];
    }
    free(A);
    free(B);
    return 0;
  }

  return 1;
}


int hyperbolic(query* const Q)
/*
  A FINIR

  Q->rep[u][0..2] = reprÃ©sentation implicite de u, ces au plus trois
  pÃ¨res. Le 3e pÃ¨re ne peut se produire que si p=3.

  rep[u][0]=successeur(u) sur le dernier niveau (cycle)
  rep[u][1]=parent(u) vers le niveau prÃ©cÃ©dant (-1 si non dÃ©fini)
  rep[u][2]=parent(u) vers le succ du niveau prÃ©cÃ©dant (-1 si non dÃ©fini)

*/
{
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);
      
  case QUERY_INIT:;
    int p=Q->param[0]; // p=taille des faces
    int k=Q->param[1]; // k=degrÃ© des sommets
    int h=Q->param[2]; // h=nombre de niveaux
    if((p<3)||(k<2)||(h<1)) RET_n(0);
    if(k==2) h=1;

    int t,c,u,i,m;
    int n2=p; // =nombre de sommets de degrÃ© 2
    int n3=0; // =nombre de sommets de degrÃ© 3
    int n4=0; // =nombre de sommets de degrÃ© 4
    Q->n=p;   // =nombre total de sommets au dÃ©part

    // calcule le nombre de sommets: pour chaque niveau, on compte le
    // nombre de sommets de degrÃ© 2, 3 ou 4. En gÃ©nÃ©ral, les sommets
    // de degrÃ© 3 sont engendrÃ©s par chaque sommet du niveau
    // prÃ©cÃ©dent. Et ceux de degrÃ© 2 sont engendrÃ©s par les chemins
    // entre les ces sommets de degrÃ© 3 consÃ©cutifs nouvellement
    // crÃ©es. Il faut ensuite tenir compte des cas particuliers. Pour
    // p=3, il faut fusionner le dernier et le premier sommets de
    // degrÃ© 3. ApparaÃ®t alors un sommet de degrÃ© 4. Pour k=3, il faut
    // tenir compte des sommets de degrÃ© 3 de la couche prÃ©cÃ©dante qui
    // doivent Ãªtre sautÃ©s. Le cas p=k=3 doit Ãªtre gÃ©rÃ© diffÃ©remment
    // encore.

    for(c=1;c<h;c++){ // pour chaque nouveau niveau
      t=n3+n2+n4; // total de la couche prÃ©cÃ©dante
      u=n3;
      PRINT(n2);
      PRINT(n3);
      PRINT(n4);
      PRINT(t);
      printf("\n");
      if(t==0){ h=c; break; } // il n'y a plus de sommets
      n3=(k-2)*n2+(k-3)*n3+(k-4)*n4;
      if(p==3){ n3-=2*t,n4=t,n2=0; if(k==3) n3=1,n4=0; h=2; }
      else{ n2=n3*(p-3) - t; if(k==3) n2-=u; }
      Q->n += n2+n3+n4;
    }
    PRINT(h);
    PRINT(Q->n);

    Q->k=3; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,3);
    for(u=0;u<Q->n;u++) Q->rep[u][1]=Q->rep[u][2]=-1; // aucun parent par dÃ©faut

    // cas trÃ¨s particulier
    if((p==3)&&(k==3)){ // cycle ou K_4
      Q->rep[0][0]=1; Q->rep[1][0]=2; Q->rep[2][0]=0;
      if(h>1){
	Q->rep[3][0]=-1;
	Q->rep[0][1]=Q->rep[1][1]=Q->rep[2][1]=3;
      }
      return 0;
    } // ici p>3 ou k>3
    
    int f=0; // f=1er sommet du nouveau niveau
    int d=0; // d=1er sommet du niveau prÃ©cÃ©dent
    t=p;     // t=nouveau sommet
    
    for(c=0;c<h;c++){ // pour chaque nouveau niveau
      for(u=d;u<f;u++){ // parcoure les sommets u du niveau prÃ©cÃ©dent
	// on dÃ©termine le parent du nouveau sommet t (par dÃ©faut il
	// n'en a pas), le cycle est dÃ©terminÃ© dans un 2e temps par un
	// parcours complet du niveau (=cycle) crÃ©e
	m=k-4+(Q->rep[u][1]<0)+(Q->rep[u][2]<0); // m=k-2 ou k-3 ou k-4 suivant deg(u)
	for(i=0;i<m;i++){ // pour chaque nouveau sommet de degrÃ© 3
	  Q->rep[t][1]=u; // t sommet de deg 3 et de parent u
	  t += p-2; // prochain sommet de deg 3 (saute les sommets de deg 2)
	}
	if(m) t--;
      }
      // parcoure le nouveau niveau pour crÃ©er le cycle, NB: ici u=f
      for(;u<t;u++) Q->rep[u][0]=u+1;
      Q->rep[u-1][0]=f; // termine le cycle (Ã©vite des calculs de modulo)
      d=f; f=t; // met Ã  jour (d,f) <- (f,t)
    }
    
    return 0;
  }

  return 1;
}


int NextDyck(int *X,int n)
/*
  Calcule le prochain mot de Dyck de longueur 2n contenu dans X (qui
  doit donc Ãªtre un tableau de taille au moins 2n). Renvoie 1 ssi le
  dernier mot de Dyck a Ã©tÃ© atteint. Si n<0, alors X est initialisÃ© au
  mot X=(10)^n.  Les mots sont Ã©numÃ©rÃ©s selon le nombre croissant de 1
  consÃ©cutifs Ã  gauche.  L'algorithme est celui de
  https://github.com/cassioneri/Dyck

  Pour n=4:

  10101010 10101100 10110010 10110100 10111000 11001010 11001100
  11010010 11010100 11011000 11100010 11100100 11101000 11110000

*/
{
  int const m=2*n-1;
  int y=0;
  int x=0;
  int i;

  if(n<0){
    x=1;
    y=-(n<<1);
    for(i=0;i<y;i++){ X[i]=x; x=1-x; }
    return (y==2);
  }

  for(i=m;i>0;i--)
    if(X[i]){
      if(X[i-1]) x++;
      else{
	X[i-1]=1;
	X[i]=0;
	for(y=y-x;y!=0;y--) X[++i]=0;
	while(i<m){
	  X[++i]=1;
	  X[++i]=0;
	}
	return 0;
      }
    }else y++;

  return 1; /* dernier mot atteint */
}


int flip(query* const Q)
/*
  Q->rep[i] = mot de Dyck = [ 1 0 1 1 0 0 ]. Ici n est le nombre de 1
  du mot = Q->param[0]-2 = 3.  Le mot reprÃ©sente un arbre binaire
  complet (chaque noeud interne Ã  deux fils exactement). On obtient le
  codage de l'arbre en mot de Dyck par un parcours DFS et en Ã©crivant
  1 si l'arÃªte parcouru mÃ¨ne Ã  un fils gauche et 0 si elle mÃ¨ne Ã  un
  fils droit. On Ã©crit rien lorsqu'on parcoure les arÃªtes vers le
  pÃ¨re.

  Rotation sur le noeud x:

               B   C                      A   B
                \ /                        \ /
             A   o           ->             o   C
              \ /                            \ /
               x                              x

  U=[... 1 A 0 1 B 0 C ...]  ->  V=[... 1 1 A 0 B 0 C ...]

  Algorithme d'adjacence entre les mots U et V:

  1. On cherche le plus grand suffixe commun. Soit i la position telle
     que U[i]<>V[i] et U[i+1...2n-1] = V[i+1...2n-1]. Pour cela on
     remonte de 2n-1 jusqu'Ã  la premiÃ¨re diffÃ©rence.

  2. On Ã©change Ã©ventuellement U et V de sorte que U[i]=1 et V[i]=0.

  3. On essaye de lire [... 1 A 0 ...] dans V Ã  partir de i et dans U
     Ã  partir de i-1, toujours selon les indices dÃ©croissant. En mÃªme
     temps qu'on vÃ©rifie que V[i]=U[i-1], on calcule la hauteur h avec
     +1 si V[i]=0 et -1 si V[i]=1 (car on lit le mot Ã  l'envers). On
     s'arÃªte dÃ¨s que h<0. Soit i l'indice dans V oÃ¹ h<0 pour la
     premiÃ¨re fois.

  4. On vÃ©rifie alors que V[j-1]=1, puis que U[0..j-2]=V[0..j-2]. On
     conclut que U et V sont adjacents.

  https://fr.wikipedia.org/wiki/Nombre_de_Catalan#Chemins_sous-diagonaux_dans_le_carr.C3.A9
  https://en.wikipedia.org/wiki/Tree_rotation
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],2*(Q->param[0]-2),"","",1,"%i");
    return 0;
  
  case QUERY_INIT:;
    int n,t,u;
    /* calcule Q->n = Catalan(n), n+2 sommets */
    n=Q->param[0]-2;
    if(n<=0) RET_n(0);
    t=(n<<1); /* t=2n */
    SET_n(Binom(t,n)/(n+1));
    ALLOC2(Q->rep,Q->n,t);
    NALLOC(int,X,t);
    NextDyck(X,-n); /* initialise le 1er mot */
    t *= sizeof(int);
    for(u=0;u<Q->n;u++){
      bcopy(X,Q->rep[u],t); /* copie le mot de Dyck courrant vers Q->rep[u] */
      NextDyck(X,n); /* calcule le mot suivant */
    }
    free(X);
    return 0;
  
  case QUERY_ADJ:;
    int* U=Q->rep[Q->i];
    int* V=Q->rep[Q->j];
    int h=1; /* hauteur de [ ... 1 A 0 ... ] */
    int k=2*(Q->param[0]-2)-1; /* k=derniÃ¨re position de U ou V */

    /* calcule suffixe */
    while(U[k]==V[k]) k--;
    if(V[k]) SWAP(U,V); /* Ã©change U et V */

    /* k=position du 0 dans V */
    while((V[k]==U[k-1])&&(h>0)) h += 1-2*V[--k];

    /* problÃ¨me ? */
    if((V[--k]==0)||(h>0)) RET_a(0);

    /* prÃ©fixe */
    while((k>=0)&&(U[k]==V[k])) k--;
    RET_a((k<0));
  }

  return 1;
}


int linegraph(query* const Q)
/*
  Chaque sommet i possÃ¨de 2 couleurs prises dans [0,k[. Les sommets i
  et j sont adjacents si une couleur de l'un est une couleur de
  l'autre.  Utilise Q->rep[u][0..1] pour la reprÃ©sentation des
  couleurs du sommet u.
*/
{
  switch(Q->code){

  case QUERY_NAME:
    sprintf(Q->name,"(%i,%i)",Q->rep[Q->i][0],Q->rep[Q->j][1]);
    return 0;

  case QUERY_END:
    return free_rep(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int u,k=Q->param[1];
    if(k<=0) RET_n(0);
    ALLOC2(Q->rep,Q->n,2);
    for(u=0;u<Q->n;u++){
      Q->rep[u][0]=randomu(k);
      Q->rep[u][1]=randomu(k);
    }
    return 0;
    
  case QUERY_ADJ:
    RET_a(
	  ((Q->rep[Q->i][0]==Q->rep[Q->j][0])||(Q->rep[Q->i][0]==Q->rep[Q->j][1])||
	   (Q->rep[Q->i][1]==Q->rep[Q->j][0])||(Q->rep[Q->i][1]==Q->rep[Q->j][1]))
	  );
  }
  
  return 1;
}


int ringarytree(query* const Q)
/*
  Utilise i<j.

  Le sommet i est un chemin P(i)=x_1,x_2,... allant de la racine (=0)
  Ã  i, chaque lettre x_t est le numÃ©ro du fils, numÃ©ro dans [0,r[ pour
  la racine et dans [0,k[ pour les autres noeuds internes. P(0)={} est
  vide.

  Principe: on calcule P(i) et P(j) en parallÃ¨le, lettre par lettre,
  et on dÃ©cide de l'adjacence de i-j, avec i<j, si:

  - si p=0 (seulement connexion dans l'arbre): P(j)=P(i),x
    c'est-Ã -dire que P(j) contient une seule lettre supplÃ©mentaire par
    rapport Ã  P(i)

  - si p=1 (p=0 et chemin entre noeuds de mÃªme niveau):
    P(i)=C,x_1,...,x_k
    P(j)=C,y_1,...,y_k
    y_1=1+x_1
    y_t=0 et x_t=k-1 pour tout t>1

  - si p=2 (p=1 et cycle entre noeuds de mÃªme niveau):
    x_1=0 et y_1=r-1 ou k-1 (suivant si C={} ou pas)
    x_t=0 et y_t=k-1 pour tout t>1

*/
{
  int h=Q->param[0];
  int const k=Q->param[1];
  int const r=Q->param[2];

  switch(Q->code){

  case QUERY_INIT: /* calcule Q->n */
    if((h<0)||(k<0)||(r<0)) RET_n(0);
    if((h==0)||(r==0)){ Q->n=1; return 0; }
    if(k==0) Q->param[0]=h=1; /* si k=0 et h>0, alors la hauteur est 1 */
    if(k<=1) Q->n=1+r*h;
    else{ /* ici k>1, h>0 et r>0 */
      int t; /* taille sous-arbre = 1+r+r^2+...+r^{h-1} = (r^h-1)/(k-1) */
      for(t=0,Q->n=1;t<h;t++) Q->n *= k; /* aprÃ¨s cette boucle, Q->n=r^h */
      Q->n=1+r*(Q->n-1)/(k-1); /* Q->n=racine + r x (taille sous-arbre) */
    }
    return 0;

  case QUERY_NAME:
    if(Q->i==0) strcpy(Q->name,"Îµ"); // pour la racine
    else{
      int d=r;    /* d=nombre de fils de la racine de T, d=r puis d=k */
      int t=Q->n; /* t=taille des fils */
      int p=0;    /* nombre de caractÃ¨res Ã©crit dans Q->name */
      int const v=(max(k,r)>10); /* vrai ssi il faut une virgule */
      VIDE(Q->name);
      while(Q->i>0){
	Q->i--;t--;t/=d;
	p+=sprintf(Q->name+p,"%i",Q->i/t); /* Ã©crit le numÃ©ro du fils f=i/t */
	if(p>NAMEMAX) Erreur(17);
	Q->i%=t;
	d=k; /* maintenant d=k */
	if((v)&&(Q->i>0)) p+=sprintf(Q->name+p,","); /* ajoute une "," */
      }
    }
    return 0;
    
  case QUERY_ADJ:;
    /* calcule le prÃ©fixe commun de P(i) et P(j) */
    int const p=Q->param[3]; /* p=0,1,2 */
    int x=Q->i,y=Q->j;/* copies de i et j, NB: x<y */
    int t=Q->n; /* t=taille de T, l'arbre oÃ¹ x et y sont */
    int d=r;    /* d=nombre de fils de la racine de T, d=r puis d=k */
    int fx,fy;
    
    /* ici x et y sont dans un arbre T de taille t ayant d fils. La
       taille des fils de T est t'=(t-1)/d. Pour trouver le fils fx de
       T contenant x il suffit de faire (x-1)/t'. Le suffixe de x dans
       ce nouveau sous-arbre est alors (x-1)%t'. */
    
    fx=fy=0;
    
    /* tant que x et y sont tout deux dans T (ils viennent du mÃªme
       fils fx=fy, et aucun d'eux n'est la racine de T: on calcule
       alors leurs fils fx et fy et met Ã  jour la nouvelle taille de T
       ainsi que le suffixe de x et y dans ce nouveau T. */
    
    while((fx==fy)&&(x>0)&&(y>0)){
      x--,y--,t--; /* on enlÃ¨ve la racine de T */
      t/=d; /* t=taille des fils de T */
      fx=x/t,fy=y/t; /* fils des sous-arbres de x et de y */
      x%=t,y%=t; /* suppression du prÃ©fixe de x et de y */
      d=k; /* maintenant d=k */
    }
    
    /*
      Ici la situation est la suivante: x et y sont dans des
      sous-arbres isomorphes Ã  T de taille t, dans les sous-arbres des
      fils fx et fy. Chacun de ces sous-arbres a d fils. Si fx=fy
      c'est que x est la racine (x=0) car y=0 est impossible puisque
      x<y.
	
                       o
                   fx / \ fy
                     o   o
                    / \ / \
                     x   y

    */

    /* fx=fy: x et y sont dans le mÃªme arbre T */
    if(fx==fy) RET_a(((y-1)%((t-1)/d)==0)); /* y fils de x ? */
    if(p<=0) RET_a(0);
    
    /* fx<>fy: x et y sont dans des sous-arbres diffÃ©rents */
    
    /* si 1er et dernier sous-arbres (voisins dans le cycle), on Ã©change x et y */
    if((p==2)&&(fx==0)){
      int const b=(t==(Q->n-1)/r)? r:k; /* fx,fy fils de niveau 1 ou > 1 ? */
      if(fy==b-1){ SWAP(x,y); fy=1; }
    }
    
    if(fy-fx>1) RET_a(0); /* sous-arbres qui ne sont pas voisins */
    
    /* x et y sont dans des sous-arbres voisins, chacun de taille t */
    
    for(;;){
      if((x==0)&&(y==0)) RET_a(1); /* racines voisines */
      if((x==0)||(y==0)) RET_a(0); /* pas mÃªme niveau */
      x--,y--,t--,t/=k,fx=x/t,fy=y/t,x%=t,y%=t; /* met Ã  jour x,y,t,fx,fy */
      if((fx!=0)&&(fy!=k-1)) RET_a(0); /* il faut x fils 0 et y fils d-1 */
    }
    return 0;
  }

  return 1;
}


int n_rectree(int const h,int const *f,int const d){
/*
  Fonction auxiliaire pour rectree() renvoyant le nombre t[h] de
  sommets d'un arbre de hauteur h, chaque fils f_i Ã©tant un arbre
  hauteur h-f[i], pour i=0..d-1.

  Cette fonction utilise une technique de mÃ©morisation qui fait que le
  temps de calcul est en O(h+d), ce qui en gÃ©nÃ©ral est bien plus petit
  que le nombre final de sommets de l'arbre, soit t[h]. Si f=NULL,
  alors le tableau static interne t[] est libÃ©rÃ©. Et dans ce cas la
  fonction renvoie 0.
*/
  static int *t=NULL;
  static int hmax=-1; // hauteur max = taille-1 de t[]
  int i;

  if(f==NULL){ free(t); hmax=-1; return 0; }
  if(h<0) return 1;

  if(h>hmax){ // il faut agrandir le tableau existant
    REALLOC(t,h+1);
    for(i=hmax+1;i<=h;i++) t[i]=1;
    hmax=h;
  }
  if((h>0)&&(t[h]==1)) for(i=0;i<d;i++) t[h] += n_rectree(h-f[i],f,d);
  return t[h];
}


int rectree(query* const Q)
/*
  Utilise i<j.

  param[]={ d+2, h, f_1, ..., f_d }.

  Principe: on route vers i dans l'arbre Ã  partir de la racine, puis
  une fois atteint (c'est-Ã -dire que i se retrouve Ã  la racine du
  sous-arbre le contenant), on vÃ©rifie si j est un fils de i. Pour
  cela, on s'appuye sur n_rectree(p-f[k],...) donnant le nombre de
  sommets du sous-arbre k de i, et oÃ¹ p est la profondeur de l'arbre
  moins celle de i. La complexitÃ© pour l'adjacence est de O(hd).

  Le nom du sommet i est un mot reprÃ©sentant le chemin P(i) =
  a_1,a_2,...,a_h allant de la racine (=0) Ã  i, chaque lettre a_t est
  le numÃ©ro du fils que l'on prend dans l'ordre pour atteindre
  i. C'est un numÃ©ro dans [0,d[. Si i=0, alors P(i)={} est vide.
*/
{
  int const d=Q->param[0]-1;
  int const h=Q->param[1];
  int const *f=Q->param+2;
  int i,s,x,y,n,p=h;

  switch(Q->code){

  case QUERY_INIT: /* calcule Q->n */
    if((h<0)||(d<0)) RET_error(6); // il faut h>=0 et d>=0
    for(i=0;i<d;i++) if(f[i]<=0) RET_error(6); // il faut f[i]>0
    RET_n(n_rectree(h,f,d));

  case QUERY_NAME: // nom du sommet i
    if(Q->i==0) strcpy(Q->name,"Îµ"); // pour la racine
    else{
      VIDE(Q->name); // name=""
      int t=0; // taille du nom de i
      y=(d>10); // vrai ssi il faut des virgules pour name
      x=Q->i; // x = copie du sommet i
      while(x>0){
	x--;
	for(i=s=0;i<d;i++){ // parcoure les fils
	  n=n_rectree(p-f[i],f,d); // taille du fils i
	  if((s<=x)&&(x<s+n)){ // x est dans le le fils i
	    t += sprintf(Q->name+t,"%i",i); // Ã©crit le numÃ©ro du fils i
	    if(t>NAMEMAX) Erreur(17); // nom trop long
	    x -= s;
	    p -= f[i]; // hauteur du fils i
	    if((y)&&(x>0)) t += sprintf(Q->name+t,","); // ajoute une ","
	    break; // sort du for(i=...)
	  }
	  s += n; // mise Ã  jour de la somme des tailles des fils
	}
      }
    }
    return 0;
    
  case QUERY_ADJ:;
    /* calcule le prÃ©fixe commun de P(i) et P(j) */
    x=Q->i, y=Q->j;/* copies de i et j, NB: x<y */
    
    for(;;){ // p = profondeur de l'arbre
      y--; // y est dans l'un des sous-arbres
      if(x==0){ // x est la racine
	for(i=s=0;i<d;i++){
	  if(y==s) RET_a(1); // y fils i de x
	  s += n_rectree(p-f[i],f,d); // taille du fils i
	}
	RET_a(0); // y n'est pas un fils de x
      }
      x--; // x est dans l'un des sous-arbres
      for(i=s=0;i<d;i++){ // trouver le sous-arbre contenant x
	n=n_rectree(p-f[i],f,d); // taille du fils i
	if((s<=x)&&(x<s+n)){ // x est dans le le fils i
	  x -= s, y -= s; // renumÃ©rote x,y Ã  partir de s
	  if(y>=n) RET_a(0); // y est dans un sous-arbre diffÃ©rent
	  // x et y sont dans le mÃªme sous-arbre, le fils i
	  p -= f[i]; // hauteur du fils i
	  break; // sort du for(i=...)
	}
	s += n; // mise Ã  jour de la somme des tailles des fils
      }
    }
  }

  return 1;
}


int rarytree(query* const Q)
/*
  Utilise Q->rep[0] pour le pÃ¨re. Les sommets sont numÃ©rotÃ©s selon un
  DFS modifiÃ©: on pose les fils avant la rÃ©cursivitÃ© (voir Dyck()).
*/
{
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:;
    int const n=Q->param[0]; if(n<=0) RET_n(0);
    int const z=Q->param[2]; if((z!=0)&&(z!=1)) RET_n(0); // z=0 ou 1
    int b=Q->param[1]; if(b<2) RET_n(0);
    int *B=Dyck(NULL,n,b-1,DYCK_KTREE); // B=arbre b-aire alÃ©atoire avec n noeuds internes
    SET_n(b*n+1+z);
    Q->k=1; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,1); // reprÃ©sentation implicite
    for(b=0;b<Q->n-z;b++) Q->rep[b][0]=B[b]; // copie l'arbre B dans Q->rep, |B|=N-z
    free(B);
    if(z) Q->rep[Q->n-1][0]=0; // le dernier sommet pointe vers la racine 0
    return 0;
  }
  
  return 1;
}


int ktree(query* const Q)
/*
  Q->rep[i][0...k[ sont les k pÃ¨res du sommet i. Cette fonction
  utilise un 3e paramÃ¨tre cachÃ©, Q->param[2] qui vaut 0, 1 ou 2
  suivant s'il faut gÃ©nÃ©rer un arbre alÃ©atoire (ktree), un chemin
  (kpath) ou un Ã©toile (kstar).
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:; /* calcule Q->n et Q->rep[i][0..k[ */
    int const n=Q->param[0];
    int const k=Q->param[1];
    if((k<0)||(n<=k)) RET_n(0);
    SET_n(n);

    int t,p,w,x,y,*T;
    if(Q->param[2]==0) T=Dyck(NULL,Q->n-k-1,1,DYCK_TREE); /* arbre de n-k noeuds */
    if(Q->param[2]==1) ALLOCZ(T,Q->n-k,_i-1); /* chemin de n-k noeuds */
    if(Q->param[2]==2) ALLOCZ(T,Q->n-k,0); /* une Ã©toile Ã  n-k noeuds */
    Q->k=k; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,k); // reprÃ©sentation implicite

    /* Chacun des k+1 sommets de la racine (numÃ©ros de 0 Ã  k) ont pour
       pÃ¨res tous les autres sommets (=> clique) */

    for(t=0;t<=k;t++) /* pour les k+1 sommets */
      for(p=0;p<k;p++) /* pour les k pÃ¨res 0..k-1 */
	Q->rep[t][p]=(t+p+1)%(k+1); /* il faut sauter t */

    /* On utilise le fait que les noeuds de T forment un DFS. En
       traitant les sommets dans l'ordre on est sÃ»r que le pÃ¨re est
       dÃ©jÃ  traitÃ© */
 
    for(t=k+1;t<Q->n;t++){ /* on dÃ©marre Ã  k+1, les k+1 sommets de la racine sont
			      dÃ©jÃ  traitÃ©s, tout sommet a donc une racine */
      p=T[t-k]; /* p=noeud pÃ¨re du sommet t du graphe */
      w=randomu(k+1); /* indice d'un des sommets du noeud pÃ¨re qui ne sera pas choisi */
      p += k; /* p=nom rÃ©el du pÃ¨re dans le graphe */
      Q->rep[t][0]=p; /* remplit avec le pÃ¨re lui-mÃªme */
      x=0; /* x=indice des pÃ¨res pour Q->rep[p], x=0..k-1 */
      y=(w>0); /* y=prochain indice des pÃ¨res pour Q->rep[t]. Si w=0 on
		  saute le pÃ¨re */
      while(x<k){
	Q->rep[t][y++]=Q->rep[p][x++];
	if(w==x) y--;
      }
    }
    
    free(T);  /* libÃ¨re l'arbre */
    return 0;
  }
  
  return 1;
}


int apollonian(query* const Q)
/*
  Utilise Q->rep pour la reprÃ©sentation implicite, Q->rep[i][0..2]
  sont les 3 pÃ¨res du sommet i. On pourrait factoriser avec polygon().
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
      
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT: /* calcule Q->n et Q->rep[i][0..2] */
    Q->n=Q->param[0]; if(Q->n<4) RET_n(0); /* il faut au moins 4 sommets */
    int const n=Q->n-3; /* n=nombre de sommets internes */
    int const m=3*n+1; /* nombre de sommets de l'arbre ternaire */
    int *P=Dyck(NULL,n,2,DYCK_KTREE); /* arbre ternaire Ã  n sommets internes */

    /*
      Principe de la construction. On part d'un arbre ternaire Ã  n
      sommets internes (dont la racine), comme dans rarytree(). La
      racine correspond Ã  un K_4 dont le centre est le sommet 3. Puis
      le k-iÃ¨me noeud interne de l'arbre (donc qui n'est pas une
      feuille) correspond Ã  un nouveau K_4 dont le centre est un
      nouveau sommet numÃ©rotÃ© k et qui est connectÃ© Ã  un triangle
      parent. Il y en a trois possibles suivant que le numÃ©ro du fils
      oÃ¹ l'on est.

                             0          3    (=triangle:012, centre:3)
                            /|\        /|\
                           1 2 3      4 . .  (=triangle:301, centre:4)
                          /|\        /|\
                         4 5 6      . 5 .    (=triangle:401, centre:5)
                          /|\        /|\
                         7 8 9      . . .

	P = [-,0,0,0,1,1,1,5,5,5] (= pÃ¨re dans l'arbre)
	C = [3,4,-,-,-,5,-,-,-,-] (= centre des triangles des sommets internes de l'arbre)
	
    */

    Q->k=3; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,3); /* reprÃ©sentation implicite */
    int u,p,c;

    /* calcule C[u]=numÃ©ro du centre du triangle ou 0 si feuille, pour
       tout noeud u de l'arbre */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[P[u]]=1; /* par dÃ©faut u est une feuille, et son pÃ¨re non */
    C[0]=c=3; /* racine = centre 3 */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numÃ©ro de centre */
    
    for(u=0;u<4;u++)   /* pour le premier tirangle */
      for(c=0;c<3;c++) /* pour les 3 fils de chaque sommet du K_4 */
	Q->rep[u][c]=(u+c+1)%4;

    /* on calcule le triangle Q->rep[c], pour chaque centre c=3...N.
       Q->rep[c][0..2] reprÃ©sente le triangle dont le centre est c */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[P[u]]; /* p=centre du pÃ¨re de u dans l'arbre */
	Q->rep[c][0]=p; /* un sommet du triangle est le centre */
	Q->rep[c][1]=Q->rep[p][u%3]; /* on en prend deux autres parmis le triangle du pÃ¨re */
	Q->rep[c][2]=Q->rep[p][(u+1)%3]; /* en fonction du numÃ©ro du fils -> u%3 */
      }
    }

    free(C);
    free(P);
    return 0;
  }
  
  return 1;
}


int polygon(query* const Q)
/*
  Utilise Q->rep pour la reprÃ©sentation implicite. Chaque sommet u
  possÃ¨de 2 pÃ¨res: Q->rep[u][0] et Q->rep[u][1].

  La construction est similaire Ã  apollonian(), avec un arbre binaire
  au lieu de ternaire, les K_3 remplaÃ§ant les K_4. On pourrait
  factoriser polygon() et apollonian() en passant en paramÃ¨tre un
  boolÃ©en, via Q->param[2] par exemple: 0 pour polygon() et 1 pour
  apollonian(). Pas clair qu'on puisse encore gÃ©nÃ©raliser cette
  construction oÃ¹ apollonian() et polygon() seraient des instances.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT: /* calcule Q->n et Q->rep[i][0..1] */
    Q->n=Q->param[0]; if(Q->n<3) RET_n(0);
    Q->k=2; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,2); /* chaque sommet a deux pÃ¨res */
    int u,c,p;
    int const n=Q->n-2; /* n=nombre de sommets internes de l'arbre binaire */
    int const m=2*n+1; /* m=nombre total de sommets de l'arbre */
    int *T=Dyck(NULL,n,1,DYCK_KTREE); /* calcule un arbre binaire alÃ©atoire T */

    /* Principe: on parcourt T selon un parcours en profondeur
       modifiÃ©, on pose les deux fils avant la rÃ©cursion.

       Dans l'exemple ci-dessous, la racine correspond au triangle
       (0,a,b). Puis, le fils gauche interne (=1) correspond au
       triangle (1,0,a). Le fils gauche interne 4 correspond au
       triangle (4,1,0). Les autres sommets feuilles ne correspondent
       Ã  aucun triangle. Dans le graphe, les numÃ©ros des sommets sont
       dÃ©calÃ©s de +2.

                             0           
			    / \            aâ”€â”€â”€b      0â”€â”€â”€1
                           1   2          / \ /      / \ /
                          / \            1â”€â”€â”€0      3â”€â”€â”€2
                         3   4            \ /        \ /
                            / \            4          4
                           5   6

       Plus prÃ©cisÃ©ment, Ã  chaque nouveau fils u>0 de T qui est un
       fils interne on associe un nouveau sommet c du graphe, un coin
       du triangle. Si u est un fils gauche de v=T[u] alors on
       connecte c aux sommets v et Q->rep[v][0] de G. Si u est un fils
       droit, on connecte c aux sommets v et Q->rep[v][1]. La paritÃ© d'un
       noeud dÃ©termine s'il s'agit d'un fils droit ou gauche.

	T = [-,0,0,1,1,4,4] (= pÃ¨re dans l'arbre)
	C = [2,3,-,-,4,-,-] (= derniers sommets des triangles des noeuds
                               internes de l'arbre, "-" si c'est une feuille)
    */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[T[u]]=1; /* par dÃ©faut u est une feuille, et son pÃ¨re non */
    C[0]=c=2; /* dernier sommet du triangle de la racine */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numÃ©ro de centre */
    
    Q->rep[0][0]=Q->rep[0][1]=0;    /* le sommet 0 n'a aucun pÃ¨re */
    Q->rep[1][0]=-1,Q->rep[1][1]=0; /* le sommet 1 a un seul pÃ¨re, 0 */
    Q->rep[2][0]=0,Q->rep[2][1]=1;  /* le sommet 2 a pour pÃ¨res 0 et 1 */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[T[u]]; /* p=centre du pÃ¨re de u dans l'arbre */
	Q->rep[c][0]=p;
	Q->rep[c][1]=Q->rep[p][u%2];
      }
    }

    free(C);
    free(T);
    return 0;
  }
  
  return 1;
}


int treep(query* const Q)
/*
  Utilise Q->wrap. Q->wrap[i] = pÃ¨re de i dans l'arbre.  Les feuilles
  de cet arbre sans sommet de degrÃ© deux sont les sommets de numÃ©ro <
  p.
*/
{
  int *P=Q->wrap; /* raccourci pour Q->wrap */
  switch(Q->code){

  case QUERY_END:
    if(Q->n>0) free(P);
    return 0;

  case QUERY_INIT:;
    int const p=Q->param[0]; /* p=#feuilles */
    if(p<3) RET_n(0);
    ALLOC(P,(p<<1)-1); /* wrap[] = 1 tableau d'au plus 2p-2 entiers */

    /*
      Principe du calcul d'un arbre alÃ©atoire Ã  p feuilles et sans
      sommets de degrÃ© deux.

      L'idÃ©e est de construire l'arbre Ã  partir des feuilles en
      dÃ©terminant le pÃ¨re de ces sommets, tout en garantissant qu'un
      sommet interne soit bien de degrÃ© au moins deux (sans son pÃ¨re).

      Soit A la liste des sommets qui n'ont encore pas de pÃ¨re. Cette
      liste de taille au plus p contient initialement toutes les
      feuilles: A={0,1,...,p-1}. Tant que |A|>2 on rÃ©pÃ¨te la procÃ©dure
      suivante:

      1. On tire un tableau alÃ©atoire R de taille |A| entiers >=0 dont
         la somme fait |A| et possÃ©dant au moins une valeur > 1. Pour
         cela on fixe une case alÃ©atoire Ã  2 puis on rÃ©pÃ¨te |A|-2 fois
         R[randomu(|A|)]++.

      2. Puis on dÃ©termine les pÃ¨res de certains sommets de A en
         fonction de R. L'idÃ©e est que si R[k]>1, alors les sommets
         A[u]...A[u+R[k]-1] vont recevoir le mÃªme pÃ¨re, un nouveau
         sommets crÃ©e et rÃ©injectÃ© Ã  la liste A. Plus prÃ©cisÃ©ment, on
         parcoure sÃ©quentiellement R. Si R[k]=0, alors on passe
         simplement Ã  k+1 sans rien faire d'autre. Si R[k]=1, on
         maintient A[u] dans la liste A on passe au prochain sommet
         u+1 de A. Le pÃ¨re du sommet A[u] n'est alors toujours pas
         fixÃ©. Par contre, si R[k]>1, alors on crÃ©e un nouveau sommet
         v que l'on ajoute Ã  la fin de la liste A, et les R[k] sommets
         A[u]...A[u+R[k]-1] ont alors tous pour pÃ¨re le sommet v.

      Lorsque |A|<=2 on ne rajoute plus aucun sommet, et le nombre de
      sommets N est alors dÃ©terminÃ©. On fixe alors que le sommet A[0]
      est racine. Si |A|=2 alors le pÃ¨re de A[1] est fixÃ© Ã  A[0]. Il
      n'est pas possible d'avoir |A|=0 car tout regroupement de sommet
      crÃ©e au moins un nouveau sommet dans pÃ¨re (et donc ajoutÃ© Ã  A).

      Une complexitÃ© de O(p^2) est possible car |A| diminue seulement
      d'une unitÃ© si Ã  chaque Ã©tape R ne contient qu'une seule valeur
      > 1. Cependant, en moyenne O(log(p)) Ã©tapes suffisent, soit
      O(p*log(p)) en tout, car il est facile de voir que R contient
      une fraction de |A| valeurs > 2. (Dans la reprÃ©sentation unaire
      de R il y a n/2 blocks et la moitiÃ© sont de longueur > 2.) Et
      pour chacun de tels blocks tous sauf 1 seront enlevÃ© de A.
    */

    int t;
    NALLOCZ(int,A,p,_i); /* tableau des sommets actifs */
    NALLOC(int,R,p); /* tableau des valeurs random */
    int u; /* u=indice dans A du prochain sommet actif */
    int q; /* q=indice dans A du prochain nouveau sommet actif, q<=u */
    int k; /* k=indice dans R, k>=u */
    int a; /* a=taille de A, a<=p */
    Q->n=a=p; /* Q->n=nb courant de sommets dans le graphe, Q->n>=p */

    while(a>2){
      for(k=0;k<a;R[k++]=0); /* tableau R Ã  zÃ©ro */
      R[randomu(a)]=2; /* met un "2" quelque part */
      for(k=2;k<a;k++) R[randomu(a)]++; /* incrÃ©mente a-2 fois des positions de R */
      for(k=u=q=0;k<a;k++){ /* parcoure les valeurs de R */
	if(R[k]==0) continue;
	if(R[k]==1) { A[q++]=A[u++]; continue; }
	t=u+R[k]; /* ici t>=2 */
	for(;u<t;u++) P[A[u]]=Q->n; /* P[A[u]]=pÃ¨re de A[u]=nouveau sommet */
	A[q++]=Q->n++; /* un sommet de plus, et un nouveau actif de plus */
      }
      a=q; /* nouvelle taille de A=nb de nouveaux sommets actifs */
    }

    P[A[0]]=-1;
    if(a==2) P[A[1]]=A[0];
    
    free(A);
    free(R);
    REALLOC(P,Q->n); /* recalibrage du tableau */
    Q->wrap=P; /* met Ã  jour Q->wrap */
    return 0;

  case QUERY_ADJ:
    RET_a(((P[Q->i]==Q->j)||(P[Q->j]==Q->i))); /* arbre */
    return 0;
  }

  return 1;
}


int halin(query* const Q)
/*
  Utilise treep().
*/
{
  int const p=Q->param[0]; /* p=#feuilles */
  if((Q->code==QUERY_ADJ)&&(Q->i<p)&&(Q->j<p)) /* cycle */
    RET_a(((Q->j==((Q->i+1)%p))||(Q->i==((Q->j+1)%p))));
  
  return treep(Q); /* arbre */
}


int permutation(query* const Q)
/*
  Utilise wrap, wrap[i] = permutation du sommet i.
*/
{
  switch(Q->code){

  case QUERY_END:
    if(Q->n>0) free(Q->wrap);
    return 0;
    
  case QUERY_NAME:
    sprintf(Q->name,"(%i,%i)",Q->i,Q->wrap[Q->i]);
    return 0;

  case QUERY_INIT: /* permutation alÃ©atoire de [0,n[ dans Q->wrap */
    SET_n(Q->param[0]);
    ALLOCZ(Q->wrap,Q->n,_i);
    Permute(Q->wrap,Q->n);
    return 0;
    
  case QUERY_ADJ:
    RET_a(((Q->i-Q->wrap[Q->i])*(Q->j-Q->wrap[Q->j])<0));
  }

  return 1;
}


int interval(query* const Q)
/*
  A chaque sommet i correspond un intervalle [a,b] de [0,2n[, avec
  a=Q->rep[i][0] et b=Q->rep[i][1].
*/
{
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    sprintf(Q->name,"[%i,%i]",Q->rep[Q->i][0],Q->rep[Q->i][1]);
    return 0;

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int const m=(Q->n<<1);
    int k,x;

    /* gÃ©nÃ¨re un intervalle Q->rep[k] pour k, [a,b] dans [0,2n[ avec a<=b */
    ALLOC2(Q->rep,Q->n,2);
    for(k=0;k<Q->n;k++){
      x=randomu(m);
      Q->rep[k][0]=x;
      Q->rep[k][1]=x+randomu(m-x);
    }
    return 0;

  case QUERY_ADJ:
    RET_a(
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][0])&&(Q->rep[Q->j][0]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][1])&&(Q->rep[Q->j][1]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->j][0]<=Q->rep[Q->i][0])&&(Q->rep[Q->i][1]<=Q->rep[Q->j][1]))
	  );
  }
  
  return 1;
}


int circle(query* const Q)
/*
  Graphe d'inclusion d'intervalle de [0,2n[, avec a=Q->rep[i][0] et
  b=Q->rep[i][1].
*/
{
  if(Q->code==QUERY_ADJ)
    RET_a(
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][0])&&(Q->rep[Q->j][1]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->j][0]<=Q->rep[Q->i][0])&&(Q->rep[Q->i][1]<=Q->rep[Q->j][1]))
	  );
  return interval(Q);
}


int sat(query* const Q)
/*
  Utilise: i<j.

  [0,2n[: les variables positives et nÃ©gatives
  [2n+t*k,2n+(t+1)*k[: t-Ã¨me clique de taille k

  Chaque sommet-clause t est connectÃ© Ã  une variable (positive ou
  nÃ©gative) alÃ©atoire stockÃ©e dans Q->rep[t][0].
*/
{
  int const n=Q->param[0];
  int const m=Q->param[1];
  int const k=Q->param[2];
  int const n2=(n<<1);
  int t; 

  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_INIT:
    if((n<1)||(m<1)||(k<1)) RET_n(0);
    SET_n(n2+m*k);
    ALLOC2(Q->rep,Q->n,1);    
    for(t=n2;t<Q->n;t++) Q->rep[t][0]=randomu(n2);
    return 0;
    
  case QUERY_ADJ:
    if(Q->j<n2) RET_a(((Q->j==Q->i+1)&&(Q->j&1))); /* i-j et j impaire */
    if(Q->i>=n2) RET_a((Q->j-Q->i<=k)); /* i et j dans la mÃªme clique ? */
    RET_a((Q->rep[Q->j][0]==Q->i)); /* i dans le matching et j dans une clique */
  }
  
  return 1;
}


int gpetersen(query* const Q)
/*
  Utilise i<j.
  u_i dans [0,n[ et v_i dans [n,2n[, Q->n=2n.
*/
{
  int const n=Q->param[0];
  int const r=Q->param[1];
  switch(Q->code){

    case QUERY_NAME:
      sprintf(Q->name,"%c%i",(Q->i<n)?'u':'v',Q->i%n);
      return 0;

    case QUERY_INIT:
      RET_n(2*n);
      
    case QUERY_ADJ:
      /* u_i-v_i */
      if(Q->j==Q->i+n) RET_a(1);

      /* sinon, cas frÃ©quent, pas d'arÃªte entre u_i et v_j */
      if((Q->i<n)&&(n<=Q->j)) RET_a(0);

      /* u_i-u_{i+1 mod n}, ici i<j<n */
      if(Q->i<n) RET_a((Q->j==(Q->i+1)%n)||(Q->i==(Q->j+1)%n));
      
      /* v_i-v_{i+r mod n} */
      /* ici i,j<n mais j<i possible*/
      RET_a((Q->j-n==(Q->i+r)%n)||(Q->i-n==(Q->j+r)%n));
  }

  return 1;
}


int antiprism(query* const Q)
{
  if(Q->code==QUERY_ADJ)
    if(Q->j==Q->param[0]+(Q->i+1)%Q->param[0]) RET_a(1);
  return gpetersen(Q);
}


int deltohedron(query* const Q)
/*
  Utilise i<j.
  Les sommets de [0,n[ forme le cycle avec n=Q->param[0].
*/
{
  int const n=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    if(n<=0) RET_n(0);
    RET_n(n+2);

  case QUERY_ADJ:
    RET_a(( (Q->j==(Q->i+1)%n) ||
	      ((Q->j==n-1)&&(Q->i==0)) ||
	      ((Q->j==n)&&((Q->i&1)==0)) ||
	      ((Q->j==n+1)&&((Q->i&1)==1)) )
	    );
  }

  return 1;
}


int helm(query* const Q)
/*
  Utilise i<j.
  On pourrait faire "wheel n -star -1".
*/
{
  int n=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    if(n<3) RET_n(0);
    RET_n(2*n+1);

  case QUERY_ADJ:
    if(Q->j-Q->i==n) RET_a(1); // branche
    if(Q->j>n) RET_a(0);
    if(Q->i==0) RET_a(1); // roue Ã  n rayons
    RET_a((Q->j-Q->i==1)||((Q->j==n)&&(Q->i==1)));
  }
  
  return 1;
}


int haar(query* const Q)
/*
  Utilise i<j.
  Utilise un paramÃ¨tre auxiliaire Q->param[1].
*/
{
  unsigned const n=Q->param[0];
  int const k=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    Q->param[1]=lg(n); // paramÃ¨tre auxiliaire
    RET_n(2*Q->param[1]); // n=0 si k=0

  case QUERY_ADJ:
    if((Q->i>=k)||(Q->j<k)) RET_a(0); // 0 si i,j dans la mÃªme part
    int const j=(Q->i-(Q->j-k)+k)%k; // u_i adjacent Ã  v_{i+j mod k} ?
    RET_a((n>>j)&1); // vrai ssi bit-j de n est Ã  1 ?
  }
  
  return 1;
}


int turan(query* const Q)
/*
  Utilise i<j.
*/
{
  int const n=Q->param[0];
  int const r=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    if((n<0)||(r<=0)||(n<r)) RET_error(6); /* il faut n>=r>0 */
    RET_n(n);

  case QUERY_ADJ:
    RET_a((Q->j-Q->i)%r);
  }
  
  return 1;
}


int kneser(query* const Q)
/*
  Q->rep[i][0..k[ sont les ensembles reprÃ©sentant les sommets.
*/
{
  int v,x,y;
  int const n=Q->param[0];
  int const k=Q->param[1];
  int const r=Q->param[2];
  
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],k,",","{}",1,"%i");
    return 0;

  case QUERY_INIT:
    if((n<0)||(k<0)||(k>n)) RET_n(0);
    SET_n(Binom(n,k)); /* on a besoin de Q->n pour allouer Q->rep */
    ALLOC2(Q->rep,Q->n,k); /* on connaÃ®t Q->n */
    if(Q->n==1) return 0; /* si Q->n=1, fin: graphe Ã  1 sommet */
    NALLOC(int,S,k);
    NextSet(S,-1,k); /* premier sous-ensemble */
    for(x=0;x<Q->n;x++){ /* pour tous les sommets x du graphe */
      for(y=0;y<k;y++) Q->rep[x][y]=S[y]; /* copie dans Q->rep[x] */
      NextSet(S,n,k); /* sous-ensemble suivant */
    }
    free(S);
    return 0;

  case QUERY_ADJ:
    /*
      Calcule si l'intersection possÃ¨de au plus r Ã©lÃ©ments.
      L'algorithme ici est en O(k) en utilisant le fait que les
      Ã©lÃ©ments de Q->rep sont rangÃ©s dans l'ordre croissant.
    */
    v=x=y=0; /* indices pour i et j, v=nb d'Ã©lements commun */

    while((x<k)&&(y<k)&&(v<=r))
      if(Q->rep[Q->i][x]==Q->rep[Q->j][y]) v++,x++,y++;
      else if(Q->rep[Q->i][x]<Q->rep[Q->j][y]) x++; else y++;

    RET_a((v<=r));
  }

  return 1;
}


int rig(query* const Q)
/*
  Q->rep[i][0] = taille t_i de l'ensemble associÃ© au sommet i
  Q->rep[i][1...t_i] = ensemble associÃ© au sommet i
*/
{
  int x,y,k,t;
  switch(Q->code){
    
  case QUERY_END:
    FREE2(Q->rep,Q->n); // car pas allouÃ© avec ALLOC2()
    return 0;
    
  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i]+1,Q->rep[Q->i][0],",","{}",1,"%i");
    return 0;
    
  case QUERY_INIT:
    SET_n(Q->param[0]);
    double const p=Q->dparam[0];
    k=Q->param[1];
    NALLOC(int,S,k+1); /* ensemble S[1...k] temporaire pour un sommet */
    ALLOC(Q->rep,Q->n);
    for(x=0;x<Q->n;x++){ /* pour chaque sommet x */
      t=0; for(y=1;y<=k;y++) if(RAND01<p) S[++t]=y; /* t=taille de S */
      ALLOC(Q->rep[x],t+1); /* espace pour le sommet x */
      Q->rep[x][0]=t; /* Ã©crit S dans Q->rep[x][1...t] */
      for(y=1;y<=t;y++) Q->rep[x][y]=S[y];
    }
    free(S);
    return 0;

  case QUERY_ADJ:
    /*
      DÃ©termine si l'intersection de Q->rep[i][1...] et
      Q->rep[j][1...]  est vide ou pas.  L'algorithme utilise le fait
      que les Ã©lÃ©ments de Q->rep sont rangÃ©s dans un ordre croissant.
    */

    x=y=1; /* indices pour les ensemble de i et j */
    k=Q->rep[Q->i][0]; /* taille de l'ensemble de i */
    t=Q->rep[Q->i][0]; /* taille de l'ensemble de j */

    while((x<=k)&&(y<=t)){
      if(Q->rep[Q->i][x]==Q->rep[Q->j][y]) RET_a(1);
      if(Q->rep[Q->i][x]<Q->rep[Q->j][y]) x++; else y++;
    }
    RET_a(0);
  }

  return 1;
}


int bdrg(query* const Q)
/*
  Utilise load().

  param[]={2t,n_1,d_1, ... n_t,d_t}, donc param[0]=|param|-1.  CrÃ©e
  dans Q->G un graphe dont la distribution des degrÃ©s est donnÃ©s par
  (n_i,d_i). La construction est basÃ©e sur un matching alÃ©atoire des
  demi-arÃªtes. Pour cela on construit un tableau T oÃ¹ chaque sommet u
  est rÃ©pÃ©tÃ© deg(u) fois dans T (cela forme les demi-arÃªtes sortantes
  de u). Les arÃªtes du graphe sont alors les arÃªtes simples entre les
  sommets T[2*i] et T[2*i+1] (on supprime boucle et arÃªte-multiple).
  La taille de T est Ã©gale Ã  âˆ‘d_i. Si elle n'est pas paire, on
  supprime une demi-arÃªte Ã  un sommet ayant un d_i>0. C'est toujours
  correct, car si d_i=0 pour tous les i, c'est que âˆ‘d_i est paire.

  ConsidÃ©rons l'exemple suivant:

    - 1 sommet  de degrÃ© 1 (le sommet 0)
    - 2 sommets de degrÃ© 2 (les sommets 1 et 2)
    - 1 sommets de degrÃ© 3 (le sommet 3)

    => param[]={6,1,1,2,2,1,3}
    => T[]=[ 0 1 1 2 2 3 3 3 ] (on rÃ©pÃ¨te le sommet u deg(u) fois)
    => permute(T) = [ 2 3 2 0 3 1 3 3 ]
    => arÃªtes:        2-3,2-0,3-1,3-3
    => arÃªtes simples: 2-3, 2-0, 3-1
*/
{
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return load(Q);

  case QUERY_INIT:;
    int c,u,u0,v,k,p;
    c=Q->param[0]; /* c=2*t=nombre de valeurs dans param[] sans le premier */
    if(c&1) RET_n(0); /* il faut un nombre paires de valeurs */
    for(Q->n=0,k=1;k<=c;k+=2) Q->n+=Q->param[k]; /* Q->n=nombre de sommets */
    TEST_n; // fin si graphe vide
    Q->G=new_graph(Q->n); /* crÃ©e le graphe (non vide) */
    for(k=u=0;k<c;k+=2){ /* alloue les listes */
      if(Q->param[k+2]>0) u0=u; /* mÃ©morise un sommet avec un degrÃ© d_i>0 */
      for(p=0;p<Q->param[k+1];p++,u++){
	Q->G->d[u]=Q->param[k+2]; /* degrÃ© du sommet u */
	ALLOC(Q->G->L[u],Q->G->d[u]); /* liste d'adjacence de u */
      }
    }
    for(v=k=0;k<c;k+=2) v+=Q->param[k+1]*Q->param[k+2]; /* v=âˆ‘ d_i*n_i=2*|E| */
    if(v&1){ /* si la somme est impaire, on enlÃ¨ve une demi-arÃªte au sommet u0 */
      Q->G->d[u0]--; /* u0 est nÃ©cessairement dÃ©fini si la somme Ã©tait impaire */
      REALLOC(Q->G->L[u0],Q->G->d[u0]); /* raccourcie la liste de u0 */
      v--; /* une demi-arÃªte de moins */
    }
    c=v; /* c=nombre de demi-arÃªtes */
    NALLOC(int,T,c); /* T=tableau des demi-arÃªtes */
    for(u=k=0;u<Q->n;u++) /* parcoure tous les arcs et remplit T */
      for(p=0;p<Q->G->d[u];p++) T[k++]=u;
    Permute(T,c); /* mÃ©lange alÃ©atoire les demi-arÃªtes */
    degres_zero(Q->G); /* pour ADD_EDGE */
    for(p=0;p<c;){ /* parcoure T et ajoute les arÃªtes simples */
      u=T[p++]; v=T[p++]; if(u==v) continue; /* pas de boucle */
      if(Q->G->d[u]>Q->G->d[v]) SWAP(u,v); /* pour une recherche plus rapide */
      if(SetSearch(v,Q->G->L[u],Q->G->d[u],0)<0) /* v dans Q->G->L[u] ? */
	ADD_EDGE(Q->G,u,v); /* non: ajoute u-v et met Ã  jour la taille des listes */
    }
    free(T);
    GraphRealloc(Q->G,Q->G->d);
    return 0;
  }

  return 1;
}

#undef SWAP
#define SWAP(x,y)						\
  do{								\
    (*(typeof(x)*)_SWAP)=(x);(x)=(y);(y)=(*(typeof(x)*)_SWAP);	\
  }while(0)

// Deux routines qui servent plusieurs fois dans fdrg().

// Ã‰change dans list[] les mini-sommets u1 et u2 tout en mettant Ã 
// jour pos[].
#define SWAP_MINI_SOMMET(u1,u2)		\
  SWAP2(list[pos[u1]],list[pos[u2]]);	\
  SWAP(pos[u1],pos[u2]);

#define SWAP2(x,y) do{ int z=(x);(x)=(y);(y)=z;}while(0)

// VÃ©rifie si les l'arÃªtes entre les sommets vi et vj formeraient une
// boucle ou une multi-arÃªte. Si c'est le cas, on met t=1, sinon t=0.
// Pour savoir si vi et vj sont dÃ©jÃ  voisins on cherche s'il existe un
// mini-sommet de vi qui a Ã©tÃ© pris (d'indice >= e) et qui a pour
// parent vj. Pour avoir un temps de recherche en
// O(min{deg(vi),deg(vj)}), cherche aussi dans les voisins de vj si
// deg[vj]<deg[vi].  NB: k^1 reprÃ©sente l'indice suivant k si k est
// pair, et prÃ©cÃ©dant si k est impair.
#define	NOT_GOOD_PAIR(vi,vj,t)				\
  if(vi==vj) t=1;					\
  else{							\
    if(deg[vj]<deg[vi]) u=vj,v=vi; else u=vi,v=vj;	\
    int d=deg[u],f=first[u],k;				\
    for(t=0;t<d;t++){				      	\
      k=pos[f+t];					\
      if((k>=e)&&(parent[list[k^1]]==v)) t=d;		\
    }							\
    t=(t>d);						\
  }


int fdrg(query* const Q)
/*
  Utilise load().

  param[] = {2t, n_1,d_1,..,n_t,d_t}.  CrÃ©e dans Q->G un graphe dont
  la distribution des degrÃ©s est donnÃ©s par param[]. La construction
  est basÃ©e sur l'algorithme (avec rejet) de [BKS10] (section 6.1)
  basÃ© sur l'implÃ©mentation de [SW97]. Une prÃ©-version a Ã©tÃ©
  programmÃ©e par AliÃ©nor Brabant, en stage de L2 en juin 2016.

  L'algorithme est le suivant (ici deg(i)=degrÃ© du sommet i souhaitÃ©,
  D[i]=degrÃ© restant du sommet i, M=sum_i deg(i)=2 fois le nombre
  total d'arÃªtes souhaitÃ©es):

  (1) E={}, D[i]=deg(i) pour i=0..n-1, V={0,...,n-1}

  (2) Choisir i,j de V avec proba p_ij = D[i]*D[j] *
      (1-deg(i)*deg(j)/2M) parmi toutes paires avec (i,j) tq i<>j et
      {i,j} pas dans E. Ajouter {i,j} Ã  E puis dÃ©crÃ©menter D[i] et
      D[j].

  (3) RÃ©pÃ©ter (2) jusqu'Ã  ce qu'on ne puisse plus rajouter d'arÃªte Ã  E

  (4) Si |E|<M/2, recommencer sinon G=(V,E).

  L'algorithme termine Ã  l'Ã©tape (4) la premiÃ¨re fois avec probabilitÃ©
  1-o(1). Pour les graphes rÃ©guliers, on peut prendre plus simplement
  p_ij = D[i]*D[j]. Pour gÃ©nÃ©rer un biparti alÃ©atoire (de sÃ©quence
  fixÃ©e) on peut remplacer (2) par p_ij = D[i]*D[j] *
  (1-deg(i)*deg(j)/M) et ajouter seulement les arÃªtes qui ne sont pas
  dans une mÃªme part.

  On implÃ©mente l'algorithme en gÃ©rant une liste de "mini-sommets",
  qui sont les brins incidents aux sommets et qu'on doit connectÃ©s.
  Chaque sommet i possÃ¨de deg(i) mini-sommets. On choisit donc deux
  mini-sommets avec la bonne probabilitÃ© p_ij tels que si connectÃ©s,
  ils ne forment ni de boucles ni de multi-arÃªtes. Si c'est le cas, on
  sÃ©lectionne cette paire, pour former une nouvelle arÃªte, et l'on
  dÃ©place ces deux mini-sommets dans la 2e partie de la liste. Ã€ la
  fin, la liste reprÃ©sente toutes les paires de mini-sommets
  connectÃ©s, et donc la liste des arÃªtes.

  Pour obtenir une complexitÃ© en O(M*dmax+dmax^4) il faut dÃ©couper
  l'Ã©tape (2) en 3 phases. Il faut aussi faire attention Ã  l'Ã©tape
  (3).  Certes, on ne peut ajouter que O(M) arÃªtes, mais dÃ©tecter si
  l'on peut ou non ajouter encore une arÃªte n'est pas Ã©vident. En
  effet, il peut rester des paires de mini-sommets non connectÃ©s sans
  pourtant pouvoir ajouter une seule arÃªtes si toutes les paires de
  mini-sommets ont des sommets dÃ©jÃ  connectÃ©s. Voici un exemple de
  blocage pour n=5 sommets et une sÃ©quence de degrÃ© 2,2,2,3,3.

                         sommets={0,1,...,4}
                    mini-sommets={a,b,...,l}
  
     0  1  2  3   4                            0  1  2  3   4
    /â”‚ /â”‚ /â”‚ /â”‚\ /â”‚\                          /â”‚ /â”‚ /â”‚ /â”‚\ /â”‚\  
    ab cd ef ghi jkl                          ab cd ef ghi jkl 
    â”‚â”‚ â”‚â”‚ â”‚â”‚ â”‚â”‚* â”‚â”‚*                          â”‚â”‚ â”‚â”‚ â”‚â”‚ â”‚â”‚â”‚ â”‚â”‚â”‚
    â”‚â””â”€â”˜â”‚ â”‚â”‚ â”‚â””â”€â”€â”˜â”‚                           â”‚â””â”€â”˜â””â”€â”¼â”¼â”€â”˜â”‚â””â”€â”˜â”‚â”‚
    â””â”€â”€â”€â”¼â”€â”˜â”‚ â”‚    â”‚                           â””â”€â”€â”€â”€â”€â”¼â”¼â”€â”€â”˜   â”‚â”‚
        â””â”€â”€â”¼â”€â”˜    â”‚                                 â””â”¼â”€â”€â”€â”€â”€â”€â”¼â”˜
           â””â”€â”€â”€â”€â”€â”€â”˜                                  â””â”€â”€â”€â”€â”€â”€â”˜

  choix: ae,bc,dg,fk,hj                    choix: ah,bc,dg,el,fk,ij
  dernier choix impossible

  Cependant, pour Ãªtre bloquÃ© il doit exister dans le graphe courant
  une clique de k sommets ayant chacun encore au moins un mini-sommet
  libre. Une clique de k sommets ayant encore un mini-sommet de libre
  implique k <= dmax: un sommet i de la clique doit Ãªtre connectÃ© Ã 
  k-1 autre sommet et ce degrÃ© est <= dmax-1 puisque un mini-sommet
  est encore libre. Autrement dit, si le nombre d'arÃªtes restantes M -
  âˆ‘ D[i] = M-e > dmax^2 on ne peut Ãªtre bloquÃ© puisqu'il reste plus de
  dmax sommets avec au plus dmax-1 mini-sommets de libres.
  
*/
{
  switch(Q->code){

  case QUERY_ADJ:
    return load(Q);

  case QUERY_INIT:;
    int c=Q->param[0]; // c=2*t=nombre de valeurs dans param[] sans le premier
    if((c&1)||(c<0)) RET_error(6); // il faut un nombre paires de valeurs
    Q->n=graphical(Q->param+1,c/2); // Q->n=nombre de sommets
    if(Q->n<0) RET_error(34); // erreur: sÃ©quence non graphique
    TEST_n; // fin si graphe vide (si n=0)
    int u,v,i,j,k,t,d,vi,vj,ui,uj,e,p,nv,nh;

    /* calcule le nombre M de mini-sommets, pour les ALLOC() suivants */
    for(t=0,i=1;i<=c;i+=2) t+=Q->param[i]*Q->param[i+1];
    int const M=t; // M=nombre total de mini-sommets=âˆ‘d_i

    /* Conventions:
         i = indice ou position, entier de [0,M[
         u = mini-sommet, entier de [0,M[
         v = sommet, entier de [0,Q->n[
    */
    
    NALLOC(int,list,M); // list[i]=u, liste des mini-sommets (au dÃ©part list[i]=i)
    NALLOC(int,pos,M); // pos[u]=i, position i du mini-sommets u dans list[]
    NALLOC(int,parent,M); // parent[u]=v, sommets correspondant au mini-sommets u
    NALLOC(int,listv,Q->n); // listv[i]=v, liste des sommets libres
    NALLOC(int,dfree,Q->n); // dfree[v]=nombre de mini-sommets libre pour le sommet v
    NALLOC(int,first,Q->n); // first[v]=1er mini-sommet du sommet v;
   // les mini-sommets de v sont donc first[v], first[v]+1,..., first[v]+deg[v]-1

    /* On parcoure les (n_i,d_i) et calcule:
       - le degrÃ© max
       - la taille des listes d'adjacence du graphe final Q->G
       - le tableau first[]
       - le tableau parent[]
    */

    Q->G=new_graph(Q->n); // crÃ©e le graphe, alloue Q->G->d
    t=u=v=0; // t=degrÃ© max, u=mini-sommet, v=sommet

    for(i=1;i<=c;i+=2){
      d=Q->param[i+1]; // degrÃ© courant
      t=max(t,d); // t=max{d_i}
      for(j=0;j<Q->param[i];j++,v++){
	first[v]=u; // premier mini-sommet du sommet v
	for(k=0;k<d;k++) parent[u++]=v; // tous le mÃªme parent v
	Q->G->d[v]=d; // degrÃ© du sommet v
	ALLOC(Q->G->L[v],d); // liste d'adjacence de u
      }
    }

    int const M2=2*M; // constante pour les probas
    int const dmax=t; // dmax=degrÃ© max du graphe final
    int const ddmax=t*t; // constante dmax^2
    int const pair=(~1); // mask permettant (avec &) de rendre pair un entier 
    int const noreg=(c>2); // noreg=vrai ssi le graphe n'est pas rÃ©gulier
    int * const deg=Q->G->d; // raccourci: deg[v]=degrÃ© de v dans le graphe final
    NALLOC(int,listh,4*ddmax); // listh[]=liste d'arÃªtes possible pour la PHASE 3

    for(;;){ // forever ... sauf si la PHASE 3 rÃ©ussit. En principe
	     // elle rÃ©ussit dÃ¨s le 1er coup avec proba 1-o(1)

      e=M; // e=nombre de mini-sommets encore disponibles

      // PHASE 1
      //
      // On titre uniformÃ©ment deux mini-sommets pris dans list[0..e[,
      // e Ã©tant le nombre de mini-sommets disponibles (non
      // connectÃ©s). La probabilitÃ© d'obtenir la paire de sommets
      // (vi,vj) est ainsi proportionnelle Ã  D[vi]*D[vj]. On garde
      // cette paire avec proba 1-deg(vi)*deg(vj)/2M (ou proba 1 si
      // c'est un graphe rÃ©gulier). Ainsi, la proba de succÃ¨s (dans le
      // cas non-rÃ©gulier) est proprotionnelle au produit (D[i]*D[j])
      // * (1-deg(vi)*deg(vj)/2M). On rejette si c'est une boucle ou
      // si vi et vj Ã©taient dÃ©jÃ  voisins, ce qui prend un temps de
      // O(dmax). Le nombre de rÃ©pÃ©titions est 2 en moyenne. La
      // complexitÃ© totale de la PHASE 1 est O(m*dmax) en moyenne.
      
      // initialise list[] et pos[]
      for(i=0;i<M;i++) list[i]=pos[i]=i;

      while(e>2*ddmax){ // on ne peut pas Ãªtre bloquÃ© avec ce test

	// choisit deux mini-sommets: list[i] et list[j]
	i=randomu(e), vi=parent[list[i]];
	j=randomu(e), vj=parent[list[j]];
	if(noreg) // NB: si le graphe est rÃ©gulier, on accÃ¨pte toujours vi,vj
	  if(M2*RAND01>M2-deg[vi]*deg[vj]) continue; // rejet ?
	NOT_GOOD_PAIR(vi,vj,t); // vi=vj ou si vi-vj ?
	if(t) continue; // recommence si pas bon
	
	// on ajoute l'arÃªte vi-vj en dÃ©plaÃ§ant les mini-sommets
	// list[i] et list[j] en position e-1 et e-2 dans list[] tout en
	// mettant Ã  jour pos[] et e.
	ui=list[i],uj=list[j];
	k=list[--e]; SWAP_MINI_SOMMET(ui,k); // Ã©change ui et k
	k=list[--e]; SWAP_MINI_SOMMET(uj,k); // Ã©change uj et k

	// NB: AprÃ¨s le premier SWAP_MINI_SOMMET(), list[] est
	// modifiÃ©e, et donc il est possible que le mini-sommet uj ne
	// soit plus en position j dans list[] (si j=e-1 par exemple).
	// Il faut donc fixer ui=list[i] et uj=list[j] AVANT de faire
	// le premier SWAP_MINI_SOMMET().
      }
    
      // PHASE 2
      //
      // Cette phase est similaire Ã  la PHASE 1, sauf qu'on considÃ¨re
      // plutÃ´t les sommets libres que les mini-sommets libres. On
      // reste dans cette phase tant que le nombre de sommets encore
      // libres (nv) est assez grand: nv > 2dmax. On choisit deux
      // sommets libres, et on rÃ©pÃ¨te jusqu'Ã  en trouver deux qui ne
      // soient pas dÃ©jÃ  voisins (2 rÃ©pÃ©titions en moyenne). On ne
      // peut pas Ãªtre bloquÃ© puisqu'on ne peut pas avoir de clique si
      // le nombre de sommets restant est > dmax. Puis on choisit
      // uniformÃ©ment un mini-sommet pour chacun des sommets, et on
      // rÃ©pÃ¨te jusqu'Ã  en avoir deux non connectÃ©s (O(dmax^2)
      // rÃ©pÃ©titions en moyenne). La complexitÃ© totale de la PHASE 2
      // est O(dmax^4) en moyenne.
      
      // initialise (1) listv[0..nv[, la liste des sommets libres,
      // c'est-Ã -dire ayant encore au moins un mini-sommet de libres;
      // (2) dfree[v], le nombre de mini-sommets de libres pour v si v
      // est dans listv[] (non dÃ©fini sinon).
      for(v=nv=0;v<Q->n;v++){
	u=first[v]; // u=premier mini-sommet de v
	d=k=deg[v]; // d=dfree[v]=deg[v] au dÃ©part
	for(i=0;i<k;i++) if(pos[u+i]>=e) d--;
	if(d) dfree[v]=d,listv[nv++]=v; // on ajoute v Ã  listv[] si d>0
      }

      while(nv>2*dmax){ // ne peut pas Ãªtre bloquÃ© si nv>dmax
	
	// choisit deux sommets de listv[] non voisins
	i=randomu(nv), vi=listv[i];
	j=randomu(nv), vj=listv[j];
	NOT_GOOD_PAIR(vi,vj,t); // vi=vj ou si vi-vj ?
	if(t) continue; // recommence si pas bon
	
	// met Ã  jour dfree[], puis Ã©ventuellement listv[] et nv
	dfree[vi]--; if(dfree[vi]==0) listv[i]=listv[--nv]; // dÃ©place en supprimant vi
	if(j==nv) j=i; // NB: dans ce cas vj a pu Ãªtre dÃ©placÃ© en i
	dfree[vj]--; if(dfree[vj]==0) listv[j]=listv[--nv]; // dÃ©place en supprimant vj

	// choisit un mini-sommet pour vi et pour vj et vÃ©rifie qu'ils
	// ne soient pas dÃ©jÃ  connectÃ©s; on ne peut pas Ãªtre bloquÃ©
	// car vi et vj ne sont pas connectÃ©s et ils ont chacun un
	// mini-sommet de libre. Le nombre de rÃ©pÃ©titions est
	// O(dmax^2) en moyenne.
	do{
	  ui=first[vi]+randomu(deg[vi]);
	  uj=first[vj]+randomu(deg[vj]);
	}while((pos[ui]>=e)||(pos[uj]>=e));
	
	// on ajoute l'arÃªte entre les mini-sommets ui et uj, et met Ã 
	// jour list[], pos[] et e
	k=list[--e]; SWAP_MINI_SOMMET(ui,k); // Ã©change ui et k
	k=list[--e]; SWAP_MINI_SOMMET(uj,k); // Ã©change uj et k
      }
    
      // PHASE 3
      //
      // On considÃ¨re le graphe H des arÃªtes encore possibles, donc
      // conctant les sommets de listv[]. Il possÃ¨de nv <= 2dmax
      // sommets et <= nv*(nv-1)/2 < 2dmax^2 arÃªtes. On choisit une
      // arÃªte alÃ©atoire uniforme de H qu'on accepte ensuite avec une
      // probabilitÃ© D[i]*D[j]/dmax^2. Le nombre de rÃ©pÃ©titions est
      // donc O(dmax^2) en moyenne. On met Ã  jour H et on rÃ©pÃ¨te
      // jusqu'Ã  ce que H n'ait plus de sommets libres.  Pour mettre Ã 
      // jour H, il faut supprimer l'arÃªte i-j choisie mais aussi
      // supprimer toutes les arÃªtes de H oÃ¹ l'un des sommets i ou j
      // apparaÃ®t si jamais D[i] ou D[j] passe Ã  0. Si, Ã  l'issue de
      // la PHASE 3 le nombre d'arÃªtes ajoutÃ©es au graphe final est <
      // M/2 (H devient vide trop tÃ´t), il faut alors recommencer les
      // trois phases. La complexitÃ© totale de la PHASE 3 est
      // O(dmax^4) en moyenne.
      
      // On considÃ¨re listh[0..nh[, la liste des arÃªtes possibles de H
      // (graphe qu'on ne construit pas en fait) oÃ¹ l'on stocke chaque
      // arÃªte possible dans des cases consÃ©cutives de listh[]. La
      // taille de ce tableau doit Ãªtre nv*(nv-1) < 4dmax^2. On
      // l'initialise en balayant toutes les paires possibles de
      // sommets de listv[] et en testant leur adjacence.  Aussi, on
      // va stocker dans list[0..p[ (au dÃ©but donc) les arÃªtes
      // gÃ©nÃ©rÃ©es Ã  la PHASE 3. On a la place car celles gÃ©nÃ©rÃ©es aux
      // phases prÃ©cÃ©dantes ont Ã©tÃ© stockÃ©es dans list[e..M[ (Ã  la fin
      // donc) et p<=e.

      for(i=nh=0;i<nv;i++)
	for(j=i+1;j<nv;j++){ // NB: i<j
	  vi=listv[i],vj=listv[j]; // vi,vj=sommets encore libre
	  NOT_GOOD_PAIR(vi,vj,t); // t=(vi=vj ou vi-vj)
	  if(t==0) listh[nh++]=vi,listh[nh++]=vj; // ajoute vi-vj
	}
      
      p=0; // p=position de la prochaine entrÃ©e libre dans list[]
      while(nh>0){ // tant qu'il reste une arÃªte dans listh[]
      
	// choisit une arÃªte i-j de listh[] uniformÃ©ment
	t=(randomu(nh)&pair); // position (paire) alÃ©atoire
	vi=listh[t],vj=listh[t+1];
	if(RAND01*ddmax>dfree[vi]*dfree[vj]) continue; // rejet ?

	// ajoute l'arÃªte vi-vj dans list[] et l'enlÃ¨ve de listh[]
	list[p++]=vi,list[p++]=vj;
	listh[t+1]=listh[--nh],listh[t]=listh[--nh];

	// met Ã  jour dfree[]
	i=(--dfree[vi]==0); // i=vrai ssi il faut supprimer vi dans listh[]
	j=(--dfree[vj]==0); // j=vrai ssi il faut supprimer vj dans listh[]

	if(i||j){
	  // ici il faut supprimer vi ou vj (ou les deux) dans tout
	  // listh[]. On parcoure listh[] une 1Ã¨re fois et on met des
	  // paires (-1,-1) si vi ou vj apparaÃ®t, puis une 2e fois
	  // pour enlever les (-1,-1) en dÃ©laÃ§ant la derniÃ¨re de
	  // listh[]. NB: si v est en position t dans listh[], alors
	  // son voisin est en position t^1. Cela prend un temps
	  // |listh|=O(dmax^2), le mÃªme temps pour accepter vi-vj en
	  // moyenne.
	  for(t=0;t<nh;t++)
	    if((i&&(listh[t]==vi))||(j&&(listh[t]==vj)))
	      listh[t]=listh[t^1]=-1; // efface l'arÃªte contenant vi ou vj
	  t=0;
	  while(t<nh)
	    if(listh[t]<0) listh[t+1]=listh[--nh],listh[t]=listh[--nh];
	    else t+=2;
	}
      }
      
      if(p==e) break; // fin: list[] est complÃ¨tement remplie
    }
    
    // on Ã©crit les arÃªtes de list[] dans Q->G, le codage des arÃªtes
    // n'est pas le mÃªme avant et aprÃ¨s la position e
    degres_zero(Q->G); // pour ADD_EDGE(Q->G,...)
    for(i=0;i<e;i+=2) ADD_EDGE(Q->G,list[i],list[i+1]);
    for(;i<M;i+=2) ADD_EDGE(Q->G,parent[list[i]],parent[list[i+1]]);
    
    free(list);
    free(pos);
    free(parent);
    free(first);
    free(listv);
    free(listh);
    free(dfree);

    return 0;
  }

  return 1;
}
#undef SWAP_MINI_SOMMET
#undef NOT_GOOD_PAIR


int seg_intersection(point p1,point q1,point p2,point q2){
/*
  Renvoie vrai ssi le segment ]p1q1[ intersecte le segment
  ]p2q2[. Pour s'intersecter il faut que: (1) les points p2,q2 soient
  de part et d'autre de la droite portÃ©e par (p1q1); et (2) les points
  p1,q1 soient de part et d'autre de la droite portÃ©e par
  [p2q2]. Sinon, si les points sont colinÃ©aires et il faut tester si
  l'un des points appartient au segment de l'autre.

  Pour dÃ©terminer si point X est au-dessus, en-dessous ou sur la
  droite portÃ©e par (AB) il suffit de calculer le signe de
  det(X-A,B-A).
*/
  double const dp2=det(q1.x-p1.x,q1.y-p1.y,p2.x-p1.x,p2.y-p1.y); // p2 au-dessus de (p1q1) ?
  double const dq2=det(q1.x-p1.x,q1.y-p1.y,q2.x-p1.x,q2.y-p1.y); // q2 au-dessus de (p1q1) ? 
  double const dp1=det(q2.x-p2.x,q2.y-p2.y,p1.x-p2.x,p1.y-p2.y); // p1 au-dessus de (p2q2) ?
  double const dq1=det(q2.x-p2.x,q2.y-p2.y,q1.x-p2.x,q1.y-p2.y); // q1 au-dessus de (p2q2) ? 
  
  if((dp2*dq2<0)&&(dp1*dq1<0)) return 1; // cas gÃ©nÃ©ral

// est-ce que C est dans ]AB[ ?
#define IS_IN(C,A,B)				\
  ( (min(A.x,B.x)<C.x)&&(C.x<max(A.x,B.x))&&	\
    (min(A.y,B.y)<C.y)&&(C.y<max(A.y,B.y)) )
 
  if((dp2==0)&&(IS_IN(p2,p1,q1))) return 1; // p2 dans ]p1q1[
  if((dq2==0)&&(IS_IN(q2,p1,q1))) return 1; // q2 dans ]p1q1[
  if((dp1==0)&&(IS_IN(p1,p2,q2))) return 1; // p1 dans ]p2q2[
  if((dq1==0)&&(IS_IN(q1,p2,q2))) return 1; // q1 dans ]p2q2[
  return 0;
}
 

int rlt(query* const Q)
/*
  Utilise load().

  Pour simplifier, on prend V = [0,2p[ x [0,2q[ comme ensemble de
  points de la grille, c'est-Ã -dire les couples (i,j) avec i et j
  pairs. On en dÃ©duit que l'ensemble M des points milieu est [0,2p-1[
  x [0,2q-1[ \ V, c'est-Ã -dire des points qui ont au moins une
  coordonnÃ©es impaires.

  Pour chaque point milieu R, il faut prendre tous les points A
  possibles de la grille et vÃ©rifier que le point B symÃ©trique de A
  selon R est dans la grille, que l'arÃªte rÃ©sultante (le segment [AB])
  soit de la bonne longueur (<=d) et n'intersectent aucune autre arÃªte
  dÃ©jÃ  en place dans E.
*/
{
  switch(Q->code){

  case QUERY_ADJ:
    return load(Q);

  case QUERY_END:
    return free_pos(Q);

  case QUERY_INIT:;
    int const p=Q->param[0];
    int const q=Q->param[1];
    if((p<=0)||(q<=0)) RET_n(0);
    Q->n=p*q;

    int d=Q->param[2];
    if(d<0) d=2*Q->n; // d=+âˆž
    d=2*d-1; // d=1->1, d=2->3, d=3->5, etc.
    
    int const m=(2*p-1)*(2*q-1)-p*q; // m=nombre d'arÃªtes=nombre de points dans M
    int const mod2=(~1); // mask (avec &) pour entier pair infÃ©rieur

    // t=nombre de points maximum que l'on peut mettre dans L
    int t=(p&1); t=(p==1)? 1 : p*(q-t)/2+t;

    NALLOC(int,M,2*m); // M=liste de points milieu de chaque arÃªte
    NALLOC(int,E,4*m); // E=liste des segments [AB] qui sont des arÃªtes
    NALLOC(int,L,2*t); // L=liste des points A possibles pour un point R de M
    Q->G=new_graph(Q->n); // NB: les degrÃ©s sont nuls

    int i,j,k,u,v;
    int Ri,Rj,i0,i1,j0,j1,di,dj;
    point A,B,E0,E1;
    int nM=0;
    int nE=0;
    int nL;

    // initialise les points de M
    if(d) // si d=0, il faut nM=0
      for(i=0;i<2*p-1;i++)
	for(j=0;j<2*q-1;j++)
	  if((i&1)||(j&1)) M[nM++]=i,M[nM++]=j; // au moins une coord. impaire

    while(nM){ // tant que M n'est pas vide (et d<>0)
      t=randomu(nM)&mod2; // position (paire) alÃ©atoire
      Ri=M[t],Rj=M[t+1]; // R=(Ri,Rj) points alÃ©atoire de M
      M[t+1]=M[--nM],M[t]=M[--nM]; // supprime R de M

      // dÃ©termine la liste L des points A=(i,j) possibles pour R. Il
      // faut R milieu de [AB] donc B=2R-A, A et B dans [0,2p-1[ x
      // [0,2q-1[, et aussi que ]AB[ ne soit pas trop grand et ne
      // contienne aucun point entier Ã  part R. Pour cela, il suffit
      // que les coordonnÃ©es de (R-A) soient premiÃ¨res entre elles. A
      // cause des symÃ©tries, seuls certains points A doivent Ãªtre
      // testÃ©s dans le rectangle [i0,i1]x[j0,j1] comme suit:
      if(Rj<q) j0=0,j1=Rj&mod2; else j0=(Rj+1)&mod2,j1=2*q-2;
      if(Ri<p) i0=0,i1=2*Ri; else i0=2*Ri-2*p+2,i1=2*p-2;

      nL=0;
      for(i=i0;i<=i1;i+=2){
	A.y=i;B.y=2*Ri-i;di=abs(i-Ri);
	for(j=j0;j<=j1;j+=2){
	  dj=abs(j-Rj);
	  if((j==Rj)&&(i>Ri)) continue; // AR mÃªme colonne et A au-dessus
	  if((j==Rj)&&(di>1)) continue; // AR mÃªme colonne et A,R pas voisins
	  if((i==Ri)&&(dj>1)) continue; // AR mÃªme ligne et A,R pas voisins
	  if(max(di,dj)>d) continue; // [AR] est trop long
	  if(pgcd(di,dj)!=1) continue; // ]AR[ contient un point entier
	  A.x=j;B.x=2*Rj-j;
	  // est-ce que ]AB[ intersecte une arÃªte de E ?
	  for(k=0;k<nE;k+=4){ // pour toutes les arÃªtes e=(E0,E1) de E
	    E0.y=E[k+0],E0.x=E[k+1];
	    E1.y=E[k+2],E1.x=E[k+3];
	    if(seg_intersection(A,B,E0,E1)) break; // si intersection
	  }
	  if(k<nE) continue; // A n'est pas bon: suivant
	  L[nL++]=i,L[nL++]=j; // ajoute A=(i,j) Ã  L
	}
      }

      // choisit un point A de L, met Ã  jour le degrÃ© des sommets A et
      // B, et met l'arÃªte correspondante dans E
      t=randomu(nL)&mod2; // position (paire) alÃ©atoire
      u=(   L[t]/2)*q + (   L[t+1]/2); // u=numÃ©ro du sommet de A
      v=(Ri-L[t]/2)*q + (Rj-L[t+1]/2); // v=numÃ©ro du sommet de B
      Q->G->d[u]++;
      Q->G->d[v]++;
      E[nE++]=L[t],E[nE++]=L[t+1]; // A
      E[nE++]=2*Ri-L[t],E[nE++]=2*Rj-L[t+1]; // B
    }

    // copie les arÃªtes de E dans Q->G
    for(u=0;u<Q->G->n;u++) ALLOC(Q->G->L[u],Q->G->d[u]); // alloue les listes
    degres_zero(Q->G); // pour ADD_EDGE()
    for(k=0;k<nE;k+=4){
      u=(E[k+0]/2)*q + (E[k+1]/2); // u=numÃ©ro du sommet de A
      v=(E[k+2]/2)*q + (E[k+3]/2); // v=numÃ©ro du sommet de B
      ADD_EDGE(Q->G,u,v); // ajoute l'arÃªte u-v
    }

    // fin
    free(M);
    free(E);
    free(L);
    return 0;
  }

  return 1;
}


int kout(query* const Q)
/*
  Modifie Q->i et Q->j.
  Marche aussi en orientÃ©.

  Q->rep[i]=tableau des voisins de i. Si Q->rep[i][j]<0, alors c'est
  la fin prÃ©maturÃ©e de la liste, par dÃ©faut on ne lit qu'au plus les k
  premiÃ¨res cases. Attention ! Si i<j, alors Q->rep[i] ne peut pas
  contenir le sommet j (qui a Ã©tÃ© crÃ©e aprÃ¨s i).
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:;
    int const n=Q->param[0];
    int const k=Q->param[1];
    if((k<1)||(k>=n)) Erreur(6); // paramÃ¨tre incorrect
    SET_n(n);
    Q->k=-k; // pour adjacency_rep()
    ALLOC2(Q->rep,n,k);
    NALLOC(int,T,n);
    int r,d,z,x,y;

    Q->rep[0][0]=-1;  // le sommet 0 est seul !
    for(x=1;x<n;x++){ // x=prochain sommet Ã  rajouter
      r=min(x,k);     // le degrÃ© de x sera au plus r=min{x,k}>0
      d=1+randomu(r); // choisir un degrÃ© d pour x: dâˆˆ[1,...,r]
      for(y=0;y<x;y++) T[y]=y; // tableau des voisins possibles
      r=x;              // r-1=index du dernier Ã©lÃ©ment de T
      for(y=0;y<d;y++){ // choisir d voisins dans T
	z=randomu(r);   // on choisit un voisin parmi ceux qui restent
        Q->rep[x][y]=T[z]; // le y-Ã¨me voisin de x est T[z]
        T[z]=T[--r];       // supprime T[z] et met le dernier Ã©lÃ©ment de T Ã  sa place
      }
      if(d<k) Q->rep[x][d]=-1; // arrÃªte la liste des voisins de x
    }

    free(T);
    return 0;

  }
  
  return 1;
}


int expander(query* const Q)
/*
  Q->rep[i]=tableau des k>0 successeurs de i. Il est possible que le mÃªme
  voisin apparaisse plusieurs fois dans Q->rep[i].

  Algorithme:
   1. on part du tableau T[i]=i, pour i=0..n-1
   2. on forme le cycle T[0]-T[1]-...-T[n-1]-T[0]
   3. on permute seulement les n-1 premiÃ¨res valeurs de T
   4. on recommence en 2.

  Il faut rÃ©pÃ©ter k fois la ligne 2, cependant on a pas besoin
  d'effectuer la derniÃ¨re permutation de la ligne 3. Rem: il est
  inutile de permuter les n cases de T, seule les n-1 premiÃ¨re cases
  suffisent car il s'agit de permutations circulaires.
*/
{
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_ADJ:
    return adjacency_rep(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    int k=Q->param[1]; if(k<1) RET_n(0); /* k>0 */
    Q->k=k; // pour adjacency_rep()
    ALLOC2(Q->rep,Q->n,k);
    NALLOCZ(int,T,Q->n,_i); /* tableau pour une permutation alÃ©atoire */

    int c,t;
    for(c=0;c<k;){ /* rÃ©pÃ¨te k fois, pour chaque cycle */
      for(t=0;t<Q->n;t++) Q->rep[T[t]][c]=T[(t+1)%Q->n]; /* suit et copie le cycle c */
      if(++c<k) Permute(T,Q->n-1); /* permutation alÃ©atoire (inutile le dernier coup) */
      /* seul Q->n-1 Ã©lÃ©ments ont besoin d'Ãªtre permutÃ©s (permutation circulaire) */
    }
    
    free(T);
    return 0;
  }
  
  return 1;
}


int margulis(query* const Q)
{
  int const n=Q->param[0]; 
  switch(Q->code){
    
  case QUERY_NAME:
    name_base(Q->name,Q->i,Q->param[0],2,",","()",1);
    return 0;

  case QUERY_INIT:;
    if(n<1) RET_n(0); /* n>0 */
    SET_n(n*n);
    return 0;
  
  case QUERY_ADJ:;
    int const x=Q->i%n;
    int const y=Q->i/n;
    int const u=Q->j%n;
    int const v=Q->j/n;
    if((u==(x+y)%n)  &&(v==y)) RET_a(1);
    if((u==(x-y+n)%n)&&(v==y)) RET_a(1);
    if((u==x)&&(v==(y+x)%n))   RET_a(1);
    if((u==x)&&(v==(y-x+n)%n)) RET_a(1);
    if((u==(x+y+1)%n)  &&(v==y)) RET_a(1);
    if((u==(x-y+1+n)%n)&&(v==y)) RET_a(1);
    if((u==x)&&(v==(y+x+1)%n))   RET_a(1);
    if((u==x)&&(v==(y-x+1+n)%n)) RET_a(1);
    RET_a(0);
  }

  return 1;
}


/***********************************

           GRAPHES DEFINIS
        A PARTIR D'UN TABLEAU
          (DE TAILLE FIXEE)

***********************************/


/* codes<0 pour GraphFromArray() */
enum{
  GFA_END   =-1, // fin du graphe
  GFA_CUT   =-2, // fin d'une sÃ©quence
  GFA_PATH  =-3, // chemin de sommets consÃ©cutifs
  GFA_HAM   =-4, // cycle Hamiltonien
  GFA_STAR  =-5, // Ã©toile
  GFA_WHEEL =-6, // camembert
};


int GraphFromArray(query* const Q,const int* const T)
/*
  Fonction d'adjacence gÃ©nÃ©rique permettant pour des graphes de petite
  taille et sans paramÃ¨tre. PlutÃ´t que d'utiliser une matrice ou une
  liste, c'est un tableau T qui dÃ©finit les adjacences du graphe.

  En gros, GraphFromArray(Q,T) fait Q->a=1 ssi Q->i et Q->j se suivent
  dans le tableau T, c'est-Ã -dire s'il existe un indice k tel que
  T[k]=Q->i et T[k+1]=Q->j ou le contraire.  Cette fonction suppose
  que Q->i < Q->j. Les arÃªtes du graphe sont ainsi couvertes par des
  chemins Ã©ventuellement non Ã©lÃ©mentaires. La fonction renvoie
  toujours 0 si tout c'est bien passÃ© et 1 sinon, comme une fonction
  d'ajacence standard.

  Les valeurs nÃ©gatives du tableau ont des significations
  particuliÃ¨res (voir ci-dessous la signification des codes GFA_xxx).
  Le principe de l'algorithme de dÃ©codage est le suivant: on lit le
  tableau T progressivement valeur aprÃ¨s valeur, et soit x le dernier
  sommet mÃ©morisÃ©, c'est-Ã -dire une valeur entiÃ¨re >=0. Si on lit une
  valeur y>=0 alors y est un sommet et le graphe possÃ¨de l'arÃªte
  x-y. On remplace alors x par y et on continue la lecture de T. Si
  y<0, c'est un code GFA_xxx et le comportement est alors spÃ©cifique
  comme dÃ©crit ci-dessous:

  GFA_END: fin du tableau et donc du graphe.

  GFA_CUT: fin d'une sÃ©quence. Le prochain sommet lu ne sera pas
     connectÃ© au dernier sommet mÃ©morisÃ© du tableau. Sans ce code, la
     prochaine valeur >=0 lue est toujours connectÃ©e au dernier sommet
     mÃ©morisÃ©. En fait, toute valeur <0 qui n'est pas un des codes
     reconnus Ã  le mÃªme effet que GFA_CUT.

  GFA_PATH: chemin de sommets consÃ©cutifs, "a,GFA_PATH,n" reprÃ©sente
      un chemin de n arÃªtes suivant le sommet a. C'est Ã©quivalent Ã 
      "a,a+1,...,a+n,GFA_CUT".

  GFA_HAM: cycle Hamiltonien, Ã©quivalent Ã 
  "0,GFA_PATH,Q->n-1,Q->n-1,0,GFA_CUT".

  GFA_STAR: Ã©toile. La sÃ©quence "a,GFA_STAR,a_1,...,a_n,GFA_CUT"
      reprÃ©sente une Ã©toile de centre a et de feuilles a_1,...,a_n.
      C'est Ã©quivalent Ã  "a,a_1,a,a_2,a,...,a,a_n,GFA_CUT". On peut
      remplacer le GFA_CUT par n'importe qu'elle valeur nÃ©gative qui
      peut ainsi Ãªtre enchaÃ®nÃ©e.

  GFA_WHEEL: comme GFA_STAR, sauf qu'en plus les sommets a_i,a_{i+1}
      sont adjacents. Attention ! a_n et a_1 ne sont pas adjacents.

  Une bonne stratÃ©gie pour trouver un code assez court pour un graphe
  (c'est-Ã -dire un tableau T de petite taille) est de l'Ã©plucher par
  degrÃ© maximum jusqu'Ã  obtenir des chemins ou cycles (non forcÃ©ment
  simples). Les sommets ainsi Ã©pluchÃ©s peuvent Ãªtre codÃ©s par GFA_STAR
  ou GFA_WHEEL, les chemins et de cycles par GFA_PATH.

    Ex:        5-6-7
               |  /      -> code = 0,GFA_PATH,10,GFA_CUT,4,8,GFA_END
       0-1-2-3-4-8-9-10

  Les sommets isolÃ©s n'ont pas Ã  Ãªtre codÃ©s. Cependant, pour
  dÃ©terminer Q->n il est nÃ©cessaire que le plus grand sommet du graphe
  figure dans T. Si cela n'est pas le cas (c'est-Ã -dire que le sommet
  n-1 se trouve Ãªtre un sommet isolÃ©), on peut ajouter la sÃ©quence
  "n-1,GFA_CUT," par exemple en dÃ©but du tableau T.
*/
{
  int const i=Q->i;
  int const j=Q->j;
  int k,a,n;
  switch(Q->code){

  case QUERY_INIT:
    for(k=n=0;T[k]!=GFA_END;k++) n=max(n,T[k]);
    RET_n(n+1);

  case QUERY_ADJ:
    for(k=0;;){
      switch(T[k]){
	  
      case GFA_END:
	RET_a(0); /* adjacence non trouvÃ©e */
	
      case GFA_PATH:
	a=T[k-1];
	n=T[++k];
	if((a<=i)&&(i<a+n)&&(j==i+1)) RET_a(1); /* suppose que i<j */
	k++;
	continue;
	
      case GFA_HAM:
	if(j==i+1) RET_a(1); /* suppose que i<j */
	if((i==0)&&(j==(Q->n-1))) RET_a(1);
	k++;
	continue;
	  
      case GFA_STAR:
	a=T[k-1];
	while(T[++k]>=0){
	  if((i==a)&&(j==T[k])) RET_a(1);
	  if((j==a)&&(i==T[k])) RET_a(1);
	}
	continue;
	  
      case GFA_WHEEL:
	a=n=T[k-1];
	while(T[++k]>=0){
	  if(((i==a)||(i==n))&&(j==T[k])) RET_a(1);
	  if(((j==a)||(j==n))&&(i==T[k])) RET_a(1);
	  n=T[k];
	}
	continue;
	
      default:
	/* ici T[k+1] existe car T[k]<>GFA_END. Si T[k]<0 alors cela
	   aura le mÃªme effet que de lire GFA_CUT, Ã  savoir le lien
	   T[k]-T[k+1] est coupÃ© puisque i et j sont >0. */
	if((i==T[k])&&(j==T[k+1])) RET_a(1);
	if((j==T[k])&&(i==T[k+1])) RET_a(1);
	k++;
	continue;
      }
    }
  }

  return 1;
}


int tutte(query* const Q){
  int const T[]={
    0,GFA_PATH,8,
    4,8,9,3,9,10,11,2,
    11,12,1,12,13,14,10,
    14,7,6,15,13,GFA_CUT,
    0,16,GFA_PATH,7,
    23,19,23,24,18,24,25,26,17,
    26,27,16,27,28,29,25,
    29,22,21,30,28,GFA_CUT,
    0,31,GFA_PATH,7,
    38,34,38,39,33,39,40,41,32,
    41,42,31,42,43,44,40,
    44,37,36,45,43,GFA_CUT,
    15,20,GFA_CUT,
    30,35,GFA_CUT,45,5,GFA_END};
  return GraphFromArray(Q,T);
}


int icosahedron(query* const Q){
  int const T[]={
    5,0,1,2,3,4,5,6,7,8,
    9,10,11,6,4,7,3,8,2,
    9,1,10,5,1,2,0,3,4,0,
    5,6,10,9,11,9,11,8,11,7,GFA_END};
  return GraphFromArray(Q,T);
}


int rdodecahedron(query* const Q){
  int const T[]={
    0,1,2,3,4,5,0,6,
    7,1,7,8,2,8,9,3,9,
    10,4,10,11,5,11,6,
    7,12,9,12,11,GFA_CUT,
    0,13,2,13,4,GFA_END};
  return GraphFromArray(Q,T);
}


int herschel(query* const Q){
  int const T[]={
    0,GFA_PATH,10,
    10,3,2,9,0,7,6,1,4,GFA_CUT,
    8,5,10,GFA_END};
  return GraphFromArray(Q,T);
}


int goldner_harary(query* const Q){
  int const T[]={
    0,1,2,0,3,4,1,5,6,7,2,8,3,
    1,6,2,3,9,6,10,1,9,2,10,GFA_CUT,
    4,9,5,GFA_CUT,
    8,9,7,GFA_END};
  return GraphFromArray(Q,T);
}


int triplex(query* const Q){
/* proche d'un flower_snark 3 */
  int const T[]={
    0,GFA_PATH,8,8,0,
    9,GFA_STAR,0,3,6,GFA_CUT,
    10,GFA_STAR,1,4,7,GFA_CUT,
    11,GFA_STAR,2,5,8,GFA_END};
  return GraphFromArray(Q,T);
}


int jaws(query* const Q){
  int const T[]={
    0,GFA_PATH,15,15,0,5,6,7,2,1,16,14,15,10,
    9,8,13,12,17,3,4,18,11,10,9,19,6,
    19,18,GFA_CUT,16,17,GFA_END};
  return GraphFromArray(Q,T);
}


int starfish(query* const Q){
/* proche d'un flower_snark 5 */
  int const T[]={
    0,1,2,3,4,0,GFA_CUT,
    5,6,7,8,9,5,GFA_CUT,
    10,11,12,13,14,10,GFA_CUT,
    15,GFA_STAR,0,5,10,GFA_CUT, 
    16,GFA_STAR,1,6,11,GFA_CUT, 
    17,GFA_STAR,2,7,12,GFA_CUT,
    18,GFA_STAR,3,8,13,GFA_CUT,
    19,GFA_STAR,4,9,14,GFA_END};
  return GraphFromArray(Q,T);
}


int fritsch(query* const Q){
  int const T[]={GFA_HAM,0,4,6,3,7,2,4,GFA_CUT,1,8,1,7,6,8,5,0,2,GFA_END};
  return GraphFromArray(Q,T);
}


int soifer(query* const Q){
  int const T[]={GFA_HAM,0,6,4,2,6,0,7,5,8,1,3,8,5,3,GFA_END};
  return GraphFromArray(Q,T);
}


int poussin(query* const Q){
  int const T[]={
    GFA_HAM,0,2,4,6,8,10,12,14,2,
    0,13,11,9,5,8,12,7,3,14,7,GFA_CUT,
    4,1,9,0,11,GFA_CUT,1,5,GFA_CUT,3,6,GFA_END};
  return GraphFromArray(Q,T);
}


int errera(query* const Q){
  int const T[]={
    0,GFA_STAR,10,9,5,15,11,16,GFA_CUT,
    1,GFA_STAR,6,16,12,13,4,7,GFA_CUT,
    2,GFA_WHEEL,6,7,8,9,10,GFA_CUT,
    3,GFA_WHEEL,11,12,13,14,15,11,GFA_CUT,
    4,GFA_STAR,7,8,5,14,13,GFA_CUT,
    5,GFA_STAR,8,9,14,15,GFA_CUT,
    16,GFA_STAR,6,11,12,10,GFA_END};
  return GraphFromArray(Q,T);
}


int kittell(query* const Q){
  int const T[]={
    0,GFA_PATH,22,
    10,12,14,5,13,11,3,10,2,9,1,6,0,
    5,16,6,17,7,1,8,22,20,18,7,GFA_CUT,
    5,15,21,14,22,9,12,GFA_CUT,
    9,14,GFA_CUT,18,8,20,GFA_CUT,
    19,GFA_STAR,15,16,17,21,GFA_CUT,
    2,0,3,0,4,13,4,11,GFA_END};
  return GraphFromArray(Q,T);
}


int frucht(query* const Q){
  int const T[]={GFA_HAM,0,4,5,3,2,7,6,8,9,11,10,1,GFA_END};
  return GraphFromArray(Q,T);
}


int moser(query* const Q){
  int const T[]={0,1,2,3,0,4,5,6,0,GFA_CUT,1,3,2,5,4,6,GFA_END};
  return GraphFromArray(Q,T);
}


int markstrom(query* const Q){
  int const T[]={
    8,0,GFA_PATH,12,5,13,GFA_PATH,7,
    3,21,22,23,0,23,1,2,21,22,18,GFA_CUT,
    15,20,14,13,4,GFA_CUT,
    6,12,7,GFA_CUT,17,19,GFA_CUT,
    9,11,10,16,GFA_END};
  return GraphFromArray(Q,T);
}


int robertson(query* const Q){
  int const T[]={GFA_HAM,0,4,8,13,17,2,6,10,14,0,
    1,9,16,5,12,1,GFA_CUT,3,11,18,7,15,3,GFA_END};
  return GraphFromArray(Q,T);
}


int heawood4(query* const Q){
  int const T[]={
    0,GFA_WHEEL,1,5,4,3,13,12,2,GFA_CUT,
    1,GFA_WHEEL,5,6,7,8,2,GFA_CUT,
    2,GFA_WHEEL,8,9,10,11,GFA_CUT,
    16,GFA_WHEEL,4,15,17,19,5,GFA_CUT,
    12,GFA_WHEEL,11,20,19,18,GFA_CUT,
    14,15,3,14,13,18,17,14,18,GFA_CUT,
    22,GFA_WHEEL,6,21,24,23,7,GFA_CUT,
    5,20,6,11,21,10,24,9,23,8,GFA_END};
  return GraphFromArray(Q,T);
}


int wiener_araya(query* const Q){
  int const T[]={
    0,GFA_STAR,1,4,12,15,GFA_CUT,
    1,GFA_PATH,39,
    41,GFA_STAR,20,23,36,GFA_CUT,
    18,36,35,16,17,1,2,19,GFA_CUT,
    3,21,22,5,4,8,7,26,27,9,10,29,28,39,25,24,6,GFA_CUT,
    8,12,11,31,32,13,14,34,33,37,38,23,GFA_CUT,
    30,40,33,GFA_END};
  return GraphFromArray(Q,T);
}


int zamfirescu(query* const Q){
  int const T[]={
    0,GFA_PATH,47,
    0,9,5,1,17,1,2,19,18,22,21,20,39,3,4,41,
    40,47,43,42,6,7,44,31,32,28,27,34,33,
    45,46,38,37,21,GFA_CUT,
    8,30,29,10,11,27,26,13,14,24,25,35,
    36,23,16,15,0,12,GFA_END};
  return GraphFromArray(Q,T);
}


int hatzel(query* const Q){
  int const T[]={
    0,GFA_PATH,55,
    46,30,31,11,12,8,7,26,27,9,10,29,28,
    47,48,44,45,38,39,43,56,52,51,23,22,
    5,6,24,25,50,49,56,GFA_CUT,
    8,4,0,12,13,32,33,37,36,40,41,42,53,
    54,20,21,3,2,19,18,55,41,GFA_CUT,
    14,34,35,16,17,1,0,15,GFA_END};
  return GraphFromArray(Q,T);
}


int harborth(query* const Q){
  int const T[]={
    0,GFA_PATH,19,GFA_CUT,19,
    0,20,1,21,2,22,3,23,4,24,5,
    25,6,26,7,27,8,28,9,29,10,
    30,11,31,12,32,13,33,14,34,15,
    35,16,36,17,37,18,38,19,39,0,
    39,38,40,41,37,36,35,42,41,43,
    44,45,21,20,45,43,40,39,GFA_CUT,
    46,44,22,23,24,46,25,26,27,
    47,48,28,29,48,49,47,46,GFA_CUT,
    50,49,51,30,31,51,50,32,33,34,
    42,50,GFA_END};
  return GraphFromArray(Q,T);
}


int doily(query* const Q){
  int const T[]={
    9,0,GFA_PATH,9,GFA_CUT,
    9,11,3,13,7,10,1,12,5,14,9,GFA_CUT,
    0,13,5,GFA_CUT,2,14,7,GFA_CUT,
    4,10,9,GFA_CUT,6,11,1,GFA_CUT,
    8,12,3,GFA_END};
  return GraphFromArray(Q,T);
}


int cricket(query* const Q){
  int const T[]={0,1,2,3,1,4,GFA_END};
  return GraphFromArray(Q,T);
}


int moth(query* const Q){
  int const T[]={0,1,2,3,4,1,5,1,3,GFA_END};
  return GraphFromArray(Q,T);
}


int dart(query* const Q){
  int const T[]={0,1,2,3,4,1,3,GFA_END};
  return GraphFromArray(Q,T);
}


int antenna(query* const Q){
  int const T[]={0,1,2,3,4,5,1,2,5,GFA_END};
  return GraphFromArray(Q,T);
}


int suzuki(query* const Q){
  int const T[]={
    0,GFA_WHEEL,1,2,3,4,5,6,7,8,GFA_CUT,
    8,2,4,6,8,
    1,3,5,7,1,
    5,6,10,1,9,
    4,3,9,10,7,
    8,10,5,9,2,
    GFA_END};
  return GraphFromArray(Q,T);
}


int bull(query* const Q){
  int const T[]={0,1,2,3,4,3,1,GFA_END};
  return GraphFromArray(Q,T);
}


int hgraph(query* const Q){
  int const T[]={0,1,2,3,GFA_CUT,4,1,2,5,GFA_END};
  return GraphFromArray(Q,T);
}


int rgraph(query* const Q){
  int const T[]={0,1,2,3,4,1,5,GFA_END};
  return GraphFromArray(Q,T);
}


/***********************************

        FONCTIONS DU PROGRAMME
              PRINCIPAL

***********************************/

int InitVertex(int const n,double const p)
/*
  Remplit le tableau V[i] donnant l'Ã©tiquette finale du sommet i et
  supprime certains des n sommets suivant la valeur p. Si p<0, alors
  on supprime exactement |p| sommets. Utilise aussi SHIFT.  Si PERMUTE
  est vrai V[] est remplit d'une permutation alÃ©atoire de
  SHIFT+[0,n[. Si V[i]=-1 cela signifie que i a Ã©tÃ© supprimÃ© (p).  La
  fonction retourne le nombre de sommets final du graphe, c'est-Ã -dire
  le nombre d'Ã©tiquettes >=0. Si k sommets ont Ã©tÃ© supprimÃ©s, alors
  les valeurs de V[] sont dans SHIFT+[0,n-k[.

  Initialise aussi le tableau VF[j], avec j=0...n-k, de sorte que
  VF[j]=i si VF[i] est le j-Ã¨me sommet non supprimÃ©. Il est important
  que VF[] ait une taille de Q->n au dÃ©part. Un realloc() le
  redimensionne plus tard dans la fonction.
*/
{
  int i,j,r; // r=n-(nb de sommets supprimÃ©s)

  /* supprime les sommets */

  if(p<0){ // ici on en supprime exactement |p| sommets
    for(i=0;i<n;i++) VF[i]=i; // on va se servir temporairement de VF
    r=-(int)p;  // les r premiÃ¨res valeurs de VF seront les sommets Ã  supprimer
    r=min(r,n); // r ne doit pas dÃ©passer n
    for(i=0;i<r;i++){
      j=i+randomu(n-i); // ressemble Ã  Permute() mais pas tout Ã  fait !
      SWAP(VF[i],VF[j]);
    }
    for(i=0;i<r;i++) V[VF[i]]=-1; // on supprime ces r sommets
    for(i=r=0;i<n;i++) // on remplit V et VF
      if(V[i]>=0) { VF[r]=i; V[i]=r++; }
  }
  else{ /* ici on supprime chacun des n sommets avec proba p */
    long const seuil=(double)p*(double)RAND_MAX;
    for(i=r=0;i<n;i++)
      if(random()<=seuil) V[i]=-1;
      else { VF[r]=i; V[i]=r++; }
  } /* dans les deux cas r = nombre de sommets restant */

  /* rÃ©ajuste le tableau VF Ã  la taille minimum */

  REALLOC(VF,r);
  if(PERMUTE) Permute(V,n);

  /* ne rien faire si SHIFT=0 */

  if(SHIFT>0)
    for(i=0;i<r;i++)
      V[VF[i]] += SHIFT;

  return r;
}


void ScanINC(int* const dmax,int* const dmin,int* const m)
/*
  Calcule, en fonction des tableaux INC[] et VF[], le degrÃ© max, le
  degrÃ© min et le nombre d'arÃªtes du graphe final.
*/
{
  int d,i;
  *m=*dmax=0;
  *dmin=INT_MAX;
  if(NF<1) *dmin=0;
  for(i=0;i<NF;i++){
    d=INC[VF[i]]; /* d=degrÃ© du i-Ã¨me sommet existant */
    *dmax=max(*dmax,d);
    *dmin=min(*dmin,d);
    *m += d;
  }
  *m >>= 1;
  return;
}


string MakeCMD(string s,int const deb,int const fin)
/*
  Routine permettant de recomposer la ligne de commande, en
  particulier pour obtenir le nom du graphe gÃ©nÃ©rÃ© (avec les
  options). On ajoute Ã  la fin de la chaÃ®ne s les arguments ARGV[i]
  pour iâˆˆ[deb,fin[. Si s=NULL alors un pointeur statique sur la chaÃ®ne
  est renvoyÃ©e, pointeur qui n'a donc pas besoin d'Ãªtre libÃ©rÃ©.

  Si un argument comprend un espace, il est alors parenthÃ©sÃ© par '...'
  de faÃ§on Ã  Ãªtre interprÃ©tÃ© comme un seul argument. Les arguments
  sont sÃ©parÃ©s par un espace. Le dernier argument est toujours suivi
  d'un espace, Ã©ventuellement inutile.
*/
{
  static char r[CMDMAX];
  if(s==NULL){ s=r; VIDE(s); }
  
  int i;
  for(i=deb;i<fin;i++)
    if(index(ARGV[i],' ')) /* si l'argument est en plusieurs mots */
      strcat(strcat(strcat(s,"'"),ARGV[i]),"' ");
    else strcat(strcat(s,ARGV[i])," ");

  return s;
}


string DateHeure(void)
/*
  Renvoie la date et l'heure courante sous forme de chaÃ®ne de
  caractÃ¨res. Il ne faut pas faire de free() sur le pointeur retournÃ©.
*/
#define DATE_FORMAT "%d/%m/%Y - %Hh%M'%S"
#define SIZE_DATE 22 // avec le 0 final
{
  static char date[SIZE_DATE];
  time_t t=time(NULL);
  struct tm *tm=localtime(&t);
  strftime(date,SIZE_DATE*sizeof(char),DATE_FORMAT,tm);
  date[SIZE_DATE-1]=0;
  return date;
}
#undef DATE_FORMAT
#undef SIZE_DATE


void Header(int const code)
/*
  Affiche et calcule un prÃ©ambule (date, nom du graphe, nombre de
  sommets). Suivant la valeur de code, le nombre d'arÃªtes est donnÃ©e.
  Si le bit-0 de code est Ã  1, alors l'en-tÃªte de base est affichÃ©e.
  Si le bit-1 de code est Ã  1, alors on affiche le nombre d'arÃªtes, le
  degrÃ© min et max.

*/
{
  /* affichage de la date, nom du graphe, ligne de commande et de n */
  if(code&1){
    printf("//\n");
    printf("// %s - seed=%u\n",DateHeure(),SEED);
    printf("// %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libÃ©rer ce pointeur
    printf("// n=%i",NF);
  }

  /* affichage du nombre d'arÃªtes, maxdeg et mindeg */
  if(code&2){
    int maxdeg,mindeg,nbedges;
    ScanINC(&maxdeg,&mindeg,&nbedges);
    if(!(code&1)) printf("//\n//");
    printf(" m=%i",nbedges);
    printf(" maxdeg=%i",maxdeg);
    printf(" mindeg=%i",mindeg);
  }

  if(code) printf("\n//\n");
  return;
}


static inline string ComputeName(query* const Q, int const i)
/*
  Routine servant plusieurs fois dans Out(Q) pour l'affichage Ã  la
  volÃ©e du sommet i. On renvoie le nom du sommet Ã  afficher dans
  Q->name. Par dÃ©faut c'est V[i]. Dans le cas du format standard on
  affiche le nom donnÃ© par Q->adj(Q) et QUERY_NAME, ou bien les
  coordonnÃ©es Q->xpos,Q->ypos suivant le cas (|LABEL|==3).

  NB: Dans le cas du format dot, il faut laisser les indices (donc
  V[i]) car les "label"s des sommets sont codÃ©s seulement Ã  la fin
  dans l'attribut "label".

  Modifie: Q->i, Q->code, Q->name
*/
{
  sprintf(Q->name,"%i",(V==NULL)?i:V[i]); /* par dÃ©faut name="V[i]" */
  if(FORMAT==F_standard){
    if(abs(LABEL)==1){ Q->code=QUERY_NAME; Q->i=i; Q->adj(Q); }
    if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",Q->xpos[i],Q->ypos[i]);
  }
  return Q->name;
}


void Out(query* const Q)
/*
  Affiche l'arÃªte i-j suivant le format FORMAT.

  Si HEADER=1, alors Out() doit se charger d'afficher un en-tÃªte.

  Si CHECK>0 alors Out() doit se charger de crÃ©er et de remplir la
  liste d'adjacence du graphe GF et de dÃ©terminer son nombre de
  sommets NF. Pour cela une liste (type "list") est crÃ©ee et
  progressivement rempli. Ã€ la fin, on calcule GF avec List2Graph().

  Variables globales Ã©ventuellement modifiÃ©es:
  - N, GF, NF, VF, NAME
  - CAPTION, NPAL, PALETTE
  - VERTEX0, VERTEXN
  - USERDOT, ROUND

  Autres variables globales utilisÃ©es:
  - CHECK, FORMAT, WIDTH
  - HEADER, VCOLOR
  - XMIN, XMAX, YMIN, YMAX
  - VSIZEK, VSIZESTD, VSIZEXY
  - POS, LABEL, LEN,
  - PARAM_PAL, CPARAM
  - COLOR_CHAR, COLOR_RGB
*/
{
  int x,y,z;
  double w;
  
  /* variables qui conservent leurs valeurs d'un appel Ã  l'autre */
  static list *L0,*L1,*L2; // pour -check: tÃªte, dernier, avant-dernier
  static int cpt;  // compteur d'arÃªtes ou de sommets isolÃ©s par ligne
  static int last; // extrÃ©mitÃ© de la derniÃ¨re arÃªte affichÃ©e
  static string sep1; // sÃ©parateur d'arÃªte: "" (pour standard) ou ";" (dot)
  static string sep2; // autre sÃ©parateur
  static string sepe; // caractÃ¨re pour une arÃªte: "-" (pour standard) ou "--" (dot)

  int const i=Q->i;
  int j=Q->j;

  switch(Q->code){


  case QUERY_INIT:
    //------------------------------------
    // initialise la fonction
    //------------------------------------

    cpt=0;
    last=-1;
    if(CHECK) L0=L1=L2=new_list(); /* initialise la liste */

    switch(FORMAT){

    case F_standard:
      if(HEADER) Header(1);
      sep1="";
      sep2=" ";
      sepe=(Q->directed)?"->":"-";
      return;

    case F_userdot:
      USERDOT.i=USERDOT.j=-1;
      USERDOT.adj=NULL;
      USERDOT.ptr=NULL;
      // puis pareil que F_dot
    case F_dot:
      if(HEADER) Header(1);
      printf("%sgraph {\n",(Q->directed)?"di":"");
      if(CAPTION){
	printf("\tgraph [label=\"%s\"];\n",CAPTION);
	free(CAPTION);
	CAPTION=NULL;
      }
      if(DOTSCALE){ // si "auto" on calule une Ã©chelle en fontion de âˆšn ou des BB
	int b=(strcmp(DOTSCALE,"auto")==0); // =1 ssi on fait un alloc()
	if(b){
	  ALLOC(DOTSCALE,20); // alloue de la place pour Ã©crire le facteur
	  if(POS) sprintf(DOTSCALE,"%g",1.0/(1+max(XMAX-XMIN,YMAX-YMIN)));
	  else sprintf(DOTSCALE,"%g",1.0/(1+(int)(sqrt(Q->n))));
	}
	printf("\tgraph [scale=\"%s\"];\n",DOTSCALE);
	if(b){ free(DOTSCALE); DOTSCALE=NULL; }
      }
      if(VCOLOR&0x10){ /* "list" */
	printf("\tgraph [layout=nop, splines=line];\n");
	printf("\tnode [height=1.0, width=0.4, style=filled, shape=rectangle];\n");
	return;
      }
      if(POS){
	w=PosAspect(Q);
	/* layout=nop: pour dire Ã  dot de ne pas re-calculer les positions */
	printf("\tgraph [layout=nop, splines=line, bb=\"%.2lf,%.2lf,%.2lf,%.2lf\"];\n",
	       w*XMIN-2*VSIZEK*VSIZEXY,w*YMIN-2*VSIZEK*VSIZEXY,
	       w*XMAX+2*VSIZEK*VSIZEXY,w*YMAX+2*VSIZEK*VSIZEXY);
	/* si on ne met pas le "2*" on a des sommets tronquÃ©s ... */

	/* affiche Ã©ventuellement une sous-grille carrÃ©e, avant le
	   graphe pour qu'elle apparaisse en arriÃ¨re-plan. Le nombre
	   de points vaut XYgrid, Q->n ou sqrt(Q->n) */
	z=Q->n; // 1er sommet disponible aprÃ¨s ceux de G
	if(XYgrid<0) XYgrid=(XYtype==XY_PERM)? Q->n : 1+(int)(sqrt(Q->n));
	if(XYgrid>1){ // il faut au moins deux lignes et colonnes
	  printf("\n\tsubgraph {\n");
	  printf("\t\tnode [label=\"\", height=0, width=0, shape=point, style=invis];\n");
	  printf("\t\tedge [color=gray, penwidth=0.6];");
	  double const rx=(double)(XMAX-XMIN)/(XYgrid-1); // pas de la grille en X
	  double const ry=(double)(YMAX-YMIN)/(XYgrid-1); // pas de la grille en Y
	  for(y=0;y<XYgrid;y++){
	      for(x=0;x<XYgrid;x++){
		/* affiche le sommet courant (x,y) ainsi que deux
		   arÃªtes incidentes vers les voisins (x+1,y) et
		   (x,y+1), s'ils existent */
		printf("\n\t\t%i [pos=\"%lf,%lf\"];",z,
		       w*(XMIN+(double)x*rx),
		       w*(YMIN+(double)y*ry));
		if(x+1<XYgrid) printf("\t%i--%i;",z,z+1); // arÃªte vers (x+1,y)
		if(y+1<XYgrid) printf("\t%i--%i;",z,z+XYgrid); // arÃªte vers (x,y+1)
		z++; // prochain sommet
	      }
	    }
	  printf("\n\t}\n\n");
	}
	if(XYzero){ /* ajoute l'origine (0,0) */
	  printf("%s\tsubgraph {\n",(XYgrid>1)?"":"\n");
	  printf("\t\tzero [label=\"\", pos=\"0,0\", shape=circle, ");
	  printf("height=0.1, width=0.1, color=red, fixedsize=true];\n");
	  printf("\n\t}\n\n");
	}
	if(XYborder){ /* ajoute un bord */
	  string s="rectangle";
	  if((XYtype=XY_DISK)||(XYtype==XY_CIRCLE)||(XYtype==XY_HYPER)) s="circle";
	  printf("%s\tsubgraph {\n",(XYgrid>1)||(XYzero)?"":"\n");
	  printf("\t\tbord [label=\"\", pos=\"0,0\", shape=%s,",s);
	  printf(" height=1.0, width=1.0, color=blue, fixedsize=true];\n");
	  printf("\n\t}\n\n");
	}
      }
      printf("\tnode [");
      if(LABEL<=0) printf("label=\"\", shape=point, "); /* sans label centrÃ© */
      w=POS?VSIZESTD:VSIZEXY; /* taille des sommets */
      printf("height=%.2lf, width=%.2lf];\n",w,w);
      if(strcmp(DOTFILTER,"neato")==0) printf("\tedge [len=%.2lf];\n",LEN);
      sep1=";";
      sep2="; ";
      sepe=(Q->directed)?"->":"--";
      return;

    case F_html:
      printf("<!DOCTYPE html>\n");
      printf("<html lang=\"fr\">\n");
      printf("<head>\n");
      printf("<meta charset=\"utf-8\" />\n");
      printf("<title>%s</title>\n\n",MakeCMD(NULL,0,ARGC));
      printf("<script>\n");
      printf("  function load_script(url){ document.write('<script src=\"'+url+'\"><\\/script>'); }\n");
      printf("  function load_style(url){\n");
      printf("    var d=document;\n");
      printf("    var c=d.createElement(\"link\");\n");
      printf("    c.rel=\"stylesheet\";\n");
      printf("    c.href=url;\n");
      printf("    d=d.head;\n");
      printf("    d.insertBefore(c,d.childNodes[d.childNodes.length-1].nextSibling);\n");
      printf("    }\n");
      printf("</script>\n\n");
      printf("<script src=\"%s\"",URL_vis_js1);
      printf(" onerror=\"load_script('%s');\"></script>\n\n",URL_vis_js2);
      printf("<link href=\"%s\" rel=\"stylesheet\" type=\"text/css\"",URL_vis_css1);
      printf(" onerror=\"load_style('%s');\" />\n\n",URL_vis_css2);
      printf("<style type=\"text/css\">\n");
      printf("  html, body { padding: 0; margin: 0; width: 100%%; height: 100%% }\n");
      printf("  #G { width: 100%%; height: 100%% }\n");
      printf("</style>\n\n");
      printf("</head>\n\n");
      printf("<body>\n");
      printf("<div id=G></div>\n");
      printf("<script>\n\n");
      printf("var E = new vis.DataSet([\n"); // pour les arÃªtes
      return;

    }
    return; // termine QUERY_INIT


  case QUERY_END:
    //----------------------------------
    // termine la fonction
    //----------------------------------

    if(CHECK){ /* on crÃ©e GF en fonction de la liste L0 */
      free(L1); /* supprime la sentienelle (dernier Ã©lÃ©ment) de L0 */
      if(L0==L1) GF=new_graph(0); /* si premier = dernier alors graphe vide */
      else{
	L2->next=NULL; /* coupe la liste Ã  l'avant dernier Ã©lÃ©ment qui a Ã©tÃ© supprimer */
	GF=List2Graph(L0,4); /* initialise GF et NF */
      }
      NF=GF->n;
    }
    
    switch(FORMAT){

    case F_standard:
    case F_userdot:
    case F_dot:
      if((VCOLOR&0x10)==0){ /* court-circuite l'affichage du graphe si "-vcolor list" */
	if(cpt>0) printf("%s\n",sep1); /* newline si fini avant la fin de ligne */
	if(FORMAT==F_standard){ /* fin si format standard */
	  if(HEADER) Header(2);
	  return;
	}

	if(POS||LABEL){ // affiche les labels des sommets (voir aussi ComputeName)
	  // NB: si POS=1, alors Q->xpos/Q->ypos existent forcÃ©ment
	  w=PosAspect(Q);
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    printf("%i [",V[x]);
	    if(POS) printf("pos=\"%lf,%lf\"",w*Q->xpos[x],w*Q->ypos[x]);
	    if(LABEL){
	      VIDE(Q->name);
	      if(abs(LABEL)==1){
		Q->code=QUERY_NAME;
		Q->i=x;
		Q->adj(Q); /* calcule Q->name */
	      }
	      if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",Q->xpos[x],Q->ypos[x]);
	      printf("%s%slabel=\"%s\"",
		     POS?", ":"",
		     (LABEL<0)? "x": "",
		     NONVIDE(Q->name)? Q->name : "\\N");
	    }
	    printf("];\n");
	  }
	}
	
	if(VSIZE&&(NF>0)){ /* taille en fonction du degrÃ© */
	  double alpha,smin;
	  ScanINC(&x,&z,&y); /* x=degmax, z=degmin */
	  smin=POS?VSIZESTD:VSIZEXY;
	  alpha=(x==z)? 0 : smin*(VSIZEK-1)/((double)(x-z));
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    w=smin + alpha*(INC[x]-z);
	    printf("%i [height=%.2lf, width=%.2lf];\n",V[x],w,w);
	  }
	}
      } /* fin du if((VCOLOR&0x10)==0) ... */

      if(VCOLOR&&(NF>0)){ /* couleur des sommets */
	color c={0,0,0},*P; /* couleur noire par dÃ©faut */
	int *D;

	if(VCOLOR&0x8){ /* option "pal" on initialise la PALETTE */
	  NPAL=(int)strlen(PARAM_PAL);
	  if(NPAL==1) { /* PARAM_PAL="x" alors PARAM_PAL="xx" */
	    PARAM_PAL[1]=PARAM_PAL[0];
	    PARAM_PAL[NPAL=2]='\0';
	  }
	  ALLOC(PALETTE,NPAL); /* PALETTE = tableau de NPAL "color" */
	  for(y=z=0;y<NPAL;y++){ /* z=prochaine entrÃ©e libre dans PALETTE */
	    x=(int)(index(COLOR_CHAR,PARAM_PAL[y])-COLOR_CHAR); /* x=index de la couleur */
	    x/=sizeof(*COLOR_CHAR); /* normalement inutile */
	    if((0<=x)&&(x<COLOR_NB)) PALETTE[z++]=COLOR_RGB[x]; /* on a trouvÃ© la couleur */
	  }
	  if(z<2){ PALETTE[0]=PALETTE[1]=c; z=2; } /* si pas trouvÃ© au moins 2 couleurs */
	  NPAL=z; /* NPAL=nb de couleurs trouvÃ©es */
	}
	
	if(VCOLOR&0x17){ /* fonction de couleur: 1,2,3,4,5 ou "-vcolor list" */
	  if((VCOLOR&0x7)>2){ /* si 3="degm", 4="randg"  ou 5="kcolor" */
	    if((VCOLOR&0x7)==3){ /* si "degm" */
	      int *T=Prune(GF,NULL);
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==4){ /* si "randg" */
	      NALLOCZ(int,T,NF,_i);
	      Permute(T,NF); /* T=permutation alÃ©atoire */
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==5){ /* si "kcolor" */
              y=MEM(CPARAM,0,int); /* y=nb de couleur */
	      D=kColor(GF,y);
              if(D==NULL){ ALLOCZ(D,NF,0); y=1; } /* pas trouvÃ© -> une seule couleur */
	    }
	  }
	  else{ /* si 1="deg" ou 2="degr" */
	    ScanINC(&x,&z,&y); /* calcule x=degmax, z=degmin */
	    y=x-z+1; /* y=nb a priori de couleurs nÃ©cessaires */
	    ALLOCZ(D,NF,INC[VF[_i]]-z);
	    if((VCOLOR&0x7)==2){ /* si "degr" */
	      int *R=SortInt(D,NULL,Q->n,0,&y,SORT_INC_RANK);
	      /* aprÃ¨s SortInt(), y=nb de degrÃ© diffÃ©rents */
	      free(D);
	      D=R; /* on substitue R Ã  D */
	    }
	  }
	  /* ici D[x]=indice de la couleur du sommet x, et y=nb de couleurs */
	  P=GradColor(PALETTE,NPAL,y); /* calcule une palette P de y couleurs. NB: ici NPAL>1 */
	  printf("\n");
	  if(VCOLOR&0x10){ /* si "-vcolor list" */
	    for(x=0;x<y;x++){
	      c=P[x];
	      for(z=0;z<COLOR_NB;z++) /* on cherche c dans COLOR_RGB */
		if((COLOR_RGB[z].r==c.r)&&(COLOR_RGB[z].g==c.g)&&(COLOR_RGB[z].b==c.b)) break;
	      printf("\t%i [pos=\"%i,0\", color=\"#%.02x%.02x%.02x\", label=\"%c\", fontcolor=%s];\n",
		     x,10+28*x,c.r,c.g,c.b,(z<COLOR_NB)?COLOR_CHAR[z]:' ',
		     (c.r+c.g+c.b<150)?"white":"black");
	    }
	  }else{ /* si pas "-vcolor list" */
	    for(x=0;x<NF;x++){
	      c=P[D[x]]; /* c=couleur du degrÃ© (ou du rang) de x */
	      printf("%i [style=filled, fillcolor=\"#%02x%02x%02x\"];\n",V[VF[x]],c.r,c.g,c.b);
	    }
	  }

	  free(D);
	  free(P);
	}
	
	if(VCOLOR&0x8){
	  free(PALETTE); /* la PALETTE ne sert plus Ã  rien */
	  PALETTE=COLOR_RGB; /* remet Ã  la valeur par dÃ©faut, au cas oÃ¹ */
	  NPAL=COLOR_NB; /* taille de la palette par dÃ©faut */
	}
      }

      printf("}\n"); /* si F_dot on ferme le "}" */
      if(HEADER) Header(2); /* affiche les arÃªtes */
      return;
      
    case F_no:
      if(HEADER) Header(3);
      return;

    case F_list:
      if(HEADER) Header(3);
      PrintGraphList(GF);
      return;
      
    case F_matrix:
    case F_smatrix:
      if(HEADER) Header(3);
      PrintGraphMatrix(GF);
      return;

    case F_xy:;
      char fmt[17]="%lf %lf\n"; // format d'affichage Q->xpos/Q->ypos par dÃ©faut , soit 6 digits
      if(HEADER) Header(3);
      if((Q->xpos==NULL)||(Q->ypos==NULL)) Erreur(8);
      printf("%i\n",NF); // nombre de sommets final
      if(ROUND>=DBL_DIG){ // par dÃ©faut vÃ©rifie si que coordonnÃ©es entiÃ¨res
	x=0; while(x<NF && Q->xpos[x]==floor(Q->xpos[x]) && Q->ypos[x]==floor(Q->ypos[x])) x++;
	if(x==NF) ROUND=0; // si que des coordonnÃ©es entiÃ¨res, on fait comme si ROUND=0
      }
      if(ROUND<DBL_DIG){ // construit le nouveau format
	int r=max(ROUND,0); // si ROUND<0 -> 0
	if(r<10){
	  strcpy(fmt,"%.0*lf %.0*lf\n"); // remplace '*' par r
	  fmt[3]=fmt[10]=(char)('0'+r);
	}else{
	  strcpy(fmt,"%.0**lf %.0**lf\n"); // remplace '**' par r
	  fmt[3]=fmt[11]=(char)('0'+(r/10));
	  fmt[4]=fmt[12]=(char)('0'+(r%10));
	}
      }
      for(y=0;y<NF;y++){ // affiche les coordonnÃ©es au format fmt
	x=VF[y]; // le sommet x existe
	printf(fmt,Q->xpos[x],Q->ypos[x]);
      }
      return;
      
    case F_html:
      printf("]);\n\n"); // fin de E = {...}
      printf("var V = new vis.DataSet([\n");
      for(y=0;y<NF;y++){
	x=VF[y]; /* le sommet x existe */
	VIDE(Q->name);
	if(LABEL){
	  if(abs(LABEL)==1){
	    Q->code=QUERY_NAME;
	    Q->i=x;
	    Q->adj(Q); /* calcule Q->name */
	  }
	  if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",Q->xpos[x],Q->ypos[x]);
	}
	if(!NONVIDE(Q->name)) sprintf(Q->name,"%i",x);
	printf("{id: %i, %s: \"%s\"},\n",
	       x,
	       (LABEL)?"label":"title",
	       Q->name);
      }
      printf("]);\n\n");
      printf("var container = document.getElementById('G');\n");
      printf("var data = { nodes: V, edges: E };\n");
      printf("var options = {\n");
      printf("  interaction: {\n");
      printf("    hover: true,\n");
      printf("    multiselect: true,\n");
      printf("    selectable: true\n");
      printf("    },\n");
      printf("  manipulation: { enabled: true },\n"); // pour l'edition
      printf("  layout: { randomSeed: %u },\n",SEED); // pour un dessin dÃ©pendant de la seed de gengraph
      printf("  nodes: {\n");
      printf("    shape: '%s',\n",(LABEL>0)?"ellipse":"dot");
      printf("    size: 11,\n"); // taille si pas de label interne
      printf("    color: {\n");
      printf("      border: 'black', background: '%s',\n",(LABEL>0)?"#dddddd":"black");
      printf("      highlight: { border: 'red', background: '#ff7f7f'},\n");
      printf("      hover: { border: 'blue', background: 'lightblue' }\n");
      printf("      }\n");
      printf("    },\n");
      printf("  edges: {\n");
      printf("    width: 2,\n");
      printf("    color: { color: 'black', highlight: 'red', hover: 'blue' }\n");
      printf("    }\n");
      printf("  };\n");
      printf("new vis.Network(container, data, options);\n");
      printf("\n</script>\n");
      printf("</body>\n");
      printf("</html>\n");
      return;

    default: Erreur(5); /* normalement sert Ã  rien */

    }
    return; // termine QUERY_END 

    
  case QUERY_ISOL: j=-1;
  case QUERY_ADJ:
    //-----------------------------------------
    // affichage Ã  la volÃ©e: "i-j", "i" ou "-j"
    //-----------------------------------------

    if(CHECK){
      L1=Insert(L2=L1,i,T_NODE); // on ajoute i
      if(j>=0) L1=Insert(L2=L1,j,(Q->directed)?T_ARC:T_EDGE); // on ajoute ->j ou -j
      if(VCOLOR&0x10) return; /* ne rien faire d'autre si "-vcolor list".
				 NB: CHECK>0 dans ce cas */
    }

    switch(FORMAT){

    case F_standard:
    case F_dot:
      if(j<0) last=-1; /* sommet isolÃ©, donc last sera diffÃ©rent de i */
      if(i!=last) printf("%s%s",(cpt==0)?"":sep2,ComputeName(Q,i));
      last=j; /* si i=last, alors affiche -j ou ->j. Si j<0 alors last<0 */
      if(j>=0) printf("%s%s",sepe,ComputeName(Q,j));
      if(++cpt==WIDTH){
	cpt=0; last=-1;
	printf("%s\n",sep1);
      }
      return;
      /* si format matrix, smatrix etc., ne rien faire, c'est fait par QUERY_END */
    
    case F_userdot: // NB: ici on suppose que WIDTH=1
      if(j<0) printf("%s%s\n",ComputeName(Q,i),sep1); // sommet isolÃ©
      else{
	Q->code=QUERY_DOT;
	if((USERDOT.i!=i)||(USERDOT.j!=j)) USERDOT.adj(Q); // l'arÃªte n'a pas Ã©tÃ© calculÃ©e
	// ici l'arÃªte i-j vient d'Ãªtre calculÃ©e, il reste Ã  la dessiner
	USERDOT.adj(Q);
      }
      return;

    case F_html: // NB: ici on suppose que WIDTH=1
      printf("{from: %i, to: %i},\n",i,j);
      return;

    }
    return; // termine QUERY_ADJ

  }
  return; // termine Out()
}


static inline double CheckProba(double const p)
/* VÃ©rifie que p est bien une probabilitÃ© */
{
  if((p<0)||(p>1)) Erreur(38);
  return p;
}


/***************************************

           FONCTIONS LIÃ‰ES
  A L'ANALYSE DE LA LIGNE DE COMMANDE

***************************************/


void Grep(int i)
/*
  Cherche le mot ARGV[i] dans l'aide contenu dans le source du
  programme, puis termine par exit(). Utilise les commandes: sed,
  tail, awk.

  Plusieurs cas peuvent se prÃ©senter:

  Cas 0: si un ARGV[j]="-", pour j=0..i, on le saute car cela ne peut
   Ãªtre ni une option ni un graphe dont on peut trouver l'aide.
  
  Cas 1: ARGV[i] est une option, ou aucun ARGV[j] avec j<i n'est pas
   une option.  Alors on affiche l'aide allant de "^....ARGV[i]($| )"
   Ã  "^$".

  Cas 2: ARGV[j] est une option mais pas ARGV[i] avec j<i. Dans ce
   cas, on pose mot=ARGV[j]" "ARGV[j+1]" "..." "ARGV[i]. Alors on
   affiche l'aide allant de "[ ]{7}mot($| )" Ã  "^$" ou "^[ ]{7}-fin"
   avec fin dÃ©fini de sorte que mot ne soit pas un prÃ©fixe.

  En fait on bride la recherche de l'option Ã  ARGV[i-1] ou ARGV[i-2]
  seulement, si bien que j=i, i-1 ou i-2.
*/
{
  int j,k,t;

  // enlÃ¨ve les arguments "-" de la liste des arguments

  for(t=j=0;t<=i;t++)
    if(strcmp(ARGV[j],"-")) j++;
    else ARGV[j]=ARGV[t],i--;
  j=i;

  DEBUG(
	for(t=0;t<=i;t++) printf("%s ",ARGV[t]);
	printf("\n");
	fflush(stdout);
	);
  
  do{ // calcule j=i, i-1 ou i-2
    if((i>0)&&(*ARGV[i-0]=='-')){ j=i-0; break; }
    if((i>1)&&(*ARGV[i-1]=='-')){ j=i-1; break; }
    if((i>2)&&(*ARGV[i-2]=='-')){ j=i-2; break; }
  }while(0);

  // construit mot

  NALLOC(char,mot,CMDMAX); VIDE(mot);
  for(t=j;t<=i;t++){
    strcat(mot,ARGV[t]);
    if(t<i) strcat(mot," ");
  }

  // construit fin (Ã  faire seulement si j<i)
  // si ARGV[j]="-option abc xy"
  // alors fin="-option ([^a]|a[^b]|ab[^c]abc) ([^x]|x[^y])"

  NALLOC(char,fin,CMDMAX); VIDE(fin);
  if(j<i){
    strcat(fin,ARGV[j]);
    strcat(fin," (");
    t=j+1; k=0;
    while(ARGV[t][k]){
      strncat(fin,ARGV[t],k);
      strcat(fin,"[^");
      strncat(fin,ARGV[t]+k,1);
      strcat(fin,"]|");
      k++;
      if(ARGV[t][k]==0){
	if(t==i) break;
	strcat(fin,ARGV[t]);
	strcat(fin,") (");
	t++; k=0;
      }
    }
    fin[strlen(fin)-1]=')';
  }

  // construit la commande s

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c|");
  /* rem: sed -E n'est pas standard */

  if(j==i)
    strcat(strcat(strcat(s,"sed -nE '/^[.]{4}"),mot),"($|[ ])/,/^$/p'|");
  else{
    strcat(strcat(strcat(s,"sed -nE '/^[ ]{7}"),mot),"($|[ ])/,/(^$)|(^[ ]{7}");
    strcat(s,strcat(fin,")/p'|tail -r|sed -n '2,$p'|"));
    // les tail -r permettent de supprimer la derniÃ¨re ligne
    // le awk est pour enlever Ã©ventuellement l'avant derniÃ¨re ligne "...."
    strcat(s,"awk '{if(NR>1)print $0;else if(!match($0,/^[.]{4}/))print $0;}'|");
    strcat(s,"tail -r|");
  }

  strcat(s,"sed -E 's/^[.]{4}/    /g'|sed 's/!!!/   /g'");
  strcat(s,"|awk '{n++;print $0;}END{if(!n) ");
  strcat(s,"print\"Erreur : argument incorrect.\\nAide sur ");
  strcat(s,mot);
  strcat(s," non trouvÃ©e.\\n\";}'");
  printf("\n");
  system(s);

  if(j<i) printf("\n");
  free(s);
  free(mot);
  free(fin);
  exit(EXIT_SUCCESS);  
}


string GetArgInc(int *i)
/*
  Retourne ARGV[i], s'il existe, puis incrÃ©mente i.  Si l'argument
  n'existe pas ou si ARGV[i]="?", on affiche l'aide en ligne sur le
  dernier argument demandÃ© et l'on s'arrÃªte.
*/
{
  if(((*i)==ARGC)||(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return ARGV[(*i)++];
}


void NextArg(int *i)
/*
  VÃ©rifie si le prochain argument existe et incrÃ©mente i. S'il
  n'existe pas ou si c'est "?", une aide sur le dernier argument est
  affichÃ©e comme l'aurait fait GetArgInc(). Cette fonction devrait
  Ãªtre appelÃ©e chaque fois que le prochain argument doit Ãªtre testÃ©
  avec: if EQUAL("..."). Si on ne le fait pas, alors c'est un
  "Segmentation fault".
*/
{
  (*i)++;
  GetArgInc(i);
  (*i)--;
  return;
}

void CheckHelp(int *i)
/*
  IncrÃ©mente i puis vÃ©rifie si l'argument est "?". Si tel est le cas,
  une aide est affichÃ©e sur ARGV[i] (avant incrÃ©mentation).  Cette
  fonction devrait Ãªtre typiquement appellÃ©e lorsqu'on vÃ©rifie l'aide
  pour un graphe sans paramÃ¨tre. Sinon c'est fait par GetArgInc().
*/
{
  (*i)++;
  if(((*i)!=ARGC)&&(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return;
}


void Help(int i)
/*
  Affiche:
  - l'aide complÃ¨te si ARGV[i] est "-help" ou "?", ou
  - les paragraphes contenant ARGV[i].
  Puis sort par exit().

  Utilise les commandes: sed, more, awk, sort.
*/
{
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' "); /* filtre l'aide */
  strcat(strcat(s,*ARGV),".c | "); /* enlÃ¨ve 1Ã¨re et derniÃ¨re ligne */
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g' ");
  i++;
  if((i==ARGC)||(ARGV[i][0]=='?'))
    strcat(s,"-e 's/^[.][.][.][.][.]*/    /g'|more"); /* aide complÃ¨te */
  else{
    strcat(s,"|awk 'BEGIN{x=\".....\"}/^[.]{4}./{x=$0} /");
    strcat(strcat(s,ARGV[i]),"/{if(x!=\".....\")print substr(x,5)}'|sort -u");
  }
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void ListGraphs(void)
/*
  Affiche les graphes disponibles et puis termine par exit(). Utilise
  les commandes: sed, grep.
*/
{
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g'| ");
  strcat(s,"grep '^[.][.][.][.][^-.]'| sed 's/^[.][.][.][.]//g'");
  //printf("%s\n",s);
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void Version(void)
/*
  Affiche la version du programme et puis termine par exit(). Utilise
  la commande: sed.
*/
{
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -n 's/.*[-] v\\([0-9]*[.][0-9]*\\) [-].*/\\1/p'");
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void PipeDot(int j)
/*
  GÃ¨re l'option "-format dot<type>".

  Rem: on pourrait utiliser popen() au lieu de rÃ©Ã©crire la ligne de
  commande et de lancer system(). Utilise la commande: dot.
*/
{
  CheckHelp(&j);j--; /* vÃ©rifie l'aide, ici ARGV[j]="dot<type>" */
  string type=strdup(ARGV[j]+3); /* type=<type> */
  ARGV[j][3]='\0'; /* ARGV[j]="dot" plutÃ´t que "dot<type>" */

  /*
    On rÃ©Ã©crit la ligne de commande:
    1. en remplaÃ§ant "-format dot<type>" par "-format dot"
    2. puis en ajoutant Ã  la fin: "| dot -T<type> -K <filter>"
  */
  
  string s=MakeCMD(NULL,0,ARGC); // ne pas libÃ©rer ce pointeur
  strcat(strcat(strcat(strcat(s,"| dot -T"),type)," -K "),DOTFILTER);
  free(type);
  system(s);
  exit(EXIT_SUCCESS);
}


void Visu(int j)
/*
  GÃ¨re l'option "-visu". Utilise la commande: dot.

  On rÃ©Ã©crit la ligne de commande:
  1. en remplaÃ§ant "-visu" par "-format userdot" si le FORMAT
     est "userdot" et par "-format dot" sinon
  2. puis en ajoutant Ã  la fin: "| dot -Tpdf -K <filter> -o g.pdf"

  Si le FORMAT est F_no (c'est le cas si l'on a fait "-check maincc"
  ou "loadc" par exemple), alors il y a un problÃ¨me puisqu'il faut
  qu'un graphe soit gÃ©nÃ©rÃ©.
*/
{
  CheckHelp(&j);j--; /* vÃ©rifie l'aide, ici ARGV[j]="-visu" */
  if(FORMAT==F_no) Erreur(24);

  /* on reconstruit dans s la ligne de commande */
  string s=MakeCMD(NULL,0,j); // ne pas libÃ©rer ce pointeur
  strcat(s,"-format ");
  if(FORMAT==F_userdot) strcat(s,"user");
  strcat(s,"dot ");
  MakeCMD(s,j+1,ARGC);
  strcat(strcat(strcat(s,"| dot -Tpdf -K "),DOTFILTER)," -o "xstr(GRAPH_PDF));
  system(s);
  exit(EXIT_SUCCESS);
}


void Visuh(int j)
/*
  GÃ¨re l'option "-visuh", de maniÃ¨re similaire Ã  "-visu".

  On rÃ©Ã©crit la ligne de commande:
  1. en remplaÃ§ant "-visu" par "-format html"
  2. puis en ajoutant Ã  la fin: "> g.html"

  Si le FORMAT est F_no (c'est le cas si l'on a fait "-check maincc"
  ou "loadc" par exemple), alors il y a un problÃ¨me puisqu'il faut
  qu'un graphe soit gÃ©nÃ©rÃ©.
*/
{
  CheckHelp(&j);j--; /* vÃ©rifie l'aide, ici ARGV[j]="-visuh" */
  if(FORMAT==F_no) Erreur(24);

  /* on reconstruit dans s la ligne de commande */
  string s=MakeCMD(NULL,0,j); // ne pas libÃ©rer ce pointeur
  strcat(s,"-format html ");
  MakeCMD(s,j+1,ARGC);
  strcat(s,"> "xstr(GRAPH_HTML));
  system(s);
  exit(EXIT_SUCCESS);
}


void MainCC(int j)
/*
  GÃ¨re l'option "-maincc".

  On rÃ©Ã©crit la ligne de commande en remplaÃ§ant "-maincc" par "-check
  maincc | ./gengraph load - -fast". Utilise la commande: gengraph.
*/
{
  CheckHelp(&j);j--; /* vÃ©rifie l'aide, ici ARGV[j]="-maincc" */

  /* on reconstruit dans s la nouvelle ligne de commande */
  string s=MakeCMD(NULL,0,j); // ne pas libÃ©rer ce pointeur
  strcat(s,"-check maincc | ./gengraph load - -fast ");
  MakeCMD(s,j+1,ARGC);
  system(s);
  exit(EXIT_SUCCESS);
}


int gsub(string s,string const t,string const r)
/*
  Remplace, dans la chaÃ®ne s, toutes les occurences de t par r et
  renvoie le nombre de remplacements. La chaÃ®ne s est modifiÃ©e
  directement. Il est nÃ©cessaire que s ait suffisamment de place.
*/
{
  int const lr=strlen(r); // longueur de t
  int const lt=strlen(t); // longueur de r
  int const d=lt-lr;

  string p=s+strlen(s)+1; // pointeur sur la fin de s avec son '\0'
  int n=0; // nombre de remplacements effectuÃ©s

  while((s=strstr(s,t))){
    n++; // une occurrence de plus
    p -= d; // met Ã  jour le pointeur de fin de s
    memmove(s+lr,s+lt,p-s); // dÃ©place la fin
    memcpy(s,r,lr); // copie le remplacement
  }

  return n;
}


/***********************************

           OPTIONS -ALGO

***********************************/


void PrintMorphism(string const s,const int* const P,int const n)
/*
  Affiche le tableau P de n Ã©lÃ©ments sous la forme:
  0->i0   1->i1 ... 7->i7
  8->i8   9->i9 ...
  ...
  oÃ¹ i_j=P[i].
*/
{
  int const k=8; /* nombre de "->" affichÃ©s par ligne */
  printf("%s",s); /* normalement printf(s) est ok, mais Warning sur certains systÃ¨mes */
  for(int i=0;i<n;i++)
    printf("%i->%i%s",i,P[i],(((i%k)<(k-1))&&(i<n-1))?"\t":"\n");
  return;
}


/***********************************

           FONCTIONS DE TESTS
              POUR -FILTER

***********************************/


int ftest_minus(graph* const G)
/*
  Retourne VRAI ssi G n'est pas dans la famille F, c'est-Ã -dire G
  n'est isomorphe Ã  aucun graphe de F.
*/
{
  graph* F=MEM(FPARAM,0,graph*);
  if(F==NULL) return (G==NULL);

  int i,*P;
  for(i=0;i<F->f;i++){
    P=Isomorphism(G,F->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }
  
  return 1;
}


int ftest_minus_id(graph* const G)
/*
  Retourne VRAI ssi F2 ne contient aucun graphe d'identifiant Ã©gale Ã 
  celui de G. La complexitÃ© est en O(log|F2|). Il est important que F2
  soit triÃ©e par ordre croissant des ID.
*/
{
  graph* F=MEM(FPARAM,0,graph*);
  if(F==NULL) return 0;
  return (bsearch(&G,F->G,F->f,sizeof(graph*),fcmp_graphid)==NULL);
}


int ftest_unique(graph* const G)
/*
  Retourne VRAI ssi la sous-famille F allant des indices i+1 Ã  F->f ne contient
  pas G, oÃ¹ i=MEM(FPARAM,0,int).

  Effet de bord: MEM(FPARAM,0,int) est incrÃ©mentÃ©.
*/
{
  int i = (MEM(FPARAM,0,int) += 1);
  int *P;

  for(;i<FAMILY->f;i++){
    P=Isomorphism(G,FAMILY->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }

  return 1;
}


int ftest_minor(graph* const G)
/*
  Retourne VRAI ssi H est mineur de G.
*/
{
  int *T=Minor(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_minor_inv(graph* const G)
{
  int *T=Minor(G,MEM(FPARAM,0,graph*));
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_sub(graph* const G)
/*
  Retourne VRAI ssi H est sous-graphe de G avec mÃªme nb de sommets.
*/
{
  graph* C=Subgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_sub_inv(graph* const G)
{
  graph* C=Subgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_isub(graph* const G)
/*
  Retourne VRAI ssi H est sous-graphe induit de G.
*/
{
  int *C=InducedSubgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_isub_inv(graph* const G)
{
  int *C=InducedSubgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_iso(graph* const G)
/*
  Retourne VRAI ssi H est isomorphe Ã  G. Aucun intÃ©rÃªt de faire
  programmer iso-inv.
*/
{
  int *T=Isomorphism(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_id(graph* const G)
{ return InRange(G->id,FPARAM); }


int ftest_vertex(graph* const G)
{ return InRange(G->n,FPARAM); }


int ftest_edge(graph* const G)
{ return InRange(nb_edges(G),FPARAM); }


int ftest_degmax(graph* const G)
{ return InRange(Degree(G,1),FPARAM); }


int ftest_degmin(graph* const G)
{ return InRange(Degree(G,0),FPARAM); }


int ftest_deg(graph* const G)
{ int u,b=1,n=G->n;
  for(u=0;(u<n)&&(b);u++) b=InRange(G->d[u],FPARAM);
  return b;
}


int ftest_degenerate(graph* const G)
{
  int x;
  int *T=Prune(G,&x);
  free(T);
  return InRange(x,FPARAM);
}


int ftest_gcolor(graph* const G)
{
  int *T=Prune(G,NULL);
  int *C=GreedyColor(G,T);
  free(T);
  free(C);
  return InRange(1+G->int1,FPARAM);
}


int ftest_component(graph* const G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM);
}


int ftest_forest(graph* const G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM)&&(nb_edges(G)==G->n-c);
}


int ftest_cutvertex(graph* const G)
{
  param_dfs *p=dfs(G,0,NULL);
  int x=p->na;
  free_param_dfs(p);
  return InRange(x,FPARAM);
}


int ftest_biconnected(graph* const G)
{
  param_dfs *p=dfs(G,0,NULL);
  int b=(p->nc==1)&&(p->na==0)&&(G->n>2);
  free_param_dfs(p);
  return b;
}


int ftest_ps1xx(graph* const G, int version)
{
  path *P=new_path(G,NULL,G->n);
  int v=PS1(G,P,version);
  free_path(P);
  return v;
}


int ftest_ps1(graph* const G) { return ftest_ps1xx(G,0); }
int ftest_ps1b(graph* const G) { return ftest_ps1xx(G,1); }
int ftest_ps1c(graph* const G) { return ftest_ps1xx(G,2); }
int ftest_ps1x(graph* const G) { return ftest_ps1xx(G,3); }

int ftest_radius(graph* const G)
{
  param_bfs *p;
  int const n=G->n;
  int x=n,u;

  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(x>=0){
      if(p->n<n) x=-1;
      else x=min(x,p->radius);
    }
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_girth(graph* const G)
{
  param_bfs *p;
  int x=1+G->n,u;
  int const n=G->n;
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(p->cycle>0) x=min(x,p->cycle);
    free_param_bfs(p);
  }
  if(x>n) x=-1;
  return InRange(x,FPARAM);
}


int ftest_diameter(graph* const G)
{
  param_bfs *p;
  int x=-1,u;
  int const n=G->n;
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(p->n==n) x=max(x,p->radius);
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_hyper(graph* const G)
{
  param_bfs *p;
  int const n=G->n;
  int u,v,x,y,d1,d2,d3,h=0;
  NALLOC(int*,D,n);

  /* calcule la matrice de distance D[][] */
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL); // Dijkstra
    D[u]=p->D;
    if(p->n<n){ /* remplace -1 par +âˆž */
      for(v=0;v<n;v++) if(p->D[v]<0) p->D[v]=INT_MAX;
    }
    p->D=NULL;
    free(p); /* efface p, mais pas p->D */
  }

  /* pour tous les quadruplets {u,v,x,y} */
  for(u=0;u<n;u++)
    for(v=u+1;v<n;v++)
      for(x=v+1;x<n;x++)
	for(y=x+1;y<n;y++){
	  d1=D[u][v]+D[x][y];
	  d2=D[u][x]+D[v][y];
	  d3=D[u][y]+D[v][x];
	  if(d1<d2) SWAP(d1,d2);
	  if(d1<d3) SWAP(d1,d3);
	  if(d2<d3) d2=d3; /* on se fiche de d3 */
	  if(d1-d2>h) h=d1-d2;
	}

  FREE2(D,n); /* efface la matrice de distances */
  if(h==INT_MAX) h=-1; /* cela ne peut jamais arriver */
  return InRange(h,FPARAM);
}


int ftest_tw(graph* const G)
{ return InRange(Treewidth(G,1),FPARAM); }


int ftest_tw2(graph* const G)
{ return (Treewidth(G,0)<=2); }


int ftest_rename(graph* const G)
{
  G->id=SHIFT++;
  return 1;
}


graph* Filter(const graph* const F,test* const f,int const code){
/*
  Etant donnÃ©e une famille de graphes et une fonction de test f,
  renvoie une sous-famille de graphes G de F telle que f(G) est
  vraie (si code=0) ou faux (si code=1). Attention! si on libÃ¨re F,
  cela dÃ©truit la sous-famille renvoyÃ©e.

  Effet de bord: si PVALUE est vrai, alors dans les graphes filtrÃ©s
  G on met dans G->int1 la valeur du paramÃ¨tre, CVALUE.
*/
  if((F==NULL)||(F->f==0)) return NULL;
  int i,j,n=F->f;

  graph* const R=new_graph(0);
  ALLOC(R->G,n); /* a priori R est de mÃªme taille que F */

  for(i=j=0;i<n;i++){
    if(f(F->G[i])^code){
      R->G[j++]=F->G[i];
      if(PVALUE) F->G[i]->int1=CVALUE;
    }
  }

  REALLOC(R->G,j);
  R->f=j;
  return R;
}


graph* Graph2Family(graph* const G)
/*
  Renvoie une famille composÃ©e d'un seul graphe G.
  Effet de bord: met G->id=0.
*/
{
  graph* const F=new_graph(0);
  ALLOC(F->G,1);
  F->G[0]=G;
  F->f=1;
  G->id=0;
  return F;
}


void ApplyFilter(int const code,int const index)
/*
  Applique le filtre FTEST (avec le code=0 ou 1) Ã  la famille de
  graphes FAMILY (voir un graphe seul), et affiche la famille
  rÃ©sultante. On affiche aussi le nombre de graphes obtenus et la
  ligne de commande (sous forme de commentaire). Si index>=0, alors
  ARGV[index] donne le paramÃ¨tre.

  Effet de bord: FAMILY est libÃ©rÃ©e.
*/
{
  graph* R;
  int i;

  if(FAMILY->f==0) FAMILY=Graph2Family(FAMILY); /* transforme en famille si graphe simple */
  R=Filter(FAMILY,FTEST,code); /* calcule la sous-famille */

  printf("// #graphs: %i\n// generated by:",(R==NULL)?0:R->f);
  for(i=0;i<ARGC;i++) printf(" %s",ARGV[i]);
  printf("\n");
  if((index>=0)&&(PVALUE)) /* on affiche la valeur de chaque graphe */
    for(i=0;i<R->f;i++)
      printf("[%i] %s: %i\n",R->G[i]->id,ARGV[index],R->G[i]->int1);
  else PrintGraph(R); /* ou bien on affiche le graphe */

  /* ! aux free() de famille de graphes ! */
  free_graph(FAMILY); /* libÃ¨re en premier la famille */
  free(R->G); /* libÃ¨re la sous-famille */
  free(R);
  return;
}


void RS_Start(string const nom,string const type,graph* const G)
/*
  Partie commune Ã  tous les schÃ©mas de routage.
  - nom: est le nom du schÃ©ma
  - type: sa catÃ©gorie (name-independent ...)
  - G: le graphe sur lequel le schÃ©ma doit Ãªtre appliquÃ©
*/
{
  printf("\nROUTING SCHEME\n");
  BARRE;
  printf("- name: %s\n",nom);
  printf("- type: %s\n",type);
  printf("- command: %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libÃ©rer ce pointeur
  printf("- date: %s\n",DateHeure());
  printf("- seed: %u\n",SEED);
  printf("- time for loading/generating the graph: %s\n",TopChrono(1));
  if(nb_edges(G)<1) Erreur(36); // il faut au moins 1 arÃªte
  param_dfs *X=dfs(G,0,NULL);
  int c=X->nc;
  free_param_dfs(X);
  if(c!=1) Erreur(11); // il doit Ãªtre connexe
  printf("- checking connectivity: Ok (%s)\n",TopChrono(1));
  int *R=SortGraph(GF,1);
  printf("- checking graph type: simple and undirected (%s)\n",TopChrono(1));
  printf("- #nodes: %i\n",G->n);
  printf("- #edges: %i\n",G->m); // G->m est a jour
  printf("- average degree: %.2lf\n",(double)(G->m<<1)/(double)G->n);
  if(!R[6]) Erreur(31); /* le graphe doit Ãªtre simple */
  string s;
  char t[128];
  s="?";
  if(HASH==H_MIX)     s="mix";
  if(HASH==H_PRIME)   s="prime";
  if(HASH==H_SHUFFLE) s="shuffle";
  if(HASH==H_MOD)     s="mod";
  printf("- hash: %s\n",s);
  s="?";
  if(SCENARIO.mode==SC_NONE)   s="none";
  if(SCENARIO.mode==SC_ALL)    s="all";
  if(SCENARIO.mode==SC_NPAIRS) s="n pairs";
  if(SCENARIO.mode==SC_PAIR)   s="pair";
  if(SCENARIO.mode==SC_EDGES)  s="edges";
  if(SCENARIO.mode==SC_ONE)    s="one-to-all";
  if(SCENARIO.mode==SC_UV)     s="u->v";
  if(SCENARIO.mode==SC_UNTIL){ sprintf(t,"until stretch â‰¥ %g",SCENARIO.stretch);s=t;}
  printf("- scenario: %s\n",s);
  return;
}


/***********************************

               MAIN

***********************************/


int main(int const argc,string argv[]){
  
  DEBUG(

	// prÃ©caution utile car sinon on peut chercher longtemps une
	// erreur alors qu'elle provient simplement d'effets de bord
	// liÃ©e au dÃ©bugage potentiellement Ã©trange dans certaines
	// situations
	
	fprintf(stderr,"//////////////////////////////////////\n");
	fprintf(stderr,"//!!! WARNING !!! DEBUGGING MODE !!!//\n");
	fprintf(stderr,"//////////////////////////////////////\n");
	
	);

  ARGC=argc;
  ARGV=argv;

  if(ARGC==1) Help(0);   /* aide si aucun argument */

  /* initialisations */

  TopChrono(0); /* intialise tous les chronomÃ¨tres internes */
  VIDE(PARAM_PAL);
  SEED=arc4random()&0xFFFF; // initialisation vraiment alÃ©atoire sur 16 bits
  srandom(SEED); // seed sur 16 bits plus facile Ã  recopier Ã  la main
  query *Q=new_query(); /* pour les entrÃ©es/sorties d'une fonction d'adjacence */
  int i,j,k,t;

  /******************************************************************

                 ANALYSE DE LA LIGNE DE COMMANDE

    o Il faut Ã©viter de faire des traitements trop coÃ»teux en
      temps/mÃ©moire dans l'analyse de la ligne de commande car on peut
      Ãªtre amenÃ© Ã  la refaire une deuxiÃ¨me fois Ã  cause des alias, les
      options comme -visu, -maincc, -format dot<type> qui causent la
      rÃ©Ã©criture puis la rÃ©-analyse de la nouvelle ligne de commande.

    o Il faut Ã©viter d'utiliser random() dans l'analyse de la ligne de
      commande, car si l'option -seed est prÃ©sente, le comportement
      dÃ©pendra de sa position dans la ligne de commande.  Cependant,
      dans certains cas comme "caterpillar" cela est inÃ©vitable.

    o Il y a essentiellement quatre cas de figure pour la lecture
      d'arguments ou d'options de la ligne de commande. Il faut
      respecter la structure suivante pour que l'aide en ligne
      fonctionne bien et pour Ã©viter les "Segmentation fault" lors de
      lecture d'argument. Cela arrive dÃ¨s qu'on essaye de faire
      EQUAL("...") avec une valeur de i en dehors de ARGV[]. En
      gÃ©nÃ©ral, il faut Ã©viter de faire un "i++" ou un "GetArgInc(&i)"
      suivit d'un EQUAL("...").

      1. option ou graphe avec au moins un argument:

         if EQUAL("kpage"){ i++;
           Q->adj=kpage;
           Q->param[0]=STRTOI(GetArgInc(&i));
           Q->param[1]=STRTOI(GetArgInc(&i));
           ...
           goto fin;
         }

      2. option ou graphe sans argument attendu:

         if EQUAL("tutte"){ CheckHelp(&i);
           Q->adj=tutte;
           goto fin;
         }

      3. option avec une ou plusieurs variantes (de type 1, 2 ou 3):

         if EQUAL("-xy"){ NextArg(&i);
           if EQUAL("load"){ ...; goto fin; }
           if EQUAL("unif"){ ...; goto fin; }
           ...
           Erreur(...);
         }

      4. option avec plusieurs variantes et argument:

         if EQUAL("-xy"){ i++;
           Q->param[0]=STRTOI(GetArgInc(&i));
           i--;NextArg(&i);
           if EQUAL("load"){ ...; goto fin; }
           if EQUAL("unif"){ ...; goto fin; }
           ...
           Erreur(...);
         }

  ******************************************************************/

  i=1; /* on dÃ©marre avec le 1er argument */
  while(i<ARGC){

    /*****************/
    /* les aides ... */
    /*****************/

    if EQUAL("-help"){ CheckHelp(&i); i--; }
    if(EQUAL("-help")||EQUAL("?")) Help(i);
    if EQUAL("-list"){ CheckHelp(&i); ListGraphs(); }
    if EQUAL("-version"){ CheckHelp(&i); Version(); }
    
    j=i; /* mÃ©morise i pour savoir si on a rÃ©ussit Ã  lire au moins une
	    option ou un graphe */

    /***********************/
    /* les options -xxx... */
    /***********************/

    if EQUAL("-visu") Visu(i);     /* se termine par system() & exit() */
    if EQUAL("-visuh") Visuh(i);   /* se termine par system() & exit() */
    if EQUAL("-maincc") MainCC(i); /* se termine par system() & exit() */
    if EQUAL("-not"){ Q->not=!Q->not; goto param0; }
    if EQUAL("-permute"){ PERMUTE=1; goto param0; }
    if EQUAL("-undirected"){ Q->directed=Q->loop=0; goto param0; }
    if EQUAL("-directed")  { Q->directed=Q->loop=1; goto param0; }
    if EQUAL("-header"){ HEADER=1; goto param0; }
    if EQUAL("-vsize"){ VSIZE=1; goto param0; }
    if EQUAL("-fast"){ FAST=1; goto param0; }
    if EQUAL("-loop"){ i++;
	Q->loop=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-seed"){ i++;
	srandom(SEED=STRTOI(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-width"){ i++;
	WIDTH=max(STRTOI(GetArgInc(&i)),0);
	goto fin;
      }
    if EQUAL("-shift"){ i++;
	SHIFT=STRTOI(GetArgInc(&i));
	if(SHIFT<0) Erreur(6);
	goto fin;
      }
    if EQUAL("-pos"){ i++;
	POS=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-label"){ i++;
	LABEL=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-dot"){ NextArg(&i);
	if EQUAL("len"){ i++;
	    LEN=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("filter"){ i++;
	    DOTFILTER=GetArgInc(&i); /* pointe sur le nom du filtre */
	    goto fin;
	  }
	if EQUAL("scale"){ i++;
	    DOTSCALE=GetArgInc(&i); /* pointe sur la/les valeurs */
	    goto fin;
	  }
	Erreur(49);
	goto fin;
      }
    if EQUAL("-norm"){ NextArg(&i);
	NORM=NORM_FAIL;
	if EQUAL("L1")    NORM=NORM_L1;
	if EQUAL("L2")    NORM=NORM_L2;
	if EQUAL("Lmax")  NORM=NORM_LMAX;
	if EQUAL("Lmin")  NORM=NORM_LMIN;
	if EQUAL("hyper") NORM=NORM_HYPER;
	if EQUAL("poly"){ i++;
	    NORM=NORM_POLY; 
	    NORM_poly=STRTOI(GetArgInc(&i));
	    if(NORM_poly<3) NORM=NORM_L2;
	    goto fin;
	  }
	if(NORM==NORM_FAIL) Erreur(43);
	i++;
	goto fin;
      }
    if EQUAL("-delv"){ i++;
	DELV=STRTOD(GetArgInc(&i));
	if(DELV>1) DELV=1;
	goto fin;
      }
    if EQUAL("-dele"){ i++;
	DELE=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-redirect"){ i++;
	REDIRECT=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-caption"){ i++;
	string c=GetArgInc(&i); // lecture de la lÃ©gende
	if(CAPTION) free(CAPTION); // si CAPTION dÃ©jÃ  allouÃ©e
	string s=strdup(c);
	k=gsub(s,"%SEED","%u");
	if(k>1) Erreur(35);
	if(k==1) asprintf(&CAPTION,s,SEED);
	if(k==0) CAPTION=s;
	goto fin;
      }
    if EQUAL("-variant"){ i++;
	VARIANT=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-vcolor"){ NextArg(&i);
	/* bits 0-2: codage de la fonction de couleur=1,2,3,4,5
	   bit 3: "pal"
	   bit 4: "list" */
	if EQUAL("deg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7; /* efface les 3 derniers bits */
	    VCOLOR |= 1; goto fin; }
	if EQUAL("degr"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 2; goto fin; }
	if EQUAL("degm"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 3; CHECK=max(CHECK,CHECK_ON); goto fin; }
	if EQUAL("randg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 4; CHECK=max(CHECK,CHECK_ON); goto fin; }
	if EQUAL("kcolor"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 5; CHECK=max(CHECK,CHECK_ON);
	    if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE);
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i)); /* CPARAM[0]=ARGV[i] */
	    goto fin; }
	if EQUAL("pal"){ NextArg(&i); /* teste si arg aprÃ¨s pal existe bien */
	    VCOLOR |= 0x8; /* set bit-3 */
	    if(strlen(ARGV[i])>=(int)(sizeof(PARAM_PAL)/sizeof(char))) Erreur(20);
	    strcpy(PARAM_PAL,GetArgInc(&i)); /* PARAM_PAL=ARGV[i] */
	    goto fin; }
	if EQUAL("list"){ i++;
	    VCOLOR |= 0x10;
	    CHECK=max(CHECK,CHECK_ON);
	    FORMAT = F_dot;
	    goto fin; }
	Erreur(9);
      }
    if EQUAL("-format"){ NextArg(&i);
	FORMAT=-1; /* sentiennelle pour savoir si on a trouvÃ© le FORMAT */
	if EQUAL("standard") FORMAT=F_standard;
	if EQUAL("xy")       FORMAT=F_xy;
	if EQUAL("no")       FORMAT=F_no;
	if EQUAL("userdot"){ FORMAT=F_userdot; WIDTH=1; }
	if EQUAL("html")   { FORMAT=F_html; WIDTH=1; }
	if EQUAL("matrix") { FORMAT=F_matrix; CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("smatrix"){ FORMAT=F_smatrix;CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("list")   { FORMAT=F_list;   CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("vertex"){ i++;
	    FORMAT=F_list;
	    CHECK=max(CHECK,CHECK_ON);
	    VERTEX0=STRTOI(GetArgInc(&i));
	    i--;
	  }
	if PREFIX("dot"){ /* si "dot" ou "dot<type>" */
	    if EQUAL("dot") FORMAT=F_dot; /* si "dot" seul */
	    else PipeDot(i); /* se termine par system() & exit() */
	  }
	if(FORMAT<0) Erreur(5); /* le format n'a pas Ã©tÃ© trouvÃ© */
	i++;
	goto fin;
      }
    if EQUAL("-xy"){ NextArg(&i);
	POS=1; /* il faut POS=1 */
	if EQUAL("unif"){ CheckHelp(&i);
	    XYtype=XY_UNIF;
	    goto fin;
	  }
	if EQUAL("load"){ i++;
	    XYtype=XY_FILE;
	    FILEXY=GetArgInc(&i); /* pointe sur le nom du fichier */
	    goto fin;
	  }
	if EQUAL("noise"){ i++;
	    XYnoiser=STRTOD(GetArgInc(&i));
	    XYnoisep=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("box"){ i++;
	    BOXX=STRTOD(GetArgInc(&i));
	    BOXY=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("seed"){ i++;
	    XYseedk=STRTOI(GetArgInc(&i));
	    XYpower=STRTOD(GetArgInc(&i));
	    XYtype=XY_PLAW;
	    goto fin;
	  }
	if EQUAL("hyper"){ i++;
	    XYpower=STRTOD(GetArgInc(&i));
	    XYtype=XY_HYPER;
	    goto fin;
	  }
	if EQUAL("round"){ i++;
	    ROUND=min(STRTOI(GetArgInc(&i)),DBL_DIG);
	    goto fin;
	  }
	if EQUAL("mesh"){ i++;
	    ROUND=0; /* a priori coordonnÃ©es entiÃ¨res dans ce cas */
	    XYtype=XY_MESH;
	    Xmesh=STRTOI(GetArgInc(&i));
	    Ymesh=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("grid"){ i++;
	    XYgrid=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("ratio"){ i++;
	    XYratio=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("vsize"){ i++;
	    XYvsize=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("polygon"){ i++;
	    XYtype=XY_RPOLY;
	    XYpoly=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("surface"){ NextArg(&i);
	    string s=NULL; // s=signature
	    if EQUAL("square")     s="bb";
	    if EQUAL("plane")      s="bb";
	    if EQUAL("cylinder")   s="hb";
	    if EQUAL("mobius")     s="cb";
	    if EQUAL("torus")      s="hh";
	    if EQUAL("klein")      s="ch";
	    if EQUAL("projective") s="cc";
	    if(s) i++; else s=GetArgInc(&i);
	    int g=0,v; // g=taille de la signature
	    while(s[g]){
	      v=2; // sentinelle
	      if(s[g]=='b') v=0;  // border
	      if(s[g]=='h') v=+1; // handle
	      if(s[g]=='c') v=-1; // crosscap
	      if(v==2) Erreur(45); // caractÃ¨re non reconnu
	      XYsurface[g++]=v;
	      if(g>SURFACEMAX) Erreur(45); // trop grand
	    }
	    if(g%2) Erreur(45); // problÃ¨me si g est impaire
	    XYsurfacesize=g;
	    XYtype=XY_RPOLY; // gÃ©nÃ©ration dans un polygone
	    XYpoly=2*g; // nombre de cotÃ©s du polygone = genre
	    XYratio=1; // polygone rÃ©gulier
	    goto fin;
	  }
	if EQUAL("zero"){ CheckHelp(&i);
	    XYzero=1;
	    goto fin;
	  }
	if EQUAL("border"){ CheckHelp(&i);
	    XYborder=1;
	    goto fin;
	  }
	if EQUAL("permutation"){ CheckHelp(&i);
	    ROUND=0; /* a priori coordonnÃ©es entiÃ¨res dans ce cas */
	    XYtype=XY_PERM;
	    goto fin;
	  }
	if EQUAL("circle"){ CheckHelp(&i);
	    XYtype=XY_CIRCLE;
	    goto fin;
	  }
	if EQUAL("cycle"){ CheckHelp(&i);
	    XYtype=XY_CYCLE;
	    goto fin;
	  }
	if EQUAL("disk") { CheckHelp(&i);
	    XYtype=XY_DISK;
	    goto fin;
	  }
	if EQUAL("convex"){ CheckHelp(&i);
	    XYtype=XY_CONVEX;
	    goto fin;
	  }
	if EQUAL("convex2"){ CheckHelp(&i);
	    XYtype=XY_CONVEX2;
	    goto fin;
	  }
	if EQUAL("unique"){ CheckHelp(&i);
	    XYunique=1;
	    goto fin;
	  }
	Erreur(1); /* l'option aprÃ¨s -xy n'a pas Ã©tÃ© trouvÃ©e */
      }
    if EQUAL("-filter"){ NextArg(&i);
	FAMILY=File2Graph(ARGV[i],2); /* lit une famille ou un graphe */
	if(FPARAM==NULL) ALLOC(FPARAM,PARAMSIZE);
	PVALUE=0; /* par dÃ©faut, on affiche pas "value" mais les graphes */
	NextArg(&i);
	if EQUAL("not"){ k=1; NextArg(&i); } else k=0;
	 /* vÃ©rifie s'il y a bien un autre argument */
	if EQUAL("rename"){
	    FTEST=ftest_rename;
	    i++;
	    SHIFT=STRTOI(GetArgInc(&i));
	    ApplyFilter(0,-1);
	    goto fin;
	  }
	if EQUAL("biconnected"){
	    FTEST=ftest_biconnected;
	  filter0:
	    i++;
	    ApplyFilter(k,-1);
	    goto fin;
	  }
	if EQUAL("id"){
	    FTEST=ftest_id;
	  filter1:
	    i++;
	    ReadRange(GetArgInc(&i),FPARAM);
	    ApplyFilter(k,i-2);
	    goto fin;
	  }
	int c; /* code pour File2Graph() */
	if EQUAL("minor"){
	    FTEST=ftest_minor;
	    c=34; /* dÃ©tection du shift et charge toujours un graphe */
	  filter2:
	    i++;
	    MEM(FPARAM,0,graph*)=File2Graph(GetArgInc(&i),c);
	    ApplyFilter(k,-1);
	    free_graph(MEM(FPARAM,0,graph*));
	    goto fin;
	  }

	if EQUAL("ps1x"){ i++;
	    FTEST=ftest_ps1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(c=MEM(CPARAM,0,int);c>=1;c--){ /* met les arguments Ã  l'envers */
	      MEM(CPARAM,(2*c-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*c)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    i--;
	    goto filter0;
	  }
	if EQUAL("ps1"){ FTEST=ftest_ps1; goto filter0; }
	if EQUAL("ps1b"){ FTEST=ftest_ps1b; goto filter0; }
	if EQUAL("ps1c"){ FTEST=ftest_ps1c; goto filter0; }
	if EQUAL("tw2"){ FTEST=ftest_tw2; goto filter0; }
	if EQUAL("unique"){ FTEST=ftest_unique; MEM(FPARAM,0,int)=0; goto filter0; }
	if EQUAL("vertex"){ FTEST=ftest_vertex; goto filter1; }
	if(EQUAL("edge")||EQUAL("edges")){ FTEST=ftest_edge; goto filter1; }
	if EQUAL("deg"){ FTEST=ftest_deg; goto filter1; }
	if EQUAL("degenerate"){ FTEST=ftest_degenerate; goto filter1; }
	if EQUAL("degmax"){ FTEST=ftest_degmax; goto filter1; }
	if EQUAL("degmin"){ FTEST=ftest_degmin; goto filter1; }
	if EQUAL("gcolor"){ FTEST=ftest_gcolor; goto filter1; }
	if EQUAL("component"){ FTEST=ftest_component; goto filter1; }
	if EQUAL("cut-vertex"){ FTEST=ftest_cutvertex; goto filter1; }
	if EQUAL("radius"){ FTEST=ftest_radius; goto filter1; }
	if EQUAL("girth"){ FTEST=ftest_girth; goto filter1; }
	if EQUAL("diameter"){ FTEST=ftest_diameter; goto filter1; }
	if EQUAL("hyper"){ FTEST=ftest_hyper; goto filter1; }
	if EQUAL("tw"){ FTEST=ftest_tw; goto filter1; }
	if EQUAL("forest"){ FTEST=ftest_forest; goto filter1; }
	if EQUAL("minor-inv"){ FTEST=ftest_minor_inv; c=34; goto filter2; }
	if EQUAL("sub"){ FTEST=ftest_sub; c=34; goto filter2; }
	if EQUAL("sub-inv"){ FTEST=ftest_sub_inv; c=34; goto filter2; }
	if EQUAL("isub"){ FTEST=ftest_isub; c=34; goto filter2; }
	if EQUAL("isub-inv"){ FTEST=ftest_isub_inv; c=34; goto filter2; }
	if EQUAL("iso"){ FTEST=ftest_iso; c=34; goto filter2; }
	if EQUAL("minus"){ FTEST=ftest_minus; c=2; goto filter2; }
	if EQUAL("minus-id"){ FTEST=ftest_minus_id; c=10; goto filter2; }
	/* alias */
	if EQUAL("connected"){ FTEST=ftest_component;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("bipartite"){ FTEST=ftest_gcolor;
	    ReadRange("<3",FPARAM); goto filter0; }
	if EQUAL("isforest"){ FTEST=ftest_forest;
	    ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("istree"){ FTEST=ftest_forest;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("cycle"){ FTEST=ftest_forest;
	    k=1-k; ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("all"){ FTEST=ftest_vertex;
	    ReadRange("t",FPARAM); goto filter0; }
	Erreur(14); /* l'option aprÃ¨s -filter n'a pas Ã©tÃ© trouvÃ©e */
      }
    if EQUAL("-check"){ NextArg(&i);
	if(CHECK>CHECK_ON) Erreur(27); /* on ne devrait jamais avoir deux fois -check */
	if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE); /* alloue les paramÃ¨tres */
	if EQUAL("bfs"){
	    CHECK=CHECK_BFS;
	  check0:
	    i++;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("iso"){
	    CHECK=CHECK_ISO;
	  check1:
	    MEM(CPARAM,0,int)=++i;
	    GetArgInc(&i); /* pour vÃ©rifier si i existe */
	    goto fin;
	  }
	if(EQUAL("deg")||EQUAL("edge")||EQUAL("edges")){
	  CHECK=CHECK_DEG;
	check_fin:
	  i++;
	  goto fin;
	}
	if EQUAL("paths"){ i++;
	    CHECK=CHECK_PATHS;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    MEM(CPARAM,sizeof(int),int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("ps1x"){ i++;
	    CHECK=CHECK_PS1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(k=MEM(CPARAM,0,int);k>=1;k--){ /* met les arguments Ã  l'envers */
	      MEM(CPARAM,(2*k-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*k)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    goto fin;
	  }
	if EQUAL("routing"){ NextArg(&i);
	    FORMAT=F_no; /* pour tous les routing schemes */
	    if EQUAL("hash"){ NextArg(&i);
		for(;;){ /* pour faire break (importants Ã  cause des i++) */
		  if EQUAL("mix")    { NextArg(&i); HASH=H_MIX; break; }
		  if EQUAL("prime")  { NextArg(&i); HASH=H_PRIME; break; }
		  if EQUAL("shuffle"){ NextArg(&i); HASH=H_SHUFFLE; break; }
		  if EQUAL("mod")    { NextArg(&i); HASH=H_MOD; break; }
		  Erreur(40); /* option aprÃ¨s "hash" non trouvÃ© */
		}
	      }
	    SCENARIO.mode=SC_NONE; /* par dÃ©faut aucun scenario */
	    SCENARIO.dist=1; /* par dÃ©faut on stocke les distance */
	    if EQUAL("scenario"){ NextArg(&i);
		if EQUAL("nomem"){ NextArg(&i); SCENARIO.dist=0; }
		for(;;){ /* pour faire break (importants Ã  cause des i++) */
		  if EQUAL("none")  { NextArg(&i); SCENARIO.mode=SC_NONE; break; }
		  if EQUAL("all")   { NextArg(&i); SCENARIO.mode=SC_ALL; break; }
		  if EQUAL("npairs"){ NextArg(&i); SCENARIO.mode=SC_NPAIRS; break; }
		  if EQUAL("edges") { NextArg(&i); SCENARIO.mode=SC_EDGES; break; }
		  if EQUAL("one")   { i++;
		      SCENARIO.mode=SC_ONE;
		      SCENARIO.u=STRTOI(GetArgInc(&i));
		      i--;NextArg(&i);
		      break;
		    }
		  if EQUAL("until"){ i++;
		      SCENARIO.mode=SC_UNTIL;
		      SCENARIO.stretch=STRTOD(GetArgInc(&i));
		      i--;NextArg(&i);
		      break;
		    }
		  if EQUAL("pair"){ i++;
		      long p=STRTOL(GetArgInc(&i));
		      SCENARIO.u=abs((int)p);
		      if(p<0L){
			if(-p>(long)INT_MAX) Erreur(46);
			SCENARIO.mode=SC_PAIR;
		      }else{
			SCENARIO.mode=SC_UV;
			SCENARIO.v=STRTOI(GetArgInc(&i));
		      }
		      i--;NextArg(&i);
		      break;
		    }
		  Erreur(41);
		}
	      }
	    if EQUAL("cluster"){ CHECK=CHECK_RS_CLUSTER; goto check0; }
	    if EQUAL("dcr"){ CHECK=CHECK_RS_DCR; goto check0; }
	    if EQUAL("agmnt"){ CHECK=CHECK_RS_AGMNT; goto check0; }
	    if EQUAL("bc"){ CHECK=CHECK_RS_BC; goto check0; }
	    if EQUAL("hdlbr"){ CHECK=CHECK_RS_HDLBR; goto check0; }
	    if EQUAL("tzrplg"){ i++;
		CHECK=CHECK_RS_TZRPLG;
		MEM(CPARAM,0,double)=STRTOD(GetArgInc(&i));
		goto fin;
	      }
	    Erreur(39);
	  }
	if EQUAL("dfs"){ CHECK=CHECK_DFS; goto check0; }
	if EQUAL("bellman"){ CHECK=CHECK_BELLMAN; goto check0; }
	if EQUAL("kcolor"){ CHECK=CHECK_KCOLOR; goto check0; }
	if EQUAL("kcolorsat"){ CHECK=CHECK_KCOLORSAT; FORMAT=F_no; goto check0; }
	if EQUAL("kindepsat"){ CHECK=CHECK_KINDEPSAT; FORMAT=F_no; goto check0; }
	if EQUAL("subdiv"){ CHECK=CHECK_SUBDIV; FORMAT=F_no; goto check0; }
	//
	if EQUAL("sub"){ CHECK=CHECK_SUB; goto check1; }
	if EQUAL("isub"){ CHECK=CHECK_ISUB; goto check1; }
	if EQUAL("minor"){ CHECK=CHECK_MINOR; goto check1; }
	//
	if EQUAL("degenerate"){ CHECK=CHECK_DEGENERATE; goto check_fin; }
	if EQUAL("gcolor"){ CHECK=CHECK_GCOLOR; goto check_fin; }
	if EQUAL("ps1"){ CHECK=CHECK_PS1; goto check_fin; }
	if EQUAL("ps1b"){ CHECK=CHECK_PS1b; goto check_fin; }
	if EQUAL("ps1c"){ CHECK=CHECK_PS1c; goto check_fin; }
	if EQUAL("twdeg"){ CHECK=CHECK_TWDEG; goto check_fin; }
	if EQUAL("tw"){ CHECK=CHECK_TW; goto check_fin; }
	if EQUAL("girth"){ CHECK=CHECK_GIRTH; goto check_fin; }
	if EQUAL("info"){ CHECK=CHECK_INFO; FORMAT=F_no; goto check_fin; }
	if(EQUAL("simplify")){ CHECK=CHECK_SIMPLIFY; FORMAT=F_no; goto check_fin; }
	if EQUAL("stretch"){ CHECK=CHECK_STRETCH; goto check_fin; }
	if EQUAL("maincc"){ CHECK=CHECK_MAINCC; FORMAT=F_no; goto check_fin; }
	if(EQUAL("ncc")||EQUAL("connected")){ CHECK=CHECK_NCC; goto check_fin; }
	if EQUAL("diameter"){ CHECK=CHECK_DIAMETER; goto check_fin; }
	if EQUAL("volm"){ CHECK=CHECK_VOLM; goto check_fin; }
	if EQUAL("radius"){ CHECK=CHECK_RADIUS; goto check_fin; }
	Erreur(12); /* l'option aprÃ¨s -check n'a pas Ã©tÃ© trouvÃ©e */
      }

    /*******************/
    /* graphes de base */
    /*******************/

    if EQUAL("tutte"){
	Q->adj=tutte;
      param0:
	CheckHelp(&i); /* CheckHelp() au lieu de i++ car aucun paramÃ¨tre */
	goto fin;
      }
    if EQUAL("prime"){
	Q->adj=prime;
      param1:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("arboricity"){
	Q->adj=arboricity;
      param2:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("sat"){
	Q->adj=sat;
      param3:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ringarytree"){
	Q->adj=ringarytree;
      param4:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("unok"){
	Q->adj=unok; POS=1;
      //param5: // non utilisÃ©
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=STRTOI(GetArgInc(&i));
	Q->param[4]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("grid"){
	Q->adj=grid;
      param_grid: // "p_1 ... p_n ." -> param[0]=n, param[1..n] = p_i
	k=++i; while(strcmp(GetArgInc(&i),".")); // NB: i-k = taille sÃ©quence avec "."
	if(i-k>PARMAX) REALLOC(Q->param,i-k); // extension du tableau
	for(t=k;t<i;t++) Q->param[t-k+1]=STRTOI(ARGV[t]); // lit les valeurs
	Q->param[0]=i-k-1; // nombre de valeurs sans le "."
	goto fin;
      }
    if EQUAL("hypercube"){
	Q->adj=grid;
      param_hypercube:
	i++;
	int d=STRTOI(GetArgInc(&i));
	if(d+1>PARMAX){ free(Q->param); ALLOC(Q->param,d+1); }
	for(k=1;k<=d;k++) Q->param[k]=2;
	Q->param[0]=d;
	goto fin;
      }
    if EQUAL("udg"){ i++;
	Q->adj=udg; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0) Q->dparam[0]=sqrt(log((double)Q->param[0])/(double)Q->param[0]);
	/* Threshold thÃ©orique rc (cf. [EMY07], Theorem 6 avec d=p=2):
	   pour n=10,    rc=0.4798
	   pour n=100,   rc=0.2145
	   pour n=1000,  rc=0.08311
	   pour n=10000, rc=0.03034
	 */
	goto fin;
      }
    if EQUAL("squashed"){ i++;
	Q->adj=squashed;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0) Q->dparam[0]=1.0/3; /* proba par dÃ©faut */
	Q->dparam[0]=CheckProba(Q->dparam[0]);
	goto fin;
      }
    if EQUAL("rig"){ i++;
	Q->adj=rig;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[0]=max(Q->param[0],1); /* n>0 */
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[1]=max(Q->param[1],1); /* k>0 */
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0){
	  Q->dparam[0]=log((double)Q->param[0])/(double)Q->param[1];
	  if(Q->param[1]>Q->param[0]) Q->dparam[0]=sqrt(Q->dparam[0]/(double)Q->param[0]);
	}
	Q->dparam[0]=CheckProba(Q->dparam[0]); /* proba */
	goto fin;
      }
    if EQUAL("rplg"){ i++;
	Q->adj=rplg;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("thetagone"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if(EQUAL("deltohedron")||EQUAL("trapezohedron")){ i++;
	Q->adj=deltohedron;
	Q->param[0]=(STRTOI(GetArgInc(&i))<<1);
	goto fin;
      }
    if EQUAL("load"){
      param_load:
	i++;
	Q->adj=load;
	Q->sparam=strdup(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ngon2"){ i++;
	POS=1;
	Q->adj=ngon2;
	Q->sparam=strdup(GetArgInc(&i));
	Q->param[0]=strlen(Q->sparam)/2+1;
	if(t>PARMAX){ free(Q->param); ALLOC(Q->param,t); }
	for(k=1;k<=Q->param[0];k++) Q->param[k]=STRTOI(GetArgInc(&i)); // lit les valeurs
	goto fin;
      }
    if EQUAL("rlt"){ i++;
	Q->adj=rlt;
	POS=1; ROUND=0; XYtype=XY_MESH;
	Q->param[0]=Ymesh=STRTOI(GetArgInc(&i)); // nb de colonnes
	Q->param[1]=Xmesh=STRTOI(GetArgInc(&i)); // nb de lignes
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("alkane"){ NextArg(&i);
	Q->adj=alkane; LABEL=1;
	Q->param[0]=-1;
	if(EQUAL("normal")||EQUAL("no")) Q->param[0]=ALK_NOR;
	if(EQUAL("cyclo") ||EQUAL("cy")) Q->param[0]=ALK_CYC;
	if(EQUAL("iso")   ||EQUAL("is")) Q->param[0]=ALK_ISO;
	if(EQUAL("neo")   ||EQUAL("ne")) Q->param[0]=ALK_NEO;
	if(EQUAL("sec")   ||EQUAL("se")) Q->param[0]=ALK_SEC;
	if(EQUAL("tert")  ||EQUAL("te")) Q->param[0]=ALK_TER;
	if(Q->param[0]<0) Erreur(6);
	i++;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if SUFFIX("ane"){
	Q->param[0]=ALK_NOR; Q->param[1]=-1;
	if PREFIX("cyclo-")  Q->param[0]=ALK_CYC;
	if PREFIX("iso-")    Q->param[0]=ALK_ISO;
	if PREFIX("neo-")    Q->param[0]=ALK_NEO;
	if PREFIX("sec-")    Q->param[0]=ALK_SEC;
	if PREFIX("tert-")   Q->param[0]=ALK_TER;
	if SUFFIX("ethane")  Q->param[1]=2; // "ethane" prefix de "methane"
	if SUFFIX("methane") Q->param[1]=1; // attention Ã  l'ordre
	if SUFFIX("propane") Q->param[1]=3;
	if SUFFIX("butane")  Q->param[1]=4;
	if SUFFIX("pentane") Q->param[1]=5;
	if SUFFIX("hexane")  Q->param[1]=6;
	if SUFFIX("heptane") Q->param[1]=7;
	if SUFFIX("octane")  Q->param[1]=8;
	if SUFFIX("nonane")  Q->param[1]=9;
	if SUFFIX("alkane"){ i++; Q->param[1]=STRTOI(GetArgInc(&i)); }
	if(Q->param[1]>=0){ Q->adj=alkane; LABEL=1; goto param0; }
      }

    /* graphes de base avec type de paramÃ¨tres dÃ©jÃ  rencontrÃ©s */

    if EQUAL("clebsch"){ Q->adj=clebsch; goto param_hypercube; }
    if EQUAL("loadc"){ LOADC=1; FORMAT=F_no; goto param_load; }
    //
    if EQUAL("ring"){ Q->adj=ring; goto param_grid; }
    if EQUAL("cage"){ Q->adj=cage; goto param_grid; }
    if EQUAL("bdrg"){ Q->adj=bdrg; goto param_grid; }
    if EQUAL("fdrg"){ Q->adj=fdrg; goto param_grid; }
    if EQUAL("ggosset"){ Q->adj=ggosset; goto param_grid; }
    if EQUAL("rectree"){ Q->adj=rectree; goto param_grid; }
    if EQUAL("rpartite"){ Q->adj=rpartite; goto param_grid; }
    if EQUAL("aqua"){ Q->adj=aqua; Q->directed=1; goto param_grid; }
    if EQUAL("collatz"){ Q->adj=collatz; Q->directed=Q->loop=1; goto param_grid; }
    //
    if EQUAL("icosahedron"){ Q->adj=icosahedron; goto param0; }
    if EQUAL("rdodecahedron"){ Q->adj=rdodecahedron; goto param0; }
    if EQUAL("herschel"){ Q->adj=herschel; goto param0; }
    if EQUAL("goldner-harary"){ Q->adj=goldner_harary; goto param0; }
    if EQUAL("triplex"){ Q->adj=triplex; goto param0; }
    if EQUAL("jaws"){ Q->adj=jaws; goto param0; }
    if EQUAL("starfish"){ Q->adj=starfish; goto param0; }
    if EQUAL("fritsch"){ Q->adj=fritsch; goto param0; }
    if EQUAL("zamfirescu"){ Q->adj=zamfirescu; goto param0; }
    if EQUAL("hatzel"){ Q->adj=hatzel; goto param0; }
    if EQUAL("soifer"){ Q->adj=soifer; goto param0; }
    if EQUAL("poussin"){ Q->adj=poussin; goto param0; }
    if EQUAL("errera"){ Q->adj=errera; goto param0; }
    if EQUAL("kittell"){ Q->adj=kittell; goto param0; }
    if EQUAL("frucht"){ Q->adj=frucht; goto param0; }
    if EQUAL("moser"){ Q->adj=moser; goto param0; }
    if EQUAL("markstrom"){ Q->adj=markstrom; goto param0; }
    if EQUAL("robertson"){ Q->adj=robertson; goto param0; }
    if EQUAL("heawood4"){ Q->adj=heawood4; goto param0; }
    if EQUAL("wiener-araya"){ Q->adj=wiener_araya; goto param0; }
    if EQUAL("hgraph"){ Q->adj=hgraph; goto param0; }
    if(EQUAL("rgraph")||EQUAL("fish")){ Q->adj=rgraph; goto param0; }
    if EQUAL("cricket"){ Q->adj=cricket; goto param0; }
    if EQUAL("moth"){ Q->adj=moth; goto param0; }
    if EQUAL("dart"){ Q->adj=dart; goto param0; }
    if EQUAL("antenna"){ Q->adj=antenna; goto param0; }
    if EQUAL("suzuki"){ Q->adj=suzuki; goto param0; }
    if EQUAL("bull"){ Q->adj=bull; goto param0; }
    if EQUAL("harborth"){ Q->adj=harborth; goto param0; }
    if EQUAL("doily"){ Q->adj=doily; goto param0; }
    if EQUAL("schlafli"){ Q->adj=schlafli; goto param0; }
    //
    if EQUAL("gear"){ Q->adj=gear; goto param1; }
    if EQUAL("pstar"){ Q->adj=pstar; Q->param[1]=2; goto param1; }
    if EQUAL("paley"){ Q->adj=paley; goto param1; }
    if(EQUAL("comb")||EQUAL("centipede")){ Q->adj=comb; goto param1; }
    if EQUAL("sunlet"){ Q->adj=sunlet; goto param1; }
    if EQUAL("mycielski"){ Q->adj=mycielski; goto param1; }
    if EQUAL("treep"){ Q->adj=treep; goto param1; }
    if EQUAL("halin"){ Q->adj=halin; goto param1; }
    if EQUAL("windmill"){ Q->adj=windmill; goto param1; }
    if EQUAL("interval"){ Q->adj=interval; goto param1; }
    if EQUAL("circle"){ Q->adj=circle; goto param1; }
    if EQUAL("permutation"){ Q->adj=permutation; goto param1; }
    if EQUAL("pancake"){ Q->adj=pancake; goto param1; }
    if EQUAL("bpancake"){ Q->adj=bpancake; goto param1; }
    if EQUAL("crown"){ Q->adj=crown; goto param1; }
    if EQUAL("shuffle"){ Q->adj=shuffle; goto param1; }
    if EQUAL("flip"){ Q->adj=flip; goto param1; }
    if EQUAL("apollonian"){ Q->adj=apollonian; goto param1; }
    if EQUAL("flower_snark"){ Q->adj=flower_snark; goto param1; }
    if EQUAL("gabriel"){ Q->adj=gabriel; POS=1; goto param1; }
    if EQUAL("sgabriel"){ Q->adj=sgabriel; POS=1; /*FORMAT=F_userdot;*/ goto param1; }
    if EQUAL("rng"){ Q->adj=rng; POS=1; goto param1; }
    if EQUAL("mst"){ Q->adj=mst; POS=1; goto param1; }
    if EQUAL("antiprism"){ Q->adj=antiprism; Q->param[1]=1; goto param1; }
    if EQUAL("butterfly"){ Q->adj=butterfly; goto param1; }
    if EQUAL("matching"){ Q->adj=matching; goto param1; }
    if EQUAL("polygon"){ Q->adj=polygon; goto param1; }
    if EQUAL("cactus"){ Q->adj=cactus; goto param1; }
    if EQUAL("helm"){ Q->adj=helm; goto param1; }
    if EQUAL("haar"){ Q->adj=haar; goto param1; }
    if EQUAL("margulis"){ Q->adj=margulis; goto param1; }
    if EQUAL("parachute"){ Q->adj=parachute; goto param1; }
    //
    if EQUAL("kout"){ Q->adj=kout; goto param2; }
    if EQUAL("ktree"){ Q->adj=ktree; Q->param[2]=0; goto param2; }
    if EQUAL("gpetersen"){ Q->adj=gpetersen; goto param2; }
    if EQUAL("debruijn"){ Q->adj=debruijn; goto param2; }
    if EQUAL("kautz"){ Q->adj=kautz; goto param2; }
    if EQUAL("gpstar"){ Q->adj=gpstar; goto param2; }
    if EQUAL("hexagon"){ Q->adj=hexagon; goto param2; }
    if EQUAL("whexagon"){ Q->adj=whexagon; goto param2; }
    if EQUAL("hanoi"){ Q->adj=hanoi; goto param2; }
    if EQUAL("sierpinski"){ Q->adj=sierpinski; goto param2; }
    if EQUAL("banana"){ Q->adj=banana; goto param2; }
    if EQUAL("kpage"){ Q->adj=kpage; goto param2; }
    if EQUAL("line-graph"){ Q->adj=linegraph; goto param2; }
    if EQUAL("linial"){ Q->adj=linial; goto param2; }
    if EQUAL("linialc"){ Q->adj=linialc; goto param2; }
    if EQUAL("expander"){ Q->adj=expander; goto param2; }
    if EQUAL("fan"){ Q->adj=fan; goto param2; }
    if EQUAL("split"){ Q->adj=split; goto param2; }
    if EQUAL("behrend"){ Q->adj=behrend; goto param2; }
    if EQUAL("turan"){ Q->adj=turan; goto param2; }
    if EQUAL("klein"){ Q->adj=klein; POS=1; goto param2; }
    if EQUAL("knng"){ Q->adj=knng; POS=1; goto param2; }
    //
    if EQUAL("rarytree"){ Q->adj=rarytree; goto param3; }
    if EQUAL("barbell"){ Q->adj=barbell; goto param3; }
    if EQUAL("planar"){ Q->adj=planar; goto param3; }
    if EQUAL("hyperbolic"){ Q->adj=hyperbolic; goto param3; }
    if EQUAL("kneser"){ Q->adj=kneser; goto param3; }
    if EQUAL("pat"){ Q->adj=pat; POS=1; goto param3; }
    if EQUAL("uno"){ Q->adj=uno; POS=1; goto param3; }
    if EQUAL("ngon"){ Q->adj=ngon; POS=1; goto param3; }
    if EQUAL("wpsl") { Q->adj=wpsl; POS=1; Q->param[3]=0; goto param3; }
    if EQUAL("wpsld"){ Q->adj=wpsl; POS=1; Q->param[3]=1; goto param3; }
    if EQUAL("upsl") { Q->adj=wpsl; POS=1; Q->param[3]=2; goto param3; }
    if EQUAL("upsld"){ Q->adj=wpsl; POS=1; Q->param[3]=3; goto param3; }
    if EQUAL("wdis") { Q->adj=wpsl; POS=1; Q->param[3]=4; goto param3; }
    if EQUAL("wdisd"){ Q->adj=wpsl; POS=1; Q->param[3]=5; goto param3; }
    if EQUAL("udis") { Q->adj=wpsl; POS=1; Q->param[3]=6; goto param3; }
    if EQUAL("udisd"){ Q->adj=wpsl; POS=1; Q->param[3]=7; goto param3; }
    //
    if EQUAL("chess"){ Q->adj=chess; goto param4; }

    /********************/
    /* graphes composÃ©s */
    /********************/

    if EQUAL("theta"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=3;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[2]=max(1,Q->param[2]);
	Q->dparam[0]=6.0/Q->param[2];
	goto fin;
      }
    if EQUAL("dtheta"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=3;
	Q->param[2]=STRTOI(GetArgInc(&i))/2;
	Q->param[2]=max(1,Q->param[2]);
	Q->dparam[0]=3.0/Q->param[2];
	goto fin;
      }
    if EQUAL("yao"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=0;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[2]=max(1,Q->param[2]);
	Q->dparam[0]=2.0/Q->param[2];
	goto fin;
      }
    if EQUAL("path"){ i++;
	Q->adj=grid;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("torus"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=-STRTOI(GetArgInc(&i));
	Q->param[2]=-STRTOI(GetArgInc(&i));
	goto fin;
    }
    if EQUAL("mesh"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("mobius"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	Q->param[0]=3;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1;
	Q->param[3]=Q->param[1]/2;
	goto fin;
      }
    if EQUAL("ladder"){ i++;
	Q->adj=grid;
	Q->param[0]=Q->param[1]=2;
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("johnson"){ i++;
	Q->adj=kneser;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=Q->param[1]-2;
	Q->not=!Q->not;
	goto fin;
      }
    if EQUAL("star"){ i++;
	Q->adj=rpartite;
	Q->param[0]=2;
	Q->param[1]=1;
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("bipartite"){ i++;
	Q->adj=rpartite;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("clique"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	Q->not=!Q->not;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("cycle"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1;
	goto fin;
      }
    if EQUAL("star-polygon"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	POS=1;
	XYtype=XY_DISK;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1;
	goto fin;
      }
    if EQUAL("convex-polygon"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	POS=1;
	XYtype=XY_CONVEX;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1;
	goto fin;
      }
    if(EQUAL("stable")||EQUAL("empty")){ i++;
      Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
      Q->param[0]=1;
      Q->param[1]=STRTOI(GetArgInc(&i));
      goto fin;
    } 
    if EQUAL("point"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	POS=1;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      } 
    if EQUAL("random"){ i++;
	Q->adj=ring; // attention ! paramÃ¨tres de type "grid"
	Q->not=!Q->not;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i));
	DELE=1-STRTOD(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("cylinder"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=-STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("caterpillar"){ i++;
	//Q->adj=grid;
	//k=STRTOI(GetArgInc(&i)); /* nb de sommets total */
	//STAR=randomu(k); /* entre 0...k-1 sommets de deg=1. Active l'opÃ©ration star() */
	//Q->param[0]=1;
	//Q->param[1]=k-STAR; /* Q->param[0]=nb de sommets du chemin, qui est >=1 */
	Q->adj=NULL; // graphe dÃ©sactivÃ©
	goto fin;
      }
    if EQUAL("sunflower"){ i++;
	Q->adj=cage;
	Q->param[0]=3;
	Q->param[1]=2*STRTOI(GetArgInc(&i));
	Q->param[2]=2;
	Q->param[3]=0;
	goto fin;
      }
    if EQUAL("wheel"){ i++;
	Q->adj=ringarytree;
	Q->param[0]=1;
	Q->param[1]=0;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=2;
	goto fin;
      }
    if(EQUAL("tadpole")||EQUAL("dragon")){ i++;
      Q->adj=barbell;
      Q->param[0]=-STRTOI(GetArgInc(&i));
      Q->param[1]=1;
      Q->param[2]=STRTOI(GetArgInc(&i));
      goto fin;
    }
    if EQUAL("pan"){ i++;
	Q->adj=barbell;
	Q->param[0]=-STRTOI(GetArgInc(&i));
	Q->param[1]=1;
	Q->param[2]=1;
	goto fin;
      }
    if EQUAL("web"){ i++;
	Q->adj=ringarytree;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=1;
	Q->param[3]=2;
	goto fin;
      }
    if EQUAL("percolation"){ i++;
	Q->adj=udg; POS=1; XYtype=XY_MESH; ROUND=0;
	Xmesh=STRTOI(GetArgInc(&i));
	Ymesh=STRTOI(GetArgInc(&i));
	Q->param[0]=Xmesh*Ymesh; /* nombre de sommets */
	Q->dparam[0]=1; /* rayon */
	DELE=1-CheckProba(STRTOD(GetArgInc(&i))); /* proba existence arÃªte */
	NORM=NORM_L1;
	goto fin;
      }
    if EQUAL("hudg"){ i++; /* paramÃ©trage Ã  revoir */
	Q->adj=udg; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	XYtype=XY_HYPER;
	XYpower=Q->dparam[0];
	NORM=NORM_HYPER;
	goto fin;
      }
    if EQUAL("plrg"){ i++;
	Q->adj=bdrg;
	int *S=power_law_seq(STRTOI(GetArgInc(&i)),STRTOD(GetArgInc(&i)),NULL);
	if(S==NULL) Erreur(6);
	// attention ! S pourrait Ãªtre de taille < PARMAX
	if(S[0]+1<PARMAX) REALLOC(S,PARMAX);
	free(Q->param); Q->param=S;
	goto fin;
      }
    if EQUAL("cubic"){ i++;
	Q->adj=fdrg;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=3;
	goto fin;
      }
    if EQUAL("regular"){ i++;
	Q->adj=fdrg;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("syracuse"){ i++;
	Q->adj=collatz;
	Q->directed=Q->loop=1;
	// attention ! paramÃ¨tres de type "grid"
	Q->param[0]=5;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1; Q->param[3]=0;
	Q->param[4]=3; Q->param[5]=1;
	goto fin;
      }
    if EQUAL("kakutami_3x+1"){ i++;
	Q->adj=collatz;
	Q->directed=Q->loop=1;
	// attention ! paramÃ¨tres de type "grid"
	Q->param[0]=5;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=1; Q->param[3]=0;
	Q->param[4]=6; Q->param[5]=2;
	goto fin;
      }
    if EQUAL("kakutami_5x+1"){ i++;
	Q->adj=collatz;
	Q->directed=Q->loop=1;
	// attention ! paramÃ¨tres de type "grid"
	Q->param[0]=13;
	Q->param[1]=STRTOI(GetArgInc(&i));
	//                                 mod6  relation
	Q->param[ 2]= 3; Q->param[3]= 0; //  0    div2
	Q->param[ 4]=30; Q->param[5]= 6; //  1    5x+1
	Q->param[ 6]= 3; Q->param[7]= 0; //  2    div2
	Q->param[ 8]= 2; Q->param[9]= 0; //  3    div3
	Q->param[10]= 3; Q->param[11]=0; //  4    div2
	Q->param[12]=30; Q->param[13]=6; //  5    5x+1
	goto fin;
      }
    if EQUAL("kakutami_7x+1"){ i++;
	Q->adj=collatz;
	Q->directed=Q->loop=1;
	// attention ! paramÃ¨tres de type "grid"
	Q->param[0]=61;
	Q->param[1]=STRTOI(GetArgInc(&i));	
	//                                   mod30  relation
	Q->param[ 2]= 15; Q->param[3]=  0; //  0     div2
	Q->param[ 4]=210; Q->param[5]= 30; //  1     7x+1
	Q->param[ 6]= 15; Q->param[7]=  0; //  2     div2
	Q->param[ 8]= 10; Q->param[9]=  0; //  3     div3
	Q->param[10]= 15; Q->param[11]= 0; //  4     div2
	Q->param[12]=  6; Q->param[13]= 0; //  5     div5
	Q->param[14]= 15; Q->param[15]= 0; //  6     div2
	Q->param[16]=210; Q->param[17]=30; //  7     7x+1
	Q->param[18]= 15; Q->param[19]= 0; //  8     div2
	Q->param[20]= 10; Q->param[21]= 0; //  9     div3
	Q->param[22]= 15; Q->param[23]= 0; // 10     div2
	Q->param[24]=210; Q->param[25]=30; // 11     7x+1
	Q->param[26]= 15; Q->param[27]= 0; // 12     div2
	Q->param[28]=210; Q->param[29]=30; // 13     7x+1
	Q->param[30]= 15; Q->param[31]= 0; // 14     div2
	Q->param[32]= 10; Q->param[33]= 0; // 15     div3
	Q->param[34]= 15; Q->param[35]= 0; // 16     div2
	Q->param[36]=210; Q->param[37]=30; // 17     7x+1
	Q->param[38]= 15; Q->param[39]= 0; // 18     div2
	Q->param[40]=210; Q->param[41]=30; // 19     7x+1
	Q->param[42]= 15; Q->param[43]= 0; // 20     div2
	Q->param[44]= 10; Q->param[45]= 0; // 21     div3
	Q->param[46]= 15; Q->param[47]= 0; // 22     div2
	Q->param[48]=210; Q->param[49]=30; // 23     7x+1
	Q->param[50]= 15; Q->param[51]= 0; // 24     div2
	Q->param[52]=  6; Q->param[53]= 0; // 25     div5
	Q->param[54]= 15; Q->param[55]= 0; // 26     div2
	Q->param[56]= 10; Q->param[57]= 0; // 27     div3
	Q->param[58]= 15; Q->param[59]= 0; // 28     div2
	Q->param[60]=210; Q->param[61]=30; // 29     7x+1
	goto fin;
      }
    if EQUAL("farkas"){ i++;
	Q->adj=collatz;
	Q->directed=Q->loop=1;
	// attention ! paramÃ¨tres de type "grid"
	Q->param[0]=25;
	Q->param[1]=STRTOI(GetArgInc(&i));
	//                                 mod12 relation
	Q->param[ 2]= 6; Q->param[3]= 0; //  0    div2
	Q->param[ 4]= 6; Q->param[5]= 6; //  1    (x+1)/2
	Q->param[ 6]= 6; Q->param[7]= 0; //  2    div2
	Q->param[ 8]= 4; Q->param[9]= 0; //  3    div3
	Q->param[10]= 6; Q->param[11]=0; //  4    div2
	Q->param[12]= 6; Q->param[13]=6; //  5    (x+1)/2
	Q->param[14]= 6; Q->param[15]=0; //  6    div2
	Q->param[16]=18; Q->param[17]=6; //  7    (3x+1)/2
	Q->param[18]= 6; Q->param[19]=0; //  8    div2
	Q->param[20]= 4; Q->param[21]=0; //  9    div3
	Q->param[22]= 6; Q->param[23]=0; // 10    div2    
	Q->param[24]=18; Q->param[25]=6; // 11    (3x+1)/2
	goto fin;
      }

    /* graphes composÃ©s avec type de paramÃ¨tres dÃ©jÃ  rencontrÃ©s */

    if EQUAL("tree"){ Q->adj=arboricity; Q->param[1]=1; goto param1; }
    if EQUAL("rbinary"){ Q->adj=rarytree; Q->param[1]=2; Q->param[2]=0; goto param1; }
    if EQUAL("rbinaryz"){ Q->adj=rarytree; Q->param[1]=2; Q->param[2]=1; goto param1; }
    if EQUAL("outerplanar"){ Q->adj=kpage; Q->param[1]=1; goto param1; }
    if EQUAL("squaregraph"){ Q->adj=planar; Q->param[1]=Q->param[2]=4; goto param1; }
    if EQUAL("prism"){ Q->adj=gpetersen; Q->param[1]=1; goto param1; }
    if EQUAL("d-octahedron"){ Q->adj=matching; Q->not=!Q->not; goto param1; }
    if EQUAL("nng"){ Q->adj=knng; Q->param[1]=1; POS=1; goto param1; }
    if EQUAL("binary"){ Q->adj=ringarytree;
	Q->param[1]=Q->param[2]=2;Q->param[3]=0; goto param1; }
    if EQUAL("td-delaunay"){ Q->adj=thetagone;
	POS=1; Q->param[1]=Q->param[2]=3; Q->dparam[0]=1; goto param1; }
    //
    if EQUAL("kpath"){ Q->adj=ktree; Q->param[2]=1; goto param2; }
    if EQUAL("kstar"){ Q->adj=ktree; Q->param[2]=2; goto param2; }
    if EQUAL("tw"){ Q->adj=ktree; Q->param[2]=0; DELE=.5; goto param2; }
    if EQUAL("pw"){ Q->adj=ktree; Q->param[2]=1; DELE=.5; goto param2; }
    if EQUAL("knight"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=2; goto param2; }
    if EQUAL("camel"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=3; goto param2; }
    if EQUAL("giraffe"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=4; goto param2; }
    if EQUAL("zebra"){ Q->adj=chess; Q->param[2]=2;Q->param[3]=3; goto param2; }
    if EQUAL("antelope"){ Q->adj=chess; Q->param[2]=2;Q->param[3]=4; goto param2; }
    if EQUAL("lollipop"){ Q->adj=barbell; Q->param[2]=0; goto param2; }
    //
    if EQUAL("arytree"){ Q->adj=ringarytree; Q->param[3]=0; goto param3; }

    /* graphes sans paramÃ¨tres mais composÃ©s d'un graphe de base
       contenant des paramÃ¨tres, on doit donc passer par CheckHelp()
       avec goto param0; */

    // 1 paramÃ¨tre
    if(EQUAL("cube")||EQUAL("hexahedron")){ Q->adj=crown; Q->param[0]=4; goto param0; }
    if EQUAL("octahedron"){ Q->adj=antiprism; Q->param[0]=3; goto param0; }
    if EQUAL("dodecahedron"){ Q->adj=gpetersen; Q->param[0]=10; Q->param[1]=2; goto param0; }
    if EQUAL("associahedron"){ Q->adj=flip; Q->param[0]=6; goto param0; }
    if EQUAL("tietze"){ Q->adj=flower_snark; Q->param[0]=3; goto param0; }
    if EQUAL("grotzsch"){ Q->adj=mycielski; Q->param[0]=4; goto param0; }
    if EQUAL("egraph"){ Q->adj=comb; Q->param[0]=3; goto param0; }
    if EQUAL("parapluie"){ Q->adj=parachute; Q->not=!Q->not; goto param0; }
    // 2 paramÃ¨tres
    if EQUAL("nauru"){ Q->adj=pstar; Q->param[0]=4;Q->param[1]=2; goto param0; }
    if EQUAL("hajos"){ Q->adj=sierpinski; Q->param[0]=2;Q->param[1]=3; goto param0; }
    if EQUAL("netgraph"){ Q->adj=sierpinski;Q->not=!Q->not;Q->param[0]=2;Q->param[1]=3; goto param0; }
    if EQUAL("house"){ Q->adj=grid; Q->not=!Q->not; Q->param[0]=1;Q->param[1]=5; goto param0; }
    if EQUAL("tetrahedron"){ Q->adj=ring; Q->not=!Q->not; Q->param[0]=4; Q->param[1]=0; goto param0; }
    if EQUAL("claw"){ Q->adj=rpartite; Q->param[0]=2;Q->param[1]=1;Q->param[2]=3; goto param0; }
    if EQUAL("desargues"){ Q->adj=gpetersen; Q->param[0]=10;Q->param[1]=3; goto param0; }
    if EQUAL("durer"){ Q->adj=gpetersen; Q->param[0]=6;Q->param[1]=2; goto param0; }
    if EQUAL("mobius-kantor"){ Q->adj=gpetersen; Q->param[0]=8;Q->param[1]=3; goto param0; }
    if EQUAL("gem"){ Q->adj=fan; Q->param[0]=4;Q->param[1]=1; goto param0; }
    if EQUAL("diamond"){ Q->adj=fan; Q->param[0]=Q->param[1]=2; goto param0; }
    if EQUAL("cross"){ Q->adj=banana; Q->param[0]=1; Q->param[1]=4; goto param0; }
    if(EQUAL("tgraph")||EQUAL("fork")){ Q->adj=banana;Q->param[0]=1;Q->param[1]=3; goto param0; }
    if EQUAL("ygraph"){ Q->adj=banana; Q->param[0]=3;Q->param[1]=1; goto param0; }
    if EQUAL("cuboctahedron"){ Q->adj=linial; Q->param[0]=4;Q->param[1]=2; goto param0; }
    // 3 paramÃ¨tres
    if EQUAL("petersen"){ Q->adj=kneser; Q->param[0]=5;Q->param[1]=2;Q->param[2]=0; goto param0; }
    if EQUAL("banner"){ Q->adj=barbell; Q->param[0]=-4;Q->param[1]=1;Q->param[2]=1; goto param0; }
    if EQUAL("paw"){ Q->adj=barbell; Q->param[0]=-3;Q->param[1]=1;Q->param[2]=1; goto param0; }
    if EQUAL("theta0"){ Q->adj=barbell; Q->param[0]=Q->param[1]=-5;Q->param[2]=-2; goto param0; }
    if EQUAL("utility"){ Q->adj=rpartite; Q->param[0]=2;Q->param[1]=Q->param[2]=3; goto param0; }
    if EQUAL("domino"){ Q->adj=grid; Q->param[0]=Q->param[1]=2;Q->param[2]=3; goto param0; }
    if EQUAL("kite"){ Q->adj=banana; Q->not=!Q->not;Q->param[0]=1;Q->param[1]=3; goto param0; }
    if EQUAL("hourglass"){ Q->adj=barbell; Q->param[0]=Q->param[1]=3;Q->param[2]=0; goto param0; }
    // 4 paramÃ¨tres
    if EQUAL("wagner"){ Q->adj=ring;
	Q->param[0]=3;
	Q->param[1]=8;
	Q->param[2]=1; Q->param[3]=4;
	goto param0;
      }
    if EQUAL("heawood"){ Q->adj=cage;
	Q->param[0]=3; 
	Q->param[1]=14;
	Q->param[2]=5; Q->param[3]=-5;
	goto param0;
      }
    if EQUAL("franklin"){ Q->adj=cage;
	Q->param[0]=3; 
	Q->param[1]=12;
	Q->param[2]=5; Q->param[3]=-5;
	goto param0;
      }
    // 5 paramÃ¨tres
    if EQUAL("mcgee"){ Q->adj=cage;
	Q->param[0]=4; 
	Q->param[1]=24;
	Q->param[2]=12; Q->param[3]=7; Q->param[4]=-7;
	goto param0;
      }
    if EQUAL("bidiakis"){ Q->adj=cage;
	Q->param[0]=4;
	Q->param[1]=12;
	Q->param[2]=-4; Q->param[3]=6; Q->param[4]=4;
	goto param0;
      }
    // 6 paramÃ¨tres
    if EQUAL("dyck"){ Q->adj=cage;
	Q->param[0]=5; 
	Q->param[1]=32;
	Q->param[2]=5; Q->param[3]=0; Q->param[4]=13; Q->param[5]=-13;
	goto param0;
      }
    if EQUAL("gosset"){ Q->adj=ggosset;
	Q->param[0]=5;
	Q->param[1]=8;
	Q->param[2]=2; Q->param[3]=3; Q->param[4]=6; Q->param[5]=-1;
	goto param0;
      }
    // 8 paramÃ¨tres
    if EQUAL("pappus"){ Q->adj=cage;
	Q->param[0]=7; 
	Q->param[1]=18;
        Q->param[2]=5; Q->param[3]=7; Q->param[4]=-7;
	Q->param[5]=7; Q->param[6]=-7; Q->param[7]=-5;
	goto param0;
      }
    if EQUAL("tutte-coexter"){ Q->adj=cage;
	Q->param[0]=7; 
	Q->param[1]=30;
	Q->param[2]=-7; Q->param[3]=9; Q->param[4]=13;
	Q->param[5]=-13; Q->param[6]=-9; Q->param[7]=7;
	goto param0;
      }
    if EQUAL("gray"){ Q->adj=cage;
	Q->param[0]=7;
	Q->param[1]=54; 
	Q->param[2]=7; Q->param[3]=-7;
	Q->param[4]=25; Q->param[5]=-25; Q->param[6]=13; Q->param[7]=-13;
	goto param0;
      }
    // 14 paramÃ¨tres
    if EQUAL("chvatal"){ Q->adj=cage;
	Q->param[0]=13;
	Q->param[1]=12;
	Q->param[2]=Q->param[4]=Q->param[7]=Q->param[10]=Q->param[12]=Q->param[13]=3;
	Q->param[3]=Q->param[5]=Q->param[6]=Q->param[8]=6;
	Q->param[9]=Q->param[11]=-3;
	goto param0;
      }

  fin:
    if(j==i){
      if PREFIX("-") Erreur(2); /* option non trouvÃ©e */
      Erreur(10); /* graphe non trouvÃ© */
    }
    
  } /* fin du while(i<ARGC) ... */

  if((Q->adj==NULL)&&(FAMILY==NULL)) Erreur(10); // graphe non trouvÃ©
  if(FAMILY) goto check; // cas d'une option -filter, on saute Ã  la fin

  /* options qui ne vont pas ensemble */

  if((LOADC)&&(PERMUTE||Q->not)) Erreur(29); // options incompatibles
  if((LOADC)&&(CHECK<=CHECK_ON)) Erreur(30); // manque -check


  /***********************************

           COEUR DU GENERATEUR

  ***********************************/

  /* initialisation du graphe, calcule Q->n avant la suppression
     Ã©ventuelle de sommets, lit le graphe Q->G si adj=load, dÃ©termine
     (Q->xpos,Q->ypos) si graphe gÃ©omÃ©trique */  

  Q->code=QUERY_INIT; Q->seed=SEED; // seed au moment de l'initialisation du graphe
  if(Q->adj(Q)&&(Q->error)) Erreur(Q->error); // exÃ©cute l'initialisation, vÃ©rifie code erreur
  if(Q->n<0) Q->n=0; // ne devrait jamais arriver
 
  if(LOADC){ GF=Q->G; NF=Q->n; Q->adj=NULL; goto check; } // saute la partie gÃ©nÃ©ration d'arÃªtes
  if(POS && Q->xpos==NULL) InitXY(Q); // il faut dÃ©terminer les positions (si pas dÃ©jÃ  fait)
  if(abs(LABEL)==1) PERMUTE=0; // si on souhaite les labels d'origine, on ne permute rien

  ALLOC(V,Q->n);       // V[i]=Ã©tiquette du sommet i, -1 si i est supprimÃ©
  ALLOC(VF,Q->n);      // VF[j]=indice du j-Ã¨me sommet non supprimÃ©
  NF=InitVertex(Q->n,DELV); // initialise V[i], VF[i] et renvoie NF=#sommets final
  ALLOC(INC,Q->n);     // INC[i]=1 ssi i possÃ¨de un voisin, 0 si sommet isolÃ©
  for(j=0;j<NF;j++) INC[VF[j]]=0; // au dÃ©part que des sommets isolÃ©s (=deg 0)
  
  /* constantes pour accÃ©lÃ©rer les tests de la boucle principale */
  long const seuil_edge=(1.0-DELE)*RAND_MAX; // dans [0,RAND_MAX]
  long const seuil_redirect=(double)(REDIRECT)*RAND_MAX; // dans [0,RAND_MAX]

  /*
    GÃ©nÃ¨re les adjacences i-j en tenant compte des sommets isolÃ©s et
    des sommets supprimÃ©s. Les sommets isolÃ©s sont affichÃ©s en dernier
    Ã  cause de l'option -redirect. On a toujours i<j lorsque l'arÃªte
    i-j doit Ãªtre sortie et que DIRECTED=0. Si on a FAST=1, alors on
    gÃ©nÃ¨re le graphe Ã  partir de Q->G, s'il existe.
  */

  /* initialise le format d'affichage */
  Q->code=QUERY_INIT; Out(Q);
  
  /* si G existe, affiche que les arÃªtes de G en O(m+n) */
  if(FAST){
    if(Q->G){
      int const redirect=(REDIRECT!=0);
      int d,t;

      /* on ne teste que les arÃªtes de Q->G */
      for(i=0;i<Q->n;i++)
	if(V[i]>=0) // si le sommet i existe
	  for(t=0,d=Q->G->d[i];t<d;t++){
	    j=Q->G->L[i][t]; if((Q->G->sym)&&(j<i)) continue; // il faut i<j si Q->G symÃ©trique
	    if((V[j]>=0)&&((!Q->directed)||(i!=j)||(Q->loop))) // si j existe ...
	      if(random()<=seuil_edge){
		if(redirect){ // si redirection
		  j=(random()<=seuil_redirect)? randomu(Q->n) : j;
		  if((V[j]<0)||(j==i)) continue; // prochain voisin
		}
		INC[i]++; // un voisin de plus pour i
		INC[j]++; // un voisin de plus pour j
		Q->code=QUERY_ADJ,Q->i=i,Q->j=j;
		Out(Q); // sort l'arÃªte i-j avec i<j
	      }
	  }
    }else Erreur(25); // -fast alors que Q->G n'existe pas
  }else{
    int const noredirect=(REDIRECT==0);

    /* teste et affiche les O(n^2) arcs ou arÃªtes possibles */

    for(i=0;i<Q->n;i++)    // pour tous les
      if(V[i]>=0){         // si i existe
	for(j=(Q->directed)?0:i+(Q->loop==0);j<Q->n;j++) // pour tous les j>i
	  if((V[j]>=0)&&((!Q->directed)||(i!=j)||(Q->loop))) // si j existe ...
	    if(random()<=seuil_edge){
	      Q->code=QUERY_ADJ,Q->i=i,Q->j=j; Q->adj(Q);
	      Q->a |= (Q->loop==2)&&(i==j); // si boucle forcÃ©e
	      if((Q->a)^(Q->not)){
		/* ici l'arÃªte i-j devrait Ãªtre sortie */
		if(noredirect){ // si pas de redirection d'arÃªte
		  INC[i]++; // un voisin de plus pour i
		  INC[j]++; // un voisin de plus pour j
		  Q->i=i,Q->j=j; // ici Q->code=QUERY_ADJ
		  Out(Q); // sort l'arÃªte i-j
		}else{ // on redirige l'arÃªte i-j vers i-k
		  k=(random()<=seuil_redirect)? randomu(Q->n) : j;
		  if((V[k]>=0)&&(k!=i)){
		    /* on affiche l'arÃªte que si k existe et si k<>i.
		       Attention ! il ne faut toucher ni Ã  i ni Ã  j */
		    INC[i]++; // un voisin de plus pour i
		    INC[k]++; // un voisin de plus pour k
		    if(k<i){
		      Q->i=k,Q->j=i;
		      Out(Q);
		    }else{
		      Q->i=i,Q->j=k;
		      Out(Q); // pour avoir i<j
		    }
		  }
		}
	      }
	    }
      }

  }

  /* affiche les sommets isolÃ©s */
  for(Q->i=0;Q->i<Q->n;Q->i++)
    if((V[Q->i]>=0)&&(!INC[Q->i])){ Q->code=QUERY_ISOL; Out(Q); }
  
  /* fin de l'affichage, doit Ãªtre fait avant adj() avec QUERY_END */
  Q->code=QUERY_END; Out(Q); /* NB: calcule GF si CHECK */

  if((CHECK)&&(POS)){ // mÃ©morise Q->xpos,Q->ypos qui vont Ãªtre supprimÃ©s par adj() avec QUERY_END
    ALLOCZ(GF->xpos,NF,Q->xpos[VF[_i]]);
    ALLOCZ(GF->ypos,NF,Q->ypos[VF[_i]]);
  }

  free(V);
  free(VF);
  free(INC);
  
 /* termine la fonction Q->adj() */
  Q->code=QUERY_END;
  Q->adj(Q);
  Q->adj=NULL;

  /***********************************

           FIN DU GENERATEUR

         Si on a CHECK<>0 alors:
         - le graphe gÃ©nÃ©rÃ© est GF
         - son nombre de sommets est NF

  ***********************************/

 check:
  /* NB: dans le cas LOADC on ne fait pas adj(Q) avec QUERY_END,
     c'est-Ã -dire load(Q) car sinon load(Q) va libÃ¨rer le graphe
     Q->G. Or on veut conserver GF=Q->G prÃ©cisÃ©ment. */

  if(CHECK){
    if((GF==NULL)||(GF->n!=NF)) Erreur(32); /* ne devrait jamais arriver */

    switch(CHECK){

    case CHECK_MAINCC:{
      param_dfs *p=dfs(GF,MEM(CPARAM,0,int),NULL);
      int d0,d1,n,c;
      for(i=c=n=d0=d1=0;i<p->nc;i++){ /* dÃ©termine c=la plus grosse cc */
	/* d1-d0=nb de sommets de la cc numÃ©ro i */
	d1=(i+1<p->nc)? p->d[p->R[i+1]] : NF;
	if(d1-d0>n){ n=d1-d0; c=i; } /* n=nb de sommets de cc max */ 
	d0=d1;
      }
      c=p->C[p->R[c]]; /* c=couleur de la composante max */
      NALLOC(int,T,n); /* T=sommet de GF Ã  garder */
      for(i=j=0;i<NF;i++) if(p->C[i]==c) T[j++]=i;
      free_param_dfs(p); /* p ne sert plus Ã  rien */
      graph* C=ExtractSubgraph(GF,T,n,1); /* construit la cc max */
      free(T); /* T ne sert plus Ã  rien */
      CHECK=CHECK_OFF; /* pour ne pas restocker le graphe */
      PrintGraph(C);
      free_graph(C);
    }break;

    case CHECK_SUBDIV:{
      int n=MEM(CPARAM,0,int);
      if(n<0) n=-n*nb_edges(GF);
      graph* C=subdivision(GF,n,VARIANT); // VARIANT=0 par dÃ©faut
      CHECK=CHECK_OFF; /* pour ne pas restocker le graphe */
      PrintGraph(C);
      free_graph(C);
    }break;

    case CHECK_BFS:{
      param_bfs *p=bfs(GF,MEM(CPARAM,0,int),NULL);
      int t,c;
      printf("root=%i\n",p->root);
      printf("rad[%i]=%i\n",p->root,p->radius);
      printf("cycle[%i]=%i%s\n",p->root,p->cycle,(p->cycle<0)?" (undefined)":"");
      printf("#vertices traversed=%i\n",p->n);
      printf("distance:\n");
      for(k=0;k<=p->radius;k++){
	printf(" d=%i:",k);
	for(c=t=0;t<NF;t++) if(p->D[t]==k){ c++; printf(" %i",t); }
	printf("  (Ã—%i)\n",c);
      }
      printf("vertices not connected to the source:");
      for(t=c=0;t<NF;t++) if(p->D[t]<0) printf(" %i",t), c++;
      if(c==0) printf(" none\n");
      else printf("  (Ã—%i)\n",c);
      printf("parent:");
      for(k=0;k<NF;k++) printf(" %i",p->P[k]);
      printf("\n");
      free_param_bfs(p);
    }break;

    case CHECK_DFS:
    case CHECK_NCC:{
      if(NF<=0) printf("Empty graph\n");
      int m,s=(CHECK==CHECK_DFS)? MEM(CPARAM,0,int) : 0;
      if(CHECK==CHECK_DFS) printf("source: %i\n",s);
      param_dfs *p=dfs(GF,s,NULL);
      printf("#component: %i%s\n",p->nc,(p->nc==1)?" (connected)":"");
      printf("#cut-vertex: %i%s\n",p->na,
	     ((p->nc==1)&&(p->na==0)&&(NF>2))?" (biconnected)":"");
      if(p->nc>1) PrintDistribution(p->C,NF,-2,"- composant");
      if(CHECK==CHECK_NCC) goto check_ncc;

      printf("root:");
      for(i=0;i<p->nc;i++) printf(" %i",p->R[i]);
      if(p->na) printf("\ncut-vertex:");
      for(i=0;i<NF;i++) if(p->A[i]) printf(" %i",i);
      if(p->nc>1){
	printf("\ncomponent:");
	for(i=0;i<NF;i++) printf(" %i",p->C[i]);
      }

      printf("\nparent:");
      for(i=m=0;i<NF;i++){
	if(p->P[i]<0) printf(" -");
	else printf(" %i",p->P[i]);
	m=max(m,p->H[i]); // max depth
      }

      printf("\ndepth: %i",m);
      if(m==NF-1) printf(" (hamiltonian path!)");
      PRINTN;
      PrintDistribution(p->H,NF,-2,"- depth");
      
      /*
      // dÃ©tail des hauteurs
      int h,c;
      for(h=0;h<=m;h++){
	printf(" h=%i:",h);
	for(i=c=0;i<NF;i++) if(p->H[i]==h) printf(" %i",i), c++;
	printf("  (Ã—%i)\n",c);
      }
      */

    check_ncc:
      free_param_dfs(p);
    }break;

    case CHECK_BELLMAN:{
      if(!InitWeights(GF,IW_POS)) Erreur(6);
      int u;
      double dk,sk,d,t,s=0;
      param_bellman *p=Bellman_Ford(GF,MEM(CPARAM,0,int),NULL);
      printf("source=%i\n",p->source);
      for(k=0;k<NF;k++){
	if(p->dist[k]==DBL_MAX) p->dist[k]=INFINITY; // p->dist[k]=inf
	printf("dist[%i]=%lf \tparent=%i",k,p->dist[k],p->parent[k]);
	if(POS){ // compare Ã  distance euclienne si graphe gÃ©omÃ©trique
	  dk=hypot(GF->xpos[p->source]-GF->xpos[k],GF->ypos[p->source]-GF->ypos[k]);
	  sk=(dk==0)? 1 : p->dist[k]/dk;
	  if(sk>s) s=sk,u=k,d=dk; // vrai une fois car s=0 au dÃ©part
	  printf(" \tstretch=%lf (=%lf/%lf)",sk,dk,p->dist[k]);
	}
	printf("\n");
      }
      t=p->dist[u]; // ici stretch s = s(k->u) = d/t 
      if(POS){ // compare Ã  distance euclienne si graphe gÃ©omÃ©trique
	printf("maximum stretch: %lf",s);
	if(d==0) printf("\n"); else printf(" (=%lf/%lf)\n",t,d);
	printf("maximum stretch pair: %i->%i\n",u,p->source);
	printf("maximum stretch path: %i",u);
	u=p->parent[u];
	while(u>=0){
	  printf("->%i",u);
	  u=p->parent[u];
	}
      }
      printf("\n");
      free_param_bellman(p);
    }break;

    case CHECK_STRETCH:
      /* Les appels multiples Ã  Bellman-Ford sont optimisÃ©s et
         nÃ©cessitent que les distances soient symÃ©triques ce qui est
         le cas des graphes gÃ©omÃ©triques. Notamment, on calcule un DFS
         de faÃ§on Ã  minimiser les distances entre le sommet visitÃ© au
         temps t et t+1. Le gain est d'environ 20% pour un graphe
         gÃ©omÃ©trique Ã  n=2000. TODO: Il faudrait plutÃ´t faire un
         premier BFS puis de transformer parent[] en un tableau
         order[] de visite grÃ¢ce Ã  MakeTree().

	 Pour trouver un exemple de graphe (gÃ©omÃ©trique alÃ©atoire) Ã  n
	 sommet dont le stretch maximum dÃ©passe sm, on peut faire:

	 n=20; sm=1.97; s=0; i=0; while [[ $(echo "$s <= $sm" | bc) -eq 1 ]]; do s=$(./gengraph dtheta $n 6 -format no -seed $i -check stretch | grep -e "maximum stretch:" | awk '{print $3}'); echo "n=$n, stretch=$s, seed=$i"; i=$((i+1)); done

      */
      if(!InitWeights(GF,IW_GEO)) Erreur(44);
      else{
	//TopChrono(0);
	int u,v,w,a,b,x,um,vm,z;
	double sv,sw,dw,dv,tv,t,d,s=0;
	double dm,tm,sm=-1;
	NALLOC(int,T,NF+1); // T[0..]=chemin du stretch_max (se termine par -1)
	NALLOC(int,C,NF+1); // C[0..]=chemin du stretch_max min (se termine par -1)
	param_bellman *p=new_param_bellman(NF); p->multiple=1; // pour appels multiples

	// pour le DFS
	param_dfs *q=dfs(GF,0,NULL); // DFS dans G depuis 0
	NALLOC(int,order,NF);
	for(u=0;u<NF;u++) order[q->d[u]]=u; // order[t]=u ssi dfs[u]=t
	free_param_dfs(q);
  
	for(z=0;z<NF;z++){ // pour tous les sommets dans l'ordre DFS
	  v=order[z]; // pour toutes les sources v
	  Bellman_Ford(GF,v,p);	  
	  sv=0; // sv=stretch max pour la source v
	  for(w=0;w<NF;w++){
	    dw=hypot(GF->xpos[v]-GF->xpos[w],GF->ypos[v]-GF->ypos[w]);
	    if(p->dist[w]==DBL_MAX) p->dist[w]=INFINITY; // p->dist[w]=inf
	    sw=(dw==0)? 1 : p->dist[w]/dw; // sw=stretch(w,v)
	    if(sw>sv) sv=sw,u=w,dv=dw; // vrai une fois car sv=0 au dÃ©part
	  }
	  tv=p->dist[u]; // ici stretch_max pour v: sv = sv(v->u) = dv/tv
	  if(sv>s){ // on a trouvÃ© un stretch supÃ©rieur -> recopie dans T
	    s=sv,a=u,b=v,d=dv,t=tv; // ici s = s(a,b) = d/t
	    T[w=0]=x=u; while(x>=0) T[++w]=x=p->parent[x]; // termine par -1
	  }
	  if((sm<0)||(sv<sm)){ // stretch_max minimum
	    sm=sv,um=u,vm=v,dm=dv,tm=tv; // sm = sm(um->vm) = dm/tm
	    C[w=0]=x=u; while(x>=0) C[++w]=x=p->parent[x]; // termine par -1
	  }
	}
	free(order);
	//printf("%s\n",TopChrono(0));
	free_param_bellman(p);
	
	printf("maximum stretch: %lf",s);
	if(NF==0){ printf(" (empty graph)\n"); free(T); free(C); break; }
	if(d==0) printf("\n"); else printf(" (=%lf/%lf)\n",t,d);
	printf("maximum stretch pair: %i->%i\n",a,b);
	printf("maximum stretch path: %i",T[w=0]);
	while(T[++w]>=0) printf("->%i",T[w]);
	printf("\n");
	free(T);
	
	printf("minimum max. stretch: %lf",sm);
	if(dm==0) printf("\n"); else printf(" (=%lf/%lf)\n",tm,dm);
	printf("minimum max. stretch pair: %i->%i\n",um,vm);
	printf("minimum max. stretch path: %i",C[w=0]);
	while(C[++w]>=0) printf("->%i",C[w]);
	printf("\n");
	free(C);
      }
      break;

    case CHECK_DEG:
      printf("#edges: %i\n",nb_edges(GF));
      PrintDistribution(GF->d,NF,-1,"- degree");
      break;

    case CHECK_DEGENERATE:{
      int *T=Prune(GF,&k);
      printf("Degenerate: %i\n",k);
      for(k=0;k<NF;k++) printf("%i ",T[k]);
      printf("\n");
      free(T);
    }break;

    case CHECK_GCOLOR:{
      int *T=Prune(GF,NULL);
      int *C=GreedyColor(GF,T);
      printf("#colors: %i\n",1+GF->int1);
      PrintMorphism("Coloring (node->color):\n",C,GF->n);
      free(C);
      free(T);
    }break;

    case CHECK_KCOLOR:{
      int k=MEM(CPARAM,0,int);
      int *C=kColor(GF,k);
      if(C==NULL) printf("There is no %i-coloration for this graph.\n",k);
      else{
	printf("#colors: %i\n",1+GF->int1);
	PrintMorphism("Coloring (node->color):\n",C,GF->n);
	free(C);
      }
    }break;

    case CHECK_KCOLORSAT:
      kColorSat(GF,MEM(CPARAM,0,int));
      break;

    case CHECK_KINDEPSAT:
      kIndepSat(GF,MEM(CPARAM,0,int));
      break;

    case CHECK_PS1:  k=0; goto check_ps;
    case CHECK_PS1b: k=1; goto check_ps;
    case CHECK_PS1c: k=2; goto check_ps;
    case CHECK_PS1x: k=3;
    check_ps:;
      path *P=new_path(GF,NULL,NF);
      int v=PS1(GF,P,k);
      printf("#tests: %i\nPS1: %s\n",GF->int1,v?"yes (PS1 for sure)":"no (probably not PS1)");
      free_path(P);
      break;
      
    case CHECK_TWDEG:{
      int *T=Prune(GF,&k);
      printf("treewidth <= %i\n",Treewidth(GF,0));
      printf("treewidth >= %i\n",k);
      free(T);
    }break;

    case CHECK_TW:
      k=Treewidth(GF,1);
      printf("#tests: %i\ntreewidth: %i\n",GF->int1,k);
      break;

    case CHECK_DIAMETER:{
      param_bfs *p;
      int u,x=-1;
      for(u=0;u<NF;u++){
	p=bfs(GF,u,NULL);
	if(p->n<NF) break; /* non connexe */
	x=max(x,p->radius);
      }
      free_param_bfs(p);
      printf("diameter: ");
      if(x<0) printf("%s","+âˆž");
      else printf("%i",x);
      printf("\n");
    }break;

    case CHECK_VOLM:{
      // on supprime tous les arcs dÃ©croissant
      // on aurait pu ensuite faire un simple DFS
      param_bellman *p=new_param_bellman(NF);
      if(!InitWeights(GF,IW_POS)) Erreur(6);
      HalfGraph(GF,1); // u->v avec u<v, attention! GF est modifiÃ©
      //PrintGraphList(GF);
      //PrintGraph(GF);
      NALLOC(int,M,NF); // M[u]=volume monotone de u
      int u,v;
      for(u=0;u<NF;u++){
	Bellman_Ford(GF,u,p);
	// calcule la densitÃ© (nombre d'arcs) de la boule couverte
	// par l'arbre de racine u
	M[u]=GF->d[u];
	for(v=0;v<NF;v++) if(p->parent[v]>=0) M[u] += GF->d[v];
      }
      free_param_bellman(p);
      //PRINTT(M,NF);
      PrintDistribution(M,NF,10,"monotonic volume");
    }break;

    case CHECK_RADIUS:{
      param_bfs *p=new_param_bfs();
      int u,x=NF-1;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->n<NF){ x=-1; break; } /* non connexe */
	x=min(x,p->radius);
      }
      free_param_bfs(p);
      printf("radius: ");
      if(x<0) printf("%s","+âˆž");
      else printf("%i",x);
      printf("\n");
    }break;

    case CHECK_GIRTH:{
      param_bfs *p=new_param_bfs();
      int u,x=1+NF;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->cycle>0) x=min(x,p->cycle);
      }
      free_param_bfs(p);
      if(x>NF) x=-1;
      printf("girth: %i%s\n",x,(x<0)?" (undefined)":"");
    }break;

    case CHECK_PATHS:{ /*  sort tous les chemins de x Ã  y */
      path *P=new_path(GF,NULL,NF); /* chemin vide d'au plus NF sommets */
      P->P[0]=MEM(CPARAM,0,int);           /* sommet dÃ©but */
      P->P[1]=MEM(CPARAM,sizeof(int),int); /* sommet fin */
      if((P->P[0]<0)||(P->P[1]<0)||(P->P[0]>=NF)||(P->P[1]>=NF)) Erreur(37);
      int u,v=NextPath(GF,P,-1); /* initialise le premier chemin */
      while(v){
	for(u=0;u<P->n;u++) printf("%i ",P->P[u]);
	printf("\n");
	v=NextPath(GF,P,0);
      }
      free_path(P);
    }break;

    case CHECK_ISO:{
      string s=ARGV[MEM(CPARAM,0,int)];
      graph* H=File2Graph(s,34);
      int *P=Isomorphism(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(P==NULL) printf("Non-isomorphic.\n");
      else{
	PrintMorphism("Isomorphism G->H:\n",P,NF);
	free(P);
      }
      free_graph(H);
    }break;

    case CHECK_SUB:{
      string s=ARGV[MEM(CPARAM,0,int)];
      graph* H=File2Graph(s,34);
      graph* S=Subgraph(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(S==NULL) printf("G is not a subgraph of H.\n");
      else{
	printf("Subgraph S of H isomorphic to G:\n");
	PrintGraph(S);
	PrintMorphism("Isomorphism S->G:\n",S->pint1,S->n);
	free_graph(S);
      }
      free_graph(H);
    }break;

    case CHECK_MINOR:{
      string s=ARGV[MEM(CPARAM,0,int)];
      graph* H=File2Graph(s,34);
      int *C=Minor(H,GF);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(C==NULL) printf("H is not a minor of G.\n");
      else{
	int c,u;
	printf("Model of H in G:\n");
	for(c=0;c<H->n;c++){ /* pour chaque couleur c */
	  printf("%i -> {",c);
	  for(u=0;u<NF;u++) /* on affiche les sommets de la couleur c */
	    if(C[u]==c) printf(" %i",u);
	  printf(" }\n");
	}
	free(C);
      }
      free_graph(H);
      }break;

    case CHECK_ISUB:{
      string s=ARGV[MEM(CPARAM,0,int)];
      graph* H=File2Graph(s,34);
      int *X=InducedSubgraph(H,GF);
      printf("H: %s\n#tests: %i\n",s,GF->int1);
      if(X==NULL) printf("H is not an induced subgraph of G.\n");
      else{
	int u;
	printf("Vertices of the induced subgraph S:");
	for(u=0;u<H->n;u++) printf(" %i",X[u]);
	for(u=0;u<H->n;u++) GF->pint1[u]=X[u];
	PrintMorphism("\nIsomorphism H->S:\n",GF->pint1,H->n);
	free(X);
      }
      free_graph(H);
      }break;

    case CHECK_INFO:{
      printf("- command: %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libÃ©rer ce pointeur
      printf("- seed: %u\n",SEED);
      printf("- time to generate or load the graph: %s\n",TopChrono(0));
      int *R=SortGraph(GF,1);
      printf("- time to traverse and sort the graph: %s\n",TopChrono(0));
      if(R){
	printf("- simple and undirected: %s\n",R[6]?"yes":"no");
	printf("- geometric (X,Y positions): %s\n",GF->xpos?"yes":"no");
	printf("- edge-weights: %s\n",GF->W?"yes":"no");
	printf("- #nodes: %s\n",millier(NF));
	printf("- #arcs: %s\n",millier(R[2]));
	printf("- #self-loops: %s\n",millier(R[0]));
	printf("- #multi-arcs: %s\n",millier(R[1]));
	printf("- #asymmetric relations: %s\n",millier(R[3])); 
	printf("- #node IDs < 0: %s\n",millier(R[4])); 
	printf("- #node IDs â‰¥ n: %s\n",millier(R[5])); 
	printf("- maximum degree: %s\n",millier(R[7]));
	printf("- minimum degree: %s\n",millier(R[8]));
	printf("- #isolated nodes: %s\n",millier(R[9]));
      }else printf("- empty graph\n");
      printf("- memory space: %s bytes\n",millier(SizeOfGraph(GF)));
    }break;

    case CHECK_SIMPLIFY:{
      int u,v,d,i;
      SortGraph(GF,0);
      for(u=0;u<NF;u++){
	d=GF->d[u];
	for(i=0;i<d;i++){
	  v=GF->L[u][i];
	  if(((v>u)||((Q->loop)&&(v==u)))&&
	     ((i==0)||(v!=GF->L[u][i-1]))) printf("%i-%i\n",u,v);
	}
      }
      }break;

    case CHECK_RS_CLUSTER:
      RS_Start("cluster",RS_NI_FP,GF);
      
      /* paramÃ¨tre */
      k=MEM(CPARAM,0,int); /* k=paramÃ¨tre */
      if(k==-1) k=ceil(sqrt((double)NF)); /* ici k>=1, toujours */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* k=0: valeur impossible */
      if(VARIANT) printf("- variant: %i\n",VARIANT);
      BARRE;
      
      /* construit, teste, puis libÃ¨re les tables */
      rs_cluster_tables *RT=rs_cluster(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_cluster_length,-1,NULL); /* teste */
      free_rs_cluster_tables(RT); /* libÃ¨re */
      break;

    case CHECK_RS_DCR:{
      RS_Start("dcr",RS_NI_FP,GF);

      /* paramÃ¨tre */
      k=MEM(CPARAM,0,int); /* k=paramÃ¨tre */
      if(k==-1) k=ceil(Minimize(func1,&NF,1,NF,0)); /* ici k>=1 */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* il faut k>0 */
      if(VARIANT) printf("- variant: %i\n",VARIANT);
      BARRE;

      /* construit, teste, puis libÃ¨re les tables */
      rs_dcr_tables *RT=rs_dcr(GF,k);
      routing_test(GF,RT,(rt_length)rs_dcr_length,-1,RT->dist);
      free_rs_dcr_tables(RT);
    }break;

    case CHECK_RS_AGMNT:{
      RS_Start("agmnt",RS_NI_FP,GF);

      /* paramÃ¨tre */
      VARIANT |= 4; // bit-2 a 1 ssi AGMNT
      k=MEM(CPARAM,0,int); /* k=paramÃ¨tre */
      if(k==-1) k=ceil(Minimize(func1,&NF,1,NF,0)); /* ici k>=1 */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* il faut k>0 */
      BARRE;

      /* construit, teste, puis libÃ¨re les tables */
      rs_dcr_tables *RT=rs_dcr(GF,k);
      routing_test(GF,RT,(rt_length)rs_agmnt_length,-1,RT->dist);
      free_rs_dcr_tables(RT);
    }break;

    case CHECK_RS_TZRPLG:{
      RS_Start("tz rplg",RS_L_FP,GF);

      /* paramÃ¨tre */
      double t=MEM(CPARAM,0,double); /* t=paramÃ¨tre (exposant du RPLG) */
      if((VARIANT<2)&&(0<t)&&(t<2)) Erreur(6); /* valeurs impossibles */
      if((VARIANT==2)&&(t<1)) Erreur(6); /* valeurs impossibles */
      printf("- parameter: %g\n",t);
      BARRE;

      /* construit, teste, puis libÃ¨re les tables */
      rs_tzrplg_tables *RT=rs_tzrplg(GF,t);
      routing_test(GF,RT,(rt_length)rs_tzrplg_length,-1,NULL); /* routage */
      free_rs_tzrplg_tables(RT); /* libÃ¨re les tables */
    }break;

    case CHECK_RS_BC:{
      RS_Start("bc",RS_L_FP,GF);

      /* paramÃ¨tre */
      k=MEM(CPARAM,0,int); /* k=paramÃ¨tre */
      printf("- parameter: %i\n",k); /* ici k>=0, toujours */
      if(k<0) Erreur(6); /* k=0 est possible */
      BARRE;
      
      /* construit, teste, puis libÃ¨re les tables */
      rs_bc_tables *RT=rs_bc(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_bc_length,-1,RT->dist); /* teste */
      free_rs_bc_tables(RT); /* libÃ¨re les tables */
    }break;

    case CHECK_RS_HDLBR:{
      RS_Start("hdlbr",RS_NI_FP,GF);
      
      /* paramÃ¨tre */
      k=MEM(CPARAM,0,int); /* k=paramÃ¨tre */
      if(k==-1) k=ceil(sqrt((double)NF)); /* ici k>=1, toujours */
      k=min(k,NF); /* pas plus que le nombre de sommets */
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<=0) Erreur(6); /* k<=0: valeur impossible */
      BARRE;
      
      /* construit, teste, puis libÃ¨re les tables */
      rs_hdlbr_tables *RT=rs_hdlbr(GF,k);
      routing_test(GF,RT,(rt_length)rs_hdlbr_length,-1,NULL); /* routage */
      free_rs_hdlbr_tables(RT); /* libÃ¨re les tables */
    }break;

    default: if(CHECK==CHECK_ON) break;
      Erreur(47); // on ne devrait jamais avoir cette erreur
    }// fin du switch(CHECK)
    
    free(CPARAM),CPARAM=NULL; /* supprime les paramÃ¨tres pour CHECK */
    free(FPARAM),FPARAM=NULL; /* supprime les paramÃ¨tres pour FILTER */
    if(GF!=Q->G) free_graph(GF),GF=NULL; /* supprime le graphe (Ã©vite le double-free) */
  }
  /* fin du "if(CHECK)" */

  free_query(Q),Q=NULL; /* supprime la requÃªte */
  TopChrono(-1); /* libÃ¨re tous les chronos */
  return 0; /* fin de gengraph */
}


/*# ###
GÃ©nÃ©rateur de graphes - v5.3 - Â© Cyril Gavoille - AoÃ»t 2019

USAGE

       gengraph [-options] graph [parameters]


DESCRIPTION

       GÃ©nÃ¨re sur la sortie standard un graphe. Par dÃ©faut le graphe
       est non orientÃ© et affichÃ© selon une liste d'arÃªtes (en texte),
       mais d'autres formats sont possibles: liste d'adjacence, format
       dot de GraphViz, xfig ou pdf par exemple. En paramÃ¨tre figure
       le nom du graphe ainsi que ses paramÃ¨tres Ã©ventuels,
       typiquement le nombre de sommets. La commande appelÃ©e seule
       affiche l'aide sur les options du gÃ©nÃ©rateur. Si les paramÃ¨tres
       d'une option ou d'un graphe sont absents ou remplacÃ©s par "?",
       une aide spÃ©cifique est affichÃ©e. Une console supportant l'UTF8
       est prÃ©fÃ©rable.

       Ex: gengraph -help
	   gengraph -list | sort
	   gengraph -not ?
	   gengraph -xy unique ?
	   gengraph tree ?
	   gengraph ? arbre
	   gengraph tutte
           gengraph hypercube 8
           gengraph mesh 7 3 -not
	   gengraph mesh 50 50 -dele .5 -maincc -visu
	   gengraph rdodecahedron -visu
           gengraph tree 100 -visu
	   gengraph web 10 3 -visu
	   gengraph gabriel 50 -caption "Grabriel with n=50" -visu
	   gengraph gabriel 2000 -xy seed 1 0.15 -visu
	   gengraph gabriel 700 -xy seed 1 -0.3 -visu
	   gengraph sierpinski 7 3 -visu
	   gengraph udg 400 .1 -visu
	   gengraph udg 400 .1 -xy seed 3 1.5 -visu
	   gengraph udg 400 -1 -vsize -vcolor deg -visu
	   gengraph arytree 6 3 3 -dot filter circo -visu
	   gengraph dyck -dot filter circo -visu
	   gengraph ringarytree 4 2 3 0 -label 1 -visu
	   gengraph arboricity 100 2 -vcolor degr -visu
	   gengraph prime 6 -directed -loop 0 -visu
	   gengraph aqua 3 2 1 . -label 1 -dot filter dot -visu
	   echo "0->1->2->0" | gengraph load - -check bfs 0
	   gengraph tutte | gengraph -filter - diameter p
           gengraph rplg 300 3 -maincc -vcolor degr -vcolor pal wz -vsize -visu
           gengraph -xy box 15 15 -xy round 0 -xy grid 16 rng 30 -visu
	   gengraph linial 7 3 -check kcolorsat 3 | ./glucose -model


   LE FORMAT STANDARD

       Le format par dÃ©faut (ou standard) est une liste d'arÃªtes ou de
       chemins Ã©crits en texte simple. Ce format minimaliste est trÃ¨s
       proche de celui du format "dot" de GraphViz.  D'autres formats
       de sortie sont possibles, notamment le format "dot" (voir
       l'option -format). Les sommets sont numÃ©rotÃ©s consÃ©cutivement
       de 0 Ã  n-1 oÃ¹ n est le nombre de sommets prÃ©sents dans le
       graphe (en fait cela peut Ãªtre changÃ© avec l'option -shift).
       Une arÃªte entre i et j est reprÃ©sentÃ©e par i-j, un arc de i
       vers j par i->j. Les sommets isolÃ©s sont simplement reprÃ©sentÃ©s
       par le numÃ©ro du sommet suivit d'un espace ou d'un retour de
       ligne. Le nombre de sommets du graphe est l'entier le plus
       grand + 1, et s'il y a i-j (ou i->j), alors il existe une arÃªte
       (ou un arc) entre les sommets i et j.

       Pour une reprÃ©sentation plus compacte, les arÃªtes (ou arcs)
       consÃ©cutives d'un chemin du graphe peuvent Ãªtre regroupÃ©es en
       blocs i-j-k-â€¦. Par exemple, les deux arÃªtes 3-5 et 5-8
       peuvent Ãªtre regroupÃ©es en 3-5-8. Mais ce n'est pas
       obligatoire.  Ã‰galement, les arÃªtes (ou arcs) d'une Ã©toile
       peuvent Ãªtre groupÃ©es avec i-(j k â€¦). Par exemple, 3-(5 7 8)
       reprÃ©sente les arÃªtes 3-5, 3-7 et 3-8. Il n'est pas possible
       cependant de combiner chemins et Ã©toiles, comme 3-(5-7-8) ou
       3-(5-(7 8)). Toutefois 3-5-(7 â€¦) est correct, mais pas 3-(5
       6)-7 ni (3 5)-6. Les sommets isolÃ©s et les arÃªtes (ou les blocs
       d'arÃªtes) sont sÃ©parÃ©s par des espaces ou des sauts de ligne.
       Une boucle sur un sommet i est codÃ©e par i-i. Les arÃªtes
       multiples sont codÃ©es par la rÃ©pÃ©tition d'une mÃªme arÃªte, comme
       par exemple i-j i-j, ou encore i-j-i (mÃªme convention pour les
       arcs i->j->i).

       Ex: 0 1-2-3-1            0   1
                                   / \
                                  3â”€â”€â”€2

       reprÃ©sente un graphe Ã  4 sommets, composÃ© d'un sommet isolÃ© (0)
       et d'un cycle Ã  trois sommets (1,2,3). Une reprÃ©sentation
       graphique possible est donnÃ©e Ã  droite.

       Ex: 4-2-1-0-3-2-5

       reprÃ©sente un graphe Ã  6 sommets composÃ© d'un cycle de longueur
       4 et de deux sommets de degrÃ© 1 attachÃ© Ã  2. On aurait pu coder
       le mÃªme graphe avec l'expression 2-(1 3 4 5) 1-0-3. En voici
       une reprÃ©sentation graphique possible:

       Ex:        1
                 / \
              4-2   0
               / \ /
              5   3

       Plus gÃ©nÃ©ralement, une famille de graphes peut Ãªtre dÃ©finie en
       prÃ©cÃ©dant chaque graphe par "[n]" oÃ¹ n est un entier unique
       reprÃ©sentant l'identifiant du graphe.

       Ex: [17] 0-1 [22] 0->1->2->0

       reprÃ©sente une famille composÃ©e de deux graphes: un chemin Ã 
       deux sommets (d'identifiant 17) ainsi d'un cycle orientÃ© Ã 
       trois sommets (d'identifiant 22).

   COMMENT FONCTIONNE LE GÃ‰NÃ‰RATEUR ?

       Pour chaque graphe une fonction adj(i,j) est dÃ©finie. Elle
       fournit d'adjacence (0 ou 1) entre les sommets i et j, des
       entiers entre 0 et n-1. Le graphe est affichÃ© en gÃ©nÃ©rant
       toutes les paires {i,j} possibles et en appelant adj(i,j) (ou
       tous les couples (i,j) possibles dans le cas orientÃ©). Les
       graphes sont ainsi gÃ©nÃ©rÃ©s de maniÃ¨re implicite. Les arÃªtes du
       graphe ne sont pas stockÃ©es en mÃ©moire, mais affichÃ©es Ã  la
       volÃ©e. Ceci permet de gÃ©nÃ©rer des graphes de trÃ¨s grande taille
       sans nÃ©cessiter O(nÂ²) espace de mÃ©moire centrale. Pour certains
       graphes cependant, comme les arbres, certains graphes
       d'intersections ou les graphes gÃ©omÃ©triques, une structure de
       donnÃ©es en O(n) peut Ãªtre utilisÃ©e. Pour les formats
       d'affichage liste, matrix, et smatrix une structure de donnÃ©es
       de taille linÃ©aire (en O(n+m) oÃ¹ m est le nombre d'arÃªtes) est
       utilisÃ©e en interne. Ces trois derniers formats sont donc Ã 
       Ã©viter. Pour la gÃ©nÃ©ration de trÃ¨s grand graphe, le format
       standard ou dot doit Ãªtre privilÃ©giÃ©.

   COMMANDES EXTERNES

       Le programme fait appel, pour certaines fonctions, aux
       commandes systÃ¨mes suivantes qui doivent Ãªtre installÃ©es: sed,
       grep, awk, more, sort, dot.


OPTIONS


....-help [word]
....? [word]
....[option|graph] ?

       Affiche l'aide en ligne qui est contenue dans le fichier source
       du gÃ©nÃ©rateur. Pour cela, le code source .c doit Ãªtre dans le
       mÃªme rÃ©pertoire que l'exÃ©cutable. Si "word" est prÃ©cisÃ©, alors
       les options et noms de graphe contenant "word" sont affichÃ©s.
       La variante "[option|graph] ?" affiche une aide dÃ©taillÃ©e sur
       une option ou un graphe prÃ©cis.
....
       Ex: gengraph ? arbre
           gengraph ktree ?
	   gengraph ? hedron
	   gengraph ? planaire
....
       La forme ? peut ne pas fonctionner correctement si un fichier
       d'un seul caractÃ¨re existe dans le rÃ©pertoire courant (Ã  cause
       de l'interprÃ©tation du shell). Il faut alors utiliser '?' au
       lieu de ?.

....-list
....
       Affiche la liste des graphes et leurs paramÃ¨tres qu'il est
       possible de gÃ©nÃ©rer. Sont listÃ©s d'abord les graphes de bases
       puis les composÃ©s. On obtient une aide sur un graphe
       particulier si son nom est suivit de " ?" ou si ses paramÃ¨tres
       Ã©ventuelles sont absents (dans ce cas il doit Ãªtre le dernier
       mot de la commande).
....
       Ex: gengraph -list | sort
	   gengraph gabriel

....-version
....
       Affiche la version courante du gÃ©nÃ©rateur (en fait du programme
       source), un rÃ©el > 1. Pour cela, le code source .c doit Ãªtre
       dans le mÃªme rÃ©pertoire que l'exÃ©cutable.

....-directed
....-undirected
....
       L'option -directed produit le graphe orientÃ© en testant les nÂ²
       arcs possibles, l'option -undirected permettant de revenir Ã  la
       situation par dÃ©faut, soit en testant les n(n-1)/2 arÃªtes. Les
       boucles peuvent Ãªtre testÃ©es ou pas grÃ¢ce Ã  l'option -loop. En
       format standard ou dot, un arc apparaÃ®t comme i->j au lieu de
       i-j ou i--j pour une arÃªte. Tous les graphes ne sont pas
       forcÃ©ment dÃ©finis pour fonctionner comme espÃ©rÃ© avec l'option
       -directed car certaines fonctions d'adjacence supposent que i<j
       et l'option par dÃ©faut (-undirected). Avec -directed, la
       plupart des graphes vont apparaÃ®tre comme orientÃ©s symÃ©triques.
....
       Ex: gengraph clique 5 -directed
           gengraph cycle 5 -directed -visu

....-not
....
       Inverse la fonction d'adjacence, et donc affiche le complÃ©ment
       du graphe. Cette option est prioritaire sur l'option -redirect.

....-loop p
....
       Supprime (p=0), autorise (p=1) ou force (p=2) les boucles du
       graphe gÃ©nÃ©rÃ©. L'option par dÃ©faut pour les graphes
       non-orientÃ©s est p=0, et pour les graphes orientÃ©s c'est
       p=1. Ces options doivent Ãªtre placÃ©es aprÃ¨s -(un)directed car
       modifiÃ©es par celles-ci.

....-dele p
....
       Permet de supprimer chaque arÃªte du graphe gÃ©nÃ©rÃ©e avec
       probabilitÃ© p.

....-delv p
....
       Similaire Ã  -dele p mais concerne les sommets. Le sommet et ses
       arÃªtes incidentes sont alors supprimÃ©s. Si p est un entier <0,
       alors exactement |p| sommets sont supprimÃ©s. Si k sommets sont
       supprimÃ©s, alors le nom des sommets restant est dans
       l'intervalle [0,n-k[ oÃ¹ n est le nombre de sommets initial du
       graphe. Les noms des sommets sont donc Ã©ventuellement
       renumÃ©rotÃ©s. Voir aussi les options -permute et -shift. Bien
       sÃ»r la fonction d'adjacence est appliquÃ©e sur les noms (i,j)
       originaux.

....-redirect p
....
       Redirige chaque arÃªte uniformÃ©ment avec probabilitÃ© p. Plus
       prÃ©cisÃ©ment, si {i,j} est une arÃªte du graphe original G, alors
       avec probabilitÃ© p l'arÃªte affichÃ©e est {i,k} au lieu de {i,j}
       oÃ¹ k est un sommet choisi uniformÃ©ment parmi les sommets du
       graphe G. Si l'arÃªte {i,j} est supprimÃ©e par -dele ou si le
       sommet i est supprimÃ© par -delv, la redirection n'a pas lieu.
       Cette option est appliquÃ©e aprÃ¨s l'option -not. Le graphe G
       tient donc compte de -not avant de rediriger ses arÃªtes.

....-star n
....
       Ajoute n sommets pendant (degrÃ© 1) aux sommets du graphe. Si
       n>0, alors n reprÃ©sente le nombre total de sommets ajoutÃ©s,
       chacun des n sommets Ã©tant connectÃ©s alÃ©atoirement uniformÃ©ment
       aux sommets du graphe original. Si n<0, alors |n| sommets sont
       ajoutÃ©s Ã  chacun des sommets du graphe. [Cette option n'est
       plus effective depuis v4.5. Elle le sera de nouveau dans une
       prochaine version.]

....-apex n
....
       Ajoute n sommets universels, donc connectÃ©s Ã  tous les sommets
       du graphe. [Cette option n'est plus effective depuis v4.5. Elle
       le sera de nouveau dans une prochaine version.]

....-seed s
....
       Permet d'initialiser le gÃ©nÃ©rateur alÃ©atoire avec la graine s,
       permettant de gÃ©nÃ©rer plusieurs fois la mÃªme suite alÃ©atoire.
       Par dÃ©faut, la graine est initialisÃ©e par une combinaison du
       numÃ©ro de processus de la commande et le temps, donc gÃ©nÃ¨re par
       dÃ©faut des suites diffÃ©rentes Ã  chaque lancement. Le gÃ©nÃ©rateur
       est initialisÃ© lorsque l'option est lue sur la ligne de
       commande. Le comportement du programme peut donc Ãªtre affectÃ©
       suivant son ordre d'apparition. Cependant le graphe est gÃ©nÃ©rÃ©
       aprÃ¨s l'analyse de la ligne de commande.

....-width m
....
       Limite Ã  m le nombre d'arÃªtes et de sommets isolÃ©s affichÃ©s par
       ligne. Cette option n'a pas de signification particuliÃ¨re en
       dehors des formats standard et dot. Par exemple, -width 1
       affiche une arrÃªte (ou un sommet isolÃ©) par ligne. L'option
       -width 0 affiche tout sur une seule ligne. La valeur par dÃ©faut
       est 12.

....-shift s
....
       Permet de renumÃ©roter les sommets Ã  partir de l'entier s
       positif. La valeur par dÃ©faut est -shift 0.  L'intÃ©rÃªt de cette
       option est de pouvoir rÃ©aliser des unions de graphes simplement
       en renumÃ©rotant les sommets et en concatÃ©nant les fichiers aux
       formats standard ou list. Cette option n'a pas d'effets pour
       les formats de sortie de type matrice.

....-permute
....
       Permute alÃ©atoirement uniformÃ©ment le nom des sommets
       lorsqu'ils sont affichÃ©s (ou chargÃ©s en mÃ©moire dans le cas
       d'options -check), et donc aprÃ¨s la gÃ©nÃ©ration du graphe. Les
       numÃ©ros restent dans l'intervalle initial qui, sauf si l'option
       -shift a Ã©tÃ© utilisÃ©e, est [0,n[ oÃ¹ n est le nombre de sommets
       du graphe rÃ©ellement gÃ©nÃ©rÃ©. Voir aussi l'option -label.

....-header
....
       Affiche un prÃ©ambule donnant certaines informations sur le
       graphe, sous forme de commentaire Ã  la C++ (//). Par dÃ©faut
       aucun prÃ©ambule n'est affichÃ©. Les informations affichÃ©es sont:
       - l'heure, la date et la graine du gÃ©nÃ©rateur alÃ©atoire
       - la ligne de commande qui a produit la gÃ©nÃ©ration du graphe
       - le nombre de sommets, d'arÃªtes, le degrÃ©s maximum et minimum
       Pour les formats standard et dot, le nombre d'arÃªtes (et les
       degrÃ©s min et max) n'est pas dÃ©terminÃ© avant l'affichage du
       graphe. Pour cette raison ces nombres ne sont affichÃ©s qu'aprÃ¨s
       le graphe. Pour n'avoir que les informations sur le graphe,
       utiliser -header avec l'option -format no. Voir aussi -check info.

....-caption title
....
       Permet d'ajouter une lÃ©gende Ã  un graphe. Cette option n'a
       d'effet qu'avec le format dot et ces variantes. Il est possible
       d'affiche la "seed" avec le format %SEED. On ne peut avoir plus
       d'une occurrence du mÃªme format (%SEED) dans cette option. 
....
       Ex: gengraph gabriel 30 -caption ex1 -visu
           gengraph gabriel 30 -caption "Exemple 2" -visu
           gengraph gabriel 30 -caption "graph with seed=%SEED" -visu

....-fast
....
       GÃ©nÃ¨re le graphe sans tester les O(nÂ²) arÃªtes possibles, mais
       en utilisant directement la liste d'adjacence du graphe
       prÃ©alablement gÃ©nÃ©rÃ©e lors de l'initialisation du graphe, comme
       c'est le cas pour le graphe "load". Le rÃ©sultat est une
       gÃ©nÃ©ration du graphe en temps O(n+m) au lieu de O(nÂ²).
       L'utilisation typique est (voir aussi "loadc"):
....
       Ex: gengraph load file -fast -delv 0.3 -check ncc
....
       Cet exemple permet de calculer le nombre de composantes
       connexes sur un sous-graphe contenu dans un fichier le tout en
       temps linÃ©aire. D'autres graphes peuvent supporter une
       gÃ©nÃ©ration rapide si elle est implantÃ©e dans l'initialisation
       de la fonction d'adjacence (pour l'instant seul le graphe
       "load" le supporte).  Certaines options, comme -not et -loop 2,
       n'ont alors plus d'effet en prÃ©sence de -fast. Cependant,
       -permute, -delv, -dele, et d'autres fonctionnent normalement.

....-variant v
....
       Permet de passer un entier v pour contrÃ´ler certaines
       fonctionnalitÃ©s du gÃ©nÃ©rateur. Par exemple, "-variant 1 -check
       routing cluster -1" permettra de calculer une variante du
       schÃ©ma de routage "cluster".

....-check [parameters]
....
       Stocke en mÃ©moire le graphe gÃ©nÃ©rÃ© sous la forme d'une liste
       d'adjacence, et lui applique un algorithme. Le graphe est
       gÃ©nÃ©ralement affichÃ©, sauf pour certaine options comme -check
       maincc ou -check routing. Utiliser "-format no" pour ne pas
       afficher le graphe gÃ©nÃ©rÃ©. Cette option nÃ©cessite un espace
       supplÃ©mentaire en O(n+m) pour le stockage du graphe.
....
       -check info
....
          Affiche quelques caractÃ©ristiques du graphe, aprÃ¨s avoir
          effectuÃ© un tri puis un simple parcours de ses listes
          d'adjacences. Indique, par exemple, si le graphe est
          orientÃ©, s'il contient des boucles, des multi-arÃªtes,
          etc. Le graphe lui-mÃªme n'est pas affichÃ©. Indique aussi
          l'occupation mÃ©moire du graphe, le temps de gÃ©nÃ©ration ou
          de chargement et le temps de parcours.
....
       -check simplify
....
          Permet de simplifier un graphe en supprimant les boucles et
          multi-arÃªtes, ceci par un simple parcours du graphe aprÃ¨s un
          tri de ses listes d'adjacences. Les arÃªtes sont affichÃ©es au
          format standard une par ligne, comme avec -width 1 -format
          standard. Cette option est un moyen rapide de sortir un
          graphe au "bon" format. Elle est insensible aux options
          d'affichage, en particulier -format. Il est possible de ne
          pas supprimer les boucles en utilisant -loop 0.
....
          Ex: gengraph loadc G -check simplify > H
	      gengraph load G -fast -permute -check simplify > H
....
	  La diffÃ©rence avec "gengraph load G > H", qui produit aussi
	  une sortie valide, est qu'avec "gengraph loadc G -check
	  simplify > H" le graphe G est lu et affichÃ© en temps
	  quasi-linÃ©aire sans la (re)gÃ©nÃ©ration de ses arÃªtes en temps
	  O(nÂ²).
....
       -check bfs s
....
          Effectue un parcours en largeur d'abord sur le graphe gÃ©nÃ©rÃ©
          depuis le sommet s. La distribution des distances depuis s
          est affichÃ©e, ainsi que l'arborescence (-1 indique que le
          sommet n'a pas de pÃ¨re). La longueur du plus petit cycle
          passant par s est aussi donnÃ©e. Elle vaut -1 s'il n'existe
          pas.
....
       -check bellman s
....
          Calcule les plus courts chemins depuis le sommet s par
          l'algorithme de Bellman-Ford. Si le graphe est gÃ©omÃ©trique,
          le poids de chaque arÃªte correspond Ã  la distance
          euclidienne entre ses extrÃ©mitÃ©s, sinon il vaut 1 et le
          rÃ©sultat sera similaire Ã  un bfs. Dans le cas gÃ©omÃ©trique,
          l'Ã©tirement maximum depuis s est calculÃ©, ainsi qu'un chemin
          le rÃ©alisant. L'implÃ©mentation Ã  l'aide d'une file prend un
          temps linÃ©aire en pratique.
....
       -check stretch
....
          Calcule, comme -check bellman s, l'Ã©tirement d'un graphe
          gÃ©omÃ©trique depuis chaque source s. On affiche une source
          atteignant l'Ã©tirement maximum, mais aussi une source
          atteignant l'Ã©tirement minimum, ainsi qu'un chemin rÃ©alisant
          ces Ã©tirements. Utilisez l'option "-format no" pour ne pas
          avoir l'affichage de la gÃ©nÃ©ration du graphe.
....
       -check volm
....
          Calcule la distribution du volume monotone des sommets. Le
          volume monotone d'un sommet u est le nombre d'arcs du
          sous-graphe des sommets accessibles depuis u dans le graphe
          oÃ¹ seuls les arcs u->v avec u<v ont Ã©tÃ© gardÃ©s. Cette mesure
          dÃ©pend de la numÃ©rotation des sommets. Ã€ utiliser en
          combinaison dans l'option -permute.
....
       -check ncc
       -check connected
....
          Donne le nombre de composantes connexes, leur taille s'il y
          en a plusieurs, ainsi que le nombre de sommets
          d'articulations (cut-vertex) du graphe. Ces informations
          sont aussi affichÃ©es par -check dfs 0.
....
       -check dfs s
....
          Effectue un parcours en profondeur d'abord de toutes les
          composantes connexes du graphe gÃ©nÃ©rÃ© depuis le sommet
          s. ComplÃ¨te l'affichage de -check ncc en affichant en plus
          l'arborescence (-1 indique une racine) ainsi que la
          distribution de profondeur des sommets.
....
       -check deg
       -check edge
       -check edges
....
          Affiche la distribution des degrÃ©s et le nombre d'arÃªtes du
          graphe.
....
       -check degenerate
....
          Donne la dÃ©gÃ©nÃ©rescence du graphe, ainsi que l'ordre
          d'Ã©limination correspondant des sommets.
....
       -check girth
....
          Donne la maille du graphe dans le cas non orientÃ©. La valeur
          -1 est renvoyÃ©e si le graphe est acyclique, et la valeur 0
          dans le cas orientÃ©.
....
       -check diameter
....
          Calcule le diamÃ¨tre du graphe gÃ©nÃ©rÃ©. Affiche +âˆž pour un
          graphe non connexe.
....
       -check radius
....
          Calcule le rayon du graphe gÃ©nÃ©rÃ©, soit la hauteur minimum
          d'un arbre couvrant. Affiche +âˆž pour un graphe non connexe.
....
       -check gcolor
....
          Donne une borne supÃ©rieure sur le nombre chromatique du
          graphe en utilisant l'heuristique du degrÃ© minimum.
....
       -check kcolor k
....
          Donne une k-coloration du graphe (et la couleur pour chaque
          sommet), si c'est possible. Pour cela une recherche
          exhaustive de toutes les k-colorations est effectuÃ©e. Le
          temps est raisonnable si k=3 et n<20.
....
       -check kcolorsat k
....
          Donne une formulation SAT de la k-coloration du graphe. Il
          s'agit de la formulation multivaluÃ©e classique, un sommet
          pouvant avoir plusieurs couleurs sans que cela nuise Ã  la
          validitÃ© du rÃ©sultat. Les contraintes sont dÃ©crites au
          format Dimacs CNF. On peut alors envoyer le rÃ©sultat Ã  un
          solveur SAT comme MiniSat ou Glucose. Le graphe n'est pas
          affichÃ©, et donc "-format no" n'est pas nÃ©cessaire.
....
          Ex: gengraph linial 6 3 -check kcolorsat 3 | ./glucose -model
....
       -check kindepsat k
....
          Donne une formulation SAT d'un ensemble indÃ©pendant de
          taille k du graphe. Les variables i=1 Ã  n indiquent si le
          sommet numÃ©rotÃ© i-1 est dans la solution ou pas. Les
          contraintes sont dÃ©crites au format Dimacs CNF. On peut
          alors envoyer le rÃ©sultat Ã  un solveur SAT comme MiniSat ou
          Glucose. Le graphe n'est pas affichÃ©, et donc "-format no"
          n'est pas nÃ©cessaire.
....
          Pour le problÃ¨me clique de taille k, il suffit de chercher
          un ensemble indÃ©pendant de taille k pour le complÃ©ment du
          graphe. Et pour le problÃ¨me "vertex cover" de taille k,
          c'est un ensemble indÃ©pendant de taille n-k sur le
          complÃ©mentaire qu'il suffit de chercher.
....
       -check ps1
       -check ps1b
       -check ps1c
       -check ps1x n u_1 v_1 â€¦ u_n v_n
....
          Applique le test ps1 ou l'une de ses variantes (voir -filter
          ps1 pour plus de dÃ©tail sur ce test). Affiche aussi le
          nombre de tests rÃ©alisÃ©s (nombre de paires de sommets et de
          chemins testÃ©s).
....
       -check paths x y
....
          Liste tous les chemins simples entre les sommets x et
          y. N'affiche rien si x et y ne sont pas connectÃ©s. L'ordre
          est dÃ©fini suivant le premier plus court chemins dans
          l'ordre des sommets depuis le sommet x.
....
       -check iso H
....
          Teste si le graphe gÃ©nÃ©rÃ© G est isomorphe Ã  H. Si oui,
          l'isomorphisme de G Ã  H est donnÃ©. Le nombre de tests
          affichÃ©s est le nombre de fois oÃ¹ les graphes sont comparÃ©s,
          la comparaison prenant un temps linÃ©aire en la taille des
          graphes. Plus les graphes sont symÃ©triques (comme un cycle
          ou un hypercube), plus le nombre de tests sera important.
....
          Ex: gengraph linial 4 2 | ./gengraph linialc 4 2 -check iso -
              (prend ~ 500 000 tests pour ce 4-rÃ©gulier Ã  12 sommets)
....
          Tester l'isomorphisme entre deux cycles de 8 sommets
          Ã©tiquetÃ©s alÃ©atoirement prend environ 4 mille tests, et
          entre deux cycles de 12 sommets, 30 millions de tests soit
          9" environ. Pour deux arbres Ã  75 sommets (alÃ©atoires mais
          isomorphes), moins de 20 tests suffisent. En revanche, le
          teste pour des graphes arÃªtes et sommets transitifs Ã  16
          sommets, comme "gpetersen 8 3" et "haar 133", est hors de
          portÃ©e.
....
       -check sub H
....
          Teste si le graphe gÃ©nÃ©rÃ© G est un sous-graphe couvrant de H
          (donc avec le mÃªme nombre de sommets). S'ils ont le mÃªme
          nombre d'arÃªtes, le test est Ã©quivalent Ã  l'isomorphisme. Le
          nombre de tests est le nombre total de fois oÃ¹ deux graphes
          sont comparÃ©s.  On peut tester si H est Hamiltonien en
          prennant pour G un cycle.
....
          Tester un cycle de longueur 12 dans une grille 3x4 prend
          jusqu'Ã  environ 32 millions de tests (parfois bien moins),
          soit au plus 10".
....
       -check minor H
....
          Teste si le graphe G gÃ©nÃ©rÃ© contient H comme mineur. Les
	  graphes peuvent Ãªtre non connexes. S'ils ont le mÃªme nombre
	  de sommets le test est Ã©quivalent Ã  celui du sous-graphe
	  (voir -check sub). Dans le cas positif, un modÃ¨le de H dans
	  G est fourni.
....
          Le principe consiste Ã  contracter des arÃªtes de G, de toutes
	  les maniÃ¨res possibles, et Ã  tester si H est un sous-graphe
	  du graphe contractÃ©. Le nombre de tests affichÃ©s est le
	  nombre de contractions plus le nombre total de tests
	  rÃ©alisÃ©s par les tests de sous-graphe. Pour H=Kâ‚„ il est
	  prÃ©fÃ©rable d'utiliser -check twdeg qui donne < 3 ssi le
	  graphe ne contient pas Kâ‚„ comme mineur.
....
       -check twdeg
....
          Donne une borne supÃ©rieure et infÃ©rieure sur la treewidth du
          graphe. Pour la borne supÃ©rieure, on utilise l'heuristique
          du sommet de degrÃ© minimum que l'on supprime et dont on
          complÃ¨te le voisinage par une clique. En cas d'Ã©galitÃ© (mÃªme
          degrÃ©) on sÃ©lectionne le sommet dont il faut rajouter le
          moins d'arÃªtes. La borne infÃ©rieure qui est donnÃ©e provient
          de la dÃ©gÃ©nÃ©rescence. La treewidth est exacte si 0,1 ou 2
          est retournÃ©. L'algorithme est en O(nÂ²).
....
       -check tw
....
          Calcule la treewidth du graphe en analysant tous les ordres
          d'Ã©liminations. La complexitÃ© est donc en n!. Il ne faut
          l'utiliser que si le nombre de sommets est < 12 (Ex:
          gengraph random 12 .5 -check tw donne 5 en environ 750
          millions de tests). Parfois, l'utilisation de -permute peut
          accÃ©lÃ©rer le traitement, car partir d'un ordre d'Ã©limination
          dÃ©jÃ  bon permet d'Ã©liminer rapidement beaucoup d'ordres
          possibles.
....
       -check maincc
....
          Affiche, dans le mode standard seulement, le graphe
          correspondant Ã  la composante connexe ayant le plus grand
          nombre de sommets. Le graphe initial n'est pas affichÃ©. Les
          sommets sont renumÃ©rotÃ©s si le graphe initial n'Ã©tait pas
          connexe. Attention ! l'affichage de la composante n'est
          sensible qu'Ã  l'option -width. En particulier il n'est pas
          possible d'afficher la composante dans un autre format
          (-format) ou avec les noms originaux (-label). Cependant,
          avec "-check maincc | ./gengraph load -" on peut afficher le
          graphe dans le format souhaitÃ©, ou ajouter -visu. (Voir
          aussi le raccourcis -maincc.)
....
	  Notez que "-check maincc -visu" provoque une erreur, car
          -visu applique l'option "-format dot" incompatible avec
	  -check maincc.
....
       -check subdiv n
....
          Affiche, dans le mode standard seulement, une subdivision
          uniforme du graphe telle que chaque arÃªte possÃ¨de n nouveaux
          sommets. Les variantes suivantes sont possibles,
          sÃ©lectionnables via l'option -variant v, m Ã©tant le nombre
          d'arÃªtes du graphe initial:
....
	  â€¢ v=0: subdivision uniforme, chaque arÃªte comprenant n
	    nouveaux sommets. C'est la valeur par dÃ©faut.
....
	  â€¢ v=1: subdivision comprenant un total de n nouveaux sommets
	    rÃ©partis alÃ©atoirement uniformÃ©ment sur les m arÃªtes du
	    graphe.
....
	  â€¢ v=2: comme ci-dessus sauf que chaque arÃªte est subdivisÃ©e
	    au moins une fois (il faut donc que n>=m).
....
	  Si n<0, alors c'est Ã©quivalent Ã  mettre |nÂ·m| comme
          paramÃ¨tre. Notez que "-check subdiv n -visu" provoque une
          erreur, car -visu applique l'option "-format dot"
          incompatible avec -check subdiv.
....
	  Ex: gengraph clique 4 -check subdiv -30 | gengraph load - -visu
....
       -check routing [hash h] [scenario [nomem] s] scheme [parameters]
....
          Construit les tables de routage pour le graphe selon le
          schÃ©ma de routage "scheme", ce schÃ©ma pouvant comporter des
          paramÃ¨tres spÃ©cifiques. La sortie consiste en statistiques
          sur les tables (taille, temps de calcul) et le graphe.
          L'option "scenario" permet en plus de tester certains types
          (s) de routage sur le graphe et d'afficher des statistiques
          sur les longueurs de routes gÃ©nÃ©rÃ©es (dont l'Ã©tirement).
          L'option "hash" permet de prÃ©ciser la fonction de hashage
          (h) appliquÃ©e le cas Ã©chÃ©ant aux sommets souvent utilisÃ©
          dans les schÃ©mas de type "name-independent". Le graphe doit
          Ãªtre connexe et comporter au moins une arÃªte, propriÃ©tÃ©s qui
          sont toujours testÃ©es.
....
          Ex: gengraph -permute rplg 200 2.3 -maincc > G1
	      gengraph loadc G1 -check routing scenario all cluster -1
....
          L'option "nomem" aprÃ¨s "scenario" permet d'optimiser la
          mÃ©moire en ne stockant pas les distances calculÃ©es pour
          Ã©tablir l'Ã©tirement (par dÃ©faut elles le sont). En
          contrepartie le temps de calcul est allongÃ©. Cela peut avoir
          un intÃ©rÃªt si le graphe est trÃ¨s grand (n â‰ƒ 1.000.000
          sommets) et qu'un scenario comme "pair -1000" est testÃ©s.
          Les scenarii possibles sont (n=nombre de sommets du graphe):
....
          â€¢ scenario none     â†’ aucun routage (scenario par dÃ©faut)
          â€¢ scenario all      â†’ n(n-1) routage possibles
          â€¢ scenario edges    â†’ tous les routages entre voisins
          â€¢ scenario npairs   â†’ n paires alÃ©atoires de sommets diffÃ©rents
          â€¢ scenario one u    â†’ n-1 routages depuis u (choix alÃ©atoire si -1)
          â€¢ scenario pair u v â†’ routage de u Ã  v (choix alÃ©atoire si -1)
          â€¢ scenario pair -p  â†’ routage depuis p>0 paires alÃ©atoires
          â€¢ scenario until s  â†’ routage jusqu'Ã  l'Ã©tirement s ou plus
....
          Les fonctions de hachages h:[0,n[->[0,k[ possibles sont:
	  (shuffle et mod atteignent le nombre de collisions minimum de âŽ¡ n/kâŽ¤)
....
          â€¢ hash prime   â†’ h(x)=((aÂ·x+b)%p)%k oÃ¹ 0<a,b<p sont alÃ©atoires
	                   et p=2^31-1 est premier (hash par dÃ©faut)
          â€¢ hash mix     â†’ h(x)=mix(a,b,x)%k oÃ¹ a,b sont alÃ©atoires sur 32 bits
	                   et mix() la fonction mÃ©lange de Bob Jenkins (2006)
          â€¢ hash shuffle â†’ h(x)=ðœ‹(x)%k oÃ¹ ðœ‹(x) est une permutation de [0,n[
	                   basÃ©e sur deux entiers alÃ©atoires de [0,n[.
          â€¢ hash mod     â†’ h(x)=(x+a)%k oÃ¹ a est alÃ©atoire dans [0,k[.
....
       -check routing cluster k
....
          SchÃ©ma de routage name-independent "cluster" de paramÃ¨tre k.
          Un sommet de degrÃ© maximum est choisi comme "centre", puis
          k-1 de ces voisins de plus haut degrÃ©s sont choisis pour
          former un cluster de taille au plus k. Un arbre BFS est
          enracinÃ© depuis le centre. Chaque sommet possÃ¨de une boule
          par rayon croissant qui s'arrÃªte avant de toucher un sommet
          du cluster. On route de u vers v d'abord dans la boule de
          u. Sinon on va dans le cluster pour chercher un sommet du
          cluster responsable du hash de v. Une fois atteint on route
          selon l'arbre BFS ou selon un plus court chemin si la
          distance Ã  la destination est â‰¤ logn/loglogn. Les sommets
          ayant un voisin dans le cluster et qui ne sont pas eux-mÃªmes
          dans le cluster possÃ¨dent dans leur table tout leur
          voisinage. Les boules sont optimisÃ©es par l'usage d'un
          voisin par dÃ©faut.
....
	  Si k=-1, alors k est initialisÃ© Ã  âŽ¡ âˆšnâŽ¤. Si k=-2 il est
          initialisÃ© Ã  n si bien que le cluster est fixÃ© Ã  tout le
          voisinage du centre. Dans tous les cas, l'Ã©tirement est
          toujours â‰¤ 5. Il est mÃªme â‰¤ 3 si k=1. La taille des tables
          est rÃ©duit dans le cas de graphes power-law (comme RPLG). Il
          existe plusieurs variantes (option -variant v) oÃ¹ chacun des
          bits de v mis Ã  1 est interprÃ©tÃ© comme suit:
....
          â€¢ bit-0: le routage est rÃ©alisÃ© sans les boules de
            voisinages ce qui rÃ©duit la taille moyenne mais augmente
            l'Ã©tirement maximum.
....
          â€¢ bit-1: le routage dans le cluster est rÃ©alisÃ©e dans
            l'Ã©toile couvrant le cluster, sans les autres arÃªtes du
            cluster.  Cela rÃ©duit la taille des tables sans modifier
            l'Ã©tirement maximum. Si de plus k=1, le routage est alors
            rÃ©alisÃ© via la racine de l'arbre BFS ce qui rÃ©duit au
            minimum la taille moyenne des tables (2 en moyenne).
....
          â€¢ bit-2: les tables des sommets voisins du cluster sont
	    vides.  L'Ã©tirement devient â‰¤ 7 au lieu de 5 au maximum,
	    mais la taille des tables est rÃ©duite.
....
          â€¢ bit-3: les tables des sommets voisins du cluster ne
            contiennent que des sommets qui ne sont pas de le cluster.
            L'Ã©tirement â‰¤ 5 est prÃ©servÃ©, mais gÃ©nÃ©ralement la taille
            maximum des tables est rÃ©duite, l'Ã©tirement moyen
            augmentant que trÃ¨s lÃ©gÃ¨rement.
....
       -check routing dcr k
....
          SchÃ©ma de routage name-independent "dcr" de paramÃ¨tre k>0
          reprÃ©sentant le nombre de couleurs. C'est une simplification
          du schÃ©ma "agmnt". L'Ã©tirement est toujours â‰¤ 5 et le nombre
          d'entrÃ©es des tables est en moyenne f(k,n) = 2n/k +
          (k-1)Â·(H(k)+1) oÃ¹ H(k) ~ ln(k)+0.577â€¦ est le k-iÃ¨me nombre
          harmonique. Le principe du schÃ©ma est le suivant.  Chaque
          sommet possÃ¨de une couleur, un entier alÃ©atoire de [0,k[,
          les sommets landmarks Ã©tant ceux de couleur 0. Les boules de
          voisinages des sommets sont dÃ©finies par volume comme la
          plus petite boule contenant au moins chacune des couleurs
          (ou tous les sommets s'il manque une couleur n'est pas
          reprÃ©sentÃ©e), les sommets du dernier niveau Ã©tant ordonnÃ©s
          par identifiant croissant. Le routage sâ†’t s'effectue selon
          un plus court chemin si t est dans la boule de s ou si t est
          un landmark.  Sinon, on route vers le sommet w de la boule
          de s dont la couleur est Ã©gale au hash de t, une valeur
          aussi dans [0,k[. Puis le routage wâ†’t s'effectue dans
          l'arbre BFS enracinÃ© dans le plus proche landmark de s ou de
          t, celui minimisant la distance de w Ã  t.
....
	  Si k=-1, alors k est initialisÃ© Ã  sa valeur optimale
          thÃ©orique, celle qui minimise le nombre moyen d'entrÃ©es
          f(k,n), valeur calculÃ©e numÃ©riquement et qui vaut environ k
          â‰ƒ âˆš(n/ln(n))/2, ce qui donne environ 2âˆš(nÂ·ln(nÂ·ln(n)))
          entrÃ©es en moyenne.  Si k=-2, le nombre de couleur est
          initialisÃ© Ã  n. Les valeurs de k>n sont possibles, il s'agit
          alors d'un routage de plus court chemins comme pour le cas
          k=n ou k=1. Il existe une variante (-variant 1) lorsque k>1
          qui a pour effet de choisir les landmarks comme les sommets
          de plus haut degrÃ©. Plus prÃ©cisÃ©ment, les âŽ¡ n/kâŽ¤ sommets de
          plus haut degrÃ© sont coloriÃ©s 0 et les autres coloriÃ©s
          alÃ©atoirement dans [1,k[. Dans cette variante, la borne sur
          l'Ã©tirement est toujours garantie mais plus sur le nombre
          maximum d'entrÃ©es.  Cependant, pour certains graphes
          l'Ã©tirement moyen est amÃ©liorÃ©.
....
       -check routing agmnt k
....
          SchÃ©ma de routage name-independent dit "agmnt" du nom de ces
	  auteurs Abraham et al. (2008). C'est la version originale du
	  schÃ©ma "dcr" qui diffÃ¨re par l'algorithme de routage. Comme
	  dans "dcr", le routage sâ†’t s'effectue directement vers t si
	  t est dans la boule de s ou est un landmark. Sinon, vers le
	  sommet w de la boule de s dont la couleur est Ã©gale au hash
	  de t. Le routage wâ†’t s'effectue suivant la meilleure des
	  options suivantes: router via un arbre BFS d'un des
	  landmarks; ou bien router via un arbre couvrant la boule
	  d'un sommet s' contenant Ã  la fois w et un sommet x voisin
	  d'un sommet y contenu dans la boule de t (les boules s' et
	  de t, si elles existent, sont dites "contiguÃ«s via l'arÃªte
	  x-y"). L'Ã©tirement est toujours â‰¤ 3 et les tables ont le
	  mÃªme nombre d'entrÃ©es que "dcr", bien que plus complexes. Le
	  temps de calcul des tables est plus important que pour
	  "dcr". Toutes les variantes de "dcr" (-variant, k<0)
	  s'appliquent aussi Ã  "agmnt".
....
       -check routing tzrplg t
....
          SchÃ©ma de routage Ã©tiquetÃ© inspirÃ© de celui de Thorup &
          Zwick (2001) optimisÃ© pour les Random Power-Law Graphs (voir
          prlg) de paramÃ¨tre rÃ©el t (power-law exponent) et proposÃ©
          par Sommer et al. (2012). L'Ã©tirement est toujours â‰¤ 5. Les
          valeurs de t entre ]0,1.5] sont interdites.  Le schÃ©ma
          utilise des sommets landmarks oÃ¹ des arbres BFS sont
          enracinÃ©s, ainsi que des boules (de voisinage) dÃ©finies par
          rayons croissant qui s'arrÃªte avant de toucher un
          landmark. Le routage s'effectue alors en prioritÃ© via les
          boules ou alors via le landmark le plus proche de la
          destination (sans raccourcis), information prÃ©cisÃ©e dans
          l'Ã©tiquette de la destination. Les landmarks sont les
          sommets de plus haut degrÃ©. Par dÃ©faut (-variant 0) leur
          nombre vaut:
....
            â€¢ si t>1.5, âŽ¡ n^((t-2)/(2t-3))âŽ¤
            â€¢ si t=0,   âŽ¡ âˆšnâŽ¤
            â€¢ si t<0,   |t|
....
	  Si -variant 1 et t>1.5, alors les landmarks sont tous les
	  sommets de degrÃ© > n^1/(2t-3). Si -variant 2 et t>0, alors
	  les landmarks sont t sommets choisis alÃ©atoirement.
	  L'Ã©tirement est â‰¤ 3 si un seul landmark est choisi.
....
       -check routing hdlbr k
....
          SchÃ©ma de routage name-independent HDLBR selon Tang et
          al. (2013) avec k>0 landmarks qui sont les sommets de plus
          haut degrÃ©.  Si k<0, alors k est initialisÃ© Ã  âŽ¡ âˆšnâŽ¤. Chaque
          sommet qui n'est pas un landmark possÃ¨de une boule dont le
          rayon est juste infÃ©rieure au plus proche landmark. Chaque
          sommet possÃ¨de sa boule boule inverse (qui peut Ãªtre vide),
          dÃ©finie comme l'ensemble des sommets le contenant dans leur
          boule. Chaque landmark a une couleur unique. On route
          directement de u Ã  v si v est dans la boule de u, la boule
          inverse de u, ou si v est un landmark. Sinon, on route vers
          le landmark (selon un plus court chemin) dont sa couleur
          vaut le hash de v. DelÃ  on route suivant un plus court
          chemin vers le plus proche landmark de v, l(v). On utilise
          ensuite le next-hop de l(v) vers v, un sommet nÃ©cessairement
          contenant v dans sa boule inverse. Ã€ chaque Ã©tape du
          routage, si v est dans la boule ou boule inverse du sommet
          courant, on route directement vers celui-ci. La longueur de
          route entre u et v est au plus 2d(u,v)+2r, oÃ¹ r est la
          distance maximum entre deux landmarks. La valeur de r est
          bornÃ©e par une constante pour k â‰ƒ âˆšn et pour les Random
          Power-Law Graphs (voir rplg).
....
       -check routing bc k
....
          SchÃ©ma de routage Ã©tiquetÃ© selon Brady-Cowen (2006).
          L'Ã©tirement est additif â‰¤ 2k, et donc multiplicativement â‰¤
          2k+1. En particulier, si k=0, il s'agit d'un routage de plus
          court chemin. Il est adaptÃ© aux Power-Law Random Graphs
          (voir plrg). Le principe est le suivant: on construit un
          arbre BFS (=T) enracinÃ© dans un sommet (=r) de plus haut
          degrÃ©. Le coeur (=C) est la boule de rayon k depuis r. On
          construit une liste (=L) de BFS couvrant G ainsi qu'une
          forÃªt (=H) de BFS de G comme suit. Au dÃ©part, L={T}, et H
          est la forÃªt T\C. Puis, pour chaque arÃªte {u,v} de G\C\T, on
          vÃ©rifie si l'ajout de {u,v} Ã  H crÃ©e un cycle ou pas. Si
          c'est non, on met Ã  jour la forÃªt H en lui ajoutant
          {u,v}. Si c'est oui, on calcule un BFS de G de racine u (ou
          v) qu'on ajoute Ã  L.  Une fois le calcul de H terminÃ©, on
          calcule une forÃªt BFS couvrante de H qu'on ajoute Ã  L.
          L'algorithme de routage de u Ã  v consiste simplement Ã 
          router dans l'arbre de L qui minimise la distance de u Ã  v.

....-filter family[:range] [not] test [parameters]
....
       Affiche les graphes d'une famille pour lesquels le test est
       vrai (ou faux si "test" est prÃ©cÃ©dÃ© de "not"). Le paramÃ¨tre
       "family" est un nom de fichier ou "-" pour l'entrÃ©e standard.
       La lecture de la famille se fait en temps linÃ©aire.  Il
       contient la famille de graphes (ou un graphe seul) au format
       standard. La premiÃ¨re ligne affichÃ©e contient le nombre de
       graphe de la famille rÃ©sultante.  L'affichage de chaque graphe
       est influencÃ© par l'option -width qui doit Ãªtre placÃ©e avant
       -filter. La variante "family:range" permet de sÃ©lectionner les
       graphes de la famille dont les identifiants sont spÃ©cifiÃ©s par
       "range", comme par exemple "family:5-8" qui sÃ©lectionne les
       graphes d'identifiant 5,6,7,8. De maniÃ¨re gÃ©nÃ©rale, "range" est
       un ensemble de valeurs selon le format "value" dÃ©crit ci-aprÃ¨s
       (voir aussi -filter F id value). La variante "-:range" est
       possible.
....
       Dans la suite, la prÃ©sence de "value" dans les paramÃ¨tres d'une
       option reprÃ©sente un ensemble de valeurs entiÃ¨res possibles.
       Par exemple, -filter F vertex '>5' filtre les graphes de la
       famille F comportant plus de 5 sommets. De maniÃ¨re gÃ©nÃ©rale,
       "value" est une suite d'intervalles d'entiers sÃ©parÃ©s par des
       "," (interprÃ©tÃ©es comme "ou"), chaque intervalle Ã©tant codÃ©
       comme suit:
....
          â€¢ <x       â†’  valeur infÃ©rieure Ã  l'entier x
          â€¢ >x       â†’  valeur supÃ©rieure Ã  l'entier x
          â€¢ x ou =x  â†’  valeur Ã©gale Ã  l'entier x
          â€¢ x-y      â†’  valeur dans l'ensemble d'entiers {x,â€¦,y}
          â€¢ t        â†’  toujours vrai (intervalle infini)
          â€¢ p        â†’  affiche la valeur plutÃ´t que le graphe
....
       Ex: -filter F vertex '5,7-13,>100'
           -filter F vertex '5-10,p'
           -filter F edge p
           -filter F id 5,7
           -filter F vertex p | head -1
           -filter â€¦ p | grep '^\[' | sort -rnk 3 | head
           -filter â€¦ p | grep '^\[' | sort -nk 3 | head
....
       Le premier exemple filtre les graphes de la famille F ayant un
       nombre de sommets n vÃ©rifiant soit n=5, soit 7 â‰¤ n â‰¤ 13, ou
       soit n>100. L'exemple 2 affiche le nombre de sommets des
       graphes ayant entre 5 et 10 sommets.  L'exemple 3 affiche le
       nombre d'arÃªtes de chaque graphe.  L'exemple 4 affiche les
       graphes d'identifiant 5 et 7 de la famille F. L'exemple 5
       affiche le nombre de graphes de la famille (ce qui correspond Ã 
       la premiÃ¨re ligne de commentaires). Les deux derniers exemples
       permettent d'avoir le maximum/minimum de "value".
....
       Si "value" contient le symbole > ou < il est alors prÃ©fÃ©rable
       de mettre des quotes ('>14' par exemple) pour que la commande
       soit correctement interprÃ©tÃ©e par le shell.
....
       La diffÃ©rence principale avec -check est que le rÃ©sultat de
       -filter est non verbeux alors que -check, qui ne s'applique pas
       a priori sur des familles de graphes mais sur un graphe seul,
       peut donner des explications sur l'exÃ©cution de l'algorithme
       dont la sortie n'est pas forcÃ©ment un ou une liste de
       graphes. Aussi, avec -check l'algorithme s'applique au graphe
       gÃ©nÃ©rÃ©, donc aprÃ¨s un temps a priori en O(nÂ²), alors qu'avec
       -filter c'est toujours Ã  partir d'un fichier (ou de l'entrÃ©e
       standard), lu en temps linÃ©aire.
....
       -filter F id value
....
          Filtre les graphes de F dont l'identifiant est dÃ©terminÃ© par
          value. Cela permet d'extraire un ou plusieurs graphes
          donnÃ©s. C'est Ã©quivalent Ã  "-filter F:value all".
....
       -filter F rename shift
....
          Affiche tous les graphes de la famille en renumÃ©rotant les
          graphes Ã  partir de l'entier "shift".
....
       -filter F vertex value
....
          Filtre les graphes de F ayant un nombre de sommets dÃ©terminÃ©
          par value.
....
       -filter F edge value
       -filter F edges value
....
          Filtre les graphes de F d'un nombre d'arÃªtes dÃ©terminÃ© par
          value.
....
       -filter F all (= vertex t)
....
          Affiche tous les graphes de F ce qui permet en particulier
          de savoir combien il y en a en examinant la premiÃ¨re ligne
          affichÃ©e.
....
       -filter F1 minus F2
....
          Affiche F1\F2, c'est-Ã -dire tous les graphes de F1 qui ne
          sont pas isomorphes Ã  F2 (si F2 est un graphe) ou Ã  l'un des
          graphes de F2 (dans le cas d'une famille).
....
       -filter F1 minus-id F2
....
          Comme "minus" mais concerne les identifiants: supprime de F1
          les graphes dont l'identifiant existe dans F2 (qui peut Ãªtre
          un graphe ou une famille de graphes). La complexitÃ© est
          environ (|F1|+|F2|)Â·log|F2|, alors que pour "minus" elle est
          en |F1|Â·|F2|Â·T oÃ¹ T est le temps pour dÃ©cider si deux
          graphes pris dans F1 et F2 sont isomorphes.
....
       -filter F minor[-inv] H
....
          Filtre les graphes de F contenant H comme mineur. La
          variante minor-inv filtre les graphes de F qui sont mineurs
          de H. Si H=Kâ‚„, il est prÃ©fÃ©rable d'utiliser -filter tw2.
....
       -filter F sub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe,
          chaque graphe de F devant avoir le mÃªme nombre de sommets
          que H. La variante sub-inv filtre les graphes de F qui sont
          un sur-graphe de H.
....
       -filter F isub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe
          induit. La variante isub-inv filtre les graphes de F qui
          sont sous-graphes induits de H.
....
       -filter F iso H
....
          Filtre les graphes de F isomorphes Ã  H.
....
       -filter F degenerate value
....
          Filtre les graphes de F de dÃ©gÃ©nÃ©rescence dÃ©terminÃ©e par
          value.
....
       -filter F forest value
....
          Filtre les graphes de F qui sont des forÃªts dont le nombre
          d'arbres est dÃ©terminÃ© par value.
....
       -filter F isforest (= forest t)
....
          Filtre les graphes de F qui sont des forÃªts.
....
       -filter F istree (= forest '=1')
....
          Filtre les graphes de F qui sont des arbres.
....
       -filter F cycle (= not forest t)
....
          Filtre les graphes de F contenant au moins un cycle.
....
       -filter F degmax/degmin value
....
          Filtre les graphes de F de degrÃ© maximum (ou minimum)
          dÃ©terminÃ© par value.
....
       -filter F deg value
....
          Filtre les graphes de F oÃ¹ tous les sommets ont un degrÃ©
          dÃ©terminÃ© par value. Ainsi -filter deg 4-7 filtre les
          graphes avec un degrÃ© minimum au moins 4 et un degrÃ© maximum
          au plus 7.
....
       -filter F gcolor value
....
          Filtre les graphes de F dont le nombre de couleurs obtenu
          selon l'heuristique du degrÃ© minimum est dÃ©terminÃ© par
          value.
....
       -filter F bipartite (= gcolor <3)
....
          Filtre les graphes de F qui sont bipartis.
....
       -filter F component value
....
          Filtre les graphes de F dont le nombre de composantes
          connexes est dÃ©terminÃ© par value.
....
       -filter F connected (= component 1)
....
          Filtre les graphes de F qui sont connexes.
....
       -filter F biconnected
....
          Filtre les graphes de F qui sont 2-connexes. Un graphe G est
          k-connexe s'il n'y a pas d'ensemble avec <k sommets qui
          dÃ©connecte G ou laisse G avec 1 sommet. Un graphe est
          2-connexe s'il est connexe, ne possÃ¨de pas de sommet
          d'articulation et a plus de 2 sommets. Les cliques de taille
          k+1 sont k-connexes.
....
       -filter F radius value
....
          Filtre les graphes de F dont le rayon est dÃ©terminÃ© par
          value. Le rayon est la profondeur du plus petit arbre
          couvrant le graphe. Il vaut -1 si le graphe n'est pas
          connexe.
....
       -filter F girth value
....
          Filtre les graphes de F dont la maille est dÃ©terminÃ©e par
          value. La maille est la taille du plus petit cycle. Elle
          vaut -1 si le graphe n'a pas de cycle. Elle n'est dÃ©finie
          que si les graphes sont orientÃ©s.
....
       -filter F diameter value
....
          Filtre les graphes de F dont le diamÃ¨tre est dÃ©terminÃ© par
          value. Le diamÃ¨tre vaut -1 si le graphe n'est pas connexe.
....
       -filter F cut-vertex value
....
          Filtre les graphes de F dont le nombre de sommets
          d'articulations est dÃ©terminÃ© par value. Un sommet est un
          point d'articulation si sa suppression augmente le nombre de
          composante connexe. Les sommets de degrÃ© 1 ne sont pas des
          point d'articulation. Le graphe est biconnexe ssi value<1 ou
          si le graphe est une clique avec au moins deux sommets. On
          peut tester si un graphe est une clique avec -filter degmin
          ou -filter deg.
....
       -filter F ps1
       -filter F ps1b
       -filter F ps1c
       -filter F ps1x n u_1 v_1 â€¦ u_n v_n
....
          Filtre les graphes G de la famille F dont le test ps1 est
          vrai, c'est-Ã -dire si l'Ã©valuation de la fonction f(G,{})
          dÃ©crite ci-aprÃ¨s est vraie.
....
	  Soit P un chemin d'un graphe G tel que G\P est connexe. La
          fonction f(G,P) est vraie ssi G est vide (en pratique
          |G|-|P|<3 suffit) ou s'il existe deux sommets x,y de G oÃ¹ y
          n'est pas dans P tels que pour tout chemin Q entre x et y
          dans G "compatible" avec P (c'est-Ã -dire P et Q
          s'intersectent en exactement un segment) on a les deux
          conditions suivantes: (1) il n'y a pas d'arÃªte entre les
          sommets de P\Q et de G\(QâˆªP); et (2) pour toute composante
          connexe C de G\(QâˆªP), f(CâˆªQ,Q) est vraie. Le test est
          optimisÃ© dans un certain nombre de cas, en particulier: les
          arbres (toujours vrai), les cliques (vrai ssi n<5).
....
	  La variante ps1b calcule et affiche de plus un graphe des
          conflits (affichage modifiable par -width), chaque noeud de
          ce graphe correspondant Ã  un argument (CâˆªQ,Q) Ã©valuÃ© Ã  faux
          par f. La valeur (ou code) d'un noeud est 0 (=lourd ou
          faux), 1 (=lÃ©ger ou vrai) ou - (indÃ©terminÃ©e). Suivant
          certaines rÃ¨gles, les valeurs 0 ou 1 sont propagÃ©es selon le
          type des arÃªtes du graphes des conflits.  RÃ©soudre le graphe
          des conflits revient Ã  trouver une affectation des valeurs 0
          ou 1 aux noeuds qui respecte (sans contradiction) toutes les
          rÃ¨gles.
....
	  La fonction f(G,{}) est Ã©valuÃ©e Ã  vraie si le graphe des
          conflits n'a pas de solution, c'est-Ã -dire si une
          contradiction a Ã©tÃ© dÃ©couverte ou si pour une paire de
          sommets (x,y) tous ses noeuds sont Ã  1.
....
	  On affiche le code d'un noeud (0,1,-) ainsi que les sommets
          de sa composante (par ex: [237]).  Les noeuds du graphe des
          conflits sont reliÃ©es par des arÃªtes typÃ©es. Les voisins v
          d'un noeud u sont listÃ©s avec le type de l'arÃªte, si l'un
          des 4 cas suivants se produit (il n'y a pas d'arÃªte entre u
          et v dans les autres cas):
....
	     v<  (la composante de v est incluse dans celle de u)
	     v>  (la composante de v contient celle de u)
	     v=  (les composantes de u et v sont les mÃªmes) 
	     v|  (les composantes de u et v sont disjointes) 
....
	  Parmi les rÃ¨gles on trouve par exemple: si deux noeuds du
          graphe des conflits u=(CâˆªQ,Q) et v=(C'âˆªQ',Q') sont
          disjoints, c'est-Ã -dire C n'intersecte pas C', alors seule
          une des expressions f(CâˆªQ,Q) ou f(C'âˆªQ',Q') peut Ãªtre
          fausse, pas les deux. Dit autrement, les composantes de u et
          v ne peuvent pas Ãªtre "lourdes" (=0) toutes les deux en mÃªme
          temps. Et donc, si le code de u est 0, celui de v est
          1. Notons que le code de u et v Ã©gale Ã  1 est compatible
          avec cette rÃ¨gle.
....
	  La variante ps1c est similaire Ã  ps1b sauf que rÃ©cursivement
          seul le test ps1 est appliquÃ©, et pas ps1b. Le test ps1c est
          plus long que ps1 mais plus rapide que ps1b. La variante
          ps1x est similaire Ã  ps1b sauf que les valeurs v_i sont
          Ã©crites dans le noeuds u_i du graphe des conflits principal
          (pas ceux gÃ©nÃ©rÃ©s lors des appels rÃ©cursifs). Plus
          prÃ©cisÃ©ment, v_1 (0 ou 1) est Ã©crit dans le noeud u_1, puis
          sa valeur est propagÃ©e. Ensuite v_2 est Ã©crit puis propagÃ©e,
          etc.
....
          Dans tous les cas, si G n'est pas connexe, le rÃ©sultat n'est
          pas dÃ©terminÃ©.
....
       -filter F tw value
....
          Filtre les graphes de F selon leur treewidth. L'algorithme
          pour le calcul de la treewidth est assez lent. Pour les
          petites valeurs de tw, des alternatives sont possibles (voir
          -check tw et -filter tw2). Pour savoir si un graphe G est de
          treewidth 3 il suffit de savoir si G contient l'un des 4
          mineurs suivants:
....
          !!! echo "[0]"  > F; ./gengraph clique 5 >> F
	      echo "[1]" >> F; ./gengraph wagner >> F
	      echo "[2]" >> F; ./gengraph prism 5 >> F
	      echo "[3]" >> F; ./gengraph hajos >> F ; echo "0-1-2-0" >> F
	      cat G |./gengraph -filter F minor-inv - -format no
....
       -filter F tw2
....
          Affiche les graphes de F de treewidth â‰¤ 2. L'algorithme est
          en O(nÂ²). Ce test peut Ãªtre utilisÃ© pour tester (plus
          rapidement qu'avec -filter minor) les graphes sans mineur
          Kâ‚„.
....
       -filter F hyper value
....
          Filtre les graphes de F selon leur hyperbolicitÃ©. Il s'agit
          de la valeur (entiÃ¨re) maximum, sur tous les quadruplets de
          sommets {u,v,x,y}, de la diffÃ©rence des deux plus grandes
          sommes parmi les sommes de distance : uv+xy, ux+vy et
          uy+vx. La complexitÃ© est en O(nâ´).

....-format type
....
       SpÃ©cifie le format de sortie. Il est prÃ©fÃ©rable d'utiliser
       cette option en dernier. Les valeurs possibles pour "type"
       sont:
....
       â€¢ standard: format standard (liste d'arÃªtes), c'est le plus compact.
       â€¢ list: liste d'adjacence.
       â€¢ matrix: matrice d'adjacence.
       â€¢ smatrix: matrice supÃ©rieure, diagonale comprise.
       â€¢ vertex i: liste des voisins du sommet i.
       â€¢ dot: format de GraphViz qui est trÃ¨s proche du format standard.
       â€¢ dot<type>: dessine le graphe avec GraphViz et convertit au format <type>.
       â€¢ html: dessin dynamique au format html et vis.js (cf. http://visjs.org).
       â€¢ xy: positions X,Y qui ont Ã©tÃ© utilisÃ©es pour le graphe gÃ©omÃ©trique.
       â€¢ no: n'affiche rien, Ã  utiliser en combinaison avec -header ou -check.
....
       Les formats matrix/smatrix/list/vertex nÃ©cessitent de stocker
       le graphe en mÃ©moire, donc nÃ©cessite un espace en O(n+m), alors
       que le graphe est gÃ©nÃ©rÃ© Ã  la volÃ©e pour les formats standard,
       dot ou html. Pour html, nÃ©cessite le script vis.min.js qui est
       chargÃ© soit localement s'il existe soit sur un dÃ©pÃ´t externe
       officiel sinon. Le script peut Ãªtre lent pour le calcul des
       coordonnÃ©es Ã  partir d'une centaine de sommets. Les valeurs les
       plus utilisÃ©es de <type> pour le format dot<type> sont: pdf,
       fig, svg, ps, jpg, gif, png (voir man dot ou faire dot -T.).
....
       L'option -format dot<type> est Ã©quivalent Ã  "-format dot | dot
       -T<type>". Par consÃ©quent, elle doit donc Ãªtre utilisÃ©e en
       dernier. Le filtre dot utilisÃ© pour dessiner le graphe peut
       Ãªtre spÃ©cifiÃ© par l'option -dot filter.  L'affichage des noms
       de sommets est contrÃ´lÃ© par l'option -label.
....
       Les positions affichÃ©es dans le format dot ([pos="â€¦"])
       diffÃ¨rent d'un facteur proportionnel Ã  âˆšn par rapport aux
       positions originales du graphe (qui peuvent Ãªtre affichÃ©es par
       -format xy ou -label -3). Ce facteur permet de garder une
       taille raisonable pour les sommets car sous dot les sommets ont
       une taille fixe minimale.

....-vcolor option [parameters]
....
       Ces options permettent de modifier la couleur des sommets. Ces
       options n'ont d'effets qu'avec le format dot (et ses variantes
       y compris -visu). Par dÃ©faut les sommets sont de couleur
       noire. Notez que les attributs par dÃ©faut des sommets
       (couleurs, formes, etc.) peuvent Ãªtre modifiÃ©s directement par
       dot (voir l'option -N de dot). Cependant l'option -vcolor
       permet d'individualiser la couleur d'un sommet, en fonction de
       son degrÃ© par exemple. Ici le degrÃ© est le degrÃ© non-orientÃ©.
       Il peut avoir plusieurs options -vcolor pour une mÃªme commande.
....
       -vcolor deg[r]
....
          La couleur dÃ©pend du degrÃ© du sommet (deg) ou du rang du
          degrÃ© du sommet (degr). Ainsi, les sommets de plus petit
          degrÃ© obtiennent la premiÃ¨re couleur de la palette, les
          sommets de plus grand degrÃ© la derniÃ¨re couleur de la
          palette, et les autres sommets une couleur intermÃ©diaire de
          la palette. Donc une seule couleur est utilisÃ©e si le graphe
          est rÃ©gulier.
....
       -vcolor degm
....
          Effectue une coloration propre (deux sommets voisins ont des
          couleurs diffÃ©rentes) suivant l'heuristique du degrÃ©
          minimum: rÃ©cursivement, le sommet de degrÃ© minimum obtient
          la plus petite couleur qui n'est pas utilisÃ©e par ses
          voisins. Cela donne des colorations avec assez peu de
          couleurs pour les graphes de faible arboricitÃ© (planaire,
          tw, pw, kout, expander, â€¦) ou de faible degrÃ©. Avec cette
          technique, les graphes bipartis (tree, crown, â€¦) sont
          coloriÃ©s avec deux couleurs. Cette option nÃ©cessite un
          espace et un temps en O(n+m).
....
       -vcolor randg
....
          Effectue une coloration propre en utilisant un algorithme
          glouton sur un ordre alÃ©atoire des sommets: rÃ©cursivement,
          le sommet d'indice i obtient la plus petite couleur qui
          n'est pas utilisÃ©e par ses voisins d'indice j<i. Cette
          option nÃ©cessite un espace et un temps en O(n+m).
....
       -vcolor kcolor k
....
          Effectue une k-coloration propre du graphe, si c'est
          possible. Si cela n'est pas possible, la premiÃ¨re couleur
          est appliquÃ©e Ã  tous les sommets. L'algorithme (exponentiel)
          est le mÃªme que celui utilisÃ© pour -check kcolor.
....
       -vcolor pal grad
....
          Permet de fixer la palette de couleurs utilisÃ©e par les
          sommets. Le paramÃ¨tre "grad" est un mot sur l'alphabet [a-z]
          (sans les guillemets). Les caractÃ¨res en dehors de cet
          alphabet sont ignorÃ©s. Chaque lettre correspond Ã  une
          couleur de base:
....
       !!! a=aquamarine     h=hotpink      o=olive         v=violet
	   b=blue           i=indigo       p=purple        w=white
	   c=cyan           j=orange       q=pink          x=gray
	   d=darkorange     k=khaki        r=red           y=yellow
	   e=chocolate      l=lavender     s=salmon        z=black
	   f=forestgreen    m=magenta      t=teal
	   g=green (lime)   n=navy         u=yellowgreen
....
          La palette est calculÃ©e selon une interpolation linÃ©aire
          entre les points dÃ©finis par le mot "grad". Par exemple, si
          "grad" vaut rb, la palette sera composÃ©e d'un dÃ©gradÃ© allant
          du rouge (r) au bleu (b). Si "grad" vaut rgbr, le dÃ©gradÃ©
          ira du rouge au vert puis au bleu et enfin au rouge. Pour
          avoir une couleur (de base) unique, disons w, sur tous les
          sommets, poser "grad" Ã©gale Ã  w. Par exemple, pour avoir
          tous les sommets blancs, on peut faire:
....
          !!! gengraph gabriel 30 -vcolor deg -vcolor pal w -visu
....
          La palette par dÃ©faut correspond au mot "grad" suivant:
          redjykugfocatbhsqvmpinzxlw. On peut visualiser la palette
          avec l'option "-vcolor list".
....
       -vcolor list
....
          Produit l'affichage de la palette des couleurs utilisÃ©es
          pour un graphe plutÃ´t que le graphe lui-mÃªme. Cela permet en
          particulier de savoir combien de couleur ont Ã©tÃ© utilisÃ©es.
          La palette est gÃ©nÃ©rÃ©e en affichant au format dot un graphe
          particulier oÃ¹ les sommets (reprÃ©sentÃ©s par un rectangle)
          sont les couleurs utilisÃ©es. Utilisez -visu pour visualiser
          la palette sous forme pdf. Le nom des sommets correspond Ã 
          la lettre de la couleur de base comme spÃ©cifiÃ© par -vcolor
          pal.
....
          Ex1: gengraph gabriel 50 -vcolor degm -vcolor list
	  (gÃ©nÃ¨re la palette utilisÃ©e pour ce graphe de Gabriel)
....
          Ex2: gengraph prime 53 -vcolor list
	  (un moyen simple de gÃ©nÃ©rer la palette par dÃ©faut)
....
          Ex3: gengraph clique 100 -vcolor degm -vcolor pal rb -vcolor list
          (gÃ©nÃ¨re un dÃ©gradÃ© de 100 couleurs allant du rouge au bleu)

....-vsize
....
       La taille des sommets est proportionnelle Ã  son degrÃ©, alors
       que par dÃ©faut elle est fixe. Cette option n'a d'effet qu'avec
       le format dot (et ses variantes). Elle peut Ãªtre combinÃ©e avec
       -vcolor.

....-visu
....-visuh
....
       CrÃ©e un fichier "g.pdf" (ou "g.html" pour -visuh) permettant de
       visualiser le graphe. Il s'agit d'un raccourci de l'option
       "-format dotpdf" ("-format html") qui rajoute Ã©galement la
       redirection "> g.pdf" ("> g.html") en fin de la ligne de
       commande.

....-maincc
....
       Affiche la composante connexe principale du graphe, les sommets
       Ã©tant Ã©ventuellement renumÃ©rotÃ©s si le graphe n'est pas
       connexe. C'est un raccourci pour "-check maincc | ./gengraph
       load - -fast". (Voir aussi -check maincc.) Cet affichage est
       rÃ©alisÃ© en temps linÃ©aire grÃ¢ce Ã  l'option -fast. Les options
       placÃ©es avant -maincc affectent le graphe initial alors que
       celles placÃ©es aprÃ¨s affectent la composante principale.  Les
       options ayant un effet pour les formats hors standard (comme
       -vsize ou -visu) ne devraient Ãªtre placÃ©es qu'aprÃ¨s cette
       option.

....-dot option [parameters]
....
       Cette option permet de controler la sortie au format dot. Elle
       permet par exemple de modifier le filtre, la longueur des
       arÃªtes ou l'Ã©chelle du dessin.
....
       -dot scale s
....
          SpÃ©cifie le facteur d'Ã©chelle pour le format dot. Cela
          affecte les coordonnÃ©es des sommets et des arÃªtes, pas des
          Ã©tiquettes (sommets ou arÃªtes). Cela permet d'Ã©carter les
          sommets les uns des autres si nÃ©cessaires. Le format pour s
          prend plusieurs formes: x ou x,y pour un facteur d'Ã©chelle
          identique ou pas en X et Y. La valeur par dÃ©faut est s=1. On
          peut aussi mettre "auto" qui calcule automatiquement un
          facteur d'Ã©chelle (symÃ©trique en X et Y) qui vaut 1/âˆšn ou
          1/max(Î”X,Î”Y) dans le cas de graphes gÃ©omÃ©triques.
....
          Ex: gengraph gabriel 10 -label -3 -dot scale 3,2 -visu
....
       -dot len p
....
          SpÃ©cifie la longueur des arÃªtes pour le format dot et le
          filtre "neato". La valeur par dÃ©faut est 1, et une valeur
          plus grande (comme 2.5 ou 3) allonge les arÃªtes et permet
          dans certains cas de mieux visualiser le graphe. C'est
          parfois nÃ©cessaire pour Ã©viter l'intersection des sommets
          lorsqu'on utilise -label 1. On peut obtenir le mÃªme genre
          d'effet avec -dot scale.
....
       -dot filter f
....
         SpÃ©cifie le filtre de GraphViz, c'est-Ã -dire l'algorithme de
         dessin utilisÃ© par dot. Par dÃ©faut, le filtre est
         "neato". Les filtres principaux sont: dot, neato, twopi,
         circo, fdp, sfdp, â€¦ Faire "dot -K ." pour afficher les
         filtres disponibles.

....-pos b
....
       Active (b=1) ou dÃ©sactive (b=0) la gÃ©nÃ©ration des positions des
       sommets pour le format dot. Cela sert Ã  indiquer Ã  l'algorithme
       de dessin dot de respecter (b=1) ou pas (b=0) les coordonnÃ©es
       des sommets. L'option par dÃ©faut est -pos 0, mais cette option
       est activÃ©e pour tous les graphes gÃ©omÃ©triques (udg, gabriel,
       thetagone, â€¦).

....-label b
....
       Active (bâ‰ 0) ou dÃ©sactive (b=0) l'affichage du nom des sommets
       pour les formats dot et standard. Les valeurs possibles de
       l'entier b sont bâˆˆ[-3,3]. Si b=1, il s'agit du nom original du
       sommet, par exemple un mot binaire pour l'hypercube. Cette
       fonctionnalitÃ© n'est pas implÃ©mentÃ©e pour tous les graphes, le
       nom par dÃ©faut Ã©tant les entiers de [0,n[ oÃ¹ n est le nombre de
       sommets du graphe gÃ©nÃ©rÃ©. L'option -label 1 -visu permet alors
       d'afficher sur le dessin du graphe le nom des sommets. Ils ne
       le sont pas par dÃ©faut (b=0).  L'option -label 2 -visu force
       l'affichage des noms sous forme d'entiers de [0,n[ et non pas
       du nom original (s'il Ã©tait dÃ©fini). L'option -label 3 permet,
       dans le cas de graphe gÃ©omÃ©trique (ou si -pos 1), d'afficher
       les coordonnÃ©es des points. Si b<0, alors l'effet est similaire
       Ã  -label |b| sauf que le nom du sommet est affichÃ© Ã  cotÃ© du
       sommet et non au centre. L'option -label 1 annule l'option
       -permute, mais -label 2 ne le fait pas. Comme l'option -label
       influence -format dot<type>, -label devrait Ãªtre placÃ©e avant
       -format.
....
       Ex: gengraph petersen -label 1 -width 1
           gengraph petersen -label 1 -format dot | grep label
           gengraph petersen -label 1 -dot len 2 -visu
           gengraph gabriel 30 -pos 0 -label 1 -visu
	   gengraph gabriel 30 -label -3 -dot scale 4 -xy round 2 -visu

....-norm â„“ [parameters]
....
       Fixe la norme d'un vecteur (x,y) du plan (ou la fonction de
       distance entre deux points du plan) pour l'adjacence de
       certains graphes gÃ©omÃ©triques (dont udg, gabriel, rng, nng, â€¦).
       Par dÃ©faut c'est la norme Euclidienne qui est utilisÃ©. Les
       valeurs possibles pour â„“ sont:
....
         â€¢ L1      â†’  |x|+|y|, distance de Manhattan
         â€¢ L2      â†’  âˆš(xÂ²+yÂ²), norme Euclidienne
         â€¢ Lmax    â†’  max{|x|,|y|}
         â€¢ Lmin    â†’  min{|x|,|y|}
         â€¢ poly p  â†’  distance polygonale de paramÃ¨tre p
         â€¢ hyper   â†’  distance hyperbolique
....
       Il s'agit de pseudo-norme (ou pseudo-distance) puisque par
       exemple la norme "Lmin" ne vÃ©rifie pas l'inÃ©galitÃ©
       triangulaire. La norme polygonale est le rayon du cercle
       inscrit dans le polygone rÃ©gulier convexe Ã  p cotÃ©s contenant
       (x,y), le polygone Ã©tant centrÃ© en (0,0) et orientÃ© de faÃ§on Ã 
       avoir son cotÃ© le plus Ã  droit vertical. Ainsi "poly 4"
       correspond Ã  la norme "Lmax". Une valeur de p<3 est interprÃ©tÃ©e
       comme p=+âˆž ce qui correspond Ã  la norme Euclidienne. Attention
       ! la norme "poly p" n'est pas toujours symÃ©trique, lorsque p
       est impair par exemple. La norme (ou distance) hyperbolique
       n'est dÃ©finie que pour des points du disque ouvert unitÃ© centrÃ©
       en (0,0).

....-xy option [parameters]
....
       Cette option contrÃ´le la faÃ§on dont sont gÃ©nÃ©rÃ©es les
       coordonnÃ©es des sommets d'un graphe gÃ©omÃ©trique. Par dÃ©faut les
       positions sont tirÃ©es alÃ©atoirement uniformÃ©ment dans le carrÃ©
       [0,1[Â², mais cela peut Ãªtre changÃ© grÃ¢ce Ã  -xy. Notez bien
       que, mÃªme si c'est improbable, deux sommets peuvent avoir les
       mÃªmes positions (voir l'option -xy unique). Il est possible de
       visualiser les points issus des options -xy (voir le graphe
       "point n"). La prÃ©sence d'une option -xy active l'option -pos
       et rend le graphe gÃ©omÃ©trique comme:
....
       Ex: gengraph path 20 -xy unique -visu
....
       -xy load file
....
          Charge les positions Ã  partir du fichier "file" ou de
          l'entrÃ©e standard si file=-. Cela permet de tester les
          adjacences d'un graphe gÃ©omÃ©trique Ã  partir de positions
          prÃ©-dÃ©terminÃ©es. Le format est celui de -format xy.
....
          Ex: gengraph gabriel 10 -xy load file.pos
....
	  Le nombre de sommets du graphe est dÃ©terminÃ© par le fichier
          et non par les paramÃ¨tres du graphe. Cette option n'a
          d'effet que pour les graphes gÃ©omÃ©triques. La structure du
          fichier texte doit Ãªtre:
....
          !!!    n
		 x_1 y_1
		 x_2 y_2
		 â€¦
		 x_n y_n
....
	  oÃ¹ n est le nombre de positions. Les positions x_i y_i ne
	  sont pas forcÃ©ment dans l'intervalle [0,1[. Notez qu'avec
	  l'option -format xy, il est possible d'effectuer la
	  transformation d'un fichier de positions. L'exemple suivant
	  normalise les coordonnÃ©es du fichier g.pos dans le carrÃ©
	  unitÃ©:
....
          Ex: gengraph -xy load g.pos -xy box 1 1 -format xy
....
       -xy box a b
....
          Effectue un redimensionement des positions de sorte quelles
          se situent dans le rectangle [0,a[ Ã— [0,b[. En prenant
          a=b=1, les coordonnÃ©es seront renormalisÃ©es dans le carrÃ©
          [0,1[Â². Cette opÃ©ration est effectuÃ©e juste avant la
          gÃ©nÃ©ration des arÃªtes, mais aprÃ¨s avoir effectuÃ© l'opÃ©ration
          -xy noise (voir ci-aprÃ¨s) et/ou -xy load.
....
       -xy grid p
....
          Ajoute une grille p Ã— p au graphe gÃ©nÃ©rÃ©, ce qui est utile
          lorsque les coordonnÃ©es des points sont entiers.
          Techniquement, on ajoute au format de sortie dot un
          sous-graphe reprÃ©sentant la grille oÃ¹ les sommets et les
          arÃªtes sont de couleur grise. Si p<0, alors le paramÃ¨tre est
          initialisÃ© Ã  1+âŽ£ âˆšnâŽ¦ ou bien Ã  n si l'option "-xy
          permutation" est prÃ©sente, n Ã©tant le nombre de sommets du
          graphe. Pour Ãªtre visible, le nombre de lignes (et de
          colonnes) de la grille gÃ©nÃ©rÃ©e doit Ãªtre au moins 2.
....
       -xy zero
....
          Ajoute l'origine (0,0) au dessin qui est reprÃ©sentÃ© par un
          cercle rouge.
....
       -xy vsize f
....
          Facteur de grossissement des sommets pour le format dot. Par
          dÃ©faut f=1.
....
       -xy noise r p
....
          Effectue une perturbation alÃ©atoire sur les positions des
	  sommets. Le dÃ©placement de chaque sommet est effectuÃ© dans
	  sa boule de rayon r (pour p>0) selon une loi en puissance de
	  paramÃ¨tre p. Prendre p=0.5 pour une perturbation uniforme
	  dans cette boule, p>0.5 pour une concentration des valeurs
	  vers le centre et p<0.5 pour un Ã©cartement du centre. Les
	  valeurs <0 de p donne des Ã©cartements au delÃ  du rayon r.
....
          Plus prÃ©cisÃ©ment, une direction (angle de 0 Ã  2ðœ‹) est
	  choisie alÃ©atoirement uniformÃ©ment, puis, selon cette
	  direction, un dÃ©calage alÃ©atoire est effectuÃ© selon une loi
	  en puissance: si x est uniforme dans [0,1[, le dÃ©calage sera
	  d(x)=rÂ·x^p.  AprÃ¨s cette opÃ©ration, il est possible que les
	  points ne soient plus dans le rectangle d'origine, ce qui
	  peut bien sÃ»r Ãªtre corrigÃ© par -xy box.
....
       -xy seed k p
....
          GÃ©nÃ¨re les points Ã  partir (ou autour) de k>0 graines. Les
	  graines sont choisies uniformÃ©ment dans le carrÃ© [0,1[Â² puis
	  centrÃ©es par rapport Ã  leur barycentre. Chaque point est
	  alors tirÃ© alÃ©atoirement autour d'une des graines et Ã  une
	  distance variant selon une loi en puissance (voir -xy noise)
	  de paramÃ¨tre p et de rayon r â‰ƒ âˆš(ln(k+1)/k). Ce rayon
	  correspond au seuil de connectivitÃ© pour un Unit Disk Graph
	  Ã  k sommets dans le carrÃ© [0,1[Â² (voir udg n r). On peut
	  obtenir une distribution uniforme dans un disque avec -xy
	  seed 1 0.5 (voir aussi -xy disk) sauf que le centre est en
	  (1/2,1/2) au lieu de (0,0) comme avec -xy disk.
....
          Ex: gengraph point 1000 -xy seed 1 1
	      gengraph point 1000 -xy seed 1 0.5
....
       -xy permutation
....
          GÃ©nÃ¨re les points correspondant Ã  une permutation ðœ‹
          alÃ©atoire uniforme. Le point i aura pour position (i,ðœ‹(i)).
....
       -xy mesh x y
....
          GÃ©nÃ¨re tous les points de coordonnÃ©es entiÃ¨res correspondant
          aux sommets d'une grille de x colonnes et de y lignes.
....
       -xy cycle
....
          GÃ©nÃ¨re les points rÃ©guliÃ¨rement espacÃ©s le long d'un cercle
	  de centre (0,0) et de rayon 1. Les points sont ordonnÃ©es
	  selon l'angle de leurs coordonnÃ©es polaires.
....
          Ex: gengraph cycle 10 -xy cycle -visu
....
       -xy unif
....
          GÃ©nÃ¨re les points alÃ©atoirement uniformÃ©ment dans le carrÃ©
          [0,1[Â². C'est la distribution par dÃ©faut.
....
       -xy circle
....
          GÃ©nÃ¨re les points alÃ©atoirement uniforme le long d'un cercle
	  de centre (0,0) et de rayon 1. Les points sont ordonnÃ©es
	  selon l'angle de leurs coordonnÃ©es polaires.
....
       -xy disk
....
          GÃ©nÃ¨re les points alÃ©atoirement uniforme dans le disque
          unitÃ© de centre (0,0) triÃ©s selon l'angle de leurs
          coordonnÃ©es polaires. Cette distribution permet de gÃ©nÃ©rer,
          par exemple, un polygone "star-shaped". La distribution est
          similaire Ã  l'option "-xy seed 1 0.5" sauf que les points
          sont ordonnÃ©es.
....
          Ex: gengraph cycle 25 -xy disk -visu
....
          Remarque: les points sont gÃ©nÃ©rÃ©s avant l'application des
          options comme -xy round, -xy noise, ou -xy unique qui
          modifient les coordonnÃ©es et qui peuvent donc produire des
          croisements avec le graphe cycle par exemple.
....
       -xy hyper p
....
          GÃ©nÃ¨re les points alÃ©atoires selon une loi exponentielle de
          paramÃ¨tre p dans le disque unitÃ© de centre (0,0). Les points
          sont triÃ©s selon l'angle de leurs coordonnÃ©es polaires.
....
       -xy convex
       -xy convex2
....
          GÃ©nÃ¨re les points alÃ©atoirement en position convexe Ã 
          l'intÃ©rieur d'un cercle de rayon 1 et de centre (0,0), ce
          qui peut Ãªtre modifiÃ© par -xy ratio. Ils sont numÃ©rotÃ©s
          consÃ©cutivement selon le parcours de l'enveloppe convexe. On
          les gÃ©nÃ¨re comme suit, n Ã©tant le nombre points Ã 
          gÃ©nÃ©rer. Inductivement, une fois que n-1 points en position
          convexe ont Ã©tÃ© gÃ©nÃ©rÃ©s, on choisit un angle ð›¼ alÃ©atoire du
          cercle de rayon 1 et de centre (0,0) supposÃ© Ã  l'intÃ©rieur
          du convexe. On dÃ©termine ensuite la partie S du segment
          d'angle ð›¼ oÃ¹ chacun des points de S=[a,b[ forment avec les
          n-1 points prÃ©cÃ©dant un convexe. Enfin, on choisit
          alÃ©atoirement un point de S selon la probabilitÃ© âˆš|b-a| pour
          obtenir n points en position convexe. Les angles des trois
          premiers points sont choisis parmi trois secteurs non
          adjacents d'angle ðœ‹/3 si bien que l'origine est toujours Ã 
          l'intÃ©rieur de l'ensemble convexe.
....
          La variante -xy convex2 gÃ©nÃ¨re Ã©galement des points
          alÃ©atoires en position convexe, selon la mÃ©thode suivante.
          On gÃ©nÃ¨re n points alÃ©atoires u_i du carrÃ© [0,1[Â² puis on
          calcule les vecteurs diffÃ©rences v_i â‰¡ u_{i+1}-u_i (mod
          n). Les vecteurs (dont la somme est nulle) sont ensuite
          triÃ©s par angle croissant, puis les points en position
          convexe sont obtenus de proche en proche en ajoutant chacun
          des vecteurs v_i. Cette mÃ©thode tend Ã  gÃ©nÃ©rer des points
          proches d'un cercle, chaque angle et chaque longueur entre
          deux points consÃ©cutifs suivant une loi normale.
....
          Ex: gengraph cycle 25 -xy convex -visu
	      gengraph dtheta 100 6 -xy convex -visu
....
          L'ordre des sommets peut Ãªtre modifiÃ© par certaines options
          (voir la remarque de l'option -xy disk).
....
       -xy polygon p
....
          GÃ©nÃ¨re des points alÃ©atoires uniformÃ©ment dans un polygone
	  convexe rÃ©gulier Ã  pâ‰¥3 cotÃ©s inscrit dans le cercle de
	  centre (0,0) et de rayon 1 de sorte qu'un des cotÃ©s du
	  polygone soit vertical. Les sommets ne sont pas
	  spÃ©cifiquement ordonnÃ©s. Pour une distribution uniforme dans
	  un disque, soit lorsque p=+âˆž, utiliser -xy disk. L'option
	  pour p=4 est similaire Ã  -xy unif, sauf que pour p=4 la
	  distribution est dans le carrÃ© [-c,+c[ Ã— [-c,+c[ oÃ¹ c =
	  cos(ðœ‹/4) = Â½âˆš2 â‰ƒ 0.707â€¦ au lieu du carrÃ© [0,1[Â².
....
       -xy ratio r
....
          Modifie les distributions de points faisant intervenir une
          forme de largeur 1 et de hauteur r, comme: -xy unif, -xy
          circle, -xy cycle, -xy convex, -xy disk, -xy seed. La valeur
          par dÃ©faut est r=1. Le rÃ©el r>0 est donc le ratio de la
          hauteur par la largeur de la forme. Par exemple, pour la
          distribution par dÃ©faut (-xy unif), les points seront
          alÃ©atoires uniformes dans le rectangle [0,1[ Ã— [0,r[. Si la
          forme est un cercle (-xy circle ou -xy disk), alors la forme
          devient une ellipse dont le rayon horizontal est 1 et celui
          vertical r. Dans le cas de -xy seed, les graines sont alors
          gÃ©nÃ©rÃ©s dans le rectangle [0,1[ Ã— [0,r[.
....
       -xy surface s
....
          DÃ©finit la signature s de la surface sur laquelle va Ãªtre
          construit le graphe gÃ©omÃ©trique. La surface peut-Ãªtre
          orientable ou non, avec ou sans bord. Elle est reprÃ©sentÃ©e
          par un polygone convexe rÃ©gulier inscrit dans un cercle de
          rayon 1 et dont les 2|s| cotÃ©s sont appariÃ©s. Cette option
          se charge Ã©galement de gÃ©nÃ©rer des points alÃ©atoirement
          uniformes sur la surface. La signature est un mot s sur
          l'alphabet {h,c,b} de longueur |s|=2g, oÃ¹ g est le genre de
          la surface, indiquant comment sont appariÃ©s les 4g cotÃ©s du
          polygone. Chaque cotÃ© est appariÃ© avec le cotÃ© +2 (le
          suivant du suivant) ou le cotÃ© -2 selon l'une des trois
          coutures suivantes:
....
          â€¢ h = handle   = couture orientÃ©e ou anse
          â€¢ c = crosscap = couture non-orientiÃ©
          â€¢ b = border   = aucune couture
....
          La caractÃ©ristique d'Euler de la surface (ou sa courbure)
          vaut 2-|s| = 2-2g. Certaines signatures ont des synonymes.
          Par exemple, "-xy surface torus" est synonyme de "-xy
          surface hh" (voir ci-dessous leurs listes).
....
          Ex: -xy surface bb (ou plane ou square)  â†’ plan rÃ©el
              -xy surface hb (ou cylinder) ....... â†’ cylindre
              -xy surface cb (ou mobius) ......... â†’ ruban de MÃ¶bius
              -xy surface hh (ou torus) .......... â†’ tore
	      -xy surface ch (ou klein) .......... â†’ bouteille de Klein
              -xy surface cc (ou projective) ..... â†’ plan projectif
	      -xy surface hhhh ................... â†’ double tore
....
	  Cette option active Ã©galement -xy polygon 4g et -xy ratio 1
	  pour gÃ©nÃ©rer des points alÃ©atoires uniformÃ©ment sur la
	  surface.
....
       -xy round p
....
          Arrondi les coordonnÃ©es Ã  10^-p prÃ¨s. Il faut que p soit un
          entier < DBL_DIG, soit p<15 en gÃ©nÃ©ral. Donc p=0 arrondi Ã 
          l'entier le plus proche. Cet opÃ©rateur est appliquÃ© aprÃ¨s
          -xy box. Il sert aussi Ã  prÃ©ciser le nombre de dÃ©cimales Ã 
          afficher pour l'option -format xy (par dÃ©faut p=6). Par
          exemple, la combinaison -xy box 100 100 -xy round -1 permet
          d'avoir des coordonnÃ©es multiples de 10.
....
       -xy unique
....
          Supprime les sommets en double, correspondant aux mÃªmes
          positions. Cela peut Ãªtre utile lorsqu'on utilise -xy round
          par exemple. Cette opÃ©ration est appliquÃ©e aprÃ¨s toutes les
          autres, notamment aprÃ¨s -xy box et -xy round. Ceci est
          rÃ©alisÃ© Ã  l'aide d'un tri des points, l'ordre n'est donc pas
          prÃ©servÃ©).


   GRAPHES

       Deux types de graphes sont possibles : les graphes de base et
       les graphes composÃ©s. Ces derniers sont obtenus en paramÃ©trant
       un graphe de base. Une catÃ©gorie importante de graphes sont les
       graphes gÃ©omÃ©triques (qui peuvent Ãªtre composÃ©s ou de bases).
       L'adjacence est dÃ©terminÃ©e par les coordonnÃ©es associÃ©es aux
       sommets. De nombreuses options s'y rÃ©fÃ¨rent.  Ils activent tous
       par dÃ©faut l'option -pos. Les graphes orientÃ©s activent quant Ã 
       eux tous l'option -directed.
       


   GRAPHES DE BASE :

....grid n_1 â€¦ n_k .
....
       Grille Ã  k dimensions de taille n_1 Ã— â‹¯ Ã— n_k. Si la taille n_i
       est nÃ©gative, alors cette dimension est cyclique. Par exemple,
       "grid -10 ." donnera un cycle Ã  10 sommets.

....ring n c_1 â€¦ c_k .
....
       Anneaux de cordes Ã  n sommets chacun ayant k cordes de longueur
       c_1,â€¦,c_k. Par exemple, "ring 10 1 ." donnera un cycle Ã  10
       sommets. Chaque cáµ¢ peut Ãªtre positif ou nÃ©gatif, mais il faut
       cáµ¢â‰¥-n. Si k=0, il s'agit d'un stable Ã  n sommets.  L'option
       -directed permet d'obtenir une k-orientation.

....cage n c_1 â€¦ c_k .
....
       Graphe pouvant servir Ã  la construction de graphes n-cage,
       c'est-Ã -dire aux plus petits graphes cubiques Ã  n sommets de
       maille donnÃ©e. Ils sont toujours Hamiltoniens. Ils peuvent Ãªtre
       vus comme des anneaux de cordes irrÃ©guliers. Ils sont
       construits Ã  partir d'un cycle de longueur n dÃ©coupÃ© en n/k
       intervalles de k>0 sommets. Le i-Ã¨me sommet de chaque
       intervalle, disons le sommet numÃ©ro j du cycle, est adjacent au
       sommet numÃ©ro j+cáµ¢ du cycle (modulo n). Chaque cáµ¢ peut Ãªtre
       positif ou nÃ©gatif, mais il faut cáµ¢â‰¥-n. Il est aussi possible
       de construire des graphes avec des sommets de degrÃ© 4 comme
       "cage 8 0 2 .", le graphe de ChvÃ¡tal, ou avec des sommets de
       degrÃ© 2 comme "cage 4 2 0 .".

....arboricity n k
....
       Graphe d'arboricitÃ© k Ã  n sommets alÃ©atoire. Ce graphe est
       composÃ© de l'union de k>0 arbres alÃ©atoires. Il est donc
       toujours connexe. Chacun des arbres est un arbre plan enracinÃ©
       alÃ©atoire uniforme dont les sommets sont permutÃ©s
       alÃ©atoirement, sauf le premier arbre dont les sommets sont
       numÃ©rotÃ©s selon un parcours en profondeur. Ces graphes
       possÃ¨dent au plus kÂ·(n-1) arÃªtes, et pour k=1 il s'agit d'un
       arbre. L'option -directed permet d'obtenir une k-orientation.

....rarytree n b z
....
       Arbre b-aire plan alÃ©atoire uniforme Ã  n noeuds internes. Il
       faut bâ‰¥2. Il possÃ¨de bn+1+z sommets, z Ã©tant un paramÃ¨tre
       valant 0 ou 1. La racine est de degrÃ© b+z, les autres sommets
       sont de degrÃ© b+1 (soit b fils) ou 1 (=feuille). Les sommets
       sont numÃ©rotÃ©s selon un parcours en profondeur modifiÃ©: tous
       les fils du sommet courant sont numÃ©rotÃ©s avant l'Ã©tape de
       rÃ©cursivitÃ©. Si n=1, alors le graphe est une Ã©toile Ã  b+z
       feuilles. Le dessin avec dot (-visu) ne respecte pas le
       plongement de l'arbre. L'option -directed permet d'obtenir une
       1-orientation.

....ringarytree h k r p
....
       Arbre de hauteur h oÃ¹ chaque noeud de profondeur < h Ã 
       exactement k fils, sauf la racine qui en possÃ¨de r. Lorsque
       p>0, un chemin (si p=1) ou un cycle (si p=2) est ajoutÃ© entre
       les sommets de mÃªme profondeur. Notez que "ringarytree h 1 r 0"
       gÃ©nÃ¨re une Ã©toile de degrÃ© r oÃ¹ chaque branche est de longueur
       h. NumÃ©rotÃ©s selon un parcours en profondeur, le nom des
       sommets est un mot correspondant au chemin depuis la racine.

....rectree h f_1 f_2 ... f_d .
....
       Arbre rÃ©cursif de hauteur h oÃ¹ chaque noeud profondeur < h Ã 
       exactement d fils, le i-Ã¨me fils Ã©tant la racine d'un arbre
       similaire de profondeur h-fáµ¢. Il faut fáµ¢>0. Ainsi "rectree h 1
       1 ." est un arbre binaire complet de hauteur h. Le graphe
       "rectree h 1 1 ... 1 ." est un arbre complet de hauteur h
       identique Ã  "ringarytree h d d 0". L'arbre est composÃ© d'un
       seul sommet si d=0 ou h<=0. C'est une Ã©toile si h>0 et d>=h. Le
       nombre de sommets est exponentiel en h dÃ¨s que dâ‰¥2, plus
       prÃ©cisÃ©ment de la forme ð›¼Â·Ï^h-1 pour des constantes ð›¼>0 et Ï>1
       dÃ©pendant des fáµ¢. Pour fâ‚=fâ‚‚=1, ð›¼=Ï=2. Pour fâ‚=1 et fâ‚‚=2,
       ð›¼â‰ƒ1.23 et Ïâ‰ƒ1.61. Pour fâ‚=1 et fâ‚‚=3, ð›¼â‰ƒ1.34 et Ïâ‰ƒ1.47.
       NumÃ©rotÃ©s selon un parcours en profondeur, le nom des sommets
       est un mot correspondant au chemin depuis la racine.

....kpage n k
....
       Graphe k-pages connexe alÃ©atoire. Un graphe k-page peut Ãªtre
       reprÃ©senter en plaÃ§ant les sommets le long d'un cercle, en
       dessinant les arÃªtes comme des segments de droites, et en
       coloriant les arÃªtes en k>0 couleurs de faÃ§on Ã  ce que les
       arÃªtes de chaque couleur induisent le dessin d'un graphe
       planaire-extÃ©rieur. La numÃ©rotation des sommets est faite le
       long du cercle. Les graphes 1-page sont les graphes
       planaires-extÃ©rieurs, les 2-pages sont les sous-graphes de
       graphes planaires Hamiltoniens. Les graphes planaires de degrÃ©
       au plus 4 sont 2-pages, les 3-arbres planaires (ou graphes
       Apolloniens) sont 3-pages, et les cliques avec 2k-1 ou 2k
       sommets des k-pages. L'option -directed permet d'obtenir une
       2k-orientation.
....
       Ces graphes sont construits par le processus alÃ©atoire suivant.
       On gÃ©nÃ¨re k graphes planaires-extÃ©rieurs alÃ©atoires uniformes
       connexes Ã  n sommets (plan et enracinÃ©) grÃ¢ce Ã  une bijection
       avec les arbres plans enracinÃ©s dont tous les sommets, sauf
       ceux de la derniÃ¨re branche, sont bicoloriÃ©s. On fait ensuite
       l'union de ces k graphes en choisissant alÃ©atoirement la racine
       des arbres, sauf celui du premier planaire-extÃ©rieur, ce qui
       correspond Ã  une permutation circulaire des sommets sur la face
       extÃ©rieure.

....cactus n
....
       Graphe cactus alÃ©atoire Ã  n sommets. Il s'agit d'arbres de
       cycles, c'est-Ã -dire de graphes connexes oÃ¹ chaque arÃªte
       appartient Ã  au plus un cycle. Ce sont aussi les graphes
       planaires-extÃ©rieurs connexes sans cordes. Ils sont gÃ©nÃ©rÃ©s Ã 
       partir d'un "outerplanar n" dans lequel les arÃªtes internes (ou
       cordes) des composantes biconnexes ont Ã©tÃ© supprimÃ©s. L'option
       -directed permet d'obtenir une 2-orientation.

....ktree n k
....
       k-arbre alÃ©atoire Ã  n sommets. Il faut n>kâ‰¥0. C'est un graphe
       chordal appelÃ© aussi graphe triangulÃ© (triangulated). Il est
       gÃ©nÃ©rÃ© Ã  partir d'un arbre enracinÃ© alÃ©atoire uniforme Ã  n-k
       noeuds de maniÃ¨re similaire Ã  "tree n-k". Cela constitue les
       "sacs" que l'on remplit avec les n sommets comme suit: on met
       k+1 sommets dans le sac racine connectÃ© en clique, puis, selon
       un parcours en profondeur de l'arbre, on met un sommet
       diffÃ©rent pour chacun des autres sacs. Ce sommet est alors
       connectÃ©s Ã  exactement k sommets choisis alÃ©atoirement dans le
       sac parent et sont ajoutÃ©s Ã  son sac. Lorsque k=1, c'est un
       arbre, et lorsque k=0, c'est un stable. L'option -directed
       permet d'obtenir une k-orientation.

....kpath n k
....
       k-chemin alÃ©atoire Ã  n sommets. La construction est similaire Ã 
       celle utilisÃ©e pour ktree, sauf que l'arbre est un chemin. Ces
       graphes sont des graphes d'intervalles particuliers (voir
       "interval n"). L'option -directed permet d'obtenir une
       k-orientation.

....kstar n k
....
       k-star alÃ©atoire Ã  n sommets. La construction est similaire Ã 
       celle utilisÃ©e pour ktree, sauf que l'arbre est une Ã©toile. Ces
       graphes, qui sont des "split graphs", sont composÃ©s d'une
       clique Ã  k+1 sommets et de n-k-1 sommets indÃ©pendants connectÃ©s
       Ã  k sommets alÃ©atoire de la clique. Il est possible d'obtenir
       le graphe "split n k" si Ã  chaque fois les k sommets de la
       clique tirÃ©s alÃ©atoires sont toujours les mÃªmes. L'option
       -directed permet d'obtenir une k-orientation.

....rig n k p
....
       Graphe d'intersections alÃ©atoire (Uniform Random Intersection
       Graph). Il possÃ¨de n sommets, chaque sommet u Ã©tant reprÃ©sentÃ©
       par un sous-ensemble S(u) alÃ©atoire de {1,â€¦,k} tel que chaque
       Ã©lÃ©ment appartient Ã  S(u) avec probabilitÃ© p. Il y a une arÃªte
       entre u et v ssi S(u) et S(v) s'intersectent. La probabilitÃ©
       d'avoir une arÃªte entre u et v est donc Pâ‚‘=1-(1-pÂ²)^k, mais les
       arÃªtes ne sont pas indÃ©pendantes (Pr(uv|uw)>Pr(uv)). En
       gÃ©nÃ©ral, pour ne pas avoir Pâ‚‘ qui tend vers 1, on choisit les
       paramÃ¨tres de faÃ§on Ã  ce que kpÂ²<cste.  Lorsque kâ‰¥nÂ³, ce modÃ¨le
       est Ã©quivalent au modÃ¨le des graphes alÃ©atoires d'ErdÃ¶s-Reny
       (voir random n p). Si p<0, alors p est fixÃ©e au seuil thÃ©orique
       de connectivitÃ©, Ã  savoir p=âˆš(ln(n)/(nk)) si k>n et p=ln(n)/k
       sinon.

....apollonian n
....
       Graphe Apollonien alÃ©atoire uniforme Ã  nâ‰¥4 sommets. Les graphes
       Apolloniens sont les 3-arbres planaires ou encore les graphes
       planaires maximaux chordaux. Ils sont obtenus en subdivisant
       rÃ©cursivement un triangle en trois autres. Ils sont
       3-dÃ©gÃ©nÃ©rÃ©s, de treewidth 3, et de nombre chromatique 4. La
       distance moyenne est Ï´(logn). Ils sont en bijection avec les
       arbres ternaires Ã  n-3 noeuds internes. Pour n=5, il s'agit
       d'un Kâ‚… moins une arÃªte qu'on peut obtenir aussi avec "split 5
       3". L'option -directed permet d'obtenir une 3-orientation.

....polygon n
....
       Triangulation alÃ©atoire uniforme d'un polygone convexe Ã  nâ‰¥3
       cotÃ©s. Ce sont aussi des graphes planaires-extÃ©rieurs maximaux
       alÃ©atoires. Ils sont Hamiltoniens, 2-dÃ©gÃ©nÃ©rÃ©s, de treewidth 2,
       et de nombre chromatique 3. Ils sont en bijection avec les
       arbres binaires (complets) Ã  n-2 noeuds internes ou encore les
       mots de Dyck de longueur 2(n-2). La numÃ©rotation des sommets
       n'est pas cyclique le long du polygone. Ce graphe n'est pas un
       graphe gÃ©omÃ©trique contrairement Ã  ses variantes utilisant -xy
       convex comme dans l'exemple ci-aprÃ¨s.
....
       Ex: gengraph polygon 20 -dot filter circo -visu
           gengraph td-delaunay 20 -xy convex2 -visu

....planar n f d
....
       Graphe planaire alÃ©atoire composÃ© de n faces internes de
       longueur fâ‰¥3, les sommets internes Ã©tant de degrÃ© au moins d et
       ceux de la face externe au moins 2. Ils possÃ¨dent entre n+f-1
       et n(f-2)+2 sommets, sont 2-connexes, 2-dÃ©gÃ©nÃ©rÃ©s, de maille
       f. Si d>4 alors ils sont d'hyperbolicitÃ© O(f). Ils sont
       construits en ajoutant itÃ©rativement les faces par le processus
       alÃ©atoire suivant. Au dÃ©part, il s'agit d'un cycle de longueur
       f. Pour chaque nouvelle face on ajoute un sommet u que l'on
       connecte Ã  un sommet quelconque du cycle C formant le bord de
       la face extÃ©rieure du graphe courant. Puis on ajoute un chemin
       allant de u Ã  un sommet v de C de faÃ§on Ã  respecter: 1) la
       contrainte des degrÃ©s des sommets qui vont devenir internes; et
       2) la contrainte sur la longueur de la nouvelle face crÃ©Ã©e. Le
       sommet v est choisit uniformÃ©ment parmi tous les sommets
       possibles de C respectant les deux contraintes. Si d<0, alors
       on fait comme si d=+âˆž (aucun sommet ne pouvant alors Ãªtre
       interne) et le rÃ©sultat est un graphe planaire-extÃ©rieur
       Hamiltonien, c'est-Ã -dire 2-connexe. Si f<0, alors chaque face
       crÃ©Ã©e est de longueur alÃ©atoire uniforme prise dans [3,|f|] au
       lieu d'Ãªtre de longueur exactement |f|. Si f=d=4, il s'agit
       d'un "squaregraph". Les valeurs d=0,1,2 sont Ã©quivalentes.
       L'option -directed permet d'obtenir une 2-orientation.

....hyperbolic p k h
....
       Graphe issu du pavage du plan hyperbolique ou euclidien par des
       polygones rÃ©guliers Ã  pâ‰¥3 cotÃ©s oÃ¹ chaque sommets est de degrÃ©
       kâ‰¥2. Le graphe est construit par couches successives de
       polygones, le paramÃ¨tre hâ‰¥1 reprÃ©sentant le nombre de couches.
       Lorque h=1, il s'agit d'un seul polygone, un cycle Ã  p sommets.
       Dans tous les cas les graphes sont planaires avec O((pk)^h)
       sommets, ils sont 2-connexes et d'arboricitÃ© 2 pour p>3.
       Lorsque p=3, ils sont 3-connexes et d'arboricitÃ© 3. L'option
       -directed permet d'obtenir une 3-orientation. Sans Ãªtre les
       mÃªmes graphes, il y a des similaritÃ©s avec les graphes "planar
       n f d". Pour paver le plan hyperbolique, reprÃ©sentable sur le
       disque de PoincarÃ©, il faut 1/p + 1/k < 1/2. Dans ce cas, le
       graphe est d'hyperbolictÃ© O(p). Pour paver le plan euclidien il
       faut prendre p=k=4 (grille carrÃ©e) ou p=3, k=6 (grille
       triangulaire) ou p=6, k=3 (grille hexagonale). Si kâ‰¤3 et pâ‰¤5,
       alors le graphe existe que pour certaines valeurs de hâ‰¤3.  Le
       cas k=3, p=4, h=2 correspond au cube, et k=3, p=5, h=3 est le
       graphe dodÃ©cahÃ¨dre ("dodecahedron"). Si k=2, alors h=1 et le
       graphe est un cycle Ã  p sommets.

....rlt p q d
....
       Random Lattice Triangulation. Il s'agit d'un graphe planaire
       alÃ©atoire construit Ã  partir d'une triangulation de la grille p
       Ã— q. Toutes les faces, sauf celle extÃ©rieure, sont des
       triangles. Il possÃ¨de pq sommets et (2p-1)(2q-1)-pq arÃªtes,
       dont les 2(p+q-2) qui sont sur le bord de la grille et le reste
       est Ã  l'intÃ©rieur. Le paramÃ¨tre d contrÃ´le la longueur des
       arÃªtes suivant la norme Lmax. Par exemple, si d=1, les arÃªtes
       seront soit celles de la grille (verticales ou horizontales) ou
       diagonales. Si d<0, alors l'effet est similaire Ã  d=+âˆž. Si d=0,
       on obtient un stable. Si p=1 ou q=1 (et dâ‰ 0), on obtient un
       chemin.
....
       Ex: gengraph rlt 8 14 2 -dot scale 0.1 -visu
....
       Il est difficile de gÃ©nÃ©rer de telles grilles alÃ©atoirement
       uniformÃ©ment. Il faut pour cela utiliser une technique de flips
       avec une chaÃ®ne de Markov dont le mixing time n'est pas
       connu. Il est cependant bien connu que le milieu de chaque
       arÃªte e d'une triangulation T de grille a pour coordonnÃ©es
       (i+1/2,j) ou (i+1/2,j+1/2) oÃ¹ i,j sont des entiers. Et
       inversement, chaque point "milieu" de la grille est coupÃ© par
       exactement une arÃªte de T. Il est aussi connu que si l'on
       parcoure les points milieux de la grille de bas en haut et de
       gauche Ã  droit, alors il n'y qu'au plus deux choix possibles
       pour l'arÃªte ayant ce milieu. Malheureusement, suivre ce
       parcours et choisir alÃ©atoirement l'une ou l'autre de ces
       arÃªtes ne donne pas une distribution alÃ©atoire uniforme. On
       propose ici la construction d'une triangulation T de maniÃ¨re
       alÃ©atoire comme suit:
....
       Tant qu'il reste au moins un point milieu faire:
       1. choisir uniformÃ©ment un point milieu R parmi les points restant
       2. dÃ©terminer la liste L des arÃªtes possibles ayant pour milieu R
          (respectant la planaritÃ© et le critÃ¨re de longueur)
       3. choisir uniformÃ©ment une arÃªte de L et l'ajouter au graphe

....kneser n k r
....
       Graphe de Kneser gÃ©nÃ©ralisÃ©. Le graphe de Kneser K(n,k)
       classique est obtenu avec r=0. Les sommets sont tous les
       sous-ensembles Ã  k Ã©lÃ©ments de [0,n[ (il faut donc 0â‰¤kâ‰¤n).
       Deux sommets sont adjacents ssi leurs ensembles correspondant
       ont au plus r Ã©lÃ©ments en commun. Le nombre chromatique de
       K(n,k), Ã©tablit par LovÃ¡sz, vaut n-2k+2 pour tout nâ‰¥2k-1>0. Le
       graphe de Petersen est le graphe K(5,2). Ils ont un lien avec
       les graphes de Johnson J(n,k).

....gpetersen n r
....
       Graphe de Petersen gÃ©nÃ©ralisÃ© P(n,r), 0â‰¤r<n/2. Ce graphe
       cubique possÃ¨de 2n sommets qui sont u_1,â€¦,u_n,v_1,â€¦,v_n.  Les
       arÃªtes sont, pour tout i: u_i-u_{i+1}, u_i-v_i et v_i-v_{i+r}
       (indice modulo n). Il peut Ãªtre dessinÃ© tel que toute ses
       arÃªtes sont de mÃªme longueur (unit distance graph). Ce graphe
       est biparti ssi n est pair et r est impair. C'est un graphe de
       Cayley ssi rÂ²=1 (modulo n). P(n,r) est Hamiltonien ssi râ‰ 2 ou
       nâ‰ 5 (modulo 6). P(n,r) est isomorphe Ã  P(n,(n-2r+3)/2)).
       P(4,1) est le cube, P(5,2) le graphe de Petersen, P(6,2) le
       graphe de DÃ¼rer, P(8,3) le graphe de MÃ¶bius-Kantor, P(10,2) le
       dodÃ©caÃ¨dre, P(10,3) le graphe de Desargues, P(12,5) le graphe
       de Nauru, P(n,1) un prisme.

....squashed n k p
....
       Squashed Cube alÃ©atoire Ã  n sommets. Il faut 0<k<n et pâˆˆ[0,1].
       Les sommets sont des mots alÃ©atoires de k lettres sur {0,1,'*'}
       oÃ¹ p est la probabilitÃ© d'obtenir '*'. La probabilitÃ© d'obtenir
       0 est la mÃªme que celle d'obtenir 1, soit (1-p)/2. Deux sommets
       sont adjacents si la distance de Hamming entre leur mots vaut 1
       avec la convention que la distance Ã  la lettre '*' est nulle.
       Lorsque que p=0, le graphe gÃ©nÃ©rÃ© est un sous-graphe
       isomÃ©trique de l'hypercube oÃ¹ certains sommets sont dupliquÃ©s
       en sommets jumeaux non-adjacents (ce sont les sommets
       correspondant au mÃªme mot). En particulier, le graphe est
       biparti et la distance entre deux sommets est donnÃ©es par la
       distance de Hamming entre leur mot, sauf s'ils ont le mÃªme
       mot. Si p=-1 alors p est fixÃ©e Ã  l'Ã©quiprobabilitÃ© de chacune
       des lettres, soit p=1/3.

....antiprism n
....
       Graphe composÃ© de deux cycles Ã  n sommets connectÃ©s par 2n
       triangles. Le prisme est similaire sauf que pour ce dernier les
       deux cycles sont connectÃ©s par des carrÃ©s. Il est ainsi
       planaire, 4-rÃ©gulier, possÃ¨de 2n sommets et 4n arÃªtes. C'est
       aussi le dual du trapÃ©zoÃ¨dre n-gonal.

....rpartite a_1 â€¦ a_k .
....
       Graphe k-parti complet K_{a_1,â€¦,a_k}. Ce graphe possÃ¨de a_1 + â‹¯
       + a_k sommets partitionnÃ©s en k parts. La i-Ã¨me part contient
       a_i sommets numÃ©rotÃ©s consÃ©cutivement dans l'intervalle [ a_1 +
       â‹¯ + a_{i-1}, a_1 + â‹¯ + a_i [. Les sommets i et j sont adjacents
       ssi i et j appartiennent Ã  des parts diffÃ©rentes.

....ggosset p d_1 v_1 â€¦ d_k v_k .
....
       Graphe de Gosset gÃ©nÃ©ralisÃ©. Les sommets sont tous les vecteurs
       entiers (et leurs opposÃ©s) de dimension d = d_1 + â‹¯ + d_k dont
       les coordonnÃ©es comprennent exactement dáµ¢â‰¥1 fois la valeur
       váµ¢. Le nombre de sommets est donc n = 2Â·âˆ_i
       binomial(d-(âˆ‘_{j<i}d_i), d_i). Il existe une arÃªte entre les
       vecteurs u et v si et seulement le produit scalaire entre u et
       v vaut l'entier p. Des valeurs intÃ©ressantes sont par exemple
       "ggosset 1 2 -1 2 0 ."  ou "ggosset 8 2 3 6 -1 ."  (le graphe
       de Gosset).

....schlafli
....
       Graphe de SchlÃ¤fli. Il s'agit du sous-graphe induit par les
       voisins d'un quelconque sommet du graphe de Gosset. Il possÃ¨de
       27 sommets, 216 arÃªtes et est 16-rÃ©gulier. Il est sans-griffe,
       Hamiltonien, de diamÃ¨tre 2, de maille 3 et de nombre
       chromatique 9.

....crown n
....
       Graphe biparti rÃ©gulier Ã  2n sommets oÃ¹ le i-Ã¨me sommet de la
       premiÃ¨re partie de taille n est voisin au j-Ã¨me sommet de la
       seconde partie ssi iâ‰ j. Pour n=3, il s'agit du cycle Ã  6
       sommets, pour n=4, il s'agit du cube (Ã  8 sommets).

....split n k
....
       Graphe split (ou fendu) Ã  n sommets et de clique maximum k. Il
       s'agit d'un graphe Ã  n sommets composÃ© d'une clique Ã  k sommets
       et d'un ensemble indÃ©pendant de n-k sommets connectÃ©s chacun Ã 
       tous ceux de la clique. C'est un graphe triangulÃ© (chordal) et
       un cas particulier de "kstar n k". On peut montrer que presque
       tous les graphes triangulÃ©s Ã  n sommets sont des graphes split
       (Bender et al. 1985). Si kâ‰¥n-1, alors il s'agit d'une clique,
       et si k=n-2, il s'agit d'une clique moins une arÃªte. Si k=1 il
       s'agit d'une Ã©toile Ã  n-1 branches.

....fan p q
....
       Graphe de p+q sommets composÃ© d'un chemin Ã  p sommets et de q
       sommets, chacun connectÃ©s Ã  tous ceux du chemin. Le graphe
       classique "fan n" correspond Ã  p=n et q=1.

....flip n
....
       Graphe des flips des triangulations d'un polygone convexe Ã  n>2
       sommets. Les sommets, qui sont les triangulations, sont en
       bijection avec des arbres binaires (complets) Ã  m=n-2 noeuds
       internes qui sont codÃ©s par les mots de Dyck de longueur 2m
       (mots que l'on peut afficher avec -label 1). Le nombre de
       sommets est donc C(m) = binom(2m,m)/(m+1), le nombre de Catalan
       d'ordre m. Les adjacences peuvent Ãªtre vues aussi comme des
       rotations d'arbres. Le diamÃ¨tre est 2n-10 pour n>12. Le nombre
       chromatique n'est pas connu. On sait pas s'il est constant ou
       pas. Il vaut 3 pour n=5..9, et 4 pour n=10 et 11.

....interval n
....
       Graphe d'intersection de n intervalles d'entiers alÃ©atoires
       uniformes pris dans [0,2n[. Des graphes d'intervalles peuvent
       aussi Ãªtre gÃ©nÃ©rÃ©s par "kpath n k".

....circle n
....
       Graphe d'intersection de n cordes alÃ©atoires d'un cercle. Il
       est rÃ©alisÃ© par le graphe d'inclusion de n intervalles
       d'entiers alÃ©atoires uniformes pris dans [0,2n[. Les graphes de
       permutation et les planaires extÃ©rieurs sont des exemples de
       circle graphs.

....permutation n
....
       Graphe de permutation sur une permutation alÃ©atoire uniforme
       des entiers de [0,n[.

....prime n
....
       Graphe Ã  n sommets tel que i est adjacent Ã  j ssi i>1 et j
       divisible par i.

....paley n
....
       Graphe de Paley Ã  n sommets. Deux sommets sont adjacents ssi
       leur diffÃ©rence est un carrÃ© modulo n. Il faut que n soit la
       puissance d'un nombre premier et que nâ‰¡1 (mod 4), mais le
       graphe est aussi dÃ©fini pour les autres valeurs. Les premiÃ¨res
       valeurs possibles pour n sont: 5, 9, 13, 17, 25, 29, 37, 41,
       49, â€¦ Ces graphes sont Hamiltoniens. Si n est simplement
       premier, alors ils sont de plus auto-complÃ©mentaires et
       rÃ©guliers.  Paley 17 est le plus grand graphe G oÃ¹ ni G ni son
       complÃ©mentaire ne contient Kâ‚„, d'oÃ¹ Ramsey(4)=18.

....mycielski k
....
       Graphe de Mycielski de paramÃ¨tre (nombre chromatique) k. C'est
       un graphe sans triangle, k-1 (sommets) connexe, et de nombre
       chromatique k. Le premier graphe de la sÃ©rie est M2 = K2, puis
       on trouve M3=C5, M4 est le graphe de GrÃ¶tzsch Ã  11 sommmets.

....windmill n
....
       Graphe composÃ© de n cycles de longueur trois ayant un sommet
       commun.

....barbell n1 n2 p
....
       Graphe des haltÃ¨res (Barbell Graph) composÃ© de deux cliques de
       n1 et n2 sommets reliÃ©es par un chemin de longueur p. Il
       possÃ¨de n1+n2+p-1 sommets. Si p=0 (p=-1), le graphe est composÃ©
       de deux cliques ayant un sommet (une arÃªte) en commun. Plus
       gÃ©nÃ©ralement, si pâ‰¤0, le graphe est composÃ© de deux cliques
       s'intersectant sur 1-p sommets.

....chess p q x y
....
       Graphe composÃ© de p x q sommets reprÃ©sentant les cases d'un
       Ã©chiquier p x q, deux cases Ã©tant connectÃ©e s'il existe un
       dÃ©placement d'une case vers l'autres avec un saut de x cases
       selon un axe et y selon un autre. Le "knight graph" classique
       est donc un "chess 8 8 2 3", et "chess n 2 1 0" correspond Ã 
       "ladder n".

....sat n m k
....
       Graphe alÃ©atoire issu de la rÃ©duction du problÃ¨me k-SAT Ã 
       Vertex Cover. Le calcul d'un Vertex Cover de taille minimum
       pour ce graphe est donc difficile pour k>2. Soit F une formule
       de k-SAT avec n>0 variables x_i et m>0 clauses CNF de k>0
       termes.  Le graphe gÃ©nÃ©rÃ© par "sat n m k" possÃ¨de un Vertex
       Cover de taille n+(k-1)m si et seulement si F est satisfiable.
....
       Ce graphe est composÃ© d'une union de n arÃªtes indÃ©pendantes et
       de m cliques Ã  k sommets, plus des arÃªtes dÃ©pendant de F
       connectant certains sommets des cliques aux n arÃªtes. Les n
       arÃªtes reprÃ©sentent les n variables, une extrÃ©mitÃ© pour x_i,
       l'autre pour Â¬(x_i). Ces sommets ont des numÃ©ros dans [0,2n[,
       x_i correspond au sommet 2i-2 et Â¬(x_i) au sommet 2i-1,
       i=1â€¦n. Les sommets des cliques ont des numÃ©ros consÃ©cutifs â‰¥ 2n
       et correspondent aux clauses. Le p-Ã¨me sommet de la q-Ã¨me
       clique (pour p=1â€¦k et q=1â€¦m) est connectÃ© Ã  l'une des
       extrÃ©mitÃ©s de la i-Ã¨me arÃªte (pour i=1â€¦n) ssi la p-Ã¨me variable
       de la q-Ã¨me clause est x_i ou Â¬(x_i).
....
       La formule F est construite en choisissant indÃ©pendamment et
       uniformÃ©ment pour chacune des m clauses et chacun des k termes
       une des variables parmi x_1,â€¦,x_n,Â¬(x_1),â€¦,Â¬(x_n).  Ainsi
       chaque sommet d'une clique possÃ¨de exactement un voisin (choisi
       alÃ©atoirement uniforme) parmi les 2n extrÃ©mitÃ©s d'arÃªtes.

....kout n k
....
       Graphe Ã  n sommets k-dÃ©gÃ©nÃ©rÃ© crÃ©e par le processus alÃ©atoire
       suivant: les sommets sont ajoutÃ©s dans l'ordre croissant de
       leur numÃ©ro, i=0,1,â€¦,n-1. Le sommet i est connectÃ© Ã  d voisins
       qui sont pris alÃ©atoirement uniformÃ©ment parmi les sommets dont
       le numÃ©ro est < i. La valeur d est choisie alÃ©atoirement
       uniformÃ©ment entre 1 et min{i,k}. Il faut 0<k<n. Le graphe est
       connexe, et pour k=1, il s'agit d'un arbre. L'option -directed
       permet d'obtenir une k-orientation.
....
       Ex: gengraph kout 50 2 -directed -vcolor deg -vcolor pal wbn
                    -vsize -visu

....expander n k
....
       Graphe Ã  n sommets composÃ© de k>0 cycles Hamiltoniens
       alÃ©atoires. Le degrÃ© des sommets varie entre 2 et 2k. Il
       possÃ¨de le cycle 0,1,â€¦,n-1,0 comme cycle Hamiltonien, et a la
       propriÃ©tÃ© d'expansion Ã  partir de kâ‰¥4. Plus prÃ©cisÃ©ment, avec
       grande probabilitÃ©, les valeurs propres de la matrice
       d'adjacence du graphe sont â‰¤ 2âˆš(2k). On rappelle que la
       constante d'expansion, isopÃ©rimÃ©trique, ou de Cheeger h(G) d'un
       graphe d-rÃ©gulier G est toujours comprise entre (d-Î»â‚‚)/2 â‰¤ h(G)
       â‰¤ âˆš(2d(d-Î»â‚‚)) oÃ¹ Î»â‚‚ est la deuxiÃ¨me plus grande valeur propre
       de la matrice d'adjacence de G. L'option -directed permet
       d'obtenir une k-orientation.

....margulis n
....
       Graphe de Margulis Ã  n^2 sommets. Il s'agit d'un expandeur avec
       Î»â‚‚ â‰¤ 5âˆš2 (cf. graphe "expander n k") de degrÃ© maximum 8 et de
       degrÃ© minimum 2. Les sommets sont les paires d'entiers (x,y)
       avec x,y âˆˆ [0,n[ avec les 8 adjacences suivantes: (x+y,y),
       (x-y,y), (x,y+x), (x,y-x), (x+y+1,y), (xâˆ’y+1,y), (x,y+x+1) et
       (x,yâˆ’x+1), toutes ces opÃ©rations Ã©tant modulo n.

....comb n
....centipede n
....
       Arbre de 2n sommets en forme de peigne, composÃ© d'un chemin Ã  n
       sommets avec un sommet pendant Ã  chacun d'eux. On peut
       l'obtenir en supprimant une arÃªte d'un sunlet n.

....sunlet n
....
       Cycle Ã  n sommets avec un sommet pendant Ã  chacun d'eux. Un
       sunlet 3 est parfois appelÃ© netgraph.

....parachute n
....
       Graphe Parachute. Il est planaire Ã  n+3 sommets composÃ©s du
       graphe "fan n 2" dont un des deux sommets de degrÃ© n possÃ¨de un
       sommet pendant. Le graphe classique correspond Ã  n=4. C'est le
       complÃ©menataire du graphe parapluie n.

....alkane t n
....
       Graphe planaires dont les sommets sont de degrÃ© 1 ou 4
       reprÃ©sentant la structure molÃ©culaire d'hydrocarbure alkalin Ã 
       n atomes de carbones. Le paramÃ¨tre t (voir ci-dessous les six
       types possibles) contrÃ´le la topologie des liaisons simples
       entre atomes de carbone (C), les atomes d'hydrogÃ¨nes (H) Ã©tant
       des sommets pendants de sorte que chaque atome C soit de degrÃ©
       4. Certaines topologies ne sont dÃ©finies que pour certaines
       valeurs de n. Les alkalins, de formule C_n H_{2n+2} si
       typeâ‰ "cyclo", sont des arbres de 3n+2 sommets alors que les
       cycloalkalin, de formule C_n H_{2n} si type="cyclo", ont 3n
       sommets et possÃ¨de un cycle. Chaque type peut abrÃ©gÃ© par ses 2
       premiÃ¨res lettres. L'option "-label 1" activÃ©e par dÃ©faut
       permet de distinguer les atomes C et H.
....
       !!! t      topologie     n             t     topologie     n
                                                    Câ”€â”
          normal  Câ”€ â‹¯ â”€C      nâ‰¥1           neo    Câ”€Câ”€ â‹¯ â”€C    nâ‰¥5
                                                    Câ”€â”˜
                  â”Œâ”€Câ”€ â‹¯ â”€C                           â”Œâ”€Câ”€ â‹¯ â”€C
          cyclo   C       â”‚    nâ‰¥3           sec    Câ”€C          nâ‰¥6
                  â””â”€Câ”€ â‹¯ â”€C                           â””â”€Câ”€ â‹¯ â”€C
                  Câ”€â”                                 â”Œâ”€Câ”€ â‹¯ â”€C
          iso       Câ”€ â‹¯ â”€C    nâ‰¥4           tret   Câ”€Câ”€Câ”€ â‹¯ â”€C  nâ‰¥7
                  Câ”€â”˜                                 â””â”€Câ”€ â‹¯ â”€C
....
       Il est possible d'utiliser les alias suivants:
....
       !!!  n-alkane n ......... (= alkane normal n) 
            cyclo-alkane n ..... (= alkane cyclo n)
	    iso-alkane n ....... (= alkane iso n)
	    neo-alkane n ....... (= alkane neo n)
	    sec-alkane n ....... (= alkane sec n)
	    tret-alkane n ...... (= alkane tret n)
	    methane ............ (= alkane normal 1)
	    ethane ............. (= alkane normal 2)
	    propane ............ (= alkane normal 3)
	    butane ............. (= alkane normal 4)
	    pentane ............ (= alkane normal 5)
	    hexane ............. (= alkane normal 6)
	    heptane ............ (= alkane normal 7)
	    octane ............. (= alkane normal 8)
	    nonane ............. (= alkane normal 9)
....
      Il est possible aussi de combiner les prÃ©fixes cyclo-, iso-,
      neo-, sec-, tret- avec les radicaux meth, eth, prop, but, hex,
      hept, oct, non, lorsque la condition sur n est satisfaite. Par
      exemple, cyclo-pentane (= alkane cyclo 5) et iso-butane (=
      alkane iso 4).

....icosahedron
....
       IsocahÃ¨dre: graphe planaire 5-rÃ©gulier Ã  12 sommets. Il possÃ¨de
       30 arÃªtes et 20 faces qui sont des triangles. C'est le dual du
       dodÃ©cahÃ¨dre.

....rdodecahedron
....
       Rhombic-dodÃ©caÃ¨dre: graphe planaire Ã  14 sommets avec des
       sommets de degrÃ© 3 ou 4. Il possÃ¨de 21 arÃªtes et 12 faces qui
       sont des carrÃ©s. C'est le dual du cuboctaÃ¨dre.

....deltohedron n
....trapezohedron n
....
       DeltoÃ¨dre ou trapÃ©zoÃ¨dre n-gonal: graphe composÃ© de 2n faces en
       forme de cerf-volant (deltoÃ¯des) dÃ©calÃ©es symÃ©triquement. C'est
       un graphe planaire de 2n+2 sommets et 4n arÃªtes oÃ¹ toutes les
       faces sont des carrÃ©es. C'est aussi le dual de l'antiprisme
       n-gonal. Il s'agit d'un cube si n=3.

....tutte
....
       Graphe de Tutte. C'est un graphe planaire cubique 3-connexe Ã 
       46 sommets qui n'est pas Hamiltonien.

....hgraph
....
       Arbre Ã  six sommets dont quatre feuilles en forme de H.

....rgraph
....fish
....
       Fish Graph, graphe Ã  six sommets en forme de R ou de
       poisson. Il est composÃ© d'un cycle Ã  quatre sommets dont un
       ayant deux sommets pendants.

....cricket
....
       Cricket Graph, graphe Ã  cinq sommets composÃ© d'un triangle oÃ¹ Ã 
       l'un des sommets est attachÃ© deux sommets pendant (de degrÃ© 1).

....moth
....
       Moth Graph, graphe Ã  six sommets composÃ© de deux triangles
       partageant une arÃªte et de deux sommets pendant (degrÃ© 1)
       attachÃ©s Ã  un sommet de degrÃ© trois.

....dart
....
       Dart Graph, graphe Ã  cinq sommets composÃ© de deux triangles
       partageant une arÃªte et d'un sommet pendant (degrÃ© 1) attachÃ© Ã 
       un sommet de degrÃ© trois. Il peut Ãªtre obtenu Ã  partir du
       graphe moth en supprimant un sommet pendant.

....bull
....
       Bull Graph, graphe Ã  cinq sommets auto-complÃ©mentaire en forme
       de A.

....antenna
....
       Antenna Graph, graphe planaire Ã  six sommets formÃ© du graphe
       house et d'un sommet pendant attachÃ© Ã  son toit. Plus
       prÃ©cisÃ©ment, il est composÃ© d'un carrÃ©, d'un triangle
       partageant une arrÃªte et d'un sommet pendant au sommet de degrÃ©
       2 du triangle. C'est le complÃ©mentaire du graphe formÃ© d'un
       carrÃ© et de deux triangles partageant deux arÃªtes consÃ©cutive
       du carrÃ©.

....suzuki
....
       Graphe de Suzuki (2010). C'est l'unique graphe 1-planaire Ã 
       n=11 sommets et ayant le nombre optimal d'arÃªtes, soit 4n-8
       arÃªtes (ici 36 donc).

....harborth
....
       Graphe de Harborth. C'est un graphe planaire 4-rÃ©gulier Ã  52
       sommets qui est distance unitaire aussi appelÃ© graphe allumette
       (voir theta0 et diamond). Il peut ainsi Ãªtre dessinÃ© sans
       croisement d'arÃªte qui ont toutes la mÃªme longueur.

....doily
....
       Graphe Doily (de Payne). C'est un graphe de 15 sommets qui est
       un carrÃ© gÃ©nÃ©ralisÃ© pouvant Ãªtre reprÃ©sentÃ© par 15 points et 15
       lignes, avec 3 points par ligne et 3 lignes par point, et sans
       triangle.

....herschel
....
       Graphe de Herschel. C'est le plus petit graphe planaire
       3-connexe qui ne soit pas Hamiltonien. Il est biparti, possÃ¨de
       11 sommets et 18 arÃªtes.

....goldner-harary
....
       Graphe de Goldner-Haray. C'est le plus petit graphe planaire
       maximal qui ne soit pas Hamiltonien. Il possÃ¨de 11 sommets et
       donc 27 arÃªtes (voir aussi Herchel). C'est un 3-arbre planaire
       (voir apollonian).

....fritsch
....
       Graphe de Fritsch. Il est planaire maximal Ã  9 sommets qui peut
       Ãªtre vu comme un graphe de HajÃ³s dans un triangle. C'est, avec
       le graphe de Soifer, le plus petit contre-exemple Ã  la
       procÃ©dure de coloration de Kempe. C'est le plus petit graphe oÃ¹
       l'heuristique de degrÃ© minimum donne cinq couleurs.

....triplex
....
       Graphe cubique de maille 5 Ã  12 sommets 1-planaire pouvant Ãªtre
       dessinÃ© avec seulement deux croisements d'arÃªte. Un des cinq
       graphes (avec le Petersen) a Ãªtre cycliquement-5-connexe
       (McCuaig 1992).

....jaws
....
       Graphe cubique de maille 5 Ã  20 sommets qui est un doublecross,
       c'est-Ã -dire dessinable sur le plan avec deux paires d'arÃªtes
       se croisant sur la face extÃ©rieure. Il est donc 1-planaire.
       Tout graphe theta-connectÃ© sans Petersen mais avec Jaws comme
       mineur est un doublecross.

....starfish
....
       Graphe cubique de maille 5 Ã  20 sommets non-planaire, mais
       peut-Ãªtre dessinÃ© comme une Ã©toile Ã  cinq branches avec une
       couronne centrale Ã  15 sommets formant un circulant avec une
       corde de longueur 3. Un graphe theta-connectÃ© (cf. Seymour et
       al. 2015) ssi il ne contient pas de Petersen comme mineur, si
       c'est un graphe apex (planaire plus un sommet), un doublecross
       (voir jaws) ou un starfish.

....soifer
....
       Graphe de Soifer. Il est planaire maximal Ã  9 sommets. C'est,
       avec le graphe de Fritsch, le plus petit contre-exemple Ã  la
       procÃ©dure de coloration de Kempe. C'est le plus petit graphe oÃ¹
       l'heuristique de degrÃ© minimum donne cinq couleurs.

....poussin
....
       Graphe de Poussin. Il est planaire maximal Ã  15 sommets. C'est
       un contre-exemple Ã  la procÃ©dure de coloration de Kempe.

....heawood4
....
       Graphe de Heawood pour la conjecture des 4 couleurs,
       contre-exemple de la preuve de Kempe. Il est planaire maximal
       avec 25 sommets, est de nombre chromatique 4, de diamÃ¨tre 5, de
       rayon 3 et Hamiltonien.

....errera
....
       Graphe d'Errera. Il est planaire maximal Ã  17 sommets. C'est un
       contre-exemple Ã  la procÃ©dure de coloration de Kempe.

....kittell
....
       Graphe de Kittell. Il est planaire maximal Ã  23 sommets. C'est
       un contre-exemple Ã  la procÃ©dure de coloration de Kempe.

....frucht
....
       Graphe de Frucht. Il est planaire cubique Ã  12 sommets. Il n'a
       pas de symÃ©trie non-triviale. C'est un graphe de Halin de
       nombre chromatique 3, de diamÃ¨tre 4 et de rayon 3.

....treep p
....
       Arbre alÃ©atoire Ã  p>2 feuilles sans sommets internes de degrÃ©
       deux. Il possÃ¨de entre p+1 et 2p-2 sommets. Ce graphe est Ã  la
       base de la construction des graphes de Halin.

....halin p
....
       Graphe de Halin alÃ©atoire basÃ© sur un arbre Ã  p>2 feuilles. Il
       possÃ¨de entre p+1 et 2p-2 sommets. Il est constituÃ© d'un arbre
       sans sommets de degrÃ© deux dont les p feuilles sont connectÃ©s
       par un cycle (de p arÃªtes). Ces graphes planaires de degrÃ©
       minimum au moins trois sont aussi arÃªte-minimale 3-connexes,
       Hamiltonien (et le reste aprÃ¨s la suppression de n'importe quel
       sommet), de treewidth exactement 3 (ils contiennent Kâ‚„ comme
       mineur). Ils contiennent toujours au moins trois triangles et
       sont de nombre chromatique 3 ou 4.

....butterfly d
....
       Graphe Butterfly de dimension d. Les sommets sont les paires
       (x,i) oÃ¹ x est un mot binaire de d bits et i un entier de
       [0,d]. Les sommets peuvent Ãªtre reprÃ©sentÃ©s en d+1 niveaux
       chacun de 2^d sommets, les arÃªtes connectant les niveaux
       consÃ©cutifs. Le sommet (x,i) est adjacent Ã  (y,i+1) ssi les
       bits de x sont identiques Ã  ceux de y sauf pour celui de numÃ©ro
       i+1 (le bit 1 Ã©tant le bit de poids le plus faible). Il possÃ¨de
       (d+1)Â·2^d sommets et dÂ·2^(d+1) arÃªtes, les sommets de niveau 0
       et d Ã©tant de degrÃ© 2 les autres de degrÃ© 4.

....shuffle d
....
       Graphe Shuffle-Exchange de dimension d. Les sommets sont les
       mots binaires de d lettres. Le sommet w et w' sont voisins si w
       et w' diffÃ¨rent du dernier bit, ou bien si w' peut Ãªtre obtenu
       par dÃ©calage cyclique Ã  droite ou Ã  gauche de w.

....debruijn d b
....
       Graphe de De Bruijn de dimension dâ‰¥0 et de base b>0. Il a b^d
       sommets qui sont tous les mots de d lettres sur un alphabet de
       b lettres. Le sommet (x_1,â€¦,x_d) est voisin des sommets
       (x_2,â€¦,x_d,*). Ce graphe est Hamiltonien, de diamÃ¨tre d et le
       degrÃ© de chaque sommet est 2b, 2b-1 ou 2b-2. Pour d=3 et b=2,
       le graphe est planaire.

....kautz d b
....
       Graphe de Kautz de dimension d>0 et de base b>1. Il a
       bÂ·(b-1)^(d-1) sommets qui sont tous les mots de d lettres sur
       un alphabet de b lettres avec la contrainte que deux lettres
       consÃ©cutives doivent Ãªtre diffÃ©rentes. L'adjacence est celle du
       graphe de De Bruijn. C'est donc un sous-graphe induit de De
       Bruijn (debruijn d b). Il est Hamiltonien, de diamÃ¨tre d et le
       degrÃ© de chaque sommet est 2b-2 ou 2b-3. Pour d=b=3 le graphe
       est planaire.

....linial n t
....
       Neighborhood graph des cycles introduit par Linial. C'est le
       graphe de voisinage des vues de taille t d'un cycle orientÃ©
       symÃ©trique Ã  n sommets ayant des identifiants uniques de [0,n[.
       Il faut nâ‰¥t>0 et nâ‰¥2. Les sommets sont les t-uplets d'entiers
       distincts de [0,n[. Le sommet (x_1,â€¦,x_t) est voisin des
       sommets (x_2,â€¦,x_t,y) oÃ¹ yâ‰ x_1 si n>t et y=x_1 si n=t. Le
       nombre chromatique de ce graphe est k ssi il existe un
       algorithme distribuÃ© qui en temps t-1 (resp.  en temps (t-1)/2
       avec t impair) peut colorier en k couleurs tout cycle orientÃ©
       (resp. orientÃ© symÃ©trique) Ã  n sommets ayant des identifiants
       uniques et entiers de [0,n[. C'est un sous-graphe induit de
       "linialc n t", et donc du graphe de Kautz (kautz t n) et de De
       Bruijn (debruijn t). Le nombre de sommets est nÂ·(n-1)â”…(n-t+1).
       Certaines propriÃ©tÃ©s se dÃ©duisent du graphe linialc n t. Pour
       n=4 et t=2, il s'agit du cuboctaÃ¨dre.

....linialc m t
....
       Neighborhood graph des cycles colorÃ©s.  Il s'agit d'une
       variante du graphe linial n t. La diffÃ©rence est que les
       sommets du cycle n'ont plus forcÃ©ment des identitÃ©s uniques,
       mais seulement une m-coloration avec mâ‰¤n. Il faut mâ‰¥tâ‰¥0 et
       mâ‰¥2. L'adjacence est identique, mais les sommets sont les
       t-uplets (x_1,â€¦,x_t) d'entiers de [0,m[ tels que x_iâ‰ x_{i+1}.
       Il s'agit donc d'un sous-graphe induit de linialc m t, lui-mÃªme
       sous-graphe induit du graphe de Kautz (kautz t m) et donc de De
       Bruijn (debruijn t m). Le nombre de sommets est mÂ·(m-1)^{t-1}
       et son degrÃ© est â‰¤ 2Â·(m-1). La taille de la clique maximum est
       3 si m>2 et t>1. Le nombre chromatique de ce graphe pour t=3
       est 3 pour m=4, 4 pour 5â‰¤mâ‰¤24. Pour 25â‰¤mâ‰¤70 c'est au moins 4 et
       au plus 5, la valeur exacte n'Ã©tant pas connue. Tout comme
       "linial 4 2", pour m=4 et t=2, il s'agit du cuboctaÃ¨dre.

....pancake n
....
       Graphe "pancake" de dimension n. Il a n! sommets qui sont les
       permutations de {1,â€¦,n} et (n-1)-rÃ©gulier. Une permutation,
       c'est-Ã -dire un sommet, est voisine de toutes celles obtenues
       en retournant un de ces prÃ©fixes. Plus prÃ©cisÃ©ment, les sommets
       x=(x_1,â€¦,x_n) et y=(y_1,â€¦,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=x_{k-i} sinon.
       Son diamÃ¨tre, qui est linÃ©aire en n, n'est pas connu
       prÃ©cisÃ©ment. Les premiÃ¨res valeurs connues, pour n=1â€¦17,
       sont: 0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18,
       19. Donc les diamÃ¨tres 2,6,12 n'existent pas.

....bpancake n
....
       Graphe "burn pancake" de dimension n. Il a n!Â·2^n sommets qui
       sont les permutations signÃ©es de {1,â€¦,n}. Les sommets
       x=(x_1,â€¦,x_n) et y=(y_1,â€¦,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=-x_{k-i}
       sinon. Dit autrement la permutation de y doit Ãªtre obtenue en
       retournant un prÃ©fixe de x et en inversant les signes. Par
       exemple, le sommet (+2,-1,-5,+4) est voisin du sommet
       (+5,+1,-2,+4). Comme le graphe pancake, c'est un graphe
       (n-1)-rÃ©gulier de diamÃ¨tre linÃ©aire en n.

....gpstar n d
....
       Graphe "permutation star" gÃ©nÃ©ralisÃ© de dimension n. Il a n!
       sommets qui sont les permutations de {1,â€¦,n}. Deux sommets sont
       adjacents si leurs permutations diffÃ¨rent par d positions. Si
       d<2, il s'agit d'un stable. C'est un graphe rÃ©gulier.

....pstar n
....
       Graphe "permutation star" de dimension n. Il a n! sommets qui
       sont les permutations de {1,â€¦,n}. Deux sommets sont adjacents
       si une permutation est obtenue en Ã©changeant le premier Ã©lÃ©ment
       avec un autre. Le graphe est (n-1)-rÃ©gulier. Le graphe est
       biparti et de diamÃ¨tre âŽ£ 3(n-1)/2)âŽ¦. C'est un sous-graphe induit
       d'un "gpstar n 1".

....hexagon p q
....
       Grille hexagonale p x q. C'est un planaire composÃ© de p rangÃ©es
       de q hexagones, le tout arrangÃ© comme un nid d'abeille. Ce
       graphe peut aussi Ãªtre vu comme un mur de p rangÃ©es de q
       briques, chaque brique Ã©tant reprÃ©sentÃ©e par un cycle de
       longueur 6. Il possÃ¨de (p+1)Â·(2p+2)-2 sommets et est de degrÃ©
       maximum 3. Sont dual est le graphe whexagon.
....
       Ex: gengraph hexagon 20 20 -dele 0.2 -maincc -visu

....whexagon p q
....
       Comme le graphe hexagon p q sauf que chaque hexagone est
       remplacÃ© par une roue de taille 6 (chaque hexagone possÃ¨de un
       sommet connectÃ© Ã  ses 6 sommets). C'est le dual de l'hexagone.
       Il possÃ¨de pÂ·q sommets de plus que l'hexagone p q.

....hanoi n b
....
       Graphe de HanoÃ¯ gÃ©nÃ©ralisÃ©, le graphe classique est obtenu avec
       b=3. Il est planaire avec b^n sommets et est dÃ©fini de maniÃ¨re
       rÃ©cursive comme suit. Le niveau n>0 est obtenu en faisant b
       copies du niveau n-1 qui sont connectÃ©s comme un cycle par une
       arÃªte, le niveau 0 Ã©tant le graphe Ã  un sommet. On obtient le
       graphe de sierpinski n b en contractant ces arÃªtes lÃ . Il faut
       bâ‰¥2 et nâ‰¥0. Lorsque n=2, on obtient un sorte de fleur, pour
       n=1 c'est un cycle et pour b=2 il s'agit d'un chemin.

....sierpinski n b
....
       Graphe de SierpiÅ„ski gÃ©nÃ©ralisÃ©, le graphe classique, le
       triangle SierpiÅ„ski qui est planaire, est obtenu avec b=3. Il a
       ((b-2)Â·b^n+b)/(b-1) sommets et est dÃ©fini de maniÃ¨re rÃ©cursive
       comme suit.  Le niveau n est obtenu en faisant b copies du
       niveau n-1 qui sont connectÃ©s comme un cycle, le niveau 1 Ã©tant
       un cycle de b sommets. Il faut bâ‰¥3 et nâ‰¥1. Ã€ la diffÃ©rence du
       graphe d'HanoÃ¯, les arÃªtes du cycle sont contractÃ©es. Le graphe
       de HajÃ³s est obtenu avec n=2 et b=3. Pour n=1 il s'agit d'un
       cycle.

....banana n k
....
       Arbre Ã  nÂ·(k+1)+1 sommets composÃ©s de n>0 copies d'Ã©toiles Ã  k
       branches connectÃ©es, par une feuille, Ã  un unique sommet. Si
       k=0, il s'agit d'un stable Ã  k+1 sommets.

....moser
....
       Graphe "Moser spindle" dÃ©couvert par les frÃ¨res Moser. C'est un
       "unit distance graph" du plan (deux points sont adjacents s'ils
       sont Ã  distance exactement 1) de nombre chromatique 4. Il est
       planaire et possÃ¨de 7 sommets. On ne connaÃ®t pas d'unit
       distance graphe avec un nombre chromatique supÃ©rieur. C'est
       aussi le complÃ©mentaire du graphe K_{3,3} dont une arÃªte a Ã©tÃ©
       subdivisÃ©e.

....markstrom
....
       Graphe de MarkstrÃ¶m. Il est cubique planaire Ã  24 sommets. Il
       n'a pas de cycle de longueur 4 et 8.

....robertson
....
       Graphe de Robertson. C'est le plus petit graphe 4-rÃ©gulier de
       maille 5. Il a 19 sommets, est 3-coloriable et de diamÃ¨tre 3.

....wiener-araya
....
       Graphe dÃ©couvert en 2009 par Wiener & Araya. C'est le plus
       petit graphe hypo-Hamiltonien planaire connu, c'est-Ã -dire
       qu'il n'a pas de cycle Hamiltonien mais la suppression de
       n'importe quel sommet le rend Hamiltonien. Il possÃ¨de 42
       sommets, 67 arÃªtes, et est de diamÃ¨tre 7.

....zamfirescu
....
       Graphe de Zamfirescu Ã  48 sommets dÃ©couvert en 2007. Il est
       planaire et hypo-Hamiltonien. C'est le second plus petit (voir
       wiener-araya). Il possÃ¨de 76 arÃªtes et a un diamÃ¨tre de 7.

....hatzel
....
       Graphe de Hatzel. Il est planaire, de diamÃ¨tre 8, possÃ¨de 57
       sommets et 88 arÃªtes, et est hypo-Hamiltonien (voir
       wiener-araya). C'Ã©tait le plus petit planaire hypo-Hamiltonien
       connu avant le graphe de Zamfirescu.

....clebsch n
....
       Graphe de Clebsch d'ordre n. Il est construit Ã  partir d'un
       hypercube de dimension n en ajoutant une arÃªte entre chaque
       paire de sommets opposÃ©s, c'est-Ã -dire Ã  distance n. Le graphe
       classique de Clebsch est rÃ©alisÃ© pour n=4 dont le diamÃ¨tre est
       deux.

....gear n
....
       Graphe planaire Ã  2n+1 sommets composÃ© d'une roue Ã  n rayons
       ("wheel n") et de n sommets chacun connectÃ©s deux sommets
       voisins consÃ©cutifs du bord de la roue. Il est construit Ã 
       partir du graphe "cage 2n 2 0 ." auquel on ajoute un sommet
       central. Pour n=3, c'est le complÃ©mentaire de "helm 3".

....helm n
....
       Graphe planaire Ã  2n+1 sommets composÃ© d'une roue Ã  nâ‰¥3 rayons
       ("wheel n") et de n sommets pendants connectÃ©s au bord de la
       roue. Pour n=3, c'est le complÃ©mentaire de "gear 3".

....haar n
....
       Graphe de Haar H(n) pour l'entier n>0. C'est un graphe biparti
       rÃ©gulier possÃ©dant 2k sommets oÃ¹ k=1+âŽ£ logâ‚‚(n)âŽ¦ est le nombre
       de bits dans l'Ã©criture binaire de n. Ses sommets sont les u_i
       et v_i pour i=0,â€¦,k-1. Le sommet u_i est adjacent Ã  v_{i+j mod
       k} ssi le bit j de n vaut 1 (j=0,â€¦,k-1). Si n est impair, H(n)
       est connexe et de maille 4 ou 6. La valeur maximale de n est
       2^32-1 = 4,294,967,295 correspondant Ã  un graphe de 64
       sommets. On retrouve respectivement les graphes de Franklin
       (n=37), de Heawood (n=69) et de MÃ¶bius-Kantor (n=133). On a
       aussi que H(2^n-1) est le biparti K_{n,n}, H(2^n) est "matching
       n+1", H(2^n+1) est un cycle Ã  n+1 sommets, H(2^n+3) est le
       "mobius n+1" si n est pair et le "prism n+1" si n est impair,
       et H(3Â·2^n-1) est le "crown n+2".

....turan n r
....
       Graphe de TurÃ¡n Ã  n sommets et r parts. Il s'agit d'un graphe
       r-parti complet de n sommets avec r-(n mod r) parts de âŽ£ n/râŽ¦
       sommets et (n mod r) parts de âŽ¡ n/râŽ¤ sommets. Il faut nâ‰¥r>0. On
       peut aussi le dÃ©finir comme le graphe ayant une arÃªte entre i
       et j ssi |i-j| mod r â‰  0. C'est la dÃ©finition utilisÃ©e pour
       gÃ©nÃ©rer ce graphe. Il est rÃ©gulier lorsque r divise n. Il
       possÃ¨de âŽ£(r-1)n^2/(2r)âŽ¦ arÃªtes et est de nombre chromatique
       r. C'est le graphe sans clique de taille r+1 ayant le plus
       grand nombre d'arÃªtes. Lorsque n=r, il s'agit d'une
       clique. Lorsque n=r+1, il s'agit d'une clique moins une
       arÃªte. Lorsque n=2r, il s'agit du "cocktail party graph" et
       pour n=8 et r=4 le graphe est 1-planar. Lorsque n=3r, il s'agit
       du graphe de Moon-Moser, le graphe possÃ©dant le plus grand
       nombre de cliques maximales (soit 3^(n/3)). Lorsque n=6 et r=3,
       c'est l'octaÃ¨dre.

....klein p q
....
       Maillage quadrangulaire de genre un, incluant le tore, la
       bouteille de Klein et le plan projectif. C'est un graphe
       4-rÃ©gulier composÃ© d'une grille |p|Ã—|q| augmentÃ©e d'arÃªtes
       connectant les bords opposÃ©s. La connexion est torique ou en
       twist suivant le signe de chaque dimension (>0 pour torique et
       <0 pour twist). Pour la bouteille de Klein, il faut pÂ·q<0, et
       pour le plan projectif il faut p<0 et q<0. Le nombre
       chromatique est 4 si |p|â‰ 1 et |q|â‰ 1 (cas d'un cycle) et qu'il y
       a une dimension <0 impaire, ou que pÂ·q=-4 (Kâ‚„). Sinon il est <
       4. Les plus petits exemples avec un nombre chromatique 4, Ã 
       part Kâ‚„, sont "klein 3 -3" et "klein -3 -3" (qui de plus est
       sans Kâ‚ƒ). Par dÃ©faut, les sommets sont dessinÃ©s selon une
       grille |p|Ã—|q|.
....
       Ex: gengraph klein -3 -4 -label -1 -xy noise .3 .7 -dot len 1 -visu

....flower_snark n
....
       Graphe cubique Ã  4n sommets construit de la maniÃ¨re suivante:
       1) on part de n Ã©toiles disjointes Ã  3 feuilles, la i-Ã¨me ayant
       pour feuilles les sommets notÃ©s u_i,v_i,w_i, i=1â€¦n; 2) pour
       chaque xâˆˆ{u,v,w}, x_1-â‹¯-x_n induit un chemin; et enfin 3) sont
       adjacents: u_0-u_n, v_0-w_n et w_0-v_n. Pour n>1, ces graphes
       sont non planaires, non Hamiltoniens, 3-coloriables et de
       maille au plus 6. Pour n=1, il s'agit d'un K_{1,3} (claw).

....udg n r
....
       Graphe gÃ©omÃ©trique alÃ©atoire (random geometric graph) sur n
       points du carrÃ© [0,1[Â² (distribution par dÃ©faut). Deux sommets
       sont adjacents si leurs points sont Ã  distance â‰¤ r. Il s'agit
       de la distance selon la norme L2 (par dÃ©faut), mais cela peut
       Ãªtre changÃ©e par l'option -norm. Le graphe devient connexe avec
       grande probabilitÃ© lorsque r=rc ~ âˆš(ln(n)/n). Si r<0, alors le
       rayon est initialisÃ© Ã  rc. Un UDG (unit disk graph) est
       normalement un graphe d'intersection de disques fermÃ©s de rayon
       1.

....gabriel n
....
       Graphe de Gabriel. Graphe gÃ©omÃ©trique dÃ©fini Ã  partir d'un
       ensemble de n points du carrÃ© [0,1[Â² (distribution par
       dÃ©faut). Les points i et j sont adjacents ssi le plus petit
       disque (voir -norm) passant par i et j ne contient aucun autre
       point. Ce graphe est connexe et planaire Ã  condition toutefois
       qu'ils n'existe pas 4 points co-cycliques et par paires
       diamÃ©tralement opposÃ©es (dans ce cas une clique et des
       croisements d'arÃªtes apparaissent). C'est un sous-graphe du
       graphe de Delaunay. Son Ã©tirement est non bornÃ©.

....rng n
....
       Graphe du proche voisinage (Relative Neighborhood Graph).
       Graphe gÃ©omÃ©trique dÃ©fini Ã  partir d'un ensemble de n points du
       carrÃ© [0,1[Â². Les points i et j sont adjacents ssi il n'existe
       aucun point k tel que max{d(k,i),d(k,j)} < d(i,j) oÃ¹ d est la
       distance (L2 par dÃ©faut, voir -norm). Dit autrement, la "lune"
       dÃ©finie par i et j doit Ãªtre vide. Ce graphe est planaire et
       connexe. C'est un sous-graphe du graphe de Gabriel.

....knng n k
....
       Graphe des k plus proches voisins (k-Nearest Neighbor
       Graph). Graphe gÃ©omÃ©trique dÃ©fini Ã  partir d'un ensemble de n
       points du carrÃ© [0,1[Â² (distribution par dÃ©faut). Chaque point
       i est connectÃ© aux k plus proches autres points (par dÃ©faut
       selon la norme L2, voir -norm). Si la norme est L2, le degrÃ©
       des sommets est â‰¤ 6k. Il peut comporter des croisements
       d'arÃªtes dÃ¨s que k>1. Cependant, chaque arÃªte ne peut Ãªtre
       coupÃ©e que O(kÂ²) fois. Plus prÃ©cisÃ©ment, c'est un t-planaire
       pour tâ‰¤78k^2-6k (cf. [DMW19]). Ã€ partir de k=3 apparaÃ®t une
       composante connexe de taille linÃ©aire.

....mst n
....
       Graphe alÃ©atoire gÃ©omÃ©trique dÃ©finissant un arbre couvrant de
       poids minimum du graphe complet Euclidien sur n points tirÃ©s
       alÃ©atoirement uniformÃ©nent dans le carrÃ© [0,1[Â² (distribution
       par dÃ©faut). Par dÃ©faut la distance est la norme L2 (voir
       -norm).
....
       Ex: gengraph mst 2000 -xy seed 1 .3 -visu

....thetagone n p k w
....
       Graphe gÃ©omÃ©trique dÃ©fini Ã  partir d'un ensemble de n points du
       carrÃ© [0,1[Â² (distribution par dÃ©faut). En gÃ©nÃ©ral le graphe
       est planaire et connexe avec des faces internes de longueur au
       plus p (pour k diviseur de p et w=1). On peut interprÃ©ter les
       paramÃ¨tres comme suit: pâ‰¥3 est le nombre de cotÃ©s d'un polygone
       rÃ©gulier, kâ‰¥1 le nombre d'axes (ou de direction), et wâˆˆ[0,1] le
       cÃ´ne de visibilitÃ©. Toute valeur de p<3 est interprÃ©tÃ©e comme
       une valeur infinie, et le polygone rÃ©gulier correspondant
       interprÃ©tÃ© comme un cercle. L'adjacence entre une paire de
       sommets est dÃ©terminÃ©e en temps O(kn).
....
       Plus formellement, pour tout point u et v, et entier i, on note
       P_i(u,v) le plus petit p-gone (polygone convexe rÃ©gulier Ã  p
       cotÃ©s) passant par u et v dont u est un sommet, et dont le
       vecteur allant de u vers son centre forme un angle de iÂ·2ðœ‹/k
       avec l'axe des abscisses, intersectÃ© avec un cÃ´ne de sommet u
       et d'angle wÂ·(p-2)Â·ðœ‹/p (wÂ·ðœ‹ si p est infini) et dont la
       bissectrice passe par le centre du p-gone. Alors, u est voisin
       de v s'il un existe au moins un entier iâˆˆ[0,k[ tel que
       l'intÃ©rieur de P_i(u,v) est vide. La distance entre u et le
       centre du p-gone dÃ©finit alors une distance (non symÃ©trique) de
       u Ã  v.
....
       Si w=1 (visibilitÃ© maximale), P_i est prÃ©cisÃ©ment un p-gone. Si
       w=0 (visibilitÃ© minimale), P_i se rÃ©duit Ã  l'axe d'angle iÂ·2ðœ‹/k
       pour un entier i. Si w=.5, P_i est un cÃ´ne formant un angle
       Ã©gale Ã  50% de l'angle dÃ©fini par deux cotÃ©s consÃ©cutifs du
       p-gone, ce dernier angle valant (p-2)ðœ‹/p. Si w=2p/((p-2)k) (ou
       simplement 2/k si p est infini) alors la visibilitÃ© correspond
       Ã  un cÃ´ne d'angle 2ðœ‹/k, l'angle entre deux axes. Comme il faut
       wâ‰¤1, cela implique que kâ‰¥2p/(p-2) (kâ‰¥2 si p infini). On
       retrouve le Theta_k-Graph pour chaque kâ‰¥6 en prenant p=3 et
       w=6/k, le demi-Theta-Graph pour tout kâ‰¥3 en prenant p=3 et
       w=3/k, le Yao_k-Graph pour chaque kâ‰¥2 en prenant p=0 (infini)
       et w=2/k, et la triangulation de Delaunay si p=0 (infini), k
       trÃ¨s grand et w=1. En fait, ce n'est pas tout-Ã -fait le graphe
       Yao_k, pour cela il faudrait que u soit le centre du polygone
       (c'est-Ã -dire du cercle).

....pat p q r
....
       Graphe possÃ©dant pqr sommets, issu d'un jeu Ã  un joueur proposÃ©
       par Pat Morin (Barbade, mars 2016). Le jeu se dÃ©roule sur une
       grille pÃ—q et comprend r coups. Un coup est un ensemble de
       positions de la grille strictement croissantes (coordonnÃ©es en
       x et en y strictement croissantes). De plus, si la position
       (x,y) est jouÃ©e alors toutes les positions situÃ©es sur la mÃªme
       ligne mais avec une abscisse au moins x ou sur la mÃªme colonne
       mais avec une ordonnÃ©es au moins y sont interdites pour tous
       les coups suivants. Le score est le nombre total de positions
       jouÃ©es en r coups. Il s'agit de trouver le score maximum.
       Lorsque r=1, le score maximum vaut min(p,q). Lorsque p=q=n et
       r=2, alors le score maximum vaut âŽ£ 4n/3âŽ¦. La question est
       ouverte lorsque r>2, c'est au moins n^1.516 pour r=n oÃ¹ la
       constante vaut log_9(28).
....
       Les sommets du graphes sont les positions dans les r grilles
       pÃ—q et deux sommets sont adjacents les positions sont en
       conflits. Le score du jeu est alors un ensemble indÃ©pendant du
       graphe. Si r=1, le graphe est une grille pÃ—q. Ce graphe active
       l'option -pos car un dessin de ce graphe (sous forme de
       grilles) est proposÃ©.
....
       Ex: gengraph pat 4 4 4 -check kindepsat 8 | ./glucose -model

....line-graph n k
....
       Line-graphe alÃ©atoire Ã  n sommets et de paramÃ¨tre k>0 entier
       Plus k est petit, plus le graphe est dense, le nombre d'arÃªtes
       Ã©tant proportionnel Ã  (n/k)Â². Si k=1, il s'agit d'une clique Ã 
       n sommets. Ces graphes sont obtenus en choisissant, pour chaque
       sommet, deux couleurs de [0,k[. Deux sommets sont adjacents ssi
       ils possÃ¨dent la mÃªme couleur. Il contient le graphe "uno n k
       k". Ces graphes sont claw-free (sans K_{1,3} induit). Tout
       line-graphe est claw-free, et les line-graphes connexes avec un
       nombre pair de sommets possÃ¨dent toujours un couplage
       parfait. On rappel qu'un graphe G est le line-graphe d'un
       graphe H si les sommets de G correspondent aux arÃªtes de H et
       oÃ¹ deux sommets de G sont adjacents ssi les arÃªtes
       correspondantes dans H sont incidentes. On parle parfois de
       graphe adjoint.

....uno n p q
....
       Line-graphe alÃ©atoire Ã  n sommets issu d'un graphe biparti de
       parts de taille p>0 et q>0. Plus prÃ©cisÃ©ment, les sommets sont
       des paires (i,j) d'entiers alÃ©atoires de [0,p[ Ã— [0,q[, pas
       nÃ©cessairement distinctes. Les sommets (i,j) et (i',j') sont
       adjacents ssi i=i' ou j=j'. C'est un sous-graphe induit du
       produit cartÃ©sien de deux cliques, K_p Ã— K_q. Ce graphe est
       gÃ©omÃ©trique, les sommets Ã©tant des points de la grille pÃ—q. Les
       sommets reprÃ©sentent aussi des cartes du jeu de UNO et les
       arÃªtes indiquent si un carte peut Ãªtre jouÃ©e consÃ©cutivement Ã 
       une autre. Le graphe "uno n k k" est un sous-graphe de
       "line-graph n k".

....unok n p q k_p k_q
....
       Graphe "uno n p q" particulier oÃ¹ les n points correspondant
       aux sommets sont pris uniformÃ©ment parmi les ensembles de n
       points distincts de [0,p[ Ã— [0,q[ ayant au plus k_p sommets par
       ligne et k_q par colonne. Il faut n â‰¤ min{pÂ·k_p,qÂ·k_q} et p, q,
       k_p, k_q>0. Si k_p<0, alors on fait comme si k_p=p, de mÃªme
       pour k_q=q si k_q<0. Contrairement Ã  uno, deux sommets ont
       toujours des coordonnÃ©es distinctes. Le graphe rÃ©sultant est de
       degrÃ© au plus k_p+k_q-2, et est de path-width (et aussi de
       tree-width) au plus celle du produit de clique K_{k_p} Ã—
       K_{k_q} soit environ k_pÂ·k_q/2. Le temps de gÃ©nÃ©ration des n
       points est en O(npq) contre O(n) pour uno, mais une
       optimisation (algorithme par rejets) fait qu'il est trÃ¨s
       souvent en O(n+p+q), dans les cas peu dense par exemple. Si k_p
       ou k_q=1, le graphe est une union de cliques, et si k_p=k_q=2
       et n=2p=2q, c'est une union de cycles.
....
       Ex: gengraph unok 200 100 100 3 2 -visu

....wpsl  n p q
....upsl  n p q
....wpsld n p q
....upsld n p q
....
       Weighted/Uniform Planar Stochastic Lattice. Graphe planaire
       alÃ©atoire connexe dont les sommets correspondent Ã  certains
       points d'une grille pÃ—q et les arÃªtes Ã  des lignes horizontales
       ou verticales. Il est gÃ©nÃ©rÃ© selon un processus en nâ‰¥0 Ã©tapes
       dÃ©crit ci-aprÃ¨s. Il faut pâ‰¥2, qâ‰¥2. Il comprend au plus 3n+1
       faces internes rectangulaires, appelÃ©s blocs, qui forment une
       partition des cases de la grille pÃ—q. C'est un graphe
       2-dÃ©gÃ©nÃ©rÃ© qui possÃ¨de au plus min{5n+4,pÂ·q} sommets dont 4
       sont de degrÃ© 2, les autres Ã©tant de degrÃ© 3 ou 4 selon une
       rÃ©partition moyenne 80%-20%. La variante wpsld (ou upsld)
       reprÃ©sente le graphe dual qui est 4-dÃ©gÃ©nÃ©rÃ© et connexe. Les
       sommets (au plus 3n+1) correspondent aux blocs (positionnÃ©s en
       leurs centres), deux blocs Ã©tant adjacents s'ils ont un bord en
       commun. Pour le dual les coordonnÃ©es des centres des blocs sont
       doublÃ©es pour Ãªtre entiÃ¨res, et le dessin n'est pas forcÃ©ment
       planaire. L'option -directed permet d'obtenir une 2-orientation
       et une 4-orientation pour le dual.
....
       Ex: gengraph -seed 0 wpsl 100 400 500 -dot scale auto -visu
           gengraph -seed 0 wpsld 100 400 500 -dot scale auto -visu
	   gengraph wpsl 2500 70000 70000 -dot scale auto -visu
	   gengraph wpsl 2 10 10 -dot scale auto -xy grid 10 -visu
....
       Le graphe est construit en n Ã©tapes. Au dÃ©part il y a un seul
       bloc contenant toutes les cases d'une grille pÃ—q. Le bord de ce
       bloc correspond aux 4 coins de la grille et forme un cycle de
       longueur 4. Ã€ chacune des n Ã©tapes ont sÃ©lectionne un bloc B
       parmi ceux dÃ©jÃ  construits selon une probabilitÃ©
       proportionnelle de sa surface (pour wpsl) ou uniformÃ©ment parmi
       tous les blocs (pour upsl). Ici la surface d'un bloc est le
       nombre de cases -- et non de points -- de la grille qu'il
       contient. Puis B est dÃ©coupÃ© selon une croix dont le centre est
       un point de la grille interne de B (pas sur un bord) choisi
       alÃ©atoirement uniformÃ©ment. (La variante consistant Ã  dÃ©couper
       un bloc en deux correspond au graphe wdis et ses variantes.) Le
       bloc n'est pas dÃ©coupÃ© s'il ne contient pas de points de la
       grille. Une autre faÃ§on de concevoir le processus pour wpsl est
       de choisir n points de la grille pÃ—q alÃ©atoirement et
       uniformÃ©ment. Puis, depuis chaque point et dans un ordre
       quelconque, faire pousser une croix jusqu'Ã  atteindre un bord
       ou une croix prÃ©cÃ©dante, les points se trouvant sur le passage
       d'une croix Ã©tant supprimÃ©s. La crÃ©ation du graphe prend un
       temps O(nlogn) en moyenne, indÃ©pendant des dimensions p et q
       qui peuvent Ãªtre donc relativement grandes. Ensuite, le test
       d'adjacence, liÃ© Ã  sa k-orientation, est constant.

....wdis  n p q
....udis  n p q
....wdisd n p q
....udisd n p q
....
       Rectangular Dissection. Graphe planaire alÃ©atoire connexe dont
       les sommets correspondent Ã  certains points d'une grille pÃ—q et
       les arÃªtes Ã  des lignes horizontales ou verticales. Il est
       gÃ©nÃ©rÃ© selon un processus en nâ‰¥0 Ã©tapes similaires Ã  wpsl (et
       ses variantes). Les variantes wdisd et udisd correspondent au
       graphe dual. Ã€ chaque Ã©tape du processus on sÃ©lectionne
       alÃ©toirement un bloc, soit proportionnellement sa surface
       (wdis) soit uniformÃ©ment (udis), que l'on le coupe en deux
       sous-blocs. Le sens de la dÃ©coupe (verticalement ou
       horizontalement) est soit alÃ©atoire uniforme (udis) soit selon
       une probabilitÃ© porprotionnelle Ã  la longueur des cotÃ©s (wdis),
       prÃ©fÃ©rant dÃ©couper le plus grand cotÃ©. Il possÃ¨de au plus
       min{2n+4,pÂ·q} sommets: 4 de degrÃ© 2, presque tous les autres de
       degrÃ© 3 sauf quelqu'uns de degrÃ© 4. Le dual possÃ¨de n+1
       sommets. Il partage un grand nombres de propriÃ©tÃ©s communes
       avec wpsl, en particulier l'orientation.  Tous les graphes wpsl
       peuvent Ãªtre gÃ©nÃ©rÃ©s par un wdis.

....ngon p c x
....
       Triangulation d'un polygone rÃ©gulier. Plusieurs types de
       triangulations sont produites suivant la valeur des
       paramÃ¨tres. Si x&4=0, alors la triangulation a 3p sommets (avec
       p>0) et est composÃ©e d'un triangle Ã©quilatÃ©ral central. Si
       x&4=1, alors la triangulation a 4p sommets et est composÃ©e d'un
       carrÃ© central avec une diagonale. La triangulation est
       symÃ©trique dans chacun des 3 ou 4 croissants dÃ©limitÃ©s par
       chacune des arÃªtes du polygone central. Si AB est l'une de ces
       arÃªtes (A avant B dans le sens direct), alors on note C le
       point de l'arc de cercle de A Ã  B Ã  distance c de A. On doit
       avoir câˆˆ[0,p/2]. Les deux bits de poids faible de x dÃ©finissent
       comment sont construits les triangulations de l'arc AC et
       CB. Tous les points de AC sont connectÃ©s Ã  A si x&1=1 et Ã  C
       sinon. Et, tous les points de BC sont connectÃ©s Ã  B si x&2=1 et
       Ã  C sinon. Si x&8=1 alors la triangulation est asymÃ©trique. On
       remplace c par p-c pour un arc AB sur deux.
....
       Si x=-1, alors il s'agit d'une autre triangulation. Elle a p
       sommets et est symÃ©trique par rapport Ã  un axe horizontal
       comprenant trois "fan": un depuis le point 0 vers tous ceux de
       [c,n-c], un depuis c vers tous ceux de [0,c], et enfin un
       depuis n-c vers tous ceux de [n-c,n].
....
       Si x=-2, alors il s'agit de la triangulation rÃ©cursive Ã  3p
       sommets. Le paramÃ¨tre c n'a pas de rÃ´le. Chacun des trois arc
       est coupÃ© en deux rÃ©cursivement. Si p n'est pas une puissance
       de deux, alors le graphe peut ne pas Ãªtre une triangulation
       complÃ¨te, mais le graphe reste cependant planaire.
....
       La triangulation qui expÃ©rimentalement minimise le stretch
       maximum (voir -check stretch) est obtenue avec "ngon p ð›¼p 3" oÃ¹
       ð›¼ = 231/512 â‰ƒ 45%. Le stretch maximum est environ 1.455 rÃ©alisÃ©
       entre les sommets u=p+30% et v=3p-20%.

....behrend p k
....
       Graphe rÃ©gulier de pÂ·k sommets possÃ©dant un trÃ¨s grands nombre
       de cycles de longueur k arÃªte-disjoints oÃ¹ p,k â‰¥ 2. Si p est
       premier, il en possÃ¨de exactement pÂ·c! = p^{2-o(1)} oÃ¹ c ~
       log(p)/loglog(p). Son degrÃ© est 2c! si k>2 ou c! si k=2. Le
       graphe est dÃ©fini que p soit premier ou pas. Il est construit Ã 
       partir de k stables S_0,â€¦,S_{k-1} de chacun p sommets. Chacun
       des cycles de longueur k contient exactement un Ã©lÃ©ment de
       chaque S_j qui est le sommet d'indice i+jÂ·x (mod p) dans S_j
       avec iâˆˆ[0,p[ et xâˆˆX oÃ¹ XâŠ‚[0,p[ est un ensemble oÃ¹ k entiers
       quelconques ne sont jamais en progression arithmÃ©tique. On
       construit X comme l'ensemble de tous les entiers < p/(k-1)
       s'Ã©crivant sur c chiffres distincts pris dans [0,c[ en base
       ck+1 avec c maximum. Donc |X|=c!. Lorsque p est petit, le
       graphe peut ne pas Ãªtre connexe.
....
       Par exemple, pour p=421 et k=3 on obtient c=3 et X = { 012,
       021, 102, 120, 201, 210 } (nombres Ã©crits en base ck+1=10). On
       vÃ©rifie qu'on a bien p > (k-1)Â·max{X} = 420. Ce graphe et donc
       12-rÃ©gulier possÃ¨de pÂ·k = 1263 sommets et pÂ·c! = 5052 triangles
       arÃªte-disjoints car 421 est premier. La table ci-dessous donne
       en fonction de k et du degrÃ© souhaitÃ© la plus petite valeur de
       p=p(k) possible. Si p est plus petit que p(k), alors le degrÃ©
       sera moindre. Lorsque k=2, le degrÃ© est c! au lieu de 2Â·c!.
....
       !!!           2Â·2!   2Â·3!    2Â·4!       2Â·5!
              degrÃ©    4     12      48        240
    	     â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
               p(2)    6    106   2,359     62,811
               p(3)   15    421  13,885    549,921
               p(4)   28  1,054  46,003  2,419,831
    	     â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
....
       Ces graphes sont utilisÃ©s en "property testing" pour montrer
       qu'il est difficile de dÃ©terminer si un graphe dense possÃ¨de ou
       pas un cycle de longueur k.

....rplg n t
....
       Random Power-Law Graph. Graphe alÃ©atoire Ã  n sommets oÃ¹ les
       degrÃ©s des sommets suivent une loi de puissance d'exposant t>1
       (typiquement un rÃ©el tâˆˆ]2,3[). L'espÃ©rance du degrÃ© du sommet
       i=0â€¦n-1 est w_i = (n/(i+1))^(1/(t-1)). La probabilitÃ© d'avoir
       l'arÃªte i-j est min{w_iÂ·w_j/S,1} avec S=âˆ‘_k w_k. La valeur
       communÃ©ment observÃ©e pour le rÃ©seau Internet Ã©tant t=2.1.

....bdrg n_1 d_1 â€¦ n_k d_k .
....
       Bounded Degree Random Graph. Graphe alÃ©atoire dont la
       distribution des degrÃ©s des sommets est fixÃ©e par les paires
       (náµ¢,dáµ¢) signifiant qu'il y a náµ¢ sommets de degrÃ© au plus
       dáµ¢. Ainsi "bdrg n 3 ." gÃ©nÃ¨re un graphe sous-cubique alÃ©atoire
       Ã  n sommets, si n est pair. Les sommets sont dupliquÃ©s selon
       leur distribution de degrÃ© puis un couplage alÃ©atoire dÃ©termine
       les arÃªtes. Les boucles et les arÃªtes multiples sont supprimer.
       Il suit que le degrÃ© des sommets ne dÃ©passe pas dáµ¢. Ils peuvent
       cependant Ãªtre infÃ©rieurs. Le nombre de sommets est n=âˆ‘áµ¢ náµ¢ et
       le nombre d'arÃªtes au plus m = Â½âˆ‘áµ¢ (náµ¢Â·dáµ¢). Si cette somme
       n'est pas entiÃ¨re, alors le degrÃ© d'un des sommets ayant dáµ¢>0
       est diminuÃ© d'un. (C'est un sommet avec dáµ¢>0 avec le plus grand
       i qui est choisi.)

....fdrg n_1 d_1 â€¦ n_k d_k .
....
       Fixed Degree Random Graph. Graphe alÃ©atoire assymptotiquement
       uniforme dont les degrÃ©s des sommets sont fixÃ©es par les paires
       (náµ¢,dáµ¢) signifiant qu'il y a náµ¢ sommets de degrÃ© dáµ¢. La suite
       des degrÃ©s doit Ãªtre graphique, Ã  savoir qu'il existe au moins
       un graphe simple ayant ces degrÃ©s (sinon une erreur est
       affichÃ©e). Ainsi "fdrg n 3 ." gÃ©nÃ¨re un graphe cubique
       alÃ©atoire asymptotiquement uniforme, Ã  condition que n soit
       pair. Il est possible d'obtenir des graphes non connexe, comme
       avec "fdrg 3 2 1 0 ." composÃ© d'un triangle et d'un sommet
       isolÃ©. La complexitÃ© est en moyenne O(mÎ”+Î”â´) oÃ¹ m=âˆ‘ náµ¢dáµ¢ et
       Î”=max{dáµ¢}, et pour Ãªtre asymptotiquement uniforme, il faut
       Î”=o(m^Â¼) ou Î”=o(âˆšn) pour les graphes rÃ©guliers (tous les dáµ¢
       Ã©gaux ou k=1).

....matching n
....
       Graphe composÃ© de n arÃªtes indÃ©pendantes, c'est-Ã -dire de n
       copies de Kâ‚‚. L'option -directed permet d'obtenir une
       1-orientation.

....load file[:range]
....loadc file[:range]
....
       Graphe dÃ©fini Ã  partir du fichier "file" ou de l'entrÃ©e
       standard si file vaut "-". Si "file" est une famille de
       graphes, alors il est possible d'utiliser la variante
       "file:range" pour prÃ©ciser l'identifiant du graphe souhaitÃ©
       (sinon c'est le premier graphe de la famille qui sera
       considÃ©rÃ©). Le graphe (ou la famille) doit Ãªtre au format
       standard, les sommets numÃ©rotÃ©s par des entiers positifs. Les
       caractÃ¨res situÃ©s sur une ligne aprÃ¨s "//" sont ignorÃ©s, ce qui
       permet de mettre des commentaires.
....
       Le temps et l'espace nÃ©cessaire au chargement du graphe sont
       linÃ©aires en la taille du fichier (si "file" est une famille de
       graphes, le fichier est entiÃ¨rement lu).  Cependant, pour la
       gÃ©nÃ©ration Ã  proprement parlÃ©e du graphe final, qui peut
       comprendre l'option -not par exemple, toutes les arÃªtes
       potentielles, soit O(nÂ²), sont passÃ©es en revue pour Ãªtre
       testÃ©es. La variante "loadc" (pour "load & check") permet une
       gÃ©nÃ©ration plus rapide lorsqu'utilisÃ©e avec -check (ou les
       alias utilisant -check, comme -maincc par exemple). Elle permet
       de passer directement de l'Ã©tape de chargement du graphe Ã 
       l'Ã©tape du test de l'algorithme en sautant la phase de
       gÃ©nÃ©ration des arÃªtes. En contre-partie, le graphe n'est pas
       affichÃ© et les options comme -not, -permute, -delv, -dele,
       etc. n'ont plus d'effet. La variante "loadc file" est environ
       20% plus rapide que "load file -fast".
....
       Pour charger un graphe au format dot on peut utiliser le script
       dot2gen.awk en amont, comme dans l'exemple suivant:
....
       !!! nop file.dot | awk -f dot2gen.awk | ./gengraph load -
....
       Le filtre nop de GraphViz, qui est recommandÃ© mais pas
       nÃ©cessaire, permet de standardiser le format dot initial. Il
       transforme par exemple les expressions du type "a--{b;c;}" en
       "a--b;a--c;".
....
       Notez que la suite d'options "load file -fast -format dot<type>"
       permet de convertir "file" au format <type> souhaitÃ©. Ce graphe
       active l'option -directed si "file" contient au moins un
       arc. Dans ce cas l'option -undirected n'aura pas d'effet.


   GRAPHES ORIENTÃ‰S :

....aqua c_1 â€¦ c_n .
....
       Graphe orientÃ© dont les sommets sont les suites de n entiers
       positifs dont la somme fait c_1 et dont le i-Ã¨me Ã©lÃ©ment est au
       plus c_i. Ils reprÃ©sentent les faÃ§ons de rÃ©partir une quantitÃ©
       c_1 de liquide dans n rÃ©cipients de capacitÃ© c_1 â€¦  c_n. Il y a
       un arc u->v s'ils existent i et j tels que v est le rÃ©sultat du
       versement du rÃ©cipient c_i vers le rÃ©cipient c_j. Le graphe est
       isomorphe au graphe oÃ¹ les c_i=0 ont Ã©tÃ© supprimÃ©s, les c_i ont
       Ã©tÃ© triÃ©s par ordre dÃ©croissant et oÃ¹ c_1 a Ã©tÃ© remplacÃ© par
       min{c_1,c_2+â€¦+c_n}. Par exemple, "aqua 4 1 0 2 ." est isomorphe
       Ã  "aqua 3 2 1 .". Le nombre de sommets ne peut pas dÃ©passer
       binom{n+c_1}{n}. Le graphe peut Ãªtre connexe mais non fortement
       connexe comme "aqua 2 2 .".
....
       Ex: gengraph aqua 3 2 1 . -label 1 -dot filter dot -visu

....collatz n a_0 b_0 â€¦ a_{k-1} b_{k-1} .
....
       Graphe de Collatz gÃ©nÃ©ralisÃ©. Il est basÃ© sur la relation C: x
       â†¦ (aáµ¢Â·x + báµ¢)/k, dÃ©finie pour tout entiers x>0, oÃ¹ i = x%k et
       oÃ¹ aáµ¢,báµ¢ sont entiers (pas forcÃ©ment positifs). Il faut aáµ¢Â·i +
       báµ¢ â‰¡ 0 (mod k) pour tout i pour que C(x) soit entier, sinon
       âŽ£C(x)âŽ¦ est considÃ©rÃ©e et Câ»Â¹ n'est plus forcÃ©ment injective.
       Les entiers gÃ©nÃ©rÃ©s forment les sommets du graphe, les arcs
       Ã©tant les relations xâ†’C(x).
....
       Le graphe pour le problÃ¨me "3x+1", dÃ©fini par la relation xâ†¦x/2
       si xâ‰¡0 (mod 2) et xâ†¦x/2 si xâ‰¡1 (mod 2), est donc le graphe
       "collatz n 1 0 3 1 .".
....
       Si n>0, la boule de volume n est gÃ©nÃ©rÃ©e en itÃ©rant la relation
       inverse depuis x=1, soit Câ»Â¹: x â†¦ (kÂ·x-báµ¢)/aáµ¢ pour chaque i tel
       que aáµ¢ divise kÂ·x-báµ¢. Il est important que C(x) soit entier
       pour que Câ»Â¹ soit injective. Si n<0, la relation est itÃ©rÃ©e
       depuis chaque entier xâˆˆ[1,|n|].  Dans le cas n>0, le graphe est
       connexe et comprend au plus n sommets, alors que pour n<0, il
       peut ne pas Ãªtre connexe et contenir plus de |n| sommets. Il
       s'agit de deux sous-graphes (n>0 et n<0) induit du mÃªme graphe
       infini. Dans tous les cas c'est un graphe orientÃ© avec au plus
       un successeur (arc sortant) et k prÃ©dÃ©cesseurs (arcs entrant),
       et peut contenir des cycles dont des boucles et des arcs
       symÃ©triques.
....
       Ex: gengraph collatz -65 1 0 5 -1 5 1 3 1 . -label 1 -visu
....
       La fameuse conjecture de Collatz pour le problÃ¨me "3x+1"
       affirme que le graphe "collatz n 1 0 3 1 ." est connexe quel
       que soit n<0 (voir aussi "syracuse n"). Pour de grandes
       familles de coefficients aáµ¢,báµ¢, il est conjecturÃ© que le graphe
       possÃ¨de un nombre constant de composantes connexes, comme par
       exemple 3 pour "collatz n 1 0 5 1 .". Savoir si le graphe
       gÃ©nÃ©ralisÃ© de Collatz est connexe est indÃ©cidable, mÃªme si tous
       les báµ¢=0 [Conway'72].
....
       La gÃ©nÃ©ration du graphe proprement dÃ®te est basÃ©e sur une file
       initialisÃ©e aux valeurs 1..|n| pour lesquelles la relation C
       est successivement appliquÃ©e (si n<0). Il s'agit donc d'une
       sorte de parcours en largeur du graphe depuis n sources en
       parallÃ¨le.  La construction est limitÃ©e arbitrairement Ã  nÂ²
       sommets car, suivant les coefficients aáµ¢,báµ¢, ce parcours peut
       ne pas converger pour certain xâˆˆ[1,|n|]. Attention! Les valeurs
       C(x) nÃ©gatives sont exclues du graphe, ce qui peut aussi
       exclure les valeurs dÃ©passant 2Â³Â¹ et donc dÃ©connecter le
       graphe. Si n>0, on initialise la file avec seulement la valeur
       1 et on itÃ¨re les k relations inverses, si elles s'appliquent,
       jusqu'Ã  produire exactement n sommets. Il s'agit donc d'un
       simple parcours en largeur depuis 1. Bien que le nombre d'arcs
       soit linÃ©aire, la gÃ©nÃ©ration de tous les successeurs prend un
       temps quadratique en le nombre de sommets final du graphe si
       n<0, et O(knÂ²) si n>0.


   GRAPHES COMPOSÃ‰S :

....mesh p q (= grid p q .)
....
       Grille 2D de p x q sommets.

....hypercube d (= grid 2 â€¦ 2 .)
....
       Hypercube de dimension d.

....path n (= grid n .)
....
       Chemin Ã  n sommets.

....cycle n (= ring n 1 .)
....
       Cycle Ã  n sommets.

....torus p q (= grid -p -q .)
....
       Tore Ã  p x q sommets.

....stable n (= ring n .)
....empty n (= ring n .)
....
       Stable Ã  n sommets.

....clique n (= -not ring n .)
....
       Graphe complet Ã  n sommets.

....bipartite p q (= rpartite p q .)
....
       Graphe biparti complet K_{p,q}.

....utility (= rpartite 3 3 .)
....
       Graphe biparti complet K_{3,3} qui doit son nom au problÃ¨me de
       la connexion planaire de trois maisons Ã  trois stations (eau,
       gaz, Ã©lectricitÃ©). C'est aussi le graphe de Haar H(7).

....domino (= grid 2 3 .)
....
       Graphe planaire Ã  6 sommets composÃ© de deux carrÃ©s partageant
       une arÃªte.

....kite (= -not banana 1 3)
....
       Kite Graph, graphe Ã  cinq sommets composÃ© de deux triangles
       partageant une arÃªte et d'un sommet pendant (degrÃ© 1) attachÃ© Ã 
       un sommet de degrÃ© deux. C'est aussi un "parachute 2".

....parapluie n (= -not parachute n)
....
       Graphe planaire Ã  n+3 sommets, complÃ©mentaire du parachute. Le
       graphe "parapluie" classique correspond Ã  n=4.

....hourglass (= barbell 3 3 0)
....
       Graphe Papillon ou encore Butterfly. Il a 5 sommets et est
       composÃ© de deux triangles partageant un sommet.

....cuboctahedron (= linial 4 2)
....
       CuboctaÃ¨dre: graphe planaire 4-rÃ©gulier Ã  12 sommets. Il
       possÃ¨de 24 arÃªtes et 14 faces qui sont des triangles ou des
       carrÃ©s. C'est le dual du rhombic-dodÃ©caÃ¨dre.

....octahedron (= antiprism 3)
....
       OctaÃ¨dre: graphe 4-rÃ©gulier planaire Ã  6 sommets ayant 8 faces
       triangulaires. Il s'agit de deux pyramides dont la base Ã  4
       sommets est commune. C'est aussi le graphe de Johnson J(4,2).

....d-octahedron d (= -not matching d)
....
       OctaÃ¨dre de dimension d: obtenu Ã  partir d'un octaÃ¨dre de
       dimension d-1 auquel on ajoute deux sommets universels,
       l'octaÃ¨dre de dimension 1 Ã©tant composÃ© d'un stable de deux
       sommets.  L'octaÃ¨dre classique est obtenu avec d=3, pour d=2 il
       s'agit d'un carrÃ©.

....tetrahedron (= -not ring 4 .)
....
       TÃ©traÃ¨dre: pyramide composÃ©e de 4 faces triangulaires. C'est
       aussi une clique Ã  4 sommets.

....cube (= crown 4)
....hexahedron (= crown 4)
....
       Hypercube de dimension 3, graphe planaire cubique Ã  8 sommets
       oÃ¹ toutes les faces sont des rectangles. C'est aussi un
       hexaÃ¨dre (6 faces carrÃ©es) ou encore le graphe de Haar H(11).

....associahedron (= flip 6)
....
       AssociaÃ¨dre (3D): graphe planaire cubique Ã  14 sommets composÃ©
       de 3 faces carrÃ©es et 6 faces pentagonales.

....johnson n k (= -not kneser n k k-2)
....
       Graphe de Johnson J(n,k). Les sommets sont tous les
       sous-ensembles Ã  k Ã©lÃ©ments de [0,n[ (il faut donc 0â‰¤kâ‰¤n). Deux
       sommets sont adjacents ssi leurs ensembles correspondant ont
       k-1 Ã©lÃ©ments en commun. La distance entre deux sommets est la
       distance de Hamming entre les ensembles correspondant. Ils sont
       rÃ©guliers de degrÃ© k(n-k), de diamÃ¨tre min{k,n-k}, de
       sommet-connectivitÃ© k(n-k). Ils sont aussi distance
       rÃ©guliers. J(n,1) est la clique K_n, J(n,2) est le complÃ©ment
       du graphe de Kneser K(n,2) et le line-graphe de K_n. En fait,
       tout sous-graphe induit de J(n,2) est un line-graphe. J(4,2)
       est l'octaÃ¨dre, J(5,2) le complÃ©ment du graphe de Petersen.

....claw (= rpartite 1 3 .)
....
       Graphe biparti complet K_{1,3}.

....star n (= rpartite 1 n .)
....
       Arbre (Ã©toile) Ã  n feuilles et de hauteur 1.

....tree n (= arboricity n 1)
....
       Arbre plan enracinÃ© alÃ©atoire uniforme Ã  n sommets. Les sommets
       sont numÃ©rotÃ©s selon un parcours en profondeur depuis la racine
       et le long de la face extÃ©rieure.

....caterpillar n (= grid n-r . -star r)
....
       Arbre Ã  n sommets dont les sommets internes (de degrÃ© > 1)
       induisent un chemin. Il est obtenu Ã  partir d'un chemin de
       longueur n-r (oÃ¹ r est un nombre alÃ©atoire entre 0 et n-1) et
       en appliquant l'option -star r. Si l'option -seed est prÃ©sente,
       (pour intervenir sur la valeur "r"), il est important qu'elle
       figure avant caterpillar. [Ce graphe est dÃ©sactivÃ©, voir
       l'option -star.]

....outerplanar n (= kpage n 1)
....
       Graphe planaire-extÃ©rieur alÃ©atoire connexe Ã  n sommets (plan
       et enracinÃ©). Ils sont en bijection avec les arbres plans
       enracinÃ©s dont tous les sommets, sauf ceux de la derniÃ¨re
       branche, sont bicoloriÃ©s. Les sommets sont numÃ©rotÃ©s le long de
       la face extÃ©rieure. C'est aussi une numÃ©rotation selon un
       parcours en profondeur depuis la racine de l'arbre bicoloriÃ©.
       Il est aussi possible de gÃ©nÃ©rer des graphes
       planaires-extÃ©rieurs alÃ©atoires Hamiltoniens, donc 2-connexes,
       avec "planar n f -1" ou "polygon n". L'option -directed permet
       d'obtenir une 2-orientation.

....squaregraph n (= planar n 4 4)
....
       Squaregraph alÃ©atoire Ã  n faces. Ce sont des graphes planaires
       2-connexes dont toutes les faces (sauf l'extÃ©rieure) sont des
       carrÃ©es. De plus, les sommets des faces internes sont de degrÃ©
       au moins 4. Ce sont des sous-graphes de quadrangulations et
       donc des 2-pages. L'option -directed permet d'obtenir une
       2-orientation.

....random n p (= -not ring n . -dele 1-p)
....
       Graphe alÃ©atoire Ã  n sommets et dont la probabilitÃ© d'avoir une
       arÃªte entre chaque paire de sommets est p. L'option -dele Ã©tant
       dÃ©jÃ  prÃ©sente, il n'est pas conseillÃ© de la rÃ©utiliser pour ce
       graphe.

....netgraph (= sierpinski 2 3 -not)
....
       Graphe Ã  6 sommets composÃ© d'un triangle avec un sommet pendant
       Ã  chacun d'eux. C'est le complÃ©mentaire du graphe de HajÃ³s. On
       peut aussi le gÃ©nÃ©rer en utilisant "fdrg 3 3 3 1 .".

....sunflower n (= cage 2n 2 2 .)
....
       Tournesol Ã  n pÃ©tales. C'est un graphe planaire-extÃ©rieur Ã  2n
       sommets composÃ© d'un cycle de longueur nâ‰¥3 oÃ¹ chaque arÃªte
       partage le cotÃ© d'un triangle. C'est le graphe "gear n" sans le
       sommet central. Pour n=3, c'est le graphe de HajÃ³s.

....gem (= fan 4 1)
....
       Graphe Ã  5 sommets composÃ© d'un chemin et d'un sommet universel.

....egraph (= comb 3)
....
       Arbre Ã  6 sommets et 3 feuilles en forme de E.

....tgraph (= banana 1 3)
....fork (= banana 1 3)
....
       Fork Graph, arbre en forme de T Ã  5 sommets dont 3 feuilles.

....ygraph (= banana 3 1)
....
       Arbre Ã  7 sommets composÃ© d'une Ã©toile Ã  trois branches.

....cross (= banana 1 4)
....
       Cross Graph, arbre Ã  six sommets en forme de croix chrÃ©tienne.

....knight p q (= chess p q 1 2)
....
       Graphe des dÃ©placements possible du chevalier dans un Ã©chiquier
       p q.

....antelope p q (= chess p q 3 4)
....
       Graphe des dÃ©placements possible d'une antilope dans un
       Ã©chiquier p q, une antilope Ã©tant une piÃ¨ce hypothÃ©tique se
       dÃ©plaÃ§ant de 3 cases selon un axe et de 4 selon l'autre.

....camel p q (= chess p q 1 3)
....
       Graphe des dÃ©placements possible d'un chameau dans un Ã©chiquier
       p q, un chameau Ã©tant une piÃ¨ce hypothÃ©tique se dÃ©plaÃ§ant de 1
       case selon un axe et 3 de selon l'autre.

....giraffe p q (= chess p q 1 4)
....
       Graphe des dÃ©placements possible d'une giraffe dans un
       Ã©chiquier p q, une giraffe Ã©tant une piÃ¨ce hypothÃ©tique se
       dÃ©plaÃ§ant de 1 case selon un axe et de 4 selon l'autre.

....zebra p q (= chess p q 2 3)
....
       Graphe des dÃ©placements possible d'un zÃ¨bre dans un Ã©chiquier p
       q, un zÃ©bre Ã©tant une piÃ¨ce hypothÃ©tique se dÃ©plaÃ§ant de 2
       cases selon un axe et de 3 selon l'autre.

....petersen (= kneser 5 2 0)
....
       Graphe de Kneser particulier. Il est cubique et possÃ¨de 10
       sommets. Il n'est pas Hamiltonien et c'est le plus petit graphe
       dont le nombre de croisements (crossing number) est 2. C'est le
       complÃ©ment du line-graphe de Kâ‚….

....tietze (= flower_snark 3)
....
       Graphe de Tietze. Il est cubique avec 12 sommets. Il possÃ¨de un
       chemin Hamiltonien, mais pas de cycle. Il peut Ãªtre plongÃ© sur
       un ruban de MÃ¶bius, a un diamÃ¨tre et une maille de 3. Il peut
       Ãªtre obtenu Ã  partir du graphe de Petersen en appliquant une
       opÃ©ration Y-Delta.

....mobius-kantor (= gpetersen 8 3)
....
       Graphe de MÃ¶bius-Kantor. Graphe cubique Ã  16 sommets de genre
       1. Il est Hamiltonien, de diamÃ¨tre 4 et de maille 6. C'est
       aussi le graphe de Haar H(133).

....dodecahedron (= gpetersen 10 2)
....
       DodÃ©caÃ¨dre: graphe planaire cubique Ã  20 sommets. Il possÃ¨de 30
       arÃªtes et 12 faces qui sont des pentagones. C'est le dual de
       l'icosaÃ¨dre.

....desargues (= gpetersen 10 3)
....
       Graphe de Desargues. Il est cubique Ã  20 sommets. Il est
       Hamiltonien, de diamÃ¨tre 5 et de maille 6.

....durer (= gpetersen 6 2)
....
       Graphe de DÃ¼rer. Graphe cubique planaire Ã  12 sommets de
       diamÃ¨tre 4 et de maille 3. Il peut Ãªtre vu comme un cube avec
       deux sommets opposÃ©s tronquÃ©s (remplacÃ©s par un cycle de
       longueur 3).

....prism n (= gpetersen n 1)
....
       Prisme, c'est-Ã -dire le produit cartÃ©sien d'un cycle Ã  n
       sommets et d'un chemin Ã  deux sommets. Pour n=3, c'est un
       graphe de Halin et aussi le complÃ©mentaire d'un cycle de
       longueur 6, et pour n=4 il s'agit du cube.

....cylinder p q (= grid p -q .)
....
       Produit cartÃ©sien d'un chemin Ã  p sommets et d'un cycle Ã  q
       sommets. Cela gÃ©nÃ©ralise le prisme (prism n = cylinder n 3). Un
       cube est un "cylinder 2 4".

....nauru (= pstar 4)
....
       Graphe de Nauru. C'est un graphe cubique Ã  24 sommets. Il
       s'agit d'un graphe "permutation star" de dimension 4. C'est
       aussi un graphe de Petersen gÃ©nÃ©ralisÃ© P(12,5).

....heawood (= cage 14 5 -5 .)
....
       Graphe de Heawood. C'est un graphe cubique bipartis Ã  14
       sommets, de maille 6 et de diamÃ¨tre 3. C'est le graphe
       d'incidence du plan projectif d'ordre 2 (plan de Fano). Il est
       1-planar. C'est le plus petit graphe dont le nombre de
       croisements (crossing number) est 3. C'est aussi le graphe de
       Haar H(69).

....franklin (= cage 12 5 -5 .)
....
       Graphe de Franklin. C'est un graphe cubique Ã  12 sommets, de
       maille 4 et de diamÃ¨tre 3. C'est aussi le graphe de Haar H(37).

....mcgee (= cage 24 12 7 -7 .)
....
       Graphe de McGee. C'est un graphe cubique Ã  24 sommets, de
       maille 7 et de diamÃ¨tre 4.

....bidiakis (= cage 12 -4 6 4 .)
....
       Graphe ou cube de Bidiakis. C'est un graphe planaire cubique Ã 
       12 sommets. Il est Hamiltonien et son nombre chromatique est
       3. On peut le reprÃ©senter comme un cube oÃ¹ deux faces opposÃ©es
       comportent une arÃªte supplÃ©mentaire perpendiculaire joignant
       deux bord opposÃ©s. On peut aussi le reprÃ©senter comme un cycle
       avec 3 colonnes et 3 lignes parallÃ¨les joingnant des sommets
       opposÃ©s (comme une raquette de tennis).

....dyck (= cage 32 5 0 13 -13 .)
....
       Graphe de Dyck. C'est un graphe cubique 3-connexe biparti Ã  32
       sommets. C'est le seul graphe cubique Ã  32 sommets Ã  Ãªtre
       symÃ©trique, c'est-Ã -dire qui est Ã  la fois arÃªte et sommet
       transitif. Il est aussi torique, c'est-Ã -dire de genre 1.

....pappus (= cage 18 5 7 -7 7 -7 5 .)
....
       Graphe de Pappus. C'est un graphe cubique Ã  18 sommets, de
       maille 6 et de diamÃ¨tre 4.

....tutte-coexter (= cage 30 -7 9 13 -13 -9 7 .)
....
       Graphe de Tutte-Coexter appelÃ© aussi 8-cage de Tutte. C'est un
       graphe cubique Ã  30 sommets, de maille 8 et de diamÃ¨tre 4.
       C'est un graphe de Levi mais surtout un graphe de Moore,
       c'est-Ã -dire un graphe d-rÃ©gulier de diamÃ¨tre k dont le nombre
       de sommets est 1+dÂ·S(d,k) (si d impair) ou 2Â·S(d,k) (si d pair)
       avec S(d,k)=âˆ‘_{i=0}^{k-1} (d-1)^i.

....gray (= cage 54 7 -7 25 -25 13 -13 .)
....
       Graphe de Gray. C'est un graphe cubique Ã  54 sommets qui peut
       Ãªtre vu comme le graphe d'incidence entre les sommets d'une
       grille 3Ã—3Ã—3 et les 27 lignes droites de la grille. Il est
       Hamiltonien, de diamÃ¨tre 6, de maille 8, et de genre 7. Il est
       arÃªte-transitif et rÃ©gulier sans Ãªtre sommet-transitif.

....chvatal (= cage 12 3 6 3 6 6 3 6 -3 3 -3 3 3 .)
....
       Graphe de ChvÃ¡tal (1970). C'est le plus petit graphe rÃ©gulier
       sans triangle de nombre chromatique 4. Il est 4-rÃ©gulier,
       possÃ¨de 12 sommets, est non-planaire, Hamiltonien et de
       diamÃ¨tre 2.

....grotzsch (= mycielski 4)
....
       Graphe de GrÃ¶tzsch. C'est le plus petit graphe sans triangle de
       nombre chromatique 4 (sans Ãªtre rÃ©gulier contrairement au
       graphe de ChvÃ¡tal). Il possÃ¨de 11 sommets et 20 arÃªtes. Comme
       le graphe de ChvÃ¡tal, il est non-planaire de diamÃ¨tre 2, de
       maille 4 et Hamiltonien. C'est le graphe de Mycielskian du
       cycle Ã  5 sommets.

....hajos (= sierpinski 2 3)
....
       Graphe de HajÃ³s. Il est composÃ© de trois triangles deux Ã  deux
       partageant un sommet distinct. On peut le dessiner comme un
       triangle dans un triangle plus grand. Il est planaire et
       possÃ¨de 6 sommets. C'est un graphe de Sierpinski ou encore le
       complÃ©mentaire d'un "sunlet 3", complÃ©mentaire du "netgraph",
       un "sunflower 3" ou encore "cage 6 2 0 .".

....house (= -not grid 5 .)
....
       Graphe planaire Ã  5 sommets en forme de maison. C'est le
       complÃ©mentaire d'un chemin Ã  5 sommets.

....wagner (= ring 8 1 4 .)
....
       Graphe de Wagner appelÃ© aussi graphe Wâ‚ˆ, un cycle Ã  8 sommets
       oÃ¹ les sommets antipodaux sont adjacents. C'est un graphe
       cubique Ã  8 sommets qui n'est pas planaire mais sans Kâ‚…. C'est
       aussi une Ã©chelle de MÃ¶bius.

....mobius n (= ring n 1 n/2 .)
....
       Ã‰chelle de MÃ¶bius, graphe cubique Ã  n sommets obtenu Ã  partir
       d'un cycle Ã  n sommets dont les sommets opposÃ©s sont
       adjacents. Lorsque n est pair, il s'agit d'un ruban de MÃ¶bius,
       c'est-Ã -dire d'une Ã©chelle dont le premier et dernier barreau
       sont recollÃ©s en sens opposÃ©. Pour nâ‰¤5, il s'agit d'une clique
       Ã  n sommets. Il est donc cubique sauf pour n=1,2,3,5. Lorsque
       nâ‰¥5, le graphe n'est plus planaire, et pour n=8, il s'agit du
       graphe de Wagner.

....ladder n (= grid 2 n .)
....
       Graphe Ã©chelle Ã  n barreaux, soit une grille Ã  2 x n sommets.

....diamond (= fan 2 2)
....
       Clique Ã  quatre sommets moins une arÃªte. C'est un graphe
       allumette, c'est-Ã -dire planaire et distance unitaire.

....gosset (= ggosset 8 2 3 6 -1 .)
....
       Graphe de Gosset. Il est 27-rÃ©gulier avec 56 sommets et 756
       arÃªtes, de diamÃ¨tre, de rayon et de maille 3. Il est
       27-arÃªte-connexe, 27-sommet-connexe et Hamiltonien. C'est
       localement un graphe de SchlÃ¤fli, c'est-Ã -dire que pour tout
       sommet le sous-graphe induit par ses voisins est isomorphe au
       graphe de SchlÃ¤fli, qui est lui-mÃªme localement un graphe de
       Clebsch.

....wheel n (=ringarytree 1 0 n 2)
....
       Roue Ã  n rayons. Graphe planaire Ã  n+1 sommets composÃ© d'un
       cycle Ã  n sommets et d'un sommet universel, donc connectÃ© Ã 
       tous les autres.

....web n r (=ringarytree r 1 n 2)
....
       Graphe planaire Ã  1+nÂ·r sommets composÃ© d'une Ã©toile Ã  n
       branches de longueur r, les sommets de mÃªme niveau Ã©tant
       connectÃ©s par un cycle. Il gÃ©nÃ©ralise "wheel n" (r=1).

....binary h (= ringarytree h 2 2 0)
....
       Arbre binaire complet de hauteur h. Il possÃ¨de 2^(h+1)-1
       sommets et la racine est de degrÃ© deux.

....arytree h k r (= ringarytree h k r 0)
....
       Arbre complet de hauteur h oÃ¹ chaque noeud interne Ã  exactement
       k fils, la racine Ã©tant de degrÃ© r.

....rbinary n (= rarytree n 2 0)
....rbinaryz n (= rarytree n 2 1)
....
       Arbre binaire plan alÃ©atoire uniforme Ã  n noeuds internes. Il
       possÃ¨de 2n-1 sommets (2n pour la variante rbinaryz) numÃ©rotÃ©s
       selon un parcours en profondeur modifiÃ©: tous les fils du
       sommet courant sont numÃ©rotÃ©s avant l'Ã©tape de rÃ©cursivitÃ©. La
       racine est de degrÃ© 2 (=rbinary) ou 1 (=rbinaryz). Le dessin
       avec dot (-visu) ne respecte pas le plongement de l'arbre.
       L'option -directed permet d'obtenir une 1-orientation.

....tw n k (= ktree n k -dele .5)
....
       Graphe de largeur arborescente au plus k alÃ©atoire Ã  n
       sommets. Il s'agit d'un k-arbre partiel alÃ©atoire dont la
       probabilitÃ© d'avoir une arÃªte est 1/2. L'option -dele Ã©tant
       dÃ©jÃ  prÃ©sente, il n'est pas conseillÃ© de la rÃ©utiliser pour ce
       graphe. L'option -directed permet d'obtenir une k-orientation.

....pw n k (= kpath n k -dele .5)
....
       Graphe de pathwidth au plus k, alÃ©atoire et avec n sommets.

....tadpole n p (= barbell -n 1 p)
....dragon n p (= barbell -n 1 p)
....
       Graphe Ã  n+p sommets composÃ© d'un cycle Ã  n sommets reliÃ© Ã  un
       chemin Ã  p sommets.

....lollipop n p (= barbell n p 0)
....
       Graphe "tapette Ã  mouches" (Lollipop Graph) composÃ© d'une
       clique Ã  n sommets reliÃ©e Ã  un chemin de longueur p. Il a n+p
       sommets.

....pan n (= barbell -n 1 1)
....
       Graphe Ã  n+1 sommets composÃ© d'un cycle Ã  n sommets et d'un
       seul sommet pendant.

....banner (= barbell -4 1 1)
....
       Graphe Ã  5 sommets composÃ© d'un carrÃ© et d'un sommet pendant.

....paw (= barbell -3 1 1)
....
       Graphe Ã  4 sommets composÃ© d'un triangle et d'un sommet
       pendant.

....theta0 (=barbell -5 -5 -2)
....
       Graphe Theta_0. C'est un graphe Ã  7 sommets sÃ©rie-parallÃ¨le
       obtenu Ã  partir d'un cycle de longueur 6 et en connectant deux
       sommets antipodaux par un chemin de longueur 2. C'est un graphe
       allumette, c'est-Ã -dire planaire et distance unitaire.

....nng n (= knng n 1)
....
       Graphe du plus proche voisin (Nearest Neighbor Graph). Graphe
       gÃ©omÃ©trique dÃ©fini Ã  partir d'un ensemble de n points du carrÃ©
       [0,1[Â² (distribution par dÃ©faut). Le point i est connectÃ© au
       plus proche autre point (par dÃ©faut selon la norme L2, voir
       -norm). Ce graphe est une forÃªt couvrante du graphe rng de
       degrÃ© au plus 6 (si la norme est L2).

....td-delaunay n (= thetagone n 3 3 1)
....
       Triangulation de Delaunay utilisant la distance triangulaire
       (TD=Triangular Distance). Ce n'est malheureusement pas toujours
       une triangulation, les sommets du bord pouvant Ãªtre de degrÃ©
       un. Il s'agit d'un graphe planaire dÃ©fini Ã  partir d'un
       ensemble de n points alÃ©atoires du carrÃ© [0,1[Â² (distribution
       par dÃ©faut). Ce graphe a un Ã©tirement de 2 par rapport Ã  la
       distance euclidienne entre deux sommets du graphe. Ce graphe,
       introduit par Chew en 1986, est le mÃªme que le graphe
       "demi-theta_6", qui est un "theta-graph" utilisant 3 des 6
       cÃ´nes. La dissymÃ©trie qui peut apparaÃ®tre entre le bord droit
       et gauche du dessin est liÃ© au fait que chaque sommet n'a
       qu'une seule bissectrice de cÃ´ne dirigÃ©e vers la droite, alors
       qu'il y en a deux obliques vers la gauche.

....theta n k (= thetagone n 3 k 6/k)
....
       Theta-graphe Ã  k>0 secteurs rÃ©guliers dÃ©fini Ã  partir d'un
       ensemble de n points du carrÃ© [0,1[Â². Les sommets u et v sont
       adjacents si le projetÃ© de v sur la bissectrice de son secteur
       est le sommet le plus proche de u. Ce graphe n'est pas planaire
       en gÃ©nÃ©ral (sauf pour k<3), mais c'est un spanner du graphe
       complet euclidien si kâ‰¥6.

....dtheta n k (= thetagone n 3 âŽ£ k/2âŽ¦ 6/k)
....
       Demi-Theta-graphe Ã  kâ‰¥2 secteurs rÃ©guliers dÃ©fini Ã  partir d'un
       ensemble de n points du carrÃ© [0,1[Â² (distribution par
       dÃ©faut). La dÃ©finition est similaire au Theta-graphe exceptÃ©
       que seul 1 secteur sur 2 est considÃ©rÃ©. Il faut k pair. Pour
       k=2, il s'agit d'un arbre, pour k=4, le graphe est de faible
       tree-width pas toujours connexe.  Pour k=6, ce graphe coÃ¯ncide
       avec le graphe td-delaunay.
....
       Ex: gengraph dtheta 500 6 -visu
           gengraph dtheta 500 4 -pos 0 -visu
           gengraph dtheta 500 2 -pos 0 -visu

....yao n k (= thetagone n 0 k 2/k)
....
       Graphe de Yao Ã  k>0 secteurs rÃ©guliers dÃ©fini Ã  partir d'un
       ensemble de n points du carrÃ© [0,1[Â² (distribution par
       dÃ©faut). Les sommets u et v sont adjacents si v est le sommet
       le plus proche de u (selon la distance euclidienne) de son
       secteur. Ce graphe n'est pas planaire en gÃ©nÃ©ral, mais c'est un
       spanner du graphe complet euclidien. Le rÃ©sultat est valide
       seulement si kâ‰¥2. En fait, ce n'est pas tout Ã  fait le graphe
       de Yao (voir thetagone).

....percolation a b p (= udg aÂ·b 1 -norm L1 -xy mesh a b -dele 1-p)
....
       Grille de percolation Ã  coordonnÃ©es entiÃ¨res (i,j) de [0,a[ Ã—
       [0,b[ oÃ¹ p reprÃ©sente la probabilitÃ© d'existence de chaque
       arÃªte. La diffÃ©rence avec le graphe "mesh a b -dele 1-p" est
       qu'ici le graphe est gÃ©omÃ©trique, donc dessinÃ© selon une grille
       si l'option -visu est ajoutÃ©e.

....hudg n r (= udg n r -norm hyper -xy hyper r)
....
       Graphe gÃ©omÃ©trique alÃ©atoire hyperbolique sur n points du
       disque unitÃ©. Deux sommets sont adjacents si leurs points sont
       Ã  distance hyperbolique â‰¤ r. [A FINIR]

....point n (= ring n . -pos 1)
....
       Graphe gÃ©omÃ©trique composÃ© de n points du plan mais sans aucune
       arÃªte (voir "stable"). Ce graphe permet de visualiser la
       distribution des points, par dÃ©faut uniforme sur [0,1[Â², mais
       qui peut Ãªtre modifiÃ©e avec l'option -xy.
....
       Ex: gengraph point 500 -xy seed 3 2.1 -visu
           gengraph point 1000 -xy seed 3 -0.1 -visu
           gengraph point 1000 -xy disk -visu

....star-polygon n (= ring n 1 . -xy disk)
....
       Polygone "star-shaped" alÃ©atoire Ã  n cotÃ©s contenu dans le
       disque unitÃ© (voir -xy ratio).

....convex-polygon n (= ring n 1 . -xy convex)
....
       Polygone convexe alÃ©atoire Ã  n cotÃ©s contenu dans le disque
       unitÃ© (voir -xy ratio). Voir aussi le graphe polygon n.


....regular n d (= fdrg n d .)
....
       Graphe d-rÃ©gulier alÃ©atoire Ã  n sommets asymptotiquement
       uniforme. Il faut que nd soit pair. L'algorithme est de
       complexitÃ© O(ndÂ²) et pour Ãªtre asymptotiquement uniforme, il
       faut d=o(âˆšn). On obtient nÃ©cessairement un matching si d=1, un
       stable si d=0, un cycle si d=2 et n<6.

....cubic n (= fdrg n 3 .)
....
       Graphe cubique alÃ©atoire Ã  n sommets asymptotiquement
       uniforme. Il faut que n soit pair.

....plrg n t (= bdrg n_1 d_1 â€¦ n_k d_k .)
....
       Power-Law Random Graph (ou scale-free graph). Graphe alÃ©atoire
       Ã  n sommets dont la distribution des degrÃ©s suit une loi en
       puissance d'exposant t>0 (typiquement un rÃ©el tâˆˆ]2,3[), la
       probabilitÃ© qu'un sommet soit de degrÃ© i>0 Ã©tant
       proportionnelle Ã  1/i^t. Plus prÃ©cisÃ©ment, la distribution est
       la suivante:
....
         â€¢ d_1=1, n_1=âŽ£ exp(ð›¼)âŽ¦ + n-s(ð›¼)
         â€¢ d_i=i, n_i=âŽ£ exp(ð›¼)/i^tâŽ¦ pour 2â‰¤iâ‰¤p(ð›¼)
....
       oÃ¹ a est un rÃ©el minimisant |n-s(ð›¼)| avec p(ð›¼)=âŽ£ exp(ð›¼/t)âŽ¦ et
       s(ð›¼)=âˆ‘_{i=1}^{p(ð›¼)} âŽ£ exp(ð›¼)/i^tâŽ¦. Ce sont les mÃªmes graphes
       que ceux gÃ©nÃ©rÃ©s par Brady-Cowen'06 ou ceux Ã©tudiÃ©s par Lu'01.

....syracuse n (= collatz n 1 0 3 1 .)
....
       Graphe orientÃ© issu du problÃ¨me "3x+1" de Collatz. D'aprÃ¨s la
       conjecture, il s'agit d'un graphe connexe. Les arcs sont
       dÃ©finis par la relation xâ†¦x/2 si x est pair, et xâ†¦(3x+1)/2
       sinon. La boule de volume n est gÃ©nÃ©rÃ©e depuis 1 si n>0 en
       itÃ©rant la relation inverse, et dans ce cas le graphe a
       prÃ©cisÃ©ment n sommets. Si n<0, la relation est itÃ©rÃ©e pour tous
       les entiers de [1,|n|] dans la limite de nÂ² sommets.
....
       Ex: gengraph syracuse 18 -label 1 -visu 
           gengraph syracuse -18 -label 1 -visu 
....
       Le rayon de la boule de volume n (excluant les valeurs
       multiples de 3 qui n'ont toujours qu'un seul prÃ©dÃ©cesseur) est
       expÃ©rimentalement < ln(n)/ln(1.36) pour tout n<10^18 et
       conjecturÃ©e en ln(n)/ln(4/3). La hauteur maximum atteinte Ã 
       partir d'un entier de [1,n] est expÃ©rimentalement < nÂ², avec
       seulement 7 exceptions pour n<10^18 (et dans ces cas lÃ , la
       hauteur est < 9nÂ²).  Pour n=27, on obtient une hauteur record
       de 4,616 â‰ƒ 7nÂ², le second record pour n<10^18. Selon les
       conjectures, les trajectoires de hauteur maximum pour n ont une
       forme gÃ©nÃ©rale qui consiste Ã  une montÃ©e â‰ƒ 8Â·ln(n) Ã©tapes pour
       atteindre le maximum, puis une descente â‰ƒ 24Â·ln(n) Ã©tapes, et
       les trajectoires extrÃªmes ont â‰ƒ 42Â·ln(n) Ã©tapes.

....kakutami_3x+1 n (= collatz n 1 0 6 2 .)
....
       Variante de syracuse "non-compressÃ©e" qui est dÃ©finie par la
       relation xâ†¦x/2 si x est pair et et xâ†¦3x+1 sinon. Comme le
       graphe de Collatz, il est possible d'avoir n<0.

....kakutami_5x+1 n (= collatz n 3 0 30 6 3 0 2 0 3 0 30 6 .)
....
       Graphe orientÃ© issu du problÃ¨me "5x+1" de Kakutami
       (cf. syracuse) dÃ©fini par la relation xâ†¦x/2 si x est pair,
       xâ†¦x/3 si x est divisible par 3 mais pas par 2, et xâ†¦5x+1
       sinon. D'aprÃ¨s la conjecture de Kakutami, il s'agit d'un graphe
       connexe. On peut se ramener Ã  un graphe de Collatz en
       considÃ©rant 2Ã—3=6 paires de coefficients. Comme le graphe de
       Collatz, il est possible d'avoir n<0.
....
       Ex: gengraph kakutami_5x+1 91 -undirected -visu

....kakutami_7x+1 n (= collatz n 15 0 210 30 15 0 10 0 â€¦ 210 30 .)
....
       Graphe orientÃ© issu du problÃ¨me "7x+1" de Kakutami
       (cf. syracuse) dÃ©fini par la relation xâ†¦x/2 si x est pair,
       xâ†¦x/3 si x est divisible par 3 mais pas par 2, xâ†¦x/5 si x est
       divisible par 5 mais ni par 2 ni par 3, et xâ†¦7x+1 sinon.
       D'aprÃ¨s la conjecture de Kakutami, il s'agit d'un graphe
       connexe. On peut se ramener Ã  un graphe de Collatz en
       considÃ©rant 2Ã—3Ã—5=30 paires de coefficients. Comme le graphe de
       Collatz, il est possible d'avoir n<0.

....farkas n (= collatz n 6 0 6 6 6 0 4 0 6 0 â€¦ 18 6 .)
....
       Graphe orientÃ© issu de la relation introduite par Farkas en
       2005 qui a prouvÃ© qu'elle dÃ©finissait un arbre de racine 1.
       Proche du problÃ¨me "3x+1" elle est dÃ©finie par xâ†¦x/2 si x est
       pair, xâ†¦x/3 si x est divisible par 3 mais pas par 2, (3x+1)/2
       si xâ‰¡3 (mod 4) et (3x+1)/2 si xâ‰¡1 (mod 4). On peut se ramener Ã 
       un graphe de Collatz en considÃ©rant 12 paires de coefficients.
       Comme le graphe de Collatz, il est possible d'avoir n<0.
....
       Ex: gengraph farkas -34 -label 1 -visu

.....
HISTORIQUE

       v1.2 octobre 2007:
            - premiÃ¨re version

       v1.3 octobre 2008:
            - options: -shift, -width
            - correction d'un bug pour les graphes de permutation
	    - accÃ©lÃ©ration du test d'ajacence pour les arbres, de O(n) Ã  O(1),
              grÃ¢ce Ã  la reprÃ©sentation implicite
	    - nouveau graphes: outerplanar, sat

       v1.4 novembre 2008:
            - format de sortie: matrix, smatrix, list
            - nouveau graphe: kout
            - correction d'un bug dans l'option -width
	    - correction d'un bug dans la combinaison -format/shift/delv

       v1.5 dÃ©cembre 2008:
            - correction d'un bug dans tree lorsque n=1

       v1.6 dÃ©cembre 2009:
            - nouveaux graphes: rpartite, bipartite

       v1.7 janvier 2010:
            - nouveaux graphes: icosa, dodeca, rdodeca, cubocta, geo,
	      wheel, cage, headwood, pappus, mcgee, levi, butterfly,
	      hexagon, whexagone, arytree, binary, ktree, tw, kpath,
	      pw, arboricity, wagner, mobius, tutte-coexter, paley
            - nouveau format de sortie: -format dot
	    - nouvelles options: -header, -h, -redirect, -dotpdf
            - correction d'un bug dans kout, et dans tree lorsque n=0
	    - tree devient un cas particulier d'arboricity.
	    - aide en ligne pour les paramÃ¨tres des graphes.

       v1.8 juillet 2010:
            - nouveaux graphes: chvatal, grotzsch, debruijn, kautz
	      gpstar, pstar, pancake, nauru, star, udg, gpetersen,
              mobius-kantor, desargues, durer, prism, franklin,
	      gabriel, thetagone, td-delaunay, yao, theta, dtheta
            - suppression du graphe geo (remplacÃ© par udg)
            - nouvelles options: -pos, -norm, -label, -dotfilter
	    - nouvelle famille d'options: -xy file/noise/scale/seed
	    - dÃ©finition plus compacte dodeca (non explicite)
	    - utilisation du gÃ©nÃ©rateur random() plutÃ´t que rand().
	    - correction d'un bug dans "-format standard" qui provoquait une erreur.
	    - correction d'un bug dans kneser pour k=0, n=0 ou k>n/2.
	    - nouveaux formats: -format dot<type>, -format xy
	    - suppression de -dotpdf (qui est maintenant: -format dotpdf)
	    - labeling pour: gpetersen, gpstar, pstar, pancake, interval,
	      permutation

       v1.9 aoÃ»t 2010:
            - renome -h en -list
	    - renome -xy file en -xy load
	    - centrage des positions sur le barycentre des graines (-xy seed)
	    - nouvelles options: -star, -visu, -xy round
	    - les graphes peuvent Ãªtre stockÃ©s en mÃ©moire, sous la forme d'une liste
	      d'adjacence grÃ¢ce Ã  l'option -check.
	    - gÃ©nÃ©ralisation de -delv p avec p<0
	    - nouveaux graphes: caterpillar, hajos, hanoi, sierpinski, sunlet, load
	    - labeling pour: hanoi, sierpinski
	    - aide sur toutes les options (nÃ©cessitant au moins un paramÃ¨tre)
              et non plus seulement pour les graphes
	    - nouvelle famille d'options: -vcolor deg/degr/pal
	    - correction d'un bug pour l'aide dans le cas de commande
	      prÃ©fixe (ex: pal & paley)

       v2.0 septembre 2010:
	    - nouvelles options: -vcolor degm/list/randg, -xy unique/permutation,
	      -check bfs, -algo iso/sub
	    - l'option -xy round p admet des valeurs nÃ©gatives pour p.
	    - les options "load file" et "-xy load file" permettent la
              lecture Ã  partir de l'entrÃ©e standard en mettant
              file="-", la lecture de famille de graphes, et supporte les commentaires.
	    - les formats list/matrix/smatrix utilisent un espace
	      linÃ©aire O(n+m) contre O(nÂ²) auparavant.
	    - les sommets sur le bord (graphes gÃ©omÃ©triques) ne sont plus coupÃ©s
	      (bounding-box (bb) plus grandes).
	    - nouveaux graphes: kpage, outerplanar n (=kpage n 1), rng, nng
	      fritsch, soifer, gray, hajos (qui avait Ã©tÃ© dÃ©finit mais non
	      implÃ©mentÃ© !), crown, moser, tietze, flower_snark, markstrom,
	      clebsch, robertson, kittell, rarytree, rbinary, poussin, errera
	    - les graphes de gabriel (et rng,nng) dÃ©pendent maintenant de -norm.
	    - "wheel n" a maintenant n+1 sommets, et non plus n.
	    - aide en ligne amÃ©liorÃ©e avec "?". Ex: gengraph tutte ? / -visu ?
	    - les options -help et ? permettent la recherche d'un mot clef.
	      Ex: gengraph -help planaire / ? arbre
	    - description plus compacte de tutte (et des graphes Ã  partir d'un tableau)
	    - correction d'un bug pour rpartite (qui ne marchait pas)

       v2.1 octobre 2010:
	    - nouvelles options:
	      -check degenerate/gcolor/edge/dfs/ps1/paths/paths2/iso/sub/minor/isub
	      -filter minor/sub/iso/vertex/edge/degenerate/ps1
	      -filter degmax/degmin/deg/gcolor/component/radius/girth/diameter
	      -filter cut-vertex/biconnected/isub/all/minor-inv/isub-inv/sub-inv
            - suppression de -algo iso/sub: l'option -algo est rÃ©servÃ©e Ã  la mis
	      au point de -check
	    - extension de -label b Ã  b=2 qui force l'affiche des noms
              sous forme d'entiers mÃªme avec -permute.
	    - correction d'un bug pour house (qui ne marchait pas)
	    - nouveau graphe: windmill

       v2.2 novembre 2010:
            - gestion des graphes orientÃ©s: lecture d'un fichier de
              graphe (ou d'une famille avec arcs et arÃªtes)
	    - nouvelles options: -(un)directed, -(no)loop, -check twdeg/tw,
	      -filter tw/id/hyperbol/rename
	    - permet l'affichage de la "value" (p) dans l'option -filter
	    - nouveau graphe: aqua
	    - correction du graphe tutte-coexter et suppression du
              graphe levi (qui en fait Ã©tait le graphe de tutte-coexter).
	    - gÃ©nÃ©ralisation de l'option "load" Ã  load:id family

       v2.3 dÃ©cembre 2010:
            - nouvelles options: -check ps1bis/edge, -filter ps1bis/tw2
	      -filter minus/minus-id/unique/connected/bipartite/forest
	      -check ps1ter
	    - remplacement de atof()/atoi() par strtod()/strtol() qui
	      sont plus standards.
	    - remplacement de LONG_MAX par RAND_MAX (=2^31-1) dans les
              expressions faisant intervenir random() qui est de type
              long mais qui est toujours dans [0,2^31[, mÃªme si
              sizeof(long)>4. Il y avait un bug pour les architectures
              avec sizeof(long)=8.
	    - nouveau graphe: cylinder
	    - suppression de la variante "load:id" au profit de la
              forme plus gÃ©nÃ©rale "file:range" valable pour load, -filter, etc.

       v2.4 janvier 2011:
            - correction d'un bug dans -filter minus-id
	    - correction d'un bug dans rpartite (incorrect Ã  partir de r>5 parts)
	    - correction d'un bug dans whexagon (nb de sommets incorrects)
	    - nouvelles options: -check ps1x/girth, -filter ps1c/ps1x
	    - renomage: ps1bis -> ps1b, ps1ter -> ps1c
	    - nouveau graphe: mycielski
	    - la graphe grotzsch est maintenant dÃ©fini Ã  partir du graphe
	      mycielski (la dÃ©finition prÃ©cÃ©dante Ã©tait fausse)
	    - bug dÃ©tectÃ©: td-delaunay 500 -check gcolor -format no -seed
              7 | grep '>6' qui donne jusqu'Ã  7 couleurs; le nb de
              couleurs affichÃ©es dans -check gcolor est erronÃ©

       v2.5 mars 2011:
	    - nouveaux graphes: line-graph, claw

       v2.6 juin 2011:
	    - amÃ©lioration du test -filter ps1: dÃ©tection de cliques et d'arbres

       v2.7 octobre 2011:
	    - nouvelle option: -check bellman (pour les gÃ©omÃ©triques seulement)
	    - ajout des champs xpos,ypos Ã  la structure "graph".
	    - nouveaux graphes: linial, linialc, cube, diamond, theta0,

       v2.8 novembre 2011:
	    - nouveaux graphes: ggosset, gosset, rplg, wiener-araya, headwood4
	    - correction d'un bug pour "-xy seed k n" lorsque k=1.
	    - nouvelles options: -check maincc, -maincc (non documentÃ©e)

       v2.9 fÃ©vrier 2013:
	    - nouveau graphe: frucht, halin
	    - correction d'un bug pour "-check gcolor" qui ne
	      renvoyait pas le nombre correct de couleurs, et qui de
	      plus n'utilisait pas l'heuristique du degrÃ© minimum.
	    - correction d'un bug pour "permutation -label 1"

       v3.0 octobre 2013:
	    - nouveaux graphes: rig, barbell, lollipop
	    - gÃ©nÃ©ralisation de l'option -filter forest
	    - nouvelles options: -apex, -filter isforest, -filter istree, -filter cycle
	    - correction d'un bug dans -filter vertex
	    - amÃ©lioration de l'aide lors d'erreurs de paramÃ¨tre comme:
              "-filter F vertex" au lieu de "-filter F vertex n"
	    - amÃ©lioration de l'option -header

       v3.1 dÃ©cembre 2013:
	    - nouveaux graphes: bpancake
	    - lÃ©gÃ¨re modification des labels des sommets des graphes pancake, 
	      gpstar et pstar
	    - nouvelles options: -xy grid, -xy vsize
	    - modification de la taille des sommets pour dot permettant de tenir
	      compte de -xy scale.

       v3.2 mai 2014:
            - amÃ©lioration du test ps1b (ajoÃ»t de rÃ¨gle et rÃ©duction
	      du nombre d'indÃ©terminÃ©es dans graphes des conflits)

       v3.3 juillet 2014:
            - modification importante du code pour -check ps1
            - modification des graphes linial et linialc
	    - nouvelles options: -check kcolor, -vcolor kcolor, -len, -check kcolorsat

       v3.4 fÃ©vrier 2015:
            - documentation et mise au point de l'option -maincc
	    - correction d'un bug lors de la combinaison de "load file" et de "-vcolor pal grad"
	    - correction d'un bug dans la fonction SortInt() qui affectait "-vcolor deg"
	    - correction d'un bug avec l'option -label 1 pour certains graphes (outerplanar ...)
            - crÃ©ation du script dot2gen.awk pour convertir le format dot en format standard
	    - nouvelles options: -fast, -caption
	    - introduction du groupement d'arÃªtes/arcs i-(j k ...) dans le format standard
	      (il reste un bug si ce format est lu depuis l'entrÃ©e standard)

       v3.5 mars 2015:
	    - correction d'un bug pour le dÃ©gradÃ© de couleurs avec "-vcolor pal grad"
	    - nouvelles options: -check ncc/connected, -check routing scenario [...] cluster
	    - nouvelle variante de graphe: loadc file
	    - nouveaux graphes: octahedron, turan, hexahedron, tetrahedron, deltohedron,
	      trapezohedron, antiprism, flip, associahedron, shuffle
	    - changement de nom (ajoÃ»t du suffixe "hedron") pour: isoca, dodeca, cubocta,
	      rdodeca

       v3.6 juin 2015:
            - nouvelles options: -version, -variant, -check info, -xy mesh,
	      -check routing tzrplg, -check routing dcr
	    - nouveaux graphes: percolation, hgraph, cricket, moth, bull, expander
	    - correction de bug dans "-help ?", dans "-check girth" qui donnait -1 pour
	      un cycle de taille n, dans butterfly (mauvais nombre de paramÃ¨tres)
	    - vÃ©rification que le nombre de sommets n'est pas nÃ©gatif
	      pour plusieurs graphes dont tree, kout, etc. ce qui
	      pouvait provoquer des erreurs
	    - vÃ©rification du format d'entrÃ©e pour load/loadc
	    - l'option -check maincc n'affiche plus le graphe initial
	    - description du format des familles de graphes dans l'aide
	    - affiche plus de messages d'erreurs: lecture d'un graphe au mauvais
	      format, mauvaise combinaison de d'options, ...

       v3.7 juillet 2015:
            - nouvelle implÃ©mentation des mots de Dyck, qui sont
              maintenant uniformes. En consÃ©quence les graphes
              alÃ©atoires tree, arboricity, rarytree, rbinary, outerplanar,
              kpage â€¦ sont gÃ©nÃ©rÃ©s de maniÃ¨re uniformes.
	    - nouveaux graphes: treep (et simplification de halin),
	      ygraph, ringarytree (et simplification de arytree,
	      binary, wheel), netgraph, egraph, rgraph (=fish),
	      gÃ©nÃ©ralisation de barbell, tadpole (=dragon), pan,
	      banner, paw, theta0 (redÃ©finition), fan, gem, chess,
	      knight, antelope, camel, giraffe, zebra, utility,
	      zamfirescu, hatzel, web, bdrg, fdrg, regular, cubic, plrg
	    - nouvelles options: -check radius, -check diameter
	    - correction d'un bug si combinaison d'options -check ncc -format list
	      (mauvaise gestion de la variable CHECK)
	    - introduction de caractÃ¨res UTF8 mathÃ©matiques dans l'aide

       v3.8 novembre 2015:
            - nouveaux graphes: ladder, matching, d-octahedron, johnson
	    - correction d'un bug pour le graphe kpage (introduit Ã  v3.7)

       v3.9 fÃ©vrier 2016:
	    - renomage de l'option "-xy scale" en "-xy box"
	    - correction d'un bug dans rbinary (mauvais alias/dÃ©finition)
	    - correction d'un bug (introduit Ã  v3.6) dans la fonction
              max(double,double) au lieu de fmax() qui affectait certains
              graphes et options gÃ©omÃ©triques (rng, -xy grid, percolation, ...)
            - correction d'un bug concernant rarytree lorsque b>2
	    - correction d'un bug dans les statistiques pour -check routing
	    - gÃ©nÃ©ralisation du graphe rarytree (augmentation du degrÃ© de la racine)
	    - nouveaux graphes: herschel, goldner-harary, rbinaryz,
	      point, empty, apollonian, dyck, cross, tgraph, bidiakis, gear,
	      centipede, harborth, sunflower, planar, squaregraph, polygon
	    - modification de certains graphes (cage, ring, gabriel,
              nng, rng) afin que le test d'adjacence ne dÃ©pende plus de N
              mais de PARAM[0] par exemple
	    - introduction du format %SEED pour l'option -caption

       v4.0 mars 2016:
	    - nouveaux graphes: pat, star-polygon, convex-polygon
	    - nouvelles options: -check kindepsat, -xy circle, -xy polar,
	      -xy unif, -xy convex, -xy ratio, -xy zero
	    - aide permettant les prÃ©fixes comme dans -check routing cluster

       v4.1 avril 2016:
	    - nouvelles options: -xy convex2, -check routing hdlbr/bc
	    - nouveaux graphes: kstar, split, cactus, squashed
	    - suppression de -check paths2 qui donnait toujours comme -check paths
	    - option -check bellman apparaÃ®t dans l'aide

       v4.2 mai 2016:
	    - nouvelles options: -check routing agmnt, -format vertex,
	      -check routing hash mix, -xy unif, -xy polygon
	    - amÃ©lioration de l'option -check info
	    - prise en compte de -xy ratio dans -xy unif et -xy seed
	    - affichage de l'Ã©cart type dans les stats affichÃ© par -check routing
	    - affichage des erreurs sur stderr plutÃ´t que stdout
	    - nouveau graphe: suzuki

       v4.3 juin 2016:
	    - nouveaux graphes: triplex, jaws, starfish
	    - suppression d'un bug (Abort trap: 6) pour -format dot<type>
	    - modification de l'initialisation du gÃ©nÃ©rateur alÃ©atoire qui pouvait
              faire comme -seed 0 sur certains systÃ¨mes oÃ¹ clock() est toujours nulle
	    - ajoÃ»t de variantes pour -check routing cluster
	    - nouvelles options: -xy cycle, -xy disk, -check stretch, -xy surface
	    - modification de l'option -norm (-norm 2 -> -norm L2, etc.)
	    - renomage de -xy polar en -xy disk

       v4.4 juillet 2016:
	    - nouvelle option: -check simplify
	    - nouveaux graphes: schlafli, doily, fork (=tgraph)
	    - suppression d'un bug pour gosset (mauvais paramÃ©trage)

       v4.5 aoÃ»t 2016:
	    - modification du terminateur de sÃ©quence pour bdrg et fdrg: -1 devient '.'
	    - refonte du prototype des fonctions d'adjacence: adj(i,j) -> adj(Q)
	    - dÃ©sactivation des options -apex, -star (et donc de caterpillar)
	    - finalisation de l'implÃ©mentation de fdrg (et donc de regular et cubic)
	    - redÃ©finition des graphes: netgraph, egraph, sunlet, cross, tgraph, ygraph
	    - nouveaux graphes: comb, alkane (et ses nombreuses variantes), banana
	      rlt, circle
	    - l'option -fast est effective pour fdrg et bdrg
	    - suppression d'un bug pour -xy mesh
	    - correction du graphe "interval" qui renvoyait en fait un "circle"
	    - gÃ©nÃ©ration d'une aide html Ã  partir du source: gengraph.html
	    - nouvelle option: -norm poly p

       v4.6 septembre 2016:
	    - nouvelle option: -check routing scenario until s
	    - nouveaux graphes: uno, unok, behrend

       v4.7 janvier 2017:
	    - dans -check dfs/bfs, dÃ©tection d'une source invalide
	    - dans -check dfs, calcule et affichage la hauteur des sommets
	    - ajoÃ»t d'un paramÃ¨tre pour unok

       v4.8 mars 2017:
	    - calcul dans -check stretch du stretch max. minimum
	    - correction d'un bug dans -xy convex qui pouvait produire des ensembles
	      non convexes et des points en dehors du cercle de rayon 1
	    - nouveaux graphes: ngon

       v4.9 juin 2017:
	    - correction d'un bug pour l'option -xy box et amÃ©lioration
	      du dessin de la grille pour -xy grid
	    - gÃ©nÃ©ralisation de -check bellman aux graphes non valuÃ©s
            - renomage des options: -dotfilter -> -dot filter,
	      -len -> -dot len, -filter hyperbol -> -filter hyper
	    - nouvelles options: -norm hyper, -xy hyper, -dot scale,
	      -label b pour b=3 et b<0, -check volm
	    - nouveaux graphes: hudg, hyperbolic

       v5.0 aoÃ»t 2017:
	    - nouveaux graphes: helm, haar
	    - nouveau format html

       v5.1 aoÃ»t 2018:
            - redÃ©finition (plus simple) du graphe de TurÃ¡n
	    - redÃ©finition de bidiakis comme graphe cage
	    - redÃ©finition de centipede comme graphe comb (et sans -star)
	    - simplification de bellman avec optimisation pour
              utilisation multiples comme dans -check stretch
	    - amÃ©lioration de dÃ©tails dans l'aide html (blocs "!!!" et "Ex:")
	    - changement des paramÃ¨tres (maintenant avec un . final)
	      pour: grid, aqua, rpartite, cage, ring, ggosset
	    - orientation corrigÃ©e avec -directed pour les graphes (et
              leurs composÃ©s) utilisant une reprÃ©sentation implicite:
              kout, ktree, kpath, kpage, hyperbolic, cactus, planar,
              rarytree, apollonian, expander, polygon
            - orientation supportÃ©e pour cage et ring (possibilitÃ© de cordes<0)
            - nouveaux graphes: margulis, collatz, syracuse, mst, farkas,
	      kakutami_[3|5|7]x+1 (3 graphes), [w|u][psl|dis][d] (8 graphes)
            - redÃ©finition de l'option -loop avec suppression de -noloop
	    - inversion des lignes/colonnes dans uno et unok
	    - redÃ©finition du cuboctahedron Ã  partir de linial
	    - aide sur le graphe chvatal qui Ã©tait abscente
	    - format xy: affichage coordonnÃ©es entiÃ¨res si elles le sont
	    - bug corrigÃ© pour -label b avec b<0: affiche exclusif centrÃ© ou pas
	    - bug corrigÃ© pour -xy zero et -xy bord apparu Ã  la v49
	    - dÃ©finitions des graphes oubliÃ©s mobius-kantor et dragon
	    - utilisation d'une version uniforme pour remplacer random()%k
	    - nouvelle option: -dot scale auto

       v5.2 juin 2019:
            - nouveaux graphes: dart, kite, domino, hourglass, antenna,
	      parachute, parapluie, klein
	    - nouvelles abrÃ©viations pour les graphes alkane et option
	      -label 1 activÃ©e par dÃ©faut, correction d'un bug pour le
	      graphe methane
	    - chargement conditionnel des scripts js dans le format html
	    - nouvelle option: -visuh, -check subdiv
	    - correction d'un bug pour unok, introduit Ã  la version v5.1
	    - correction d'un bug pour -check stretch (Bellman-Ford
              avec appels multiples)
	    - correction d'un bug pour nb_edges() introduit Ã  la
	      version v5.1 et qui affectait par exemple -check edges
	    - correction bug dans la dÃ©tection de sÃ©quence graphique
	      qui pouvait affecter fdrg, regular et cubic

       v5.3 aoÃ»t 2019:
            - nouveaux graphes: knng, rectree
	    - redÃ©finition de nng Ã  partir de knng
	    - renommage de headwood et headwood4 en heawood et heawood
### #*/
