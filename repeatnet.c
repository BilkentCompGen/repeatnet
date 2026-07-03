/*
	RepeatNet: an ab initio centromeric sequence detection algorithm
  
*/




/*

TODO : 

 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>

#define SEQ_LENGTH 150

#define FORWARD 0
#define REVERSE 1

#define LOAD 0
#define CREATE 1

#define ASCENDING 1
#define DESCENDING 0

typedef struct imatrix{
  int i;
  int j;
  int eweight;
}_imatrix;

typedef struct wpair{
  int id;
  struct wpair *next;
}_wpair;


typedef struct Node{
  int v1,v2,w;
  struct Node *next;
}_Node;

typedef struct Component{
  struct Node *node;
  struct Component *next;
}_Component;

void my_fgets(char *, int, FILE *);
int readSingleFasta(FILE *, int);
int readFastaRecord(FILE *, char **, char **);
int readPairedFasta(FILE *, FILE *);
void usage(char *);
int hash(char *, int, int);
void end2clone(char *, char *);
void rcomp(char *, char *);
int isReverse(char *);
int encode_window(char *);
void decode_window(int, char *);
int getWindow(char *, int, char *, int);
void do_end(int, int);
void insert_wloc(struct wpair **, int, int, int);
void compare_wloc(struct wpair *, struct wpair *, int *);
int hashTest(FILE *, int);
void dumpwins(FILE *, int);
void do_loaded(FILE *, char *);
void loadwins(FILE *);
char reverseindex(int);
void compare(char *, int, float, int, int, int);
unsigned int cindex(char);
int encode_window2(char *);
int revcomp_encoded(unsigned int);
void freewins(struct wpair *);
struct wpair *copyWins(struct wpair *);
void compactwins(void);
void merge(int, int);
FILE *myfopen(char *, char *);
int loadMatrix(char *, char *);
void decode_wins(char *);
void loadNames(char *);
void dumpclones(char *);
void trimwins(char [], int, char []);


/* component stuff */
void connectedComponent(int, int);
void alloc_node(struct Node **);
void alloc_component(struct Component **);
void insert_node_to_component(struct Component **, int, int, int);
void makecolors(void);
void connComp(int, int);


char **names;
char **forwards;
char **reverses;
int nseq;
int WS;
int *windowflag;
int wincnt;
struct wpair **winlocs;
int collisions;
int *windowcnt;

char **colors;

struct imatrix *interactions;

int EDGEHIST;
int WINORDER;


int main(int argc, char **argv){
  
  FILE *in = NULL;
  int i,j,k;
  FILE *winlog;
  FILE *namelog;
  int npairs;
  struct timeval start, end;
  struct timezone tz;
  char fname[SEQ_LENGTH];
  char fname2[SEQ_LENGTH];
  char namefile[SEQ_LENGTH];
  char forwardfile[SEQ_LENGTH];
  char reversefile[SEQ_LENGTH];
  FILE *fwdfp = NULL;
  FILE *revfp = NULL;
  int PAIRED = 0;
  char *encodearg = NULL;
  int hashtype;
  int hashtest=0;
  int mintime;
  int minhash;
  float ecut;
  int acut;

  long inserttime;
  long compacttime;
  long readtime;
  long inittime;
  int mincollide=999999999;
  int mincolhash;
  int ccut; // compare cut; do not compare any kmer with cnt<ccut
  int LOADWIN = 0;
  int COMPARE = 0;
  int MERGE = 0;
  int LOADMATRIX = 0;
  int DECODE = 0;
  int DUMPCLONES = 0;
  int LOADNAMES = 0;
  int COMPONENTS = 0;
  int TRIMWIN = 0;
  char trimwinfile[SEQ_LENGTH];

  int matrixcnt;

  EDGEHIST = 0;

  collisions=0;
  
  WS = 12;
  ccut = 1;
  acut = 0;

  if (argc < 2){
    usage(argv[0]);
    return 0;
  }

  fname[0]=0;
  fname2[0]=0;
  namefile[0]=0;
  forwardfile[0]=0;
  reversefile[0]=0;
  hashtype=11;
  mintime=999999999;

  inserttime=0;
  compacttime=0;
  readtime=0;
  inittime=0;
  ecut = 0.0;

  /* long-only option identifiers (values outside the printable-char range) */
  enum {
    OPT_LOADWIN = 256, OPT_LOADNAMES, OPT_LOADMATRIX, OPT_GETN,
    OPT_DECODE, OPT_ENCODE, OPT_COMPONENTS, OPT_COMPARE,
    OPT_EDGEHIST, OPT_CLONES, OPT_TRIM, OPT_HELP
  };

  static struct option long_opts[] = {
    {"input",      required_argument, 0, 'i'},
    {"forward",    required_argument, 0, 'f'},
    {"reverse",    required_argument, 0, 'r'},
    {"hash",       required_argument, 0, 'H'},
    {"hashtest",   required_argument, 0, 't'},
    {"ecut",       required_argument, 0, 'e'},
    {"acut",       required_argument, 0, 'a'},
    {"ccut",       required_argument, 0, 'c'},
    {"merge",      no_argument,       0, 'm'},
    {"loadwin",    required_argument, 0, OPT_LOADWIN},
    {"loadnames",  required_argument, 0, OPT_LOADNAMES},
    {"loadmatrix", required_argument, 0, OPT_LOADMATRIX},
    {"getn",       required_argument, 0, OPT_GETN},
    {"decode",     required_argument, 0, OPT_DECODE},
    {"encode",     required_argument, 0, OPT_ENCODE},
    {"components", no_argument,       0, OPT_COMPONENTS},
    {"compare",    no_argument,       0, OPT_COMPARE},
    {"edgehist",   no_argument,       0, OPT_EDGEHIST},
    {"clones",     required_argument, 0, OPT_CLONES},
    {"trim",       required_argument, 0, OPT_TRIM},
    {"help",       no_argument,       0, OPT_HELP},
    {0, 0, 0, 0}
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "i:f:r:H:t:e:a:c:m",
			    long_opts, NULL)) != -1){
    switch (opt){
    case 'i': strcpy(fname, optarg); break;
    case 'f': strcpy(forwardfile, optarg); break;
    case 'r': strcpy(reversefile, optarg); break;
    case 'H': hashtype = atoi(optarg); break;
    case 't': hashtest = atoi(optarg); break;
    case 'e': ecut = atof(optarg); break;
    case 'a': acut = atoi(optarg); break;
    case 'c': ccut = atof(optarg); break;
    case 'm': MERGE = 1; break;
    case OPT_LOADWIN:    strcpy(fname, optarg); LOADWIN = 1; break;
    case OPT_LOADNAMES:  strcpy(namefile, optarg); LOADNAMES = 1; break;
    case OPT_LOADMATRIX: strcpy(fname, optarg); LOADMATRIX = 1; break;
    case OPT_GETN:       strcpy(fname2, optarg); break;
    case OPT_DECODE:     strcpy(fname, optarg); DECODE = 1; break;
    case OPT_ENCODE:     encodearg = optarg; break;
    case OPT_COMPONENTS: COMPONENTS = 1; break;
    case OPT_COMPARE:    COMPARE = 1; break;
    case OPT_EDGEHIST:   EDGEHIST = 1; break;
    case OPT_CLONES:     DUMPCLONES = 1; strcpy(fname2, optarg); break;
    case OPT_TRIM:       TRIMWIN = 1; strcpy(trimwinfile, optarg); break;
    case OPT_HELP:       usage(argv[0]); return 0;
    case '?':            usage(argv[0]); return 1;
    }
  }

  /* --encode just converts a k-mer string to its integer code and exits. */
  if (encodearg != NULL){
    printf("%d\n", encode_window(encodearg));
    return 1;
  }

  /* Two-file (separate forward/reverse) input requires both files. */
  if (forwardfile[0] || reversefile[0]){
    if (!forwardfile[0] || !reversefile[0]){
      fprintf(stderr, "Both --forward and --reverse are required for two-file input.\n");
      return 0;
    }
    PAIRED = 1;
  }

  if (PAIRED){
    fwdfp = myfopen(forwardfile, "r");
    revfp = myfopen(reversefile, "r");
    /* base the output file names on the forward reads file */
    strcpy(fname, forwardfile);
  }
  else{
    if (fname[0] == 0){
      fprintf(stderr, "input file?\n");
      return 0;
    }
    in = myfopen(fname, "r");
  }


  if (LOADNAMES){
    loadNames(namefile);
  }

  if (DECODE){
    decode_wins(fname);
    return 1;
  }
  
  if (LOADMATRIX){
    matrixcnt = loadMatrix(fname, fname2);
    if (COMPONENTS)
      connComp(matrixcnt, wincnt);
    return 1;
  }

  if (LOADWIN){
    do_loaded(in, fname);
    if (COMPARE)
      compare(fname, hashtype, ecut, ccut, acut, MERGE);
    if (DUMPCLONES)
      dumpclones(fname2);
    if (TRIMWIN)
      trimwins(fname, hashtype, trimwinfile);
    return 1;
  }

  if (hashtest){
    for (i=3;i<hashtest;i++){
      rewind(in);
      for (j=0;j<nseq;j++)
	free(names[j]);
      free(names);
      gettimeofday(&start, &tz);
      collisions = 0;
      nseq = hashTest(in, i);
      gettimeofday(&end, &tz);
      fprintf(stdout, "hash type %d read time %d microseconds.\tCollisions:%d\n", i, ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec)), collisions);
      if (((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec))<mintime){
	mintime = ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec));
	minhash=i;
      }
      if (collisions<mincollide){
	mincollide=collisions;
	mincolhash=i;
      }
    }
    fprintf(stdout, "\n\nmintime: %d microseconds\nbesthash:%d\n",mintime, minhash);
    fprintf(stdout, "\n\nmincollide: %d times\nbesthash:%d\n",mincollide, mincolhash);
    return 0;
  }

  gettimeofday(&start, &tz);
  if (PAIRED)
    nseq = readPairedFasta(fwdfp, revfp);
  else
    nseq = readSingleFasta(in, hashtype);
  gettimeofday(&end, &tz);
  readtime = ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec));
	
  npairs = 0;
  /* compact seqs */

  gettimeofday(&start, &tz);
  sprintf(fname2, "%s.h%d.names", fname, hashtype);
  namelog = fopen(fname2, "w");
  fprintf(namelog, "%d\n",nseq);
  for (i=0;i<nseq;i++){
    fprintf(stderr, "\rCompacting %d\tof\t%d",(i+1),nseq);
    if (names[i]!=NULL){
      if (forwards[i]==NULL){
	free(reverses[i]);
	free(names[i]);
	names[i] = NULL;
      }
      else if (reverses[i]==NULL){
	free(forwards[i]);
	free(names[i]);
	names[i] = NULL;
      }
      else{
	npairs++;
	fprintf(namelog, "%d\t%s\n", i, names[i]);
      }
    }
  }
  fclose(namelog);
  fprintf(stderr, "\n");
  gettimeofday(&end, &tz);
  compacttime = ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec));
  
  
  fprintf(stderr, "Remaining pairs: %d\n", npairs);
  wincnt = pow(4, WS);
  
  fprintf(stdout, "Number of windows:%d\n", wincnt);

  winlocs = (struct wpair **)malloc(sizeof(struct wpair *)*wincnt);
  windowflag = (int *) malloc(sizeof(int)*wincnt);
  windowcnt = (int *) malloc(sizeof(int)*wincnt);
  
  gettimeofday(&start, &tz);  
  fprintf(stdout, "Initializing arrays...");
  for (i=0;i<wincnt;i++){
    winlocs[i]=NULL;
    windowflag[i] = -1;
    windowcnt[i] = 0;
  }
  fprintf(stdout, "\n");
  fprintf(stderr, "\n");
  gettimeofday(&end, &tz);
  inittime = ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec));

  k=0;
  for (i=0; i<nseq; i++){
    if (names[i] == NULL)
      continue;
    
    // else, both ends are here.
    k++;
    fprintf(stderr, "\rInserting %d\tof\t%d",k,npairs);
    
    gettimeofday(&start, &tz);  
    do_end(i, FORWARD);
    do_end(i, REVERSE);
    gettimeofday(&end, &tz);
    inserttime += ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec));
    

    // remove index i, I don't need them anymore.
    // but keep the names for now. 
    free(forwards[i]);
    free(reverses[i]);
    forwards[i]=NULL;
    reverses[i]=NULL;

  }
  fprintf(stderr, "\n\n");
  
  sprintf(fname2, "%s.h%d.winlog", fname, hashtype);
  winlog = fopen(fname2, "w");
  

  fprintf(stderr, "Comparing...");
  fflush(stderr);

  fprintf(winlog,"%d\n",WS);

  sprintf(fname2, "%s.h%d.dump", fname, hashtype);
  FILE *dumphist = fopen(fname2, "w");


  for(i=0;i<wincnt;i++){
    if (winlocs[i]==NULL)
      continue;
    dumpwins(winlog, i);
    fprintf(dumphist, "%d\t%d\n", i, windowcnt[i]);
  }
  fclose(winlog);
  fclose(dumphist);
  fprintf(stderr, "\n");

  fprintf(stdout,"\n-------------------------------------\n");
  fprintf(stdout,"\tread time:\t%ld usec\n", readtime);
  fprintf(stdout,"\tcompact time:\t%ld usec\n", compacttime);
  fprintf(stdout,"\tinit time:\t%ld usec\n", inittime);
  fprintf(stdout,"\tinsert time:\t%ld usec\n", inserttime);
  fprintf(stdout,"-------------------------------------\n");

  return 1;

}


int hashTest(FILE *fastaFile, int hashtype){
  int cnt;
  char myclone[SEQ_LENGTH];
  int index;
  char dummy[SEQ_LENGTH];
  int i;

  cnt = 0;
  fprintf(stderr, "Counting sequences.\n");

  while(fscanf(fastaFile, "%s", dummy)>0)
    cnt++;

  names = (char **) calloc(cnt, sizeof(char *));
  rewind(fastaFile);
  i=0;
  while(fscanf(fastaFile, "%s", dummy)>0){
    end2clone(myclone, dummy);
    index = hash(myclone , cnt, hashtype); // hope it works
    i++;
    fprintf(stderr, "\r%d\tof\t%d", i,cnt);
    if (names[index] == NULL)
      names[index] = (char *) malloc(sizeof(char)*(strlen(myclone)+1));
    strcpy(names[index], myclone);
  }

  fprintf(stderr, "\n");
  return cnt;
}

int readSingleFasta(FILE *fastaFile, int hashtype){
  int cnt;
  char ch;
  int i;
  int seqcnt=0, seqlen=0;
  int clonecnt;
  int maxnamelen;
  int maxlen;
  char dummy[SEQ_LENGTH];
  char myname[SEQ_LENGTH];
  char myclone[SEQ_LENGTH];
  int index;
  char *thisseq;
  char *thisrevseq;

  cnt = 0; i=0;
  fprintf(stderr, "Counting sequences.\n");
  maxlen=0; maxnamelen=0;
  rewind(fastaFile);
  while (fscanf(fastaFile, "%c", &ch) > 0){
    if (ch == '>'){
      fscanf(fastaFile, "%s", dummy);
      if (strlen(dummy)>maxnamelen)
	maxnamelen=strlen(dummy);
      if (seqlen>maxlen)
	maxlen=seqlen;
      cnt++;
      seqlen=0;
      my_fgets(dummy, SEQ_LENGTH, fastaFile);
    }
    else if (!isspace(ch))
      seqlen++;
  }

  seqcnt = cnt;  
  clonecnt = seqcnt;
    //seqcnt/2+1; // divide fwds & reverses 
  if (seqlen>maxlen)
    maxlen=seqlen;

  //length = (int *) malloc((clonecnt) * sizeof(int));


  //int hash(char *this_name, int nseq, int mode)

  cnt = 0; i=0;

  fprintf(stderr, "Allocating memory for %d sequences with max length %d, maxnamelength %d.\n", clonecnt, maxlen, maxnamelen);
  
  
  forwards = (char **) malloc((clonecnt) * sizeof(char *));
  
  /*
  for (i=0; i<clonecnt; i++)
    forwards[i] = (char *) malloc(maxlen);
  */

  reverses = (char **) malloc((clonecnt) * sizeof(char *));
  
  /*
  for (i=0; i<clonecnt; i++)
    reverses[i] = (char *) malloc(maxlen);
  */

  names = (char **) malloc((clonecnt) * sizeof(char *));

  /*
  for (i=0; i<clonecnt; i++)
    names[i] = (char *) malloc(SEQ_LENGTH);
  */
  

  thisseq = (char *) malloc(maxlen+1);
  thisrevseq = (char *) malloc(maxlen+1);
  
  

  for (i=0; i<clonecnt; i++){
    forwards[i] = NULL;
    reverses[i] = NULL;
    names[i] = NULL;
    /*
    forwards[i][0] = 0;
    reverses[i][0] = 0;
    names[i][0] = 0; */
  }
   
  // FIX THE REST for FWD/REV

  fprintf(stderr, "Reading sequences.\n");
  rewind(fastaFile);
  cnt = -1; 
  while (fscanf(fastaFile, "%c", &ch) > 0){
    if (ch == '>'){
      cnt++;
      my_fgets(myname, SEQ_LENGTH, fastaFile);
      myname[strlen(myname)-1] = 0;

      end2clone(myclone, myname);
      
      index = hash(myclone , clonecnt, hashtype);
      /*
	fprintf(stderr, "\r%d\tof\t%d\tindex\t%d\tname\t%s\tclone\t%s", (cnt+1), clonecnt*2, index, myname, myclone);
      */
      fprintf(stderr, "\r%d\tof\t%d", (cnt+1),seqcnt);
      if (names[index] == NULL)
	names[index] = (char *) malloc(sizeof(char)*(strlen(myclone)+1));
      strcpy(names[index], myclone);

    }
    i = 0;
    if (cnt != 0){
      thisseq[i++] = ch;
      /*
      if (isReverse(myname))
	reverses[index][i++] = ch;
      else
	forwards[index][i++] = ch;
      */
    }
    do{
      if (!(fscanf(fastaFile, "%c", &ch) > 0))
	break;
      if (ch!='>' && ch!='\r' && ch!='\n'){
	thisseq[i++] = ch;
	/*
	if (isReverse(myname))
	  reverses[index][i++] = ch;
	else
	  forwards[index][i++] = ch;
	*/
      }
    } while (ch != '>');
    
    thisseq[i++] = 0;

    if (isReverse(myname)){
      rcomp(thisseq, thisrevseq);
      reverses[index] = (char *)malloc(sizeof(char)*(strlen(thisseq)+1));
      strcpy(reverses[index], thisrevseq);
    }
    else{
      forwards[index] = (char *)malloc(sizeof(char)*(strlen(thisseq)+1));
      strcpy(forwards[index], thisseq);
    }
    
    if (ch == '>'){
      cnt++;
      fprintf(stderr, "\r%d\tof\t%d", (cnt+1), seqcnt );
      if (cnt != seqcnt){
	  my_fgets(myname, SEQ_LENGTH, fastaFile);
	  myname[strlen(myname)-1] = 0;
	  
	  end2clone(myclone, myname);
	  
	  index = hash(myclone , clonecnt, hashtype); 
	  if (names[index] == NULL)
	    names[index] = (char *) malloc(sizeof(char)*(strlen(myclone)+1));
	  strcpy(names[index], myclone);
	  
      }
    } // if
  } // while
  
  
  fprintf(stderr, "\n[OK] %d sequences read from fasta file.\n",clonecnt);

  free(thisseq); free(thisrevseq);

  return clonecnt;

}

/* Read the next FASTA record from f.
   On success returns 1 and stores freshly malloc'd strings in *name_out
   (header line, without the leading '>') and *seq_out (sequence with all
   whitespace stripped). Returns 0 at end of file. The caller frees both. */
int readFastaRecord(FILE *f, char **name_out, char **seq_out){
  int c;
  int ncap = 0, nlen = 0;
  int scap = 0, slen = 0;
  char *name = NULL;
  char *seq = NULL;

  /* find the start of the next record */
  while ((c = fgetc(f)) != EOF && c != '>')
    ;
  if (c == EOF)
    return 0;

  /* header line (up to end of line) */
  while ((c = fgetc(f)) != EOF && c != '\n' && c != '\r'){
    if (nlen + 1 >= ncap){
      ncap = ncap ? ncap * 2 : 128;
      name = (char *) realloc(name, ncap);
    }
    name[nlen++] = c;
  }
  if (name == NULL)
    name = (char *) malloc(1);
  name[nlen] = 0;

  /* sequence: everything up to the next '>' or EOF, whitespace removed */
  while ((c = fgetc(f)) != EOF){
    if (c == '>'){
      ungetc(c, f);
      break;
    }
    if (isspace(c))
      continue;
    if (slen + 1 >= scap){
      scap = scap ? scap * 2 : 256;
      seq = (char *) realloc(seq, scap);
    }
    seq[slen++] = c;
  }
  if (seq == NULL)
    seq = (char *) malloc(1);
  seq[slen] = 0;

  *name_out = name;
  *seq_out = seq;
  return 1;
}

/* Two-file input: forward reads in one file, reverse reads in another, paired
   by position (the i-th record of each file form a forward/reverse pair, so no
   name matching is needed). Fills the global forwards/reverses/names arrays and
   returns the number of pairs. Reverse reads are reverse-complemented, matching
   how readSingleFasta stores them. */
int readPairedFasta(FILE *fwd, FILE *rev){
  int cap = 0, n = 0;
  char *fwdname, *fwdseq;
  char *revname, *revseq;

  forwards = NULL;
  reverses = NULL;
  names = NULL;

  fprintf(stderr, "Reading paired forward/reverse sequences.\n");

  while (readFastaRecord(fwd, &fwdname, &fwdseq)){
    char *clone;
    char *revcomp;

    if (!readFastaRecord(rev, &revname, &revseq)){
      fprintf(stderr,
	      "Warning: forward file has more reads than reverse; extra forward reads ignored.\n");
      free(fwdname); free(fwdseq);
      break;
    }

    if (n >= cap){
      cap = cap ? cap * 2 : 1024;
      forwards = (char **) realloc(forwards, cap * sizeof(char *));
      reverses = (char **) realloc(reverses, cap * sizeof(char *));
      names    = (char **) realloc(names,    cap * sizeof(char *));
    }

    /* forward read stored as-is */
    forwards[n] = fwdseq;

    /* reverse read stored reverse-complemented */
    revcomp = (char *) malloc(strlen(revseq) + 1);
    rcomp(revseq, revcomp);
    reverses[n] = revcomp;
    free(revseq);

    /* clone identifier for the .names file (strip any FORWARD/REVERSE tag) */
    clone = (char *) malloc(strlen(fwdname) + 1);
    end2clone(clone, fwdname);
    names[n] = clone;

    free(fwdname);
    free(revname);
    n++;

    fprintf(stderr, "\r%d pairs", n);
  }

  /* drain any leftover reverse reads just to report the mismatch */
  if (readFastaRecord(rev, &revname, &revseq)){
    fprintf(stderr,
	    "\nWarning: reverse file has more reads than forward; extra reverse reads ignored.\n");
    free(revname); free(revseq);
  }

  fprintf(stderr, "\n[OK] %d forward/reverse pairs read.\n", n);
  return n;
}

void usage(char *prog){
  fprintf(stderr,
"RepeatNet: an ab initio centromeric sequence detection algorithm\n"
"\n"
"Usage: %s [options]\n"
"\n"
"Input (choose one):\n"
"  -i, --input FILE         interleaved FASTA (forward/reverse reads in one file)\n"
"  -f, --forward FILE       forward reads FASTA (use together with --reverse)\n"
"  -r, --reverse FILE       reverse reads FASTA (paired by position with --forward)\n"
"\n"
"Parameters:\n"
"  -H, --hash N             hash function type (default 11)\n"
"  -e, --ecut F             edge-weight cutoff (default 0.0)\n"
"  -a, --acut N             absolute co-occurrence cutoff (0 = use --ecut)\n"
"  -c, --ccut N             vertex count cutoff (default 1)\n"
"  -m, --merge              merge reverse-complement k-mers before comparing\n"
"  -t, --hashtest N         benchmark hash functions 3..N-1 and exit\n"
"\n"
"Pipeline stages:\n"
"      --loadwin FILE       load a .winlog produced by an earlier run\n"
"      --loadnames FILE     load a .names file\n"
"      --loadmatrix FILE    load a binary .matrix file\n"
"      --getn FILE          list of vertices whose neighbors to report (with --loadmatrix)\n"
"      --compare            build the co-occurrence graph/matrix\n"
"      --components         split a loaded matrix into connected components\n"
"      --clones FILE        dump clone names for the vertices listed in FILE\n"
"      --trim FILE          drop the windows listed in FILE from a loaded .winlog\n"
"      --edgehist           also write an edge-weight histogram\n"
"      --decode FILE        decode integer window codes back to k-mer strings\n"
"      --encode KMER        print the integer code for a k-mer and exit\n"
"      --help               show this message\n",
	  prog);
}

void my_fgets(char *str, int length, FILE *in){
  char ch;
  int i=0;
  while (i<length && fscanf(in, "%c", &ch)){
    if ((ch==' '  || ch=='\t')  && i!=0){
      if (str[i-1] == ' ' || str[i] == '\t')
	;
      else
	str[i++] = ch;
    }
    else
      str[i++] = ch;
    if (ch == '\n' || ch=='\r')
      break;
  }
  str[i] = 0;
}


int hash(char *this_name, int nseq, int type){
  int i;
  int len = strlen(this_name);
  unsigned long sum;
  int index;
  unsigned int a,b;

  sum = 0;
  a=378551;
  b=63689;

  if (type == 0){
    for (i=0;i<len;i++){
      if (this_name[i]!=' ')
	sum += i*this_name[i]+i;
    }
  }
  else if (type == 1){
    for (i=0;i<len;i++){
      if (this_name[i]!=' ')
	sum += this_name[i];
    }  
  }
  else if (type == 2){
    for (i=0;i<len;i++){
      if (this_name[i]!=' ')
	sum += i*this_name[i];
    }
  }
  else if (type == 3){
    for (i=0;i<len;i++){
      if (this_name[i]!=' ')
	sum += this_name[i] * pow(7,i); // rabin-karp style
    }
  }
  else if (type == 4){
    for (i=0;i<len;i++){
      if (this_name[i]!=' ')
	sum += this_name[i] * pow(11,i); // rabin-karp style, part 2
    }
  }
  else if (type == 5){
   for (i=0;i<len;i++){
     sum += sum*a + this_name[i];
     a = a*b; //RShash
   }
  }
  else if (type == 6){
    sum=1315423911;
   for (i=0;i<len;i++){
     sum ^= ((sum<<5) + this_name[i] + (sum>>2)); // JShash
   }
  }
  else if (type == 7){
   for (i=0;i<len;i++){
     sum = (sum << 4) + this_name[i];
     if ((a = sum & 0XFF0000000L) != 0)
       sum ^= a>>24;
     sum &= ~a; // ELFhash
   }  
  }
  else if (type == 8){
    a=131;
    for (i=0;i<len;i++){
      sum=(sum*a)+this_name[i]; //bkdrhash
    }
  }
  else if (type == 9){
    for (i=0;i<len;i++){
      sum = this_name[i] + (sum<<6) + (sum<<16) - sum; //sdbmhash
    }
  }
  else if (type == 10){
    for (i=0;i<len;i++){
      sum = ((sum<<5) + sum) + this_name[i]; // djb
    }
    
  }
  else if (type == 11){
    for (i=0;i<len;i++){
      sum=((sum<<5)^(sum>>27))^this_name[i]; //DEK
    }
  }
  else if (type == 12){
    for (i=0;i<len;i++){
      sum = sum<<7 ^ this_name[i];
    }
  }
  else if (type == 13){
    for (i=0;i<len;i++){
      a=0x811C9DC5;
      sum*=a;
      sum^=this_name[i]; // fnv
    }
  }
  else if (type == 14){
    for (i=0;i<len;i++){
      sum ^= ((i&1) == 0)?((sum<<7)^this_name[i]^(sum>>3)):
	(~((sum<<11)^this_name[i]^(sum>>5))); //ap
    }
  }    

  
  index = sum % nseq;

  
  while (names[index] != NULL && strcmp(names[index], this_name)){
    collisions++;
    index++;
    if (index == nseq)
      index = 0;
  }

  return index;

}

void end2clone(char *clone, char *end){
  int i, len, j;
  int flag = 0;
  char *rloc, *floc;
  char *pass;
  
  len = strlen(end);

  rloc = strstr(end, "REVERSE");
  floc = strstr(end, "FORWARD");
  
  if (floc!=NULL || rloc!=NULL){
    
    if (rloc != NULL)
      pass = rloc;
    else
      pass = floc;
    
    i = 0; j=0;
    while (end+i+1 != pass){
      clone[j] = end[i];
      i++; j++;
    }
    // we're at the '.' before REVERSE/FORWARD
    // jump over that
    i++;
    // find the '.' after REVERSE/FORWARD
    while (end[i]!='.') i++;
    // don't skip that
    // copy the rest
    while(i<len)
      clone[j++] = end[i++];
    // terminate
    clone[j] = 0;
    return;
  }


  j=0;
  for (i=0; i<len; i++){
    if ((end[i] != 'F' && end[i] != 'R') || flag){
      clone[j++] = end[i];  
    }
    else if (end[i] != 'F' || end[i] != 'R'){
      if (flag)
	clone[j++] = end[i];  
      else{
	flag = 1;
      }
    } 
    
  }
  clone[j] = 0;

  //fprintf(stderr, "name:\t%s\tclone:\t%s\n", end, clone);
  
}


void rcomp(char *window, char *rcomp){
  /* reverse complement */
  int i;
  int len = strlen(window);
  for (i=0;i<len;i++)
    switch (toupper(window[i])){
    case 'A':
      rcomp[len-i-1] = 'T';
      break;
    case 'C':
      rcomp[len-i-1] = 'G';
      break;
    case 'G':
      rcomp[len-i-1] = 'C';
      break;
    case 'T':
      rcomp[len-i-1] = 'A';
      break;
    default: // may be N
      rcomp[len-i-1] = window[i];
      break;
    }    
  rcomp[len] = 0;
}


int isReverse(char *name){
  int i, len;
  int isNum = 0;

  if(strstr(name, "REVERSE"))
    return 1;
  if(strstr(name, "FORWARD"))
    return 0;


  len = strlen(name);
  

  for (i=0;i<len;i++){
    if (isdigit(name[i]))
      isNum = 1;
    if (isNum && name[i] == 'R')
      return 1;
    if (isNum && name[i] == 'F')
      return 0;
  }
  
  return 0;
}

int encode_window(char *window){
  /* 2-bit encoding: A:0, C:1, G:2, T:3, packed MSB-first. */
  return encode_window2(window);
}

int encode_window2(char *window){
  int len = strlen(window);
  unsigned int code;
  //int totalbits = sizeof(int) * 8;
  unsigned int this;
  int shiftleft;
  //int shiftright;
  int i;
  unsigned int left; //, right;

  code = 0;
  
  for (i=0;i<len;i++){
    this = cindex(window[i]);
    shiftleft = 2*(len-i-1);


    left = this << shiftleft;
    code = code | left;
  }
  return code;
}

int revcomp_encoded(unsigned int code){
  /* Decode to text, reverse-complement, then re-encode. */
  char thiswin[20];
  char reverse[20];

  decode_window(code, thiswin);
  rcomp(thiswin, reverse);
  return encode_window2(reverse);
}


unsigned int cindex(char ch){
  unsigned int this;
  switch(toupper(ch)){
  case 'A':
    this = 0;
    break;
  case 'C':
    this = 1;
    break;
  case 'G':
    this = 2;
    break;
  case 'T':
    this = 3;
    break;
  default:
    assert(0);
    break;
  }
  return this;
}

char reverseindex(int code){
  char ch;
  switch(code){
  case 0:
    ch = 'A';
    break;
  case 1:
    ch = 'C';
    break;
  case 2:
    ch = 'G';
    break;
  case 3:
    ch = 'T';
    break;
  default:
    assert(0);
    break;
  }
  return ch;
}

int getWindow(char *genome, int glen, char *window, int win_start){
  int i;
  int win_end;
  
  win_end = win_start+WS;
  if (win_start >= glen)
    return 0;
  if (win_end>=glen)
    return 0;
  //    win_end = glen;
  for (i=win_start;i<win_end;i++){
    if (genome[i] == 'N' || genome[i] == 'O'){
      return 0;
    }
    window[i-win_start] = genome[i];
    
  }

  window[i-win_start] = 0;


  return 1;
}

void do_end(int index, int orient){
  char *window;
  int i;
  int willdo;
  int flen;
  int code;
  if (orient == FORWARD)
    flen = strlen(forwards[index]);
  else
    flen = strlen(reverses[index]);

  /* getWindow writes WS characters plus a terminating '\0' at window[WS],
     so the buffer needs WS+1 bytes. */
  window = (char *) malloc(sizeof(char)*(WS+1));

  i=0;
  if (orient == FORWARD){
    while (i<flen-WS+1){
      willdo = getWindow(forwards[index], flen, window, i);  
      
      if (willdo){
	code = encode_window(window);
	windowflag[code] = index;
      }
      i++;
    }
  }

  else{
    while (i<flen-WS+1){
      willdo = getWindow(reverses[index], flen, window, i);  
      if (willdo){
	code = encode_window(window);
	if (windowflag[code] == index){
	  // index is the sequence number 0...nseqs-1
	  // code is the window code AAA..AA=0, etc.
	  //printf("DOUBLE SIDED %s\n", window);
	  insert_wloc(&(winlocs[code]), index, code, CREATE);
	}
      }
      i++;
    }
  }

  free(window);
}

void insert_wloc(struct wpair **winloc, int index, int code, int mode){
  //printf("Code: %d\t index: %d\n", code, index);
  struct wpair *new;
  
  if ((*winloc) == NULL){
    (*winloc) = (struct wpair *) malloc(sizeof(struct wpair));
    (*winloc)->id = index;
    (*winloc)->next = NULL;
    if (mode == CREATE)
      windowcnt[code]++;
    //printf("insert %d\tcnt %d\t code %d\n", index, windowcnt[code], code);
    return;
  }
  
  if ((*winloc)->id == index) // no double insertions
    return;
  
  new = (struct wpair *) malloc(sizeof(struct wpair));
  new->id = index;
  new->next = (*winloc);
  if (new->id > (*winloc)->id)
    WINORDER = DESCENDING;
  else
    WINORDER = ASCENDING;


  (*winloc) = new;
  if (mode == CREATE)
    windowcnt[code]++;
  //printf("insert %d\tcnt %d\t code %d\n", index, windowcnt[code], code);

}


void compare_wloc(struct wpair *wloc1, struct wpair *wloc2, int *eboth){
  struct wpair *w1, *w2;
  int eb;

  eb = 0;

  w1 = wloc1;
  w2 = wloc2;

  // they are both sorted in DESCENDING ORDER if mode == 0
  
  //printf("Mode : %s\n", WINORDER==ASCENDING ? "ascending" : "descending");

  if (WINORDER == ASCENDING){

    while (w1!=NULL && w2!=NULL){
      //printf("Check %d with %d\n", w1->id, w2->id);
      //getchar();
      if (w1->id == w2->id){
	eb++;
	w1=w1->next;
	w2=w2->next;
      }
      
      else if (w1->id < w2->id){
	// we wont find w1.id in the rest of w2 anyway
	w1=w1->next;
      }
      else{
	// we wont find w2.id in the rest of w1 anyway
	w2=w2->next;
      }
      
    }
  }

  else{ // if descending
    while (w1!=NULL && w2!=NULL){
      //printf("Check %d with %d\n", w1->id, w2->id);
      //getchar();
      if (w1->id == w2->id){
	eb++;
	w1=w1->next;
	w2=w2->next;
      }
      
      else if (w1->id < w2->id){
	// we wont find w2.id in the rest of w1 anyway
	w2=w2->next;
      }
      else{
	// we wont find w1.id in the rest of w2 anyway
	w1=w1->next;
      }
      
    }
  }


  (*eboth) = eb;
}

void dumpwins(FILE *winlog, int index){
  struct wpair *w1;
  w1 = winlocs[index];

  fwrite(&index, sizeof(int), 1, winlog);
  fwrite(&(windowcnt[index]), sizeof(int), 1, winlog);
  while (w1!=NULL){
    fwrite(&(w1->id), sizeof(int), 1, winlog);
    w1=w1->next;
  }
  //fwrite(&separator, sizeof(int), 1, winlog);
}

void loadwins(FILE *winlog){
  int cnt;
  int index;
  int myid;
  char ch;
  int i;
  int status;

  fscanf(winlog, "%d", &WS);
  // pass newline
  fscanf(winlog, "%c", &ch);

  wincnt = pow(4, WS);
  printf("read WS %d\twincnt %d\n", WS, wincnt);
  winlocs = (struct wpair **)malloc(sizeof(struct wpair *)*wincnt);
  windowflag = (int *) malloc(sizeof(int)*wincnt);
  windowcnt = (int *) malloc(sizeof(int)*wincnt);


  for (i=0;i<wincnt;i++){
    winlocs[i]=NULL;
    windowflag[i] = -1;
    windowcnt[i] = 0;
  }
  
  

  while (!feof(winlog)){
    status = fread(&index, sizeof(int), 1, winlog);
    if (status == 0)
      break;
    fread(&cnt, sizeof(int), 1, winlog);
    //printf("read %d\t%d\n", index,cnt);
    windowcnt[index] = cnt;
    //printf("read %d\t%d\n", index, windowcnt[index]);
    //getchar();
    for (i=0;i<cnt;i++){
      fread(&myid, sizeof(int), 1, winlog);
      insert_wloc(&(winlocs[index]), myid, index, LOAD);      
    }   
  }

}


void do_loaded(FILE *in, char *fname){
  FILE *histogram;
  int i;
  char fname2[SEQ_LENGTH];

  // CALKAN: Jun 13, 12:37pm.
  // something's wrong with counter
  
  fprintf(stdout, "Loading windows...");
  fflush(stdout);

  loadwins(in);
  
  sprintf(fname2, "%s.histogram", fname);
  histogram = fopen(fname2, "w");
  
  for (i=0;i<wincnt;i++){
    fprintf(histogram, "%d\n", windowcnt[i]);
  }
  fclose(histogram);
  fprintf(stdout, "\nWindows loaded.\n");
  
}

void decode_window(int code, char *window){
  int i;
  int shiftleft, shiftright;
  unsigned int left, right;
  int totalbits = sizeof(int)*8;
  for (i=0;i<WS;i++){
    shiftleft = totalbits - 2*WS + 2*i;
    shiftright = totalbits - 2;
    /* shift on an unsigned value: shifting into the sign bit of a signed
       int is undefined behavior. */
    left = (unsigned int)code << shiftleft;
    right = left >> shiftright;
    window[i] = reverseindex(right);
  }
  window[WS]=0;
}



void compare(char *fname, int hashtype, float ecut, int ccut, int acut, int merge){
  int i,j ;
  char fname2[SEQ_LENGTH];
  FILE *graph, *matlog, *edgelog;
  int e1, e2, eboth;
  float eweight;
  float thiscut=ecut;
  int dumped=-12345;

  if (merge)
    compactwins();
  
  if (acut!=0) thiscut=acut;
  
  fprintf(stdout, "acut: %d\nthiscut: %d\n", acut, (int)thiscut);

  if (!merge)
    sprintf(fname2, "%s.h%d.cut%d.e%d.viz", fname, hashtype, ccut, (int)(thiscut*100));
  else
    sprintf(fname2, "%s.h%d.cut%d.e%d.merged.viz", fname, hashtype, ccut, (int)(thiscut*100));
  graph = fopen(fname2, "w");
  if (!merge)
    sprintf(fname2, "%s.h%d.cut%d.e%d.matrix", fname, hashtype, ccut, (int)(thiscut*100));
  else
    sprintf(fname2, "%s.h%d.cut%d.e%d.merged.matrix", fname, hashtype, ccut, (int)(thiscut*100));
  matlog = fopen(fname2, "w");


  if (EDGEHIST){
    if (!merge)
      sprintf(fname2, "%s.h%d.cut%d.e%d.edges", fname, hashtype, ccut, (int)(thiscut*100));
    else
      sprintf(fname2, "%s.h%d.cut%d.e%d.merged.edges", fname, hashtype, ccut, (int)(thiscut*100));
    edgelog = fopen(fname2, "w");
  }


  fprintf(graph, "graph G{\n");
  //fprintf(graph, "\tnode[style=filled,label=\"\"];\n");

  // NOTE: THE FILE HEADER CHANGED BY 4 BYTES 
  // ALSO WS IS NOW WRITTEN BINARY

  fwrite(&dumped, sizeof(int), 1 , matlog);
  dumped = 0;

  fwrite(&WS, sizeof(int), 1 , matlog);
  //fprintf(matlog,"%d\n",WS);

  printf("Mode : %s\n", WINORDER==ASCENDING ? "ascending" : "descending");


  fprintf(stdout, "Comparing...");
  fflush(stdout);

  /* Only windows with an occupied location list and count >= ccut can
     produce an edge, and there are far fewer of those than the 4^WS slots.
     Collect them once (in ascending index order, so i<j is preserved) and
     compare only those pairs instead of scanning the whole space. */
  int *active = (int *) malloc(sizeof(int) * wincnt);
  int nactive = 0;
  for (i=0;i<wincnt;i++)
    if (windowcnt[i] >= ccut && winlocs[i] != NULL)
      active[nactive++] = i;

  int ai, aj;
  for (ai=0; ai<nactive; ai++){
    i = active[ai];
    for (aj=ai+1; aj<nactive; aj++){
      j = active[aj];

      compare_wloc(winlocs[i], winlocs[j], &eboth);

      e1 = windowcnt[i];
      e2 = windowcnt[j];

      eweight = 2*(float)eboth/(float)(e1+e2);

      if (acut!=0) eweight=eboth;

      if (eweight>thiscut){
	// if ecut=0.0, then return any overlap (1)
	fprintf(graph, "\t%d -- %d;\n", i, j);
	fwrite(&i,sizeof(int),1,matlog);
	fwrite(&j,sizeof(int),1,matlog);
	fwrite(&eboth,sizeof(int),1,matlog);
	dumped++;
	if(EDGEHIST)
	  fprintf(edgelog, "%d\n", eboth);
      }
    }
  }

  free(active);

  fprintf(graph, "}\n");
  fclose(graph);

  // replace dump count
  rewind(matlog);
  fwrite(&dumped, sizeof(int), 1 , matlog);

  fclose(matlog);
  if (EDGEHIST)
    fclose(edgelog);
  printf("matsize: %d\nWS: %d\n", dumped, WS);

  fprintf(stdout, "\nCompared...\n");

}

void compactwins(void){
  // concatenate linked lists of kmers that are revcomp'ed
  // use concatlocations from pfast
  // do this for loaded wins for now.
  // the numbers are in ASCENDING ORDER

  int i;
  int rci;

  fprintf(stdout, "Merging...");

  for (i=0;i<wincnt;i++){
    if (windowcnt[i]==0)
      continue;
    rci = revcomp_encoded(i);
    if (windowcnt[rci]==0)
      continue;
    if (rci<=i)
      continue;
    merge(i, rci);
  }
  
  fprintf(stdout, "\nMerged.\n");
}

void merge(int i, int rci){
  struct wpair *w1, *w2, *new, *w1p;
  // winloc_f = winloc_f+winloc_r
  w1 = winlocs[i];
  w1p = w1;
  w2 = winlocs[rci];
  

  if (WINORDER == ASCENDING){
    while (w2 != NULL){
      if (w2->id < w1->id){
	new = (struct wpair *) malloc(sizeof(struct wpair));
	new->id = w2->id;
	new->next = w1;
	if (w1p == winlocs[i] && w1p->id > new->id){
	  winlocs[i] = new;
	  w1p = winlocs[i];
	  w1 = winlocs[i];
	}
	else{
	  w1p->next = new;
	  w1 = new;
	}
	w2 = w2->next;
      } 
      else if (w2->id > w1->id){
	// no copy here, w2 stays the same, w1 moves
	w1p = w1;
	w1 = w1->next;
      }
      else{
	// they are equal; move along
	w1p = w1;
	w1 = w1->next;
	w2 = w2->next;
      }
      if (w1p->next==NULL && w2!=NULL){
	w1p->next = copyWins(w2);
	break;
      }
    }
  }
  
  else{
    while (w2 != NULL){
      if (w2->id > w1->id){
	new = (struct wpair *) malloc(sizeof(struct wpair));
	new->id = w2->id;
	new->next = w1;
	if (w1p == winlocs[i]  && w1p->id < new->id){
	  winlocs[i] = new;
	  w1p = winlocs[i];
	  w1 = winlocs[i];
	}
	else{
	  w1p->next = new;
	  w1 = new;
	}
	w2 = w2->next;
      } 
      else if (w2->id < w1->id){
	// no copy here, w2 stays the same, w1 moves
	w1p = w1;
	w1 = w1->next;
      }
      else{
	// they are equal; move along
	w1p = w1;
	w1 = w1->next;
	w2 = w2->next;
      }
      if (w1p->next==NULL && w2!=NULL){
	w1p->next = copyWins(w2);
	break;
      }
    }
  }
  
  freewins(winlocs[rci]);
  windowcnt[i] += windowcnt[rci];
  windowcnt[rci] = 0;

}

struct wpair *copyWins(struct wpair *w){
  struct wpair *new;
  struct wpair *ret;
  struct wpair *this;
  struct wpair *current;

  if (w==NULL)
    return NULL;
 
  new = (struct wpair *) malloc(sizeof(struct wpair));
  new->id = w->id;

  new->next = NULL;
  ret = new;
  this = ret;
  current = w->next;

  while(current != NULL){
    new = (struct wpair *) malloc(sizeof(struct wpair));
    new->id = current->id;
    new->next = NULL;
    this->next = new;
    this = this->next;
    current = current->next;
  }

  return ret;
}

void freewins(struct wpair *w){
  struct wpair *tmp;
  while (w!=NULL){
    tmp = w->next;
    free(w);
    w = tmp;
  }
}

int loadMatrix(char *fname, char *fname2){
  FILE *matlog;
  FILE *neighbor;
  FILE *nfile;
  char name[200];
  int i, j, eboth;
  int nsearch[5000]; // for now
  int cnt=0;
  int k;
  FILE *degrees;
  int *matdegrees;
  int status;
  int dumped;
  int loadcnt;

  fprintf(stdout, "Loading matrix...");
  fflush(stdout);

  matlog = myfopen(fname, "r");

  if (fname2[0]!=0){
    neighbor = myfopen(fname2, "r");
    while (fscanf(neighbor, "%d", &(nsearch[cnt])) > 0)
      printf("Look for: %d\n", nsearch[cnt++]);
    fclose(neighbor);
    sprintf(name, "%s.neighbors", fname2);
    nfile = myfopen(name, "w");
  }


  fread(&dumped, sizeof(int), 1, matlog);
  fread(&WS, sizeof(int), 1, matlog);

  printf("matsize: %d\nWS: %d\n", dumped, WS);

  interactions = (struct imatrix *) malloc(sizeof(struct imatrix) * dumped);
  loadcnt = 0;

  //fscanf(matlog, "%d", &WS);
  //fscanf(matlog, "%c", &ch); // skip newline

  wincnt = pow(4, WS);
  
  matdegrees = (int *)malloc(sizeof(int)*wincnt);


  for (i=0;i<wincnt;i++) {
    matdegrees[i]=0;
  }

  while (!feof(matlog)){
    status = fread(&i, sizeof(int), 1, matlog);
    if (status == 0) break;
    fread(&j, sizeof(int), 1, matlog);
    matdegrees[i]++; matdegrees[j]++;
    fread(&eboth, sizeof(int), 1, matlog);    
    


    interactions[loadcnt].i = i;
    interactions[loadcnt].j = j;
    interactions[loadcnt++].eweight = eboth;

    

    // NOTE: i < j ALWAYS
    for (k=0;k<cnt;k++){ 
      if (i==nsearch[k]){
	fprintf(nfile, "%d\n", j);
	break;
      }
      else if (j==nsearch[k]){
	fprintf(nfile, "%d\n", i);
	break;
      }
    }
  }
  fclose(matlog);
  if (fname2[0]!=0) fclose(nfile);
  fprintf(stdout, "\nMatrix loaded.\n");

  sprintf(name, "%s.degrees", fname);
  degrees = fopen(name, "w");
  for (i=0;i<wincnt;i++) 
    if (matdegrees[i]!=0)
      fprintf(degrees, "%d\t%d\n", i, matdegrees[i]);

  free(matdegrees);

  return dumped;

}

FILE *myfopen(char *fname, char *mode){
  FILE *ret;
  ret = fopen(fname, mode);
  if (ret == 0){
    fprintf(stderr, "File %s not found.\n", fname);
    exit(0);
  }
  return ret;
}

void decode_wins(char *fname){
  int i;
  char decoded[20];
  FILE *in;
  FILE *out;
  char name[200];
  sprintf(name, "%s.decoded", fname);
  in = myfopen(fname, "r");
  out = myfopen(name, "w");

  while (fscanf(in, "%d", &i)>0){
    decode_window(i, decoded);
    fprintf(out, ">%d\n%s\n", i, decoded);
  }
  fclose(in); fclose(out);
}

void dumpclones(char *fname){
  FILE *in;
  FILE *out;
  int i;
  char name[200];
  struct wpair *w;

  // for now, i'm assuming names[] are loaded already

  sprintf(name, "%s.clones", fname);
  
  in = myfopen(fname, "r");
  out = myfopen(name, "w");
  
  while(fscanf(in, "%d", &i) > 0) {
    printf("Clones for: %d\n", i);
    if (windowcnt[i]!=0 && winlocs[i]!=NULL){ // over cautious
      w = winlocs[i];
      while (w!=NULL){
	fprintf(out, "%s\n", names[w->id]);
	w = w->next;
      }
    }
  }
  fclose(in); fclose(out);
}


void loadNames(char *fname){
  FILE *in; 
  char thisname[SEQ_LENGTH];
  int i;
  in = fopen(fname, "r");

  fscanf(in, "%d", &nseq);
  names = (char **) malloc((nseq) * sizeof(char *));
  fprintf(stdout, "Loading %d names...", nseq);
  fflush(stdout);

  while (fscanf(in, "%d %s", &i, thisname) > 0){
    names[i] = (char*) malloc((strlen(thisname)+1) * sizeof(char));
    strcpy(names[i], thisname);
  }
  
  fclose(in);

  fprintf(stdout, "\nNames loaded.\n");
}


//gozde modifications

#define UNCLUSTERED 0 //sequence is not in a connected component
#define SEQ_EXISTS 1    //sequence exists in the graph
#define NULL_SEQ_ID -1

//sequence Node
typedef struct listNode{
  int seqId;
  struct listNode* next;
}_listNode;

//interaction Node
typedef struct interNode{
  int v1;
  int v2;
  int w;
  struct interNode* next;
}_interNode;

void connComp(int matrixcnt, int wincnt){

  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);  

  printf("Computing connected components...\n");

  //first find the number of nodes in the graph
  int* seqs = (int*) calloc(wincnt, sizeof(int));
  int numNodes = 0;
  int v1, v2, w;
  int tmp;

  for (tmp = 0; tmp < matrixcnt; tmp++) {
    v1 = interactions[tmp].i;
    v2 = interactions[tmp].j;
    w = interactions[tmp].eweight;
    seqs[v1] = SEQ_EXISTS;
    seqs[v2] = SEQ_EXISTS;
  } 

  for (tmp = 0; tmp < wincnt; tmp++) {
    if (seqs[tmp] == SEQ_EXISTS)
      numNodes++;
  }

  free(seqs);

  printf("total number of nodes in the graph:%d\ntotal number of interactions:%d\n", numNodes, matrixcnt);

  //keep an array to keep the component of each node
  int* seqComp = (int*) calloc(wincnt, sizeof(int));
  
  //keep the components in a two dimensional listNode array
  struct listNode** components = (struct listNode**) malloc((numNodes+1) * sizeof(struct listNode*));

  //keep the last elements of components in a one dimensional listNode array
  struct listNode** lastNodes = (struct listNode**) malloc((numNodes+1) * sizeof(struct listNode*));

  //keep the interaction data in a two dimensional interNode array
  struct interNode** interNodes = (struct interNode**) malloc((numNodes+1) * sizeof(struct interNode*));
  //keep the last elements of interactions in a one dimensional interNode array
  struct interNode** lastInterNodes = (struct interNode**) malloc((numNodes+1) * sizeof(struct interNode*));
  
  int numComponents = 0; //current number of components

  struct listNode *compItr, *newNode, *newNode2;
  int seq1Comp, seq2Comp;
  struct interNode *interNodesItr, *newInterNode;

  int k;
  for (k = 0; k < matrixcnt; k++) {
    
    v1 = interactions[k].i;
    v2 = interactions[k].j;
    w  = interactions[k].eweight;
    
    //create the new interaction node
    newInterNode = (struct interNode*) malloc(sizeof(struct interNode));
    newInterNode -> v1 = v1;
    newInterNode -> v2 = v2;
    newInterNode -> w = w;
    newInterNode -> next = NULL;

    //if both sequences are not clustered, put them in a new cluster
    if (seqComp[v1] == UNCLUSTERED && seqComp[v2] == UNCLUSTERED) {
      numComponents++;
      newNode = (struct listNode*) malloc(sizeof(struct listNode));
      newNode -> seqId = v1;
      components[numComponents] = newNode;
      newNode2 = (struct listNode*) malloc(sizeof(struct listNode));
      newNode2 -> seqId = v2;
      newNode2 -> next = NULL;
      components[numComponents] -> next = newNode2;
      lastNodes[numComponents] = newNode2;
      seqComp[v1] = numComponents;
      seqComp[v2] = numComponents;  
      //arrange interaction nodes
      interNodes[numComponents] = newInterNode;
      lastInterNodes[numComponents] = interNodes[numComponents];
    }

    //if one of the sequences is clustered, put the other one in that cluster
    //go to the end of the linked list where v1 resides and put v2
    else if (seqComp[v1] != UNCLUSTERED && seqComp[v2] == UNCLUSTERED) {
      seq1Comp = seqComp[v1];
      newNode = (struct listNode*) malloc(sizeof(struct listNode));
      newNode -> seqId = v2;
      newNode -> next = NULL;
      lastNodes[seq1Comp] -> next = newNode;
      lastNodes[seq1Comp] = newNode;
      seqComp[v2] = seq1Comp;
      //arrange interaction nodes
      lastInterNodes[seq1Comp] -> next = newInterNode;
      lastInterNodes[seq1Comp] = newInterNode;      
    }

    //go to the end of the linked list where v2 resides and put v1
    else if (seqComp[v1] == UNCLUSTERED && seqComp[v2] != UNCLUSTERED) {
      seq2Comp = seqComp[v2];
      newNode = (struct listNode*) malloc(sizeof(struct listNode));
      newNode -> seqId = v1;
      newNode -> next = NULL;
      lastNodes[seq2Comp] -> next = newNode;
      lastNodes[seq2Comp] = newNode;
      seqComp[v1] = seq2Comp;
      //arrange interaction nodes
      lastInterNodes[seq2Comp] -> next = newInterNode;
      lastInterNodes[seq2Comp] = newInterNode;      
    }

    //if both sequences are already in connected components,
    //combine their components
    //add the component of v2 to the end of the component of v1
    //if only they aren't in the same component already
    else if ((seqComp[v1] != UNCLUSTERED && seqComp[v2] != UNCLUSTERED) && 
        (seqComp[v1] != seqComp[v2])) {
      seq1Comp = seqComp[v1];
      seq2Comp = seqComp[v2];
      //update the components of each seq in the about-to-be-deleted component
      compItr = components[seq2Comp];
      while (compItr != NULL) {
        seqComp[compItr->seqId] = seq1Comp;
        compItr = compItr -> next;
      }
      //combine the components and update the lastNodes array
      lastNodes[seq1Comp] -> next = components[seq2Comp];
      lastNodes[seq1Comp] = lastNodes[seq2Comp];      
      components[seq2Comp] = NULL;
      lastNodes[seq2Comp] = NULL;
      //arrange interaction nodes
      lastInterNodes[seq1Comp] -> next = interNodes[seq2Comp];
      lastInterNodes[seq1Comp] = lastInterNodes[seq2Comp];
      interNodes[seq2Comp] = NULL;
      lastInterNodes[seq2Comp] = NULL;
    }

    //if both sequences are already in the same connected component,
    //record their interaction info only
    else if ((seqComp[v1] != UNCLUSTERED && seqComp[v2] != UNCLUSTERED) &&
      (seqComp[v1] == seqComp[v2])) {
      seq1Comp = seqComp[v1];
      lastInterNodes[seq1Comp] -> next = newInterNode;
      lastInterNodes[seq1Comp] = newInterNode;
    }
  }  


  //print the components info into files
  makecolors();
  
  FILE *matlog;
  FILE *matviz;
  FILE *mattxt;
  char fname[300];
  int dumped;
  int ncomponents = 1; //the actual number of components, after merges
  char window1[20];
  char window2[20];

  for (tmp = 1; tmp <= numComponents; tmp++) {
    interNodesItr = interNodes[tmp];
    if (interNodesItr != NULL) {
      sprintf(fname, "component-%d.matrix", ncomponents);
      matlog = myfopen(fname, "w");
      sprintf(fname, "component-%d.viz", ncomponents);
      matviz = myfopen(fname, "w");
      sprintf(fname, "component-%d.txt", ncomponents);
      mattxt = myfopen(fname, "w");

      fprintf(matviz, "subgraph cluster%d{\n", ncomponents);
      fprintf(matviz, "\tlabel = \"Cluster #%d\";\n", ncomponents);
      dumped = 0;
      fwrite(&dumped, sizeof(int), 1 , matlog);
      fwrite(&WS, sizeof(int), 1, matlog);
      
      while(interNodesItr != NULL){
        dumped++;
        fwrite(&(interNodesItr->v1), sizeof(int), 1, matlog);
        fwrite(&(interNodesItr->v2), sizeof(int), 1, matlog);
        fwrite(&(interNodesItr->w), sizeof(int), 1, matlog);  
   
        decode_window(interNodesItr->v1, window1);
        decode_window(interNodesItr->v2, window2);
        fprintf(mattxt, "%s\n%s\n", window1, window2);
        //fprintf(matviz, "\t%s\t--\t%s [ color=\"%s\" ];\n", window1, window2, colors[ncomponents%10]);
        fprintf(matviz, "\t%s\t--\t%s;\n", window1, window2);

        interNodesItr = interNodesItr->next;
      }

      fprintf(stdout, "Component %d\tdumped:%d\n", ncomponents, dumped);
      fprintf(matviz, "}\n");
      rewind(matlog);
      fwrite(&dumped, sizeof(int), 1 , matlog);
      fclose(matlog);
      fclose(matviz);
      fclose(mattxt);
      /*
      if (dumped < 50){ // clean
        sprintf(fname, "component-%d.matrix", ncomponents);
        remove(fname);
        sprintf(fname, "component-%d.viz", ncomponents);
        remove(fname);
        sprintf(fname, "component-%d.txt", ncomponents);
        remove(fname);
	}*/
      ncomponents++;
    }
  }

  system("for i in `ls component*txt`; do sort -u $i > tmp; mv tmp $i; done");

  gettimeofday(&end, &tz);
  printf("\n***************\nComponents computation: %d microseconds.\n", ((int)(end.tv_sec*1000000+end.tv_usec)-(int)(start.tv_sec*1000000+start.tv_usec)));
  
}
//gozde modifications end



void alloc_node(struct Node **node){
  (*node) = (struct Node *) malloc (sizeof(struct Node));
  (*node)->next = NULL;
  (*node)->v1 = -1;
  (*node)->v2 = -1;
  (*node)->w = -1;
}


void alloc_component(struct Component **cmp){
  (*cmp) = (struct Component *) malloc (sizeof(struct Component));
  (*cmp)->node = NULL;
  (*cmp)->next = NULL;
}



void trimwins(char fname[], int hashtype, char trimwinfile[]){
  char outname[SEQ_LENGTH*2];
  char fname2[SEQ_LENGTH];
  FILE *winlog;
  FILE *trimfile;
  int i,j;
  int this;
  int cnt;
  int *trim;
  int breakthis;

  printf("Trimming...\n");

  strcpy(fname2, fname);
  fname2[strlen(fname)-7] = 0;
  sprintf(outname, "%s-%s.winlog", fname2, trimwinfile);
  trimfile = myfopen(trimwinfile, "r");
  
  cnt = 0;
  while (fscanf(trimfile, "%d", &this)  > 0) cnt++;
  rewind(trimfile);
  cnt = cnt*2;
  trim = (int *) malloc(sizeof(int)*cnt);

  wincnt = pow(4, WS);

  i = 0;
  while (fscanf(trimfile, "%d", &this)  > 0){ 
    trim[i++] = this;
    trim[i++] = revcomp_encoded(this);
  }

  winlog = fopen(outname, "w");

  fprintf(winlog,"%d\n",WS);

  sprintf(fname2, "%s-%s.h%d.dump", fname2, trimwinfile, hashtype);
  FILE *dumphist = fopen(fname2, "w");


  for(i=0;i<wincnt;i++){

    if (winlocs[i]==NULL)
      continue;
    breakthis = 0;

    for (j=0;j<cnt;j++)
      if (i == trim[j]){
	printf ("Trim %d\n", i);
	breakthis = 1;
	break;
      }
    
    if (breakthis) continue;

    dumpwins(winlog, i);
    fprintf(dumphist, "%d\t%d\n", i, windowcnt[i]);
  }

  fclose(winlog); fclose(dumphist);
		  
  free(trim);

}


void makecolors(void){
  int i;
  int maxcolor=10;
  colors = (char **)malloc(sizeof(char *)*maxcolor);

  for (i=0;i<maxcolor;i++)
    colors[i] = (char *)malloc(25);

  strcpy(colors[0], "red");
  strcpy(colors[1], "green");
  strcpy(colors[2], "blue");
  strcpy(colors[3], "magenta");
  strcpy(colors[4], "orange");
  strcpy(colors[5], "gold");
  strcpy(colors[6], "deeppink");
  strcpy(colors[7], "black");
  strcpy(colors[8], "darkslategary");
  strcpy(colors[9], "plum");


}
