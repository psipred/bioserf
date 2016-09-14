/*
 * PDB Sequence Extraction
 */

/* Version 0.01     by David T. Jones, April 2010 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define max(x,y) ((x)>(y) ? (x) : (y))
#define min(x,y) ((x)<(y) ? (x) : (y))
#define ran0() ((random()&32767)/(float)32768.0)

#define FALSE              0
#define TRUE               1

#define BIG 1000000

#define MAXNRES 3000
#define MAXNSEQ 500000

#define NTOKENS            28

/* A list of common PDB record types... */

#define HEADER             1
#define COMPND             2
#define SOURCE             3
#define AUTHOR             4
#define REVDAT             5
#define REMARK             6
#define SEQRES             7
#define FTNOTE             8
#define HET                9
#define FORMUL             10
#define HELIX              11
#define CRYST1             12
#define ORIGX1             13
#define ORIGX2             14
#define ORIGX3             15
#define SCALE1             16
#define SCALE2             17
#define SCALE3             18
#define ATOM               19
#define TER                20
#define HETATM             21
#define CONECT             22
#define ENDENT             23
#define JRNL               24
#define TURN               25
#define ENDMDL             26
#define EXPDTA             27
#define TITLE              28

void           *calloc(), *malloc();

struct SEQ
{
    int             nres;
    short           calpha, badseq;
    char            chainid, brkid[5], *seq;
    float           resol, rvalue;
}
db[MAXNSEQ];

/* Record names for decoding record types */
char           *tokstr[] =
{
    "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
    "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL",
    "HELIX", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
    "SCALE1", "SCALE2", "SCALE3", "ATOM", "TER",
    "HETATM", "CONECT", "END", "JRNL", "TURN", "ENDMDL",
    "EXPDTA", "TITLE"
};

/* Non-amino acid names to detect non-protein files */
char           *badnames = "M2G H2U PSU NAG NAM";

/* Residue name to allow conversion of a.a. name into numeric code */
char           *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "UNK", NULL
};

char            resnums[] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
};

char           *rescodes = "ARNDCQEGHILKMFPSTWYVBZX";

int             nseqs;
short         **idmat;

void
                fail(errstr)
     char           *errstr;
{
    printf("\n*** %s\n\n", errstr);
    exit(-1);
}

static char    *
                instr(ct, cs, l)
     char           *ct, *cs;
     int             l;
{
    for (; *ct; ct++)
	if (!strncmp(ct, cs, l))
	    return (ct);
    return (NULL);
}

void
                procfile(char *fname, char reqchain)
{
    int             aa, i,j, nunk, token, nres, isprot, natoms, ncalpha, nchains;
    char            brkid[5], keyword[8], buf[160], chainid, seq[MAXNRES], *vp;
    float           resol, rvalue;
    FILE           *ifp;

    printf("Processing : %s\n", fname);
    ifp = fopen(fname, "r");
    if (!ifp)
    {
	fprintf(stderr, "Failed to open PDB file!");
	return;
    }

    brkid[0] = '\0';
    resol = 99.9;
    rvalue = 9.999;
    nres = ncalpha = natoms = nunk = nchains = 0;
    isprot = TRUE;
    chainid = 0;
    while (!feof(ifp))
    {
	if (!fgets(buf, 160, ifp))
	    break;
	buf[6] = ' ';
	sscanf(buf, "%s", keyword);	/* Read the record name */

	if (!keyword[0])	/* Odd - there isn't a record name! Exit. */
	    break;

	token = 0;
	for (i = 1; i <= NTOKENS; i++)	/* Decode record type */
	    if (!strcmp(keyword, tokstr[i - 1]))
		token = i;

	switch (token)
	{
	case COMPND:
	case TITLE:
	    if (strstr(buf, "THEORETICAL MODEL"))
	    {
		puts("Theoretical model!");
		fclose(ifp);
		return;
	    }
	    break;

	case EXPDTA:
	    if (strstr(buf, "THEORETICAL") || strstr(buf, "NONE"))
	    {
		puts("Theoretical model!");
		fclose(ifp);
		return;
	    }
	    break;

	case HEADER:
	    strncpy(brkid, buf + 62, 4);
	    brkid[0] = tolower(brkid[0]);
	    brkid[1] = tolower(brkid[1]);
	    brkid[2] = tolower(brkid[2]);
	    brkid[3] = tolower(brkid[3]);
	    brkid[4] = '\0';
	    break;

	case REMARK:
	    if (!strncmp(buf + 8, " 2 RESOLUTION.", 14))
		if (sscanf(buf + 22, "%f", &resol) != 1)
		    puts("Bad RESOLUTION record!");
	    if (!strncmp(buf + 8, " 3   R VALUE", 12))
	    {
		if (!strstr(buf, "NULL"))
		{
		    vp = strstr(buf, " 0.");
		    if (rvalue > 1.0 && (!vp || sscanf(vp, "%f", &rvalue) != 1))
			puts("Bad R VALUE record!");
		}
	    }
	    break;

	case ENDENT:
	case ENDMDL:
	    if (isprot && nres >= 10)
	    {
		if (chainid == reqchain || reqchain == '?')
		{
		    if (nseqs > MAXNSEQ)
			fail("Too many chains!");
		    strcpy(db[nseqs].brkid, brkid);
		    db[nseqs].chainid = chainid;
		    db[nseqs].nres = nres;
		    db[nseqs].resol = resol;
		    db[nseqs].rvalue = rvalue;
		    db[nseqs].calpha = ncalpha > natoms / 2;
		    db[nseqs].badseq = nunk > nres / 4;
		    if ((db[nseqs].seq = malloc(nres)) == NULL)
			fail("Out of memory!");
		    memcpy(db[nseqs].seq, seq, nres);
		    nseqs++;
		    nchains++;
		}
	    }
	    break;
	    
	case ATOM:
	    if (chainid && buf[21] != chainid)
	    {
		if (isprot && nres > 10 && (!nchains || db[nseqs-1].nres != nres || memcmp(db[nseqs-1].seq, seq, nres)) && (nchains < 2 || db[nseqs-2].nres != nres || memcmp(db[nseqs-2].seq, seq, nres)))
		{
		    if (chainid == reqchain || reqchain == '?')
		    {
			if (nseqs > MAXNSEQ)
			    fail("Too many chains!");
			strcpy(db[nseqs].brkid, brkid);
			db[nseqs].chainid = chainid;
			db[nseqs].nres = nres;
			db[nseqs].resol = resol;
			db[nseqs].rvalue = rvalue;
			db[nseqs].calpha = ncalpha > natoms / 2;
			db[nseqs].badseq = nunk > nres / 4;
			if ((db[nseqs].seq = malloc(nres)) == NULL)
			    fail("Out of memory!");
			memcpy(db[nseqs].seq, seq, nres);
			nseqs++;
			nchains++;
		    }
		}
		nres = ncalpha = natoms = nunk = 0;
		isprot = TRUE;
	    }
	    chainid = buf[21];
	    if (isprot && buf[17] == ' ')
	    {
		puts("Not a protein chain!");
		isprot = FALSE;
	    }
	    if (!strncmp(buf + 12, " CA ", 4) && (buf[16] == ' ' || buf[16] == 'A') && (buf[20] == ' ' || isdigit(buf[20])))
	    {
		for (aa = 0; rnames[aa] != NULL; aa++)
		    if (!strncmp(buf + 17, rnames[aa], 3))
			break;
		if (aa >= 20)
		    nunk++;
		if (rnames[aa] == NULL)
		    seq[nres++] = 22;
		else
		    seq[nres++] = aa;
		ncalpha++;
	    }		
	    natoms++;
	    break;

	default:		/* Ignore all other types in this version */
	    break;
	}
	
	if (token == ENDENT || token == ENDMDL)
	    break;
    }
    
    fclose(ifp);
}

int 
nwscore(seq1, seq2, len1, len2, gap_pen)
    char           *seq1, *seq2;
    int             len1, len2, gap_pen;
{
    int             diag, col, row, maxcol, maxrows[MAXNRES], toprows[MAXNRES], i, j, maxscore = -BIG;
    int             maxcgap, maxrgap[MAXNRES], maxid, now = 0, last = 1;
    int             mat[2][MAXNRES], gap[2][MAXNRES];

    for (i = 0; i < MAXNRES; i++)
	maxrows[i] = -BIG;

    for (j = len2 - 1; j >= 0; j--)
    {
	maxcol = -BIG;

	for (i = len1 - 1; i >= 0; i--)
	{
	    gap[now][i] = mat[now][i] = (seq1[i] == seq2[j]) ? 1 : 0;
	    if (j != len2 - 1 && i != len1 - 1)
	    {
		diag = mat[last][i + 1];
		col = maxcol - gap_pen;
		row = maxrows[i] - gap_pen;

		if (diag >= col && diag >= row)
		{
		    mat[now][i] += diag;
		    gap[now][i] += gap[last][i + 1];
		}
		else
		{
		    if (row > col)
		    {
			mat[now][i] += row;
			gap[now][i] += maxrgap[i+1];
		    }
		    else
		    {
			mat[now][i] += col;
			gap[now][i] += maxcgap;
		    }
		}

		if (mat[now][i] > maxscore)
		{
		    maxscore = mat[now][i];
		    maxid = gap[now][i];
		}

		if (diag > maxrows[i])
		{
		    maxrows[i] = diag;
		    maxrgap[i] = gap[last][i+1];
		    toprows[i] = j + 1;
		}
		if (diag > maxcol)
		{
		    maxcol = diag;
		    maxcgap = gap[last][i+1];
		}
	    }
	}
	now = !now;
	last = !last;
    }

    return maxid;
}

/* Perform random shuffle of sequence */
shuffle(s, len)
     char           *s;
     int             len;
{
    int             i, r;
    char            temp;

    for (i = 0; i < len; i++)
    {
	r = ran0() * len;
	temp = s[i];
	s[i] = s[r];
	s[r] = temp;
    }
}

/* Check significance of NW score */
int             is_sig(score, seq1, seq2, len1, len2, gap_pen)
     int             score, len1, len2, gap_pen;
     char           *seq1, *seq2;
{
    int             i, minlen;
    char            rseq1[MAXNRES], rseq2[MAXNRES];
    float           perid;

    minlen = min(len1, len2);
    if (minlen >= 50)
    {
	perid = (float) 100.0 *score / minlen;

	/* Check if we are in the Twilight Zone... */
	if (perid < (float) 20.0)
	    return TRUE;
	if (perid >= (float) 30.0)
	    return TRUE;
    }

    return ((float) minlen / max(len1, len2) > 0.1);
}

/* Compare two sequence db entries */
int
                cmpdb(const void *s1, const void *s2)
{
    int             c;

    if (((struct SEQ *)s1)->resol < ((struct SEQ *)s2)->resol)
	return (-1);
    if (((struct SEQ *)s1)->resol > ((struct SEQ *)s2)->resol)
	return (1);
    if (((struct SEQ *)s1)->rvalue < ((struct SEQ *)s2)->rvalue)
	return (-1);
    if (((struct SEQ *)s1)->rvalue > ((struct SEQ *)s2)->rvalue)
	return (1);

    c = strcmp(((struct SEQ *)s1)->brkid, ((struct SEQ *)s2)->brkid);
    return (c) ? c : ((struct SEQ *)s1)->chainid - ((struct SEQ *)s2)->chainid;
}

main(argc, argv)
     int             argc;
     char           *argv[];
{
    int             i, j, sc, hsc;
    char            buf[512], fname[160], pathname[512], chainid;
    float           prob;
    FILE           *ofp, *bfp, *lfp;

    if (argc != 4)
    {
	printf("Usage: pdbhomol pdbpath pdblist outputfile\n\n");
	exit(-1);
    }

    srandom(1);

    lfp = fopen(argv[2], "r");
    if (!lfp)
	fail("Cannot open list file!");

    while (!feof(lfp))
    {
	if (fscanf(lfp, "%s", fname) != 1)
	    break;
	chainid = '?';
	sprintf(pathname, "%s/%s", argv[1], fname);
	j = nseqs;
	procfile(pathname, chainid);
	printf("%d chains from %s.\n", nseqs-j, pathname);
    }

    fclose(lfp);

    printf("%d chains.\n", nseqs);
    fflush(stdout);

    qsort(db, nseqs, sizeof(struct SEQ), cmpdb);

    ofp = fopen(argv[3], "w");
    
    if (!ofp)
	fail("Cannot open output file!");
    
    for (i = 0; i < nseqs; i++)
	if (!db[i].calpha && !db[i].badseq && (db[i].resol < 4.0 || db[i].resol > 99.0) && db[i].nres >= 20)
	{
	    fprintf(ofp, ">%s%c0 %4d %5.2f %5.3f\n", db[i].brkid, db[i].chainid, db[i].nres, db[i].resol, db[i].rvalue);
	    for (j = 0; j < db[i].nres; j++)
		fputc(rescodes[db[i].seq[j]], ofp);
	    fprintf(ofp, "\n");
	}

    fclose(ofp);
}
