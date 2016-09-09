/* Protein Chain Shotgun Superposition Program - by David Jones, April 2004 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQLEN 5000
#define MAXNSTRUC 10000

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#define BIG ((float)1.0e10)
#define SMALL ((float)1.0e-3)
#define ZERO ((float)0.0)
#define HALF ((float)0.5)
#define ONE ((float)1.0)
#define TWO ((float)2.0)
#define THREE ((float)3.0)
#define FOUR ((float)4.0)
#define PI ((float)3.1415927)

#define NTOKENS 26
#define DEFR 1.80

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

#define MEMGRAIN	   100

#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define veczero(v) memset(v, 0, sizeof(v))
#define veccopy(a,b) (a[0]=b[0],a[1]=b[1],a[2]=b[2])
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) (a[0]=b[0]*c,a[1]=b[1]*c,a[2]=b[2]*c)
#define ran0() ((rand()&32767)/(float)32768.0)
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define SQR(x) ((x)*(x))
#define dist(a,b) sqrt(SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))

extern void    *
                calloc(), *malloc(), *realloc();

typedef float   Transform[4][4];
typedef float   Point[3];

typedef struct
{
    short           length;
    short           posn_a[MAXSEQLEN], posn_b[MAXSEQLEN];
}
ALNS;

struct pdbatm
{
    Point           pos;
    float           radius;	/* Radius of sphere */
    int             resnum, aac, nequiv;
    char            ispolar, isback, donors, acceptors, chain;
    char            ndon, nacc, donflag, accflag;
    char            atmnam[5];
}              *atoms[MAXNSTRUC];

float mqapscore[MAXNSTRUC];

char            pdbfn[80], csdfn[80], logfn[80], keyword[40], buf[8192];

/* Record names for decoding record types */
char           *tokstr[] =
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL",
 "HELIX", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM", "TER",
 "HETATM", "CONECT", "END", "JRNL", "TURN", "ENDMDL"
};

/* Residue name to allow conversion of a.a. name into numeric code */
char           *rnames[] =
{
 "ALA", "ARG", "ASN", "ASP", "CYS",
 "GLN", "GLU", "GLY", "HIS", "ILE",
 "LEU", "LYS", "MET", "PHE", "PRO",
 "SER", "THR", "TRP", "TYR", "VAL",
 "ASX", "GLX", "UNK"
};

enum RESCODE
{
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys,
    Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val,
    Asx, Glx, Unk
};

struct distrec
{
    short           a, b;
}              *distmat, **sortlist;

int             nmodels, natoms[MAXNSTRUC];
float         **dmat;

int pat[MAXSEQLEN][MAXSEQLEN];


/***************************************************************************/

void
                fail(errstr)
    char           *errstr;
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = (void **) malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	for (i = 0; i < rows; i++)
	    if ((p[i] = calloc(columns, size)) == NULL)
		fail("allocmat: calloc [][] failed!");
    }
    else
	for (i = 0; i < rows; i++)
	    if ((p[i] = malloc(columns * size)) == NULL)
		fail("allocmat: malloc [][] failed!");

    return p;
}

void            freemat(void *p, int rows)
{
    int             i;

    for (i = rows - 1; i >= 0; i--)
	free(((void **) p)[i]);
    free(p);
}

/* Apply a Transform matrix to a point */
void
                transform_point(transform, p, tp)
    Transform       transform;	/* transform to apply to the point */
    Point           p;		/* the point to transform */
    Point           tp;		/* the returned point after transformation */
{
    int             i, j;
    Point           temp;

    temp[0] = p[0] + transform[0][3];
    temp[1] = p[1] + transform[1][3];
    temp[2] = p[2] + transform[2][3];
    tp[0] = dotprod(transform[0], temp);
    tp[1] = dotprod(transform[1], temp);
    tp[2] = dotprod(transform[2], temp);

#if 0
    printf("%g %g %g\n%g %g %g\n\n", p[0], p[1], p[2], tp[0], tp[1], tp[2]);
#endif
}


#define ROTATE(a,i,j,k,l) (g=a[i][j],h=a[k][l],a[i][j]=g-s*(h+g*tau),a[k][l]=h+s*(g-h*tau))

int             jacobi(float a[6][6], float d[6], float v[6][6], int *nrot)
{
    int             j, iq, ip, i;
    float           tresh, theta, tau, t, sm, s, h, g, c, b[6], z[6];

    for (ip = 0; ip < 6; ip++)
    {
	for (iq = 0; iq < 6; iq++)
	    v[ip][iq] = ZERO;
	v[ip][ip] = 1.0;
    }

    for (ip = 0; ip < 6; ip++)
    {
	b[ip] = d[ip] = a[ip][ip];
	z[ip] = ZERO;
    }

    *nrot = 0;
    for (i = 0; i < 50; i++)
    {
	sm = ZERO;
	for (ip = 0; ip < 5; ip++)
	{
	    for (iq = ip + 1; iq < 6; iq++)
		sm += fabs(a[ip][iq]);
	}
	if (sm == ZERO)
	    return 0;
	if (i < 3)
	    tresh = 0.2 * sm / 36;
	else
	    tresh = ZERO;
	for (ip = 0; ip < 5; ip++)
	{
	    for (iq = ip + 1; iq < 6; iq++)
	    {
		g = 100.0 * fabs(a[ip][iq]);
		if (i > 3 && fabs(d[ip]) + g == fabs(d[ip])
		    && fabs(d[iq]) + g == fabs(d[iq]))
		    a[ip][iq] = ZERO;
		else
		    if (fabs(a[ip][iq]) > tresh)
		    {
			h = d[iq] - d[ip];
			if (fabs(h) + g == fabs(h))
			    t = (a[ip][iq]) / h;
			else
			{
			    theta = 0.5 * h / (a[ip][iq]);
			    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
			    if (theta < ZERO)
				t = -t;
			}
			c = 1.0 / sqrt(1 + t * t);
			s = t * c;
			tau = s / (1.0 + c);
			h = t * a[ip][iq];
			z[ip] -= h;
			z[iq] += h;
			d[ip] -= h;
			d[iq] += h;
			a[ip][iq] = ZERO;
			for (j = 0; j <= ip - 1; j++)
			{
			    ROTATE(a, j, ip, j, iq);
			}
			for (j = ip + 1; j <= iq - 1; j++)
			{
			    ROTATE(a, ip, j, j, iq);
			}
			for (j = iq + 1; j < 6; j++)
			{
			    ROTATE(a, ip, j, iq, j);
			}
			for (j = 0; j < 6; j++)
			{
			    ROTATE(v, j, ip, j, iq);
			}
			++(*nrot);
		    }
	    }
	}
	for (ip = 0; ip < 6; ip++)
	{
	    b[ip] += z[ip];
	    d[ip] = b[ip];
	    z[ip] = ZERO;
	}
    }
    return 1;
}

#undef ROTATE

/*
 * function eigsrt 
 *
 * Given the eigenvalues d[n] and eignevectors v[n][n] as output from jacobi
 * this routine sourts the eigenvalues into descending order and rearranges
 * the columns of v correspondingly 
 */
void            eigsrt(float d[6], float v[6][6])
{
    int             k, j, i;
    float           p;

    for (i = 0; i < 5; i++)
    {
	p = d[k = i];
	for (j = i + 1; j < 6; j++)
	    if (d[j] >= p)
		p = d[k = j];
	if (k != i)
	{
	    d[k] = d[i];
	    d[i] = p;
	    for (j = 0; j < 6; j++)
	    {
		p = v[j][i];
		v[j][i] = v[j][k];
		v[j][k] = p;
	    }
	}
    }
}

int
                lsq_fit(float u[3][3], Transform R)
{
    float           du, omega[6][6], vom[6][6];
    float           dom[6], root2, h[3][3], k[3][3], sign;
    int             i, j, l, rot;

    /* Constant */
    root2 = sqrt(2.0);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = ZERO;

    for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    omega[i][j] = 0;

    /* Calculate determinant of U */
    du = u[0][0] * u[1][1] * u[2][2] - u[0][0] * u[1][2] * u[2][1]
	- u[0][1] * u[1][0] * u[2][2] + u[0][1] * u[1][2] * u[2][0]
	+ u[0][2] * u[1][0] * u[2][1] - u[0][2] * u[1][1] * u[2][0];

    /* If determinant is zero return */
    if (fabs(du) < 1.0e-5)
	return TRUE;

    /* Make metric matrix omega */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    omega[i][j] = ZERO;
	    omega[i + 3][j + 3] = ZERO;
	    omega[i + 3][j] = u[j][i];
	    omega[i][j + 3] = u[i][j];
	}

    /* Diagonalise matrix */
    if (jacobi(omega, dom, vom, &rot))
	return TRUE;

    /* Sort by eigenvalues */
    eigsrt(dom, vom);

    /* Check for degeneracy */
    if (du <= ZERO && fabs(dom[2] - dom[5]) < 1.0e-5)
	return TRUE;

    /* Determine h and k */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    h[i][j] = root2 * vom[i][j];
	    k[i][j] = root2 * vom[i + 3][j];
	}

    sign = h[0][0] * h[1][1] * h[2][2] - h[0][0] * h[1][2] * h[2][1]
	- h[0][1] * h[1][0] * h[2][2] + h[0][1] * h[1][2] * h[2][0]
	+ h[0][2] * h[1][0] * h[2][1] - h[0][2] * h[1][1] * h[2][0];

    if (sign <= ZERO)
	for (i = 0; i < 3; i++)
	{
	    h[i][2] = -h[i][2];
	    k[i][2] = -k[i][2];
	}

    /* Determine rotation matrix */
    du /= fabs(du);
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = k[j][0] * h[i][0] + k[j][1] * h[i][1] + du * k[j][2] * h[i][2];
    R[0][3] = R[1][3] = R[2][3] = R[3][0] = R[3][1] = R[3][2] = ZERO;
    R[3][3] = ONE;

    return FALSE;
}


/*
 * Trace back highest scoring path 
 */
void
                trace(short *posa, short *posb, int mati, int matj,
		      int pati, int patj, int lasti, int lastj, short *n)
{
    int             pij = pat[pati][patj], i, j;

    for (i = lasti + 1; i < mati; i++)
    {
	*(++posa) = i;
	*(++posb) = 0;
	(*n)++;
    }
    for (j = lastj + 1; j < matj; j++)
    {
	*(++posa) = 0;
	*(++posb) = j;
	(*n)++;
    }
    *(++posa) = mati;
    *(++posb) = matj;
    (*n)++;

    if (!pij)
	return;

    if (pij == 1)
	trace(posa, posb, mati + 1, matj + 1, pati + 1, patj + 1,
	      mati, matj, n);
    if (pij < 1)
	trace(posa, posb, mati + 1, matj - pij, pati + 1, patj - pij,
	      mati, matj, n);
    if (pij > 1)
	trace(posa, posb, mati + pij, matj + 1, pati + pij, patj + 1,
	      mati, matj, n);
}

int
                seqscore(const struct pdbatm *atom1, const struct pdbatm *atom2, const int gappen, const int gapext, ALNS * aln, const int seq1len, const int seq2len)
{
    short          *posa, *posb;
    int             trace_back = (aln != NULL);
    int             now = 0, last = 1;
    int             pati, patj, mati, matj, i, j;
    int             toprows[MAXSEQLEN + 1], topcol, toprow;
    int             maxrows[MAXSEQLEN + 1], maxscore, maxflg = FALSE, maxcol,
                    maxrow, diag, row, col;
    int             mat[2][MAXSEQLEN + 1];

    for (i = 1; i <= seq1len; i++)
    {
	mat[0][i] = mat[1][i] = maxrows[i] = -1000000;
	pat[i][seq2len] = toprows[i] = 0;
    }

    for (j = seq2len; j > 0; j--)
    {
	maxcol = -1000000;
	topcol = 0;

	if (trace_back)
	    pat[seq1len][j] = 0;

	for (i = seq1len; i > 0; i--)
	{
	    /* Get matrix element */
	    mat[now][i] = atom2[j - 1].aac == atom1[i - 1].aac ? 1 : -1;

	    if (j != seq2len && i != seq1len)
	    {
		diag = mat[last][i + 1];
		maxrow = maxrows[i];
		toprow = toprows[i];

		switch (gapext)
		{
		case 1:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) + 1;
			row = maxrow - gappen - (toprow - j) + 1;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		    break;
		case 2:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) - (topcol - i) + 1;
			row = maxrow - gappen - (toprow - j) - (toprow - j) + 1;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		    break;
		case 0:
		    col = maxcol - gappen;
		    row = maxrow - gappen;
		    break;
		default:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) * gapext + gapext;
			row = maxrow - gappen - (toprow - j) * gapext + gapext;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		}

		if (diag >= col && diag >= row)
		{
		    mat[now][i] += diag;
		    if (trace_back)
			pat[i][j] = 1;
		}
		else
		{
		    if (row > col)
		    {
			mat[now][i] += row;
			if (trace_back)
			    pat[i][j] = -(toprow - j);
		    }
		    else
		    {
			mat[now][i] += col;
			if (trace_back)
			    pat[i][j] = topcol - i;
		    }
		}

		if (diag > maxrows[i])
		{
		    maxrows[i] = diag;
		    toprows[i] = j + 1;
		}
		if (diag > maxcol)
		{
		    maxcol = diag;
		    topcol = i + 1;
		}

		if ((i == 1 || j == 1) && (!maxflg || mat[now][i] > maxscore))
		{
		    maxflg = TRUE;
		    maxscore = mat[now][i];
		    if (trace_back)
		    {
			pati = i;
			patj = j;
			mati = matj = 1;
			if (i == 1)
			    matj = j;
			if (j == 1)
			    mati = i;
		    }
		}
	    }
	}
	now = !now;
	last = !last;
    }

    if (!trace_back)
	return (maxscore);

    posa = aln->posn_a;
    posb = aln->posn_b;
    aln->length = 0;
    trace(posa, posb, mati, matj, pati, patj, 0, 0, &aln->length);
    posa += aln->length;
    posb += aln->length;
    if (*posa == seq1len)
    {
	for (i = *(posb) + 1; i <= seq2len; i++)
	{
	    *(++posa) = 0;
	    *(++posb) = i;
	    (aln->length)++;
	}
	if (aln->length > MAXSEQLEN)
	    puts("score : max. align length exceeded!");
	return (maxscore);
    }
    if (*posb == seq2len)
	for (i = *(posa) + 1; i <= seq1len; i++)
	{
	    *(++posb) = 0;
	    *(++posa) = i;
	    (aln->length)++;
	}
    if (aln->length > MAXSEQLEN)
	puts("score : max. align length exceeded!");

    return (maxscore);
}

/* Get atomic coordinates from ATOM records */
void
                getcoord(x, y, z, chain, n, aacode, atmnam)
    float          *x, *y, *z;
    int            *n, *aacode;
    char           *atmnam, *chain;
{
    int             i, atomn;

    strncpy(atmnam, buf + 12, 4);
    atmnam[4] = '\0';
    if (atmnam[2] == ' ')
	atmnam[3] = ' ';
    sscanf(buf + 6, "%d", &atomn);
    *chain = buf[21];

    sscanf(buf + 22, "%d%f%f%f", n, x, y, z);
    for (i = 0; i < 22; i++)
	if (!strncmp(rnames[i], buf + 17, 3))
	    break;
    *aacode = i;
}

int
                read_atoms(FILE *ifp, int *natoms, struct pdbatm **atmptr)
{
    int             i, token=ENDENT, namino, aac, resnum = 0, blksize, read_end = 0;
    float           x, y, z, u[3][3];
    char            atmnam[5], chain = '?';
    struct pdbatm  *atom = *atmptr;

    blksize = 10 * MEMGRAIN;
    if (!(atom = malloc(blksize * sizeof(struct pdbatm))))
	fail("read_atoms : Out of memory!");

    while (!feof(ifp))
    {
	if (!fgets(buf, 8192, ifp))
	    break;

	if (!isalpha(buf[0]))
	    continue;

	/* Read the record name */
	if (sscanf(buf, "%s", keyword) != 1)
	    continue;

	if (!keyword[0])	/* Odd - there isn't a record name! Skip. */
	    continue;

	token = 0;
	for (i = 1; i <= NTOKENS; i++)	/* Decode record type */
	    if (!strcmp(keyword, tokstr[i - 1]))
		token = i;

	if (token == ENDENT || token == ENDMDL)
	{
	    read_end = 1;
	    break;
	}
	
	switch (token)
	{
	case ATOM:
	    if (buf[16] != ' ' && buf[16] != 'A')
		continue;
	    buf[26] = ' ';
	    if (*natoms >= blksize)
	    {
		blksize += MEMGRAIN;
		atom = realloc(atom, blksize * sizeof(struct pdbatm));
		if (!atom)
		    fail("read_atoms : Out of Memory!");
	    }
	    getcoord(&x, &y, &z, &chain, &namino, &aac, atmnam);
	    if (strncmp(atmnam, " CA", 3))
		continue;
	    resnum++;
	    atom[*natoms].pos[0] = x;
	    atom[*natoms].pos[1] = y;
	    atom[*natoms].pos[2] = z;
	    atom[*natoms].chain = chain;
	    atom[*natoms].resnum = resnum;
	    atom[*natoms].aac = aac;
	    atom[*natoms].isback = i > 0 && i < 4;
	    strcpy(atom[*natoms].atmnam, atmnam);
	    ++*natoms;
	    break;
	default:		/* Ignore all other types in this version */
	    break;
	}
    }

    *atmptr = atom;

    if (!read_end)
    {
	*natoms = 0;
	return 1;
    }
    

    return 0;
}

/* Model score */
struct mscore
{
    float             score;
    unsigned int modelnum, skipflg;
};

struct mscore scoretab[MAXNSTRUC];

int             cmpsc(const void *a, const void *b)
{
    if (((struct mscore *) a)->score < ((struct mscore *) b)->score)
	return 1;
    else if (((struct mscore *) a)->score > ((struct mscore *) b)->score)
	return -1;

    return 0;
}

/* Convert AA letter to numeric code (0-22) */
int
  aanum(int ch)
{
    static int aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}

/* This routine will read in one sequence from a database file. The
   sequence can be in any of the supported formats. Returns length
   of sequence.
*/
int
  getseq(char *dbname, char *dseq, FILE * lfil)
{
    int i, j, len;
    short badln, fformat;
    enum
    {
	unknown, embl, genbank, staden, fastp, codata, owl, intelgen, gcg
    };
    char temp[8192], split;
    int offset;

    offset = j = 0;

    if (!fgets(temp, 8192, lfil))
	return (-1);
    if (strstr(temp, "of:") != NULL && strstr(temp, "check:") != NULL)
	fformat = gcg;
    else if ((temp[0] == '<') && (temp[19] == '>'))
	fformat = staden;
    else if (strncmp(temp, "ID   ", 5) == 0)
	fformat = embl;
    else if (strncmp(temp, "LOCUS     ", 10) == 0)
	fformat = genbank;
    else if (strncmp(temp, "ENTRY", 5) == 0)
	fformat = codata;
    else if (temp[0] == ';')
	fformat = intelgen;
    else if (temp[0] == '>' && (temp[1] == '>' || temp[3] == ';'))
	fformat = owl;
    else if (temp[0] == '>')
	fformat = fastp;
    else
	fformat = unknown;

    switch (fformat)
    { 
     case gcg:
	sscanf(strstr(temp, "of:")+3, "%s", dbname);
	while (strstr(temp, "..") == NULL)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;
   case embl:
	strncpy(dbname, temp + 5, 70);
	while (temp[0] != ' ')
	    fgets(temp, 8192, lfil);
	break;

    case genbank:
	while (strncmp(temp, "ORIGIN", 6) != 0)
	{
	    fgets(temp, 8192, lfil);
	    if (strncmp(temp, "DEFINITION", 10) == 0)
		strncpy(dbname, temp + 12, 70);
	}
	fgets(temp, 8192, lfil);
	break;

    case codata:
	strncpy(dbname, temp + 6, 70);
	while (strncmp(temp, "SEQUENCE", 8) != 0)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    case owl:
	fgets(temp, 8192, lfil);
	strncpy(dbname, temp, 70);
	fgets(temp, 8192, lfil);
	break;

    case fastp:
	strncpy(dbname, temp + 1, 70);
	fgets(temp, 8192, lfil);
	break;

    case staden:
	strncpy(dbname, temp + 1, 18);
	offset = 20;
	break;

    case intelgen:
	while (*temp == ';')
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    default:
	do
	{
	    len = strlen(temp);
	    for (badln = i = 0; i < len; i++)
		if (islower(temp[i]) || temp[i] == 'J' || temp[i] == 'O' || temp[i] == 'U')
		{
		    badln = TRUE;
		    break;
		}
	    if (badln && !fgets(temp, 8192, lfil))
		return (-1);
	}
	while (badln);
	strcpy(dbname, "<NO NAME>");
	break;
    }

    if (dbname[(len = strlen(dbname)) - 1] == '\n')
	dbname[--len] = '\0';
    if (len >= 70)
	dbname[70] = '\0';

    for (;;)
    {
	if (!strncmp(temp, "//", 2))
	    break;
	len = strlen(temp);
	for (i = offset; i < len && j < MAXSEQLEN; i++)
	{
	    split = islower(temp[i]) ? toupper(temp[i]) : temp[i];
	    if (split == '@' || (fformat == owl && split == '*'))
	    {
		dseq[j] = '\0';
		while (fgets(temp, 8192, lfil));
		return (j);
	    }
	    if (isalpha(split))
		dseq[j++] = split;
	    else if (temp[i] == '\n')
		break;
	}
	if (staden)
	    offset = 0;
	if (!fgets(temp, 8192, lfil))
	    break;
    }

    if (j == MAXSEQLEN)
	printf("\nWARNING: sequence %s over %d long; truncated!\n",
	       dbname, MAXSEQLEN);

    dseq[j] = '\0';
    return (j);
}


main(argc, argv)
    int             argc;
    char          **argv;
{
    int             i, j, k, ii, jj, kk, l, n, nmax = 0, at1, at2, nequiv, **neqmat, seqlen=0, ntemplates, nmqap;
    int             blksize, hashval, moda, modb, maxclusz, repnum, nclust, naln, mqapok[MAXNSTRUC];
    float           x, y, z, d, r, rmsd, u[3][3], **tmscmat, tmsc, matchsum, cutoff, jury[MAXNSTRUC];
    int             mineqv = 56, nsup;
    int             **incompatmat, tlist[MAXNSTRUC];
    long            filepos[MAXNSTRUC];
    float           eqvdist = 5.93, maxrms = 1000.0, d0sq, rsq;
    struct pdbatm   targstr[MAXSEQLEN];
    FILE           *ifp, *ofp, *afp;
    Point           new, CG_a, CG_b;
    Transform       fr_xf;
    ALNS tplt;
    char fname[160], seq[MAXSEQLEN], id[160];
    const char *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";
    
    if (argc != 4)
	fail("usage : tmjury3d_filter target-seq ensemble.pdb output.pdb");

    ifp = fopen(argv[1], "r");	/* Open seq file in TEXT mode */
    if (!ifp)
    {
	printf("Sequence file error!\n");
	exit(-1);
    }

    seqlen = getseq(id, seq, ifp);

    fclose(ifp);

    for (i=0; i<seqlen; i++)
	targstr[i].aac = aanum(seq[i]);

    ifp = fopen(argv[2], "r");	/* Open PDB file in TEXT mode */
    if (!ifp)
    {
	printf("PDB file error!\n");
	exit(-1);
    }

    /* Read models */
    for (i=nmodels=0; nmodels<MAXNSTRUC; i++)
    {
	filepos[nmodels] = ftell(ifp);
	if (read_atoms(ifp, natoms+nmodels, atoms+nmodels))
	    break;
//	printf("%d %d\n", i, natoms[nmodels]);
	if (natoms[nmodels])
	    nmodels++;
	else
	{
	    printf("Warning - model %d has zero length!!\n", nmodels+1);
	    i--;
	}
    }

    printf("%d models read from PDB file\n", nmodels);

    /* Calculate all pairwise equivalences */

    incompatmat = allocmat(nmodels, nmodels, sizeof(int), TRUE);
	
    for (moda=0; moda<nmodels; moda++)
    {
	jury[moda] = 0;

	if (natoms[moda] > seqlen)
	    fail("Model longer than target sequence!");
	
	if (natoms[moda] > 0)
	    for (i=0; i<natoms[moda]; i++)
		atoms[moda][i].nequiv = 0;
    }

    for (moda=0; moda<nmodels; moda++)
	for (modb=moda+1; modb<nmodels; modb++)
	    if (natoms[moda] > 10 && natoms[modb] > 10)
	    {
		if (seqscore(atoms[moda], atoms[modb], 1, 0, &tplt, natoms[moda], natoms[modb]) < 10)
		    continue;
		
		/* Calculate centroids */
		veczero(CG_a);
		veczero(CG_b);
		for (naln = 0, j = 1; j <= tplt.length; j++)
		    if (tplt.posn_a[j] > 0 && tplt.posn_b[j] > 0)
		    {
			vecadd(CG_a, CG_a, atoms[moda][tplt.posn_a[j]-1].pos);
			vecadd(CG_b, CG_b, atoms[modb][tplt.posn_b[j]-1].pos);
			naln++;
		    }
		
		vecscale(CG_a, CG_a, ONE / naln);
		vecscale(CG_b, CG_b, ONE / naln);
		
		for (i = 0; i < natoms[moda]; i++)
		    vecsub(atoms[moda][i].pos, atoms[moda][i].pos, CG_a);
		for (j = 0; j < natoms[modb]; j++)
		    vecsub(atoms[modb][j].pos, atoms[modb][j].pos, CG_b);
		
		/* Calculate U */
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++)
			u[i][j] = ZERO;
		
		for (j = 1; j <= tplt.length; j++)
		    if (tplt.posn_a[j] > 0 && tplt.posn_b[j] > 0)
			for (k = 0; k < 3; k++)
			    for (l = 0; l < 3; l++)
				u[k][l] += atoms[moda][tplt.posn_a[j]-1].pos[k] * atoms[modb][tplt.posn_b[j]-1].pos[l];
		
		if (lsq_fit(u, fr_xf))
		    continue; /* LSQ fit failed - all we can do is move to the next pair! */
		
		for (j = 0; j < natoms[modb]; j++)
		{
		    transform_point(fr_xf, atoms[modb][j].pos, new);
		    veccopy(atoms[modb][j].pos, new);
		}
		
		for (nequiv = 0, i = 1; i <= tplt.length; i++)
		    if (tplt.posn_a[i] > 0 && tplt.posn_b[i] > 0)
			nequiv++;

		d0sq = SQR(1.24 * pow(nequiv-15.0, 1.0/3.0) - 1.8);
		
		for (nsup = tmsc = 0, i = 1; i <= tplt.length; i++)
		    if (tplt.posn_a[i] > 0 && tplt.posn_b[i] > 0)
		    {
			rsq = distsq(atoms[moda][tplt.posn_a[i]-1].pos, atoms[modb][tplt.posn_b[i]-1].pos);
			if (rsq < 25.0F)
			    nsup++;
			tmsc += 1.0F / (1.0F + rsq/d0sq);
		    }

		tmsc /= nequiv;

		if (nequiv > 50 && tmsc < 0.6)
		    incompatmat[moda][modb] = incompatmat[modb][moda] = 1;
		
		if (tmsc < 0.99)
		{
		    jury[moda] += tmsc * natoms[moda] / seqlen;
		    jury[modb] += tmsc * natoms[modb] / seqlen;
		}
	    }

    printf("nmodels = %d\n", nmodels);

    ofp = fopen(argv[3], "w");

    if (ofp == NULL)
	fail("Cannot output PDB file!");

    for (moda=0; moda<nmodels; moda++)
    {
	for (modb=moda-1; modb>=0; modb--)
	    if (incompatmat[moda][modb])
		break;
	
	if (modb < 0)
	{	
	    fseek(ifp, filepos[moda], SEEK_SET);
	    while (!feof(ifp))
	    {
		if (!fgets(buf, 1024, ifp))
		    break;
		if (strncmp(buf, "ATOM", 4) && strncmp(buf, "END", 3))
		    continue;
		fputs(buf, ofp);
		if (!strncmp(buf, "END", 3))
		    break;
	    }
	}
    }

    fclose(ofp);

    return 0;
}
