
#ifndef _SS_C50Classifier_H
#define _SS_C50Interface_H


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <vector>
#include <string.h>
#include <string>
#ifdef WIN32
#include <windows.h>
#endif


/*************************************************************************/
/*									 */
/*		Constants, macros etc.					 */
/*									 */
/*************************************************************************/


#define	 SEE5

#define	 Nil	   0		/* null pointer */
#define	 false	   0
#define	 true	   1
#define	 None	   -1
#define	 Epsilon   1E-4

#define  EXCLUDE   1		/* special attribute status: do not use */
#define  SKIP	   2		/* do not use in classifiers */
#define  DISCRETE  4		/* ditto: collect values as data read */
#define  ORDERED   8		/* ditto: ordered discrete values */
#define  DATEVAL   16		/* ditto: YYYY/MM/DD or YYYY-MM-DD */
#define  STIMEVAL  32		/* ditto: HH:MM:SS */
#define	 TSTMPVAL  64		/* date time */

/* unknown and N/A values are represented by
   unlikely floating-point numbers
   (octal 01600000000 and 01) */
#define	 UNKNOWN   01600000000	/* 1.5777218104420236e-30 */
#define	 NA	   01		/* 1.4012984643248171e-45 */

#define	 BrDiscr   1
#define	 BrThresh  2
#define	 BrSubset  3

#define  Alloc(N,T)		(T *) Pmalloc((N)*sizeof(T))
#define  AllocZero(N,T)		(T *) Pcalloc(N, sizeof(T))
#define  Realloc(V,N,T)		V = (T *) Prealloc(V, (N)*sizeof(T))

#define	 Max(a,b)		((a)>(b) ? (a) : (b))
#define	 Min(a,b)		((a)<(b) ? (a) : (b))

#define	 Bit(b)			(1 << (b))
#define	 In(b,s)		((s[(b) >> 3]) & Bit((b) & 07))
#define	 SetBit(b,s)		(s[(b) >> 3] |= Bit((b) & 07))

#define	 ForEach(v,f,l)		for(v=f ; v<=l ; ++v)

#define	 StatBit(a,b)		(SpecialStatus[a]&(b))
#define	 Exclude(a)		StatBit(a,EXCLUDE)
#define	 Skip(a)		StatBit(a,EXCLUDE|SKIP)
#define  Discrete(a)		(MaxAttVal[a] || StatBit(a,DISCRETE))
#define  Continuous(a)		(! MaxAttVal[a] && ! StatBit(a,DISCRETE))
#define	 Ordered(a)		StatBit(a,ORDERED)
#define	 DateVal(a)		StatBit(a,DATEVAL)
#define	 TimeVal(a)		StatBit(a,STIMEVAL)
#define	 TStampVal(a)		StatBit(a,TSTMPVAL)

#define  Space(s)		(s==' ' || s=='\n' || s=='\r' || s=='\t')
#define  SkipComment		while ( ( c = InChar(f) ) != '\n' && c != EOF )

#define	 FreeUnlessNil(p)	if((p)!=Nil) free(p)
#define	 Free(x)	 	free(x)

#ifdef WIN32
#define  rint(x)		((x)<0 ? (double)((int)((x)-0.5)) :\
					 (double)((int)((x)+0.5)))
#define	 finite(x)		_finite(x)
#define	 strdup(x)		_strdup(x)
#endif

#define	 P1(x)			(rint((x)*10) / 10)
#define	 Of			stdout
#define	 Goodbye(x)		exit(x)
#define	 CharWidth(s)		((int) strlen(s))


#define	 NOFILE		 0
#define	 BADCLASSTHRESH	 1
#define	 LEQCLASSTHRESH	 2
#define	 BADATTNAME	 3
#define	 EOFINATT	 4
#define	 SINGLEATTVAL	 5
#define	 BADATTVAL	 6
#define	 BADNUMBER	 7
#define	 BADCLASS	 8
#define	 BADCOSTCLASS	 9
#define	 BADCOST	10
#define	 NOMEM		11
#define	 TOOMANYVALS	12
#define	 BADDISCRETE	13
#define	 NOTARGET	14
#define	 BADCTARGET	15
#define	 BADDTARGET	16
#define	 LONGNAME	17
#define	 HITEOF		18
#define	 MISSNAME	19
#define	 BADDATE	20
#define	 BADTIME	21
#define	 BADTSTMP	22
#define	 DUPATTNAME	23
#define	 UNKNOWNATT	24
#define	 BADDEF1	25
#define	 BADDEF2	26
#define	 BADDEF3	27
#define	 BADDEF4	28
#define	 SAMEATT	29
#define	 MODELFILE	30
#define	 CWTATTERR	31


/*************************************************************************/
/*									 */
/*		Type definitions					 */
/*									 */
/*************************************************************************/


typedef  unsigned char	Boolean, BranchType, *Set, Byte;
typedef	 char	const *String;

typedef  int	CaseNo;			/* data item number */
typedef  double	CaseCount;		/* count of (partial) items */

typedef  int	ClassNo,		/* class number, 1..MaxClass */
            DiscrValue,		/* discrete attribute value (0 = ?) */
            Attribute;		/* attribute number, 1..MaxAtt */

typedef	 float	ContValue;		/* continuous attribute value */
#define	 PREC	 7			/* precision */


typedef  union	 _def_val
{
    String	_s_val;		/* att val for comparison */
    ContValue	_n_val;		/* number for arith */
}
DefVal;


typedef  struct  _def_elt
{
    short	_op_code;	/* type of element */
    DefVal	_operand;	/* string or numeric value */
}
DefElt, *Definition;


typedef  struct  _elt_rec
{
    int		Fi,		/* index of first char of element */
            Li;		/* last ditto */
    char	Type;		/* 'B', 'S', or 'N' */
}
EltRec;


#define	 DefOp(DE)	DE._op_code
#define	 DefSVal(DE)	DE._operand._s_val
#define	 DefNVal(DE)	DE._operand._n_val

#define	 OP_ATT			 0	/* opcodes */
#define	 OP_NUM			 1
#define	 OP_STR			 2
#define	 OP_MISS		 3
#define	 OP_AND			10
#define	 OP_OR			11
#define	 OP_EQ			20
#define	 OP_NE			21
#define	 OP_GT			22
#define	 OP_GE			23
#define	 OP_LT			24
#define	 OP_LE			25
#define	 OP_SEQ			26
#define	 OP_SNE			27
#define	 OP_PLUS		30
#define	 OP_MINUS		31
#define	 OP_UMINUS		32
#define	 OP_MULT		33
#define	 OP_DIV			34
#define	 OP_MOD			35
#define	 OP_POW			36
#define	 OP_SIN			40
#define	 OP_COS			41
#define	 OP_TAN			42
#define	 OP_LOG			43
#define	 OP_EXP			44
#define	 OP_INT			45
#define	 OP_END			99


typedef  union  _attribute_value
{
    DiscrValue	_discr_val;
    ContValue	_cont_val;
}
AttValue, *DataRec;

#define  CVal(Case,Attribute)   Case[Attribute]._cont_val
#define  DVal(Case,Attribute)   Case[Attribute]._discr_val
#define  XDVal(Case,Att)	(Case[Att]._discr_val & 077777777)
#define  SVal(Case,Attribute)   Case[Attribute]._discr_val
#define  Class(Case)		(*Case)._discr_val
#define  Weight(Case)		(*(Case-1))._cont_val

#define	 Unknown(Case,Att)	(DVal(Case,Att)==UNKNOWN)
#define	 UnknownVal(AV)		(AV._discr_val==UNKNOWN)
#define	 NotApplic(Case,Att)	(DVal(Case,Att)==NA)
#define	 NotApplicVal(AV)	(AV._discr_val==NA)


typedef  struct _treerec	*Tree;
typedef  struct _treerec
{
    BranchType	NodeType;
    ClassNo	Leaf;		/* best class at this node */
    CaseCount	Cases,		/* no of items at this node */
                *ClassDist,	/* class distribution of items */
                Errors;		/* no of errors at this node */
    Attribute	Tested; 	/* attribute referenced in test */
    int		Forks;		/* number of branches at this node */
    ContValue	Cut,		/* threshold for continuous attribute */
                Lower,		/* lower limit of soft threshold */
                Upper,		/* upper limit ditto */
                Mid;		/* 50% point */
    Set         *Subset;	/* subsets of discrete values  */
    Tree	*Branch;	/* Branch[x] = subtree for outcome x */
}
TreeRec;


typedef  int	RuleNo;			/* rule number */

typedef  struct _condrec
{
    BranchType	NodeType;	/* test type (see tree nodes) */
    Attribute	Tested;		/* attribute tested */
    int		Forks;		/* possible branches */
    ContValue	Cut;		/* threshold (if relevant) */
    Set		Subset;		/* subset (if relevant) */
    int		TestValue,	/* specified outcome of test */
            TestI;		/* rule tree index of this test */
}
CondRec, *Condition;


typedef  struct _rulerec
{
    RuleNo	RNo;		/* rule number */
    int		TNo,		/* trial number */
            Size;		/* number of conditions */
    Condition	*Lhs;		/* conditions themselves */
    ClassNo	Rhs;		/* class given by rule */
    CaseCount	Cover,		/* number of cases covered by rule */
                Correct;	/* number on which correct */
    float	Prior;		/* prior probability of RHS */
    int		Vote;		/* unit = 0.001 */
}
RuleRec, *CRule;


typedef  struct _ruletreerec *RuleTree;
typedef  struct _ruletreerec
{
    RuleNo	*Fire;		/* rules matched at this node */
    Condition	CondTest;	/* new test */
    int		Forks;		/* number of branches */
    RuleTree	*Branch;	/* subtrees */
}
RuleTreeRec;


typedef struct _rulesetrec
{
    RuleNo	SNRules;	/* number of rules */
    CRule	*SRule;		/* rules */
    ClassNo	SDefault;	/* default class for this ruleset */
    RuleTree	RT;		/* rule tree for this ruleset */
}
RuleSetRec, *CRuleSet;


typedef	 struct _classify_environment
{
    CaseNo	Fp;		/* for SMP */
    double	*ClassWt;	/* total class votes */
    float	*Vote,		/* class boost votes */
            Confidence;	/* prediction CF */
    RuleNo	*Active,	/* active rules */
            ActiveSpace,	/* space for same */
            NActive;	/* number of same */
    CRule	*MostSpec;	/* most specific active rules */
    Boolean	*AttUsed;	/* reserved for possible later use */
    RuleNo	*RulesUsed,	/* all applicable rules */
            NRulesUsed;	/* number of same */
}
CEnvRec, *CEnv;


/*************************************************************************/
/*									 */
/*		Function prototypes					 */
/*									 */
/*************************************************************************/

/*************************************************************************/
/*									 */
/*		Text strings						 */
/*									 */
/*************************************************************************/


#define	 TX_Line(l,f)		"\n*** line %d of `%s': ", l, f
#define	 E_NOFILE(f,e)		"cannot open file %s%s\n", f, e
#define	 E_BADATTNAME		"`:' or `:=' expected after attribute name"\
					" `%s'\n"
#define	 E_EOFINATT		"unexpected eof while reading attribute `%s'\n"
#define	 E_SINGLEATTVAL(a,v)	"attribute `%s' has only one value `%s'\n",\
					a, v
#define	 E_DUPATTNAME		"multiple attributes with name `%s'\n"
#define	 E_CWTATTERR		"case weight attribute must be continuous\n"
#define	 E_BADATTVAL(v,a)	"bad value of `%s' for attribute `%s'\n", v, a
#define	 E_BADNUMBER(a)		"value of `%s' changed to `?'\n", a
#define	 E_BADCLASS		"bad class value `%s'l\n"
#define	 E_BADCLASSTHRESH	"bad class threshold `%s'\n"
#define	 E_LEQCLASSTHRESH	"class threshold `%s' <= previous threshold\n"
#define	 E_BADCOSTCLASS		"bad class `%s'\n"
#define	 E_BADCOST		"bad cost value `%s'\n"
#define	 E_NOMEM		"unable to allocate sufficient memory\n"
#define	 E_TOOMANYVALS(a,n)	"too many values for attribute `%s'"\
					" (max %d)\n", a, n
#define	 E_BADDISCRETE		"bad number of discrete values for attribute"\
					" `%s'\n"
#define	 E_NOTARGET		"target attribute `%s' not found\n"
#define	 E_BADCTARGET		"target attribute `%s' must be"\
					" type `continuous'\n"
#define	 E_BADDTARGET		"target attribute `%s' must be specified by"\
					" a list of discrete values\n"
#define	 E_LONGNAME		"overlength name: check data file formats\n"
#define	 E_HITEOF		"unexpected end of file\n"
#define	 E_MISSNAME		"missing name or value before `%s'\n"
#define	 E_BADTSTMP(d,a)	"bad timestamp `%s' for attribute `%s'\n", d, a
#define	 E_BADDATE(d,a)		"bad date `%s' for attribute `%s'\n", d, a
#define	 E_BADTIME(d,a)		"bad time `%s' for attribute `%s'\n", d, a
#define	 E_UNKNOWNATT		"unknown attribute name `%s'\n"
#define	 E_BADDEF1(a,s,x)	"in definition of attribute `%s':\n"\
					"\tat `%.12s': expect %s\n", a, s, x
#define	 E_BADDEF2(a,s,x)	"in definition of attribute `%s':\n"\
					"\t`%s': %s\n", a, s, x
#define	 E_BADDEF3		"cannot define target attribute `%s'\n"
#define	 E_BADDEF4		"[warning] target attribute appears in"\
					" definition of attribute `%s'\n"
#define	 E_SAMEATT(a,b)		"[warning] attribute `%s' is identical to"\
					" attribute `%s'\n", a, b
#define	 EX_MODELFILE(f)	"file %s incompatible with .names file\n", f
#define	 E_MFATT		"undefined or excluded attribute"
#define	 E_MFATTVAL		"undefined attribute value"
#define	 E_MFCLASS		"undefined class"
#define	 E_MFEOF		"unexpected eof"
#define	 T_ErrorLimit		"Error limit exceeded\n"

#define	MAXLINEBUFFER	10000

#define FailSyn(Msg)	 {DefSyntaxError(Msg); return false;}
#define FailSem(Msg)	 {DefSemanticsError(Fi, Msg, OpCode); return false;}

typedef  union  _xstack_elt
         {
            DiscrValue  _discr_val;
            ContValue   _cont_val;
            String      _string_val;
         }
	 XStackElt;

#define	cval		_cont_val
#define	sval		_string_val
#define	dval		_discr_val
#define XError(a,b,c)	Error(a,b,c)
#define	CUnknownVal(AV)		(AV.cval==_UNK.cval)
#define	DUnknownVal(AV)		(AV.dval==_UNK.dval)
#define DUNA(a)	(DUnknownVal(XStack[a]) || NotApplicVal(XStack[a]))
#define CUNA(a)	(CUnknownVal(XStack[a]) || NotApplicVal(XStack[a]))
#define	C1(x)	(CUNA(XSN-1) ? _UNK.cval : (x))
#define	C2(x)	(CUNA(XSN-1) || CUNA(XSN-2) ? _UNK.cval : (x))
#define	CD2(x)	(CUNA(XSN-1) || CUNA(XSN-2) ? _UNK.dval : (x))
#define	D2(x)	(DUNA(XSN-1) || DUNA(XSN-2) ? _UNK.dval : (x))



#define	PROPS 23

#define	ERRORP		0
#define ATTP		1
#define CLASSP		2
#define CUTP		3
#define	CONDSP		4
#define ELTSP		5
#define ENTRIESP	6
#define FORKSP		7
#define FREQP		8
#define IDP		9
#define TYPEP		10
#define LOWP		11
#define MIDP		12
#define HIGHP		13
#define RESULTP		14
#define RULESP		15
#define VALP		16
#define LIFTP		17
#define COVERP		18
#define OKP		19
#define DEFAULTP	20
#define COSTSP		21

class C50Classifier  {

public :
    String	OptArg, Option;
int	Entry;
Condition	*Test=Nil;
int		NTest,
        TestSpace,
        *TestOccur=Nil,
         *RuleCondOK=Nil;
Boolean		*TestUsed=Nil;
char	PropName[20],
        *PropVal=Nil,
         *Unquoted;
int	PropValSize=0;


    Boolean	    ReadName(FILE *f, String s, int n, char ColonOpt);
    void	    GetNames(FILE *Nf);
    void	    ExplicitAtt(FILE *Nf);
    int	    Which(String Val, String *List, int First, int Last);
    int	    InChar(FILE *f);

    DataRec GetDataRecFromVec(std::vector<std::string> &datavec, Boolean Train);
    DataRec	GetDataRec(FILE *Df, Boolean Train);
    
    int	    StoreIVal(String S);
    void	    CheckValue(DataRec DVec, Attribute Att);

    void	    ImplicitAtt(FILE *Nf);
    void	    ReadDefinition(FILE *f);
    void	    Append(char c);
    Boolean	    Expression();
    Boolean	    Conjunct();
    Boolean	    SExpression();
    Boolean	    AExpression();
    Boolean	    Term();
    Boolean	    Factor();
    Boolean	    Primary();
    Boolean	    Atom();
    Boolean	    Find(String S);
    int	    FindOne(String *Alt);
    Attribute   FindAttName();
    void	    DefSyntaxError(String Msg);
    void	    DefSemanticsError(int Fi, String Msg, int OpCode);
    void	    Dump(char OpCode, ContValue F, String S, int Fi);
    void	    DumpOp(char OpCode, int Fi);
    Boolean	    UpdateTStack(char OpCode, ContValue F, String S, int Fi);
    AttValue    EvaluateDef(Definition D, DataRec Case);

    void	    ReadFilePrefix(String Extension);
    void	    ReadHeader();
    Tree	    GetTree(String Extension);
    Tree	    InTree();
    CRuleSet    GetRules(String Extension);
    CRuleSet    InRules();
    CRule	    InRule();
    Condition   InCondition();
    void	    ConstructRuleTree(CRuleSet RS);
    void	    SetTestIndex(Condition C);
    RuleTree    GrowRT(RuleNo *RR, int RRN, CRule *Rule);
    int	    DesiredOutcome(CRule R, int TI);
    int	    SelectTest(RuleNo *RR, int RRN, CRule *Rule);
    int	    ReadProp(char *Delim);
    char*    RemoveQuotes(String S);
    Set	    MakeSubset(Attribute Att);
    Tree	    Leaf(double *Freq, ClassNo NodeClass, CaseCount Cases);

    void	    GetMCosts(FILE *f);

    ClassNo	    TreeClassify(DataRec CaseDesc, Tree DecisionTree, CEnv E);
    void	    FollowAllBranches(DataRec CaseDesc, Tree T, float Fraction,
                                  double *Prob, Boolean *AttUsed);
    void	    FindLeaf(DataRec CaseDesc, Tree T, Tree PT, float Wt, double *Prob,
                         Boolean *AttUsed);
    ClassNo	    RuleClassify(DataRec CaseDesc, CRuleSet RS, CEnv E);
    int	    FindOutcome(DataRec Case, Condition OneCond);
    Boolean	    Satisfies(DataRec CaseDesc, Condition OneCond);
    Boolean	    Matches(CRule R, DataRec Case);
    void	    CheckActiveSpace(int N, CEnv E);
    void	    MarkActive(RuleTree RT, DataRec Case, CEnv E);
    ClassNo	    BoostClassify(DataRec CaseDesc, int MaxTrial, CEnv E);
    ClassNo	    SelectClass(ClassNo Default, Boolean UseCosts, double *Prob);
    double	    MisclassCost(double *LocalFreq, ClassNo C);
    ClassNo	    Classify(DataRec CaseDesc, CEnv E);
    float	    Interpolate(Tree T, ContValue Val);

    FILE *	    GetFile(String Extension, String RW);
    void	    CheckFile(String Extension, Boolean Write);

    char	    ProcessOption(int Argc, char *Argv[], char *Options);
    void	    *Pmalloc(size_t Bytes);
    void	    *Prealloc(void *Present, size_t Bytes);
    void	    *Pcalloc(size_t Number, unsigned Size);
    void	    Error(int ErrNo, String S1, String S2);
    int	    Denominator(ContValue Val);
    int	    GetInt(String S, int N);
    int	    DateToDay(String DS);
    int	    TimeToSecs(String TS);
    void	    SetTSBase(int y);
    int	    TStampToMins(String TS);

    void	    FreeGlobals();
    void	    FreeCosts();
    void	    FreeNames();
    void	    FreeTree(Tree T);
    void	    FreeRule(CRule R);
    void	    FreeRuleTree(RuleTree RT);
    void	    FreeRules(CRuleSet RS);
    void	    FreeLastCase(DataRec DVec);
    void	    FreeVector(void **V, int First, int Last);
    
    void      loadModel();
    
    ~C50Classifier();


    int		TRIALS=1,	/* number of trees to be grown */
            Trial;		/* trial number for boosting */

    Boolean		RULES=0,	/* rule-based classifiers */
                RULESUSED=0;	/* list applicable rules */


    /*************************************************************************/
    /*									 */
    /*		Attributes and data					 */
    /*									 */
    /*************************************************************************/

    Attribute	ClassAtt=0,	/* attribute to use as class */
                LabelAtt,	/* attribute to use as case ID */
                CWtAtt;		/* attribute to use for case weight */

    String		*ClassName=0,	/* class names */
                 *AttName=0,	/* att names */
                  **AttValName=0;	/* att value names */

    char		*IgnoredVals=0;	/* values of labels and atts marked ignore */
    int		IValsSize=0,	/* size of above */
            IValsOffset=0;	/* index of first free char */

    int		MaxAtt,		/* max att number */
            MaxClass=0,	/* max class number */
            AttExIn=0,	/* attribute exclusions/inclusions */
            LineNo=0,	/* input line number */
            ErrMsgs=0,	/* errors found */
            Delimiter,	/* character at end of name */
            TSBase=0;	/* base day for time stamps */

    DiscrValue	*MaxAttVal=0;	/* number of values for each att */

    ContValue	*ClassThresh=0;	/* thresholded class attribute */

    char		*SpecialStatus=0;/* special att treatment */

    Definition	*AttDef=0;	/* definitions of implicit atts */

    Boolean		*SomeMiss=Nil,	/* att has missing values */
                 *SomeNA=Nil;	/* att has N/A values */

    String		FileStem="undefined";

    /*************************************************************************/
    /*									 */
    /*		Classification environment				 */
    /*									 */
    /*************************************************************************/

    CEnvRec		*GCEnv;		/* used by classification routines */

    /*************************************************************************/
    /*									 */
    /*		Trees							 */
    /*									 */
    /*************************************************************************/

    Tree		*Pruned=0;	/* decision trees */

    ClassNo		*TrialPred=0;	/* predictions for each boost trial */

    float		**MCost=0;	/* misclass cost [pred][real] */

    /*************************************************************************/
    /*									 */
    /*		Rules							 */
    /*									 */
    /*************************************************************************/


    CRuleSet	*RuleSet=0;	/* rulesets */

    ClassNo		Default;	/* default class associated with ruleset or
				   boosted classifier */


    /*************************************************************************/
    /*									 */
    /*		Misc							 */
    /*									 */
    /*************************************************************************/

    FILE		*TRf=0;		/* file pointer for tree and rule i/o */
    char		Fn[500];	/* file name */
    char	LineBuffer[MAXLINEBUFFER], *LBp=LineBuffer;

char	*Buff;			/* buffer for input characters */
int	BuffSize, BN;		/* size and index of next character */

EltRec	*TStack;		/* expression stack model */
int	TStackSize, TSN;	/* size of stack and index of next entry */

int	DefSize, DN;		/* size of definition and next element */

Boolean PreviousError;		/* to avoid parasytic errors */

AttValue _UNK,			/* quasi-constant for unknown value */
	 _NA;			/* ditto for not applicable */


};
