/* Gap Weighted Subsequence Kernel 
 * 
 * This implementation is based on the gapped substring kernel algorithm by 
 * Juho Rousu and John Shawe-Taylor.
 * 
 * Citation:
 * Juho Rousu, John Shawe-Taylor. Efficient Computation of Gapped Substring 
 * Kernels on Large Alphabets. Journal of Machine Learning Research 6 (2005), 
 * 1323-1344
 * 
 */

#include <algorithm>
#include <memory>
#include <vector>
#include <cmath>
//#include <iostream>
#include <assert.h>
//#include <string>

using namespace std;
/*typedef unsigned int SEQelem;
typedef struct {
  int len;
  SEQelem *S;
} SEQ;
*/
typedef unsigned int uint;
typedef uint SEQ;

/*
class MatchListElem {
  public:
  uint len;
  uint* p; // position
  double* w; //weight
  MatchListElem() {
	p=NULL; w=NULL; len=0;
  }
  MatchListElem(uint* position, double* weight, uint length):
    p(position), w(weight), len(length) 
  {}
  MatchListElem(int size) : len(size), p(new uint[size]), w(new double[size]) {}
  MatchListElem(const MatchListElem& mle) : 
    len(mle.len), p(new uint[mle.len]), w(new double[mle.len]) {
      copy(mle.p, mle.p + mle.len, p);
      copy(mle.w, mle.w + mle.len, w);
    }
  void setElem(uint* position, double* weight, uint length) {
    assert(len==0);
    len = length; p=position; w=weight;
  }
  ~MatchListElem() {delete[] p; delete[] w;}
};
*/

typedef struct { int p; double w; } PosWeight;
typedef vector<PosWeight> MatchListElem;
typedef vector<MatchListElem> MatchList;
typedef vector< vector<int> > multivec;

MatchList* createMatchListFrom(const SEQ* s, const SEQ* t, int m, int n, int alphaSize, double lambda, int HSH, int VSW) {
    
    multivec index;
    index.resize(alphaSize);
    for (int i=0;i<n;i++) {
        index.at((int) t[i]).push_back(i);
    }
    
    MatchList* L = new MatchList();
    
    MatchListElem Buff;
    for (int i=0;i<m;i++) {
        int h_strip_max = ceil(((double) i+1)/HSH)*HSH;
        vector<int>& ind = index.at((int) s[i]);
        
        Buff.resize(ind.size());
        for (int j=0;j<ind.size();j++) {
            int v_strip_max = ceil(((double) ind[j]+1)/VSW)*VSW;
            Buff[j].p = ind[j];
            Buff[j].w = pow(lambda, 2+h_strip_max-(i+1)+v_strip_max-(ind[j]+1));
        }
        L->push_back(Buff);
    }
    
    return L;
}


void precompute_rangepaths_opt(int n, multivec& Query, multivec& Update) {
    int highbit = ceil(log(n+1.)/log(2.));
    //int* powersof2 = new int[highbit];
    //for (int i=0; i<highbit; i++) powersof2[i] = power(2, highbit-1-i);
    // Note: element 0 is left empty for efficiency: The tree structrure depends on this fact.
    Query.resize(n);
    Update.resize(n);
    int rootindex = pow(2.0,highbit);
    
    for (int j=1; j<n; j++) {
        int* path = new int[highbit];
        //int bmask = pow(2,highbit+1)-1;
        path[0] = rootindex>>1;
        for (int b=1; b<highbit; b++) {
            int xx = 1<<(highbit-b);
            path[b] = path[b-1] + (xx>>1) * ( xx&j ? 1 : -1);
        }
        
        /* confusing code...
                highbit = most significant bit
                path high to low
                bits (in matlab) = high to low
                querypath = path[bits==1]
        
                    */
        int lsb;
        for (lsb=0;lsb<highbit;lsb++) {
            if (j & (1<<lsb)) break;
        }
        lsb = highbit - lsb - 1;
        vector<int>& Q = Query.at(j);
        vector<int>& U = Update.at(j-1);
        for (int b=0; b<lsb; b++) {
            int xx = 1<<(highbit-b-1);
            if (xx&j) 
                Q.push_back(path[b]);
            else 
                U.push_back(path[b]);
        }
        Q.push_back(path[lsb]);
        U.push_back(path[lsb]);
        delete[] path;
    }
}
    

//extern "C" {

void GapWeightedSubsequences_Sparse(double* kvalues, const SEQ* s, const SEQ* t, int m, int n, int seqlen, int alphaSize, double lambda, int minseqlen = -1) {
    
    
    // Striping definitions
    int VSW; // V_STRIP_WIDTH;
    int HSH; // H_STRIP_HEIGHT;
    
    if (lambda < 1 && 1) {
        VSW = (int) pow(2, floor(log(min(-150/log10(lambda),(double) n))/log(2.)));
        HSH = (int) pow(2, floor(log(min(-150/log10(lambda),(double) m))/log(2.)));
    } else {
        VSW = n;
        HSH = m;
    }
    
    MatchList* MatchListInit = createMatchListFrom(s, t, m, n, alphaSize, lambda, HSH, VSW);
    
    // should we compute a range of kernels?
    if (minseqlen==-1) minseqlen = seqlen;
    assert(minseqlen<=seqlen && minseqlen>=1);
    //double[] kvalues = new double[seqlen-minseqlen+1];
    
    if (minseqlen==1) {
        kvalues[0] = 0;
        // compute value for length 1 = number of matches in MatchListInit
        for (MatchList::iterator i = MatchListInit->begin(); i!=MatchListInit->end();i++) {
            kvalues[0] += i->size();
        }
    }

        
    
    
    vector<double> Weight;
	Weight.resize((int) ceil((double) n/VSW)*VSW);
    int Wsize = ceil((double)n/VSW)*VSW;
    int Wstrip = VSW; //ceil((double)n/VSW);
    
    multivec QueryPath; 
    multivec UpdatePath;
    precompute_rangepaths_opt(min(n,VSW+1), QueryPath, UpdatePath); // or VSW+1?
    
    // create a buffer for new matchlists
    MatchListElem NewL;
    
    for(int l=2; l <= seqlen; l++) {
        // more strips
        int h_strip_number = 0;
        int h_strip_max = HSH;
        // MatchList& OldMatchList = *MatchListInit;
        MatchList& NewMatchList = *(new MatchList(m));
        
		int firstrow=1;
        
        // (re-)set weights
        for (int ii=0;ii<Wsize;ii++) {
            Weight[ii] = 0;
        }
        
        //process match lists
        for (int i=0; i<m; i++) {
            MatchListElem& L = MatchListInit->at(i);
            
            if(L.size()>0) {
                NewL.resize(0); // resetting length to zero
                if (firstrow) {
					firstrow = 0; // ignore the first row.
                    // NewMatchList.push_back(MatchListElem());
                } 
                else { 
                    if (i>=h_strip_max) {
                        int strips_passed = ceil((double)(i-h_strip_max+1)/HSH);
                        h_strip_number = h_strip_number + strips_passed;
                        h_strip_max = h_strip_max + strips_passed*HSH;
                        
                        // rescale Weights
                        for (int ii=0;ii<Wsize;ii++) {
                          Weight[ii]=Weight[ii]*pow(lambda,strips_passed*HSH);
                        }
                    }
                    
                    // more striping
                    int v_strip_number = 0;
                    int v_strip_max = VSW;
                    double v_strip_sum = 0;
                    
                    // TODO check for 0/1 array start probs
                    for (int j = 0; j<L.size(); j++) {
                        int node = L[j].p;
                        while(node >= v_strip_max) {
                            v_strip_sum = pow(lambda, VSW)*
                                (v_strip_sum+Weight[v_strip_number*Wstrip+VSW-1]);
                            v_strip_number++;
                            v_strip_max += VSW;
                        }
                        
                        double sumOfWeights = v_strip_sum;
                        int qnode = node - v_strip_max + VSW;
                        
                        vector<int>& Q = QueryPath.at(qnode);
                     
                        for (int i = 0;i<Q.size();i++) {
                            int iq = Q[i]-1;
                            sumOfWeights += Weight[v_strip_number*Wstrip+iq];
                        }
                        
                        if (sumOfWeights>0) {

							PosWeight x = {node, sumOfWeights};
                            NewL.push_back(x);
                        }
                    }
                    if (NewL.size()>0) {
                        NewMatchList.at(i) = NewL; // will copy current size of NewL
                    } else {
                        //NewMatchList.push_back(MatchListElem());
                    }
                }    
                
                int v_strip_number = 0;
                int v_strip_max = VSW;
                
                for (int j=0; j < L.size(); j++) {
                    uint node = L[j].p;
                    while(node >= v_strip_max) {
                        v_strip_number++;
                        v_strip_max += VSW;
                    }
                    
                    int w_node = node - v_strip_max + VSW;
                    
                    vector<int>& U = UpdatePath.at(w_node);
                     
                    for (int i = 0;i<U.size();i++) {
                        int iw = U[i]-1;
                        Weight.at(v_strip_number*Wstrip+iw) += L[j].w;
                    } 
                }    
            }
        }
        
        delete MatchListInit;
        MatchListInit = &NewMatchList;
        
        if (minseqlen<=l) {
            //compute kernel value
            double Kappa = 0;
            
            for (int i=0;i<m;i++) {
                int h_strip_max = ceil(((double)i+1)/HSH)*HSH;
                MatchListElem& L = MatchListInit->at(i);
                
                for (int j=0;j<L.size();j++) {
                    int v_strip_max = ceil(((double)L[j].p+1)/VSW)*VSW;
        			double expo = (i+1)-h_strip_max + (L[j].p+1)-v_strip_max;
        			double lpow = pow(lambda, expo);
                    Kappa += L[j].w*lpow;
                }
            }
            
            kvalues[l-minseqlen] = Kappa/pow(lambda, 2*l);
        }
            
        
    }    
    
    // Compute Kappa
 
    delete MatchListInit;

    
}

//}
/*
int main(void) {

    cout << "Test";
	// the cat was chased by the fat dog
	double lambda = 0.8; 
    SEQ a[] = {0,1,2,3,4,0,5,6}; int m = 8;
	// the fat cat bit the dog
    SEQ b[] = {0,5,1,7,0,6}; int n = 6;
    double k  = GapWeightedSubsequences_Sparse(&a[0], &b[0], m, n, 2, 8, lambda);
    double k2 = GapWeightedSubsequences_Sparse(&b[0], &a[0], n, m, 2, 8, lambda);
    cout << "GapWeighted " << k << " "; // << k2;
    return 0;

}

*/
extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

  SEXP gapkernel_range(SEXP rtext, // text document
		SEXP ltext, // list or vector of text documents to compute kvalues against
		SEXP rntok,  // number of tokens in text
		SEXP lntok, // tokens per document in ltext
		SEXP alphasize,   // type of kernel
		SEXP len, // len of sequence
        SEXP minlen, // minimum length for range
		SEXP lamb)   // parameter for kernel
  {
    // R interface for text and list of text computation. Should return a vector of computed kernel values.
    // Construct ESASK
    int r_size = *INTEGER(rntok);
    int l_size = *INTEGER(lntok);
    const SEQ *rtxt = (SEQ*) INTEGER(rtext);
	const SEQ *ltxt = (SEQ*) INTEGER(ltext);
    double lambda = *REAL(lamb); 
    int as = *INTEGER(alphasize);
	int seqlen = *INTEGER(len);
    int minseqlen = *INTEGER(minlen);
    SEXP kernel;
    
    PROTECT(kernel = allocVector(REALSXP, seqlen-minseqlen+1));
    
	GapWeightedSubsequences_Sparse(REAL(kernel), rtxt, ltxt, r_size, l_size, seqlen, as, lambda, minseqlen);

	//symmetry check

    UNPROTECT(1); 
    return kernel;

  }
}
