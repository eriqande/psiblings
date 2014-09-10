// familyFinder.cc

// Partitions a single generation of individuals into full-sib families 
// based on genetic data.

// input: genotype file =  tab-delimited text file
//        comment lines, which begin with #, are ignored
//        first non-comment line is number of individuals (individuals=60)
//        second non-comment line is number of loci (loci=11)
//        third non-comment line is column headings
//        rest of lines are data:
//          column 1 is the individual id
//          each add'l column is the genotype for a given locus, in the format
//          ... allele1-allele2

// output: a file containing a list of the full-sib families followed by the
//         original data set organized by family.

#include <R.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "list.h"
#include "allele.h"
#include "intList.h"
#include "intNode.h"
#include "node.h"
#include "edge.h"
#include "graph.h"

#define MINSIZE 5

int parse_int( string str );
string parse_allele( string geno, int al_num );
int update_locus( Allele * ap, List l_arr[], int loc, string al_name );
void update_genotype( int * g_arr[], int indiv, int loc, int al_num, int cv );
void check_maxAlleles( int & max, List l_arr[], int loc );
void update_freqs( List l_arr[], double * f_arr[], int nLoc );
void set_permutation( int p, int & Xm, int & Xp, int & Ym, int & Yp,
                      int ind1, int ind2, int loc, int * g_arr[] );
char get_equation( int Xm, int Xp, int Ym, int Yp );
double mSame_pSame( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp );
double mSame_pDiff( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp, const double & Pyp );
double mDiff_pSame( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp, const double & Pym );
double mDiff_pDiff( const double & Rm, const double & Rp, const double & Pxm,
                   const double & Pxp, const double & Pym, const double & Pyp);
int draw_from_pop( double * f_arr[], int loc );
void merge_sort( double A[], int p, int r );
void merge( double A[], int p, int q, int r );
double get_t2_error( double ratios_arr[], double x_val ); 

int connected_components( intList adjVert[], int n, int cc[] );
void ccDFS( intList adjVert[], int color[], int current, int component, int cc[] );
int countSibs( int famNum, int nIndiv, int cc[] );
int countEdges( int famNum, int nIndiv, int cc[], intList adjVert[] );
double getScore( int famNum, int nIndiv, int cc[], intList adjVert[] );
bool scores_higher( double s1best, double s2best, double s1, double s2 );
void number_for_ghc( int nIndiv, int famNum, int cc[], int ghc_num[] );
void convert_graph( Graph * gr, intList adjVert[], int famNum, int n,
                    int m, int nIndiv, int cc[], int ghc_num[], string id_arr[] );
bool ghc_tree( Graph * gr, Node * active[], int number[],
               int & max_dist, int & bound, bool & co_check );
int maxflow( Graph * gr, Node * sptr, Node * tptr, Node * active[],
             int number[], int & max_dist, int & bound, bool & co_check );
void initialize_maxflow( Graph * gr, Node * active[], int number[] );
void global_relabel( Graph * gr, Node * tptr, Node * active[],
                     int number[], int & max_dist, int & bound );
void compute_subtreeSizes( Graph * gr );
void percolate_increment( Node * nptr, Graph * gr );
Node * find_mincut( Graph * gr );
bool size_ok( Node * cut_ptr, Graph * gr );
void label_T_nodes( Node * cut_ptr, Graph * gr, int T_size );
void update_components( int nIndiv, int famNum, int nFams,
                        int cc[], int ghc_nodes[], Graph * gr );
void mark_ST_edges( int nIndiv, int cc[], int T_id, int S_id, intList adjVert[] );
void restore_components( int nIndiv, int cc[], int nFams, int old_famNum  );
void print_tree( Graph * gr );




/* now L is the length of e1 and e2.  They are the 
vectors giving indices of the pairwise adjacencies */
extern "C" {  /* I wrap this in this extern "C" block so that R can find all right in the shared object file */
	 void famfind(int *NumIndivs, int *e1, int *e2, int *L,   int *FamVec,  int *FamStart, int *FamEnd, int *NumFams ) {
		 ifstream infile;
		 ofstream outfile;
		 int nIndiv, nLoci, relation, edges;
		 int max_dist, bound, T_id, T_size, famNum, locusNum, conv;
		 double fullscore, best_score1, best_score2, score1, score2;
		 string infile_name, outfile_name, str, allele1, allele2;
		 string node1, node2, node1fam, node2fam, allele_name;
		 bool tree_ok, done, co_check;
		 Allele * ap;
		 IntNode * iPtr1, * iPtr2, * ip;
		 Node * cut_ptr, * best_ptr;



	 //-------------------  done calculating pairwise likelihoods  ----------------------------

		 // Eric put this here to do a simple test example:
		 nIndiv = *NumIndivs;
		 nLoci = 1;

		 // these are here for now just to stop the compiler errors after deleting all the initial stuff	
		 // mostly these are things that are stored so that the original data file can be written out in
		 // the end, or something like that
		 int ** results = new int * [nIndiv];                     // declare results matrix 
		 for( int j = 0; j < nIndiv; j++ ) {                      // stores pairwise relationships
			 results[j] = new int[nIndiv];                          // represented as ints
		 }
	
		 for( int i = 0; i < nIndiv; i++ ) {                      // initialize results matrix
			 for( int j = 0; j < nIndiv; j++ ) {
				 results[i][j] = 1;
			 }
		 }
	
		 string * id_arr = new string[nIndiv];
		 string * headings = new string[nLoci + 1];
		 int ** genotypes = new int * [nIndiv];                  // declare genotypes matrix
		 List * locus_arr = new List[nLoci + 1];                 // declare allele lists
	
	
		 // now we make the adjacency matrix
		 for( int i=0; i<(*L); i++) {
			 results[ e1[i]-1 ][ e2[i]-1 ] = 11;
			 results[ e2[i]-1 ][ e1[i]-1 ] = 11;
		 }
	
	 /*  // right here, Eric is going to give values to some variables to test some 
		 // stuff.  Let's imagine a cluster of size 7 (or 10, etc) with the first three fully connected and the second
		 // four fully connected, and a single edge connected the first three to the second four
		 for( int indiv1 = 1; indiv1 < nIndiv; indiv1++ ) {
			 for( int indiv2 = 0; indiv2 < indiv1; indiv2++ ) {
				 if(indiv1<5 && indiv2<5) {
					 results[indiv1][indiv2] = 11;
				 } 
				 else if(indiv1>=5 && indiv2>=5) {
					 results[indiv1][indiv2] = 11;
				 }
			 }
		 }
		 results[6][3] = 11;  // this is the one edge that connects the two cliques
		 results[6][2] = 11;
		 results[6][1] = 11;
		 results[6][0] = 11;
	 */
	
		 for(int i=0;i<nIndiv; i++) {
			 char buff[100];
			 sprintf(buff, "%d", i+1);
			 id_arr[i] = buff;
		 }
	 /*  std::string buffAsStdStr = buff;
	
		 id_arr[5] = buff; 
	
		 id_arr[0] = "dude1";
		 id_arr[1] = "dude2";
		 id_arr[2] = "dude3";
		 id_arr[3] = "dude4";
		 id_arr[4] = "dude5";
		 id_arr[5] = "dude6";
		 id_arr[6] = "dude7";
		 id_arr[7] = "dude8";
		 id_arr[8] = "dude9";
		 id_arr[9] = "dude10";
	 //  id_arr[6] = "dude7";
	 */

	
	
		 outfile.open( "boing.txt" );
		 if( !outfile ) {
			 cout << "Unable to open output file: " << outfile_name << endl;
			 cout << "Program aborted" << endl  << endl;
			 exit(-3);
		 }
	
	
	
	

		 intList * adj_list = new intList[nIndiv];           // declare adjacency list
		 int * components = new int[nIndiv];

		 // build adjacency list
		 for( int indiv1 = 1; indiv1 < nIndiv; indiv1++ ) {
			 for( int indiv2 = 0; indiv2 < indiv1; indiv2++ ) {
				 relation = results[indiv1][indiv2];

				 if( (relation == 11) || (relation == 12) || (relation == 13) ) {
					 iPtr1 = new IntNode( indiv1 );
					 iPtr2 = new IntNode( indiv2 );
					 adj_list[indiv1].insert( iPtr2 );
					 adj_list[indiv2].insert( iPtr1 );
				 }
			 }
		 }
	
		 for( int i = 0; i < nIndiv; i++ )
			 delete [] results[i];
		 delete [] results;

		 // find connected components
		 int nFams = connected_components( adj_list, nIndiv, components );

		 int * ghc_nodes = new int[nIndiv];           // holds within-family id's (1..n)
		 intList famQueue;

		 // enqueue all families found by connected_components()
		 for( int i = 1; i <= nFams; i++ ) {
			 ip = new IntNode( i );
			 famQueue.insert( ip );
		 }

		 while( !famQueue.isEmpty() ) {
			 ip = famQueue.get_first();
			 famNum = ip->get_id_num();
			 famQueue.delete_first();

			 int famSize = countSibs( famNum, nIndiv, components );
			 edges = countEdges( famNum, nIndiv, components, adj_list );
			 number_for_ghc( nIndiv, famNum, components, ghc_nodes );
			 fullscore = getScore( famNum, nIndiv, components, adj_list ); // score the original family

			 if( famSize >= (2 * MINSIZE) ) {
				 // convert family to Graph data structure
				 Graph * graph_ptr = new Graph( famSize, edges );
				 convert_graph( graph_ptr, adj_list, famNum, famSize, edges,
												nIndiv, components, ghc_nodes, id_arr );

				 // data structures for maxflow()
				 Node ** active = new Node *[famSize+1];
				 int * number = new int[famSize+1];

				 // build a Gomory-Hu cut tree
				 tree_ok = ghc_tree( graph_ptr, active, number, max_dist, bound, co_check );

				 delete [] active;
				 delete [] number;

				 if( tree_ok ) {
					 done = false;
					 best_ptr = 0;
					 best_score1 = 0;
					 best_score2 = 0;

					 while( 1 ) {
						 if( !done )
							 cut_ptr = find_mincut( graph_ptr );

						 if( cut_ptr ) {
							 // label nodes that are in T
							 T_size = cut_ptr->get_subtreeSize();
							 label_T_nodes( cut_ptr, graph_ptr, T_size );

							 // modify cc[], T nodes get id nFams+1
							 T_id = nFams+1;
							 update_components( nIndiv, famNum, nFams, components, ghc_nodes, graph_ptr );

							 // modify adjacency list, mark edges connecting S and T
							 mark_ST_edges( nIndiv, components, T_id, famNum, adj_list );

							 // score the two new families
							 score1 = getScore( famNum, nIndiv, components, adj_list );
							 score2 = getScore( T_id, nIndiv, components, adj_list );

							 if( (score1 > fullscore) && (score2 > fullscore) ) {
								 if( done ) {
									 // remove marked edges (connecting S & T) from adj_list
									 for( int i = 0; i < nIndiv; i++ )
										 adj_list[i].delete_notInList();
									 nFams++;
									 ip = new IntNode( famNum );
									 famQueue.insert( ip );
									 ip = new IntNode( nFams );
									 famQueue.insert( ip );
									 break;
								 }
								 else if( scores_higher(best_score1, best_score2, score1, score2) ) {
									 best_score1 = score1;
									 best_score2 = score2;
									 best_ptr = cut_ptr;
								 }
							 }

							 // change all nodes with T_id back to component i
							 // and unmark all previously marked edges in adj_list
							 // then check for a cut with a better score

							 restore_components( nIndiv, components, nFams, famNum );
							 for( int i = 0; i < nIndiv; i++ )
								 adj_list[i].restore_notInList();

						 }//if cut_ptr

						 else {       // cut_ptr == null
							 if( best_ptr != 0 ) {
								 cut_ptr = best_ptr;
								 done = true;
							 }
							 else break;
						 }

					 }//while(1)
				 }//if tree_ok

				 delete graph_ptr;

			 }//if famSize >= 2MINSIZE
		 }//while(famQueue)


	 // -----------------------  write to outfile  -------------------------------------------

		 outfile << "# file: " << outfile_name << endl;
		 outfile << "# input file: " << infile_name << endl << endl;
		 outfile << "individuals = " << nIndiv << endl;
		 outfile << "loci = " << nLoci << endl;
		 outfile << "families = " << nFams << endl << endl;

		 // write original data file organized by family
		 outfile << "\n------------------------------" << endl;
		 outfile << "Input file arranged by family" << endl;
		 outfile << "------------------------------\n\n" << endl;

		 outfile << endl << endl;


	 // put this in a block so I can declare some variables
	 *NumFams=nFams;  // send this back to R
	 { int s=0; 	
	
		 for( int i = 1; i <= nFams; i++ ) {
			 FamStart[i-1] = s+1;  // the new start point is the last endpoint + 1
			 outfile << "family " << i << ":  " << endl;
			 for( int j = 0; j < nIndiv; j++ ) {
				 if( components[j] == i ) {
					 outfile << id_arr[j];
					 outfile << endl;
					 FamVec[s++]=atoi(id_arr[j].c_str());
				 }
			 }
			 FamEnd[i-1]=s;
			 outfile << endl;
		 }
	 }






		 // free dynamically allocated memory

	 /*  for( int i = 0; i < nIndiv; i++ )
			 delete [] genotypes[i];
		 delete [] genotypes;
 
		 delete [] id_arr;
		 delete [] headings;
		 delete [] locus_arr;
		 delete [] adj_list;
		 delete [] components;
		 delete [] ghc_nodes;
	 */

		 outfile.close();

		 //cout << "\n\nResults are in file: " << outfile_name << endl << endl;

		 //return 0;
	}
}




//--------------------------  function definitions  -------------------------------------



/*---------------------------  parse_int()  --------------------------------
|   Parses an integer from an expression of form: identifier=integer.      |
|   The integer is returned.                                               |
--------------------------------------------------------------------------*/

int parse_int( string str ) {
  int number;
  string s;

  int x = str.find("=");
  int lastChar = str.length() - 1;

  s = str.substr( x+1, lastChar - x );
  number = atoi( s.c_str() );

  return number;
}



/*---------------------------  parse_allele()  -----------------------------
|   Parses the specified allele from a genotype string.                    |
|   The allele string is returned.                                         |
--------------------------------------------------------------------------*/

string parse_allele( string geno, int al_num ) {
  string str;
  int x = geno.find("-");
  int lastChar = geno.length() - 1;

  if( al_num == 1 )
    str = geno.substr( 0, x );
  else
    str = geno.substr( x+1, lastChar - x );

  return str;
}



/*---------------------------  update_locus()  -----------------------------
|   Updates the locus matrix with the new allele.                          |
|   If the allele already exists, its count is incremented,                |
|   otherwise a new allele node is created.                                |
|   The integer representation of the allele is returned.                  |
--------------------------------------------------------------------------*/

int update_locus( Allele * ap, List l_arr[], int loc, string al_name ) {
  int conversion = 0;

  if( ap != 0 ) {                                     // allele found in list
    int newCount = ap->get_count() + 1;
    ap->set_count( newCount );                        // increment allele count
    conversion = ap->get_conversion();
  }

  else {
    Allele * lastNode;
    lastNode = l_arr[loc].get_last();
    int prevConv;

    if( lastNode != 0 )                               // calculate conversion for allele
      prevConv = lastNode->get_conversion();
    else 
      prevConv = 0;
    conversion = prevConv + 1;

    ap = new Allele( al_name, conversion, 1 );        // create new allele
    l_arr[loc].insert( ap );                          // add allele to list
  }

  return conversion;
}



/*------------------------  update_genotype()  -----------------------------
|   Updates the genotypes matrix with the integer representation           |
|   of the new allele.                                                     |
--------------------------------------------------------------------------*/

void update_genotype( int * g_arr[], int indiv, int loc, int al_num, int cv ) {
  if( al_num == 1)  
    g_arr[indiv][(loc * 2) - 1] = cv;
  else  
    g_arr[indiv][loc * 2] = cv;
}

 

/*-----------------------  check_maxAlleles()  -----------------------------
|   Updates the maximum number of alleles per locus.                       |
--------------------------------------------------------------------------*/

void check_maxAlleles( int & max, List l_arr[], int loc ) {
  Allele * lastNode;
  int num;

  lastNode = l_arr[loc].get_last();
  if( lastNode != 0 ) {
    num = lastNode->get_conversion();
    if( num > max )
      max = num;
  }
}



/*-------------------------  update_freqs()  -------------------------------
|   Updates the allele frequencies matrix with the frequency of each       |
|   allele at each locus.  The allele count for each locus is added to     |
|   column 0.                                                              |
--------------------------------------------------------------------------*/

void update_freqs( List l_arr[], double * f_arr[], int nLoc ) {
  for( int i = 1; i <= nLoc; i++ ) {
    double total = 0;
    int allele_num = 1;
    int nAlleles;
    Allele * current;

    for( current = l_arr[i].get_first(); current; current = current->get_next() ) {
      total += current->get_count();
    }

    for( current = l_arr[i].get_first(); current; current = current->get_next() ) {
      double freq;
      double ct = current->get_count();
      freq = ct / total;
      f_arr[i][allele_num] = freq;
      allele_num++;
    }
    nAlleles = allele_num - 1;
    f_arr[i][0] = nAlleles;               // col 0 = allele count for each locus
  }
}



/*----------------------  set_permutation()  -------------------------------
|  Designates alleles of individuals x and y as maternal or paternal       |
|  based on the value of p. Assumes that it is not known which alleles     |
|  are maternal/paternal.                                                  |
--------------------------------------------------------------------------*/

void set_permutation( int p, int & Xm, int & Xp, int & Ym, int & Yp,
                      int ind1, int ind2, int loc, int * g_arr[] ) {
  if( p < 2 ) {
    Xm = g_arr[ind1][(loc * 2) - 1];
    Xp = g_arr[ind1][loc * 2];
  }
  else {
    Xm = g_arr[ind1][loc * 2];
    Xp = g_arr[ind1][(loc * 2) - 1];
  }

  if( p == 0 || p == 2 ) {
    Ym = g_arr[ind2][(loc * 2) - 1];
    Yp = g_arr[ind2][loc * 2];
  }
  else {
    Ym = g_arr[ind2][loc * 2];
    Yp = g_arr[ind2][(loc * 2) - 1];
  }
}


 
/*----------------------  get_equation()  ----------------------------------
|  Chooses the likelihood equation based on identity between alleles of    |
|  individuals x and y.                                                    |
--------------------------------------------------------------------------*/

char get_equation( int Xm, int Xp, int Ym, int Yp ) {
  char ch;
  
  if( Xm == Ym) {
    if( Xp == Yp ) {
      ch = 'a';
    }
    else 
      ch = 'b';
    }
  else {
    if( Xp == Yp ) {
      ch = 'c';
    }
    else 
      ch = 'd';
    }
  
  return ch;
}



/*-------------------------  mSame_pSame()  --------------------------------
|   Calculates the likelihood of two individuals, x and y, being           |
|   related by the hypothesized relationship, Rm/Rp.                       |
|   This equation is used when the individuals have identical maternal     |
|   alleles and identical paternal alleles.                                |
--------------------------------------------------------------------------*/

double mSame_pSame( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp ) {
  double answer;
  answer =  (( Pxm * ( Rm + (( 1 - Rm ) * Pxm )))
           * ( Pxp * ( Rp + (( 1 - Rp ) * Pxp ))));
  return answer;
}



/*-------------------------  mSame_pDiff()  --------------------------------
|   Calculates the likelihood of two individuals, x and y, being           |
|   related by the hypothesized relationship, Rm/Rp.                       |
|   This equation is used when the individuals have identical maternal     |
|   alleles and different paternal alleles.                                |
--------------------------------------------------------------------------*/

double mSame_pDiff( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp, const double & Pyp ) {
  double answer;
  answer =   ( Pxm * ( Rm + ( 1 - Rm ) * Pxm ))
           * ( Pxp * ( 1 - Rp ) * Pyp );
  return answer;
}



/*-------------------------  mDiff_pSame()  --------------------------------
|   Calculates the likelihood of two individuals, x and y, being           |
|   related by the hypothesized relationship, Rm/Rp.                       |
|   This equation is used when the individuals have different maternal     |
|   alleles and identical paternal alleles.                                |
--------------------------------------------------------------------------*/

double mDiff_pSame( const double & Rm, const double & Rp,
                   const double & Pxm, const double & Pxp, const double & Pym ) {
  double answer;
  answer =   ( Pxm * ( 1 - Rm ) * Pym )
           * ( Pxp * ( Rp + ( 1 - Rp ) * Pxp ));
  return answer;
}



/*-------------------------  mDiff_pDiff()  --------------------------------
|   Calculates the likelihood of two individuals, x and y, being           |
|   related by the hypothesized relationship, Rm/Rp.                       |
|   This equation is used when the individuals have different maternal     |
|   alleles and different paternal alleles.                                |
--------------------------------------------------------------------------*/

double mDiff_pDiff( const double & Rm, const double & Rp, const double & Pxm,
                   const double & Pxp, const double & Pym, const double & Pyp ) {
  double answer;
  answer =   ( Pxm * ( 1 - Rm ) * Pym )
           * ( Pxp * ( 1 - Rp ) * Pyp );
  return answer;
}



/*------------------------  draw_from_pop()  -------------------------------
|   Draws an allele from the population based on the population allele     | 
|   frequencies in f_arr for locus loc.                                    |
--------------------------------------------------------------------------*/

int draw_from_pop( double * f_arr[], int loc ) {
  int result;
  double prob, freq, nAlleles;

  prob = rand() % 100 + 1;
  prob = prob * 0.01;
  freq = 0;
  nAlleles = f_arr[loc][0];                       // number of alleles at locus

  for( int i = 1; i <= nAlleles; i++ ) {          // look at each allele
    freq = freq + f_arr[loc][i];                  // cumulative freq of alleles so far
    if( (prob <= freq) || (i == nAlleles) ) {     // if prob <= cumulative freq
      result = i;                                 // .. assign allele i
      break;
    }
  }
  return result;
}



/*---------------------------  merge_sort()  -------------------------------
|   This is a recursive function that divides a double array into           |
|   sorted subarrays containing one element.  Merge_sort then calls        |
|   merge(), which combines the subarrays into a single sorted array.      |
--------------------------------------------------------------------------*/

void merge_sort( double A[], int p, int r ) {     // p = start index of array
  int q = (p + r) / 2;                            // r = end index of array
                                                  // q divides subarrays
  if( p == r )
    return;
  merge_sort( A, p, q );
  merge_sort( A, q + 1, r );
  merge( A, p, q, r );
}



/*------------------------------  merge()  ---------------------------------
|   This function merges two sorted subarrays, A[p..q] and A[q + 1..r],    |
|   into a single sorted array that replaces A[p..r].                      |
--------------------------------------------------------------------------*/

void merge( double A[], int p, int q, int r ) {
  int end1 = q + 1;
  int end2 = r + 1;
  double * T = new double[r - p + 1];
  int i = p;
  int j = q + 1;
  int k = 0;

  while( i < end1 && j < end2 ) {            // neither subarray is at end
    if( A[i] <= A[j] ) {
      T[k] = A[i];
      i++;
      k++;
    }
    else {                                   // A[i] > A[j]
      T[k] = A[j];
      j++;
      k++;
    }
  }

  if( i == end1 ) {                          // 1st subarray is at end
    while( j < end2 ) {                      // add rest of 2nd subarray to T
      T[k] = A[j];
      j++;
      k++;
    }
  }
  else {                                     // 2nd subarray is at end
    while( i < end1 ) {                      // add rest of 1st subarray to T
      T[k] = A[i];
      i++;
      k++;
    }
  }

  for( i = 0; i + p <= r; i++ ) {            // assign the sorted array into A
    A[i + p] = T[i];
  }
}



/*-------------------------  get_t2_error()  -------------------------------
|   Calculates the type II error rate for a given p value (x_val).         |
--------------------------------------------------------------------------*/

double get_t2_error( double ratios_arr[], double x_val ) {
  int count = 0;
  double t2_error;

  for( int i = 0; i < 1000; i++ ) {
    if( ratios_arr[i] <= x_val )
      count++;
    else
      break;
  }
  
  t2_error = count * 0.001;
  return t2_error;
}



/*-----------------------  connected_components()  -------------------------
|  Finds the connected components in a graph represented by the adjacency  |
|  list adjVert.  Argument n is the number of vertices in the graph, and   |
|  array cc, which is filled by the function, contains the component       |
|  number to which each vertex belongs.  The component number is the       |
|  number of some vertex in the component.                                 |
|  The total number of components is returned.                             |
--------------------------------------------------------------------------*/

int connected_components( intList adjVert[], int n, int cc[] ) {
  int * color = new int[n];
  int component = 1;

  for( int i = 0; i < n; i++ )                     // initialize array to 'white'
    color[i] = 0;
 
  for( int v = 0; v < n; v++ ) {
    if( color[v] == 0 ) {                          // 0 = white, undiscovered
      ccDFS(adjVert, color, v, component, cc);
      component++;
    }
  }

  return component - 1;
}



/*---------------------------  ccDFS()  ------------------------------------
|  This recursive function is called by connected_components() to find     |
|  connected components within a graph.                                    |
--------------------------------------------------------------------------*/

void ccDFS( intList adjVert[], int color[], int current, int component, int cc[] ) {
  int w;
  intList remAdj;
  IntNode * ip;

  color[current] = 1;                              // 1 = grey, visited but not finished
  cc[current] = component;
  remAdj = adjVert[current];

  while( remAdj.isEmpty() == false ) {
    ip = remAdj.get_first();
    w = ip->get_id_num();

    if( color[w] == 0 )
      ccDFS( adjVert, color, w, component, cc);

    remAdj.delete_first();                         // remAdj = tail(remAdj)
    color[current] = 2;                            // 2 = black, finished
  }
}



/*-------------------------  countSibs()  ----------------------------------
|  Counts the number of siblings in family ( the number of vertices in     |
|  a component).                                                           |
--------------------------------------------------------------------------*/

int countSibs( int famNum, int nIndiv, int cc[] ) {
  int count = 0;

  for( int i = 0; i < nIndiv; i++ ) {
    if( cc[i] == famNum )
      count++;
  }
  return count;
}



/*-------------------------  countEdges()  ---------------------------------
|  Counts the number of edges in family ( the number of edges in           |
|  a component).                                                           |
--------------------------------------------------------------------------*/

int countEdges( int famNum, int nIndiv, int cc[], intList adjVert[] ) {
  int nEdges = 0;

  for( int i = 0; i < nIndiv; i++ ) {
    if( cc[i] == famNum ) {
      nEdges += adjVert[i].countAdj(i);
    }
  }

  return nEdges;
}



/*-------------------------  getScore()  -----------------------------------
|  Calculates the score for a component/family.                            |
|  The score is the number of edges divided by the number of edges in a    |
|  complete graph with the same number of vertices.                        |
--------------------------------------------------------------------------*/

double getScore( int famNum, int nIndiv, int cc[], intList adjVert[] ) {
  int nVertices, nEdges;
  double complete, score = 0;

  nVertices = countSibs( famNum, nIndiv, cc );
  nEdges = countEdges( famNum, nIndiv, cc, adjVert );
  complete = (nVertices * (nVertices - 1)) / 2.0;

  if( complete != 0 )
    score = ( nEdges / complete ) * 100;

  return score;
}



/*-----------------------  scores_higher()  --------------------------------
|  Compares two scores to two other scores.  Returns true if s1 and s2     |
|  are determined to be better than s1best and s2best.                     |
--------------------------------------------------------------------------*/

bool scores_higher( double s1best, double s2best, double s1, double s2 ) {
  bool better;
  double higher, lower;

  if( s1best > s2best ) {
    higher = s1best;
    lower = s2best;
  }
  else {
    higher = s2best;
    lower = s1best;
  }

  if( s1 > s2 ) {
    if( ((s1 > higher) && (s2 > lower))
       || ((s1 > higher + 10) && (s2 > lower - 2)) )
      better = true;
    else better = false;   
  }
  else {
    if( ((s2 > higher) && (s1 > lower))
       || ((s2 > higher + 10) && (s1 > lower - 2)) )
      better = true;
    else better = false;
  }

  return better;
}



/*------------------------  number_for_ghc()  ------------------------------
|  Renumbers nodes in a family from 0 to n-1, as required for use in       |
|  ghc_tree().                                                             |
--------------------------------------------------------------------------*/

void number_for_ghc( int nIndiv, int famNum, int cc[], int ghc_num[] ) {
  int number = 0;

  for( int m = 0; m < nIndiv; m++ ) {
    if( cc[m] == famNum ) {
      ghc_num[m] = number;
      number++;
    }
  }
}



/*------------------------  convert_graph()  -------------------------------
|  Converts a family/connected component to the Graph data structure used  |
|  in creating the Gomory-Hu cut tree.                                     |
|                                                                          |
|  Undirected edges are represented by pairs of edge desriptors associated |
|  to incident nodes, the incidence list of a node is represented as a     |
|  one-way linked list of edge descriptors each of which contains a        |
|  pointer to an adjacent node and a "back" pointer to the descriptor of   |
|  the other edge in the pair contained in the incidence list of the       |
|  adjacent node.                                                          |
|                                                                          |
|  This code is a modification of code available at :                      |
|  http://elib.zib.de/pub/Packages/mathprog/mincut/all-pairs               |
--------------------------------------------------------------------------*/

void convert_graph( Graph * gr, intList adjVert[], int famNum, int n,
                    int m, int nIndiv, int cc[], int ghc_num[], string id_arr[] ) {
  int x, y, node1, node2;
  Edge * eptr1, * eptr2;
  Node * nptr1, * nptr2;
  bool adjacent;
  string node1_id, node2_id;

  // set node id's
  for( int i = 0; i < n; i++ )
    gr->nodes[i].set_id( i + 1 );

  x = 0;
  y = m;
  eptr1 = &(gr->edges[x]);
  eptr2 = &(gr->edges[y]);
  gr->set_n_edges0( m );

  // read in an edge (node1, node2)
  for( int j = 0; j < nIndiv; j++ ) {
    if( cc[j] == famNum ) {
      for( int k = j + 1; k < nIndiv; k++ ) {
        if( cc[k] == famNum ) {
          adjacent = adjVert[j].search( k );
          if( adjacent ) {
            node1 = ghc_num[j] + 1;
            node2 = ghc_num[k] + 1;
            node1_id = id_arr[j];
            node2_id = id_arr[k];

            // insert edge into data structure
            node1--;
            node2--;
            nptr1 = &(gr->nodes[node1]);
            nptr2 = &(gr->nodes[node2]);
            nptr1->set_trueID( node1_id );
            nptr2->set_trueID( node2_id );
            eptr1->set_adjacent( nptr2 );
            eptr2->set_adjacent( nptr1 );
            eptr1->set_cap( 1 );
            eptr2->set_cap( 1 );
            eptr1->set_back( eptr2 );
            eptr2->set_back( eptr1 );
            if( nptr1->get_first() == 0 ) {
              nptr1->set_first( eptr1 );
              eptr1->set_next( 0 );
            }
            else {
              eptr1->set_next( nptr1->get_first() );
              nptr1->set_first( eptr1 );
            }
            if( nptr2->get_first() == 0 ) {
              nptr2->set_first( eptr2 );
              eptr2->set_next( 0 );
            }
            else {
              eptr2->set_next( nptr2->get_first() );
              nptr2->set_first( eptr2 );
            }
            x++;
            y++;
            eptr1 = &(gr->edges[x]);
            eptr2 = &(gr->edges[y]);
          }//if adjacent
        }//if comp[k]
      }//for k
    }//if comp[j]
  }//for j
}



/*----------------------------  ghc_tree ()  -------------------------------
|  Determines Gomory/Hu cut tree for input graph with capacitated edges.   |
|  The tree structure is represented by parent pointers within the         |
|  node object.  The capacity of a tree edge is stored at the child node.  |
|  The root of the cut tree is the first node in the list of graph nodes   |
|  (&gr->nodes[0]). The implementation is described in [1].                |
|                                                                          |
|  References:                                                             |
|  ----------                                                              |
|  1) D. Gusfield: "Very Simple Algorithms and Programs for                |
|     All Pairs Network Flow Analysis", Computer Science Division,         |
|     University of California, Davis, 1987.                               |
|                                                                          |
|  2) R.E. Gomory and T.C. Hu: "Multi-Terminal Network Flows",             |
|         SIAM J. Applied Math. 9 (1961), 551-570.                         |
|                                                                          |
|  This code is a modification of code available at :                      |
|  http://elib.zib.de/pub/Packages/mathprog/mincut/all-pairs               |
--------------------------------------------------------------------------*/

bool ghc_tree( Graph * gr, Node * active[], int number[],
               int & max_dist, int & bound, bool & co_check ) {

  Node *nptr, *sptr, *tptr, *tparent;
  int n, m, maxfl, tmincap;

  n = gr->get_n_nodes();
  m = gr->get_n_edges();

  nptr = &(gr->nodes[0]);
  for( int i = 0; i < n; i++ )
    gr->nodes[i].set_parent( nptr );

  for( int j = 1; j < n; j++ ) {
    sptr = &(gr->nodes[j]);
    tptr = sptr->get_parent();
    maxfl = maxflow( gr, sptr, tptr, active, number, max_dist, bound, co_check );
    if( maxfl < 0 )
      return false;

    sptr->set_mincap( maxfl );

    for( int k = 1; k < n; k++ ) {
      nptr = &(gr->nodes[k]);
      if( (nptr != sptr) && (nptr->get_alive() == false)
                         && (nptr->get_parent() == tptr) ) {
        nptr->set_parent( sptr );
      }
    }
    if( tptr->get_parent()->get_alive() == false ) {
      tparent = tptr->get_parent();
      tmincap = tptr->get_mincap();
      sptr->set_parent( tparent );
      tptr->set_parent( sptr );
      sptr->set_mincap( tmincap );
      tptr->set_mincap( maxfl );
    }
  }

  compute_subtreeSizes( gr );

  return true;
}



/*---------------------------  maxflow()  ----------------------------------
|  Determines maximum flow and minimum cut between nodes                   |
|  s (= *s_ptr) and t (= *t_ptr) in an undirected graph.                   |
|                                                                          |
|  References:                                                             |
|  ----------                                                              |
|  A. Goldberg/ E. Tarjan: "A New Approach to the Maximum Flow Problem",   |
|  Proc. 18th ACM Symp. on Theory of Computing, 1986.                      |
|                                                                          |
|  This code is a modification of code available at :                      |
|  http://elib.zib.de/pub/Packages/mathprog/mincut/all-pairs               |
--------------------------------------------------------------------------*/

int maxflow( Graph * gr, Node * sptr, Node * tptr, Node * active[],
             int number[], int & max_dist, int & bound, bool & co_check ) {
  Node *aptr, *nptr, *q_front, *q_rear;
  Edge *eptr;
  int n, m, m0, level, n_discharge, cap, old_rcap, incre, dmin;

  // node ids range from 1 to n, node array indices range from 0 to n-1.

  n = gr->get_n_nodes();
  m = gr->get_n_edges();
  m0 = gr->get_n_edges0();
  co_check = true;

  initialize_maxflow( gr, active, number );

  // breadth first search to get exact distances from sink
  // and for test of graph connectivity

  tptr->set_dist( 0 );
  tptr->set_unmarked( false );
  q_front = tptr;
  q_rear = q_front;

 bfs_next:                                          // ---------- bfs_next ----------
  level = q_rear->get_dist() + 1;
  eptr = q_rear->get_first();
  while( eptr != 0 ) {
    if( (eptr->get_adjacent()->get_unmarked())
         && eptr->get_back()->get_rcap() > 0 ) {
      nptr = eptr->get_adjacent();
      nptr->set_unmarked( false );
      nptr->set_dist( level);
      number[level]++;
      q_front->set_bfs_link( nptr);
      q_front = nptr;
    }
    eptr = eptr->get_next();
  }
  if( q_rear == q_front )
    goto bfs_ready;

  q_rear = q_rear->get_bfs_link();
  goto bfs_next;

 bfs_ready:                                       // ---------- bfs_ready ----------
  if( co_check ) {
    co_check = false;
    for( int i = 0; i < n; i++ ) {
      nptr = &(gr->nodes[i]);
      if( nptr->get_unmarked() ) {
        cout << "Input graph not connected." << endl;
        exit (-5);
      }
    }
  }

  sptr->set_dist( n );             // number[0] and number[n] not required
  tptr->set_dist( 0 );
  tptr->set_excess( 1 );           // to be subtracted again

  // initial preflow push from source node

  max_dist = 0;                    // = max_dist of active nodes
  eptr = sptr->get_first();
  while( eptr != 0 ) {
    nptr = eptr->get_adjacent();
    cap = eptr->get_rcap();
    nptr->set_excess( nptr->get_excess() + cap );
    sptr->set_excess( sptr->get_excess() - cap );
    old_rcap = eptr->get_back()->get_rcap();
    eptr->get_back()->set_rcap( old_rcap + cap );
    eptr->set_rcap( 0 );

    if( (nptr != tptr) && (nptr->get_excess() <= cap) ) {
      // push node nptr onto stack for nptr->dist,
      // but only once in case of double edges
      nptr->set_stack_link( active[nptr->get_dist()] );
      active[nptr->get_dist()] = nptr;
      if( nptr->get_dist() > max_dist )
        max_dist = nptr->get_dist();
    }
      eptr = eptr->get_next();
  }

  sptr->set_alive( false );
  bound = n - 1;
  n_discharge = 0;

  // -----------  main loop -----------

  do {
    // get maximum distance active node
    aptr = active[max_dist];
    while( aptr != 0 ) {
      active[max_dist] = aptr->get_stack_link();
      eptr = aptr->get_scan();

      edge_scan:  // for current active node      ---------- edge_scan ----------
        nptr = eptr->get_adjacent();
        if( (nptr->get_dist() == aptr->get_dist() - 1)
            && eptr->get_rcap() > 0) {
          incre = aptr->get_excess();
          if( incre <= eptr->get_rcap() ) {
            // perform a non saturating push
            old_rcap = eptr->get_rcap();
            eptr->set_rcap( old_rcap - incre );
            old_rcap = eptr->get_back()->get_rcap();
            eptr->get_back()->set_rcap( old_rcap + incre );
            aptr->set_excess( 0 );
            nptr->set_excess( nptr->get_excess() + incre );
            if( nptr->get_excess() <= incre ) {
              // push nptr onto active stack
              nptr->set_stack_link( active[nptr->get_dist()] );
              active[nptr->get_dist()] = nptr;
            }
            aptr->set_scan( eptr);
            goto node_ready;
          }
          else {
            // perform a saturating push
            incre = eptr->get_rcap();
            old_rcap = eptr->get_back()->get_rcap();
            eptr->get_back()->set_rcap( old_rcap + incre );
            aptr->set_excess( aptr->get_excess() - incre );
            nptr->set_excess( nptr->get_excess() + incre );
            eptr->set_rcap( 0 );
            if( nptr->get_excess() <= incre ) {
              // push nptr onto active stack
              nptr->set_stack_link( active[nptr->get_dist()] );
              active[nptr->get_dist()] = nptr;
            }
            if( aptr->get_excess() <= 0 ) {
              aptr->set_scan( eptr );
              goto node_ready;
            }
          }
        }//if nptr
        if( eptr->get_next() == 0 ) {
          // all incident arcs scanned, but node still has positive
          // excess, check if for all nptr, nptr->dist != aptr->dist

          if( number[aptr->get_dist()] == 1 ) {
            // put all nodes v with dist[v] >= dist[a] into the set of
            // "dead" nodes since they are disconnected from the sink

            for( int i = 0; i < n; i++ ) {
              nptr = &(gr->nodes[i]);
              if( (nptr->get_alive()) && (nptr->get_dist() > aptr->get_dist()) ) {
                --number[nptr->get_dist()];
                active[nptr->get_dist()] = 0;
                nptr->set_alive( false );
                nptr->set_dist( n );
                --bound;
              }
            }
              --number[aptr->get_dist()];
              active[aptr->get_dist()] = 0;
              aptr->set_alive( false );
              aptr->set_dist( n );
              --bound;
              goto node_ready;
          }
          else {                          // determine new label value
            dmin = n;
            aptr->set_scan( 0 );
            eptr = aptr->get_first();
            while( eptr != 0 ) {
              if( (eptr->get_adjacent()->get_dist() < dmin)
                  && (eptr->get_rcap() > 0) ) {
                dmin = eptr->get_adjacent()->get_dist();
                if( aptr->get_scan() == 0 )
                  aptr->set_scan( eptr );
              }
              eptr = eptr->get_next();
            }
            if( ++dmin < bound ) {         // ordinary relabel operation
              --number[aptr->get_dist()];
              aptr->set_dist( dmin );
              ++number[dmin];
              max_dist = dmin;
              eptr = aptr->get_scan();
              goto edge_scan;
            }
            else {
              aptr->set_alive( false );
              --number[aptr->get_dist()];
              aptr->set_dist( n );
              --bound;
              goto node_ready;
            }
          }
        }
        else {
          eptr = eptr->get_next();
          goto edge_scan;
        }

        node_ready:                      // ---------- node_ready ----------
          ++n_discharge;
          if( n_discharge == n ) {
            n_discharge = 0;
            global_relabel( gr, tptr, active, number, max_dist, bound );
          }
          aptr = active[max_dist];
       }// while( aptr != 0 )

       // aptr != 0
       --max_dist;

     } while( max_dist > 0 );

  return( tptr->get_excess() - 1 );
}



/*-----------------------  initialize_maxflow()  ---------------------------
|  Initializes data structures and variables needed for max_flow().        |
--------------------------------------------------------------------------*/

void initialize_maxflow( Graph * gr, Node * active[], int number[] ) {
  Node * nptr;
  Edge * eptr;
  int n, m, m0;

  n = gr->get_n_nodes();
  m = gr->get_n_edges();
  m0 = gr->get_n_edges0();

  for( int i = 0; i < n; i++ ) {
    nptr = &(gr->nodes[i]);
    nptr->set_scan( nptr->get_first() );
    if( nptr->get_scan() == 0 ) {
      cout << "isolated node in input graph." << endl;
      exit (-4);
    }
    nptr->set_excess( 0 );
    nptr->set_stack_link( 0 );
    nptr->set_alive( true );
    nptr->set_unmarked( true );
  }

  for( int j = 0; j < m; j++ ) {
     eptr = &(gr->edges[j]);
     eptr->set_rcap( eptr->get_cap() );
  }

  for( int k = m0; k < (m0 + m); k++ ) {
     eptr = &(gr->edges[k]);
     eptr->set_rcap( eptr->get_cap() );
  }

  for( int p = 0; p <= n; p++ ) {
    number[p] = 0;
    active[p] = 0;
  }
}



/*------------------------  global_relabel()  ------------------------------
|  Breadth first search to get exact distance labels from sink with        |
|  reordering of stack of active nodes.                                    |
|                                                                          |
|  This code is a modification of code available at :                      |
|  http://elib.zib.de/pub/Packages/mathprog/mincut/all-pairs               |
--------------------------------------------------------------------------*/

void global_relabel( Graph * gr, Node * tptr, Node * active[],
                     int number[], int & max_dist, int & bound ) {
  Node *front, *rear, *nptr, **ptr;
  Edge *eptr;
  int n, level, count;

  n = gr->get_n_nodes();
  for( int i = 0; i < n; i++ ) {
    nptr = &(gr->nodes[i]);
    nptr->set_unmarked( true );
    nptr->set_stack_link( 0 );
    nptr->set_scan( nptr->get_first() );
  }
  tptr->set_unmarked( false );

  // initialize stack of active nodes
  for( int j = 0; j <= n; j++ ) {        // initialize stack of active nodes
    ptr = &(active[j]);
    *ptr = 0;
  }

  for( int i = 0; i <= n; i++)           // initialize number array
    number[i] = 0;

  max_dist = 0;
  count = 1;                             // number of alive nodes
  front = tptr;
  rear = front;

 bfs_next:                               // ---------- bfs_next ----------
  level = rear->get_dist() + 1;
  eptr = rear->get_first();
  while( eptr != 0 ) {
    nptr = eptr->get_adjacent();
    if( (nptr->get_alive()) && (nptr->get_unmarked())
                    && (eptr->get_back()->get_rcap() > 0) ) {
      nptr->set_unmarked( false );
      nptr->set_dist( level );
      count++;
      number[level]++;
      if( nptr->get_excess() > 0 ) {
        nptr->set_stack_link( active[level] );
        active[level] = nptr;
        max_dist = level;
      }
      front->set_bfs_link( nptr );
      front = nptr;
    }
    eptr = eptr->get_next();
  }
  if( front == rear )
    goto bfs_ready;

  rear = rear->get_bfs_link();
  goto bfs_next;

 bfs_ready:                              // ---------- bfs_ready ----------
  if( count < bound ) {
    // identify nodes that are marked alive but have
    // not been reached by BFS and mark them as dead.
    for( int i = 0; i < n; i++ ) {
      nptr = &(gr->nodes[i]);
      if( nptr->get_unmarked() && nptr->get_alive() ) {
        nptr->set_dist( n );
        nptr->set_alive( false );
      }
    }
    bound = count;
  }
}



/*------------------------  compute_subtreeSizes()  -----------------------
|  Sets the field subtreeSize for each node in the ghc tree.               |
--------------------------------------------------------------------------*/

void compute_subtreeSizes( Graph * gr ) {
  Node * nptr, * root_ptr;
  int n = gr->get_n_nodes();

  // temporarily set parent of root to null
  root_ptr = &(gr->nodes[0]);
  root_ptr->set_parent( 0 );

  for( int i = 1; i < n; i++ ) {
    nptr = &(gr->nodes[i]);
    percolate_increment( nptr, gr );
  }

  // restore parent of root to root
  root_ptr->set_parent( root_ptr );
}



/*------------------------  percolate_increment()  ------------------------
|  Helper function to compute_subtreeSize.                                 |
--------------------------------------------------------------------------*/

void percolate_increment( Node * nptr, Graph * gr ) {
  Node * parent_ptr;

  parent_ptr = nptr->get_parent();
  if( parent_ptr ) {
    parent_ptr->set_subtreeSize( parent_ptr->get_subtreeSize() + 1 );
    percolate_increment( parent_ptr, gr );
  }
}



/*--------------------------  find_mincut()  ------------------------------
|  Searches the ghc tree for the smallest feasible unchecked cut.          |
--------------------------------------------------------------------------*/

Node * find_mincut( Graph * gr ) {
  int n, cap, min;
  bool checked;
  Node * temp_ptr;
  double one_third;

  n = gr->get_n_nodes();
  min = gr->get_n_edges();
  one_third = n / 3.0;

  while( 1 ) {
    min = gr->get_n_edges();
    temp_ptr = 0;

    // search tree for smallest unchecked cut
    for( int i = 1; i < n; i++ ) {
      cap = gr->nodes[i].get_mincap();
      checked = gr->nodes[i].get_checked();
      if( (cap < min) && (!checked) ) {
        temp_ptr = &(gr->nodes[i]);
        min = cap;
      }
    }

    // no feasible unchecked cut
    if( temp_ptr == 0 ) {
      return temp_ptr;
    }

    // feasible cut found
    temp_ptr->set_checked( true );
    if( size_ok(temp_ptr, gr) ) {
      if( min < one_third )
        return temp_ptr;
    }
  }
}



/*----------------------------  size_ok()  --------------------------------
|  Returns true if sets S and T are both greater than MINSIZE.             |
--------------------------------------------------------------------------*/

bool size_ok( Node * cut_ptr, Graph * gr ) {
  int S_size, T_size;

  T_size = cut_ptr->get_subtreeSize();
  S_size = gr->get_n_nodes() - T_size;

  if( (S_size >= MINSIZE) && (T_size >= MINSIZE) )
    return true;
  else
    return false;
}



/*--------------------------  label_T_nodes()  ----------------------------
|  Labels each node in the ghc tree that belongs to set T.                 |
--------------------------------------------------------------------------*/

void label_T_nodes( Node * cut_ptr, Graph * gr, int T_size ) {
  int count;
  Node * parent_ptr, * nptr;
  int n = gr->get_n_nodes();

  for( int i = 0; i < n; i++ ) {
    nptr = &(gr->nodes[i]);
    nptr->set_memberT( false );
  }

  cut_ptr->set_memberT( true );
  count = 1;

  while( count < T_size ) {
    for( int i = 0; i < n; i++ ) {
      nptr = &(gr->nodes[i]);
      if( nptr->get_memberT() == false ) {
        parent_ptr = nptr->get_parent();
        if( parent_ptr->get_memberT() == true ) {     // if parent of i is in set T
          nptr->set_memberT( true );                  // then i is in set T
          count++;
        }
      }
    }
  }
}



/*----------------------  update_components ()  ---------------------------
|  Assigns all nodes belonging to set T to a new component.                |
--------------------------------------------------------------------------*/

void update_components( int nIndiv, int famNum, int nFams,
                        int cc[], int ghc_nodes[], Graph * gr ) {
  int ghc_num;

  for( int i = 0; i < nIndiv; i++ ) {
    if( cc[i] == famNum ) {
      ghc_num = ghc_nodes[i];
      if( gr->nodes[ghc_num].get_memberT() )
        cc[i] = nFams + 1;
    }
  }
}



/*------------------------  mark_ST_edges()  -----------------------------
|  Labels all edges belonging to the ST cut.                              |
--------------------------------------------------------------------------*/

void mark_ST_edges( int nIndiv, int cc[], int T_id, int S_id, intList adjVert[] ) {

  for( int i = 0; i < nIndiv; i++ ) {
    if( cc[i] == T_id )
      adjVert[i].mark_STcut( T_id, cc );

    else if( cc[i] == S_id )
      adjVert[i].mark_STcut( S_id, cc );
  }
}



/*----------------------  restore_components()  ---------------------------
|  Assigns all nodes belonging to set T to their old component, since it   |
|  was not split into two families.                                        |
--------------------------------------------------------------------------*/

void restore_components( int nIndiv, int cc[], int nFams, int old_famNum  ) {

  for( int i = 0; i < nIndiv; i++ ) {
    if( cc[i] == nFams + 1 )
      cc[i] = old_famNum;
  }
}



/*--------------------------  print_tree()  -------------------------------
|  Prints the ghc tree.                                                    |
--------------------------------------------------------------------------*/

void print_tree( Graph * gr ) {
  Node * nptr;
  int n = gr->get_n_nodes();

  cout << "\nGomory-Hu cut tree:" << endl;
  cout << "root  " << gr->nodes[0].get_trueID() << endl;

  for( int i = n-1; i > 0; i-- ) {
    nptr = &(gr->nodes[i]);
    // nodes numbered 1..n
    cout << "e  " << nptr->get_trueID() << " -> " << nptr->get_parent()->get_trueID()
         << "    " << nptr->get_mincap() << endl;
  }

  for( int i = 0; i < n; i++ ) {
    nptr = &(gr->nodes[i]);
    cout << "size of subtree with root " << nptr->get_trueID()
         << " : " << nptr->get_subtreeSize() << endl;
  }
}
