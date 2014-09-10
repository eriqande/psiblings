// Pseudo likelihood estimation for Eric and Carlo's Big Creek data
// Version without "truncaton" (ability to handle small pairwise likelihood ratios approximately)   

DATA_SECTION

  init_int n			// Number of individuals
  init_int S			// Number of SNP loci
  init_imatrix A(1,n,1,S)	// Allele 1
  init_imatrix B(1,n,1,S)	// Allele 2
  init_vector position(1,n) 	// Position measured from joining point of the three reaches
  init_vector reach(1,n) 	// Reach indicator 1,2,3    
  init_ivector max_allele(1,S)	// Highest allel number per locus
  init_matrix f(1,S,1,max_allele)// Allele frequencies
  init_number eps_error		// Genotyping probability in E.A Thompson simple model
  
  int N 			// Number of number of pairwise comparisons
  !!N = n*(n-1)/2;

  // Arrays containing quantities that can be pre-calculated (before parameter estimation)
  vector L1(1,N)		// Likelihood ratio for full sibs (versus unrelated)
  vector L2(1,N)		// Likelihood ratio for half sibs (versus unrelated)
  vector L3(1,N)		// Likelihood ratio for parent-offspring (versus unrelated)
  vector d(1,N)			// Physical distance (in creek) between members of pair
  number mean_d			// Average distance

  ivector I(1,N)		// Number in dataset of first pair member
  ivector J(1,N)		// Number in dataset of second pair member

  // Change the input data file
  !!USER_CODE ad_comm::change_datafile_name("phase_pseudo.dat");
  init_int fit_reg		// 1 = Include regression parameters, 0 = do not include
  init_int include_PO		// 1 = Include offspring relationships, 0 = do not include
  init_int print_posterior	// 1 = print posterior probabiliteis in the rep section

PARAMETER_SECTION

  // Full sibling parameters
  init_bounded_number a(-50,0.0,1)				// Intercept
  init_bounded_number b(-10,5.0,3*fit_reg-1)   			// Slope of distance

  // Half-sib parameters
  init_bounded_number a_hs(-50,0.0,1)				// Intercept
  init_bounded_number b_hs(-10,5.0,3*fit_reg-1)   		// Slope of distance

  // Parent offsping parameters
  init_bounded_number a_po(-50,0.0,2*include_PO-1)		// Intercept: has to be low,
  init_bounded_number b_po(-20,20.0,3*include_PO*fit_reg-1)   	// Slope of distance

  objective_function_value ll

PRELIMINARY_CALCS_SECTION

  if(!((eps_error>=0.0)&&(eps_error<1.0)))
    exit(1);			// Error rate not in valid range
    
  int i,j,s,ii=0;

  // Pairwise comparisons
  for(i=2;i<=n;i++)
    for(j=1;j<i;j++)
    {

      // Initialization of likelihood ratios (relative to unrelated)
      double tmpL1 = 1.0;  		// Full siblings
      double tmpL2 = 1.0;  		// Half siblings
      double tmpL3 = 1.0;		// Parent-offspring  

      // Loop over SNP loci
      for(s=1;s<=S;s++)
      {
        if(!((A(i,s)==0)||(A(j,s)==0))) // Checks for missing values in either individual (Assumes same zero pattern in A and B)
        {
          // Likehood ratios given Z (number of ibd alleles) realitive to Z=0   
          double R1 = 1.0/f(s,A(i,s))*((A(i,s)==A(j,s))+(A(i,s)==B(j,s)) ) 
          	      + 1.0/f(s,B(i,s))*((B(i,s)==A(j,s))+(B(i,s)==B(j,s)));		// z=1 (ie. 1 alleles ibd)  
          double R2 = 1.0/(f(s,A(i,s))*f(s,B(i,s))) * ((A(i,s)==A(j,s))*(B(i,s)==B(j,s)) 
          	      + (B(i,s)==A(j,s))*(A(i,s)==B(j,s)));				// z=2 (ie. 2 alleles ibd)

          // Simple model for genotyping error (E.A. Thompson)
          R1 = (1.0-eps_error)*R1 + eps_error;
          R2 = (1.0-eps_error)*R2 + eps_error;
          
          // Likelihood ratio given type of relationship (relative to unrelated)
          tmpL1 *= .25 + .5*R1/4.0 + .25*R2/2.0;		// Full siblings
          tmpL2 *= .5 + .5*R1/4.0;				// Half siblings
          tmpL3 *= R1/4.0;					// Parent-offspring
          //cout << R1 << endl;
        }
      }
   //exit(1);
      // Put stuff into the array now that we have accumulated over marker
      ii++;
      L1(ii) = tmpL1;
      L2(ii) = tmpL2;
      L3(ii) = tmpL3;
      I(ii) = i;
      J(ii) = j;
      if(reach(i)==reach(j))
        d(ii) = fabs(position(i)-position(j));
      else
        d(ii) = position(i)+position(j);
    }
    
    mean_d = sum(d)/N;		// Centering  d's gives better convergence
 
PROCEDURE_SECTION

  int i, ii;

  ll = 0.0;

  // Terms in the pseudo likelihood which are handled exactly (without approximation)
  for(ii=1;ii<=N;ii++)
  {
    double dd = d(ii)-mean_d;
    dvariable P = mfexp(a+b*dd);
    P /= 1+P;
    dvariable P_hs = mfexp(a_hs+b_hs*dd);
    P_hs /= 1+P_hs;
    dvariable P_po = mfexp(a_po+b_po*dd);
    P_po /= 1+P_po;

    ll += log(1.0-(P+P_hs+P_po) + P*L1(ii) + P_hs*L2(ii) + P_po*L3(ii) + 1.e-10);
  }

  ll = -ll;		// ADMB does function minimization


REPORT_SECTION

  if(print_posterior==1)
  {
    for(int ii=1;ii<=N;ii++)
    {
      dvariable P = mfexp(a+b*d(ii));
      P /= 1+P;
      dvariable P_hs = mfexp(a_hs+b_hs*d(ii));
      P_hs /= 1+P_hs;
      dvariable P_po = mfexp(a_po+b_po*d(ii));
      P_po /= 1+P_po;
      dvariable P0 = 1.0-(P+P_hs+P_po);
      dvariable L0 = 1.0;			//  Likelihood ratio for unrelated is 1
      
      dvariable tot = P0*L0 + P*L1(ii) +  P_hs*L2(ii) +  P_po*L3(ii);   // Pr(Data|unrelated)

      report << I(ii) << " "  << J(ii) << " " << P0*L0/tot << " " <<  P*L1(ii)/tot   << " " 
      				<< P_hs*L2(ii)/tot << " " << P_po*L3(ii)/tot << endl;
    }
  }

