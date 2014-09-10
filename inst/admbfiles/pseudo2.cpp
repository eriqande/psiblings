#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pseudo2.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  S.allocate("S");
  A.allocate(1,n,1,S,"A");
  B.allocate(1,n,1,S,"B");
  position.allocate(1,n,"position");
  reach.allocate(1,n,"reach");
  max_allele.allocate(1,S,"max_allele");
  f.allocate(1,S,1,max_allele,"f");
  eps_error.allocate("eps_error");
N = n*(n-1)/2;
  L1.allocate(1,N);
  L2.allocate(1,N);
  L3.allocate(1,N);
  d.allocate(1,N);
  I.allocate(1,N);
  J.allocate(1,N);
 ad_comm::change_datafile_name("phase_pseudo.dat");
  fit_reg.allocate("fit_reg");
  include_PO.allocate("include_PO");
  print_posterior.allocate("print_posterior");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  a.allocate(-50,0.0,1,"a");
  b.allocate(-10,5.0,3*fit_reg-1,"b");
  a_hs.allocate(-50,0.0,1,"a_hs");
  b_hs.allocate(-10,5.0,3*fit_reg-1,"b_hs");
  a_po.allocate(-50,0.0,2*include_PO-1,"a_po");
  b_po.allocate(-20,20.0,3*include_PO*fit_reg-1,"b_po");
  ll.allocate("ll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
 
}

void model_parameters::userfunction(void)
{
  ll =0.0;
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
