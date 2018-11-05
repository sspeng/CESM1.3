
#include <slave.h>
#include "dma_macros.h"
#include "cpe_print.h"
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var)) 
#define LOCAL(p,a,len1,b,len2) (((void*)p)+(a)*(len1)+((b)*(len2)*8))     
#define timelevels 3   
#define np 4 
void gradient_sphere(double s[np][np],double d_Dvv[np][np],double Dinv[2][2][np][np],double ds[2][np][np],double rrearth){
	int i,j,l;
	double dsdx00,dsdy00,v1[np][np],v2[np][np];
	for(j=0;j<np;j++){
		for(l=0;l<np;l++){
			dsdx00=0.0;
			dsdy00=0.0;
			for(i=0;i<np;i++){
				dsdx00 = dsdx00 +d_Dvv[l][i]*s[j][i];
				dsdy00 = dsdy00 +d_Dvv[l][i]*s[i][j];
			}
			v1[j][l]= dsdx00*rrearth;
			v2[l][j]= dsdy00*rrearth;
		}
	}
	for(j=0;j<np;j++){
		for(i=0;i<np;i++){
			ds[0][j][i]=Dinv[0][0][j][i]*v1[j][i]+Dinv[0][1][j][i]*v2[j][i];
			ds[1][j][i]=Dinv[1][0][j][i]*v1[j][i]+Dinv[1][1][j][i]*v2[j][i];
		}
	}
}
void vorticity_sphere(double v[2][np][np],double e_D[2][2][np][np],
 		double d_Dvv[np][np],double e_rmetdet[np][np],double rrearth,
 		double vort[np][np]){
 	int i,j,l;
 	double dvdx00,dudy00;
 	double vco[2][np][np],vtemp[np][np];
 	for(j=0;j<np;j++)
 		for(i=0;i<np;i++){
 			vco[0][j][i]=(e_D[0][0][j][i]*v[0][j][i]+e_D[0][1][j][i]*v[1][j][i]);
 			vco[1][j][i]=(e_D[1][0][j][i]*v[0][j][i]+e_D[1][1][j][i]*v[1][j][i]);
 		}
 	for(j=0;j<np;j++)
 		for(l=0;l<np;l++){
 			dudy00=0.0;
 			dvdx00=0.0;
 			for(i=0;i<np;i++){
 				dvdx00 = dvdx00 + d_Dvv[l][i]*vco[1][j][i];
 				dudy00 = dudy00 + d_Dvv[l][i]*vco[0][i][j];
 			}
 			vort[j][l]=dvdx00;
 			vtemp[l][j]=dudy00;
		}
	for(j=0;j<np;j++)
		for(i=0;i<np;i++){

			vort[j][i]=(vort[j][i]-vtemp[j][i])*(e_rmetdet[j][i]*rrearth);
	}
}
void divergence_sphere(double v[2][np][np],double d_Dvv[np][np],double e_metdet[np][np],
		double e_Dinv[2][2][np][np],double e_rmetdet[np][np],double div[np][np],double rrearth){
	//local variable
	int i,j,l;
	double dudx00,dvdy00,gv[2][np][np],vvtemp[np][np];
	for(j=0;j<np;j++)
		for(i=0;i<np;i++){
			gv[0][j][i]=e_metdet[j][i]*(e_Dinv[0][0][j][i]*v[0][j][i]+e_Dinv[1][0][j][i]*v[1][j][i]);
			gv[1][j][i]=e_metdet[j][i]*(e_Dinv[0][1][j][i]*v[0][j][i]+e_Dinv[1][1][j][i]*v[1][j][i]);
		}
	//compute d/dx and d/dy
	for(j=0;j<np;j++)
		for(l=0;l<np;l++){
			dudx00=0.0;
			dvdy00=0.0;
			for(i=0;i<np;i++){
				dudx00=dudx00+d_Dvv[l][i]*gv[0][j][i];
				dvdy00=dvdy00+d_Dvv[l][i]*gv[1][i][j];
			}
			div[j][l]=dudx00;
			vvtemp[l][j]=dvdy00;
		}
	for(j=0;j<np;j++)
		for(i=0;i<np;i++){
			div[j][i]=(div[j][i]+vvtemp[j][i])*(e_rmetdet[j][i]*rrearth);
		}
}
typedef struct elem_accum_t{
    double KEvert1[np][np];                           // term from continuity equ
    double KEvert2[np][np];                           // term from momentum equ
    double IEvert1[np][np];                           // term from continuity equ
    double IEvert2[np][np];                           // term from T equ
    double IEvert1_wet[np][np];                       // wet term from continuity equ
    double IEvert2_wet[np][np];                       // wet term from T equ

    double KEhorz1[np][np];                           // at time t
    double KEhorz2[np][np];                           // after calling time_advance, these will be at time t-1
    double IEhorz1[np][np];
    double IEhorz2[np][np];
    double IEhorz1_wet[np][np];
    double IEhorz2_wet[np][np];

    double T1[np][np];
    double T2[np][np];
    double T2_s[np][np];
    double S1[np][np];
    double S1_wet[np][np];
    double S2[np][np];
  }elem_accum_t;
typedef struct toal{
    // elem_state
    double *v_n0;                           //[nlev][2][np][np];
    double *v_np1;
    double *v_nm1;
    double *T_n0;                           //[nlev][np][np];
     double *T_np1;
      double *T_nm1;
    double *dp3d_n0;                        //[[nlev][np][np];
    double *dp3d_np1;
    double *dp3d_nm1;
    double *ps_v;                        //[np][np];
    double *phis;                        //[np][np];
    double *Q;                           // [1][nlev][np][np];
    double *Qdp;                         //[2][1][nlev][np][np];
            //elm_derived
    double *vn0;                         //[nlev][2][np][np];
    double *phi;                          //[nlev][np][np];
    double *omega_p;                      //[nlev][np][np];
    double *eta_dot_dpdn;                 //[nlev+1][np][np];
    double *pecnd;                        //[nlev][np][np];
    //elem_accun_t
    double *elem_accum_t_first;            //17*4*4  //17个数组一起传回，不用传入get也要放再判断里//生成17*4*4的数组；分17次给每个数组赋值 
    double *CONV;                          //nlev*2*np*np
    //elem
    double *metdet;                        //np*np
    double *rmetdet;                       //np*np
    double *D;                             //[2][2][np][np]
    double *Dinv;                          //[2][2][np][np]
    double *spheremp;                      //[np][np];
    double *rspheremp;                     //[np][np];
    double *fcor;                          //[np][np];
    //hvcoord_T
    double *hyai;                          //[nlev+1];
    double *hyam;                          //[nlev];
    double *hybi;                          //[nlev+1];
    double *hybm;                          //[nlev];
    //deriv
    double *Dvv;                           //np*np;
	double ps0;
    double dt2;
    double eta_ave_w;
    double kappa,Rgas,Rw_R,rrearth;
    int n0,np1,nm1,nlev,nets,nete,qn0;
    int compute_diagnostics;
    int rsplit,use_cpstar;
    int len_ele;                        //ele长度
    int slave_flag; 
}toal;
                               
void slave_compute_and_apply_rhs_(toal *toalll)
{
    int my_id;
    toal toall;
    elem_accum_t accum;
    my_id = athread_get_id(-1);
     dma_init();
  if(my_id==0)
  {
  bcast_get(toalll, &toall,sizeof(toal));
  dma_syn();
  }
  athread_syn(ARRAY_SCOPE,0xffff);
  
    int ie,i,j,k,l;
    int nets=toall.nets-1;
    int nete=toall.nete;
    int nlev=toall.nlev;
    int len_ele=toall.len_ele;
    double dt2 =toall.dt2;
    int nlev_temp,nlev_s,nlev_r;
    int n0=toall.n0-1;
    int np1=toall.np1-1;
    int nm1=toall.nm1-1;
    int use_cpstar=toall.use_cpstar;
    int slave_flag=toall.slave_flag;
    int compute_diagnostics=toall.compute_diagnostics;
    double eta_ave_w=toall.eta_ave_w;
    int qn0;
    if(toall.qn0>0)
        qn0=toall.qn0-1;
    else
        qn0=toall.qn0;
    int rsplist=toall.rsplit;
 
    nlev_temp=nlev%8;
    if(my_id%8<nlev_temp)  {
      nlev_s = (my_id%8)*((int)(nlev/8)+1);
      nlev_r = nlev_s+((int)(nlev/8)+1);
    }
    else {
      nlev_s=(nlev%8)*((int)(nlev/8)+1)+(my_id%8 - nlev_temp)*(int)(nlev/8);
      nlev_r = nlev_s+((int)(nlev/8));
    }
    double comm_buff[32];
    double p_dp[2][4][4];
    int pp=0;
    int nlev_range = nlev_r -nlev_s;
    //init_ const
     double Rgas=toall.Rgas;
    double Cp= 1005.0;
    double Cpwater_vapor= 1870.0;
    double kappa= toall.kappa;
	double Rw_R= toall.Rw_R;
    double sdot_sum[np][np];
    double  rearth= 6376000;
    double   rrearth= toall.rrearth;
    //int
    double ps0=toall.ps0;
    double hyai[nlev+1];
    double hybi[nlev+1];
    double hyam[nlev];
    double hybm[nlev];
    double ps_v[3][np][np];
    double grad_ps[2][np][np];
    double dp[nlev_range][np][np];
    double dp3d[3][nlev_range][np][np];
    double p[nlev_range][np][np];
    double grad_p[nlev_range][2][np][np];
    double v1,v2;
    double vgrad_p[nlev_range][np][np];
    double vtemp[2][np][np];
    double rdp[nlev_range][np][np];
    double v[3][nlev_range][2][np][np];
    double vn0[nlev_range][2][np][np];
    double T[3][nlev_range][np][np];
    double T_v[nlev_range][np][np];
    double Qdp[nlev_range][np][np];
    double kappa_star[nlev_range][np][np];
    double Qt;
    double ele_eta_dot_dpdn[nlev_range+1][np][np];
    double eta_dot_dpdn[nlev_range+1][np][np];
    double omega_p[nlev_range][np][np];
    double ele_omega_p[nlev_range][np][np];
    double T_vadv[nlev_range][np][np];
    double v_vadv[nlev_range][2][np][np];
    double divdp[nlev_range][np][np];
    double glnps1,glnps2,gpterm;
    double vtens1[nlev_range][np][np];
    double vtens2[nlev_range][np][np];
    double ttens[nlev_range][np][np];
    double pecnd[nlev_range][np][np];
    double phi[nlev_range][np][np];
    double fcor[np][np];
    double vort[nlev_range][np][np];
    double Ephi[np][np];
    double Q[nlev][np][np];
    double E;
    double phis[np][np];
    double vgrad_T[np][np];
   // elem_accum_t accum;
    double CONV[nlev_range][2][np][np];
    double spheremp[np][np];
    double rspheremp[np][np];
    double facq,facm;
    double hkk,hkl;
    double ckk,ckl,term;
    double suml[np][np];
    double phii[nlev][np][np];
    double T_reg[np][np];
    double v_reg[2][np][np];
    double metdet[np][np];
    double rmetdet[np][np];
	double D[2][2][np][np];
    double Dinv[2][2][np][np];
    double Dvv[np][np];
	double com_short[7][np][np];
    double g_s[np][np];
	double g_me[np][np];
	double g_rme[np][np];
	double g_D[2][2][np][np];
	double g_Dvv[np][np];
    double g_ds[2][np][np];
	double g_rre;
    double Q_p[nlev_range][np][np];
    double de,cp2,cp_ratio;
	int r;
	  pe_get(toall.hyai,hyai, sizeof(double)*(nlev+1));
   	  pe_get(toall.hybi,hybi, sizeof(double)*(nlev+1));
	  pe_get(toall.hyam,hyam, sizeof(double)*(nlev));
	  pe_get(toall.hybm,hybm, sizeof(double)*(nlev));
   dma_syn(); 


    ie=(int)(my_id/8);
    for(;ie<nete-nets;ie=ie+8){
    if(slave_flag==1){	
	  pe_get(LOCAL(toall.ps_v,ie,len_ele,0,0),ps_v, sizeof(double)*np*np*3);
	  
      
	  pe_get(LOCAL(toall.spheremp,ie,len_ele,0,0),spheremp, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.phis,ie,len_ele,0,0),phis, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.rspheremp,ie,len_ele,0,0),rspheremp, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.metdet,ie,len_ele,0,0),metdet, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.rmetdet,ie,len_ele,0,0),rmetdet, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.fcor,ie,len_ele,0,0),fcor, sizeof(double)*np*np);
	  pe_get(LOCAL(toall.Dvv,0,0,0,0),Dvv, sizeof(double)*np*np);
	  

	  pe_get(LOCAL(toall.D,ie,len_ele,0,0),D,sizeof(double)*np*np*2*2);
	  
	   	//把只有np*np的大小的数据放在同一个数组内，数组为k(个数)维，一起取，再分开赋值
	   pe_get(LOCAL(toall.Dinv,ie,len_ele,0,0),Dinv, sizeof(double)*2*np*np*2);
	  pe_get(LOCAL(toall.Qdp,ie,len_ele,nlev_s,np*np),&Qdp,sizeof(double)*(nlev_range)*np*np);
	  pe_get(LOCAL(toall.eta_dot_dpdn,ie,len_ele,nlev_s,np*np),&ele_eta_dot_dpdn, sizeof(double)*(nlev_range+1)*np*np);
      pe_get(LOCAL(toall.omega_p,ie,len_ele,nlev_s,np*np),&ele_omega_p, sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.phi,ie,len_ele,nlev_s,np*np),&phi, sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.pecnd,ie,len_ele,nlev_s,np*np),&pecnd, sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.Q,ie,len_ele,nlev_s,np*np),&Q, sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.vn0,ie,len_ele,nlev_s,np*np*2),&vn0, sizeof(double)*(nlev_range)*np*np*2);
    	
		pe_get(LOCAL(toall.v_n0,ie,len_ele,nlev_s,np*np*2),&v[n0][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
    	pe_get(LOCAL(toall.v_np1,ie,len_ele,nlev_s,np*np*2),&v[np1][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
    	pe_get(LOCAL(toall.v_nm1,ie,len_ele,nlev_s,np*np*2),&v[nm1][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
    	pe_get(LOCAL(toall.T_n0,ie,len_ele,nlev_s,np*np),&T[n0][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.T_np1,ie,len_ele,nlev_s,np*np),&T[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.T_nm1,ie,len_ele,nlev_s,np*np),&T[nm1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.dp3d_n0,ie,len_ele,nlev_s,np*np),&dp3d[n0][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.dp3d_np1,ie,len_ele,nlev_s,np*np),&dp3d[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.dp3d_nm1,ie,len_ele,nlev_s,np*np),&dp3d[nm1][0][0][0], sizeof(double)*(nlev_range)*np*np);
       
		dma_syn();
		
		
		if(rsplist==0)
		{
		  for(i=0;i<4;i++)
          for(j=0;j<4;j++)
	         g_s[i][j]=ps_v[n0][i][j];
		      gradient_sphere(g_s,Dvv ,Dinv,g_ds,rrearth); 
	        for(l=0;l<2;l++)
		    for(i=0;i<4;i++)
			  for(j=0;j<4;j++)
			    grad_ps[l][i][j]=g_ds[l][i][j]; 
		} 
    
		if(rsplist==0){
       for(k=0;k<nlev_range;k++) 
        for(i=0;i<np;i++)
        for(j=0;j<np;j++){
            dp[k][i][j] = hyai[nlev_s+k+1]*ps0 + hybi[nlev_s+k+1]*ps_v[n0][i][j] - hyai[nlev_s+k]*ps0 - hybi[nlev_s+k]*ps_v[n0][i][j];  
            p[k][i][j] =hyam[nlev_s+k]*ps0 + hybm[nlev_s+k+1]*ps_v[n0][i][j];
            for (l=0;l<2;l++)
              grad_p[k][l][i][j] = hybm[nlev_s+k]*grad_ps[l][i][j];
          }
		}
		else { 
      for(k=0;k<nlev_range;k++)
        for(i=0;i<np;i++)
          for(j=0;j<np;j++)
           dp[k][i][j] = dp3d[n0][k][i][j];
      if(my_id%8==0) {
        for(i=0;i<np;i++)
          for(j=0;j<np;j++)
	     p[0][i][j] = hyai[0]*ps0+dp[0][i][j]/2.0;
           
        for(k=1;k<nlev_range;k++)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              p[k][i][j]=p[k-1][i][j]+dp[k-1][i][j]/2.0+dp[k][i][j]/2.0;
    pp=0;                  
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=p[nlev_range-1][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=dp[nlev_range-1][i][j];
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)+1);   
     }
      else{           
        for(pp=0;pp<32;pp++)    
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            p_dp[0][i][j]=comm_buff[pp++];
                        
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            p_dp[1][i][j]=comm_buff[pp++];
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            p[0][i][j]=p_dp[0][i][j]+p_dp[1][i][j]/2.0+dp[0][i][j]/2.0;
        for(k=1;k<nlev_range;k++)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              p[k][i][j]=p[k-1][i][j]+dp[k-1][i][j]/2.0+dp[k][i][j]/2.0;
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=p[nlev_range-1][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=dp[nlev_range-1][i][j];    
        if(my_id%8!=7)
          {
            for(pp=0;pp<32;pp++)                
              REG_PUTR(comm_buff[pp],(my_id%8)+1);                  
          }
      }
 }
   if(rsplist!=0){ 
    for(k=0;k<nlev_range;k++)
    { 
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
	         g_s[i][j]=p[k][i][j];
		      gradient_sphere(g_s,Dvv ,Dinv,g_ds,rrearth); 
	        for(l=0;l<2;l++)
		    for(i=0;i<4;i++)
			  for(j=0;j<4;j++)
			    grad_p[k][l][i][j]=g_ds[l][i][j]; 
		    
		    }
   }
	 for(k=0;k<nlev_range;k++){
		 for(i=0;i<np;i++)
			 for(j=0;j<np;j++){
				 rdp[k][i][j]=1.0/dp[k][i][j];
				 v1 = v[n0][k][0][i][j];
				 v2 = v[n0][k][1][i][j];
				 vgrad_p[k][i][j] = v1*grad_p[k][0][i][j] + v2*grad_p[k][1][i][j];
				 vtemp[0][i][j] = v1*dp[k][i][j];
				 vtemp[1][i][j] = v2*dp[k][i][j]; 
				 for(l=0;l<2;l++)
					 vn0[k][l][i][j]=vn0[k][l][i][j]+ eta_ave_w*vtemp[l][i][j]; 

			 }  
         
		 for(i=0;i<np;i++)
			 for(j=0;j<np;j++)
			 {
			     g_me[i][j]=metdet[i][j];
				 g_rme[i][j]=rmetdet[i][j];
		       g_Dvv[i][j]=Dvv[i][j];
			 }
		 for(l=0;l<2;l++)
		    for(r=0;r<2;r++)
				for(i=0;i<np;i++)
					 for(j=0;j<np;j++)
				g_D[l][r][i][j]=Dinv[l][r][i][j];
		 g_rre=rrearth; 
			  
				 divergence_sphere(vtemp,g_Dvv,g_me,g_D,g_rme,g_s,g_rre); 
		 for(i=0;i<np;i++)
			 for(j=0;j<np;j++)
				 divdp[k][i][j]=g_s[i][j];
		 for(l=0;l<2;l++)
			 for(i=0;i<4;i++)
				 for(j=0;j<4;j++)
					 g_ds[l][i][j]=v[n0][k][l][i][j];

		 for(i=0;i<np;i++)
			 for(j=0;j<np;j++)
			 {
			     g_me[i][j]=metdet[i][j];
				 g_rme[i][j]=rmetdet[i][j];
		       g_Dvv[i][j]=Dvv[i][j];
			 }
		 for(l=0;l<2;l++)
		    for(r=0;r<2;r++)
				for(i=0;i<np;i++)
					 for(j=0;j<np;j++)
				g_D[l][r][i][j]=D[l][r][i][j];
		 g_rre=rrearth; 
		 vorticity_sphere(g_ds,g_D,g_Dvv,g_rme,g_rre,g_s);

		 for(i=0;i<np;i++)
			 for(j=0;j<np;j++)
				 vort[k][i][j]=g_s[i][j];

	 }


		
	for(k=0;k<nlev_range;k++){
      for(i=0;i<np;i++)
        for(j=0;j<np;j++){
          if(qn0==-1){
            T_v[k][i][j] = T[n0][k][i][j];
            kappa_star[k][i][j] = kappa;
          }  
          else{
            Qt = Qdp[k][i][j]/dp[k][i][j]; 
         //   Q_p[k][i][j] = 1.00+(Rw_R-1.00)*Qt; 
            T_v[k][i][j] = (T[n0][k][i][j])*((1.00) +((Rw_R-1.00)*Qt));
           //t_v_compute_(&T[n0][k][i][j], &Qt, &T_v[k][i][j]);
			if (use_cpstar==1)
              kappa_star[k][i][j] =  Rgas/Cp*(1.000 + (Cpwater_vapor/Cp - 1.000)*Qt);
            else
              kappa_star[k][i][j] = kappa;
          }
                                   
        }        
    }

//		 if(my_id==0&&ie==0) cpe_printf("qn0%.17lf\n T%.17lf\n Qdp%.17lf\n dp%.17lf\n T_v%.17lf\n",qn0,T[n0][0][0][0],Qdp[0][0][0],dp[0][0][0],T_v[0][0][0]);

    ///////////////preq_hydrostatic
    if(my_id%8==7) {
      for(i=0;i<np;i++)
        for(j=0;j<np;j++){
                
          hkk = dp[nlev_range-1][i][j]*0.5/p[nlev_range-1][i][j];
          hkl=2*hkk;
          phii[nlev_range-1][i][j] = Rgas*T_v[nlev_range-1][i][j]*hkl;
          phi[nlev_range-1][i][j]= phis[i][j] + Rgas*T_v[nlev_range-1][i][j]*hkk;
        }
      for(k=nlev_range-2;k>=0;k--)
        for(i=0;i<4;i++)
          for(j=0;j<4;j++){
            hkk = dp[k][i][j]*0.5/p[k][i][j];
            hkl=2*hkk;  
            phii[k][i][j] = phii[k+1][i][j]+Rgas*T_v[k][i][j]*hkl;
            phi[k][i][j]= phis[i][j] +phii[k+1][i][j]+ Rgas*T_v[k][i][j]*hkk;
          }
      pp=0;
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=phii[0][i][j];
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=0;
      for(pp=0;pp<32;pp++)
        REG_PUTR(comm_buff[pp],(my_id%8)-1);
    }
    else{
      for(pp=0;pp<32;pp++)
                
        REG_GETR(comm_buff[pp]);

      pp=0;
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          {
            hkk = dp[nlev_range-1][i][j]*0.5/p[nlev_range-1][i][j];
            hkl=2*hkk;  
            phii[nlev_range-1][i][j] = comm_buff[pp]+Rgas*T_v[nlev_range-1][i][j]*hkl;
            phi[nlev_range-1][i][j]= phis[i][j] +comm_buff[pp]+ Rgas*T_v[nlev_range-1][i][j]*hkk;
            pp++;
        
          }
      for(k=nlev_range-2;k>=0;k--)
        for(i=0;i<4;i++)
          for(j=0;j<4;j++){
            hkk = dp[k][i][j]*0.5/p[k][i][j];
            hkl=2*hkk;  
            phii[k][i][j] = phii[k+1][i][j]+Rgas*T_v[k][i][j]*hkl;
            phi[k][i][j]= phis[i][j] +phii[k+1][i][j]+ Rgas*T_v[k][i][j]*hkk;
          }
      pp=0;
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=phii[0][i][j];
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=0;
            
      if(my_id%8!=0)
        {for(pp=0;pp<32;pp++)
            REG_PUTR(comm_buff[pp],(my_id%8)-1);
        }

    }
    ///////////////preq_omaga_ps
    if(my_id%8==0) {
      for(i=0;i<np;i++)
        for(j=0;j<np;j++)
          {
            ckk=0.5/p[0][i][j];
            term=divdp[0][i][j];
            omega_p[0][i][j] = vgrad_p[0][i][j]/p[0][i][j];
            omega_p[0][i][j] = omega_p[0][i][j] - ckk*term;
            suml[i][j] = term;
                                
          }
      for(k=1;k<nlev_range;k++)
        for(i=0;i<np;i++)
          for(j=0;j<4;j++){
            ckk=0.5/p[k][i][j];
            ckl=2*ckk;
            term=divdp[k][i][j];
            omega_p[k][i][j] = vgrad_p[k][i][j]/p[k][i][j];
            omega_p[k][i][j] = omega_p[k][i][j]-ckl*suml[i][j] - ckk*term;
            suml[i][j] = suml[i][j]+term;
          } 

      pp=0;
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=suml[i][j];
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=0;
      for(pp=0;pp<32;pp++)
        REG_PUTR(comm_buff[pp],(my_id%8)+1);
    }
     
    else{
      for(pp=0;pp<32;pp++)
        REG_GETR(comm_buff[pp]);
      pp=0;
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          suml[i][j]=comm_buff[pp++];
      pp=0;
      for(k=0;k<nlev_range;k++)
        for(i=0;i<4;i++)
          for(j=0;j<4;j++){
            ckk=0.5/p[k][i][j];
            ckl=2*ckk;
            term=divdp[k][i][j];
            omega_p[k][i][j] = vgrad_p[k][i][j]/p[k][i][j];
            omega_p[k][i][j] = omega_p[k][i][j]-ckl*suml[i][j] - ckk*term;
            suml[i][j] = suml[i][j]+term;
          }
      pp=0;                      
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=suml[i][j];
      for(i=0;i<4;i++)
        for(j=0;j<4;j++)
          comm_buff[pp++]=0;
      if(my_id%8!=7)
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)+1);                  
    }
    
     //////////////
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        sdot_sum[i][j]=0;
    if(rsplist>0){

      for(k=0;k<nlev_range+1;k++)
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            eta_dot_dpdn[k][i][j]=0;
      for(k=0;k<nlev_range;k++)
        for(i=0;i<4;i++)
          for(j=0;j<4;j++){
            T_vadv[k][i][j]=0;
            v_vadv[k][0][i][j]=0;
            v_vadv[k][1][i][j]=0;}
      //        if(my_id==0) printf("4");
              
    }
    else{      
      if(my_id%8==0)
        {   for(i=0;i<np;i++)
            for(j=0;j<np;j++)
              eta_dot_dpdn[0][i][j]=0;
        } //first do 
      if(my_id%8==0){
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            sdot_sum[i][j]=0;

        for(k=0;k<nlev_range;k++)
          for(i=0;i<np;i++)
            for(j=0;j<np;j++){
              sdot_sum[i][j]= sdot_sum[i][j] + divdp[k][i][j];
              eta_dot_dpdn[k+1][i][j]= sdot_sum[i][j];  
                                 
              //    eta_dot_dpdn[k+1][i][j] = hybi[nlev_s+k+1]*sdot_sum[i][j] - eta_dot_dpdn[k+1][i][j]; 
            }
        //  if(my_id<8) printf("aaa%daaa%.17lf\n",my_id,sdot_sum[1][1]);
        // if(my_id<8) printf("aaa%daaa%.17lf\n",my_id,sdot_sum[3][0]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=sdot_sum[i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=eta_dot_dpdn[nlev_range][i][j];
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)+1);  
      }
                
      else {
        for(pp=0;pp<32;pp++)    
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            sdot_sum[i][j]=comm_buff[pp++];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            eta_dot_dpdn[0][i][j]=comm_buff[pp++];

    
        for(k=0;k<nlev_range;k++)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++){
              sdot_sum[i][j]= sdot_sum[i][j] + divdp[k][i][j];
              eta_dot_dpdn[k+1][i][j]= sdot_sum[i][j];  
              //    eta_dot_dpdn[k+1][i][j] = hybi[nlev_s+k+1]*sdot_sum[i][j] - eta_dot_dpdn[k+1][i][j]; 
            }
        // if(my_id<8) printf("aaa%daaa%.17lf\n",my_id,sdot_sum[1][1]);
        // if(my_id<8) printf("bbb%daaa%.17lf\naaaa%.17lf\n",my_id,divdp[0][3][3],sdot_sum[3][3]);


        if(my_id%8!=7){
                
          pp=0;
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=sdot_sum[i][j];
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=eta_dot_dpdn[nlev_range][i][j];
          for(pp=0;pp<32;pp++)
            REG_PUTR(comm_buff[pp],(my_id%8)+1);  
        }
      }
      /////second do
      // if(my_id==0) printf("aaa%.17lf\n",eta_dot_dpdn[2][0][0]);
      if(my_id%8==7){
        //              if(my_id==7)    printf("%.17lf\n",sdot_sum[3][3]);
        for(k=-1;k<nlev_range;k++)
          for(i=0;i<np;i++)
            for(j=0;j<np;j++)
              eta_dot_dpdn[k+1][i][j] = hybi[nlev_s+k+1]*sdot_sum[i][j] - eta_dot_dpdn[k+1][i][j]; 
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=sdot_sum[i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=0;
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)-1);  
      }
      else{
        for(pp=0;pp<32;pp++)    
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            sdot_sum[i][j]=comm_buff[pp++];
        for(k=-1;k<nlev_range;k++)
          for(i=0;i<np;i++)
            for(j=0;j<np;j++)
              eta_dot_dpdn[k+1][i][j] = hybi[nlev_s+k+1]*sdot_sum[i][j] - eta_dot_dpdn[k+1][i][j]; 
        if(my_id%8!=0){
          
          for(pp=0;pp<32;pp++)
            REG_PUTR(comm_buff[pp],(my_id%8)-1);  
        }          
      }
                
                
      if(my_id%8==7)
        for(i=0;i<np;i++)
          for(j=0;j<np;j++)
            eta_dot_dpdn[nlev_range][i][j] =0.0;
    ////////////////preq_vertadv
      ////first do
      if(my_id%8==0) {
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            facq=0.5*rdp[0][i][j]*eta_dot_dpdn[1][i][j];
            T_vadv[0][i][j]  = facq*(T[n0][1][i][j]- T[n0][0][i][j]);
            v_vadv[0][0][i][j] = facq*(v[n0][1][0][i][j]- v[n0][0][0][i][j]);
            v_vadv[0][1][i][j] = facq*(v[n0][1][1][i][j]- v[n0][0][1][i][j]);
          } 

        for(k=1;k<nlev_range;k++)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++){
              facm=0.5*rdp[k][i][j]*eta_dot_dpdn[k][i][j];
              T_vadv[k][i][j]  = facm*(T[n0][k][i][j]- T[n0][k-1][i][j]);
              v_vadv[k][0][i][j] = facm*(v[n0][k][0][i][j]- v[n0][k-1][0][i][j]);
              v_vadv[k][1][i][j] = facm*(v[n0][k][1][i][j]- v[n0][k-1][1][i][j]);
            } 
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=T[n0][nlev_range-1][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=0;
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)+1);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=v[n0][nlev_range-1][0][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=v[n0][nlev_range-1][1][i][j];
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)+1);
      }
      else{
        for(pp=0;pp<32;pp++)
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            T_reg[i][j]=comm_buff[pp++];
        for(pp=0;pp<32;pp++)
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            v_reg[0][i][j]=comm_buff[pp++];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            v_reg[1][i][j]=comm_buff[pp++];
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++){
            facm=0.5*rdp[0][i][j]*eta_dot_dpdn[0][i][j];
            T_vadv[0][i][j]  = facm*(T[n0][0][i][j]- T_reg[i][j]);
            v_vadv[0][0][i][j] = facm*(v[n0][0][0][i][j]- v_reg[0][i][j]);
            v_vadv[0][1][i][j] = facm*(v[n0][0][1][i][j]- v_reg[1][i][j]);
          } 
        for(k=1;k<nlev_range;k++)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++){
              facm=0.5*rdp[k][i][j]*eta_dot_dpdn[k][i][j];
              T_vadv[k][i][j]  = facm*(T[n0][k][i][j]- T[n0][k-1][i][j]);
              v_vadv[k][0][i][j] = facm*(v[n0][k][0][i][j]- v[n0][k-1][0][i][j]);
              v_vadv[k][1][i][j] = facm*(v[n0][k][1][i][j]- v[n0][k-1][1][i][j]);
            } 
        if(my_id%8!=7){
          pp=0;
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=T[n0][nlev_range-1][i][j];
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=0;
          for(pp=0;pp<32;pp++)
            REG_PUTR(comm_buff[pp],(my_id%8)+1);
          pp=0;
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=v[n0][nlev_range-1][0][i][j];
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              comm_buff[pp++]=v[n0][nlev_range-1][1][i][j];
          for(pp=0;pp<32;pp++)
            REG_PUTR(comm_buff[pp],(my_id%8)+1);                 
        }       
       // if(my_id==7) printf("uuu%.17lf\n",rdp[k-1][0][0]);
      }    
      ///second
      if(my_id%8==7) {
        k=nlev_range-1;
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            facm = 0.5*rdp[k][i][j]*eta_dot_dpdn[k][i][j];
            T_vadv[k][i][j]  = T_vadv[k][i][j]+facm*(T[n0][k][i][j]- T[n0][k-1][i][j]);
            v_vadv[k][0][i][j] = v_vadv[k][0][i][j]+facm*(v[n0][k][0][i][j]- v[n0][k-1][0][i][j]);
            v_vadv[k][1][i][j] = v_vadv[k][0][i][j]+facm*(v[n0][k][1][i][j]- v[n0][k-1][1][i][j]);
          }
        for(k=nlev_range-2;k>=0;k--)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++){
              facq=0.5*rdp[k][i][j]*eta_dot_dpdn[k+1][i][j];
              T_vadv[k][i][j]  = T_vadv[k][i][j]+facq*(T[n0][k+1][i][j]- T[n0][k][i][j]);
              v_vadv[k][0][i][j] = v_vadv[k][0][i][j]+facq*(v[n0][k+1][0][i][j]- v[n0][k][0][i][j]);
              v_vadv[k][1][i][j] = v_vadv[k][1][i][j]+facq*(v[n0][k+1][1][i][j]- v[n0][k][1][i][j]);
            }
        // if(my_id==7) printf("app%.17lf\n",T[n0][nlev_range-1][0][0]); 
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=T[n0][0][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=0;
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)-1);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=v[n0][0][0][i][j];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            comm_buff[pp++]=v[n0][0][1][i][j];
        for(pp=0;pp<32;pp++)
          REG_PUTR(comm_buff[pp],(my_id%8)-1);
      }
      else{
        for(pp=0;pp<32;pp++)
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            T_reg[i][j]=comm_buff[pp++];
        for(pp=0;pp<32;pp++)
          REG_GETR(comm_buff[pp]);
        pp=0;
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            v_reg[0][i][j]=comm_buff[pp++];
        for(i=0;i<4;i++)
          for(j=0;j<4;j++)
            v_reg[1][i][j]=comm_buff[pp++];
        pp=0;
        k=nlev_range-1;

     //   if(my_id==6)printf("T26%.17lf\n",T[n0][k][0][0]);
      //  if(my_id==6)printf("T27%.17lf\n",T_reg[0][0]);
       // if(my_id==6)printf("eta%.17lf\n",eta_dot_dpdn[k+1][0][0]);
        //if(my_id==6)printf("rdp%.17lf\n",rdp[k][0][0]);

        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            facq = 0.5*rdp[k][i][j]*eta_dot_dpdn[k+1][i][j];
            T_vadv[k][i][j]  =T_vadv[k][i][j]+ facq*(T_reg[i][j]-T[n0][k][i][j]);
            v_vadv[k][0][i][j] = v_vadv[k][0][i][j]+facq*(v_reg[0][i][j]-v[n0][k][0][i][j]);
            v_vadv[k][1][i][j] = v_vadv[k][1][i][j]+facq*(v_reg[0][i][j]-v[n0][k][1][i][j]);
          }
       // if(my_id==6)printf("T_vadv27%.17lf\n",T_vadv[k][0][0]);
        for(k=nlev_range-2;k>=0;k--)
          for(i=0;i<4;i++)
            for(j=0;j<4;j++){
              facq=0.5*rdp[k][i][j]*eta_dot_dpdn[k+1][i][j];
              T_vadv[k][i][j]  = T_vadv[k][i][j]+facq*(T[n0][k+1][i][j]- T[n0][k][i][j]);
              v_vadv[k][0][i][j] = v_vadv[k][0][i][j]+facq*(v[n0][k+1][0][i][j]- v[n0][k][0][i][j]);
              v_vadv[k][1][i][j] = v_vadv[k][1][i][j]+facq*(v[n0][k+1][1][i][j]- v[n0][k][1][i][j]);
            }

             
        if(my_id%8!=0)
          { 
            pp=0;
            for(i=0;i<4;i++)
              for(j=0;j<4;j++)
                comm_buff[pp++]=T[n0][0][i][j];
            for(i=0;i<4;i++)
              for(j=0;j<4;j++)
                comm_buff[pp++]=0;
            for(pp=0;pp<32;pp++)
              REG_PUTR(comm_buff[pp],(my_id%8)-1);
            pp=0;
            for(i=0;i<4;i++)
              for(j=0;j<4;j++)
                comm_buff[pp++]=v[n0][0][0][i][j];
            for(i=0;i<4;i++)
              for(j=0;j<4;j++)
                comm_buff[pp++]=v[n0][0][1][i][j];
            for(pp=0;pp<32;pp++)
              REG_PUTR(comm_buff[pp],(my_id%8)-1);
          }

      }
      /////////////end          preq_vertadv
    }                                
    	
		 for(k=0;k<=nlev_range;k++)
      for(i=0;i<np;i++)
        for(j=0;j<np;j++){
          ele_eta_dot_dpdn[k][i][j]=ele_eta_dot_dpdn[k][i][j]+eta_ave_w*eta_dot_dpdn[k][i][j];
                  
        }
     for(k=0;k<nlev_range;k++)
      for(i=0;i<np;i++)
        for(j=0;j<np;j++)
		ele_omega_p[k][i][j]=ele_omega_p[k][i][j]+eta_ave_w*omega_p[k][i][j]; 
		
     for(k=0;k<nlev_range;k++){       
       for(i=0;i<np;i++)
       for(j=0;j<np;j++){
       v1 = v[n0][k][0][i][j];
       v2 = v[n0][k][1][i][j];
       E = 0.5*( v1*v1 + v2*v2 );
       Ephi[i][j]=E+phi[k][i][j]+pecnd[k][i][j];
       }
	   for(i=0;i<4;i++)
          for(j=0;j<4;j++)
	         g_s[i][j]=T[n0][k][i][j];
		      gradient_sphere(g_s,Dvv ,Dinv,vtemp,rrearth);
      
       for(i=0;i<np;i++)
       for(j=0;j<np;j++){
       v1 = v[n0][k][0][i][j];
       v2 = v[n0][k][1][i][j];
       vgrad_T[i][j] = v1*vtemp[0][i][j]+v2*vtemp[1][i][j];
       }

        gradient_sphere(Ephi,Dvv ,Dinv,vtemp,rrearth);
       
		for(i=0;i<np;i++)
       for(j=0;j<np;j++){
       gpterm = T_v[k][i][j]/p[k][i][j];
       glnps1 = Rgas*gpterm*grad_p[k][0][i][j];
       glnps2 = Rgas*gpterm*grad_p[k][1][i][j];
       v1 = v[n0][k][0][i][j];
       v2 = v[n0][k][1][i][j];
                   
       vtens1[k][i][j] = - v_vadv[k][0][i][j] + v2*(fcor[i][j] + vort[k][i][j]) - vtemp[0][i][j] - glnps1;                        
       vtens2[k][i][j] = - v_vadv[k][1][i][j] - v1*(fcor[i][j] + vort[k][i][j]) - vtemp[1][i][j] - glnps2;            
       ttens[k][i][j]  = - T_vadv[k][i][j] - vgrad_T[i][j] + kappa_star[k][i][j]*T_v[k][i][j]*omega_p[k][i][j];
       }
	 }
		

 if (compute_diagnostics){
      for(i=0;i<np;i++)
        for(j=0;j<np;j++){//conv alone put 
                                
          accum.KEhorz1[i][j]=0;
          accum.KEhorz2[i][j]=0;
          accum.IEhorz1[i][j]=0;
          accum.IEhorz2[i][j]=0;
          accum.IEhorz1_wet[i][j]=0;
          accum.IEhorz2_wet[i][j]=0;
          accum.KEvert1[i][j]=0;
          accum.KEvert2[i][j]=0;
          accum.IEvert1[i][j]=0;
          accum.IEvert2[i][j]=0;
          accum.IEvert1_wet[i][j]=0;
	      accum.IEvert2_wet[i][j]=0;
          accum.T1[i][j]=0;
          accum.T2[i][j]=0;
          accum.T2_s[i][j]=0;
          accum.S1[i][j]=0;
          accum.S1_wet[i][j]=0;
          accum.S2[i][j]=0;
        }
      for(i=0;i<np;i++)
        for(j=0;j<np;j++){
          accum.S2[i][j] = accum.S2[i][j] -  sdot_sum[i][j]*phis[i][j];
        }
      for(k=0;k<nlev_range;k++){
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            v1 = v[n0][k][0][i][j];
            v2 = v[n0][k][1][i][j];
            Ephi[i][j]=0.5*( v1*v1 + v2*v2 );
          }
        gradient_sphere(Ephi,Dvv ,Dinv,vtemp,rrearth);
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            v1 = v[n0][k][0][i][j];
            v2 = v[n0][k][1][i][j];
            accum.KEhorz2[i][j] = accum.KEhorz2[i][j] + (v1*vtemp[0][i][j]  + v2*vtemp[1][i][j] )*dp[k][i][j];
            accum.KEhorz1[i][j] = accum.KEhorz1[i][j] + Ephi[i][j]*divdp[k][i][j];
            accum.IEhorz1[i][j] = accum.IEhorz1[i][j] + Cp*T[n0][k][i][j]*divdp[k][i][j];
          }
          for(i=0;i<4;i++)
          for(j=0;j<4;j++)
	         g_s[i][j]=phi[k][i][j];
        gradient_sphere(g_s,Dvv ,Dinv,vtemp,rrearth);
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            v1 = v[n0][k][0][i][j];
            v2 = v[n0][k][1][i][j];
            E = 0.5*( v1*v1 + v2*v2 );
          
        de =  eta_dot_dpdn[k+1][i][j]-eta_dot_dpdn[k][i][j];

            accum.IEvert1[i][j]=accum.IEvert1[i][j] + Cp*T[n0][k][i][j]*de;
            accum.KEvert1[i][j]=accum.KEvert1[i][j] + E*de;
            accum.IEvert2[i][j]=accum.IEvert2[i][j] + Cp*T_vadv[k][i][j]*dp[k][i][j];
            accum.KEvert2[i][j]=accum.KEvert2[i][j] + (v1*v_vadv[k][0][i][j] + v2*v_vadv[k][1][i][j]) *dp[k][i][j];
            if (use_cpstar==1)
              accum.IEvert2_wet[i][j]=accum.IEvert2_wet[i][j] +(Cpwater_vapor-Cp)*Q[k][i][j]*T_vadv[k][i][j]*dp[k][i][j];

            gpterm = T_v[k][i][j]/p[k][i][j];
            accum.T1[i][j] = accum.T1[i][j] - Rgas*gpterm*(grad_p[k][0][i][j]*v1 + grad_p[k][1][i][j]*v2)*dp[k][i][j];
            accum.T2[i][j] = accum.T2[i][j] - (vtemp[0][i][j]*v1 + vtemp[1][i][j]*v2)*dp[k][i][j];
            accum.S1[i][j] = accum.S1[i][j] + Rgas*T_v[k][i][j]*omega_p[k][i][j]*dp[k][i][j];
            if (use_cpstar==1) {
              cp2 = (Cpwater_vapor-Cp)*Q[k][i][j];
              cp_ratio = cp2/(Cp+cp2);
              accum.S1_wet[i][j] = accum.S1_wet[i][j] +cp_ratio*(Rgas*T_v[k][i][j]*omega_p[k][i][j]*dp[k][i][j]);
            }
            for(l=0;l<2;l++)
              CONV[k][l][i][j]=-Rgas*gpterm*grad_p[k][l][i][j]-vtemp[l][i][j];
          }
        gradient_sphere(phis,Dvv ,Dinv,vtemp,rrearth);
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            v1 = v[n0][k][0][i][j];
            v2 = v[n0][k][1][i][j];
            accum.T2_s[i][j] = accum.T2_s[i][j] -  (vtemp[0][i][j]*v1 + vtemp[1][i][j]*v2)*dp[k][i][j];
          }
          for(i=0;i<4;i++)
          for(j=0;j<4;j++)
	         g_s[i][j]=T[n0][k][i][j];
        gradient_sphere(g_s,Dvv ,Dinv,vtemp,rrearth);
        for(i=0;i<np;i++)
          for(j=0;j<np;j++){
            v1 = v[n0][k][0][i][j];
            v2 = v[n0][k][1][i][j];
            accum.IEhorz2[i][j] = accum.IEhorz2[i][j] + Cp*(v1*vtemp[0][i][j] + v2*vtemp[1][i][j])*dp[k][i][j];
            if (use_cpstar==1)
              accum.IEhorz2_wet[i][j] = accum.IEhorz2_wet[i][j] +  (Cpwater_vapor-Cp)*Q[k][i][j]*(v1*vtemp[0][i][j] + v2*vtemp[1][i][j])*dp[k][i][j];
          }
      }
    }
    
     if(dt2<0){
	   	    for(k=0;k<nlev_range;k++){
				for(i=0;i<np;i++)
	               for(j=0;j<np;j++){
	               	    v[np1][k][0][i][j] = spheremp[i][j]*vtens1[k][i][j];
                        v[np1][k][1][i][j]= spheremp[i][j]*vtens2[k][i][j];
                        T[np1][k][i][j] = spheremp[i][j]*ttens[k][i][j];
	                    if (rsplist>0)
                        dp3d[np1][k][i][j] = -spheremp[i][j]*(divdp[k][i][j]+ eta_dot_dpdn[k+1][i][j]-eta_dot_dpdn[k][i][j]);
				   }
	        }
	        for(i=0;i<np;i++)
	            for(j=0;j<np;j++){
				    ps_v[np1][i][j] = -spheremp[i][j]*sdot_sum[i][j];
				}
	}
	else{
	
			for(k=0;k<nlev_range;k++){
				for(i=0;i<np;i++)
	               for(j=0;j<np;j++){ 
	               	    v[np1][k][0][i][j] = spheremp[i][j]*( v[nm1][k][0][i][j] + dt2*vtens1[k][i][j]);
                        v[np1][k][1][i][j]= spheremp[i][j]*( v[nm1][k][1][i][j] + dt2*vtens2[k][i][j]);
                        T[np1][k][i][j] = spheremp[i][j]*(T[nm1][k][i][j] + dt2*ttens[k][i][j]);
                        if (rsplist>0)
                           dp3d[np1][k][i][j] = spheremp[i][j]*(dp3d[nm1][k][i][j]-dt2*(divdp[k][i][j] + eta_dot_dpdn[k+1][i][j]-eta_dot_dpdn[k][i][j]));
				   }
	        }
	        for(i=0;i<np;i++)
	            for(j=0;j<np;j++)
				    ps_v[np1][i][j] = spheremp[i][j]*(ps_v[nm1][i][j] - dt2*sdot_sum[i][j]);
	   }   
    


	
    	pe_put(LOCAL(toall.dp3d_np1,ie,len_ele,nlev_s,np*np),&dp3d[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_put(LOCAL(toall.v_np1,ie,len_ele,nlev_s,np*np*2),&v[np1][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
     	pe_put(LOCAL(toall.T_np1,ie,len_ele,nlev_s,np*np),&T[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
		dma_syn();

    		if(my_id%8==0){
		  pe_put(LOCAL(toall.eta_dot_dpdn,ie,len_ele,nlev_s,np*np),&ele_eta_dot_dpdn, sizeof(double)*(nlev_range+1)*np*np);
  		  dma_syn();
			 }
			else{
    	  pe_put(LOCAL(toall.eta_dot_dpdn,ie,len_ele,nlev_s+1,np*np),&ele_eta_dot_dpdn[1][0][0], sizeof(double)*(nlev_range)*np*np);
	  	 dma_syn();}
   	  
		 
		pe_put(LOCAL(toall.ps_v,ie,len_ele,np1,np*np),&ps_v[np1][0][0], sizeof(double)*np*np);
	    pe_put(LOCAL(toall.omega_p,ie,len_ele,nlev_s,np*np),&ele_omega_p, sizeof(double)*(nlev_range)*np*np);
      	pe_put(LOCAL(toall.vn0,ie,len_ele,nlev_s,np*np*2),&vn0, sizeof(double)*(nlev_range)*np*np*2);

      dma_syn(); 
    
    
#ifdef ENERGY_DIAGNOSTICS
   pe_put((void*)(toall.elem_accum_t_first)+(ie*8*sizeof(accum))+((my_id%8)*sizeof(accum)),&accum,sizeof(accum)); 
    pe_put(LOCAL(toall.CONV,ie,7680,nlev_s,np*np*2),&CONV[0][0][0][0], sizeof(double)*np*np*(nlev_range)*2);
	dma_syn();
#endif

   }

  		else{
  		pe_get(LOCAL(toall.rspheremp,ie,len_ele,0,0),rspheremp, sizeof(double)*np*np);
		pe_get(LOCAL(toall.v_np1,ie,len_ele,nlev_s,np*np*2),&v[np1][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
     	pe_get(LOCAL(toall.T_np1,ie,len_ele,nlev_s,np*np),&T[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_get(LOCAL(toall.dp3d_np1,ie,len_ele,nlev_s,np*np),&dp3d[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	dma_syn();  
    		for(k=0;k<nlev_range;k++){
				for(i=0;i<np;i++)
	               for(j=0;j<np;j++){ 
	               	    v[np1][k][0][i][j] = rspheremp[i][j]*v[np1][k][0][i][j];
                        v[np1][k][1][i][j]= rspheremp[i][j]*v[np1][k][1][i][j];
                        T[np1][k][i][j] = rspheremp[i][j]*T[np1][k][i][j];
                        if (rsplist>0)
                           dp3d[np1][k][i][j] = rspheremp[i][j]*dp3d[np1][k][i][j];
                        else
                            ps_v[np1][i][j]=rspheremp[i][j]*ps_v[np1][i][j];
				   }
	        }
	    pe_put(LOCAL(toall.dp3d_np1,ie,len_ele,nlev_s,np*np),&dp3d[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
    	pe_put(LOCAL(toall.v_np1,ie,len_ele,nlev_s,np*np*2),&v[np1][0][0][0][0], sizeof(double)*(nlev_range)*np*np*2);
     	pe_put(LOCAL(toall.T_np1,ie,len_ele,nlev_s,np*np),&T[np1][0][0][0], sizeof(double)*(nlev_range)*np*np);
		pe_put(LOCAL(toall.ps_v,ie,len_ele,np1,np*np),&ps_v[np1][0][0], sizeof(double)*np*np);
		dma_syn();
	    
   }
  
 }
}
