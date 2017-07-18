#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <limits>        // for std::numeric_limits<long double>::max()
#include <R.h>
#include <Rmath.h>
#include "matrix.h"
#include "rgwish.h"
#include "copula.h"

using namespace std;

extern "C" {
/*
 * birth-death MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for Bayesian model averaging
*/
void gcgm_DMH_bdmcmc_ma( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, counter;

	double Dsijj, Dsjj, Dsij, K022, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	double sum_weights = 0.0;
	long double rate, sum_rates; 

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );             // K[e, -e]
	vector<double> sigma11( 4 );            // sigma[e, e]
	vector<double> sigma12( p2x2 );         // sigma[e, -e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------|
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|		
		for( j = 1; j < dim; j++ )
		{			
			jj    = j * dim + j;
			Dsjj  = Ds[jj];
			Djj  = D[jj];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ij];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ij];
				Dijj  = - Dij * Dij / Djj;

				log_H_ij( &K[0], &sigma[0], &logH_ij, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dsijj, &Dsij, &Dsjj );

				log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dijj, &Dij, &Djj );

				rate = ( G[ij] ) ? exp( logH_ij - logI_p ) : exp( logI_p - logH_ij );
				
				rates[counter++] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
			}
		}
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}

		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ )
			{
				K_hat_Cpp[i]   += K[i] / sum_rates;
				if( G[i] ) p_links_Cpp[i] += 1.0 / sum_rates;
			}
			
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
   
/*
 * birth-death MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for maximum a posterior probability estimation (MAP)
*/
void gcgm_DMH_bdmcmc_map( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int counterallG = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, counter, size_sample_graph = *size_sample_g;
	bool this_one;

	double Dsijj, Dsjj, Dsij, sum_diag, K022, a11, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, nu_star, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];
	double sum_weights = 0.0;
	long double rate, sum_rates; 
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );        // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2xp2 );         // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings
	// ----------------------------------------

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------|
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|		
		for( j = 1; j < dim; j++ )
		{			
			jj    = j * dim + j;
			Dsjj  = Ds[jj];
			Djj  = D[jj];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ij];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ij];
				Dijj  = - Dij * Dij / Djj;

				log_H_ij( &K[0], &sigma[0], &logH_ij, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dsijj, &Dsij, &Dsjj );

				log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dijj, &Dij, &Djj );

				rate = ( G[ij] ) ? exp( logH_ij - logI_p ) : exp( logI_p - logH_ij );
				
				rates[counter] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
				
				char_g[counter] = G[ij] + '0'; 
				counter++; 
			}
		}
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}

		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) K_hat[i] += K[i] / sum_rates;	

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[counterallG] = 1.0 / sum_rates;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[counterallG];
					all_graphs[counterallG] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[counterallG];
				all_graphs[counterallG]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			counterallG++; 
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy(sample_graphs[i], qp, 0);
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;

	for( i = 0; i < pxp; i++ ) 
		K_hat[i] /= sum_weights;
}
           
/*
 * Multiple birth-death MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for Bayesian model averaging
*/
void gcgm_DMH_bdmcmc_ma_multi_update( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double D[], double Ds[], double *threshold, int *multi_update )
{
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int multi_update_C = *multi_update;
	
	int selected_edge_i, selected_edge_j, selected_edge_ij, counter;

	double Dsijj, Dsjj, Dsij, K022, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	double sum_weights = 0.0;
	long double rate, sum_rates; 

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );        // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2xp2 );         // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------|
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|		
		for( j = 1; j < dim; j++ )
		{			
			jj    = j * dim + j;
			Dsjj  = Ds[jj];
			Djj  = D[jj];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ij];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ij];
				Dijj  = - Dij * Dij / Djj;

				log_H_ij( &K[0], &sigma[0], &logH_ij, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dsijj, &Dsij, &Dsjj );

				log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dijj, &Dij, &Djj );

				rate = ( G[ij] ) ? exp( logH_ij - logI_p ) : exp( logI_p - logH_ij );
				
				rates[counter++] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
			}
		}
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}
		
		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ )
			{
				K_hat_Cpp[i] += K[i] / sum_rates;
				if( G[i] ) p_links_Cpp[i] += 1.0 / sum_rates;
			}
			
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
    
/*
 * Multiple birth-death MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for maximum a posterior probability estimation (MAP)
*/
void gcgm_DMH_bdmcmc_map_multi_update( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double D[], double Ds[], double *threshold, int *multi_update )
{
	int multi_update_C = *multi_update;
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int counterallG = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	int selected_edge_i, selected_edge_j, selected_edge_ij, counter, size_sample_graph = *size_sample_g;
	bool this_one;

	double Dsijj, Dsjj, Dsij, sum_diag, K022, a11, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, nu_star, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];
	double sum_weights = 0.0;
	long double rate, sum_rates; 
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );        // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2xp2 );         // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings
	// ----------------------------------------

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------|
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|		
		for( j = 1; j < dim; j++ )
		{			
			jj    = j * dim + j;
			Dsjj  = Ds[jj];
			Djj  = D[jj];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ij];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ij];
				Dijj  = - Dij * Dij / Djj;

				log_H_ij( &K[0], &sigma[0], &logH_ij, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dsijj, &Dsij, &Dsjj );

				log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &i, &j,
					   &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
					   &dim, &p1, &p2, &p2xp2, &jj,
					   &Dijj, &Dij, &Djj );

				rate = ( G[ij] ) ? exp( logH_ij - logI_p ) : exp( logI_p - logH_ij );
				
				rates[counter] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
				
				char_g[counter] = G[ij] + '0'; 
				counter++; 
			}
		}
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}

		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) K_hat[i] += K[i] / sum_rates;	

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[counterallG] = 1.0 / sum_rates;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[counterallG];
					all_graphs[counterallG] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[counterallG];
				all_graphs[counterallG]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			counterallG++; 
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy(sample_graphs[i], qp, 0);
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;

	for( i = 0; i < pxp; i++ ) 
		K_hat[i] /= sum_weights;
}
       
/*
 * Reversible Jump MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for Bayesian model averaging
*/

// main function of RJMCMC for GCGMs based on double Metropolis-Hastings
void gcgm_DMH_rjmcmc_ma( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 double K_hat[], int p_links[], 
			 int *b, int *b_star, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	
	int randomEdge, counter, selected_edge_i, selected_edge_j;

	double Dsijj, Dsjj, Dsij, K022, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );        // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2xp2 );         // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings

	double alpha_ij;
	
	GetRNGstate();
	// main loop for Reversible Jump MCMC sampling algorithm ------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------| 
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// STEP 2: selecting edge and calculating alpha -----------------------| 
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		randomEdge = static_cast<int>( runif( 0, 1 ) * qp );

		counter   = 0;
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				if( counter == randomEdge )    
				{
					selected_edge_i = i;
					selected_edge_j = j;
				}
				
				counter++;
			}
		
		// -------- Calculating alpha -----------------------------------------|
		ij    = selected_edge_j * dim + selected_edge_i;
		jj    = selected_edge_j * dim + selected_edge_j;
		Dsij  = Ds[ij];
		Dsjj  = Ds[jj];
		Dsijj = - Dsij * Dsij / Dsjj;
		Dij   = D[ij];
		Djj   = D[jj];
		Dijj  = - Dij * Dij / Djj;
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		log_H_ij( &K[0], &sigma[0], &logH_ij, &selected_edge_i, &selected_edge_j,
               &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
               &dim, &p1, &p2, &p2xp2, &jj,
               &Dsijj, &Dsij, &Dsjj );

		log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &selected_edge_i, &selected_edge_j,
               &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
               &dim, &p1, &p2, &p2xp2, &jj,
               &Dijj, &Dij, &Djj );

		alpha_ij = ( G[ij] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );

		// -------- End calculating alpha -------------------------------------|
		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( runif( 0, 1 ) ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selected_edge_i * dim + selected_edge_j] = G[ij];

			if( G[ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}
		}

		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
			for( i = 0; i < pxp ; i++ )
			{
				K_hat[i]   += K[i];
				p_links[i] += G[i];
			}	
	    // End of saving result -----------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();
}
    
/*
 * Reversible Jump MCMC for Gaussian copula Graphical models  
 * Based on Double Metropolis-Hastings
 * it is for maximum a posterior probability estimation (MAP)
*/
void gcgm_DMH_rjmcmc_map( int *iter, int *burnin, int G[], double Ts[], double Ti[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int counterallG = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	int randomEdge, counter, selected_edge_i, selected_edge_j, size_sample_graph = *size_sample_g;
	bool this_one;

	double Dsijj, Dsjj, Dsij, sum_diag, K022, a11, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, nu_star, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );        // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2xp2 );         // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ---- for copula ------------------------
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// ----------------------------------------
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	double logH_ij, logI_p, Dij, Djj, Dijj;   // for double Metropolis-Hastings
	// ----------------------------------------

	double alpha_ij;
	GetRNGstate();
	// main loop for Reversible Jump MCMC sampling algorithm ------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: copula -----------------------------------------------------| 
		get_Ds( K, Z, R, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
		// STEP 2: selecting edge and calculating alpha
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		randomEdge = static_cast<int>( runif( 0, 1 ) * qp );

		counter   = 0;
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				if( counter == randomEdge )    
				{
					selected_edge_i = i;
					selected_edge_j = j;
				}
				
				char_g[counter++] = G[j * dim + i] + '0'; // adjToString( G, &all_graphs[g], p );
			}
		
		// -------- Calculating alpha -----------------------------------------|
		ij    = selected_edge_j * dim + selected_edge_i;
		jj    = selected_edge_j * dim + selected_edge_j;
		Dsij  = Ds[ij];
		Dsjj  = Ds[jj];
		Dsijj = - Dsij * Dsij / Dsjj;
		Dij   = D[ij];
		Djj   = D[jj];
		Dijj  = - Dij * Dij / Djj;
		
		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		log_H_ij( &K[0], &sigma[0], &logH_ij, &selected_edge_i, &selected_edge_j,
               &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
               &dim, &p1, &p2, &p2xp2, &jj,
               &Dsijj, &Dsij, &Dsjj );

		log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &selected_edge_i, &selected_edge_j,
               &Kj22_inv[0], &Kj12[0], &Kj12xK22_inv[0], &K022, &K12[0], &K22_inv[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma11[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0], &sigma2112[0],
               &dim, &p1, &p2, &p2xp2, &jj,
               &Dijj, &Dij, &Djj );

		alpha_ij = ( G[ij] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );
		// -------- End calculating alpha -------------------------------------|		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( runif( 0, 1 ) ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selected_edge_i * dim + selected_edge_j] = G[ij];

			if( G[ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}
		}

		// STEP 3: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) K_hat[i] += K[i];	

			string_g = string( char_g.begin(), char_g.end() );	
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i]++;     // += all_weights[counterallG];
					all_graphs[counterallG] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[counterallG];
				all_graphs[counterallG]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			counterallG++; 
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy(sample_graphs[i], qp, 0);
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
    
} // End of exturn "C"
