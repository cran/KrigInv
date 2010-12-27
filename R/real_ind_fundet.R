real_ind_fundet <- function(x,T,epsilon){
	real_ind_fundet <- ( fundet(x) > T-epsilon )&( fundet(x) < T+epsilon )
	real_ind_fundet <- as.numeric(real_ind_fundet)
}

