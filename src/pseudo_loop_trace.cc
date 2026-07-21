#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

////////////////////////// Traceback ////////////////////////////////
void pseudo_loop::Trace_W(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("W at %d and %d with %d\n", i, j, e);
	if (j<=i) return;

	energy_t acc = INF;

	// this case is for j unpaired, so I have to check that.
	energy_t tmp = W[j-1];
	if (e==tmp){
		Trace_W(i,j-1,W[j-1]);
		return;
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i>1) ? W[i-1] : 0;
		base_type si1 = i>1 ? S_[i-1] : -1;
		base_type sj1 = j<n ? S_[j+1] : -1;
		tmp = acc + get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_ExtLoop(pair[S_[i]][S_[j]],si1,sj1,params_) : E_ExtLoop(pair[S_[i]][S_[j]],-1,-1,params_));
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			Trace_V(i,j,get_energy(i,j));
			return;
		}
		if(params_->model_details.dangles ==1){
			tmp = acc + get_energy(i+1,j) + E_ExtLoop(pair[S_[i+1]][S_[j]],S_[i],-1,params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j,get_energy(i+1,j));
				return;
			}
			tmp = acc + get_energy(i,j-1) + E_ExtLoop(pair[S_[i]][S_[j-1]],-1,S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i,j-1,get_energy(i,j-1));
				return;
			}
			tmp = acc + get_energy(i+1,j-1) + E_ExtLoop(pair[S_[i+1]][S_[j-1]],S_[i],S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j-1,get_energy(i+1,j-1));
				return;
			}
		}
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i-1>0) ? W[i-1]: 0;
		tmp = acc + P.get(i,j)+ PS_penalty;
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			Trace_P(i,j,P.get(i,j));
			return;
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_V(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("V at %d and %d as type: %c with %d\n", i, j,get_type(i,j), e);
	fres[i] = j;
	fres[j] = i;
	char type = get_type(i,j);

	switch(type){
		case HAIRP:{
			return;
		}
		case INTER:{
			cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
			for (cand_pos_t k = i+1; k <= max_k; ++k){
				cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
				for (cand_pos_t l = j-1; l >= min_l; --l)
				{
					energy_t tmp = compute_int(i,j,k,l) + get_energy(k,l);
					if (e == tmp)
					{
						Trace_V(k,l,get_energy(k,l));
						return;
					}
				}
				
			}
		}
		break;
		case MULTI: {
			energy_t tmp = INF;
			for (cand_pos_t k = i+1; k <= j-1; k++){
				tmp = WM.get(i+1,k-1) + std::min(WMv.get(k,j-1),WMp.get(k,j-1)) + params_->MLclosing;
				if(params_->model_details.dangles == 2){
					tmp += E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_);
				} else {
					tmp += E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_);
				}
				if (e==tmp){
					tmp -= params_->MLclosing;
					tmp -= (params_->model_details.dangles == 2 ? E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) : E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_));
					Trace_WM(i+1,k-1,WM.get(i+1,k-1));
					if(tmp == WM.get(i+1,k-1) + WMv.get(k,j-1)){
						Trace_WMv(k,j-1,WMv.get(k, j-1));
					} else {
						Trace_WMp(k,j-1,WMp.get(k, j-1));
					}
					return;
				}

				tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + WMp.get(k,j-1))+ E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing;
				if (e==tmp){
					Trace_WMp(k,j-1,WMp.get(k,j-1));
					return;
				}
				if(params_->model_details.dangles ==1){
					tmp = WM.get(i+2,k-1) + std::min(WMv.get(k,j-1),WMp.get(k,j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp)
					{
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLbase);
						Trace_WM(i+2,k-1,WM.get(i+2,k-1));
						if(tmp == WM.get(i+2,k-1) + WMv.get(k,j-1)){
							Trace_WMv(k,j-1,WMv.get(k, j-1));
						} else {
							Trace_WMp(k,j-1,WMv.get(k, j-1));
						}
						return;
					}
					tmp = WM.get(i+1,k-1) + std::min(WMv.get(k,j-2),WMp.get(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp)
					{
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLbase);
						Trace_WM(i+1,k-1,WM.get(i+1,k-1));
						if(tmp == WM.get(i+1,k-1) + WMv.get(k,j-2)){
							Trace_WMv(k,j-2,WMv.get(k, j-2));
						} else {
							Trace_WMp(k,j-2,WMv.get(k, j-2));
						}
						return;
					}
					tmp = WM.get(i+2,k-1) + std::min(WMv.get(k,j-2),WMp.get(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
					if (e==tmp)
					{
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + 2*params_->MLbase);
						Trace_WM(i+2,k-1,WM.get(i+2,k-1));
						if(tmp == WM.get(i+2,k-1) + WMv.get(k,j-2)){
							Trace_WMv(k,j-2,WMv.get(k, j-2));
						} else {
							Trace_WMp(k,j-2,WMv.get(k, j-2));
						}
						return;
					}

					if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-1) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-1,WMp.get(k, j-1));
						return;
					}
					tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-2,WMp.get(k, j-2));
						return;
					}
					
					if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-2,WMp.get(k, j-2));
						return;
					}	
				}				
			}
		}
		break;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WM at %d and %d with %d\n", i, j, e);
	energy_t tmp = INF;

	tmp = WM.get(i,j-1)+params_->MLbase;
	if(e==tmp){
		Trace_WM(i,j-1,WM.get(i,j-1));
		return;
	}
	for (cand_pos_t k=i; k <= j-TURN-1; k++){	
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + WMv.get(k,j);
		if(e==tmp){
			Trace_WMv(k,j,WMv.get(k,j));
			return;
		}
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + WMp.get(k,j);
		if(e==tmp){
			Trace_WMp(k,j,WMp.get(k,j));
			return;
		}
		tmp = WM.get(i,k-1) + WMv.get(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,WM.get(i, k-1));
			Trace_WMv(k,j,WMv.get(k,j));
			return;
		}
		tmp = WM.get(i,k-1) + WMp.get(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,WM.get(i, k-1));
			Trace_WMp(k,j,WMp.get(k,j));
			return;
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMv at %d and %d with %d\n", i, j, e);
	cand_pos_t si = S_[i];
	cand_pos_t sj = S_[j];
	cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
	cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
	pair_type tt = pair[S_[i]][S_[j]];
	energy_t tmp = get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_MLstem(tt,si1,sj1,params_) : E_MLstem(tt,-1,-1,params_));
	if(e==tmp){
		Trace_V(i,j,get_energy(i,j));
		return;
	}

	if(params_->model_details.dangles == 1){
		tt = pair[S_[i+1]][S_[j]];
		energy_t tmp = get_energy(i+1,j) + E_MLstem(tt,si,-1,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j,get_energy(i+1,j));
			return;
		}
		tt = pair[S_[i]][S_[j-1]];
		tmp = get_energy(i,j-1) + E_MLstem(tt,-1,sj,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i,j-1,get_energy(i,j-1));
			return;
		}
		tt = pair[S_[i+1]][S_[j-1]];
		tmp = get_energy(i+1,j-1) + E_MLstem(tt,si,sj,params_) + 2*params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j-1,get_energy(i+1,j-1));
			return;
		}
	}

	tmp = WMv.get(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMv(i,j-1,WMv.get(i,j-1));
		return;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMp at %d and %d with %d\n", i, j, e);
	energy_t tmp = P.get(i,j) + PSM_penalty + b_penalty;
	if(e==tmp){
		Trace_P(i,j,WMp.get(i,j));
		return;
	}
	tmp = WMp.get(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMp(i,j-1,WMp.get(i,j-1));
		return;
	}
	UNREACHABLE();
}

void pseudo_loop::Trace_WB(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WB at %d and %d with %d\n", i, l, e);
	if (i>l) return;
	if(e==cp_penalty*(l-i+1)) return;
	if(e==WBP.get(i,l)){
		Trace_WBP(i,l,WBP.get(i,l));
		return;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WBP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WBP at %d and %d with %d\n", i, l, e);
	energy_t tmp = WBP.get(i,l-1)+cp_penalty;
	if(e==tmp){
		Trace_WBP(i,l-1,WBP.get(i,l-1));
		return;
	}
	for(cand_pos_t d=i; d< l; ++d){
		tmp = calc_WB(i,d-1) + get_energy(d,l) + bp_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WB(i,d-1,calc_WB(i,d-1));
			Trace_V(d,l,get_energy(d,l));
			return;
		}
		tmp = calc_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WB(i,d-1,calc_WB(i,d-1));
			Trace_P(d,l,P.get(d,l));
			return;
		}
	}
	UNREACHABLE();
}

// PK portion

void pseudo_loop::Trace_WP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WP at %d and %d with %d\n", i, l, e);
	if (i>l) return;
	if(e==PUP_penalty*(l-i+1)) return;
	if(e==WPP.get(i,l)){
		Trace_WPP(i,l,WPP.get(i,l));
		return;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WPP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WPP at %d and %d with %d\n", i, l, e);
	energy_t tmp = WPP.get(i,l-1)+PUP_penalty;
	if(e==tmp){
		Trace_WPP(i,l-1,WPP.get(i,l-1));
		return;
	}
	for(cand_pos_t d=i; d<l; ++d){
		tmp = calc_WP(i,d-1) + get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		if(e==tmp){
			Trace_WP(i,d-1,calc_WP(i,d-1));
			Trace_V(d,l,get_energy(d,l));
			return;
		}
		tmp = calc_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WP(i,d-1,calc_WP(i,d-1));
			Trace_P(d,l,P.get(d,l));
			return;
		}
	}
	UNREACHABLE();
}

void pseudo_loop::Trace_P(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("P at %d and %d with %d\n", i, l, e);
	energy_t tmp = INF;
    recompute_slice_PK(Index4D(i,l,i,l));

	for (const candidate_PK c : PK_CL[l]) {
        if (c.d() < i) continue; // skip candidates outside current interval

        // decomposition 1212 (where 2 is candidate c)

        int e1 = PK.get(i,c.d()-1,c.j()+1,c.k()-1);
        int e2 = c.w();

        if (e1+e2 == e) {
            Trace_PK(i,c.d()-1,c.j()+1,c.k()-1,e1);
            Trace_PK(c.d(),c.j(),c.k(),l,e2);
            return;
        }
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_PK(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, energy_t e){
    int best_d = -1, best_branch=-1;
    Index4D x = Index4D(i,j,k,l);
    recompute_slice_PK(x);

    energy_t min_energy = INF;

    for(int d=i+1; d < j; d++){
        if (PK.get(i,d,k,l) + WP.get(d+1,j) == e){  // 12G1
            Trace_PK(i,d,k,l,PK.get(i,d,k,l));
            Trace_WP(d+1,j,WP.get(d+1,j));
            return;
        }
    }
    for(int d=k+1; d < l; d++){
        int temp = PK.get(i,j,d,l) + WP.get(k,d-1);
        if (PK.get(i,j,d,l) + WP.get(k,d-1) == e){ // 1G21
            Trace_PK(i,j,d,l,PK.get(i,j,d,l));
            Trace_WP(k,d-1,WP.get(k,d-1));
            return;
        }
    }
    // continue trace with one of the recursion cases to PL,PM,PR,PO
    int best_tgt_energy = INF;
    char best_tgt_type = '\0';
    min_energy = INF;

    for (auto type : {MType::L, MType::M, MType::R, MType::O} ) {
        int px_e = recompute_PX(x, type);
        int pen = penalty(x, gamma2, type) + PB_penalty;
        if (px_e + pen == e) {
            Trace_PX(i,j,k,l,type,px_e + pen);
            return;
            // best_tgt_type = pid_by_mtype(type);
        }
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_PLiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PLiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i+TURN+2<j) { 
		tmp = PX.get(i+1,j-1,k,l) + get_e_stP(i,j);
		if(e==tmp){
			Trace_PX(i+1,j-1,k,l,type,PX.get(i+1,j-1,k,l));
			return;
		}
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(i,d,dp,j) + PX.get(d,dp,k,l);
			if(e==tmp){
				Trace_PX(d,dp,k,l,type,PX.get(d,dp,k,l));
				return;
			}
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_PMiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PMiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i<j && k<l) {
		tmp = PX.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
		if(e==tmp){
			Trace_PX(i,j-1,k+1,l,type,PX.get(i,j-1,k+1,l));
			return;
		}
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(d,j,k,dp) + PX.get(i,d,dp,l);
			if(e==tmp){
				Trace_PX(i,d,dp,l,type,PX.get(i,d,dp,l));
				return;
			}
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_PRiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PRiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (k+TURN+2<l) { 
		tmp = PX.get(i,j,k+1,l-1) + get_e_stP(k,l);
		if(e==tmp){
			Trace_PX(i,j,k+1,l-1,type,PX.get(i,j,k+1,l-1));
			return;
		}
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(k,d,dp,l) + PX.get(i,j,d,dp);
			if(e==tmp){
				Trace_PX(i,j,d,dp,type,PX.get(i,j,d,dp));
				return;
			}
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_POiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "POiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i<j && k<l ) { 
		tmp = PX.get(i+1,j,k,l-1) + get_e_stP(i,l);
		if(e==tmp){
			Trace_PX(i+1,j,k,l-1,type,PX.get(i+1,j,k,l-1));
			return;
		}
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(i,d,dp,l) + PX.get(d,j,k,dp);
			if(e==tmp){
				Trace_PX(d,j,k,dp,type,PX.get(d,j,k,dp));
				return;
			}
		}
	}
	UNREACHABLE();
}
// /**
//  * 
//  * 
//  * 
//  * 
//  */
// void pseudo_loop::Trace_PfromL(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e){
// 	if (debug) std::cout << "PfromL at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
// 	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 	for(cand_pos_t d=i; d<=j; ++d){
// 		energy_t tmp = calc_WP(i,d-1) + PfromXprime.get(d,j,k,l);
// 		if(e==tmp){
// 			Trace_WP(i,d-1,calc_WP(i,d-1));
// 			Trace_PfromXprime(d,j,k,l,type,PfromXprime.get(d,j,k,l));
// 			return;
// 		}
// 	}
// 	UNREACHABLE();
// }
// void pseudo_loop::Trace_PfromM(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e){
// 	if (debug) std::cout << "PfromM at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
// 	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 	for(cand_pos_t d=i; d<=j; ++d){
// 		energy_t tmp = PfromXprime.get(i,d,k,l) + calc_WP(d+1,j);
// 		if(e==tmp){
// 			Trace_PfromXprime(i,d,k,l,type,PfromXprime.get(i,d,k,l));
// 			Trace_WP(d+1,j,calc_WP(d+1,j));
// 			return;
// 		}
// 	}
// 	UNREACHABLE();
// }
// void pseudo_loop::Trace_PfromR(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e){
// 	if (debug) std::cout << "PfromR at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
// 	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 	for(cand_pos_t d=k; d<=l; ++d){
// 		energy_t tmp = calc_WP(k,d-1) + PfromXprime.get(i,j,d,l);
// 		if(e==tmp){
// 			Trace_WP(k,d-1,calc_WP(k,d-1));
// 			Trace_PfromXprime(i,j,d,l,type,PfromXprime.get(i,j,d,l));
// 			return;
// 		}
// 	}
// 	UNREACHABLE();
// }
// void pseudo_loop::Trace_PfromO(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e){
// 	if (debug) std::cout << "PfromO at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
// 	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 	for(cand_pos_t d=i; d<=j; ++d){
// 		energy_t tmp = calc_WP(i,d-1) + PfromXprime.get(d,j,k,l);
// 		if(e==tmp){
// 			Trace_WP(i,d-1,calc_WP(i,d-1));
// 			Trace_PfromXprime(d,j,k,l,type,PfromXprime.get(d,j,k,l));
// 			return;
// 		}
// 	}
// 	UNREACHABLE();
// }
// /**
//  * 
//  * 
//  * 
//  */
void pseudo_loop::Trace_PLmloop1(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PLmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i+1; d <= j; ++d){
        tmp = PX.get(i,d,k,l) + calc_WB(d+1,j);
        if(e==tmp){
			Trace_PX(i,d,k,l,type,PX.get(i,d,k,l));
			Trace_WB(d+1,j,calc_WB(d+1,j));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_PMmloop1(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PMmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = k+1; d <= l; ++d){
        tmp = PX.get(i,j,d,l) + calc_WB(k,d-1);
        if(e==tmp){
			Trace_PX(i,d,k,l,type,PX.get(i,d,k,l));
			Trace_WB(d+1,j,calc_WB(d+1,j));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_PRmloop1(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PRmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = k+1; d <= l; ++d){
        tmp = PX.get(i,j,k,d) + calc_WB(d+1,l);
        if(e==tmp){
			Trace_PX(i,j,k,d,type,PX.get(i,j,k,d));
			Trace_WB(d+1,l,calc_WB(d+1,l));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_POmloop1(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "POmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d=k+1; d<=l;++d){
        tmp=PX.get(i,j,k,d) + calc_WB(d+1,l);
        if(e==tmp){
			Trace_PX(i,j,k,d,type,PX.get(i,j,k,d));
			Trace_WB(d+1,l,calc_WB(d+1,l));
			return;
		}
    }
	UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void pseudo_loop::Trace_PLmloop0(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PLmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp=INF;
	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        tmp = cp_penalty*(d-i) + PX.get(d,j,k,l);
		if(e==tmp){
			Trace_PX(d,j,k,l,type,PX.get(d,j,k,l));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_PMmloop0(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PMmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        tmp = cp_penalty*(j-d) + PX.get(i,d,k,l);
        if(e==tmp){
			Trace_PX(i,d,k,l,type,PX.get(i,d,k,l));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_PRmloop0(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PRmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = k; d < l; d++){
        tmp = cp_penalty*(d-k) + PX.get(i,j,d,l);
		if(e==tmp){
			Trace_PX(i,j,d,l,type,PX.get(i,j,d,l));
			return;
		}
    }
	UNREACHABLE();
}
void pseudo_loop::Trace_POmloop0(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "POmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	energy_t tmp = INF;
	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        tmp = cp_penalty*(d-i) + PX.get(d,j,k,l);
		if(e==tmp){
			Trace_PX(d,j,k,l,type,PX.get(d,j,k,l));
			return;
		}
    }
	UNREACHABLE();
}