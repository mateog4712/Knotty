#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

void pseudo_loop::compute_WMv_WMp(cand_pos_t i, cand_pos_t j){
	if(j-i+1<4) return;

	WMv.set(i,j) = std::min(E_MLStem(get_energy(i,j),get_energy(i+1,j),get_energy(i,j-1),get_energy(i+1,j-1),i,j),WMv.get(i,j-1) + params_->MLbase);
	WMp.set(i,j) = std::min(P.get(i,j)+PSM_penalty+b_penalty,WMp.get(i,j-1) + params_->MLbase);
}

void pseudo_loop::compute_energy_WM (cand_pos_t i, cand_pos_t j)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    if(j-i+1<4) return;
	energy_t m1 = INF,m2=INF,m3=INF,m4=INF,m5=INF;
	
	for (cand_pos_t k=j-TURN-1; k >= i; --k)
	{
		energy_t wm_kj = E_MLStem(get_energy(k,j),get_energy(k+1,j),get_energy(k,j-1),get_energy(k+1,j-1),k,j);
		energy_t wmb_kj = P.get(k,j)+PSM_penalty+b_penalty;
		m1 = std::min(m1,static_cast<energy_t>((k-i)*params_->MLbase) + wm_kj);
		m2 = std::min(m2,static_cast<energy_t>((k-i)*params_->MLbase) + wmb_kj);
		m3 =  std::min(m3,WM.get(i,k-1) + wm_kj);
		m4 =  std::min(m4,WM.get(i,k-1) + wmb_kj);

	}
	m5 = std::min(m5,WM.get(i,j-1) + params_->MLbase);
	WM.set(i,j) = std::min({m1,m2,m3,m4,m5});  
}

energy_t pseudo_loop::compute_energy_VM(cand_pos_t i, cand_pos_t j)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    energy_t min = INF;
    for (cand_pos_t k = i+1; k <= j-3; ++k)
    {
        energy_t WM2ij = WM.get(i+1,k-1) + WMv.get(k,j-1);
		WM2ij = std::min(WM2ij,WM.get(i+1,k-1) + WMp.get(k,j-1));
		WM2ij = std::min(WM2ij,static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-1));

        energy_t WM2ip1j = WM.get(i+2,k-1) + WMv.get(k,j-1);
		WM2ip1j = std::min(WM2ip1j,WM.get(i+2,k-1) + WMp.get(k-1,j-1));
		WM2ip1j = std::min(WM2ip1j,static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-1));

        energy_t WM2ijm1 = WM.get(i+1,k-1) + WMv.get(k,j-2);
		WM2ijm1 = std::min(WM2ijm1, WM.get(i+1,k-1) + WMp.get(k,j-2));
		WM2ijm1 = std::min(WM2ijm1,static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-2));

        energy_t WM2ip1jm1 = WM.get(i+2,k-1) + WMv.get(k,j-2);
		WM2ip1jm1 = std::min(WM2ip1jm1,WM.get(i+2,k-1) + WMp.get(k,j-2));
		WM2ip1jm1 = std::min(WM2ip1jm1,static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-2));

        min = std::min(min,E_MbLoop(WM2ij,WM2ip1j,WM2ijm1,WM2ip1jm1,i,j));
    }
    return min;
}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t pseudo_loop::HairpinE(const std::string& seq, cand_pos_t i, cand_pos_t j) {
	
	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1_[i+1],S1_[j-1],&seq.c_str()[i-1], const_cast<vrna_param_t*>(params_));
}

/**
 * @brief non-restricted version
*/
energy_t pseudo_loop::compute_internal(cand_pos_t i, cand_pos_t j){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const int ptype_closing = pair[S_[i]][S_[j]];
	for ( cand_pos_t k=i+1; k<=max_k; ++k) {
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
		for (int l=j-1; l>=min_l; --l) {
			energy_t v_iloop_kl = E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<vrna_param_t*>(params_)) + get_energy(k,l);
			v_iloop = std::min(v_iloop,v_iloop_kl);	
		} 
	}
	return v_iloop;
}

void pseudo_loop::compute_energy (cand_pos_t i, cand_pos_t j)
// compute the V(i,j) value, if the structure must be restricted
{
    energy_t min, min_en[3];
    cand_pos_t k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;

	min_en[0] = HairpinE(seq,i,j);
	min_en[1] = compute_internal(i,j);
	min_en[2] = compute_energy_VM(i,j);
    
    for (k=0; k<3; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = INTER; break;
        case  2: type = MULTI; break;
        default: type = NONE;
    }

    if (min < INF/2) {
        int ij = index[i]+j-i;
        V[ij].energy = min;
        V[ij].type = type;
    }
}

void pseudo_loop::compute_PK(const Index4D &x){
	if (impossible_case(x)) return;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	energy_t min_energy = INF;
    int best_branch = 0;
	for(cand_pos_t d=i; d < j; d++){
        energy_t tmp = PK.get(i,d,k,l) + WP.get(d+1,j);  // 12G1

        if (tmp < min_energy){
            min_energy=tmp;
            best_branch = 1;
        }
    }
	for(cand_pos_t d=k+1; d <= l; d++){
        energy_t tmp = PK.get(i,j,d,l) + WP.get(k,d-1);  //1G21
        if (tmp < min_energy){
            min_energy=tmp;
            best_branch = 2;
        }
    }

	energy_t tmp = calc_PX<MType::L>(x) + gamma2(j,i)+PB_penalty;
    if(tmp < min_energy){
        min_energy = tmp;
        best_branch = 3;
    }

    tmp = calc_PX<MType::M>(x) + gamma2(j,k)+PB_penalty;
    if(tmp < min_energy){
        min_energy = tmp;
        best_branch = 4;
    }

    tmp = calc_PX<MType::R>(x) + gamma2(l,k)+PB_penalty;
    if(tmp < min_energy){
        min_energy = tmp;
        best_branch = 5;
    }

    tmp = calc_PX<MType::O>(x) + gamma2(l,i)+PB_penalty;
    if(tmp < min_energy){
        min_energy = tmp;
        best_branch = 6;
    }
	PK.set(x, min_energy);
	if (min_energy < INF/2){
        // adding candidates
        if (best_branch > 1) {
            PK_CL.push_candidate(x, min_energy);
        }
    }
}

/**
 * These should probably just be in a matrix at this point; there are so many arrays, why have four less and recompute every time;
 * It would mean not having four more though, but 11. The space might be more important
 */
energy_t pseudo_loop::calc_PLiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t min_energy = INF, tmp = INF;
	if (i+TURN+2<j) { 
		min_energy = PX.get(i+1,j-1,k,l) + get_e_stP(i,j);
		best_d_=i+1;
        best_dp_=j-1;
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(i,d,dp,j) + PX.get(d,dp,k,l);
			if(tmp<min_energy){
				min_energy = tmp;
				best_d_ = d;
				best_dp_ = dp;
			}
		}
	}
	return min_energy;
}
energy_t pseudo_loop::calc_PMiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t min_energy = INF, tmp = INF;
	if (i<j && k<l) {
		min_energy = PX.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
		best_d_=j-1;
        best_dp_=k+1;
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(d,j,k,dp) + PX.get(i,d,dp,l);
			if(tmp<min_energy){
				min_energy = tmp;
				best_d_ = d;
				best_dp_ = dp;
			}
		}
	}
	return min_energy;
}
energy_t pseudo_loop::calc_PRiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t min_energy = INF, tmp = INF;
	if (k+TURN+2<l) { 
		min_energy = PX.get(i,j,k+1,l-1) + get_e_stP(k,l);
		best_d_=k+1;
        best_dp_=l-1;
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(k,d,dp,l) + PX.get(i,j,d,dp);
			if(tmp<min_energy){
				min_energy = tmp;
				best_d_ = d;
				best_dp_ = dp;
			}
		}
	}
	return min_energy;
}
energy_t pseudo_loop::calc_POiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PX = PX_by_mtype(type);
	energy_t min_energy = INF, tmp = INF;
	if (i<j && k<l ) { 
		min_energy = PX.get(i+1,j,k,l-1) + get_e_stP(i,l);
		best_d_=i+1;
        best_dp_=l-1;
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			tmp = get_intP(i,d,dp,l) + PX.get(d,j,k,dp);
			if(tmp<min_energy){
				min_energy = tmp;
				best_d_ = d;
				best_dp_ = dp;
			}
		}
	}
	return min_energy;
}

void pseudo_loop::compute_WB(cand_pos_t i, cand_pos_t l){
    assert(!impossible_case(i,l));
    WB.set(i,l) = (std::min(cp_penalty*(l-i+1),WBP.get(i,l)));
}
void pseudo_loop::compute_WP(cand_pos_t i, cand_pos_t l){
    assert(!impossible_case(i,l));
    WP.set(i,l) = (std::min(PUP_penalty*(l-i+1),WPP.get(i,l)));
}

energy_t pseudo_loop::Liloop_energy(const Index4D &x, cand_pos_t d, cand_pos_t dp){
    if ( d == x.i()+1 && dp == x.j()-1 ) {
        return get_e_stP(x.i(),x.j());
    } else {
        return get_intP(x.i(),d,dp,x.j());
    }
}

energy_t pseudo_loop::Miloop_energy(const Index4D &x, cand_pos_t d, cand_pos_t dp){
    if ( d == x.j()-1 && dp == x.k()+1 ) {
        return get_e_stP(d,dp);
    } else {
        return get_intP(d,x.j(),x.k(),dp);
    }
}

energy_t pseudo_loop::Riloop_energy(const Index4D &x, cand_pos_t d, cand_pos_t dp){
    if ( d == x.k()+1 && dp == x.l()-1 ) {
        return get_e_stP(x.k(),x.l());
    } else {
        return get_intP(x.k(),d,dp,x.l());
    }
}

energy_t pseudo_loop::Oiloop_energy(const Index4D &x, cand_pos_t d, cand_pos_t dp){
    if ( d == x.i()+1 && dp == x.l()-1 ) {
        return get_e_stP(x.i(),x.l());
    } else {
        return get_intP(x.i(),d,dp,x.l());
    }
}