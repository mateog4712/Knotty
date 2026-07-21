#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

pseudo_loop::pseudo_loop(std::string seq, int dangle) : seq(seq), params_(vrna_params(NULL))
{
	n = seq.length();
	params_->model_details.dangles = dangle;
	make_pair_matrix();
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
    allocate_space();
}

void pseudo_loop::allocate_space()
{

	TriangleMatrix::new_index(index,n+1);
	MatrixSlices3D::construct_index(index3D,n);

	ta = new MasterTraceArrows(n, index);

    PfromL_CL = new candidate_lists(MType::L, n);
    PfromM_CL = new candidate_lists(MType::M, n);
    PfromR_CL = new candidate_lists(MType::R, n);
    PfromO_CL = new candidate_lists(MType::O, n);
    PLmloop0_CL = new candidate_lists(MType::L, n);
	PMmloop0_CL = new candidate_lists(MType::M, n);
    PRmloop0_CL = new candidate_lists(MType::R, n);
    POmloop0_CL = new candidate_lists(MType::O, n);

    PK_CL.resize(n);

	fres.resize(n+1,-2);
	structure = std::string (n+1,'.');

	W.resize(n+1,0);
	cand_pos_t total_length = ((n+1) *(n+2))/2;
	V.resize(total_length);
	WM.init(n+1,index);
	WMv.init(n+1,index);
	WMp.init(n+1,index);
    WBP.init(n+1,index);
	WPP.init(n+1,index);
	WB.init(n+1,index);
	WP.init(n+1,index);
	P.init(n+1,index);
	
	// 4D matrix initialization
	PK.init(n,index3D);
	PL.init(n,index3D);
	PfromL.init(n,index3D);
	PLmloop0.init(n,index3D);
	PLmloop1.init(n,index3D);

	PR.init(n,index3D);
	PfromR.init(n,index3D);
	PRmloop0.init(n,index3D);
	PRmloop1.init(n,index3D);

	PM.init(n,index3D);
	PfromM.init(n,index3D);
	PMmloop0.init(n,index3D);
	PMmloop1.init(n,index3D);

	PO.init(n,index3D);
	PfromO.init(n,index3D);
	POmloop0.init(n,index3D);
	POmloop1.init(n,index3D);
}

pseudo_loop::~pseudo_loop()
{
	free(params_);
	free(S_);
	free(S1_);
	delete ta;
	delete PfromL_CL;
	delete PfromM_CL; // switch these to non-pointers
	delete PfromR_CL;
	delete PfromO_CL;
	delete PLmloop0_CL;
	delete PMmloop0_CL;
	delete PRmloop0_CL;
	delete POmloop0_CL;
}

double pseudo_loop::ccj (){
	for (cand_pos_t i = n; i>=1; --i){	
		for (cand_pos_t j =i; j<=n; ++j){
			compute_energy (i,j);
			compute_energies(i,j);
			compute_WMv_WMp(i,j);
			compute_energy_WM(i,j);
		}
	}
	for (cand_pos_t j= TURN+1; j <= n; j++){
		energy_t m1 = INF, m2 = INF, m3 = INF;
		m1 = W[j-1];
		for (cand_pos_t k=1; k<=j-TURN-1; ++k){
			energy_t acc = (k>1) ? W[k-1]: 0;
			m2 = std::min(m2,acc + E_ext_Stem(get_energy(k,j),get_energy(k+1,j),get_energy(k,j-1),get_energy(k+1,j-1),k,j));
			m3 = std::min(m3,acc + P.get(k,j) + PS_penalty);
		}
		W[j] = std::min({m1,m2,m3});
	}

    double energy = W[n]/100.0;

	// backtrack
	// backtrack();

	fill_structure(fres,structure);
	this->structure = structure.substr(1,n);
    return energy;
}

void pseudo_loop::compute_energies(cand_pos_t i, cand_pos_t l)
{

	// 1) compute all energies over region [i,l]
	compute_P(i,l);
	compute_WBP(i,l);
	compute_WPP(i,l);
	compute_WB(i,l);
	compute_WP(i,l);

	//2) compute all energies over gapped region [i,j]U[k,l]
	for(cand_pos_t j = i; j<l; ++j){
		// Hosna, July 8, 2014
		// in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
		for(cand_pos_t k = l; k>=j+2; --k){
			Index4D x(i,j,k,l);

			compute_PXmloop0(x,MType::L);
			compute_PXmloop0(x,MType::M);
			compute_PXmloop0(x,MType::R);
			compute_PXmloop0(x,MType::O);

			compute_PXmloop1(x,MType::L);
			compute_PXmloop1(x,MType::M);
			compute_PXmloop1(x,MType::R);
			compute_PXmloop1(x,MType::O);

			compute_PX(x,MType::L);
			compute_PX(x,MType::M);
			compute_PX(x,MType::R);
			compute_PX(x,MType::O);

			compute_PfromX(x,MType::L);
			compute_PfromX(x,MType::M);
			compute_PfromX(x,MType::R);
			compute_PfromX(x,MType::O);

			compute_PK(x);
		}
	}
}

void pseudo_loop::compute_WBP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy= INF, b1 = INF, b2=INF, b3 = INF,tmp =INF;

	for(cand_pos_t d=i; d< l; ++d){
		tmp = calc_WB(i,d-1) + get_energy(d,l) + bp_penalty + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = calc_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WBP.get(i,l-1)+cp_penalty;
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		WBP.set(i,l) = min_energy;
	}
}

void pseudo_loop::compute_WPP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF, b1 = INF, b2=INF, b3 =INF, tmp = INF;

	for(cand_pos_t d=i; d<l; ++d){
		tmp = calc_WP(i,d-1) + get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = calc_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WPP.get(i,l-1)+PUP_penalty;
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		WPP.set(i,l) = min_energy;
	}
}

void pseudo_loop::compute_P(cand_pos_t i, cand_pos_t l){
    int min_energy = INF, temp=INF;

    if (impossible_case(i,l) || i==l) {return;}

    for (const candidate_PK c : PK_CL[l]) {
        int w = c.w();
        Index4D x(i,c.d(),c.j(),c.k());
        temp = PK.get(x + Index4D(0, -1, 1, -1)) + w;

        if (temp < min_energy) {
            min_energy = temp;
        }
    }

    // SW: no trace arrows required, since PK can be recomputed from PK candidates
    // (note: in the case of P / PK this is slightly overoptimized, since it
    // saves only O(n^2) TAs)

    if (min_energy < INF/2){
        P.set(i,l) = min_energy;
    }
}
/**
 * This was a case where the code was almost the same for all recurrences.
 * Each time, you would take the respective recurrence, shrink the related indices
 * by one and then add the penalties. This could be solve by the Index4D shrink function
 * which shrinks based on the Mtype. And so 10 functions become one
 */
energy_t pseudo_loop::calc_PXmloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;

	Index4D xp(x);
	xp.shrink(type);
	MatrixSlices3D &PXmloop00 = PXmloop1_by_mtype(type);
	return PXmloop00.get(xp.i(),xp.j(),xp.k(),xp.l())+ ap_penalty + bp_penalty;
}

energy_t pseudo_loop::calc_PXiloop(const Index4D &x, MType type){
	switch(type) {
    case MType::L: return calc_PLiloop(x,type);
    case MType::M: return calc_PMiloop(x,type);
    case MType::R: return calc_PRiloop(x,type);
    case MType::O: return calc_POiloop(x,type);
    }
    UNREACHABLE();
}
void pseudo_loop::compute_PfromX(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PfromX = PfromX_by_mtype(type);
	candidate_lists &PfromX_CL = PfromX_CL_by_mtype(type);
	energy_t min_energy = generic_decomposition(i, j, k, l, lmro_case(type), PfromX_CL, WP, PfromX, lmro_cases_in_fromX_by_mtype(type), penalty(x,gamma2,type) + PB_penalty);

	PfromX.set(i,j,k,l,min_energy);
	if (min_energy < INF/2){
        if (!decomposing_branch_ && i < j) {
            PfromX_CL.push_candidate(x, min_energy);
        }
    }
}
void pseudo_loop::compute_PXmloop1(const Index4D &x, MType type){
	if (impossible_case(x)) return;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	MatrixSlices3D &PXmloop0 = PXmloop0_by_mtype(type);
	candidate_lists &PXmloop0_CL = PXmloop0_CL_by_mtype(type);
	energy_t min_energy = generic_decomposition(i, j, k, l, lmro_case(type), PXmloop0_CL, WBP, PXmloop0);

	MatrixSlices3D &PXmloop1 = PXmloop1_by_mtype(type);
	PXmloop1.set(i,j,k,l,min_energy);
}
void pseudo_loop::compute_PXmloop0(const Index4D &x, MType type){
	if (impossible_case(x)) return;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	
	MatrixSlices3D &PXmloop0 = PXmloop0_by_mtype(type);
	candidate_lists &PXmloop0_CL = PXmloop0_CL_by_mtype(type);
	energy_t min_energy = generic_decomposition(i, j, k, l, lmro_case(type), PXmloop0_CL, WB, PXmloop0, lmro_caseP(type), bp_penalty);
	
	PXmloop0.set(i,j,k,l,min_energy);

	if (min_energy < INF/2){
        if ( !decomposing_branch_ ){
            PXmloop0_CL.push_candidate(x, min_energy);
        }
    }
}
/**
 * As this just calls the other functions, we can reduce it to just a PX.
 * ptype closing can be done because lend and rend give the correct indices
 * for each recurrence. The functions coming after just use their respective
 * PX functions
 */
void pseudo_loop::compute_PX(const Index4D &x, MType type){
    energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;
	MatrixSlices3D &PX = PX_by_mtype(type);

	const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];

	if (ptype_closing>0){
		b1 = calc_PXiloop(x, type);
		b2 = calc_PXmloop(x,type);

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (x.difference(type)>TURN){
			if(type == MType::O){
				if(x.i()==x.j() && x.k()==x.l()){
					b3=gamma2(x.l(),x.i());
				}
			}
			Index4D xp(x);
			xp.shrink(type);
			MatrixSlices3D &PfromX = PfromX_by_mtype(type);
			b3 = PfromX.get(xp) + penalty(xp, gamma2, type);
		}
	}
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		PX.setI(x, min_energy);
	}
}
///////////////// Traceback ////////////////////////////////

// void pseudo_loop::backtrack(){
//    Trace_W(1,n,W[n]);
//    return;
// }

// void pseudo_loop::Trace_PXmloop(const Index4D &x, MType type, energy_t e){
// 	if (debug) std::cout << "PXmloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
// 	assert(impossible_case(x));

// 	Index4D xp(x);
// 	xp.shrink(type);
// 	MatrixSlices3D &PXmloop1 = PXmloop1_by_mtype(type);
// 	energy_t tmp = PXmloop1.get(xp.i(),xp.j(),xp.k(),xp.l())+ ap_penalty + bp_penalty;
// 	if(e==tmp){
// 		return Trace_PXmloop1(xp,type,PXmloop1.get(xp.i(),xp.j(),xp.k(),xp.l()));
// 	}
// 	UNREACHABLE();
// }

// void pseudo_loop::Trace_PXiloop(const Index4D &x, MType type, energy_t e){
// 	if (debug) std::cout << "PXiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
// 	switch(type) {
//     case MType::L: return Trace_PLiloop(x,type,e);
//     case MType::M: return Trace_PMiloop(x,type,e);
//     case MType::R: return Trace_PRiloop(x,type,e);
//     case MType::O: return Trace_POiloop(x,type,e);
//     }
//     UNREACHABLE();
// }

// void pseudo_loop::Trace_PfromX(const Index4D &x, MType type,energy_t e){
//     switch(type) {
//     case MType::L: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
//     case MType::M: return Trace_PfromM(x.i(),x.j(),x.k(),x.l(),type,e);
//     case MType::R: return Trace_PfromR(x.i(),x.j(),x.k(),x.l(),type,e);
//     case MType::O: return Trace_PfromO(x.i(),x.j(),x.k(),x.l(),type,e);
//     }
//     UNREACHABLE();
// }

// void pseudo_loop::Trace_PXmloop1(const Index4D &x, MType type, energy_t e){
// 	if (debug) std::cout << "PXmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
//     switch(type) {
//     case MType::L: return Trace_PLmloop1(x,type,e);
//     case MType::M: return Trace_PMmloop1(x,type,e);
//     case MType::R: return Trace_PRmloop1(x,type,e);
//     case MType::O: return Trace_POmloop1(x,type,e);
//     }
//     UNREACHABLE();
// }
// void pseudo_loop::Trace_PXmloop0(const Index4D &x, MType type, energy_t e){
// 	if (debug) std::cout << "PXmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
//     switch(type) {
//     case MType::L: return Trace_PLmloop0(x,type,e);
//     case MType::M: return Trace_PMmloop0(x,type,e);
//     case MType::R: return Trace_PRmloop0(x,type,e);
//     case MType::O: return Trace_POmloop0(x,type,e);
//     }
//     UNREACHABLE();
// }
// void pseudo_loop::Trace_PX(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
// 	if (debug) std::cout << "PX at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and en: " << e << std::endl;
// 	const Index4D x(i,j,k,l);
// 	const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];
// 	fres[x.lend(type)] = x.rend(type);
// 	fres[x.rend(type)] = x.lend(type);

// 	if (ptype_closing>0){
// 		energy_t tmp = calc_PXmloop(x,type);
// 		if(e==tmp){
// 			Trace_PXmloop(x,type,tmp);
// 			return;
// 		}
// 		if (x.difference(type)>TURN){
// 			if(type == MType::O){
// 				if(x.i()==x.j() && x.k()==x.l()){
// 					if(e==gamma2(x.l(),x.i())) return;
// 				}
// 			}
// 			Index4D xp(x);
// 			xp.shrink(type);
// 			MatrixSlices3D &PfromX = PfromX_by_mtype(type);
// 			tmp = PfromX.get(xp) + penalty(xp, gamma2, type);
// 			if(e==tmp){
// 				Trace_PfromX(xp,type,PfromX.get(xp));
// 				return;
// 			}
// 		}
// 		tmp = calc_PXiloop(x, type);
// 		if(e==tmp){
// 			Trace_PXiloop(x,type,tmp);
// 			return;
// 		}
// 	}
// 	UNREACHABLE();
// }

energy_t pseudo_loop::generic_decomposition(int i, int j, int k, int l, int decomp_cases, candidate_lists &CL, const TriangleMatrix &w, const MatrixSlices3D &PX, int LMRO_ndcases, energy_t penfun) {
    auto x = Index4D(i,j,k,l);

    int min_energy = INF;

    best_branch_ = -1;
    best_d_ = -1;
    decomposing_branch_ = 0;
    best_tgt_energy_ = INF;

    if ( decomp_cases & CASE_12G2 ) {
		// Ian Wark Jan 23 2017
		// 12G2 using candidate list
		for (const auto& c : CL.get_list(j, k, l)) {
			if ( c.first <= i ) break;
			int temp = w.get_uc(i, c.first - 1) + c.second;
			if (temp < min_energy) {
				min_energy = temp;
				best_branch_ = CASE_12G2;
				best_tgt_energy_ = c.second;
				best_d_ = c.first;
			}
		}
    }

    if (decomp_cases & CASE_12G1) {
        // case 1G21
        for (int d = i; d < j; d++) {
            int px_e = PX.get_uc(i, d, k, l);
            int temp = px_e + w.get_uc(d + 1, j);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_12G1;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G21) {
        // case 1G21
        for(int d = k+1; d <= l; d++){
            int px_e = PX.get_uc(i, j, d, l);
            int temp = w.get_uc(k, d - 1) + px_e;
            if (temp < min_energy){
                min_energy = temp;
                best_branch_ = CASE_1G21;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (decomp_cases & CASE_1G12) {
        // case 1G12
        for (int d = k; d < l; d++) {
            int px_e = PX.get_uc(i, j, k, d);
            int temp = px_e + w.get_uc(d + 1, l);
            if (temp < min_energy) {
                min_energy = temp;
                best_branch_ = CASE_1G12;
                best_tgt_energy_ = px_e;
                best_d_ = d;
            }
        }
    }

    if (LMRO_ndcases & CASE_PL) {
        // int px_e = calc_PX<MType::L>(x);
        int px_e = calc_PX<MType::L>(x);
        int temp = px_e + penfun;
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PL;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PM) {
        int px_e = calc_PX<MType::M>(x);
        int temp = px_e + penfun;
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PM;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PR) {
        int px_e = calc_PX<MType::R>(x);
        int temp = px_e + penfun;
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PR;
            best_tgt_energy_ = px_e;
        }
    }

    if (LMRO_ndcases & CASE_PO) {
        int px_e = calc_PX<MType::O>(x);
        int temp = px_e + penfun;
        if (temp < min_energy) {
            min_energy = temp;
            best_branch_ = CASE_PO;
            best_tgt_energy_ = px_e;
        }
    }

    decomposing_branch_ = best_branch_ & ( CASE_12G2 | CASE_12G1 | CASE_1G21 | CASE_1G12 );

    return min_energy;
}

////////////////// Util Functions ///////////////////////////

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t pseudo_loop::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type tt, base_type s5, base_type s3) {
        if (!valid || v == INF) return;
        e = std::min(e, v + E_ExtLoop(tt, s5, s3, params_));
    };
	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, true, pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1);
	if (dangle1) {
        consider(vi1j, j-i-1 > TURN, pair[S_[i+1]][S_[j]], si, -1);
        consider(vij1, j-1-i > TURN, pair[S_[i]][S_[j-1]], -1, sj);
        consider(vi1j1, j-1-i-1 > TURN, pair[S_[i+1]][S_[j-1]], si, sj);
    }
	return e;
}

/**
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t pseudo_loop::E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type type, base_type s5, base_type s3, int ml_count) {
        if (!valid || v == INF) return;
        e = std::min(e, v + E_MLstem(type, s5, s3, params_) + ml_count * params_->MLbase);
    };

	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, true, pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1, 0);
	if (dangle1) {
		consider(vi1j,  j-i-1 > TURN, pair[S_[i+1]][S_[j]], si, -1, 1);
        consider(vij1,  j-1-i > TURN, pair[S_[i]][S_[j-1]], -1, sj, 1);
        consider(vi1j1, j-1-i-1 > TURN, pair[S_[i+1]][S_[j-1]], si, sj, 2);
	}
    return e;
}

/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
*/
energy_t pseudo_loop::E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, cand_pos_t i, cand_pos_t j){
	energy_t e = INF;

    pair_type tt = pair[S_[j]][S_[i]];
    base_type si1 = S_[i+1];
    base_type sj1 = S_[j-1];

	auto consider = [&](energy_t v, base_type s5, base_type s3, int ml_count) {
        if (v == INF) return;
        e = std::min(e, v + E_MLstem(tt, s5, s3, params_) + params_->MLclosing + ml_count * params_->MLbase);
    };

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(WM2ij, dangle2 ? sj1 : -1, dangle2 ? si1 : -1, 0);
	if(dangle1){
		// ML pair 5 — closing (i,j) with mb part [i+2, j-1]
		consider(WM2ip1j, -1, si1, 1);
        // ML pair 3 — closing (i,j) with mb part [i+1, j-2]
        consider(WM2ijm1, sj1, -1, 1);
        // ML pair 53 — closing (i,j) with mb part [i+2, j-2]
        consider(WM2ip1jm1, sj1, si1, 2);
	}
	return e;
}

energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<vrna_param_t*>(params_));
}

energy_t pseudo_loop::get_e_stP(cand_pos_t i, cand_pos_t j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	energy_t ss = compute_int(i,j,i+1,j-1);
	return lrint(e_stP_penalty * ss);
}

energy_t pseudo_loop::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j){
	energy_t e_int = compute_int(i,j,ip,jp);
	energy_t energy = lrint(e_intP_penalty * e_int);
	return energy;
}