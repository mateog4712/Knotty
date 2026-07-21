#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "constants.hh"
#include "matrices.hh"
#include "candidate_list.hh"
#include "trace_arrow.hh"
#include <cstdio>
#include <cstring>

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"
#define debug 0

#ifdef NDEBUG
	#define UNREACHABLE() __builtin_unreachable()
#else
	#define UNREACHABLE() \
		do { \
			std::cerr << "Reached unreachable at line " << __LINE__ << " in File: " << __FILE__ << std::endl; \
			abort(); \
		} while(0)
#endif


class pseudo_loop{

public:
	// constructor
	pseudo_loop(std::string seq, int dangle);

	// destructor
	~pseudo_loop();

	double ccj ();
    void compute_energies(cand_pos_t i, cand_pos_t j);
	
	char get_type (cand_pos_t i, cand_pos_t j) { cand_pos_t ij = index[i]+j-i; return V[ij].type;}
    energy_t get_energy (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return V[ij].energy; }

	std::string structure;
private:

 	// function to allocate space for the arrays
    void allocate_space();

	cand_pos_t n;
	std::string res;
	std::string seq;

	std::vector<int> fres;
	vrna_param_t *params_;

	std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WPP and WBP
	index_offset_t3 index3D;

	short *S_;
	short *S1_;
	std::vector<energy_t> W; // size n+1 so left as non-Trianglematrix
	std::vector<free_energy_node> V;
	TriangleMatrix WM;
	TriangleMatrix WMv;
	TriangleMatrix WMp;
	TriangleMatrix WP;
	TriangleMatrix WB;
	TriangleMatrix WPP;	// similar to WP but has at least one base pair
	TriangleMatrix WBP;	// similar to WB but has at least one base pair
	
	TriangleMatrix P;					// the main loop for pseudoloops and bands
	MatrixSlices3D PK;				// MFE of a TGB structure over gapped region [i,j] U [k,l]

	MatrixSlices3D PL;		// MFE of a TGB structure s.t. i.j is paired
	MatrixSlices3D PLmloop0;
	MatrixSlices3D PLmloop1;
	MatrixSlices3D PfromL;

	MatrixSlices3D PR;		// MFE of a TGB structure s.t. k.l is paired
	MatrixSlices3D PRmloop0;
	MatrixSlices3D PRmloop1;
	MatrixSlices3D PfromR;

	MatrixSlices3D PM;		// MFE of a TGB structure s.t. j.k is paired
	MatrixSlices3D PfromM;
	MatrixSlices3D PMmloop0;
	MatrixSlices3D PMmloop1;

	MatrixSlices3D PO;	// MFE of a TGB structure s.t. i.l is paired multiple
	MatrixSlices3D POmloop0;
	MatrixSlices3D POmloop1;
	MatrixSlices3D PfromO;

	// Candidates and Trace arrows
	MasterTraceArrows *ta;

	candidate_lists *PLmloop0_CL;
    candidate_lists *POmloop0_CL;
    candidate_lists *PfromL_CL;
    candidate_lists *PfromO_CL;
    candidate_lists *PMmloop0_CL;
    candidate_lists *PRmloop0_CL;
    candidate_lists *PfromM_CL;
    candidate_lists *PfromR_CL;

    CandidateListsPK PK_CL;

	void compute_WBP(cand_pos_t i, cand_pos_t l);
	void compute_WPP(cand_pos_t i, cand_pos_t l);
	void compute_WP(cand_pos_t i, cand_pos_t l);
	void compute_WB(cand_pos_t i, cand_pos_t l);
		
	void compute_P(cand_pos_t i, cand_pos_t l);
	void compute_PK(const Index4D &x);
	void compute_PX(const Index4D &x, MType type);
	void compute_PfromX(const Index4D &x, MType type);
	void compute_PXmloop0(const Index4D &x, MType type);
	void compute_PXmloop1(const Index4D &x, MType type);
	//Get rid of these later
	template<MType type> int calc_PX_checked(const Index4D &x){
		assert(x.difference(type) > TURN);

		return PX_by_mtype(type).get(x);
	}

	template <MType type> int calc_PX(const Index4D &x){
		const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];
		if (!(ptype_closing>0)) {
			return INF;
		}
		return calc_PX_checked<type>(x);
	}

	energy_t calc_PLiloop(const Index4D &x, MType type);
	energy_t calc_PRiloop(const Index4D &x, MType type);
	energy_t calc_PMiloop(const Index4D &x, MType type);
	energy_t calc_POiloop(const Index4D &x, MType type);

	energy_t calc_WB(cand_pos_t i, cand_pos_t l);
	energy_t calc_WP(cand_pos_t i, cand_pos_t l);


	// Traceback //
	// void backtrack();
	// void Trace_W(cand_pos_t i, cand_pos_t j, energy_t e);
	// void Trace_P(cand_pos_t i, cand_pos_t l, energy_t e);
	// void Trace_V(cand_pos_t i, cand_pos_t j, energy_t e);
	// void Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e);
	// void Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e);
	// void Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e);

	// void Trace_WB(cand_pos_t i, cand_pos_t l, energy_t e);
	// void Trace_WBP(cand_pos_t i, cand_pos_t l, energy_t e);
	// void Trace_WP(cand_pos_t i, cand_pos_t l, energy_t e);
	// void Trace_WPP(cand_pos_t i, cand_pos_t l, energy_t e);

	// void Trace_PX(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);
	// void Trace_PXiloop(const Index4D &x, MType type, energy_t e);
	// void Trace_PXmloop(const Index4D &x, MType type, energy_t e);
	// void Trace_PXmloop0(const Index4D &x, MType type, energy_t e);
	// void Trace_PXmloop1(const Index4D &x, MType type, energy_t e);
	// void Trace_PfromX(const Index4D &x, MType type, energy_t e);

	// void Trace_PLiloop(const Index4D &x, MType type, energy_t e);
	// void Trace_PMiloop(const Index4D &x, MType type, energy_t e);
	// void Trace_PRiloop(const Index4D &x, MType type, energy_t e);
	// void Trace_POiloop(const Index4D &x, MType type, energy_t e);

	// void Trace_PLmloop0(const Index4D &x, MType type, energy_t e);
	// void Trace_PMmloop0(const Index4D &x, MType type, energy_t e);
	// void Trace_PRmloop0(const Index4D &x, MType type, energy_t e);
	// void Trace_POmloop0(const Index4D &x, MType type, energy_t e);

	// void Trace_PLmloop1(const Index4D &x, MType type, energy_t e);
	// void Trace_PMmloop1(const Index4D &x, MType type, energy_t e);
	// void Trace_PRmloop1(const Index4D &x, MType type, energy_t e);
	// void Trace_POmloop1(const Index4D &x, MType type, energy_t e);

	// void Trace_PfromL(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	// void Trace_PfromM(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	// void Trace_PfromR(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	// void Trace_PfromO(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);

	void compute_energy (cand_pos_t i, cand_pos_t j);
	energy_t HairpinE(const std::string& seq, cand_pos_t i, cand_pos_t j);
	energy_t compute_internal(cand_pos_t i, cand_pos_t j);
	void compute_energy_WM (cand_pos_t i, cand_pos_t j);
	energy_t compute_energy_VM (cand_pos_t i, cand_pos_t j);
	void compute_WMv_WMp(cand_pos_t i, cand_pos_t j);
	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j);
	energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j);
	energy_t E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, cand_pos_t i, cand_pos_t j);

	// penalty for closing pair i.l or l.i of a pseudoloop
	static constexpr energy_t gamma2(cand_pos_t i, cand_pos_t l){
		return 0;
	}
	template<class Penalty> energy_t penalty(const Index4D &x, Penalty p, MType type) {
		switch(type) {
		case MType::L: return p(x.j(),x.i());
		case MType::M: return p(x.j(),x.k());
		case MType::R: return p(x.l(),x.k());
		case MType::O: return p(x.l(),x.i());
		}
		UNREACHABLE();
	}
	energy_t calc_PXiloop(const Index4D &x, MType type);
	energy_t calc_PXmloop(const Index4D &x, MType type);
	energy_t generic_decomposition(int i, int j, int k, int l, int decomp_cases, candidate_lists &CL, const TriangleMatrix &w, const MatrixSlices3D &PX, int LMRO_cases = 0, energy_t penalty = 0);
	inline bool impossible_case(const Index4D &x) const {
    	return !x.is_valid(n);
	}
	inline bool impossible_case(int i, int l) const {
		return((i<=0 ||l<=0 || i>n || l >n || i> l));
	}
	// output parameters of generic_decomposition function
    // and recompute PL functions
    int best_d_; //!< best split point (end of the gap matrix or first interior loop end)
    int best_dp_; //!< best second split point (interior loop)
    int best_branch_; //!< index of best branch
    bool decomposing_branch_; //!< whether best branch is decomposing
    int best_tgt_type_;
    int best_tgt_energy_;
	// cases for the generic decomposition
    static const int CASE_12G2 = 1<<0;
    static const int CASE_12G1 = 1<<1;
    static const int CASE_1G21 = 1<<2;
    static const int CASE_1G12 = 1<<3;
    static const int CASE_L = CASE_12G2 | CASE_12G1;
    static const int CASE_M = CASE_12G1 | CASE_1G21;
    static const int CASE_R = CASE_1G21 | CASE_1G12;
    static const int CASE_O = CASE_12G2 | CASE_1G12;
    static const int CASE_PL = 1<<4;
    static const int CASE_PM = 1<<5;
    static const int CASE_PR = 1<<6;
    static const int CASE_PO = 1<<7;

	int lmro_case(MType type) const {
        return select_by_mtype<int,CASE_L,CASE_M,CASE_R,CASE_O>(type);
    }
	int lmro_caseP(MType type) const {
        return select_by_mtype<int,CASE_PL,CASE_PM,CASE_PR,CASE_PO>(type);
    }
	//! non-decomposing cases in the from recursions by type
    //! @param type
    int lmro_cases_in_fromX_by_mtype(MType type) const {
        return select_by_mtype<int,
                            CASE_PM | CASE_PR | CASE_PO, // fromL
                            CASE_PL | CASE_PR,           // fromM
                            CASE_PM | CASE_PO,           // fromR
                            CASE_PL | CASE_PR>           // fromO
            (type);
    }

	inline MatrixSlices3D& PX_by_mtype(MType type) {
		static std::array<MatrixSlices3D*,4> matrices{&PL, &PM, &PR, &PO};
		return *matrices[static_cast<int>(type)];
	}
	inline MatrixSlices3D& PfromX_by_mtype(MType type) {
		static std::array<MatrixSlices3D*,4> matrices{&PfromL, &PfromM, &PfromR, &PfromO};
		return *matrices[static_cast<int>(type)];
	}
	inline MatrixSlices3D& PXmloop0_by_mtype(MType type) {
		static std::array<MatrixSlices3D*,4> matrices{&PLmloop0, &PMmloop0, &PRmloop0, &POmloop0};
		return *matrices[static_cast<int>(type)];
	}
	inline MatrixSlices3D& PXmloop1_by_mtype(MType type) {
		static std::array<MatrixSlices3D*,4> matrices{&PLmloop1, &PMmloop1, &PRmloop1, &POmloop1};
		return *matrices[static_cast<int>(type)];
	}
	inline candidate_lists& PXmloop0_CL_by_mtype(MType type) {
		static std::array<candidate_lists*,4> matrices{PLmloop0_CL, PMmloop0_CL, PRmloop0_CL, POmloop0_CL};
		return *matrices[static_cast<int>(type)];
	}
	inline candidate_lists& PfromX_CL_by_mtype(MType type) {
		static std::array<candidate_lists*,4> matrices{PfromL_CL, PfromM_CL, PfromR_CL, PfromO_CL};
		return *matrices[static_cast<int>(type)];
	}
};
#endif /*PSEUDO_LOOP_H_*/
