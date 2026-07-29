#include "candidate_list.hh"
#include <iostream>

candidate_lists::list_t candidate_lists::empty_list = candidate_lists::list_t();

/**
 *  @brief push candidate with information w, i, to front of CL
 */
void candidate_lists::push_candidate(const Index4D &x, int w) {
    assert( w >= std::numeric_limits<energy_16t>::min() && w <= std::numeric_limits<energy_16t>::max() );
    assert( x.i() >= std::numeric_limits<index_t>::min() && x.i() <= std::numeric_limits<index_t>::max() );

    cls_[x.j()][index(x.k(),x.l())].push_sorted(x.i(), w);
}

/**
 * Find candidate in candidate list CL
 * @returns on failure returns nullptr, else candidate
 */
int candidate_lists::find_candidate(int i, int j, int k, int l) const {
    const auto it = cls_[j].find(index(k,l));
    if (it == cls_[j].end()) {
        return INF;
    }
    const auto it2 = it->second.find(i);
    if (it2 == it->second.end()) {
        return INF;
    }
    return it2->second;
}

/** @brief adds number of candidates or empty candidate lists to candidates/empty_lists
*   used in print_CL_sizes
*/
void candidate_lists::get_CL_size(int &candidates, int &capacity) const {
    for (const auto &x : cls_ ){
        for (const auto &cl : x ){
            candidates += cl.second.size();
            capacity += cl.second.capacity();
        }
    }
}

/**
 * @brief prints information on a single candidate list
 */
void candidate_lists::print_CL_size(std::string type) const {
    int candidates = 0, capacity = 0;
    get_CL_size(candidates, capacity);

    std::cout << std::endl;
    std::cout << type << std::endl;


    int num_lists=0;
    for(const auto &x:cls_)
        num_lists+=x.size();

    printf("Num lists: %d\n",num_lists);

    printf("Num candidates: %d\n", candidates);

    printf("Avg cands per list: %f\n", (float)candidates/num_lists);


    int c_size = sizeof(candidate);
    printf("Total candidate space: %d\n", candidates*c_size );

    printf("Size: %d, Capacity: %d\n", candidates, capacity);
}
