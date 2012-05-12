#ifndef FSG_H_
#define FSG_H_

#include <vector>
#include <map>
#include <list>
#include <cassert>
#include <algorithm>
#include <iterator>

#include "print_graph.hpp"
#include "lattice.hpp"
#include "fsg_candidate.hpp"

namespace graph_alg
{

    template<class CandidateGraph, class L>
    void fsg_gen_2_edge(GraphCollection<CandidateGraph>& cands,
			const GraphCollection<CandidateGraph>& freqs,
			const L& lops)
    {
	typename GraphCollection<CandidateGraph>::const_iterator i, j;
	for (i = freqs.begin(); i != freqs.end(); ++i)
	{
	    assert(num_vertices(**i) == 2 && num_edges(**i) == 1);
	    for (j = i; j != freqs.end(); ++j)
	    {
		const CandidateGraph& g1 = **i;
		const CandidateGraph& g2 = **j;
	    
		typedef typename
		    boost::graph_traits<CandidateGraph>::edge_descriptor E;
		typedef typename
		    boost::graph_traits<CandidateGraph>::vertex_descriptor V;

		E e1 = *edges(g1).first;
		E e2 = *edges(g2).first;
		V v1_s = source(e1, g1);
		V v1_t = target(e1, g1);
		V v2_s = source(e2, g2);
		V v2_t = target(e2, g2);

		if (lops.label_equal(g1[v1_s], g2[v2_s]))
		    insert_new_candidate<CandidateGraph>(cands, g1, v1_s, g2[v2_t]);
		if (lops.label_equal(g1[v1_s], g2[v2_t]))
		    insert_new_candidate<CandidateGraph>(cands, g1, v1_s, g2[v2_s]);
		if (lops.label_equal(g1[v1_t], g2[v2_s]))
		    insert_new_candidate<CandidateGraph>(cands, g1, v1_t, g2[v2_t]);
		if (lops.label_equal(g1[v1_t], g2[v2_t]))
		    insert_new_candidate<CandidateGraph>(cands, g1, v1_t, g2[v2_s]);
	    }
	}
    }


    template<class CandidateGraph, class L>
    void fsg_join(GraphCollection<CandidateGraph>& cands,
		  const CandidateGraph& g1,
		  const CandidateGraph& g2,
		  const CandidateGraph& ancestor,
		  const typename LatticeLink<CandidateGraph>::M& m1,
		  const typename LatticeLink<CandidateGraph>::M& m2,
		  const L& l)
    {
	typedef typename
	    boost::graph_traits<CandidateGraph>::vertex_descriptor V;
	const V NULV = boost::graph_traits<CandidateGraph>::null_vertex();

	if (m2.e.vertex_ancestor_2 == NULV)
	{
	    V v1 = m1[m2.e.vertex_ancestor_1];
	    insert_new_candidate<CandidateGraph>(cands, g1, v1, g2[m2.e.vertex_descendant_2]);

	    if (l.label_equal(g1[m1.e.vertex_descendant_2], g2[m2.e.vertex_descendant_2]) &&
		v1 != m1.e.vertex_descendant_1)
	    {
		insert_new_candidate<CandidateGraph>(cands, g1, v1, m1.e.vertex_descendant_2);
	    }
	}
	else
	{
	    V v1_s = m1[m2.e.vertex_ancestor_1];
	    V v1_t = m1[m2.e.vertex_ancestor_2];
	    insert_new_candidate<CandidateGraph>(cands, g1, v1_s, v1_t);
	}
    }


    template<class CandidateGraph, class L>
    void fsg_join(GraphCollection<CandidateGraph>& cands,
		  const CandidateGraph& g1,
		  const CandidateGraph& g2,
		  const LatticeLink<CandidateGraph>& lk1,
		  const LatticeLink<CandidateGraph>& lk2,
		  const L& lops)
    {
	typedef typename LatticeLink<CandidateGraph>::MI MI;
	for (MI mi_1 = lk1.mm.begin(); mi_1 != lk1.mm.end(); ++mi_1)
	    for (MI mi_2 = lk2.mm.begin(); mi_2 != lk2.mm.end(); ++mi_2)
		fsg_join(cands, g1, g2, *lk1.ancestor, *mi_1, *mi_2, lops);
    }
    

    template<class CandidateGraph, class L>
    void fsg_gen(GraphCollection<CandidateGraph>& cands,
		 const GraphCollection<CandidateGraph>& freqs,
		 const LatticeLayer<CandidateGraph>& ll,
		 const L& lops)
    {
	typename GraphCollection<CandidateGraph>::const_iterator i, j;
	for (i = freqs.begin(); i != freqs.end(); ++i)
	{
	    const CandidateGraph* g1 = *i;

	    for (j = i; j != freqs.end(); ++j)
	    {
		const CandidateGraph* g2 = *j;

		assert(num_edges(*g1) == num_edges(*g2));

		typedef LatticeLayer<CandidateGraph> LL;
		typedef LatticeLink<CandidateGraph> LK;
		typedef typename LL::const_iterator LLI;
		std::pair<LLI,LLI> lli_g1 = ll.equal_range(g1);
		std::pair<LLI,LLI> lli_g2 = ll.equal_range(g2);
		for ( ; lli_g1.first != lli_g1.second; ++lli_g1.first)
		{
		    const LK& lk1 = lli_g1.first->second;
		    for (LLI llig2 = lli_g2.first; llig2 != lli_g2.second; ++llig2)
		    {
			const LK& lk2 = llig2->second;
			if (lk1.ancestor == lk2.ancestor)
			    fsg_join(cands, *g1, *g2, lk1, lk2, lops);
		    }
		}
	    }
	}
    }


    template<class CandidateGraph, class L>
    void remove_dupls(GraphCollection<CandidateGraph>& cands,
		      const L& lops)
    {
	typename GraphCollection<CandidateGraph>::iterator i, j;
	for (i = cands.begin(); i != cands.end(); ++i)
	{
	    j = i; ++j;
	    while (j != cands.end())
	    {
		VertexCompatible_<CandidateGraph,L> vc(**i, **j, lops);
		if (isomorphism_test(**i, **j, vc))
		{
		    delete *j;
		    j = cands.erase(j);
		}
		else
		    ++j;
	    }
	}
    }


    template<class CandidateGraph, class L>
    void downward_closure(GraphCollection<CandidateGraph>& cands,
			  const GraphCollection<CandidateGraph>& freqs,
			  LatticeLayer<CandidateGraph>& ll,
			  const L& lops)
    {
	//ll.clear();
	typename GraphCollection<CandidateGraph>::iterator ci = cands.begin();
	typename GraphCollection<CandidateGraph>::iterator ci_end = cands.end();
	while (ci != ci_end)
	{  
	    if (! insert(ll, **ci, freqs, lops))
	    {
		delete *ci;
		ci = cands.erase(ci);
	    }
	    else
		++ci;
	}
    }
    


    /////////////////////////////////////////////////////
    //	fsg algorithm
    // -------------------------------------------------
    // template parameters:
    //
    // L must have
    //  label_type
    //  label_type label_value(OV, OriGraph)
    //        bool label_equal(label_type, label_type)
    //
    // label_iterator dereferenced to label_type
    // label_iterator begin()
    // label_iterator end()
    //
    template<class GraphIterator, class L>
    void fsg(Lattice<
	     Candidate<
		 typename std::iterator_traits<GraphIterator>::value_type,
		 L> >& lattice,
	     GraphIterator gi_first, GraphIterator gi_end,
	     const L& lops)
    {
	typedef typename std::iterator_traits<GraphIterator>::value_type OriGraph;
	typedef Candidate<OriGraph, L> CandidateGraph;
	typedef GraphCollection<CandidateGraph> CandidateCollection;

	CandidateCollection cands;

	// construct all frequent 1-edge subgraphs candidates
	typename L::label_iterator li, lj, li_end;
	for (li = lops.begin(), li_end = lops.end(); li != li_end; ++li)
	    for (lj = li; lj != li_end; ++lj)
	{
	    try {
		CandidateGraph* c =
		    new CandidateGraph(gi_first, gi_end, *li, *lj, lops);
		cands.push_back(c);
	    } catch (NotFrequentException) {}
	    
	}
	if (cands.empty())
	    return;
	lattice.push_back(CandidateCollection());
	lattice.back().swap(cands);
	
	// construct all frequent 2-edge subgraphs candidates
	fsg_gen_2_edge(cands, lattice.back(), lops);
	remove_dupls(cands, lops);
	downward_closure(cands, lattice.back(), lattice.ll, lops);

	// construct all other frequent subgraphs candidates
	//int max_i = 3;
	while (!cands.empty())
	{
	    lattice.push_back(CandidateCollection());
	    lattice.back().swap(cands);
	    
	    //if (max_i-- == 0) break;

	    fsg_gen(cands, lattice.back(), lattice.ll, lops);
	    remove_dupls(cands, lops);
	    downward_closure(cands, lattice.back(), lattice.ll, lops);
	}	
    }
}
#endif
