#ifndef FSG_H_
#define FSG_H_

#include <vector>
#include <map>
#include <list>
#include <cassert>
#include <algorithm>
#include <iterator>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "vf2.hpp"

#define BR asm("int $3;");

namespace graph_alg
{
    namespace detail
    {
	////////////////////////////////////////////////////
	//	class ParentChildMapping
	// -------------------------------------------------
	template<class Candidate>
	struct ParentChildMapping
	{
	    typedef typename Candidate::MV V;
	    std::vector<V> p2c;
	    V operator[] (V v) const { return p2c.at(v); }
	    
	    struct ExtEdge
	    {
		V vertex_descendant_1;
		V vertex_descendant_2;
		V vertex_ancestor_1;
		V vertex_ancestor_2; // may be null_vertex()
	    } e;
	};


	////////////////////////////////////////////////////
	//	class ParentLink
	// -------------------------------------------------
	template<class Candidate>
	struct ParentLink
	{
	    const Candidate* parent;

	    typedef typename std::list<ParentChildMapping<Candidate> > PCMS;
	    PCMS pcms;

	    void swap(ParentLink& pl)
		{
		    std::swap(parent, pl.parent);
		    pcms.swap(pl.pcms);
		}
	};


	////////////////////////////////////////////////////
	//	class Candidate
	// -------------------------------------------------
	class NotFrequentException {};

	template<class G, class L>
	class Candidate :
	    public boost::adjacency_list<boost::setS,
					 boost::vecS,
					 boost::undirectedS,
					 typename L::label_type>
	{
	public:
	    typedef boost::adjacency_list<boost::setS,
					  boost::vecS,
					  boost::undirectedS,
					  typename L::label_type> ModelGraph;

	    typedef typename boost::graph_traits<ModelGraph>::vertex_descriptor MV;
	    typedef typename boost::graph_traits<ModelGraph>::edge_descriptor ME;
	    
	    typedef G SrcGraph;
	    typedef typename boost::graph_traits<SrcGraph>::vertex_descriptor SV;
	    typedef typename boost::graph_traits<SrcGraph>::edge_descriptor SE;
	    typedef typename L::label_type LabelType;
	 
	    // ctor for one edge graph
	    template<class SrcGraphIterator>
	    Candidate(SrcGraphIterator gi_begin, SrcGraphIterator gi_end,
		      const LabelType& lab1, const LabelType& lab2,
		      const L& l);

	    // create new edge and new vertex
	    Candidate(const Candidate& parent, MV mv, const LabelType& lab);

	    // create new edge
	    Candidate(const Candidate& parent, MV mv1, MV mv2);

	    

	    //private:
	    const L& lops;

	    // ------------------------------
	    // links to SrcGraph`s subgraphs
	    typedef std::vector<SV> M2S;
	    typedef std::list<M2S> M2SList;
	    typedef std::list<std::pair<const SrcGraph*, M2SList> > SBS;
	    SBS sbs;

	    // ------------------------------
	    // links to parent candidates
	    // forms lattice structure
	    std::list<ParentLink<Candidate> > parents;
	    typedef typename std::list<ParentLink<Candidate> >::const_iterator P_CI;


	    bool label_eq(MV mv, SV sv, const SrcGraph& sg) const
		{ return lops.label_equal((*this)[mv], lops.label_value(sv, sg)); }

	    static void m2s_construct(M2SList& result, const M2S& m2s,
				      const SrcGraph& sg, const ModelGraph& mg,
				      MV mv1, MV mv2, const L& l);

	    void sbs_construct(const SBS& par_sbs, MV mv1, MV mv2);

	    int id;
	    static int id_;
	};


	template<class G, class L> int Candidate<G,L>::id_ = 0;

	template<class G, class L>
	void Candidate<G,L>::m2s_construct(M2SList& result, const M2S& m2s,
					   const SrcGraph& sg, const ModelGraph& mg,
					   MV mv1, MV mv2, const L& l)
	{
	    assert(mv1 < m2s.size());
	    SV sv1 = m2s[mv1];

	    typedef typename boost::graph_traits<G>::adjacency_iterator AdjI;
	    typedef std::pair<AdjI,AdjI> AP;
	    if (mv2 < m2s.size())
	    {
		// mv2 has already mapped to subgraph of the SrcGraph
		SV sv2 = m2s[mv2];
		for (AP ap = adjacent_vertices(sv1, sg);
		     ap.first != ap.second; ++ap.first)
		{
		    if (*ap.first == sv2)
		    {
			result.push_back(m2s);
			break;
		    }
		}
	    }
	    else
	    {
		// mv2 has not mapped to subgraph of the SrcGraph
		const LabelType& lab2 = mg[mv2];
		for (AP ap = adjacent_vertices(sv1, sg);
		     ap.first != ap.second; ++ap.first)
		{
		    SV sv2 = *ap.first;
		    if (std::find(m2s.begin(), m2s.end(), sv2) != m2s.end())
			continue;

		    if (l.label_equal(lab2, l.label_value(sv2, sg)))
		    {
			// map mv2 --> sv2
			result.push_back(m2s);
			result.back().push_back(sv2);
			assert(result.back().size() - 1 == mv2);
		    }
		}
	    }
	}

	

	template<class G, class L>
	void Candidate<G,L>::sbs_construct(const SBS& par_sbs, MV mv1, MV mv2)
	{
	    for (typename SBS::const_iterator i = par_sbs.begin(); i != par_sbs.end(); ++i)
	    {
		M2SList sblist_new;
		const SrcGraph& sg = *i->first;
		const M2SList& sblist_par = i->second;
		for (typename M2SList::const_iterator j = sblist_par.begin();
		     j != sblist_par.end(); ++j)
		{
		    m2s_construct(sblist_new, *j, sg, *this, mv1, mv2, lops);
		}

		if (!sblist_new.empty())
		{
		    sbs.push_back(std::make_pair(&sg, M2SList()));
		    sbs.back().second.swap(sblist_new);
		}
	    }

	    if (sbs.empty())
		throw NotFrequentException();
	}


	template<class G, class L>
	template<class SrcGraphIterator>
	Candidate<G,L>::Candidate(SrcGraphIterator gi_begin,
				  SrcGraphIterator gi_end,
				  const LabelType& lab1, const LabelType& lab2,
				  const L& l)
	    :lops(l)
	{
	    id = ++id_;

	    MV mv1 = add_vertex(*this);
	    MV mv2 = add_vertex(*this);
	    ME me = add_edge(mv1, mv2, *this).first;
	    (*this)[mv1] = lab1;
	    (*this)[mv2] = lab2;

	    // for each SrcGraph
	    for (; gi_begin != gi_end; ++gi_begin)
	    {
		M2SList sblist_new;
		const SrcGraph& sg = *gi_begin;

		// for each edges
		typedef typename boost::graph_traits<SrcGraph>::edge_iterator SEI;
		for (std::pair<SEI,SEI> ep = edges(sg);
		     ep.first != ep.second;
		     ++ep.first)
		{
		    SV sv1 = source(*ep.first, sg);
		    SV sv2 = target(*ep.first, sg);

		    if (label_eq(mv1, sv1, sg) && label_eq(mv2, sv2, sg))
		    {
			sblist_new.push_back(M2S());
			sblist_new.back().push_back(sv1);
			sblist_new.back().push_back(sv2);
		    }

		    if (label_eq(mv1, sv2, sg) && label_eq(mv2, sv1, sg))
		    {
			sblist_new.push_back(M2S());
			sblist_new.back().push_back(sv2);
			sblist_new.back().push_back(sv1);
		    }
		}

		if (!sblist_new.empty())
		{
		    sbs.push_back(std::make_pair(&sg, M2SList()));
		    sbs.back().second.swap(sblist_new);
		}
	    }

	    if (sbs.empty())
		throw NotFrequentException();
	}


	template<class G, class L>
	Candidate<G,L>::Candidate(const Candidate& parent,
				  MV mv, const LabelType& lab)
	    :ModelGraph(parent), lops(parent.lops)
	{
	    id = ++id_;

	    MV mv_new = add_vertex(*this);
	    std::pair<ME, bool> r = add_edge(mv, mv_new, *this);
	    (*this)[mv_new] = lab;

	    sbs_construct(parent.sbs, mv, mv_new);
	}


	template<class G, class L>
	Candidate<G,L>::Candidate(const Candidate& parent, MV mv1, MV mv2)
	    :ModelGraph(parent), lops(parent.lops)
	{
	    id = ++id_;

	    std::pair<ME, bool> r = add_edge(mv1, mv2, *this);
	    if (!r.second)
		throw NotFrequentException();

	    sbs_construct(parent.sbs, mv1, mv2);	    
	}


	////////////////////////////////////////////////////
	//		class VertexCompatible
	// -------------------------------------------------
	// adapter for vf2 functions
	//
	template<class G, class L>
	struct VertexCompatible
	{
	    VertexCompatible(const G& g1,
			     const G& g2,
			     const L& lops)
		:g1(g1),g2(g2),lops(lops) {}
	    const G& g1;
	    const G& g2;
	    const L& lops;
	    typedef typename boost::graph_traits<G>::vertex_descriptor V;
	    bool operator() (V vg1, V vg2) const
		{ return lops.label_equal(g1[vg1], g2[vg2]); }
	};


	////////////////////////////////////////////////////
	//		class ConnectedSubgraph
	// -------------------------------------------------
	class NotConnectedException {};

	template<class G>
	class ConnectedSubgraph
	{
	public:
	    typedef typename boost::graph_traits<G>::vertex_descriptor V;
	    typedef typename boost::graph_traits<G>::edge_descriptor E;

	    template <typename T>
	    struct Predicate
	    {
		Predicate() :enabled(false), v() {}
		Predicate(const T& x) :enabled(true), v(x) {}
		bool operator() (const T& x) const { return enabled ? v != x : true; }
		bool enabled;
		T v;
	    };
	
	    int num_ver;

	    typedef boost::filtered_graph<G, Predicate<E>, Predicate<V> > FilGraph;

	    ConnectedSubgraph(const ConnectedSubgraph&);
	    ConnectedSubgraph& operator= (const ConnectedSubgraph&);

	    // ************************
	    // Color Map
	    // ************************
	    typedef char ColorValue;
	    typedef boost::color_traits<ColorValue> Color;
	    typedef std::vector<ColorValue> ColorMap;
	    static void put(ColorMap& cm, V v, ColorValue c) { cm[v] = c; }
	    static ColorValue get(ColorMap& cm, V v) { return cm[v]; }
	    static bool search_other_path(const G& g, V v, E e, ColorMap& color);
	public:
	    FilGraph* fil;
	    V v1, v2;


	    ConnectedSubgraph(const G& graph) 
		:num_ver(boost::num_vertices(graph)),
		 v1(boost::graph_traits<G>::null_vertex()),
		 v2(boost::graph_traits<G>::null_vertex())
		{ fil = new FilGraph(graph,Predicate<E>()); }
	
	    ConnectedSubgraph(const G& graph, const E& removed_edge)	    
		:num_ver(boost::num_vertices(graph)),
		 fil(0),
		 v1(boost::graph_traits<G>::null_vertex()),
		 v2(boost::graph_traits<G>::null_vertex())
		{
		    ColorMap color(num_ver, Color::white());
		
		    if (degree(v2 = source(removed_edge, graph), graph) == 1)
		    {
			v1 = target(removed_edge, graph);
			--num_ver;
			fil = new FilGraph(graph,
					   Predicate<E>(removed_edge),
					   Predicate<V>(v2));
		    }
		    else if (degree(v2 = target(removed_edge, graph), graph) == 1)
		    {
			v1 = source(removed_edge, graph);
			--num_ver;
			fil = new FilGraph(graph,
					   Predicate<E>(removed_edge),
					   Predicate<V>(v2));
		    }
		    else if (search_other_path(graph, source(removed_edge, graph),
					       removed_edge, color))
		    {
			v1 = source(removed_edge, graph);
			v2 = target(removed_edge, graph);
			fil = new FilGraph(graph, Predicate<E>(removed_edge));
		    }
		    else
			throw NotConnectedException();		    
		}

	    ~ConnectedSubgraph() { delete fil; }
	
	    operator const FilGraph& () const { return *fil; }
	    int num_vertices() const { return num_ver; }
	};


	template<class G>
	bool ConnectedSubgraph<G>::search_other_path(const G& g, V v, E edge, ColorMap& color)
	{
	    ConnectedSubgraph<G>::put(color, v, Color::gray());
	    if (v == target(edge, g))
		return true;
	    typedef typename boost::graph_traits<G>::out_edge_iterator OutI;
	    for (std::pair<OutI,OutI> ap = out_edges(v, g);
		 ap.first != ap.second; ++ap.first)
	    {
		if (*ap.first == edge)
		    continue;
		V u = target(*ap.first, g);
		ColorValue c = get(color, u);
		if (c == Color::white())
		{
		    if (search_other_path(g, u, edge, color))
			return true;
		}
	    }
	    put(color, v, Color::black());
	    return false;
	}


	template<class G, class L, class TwoMPairSec>
	bool connsb_isomorphism(const ConnectedSubgraph<G>& g1,
				const ConnectedSubgraph<G>& g2,
				const L& lops,
				TwoMPairSec& sec2)
	{
	    typedef typename ConnectedSubgraph<G>::FilGraph FilGraph;
	    VertexCompatible<FilGraph, L> vc(g1, g2, lops);
	
	    return isomorphism_all(*g1.fil, *g2.fil,
				   sec2,
				   vc,
				   g1.num_vertices(),
				   g2.num_vertices());
	}


    } // namespace detail


    ////////////////////////////////////////////////////////
    //		class FSG_Lattice
    // -----------------------------------------------------
    //
    // template parameters:
    //
    // G transactional graph
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
    template<class G, class L>
    class FSG_Lattice
    {
    public:
	template<class GI>
	FSG_Lattice(GI gi_begin, GI gi_end, const L& l = L());
	
	int get_num_graphs() const { assert(num_graphs >= 0); return num_graphs; }
	int get_num_all_created_graphs() const { return num_all_created_graphs; }
	//private:
	typedef G SrcGraph;
	
	typedef detail::Candidate<SrcGraph, L> CandidateGraph;
	typedef typename boost::graph_traits<CandidateGraph>::edge_descriptor E;
	typedef typename boost::graph_traits<CandidateGraph>::vertex_descriptor V;
	typedef typename boost::graph_traits<CandidateGraph>::edge_iterator EI;

	typedef boost::ptr_vector<CandidateGraph> CandidateCollection;
	typedef typename CandidateCollection::const_iterator CG_CI;
	typedef typename CandidateCollection::iterator CG_I;

	typedef detail::ParentLink<CandidateGraph> PL;
	
	const L& l;
	std::list<CandidateCollection> lattice;
	
	typedef typename
	std::list<CandidateCollection>::const_iterator CC_CI;
	
	void join(CandidateCollection& cands,
		  const CandidateGraph& g1, const CandidateGraph& g2,
		  const CandidateGraph& parent,
		  const detail::ParentChildMapping<CandidateGraph>& m1,
		  const detail::ParentChildMapping<CandidateGraph>& m2);
	
	void gen(CandidateCollection& cands,
		 const CandidateCollection& freqs);
	
	void gen_2_edge(CandidateCollection& cands,
			const CandidateCollection& freqs);

	bool downward_closure_test(CandidateGraph& cand,
				   const CandidateCollection& freqs);

	void downward_closure(CandidateCollection& cands,
			      const CandidateCollection& freqs);

	void remove_dupls(CandidateCollection& cands);

	void insert_new_candidate(CandidateCollection& cands,
				  const CandidateGraph& ancestor,
				  V v,
				  const typename L::label_type& lab)
	    {
		try {
		    cands.push_back(new CandidateGraph(ancestor, v, lab));
		    ++num_graphs;
		    ++num_all_created_graphs;
		} catch (detail::NotFrequentException) {}
	    }

	void insert_new_candidate(CandidateCollection& cands,
				  const CandidateGraph& ancestor,
				  V v, V u)
	    {
		try {
		    cands.push_back(new CandidateGraph(ancestor, v, u));
		    ++num_graphs;
		    ++num_all_created_graphs;
		} catch (detail::NotFrequentException) {}
	    }

	CG_I del_candidate(CG_I i, CandidateCollection& coll)
	    { --num_graphs; return coll.erase(i); }

	int num_graphs;
	int num_all_created_graphs;
    };


    template<class G, class L>
    void FSG_Lattice<G,L>::join(CandidateCollection& cands,
				const CandidateGraph& g1, const CandidateGraph& g2,
				const CandidateGraph& parent,
				const detail::ParentChildMapping<CandidateGraph>& m1,
				const detail::ParentChildMapping<CandidateGraph>& m2)
    {
	const V NULV = boost::graph_traits<CandidateGraph>::null_vertex();
	
	if (m2.e.vertex_ancestor_2 == NULV)
	{
	    V v1 = m1[m2.e.vertex_ancestor_1];
	    insert_new_candidate(cands, g1, v1, g2[m2.e.vertex_descendant_2]);
	    
	    if (l.label_equal(g1[m1.e.vertex_descendant_2],
			      g2[m2.e.vertex_descendant_2]) &&
		v1 != m1.e.vertex_descendant_1)
	    {
		insert_new_candidate(cands, g1, v1, m1.e.vertex_descendant_2);
	    }
	}
	else
	{
	    V v1_s = m1[m2.e.vertex_ancestor_1];
	    V v1_t = m1[m2.e.vertex_ancestor_2];
	    insert_new_candidate(cands, g1, v1_s, v1_t);
	}
    }


    template<class G, class L>
    void FSG_Lattice<G,L>::gen(CandidateCollection& cands,
			       const CandidateCollection& freqs)
    {
	for (CG_CI i = freqs.begin(); i != freqs.end(); ++i)
	    for (CG_CI j = i; j != freqs.end(); ++j)
	    {
		const CandidateGraph& g1 = *i;
		const CandidateGraph& g2 = *j;

		assert(num_edges(g1) == num_edges(g2));
		
		typedef typename CandidateGraph::P_CI P_CI;
		for (P_CI pi = g1.parents.begin(); pi != g1.parents.end(); ++pi)
		    for (P_CI pj = g2.parents.begin(); pj != g2.parents.end(); ++pj)
		    {
			if (pi->parent == pj->parent)
			{
			    typedef typename PL::PCMS::const_iterator PCM_CI;
			    for (PCM_CI pcm_i = pi->pcms.begin();
				 pcm_i != pi->pcms.end(); ++pcm_i)
				for (PCM_CI pcm_j = pj->pcms.begin();
				     pcm_j != pj->pcms.end(); ++pcm_j)
				{
				    join(cands, g1, g2, *pi->parent, *pcm_i, *pcm_j);
				}
			}
		    }
	    }
    }


    template<class G, class L>
    void FSG_Lattice<G,L>::gen_2_edge(CandidateCollection& cands,
				      const CandidateCollection& freqs)
    {
	for (CG_CI i = freqs.begin(); i != freqs.end(); ++i)
	    for (CG_CI j = i; j != freqs.end(); ++j)
	    {
		const CandidateGraph& g1 = *i;
		const CandidateGraph& g2 = *j;
		
		E e1 = *edges(g1).first;
		E e2 = *edges(g2).first;
		V v1_s = source(e1, g1);
		V v1_t = target(e1, g1);
		V v2_s = source(e2, g2);
		V v2_t = target(e2, g2);

		if (l.label_equal(g1[v1_s], g2[v2_s]))
		    insert_new_candidate(cands, g1, v1_s, g2[v2_t]);
		if (l.label_equal(g1[v1_s], g2[v2_t]))
		    insert_new_candidate(cands, g1, v1_s, g2[v2_s]);
		if (l.label_equal(g1[v1_t], g2[v2_s]))
		    insert_new_candidate(cands, g1, v1_t, g2[v2_t]);
		if (l.label_equal(g1[v1_t], g2[v2_t]))
		    insert_new_candidate(cands, g1, v1_t, g2[v2_s]);
	    }
    }


    template<class G, class L>
    bool FSG_Lattice<G,L>::downward_closure_test(CandidateGraph& cand,
						 const CandidateCollection& freqs)
    {
	bool ret = false;
	using namespace std;
	typedef detail::ConnectedSubgraph<CandidateGraph> ConnSubGraph;

	// for each parent graph
	for (CG_CI i = freqs.begin(); i != freqs.end(); ++i)
	{
	    PL pl;
	    pl.parent = &*i;

	    ConnSubGraph par_graph(*i);
	    
	    // masking each edges
	    for (pair<EI,EI> ep = edges(cand); ep.first != ep.second; ++ep.first)
	    {
		E edge = *ep.first;
		
		try
		{
		    ConnSubGraph cand_subgraph(cand, edge);
		    
		    // find all isomorphisms
		    typedef std::vector<MPair> Sec;
		    typedef std::vector<Sec> Sec2;
		    Sec2 sec2;
		    connsb_isomorphism(par_graph, cand_subgraph, l, sec2);
		    typedef Sec2::const_iterator II;
		    for (II ii = sec2.begin(); ii != sec2.end(); ++ii)
		    {
			typedef detail::ParentChildMapping<CandidateGraph> PCM;
			pl.pcms.push_back(PCM());
			PCM& pcm = pl.pcms.back();

			const V NULV = boost::graph_traits<G>::null_vertex();
			
			pcm.e.vertex_descendant_1 = source(edge, cand);
			pcm.e.vertex_descendant_2 = target(edge, cand);
			pcm.e.vertex_ancestor_1 = NULV;
			pcm.e.vertex_ancestor_2 = NULV;
			pcm.p2c.resize(num_vertices(*i), NULV);

			for (Sec::const_iterator j = ii->begin();
			     j != ii->end(); ++j)
			{
			    pcm.p2c[j->first] = j->second;
			    
			    if (pcm.e.vertex_descendant_1 == static_cast<V>(j->second))
				pcm.e.vertex_ancestor_1 = j->first;
			    if (pcm.e.vertex_descendant_2 == static_cast<V>(j->second))
				pcm.e.vertex_ancestor_2 = j->first;
			}

			if (pcm.e.vertex_ancestor_1 == NULV)
			{
			    std::swap(pcm.e.vertex_descendant_1, pcm.e.vertex_descendant_2);
			    std::swap(pcm.e.vertex_ancestor_1, pcm.e.vertex_ancestor_2);
			}

			assert(pcm.e.vertex_descendant_1 != NULV);
			assert(pcm.e.vertex_descendant_2 != NULV);
			assert(pcm.e.vertex_ancestor_1 != NULV);
			assert(find(pcm.p2c.begin(), pcm.p2c.end(), NULV) == pcm.p2c.end());
		    }

		}
		catch (detail::NotConnectedException) {}
	    } // for: masking each edges

	    if (! pl.pcms.empty())
	    {
		cand.parents.push_back(PL());
		cand.parents.back().swap(pl);
		ret = true;
	    }
	}

	return ret;
    }


    template<class G, class L>
    void FSG_Lattice<G,L>::downward_closure(CandidateCollection& cands,
					    const CandidateCollection& freqs)
    {
	CG_I i = cands.begin();
	CG_I i_end = cands.end();
	while (i != i_end)
	{
	    if (! downward_closure_test(*i, freqs))
		i = del_candidate(i, cands);
	    else
		++i;
	}
    }


    template<class G, class L>
    void FSG_Lattice<G,L>::remove_dupls(CandidateCollection& cands)
    {
	for (CG_I i = cands.begin(); i != cands.end(); ++i)
	{
	    CG_I j = i; ++j;
	    while (j != cands.end())
	    {
		detail::VertexCompatible<CandidateGraph,L> vc(*i, *j, l);
		if (isomorphism_test(*i, *j, vc))
		   j = del_candidate(j, cands);
		else
		    ++j;
	    }
	}
    }


    template<class G, class L>
    template<class GI>   
    FSG_Lattice<G,L>::FSG_Lattice(GI gi_begin, GI gi_end, const L& l_)
	: l(l_), num_graphs(0), num_all_created_graphs(0)
    {
	CandidateCollection cands;

	// construct all frequent 1-edge subgraphs candidates
	typename L::label_iterator li, lj, li_end;
	for (li = l.begin(), li_end = l.end(); li != li_end; ++li)
	    for (lj = li; lj != li_end; ++lj)
	    {
		try
		{
		    cands.push_back(new CandidateGraph(gi_begin, gi_end,
						       *li, *lj, l));
		    ++num_graphs;
		    ++num_all_created_graphs;
		}
		catch (detail::NotFrequentException) {}
	    }
	if (cands.empty())
	    return;
	lattice.push_back(CandidateCollection());
	lattice.back().swap(cands);

	// construct all frequent 2-edge subgraphs candidates
	gen_2_edge(cands, lattice.back());
	remove_dupls(cands);
	downward_closure(cands, lattice.back());

	// construct all other frequent subgraphs candidates
	while (! cands.empty())
	{
	    lattice.push_back(CandidateCollection());
	    lattice.back().swap(cands);

	    gen(cands, lattice.back());
	    remove_dupls(cands);
	    //BR;
	    downward_closure(cands, lattice.back());
	    //BR;
	}
    }

}
#endif
