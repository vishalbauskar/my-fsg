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

#include "vf2.hpp"

namespace graph_alg
{
    template<class G>
    struct gtraits
    {
	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef typename boost::graph_traits<G>::out_edge_iterator OutEI;
	typedef std::pair<OutEI,OutEI> OutEIP;

	typedef typename boost::graph_traits<G>::adjacency_iterator AdjI;
	typedef std::pair<AdjI,AdjI> AdjIP;

	typedef typename boost::graph_traits<G>::edge_iterator EI;
	typedef std::pair<EI,EI> EIP;

	static V vnil() { return boost::graph_traits<G>::null_vertex(); }
    };


    // adapter for vf2 functions
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

	
    /////////////////////////////////////////////////////
    //	class Candidate
    // -------------------------------------------------
    // template parameters:
    //
    // L must have
    //  label_type
    //  label_type label_value(OV, OriGraph)
    //        bool label_equal(label_type, label_type)
    //
    class NotFrequentException {};

    template<class G, class L>
    class Candidate :
	public boost::adjacency_list<boost::setS,
				     boost::vecS,
				     boost::undirectedS,
				     typename L::label_type>
    {
	static int N_G;
    public:
	int ID;

	// -----------------------------------------
	typedef G OriGraph;
	    
	typedef typename L::label_type VertexLabel;

	typedef typename gtraits<OriGraph>::V OV;
	typedef typename gtraits<OriGraph>::E OE;

	typedef boost::adjacency_list<boost::setS,
				      boost::vecS,
				      boost::undirectedS,
				      VertexLabel> ModelGraph;
	typedef typename gtraits<ModelGraph>::V MV;
	typedef typename gtraits<ModelGraph>::E ME;
	
	typedef std::vector<MV> M2O;
	typedef std::list<M2O> M2OList;
	typedef std::list<std::pair<const OriGraph*, M2OList> > OGS;
	
	const L& lops;
	OGS ogs;

	static void m2o_construct(M2OList& result,
				  const M2O& m2o,
				  const OriGraph& og,
				  const ModelGraph& mg,
				  MV mv1, MV mv2,
				  const L& lops);

	void mk(const OGS& r_ogs, MV mv1, MV mv2);

	// -----------------------------------------
	template<class OriGraphIterator>
	Candidate(OriGraphIterator, OriGraphIterator,
		  const VertexLabel&, const VertexLabel&,
		  const L& lops = L());
	
	Candidate(const Candidate&, MV, const VertexLabel&);
	
	Candidate(const Candidate&, MV, MV);
    };

    template<class G, class L>
    int Candidate<G,L>::N_G = 0;

    template<class G, class L>
    void Candidate<G,L>::m2o_construct(M2OList& result,
				       const M2O& m2o,
				       const OriGraph& og,
				       const ModelGraph& mg,
				       MV mv1, MV mv2,
				       const L& lops)
    {
	assert(mv1 < m2o.size());
	OV ov1 = m2o[mv1];
    
	typedef typename gtraits<G>::AdjIP AdjIP;

	if (mv2 < m2o.size())
	{
	    OV ov2 = m2o[mv2];
	    for (AdjIP ap = adjacent_vertices(ov1, og); ap.first != ap.second; ++ap.first)
	    {
		if (*ap.first == ov2)
		{
		    result.push_back(m2o);
		    break;
		}
	    }
	}
	else
	{
	    const VertexLabel& lab2 = mg[mv2];
	    for (AdjIP ap = adjacent_vertices(ov1, og); ap.first != ap.second; ++ap.first)
	    {
		bool skip = false;
		for (typename M2O::const_iterator k = m2o.begin(); k != m2o.end(); ++k)
		    if (*ap.first == *k)
			skip = true;
		if (skip) continue;
	    
		if (lops.label_equal(lab2, lops.label_value(*ap.first, og)))
		{
		    result.push_back(m2o);
		    result.back().push_back(*ap.first);
		    assert(result.back().size() - 1 == mv2);
		}
	    }
	}
    }


    template<class G, class L>
    void Candidate<G,L>::mk(const OGS& r_ogs, MV mv1, MV mv2)
    {
	for (typename OGS::const_iterator i = r_ogs.begin(); i != r_ogs.end(); ++i)
	{
	    M2OList sblist_new;
	    const OriGraph& og = *i->first;
	    for (typename M2OList::const_iterator j = i->second.begin();
		 j != i->second.end(); ++j)
	    {
		m2o_construct(sblist_new, *j, og, *this, mv1, mv2, lops);
	    }

	    if (!sblist_new.empty())
	    {
		ogs.push_back(std::make_pair(&og, M2OList()));
		ogs.back().second.swap(sblist_new);
	    }
	}

	if (ogs.empty())
	    throw NotFrequentException();
    }


    template<class G, class L>
    template<class OriGraphIterator>
    Candidate<G,L>::Candidate(OriGraphIterator begin,
			      OriGraphIterator end,
			      const VertexLabel& lab1, const VertexLabel& lab2,
			      const L& lops)
	: lops(lops)
    {
	ID = ++N_G;

	//
	// Create model graph
	//
	MV mv1 = add_vertex(*this);
	MV mv2 = add_vertex(*this);
	ME me = add_edge(mv1, mv2, *this).first;
	(*this)[mv1] = lab1;
	(*this)[mv2] = lab2;


	//
	// 
	//
	for (; begin != end; ++begin)
	{
	    M2OList sblist_new;
	    const OriGraph& og = *begin;
	    for (typename gtraits<OriGraph>::EIP ep = edges(og);
		 ep.first != ep.second; ++ep.first)
	    {
		OV ov1 = source(*ep.first, og);
		OV ov2 = target(*ep.first, og);

		if (lops.label_equal(lab1, lops.label_value(ov1, og)) &&
		    lops.label_equal(lab2, lops.label_value(ov2, og)))
		{
		    sblist_new.push_back(M2O());
		    sblist_new.back().push_back(ov1);
		    sblist_new.back().push_back(ov2);
		}

		if (lops.label_equal(lab1, lops.label_value(ov2, og)) &&
		    lops.label_equal(lab2, lops.label_value(ov1, og)))
		{
		    sblist_new.push_back(M2O());
		    sblist_new.back().push_back(ov2);
		    sblist_new.back().push_back(ov1);
		}
	    }
	
	    if (!sblist_new.empty())
	    {
		ogs.push_back(std::make_pair(&og, M2OList()));
		ogs.back().second.swap(sblist_new);
	    }
	}
    
	if (ogs.empty())
	    throw NotFrequentException();
    }


    template<class G, class L>
    Candidate<G,L>::Candidate(const Candidate& ancestor,
			      MV mv, const VertexLabel& lab)
	: ModelGraph(ancestor), lops(ancestor.lops)
    {
	ID = ++N_G;

	//
	// Create model graph
	//
	MV mv_new = add_vertex(*this);
	std::pair<ME, bool> r = add_edge(mv, mv_new, *this);
	assert (r.second);
	ME me = r.first;
	(*this)[mv_new] = lab;


	//
	//
	//
	mk(ancestor.ogs, mv, mv_new);

	if (ogs.empty())
	    throw NotFrequentException();
    }



    template<class G, class L>
    Candidate<G,L>::Candidate(const Candidate& ancestor, MV mv1, MV mv2)
	: ModelGraph(ancestor), lops(ancestor.lops)
    {
	ID = ++N_G;

	//
	// Create model graph
	//
	std::pair<ME, bool> r = add_edge(mv1, mv2, *this);
	if (!r.second)
	    throw NotFrequentException();
	ME me = r.first;


	//
	//
	//
	mk(ancestor.ogs, mv1, mv2);

	if (ogs.empty())
	    throw NotFrequentException();
    }



    template<class Cand, class Cont>
    void insert_new_candidate(Cont& cont,
			      const Cand& ancestor,
			      typename Cand::MV v,
			      const typename Cand::VertexLabel& lab)
    {
	try { cont.push_back(new Cand(ancestor, v, lab));
	} catch (NotFrequentException) {}
    }


    template<class Cand, class Cont>
    void insert_new_candidate(Cont& cont,
			      const Cand& ancestor,
			      typename Cand::MV v,
			      typename Cand::MV u)
    {
	try { cont.push_back(new Cand(ancestor, v, u));
	} catch (NotFrequentException) {}
    }



    /////////////////////////////////////////////////////
    //		class ConnectedSubgraph
    // -------------------------------------------------
    //
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

    

    /////////////////////////////////////////////////////
    //		class GraphCollection
    // -------------------------------------------------
    // 
    template<class G>
    class GraphCollection : public std::list<G*>
    {
    public:
	~GraphCollection()
	    {
		typedef std::list<G*> Base;
		for (typename std::list<G*>::iterator i = Base::begin();
		     i != Base::end(); ++i)
		    delete *i;
	    }
    };
    

    /////////////////////////////////////////////////////
    //		class LatticeLink
    // -------------------------------------------------
    // 
    template<class G>
    class LatticeLink
    {
    public:
	typedef typename boost::graph_traits<G>::vertex_descriptor V;

	const G* ancestor;

	// isomap and additional edge
	struct M
	{
	    std::vector<V> a2d; // map ancestor vertices to descendant
	    V operator[] (V v_a) const
		{ return a2d[v_a]; }

	    struct ExtEdge
	    {
		V vertex_descendant_1;
		V vertex_descendant_2;
		V vertex_ancestor_1;
		V vertex_ancestor_2; // may be null_vertex()
	    } e;
	};
	typedef std::list<M> MList;
	typedef typename MList::const_iterator MI;
	MList mm;
	
	void swap(LatticeLink& r)
	    {
		std::swap(ancestor, r.ancestor);
		mm.swap(r.mm);
	    }
    };



    /////////////////////////////////////////////////////
    //		class LatticeLayer
    // -------------------------------------------------
    // lattice layer
    // map Descendant graph to multiple Ancestor graphs
    //
    template<class G>
    struct LatticeLayer : public std::multimap<const G*, LatticeLink<G> > {};
    

    /////////////////////////////////////////////////////
    //		construct LatticeLayer
    // -------------------------------------------------
    //
    template<class G, class L>
    bool insert(LatticeLayer<G>& lat,
		const G& descendant,
		const GraphCollection<G>& ancestors,
		const L& lops)
    {
	bool ret = false;
	typedef ConnectedSubgraph<G> ConnSubGraph;

	// for each ancestor graph
	for (typename GraphCollection<G>::const_iterator i = ancestors.begin();
	     i != ancestors.end(); ++i)
	{
	    LatticeLink<G> li;
	    li.ancestor = *i;
	    
	    ConnSubGraph ancestor_graph(**i);


	    typedef typename boost::graph_traits<G>::vertex_descriptor V;
	    typedef typename boost::graph_traits<G>::edge_descriptor E;
	    typedef typename boost::graph_traits<G>::edge_iterator EI;
	    typedef std::pair<EI,EI> EIP;
	    for (EIP eip = edges(descendant); eip.first != eip.second; ++eip.first)
	    {
		E edge = *eip.first;

		try
		{
		    ConnSubGraph descendant_subgraph(descendant, edge);

		    // find all isomorphisms
		    typedef std::vector<MPair> Sec;
		    typedef std::vector<Sec> Sec2;
		    Sec2 sec2;
		    connsb_isomorphism(ancestor_graph,
				       descendant_subgraph, lops, sec2);
		    

		    typedef Sec2::const_iterator II;
		    for (II ii = sec2.begin(); ii != sec2.end(); ++ii)
		    {

			typedef typename LatticeLink<G>::M M;
			li.mm.push_back(M());
			M& m = li.mm.back();

			const V NULV = boost::graph_traits<G>::null_vertex();
			m.e.vertex_descendant_1 = source(edge, descendant);
			m.e.vertex_descendant_2 = target(edge, descendant);
			m.e.vertex_ancestor_1 = NULV;
			m.e.vertex_ancestor_2 = NULV;
			m.a2d.resize(num_vertices(**i), NULV);
			for (Sec::const_iterator j = ii->begin(); j != ii->end(); ++j)
			{
			    m.a2d[j->first] = j->second;
			    
			    if (m.e.vertex_descendant_1 == static_cast<V>(j->second))
				m.e.vertex_ancestor_1 = j->first;
			    if (m.e.vertex_descendant_2 == static_cast<V>(j->second))
				m.e.vertex_ancestor_2 = j->first;
			}
			
			if (m.e.vertex_ancestor_1 == NULV)
			{
			    std::swap(m.e.vertex_descendant_1, m.e.vertex_descendant_2);
			    std::swap(m.e.vertex_ancestor_1, m.e.vertex_ancestor_2);
			}

			assert(m.e.vertex_descendant_1 != NULV);
			assert(m.e.vertex_descendant_2 != NULV);
			assert(m.e.vertex_ancestor_1 != NULV);
			assert(find(m.a2d.begin(), m.a2d.end(), NULV) == m.a2d.end());
		    }
		}
		catch (NotConnectedException)
		{
		    if (descendant.ID == 1398 && (*i)->ID == 1128)
		    {
			if ((source(edge, descendant) == 0 && target(edge, descendant) == 1) ||
			    (source(edge, descendant) == 1 && target(edge, descendant) == 0))
			    ;//BR;
		    }

		}
	    }

	    if (! li.mm.empty())
	    {
		lat.insert(std::make_pair(&descendant, LatticeLink<G>()))
		    ->second.swap(li);
		ret = true;
	    }
	}
	
	return ret;
    }


    /////////////////////////////////////////////////////
    //		class Lattice
    // -------------------------------------------------
    //
    template<class G>
    class Lattice : public std::list<GraphCollection<G> >
    {
    public:
	LatticeLayer<G> ll;
    };



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
		VertexCompatible<CandidateGraph,L> vc(**i, **j, lops);
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
