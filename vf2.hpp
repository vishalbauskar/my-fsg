#ifndef VF2_H_
#define VF2_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <algorithm>
#include <utility>
#include <map>

namespace graph_alg
{
    typedef int VertexID;
    const VertexID NULL_VERTEX = -1;
    
    typedef std::pair<VertexID,VertexID> MPair;

    namespace detail
    {
	// ===========================================================
	//	class EdgesMap
	// ===========================================================
	template<class G>
	class EdgesMap
	{
	public:
	    typedef typename G::edge_descriptor E;
	    typedef typename G::vertex_descriptor V;
	    EdgesMap(const G& g)
		{
		    typedef typename boost::graph_traits<G>::edge_iterator EI;
		    for (std::pair<EI,EI> i = edges(g);
			 i.first != i.second; ++i.first)
			val(source(*i.first, g), target(*i.first, g)) = *i.first;
		}

	    std::pair<E,bool> operator() (V v1, V v2) const
		{
		    typename std::map<V, std::map<V, E> >::const_iterator
			i = m.find(std::min(v1,v2));
		    if (i != m.end())
		    {
			typename std::map<V, E>::const_iterator
			    ii = i->second.find(std::max(v1,v2));
			if (ii != i->second.end())
			    return std::pair<E,bool>(ii->second, true);
		    }
		    return std::pair<E,bool>(E(), false);
		}
	private:
	    std::map<V, std::map<V, E> > m;
	    E& val(V v1, V v2) { return m[std::min(v1,v2)][std::max(v1,v2)]; }
	    const E& val(V v1, V v2) const { return m[std::min(v1)][std::max(v2)]; }
	};


	// ===========================================================
	//	class VertexEnumerator
	// ===========================================================	
	template<class G>
	class VertexEnumerator
	{
	    typedef typename boost::graph_traits<G>::vertex_descriptor Vertex;
	    typedef typename boost::graph_traits<G>::vertex_iterator VertexI;

	    const VertexID* c;
	    const int* t;
	    int n;

	    std::pair<VertexI,VertexI> viter;
	    bool has;
	    bool is_unmapped_vertex(Vertex v) const
		{ return t ?
			c[v] == NULL_VERTEX && t[v] :
			c[v] == NULL_VERTEX; }

	    bool find_next()
		{
		    //if (n)BR;
		    while (n && viter.first != viter.second)
		    {
			if (is_unmapped_vertex(*viter.first))
			    return true;
			++viter.first;
		    }
		    return false;
		}
	public:
	    VertexEnumerator(const G& g,
			      const VertexID* c = 0,
			      const int* t = 0, int n = 0)
		:c(c), t(t), n(n), viter(vertices(g)), has(find_next()) {}
	    VertexID get_v() const { assert(has); return *viter.first; }
	    bool has_vertex() const { return has; }
	    void next() { assert(has); ++viter.first; has = find_next(); }
	};
	

	
	// ===========================================================
	//	class State
	// ===========================================================
	template<class G,
		 class VertexCompatible,
		 class EdgeCompatible>
	class VF2State
	{
	    struct MapData {
		struct GData {
		    VertexID* core;
		    int* out;
		    int* in;
		    const G& graph;
		    explicit GData(const G&);
		    ~GData();
		    int num() const { return num_vertices(graph); }
		    void prt(int y) const
			{
			    std::cerr<<"t_"<<y<<": ";
			    for (int i=0; i < num(); i++)
				std::cerr<< "("<<out[i]<<","<<in[i]<<") ";
			    std::cerr<<std::endl;
			}
		};
		GData t1;
		GData t2;

		const VertexCompatible& vc;
		const EdgeCompatible& ec;

		EdgesMap<G> edgemap;

		MapData(const G& g1,const G& g2,
			const VertexCompatible& vc,
			const EdgeCompatible& ec) : t1(g1), t2(g2), vc(vc), ec(ec),
						    edgemap(g2)
		    {}

		void prt() const { t1.prt(1); t2.prt(2); }
	    };

	    boost::shared_ptr<MapData> m;
	    
	    struct TLen {
		int out;
		int in;
		TLen() : out(0), in(0) {}
	    };

	    TLen t1len;
	    TLen t2len;

	    const int orig_core_len;
	    int core_len;
	    MPair last_pair;

	    const int n1, n2;

	    static void add_pair(VertexID v1, VertexID v2,
				 typename MapData::GData& tp,
				 TLen& tlen, int core_len);
	    static void backtrack(VertexID v,
				  typename MapData::GData& tp,
				  int core_len);
	public:
	    void prt() const { m->prt(); }
	    class MPairEnumerator
	    {
		const VF2State& s;
		enum St { T_OUT, T_IN, T_V };
		St st;

		typedef VertexEnumerator<G> VEnum;
		VEnum ven1;
		VEnum ven2;
	    public:
		MPairEnumerator(const VF2State& s);
		MPair get_mpair() const { return MPair(ven1.get_v(), ven2.get_v()); }
		void next();
		bool has_mpair() const
		    { return ven1.has_vertex() && ven2.has_vertex(); }
	    };
	    friend class MPairEnumerator;

	    VF2State(const G& g1, const G& g2,
		     const VertexCompatible& vc,
		     const EdgeCompatible& ec);

	    VF2State(const G& g1, const G& g2,
		     int nv1, int nv2,
		     const VertexCompatible& vc,
		     const EdgeCompatible& ec);

	    VF2State(const MPair&, const VF2State&);

	    bool is_goal() const;
	    bool is_dead() const;
	    MPairEnumerator enumerate_pairs() const;
	    bool is_feasible(const MPair&) const;
	    void backtrack();

	    template<class OutputIterator>
	    void yield(OutputIterator iout) const;

	    template<class IsoMap>
	    void yield2(IsoMap& isom) const;
	};



	// ===========================================================
	//	class VF2State::MapData::GData
	// ===========================================================
	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::MapData::GData::GData(const G& g)
	    :graph(g)
	{
	    core = new VertexID[num()];
	    out = new int[num()];
	    in = new int[num()];

	    std::fill_n(core, num(), NULL_VERTEX);
	    std::fill_n(out, num(), 0);
	    std::fill_n(in, num(), 0);
	}

	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::MapData::GData::~GData()
	{
	    delete[] core;
	    delete[] out;
	    delete[] in;
	}

	// ===========================================================
	//	class VF2State::MPairEnumerator
	// ===========================================================
	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::MPairEnumerator::MPairEnumerator(const VF2State<G,VC,EC>& s)
	    : s(s),
	      ven1(VEnum(s.m->t1.graph)),
	      ven2(VEnum(s.m->t2.graph))
	{
	    if (s.t1len.out > s.core_len && s.t2len.out > s.core_len)
	    {
		st = T_OUT;
		ven1 = VEnum(s.m->t1.graph, s.m->t1.core, s.m->t1.out, s.m->t1.num());
		ven2 = VEnum(s.m->t2.graph, s.m->t2.core, s.m->t2.out, s.m->t2.num());
	    }
	    else if (s.t1len.in > s.core_len && s.t2len.in > s.core_len)
	    {
		st = T_IN;
		ven1 = VEnum(s.m->t1.graph, s.m->t1.core, s.m->t1.in, s.m->t1.num());
		ven2 = VEnum(s.m->t2.graph, s.m->t2.core, s.m->t2.in, s.m->t2.num());
	    }
	    else
	    {
		st = T_V;
		ven1 = VEnum(s.m->t1.graph, s.m->t1.core, 0, s.m->t1.num());
		ven2 = VEnum(s.m->t2.graph, s.m->t2.core, 0, s.m->t2.num());
	    }
	}

	template<class G, class VC, class EC>
	void VF2State<G,VC,EC>::MPairEnumerator::next()
	{
	    ven2.next();
	}


	// ===========================================================
	//	class VF2State
	// ===========================================================
	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::VF2State(const G& g1, const G& g2,
				    const VC& vc,
				    const EC& ec)
	    : m(new MapData(g1,g2, vc, ec)),
	      orig_core_len(0), core_len(0),
	      last_pair(NULL_VERTEX,NULL_VERTEX),
	      n1(num_vertices(g1)), n2(num_vertices(g2))
	{
	}


	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::VF2State(const G& g1, const G& g2,
				    int nv1, int nv2,
				    const VC& vc,
				    const EC& ec)
	    : m(new MapData(g1, g2, vc, ec)),
	      orig_core_len(0), core_len(0),
	      last_pair(NULL_VERTEX,NULL_VERTEX), n1(nv1), n2(nv2)
	{
	}


	template<class G, class VC, class EC>
	VF2State<G,VC,EC>::VF2State(const MPair& p, const VF2State<G,VC,EC>& s)
	    : m(s.m),
	      t1len(s.t1len), t2len(s.t2len),
	      orig_core_len(s.core_len), core_len(s.core_len),
	      last_pair(p), n1(s.n1), n2(s.n2)
	{
	    ++core_len;
	    add_pair(p.first, p.second, m->t1, t1len, core_len);
	    add_pair(p.second, p.first, m->t2, t2len, core_len);
	}


	template<class G, class VC, class EC>
	void VF2State<G,VC,EC>::add_pair(VertexID v1, VertexID v2,
					 typename MapData::GData& tp,
					 TLen& tlen, int core_len)
	{
	    tp.core[v1] = v2;

	    if (! tp.out[v1])
	    {
		tp.out[v1] = core_len;
		++tlen.out;
	    }

	    if (! tp.in[v1])
	    {
		tp.in[v1] = core_len;
		++tlen.in;
	    }

	    typedef typename boost::graph_traits<G>::out_edge_iterator OutI;
	    for (std::pair<OutI,OutI> ap = out_edges(v1, tp.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		VertexID out_vi = target(*ap.first, tp.graph);
		if (! tp.out[out_vi])
		{
		    tp.out[out_vi] = core_len;
		    ++tlen.out;
		}
	    }


	    typedef typename boost::graph_traits<G>::in_edge_iterator InI;
	    for (std::pair<InI,InI> ap = in_edges(v1, tp.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		VertexID in_vi = source(*ap.first, tp.graph);
		if (! tp.in[in_vi])
		{
		    tp.in[in_vi] = core_len;
		    ++tlen.in;
		}
	    }
	}


	template<class G, class VC, class EC>
	void VF2State<G,VC,EC>::backtrack(VertexID v,
					  typename MapData::GData& tp,
					  int core_len)
	{
	    if (tp.out[v] == core_len) tp.out[v] = 0;
	    if (tp.in[v] == core_len) tp.in[v] = 0;

	    typedef typename boost::graph_traits<G>::out_edge_iterator OutI;
	    for (std::pair<OutI,OutI> ap = out_edges(v, tp.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		VertexID out_vi = target(*ap.first, tp.graph);
		if (tp.out[out_vi] == core_len) tp.out[out_vi] = 0;
	    }

	    typedef typename boost::graph_traits<G>::in_edge_iterator InI;
	    for (std::pair<InI,InI> ap = in_edges(v, tp.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		VertexID in_vi = source(*ap.first, tp.graph);
		if (tp.in[in_vi] == core_len) tp.in[in_vi] = 0;
	    }

	    tp.core[v] = NULL_VERTEX;
	}


	template<class G, class VC, class EC>
	void VF2State<G,VC,EC>::backtrack()
	{
	    if (core_len > orig_core_len)
	    {
		backtrack(last_pair.first, m->t1, core_len);
		backtrack(last_pair.second, m->t2, core_len);
		last_pair = MPair(NULL_VERTEX,NULL_VERTEX);
		core_len = orig_core_len;
	    }
	}


	template<class G, class VC, class EC>
	bool VF2State<G,VC,EC>::is_goal() const
	{
	    return core_len == n1;
	}


	template<class G, class VC, class EC>
	bool VF2State<G,VC,EC>::is_dead() const
	{
	    return
		m->t1.num() > m->t2.num() ||
		t1len.out > t2len.out ||
		t1len.in > t2len.in;
	}


	template<class G, class VC, class EC>
	typename VF2State<G,VC,EC>::MPairEnumerator VF2State<G,VC,EC>::enumerate_pairs() const
	{
	    return MPairEnumerator(*this);
	}


	template<class G, class VC, class EC>
	bool VF2State<G,VC,EC>::is_feasible(const MPair& p) const
	{
	    const VertexID v1 = p.first;
	    const VertexID v2 = p.second;

	    assert(v1 < m->t1.num());
	    assert(v2 < m->t2.num());
	    assert(m->t1.core[v1] == NULL_VERTEX);
	    assert(m->t2.core[v2] == NULL_VERTEX);

	    if (! m->vc(v1, v2))
		return false;


	    int termout1 = 0, termout2 = 0, termin1 = 0, termin2 = 0;
	    unsigned int new1 = 0, new2 = 0;

	    typedef typename G::edge_descriptor E;
	    typedef typename boost::graph_traits<G>::out_edge_iterator OutI;
	    for (std::pair<OutI,OutI> ap = out_edges(v1, m->t1.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		E edge_1 = *ap.first;
		VertexID adj_v1 = target(edge_1, m->t1.graph);
		if (m->t1.core[adj_v1] != NULL_VERTEX)
		{
		    VertexID adj_v2 = m->t1.core[adj_v1];
		    E edge_2;
		    bool found;
		    boost::tie(edge_2,found) = m->edgemap(v2, adj_v2);
		    if (!found || !m->ec(edge_1, edge_2))
			return false;
		}
		else
		{
		    if (m->t1.in[adj_v1]) ++termin1;
		    if (m->t1.out[adj_v1]) ++termout1;
		    if (!m->t1.in[adj_v1] && !m->t1.out[adj_v1]) ++new1;
		}
	    }


	    for (std::pair<OutI,OutI> ap = out_edges(v2, m->t2.graph);
		 ap.first != ap.second; ++ap.first)
	    {
		E edge_2 = *ap.first;
		VertexID adj_v2 = target(edge_2, m->t1.graph);
		if (m->t2.core[adj_v2] == NULL_VERTEX)
		{
		    if (m->t2.in[adj_v2]) ++termin2;
		    if (m->t2.out[adj_v2]) ++termout2;
		    if (!m->t2.in[adj_v2] && !m->t2.out[adj_v2]) ++new2;
		}
	    }

	    return termin1 <= termin2 && termout1 <= termout2 &&
		(termin1+termout1+new1)<=(termin2+termout2+new2);
	}

	template<class G, class VC, class EC>
	template<class OutputIterator>
	void VF2State<G,VC,EC>::yield(OutputIterator iout) const
	{
	    for (VertexID i = 0; i < m->t1.num(); ++i)
	    {
		if (m->t1.core[i] != NULL_VERTEX)
		{
		    *iout = MPair(i, m->t1.core[i]);
		    ++iout;
		}
	    }
	}


	template<class G, class VC, class EC>
	template<class IsoMap>
	void VF2State<G,VC,EC>::yield2(IsoMap& isom) const
	{
	    for (VertexID i = 0; i < m->t1.num(); ++i)
		if (m->t1.core[i] != NULL_VERTEX)
		    isom[i] = m->t1.core[i];
	}


	// ===========================================================
	//	match()
	// ===========================================================
	template<class G, class VC, class EC, class OutputIterator>
	bool match_one(VF2State<G,VC,EC>& s, OutputIterator& iout)
	{
	    if (s.is_goal())
	    {
		s.yield(iout);
		return true;
	    }

	    if (s.is_dead())
		return false;

	    bool found = false;
	    typedef typename VF2State<G,VC,EC>::MPairEnumerator MPE;
	    for (MPE e = s.enumerate_pairs(); !found && e.has_mpair(); e.next())
	    {
		const MPair p = e.get_mpair();	
		if (s.is_feasible(p))
		{
		    VF2State<G,VC,EC> h(p, s);
		    found = match_one(h, iout);
		    h.backtrack();
		}
	    }
	    return found;
	}


	template<class G,
		 class VC, class EC,
		 class DoubleBackInsertionSequence>
	void match_all(VF2State<G,VC,EC>& s, DoubleBackInsertionSequence& sec)
	{

	    if (s.is_goal())
	    {
		typedef typename DoubleBackInsertionSequence::value_type OneSecT;
		sec.push_back(OneSecT());
		s.yield(back_inserter(sec.back()));
		return;
	    }

	    if (s.is_dead())
		return;

	    typedef typename VF2State<G,VC,EC>::MPairEnumerator MPE;
	    for (MPE e = s.enumerate_pairs(); e.has_mpair(); e.next())
	    {
		const MPair p = e.get_mpair();
		if (s.is_feasible(p))
		{
		    VF2State<G,VC,EC> h(p, s);
		    match_all(h, sec);
		    h.backtrack();
		}
	    }
	}

	template<class G,
		 class VC, class EC,
		 class IsoMaps>
	void match_all2(VF2State<G,VC,EC>& s, IsoMaps& isom)
	{
	    if (s.is_goal())
	    {
		isom.push_back(typename IsoMaps::value_type());
		s.yield2(isom.back());
		return;
	    }

	    if (s.is_dead())
		return;

	    typedef typename VF2State<G,VC,EC>::MPairEnumerator MPE;
	    for (MPE e = s.enumerate_pairs(); e.has_mpair(); e.next())
	    {
		const MPair p = e.get_mpair();	
		if (s.is_feasible(p))
		{
		    VF2State<G,VC,EC> h(p, s);
		    match_all2(h, isom);
		    h.backtrack();
		}
	    }
	}

    } // namespace detail



    template<class G>
    struct VertexDefaultCompatible
    {
	typedef typename boost::graph_traits<G>::vertex_descriptor VDecr;
	bool operator() (VDecr, VDecr) const { return true; }
    };

    template<class G>
    struct EdgeDefaultCompatible
    {
	typedef typename boost::graph_traits<G>::edge_descriptor EDecr;
	bool operator() (EDecr, EDecr) const { return true; }
    };

    // ===========================================================
    //	isomorphism_one()
    // ===========================================================

    template<class G,
	     class OutputIterator,
	     class VertexCompatible,
	     class EdgeCompatible>
    bool isomorphism_one(const G& g1, const G& g2, OutputIterator iout,
			const VertexCompatible& vc, const EdgeCompatible& ec)
    {
	if (num_vertices(g1) != num_vertices(g2)) return false;
	detail::VF2State<G,VertexCompatible,EdgeCompatible> s0(g1, g2, vc, ec);
	return match_one(s0, iout);
    }

    template<class G,
	     class OutputIterator,
	     class VertexCompatible>
    bool isomorphism_one(const G& g1, const G& g2, OutputIterator iout,
			 const VertexCompatible& vc)
    {
	if (num_vertices(g1) != num_vertices(g2)) return false;

	detail::VF2State<G, VertexCompatible, EdgeDefaultCompatible<G> >
	    s0(g1, g2, vc, EdgeDefaultCompatible<G>());
	return match_one(s0, iout);
    }

    template<class G,
	     class OutputIterator>
    bool isomorphism_one(const G& g1, const G& g2, OutputIterator iout)
    {
	if (num_vertices(g1) != num_vertices(g2)) return false;

	detail::VF2State<G, VertexDefaultCompatible<G>, EdgeDefaultCompatible<G> >
	    s0(g1, g2);
	return match_one(s0, iout);
    }

    // ===========================================================
    //	isomorphism_test
    // ===========================================================
    template<class G,
	     class VertexCompatible>
    bool isomorphism_test(const G& g1, const G& g2,
			  const VertexCompatible& vc)
    {
	std::list<MPair> foo;
	return isomorphism_one(g1, g2, back_inserter(foo), vc);
    }

    // ===========================================================
    //	isomorphism_all()
    // ===========================================================

    template<class G,
	     class DoubleBackInsertionSequence,
	     class VertexCompatible,
	     class EdgeCompatible>
    bool isomorphism_all(const G& g1, const G& g2,
			 DoubleBackInsertionSequence& sec,
			 const VertexCompatible& vc,
			 const EdgeCompatible& ec,
			 int nv1 = -1, int nv2 = -1)
    {
	if (nv1 < 0) nv1 = num_vertices(g1);
	if (nv2 < 0) nv2 = num_vertices(g2);
	if (nv1 != nv2) return false;

	detail::VF2State<G,VertexCompatible,EdgeCompatible>
	    s0(g1, g2, nv1, nv2, vc, ec);
	std::size_t size = sec.size();
	match_all(s0, sec);
	return size != sec.size();
    }

    template<class G,
	     class DoubleBackInsertionSequence,
	     class VertexCompatible>
    bool isomorphism_all(const G& g1, const G& g2,
			 DoubleBackInsertionSequence& sec,
			 const VertexCompatible& vc,
			 int nv1 = -1, int nv2 = -1)
    {
	if (nv1 < 0) nv1 = num_vertices(g1);
	if (nv2 < 0) nv2 = num_vertices(g2);
	if (nv1 != nv2) return false;

	detail::VF2State<G,VertexCompatible,EdgeDefaultCompatible<G> >
	    s0(g1, g2, nv1, nv2, vc, EdgeDefaultCompatible<G>());
	std::size_t size = sec.size();
	match_all(s0, sec);
	return size != sec.size();
    }


    template<class G,
	     class DoubleBackInsertionSequence>
    bool isomorphism_all(const G& g1, const G& g2,
			 DoubleBackInsertionSequence& sec,
			 int nv1 = -1, int nv2 = -1)
    {
	if (nv1 < 0) nv1 = num_vertices(g1);
	if (nv2 < 0) nv2 = num_vertices(g2);
	if (nv1 != nv2) return false;

	detail::VF2State<G,VertexDefaultCompatible<G>,EdgeDefaultCompatible<G> >
	    s0(g1, g2, nv1, nv2,
	       VertexDefaultCompatible<G>(), EdgeDefaultCompatible<G>());
	std::size_t size = sec.size();
	match_all(s0, sec);
	return size != sec.size();
    }



    // ===========================================================
    //	isomorphism_all2()
    // ===========================================================

    template<class G,
	     class IsoMaps,
	     class VertexCompatible>
    bool isomorphism_all2(const G& g1, const G& g2,
			  IsoMaps& isom,
			  const VertexCompatible& vc,
			  int nv1 = -1, int nv2 = -1)
    {
	if (nv1 < 0) nv1 = num_vertices(g1);
	if (nv2 < 0) nv2 = num_vertices(g2);

	if (nv1 != nv2) return false;

	detail::VF2State<G,VertexCompatible,EdgeDefaultCompatible<G> >
	    s0(g1, g2, nv1, nv2, vc, EdgeDefaultCompatible<G>());
	//std::size_t size = sec.size();
	match_all2(s0, isom);
	//return size != sec.size();
	return 1;
    }


} // namespace vf2

#endif
