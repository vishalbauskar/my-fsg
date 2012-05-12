#ifndef VF2_H_
#define VF2_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <algorithm>
#include <utility>
#include <map>

#define BR asm("int $3;")

namespace graph_alg
{
    typedef int VertexID;
    const VertexID NULL_VERTEX = -1;
    
    typedef std::pair<VertexID,VertexID> MPair;

    namespace detail
    {

	// ===========================================================
	//
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
	class VertexEnumerator2
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
	    VertexEnumerator2(const G& g,
			     const VertexID* c = 0,
			     const int* t = 0, int n = 0)
		:c(c), t(t), n(n), viter(vertices(g)), has(find_next()) {}
	    VertexID get_v() const { assert(has); return *viter.first; }
	    bool has_vertex() const { return has; }
	    void next() { assert(has); ++viter.first; has = find_next(); }
	};
	

	class VertexEnumerator
	{
	    const VertexID* c;
	    const int* t;
	    int n;
	    VertexID v;
	    bool has;
	    bool is_unmapped_vertex(VertexID v) const
		{ return t ?
			c[v] == NULL_VERTEX && t[v] :
			c[v] == NULL_VERTEX; }
	    bool find_next();
	public:
	    template<class G>
	    VertexEnumerator(const  G&, const VertexID* c = 0, const int* t = 0, int n = 0)
		:c(c), t(t), n(n), v(NULL_VERTEX), has(find_next()) {}
	    VertexID get_v() const { assert(has); return v; }
	    bool has_vertex() const { return has; }
	    void next() { assert(has); has = find_next(); }
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

		typedef VertexEnumerator2<G> VEnum;
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

#include "vf2_impl.hpp"
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
