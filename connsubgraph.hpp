#ifndef CONNSUBGRAPH_H_
#define CONNSUBGRAPH_H_

#include "vf2.hpp"

namespace graph_alg
{

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




    template<class G, class L>
    struct VertexCompatible_
    {
	VertexCompatible_(const G& g1,
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


    template<class G, class L, class TwoMPairSec>
    bool connsb_isomorphism(const ConnectedSubgraph<G>& g1,
			    const ConnectedSubgraph<G>& g2,
			    const L& lops,
			    TwoMPairSec& sec2)
    {
	typedef typename ConnectedSubgraph<G>::FilGraph FilGraph;
	VertexCompatible_<FilGraph, L> vc(g1, g2, lops);
	
	return isomorphism_all(*g1.fil, *g2.fil,
			       sec2,
			       vc,
			       g1.num_vertices(),
			       g2.num_vertices());
    }
}

#endif
