#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
#include <map>
#include <list>
#include <cassert>
#include <algorithm>

#include <boost/graph/filtered_graph.hpp>

#include "vf2.hpp"
#include "connsubgraph.hpp"

namespace graph_alg
{  
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

		if (descendant.ID == 1398 && (*i)->ID == 1128)
		{
		    if ((source(edge, descendant) == 0 && target(edge, descendant) == 1) ||
			(source(edge, descendant) == 1 && target(edge, descendant) == 0))
		    {
			//BR;
			print_graph(descendant);
			std::cerr << "Edge = " << edge << "\n";
		    }
		}

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
    
}
#endif
