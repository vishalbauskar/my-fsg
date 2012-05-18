#ifndef VIZLATTICE_H_
#define VIZLATTICE_H_

#include <stdlib.h>

#include <sstream>
#include <fstream>
#include <utility>
#include <string>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

template<class G>
class LatGraph
{
public:
    typedef boost::adjacency_list<boost::vecS,
				  boost::vecS,
				  boost::directedS,
				  std::pair<const G*, std::string> > LG;
    typedef typename boost::graph_traits<LG>::vertex_descriptor V;
    typedef typename boost::graph_traits<LG>::vertex_iterator VI;

    LG lat_graph;
    std::map<const G*, V> g2v;

    LatGraph(const graph_alg::Lattice<G>& l);
    
    void write_gviz() const;
};


template<class G>
LatGraph<G>::LatGraph(const graph_alg::Lattice<G>& l)
{
    // create vertices
    typedef typename graph_alg::Lattice<G>::const_iterator LI;
    typedef typename graph_alg::GraphCollection<G> GColl;
    typedef typename GColl::const_iterator GI;
    for (LI li = l.begin(); li != l.end(); ++li)
	for (GI gi = li->begin(); gi != li->end(); ++gi)
	{
	    std::ostringstream file_name;
	    file_name << "g" << (*gi)->ID;

	    V v = add_vertex(lat_graph);
	    g2v[*gi] = v;
	    lat_graph[v] = std::make_pair(*gi, file_name.str());
	}
    
    // create edges
    typedef typename graph_alg::LatticeLayer<G>::const_iterator LLI;
    for (LLI lli = l.ll.begin(); lli != l.ll.end(); ++lli)
    {
	const G* g_descen = lli->first;
	const G* g_ascen  = lli->second.ancestor;
	typename std::map<const G*, V>::const_iterator g2v_i;

	g2v_i = g2v.find(g_descen);
	assert(g2v_i != g2v.end());
	V v_descen = g2v_i->second;
	
	g2v_i = g2v.find(g_ascen);
	assert(g2v_i != g2v.end());	
	V v_ascen = g2v_i->second;
	
	bool edge_added = add_edge(v_descen, v_ascen, lat_graph).second;
	assert(edge_added);

	v_descen=v_descen;
	v_ascen=v_ascen;
    }
}

template <class Name>
class my_label_writer {
public:
    my_label_writer(Name _name) : name(_name) {}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {
	out << "[label=\"" << v << name[v] << "\"]";
    }
private:
    Name name;
};


template <class Name>
class my_image_writer {
public:
    my_image_writer(Name _name) : name(_name) {}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {
	out << "[image=\"" << name[v].second + ".png" << "\"]";
	out << "[label=" << name[v].first->ID << "]";
    }
private:
    Name name;
};

template<class G>
void LatGraph<G>::write_gviz() const
{
    for (std::pair<VI,VI> vip = vertices(lat_graph); vip.first != vip.second;
	 ++vip.first)
    {
	const std::pair<const G*, std::string>& pp = lat_graph[*vip.first];
	
	const G& g = *pp.first;
	std::filebuf dot;
	std::string dot_file = "grviz/" + pp.second;
	std::filebuf* p = dot.open((dot_file + ".dot").c_str(), std::ios::out);
	assert(p);
	std::ostream os(&dot);
	write_graphviz(os, g, my_label_writer<G>(g));

	dot.close();

	::system(("dot -Tpng " + dot_file + ".dot > " + dot_file + ".png").c_str());
    }


    std::filebuf dot;
    std::filebuf* p = dot.open("grviz/graph.dot", std::ios::out);
    assert(p);
    std::ostream os(&dot);
    write_graphviz(os, lat_graph, my_image_writer<LG>(lat_graph));
    dot.close();

    ::system("cd grviz && dot -Tpng graph.dot > graph.png");
}


template<class G>
void pr_lattice(const graph_alg::Lattice<G>& l)
{
    ::system("rm -f /home/dedal/programming/vf2/grviz/*");
    LatGraph<G> latg(l);
    latg.write_gviz();
}

#endif

