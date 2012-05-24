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

#include "fsg.hpp"


template<class G, class L>
class VizFSGLattice : public graph_alg::FSG_Lattice<G,L>
{
    typedef graph_alg::FSG_Lattice<G,L> Base;
public:
    template<class GI>
    VizFSGLattice(GI gi_begin, GI gi_end, const L& l = L());
    void grviz_write() const;
private:
    
    typedef graph_alg::detail::Candidate<G,L> CandidateGraph;

    typedef boost::adjacency_list<boost::vecS,
				  boost::vecS,
				  boost::directedS,
				  std::pair<const CandidateGraph*, std::string> > LG;
    typedef typename boost::graph_traits<LG>::vertex_descriptor V;
    typedef typename boost::graph_traits<LG>::vertex_iterator VI;

    LG lat_graph;
    std::map<const CandidateGraph*, V> g2v;
};


template<class G, class L>
template<class GI>
VizFSGLattice<G,L>::VizFSGLattice(GI gi_begin, GI gi_end, const L& l)
    : Base(gi_begin, gi_end, l)
{
    // create vertices
    for (typename Base::CC_CI i = Base::lattice.begin(); i != Base::lattice.end(); ++i)
	for (typename Base::CG_CI j = i->begin(); j != i->end(); ++j)
	{
	    std::ostringstream file_name;
	    file_name << "g" << j->id;

	    V v = add_vertex(lat_graph);
	    g2v[&*j] = v;
	    lat_graph[v] = std::make_pair(&*j, file_name.str());
	}

    // create edges
    typename std::map<const CandidateGraph*, V>::const_iterator k;
    for (k = g2v.begin(); k != g2v.end(); ++k)
    {
	typename CandidateGraph::P_CI h;
	for (h = k->first->parents.begin(); h != k->first->parents.end(); ++h)
	    add_edge(k->second, g2v[h->parent], lat_graph);
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
	out << "[label=" << name[v].first->id << "]";
    }
private:
    Name name;
};


template<class G, class L>
void VizFSGLattice<G,L>::grviz_write() const
{
    ::system("rm -f grviz/*");

    for (std::pair<VI,VI> vip = vertices(lat_graph); vip.first != vip.second;
	 ++vip.first)
    {
	const std::pair<const CandidateGraph*, std::string>& pp = lat_graph[*vip.first];
	
	const CandidateGraph& g = *pp.first;
	std::filebuf dot;
	std::string dot_file = "grviz/" + pp.second;
	std::filebuf* p = dot.open((dot_file + ".dot").c_str(), std::ios::out);
	assert(p);
	std::ostream os(&dot);
	write_graphviz(os, g, my_label_writer<CandidateGraph>(g));

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

#endif

