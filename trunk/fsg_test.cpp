
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include "vizlattice.hpp"

using namespace boost;
using namespace std;
using namespace graph_alg;
using namespace graph_alg::detail;


///////////////////////////////////////////////////////////
//  Original Graph
typedef char LabelType;
typedef adjacency_list<vecS,vecS, undirectedS, LabelType> Graph;
typedef graph_traits<Graph>::vertex_descriptor V;
typedef graph_traits<Graph>::edge_descriptor E;

///////////////////////////////////////////////////////////
//  LabelOps
LabelType labels[] = {'a', 'b', 'c'};

struct LabelOps {
    typedef LabelType label_type;
    const label_type& label_value(V v, const Graph& g) const
	{ return g[v]; }
    bool label_equal(const label_type& l1, const label_type& l2) const
	{ return l1 == l2; }
    
    typedef const label_type* label_iterator;
    label_iterator begin() const	{ return labels; }
    label_iterator end() const		{ return labels+3; }
};


typedef VizFSGLattice<Graph, LabelOps> FSGLattice;

void mk_g_1(Graph& g)
{
    V v1 = add_vertex(g); g[v1] = 'a';
    V v2 = add_vertex(g); g[v2] = 'b';
    V v3 = add_vertex(g); g[v3] = 'b';
    V v4 = add_vertex(g); g[v4] = 'a';
    V v5 = add_vertex(g); g[v5] = 'a';
    V v6 = add_vertex(g); g[v6] = 'c';
    
    add_edge(v1, v2, g);
    add_edge(v2, v3, g);
    add_edge(v3, v4, g);
    add_edge(v4, v1, g);
    add_edge(v4, v5, g);
    add_edge(v3, v6, g);
    add_edge(v5, v6, g);
    add_edge(v2, v5, g);
}

void mk_g_2(Graph& g)
{
    V v1 = add_vertex(g); g[v1] = 'c';
    V v2 = add_vertex(g); g[v2] = 'c';
    V v3 = add_vertex(g); g[v3] = 'c';
    V v4 = add_vertex(g); g[v4] = 'c';
    
/*
    add_edge(v1, v2, g);
    add_edge(v2, v3, g);
    add_edge(v1, v3, g);
    add_edge(v4, v3, g);
*/
}


int main()
{
    const int NR_GRAPHS = 2;
    Graph g[NR_GRAPHS];

    mk_g_1(g[0]);
    mk_g_2(g[1]);

    FSGLattice lattice(g, g+1);
    lattice.grviz_write();

    cout << "Num graphs = " << lattice.get_num_graphs() << endl;
    cout << "Num all graphs = " << lattice.get_num_all_created_graphs() << endl;
}
