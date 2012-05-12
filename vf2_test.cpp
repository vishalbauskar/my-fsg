
#include "fsg.hpp"
#include "vf2.hpp"

#include "vf2_ori.hpp"

#include <vector>

using namespace boost;
using namespace std;

using namespace graph_alg;
using namespace graph_alg::detail;

typedef boost::adjacency_list<vecS,vecS, undirectedS, char> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor V;
typedef boost::graph_traits<Graph>::edge_descriptor E;

bool fv(V, V) { return true; }
bool fe(E, E) { return true; }


//typedef std::vector<MPair> OneMPairSec;
//typedef std::vector<OneMPairSec> TwoMPairSec;
/*
void prt(const OneMPairSec& c) {
    for (unsigned int i = 0; i < c.size(); i++)
	cout << "("<<c[i].first << ", " << c[i].second << ") ";
    cout << endl;
}

void prt(const TwoMPairSec& c) {
    for (unsigned int i = 0; i < c.size(); i++)
	prt(c[i]);
    cout << endl;
}
*/

struct VertexCompatible
{
    const Graph *g1, *g2;
    typedef boost::graph_traits<Graph>::vertex_descriptor V;
    VertexCompatible(const Graph &g1, const Graph &g2) :g1(&g1), g2(&g2) {}
    bool operator() (V v1, V v2) const {
	return (*g1)[v1] == (*g2)[v2];
    }

    bool operator() (V v1, V v2,
		     const Graph& g1, const Graph& g2) const {
	return (g1)[v1] == (g2)[v2];
    }

    V operator() (Graph& g, const Graph& g_src, V v) const
	{
	    V u = add_vertex(g);
	    g[u] = g_src[v];
	    return u;
	}

    // -----------------------------------------------
    typedef char VertexLabel;
    typedef const VertexLabel* VertexLabelIterator;
    static VertexLabel labels[3];
    const VertexLabelIterator begin() const { return labels; }
    const VertexLabelIterator end() const { return labels+3; }


    // -----------------------------------------------
    bool operator() (V v, const VertexLabel& vl, const Graph& g) const
	{ return g[v] == vl; }
    
};

VertexCompatible::VertexLabel
VertexCompatible::labels[3] = {'a', 'b', 'c'};


/*
template <typename T>
struct Predicate
{
    Predicate() {}
    Predicate(const T& x) :v(x) {}
    bool operator() (const T& x) const { return v != x; }
    T v;
};
*/
 //typedef Predicate<E> EdgePred;
 //typedef Predicate<V> VertexPred;
 //typedef boost::filtered_graph<Graph, EdgePred, VertexPred> FilGraph;

typedef std::vector<MPair> OneMPairSec;
typedef std::vector<OneMPairSec> TwoMPairSec;


void test_g(const Graph& g1, const Graph& g2, const VertexCompatible& vc)
{
    cout << "*** vf2 **** **************\n";
    TwoMPairSec vvec1;
    isomorphism_all(g1, g2, vvec1, vc);
    //prt(vvec1);


    cout << "*** vf2_ori ***************\n";
    TwoMPairSec vvec2;
    boost::vf2_all(g1, g2, vc, fe, vvec2);
    //prt(vvec2);
    
    if (vvec1 == vvec2)
	cout << "OK: Equal\n";
    else
	cout << "ERROR: not equal\n";
}

void mk_g(Graph& g)
{
/*
    g[0] = 'a';
    g[1] = 'b';
    g[2] = 'b';
    g[3] = 'b';
    g[4] = 'b';

    add_edge(0, 1, g);
    add_edge(1, 2, g);
    add_edge(2, 3, g);
    add_edge(3, 4, g);
    add_edge(4, 1, g);
*/

    g[0] = 'a';
    g[1] = 'a';
    g[2] = 'b';
    g[3] = 'c';

    add_edge(0, 1, g);
    add_edge(1, 2, g);
    add_edge(2, 3, g);
    //add_edge(3, 0, g);
}

void mk_g(Graph& g1, Graph& g2)
{
    const int N = 4;

    for (int i = 0; i < N; i++)
    {
	add_vertex(g1);
	add_vertex(g2);
    }

    // -------------------
    g1[0] = 'a';
    g1[1] = 'a';
    g1[2] = 'b';
    g1[3] = 'c';

    add_edge(0, 1, g1);
    add_edge(1, 2, g1);
    add_edge(2, 3, g1);
    add_edge(1, 3, g1);
    // -------------------

    g2[0] = 'a';
    g2[1] = 'c';
    g2[2] = 'a';
    g2[3] = 'b';

    add_edge(0, 1, g2);
    add_edge(1, 2, g2);
    add_edge(2, 3, g2);
    add_edge(1, 3, g2);
}

int main()
{
    const int N = 2;
    Graph g[N];

    mk_g(g[0], g[1]);
    VertexCompatible vc(g[0], g[1]);

/*
    FilGraph fg1(g1, EdgePred(e1), VertexPred(-111));
    FilGraph fg2(g2, EdgePred(e2), VertexPred(-111));
    TwoMPairSec isoms;
    bool iso_found = isomorphism_all(fg1, fg2, isoms, vc, 4, 4);
    if (iso_found)
    {
	//prt(fg1, fg2, e1, e2, -1, -1, isoms);
    }
*/

    //Lattice<Graph> lattice;
    //fsg(lattice, g, g+N, vc);
    
    //fsg_t(g[0], g[1], vc);
    //test_g(g1, g2, vc);
    return 0;
}
