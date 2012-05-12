#ifndef FSG_CANDIDATE_IMPL_H_
#define FSG_CANDIDATE_IMPL_H_


#include "print_graph.hpp"

template<class G, class L>
int Candidate<G,L>::N_G = 0;

template<class G, class L>
void debug_print(const Candidate<G,L>& g)
{
    using namespace std;

    print_graph(g, "model graph");
/*
    typedef typename Candidate<G,L>::OGS OGS;
    for (typename OGS::const_iterator i = g.ogs.begin(); i != g.ogs.end(); ++i)
    {
	//std::cerr << "origraph=" << i->first;
	//print_graph(*i->first, "  ");

	typedef typename Candidate<G,L>::M2OList::const_iterator Sb_CI;
	for (Sb_CI s = i->second.begin(); s != i->second.end(); ++s)
	{
	    cerr << "  subgraph: ";
	    cerr << "\n";
	}
    }

    cerr << "\n";
*/
}



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
//    debug_print(*this);
}


template<class G, class L>
Candidate<G,L>::Candidate(const Candidate& ancestor,
			  MV mv, const VertexLabel& lab)
    : ModelGraph(ancestor), lops(ancestor.lops)
{
    ID = ++N_G;
//    std::cerr << "*******************************************************************\n";
//    std::cerr << "Ctor2: mv=" << mv << " label=" << lab << "\n";
//    std::cerr << "\nAncestor:\n";
//    debug_print(ancestor);


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

    //  std::cerr << "Result:\n";
//    debug_print(*this);
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

    if (num_edges(*this) == 5 && 0)
    {
	std::cerr << "*******************************************************************\n";
	std::cerr << "Ctor3: mv1=" << mv1 << " mv2=" << mv2 << "\n";
	std::cerr << "\nAncestor:\n";
	debug_print(ancestor);


	std::cerr << "Result:\n";
	debug_print(*this);
    }
}


#endif
