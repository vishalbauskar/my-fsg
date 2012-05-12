#ifndef FSG_CANDIDATE_H_
#define FSG_CANDIDATE_H_

#include <list>
#include <map>

#include <boost/graph/adjacency_list.hpp>
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

#include "fsg_candidate_impl.hpp"

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

} // namespace graph_alg

#endif
