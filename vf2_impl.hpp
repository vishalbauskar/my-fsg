


// ===============================================================
//	class VertexEnumerator
// ===============================================================

inline bool VertexEnumerator::find_next()
{
    assert(0);
    if (v == NULL_VERTEX)
	v = 0;
    else
	++v;
    while (v < n)
    {
	if (is_unmapped_vertex(v))
	    return true;
	++v;
    }
    return false;
}


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

/*
    std::cerr<<"Add Pair: "<<p.first<<" " <<p.second<< std::endl;
    std::cerr<<"core_len = " << core_len << std::endl;
    std::cerr<<"out1_len = " << t1len.out << std::endl;
    std::cerr<<"out2_len = " << t2len.out << std::endl;
    std::cerr<<"in1_len  = " << t1len.in << std::endl;
    std::cerr<<"in2_len  = " << t2len.in << std::endl;
    prt();
*/
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

