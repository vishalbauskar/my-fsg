#ifndef PRINT_GRAPH_H_
#define PRINT_GRAPH_H_

#include <iostream>
#include <utility>
#include <string>

template<class G>
void print_graph(const G& g, const std::string& name = "")
{
    if (name != "")
	std::cout << name << std::endl;
    else
	std::cout << "****************************************" << std::endl;	
    typedef typename boost::graph_traits<G>::vertex_iterator VI;
    for (std::pair<VI,VI> vi = vertices(g); vi.first != vi.second; ++vi.first)
	std::cout << *vi.first << "{" << g[*vi.first] << "} ";
    std::cout << std::endl;

    typedef typename boost::graph_traits<G>::edge_iterator EI;
    for (std::pair<EI,EI> ei = edges(g); ei.first != ei.second; ++ei.first)
	std::cout << *ei.first << " ";
    std::cout << std::endl;
}


#endif
