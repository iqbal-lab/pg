#include <iostream>
#include <fstream>
#include <limits>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include "create_datastructures.hpp"
#include "partial_lcp.hpp"
#include "handle_graph.hpp"

using namespace std;
using namespace sdsl;
using namespace std::chrono;

struct node
{
	uint64_t len;
	vector<uint64_t> adj_list;
	vector<uint64_t> pos_list;
	node()
	{
		len = 0;
	}
};

tuple<vector<node>, vector<uint64_t>> create_cdbg_with_lf(cache_config& config, uint64_t k)
{
	// Create WT of the BWT
	typedef wt_huff<bit_vector, rank_support_v<>, select_support_scan<1>, select_support_scan<0>> wt;
	wt wt_bwt;
	construct(wt_bwt, cache_file_name(conf::KEY_BWT, config));

	// Create C-array (needed for interval_symbols)
	vector<uint64_t> carray(256, 0);
	for(uint64_t i=0, sum=0; i<256; ++i)
	{
		carray[i] = sum;
		sum += wt_bwt.rank(wt_bwt.size(), i);
	}

	// Create bit-vectors
	bit_vector bv(wt_bwt.size(), 0);
	bit_vector bv2(wt_bwt.size(), 0);
	bit_vector bv3(wt_bwt.size(), 0);
	{
		auto start = high_resolution_clock::now();
		int_vector<2> lcp_k = construct_partial_lcp<wt>(wt_bwt, carray, k);
		auto stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for partial LCP construction (sdsl-bblaca)" << endl;
		start = high_resolution_clock::now();
		int_vector_buffer<8> bwt(cache_file_name(conf::KEY_BWT, config));
		bool open=false;
		uint64_t kvalue=0;
		uint64_t lb=0;
		uint64_t last_change=0;
		vector<uint64_t> occ(256, 0);
		vector<uint8_t> occ_list;
		occ_list.reserve(256);
		vector<uint64_t> lf = carray;
		for(uint64_t i=1; i<lcp_k.size(); ++i)
		{
			++lf[bwt[i-1]];
			if(lcp_k[i] == gt_k or lcp_k[i] == eq_k)
			{
				open = true;
				if(lcp_k[i] == eq_k)
				{
					kvalue = i;
				}
			}
			else
			{
				if(open)
				{
					if(kvalue > lb)
					{
						bv[lb] = true;
						bv[i-1] = true;
					}
					if(last_change > lb)
					{
						for(uint64_t j=lb; j<=i-1; ++j)
						{
							bv2[j] = true;
							uint8_t c = bwt[j];
							if(occ[c]==0)
							{
								occ[c] = j;
								occ_list.emplace_back(c);
							}
						}
						for(uint64_t j=0; j<occ_list.size(); ++j)
						{
							bv3[lf[occ_list[j]]-1] = true;
							occ[occ_list[j]] = 0;
						}
						occ_list.resize(0);
					}
					open = false;
				}
				lb = i;
			}
			if(bwt[i] != bwt[i-1] or bwt[i] <= 1)
			{
				last_change = i;
			}
		}
		if(open)
		{
			++lf[bwt[lcp_k.size()-1]];
			if(kvalue > lb)
			{
				bv[lb] = true;
				bv[lcp_k.size()-1] = true;
			}
			if(last_change > lb)
			{
				for(uint64_t j=lb; j<=lcp_k.size()-1; ++j)
				{
					bv2[j] = true;
					uint8_t c = bwt[j];
					if(occ[c]==0)
					{
						occ[c] = j;
						occ_list.emplace_back(c);
					}
				}
				for(const auto& c : occ_list)
				{
					bv3[lf[c]-1] = true;
					occ[c] = 0;
				}
				occ_list.resize(0);
			}
		}
		bv3[0] = 1;
		for(uint64_t i=1; i<carray[2]; ++i)
		{
			bv3[i] = 0;
		}
		open = false;
		for(uint64_t i=0; i<bv.size(); ++i)
		{
			if(open)
			{
				bv3[i] = 0;
				if(bv[i])
				{
					open = false;
				}
			}
			else if(bv[i])
			{
				open = true;
				bv3[i] = 0;
			}
		}
		stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for create bv and initial nodes" << endl;
	}

	// Init rank support and graph
	bit_vector::rank_1_type bv_rank, bv3_rank;
	util::init_support(bv_rank, &bv);
	util::init_support(bv3_rank, &bv3);
	uint64_t number_right_max_nodes = bv_rank(bv.size())/2;
	vector<node> graph(number_right_max_nodes+bv3_rank(bv3.size()));
	vector<uint64_t> start_nodes;

	// Create compressed de bruijn graph
	{
		auto start = high_resolution_clock::now();
		uint64_t lb = 0;
		uint64_t cur_node = number_right_max_nodes;
		graph[cur_node].len = 1;
		for(uint64_t i=wt_bwt.size()-1; i>0; --i)
		{
			// LF
			auto res = wt_bwt.inverse_select(lb);
			uint8_t c = res.second;
			uint64_t lb_new = carray[c] + res.first;

			// Right maximal
			uint64_t ones = bv_rank(lb_new+1);
			uint64_t next_node = numeric_limits<uint64_t>::max();
			if(ones % 2 == 1 or bv[lb_new] == 1)
			{
				next_node = (ones-1)/2;
			}

			if(c <= 1) // c == sentinal
			{
				graph[cur_node].pos_list.emplace_back(i+1);
				start_nodes.emplace_back(cur_node);
				cur_node = graph.size();
				graph.emplace_back(node());
				graph[cur_node].len = 1;
			}
			else if(next_node != numeric_limits<uint64_t>::max() or bv2[lb]) // Next node is right max or cur node is left maximal => split
			{
				if(next_node == numeric_limits<uint64_t>::max())
				{
					next_node = number_right_max_nodes + bv3_rank(lb_new);
				}
				graph[cur_node].pos_list.emplace_back(i+1);
				graph[next_node].adj_list.emplace_back(cur_node);
				graph[next_node].len = k;
				cur_node = next_node;
			}
			else
			{
				++graph[cur_node].len;
			}
			lb = lb_new;
		}
		graph[cur_node].pos_list.emplace_back(1);
		start_nodes.emplace_back(cur_node);
		reverse(begin(start_nodes), end(start_nodes));

		auto stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for creating graph" << endl;
	}

	return make_tuple(move(graph), move(start_nodes));
}

int main(int argc, char *argv[])
{
	// Check parameters
	if(argc != 4)
	{
		cerr << "Usage: " << argv[0] << " inputfile outputfile kFile" << endl;
		return 1;
	}

	// Get parameters
	string inputfile = argv[1];
	string outputfile = argv[2];
	string kfilename = argv[3];

	// Create datastructures
	cache_config config(true, ".", "tmp");
	uint64_t errors = create_datastructures(config, inputfile, kfilename, false);

	// Read k-values
	ifstream kfile(kfilename);
	uint64_t k;
	while(!errors and kfile >> k)
	{
		// Create graph
		vector<node> graph;
		vector<uint64_t> start_nodes;
		tie(graph, start_nodes) = create_cdbg_with_lf(config, k);

		// Print graph
		{
			auto start = high_resolution_clock::now();
			ofstream output(outputfile+".k"+to_string(k)+".dot");
			ofstream output_start_nodes(outputfile+".k"+to_string(k)+".start_nodes.txt");
			print_graph(graph, start_nodes, output, output_start_nodes);
			auto stop = high_resolution_clock::now();
			cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for printing graph" << endl;
		}
	}

	// Delete files
	if(config.delete_files)
	{
		util::delete_all_files(config.file_map);
	}

	return errors;
}
