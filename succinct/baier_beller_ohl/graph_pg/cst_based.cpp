#include <algorithm>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <assert.h>
#include <list>
#include <unistd.h>

#include <sdsl/construct.hpp>
#include <sdsl/cst_sada.hpp>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/bit_vectors.hpp>


using namespace std;
using namespace sdsl;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


struct cdbg_expl2
{
	private:
		struct node
		{
			uint64_t len;
			std::vector<uint64_t> adj_list;
			std::vector<uint64_t> pos_list;
			node()
			{
				len = 0;
			}
		};

		std::vector<node> G;
		uint64_t k;
	public:
		typedef uint64_t size_type;
		typedef uint64_t node_type;
		cdbg_expl2(uint64_t _k, uint64_t n=0)
		{
			k = _k;
			G = std::vector<node>(n);
		}

		void resize(uint64_t n)
		{
			G.resize(n);
		}

		void add_node(uint64_t l, uint64_t pos)
		{
			node n;
			n.len = l;
			n.pos_list.emplace_back(pos);
			G.emplace_back(n);
		}

		void add_node(uint64_t l, uint64_t pos, uint64_t edge)
		{
			node n;
			n.len = l;
			n.adj_list.emplace_back(edge);
			n.pos_list.emplace_back(pos);
			G.emplace_back(n);
		}

		//! Order of the debruijn graph
		size_type order() const {
			return k;
		}

		//! Number of nodes
		size_type size() const {
			return G.size();
		}

		//! Indicates if the graph is empty
		bool empty() const {
			return size() == 0;
		}

		//! Number of edges of node n
		uint64_t edges(node_type n) const {
			assert( n < size() );
			return G[n].adj_list.size();
		}

		//! k-th edge of node n
		uint64_t edge(node_type n, size_type k) const {
			assert( k < edges(n) );
			return G[n].adj_list[k];
		}

		//! Number of suffixes of node n
		uint64_t suffixes(node_type n) const {
			assert( n < size() );
			return G[n].pos_list.size();
		}

		//! k-th suffix of node n
		uint64_t suffix(node_type n, size_type k) const {
			assert( k < suffixes(n) );
			return G[n].pos_list[k];
		}

		//! length of node n
		size_type length(node_type n) const {
			assert(n < size());
			return G[n].len;
		}

		void set_length(node_type n, size_type l) {
			assert(n < size());
			G[n].len = l;
		}
		void insert_edge(node_type n, size_type e) {
			assert(n < size());
			G[n].adj_list.emplace_back(e);
		}
		void insert_suffix(node_type n, size_type s) {
			assert(n < size());
			G[n].pos_list.emplace_back(s);
		}
		size_type remove_edge(node_type n) {
			assert(n < size());
			size_type res = G[n].adj_list.back();
			G[n].adj_list.pop_back();
			return res;
		}
		size_type remove_suffix(node_type n) {
			assert(n < size());
			size_type res = G[n].pos_list.back();
			G[n].pos_list.pop_back();
			return res;
		}
};


//! computation of repeat nodes, where all left-maximal nodes and all nodes with depth < k get marked.
//! While nodes get marked, all left-maximal nodes whose children all are leafs are stored to a list,
//! so the cst need not to be iterated fully for second step.
/**
 * \param cst the suffix tree used for computation
 * \param k order of debruijn-graph, for which to construct repeat nodes
 * \param create_rep_node a function expecting the cst, a node and the length of the repnode
 *                        which is able to create a repeat node. Note that the string depth
 *                        of the given cst node may be larger than the given length, so the
 *                        tree may must be iterated upwards. Also note that the function must
 *                        handle duplicates, i.e., the same node could be constructed twice.
 * Used Suffix tree should use a function id, what maps a node to a number in
 * range [0...#nodecount). Also, this function should map leafs to the range [0..#leafcount),
 * and internal nodes to the range [#leafcount..#nodecount).
 *
 */
template<class CST,
         typename node_type = typename CST::node_type,
         typename size_type = typename CST::size_type>
struct rep_nodes_algorithm_3 {

	template<class CreateRepNode>
	static void run( const CST &cst, size_type k, CreateRepNode &&create_rep_node) {
		assert(k > 0);
		if (cst.size() == 0)	return;
		assert(k <= cst.size());

		//step 1: mark all left-maximal internal nodes or internal nodes with sdepth < k.
		//Also, set up a bit vector so LAQs for fixed string depth k can be computed fast,
		//and store nodes to be iterated in second step
		std::list<node_type> nodelist;
		sdsl::bit_vector marked( cst.nodes() - cst.size() );
		sdsl::bit_vector laqs_k_bv( cst.size() ); //number of leafs
		for (auto it = cst.begin_bottom_up(); it != cst.end_bottom_up(); ++it) {
			node_type v = *it;
			if (!cst.is_leaf(v)) {
				bool markv = (cst.depth(v) < k);
				if (!markv) {
					//check all children, if they have different beginning characters
					// or are marked already
					auto children = cst.children(v);
					auto child_it = children.begin();
					char c = (char)cst.csa.bwt[cst.lb(*child_it)];
					bool onlyleafs = true;
					do {
						auto u = *child_it;
						onlyleafs &= cst.is_leaf(u);
						//mark node if child node is marked or bwt characters differ
						markv |= (!cst.is_leaf(u) && marked[cst.id(u)-cst.size()])
						      || (c != cst.csa.bwt[cst.lb(u)]);

					} while (!markv && ++child_it != children.end());

					//second, if the actual node is LAQ of size k, mark it's interval in
					// the laqs_k bitvector
					if (cst.depth(cst.parent(v)) < k) {
						laqs_k_bv[cst.lb(v)] = true;
						laqs_k_bv[cst.rb(v)] = true;
					}

					//finally, if the node is left maximal and has only leafs as
					//children, store it to the list
					if (markv && onlyleafs) {
						//check if the rest of the children are leafs
						do {
							if (child_it == children.end()) {
								nodelist.push_back(v);
								break;
							}
							onlyleafs &= cst.is_leaf(*child_it);
							++child_it;
						} while (onlyleafs);
					}
				}
				//save resulting mark
				marked[cst.id(v)-cst.size()] = markv;
			}
		}

		//init rank and select support for LAQs with fixed depth k
		sdsl::bit_vector::select_1_type laqs_k_sel(&laqs_k_bv);
		sdsl::bit_vector::rank_1_type   laqs_k_rank(&laqs_k_bv);

		//step 2: create repeat nodes
		//iterate the nodes of list created in first step
		while (!nodelist.empty()) {
			node_type v = nodelist.back();
			nodelist.pop_back();

			node_type cur = v;
			size_type l = 0;
			do {
				//compute u=LAQs(v,k) efficient using rank/select supports
				//compute lb and rb of node u
				size_type r = laqs_k_rank(cst.rb(v));
				size_type lb = laqs_k_sel(r);
				size_type rb = laqs_k_sel(r+1);
				node_type u = cst.node(lb, rb); //use LCA to get u, using lb and rb

				v = cst.sl(v);
				if (marked[cst.id(u)-cst.size()]) {
					if (l > 0) {
						create_rep_node(cst,cur,l+k-1);
					}
					if (cst.depth(u) == k) {
						create_rep_node(cst,u,k);
						cur = v;
						l = 0;
					} else { //depth(u) != k
						cur = u;
						l = 1;
					}
				} else { //u is not left-maximal
					if (cst.depth(u) == k) {
						create_rep_node(cst,cur,l+k);
						cur = v;
						l = 0;
					} else { //depth(u) != k
						++l;
					}
				}
			} while (!marked[cst.id(v)-cst.size()]);
			if (l > 0) {
				create_rep_node(cst,cur,l+k-1);
			}
		}
	}
};


template<class CST,
         typename node_type = typename CST::node_type,
         typename size_type = typename CST::size_type>
using rep_nodes_algorithm = rep_nodes_algorithm_3<CST, node_type, size_type>;

template<class CST,
         typename node_type = typename CST::node_type,
         typename size_type = typename CST::size_type>
cdbg_expl2 myFunc(CST &cst, size_type k, sdsl::cache_config& config) {
	cdbg_expl2 G(k);
	uint64_t n = cst.size();
	{
		//compute repeat nodes and store them to a list
		sdsl::bit_vector B( cst.size() );
		sdsl::rank_support_v<01,2> B_rank;
		std::list<std::pair<node_type,size_type>> L;
		auto create_rep_node = [&B, &L]
			(const CST &cst,node_type v,size_type l) {
				//iterate the tree up until u = LAQs(v,l), and mark (lb(u)..rb(u)] in bitvector.
				//to ensure no duplicate nodes are created, check if rb is set before
				//iterating up the tree
				if (B[cst.rb(v)])	return; //node exists already
				for (node_type u = cst.parent(v); cst.depth(u) >= l; u = cst.parent(u)) {
					v = u;
				}
				//now, v is the searched internal node, so mark (lb(v)..rb(v)],
				// and store node and node length to the list.
				for (size_type i = cst.lb(v) + 1; i <= cst.rb(v); ++i) {
					B[i] = true;
				}
				L.push_back( std::make_pair( v, l ) );
			};
		rep_nodes_algorithm<CST,node_type,size_type>::run(cst, k, create_rep_node);
		sdsl::util::init_support( B_rank, &B );

		size_type counter = L.size();
		G.resize(counter);

		for (const auto& element : L) {
			uint64_t lb = cst.lb( element.first );
			uint64_t rb = cst.rb( element.first );
			size_type id = B_rank( lb );
			size_type l = element.second;
			G.set_length(id, l);
			G.insert_edge(id, lb);
			G.insert_suffix(id, rb);
		}
		sdsl::util::clear(cst);
	}

	// Fill Graph
	std::string filename_a = "a_array.sdsl";
	{
		uint64_t undef = n+1;
		uint8_t int_width = sdsl::bits::hi(undef+1)+1;
		{
			sdsl::int_vector<> A(n, undef, int_width);
			{
				sdsl::int_vector_buffer<> sa_buff(sdsl::cache_file_name(sdsl::conf::KEY_SA, config));
				for(uint64_t i=0; i<G.size(); ++i)
				{
					uint64_t lb = G.remove_edge(i);
					uint64_t rb = G.remove_suffix(i);
					for(uint64_t j=lb; j<=rb; ++j)
					{
						A[sa_buff[j]] = i;
					}
				}
			}
			sdsl::store_to_file(A, filename_a);
		}

		sdsl::int_vector_buffer<> A_buf(filename_a);
		uint64_t first_pos = 0;
		while(A_buf[first_pos] == undef)
		{
			++first_pos;
		}
		G.insert_suffix(A_buf[first_pos], first_pos);
		if(first_pos != 0)
		{
			G.add_node(first_pos-1+k, 0, A_buf[first_pos]);
		}

		uint64_t prev_pos = first_pos;
		uint64_t prev_id  = A_buf[first_pos];
		for(uint64_t pos = first_pos+1; pos < n; ++pos)
		{
			uint64_t id = A_buf[pos];
			if(id != undef)
			{
				uint64_t new_pos = prev_pos + G.length(prev_id) - (k-1);
				if(pos == new_pos)
				{
					G.insert_edge(prev_id, id);
				}
				else
				{
					G.insert_edge(prev_id, G.size());
					G.add_node(pos-new_pos+k-1, new_pos, id);
				}
				G.insert_suffix(id, pos);
				prev_pos = pos;
				prev_id  = id;
			}
		}
		G.insert_edge(prev_id, G.size());
		uint64_t new_pos = prev_pos + G.length(prev_id) - (k-1);
		G.add_node(n - new_pos, new_pos);
		A_buf.close(true);
	}
	return G;
}




void printUsage( char *command ) {
	cerr << command << " k INFILE OUTFILE" << endl;
}

void unregister_cache_file(const std::string& key, cache_config& config)
{
    config.file_map.erase(key);
}

typedef sdsl::cst_sct3<> CST;

int main( int argc, char **argv) {
	typedef typename CST::size_type size_type;

	size_type k;
	construct_config::byte_algo_sa = SE_SAIS; // or LIBDIVSUFSORT for less space-efficient but faster construction

	CST cst;

	if ( argc != 4 ) {
		cerr << "Illegal number of parameters" << endl;
		printUsage( argv[0] );
		return 1;
	}
	//check if file exists
	if (access(argv[2], R_OK) < 0) {
		cerr << "File " << argv[2] << " cannot be accessed" << endl;
		printUsage( argv[0] );
		return 1;
	}
	// Transform FASTA-Format
	cache_config config(false, ".", "tmp");
	string filename = sdsl::tmp_file(config);
	{
		string input = argv[2];
		filename += ".tmp";
		int_vector_buffer<8> input_buf(argv[2], std::ios::in, 1024*1024, 8, true);
		int_vector_buffer<8> output_buf(filename, std::ios::out, 1024*1024, 8, true);
		if(!output_buf.good())
		{
			cerr << "File " << output_buf.filename() << " could not be written" << endl;
			return 1;
		}
		for(uint64_t i=0, target=0, skip=0, seperator=1; i<input_buf.size(); ++i)
		{
			uint8_t input_char = input_buf[i];
			if(input_char == '>')
			{
				skip = 1;
				if(i>0)
				{
					output_buf[target] = seperator;
					++target;
				}
			}
			else if(input_char == '\n')
			{
				skip = 0;
			}
			else if(skip == 0)
			{
				output_buf[target] = input_char;
				++target;
			}
		}
	}

	//construct cst from file
	construct( cst, filename, config, 1 );
	unregister_cache_file(conf::KEY_SA, config);
        util::delete_all_files(config.file_map);
	sdsl::remove(filename);

	//fetch k
	int k_ = atoi(argv[1]);
	if (k_ < 1 || (size_type)k_ > cst.size()) {
		cerr << "order k must be greater 0 and smaller or equal to text size" << endl;
		cerr << "k=" << k << " (size_type)k_=" << ((size_type)k_) << " cst.size()=" << cst.size() << endl;
		printUsage( argv[0] );
		return 1;
	}
	k = (size_type)k_;

	auto cdbg = myFunc( cst, k, config );

	register_cache_file(conf::KEY_SA, config);
        util::delete_all_files(config.file_map);

	//and output it in .dot format to given file
	ofstream out( argv[3], ios_base::out | ios_base::trunc );
	if (!out.good()) {
		cerr << "Problem while opening file " << argv[3] << " for writing" << endl;
		printUsage( argv[0] );
		return 1;
	}
	//write the graph
	out << "digraph G {" << endl;
	for (size_type n = 0; n < cdbg.size(); n++) { //for each node
		out << "  " << n << " [label=\"" << (cdbg.suffix(n, 0)+1);
		//print suffixes
		for (size_type s = 1; s < cdbg.suffixes(n); s++) {
			out << "," << (cdbg.suffix(n, s)+1);
		}
		//print length of node
		out << ":" << cdbg.length(n) << "\"]" << endl;
		//second, print edges
		for (size_type e = 0; e < cdbg.edges(n); e++) {
			out << "  " << n << " -> " << cdbg.edge(n,e) << endl;
		}
	}
	out << "}";
	//check if everything went well
	out.close();
	if (!out.good()) {
		cerr << "Problems while writing to file " << argv[3] << endl;
		return 1;
	}
	return 0;
}
