#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

typedef std::pair<int, int> grid_point;
typedef std::pair<grid_point, grid_point> edge;

struct comp_t {
    bool operator()(const grid_point &lhs, const grid_point &rhs) const {
        if (lhs.second != rhs.second)
            return lhs.second < rhs.second;
        return lhs.first < rhs.first;
    }
};

struct Graph {
    std::set<grid_point> P_x_sorted;
    std::set<grid_point, comp_t> P_y_sorted;
    std::list<edge> E;
    bool grown = false;
    int L = 0;
};

std::vector<int> TerMapPermut(std::vector<grid_point> &terminals) {
    sort(terminals.begin(), terminals.end(),
        [](const auto &lhs, const auto &rhs) {
            if (lhs.second != rhs.second)
                return lhs.second < rhs.second;
            return lhs.first < rhs.first;
        }
    );

    std::map<grid_point, int> order;
    for (int i = 0; i < terminals.size(); ++i)
        order[terminals[i]] = i + 1;

    sort(terminals.begin(), terminals.end());

    std::vector<int> C(terminals.size());
    for (int i = 0; i < C.size(); ++i)
        C[i] = order[terminals[i]];

    return C;
}

Graph Permut(const std::vector<int> &C) {
    Graph P;
    for (int i = 0; i < C.size(); ++i) {
        P.P_x_sorted.insert(std::make_pair(C[i], i+1));
        P.P_y_sorted.insert(std::make_pair(C[i], i+1));
    }
    return P;
}

bool extreme(Graph &G) {
    if (G.P_x_sorted.size() == 2) {
        auto p1 = *G.P_x_sorted.begin(), p2 = *G.P_x_sorted.rbegin();
        G.E.push_back({p1, p2});
        G.P_x_sorted.erase(G.P_x_sorted.begin());
        G.P_y_sorted.erase(p1);
        G.L += abs(p1.first - p2.first) + abs(p2.second - p1.second);
        return false;
    }

    auto x_it1  = G.P_x_sorted.begin(),  x_it2  = G.P_x_sorted.begin();  ++x_it2;
    auto x_rit1 = G.P_x_sorted.rbegin(), x_rit2 = G.P_x_sorted.rbegin(); ++x_rit2;
    auto y_it1  = G.P_y_sorted.begin(),  y_it2  = G.P_y_sorted.begin();  ++y_it2;
    auto y_rit1 = G.P_y_sorted.rbegin(), y_rit2 = G.P_y_sorted.rbegin(); ++y_rit2;

    grid_point p, adj_p;
    if (y_it1->second != y_it2->second) {
        p = *y_it1;
        G.P_x_sorted.erase(p);
        G.P_y_sorted.erase(y_it1);
        adj_p = {p.first, p.second + 1};
    }
    else if (y_rit1->second != y_rit2->second) {
        p = *y_rit1;
        G.P_x_sorted.erase(p);
        G.P_y_sorted.erase(y_rit2.base());
        adj_p = {p.first, p.second - 1};
    }
    else if (x_rit1->first != x_rit2->first) {
        p = *x_rit1;
        G.P_x_sorted.erase(x_rit2.base());
        G.P_y_sorted.erase(p);
        adj_p = {p.first - 1, p.second};
    }
    else if (x_it1->first != x_it2->first) {
        p = *x_it1;
        G.P_x_sorted.erase(x_it1);
        G.P_y_sorted.erase(p);
        adj_p = {p.first + 1, p.second};
    }
    else return true;

    G.P_x_sorted.insert(adj_p);
    G.P_y_sorted.insert(adj_p);
    G.E.push_back({p, adj_p});
    ++G.L;

    return extreme(G);
}

bool equal_edge_lists(const std::list<edge> &lhs, const std::list<edge> &rhs) {
    if (lhs.size() != rhs.size())
        return false;

    for (auto itl = lhs.begin(), itr = rhs.begin(); itl != lhs.end(); ++itl, ++itr)
        if (*itl != *itr)
            return false;

    return true;
}

void fork(const Graph &G, std::list<Graph> &TreeList) {
    Graph G1, G2, G3;
    G1 = G2 = G3 = G;

    auto x_it1  = G.P_x_sorted.begin(),  x_it2  = G.P_x_sorted.begin();  ++x_it2;
    auto x_rit1 = G.P_x_sorted.rbegin(), x_rit2 = G.P_x_sorted.rbegin(); ++x_rit2;
    auto y_it1  = G.P_y_sorted.begin(),  y_it2  = G.P_y_sorted.begin();  ++y_it2;
    auto y_rit1 = G.P_y_sorted.rbegin(), y_rit2 = G.P_y_sorted.rbegin(); ++y_rit2;

    if (x_it1->first == x_it2->first) {
        auto eps_i_j = *x_it1, eps_i_k = *x_it2;
        auto eps_i_plus1_j = std::make_pair(x_it1->first + 1, x_it1->second);
        auto eps_i_plus1_k = std::make_pair(x_it2->first + 1, x_it2->second);

        G1.P_x_sorted.erase(eps_i_j); G1.P_x_sorted.erase(eps_i_k);
        G1.P_y_sorted.erase(eps_i_j); G1.P_y_sorted.erase(eps_i_k);
        G1.P_x_sorted.insert(eps_i_plus1_j); G1.P_x_sorted.insert(eps_i_plus1_k);
        G1.P_y_sorted.insert(eps_i_plus1_j); G1.P_y_sorted.insert(eps_i_plus1_k);
        G1.E.push_back({eps_i_j, eps_i_plus1_j});
        G1.E.push_back({eps_i_k, eps_i_plus1_k});
        G1.L += 2;
        TreeList.push_back(G1);

        G2.P_x_sorted.erase(eps_i_j); G2.P_x_sorted.erase(eps_i_k);
        G2.P_y_sorted.erase(eps_i_j); G2.P_y_sorted.erase(eps_i_k);
        G2.P_x_sorted.insert(eps_i_plus1_j);
        G2.P_y_sorted.insert(eps_i_plus1_j);
        G2.E.push_back(std::make_pair(eps_i_j, eps_i_plus1_j));
        G2.E.push_back(std::make_pair(eps_i_j, eps_i_k));
        G2.L += 1 + abs(eps_i_j.second - eps_i_k.second);
        TreeList.push_back(G2);

        G3.P_x_sorted.erase(eps_i_j); G3.P_x_sorted.erase(eps_i_k);
        G3.P_y_sorted.erase(eps_i_j); G3.P_y_sorted.erase(eps_i_k);
        G3.P_x_sorted.insert(eps_i_plus1_k);
        G3.P_y_sorted.insert(eps_i_plus1_k);
        G3.E.push_back(std::make_pair(eps_i_k, eps_i_plus1_k));
        G3.E.push_back(std::make_pair(eps_i_j, eps_i_k));
        G3.L += 1 + abs(eps_i_j.second - eps_i_k.second);
        TreeList.push_back(G3);

        for (auto it = TreeList.begin(); it != TreeList.end(); ++it)
            if (it->P_x_sorted == G.P_x_sorted && equal_edge_lists(it->E, G.E)) {
                TreeList.erase(it);
                break;
            }

        return;
    }

    if (x_rit1->first == x_it2->first) {
        auto eps_i_j = *x_rit2, eps_i_k = *x_it1;
        auto eps_i_minus1_j = std::make_pair(x_rit2->first - 1, x_it2->second);
        auto eps_i_minus1_k = std::make_pair(x_rit1->first - 1, x_it1->second);

        G1.P_x_sorted.erase(eps_i_j); G1.P_x_sorted.erase(eps_i_k);
        G1.P_y_sorted.erase(eps_i_j); G1.P_y_sorted.erase(eps_i_k);
        G1.P_x_sorted.insert(eps_i_minus1_j); G1.P_x_sorted.insert(eps_i_minus1_k);
        G1.P_y_sorted.insert(eps_i_minus1_j); G1.P_y_sorted.insert(eps_i_minus1_k);
        G1.E.push_back({eps_i_j, eps_i_minus1_j});
        G1.E.push_back({eps_i_k, eps_i_minus1_k});
        G1.L += 2;
        TreeList.push_back(G1);

        G2.P_x_sorted.erase(eps_i_j); G2.P_x_sorted.erase(eps_i_k);
        G2.P_y_sorted.erase(eps_i_j); G2.P_y_sorted.erase(eps_i_k);
        G2.P_x_sorted.insert(eps_i_minus1_j);
        G2.P_y_sorted.insert(eps_i_minus1_j);
        G2.E.push_back(std::make_pair(eps_i_j, eps_i_minus1_j));
        G2.E.push_back(std::make_pair(eps_i_j, eps_i_k));
        G2.L += 1 + abs(eps_i_j.second - eps_i_k.second);
        TreeList.push_back(G2);

        G3.P_x_sorted.erase(eps_i_j); G3.P_x_sorted.erase(eps_i_k);
        G3.P_y_sorted.erase(eps_i_j); G3.P_y_sorted.erase(eps_i_k);
        G3.P_x_sorted.insert(eps_i_minus1_k);
        G3.P_y_sorted.insert(eps_i_minus1_k);
        G3.E.push_back(std::make_pair(eps_i_k, eps_i_minus1_k));
        G3.E.push_back(std::make_pair(eps_i_j, eps_i_k));
        G3.L += 1 + abs(eps_i_j.second - eps_i_k.second);
        TreeList.push_back(G3);

        for (auto it = TreeList.begin(); it != TreeList.end(); ++it)
            if (it->P_x_sorted == G.P_x_sorted && equal_edge_lists(it->E, G.E)) {
                TreeList.erase(it);
                break;
            }

        return;
    }

    if (y_it1->second == y_it2->second) {
        auto eps_i_k = *y_it1, eps_j_k = *y_it2;
        auto eps_i_k_plus1 = std::make_pair(y_it1->first, y_it1->second + 1);
        auto eps_j_k_plus1 = std::make_pair(y_it2->first, y_it2->second + 1);

        G1.P_x_sorted.erase(eps_i_k); G1.P_x_sorted.erase(eps_j_k);
        G1.P_y_sorted.erase(eps_i_k); G1.P_y_sorted.erase(eps_j_k);
        G1.P_x_sorted.insert(eps_i_k_plus1); G1.P_x_sorted.insert(eps_j_k_plus1);
        G1.P_y_sorted.insert(eps_i_k_plus1); G1.P_y_sorted.insert(eps_j_k_plus1);
        G1.E.push_back({eps_i_k, eps_i_k_plus1});
        G1.E.push_back({eps_j_k, eps_j_k_plus1});
        G1.L += 2;
        TreeList.push_back(G1);

        G2.P_x_sorted.erase(eps_i_k); G2.P_x_sorted.erase(eps_j_k);
        G2.P_y_sorted.erase(eps_i_k); G2.P_y_sorted.erase(eps_j_k);
        G2.P_x_sorted.insert(eps_i_k_plus1);
        G2.P_y_sorted.insert(eps_i_k_plus1);
        G2.E.push_back(std::make_pair(eps_i_k, eps_i_k_plus1));
        G2.E.push_back(std::make_pair(eps_i_k, eps_j_k));
        G2.L += 1 + abs(eps_i_k.first - eps_j_k.first);
        TreeList.push_back(G2);

        G3.P_x_sorted.erase(eps_i_k); G3.P_x_sorted.erase(eps_j_k);
        G3.P_y_sorted.erase(eps_i_k); G3.P_y_sorted.erase(eps_j_k);
        G3.P_x_sorted.insert(eps_j_k_plus1);
        G3.P_y_sorted.insert(eps_j_k_plus1);
        G3.E.push_back(std::make_pair(eps_j_k, eps_j_k_plus1));
        G3.E.push_back(std::make_pair(eps_i_k, eps_j_k));
        G3.L += 1 + abs(eps_i_k.first - eps_j_k.first);
        TreeList.push_back(G3);

        for (auto it = TreeList.begin(); it != TreeList.end(); ++it)
            if (it->P_x_sorted == G.P_x_sorted && equal_edge_lists(it->E, G.E)) {
                TreeList.erase(it);
                break;
            }

        return;
    }

    if (y_rit1->first == y_rit2->first) {
        auto eps_i_k = *y_rit2, eps_j_k = *y_rit1;
        auto eps_i_k_minus1 = std::make_pair(y_rit2->first, y_rit2->second - 1);
        auto eps_j_k_minus1 = std::make_pair(y_rit1->first, y_rit1->second - 1);

        G1.P_x_sorted.erase(eps_i_k); G1.P_x_sorted.erase(eps_j_k);
        G1.P_y_sorted.erase(eps_i_k); G1.P_y_sorted.erase(eps_j_k);
        G1.P_x_sorted.insert(eps_i_k_minus1); G1.P_x_sorted.insert(eps_j_k_minus1);
        G1.P_y_sorted.insert(eps_i_k_minus1); G1.P_y_sorted.insert(eps_j_k_minus1);
        G1.E.push_back({eps_i_k, eps_i_k_minus1});
        G1.E.push_back({eps_j_k, eps_j_k_minus1});
        G1.L += 2;
        TreeList.push_back(G1);

        G2.P_x_sorted.erase(eps_i_k); G2.P_x_sorted.erase(eps_j_k);
        G2.P_y_sorted.erase(eps_i_k); G2.P_y_sorted.erase(eps_j_k);
        G2.P_x_sorted.insert(eps_i_k_minus1);
        G2.P_y_sorted.insert(eps_i_k_minus1);
        G2.E.push_back(std::make_pair(eps_i_k, eps_i_k_minus1));
        G2.E.push_back(std::make_pair(eps_i_k, eps_j_k));
        G2.L += 1 + abs(eps_i_k.first - eps_j_k.first);
        TreeList.push_back(G2);

        G3.P_x_sorted.erase(eps_i_k); G3.P_x_sorted.erase(eps_j_k);
        G3.P_y_sorted.erase(eps_i_k); G3.P_y_sorted.erase(eps_j_k);
        G3.P_x_sorted.insert(eps_j_k_minus1);
        G3.P_y_sorted.insert(eps_j_k_minus1);
        G3.E.push_back(std::make_pair(eps_j_k, eps_j_k_minus1));
        G3.E.push_back(std::make_pair(eps_i_k, eps_j_k));
        G3.L += 1 + abs(eps_i_k.first - eps_j_k.first);
        TreeList.push_back(G3);

        for (auto it = TreeList.begin(); it != TreeList.end(); ++it)
            if (it->P_x_sorted == G.P_x_sorted && equal_edge_lists(it->E, G.E)) {
                TreeList.erase(it);
                break;
            }

        return;
    }

    std::cout << "Impossible situation.\n";
    return;
}

Graph Const_optRST(Graph G) {
    std::list<Graph> TreeList;
    TreeList.emplace_back(G);

    auto it = TreeList.begin();
    while (it != TreeList.end()) {
        if (it->grown == false)
            if (extreme(*it)) {
                fork(*it, TreeList);
                it = TreeList.begin();
            }
            else {
                it->grown = true;
                ++it;
            }
        else
            ++it;
    }

    return *std::min_element(TreeList.begin(), TreeList.end(),
        [](const Graph &lhs, const Graph &rhs) {
            return lhs.L < rhs.L;
        }
    );
}

Graph RSTC(std::vector<grid_point> terminals, int district_size) {
    auto C = TerMapPermut(terminals);
    int max_size = std::min(district_size, 7);

    int n = C.size(), d_solved = 0;
    int size = std::min(n, max_size);

    auto itl = C.begin();
    auto itr = std::next(itl, size);
    auto last = std::next(C.end(), -1);

    Graph G_combined;

    std::vector<int> index(n);
    while (itl != last) {
        std::fill(index.begin(), index.end(), -1);
        auto it = itl;
        for (int i = 0; i < size; ++i, ++it)
            index[*it - 1] = i;

        std::vector<int> district(size);
        for (int i = 0, count = 1; i < n; ++i)
            if (index[i] >= 0) {
                district[index[i]] = count++;
            }
        auto G_d = Const_optRST(Permut(district));

        /* for (auto p : G_d.P_x_sorted) {
            auto p_copy = p;
            retrieve_permutation(p_copy, );
            G_combined.P_x_sorted.insert(p_copy);
        }
        for (auto p : G_d.P_y_sorted) {
            auto p_copy = p;
            retrieve_permutation(p_copy, );
            G_combined.P_y_sorted.insert(p_copy);
        }

        auto E_d = G_d.E;
        std::for_each(E_d.begin(), E_d.end(),
            [//](edge &e) {
                retrieve_permutation(e.first, );
                retrieve_permutation(e.second, );
            }
        );
        G_combined.E.merge(E_d,
            [](const auto &lhs, const auto &rhs) {
                return lhs.first == rhs.first  && lhs.second == rhs.second
                    || lhs.first == rhs.second && lhs.second == rhs.first;
            }
        ); */

        G_combined.L += G_d.L;

        ++d_solved;
        size = std::min(n - (d_solved * max_size), max_size);
        itl = std::next(itr, -1);
        std::advance(itr, size);
    }
    G_combined.grown = true;

    return G_combined;
}

int main() {
    std::vector<grid_point> terminals{{4, 3}, {6, 6}, {0, 47}, {2, 8}, {1, 6}, {3, 4}};
    int district_size = 7;

    Graph G = RSTC(terminals, district_size);

    for (auto it = G.P_x_sorted.begin(); it != G.P_x_sorted.end(); ++it)
        std::cout << "(" << it->first << ", " << it->second << ") ; ";
    std::cout << "\b\b \n";

    for (auto it = G.E.begin(); it != G.E.end(); ++it)
        std::cout << "(" << it->first.first  << ", " << it->first.second  << ") -> "
                  << "(" << it->second.first << ", " << it->second.second << ")\n";

    std::cout << G.L << std::endl;

    return 0;
}
