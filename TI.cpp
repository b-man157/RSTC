#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
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

std::vector<int> TerMapPermut(const std::vector<grid_point> &terminals) {
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

    auto x_it1  = G.P_x_sorted.begin(),  x_it2  = G.P_y_sorted.begin();  ++x_it2;
    auto x_rit1 = G.P_x_sorted.rbegin(), x_rit2 = G.P_y_sorted.rbegin(); ++x_rit2;
    auto y_it1  = G.P_y_sorted.begin(),  y_it2  = G.P_y_sorted.begin();  ++y_it2;
    auto y_rit1 = G.P_y_sorted.rbegin(), y_rit2 = G.P_y_sorted.rbegin(); ++y_rit2;

    grid_point p, adj_p;
    if (x_it1->first != x_it2->first) {
        p = *x_it1;
        G.P_x_sorted.erase(x_it1);
        G.P_y_sorted.erase(p);
        adj_p = {p.first + 1, p.second};
    }
    else if (x_rit1->first != x_rit2->first) {
        p = *x_rit1;
        G.P_x_sorted.erase(x_rit2.base());
        G.P_y_sorted.erase(p);
        adj_p = {p.first - 1, p.second};
    }
    else if (y_rit1->second != y_rit2->second) {
        p = *y_rit1;
        G.P_x_sorted.erase(p);
        G.P_y_sorted.erase(y_rit2.base());
        adj_p = {p.first, p.second - 1};
    }
    else if (y_it1->second != y_it2->second) {
        p = *y_it1;
        G.P_x_sorted.erase(p);
        G.P_y_sorted.erase(y_it1);
        adj_p = {p.first, p.second + 1};
    }
    else return true;

    G.P_x_sorted.insert(adj_p);
    G.P_y_sorted.insert(adj_p);
    G.E.push_back({p, adj_p});
    ++G.L;

    return extreme(G);
}

bool equal_lists(const std::list<edge> &lhs, const std::list<edge> &rhs) {
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

    auto x_it1  = G.P_x_sorted.begin(),  x_it2  = G.P_y_sorted.begin();  ++x_it2;
    auto x_rit1 = G.P_x_sorted.rbegin(), x_rit2 = G.P_y_sorted.rbegin(); ++x_rit2;
    auto y_it1  = G.P_y_sorted.begin(),  y_it2  = G.P_y_sorted.begin();  ++y_it2;
    auto y_rit1 = G.P_y_sorted.rbegin(), y_rit2 = G.P_y_sorted.rbegin(); ++y_rit2;

    if (x_it1->first == x_it2->first) {
        auto eps_i_j = *x_it1, eps_i_k = *x_it2;

        return;
    }

    if () {
        return;
    }

    if () {
        return;
    }

    if () {
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

int main() {
    std::vector<grid_point> terminals{{4, 3}, {6, 6}, {0, 47}, {2, 8}, {1, 6}, {3, 4}};
    auto C = TerMapPermut(terminals);

    for (auto i : C) {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    Graph G = Permut(C);

    for (auto it = G.P_x_sorted.begin(); it != G.P_x_sorted.end(); ++it)
        std::cout << "(" << it->first << ", " << it->second << ") ; ";
    std::cout << "\b\b \n";

    for (auto it = G.P_y_sorted.begin(); it != G.P_y_sorted.end(); ++it)
        std::cout << "(" << it->first << ", " << it->second << ") ; ";
    std::cout << "\b\b \n";

    Graph G_ret = Const_optRST(G);

    for (auto it = G_ret.P_x_sorted.begin(); it != G_ret.P_x_sorted.end(); ++it)
        std::cout << "(" << it->first << ", " << it->second << ") ; ";
    std::cout << "\b\b \n";

    for (auto it = G_ret.E.begin(); it != G_ret.E.end(); ++it)
        std::cout << "(" << it->first.first  << ", " << it->first.second  << ") -> "
                  << "(" << it->second.first << ", " << it->second.second << ")\n";


    std::cout << G_ret.L << std::endl;

    return 0;
}