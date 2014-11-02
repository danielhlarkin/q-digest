#if !defined QDIGEST_H
#define QDIGEST_H

#include <iostream>
#include <vector>
#include <memory>
#include <cstddef>
#include <assert.h>

// #define DEBUGERR(X)
#define DEBUGERR(X) std::cerr << X

#if !defined(nullptr)
const                        // this is a const object...
class {
 public:
  template<class T>          // convertible to any type
    operator T*() const      // of null non-member
    { return 0; }            // pointer...
  template<class C, class T> // or any type of null
    operator T C::*() const  // member pointer...
  { return 0; }
 private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};              // and whose name is nullptr
#endif

namespace qdigest {

  size_t log2Ceil(size_t n) {
    bool is_pow2 = (n & -n) == n;
    int ret = 0;
    while (n > 1) { n /= 2; ++ret; }
    return ret + (is_pow2 ? 0 : 1);
  }

  struct QDigestNode {
    QDigestNode *left, *right, *parent;
    size_t count;
    size_t lb, ub; // Range is [lb..ub] (i.e. both inclusive)
  QDigestNode(size_t _lb, size_t _ub)
  : left(nullptr), right(nullptr), parent(nullptr),
      count(0),
      lb(_lb), ub(_ub) { }
    ~QDigestNode() {
      delete this->left;
      delete this->right;
      this->left = this->right = nullptr;
      this->parent = nullptr;
    }
  };

  class QDigest {
    // A QDigest can NOT be copied
    std::unique_ptr<QDigestNode> root;
    size_t num_nodes;
    size_t N, K;
    size_t num_inserts;

    size_t node_and_sibbling_count(QDigestNode *n) {
      size_t ret = n->count;
      if (n->left) { ret += n->left->count; }
      if (n->right) { ret += n->right->count; }
      return ret;
    }

    bool delete_node_if_needed(QDigestNode *n, int level, int l_max) {
      if (n->count == 0 && ((level == l_max) || (!n->left && !n->right))) {
        if (n->parent->left == n) {
          n->parent->left = nullptr;
        } else {
          n->parent->right = nullptr;
        }
        delete n;
        --num_nodes;
        return true;
      }
      return false;
    }

    void compact(QDigestNode *n, int level, int l_max, size_t nDivk) {
      if (!n) return;
      DEBUGERR("compact [" << n->lb << ", " << n->ub << "] "
               << "level: " << level << "\n");
      compact(n->left, level + 1, l_max, nDivk);
      compact(n->right, level + 1, l_max, nDivk);

      if (level > 0) {
        bool deleted = delete_node_if_needed(n, level, l_max);
        if (!deleted && node_and_sibbling_count(n->parent) < nDivk) {
          auto par = n->parent;
          par->count = node_and_sibbling_count(par);
          if (par->left) {
            par->left->count = 0;
            delete_node_if_needed(par->left, level, l_max);
          }
          if (par->right) {
            par->right->count = 0;
            delete_node_if_needed(par->right, level, l_max);
          }
        } // if (!deleted && ...)
      } // if (level > 0)
    }

    void printTree(std::ostream &out, QDigestNode *n) const {
      if (!n) return;
      printTree(out, n->left);
      printTree(out, n->right);
      out << "[" << n->lb << ".." << n->ub << "] --> " << n->count << "\n";
    }

    void expandTree(size_t ub) {
      DEBUGERR("expandTree(" << ub << ")\n");
      assert((ub & (-ub)) == ub); // ub should be a power of 2
      --ub;
      QDigest tmp(this->K, ub);
      const bool try_compact = false;
      tmp._insert(this->root->ub, 1, try_compact);

      QDigestNode *n = tmp.root.get();
      while (n->ub != this->root->ub) {
        DEBUGERR("UB: " << ub << std::endl);
        n = n->left;
      }
      QDigestNode *par = n->parent;
      int to_remove = 0;
      while (n) {
        DEBUGERR("n (lb, ub): (" << n->lb << ", " << n->ub << ")\n");
        n = n->right; ++to_remove;
      }
      par->left = this->root.release();
      par->left->parent = par;
      DEBUGERR("to_remove: " << to_remove << "\n");
      tmp.num_nodes -= to_remove;
      tmp.num_nodes += this->num_nodes;
      tmp.N = this->N;

      this->swap(tmp);
    }

    void _insert(size_t key, unsigned int count, bool try_compact) {
      if (key > this->root->ub) {
        expandTree(1 << log2Ceil(key));
      }
      size_t lb = 0;
      size_t ub = this->root->ub;
      QDigestNode *prev = this->root.get();
      QDigestNode *curr = prev;
      while (lb != ub) {
        size_t mid = lb + (ub - lb) / 2;
        prev = curr;
        if (key <= mid) {
          // Go left
          if (!curr->left) {
            prev->left = new QDigestNode(lb, mid);
            prev->left->parent = prev;
            ++this->num_nodes;
          }
          curr = prev->left;
          ub = mid;
        } else {
          // Go right
          assert(mid + 1 <= ub);
          if (!curr->right) {
            prev->right = new QDigestNode(mid + 1, ub);
            prev->right->parent = prev;
            ++this->num_nodes;
          }
          curr = prev->right;
          lb = mid + 1;
        }
      } // while()
      curr->count += count;
      this->N += count;
      if (try_compact && (this->num_nodes >= K * 6)) {
        const size_t nDivk = (N / K) + (N % K ? 1 : 0);
        const int l_max = log2Ceil(this->root->ub + 1);
        DEBUGERR("nDivk: " << nDivk << ", l_max: " << l_max << "\n");
        this->compact(this->root.get(), 0, l_max, nDivk);
      }
    }

    size_t postorder_by_rank(QDigestNode *n,
                             size_t &curr_rank,
                             size_t req_rank) const {
      if (!n) return 0;
      size_t val = postorder_by_rank(n->left, curr_rank, req_rank);
      if (curr_rank >= req_rank) return val;
      val = postorder_by_rank(n->right, curr_rank, req_rank);
      if (curr_rank >= req_rank) return val;

      val = n->ub;
      curr_rank += n->count;
      return val;
    }

  public:
    explicit QDigest(size_t _k, size_t ub = 1)
      : root(new QDigestNode(0, ub)),
      num_nodes(1),
      N(0), K(_k)
    { }

    void swap(QDigest &other) {
      std::swap(this->root, other.root);
      std::swap(this->num_nodes, other.num_nodes);
      std::swap(this->N, other.N);
      std::swap(this->K, other.K);
    }

    void insert(size_t key, unsigned int count) {
      const bool try_compact = true;
      this->_insert(key, count, try_compact);
    }

    size_t percentile(double p) const {
      // p is in the range [0..1]
      size_t curr_rank = 0;
      const size_t req_rank = p * N;
      return postorder_by_rank(this->root.get(),
                               curr_rank,
                               req_rank);
    }

    void printTree(std::ostream &out) const {
      out << "num_nodes: " << num_nodes << ", (N, K): ("
          << N << ", " << K << ")\n";
      this->printTree(out, this->root.get());
    }
  };
} // namespace qdigest

#endif // QDIGEST_H
