#ifndef QDIGEST_H
#define QDIGEST_H

/**
 * Based on the Q-Digest algorithm documented in
 * http://www.cs.virginia.edu/~son/cs851/papers/ucsb.sensys04.pdf
 *
 * Builds a Q-Digest with a compression factor 'K' (parameter to the
 * constructor of the QDigest class). Higher values of 'K' result in
 * lesser compression, but more accuracy.
 *
 */

#include <assert.h>
#include <cstddef>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <vector>

#if defined DODEBUG
#define DEBUGERR(X) std::cerr << X
#else
#define DEBUGERR(X)
#endif

namespace qdigest {

extern size_t log2Ceil(size_t n);

struct QDigestNode {
  QDigestNode *left, *right, *parent;
  size_t count;
  size_t lb, ub; // Range is [lb..ub] (i.e. both inclusive)
  QDigestNode(size_t _lb, size_t _ub);
  ~QDigestNode();
};

std::ostream &operator<<(std::ostream &out, QDigestNode const &n);

class QDigest {
  // A QDigest can NOT be copied
  std::unique_ptr<QDigestNode> root;
  size_t num_nodes;
  size_t N, K;
  size_t num_inserts;

  /**
   * Returns the aggregated count for a node and its siblings.
   *
   */
  size_t node_and_sibbling_count(QDigestNode *n);

  /**
   * A tree node which has a count of 0 can be deleted only if it
   * has no children.
   *
   * Returns 'true' or 'false' depending on whether it deleted the
   * node 'n' from the tree.
   *
   */
  bool delete_node_if_needed(QDigestNode *n, int level, int l_max);

  /**
   * Perform compaction. Specifically, ensure that no node is too
   * small. i.e. apart from the root node, try to see if we can
   * compress counts from a set of 3 nodes (i.e. a node, its parent
   * and its sibling) and promote them to the parent node.
   *
   */
  void compact(QDigestNode *n, int level, int l_max, size_t nDivk);

  void printTree(std::ostream &out, QDigestNode *n) const;

  /**
   * Expand the range of the tree to include numbers in the range
   * [0..ub).
   *
   */
  void expandTree(size_t ub);

  /**
   * Insert the equivalent of the values present in node 'n' into
   * the current tree. This will either create new nodes along the
   * way and then create the final node or will update the count in
   * the destination node if that node is already present in the
   * tree. No compaction is attempted after the new node is inserted
   * since this function is assumed to be called by the
   * deserialization routine.
   *
   */
  void _insert_node(const QDigestNode *n);

  /**
   * Bump up the count for 'key' by 'count'.
   *
   * If 'try_compact' is 'true' then attempt compaction if
   * applicable. We don't compact when we want to build a tree which
   * has a specific shape since we assume that certain nodes will be
   * present at specific positions (for example when called by
   * expandTree().
   *
   */
  void _insert(size_t key, unsigned int count, bool try_compact);

  void compact_if_needed();

  /**
   * Perform a post-order traversal of the tree and fetch the
   * element at rank 'req_rank' starting from the smallest element
   * in the structure.
   *
   */
  size_t postorder_by_rank(QDigestNode *n, size_t &curr_rank,
                           size_t req_rank) const;

  /**
   * Perform a pre-prder traversal of the tree and serialize all the
   * nodes with a non-zero count. Separates each node with a newline
   * (\n).
   *
   */
  void preorder_toString(QDigestNode *n, std::ostream &out) const;

public:
  explicit QDigest(size_t _k, size_t ub = 1);
  QDigest(QDigest &&rhs);

  void swap(QDigest &other);

  void insert(size_t key, unsigned int count);

  size_t size() const;

  bool empty() const;

  double compression_ratio() const;

  /**
   * Returns the approximate 100p'th percentile element in the
   * structure. i.e. passing in 0.7 will return the 70th percentile
   * element (which is the 70th percentile element starting from the
   * smallest element).
   *
   */
  size_t percentile(double p) const;

  /**
   * Serialized format consists of newline separated entries which
   * are tripples of the form: (LB, UB, COUNT)
   *
   * That means that we have a node which has COUNT elements which
   * have values in the range [LB..UB]. Only non-empty ranges will
   * be serialized (i.e. the serialized tree will be sparse). Also,
   * the ranges will be serialized in pre-order tree traversal so
   * that re-construction is easy.
   *
   */
  std::string toString() const;

  /**
   * Deserialize the tree from the serialized version in the string
   * 'ser'. The serialized version is obtained by calling
   * toString().
   *
   */
  void fromString(std::string const &ser);

  void merge(QDigest const &rhs);

  QDigest &operator+=(QDigest const &rhs);

  void printTree(std::ostream &out) const;
};

extern std::ostream &operator<<(std::ostream &out, QDigest const &digest);

} // namespace qdigest

#endif // QDIGEST_H
