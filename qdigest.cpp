#include "qdigest.h"

using namespace qdigest;

size_t qdigest::log2Ceil(size_t n) {
  bool is_pow2 = (n & -n) == n;
  int ret = 0;
  while (n > 1) {
    n /= 2;
    ++ret;
  }
  return ret + (is_pow2 ? 0 : 1);
}

QDigestNode::QDigestNode(size_t _lb, size_t _ub)
    : left(nullptr), right(nullptr), parent(nullptr), count(0), lb(_lb),
      ub(_ub) {}

QDigestNode::~QDigestNode() {
  delete this->left;
  delete this->right;
  this->left = this->right = nullptr;
  this->parent = nullptr;
}

std::ostream &qdigest::operator<<(std::ostream &out, QDigestNode const &n) {
  out << "[" << n.lb << ".." << n.ub << "], count -> " << n.count;
  return out;
}

size_t QDigest::node_and_sibbling_count(QDigestNode *n) {
  size_t ret = n->count;
  if (n->left) {
    ret += n->left->count;
  }
  if (n->right) {
    ret += n->right->count;
  }
  return ret;
}

bool QDigest::delete_node_if_needed(QDigestNode *n, int level, int l_max) {
  if (n->count == 0 && (!n->left && !n->right)) {
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

void QDigest::compact(QDigestNode *n, int level, int l_max, size_t nDivk) {
  if (!n)
    return;
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
  }   // if (level > 0)
}

void QDigest::printTree(std::ostream &out, QDigestNode *n) const {
  if (!n)
    return;
  printTree(out, n->left);
  printTree(out, n->right);
  out << *n << "\n";
}

void QDigest::expandTree(size_t ub) {
  DEBUGERR("expandTree(" << ub << ")\n");
  assert(ub - 1 > this->root->ub);
  assert((ub & (-ub)) == ub); // ub should be a power of 2
  --ub;
  QDigest tmp(this->K, ub);

  if (this->N == 0) {
    this->swap(tmp);
    return;
  }

  const bool try_compact = false;
  tmp._insert(this->root->ub, 1, try_compact);

  // Intuitively, we keep going down the left child till we find
  // that no left child exists. This is the point where we
  // branched off to the right, and we need to trim the tree below
  // this node and replace this node with the root node of the
  // original tree.

  QDigestNode *n = tmp.root.get();
  while (n->ub != this->root->ub) {
    DEBUGERR("UB: " << ub << std::endl);
    n = n->left;
  }
  QDigestNode *par = n->parent;
  int to_remove = 0;
  while (n) {
    DEBUGERR("node --> (lb, ub): (" << n->lb << ", " << n->ub << ")\n");
    n = n->right;
    ++to_remove;
  }
  par->left = (QDigestNode *)0x1;
  par->left = this->root.release();
  par->left->parent = par;
  DEBUGERR("to_remove: " << to_remove << "\n");
  tmp.num_nodes -= to_remove;
  tmp.num_nodes += this->num_nodes;
  tmp.N = this->N;

  this->swap(tmp);
}

void QDigest::_insert_node(const QDigestNode *n) {
  DEBUGERR("_insert_node (" << *n << ")\n");
  auto r = this->root.get();
  assert(n->lb >= r->lb);
  assert(n->ub <= r->ub);

  QDigestNode *prev = this->root.get();
  QDigestNode *curr = prev;
  while (curr->lb != n->lb || n->ub != curr->ub) {
    size_t mid = curr->lb + (curr->ub - curr->lb) / 2;
    DEBUGERR("mid: " << mid << "\n");
    prev = curr;
    if (n->ub <= mid) {
      // Go left
      if (!prev->left) {
        prev->left = new QDigestNode(curr->lb, mid);
        prev->left->parent = prev;
        ++this->num_nodes;
      }
      curr = prev->left;
    } else {
      // Go right
      assert(mid + 1 <= curr->ub);
      if (!prev->right) {
        prev->right = new QDigestNode(mid + 1, curr->ub);
        prev->right->parent = prev;
        ++this->num_nodes;
      }
      curr = prev->right;
    }
  } // while()
  assert(curr->lb == n->lb);

  // curr should get the contents of 'n'
  curr->count += n->count;
  this->N += n->count;
  DEBUGERR("(curr, n): (" << *curr << ", " << *n << ")\n");
}

void QDigest::_insert(size_t key, unsigned int count, bool try_compact) {
  if (key > this->root->ub) {
    size_t new_ub_plus_one = ((size_t)1) << log2Ceil(key);
    if (this->root->ub + 1 == new_ub_plus_one) {
      new_ub_plus_one *= 2;
    }
    expandTree(new_ub_plus_one);
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
  if (try_compact) {
    compact_if_needed();
  }
}

void QDigest::compact_if_needed() {
  if (this->num_nodes >= K * 6) {
    const size_t nDivk = (N / K); // + (N % K ? 1 : 0);
    const int l_max = log2Ceil(this->root->ub + 1);
    DEBUGERR("nDivk: " << nDivk << ", l_max: " << l_max << "\n");
    this->compact(this->root.get(), 0, l_max, nDivk);
  }
}

size_t QDigest::postorder_by_rank(QDigestNode *n, size_t &curr_rank,
                                  size_t req_rank) const {
  if (!n)
    return 0;
  size_t val = postorder_by_rank(n->left, curr_rank, req_rank);
  if (curr_rank >= req_rank)
    return val;
  val = postorder_by_rank(n->right, curr_rank, req_rank);
  if (curr_rank >= req_rank)
    return val;

  val = n->ub;
  curr_rank += n->count;
  return val;
}

void QDigest::preorder_toString(QDigestNode *n, std::ostream &out) const {
  if (!n)
    return;
  if (n->count > 0) {
    out << n->lb << " " << n->ub << " " << n->count << "\n";
  }
  preorder_toString(n->left, out);
  preorder_toString(n->right, out);
}

QDigest::QDigest(size_t _k, size_t ub)
    : root(new QDigestNode(0, ub)), num_nodes(1), N(0), K(_k) {}

QDigest::QDigest(QDigest &&rhs) { this->swap(rhs); }

void QDigest::swap(QDigest &other) {
  std::swap(this->root, other.root);
  std::swap(this->num_nodes, other.num_nodes);
  std::swap(this->N, other.N);
  std::swap(this->K, other.K);
}

void QDigest::insert(size_t key, unsigned int count) {
  const bool try_compact = true;
  this->_insert(key, count, try_compact);
}

size_t QDigest::size() const { return this->N; }

bool QDigest::empty() const { return this->size() == 0; }

double QDigest::compression_ratio() const {
  return double(num_nodes) / (double)N;
}

size_t QDigest::percentile(double p) const {
  // p is in the range [0..1]
  size_t curr_rank = 0;
  const size_t req_rank = p * N;
  return postorder_by_rank(this->root.get(), curr_rank, req_rank);
}

std::string QDigest::toString() const {
  std::ostringstream sout;
  auto r = this->root.get();
  sout << N << " " << K << " " << r->lb << " " << r->ub << "\n";
  preorder_toString(this->root.get(), sout);
  return sout.str();
}

void QDigest::fromString(std::string const &ser) {
  std::istringstream sin(ser);
  int _n, _k, _lb, _ub;
  sin >> _n >> _k >> _lb >> _ub;
  this->K = _k;
  this->root.reset(new QDigestNode(_lb, _ub));

  while (true) {
    size_t lb, ub, count;
    sin >> lb >> ub >> count;
    if (!sin)
      break;
    QDigestNode n(lb, ub);
    n.count = count;
    DEBUGERR("Read QDigestNode: " << n << "\n");
    this->_insert_node(&n);
  }
}

void QDigest::merge(QDigest const &rhs) {
  const size_t max_k = std::max(this->K, rhs.K);
  const size_t max_ub = std::max(this->root->ub, rhs.root->ub);
  QDigest tmp(max_k, max_ub);

  std::queue<const QDigestNode *> nodesq;
  nodesq.push(this->root.get());
  nodesq.push(rhs.root.get());
  while (!nodesq.empty()) {
    auto n = nodesq.front();
    nodesq.pop();
    if (n->left) {
      nodesq.push(n->left);
    }
    if (n->right) {
      nodesq.push(n->right);
    }
    tmp._insert_node(n);
  }
  tmp.compact_if_needed();
  this->swap(tmp);
}

QDigest &QDigest::operator+=(QDigest const &rhs) {
  this->merge(rhs);
  return *this;
}

void QDigest::printTree(std::ostream &out) const {
  out << "[TREE] num_nodes: " << num_nodes << ", (N, K): (" << N << ", " << K
      << ")\n";
  this->printTree(out, this->root.get());
}

std::ostream &qdigest::operator<<(std::ostream &out, QDigest const &digest) {
  digest.printTree(out);
  return out;
}
