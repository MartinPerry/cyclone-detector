#ifndef BVH_H
#define BVH_H

#include <algorithm>

/// <summary>
/// Mostly based on:
/// http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies.html
/// </summary>
template <typename T, typename Aabb = d2::Aabb<MyMath::Vector2d>>
class BVH
{
protected:

	/// <summary>
	/// BVH tree node
	/// </summary>
	struct Node
	{
		int dataIndex; //index to orderedData array, -1 = internal node
		int dataCount; //number of items starting with [dataIndex], 0 = internal node
		Aabb aabb;	   //bounds of node

		int axis; //separating axis

		Node * left;
		Node * right;

		static Node * CreateLeaf(int index, int count,
			const Aabb & aabb)
		{
			return new Node(index, count, -1, aabb);
		}

		static Node * CreateInner(int axis, Node * left, Node * right)
		{
			Aabb aabb = Aabb::CreateUnion(left->aabb, right->aabb);
			Node * n = new Node(-1, 0, axis, aabb);
			n->left = left;
			n->right = right;
			return n;
		}


		Node(int index, int count, int axis, const Aabb & aabb) :
			dataIndex(index),
			dataCount(count),
			aabb(aabb),
			axis(axis),
			left(nullptr),
			right(nullptr)
		{
		};
	};

	/// <summary>
	/// Linearized node used for linear tree representation
	/// </summary>
	struct LinearNode
	{
		Aabb aabb;	 //bounds of node
		union
		{
			int dataIndex;			// used in leaf
			int rightChildOffset;   // used in internal - offset for right child, 
									// left child is on next index
		};
		uint16_t dataCount;		// 0 = internal node
		uint8_t axis;          // separation axis (0 - x, 1 - y, ...)
		uint8_t pad[1];        // ensure 32 byte total size
	};

	/// <summary>
	/// Wrapper around inserted data
	/// </summary>
	struct DataInfo
	{
		T val;
		Aabb aabb;
	};

public:
	enum class SplitMethod { Middle, EqualCounts };

	BVH(SplitMethod splitMethod = SplitMethod::Middle);
	~BVH();

	void Clear();

	template <typename Point>
	std::vector<T> GetPointInside(const Point & pt);

	void Add(const T & val, const Aabb & aabb);
	void Build();

protected:
	SplitMethod splitMethod;

	int nodesCount;

	std::vector<DataInfo> data;
	std::vector<DataInfo *> orderedData;

	LinearNode * nodes;

	Node * BuildRecursively(int start, int end);
	int SplitMiddle();
	void ClearTree(Node * root);

	int CreateLinearRepresentation(Node * node, int & nodeOffset);

};

/// <summary>
/// ctor
/// You can select splitting method for children bounding boxes
/// based on parent
/// todo: SAH (Surface Area Heuristic)
/// </summary>
/// <param name="splitMethod"></param>
template <typename T, typename Aabb>
BVH<T, Aabb>::BVH(SplitMethod splitMethod) :
	splitMethod(splitMethod),
	nodesCount(0),
	nodes(nullptr)
{
}


/// <summary>
/// dtor
/// </summary>
template <typename T, typename Aabb>
BVH<T, Aabb>::~BVH()
{
	this->data.clear();
	this->Clear();
}

/// <summary>
/// Clear tree and release linearized nodes
/// </summary>
template <typename T, typename Aabb>
void BVH<T, Aabb>::Clear()
{
	if (this->orderedData.size() == 0)
	{
		return;
	}

	this->orderedData.clear();

	delete[] nodes;
	nodes = nullptr;

	this->nodesCount = 0;
}

/// <summary>
/// When BVH is constructed, pointer-based tree is first created
/// -> release this tree starting from root
/// </summary>
/// <param name="root"></param>
template <typename T, typename Aabb>
void BVH<T, Aabb>::ClearTree(Node * root)
{
	std::stack<Node *> s;

	Node * tmp = root;

	do
	{
		// Move to leftmost node 
		while (tmp)
		{
			// Push root's right child and then root to stack. 
			if (tmp->right)
			{
				s.push(tmp->right);
			}
			s.push(tmp);

			// Set root as root's left child   
			tmp = tmp->left;
		}

		// Pop an item from stack and set it as root     
		tmp = s.top();
		s.pop();

		auto tmp2 = (s.empty()) ? nullptr : s.top();

		// If the popped item has a right child and the right child is not 
		// processed yet, then make sure right child is processed before root 
		if (tmp->right && tmp2 == tmp->right)
		{
			s.pop();  // remove right child from stack 
			s.push(tmp);  // push root back to stack 
			tmp = tmp->right; // change root so that the right  
								// child is processed next 
		}
		else
		{
			//process node
			delete tmp;
			tmp = nullptr;
		}
	} while (!s.empty());

	root = nullptr;
}

/// <summary>
/// Get all values which bounding box contain input point
/// </summary>
/// <param name="pt"></param>
/// <returns></returns>
template <typename T, typename Aabb>
template <typename Point>
std::vector<T> BVH<T, Aabb>::GetPointInside(const Point & pt)
{
	std::vector<T> res;

	if (this->nodes == nullptr)
	{
		return res;
	}

	int stackTop = 0;
	int currentNodeIndex = 0;
	int nodesToVisit[64];

	nodesToVisit[stackTop] = currentNodeIndex;
	stackTop++;

	while (stackTop != 0)
	{
		stackTop--;
		currentNodeIndex = nodesToVisit[stackTop];
		const LinearNode *node = &nodes[currentNodeIndex];


		if (node->aabb.IsInside(pt))
		{
			if (node->dataCount > 0)
			{
				//leaf
				for (int i = 0; i < node->dataCount; i++)
				{
					const DataInfo * di = orderedData[node->dataIndex + i];
					if (di->aabb.IsInside(pt))
					{
						res.push_back(di->val);
					}
				}
			}
			else
			{
				//internal node
				//we need to visit left and right child

				nodesToVisit[stackTop] = currentNodeIndex + 1;
				stackTop++;

				nodesToVisit[stackTop] = node->rightChildOffset;
				stackTop++;
			}
		}
	}

	return res;
}

/// <summary>
/// Add element to tree
/// Element is added, but tree is not build
/// If tree is already build, current tree is cleared
/// </summary>
/// <param name="val"></param>
/// <param name="aabb"></param>
template <typename T, typename Aabb>
void BVH<T, Aabb>::Add(const T & val, const Aabb & aabb)
{
	this->Clear();
	this->data.push_back({ val, aabb });
}

/// <summary>
/// Create tree from added elements
/// </summary>
template <typename T, typename Aabb>
void BVH<T, Aabb>::Build()
{
	//build classic binary tree
	Node * root = this->BuildRecursively(0, static_cast<int>(this->data.size()));

	if (root == nullptr)
	{
		return;
	}

	//create linear representation of tree
	//used DFS to fill array (pre-order)
	nodes = new LinearNode[this->nodesCount];
	int nodeOffset = 0;
	this->CreateLinearRepresentation(root, nodeOffset);

	//free classic binary tree from memory
	this->ClearTree(root);
}

/// <summary>
/// Build standard tree representation (with pointers)
/// - use recursive algorithm
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
/// <returns></returns>
template <typename T, typename Aabb>
typename BVH<T, Aabb>::Node * BVH<T, Aabb>::BuildRecursively(int start, int end)
{
	int itemsCount = end - start;

	if (itemsCount == 0)
	{
		//no item at all
		return nullptr;
	}

	this->nodesCount++;

	

	if (itemsCount == 1)
	{
		//only one item - insert it as leaf

		int firstDataIndex = static_cast<int>(orderedData.size());
		for (int i = start; i < end; i++)
		{
			orderedData.push_back(&data[i]);
		}

		return Node::CreateLeaf(firstDataIndex, 1, data[start].aabb);
	}


	//calculate bounding box of multiple objects

	auto cI = data[start].aabb.GetCentroid();
	Aabb bounds = data[start].aabb;
	Aabb centroidBounds(cI, cI);

	for (int i = start + 1; i < end; i++)
	{
		bounds = Aabb::CreateUnion(bounds, data[i].aabb);
		centroidBounds = Aabb::CreateUnion(centroidBounds, data[i].aabb.GetCentroid());
	}
	int dim = centroidBounds.GetMaxExtent();


	if (centroidBounds.max[dim] == centroidBounds.min[dim])
	{
		//all centroids are along the line or same
		int firstDataIndex = static_cast<int>(orderedData.size());
		for (int i = start; i < end; i++)
		{
			orderedData.push_back(&data[i]);
		}

		return Node::CreateLeaf(firstDataIndex, itemsCount, bounds);
	}
	else
	{
		int mid = (start + end) / 2;

		switch (splitMethod)
		{
		case SplitMethod::Middle:
		{
			//std::partition:
			//Rearranges the elements from the range[first, last), 
			//in such a way that all the elements for which 
			//pred returns true precede all those for which it returns false.
			//The iterator returned points to the first element of the second group.

			auto pmid = (centroidBounds.min[dim] + centroidBounds.max[dim]) / 2;
			DataInfo * midPtr = std::partition(&data[start], &data[end - 1] + 1,
				[dim, pmid](const DataInfo &pi) {
				return pi.aabb.GetCentroid()[dim] < pmid;
			});

			mid = static_cast<int>(midPtr - &data[0]);

			if (mid != start && mid != end)
			{
				break;
			}

		}
		case SplitMethod::EqualCounts:
		{
			mid = (start + end) / 2;
			std::nth_element(&data[start], &data[mid], &data[end - 1] + 1,
				[dim](const DataInfo &a, const DataInfo &b) {
				return a.aabb.GetCentroid()[dim] < b.aabb.GetCentroid()[dim];
			});
			break;
		}

		}


		Node * left = BuildRecursively(start, mid);
		Node * right = BuildRecursively(mid, end);

		return Node::CreateInner(dim, left, right);
	}
}

/// <summary>
/// Create linear representation of tree
/// Tree nodes are linearized to 1D array
/// in DFS order
/// Parents left child has index parent + 1
/// Right child index must be stored as offset
/// </summary>
/// <param name="node"></param>
/// <param name="nodeOffset"></param>
/// <returns></returns>
template <typename T, typename Aabb>
int BVH<T, Aabb>::CreateLinearRepresentation(Node * node, int & nodeOffset)
{
	LinearNode * n = &this->nodes[nodeOffset];
	n->aabb = node->aabb;

	int curOffset = nodeOffset;
	nodeOffset++;

	if (node->dataCount > 0)
	{
		//leaf
		n->dataIndex = node->dataIndex;
		n->dataCount = node->dataCount;
	}
	else
	{
		//internal node
		n->axis = n->axis;
		n->dataCount = 0;
		this->CreateLinearRepresentation(node->left, nodeOffset);
		n->rightChildOffset = this->CreateLinearRepresentation(node->right, nodeOffset);
	}

	return curOffset;
}



#endif
