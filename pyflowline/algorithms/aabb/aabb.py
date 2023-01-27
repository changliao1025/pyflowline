class AABBNode:
    def __init__(self, min_x, min_y, max_x, max_y, left=None, right=None):
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y
        self.left = left
        self.right = right


class AABBTree:
    def __init__(self):
        self.root = None

    def insert(self, min_x, min_y, max_x, max_y):
        if self.root is None:
            self.root = AABBNode(min_x, min_y, max_x, max_y)
            return

        current = self.root
        while current is not None:
            if min_x < current.min_x:
                if current.left is None:
                    current.left = AABBNode(min_x, min_y, max_x, max_y)
                    return
                current = current.left
            else:
                if current.right is None:
                    current.right = AABBNode(min_x, min_y, max_x, max_y)
                    return
                current = current.right

    def query(self, min_x, min_y, max_x, max_y):
        def intersect(a, b):
            return not (a.min_x > b.max_x or a.min_y > b.max_y or a.max_x < b.min_x or a.max_y < b.min_y)

        def traverse(node):
            if node is None:
                return []

            if intersect(node, (min_x, min_y, max_x, max_y)):
                return [node] + traverse(node.left) + traverse(node.right)
            else:
                return traverse(node.left) + traverse(node.right)

        return traverse(self.root)
