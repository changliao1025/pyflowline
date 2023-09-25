#
#    tinyr - a 2D-RTree implementation in Cython
#    Copyright (C) 2011  Matthias Simon
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


'''
Fast 2D R-Tree implementation in cython.

The implementation is based on the original paper [Gut84]. 



[Gut84] R-trees: A Dynamic Index Structure for Spatial Searching, A. Guttman,
        Proc ACM SIGMOD Int, 47-57, 1984
        
[Bec90] The R*-tree: An efficient and robust access method for points and
        rectangles, N. Beckmann, H.-P. Kriegel, R. Schneider, and B. Seeger,
        ACM SIGMOD, 322-331, 1990
'''

# cython: profile=False

import itertools


from libc.float cimport DBL_MAX as Dfloat_max, DBL_MIN as Dfloat_min

ctypedef double Dfloat
ctypedef unsigned Duint


###############################################################################
#--- forward declarations
###############################################################################

cdef class _Record
cdef class _ChildRecord(_Record)
cdef class _InnerRecord(_Record)
cdef class _Node
cdef class RTree
cdef class _Iterator
cdef class _ChildRecordIter(_Iterator)
cdef class _KVIIter(_Iterator)
cdef class _KeysIter(_KVIIter)
cdef class _ValuesIter(_KVIIter)
cdef class _ItemsIter(_KVIIter)

class RTreeInvalid(Exception):
    pass

###############################################################################

cdef inline void init_Dfloat4(Dfloat *f):
    # array initialization - compiler cries if we don't do that
    f[0] = 0
    f[1] = 0
    f[2] = 0
    f[3] = 0

cdef inline void make_order(Dfloat *coords):
    cdef Dfloat switch
    if coords[0] > coords[2]:
        switch = coords[0]
        coords[0] = coords[2]
        coords[2] = switch
    if coords[1] > coords[3]:
        switch = coords[1]
        coords[1] = coords[3]
        coords[3] = switch

cdef void tuple_to_array_interleaved(object t, Dfloat *coords):
    for i in range(4):
        coords[i] = <Dfloat>t[i]
    make_order(coords)

cdef void tuple_to_array_normal(object t, Dfloat *coords):
    coords[0] = <Dfloat>t[0]
    coords[2] = <Dfloat>t[1]
    coords[1] = <Dfloat>t[2]
    coords[3] = <Dfloat>t[3]
    make_order(coords)


cdef tuple array_to_tuple_interleaved(Dfloat *coords):
    return tuple([ coords[i] for i in range(4)])

cdef tuple array_to_tuple_normal(Dfloat *coords):
    return tuple([ coords[0], coords[2], coords[1], coords[3] ])

cdef inline _InnerRecord parent_record(_Node node):
    cdef:
        _InnerRecord e
        
    for e in node.parent.records:
        if e.child == node:
            return e

cdef inline void common_boundaries(list records, Dfloat *target):
    cdef:
        Dfloat *coords
        _Record r
        
    target[0] = target[1] = Dfloat_max
    target[2] = target[3] = -Dfloat_min
    for r in records:
        coords = r.coords
        target[0] = min(target[0], coords[0])
        target[1] = min(target[1], coords[1])
        target[2] = max(target[2], coords[2])
        target[3] = max(target[3], coords[3])

cdef inline Dfloat area(Dfloat *coords):
    return (coords[2] - coords[0]) * (coords[3] - coords[1])

cdef class _Record:
    cdef:
        Dfloat coords[4]

    cdef inline bint overlaps(self, Dfloat *rect):
        return self.coords[0] < rect[2] and self.coords[2] > rect[0] and self.coords[1] < rect[3] and self.coords[3] > rect[1]
    
    cdef inline void copy_coords_to(self, Dfloat *coords):
        for i in range(4):
            coords[i] = self.coords[i]
    
    cdef _Record copy(self, _Node newparent):
        raise NotImplemented

cdef class _ChildRecord(_Record):
    cdef:
        object identifier
    
    def __init__(self, object identifier):
        self.identifier = identifier
    
    cdef _Record copy(self, _Node newparent):
        cdef:
            _ChildRecord ret
        
        ret = _ChildRecord(self.identifier)
        self.copy_coords_to(ret.coords)
        return ret
    
cdef class _InnerRecord(_Record):
    cdef:
        _Node child
    
    def __cinit__(self, _Node child):
        self.child = child

    cdef _Record copy(self, _Node newparent):
        cdef:
            _InnerRecord ret
        
        ret = _InnerRecord(self.child.copy(newparent))
        self.copy_coords_to(ret.coords)
        return ret
        
cdef class _Node(object):
    cdef:
        _Node parent
        list records
        bint is_leaf
    
    def __cinit__(self, _Node parent):
        self.parent = parent
        self.records = list()
    
    cdef _Node copy(self, _Node parent):
        cdef:
            _Node ret
        ret = _Node(parent)
        ret.records = [ (<_Record>rec).copy(ret) for rec in self.records ]
        ret.is_leaf = self.is_leaf
        return ret
        
    cdef inline void common_boundaries(self, Dfloat *target):
        common_boundaries(self.records, target)

    cdef inline _Record choose_subtree_least_enlargement(self, _Record ir):
        cdef:
            Dfloat least_enlargement, enlagrgement, current_area, target_area
            Dfloat combined_rectangle[4]
            _Record current_record, target
            list records
        
        init_Dfloat4(combined_rectangle) 
        
        target_area = least_enlargement = Dfloat_max
        records = [ir, None]
        
        
        for current_record in self.records:
            current_area = area(current_record.coords)
            
            records[1] = current_record
            
            common_boundaries(records, combined_rectangle)
            enlagrgement = area(combined_rectangle) - current_area
            
            if enlagrgement < least_enlargement:
                target = current_record
                target_area = area(current_record.coords)
                least_enlargement = enlagrgement
            elif enlagrgement == least_enlargement and current_area < target_area:
                target = current_record
                least_enlargement = target_area = current_area
            
        return target
    
    cdef void find_overlapping_leafs_recursive(self, list result, Dfloat *coords):
        cdef:
            _InnerRecord record
    
        if self.is_leaf:
            result.append(self)
        else:
            for record in self.records:
                if record.overlaps(coords):
                    record.child.find_overlapping_leafs_recursive(result, coords)
                    
    cdef void addChild(self, _Node node):
        cdef:
            _InnerRecord ir

        node.parent = self
        ir = _InnerRecord(node)
        node.common_boundaries(ir.coords)
        self.records.append(ir)

cdef list find_overlapping_leafs_recursive(RTree rtree, Dfloat *coords):
    cdef:
        _Record rec
        list result
    
    # don't know the whole surrounding rectangle of root, ask records for overlapping
    for rec in rtree.root.records:
        if rec.overlaps(coords):
            break
    else:
        return []
    
    result = list()
    rtree.root.find_overlapping_leafs_recursive(result, coords)
    return result

# performs worse than the recursive one - so we skip it
#cdef list find_overlapping_leafs_linear(RTree rtree, Dfloat *coords):
#    cdef:
#        _Node node
#        _InnerRecord ir
#        list nodes_to_check, next_level_nodes
#        unsigned level
#
#    nodes_to_check = [rtree.root]
#    next_level_nodes = []
#    level = 0
#    
#    while level != rtree._leaf_level:
#        while nodes_to_check:
#            node = <_Node>nodes_to_check.pop()
#            for ir in node.records:
#                if (<_InnerRecord>ir).overlaps(coords):
#                    next_level_nodes.append((<_InnerRecord>ir).child)
#        
#        nodes_to_check = next_level_nodes
#        next_level_nodes = []
#        level += 1
#    
#    return nodes_to_check

cdef int PickSeedsQuadratic(_Node node, list remaining, _Node newnode) except -1:
    cdef:
        Dfloat d, d_max, E1_E2[4]
        Dfloat J, a_E1, a_E2
        _Record E1_ret, E2_ret
        list combi
    
    init_Dfloat4(E1_E2)
    
    combi = [None, None]
    d_max = -Dfloat_max
    
    for E1, E2 in itertools.combinations(remaining, 2):
        
        combi[0] = E1
        combi[1] = E2

        common_boundaries(combi, E1_E2)
        J = area(E1_E2)
        a_E1 = area((<_Record>E1).coords) 
        a_E2 = area((<_Record>E2).coords)
        d = J - a_E1 - a_E2
        
        if d > d_max:
            E1_ret = E1
            E2_ret = E2

    remaining.remove(E1_ret)
    remaining.remove(E2_ret)
    
    node.records = [E1_ret]
    newnode.records = [E2_ret]
    
    return 0

#cdef int PickSeedsLinear(_Node node, list remaining, _Node newnode) except -1:
#    # specification intended to have such a method, but
#    # too complicated, since spec lets an unclear situation left:
#    # the two seed nodes may be the same. the spec doen't say anything about how
#    # to handle this situation... my guess 'take the next best one' would be too
#    # rude.
#    raise NotImplemented

cdef class RTree(object):
    '''RTree - structure for 2D-indexed items.
    
    @param interleaved: It True, all methods that accept coordinates expect it to be in format (x1, y1, x2, y2). 
                        If False, coordinates are supposed to be (x1, x2, y1, y2)
                        Default: True
    @param max_cap: Maximum of capacity of an index node. Default: 5
    @param min_cap: Minimum of entries an index node must hold. Default: 2
    '''
    cdef:
        _Node root
        Duint _leaf_level
        Duint _node_cnt
        int M_cap, m_cap
        
        list (*find_overlapping_leafs)(RTree rtree, Dfloat *coords)
        int (*PickSeeds)(_Node node, list remaining, _Node newnode) except -1
        void (*tuple_to_array)(object t, Dfloat *coords)
        tuple (*array_to_tuple)(Dfloat *coords)
    
    property min_cap:
        def __get__(self):
            return self.m_cap
        
    property max_cap:
        def __get__(self):
            return self.M_cap
    
    def __cinit__(self):
        self.root = _Node(None)
        self.root.is_leaf = True
        self._leaf_level = 0
        self._node_cnt = 0
    
    def __init__(self, interleaved=True, max_cap=5, min_cap=2):
        if min_cap > max_cap/2.:
            raise ValueError('for min/max capacities, it\'s required that min_cap <= max_cap/2.0')
        
        if interleaved:
            self.tuple_to_array = tuple_to_array_interleaved
            self.array_to_tuple = array_to_tuple_interleaved
        else:
            self.tuple_to_array = tuple_to_array_normal
            self.array_to_tuple = array_to_tuple_normal
        
        # the iterative method performed worse
        self.find_overlapping_leafs = find_overlapping_leafs_recursive
        
        self.PickSeeds = PickSeedsQuadratic

        self.M_cap = max_cap
        self.m_cap = min_cap
    
    def __len__(self):
        return self._node_cnt
    
    def __getitem__(self, item):
        if isinstance(item, (list, tuple)):
           if len(item) == 2:
                return self.search_surrounding(item)
           elif len(item) == 4:
                return self.search(item)
        raise AttributeError('indexing RTree only possible with rectangle or point (x, y)')
    
    def iterkeys(self):
        return _KeysIter(self)
    
    def itervalues(self):
        return _ValuesIter(self)
    
    def iteritems(self):
        return _ItemsIter(self)
    
    def get_info(self):
        return RTreeInfo(self)

    def insert(self, object key, object coords):
        '''Insert an item.
        
        @param key: an object to insert, where all keys should be unique (regarding !=)
        @param coords: 2D coordinates
        '''
        
        cdef:
            _ChildRecord cr
            _Node L
        
        if not isinstance(coords, (list, tuple)):
            raise TypeError('coordinates as list or tuple expected, got %s' % coords.__class__)
        if not len(coords) == 4:
            raise TypeError('len(coords) must be 4, length is %d' % len(coords))
        
        cr = _ChildRecord(key)
        self.tuple_to_array(coords, cr.coords)
        
        L = self._ChooseLeaf(cr)
        
        self.insert_at_node(cr, L)
        
        self._node_cnt += 1


    def search(self, object coords):
        '''Search overlapping items.
        
        @param coords: list or tuple of four values that make a rectangle
        @return: a list of identifiers whose coordinates overlap with coords
        '''
        cdef:
            _Node node
            _ChildRecord cr
            Dfloat coords_a[4]
            list leafnodes, result
            
        if not isinstance(coords, (list, tuple)):
            raise TypeError('coordinates as list or tuple expected, got %s' % coords.__class__)
        
        if len(coords) != 4:
            raise TypeError('len(coords) must be 4, len is %d' % len(coords))
        
        self.tuple_to_array(coords, coords_a)
        
        leafnodes = self.find_overlapping_leafs(self, coords_a)

        result = []
        for node in leafnodes:
            assert node.is_leaf
            for cr in node.records:
                if cr.overlaps(coords_a):
                    result.append(cr.identifier)
        
        return result
    
    def copy(self):
        cdef:
            RTree ret

        ret = RTree(max_cap=self.M_cap, min_cap=self.m_cap)
        ret._leaf_level = self._leaf_level
        ret._node_cnt = self._node_cnt
        ret.root = self.root.copy(None)
        return ret
        
    def search_surrounding(self, point):
        '''Search items that surround a point.
        
        @param item: a point in form (x, y)
        @return: a list of identifiers whose coordinates surround with point
        '''
        cdef:
            _Node node
            _ChildRecord cr
            Dfloat coords_a[4]
            list leafnodes, result
        
        if not isinstance(point, (list, tuple)):
            raise TypeError('coordinates as list or tuple expected, got %s' % point.__class__)
        
        if len(point) == 2:
            point = 2 * point    
        else:
            raise TypeError('len(point) must be 2, len is %d' % len(point))
        
        self.tuple_to_array(point, coords_a)
        
        leafnodes = self.find_overlapping_leafs(self, coords_a)

        result = []
        for node in leafnodes:
            assert node.is_leaf
            for cr in node.records:
                if cr.overlaps(coords_a):
                    result.append(cr.identifier)
        
        return result

    cdef inline int insert_at_node(self, _Record r, _Node L) except -1:
        cdef:
            _Node LL = None
    
        if self.is_full(L):
            LL = self.SplitNode(L, r)
        else:
            L.records.append(r)
            
        L, LL = self._AdjustTree(L, LL)

        if LL != None:
            self.root = _Node(None)
            self.root.addChild(L)
            self.root.addChild(LL)
            self._leaf_level += 1
        
        return 0
    
    def remove(self, identifier, coords):
        '''Delete an item
        
        @param identifier: the item that should be removed
        @param coords: the coordinates of item
        '''
        cdef:
            _Node leaf
            Dfloat coords_a[4]
        
        if not isinstance(coords, (list, tuple)):
            raise TypeError('coordinates as list or tuple expected, got %s' % coords.__class__)
        if not len(coords) == 4:
            raise TypeError('len(coords) must be 4, len is %d' % len(coords))
        
        self.tuple_to_array(coords, coords_a)

        leaf = self.FindLeaf_and_remove(identifier, coords_a)
        
        if leaf is None:
            raise KeyError()
        
        self.after_remove(leaf)
        
        self._node_cnt -= 1
    
    cdef inline void after_remove(self, _Node node):
        self.CondenseTree(node)
        if not self.root.is_leaf and len(self.root.records) == 1:
            self.root = (<_InnerRecord>self.root.records.pop()).child
            self._leaf_level -= 1
    
    cdef inline _Node FindLeaf_and_remove(self, object identifier, Dfloat *coords):
        # A bit different that in spec: we find overlapping leaf nodes with the 
        # method the user wants to find it (recursive or iterative). Spec only
        # defines a recursive way.
        # It also removes the item that should be removed.
        cdef:
            list leafs
            _Node node
            _ChildRecord rec
        
        leafs = self.find_overlapping_leafs(self, coords)
        for node in leafs:
            assert node.is_leaf
            for rec in node.records:
                if rec.identifier == identifier and rec.overlaps(coords):
                    node.records.remove(rec) # D2
                    return node 
        
        return None

    cdef inline int CondenseTree(self, _Node node) except -1:
        cdef:
            _Node P, N
            _InnerRecord EN
            _Record rec
            _InnerRecord ir
            _ChildRecord crec
            list removed_childrecords, removed_innerrecords, reclist
        
        removed_innerrecords = list()
        
        N = node
        while not N is self.root:
            
            P = N.parent
            EN = parent_record(N)
            
            if len(N.records) < self.m_cap:
                assert EN is not None, str(N.is_leaf)
                
                P.records.remove(EN)
                if N.is_leaf:
                    assert removed_childrecords is None
                    removed_childrecords = N.records
                else:
                    removed_innerrecords.extend(N.records)
            
            else:
                common_boundaries(N.records, EN.coords)
            
            N = P
        
        
        for ir in removed_innerrecords:
            node = self.find_inner_node(ir)
            ir.child.parent = node
            self.insert_at_node(ir, node)
        
        if removed_childrecords is not None:
            for crec in removed_childrecords:
                node = self._ChooseLeaf(crec)
                self.insert_at_node(crec, node)
        
        return 0
    
    cdef _Node find_inner_node(self, _InnerRecord record):
        cdef:
            _Node node
            _InnerRecord rec
            Duint subtree_len, level_reverse
        
        subtree_len = self.subtree_len(record)
        level_reverse = self._leaf_level
        
        node = self.root
        
        while level_reverse != subtree_len:
            rec = node.choose_subtree_least_enlargement(record)
            node = rec.child
            level_reverse -= 1
            
        return node
        
     
    cdef Duint subtree_len(self, _InnerRecord rec):
        cdef:
            _Node node
            Duint slen
            
        slen = 1
        
        node = rec.child
        while not node.is_leaf:
            assert node.records
            node = (<_InnerRecord>node.records[0]).child
            slen += 1
        
        return slen
     

    cdef inline bint is_full(self, _Node node):
        return len(node.records) == self.M_cap

    cdef inline _Node _ChooseLeaf(self, _Record ir):
        cdef:
            _Node N
            _InnerRecord F
        
        N = self.root
        while True:
            if N.is_leaf:
                break
            r = N.choose_subtree_least_enlargement(ir)
            N = (<_InnerRecord>r).child
        return N

    cdef inline tuple _AdjustTree(self, _Node L, _Node LL):
        cdef:
            _Node N, NN, P
            _Record EN, e
        
        N = L
        NN = LL

        while N != self.root:
            P = N.parent

            # search for the entry in parent of N that holds N
            EN = parent_record(N)
            
            assert EN != None, 'no parent entry holds the child'

            N.common_boundaries(EN.coords)

            if NN != None:
                ENN = _InnerRecord(NN)
                NN.common_boundaries(ENN.coords)
                if not self.is_full(P):
                    P.records.append(ENN)
                    NN = None
                else:
                    NN = self.SplitNode(P, ENN)

            N = P
        
        return (N, NN)

    cdef _Node SplitNode(self, _Node node, _Record ir):
        cdef:
            list remaining
            _Node newnode
            unsigned remaining_count
            _InnerRecord r
        
        remaining = node.records
        remaining.append(ir)
        node.records = list()
        
        newnode = _Node(node.parent)

        self.PickSeeds(node, remaining, newnode) 

        remaining_count = len(remaining)
        while remaining_count > 0:
            if len(node.records) + remaining_count - 1 < self.m_cap:
                node.records.extend(remaining)
                break
            if len(newnode.records) + remaining_count - 1 < self.m_cap:
                newnode.records.extend(remaining)
                break
            
            self.PickNext(node, remaining, newnode) 
            remaining_count -= 1

        if node.is_leaf:
            newnode.is_leaf = True
        else:
            # child records are inner records - re-set parent of them to newnode
            for r in newnode.records:
                r.child.parent = newnode
        
        return newnode

    cdef inline int PickNext(self, _Node node, list remaining, _Node newnode) except -1:
        cdef:
            Dfloat area_group_L, area_group_LL
            Dfloat d1_next, d2_next # for QS3 decision
            Dfloat d1, d2, d_diff, d_diff_max
            Dfloat coords[4]
            _Record E_next, E, r
        
        area_group_L = 0
        area_group_LL = 0
        d1_next = 0
        d2_next = 0
        init_Dfloat4(coords)
        
        if len(remaining) == 1:
            E_next = remaining.pop()
        else:
            for r in node.records:
                area_group_L += area(r.coords)
    
            for r in newnode.records:
                area_group_LL += area(r.coords)
            
            d_diff_max = -Dfloat_min
    
            node.records.append(None)
            newnode.records.append(None)
            
            for E in remaining:
                
                # temporary add E to self.records / newnode.records as parameter for common_boundires()
                # -> don't init a new list of records
                node.records[-1] = E
                newnode.records[-1] = E
                
                common_boundaries(node.records, coords)
                d1 = area(coords) - area_group_L
                
                common_boundaries(newnode.records, coords)
                d2 = area(coords) - area_group_LL
                
                d_diff = d1 - d2
                d_diff = d_diff if d_diff >= 0 else -d_diff # abs diff
                
                if d_diff > d_diff_max:
                    d_diff_max = d_diff
                    E_next = E
    
                    d1_next = d1
                    d2_next = d2
            
            node.records.pop()
            newnode.records.pop()
            
            remaining.remove(E_next)
            
        # QS3
        if d1_next < d2_next:
            node.records.append(E_next)
        elif d1_next > d2_next:
            newnode.records.append(E_next)
        elif len(node.records) > len(newnode.records):
            newnode.records.append(E_next)
        else:
            node.records.append(E_next)
            
    
    def valid(self):
        # numbers in comments refer to invariant properties of an RTree - taken from [Gut84]
        cdef:
            Dfloat tmp_coords[4]
            _Node node, n
            Duint depth, item_count
            
        nodedepht = None
        item_count = 0
        leafnodes_and_depth = []
        
        for depth, node in _NodeIter(self):
            if node.is_leaf:
                for rec in node.records:
                    assert isinstance(rec, _ChildRecord), depth
                
                # (1)
                if not (node is self.root or len(node.records) <= self.M_cap and len(node.records) >= self.m_cap):
                    raise RTreeInvalid('non-leaf node has invalid number of entries: %d <= len(%s) <= %d == False' % (self.m_cap, str(node), self.m_cap))
                        
                # (2)
                # Always True since I _IS_ the N-dimensional representation of the data object.
                
                # prepare (6)
                if nodedepht is None:
                    nodedepht = depth
                leafnodes_and_depth.append((depth, node))
                
                item_count += len(node.records)
                
                # height is right
                if self._leaf_level != depth:
                    raise RTreeInvalid('inconsistent node depth: %d != %d' % (self._leaf_level, depth))
                
            else:
                # (3)
                if not node is self.root:
                    if not (len(node.records) <= self.M_cap and len(node.records) >= self.m_cap):
                        raise RTreeInvalid('non-leaf node has invalid number of entries: %d <= len(%s) <= %d == False' % (self.m_cap, str(node), self.m_cap))
                        
                # (4)
                for rec in node.records:
                    if not isinstance(rec, _InnerRecord) or \
                        not isinstance((<_InnerRecord>rec).child, _Node):
                            raise RTreeInvalid('misspositioned inner node')
                        
                    common_boundaries((<_InnerRecord>rec).child.records, tmp_coords)
                    
                    if area((<_Record>rec).coords) != area(tmp_coords):
                        raise RTreeInvalid('wrong combined area of records of a node')

                # double-links cover each other
                for ir in node.records:
                    if not (<_InnerRecord>ir).child.parent is node:
                        raise RTreeInvalid('double-link doesn\'t cover itself')
                
        # (5)
        if not self.root.is_leaf:
            
            if len(self.root.records) < 2:
                raise RTreeInvalid('root must have at least 2 records, not %d.' % len(self.root.records))
            
        # (6)
        for d, leaf in leafnodes_and_depth:
            if d != nodedepht:
                raise RTreeInvalid('leaf level is not equal for all leafs')
        
        # count of items is right
        if item_count != self._node_cnt:
            raise RTreeInvalid('inconsistent item count %d != %d' % (item_count, self._node_cnt))
        
        return True

###############################################################################
#--- Iterators and RTree Info
###############################################################################

cdef class _Iterator:
    cdef RTree rtree

    def __cinit__(self, RTree rtree):
        self.rtree = rtree

    def __iter__(self):
        return self

cdef class _NodeIter(_Iterator):
    cdef:
        list nodeset
    
    def __cinit__(self, RTree rtree):
        if rtree.root:
            self.nodeset = [(0, rtree.root)]
        else:
            self.nodeset = []
    
    def __next__(self):
        cdef:
            _Node next
            Duint level
            
        try:
            level, next = <_Node>self.nodeset.pop(0)
        except IndexError:
            raise StopIteration
        
        if not next.is_leaf:
            self.nodeset = [ (level+1, (<_InnerRecord>r).child) for r in next.records ] + self.nodeset
        
        return level, next

cdef  class _RectangleIter(_Iterator):
    cdef:
        _NodeIter ni
        list leafrecords
        int l_level
    
    def __cinit__(self, RTree rtree):
        self.ni = _NodeIter(rtree)

    def __next__(self):
        cdef:
            _Node nextnode
            Dfloat coords[4]
            int level
            
        if self.leafrecords:
            rec = self.leafrecords.pop()
            return (self.l_level, self.rtree.array_to_tuple((<_ChildRecord>rec).coords) )
        
        level, nextnode = self.ni.next()
        
        if nextnode.is_leaf:
           self.leafrecords = list(nextnode.records)
           self.l_level = level+1
        
        common_boundaries(nextnode.records, coords)
        return (level, self.rtree.array_to_tuple(coords))
    

cdef class _ChildRecordIter(_Iterator):
    cdef:
        list nodeset, current_crs
        object (*returner)(_ChildRecordIter this, _ChildRecord cr)
    
    def __cinit__(self, RTree rtree):
        if rtree.root:
            self.nodeset = [rtree.root]
        else:
            self.nodeset = []
        self.current_crs = []
    
    def __next__(self):
        cdef:
            _Node node
            _ChildRecord cr
        
        if not self.current_crs:
            try:
                node = self.nodeset.pop(0)
            except IndexError:
                raise StopIteration
            
            while not node.is_leaf:
                self.nodeset = [ (<_InnerRecord>r).child for r in node.records ] + self.nodeset
                try:
                    node = self.nodeset.pop(0)
                except IndexError:
                    raise StopIteration
            
            self.current_crs.extend(node.records)

        try:
            return self.current_crs.pop(0)
        except IndexError:
            # should not happen 
            raise StopIteration
        

cdef class _KVIIter(_Iterator):
    cdef:
        _ChildRecordIter cri
    
    def __cinit__(self, RTree rtree):
        self.cri = _ChildRecordIter(rtree)

cdef class _KeysIter(_KVIIter):
    def __next__(self):
        return self.rtree.array_to_tuple((<_ChildRecord>self.cri.next()).coords)

cdef class _ValuesIter(_KVIIter):
    def __next__(self):
        return (<_ChildRecord>self.cri.next()).identifier

cdef class _ItemsIter(_KVIIter):
    def __next__(self):
        cdef _ChildRecord cr
        cr = self.cri.next()
        return ( self.rtree.array_to_tuple(cr.coords), cr.identifier )


cdef class RTreeInfo:
    '''Provides some meta-info about an existing RTree
    
    @param rtree: an RTree object
    '''
    cdef:
        RTree rtree

    property levels:
        def __get__(self):
            return self.rtree._leaf_level + 1 # count leaf rectangles, too 
    
    def __init__(self, RTree rtree):
        self.rtree = rtree
    
    def iter_rectangles(self):
        '''Iterate over key coordinates and combined frames.
        
        Combined frames are 2D-indexes under which rectangle entries are stored.
        The iterator returns a tuple (level, coordinates) where level is the
        level of the index node or entry in the tree structure.
        '''
        return _RectangleIter(self.rtree)

    property common_boundary:
        def __get__(self):
            cdef:
                Dfloat coords[4]
            if not self.rtree.root:
                return 0
            common_boundaries(self.rtree.root.records, coords)
            return self.rtree.array_to_tuple(coords)

    property width:
        def __get__(self):
            cdef:
                Dfloat coords[4]
            if not self.rtree.root:
                return 0
            common_boundaries(self.rtree.root.records, coords)
            return coords[2]-coords[0]
    
    property height:
        def __get__(self):
            cdef:
                Dfloat coords[4]
            if not self.rtree.root:
                return 0
            common_boundaries(self.rtree.root.records, coords)
            return coords[3]-coords[1]
    
    def to_dot(self, filelike):
        '''Writes simple dot-code into filelike that represents the tree structure.'''
        cdef:
            _Node node
            _InnerRecord ir
            _ChildRecord cr
            int level
            unsigned nid, rid, chid, lid
            list node_record_leaf, node_record_child
        
        filelike.write('digraph RTree {\n    node [shape=\"none\"]; \n')
        
        node_record_leaf = list()
        node_record_child = list()
        
        for level, node in _NodeIter(self.rtree):
            if node.is_leaf:
                for cr in node.records:
                    node_record_leaf.append((id(node), id(cr), id(cr.identifier)))
                    
            else:
                for ir in node.records:
                    node_record_child.append((id(node), id(ir), id(ir.child)))
                
            filelike.write('    %d [label=<<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\"><TR>' % (id(node)))
            for record in node.records:
                filelike.write('<TD PORT=\"%d\" CELLSPACING=\"0\">%d</TD>' % (id(record), id(record)))
            filelike.write('</TR></TABLE>>];\n')
        
        filelike.write('\n    node [shape=\"ellipse\"]; \n')
        for nid, rid, lid in node_record_leaf:
            filelike.write('    %d;\n' % lid)
        
        filelike.write('\n')
        
        for nid, rid, chid in node_record_child:
            filelike.write('    %d->%d [tailport=\"%d:s\"]\n' % ( nid, chid, rid ))
        
        filelike.write('\n')
        
        for nid, rid, lid in node_record_leaf:
            filelike.write('    %d->%d [tailport=\"%d:s\"]\n' % ( nid, lid, rid ))
        
        filelike.write('}')
        




