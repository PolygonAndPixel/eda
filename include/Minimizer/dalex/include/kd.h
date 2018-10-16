/*
Based on description in

Bentley, J. L.
Communications of the Association for Computing Machinery
volume 18
p 509
(1975)
*/

#ifndef KD_H
#define KD_H

#include "time.h"
#include "containers.h"
#include "goto_tools.h"

class kd_tree{
    /*
    This class stores points in an n-dimensional parameter space in a
    KD-tree for easy nearest neighbor searching.
    */
 
    /*this friend declaration makes it possible for the Gaussian Process class 
    to access the data stored in the tree*/
    friend class gp;
    friend class box;
    
    public:
        index_t ktests;
        
        /*Build a tree out of the list of points provided (each row is a 
        different point in parameter space; therefore, the number of columns
        is the dimensionality of the parameter space)*/
        kd_tree(array_2d<value_t>&);
        
        /*Build a tree as above.  The array_1d<value_t> arguments are minimum
        and maximum values of each parameter in parameter space.  These
        are not bounds.  max-min is used to normalize distances in parameter
        space when searching for nearest neighbors.*/
        kd_tree(array_2d<value_t>&,array_1d<value_t>&,array_1d<value_t>&);
        
        kd_tree(const kd_tree&);
        
        kd_tree();
        
        ~kd_tree();
        
        void copy(const kd_tree&);
        
        /*
        These routines provide the back end for building the KD tree
        */
        void build_tree(array_2d<value_t>&);
        void build_tree(array_2d<value_t>&,array_1d<value_t>&,array_1d<value_t>&);
 
        /*
        These routines will set the maximum and minimum values (used for
        normalizing parameter distances; see above)
        */
        void set_max(index_t,value_t);
        void set_min(index_t,value_t);
        
        /*check to make sure the tree is properly constructed;
        if it is, set the global variable diagnostic=1;
        if not, set diagnostic=0*/
        void check_tree();
        
        /*check to make sure the part of the tree that descendes from
        the node specified by the index_t is properly constructed*/
        void check_tree(index_t);
        
        /*
        The distance routines all return parameter space distances between
        points.  These points do not have to be on the tree, but they can be.
        
        The distances are Euclidean, except that they are normalized by the
        values stored in the private global variables maxs and mins, i.e.
        
        distance = sqrt{\sum_i [(p1_i - p2_i)/(maxs_i - mins_i)]^2}
        */
        
        /*the parameter space distance between arbitrary points*/
        value_t distance(const array_1d<value_t>&, const array_1d<value_t>&);
 
        /*the parameter space distance between and arbitrary point and a node
        on the tree*/
        value_t distance(const array_1d<value_t>&,index_t);
        value_t distance(index_t,const array_1d<value_t>&);
        
        /*the parameter space distance between two nodes on the tree*/
        value_t distance(index_t,index_t);

        /*fill the array-1d<value_t> with the node specified by the index_t*/
        void get_pt(index_t,array_1d<value_t>&);
        
        /*return a node on the tree one dimension at a time; the first
        index_t specifies the node; the second index_t specifies the dimension*/
        value_t get_pt(index_t,index_t);
        
        /*return a pointer to a node on the tree*/
        array_1d<value_t> get_pt(index_t);
        
        /*return the number of points stored in the tree*/
        index_t get_pts();

        void write_tree(char*);
        
        /*add a point to the tree*/
        void add(const array_1d<value_t>&);
        
        /*removes a node from the tree and then rebuilds the tree*/
        void remove(index_t);
        
        /*counts the number of nodes descended from a given parent*/
        void count(index_t,index_t*);
        
        /*
        The nn_srch routins perform the nearest neighbor searches.
        The first argument specifies the point whose nearest neighbors
        you want to find (if it is just an index_t, you are looking for the
        nearest neighbors of a node onthe tree)
        
        The second argument specifies the number of nearest neighbors to find.
        
        The third argument will store the indices of the nearest neighbors.
        
        The fourth argument will store the (normalized) parameter space
        distances to those nearest neighbors.
        
        Remember: all distances are normalized by maxs-mins (see documentation
        of the distance() routines)
        */
        void nn_srch(const array_1d<value_t>&,index_t,array_1d<index_t>&,array_1d<value_t>&);
        void nn_srch(index_t,index_t,array_1d<index_t>&,array_1d<value_t>&);
        
        /*
        kernel_srch will do a neighbor search centered on a point (the first
        argument).  It will return all of the points located within a hyperbox 
        whose half-width is specified by the second argument.  The third
        argument stores the indices of the found neighbors.  The routine returns
        the number of found neighbors.
        
        NOTE: THIS IS NOT WELL-TESTED
        */
        index_t kernel_srch(array_1d<value_t>&,array_1d<value_t>&,array_1d<index_t>&);
        
        /*
        radial_srch performs a neighbor search centered on a point (the first
        argument).  It returns all points within a given (normalized) radius
        (the second argument).  The third argument stores the indices of the
        found neighbors.  The routine returns the number of found neighbors.
        
        NOTE: THIS IS NOT WELL-TESTED
        */
        index_t radial_srch(array_1d<value_t>&,value_t,array_1d<index_t>&);
        
        /*return the dimensionality of the parameter space*/
        index_t get_dim();
        
        /*return diagnostic, the global variable that logs whether or not
        the tree was properly constructed*/
        index_t get_diagnostic();
        
        /*return the maximum and minimum values in a dimension of parameter
        space*/
        value_t get_max(index_t);
        value_t get_min(index_t);
        
        value_t get_search_time();
        index_t get_search_ct();
        
        index_t get_search_ct_solo();
        value_t get_search_time_solo();
        
        void set_search_ct(index_t);
        void set_search_time(value_t);
        
        void set_search_ct_solo(index_t);
        void set_search_time_solo(value_t);
        
    private:
        
        value_t search_time,search_time_solo;
        index_t search_ct,search_ct_solo;
        
        /*a global variable logging whether or not the tree was properly
        constructed; diagnostic=1 if the tree is correct; diagnostic=0 if not*/
        index_t diagnostic;
        
        /*
        The array_2d<index_t> tree stores the structure of the KD tree.
        Each row corresponds to a data point stored in the tree.
        There are four columns.
        
        The 0th column stores the dimension that the tree is branching on
        at that point.
        
        The 1st column is the left-hand daughter (points whose value in the
        dimension specified in the 0th column are less than the current
        point).
        
        The 2nd column is the right-hand daughter (points whose value in
        the dimension specified in the 0th column are greater than or
        equal to the current point).
        
        The 3rd column is the parent of the current point.
        
        Columns are set to -1 if they have no meaningful answer (i.e. the parent
        of the original point or the daughter of a terminal node).
        */
        array_2d<index_t> tree;
        
        /*
        The array_2d<value_t> data contains the points that are stored in this
        tree.
        */
        array_2d<value_t> data;
        
        /*
        maxs and mins are maximum and minimum values in each dimension of
        parameter space.  Note: these are not bounds.  The difference max-min is
        used to normalize parameter space distances when trying to find nearest
        neighbors.  Setting all maxs=1 and all mins=0 will result in nearest
        neighbors reckoned by unnormalized parameter space distances.
        */
        array_1d<value_t> maxs,mins;
        
        /*masterparent is the point that is the first node of the tree*/
        index_t masterparent;
        
        /*a global variable to define the tolerance with which dimensions
        are sorted on the tree*/
        value_t tol;
        
        /*a global variable used by kernel_srch to keep track of how many
        points were found within the kernel
        
        also used by radial_srch
        */
        index_t nkernel;
        
        /*this provides the backend for check_tree;
        see source code in kd.cpp for more details*/
        void confirm(index_t,index_t,index_t,index_t);
        
        /*the iterative backend of build_tree*/
        void organize(array_1d<index_t>&,index_t,index_t,index_t,index_t,index_t);
        
        /*find the node where a new point would be added to the tree;
        this is part of the backend for the add() routine*/
        index_t find_node(const array_1d<value_t>&);
        
        /*neigh_check provides the back end for nn_srch*/
        void neigh_check(const array_1d<value_t>&,
            index_t,array_1d<index_t>&,array_1d<value_t>&,index_t,index_t);
        
        /*kernel_check provides the backend for kernel_srch*/
        void kernel_check(array_1d<value_t>&,array_1d<value_t>&,array_1d<index_t>&,
            index_t,index_t);
        
        /*radial_check provides the backend for radial_srch*/
        void radial_check(array_1d<value_t>&,value_t,array_1d<index_t>&,index_t,index_t);
        
        /*
        reassign() and descend() provide some of the backend for remove()
        */
        void reassign(index_t);
        void descend(index_t);

};


void convert_to_boundary(array_2d<value_t>&, value_t, value_t, array_2d<value_t>&);

#endif
