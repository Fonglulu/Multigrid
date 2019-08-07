# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:49:15 2013

@author: stals

This file contains the routines that may be used to communicate information
between different processors.

It currently uses the lowercase commands from mpi4py, which means that 
everything is pickeled before being passed to MPI. This is simple and a 
quick way to set-up the MPI calls, but it is NOT efficient. Before doing any
serious efficiency testing the lowercase mpi4py commands must be replaced by
appropriate uppercase mpi4py commands.

The FunctionStore class is a bit of hack to get around the fact that MPI 
can't send a function. It is not needed with python. Using python's pickle
commands will allow functions to be sent between processors. I have decided to
keep on using the FunctionStore class because the ability to pickle a function
is so unique to python, meaning the ideas taught in this course would not 
carry through to another language. If the following code is improved so that
it does not use pickle, the FunctionStore should be used to help communicate
a function between processors.

The communication model that we are using is one of host (or master or boss)
and workers. The host is expected to do very little work, it essentially just
farms the original coarse grid to the works and waits for the results to be
send back. The works do most of the computational work. Several of the
communication calls are just designed for the worker group 

"""

# Import appropriate information from other classes

from mpi4py import MPI
#from function.FunctionStore import FunctionStore
from grid.NghNodes import ngh_iterator
from grid.NodeTable import NodeTable
from grid.NodeTable import node_iterator
from grid.EdgeTable import endpt_iterator


class Communication:
    """Communication between different processors"""
    
    # MPI variables
    
    # MPI communcation class
    _comm = 0
    
    # MPI groups
    _group = 0

    # MPI rank
    _rank = 0
    
    # MPI total number of processors
    _size = 0
    
    # MPI communication tags
    _tags = 0
    
    # MPI taks ids
    _tids = []
    
    # Number of workers
    _no_workers = 0
    
    # The host (or master) processor
    _host_no = 0
    
    # The id of the current processor
    _my_no = 0
    
    # Communication group amongst the worker processors
    _worker_comm = 0
    _worker_group = 0


    def __init__(self):
        """Initialise the MPI communication protocols
        
        Set up a host/worker environment with the number of
        workers being extracted from the mpirun call
        
        """
        
        # Extract the world communication information from MPI
        self._comm = MPI.COMM_WORLD
        self._rank = self._comm.Get_rank()
        self._size = self._comm.Get_size()
        self._group = self._comm.Get_group()
        
        # Set up default tag numbers
        self._tag = 11
        self._tag_worker = 12
        
        # Set up a host/worker structure
        
        # The number workers is the number of processors - 1
        self._no_workers = self._size-1
        
        # Find an id for the host processor
        self._host_no = self._no_workers
        self._tids = [0]*self._size
        self._tids[self._host_no] = 0
        
        # Find the id for the current processor
        if (self._rank == 0):
            self._my_no = self._host_no
        else:
            self._my_no = self._rank-1
            
        # Find the ids for the worker processors
        for i in range(self._no_workers):
            self._tids[i] = i+1
            
        # Create a sub group that just communicates amongst the workers
        cell_index  = self._tids[0:self._no_workers]
        self._worker_group = self._group.Incl(cell_index)
        self._worker_comm = self._comm.Create(self._worker_group)


    def get_my_no(self):
        """Get the current processor id"""
        return self._my_no
        
    def get_host_no(self):
        """Get the id of the host processor"""
        return self._host_no
        
    def get_no_workers(self):
        """Get the number of workers"""
        return self._no_workers
        
    def send_grid(self, processor_no, grid):
        """
        Send a grid to the given processor no
        
        This applies to the global group
        
        """
        
        # Send the (pickled) grid to the neighbouring processor
        self._comm.send(grid, dest=self._tids[processor_no], tag=self._tag)
       

    def get_grid(self, processor_no):
        """
        Get a grid from the given processor no 
        
        This applies to the global group
        
        """
        
        # Get the grid from the neighbouring processor
        grid = self._comm.recv(source=self._tids[processor_no], tag=self._tag)
        
        # When python pickles the grid it does not appear to handle the
        # static variables correctly. In particular, the counter in NodeTable
        # used to determine the local ID is not set correctly. The
        # We must make sure that counter is greater than the maximum
        # local ID value so that if a new node is added it will be given
        # the correct local ID
        
        # Find the maximum counter in the full nodes
        max_counter = 0
        for node in node_iterator(grid):
            node_id = node.get_node_id()
            if (node_id.get_no() > max_counter):
                max_counter = node_id.get_no()
                
        # Find the maximum counter in the full and ghost nodes
        ghost_table = grid.reference_ghost_table()
        for node in node_iterator(ghost_table):
            node_id = node.get_node_id()
            if (node_id.get_no() > max_counter):
                max_counter = node_id.get_no()
                
        # Adjust the static variable in the NodeTable accordingly
        if max_counter > NodeTable._counter:
            NodeTable._counter = max_counter+1

        # Return the grid
        return grid
        
    def send_function(self, processor_no, function):
        """
        Send a function to the given processor 
        
        This applies to the global group

        """
        
        # Send a pickled copy of the function
        self._comm.send(function, dest=self._tids[processor_no], tag=self._tag)

        
    def get_function(self, processor_no):
        """
        Get a function from the given processor
        
        This applies to the global group

        """     
        
        # Get the copy of the function
        function = self._comm.recv(source=self._tids[processor_no], 
                                   tag=self._tag)

        return function
        
    def update_vector_table(self, grid, vector_table, ghost_vector_table):
        """Update the ghost node information
        
        This applies to the worker group

        """
        
        # Get information about where the updates have to be sent to
        full_neighbour_table = grid.reference_full_commun()
        
        # Where will the updates come from?
        ghost_neighbour_table = grid.reference_ghost_commun()
        
        # The updated information will be stored in the ghost table
        ghost_table = grid.reference_ghost_table()
        
        # Loop over the neighbouring processor
        send_list = []
        reqs = []
        for processor_no, id_list in ngh_iterator(full_neighbour_table):
            
            # Make a record of the information that needs to be sent
            send_vector_table = dict()
            for node_id in id_list:
                global_id = grid.get_global_id(node_id)
                send_vector_table[str(global_id)]=vector_table[str(node_id)]
                
            # Send that information
            reqs.append(self._comm.isend(send_vector_table, 
                             dest=self._tids[int(processor_no)],                                                                
                            tag=self._tag))
                            
            # Do not delete the information until all of the works have
            # completed
            send_list.append(send_vector_table)

        # Loop over the neighbouring processors
        for processor_no, id_list in ngh_iterator(ghost_neighbour_table):
                        
            # Get the updated information
            get_vector_table = self._comm.recv(
                source=self._tids[int(processor_no)], 
                tag=self._tag)
                            
            # Update the information in the ghost node table
            for global_id in get_vector_table:
                node_id = ghost_table.get_local_id(global_id)
                ghost_vector_table[str(node_id)]=get_vector_table[global_id]
                
        self._worker_comm.barrier()   
                

             
    def send_array(self, processor_no, array):
        """
        Send an array to the given processor
        
        This applies to the global group
        """
        self._comm.send(array, dest=self._tids[processor_no], tag=self._tag)
         
    def get_array(self, processor_no):
        """
        Get an array from the given processor
        
        This applies to the global group
        
        """
        
        array = self._comm.recv(source=self._tids[processor_no], 
                                tag=self._tag)
        return array        
        
    def global_sum(self, my_value):
        """Find a global sum
        
        This applies to the worker group

        """
        from numpy import array
        
        # Use the reduce routine to find the sum
        sendBuf = array(my_value, 'd')
        recvBuf = array(0.0, 'd')
        self._worker_comm.Allreduce([sendBuf, MPI.DOUBLE], 
                                    [recvBuf, MPI.DOUBLE], 
                                      op = MPI.SUM)
                                      
        # Return the sum
        return recvBuf
        
    def all_gather(self, my_value):
        """all to all communication
        
        This applies to the worker group
        
        """
        
        from numpy import array

        # Use the gather routine to get the information from other processor
        sendBuf = array(my_value, 'd')
        recvBuf = array([0.0]*self._no_workers, 'd')
        self._worker_comm.Allgather([sendBuf, MPI.DOUBLE], 
                                    [recvBuf, MPI.DOUBLE])
                                    
        # Return an array containing the value from each processor
        return recvBuf
        

 