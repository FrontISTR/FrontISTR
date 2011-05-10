module m_mw3_interface
   use, intrinsic :: iso_c_binding

   interface
     integer(C_INT) function mw_solve_( iter_max, tolerance, method, pre_condition ) bind(c)
        use iso_c_binding
        implicit none
        type(C_PTR), VALUE :: iter_max
        type(C_PTR), VALUE :: tolerance
        type(C_PTR), VALUE :: method
        type(C_PTR), VALUE :: pre_condition
	 END FUNCTION mw_solve_
	 
	 subroutine mw_allreduce_r_(rval, val_size, op ) bind(c)
	    use iso_c_binding
        implicit none
        type(C_PTR), VALUE :: rval
        type(C_PTR), VALUE :: val_size
        type(C_PTR), VALUE :: op
	 END subroutine mw_allreduce_r_
   end interface
   
   
   
end module