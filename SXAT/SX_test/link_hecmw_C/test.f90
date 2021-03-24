program main
    Use iso_c_binding

    interface
        Function get_aveo_version() bind(c)
            Import
            Integer(C_int) get_aveo_version
        end Function
    end interface

    Integer(C_int) :: version
    version = get_aveo_version()
    print *, version

end program main
