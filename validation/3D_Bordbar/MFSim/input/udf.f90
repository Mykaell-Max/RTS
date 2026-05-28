
subroutine lb_calculate_weights(mesh_ptr, headlev_ptr, igxb, igyb, igzb, ratio)
    use lb_udf
    use data_grid
    use data_grid_obj 
    use par_data_types
    implicit none
    ! In
    type(Mesh_Type), pointer :: mesh_ptr
    integer :: igxb, igyb, igzb, ratio
    type(level_components1), dimension(:), pointer :: headlev_ptr
    ! Local
    integer :: l

    do l = lbot+2,lbot+2 !laux,lbot+1,-1 !lbot,-1 !ltop,lbot+1,-1
        call set_wgts(l,headlev_ptr, mesh_ptr, igxb, igyb, igzb, ratio)
    end do
end subroutine lb_calculate_weights