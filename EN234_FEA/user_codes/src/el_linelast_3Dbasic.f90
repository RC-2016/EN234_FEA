!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,kk

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12             ! Material properties
    real (prec)  ::  B_correction(6, length_dof_array) ! B_correction
    real (prec)  ::  B_bar(6,length_dof_array)         ! B_bar calculation
    real (prec)  ::  el_vol                            ! elemental volume
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!10/26/15 - hw6
    !!determining inf the input file is nonlinear
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(6,*) n_properties

    IF (n_properties == 4) THEN
        write(6,*) 'hypoelastic material based on input file'
    else
        write(6,*) 'linear elastic material based on input file'
    END IF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Checking against which element type created
    !Added do loop to compute dNbar

    IF (element_identifier == 1002) THEN

        !if statement here to determine if linear or not

        write(6,*) 'element is 1002'

        write(6,*) 'enter linear elastic mode'

        dNbardx = 0.d0
        el_vol = 0.d0

        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            el_vol = el_vol + w(kint)*determinant

            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) +  dNdx(1:n_nodes,1:3) * w(kint) * determinant

        end do

        dNbardx = dNbardx/el_vol

    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    !     --  Loop over integration points
    do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Added B_correction 10/20
        !Added b_bar and if statement for element identifier
        !compute new strain, stress values based with b_bar

        B_correction = 0.d0

        if (element_identifier == 1002) then

!            write(6,*) 'enter linear elastic model b matrix model calculations'

            do kk = 1,n_nodes
                B_correction(1,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
                B_correction(2,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
                B_correction(3,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
            end do

            B_bar = B + B_correction*(1.d0/3.d0)

            strain = matmul(B_bar,dof_total)
            dstrain = matmul(B_bar,dof_increment)

            if (n_properties == 4) then

                write(6,*) 'material call, hypoelastic'

                call hypoelastic_material(strain, dstrain, n_properties, element_properties, &
                 stress, D)

            else
                stress = matmul(D,strain+dstrain)
            end if

            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(D,B_bar(1:6,1:3*n_nodes)))*w(kint)*determinant

            write(6,*) ' 1002 calcuations complete '

        else

            !if it is not a 1002 element or hypoelastic material than it is normal linear elastic behavior

            !original code
            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            stress = matmul(D,strain+dstrain)
            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

            write(6,*) 'element stiffness below'
            write(6,*) element_stiffness

            write(6,*) ' 1001 calcuations complete '

        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!my notes for changes and computations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! Changed this to for l
        !calc b_bar


!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !adjusting ouput to compute with B_bar
!
!        stress = matmul(D,strain+dstrain)
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress)*w(kint)*determinant
!
!        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
!            + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(D,B_bar(1:6,1:3*n_nodes)))*w(kint)*determinant
!
!        write(6,*) ' element stiffness complete '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !original code
!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)
!
!        stress = matmul(D,strain+dstrain)
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
!
!        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
!            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do

  
    return
end subroutine el_linelast_3dbasic


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute hypoelastic material model
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================SUBROUTINE hypoelastic_material ==============================
subroutine hypoelastic_material(strain, dstrain, n_properties, element_properties, stress, D) ! Output variables are stress and D

    use Types
    use ParamIO
    use Mesh, only : node


! verify this is the correct format
!   passed variables
    real( prec ), intent( in )    :: strain(6)
    real( prec ), intent( in )    :: dstrain(6)
    integer, intent( in )         :: n_properties
    real( prec ), intent( in )    :: element_properties(n_properties)  ! Element or material properties, stored in order listed in input file
    real( prec ), intent( out )   :: stress(6)
    real( prec ), intent( out )   :: D(6,6)


!   local variables
    real( prec )    ::  total_strain_vector(6)                     !strain vector (strain + dstrain)
    !real( prec )    ::  stress_vector(6)                          !stress vector this will be calculated
    real( prec )    ::  devi_strain(6)                             !deviatoric strain tensor
    real( prec )    ::  von_mises_strain                           !von mises strain
    real( prec )    ::  E_s                                        !
    real( prec )    ::  E_t                                        !
    real( prec )    ::  n,k                                          !element numbers - material prop
    real( prec )    ::  strain_int                                 !strain initial - material prop
    real( prec )    ::  sigma_int                                  !stress initial - material prop
    real( prec )    ::  sigma_e
    real( prec )    ::  d_sigma_e
    real( prec )    ::  d_sigma_e_bottom
    real( prec )    ::  d_sigma_e_top
    real( prec )    ::  id_matrix_1(6,6)
    real( prec )    ::  id_matrix_2(6,6)
    real( prec )    ::  D_1(6,6)
    real( prec )    ::  D_2(6,6)
    real( prec )    ::  e_dyadic_e(6,6)


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!Variable block!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

    sigma_int = element_properties(1)
    strain_int = element_properties(2)
    n = element_properties(3)
    k = element_properties(4)

!    write(6,*) sigma_int
!    write(6,*) strain_int
!    write(6,*) n
!    write(6,*) k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!Calculations!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   total strain vector
    total_strain_vector = (strain + dstrain)

    write(6,*) 'tsv',total_strain_vector

!   deviatoric strain
    devi_strain(4:6) = total_strain_vector(4:6)/2.d0
    devi_strain(1:3) = total_strain_vector(1:3) - (1.d0/3.d0)*sum(total_strain_vector(1:3))

    write(6,*) ' devi',devi_strain

!    write(6,*) devi_strain

!   von mises strain
!!!!!!!!!!!!! CHECK HERE ???? !!!!!!!!
    von_mises_strain = dot_product(devi_strain(1:3) , devi_strain(1:3)) + 2.d0*dot_product(devi_strain(4:6),devi_strain(4:6))
    von_mises_strain = dsqrt((2.d0/3.d0)*von_mises_strain)

    write(6,*) 'von mises strain below'
    write(6,*) von_mises_strain


    !solve for sigma_e
    if  (von_mises_strain < strain_int) then

        sigma_e = dsqrt( (1.d0+n**2.d0)/((n-1.d0)**2.d0) - ( n/(n-1.d0) - von_mises_strain / strain_int )**2.d0 )
        sigma_e = sigma_e - (1.d0/(n-1.d0))
        sigma_e = sigma_e * sigma_int

        write(6,*) 'von mises is less then strain initial - sigma_e below'
        write(6,*) sigma_e

    else

        sigma_e = sigma_int * (von_mises_strain / strain_int)**(1/n)

        write(6,*) 'von mises is greater than strain int - sigma_e below'
        write(6,*) sigma_e

    end if


!   solving for stress vector
    if (von_mises_strain == 0.d0) then
        stress = 0.d0
        stress(1:3) = k*sum(total_strain_vector(1:3))
    else
        stress(1:3) = (2.d0/3.d0)*sigma_e*(devi_strain(1:3) / von_mises_strain) + k * sum(total_strain_vector(1:3))
        stress(4:6) = (2.d0/3.d0)*sigma_e*(devi_strain(4:6) / von_mises_strain)
    endif

!   ??????????????????
!   E_t
    if (von_mises_strain == 0.d0) then

        d_sigma_e_top = ((n/(n-1.d0)) - (von_mises_strain/strain_int))
        d_sigma_e_bottom = dsqrt( (1.d0+n**2)/(n-1.d0)**2.d0 - (n/(n-1.d0) - von_mises_strain / strain_int)**2.d0)
        d_sigma_e_bottom = strain_int * d_sigma_e_bottom
        d_sigma_e = sigma_int*(d_sigma_e_top / d_sigma_e_bottom)

        E_t = d_sigma_e

    else if (von_mises_strain > strain_int) then

        write(6,*) 'von mises less that intial strain for derivative calc'

        d_sigma_e = sigma_e/n/von_mises_strain
        E_t = d_sigma_e

        write(6,*) 'HERE - E_t value is below - top part of code'
        write(6,*) E_t

    else

        write(6,*) 'von mises greater than initial strain for deriviatve calc'

        d_sigma_e_top = ((n/(n-1.d0)) - (von_mises_strain/strain_int))
        d_sigma_e_bottom = dsqrt( (1.d0+n**2)/(n-1.d0)**2 - (n/(n-1.d0) - von_mises_strain / strain_int)**2)
        d_sigma_e_bottom = strain_int * d_sigma_e_bottom
        d_sigma_e = sigma_int*(d_sigma_e_top / d_sigma_e_bottom)

        E_t = d_sigma_e

        write(6,*) 'HERE - E_t value is below - bottom part of code'
        write(6,*) E_t

    end if

    !!!!!calculating E_s
    if (von_mises_strain == 0) then

        write(6,*) 'von mises strain is zero - E_s = E_t'

        E_s = E_t

        write(6,*) 'E_s value below'
        write(6,*) E_s

    else

        E_s = sigma_e / von_mises_strain

        write(6,*) 'E_s value below'
        write(6,*) E_s

    end if

    !!!!calculation for the D matrix
    !matrix building
    id_matrix_1 = 0.d0
    id_matrix_1(1,1) = 2.d0
    id_matrix_1(2,2) = 2.d0
    id_matrix_1(3,3) = 2.d0
    id_matrix_1(4,4) = 1.d0
    id_matrix_1(5,5) = 1.d0
    id_matrix_1(6,6) = 1.d0

    id_matrix_2 = 0.d0
    id_matrix_2(1:3,1) = 1.d0
    id_matrix_2(1:3,2) = 1.d0
    id_matrix_2(1:3,3) = 1.d0

    !!!Building the D matrix
    !two cases et - es = 0 or doesnt
    if (von_mises_strain == 0.d0) then
        E_s = E_t

        D = (E_s / 3.d0)*id_matrix_1 + (k-((2.d0*E_s)/9.d0))*id_matrix_2
    else

        e_dyadic_e = spread(devi_strain,dim=2,ncopies=6)*spread(devi_strain,dim=1,ncopies=6)
        D_1 = ((4.d0 /(9.d0*von_mises_strain**2)) * (E_t - E_s)) * e_dyadic_e

        D_2 = (E_s / 3.d0) * id_matrix_1 + (k - ((2.d0*E_s) / 9.d0))*id_matrix_2
        D = D_1 + D_2

    end if

    !stress and D matrix calculations complete

!     write(6,*) 'haha - your code ran but your values need work'
    write(6,*) 'hypo stress and D matrix calcs complete'

    return
end subroutine hypoelastic_material











    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do

    return
end subroutine el_linelast_3dbasic_dynamic





!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,kk

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    real (prec)  ::  B_correction(6, length_dof_array) ! B_correction
    real (prec)  ::  B_bar(6,length_dof_array)         ! B_bar calculation
    real (prec)  ::  el_vol                            ! elemental volume

    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Checking against which element type created
    !Added do loop to compute dNbar

    IF (element_identifier == 1002) THEN

        write(6,*) 'element is 1002 - field vars'

        dNbardx = 0.d0
        el_vol = 0.d0

        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            el_vol = el_vol + w(kint)*determinant

            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) +  dNdx(1:n_nodes,1:3) * w(kint) * determinant

        end do

        dNbardx = dNbardx/el_vol

    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Added B_correction 10/20

        B_correction = 0.d0

        if (element_identifier == 1002) then

            do kk = 1,n_nodes
                B_correction(1,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
                B_correction(2,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
                B_correction(3,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
            end do

            !calc b_bar
            B_bar = B + B_correction*(1.d0/3.d0)

            strain = matmul(B_bar,dof_total)
            dstrain = matmul(B_bar,dof_increment)

            write(6,*) 'Field vars complete for element 1002'

        else

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            write(6,*) 'Field vars subroutine compelete for element 1001'

        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !        strain = matmul(B,dof_total)
        !        dstrain = matmul(B,dof_increment)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call to hypo subroutine to get new D matrix
        if (n_properties == 4) then
            write(6,*) 'material call, hypoelastic'
            call hypoelastic_material(strain, dstrain, n_properties, element_properties, &
                stress, D)
        else
            stress = matmul(D,strain+dstrain)
        end if
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        stress = matmul(D,strain+dstrain)
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_linelast_3dbasic





