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
    integer      :: n_points,kint,kk,i

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
    real (prec)  ::  dNdy(length_dof_array/3,3)          ! shape funciton wrt to y
    real (prec)  ::  delta_ij(3,3)                     !used to calc the new shap function
    real (prec)  ::  f_ij(3,3)                         !deformation gradient matrix for finding y shape function
    real (prec)  ::  f_ij_inv(3,3)                         !invesre of f
    real (prec)  ::  J                                !scalar value for determinant of f
    real (prec)  ::  u_i(3, length_dof_array/3)
    real (prec)  ::  B_star(9,3*n_nodes)        !maps nodal velocities to the vel gradient
    real (prec)  ::  G(6,9)
    real (prec)  ::  S(3,length_dof_array/3)
    real (prec)  ::  Pvec(length_dof_array)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  Svec(length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array)
    real (prec)  ::  Sigma(length_dof_array,length_dof_array)
!    real (prec) ::   temp(length_dof_array)

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

!    write(6,*) n_properties

!    IF (n_properties == 4) THEN
!        write(6,*) 'hypoelastic material based on input file'
!    else
!        write(6,*) 'linear elastic material based on input file'
!    END IF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Checking against which element type created
    !Added do loop to compute dNbar

!    IF (element_identifier == 1002) THEN
!
!        !if statement here to determine if linear or not
!
!        write(6,*) 'element is 1002'
!
!        write(6,*) 'enter linear elastic mode'
!
!        dNbardx = 0.d0
!        el_vol = 0.d0
!
!        do kint = 1, n_points
!            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!
!            call invert_small(dxdxi,dxidx,determinant)
!            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!
!            el_vol = el_vol + w(kint)*determinant
!
!            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) +  dNdx(1:n_nodes,1:3) * w(kint) * determinant
!
!        end do
!
!        dNbardx = dNbardx/el_vol
!
!    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    element_residual = 0.d0
    element_stiffness = 0.d0
	
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !creating u - hw7
    u_i = 0.d0
    u_i = reshape(dof_total+dof_increment,(/3,length_coord_array/3/))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !     --  Loop over integration points
    do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call hyperelastic material HW7

        !idenity matrix to compute the deformation gradient f
        delta_ij = 0.d0
        delta_ij(1,1) = 1.d0
        delta_ij(2,2) = 1.d0
        delta_ij(3,3) = 1.d0

        !creating the f matrix - u call is outside of the loop
        f_ij(1:3,1:3) = delta_ij + (matmul(u_i(1:3,1:n_nodes), dNdx(1:n_nodes,1:3)))

!        write(6,*) 'solving f_ij below'

         !compute the inverse of the deformation gradient and J
        call invert_small(f_ij,f_ij_inv,J)


        !computing the shape function to the right coordinates (y)
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),f_ij_inv(1:3,1:3))
!        write(6,*) 'solving dNdy to the right'

        !need to remake my b matrix wrt y meow (yes that's a cat joke)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

!        write(6,*) 'b matrix adjusted wrt y'


        !create the b* function
        B_star = 0.d0
        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)


!        write(6,*) 'b_star matrix had been created'
!        write(6,*) 'call to hyperelastic material'

        call hyperelastic_material(n_properties, element_properties, f_ij, G, D, stress) !d and stress

!        temp = matmul(transpose(B(1:6,1:n_nodes)),stress(1:6))
!        write(6,*) temp
!        S = reshape(temp,(/3,length_dof_array/3/))

        Sigma = 0.d0
        S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do

        Sigma = Pmat*transpose(Smat)
!        write(6,*) Sigma

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                    + matmul(transpose(B),matmul(matmul(D,G),B_star))*w(kint)*determinant &
                    - Sigma*w(kint)*determinant


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Added B_correction 10/20
        !Added b_bar and if statement for element identifier
        !compute new strain, stress values based with b_bar

!        B_correction = 0.d0
!
!        if (element_identifier == 1002) then
!
!!            write(6,*) 'enter linear elastic model b matrix model calculations'
!
!            do kk = 1,n_nodes
!                B_correction(1,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!                B_correction(2,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!                B_correction(3,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!            end do
!
!            B_bar = B + B_correction*(1.d0/3.d0)
!
!            strain = matmul(B_bar,dof_total)
!            dstrain = matmul(B_bar,dof_increment)
!
!            if (n_properties == 4) then
!
!                write(6,*) 'material call, hypoelastic'
!
!                call hypoelastic_material(strain, dstrain, n_properties, element_properties, &
!                 stress, D)
!
!            else
!                stress = matmul(D,strain+dstrain)
!            end if
!
!            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress)*w(kint)*determinant
!
!            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
!                + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(D,B_bar(1:6,1:3*n_nodes)))*w(kint)*determinant
!
!            write(6,*) ' 1002 calcuations complete '
!
!        else
!
!            !if it is not a 1002 element or hypoelastic material than it is normal linear elastic behavior
!
!            !original code
!!            strain = matmul(B,dof_total)
!!            dstrain = matmul(B,dof_increment)
!!
!!            stress = matmul(D,strain+dstrain)
!!            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
!!
!!            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
!!                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant
!!
!!            write(6,*) 'element stiffness below'
!!            write(6,*) element_stiffness
!
!            write(6,*) ' 1001 calcuations complete '
!
!        endif

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
!compute hyperelastic material model - hw7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hyperelastic_material(n_properties, element_properties, f_ij, G, D, stress)

    use Types
    use ParamIO
    use Mesh, only : node

    use Element_Utilities, only : invert_small

    implicit none
    !passed variables here
    integer, intent(in)            ::    n_properties
    real( prec ), intent( in )    ::    element_properties(n_properties)
    real( prec ), intent( in )    ::    f_ij(3,3)
    real( prec ), intent( out )   ::    G(6,9)
    real( prec ), intent( out )   ::    stress(6)
    real( prec ), intent( out )   ::    D(6,6)

    !local variables here
!    real( prec )    ::  G(6,9)                      !used in derivation of D
    real( prec )    ::  B_ij(3,3)                   ! not the right size and will need to adjust
    real( prec )    ::  calc_matrix(6,6)            !this is used in calculation of the D matrix
    real( prec )    ::  I_bar(6)                  !this is used in calculation of the D matrix
    real( prec )    ::  B_bar(6)                  !this is used in calculation of the D matrix
    real( prec )    ::  mu                          !this is from the element properties
    real( prec )    ::  k                           !this is from the element properties input file
    real( prec )    ::  D_start(6,6)
    real( prec )    ::  D_mid(6,6)
    real( prec )    ::  D_end(6,6)
    real( prec )    ::  f_ij_inv(3,3)
    real( prec )    ::  J
    real( prec )    ::  B_ij_trace
    real( prec )    ::  ib_dydaic(6,6)
    real( prec )    ::  ii_dydaic(6,6)
    real( prec )    ::  bb_dydaic(6,6)
!    real( prec )    ::  stress_left(6)
!    real( prec )    ::  stress_right(6)
    real( prec )    ::  B_kk
    real( prec )    ::  delta_ij(3,3)
    real( prec )    ::  B_ij_inv(3,3)
    real( prec )    ::  J_B_ij_inv
    real( prec )    ::  B_bar_inverse(6)


    !element properties
    mu = element_properties(1)
    k = element_properties(2)

    !idenity matrix to compute the deformation gradient f
    delta_ij = 0.d0
    delta_ij(1,1) = 1.d0
    delta_ij(2,2) = 1.d0
    delta_ij(3,3) = 1.d0


    !stating what the calc_matrix is
    calc_matrix = 0.d0
    calc_matrix(1,1) = 1.d0
    calc_matrix(2,2) = 1.d0
    calc_matrix(3,3) = 1.d0
    calc_matrix(4,4) = .5d0
    calc_matrix(5,5) = .5d0
    calc_matrix(6,6) = .5d0

    !calc I_bar and B_bar vectors for calc of D matrix
    I_bar = 0.d0
    I_bar(1) = 1.d0
    I_bar(2) = 1.d0
    I_bar(3) = 1.d0

    !solving for adjusted b matrix
    B_ij = matmul(f_ij, transpose(f_ij))
    B_bar(1) = B_ij(1,1)
    B_bar(2) = B_ij(2,2)
    B_bar(3) = B_ij(3,3)
    B_bar(4) = B_ij(1,2)
    B_bar(5) = B_ij(1,3)
    B_bar(6) = B_ij(2,3)

    !computing j
    call invert_small(f_ij,f_ij_inv,J)

    B_kk = B_ij(1,1) + B_ij(2,2) + B_ij(3,3)

    !calculating stress
!    stress = (mu/(J**(2.d0/3.d0))) * (B_bar - ((1.d0/3.d0)*B_kk*I_bar)) + (K*J*(J-1)*I_bar) !bower adjustment divide by j when plotting
!    stress = (mu/(J**(5.d0/3.d0))) * (B_bar - ((1.d0/3.d0)*B_kk*I_bar)) + (K*(J-1)*I_bar)
!    stress = J*stress
   stress(1) = (mu/(J**(5.d0/3.d0))) * (B_bar(1) - ((1.d0/3.d0)*B_kk*I_bar(1))) + (K*(J-1.d0))
   stress(2) = (mu/(J**(5.d0/3.d0))) * (B_bar(2) - ((1.d0/3.d0)*B_kk*I_bar(2))) + (K*(J-1.d0))
   stress(3) = (mu/(J**(5.d0/3.d0))) * (B_bar(3) - ((1.d0/3.d0)*B_kk*I_bar(3))) + (K*(J-1.d0))
   stress(4) = (mu/(J**(5.d0/3.d0))) * (B_bar(4))
   stress(5) = (mu/(J**(5.d0/3.d0))) * (B_bar(5))
   stress(6) = (mu/(J**(5.d0/3.d0))) * (B_bar(6))

!   stress = stress * J

!    write(6,*) 'stress calculated'

    !creating the G matrix the hard way
    G = 0.d0

    !first row
    G(1,1) = 2.d0*B_ij(1,1)
    G(1,4) = 2.d0*B_ij(1,2)
    G(1,6) = 2.d0*B_ij(1,3)
    !second row
    G(2,2) = 2.d0*B_ij(2,2)
    G(2,5) = 2.d0*B_ij(1,2)
    G(2,8) = 2.d0*B_ij(2,3)
    !third row
    G(3,3) = 2.d0*B_ij(3,3)
    G(3,7) = 2.d0*B_ij(1,3)
    G(3,9) = 2.d0*B_ij(2,3)
    !fourth row
    G(4,1) = 2.d0*B_ij(1,2)
    G(4,2) = 2.d0*B_ij(1,2)
    G(4,4) = 2.d0*B_ij(2,2)
    G(4,5) = 2.d0*B_ij(1,1)
    G(4,6) = 2.d0*B_ij(2,3)
    G(4,8) = 2.d0*B_ij(1,3)
    !fifth row
    G(5,1) = 2.d0*B_ij(1,3)
    G(5,3) = 2.d0*B_ij(1,3)
    G(5,4) = 2.d0*B_ij(2,3)
    G(5,6) = 2.d0*B_ij(3,3)
    G(5,7) = 2.d0*B_ij(1,1)
    G(5,9) = 2.d0*B_ij(1,2)
    !sixth row
    G(6,2) = 2.d0*B_ij(2,3)
    G(6,3) = 2.d0*B_ij(2,3)
    G(6,5) = 2.d0*B_ij(1,3)
    G(6,7) = 2.d0*B_ij(1,2)
    G(6,8) = 2.d0*B_ij(3,3)
    G(6,9) = 2.d0*B_ij(2,2)


    !calculating D
    call invert_small(B_ij, B_ij_inv, J_B_ij_inv)

    B_bar_inverse(1) = B_ij_inv(1,1)
    B_bar_inverse(2) = B_ij_inv(2,2)
    B_bar_inverse(3) = B_ij_inv(3,3)
    B_bar_inverse(4) = B_ij_inv(1,2)
    B_bar_inverse(5) = B_ij_inv(1,3)
    B_bar_inverse(6) = B_ij_inv(2,3)

    B_ij_trace = B_ij(1,1) + B_ij(2,2) + B_ij(3,3)

    ib_dydaic = spread(I_bar,dim=2,ncopies=6)*spread(B_Bar_inverse,dim=1,ncopies=6)
    ii_dydaic = spread(I_bar,dim=2,ncopies=6)*spread(I_bar,dim=1,ncopies=6)
    bb_dydaic = spread(B_bar,dim=2,ncopies=6)*spread(B_Bar_inverse,dim=1,ncopies=6)

    D_start = (mu/(J**(2.d0/3.d0))) * calc_matrix
    D_mid = (mu/(3*J**(2.d0/3.d0)))*(((B_ij_trace)/3)*ib_dydaic - ii_dydaic - bb_dydaic)
    D_end = K*J*(J-(1.d0/2.d0))*ib_dydaic

    D = D_start + D_mid + D_end


!    write(6,*) 'exiting hyperelastic routine yall'

    return

end subroutine hyperelastic_material




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
    real (prec)  ::  delta_ij(3,3)                     !used to calc the new shap function
    real (prec)  ::  f_ij(3,3)                         !deformation gradient matrix for finding y shape function
    real (prec)  ::  f_ij_inv(3,3)                         !invesre of f
    real (prec)  ::  J                                !scalar value for determinant of f
    real (prec)  ::  u_i(3, length_dof_array/3)
!    real (prec)  ::  B_star(9,3*n_nodes)        !maps nodal velocities to the vel gradient
    real (prec)  ::  dNdy(length_dof_array/3,3)          ! shape funciton wrt to y
    real (prec)  ::  G(6,9)


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

!    IF (element_identifier == 1002) THEN
!
!        write(6,*) 'element is 1002 - field vars'
!
!        dNbardx = 0.d0
!        el_vol = 0.d0
!
!        do kint = 1, n_points
!            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!
!            call invert_small(dxdxi,dxidx,determinant)
!            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!
!            el_vol = el_vol + w(kint)*determinant
!
!            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) +  dNdx(1:n_nodes,1:3) * w(kint) * determinant
!
!        end do
!
!        dNbardx = dNbardx/el_vol
!
!    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nodal_fieldvariables = 0.d0
!
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !creating u - hw7
    u_i = 0.d0
    u_i = reshape(dof_total+dof_increment,(/3,length_coord_array/3/))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !hw7
          !idenity matrix to compute the deformation gradient f
        delta_ij = 0.d0
        delta_ij(1,1) = 1.d0
        delta_ij(2,2) = 1.d0
        delta_ij(3,3) = 1.d0

        !creating the f matrix - u call is outside of the loop
        f_ij(1:3,1:3) = delta_ij + (matmul(u_i(1:3,1:n_nodes), dNdx(1:n_nodes,1:3)))

!        write(6,*) 'solving f_ij below'

         !compute the inverse of the deformation gradient and J
        call invert_small(f_ij,f_ij_inv,J)
!        write(6,*) 'f_ij has been inverted'

        !computing the shape function to the right coordinates (y)
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),f_ij_inv(1:3,1:3))
!        write(6,*) 'solving dNdy to the right'

        !need to remake my b matrix wrt y meow (yes that's a cat joke)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)


        call hyperelastic_material(n_properties, element_properties, f_ij, G, D, stress) !d and stress



!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Added B_correction 10/20

!        B_correction = 0.d0
!
!        if (element_identifier == 1002) then
!
!            do kk = 1,n_nodes
!                B_correction(1,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!                B_correction(2,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!                B_correction(3,3*kk-2:3*kk) = (dNbardx(kk,1:3) - dNdx(kk,1:3))
!            end do
!
!            !calc b_bar
!            B_bar = B + B_correction*(1.d0/3.d0)
!
!            strain = matmul(B_bar,dof_total)
!            dstrain = matmul(B_bar,dof_increment)
!
!            write(6,*) 'Field vars complete for element 1002'
!
!        else
!
!            strain = matmul(B,dof_total)
!            dstrain = matmul(B,dof_increment)
!
!            write(6,*) 'Field vars subroutine compelete for element 1001'
!
!        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !        strain = matmul(B,dof_total)
        !        dstrain = matmul(B,dof_increment)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call to hypo subroutine to get new D matrix
!        if (n_properties == 4) then
!            write(6,*) 'material call, hypoelastic'
!            call hypoelastic_material(strain, dstrain, n_properties, element_properties, &
!                stress, D)
!        else
!!            stress = matmul(D,strain+dstrain)
!        end if
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





