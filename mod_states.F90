module mod_states
! Modelstate definition for CICEv5.1 model
use ice_kinds_mod
use mod_dimensions
use ice_blocks, only: nx_block, ny_block
use ice_domain, only: nblocks
use ice_domain_size, only: nx_global, ny_global, max_blocks

   type states
      real var
   end type states

   type gstates
      real var(nx_global,ny_global)
   end type gstates

   type states4
      real*4 var
   end type states4

! Overloaded and generic operators
   interface operator(+)
      module procedure add_states
   end interface

   interface operator(-)
      module procedure subtract_states
   end interface

   interface operator(*)
      module procedure states_real_mult,&
                       real_states_mult,&
                       states_states_mult
   end interface

   interface assignment(=)
      module procedure assign_states
      module procedure states4to8
      module procedure states8to4
   end interface

contains

  function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       add_states%var = A%var + B%var
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       subtract_states%var = A%var - B%var
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_mult%var = B*A%var
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_mult%var = B*A%var
   end function real_states_mult

  function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_mult%var = A%var * B%var
   end function states_states_mult


   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%var = r
   end subroutine assign_states

   subroutine states4to8(A,B)
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%var=DBLE(B%var)
   end subroutine states4to8

   subroutine states8to4(A,B)
      type(states), intent(in)  :: B
      type(states4),  intent(out) :: A
      A%var=real(B%var)
   end subroutine states8to4



end module mod_states

