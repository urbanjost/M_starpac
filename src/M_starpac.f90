module M_starpac
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, M_starpac!"
  end subroutine say_hello
end module M_starpac
